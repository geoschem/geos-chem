! $Id: sulfate_mod.f,v 1.5 2003/12/05 21:14:05 bmy Exp $
      MODULE SULFATE_MOD
!
!******************************************************************************
!  Module SULFATE_MOD contains arrays and routines for performing either a
!  coupled chemistry/aerosol run or an offline sulfate aerosol simulation.
!  Original code taken from Mian Chin's GOCART model and modified accordingly.
!  (rjp, bdf, bmy, 6/22/00, 8/1/03)
!
!  Module variables:
!  ============================================================================
!  (1 ) XNUMOL_OH  (REAL*8 ) : Molecules OH  per kg OH          [molec/kg]
!  (2 ) XNUMOL_O3  (REAL*8 ) : Molecules O3  per kg O3          [molec/kg]
!  (3 ) XNUMOL_NO3 (REAL*8 ) : Molecules NO3 per kg NO3         [molec/kg]
!  (4 ) TCVV_S     (REAL*8 ) : Ratio: Molwt air / Molwt S       [unitless]
!  (5 ) DMSo       (REAL*8 ) : DMS oceanic emissions            [v/v/timestep]
!  (6 ) DRYH2O2    (INTEGER) : Pointer to H2O2 in DEPVEL array  [unitless] 
!  (7 ) DRYSO2     (INTEGER) : Pointer to SO2  in DEPVEL array  [unitless]
!  (8 ) DRYSO4     (INTEGER) : Pointer to SO4  in DEPVEL array  [unitless]
!  (9 ) DRYMSA     (INTEGER) : Pointer to MSA  in DEPVEL array  [unitless]
!  (10) DRYNH3     (INTEGER) : Pointer to NH3  in DEPVEL array  [unitless]
!  (11) DRYNH4     (INTEGER) : Pointer to NH4  in DEPVEL array  [unitless]
!  (12) DRYNIT     (INTEGER) : Pointer to NIT  in DEPVEL array  [unitless]
!  (13) ENH3_an    (REAL*8 ) : NH3 anthropogenic emissions      [kg NH3/box/s]
!  (14) ENH3_bb    (REAL*8 ) : NH3 biomass emissions            [kg NH3/box/s]
!  (15) ENH3_bf    (REAL*8 ) : NH3 biofuel emissions            [kg NH3/box/s]
!  (16) ENH3_na    (REAL*8 ) : NH3 natural source emissions     [kg NH3/box/s]
!  (17) ESO2_ac    (REAL*8 ) : SO2 aircraft emissions           [kg SO2/box/s]
!  (18) ESO2_an    (REAL*8 ) : SO2 anthropogenic emissions      [kg SO2/box/s]
!  (19) ESO2_ev    (REAL*8 ) : SO2 eruptive volcanic em.        [kg SO2/box/s]
!  (20) ESO2_nv    (REAL*8 ) : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
!  (21) ESO2_bb    (REAL*8 ) : SO2 biomass burning emissions    [kg SO2/box/s]
!  (22) ESO2_bf    (REAL*8 ) : SO2 biofuel burning emissions    [kg SO2/box/s]
!  (23) ESO4_an    (REAL*8 ) : SO4 anthropogenic emissions      [kg SO2/box/s]
!  (24) IJSURF     (INTEGER) : 1-D grid box indices             [unitless]
!  (25) JH2O2      (REAL*8 ) : Monthly mean J(H2O2) values      [s-1]
!  (26) O3m        (REAL*8 ) : Monthly mean O3 concentration    [v/v]
!  (27) PH2O2m     (REAL*8 ) : Monthly mean P(H2O2)             [molec/cm3/s]
!  (28) PMSA_DMS   (REAL*8 ) : P(MSA) from DMS                  [v/v/timestep]
!  (29) PSO2_DMS   (REAL*8 ) : P(SO2) from DMS                  [v/v/timestep]
!  (30) PSO4_SO2   (REAL*8 ) : P(SO4) from SO2                  [v/v/timestep]
!  (31) SSTEMP     (REAL*8 ) : Sea surface temperatures         [K]
!  (32) VCLDF      (REAL*8 ) : Volume cloud frac. for SO2 aq.   [unitless]
!  (33) NEV        (INTEGER) : Max # of eruptive volcanoes      [unitless]
!  (34) IEV        (INTEGER) : Longitudes of eruptive volcanoes [degrees]  
!  (35) JEV        (INTEGER) : Latitudes of eruptive volcanoes  [degrees ]
!  (36) IHGHT      (INTEGER) : Height of eruptive volcano plume [m]
!  (37) IELVe      (INTEGER) : Elevation of eruptive volcanoes  [m]
!  (38) Eev        (REAL*8 ) : SO2 em. from eruptive volcanoes  [kg SO2/box/s]
!  (39) NNV        (INTEGER) : Max # of non-eruptive volcanoes  [unitless]
!  (40) NNVOL      (INTEGER) : Number of non-eruptive volcanoes [unitless]
!  (41) INV        (INTEGER) : Longitude of non-erup volcanoes  [degrees]
!  (42) JNV        (INTEGER) : Latitude of non-erup volcanoes   [degrees]
!  (43) IELVn      (INTEGER) : Elevation of non-erup volcanoes  [m]
!  (44) Env        (INTEGER) : SO2 em. from non-erup volcanoes  [kg SO2/box/s]
!  (45) TCOSZ      (REAL*8 ) : Sum of cos(SZA) for offline run  [unitless] 
!  (46) TTDAY      (REAL*8 ) : Total daylight length at (I,J)   [minutes]
!  (47) SMALLNUM   (REAL*8 ) : Small number - prevent underflow [unitless]
!  
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_VCLDF         : Computes volume cloud fraction for SO2 chemistry 
!  (2 ) GET_LWC           : Computes liquid water content as a function of T
!  (3 ) CHEMSULFATE       : Driver routine for sulfate/aerosol chemistry
!  (4 ) CHEM_DMS          : Chemistry routine for DMS tracer
!  (5 ) CHEM_H2O2         : Chemistry routine for H2O2 tracer
!  (6 ) CHEM_SO2          : Chemistry routine for SO2 tracer
!  (7 ) AQCHEM_SO2        : Computes reaction rates for aqueous SO2 chemistry
!  (8 ) CHEM_SO4          : Chemistry routine for SO4 tracer
!  (9 ) CHEM_MSA          : Chemistry routine for MSA tracer
!  (10) CHEM_NH3          : Chemistry routine for ammonia tracer
!  (11) CHEM_NH4          : Chemistry routine for ammonium tracer
!  (12) CHEM_NIT          : Chemistry routine for nitrates tracer
!  (13) EMISSSULFATE      : Driver routine for sulfate/aerosol emissions
!  (14) SRCDMS            : Emission routine for DMS tracer
!  (15) SRCSO2            : Emission routine for SO2 tracer
!  (16) SRCSO4            : Emission routine for SO4 tracer
!  (17) SRCNH3            : Emission routine for NH3 tracer
!  (18) GET_OH            : Returns OH for coupled or offline simulations
!  (19) SET_OH            : Resets modified OH in SMVGEAR's CSPEC array
!  (20) GET_NO3           : Returns NO3 for coupled or offline simulations
!  (21) SET_NO3           : Resets modified OH in SMVGEAR's CSPEC array
!  (22) GET_O3            : Returns O3 for coupled or offline simulations
!  (23) READ_NONERUP_VOLC : Reads SO2 emissions from non-eruptive volcanoes
!  (24) READ_ERUP_VOLC    : Reads SO2 emissions from eruptive volcanoes 
!  (25) READ_ANTHRO_SOx   : Reads anthropogenic SO2 and SO4 emissions
!  (26) READ_OCEAN_DMS    : Reads biogenic DMS emissions from oceans
!  (27) READ_SST          : Reads monthly mean sea-surface temperatures
!  (28) READ_BIOMASS_SO2  : Reads SO2 emissions from biomass burning
!  (29) READ_AIRCRAFT_SO2 : Reads SO2 emissions from aircraft exhaust
!  (30) READ_ANTHRO_NH3   : Reads NH3 emissions from anthropogenic sources
!  (31) READ_NATURAL_NH3  : Reads NH3 emissions from natural sources
!  (32) READ_BIOMASS_NH3  : Reads NH3 biomass burning emissions
!  (33) READ_OXIDANT      : Reads monthly mean O3 and H2O2 for offline run
!  (34) OHNO3TIME         : Computes time arrays for scaling offline OH, NO3
!  (35) INIT_SULFATE      : Allocates & zeroes module arrays
!  (36) CLEANUP_SULFATE   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by sulfate_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f       : Module containing routines for binary pch file I/O
!  (2 ) comode_mod.f      : Module containing SMVGEAR allocatable arrays
!  (3 ) dao_mod.f         : Module containing DAO met field arrays
!  (4 ) diag_mod.f        : Module containing GEOS-CHEM diagnostic arrays
!  (5 ) drydep_mod.f      : Module containing GEOS-CHEM dry deposition routines
!  (6 ) error_mod.f       : Module containing NaN, other error check routines
!  (7 ) file_mod.f        : Module containing file unit numbers & error checks
!  (8 ) grid_mod.f        : Module containing horizontal grid information
!  (9 ) global_no3_mod.f  : Module containing routines to read 3-D NO3 field
!  (10) global_oh_mod.f   : Module containing routines to read 3-D OH field
!  (11) pressure_mod.f    : Module containing routines to compute P(I,J,L)
!  (12) tracerid_mod.f    : Module containing pointers to tracers & emissions
!  (13) transfer_mod.f    : Module containing routines to cast & resize arrays
!  (14) time_mod.f        : Module containing routines to compute time & date
!  (15) uvalbedo_mod.f    : Module containing UV albedo array and reader
!  (16) wetscav_mod.f     : Module containing routines for wetdep & scavenging
!
!  References
!  ============================================================================
!  (1 ) Andreae & Merlet, 2001
!
!  NOTES:
!  (1 ) All module variables are declared PRIVATE (i.e., they can only
!        be seen from within this module (bmy, 6/2/00)
!  (2 ) The routines in "sulfate_mod.f" assume that we are doing chemistry
!        over the global region (e.g. IIPAR=IGLOB, JJPAR=JGLOB). (bmy, 6/8/00)
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
!******************************************************************************
!
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "sulfate_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: SSTEMP,     DMSo,       ESO2_ac,    ESO2_an    
      PRIVATE :: ENH3_an,    ENH3_bb,    ENH3_bf,    ENH3_na
      PRIVATE :: ESO2_bb,    ESO2_ev,    ESO2_nv,    ESO4_an
      PRIVATE :: PSO2_DMS,   PMSA_DMS,   PSO4_SO2,   LSO2_AQ
      PRIVATE :: VCLDF,      NEV,        NEVOL,      IEV        
      PRIVATE :: JEV,        IDAYs,      IDAYe,      IHGHT      
      PRIVATE :: IELVe,      Eev,        NNV,        NNVOL      
      PRIVATE :: INV,        JNV,        IELVn,      Env
      PRIVATE :: XNUMOL_OH,  XNUMOL_O3,  XNUMOL_NO3, JH2O2
      PRIVATE :: DRYSO2,     DRYSO4,     DRYMSA,     DRYNH3
      PRIVATE :: DRYNH4,     DRYNIT,     DRYH2O2,    TCOSZ
      PRIVATE :: TTDAY,      ESO2_bf,    SMALLNUM

      ! PRIVATE module routines
      PRIVATE :: GET_VCLDF,         GET_LWC,           CHEM_DMS 
      PRIVATE :: CHEM_SO2,          CHEM_SO4,          CHEM_MSA    
      PRIVATE :: CHEM_NH3,          CHEM_NH4,          CHEM_NIT    
      PRIVATE :: SRCDMS,            SRCSO2,            SRCSO4
      PRIVATE :: SRCNH3,            READ_NONERUP_VOLC, READ_ERUP_VOLC
      PRIVATE :: READ_ANTHRO_SOx,   READ_OCEAN_DMS,    READ_SST
      PRIVATE :: READ_AIRCRAFT_SO2, READ_ANTHRO_NH3,   READ_NATURAL_NH3
      PRIVATE :: READ_BIOMASS_NH3,  AQCHEM_SO2,        GET_O3            
      PRIVATE :: GET_OH,            SET_OH,            READ_OXIDANT
      PRIVATE :: OHNO3TIME

      !=================================================================
      ! MODULE VARIABLES (see descriptions listed above)
      !=================================================================

      ! Time variable
      INTEGER              :: ELAPSED_SEC

      ! Logical Flags
      LOGICAL, PARAMETER   :: LENV = .TRUE.
      LOGICAL, PARAMETER   :: LEEV = .TRUE.
      
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
      REAL*8,  ALLOCATABLE :: ESO4_an(:,:,:) 
      REAL*8,  ALLOCATABLE :: IJSURF(:,:)
      REAL*8,  ALLOCATABLE :: JH2O2(:,:,:)
      REAL*8,  ALLOCATABLE :: LSO2_AQ(:,:,:)
      REAL*8,  ALLOCATABLE :: O3m(:,:,:)
      REAL*8,  ALLOCATABLE :: PH2O2m(:,:,:)
      REAL*8,  ALLOCATABLE :: PMSA_DMS(:,:,:)
      REAL*8,  ALLOCATABLE :: PSO2_DMS(:,:,:)
      REAL*8,  ALLOCATABLE :: PSO4_SO2(:,:,:)
      REAL*8,  ALLOCATABLE :: SOx_SCALE(:,:)
      REAL*8,  ALLOCATABLE :: SSTEMP(:,:)
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: VCLDF(:,:,:)


      ! Eruptive volcanoes
      INTEGER, PARAMETER   :: NEV=50
      INTEGER              :: NEVOL
      INTEGER, ALLOCATABLE :: IEV(:),   JEV(:)
      INTEGER, ALLOCATABLE :: IDAYs(:), IDAYe(:)
      INTEGER, ALLOCATABLE :: IHGHT(:), IELVe(:)
      REAL*8,  ALLOCATABLE :: EEV(:)

      ! Non-eruptive volcanoes 
      INTEGER, PARAMETER   :: NNV=50
      INTEGER              :: NNVOL
      INTEGER, ALLOCATABLE :: INV(:), JNV(:), IELVn(:)
      REAL*8,  ALLOCATABLE :: ENV(:)
      
      ! Pointers to drydep species w/in DEPSAV
      INTEGER              :: DRYSO2, DRYSO4, DRYMSA
      INTEGER              :: DRYNH3, DRYNH4, DRYNIT, DRYH2O2

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
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : RH
      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER              :: I,    J,    L
      REAL*8               :: PRES, PSFC, RH2, R0, B0

      ! Parameters
      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
		
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================
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
!  Function GET_LWC returns the cloud liquid water content at a GEOS-CHEM
!  grid box as a function of temperature. (rjp, bmy, 10/31/02, 1/14/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T (REAL*8) : Temperature value at a GEOS-CHEM grid box [K]
!
!  NOTES:
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
!  (rjp, bdf, bmy, 5/31/00, 3/27/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD,     AIRDEN,  CLDF, 
     &                           SUNCOS, MAKE_RH, CONVERT_UNITS
      USE DRYDEP_MOD,     ONLY : DEPSAV, DEPNAME, NUMDEP
      USE ERROR_MOD,      ONLY : DEBUG_MSG, CHECK_STT
      USE GLOBAL_OH_MOD,  ONLY : GET_GLOBAL_OH
      USE GLOBAL_NO3_MOD, ONLY : GET_GLOBAL_NO3
      USE TIME_MOD,       ONLY : GET_MONTH, GET_TS_CHEM, GET_ELAPSED_SEC
      USE TRACERID_MOD
 
#     include "CMN_SIZE"     ! Size parameters 
#     include "CMN"          ! STT, TCVV, NSRCX

      ! Local variables
      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE       :: LASTMONTH = -99
      INTEGER             :: I, J, L, N
      REAL*8              :: DTCHEM

      ! External functions   
      REAL*8,  EXTERNAL   :: BOXVL

      !=================================================================
      ! CHEMSULFATE begins here!
      !=================================================================

      ! Establish indices w/in DEPSAV array
      IF ( FIRSTCHEM ) THEN

         ! Initialize arrays (if not already done before)
         CALL INIT_SULFATE

         ! Find drydep species in DEPSAV
         DO N = 1, NUMDEP
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE ( 'H2O2' )
                  DRYH2O2 = N
               CASE ( 'SO2' )
                  DRYSO2 = N
               CASE ( 'SO4' )
                  DRYSO4 = N
               CASE ( 'MSA' )
                  DRYMSA = N
               CASE ( 'NH3' )
                  DRYNH3 = N
               CASE ( 'NH4' )
                  DRYNH4 = N
               CASE ( 'NIT' )
                  DRYNIT = N
               CASE DEFAULT
                  ! Nothing
            END SELECT        
         ENDDO
         
         ! Reset first-time flag
         FIRSTCHEM = .FALSE.
      ENDIF

      ! Read monthly mean fields for offline runs 
      IF ( NSRCX == 10 .and. GET_MONTH() /= LASTMONTH ) THEN
         CALL GET_GLOBAL_OH( GET_MONTH() )
         CALL GET_GLOBAL_NO3( GET_MONTH() )
         LASTMONTH = GET_MONTH()
      ENDIF

      ! Compute time scaling arrays for offline OH, NO3
      IF ( NSRCX == 10 ) CALL OHNO3TIME

      ! Store NTIME in a shadow variable
      ELAPSED_SEC = GET_ELAPSED_SEC()

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize module arrays
      PSO2_DMS = 0d0
      PMSA_DMS = 0d0
      PSO4_SO2 = 0d0
                  
      !================================================================= 
      ! Call individual chemistry routines for sulfate/aerosol tracers
      !=================================================================

      ! Convert STT from [kg] -> [v/v] 
      CALL CONVERT_UNITS( 1, NTRACE, TCVV(1:NTRACE), 
     &                    AD,  STT(:,:,:,1:NTRACE) )

      ! DMS (offline only)
      IF ( NSRCX == 10 ) THEN
         CALL CHEM_DMS
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_DMS' ) 
      ENDIF

      ! H2O2 (offline only)
      IF ( NSRCX == 10 ) THEN
         CALL CHEM_H2O2
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: CHEM_H2O2' )
      ENDIF

      ! SO2 
      CALL GET_VCLDF
      CALL CHEM_SO2
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_SO2' )

      ! SO4 
      CALL CHEM_SO4
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_SO4' )

      ! MSA 
      CALL CHEM_MSA
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_MSA' )

      ! NH3 
      CALL CHEM_NH3
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_NH3' )

      ! NH4 
      CALL CHEM_NH4
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_NH4' )

      ! SULFUR NITRATE 
      CALL CHEM_NIT
      IF ( LPRT ) CALL DEBUG_MSG( '### CHEMSULFATE: after CHEM_NIT' )
      
      ! Convert STT from [v/v] -> [kg]
      CALL CONVERT_UNITS( 2, NTRACE, TCVV(1:NTRACE), 
     &                    AD, STT(:,:,:,1:NTRACE) )

      ! We have already gone thru one chemistry iteration
      FIRSTCHEM = .FALSE. 
         
      ! Return to calling program
      END SUBROUTINE CHEMSULFATE

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_DMS
!
!******************************************************************************
!  Subroutine CHEM_DMS is the DMS chemistry subroutine from Mian Chin's    
!  GOCART model, modified for use with the GEOS-CHEM model.
!  (rjp, bdf, bmy, 5/31/00, 3/27/03)  
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
!******************************************************************************
!
      ! Reference to F90 modules
      USE DAO_MOD,      ONLY : AD, AIRDEN, SUNCOS, T
      USE DIAG_MOD,     ONLY : AD05
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTDMS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NSRCX, LPAUSE
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_GCTM"     ! AIRMW
#     include "CMN_DIAG"     ! ND05, LD05

      ! Local variables
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
         IF ( L >= LPAUSE(I,J) ) CYCLE

         ! IJLOOP is the 1-D grid box index for SUNCOS
         IJLOOP = IJSURF(I,J)

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

            RK2 = 1.2d-11 * EXP( -260.d0 / TK ) * OH 
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
         !###DMS    = MAX( DMS, 1d-32 )
         IF ( DMS < SMALLNUM ) DMS = 0d0

         ! Archive initial OH and NO3 for diagnostics
         OH0    = OH
         XNO30  = XNO3

         IF ( NSRCX == 3 ) THEN

            ! Update OH after rxn w/ DMS (coupled runs only)
            OH    = OH0 - ( ( DMS0 - DMS_OH ) * AIRDEN(L,I,J) * f )
            !###OH    = MAX( OH, 1D-32 )
            IF ( OH < SMALLNUM ) OH = 0d0

            ! Update NO3 after rxn w/ DMS (coupled runs only)
            XNO3  = XNO30 - ( ( DMS_OH - DMS ) * AIRDEN(L,I,J) * f )
            !###XNO3  = MAX( XNO3, 1D-32 )
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
            XOH  = ( DMS0   - DMS_OH ) / Fx * AD(I,J,L) / TCVV_S
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
            AD05(I,J,L,8) = AD05(I,J,L,8) + LOH

            ! Store L(NO3) by DMS [kg NO3/timestep]
            AD05(I,J,L,9) = AD05(I,J,L,9) + LNO3

         ENDIF

         !==============================================================
         ! For a coupled fullchem/aerosol run, save OH [molec/cm3] 
         ! and NO3 [molec/cm3] back into the CSPEC array of SMVGEAR
         !==============================================================
         IF ( NSRCX == 3 ) THEN
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
!  the SMVGEAR module. (rjp, bmy, 11/26/02, 8/1/03)
!                                                                           
!  NOTES:
!  (1 ) Bug fix: need to multiply DXYP by 1d4 for cm2 (bmy, 3/23/03)
!  (2 ) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f"
!        Now use functions GET_MONTH and GET_TS_CHEM from "time_mod.f".
!        (bmy, 3/27/03)
!  (3 ) Now references PBLFRAC from "drydep_mod.f".  Now apply dry deposition 
!        throughout the entire PBL.  Added FREQ variable. (bmy, 8/1/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,      ONLY : AD, AIRDEN, OPTD, SUNCOS, T
      USE DIAG_MOD,     ONLY : AD44 
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTH2O2
      USE TRANSFER_MOD, ONLY : TRANSFER_3D_TROP
      USE UVALBEDO_MOD, ONLY : UVALBEDO
 
#     include "cmn_fj.h"     ! IPAR, JPAR, LPAR, CMN_SIZE
#     include "CMN"          ! TCVV, LPAUSE
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_GCTM"     ! AIRMW
#     include "CMN_SETUP"    ! DATA_DIR
#     include "CMN_DIAG"     ! ND44
      
      ! Local variables
      LOGICAL                :: FIRST     = .TRUE.
      INTEGER, SAVE          :: LASTMONTH = -99
      INTEGER                :: I, J, L
      REAL*4                 :: ARRAY(IGLOB,JGLOB,LLTROP)
      REAL*8                 :: DT,   Koh,   DH2O2, M,     F   
      REAL*8                 :: XTAU, H2O20, H2O2,  ALPHA, FLUX, FREQ
      REAL*8,  PARAMETER     :: A = 2.9d-12
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! CHEM_H2O2 begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. DRYH2O2 == 0 ) RETURN 

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
      F  = 1000.d0 / AIRMW * 6.022d23 * 1.d-6

      !=================================================================
      ! For offline run: read J(H2O2) from disk below
      !=================================================================
      IF ( GET_MONTH() /= LASTMONTH ) THEN 

         ! File name to read data 
         FILENAME = TRIM( DATA_DIR )           // 
     &              'sulfate_sim_200210/JH2O2.'// GET_NAME_EXT() //
     &              '.'                        // GET_RES_EXT()
           
         ! Print filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - CHEM_H2O2: Reading ', a )

         ! Get TAU0 value for 1998
         XTAU = GET_TAU0( GET_MONTH(), 1, 1998 )
	
         ! Read J(H2O2) [s-1]  from disk (only up to tropopause)
         CALL READ_BPCH2( FILENAME, 'JV-MAP-$', 3,      XTAU, 
     &                    IGLOB,    JGLOB,      LLTROP, ARRAY )

         ! Cast to REAL*8 and resize if necessary
         CALL TRANSFER_3D_TROP( ARRAY, JH2O2 )
            
         ! Reset LASTMONTH
         LASTMONTH = GET_MONTH()
      ENDIF

      !=================================================================
      ! Loop over tropopsheric grid boxes and do chemistry
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, M, H2O20, KOH, FREQ, ALPHA, DH2O2, H2O2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L  = 1, LLTROP
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR

         ! Initialize for safety's sake 
         FLUX = 0d0
         FREQ = 0d0

         ! Skip stratospheric boxes
         IF ( L >= LPAUSE(I,J) ) CYCLE

         ! Density of air [molec/cm3]
         M     = AIRDEN(L,I,J) * f  

         ! Initial H2O2 [v/v]
         H2O20 = STT(I,J,L,IDTH2O2)

         ! Loss frequenty due to OH oxidation [s-1]
         KOH   = A * EXP( -160.d0 / T(I,J,L) ) * GET_OH(I,J,L)  

         ! H2O2 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
         ! of the grid box (I,J,L) that is located beneath the PBL top
         FREQ  = DEPSAV(I,J,DRYH2O2) * PBLFRAC(I,J,L)

!------------------------------------------------------------------------------
! Prior to 8/1/03:
! Now multiply DEPSAV by PBLFRAC, in order to do dry deposition below the 
! PBLTOP.  PBLFRAC is the fraction of each level below the PBL top.
! (rjp, bmy, 8/1/03)
!         ! Compute loss fraction [unitless].  
!         ! At surface, also include the loss by dry deposition.
!         IF ( L == 1 ) THEN 
!            ALPHA = 1.D0 + ( KOH + JH2O2(I,J,L) + 
!     &                       DEPSAV(I,J,DRYH2O2) ) * DT 
!         ELSE 
!            ALPHA = 1.D0 + ( KOH + JH2O2(I,J,L) ) * DT               
!         ENDIF
!------------------------------------------------------------------------------

         ! Compute loss fraction from OH, photolysis, drydep [unitless].  
         ALPHA = 1.D0 + ( KOH + JH2O2(I,J,L) + FREQ ) * DT 

         ! Delta H2O2 [v/v]
         DH2O2 = PH2O2m(I,J,L) * DT / ( ALPHA * M )
         
         ! Final H2O2 [v/v]
         !###H2O2  = MAX( ( H2O20 / ALPHA + DH2O2 ), 1.D-32 )
         H2O2  = ( H2O20 / ALPHA + DH2O2 )
         IF ( H2O2 < SMALLNUM ) H2O2 = 0d0

         ! Store final H2O2 in STT
         STT(I,J,L,IDTH2O2) = H2O2

         !==============================================================
         ! ND44 diagnostics (make sure this is correct!)
         !==============================================================
         !--------------------------------------------------------------
         ! Prior to 8/1/03:
         ! Now apply drydep to entire PBL (bmy, 8/1/03)
         !IF ( ND44 > 0 .and. L == 1 ) THEN
         !--------------------------------------------------------------
         IF ( ND44 > 0 .AND. FREQ > 0d0 ) THEN

            ! Convert H2O2 from [v/v] to H2O2 [molec/cm2]
            !-----------------------------------------------------------
            ! Prior to 8/1/03:
            ! AD(I,J,1) needs to be AD(I,J,L) (bmy, 8/1/03)
            !FLUX = ( H2O20 - H2O2 ) * AD(I,J,1) / TCVV(IDTH2O2)
            !-----------------------------------------------------------
            FLUX = ( H2O20 - H2O2 ) * AD(I,J,L) / TCVV(IDTH2O2)
            FLUX = FLUX * XNUMOL(IDTH2O2) / GET_AREA_CM2( J )

            ! Multiply by drydep freq [s-1] to convert to [molec/cm2/s]
            !------------------------------------------------------------
            ! Prior to 8/1/03:
            ! Need to multiply DEPSAV by PBLFRAC (bmy, 8/1/03)
            !FLUX = FLUX * DEPSAV(I,J,DRYH2O2) 
            !------------------------------------------------------------
            FLUX = FLUX * FREQ

            ! Save drydep loss in AD44 [molec/cm2/s]
            AD44(I,J,DRYH2O2,1) = AD44(I,J,DRYH2O2,1) + FLUX 
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE CHEM_H2O2

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_SO2
!
!******************************************************************************
!  Subroutine CHEM_SO2 is the SO2 chemistry subroutine 
!  (rjp, bmy, 11/26/02, 8/1/03) 
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
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,      ONLY : AD,      AIRDEN, T
      USE DIAG_MOD,     ONLY : AD05,    AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV,  PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTH2O2, IDTSO2
      USE WETSCAV_MOD,  ONLY : H2O2s,   SO2s

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN"         ! NCHEM, STT, LPAUSE
#     include "CMN_GCTM"    ! AIRMW
#     include "CMN_O3"      ! XNUMOL
#     include "CMN_DIAG"    ! LD05, ND05, ND44
#     include "CMN_SETUP"   ! DATA_DIR

      ! Local variables
      INTEGER               :: I,      J,       L,      I1,   I2
      INTEGER               :: II,     NSTEP
      REAL*8                :: K0,     Ki,      KK,     M,    L1
      REAL*8                :: L2,     L3,      Ld,     F,    Fc
      REAL*8                :: RK,     RKT,     DTCHEM, DT_T, TK
      REAL*8                :: F1,     RK1,     RK2,    RK3,  SO20
      REAL*8                :: SO2_cd, H2O20,   O3,     L2S,  L3S
      REAL*8                :: LWC,    KaqH2O2, KaqO3,  PATM, FLUX
      REAL*8                :: AREA_CM2

      ! Parameters
      REAL*8,  PARAMETER    :: HPLUS  = 3.16227766016837953d-5  !pH = 4.5
      REAL*8,  PARAMETER    :: MINDAT = 1.d-20

      !=================================================================
      ! CHEM_SO2 begins here!
      !=================================================================
      IF ( IDTH2O2 == 0 .or. IDTSO2 == 0 .or. DRYSO2 == 0 ) RETURN

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
      F      = 1000.d0 / AIRMW * 6.022d23 * 1.d-6
      Ki     = 1.5d-12

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO20, H2O20, O3, PATM, TK, K0, M, KK, F1, RK1  )
!$OMP+PRIVATE( RK2, RK, RKT, SO2_cd, L1, Ld, L2, L2S, L3, L3S, FC, LWC )
!$OMP+PRIVATE( KaqH2O2, KaqO3, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP  
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake 
         AREA_CM2 = 0d0
         FLUX     = 0d0
         Ld       = 0d0

         ! Skip stratospheric boxes
         IF ( L >= LPAUSE(I,J) ) CYCLE

         ! Initial SO2, H2O2 and O3 [v/v]
         SO20   = STT(I,J,L,IDTSO2)         
         H2O20  = STT(I,J,L,IDTH2O2)
         O3     = GET_O3(I,J,L)

         ! PATM  : Atmospheric pressure in atm
         PATM   = GET_PCENTER( I, J, L ) / 1013.25d0

         ! TK : Temperature [K]
         TK     = T(I,J,L)

         IF ( NSRCX == 10 ) THEN

            ! Gas phase SO4 production is done here in offline run only 
            ! RK1: SO2 + OH(g) [s-1]  (rjp, bmy, 3/23/03)
            K0  = 3.0d-31 * ( 300.d0 / TK )**3.3d0
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

         ! SO2 drydep frequency [s-1]  
         !-------------------------------------------------------------
         ! Prior to 8/1/03:
         ! Multiply DEPSAV by PBLFRAC to account for the fraction of
         ! each grid box (I,J,L) beneath the PBL top. (bmy, 8/1/03)
         !IF ( L == 1 ) THEN
         !   RK2 = DEPSAV(I,J,DRYSO2)
         !ELSE
         !   RK2 = 0.D0
         !ENDIF
         !-------------------------------------------------------------

         ! SO2 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
         ! of grid box (I,J,L) that is located beneath the PBL top
         RK2    = DEPSAV(I,J,DRYSO2) * PBLFRAC(I,J,L)

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
             
            !-----------------------------------------------------------
            ! Prior to 8/1/03:
            ! We no longer need to restrict drydep to the first level
            ! (rjp, bmy, 8/1/03)
            !IF ( L == 1 ) THEN
            !   Ld  = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK2/RK
            !ELSE
            !   Ld  = 0.d0
            !ENDIF
            !-----------------------------------------------------------
            Ld  = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK2/RK
            
         ELSE
            SO2_cd = SO20
            L1     = 0.d0
         ENDIF

         !==============================================================
         ! Update SO2 concentration after cloud chemistry          
         ! SO2 chemical loss rate = SO4 production rate [v/v/timestep]
         !==============================================================
      
         ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
         FC      = VCLDF(I,J,L)

         ! Liquid water content in cloudy area of grid box [m3/m3]
         LWC     = GET_LWC( TK ) * FC

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
     &        ( SO2_cd > MINDAT )  .AND. 
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
            CALL AQCHEM_SO2( LWC, TK,    PATM,    SO2_cd, H2O20, 
     &                       O3,  HPLUS, KaqH2O2, KaqO3 ) 

            ! Aqueous phase SO2 loss rate (v/v/timestep): 
            L2  = EXP( ( SO2_cd - H2O20 ) * KaqH2O2 * DTCHEM )  
            L3  = EXP( ( SO2_cd - O3    ) * KaqO3   * DTCHEM )       

            ! Loss by H2O2
            L2S = SO2_cd * H2O20 * (L2 - 1.D0) / ((SO2_cd * L2) - H2O20)  

            ! Loss by O3
            L3S = SO2_cd * O3    * (L3 - 1.D0) / ((SO2_cd * L3) - O3)     
          
            SO2_cd = MAX( SO2_cd - ( L2S + L3S ), MINDAT )
            H2O20  = MAX( H2O20  - L2S,           MINDAT )

            ! Update SO2 level, save SO2[ppv], H2O2[ppv] for WETDEP
            SO2s( I,J,L) = SO2_cd
            H2O2s(I,J,L) = H2O20

         ELSE

            ! Otherwise, don't do aqueous chemistry, and
            ! save the original concentrations into SO2 and H2O2
            H2O2s(I,J,L) = MAX( H2O20,  1.0d-32 )
            SO2s(I,J,L ) = MAX( SO2_cd, 1.0d-32 )
            L2S          = 0.d0
            L3S          = 0.d0

         ENDIF

         ! Store updated SO2, H2O2 back to the tracer arrays 
         STT(I,J,L,IDTSO2)  = SO2s( I,J,L)
         STT(I,J,L,IDTH2O2) = H2O2s(I,J,L)

         ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
         PSO4_SO2(I,J,L) = L1 + L2S + L3S

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
         ENDIF

         !=================================================================
         ! ND44 Diagnostic: Drydep flux of SO2 [molec/cm2/s]
         !=================================================================
         !------------------------------------------------------------
         ! Prior to 8/1/03:
         ! We no longer need to restrict drydep to the first layer
         ! (rjp, bmy, 8/1/03)
         !IF ( ND44 > 0 .AND. L == 1 ) THEN
         !------------------------------------------------------------
         IF ( ND44 > 0 .AND. Ld > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert [v/v/timestep] to [molec/cm2/s]
            !--------------------------------------------------------
            ! Prior to 8/1/03:
            ! AD(I,J,1) is now AD(I,J,L) (bmy, 8/1/03)
            !FLUX = Ld   * AD(I,J,1)      / TCVV(IDTSO2)
            !--------------------------------------------------------
            FLUX = Ld   * AD(I,J,L)      / TCVV(IDTSO2)
            FLUX = FLUX * XNUMOL(IDTSO2) / AREA_CM2 / DTCHEM

            ! Store in AD44
            AD44(I,J,DRYSO2,1) = AD44(I,J,DRYSO2,1) + FLUX
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_SO2

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
!  model, modified for the GEOS-CHEM model. (rjp, bdf, bmy, 5/31/00, 8/1/03) 
!                                                                          
!  Module Variables Used:
!  ============================================================================
!  (1 ) PSO4_SO2 (REAL*8 ) : Array for P(SO4) from SO2 [v/v/timestep]
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTSO4

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER                :: I, J, L
      REAL*8                 :: SO4, SO40, RKT, DTCHEM, FLUX, AREA_CM2

      !=================================================================
      ! CHEM_SO4 begins here!
      !=================================================================
      IF ( IDTSO4 == 0 .or. DRYSO4 == 0 ) RETURN

      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, SO40, RKT, SO4, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize for safety's sake
         AREA_CM2 = 0d0
         RKT      = 0d0
         FLUX     = 0d0

         ! Skip stratospheric boxes
         IF ( L >= LPAUSE(I,J) ) CYCLE

         ! Initial SO4 [v/v]
         SO40 = STT(I,J,L,IDTSO4)
         
         !==============================================================
         ! At surface, we have L(SO4) due to drydep and P(SO4) from SO2
         ! At altitude, we only have P(SO4) from SO2
         !==============================================================
!-----------------------------------------------------------------------------
! Prior to 8/1/03:
! Now apply drydep throughout the entire PBL (bmy, 8/1/03)
!         IF ( L == 1 ) THEN
!
!            ! RKT: Fraction of SO4 lost to drydep [unitless]
!            RKT = DEPSAV(I,J,DRYSO4) * DTCHEM
!            
!            ! SO4 concentration [v/v]
!            SO4 = ( SO40 * EXP( -RKT )                           ) + 
!     &            ( PSO4_SO2(I,J,L)/RKT * ( 1.d0 - EXP( -RKT ) ) )
!         ELSE
!
!            ! SO4 production from SO2 [v/v/timestep]
!            SO4 = SO40 + PSO4_SO2(I,J,L)
!
!         ENDIF
!-----------------------------------------------------------------------------

         ! SO4 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
         ! of each vertical level that is located below the PBL top
         RKT = DEPSAV(I,J,DRYSO4) * PBLFRAC(I,J,L)

         ! RKT > 0 denotes that we have drydep occurring
         IF ( RKT > 0d0 ) THEN
            
            ! Fraction of SO4 lost to drydep [unitless]
            RKT = RKT * DTCHEM

            ! SO4 concentration [v/v]
            SO4 = ( SO40 * EXP( -RKT )                           ) + 
     &            ( PSO4_SO2(I,J,L)/RKT * ( 1.d0 - EXP( -RKT ) ) )
         ELSE

            ! SO4 production from SO2 [v/v/timestep]
            SO4 = SO40 + PSO4_SO2(I,J,L)

         ENDIF

         ! Final SO4 [v/v]
         !###SO4               = MAX( SO4, 1.0d-32 )
         IF ( SO4 < SMALLNUM ) SO4 = 0d0
         STT(I,J,L,IDTSO4) = SO4

         !==============================================================
         ! ND44 Diagnostic: Drydep flux of SO4 [molec/cm2/s]
         !==============================================================
         !--------------------------------------------------------------
         ! Prior to 8/1/03:
         ! We no longer need to restrict drydep to the first layer
         ! (rjp, bmy, 8/1/03)
         !IF ( ND44 > 0 .and. L == 1 ) THEN
         !--------------------------------------------------------------
         IF ( ND44 > 0 .AND. RKT > 0d0 ) THEN

            ! Surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert from [v/v/timestep] to [molec/cm2/s]
            FLUX = SO40 - SO4 + PSO4_SO2(I,J,L) 
            FLUX = FLUX * AD(I,J,L)      / TCVV(IDTSO4)
            FLUX = FLUX * XNUMOL(IDTSO4) / AREA_CM2 / DTCHEM

            ! Store in AD44
            AD44(I,J,DRYSO4,1) = AD44(I,J,DRYSO4,1) + FLUX
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_SO4

!-----------------------------------------------------------------------------

      SUBROUTINE CHEM_MSA
!
!******************************************************************************
!  Subroutine CHEM_MSA is the SO4 chemistry subroutine from Mian Chin's GOCART
!  model, modified for the GEOS-CHEM model. (rjp, bdf, bmy, 5/31/00, 8/1/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTMSA

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! NCHEM
#     include "CMN_GCTM"     ! AIRMW
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER  :: I,      J,    L
      REAL*8   :: DTCHEM, MSA0, MSA, RK, RKT, FLUX, AREA_CM2

      !=================================================================
      ! CHEM_MSA begins here!
      !=================================================================
      IF ( IDTMSA == 0 .or. DRYMSA == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      ! Loop over tropospheric grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, MSA0, RKT, MSA, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP 
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only apply drydep loss to boxes w/in the PBL
         IF ( PBLFRAC(I,J,L) > 0 ) THEN
         
            ! Initial MSA [v/v]
            MSA0 = STT(I,J,L,IDTMSA) 

!------------------------------------------------------------------------------
! Prior to 8/1/03:
! Now apply drydep throughout the entire PBL (bmy, 8/1/03)
!         !==============================================================
!         ! At surface, we have L(MSA) due to drydep and P(MSA) from DMS
!         ! At altitude, we only have P(MSA) from DMS
!         !==============================================================
!         IF ( L == 1 ) THEN
!
!            ! RKT: Fraction of MSA lost to drydep [unitless]
!            RKT = DEPSAV(I,J,DRYMSA) * DTCHEM
!
!            ! Modified MSA concentration 
!            MSA = ( MSA0 * EXP( -RKT )                        ) +
!     &            ( PMSA_DMS(I,J,L)/RKT * ( 1d0 - EXP( -RKT ) ) )
!
!         ELSE
!
!            ! MSA production from DMS [v/v/timestep]
!            MSA = MSA0 + PMSA_DMS(I,J,L)
!
!         ENDIF
!------------------------------------------------------------------------------



            ! MSA drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            RKT = DEPSAV(I,J,DRYMSA) * PBLFRAC(I,J,L)

            ! RKT > 0 denotes that we have drydep occurring
            IF ( RKT > 0.d0 ) THEN

               ! Fraction of MSA lost to drydep [unitless]
               RKT = RKT * DTCHEM

               ! Modified MSA concentration 
               MSA = ( MSA0 * EXP( -RKT )                        ) +
     &               ( PMSA_DMS(I,J,L)/RKT * ( 1d0 - EXP( -RKT ) ) )

            ELSE

               ! MSA production from DMS [v/v/timestep]
               MSA = MSA0 + PMSA_DMS(I,J,L)

            ENDIF

            ! Final MSA [v/v]
            !###MSA               = MAX( MSA, 1.0d-32 )
            IF ( MSA < SMALLNUM ) MSA = 0d0
            STT(I,J,L,IDTMSA) = MSA

            !===========================================================
            ! ND44 Diagnostic: Drydep flux of MSA [molec/cm2/s]
            !===========================================================
            !-----------------------------------------------------------
            ! Prior to 8/1/03:
            ! We are no longer restricted to doing drydep in the first 
            ! level. (rjp, bmy, 8/1/03)
            !IF ( ND44 > 0 .and. L == 1 ) THEN
            !-----------------------------------------------------------
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Convert [v/v/timestep] to [molec/cm2/s]
               FLUX = MSA0 - MSA + PMSA_DMS(I,J,L)                    
               FLUX = FLUX * AD(I,J,L)      / TCVV(IDTMSA)            
               FLUX = FLUX * XNUMOL(IDTMSA) / AREA_CM2 / DTCHEM    
               
               ! Store in AD44
               AD44(I,J,DRYMSA,1) = AD44(I,J,DRYMSA,1) + FLUX
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_MSA

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH3
!
!******************************************************************************
!  Subroutine CHEM_NH3 removes NH3 from the surface via dry deposition.
!  (rjp, bdf, bmy, 1/2/02, 8/1/03)  
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTNH3

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! TCVV, LPAUSE
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER :: I,      J,        L
      REAL*8  :: DTCHEM, NH3,      NH30
      REAL*8  :: FREQ,   AREA_CM2, FLUX

      !=================================================================
      ! CHEM_NH3 begins here!
      !=================================================================
      IF ( IDTNH3 == 0 .or. DRYNH3 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

!-----------------------------------------------------------------------------
! Prior to 8/1/03:
! We now need to apply drydep losses throughout the PBL (rjp, bmy, 8/1/03)
!      ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, NH30, NH3, AREA_CM2, FLUX )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Initial NH3 [v/v]
!         NH30 = STT(I,J,1,IDTNH3)
!         
!         NH3  = NH30 * EXP( -DEPSAV(I,J,DRYNH3) * DTCHEM )
!         NH3  = MAX( NH3, 1d-32 )
!
!         !==============================================================
!         ! ND44 diagnostic: Drydep flux of NH3 [molec/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 ) THEN
!
!            ! Surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!
!            ! Convert [v/v/timestep] to [molec/cm2/s]
!            FLUX = ( NH30 - NH3 ) * AD(I,J,1) / TCVV(IDTNH3)
!            FLUX = FLUX * XNUMOL(IDTNH3) / AREA_CM2 / DTCHEM
!
!            ! Store in AD44
!            AD44(I,J,DRYNH3,1) = AD44(I,J,DRYNH3,1) + FLUX
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!-----------------------------------------------------------------------------

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FREQ, NH30, NH3, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only apply drydep to boxes w/in the PBL
         IF ( PBLFRAC(I,J,L) > 0d0 ) THEN

            ! NH3 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH3) * PBLFRAC(I,J,L)

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
               IF ( ND44 > 0 ) THEN

                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )

                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH3  * AD(I,J,L)      / TCVV(IDTNH3)
                  FLUX = FLUX * XNUMOL(IDTNH3) / AREA_CM2 / DTCHEM

                  ! Store in AD44
                  AD44(I,J,DRYNH3,1) = AD44(I,J,DRYNH3,1) + FLUX
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_NH3

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NH4
!
!******************************************************************************
!  Subroutine CHEM_NH4 removes NH4 from the surface via dry deposition.
!  (rjp, bdf, bmy, 1/2/02, 8/1/03)  
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTNH4

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, LPAUSE
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER                :: I,      J,    L
      REAL*8                 :: DTCHEM, NH4,  NH40
      REAL*8                 :: FREQ,   FLUX, AREA_CM2

      !=================================================================
      ! CHEM_NH4 begins here!
      !=================================================================
      IF ( IDTNH4 == 0 .or. DRYNH4 == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

!------------------------------------------------------------------------------
! Prior to 8/1/03:
! We now need to apply drydep losses throughout the PBL (rjp, bmy, 8/1/03)
!      ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, NH40, NH4, AREA_CM2, FLUX )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Initial NH4 [v/v]
!         NH40 = STT(I,J,1,IDTNH4)
!         
!         ! Amount of NH4 left after drydep [v/v]
!         NH4 = NH40 * EXP( -DEPSAV(I,J,DRYNH4) * DTCHEM )
!         NH4 = MAX( NH4, 1.0d-32 )
!
!         ! Final NH4 [v/v]
!         STT(I,J,1,IDTNH4) = NH4
!
!         !==============================================================
!         ! ND44 diagnostic: Drydep flux of NH4 [molec/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 ) THEN
!         
!            ! Surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!
!            ! Convert from [v/v/timestep] to [molec/cm2/s]
!            FLUX = ( NH40 - NH4 ) * AD(I,J,1) / TCVV(IDTNH4)
!            FLUX = FLUX * XNUMOL(IDTNH4) / AREA_CM2 / DTCHEM
!
!            ! Store in AD44
!            AD44(I,J,DRYNH4,1) = AD44(I,J,DRYNH4,1) + FLUX
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!------------------------------------------------------------------------------

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FREQ, NH40, NH4, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only apply drydep to boxes w/in the PBL
         IF ( PBLFRAC(I,J,L) > 0d0 ) THEN

            ! NH4 drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each grid box (I,J,L) that is located beneath the PBL top
            FREQ = DEPSAV(I,J,DRYNH4) * PBLFRAC(I,J,L)

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
               IF ( ND44 > 0 ) THEN
         
                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )

                  ! Convert drydep loss from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NH4  * AD(I,J,L)      / TCVV(IDTNH4)
                  FLUX = FLUX * XNUMOL(IDTNH4) / AREA_CM2 / DTCHEM

                  ! Store in AD44
                  AD44(I,J,DRYNH4,1) = AD44(I,J,DRYNH4,1) + FLUX
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_NH4

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_NIT
!
!******************************************************************************
!  Subroutine CHEM_NIT removes SULFUR NITRATES (NIT) from the surface 
!  via dry deposition. (rjp, bdf, bmy, 1/2/02, 8/1/03)  
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV, PBLFRAC
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACERID_MOD, ONLY : IDTNIT

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, LPAUSE
#     include "CMN_GCTM"     ! AIRMW
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND44

      ! Local variables
      INTEGER :: I,      J,   L
      REAL*8  :: DTCHEM, NIT, NIT0, FREQ, AREA_CM2, FLUX

      !=================================================================
      ! CHEM_NIT begins here!
      !=================================================================
      IF ( IDTNIT == 0 .or. DRYNIT == 0 ) RETURN

      ! DTCHEM is the chemistry interval in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

!------------------------------------------------------------------------------
! Prior to 8/1/03:
! We now have to apply drydep throughout the PBL (rjp, bmy, 8/1/03)
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, NIT0, NIT, AREA_CM2, FLUX )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Initial NITRATE [v/v]
!         NIT0 = STT(I,J,1,IDTNIT)
!         
!         ! Amount of NITRATE left after dry drydep [v/v]
!         NIT  = NIT0 * EXP( -DEPSAV(I,J,DRYNIT) * DTCHEM )
!         NIT  = MAX( NIT, 1.0d-32 )
!         
!         ! Final NITRATE [v/v]
!         STT(I,J,1,IDTNIT) = NIT
!
!         !==============================================================
!         ! ND44 Diagnostic: Drydep flux of NIT [molec/cm2/s]
!         !==============================================================
!         IF ( ND44 > 0 ) THEN
!         
!            ! Surface area [cm2]
!            AREA_CM2 = GET_AREA_CM2( J )
!            
!            ! Convert from [v/v/timestep] to [molec/cm2/s]
!            FLUX = ( NIT0 - NIT ) * AD(I,J,1) / TCVV(IDTNIT)
!            FLUX = FLUX * XNUMOL(IDTNIT) / AREA_CM2 / DTCHEM
!         
!            ! Store in ND44
!            AD44(I,J,DRYNIT,1) = AD44(I,J,DRYNIT,1) + FLUX
!         ENDIF
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!------------------------------------------------------------------------------

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FREQ, NIT0, NIT, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only apply drydep to boxes w/in the PBL
         IF ( PBLFRAC(I,J,L) > 0d0 ) THEN 

            ! NIT drydep frequency [1/s] -- PBLFRAC accounts for the fraction
            ! of each vertical level that is located below the PBL top
            FREQ = DEPSAV(I,J,DRYNIT) * PBLFRAC(I,J,L)

            ! Only apply drydep loss if FREQ is nonzero
            IF ( FREQ > 0d0 ) THEN

               ! Initial NITRATE [v/v]
               NIT0 = STT(I,J,L,IDTNIT)

               ! Amount of NITRATE lost to drydep [v/v]
               NIT = NIT0 * ( 1d0 - EXP( -FREQ * DTCHEM ) )

               ! Prevent underflow condition
               IF ( NIT < SMALLNUM ) NIT = 0d0

               ! Subtract NITRATE lost to drydep from initial NITRATE [v/v]
               STT(I,J,L,IDTNIT) = NIT0 - NIT

               !========================================================
               ! ND44 Diagnostic: Drydep flux of NIT [molec/cm2/s]
               !========================================================
               IF ( ND44 > 0 ) THEN
         
                  ! Surface area [cm2]
                  AREA_CM2 = GET_AREA_CM2( J )
            
                  ! Convert from [v/v/timestep] to [molec/cm2/s]
                  FLUX = NIT  * AD(I,J,L)      / TCVV(IDTNIT)
                  FLUX = FLUX * XNUMOL(IDTNIT) / AREA_CM2 / DTCHEM
         
                  ! Store flux for ND44
                  AD44(I,J,DRYNIT,1) = AD44(I,J,DRYNIT,1) + FLUX
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_NIT

!-----------------------------------------------------------------------------

      SUBROUTINE EMISSSULFATE
!
!******************************************************************************
!  Subroutine EMISSSULFATE is the interface between the GEOS-CHEM model and
!  the sulfate emissions routines in "sulfate_mod.f" (bmy, 6/7/00, 11/6/02)
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : CHECK_STT
      USE TIME_MOD,    ONLY : GET_SEASON, GET_MONTH
      USE TRACERID_MOD

#     include "CMN_SIZE" ! Size parameters
#     include "CMN"      ! STT, NSRCX

      ! Local variables
      LOGICAL, SAVE      :: FIRSTEMISS = .TRUE. 
      INTEGER, SAVE      :: LASTMONTH  = -99
      INTEGER            :: NSEASON, MONTH   

      !=================================================================
      ! EMISSSULFATE begins here!
      !=================================================================

      ! Do only on the first timestep
      IF ( FIRSTEMISS ) THEN

         ! Initialize arrays
         CALL INIT_SULFATE

         ! Read emissions from volcanoes
         IF ( LENV ) CALL READ_NONERUP_VOLC
         IF ( LEEV ) CALL READ_ERUP_VOLC

         ! We have now gone thru the first timestep
         FIRSTEMISS  = .FALSE.
      ENDIF

      ! Get the season and month
      NSEASON = GET_SEASON()
      MONTH   = GET_MONTH()

      !=================================================================
      ! If this is a new month, read in the monthly mean quantities
      !=================================================================
      IF ( MONTH /= LASTMONTH ) THEN 
      
         ! Read monthly mean data
         CALL READ_SST( MONTH )
         CALL READ_OCEAN_DMS( MONTH )
         CALL READ_AIRCRAFT_SO2( MONTH )
         CALL READ_BIOMASS_SO2( MONTH )
         CALL READ_ANTHRO_SOx( MONTH, NSEASON )
         CALL READ_ANTHRO_NH3( MONTH )
         CALL READ_BIOMASS_NH3( MONTH )
         CALL READ_NATURAL_NH3( MONTH )
      
         ! Read oxidants for the offline simulation only
         IF ( NSRCX == 10 ) CALL READ_OXIDANT( MONTH )
      
         ! Save for next month
         LASTMONTH = MONTH
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

!-----------------------------------------------------------------------------

      SUBROUTINE SRCDMS( TC )
!
!***************************************************************************** 
!  Subroutine SRCDMS, from Mian Chin's GOCART model, add DMS emissions 
!  to the tracer array.  Modified for use with the GEOS-CHEM model.
!  (bmy, 6/2/00, 3/27/03)
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
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG_MOD,     ONLY : AD13_DMS
      USE DAO_MOD,      ONLY : IS_WATER, LWI, PBL
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_TS_EMIS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, XTRA2
#     include "CMN_DIAG"     ! ND13 (for now)

      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I,      J,    L,     NTOP
      REAL*8                 :: DTSRCE, SST,  Sc,    CONC,   W10 
      REAL*8                 :: ScCO2,  AKw,  ERATE, DMSSRC, BLTOP
      REAL*8                 :: P1,     P2,   DELP,  FEMIS

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
!$OMP+PRIVATE( I, J, SST, Sc, CONC, W10, ScCO2, AKw, ERATE )
!$OMP+PRIVATE( DMSSRC, NTOP, BLTOP, L, P1, P2, DELP, FEMIS )
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
            ! Liss and Merlivat (1986) 
            ScCO2 = 600.d0
            IF ( W10 <= 3.6d0 ) then
               AKw = 0.17d0 * W10
               
            ELSE IF ( W10 <= 13.d0 ) THEN
               AKw = 2.85d0 * W10 - 9.65d0
               
            ELSE
               AKw = 5.90d0 * W10 - 49.3d0
               
            ENDIF
            !-----------------------------------------------------------

            IF ( W10 <= 3.6d0 ) THEN
               AKw = AKw * ( (ScCO2/Sc)**0.667 )
            ELSE
               AKw = AKw * SQRT(ScCO2/Sc)
            ENDIF

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
            NTOP = CEILING( XTRA2(I,J) )
            
            ! PBL height is in the 3rd model layer or higher
            IF ( NTOP >= 2 ) THEN

               ! PBL top pressure [hPa]
               BLTOP = GET_PEDGE(I,J,1) - PBL(I,J)

               ! Loop thru the boundary layer
               DO L = 1, NTOP

                  ! DELP is the pressure thickness of level L [hPa]
                  P1   = GET_PEDGE(I,J,L) 
                  P2   = GET_PEDGE(I,J,L+1)
                  DELP = P1 - P2

                  ! PBL top occurs at level L
                  IF ( BLTOP <= P2 )  THEN
                     FEMIS = DELP / PBL(I,J)
   
                  ! Level L lies completely w/in the PBL
                  ELSE IF ( BLTOP > P2 .AND. BLTOP < P1 ) THEN
                     FEMIS = ( P1 - BLTOP ) / PBL(I,J)               

                  ! Level L lies completely out of the PBL
                  ELSE IF ( BLTOP > P1 ) THEN
                     CYCLE       

                  ENDIF

                  ! Fraction of total DMS in level L
                  TC(I,J,L) = TC(I,J,L) + ( FEMIS * DMSSRC )
               ENDDO

            ELSE

               ! If PBL height and lower or similar to the second model layer
               !then surface emission is emitted to the first model layer.
               TC(I,J,1) = TC(I,J,1) + DMSSRC 

            ENDIF 

         ELSE                   

            ! If we are not over water, then there is no DMS source
            DMSSRC = 0.d0

         ENDIF                  

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
!  aircraft, biomass, and anthro sources. (rjp, bdf, bmy, 6/2/00, 3/27/03)
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
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG_MOD,     ONLY : AD13_SO2_an, AD13_SO2_ac, AD13_SO2_bb,
     &                         AD13_SO2_nv, AD13_SO2_ev, AD13_SO2_bf
      USE DAO_MOD,      ONLY : BXHEIGHT, PBL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_TS_EMIS, GET_DAY_OF_YEAR

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! XTRA2
#     include "CMN_DIAG"     ! ND13, LD13 (for now)
  
      ! Arguments
      INTEGER, INTENT(IN)    :: NSEASON
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, K, L, LV1, LV2, NTOP, JDAY
      REAL*8                 :: ZH(0:LLPAR), DZ(LLPAR), SO2(LLPAR)
      REAL*8                 :: DTSRCE,      HGHT,      SO2SRC
      REAL*8                 :: SLAB,        SLAB1,     BLTOP
      REAL*8                 :: TSO2,        P1,        P2
      REAL*8                 :: DELP,        FEMIS

      ! Ratio of molecular weights: S/SO2
      REAL*8,  PARAMETER     :: S_SO2 = 32d0 / 64d0

      !=================================================================
      ! SRCSO2 begins here!
      !================================================================

      ! DTSRCE is the emission timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      ! JDAY is the day of year (0-365 or 0-366)
      JDAY   = GET_DAY_OF_YEAR()

      !=================================================================
      ! SO2 emissions from non-eruptive volcanoes [kg SO2/box/s].
      ! Assume that emission only occurs at the crater altitude.
      !=================================================================
      IF ( LENV ) THEN

         ! Initialize
         ESO2_nv = 0.d0

         ! Loop thru each non-erupting volcano
         DO K = 1, NNVOL

            ! Elevation of volcano crater
            HGHT  = DBLE( IELVn(k) )

            ! Altitude of crater from the ground
            ZH(0) = 0.d0

            ! Loop over levels
            DO L = 1, LLPAR

               ! Thickness of layer [m] w/ crater
               DZ(L) = BXHEIGHT(INV(K),JNV(K),L)

               ! Increment altitude
               ZH(L) = ZH(L-1) + DZ(L)

               ! If we are at the crater altitude, add emissions and exit
               IF ( ZH(L-1) <= HGHT .AND. ZH(L) > HGHT ) THEN 
                  ESO2_nv(INV(K),JNV(K),L) = 
     &                 ESO2_nv(INV(K),JNV(K),L) + Env(K)
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDIF  

      !=================================================================
      ! Calculate eruptive volcanic emission of SO2.
      !=================================================================
      IF ( LEEV ) THEN

         ! Initialize
         ESO2_ev = 0.D0

         ! Loop thru each erupting volcano
         DO K = 1, NEVOL

            ! Test to see if the volcano is erupting
            IF ( JDAY < IDAYS(K) .OR. JDAY > IDAYe(K) ) GOTO 20

            !===========================================================
            ! Define a slab at the top 1/3 of the volcano plume.
            !===========================================================
            HGHT  = DBLE( IHGHT(K) )

            ! slab bottom height
            SLAB1 = HGHT - ( HGHT - DBLE ( IELVe(K) ) ) / 3.d0 

            ! Slab thickness
            SLAB  = HGHT - SLAB1 
            ZH(0) = 0.d0 
        
            ! Loop over each level
            DO L = 1, LLPAR

               ! DZ is the thickness of level L [m]
               DZ(L) = BXHEIGHT(IEV(K),JEV(K),L)

               ! ZH is the height of the top edge of 
               ! level L, measured from the ground up [m]
               ZH(L) = ZH(L-1) + DZ(L) 

               ! max model erup.height
               IF ( L == LLPAR .AND. HGHT > ZH(L) ) THEN 
                  LV2 = LLPAR
                  !HGHT = ZH(L)
                  !SLAB1 = SLAB1 - ( HGHT - ZH(L) )
               ENDIF

               !========================================================
               ! If Zh(l) <= bottom of the slab, go to next level.
               !========================================================
               IF ( ZH(L) <= SLAB1 ) GOTO 22

               !========================================================
               ! If the slab is only in current level: CASE 1
               !       ---------------------------------- ZH(L)
               ! HGHT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! SLAB1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
               !       ---------------------------------- ZH(L-1)
               !========================================================
               IF ( ZH(L-1) <= SLAB1 .AND. ZH(L) >= HGHT ) THEN
                  LV1   = L
                  LV2   = L
                  DZ(L) = SLAB

               !========================================================
               ! The slab extends more then one level.  Find the 
               ! lowest (lv1) and the highest (lv2) levels:
               !       --------------------------------- ZH(L)
               ! HGHT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !       --------------------------------- ZH(L-1)
               ! 
               !       --------------------------------- ZH(L)
               ! SLAB1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               !       --------------------------------- ZH(L-1)
               !========================================================
               ELSE IF (ZH(L-1) <= SLAB1 .AND. ZH(L) > SLAB1)  THEN
                  LV1   = L
                  DZ(L) = ZH(L) - SLAB1
                  
               ELSE IF (ZH(L-1) < HGHT  .AND. ZH(L) > HGHT )  THEN
                  LV2   = L
                  DZ(L) = HGHT - ZH(L-1)
                  EXIT          ! do 20 
        
               ENDIF
               
               ! Go to next level
 22            CONTINUE         
            ENDDO

            !===========================================================
            ! Calculate SO2 emission in the levels between LV1 and LV2.  
            ! Convert Eev from [kg SO2/box/event] to [kg SO2/box/s].  
            ! ESO2_ev is distributed evenly with altitude among the slab.
            !===========================================================
            DO L = LV1, LV2
               ESO2_ev(IEV(K),JEV(K),L) = ESO2_ev(IEV(K),JEV(K),L) + 
     &              EEV(K) / ( (IDAYe(K)-IDAYs(K)+1) * 24.d0 * 3600.d0 )
     &              * DZ(L) / SLAB
            ENDDO

            ! Go to next volcano
 20         CONTINUE 
         ENDDO

      ENDIF  

      !=================================================================
      ! Add SO2 emissions into model levels
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,  NTOP, L,    SO2,   TSO2   ) 
!$OMP+PRIVATE( BLTOP, P1, P2,   DELP, FEMIS, SO2SRC )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top of the boundary layer
         NTOP  = CEILING( XTRA2(I,J) )

         ! Zero SO2 array
         DO L = 1, LLPAR
            SO2(L) = 0d0
         ENDDO

         ! Sum of anthro (surface + 100m), biomass, biofuel SO2 at (I,J)
         TSO2  = SUM( ESO2_an(I,J,:) ) + ESO2_bb(I,J) + ESO2_bf(I,J)

         ! Zero SO2SRC
         SO2SRC = 0d0

         !===============================================================
         ! Partition the total anthro and biomass SO2 emissions thru
         ! the entire boundary layer (if PBL top is higher than level 2)
         !===============================================================
         IF ( NTOP > 2 ) THEN

            ! BLTOP = pressure of PBL top [hPa]
            BLTOP = GET_PEDGE(I,J,1) - PBL(I,J)

            ! Loop thru levels in the PBL
            DO L  = 1, NTOP

               ! DELP is the pressure thickness of level K
               P1   = GET_PEDGE(I,J,L)
               P2   = GET_PEDGE(I,J,L+1)
               DELP = P1 - P2

               ! PBL top occurs w/in level L
               IF ( BLTOP <= P2 )  THEN
                  FEMIS = DELP / PBL(I,J)

               ! Level L lies completely w/in the PBL
               ELSE IF ( BLTOP >  P2 .AND. BLTOP < P1 ) THEN
                  FEMIS = ( P1 - BLTOP ) / PBL(I,J)

               ! Level L lies completely out of the PBL
               ELSE IF ( BLTOP > P1 ) THEN
                  CYCLE

               ENDIF

               ! Partition total SO2 into level K
               SO2(L) = FEMIS * TSO2
            ENDDO

         !===============================================================
         ! If PBL top occurs lower than or close to the top of level 2,
         ! then then surface SO2 goes into level 1 and the smokestack
         ! stack SO2 goes into level 2.
         !===============================================================
         ELSE
            SO2(1) = ESO2_an(I,J,1) + ESO2_bb(I,J) + ESO2_bf(I,J)
            SO2(2) = ESO2_an(I,J,2) 

         ENDIF 

         ! Error check
         IF ( ABS( SUM( SO2 ) - TSO2 ) > 1.D-5 ) THEN
!$OMP CRITICAL
            PRINT*, '### ERROR in SRCSO2!'
            PRINT*, '### I, J, L, : ', I, J, L
            PRINT*, '### SUM(SO2) : ', SUM( SO2 )
            PRINT*, '### TSO2     : ', TSO2
!$OMP END CRITICAL
            CALL ERROR_STOP( 'Check SO2 redistribution!',
     &                       'SRCSO2 (sulfate_mod.f)' )
         ENDIF
              
         !==============================================================
         ! Add anthro SO2, aircraft SO2, volcano SO2, and biomass SO2
         ! Convert from [kg SO2/box/s] -> [kg SO2/box/timestep]
         !==============================================================
         DO L = 1, LLPAR
            
            ! SO2 emissions [kg/box/s]
            SO2SRC = SO2(L)         + ESO2_ac(I,J,L) + 
     &               ESO2_nv(I,J,L) + ESO2_ev(I,J,L) 

            ! Add SO2 to TC array [kg/box/timestep]
            TC(I,J,L) = TC(I,J,L) + ( SO2SRC * DTSRCE )

         ENDDO

         !==============================================================
         ! ND13 Diagnostic: SO2 emissions in [kg S/box/timestep]
         !==============================================================
         IF ( ND13 > 0 ) THEN 

            ! Anthropogenic SO2 -- Levels 1-2
            DO L = 1, 2
               AD13_SO2_an(I,J,L) = AD13_SO2_an(I,J,L) +
     &                              ( ESO2_an(I,J,L) * S_SO2 * DTSRCE )
            ENDDO

            ! SO2 from biomass burning
            AD13_SO2_bb(I,J)      = AD13_SO2_bb(I,J) +
     &                              ( ESO2_bb(I,J) * S_SO2 * DTSRCE )
 
            ! SO2 from biofuel burning
            AD13_SO2_bf(I,J)      = AD13_SO2_bf(I,J) +
     &                              ( ESO2_bf(I,J) * S_SO2 * DTSRCE )

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
!***************************************************************************** 
!  Subroutine SRCSO4 (originally from Mian Chin) computes SO4 emissions from 
!  anthropogenic sources (rjp, bdf, bmy, 6/2/00, 3/27/03)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  (2) TC     (REAL*8 ) : Array for SO4 mass [kg]
! 
!  NOTES:
!  (1 ) Emission of SO4 is read in SULFATE_READYR, in [kg/box/s]. 
!        It is converted to [kg/box/timestep] here.  
!  (2 ) Now use routine GET_TS_EMIS from the new "time_mod.f" (bmy, 3/27/03)
!*****************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,      ONLY : PBL
      USE DIAG_MOD,     ONLY : AD13_SO4_an
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_TS_EMIS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! XTRA2
#     include "CMN_DIAG"     ! ND13 (for now)

      ! Arguments      
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, K, L, NTOP
      REAL*8                 :: SO4(LLPAR), DTSRCE, BLTOP, TSO4
      REAL*8                 :: P1,         P2,     DELP,  FEMIS

      ! Ratio of molecular weights: S/SO4
      REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0

      !=================================================================
      ! SRCSO4 begins here!
      !=================================================================

      ! DTSRCE is the emission timestep in seconds
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Compute SO4 emissions 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, SO4, TSO4, L, BLTOP, P1, P2, DELP, FEMIS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Top level of boundary layer at (I,J)
         NTOP = CEILING( XTRA2(I,J) )

         ! Zero SO4 array at all levels 
         DO L = 1, LLPAR
            SO4(L) = 0.0
         ENDDO

         ! Compute total anthropogenic SO4 (surface + 100m)
         TSO4 = SUM( ESO4_an(I,J,:) )

         !==============================================================
         ! Partition the total anthro SO4 emissions thru the entire 
         ! boundary layer (if PBL top is higher than level 2)
         !==============================================================
         IF ( NTOP > 2 ) THEN

            ! Boundary layer top
            BLTOP = GET_PEDGE(I,J,1) - PBL(I,J)
           
            ! Loop thru boundary layer
            DO L = 1, NTOP

               ! DELP is the pressure thickness of level L [hPa]
               P1   = GET_PEDGE(I,J,L)
               P2   = GET_PEDGE(I,J,L+1)
               DELP = P1 - P2

               ! PBL top occurs w/in level L
               IF ( BLTOP <= P2 )  THEN
                  FEMIS = DELP / PBL(I,J)

               ! Level L lies completely w/in the PBL
               ELSE IF ( BLTOP >  P2 .AND. BLTOP < P1 ) THEN
                  FEMIS = ( P1 - BLTOP ) / PBL(I,J)               

               ! Level L lies completely out of the PBL
               ELSE IF ( BLTOP > P1 ) THEN
                  CYCLE 

               ENDIF

               ! Fraction of total SO4 in layer L
               SO4(L) = FEMIS * TSO4
            ENDDO

         !==============================================================
         ! If PBL height is low and lower or similar to the second 
         ! model layer then surface emission is emitted to the first
         ! model layer and the stack emission goes to the second model
         ! layer
         !==============================================================
         ELSE

            SO4(1) = ESO4_an(I,J,1)
            SO4(2) = ESO4_an(I,J,2) 
            
         ENDIF 

         IF ( ABS( SUM( SO4 ) - TSO4 ) > 1.D-5 ) THEN
!$OMP CRITICAL
            PRINT*, '### ERROR in SRCSO4!'
            PRINT*, '### I, J, L, : ', I, J, L
            PRINT*, '### SUM(SO4) : ', SUM( SO4 )
            PRINT*, '### TSO4     : ', TSO4
!$OMP END CRITICAL
            CALL ERROR_STOP( 'Check SO4 redistribution',
     &                       'SRCSO4 (sulfate_mod.f)' )
         ENDIF

         !=============================================================
         ! Add SO4 emissions to tracer array 
         ! Convert from [kg SO4/box/s] -> [kg SO4/box/timestep]
         !=============================================================
         DO L = 1, LLPAR
            TC(I,J,L) = TC(I,J,L) + ( SO4(L) * DTSRCE )
         ENDDO

         !==============================================================
         ! ND13 Diagnostic: SO4 emission in [kg S/box/timestep]       
         !==============================================================
         IF ( ND13 > 0 ) THEN 
            DO L = 1, 2      
               AD13_SO4_an(I,J,L) = AD13_SO4_an(I,J,L) + 
     &                              ( ESO4_an(I,J,L) * S_SO4 * DTSRCE )
            ENDDO
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
!  (rjp, bmy, 12/17/01, 3/27/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD13_NH3_an, AD13_NH3_bb,
     &                         AD13_NH3_bf, AD13_NH3_na
      USE DAO_MOD,      ONLY : PBL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_TS_EMIS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT, XTRA2
#     include "CMN_DIAG"     ! ND13
      
      ! Argumetns
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)

      ! Local variables
      INTEGER                :: I, J, L, K, NTOP
      REAL*8                 :: BLTOP, P1, P2, DELP 
      REAL*8                 :: FEMIS, DTSRCE, NH3SRC, TNH3

      !=================================================================
      ! SRCNH3 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0
 
      ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NTOP, NH3SRC, TNH3, BLTOP, L, P1, P2, DELP, FEMIS )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
          
         ! Layer where the PBL top happens
         NTOP   = CEILING( XTRA2(I,J) )

         ! Initialize
         NH3SRC = 0.d0

         ! Sum all types of NH3 emission [kg/box/s]
         TNH3   = ENH3_an(I,J) + ENH3_bb(I,J) + 
     &            ENH3_bf(I,J) + ENH3_na(I,J)

         !==============================================================
         ! Add NH3 emissions [kg NH3/box] into the tracer array
         ! Partition total NH3 throughout the entire boundary layer
         !==============================================================
         IF ( NTOP >= 2 ) THEN

            ! Boundary layer top pressure [hPa]
            BLTOP = GET_PEDGE(I,J,1) - PBL(I,J)

            ! Loop over all levels in the boundary layer
            DO L = 1, NTOP

               ! DELP is the pressure thickness of level K [hPa]
               P1   = GET_PEDGE(I,J,L)
               P2   = GET_PEDGE(I,J,L+1)
               DELP = P1 - P2

               ! Case of model grid is lower than PBL
               IF ( BLTOP <= P2 )  THEN
                  FEMIS = DELP / PBL(I,J)

               ELSE IF ( BLTOP >  P2 .AND. BLTOP < P1 ) THEN
                  FEMIS = ( P1 - BLTOP ) / PBL(I,J)               

               ELSE IF ( BLTOP > P1 ) THEN
                  CYCLE       
               ENDIF

               ! Partition total NH3 into level K [kg NH3/s]
               ! This is just for error checking
               NH3SRC    = NH3SRC  + ( FEMIS * TNH3 )

               ! Add NH3 emissions into tracer array [kg NH3/timestep]
               TC(I,J,L) = TC(I,J,L) + ( TNH3 * FEMIS * DTSRCE )
            ENDDO

            ! Error check
            IF ( ABS( NH3SRC - TNH3 ) > 1.D-5 ) THEN
!$OMP CRITICAL
               PRINT*, '### ERROR in SRCNH3!'
               PRINT*, '### I, J         : ', I, J
               PRINT*, '### NH3SRC       : ', NH3SRC
               PRINT*, '### ENH3_an(I,J) : ', ENH3_an(I,J) 
               PRINT*, '### ENH3_bb(I,J) : ', ENH3_bb(I,J)
               PRINT*, '### ENH3_bf(I,J) : ', ENH3_bf(I,J)
!$OMP END CRITICAL
               CALL ERROR_STOP( 'Check NH3 redistribution', 
     &                          'SRCNH3 (sulfate_mod.f)' )
            ENDIF

         !============================================================
         ! If PBL height close to the top of level 2, then put all of
         ! the emission into the surface layer [kg NH3/box/timestep]
         !============================================================
         ELSE
            TC(I,J,1) = TC(I,J,1) + ( TNH3 * DTSRCE )
         ENDIF

         !============================================================
         ! ND13 diagnostics: NH3 emissions [kg NH3/box/timestep]
         !============================================================
         IF ( ND13 > 0 ) THEN

            ! Anthro NH3
            AD13_NH3_an(I,J) = AD13_NH3_an(I,J) + 
     &                         ( ENH3_an(I,J) * DTSRCE )
            
            ! Biomass NH3
            AD13_NH3_bb(I,J) = AD13_NH3_bb(I,J) + 
     &                         ( ENH3_bb(I,J) * DTSRCE )
                  
            ! Biofuel NH3
            AD13_NH3_bf(I,J) = AD13_NH3_bf(I,J) +
     &                         ( ENH3_bf(I,J) * DTSRCE )

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
!  OH for offline simulations. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDOH (bmy, 11/1/02)
!  (2 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DAO_MOD,       ONLY : SUNCOS
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACERID_MOD,  ONLY : IDOH

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================
      SELECT CASE ( NSRCX ) 
      
         !---------------------
         ! Coupled simulation
         !---------------------
         CASE ( 3 )

            ! JLOOP = SMVGEAR 1-D grid box index
            JLOOP = JLOP(I,J,L)

            ! Take OH from the SMVGEAR array CSPEC
            ! OH is defined only in the troposphere
            IF ( JLOOP > 0 ) THEN
               OH_MOLEC_CM3 = CSPEC(JLOOP,IDOH)
            ELSE
               OH_MOLEC_CM3 = 0d0
            ENDIF

         !---------------------
         ! Offline simulation
         !---------------------
         CASE ( 10 )

            ! 1-D grid box index for SUNCOS
            JLOOP = IJSURF(I,J)

            ! Test for sunlight...
            IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

               ! Impose a diurnal variation on OH during the day
               OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                        ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                        ( 1440d0        / GET_TS_CHEM() )

               ! Make sure OH is not negative
               OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
            ELSE

               ! At night, OH goes to zero
               OH_MOLEC_CM3 = 0d0

            ENDIF

         ! Error check
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      END SELECT

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
 
#     include "CMN_SIZE"  ! Size parameters

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

      ! Replace OH into CSPEC (troposphere only)
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
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : CSPEC, JLOP
      USE DAO_MOD,        ONLY : AD,    SUNCOS
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE GLOBAL_NO3_MOD, ONLY : NO3
      USE TRACERID_MOD,   ONLY : IDNO3

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX

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
      SELECT CASE ( NSRCX )

         ! Coupled chemistry/aerosol simulation
         CASE ( 3 )
            
            ! 1-D SMVGEAR grid box index
            JLOOP = JLOP(I,J,L)

            ! Take NO3 from the SMVGEAR array CSPEC
            ! NO3 is defined only in the troposphere
            IF ( JLOOP > 0 ) THEN
               NO3_MOLEC_CM3 = CSPEC(JLOOP,IDNO3)
            ELSE
               NO3_MOLEC_CM3 = 0d0
            ENDIF

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
         CASE ( 10 )

            ! 1-D grid box index for SUNCOS
            JLOOP = IJSURF(I,J)

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

         ! Error check
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid NSRCX!','GET_NO3 (sulfate_mod.f)')

      END SELECT

      ! Return to calling program
      END FUNCTION GET_NO3

!------------------------------------------------------------------------------

      SUBROUTINE SET_NO3( I, J, L, NO3 ) 
!
!******************************************************************************
!  Function SET_NO3 saves the modified NO3 value back to SMVGEAR's CSPEC array
!  for coupled lfate/aerosol simulations. (rjp, bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!  (4  ) NO3     (REAL*8 ) : OH at grid box (I,J,L) to be saved into CSPEC
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDNO3. (bmy, 12/16/02)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACERID_MOD, ONLY : IDNO3
      
#     include "CMN_SIZE"  ! Size parameters 
#     include "CMN"       ! NSRCX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: NO3

      ! Local variables
      INTEGER             :: JLOOP

      !=================================================================
      ! SET_NO3 begins here!
      !=================================================================

      SELECT CASE ( NSRCX )

         !--------------------
         ! Coupled simulation
         !--------------------
         CASE ( 3 ) 

            ! 1-D grid box index for CSPEC
            JLOOP = JLOP(I,J,L) 

            ! Replace OH into CSPEC (troposphere only)
            IF ( JLOOP > 0 ) THEN
               CSPEC(JLOOP,IDNO3) = NO3
            ENDIF

         ! Error check
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid NSRCX!','SET_NO3 (sulfate_mod.f)')

      END SELECT

      ! Return to calling program
      END SUBROUTINE SET_NO3
      
!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_VV )
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
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC, JLOP, VOLUME
      USE DAO_MOD,      ONLY : AIRDEN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACERID_MOD, ONLY : IDO3

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX, LPAUSE
#     include "CMN_O3"    ! XNUMOLAIR

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: O3_VV
 
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================
      SELECT CASE ( NSRCX ) 

         !--------------------
         ! Coupled simulation
         !--------------------
         CASE ( 3 )

            ! JLOOP = SMVGEAR 1-D grid box index
            JLOOP = JLOP(I,J,L)

            ! Get O3 from CSPEC [molec/cm3] and convert it to [v/v]
            ! O3 data will only be defined below the tropopause
            IF ( JLOOP  > 0 ) THEN
               O3_VV = ( CSPEC(JLOOP,IDO3) * 1d6       ) / 
     &                 ( AIRDEN(L,I,J)     * XNUMOLAIR )
            ELSE
               O3_VV = 0d0
            ENDIF
         
         !--------------------
         ! Offline simulation
         !--------------------
         CASE ( 10 )

            ! Get O3 [v/v] for this gridbox & month
            ! O3 data will only be defined below the tropopause
            IF ( L <= LLTROP ) THEN
               O3_VV = O3m(I,J,L)
            ELSE
               O3_VV = 0d0
            ENDIF

         ! Error check
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid NSRCX!', 'GET_OH (sulfate_mod.f)')

      END SELECT

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_NONERUP_VOLC
!
!******************************************************************************
!  Subroutine READ_NONERUP_VOLC reads SO2 emissions from non-eruptive
!  volcanoes. (rjp, bdf, bmy, 9/19/02)
!
!  NOTES:
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!******************************************************************************
! 
      ! References to F90 modules
      USE BPCH2_MOD
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

      ! Local variables
      INTEGER              :: I, IOS, J, K, L
      REAL*8               :: EE
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_NONERUP_VOLC begins here!
      !=================================================================

      ! Initialize
      K        = 1
      Env      = 0.d0
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'sulfate_sim_200210/volcano.con.' // GET_RES_EXT()

      !=================================================================
      ! Read NON-eruptive volcanic SO2 emission (GEIA) into Env.  
      ! Convert Env from [Mg SO2/box/day] to [kg SO2/box/s].
      !=================================================================

      ! Fancy output
      WRITE( 6, 100 ) TRIM( FILENAME ) 
 100  FORMAT( '     - READ_NONERUP_VOLC: Reading ', a )
     
      ! Open file 
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:1' )

      ! Read header lines
      DO L = 1, 2
         READ( IU_FILE, '(a)', IOSTAT=IOS )             
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:2' )
         ENDIF
      ENDDO

      ! Read data values
      DO 
         READ( IU_FILE, '(49x,i4,e11.3,1x,2i4)', IOSTAT=IOS ) 
     &        IELVn(k), EE, INV(K), JNV(k) 
         
         ! Check for EOF
         IF ( IOS < 0 ) EXIT

         ! Trap I/O error
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'read_nonerup_volc:3' )
         ENDIF

         ! Unit conversion: [Mg SO2/box/day] -> [kg SO2/box/s]
         Env(k) = EE * 1000.d0 / ( 24.d0 * 3600.d0 )

         ! Increment counter
         K = K + 1 
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! NNVOL = Number of non-eruptive volcanoes
      NNVOL = K - 1

      ! Return to calling program
      END SUBROUTINE READ_NONERUP_VOLC

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_ERUP_VOLC
!
!***************************************************************************** 
!  Subroutine READ_ERUP_VOLC reads SO2 emissions from eruptive
!  volcanoes. (rjp, bdf, bmy, 9/19/02)
!
!  NOTES:
!  (1 ) Split off from old module routine "sulfate_readyr" (bmy, 9/19/02)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE FILE_MOD, ONLY : IU_FILE, IOERROR
      
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_SETUP"  ! DATA_DIR
 
      ! Local variables
      INTEGER              :: I, IOS, IUNIT, J, K, L, M
      REAL*8               :: A, B, Fe, X, EE
      CHARACTER(LEN=255)   :: FILENAME

      !==================================================================
      ! READ_ERUP_VOLC begins here
      !==================================================================

      ! Initialize
      K        = 1
      Eev(:)   = 0.d0
      FILENAME = TRIM( DATA_DIR )                        // 
     &           'sulfate_sim_200210/volcano.erup.1990.' // 
     &           GET_RES_EXT()
      
      !==================================================================
      ! Read eruptive volcanic SO2 emission (based on Smithsonian data 
      ! base, SO2 emission and cloud height are a function of VEI.  
      ! Data are over-written if TOMS observations are available.  
      ! Also define a slab with a thickness of 1/3 of the cloud column, 
      ! and SO2 are emitted uniformely within the slab.  
      !
      ! Convert Ee from [kton SO2] to [kg SO2/box] and store in Eev.
      ! ESO2_ev(i,j,l) in [kg SO2/box/s] will be calculated in SRCSO2.  
      !==================================================================
      
      ! Fancy output
      WRITE( 6, 100 ) TRIM( FILENAME ) 
 100  FORMAT( '     - READ_ERUP_VOLC: Reading ', a )
   
      ! Open file 
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_erup_volc:1' )

      ! Read header lines
      DO L = 1, 2
         READ( IU_FILE, '(a)', IOSTAT=IOS )
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'read_erup_volc:2' )
         ENDIF
      ENDDO

         ! Read data values
      DO 
         READ( IU_FILE, '(47x,3i6,6x,i6,es11.3,1x,2i4)', IOSTAT=IOS )
     &        IELVe(K), IDAYs(K), IDAYe(K), IHGHT(K), 
     &        Ee,       IEV(K),   JEV(K)
         
         ! Check for EOF
         IF ( IOS < 0 ) EXIT

          ! Trap I/O error
         IF ( IOS > 0 ) THEN
            CALL IOERROR( IOS, IU_FILE, 'sulfate_readyr:6' )
         ENDIF

         ! Unit conversion: [kton SO2/box/event] -> [kg SO2/box/event]
         Eev(k) = Ee * 1.d6

         ! Increment count
         K = K + 1
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      ! NEVOL = Number of eruptive volcanoes
      NEVOL = K - 1

      ! Return to calling program
      END SUBROUTINE READ_ERUP_VOLC
      
!------------------------------------------------------------------------------

      SUBROUTINE READ_ANTHRO_SOx( THISMONTH, NSEASON )
!
!******************************************************************************
!  Suborutine READ_ANTHRO_SOx reads the anthropogenic SOx from disk, 
!  and partitions it into anthropogenic SO2 and SO4. 
!  (rjp, bdf, bmy, 9/20/02, 3/27/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE GRID_MOD,     ONLY : GET_XMID, GET_YMID, GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_YEAR
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
      USE TRACERID_MOD, ONLY : IDTSO2, IDTSO4

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSEASON, FSCALYR
#     include "CMN_O3"    ! XNUMOL
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH, NSEASON

      ! Local variables
      INTEGER             :: I, J, L, IX, JX, IOS
      INTEGER, SAVE       :: LASTYEAR = -99
      REAL*4              :: E_SOx(IGLOB,JGLOB,2)
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, Fe, X, Y, AREA_CM2
      CHARACTER(LEN=4)    :: SYEAR
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! Read SOx emissions [molec SOx/cm2/s] 
      !=================================================================

      ! Define filename
      FILENAME = TRIM( DATA_DIR )                       //
     &           'fossil_200104/merge_nobiofuels.geos.' // 
     &           GET_RES_EXT() 
     
      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANTHRO_SOx: Reading ', a )

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
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 2,         E_SOx,     QUIET=.TRUE. )

      !=================================================================
      ! Read in yearly SO2 scale factors here
      ! (For now we only have 1998, deal w/ other years later)
      !=================================================================
      IF ( LASTYEAR < 0 ) THEN

         ! put in SOX scale year here (hardwired to 1998 for now)
         SYEAR    = '1998'
         FILENAME = TRIM( DATA_DIR )                    // 
     &              'sulfate_sim_200210/scalefoss.SOx.' //  
     &              GET_RES_EXT()  // '.' // SYEAR
        
         ! Echo output
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Get TAU value (use Jan 1, 1998 for scale factors)
         XTAU = GET_TAU0( 1, 1, 1998 )

         ! Read anthropogenic SOx [molec/cm2/s] 
         CALL READ_BPCH2( FILENAME, 'SCALFOSS', 3, 
     &                    XTAU,      IGLOB,     JGLOB,     
     &                    1,         ARRAY,     QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8
         CALL TRANSFER_2D( ARRAY(:,:,1), SOx_SCALE )
         
         ! Reset LASTYEAR
         LASTYEAR = GET_YEAR()
      ENDIF

      !=================================================================
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
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_CM2, Y, X, Fe )
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

            ! Compute SO4 (tracer #3) from SOx
            ! Convert from [molec SOx/cm2/s] to [kg SO4/box/s]
            ESO4_an(I,J,L) = E_SOx(I,J,L) * Fe *
     &                       AREA_CM2 / XNUMOL(IDTSO4)
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE READ_ANTHRO_SOx

!------------------------------------------------------------------------------

      SUBROUTINE READ_OCEAN_DMS( THISMONTH )
!
!***************************************************************************** 
!  Subroutine READ_OCEAN_DMS reads seawater concentrations of DMS (nmol/L).
!  (rjp, bdf, bmy, 9/20/02, 3/27/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters 
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH
      
      ! Local variables
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      !==================================================================
      ! READ_OCEAN_DMS begins here!
      !==================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                        // 
     &           'sulfate_sim_200210/DMS_seawater.geos.' //
     &           GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )  
 100  FORMAT( '     - READ_OCEAN_DMS: Reading ', a ) 

      ! TAU value (use generic year 1985)
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read ocean DMS [nmol/L]
      CALL READ_BPCH2( FILENAME, 'BIOGSRCE',    25, 
     &                 XTAU,      IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8 
      CALL TRANSFER_2D( ARRAY(:,:,1), DMSo )
      
      ! Return to calling program
      END SUBROUTINE READ_OCEAN_DMS

!------------------------------------------------------------------------------

      SUBROUTINE READ_SST( THISMONTH )
!
!***************************************************************************** 
!  Subroutine READ_SST reads monthly mean sea surface temperatures.
!  (rjp, bdf, bmy, 9/18/02, 3/27/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!  (2 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters 
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH
      
      ! Local variables
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      !==================================================================
      ! READ_SST begins here!
      !==================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )               // 
     &           'sulfate_sim_200210/SST.geos.' // GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )  
 100  FORMAT( '     - READ_SST: Reading ', a ) 

      ! TAU value (use generic year 1985)
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read sea surface temperature [K]
      CALL READ_BPCH2( FILENAME, 'DAO-FLDS',    5, 
     &                 XTAU,      IGLOB,        JGLOB,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8 
      CALL TRANSFER_2D( ARRAY(:,:,1), SSTEMP )
      
      ! Return to calling program
      END SUBROUTINE READ_SST

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOMASS_SO2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_BIOMASS_SOx reads monthly mean biomass burning 
!  emissions for SO2.  SOx is read in, and converted to SO2. 
!  (rjp, bdf, bmy, 1/16/03, 3/27/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRANSFER_MOD, ONLY : TRANSFER_2D
      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! MONTH
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: BIOCO(IIPAR,JJPAR)
      REAL*8              :: XTAU, AREA_CM2
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_BIOMASS_SO2 begins here!
      !
      ! Compute biomass SO2 from biomass CO.  Use a molar ratio 
      ! of 0.0026 moles SO2/mole CO. (rjp, bmy, 1/16/03)
      !=================================================================
      
      ! File name for climatological biomass burning
      FILENAME = TRIM( DATA_DIR )                       // 
     &           'biomass_200110/bioburn.seasonal.geos.'// GET_RES_EXT()

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BIOMASS_SO2: Reading ', a )

      ! Get TAU0 value (use generic year 1985)
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read Biomass burning of CO [molec/cm2/month]
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE',    4, 
     &                 XTAU,      IGLOB,        JGLOB,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), BIOCO )

      ! Convert from [molec CO/cm2/month] to [kg SO2/box/s]
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )
         
         DO I = 1, IIPAR
            ESO2_bb(I,J) = BIOCO(I,J) * AREA_CM2 * 
     &                     0.0026d0 * 64d-3 * 12.d0 /
     &                    ( 6.022d23 * 86400.d0 * 365.25d0 )
         ENDDO
      ENDDO

      !=================================================================
      ! Compute biofuel SO2 from biofuel CO.  Use a molar 
      ! ratio of 0.0015 moles SO2/mole CO from biofuel burning. 
      ! (Table 2, [Andreae and Merlet, 2001])
      !=================================================================
      
      ! File name for biofuel burning file
      FILENAME = TRIM( DATA_DIR )               // 
     &           'biofuel_200202/biofuel.geos.' // GET_RES_EXT()

      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Get TAU0 value (use generic year 1985)
      XTAU = GET_TAU0( 1, 1, 1985 )

      ! Read Biofuel burning of CO [kg/yr]
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    4, 
     &                 XTAU,      IGLOB,        JGLOB,     
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE.  ) 

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), BIOCO )

      ! Convert from [kg CO/box/yr] to [kg SO2/box/s]
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ESO2_bf(I,J) = ( BIOCO(I,J) * 64d-3 * 0.0015d0 /
     &                    ( 28D-3 * 86400.d0 * 365.25d0 ) )
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_BIOMASS_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_AIRCRAFT_SO2( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_AIRCRAFT_SO2 reads monthly mean aircraft fuel emissions 
!  and converts them to SO2 emissions. (rjp, bdf, bmy, 9/18/02)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Extracted from old module routine SULFATE_READMON (bmy, 9/18/02)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,  ONLY : BXHEIGHT
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, IOS, J, K, L
      REAL*8              :: ACSO2(IGLOB,JGLOB,20)
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
     &           'sulfate_sim_200210/aircraft.' // GET_RES_EXT() //
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

      SUBROUTINE READ_ANTHRO_NH3( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_ANTHRO_NH3 reads the monthly mean anthropogenic 
!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 3/27/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_ANTHRO_NH3 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                        //
     &           'sulfate_sim_200210/NH3_anthsrce.geos.' //
     &           GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANTHRO_NH3: Reading ', a )
      
      ! Get TAU value (use year 1990, the year of the data!)
      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
	
      ! Read 1990 NH3 emissions [kg N/box/mon]
      CALL READ_BPCH2( FILENAME, 'NH3-ANTH',    29,  
     &                 XTAU,      IGLOB,        JGLOB,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_an )

      ! Convert from [kg N/box/mon] to [kg NH3/box/s]
      ENH3_an = ENH3_an * ( 17.d0 / 14.d0 ) 
     &        / ( NMDAY(THISMONTH) * 86400.d0 ) 
 
      ! Return to calling program
      END SUBROUTINE READ_ANTHRO_NH3

!------------------------------------------------------------------------------

      SUBROUTINE READ_NATURAL_NH3( THISMONTH )
!
!******************************************************************************
!  Subroutine READ_NATURAL_NH3 reads the monthly mean natural 
!  NH3 emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 4/8/03)
!
!  Arguments as input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Updated FORMAT string.  Now also call READ_BPCH2 with QUIET=.TRUE.
!        (bmy, 4/8/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_NATURAL_NH3 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                        //
     &           'sulfate_sim_200210/NH3_natusrce.geos.' //
     &           GET_RES_EXT()

      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_NATURAL_NH3: Reading ', a )
      
      ! Get TAU value (use year 1990, the year of the data!)
      XTAU = GET_TAU0( THISMONTH, 1, 1990 )
	
      ! Read 1990 NH3 emissions [kg N/box/mon]
      CALL READ_BPCH2( FILENAME, 'NH3-NATU',    29,  
     &                 XTAU,      IGLOB,        JGLOB,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_na )

      ! Convert from [kg N/box/mon] to [kg NH3/box/s]
      ENH3_na = ENH3_na * ( 17.d0 / 14.d0 ) /
     &          ( NMDAY(THISMONTH) * 86400.d0 ) 
 
      ! Return to calling program
      END SUBROUTINE READ_NATURAL_NH3

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOMASS_NH3( THISMONTH ) 
!
!******************************************************************************
!  Subroutine READ_BIOMASS_NH3 reads the monthly mean biomass NH3 
!  and biofuel emissions from disk and converts to [kg NH3/box/s]. 
!  (rjp, bdf, bmy, 9/20/02, 12/2/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TIME_MOD,     ONLY : GET_YEAR, GET_TAU
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! JYEAR
#     include "CMN_SETUP" ! DATA_DIR, LBBSEA

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J, YEAR
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU, DMASS, AREA_M2, TAU
      REAL*8              :: NMDAY(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                      31d0, 30d0, 31d0, 31d0, 
     &                                      30d0, 31d0, 30d0, 31d0 /)
      CHARACTER(LEN=4  )  :: CYEAR
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_BIOMASS_NH3 begins here!
      !=================================================================

      ! Current TAU value
      TAU = GET_TAU()

      ! Use seasonal or interannual emisisons?
      IF ( LBBSEA ) THEN

         !-----------------------------------------
         ! Use seasonal biomass emissions
         !-----------------------------------------

         ! File name for seasonal BB emissions
         FILENAME = TRIM( DATA_DIR )                        //
     &              'biomass_200110/bioburn.seasonal.geos.' //
     &              GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )      
         
      ELSE IF ( ( .not. LBBSEA ) .AND. 
     &          ( TAU < 101520d0 .or. TAU > 140256d0 ) ) THEN

         !-----------------------------------------
         ! Use seasonal biomass emissions as a
         ! proxy for missing interannual emissions
         !-----------------------------------------

         ! File name for seasonal BB emissions
         FILENAME = TRIM( DATA_DIR )                        //
     &              'biomass_200110/bioburn.seasonal.geos.' //
     &              GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( THISMONTH, 1, 1985 )   

      ELSE

         !-----------------------------------------
         ! Use interannual biomass emissions for
         ! years between 1996 and 2000 
         !-----------------------------------------

         ! Get year for interannual biomass emissions
         YEAR = MAX( MIN( GET_YEAR(), 2000 ), 1996 )

         ! Convert YEAR to a string
         ! Now use ENCODE to define CYEAR string for PGI/Linux (bmy, 9/29/03)
#if   defined ( LINUX_PGI ) 
         ENCODE( 4, '(i4)', CYEAR ) YEAR
#else
         WRITE( CYEAR, '(i4)' ) YEAR
#endif
    
         ! File name for interannual biomass burning emissions
         FILENAME = TRIM( DATA_DIR )                           //
     &              'biomass_200110/bioburn.interannual.geos.' //
     &               GET_RES_EXT() // '.' // CYEAR

         ! Get TAU0 value for the given year
         XTAU = GET_TAU0( THISMONTH, 1, YEAR ) 

      ENDIF

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BIOMASS_NH3: Reading ', a )

      ! Read TOTAL Drymass burned [g/cm2] as tracer #33
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE',    33, 
     &                 XTAU,      IGLOB,        JGLOB,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_bb )

      ! Loop over grid boxes
      DO J = 1, JJPAR
         
         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )
         
         DO I = 1, IIPAR

            ! Convert [g/cm2/month] to [kg/box/s] dry biomass burned
            DMASS = ENH3_bb(I,J) * AREA_M2 * 10.d0 /
     &              ( NMDAY(THISMONTH) * 86400d0 )

            ! Convert [kg/box/s] dry biomass to [kg NH3/box/s]
            ! Using emission factor of 0.0013 from Andreae & Merlet 2001
            ENH3_bb(I,J) = DMASS * 0.0013d0          
         ENDDO
      ENDDO
 
      !=================================================================
      ! Read NH3 biofuel emissions
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR )                       // 
     &           'sulfate_sim_200210/NH3_biofuel.geos.' // GET_RES_EXT()
   
      ! Echo output
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Get TAU0 value for 1998
      XTAU  = GET_TAU0( THISMONTH, 1, 1998 )

      ! Read NH3 biofuel data [kg NH3/box/month]
      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',    29, 
     &                 XTAU,      IGLOB,        JGLOB,       
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )
 
      ! Cast from REAL*4 to REAL*8 and resize if necesary
      CALL TRANSFER_2D( ARRAY(:,:,1), ENH3_bf )

      ! Store NH3 in ENH3_bf array [kg NH3/box/s]
      ENH3_bf = ENH3_bf / ( NMDAY(THISMONTH) * 86400.d0 )
 
      ! Return to calling program
      END SUBROUTINE READ_BIOMASS_NH3

!------------------------------------------------------------------------------

      SUBROUTINE READ_OXIDANT( MONTH )
!
!******************************************************************************
!  Subroutine READ_OXIDANT reads in monthly mean H2O2 and O3 fields for the
!  offline sulfate simulation. (rjp, bdf, bmy, 11/1/02, 3/27/03)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) MONTH    (INTEGER  ) : Emission timestep in minutes
!
!  NOTES:
!  (1 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/27/03)
!******************************************************************************
!  
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_3D_TROP

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: MONTH

      ! Local variables 
      INTEGER             :: I, J, L, K      
      REAL*4              :: ARRAY(IGLOB,JGLOB,LLTROP)
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
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/PH2O2.' //
     &           GET_NAME_EXT()    // '.' // GET_RES_EXT()

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_OXIDANT: Reading ', a ) 

      ! Read data
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 20,     
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_3D_TROP( ARRAY, PH2O2m )
            
      !=================================================================
      ! Read monthly mean O3 [v/v]
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/O3.' //
     &           GET_NAME_EXT()    // '.' // GET_RES_EXT()

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 32,     
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

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
!  (rjp, bmy, 12/16/02, 3/27/03)
!
!  NOTES:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC, 
     &                     GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"

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
!  declared in "sulfate_mod.f" (bmy, 6/2/00, 3/23/03)
!
!  NOTES:
!  (1 ) Only allocate some arrays for the standalone simulation (NSRCX==10).
!        Also reference NSRCX from "CMN".  Now eferences routine ALLOC_ERR
!        from "error_mod.f" ((rjp, bdf, bmy, 10/15/02)
!  (2 ) Now also allocate the IJSURF array to keep the 1-D grid box indices
!        for SUNCOS (for both coupled & offline runs).  Now allocate PH2O2m 
!        and O3m for offline runs.  Also allocate ESO2_bf (bmy, 1/16/03)
!  (3 ) Now allocate ENH3_na array (rjp, bmy, 3/23/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" ! Size parameters
#     include "CMN"      ! NSRCX

      ! Local variables
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS, I, J, IJLOOP

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

      ALLOCATE( EEV( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Eev' )
      Eev = 0d0

      ALLOCATE( ENV( NNV ), STAT=AS )
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

      ALLOCATE( ESO2_an( IIPAR, JJPAR, 2 ), STAT=AS )
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

      ALLOCATE( ESO4_an( IIPAR, JJPAR, 2  ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ESO4_an' )
      ESO4_an = 0d0
  
      ALLOCATE( IDAYe( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IDAYe' )
      IDAYe = 0d0

      ALLOCATE( IDAYs( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IDAYs' )
      IDAYs = 0d0

      ALLOCATE( IELVe( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IELVe' )
      IELVe = 0d0

      ALLOCATE( IELVn( NNV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IELVn' )
      IELVn = 0d0

      ALLOCATE( IEV( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IEV' )
      IEV = 0d0

      ALLOCATE( IHGHT( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IHGHT' )
      IHGHT = 0d0

      ALLOCATE( INV( NNV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'INV' )
      INV = 0d0

      ALLOCATE( JEV( NEV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JEV' )
      JEV = 0d0

      ALLOCATE( JNV( NNV ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JNV' )
      JNV = 0d0
         
      ALLOCATE( PMSA_DMS( IIPAR, JJPAR, LLTROP ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PMSA_DMS' )
      PMSA_DMS = 0d0

      ALLOCATE( PSO2_DMS( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO2_DMS' )
      PSO2_DMS = 0d0

      ALLOCATE( PSO4_SO2( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSO4_SO2' )
      PSO4_SO2 = 0d0

      ALLOCATE( SOx_SCALE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOx_SCALE' )
      SOx_SCALE = 0d0

      ALLOCATE( VCLDF( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VCLDF' )
      VCLDF = 0d0

      !=================================================================
      ! Only initialize the following for offline runs (NSRCX == 10)
      !=================================================================
      IF ( NSRCX == 10 ) THEN
         ALLOCATE( PH2O2m( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PH2O2m' )
         PH2O2m = 0d0

         ALLOCATE( O3m( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3m' )
         O3m = 0d0

         ALLOCATE( JH2O2( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3m' )
         JH2O2 = 0d0

         ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
         TCOSZ = 0d0

         ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
         TTDAY = 0d0
      ENDIF

      !=================================================================
      ! Initialize IJSURF, which indexes the 1-D SUNCOS array
      !=================================================================
      ALLOCATE( IJSURF( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IJSURF' )

      IJLOOP = 0
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IJLOOP      = IJLOOP +1 
         IJSURF(I,J) = IJLOOP
      ENDDO
      ENDDO

      ! Set IS_INIT
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_SULFATE

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SULFATE
!
!******************************************************************************
!  Subroutine CLEANUP_SULFATE deallocates all previously allocated arrays 
!  for sulfate emissions -- call at the end of the run (bmy, 6/1/00, 3/23/03)
! 
!  NOTES:
!  (1 ) Now also deallocates IJSURF. (bmy, 11/12/02)
!  (2 ) Now also deallocates ENH3_na (rjp, bmy, 3/23/03)
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
      IF ( ALLOCATED( ESO4_an   ) ) DEALLOCATE( ESO4_an   )
      IF ( ALLOCATED( IDAYs     ) ) DEALLOCATE( IDAYs     )
      IF ( ALLOCATED( IDAYe     ) ) DEALLOCATE( IDAYe     )
      IF ( ALLOCATED( IELVe     ) ) DEALLOCATE( IELVe     )
      IF ( ALLOCATED( IELVn     ) ) DEALLOCATE( IELVn     )
      IF ( ALLOCATED( IEV       ) ) DEALLOCATE( IEV       )
      IF ( ALLOCATED( IHGHT     ) ) DEALLOCATE( IHGHT     )
      IF ( ALLOCATED( IJSURF    ) ) DEALLOCATE( IJSURF    )
      IF ( ALLOCATED( INV       ) ) DEALLOCATE( INV       )
      IF ( ALLOCATED( JEV       ) ) DEALLOCATE( JEV       )
      IF ( ALLOCATED( JH2O2     ) ) DEALLOCATE( JH2O2     )
      IF ( ALLOCATED( JNV       ) ) DEALLOCATE( JNV       )
      IF ( ALLOCATED( O3m       ) ) DEALLOCATE( O3m       )
      IF ( ALLOCATED( PH2O2m    ) ) DEALLOCATE( PH2O2m    )
      IF ( ALLOCATED( PMSA_DMS  ) ) DEALLOCATE( PMSA_DMS  )
      IF ( ALLOCATED( PSO2_DMS  ) ) DEALLOCATE( PSO2_DMS  )
      IF ( ALLOCATED( PSO4_SO2  ) ) DEALLOCATE( PSO4_SO2  )
      IF ( ALLOCATED( SOx_SCALE ) ) DEALLOCATE( SOx_SCALE )
      IF ( ALLOCATED( SSTEMP    ) ) DEALLOCATE( SSTEMP    )
      IF ( ALLOCATED( TCOSZ     ) ) DEALLOCATE( TCOSZ     )
      IF ( ALLOCATED( TTDAY     ) ) DEALLOCATE( TTDAY     )          
      IF ( ALLOCATED( VCLDF     ) ) DEALLOCATE( VCLDF     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_SULFATE

!------------------------------------------------------------------------------

      END MODULE SULFATE_MOD
