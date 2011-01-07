! $Id: drydep_mod.f,v 1.5 2010/03/15 19:33:24 ccarouge Exp $
      MODULE DRYDEP_MOD
!
!******************************************************************************
!  Module DRYDEP_MOD contains variables and routines for the GEOS-CHEM dry
!  deposition scheme. (bmy, 1/27/03, 10/19/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) MAXDEP   (INTEGER) : Maximum number of drydep species
!  (2 ) NNTYPE   (INTEGER) : Max # of landtypes / grid box
!  (3 ) NNPOLY   (INTEGER) : Number of drydep polynomial coefficients
!  (4 ) NNVEGTYPE(INTEGER) : Number of Olson land types
!  (5 ) XCKMAN   (REAL*8 ) : Von Karman constant?
!  (6 ) DRYDHNO3 (INTEGER) : Internal flag for location of HNO3 in DEPVEL
!  (7 ) DRYDNO2  (INTEGER) : Internal flag for location of NO2  in DEPVEL
!  (8 ) DRYDPAN  (INTEGER) : Internal flag for location of PAN  in DEPVEL
!  (9 ) NUMDEP   (INTEGER) : Actual number of drydep species
!  (10) NWATER   (INTEGER) : Number of Olson's surface types that are water
!  (11) AIROSOL  (LOGICAL) : Array flags to denote aerosol drydep species
!  (12) IDEP     (INTEGER) : ID #'s for dry deposition surface types 
!  (13) IRAC     (INTEGER) : ???       resistance for drydep land type
!  (14) IRCLO    (INTEGER) : ???       resistance for drydep land type
!  (15) IRCLS    (INTEGER) : ???       resistance for drydep land type
!  (16) IRGSO    (INTEGER) : ???       resistance for drydep land type
!  (17) IRGSS    (INTEGER) : ???       resistance for drydep land type
!  (18) IRI      (INTEGER) : Internal  resistance for drydep land types
!  (19) IRLU     (INTEGER) : Cuticular resistance for drydep land types
!  (20) IVSMAX   (INTEGER) : ???       resistance for drydep land type
!  (21) IWATER   (INTEGER) : ID #'s for Olson surface types that are water 
!  (22) IZO      (INTEGER) : Roughness heights for each Olson surface type
!  (23) NDVZIND  (INTEGER) : Index array for ordering drydep species in DEPVEL
!  (24) NTRAIND  (INTEGER) : Stores tracer numbers of drydep species
!  (25) DEPSAV   (REAL*8 ) : Array containing dry deposition frequencies [s-1]
!  (26) PBLFRAC  (REAL*8 ) : Array for multiplicative factor for drydep freq
!  (27) DRYCOEFF (REAL*8 ) : Polynomial coefficients for dry deposition
!  (28) HSTAR    (REAL*8 ) : Henry's law constant
!  (29) F0       (REAL*8 ) : Reactivity factor for biological oxidation
!  (30) XMW      (REAL*8 ) : Molecular weight of drydep species [kg]
!  (32) A_RADI   (REAL*8 ) : Radius of aerosol for size-resolved drydep [um]
!  (33) A_DEN    (REAL*8 ) : Density of aerosol for size-res'd drydep [kg/m3]
!  (33) DEPNAME  (CHAR*14) : Names of dry deposition species
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_DRYDEP          : Dry deposition driver routine
!  (2 ) DVZ_MINVAL         : Sets minimum drydep velocities for SULFATE tracers
!  (3 ) METERO             : Computes meterological fields for dry deposition
!  (4 ) DRYFLX             : Applies drydep losses from SMVGEAR to tracer array
!  (5 ) DRYFLXRnPbBe       : Applies drydep losses to 210Pb and 7Be 
!  (6 ) DRYFLXH2HD         : Applies drydep losses to H2 and HD
!  (7 ) DEPVEL             : Computes dry deposition velocities (by D. Jacob)
!  (8 ) DIFFG              : Computes diffusion coefficient for a gas
!  (9 ) MODIN              : Reads inputs for DEPVEL from "drydep.table"
!  (10) RDDRYCF            : Reads drydep polynomial coeffs from "drydep.coef"
!  (11) AERO_SFCRSI        : Computes dust sfc resistance ff Seinfeld et al 86
!  (12) AERO_SFCRSII       : Conputes dust sfc resistance ff Zhang et al 2001
!  (13) INIT_DRYDEP        : Initializes and allocates module arrays
!  (14) CLEANUP_DRYDEP     : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "drydep_mod.f":
!  ============================================================================
!  (1 ) comode_mod.f       : Module w/ SMVGEAR allocatable arrays
!  (2 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f    : Module w/ GEOS-CHEM data & met field dirs
!  (4 ) error_mod.f        : Module w/ NaN, other error check routines
!  (5 ) file_mod.f         : Module w/ file unit #'s and error checks
!  (6 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (7 ) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (8 ) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (9 ) tracer_mod.f       : Module w/ GEOS-CHEM tracer array etc.
!  (10) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  References:
!  ============================================================================
!  (1 ) Baldocchi, D.D., B.B. Hicks, and P. Camara, "A canopy stomatal
!        resistance model for gaseous deposition to vegetated surfaces",
!        Atmos. Environ. 21, 91-101, 1987.
!  (2 ) Brutsaert, W., "Evaporation into the Atmosphere", Reidel, 1982.
!  (3 ) Businger, J.A., et al., "Flux-profile relationships in the atmospheric 
!        surface layer", J. Atmos. Sci., 28, 181-189, 1971.
!  (4 ) Dwight, H.B., "Tables of integrals and other mathematical data",
!        MacMillan, 1957.
!  (5 ) Guenther, A., and 15 others, A global model of natural volatile
!         organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!  (6 ) Hicks, B.B., and P.S. Liss, "Transfer of SO2 and other reactive
!        gases across the air-sea interface", Tellus, 28, 348-354, 1976.
!  (7 ) Jacob, D.J., and S.C. Wofsy, "Budgets of reactive nitrogen,
!        hydrocarbons, and ozone over the Amazon forest during the wet season",
!        J.  Geophys. Res., 95, 16737-16754, 1990.
!  (8 ) Jacob, D.J., et al, "Deposition of ozone to tundra", J. Geophys. Res., 
!        97, 16473-16479, 1992.
!  (9 ) Levine, I.N., "Physical Chemistry, 3rd ed.", McGraw-Hill, 
!        New York, 1988.
!  (10) Munger, J.W., et al, "Atmospheric deposition of reactive nitrogen 
!        oxides and ozone in a temperate deciduous forest and a sub-arctic 
!        woodland", J. Geophys. Res., in press, 1996.
!  (11) Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, "SO2, sulfate, 
!        and HNO3 deposition velocities computed using regional landuse and
!        meteorological data", Atmos. Environ., 20, 949-964, 1986.
!  (12) Wang, Y.H., paper in preparation, 1996.
!  (13) Wesely, M.L, "Improved parameterizations for surface resistance to
!        gaseous dry deposition in regional-scale numerical models", 
!        Environmental Protection Agency Report EPA/600/3-88/025,
!        Research Triangle Park (NC), 1988.
!  (14) Wesely, M. L., Parameterization of surface resistance to gaseous dry 
!        deposition in regional-scale numerical models.  Atmos. Environ., 23
!        1293-1304, 1989. 
!  (15) Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!        Global Budget of Molecular Hydrogen and its Deuterium Content: 
!        Constraints from Ground Station, Cruise, and Aircraft Observations,
!        submitted to J. Geophys. Res., 2007.
!!
!  NOTES:
!  (1 ) Bug fix: Do not assume NO2 is the 2nd drydep species.  This causes
!        a mis-indexing for CANOPYNOX.  Now archive ND44 diagnostic in kg for
!        Radon runs in routine DRYFLXRnPbBe; convert to kg/s in diag3.f
!        (bmy, 1/27/03)
!  (2 ) Now references "grid_mod.f" and the new "time_mod.f".  Renamed DRYDEP
!        routine to DO_DRYDEP for consistency w/ other drivers called from
!        the MAIN program. (bmy, 2/11/03)
!  (3 ) Added error check in DRYFLX for SMVGEAR II (bmy, 4/28/03)
!  (4 ) Added drydep of N2O5.  Now added PBLFRAC array, which is the fraction
!        of each level below the PBL top.  Also now compute drydep throughout 
!        the entire PBL, in order to prevent short-lived species such as HNO3 
!        from being depleted in the shallow GEOS-3 surface layer.  
!        (rjp, bmy, 7/21/03)
!  (5 ) Bug fix for GEOS-4 in DRYFLXRnPbBe (bmy, 12/2/03)
!  (6 ) Now made CFRAC, RADIAT local variables in DO_DRYDEP (bmy, 12/9/03)
!  (7 ) Now enclose AD44 in !$OMP CRITICAL block for drydep flux (bmy, 3/24/04)
!  (8 ) Now handle extra carbon & dust tracers (rjp, tdf, bmy, 4/1/04)
!  (9 ) Added routines AERO_SFCRS1, AERO_SFCRSII.  Increased MAXDEP to 25.
!        Now handles extra carbon & dust tracers. (rjp, tdf, bmy, 4/1/04)
!  (10) Increased MAXDEP to 26.  Added A_RADI and A_DEN module variables.
!        Other modifications for size-resolved drydep. (rjp, bec, bmy, 4/20/04)
!  (11) Increased MAXDEP to 35 and handle extra SOA tracers (rjp, bmy, 7/13/04)
!  (12) Now references "logical_mod.f", "directory_mod.f", and "tracer_mod.f"
!        (bmy, 7/20/04)
!  (13) Add Hg2, HgP as drydep tracers (eck, bmy, 12/8/04)
!  (14) Updated for AS, AHS, LET, NH4aq, SO4aq (cas, bmy, 1/6/05)
!  (15) Now references "pbl_mix_mod.f".  Removed PBLFRAC array. (bmy, 2/22/05)
!  (16) Now include SO4s, NITs tracers.  Now accounts for hygroscopic growth
!        of seasalt aerosols when computing aerodynamic resistances.
!        (bec, bmy, 4/13/05)
!  (17) Now modified for GEOS-5 and GCAP met fields (bmy, 5/25/05)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Now change Reynold's # criterion from 1 to 0.1 in DEPVEL.  Also 
!        change Henry's law constant for Hg2.  Also increase MAXDEP from
!        35 to 37. (eck, djj, bmy, 2/1/06)
!  (20) Bug fix in INIT_DRYDEP (bmy, 4/17/06)
!  (21) Now bundle function DIFFG into "drydep_mod.f".  Also updated for SOG4
!        and SOA4 tracers.  Bug fix in INIT_DRYDEP. (dkh, bmy, 5/24/06)
!  (22) Fix typo in INIT_DRYDEP (dkh, bmy, 6/23/06)
!  (23) Add H2 and HD as drydep tracers. Added subroutine DRYFLXH2HD for H2HD
!        offline sim (phs, 9/18/07)
!  (24) Extra error check for small RH in AERO_SFCRII (phs, 6/11/08)
!  (25) Added 15 more dry deposition species (tmf, 7/31/08)
!  (26) Modify dry depostion to follow the non-local PBL scheme.
!        (lin, ccc, 5/29/09)
!  (27) Minor bug fix in mol wts for ALPH, LIMO (bmy, 10/19/09)
!  (28) modified to use Zhang 2001 for all non-size resolved aerosols (hotp)
!  (29) Add aromatics SOA (dkh)
!  (30) Add new species. Some tracers give 2 deposition species: ISOPN-> ISOPNB
!       and ISOPND. (fp)
!  (31) Updates for mercury simulation (ccc, 5/17/10)
!  (32) Add POPs (eck, 9/20/10)
!******************************************************************************
!
      USE LOGICAL_MOD,     ONLY : LNLPBL ! (Lin, 03/31/09)

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      PUBLIC :: DEPNAME
      PUBLIC :: DEPSAV
      PUBLIC :: MAXDEP
      PUBLIC :: NUMDEP
      PUBLIC :: NTRAIND
      PUBLIC :: DRYHg0, DRYHg2, DryHgP !CDH
      PUBLIC :: DRYPOPG, DRYPOPP
      
      ! ... and these routines
      PUBLIC :: CLEANUP_DRYDEP     
      PUBLIC :: DO_DRYDEP
      PUBLIC :: DRYFLX   
      PUBLIC :: DRYFLXH2HD
      PUBLIC :: DRYFLXRnPbBe       
      PUBLIC :: DVZ_MINVAL         
      PUBLIC :: INIT_DRYDEP        

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER   :: MAXDEP    = 50
      INTEGER, PARAMETER   :: NNTYPE    = 15     ! NTYPE    from "CMN_SIZE"
      INTEGER, PARAMETER   :: NNPOLY    = 20     ! NPOLY    from "CMN_SIZE"
      INTEGER, PARAMETER   :: NNVEGTYPE = 74     ! NVEGTYPE from "CMN_SIZE"
      REAL*8,  PARAMETER   :: XCKMAN    = 0.4d0
 
      ! Scalars
      INTEGER              :: DRYDHNO3, DRYDNO2, DRYDPAN
      !FP_ISOP (6/2009)
      INTEGER              :: DRYDH2O2
      INTEGER              :: NUMDEP,   NWATER

      ! Arrays
      LOGICAL              :: AIROSOL(MAXDEP)
      INTEGER              :: IDEP(NNVEGTYPE)
      INTEGER              :: IRAC(NNTYPE)
      INTEGER              :: IRCLO(NNTYPE)
      INTEGER              :: IRCLS(NNTYPE)
      INTEGER              :: IRGSS(NNTYPE)
      INTEGER              :: IRGSO(NNTYPE)
      INTEGER              :: IRI(NNTYPE)
      INTEGER              :: IRLU(NNTYPE)
      INTEGER              :: IVSMAX(NNTYPE)
      INTEGER              :: IZO(NNVEGTYPE)
      INTEGER              :: IWATER(NNVEGTYPE)
      INTEGER              :: NDVZIND(MAXDEP)
      INTEGER              :: NTRAIND(MAXDEP)
      REAL*8,  ALLOCATABLE :: DEPSAV(:,:,:)
      REAL*8               :: DRYCOEFF(NNPOLY)
      REAL*8               :: HSTAR(MAXDEP)
      REAL*8               :: KOA(MAXDEP)
      REAL*8               :: F0(MAXDEP)
      REAL*8               :: XMW(MAXDEP)
      REAL*8               :: A_RADI(MAXDEP)
      REAL*8               :: A_DEN(MAXDEP)
      CHARACTER(LEN=14)    :: DEPNAME(MAXDEP)

      ! Additional variables for mercury sim (cdh, 9/1/09)
      INTEGER              :: DRYHg0, DRYHg2, DryHgP

      ! Additional variables for POPs sim (eck, 9/20/10)
      INTEGER              :: DRYPOPP, DRYPOPG

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_DRYDEP
!
!******************************************************************************
!  Subroutine DO_DRYDEP is the driver for the GEOS-CHEM dry deposition scheme.
!  DO_DRYDEP calls DEPVEL to compute deposition velocities [m/s], which are 
!  then converted to [cm/s].  Drydep frequencies are also computed.  
!  (lwh, gmg, djj, 1989, 1994; bmy, 2/11/03, 5/25/05)
!
!  DAO met fields passed via "dao_mod.f":
!  ============================================================================
!  (1 ) AD      (REAL*8 ) : Array for dry air mass at each grid box [kg]
!  (2 ) AZO     (REAL*8 ) : Array for surface roughness heights     [m]
!  (3 ) ALBD    (REAL*8 ) : Array for surface albedo                [m]
!  (5 ) SUNCOS  (REAL*8 ) : Array for COSINE( solar zenith angle )  [unitless]
!  (6 ) T       (REAL*8 ) : Array for grid box temperature          [K]
!  (7 ) USTAR   (REAL*8 ) : Array for grid box friction velocity    [m/s]
!
!  Other important quantities:
!  ============================================================================
!  (1 ) LSNOW   (LOGICAL) : Array to flag whether there is snow/ice on the sfc.
!  (2 ) CZ1     (REAL*8 ) : Midpoint height of first model level     [m]
!  (3 ) OBK     (REAL*8 ) : Array for Monin-Obhukov Length           [m]
!  (4 ) TC0     (REAL*8 ) : Array for grid box surface temperature   [K]
!  (5 ) ZH      (REAL*8 ) : Array for PBL heights at each grid box   [m]
!  (6 ) DVEL    (REAL*8 ) : Array containing drydep velocities       [m/s]
!  (7 ) CFRAC   (REAL*8 ) : Array containing column cloud frac       [unitless]
!  (8 ) RADIAT  (REAL*8 ) : Array containing solar radiation         [W/m2]
!  (9 ) RHB     (REAL*8 ) : Array containing relative humidity       [unitless]
!
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989 
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
!  NOTES:
!  (1 ) Remove SUNCOS, USTAR, AZO, OBK from the arg list; now reference these
!        as well as AD and T from "dao_mod.f".  Cleaned up code and updated
!        comments.  Now only order tracer numbers into NTRAIND on the first
!        call.  Now force double-precision with "D" exponents.  Now also 
!        reference IDTNOX, IDTOX, etc. from "tracerid_mod.f".  Bundled into
!        "drydep_mod.f" (bmy, 11/19/02) 
!  (2 ) Now make sure that the PBL depth (THIK) is greater than or equal to 
!        the thickness of the first layer.  Now initialize PBLFRAC array on
!        each call. (rjp, bmy, 7/21/03)
!  (3 ) Now declare CFRAC, RADIAT, AZO, USTAR as local variables, which are 
!        returned by METERO.  CFRAC and RADIAT have also been deleted from 
!        "CMN_DEP". (bmy, 12/9/03)
!  (4 ) Now use explicit formula for IJLOOP to allow parallelization.
!        Also reference LPRT from "logical_mod.f" (bmy, 7/20/04)
!  (5 ) Now use routines from "pbl_mix_mod.f" to get PBL quantities, instead
!        of re-computing them here.  Removed PBLFRAC array.  Removed reference
!        to "pressure_mod.f".  Removed reference to header file CMN.
!        Parallelize DO-loops. (bmy, 2/22/05)
!  (6 ) Now define RHB as a local array, which is defined in METERO and then
!        passed to DEPVEL. (bec, bmy, 4/13/05)
!  (7 ) Now dimension AZO for GEOS or GCAP met fields.  Remove obsolete
!        variables. (swu, bmy, 5/25/05)
!  (8 ) Remove reference to TRACERID_MOD, it's not needed (bmy, 10/3/05)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DAO_MOD,      ONLY : AD, ALBD, BXHEIGHT, SUNCOS
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! ND44
#     include "CMN_DEP"  ! IREG, ILAND, IUSE, etc.
#     include "CMN_GCTM" ! Physical constants

      ! Local variables
      LOGICAL, SAVE     :: FIRST = .TRUE.
      LOGICAL           :: LSNOW(MAXIJ)
      INTEGER           :: I, J, L, N, IJLOOP, NN, NDVZ
      REAL*8            :: THIK,              DVZ
      REAL*8            :: CZ1(MAXIJ),        TC0(MAXIJ)
      REAL*8            :: ZH(MAXIJ),         OBK(MAXIJ)
      REAL*8            :: CFRAC(MAXIJ),      RADIAT(MAXIJ)
      REAL*8            :: USTAR(MAXIJ),      RHB(MAXIJ)
      REAL*8            :: DVEL(MAXIJ,MAXDEP)

      ! Dimension AZO for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
      REAL*8            :: AZO(NTYPE)
#else
      REAL*8            :: AZO(MAXIJ)
#endif

      !=================================================================
      ! DO_DRYDEP begins here!
      !=================================================================
      
      ! Read drydep coeff's and land types on first call
      IF ( FIRST ) THEN
         CALL RDDRYCF
         CALL MODIN
         FIRST = .FALSE.
      ENDIF
 
      ! Call METERO to obtain meterological fields (all 1-D arrays)
      CALL METERO( CZ1, TC0,   OBK, CFRAC, RADIAT, 
     &             AZO, USTAR, ZH,  LSNOW, RHB )

      !=================================================================
      ! Call DEPVEL to compute dry deposition velocities [m/s]
      !=================================================================
      CALL DEPVEL( MAXIJ, RADIAT,  TC0,   SUNCOS, F0,  HSTAR, 
     &             XMW,   AIROSOL, USTAR, CZ1,    OBK, CFRAC,  
     &             ZH,    LSNOW,   DVEL,  AZO,    RHB ) 

      !=================================================================
      ! Compute dry deposition frequencies; archive diagnostics
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, THIK, N, NN, NDVZ, DVZ )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D grid box index
         IJLOOP  = ( (J-1) * IIPAR ) + I

         ! THIK = thickness of surface layer [m]
         THIK    = BXHEIGHT(I,J,1)

         ! Now we calculate drydep throughout the entire PBL.  
         ! Make sure that the PBL depth is greater than or equal 
         ! to the thickness of the 1st layer (rjp, bmy, 7/21/03)
         ! Add option for non-local PBL mixing scheme: THIK must
         ! be the first box height. (Lin, 03/31/09) 
         IF (.NOT. LNLPBL) THIK    = MAX( ZH(IJLOOP), THIK )

         ! Loop over drydep species
         DO N = 1, NUMDEP

            ! GEOS-CHEM tracer number
            NN   = NTRAIND(N)

            ! Index of drydep species in the DVEL array 
            ! as passed back from subroutine DEPVEL
            NDVZ = NDVZIND(N)
                   
            ! Dry deposition velocity [cm/s]
            DVZ  = DVEL(IJLOOP,NDVZ) * 100.d0

            ! Set minimum velocity for sulfate tracers
            DVZ  = DVZ_MINVAL( NN, LSNOW(IJLOOP), DVZ )

            ! Dry deposition frequency [1/s]
            DEPSAV(I,J,N) = ( DVZ / 100.d0 ) / THIK

            ! ND44 diagnostic: drydep velocity [cm/s]
            IF ( ND44 > 0 ) THEN 
               AD44(I,J,N,2) = AD44(I,J,N,2) + DVZ
            ENDIF
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )

      ! Return to calling program
      END SUBROUTINE DO_DRYDEP

!------------------------------------------------------------------------------

      FUNCTION DVZ_MINVAL( N, LSNOW, DVZ ) RESULT( NEWDVZ )
!
!******************************************************************************
!  Function DVZ_MINVAL sets minimum values for drydep velocities for 
!  SULFATE TRACERS, according to Mian Chin's GOCART model. 
!  (rjp, bmy, 11/21/02, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N     (INTEGER) : Tracer number
!  (2 ) LSNOW (LOGICAL) : Flag for denoting snow/ice 
!  (3 ) DVZ   (REAL*8 ) : Deposition velocity [cm/s]
!
!  NOTES:
!  (1 ) Don't put a min drydep value on H2O2 for offline run (rjp, bmy,3/31/03)
!  (2 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTMSA, IDTNH3, IDTNH4
      USE TRACERID_MOD, ONLY : IDTNIT, IDTSO2, IDTSO4

#     include "CMN_SIZE"   ! Size parameters!

      ! Arguments
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: LSNOW
      REAL*8,  INTENT(IN) :: DVZ

      ! Function value
      REAL*8              :: NEWDVZ

      !=================================================================
      ! DVZ_MINVAL begins here!
      !=================================================================

      !---------------------------------------
      ! SO2, NH3, offline H2O2
      ! Min Vd = 2.0e-1 [cm/s] over ice/snow
      !        = 3.0e-1 [cm/s] over land
      !---------------------------------------
      IF ( N == IDTSO2 .or. N == IDTNH3 ) THEN 

         IF ( LSNOW ) THEN
            NEWDVZ = MAX( DVZ, 2.0d-1 )       
         ELSE
            NEWDVZ = MAX( DVZ, 3.0d-1 )
         ENDIF

      !---------------------------------------
      ! SO4, MSA, NH4, NIT
      ! Min Vd = 1.0e-2 [cm/s] 
      !---------------------------------------         
      ELSE IF ( N == IDTSO4 .or. N == IDTMSA  .or.
     &          N == IDTNH4 .or. N == IDTNIT ) THEN

         NEWDVZ = MAX( DVZ, 1.0d-2 )

      !---------------------------------------
      ! Other drydep species: do nothing
      !---------------------------------------         
      ELSE
         NEWDVZ = DVZ

      ENDIF

      ! Return to calling program
      END FUNCTION DVZ_MINVAL

!------------------------------------------------------------------------------

      SUBROUTINE METERO( CZ1, TC0,  OBK, CFRAC, RADIAT, 
     &                   AZO, USTR, ZH,  LSNOW, RHB )
!
!******************************************************************************
!  Subroutine METERO calculates meteorological constants needed for the      
!  dry deposition velocity module. (lwh, gmg, djj, 1989, 1994; bmy, 10/3/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) CZ1    (REAL*8 ) : Midpoint height of first model level     [m]
!  (2 ) TC0    (REAL*8 ) : Array for grid box surface temperature   [K]
!  (3 ) OBK    (REAL*8 ) : Array for the Monin-Obhukov length       [m]
!  (4 ) CFRAC  (REAL*8 ) : Array for the column cloud fraction      [unitless]
!  (5 ) RADIAT (REAL*8 ) : Array for the solar radiation @ ground   [W/m2]
!  (6 ) AZO    (REAL*8 ) : Array for the roughness heights          [m]
!  (7 ) USTR   (REAL*8 ) : Array for the friction velocity          [m/s]
!  (8 ) ZH     (REAL*8 ) : Height of the mixed layer (aka PBL)      [m]
!  (9 ) LSNOW  (LOGICAL) : Flag to denote ice & snow (ALBEDO < 0.4)
!  (10) RHB    (REAL*8 ) : Relative humidity at surface             [unitless]
!
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989. 
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
!  NOTES: 
!  (1 ) Now reference GET_PEDGE from "pressure_mod.f".  Now reference T from 
!        "dao_mod.f".  Removed obsolete code & comments, and added new 
!         documentation header.  Now force double precision with "D" 
!         exponents.  Now compute OBK here as well.  Bundled into F90 module
!         "drydep_mod.f" (bmy, 11/20/02)
!  (2 ) Now reference CLDFRC, RADSWG, ZO, USTAR from "dao_mod.f".  Also now 
!         pass CFRAC, RADIAT, AZO, USTR back to the calling routine 
!         via the arg list. (bmy, 12/9/03)
!  (3 ) Now use explicit formula for IJLOOP to allow parallelization
!        (bmy, 7/20/04)
!  (4 ) Now compute ZH and LSNOW here instead of w/in DO_DRYDEP.  Parallelize
!        DO-loops.  Now use BXHEIGHT from "dao_mod.f" instead of computing 
!        the thickness of the 1st level here.  Remove reference to 
!        "pressure_mod.f".  Remove reference to T from "dao_mod.f".  Now
!        reference ALBD from "dao_mod.f" (bmy, 2/22/05)
!  (5 ) Now references RH from "dao_mod.f".  Now passes relative humidity
!        from the surface layer back via RHB argument. (bec, bmy, 4/13/05)
!  (6 ) Now call GET_OBK from "dao_mod.f" to get the M-O length for both
!        GEOS or GCAP met fields.  Remove local computation of M-O length
!        here.  Also now dimension AZO appropriately for GCAP or GEOS met
!        fields.  Remove obsolete variables. (swu, bmy, 5/25/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Move XLTMMP function to module MEGANUT_MOD. (ccc, 11/20/09)
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : ALBD,   BXHEIGHT, CLDFRC, GET_OBK
      USE DAO_MOD,      ONLY : RADSWG, RH,       TS,     USTAR,   Z0
      USE PBL_MIX_MOD,  ONLY : GET_PBL_TOP_m
      USE MEGANUT_MOD,  ONLY : XLTMMP
                                  
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! Physical constants

      ! Arguments
      LOGICAL, INTENT(OUT)  :: LSNOW(MAXIJ)
      REAL*8,  INTENT(OUT)  :: CZ1(MAXIJ)
      REAL*8,  INTENT(OUT)  :: TC0(MAXIJ)
      REAL*8,  INTENT(OUT)  :: OBK(MAXIJ)
      REAL*8,  INTENT(OUT)  :: CFRAC(MAXIJ)
      REAL*8,  INTENT(OUT)  :: RADIAT(MAXIJ)
      REAL*8,  INTENT(OUT)  :: RHB(MAXIJ)
      REAL*8,  INTENT(OUT)  :: USTR(MAXIJ)
      REAL*8,  INTENT(OUT)  :: ZH(MAXIJ)

      ! Dimension AZO for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
      REAL*8, INTENT(OUT)   :: AZO(NTYPE)
#else
      REAL*8, INTENT(OUT)   :: AZO(MAXIJ)
#endif

      ! Local variables
      INTEGER               :: I,  J,  IJLOOP
      REAL*8                :: THIK

      ! External functions
!-- Moved to megan_mod.f. (ccc, 11/20/09)
!      REAL*8, EXTERNAL      :: XLTMMP

      !=================================================================
      ! METERO begins here!
      !=================================================================

#if   defined( GCAP )
      ! For GCAP: AZO (roughness ht) is a function of Olson land type 
      ! instead of lat/lon location.  Zero AZO here; AZO will be 
      ! computed internally w/in routine DEPVEL (swu, bmy, 5/25/05)
      AZO(:) = 0d0
#endif

      ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, IJLOOP, THIK )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D grid box index
         IJLOOP      = ( (J-1) * IIPAR ) + I

         ! THIK = thickness of layer 1 [m]
         THIK        = BXHEIGHT(I,J,1)

         ! Midpoint height of first model level [m]
         CZ1(IJLOOP) = THIK / 2.0d0

         !==============================================================
         ! Return meterological quantities as 1-D arrays for DEPVEL
         !==============================================================

#if   !defined( GCAP )
         ! For GEOS: Roughness height [m] is a function of lat/lon
         AZO(IJLOOP)    = Z0(I,J)
#endif

         ! Column cloud fraction [unitless]
         CFRAC(IJLOOP)  = CLDFRC(I,J)

         ! Set logical LSNOW if snow and sea ice (ALBEDO > 0.4)
         LSNOW(IJLOOP)  = ( ALBD(I,J) > 0.4 )

         ! Monin-Obhukov length [m]
         OBK(IJLOOP)    = GET_OBK( I, J )

         ! Solar insolation @ ground [W/m2]
         RADIAT(IJLOOP) = RADSWG(I,J) 

         ! Surface temperature [K]
         TC0(IJLOOP)    = TS(I,J)
         
         ! Friction velocity [m/s]
         USTR(IJLOOP)   = USTAR(I,J)

         ! Mixed layer depth [m]
         ZH(IJLOOP)     = GET_PBL_TOP_m( I, J )

         ! Relative humidity @ surface [unitless] (bec, bmy, 4/13/05)
         RHB(IJLOOP)    = MIN( 0.99d0, RH(I,J,1) * 1.d-2 ) 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE METERO

!------------------------------------------------------------------------------

      SUBROUTINE DRYFLX
!
!******************************************************************************
!  Subroutine DRYFLX sets up the dry deposition flux diagnostic for tracers
!  which are part of the SMVGEAR mechanism. (bmy, bdf, 4/20/99, 3/24/04)
!
!  NOTES:
!  (1 ) Bug fix -- now skip tracers for which NTDEP(N) is zero, in order
!        to avoid array-out-of-bounds errors. (bmy, 5/2/00)
!  (2 ) Now reference the CSPEC array from "comode_mod.f" instead of from
!        common block header "comode.h". (bmy, 7/11/00)
!  (3 ) Also reference JLOP and VOLUME from "comode_mod.f" (bmy, 10/19/00)
!  (4 ) Updated comments, cosmetic changes (bmy, 3/14/02)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Removed reference to "comtrid.h", "CMN_SAV", "CMN_DEP", and "CMN_O3",
!        these are not used in this routine.  Also bundled into "drydep_mod.f"
!        for more convenient packaging. (bmy, 11/19/02)
!  (7 ) Replaced DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid"mod.f".
!        Also removed references to JREF and FLUXRUL.  Now use function
!        GET_TS_CHEM from "time_mod.f". (bmy, 2/11/03)
!  (8 ) Now references ERROR_STOP from "error_mod.f" (bmy, 4/28/03)
!  (9 ) Now sum drydep fluxes throughout the entire PBL.  Added L variable.
!        AREA_CM2 has now been made into a lookup table. Now implement a 
!        parallel DO loop for efficiency. (rjp, bmy, 7/21/03)
!  (10) Now bracket AD44 with a !$OMP CRITICAL block in order to avoid
!        multiple threads writing to the same element (bmy, 3/24/04)
!  (11) Now reference GET_FRAC_UNDER_PBLTOP and GET_PBL_MAX_L from 
!        "pbl_mix_mod.f".  Remove reference to CMN. (bmy, 2/22/05)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,  ONLY : CSPEC, JLOP, VOLUME
      USE DIAG_MOD,    ONLY : AD44
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD, ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,    ONLY : GET_TS_CHEM

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! Diagnostic switches & arrays
#     include "comode.h"   ! CSPEC

      ! Local variables
      INTEGER :: I, J, JJ, JLOOP, L, L_PBLTOP, N, NK, NN
      REAL*8  :: DTCHEM, PBL_MAX, TDRYFX, AREA_CM2(JJPAR)

      !=================================================================
      ! DRYFLX begins here!
      !=================================================================

      ! Return unless we have turned on ND44 drydep diagnostic
      IF ( ND44 == 0 ) RETURN

      ! There is only drydep in the surface layer, which
      ! is accounted for in the "URBAN" chemistry slot
      NCS     = NCSURBAN

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Highest extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L()

      !=================================================================
      ! ND44 diagnostic: Dry deposition flux [molec/cm2/s]
      !
      ! NOTE: DRYFLX will only archive the dry deposition fluxes for
      ! tracers which are SMVGEAR species.  Fluxes for sulfate tracers
      ! will be updated in "sulfate_mod.f". (bmy, 11/19/02)
      !=================================================================

      ! Save grid box surface area [cm2] in a lookup table (bmy, 7/23/03)
      DO J = 1, JJPAR
         AREA_CM2(J) = GET_AREA_CM2(J)
      ENDDO

      ! Loop over dry deposition species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NK, JJ, JLOOP, TDRYFX )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NUMDEP

         ! Index for drydep species #N, from SMVGEAR
         NK = NTDEP(N)

         ! If NK <= 0, then skip to the next tracer.  
         ! This avoids array-out-of-bounds errors (bmy, 5/2/00)
         IF ( NK <= 0 ) CYCLE

         ! Index for drydep flux in CSPEC array 
         JJ = IRM(NPRODLO+1,NK,NCS)

         ! Error check JJ -- can't be zero
         IF ( JJ <= 0 ) THEN 
            CALL ERROR_STOP( 'Drydep species mis-indexing!', 
     &                       'DRYFLX ("error_mod.f")' )
         ENDIF

         ! Loop over grid boxes
         DO L = 1, PBL_MAX
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         
            ! Only deal w/ boxes w/in the boundary layer
            IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

               ! 1-D grid box index for CSPEC & VOLUME
               JLOOP = JLOP(I,J,L)

               ! Dry dep flux [molec] for species N = 
               !  CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
               !  [molec/cm3]     * [cm3]  
               TDRYFX = CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
                    
               ! Convert TDRYFX from [molec] to [molec/cm2/s]        
               TDRYFX = TDRYFX / ( AREA_CM2(J) * DTCHEM ) 
                    
!$OMP CRITICAL
               ! Save into AD44 diagnostic array
               AD44(I,J,N,1) = AD44(I,J,N,1) + TDRYFX
!$OMP END CRITICAL

            ENDIF
         ENDDO           
         ENDDO               
         ENDDO  
      ENDDO                     
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRYFLX

!------------------------------------------------------------------------------

      SUBROUTINE DRYFLXRnPbBe
!
!******************************************************************************
!  Subroutine DRYFLXRnPbBe removes dry deposition losses from the STT tracer
!  array and archives deposition fluxes to the ND44 diagnostic. 
!  (hyl, bmy, bdf, 4/2/99, 5/25/05)
!
!  NOTES:
!  (1 ) Now eliminate DEPFLUX from CMN_SAV, in order to save memory.
!        DEPFLUX is now a local variable (bdf, 4/2/99)
!  (2 ) Now make DEPFLUX of dimension (IIPAR,JJPAR,MAXDEP) (bmy, 4/2/99)
!  (3 ) Now use an allocatable array for the ND44 diagnostic.
!        Also made cosmetic changes, updated comments. (bmy, 3/16/00)
!  (4 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (5 ) Added to module "RnPbBe_mod.f".  Also made cosmetic changes
!        and updated comments (bmy, 6/14/01)
!  (6 ) Updated comments (bmy, 3/29/02)
!  (7 ) Replace all instances of IM, JM, IMX, JMX, with IIPAR, JJPAR, IGLOB,
!        and JGLOB.  Now replaced DEPFLUX array w/ AMT_LOST scalar
!        variable.  Also make sure that the amount of tracer lost to drydep
!        is now accurately accounted in the ND44 diagnostic. (bmy, 8/7/02)
!  (8 ) Now call GEOS_CHEM_STOP or ERROR_STOP (from "error_mod.f") when 
!        stopping the run w/ an error condition. (bmy, 10/15/02)
!  (9 ) Now moved from "RnPbBe_mod.f" to "drydep_mod.f".  (bmy, 1/27/03)
!  (10) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 2/11/03)
!  (11) Now compute drydep fluxes throughout the entire PBL.  Now references
!        PBLFRAC.  Added L_PBLTOP variable. (bmy, 7/21/03)
!  (12) Now follow GEOS-3 algorithm for GEOS-4 model (bmy, 12/2/03)
!  (13) Now reference STT from "tracer_mod.f" and LDRYD from "logical_mod.f"
!        (bmy, 7/20/04)
!  (14) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,    ONLY : AD44
      USE ERROR_MOD,   ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE LOGICAL_MOD, ONLY : LDRYD
      USE PBL_MIX_MOD, ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,    ONLY : GET_TS_CHEM
      USE TRACER_MOD,  ONLY : STT

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND44
#     include "CMN_DEP"   ! Dry deposition variables

      ! Local variables
      INTEGER             :: I, J, L, PBL_MAX, N, NN
      REAL*8              :: DTCHEM, FRACLOST, F_UNDER_TOP, AMT_LOST

      !=================================================================
      ! DRYFLXRnPbBe begins here!!
      !=================================================================

      ! Return if drydep is turned off
      IF ( .not. LDRYD ) RETURN

      ! Chemistry timestep in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Maximum extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L() 

      ! Add option for non-local PBL mixing scheme: only done at the surface
      ! (Lin, 03/31/09) 
      IF (LNLPBL) PBL_MAX = 1

      ! Loop over drydep species
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NN, F_UNDER_TOP, FRACLOST, AMT_LOST )
      DO N = 1, NUMDEP

         ! Tracer index in STT that corresponds to drydep species N
         ! If invalid, then cycle
         NN = NTRAIND(N)
         IF ( NN == 0 ) CYCLE

         ! Loop over grid boxes
         DO L = 1, PBL_MAX
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            
            ! Fraction of box (I,J,L) under PBL top [unitless]
            F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

            ! FRACLOST is the fraction of tracer lost.  PBLFRAC is 
            ! the fraction of layer L located totally w/in the PBL.
            FRACLOST = DEPSAV(I,J,N) * F_UNDER_TOP * DTCHEM

            ! add option for non-local PBL mixing scheme: only the surface
            ! (Lin, 03/31/09) 
            IF (LNLPBL) FRACLOST = DEPSAV(I,J,N) * DTCHEM

            !===========================================================
            ! Proceed as follows:
            ! --------------------------------
            ! (a) If FRACLOST < 0, then stop the run.
            !
            ! (b) If FRACLOST > 1, use an exponential loss to 
            !     avoid negative tracer
            !
            ! (c) If FRACLOST is in the range (0-1), then use the
            !     the regular formula (STT * FRACLOST) to compute
            !     loss from dry deposition.
            !=====================================================

            ! Stop the run on negative FRACLOST!
            IF ( FRACLOST < 0 ) THEN
               CALL ERROR_STOP( 'FRACLOST < 0', 'dryflxRnPbBe' )
            ENDIF

            ! AMT_LOST = amount of tracer lost to drydep [kg]
            IF ( FRACLOST > 1 ) THEN
               AMT_LOST = STT(I,J,L,NN) * ( 1d0 - EXP(-FRACLOST) )
            ELSE
               AMT_LOST = STT(I,J,L,NN) * FRACLOST
            ENDIF

            ! ND44 diagnostic: drydep flux [kg/s]
            IF ( ND44 > 0 ) THEN
!$OMP CRITICAL
               AD44(I,J,N,1) = AD44(I,J,N,1) + ( AMT_LOST/DTCHEM ) 
!$OMP END CRITICAL
            ENDIF

            ! Subtract AMT_LOST from the STT array [kg]
            STT(I,J,L,NN)  = STT(I,J,L,NN) - AMT_LOST

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DRYFLXRnPbBe

!------------------------------------------------------------------------------

      SUBROUTINE DRYFLXH2HD
!
!******************************************************************************
!  Subroutine DRYFLXH2HD removes dry deposition losses from the tracer
!  array and archives deposition fluxes AND VELOCITY to the ND44 diagnostic. 
!  (adapted from DRYFLX v5-05, jaegle 11/02/2005). 
!
!  NOTES:
!  (1) Now deposit through the PBL. Commented but kept code related to soil
!       temperature (phs, 5/16/07)
!  (2) Move XLTMMP to module MEGANUT_MOD (ccc, 11/20/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE ERROR_MOD,    ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE GRID_MOD,     ONLY : GET_AREA_CM2, GET_XOFFSET, GET_YOFFSET
      USE DAO_MOD,      ONLY : T, TS, ALBD
      USE TRACER_MOD,   ONLY : STT
      USE LOGICAL_MOD,  ONLY : LDRYD
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE PBL_MIX_MOD,  ONLY : GET_PBL_TOP_m
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE MEGANUT_MOD,  ONLY : XLTMMP
      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! Diagnostic switches & arrays
#     include "CMN_VEL"   ! IJLAND
#     include "CMN_DEP"   ! Dry deposition variables
#     include "commsoil.h" ! Soil pulsing & wetness variables

      ! Local variables
      INTEGER     :: I, J, L, N, NN, M, PBL_MAX
      INTEGER     :: IJLOOP, I0, J0, IREF, JREF, K, STYP
      INTEGER     :: JLOP(IIPAR,JJPAR,1), NTYP(IIPAR,JJPAR)
      REAL*8      :: DTCHEM, FRACLOST, AMT_LOST
      REAL*8      :: THIK, DRYF, SVEL, FSOIL, AREA_CM2
      REAL*8      :: SOIL_H2, SOIL_HD, TMMP, STEMP(IIPAR,JJPAR)
      REAL*8      :: MLD
      REAL*8      :: F_UNDER_TOP

      ! External functions, for calculating soil temperature
!-- XLTMMP moved to meganut_mod.f (ccc, 11/20/09)
!      REAL*8, EXTERNAL   :: SOILTEMP, XLTMMP
      REAL*8, EXTERNAL   :: SOILTEMP

      !=================================================================
      ! DRYFLXH2HD begins here!!
      !=================================================================

      ! Chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Call soiltype to determine whether soil is dry or
      ! wet for all land grid-boxes
      CALL SOILTYPE

      ! Only do the following if DRYDEP is turned on
      IF ( .not. LDRYD ) RETURN

      ! Maximum extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L() 

      ! Add option for non-local PBL mixing scheme: only the surface
      ! (Lin, 03/31/09) 
      IF (LNLPBL) PBL_MAX = 1

      ! Need nested-grid offsets for soiltemp code
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()

      ! Initalize
      IJLOOP = 0
      DO J = 1, JJPAR
         DO I = 1, IIPAR
            IJLOOP        = IJLOOP + 1
            JLOP(I,J,1)   = IJLOOP
         ENDDO
      ENDDO


      ! Loop over drydep species
      DO N = 1, NUMDEP

         ! Tracer index in STT that corresponds to drydep species N
         ! If invalid, then cycle
         NN = NTRAIND(N)
         IF ( NN == 0 ) CYCLE

         ! Loop over layers (most efficient if moved below?)
         DO L = 1, PBL_MAX

           ! reset STEMP (could be a scalar -depends on future usage-)
           ! STEMP = 0  ! Not use yet

           ! Loop over each land grid-box
            DO M      = 1, NLAND
               IREF   = INDEXSOIL(1,M)
               JREF   = INDEXSOIL(2,M)
               I      = IREF - I0
               J      = JREF - J0
               IJLOOP = JLOP(I,J,1)

               ! Fraction of grid box that is ocean
               FSOIL = FRCLND(I,J)
 

               ! Only apply dry deposition over land surfaces which
                    ! are not covered with ice or desert (albedo < 0.4)
               ! and if we are in the window simulation.
               IF ( (I.GE.1)       .AND. (I.LE.IIPAR)  .AND.
     &              (J.GE.1)       .AND. (J.LE.JJPAR)  .AND.
     &              (FSOIL > 0.d0) .AND. (ALBD(I,J) < 0.4d0) ) THEN

                  ! Grid box area in cm2
                  AREA_CM2 = GET_AREA_CM2( J )

!                 !===========================================================
!                 ! Get SOIL TEMPerature from function SOILTEMP(I,J,M,NTYP)
!                 ! Right now Surface Temp is used instead.
!                 ! So commented for now (phs, 3/5/07)  
!                 !===========================================================
!                 TMMP = XLTMMP(I,J,IJLOOP) - 273.15
!  
!                 ! Loop over landtype
!                 DO K = 1, IREG(IREF,JREF)
!                      
!                    ! NCONSOIL Converts from Olson type  -> soil type           
!                    STYP = NCONSOIL(ILAND(IREF,JREF,K)+1)
!                          
!                    ! Temperature factor
!                    ! STEMP(I,J) is the weighted soil temperature in
!                    ! gridbox i,j for all the soil types
!                    ! IUSE is the fraction ((per mil) of box covered by land types
!                    STEMP(I,J) = STEMP(I,J) + 
!     &                           SOILTEMP(I,J,M,STYP,TMMP)*
!     &                           DBLE(IUSE(IREF,JREF,K))/1000.D0
!    
!                      !write(*,*)'TEst',TMMP, TS(I,J)-273.15, STEMP(I,J)
!                 ENDDO
!

                  ! SVEL [cm/s] is the air-to-soil transfer velocity
                  ! Use uniform value of 3.94d-2 cm/s over land
                  ! not covered by snow or desert.
                  SVEL = 3.94d-2
		  
            	
                  ! if soil temperature is below freezing reduce dep vel
                  ! by 1/2, and additional 1/2 below -15C(hup, 6/21/2005)
                  ! for now use surface temperature (TS) instead of
                  ! air temperature (T) jaegle, 12/12/2005
		  IF (TS(I,J) <= 273.15d0) SVEL = SVEL / 2.0d0
		  IF (TS(I,J) <= 258.15d0) SVEL =  SVEL / 2.0d0
                  
		  ! if desert, set deposition velocity to zero by multiplying
		  ! dep vel by the fraction covered by desert(hup, 5/1/2006)
		  ! IJLAND+1 is the Olson Land type index 
            	  ! 51: desert	52: desert set SVEL = 0
		  ! IJUSE is the fraction of the grid square occupied by surface K
		  ! in units of per mil (IJUSE=500 -> 50% of the grid square). 
		  DO K = 1, IREG(IREF,JREF)
	          
                     NTYP(I,J) = IJLAND(IJLOOP, K) + 1
                     
                     IF (NTYP(I,J) .eq. 52 .or. NTYP(I,J) .eq. 51) THEN 
		  	SVEL = SVEL*(1-(IJUSE(IJLOOP,K)/1.d3)) 
                     ENDIF
                     
		  ENDDO
		  
		  ! For HD add soil fractionation with 
                  ! an alpha coefficient of 0.943 Gerst & Quay, 2001
		  IF (N .eq. 2) SVEL = SVEL * 0.943


                  ! Get THIK  (cannot use ZH variable, since
                  ! DO_DRYDEP, METERO, and DEPVEL are not called in H2/HD sims)

                  ! Fraction of box (I,J,L) under PBL top [unitless]
                  F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )
      
                  ! Mixed layer depth [m]
                  MLD = GET_PBL_TOP_m( I, J )
      
                  ! THIK = thickness of surface layer [m]
                  THIK    = BXHEIGHT(I,J,1)
                  ! Add option for non-local PBL mixing scheme: 
                  ! only the surface (Lin, 03/31/09) 
                  IF (.NOT. LNLPBL) THIK    = MAX( MLD, THIK ) 
                  
                  ! Dry deposition frequency [1/s]
                  DRYF    = ( SVEL / 100.d0 ) / THIK
		  
		  ! FRACLOST = Fraction of species lost to drydep [unitless]
                  FRACLOST = DRYF * DTCHEM * F_UNDER_TOP
                  ! Add option for non-local PBL mixing scheme:
                  ! only the surface (Lin, 03/31/29)
                  IF (LNLPBL) FRACLOST = DRYF * DTCHEM

                  !========================================================
                  ! Proceed as follows:
                  ! -------
                  ! (a) If FRACLOST < 0, then stop the run.
                  !
                  ! (b) If FRACLOST > 1, use an exponential loss to 
                  !     avoid negative tracer
                  !
                  ! (c) If FRACLOST is in the range (0-1), then use the 
                  !     regular formula (STT * FRACLOST) to compute
                  !     the loss from dry deposition.
                  !========================================================

                  ! Stop the run on negative FRACLOST!
                  IF ( FRACLOST < 0 ) THEN
                     CALL ERROR_STOP( 'FRACLOST < 0', 'dryflxH2HD' )
                  ENDIF


                  ! AMT_LOST = amount of tracer lost to drydep [kg]
                  IF ( FRACLOST > 1 ) THEN
                     AMT_LOST = STT(I,J,L,NN) * ( 1d0 - EXP(-FRACLOST) )
     &                    * FSOIL
                  ELSE
                     AMT_LOST = STT(I,J,L,NN) * FRACLOST * FSOIL
                  ENDIF


                  ! ND44 diagnostic: drydep flux [kg/s]
                  ! ND44 diagnostic: drydep velocity [cm/s]
                  IF ( ND44 > 0 ) THEN
                     AD44(I,J,N,1) = AD44(I,J,N,1) + ( AMT_LOST/DTCHEM ) 
                     AD44(I,J,N,2) = AD44(I,J,N,2) + SVEL * FSOIL
                  ENDIF


                  ! Subtract AMT_LOST from the STT array [kg]
                  STT(I,J,L,NN)  = STT(I,J,L,NN) - AMT_LOST
		
                    
            ENDIF               ! I and J within bounds, ALBD<0.4 and FSOIL>0
         ENDDO                  ! M = LAND GRID BOXES
         ENDDO                  ! PBL layers
      ENDDO                     ! NUMDEP = Number of species that drydep


      ! Return to calling program
      END SUBROUTINE DRYFLXH2HD

!------------------------------------------------------------------------------

      SUBROUTINE DEPVEL( NPTS, RADIAT,  TEMP,  SUNCOS, F0,  HSTAR,
     &                   XMW,  AIROSOL, USTAR, CZ1,    OBK, CFRAC,
     &                   ZH,   LSNOW,   DVEL,  ZO,     RHB )

      ! References to F90 modules (bmy, 3/8/01)
      USE ERROR_MOD, ONLY : IT_IS_NAN
      USE TRACER_MOD, ONLY : ITS_A_POPS_SIM ! (clf, 1/3/2011)
      
                        
C     Subroutine computes the dry deposition velocities using 
C      a resistance-in-series model.
C
C** Contact: D.J. Jacob, Harvard U. (djj@io.harvard.edu)
C** Modularized by G.M. Gardner, Harvard U.
C** Version 3.2:   5/27/97
C** Version 3.2.1: 3/4/99   -- bug fix in expression for RT 
C** Version 3.2.2: 3/26/99  -- bug fix: specify a large Ra for aerosols
C** Version 3.2.3: 11/12/99 -- change Reynolds # criterion from 10 to 1
C                           -- force double precision w/ "D" exponents
C** Version 3.3:   5/8/00   -- bug fixes, cleanup, updated comments.
C** Version 3.4:   1/22/03  -- remove hardwire for CANOPYNOX
C** Version 3.5    7/21/03  -- Remove cap of surface resistance in RLUXX
C** Version 3.6    4/01/04  -- Now do drydep of DUST aerosol tracers
C** Version 3.7    4/20/04  -- Now also do drydep of SEASALT aerosol tracers
C** Version 3.8    4/13/05  -- Accounts for hygroscopic growth of SEASALT
C**                             aerosol tracers.  DUST aerosol tracers do
C**                             not grow hygroscopically.  Added RHB as
C**                             an input argument.
C** Version 3.9    5/25/05  -- Now restore GISS-specific code for GCAP model
C** Version 3.9.1  11/17/05 -- change Reynolds # criterion from 1 to 0.1
C
C***********************************************************************
C   Changes from Version 3.2 to Version 3.3:                         ***
C   * We now suppress dry deposition over aerodynamically smooth     ***
C     surfaces.  The previous algorithm yielded negative numbers     ***
C     when u* was very small (due to the logarithm going negative).  ***
C     See the comments below for more information.                   ***
C   * Now eliminate obsolete variables ZLMO and SIH from the code.   ***
C   * Obsolete comments have been updated or removed.                ***
C***********************************************************************
C   Changes from version 3.1 to version 3.2:                         ***
C   * In unstable atmospheres with |ZLMO| < ZO, as can happen        ***
C    occasionally under very low wind conditions with tall canopies, ***
C    application of Monin-Obukhov similarity yields negative values  ***
C    for RA.  This was a problem in version 3.1.  In fact,           ***
C    Monin-Obukhov similarity does not apply under such conditions,  ***
C    so we now set RA to zero and let the boundary                   ***
C    resistance RB define the overall aerodynamic resistance.  Since *** 
C    RB varies inversely with U* it will impose a large aerodynamic  ***
C    resistance under very low wind conditions.                      ***
C   * The range of applicability of stability correction functions   ***
C    to Monin-Obukhov similarity has been extended to                ***
C    -2.5 < z/zMO < 1.5, based on Figure 2 of Businger et al. [1971].***
C    The range used to be -1 < z/zMO < 1 in version 3.1.             ***
C***********************************************************************
C                                  
C  Literature cited:
C     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
C       resistance model for gaseous deposition to vegetated surfaces,
C       Atmos. Environ. 21, 91-101, 1987.
C     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
C     Businger, J.A., et al., Flux-profile relationships in the atmospheric 
C       surface layer, J. Atmos. Sci., 28, 181-189, 1971.
C     Dwight, H.B., Tables of integrals and other mathematical data,
C       MacMillan, 1957.
C     Guenther, A., and 15 others, A global model of natural volatile
C       organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
C     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
C       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
C     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
C       hydrocarbons, and ozone over the Amazon forest during the wet season,
C       J.  Geophys. Res., 95, 16737-16754, 1990.
C     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
C       J. Geophys. Res., 97, 16473-16479, 1992.
C     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York, 1988.
C     Munger, J.W., and 8 others, Atmospheric deposition of reactive
C       nitrogen oxides and ozone in a temperate deciduous forest and a
C       sub-arctic woodland, J. Geophys. Res., in press, 1996.
C     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2, sulfate, and
C       HNO3 deposition velocities computed using regional landuse and
C       meteorological data, Atmos. Environ., 20, 949-964, 1986.
C     Wang, Y.H., paper in preparation, 1996.
C     Wesely, M.L, Improved parameterizations for surface resistance to
C       gaseous dry deposition in regional-scale numerical models, 
C       Environmental Protection Agency Report EPA/600/3-88/025,
C       Research Triangle Park (NC), 1988.
C     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
C
C***********************************************************************       
C
C Need as landtype input for each grid square (I,J) see (RDLAND & CMN_VEL):
C     IJREG(JLOOP)       - # of landtypes in grid square
C     IJLAND(IJLOOP,LDT) - Land type ID for element LDT =1, IJREG(IJLOOP)
C                          (could be from any source - mapped to deposition 
C                          surface ID in input unit 65)
C     IJUSE(IJLOOP,LDT)  - Fraction ((per mil) of gridbox area occupied by
C                          land type element LDT
C
C Need as leaf area index see (RDLAI & CMN_VEL):
C     XYLAI(IJLOOP,LDT)  - Leaf Area Index of land type element LDT
C
C Need as meteorological input for each grid square(I,J) (passed):
C     RADIAT(IJLOOP) - Solar radiation in W m-2
C     TEMP(IJLOOP)   - Surface air temperature in K
C     SUNCOS(IJLOOP) - Cosine of solar zenith angle
C     LSNOW(IJLOOP)  - Logical for snow and sea ice
C     RHB(IJLOOP)    - Relative humidity at the surface
C
C Need as input for each species K (passed):
C     F0(K)          - reactivity factor for oxidation of biological substances
C     HSTAR(K)       - Henry's Law constant
C     XMW(K)         - Molecular weight (kg/mole) of species K
C                      (used to calculate molecular diffusivities)
C     AIROSOL(K)     - LOGICAL flag (T = aerosol species; 
C                                    F = gas-phase species)
C
C Also need to call the following subroutines to read drydep input data:
C     "modin.f"      - reads Olson land types, dry deposition land types,
C                      and roughness heights from "drydep.table". 
C                      (NOTE: For GEOS model, roughness heights are taken
C                       from met field input instead of from "drydep.table").
C     "rddrycf.f     - reads drydep polynomial coeff's from file "drydep.coef"
C     "rdlai.f"      - reads Leaf Area Indices from files "lai**.global"
C     "rdland.f"     - reads Olson land types from file "vegtype.global"
C
C Some variables used in the subroutine (passed):
C     LRGERA(IJLOOP) T -> stable atmosphere; a high aerodynamic resistance
C                        (RA=1.E4 m s-1) is imposed; else RA is calculated
C     USTAR(IJLOOP)  - Friction velocity (m s-1)
C     CZ1(IJLOOP)    - Altitude (m) at which deposition velocity is computed
C     OBK(IJLOOP)    - Monin-Obukhov length (m): set to 1.E5 m under neutral 
C                      conditions
C     CFRAC(IJLOOP)  - Fractional cloud cover
C     ZH(IJLOOP)     - Mixing depth (m)
C
C Some variables used in the subroutine:
C     MAXDEP         - the maximum number of species for which the dry 
C                      deposition calculation is done
C     ZO(LDT)        - Roughness height (m) for specific surface type indexed 
C                      by LDT
C     RSURFC(K,LDT)  - Bulk surface resistance (s m-1) for species K to 
C                      surface LDT
C     C1X(K)         - Total resistance to deposition (s m-1) for species K
C
C Returned:
C     DVEL(IJLOOP,K) - Deposition velocity (m s-1) of species K
C***********************************************************************       

#     include "CMN_SIZE"
#     include "CMN_VEL"
#     include "commsoil.h"

      INTEGER NPTS
      REAL*8  RADIAT(MAXIJ),TEMP(MAXIJ),SUNCOS(MAXIJ)
      REAL*8  USTAR(MAXIJ),CZ1(MAXIJ)
      REAL*8  OBK(MAXIJ),CFRAC(MAXIJ),ZH(MAXIJ)
      REAL*8  DVEL(MAXIJ,MAXDEP)

      ! Added relative humidity array (bec, bmy, 4/13/05)
      REAL*8 :: RHB(MAXIJ)

      REAL*8  RI(NTYPE),RLU(NTYPE),RAC(NTYPE),RGSS(NTYPE),
     1        RGSO(NTYPE),RCLS(NTYPE),RCLO(NTYPE),
     2        RSURFC(MAXDEP,NTYPE)       
                                   
      REAL*8  C1X(MAXDEP),VD(MAXDEP),VK(MAXDEP)              

#if   defined( GCAP )
      ! For the GISS/GCAP model, ZO is a function of land type 
      ! and is of dimension NTYPE (swu, bmy, 5/25/05)
      REAL*8  ZO(NTYPE)           
#else
      ! For GEOS-CTM, ZO is now of size MAXIJ and is passed via 
      ! the argument list, since it is a DAO met field. (bmy, 11/10/99)
      REAL*8  ZO(MAXIJ)           
#endif

      LOGICAL LDEP(MAXDEP)
      LOGICAL AIROSOL(MAXDEP)
      REAL*8  F0(MAXDEP),HSTAR(MAXDEP),XMW(MAXDEP)

      LOGICAL LRGERA(MAXIJ)

      REAL*8  VDS
      REAL*8  CZ,C1,RT,XNU,RAD0,RIX,GFACT,GFACI
      REAL*8  RDC,RLUXX,RGSX,RCL,DTMP1,DTMP2,DTMP3,DTMP4
      REAL*8  CZH,CKUSTR,REYNO,CORR1,CORR2,Z0OBK
      REAL*8  RA,RB,DUMMY1,DUMMY2,DUMMY3,DUMMY4
      REAL*8  XMWH2O,DAIR,TEMPK,TEMPC
      INTEGER IOLSON,II,IW
      INTEGER K,IJLOOP,LDT
      REAL*8  RCLX,RIXX,BIOFIT
      REAL*8  PRESS
      DATA PRESS /1.5D5/
C
C Logical for snow and sea ice
C

      LOGICAL LSNOW(MAXIJ)

      ! Adding logical switch for POPs in order to use octanol-air partitioning
      ! instead of Henry's law (water-air) for scaling of cuticular resistances
      ! only (since POPs accumulate in waxy part of cuticules), clf 1/3/2011

      LOGICAL :: IS_POPS

      


C***********************************************************************       

      ! Is this a POPs simulation?
      IS_POPS = ITS_A_POPS_SIM()  ! clf, 1/3/2011

C***********************************************************************
C
C
C** If LDEP(K)=F, species does not deposit.
C** Deposition is applied only to species with LDEP=T.

      DO K = 1,NUMDEP
         LDEP(K) = (HSTAR(K).GT.0.D0 .OR. F0(K).GT.0.D0 
     &              .OR. AIROSOL(K)) 
      ENDDO
                                   
      DO K = 1,NUMDEP    
         DO IJLOOP =1,NPTS
            DVEL(IJLOOP,K) = 0.0D0
         ENDDO
      ENDDO
C***********************************************************************
C*                                 
C*    Begin section for computing deposition velocities           
C*                                 
C*                                 
      ! Add parallel DO-loop (bmy, 2/22/05)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( IJLOOP, CZ,     TEMPK,  TEMPC,  K,      VD     )
!$OMP+PRIVATE( LDT,    RSURFC, C1,     XNU,    RT,     IOLSON )
!$OMP+PRIVATE( II,     RI,     RLU,    RAC,    RGSS,   RGSO   )
!$OMP+PRIVATE( RCLS,   RCLO,   RAD0,   RIX,    GFACT,  GFACI  )
!$OMP+PRIVATE( RDC,    XMWH2O, RIXX,   RLUXX,  RGSX,   RCLX   )
!$OMP+PRIVATE( DTMP1,  DTMP2,  DTMP3,  DTMP4,  VDS,    CZH    )
!$OMP+PRIVATE( CKUSTR, REYNO,  CORR1,  CORR2,  Z0OBK,  RA     )
!$OMP+PRIVATE( DUMMY1, DUMMY2, DUMMY3, DUMMY4, DAIR,   RB     )
!$OMP+PRIVATE( C1X,    VK                                     )
      DO 560 IJLOOP =1,NPTS

C** CZ is Altitude (m) at which deposition velocity is computed
         CZ = CZ1(IJLOOP)
C** TEMPK and TEMPC are surface air temperatures in K and in C
         TEMPK = TEMP(IJLOOP)
         TEMPC = TEMP(IJLOOP)-273.15D0
C* Initialize variables
         DO K = 1,NUMDEP
            VD(K) = 0.0D0
            DO LDT = 1,NTYPE
               RSURFC(K,LDT) = 0.D0
            END DO
         END DO

C** Calculate the kinematic viscosity XNU (m2 s-1) of air
C** as a function of temperature.
C** The kinematic viscosity is used to calculate the roughness heights over
C** water surfaces and to diagnose whether such surfaces are aerodynamically
C** rough or smooth using a Reynolds number criterion.
C** The expression for the temperature dependence of XNU
C** is from the FORTRAN code in Appendix II of Wesely [1988];
C** I wasn't able to find an original reference but it seems benign enough.
         C1 = TEMPK/273.15D0
         XNU = 0.151D0*(C1**1.77D0)*1.0D-04

C* Compute bulk surface resistance for gases.    
C*                                 
C* Adjust external surface resistances for temperature; 
C* from Wesely [1989], expression given in text on p. 1296.        
C*                                 
C* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
C* so the inner parentheses are not needed (bmy, 3/4/99)
C*        RT = 1000.0*EXP(-(TEMPC-4.0))
         RT = 1000.0D0*EXP(-TEMPC-4.0D0)
C*
C    Get surface resistances - loop over land types LDT
C***************************************************************************
C* The land types within each grid square are defined using the Olson
C* land-type database.  Each of the Olson land types is assigned a 
C* corresponding "deposition land type" with characteristic values of surface
C* resistance components.  There are 74 Olson land-types but only 11 deposition
C* land-types (i.e., many of the Olson land types share the same deposition
C* characteristics).  Surface resistance components for the "deposition land 
C* types" are from Wesely [1989] except for tropical forests [Jacob and Wofsy,
C* 1990] and for tundra [Jacob et al., 1992].  All surface resistance 
C* components are normalized to a leaf area index of unity.
C*
C* Olson land types, deposition land types, and surface resistance components
C* are read from file 'drydep.table'; check that file for further details.
C****************************************************************************
         DO 170 LDT = 1,IJREG(IJLOOP)
            IF (IJUSE(IJLOOP,LDT) .EQ. 0) GOTO 170
            IOLSON = IJLAND(IJLOOP,LDT)+1
            II = IDEP(IOLSON)
C
C** If the surface to be snow or ice;
C** set II to 1 instead.
C
            IF(LSNOW(IJLOOP)) II=1

C* Read the internal resistance RI (minimum stomatal resistance for water 
C* vapor,per unit area of leaf) from the IRI array; a '9999' value means no 
C* deposition to stomata so we impose a very large value for RI.

            RI(LDT) = DBLE(IRI(II))
            IF (RI(LDT)   .GE. 9999.D0) RI(LDT)   = 1.D12

C** Cuticular resistances IRLU read in from 'drydep.table'
C** are per unit area of leaf;
C** divide them by the leaf area index to get a cuticular resistance for the
C** bulk canopy.  If IRLU is '9999' it means there are no cuticular
C** surfaces on which to deposit so we impose a very large value for RLU.
            IF (IRLU(II) .GE. 9999 .OR. 
     &           XYLAI(IJLOOP,LDT).LE.0.D0) THEN
               RLU(LDT)  = 1.D6
            ELSE
               RLU(LDT)= DBLE(IRLU(II))/XYLAI(IJLOOP,LDT)+RT
            ENDIF
C** The following are the remaining resistances for the Wesely
C** resistance-in-series model for a surface canopy
C** (see Atmos. Environ. paper, Fig.1).  
            RAC(LDT)  = MAX(DBLE(IRAC(II)), 1.D0)
            IF (RAC(LDT)  .GE. 9999.D0) RAC(LDT)  = 1.D12
            RGSS(LDT) = MAX(DBLE(IRGSS(II)) + RT ,1.D0)
            IF (RGSS(LDT) .GE. 9999.D0) RGSS(LDT) = 1.D12
            RGSO(LDT) = MAX(DBLE(IRGSO(II)) + RT ,1.D0) 
            IF (RGSO(LDT) .GE. 9999.D0) RGSO(LDT) = 1.D12
            RCLS(LDT) = DBLE(IRCLS(II)) + RT           
            IF (RCLS(LDT) .GE. 9999.D0) RCLS(LDT) = 1.D12
            RCLO(LDT) = DBLE(IRCLO(II)) + RT          
            IF (RCLO(LDT) .GE. 9999.D0) RCLO(LDT) = 1.D12
C***************************************************************************
C*                                 
C*    Adjust stomatal resistances for insolation and temperature:  
C*
C*     Temperature adjustment is from Wesely [1989], equation (3).
C*       
C*     Light adjustment by the function BIOFIT is described by Wang [1996].
C*     It combines
C*       - Local dependence of stomal resistance on the intensity I of light 
C*         impinging the leaf; this is expressed as a mutliplicative 
C*         factor I/(I+b) to the stomatal resistance where b = 50 W m-2
C*         (equation (7) of Baldocchi et al. [1987])
C*       - radiative transfer of direct and diffuse radiation in the 
C*         canopy using equations (12)-(16) from Guenther et al. [1995]
C*       - separate accounting of sunlit and shaded leaves using
C*         equation (12) of Guenther et al. [1995]
C*       - partitioning of the radiation at the top of the canopy into direct
C*         and diffuse components using a parameterization to results from
C*         an atmospheric radiative transfer model [Wang, 1996]
C*     The dependent variables of the function BIOFIT are the leaf area
C*     index (XYLAI), the cosine of zenith angle (SUNCOS) and the fractional
C*     cloud cover (CFRAC).  The factor GFACI integrates the light
C*     dependence over the canopy depth; sp even though RI is input per
C*     unit area of leaf it need not be scaled by LAI to yield a bulk
C*     canopy value because that's already done in the GFACI formulation.
C***************************************************************************

            RAD0 = RADIAT(IJLOOP)
            RIX = RI(LDT)
            IF (RIX .GE. 9999.D0) GO TO 150
            GFACT = 100.0D0
            IF (TEMPC .GT. 0.D0 .AND. TEMPC .LT. 40.D0)
     *           GFACT = 400.D0/TEMPC/(40.0D0-TEMPC)
            GFACI = 100.D0
            IF (RAD0.GT.0.D0 .AND. XYLAI(IJLOOP,LDT).GT.0.D0) THEN
               GFACI=1.D0/BIOFIT(DRYCOEFF,XYLAI(IJLOOP,LDT),
     *              SUNCOS(IJLOOP),CFRAC(IJLOOP))
            ENDIF
            
            RIX = RIX*GFACT*GFACI
 150        CONTINUE
C*                                 
C*    Compute aerodynamic resistance to lower elements in lower part           
C*    of the canopy or structure, assuming level terrain - 
C*    equation (5) of Wesely [1989].
C*                                 
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))
C*
C*    Loop over species; species-dependent corrections to resistances
C*    are from equations (6)-(9) of Wesely [1989].
C*
            DO 160  K = 1,NUMDEP
C** exit for non-depositing species or aerosols.
               IF (.NOT. LDEP(K) .OR. AIROSOL(K)) GOTO 155
               XMWH2O = 18.D-3
               RIXX = RIX*DIFFG(TEMPK,PRESS,XMWH2O)/
     C              DIFFG(TEMPK,PRESS,XMW(K))
     C              + 1.D0/(HSTAR(K)/3000.D0+100.D0*F0(K))
               RLUXX = 1.D12
               IF (RLU(LDT).LT.9999.D0)
     C              RLUXX = RLU(LDT)/(HSTAR(K)/1.0D+05 + F0(K))

               ! If a POPs simulation, scale cuticular resistances with 
               ! octanol-air partition coefficient (Koa) instead of HSTAR
               ! (clf, 1/3/2011)

               IF (IS_POPS) THEN
                    RLUXX = RLU(LDT)/(KOA(K)/1.0D+05 + F0(K))
               ENDIF
C*
C* To prevent virtually zero resistance to species with huge HSTAR, such
C* as HNO3, a minimum value of RLUXX needs to be set. The rationality
C* of the existence of such a minimum is demonstrated by the observed
C* relationship between Vd(NOy-NOx) and Ustar in Munger et al.[1996];
C* Vd(HNO3) never exceeds 2 cm s-1 in observations. The
C* corresponding minimum resistance is 50 s m-1.  This correction
C* was introduced by J.Y. Liang on 7/9/95.
C*
               !-----------------------------------------------------------
               ! Prior to 7/21/03:
               ! Remove the cap of surface resistance (rjp, bmy, 7/21/03)
               !IF(RLUXX.LT. 50.D0) RLUXX= 50.D0
               !-----------------------------------------------------------
C     
               RGSX = 1.D0/(HSTAR(K)/1.0D+05/RGSS(LDT) + 
     1              F0(K)/RGSO(LDT))
               RCLX = 1.D0/(HSTAR(K)/1.0D+05/RCLS(LDT) + 
     1              F0(K)/RCLO(LDT))
C*
C** Get the bulk surface resistance of the canopy, RSURFC, from the network
C** of resistances in parallel and in series (Fig. 1 of Wesely [1989])
               DTMP1=1.D0/RIXX
               DTMP2=1.D0/RLUXX
               DTMP3=1.D0/(RAC(LDT)+RGSX)
               DTMP4=1.D0/(RDC+RCLX)
               RSURFC(K,LDT) = 1.D0/(DTMP1 + DTMP2 + DTMP3 + DTMP4)
C  Save the within canopy depvel of NOx, used in calculating the 
C  canopy reduction factor for soil emissions.
               ! Remove hardwire for CANOPYNOX (bmy, 1/24/03)
               IF ( K == DRYDNO2 ) THEN
                  CANOPYNOX(IJLOOP,LDT)=DTMP1+DTMP2+DTMP3+DTMP4
               ENDIF
C** get surface deposition velocity for aerosols if needed;
C** equations (15)-(17) of Walcek et al. [1986]
 155           IF (.NOT. AIROSOL(K)) GOTO 160

              !===========================================================
              ! The difference between sea-salt and dust tracers below
              ! is whether or not we account for hygroscopic growth.
	      ! Seasalt (yes), Dust (no)  (bec, bmy, 4/13/05 )
              !===========================================================

               IF ( ( DEPNAME(K) == 'SALA' )  .OR. 
     &              ( DEPNAME(K) == 'SALC' )  .OR. 
     &              ( DEPNAME(K) == 'SO4S' )  .OR. 
     &              ( DEPNAME(K) == 'NITS' ) ) THEN 

                  !=====================================================
                  ! Use size-resolved dry deposition calculations for 
                  ! seasalt aerosols.  We need to account for the
                  ! hygroscopic growth of the aerosol particles.
                  ! (rjp, bec, bmy, 4/13/05)
                  !=====================================================

!---------------------------------------------------------------------------
! NOTE: We need to add a new subroutine if you want to use the
!       Seinfeld 1986 mechanism (bec, bmy, 4/13/05)
!                  ! [Seinfeld, 1986] 
!                  RSURFC(K,LDT) = 
!     &             AERO_sfcRsI(K, II, PRESS*1D-3, TEMPK, USTAR(IJLOOP))
!---------------------------------------------------------------------------

                  ! [Zhang et al., 2001]
                  RSURFC(K,LDT) = 
     &               AERO_SFCRSII( K,     II,            PRESS*1D-3, 
     &                             TEMPK, USTAR(IJLOOP), RHB(IJLOOP) )

               ELSE IF ( ( DEPNAME(K) == 'DST1' )  .OR. 
     &                   ( DEPNAME(K) == 'DST2' )  .OR. 
     &                   ( DEPNAME(K) == 'DST3' )  .OR. 
     &                   ( DEPNAME(K) == 'DST4' ) ) THEN 

                  !=====================================================
                  ! Use size-resolved dry deposition calculations for 
                  ! dust aerosols only.  Do not account for hygroscopic
                  ! growth of the dust aerosol particles.
                  ! (rjp, bec, bmy, 4/13/05)
                  !=====================================================     

!                  ! [Seinfeld, 1986] 
!                  RSURFC(K,LDT) = 
!     &             DUST_sfcRsI(K, II, PRESS*1D-3, TEMPK, USTAR(IJLOOP))

                  ! [Zhang et al., 2001]
                  RSURFC(K,LDT) = 
     &             DUST_SFCRSII(K, II, PRESS*1D-3, TEMPK, USTAR(IJLOOP))

               ELSE 

                  !=====================================================
                  ! Replace original code to statement 160 here: only
                  ! do this for non-size-resolved tracers where 
                  ! AIROSOL(K)=T. (rjp, tdf, bec, bmy, 4/20/04)
                  !=====================================================
!                  ! use Zhang et al for all aerosols (hotp 10/26/07)
!                  VDS = 0.002D0*USTAR(IJLOOP)
!                  IF (OBK(IJLOOP) .LT. 0.0D0) THEN
!                     VDS = VDS*(1.D0+(-300.D0/OBK(IJLOOP))**0.6667D0)
!                  ENDIF
!C***                               
!                  IF ( OBK(IJLOOP) .EQ. 0.0D0 )
!     c                 WRITE(6,156) OBK(IJLOOP),IJLOOP,LDT
! 156              FORMAT(1X,'OBK(IJLOOP)=',E11.2,1X,' IJLOOP =',I4,
!     c                   1X,'LDT=',I3/) 
!                  CZH  = ZH(IJLOOP)/OBK(IJLOOP)
!                  IF (CZH.LT.-30.0D0) VDS = 0.0009D0*USTAR(IJLOOP)*
!     x                                 (-CZH)**0.6667D0

                  RSURFC(K, LDT) =
     &                ADUST_SFCRSII(K, II, PRESS*1D-3, TEMPK, 
     &                              USTAR(IJLOOP))
                  
C*                                 
C*    Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
C*    VDSMAX is taken from Table 2 of Walcek et al. [1986].
C*    Invert to get corresponding R

!                  RSURFC(K,LDT) = 1.D0/MIN(VDS, DBLE(IVSMAX(II))/1.D4)
               ENDIF
 160        CONTINUE
C*
 170     CONTINUE
C*
C*    Set max and min values for bulk surface resistances         
C*                                 
         DO 190 K = 1,NUMDEP 
            IF (.NOT.LDEP(K)) GOTO 190
            DO 180 LDT = 1,IJREG(IJLOOP)
               IF (IJUSE(IJLOOP,LDT) .EQ. 0) GOTO 180
               RSURFC(K,LDT)= MAX(1.D0, MIN(RSURFC(K,LDT), 9999.D0))
 180        CONTINUE
 190     CONTINUE
C*                                 
C*    Loop through the different landuse types present in the grid square     
C*                                 
         DO 500 LDT=1, IJREG(IJLOOP)
            IF (IJUSE(IJLOOP,LDT) .EQ. 0) GOTO 500
            IOLSON = IJLAND(IJLOOP,LDT)+1

#if   defined( GCAP )
! NOTE: This section only applies to the GCAP/GISS model (swu, bmy, 5/25/05)
!** Get roughness heights; they are specified constants for each surface
!** type except over water where zo = f(u*).  The latter dependence
!** is from equation (6) of Hicks and Liss [1976]. 
            DO 200 IW=1,NWATER
               IF (IOLSON .NE. IWATER(IW)) GOTO 200
               ZO(LDT) = 1.4D-02*USTAR(IJLOOP)*USTAR(IJLOOP)/9.8D0   
     1                 + 1.1D-01*XNU/USTAR(IJLOOP)        
               GOTO 210
 200        CONTINUE
            ZO(LDT) = DBLE(IZO(IOLSON))*1.D-4
 210        CONTINUE
#endif

C***** Get aerodynamic resistances Ra and Rb. ***********************
C   The aerodynamic resistance Ra is integrated from altitude z0+d up to the
C   altitude z1 at which the dry deposition velocity is to be referenced.
C   The integration corrects for stability using Monin-Obukhov similarity 
C   formulas from Businger et al. [1971] which apply over the range 
C   -2.5 < z/zMO < 1.5 (see their Figure 2).
C   Under very unstable conditions when z1 > -2.5 zMO, we assume that there is
C   no resistance to transfer in the convective column between zMO and z1.
C   Under very stable conditions when z1 > 1.5 zMO, we assume that vertical
C   transfer in the column between zMO and z1 is strongly suppressed so
C   that the deposition velocity at altitude z1 is very low.  Under these
C   conditions we just specify a very large Ra=1.E4 s m-1 (LRGERA = T).
C**
C   The Reynolds number REYNO diagnoses whether a surface is
C   aerodynamically rough (REYNO > 1) or smooth.  
C
C   NOTE: The criterion "REYNO > 1" was originally "REYNO > 10".
C         See below for an explanation of why it was changed (hyl, 10/15/99)
C
C   Surface is rough in all cases except over water with low wind speeds.  
C   In the smooth case, vertical transport IN THE SUBLAYER near the surface 
C   is limited by molecular diffusion and is therefore very slow; we assign 
C   a large value we assign a large value of Ra + Rb to account for this 
C   effect.  [In Versions 3.2 and earlier we used the formulation for Ra + Rb 
C   given in Equation (12) of Walcek et al [1986] to calculate the aerodynamic 
C   resistance over smooth surfaces.  However, that expression fails when 
C   u* is very small, as it yields negative values of Ra + Rb].  
C   (djj, hyl, bmy, 5/8/00)
C**
C   In the aerodynamically rough case, the expression for Ra is as
C   given in equation (5) of Jacob et al. [1992]:
C
C          Ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
C
C    where x = (z-D)/zMO, z is the height above ground, and D is the
C    displacement height which is typically 70-80% of the canopy height
C    [Brutsaert, 1982].  We change the vertical coordinate so that z=0 at
C    the displacement height; that's OK since for all practical applications
C    z1 >> D.  In this manner we don't need to assume any specific value for
C    the displacement height.  Applying the variable transformation 
C    z -> x = z/zMO, the equation above becomes
C
C          Ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx   with x=z/zMO
C
C    Here phi is a stability correction function originally formulated by
C    Businger et al. [1971] and given in eqns 5a and 5b of Jacob et al. [1992].
C    For unstable conditions,
C
C          phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
C
C    The analytical solution to the integral is 
C    [Dwight, 1957, integral 192.11]:
C
C          int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)/(sqrt(1-bx)+1)))
C
C    which yields the expression for Ra used in the code for unstable 
C    conditions.  For stable conditions,
C
C          phi(x) = a + bx        where a=0.74, b = 4.7
C
C    and the analytical solution to the integral is
C
C          int((a/x)+b)dx = a*ln(x) + bx
C
C    which yields the expression of Ra used in the code for stable conditions.
C**
C   The formulation of RB for gases is equation (12) of 
C   Walcek et al. [1986].  The parameterization for deposition of
C   aerosols does not include an RB term so RB for aerosols is set
C   to zero.
C   Modify phi(x) according to the non-local mixing scheme by Holtslag 
C   and Boville [1993] ( Lin, 07/18/08 )
C   For unstable conditions,
C          phi(x) = a/sqrt(1-bx)  where a=1.0, b=15.0
C
C   For stable conditions,
C          phi(x) = a + bx
C              where a=1.0, b=5.0 for 0 <= x <= 1, and
C                    a=5.0, b=1.0 for x > 1.0 
C*********************************************************************
            CKUSTR = XCKMAN*USTAR(IJLOOP)

            ! Define REYNO for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
            REYNO = USTAR(IJLOOP)*ZO(LDT)/XNU 
#else
            REYNO = USTAR(IJLOOP)*ZO(IJLOOP)/XNU  
#endif

            IF ( OBK(IJLOOP) .EQ. 0.0D0 )
     c           WRITE(6,211) OBK(IJLOOP),IJLOOP,LDT                
 211        FORMAT(1X,'OBK(IJLOOP)=',E11.2,1X,' IJLOOP = ',I4,1X,
     c           'LDT=',I3/) 
            CORR1 = CZ/OBK(IJLOOP)

            ! Define Z0OBK for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
            Z0OBK = ZO(LDT)/OBK(IJLOOP)
#else
            Z0OBK = ZO(IJLOOP)/OBK(IJLOOP) 
#endif

            LRGERA(IJLOOP) = .FALSE.
            ! Add option for non-local PBL (Lin, 03/31/09) 
            IF (.NOT. LNLPBL) THEN
               IF (CORR1 .GT. 0.D0) THEN
                  IF (CORR1 .GT.  1.5D0) LRGERA(IJLOOP) = .TRUE.
               ELSEIF(CORR1 .LE. 0.D0) THEN
                  IF (CORR1 .LE. -2.5D0) CORR1 = -2.5D0
                  CORR2 = LOG(-CORR1)
               ENDIF
            ENDIF
C*                                 
            IF (CKUSTR.EQ.0.0D0) THEN
               WRITE(6,212) IJLOOP,CKUSTR,XCKMAN,USTAR(IJLOOP)
 212           FORMAT(1X,'IJLOOP= ',I4,1X,'CKUSTR=',E10.1,1X,
     x              'XCKMAN= ',E12.4,1X,'USTAR(IJLOOP)= ',
     x              E12.4)
               CLOSE(98)
               STOP             ! debug
            ENDIF
C
C
C...aerodynamically rough or smooth surface
C "In the classic study by Nikuradse (1933) the transition from smooth
C  to rough was examined in pipe flow. He introduced a roughness Reynolds
C  number Rr = U* Z0 / Nu and found the flow to be smooth for Rr < 0.13
C  and rough for Rr > 2.5 with a transition regime in between."
C  (E.B. Kraus and J.A. Businger, Atmosphere-Ocean Interaction, second
C  edition, P.144-145, 1994). Similar statements can be found in the books:
C  Evaporation into the atmosphere, by Wilfried Brutsaert, P.59,89, 1982;
C  or Seinfeld & Pandis, P.858, 1998. Here we assume a sudden transition
C  point Rr = 1 from smooth to rough, following L. Merlivat (1978, The
C  dependence of bulk evaporation coefficients on air-water interfacial
C  conditions as determined by the isotopic method, J. Geophys. Res.,
C  Oceans & Atmos., 83, C6, 2977-2980). Also refer to Brutsaert's book,
C  P.125. We used to use the criterion "REYNO > 10" for aerodynamically
C  rough surface and now change to "REYNO > 1". (hyl, 10/15/99)
C  
C  11/17/05: D. J. Jacob says to change the criterion for aerodynamically
C  rough surface to REYNO > 0.1 (eck, djj, bmy, 11/17/05)
            IF ( REYNO < 0.1d0 ) GOTO 220

            ! Add option for non-local PBL (Lin, 03/31/09) 
            IF (.NOT. LNLPBL) THEN

C...aerodynamically rough surface.
C*                                 
               IF (CORR1.LE.0.0D0 .AND. Z0OBK .LT. -1.D0)THEN
C*... unstable condition; set RA to zero. (first implemented in V. 3.2)
                  RA = 0.D0
               ELSEIF (CORR1.LE.0.0D0 .AND. Z0OBK .GE. -1.D0) THEN
C*... unstable conditions; compute Ra as described above.
                  DUMMY1 = (1.D0 - 9D0*CORR1)**0.5D0
                  DUMMY2 = (1.D0 - 9D0*Z0OBK)**0.5D0
                  DUMMY3 = ABS((DUMMY1 - 1.D0)/(DUMMY1 + 1.D0))
                  DUMMY4 = ABS((DUMMY2 - 1.D0)/(DUMMY2 + 1.D0))
                  RA = 0.74D0* (1.D0/CKUSTR) * LOG(DUMMY3/DUMMY4)
                                 
               ELSEIF((CORR1.GT.0.0D0).AND.(.NOT.LRGERA(IJLOOP)))THEN
C*...moderately stable conditions (z/zMO <1); compute Ra as described above
                  RA = (1D0/CKUSTR) * 
     &                 (.74D0*LOG(CORR1/Z0OBK) + 4.7D0*(CORR1-Z0OBK))
               ELSEIF(LRGERA(IJLOOP)) THEN
C*... very stable conditions
                  RA = 1.D+04
               ENDIF
C* check that RA is positive; if RA is negative (as occasionally
C* happened in version 3.1) send a warning message.

            ELSE

               IF (CORR1.LT.0.0D0) THEN
C*... unstable conditions; compute Ra as described above.
                  !coef_a=1.d0
                  !coef_b=15.d0
                  DUMMY1 = (1.D0 - 15.D0*CORR1)**0.5D0
                  DUMMY2 = (1.D0 - 15.D0*Z0OBK)**0.5D0
                  DUMMY3 = ABS((DUMMY1 - 1.D0)/(DUMMY1 + 1.D0))
                  DUMMY4 = ABS((DUMMY2 - 1.D0)/(DUMMY2 + 1.D0))
                  RA = 1.D0 * (1.D0/CKUSTR) * LOG(DUMMY3/DUMMY4)
               ELSEIF((CORR1.GE.0.0D0).AND.(CORR1.LE.1.0D0)) THEN
                  !coef_a=1.d0
                  !coef_b=5.d0
                  RA = (1D0/CKUSTR) *
     &                 (1.D0*LOG(CORR1/Z0OBK) + 5.D0*(CORR1-Z0OBK))
               ELSE ! CORR1 .GT. 1.0D0
                  !coef_a=5d0
                  !coef_b=1.d0
                  RA = (1D0/CKUSTR) *
     &                 (5.D0*LOG(CORR1/Z0OBK) + 1.D0*(CORR1-Z0OBK))
               ENDIF

            ENDIF

            RA = MIN(RA,1.D4)

#if   defined( GCAP )
            ! Debug output for GISS/GCAP model (swu, bmy, 5/25/05)
            IF (RA .LT. 0.) THEN
               WRITE (6,1001) IJLOOP,RA,CZ,ZO(LDT),OBK(IJLOOP)
               RA = 0.0D0
            ENDIF
#else
            ! For GEOS-CTM, We use ZO(MAXIJ), and IJLOOP is the index.
            ! Also, if RA is < 0, set RA = 0 (bmy, 11/12/99)
            IF (RA .LT. 0.D0) THEN
               WRITE (6,1001) IJLOOP,RA,CZ,ZO(IJLOOP),OBK(IJLOOP)  
               RA = 0.0D0
            ENDIF
#endif
 1001       FORMAT('WARNING: RA < 0 IN SUBROUTINE DEPVEL',
     &             I10,4(1X,E12.5))
C* Get total resistance for deposition - loop over species.
            DO 215 K = 1,NUMDEP 
               IF (.NOT.LDEP(K)) GOTO 215
C** DAIR is the thermal diffusivity of air; value of 0.2*1.E-4 m2 s-1 
C** cited on p. 16,476 of Jacob et al. [1992]
               DAIR = 0.2D0*1.D-4
               RB = (2.D0/CKUSTR)*
     x              (DAIR/DIFFG(TEMPK,PRESS,XMW(K)))**0.667D0
               IF (AIROSOL(K)) RB=0.D0
               C1X(K) = RA + RB + RSURFC(K,LDT)
 215        CONTINUE
            GOTO 240
 220        CONTINUE 
C** ... aerodynamically smooth surface
C** BUG FIX -- suppress drydep over smooth surfaces by setting Ra to a large
C** value (1e4).  This prevents negative dry deposition velocities when u*
C** is very small (djj, bmy, 5/8/00)
            DO 230 K = 1,NUMDEP 
               IF ( LDEP(K) ) THEN
                  RA     = 1.0D4
                  C1X(K) = RA + RSURFC(K,LDT)
               ENDIF
 230        CONTINUE

 240        CONTINUE
C*                                 
C* IJUSE is the fraction of the grid square occupied by surface LDT
C* in units of per mil (IJUSE=500 -> 50% of the grid square).  Add
C* the contribution of surface type LDT to the deposition velocity;
C* this is a loop over all surface types in the gridbox.
C*
            DO 400 K = 1,NUMDEP
               IF (.NOT.LDEP(K)) GOTO 400
               VK(K) = VD(K)
               VD(K) = VK(K) +.001D0*DBLE(IJUSE(IJLOOP,LDT))/C1X(K)
 400        CONTINUE
 500     CONTINUE

C** Load array DVEL
         DO 550 K=1,NUMDEP
            IF (.NOT.LDEP(K)) GOTO 550
            DVEL(IJLOOP,K) = VD(K)
            
            ! Now check for negative deposition velocity 
            ! before returning to calling program (bmy, 4/16/00)
            ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
            IF ( DVEL(IJLOOP,K) < 0d0 ) THEN
!$OMP CRITICAL
               PRINT*, 'DEPVEL: Deposition velocity is negative!'
               PRINT*, 'Dep. Vel = ', DVEL(IJLOOP,K)
               PRINT*, 'Species  = ', K
               PRINT*, 'IJLOOP   = ', IJLOOP
               PRINT*, 'RADIAT   = ', RADIAT(IJLOOP)
               PRINT*, 'TEMP     = ', TEMP(IJLOOP)
               PRINT*, 'SUNCOS   = ', SUNCOS(IJLOOP)
               PRINT*, 'USTAR    = ', USTAR(IJLOOP)
               PRINT*, 'CZ1      = ', CZ1(IJLOOP)
               PRINT*, 'OBK      = ', OBK(IJLOOP)
               PRINT*, 'CFRAC    = ', CFRAC(IJLOOP)
               PRINT*, 'ZH       = ', ZH(IJLOOP)
               PRINT*, 'LRGERA   = ', LRGERA(IJLOOP)
               PRINT*, 'ZO       = ', ZO(IJLOOP)
               PRINT*, 'STOP in depvel.f!'
               CALL CLEANUP
               STOP
!$OMP END CRITICAL
            ENDIF

            ! Now check for IEEE NaN (not-a-number) condition 
            ! before returning to calling program (bmy, 4/16/00)
            ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
            IF ( IT_IS_NAN( DVEL(IJLOOP,K) ) ) THEN
!$OMP CRITICAL
               PRINT*, 'DEPVEL: Deposition velocity is NaN!'
               PRINT*, 'Dep. Vel = ', DVEL(IJLOOP,K)
               PRINT*, 'Species  = ', K
               PRINT*, 'IJLOOP   = ', IJLOOP
               PRINT*, 'RADIAT   = ', RADIAT(IJLOOP)
               PRINT*, 'TEMP     = ', TEMP(IJLOOP)
               PRINT*, 'SUNCOS   = ', SUNCOS(IJLOOP)
               PRINT*, 'USTAR    = ', USTAR(IJLOOP)
               PRINT*, 'CZ1      = ', CZ1(IJLOOP)
               PRINT*, 'OBK      = ', OBK(IJLOOP)
               PRINT*, 'CFRAC    = ', CFRAC(IJLOOP)
               PRINT*, 'ZH       = ', ZH(IJLOOP)
               PRINT*, 'LRGERA   = ', LRGERA(IJLOOP)
               PRINT*, 'ZO       = ', ZO(IJLOOP)
               CALL CLEANUP
               STOP
!$OMP END CRITICAL
            ENDIF
 550     CONTINUE
 560  CONTINUE
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DEPVEL

!------------------------------------------------------------------------------

      FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
!
!******************************************************************************
!  Function DIFFG calculates the molecular diffusivity [m2/s] in air for a
!  gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].  (bmy, 5/16/06)
!
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TK    (REAL*8) : Temperature [K]
!  (2 ) PRESS (REAL*8) : Pressure [Pa]
!  (3 ) XM    (REAL*8) : Molecular weight of gas [kg]
!
!  NOTES:
!  (1 ) Originally was a standalone function; now bundled into drydep_mod.f.
!        Also now force REAL*8 precision with D exponents.  Now use F90
!        style syntax and updated comments. (bmy, 5/16/06)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TK
      REAL*8, INTENT(IN) :: PRESS
      REAL*8, INTENT(IN) :: XM

      ! Local variables
      REAL*8             :: AIRDEN, Z, DIAM, FRPATH, SPEED, DIFF_G    
      REAL*8, PARAMETER  :: XMAIR  = 28.8d-3 
      REAL*8, PARAMETER  :: RADAIR = 1.2d-10
      REAL*8, PARAMETER  :: PI     = 3.1415926535897932d0
      REAL*8, PARAMETER  :: RADX   = 1.5d-10
      REAL*8, PARAMETER  :: RGAS   = 8.32d0
      REAL*8, PARAMETER  :: AVOGAD = 6.023d23

      !=================================================================
      ! DIFFG begins here!
      !=================================================================

      ! Air density
      AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )

      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR

      ! Calculate the mean free path for gas X in air: 
      ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1d0 /( PI * SQRT( 1d0 + Z ) * AIRDEN*( DIAM**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8d0 * RGAS * TK / ( PI * XM ) )

      ! Calculate diffusion coefficient of gas X in air; 
      ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3d0 * PI / 32d0 ) * ( 1d0 + Z ) * FRPATH * SPEED

      ! Return to calling program
      END FUNCTION DIFFG

!------------------------------------------------------------------------------

      SUBROUTINE MODIN
!
!******************************************************************************
!  Subroutine MODIN reads Olson's data from the file "drydep.table".
!  (bmy, 4/1/02, 7/20/04)
!
!  NOTE: The roughness heights (IZO) from "drydep.table" are supplanted by
!  the Z0 field from the DAO met field archive.  The old GISS-II routines did
!  not archive Z0 as a met field, so roughness heights for each land type
!  were specified in this file.  This is historical baggage, but we still
!  need to keep IZO for compatibility w/ existing routine "depvel.f".  
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Wesely, M.L., 1988.
!  (2 ) Wesely, M.L., 1989.
!
!  NOTES:
!  (1 ) MODIN is one of the original GEOS-CHEM subroutines, that go back
!        to the days of the GISS-II code.  This has been cleaned up and
!        new comments added.  Also use subroutine "ioerror.f" to trap
!        I/O errors across all platforms . Now read the "drydep.table" file 
!        from the DATA_DIR/drydep_200203/ directory. (bmy, 4/1/02) 
!  (2 ) Remove obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE as the file unit
!        number instead of IUNIT. (bmy, 6/27/02)
!  (3 ) Now bundled into "drydep_mod.f".  Changed NVEGTYPE to NNVEGTYPE.
!        (bmy, 11/21/02)
!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER             :: L, IOLSON, I, IOS, IUNIT
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! MODIN begins here!
      !=================================================================

      ! Logical unit number
      IUNIT    = IU_FILE

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'drydep_200203/drydep.table'

      WRITE( 6, 50 ) TRIM( FILENAME )
 50   FORMAT( '     - MODIN: Reading ', a )

      ! Open file
      OPEN( IUNIT, FILE=TRIM( FILENAME ), FORM='FORMATTED', 
     &             STATUS='OLD',          IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:1' ) 

      ! Read 5 header comment lines
      DO L = 1, 5
         READ( IUNIT, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:2' )
      ENDDO

      !=================================================================
      ! For each of the NVEGTYPE Olson land types, read:
      !
      ! IOLSON (INTEGER) : Olson surface type ID #
      ! IDEP   (INTEGER) : Drydep ID # corresponding to IOLSON
      ! IZO    (INTEGER) : Roughness height [1e-4 m]
      !=================================================================
      DO L = 1, NNVEGTYPE
         READ( IUNIT, '(3i6)', IOSTAT=IOS )
     &      IOLSON, IDEP(IOLSON), IZO(IOLSON)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:3' )
      ENDDO

      ! Read comment line
      READ( IUNIT, '(a)', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:4' )
      
      !=================================================================
      ! For the water surface types, zO is input as 1.E-4 m but is 
      ! recalculated elsewhere as function of wind speed.  Read the # 
      ! of Olson's surface types that are water (NWATER) and the 
      ! corresponding ID's (IWATER)
      !=================================================================
      READ( IUNIT, '(10i3)', IOSTAT=IOS ) NWATER, (IWATER(I),I=1,NWATER)
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:5' )

      ! Read 3 lines of comments
      DO L = 1, 3
         READ( IUNIT, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'modin:6' )
      ENDDO

      !=================================================================
      ! Read in resistances for each surface type (see "depvel.f")
      ! IRI,IRLU,IRAC,IRGSS,IRGSO,IRCLS,IRCLO,IVSMAX 
      !=================================================================
      DO L = 1, NNVEGTYPE
         READ( IUNIT, '(9i5)', IOSTAT=IOS ) 
     &        I, IRI(I),   IRLU(I),  IRAC(I),  IRGSS(I),
     &           IRGSO(I), IRCLS(I), IRCLO(I), IVSMAX(I)

         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IUNIT, 'modin:7' )
      ENDDO

      ! Close the file
      CLOSE( IUNIT )

      ! Return to calling program
      END SUBROUTINE MODIN

!------------------------------------------------------------------------------

      SUBROUTINE RDDRYCF
!
!******************************************************************************
!  Subroutine RDDRYCF read polynomial coefficients from the "drydep.coef"
!  file in the data directory (bmy, 7/6/01, 7/20/04)
!
!  NOTES:
!  (1 ) Use F90 syntax.  Now read "drydep.coef" directly from DATA_DIR.
!        Now use IOERROR to trap I/O errors.  Updated comments and made
!        cosmetic changes (bmy, 7/6/01)
!  (2 ) Removed obsolete code from ages past (bmy, 9/4/01)
!  (3 ) Now read the "drydep.coef" file from the DATA_DIR/drydep_200203/ 
!        directory.  Make IUNIT a dynamic variable and not a parameter. 
!        (bmy, 3/29/02)
!  (4 ) Removed obsolete code from March 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE as the logical unit
!        number. (bmy, 6/27/02)
!  (5 ) Bundled into "drydep_mod.f" (bmy, 11/21/02) 
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

      IMPLICIT NONE
      
#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER            :: I, IOS
      CHARACTER(LEN=80)  :: DUM
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! RDDRYCF begins here!
      !=================================================================

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'drydep_200203/drydep.coef'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - RDDRYCF: Reading ', a ) 

      ! Open file
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &               FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rddrycf:1' )

      ! Read header line
      READ( IU_FILE, '(a80)', IOSTAT=IOS ) DUM
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rddrycf:2' )

      ! Read polynomial coefficients
      READ( IU_FILE,'(8(1pe10.2))', IOSTAT=IOS) (DRYCOEFF(I),I=1,NNPOLY)
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rddrycf:3' )

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE RDDRYCF

!------------------------------------------------------------------------------

      FUNCTION AERO_SFCRSII( K, II, PRESS, TEMP, USTAR, RHB ) RESULT(RS)
!
!******************************************************************************
!  Function AERO_SFCRSII computes the aerodynamic resistance of seasalt aerosol
!  tracers according to Zhang et al 2001.  We account for hygroscopic growth
!  of the seasalt aerosol particles (rjp, tdf, bec, bmy, 4/1/04, 6/11/08)
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) K     (INTEGER) : Dry deposition tracer index (range: 1-NUMDEP)
!  (2 ) II    (INTEGER) : GEOS-CHEM surface type index
!  (3 ) PRESS (REAL*8 ) : Pressure [kPa] (where 1 Kpa = 0.1 mb)
!  (4 ) TEMP  (REAL*8 ) : Temperature [K]
!  (5 ) USTAR (REAL*8 ) : Friction Velocity [m/s]
!  (6 ) RHB   (REAL*8)  : Relative humidity (fraction)
!
!  Function Value
!  ============================================================================
!  (6 ) Rs    (REAL*8 ) : Surface resistance for dust particles [s/m]
!
!  NOTES
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 6/11/08)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: K     ! INDEX OF NUMDEP
      INTEGER, INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
      REAL*8,  INTENT(IN) :: TEMP  ! Temperature (K)    
      REAL*8,  INTENT(IN) :: USTAR ! Friction velocity (m/s)
      REAL*8,  INTENT(IN) :: RHB   ! Relative humidity (fraction)

      ! Function value
      REAL*8              :: RS    ! Surface resistance for particles [s/m]

      ! Local variables
      INTEGER             :: N
      REAL*8, PARAMETER   :: C1 = 0.7674d0,  C2 = 3.079d0, 
     &                       C3 = 2.573d-11, C4 = -1.424d0

      REAL*8, PARAMETER   :: G0 = 9.8D0
      REAL*8, PARAMETER   :: BETA  = 2.d0
      REAL*8, PARAMETER   :: BOLTZ = 1.381d-23  ! Boltzmann constant (J/K)
      REAL*8, PARAMETER   :: E0 = 3.d0
      REAL*8  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL*8  :: DP          ! Diameter of aerosol [um]
      REAL*8  :: PDP         ! Press * Dp      
      REAL*8  :: CONST       ! Constant for settling velocity calculations
      REAL*8  :: SLIP        ! Slip correction factor
      REAL*8  :: VISC        ! Viscosity of air (Pa s)
      REAL*8  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL*8  :: SC, ST      ! Schmidt and Stokes number (nondim)
      REAL*8  :: RHBL        ! Relative humidity local

      REAL*8  :: DIAM, DEN, RATIO_R, RWET, RCM
      REAL*8  :: FAC1, FAC2
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS

!=======================================================================
!   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
!-----------------------------------------------------------------------
!   1 - Evergreen needleleaf trees             Snow/Ice          (12)
!   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
!   3 - Deciduous needleleaf trees             Coniferous forest ( 1)  
!   4 - Deciduous broadleaf trees              Agricultural land ( 7)    
!   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
!   6 - Grass                                  Amazon forest     ( 2)
!   7 - Crops and mixed farming                Tundra            ( 9)
!   8 - Desert                                 Desert            ( 8)
!   9 - Tundra                                 Wetland           (11)
!  10 - Shrubs and interrupted woodlands       Urban             (15)
!  11 - Wet land with plants                   Water             (14)
!  12 - Ice cap and glacier                    
!  13 - Inland water                           
!  14 - Ocean                                  
!  15 - Urban                                  
!=======================================================================      
!     GEOS-CHEM LUC                1, 2, 3, 4, 5, 6, 7  8, 9,10,11
      INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
      INTEGER :: LUC

      !=================================================================
      !   LUC       1,    2,    3,    4,    5,    6,    7,    8,    
      !   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0, 
      !   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54 
      !
      !   LUC       9,   10,   11,   12,   13,   14,   15
      !   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
      !   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
      !=================================================================

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  :: 
     & ALPHA(15) = (/ 1.0d0,  0.6d0,   1.1d0,   0.8d0, 0.8d0,  
     &                1.2d0,  1.2d0,  50.0d0,  50.0d0, 1.3d0, 
     &                2.0d0, 50.0d0, 100.0d0, 100.0d0, 1.5d0  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  ::
     & GAMMA(15) = (/ 0.56d0, 0.58d0, 0.56d0, 0.56d0, 0.56d0, 
     &                0.54d0, 0.54d0, 0.54d0, 0.54d0, 0.54d0, 
     &                0.54d0, 0.54d0, 0.50d0, 0.50d0, 0.56d0 /)

!...A unit is (mm) so multiply by 1.D-3 to (m)
!   LUC       1,    2,    3,    4,    5,    6,    7,    8,     
!   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
!   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,

!   LUC       9,   10,   11,   12,   13,   14,   15
!   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0

      REAL*8  :: A(15,5)

      REAL*8  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,  
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.D-3

      !=================================================================
      !...Ref. Zhang et al., AE 35(2001) 549-560
      !. 
      !...Model theroy
      !    Vd = Vs + 1./(Ra+Rs)
      !      where Vs is the gravitational settling velocity, 
      !      Ra is the aerodynamic resistance above the canopy
      !      Rs  is the surface resistance
      !    Here we calculate Rs only..
      !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
      !      where Eo is an empirical constant ( = 3.)
      !      Ustar is the friction velocity
      !      Collection efficiency from 
      !        Eb,  [Brownian diffusion]
      !        Eim, [Impaction]
      !        Ein, [Interception]
      !      R1 is the correction factor representing the fraction 
      !         of particles that stick to the surface.
      !=======================================================================
      !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
      !         Sc = v/D, v (the kinematic viscosity of air) 
      !                   D (particle brownian diffusivity)
      !         r usually lies between 1/2 and 2/3
      !      Eim is a function of Stokes number, St
      !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
      !          St = Vs * Ustar * Ustar / v  for smooth surface
      !          A is the characteristic radius of collectors.
      !        
      !       1) Slinn (1982)
      !           Eim = 10^(-3/St)          for smooth surface      
      !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
      !       2) Peters and Eiden (1992)
      !           Eim = ( St / ( alpha + St ) )^(beta)
      !                alpha(=0.8) and beta(=2) are constants
      !       3) Giorgi (1986)
      !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
      !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
      !       4) Davidson et al.(1982)
      !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
      !       5) Zhang et al.(2001) used 2) method with alpha varying with
      !          vegetation type and beta equal to 2
      !
      !      Ein = 0.5 * ( Dp / A )^2
      !
      !      R1 (Particle rebound)  = exp(-St^0.5)
      !=================================================================

      ! Particle radius [cm]
      RCM  = A_RADI(K) * 1.d2
      
      ! Exponential factors used for hygroscopic growth
      FAC1 = C1 * ( RCM**C2 )
      FAC2 = C3 * ( RCM**C4 )

      ! Aerosol growth with relative humidity in radius [m] 
      ! (Gerber, 1985) (bec, 12/8/04)
      ! Added safety check for LOG (phs, 6/11/08)
      RHBL    = MAX( TINY(RHB), RHB )
      RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.d0)**0.33d0

      ! Ratio dry over wet radii at the cubic power
      RATIO_R = ( A_RADI(K) / RWET )**3.d0

      ! Diameter of the wet aerosol [m]
      DIAM  = RWET * 2.d0  

      ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
      DEN   = RATIO_R * A_DEN(K) + ( 1.d0 - RATIO_R ) * 1000.d0 

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.d6 
 
      ! Constant for settling velocity calculation       
      CONST = DEN * DIAM**2 * G0 / 18.d0
       
      !=================================================================
      !   # air molecule number density
      !   num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
      !   # gas mean free path
      !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
      !   # Slip correction
      !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
      ! &     / (2. * lamda))) / Dp
      !=================================================================
      ! Note, Slip correction factor calculations following Seinfeld, 
      ! pp464 which is thought to be more accurate but more computation 
      ! required.
      !=================================================================

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP( -0.059d0 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458d-6 * (TEMP)**(1.5d0) / (TEMP + 110.4d0)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928d0  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP 
     &      / (3.d0 * 3.141592d0 * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.D0/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0d0 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * G0 ) ! for smooth surface 
         EIN  = 0D0
      ELSE
         ST   = VTS   * USTAR / ( G0 * AA )          ! for vegetated surfaces
         EIN  = 0.5d0 * ( DIAM / AA )**2
      ENDIF

      EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

      EIM  = MIN( EIM, 0.6D0 )

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.D0
      ELSE
         R1 = EXP( -1D0 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      RS   = 1.D0 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

      ! Return to calling program
      END FUNCTION AERO_SFCRSII

!------------------------------------------------------------------------------

      FUNCTION DUST_SFCRSI( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
!******************************************************************************
!  Function DUST_SFCRSI computes the aerodynamic resistance of dust aerosol
!  tracers according to Seinfeld et al 96.  We do not consider hygroscopic
!  growth of the dust aerosol particles. (rjp, tdf, bmy, bec, 4/1/04, 4/15/05)
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) K     (INTEGER) : Dry deposition tracer index (range: 1-NUMDEP)
!  (2 ) II    (INTEGER) : GEOS-CHEM surface type index
!  (3 ) PRESS (REAL*8 ) : Pressure [kPa] (where 1 Kpa = 0.1 mb)
!  (4 ) TEMP  (REAL*8 ) : Temperature [K]
!  (5 ) USTAR (REAL*8 ) : Friction Velocity [m/s]
!
!  Function Value
!  ============================================================================
!  (6 ) Rs    (REAL*8 ) : Surface resistance for dust particles [s/m]
!
!  NOTES
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!******************************************************************************
!    
      INTEGER, INTENT(IN) :: K     ! INDEX OF NUMDEP
      INTEGER, INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8, INTENT(IN)  :: PRESS ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
      REAL*8, INTENT(IN)  :: TEMP  ! Temperature (K)    
      REAL*8, INTENT(IN)  :: USTAR ! Friction velocity (m/s)

      ! Function value
      REAL*8              :: RS    ! Surface resistance for particles [s/m]

      ! Local variables
      INTEGER             :: N
      REAL*8, PARAMETER   :: C1 = 0.7674d0,  C2 = 3.079d0, 
     &                       C3 = 2.573d-11, C4 = -1.424d0

      REAL*8, PARAMETER   :: G0 = 9.8d0
      REAL*8, PARAMETER   :: BETA  = 2.d0
      REAL*8, PARAMETER   :: BOLTZ = 1.381D-23  ! Baltzmann constant (J/K)
      rEAL*8, PARAMETER   :: E0 = 1.d0
      REAL*8  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL*8  :: DP          ! Diameter of aerosol [um]
      REAL*8  :: PDP         ! Press * Dp      
      REAL*8  :: CONST       ! Constant for settling velocity calculations
      REAL*8  :: SLIP        ! Slip correction factor
      REAL*8  :: VISC        ! Viscosity of air (Pa s)
      REAL*8  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL*8  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL*8  :: DIAM, DEN
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS

      !=================================================================
      ! Ref. Zhang et al., AE 35(2001) 549-560 and Seinfeld(1986)
      ! 
      ! Model theory
      !    Vd = Vs + 1./(Ra+Rs)
      !      where Vs is the gravitational settling velocity, 
      !      Ra is the aerodynamic resistance above the canopy
      !      Rs  is the surface resistance
      !    Here we calculate Rs only..
      !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
      !      where Eo is an empirical constant ( = 3.)
      !      Ustar is the friction velocity
      !      Collection efficiency from 
      !        Eb,  [Brownian diffusion]
      !        Eim, [Impaction]
      !        Ein, [Interception]
      !      R1 is the correction factor representing the fraction 
      !         of particles that stick to the surface.
      !=================================================================
      !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
      !         Sc = v/D, v (the kinematic viscosity of air) 
      !                   D (particle brownian diffusivity)
      !         r usually lies between 1/2 and 2/3
      !      Eim is a function of Stokes number, St
      !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
      !          St = Vs * Ustar * Ustar / v  for smooth surface
      !          A is the characteristic radius of collectors.
      !        
      !       1) Slinn (1982)
      !           Eim = 10^(-3/St)          for smooth surface      
      !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
      !       2) Peters and Eiden (1992)
      !           Eim = ( St / ( alpha + St ) )^(beta)
      !                alpha(=0.8) and beta(=2) are constants
      !       3) Giorgi (1986)
      !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
      !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
      !       4) Davidson et al.(1982)
      !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
      !       5) Zhang et al.(2001) used 2) method with alpha varying with
      !          vegetation type and beta equal to 2
      !
      !      Ein = 0.5 * ( Dp / A )^2
      !
      !      R1 (Particle rebound)  = exp(-St^0.5)
      !=================================================================

      ! Particle diameter [m]
      DIAM  = A_RADI(K) * 2.d0 

      ! Particle density [kg/m3]
      DEN   = A_DEN(K)          

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.d6 
 
      ! Constant for settling velocity calculation       
      CONST = DEN * DIAM**2 * G0 / 18.d0
       
      !=================================================================
      !   # air molecule number density
      !   num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
      !   # gas mean free path
      !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
      !   # Slip correction
      !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
      ! &     / (2. * lamda))) / Dp
      !================================================================
      ! Note, Slip correction factor calculations following Seinfeld, 
      ! pp464 which is thought to be more accurate but more computation 
      ! required.
      !=================================================================

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP( -0.059d0 * PDP ) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458d-6 * (TEMP)**(1.5d0) / (TEMP + 110.4d0)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928d0  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP 
     &     / (3.d0 * 3.141592d0 * VISC * DIAM)  

      ! Schmidt number and Diffusion term
      SC   = AIRVS / DIFF                            
      EB   = SC**(-0.666667d0)

      ! Stokes number and impaction term
      ST   = VTS * USTAR * USTAR / ( AIRVS * G0 )
      EIM  = 10.d0**(-3.d0 / ST) 

      ! surface resistance for particle
      RS   = 1.D0 / ( E0 * USTAR * (EB + EIM) )
      
      ! Return to calling program
      END FUNCTION DUST_SFCRSI

!------------------------------------------------------------------------------
! added MERGE1 (hotp 5/25/09)

      FUNCTION ADUST_SFCRSII( K, II, PRESS, TEMP, USTAR ) RESULT( RS )

! This routine is used for all aerosols except dust, sulfate, and seasalt (hotp 7/31/09)
! modified hotp 7/12/07 for non size resolved aerosols
! this is just DUST_SFCRSII rename and the diameter and density fixed
!
!******************************************************************************
!  Function ADUST_SFCRSII computes the aerodynamic resistance of non-size
!  resolved aerosol according to Zhang et al 2001.  We do not consider the hygroscopic
!  growth of the aerosol particles. (rjp, tdf, bec, bmy, 4/1/04, 4/15/05)
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) K     (INTEGER) : Dry deposition tracer index (range: 1-NUMDEP)
!  (2 ) II    (INTEGER) : GEOS-CHEM surface type index
!  (3 ) PRESS (REAL*8 ) : Pressure [kPa] (where 1 Kpa = 0.1 mb)
!  (4 ) TEMP  (REAL*8 ) : Temperature [K]
!  (5 ) USTAR (REAL*8 ) : Friction Velocity [m/s]
!
!  Function Value
!  ============================================================================
!  (6 ) Rs    (REAL*8 ) : Surface resistance for dust particles [s/m]
!
!  NOTES
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: K     ! INDEX OF NUMDEP
      INTEGER, INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
      REAL*8,  INTENT(IN) :: TEMP  ! Temperature (K)    
      REAL*8,  INTENT(IN) :: USTAR ! Friction velocity (m/s)

      ! Function value
      REAL*8              :: RS    ! Surface resistance for particles [s/m]

      ! Local variables
      INTEGER             :: N
      REAL*8, PARAMETER   :: C1 = 0.7674d0,  C2 = 3.079d0, 
     &                       C3 = 2.573d-11, C4 = -1.424d0

      REAL*8, PARAMETER   :: G0 = 9.8D0
      REAL*8, PARAMETER   :: BETA  = 2.d0
      REAL*8, PARAMETER   :: BOLTZ = 1.381d-23  ! Boltzmann constant (J/K)
      REAL*8, PARAMETER   :: E0 = 3.d0
      REAL*8  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL*8  :: DP          ! Diameter of aerosol [um]
      REAL*8  :: PDP         ! Press * Dp      
      REAL*8  :: CONST       ! Constant for settling velocity calculations
      REAL*8  :: SLIP        ! Slip correction factor
      REAL*8  :: VISC        ! Viscosity of air (Pa s)
      REAL*8  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL*8  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL*8  :: DIAM, DEN
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS

!=======================================================================
!   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
!-----------------------------------------------------------------------
!   1 - Evergreen needleleaf trees             Snow/Ice          (12)
!   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
!   3 - Deciduous needleleaf trees             Coniferous forest ( 1)  
!   4 - Deciduous broadleaf trees              Agricultural land ( 7)    
!   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
!   6 - Grass                                  Amazon forest     ( 2)
!   7 - Crops and mixed farming                Tundra            ( 9)
!   8 - Desert                                 Desert            ( 8)
!   9 - Tundra                                 Wetland           (11)
!  10 - Shrubs and interrupted woodlands       Urban             (15)
!  11 - Wet land with plants                   Water             (14)
!  12 - Ice cap and glacier                    
!  13 - Inland water                           
!  14 - Ocean                                  
!  15 - Urban                                  
!=======================================================================      
!     GEOS-CHEM LUC                1, 2, 3, 4, 5, 6, 7  8, 9,10,11
      INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
      INTEGER :: LUC

!=======================================================================
!   LUC       1,    2,    3,    4,    5,    6,    7,    8,    
!   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0, 
!   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54

!   LUC       9,   10,   11,   12,   13,   14,   15
!   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
!   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
!=======================================================================

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  :: 
     & ALPHA(15) = (/ 1.0d0,  0.6d0,   1.1d0,   0.8d0, 0.8d0,  
     &                1.2d0,  1.2d0,  50.0d0,  50.0d0, 1.3d0, 
     &                2.0d0, 50.0d0, 100.0d0, 100.0d0, 1.5d0  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  ::
     & GAMMA(15) = (/ 0.56d0, 0.58d0, 0.56d0, 0.56d0, 0.56d0, 
     &                0.54d0, 0.54d0, 0.54d0, 0.54d0, 0.54d0, 
     &                0.54d0, 0.54d0, 0.50d0, 0.50d0, 0.56d0 /)

!...A unit is (mm) so multiply by 1.D-3 to (m)
!   LUC       1,    2,    3,    4,    5,    6,    7,    8,     
!   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
!   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,

!   LUC       9,   10,   11,   12,   13,   14,   15
!   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0

      REAL*8  :: A(15,5)

      REAL*8  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,  
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.D-3

      !=================================================================
      !...Ref. Zhang et al., AE 35(2001) 549-560
      !. 
      !...Model theroy
      !    Vd = Vs + 1./(Ra+Rs)
      !      where Vs is the gravitational settling velocity, 
      !      Ra is the aerodynamic resistance above the canopy
      !      Rs  is the surface resistance
      !    Here we calculate Rs only..
      !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
      !      where Eo is an empirical constant ( = 3.)
      !      Ustar is the friction velocity
      !      Collection efficiency from 
      !        Eb,  [Brownian diffusion]
      !        Eim, [Impaction]
      !        Ein, [Interception]
      !      R1 is the correction factor representing the fraction 
      !         of particles that stick to the surface.
      !=======================================================================
      !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
      !         Sc = v/D, v (the kinematic viscosity of air) 
      !                   D (particle brownian diffusivity)
      !         r usually lies between 1/2 and 2/3
      !      Eim is a function of Stokes number, St
      !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
      !          St = Vs * Ustar * Ustar / v  for smooth surface
      !          A is the characteristic radius of collectors.
      !        
      !       1) Slinn (1982)
      !           Eim = 10^(-3/St)          for smooth surface      
      !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
      !       2) Peters and Eiden (1992)
      !           Eim = ( St / ( alpha + St ) )^(beta)
      !                alpha(=0.8) and beta(=2) are constants
      !       3) Giorgi (1986)
      !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
      !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
      !       4) Davidson et al.(1982)
      !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
      !       5) Zhang et al.(2001) used 2) method with alpha varying with
      !          vegetation type and beta equal to 2
      !
      !      Ein = 0.5 * ( Dp / A )^2
      !
      !      R1 (Particle rebound)  = exp(-St^0.5)
      !=================================================================
      
      ! Particle diameter [m] hotp 10/26/07
      DIAM  = 0.5d-6  

      ! Particle density [kg/m3] hotp 10/26/07
      DEN   = 1500          

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.d6 
 
      ! Constant for settling velocity calculation       
      CONST = DEN * DIAM**2 * G0 / 18.d0
       
      !=================================================================
      !   # air molecule number density
      !   num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
      !   # gas mean free path
      !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
      !   # Slip correction
      !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
      ! &     / (2. * lamda))) / Dp
      !=================================================================
      ! Note, Slip correction factor calculations following Seinfeld, 
      ! pp464 which is thought to be more accurate but more computation 
      ! required.
      !=================================================================

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP( -0.059d0 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458d-6 * (TEMP)**(1.5d0) / (TEMP + 110.4d0)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928d0  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP 
     &      / (3.d0 * 3.141592d0 * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.D0/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0d0 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * G0 ) ! for smooth surface 
         EIN  = 0D0
      ELSE
         ST   = VTS   * USTAR / ( G0 * AA )          ! for vegetated surfaces
         EIN  = 0.5d0 * ( DIAM / AA )**2
      ENDIF

      EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

      EIM  = MIN( EIM, 0.6D0 )

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.D0
      ELSE
         R1 = EXP( -1D0 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      RS   = 1.D0 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

      ! Return to calling program
      END FUNCTION ADUST_SFCRSII

!------------------------------------------------------------------------------

      FUNCTION DUST_SFCRSII( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
! NOW ONLY CALLED FOR DUST
!

!******************************************************************************
!  Function DUST_SFCRSII computes the aerodynamic resistance of dust aerosol
!  tracers according to Zhang et al 2001.  We do not consider the hygroscopic
!  growth of the aerosol particles. (rjp, tdf, bec, bmy, 4/1/04, 4/15/05)
!
!  Arguments as Input: 
!  ============================================================================
!  (1 ) K     (INTEGER) : Dry deposition tracer index (range: 1-NUMDEP)
!  (2 ) II    (INTEGER) : GEOS-CHEM surface type index
!  (3 ) PRESS (REAL*8 ) : Pressure [kPa] (where 1 Kpa = 0.1 mb)
!  (4 ) TEMP  (REAL*8 ) : Temperature [K]
!  (5 ) USTAR (REAL*8 ) : Friction Velocity [m/s]
!
!  Function Value
!  ============================================================================
!  (6 ) Rs    (REAL*8 ) : Surface resistance for dust particles [s/m]
!
!  NOTES
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: K     ! INDEX OF NUMDEP
      INTEGER, INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
      REAL*8,  INTENT(IN) :: TEMP  ! Temperature (K)    
      REAL*8,  INTENT(IN) :: USTAR ! Friction velocity (m/s)

      ! Function value
      REAL*8              :: RS    ! Surface resistance for particles [s/m]

      ! Local variables
      INTEGER             :: N
      REAL*8, PARAMETER   :: C1 = 0.7674d0,  C2 = 3.079d0, 
     &                       C3 = 2.573d-11, C4 = -1.424d0

      REAL*8, PARAMETER   :: G0 = 9.8D0
      REAL*8, PARAMETER   :: BETA  = 2.d0
      REAL*8, PARAMETER   :: BOLTZ = 1.381d-23  ! Boltzmann constant (J/K)
      REAL*8, PARAMETER   :: E0 = 3.d0
      REAL*8  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL*8  :: DP          ! Diameter of aerosol [um]
      REAL*8  :: PDP         ! Press * Dp      
      REAL*8  :: CONST       ! Constant for settling velocity calculations
      REAL*8  :: SLIP        ! Slip correction factor
      REAL*8  :: VISC        ! Viscosity of air (Pa s)
      REAL*8  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL*8  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL*8  :: DIAM, DEN
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS

!=======================================================================
!   #  LUC [Zhang et al., 2001]                GEOS-CHEM LUC (Corr. #)
!-----------------------------------------------------------------------
!   1 - Evergreen needleleaf trees             Snow/Ice          (12)
!   2 - Evergreen broadleaf trees              Deciduous forest  ( 4)
!   3 - Deciduous needleleaf trees             Coniferous forest ( 1)  
!   4 - Deciduous broadleaf trees              Agricultural land ( 7)    
!   5 - Mixed broadleaf and needleleaf trees   Shrub/grassland   (10)
!   6 - Grass                                  Amazon forest     ( 2)
!   7 - Crops and mixed farming                Tundra            ( 9)
!   8 - Desert                                 Desert            ( 8)
!   9 - Tundra                                 Wetland           (11)
!  10 - Shrubs and interrupted woodlands       Urban             (15)
!  11 - Wet land with plants                   Water             (14)
!  12 - Ice cap and glacier                    
!  13 - Inland water                           
!  14 - Ocean                                  
!  15 - Urban                                  
!=======================================================================      
!     GEOS-CHEM LUC                1, 2, 3, 4, 5, 6, 7  8, 9,10,11
      INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
      INTEGER :: LUC

!=======================================================================
!   LUC       1,    2,    3,    4,    5,    6,    7,    8,    
!   alpha   1.0,  0.6,  1.1,  0.8,  0.8,  1.2,  1.2, 50.0, 
!   gamma  0.56, 0.58, 0.56, 0.56, 0.56, 0.54, 0.54, 0.54

!   LUC       9,   10,   11,   12,   13,   14,   15
!   alpha  50.0,  1,3,  2.0, 50.0,100.0,100.0,  1.5
!   gamma  0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56
!=======================================================================

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  :: 
     & ALPHA(15) = (/ 1.0d0,  0.6d0,   1.1d0,   0.8d0, 0.8d0,  
     &                1.2d0,  1.2d0,  50.0d0,  50.0d0, 1.3d0, 
     &                2.0d0, 50.0d0, 100.0d0, 100.0d0, 1.5d0  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  ::
     & GAMMA(15) = (/ 0.56d0, 0.58d0, 0.56d0, 0.56d0, 0.56d0, 
     &                0.54d0, 0.54d0, 0.54d0, 0.54d0, 0.54d0, 
     &                0.54d0, 0.54d0, 0.50d0, 0.50d0, 0.56d0 /)

!...A unit is (mm) so multiply by 1.D-3 to (m)
!   LUC       1,    2,    3,    4,    5,    6,    7,    8,     
!   SC1     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
!   SC2     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,
! A SC3     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC4     2.0,  5.0,  5.0, 10.0,  5.0,  5.0,  5.0,-999.,
!   SC5     2.0,  5.0,  2.0,  5.0,  5.0,  2.0,  2.0,-999.,

!   LUC       9,   10,   11,   12,   13,   14,   15
!   SC1   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC2   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
! A SC3   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC4   -999., 10.0, 10.0,-999.,-999.,-999., 10.0
!   SC5   -999., 10.0, 10.0,-999.,-999.,-999., 10.0

      REAL*8  :: A(15,5)

      REAL*8  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   5.0d0,  10.0d0,  5.0d0,  
     &            5.0d0,   5.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0,
     &
     &            2.0d0,   5.0d0,   2.0d0,   5.0d0,  5.0d0,  
     &            2.0d0,   2.0d0, -999.d0, -999.d0, 10.0d0, 
     &           10.0d0, -999.d0, -999.d0, -999.d0, 10.0d0  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.D-3

      !=================================================================
      !...Ref. Zhang et al., AE 35(2001) 549-560
      !. 
      !...Model theroy
      !    Vd = Vs + 1./(Ra+Rs)
      !      where Vs is the gravitational settling velocity, 
      !      Ra is the aerodynamic resistance above the canopy
      !      Rs  is the surface resistance
      !    Here we calculate Rs only..
      !    Rs = 1 / (Eo*Ustar*(Eb+Eim+Ein)*R1)
      !      where Eo is an empirical constant ( = 3.)
      !      Ustar is the friction velocity
      !      Collection efficiency from 
      !        Eb,  [Brownian diffusion]
      !        Eim, [Impaction]
      !        Ein, [Interception]
      !      R1 is the correction factor representing the fraction 
      !         of particles that stick to the surface.
      !=======================================================================
      !      Eb is a funciont of Schmidt number, Eb = Sc^(-gamma)
      !         Sc = v/D, v (the kinematic viscosity of air) 
      !                   D (particle brownian diffusivity)
      !         r usually lies between 1/2 and 2/3
      !      Eim is a function of Stokes number, St
      !          St = Vs * Ustar / (g0 * A)   for vegetated surfaces
      !          St = Vs * Ustar * Ustar / v  for smooth surface
      !          A is the characteristic radius of collectors.
      !        
      !       1) Slinn (1982)
      !           Eim = 10^(-3/St)          for smooth surface      
      !           Eim = St^2 / ( 1 + St^2 ) for vegetative canopies
      !       2) Peters and Eiden (1992)
      !           Eim = ( St / ( alpha + St ) )^(beta)
      !                alpha(=0.8) and beta(=2) are constants
      !       3) Giorgi (1986)
      !           Eim = St^2 / ( 400 + St^2 )     for smooth surface
      !           Eim = ( St / (0.6 + St) )^(3.2) for vegetative surface
      !       4) Davidson et al.(1982)
      !           Eim = St^3 / (St^3+0.753*St^2+2.796St-0.202) for grassland
      !       5) Zhang et al.(2001) used 2) method with alpha varying with
      !          vegetation type and beta equal to 2
      !
      !      Ein = 0.5 * ( Dp / A )^2
      !
      !      R1 (Particle rebound)  = exp(-St^0.5)
      !=================================================================
      
      ! Particle diameter [m]
      DIAM  = A_RADI(K) * 2.d0  

      ! Particle density [kg/m3]
      DEN   = A_DEN(K)          

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.d6 
 
      ! Constant for settling velocity calculation       
      CONST = DEN * DIAM**2 * G0 / 18.d0
       
      !=================================================================
      !   # air molecule number density
      !   num = P * 1d3 * 6.023d23 / (8.314 * Temp) 
      !   # gas mean free path
      !   lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 ) 
      !   # Slip correction
      !   Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp     
      ! &     / (2. * lamda))) / Dp
      !=================================================================
      ! Note, Slip correction factor calculations following Seinfeld, 
      ! pp464 which is thought to be more accurate but more computation 
      ! required.
      !=================================================================

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP( -0.059d0 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458d-6 * (TEMP)**(1.5d0) / (TEMP + 110.4d0)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928d0  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP 
     &      / (3.d0 * 3.141592d0 * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.D0/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0d0 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * G0 ) ! for smooth surface 
         EIN  = 0D0
      ELSE
         ST   = VTS   * USTAR / ( G0 * AA )          ! for vegetated surfaces
         EIN  = 0.5d0 * ( DIAM / AA )**2
      ENDIF

      EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

      EIM  = MIN( EIM, 0.6D0 )

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.D0
      ELSE
         R1 = EXP( -1D0 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      RS   = 1.D0 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

      ! Return to calling program
      END FUNCTION DUST_SFCRSII

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DRYDEP
!
!******************************************************************************
!  Subroutine INIT_DRYDEP initializes certain variables for the GEOS-CHEM
!  dry deposition subroutines. (bmy, 11/19/02, 10/19/09)
!
!  NOTES:
!  (1 ) Added N2O5 as a drydep tracer, w/ the same drydep velocity as
!        HNO3.  Now initialize PBLFRAC array. (rjp, bmy, 7/21/03)
!  (2 ) Added extra carbon & dust aerosol tracers (rjp, tdf, bmy, 4/1/04)
!  (3 ) Added seasalt aerosol tracers.  Now use A_RADI and A_DEN to store
!        radius & density of size-resolved tracers.  Also added fancy
!        output. (bec, rjp, bmy, 4/26/04)
!  (3 ) Now handles extra SOA tracers (rjp, bmy, 7/13/04)
!  (4 ) Now references LDRYD from "logical_mod.f" and N_TRACERS, 
!        SALA_REDGE_um, and SALC_REDGE_um from "tracer_mod.f" (bmy, 7/20/04)
!  (5 ) Included Hg2, HgP tracers (eck, bmy, 12/14/04)
!  (6 ) Included AS, AHS, LET, NH4aq, SO4aq tracers (cas, bmy, 1/6/05)
!  (7 ) Remove reference to PBLFRAC array -- it's obsolete (bmy, 2/22/05)
!  (8 ) Included SO4s, NITs tracers (bec, bmy, 4/13/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now set Henry's law constant to 1.0d+14 for Hg2.  Now use ID_Hg2, 
!        ID_HgP, and ID_Hg_tot from "tracerid_mod.f".  Bug fix: split up
!        compound IF statements into separate 2 IF statements for ID_Hg2, 
!        ID_HgP to avoid seg faults. (eck, cdh, bmy, 4/17/06)
!  (11) Now also initialize SOG4, SOA4 drydep species.  Bug fix: Remove 2nd
!        "IF ( IS_Hg ) THEN" statement. (dkh, bmy, 5/24/06)
!  (12) Bug fix: fix TYPO in IF block for IDTSOA4 (dkh, bmy, 6/23/06)
!  (13) Included H2/HD tracers for offline H2-HD sim (phs, 9/18/07)
!  (14) Add dicarbonyl chemistry species (tmf, ccc, 3/6/09)
!  (15) Minor bug fix: ALPH, LIMO should have molwt = 136.23, not 136 even
!        (bmy, 10/19/09)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE LOGICAL_MOD,  ONLY : LDRYD
      USE TRACER_MOD,   ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,   ONLY : ITS_A_POPS_SIM
      USE TRACER_MOD,   ONLY : N_TRACERS, SALA_REDGE_um, SALC_REDGE_um
      USE TRACERID_MOD, ONLY : IDTPB,     IDTBE7,        IDTNOX
      USE TRACERID_MOD, ONLY : IDTOX,     IDTPAN,        IDTHNO3 
      USE TRACERID_MOD, ONLY : IDTH2O2,   IDTPMN,        IDTPPN  
      USE TRACERID_MOD, ONLY : IDTISN2,   IDTR4N2,       IDTCH2O 
      USE TRACERID_MOD, ONLY : IDTN2O5,   IDTSO2,        IDTSO4  
      USE TRACERID_MOD, ONLY : IDTSO4S,   IDTMSA,        IDTNH3  
      USE TRACERID_MOD, ONLY : IDTNH4,    IDTNIT,        IDTNITS 
      USE TRACERID_MOD, ONLY : IDTAS,     IDTAHS,        IDTLET  
      USE TRACERID_MOD, ONLY : IDTSO4aq,  IDTNH4aq,      IDTBCPI 
      USE TRACERID_MOD, ONLY : IDTOCPI,   IDTBCPO,       IDTOCPO 
      USE TRACERID_MOD, ONLY : IDTALPH,   IDTLIMO,       IDTALCO 
      USE TRACERID_MOD, ONLY : IDTSOG1,   IDTSOG2,       IDTSOG3 
      USE TRACERID_MOD, ONLY : IDTSOG4,   IDTSOA1,       IDTSOA2       
      USE TRACERID_MOD, ONLY : IDTSOA3,   IDTSOA4,       IDTDST1
      ! (hotp 5/25/09)
      USE TRACERID_MOD, ONLY : IDTSOA5,   IDTSOG5
      USE TRACERID_MOD, ONLY : IDTDST2,   IDTDST3,       IDTDST4
      USE TRACERID_MOD, ONLY : IDTSALA,   IDTSALC,       Id_Hg2
      USE TRACERID_MOD, ONLY : ID_HgP,    ID_Hg_tot
      USE TRACERID_MOD, ONLY : ID_Hg0 !eck added 19jul06
      USE TRACERID_MOD, ONLY : IDTPOPG, IDTPOPP
      USE TRACERID_MOD, ONLY : IDTH2,     IDTHD
      USE TRACERID_MOD, ONLY : IDTGLYX,   IDTMGLY
      USE TRACERID_MOD, ONLY : IDTSOAG,   IDTSOAM
      USE TRACERID_MOD, ONLY : IDTGLYC
      USE TRACERID_MOD, ONLY : IDTAPAN, IDTENPAN, IDTGLPAN
      USE TRACERID_MOD, ONLY : IDTGPAN, IDTMPAN, IDTNIPAN
      !add some species (fp, 6/2009)
      USE TRACERID_MOD, ONLY : IDTPROPNN
      USE TRACERID_MOD, ONLY : IDTISOPN    
      USE TRACERID_MOD, ONLY : IDTMMN
      USE TRACERID_MOD, ONLY : IDTRIP, IDTIEPOX, IDTPYPAN
      USE TRACERID_MOD, ONLY : IDTMAP

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      LOGICAL :: IS_Hg
      LOGICAL :: IS_POPS
      INTEGER :: AS, N

      !=================================================================
      ! INIT_DRYDEP begins here!
      !=================================================================

      ! Is this a mercury simulation?
      IS_Hg      = ITS_A_MERCURY_SIM()
      ! Is this a pops simulation?
      IS_POPS     = ITS_A_POPS_SIM()
      ! Zero variables
      DRYDNO2    = 0
      DRYDPAN    = 0
      DRYDHNO3   = 0 
      !(fp, 06/09)
      DRYDH2O2   = 0 
 
      NUMDEP     = 0
      NTRAIND(:) = 0
      NDVZIND(:) = 0
      HSTAR(:)   = 0d0
      F0(:)      = 0d0
      XMW(:)     = 0d0
      A_RADI(:)  = 0d0
      A_DEN(:)   = 0d0
      AIROSOL(:) = .FALSE.

      !=================================================================
      ! First identify tracers that dry deposit and then initialize 
      ! DEPNAME, NDVZIND, HSTAR, F0, XMW and AIROSOL accordingly
      !=================================================================
      DO N = 1, N_TRACERS

         !----------------------------------
         ! Regular full-chemistry tracers
         !----------------------------------

         ! 210Pb (aerosol)
         IF ( N == IDTPB ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPB
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = '210Pb'
            HSTAR(NUMDEP)   = 0.0d+3
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 210d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! 7Be (aerosol)
         ELSE IF ( N == IDTBE7 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTBE7
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = '7Be'
            HSTAR(NUMDEP)   = 0.0d+3
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 7d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! NO2 (as part of NOx)
         ELSE IF ( N == IDTNOX ) THEN
            NUMDEP          = NUMDEP + 1
            DRYDNO2         = NUMDEP
            NTRAIND(NUMDEP) = IDTNOX
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NO2'
            HSTAR(NUMDEP)   = 0.01d0
            F0(NUMDEP)      = 0.1d0
            XMW(NUMDEP)     = 46d-3
            AIROSOL(NUMDEP) = .FALSE.
            
         ! O3 (as part of Ox)
         ELSE IF ( N == IDTOX ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOX
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'O3'
            HSTAR(NUMDEP)   = 0.01d0
            F0(NUMDEP)      = 1.0d0
            XMW(NUMDEP)     = 48d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! PAN 
         ELSE IF ( N == IDTPAN ) THEN
            NUMDEP          = NUMDEP + 1       
            DRYDPAN         = NUMDEP
            NTRAIND(NUMDEP) = IDTPAN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'PAN'
            HSTAR(NUMDEP)   = 3.6d0
            F0(NUMDEP)      = 0.1d0
            XMW(NUMDEP)     = 121d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! HNO3
         ELSE IF ( N == IDTHNO3 ) THEN
            NUMDEP          = NUMDEP + 1
            DRYDHNO3        = NUMDEP
            NTRAIND(NUMDEP) = IDTHNO3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'HNO3'
            HST AR(NUMDEP)  = 1.0d+14
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 63d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! H2O2
         ELSE IF ( N == IDTH2O2 ) THEN
            NUMDEP          = NUMDEP + 1
            ! FP (6/2009)
            DRYDH2O2        = NUMDEP
            NDVZIND(NUMDEP) = NUMDEP
            NTRAIND(NUMDEP) = IDTH2O2
            DEPNAME(NUMDEP) = 'H2O2'
            HSTAR(NUMDEP)   = 1.0d+5
            F0(NUMDEP)      = 1.0d0
            XMW(NUMDEP)     = 34d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! PMN (uses same dep vel as PAN)
         ELSE IF ( N == IDTPMN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPMN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'PMN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! PPN (uses same dep vel as PAN)
         ELSE IF ( N == IDTPPN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPPN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'PPN'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.
         
         ! PYPAN (uses same dep vel as PAN)
         !FP_ISOP
         ELSE IF ( N == IDTPYPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPYPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'PYPAN'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! ISN2 (uses same dep vel as HNO3)
         ELSE IF ( N == IDTISN2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISN2
            NDVZIND(NUMDEP) = DRYDHNO3
            DEPNAME(NUMDEP) = 'ISN2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! R4N2 (uses same dep vel as PAN)
         ELSE IF ( N == IDTR4N2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTR4N2
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'R4N2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! CH2O 
         ELSE IF ( N == IDTCH2O ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTCH2O
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'CH2O'
            HSTAR(NUMDEP)   = 6.0d+3
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 30d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! Add GLYX and MGLY dry deposition, 
         ! using same algorithm as CH2O. (tmf, 5/25/06) 
         ! GLYX 
         ELSE IF ( N == IDTGLYX ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTGLYX
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'GLYX'
            HSTAR(NUMDEP)   = 3.6d+5
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 58d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! MGLY 
         ELSE IF ( N == IDTMGLY ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMGLY
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MGLY'
            HSTAR(NUMDEP)   = 3.7d+3
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 72d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! GLYC
         ELSE IF ( N == IDTGLYC ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTGLYC
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'GLYC'
            HSTAR(NUMDEP)   = 4.1d+4
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 60d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! APAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTAPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTAPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'APAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! ENPAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTENPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTENPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'ENPAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! GLPAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTGLPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTGLPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'GLPAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! GPAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTGPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTGPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'GPAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! MPAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTMPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'MPAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! NIPAN (uses same dep vel as PAN)
         ELSE IF ( N == IDTNIPAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNIPAN
            NDVZIND(NUMDEP) = DRYDPAN
            DEPNAME(NUMDEP) = 'NIPAN'           
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         ! N2O5  (uses same dep vel as HNO3) 
         ELSE IF ( N == IDTN2O5 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTN2O5
            NDVZIND(NUMDEP) = DRYDHNO3
            DEPNAME(NUMDEP) = 'N2O5'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.   

         !FP_ISOP
         ! !ISOPN=ISOPNDISOPNB
         ELSE IF ( N == IDTISOPN ) THEN
         !ISOPNNR
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOPN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOPND'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 0d0  
            XMW(NUMDEP)     = 147d-3
            AIROSOL(NUMDEP) = .FALSE.

            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOPN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOPNB'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 0d0  
            XMW(NUMDEP)     = 147d-3
            AIROSOL(NUMDEP) = .FALSE.

         !MMN=MACRN+MVKN
         ELSE IF ( N == IDTMMN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMMN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MACRN'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 0d0  !TO CHECK
            XMW(NUMDEP)     = 149d-3
            AIROSOL(NUMDEP) = .FALSE.
       
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMMN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MVKN'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 0d0  
            XMW(NUMDEP)     = 149d-3
            AIROSOL(NUMDEP) = .FALSE.   

         !ANIT
         ELSE IF ( N == IDTPROPNN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPROPNN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'PROPNN'  !NITROOXYACETONE IN SANDER TABLE
            HSTAR(NUMDEP)   = 1d3       
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 119d-3
            AIROSOL(NUMDEP) = .FALSE.   

         !RIP
         ELSE IF ( N == IDTRIP ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTRIP
            NDVZIND(NUMDEP) = DRYDH2O2 !USE H2O2
            DEPNAME(NUMDEP) = 'RIP'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         !IEPOX
         ELSE IF ( N == IDTIEPOX ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTIEPOX
            NDVZIND(NUMDEP) = DRYDH2O2 !USE H2O2
            DEPNAME(NUMDEP) = 'IEPOX'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 0d0
            AIROSOL(NUMDEP) = .FALSE.

         !MAP
         ELSE IF ( N == IDTMAP ) THEN
            NUMDEP          = NUMDEP + 1
            NDVZIND(NUMDEP) = NUMDEP
            NTRAIND(NUMDEP) = IDTMAP
            DEPNAME(NUMDEP) = 'MAP'
            HSTAR(NUMDEP)   = 8.4d+2 !FROM R. Sander
            F0(NUMDEP)      = 1.0d0  !Assume reactive
            XMW(NUMDEP)     = 76d-3
            AIROSOL(NUMDEP) = .FALSE.

         !----------------------------------
         ! Sulfur & Nitrate aerosol tracers
         !----------------------------------

         ! SO2
         ELSE IF ( N == IDTSO2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSO2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SO2'
            HSTAR(NUMDEP)   = 1.0d+5
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 64d-3
            AIROSOL(NUMDEP) = .FALSE. 

         ! SO4 (aerosol)
         ELSE IF ( N == IDTSO4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSO4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SO4'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 96d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! SO4 in seasalt aerosol (bec, bmy, 4/13/05)
         ELSE IF ( N == IDTSO4s ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSO4s
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SO4S'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 36d-3                   ! MW of seasalt
            A_RADI(NUMDEP)  = ( SALC_REDGE_um(1) + 
     &                          SALC_REDGE_um(2) ) * 0.5d-6
            A_DEN(NUMDEP)   = 2200.d0 
            AIROSOL(NUMDEP) = .TRUE. 

         ! MSA (aerosol)
         ELSE IF ( N == IDTMSA ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMSA
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MSA'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 96d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! NH3 
         ELSE IF ( N == IDTNH3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNH3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NH3'
            HSTAR(NUMDEP)   = 2.0d+4
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 17d-3
            AIROSOL(NUMDEP) = .FALSE. 

         ! NH4 (aerosol)
         ELSE IF ( N == IDTNH4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNH4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NH4'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 18d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! NIT (aerosol)
         ELSE IF ( N == IDTNIT ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNIT
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NIT'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 62d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! NIT in seasalt aerosol (bec, bmy, 4/13/05)
         ELSE IF ( N == IDTNITs ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNITs
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NITS'
            HSTAR(NUMDEP)   = 0.0d0 
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 36d-3                   ! MW of seasalt
            A_RADI(NUMDEP)  = ( SALC_REDGE_um(1) + 
     &                          SALC_REDGE_um(2) ) * 0.5d-6
            A_DEN(NUMDEP)   = 2200.d0 
            AIROSOL(NUMDEP) = .TRUE. 

         !----------------------------------
         ! Crystalline & aqueous aerosols
         !----------------------------------            

         ! AS (crystalline ammonium sulfate)
         ELSE IF ( N == IDTAS ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTAS
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'AS'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 132d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! AHS (crystaline ammonium bisulfite) 
         ELSE IF ( N == IDTAHS ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTAHS
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'AHS'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 115d-3
            AIROSOL(NUMDEP) = .TRUE. 
        
         ! LET (crystaline LETOVOCITE)
         ELSE IF ( N == IDTLET  ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTLET
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'LET'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 248.0d-3 
            AIROSOL(NUMDEP) = .TRUE. 

         ! SO4aq (aqueous sulfate aerosol) 
         ELSE IF ( N == IDTSO4aq ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSO4aq
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SO4aq'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 96.0d-3
            AIROSOL(NUMDEP) = .TRUE. 

         ! NH4aq (aqueous NH4 aerosol)
         ELSE IF ( N == IDTNH4aq ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNH4aq
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NH4aq'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 18d-3
            AIROSOL(NUMDEP) = .TRUE. 

         !----------------------------------
         ! Carbon & SOA aerosol tracers
         !----------------------------------

         ! Hydrophilic BC (aerosol)
         ELSE IF ( N == IDTBCPI ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTBCPI
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'BCPI'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! Hydrophilic OC (aerosol)
         ELSE IF ( N == IDTOCPI ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOCPI
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OCPI'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! Hydrophobic BC (aerosol)
         ELSE IF ( N == IDTBCPO ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTBCPO
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'BCPO'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! Hydrophobic OC (aerosol)
         ELSE IF ( N == IDTOCPO ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOCPO
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OCPO'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! ALPH (Alpha-pinene)
         ELSE IF ( N == IDTALPH ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTALPH
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ALPH'
            HSTAR(NUMDEP)   = 0.023d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 136.23d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! LIMO (Limonene)
         ELSE IF ( N == IDTLIMO ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTLIMO
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'LIMO'
            HSTAR(NUMDEP)   = 0.07d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 136.23d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! ALCO (Alcohols)
         ELSE IF ( N == IDTALCO ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTALCO
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ALCO'
            HSTAR(NUMDEP)   = 54.d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 142d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOG1
         ELSE IF ( N == IDTSOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOG1'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOG2
         ELSE IF ( N == IDTSOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOG2'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 160d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOG3
         ELSE IF ( N == IDTSOG3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOG3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOG3'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 220d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOG4
         ELSE IF ( N == IDTSOG4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOG4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOG4'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 130d-3
            AIROSOL(NUMDEP) = .FALSE.

        ! SOG5  (dkh, 03/27/07)  
         ELSE IF ( N == IDTSOG5 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOG5
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOG5'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            ! MWT is 150 not 130 hotp
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOA1
         ELSE IF ( N == IDTSOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOA1'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOA2
         ELSE IF ( N == IDTSOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOA2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 160d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOA3
         ELSE IF ( N == IDTSOA3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOA3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOA3'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 220d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOA4
         ELSE IF ( N == IDTSOA4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOA4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOA4'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 130d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOA5 (dkh, 03/27/07)  
         ELSE IF ( N == IDTSOA5 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOA5
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOA5'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            ! hotp 5/24/09 MWT fix
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
 
         ! SOAG
         ELSE IF ( N == IDTSOAG ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOAG
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOAG'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 58d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOAM
         ELSE IF ( N == IDTSOAM ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSOAM
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SOAM'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 72d-3
            AIROSOL(NUMDEP) = .TRUE.

         !----------------------------------
         ! Dust aerosol tracers
         !----------------------------------

         ! DUST1 (aerosol)
         ELSE IF ( N == IDTDST1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTDST1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'DST1'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 29d-3
            A_RADI(NUMDEP)  = 0.73d-6
            A_DEN(NUMDEP)   = 2500.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! DUST2 (aerosol)
         ELSE IF ( N == IDTDST2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTDST2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'DST2'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 29d-3
            A_RADI(NUMDEP)  = 1.4d-6
            A_DEN(NUMDEP)   = 2650.d0   
            AIROSOL(NUMDEP) = .TRUE.

         ! DUST3 (aerosol)
         ELSE IF ( N == IDTDST3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTDST3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'DST3'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 29d-3
            A_RADI(NUMDEP)  = 2.4d-6
            A_DEN(NUMDEP)   = 2650.d0  
            AIROSOL(NUMDEP) = .TRUE.

         ! DUST4 (aerosol)
         ELSE IF ( N == IDTDST4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTDST4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'DST4'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 29d-3
            A_RADI(NUMDEP)  = 4.5d-6
            A_DEN(NUMDEP)   = 2650.d0   
            AIROSOL(NUMDEP) = .TRUE.

         !----------------------------------
         ! Sea salt aerosol tracers
         !----------------------------------

         ! Accum mode seasalt (aerosol) 
         ELSE IF ( N == IDTSALA ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSALA
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SALA'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 36d-3     
            A_RADI(NUMDEP)  = ( SALA_REDGE_um(1) + 
     &                          SALA_REDGE_um(2) ) * 0.5d-6
            A_DEN(NUMDEP)   = 2200.d0         
            AIROSOL(NUMDEP) = .TRUE. 

         ! Coarse mode seasalt (aerosol) 
         ELSE IF ( N == IDTSALC ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTSALC
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'SALC'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 36d-3 
            A_RADI(NUMDEP)  = ( SALC_REDGE_um(1) + 
     &                          SALC_REDGE_um(2) ) * 0.5d-6
            A_DEN(NUMDEP)   = 2200.d0         
            AIROSOL(NUMDEP) = .TRUE. 
         !----------------------------------
	 ! H2/HD tracers 
         ! (hup, jaegle, phs, 9/17/08)
         !----------------------------------
	 ELSE IF ( N == IDTH2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTH2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'H2'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 0d-3
            AIROSOL(NUMDEP) = .FALSE.

	 ELSE IF ( N == IDTHD ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTHD
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'HD'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 0d-3
            AIROSOL(NUMDEP) = .FALSE.

         !----------------------------------
         ! Mercury tracers
         !----------------------------------

         ELSE IF ( IS_Hg ) THEN
             ! add dry dep of Hg0 (eck, 19jul06)
             ! Hg0 -- Elemental Mercury
             IF ( N == ID_Hg0(ID_Hg_tot) ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = ID_Hg0(ID_Hg_tot)
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'Hg0'
               HSTAR(NUMDEP)   = 0.11
               ! F0 consistent with Lin et al (2006)
               F0(NUMDEP)      = 1.0d-5
               XMW(NUMDEP)     = 201d-3
               AIROSOL(NUMDEP) = .FALSE. 
            ENDIF

            ! Hg2 -- Divalent Mercury
            IF ( N == ID_Hg2(ID_Hg_tot) ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = ID_Hg2(ID_Hg_tot)
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'Hg2'
               HSTAR(NUMDEP)   = 1.0d+6
               F0(NUMDEP)      = 0.0d0
               XMW(NUMDEP)     = 201d-3
               AIROSOL(NUMDEP) = .FALSE. 
            ENDIF

            ! HgP -- Particulate Mercury
            IF ( N == ID_HgP(ID_Hg_tot) ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = ID_HgP(ID_Hg_tot)
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'HgP'
               HSTAR(NUMDEP)   = 0.0d0
               F0(NUMDEP)      = 0.0d0
               XMW(NUMDEP)     = 201d-3
               AIROSOL(NUMDEP) = .TRUE. 
            ENDIF
         
         
         !----------------------------------
         ! POPS tracers
         !----------------------------------

          ELSE IF ( IS_POPS ) THEN
             ! POPs GASEOUS
             IF ( N == IDTPOPG) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = IDTPOPG
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'POPG'
               HSTAR(NUMDEP)   = 2.352d1
               ! Adding octanol-air partition coefficient for POPs            clf 1/3/2011
               ! to account for accumulation in leaf cuticles
               ! Needs to be in units of mol/liter/atm as with HSTAR
               ! Divide unitless Koa at 298 K by product of R (0.0821 atm/M/K) and T (298 K)
               ! Currently set for phenanthrene               
               KOA(NUMDEP)     = 1.78d6
               ! F0 zero for now   clf, 1/3/2011
               F0(NUMDEP)      = 0.0d0
               ! Need to change molecular weight for different POPs
               ! Currently set for phenanthrene (kg/mol)
               XMW(NUMDEP)     = 178.23d-3
               AIROSOL(NUMDEP) = .FALSE. 
            ENDIF

             ! POPs PARTICLES
            IF ( N == IDTPOPP ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = IDTPOPP
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'POPP'
               HSTAR(NUMDEP)   = 0.0d0
               ! Koa for particulate POPs is just set to the equivalent of Henry's Law
               ! so that cuticular accumulation is not considered 
               KOA(NUMDEP)     = 0.0d0
               F0(NUMDEP)      = 0.0d0
               ! Need to change molecular weight for different POPs
               ! Currently set for phenanthrene (kg/mol)
               XMW(NUMDEP)     = 178.23d-3
               AIROSOL(NUMDEP) = .TRUE. 
            ENDIF
          ENDIF

         
      ENDDO
      
      !=================================================================
      ! Additional variables required for pops simultion (following hg)
      ! Locate the drydep species w/in the DEPSAV array
      !=================================================================
      IF ( IS_POPS ) THEN

         ! Initialize flags
         ! add dry dep of POPG and POPP
         DRYPOPG = 0
         DRYPOPP = 0
         
         ! If drydep is turned on ...
         IF ( LDRYD ) THEN
         
            ! Loop over drydep species
            DO N = 1, NUMDEP

               ! Locate by DEPNAME
               SELECT CASE ( TRIM( DEPNAME(N) ) )
                  CASE( 'POPG' )
                     DRYPOPG = N
                  CASE( 'POPP' )
                     DRYPOPP = N
                  CASE DEFAULT
                     ! nothing
                END SELECT
             ENDDO

          ENDIF

       ENDIF
       
      !=================================================================
      ! Additional variables required for mercury simultion
      ! Locate the drydep species w/in the DEPSAV array
      !=================================================================
      IF ( IS_HG ) THEN

         ! Initialize flags
         ! add dry dep of Hg0
         DRYHg0 = 0
         DRYHg2 = 0
         DRYHgP = 0
         
         ! If drydep is turned on ...
         IF ( LDRYD ) THEN
         
            ! Loop over drydep species
            DO N = 1, NUMDEP

               ! Locate by DEPNAME
               SELECT CASE ( TRIM( DEPNAME(N) ) )
                  ! add dry dep of Hg(0)
                  CASE( 'Hg0' )
                     DRYHg0 = N
                  CASE( 'Hg2' )
                     DRYHg2 = N
                  CASE( 'HgP' )
                     DRYHgP = N
                  CASE DEFAULT
                     ! nothing
                END SELECT
             ENDDO

          ENDIF

      ENDIF



      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( DEPSAV( IIPAR, JJPAR, NUMDEP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEPSAV' )
      DEPSAV = 0d0

      !=================================================================
      ! Echo information to stdout
      !=================================================================
      WRITE( 6, '(/,a)' ) 'INIT_DRYDEP: List of dry deposition species:'
      WRITE( 6, '(/,a)'   )
     & '  #   Name  Tracer DEPVEL Henry''s    React.   Molec.  Aerosol?'
      WRITE( 6, '(a)'   )
     & '            Number Index  Law Const  Factor   Weight  (T or F)'
      WRITE( 6, '(a)'   ) REPEAT( '-', 65 )

      DO N = 1, NUMDEP
         WRITE( 6, 100 ) N, TRIM( DEPNAME(N) ), NTRAIND(N), NDVZIND(N), 
     &                   HSTAR(N),  F0(N),      XMW(N),     AIROSOL(N)
      ENDDO
 100  FORMAT( i3, 3x, a4, 2(3x,i3), 4x, es8.1, 2(3x,f6.3), 3x, L3 )

      ! Return to calling program
      END SUBROUTINE INIT_DRYDEP

!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_DRYDEP
!
!******************************************************************************
!  Subroutine CLEANUP_DRYDEP deallocates all module arrays.
!  (bmy, 2/27/03, 2/22/05)
! 
!  NOTES:
!  (1 ) Remove reference to PBLFRAC array; it's obsolete (bmy, 2/22/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DRYDEP begins here!
      !=================================================================
      IF ( ALLOCATED( DEPSAV  ) ) DEALLOCATE( DEPSAV  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DRYDEP

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DRYDEP_MOD

      
