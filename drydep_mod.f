! $Id: drydep_mod.f,v 1.6 2003/12/11 21:54:09 bmy Exp $
      MODULE DRYDEP_MOD
!
!******************************************************************************
!  Module DRYDEP_MOD contains variables and routines for the GEOS-CHEM dry
!  deposition scheme. (bmy, 1/27/03, 12/9/03)
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
!  (31) DEPNAME  (CHAR*14) : Names of dry deposition species
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_DRYDEP          : Dry deposition driver routine
!  (2 ) DVZ_MINVAL         : Sets minimum drydep velocities for SULFATE tracers
!  (3 ) METERO             : Computes meterological fields for dry deposition
!  (4 ) DRYFLX             : Applies drydep losses from SMVGEAR to tracer array
!  (5 ) DRYFLXRnPbBe       : Applies drydep losses to 210Pb and 7Be 
!  (6 ) DEPVEL             : Computes dry deposition velocities (by D. Jacob)
!  (7 ) MODIN              : Reads inputs for DEPVEL from "drydep.table"
!  (8 ) RDDRYCF            : Reads drydep polynomial coeffs from "drydep.coef"
!  (9 ) INIT_DRYDEP        : Initializes and allocates module arrays
!  (10) CLEANUP_DRYDEP     : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "drydep_mod.f":
!  ============================================================================
!  (1 ) comode_mod.f       : Module containing SMVGEAR allocatable arrays
!  (2 ) dao_mod.f          : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f         : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) error_mod.f        : Module containing NaN, other error check routines
!  (5 ) file_mod.f         : Module containing file unit #'s and error checks
!  (6 ) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!  (7 ) tracerid_mod.f     : Module containing pointers to tracers & emissions
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
!
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS
      !=================================================================
      
      ! Module variables
      PRIVATE :: NNTYPE,   NNPOLY,   NNVEGTYPE, XCKMAN
      PRIVATE :: DRYDPAN,  DRYDHNO3, DRYDNO2,   NWATER
      PRIVATE :: AIROSOL,  NDVZIND,  IDEP,      IZO       
      PRIVATE :: IWATER,   IRI,      IRLU,      IRAC      
      PRIVATE :: IRGSS,    IRGSO,    IRCLS,     IRCLO
      PRIVATE :: IVSMAX,   DRYCOEFF, HSTAR,     F0
      PRIVATE :: XMW

      ! Module Routines
      PRIVATE :: DEPVEL,   METERO,   MODIN,     RDDRYCF

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER   :: MAXDEP    = 16     
      INTEGER, PARAMETER   :: NNTYPE    = 15     ! NTYPE    from "CMN_SIZE"
      INTEGER, PARAMETER   :: NNPOLY    = 20     ! NPOLY    from "CMN_SIZE"
      INTEGER, PARAMETER   :: NNVEGTYPE = 74     ! NVEGTYPE from "CMN_SIZE"
      REAL*8,  PARAMETER   :: XCKMAN    = 0.4d0
 
      ! Scalars
      INTEGER              :: DRYDHNO3, DRYDNO2, DRYDPAN
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
      REAL*8,  ALLOCATABLE :: PBLFRAC(:,:,:)
      REAL*8               :: DRYCOEFF(NNPOLY)
      REAL*8               :: HSTAR(MAXDEP)
      REAL*8               :: F0(MAXDEP)
      REAL*8               :: XMW(MAXDEP)
      CHARACTER(LEN=14)    :: DEPNAME(MAXDEP)

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
!  (lwh, gmg, djj, 1989, 1994; bmy, 2/11/03, 12/9/03)
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
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      !--------------------------------------------------------------
      ! Prior to 12/9/03:
      ! Now declare AZO and USTAR as local variables
      !USE DAO_MOD,      ONLY : AD, ALBD, AZO=>Z0, SUNCOS, T, USTAR
      !--------------------------------------------------------------
      USE DAO_MOD,      ONLY : AD, ALBD, SUNCOS, T
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TRACERID_MOD
      
#     include "CMN_SIZE" ! Size parameters
#     include "CMN"      ! XTRA2
#     include "CMN_DIAG" ! ND44
#     include "CMN_DEP"  ! IREG, ILAND, IUSE, etc.

      ! Local variables
      LOGICAL, SAVE     :: FIRST = .TRUE.
      LOGICAL           :: LSNOW(MAXIJ)
      INTEGER           :: I, J, L, N, M, IJLOOP, LAYBOT, NN, NDVZ
      REAL*8            :: PS, X25, RESIDU, P1, P2, THIK, DVZ, PL1, PL2 
      REAL*8            :: HEIGHT(LLPAR),  CZ1(MAXIJ) 
      REAL*8            :: TC0(MAXIJ),     DVEL(MAXIJ,MAXDEP)
      REAL*8            :: ZH(MAXIJ),      OBK(MAXIJ)
      REAL*8            :: CFRAC(MAXIJ),   RADIAT(MAXIJ)
      REAL*8            :: USTAR(MAXIJ),   AZO(MAXIJ)

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
      CALL METERO( CZ1, TC0, OBK, CFRAC, RADIAT, AZO, USTAR )

      !=================================================================
      ! Compute mixing heights
      !=================================================================

      ! 1-D grid box index
      IJLOOP = 0

      ! Compute mixing heights
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Increment IJLOOP
         IJLOOP        = IJLOOP + 1

         ! Set logical LSNOW if snow and sea ice (ALBEDO > 0.4)
         LSNOW(IJLOOP) = ( ALBD(I,J) > 0.4 )

         ! Calculate height of mixed layer [m]
         ZH(IJLOOP)    = 0.0d0
         X25           = XTRA2(I,J)
         LAYBOT        = INT(X25)
           
         ! Compute height of each level that is fully
         ! included w/in the mixed layer [m]
         DO L = 1, LAYBOT + 1
            PL1        = GET_PEDGE(I,J,L)
            PL2        = GET_PEDGE(I,J,L+1)
            HEIGHT(L)  = 2.9271d+01 * T(I,J,L) * LOG( PL1 / PL2 )
         ENDDO

         ! Sum the mixed layer height up to level LAYBOT [m]
         DO M = 1, LAYBOT
            ZH(IJLOOP) = ZH(IJLOOP) + HEIGHT(M)
         ENDDO

         ! For the level where the mixed layer ends, compute the
         ! fractional height of that level w/in the mixed layer [m]
         RESIDU        = X25 - LAYBOT
         ZH(IJLOOP)    = ZH(IJLOOP) + RESIDU * HEIGHT(LAYBOT+1)

         ! If the PBL top occurs in the first level, then set 
         ! RESIDU=1 to apply drydep throughout the surface layer
         IF ( LAYBOT == 0 ) RESIDU = 1d0

         ! Loop up to tropopause
         DO L = 1, LLTROP

            ! Test for level ...
            IF ( L <= LAYBOT ) THEN

               ! Below PBL top: apply drydep frequency to entire level
               PBLFRAC(I,J,L) = 1d0

            ELSE IF ( L == LAYBOT+1 ) THEN

               ! Layer LAYBOT+1 is where PBL top occurs.  Only apply drydep
               ! frequency to the fraction of this level below PBL top.
               PBLFRAC(I,J,L) = RESIDU

            ELSE

               ! Above PBL top: there is no drydep
               PBLFRAC(I,J,L) = 0d0

            ENDIF
         ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! Call DEPVEL to compute dry deposition velocities [m/s]
      !=================================================================
      CALL DEPVEL( MAXIJ, RADIAT, TC0, SUNCOS, F0, HSTAR, XMW,  AIROSOL,
     &             USTAR, CZ1,    OBK, CFRAC,  ZH, LSNOW, DVEL, AZO ) 
      
      !=================================================================
      ! Compute dry deposition frequencies; archive diagnostics
      !=================================================================

      ! 1-D grid box index
      IJLOOP = 0

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Increment IJLOOP
         IJLOOP  = IJLOOP + 1

         ! THIK = thickness of surface layer [m]
         P1      = GET_PEDGE(I,J,1)
         P2      = GET_PEDGE(I,J,2)
         THIK    = 2.9271d+01 * T(I,J,1) * LOG( P1 / P2 )    

         ! Now we calculate drydep throughout the entire PBL.  
         ! Make sure that the PBL depth is greater than or equal 
         ! to the thickness of the 1st layer (rjp, bmy, 7/21/03)
         THIK    = MAX( ZH(IJLOOP), THIK )

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

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )

      ! Return to calling program
      END SUBROUTINE DO_DRYDEP

!------------------------------------------------------------------------------

      FUNCTION DVZ_MINVAL( N, LSNOW, DVZ ) RESULT( NEWDVZ )
!
!******************************************************************************
!  Function DVZ_MINVAL sets minimum values for drydep velocities for SULFATE
!  TRACERS, according to Mian Chin's GOCART model (rjp, bmy, 11/21/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N     (INTEGER) : Tracer number
!  (2 ) LSNOW (LOGICAL) : Flag for denoting snow/ice 
!  (3 ) DVZ   (REAL*8 ) : Deposition velocity [cm/s]
!
!  NOTES:
!  (1 ) Don't put a min drydep value on H2O2 for offline run (rjp, bmy,3/31/03)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! NSRCX

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

      SUBROUTINE METERO( CZ1, TC0, OBK, CFRAC, RADIAT, AZO, USTR )
!
!******************************************************************************
!  Subroutine METERO calculates meteorological constants needed for the      
!  dry deposition velocity module. (lwh, gmg, djj, 1989, 1994; bmy, 12/9/03)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) CZ1    (REAL*8) : Midpoint height of first model level   [m]
!  (2 ) TC0    (REAL*8) : Array for grid box surface temperature [K]
!  (3 ) OBK    (REAL*8) : Array for the Monin-Obhukov length     [m]
!  (4 ) CFRAC  (REAL*8) : Array for the column cloud fraction    [unitless]
!  (5 ) RADIAT (REAL*8) : Array for the solar radiation @ ground [W/mw]
!  (6 ) AZO    (REAL*8) : Array for the roughness heights        [m]
!  (7 ) USTR   (REAL*8) : Array for the friction velocity        [m/s]
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
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : AIRDEN, CLDFRC, HFLUX, RADSWG,
     &                         T,      TS,     USTAR, Z0
      USE PRESSURE_MOD, ONLY : GET_PEDGE
                                  
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants

      ! Arguments
      REAL*8, INTENT(OUT) :: CZ1(MAXIJ)
      REAL*8, INTENT(OUT) :: TC0(MAXIJ)
      REAL*8, INTENT(OUT) :: OBK(MAXIJ)
      REAL*8, INTENT(OUT) :: CFRAC(MAXIJ)
      REAL*8, INTENT(OUT) :: RADIAT(MAXIJ)
      REAL*8, INTENT(OUT) :: AZO(MAXIJ)
      REAL*8, INTENT(OUT) :: USTR(MAXIJ)

      ! Local variables
      INTEGER             :: I,  J,  IJLOOP
      REAL*8              :: P1, P2, THIK, NUM, DEN
      REAL*8, PARAMETER   :: KAPPA = 0.4d0 
      REAL*8, PARAMETER   :: CP    = 1000.0d0

      ! External functions
      REAL*8, EXTERNAL    :: XLTMMP

      !=================================================================
      ! METERO begins here!
      !=================================================================
      IJLOOP = 0

      ! Loop over surface grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D grid box index
         IJLOOP = IJLOOP + 1

         ! THIK = thickness of layer 1 [m]
         P1     = GET_PEDGE(I,J,1)
         P2     = GET_PEDGE(I,J,2)
         THIK   = 2.9271d+01 * T(I,J,1) * LOG( P1 / P2 )

         ! Midpoint height of first model level [m]
         CZ1(IJLOOP) = THIK / 2.0d0

         ! Surface temperature [K]
         TC0(IJLOOP) = TS(I,J)

         !==============================================================
         ! The direct computation of the Monin-Obhukov length is:
         !
         !            - Air density * Cp * T(surface air) * Ustar^3 
         !    OBK =  -----------------------------------------------
         !              Kappa       * g  * Sensible Heat flux
         !
         ! Cp    = 1000 J / kg / K = specific heat of air at constant P
         ! Kappa = 0.4             = Von Karman's constant
         !
         !
         !  Also test the denominator in order to prevent div by zero.
         !=================================================================

         ! Numerator
         NUM = -AIRDEN(1,I,J) *  CP            * TS(I,J) *
!-------------------------------------------------------------------------- 
! Prior to 12/9/03:
! USTAR is now a 2-D array (bmy, 12/9/03)
!     &          USTAR(IJLOOP) *  USTAR(IJLOOP) * USTAR(IJLOOP)
!-------------------------------------------------------------------------- 
     &          USTAR(I,J)    *  USTAR(I,J)    * USTAR(I,J)

         ! Denominator
         DEN =  KAPPA * g0 * HFLUX(I,J) 

         ! Prevent div by zero
         IF ( ABS( DEN ) > 0d0 ) THEN
            OBK(IJLOOP) = NUM / DEN
         ELSE
            OBK(IJLOOP) = 1.0d5
         ENDIF

         !=================================================================
         ! Return meterological quantities as 1-D arrays for DEPVEL
         !=================================================================

         ! Roughness height [m]
         AZO(IJLOOP)    = Z0(I,J)

         ! Column cloud fraction [unitless]
         CFRAC(IJLOOP)  = CLDFRC(I,J)

         ! Solar insolation @ ground [W/m2]
         RADIAT(IJLOOP) = RADSWG(I,J) 
         
         ! Friction velocity [m/s]
         USTR(IJLOOP)   = USTAR(I,J)

      ENDDO
      ENDDO
      
      ! Return to calling program
      END SUBROUTINE METERO

!------------------------------------------------------------------------------

      SUBROUTINE DRYFLX
!
!******************************************************************************
!  Subroutine DRYFLX sets up the dry deposition flux diagnostic for tracers
!  which are part of the SMVGEAR mechanism. (bmy, bdf, 4/20/99, 7/21/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : CSPEC, JLOP, VOLUME
      USE DIAG_MOD,   ONLY : AD44
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE GRID_MOD,   ONLY : GET_AREA_CM2
      USE TIME_MOD,   ONLY : GET_TS_CHEM

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! STT, many other variables
#     include "CMN_DIAG"   ! Diagnostic switches & arrays
#     include "comode.h"   ! CSPEC

      ! Local variables
      INTEGER :: I, J, JJ, JLOOP, L, L_PBLTOP, N, NK, NN
      REAL*8  :: DTCHEM, TDRYFX, AREA_CM2(JJPAR)

      !=================================================================
      ! DRYFLX begins here!
      !=================================================================

      ! Return unless we have turned on ND44 drydep diagnostic
      IF ( ND44 == 0 ) RETURN

      ! There is only drydep in the surface layer, which
      ! is accounted for in the "URBAN" chemistry slot
      NCS    = NCSURBAN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

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
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         
            ! Only deal w/ boxes w/in the boundary layer
            IF ( PBLFRAC(I,J,L) > 0d0 ) THEN

               ! 1-D grid box index for CSPEC & VOLUME
               JLOOP = JLOP(I,J,L)

               ! Dry dep flux [molec] for species N = 
               !  CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
               !  [molec/cm3]     * [cm3]  
               TDRYFX = CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
                    
               ! Convert TDRYFX from [molec] to [molec/cm2/s]        
               TDRYFX = TDRYFX / ( AREA_CM2(J) * DTCHEM ) 
                    
               ! Save into AD44 diagnostic array
               AD44(I,J,N,1) = AD44(I,J,N,1) + TDRYFX
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
!  (hyl, bmy, bdf, 4/2/99, 12/2/03)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,  ONLY : AD44
      USE ERROR_MOD, ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE TIME_MOD,  ONLY : GET_TS_CHEM
      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT
#     include "CMN_DIAG"  ! Diagnostic switches & arrays
#     include "CMN_DEP"   ! Dry deposition variables
#     include "CMN_SETUP" ! LDRYD

      ! Local variables
      INTEGER             :: I, J, L, L_PBLTOP, N, NN
      REAL*8              :: DTCHEM, FRACLOST, AMT_LOST, RESIDU

      !=================================================================
      ! DRYFLXRnPbBe begins here!!
      !=================================================================

      ! Chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Only do the following if DRYDEP is turned on
      IF ( LDRYD ) THEN

         ! Loop over drydep species
         DO N = 1, NUMDEP

            ! Tracer index in STT that corresponds to drydep species N
            ! If invalid, then cycle
            NN = NTRAIND(N)
            IF ( NN == 0 ) CYCLE

            ! Loop over surface grid boxes
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! L_PBLTOP is the level where the PBL top occurs
               L_PBLTOP = INT( XTRA2(I,J) ) + 1

               ! Loop up to the PBL top
               DO L = 1, L_PBLTOP

                  ! FRACLOST is the fraction of tracer lost.  PBLFRAC is 
                  ! the fraction of layer L located totally w/in the PBL.
                  FRACLOST = DEPSAV(I,J,N) * PBLFRAC(I,J,L) * DTCHEM

#if defined( GEOS_1 ) || defined( GEOS_STRAT )
                  !=====================================================
                  ! GEOS-1 or GEOS-STRAT:
                  ! ---------------------
                  ! (a) If FRACLOST >= 1 or FRACLOST < 0, stop the run
                  ! (b) If not, then subtract drydep losses from STT
                  !=====================================================
                  IF ( FRACLOST >= 1 .or. FRACLOST < 0 ) THEN
                     WRITE(6,*) 'DEPSAV*DTCHEM >=1 or <0 : '
                     WRITE(6,*) 'DEPSAV*DTCHEM   = ', FRACLOST
                     WRITE(6,*) 'DEPSAV(I,J,1,N) = ', DEPSAV(I,J,N)
                     WRITE(6,*) 'DTCHEM (s)      = ', DTCHEM
                     WRITE(6,*) 'I,J             = ', I, J
                     WRITE(6,*) 'STOP in dryflxRnPbBe.f'
                     CALL GEOS_CHEM_STOP
                  ENDIF

                  ! AMT_LOST = amount of species lost to drydep [kg]
                  AMT_LOST = STT(I,J,L,NN) * FRACLOST

                  ! ND44 diagnostic: drydep flux [kg/s]
                  IF ( ND44 > 0 ) THEN
                     AD44(I,J,N,1) = AD44(I,J,N,1) + ( AMT_LOST/DTCHEM )
                  ENDIF

                  ! Subtract AMT_LOST from the STT array [kg]
                  STT(I,J,L,NN) = STT(I,J,L,NN) - AMT_LOST

#elif defined( GEOS_3 ) || defined( GEOS_4 )
                  !=====================================================
                  ! GEOS-3 or GEOS-4:
                  ! -----------------
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
                     AD44(I,J,N,1) = AD44(I,J,N,1) + ( AMT_LOST/DTCHEM ) 
                  ENDIF

                  ! Subtract AMT_LOST from the STT array [kg]
                  STT(I,J,L,NN)  = STT(I,J,L,NN) - AMT_LOST
#endif
               ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE DRYFLXRnPbBe

!------------------------------------------------------------------------------

      SUBROUTINE DEPVEL(NPTS,RADIAT,TEMP,SUNCOS,F0,HSTAR,XMW,
     1                  AIROSOL,USTAR,CZ1,OBK,CFRAC,ZH,
     2                  LSNOW,DVEL,ZO)

      ! References to F90 modules (bmy, 3/8/01)
      USE ERROR_MOD, ONLY : IT_IS_NAN
                        
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

      REAL*8  RI(NTYPE),RLU(NTYPE),RAC(NTYPE),RGSS(NTYPE),
     1        RGSO(NTYPE),RCLS(NTYPE),RCLO(NTYPE),
     2        RSURFC(MAXDEP,NTYPE)       
                                   
      REAL*8  C1X(MAXDEP),VD(MAXDEP),VK(MAXDEP)              

      !----------------------------------------------------------------------- 
      ! This only applies to the GISS-CTM...comment out (bmy, 11/12/99)
      !REAL*8  ZO(NTYPE)           
      !-----------------------------------------------------------------------

      ! For GEOS-CTM, ZO is now of size MAXIJ and is passed via 
      ! the argument list, since it is a DAO met field. (bmy, 11/10/99)
      REAL*8  ZO(MAXIJ)           

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
      REAL*8  RCLX,RIXX,BIOFIT,DIFFG
      REAL*8  PRESS
      DATA PRESS /1.5D5/
C
C Logical for snow and sea ice
C

      LOGICAL LSNOW(MAXIJ)
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

               VDS = 0.002D0*USTAR(IJLOOP)
               IF (OBK(IJLOOP) .LT. 0.0D0) THEN
                  VDS = VDS*(1.D0+(-300.D0/OBK(IJLOOP))**0.6667D0)
               ENDIF
C***                               
               IF ( OBK(IJLOOP) .EQ. 0.0D0 )
     c              WRITE(6,156) OBK(IJLOOP),IJLOOP,LDT
 156           FORMAT(1X,'OBK(IJLOOP)=',E11.2,1X,' IJLOOP =',I4,
     c              1X,'LDT=',I3/) 
               CZH  = ZH(IJLOOP)/OBK(IJLOOP)
               IF (CZH.LT.-30.0D0) VDS = 0.0009D0*USTAR(IJLOOP)*
     x                             (-CZH)**0.6667D0
C*                                 
C*    Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
C*    VDSMAX is taken from Table 2 of Walcek et al. [1986].
C*    Invert to get corresponding R

               RSURFC(K,LDT) = 1.D0/MIN(VDS, DBLE(IVSMAX(II))/1.D4)
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

!-----------------------------------------------------------------------------
! This section only applies to the GISS-CTM...comment out (bmy, 11/12/99)
!
!C** Get roughness heights; they are specified constants for each surface
!C** type except over water where zo = f(u*).  The latter dependence
!C** is from equation (6) of Hicks and Liss [1976]. 
!              DO 200 IW=1,NWATER
!                  IF (IOLSON .NE. IWATER(IW)) GOTO 200
!                  ZO(LDT) = 1.4D-02*USTAR(IJLOOP)*USTAR(IJLOOP)/9.8D0   
!     1                               + 1.1D-01*XNU/USTAR(IJLOOP)        
!                  GOTO 210
! 200          CONTINUE
!              ZO(LDT) = DBLE(IZO(IOLSON))*1.D-4
! 210          CONTINUE
!-----------------------------------------------------------------------------

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
C*********************************************************************
            CKUSTR = XCKMAN*USTAR(IJLOOP)

            !-----------------------------------------------------------------
            ! This applies to the GISS-CTM...comment out (bmy, 11/12/99)
            !REYNO = USTAR(IJLOOP)*ZO(LDT)/XNU 
            !-----------------------------------------------------------------

            ! For GEOS-CTM, Z0 is now of dimension MAXIJ (bmy, 11/12/99)
            REYNO = USTAR(IJLOOP)*ZO(IJLOOP)/XNU  
            IF ( OBK(IJLOOP) .EQ. 0.0D0 )
     c           WRITE(6,211) OBK(IJLOOP),IJLOOP,LDT                
 211        FORMAT(1X,'OBK(IJLOOP)=',E11.2,1X,' IJLOOP = ',I4,1X,
     c           'LDT=',I3/) 
            CORR1 = CZ/OBK(IJLOOP)

            !-----------------------------------------------------------------
            ! This applies to the GISS-CTM...comment out (bmy, 11/12/99)
            !Z0OBK = ZO(LDT)/OBK(IJLOOP)
            !-----------------------------------------------------------------

            ! For GEOS-CTM, Z0 is now of dimension MAXIJ (bmy, 11/12/99)
            Z0OBK = ZO(IJLOOP)/OBK(IJLOOP) 
            LRGERA(IJLOOP) = .FALSE.
            IF (CORR1 .GT. 0.D0) THEN
               IF (CORR1 .GT.  1.5D0) LRGERA(IJLOOP) = .TRUE.
            ELSEIF(CORR1 .LE. 0.D0) THEN
               IF (CORR1 .LE. -2.5D0) CORR1 = -2.5D0
               CORR2 = LOG(-CORR1)
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
            IF (REYNO .LT. 1.0D0) GOTO 220

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
     &              (.74D0*LOG(CORR1/Z0OBK) + 4.7D0*(CORR1-Z0OBK))
            ELSEIF(LRGERA(IJLOOP)) THEN
C*... very stable conditions
               RA = 1.D+04
            ENDIF
C* check that RA is positive; if RA is negative (as occasionally
C* happened in version 3.1) send a warning message.
            !-----------------------------------------------------------------
            ! Debug output for GISS-CTM...comment out (bmy, 5/8/00)
            !IF (RA .LT. 0.) WRITE (6,1001) RA,CZ,ZO(LDT),OBK(IJLOOP)
            !-----------------------------------------------------------------

            ! For GEOS-CTM, We use ZO(MAXIJ), and IJLOOP is the index.
            ! Also, if RA is < 0, set RA = 0 (bmy, 11/12/99)
            IF (RA .LT. 0.D0) THEN
               WRITE (6,1001) IJLOOP,RA,CZ,ZO(IJLOOP),OBK(IJLOOP)  
               RA = 0.0D0
            ENDIF
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
!-----------------------------------------------------------------------------
! Prior to 5/8/00:
! 220          CONTINUE 
!C** ... aerodynamically smooth surface
!              DO 230 K = 1,NUMDEP 
!                 IF (.NOT.LDEP(K)) GOTO 230
!
!C...Simply specify a large C1X(K) for aerosols because C1X(K) = Ra + (Rb+Rs)
!C   where Rb is large over aerodynamically smooth surface. Here do not 
!C   calculate Ra and Rs for aerosols since Rb is large anyway. (hyl,10/15/99)
!                  IF ( AIROSOL(K) ) THEN
!                     C1X(K) = 1.D4 
!                  ELSE
!                     RA = (1.D0/CKUSTR)
!     *               *(LOG(CZ*CKUSTR/DIFFG(TEMPK,PRESS,XMW(K))) - SIH)
!                     C1X(K) = RA + RSURFC(K,LDT)
!                  ENDIF 
!
! 230          CONTINUE
!-----------------------------------------------------------------------------
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
            ENDIF

            ! Now check for IEEE NaN (not-a-number) condition 
            ! before returning to calling program (bmy, 4/16/00)
            ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
            IF ( IT_IS_NAN( DVEL(IJLOOP,K) ) ) THEN
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
            ENDIF
 550     CONTINUE
 560  CONTINUE

      ! Return to calling program
      END SUBROUTINE DEPVEL

!------------------------------------------------------------------------------

      SUBROUTINE MODIN
!
!******************************************************************************
!  Subroutine MODIN reads Olson's data from the file "drydep.table".
!  (bmy, 4/1/02, 11/21/02)
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
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! for DATA_DIR (bmy, 7/6/01)

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
!  file in the data directory (bmy, 7/6/01, 11/21/02)
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
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

      IMPLICIT NONE
      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_SETUP" ! DATA_DIR

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

      SUBROUTINE INIT_DRYDEP
!
!******************************************************************************
!  Subroutine INIT_DRYDEP initializes certain variables for the GEOS-CHEM
!  dry deposition subroutines. (bmy, 11/19/02, 7/21/03)
!
!  NOTES:
!  (1 ) Added N2O5 as a drydep tracer, w/ the same drydep velocity as
!        HNO3.  Now initialize PBLFRAC array. (rjp, bmy, 7/21/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TRACERID_MOD

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX, NTRACE
#     include "CMN_SETUP" ! LDRYD

      ! Local variables
      INTEGER :: AS, N

      !=================================================================
      ! INIT_DRYDEP begins here!
      !=================================================================

      ! Zero variables
      DRYDNO2    = 0
      DRYDPAN    = 0
      DRYDHNO3   = 0 
      NUMDEP     = 0
      NTRAIND(:) = 0
      NDVZIND(:) = 0
      HSTAR(:)   = 0d0
      F0(:)      = 0d0
      XMW(:)     = 0d0
      AIROSOL(:) = .FALSE.

      !=================================================================
      ! First identify tracers that dry deposit and then initialize 
      ! DEPNAME, NDVZIND, HSTAR, F0, XMW and AIROSOL accordingly
      !=================================================================
      DO N = 1, NTRACE

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

         ENDIF
      ENDDO

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( DEPSAV( IIPAR, JJPAR, NUMDEP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEPSAV' )
      DEPSAV = 0d0

      ALLOCATE( PBLFRAC( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBLFRAC' )
      PBLFRAC = 0d0

      !=================================================================
      ! Echo information to stdout
      !=================================================================
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'INIT_DRYDEP: List of dry deposition species:'
      WRITE( 6, '(a)'   )
     & '  #   Name  Tracer DEPVEL Henry''s    React.   Molec.  Aerosol?'
      WRITE( 6, '(a)'   )
     & '            Number Index  Law Const  Factor   Weight  (T or F)'
      WRITE( 6, '(a)'   ) REPEAT( '-', 65 )

      DO N = 1, NUMDEP
         WRITE( 6, 100 ) N, TRIM( DEPNAME(N) ), NTRAIND(N), NDVZIND(N), 
     &                   HSTAR(N),  F0(N),      XMW(N),     AIROSOL(N)
      ENDDO
 100  FORMAT( i3, 3x, a4, 2(3x,i3), 4x, es8.1, 2(3x,f6.3), 3x, L3 )

      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE INIT_DRYDEP

!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_DRYDEP

      !=================================================================
      ! CLEANUP_DRYDEP begins here!
      !=================================================================
      IF ( ALLOCATED( DEPSAV  ) ) DEALLOCATE( DEPSAV  )
      IF ( ALLOCATED( PBLFRAC ) ) DEALLOCATE( PBLFRAC )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DRYDEP

!------------------------------------------------------------------------------

      END MODULE DRYDEP_MOD

      
