!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: drydep_mod
!
! !DESCRIPTION: Module DRYDEP\_MOD contains variables and routines for the
!  GEOS-Chem dry deposition scheme.
!\\
!\\
! !INTERFACE: 
!
      MODULE DRYDEP_MOD
! 
! !USES:
!
      USE CMN_SIZE_MOD                          ! Size parameters
      USE CMN_DIAG_MOD                          ! Diag counters & flags
      USE CMN_GCTM_MOD                          ! Physical constants
      USE COMMSOIL_MOD                          ! Soil wetness variables
      USE COMODE_MOD                            ! Large arrays for SMVGEAR
      USE COMODE_LOOP_MOD                       ! Formerly "comode.h"
      USE DAO_MOD                               ! Met field subroutines
      USE DIAG_MOD,       ONLY : AD44           ! Diagnostic arrays
      USE DIRECTORY_MOD                         ! Data directory paths
      USE ERROR_MOD                             ! Error handling routines
      USE GET_POPSINFO_MOD                      ! For POPs simulation
      USE GRID_MOD,       ONLY : GET_AREA_CM2   ! Grid box surface areas [cm2]
      USE LOGICAL_MOD                           ! GEOS-Chem logical switches
      USE PBL_MIX_MOD                           ! Boundary layer quantities
      USE PRESSURE_MOD,   ONLY : GET_PEDGE      ! Pressure @ level edges
      USE TIME_MOD,       ONLY : GET_TS_CHEM    ! Chemistry timestep
      USE TRACERID_MOD                          ! Tracer ID flags
#if   defined( TOMAS )
      USE TOMAS_MOD                             ! For TOMAS microphysics
#endif

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_DRYDEP     
      PUBLIC :: DO_DRYDEP
      PUBLIC :: DRYFLX   
      PUBLIC :: DRYFLXH2HD
      PUBLIC :: DRYFLXRnPbBe       
      PUBLIC :: DVZ_MINVAL         
      PUBLIC :: INIT_DRYDEP
      PUBLIC :: INIT_WEIGHTSS
#if   defined( APM )
      PUBLIC :: DEPSAVSIZE
#endif
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: DEPNAME
      PUBLIC :: DEPSAV
      PUBLIC :: SHIPO3DEP
      PUBLIC :: MAXDEP
      PUBLIC :: NUMDEP
      PUBLIC :: NTRAIND
      PUBLIC :: DRYHg0, DRYHg2, DryHgP    !CDH
      PUBLIC :: DRYPOPG, DRYPOPP_OC, DRYPOPP_BC
      PUBLIC :: IDEP,   IRGSS,  IRAC, IRCLS
      PUBLIC :: IRGSO,  IRLU,   IRI,  IRCLO, DRYCOEFF
      PUBLIC :: NDVZIND ! MSL -> For MPI broadcasting in GIGC
!
! !REMARKS:
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
!  (16) Karl, T., Harley, P., Emmons, L., Thornton, B., Guenther, A., Basu, C.,
!        Turnipseed, A., and Jardine, K.: Efficient Atmospheric Cleansing of
!        Oxidized Organic Trace Gases by Vegetation, Science, 330, 816-819,
!        10.1126/science.1192534, 2010. 
!
! !REVISION HISTORY:
!  27 Jan 2003 - R. Yantosca - Moved standalone routines into this module
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
!  (28) Change MAXDEP from 50 to 81 (win, 7/14/09)
!  (28a) modified to use Zhang 2001 for all non-size resolved aerosols (hotp)
!  (29) Add aromatics SOA (dkh)
!  (30) Add new species. Some tracers give 2 deposition species: ISOPN-> ISOPNB
!       and ISOPND. (fp)
!  (31) Updates for mercury simulation (ccc, 5/17/10)
!  (32) Add POPs (eck, 9/20/10)
!  (33) Increase MAXDEP to 51 for dicarbonyls simulation. (ccc, 10/8/10)
!  15 Feb 2011 - R. Yantosca - Added modifications for APM from G. Luo
!  01 Aug 2011 - J. Fisher - Set aerosol dry deposition velocity to 0.03 cm/s
!                            over snow and ice based on Nilsson & Rannik, 2001
!  21 Dec 2011 - M. Payer  - Updates for sea salt (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer  - Added ProTeX headers
!  10 Jan 2012 - M. Payer  - Update to use local surface pressure
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  26 Mar 2012 - R. Yantosca - Now reference CMN_SIZE_MOD at the top of module
!  26 Mar 2012 - R. Yantosca - Replace NNTYPE, NNPOLY, NNVEGTYPE w/ the 
!                              values NTYPE, NPOLY, NVEGTYPE from CMN_SIZE
!  26 Mar 2012 - R. Yantosca - Now retire MODIN and RDDRYCF; read drydep inputs
!                              from a netCDF file w/ routine READ_DRYDEP_INPUTS
!  26 Mar 2012 - R. Yantosca - Reorganize module USE statements for clarity
!  09 Apr 2012 - R. Yantosca - Now replace IJREG, IJLAND, IJUSE, XYLAI arrays
!                              with IREG, ILAND, IUSE, XLAI.
!  31 Jul 2012 - R. Yantosca - Modifications for grid-independence
!  11 Dec 2012 - M. Long     - Now call READ_DRYDEP_INPUTS from INIT_DRYDEP
!  11 Dec 2012 - R. Yantosca - Now call INIT_WEIGHTSS from INIT_DRYDEP
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  26 Feb 2013 - R. Yantosca - Now use Input_Opt fields where possible
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  29 Jan 2014 - R. Yantosca - Now set MAXDEP=105 for all simulations.  For
!                              TOMAS we had MAXDEP=100; this is close enough.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SOAupdate: 55 to 105 jje (8/19/2010)
      INTEGER, PARAMETER :: MAXDEP    = 105     ! Max # of depositing species
      REAL*8,  PARAMETER :: XCKMAN    = 0.4d0   ! Von Karman's constant
      INTEGER, PARAMETER :: NR_MAX    = 200     ! # of seasalt bins 
      INTEGER, PARAMETER :: NDRYDTYPE = 11      ! # of drydep land types
!
! PRIVATE TYPES:
!
      !========================================================================
      !  MODULE VARIABLES:
      !
      !  MAXDEP    : Maximum number of drydep species
      !  XCKMAN    : Von Karman constant?
      !  DRYDHNO3  : Internal flag for location of HNO3 in DEPVEL
      !  DRYDNO2   : Internal flag for location of NO2  in DEPVEL
      !  DRYDPAN   : Internal flag for location of PAN  in DEPVEL
      !  NUMDEP    : Actual number of drydep species
      !  NWATER    : Number of Olson's surface types that are water
      !  AIROSOL   : Array flags to denote aerosol drydep species
      !  IDEP      : ID #'s for dry deposition surface types 
      !  IRAC      : ???       resistance for drydep land type
      !  IRCLO     : ???       resistance for drydep land type
      !  IRCLS     : ???       resistance for drydep land type
      !  IRGSO     : ???       resistance for drydep land type
      !  IRGSS     : ???       resistance for drydep land type
      !  IRI       : Internal  resistance for drydep land types
      !  IRLU      : Cuticular resistance for drydep land types
      !  IVSMAX    : ???       resistance for drydep land type
      !  IWATER    : ID #'s for Olson surface types that are water 
      !  IZO       : Roughness heights for each Olson surface type
      !  NDVZIND   : Index array for ordering drydep species in DEPVEL
      !  NTRAIND   : Stores tracer numbers of drydep species
      !  DEPSAV    : Array containing dry deposition frequencies [s-1]
      !  PBLFRAC   : Array for multiplicative factor for drydep freq
      !  DRYCOEFF  : Polynomial coefficients for dry deposition
      !  HSTAR     : Henry's law constant
      !  F0        : Reactivity factor for biological oxidation
      !  XMW       : Molecular weight of drydep species [kg]
      !  A_RADI    : Radius of aerosol for size-resolved drydep [um]
      !  A_DEN     : Density of aerosol for size-res'd drydep [kg/m3]
      !  DEPNAME   : Names of dry deposition species 
      !
      !  NOTE: these variables are defined in CMN_SIZE_mod.F
      !    NTYPE    : Max # of landtypes / grid box
      !    NPOLY    : Number of drydep polynomial coefficients
      !    NVEGTYPE : Number of Olson land types
      !========================================================================

      ! Scalars
      INTEGER              :: NUMDEP,   NWATER
      INTEGER              :: DRYDHNO3, DRYDNO2, DRYDPAN, DRYDH2O2
      INTEGER              :: DRYHg0,   DRYHg2,  DryHgP
      ! Additional variables for POPs sim (eck, 9/20/10)
      INTEGER              :: DRYPOPP_OC
      INTEGER              :: DRYPOPP_BC
      INTEGER              :: DRYPOPG

      ! Arrays for Baldocchi drydep polynomial coefficients
      REAL*8               :: DRYCOEFF(NPOLY    )

      ! Arrays that hold information for each of the 74 Olson land types
      INTEGER              :: INDOLSON(NVEGTYPE )
      INTEGER              :: IDEP    (NVEGTYPE )
      INTEGER              :: IZO     (NVEGTYPE )
      INTEGER              :: IWATER  (NVEGTYPE )

      ! Arrays that hold information for each of the 11 drydep land types
      INTEGER              :: IDRYDEP (NDRYDTYPE)
      INTEGER              :: IRAC    (NDRYDTYPE)
      INTEGER              :: IRCLO   (NDRYDTYPE)
      INTEGER              :: IRCLS   (NDRYDTYPE)
      INTEGER              :: IRGSS   (NDRYDTYPE)
      INTEGER              :: IRGSO   (NDRYDTYPE)
      INTEGER              :: IRI     (NDRYDTYPE)
      INTEGER              :: IRLU    (NDRYDTYPE)
      INTEGER              :: IVSMAX  (NDRYDTYPE)

      ! Arrays that hold information about the dry-depositing species
      LOGICAL              :: AIROSOL (MAXDEP   )
      INTEGER              :: NDVZIND (MAXDEP   )
      INTEGER              :: NTRAIND (MAXDEP   )
      REAL*8               :: HSTAR   (MAXDEP   )
      REAL*8               :: KOA     (MAXDEP   )
      REAL*8               :: F0      (MAXDEP   )
      REAL*8               :: XMW     (MAXDEP   )
      REAL*8               :: A_RADI  (MAXDEP   )
      REAL*8               :: A_DEN   (MAXDEP   )
      CHARACTER(LEN=14)    :: DEPNAME (MAXDEP   )

      ! Allocatable arrays
      REAL*8,  ALLOCATABLE :: DEPSAV  (:,:,:    )   ! Drydep frequencies
      REAL*8,  ALLOCATABLE :: DMID    (:        )   ! 
      REAL*8,  ALLOCATABLE :: SALT_V  (:        )
      REAL*8,  ALLOCATABLE :: SHIPO3DEP(:,:)
#if   defined( APM )
      REAL*8,  ALLOCATABLE :: DEPSAVSIZE(:,:,:)
#endif

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
! !IROUTINE: do_drydep
!
! !DESCRIPTION: Subroutine DO\_DRYDEP is the driver for the GEOS-CHEM dry
!  deposition scheme. DO\_DRYDEP calls DEPVEL to compute deposition velocities
!  [m/s], which are then converted to [cm/s].  Drydep frequencies are also
!  computed. (lwh, gmg, djj, 1989, 1994; bmy, 2/11/03, 5/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_DRYDEP( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
#if   defined( APM ) 
      USE APM_INIT_MOD,       ONLY : NSO4, NSEA, NDSTB
#endif
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
!  NOTE: Modeled aerosol dry deposition velocities over snow and ice surfaces 
!  in the Arctic are much higher than estimated from measured values (e.g., 
!  Ibrahim et al. [1983]; Duan et al. [1988]; Nilsson and Rannik [2001]).  
!  We will impose a dry deposition velocity of 0.03 cm/s for all aerosols 
!  over snow and ice surfaces. (Jenny Fisher, 01 Aug 2011)
!
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989 
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
! !REVISION HISTORY: 
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
!  15 Feb 2011 - R. Yantosca - Added modifications for APM from G. Luo
!  01 Aug 2011 - J. Fisher - Set aerosol dry deposition velocity to 0.03 cm/s
!                            over snow and ice based on Nilsson & Rannik, 2001
!  15 Aug 2011 - R. Yantosca - Now reference IDTxxx flags from tracerid_mod.f
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  10 Jan 2012 - M. Payer    - Added local surface pressure
!  26 Mar 2012 - R. Yantosca - Now read drydep inputs from a netCDF file
!                              via routine READ_DRYDEP_INPUTS
!  26 Mar 2012 - R. Yantosca - Remove calls to obsolete MODIN, RDDRYCF routines
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Now make SUNCOS_MID a local array of size
!                              MAXIJ, populated from State_Met%SUNCOSmid
!  11 Dec 2012 - R. Yantosca - Now do not call READ_DRYDEP_INPUTS and
!                              INIT_WEIGHTSS when using the ESMF environment
!  11 Dec 2012 - R. Yantosca - Remove FIRST variable, as we now read inputs
!                              from disk in routine INIT_DRYDEP
!  12 Dec 2012 - R. Yantosca - Now pass State_Met to DEPVEL
!  26 Feb 2013 - R. Yantosca - Now use Input_Opt fields where possible.  This
!                              facilitates connection to the GEOS-5 GCM.   
!  31 May 2013 - R. Yantosca - Now pass Input_Opt & State_Chm to DEPVEL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER       :: I,   J,   L,  N, IJLOOP, NN, NDVZ
      REAL*8        :: DVZ, THIK

      ! Arrays
      LOGICAL       :: LSNOW (MAXIJ       )   ! Flag for snow/ice on the sfc.
      REAL*8        :: CZ1   (MAXIJ       )   ! Midpt ht of 1st model level [m]
      REAL*8        :: TC0   (MAXIJ       )   ! Grid box sfc temperature [K]
      REAL*8        :: ZH    (MAXIJ       )   ! PBL height [m]
      REAL*8        :: OBK   (MAXIJ       )   ! Monin-Obhukov Length [m]
      REAL*8        :: CFRAC (MAXIJ       )   ! Column cloud frac [unitless]
      REAL*8        :: RADIAT(MAXIJ       )   ! Solar radiation [W/m2]
      REAL*8        :: USTAR (MAXIJ       )   ! Grid box friction velocity[m/s]
      REAL*8        :: RHB   (MAXIJ       )   ! Relative humidity [unitless]
      REAL*8        :: DVEL  (MAXIJ,MAXDEP)   ! Drydep velocities [m/s]
      REAL*8        :: PRESSU(MAXIJ       )   ! Local surface pressure [Pa]
      REAL*8        :: W10   (MAXIJ       )   ! 10m windspeed [m/s]

      ! Now make SUNCOS_MID a local array (bmy, 11/28/12)
      REAL*8        :: SUNCOS_MID(MAXIJ   )   ! COS(SZA) @ midpoint of the
                                              !  current chemistry timestep

      ! Dimension AZO for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
      REAL*8        :: AZO   (NTYPE       )   ! Z0, per Olson land type
#else
      REAL*8        :: AZO   (MAXIJ       )   ! Z0, per (I,J) square
#endif

      ! For ESMF, need to assign these from Input_Opt
      LOGICAL       :: LNLPBL
      LOGICAL       :: LPRT

#if   defined( APM )
      ! size resolved dry velocities for APM simulation
      REAL*8            :: DVELSIZE(MAXIJ,NSO4+NSEA+NDSTB)
#endif

      !=================================================================
      ! DO_DRYDEP begins here!
      !=================================================================

      ! Assume success
      RC     = GIGC_SUCCESS
      
      ! Copy values from the Input Options object to local variables
      LNLPBL  = Input_Opt%LNLPBL
      LPRT    = Input_Opt%LPRT

      ! Be careful only to copy NUMDEP values from Input_Opt%NDVZIND
      ! to the NDVZIND local array.  This will prevent an OOB error.
      DO N = 1, Input_Opt%NUMDEP
         NDVZIND(N) = Input_Opt%NDVZIND(N)
      ENDDO

      ! Call METERO to obtain meterological fields (all 1-D arrays)
      ! Added sfc pressure as PRESSU and 10m windspeed as W10 
      !  (jaegle 5/11/11, mpayer 1/10/12)
      CALL METERO( State_Met, CZ1,    TC0,     OBK,       CFRAC,  
     &             RADIAT,    AZO,    USTAR,   ZH,        LSNOW, 
     &             RHB,       PRESSU, W10,     SUNCOS_MID             )

      ! Call DEPVEL to compute dry deposition velocities [m/s]
      CALL DEPVEL( Input_Opt, State_Met, State_Chm,   MAXIJ,  
     &             RADIAT,    TC0,       SUNCOS_MID,  F0,        
     &             HSTAR,     XMW,       AIROSOL,     USTAR,      
     &             CZ1,       OBK,       CFRAC,       ZH,        
     &             LSNOW,     DVEL,      AZO,         RHB,     
     &             PRESSU,    W10                                     )

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
         THIK    = State_Met%BXHEIGHT(I,J,1)
         
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
            
            ! For ACET, we need to only do drydep over the land
            ! and not over the oceans.
            IF ( NN == IDTACET ) THEN

               IF ( IS_LAND( I, J, State_Met ) ) THEN
                  DVZ = 0.1d0
               ELSE
                  DVZ = 0d0
               ENDIF

            ENDIF

            ! Set Vd for aerosol tracers over snow/ice (J. Fisher, 8/1/11)
            IF ( ( AIROSOL(N)    ) .and. ( LSNOW(IJLOOP) )  .and. 
     &           ( NN /= IDTSALA ) .and. ( NN /= IDTSALC )  .and.
     &           ( NN /= IDTSO4s ) .and. ( NN /= IDTNITs )  .and.
     &           ( NN /= IDTDST1 ) .and. ( NN /= IDTDST2 )  .and.
     &           ( NN /= IDTDST3 ) .and. ( NN /= IDTDST4 ) ) THEN
               DVZ = 3.0d-2
            ENDIF

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

! Activate this shortly (bmy, 3/5/13)
!      IF ( Input_Opt%LDRYD ) THEN
!         DO N = 1, Input_Opt%NUMDEP
!            State_Chm%DEPSAV(:,:,N) = DEPSAV(:,:,N)
!         ENDDO
!      ENDIF

      !### Debug
      IF ( LPRT .and. am_I_Root ) THEN
         CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )
      ENDIF

!APM_GanLuo+
#if   defined( APM )

      !=================================================================
      ! Call DEPVELBIN to compute size resolved dry velocities
      ! for APM simulation [m/s]
      !=================================================================
      !GanLuo for newbins
      CALL DEPVELBIN( MAXIJ, RADIAT,   TC0,        SUNCOS_MID,
     &                USTAR, CZ1,      OBK,        CFRAC,  
     &                ZH,    LSNOW,    DVELSIZE,   AZO,
     &                RHB,   State_Met )

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
         THIK    = State_Met%BXHEIGHT(I,J,1)

         ! Now we calculate drydep throughout the entire PBL.  
         ! Make sure that the PBL depth is greater than or equal 
         ! to the thickness of the 1st layer (rjp, bmy, 7/21/03)
         THIK    = MAX( ZH(IJLOOP), THIK )

         ! Loop over drydep species
         DO N = 1, NSO4+NSEA+NDSTB

            ! Dry deposition velocity [cm/s]
            DVZ  = DVELSIZE(IJLOOP,N) * 100.d0

            ! Set minimum velocity for sulfate tracers
            DVZ  = MAX( 1.0d-2, DVZ )

            ! Dry deposition frequency [1/s]
            DEPSAVSIZE(I,J,N) = ( DVZ / 100.d0 ) / THIK

         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_DRYDEP: after DEPVELBIN' )
#endif

      END SUBROUTINE DO_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dvz_minval
!
! !DESCRIPTION: Function DVZ\_MINVAL sets minimum values for drydep velocities 
!  for SULFATE TRACERS, according to Mian Chin's GOCART model. (rjp, bmy, 
!  11/21/02, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      FUNCTION DVZ_MINVAL( N, LSNOW, DVZ ) RESULT( NEWDVZ )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: N      ! Tracer number
      LOGICAL, INTENT(IN) :: LSNOW  ! Flag for denoting snow/ice
      REAL*8,  INTENT(IN) :: DVZ    ! Deposition velocity [cm/s]
!
! !RETURN VALUE:
!
      REAL*8              :: NEWDVZ
!
! !REVISION HISTORY: 
!  (1 ) Don't put a min drydep value on H2O2 for offline run (rjp, bmy,3/31/03)
!  (2 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
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

      END FUNCTION DVZ_MINVAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: metero
!
! !DESCRIPTION: Subroutine METERO calculates meteorological constants needed
!  for the dry deposition velocity module. (lwh, gmg, djj, 1989, 1994; bmy,
!  10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE METERO( State_Met, CZ1,    TC0,  OBK,       CFRAC,  
     &                   RADIAT,    AZO,    USTR, ZH,        LSNOW, 
     &                   RHB,       PRESSU, W10,  SUNCOS_MID        )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      LOGICAL, INTENT(OUT) :: LSNOW (MAXIJ)  ! Flag for denoting snow/ice
      REAL*8,  INTENT(OUT) :: CZ1   (MAXIJ)  ! Midpt ht of 1st model level [m]
      REAL*8,  INTENT(OUT) :: TC0   (MAXIJ)  ! Grid box sfc temperature [K]
      REAL*8,  INTENT(OUT) :: OBK   (MAXIJ)  ! Monin-Obhukov length [m]
      REAL*8,  INTENT(OUT) :: CFRAC (MAXIJ)  ! Column cloud fraction [unitless]
      REAL*8,  INTENT(OUT) :: RADIAT(MAXIJ)  ! Solar radiation @ ground [W/m2]
      REAL*8,  INTENT(OUT) :: RHB   (MAXIJ)  ! Rel humidity at sfc [unitless]
      REAL*8,  INTENT(OUT) :: USTR  (MAXIJ)  ! Friction velocity [m/s]
      REAL*8,  INTENT(OUT) :: ZH    (MAXIJ)  ! PBL height [m]
      REAL*8,  INTENT(OUT) :: PRESSU(MAXIJ)  ! Local surface pressure [Pa]
      REAL*8,  INTENT(OUT) :: W10   (MAXIJ)  ! 10 meter windspeed [m/s]
      REAL*8,  INTENT(OUT) :: SUNCOS_MID(MAXIJ)  ! COS(SZA) @ midpt of current
                                                 !  chemistry timestep

      ! Dimension AZO for GCAP or GEOS met fields (swu, bmy, 5/25/05)
#if   defined( GCAP )
      REAL*8, INTENT(OUT)  :: AZO(NTYPE)     ! Roughness heights, by landtype
#else
      REAL*8, INTENT(OUT)  :: AZO(MAXIJ)     ! Roughness heights, by grid box
#endif
!
! !REMARKS:
!  NOTE: We save into arrays of dimension MAXIJ=IIPAR*JJPAR for compatibility
!  with the legacy drydep routine DEPVEL.
!                                                                             .
!  References (see full citations above):
!  ============================================================================
!  (1 ) Wesely, M. L., 1989. 
!  (2 ) Jacob, D.J., and S.C. Wofsy, 1990
!
! !REVISION HISTORY: 
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
!  (9 ) Add sea level pressure and 10m windspeed as arguments (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  10 Jan 2012 - M. Payer    - Added local surface pressure
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Add SUNCOS_MID to the argument list and 
!                              populate that with State_Met%SUNCOSmid
!  21 Oct 2013 - R. Yantosca - Bug fix: need to hold SP private in OMP loop
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I,  J,  IJLOOP
      REAL*8                :: THIK
      REAL*8                :: SP

      ! External functions
      REAL*8, EXTERNAL      :: SFCWINDSQR ! Surface wind speed (jaegle 5/11/11)

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
!$OMP+PRIVATE( I, J, IJLOOP, THIK, SP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! 1-D grid box index
         IJLOOP         = ( (J-1) * IIPAR ) + I

         ! THIK = thickness of layer 1 [m]
         THIK           = State_Met%BXHEIGHT(I,J,1)

         ! Midpoint height of first model level [m]
         CZ1(IJLOOP)    = THIK / 2.0d0

         ! Local surface pressure [hPa] (mpayer, 1/10/12)
         SP             = GET_PEDGE(I,J,1)

         ! Convert from hPa to Pa for SFCPRESS
         PRESSU(IJLOOP) = SP * 1.d2

         !==============================================================
         ! Return meterological quantities as 1-D arrays for DEPVEL
         !==============================================================

#if   !defined( GCAP )
         ! For GEOS: Roughness height [m] is a function of lat/lon
         AZO(IJLOOP)    = State_Met%Z0(I,J)
#endif

         ! Column cloud fraction [unitless]
         CFRAC(IJLOOP)  = State_Met%CLDFRC(I,J)

         ! Set logical LSNOW if snow and sea ice (ALBEDO > 0.4)
         LSNOW(IJLOOP)  = ( State_Met%ALBD(I,J) > 0.4 )

         ! Monin-Obhukov length [m]
         OBK(IJLOOP)    = GET_OBK( I, J, State_Met )

         ! Solar insolation @ ground [W/m2]
         RADIAT(IJLOOP) = State_Met%RADSWG(I,J) 

         ! Surface temperature [K]
         TC0(IJLOOP)    = State_Met%TS(I,J)

         ! Friction velocity [m/s]
         USTR(IJLOOP)   = State_Met%USTAR(I,J)

         ! Mixed layer depth [m]
         ZH(IJLOOP)     = GET_PBL_TOP_m( I, J )

         ! Relative humidity @ surface [unitless] (bec, bmy, 4/13/05)
         !RHB(IJLOOP)    = MIN( 0.99d0, RH(I,J,1) * 1.d-2 ) 
         !  changed to 98% due to vapor pressure lowering above sea water
         !  (Lewis & Schwartz, 2004)
         !  jaegle (5/11/11)
         RHB(IJLOOP)    = MIN( 0.98d0, State_Met%RH(I,J,1) * 1.d-2 )

         ! 10m windspeed (jaegle 5/11/11)
         W10(IJLOOP)    =  SQRT( SFCWINDSQR( I, J,
     &                                       State_Met%U10M,
     &                                       State_Met%V10M ) )

         ! Cosine of solar zenith angle at midpoint 
         ! of the current chemistry timestep.
         SUNCOS_MID(IJLOOP) = State_Met%SUNCOSmid(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      END SUBROUTINE METERO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dryflx
!
! !DESCRIPTION: Subroutine DRYFLX sets up the dry deposition flux diagnostic
!  for tracers which are part of the SMVGEAR mechanism. (bmy, bdf, 4/20/99,
!  3/24/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRYFLX( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE COMMSOIL_MOD
      USE COMODE_MOD
      USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD,       ONLY : IDTHNO3
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REMARKS:
!  DRYFLX is not called if non-local PBL mixing is done (i.e. when LNLPBL=F).
!  If using non-local PBL mixing, the drydep fluxes are archived in routine
!  VDIFF (GeosCore/vdiff_mod.F90).
!
! !REVISION HISTORY: 
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
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  31 Jul 2012 - R. Yantosca - Now loop from 1..LLPAR for GIGC compatibility
!  05 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%LSOILNOX and Input_Opt%NUMDEP
!  30 Oct 2013 - M. Sulprizio- Bug fix: Hold AREA_CM2 private in !$OMP loop
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, JJ, JLOOP, L, L_PBLTOP, N, NK, NN
      REAL*8  :: DTCHEM, PBL_MAX, TDRYFX, AREA_CM2

      !=================================================================
      ! DRYFLX begins here!
      !=================================================================

      ! Assume success
      RC      = GIGC_SUCCESS

      ! Return unless we have turned on ND44 drydep diagnostic
      !IF ( ND44 == 0 ) RETURN

      ! There is only drydep in the surface layer, which
      ! is accounted for in the "URBAN" chemistry slot
      NCS     = NCSURBAN

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0

      !=================================================================
      ! ND44 diagnostic: Dry deposition flux [molec/cm2/s]
      !
      ! NOTE: DRYFLX will only archive the dry deposition fluxes for
      ! tracers which are SMVGEAR species.  Fluxes for sulfate tracers
      ! will be updated in "sulfate_mod.f". (bmy, 11/19/02)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NN, NK, JJ, JLOOP, TDRYFX, AREA_CM2 )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, Input_Opt%NUMDEP

         ! Index for drydep species #N, from SMVGEAR
         NK = NTDEP(N)

         ! Dry deposition tracer id
         NN =  NTRAIND(N) 

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
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         
            ! Only deal w/ boxes w/in the boundary layer
            IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

               ! Grid box surface area [cm2]
               AREA_CM2      = GET_AREA_CM2( I, J, L )

               ! 1-D grid box index for CSPEC & VOLUME
               JLOOP         = JLOP(I,J,L)

               ! Dry dep flux [molec] for species N = 
               !  CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
               !  [molec/cm3]     * [cm3]  
               TDRYFX        = CSPEC(JLOOP,JJ) * VOLUME(JLOOP)
                    
               ! Convert TDRYFX from [molec] to [molec/cm2/s]        
               TDRYFX        = TDRYFX / ( AREA_CM2 * DTCHEM ) 
                    
!$OMP CRITICAL
               ! Save into AD44 diagnostic array
               IF ( ND44 > 0 ) THEN 
                  AD44(I,J,N,1) = AD44(I,J,N,1) + TDRYFX
               ENDIF 
!$OMP END CRITICAL

               ! Archive dry deposited N species for use in soil_nox_mod
               IF ( Input_Opt%LSOILNOX ) THEN
                  CALL SOIL_DRYDEP( I, J, L, NN, TDRYFX )
               ENDIF

            ENDIF
         ENDDO           
         ENDDO               
         ENDDO  
      ENDDO                     
!$OMP END PARALLEL DO

      END SUBROUTINE DRYFLX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dryflxRnPbBe
!
! !DESCRIPTION: Subroutine DRYFLXRnPbBe removes dry deposition losses from the
!  STT tracer array and archives deposition fluxes to the ND44 diagnostic. 
!  (hyl, bmy, bdf, 4/2/99, 5/25/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRYFLXRnPbBe( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
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
!  (1 ) Now eliminate DEPFLUX from CMN_SAV, in order to save memory.
!        DEPFLUX is now a local variable (bdf, 4/2/99)
!  (2 ) Now make DEPFLUX of dimension (IIPAR,JJPAR,MAXDEP) (bmy, 4/2/99)
!  (3 ) Now use an allocatable array for the ND44 diagnostic.
!        Also made cosmetic changes, updated comments. (bmy, 3/16/00)
!  (4 ) Eliminate obsolete code and ND63 diagnostic (bmy, 4/12/00)
!  (5 ) Added to module "RnPbBe_mod.f".  Also made cosmetic changes
!        and updated comments (bmy, 6/14/01)
!  (6 ) Updated comments (bmy, 3/29/02)
!  (7 ) Replace all instances of IM, JM, IMX, JMX, with IIPAR, JJPAR, IIPAR,
!        and JJPAR.  Now replaced DEPFLUX array w/ AMT_LOST scalar
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
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, PBL_MAX, N, NN
      REAL*8              :: DTCHEM, FRACLOST, F_UNDER_TOP, AMT_LOST

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! DRYFLXRnPbBe begins here!!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Return if drydep is turned off
      IF ( .not. Input_Opt%LDRYD ) RETURN

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Chemistry timestep in seconds
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Maximum extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L() 

      ! Add option for non-local PBL mixing scheme: only done at the surface
      ! (Lin, 03/31/09) 
      IF ( Input_Opt%LNLPBL ) PBL_MAX = 1

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

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE DRYFLXRnPbBe
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dryflxh2hd
!
! !DESCRIPTION: Subroutine DRYFLXH2HD removes dry deposition losses from the
!  tracer array and archives deposition fluxes AND VELOCITY to the ND44
!  diagnostic. (adapted from DRYFLX v5-05, jaegle 11/02/2005). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DRYFLXH2HD( State_Met, State_Chm )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD44
      USE ERROR_MOD,          ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,        ONLY : LDRYD
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_m
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches & arrays
      USE COMMSOIL_MOD             ! Soil wetness variables
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% NOTE: THIS ROUTINE SHOULD BE REWRITTEN FOR THE NEW SOIL NOX SCHEME. %%%
!  %%% We have commented out the obsolete INDEXSOIL and NLAND variables,   %%%
!  %%% which were different depending on which horizontal grid was used.   %%%
!  %%% The grid-independent GEOS-Chem cannot know a-priori what the size   %%%
!  %%% of the horizontal grid is, since it will be told that by the        %%%
!  %%% interface to the external GEOS-Chem.  Therefore, we cannot use      %%%
!  %%% fixed parameters to define the horizontal (and vertical) grids that %%%
!  %%% are used in GEOS-Chem.  (bmy, 10/30/12)                             %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY: 
!  (1) Now deposit through the PBL. Commented but kept code related to soil
!       temperature (phs, 5/16/07)
!  (2) Move XLTMMP to module MEGANUT_MOD (ccc, 11/20/09)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  09 Apr 2012 - R. Yantosca - Replace IJLAND, IJUSE arrays w/ ILAND, IUSE
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER     :: I, J, L, N, NN, M, PBL_MAX
      INTEGER     :: IJLOOP, I0, J0, IREF, JREF, K, STYP
      INTEGER     :: JLOP(IIPAR,JJPAR,1), NTYP(IIPAR,JJPAR)
      REAL*8      :: DTCHEM, FRACLOST, AMT_LOST
      REAL*8      :: THIK, DRYF, SVEL, FSOIL, AREA_CM2
      REAL*8      :: SOIL_H2, SOIL_HD, TMMP, STEMP(IIPAR,JJPAR)
      REAL*8      :: MLD
      REAL*8      :: F_UNDER_TOP

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      ! External functions, for calculating soil temperature
      REAL*8, EXTERNAL   :: SOILTEMP

      !=================================================================
      ! DRYFLXH2HD begins here!!
      !=================================================================

      ! Chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Obsolete - Call soiltype to determine whether soil is dry or
      ! wet for all land grid-boxes
      !CALL SOILTYPE

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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

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

!           ! Loop over each land grid-box
!            DO M      = 1, NLAND
!
!               !IREF   = INDEXSOIL(1,M)
!               !JREF   = INDEXSOIL(2,M)
!               I      = IREF - I0
!               J      = JREF - J0
!               IJLOOP = JLOP(I,J,1)
!
!               ! Fraction of grid box that is ocean
!               FSOIL = State_Met%FRCLND(I,J)
! 
!
!               ! Only apply dry deposition over land surfaces which
!                    ! are not covered with ice or desert (albedo < 0.4)
!               ! and if we are in the window simulation.
!               IF ( (I.GE.1)       .AND. (I.LE.IIPAR)  .AND.
!     &              (J.GE.1)       .AND. (J.LE.JJPAR)  .AND.
!     &              (FSOIL > 0.d0) .AND. 
!     &              (State_Met%ALBD(I,J) < 0.4d0) ) THEN
!
!                  ! Grid box area in cm2
!                  AREA_CM2 = GET_AREA_CM2( I, J, 1 )
!
!                  ! SVEL [cm/s] is the air-to-soil transfer velocity
!                  ! Use uniform value of 3.94d-2 cm/s over land
!                  ! not covered by snow or desert.
!                  SVEL = 3.94d-2
!		  
!                  ! if soil temperature is below freezing reduce dep vel
!                  ! by 1/2, and additional 1/2 below -15C(hup, 6/21/2005)
!                  ! for now use surface temperature (TS) instead of
!                  ! air temperature (T) jaegle, 12/12/2005
!		  IF (State_Met%TS(I,J) <= 273.15d0) SVEL = SVEL / 2.0d0
!		  IF (State_Met%TS(I,J) <= 258.15d0) SVEL = SVEL / 2.0d0
!                  
!		  ! if desert, set deposition velocity to zero by multiplying
!		  ! dep vel by the fraction covered by desert(hup, 5/1/2006)
!		  ! IJLAND+1 is the Olson Land type index 
!            	  ! 51: desert	52: desert set SVEL = 0
!		  ! IJUSE is the fraction of the grid square occupied by surface K
!		  ! in units of per mil (IJUSE=500 -> 50% of the grid square). 
!		  DO K = 1, IREG(IREF,JREF)
!	          
!                     NTYP(I,J) = ILAND(I,J,K) + 1
!                     
!                     IF (NTYP(I,J) .eq. 52 .or. NTYP(I,J) .eq. 51) THEN 
!		  	SVEL = SVEL*(1-(IUSE(I,J,K)/1.d3)) 
!                     ENDIF
!                     
!		  ENDDO
!		  
!		  ! For HD add soil fractionation with 
!                  ! an alpha coefficient of 0.943 Gerst & Quay, 2001
!		  IF (N .eq. 2) SVEL = SVEL * 0.943
!
!
!                  ! Get THIK  (cannot use ZH variable, since
!                  ! DO_DRYDEP, METERO, and DEPVEL are not called in H2/HD sims)
!
!                  ! Fraction of box (I,J,L) under PBL top [unitless]
!                  F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )
!      
!                  ! Mixed layer depth [m]
!                  MLD = GET_PBL_TOP_m( I, J )
!      
!                  ! THIK = thickness of surface layer [m]
!                  THIK    = State_Met%BXHEIGHT(I,J,1)
!
!                  ! Add option for non-local PBL mixing scheme: 
!                  ! only the surface (Lin, 03/31/09) 
!                  IF (.NOT. LNLPBL) THIK    = MAX( MLD, THIK ) 
!                  
!                  ! Dry deposition frequency [1/s]
!                  DRYF    = ( SVEL / 100.d0 ) / THIK
!		  
!		  ! FRACLOST = Fraction of species lost to drydep [unitless]
!                  FRACLOST = DRYF * DTCHEM * F_UNDER_TOP
!                  ! Add option for non-local PBL mixing scheme:
!                  ! only the surface (Lin, 03/31/29)
!                  IF (LNLPBL) FRACLOST = DRYF * DTCHEM
!
!                  !========================================================
!                  ! Proceed as follows:
!                  ! -------
!                  ! (a) If FRACLOST < 0, then stop the run.
!                  !
!                  ! (b) If FRACLOST > 1, use an exponential loss to 
!                  !     avoid negative tracer
!                  !
!                  ! (c) If FRACLOST is in the range (0-1), then use the 
!                  !     regular formula (STT * FRACLOST) to compute
!                  !     the loss from dry deposition.
!                  !========================================================
!
!                  ! Stop the run on negative FRACLOST!
!                  IF ( FRACLOST < 0 ) THEN
!                     CALL ERROR_STOP( 'FRACLOST < 0', 'dryflxH2HD' )
!                  ENDIF
!
!
!                  ! AMT_LOST = amount of tracer lost to drydep [kg]
!                  IF ( FRACLOST > 1 ) THEN
!                     AMT_LOST = STT(I,J,L,NN) * ( 1d0 - EXP(-FRACLOST) )
!     &                    * FSOIL
!                  ELSE
!                     AMT_LOST = STT(I,J,L,NN) * FRACLOST * FSOIL
!                  ENDIF
!
!
!                  ! ND44 diagnostic: drydep flux [kg/s]
!                  ! ND44 diagnostic: drydep velocity [cm/s]
!                  IF ( ND44 > 0 ) THEN
!                     AD44(I,J,N,1) = AD44(I,J,N,1) + ( AMT_LOST/DTCHEM ) 
!                     AD44(I,J,N,2) = AD44(I,J,N,2) + SVEL * FSOIL
!                  ENDIF
!
!
!                  ! Subtract AMT_LOST from the STT array [kg]
!                  STT(I,J,L,NN)  = STT(I,J,L,NN) - AMT_LOST
!		
!                    
!            ENDIF               ! I and J within bounds, ALBD<0.4 and FSOIL>0
!         ENDDO                  ! M = LAND GRID BOXES
         ENDDO                  ! PBL layers
      ENDDO                     ! NUMDEP = Number of species that drydep

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE DRYFLXH2HD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: depvel
!
! !DESCRIPTION: Subroutine DEPVEL computes the dry deposition velocities using 
!  a resistance-in-series model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEPVEL( Input_Opt, State_Met, State_Chm, NPTS,   
     &                   RADIAT,    TEMP,      SUNCOS,    F0,        
     &                   HSTAR,     XMW,       AIROSOL,   USTAR,  
     &                   CZ1,       OBK,       CFRAC,     ZH,      
     &                   LSNOW,     DVEL,      ZO,        RHB,    
     &                   PRESSU,    W10                          )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input Options object
      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology state object

      INTEGER, INTENT(IN) :: NPTS             ! # of grid boxes = IIPAR*JJPAR
      REAL*8,  INTENT(IN) :: RADIAT (MAXIJ )  ! Solar radiation [W/m2]
      REAL*8,  INTENT(IN) :: TEMP   (MAXIJ )  ! Temperature [K]
      REAL*8,  INTENT(IN) :: SUNCOS (MAXIJ )  ! Cosine of solar zenith angle
      LOGICAL, INTENT(IN) :: AIROSOL(MAXDEP)  ! =T denotes aerosol species
      REAL*8,  INTENT(IN) :: F0     (MAXDEP)  ! Reactivity factor for oxidation
                                              !  of biological substances
      REAL*8,  INTENT(IN) :: HSTAR  (MAXDEP)  ! Henry's law constant
      REAL*8,  INTENT(IN) :: XMW    (MAXDEP)  ! Molecular weight [kg/mol]
      REAL*8,  INTENT(IN) :: USTAR  (MAXIJ )  ! Friction velocity [m/s]
      REAL*8,  INTENT(IN) :: CZ1    (MAXIJ )  ! Alt @ which Vd is computed [m]
      REAL*8,  INTENT(IN) :: OBK    (MAXIJ )  ! Monin-Obhukov length [m]
      REAL*8,  INTENT(IN) :: CFRAC  (MAXIJ )  ! Surface cloud fraction
      REAL*8,  INTENT(IN) :: ZH     (MAXIJ )  ! Roughness height [m]
      REAL*8,  INTENT(IN) :: RHB    (MAXIJ )  ! Relative humidity [%]
      REAL*8,  INTENT(IN) :: PRESSU (MAXIJ )  ! Surface pressure [hPa]
      REAL*8,  INTENT(IN) :: W10    (MAXIJ )  ! Wind speed @ 10m altitude [m/s]
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
! 
      REAL*8, INTENT(OUT) :: DVEL(MAXIJ,MAXDEP)   ! Drydep velocity [m/s]
!
! !REMARKS:
!  Need as landtype input for each grid square (I,J); see CMN_DEP_mod.F
!     IREG(I,J)       - # of landtypes in grid square
!     ILAND(I,J,LDT)  - Land type ID for element LDT =1, IREG(I,J)
!                        (could be from any source - mapped to deposition 
!                        surface ID in input unit 65)
!     IJUSE(I,J,LDT) - Fraction ((per mil) of gridbox area occupied by
!                       land type element LDT
!                                                                             .
!  Need as leaf area index; see CMN_DEP_mod.F
!     XLAI(I,J,LDT)  - Leaf Area Index of land type element LDT
!                                                                             .
!  Need as meteorological input for each grid square(I,J) (passed):
!     RADIAT(IJLOOP) - Solar radiation in W m-2
!     TEMP(IJLOOP)   - Surface air temperature in K
!     SUNCOS(IJLOOP) - Cosine of solar zenith angle
!     LSNOW(IJLOOP)  - Logical for snow and sea ice
!     RHB(IJLOOP)    - Relative humidity at the surface
!     PRESSU(IJLOOP) - Local surface pressure
!     W10(IJLOOP)    - 10m wind speed
!                                                                             .
!  Need as input for each species K (passed):
!     F0(K)          - reactivity factor for oxidation of biological substances
!     HSTAR(K)       - Henry's Law constant
!     XMW(K)         - Molecular weight (kg/mole) of species K
!                      (used to calculate molecular diffusivities)
!     AIROSOL(K)     - LOGICAL flag (T = aerosol species; 
!                                    F = gas-phase species)
!                                                                             .
!  Also need to call the following subroutines to read drydep input data:
!     READ_DRYDEP_INPUTS    - (in this module) Reads in Olson land type 
!                             indices, dry deposition land type indices, 
!                             default roughness heights, and polynomial 
!                             coefficients.  (This supersedes MODIN, RDDRYCF)
!     COMPUTE_OLSON_LANDMAP - (in olson_landmap_mod.F90).  Reads in the
!                             Olson land types at native resolution and re-bins
!                             them on-the-fly to the GEOS-Chem grid resolution.
!                             (This supersedes RDLAND)
!     "rdlai.f"             - reads Leaf Area Indices from files "lai**.global"
!                                                                             .
!  Some variables used in the subroutine (passed):
!     LRGERA(IJLOOP) T -> stable atmosphere; a high aerodynamic resistance
!                        (RA=1.E4 m s-1) is imposed; else RA is calculated
!     USTAR(IJLOOP)  - Friction velocity (m s-1)
!     CZ1(IJLOOP)    - Altitude (m) at which deposition velocity is computed
!     OBK(IJLOOP)    - Monin-Obukhov length (m): set to 1.E5 m under neutral 
!                      conditions
!     CFRAC(IJLOOP)  - Fractional cloud cover
!     ZH(IJLOOP)     - Mixing depth (m)
!                                                                             .
!  Some variables used in the subroutine:
!     MAXDEP         - the maximum number of species for which the dry 
!                      deposition calculation is done
!     ZO(LDT)        - Roughness height (m) for specific surface type indexed 
!                      by LDT
!     RSURFC(K,LDT)  - Bulk surface resistance (s m-1) for species K to 
!                      surface LDT
!     C1X(K)         - Total resistance to deposition (s m-1) for species K
!                                                                             .
!  Returned:
!     DVEL(IJLOOP,K) - Deposition velocity (m s-1) of species K
!                                                                             .
!  References:
!  ============================================================================
!     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
!       resistance model for gaseous deposition to vegetated surfaces,
!       Atmos. Environ. 21, 91-101, 1987.
!     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
!     Businger, J.A., et al., Flux-profile relationships in the atmospheric 
!       surface layer, J. Atmos. Sci., 28, 181-189, 1971.
!     Dwight, H.B., Tables of integrals and other mathematical data,
!       MacMillan, 1957.
!     Guenther, A., and 15 others, A global model of natural volatile
!       organic compound emissions, J. Geophys. Res., 100, 8873-8892, 1995.
!     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
!       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
!     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
!       hydrocarbons, and ozone over the Amazon forest during the wet season,
!       J.  Geophys. Res., 95, 16737-16754, 1990.
!     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
!       J. Geophys. Res., 97, 16473-16479, 1992.
!     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York, 1988.
!     Munger, J.W., and 8 others, Atmospheric deposition of reactive
!       nitrogen oxides and ozone in a temperate deciduous forest and a
!       sub-arctic woodland, J. Geophys. Res., in press, 1996.
!     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2, sulfate, and
!       HNO3 deposition velocities computed using regional landuse and
!       meteorological data, Atmos. Environ., 20, 949-964, 1986.
!     Wang, Y.H., paper in preparation, 1996.
!     Wesely, M.L, Improved parameterizations for surface resistance to
!       gaseous dry deposition in regional-scale numerical models, 
!       Environmental Protection Agency Report EPA/600/3-88/025,
!       Research Triangle Park (NC), 1988.
!     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
!
! !REVISION HISTORY: 
!** Contact: D.J. Jacob, Harvard U. (djj@io.harvard.edu)
!** Modularized by G.M. Gardner, Harvard U.
!** Version 3.2:   5/27/97
!** Version 3.2.1: 3/4/99   -- bug fix in expression for RT 
!** Version 3.2.2: 3/26/99  -- bug fix: specify a large Ra for aerosols
!** Version 3.2.3: 11/12/99 -- change Reynolds # criterion from 10 to 1
!                           -- force double precision w/ "D" exponents
!** Version 3.3:   5/8/00   -- bug fixes, cleanup, updated comments.
!** Version 3.4:   1/22/03  -- remove hardwire for CANOPYNOX
!** Version 3.5    7/21/03  -- Remove cap of surface resistance in RLUXX
!** Version 3.6    4/01/04  -- Now do drydep of DUST aerosol tracers
!** Version 3.7    4/20/04  -- Now also do drydep of SEASALT aerosol tracers
!** Version 3.8    4/13/05  -- Accounts for hygroscopic growth of SEASALT
!**                             aerosol tracers.  DUST aerosol tracers do
!**                             not grow hygroscopically.  Added RHB as
!**                             an input argument.
!** Version 3.9    5/25/05  -- Now restore GISS-specific code for GCAP model
!** Version 3.9.1  11/17/05 -- change Reynolds # criterion from 1 to 0.1
!  11 May 2011 - L. Jaegle   - Updated to use actual Sea level pressure instead
!                              of 1000 hPa
!                            - Modified to used Slinn & Slinn (1980) over Ocean
!                              surfaces
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  10 Jan 2012 - M. Payer    - Updated to use local surface pressure
!  09 Apr 2012 - R. Yantosca - Remove IJREG, IJLAND, IJUSE, XYLAI arrays and
!                              replace w/ IREG, ILAND, IUSE, XLAI
!  09 Apr 2012 - R. Yantosca - Remove reference to CMN_VEL_mod.F
!  09 Apr 2012 - R. Yantosca - Now use INTENT(IN), INTENT(OUT) for arguments
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  12 Dec 2012 - R. Yantosca - Now get ILAND, IUSE, IREG from State_Met
!  13 Dec 2012 - R. Yantosca - Now get XLAI from State_Met
!  31 May 2013 - R. Yantosca - Now pass State_Chm, for TOMAS
!  14 Jun 2013 - R. Yantosca - Now use Input_Opt%ITS_A_POPS_SIM
!  29 Aug 2013 - R. Yantosca - Bug fix: Skip to the next species if unless
!                              HSTAR>0 and XMW>0, or AIROSOL=t.  This avoids
!                              a floating-point invalid condition.
!  28 Jan 2014 - R. Yantosca - For TOMAS, don't hold A_RADI and A_DEN PRIVATE
!
!***********************************************************************
!   Changes from Version 3.2 to Version 3.3:                         ***
!   * We now suppress dry deposition over aerodynamically smooth     ***
!     surfaces.  The previous algorithm yielded negative numbers     ***
!     when u* was very small (due to the logarithm going negative).  ***
!     See the comments below for more information.                   ***
!   * Now eliminate obsolete variables ZLMO and SIH from the code.   ***
!   * Obsolete comments have been updated or removed.                ***
!***********************************************************************
!   Changes from version 3.1 to version 3.2:                         ***
!   * In unstable atmospheres with |ZLMO| < ZO, as can happen        ***
!    occasionally under very low wind conditions with tall canopies, ***
!    application of Monin-Obukhov similarity yields negative values  ***
!    for RA.  This was a problem in version 3.1.  In fact,           ***
!    Monin-Obukhov similarity does not apply under such conditions,  ***
!    so we now set RA to zero and let the boundary                   ***
!    resistance RB define the overall aerodynamic resistance.  Since *** 
!    RB varies inversely with U* it will impose a large aerodynamic  ***
!    resistance under very low wind conditions.                      ***
!   * The range of applicability of stability correction functions   ***
!    to Monin-Obukhov similarity has been extended to                ***
!    -2.5 < z/zMO < 1.5, based on Figure 2 of Businger et al. [1971].***
!    The range used to be -1 < z/zMO < 1 in version 3.1.             ***
!***********************************************************************
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

#if   defined( TOMAS )
      ! For TOMAS aerosol (win, 7/15/09)
      INTEGER  BIN
      REAL*8   SIZ_DIA(MAXIJ, IBINS), SIZ_DEN(MAXIJ, IBINS)
#endif

      ! Logical for snow and sea ice
      LOGICAL LSNOW(MAXIJ)

      ! Loop indices (bmy, 3/29/12)
      INTEGER :: I, J

      ! Added these to pass particle density & number to DUST_SFCRSII
      ! so as to avoid parallelization errors with TOMAS (bmy, 1/31/14)
      REAL*8  :: DIAM, DEN

      ! Logical switch for POPS in order to use octanol-air partitioning instead
      ! of Henry's Law for scaling of cuticular resitances (clf, 1/3/2011)
      LOGICAL :: IS_POPS

      ! Pointers
      INTEGER, POINTER :: IREG(:,:)
      INTEGER, POINTER :: ILAND(:,:,:)
      INTEGER, POINTER :: IUSE(:,:,:)
      REAL*8,  POINTER :: XLAI(:,:,:)

      !=================================================================
      ! DEPVEL begins here!
      !=================================================================

      ! Initialize pointers
      IREG  => State_Met%IREG
      ILAND => State_Met%ILAND
      IUSE  => State_Met%IUSE
      XLAI  => State_Met%XLAI

      ! Is this a POPs simmulation?
      IS_POPS = Input_Opt%ITS_A_POPS_SIM  ! clf, 1/3/2011

!***********************************************************************       
!
!
!** If LDEP(K)=F, species does not deposit.
!** Deposition is applied only to species with LDEP=T.
      DO K = 1,NUMDEP

         ! Better test for depositing species K: We need both HSTAR and XMW
         ! to be nonzero, OR the value of AIROSOL to be true.  This should
         ! avoid any futher floating point invalid issues caused by putting
         ! a zero value in a denominator. (bmy, 8/29/13)
         ! (bmy, 8/29/13)
         IF ( ( HSTAR(K) > 0d0   .and. 
     &          XMW  (K) > 0d0 ) .or.  AIROSOL(K) ) THEN
            LDEP(K) = .TRUE.
         ELSE
            LDEP(K) = .FALSE.
         ENDIF
      ENDDO

      DO K = 1,NUMDEP    
         DO IJLOOP =1,NPTS
            DVEL(IJLOOP,K) = 0.0D0
         ENDDO
      ENDDO
!***********************************************************************
!*                                 
!*    Begin section for computing deposition velocities           
!*                                 
*                                 
#if   defined( TOMAS )
      !=================================================================
      ! FOR TOMAS MICROPHYSICS:
      !
      ! Calculate 30-bin aerosol size and diameter at each grid box
      ! This is done differently than 2-bin seasalt and 4-bin dust
      ! because at each box the aerosol size&density is different
      ! depending on the internally-mixed composition (win, 7/15/09)
      !=================================================================
      IF ( IDTNK1 > 0 ) THEN
         CALL AERO_DIADEN( SIZ_DIA,   SIZ_DEN,   1, 
     &                     Input_Opt, State_Met, State_Chm )
      ENDIF
#endif

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
!$OMP+PRIVATE( C1X,    VK,     I,      J,      IW             )
!$OMP+PRIVATE( DIAM,   DEN                                    )
#if defined( TOMAS )
!$OMP+PRIVATE( BIN                                            )
#endif
#if defined( GCAP )
!$OMP+PRIVATE( ZO                                             )
#endif

      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Convert (I,J) to 1-D index IJLOOP
         IJLOOP = ( ( J - 1 ) * IIPAR ) + I 

!** CZ is Altitude (m) at which deposition velocity is computed
         CZ = CZ1(IJLOOP)
!** TEMPK and TEMPC are surface air temperatures in K and in C
         TEMPK = TEMP(IJLOOP)
         TEMPC = TEMP(IJLOOP)-273.15D0
!* Initialize variables
         DO K = 1,NUMDEP
            VD(K) = 0.0D0
            DO LDT = 1,NTYPE
               RSURFC(K,LDT) = 0.D0
            END DO
         END DO

!** Calculate the kinematic viscosity XNU (m2 s-1) of air
!** as a function of temperature.
!** The kinematic viscosity is used to calculate the roughness heights over
!** water surfaces and to diagnose whether such surfaces are aerodynamically
!** rough or smooth using a Reynolds number criterion.
!** The expression for the temperature dependence of XNU
!** is from the FORTRAN code in Appendix II of Wesely [1988];
!** I wasn't able to find an original reference but it seems benign enough.
         C1 = TEMPK/273.15D0
         XNU = 0.151D0*(C1**1.77D0)*1.0D-04

!* Compute bulk surface resistance for gases.    
!*                                 
!* Adjust external surface resistances for temperature; 
!* from Wesely [1989], expression given in text on p. 1296.        
!*                                 
!* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
!* so the inner parentheses are not needed (bmy, 3/4/99)
!*        RT = 1000.0*EXP(-(TEMPC-4.0))
         RT = 1000.0D0*EXP(-TEMPC-4.0D0)
!*
!    Get surface resistances - loop over land types LDT
!***************************************************************************
!* The land types within each grid square are defined using the Olson
!* land-type database.  Each of the Olson land types is assigned a 
!* corresponding "deposition land type" with characteristic values of surface
!* resistance components.  There are 74 Olson land-types but only 11 deposition
!* land-types (i.e., many of the Olson land types share the same deposition
!* characteristics).  Surface resistance components for the "deposition land 
!* types" are from Wesely [1989] except for tropical forests [Jacob and Wofsy,
!* 1990] and for tundra [Jacob et al., 1992].  All surface resistance 
!* components are normalized to a leaf area index of unity.
!*
!* Olson land types, deposition land types, and surface resistance components
!* are read from file 'drydep.table'; check that file for further details.
!****************************************************************************

         ! Loop over the # of Olson land types in this grid box (I,J)
         DO 170 LDT = 1, IREG(I,J)

            ! If the land type is not represented in grid 
            ! box  (I,J), then skip to the next land type
            IF ( IUSE(I,J,LDT) == 0 ) GOTO 170

            ! Olson land type index + 1
            IOLSON = ILAND(I,J,LDT)+1

            ! Dry deposition land type index
            II     = IDEP(IOLSON)
!
!** If the surface to be snow or ice;
!** set II to 1 instead.
!
            IF(LSNOW(IJLOOP)) II=1

!* Read the internal resistance RI (minimum stomatal resistance for water 
!* vapor,per unit area of leaf) from the IRI array; a '9999' value means no 
!* deposition to stomata so we impose a very large value for RI.

            RI(LDT) = DBLE(IRI(II))
            IF (RI(LDT)   .GE. 9999.D0) RI(LDT)   = 1.D12

!** Cuticular resistances IRLU read in from 'drydep.table'
!** are per unit area of leaf;
!** divide them by the leaf area index to get a cuticular resistance for the
!** bulk canopy.  If IRLU is '9999' it means there are no cuticular
!** surfaces on which to deposit so we impose a very large value for RLU.
            IF ( IRLU(II) >= 9999 .or. XLAI(I,J,LDT) <= 0.D0 ) THEN
               RLU(LDT) = 1.D6
            ELSE
               RLU(LDT) = DBLE( IRLU(II) ) / XLAI(I,J,LDT) + RT
            ENDIF
!** The following are the remaining resistances for the Wesely
!** resistance-in-series model for a surface canopy
!** (see Atmos. Environ. paper, Fig.1).  
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
!***************************************************************************
!*                                 
!*    Adjust stomatal resistances for insolation and temperature:  
!*
!*     Temperature adjustment is from Wesely [1989], equation (3).
!*       
!*     Light adjustment by the function BIOFIT is described by Wang [1996].
!*     It combines
!*       - Local dependence of stomal resistance on the intensity I of light 
!*         impinging the leaf; this is expressed as a mutliplicative 
!*         factor I/(I+b) to the stomatal resistance where b = 50 W m-2
!*         (equation (7) of Baldocchi et al. [1987])
!*       - radiative transfer of direct and diffuse radiation in the 
!*         canopy using equations (12)-(16) from Guenther et al. [1995]
!*       - separate accounting of sunlit and shaded leaves using
!*         equation (12) of Guenther et al. [1995]
!*       - partitioning of the radiation at the top of the canopy into direct
!*         and diffuse components using a parameterization to results from
!*         an atmospheric radiative transfer model [Wang, 1996]
!*     The dependent variables of the function BIOFIT are the leaf area
!*     index (XYLAI), the cosine of zenith angle (SUNCOS) and the fractional
!*     cloud cover (CFRAC).  The factor GFACI integrates the light
!*     dependence over the canopy depth; sp even though RI is input per
!*     unit area of leaf it need not be scaled by LAI to yield a bulk
!*     canopy value because that's already done in the GFACI formulation.
!***************************************************************************

            RAD0 = RADIAT(IJLOOP)
            RIX = RI(LDT)
            IF (RIX .GE. 9999.D0) GO TO 150
            GFACT = 100.0D0
            IF (TEMPC .GT. 0.D0 .AND. TEMPC .LT. 40.D0)
     *           GFACT = 400.D0/TEMPC/(40.0D0-TEMPC)
            GFACI = 100.D0
            IF ( RAD0 > 0.D0 .and. XLAI(I,J,LDT) > 0.D0 ) THEN
               GFACI = 1.D0 / BIOFIT( DRYCOEFF,       XLAI(I,J,LDT),
     &                                SUNCOS(IJLOOP), CFRAC(IJLOOP) )
            ENDIF
            
            RIX = RIX*GFACT*GFACI
 150        CONTINUE
!*                                 
!*    Compute aerodynamic resistance to lower elements in lower part           
!*    of the canopy or structure, assuming level terrain - 
!*    equation (5) of Wesely [1989].
!*                                 
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))
!*
!*    Loop over species; species-dependent corrections to resistances
!*    are from equations (6)-(9) of Wesely [1989].
!*
            DO 160  K = 1,NUMDEP
!** exit for non-depositing species or aerosols.
               IF (.NOT. LDEP(K) .OR. AIROSOL(K)) GOTO 155

               XMWH2O = 18.D-3
               RIXX = RIX*DIFFG(TEMPK,PRESSU(IJLOOP),XMWH2O)/
     C              DIFFG(TEMPK,PRESSU(IJLOOP),XMW(K))
     C              + 1.D0/(HSTAR(K)/3000.D0+100.D0*F0(K))
               RLUXX = 1.D12
               IF (RLU(LDT).LT.9999.D0)
     C              RLUXX = RLU(LDT)/(HSTAR(K)/1.0D+05 + F0(K))

        ! If a POPs simulation, scale cuticular resistances with octanol-air
        ! partition coefficient (Koa) instead of HSTAR (clf, 1/3/2011)

               IF (IS_POPS)
     C              RLUXX = RLU(LDT)/(KOA(K)/1.0D+05 + F0(K))

!*
!* To prevent virtually zero resistance to species with huge HSTAR, such
!* as HNO3, a minimum value of RLUXX needs to be set. The rationality
!* of the existence of such a minimum is demonstrated by the observed
!* relationship between Vd(NOy-NOx) and Ustar in Munger et al.[1996];
!* Vd(HNO3) never exceeds 2 cm s-1 in observations. The
!* corresponding minimum resistance is 50 s m-1.  This correction
!* was introduced by J.Y. Liang on 7/9/95.
!*

               RGSX = 1.D0/(HSTAR(K)/1.0D+05/RGSS(LDT) + 
     1              F0(K)/RGSO(LDT))
               RCLX = 1.D0/(HSTAR(K)/1.0D+05/RCLS(LDT) + 
     1              F0(K)/RCLO(LDT))
!*
!** Get the bulk surface resistance of the canopy, RSURFC, from the network
!** of resistances in parallel and in series (Fig. 1 of Wesely [1989])
               DTMP1=1.D0/RIXX
               DTMP2=1.D0/RLUXX
               DTMP3=1.D0/(RAC(LDT)+RGSX)
               DTMP4=1.D0/(RDC+RCLX)
               RSURFC(K,LDT) = 1.D0/(DTMP1 + DTMP2 + DTMP3 + DTMP4)
!** get surface deposition velocity for aerosols if needed;
!** equations (15)-(17) of Walcek et al. [1986]
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
     &               AERO_SFCRSII( K,     II,   PRESSU(IJLOOP)*1D-3, 
     &                             TEMPK, USTAR(IJLOOP), RHB(IJLOOP),
     &                             W10(IJLOOP), Input_Opt )

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

                  ! Particle diameter [m]
                  DIAM  = A_RADI(K) * 2.d0 

                  ! Particle density [kg/m3]
                  DEN   = A_DEN(K) 

                  ! [Zhang et al., 2001]
                  ! Now pass DIAM, DEN as arguments (bmy, 1/31/14)
                  RSURFC(K,LDT) = DUST_SFCRSII( K, II, 
     &                                          PRESSU(IJLOOP)*1D-3, 
     &                                          TEMPK, 
     &                                          USTAR(IJLOOP),
     &                                          DIAM,  
     &                                          DEN )

#if   defined( TOMAS )
               ! (win, 7/15/09)
               ELSE IF ( ( NTRAIND(K) >= IDTNK1  )   .AND.
     &                   ( NTRAIND(K) <  IDTNK1+IBINS ) ) THEN

                  !=====================================================
                  ! Use size-resolved dry deposition calculations for 
                  ! size-resolved aerosols.  Choose to use DUST_SFCRSII
                  ! because these aerosol have diameter and density 
                  ! calculation that takes hygroscopic growth into 
                  ! account already, thus no need for such calculation
                  ! like in AERO_SFCRSII.  (win, 7/19/07)
                  !=====================================================     

                  ! Get current aerosol diameter and density
                  ! NOTE: In TOMAS the aerosol diameter and density
                  ! evolves with time.  We have to save these in the 
                  ! DIAM and DEN variables so that we can hold these
                  ! PRIVATE for the parallel loop.
                  BIN  = NTRAIND(K) - IDTNK1 + 1
                  DIAM = SIZ_DIA( IJLOOP, BIN )     ! Diameter [m]
                  DEN  = SIZ_DEN( IJLOOP, BIN )     ! Density  [kg/m3]

                  ! [Zhang et al., 2001]
                  ! Now pass DIAM, DEN as arguments (bmy, 1/31/14)
                  RSURFC(K,LDT) = DUST_SFCRSII( K, II, 
     &                                          PRESSU(IJLOOP)*1D-3,
     &                                          TEMPK, 
     &                                          USTAR(IJLOOP),   
     &                                          DIAM, 
     &                                          DEN )
#endif

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
     &                ADUST_SFCRSII(K, II, PRESSU(IJLOOP)*1D-3, TEMPK,
     &                              USTAR(IJLOOP))
                  
!*                                 
!*    Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
!*    VDSMAX is taken from Table 2 of Walcek et al. [1986].
!*    Invert to get corresponding R

!                  RSURFC(K,LDT) = 1.D0/MIN(VDS, DBLE(IVSMAX(II))/1.D4)
               ENDIF
 160        CONTINUE
!*
 170     CONTINUE
!*
!*    Set max and min values for bulk surface resistances         
!*                                 
         DO 190 K = 1,NUMDEP 
            IF (.NOT.LDEP(K)) GOTO 190
            DO 180 LDT = 1, IREG(I,J)
               IF ( IUSE(I,J,LDT) == 0 ) GOTO 180
               RSURFC(K,LDT)= MAX(1.D0, MIN(RSURFC(K,LDT), 9999.D0))
 180        CONTINUE
 190     CONTINUE
!*                                 
!*    Loop through the different landuse types present in the grid square     
!*           
         DO 500 LDT = 1, IREG(I,J)
            IF ( IUSE(I,J,LDT) == 0 ) GOTO 500
            IOLSON = ILAND(I,J,LDT) + 1

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

!***** Get aerodynamic resistances Ra and Rb. ***********************
!   The aerodynamic resistance Ra is integrated from altitude z0+d up to the
!   altitude z1 at which the dry deposition velocity is to be referenced.
!   The integration corrects for stability using Monin-Obukhov similarity 
!   formulas from Businger et al. [1971] which apply over the range 
!   -2.5 < z/zMO < 1.5 (see their Figure 2).
!   Under very unstable conditions when z1 > -2.5 zMO, we assume that there is
!   no resistance to transfer in the convective column between zMO and z1.
!   Under very stable conditions when z1 > 1.5 zMO, we assume that vertical
!   transfer in the column between zMO and z1 is strongly suppressed so
!   that the deposition velocity at altitude z1 is very low.  Under these
!   conditions we just specify a very large Ra=1.E4 s m-1 (LRGERA = T).
!**
!   The Reynolds number REYNO diagnoses whether a surface is
!   aerodynamically rough (REYNO > 1) or smooth.  
!
!   NOTE: The criterion "REYNO > 1" was originally "REYNO > 10".
!         See below for an explanation of why it was changed (hyl, 10/15/99)
!
!   Surface is rough in all cases except over water with low wind speeds.  
!   In the smooth case, vertical transport IN THE SUBLAYER near the surface 
!   is limited by molecular diffusion and is therefore very slow; we assign 
!   a large value we assign a large value of Ra + Rb to account for this 
!   effect.  [In Versions 3.2 and earlier we used the formulation for Ra + Rb 
!   given in Equation (12) of Walcek et al [1986] to calculate the aerodynamic 
!   resistance over smooth surfaces.  However, that expression fails when 
!   u* is very small, as it yields negative values of Ra + Rb].  
!   (djj, hyl, bmy, 5/8/00)
!**
!   In the aerodynamically rough case, the expression for Ra is as
!   given in equation (5) of Jacob et al. [1992]:
!
!          Ra = (1/ku*)*int(from z0 to z1) (phi(x)/z)dz
!
!    where x = (z-D)/zMO, z is the height above ground, and D is the
!    displacement height which is typically 70-80% of the canopy height
!    [Brutsaert, 1982].  We change the vertical coordinate so that z=0 at
!    the displacement height; that's OK since for all practical applications
!    z1 >> D.  In this manner we don't need to assume any specific value for
!    the displacement height.  Applying the variable transformation 
!    z -> x = z/zMO, the equation above becomes
!
!          Ra = (1/ku*)*int(from x0 to x1) (phi(x)/x)dx   with x=z/zMO
!
!    Here phi is a stability correction function originally formulated by
!    Businger et al. [1971] and given in eqns 5a and 5b of Jacob et al. [1992].
!    For unstable conditions,
!
!          phi(x) = a/sqrt(1-bx)  where a=0.74, b = 9
!
!    The analytical solution to the integral is 
!    [Dwight, 1957, integral 192.11]:
!
!          int(dx/(x*sqrt(1-bx))) = log(abs((sqrt(1-bx)-1)/(sqrt(1-bx)+1)))
!
!    which yields the expression for Ra used in the code for unstable 
!    conditions.  For stable conditions,
!
!          phi(x) = a + bx        where a=0.74, b = 4.7
!
!    and the analytical solution to the integral is
!
!          int((a/x)+b)dx = a*ln(x) + bx
!
!    which yields the expression of Ra used in the code for stable conditions.
!**
!   The formulation of RB for gases is equation (12) of 
!   Walcek et al. [1986].  The parameterization for deposition of
!   aerosols does not include an RB term so RB for aerosols is set
!   to zero.
!   Modify phi(x) according to the non-local mixing scheme by Holtslag 
!   and Boville [1993] ( Lin, 07/18/08 )
!   For unstable conditions,
!          phi(x) = a/sqrt(1-bx)  where a=1.0, b=15.0
!
!   For stable conditions,
!          phi(x) = a + bx
!              where a=1.0, b=5.0 for 0 <= x <= 1, and
!                    a=5.0, b=1.0 for x > 1.0 
!*********************************************************************
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
!*                                 
            IF (CKUSTR.EQ.0.0D0) THEN
               WRITE(6,212) IJLOOP,CKUSTR,XCKMAN,USTAR(IJLOOP)
 212           FORMAT(1X,'IJLOOP= ',I4,1X,'CKUSTR=',E10.1,1X,
     x              'XCKMAN= ',E12.4,1X,'USTAR(IJLOOP)= ',
     x              E12.4)
               CLOSE(98)
               STOP             ! debug
            ENDIF
!
!
!...aerodynamically rough or smooth surface
! "In the classic study by Nikuradse (1933) the transition from smooth
!  to rough was examined in pipe flow. He introduced a roughness Reynolds
!  number Rr = U* Z0 / Nu and found the flow to be smooth for Rr < 0.13
!  and rough for Rr > 2.5 with a transition regime in between."
!  (E.B. Kraus and J.A. Businger, Atmosphere-Ocean Interaction, second
!  edition, P.144-145, 1994). Similar statements can be found in the books:
!  Evaporation into the atmosphere, by Wilfried Brutsaert, P.59,89, 1982;
!  or Seinfeld & Pandis, P.858, 1998. Here we assume a sudden transition
!  point Rr = 1 from smooth to rough, following L. Merlivat (1978, The
!  dependence of bulk evaporation coefficients on air-water interfacial
!  conditions as determined by the isotopic method, J. Geophys. Res.,
!  Oceans & Atmos., 83, C6, 2977-2980). Also refer to Brutsaert's book,
!  P.125. We used to use the criterion "REYNO > 10" for aerodynamically
!  rough surface and now change to "REYNO > 1". (hyl, 10/15/99)
!  
!  11/17/05: D. J. Jacob says to change the criterion for aerodynamically
!  rough surface to REYNO > 0.1 (eck, djj, bmy, 11/17/05)
            IF ( REYNO < 0.1d0 ) GOTO 220

            ! Add option for non-local PBL (Lin, 03/31/09) 
            IF (.NOT. LNLPBL) THEN

!...aerodynamically rough surface.
!*                                 
               IF (CORR1.LE.0.0D0 .AND. Z0OBK .LT. -1.D0)THEN
!*... unstable condition; set RA to zero. (first implemented in V. 3.2)
                  RA = 0.D0
               ELSEIF (CORR1.LE.0.0D0 .AND. Z0OBK .GE. -1.D0) THEN
!*... unstable conditions; compute Ra as described above.
                  DUMMY1 = (1.D0 - 9D0*CORR1)**0.5D0
                  DUMMY2 = (1.D0 - 9D0*Z0OBK)**0.5D0
                  DUMMY3 = ABS((DUMMY1 - 1.D0)/(DUMMY1 + 1.D0))
                  DUMMY4 = ABS((DUMMY2 - 1.D0)/(DUMMY2 + 1.D0))
                  RA = 0.74D0* (1.D0/CKUSTR) * LOG(DUMMY3/DUMMY4)
                                 
               ELSEIF((CORR1.GT.0.0D0).AND.(.NOT.LRGERA(IJLOOP)))THEN
!*...moderately stable conditions (z/zMO <1); compute Ra as described above
                  RA = (1D0/CKUSTR) * 
     &                 (.74D0*LOG(CORR1/Z0OBK) + 4.7D0*(CORR1-Z0OBK))
               ELSEIF(LRGERA(IJLOOP)) THEN
!*... very stable conditions
                  RA = 1.D+04
               ENDIF
!* check that RA is positive; if RA is negative (as occasionally
!* happened in version 3.1) send a warning message.

            ELSE

               IF (CORR1.LT.0.0D0) THEN
!*... unstable conditions; compute Ra as described above.
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
!* Get total resistance for deposition - loop over species.
            DO 215 K = 1,NUMDEP 
               IF (.NOT.LDEP(K)) GOTO 215
!** DAIR is the thermal diffusivity of air; value of 0.2*1.E-4 m2 s-1 
!** cited on p. 16,476 of Jacob et al. [1992]
               DAIR = 0.2D0*1.D-4
               RB = (2.D0/CKUSTR)*
     x              (DAIR/DIFFG(TEMPK,PRESSU(IJLOOP),XMW(K)))**0.667D0
               IF (AIROSOL(K)) RB=0.D0
               C1X(K) = RA + RB + RSURFC(K,LDT)
 215        CONTINUE
            GOTO 240
 220        CONTINUE 
!** ... aerodynamically smooth surface
!** BUG FIX -- suppress drydep over smooth surfaces by setting Ra to a large
!** value (1e4).  This prevents negative dry deposition velocities when u*
!** is very small (djj, bmy, 5/8/00)
            DO 230 K = 1,NUMDEP 
               IF ( LDEP(K) ) THEN
                  RA     = 1.0D4
                  C1X(K) = RA + RSURFC(K,LDT)
               ENDIF
 230        CONTINUE

 240        CONTINUE
!*                                 
!* IJUSE is the fraction of the grid square occupied by surface LDT
!* in units of per mil (IJUSE=500 -> 50% of the grid square).  Add
!* the contribution of surface type LDT to the deposition velocity;
!* this is a loop over all surface types in the gridbox.
!*
            DO 400 K = 1,NUMDEP
               IF (.NOT.LDEP(K)) GOTO 400
               VK(K) = VD(K)
               VD(K) = VK(K) +.001D0* DBLE( IUSE(I,J,LDT) )/C1X(K)
 400        CONTINUE
 500     CONTINUE

!** Load array DVEL
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
#if !defined( GCAP )
               PRINT*, 'ZO       = ', ZO(IJLOOP)
#endif
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
#if !defined( GCAP )
               PRINT*, 'ZO       = ', ZO(IJLOOP)
#endif
               CALL CLEANUP
               STOP
!$OMP END CRITICAL
            ENDIF
 550     CONTINUE
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      NULLIFY( IREG  )
      NULLIFY( ILAND )
      NULLIFY( IUSE  )
      NULLIFY( XLAI  )

      END SUBROUTINE DEPVEL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diffg
!
! !DESCRIPTION: Subroutine DIFFG calculates the molecular diffusivity [m2/s] in
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].  (bmy, 5/16/06)
!\\
!\\
! !INTERFACE:
!
      FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
!
! !USES:
!

!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: TK     ! Temperature [K]
      REAL*8, INTENT(IN) :: PRESS  ! Pressure [Pa]
      REAL*8, INTENT(IN) :: XM     ! Molecular weight of gas [kg]
!
! !REMARKS:
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
! 
! !REVISION HISTORY: 
!  (1 ) Originally was a standalone function; now bundled into drydep_mod.f.
!        Also now force REAL*8 precision with D exponents.  Now use F90
!        style syntax and updated comments. (bmy, 5/16/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      ! Air density [molec/m3]
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

      END FUNCTION DIFFG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_drydep_inputs
!
! !DESCRIPTION: Subroutine READ\_DRYDEP\_INPUTS reads inputs for the dry
!  deposition module corresponding to either the Olson 1992 (GEOS-Chem default)
!  or Olson 2001 (planned replacement for Olson 1992) land map.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_DRYDEP_INPUTS( am_I_Root, DRYCOEFF, IOLSON, 
     &                               IDEP,      IWATER,   NWATER, 
     &                               IZO,       IDRYDEP,  IRI,    
     &                               IRLU,      IRAC,     IRGSS,  
     &                               IRGSO,     IRCLS,    IRCLO,  
     &                               IVSMAX,    Data_Dir_1x1      )
!
! !USES:
!
      ! Modules for netCDF read
      USE m_netcdf_io_open
      USE m_netcdf_io_get_dimlen
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
      CHARACTER(LEN=255), INTENT(IN) :: Data_Dir_1x1
!
! !OUTPUT PARAMETERS:
!  
      !-----------------------------------------------------------------------
      ! DRYCOEFF : Baldocchi polynomial coeffs
      ! IOLSON   : Olson land type indices (+1) 
      ! IDEP     : Mapping: Olson ==> drydep ID
      ! IWATER   : Olson types that represent water
      ! NWATER   : Number of Olson types that are water
      ! IZO      : Default Z0 (routgness height) for each Olson land type
      ! IDRYDEP  : Dry deposition land type indices
      ! IRI      : RI   resistance for drydep
      ! IRLU     : RLU  resistance for drydep
      ! IRAC     : RAC  resistance for drydep
      ! IRGSS    : RGSS resistance for drydep
      ! IRGSO    : RGSO resistance for drydep
      ! IRCLS    : RCLS resistance for drydep
      ! IRCLO    : RCLO resistance for drydep
      ! IVSMAX   : Max drydep velocity (for aerosol) perr drydep land type
      !-----------------------------------------------------------------------
      REAL*8,  INTENT(OUT) :: DRYCOEFF(NPOLY    ) 
      INTEGER, INTENT(OUT) :: IOLSON  (NVEGTYPE )       
      INTEGER, INTENT(OUT) :: IDEP    (NVEGTYPE )      
      INTEGER, INTENT(OUT) :: IWATER  (NVEGTYPE )    
      INTEGER, INTENT(OUT) :: NWATER        
      INTEGER, INTENT(OUT) :: IZO     (NVEGTYPE )       
      INTEGER, INTENT(OUT) :: IDRYDEP (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRI     (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRLU    (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRAC    (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRGSS   (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRGSO   (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRCLS   (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IRCLO   (NDRYDTYPE)
      INTEGER, INTENT(OUT) :: IVSMAX  (NDRYDTYPE)
!
! !REMARKS:
!  Routine READ_DRYDEP_INPUTS replaces routines MODIN (which read the ASCII 
!  file "drydep.table") and RDDRYCF (which read the ASCII file "drydep.coef").
!                                                                             .
!  READ_DRYDEP_INPUTS was generated from the Perl script "ncCodeRead", which
!  is part of the NcdfUtilities package (with subsequent hand-editing).
!                                                                             .
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!
! !REVISION HISTORY:
!  26 Mar 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! Variable declarations
      !=================================================================

      ! Scalars
      INTEGER            :: fId                 ! netCDF file ID
      CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
      CHARACTER(LEN=255) :: nc_file             ! netCDF file name
      CHARACTER(LEN=255) :: nc_path             ! netCDF path name
      CHARACTER(LEN=255) :: v_name              ! netCDF variable name 
      CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
      CHARACTER(LEN=255) :: a_val               ! netCDF attribute value

      ! Arrays for netCDF start and count values
      INTEGER            :: st1d(1), ct1d(1)    ! For 1D arrays 

      !=================================================================
      ! Open and read data from the netCDF file
      !=================================================================

      ! Pick the proper input file for Olson 1992 or Olson 2001
      IF ( USE_OLSON_2001 ) THEN
         nc_file = 'Olson_2001_Drydep_Inputs.nc'
      ELSE
         nc_file = 'Olson_1992_Drydep_Inputs.nc'
      ENDIF

      ! Open file for reading (construct full data path)
      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'Olson_Land_Map_201203/'
      nc_path = TRIM( nc_dir )       // TRIM( nc_file )
      CALL Ncop_Rd( fId, TRIM(nc_path) )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) REPEAT( '%', 79 )
         WRITE( 6, 110 ) TRIM(nc_file)
         WRITE( 6, 120 ) TRIM(nc_dir)
      ENDIF

      !----------------------------------------
      ! VARIABLE: DRYCOEFF
      !----------------------------------------

      ! Variable name
      v_name = "DRYCOEFF"

      ! Read DRYCOEFF from file
      st1d   = (/ 1      /)
      ct1d   = (/ NPOLY /)
      CALL NcRd( DRYCOEFF, fId, TRIM(v_name), st1d, ct1d )

      ! Read the DRYCOEFF:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IOLSON
      !----------------------------------------

      ! Variable name
      v_name = "IOLSON"

      ! Read IOLSON from file
      st1d   = (/ 1         /)
      ct1d   = (/ NVEGTYPE /)
      CALL NcRd( IOLSON, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IOLSON:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IDEP
      !----------------------------------------

      ! Variable name
      v_name = "IDEP"

      ! Read IDEP from file
      st1d   = (/ 1         /)
      ct1d   = (/ NVEGTYPE /)
      CALL NcRd( IDEP, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IDEP:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IWATER
      !----------------------------------------

      ! Variable name
      v_name = "IWATER"

      ! Get the # of Olson types that are water
      ! (NOTE: IWATER is an index array, dimension name = variable name)
      CALL NcGet_DimLen( fId, TRIM(v_name), NWATER )

      ! Initialize
      IWATER = 0

      ! Read IWATER from file
      ! NOTE: IWATER is declared with NNVEGTYPE, but has NWATER values
      ! The rest can be zeroed out
      st1d   = (/ 1      /)
      ct1d   = (/ NWATER /)
      CALL NcRd( IWATER(1:NWATER), fId, TRIM(v_name), st1d, ct1d )

      ! Read the IWATER:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IZO
      !----------------------------------------

      ! Variable name
      v_name = "IZO"

      ! Read IZO from file
      st1d   = (/ 1         /)
      ct1d   = (/ NVEGTYPE /)
      CALL NcRd( IZO, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IZO:long_name attribute
      a_name = "long_name"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IDRYDEP
      !----------------------------------------

      ! Variable name
      v_name = "IDRYDEP"

      ! Read IDRYDEP from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IDRYDEP, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IDRYDEP:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRI
      !----------------------------------------

      ! Variable name
      v_name = "IRI"

      ! Read IRI from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRI, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRI:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRLU
      !----------------------------------------

      ! Variable name
      v_name = "IRLU"

      ! Read IRLU from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRLU, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRLU:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRAC
      !----------------------------------------

      ! Variable name
      v_name = "IRAC"

      ! Read IRAC from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRAC, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRAC:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRGSS
      !----------------------------------------

      ! Variable name
      v_name = "IRGSS"

      ! Read IRGSS from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRGSS, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRGSS:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRGSO
      !----------------------------------------

      ! Variable name
      v_name = "IRGSO"

      ! Read IRGSO from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRGSO, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRGSO:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRCLS
      !----------------------------------------

      ! Variable name
      v_name = "IRCLS"

      ! Read IRCLS from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRCLS, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRCLS:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IRCLO
      !----------------------------------------

      ! Variable name
      v_name = "IRCLO"

      ! Read IRCLO from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IRCLO, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IRCLO:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: IVSMAX
      !----------------------------------------

      ! Variable name
      v_name = "IVSMAX"

      ! Read IVSMAX from file
      st1d   = (/ 1         /)
      ct1d   = (/ NDRYDTYPE /)
      CALL NcRd( IVSMAX, fId, TRIM(v_name), st1d, ct1d )

      ! Read the IVSMAX:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Close netCDF file
      CALL NcCl( fId )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 140 )
         WRITE( 6, 100 ) REPEAT( '%', 79 )
      ENDIF

      ! FORMAT statements
 100  FORMAT( a                                              )
 110  FORMAT( '%% Opening file  : ',         a               )
 120  FORMAT( '%%  in directory : ',         a, / , '%%'     )
 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140  FORMAT( '%% Successfully closed file!'                 )

      END SUBROUTINE READ_DRYDEP_INPUTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aero_sfcrsii
!
! !DESCRIPTION: Function AERO\_SFCRSII computes the aerodynamic resistance of
!  seasalt aerosol tracers according to Zhang et al 2001.  We account for
!  hygroscopic growth of the seasalt aerosol particles. (rjp, tdf, bec, bmy,
!  4/1/04, 6/11/08)
!\\
!\\
! !INTERFACE:
!
      FUNCTION AERO_SFCRSII( K, II, PRESS, TEMP, USTAR, RHB, 
     &                       W10, Input_Opt ) RESULT(RS)
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: K     ! Drydep tracer index (range: 1-NUMDEP)
      INTEGER,        INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8,         INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb=100 Pa=0.1 kPa)
      REAL*8,         INTENT(IN) :: TEMP  ! Temperature [K]   
      REAL*8,         INTENT(IN) :: USTAR ! Friction velocity [m/s]
      REAL*8,         INTENT(IN) :: RHB   ! Relative humidity (fraction)
      REAL*8,         INTENT(IN) :: W10   ! 10 m windspeed [m/s]
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL*8              :: RS    ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY: 
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Now limit relative humidity to [tiny(real*8),0.99] range for DLOG
!         argument (phs, 6/11/08)
!  (3 ) Bug fixes to the Gerber (1985) growth function (jaegle 5/11/11)
!  (4)  Update growth function to Lewis and Schwartz (2006) and density
!       calculation based on Tang et al. (1997) (bec, jaegle 5/11/11)
!  (5 ) Updates of sea salt deposition over water to follow the Slinn & Slinn
!       (1980) formulation over water surface. Described in Jaegle et al. (ACP,
!       11, 2011) (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  14 Jun 2013 - R. Yantosca - Now pass Input_Opt via the arg list

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      ! replace RCM with RUM (radius in microns instead of cm) - jaegle 5/11/11
      !REAL*8  :: DIAM, DEN, RATIO_R, RWET, RCM
      REAL*8  :: DIAM, DEN, RATIO_R, RWET, RUM
      REAL*8  :: FAC1, FAC2
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS
      ! New variables added (jaegle 5/11/11)
      REAL*8  :: SW
      REAL*8  :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW ! for weighting the settling velocity
      REAL*8  :: D0, D1  !lower and upper bounds of sea-salt dry diameter bins 
      REAL*8  :: DEDGE
      REAL*8  :: DEN1, WTP 
      INTEGER :: ID,NR
      LOGICAL, SAVE          :: FIRST = .TRUE.

      !increment of radius for integration of settling velocity (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2

      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) - jaegle 5/11/11
      REAL*8,  PARAMETER     :: A1 =  7.93d-3
      REAL*8,  PARAMETER     :: A2 = -4.28d-5
      REAL*8,  PARAMETER     :: A3 =  2.52d-6
      REAL*8,  PARAMETER     :: A4 = -2.35d-8
      REAL*8,  PARAMETER     :: EPSI = 1.0D-4

      ! parameters for assumed size distribution of accumulation and coarse
      ! mode sea salt aerosols, as described in Jaegle et al. (ACP, 11, 2011)
      ! (jaegle, 5/11/11)
      ! 1) geometric dry mean diameters (microns)
      REAL*8,  PARAMETER     ::   RG_A = 0.085d0
      REAL*8,  PARAMETER     ::   RG_C = 0.4d0
      ! 2) sigma of the size distribution
      REAL*8,  PARAMETER     ::   SIG_A = 1.5d0
      REAL*8,  PARAMETER     ::   SIG_C = 1.8d0
      REAL*8,  PARAMETER     :: PI =3.14159D0


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
      !      Update (jaegle 5/11/2011): The above formulation of Zhang et al
      !      (2001) is valid for land surfaces and was originally based on the
      !      work of Slinn (1982). Over water surfaces, the work of reference
      !      is that of Slinn and Slinn (1980) who use the term "viscous
      !      sublayer" to refer to the thin layer extending 0.1-1mm above the 
      !      water surface. Due to the proximity of the water, the RH in this
      !      layer is much higher than the ambient RH in the surface layer.
      !      According to Lewis and Schwartz (2004): "Relative humidities of
      !      99% and 100% were considered by Slinn and Slinn for the viscous
      !      sublayer, however near the ocean surface RH would be limited to
      !      near 98% because of the vapor pressure lowering of water over
      !      seawater due to the salt content". We will thus use a constant
      !      value RH=98% over all ocean boxes. This affects the growth of
      !      particles (the wet radius at RH=98% is x4 the dry radius) and thus
      !      affects all the terms depending on particle size.
      !      
      !      Other updates for ocean surfaces:
      !         a)   Over ocean surfaces the formulation from Slinn & Slinn for
      !              the resistance in the viscous layer is 
      !                Rs = 1 / (Cd/XCKMAN*U10m*(Eb+Eim)+VTS)
      !              with  Cd=(Ustar/U10m)**2, and VTS is the gravitational
      !              settling in the viscous layer. Note that the gravitational
      !              settling calculated here for the viscous layer is >> than
      !              the one calculated for the surface layer in seasalt_mod.f
      !              because of the higher RH.
      !         b)   Eim = 10^(-3/St) based on Slinn and Slinn (1980)
      !
      ! References:
      !  LEWIS and SCHWARTZ (2004), "SEA SALT AEROSOL PRODUCTION, MECHANISMS,
      !    METHODS AND MODELS" AGU monograph 152.
      !  SLINN and SLINN (1980), "PREDICTIONS FOR PARTICLE DEPOSITION ON
      !    NATURAL-WATERS" Atmos Environ (1980) vol. 14 (9) pp. 1013-1016.
      !  SLINN (1982), "PREDICTIONS FOR PARTICLE DEPOSITION TO VEGETATIVE
      !    CANOPIES" Atmos Environ (1982) vol. 16 (7) pp. 1785-1794.
      !==================================================================

      ! Number of bins for sea salt size distribution
      NR =INT((( Input_Opt%SALC_REDGE_um(2) - 
     &           Input_Opt%SALA_REDGE_um(1) ) / DR ) 
     &          + 0.5d0 )

      ! Particle radius [cm]
      ! Bug fix: The Gerber [1985] growth should use the dry radius 
      ! in micromenters and not cm. Replace RCM with RUM (jaegle 5/11/11)
      !RCM  = A_RADI(K) * 1.d2
      RUM  = A_RADI(K) * 1.d6
      
      ! Exponential factors used for hygroscopic growth
      ! Replace RCM with RUM (jaegle 5/11/11)
      !FAC1 = C1 * ( RCM**C2 )
      !FAC2 = C3 * ( RCM**C4 )
      FAC1 = C1 * ( RUM**C2 )
      FAC2 = C3 * ( RUM**C4 )

      ! Aerosol growth with relative humidity in radius [m] 
      ! (Gerber, 1985) (bec, 12/8/04)
      ! Added safety check for LOG (phs, 6/11/08)
      RHBL    = MAX( TINY(RHB), RHB )

      ! Check whether we are over the oceans or not:
      ! Over oceans the RH in the viscous sublayer is set to 98%, following
      ! Lewis and Schwartz (2004), see discussion above (jaegle 5/11/11)
      IF (LUC == 14) THEN
          RHBL = 0.98
      ENDIF
      ! Corrected bug in Gerber formulation: use of LOG10  (jaegle 5/11/11)
      !RWET    = 0.01d0*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.d0)**0.33d0
      !RWET = 1.d-6*(FAC1/(FAC2-LOG10(RHBL))+RUM**3.d0)**0.33333d0

      ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth [m] (jaegle 5/11/11)
      RWET = A_RADI(K) * (4.d0 / 3.7d0) *
     &          ( (2.d0 - RHBL)/(1.d0 - RHBL) )**(1.d0/3.d0)

      ! Ratio dry over wet radii at the cubic power
      !RATIO_R = ( A_RADI(K) / RWET )**3.d0

      ! Diameter of the wet aerosol [m]
      DIAM  = RWET * 2.d0  

      ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
      !DEN   = RATIO_R * A_DEN(K) + ( 1.d0 - RATIO_R ) * 1000.d0 

      ! Above density calculation is chemically unsound because it ignores
      ! chemical solvation.  
      ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of
      ! wet aerosol (kg/m3) 
      ! (bec, 6/17/10, jaegle 5/11/11)
      ! Redefine RATIO_R
      RATIO_R = A_RADI(K) / RWET

      ! Assume an initial density of 1000 kg/m3
      DEN  = 1000.D0
      DEN1 = 0.d0 !initialize (bec, 6/21/10)
      DO WHILE ( ABS( DEN1-DEN ) .gt. EPSI )
         ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
         WTP    = 100.d0 * A_DEN(K)/DEN * RATIO_R**3.d0
         ! Then calculate density of wet aerosol using equation 5 
         ! in Tang et al., 1997 [kg/m3]
         DEN1   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2) + 
     $               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.d0
         ! Now calculate new weight percent using above density calculation
         WTP    = 100.d0 * A_DEN(K)/DEN1 * RATIO_R**3
         ! Now recalculate new wet density [kg/m3]
         DEN   = ( 0.9971d0 + (A1 * WTP) + (A2 * WTP**2) + 
     $              (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.d0
      ENDDO

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

      ! This settling velocity is for the mid-point of the size bin.
      ! Need to integrate over the size bin, taking into account the
      ! mass distribution of sea-salt and the dependence of VTS on aerosol
      ! size. See WET_SETTLING in SEASALT_MOD.f for more details.
      ! (jaegle 5/11/11)
      SALT_MASS_TOTAL = 0d0
      VTS_WEIGHT      = 0d0
      ! Check what the min/max range of the SS size bins are
      IF ( RUM .le. Input_Opt%SALA_REDGE_um(2) ) THEN
        D0 = Input_Opt%SALA_REDGE_um(1)*2d0
        D1 = Input_Opt%SALA_REDGE_um(2)*2d0
      ELSE 
        D0 = Input_Opt%SALC_REDGE_um(1)*2d0
        D1 = Input_Opt%SALC_REDGE_um(2)*2d0
      ENDIF
     

      DO ID = 1, NR
      ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diamter):
      ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * DEN /Rw
        IF (DMID(ID) .ge. D0 .and. DMID(ID) .le. D1 ) THEN
           DMIDW = DMID(ID) * RWET/A_RADI(K)   ! wet radius [um]
           SALT_MASS   = SALT_V(ID) * RWET/A_RADI(K) * DEN / 
     &                    (DMIDW*0.5d0)
           VTS_WEIGHT  = VTS_WEIGHT + 
     &              SALT_MASS * VTS * (DMIDW/(RWET*1d6*2d0) )**2d0 *
     &                            (2d0 * DR *  RWET/A_RADI(K))
           SALT_MASS_TOTAL=SALT_MASS_TOTAL+SALT_MASS *
     &                            (2d0 * DR *  RWET/A_RADI(K))
        ENDIF

      ENDDO

      ! Final mass weighted setting velocity:
      VTS = VTS_WEIGHT/SALT_MASS_TOTAL


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

      ! Use the formulation of Slinn and Slinn (1980) for the impaction over
      ! water surfaces (jaegle 5/11/11)
      IF (LUC == 14) THEN
         EIM  = 10.d0**( -3.d0/ ST )   ! for water surfaces
      ELSE
         EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)
         EIM  = MIN( EIM, 0.6D0 )
      ENDIF

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.D0
      ELSE
         R1 = EXP( -1D0 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      ! Use the formulation of Slinn and Slinn (1980) for the impaction over
      ! water surfaces (jaegle 5/11/11)
      IF (LUC == 14) THEN
         RS   = 1.D0 / (USTAR**2.d0/ (W10*XCKMAN) * (EB + EIM ) + VTS)
      ELSE
         RS   = 1.D0 / (E0 * USTAR * (EB + EIM + EIN) * R1 )
      ENDIF

      END FUNCTION AERO_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_weightss
!
! !DESCRIPTION: Subroutine INIT\_WEIGHTSS calculates the volume size
!  distribution of sea-salt. This only has to be done once. We assume that
!  sea-salt is the combination of a coarse mode and accumulation model
!  log-normal distribution functions. The resulting arrays are: DMID = diameter
!  of bin and SALT\_V = dV/dln(D) [in um3]. (jaegle 5/11/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_WEIGHTSS( Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt
! 
! !REVISION HISTORY: 
!  11 May 2011 - L. Jaegle   - Initial version
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  14 Jun 2013 - R. Yantosca - Now accept Input_Opt via tha argument list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: N

      REAL*8  :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW ! jaegle, for weighting the settling vel.
      REAL*8  :: DEDGE
      INTEGER :: ID,NR

      ! increment of radius for integration of settling velocity (um)
      REAL*8, PARAMETER      :: DR    = 5.d-2

      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols
      ! geometric dry mean diameters (microns)
      REAL*8,  PARAMETER     ::   RG_A = 0.085d0
      REAL*8,  PARAMETER     ::   RG_C = 0.4d0
      ! sigma of the size distribution
      REAL*8,  PARAMETER     ::   SIG_A = 1.5d0
      REAL*8,  PARAMETER     ::   SIG_C = 1.8d0
      REAL*8, PARAMETER   :: PI =3.14159D0


      ! Number of bins between the lowest bound of of the accumulation mode
      ! sea salt and the upper bound of the coarse mode sea salt.
      NR =INT((( Input_Opt%SALC_REDGE_um(2) - 
     &           Input_Opt%SALA_REDGE_um(1) ) / DR )
     &                  + 0.5d0 )

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a
      ! coarse mode and accumulation model log-normal distribution functions
      !=================================================================

       ! Lower edge of 0th bin diameter [um]
       DEDGE=Input_Opt%SALA_REDGE_um(1) * 2d0

       ! Loop over diameters
       DO ID = 1, NR

           ! Diameter of mid-point in microns
           DMID(ID)  = DEDGE + ( DR )

           ! Calculate the dry volume size distribution as the sum of two
           ! log-normal size distributions. The parameters for the size
           ! distribution are based on Reid et al. and Quinn et al.
           ! The scaling factors 13. and 0.8 for acc and coarse mode aerosols
           ! are chosen to obtain a realistic distribution  
           ! SALT_V (D) = dV/dln(D) [um3]
           SALT_V(ID) = PI / 6d0* (DMID(ID)**3) * (
     &         13d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_A*2d0) )**2d0/
     &                   LOG(SIG_A)**2d0 )
     &         /( sqrt(2d0 * PI) * LOG(SIG_A) )  +
     &         0.8d0*exp(-0.5*( LOG(DMID(ID))-LOG(RG_C*2d0) )**2d0/
     &                   LOG(SIG_C)**2d0)
     &         /( sqrt(2d0 * PI) * LOG(SIG_C) )  )
           ! update the next edge
           DEDGE = DEDGE + DR*2d0
        ENDDO

      END SUBROUTINE INIT_WEIGHTSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dust_sfcrsi
!
! !DESCRIPTION: Function DUST\_SFCRSI computes the aerodynamic resistance of
!  dust aerosol tracers according to Seinfeld et al 96.  We do not consider
!  hygroscopic growth of the dust aerosol particles. (rjp, tdf, bmy, bec, 
!  4/1/04, 4/15/05)
!\\
!\\
! !INTERFACE:
!
      FUNCTION DUST_SFCRSI( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K       ! Drydep species (range: 1-NUMDEP)
      INTEGER, INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS   ! Pressure [kPa]
      REAL*8,  INTENT(IN) :: TEMP    ! Temperature [K]    
      REAL*8,  INTENT(IN) :: USTAR   ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
      REAL*8              :: RS      ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY: 
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
      
      END FUNCTION DUST_SFCRSI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: adust_sfcrsii
!
! !DESCRIPTION: Function ADUST\_SFCRSII computes the aerodynamic resistance of
!  non-size resolved aerosol according to Zhang et al 2001.  We do not consider
!  the hygroscopic growth of the aerosol particles. (rjp, tdf, bec, bmy, 
!  4/1/04, 4/15/05) 
!\\
!\\
!  This routine is used for all aerosols except dust, sulfate, and seasalt
!  (hotp 7/31/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION ADUST_SFCRSII( K, II, PRESS, TEMP, USTAR ) RESULT( RS )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K     ! Drydep tracer index (range: 1-NUMDEP)
      INTEGER, INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb = 100 Pa = 0.1 kPa)
      REAL*8,  INTENT(IN) :: TEMP  ! Temperature [K]    
      REAL*8,  INTENT(IN) :: USTAR ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
      REAL*8              :: RS    ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY:
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!  (3 ) Modified hotp for non size resolved aerosols. This is just DUST_SFCRSII
!        renamed and the diameter and density fixed. (hotp 7/12/07)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      END FUNCTION ADUST_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dust_sfcrsii
!
! !DESCRIPTION: Function DUST\_SFCRSII computes the aerodynamic resistance of
!  dust aerosol tracers according to Zhang et al 2001.  We do not consider the
!  hygroscopic growth of the aerosol particles. (rjp, tdf, bec, bmy, 4/1/04,
!  4/15/05)
!\\
!\\
! !INTERFACE:
!

      FUNCTION DUST_SFCRSII( K, II, PRESS, TEMP, USTAR, DIAM, DEN ) 
     & RESULT( RS )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K       ! Drydep tracer index (range: 1-NUMDEP)
      INTEGER, INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS   ! Pressure [kPa]
      REAL*8,  INTENT(IN) :: TEMP    ! Temperature [K]    
      REAL*8,  INTENT(IN) :: USTAR   ! Friction velocity [m/s]
      REAL*8,  INTENT(IN) :: DIAM    ! Particle diameter [m]
      REAL*8,  INTENT(IN) :: DEN     ! Particle density [kg/m3]
!
! !RETURN VALUE:
!
      REAL*8              :: RS      ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY: 
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  31 Jan 2014 - R. Yantosca - Now pass DIAM and DEN as arguments so as to
!                              avoid parallelization errors when using
!                              the TOMAS microphysics package.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

!------------------------------------------------------------------------------
! Prior to 1/31/14:
! Now pass these as arguments to avoid parallel errors w/ TOMAS (bmy, 1/31/14)
!      REAL*8  :: DIAM, DEN
!------------------------------------------------------------------------------
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
 
!-----------------------------------------------------------------------------  
! Prior to 1/31/14:
! Now pass these as arguments to avoid parallel errors w/ TOMAS (bmy, 1/31/14)
!      ! Particle diameter [m]
!      DIAM  = A_RADI(K) * 2.d0  
!
!      ! Particle density [kg/m3]
!      DEN   = A_DEN(K)          
!-----------------------------------------------------------------------------  

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

      END FUNCTION DUST_SFCRSII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_drydep
!
! !DESCRIPTION: Subroutine INIT\_DRYDEP initializes certain variables for the
!  GEOS-CHEM dry deposition subroutines. (bmy, 11/19/02, 10/19/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DRYDEP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
#if   defined( APM ) 
      USE APM_INIT_MOD, ONLY : NSO4, NSEA, NDSTB
#endif
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
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
!  (16) Add TOMAS aerosol NK1-NK30 and H2SO4 to drydep list (win, 7/14/09)
!  15 Feb 2011 - R. Yantosca - Added modifications for APM from G. Luo
!  15 Dec 2011 - M. Payer    - Update OVOC drydep according to Karl et al. 2010
!                              and add drydep for MVK and MACR. (J. Mao)
!  21 Dec 2011 - M. Payer    - Add allocation for size distribution of sea salt
!                              SALT_V and DMID (jaegle, 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO2 and O3 as part
!                              of removal of NOx-Ox partitioning
!  12 Jun 2013 - R. Yantosca - Bug fix: now only copy NUMDEP values to
!                              Input_Opt%NDVZIND and Input_Opt%DEPNAME
!  14 Jun 2013 - R. Yantosca - Now replace fields from tracer_mod.F
!                              with fields from Input_Opt
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  29 Aug 2013 - R. Yantosca - Assign XMW=118d-3 to RIP and IEPOX.  This now
!                              prevents XMW=0d0 from being passed to function
!                              DIFFG, where it is in the denominator.
!  04 Sep 2013 - R. Yantosca - Improve printout of drydep species
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_Hg
      LOGICAL :: IS_POPS
      INTEGER :: AS, N

      REAL*8, SAVE :: POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH
      REAL*8, SAVE :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL*8, SAVE :: POP_HSTAR, POP_DEL_H, POP_DEL_Hw
      REAL*8       :: DUM

      CHARACTER(LEN=255) :: MSG

      !=================================================================
      ! INIT_DRYDEP begins here!
      !=================================================================
      
      ! Assume success
      RC = GIGC_SUCCESS
     
      ! For GEOS-5, we need to do this only on the root CPU
      IF ( am_I_Root ) THEN
      DUM = 1.0

      ! Is this a mercury simulation?
      IS_Hg      = Input_Opt%ITS_A_MERCURY_SIM

      ! Is this a pops simulation?
      IS_POPS    = Input_Opt%ITS_A_POPS_SIM

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
      DO N = 1, Input_Opt%N_TRACERS

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

         ! NO2
         ! Replaced NOx with NO2  (mpayer, 3/14/13)
         ELSE IF ( N == IDTNO2 ) THEN
            NUMDEP          = NUMDEP + 1
            DRYDNO2         = NUMDEP
            NTRAIND(NUMDEP) = IDTNO2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NO2'
            HSTAR(NUMDEP)   = 0.01d0
            F0(NUMDEP)      = 0.1d0
            XMW(NUMDEP)     = 46d-3
            AIROSOL(NUMDEP) = .FALSE.
            
         ! O3
         ! Replaced Ox with O3  (mpayer, 3/14/13)
         ELSE IF ( N == IDTO3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTO3
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
            F0(NUMDEP)      = 1.0d0
            XMW(NUMDEP)     = 121d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! HNO3
         ELSE IF ( N == IDTHNO3 ) THEN
            NUMDEP          = NUMDEP + 1
            DRYDHNO3        = NUMDEP
            NTRAIND(NUMDEP) = IDTHNO3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'HNO3'
            HSTAR(NUMDEP)   = 1.0d+14
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

         ! ACET
         ELSE IF ( N == IDTACET ) THEN
            NUMDEP          = NUMDEP + 1
            NDVZIND(NUMDEP) = NUMDEP
            NTRAIND(NUMDEP) = IDTACET
            DEPNAME(NUMDEP) = 'ACET'
            HSTAR(NUMDEP)   = 1.0d+5
            F0(NUMDEP)      = 1.0d0
            XMW(NUMDEP)     = 58d-3
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
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
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
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
            XMW(NUMDEP)     = 58d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! MGLY 
         ELSE IF ( N == IDTMGLY ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMGLY
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MGLY'
            HSTAR(NUMDEP)   = 3.7d+3
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
            XMW(NUMDEP)     = 72d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! GLYC
         ELSE IF ( N == IDTGLYC ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTGLYC
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'GLYC'
            HSTAR(NUMDEP)   = 4.1d+4
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
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

         !--------------------------------------------
         ! jpp, 2/27/08, adding HOBr and HBr dry dep
         ! HOBr
         ELSE IF ( N == IDTHOBR ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTHOBr
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'HOBr'
            ! from R. Sander's compilation...
            ! jpp 4/27/09: This adjusted Henry's Law
            !              coefficient is calculated as the avg
            !              of the following values (as found in
            !              R. Sander's document).
            !  1. Blatchley et al. 1992 (they said > 1.9E3, so I make it 2E3
            !  2. Frenzel et al. 1996 (this is an estimate, 6.1E3)
            !  3. Mozurkewich et al. 1995 (thermo calculation, 1.8E0)
            !
            ! jpp 1/14/2011: changing Hstar from 2.8d3 to 6.1d3, consistent
            !                with p-TOMCAT. ref: Freznel et al. 1998
            HSTAR(NUMDEP)   = 6.1d+3
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 97d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! HBr
         ELSE IF ( N == IDTHBR ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTHBr
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'HBr'
            ! R. Sander's compilation sheet...
            ! Chameides [1992]
            !HSTAR(NUMDEP)   = 7.2d-1
            !HSTAR(NUMDEP)   = 7.1d+8
            !HSTAR(NUMDEP)   = 2.9d+7
            !
            ! jpp 4/27/09: This adjusted Henry's Law
            !           coefficient is calculated as
            !           the average of the adjusted
            !           quantities from the following sources
            !           (all listed in R. Sander's Henry's law document)
            !   1. Brimblecombe and Clegg 1989
            !   2. Wagman et al. 1982
            !   3. Chameides and Stelson 1992
            !HSTAR(NUMDEP)   = 2.97d+16
            !
            ! jpp, 1/14/2011: made Hstar consistent with
            !                 p-TOMCAT. Ref: Dean 1992.
            !                 Assuming a pH of 7 for the plant stomata.
            HSTAR(NUMDEP)   = 7.1d+15
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 81d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! BrNO3
         ELSE IF ( N == IDTBRNO3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTBrNO3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'BrNO3'
            ! see R. Sander's compilation...
            ! infininte Hstar
            HSTAR(NUMDEP)   = 1.d+20
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 142.d-3
            AIROSOL(NUMDEP) = .FALSE.

         ELSE IF ( N == IDTBR2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTBr2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'Br2'
            ! taken from Table 2 of Yang et al. 2005, assuming
            ! that the temperature = 298 K (jpp, 1/13/2011)
            HSTAR(NUMDEP)   = 0.76d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 160.d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! NOTE: ISOPN is composed of 2 species (ISOPND + ISOPNB)
         ELSE IF ( N == IDTISOPN ) THEN

            ! First species; ISOPND
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOPN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOPND'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 1d0  !Karl et al. 2010(science)
            XMW(NUMDEP)     = 147d-3
            AIROSOL(NUMDEP) = .FALSE.

            ! Second species: ISOPNB
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOPN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOPNB'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 1d0  !Karl et al. 2010(science)
            XMW(NUMDEP)     = 147d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! NOTE: MMN is composed of 2 species (MACRN + MVKN)
         ELSE IF ( N == IDTMMN ) THEN

            ! First species; MACRN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMMN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MACRN'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 1d0  !Karl et al. 2010(science)
            XMW(NUMDEP)     = 149d-3
            AIROSOL(NUMDEP) = .FALSE.
       
            ! Second species: MVKN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMMN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MVKN'
            HSTAR(NUMDEP)   = 17d3 !ITO 2007
            F0(NUMDEP)      = 1d0  !Karl et al. 2010(science) 
            XMW(NUMDEP)     = 149d-3
            AIROSOL(NUMDEP) = .FALSE.   

         !ANIT
         ELSE IF ( N == IDTPROPNN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPROPNN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'PROPNN'  !NITROOXYACETONE IN SANDER TABLE
            HSTAR(NUMDEP)   = 1d3       
            F0(NUMDEP)      = 1d0 !Karl et al. 2010(science)
            XMW(NUMDEP)     = 119d-3
            AIROSOL(NUMDEP) = .FALSE.   

         !RIP
         ELSE IF ( N == IDTRIP ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTRIP
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'RIP'
            HSTAR(NUMDEP)   = 1.70d6    ! US EPA 2011
            F0(NUMDEP)      = 1.0d0     ! Karl et al. 2010(science)
            XMW(NUMDEP)     = 118d-3
            AIROSOL(NUMDEP) = .FALSE.

         !IEPOX
         ELSE IF ( N == IDTIEPOX ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTIEPOX
            NDVZIND(NUMDEP) = NUMDEP    !USE H2O2
            DEPNAME(NUMDEP) = 'IEPOX'
            HSTAR(NUMDEP)   = 1.30d8
            F0(NUMDEP)      = 1.0d0     ! Karl et al. 2010(science)
            XMW(NUMDEP)     = 118d-3    ! (Eddingsaas et al., 2010)
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
            
         !MVK
         ELSE IF ( N == IDTMVK ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMVK
            NDVZIND(NUMDEP) = NUMDEP !USE H2O2
            DEPNAME(NUMDEP) = 'MVK'
            HSTAR(NUMDEP)   = 4.4d1 !from R. Sander
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
            XMW(NUMDEP)     = 70d-3
            AIROSOL(NUMDEP) = .FALSE.

         !MACR
         ELSE IF ( N == IDTMACR ) THEN
            NUMDEP          = NUMDEP + 1
            NDVZIND(NUMDEP) = NUMDEP
            NTRAIND(NUMDEP) = IDTMACR
            DEPNAME(NUMDEP) = 'MACR'
            HSTAR(NUMDEP)   = 6.5d0 !from R. Sander
            F0(NUMDEP)      = 1.0d0 !Assume reactive
            XMW(NUMDEP)     = 70d-3
            AIROSOL(NUMDEP) = .FALSE.

         ELSE IF ( N == IDTHAC ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTHAC
            NDVZIND(NUMDEP) = NUMDEP 
            DEPNAME(NUMDEP) = 'HAC'
            HSTAR(NUMDEP)   = 2.9d3 !from Ito et al. (2007, JGR)
            F0(NUMDEP)      = 1.0d0 !Karl et al. 2010(science)
            XMW(NUMDEP)     = 74d-3
            AIROSOL(NUMDEP) = .FALSE.

         !MACR
         ELSE IF ( N == IDTALD2 ) THEN
            NUMDEP          = NUMDEP + 1
            NDVZIND(NUMDEP) = NUMDEP
            NTRAIND(NUMDEP) = IDTALD2
            DEPNAME(NUMDEP) = 'ALD2'
            HSTAR(NUMDEP)   = 1.5d1 !from R. Sander
            F0(NUMDEP)      = 1.0d0 !Assume reactive
            XMW(NUMDEP)     = 44d-3
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
            A_RADI(NUMDEP)  = ( Input_Opt%SALC_REDGE_um(1) + 
     &                          Input_Opt%SALC_REDGE_um(2) ) * 0.5d-6
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
            A_RADI(NUMDEP)  = ( Input_Opt%SALC_REDGE_um(1) + 
     &                          Input_Opt%SALC_REDGE_um(2) ) * 0.5d-6
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

         ! SOAupdate (hotp 6/12/09)
         ! semivolpoa: add POA (hotp 2/17/09)
         ELSE IF ( N == IDTPOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'POA1'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTPOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'POA2'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOAupdate: MTPA, MTPO use H-law = 0.049 (hotp 5/22/10)
         ELSE IF ( N == IDTMTPA ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMTPA
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MTPA'
            HSTAR(NUMDEP)   = 0.049d0
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

         ! SOAupdate: now MTPO for new mtp (hotp 5/22/10)
         ELSE IF ( N == IDTMTPO ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTMTPO
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'MTPO'
            HSTAR(NUMDEP)   = 0.049d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 136d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOAupdate: mono+sesq terpene in new mtp (hotp 5/22/10)
         ELSE IF ( N == IDTTSOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOG1'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTTSOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOG2'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTTSOG3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOG3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOG3'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTTSOG0 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOG0
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOG0'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTTSOA0 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOA0
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOA0'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTTSOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOA1'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTTSOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOA2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTTSOA3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTTSOA3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'TSOA3'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOAupdate: isoprene in net mtp formulation (hotp 5/22/10)
         ELSE IF ( N == IDTISOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOG1'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTISOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOG2'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTISOG3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOG3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOG3'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTISOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOA1'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTISOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOA2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTISOA3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTISOA3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ISOA3'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOAupdate: lumped aromatic/IVOC aerosol (hotp 5/12/10)
         ! LUMPAROMIVOC
         ELSE IF ( N == IDTASOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOG1'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTASOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOG2'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTASOG3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOG3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOG3'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTASOAN ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOAN
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOAN'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTASOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOA1'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTASOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOA2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTASOA3 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTASOA3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'ASOA3'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 150d-3
            AIROSOL(NUMDEP) = .TRUE.

         ! SOAupdate: (hotp 6/12/09)
         ! semivolpoa2: POG (hotp 3/2/09)  
         ELSE IF ( N == IDTPOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'POG1'
            ! make POG hydrophobic (hotp 9/24/09)
            !HSTAR(NUMDEP)   = 1d5
            HSTAR(NUMDEP)   = 9.5d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTPOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTPOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'POG2'
            HSTAR(NUMDEP)   = 9.5d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOAupdate: (hotp 6/12/09)
         ! semivolpoa4opoa: OPOG (hotp 3/18/09)  
         ELSE IF ( N == IDTOPOG1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOPOG1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OPOG1'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .FALSE.
         ELSE IF ( N == IDTOPOG2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOPOG2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OPOG2'
            HSTAR(NUMDEP)   = 1d5
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .FALSE.

         ! SOAupdate: (hotp 6/12/09)
         ! semivolpoa4opoa: add OPOA (hotp 3/18/09)
         ELSE IF ( N == IDTOPOA1 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOPOA1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OPOA1'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.
         ELSE IF ( N == IDTOPOA2 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTOPOA2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'OPOA2'
            HSTAR(NUMDEP)   = 0d0
            F0(NUMDEP)      = 0d0
            XMW(NUMDEP)     = 12d-3
            AIROSOL(NUMDEP) = .TRUE.
         ! End SOAupdate

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
            A_RADI(NUMDEP)  = ( Input_Opt%SALA_REDGE_um(1) + 
     &                          Input_Opt%SALA_REDGE_um(2) ) * 0.5d-6
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
            A_RADI(NUMDEP)  = ( Input_Opt%SALC_REDGE_um(1) + 
     &                          Input_Opt%SALC_REDGE_um(2) ) * 0.5d-6
            A_DEN(NUMDEP)   = 2200.d0         
            AIROSOL(NUMDEP) = .TRUE. 

#if   defined( TOMAS )
         !----------------------------------
         ! Sulfuric acid gas tracer   !(win, 7/14/09)
         !----------------------------------
         ! H2SO4
         ELSE IF ( N == IDTH2SO4 ) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTH2SO4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'H2SO4'
            HSTAR(NUMDEP)   = 1.0d+5
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 98d-3
            AIROSOL(NUMDEP) = .FALSE. 
                     
         !----------------------------------
         ! Aerosol number tracers     !(win, 7/14/09)
         !----------------------------------

         ! NK1 (aerosol)
         ELSE IF  ( N == IDTNK1) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK1
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK1'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

         ! NK2 (aerosol)
         ELSE IF  ( N == IDTNK2) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK2
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK2'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

         ! NK3 (aerosol)
         ELSE IF  ( N == IDTNK3) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK3
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK3'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK4 (aerosol)
         ELSE IF  ( N == IDTNK4) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK4
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK4'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK5 (aerosol)
         ELSE IF  ( N == IDTNK5) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK5
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK5'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK6 (aerosol)
         ELSE IF  ( N == IDTNK6) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK6
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK6'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK7 (aerosol)
         ELSE IF  ( N == IDTNK7) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK7
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK7'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK8 (aerosol)
         ELSE IF  ( N == IDTNK8) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK8
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK8'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK9 (aerosol)
         ELSE IF  ( N == IDTNK9) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK9
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK9'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK10 (aerosol)
         ELSE IF  ( N == IDTNK10) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK10
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK10'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK11 (aerosol)
         ELSE IF  ( N == IDTNK11) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK11
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK11'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

         ! NK12 (aerosol)
         ELSE IF  ( N == IDTNK12) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK12
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK12'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

# if  !defined( TOMAS12 )
         ! NK13 (aerosol)
         ELSE IF  ( N == IDTNK13) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK13
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK13'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK14 (aerosol)
         ELSE IF  ( N == IDTNK14) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK14
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK14'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK15 (aerosol)
         ELSE IF  ( N == IDTNK15) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK15
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK15'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.
 
#  if !defined( TOMAS15 )
         ! NK16 (aerosol)
         ELSE IF  ( N == IDTNK16) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK16
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK16'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK17 (aerosol)
         ELSE IF  ( N == IDTNK17) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK17
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK17'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK18 (aerosol)
         ELSE IF  ( N == IDTNK18) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK18
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK18'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK19 (aerosol)
         ELSE IF  ( N == IDTNK19) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK19
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK19'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK20 (aerosol)
         ELSE IF  ( N == IDTNK20) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK20
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK20'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK21 (aerosol)
         ELSE IF  ( N == IDTNK21) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK21
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK21'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

         ! NK22 (aerosol)
         ELSE IF  ( N == IDTNK22) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK22
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK22'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE. 

         ! NK23 (aerosol)
         ELSE IF  ( N == IDTNK23) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK23
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK23'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK24 (aerosol)
         ELSE IF  ( N == IDTNK24) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK24
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK24'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK25 (aerosol)
         ELSE IF  ( N == IDTNK25) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK25
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK25'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK26 (aerosol)
         ELSE IF  ( N == IDTNK26) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK26
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK26'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK27 (aerosol)
         ELSE IF  ( N == IDTNK27) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK27
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK27'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK28 (aerosol)
         ELSE IF  ( N == IDTNK28) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK28
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK28'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK29 (aerosol)
         ELSE IF  ( N == IDTNK29) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK29
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK29'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK30 (aerosol)
         ELSE IF  ( N == IDTNK30) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK30
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK30'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.
#  endif
!  endif not 15
# endif
! endif not 12

# if  defined( TOMAS40 )
         ! NK31 (aerosol)
         ELSE IF  ( N == IDTNK31) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK31
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK31'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK32 (aerosol)
         ELSE IF  ( N == IDTNK32) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK32
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK32'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK33 (aerosol)
         ELSE IF  ( N == IDTNK33) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK33
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK33'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK34 (aerosol)
         ELSE IF  ( N == IDTNK34) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK34
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK34'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK35 (aerosol)
         ELSE IF  ( N == IDTNK35) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK35
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK35'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK36 (aerosol)
         ELSE IF  ( N == IDTNK36) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK36
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK36'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK37 (aerosol)
         ELSE IF  ( N == IDTNK37) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK37
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK37'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK38 (aerosol)
         ELSE IF  ( N == IDTNK38) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK38
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK38'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK39 (aerosol)
         ELSE IF  ( N == IDTNK39) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK39
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK39'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.

         ! NK40 (aerosol)
         ELSE IF  ( N == IDTNK40) THEN
            NUMDEP          = NUMDEP + 1
            NTRAIND(NUMDEP) = IDTNK40
            NDVZIND(NUMDEP) = NUMDEP
            DEPNAME(NUMDEP) = 'NK40'
            HSTAR(NUMDEP)   = 0.0d0
            F0(NUMDEP)      = 0.0d0
            XMW(NUMDEP)     = 1.d0
            AIROSOL(NUMDEP) = .TRUE.
# endif
#endif

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
               HSTAR(NUMDEP)   = 1.0d+14
               ! Change from 1d6 to 1d14 per HAmos 9/23/11
               ! implemented 10/19/11 eck
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
            ! POPs gas phase
            IF ( N == IDTPOPG) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = IDTPOPG
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'POPG'
               ! HSTAR is Henry's Law in mol/L/atm. 
               ! For PHENANTHRENE, log Kaw = -2.76  
               ! so unitless Kaw = 1.73*10^-3 and Kwa = 1/Kaw
               ! Divide by R (0.0821 atm/M/K) and T (298 K) and get HSTAR = 23.5 M/atm
               ! For PYRENE, log Kaw = -3.27
               ! Using the same conversion, HSTAR = 76.1 M/atm
               ! For BENZO[a]PYRENE, log Kaw = -4.51
               ! Using the same conversion, HSTAR = 1.32d3 M/atm
               ! All log Kaws from Ma et al., J Chem Eng Data 2010, 55:819 
               HSTAR(NUMDEP)   = GET_POP_HSTAR(DUM)
               ! Adding Koa (octanol-ar partition coefficient) for POPs to
               ! account for accumulation in leaf cuticles
               ! Needs to be in units of mol/liter/atm as with HSTAR
               ! Divide unitless Koa at 298 K by product of R (0.0821 atm/M/K)
               ! and T (298 K)
               ! For PHENANTHRENE, log Koa = 7.64  
               ! use same conversion as for HSTAR to get 1.78d6 M/atm
               ! For PYRENE, log Koa = 8.86
               ! use same conversion to get 2.96d7 M/atm
               ! For BENZO[a]PYRENE, log Koa = 11.48
               ! use same conversion to get 1.23d10 M/atm 
               ! All log Koas from Ma et al., J Chem Eng Data 2010, 55:819 
               KOA(NUMDEP)     = GET_POP_KOA(DUM) *0.0409
               F0(NUMDEP)      = 0.0d0
               ! Need to change molecular weight for different POPs
               ! For PHENANTHRENE, MW = 178d-3 (kg/mol)
               ! For PYRENE, MW = 202d-3 (kg/mol)
               ! for BENZO[a]PYRENE, MW = 252d-3 (kg/mol)
               XMW(NUMDEP)     = GET_POP_XMW(DUM)
               AIROSOL(NUMDEP) = .FALSE. 
            ENDIF

            ! POPs PARTICLES - OC 
            IF ( N == IDTPOPPOC ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = IDTPOPPOC
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'POPPOC'
               HSTAR(NUMDEP)   = 0.0d0
               ! Koa for particulate POPs is set to equivalent of Henry's Law
               ! so cuticular accumulation is not considered
               KOA(NUMDEP)     = 0.0d0
               F0(NUMDEP)      = 0.0d0
               ! Need to change molecular weight for different POPs
               ! For PHENANTHRENE, MW = 178d-3 (kg/mol)
               ! For PYRENE, MW = 202d-3 (kg/mol)
               ! for BENZO[a]PYRENE, MW = 252d-3 (kg/mol)
               XMW(NUMDEP)     = GET_POP_XMW(DUM)
               AIROSOL(NUMDEP) = .TRUE. 
            ENDIF

            ! POPs PARTICLES - BC 
            IF ( N == IDTPOPPBC ) THEN
               NUMDEP          = NUMDEP + 1
               NTRAIND(NUMDEP) = IDTPOPPBC
               NDVZIND(NUMDEP) = NUMDEP
               DEPNAME(NUMDEP) = 'POPPBC'
               HSTAR(NUMDEP)   = 0.0d0
               ! Koa for particulate POPs is set to equivalent of Henry's Law
               ! so cuticular accumulation is not considered
               KOA(NUMDEP)     = 0.0d0
               F0(NUMDEP)      = 0.0d0
               ! Need to change molecular weight for different POPs
               ! For PHENANTHRENE, MW = 178d-3 (kg/mol)
               ! For PYRENE, MW = 202d-3 (kg/mol)
               ! for BENZO[a]PYRENE, MW = 252d-3 (kg/mol)
               XMW(NUMDEP)     = GET_POP_XMW(DUM)
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
         ! add dry dep of Hg0
         DRYPOPG = 0
         DRYPOPP_OC = 0
         DRYPOPP_BC = 0
         
         ! If drydep is turned on ...
         IF ( LDRYD ) THEN
         
            ! Loop over drydep species
            DO N = 1, NUMDEP

               ! Locate by DEPNAME
               SELECT CASE ( TRIM( DEPNAME(N) ) )
                  CASE( 'POPG' )
                     DRYPOPG = N
                  CASE( 'POPPOC' )
                     DRYPOPP_OC = N
                  CASE( 'POPPBC' )
                     DRYPOPP_BC = N
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
      ELSE
         NUMDEP  = Input_Opt%NUMDEP
      ENDIF

      !=================================================================
      ! Allocate arrays
      ! add allocation for SALT_V and DMID (jaegle 5/11/11)
      !=================================================================
      ALLOCATE( DEPSAV( IIPAR, JJPAR, NUMDEP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEPSAV' )
      DEPSAV = 0d0

      ALLOCATE( SALT_V( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SALT_V' )
      SALT_V = 0d0

      ALLOCATE( DMID( NR_MAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DMID' )
      DMID = 0d0

      ALLOCATE( SHIPO3DEP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SHIPO3DEP' )
      SHIPO3DEP = 0d0

      IF ( am_I_Root ) THEN
         Input_Opt%NUMDEP            = NUMDEP
         Input_Opt%NDVZIND(1:NUMDEP) = NDVZIND(1:NUMDEP)
         Input_Opt%DEPNAME(1:NUMDEP) = DEPNAME(1:NUMDEP)
      ENDIF

#if   defined( APM )
      ALLOCATE( DEPSAVSIZE( IIPAR, JJPAR, NSO4+NSEA+NDSTB ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DEPSAVSIZE' )
      DEPSAVSIZE = 0d0
#endif

      !=================================================================
      ! Echo information to stdout
      !=================================================================
      IF ( am_I_Root ) THEN

         ! Line 1
         MSG = 'INIT_DRYDEP: List of dry deposition species:'
         WRITE( 6, '(/,a)' ) TRIM( MSG )

         ! Line 2
         MSG =  '  #   Name      Tracer DEPVEL Henry''s    React.'
     &       // '   Molec.   Aerosol?'
         WRITE( 6, '(/,a)'   ) TRIM( MSG )

         ! Line 3
         MSG =  '                Number Index  Law Const  Factor'
     &       // '   Weight   (T or F)'
         WRITE( 6, '(a)'   ) TRIM( MSG )

         ! Separator
         WRITE( 6, '(a)'   ) REPEAT( '-', 70 )

         ! Output
         DO N = 1, NUMDEP
            WRITE( 6, 100 ) N,          ADJUSTL( DEPNAME(N) ), 
     &                      NTRAIND(N), NDVZIND(N), 
     &                      HSTAR(N),   F0(N),      
     &                      XMW(N),     AIROSOL(N)

         ENDDO
 100     FORMAT( i3, 3x, a8, 2(3x,i3), 4x, es8.1, 2(3x,f6.3), 3x, L3 )
      ENDIF

      !=================================================================
      ! Get input information 
      !=================================================================

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NOTE: Because READ_DRYDEP_INPUTS reads info from a netCDF %%%
      !%%% file, we may have to broadcast these.  However, the file  %%%
      !%%% dimensions are not very great (10 or 20 indices each)     %%%
      !%%% (bmy, 12/11/12)                                           %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! Read drydep inputs from the netCDF file
      ! Save Olson indices in INDOLSON array, in order to avoid
      ! confusion w/ previously-assinged variable name IOLSON
      CALL READ_DRYDEP_INPUTS( am_I_Root,
     &                         DRYCOEFF, INDOLSON, IDEP,    
     &                         IWATER,   NWATER,   IZO,      
     &                         IDRYDEP,  IRI,      IRLU,     
     &                         IRAC,     IRGSS,    IRGSO, 
     &                         IRCLS,    IRCLO,    IVSMAX, 
     &                         Input_Opt%Data_Dir_1x1    )


      ! Calls INIT_WEIGHTSS to calculate the volume distribution of 
      ! sea salt aerosols (jaegle 5/11/11)
      CALL INIT_WEIGHTSS( Input_Opt )
      
      ! Need to save NUMDEP into Input_Opt
      Input_Opt%NUMDEP = NUMDEP

      END SUBROUTINE INIT_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_drydep
!
! !DESCRIPTION: Subroutine CLEANUP\_DRYDEP deallocates all module arrays.
!  (bmy, 2/27/03, 2/22/05)
!\\
!\\
! !INTERFACE:
!   
      SUBROUTINE CLEANUP_DRYDEP
! 
! !REVISION HISTORY: 
!  (1 ) Remove reference to PBLFRAC array; it's obsolete (bmy, 2/22/05)
!  (2 ) Added SALT_V and DMID (jaegle, 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DRYDEP begins here!
      !=================================================================
      IF ( ALLOCATED( DEPSAV   ) ) DEALLOCATE( DEPSAV   )
      IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )
      IF ( ALLOCATED( SHIPO3DEP) ) DEALLOCATE( SHIPO3DEP)

#if   defined( APM )
      IF ( ALLOCATED( DEPSAVSIZE  ) ) DEALLOCATE( DEPSAVSIZE  )
#endif

      END SUBROUTINE CLEANUP_DRYDEP
!EOC
!------------------------------------------------------------------------------
#if   defined( APM )

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: depvelbin
!
! !DESCRIPTION: Computes size resolved dry velocities for APM simulation
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEPVELBIN( NPTS, RADIAT, TEMP,  SUNCOS, USTAR,  
     &                      CZ1,  OBK,    CFRAC, ZH,     LSNOW,  
     &                      DVEL, ZO,     RHB,   State_Met      )
!
! !USES:
!
      USE APM_INIT_MOD,       ONLY : NSO4,NSEA,NDSTB
      USE APM_INIT_MOD,       ONLY : RDRY,RSALT,RDST,DENDST
      USE APM_DRIV_MOD,       ONLY : GFTOT3D,DENWET3D
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: NPTS            ! # of grid boxes = IIPAR*JJPAR
      REAL*8,  INTENT(IN)  :: RADIAT(MAXIJ)   ! Solar radiation [W/m2]
      REAL*8,  INTENT(IN)  :: TEMP  (MAXIJ)   ! Temperature [K]
      REAL*8,  INTENT(IN)  :: SUNCOS(MAXIJ)   ! Cosine of solar zenith angle
      REAL*8,  INTENT(IN)  :: USTAR (MAXIJ)   ! Friction velocity [m/s]
      REAL*8,  INTENT(IN)  :: CZ1   (MAXIJ)   ! Alt @ which Vd is computed [m]
      REAL*8,  INTENT(IN)  :: OBK   (MAXIJ)   ! Monin-Obhukov length [m]
      REAL*8,  INTENT(IN)  :: CFRAC (MAXIJ)   ! Surface cloud fraction
      REAL*8,  INTENT(IN)  :: ZH    (MAXIJ)   ! Roughness height [m]
      REAL*8,  INTENT(IN)  :: RHB   (MAXIJ)   ! Relative humidity [%]
      TYPE(MetState), INTENT(IN) :: State_Met ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: DVEL(MAXIJ,NSO4+NSEA+NDSTB)   ! Drydep vel. [m/s]
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: As of GEOS-Chem v9-01-03d (12 Jan 2012), the addtional     %%%
!  %%%  arguments PRESSU and W10 have been passed.  PRESSU replaces the  %%%
!  %%%  PRESS=1.5d5 parameter and W10 is needed for seasalt.  These      %%%
!  %%%  would need to be built into the DEPVELBIN routine by the APM     %%%
!  %%%  programming team.  We have left the code as-is for now.          %%%
!  %%%     -- Bob Yantosca (09 Apr 2012)                                 %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Need as landtype input for each grid square (I,J); see CMN_DEP_mod.F
!     IREG(I,J)       - # of landtypes in grid square
!     ILAND(I,J,LDT)  - Land type ID for element LDT =1, IREG(I,J)
!                        (could be from any source - mapped to deposition 
!                        surface ID in input unit 65)
!     IJUSE(I,J,LDT) - Fraction ((per mil) of gridbox area occupied by
!                       land type element LDT
!                                                                             .
!  Need as leaf area index; see CMN_DEP_mod.F:
!     XLAI(I,J,LDT)  - Leaf Area Index of land type element LDT
!
! Need as meteorological input for each grid square(I,J) (passed):
!     RADIAT(IJLOOP) - Solar radiation in W m-2
!     TEMP(IJLOOP)   - Surface air temperature in K
!     SUNCOS(IJLOOP) - Cosine of solar zenith angle
!     LSNOW(IJLOOP)  - Logical for snow and sea ice
!     RHB(IJLOOP)    - Relative humidity at the surface
!
! Need as input for each species K (passed):
!     F0(K)          - reactivity factor for oxidation of biological
!     substances
!     HSTAR(K)       - Henry's Law constant
!     XMW(K)         - Molecular weight (kg/mole) of species K
!                      (used to calculate molecular diffusivities)
!
! Some variables used in the subroutine (passed):
!     LRGERA(IJLOOP) T -> stable atmosphere; a high aerodynamic
!     resistance
!                        (RA=1.E4 m s-1) is imposed; else RA is
!                        calculated
!     USTAR(IJLOOP)  - Friction velocity (m s-1)
!     CZ1(IJLOOP)    - Altitude (m) at which deposition velocity is
!     computed
!     OBK(IJLOOP)    - Monin-Obukhov length (m): set to 1.E5 m under
!     neutral 
!                      conditions
!     CFRAC(IJLOOP)  - Fractional cloud cover
!     ZH(IJLOOP)     - Mixing depth (m)
!
! Some variables used in the subroutine:
!     MAXDEP         - the maximum number of species for which the dry 
!                      deposition calculation is done
!     ZO(LDT)        - Roughness height (m) for specific surface type
!     indexed 
!                      by LDT
!     RSURFC(K,LDT)  - Bulk surface resistance (s m-1) for species K to 
!                      surface LDT
!     C1X(K)         - Total resistance to deposition (s m-1) for
!     species K
!
! Returned:
!     DVEL(IJLOOP,K) - Deposition velocity (m s-1) of species K
!
! !REVISION HISTORY: 
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  09 Apr 2012 - R. Yantosca - Remove reference to CMN_VEL_mod.F
!  09 Apr 2012 - R. Yantosca - Now use INTENT(IN), INTENT(OUT) for all args
!  14 Dec 2012 - M. Payer    - Now get ILAND, IUSE, IREG from State_Met

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  RI(NTYPE),RLU(NTYPE),RAC(NTYPE),RGSS(NTYPE),
     1        RGSO(NTYPE),RCLS(NTYPE),RCLO(NTYPE),
     2        RSURFC(NSO4+NSEA+NDSTB,NTYPE)       
                                   
      REAL*8  C1X(NSO4+NSEA+NDSTB),VD(NSO4+NSEA+NDSTB),
     &        VK(NSO4+NSEA+NDSTB)

#if   defined( GCAP )
      ! For the GISS/GCAP model, ZO is a function of land type 
      ! and is of dimension NTYPE (swu, bmy, 5/25/05)
      REAL*8  ZO(NTYPE)           
#else
      ! For GEOS-CTM, ZO is now of size MAXIJ and is passed via 
      ! the argument list, since it is a DAO met field. (bmy, 11/10/99)
      REAL*8  ZO(MAXIJ)           
#endif

      LOGICAL LRGERA(MAXIJ)

      REAL*8  VDS
      REAL*8  CZ,C1,RT,XNU,RAD0,RIX,GFACT,GFACI
      REAL*8  RDC,RLUXX,RGSX,RCL,DTMP1,DTMP2,DTMP3,DTMP4
      REAL*8  CZH,CKUSTR,REYNO,CORR1,CORR2,Z0OBK
      REAL*8  RA,RB,DUMMY1,DUMMY2,DUMMY3,DUMMY4
      REAL*8  XMWH2O,DAIR,TEMPK,TEMPC
      INTEGER IOLSON,II,IW
      INTEGER I,J,K,IJLOOP,LDT
      REAL*8  RCLX,RIXX,BIOFIT
      REAL*8  PRESS
      DATA PRESS /1.5D5/
!
! Logical for snow and sea ice
!
      LOGICAL LSNOW(MAXIJ)
!***********************************************************************       
!
      DO K = 1,NSO4+NSEA+NDSTB 
         DO IJLOOP =1,NPTS
            DVEL(IJLOOP,K) = 0.0D0
         ENDDO
      ENDDO
!***********************************************************************
!*                                 
!*    Begin section for computing deposition velocities           
!*                                 
!*                                 
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
!$OMP+PRIVATE( C1X,    VK, I, J                               )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      
        IJLOOP = (J-1)*IIPAR+I

!** CZ is Altitude (m) at which deposition velocity is computed
         CZ = CZ1(IJLOOP)
!** TEMPK and TEMPC are surface air temperatures in K and in C
         TEMPK = TEMP(IJLOOP)
         TEMPC = TEMP(IJLOOP)-273.15D0

!* Initialize variables
         DO K = 1,NSO4+NSEA+NDSTB
            VD(K) = 0.0D0
            DO LDT = 1,NTYPE
               RSURFC(K,LDT) = 0.D0
            END DO
         END DO

!** Calculate the kinematic viscosity XNU (m2 s-1) of air
         C1 = TEMPK/273.15D0
         XNU = 0.151D0*(C1**1.77D0)*1.0D-04
         RT = 1000.0D0*EXP(-TEMPC-4.0D0)

!****************************************************************************

         ! Loop over the # of Olson land types in this grid box (I,J)
         DO 175 LDT = 1, State_Met%IREG(I,J)

            ! If the land type is not represented in grid 
            ! box  (I,J), then skip to the next land type
            IF ( State_Met%IUSE(I,J,LDT) == 0 ) GOTO 175

            ! Olson land type index + 1
            IOLSON = State_Met%ILAND(I,J,LDT)+1

            ! Dry deposition land type index
            II     = IDEP(IOLSON)
C
            IF(LSNOW(IJLOOP)) II=1

            RI(LDT) = DBLE(IRI(II))
            IF (RI(LDT)   .GE. 9999.D0) RI(LDT)   = 1.D12

            IF ( IRLU(II) >= 9999 .or.
     &           State_Met%XLAI(I,J,LDT) <= 0.D0 ) THEN
               RLU(LDT) = 1.D6
            ELSE
               RLU(LDT) = DBLE( IRLU(II) )
     &                  / State_Met%XLAI(I,J,LDT) + RT
            ENDIF

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
!***************************************************************************

            RAD0 = RADIAT(IJLOOP)
            RIX = RI(LDT)
            IF (RIX .GE. 9999.D0) GOTO 150
            GFACT = 100.0D0
            IF (TEMPC .GT. 0.D0 .AND. TEMPC .LT. 40.D0)
     *           GFACT = 400.D0/TEMPC/(40.0D0-TEMPC)
            GFACI = 100.D0

            IF ( RAD0 > 0.D0 .and. State_Met%XLAI(I,J,LDT) > 0.D0 ) THEN
               GFACI=1.D0 / BIOFIT( DRYCOEFF,
     &                              State_Met%XLAI(I,J,LDT),
     &                              SUNCOS(IJLOOP),
     &                              CFRAC(IJLOOP)  )
            ENDIF
            
            RIX = RIX*GFACT*GFACI
 150        CONTINUE
!*                                 
            RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))
!*
            DO K = 1, NSO4

               ! [Zhang et al., 2001]
               RSURFC(K,LDT) = 
     &         PWET_SFCRSII( II, PRESS*1D-3, TEMPK, USTAR(IJLOOP),
     &         DENWET3D(I,J,1,1), GFTOT3D(I,J,1,1), RDRY(K))  
            ENDDO
 170        CONTINUE

!luogan add new dust
            DO K = 1, NDSTB
               ! [Zhang et al., 2001]
               RSURFC((NSO4+NSEA+K),LDT) =
     &         PWET_SFCRSII( II, PRESS*1D-3, TEMPK, USTAR(IJLOOP),
     &         DENDST(K), 1.d0, RDST(K))

            ENDDO

!luogan add new seasalt
            DO K = 1, NSEA
               ! [Zhang et al., 2001]
               RSURFC((NSO4+K),LDT) =
     &         PWET_SFCRSII( II, PRESS*1D-3, TEMPK, USTAR(IJLOOP),
     &         DENWET3D(I,J,1,2), GFTOT3D(I,J,1,2), RSALT(K))
            ENDDO
 175        CONTINUE
!*
!*    Set max and min values for bulk surface resistances         
!*                                 
         DO 190 K = 1,NSO4+NSEA+NDSTB
            DO 180 LDT = 1, State_Met%IREG(I,J)
               IF ( State_Met%IUSE(I,J,LDT) == 0 ) GOTO 180
               RSURFC(K,LDT)= MAX(1.D0, MIN(RSURFC(K,LDT), 9999.D0))
 180        CONTINUE
 190     CONTINUE
!*                                 
!*    Loop through the different landuse types present in the grid
!square     
!*                                 
         DO 500 LDT = 1, State_Met%IREG(I,J)
            IF ( State_Met%IUSE(I,J,LDT) == 0 ) GOTO 500
            IOLSON = State_Met%ILAND(I,J,LDT)+1

#if   defined( GCAP )
            DO 200 IW=1,NWATER
               IF (IOLSON .NE. IWATER(IW)) GOTO 200
               ZO(LDT) = 1.4D-02*USTAR(IJLOOP)*USTAR(IJLOOP)/9.8D0   
     1                 + 1.1D-01*XNU/USTAR(IJLOOP)        
               GOTO 210
 200        CONTINUE
            ZO(LDT) = DBLE(IZO(IOLSON))*1.D-4
 210        CONTINUE
#endif

!*********************************************************************
            CKUSTR = XCKMAN*USTAR(IJLOOP)

            ! Define REYNO for GCAP or GEOS met fields (swu, bmy,
            ! 5/25/05)
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

            ! Define Z0OBK for GCAP or GEOS met fields (swu, bmy,
            ! 5/25/05)
#if   defined( GCAP )
            Z0OBK = ZO(LDT)/OBK(IJLOOP)
#else
            Z0OBK = ZO(IJLOOP)/OBK(IJLOOP) 
#endif

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
!
!...aerodynamically rough or smooth surface
            IF ( REYNO < 0.1d0 ) GOTO 220

!...aerodynamically rough surface.
!*                                 
            IF (CORR1.LE.0.0D0 .AND. Z0OBK .LT. -1.D0)THEN
!*... unstable condition; set RA to zero. (first implemented in V. 3.2)
               RA = 0.D0
            ELSEIF (CORR1.LE.0.0D0 .AND. Z0OBK .GE. -1.D0) THEN
!*... unstable conditions; compute Ra as described above.
               DUMMY1 = (1.D0 - 9D0*CORR1)**0.5D0
               DUMMY2 = (1.D0 - 9D0*Z0OBK)**0.5D0
               DUMMY3 = ABS((DUMMY1 - 1.D0)/(DUMMY1 + 1.D0))
               DUMMY4 = ABS((DUMMY2 - 1.D0)/(DUMMY2 + 1.D0))
               RA = 0.74D0* (1.D0/CKUSTR) * LOG(DUMMY3/DUMMY4)
                                 
            ELSEIF((CORR1.GT.0.0D0).AND.(.NOT.LRGERA(IJLOOP)))THEN
!*...moderately stable conditions (z/zMO <1); compute Ra as described
!above
               RA = (1D0/CKUSTR) * 
     &              (.74D0*LOG(CORR1/Z0OBK) + 4.7D0*(CORR1-Z0OBK))
            ELSEIF(LRGERA(IJLOOP)) THEN
!*... very stable conditions
               RA = 1.D+04
            ENDIF
!* check that RA is positive; if RA is negative (as occasionally
!* happened in version 3.1) send a warning message.
#if   defined( GCAP )
            ! Debug output for GISS/GCAP model (swu, bmy, 5/25/05)
            IF (RA .LT. 0.) THEN
               WRITE (6,1001) IJLOOP,RA,CZ,ZO(LDT),OBK(IJLOOP)
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
!* Get total resistance for deposition - loop over species.
            DO 215 K = 1,NSO4+NSEA+NDSTB
               C1X(K) = RA + RSURFC(K,LDT)
 215        CONTINUE
            GOTO 240
 220        CONTINUE 
!** ... aerodynamically smooth surface
!** BUG FIX -- suppress drydep over smooth surfaces by setting Ra to a
!large
!** value (1e4).  This prevents negative dry deposition velocities when
!u*
!** is very small (djj, bmy, 5/8/00)
            DO 230 K = 1,NSO4+NSEA+NDSTB
                  RA     = 1.0D4
                  C1X(K) = RA + RSURFC(K,LDT)
 230        CONTINUE

 240        CONTINUE
C*                                 
            DO 400 K = 1,NSO4+NSEA+NDSTB
               VK(K) = VD(K)
               VD(K) = VK(K) + .001D0 * DBLE( State_Met%IUSE(I,J,LDT) )
     &               / C1X(K)
 400        CONTINUE
 500     CONTINUE

!**********************************************************************
!** Load array DVELBIN
         DO 550 K=1,NSO4+NSEA+NDSTB
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
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE DEPVELBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: pwet_sfcrsii
!
! !DESCRIPTION: For APM 
!\\
!\\
! !INTERFACE:
!
      FUNCTION PWET_SFCRSII( II, PRESS, TEMP, USTAR,
     &                       AERODEN, HGROW, A_RADIDRY) RESULT(RS)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: II         ! Surface type index of GEOS-CHEM
      REAL*8,  INTENT(IN) :: PRESS      ! Pressure in Kpa 1 mb = 100 pa = 0.1kPa
      REAL*8,  INTENT(IN) :: TEMP       ! Temperature (K)    
      REAL*8,  INTENT(IN) :: USTAR      ! Friction velocity (m/s)
      REAL*8,  INTENT(IN) :: AERODEN    ! Relative humidity (fraction)
      REAL*8,  INTENT(IN) :: HGROW      ! Relative humidity (fraction)
      REAL*8,  INTENT(IN) :: A_RADIDRY  ! Sea salt radius (m)
!
! !RETURN VALUE
!
      REAL*8              :: RS    ! Surface resistance for particles [s/m]

! !REVISION HISTORY: 
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: N

      REAL*8, PARAMETER   :: G0 = 9.8D0
      REAL*8, PARAMETER   :: BETA  = 2.d0
      REAL*8, PARAMETER   :: BOLTZ = 1.381d-23  ! Boltzmann constant (J/K)
      REAL*8, PARAMETER   :: E0 = 3.d0
      REAL*8  :: AIRVS     ! kinematic viscosity of Air (m^2/s)
      REAL*8  :: DP        ! Diameter of aerosol [um]
      REAL*8  :: PDP       ! Press * Dp      
      REAL*8  :: CONST     ! Constant for settling velocity calculations
      REAL*8  :: SLIP      ! Slip correction factor
      REAL*8  :: VISC      ! Viscosity of air (Pa s)
      REAL*8  :: DIFF      ! Brownian Diffusion constant for particles (m2/s)
      REAL*8  :: SC, ST    ! Schmidt and Stokes number (nondim)

      REAL*8  :: DIAM, DEN, RATIO_R, RWET
      REAL*8  :: EB, EIM, EIN, R1, AA, VTS

      INTEGER :: LUCINDEX(11) = (/12, 4, 1, 7,10, 2, 9, 8,11,15,14/)
      INTEGER :: LUC

      REAL*8  :: 
     & ALPHA(15) = (/ 1.0d0,  0.6d0,   1.1d0,   0.8d0, 0.8d0,  
     &                1.2d0,  1.2d0,  50.0d0,  50.0d0, 1.3d0, 
     &                2.0d0, 50.0d0, 100.0d0, 100.0d0, 1.5d0  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL*8  ::
     & GAMMA(15) = (/ 0.56d0, 0.58d0, 0.56d0, 0.56d0, 0.56d0, 
     &                0.54d0, 0.54d0, 0.54d0, 0.54d0, 0.54d0, 
     &                0.54d0, 0.54d0, 0.50d0, 0.50d0, 0.56d0 /)

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

      ! Aerosol growth with relative humidity in radius [m] 
      RWET = A_RADIDRY * HGROW 

      ! Diameter of the wet aerosol [m]
      DIAM  = RWET * 2.d0  

      ! Density of the wet aerosol [kg/m3]
      DEN   = AERODEN *1000.d0 !luo change g/cm3 to kg/m3 

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.d6 
 
      ! Constant for settling velocity calculation       
      CONST = DEN * DIAM**2 * G0 / 18.d0

      ! Slip correction factor as function of (P*dp)
      PDP  = PRESS * DP
      SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP( -0.059d0 * PDP) ) / PDP

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
      END FUNCTION PWET_SFCRSII
#endif
      END MODULE DRYDEP_MOD

      
