!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: drydep_mod.F
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
      USE DAO_MOD                               ! Met field subroutines
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD                          ! Diag counters & flags
      USE DIAG_MOD,       ONLY : AD44           ! Diagnostic arrays
#endif
      USE ERROR_MOD                             ! Error handling routines
      USE GC_GRID_MOD,    ONLY : GET_AREA_CM2   ! Grid box surface areas [cm2]
      USE PBL_MIX_MOD                           ! Boundary layer quantities
#if   defined( TOMAS )
      USE TOMAS_MOD                             ! For TOMAS microphysics
#endif
      USE PhysConstants                         ! Physical constants
      USE PRECISION_MOD                         ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_DRYDEP     
      PUBLIC :: DO_DRYDEP
      PUBLIC :: INIT_DRYDEP
      PUBLIC :: INIT_WEIGHTSS
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: DEPNAME
      PUBLIC :: NUMDEP
      PUBLIC :: NTRAIND
      PUBLIC :: IDEP,   IRGSS,  IRAC, IRCLS
      PUBLIC :: IRGSO,  IRLU,   IRI,  IRCLO, DRYCOEFF
!      PUBLIC :: NDVZIND ! MSL -> For MPI broadcasting in GCHP
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
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  29 Jan 2014 - R. Yantosca - Now set MAXDEP=105 for all simulations.  For
!                              TOMAS we had MAXDEP=100; this is close enough.
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jul 2014 - R. Yantosca - Removed reference to commsoil_mod.F
!  12 Nov 2014 - M. Yannetti - Added PRECISION_MOD, changed REAL*8 to REAL(fp)
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE
!  03 Mar 2015 - C. Keller   - Disabled DRYFLX. Now done in mixing_mod.F90
!  15 Jun 2015 - R. Yantosca - Disabled DRYFLXRnPbBe, also done in mixing_mod
!  01 Oct 2015 - R. Yantosca - Remove DVZ_MINVAL function
!  01 Oct 2015 - R. Yantosca - Added PRIVATE flag ID_ACET
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!  13 Jul 2016 - R. Yantosca - Remove DRYHg0, DRYHg2, DRYHgP ID flags
!  13 Jul 2016 - R. Yantosca - Also make NDVZIND a local variable for now
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  10 Mar 2017 - M. Sulprizio- Add special treatment for ALD2, following what is
!                              done for ACET. Drydep over the oceans is handled
!                              via the HEMCO SeaFlux extension. 
!  16 Mar 2017 - R. Yantosca - Remove references to obsolete vars in Input_Opt
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  01 Dec 2018 - H.P. Lin    - Move DEPSAV to State_Chm%DryDepSav
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER,  PARAMETER :: NR_MAX    = 200       ! # of seasalt bins 
      INTEGER,  PARAMETER :: NDRYDTYPE = 11        ! # of drydep land types
!
! PRIVATE TYPES:
!
      !========================================================================
      !  MODULE VARIABLES:
      !
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
      !  NTRAIND   : Stores species numbers of drydep species
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
      !    NTYPE     : Max # of landtypes / grid box
      !    NPOLY     : Number of drydep polynomial coefficients
      !    NSURFTYPE : Number of Olson land types
      !
      !  NOTE: these grid-dependent variables are defined in State_Chm_Mod.F90
      !    DEPSAV    : Array containing dry deposition frequencies [s-1]
      !========================================================================

      ! Scalars
      INTEGER                        :: NUMDEP,   NWATER
      INTEGER                        :: DRYHg0,   DRYHg2,   DryHgP
      INTEGER                        :: id_ACET,  id_ALD2
      INTEGER                        :: id_NK1
      INTEGER                        :: id_HNO3,  id_PAN,   id_ISOPND

      ! Arrays for Baldocchi drydep polynomial coefficients
      REAL(fp), TARGET               :: DRYCOEFF(NPOLY    )

      ! Arrays that hold information for each of the 74 Olson land types
      INTEGER                        :: INDOLSON(NSURFTYPE )
      INTEGER                        :: IDEP    (NSURFTYPE )
      INTEGER                        :: IZO     (NSURFTYPE )
      INTEGER                        :: IWATER  (NSURFTYPE )

      ! Arrays that hold information for each of the 11 drydep land types
      INTEGER                        :: IDRYDEP (NDRYDTYPE)
      INTEGER                        :: IRAC    (NDRYDTYPE)
      INTEGER                        :: IRCLO   (NDRYDTYPE)
      INTEGER                        :: IRCLS   (NDRYDTYPE)
      INTEGER                        :: IRGSS   (NDRYDTYPE)
      INTEGER                        :: IRGSO   (NDRYDTYPE)
      INTEGER                        :: IRI     (NDRYDTYPE)
      INTEGER                        :: IRLU    (NDRYDTYPE)
      INTEGER                        :: IVSMAX  (NDRYDTYPE)

      ! Arrays that hold information about the dry-depositing species
      LOGICAL,           ALLOCATABLE :: AIROSOL (:    ) ! Is Aerosol? (T/F)
      INTEGER,           ALLOCATABLE :: NDVZIND (:    ) ! Drydep index
      INTEGER,           ALLOCATABLE :: FLAG    (:    ) ! Drydep scaling flag
      INTEGER,           ALLOCATABLE :: NTRAIND (:    ) ! Species index
      REAL(f8),          ALLOCATABLE :: HSTAR   (:    ) ! Henry's K0 [M/atm]
      REAL(f8),          ALLOCATABLE :: KOA     (:    ) ! POP's KOA
      REAL(f8),          ALLOCATABLE :: F0      (:    ) ! Reactivity factor [1]
      REAL(f8),          ALLOCATABLE :: XMW     (:    ) ! Mol wt. [kg/mol]
      REAL(f8),          ALLOCATABLE :: A_RADI  (:    ) ! Aer radius [m]
      REAL(f8),          ALLOCATABLE :: A_DEN   (:    ) ! Aer density [kg/m3]
      CHARACTER(LEN=14), ALLOCATABLE :: DEPNAME (:    ) ! Species name

      ! Allocatable arrays
      REAL(f8),          ALLOCATABLE :: DMID    (:    )
      REAL(f8),          ALLOCATABLE :: SALT_V  (:    )
      
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
      SUBROUTINE DO_DRYDEP( am_I_Root, Input_Opt,  State_Met,
     &                      State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE ErrCode_Mod
      USE HCO_ERROR_MOD
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
      USE State_Met_Mod,      ONLY : MetState
      USE Time_Mod,           ONLY : Get_Ts_Chem
      USE UnitConv_Mod
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
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
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
!  15 Jan 2015 - R. Yantosca - Add new netCDF diagnostics structure
!  12 Aug 2015 - E. Lundgren - Tracer units are now [kg/kg] 
!  01 Oct 2015 - R. Yantosca - Remove DVZ_MINVAL; now use species database
!                              to set the minimum drydep velocity for aerosols
!  01 Oct 2015 - R. Yantosca - Now use species database to set Vd for aerosols
!                              over snow/ice; combine with the prior update to
!                              enforce a minimum Vd for aerosol species
!  01 Oct 2015 - R. Yantosca - Remove references to tracerid_mod.F (except
!                              for TOMAS for now)
!  21 Jan 2016 - E. Lundgren - Update netcdf diagnostics for drydep trcr names
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  17 May 2016 - M. Sulprizio- Remove IJLOOP and change dimension of arrays
!                              from (MAXIJ) to (IIPAR,JJPAR)
!  18 Jul 2016 - M. Sulprizio- Remove special handling for ISOPN and MVKN for
!                              ND44 netCDF diagnostic. Family tracers have been
!                              eliminated.
!  13 Jul 2017 - M. Sulprizio- Add scaling of species drydep relative to HNO3,
!                              PAN, and ISOPN (from K. Travis, J. Fisher)
!  05 Oct 2017 - R. Yantosca - Now update State_Diag%DryDepVel diagnostic
!  27 Aug 2018 - E. Lundgren - Implement budget diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: I,   J,   L,   D,   N, NDVZ
      REAL(f8)           :: DVZ, THIK
      CHARACTER(LEN=255) :: ErrMsg,  ThisLoc


      ! Arrays
      LOGICAL       :: LSNOW (IIPAR,JJPAR       ) ! Flag for snow/ice on the sfc
      REAL(f8)      :: CZ1   (IIPAR,JJPAR       ) ! Midpt ht of 1st model lev[m]
      REAL(f8)      :: TC0   (IIPAR,JJPAR       ) ! Grid box sfc temperature [K]
      REAL(f8)      :: ZH    (IIPAR,JJPAR       ) ! PBL height [m]
      REAL(f8)      :: OBK   (IIPAR,JJPAR       ) ! Monin-Obhukov Length [m]
      REAL(f8)      :: CFRAC (IIPAR,JJPAR       ) ! Column cloud frac [unitless]
      REAL(f8)      :: RADIAT(IIPAR,JJPAR       ) ! Solar radiation [W/m2]
      REAL(f8)      :: USTAR (IIPAR,JJPAR       ) ! Grid box friction vel [m/s]
      REAL(f8)      :: RHB   (IIPAR,JJPAR       ) ! Relative humidity [unitless]
      REAL(f8)      :: DVEL  (IIPAR,JJPAR,NUMDEP) ! Drydep velocities [m/s]
      REAL(f8)      :: PRESSU(IIPAR,JJPAR       ) ! Local surface pressure [Pa]
      REAL(f8)      :: W10   (IIPAR,JJPAR       ) ! 10m windspeed [m/s]
      REAL(f8)      :: AZO   (IIPAR,JJPAR)        ! Z0, per (I,J) square
      REAL(f8)      :: SUNCOS_MID(IIPAR,JJPAR)    ! COS(SZA) @ midpoint of the
                                                  !  current chemistry timestep

      ! Pointers
      REAL(fp), POINTER :: DEPSAV (:,:,:   )      ! Dry deposition frequencies [s-1]

      ! For ESMF, need to assign these from Input_Opt
      LOGICAL       :: PBL_DRYDEP 
      LOGICAL       :: LPRT

      ! Objects
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! DO_DRYDEP begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Initialize
      SpcInfo => NULL()
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Do_DryDep  (in module GeosCore/drydep_mod.F)'    

      ! Point to columns of derived-type object fields
      DEPSAV     => State_Chm%DryDepSav  
      
      ! Copy values from the Input Options object to local variables
      PBL_DRYDEP = Input_Opt%PBL_DRYDEP
      LPRT       = Input_Opt%LPRT

      ! Call METERO to obtain meterological fields (all 1-D arrays)
      ! Added sfc pressure as PRESSU and 10m windspeed as W10 
      !  (jaegle 5/11/11, mpayer 1/10/12)
      CALL METERO( State_Met, CZ1,    TC0,     OBK,       CFRAC,  
     &             RADIAT,    AZO,    USTAR,   ZH,        LSNOW, 
     &             RHB,       PRESSU, W10,     SUNCOS_MID          )

      ! Call DEPVEL to compute dry deposition velocities [m/s]
      CALL DEPVEL( am_I_Root, Input_Opt, State_Met,  State_Chm,     
     &             State_Diag, RADIAT,   TC0,        SUNCOS_MID,
     &             F0,        HSTAR,     XMW,        AIROSOL,           
     &             USTAR,     CZ1,       OBK,        CFRAC,
     &             ZH,        LSNOW,     DVEL,       AZO,              
     &             RHB,       PRESSU,    W10,        RC            )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in call to "DEPVEL!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Compute dry deposition frequencies; archive diagnostics
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, THIK, D, N, NDVZ, DVZ, SpcInfo )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! THIK = thickness of surface layer [m]
         THIK   = State_Met%BXHEIGHT(I,J,1)

         ! Now we calculate drydep throughout the entire PBL.  
         ! Make sure that the PBL depth is greater than or equal 
         ! to the thickness of the 1st layer (rjp, bmy, 7/21/03)
         ! Add option for non-local PBL mixing scheme: THIK must
         ! be the first box height. (Lin, 03/31/09)
         ! Now use PBL_DRYDEP instead of LNLPBL (ckeller, 3/5/15). 
         IF (PBL_DRYDEP) THIK = MAX( ZH(I,J), THIK )

         ! Loop over drydep species
         DO D = 1, State_Chm%nDryDep

            ! GEOS-CHEM species number
            N       =  State_Chm%Map_DryDep(D)

            ! Get info about this species from the database
            SpcInfo => State_Chm%SpcData(N)%Info

            ! Index of drydep species in the DVEL array 
            ! as passed back from subroutine DEPVEL
            NDVZ    =  NDVZIND(D)

            ! Dry deposition velocity [cm/s]
            DVZ     =  DVEL(I,J,NDVZ) * 100.e+0_f8

            ! Scale relative to specified species (krt, 3/1/15)
            IF ( FLAG(D) .eq. 1 )  THEN

               ! Scale species to HNO3 
               DVZ = DVZ * sqrt(State_Chm%SpcData(id_HNO3)%Info%MW_g)
     &                   / sqrt(SpcInfo%MW_g)

            ELSE IF ( FLAG(D) .eq. 2 ) THEN

               ! Scale species to PAN
               DVZ = DVZ * sqrt(State_Chm%SpcData(id_PAN)%Info%MW_g)
     &                   / sqrt(SpcInfo%MW_g)         

            ELSE IF ( FLAG(D) .eq. 3 ) THEN

               ! Scale species to ISOPN
               DVZ = DVZ * sqrt(State_Chm%SpcData(id_ISOPND)%Info%MW_g)
     &                   / sqrt(SpcInfo%MW_g)            

            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for snow vs. ice
            !-----------------------------------------------------------
            IF ( LSNOW(I,J) ) THEN 
               
               !-------------------------------------
               ! %%% SURFACE IS SNOW OR ICE %%%
               !-------------------------------------
               IF ( SpcInfo%DD_DvzAerSnow > 0.0_fp ) THEN

                  ! For most aerosol species (basically everything
                  ! except sea salt and dust species), we just set
                  ! the deposition velocity over snow to a fixed value.
                  ! (Modification by Jenny Fisher, dated 8/1/11)
                  DVZ = DBLE( SpcInfo%DD_DvzAerSnow )

               ELSE
                  
                  ! Otherwise, enforce a minimum drydep velocity over snow
                  ! (cf. the GOCART model).  NOTE: In practice this will 
                  ! only apply to the species SO2, SO4, MSA, NH3, NH4, NIT.
                  DVZ = MAX( DVZ, DBLE( SpcInfo%DD_DvzMinVal(1) ) )

               ENDIF

            ELSE

               !-------------------------------------
               ! %%% SURFACE IS NOT SNOW OR ICE %%%
               !-------------------------------------

               ! Enforce a minimum drydep velocity over land (cf. the 
               ! GOCART model).  NOTE: In practice this will only apply 
               ! to the species SO2, SO4, MSA, NH3, NH4, NIT.
               DVZ = MAX( DVZ, DBLE( SpcInfo%DD_DvzMinVal(2) ) ) 

            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for ACETONE
            !-----------------------------------------------------------

            ! For ACET, we need to only do drydep over the land
            ! and not over the oceans.
            IF ( N == id_ACET ) THEN
               IF ( IS_LAND( I, J, State_Met ) ) THEN
                  DVZ = 0.1e+0_f8
               ELSE
                  DVZ = 0e+0_f8
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Special treatment for ALD2
            !-----------------------------------------------------------

            ! For ALD2, we need to only do drydep over the land
            ! and not over the oceans.
            IF ( N == id_ALD2 ) THEN
               IF ( .not. IS_LAND( I, J, State_Met ) ) THEN
                  DVZ = 0e+0_f8
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! Compute drydep frequency and update diagnostics
            !-----------------------------------------------------------

            ! Dry deposition frequency [1/s]
            DEPSAV(I,J,D) = ( DVZ / 100.e+0_f8 ) / THIK

#if defined( BPCH_DIAG )
            ! ND44 diagnostic
            IF ( ND44 > 0 ) THEN 
               ! Drydep velocity [cm/s]
               AD44(I,J,D,2) = AD44(I,J,D,2) + DVZ

!------------------------------------------------------------------------------
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!               ! Instantaneous drydep velocity for ND49 [cm/s]
!               AD44b(I,J,D)  = DVZ
!------------------------------------------------------------------------------
            ENDIF
#endif

            ! Archive dry dep velocity [cm/s]
            IF ( State_Diag%Archive_DryDepVel ) THEN
               State_Diag%DryDepVel(I,J,D) = DVZ
            ENDIF

            ! Free pointer
            SpcInfo => NULL()
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT .and. am_I_Root ) THEN
         CALL DEBUG_MSG( '### DO_DRYDEP: after dry dep' )
      ENDIF

      ! Nullify pointers
      NULLIFY( DEPSAV )

      END SUBROUTINE DO_DRYDEP
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
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      LOGICAL,  INTENT(OUT) :: LSNOW (IIPAR,JJPAR)  ! Flag for denoting snow/ice
      REAL(f8), INTENT(OUT) :: CZ1   (IIPAR,JJPAR)  ! Midpt ht of 1st model lev
                                                    !  [m]
      REAL(f8), INTENT(OUT) :: TC0   (IIPAR,JJPAR)  ! Grid box sfc temp [K]
      REAL(f8), INTENT(OUT) :: OBK   (IIPAR,JJPAR)  ! Monin-Obhukov length [m]
      REAL(f8), INTENT(OUT) :: CFRAC (IIPAR,JJPAR)  ! Column cloud fraction
                                                    !  [unitless]
      REAL(f8), INTENT(OUT) :: RADIAT(IIPAR,JJPAR)  ! Solar radiation @ ground
                                                    !  [W/m2]
      REAL(f8), INTENT(OUT) :: RHB   (IIPAR,JJPAR)  ! Rel humidity at sfc
                                                    !  [unitless]
      REAL(f8), INTENT(OUT) :: USTR  (IIPAR,JJPAR)  ! Friction velocity [m/s]
      REAL(f8), INTENT(OUT) :: ZH    (IIPAR,JJPAR)  ! PBL height [m]
      REAL(f8), INTENT(OUT) :: PRESSU(IIPAR,JJPAR)  ! Local surface press [Pa]
      REAL(f8), INTENT(OUT) :: W10   (IIPAR,JJPAR)  ! 10 meter windspeed [m/s]
      REAL(f8), INTENT(OUT) :: SUNCOS_MID(IIPAR,JJPAR) ! COS(SZA) @ midpt of
                                                       !  current chem timestep
      REAL(f8), INTENT(OUT) :: AZO(IIPAR,JJPAR) ! Roughness heights, by grid box
!
! !REMARKS:
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
!  25 Jul 2014 - R. Yantosca - Now remove reference to function SFCWINDSQR
!  17 May 2016 - M. Sulprizio- Remove IJLOOP and change dimension of output
!                              arrays from (MAXIJ) to (IIPAR,JJPAR)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: I,  J
      REAL(f8) :: THIK
      REAL(f8) :: SP
      REAL(f8) :: SFCWINDSQR

      !=================================================================
      ! METERO begins here!
      !=================================================================

      ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, THIK, SP, SFCWINDSQR )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! THIK = thickness of layer 1 [m]
         THIK        = State_Met%BXHEIGHT(I,J,1)

         ! Midpoint height of first model level [m]
         CZ1(I,J)    = THIK / 2.0e+0_f8

         ! Local surface pressure [hPa] (mpayer, 1/10/12)
         ! Use moist air pressure for mean free path  (ewl, 3/2/15)
         SP          = State_Met%PEDGE(I,J,1)

         ! Convert from hPa to Pa for SFCPRESS
         PRESSU(I,J) = SP * 1.e+2_f8

         !==============================================================
         ! Return meterological quantities as 1-D arrays for DEPVEL
         !==============================================================

         ! Roughness height [m]
         AZO(I,J)    = State_Met%Z0(I,J)

         ! Column cloud fraction [unitless]
         CFRAC(I,J)  = State_Met%CLDFRC(I,J)

         ! Set logical LSNOW if snow and sea ice (ALBEDO > 0.4)
         LSNOW(I,J)  = ( State_Met%ALBD(I,J) > 0.4 )

         ! Monin-Obhukov length [m]
         OBK(I,J)    = GET_OBK( I, J, State_Met )

         ! Solar insolation @ ground [W/m2]
         RADIAT(I,J) = State_Met%SWGDN(I,J)

         ! Surface temperature [K]
         TC0(I,J)    = State_Met%TS(I,J)

         ! Friction velocity [m/s]
         USTR(I,J)   = State_Met%USTAR(I,J)

         ! Mixed layer depth [m]
         ZH(I,J)     = GET_PBL_TOP_m( I, J )

         ! Relative humidity @ surface [unitless] (bec, bmy, 4/13/05)
         !RHB(I,J)    = MIN( 0.99e+0_f8, RH(I,J,1) * 1.d-2 ) 
         !  changed to 98% due to vapor pressure lowering above sea water
         !  (Lewis & Schwartz, 2004)
         !  jaegle (5/11/11)
         RHB(I,J)    = MIN( 0.98e+0_f8, State_Met%RH(I,J,1) *
     &                    1.e-2_f8 )

         ! 10m windspeed [m/s] (jaegle 5/11/11)
         SFCWINDSQR  = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2
         W10(I,J)    = SQRT( SFCWINDSQR )

         ! Cosine of solar zenith angle at midpoint 
         ! of the current chemistry timestep.
         SUNCOS_MID(I,J) = State_Met%SUNCOSmid(I,J)

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
! !IROUTINE: depvel
!
! !DESCRIPTION: Subroutine DEPVEL computes the dry deposition velocities using 
!  a resistance-in-series model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DEPVEL( am_I_Root, Input_Opt, State_Met, State_Chm,    
     &                   State_Diag, RADIAT,   TEMP,      SUNCOS,
     &                   F0,        HSTAR,     XMW,       AIROSOL,     
     &                   USTAR,     CZ1,       OBK,       CFRAC,           
     &                   ZH,        LSNOW,     DVEL,      ZO,            
     &                   RHB,       PRESSU,    W10,       RC          )
!
! !USES:
!
      USE Drydep_Toolbox_Mod, ONLY : BioFit
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root      ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt      ! Input Options object
      TYPE(MetState), INTENT(IN) :: State_Met      ! Meteorology state object
      TYPE(DgnState), INTENT(IN) :: State_Diag     ! Diagnostics state object

      REAL(f8), INTENT(IN) :: RADIAT (IIPAR,JJPAR) ! Solar radiation [W/m2]
      REAL(f8), INTENT(IN) :: TEMP   (IIPAR,JJPAR) ! Temperature [K]
      REAL(f8), INTENT(IN) :: SUNCOS (IIPAR,JJPAR) ! Cos of solar zenith angle
      LOGICAL,  INTENT(IN) :: AIROSOL(NUMDEP)      ! =T denotes aerosol species
      REAL(f8), INTENT(IN) :: F0     (NUMDEP)      ! React. factor for oxidation
                                                   !  of biological substances
      REAL(f8), INTENT(IN) :: HSTAR  (NUMDEP)      ! Henry's law constant
      REAL(f8), INTENT(IN) :: XMW    (NUMDEP)      ! Molecular weight [kg/mol]
      REAL(f8), INTENT(IN) :: USTAR  (IIPAR,JJPAR) ! Friction velocity [m/s]
      REAL(f8), INTENT(IN) :: CZ1    (IIPAR,JJPAR) ! Alt @ which Vd is computed
                                                   !  [m]
      REAL(f8), INTENT(IN) :: OBK    (IIPAR,JJPAR) ! Monin-Obhukov length [m]
      REAL(f8), INTENT(IN) :: CFRAC  (IIPAR,JJPAR) ! Surface cloud fraction
      REAL(f8), INTENT(IN) :: ZH     (IIPAR,JJPAR) ! Roughness height [m]
      REAL(f8), INTENT(IN) :: RHB    (IIPAR,JJPAR) ! Relative humidity [%]
      REAL(f8), INTENT(IN) :: PRESSU (IIPAR,JJPAR) ! Surface pressure [hPa]
      REAL(f8), INTENT(IN) :: W10    (IIPAR,JJPAR) ! Wind speed @ 10m altitude
                                                   !  [m/s]
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,  INTENT(OUT) :: RC                  ! Success or failure?
      REAL(f8), INTENT(OUT) :: DVEL(IIPAR,JJPAR,NUMDEP) ! Drydep velocity [m/s]
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
!     RADIAT(I,J)    - Solar radiation in W m-2
!     TEMP(I,J)      - Surface air temperature in K
!     SUNCOS(I,J)    - Cosine of solar zenith angle
!     LSNOW(I,J)     - Logical for snow and sea ice
!     RHB(I,J)       - Relative humidity at the surface
!     PRESSU(I,J)    - Local surface pressure
!     W10(I,J)       - 10m wind speed
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
!     LRGERA(I,J)    T -> stable atmosphere; a high aerodynamic resistance
!                        (RA=1.E4 m s-1) is imposed; else RA is calculated
!     USTAR(I,J)     - Friction velocity (m s-1)
!     CZ1(I,J)       - Altitude (m) at which deposition velocity is computed
!     OBK(I,J)       - Monin-Obukhov length (m): set to 1.E5 m under neutral 
!                      conditions
!     CFRAC(I,J)     - Fractional cloud cover
!     ZH(I,J)        - Mixing depth (m)
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
!     DVEL(I,J,K) - Deposition velocity (m s-1) of species K
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
!  06 Mar 2013 - H. Amos     - Merge C. Friedman's PAH code
!  31 May 2013 - R. Yantosca - Now pass State_Chm, for TOMAS
!  14 Jun 2013 - R. Yantosca - Now use Input_Opt%ITS_A_POPS_SIM
!  29 Aug 2013 - R. Yantosca - Bug fix: Skip to the next species if unless
!                              HSTAR>0 and XMW>0, or AIROSOL=t.  This avoids
!                              a floating-point invalid condition.
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  28 Jan 2014 - R. Yantosca - For TOMAS, don't hold A_RADI and A_DEN PRIVATE
!  19 May 2014 - C. Keller   - Now call BIOFIT from drydep_toolbox_mod.F\
!  12 Aug 2015 - E. Lundgren - Now accept am_I_Root and RC as input args
!  12 Aug 2015 - E. Lundgren - Now pass am_I_Root and RC to AERO_DIADEN 
!                              to enable unit conversion in that routine
!  22 Sep 2015 - R. Yantosca - Now use NUMDEP instead of MAXDEP for arrays
!  15 Mar 2016 - C. Keller   - Prevent very small numbers of CORR1 and Z0OBK
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  17 May 2016 - M. Sulprizio- Remove IJLOOP and change dimension of arrays
!                              from (MAXIJ) to (IIPAR,JJPAR); Also remove NPTS
!                              input argument
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
      REAL(f8)  RI(NTYPE),RLU(NTYPE),RAC(NTYPE),RGSS(NTYPE),
     1        RGSO(NTYPE),RCLS(NTYPE),RCLO(NTYPE),
     2        RSURFC(NUMDEP,NTYPE)       
                                   
      REAL(f8)  C1X(NUMDEP),VD(NUMDEP),VK(NUMDEP)              
      REAL(f8)  ZO(IIPAR,JJPAR)           

!      LOGICAL LDEP(MAXDEP)
      LOGICAL LDEP(NUMDEP)
      LOGICAL LRGERA(IIPAR,JJPAR)

      REAL(f8)  VDS
      REAL(f8)  CZ,C1,RT,XNU,RAD0,RIX,GFACT,GFACI
      REAL(f8)  RDC,RLUXX,RGSX,RCL,DTMP1,DTMP2,DTMP3,DTMP4
      REAL(f8)  CZH,CKUSTR,REYNO,CORR1,CORR2,Z0OBK
      REAL(f8)  RA,RB,DUMMY1,DUMMY2,DUMMY3,DUMMY4
      REAL(f8)  XMWH2O,DAIR,TEMPK,TEMPC
      REAL(f8)  RAKT,DUMMY2KT, DUMMY4KT, Z0OBKKT ! for corrected O3, krt,11/2017
      INTEGER IOLSON,II,IW
      INTEGER K,LDT
      REAL(f8)  RCLX,RIXX    !,BIOFIT

#if   defined( TOMAS )
      ! For TOMAS aerosol (win, 7/15/09)
      INTEGER  BIN
      REAL(f8)   SIZ_DIA(IIPAR,JJPAR, IBINS), SIZ_DEN(IIPAR,JJPAR,IBINS)
#endif

      ! Logical for snow and sea ice
      LOGICAL LSNOW(IIPAR,JJPAR)

      ! Loop indices (bmy, 3/29/12)
      INTEGER :: I, J

      ! Size of DRYCOEFF (ckeller, 05/19/14)
      INTEGER :: NN

      ! Added these to pass particle density & number to DUST_SFCRSII
      ! so as to avoid parallelization errors with TOMAS (bmy, 1/31/14)
      REAL(f8)  :: DIAM, DEN

      ! Logical switch for POPS in order to use octanol-air partitioning instead
      ! of Henry's Law for scaling of cuticular resitances (clf, 1/3/2011)
      LOGICAL :: IS_POPS

      ! Pointers to fields in State_Met
      INTEGER,  POINTER :: IREG(:,:)
      INTEGER,  POINTER :: ILAND(:,:,:)
      INTEGER,  POINTER :: IUSE(:,:,:)
      REAL(fp), POINTER :: XLAI(:,:,:)

      ! For making sure that all inputs to BIOFIT are of the same type
      REAL(fp)          :: XLAI_FP
      REAL(fp)          :: SUNCOS_FP
      REAL(fp)          :: CFRAC_FP

#if defined( MODEL_GEOS )
      ! Archive Ra?
      REAL(fp)          :: RA2M,  Z0OBK2M
      REAL(fp)          :: RA10M, Z0OBK10M
      !LOGICAL           :: WriteRa2m
      !LOGICAL           :: WriteRa10m
#endif

      ! For the species database
      INTEGER                :: SpcId
      TYPE(Species), POINTER :: SpcInfo

      ! Strings
      CHARACTER(LEN=255)  :: ThisLoc
      CHARACTER(LEN=512)  :: ErrMsg
!
! !DEFINED PARAMETERS:
!      
      ! Small number
      REAL(fp), PARAMETER :: SMALL = 1.0e-10_f8

      !=================================================================
      ! DEPVEL begins here!
      !=================================================================

      ! Assume success
      RC      =  GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at DEPVEL (in module GeosCore/drydep_mod.F)'

      ! Error check: PBL_DRYDEP must not be used together with non-local
      ! PBL mixing. 
      IF ( Input_Opt%LNLPBL .AND. Input_Opt%PBL_DRYDEP ) THEN 
         ErrMsg = 'PBL_DRYDEP must be disabled if LNLPBL=T!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Initialize pointers
      IREG    => State_Met%IREG
      ILAND   => State_Met%ILAND
      IUSE    => State_Met%IUSE
      XLAI    => State_Met%XLAI
      SpcInfo => NULL()

      ! Is this a POPs simmulation?
      IS_POPS = Input_Opt%ITS_A_POPS_SIM  ! clf, 1/3/2011

      ! Size of drycoeff (ckeller, 05/19/14)
      NN = SIZE(DRYCOEFF)

#if defined( MODEL_GEOS )
      ! Logical flag for Ra (ckeller, 12/29/17)
      State_Chm%DryDepRa2m  = 0.0_fp
      State_Chm%DryDepRa10m = 0.0_fp
      !WriteRa2m = ASSOCIATED ( State_Diag%DryDepRa2m )
      !IF ( WriteRa2m ) THEN
      !   State_Diag%DryDepRa2m = 0.0_fp
      !ENDIF
      !WriteRa10m = ASSOCIATED ( State_Diag%DryDepRa10m )
      !IF ( WriteRa10m ) THEN
      !   State_Diag%DryDepRa10m = 0.0_fp
      !ENDIF
#endif

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
         IF ( ( HSTAR(K) > 0e+0_f8   .and. 
     &          XMW  (K) > 0e+0_f8 ) .or.  AIROSOL(K) ) THEN
            LDEP(K) = .TRUE.
         ELSE
            LDEP(K) = .FALSE.
         ENDIF
      ENDDO

      ! Initialize DVEL
      DVEL = 0.0e+0_f8

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
      !
      ! Now pass am_I_Root and RC to enable species unit conversion to 
      ! [kg] in TOMAS code (ewl, 8/12/15)
      !=================================================================

      IF ( id_NK1 > 0 ) THEN
         CALL AERO_DIADEN( am_I_Root, 1,       Input_Opt, State_Met,
     &                     State_Chm, SIZ_DIA, SIZ_DEN,   RC        )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in call to "AERO_DIADEN"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF
#endif

      ! Add parallel DO-loop (bmy, 2/22/05)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( CZ,     TEMPK,  TEMPC,  K,      VD             )
!$OMP+PRIVATE( LDT,    RSURFC, C1,     XNU,    RT,     IOLSON )
!$OMP+PRIVATE( II,     RI,     RLU,    RAC,    RGSS,   RGSO   )
!$OMP+PRIVATE( RCLS,   RCLO,   RAD0,   RIX,    GFACT,  GFACI  )
!$OMP+PRIVATE( RDC,    XMWH2O, RIXX,   RLUXX,  RGSX,   RCLX   )
!$OMP+PRIVATE( DTMP1,  DTMP2,  DTMP3,  DTMP4,  VDS,    CZH    )
!$OMP+PRIVATE( CKUSTR, REYNO,  CORR1,  CORR2,  Z0OBK,  RA     )
!$OMP+PRIVATE( Z0OBKKT,RAKT,   DUMMY2KT,       DUMMY4KT       )
!$OMP+PRIVATE( DUMMY1, DUMMY2, DUMMY3, DUMMY4, DAIR,   RB     )
!$OMP+PRIVATE( C1X,    VK,     I,      J,      IW             )
!$OMP+PRIVATE( DIAM,   DEN,    XLAI_FP, SUNCOS_FP, CFRAC_FP   )
#if defined( MODEL_GEOS )
!$OMP+PRIVATE( RA2M,   Z0OBK2M, RA10M, Z0OBK10M               )
#endif
#if defined( TOMAS )
!$OMP+PRIVATE( BIN                                            )
#endif
!$OMP+PRIVATE( SpcId,  SpcInfo                                )

      DO J = 1, JJPAR
      DO I = 1, IIPAR

!** CZ is Altitude (m) at which deposition velocity is computed
         CZ = CZ1(I,J)
!** TEMPK and TEMPC are surface air temperatures in K and in C
         TEMPK = TEMP(I,J)
         TEMPC = TEMP(I,J)-273.15e+0_f8
!* Initialize variables
         DO K = 1,NUMDEP
            VD(K) = 0.0e+0_f8
            DO LDT = 1,NTYPE
               RSURFC(K,LDT) = 0.e+0_f8
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
         C1 = TEMPK/273.15e+0_f8
         XNU = 0.151e+0_f8*(C1**1.77e+0_f8)*1.0e-04_f8

!* Compute bulk surface resistance for gases.    
!*                                 
!* Adjust external surface resistances for temperature; 
!* from Wesely [1989], expression given in text on p. 1296.        
!*                                 
!* BUG FIX!  Wesely [1989] gives RT = 1000.0*EXP(-TEMPC-4.0)
!* so the inner parentheses are not needed (bmy, 3/4/99)
!*        RT = 1000.0*EXP(-(TEMPC-4.0))
         RT = 1000.0e+0_f8*EXP(-TEMPC-4.0e+0_f8)
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
            IF(LSNOW(I,J)) II=1

!* Read the internal resistance RI (minimum stomatal resistance for water 
!* vapor,per unit area of leaf) from the IRI array; a '9999' value means no 
!* deposition to stomata so we impose a very large value for RI.

            RI(LDT) = DBLE(IRI(II))
            IF (RI(LDT)   .GE. 9999.e+0_f8) RI(LDT)   = 1.e+12_f8

!** Cuticular resistances IRLU read in from 'drydep.table'
!** are per unit area of leaf;
!** divide them by the leaf area index to get a cuticular resistance for the
!** bulk canopy.  If IRLU is '9999' it means there are no cuticular
!** surfaces on which to deposit so we impose a very large value for RLU.
            IF ( IRLU(II) >= 9999 .or. XLAI(I,J,LDT) <= 0.e+0_f8 ) THEN
               RLU(LDT) = 1.e+6_f8
            ELSE
               RLU(LDT) = DBLE( IRLU(II) ) / XLAI(I,J,LDT) + RT
            ENDIF
!** The following are the remaining resistances for the Wesely
!** resistance-in-series model for a surface canopy
!** (see Atmos. Environ. paper, Fig.1).  
            RAC(LDT)  = MAX(DBLE(IRAC(II)), 1.e+0_f8)
            IF (RAC(LDT)  .GE. 9999.e+0_f8) RAC(LDT)  = 1.e+12_f8
            RGSS(LDT) = MAX(DBLE(IRGSS(II)) + RT ,1.e+0_f8)
            IF (RGSS(LDT) .GE. 9999.e+0_f8) RGSS(LDT) = 1.e12_f8
            RGSO(LDT) = MAX(DBLE(IRGSO(II)) + RT ,1.e+0_f8) 
            IF (RGSO(LDT) .GE. 9999.e+0_f8) RGSO(LDT) = 1.e+12_f8
            RCLS(LDT) = DBLE(IRCLS(II)) + RT           
            IF (RCLS(LDT) .GE. 9999.e+0_f8) RCLS(LDT) = 1.e+12_f8
            RCLO(LDT) = DBLE(IRCLO(II)) + RT          
            IF (RCLO(LDT) .GE. 9999.e+0_f8) RCLO(LDT) = 1.e+12_f8
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

            RAD0 = RADIAT(I,J)
            RIX = RI(LDT)
            IF (RIX .GE. 9999.e+0_f8) GO TO 150
            GFACT = 100.0e+0_f8
            IF (TEMPC .GT. 0.e+0_f8 .AND. TEMPC .LT. 40.e+0_f8)
     *           GFACT = 400.e+0_f8/TEMPC/(40.0e+0_f8-TEMPC)
            GFACI = 100.e+0_f8
            IF ( RAD0 > 0.e+0_f8 .and. XLAI(I,J,LDT) > 0.e+0_f8 ) THEN
               ! Now make sure all inputs to BIOFIT are flexible precision
               ! so that the code will compile properly (bmy, 12/18/14)
               XLAI_FP   = XLAI(I,J,LDT)
               SUNCOS_FP = SUNCOS(I,J)
               CFRAC_FP  = CFRAC(I,J)
               GFACI     = 1.e+0_f8 / BIOFIT( DRYCOEFF,  XLAI_FP, 
     &                                        SUNCOS_FP, CFRAC_FP, NN )
            ENDIF
            
            RIX = RIX*GFACT*GFACI
 150        CONTINUE
!*                                 
!*    Compute aerodynamic resistance to lower elements in lower part           
!*    of the canopy or structure, assuming level terrain - 
!*    equation (5) of Wesely [1989].
!*                                 
            RDC = 100.e+0_f8*(1.0e+0_f8+1000.0e+0_f8/(RAD0+10.e+0_f8))
!*
!*    Loop over species; species-dependent corrections to resistances
!*    are from equations (6)-(9) of Wesely [1989].
!*
            DO 160  K = 1,NUMDEP

!** exit for non-depositing species or aerosols.
               IF (.NOT. LDEP(K) .OR. AIROSOL(K)) GOTO 155

               !XMWH2O = 18.e-3_f8 ! Use global H2OMW (ewl, 1/6/16)
               XMWH2O = H2OMW * 1.e-3_f8
               RIXX = RIX*DIFFG(TEMPK,PRESSU(I,J),XMWH2O)/
     C              DIFFG(TEMPK,PRESSU(I,J),XMW(K))
     C              + 1.e+0_f8/(HSTAR(K)/3000.e+0_f8+100.e+0_f8*F0(K))
               RLUXX = 1.e+12_f8
               IF (RLU(LDT).LT.9999.e+0_f8)
     C              RLUXX = RLU(LDT)/(HSTAR(K)/1.0e+05_f8 + F0(K))

        ! If a POPs simulation, scale cuticular resistances with octanol-air
        ! partition coefficient (Koa) instead of HSTAR (clf, 1/3/2011)

               IF (IS_POPS)
     C              RLUXX = RLU(LDT)/(KOA(K)/1.0e+05_f8 + F0(K))

!*
!* To prevent virtually zero resistance to species with huge HSTAR, such
!* as HNO3, a minimum value of RLUXX needs to be set. The rationality
!* of the existence of such a minimum is demonstrated by the observed
!* relationship between Vd(NOy-NOx) and Ustar in Munger et al.[1996];
!* Vd(HNO3) never exceeds 2 cm s-1 in observations. The
!* corresponding minimum resistance is 50 s m-1.  This correction
!* was introduced by J.Y. Liang on 7/9/95.
!*

               RGSX = 1.e+0_f8/(HSTAR(K)/1.0e+05_f8/RGSS(LDT) + 
     1              F0(K)/RGSO(LDT))
               RCLX = 1.e+0_f8/(HSTAR(K)/1.0e+05_f8/RCLS(LDT) + 
     1              F0(K)/RCLO(LDT))
!*
!** Get the bulk surface resistance of the canopy, RSURFC, from the network
!** of resistances in parallel and in series (Fig. 1 of Wesely [1989])
               DTMP1=1.e+0_f8/RIXX
               DTMP2=1.e+0_f8/RLUXX
               DTMP3=1.e+0_f8/(RAC(LDT)+RGSX)
               DTMP4=1.e+0_f8/(RDC+RCLX)
               RSURFC(K,LDT) = 1.e+0_f8/(DTMP1 + DTMP2 + DTMP3 + DTMP4)
!** get surface deposition velocity for aerosols if needed;
!** equations (15)-(17) of Walcek et al. [1986]
 155           IF (.NOT. AIROSOL(K)) GOTO 160

!-----------------------------------------------------------------------------
! Prior to 9/30/15:
! Use fields from the species database instead of string tests (bmy, 9/30/15) 
! Leave the TOMAS code here for reference
!
!#if   defined( TOMAS )
!               ! (win, 7/15/09)
!               ELSE IF ( ( NTRAIND(K) >= id_NK1  )   .AND.
!     &                   ( NTRAIND(K) <  id_NK1+IBINS ) ) THEN
!
!                  !=====================================================
!                  ! Use size-resolved dry deposition calculations for 
!                  ! size-resolved aerosols.  Choose to use DUST_SFCRSII
!                  ! because these aerosol have diameter and density 
!                  ! calculation that takes hygroscopic growth into 
!                  ! account already, thus no need for such calculation
!                  ! like in AERO_SFCRSII.  (win, 7/19/07)
!                  !=====================================================     
!
!                  ! Get current aerosol diameter and density
!                  ! NOTE: In TOMAS the aerosol diameter and density
!                  ! evolves with time.  We have to save these in the 
!                  ! DIAM and DEN variables so that we can hold these
!                  ! PRIVATE for the parallel loop.
!                  BIN  = NTRAIND(K) - id_NK1 + 1
!                  DIAM = SIZ_DIA( I, J, BIN )     ! Diameter [m]
!                  DEN  = SIZ_DEN( I, J, BIN )     ! Density  [kg/m3]
!
!                  ! [Zhang et al., 2001]
!                  ! Now pass DIAM, DEN as arguments (bmy, 1/31/14)
!                  RSURFC(K,LDT) = DUST_SFCRSII( K, II, 
!     &                                     PRESSU(I,J)*1e-3_f8,
!     &                                          TEMPK, 
!     &                                          USTAR(I,J),   
!     &                                          DIAM, 
!     &                                          DEN )
!#endif
!-----------------------------------------------------------------------------

               ! Get information about this species from the database
               SpcId   =  NTRAIND(K)
               SpcInfo => State_Chm%SpcData(SpcId)%Info

               IF ( SpcInfo%DD_AeroDryDep ) THEN

                  !=====================================================
                  ! Use size-resolved dry deposition calculations for 
                  ! seasalt aerosols.  We need to account for the
                  ! hygroscopic growth of the aerosol particles.
                  !=====================================================

                  ! [Zhang et al., 2001]
                  RSURFC(K,LDT) = AERO_SFCRSII( K,     
     &                                          II,   
     &                                          PRESSU(I,J)*1e-3_f8, 
     &                                          TEMPK, 
     &                                          USTAR(I,J), 
     &                                          RHB(I,J),
     &                                          W10(I,J), 
     &                                          Input_Opt )

               ELSE IF ( SpcInfo%DD_DustDryDep ) THEN

                  !=====================================================
                  ! Use size-resolved dry deposition calculations for 
                  ! dust aerosols only.  Do not account for hygroscopic
                  ! growth of the dust aerosol particles.
                  !=====================================================  
#if defined( TOMAS ) 

                  !-------------------------------
                  !%%% TOMAS SIMULATIONS %%%
                  !-------------------------------

                  ! Get current aerosol diameter and density
                  ! NOTE: In TOMAS the aerosol diameter and density
                  ! evolves with time.  We have to save these in the 
                  ! DIAM and DEN variables so that we can hold these
                  ! PRIVATE for the parallel loop.
                  BIN  = NTRAIND(K) - id_NK1 + 1
                  DIAM = SIZ_DIA( I, J, BIN )     ! Diameter [m]
                  DEN  = SIZ_DEN( I, J, BIN )     ! Density  [kg/m3]

#else
                  !-------------------------------
                  !%%% NON-TOMAS SIMULATIONS %%%
                  !-------------------------------

                  ! Particle diameter, convert [m] -> [um]
                  DIAM  = A_RADI(K) * 2.e+0_f8 

                  ! Particle density [kg/m3]
                  DEN   = A_DEN(K)             

#endif

                  ! [Zhang et al., 2001]
                  RSURFC(K,LDT) = DUST_SFCRSII( K, 
     &                                          II,
     &                                          PRESSU(I,J)*1e-3_f8,
     &                                          TEMPK,
     &                                          USTAR(I,J),
     &                                          DIAM,  
     &                                          DEN )

               ELSE 

                  !=====================================================
                  ! Replace original code to statement 160 here: only
                  ! do this for non-size-resolved species where 
                  ! AIROSOL(K)=T. (rjp, tdf, bec, bmy, 4/20/04)
                  !=====================================================
!                  ! use Zhang et al for all aerosols (hotp 10/26/07)
!                  VDS = 0.002D0*USTAR(I,J)
!                  IF (OBK(I,J) .LT. 0.0D0) THEN
!                     VDS = VDS*(1.D0+(-300.D0/OBK(I,J))**0.6667D0)
!                  ENDIF
!C***                               
!                  IF ( OBK(I,J) .EQ. 0.0D0 )
!     c                 WRITE(6,156) OBK(I,J),I,J,LDT
! 156              FORMAT(1X,'OBK(I,J)=',E11.2,1X,' I,J =',2I4,
!     c                   1X,'LDT=',I3/) 
!                  CZH  = ZH(I,J)/OBK(I,J)
!                  IF (CZH.LT.-30.0D0) VDS = 0.0009D0*USTAR(I,J)*
!     x                                 (-CZH)**0.6667D0

                  RSURFC(K, LDT) =
     &                ADUST_SFCRSII(K, II, PRESSU(I,J)*1e-3_f8, 
     &                  TEMPK, USTAR(I,J))
                  
!*                                 
!*    Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
!*    VDSMAX is taken from Table 2 of Walcek et al. [1986].
!*    Invert to get corresponding R

!                  RSURFC(K,LDT) = 1.D0/MIN(VDS, DBLE(IVSMAX(II))/1.D4)
               ENDIF

               ! Free pointer
               SpcInfo => NULL()

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
               RSURFC(K,LDT)= MAX(1.e+0_f8, MIN(RSURFC(K,LDT), 
     &                        9999.e+0_f8))
 180        CONTINUE
 190     CONTINUE
!*                                 
!*    Loop through the different landuse types present in the grid square     
!*           
         DO 500 LDT = 1, IREG(I,J)
            IF ( IUSE(I,J,LDT) == 0 ) GOTO 500
            IOLSON = ILAND(I,J,LDT) + 1

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

            CKUSTR = VON_KARMAN * USTAR(I,J)

            REYNO = USTAR(I,J)*ZO(I,J)/XNU  

            IF ( OBK(I,J) .EQ. 0.0e+0_f8 )
     c           WRITE(6,211) OBK(I,J),I,J,LDT                
 211        FORMAT(1X,'OBK(I,J)=',E11.2,1X,' I,J = ',2I4,1X,
     c           'LDT=',I3/) 
            CORR1 = CZ/OBK(I,J)

            ! Define Z0OBK
            Z0OBK = ZO(I,J)/OBK(I,J) 
!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!            ! Set CASTNET height at 10m (krt, 11/2017)
#if defined( MODEL_GEOS )
            Z0OBK2M  = MAX(Z0OBK,  2e+0_fp/OBK(I,J) )
            Z0OBK10M = MAX(Z0OBK, 10e+0_fp/OBK(I,J) )
#else
!            Z0OBKKT = 10e+0_fp/OBK(I,J)
#endif
!==============================================================================

            LRGERA(I,J) = .FALSE.
            ! Add option for non-local PBL (Lin, 03/31/09) 
            IF (.NOT. Input_Opt%LNLPBL) THEN
               IF (CORR1 .GT. 0.e+0_f8) THEN
                  IF (CORR1 .GT.  1.5e+0_f8) LRGERA(I,J) = .TRUE.
               ELSEIF(CORR1 .LE. 0.e+0_f8) THEN
                  IF (CORR1 .LE. -2.5e+0_f8) CORR1 = -2.5e+0_f8
                  CORR2 = LOG(-CORR1)
               ENDIF
            ENDIF
!*                                 
            IF (CKUSTR.EQ.0.0e+0_f8) THEN

               WRITE(6,212) I,J,CKUSTR,VON_KARMAN,USTAR(I,J)

 212           FORMAT(1X,'I,J= ',2I4,1X,'CKUSTR=',E10.1,1X,
     x              'VON_KARMAN= ',E12.4,1X,'USTAR(I,J)= ',

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
            IF ( REYNO < 0.1e+0_f8 ) GOTO 220

            ! Add option for non-local PBL (Lin, 03/31/09) 
            IF (.NOT. Input_Opt%LNLPBL) THEN

!...aerodynamically rough surface.
!*                                 
               IF (CORR1.LE.0.0e+0_f8 .AND. Z0OBK .LT. -1.e+0_f8)THEN
!*... unstable condition; set RA to zero. (first implemented in V. 3.2)
                  RA    = 0.e+0_f8
#if defined( MODEL_GEOS )
                  !RAKT  = 0.e+0_f8
                  RA2M  = 0.e+0_f8
                  RA10M = 0.e+0_f8
#endif

!*... error trap: prevent CORR1 or Z0OBK from being zero or close to
!*... zero (ckeller, 3/15/16)
               ELSEIF ( ABS(CORR1)<=SMALL .OR. ABS(Z0OBK)<=SMALL ) THEN
                  RA = 0.e+0_f8
#if defined( MODEL_GEOS )
                  !RAKT  = 0.e+0_f8
                  RA2M  = 0.e+0_f8
                  RA10M = 0.e+0_f8
#endif

               ELSEIF (CORR1.LE.0.0e+0_f8 .AND. Z0OBK .GE. -1.e+0_f8) 
     &          THEN
!*... unstable conditions; compute Ra as described above.
                  DUMMY1 = (1.e+0_f8 - 9e+0_f8*CORR1)**0.5e+0_f8
                  DUMMY2 = (1.e+0_f8 - 9e+0_f8*Z0OBK)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  DUMMY4 = ABS((DUMMY2 - 1.e+0_f8)/(DUMMY2 + 1.e+0_f8))
                  RA = 0.74e+0_f8* (1.e+0_f8/CKUSTR) * 
     &                 LOG(DUMMY3/DUMMY4)
#if defined( MODEL_GEOS )
!*... correct ozone to 10m for diagnostic (krt, 11/2017)
                  !DUMMY2KT=(1.e+0_f8 - 9e+0_f8*Z0OBKKT)**0.5e+0_f8
                  !DUMMY4KT=ABS((DUMMY2KT-1.e+0_f8)/(DUMMY2KT+1.e+0_f8))
                  !RAKT = 0.74e+0_f8* (1.e+0_f8/CKUSTR) *
     &            !     LOG(DUMMY3/DUMMY4KT)
                  ! 2M
                  DUMMY1 = (1.e+0_f8 - 9e+0_f8*Z0OBK2M)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  RA2M   = 0.74e+0_f8* (1.e+0_f8/CKUSTR) *
     &                     LOG(DUMMY3/DUMMY4)
                  DUMMY1 = (1.e+0_f8 - 9e+0_f8*Z0OBK10M)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  RA10M  = 0.74e+0_f8* (1.e+0_f8/CKUSTR) *
     &                     LOG(DUMMY3/DUMMY4)
#endif                       
          
               ELSEIF((CORR1.GT.0.0e+0_f8).AND.(.NOT.LRGERA(I,J)))
     &          THEN
!*...moderately stable conditions (z/zMO <1); compute Ra as described above
                  RA = (1e+0_f8/CKUSTR) * 
     &            (.74e+0_f8*LOG(CORR1/Z0OBK) + 4.7e+0_f8*(CORR1-Z0OBK))
#if defined( MODEL_GEOS )
     &            !(.74_f8*LOG(CORR1/Z0OBKKT)+4.7_f8*(CORR1-Z0OBKKT))
                  RA2M = (1e+0_f8/CKUSTR) *
     &            (0.74_f8*LOG(Z0OBK2M/Z0OBK)+4.7_f8*(Z0OBK2M-Z0OBK))
                  RA10M = (1e+0_f8/CKUSTR) *
     &            (0.74_f8*LOG(Z0OBK10M/Z0OBK)+4.7_f8*(Z0OBK10M-Z0OBK))
#endif

               ELSEIF(LRGERA(I,J)) THEN
!*... very stable conditions
                  RA    = 1.e+04_f8
#if defined( MODEL_GEOS )
                  !RAKT = 1.e+04_f8
                  RA2M  = 1.e+04_f8
                  RA10M = 1.e+04_f8
#endif
               ENDIF
!* check that RA is positive; if RA is negative (as occasionally
!* happened in version 3.1) send a warning message.

            ELSE

               IF (CORR1.LT.0.0e+0_f8) THEN
!*... unstable conditions; compute Ra as described above.
                  !coef_a=1.e+0_f8
                  !coef_b=15.e+0_f8
                  DUMMY1 = (1.D0 - 15.e+0_f8*CORR1)**0.5e+0_f8
                  DUMMY2 = (1.D0 - 15.e+0_f8*Z0OBK)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  DUMMY4 = ABS((DUMMY2 - 1.e+0_f8)/(DUMMY2 + 1.e+0_f8))
                  RA = 1.e+0_f8 * (1.e+0_f8/CKUSTR) * LOG(DUMMY3/DUMMY4)
!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!!*... correct ozone to 10m for diagnostic (krt, 11/2017)
!                  DUMMY2KT = (1.D0 - 15.e+0_f8*Z0OBKKT)**0.5e+0_f8
!                  DUMMY4KT = ABS((DUMMY2KT - 1.e+0_f8)/
!     &                 (DUMMY2KT + 1.e+0_f8))
!                  RAKT = 1.e+0_f8 * (1.e+0_f8/CKUSTR) *
!     &                 LOG(DUMMY3/DUMMY4KT)
!==============================================================================
#if defined( MODEL_GEOS )
                  ! 2M
                  DUMMY1 = (1.D0 - 15.e+0_f8*Z0OBK2M)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  RA2M = 1.e+0_f8 * (1.e+0_f8/CKUSTR)*LOG(DUMMY3/DUMMY4)
                  ! 10M
                  DUMMY1 = (1.D0 - 15.e+0_f8*Z0OBK10M)**0.5e+0_f8
                  DUMMY3 = ABS((DUMMY1 - 1.e+0_f8)/(DUMMY1 + 1.e+0_f8))
                  RA10M= 1.e+0_f8 * (1.e+0_f8/CKUSTR)*LOG(DUMMY3/DUMMY4)
#endif

               ELSEIF((CORR1.GE.0.0e+0_f8).AND.(CORR1.LE.1.0e+0_f8)) 
     &          THEN
                  !coef_a=1.e+0_f8
                  !coef_b=5.e+0_f8
                  RA = (1D0/CKUSTR) *
     &          (1.e+0_f8*LOG(CORR1/Z0OBK) + 5.e+0_f8*(CORR1-Z0OBK))
!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!                  RAKT = (1D0/CKUSTR) *
!     &          (1.e+0_f8*LOG(CORR1/Z0OBKKT) + 5.e+0_f8*(CORR1-Z0OBKKT))
!==============================================================================
#if defined( MODEL_GEOS )
                  RA2M = (1D0/CKUSTR) *
     &          (1.e+0_f8*LOG(Z0OBK2M/Z0OBK)+5.e+0_f8*(Z0OBK2M-Z0OBK))
                  RA10M = (1D0/CKUSTR) *
     &          (1.e+0_f8*LOG(Z0OBK10M/Z0OBK)+5.e+0_f8*(Z0OBK10M-Z0OBK))
#endif

               ELSE ! CORR1 .GT. 1.0D0
                  !coef_a=5e+0_f8
                  !coef_b=1.e+0_f8
                  RA = (1D0/CKUSTR) *
     &          (5.e+0_f8*LOG(CORR1/Z0OBK) + 1.e+0_f8*(CORR1-Z0OBK))
!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!                  RAKT = (1D0/CKUSTR) *
!     &          (5.e+0_f8*LOG(CORR1/Z0OBKKT) +1.e+0_f8*(CORR1-Z0OBKKT))
!==============================================================================
#if defined( MODEL_GEOS )
                  RA2M = (1D0/CKUSTR) *
     &          (5.e+0_f8*LOG(Z0OBK2M/Z0OBK)+1.e+0_f8*(Z0OBK2M-Z0OBK))
                  RA10M = (1D0/CKUSTR) *
     &          (5.e+0_f8*LOG(Z0OBK10M/z0OBK)+1.e+0_f8*(Z0OBK10M-Z0OBK))
#endif

               ENDIF

            ENDIF

            RA   = MIN(RA,1.e+4_f8)
!==============================================================================
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!            RAKT = MIN(RAKT,1.e+4_f8)
!==============================================================================
#if defined( MODEL_GEOS )
            RA2M   = MIN(RA2M,  1.e+4_f8)
            RA10M  = MIN(RA10M, 1.e+4_f8)
#endif

            ! If RA is < 0, set RA = 0 (bmy, 11/12/99)
            IF (RA .LT. 0.e+0_f8) THEN
               WRITE (6,1001) I,J,RA,CZ,ZO(I,J),OBK(I,J)  
               RA = 0.0e+0_f8
            ENDIF

#if defined( MODEL_GEOS )
            ! Adjust 2M Ra if needed
            IF ( RA2M  < 0.e+0_f8 ) RA2M  = 0.e+0_f8
            IF ( RA10M < 0.e+0_f8 ) RA10M = 0.e+0_f8
#endif

!==============================================================================
!#if defined ( BPCH_DIAG )
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!            IF (RAKT .LT. 0.e+0_f8) THEN
!               WRITE (6,1001) I,J,RAKT,CZ,ZO(I,J),OBK(I,J)  
!               RAKT = 0.0e+0_f8
!            ENDIF
!            ! Store Ra for diagnostic
!            AD_RA(I,J) = RAKT
!            ! Store z/L to determine if M/O theory applies
!            AD_ZL(I,J) = CZ1(I,J)/OBK(I,J)
!#endif
!==============================================================================
            
 1001       FORMAT('WARNING: RA < 0 IN SUBROUTINE DEPVEL',
     &             2I4,4(1X,E12.5))
!* Get total resistance for deposition - loop over species.
            DO 215 K = 1,NUMDEP 
               IF (.NOT.LDEP(K)) GOTO 215
!** DAIR is the thermal diffusivity of air; value of 0.2*1.E-4 m2 s-1 
!** cited on p. 16,476 of Jacob et al. [1992]
               DAIR = 0.2e0_f8*1.e-4_f8
               RB = (2.e+0_f8/CKUSTR)*
     &              (DAIR/DIFFG(TEMPK,PRESSU(I,J),XMW(K)))
     &              **0.667e+0_f8
               IF (AIROSOL(K)) RB=0.e+0_f8
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
#if defined( MODEL_GEOS )
                  RA2M   = 1.0D4
                  RA10M  = 1.0D4
#endif
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

#if defined( MODEL_GEOS )
!---------- Eventually archive aerodynamic resistance. 
!---------- Convert s m-1 to s cm-1. 
!---------- Ra is set to an arbitrary large value of 1.0e4
!---------- in stable conditions. Adjust archived Ra's 
!---------- downward to avoid excessive surface concentration
!---------- adjustments when using C'=(1-Ra*vd)*C.
!---------- (ckeller, 2/2/18)
            !IF ( RA2M >= 1.0d4 ) RA2M = 0.0
            State_Chm%DryDepRa2m(I,J) = State_Chm%DryDepRa2m(I,J) +
     &                                  0.001d0 * DBLE(IUSE(I,J,LDT)) *
     &                                   ( MAX(0.0d0,RA-RA2M)/100.0d0 )
            !IF ( RAKT >= 1.0d4 ) RAKT = 0.0
            State_Chm%DryDepRa10m(I,J) = State_Chm%DryDepRa10m(I,J) +
     &                                   0.001d0 * DBLE(IUSE(I,J,LDT)) *
     &                                   ( MAX(0.0d0,RA-RA10M)/100.0d0 )
#endif

 500     CONTINUE

!** Load array DVEL
         DO 550 K=1,NUMDEP
            IF (.NOT.LDEP(K)) GOTO 550
            DVEL(I,J,K) = VD(K)
            
            ! Now check for negative deposition velocity 
            ! before returning to calling program (bmy, 4/16/00)
            ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
            IF ( DVEL(I,J,K) < 0e+0_f8 ) THEN
!$OMP CRITICAL
               PRINT*, 'DEPVEL: Deposition velocity is negative!'
               PRINT*, 'Dep. Vel = ', DVEL(I,J,K)
               PRINT*, 'Species  = ', K
               PRINT*, 'I, J     = ', I, J
               PRINT*, 'RADIAT   = ', RADIAT(I,J)
               PRINT*, 'TEMP     = ', TEMP(I,J)
               PRINT*, 'SUNCOS   = ', SUNCOS(I,J)
               PRINT*, 'USTAR    = ', USTAR(I,J)
               PRINT*, 'CZ1      = ', CZ1(I,J)
               PRINT*, 'OBK      = ', OBK(I,J)
               PRINT*, 'CFRAC    = ', CFRAC(I,J)
               PRINT*, 'ZH       = ', ZH(I,J)
               PRINT*, 'LRGERA   = ', LRGERA(I,J)
               PRINT*, 'ZO       = ', ZO(I,J)
               PRINT*, 'STOP in depvel.f!'
               CALL CLEANUP
               STOP
!$OMP END CRITICAL
            ENDIF

            ! Now check for IEEE NaN (not-a-number) condition 
            ! before returning to calling program (bmy, 4/16/00)
            ! Also call CLEANUP to deallocate arrays (bmy, 10/15/02)
            IF ( IT_IS_NAN( DVEL(I,J,K) ) ) THEN
!$OMP CRITICAL
               PRINT*, 'DEPVEL: Deposition velocity is NaN!'
               PRINT*, 'Dep. Vel = ', DVEL(I,J,K)
               PRINT*, 'Species  = ', K
               PRINT*, 'I, J     = ', I, J
               PRINT*, 'RADIAT   = ', RADIAT(I,J)
               PRINT*, 'TEMP     = ', TEMP(I,J)
               PRINT*, 'SUNCOS   = ', SUNCOS(I,J)
               PRINT*, 'USTAR    = ', USTAR(I,J)
               PRINT*, 'CZ1      = ', CZ1(I,J)
               PRINT*, 'OBK      = ', OBK(I,J)
               PRINT*, 'CFRAC    = ', CFRAC(I,J)
               PRINT*, 'ZH       = ', ZH(I,J)
               PRINT*, 'LRGERA   = ', LRGERA(I,J)
               PRINT*, 'ZO       = ', ZO(I,J)
               CALL CLEANUP
               STOP
!$OMP END CRITICAL
            ENDIF
 550     CONTINUE
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#if defined( MODEL_GEOS )
      ! Pass Monin Obukhov diagnostics
      IF ( ASSOCIATED ( State_Diag%MoninObukhov ) ) THEN
         State_Diag%MoninObukhov(:,:) = OBK(:,:)
      ENDIF
#endif

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
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(IN) :: TK     ! Temperature [K]
      REAL(f8), INTENT(IN) :: PRESS  ! Pressure [Pa]
      REAL(f8), INTENT(IN) :: XM     ! Molecular weight of gas [kg]
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
!        Also now force REAL(f8) precision with D exponents.  Now use F90
!        style syntax and updated comments. (bmy, 5/16/06)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  06 Jan 2016 - E. Lundgren - Now use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(f8)             :: AIRDEN, Z, DIAM, FRPATH, SPEED, DIFF_G    
!
! !DEFINED PARAMETERS:
!
      REAL(f8), PARAMETER  :: XMAIR  = 28.8e-3_f8 ! Moist air molec wt? (ewl)
      REAL(f8), PARAMETER  :: RADAIR = 1.2e-10_f8
      !REAL(f8), PARAMETER  :: PI     = 3.1415926535897932e+0_f8 ! Use global PI
      REAL(f8), PARAMETER  :: RADX   = 1.5e-10_f8
      !REAL(f8), PARAMETER  :: RGAS   = 8.32e+0_f8     ! Now use global RSTARG
      !REAL(f8), PARAMETER  :: AVOGAD = 6.023e+23_f8   ! Now use global AVO

      !=================================================================
      ! DIFFG begins here!
      !=================================================================

      ! Air density [molec/m3]
      AIRDEN = ( PRESS * AVO ) / ( RSTARG * TK )

      ! DIAM is the collision diameter for gas X with air.
      DIAM   = RADX + RADAIR

      ! Calculate the mean free path for gas X in air: 
      ! eq. 8.5 of Seinfeld [1986];
      Z      = XM  / XMAIR
      FRPATH = 1e+0_f8 /( PI * SQRT( 1e+0_f8 + Z ) * AIRDEN*
     &         ( DIAM**2 ) )

      ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED  = SQRT( 8e+0_f8 * RSTARG * TK / ( PI * XM ) )

      ! Calculate diffusion coefficient of gas X in air; 
      ! eq. 8.9 of Seinfeld [1986]
      DIFF_G = ( 3e+0_f8 * PI / 32e+0_f8 ) * ( 1e+0_f8 + Z ) * 
     &         FRPATH * SPEED

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
      SUBROUTINE READ_DRYDEP_INPUTS( am_I_Root, Input_Opt,
     &                               DRYCOEFF,  IOLSON,   IDEP,      
     &                               IWATER,    NWATER,   IZO,       
     &                               IDRYDEP,   IRI,      IRLU,      
     &                               IRAC,      IRGSS,    IRGSO,     
     &                               IRCLS,     IRCLO,    IVSMAX,    
     &                               RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput

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
      LOGICAL,        INTENT(IN) :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
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
      REAL(fp), INTENT(OUT) :: DRYCOEFF(NPOLY    ) 
      INTEGER,  INTENT(OUT) :: IOLSON  (NSURFTYPE)       
      INTEGER,  INTENT(OUT) :: IDEP    (NSURFTYPE)      
      INTEGER,  INTENT(OUT) :: IWATER  (NSURFTYPE)    
      INTEGER,  INTENT(OUT) :: NWATER        
      INTEGER,  INTENT(OUT) :: IZO     (NSURFTYPE)       
      INTEGER,  INTENT(OUT) :: IDRYDEP (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRI     (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRLU    (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRAC    (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRGSS   (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRGSO   (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRCLS   (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IRCLO   (NDRYDTYPE)
      INTEGER,  INTENT(OUT) :: IVSMAX  (NDRYDTYPE)

      ! Success or failure flag
      INTEGER,  INTENT(OUT) :: RC
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
!  03 Feb 2014 - M. Sulprizio- Change the internal resistance for coniferous
!                              forests to match the internal resistance for
!                              deciduous forests when using the Olson 2001
!                              land map (skim, 1/27/14)
!  18 Dec 2014 - R. Yantosca - Now read DRYCOEFF at REAL*8 precision and
!                              then cast to flexible precision
!  13 Mar 2015 - R. Yantosca - Replace DATA_DIR_1x1 w/ CHEM_INPUTS_DIR
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
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

      ! Shadow variable for reading in data at REAL*8
      REAL(f8)           :: DRYCOEFF_R8(NPOLY)

      !=================================================================
      ! Open and read data from the netCDF file
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Input file for Olson 2001
      nc_file = 'Olson_2001_Drydep_Inputs.nc'

      ! Open file for reading (construct full data path)
      nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 
     &          'Olson_Land_Map_201203/'
      nc_path = TRIM( nc_dir ) // TRIM( nc_file )
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
      st1d   = (/ 1     /)
      ct1d   = (/ NPOLY /)
      CALL NcRd( DRYCOEFF_R8, fId, TRIM(v_name), st1d, ct1d )

      ! Read the DRYCOEFF:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      ! Cast from REAL*8 to flexible precision
      DRYCOEFF = DRYCOEFF_R8

      !----------------------------------------
      ! VARIABLE: IOLSON
      !----------------------------------------

      ! Variable name
      v_name = "IOLSON"

      ! Read IOLSON from file
      st1d   = (/ 1         /)
      ct1d   = (/ NSURFTYPE /)
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
      ct1d   = (/ NSURFTYPE /)
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
      ! NOTE: IWATER is declared with NNSURFTYPE, but has NWATER values
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
      ct1d   = (/ NSURFTYPE /)
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

      ! For Olson 2001 land map, change IRI for coniferous forests
      ! to match IRI for deciduous forests (skim, mps, 2/3/14)
      IRI(3) = 200

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
!  seasalt aerosol species according to Zhang et al 2001.  We account for
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
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: K     ! Drydep species index (range: 1-NUMDEP)
      INTEGER,        INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL(f8),       INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb=100 Pa=0.1 kPa)
      REAL(f8),       INTENT(IN) :: TEMP  ! Temperature [K]   
      REAL(f8),       INTENT(IN) :: USTAR ! Friction velocity [m/s]
      REAL(f8),       INTENT(IN) :: RHB   ! Relative humidity (fraction)
      REAL(f8),       INTENT(IN) :: W10   ! 10 m windspeed [m/s]
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(f8)                   :: RS    ! Surface resistance for particles [s/m]
! 
! !REMARKS:
!  Do computations internally with REAL*8 (8-byte) floating-point precision,
!  in order to avoid a loss of precision.
!
! !REVISION HISTORY: 
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Now limit relative humidity to [tiny(real(f8)),0.99] range for DLOG
!         argument (phs, 6/11/08)
!  (3 ) Bug fixes to the Gerber (1985) growth function (jaegle 5/11/11)
!  (4)  Update growth function to Lewis and Schwartz (2006) and density
!       calculation based on Tang et al. (1997) (bec, jaegle 5/11/11)
!  (5 ) Updates of sea salt deposition over water to follow the Slinn & Slinn
!       (1980) formulation over water surface. Described in Jaegle et al. (ACP,
!       11, 2011) (jaegle 5/11/11)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  14 Jun 2013 - R. Yantosca - Now pass Input_Opt via the arg list
!  06 Jan 2015 - E. Lundgren - Use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: N
      REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8,  C2 = 3.079e+0_f8, 
     &                         C3 = 2.573e-11_f8, C4 = -1.424e+0_f8
      !REAL(f8), PARAMETER   :: G0 = 9.8e+0_f8  ! Now use global g0
      REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
      !REAL(f8), PARAMETER   :: BOLTZ = 1.381e-23_f8  ! Boltzmann constant (J/K)
                                                      ! Now use global BOLTZ
      REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
      REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL(f8)  :: DP          ! Diameter of aerosol [um]
      REAL(f8)  :: PDP         ! Press * Dp      
      REAL(f8)  :: CONST       ! Constant for settling velocity calculations
      REAL(f8)  :: SLIP        ! Slip correction factor
      REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
      REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)
      REAL(f8)  :: RHBL        ! Relative humidity local

      ! replace RCM with RUM (radius in microns instead of cm) - jaegle 5/11/11
      !REAL(f8)  :: DIAM, DEN, RATIO_R, RWET, RCM
      REAL(f8)  :: DIAM, DEN, RATIO_R, RWET, RUM
      REAL(f8)  :: FAC1, FAC2
      REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS
      ! New variables added (jaegle 5/11/11)
      REAL(f8)  :: SW
      REAL(f8)  :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT, DMIDW ! for weighting the settling velocity
      REAL(f8)  :: D0, D1  !lower and upper bounds of sea-salt dry diameter bins 
      REAL(f8)  :: DEDGE
      REAL(f8)  :: DEN1, WTP 
      INTEGER   :: ID,NR
      LOGICAL, SAVE          :: FIRST = .TRUE.

      !increment of radius for integration of settling velocity (um)
      REAL(f8), PARAMETER      :: DR    = 5.e-2_f8

      ! Parameters for polynomial coefficients to derive seawater
      ! density. From Tang et al. (1997) - jaegle 5/11/11
      REAL(f8),  PARAMETER     :: A1 =  7.93e-3_f8
      REAL(f8),  PARAMETER     :: A2 = -4.28e-5_f8
      REAL(f8),  PARAMETER     :: A3 =  2.52e-6_f8
      REAL(f8),  PARAMETER     :: A4 = -2.35e-8_f8
      REAL(f8),  PARAMETER     :: EPSI = 1.0e-4_f8

      ! parameters for assumed size distribution of accumulation and coarse
      ! mode sea salt aerosols, as described in Jaegle et al. (ACP, 11, 2011)
      ! (jaegle, 5/11/11)
      ! 1) geometric dry mean diameters (microns)
      REAL(f8),  PARAMETER     ::   RG_A = 0.085e+0_f8
      REAL(f8),  PARAMETER     ::   RG_C = 0.4e+0_f8
      ! 2) sigma of the size distribution
      REAL(f8),  PARAMETER     ::   SIG_A = 1.5e+0_f8
      REAL(f8),  PARAMETER     ::   SIG_C = 1.8e+0_f8
      !REAL(f8),  PARAMETER     :: PI =3.14159e+0_f8 ! Now use global PI

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
      REAL(f8)  :: 
     & ALPHA(15) = (/ 
     &  1.0e+0_f8,  0.6e+0_f8,   1.1e+0_f8,   0.8e+0_f8, 0.8e+0_f8,  
     &  1.2e+0_f8,  1.2e+0_f8,  50.0e+0_f8,  50.0e+0_f8, 1.3e+0_f8, 
     &  2.0e+0_f8, 50.0e+0_f8, 100.0e+0_f8, 100.0e+0_f8, 1.5e+0_f8  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL(f8)  ::
     & GAMMA(15) = (/ 
     &  0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 
     &  0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 
     &  0.54e+0_f8, 0.54e+0_f8, 0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8 /)

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

      REAL(f8)  :: A(15,5)

      REAL(f8)  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,  
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.e-3_f8

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
      !                Rs = 1 / (Cd/VON_KARMAN*U10m*(Eb+Eim)+VTS)
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
     &          + 0.5e+0_f8 )

      ! Particle radius [cm]
      ! Bug fix: The Gerber [1985] growth should use the dry radius 
      ! in micromenters and not cm. Replace RCM with RUM (jaegle 5/11/11)
      !RCM  = A_RADI(K) * 1.d2
      RUM  = A_RADI(K) * 1.e+6_f8

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
      !RWET    = 0.01e+0_f8*(FAC1/(FAC2-DLOG(RHBL))+RCM**3.e+0_f8)**0.33e+0_f8
      !RWET = 1.d-6*(FAC1/(FAC2-LOG10(RHBL))+RUM**3.e+0_f8)**0.33333e+0_f8

      ! Use equation 5 in Lewis and Schwartz (2006) for sea salt growth [m] (jaegle 5/11/11)
      RWET = A_RADI(K) * (4.e+0_f8 / 3.7e+0_f8) *
     & ( (2.e+0_f8 - RHBL)/(1.e+0_f8 - RHBL) )**(1.e+0_f8/3.e+0_f8)

      ! Ratio dry over wet radii at the cubic power
      !RATIO_R = ( A_RADI(K) / RWET )**3.e+0_f8

      ! Diameter of the wet aerosol [m]
      DIAM  = RWET * 2.e+0_f8  

      ! Density of the wet aerosol [kg/m3] (bec, 12/8/04)
      !DEN   = RATIO_R * A_DEN(K) + ( 1.e+0_f8 - RATIO_R ) * 1000.e+0_f8 

      ! Above density calculation is chemically unsound because it ignores
      ! chemical solvation.  
      ! Iteratively solve Tang et al., 1997 equation 5 to calculate density of
      ! wet aerosol (kg/m3) 
      ! (bec, 6/17/10, jaegle 5/11/11)
      ! Redefine RATIO_R
      RATIO_R = A_RADI(K) / RWET

      ! Assume an initial density of 1000 kg/m3
      DEN  = 1000.e+0_f8
      DEN1 = 0.e+0_f8 !initialize (bec, 6/21/10)
      DO WHILE ( ABS( DEN1-DEN ) .gt. EPSI )
         ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet) 
         WTP    = 100.e+0_f8 * A_DEN(K)/DEN * RATIO_R**3.e+0_f8
         ! Then calculate density of wet aerosol using equation 5 
         ! in Tang et al., 1997 [kg/m3]
         DEN1   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + 
     $               (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8

         ! Now calculate new weight percent using above density calculation
         WTP    = 100.e+0_f8 * A_DEN(K)/DEN1 * RATIO_R**3
         ! Now recalculate new wet density [kg/m3]
         DEN   = ( 0.9971e+0_f8 + (A1 * WTP) + (A2 * WTP**2) + 
     $              (A3 * WTP**3) + (A4 * WTP**4) ) * 1000.e+0_f8
      ENDDO

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.e+6_f8 
 
      ! Constant for settling velocity calculation
      CONST = DEN * DIAM**2 * g0 / 18.e+0_f8

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
      SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * 
     &       EXP( -0.059e+0_f8 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928e+0_f8  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! This settling velocity is for the mid-point of the size bin.
      ! Need to integrate over the size bin, taking into account the
      ! mass distribution of sea-salt and the dependence of VTS on aerosol
      ! size. See WET_SETTLING in SEASALT_MOD.f for more details.
      ! (jaegle 5/11/11)
      SALT_MASS_TOTAL = 0e+0_f8
      VTS_WEIGHT      = 0e+0_f8
      ! Check what the min/max range of the SS size bins are
      IF ( RUM .le. Input_Opt%SALA_REDGE_um(2) ) THEN
        D0 = Input_Opt%SALA_REDGE_um(1)*2e+0_f8
        D1 = Input_Opt%SALA_REDGE_um(2)*2e+0_f8
      ELSE 
        D0 = Input_Opt%SALC_REDGE_um(1)*2e+0_f8
        D1 = Input_Opt%SALC_REDGE_um(2)*2e+0_f8
      ENDIF
     

      DO ID = 1, NR
      ! Calculate mass of wet aerosol (Dw = wet diameter, D = dry diamter):
      ! Overall = dM/dDw = dV/dlnD * Rwet/Rdry * DEN /Rw
        IF (DMID(ID) .ge. D0 .and. DMID(ID) .le. D1 ) THEN
           DMIDW = DMID(ID) * RWET/A_RADI(K)   ! wet radius [um]
           SALT_MASS   = SALT_V(ID) * RWET/A_RADI(K) * DEN / 
     &                    (DMIDW*0.5e+0_f8)
           VTS_WEIGHT  = VTS_WEIGHT + 
     &              SALT_MASS * VTS * (DMIDW/(RWET*1d6*2e+0_f8) )**
     &              2e+0_f8 * (2e+0_f8 * DR *  RWET/A_RADI(K))
           SALT_MASS_TOTAL=SALT_MASS_TOTAL+SALT_MASS *
     &                            (2e+0_f8 * DR *  RWET/A_RADI(K))
        ENDIF

      ENDDO

      ! Final mass weighted setting velocity:
      VTS = VTS_WEIGHT/SALT_MASS_TOTAL


      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP 
     &      / (3.e+0_f8 * PI * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.e+0_f8/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0e+0_f8 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface 
         EIN  = 0e+0_f8
      ELSE
         ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
         EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
      ENDIF

      ! Use the formulation of Slinn and Slinn (1980) for the impaction over
      ! water surfaces (jaegle 5/11/11)
      IF (LUC == 14) THEN
         EIM  = 10.e+0_f8**( -3.e+0_f8/ ST )   ! for water surfaces
      ELSE
         EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)
         EIM  = MIN( EIM, 0.6e+0_f8 )
      ENDIF

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.e+0_f8
      ELSE
         R1 = EXP( -1e+0_f8 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      ! Use the formulation of Slinn and Slinn (1980) for the impaction over
      ! water surfaces (jaegle 5/11/11)
      IF (LUC == 14) THEN
         RS   = 1.e+0_f8 / (USTAR**2.e+0_f8/ (W10*VON_KARMAN) * 
     &          (EB + EIM ) + VTS)
      ELSE
         RS   = 1.e+0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )
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
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt
! 
! !REVISION HISTORY: 
!  11 May 2011 - L. Jaegle   - Initial version
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  14 Jun 2013 - R. Yantosca - Now accept Input_Opt via tha argument list
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: N
      REAL(f8)            :: SALT_MASS, SALT_MASS_TOTAL, VTS_WEIGHT
      REAL(f8)            :: DEDGE
      INTEGER             :: ID,NR
!
! !DEFINED PARAMETERS:
!
      ! increment of radius for integration of settling velocity (um)
      REAL(f8), PARAMETER :: DR    = 5.e-2_f8

      ! parameters for assumed size distribution of acc and coarse mode
      ! sea salt aerosols
      ! geometric dry mean diameters (microns)
      REAL(f8), PARAMETER :: RG_A  = 0.085e+0_f8
      REAL(f8), PARAMETER :: RG_C  = 0.4e+0_f8
      ! sigma of the size distribution
      REAL(f8), PARAMETER :: SIG_A = 1.5e+0_f8
      REAL(f8), PARAMETER :: SIG_C = 1.8e+0_f8
      !REAL(f8), PARAMETER :: PI    = 3.14159e+0_f8  ! Use global PI

      ! Number of bins between the lowest bound of of the accumulation mode
      ! sea salt and the upper bound of the coarse mode sea salt.
      NR =INT((( Input_Opt%SALC_REDGE_um(2) - 
     &           Input_Opt%SALA_REDGE_um(1) ) / DR )
     &                  + 0.5e+0_f8 )

      !=================================================================
      ! Define the volume size distribution of sea-salt. This only has
      ! to be done once. We assume that sea-salt is the combination of a
      ! coarse mode and accumulation model log-normal distribution functions
      !=================================================================

       ! Lower edge of 0th bin diameter [um]
       DEDGE=Input_Opt%SALA_REDGE_um(1) * 2e+0_f8

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
           SALT_V(ID) = PI / 6e+0_f8* (DMID(ID)**3) * (
     &         13e+0_f8*exp(-0.5*( LOG(DMID(ID))-
     &         LOG(RG_A*2e+0_f8) )**2e+0_f8/
     &                   LOG(SIG_A)**2e+0_f8 )
     &         /( sqrt(2e+0_f8 * PI) * LOG(SIG_A) )  +
     &         0.8e+0_f8*exp(-0.5*( LOG(DMID(ID))-
     &         LOG(RG_C*2e+0_f8) )**2e+0_f8/
     &                   LOG(SIG_C)**2e+0_f8)
     &         /( sqrt(2e+0_f8 * PI) * LOG(SIG_C) )  )
           ! update the next edge
           DEDGE = DEDGE + DR*2e+0_f8
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
!  dust aerosol species according to Seinfeld et al 96.  We do not consider
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
      INTEGER,  INTENT(IN) :: K       ! Drydep species (range: 1-NUMDEP)
      INTEGER,  INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
      REAL(f8), INTENT(IN) :: PRESS   ! Pressure [kPa]
      REAL(f8), INTENT(IN) :: TEMP    ! Temperature [K]    
      REAL(f8), INTENT(IN) :: USTAR   ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
      REAL(f8)             :: RS      ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY: 
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  06 Jan 2016 - E. Lundgren - Use global physical paramters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: N
      REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8,  C2 = 3.079e+0_f8, 
     &                         C3 = 2.573e-11_f8, C4 = -1.424e+0_f8
      !REAL(f8), PARAMETER   :: G0 = 9.8e+0_f8 ! Use global g0
      REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
      !REAL(f8), PARAMETER   :: BOLTZ = 1.381e-23_f8  ! Boltzmann constant (J/K)
                                                      ! Use global BOLTZ
      REAL(f8), PARAMETER   :: E0 = 1.e+0_f8
      REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL(f8)  :: DP          ! Diameter of aerosol [um]
      REAL(f8)  :: PDP         ! Press * Dp      
      REAL(f8)  :: CONST       ! Constant for settling velocity calculations
      REAL(f8)  :: SLIP        ! Slip correction factor
      REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
      REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL(f8)  :: DIAM, DEN
      REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

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
      DIAM  = A_RADI(K) * 2.e+0_f8 

      ! Particle density [kg/m3]
      DEN   = A_DEN(K)          

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.e+6_f8 
 
      ! Constant for settling velocity calculation
      CONST = DEN * DIAM**2 * g0 / 18.e+0_f8
       
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
      SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * 
     &       EXP( -0.059e+0_f8 * PDP ) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928e+0_f8  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP
     &     / (3.e+0_f8 * Pi * VISC * DIAM)  

      ! Schmidt number and Diffusion term
      SC   = AIRVS / DIFF                            
      EB   = SC**(-0.666667e+0_f8)

      ! Stokes number and impaction term
      ST   = VTS * USTAR * USTAR / ( AIRVS * g0 )
      EIM  = 10.e+0_f8**(-3.e+0_f8 / ST) 

      ! surface resistance for particle
      RS   = 1.e+0_f8 / ( E0 * USTAR * (EB + EIM) )
      
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
      INTEGER,  INTENT(IN) :: K     ! Drydep species index (range: 1-NUMDEP)
      INTEGER,  INTENT(IN) :: II    ! Surface type index of GEOS-CHEM
      REAL(f8), INTENT(IN) :: PRESS ! Pressure [kPa] (1 mb = 100 Pa = 0.1 kPa)
      REAL(f8), INTENT(IN) :: TEMP  ! Temperature [K]    
      REAL(f8), INTENT(IN) :: USTAR ! Friction velocity [m/s]
!
! !RETURN VALUE:
!
      REAL(f8)             :: RS    ! Surface resistance for particles [s/m]
! 
! !REVISION HISTORY:
!  (1 ) Updated comments.  Also now force double precision w/ "D" exponents.
!        (bmy, 4/1/04)
!  (2 ) Renamed to DUST_SFCRSII, since this will only be used to compute
!        aerodynamic resistance of dust aerosols.  (bec, bmy, 4/15/05)
!  (3 ) Modified hotp for non size resolved aerosols. This is just DUST_SFCRSII
!        renamed and the diameter and density fixed. (hotp 7/12/07)
!  22 Dec 2011 - M. Payer    - Added ProTeX headers
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: N
      REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8,  C2 = 3.079e+0_f8, 
     &                         C3 = 2.573e-11_f8, C4 = -1.424e+0_f8
      !REAL(f8), PARAMETER   :: G0 = 9.8e+0_f8 ! Use global g0
      REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
      !REAL(f8), PARAMETER   :: BOLTZ = 1.381e-23_f8  ! Boltzmann constant (J/K)
                                                      ! Use global BOLTZ
      REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
      REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL(f8)  :: DP          ! Diameter of aerosol [um]
      REAL(f8)  :: PDP         ! Press * Dp      
      REAL(f8)  :: CONST       ! Constant for settling velocity calculations
      REAL(f8)  :: SLIP        ! Slip correction factor
      REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
      REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL(f8)  :: DIAM, DEN
      REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

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
      REAL(f8)  :: 
     & ALPHA(15) = (/ 
     & 1.0e+0_f8,  0.6e+0_f8,   1.1e+0_f8,   0.8e+0_f8, 0.8e+0_f8,  
     & 1.2e+0_f8,  1.2e+0_f8,  50.0e+0_f8,  50.0e+0_f8, 1.3e+0_f8, 
     & 2.0e+0_f8, 50.0e+0_f8, 100.0e+0_f8, 100.0e+0_f8, 1.5e+0_f8  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL(f8)  ::
     & GAMMA(15) = (/ 
     & 0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 
     & 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 
     & 0.54e+0_f8, 0.54e+0_f8, 0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8 /)

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

      REAL(f8)  :: A(15,5)

      REAL(f8)  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,  
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     & 10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.e-3_f8

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
      DIAM  = 0.5e-6_f8  

      ! Particle density [kg/m3] hotp 10/26/07
      DEN   = 1500          

      ! Dp [um] = particle diameter
      DP    = DIAM * 1.e+6_f8 
 
      ! Constant for settling velocity calculation
      CONST = DEN * DIAM**2 * g0 / 18.e+0_f8
       
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
      SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * 
     &       EXP( -0.059e+0_f8 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928e+0_f8  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP
     &      / (3.e+0_f8 * Pi * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.e+0_f8/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0e+0_f8 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface 
         EIN  = 0e+0_f8
      ELSE
         ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
         EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
      ENDIF

      EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

      EIM  = MIN( EIM, 0.6e+0_f8 )

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.e+0_f8
      ELSE
         R1 = EXP( -1e+0_f8 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      RS   = 1.e0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

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
!  dust aerosol species according to Zhang et al 2001.  We do not consider the
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
      INTEGER,  INTENT(IN) :: K       ! Drydep species index (range: 1-NUMDEP)
      INTEGER,  INTENT(IN) :: II      ! Surface type index of GEOS-CHEM
      REAL(f8), INTENT(IN) :: PRESS   ! Pressure [kPa]
      REAL(f8), INTENT(IN) :: TEMP    ! Temperature [K]    
      REAL(f8), INTENT(IN) :: USTAR   ! Friction velocity [m/s]
      REAL(f8), INTENT(IN) :: DIAM    ! Particle diameter [m]
      REAL(f8), INTENT(IN) :: DEN     ! Particle density [kg/m3]
!
! !RETURN VALUE:
!
      REAL(f8)             :: RS      ! Surface resistance for particles [s/m]
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
!  06 Jan 2016 - E. Lundgren - Use global physical parameters
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: N
      REAL(f8), PARAMETER   :: C1 = 0.7674e+0_f8,  C2 = 3.079e+0_f8, 
     &                         C3 = 2.573e-11_f8, C4 = -1.424e+0_f8
      !REAL(f8), PARAMETER   :: G0 = 9.8e+0_f8 ! Use global g0
      REAL(f8), PARAMETER   :: BETA  = 2.e+0_f8
      !REAL(f8), PARAMETER   :: BOLTZ = 1.381e-23_f8  ! Boltzmann constant (J/K)
                                                      ! Use global BOLTZ
      REAL(f8), PARAMETER   :: E0 = 3.e+0_f8
      REAL(f8)  :: AIRVS       ! kinematic viscosity of Air (m^2/s)
      REAL(f8)  :: DP          ! Diameter of aerosol [um]
      REAL(f8)  :: PDP         ! Press * Dp      
      REAL(f8)  :: CONST       ! Constant for settling velocity calculations
      REAL(f8)  :: SLIP        ! Slip correction factor
      REAL(f8)  :: VISC        ! Viscosity of air (Pa s)
      REAL(f8)  :: DIFF        ! Brownian Diffusion constant for particles (m2/s)
      REAL(f8)  :: SC, ST      ! Schmidt and Stokes number (nondim)

      REAL(f8)  :: EB, EIM, EIN, R1, AA, VTS

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
      REAL(f8)  :: 
     & ALPHA(15) = (/ 
     & 1.0e+0_f8,  0.6e+0_f8,   1.1e+0_f8,   0.8e+0_f8, 0.8e+0_f8,  
     & 1.2e+0_f8,  1.2e+0_f8,  50.0e+0_f8,  50.0e+0_f8, 1.3e+0_f8, 
     & 2.0e+0_f8, 50.0e+0_f8, 100.0e+0_f8, 100.0e+0_f8, 1.5e+0_f8  /)

      ! Now force to double precision (bmy, 4/1/04)
      REAL(f8)  ::
     & GAMMA(15) = (/ 
     & 0.56e+0_f8, 0.58e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 0.56e+0_f8, 
     & 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 0.54e+0_f8, 
     & 0.54e+0_f8, 0.54e+0_f8, 0.50e+0_f8, 0.50e+0_f8, 0.56e+0_f8 /)

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

      REAL(f8)  :: A(15,5)

      REAL(f8)  :: Aavg(15)

      ! Now force to double precision (bmy, 4/1/04)
      DATA   A /  
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     &  10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     &  10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     &  10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   5.0e+0_f8,  10.0e+0_f8,  5.0e+0_f8,  
     &  5.0e+0_f8,   5.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     &  10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8,
     &
     &  2.0e+0_f8,   5.0e+0_f8,   2.0e+0_f8,   5.0e+0_f8,  5.0e+0_f8,  
     &  2.0e+0_f8,   2.0e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8, 
     &  10.0e+0_f8, -999.e+0_f8, -999.e+0_f8, -999.e+0_f8, 10.0e+0_f8  /

      ! Annual average of A
      Aavg(:) = (A(:,1)+A(:,2)+A(:,3)+A(:,4)+A(:,5))/5.
      LUC     = LUCINDEX(II)
      AA      = Aavg(LUC) * 1.e-3_f8

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
 
      ! Dp [um] = particle diameter
      DP    = DIAM * 1.e+6_f8 
 
      ! Constant for settling velocity calculation
      CONST = DEN * DIAM**2 * g0 / 18.e+0_f8
       
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
      SLIP = 1e+0_f8 + ( 15.60e+0_f8 + 7.0e+0_f8 * 
     &       EXP( -0.059e+0_f8 * PDP) ) / PDP

      !=================================================================
      ! Note, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
      ! which produce slip correction factore with small error
      ! compared to the above with less computation.
      !=================================================================

      ! Viscosity [Pa s] of air as a function of temp (K)
      VISC = 1.458e-6_f8 * (TEMP)**(1.5e+0_f8) / (TEMP + 110.4e+0_f8)

      ! Kinematic viscosity (Dynamic viscosity/Density)
      AIRVS= VISC / 1.2928e+0_f8  

      ! Settling velocity [m/s]
      VTS  = CONST * SLIP / VISC

      ! Brownian diffusion constant for particle (m2/s)
      DIFF = BOLTZ * TEMP * SLIP
     &      / (3.e+0_f8 * Pi * VISC * DIAM)  

      ! Schmidt number 
      SC   = AIRVS / DIFF                            
      EB   = 1.e+0_f8/SC**(gamma(LUC))

       ! Stokes number  
      IF ( AA < 0e+0_f8 ) then
         ST   = VTS * USTAR * USTAR / ( AIRVS * g0 ) ! for smooth surface 
         EIN  = 0e+0_f8
      ELSE
         ST   = VTS   * USTAR / ( g0 * AA )          ! for vegetated surfaces
         EIN  = 0.5e+0_f8 * ( DIAM / AA )**2
      ENDIF

      EIM  = ( ST / ( ALPHA(LUC) + ST ) )**(BETA)

      EIM  = MIN( EIM, 0.6e+0_f8 )

      IF (LUC == 11 .OR. LUC == 13 .OR. LUC == 14) THEN
         R1 = 1.D0
      ELSE
         R1 = EXP( -1e+0_f8 * SQRT( ST ) )
      ENDIF

      ! surface resistance for particle
      RS   = 1.e+0_f8 / (E0 * USTAR * (EB + EIM + EIN) * R1 )

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
      SUBROUTINE INIT_DRYDEP( am_I_Root, Input_Opt, 
     &                        State_Chm, State_Diag, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE Species_Mod,    ONLY : Species
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REMARKS:
!  We now know how many drydep species there are before INIT_DRYDEP is 
!  called.  This allows us to get rid of MAXDEP.  NUMDEP should be
!  equal to State_Chm%nDryDep, otherwise there is an error.
!
!  Also note: we need to use the actual molecular weights instead of
!  the emitted molecular weights.  These are necessary for the Schmidt #
!  computation.
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
!                              prevents XMW=0e+0_f8 from being passed to function
!                              DIFFG, where it is in the denominator.
!  04 Sep 2013 - R. Yantosca - Improve printout of drydep species
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  15 Jan 2015 - R. Yantosca - Now save NTRAIND to Input_Opt%NTRAIND
!  08 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  22 Sep 2015 - R. Yantosca - Now allocate internal variables that used 
!                              to be of size MAXDEP.
!  30 Sep 2015 - R. Yantosca - DD_A_Density is renamed to Density
!  30 Sep 2015 - R. Yantosca - DD_A_Radius is renamed to Radius
!  14 Dec 2015 - R. Yantosca - Now use actual MW (MW_g) instead of emitted MW
!                              (EmMW_g) to define the XMW array
!  20 Jun 2016 - R. Yantosca - Now define a couple of species ID flags here
!                              as module variables on the first call
!  18 Jul 2016 - M. Sulprizio- Remove special handling for ISOPN and MVKN.
!                              Family tracers have been eliminated.
!  13 Jul 2017 - M. Sulprizio- Add scaling of species drydep relative to HNO3,
!                              PAN, and ISOPN (from K. Travis, J. Fisher)
!  05 Oct 2017 - R. Yantosca - Now accept State_Diag as an argument
!  07 Nov 2017 - R. Yantosca - Now return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: LDRYD
      LOGICAL                :: IS_Hg
      INTEGER                :: N

      ! Strings
      CHARACTER(LEN=255)     :: Msg, ErrMsg, ThisLoc

      ! Objects
      TYPE(Species), POINTER :: SpcInfo 

      !=================================================================
      ! INIT_DRYDEP begins here!
      !=================================================================
      
      ! Initialize
      RC        = GC_SUCCESS
      IS_Hg     = Input_Opt%ITS_A_MERCURY_SIM
      LDRYD     = Input_Opt%LDRYD
      NUMDEP    = 0
      id_ACET   = 0
      id_ALD2   = 0
      id_HNO3   = IND_('HNO3'  )
      id_PAN    = IND_('PAN'   )
      id_ISOPND = IND_('ISOPND')
      ErrMsg    = ''
      ThisLoc   = ' -> at Init_Drydep (in module GeosCore/drydep_mod.F)'

      ! If the dry deposition module has already been initialized,
      ! the arrays do not need to be allocated again, as they are only
      ! dependent on the chemistry configuration (State_Chm%nDryDep)
      !
      ! This is necessary for integrating GEOS-Chem with a variable
      ! domain model like WRF-GC, where multiple instances of GEOS-Chem
      ! run in the same CPU. (hplin, 2/16/2019)
      IF ( ALLOCATED( A_DEN ) ) RETURN

      !===================================================================
      ! Arrays that hold information about dry-depositing species
      ! Only allocate these if dry deposition is activated
      !===================================================================
      IF ( State_Chm%nDryDep > 0 ) THEN

         ! Aerosol density [kg/m3]
         ALLOCATE( A_DEN( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:A_DEN', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         A_DEN(:)   = 0e+0_f8

         ! Aerosol radius [um]
         ALLOCATE( A_RADI( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:A_RADI', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         A_RADI = 0e+0_f8

         ! Is the species an aerosol? (T/F)
         ALLOCATE( AIROSOL( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:AIROSOL', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         AIROSOL = .FALSE.

         ! Drydep species name
         ALLOCATE( DEPNAME( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:DEPNAME', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         DEPNAME = ''

         ! Reactivity factor
         ALLOCATE( F0( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:F0', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         F0 = 0e+0_f8

         ! Henry's law K0
         ALLOCATE( HSTAR( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:HSTAR', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         HSTAR = 0e+0_f8

         ! POPs KOA
         ALLOCATE( KOA( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:KOA', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         KOA  = 0e+0_f8

         ! Drydep species indicies
         ALLOCATE( NDVZIND( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:NDVZIND', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         NDVZIND = 0

         ! Drydep scaling flag
         ALLOCATE( FLAG( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:FLAG', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         FLAG = 0

         ! Species indices
         ALLOCATE( NTRAIND( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:NTRAIND', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         NTRAIND = 0

         ! Molecular weight [kg/mol]
         ALLOCATE( XMW( State_Chm%nDryDep ), STAT=RC )
         CALL GC_CheckVar( 'drydep_mod:XMW', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         XMW = 0e+0_f8

      ENDIF

      !=================================================================
      ! First identify species that dry deposit and then initialize 
      ! the dry deposition quantities accordingly:
      !
      ! Quantity  Description
      ! ----------------------------------------------------------------
      ! NUMDEP    Number of dry depositing species
      ! NTRAIND   GEOS-Chem species ID number (advected index)
      ! NDVZIND   Coresponding index in the DVEL drydep velocity array
      ! HSTAR     Henry's law solubility constant [M atm-1]
      ! F0        Reactivity (0.0 = not reactive, 1.0=very reactive)
      ! XMW       Molecular weight of species [kg mol-1]
      ! AIROSOL    = T if the species is aerosol, =F if gas
      ! A_DEN     Aerosol density [kg m-3]
      ! A_RADIUS  Aerosol radius [um]
      ! KOA       POPs KOA parameter
      ! 
      ! NOTES: 
      ! (1) For XMW, we take the species emitted molecular weight
      !      (EmMW_g * 1e-3).  Several hydrocarbons are transported
      !      as equivalent molecules of carbon.  In this case the
      !      EmMw_g will be 12.0.
      ! (2) We have edited the molecular weights of some species to
      !      match the prior code.  Some of these definitions are
      !      inconsistent with the listed molecular weights. (e.g.
      !      ACET is transported as 3 carbons, but we give it a
      !      molecular weight of 58 instead of 12).  These will need
      !      to be researched further.
      ! (3) The deposition names of SO4s and NITs need to be in
      !      uppercase.  Therefore, we overwrite the values from
      !      the species database with SO4S, NITS.
      !=================================================================
      DO N = 1, State_Chm%nAdvect

         ! Point to the Nth species in the species database
         SpcInfo => State_Chm%SpcData(N)%Info

         ! Only proceed if the species dry deposits
         IF ( SpcInfo%Is_Drydep ) THEN

            ! Initialize dry deposition quantities
            NUMDEP            = NUMDEP + 1
            NTRAIND(NUMDEP)   = SpcInfo%ModelID
            NDVZIND(NUMDEP)   = SpcInfo%DryDepID
            DEPNAME(NUMDEP)   = TRIM( SpcInfo%Name )
            HSTAR(NUMDEP)     = DBLE( SpcInfo%DD_Hstar_old   )
            F0(NUMDEP)        = DBLE( SpcInfo%DD_F0          )
            KOA(NUMDEP)       = DBLE( SpcInfo%DD_KOA         )
            XMW(NUMDEP)       = DBLE( SpcInfo%MW_g * 1e-3_fp )
            AIROSOL(NUMDEP)   = ( .not. SpcInfo%Is_Gas       )

            ! Only copy DENSITY if it's not a missing value
            IF ( SpcInfo%Density > 0.0_fp ) THEN
               A_DEN(NUMDEP)  = DBLE( SpcInfo%Density        )
            ENDIF

            ! Only copy RADIUS if it's not a missing value
            IF ( SpcInfo%Radius > 0.0_fp ) THEN
               A_RADI(NUMDEP) = DBLE( SpcInfo%Radius         )
            ENDIF

            !-----------------------------------------------------
            ! Kludges to match behavior of older code
            !-----------------------------------------------------
            SELECT CASE ( TRIM( SpcInfo%Name ) ) 

               CASE( 'ACET' )
                  ! Flag the species ID of ACET for use above.
                  id_ACET = SpcInfo%ModelId
                        
               CASE( 'ALD2' )
                  ! Flag the species ID of ALD2 for use above.
                  id_ALD2 = SpcInfo%ModelId

               CASE( 'NITs', 'NITS' )
                  ! DEPNAME for NITs has to be in all caps, for
                  ! backwards compatibility with older code.
                  DEPNAME(NUMDEP) = 'NITS' 


               CASE( 'N2O5', 'HC187' )
                  ! These species scale to the Vd of HNO3. We will
                  ! explicitly compute the Vd of these species instead
                  ! of assigning the Vd of HNO3 from the DVEL array.
                  ! The scaling is applied in DO_DRYDEP using FLAG=1.
                  !
                  ! Make sure to set XMW to the MW of HNO3
                  ! for the computation of Vd to work properly.
                  XMW(NUMDEP)  = State_Chm%SpcData(id_HNO3)%Info%MW_g
     &                           * 1e-3_fp
                  FLAG(NUMDEP) = 1

               CASE(  'NPMN', 'IPMN', 'PPN', 'R4N2' )
                  ! These specied scale to the Vd of PAN.  We will 
                  ! explicitly compute the Vd of these species instead 
                  ! of assigning the Vd of PAN from the DVEL array.
                  ! The scaling is applied in DO_DRYDEP using FLAG=2.
                  !
                  ! Make sure to set XMW to the MW of PAN
                  ! for the computation of Vd to work properly.
                  XMW(NUMDEP)  = State_Chm%SpcData(id_PAN)%Info%MW_g
     &                           * 1e-3_fp
                  FLAG(NUMDEP) = 2

               CASE( 'DHDN', 'ETHLN', 'ISN1', 'MONITS', 'MONITU', 
     &               'HONIT' )
                  ! These species scale to the Vd of ISOPN. We will
                  ! explicitly compute the Vd of these species instead
                  ! of assigning the Vd of ISOPN from the DVEL array.
                  ! The scaling is applied in DO_DRYDEP using FLAG=3.
                  !
                  ! Make sure to set XMW to the MW of ISOPN
                  ! for the computation of Vd to work properly.
                  XMW(NUMDEP)  = State_Chm%SpcData(id_ISOPND)%Info%MW_g
     &                           * 1e-3_fp
                  FLAG(NUMDEP) = 3

               CASE( 'SO4s', 'SO4S' )
                  ! DEPNAME for SO4s has to be in all caps, for 
                  ! backwards compatibility with older code
                  DEPNAME(NUMDEP) = 'SO4S'

               CASE DEFAULT
                  ! Do nothing

            END SELECT

         ENDIF
            
         ! Free pointer
         SpcInfo => NULL()
      ENDDO

      ! For TOMAS 
      id_NK1   = Ind_('NK1'  )

      !=================================================================
      ! Allocate arrays
      ! add allocation for SALT_V and DMID (jaegle 5/11/11)
      !=================================================================
      ALLOCATE( SALT_V( NR_MAX ), STAT=RC )
      CALL GC_CheckVar( 'drydep_mod:SALT_V', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SALT_V = 0e+0_f8

      ALLOCATE( DMID( NR_MAX ), STAT=RC )
      CALL GC_CheckVar( 'drydep_mod:DMID', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DMID = 0e+0_f8

      ! SDE 2016-03-27: This is now called after Input_Opt has already
      ! been broadcast, so each CPU will do this independently anyway
      ! ckeller: moved to end of routine (see below)
!      Input_Opt%NUMDEP            = NUMDEP
!      Input_Opt%DEPNAME(1:NUMDEP) = DEPNAME(1:NUMDEP)
!      Input_Opt%NTRAIND(1:NUMDEP) = NTRAIND(1:NUMDEP)

      !=================================================================
      ! Echo information to stdout
      !=================================================================
      IF ( am_I_Root ) THEN

         ! Line 1
         MSG = 'INIT_DRYDEP: List of dry deposition species:'
         WRITE( 6, '(/,a)' ) TRIM( MSG )

         ! Line 2
         MSG =  '  #   Name      Species DEPVEL Henry''s    React.'
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
      CALL READ_DRYDEP_INPUTS( am_I_Root, Input_Opt, 
     &                         DRYCOEFF,  INDOLSON, IDEP,
     &                         IWATER,    NWATER,   IZO,      
     &                         IDRYDEP,   IRI,      IRLU,     
     &                         IRAC,      IRGSS,    IRGSO, 
     &                         IRCLS,     IRCLO,    IVSMAX, 
     &                         RC )

       ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Read_Drydep_Inputs"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Calls INIT_WEIGHTSS to calculate the volume distribution of 
      ! sea salt aerosols (jaegle 5/11/11)
      CALL INIT_WEIGHTSS( Input_Opt )
    
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
!  22 Sep 2015 - R. Yantosca - Now deallocate arrays that were formerly
!                              sized to MAXDEP.
!  05 Jul 2018 - H. Lin      - Bug fix: deallocate FLAG
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DRYDEP begins here!
      !=================================================================
      IF ( ALLOCATED( A_DEN    ) ) DEALLOCATE( A_DEN    )
      IF ( ALLOCATED( A_RADI   ) ) DEALLOCATE( A_RADI   )
      IF ( ALLOCATED( AIROSOL  ) ) DEALLOCATE( AIROSOL  )
      IF ( ALLOCATED( DEPNAME  ) ) DEALLOCATE( DEPNAME  )
      IF ( ALLOCATED( DMID     ) ) DEALLOCATE( DMID     )
      IF ( ALLOCATED( F0       ) ) DEALLOCATE( F0       )
      IF ( ALLOCATED( FLAG     ) ) DEALLOCATE( FLAG     )
      IF ( ALLOCATED( HSTAR    ) ) DEALLOCATE( HSTAR    )
      IF ( ALLOCATED( KOA      ) ) DEALLOCATE( KOA      )
      IF ( ALLOCATED( NDVZIND  ) ) DEALLOCATE( NDVZIND  )
      IF ( ALLOCATED( NTRAIND  ) ) DEALLOCATE( NTRAIND  )
      IF ( ALLOCATED( SALT_V   ) ) DEALLOCATE( SALT_V   )
      IF ( ALLOCATED( XMW      ) ) DEALLOCATE( XMW      )

      END SUBROUTINE CLEANUP_DRYDEP
!EOC
      END MODULE DRYDEP_MOD

      
