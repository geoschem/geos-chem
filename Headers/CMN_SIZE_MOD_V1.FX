!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_SIZE_mod.F
!
! !DESCRIPTION: CMN\_SIZE contains size parameters for GEOS-Chem arrays.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_SIZE_MOD
!
! !USES:
!
      ! C Preprocessor #define statements for conditional compilation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% This is already included in CMN_GRID_MOD (bmy, 8/5/11)
!#     include "define.h"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#     include "define.h"
      USE CMN_GRID_MOD

      IMPLICIT NONE
!      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! TRACER & EMISSION SPECIES PARAMETERS
      !
      ! NNPAR   = max number of tracers
      ! NEMPARA = max number of anthropogenic emission species
      ! NEMPARB = max number of biogenic      emission species
      !=================================================================
      ! increase NNPAR and NEMPARA an extra amount (hotp 7/31/09)
#if   defined( TOMAS )
      INTEGER, PARAMETER :: NNPAR   = 320   ! For TOMAS (win, bmy, 1/25/10)
#else
      ! increase NNPAR to 100 (FP 8/2009)
      INTEGER, PARAMETER :: NNPAR   = 75    ! For non-TOMAS simulations
      !INTEGER, PARAMETER :: NNPAR   = 100
#endif

      ! Nempara increased to 26. (fp, 2/8/10)
      ! new emissions HNO3 and O3 (phs)
      !  Add non-biogenic emission species:
      !  BENZ, TOLU, XYLE, C2H2, C2H4, GLYX, MGLY, GLYC, HAC. (tmf, 1/7/09)
      !INTEGER, PARAMETER :: NEMPARA = 21 
      INTEGER, PARAMETER :: NEMPARA = 26 

      ! Add biogenic emissions: MBO, MONX. (tmf, 1/7/09) 
!      INTEGER, PARAMETER :: NEMPARB = 3
      INTEGER, PARAMETER :: NEMPARB = 17

      !=================================================================
      ! OTHER PARAMETERS 
      !=================================================================

      ! NVEGTYPE - Maximum number of surface types: 74 olson
      ! NTYPE    - Maximum number of veg types in a CTM grid box
      ! NPOLY    - Number of coefficients for polynomial fits
      INTEGER, PARAMETER :: NVEGTYPE = 74
      INTEGER, PARAMETER :: NTYPE    = 15
      INTEGER, PARAMETER :: NPOLY    = 20

      ! NNSTA = max number of time series stations (in inptr.ctm)
      INTEGER, PARAMETER :: NNSTA = 800

      ! LLCONVM - Max number of layers for convection
      INTEGER, PARAMETER :: LLCONVM = LLPAR - 1

      ! NOXLEVELS = Number of levels of anthro NOx emission 
      !             (e.g. surface and 100m)
      ! NOXEXTENT = Highest sigma level that receives anthro NOx emission 
      INTEGER, PARAMETER :: NOXLEVELS = 2
      INTEGER, PARAMETER :: NOXEXTENT = 2 

      ! MAXFAM -- Max number of families for prod and loss output
      INTEGER, PARAMETER :: MAXFAM = 40

      ! MAXMEM is maximum number of families of prod and loss
      ! moved from input_mod and diag_pl_mod to here (hotp 7/31/09)
      ! MAXMEM also increased from 10 to 20 by FP
      INTEGER, PARAMETER :: MAXMEM = 20

      ! MAXPL increased from 100 to 500 and moved from diag_pl_mod
      ! to here by FP (hotp 7/31/09)
      INTEGER, PARAMETER :: MAXPL = 500

      ! MAXIJ - Maximum number of 1st level grid boxes
      INTEGER, PARAMETER :: MAXIJ = IIPAR * JJPAR

      ! NDUST -- Number of FAST-J aerosol size bins (rvm, bmy, 11/15/01)
      INTEGER, PARAMETER :: NDUST = 7

      ! NAER -- number of other aerosol categories (rvm, bmy, 2/27/02)
      INTEGER, PARAMETER :: NAER = 5

      ! NRH -- number of relative humidity bins (rvm, bmy, 2/27/02)
      INTEGER, PARAMETER :: NRH = 5

      !NBIOMAX -- biomass burning
      ! increase NBIOMAX to 20 (hotp 7/31/09)
      ! increase NBIOMAX to 24 for dicarbonyls (ccc, 2/02/10)
      INTEGER, PARAMETER :: NBIOMAX = 24

#if   defined( TOMAS )

      ! NDSTBIN -- redimensioned for TOMAS (dwest, bmy, 2/1/10)
      INTEGER, PARAMETER :: NDSTBIN   = 30

      ! Number of TOMAS bins 
      INTEGER, PARAMETER :: TOMASBIN  = 30
      INTEGER, PARAMETER :: TOMASSPEC = 8

#else

      ! NDSTBIN -- number of dust size bins for transport (tdf, bmy, 3/31/04)
      INTEGER, PARAMETER :: NDSTBIN   = 4

#endif
!
! !REVISION HISTORY:
!  (1 ) Now set LLTROP = 20 for GEOS-3 (bmy, 4/12/01)
!  (2 ) Eliminated obsolete commented-out code (bmy, 4/20/01)
!  (3 ) Now set MAXFAM = 12 for more P-L families (bmy, 6/28/01)  
!  (4 ) Comment out {IJL}GCMPAR -- these are obosolete (bmy, 9/24/01)
!  (5 ) Also set LLPAR = 30 for GEOS-3, will regrid online (bmy, 9/24/01) 
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (7 ) Removed NAIR, LAIREMS, these are now defined 
!        in "aircraft_nox_mod.f" (bmy, 2/14/02)
!  (8 ) Eliminated commented-out code from 2/14/02.  Also added NAER 
!        and NRH parameters for aerosols. (rvm, bmy, 2/27/02)
!  (9 ) Removed IM, JM, IMX, JMX to avoid namespace pollution.  This
!        is needed to get the new TPCORE to work.  Also changed RCS
!        ID tag comment character from "C" to "!" to allow freeform
!        compilation. (bmy, 6/25/02)
!  (10) Removed obsolete code from 6/02 (bmy, 8/26/02)
!  (11) Added NUMDEP_SULF in a common block for sulfate dry deposition. 
!        Also set MAXDEP=31 and NNPAR=31 for coupled fullchem/sulfate
!        simulations. (rjp, bdf, bmy, 11/15/02)
!  (12) Removed I0, J0; these are now superseded by "grid_mod.f"
!        (bmy, 2/11/03)
!  (13) Added parameters for GEOS-4 (bmy, 6/18/03)
!  (14) Now defines both 55 level and 30 level GEOS-4 grids.  Also
!        define LLTROP=19 for GEOS-4 grids.  Also remove obsolete
!        GEOS-2 grid declarations. (bmy, 10/31/03)
!  (15) LLTROP should be 17 for GEOS-4...based on the ND55 diagnostic
!        when computed for 2003 met fields (bmy, 2/18/04)
!  (16) Increase NNPAR from 31 to 39 for carbon & dust tracers.  Also
!        declare NDSTBIN as # of dust bins. (rvm, tdf, bmy, 4/1/04)
!  (17) Increase NNPAR to 41 for seasalt tracers (rjp, bec, bmy, 4/20/04)
!  (18) Increase NNPAR to 50 for SOA tracers (rjp, bmy, 7/15/04)
!  (19) Now use NESTED_CH and NESTED_NA cpp switches to define 
!        parameters for 1x1 nested grids.  Also add parameters for
!        the 1 x 1.25 global grid. (bmy, 12/1/04)
!  (20) Now add parameters for GCAP and GEOS-5 grids.  Remove references
!        to obsolete LGEOSCO and FULLCHEM Cpp switches (bmy, 6/24/05)
!  (21) Now add I1x1 and J1x1 parameters for data on the 1x1 GEOS
!        grid. (bmy, 10/24/05)
!  (22) Increase NNPAR to 52 (bmy, 12/6/05)
!  (23) Increase NNPAR to 54 (dkh, bmy, 5/22/06)
!  (24) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (25) Added variable tropopause case (bmy, phs, bdf, 8/21/06)
!  (26) Set LLTROP to level of first box entirely above 20km for GEOS-3 
!       and GEOS-4 (phS, 9/14/06)
!  (27) Bug fix: set LLTROP_FIX = LLPAR for GCAP (bmy, 11/29/06)
!  (28) Reset vertical coordinates for GEOS-5.  Also renamed GRID30LEV
!        to GRIDREDUCED (bmy, 4/3/07)
!  (29) New parameters for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (30) NEMPARA set to 12 to emit O3 and HNO3 (phs, 4/3/08)
!  (31) Add tracers to NNPAR = 73. (tmf, 1/7/09)
!  (32) NEMPARA set to 21 to emit new tracers for GLYX chemistry 
!       (tmf, ccc, 3/2/09)
!  (33) NEMPARB set to 3 to emit MBO, MONX (tmf, ccc, 3/2/09)
!  (34) Added EUROPE grid parameters (amv, 10/19/09)
!  18 Dec 2009 - Aaron van D - Added NESTED_EU grid parameters
!  18 Dec 2009 - R. Yantosca - Added ProTeX headers
!  25 Jan 2010 - R. Yantosca - Set NNPAR=320 for TOMAS simulations
!  25 Jan 2010 - R. Yantosca - Define TOMASBIN and TOMASSPEC for TOMAS sims
!  08 Feb 2010 - F. Paulot   - Increase NNPAR, NEMPARA and NEMPARB
!  08 Feb 2010 - F. Paulot   - Move MAXMEM and MAXPL from diag_pl_mod.
!  30 Nov 2010 - R. Yantosca - Increase LLTROP (from 38 to 40) for GEOS-5 and 
!                              MERRA for the full 72-layer grids (i.e. when
!                              the Cpp switch GRIDREDUCED is not set).
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!EOP
	END MODULE CMN_SIZE_MOD
!------------------------------------------------------------------------------
