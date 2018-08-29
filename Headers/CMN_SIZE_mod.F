!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      USE PRECISION_MOD, ONLY : fpp => fp    ! For GEOS-Chem Precision (fpp)

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! DISIZE = size (in degrees) of a longitude grid box
      ! DJSIZE = size (in degrees) of a latitude  grid box
      !
      ! NOTE: DISIZE and DJSIZE are only needed for "GEOS-Chem Classic",
      ! as these are used to define the Delta-Lon and Delta-Lat values
      ! that need to get saved to the binary punch files.  If you are
      ! using "GEOS-Chem HP", then you can set this to a dummy value.
      ! (mlong, bmy, 5/22/15)
      !=================================================================
#if   defined( GRID4x5 ) 
      REAL(fpp), PARAMETER :: DISIZE = 5.0e+0_fpp
      REAL(fpp), PARAMETER :: DJSIZE = 4.0e+0_fpp
#elif defined( GRID2x25 )
      REAL(fpp), PARAMETER :: DISIZE = 2.5e+0_fpp 
      REAL(fpp), PARAMETER :: DJSIZE = 2.0e+0_fpp
#elif defined( GRID05x0625 )
      REAL(fpp), PARAMETER :: DISIZE = 0.625e+0_fpp
      REAL(fpp), PARAMETER :: DJSIZE = 0.5e+0_fpp
#elif defined( GRID025x03125)
      REAL(fpp), PARAMETER :: DISIZE = 0.3125e+0_fpp
      REAL(fpp), PARAMETER :: DJSIZE = 0.25e+0_fpp
#elif defined( EXTERNAL_GRID )
      REAL(fpp)            :: DISIZE
      REAL(fpp)            :: DJSIZE
#endif

      !=================================================================
      ! GRID SETTINGS #1: Necessary for the grid-independent GEOS-Chem
      !=================================================================
      INTEGER                :: I_LO         ! Minimum lon index on this CPU
      INTEGER                :: J_LO         ! Minimum lat index on this CPU
      INTEGER                :: I_HI         ! Maximum lon index on this CPU
      INTEGER                :: J_HI         ! Maximum lat index on this CPU
      INTEGER                :: IM_WORLD     ! # of lons in the global grid
      INTEGER                :: JM_WORLD     ! # of lats in the global grid
      INTEGER                :: LM_WORLD     ! # of levs in the global grid
      REAL(fpp), ALLOCATABLE :: DLON(:,:,:)  ! Array of delta-lon [degrees]
      REAL(fpp), ALLOCATABLE :: DLAT(:,:,:)  ! Array of delta-lat [degrees]

      !=================================================================
      ! GRID SETTINGS #2: Mostly historical declarations (keep for now)
      !
      ! IGLOB      = global longitude dimension
      ! JGLOB      = global latitude dimension
      ! LGLOB      = max number of sigma levels 
      ! IIPAR      = window longitude dimension
      ! JJPAR      = window latitude dimension
      ! LLPAR      = window vertical dimension
      ! LLTROP     = maximum number of tropospheric levels for variable
      !              tropopause
      ! LLTROP_FIX = number of tropospheric levels for offline simulations
      ! LLSTRAT    = maximum number of levels below stratopause
      ! LLCHEM     = maximum number of levels included in chemistry grid
      ! LLCHEM_FIX = number of chemgrid levels for offline simulations
      ! PTOP       = model top pressure (mb)
      !
      ! Most of the time, GEOS-CHEM is used for global simulations.
      ! In this case, then IIPAR=IGLOB, JJPAR=JGLOB, LLPAR=LGLOB.
      !
      ! For nested grids, then IIPAR<IGLOB, JJPAR<JGLOB, LLPAR<LGLOB.
      !=================================================================

      ! IIPAR and JJPAR are always set at runtime
      ! for both GEOS-Chem HP and GEOS-Chem Classic
      INTEGER               :: IIPAR
      INTEGER               :: JJPAR

#if defined ( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) || defined( ESMF_ )
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! We need to make all of these be variables rather than 
      ! parameters, so that they can be defined directly from the
      ! ESMF interface. (bmy, 10/11/12)
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB
      INTEGER               :: JGLOB
      INTEGER               :: LGLOB
      INTEGER               :: LLPAR                
      INTEGER               :: LLTROP_FIX
      INTEGER               :: LLTROP
      INTEGER               :: LLSTRAT
      REAL*8,     PARAMETER :: PTOP       = 0.01d0
#else

      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! For GEOS-Chem "Classic", we can preset the values of IGLOB,
      ! JGLOB, LGLOB, LLTROP, LLTROP_FIX, LLSTRAT, (bmy, 5/22/15)
      !-----------------------------------------------------------------

#if defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_CH )

      !-----------------------------------------------------------------
      ! GEOS-FP: Nested China Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 225
      INTEGER               :: JGLOB      = 161
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else 
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_NA)

      !-----------------------------------------------------------------
      ! GEOS-FP Nested NA Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 225
      INTEGER               :: JGLOB      = 202
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

      !Anna Prot added the following sizes for nesting EU, 7 May 2015
#elif defined( GEOS_FP ) && defined( GRID025x03125 ) && defined( NESTED_EU)

      !-----------------------------------------------------------------
      ! GEOS-FP Nested EU Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 177
      INTEGER               :: JGLOB      = 115
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( GEOS_FP ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-FP: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 144
      INTEGER               :: JGLOB      = 91
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else 
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( GEOS_FP ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-FP: 4 x 5
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 72
      INTEGER               :: JGLOB      = 46
      INTEGER               :: LGLOB      = 72         
#if   defined( GRIDREDUCED ) 
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels 
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( MERRA2 ) && defined( GRID05x0625 ) && defined( NESTED_AS )

      !-----------------------------------------------------------------
      ! MERRA2: Nested China Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 145
      INTEGER               :: JGLOB      = 133
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else 
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( MERRA2 ) && defined( GRID05x0625 ) && defined( NESTED_NA)

      !-----------------------------------------------------------------
      ! MERRA2: Nested NA Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 161
      INTEGER               :: JGLOB      = 121
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( MERRA2 ) && defined( GRID05x0625 ) && defined( NESTED_EU)

      !-----------------------------------------------------------------
      ! MERRA2: Nested EU Grid
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 129
      INTEGER               :: JGLOB      = 81
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( MERRA2 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! MERRA2: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 144
      INTEGER               :: JGLOB      = 91
      INTEGER               :: LGLOB      = 72
#if   defined( GRIDREDUCED )
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else 
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

#elif defined( MERRA2 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! MERRA2: 4 x 5
      !-----------------------------------------------------------------
      INTEGER               :: IGLOB      = 72
      INTEGER               :: JGLOB      = 46
      INTEGER               :: LGLOB      = 72         
#if   defined( GRIDREDUCED ) 
      INTEGER               :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER,    PARAMETER :: LLTROP     = 38  
      INTEGER,    PARAMETER :: LLSTRAT    = 44
#else
      INTEGER               :: LLPAR                  ! Full vertical grid
      INTEGER,    PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels 
      INTEGER,    PARAMETER :: LLTROP     = 40  
      INTEGER,    PARAMETER :: LLSTRAT    = 59
#endif
      REAL(fpp),  PARAMETER :: PTOP       = 0.01e+0_fpp

      ! This #endif below is to close (metfield,grid) choices 
      ! that specified by many #elif statements above
#endif

      ! This #endif below is to distinguish between classic and HP
#endif

      !-----------------------------------------------------------------
      ! For stratospheric chemistry
      !-----------------------------------------------------------------
      INTEGER               :: LLCHEM
      INTEGER               :: LLCHEM_FIX

      !=================================================================
      ! OTHER PARAMETERS 
      !=================================================================

      ! NSURFTYPE - Maximum number of surface types: 73 olson
      ! NTYPE     - Maximum number of veg types in a CTM grid box
      ! NPOLY     - Number of coefficients for polynomial fits
      INTEGER,    PARAMETER :: NSURFTYPE = 73
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_TYPE )
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! For testing the ESMF interface to GEOS-Chem with a grid that 
      ! is smaller than the usual 72x46, increase NTYPE (bmy, 12/4/12)
      !-----------------------------------------------------------------
      INTEGER,    PARAMETER :: NTYPE    = 50
#else
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! Current practice in the std GEOS-Chem is to set NTYPE to 25,
      ! which is large enough if using the Olson 2001 land map at
      ! 0.25 x 0.25 resolution. (bmy, 12/4/12)
      !-----------------------------------------------------------------
      INTEGER,    PARAMETER :: NTYPE    = 25
#endif
      INTEGER,    PARAMETER :: NPOLY    = 20

      ! MAXFAM -- Max number of families for prod and loss output
      INTEGER,    PARAMETER :: MAXFAM = 40

      ! MAXMEM is maximum number of families of prod and loss
      ! moved from input_mod and diag_pl_mod to here (hotp 7/31/09)
      ! MAXMEM also increased from 10 to 20 by FP
      ! MAXMEM increased from 20 to 22 for bromine (jpp, mpayer, 12/28/11)
      INTEGER,    PARAMETER :: MAXMEM = 22

      ! MAXPL increased from 100 to 500 and moved from diag_pl_mod
      ! to here by FP (hotp 7/31/09)
      INTEGER,    PARAMETER :: MAXPL = 500

      ! NDUST -- Number of FAST-J aerosol size bins (rvm, bmy, 11/15/01)
      INTEGER,    PARAMETER :: NDUST = 7

      ! NRHAER -- Number of aerosols undergoing hygroscopic growth
      INTEGER,    PARAMETER :: NRHAER = 5

      ! NSTRATAER -- Number of stratospheric aerosols (SDE 04/17/13)
      INTEGER,    PARAMETER :: NSTRATAER = 2

      ! NAER -- number of other aerosol categories (rvm, bmy, 2/27/02)
      ! Now set in Init_CMN_SIZE below (mps, 1/3/18)
      INTEGER               :: NAER

      ! NRH -- number of relative humidity bins (rvm, bmy, 2/27/02)
      INTEGER,    PARAMETER :: NRH = 5

#if   defined( TOMAS )
      INTEGER,    PARAMETER :: TOMASSPEC = 8

# if  defined( TOMAS40 )
      INTEGER,    PARAMETER :: NDSTBIN   = 40
      INTEGER,    PARAMETER :: TOMASBIN  = 40
# elif  defined( TOMAS15 )
      INTEGER,    PARAMETER :: NDSTBIN   = 15
      INTEGER,    PARAMETER :: TOMASBIN  = 15
# elif  defined( TOMAS12 )
      INTEGER,    PARAMETER :: NDSTBIN   = 12
      INTEGER,    PARAMETER :: TOMASBIN  = 12
# else
      ! NDSTBIN -- redimensioned for TOMAS (dwest, bmy, 2/1/10)
      INTEGER,    PARAMETER :: NDSTBIN   = 30
      INTEGER,    PARAMETER :: TOMASBIN  = 30 ! Number of TOMAS bins 
# endif

#else
      ! NDSTBIN -- number of dust size bins for transport (tdf, bmy, 3/31/04)
      INTEGER,    PARAMETER :: NDSTBIN   = 4
#endif

      ! NJDAY is the array that is used to determine when it is time
      ! to end a GEOS-Chem run.  Moved this here from CMN_DIAG_mod.F,
      ! so that we can phase that module out. (bmy, 1/16/18)
      INTEGER               :: NJDAY(366)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%% NOTE: THIS WAS MADE INTO A MODULE IN ORDER TO REMOVE COMMON BLOCKS  %%%
!  %%% WE WILL KEEP THIS FOR NOW.  EVENTUALLY WE MIGHT MIGRATE THESE DATA  %%%
!  %%% INTO A DERIVED TYPE OBJECT. (bmy, 12/3/12)                          %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  09 Mar 2011 - R. Yantosca - Updated NNPAR for APM (G. Luo)
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  27 Dec 2011 - M. Payer    - Updated NNPAR, NEMPARB, MAXMEM for bromine
!                              chemistry (J. Parrella) 
!  10 Feb 2012 - R. Yantosca - Added #if blocks for GEOS-5.7.x nested CH grid
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  27 Mar 2012 - R. Yantosca - Increase NTYPE from 15 to 25 for Olson 2001 map
!  22 Oct 2012 - M. Payer    - Increase NNPAR to 100 for tagged Hg simulation
!                              (E. Corbitt)
!  25 Oct 2012 - R. Yantosca - Now also set LLPAR, LLTROP, LLTROP_FIX to
!                              LGLOB for grid-independent simulation
!  19 Nov 2012 - R. Yantosca - Renamed to INIT_CMN_SIZE, to better follow
!                              adopted GEOS-Chem naming convention
!  27 Nov 2012 - R. Yantosca - Removed commented out code
!  05 Jun 2013 - R. Yantosca - Now define GEOS-5 0.25x0.3125 nested NA grid
!  13 Aug 2013 - M. Sulprizio- Increase NNPAR & NBIOMAX for updated SOA and
!                              SOA + semivolatile POA simulations (H. Pye) 
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  22 May 2015 - R. Yantosca - Updated comments, cosmetic changes
!  22 May 2015 - R. Yantosca - Removed variables made obsolete by HEMCO
!  28 Jan 2016 - M. Sulprizio- Remove NBIOMAX made obsolete by HEMCO
!  18 May 2016 - M. Sulprizio- Remove MAXIJ made obsolete by FlexChem
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  16 Mar 2017 - R. Yantosca - Remove obsolete NNPAR parameter
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  16 Jan 2018 - R. Yantosca - Moved NJDAY here from obsolete CMN_DIAG_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_Size
!
! !DESCRIPTION: Routine INIT\_CMN\_SIZE initializes the grid dimension values
!  in module CMN\_SIZE\_mod.F.
!\\
!\\
! !INTERFACE:

      SUBROUTINE Init_CMN_SIZE( am_I_Root,      Input_Opt,     RC,
     &                          value_I_LO,     value_J_LO,
     &                          value_I_HI,     value_J_HI,
     &                          value_IM,       value_JM,  
     &                          value_LM,       value_IM_WORLD, 
     &                          value_JM_WORLD, value_LM_WORLD,
     &                          value_LLTROP,   value_LLSTRAT  )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
! 
      LOGICAL,        INTENT(IN) :: am_I_Root  ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
      INTEGER,        OPTIONAL   :: value_I_LO    ! Lower lon index on this CPU
      INTEGER,        OPTIONAL   :: value_J_LO    ! Lower lat index on this CPU
      INTEGER,        OPTIONAL   :: value_I_HI    ! Upper lon index on this CPU
      INTEGER,        OPTIONAL   :: value_J_HI    ! Upper lat index on this CPU
      INTEGER,        OPTIONAL   :: value_IM      ! # of lons   on this CPU
      INTEGER,        OPTIONAL   :: value_JM      ! # of lats   on this CPU
      INTEGER,        OPTIONAL   :: value_LM      ! # of levels on this CPU
      INTEGER,        OPTIONAL   :: value_IM_WORLD! # of lons in the global grid
      INTEGER,        OPTIONAL   :: value_JM_WORLD! # of lats in the global grid
      INTEGER,        OPTIONAL   :: value_LM_WORLD! # of levs in the global grid
      INTEGER,        OPTIONAL   :: value_LLTROP  ! # of levs in troposphere
      INTEGER,        OPTIONAL   :: value_LLSTRAT ! # of levs in stratosphere
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_Dimensions
!  03 Dec 2012 - R. Yantosca - Rewritten for clarity.  Also pass optional 
!                              arguments to carry values from ESMF environment
!  03 Dec 2012 - R. Yantosca - Now allocate DLON, DLAT arrays here.
!  04 Aug 2015 - C. Keller   - Added value_LLTROP and value_LLSTRAT. 
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  05 Jul 2018 - H. Lin      - Initialize DLAT to 0 (previously missing)
!EOP
!------------------------------------------------------------------------------
!BOC      

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING ) || defined( ESMF_ )
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! If we are connecting to the ESMF interface, then define all
      ! size parameters with values passed from the ESMF environment
      ! via the optional arguments.  
      !
      ! We have introduced variables for the lower and upper longitude
      ! and latitude indices on the current CPU (I_LO, J_LO, I_HI, J_HI)
      ! as well as the # of longitudes, latitudes, and levels in the 
      ! entire global grid (IM_WORLD, JM_WORLD, LM_WORLD).  Knowing 
      ! these values will facilitate computation of certain quantities
      ! when running with  MPI parallelization in the ESMF environment.  
      !
      ! NOTE: We don't need to define min & max level variables L_LO, 
      ! L_HI, because we will always be decomposing the data with MPI 
      ! in the  horizontal.  Thus L_LO=1 and L_HI=LLPAR.
      !
      ! We must also fix LLTROP = LLPAR and LLTROP_FIX = LLPAR, 
      ! because we can only have one vertical dimension within the
      ! ESMF environment. (bmy, 12/3/12)
      !-----------------------------------------------------------------

      !=================================================================
      ! Set values for the upper & lower longitude and latitude
      ! indices that will be processed by the current CPU.  
      !=================================================================
      I_LO       = value_I_LO       ! Min longitude index on this CPU
      J_LO       = value_J_LO       ! Min latitude  index on this CPU
      I_HI       = value_I_HI       ! Max longitude index on this CPU
      J_HI       = value_J_HI       ! Max latitude  index on this CPU

      !=================================================================
      ! Set values for the number of longitudes, latitudes, and 
      ! levels that will be processed by the current CPU. 
      !=================================================================

      IIPAR      = value_IM         ! # of longitudes on this CPU
      JJPAR      = value_JM         ! # of latitudes  on this CPU
      LLPAR      = value_LM         ! # of levels     on this CPU

      IGLOB      = IIPAR            ! For backwards compatibility, set = IIPAR
      JGLOB      = JJPAR            ! For backwards compatibility, set = JJPAR
      LGLOB      = LLPAR            ! For backwards compatibility, set = LLPAR

      IF ( PRESENT(value_LLTROP) ) THEN
         LLTROP  = value_LLTROP
      ELSE
         LLTROP  = 40
      ENDIF
      LLTROP_FIX = LLTROP

      IF ( PRESENT(value_LLSTRAT) ) THEN
         LLSTRAT = value_LLSTRAT
      ELSE
         LLSTRAT = 59
      ENDIF

      IF ( Input_Opt%LUCX ) THEN
         LLCHEM     = LLSTRAT
         LLCHEM_FIX = LLSTRAT
      ELSE
         LLCHEM     = LLTROP
         LLCHEM_FIX = LLTROP_FIX
      ENDIF

      ! For GEOS-Chem HP, the DISIZE and DJSIZE variables can just be set 
      ! to dummy variables.  These are not used in the ESMF environment.
      ! (mlong, bmy, 5/22/15)
      DISIZE     = 1
      DJSIZE     = 1

      !=================================================================
      ! Set values for the number of longitudes, latitudes, and levels 
      ! in the whole global grid.  We need these values to compute 
      ! certain quantities (i.e. grid box areas) when connecting to 
      ! the GEOS-5 GCM via the ESMF environment.
      !=================================================================
      IM_WORLD   = value_IM_WORLD   ! # of longitudes in the global grid
      JM_WORLD   = value_JM_WORLD   ! # of latitudes  in the global grid
      LM_WORLD   = value_LM_WORLD   ! # of levels     in the global grid

#else
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! Current practice in the standard GEOS-Chem is to set the IIPAR,
      ! JJPAR, and LLPAR variables from the IGLOB, JGLOB, and LGLOB
      ! parameters defined at the top of the module, above the CONTAINS
      ! statement.  The values of IGLOB, JGLOB, LGLOB are pegged to the 
      ! particular met field that is being used.
      !
      ! Also, only the GMAO met fields (GEOS-4, GEOS-5, GEOS-5.7, and
      ! MERRA) have a reduced-grid option.  Therefore, set LLPAR to
      ! LGLOB for GCAP here.  (bmy, 12/3/12)
      !-----------------------------------------------------------------

      !=================================================================
      ! Set values for for IIPAR, JJPAR, LLPAR variables from the 
      ! parameter settings IGLOB, JGLOB, LGLOB defined above
      !=================================================================
      IIPAR      = IGLOB            ! Set IIPAR to IGLOB parameter
      JJPAR      = JGLOB            ! Set JJPAR to JGLOB parameter
#if !defined( GRIDREDUCED )
      LLPAR      = LGLOB            ! If GRIDREDUCED is chosen, then LLPAR
                                    !  and LLTROP are defined above, so no
                                    !  need to reset them here.
#endif

      IF ( Input_Opt%LUCX ) THEN
         LLCHEM     = LLSTRAT
         LLCHEM_FIX = LLSTRAT
      ELSE
         LLCHEM     = LLTROP
         LLCHEM_FIX = LLTROP_FIX
      ENDIF

      !=================================================================
      ! In the traditional GEOS-Chem, we use the OpenMP (or loop-level)
      ! paralellization instead of MPI parallelization.  The following
      ! variables may be therefore be set to global values.
      !=================================================================
      I_LO       = 1                ! Min longitude index on this CPU
      J_LO       = 1                ! Min latitude  index on this CPU
      I_HI       = IIPAR            ! Max longitude index on this CPU
      J_HI       = JJPAR            ! Max latitude  index on this CPU
      IM_WORLD   = IIPAR            ! # of longitudes in the global grid
      JM_WORLD   = JJPAR            ! # of latitudes  in the global grid
      LM_WORLD   = LLPAR            ! # of levels     in the global grid

#endif   

      !=================================================================
      ! Set values for former variables declared as parameters
      !=================================================================

      ! Number of aerosol categories
      IF ( Input_Opt%LUCX ) THEN
         ! UCX-based mechanisms include stratospheric aerosols
         NAER = NRHAER + NSTRATAER
      ELSE
         NAER = NRHAER
      ENDIF

      !================================================================
      ! Allocate arrays here
      !================================================================

      ! Delta-longitudes [degrees]
      ALLOCATE( DLON( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DLON = 0e+0_fpp

      ! Delta-latitudes [degrees]
      ALLOCATE( DLAT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DLAT = 0e+0_fpp

      ! Return w/ success
      RC = GC_SUCCESS
      
      END SUBROUTINE Init_CMN_SIZE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_Size
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_SIZE deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_SIZE( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!   3 Dec 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Deallocate arrays
      IF ( ALLOCATED( DLON ) ) DEALLOCATE( DLON )
      IF ( ALLOCATED( DLAT ) ) DEALLOCATE( DLAT )

      ! Return successfully
      RC = GC_SUCCESS

      END SUBROUTINE Cleanup_CMN_SIZE
!EOC
      END MODULE CMN_SIZE_MOD
