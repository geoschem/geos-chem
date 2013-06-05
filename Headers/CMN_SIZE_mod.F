!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_SIZE
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
      IMPLICIT NONE
      PUBLIC

#     include "define.h"
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! DISIZE = size (in degrees) of a longitude grid box
      ! DJSIZE = size (in degrees) of a latitude  grid box
      !=================================================================
#if   defined( GRID4x5 ) 
      REAL*8, PARAMETER :: DISIZE = 5.0d0
      REAL*8, PARAMETER :: DJSIZE = 4.0d0
#elif defined( GRID2x25 )
      REAL*8, PARAMETER :: DISIZE = 2.5d0 
      REAL*8, PARAMETER :: DJSIZE = 2.0d0
#elif defined( GRID1x125 )
      REAL*8, PARAMETER :: DISIZE = 1.25d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0
#elif defined( GRID1x1 )
      REAL*8, PARAMETER :: DISIZE = 1.0d0 
      REAL*8, PARAMETER :: DJSIZE = 1.0d0
#elif defined( GRID05x0666 )
      REAL*8, PARAMETER :: DISIZE = 2d0/3d0
      REAL*8, PARAMETER :: DJSIZE = 0.5d0
#elif defined( GRID025x03125)
      REAL*8, PARAMETER :: DISIZE = 0.3125d0
      REAL*8, PARAMETER :: DJSIZE = 0.25d0
#elif defined( EXTERNAL_GRID )
      REAL*8            :: DISIZE
      REAL*8            :: DJSIZE
#endif

      !=================================================================
      ! GRID SETTINGS #1: Necessary for the grid-independent GEOS-Chem
      !=================================================================
      INTEGER             :: I_LO         ! Minimum lon index on this CPU
      INTEGER             :: J_LO         ! Minimum lat index on this CPU
      INTEGER             :: I_HI         ! Maximum lon index on this CPU
      INTEGER             :: J_HI         ! Maximum lat index on this CPU
      INTEGER             :: IM_WORLD     ! # of lons in the whole global grid
      INTEGER             :: JM_WORLD     ! # of lats in the whole global grid
      INTEGER             :: LM_WORLD     ! # of levs in the whole global grid
      REAL*8, ALLOCATABLE :: DLON(:,:,:)  ! Array of delta-longitude [degrees]
      REAL*8, ALLOCATABLE :: DLAT(:,:,:)  ! Array of delta-latitude  [degrees]

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
      ! PTOP       = model top pressure (mb)
      !
      ! Most of the time, GEOS-CHEM is used for global simulations.
      ! In this case, then IIPAR=IGLOB, JJPAR=JGLOB, LLPAR=LGLOB.
      !
      ! For nested grids, then IIPAR<IGLOB, JJPAR<JGLOB, LLPAR<LGLOB.
      !=================================================================
#if   defined( GCAP ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GCAP: 4 x 5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 72
      INTEGER            :: JGLOB      = 45
      INTEGER            :: LGLOB      = 23
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
      INTEGER            :: LLPAR
      INTEGER, PARAMETER :: LLTROP     = 12   
      INTEGER, PARAMETER :: LLTROP_FIX = LLTROP     
      REAL*8,  PARAMETER :: PTOP       = 0.002d0

#elif defined( GEOS_4 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-4: 4 x 5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 72
      INTEGER            :: JGLOB      = 46
      INTEGER            :: LGLOB      = 55           
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED ) 
      INTEGER            :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-4: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 144
      INTEGER            :: JGLOB      = 91
      INTEGER            :: LGLOB      = 55
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_4 ) && defined( GRID1x125 )

      !-----------------------------------------------------------------
      ! GEOS-4: 1 x 1.2.5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 288
      INTEGER            :: JGLOB      = 181
      INTEGER            :: LGLOB      = 55
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 30        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP     = 22        !  -- 30 levels
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP     = 23        !  -- 55 levels
#endif
      INTEGER, PARAMETER :: LLTROP_FIX = 17
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-5: 4 x 5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 72
      INTEGER            :: JGLOB      = 46
      INTEGER            :: LGLOB      = 72         
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED ) 
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-5: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 144
      INTEGER            :: JGLOB      = 91
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID1x125 )
      !-----------------------------------------------------------------
      ! GEOS-5: 1 x 1.25
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 288
      INTEGER            :: JGLOB      = 181
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels    
      INTEGER, PARAMETER :: LLTROP     = 40
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_5 ) && defined( GRID05x0666 )

      !-----------------------------------------------------------------
      ! GEOS-5: 0.5 x 0.666
      !-----------------------------------------------------------------
#if   defined( NESTED_CH )
      INTEGER            :: IGLOB      = 121       ! NESTED CHINA   0.5x0.666
      INTEGER            :: JGLOB      = 133
      INTEGER            :: LGLOB      = 72
#elif defined( NESTED_NA )
      INTEGER            :: IGLOB      = 151       ! NESTED N.AMER. 0.5x0.666
      INTEGER            :: JGLOB      = 121
      INTEGER            :: LGLOB      = 72
#elif defined( NESTED_EU )
      INTEGER            :: IGLOB      = 121       ! NESTED EUROPE  0.5x0.666
      INTEGER            :: JGLOB      = 81
      INTEGER            :: LGLOB      = 72
#endif
      INTEGER            :: IIPAR  
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( MERRA ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! MERRA: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 144
      INTEGER            :: JGLOB      = 91
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( MERRA ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! MERRA: 4 x 5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 72
      INTEGER            :: JGLOB      = 46
      INTEGER            :: LGLOB      = 72         
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED ) 
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels 
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_57 ) && defined( GRID025x03125 ) && defined( SEAC4RS )

      !-----------------------------------------------------------------
      ! GEOS-5.7.x: SEA4CRS Grid
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 177
      INTEGER            :: JGLOB      = 161
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_57 ) && defined( GRID025x03125 ) && defined( NESTED_CH )

      !-----------------------------------------------------------------
      ! GEOS-5.7.x: Nested China Grid
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 225
      INTEGER            :: JGLOB      = 161
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_57 ) && defined( GRID025x03125 ) && defined( NESTED_NA)

      !-----------------------------------------------------------------
      ! GEOS-5.7.x: Nested NA Grid
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 225
      INTEGER            :: JGLOB      = 202
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_57 ) && defined( GRID2x25 )

      !-----------------------------------------------------------------
      ! GEOS-5.7.x: 2 x 2.5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 144
      INTEGER            :: JGLOB      = 91
      INTEGER            :: LGLOB      = 72
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED )
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else 
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined( GEOS_57 ) && defined( GRID4x5 )

      !-----------------------------------------------------------------
      ! GEOS-5.7.x: 4 x 5
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB      = 72
      INTEGER            :: JGLOB      = 46
      INTEGER            :: LGLOB      = 72         
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
#if   defined( GRIDREDUCED ) 
      INTEGER            :: LLPAR      = 47        ! Reduced vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 38        !  -- 47 levels
      INTEGER, PARAMETER :: LLTROP     = 38  
#else
      INTEGER            :: LLPAR                  ! Full vertical grid
      INTEGER, PARAMETER :: LLTROP_FIX = 40        !  -- 72 levels 
      INTEGER, PARAMETER :: LLTROP     = 40  
#endif
      REAL*8,  PARAMETER :: PTOP       = 0.01d0

#elif defined ( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! We need to make all of these be variables rather than 
      ! parameters, so that they can be defined directly from the
      ! ESMF interface. (bmy, 10/11/12)
      !-----------------------------------------------------------------
      INTEGER            :: IGLOB
      INTEGER            :: JGLOB
      INTEGER            :: LGLOB
      INTEGER            :: IIPAR
      INTEGER            :: JJPAR
      INTEGER            :: LLPAR                
      INTEGER            :: LLTROP_FIX
      INTEGER            :: LLTROP
      REAL*8,  PARAMETER :: PTOP       = 0.01d0
#endif

      !-----------------------------------------------------------------
      ! For GEOS 1x1 files
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: I1x1       = 360
      INTEGER, PARAMETER :: J1x1       = 181

      !-----------------------------------------------------------------
      ! For GEOS 05x0666 files
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: I05x0666   = 540 
      INTEGER, PARAMETER :: J05x0666   = 361

      !-----------------------------------------------------------------
      ! For GFED3
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: IGFED3     = 720 
      INTEGER, PARAMETER :: JGFED3     = 360

      !=================================================================
      ! TRACER & EMISSION SPECIES PARAMETERS
      !
      ! NNPAR   = max number of tracers
      ! NEMPARA = max number of anthropogenic emission species
      ! NEMPARB = max number of biogenic      emission species
      !=================================================================
      ! increase NNPAR and NEMPARA an extra amount (hotp 7/31/09)
#if   defined( TOMAS )
# if  defined( TOMAS40 )
      INTEGER, PARAMETER :: NNPAR   = 410   ! For TOMAS40 (jeff 2/25/10)
# else
      INTEGER, PARAMETER :: NNPAR   = 320   ! For TOMAS (win, bmy, 1/25/10)
# endif
#elif defined( APM )
      INTEGER, PARAMETER :: NNPAR   = 154   ! For APM (G. Luo, 3/8/11)
#else
      INTEGER, PARAMETER :: NNPAR   = 100   ! For non-TOMAS simulations
#endif

      ! Nempara increased to 26. (fp, 2/8/10)
      ! new emissions HNO3 and O3 (phs)
      !  Add non-biogenic emission species:
      !  BENZ, TOLU, XYLE, C2H2, C2H4, GLYX, MGLY, GLYC, HAC. (tmf, 1/7/09)
      !INTEGER, PARAMETER :: NEMPARA = 21 
      !  Add RCHO, HCOOH, ACTA
      !INTEGER, PARAMETER :: NEMPARA = 26 
      INTEGER, PARAMETER :: NEMPARA = 29
      INTEGER, PARAMETER :: NEMPARB = 20

      !=================================================================
      ! OTHER PARAMETERS 
      !=================================================================

      ! NVEGTYPE - Maximum number of surface types: 74 olson
      ! NTYPE    - Maximum number of veg types in a CTM grid box
      ! NPOLY    - Number of coefficients for polynomial fits
      INTEGER, PARAMETER :: NVEGTYPE = 74
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_TYPE )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! For testing the ESMF interface to GEOS-Chem with a grid that 
      ! is smaller than the usual 72x46, increase NTYPE (bmy, 12/4/12)
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: NTYPE    = 50
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the std GEOS-Chem is to set NTYPE to 25,
      ! which is large enough if using the Olson 2001 land map at
      ! 0.25 x 0.25 resolution. (bmy, 12/4/12)
      !-----------------------------------------------------------------
      INTEGER, PARAMETER :: NTYPE    = 25
#endif
      INTEGER, PARAMETER :: NPOLY    = 20

      ! NNSTA = max number of time series stations (in inptr.ctm)
      INTEGER, PARAMETER :: NNSTA = 800

      ! MAXIJ - Maximum number of 1st level grid boxes
      INTEGER :: MAXIJ

      ! LLCONVM - Max number of layers for convection
      INTEGER :: LLCONVM

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
      ! MAXMEM increased from 20 to 22 for bromine (jpp, mpayer, 12/28/11)
      INTEGER, PARAMETER :: MAXMEM = 22

      ! MAXPL increased from 100 to 500 and moved from diag_pl_mod
      ! to here by FP (hotp 7/31/09)
      INTEGER, PARAMETER :: MAXPL = 500


      ! NDUST -- Number of FAST-J aerosol size bins (rvm, bmy, 11/15/01)
      INTEGER, PARAMETER :: NDUST = 7

      ! NAER -- number of other aerosol categories (rvm, bmy, 2/27/02)
      INTEGER, PARAMETER :: NAER = 5

      ! NRH -- number of relative humidity bins (rvm, bmy, 2/27/02)
      INTEGER, PARAMETER :: NRH = 5

      !NBIOMAX -- biomass burning
      ! increase NBIOMAX to 20 (hotp 7/31/09)
      ! increase NBIOMAX to 24 for dicarbonyls (ccc, 2/02/10)
      ! increase NBIOMAX to 25 fpr CH4 (kjw)
      INTEGER, PARAMETER :: NBIOMAX = 25


#if   defined( TOMAS )
      INTEGER, PARAMETER :: TOMASSPEC = 8

# if  defined( TOMAS40 )
      INTEGER, PARAMETER :: NDSTBIN   = 40
      INTEGER, PARAMETER :: TOMASBIN  = 40
# else
      ! NDSTBIN -- redimensioned for TOMAS (dwest, bmy, 2/1/10)
      INTEGER, PARAMETER :: NDSTBIN   = 30
      INTEGER, PARAMETER :: TOMASBIN  = 30 ! Number of TOMAS bins 
# endif

#else
      ! NDSTBIN -- number of dust size bins for transport (tdf, bmy, 3/31/04)
      INTEGER, PARAMETER :: NDSTBIN   = 4
#endif
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
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn_size
!
! !DESCRIPTION: Routine INIT\_CMN\_SIZE initializes the grid dimension values
!  in module CMN_SIZE_mod.F.
!\\
!\\
! !INTERFACE:

      SUBROUTINE Init_CMN_SIZE( am_I_Root,      RC,                    
     &                          value_I_LO,     value_J_LO,
     &                          value_I_HI,     value_J_HI,
     &                          value_IM,       value_JM,  
     &                          value_LM,       value_IM_WORLD, 
     &                          value_JM_WORLD, value_LM_WORLD  )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
! 
      LOGICAL, INTENT(IN)  :: am_I_Root        ! Are we on the root CPU?
      INTEGER, OPTIONAL    :: value_I_LO       ! Lower lon index on this CPU
      INTEGER, OPTIONAL    :: value_J_LO       ! Lower lat index on this CPU
      INTEGER, OPTIONAL    :: value_I_HI       ! Upper lon index on this CPU
      INTEGER, OPTIONAL    :: value_J_HI       ! Upper lat index on this CPU
      INTEGER, OPTIONAL    :: value_IM         ! # of lons   on this CPU
      INTEGER, OPTIONAL    :: value_JM         ! # of lats   on this CPU
      INTEGER, OPTIONAL    :: value_LM         ! # of levels on this CPU
      INTEGER, OPTIONAL    :: value_IM_WORLD   ! # of lons in the global grid
      INTEGER, OPTIONAL    :: value_JM_WORLD   ! # of lats in the global grid
      INTEGER, OPTIONAL    :: value_LM_WORLD   ! # of levs in the global grid
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
!EOP
!------------------------------------------------------------------------------
!BOC      

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
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
      LLTROP     = LLPAR            ! For backwards compatibility, set = LLPAR
      LLTROP_FIX = LLPAR            ! For backwards compatibility, set = LLPAR

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
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
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

#if defined( GCAP )
      LLPAR      = LGLOB            ! GCAP does not have a reduced grid option
#else 

#if !defined( GRIDREDUCED )
      LLPAR      = LGLOB            ! If GRIDREDUCED is chosen, then LLPAR
                                    !  and LLTROP are defined above, so no
                                    !  need to reset them here.
#endif
#endif

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

      ! Set values for former variables declared as parameters
      MAXIJ   = IIPAR * JJPAR      
      LLCONVM = LLPAR - 1

      !================================================================
      ! Allocate arrays here
      !================================================================

      ! Delta-longitudes [degrees]
      ALLOCATE( DLON( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      DLON = 0d0

      ! Delta-latitudes [degrees]
      ALLOCATE( DLAT( IIPAR, JJPAR, LLPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      DLON = 0d0

      ! Return w/ success
      RC = GIGC_SUCCESS
      
      END SUBROUTINE Init_CMN_SIZE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_size
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
      USE GIGC_ErrCode_Mod
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
      RC = GIGC_SUCCESS

      END SUBROUTINE Cleanup_CMN_SIZE
!EOC
      END MODULE CMN_SIZE_MOD
