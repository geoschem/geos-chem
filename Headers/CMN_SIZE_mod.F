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
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! Model top pressure (mb)
      REAL*8,     PARAMETER :: PTOP       = 0.01d0

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

#if   defined( APM )
      INTEGER,    PARAMETER :: APMSPEC    = 119

      INTEGER,    PARAMETER :: NAPMLVGs   = 2
      INTEGER,    PARAMETER :: NAPMSO4BIN = 40
      INTEGER,    PARAMETER :: NAPMLVA    = 1
      INTEGER,    PARAMETER :: NAPMCTBC   = 2
      INTEGER,    PARAMETER :: NAPMCTOC   = 2
      INTEGER,    PARAMETER :: NAPMCTDST  = 2
      INTEGER,    PARAMETER :: NAPMCTSEA  = 2
      INTEGER,    PARAMETER :: NAPMSEABIN = 20
      INTEGER,    PARAMETER :: NAPMDSTBIN = 15
      INTEGER,    PARAMETER :: NAPMBCBIN  = 15
      INTEGER,    PARAMETER :: NAPMOCBIN  = 15
      INTEGER,    PARAMETER :: NAPMAMINE  = 3
#endif

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

      SUBROUTINE Init_CMN_SIZE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root     ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt     ! Input Options object
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

      ! Return successfully
      RC = GC_SUCCESS

      END SUBROUTINE Cleanup_CMN_SIZE
!EOC
      END MODULE CMN_SIZE_MOD
