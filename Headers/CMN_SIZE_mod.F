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

#ifdef APM
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

#ifdef TOMAS
      INTEGER,    PARAMETER :: TOMASSPEC = 8

# if defined( TOMAS40 )
      INTEGER,    PARAMETER :: NDSTBIN   = 40
      INTEGER,    PARAMETER :: TOMASBIN  = 40
# elif defined( TOMAS15 )
      INTEGER,    PARAMETER :: NDSTBIN   = 15
      INTEGER,    PARAMETER :: TOMASBIN  = 15
# elif defined( TOMAS12 )
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

#ifdef BPCH_DIAG
      ! NJDAY is the array that is used to determine when it is time
      ! to end a GEOS-Chem run.  Moved this here from CMN_DIAG_mod.F,
      ! so that we can phase that module out. (bmy, 1/16/18)
      INTEGER               :: NJDAY(366)
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
!  See the Git history with the gitk browser!
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
!  See the Git history with the gitk browser!
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
