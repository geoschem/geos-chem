!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: smv_dimension_mod
!
! !DESCRIPTION: This include file contains the various placeholder parameters
!  that are required to replace references to GEOS-Chem grid parameters.  
!  This is necessary because several quantities in the FAST-J and SMVGEAR
!  codes are contained in common blocks, and we need to have these parameters
!  for sizing those arrays properly.
!\\
!\\
! !INTERFACE: 
!
      MODULE SMV_DIMENSION_MOD
!
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! Locally defined replacement for GEOS-Chem parameter "LLPAR"
      ! This is used to set the "NB" value in jv_cmn_mod.F
#if defined( GCAP )

      !----------------------------------------
      ! GCAP vertical grid
      !----------------------------------------
      INTEGER, PARAMETER :: MAX_COLUMN   = 23

#elif defined( GEOS_4 ) 

      !----------------------------------------
      ! GEOS-4 vertical grid
      !----------------------------------------
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: MAX_COLUMN   = 30   ! Reduced GEOS-4 vertical grid
#else
      INTEGER, PARAMETER :: MAX_COLUMN   = 55   ! Full GEOS-4 vertical grid
#endif

#else

      !----------------------------------------
      ! GEOS-5, MERRA, GEOS-FP vertical grids
      !----------------------------------------
#if   defined( GRIDREDUCED )
      INTEGER, PARAMETER :: MAX_COLUMN   = 47   ! Reduced GEOS-5 vertical grid
#else
      INTEGER, PARAMETER :: MAX_COLUMN   = 72   ! Full GEOS-5 vertical grid
#endif

#endif

      ! Locally defined replacement for GEOS-Chem parameters "IIPAR" & "JJPAR"
      INTEGER, PARAMETER :: MAX_LAT      = 47   !
      INTEGER, PARAMETER :: MAX_LON      = 47   !

      ! Locally defined replacement for GEOS-Chem parameter "NNPAR"
      INTEGER, PARAMETER :: MAX_TRACERS = 100

      ! Locally defined replacement for "comode.h" parameter "IGAS"
      INTEGER, PARAMETER :: MAX_SPECIES = 125
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WAS ORIGINALLY DEVELOPED FOR THE COLUMN CODE.  %%%
!  %%%  THE PARAMETERS HERE CAN BE COMBINED INTO A SINGLE MODULE WITH    %%%
!  %%%  DECLARED SIZE PARAMETERS (IN THE NEAR FUTURE).                   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY: 
!  24 Mar 2009 - R. Yantosca - Initial version
!  16 Apr 2010 - R. Yantosca - Added MAX_SPECIES = 125
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!  10 Aug 2012 - R. Yantosca - Now define MAX_COLUMN=47 for GRIDREDUCED or
!                              =72 otherwise.  Kludge for testing DEVEL code.
!  02 Jul 2013 - R. Yantosca - Set MAX_COLUMN properly for GCAP grid
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  28 Oct 2013 - R. Yantosca - Bug fix: add GEOS-4 info
!EOP
      END MODULE SMV_DIMENSION_MOD
!------------------------------------------------------------------------------
