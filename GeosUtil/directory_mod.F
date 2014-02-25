!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: directory_mod
!
! !DESCRIPTION: Module DIRECTORY\_MOD contains the directory path variables
!  used by GEOS-Chem.
!\\
!\\
! !INTERFACE: 
!
      MODULE DIRECTORY_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      CHARACTER(LEN=255) :: DATA_DIR      ! Main DAO met field directory
      CHARACTER(LEN=255) :: DATA_DIR_1x1  ! Root data dir for 1x1 emissions
      CHARACTER(LEN=255) :: GCAP_DIR      ! Subdir for GCAP met data
      CHARACTER(LEN=255) :: GEOS_1_DIR    ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: GEOS_S_DIR    ! !%%% OBSOLETE 
      CHARACTER(LEN=255) :: GEOS_4_DIR    ! Subdir for GEOS-4 met data
      CHARACTER(LEN=255) :: GEOS_5_DIR    ! Subdir for GEOS-5 met data
      CHARACTER(LEN=255) :: GEOS_FP_DIR   ! Subdir for GEOS-FP met data
      CHARACTER(LEN=255) :: MERRA_DIR     ! Subdir for MERRA  met data
      CHARACTER(LEN=255) :: TEMP_DIR      ! Temp dir for unzipping met data
      CHARACTER(LEN=255) :: RUN_DIR       ! Run directory for GEOS-Chem
      CHARACTER(LEN=255) :: OH_DIR        ! Dir w/ mean OH files
      CHARACTER(LEN=225) :: BC_DIR        ! Dir w/ mean BC files
      CHARACTER(LEN=255) :: OC_DIR        ! Dir w/ mean OC files
      CHARACTER(LEN=255) :: O3PL_DIR      ! Dir w/ archived O3 P/L rate files 
      CHARACTER(LEN=255) :: TPBC_DIR      ! Dir w/ TPCORE boundary conditions
      CHARACTER(LEN=255) :: TPBC_DIR_NA   ! Dir w/ TPCORE BC's for NA nest grid
      CHARACTER(LEN=255) :: TPBC_DIR_EU   ! Dir w/ TPCORE BC's for EU nest grid
      CHARACTER(LEN=255) :: TPBC_DIR_CH   ! Dir w/ TPCORE BC's for CH nest grid
      CHARACTER(LEN=255) :: TPBC_DIR_SE   ! Dir w/ TPCORE BC's for SEAC4RS grid
      CHARACTER(LEN=255) :: POP_EMISDIR   ! Dir w/ POPS emissions files
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS NOW OBSOLETE AND IS SLATED TO BE REMOVED!   %%%
!  %%%  WE NOW USE THE DIR NAMES FROM Input_Opt INSTEAD OF FROM HERE.    %%%
!  %%%     -- Bob Yantosca, 20 Aug 2013                                  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  25 May 2005 - R. Yantosca - Added variables GCAP_DIR and GEOS_5_DIR
!  24 Oct 2005 - R. Yantosca - Added DATA_DIR_1x1
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  18 Dec 2009 - Aaron van D - Added TPBC_DIR_NA, TPBC_DIR_EU, TPBC_DIR_CH
!  13 Aug 2010 - R. Yantosca - Added MERRA_DIR for MERRA met fields
!  20 Jan 2011 - CL Friedman - Added OC_DIR and BC_DIR
!  01 Feb 2012 - R. Yantosca - Added GEOS_57_DIR for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  29 May 2012 - S. Kim      - Added directory for SEAC4RS BCs
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57_DIR to GEOS_FP_DIR
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE DIRECTORY_MOD
!EOC
