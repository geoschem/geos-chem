! $Id: directory_mod.f,v 1.1 2004/09/21 18:04:12 bmy Exp $
      MODULE DIRECTORY_MOD
!
!******************************************************************************
!  Module DIRECTORY_MOD contains the directory path variables used by 
!  GEOS-CHEM. (bmy, 7/20/04)
!     
!  Module Variables:
!  ============================================================================
!  (1 ) DATA_DIR   (CHAR*255) : Main DAO met field directory
!  (2 ) GEOS_1_DIR (CHAR*255) : Subdir where GEOS-1     met data are stored
!  (3 ) GEOS_S_DIR (CHAR*255) : Subdir where GEOS-STRAT met data are stored
!  (4 ) GEOS_3_DIR (CHAR*255) : Subdir where GEOS-3     met data are stored
!  (5 ) GEOS_4_DIR (CHAR*255) : Subdir where GEOS-4     met data are stored
!  (6 ) TEMP_DIR   (CHAR*255) : Temporary directory for unzipping met dat
!  (7 ) RUN_DIR    (CHAR*255) : Run directory for GEOS-CHEM
!  (8 ) OH_DIR     (CHAR*255) : Dir containing OH files are stored
!  (9 ) O3PL_DIR   (CHAR*255) : Dir containing archived O3 P/L rate files 
!  (10) TPBC_DIR   (CHAR*255) : TPCORE boundary conditions dir (nested grid)
!
!  NOTES:
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER(LEN=255) :: DATA_DIR
      CHARACTER(LEN=255) :: GEOS_1_DIR
      CHARACTER(LEN=255) :: GEOS_S_DIR
      CHARACTER(LEN=255) :: GEOS_3_DIR
      CHARACTER(LEN=255) :: GEOS_4_DIR
      CHARACTER(LEN=255) :: TEMP_DIR
      CHARACTER(LEN=255) :: RUN_DIR
      CHARACTER(LEN=255) :: OH_DIR
      CHARACTER(LEN=255) :: O3PL_DIR
      CHARACTER(LEN=255) :: TPBC_DIR

      ! End of module
      END MODULE DIRECTORY_MOD
