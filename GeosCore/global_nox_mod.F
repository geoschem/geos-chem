!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_NOx_mod
!
! !DESCRIPTION: Module GLOBAL\_NOx\_MOD contains variables and routines for 
!  reading the global monthly mean NOx concentration from disk.
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_NOX_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean BNOX field
      REAL*8, PUBLIC, ALLOCATABLE :: BNOX(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_GLOBAL_NOx
      PUBLIC :: GET_GLOBAL_NOx
      PUBLIC :: INIT_GLOBAL_NOx
!
! !REVISION HISTORY:
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Updated comments, made cosmetic changes (bmy, 6/13/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Now regrid BNOX array from 48L to 30L for GEOS-3 if necessary.
!        (bmy, 1/14/02)
!  (4 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (6 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (7 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (8 ) Cosmetic changes to improve output (bmy, 3/27/03)
!  (9 ) Now references "directory_mod.f" and "unix_cmds_mod.f" (bmy, 7/20/04)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: get_global_nox
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_NOx reads global NOx from binary 
!  punch files from a full chemistry run.  This NOx data is needed to 
!  calculate the CO yield from isoprene oxidation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_NOx( THISMONTH, Input_Opt )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT
      USE BPCH2_MOD,          ONLY : GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRANSFER_MOD,       ONLY : TRANSFER_3D
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Now use version of GET_TAU0 with 3 arguments.  Now call READ_BPCH2 
!        with IIPAR,JJPAR,LGLOB.  Call TRANSFER_3D to cast from REAL*4 to 
!        REAL*8 and to regrid to 30 levels for GEOS-3 (if necessary).  ARRAY 
!        should now be of size (IIPAR,JJPAR,LGLOB). (bmy, 1/14/02)
!  (2 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas.  Also make sure
!        to define FILENAME before printing it (bmy, 4/28/03)
!  (4 ) Now references TEMP_DIR, DATA_DIR from "directory_mod.f".  Also
!        references Unix unzipping commands from "unix_cmds_mod.f".
!        (bmy, 7/20/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  11 Apr 2013 - R. Yantosca - Now pass directory fields etc with Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L
      REAL*4             :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255) :: FIELD_DIR, RGNAME, TEMPO, CHAROP
      CHARACTER(LEN=3)   :: BMONTH(12) = (/ 'jan', 'feb', 'mar', 
     &                                      'apr', 'may', 'jun', 
     &                                      'jul', 'aug', 'sep', 
     &                                      'oct', 'nov', 'dec' /)

      ! First time flag
      LOGICAL, SAVE      :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_NOX begins here!
      !=================================================================

      ! Allocate NOx array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_NOx
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Construct file names and uncompress commands
      !=================================================================

      ! Name of unzipped file in TEMP_DIR
      TEMPO = 'tempo'
      
      ! Directory where the NOx files reside
      FIELD_DIR = '/data/ctm/GEOS_MEAN/OHparam/'

      ! Name of the zipped punch file w/ NOx in FIELD_DIR
      RGNAME = TRIM( FIELD_DIR )   // 'ctm.bpch.'                  // 
     &         BMONTH( THISMONTH ) // '.'                          // 
     &         GET_NAME_EXT()      // TRIM( Input_Opt%ZIP_SUFFIX )

      ! Construct the command to unzip the file & copy to TEMP_DIR
      CHAROP = TRIM( Input_Opt%UNZIP_CMD ) // ' ' //
     &         TRIM( RGNAME              ) // 
     &         TRIM( Input_Opt%REDIRECT  ) // ' ' //
     &         TRIM( Input_Opt%TEMP_DIR  ) //
     &         TRIM( TEMPO               )

      ! Uncompress the file and store in TEMP_DIR
      CALL SYSTEM( TRIM( CHAROP ) )

      !=================================================================
      ! Read NOx data from the punch file
      !=================================================================

      ! Read 1997 NOx data for Jan-Aug; Read 1996 NOx data for Sep-Dec 
      ! This avoids the 1997 El Nino signal in the NOx data
      IF ( THISMONTH >= 9 ) THEN
         XTAU = GET_TAU0( THISMONTH, 1, 1996 )
      ELSE
         XTAU = GET_TAU0( THISMONTH, 1, 1997 )
      ENDIF

      ! Name of unzipped file in TEMP_DIR
      FILENAME = TRIM( Input_Opt%TEMP_DIR ) // TRIM( TEMPO )

      ! Echo info
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_NOX: Reading NOX from: ', a )
      
      ! Read NOX data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,     XTAU,  
     &                 IIPAR,    JJPAR,      LGLOB, ARRAY )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_3D( ARRAY, BNOX )

      END SUBROUTINE GET_GLOBAL_NOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_NOx
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_NOx allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_NOX
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) BNOX now needs to be sized (IIPAR,JJPAR,LLPAR) (bmy, 1/14/02)
!  (2 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      ! Allocate BNOX array
      ALLOCATE( BNOX( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNOX' )

      ! Zero BNOX array
      BNOX = 0d0

      END SUBROUTINE INIT_GLOBAL_NOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_nox
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_NOx deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_NOX
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC      
      IF ( ALLOCATED( BNOX ) ) DEALLOCATE( BNOX ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_NOX
!EOC
      END MODULE GLOBAL_NOX_MOD
