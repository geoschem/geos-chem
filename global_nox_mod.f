! $Id: global_nox_mod.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      MODULE GLOBAL_NOX_MOD
!
!******************************************************************************
!  Module GLOBAL_NOX_MOD contains variables and routines for reading the
!  global monthly mean NOX concentration from disk. (bmy, 7/28/00, 4/28/03)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) BNOX (REAL*8)      : stores global monthly mean NOx field [ppbv]
!
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_GLOBAL_NOX     : reads global monthly mean NOx from disk
!  (2 ) INIT_GLOBAL_NOX    : allocates & initializes the NOx array
!  (3 ) CLEANUP_GLOBAL_NOX : deallocates the NOx array
!
!  GEOS-CHEM modules referenced by global_nox_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f  : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f  : Module containing NaN and other error check routines
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean BNOX field
      REAL*8, ALLOCATABLE :: BNOX(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_NOX( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_NOX reads global NOX from binary punch files from a 
!  a full chemistry run.  This NOx data is needed to calculate the CO yield
!  from isoprene oxidation. (bmy, 7/28/00, 4/28/03)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Now use version of GET_TAU0 with 3 arguments.  Now call READ_BPCH2 
!        with IGLOB,JGLOB,LGLOB.  Call TRANSFER_3D to cast from REAL*4 to 
!        REAL*8 and to regrid to 30 levels for GEOS-3 (if necessary).  ARRAY 
!        should now be of size (IGLOB,JGLOB,LGLOB). (bmy, 1/14/02)
!  (2 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas.  Also make sure
!        to define FILENAME before printing it (bmy, 4/28/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_3D

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_SETUP"   ! TEMP_DIR, DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME
      CHARACTER(LEN=255)   :: FIELD_DIR, RGNAME, TEMPO, CHAROP
      CHARACTER(LEN=3)     :: BMONTH(12) = (/ 'jan', 'feb', 'mar', 
     &                                        'apr', 'may', 'jun', 
     &                                        'jul', 'aug', 'sep', 
     &                                        'oct', 'nov', 'dec' /)

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

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
      RGNAME = TRIM( FIELD_DIR )   // 'ctm.bpch.'         // 
     &         BMONTH( THISMONTH ) // '.'                 // 
     &         GET_NAME_EXT()      // TRIM( ZIP_SUFFIX )

      ! Construct the command to unzip the file & copy to TEMP_DIR
      CHAROP = TRIM( 'gzcat' )   // SPACE(1:1)          //
     &         TRIM( RGNAME  )   // TRIM( REDIRECT  )   //
     &         SPACE(1:1)        // TRIM( TEMP_DIR  )   //
     &         TRIM( TEMPO   )

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
      FILENAME = TRIM( TEMP_DIR ) // TRIM( TEMPO )

      ! Echo info
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_NOX: Reading NOX from: ', a )
      
      ! Read NOX data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,     XTAU,  
     &                 IGLOB,    JGLOB,      LGLOB, ARRAY )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_3D( ARRAY, BNOX )

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_NOX

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_NOX
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_NOX allocates and zeroes the NOX array, which 
!  holds global monthly mean NOX concentrations. (bmy, 7/28/00, 10/15/02)
!
!  NOTES:
!  (1 ) BNOX now needs to be sized (IIPAR,JJPAR,LLPAR) (bmy, 1/14/02)
!  (2 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS

      ! Allocate NOX array
      ALLOCATE( BNOX( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BNOX' )

      ! Zero BNOX array
      BNOX = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_NOX
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_NOX
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_NOX deallocates the NOX array.
!******************************************************************************
!                               
      IF ( ALLOCATED( BNOX ) ) DEALLOCATE( BNOX ) 
     
      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_NOX

!------------------------------------------------------------------------------

      END MODULE GLOBAL_NOX_MOD
