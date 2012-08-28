!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdsoil
!
! !DESCRIPTION: Subroutine RDSOIL reads in soiltype data, fertilizer data, 
!  and monthly soil precipitation data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDSOIL
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_MONTH

!maasa
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      !USE FILE_MOD,    ONLY : IU_SOIL
!maasa end
      
      USE CMN_SIZE_MOD   ! Size parameters
      USE COMMSOIL_MOD ! Soil variables

      USE inquireMod,ONLY : findFreeLUN

      IMPLICIT NONE
#     include "define.h"

!
! !REMARKS:
!  RDSOIL is one of the original GEOS-CHEM subroutines, and has its origins
!  from the GISS-II model that was used at Harvard in the early 90's.  This
!  was cleaned up and improved error checking was added. (bmy, 4/2/02)
!                                                                             .
!  Variables from "commsoil.h" header file:
!  ============================================================================
!  (1 ) NCONSOIL  (INTEGER) : Olson -> soil type mapping index                
!  (2 ) INDEXSOIL (INTEGER) : Array containing grid box indices (I,J)
!  (3 ) SOILFERT  (REAL*8 ) : Array containing fertilizer NOx [ng N/m2/s]  
!  (4 ) SOILPREP  (REAL*8 ) : Array containing 2 months of observed
!                              soil precipitation [mm/day]    
!                                                                             .
!  Files read in by "rdsoil.f":
!  ============================================================================
!  (1 ) DATA_DIR/soil_NOx_200203/soiltype.dat       : Olson and soil land types
!  (2 ) DATA_DIR/soil_NOx_200203/fert_scale.dat     : NOx from fertilizers
!  (3 ) DATA_DIR/soil_NOx_200203/climatprep4x5.dat  : 1x1   monthly soil precip
!                                climatprep2x25.dat : 2x2.5 monthly soil precip
!                                climatprep1x1.dat  : 4x5   monthly soil precip
! 
! !REVISION HISTORY: 
!  05 Jan 1994 - Y. H. Wang, G. M. Gardner, - Initial version
!  (1 ) Be sure to force double precision with the DBLE function and the "D" 
!        exponent, wherever necessary (bmy, 10/6/99)            *
!  (2 ) Now read soil data files directly from the from 
!        DATA_DIR/soil_NOx_200203/ subdirectory.  Now use IOERROR to trap
!        I/O errors across all platforms.  Added comment header.  Updated
!        comments, cosmetic changes. (bmy, 4/2/02)
!  (3 ) Removed obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE as the file unit number,
!        assign it to IUNIT. (bmy, 6/27/02)
!  (4 ) Now reference GEOS_CHEM_STOP from "error_mod.f".  Bug fix: remove 
!        duplicate declaration of IOS.  This causes compile errors for the 
!        ALPHA platform. (gcc, bmy, 11/6/02)
!  (5 ) Now use function GET_MONTH from "time_mod.f".  Now make MONTH a local
!        variable. (bmy, 2/11/03)
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST   = .TRUE.
      INTEGER, SAVE      :: MONSAVE = 0 
      INTEGER            :: I, IUNIT, IOS, J, K, KK, M, M1, MONTH
      INTEGER            :: IU_SOIL
      REAL*8             :: TMP(12)

      ! Name of file to read in
      CHARACTER(LEN=255) :: FILENAME

!Maasa
      INTEGER             :: NI,     NJ, NL, N
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      REAL*4              :: TMPLND(IIPAR,JJPAR)
!Maasa end

      !=================================================================
      ! RDSOIL begins here!
      !=================================================================

      ! Define the file unit
      IUNIT = IU_FILE
      
      ! Get the current month
      MONTH = GET_MONTH()

      ! First-time only initialization
      IF ( FIRST ) THEN
         
         ! Reset First-time flag
         FIRST = .FALSE.

 
         
         !==============================================================
         !Initialize deposition arrays
         !==============================================================          
$
          DRY_NO2  = 0.d0
          DRY_PAN  = 0.d0
          DRY_HNO3 = 0.d0
          DRY_NH3  = 0.d0
          DRY_NH4  = 0.d0
          DRY_NIT  = 0.d0
          WET_HNO3 = 0.d0
          WET_NH3  = 0.d0
          WET_NH4  = 0.d0
          WET_NIT  = 0.d0

          TOTAL_WETOLD = 0.0 ! Testing
          TOTAL_DRYOLD = 0.0 ! Testing

        !==============================================================
         ! Read in soil type data (first pass only)
         !==============================================================

         ! Define soiltype file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/soiltype.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - RDSOIL: Reading ', a )

         ! Find a free file LUN
         IUNIT = findFreeLUN()

         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:1' )

         ! Read header line
         READ( IUNIT, '(a)', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:2' )

         ! Read data
         DO K = 1, NVEGTYPE
	    READ( IUNIT, *, IOSTAT=IOS ) KK, NCONSOIL(KK+1) 
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:3' )
         ENDDO

         ! Close file
         CLOSE( IUNIT )

      ENDIF  

      END SUBROUTINE RDSOIL






