! $Id: rdsoil.f,v 1.2 2004/09/21 18:04:17 bmy Exp $
      SUBROUTINE RDSOIL
!
!******************************************************************************
!  Subroutine RDSOIL reads in soiltype data, fertilizer data, and monthly
!  soil precipitation data. (yhw, gmg, djj, bmy, 1994, 7/20/04)
!
!  RDSOIL is one of the original GEOS-CHEM subroutines, and has its origins
!  from the GISS-II model that was used at Harvard in the early 90's.  This
!  was cleaned up and improved error checking was added. (bmy, 4/2/02)
!
!  Variables from "commsoil.h" header file:
!  ============================================================================
!  (1 ) NCONSOIL  (INTEGER) : Olson -> soil type mapping index                
!  (2 ) INDEXSOIL (INTEGER) : Array containing grid box indices (I,J)
!  (3 ) SOILFERT  (REAL*8 ) : Array containing fertilizer NOx [ng N/m2/s]  
!  (4 ) SOILPREP  (REAL*8 ) : Array containing 2 months of observed
!                              soil precipitation [mm/day]    
!
!  Files read in by "rdsoil.f":
!  ============================================================================
!  (1 ) DATA_DIR/soil_NOx_200203/soiltype.dat       : Olson and soil land types
!  (2 ) DATA_DIR/soil_NOx_200203/fert_scale.dat     : NOx from fertilizers
!  (3 ) DATA_DIR/soil_NOx_200203/climatprep4x5.dat  : 1x1   monthly soil precip
!                                climatprep2x25.dat : 2x2.5 monthly soil precip
!                                climatprep1x1.dat  : 4x5   monthly soil precip
! 
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,      ONLY : GET_MONTH
      
      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
!------------------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN"        ! MONTH
!#     include "CMN_SETUP"  ! for DATA_DIR (bmy, 7/6/01)
!------------------------------------------------------------
#     include "commsoil.h" ! Soil variables

      ! Local variables
      LOGICAL, SAVE      :: FIRST   = .TRUE.
      INTEGER, SAVE      :: MONSAVE = 0 
      INTEGER            :: I, IUNIT, IOS, J, K, KK, M, M1, MONTH
      REAL*8             :: TMP(12)

      ! Name of file to read in
      CHARACTER(LEN=255) :: FILENAME

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
         ! Read in soil type data (first pass only)
         !==============================================================

         ! Define soiltype file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/soiltype.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - RDSOIL: Reading ', a )

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

         !==============================================================
         ! Read in fertilizer data (first pass only)
         ! Units are [ng N/m2/s] 
         !==============================================================

         ! Define fertilizer file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/fert_scale.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:4' )

         ! Read data -- save (I,J) pairs into INDEXSOIL array
         DO M = 1, NLAND
	    READ( IUNIT,*, IOSTAT=IOS ) 
     &           ( INDEXSOIL(I,M), I=1,2 ), SOILFERT(M)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:5' )
         ENDDO

         ! Close file
         CLOSE( IUNIT )

      ENDIF
         
      !=================================================================
      ! Read in monthly soil precipitation data
      !=================================================================

      ! Only read data when we have entered a new month...
      IF ( MONSAVE /= MONTH ) THEN
      
         ! Save the current month 
         MONSAVE = MONTH
      
         ! M1 is the previous month
         IF ( MONTH == 1 ) THEN
            M1 = 12
         ELSE
            M1 = MONTH - 1
         END IF

         ! Define soil precip file name
         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/climatprep' //
     &              GET_RES_EXT()    // '.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Open soil precip file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:6' )

         ! Loop over Olson land types
         DO M = 1, NLAND

            ! Read monthly soil precip data for each (I,J) box
            READ( IUNIT, *, IOSTAT=IOS ) I, J, ( TMP(K), K=1,12 )
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:7' )

            ! Error check -- make sure that each (I,J) pair has both
            ! soil precip data and fertilizer data defined
            IF ( INDEXSOIL(1,M) /= I .OR. INDEXSOIL(2,M) /= J ) THEN
               WRITE(6,*) 'CORRUPTED TEMPCORR OR CLMATPRECIP DATA'
               WRITE(6,*) 'CHECK (I,J)',I,J
               CALL GEOS_CHEM_STOP
            ELSE
               SOILPREP(1,M) = TMP(M1)
               SOILPREP(2,M) = TMP(MONTH)
            ENDIF
         ENDDO

         ! Close file
         CLOSE( IUNIT )

      ENDIF

      ! Return to calling program
      END SUBROUTINE RDSOIL






