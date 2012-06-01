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
      USE FILE_MOD,    ONLY : IU_SOIL
!maasa end
      
      USE CMN_SIZE_MOD   ! Size parameters
      USE COMMSOIL_MOD ! Soil variables

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST   = .TRUE.
      INTEGER, SAVE      :: MONSAVE = 0 
      INTEGER            :: I, IUNIT, IOS, J, K, KK, M, M1, MONTH
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

!||||||||||||||||||||||||||||||||||||||
!Read landtype module, maasa
!||||||||||||||||||||||||||||||||||||||

 !FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_201004/Landtype.bpch' !Maasa
	FILENAME = '/as/home/m/maasakkers/data/soil_NOx_201004/Landtype.bpch'

      WRITE(6,*) FILENAME

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_SOIL, FILENAME )
      
      
      !Loop over 24 Landtypes

      N=1

      DO
      

         READ( IU_SOIL, IOSTAT=IOS )       
     &         MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
      
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_Aws:1')

         READ( IU_SOIL, IOSTAT=IOS ) 
     &         CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &         NI,       NJ,    NL, IFIRST, JFIRST, LFIRST,
     &         NSKIP

         IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_Aws:2')


         ! Put land fraction into temporary array
         READ( IU_SOIL, IOSTAT=IOS ) 
     &       ( (   TMPLND(I,J), I=1,NI ), J=1,NJ  )

         IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_Aws:3')

         !WRITE(6,*) 'I = 100 J =61 ', N, ':', TMPLND(100,61) 


         ! Put into global array
         LAND2(:,:,N) = TMPLND

         N=N+1
  
      ENDDO
         CLOSE( IU_SOIL )

         !WRITE(6,*) 'I = 100 J =61 ', land2(100,61,:)

!||||||||||||||||||||||||||||||||||||||
!End of landtype module, maasa
!||||||||||||||||||||||||||||||||||||||

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
         !FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_201004/fert_res.dat' !Maasatemp
	 FILENAME = '/as/home/m/maasakkers/data/soil_NOx_201004/fert_res.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:4' )

         ! Read data -- save (I,J) pairs into INDEXSOIL array
         DO M = 1, NFERT
	    READ( IUNIT,*, IOSTAT=IOS ) 
     &           ( INDEX_INT(I,M), I=1,2 ), SOIL_INT(M) !maasatemp
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:5' )
         ENDDO

         DO M = 1, NLAND
                INDEXSOIL(1,M) = INDEX_INT(1,(1 + (M-1)*365))
                INDEXSOIL(2,M) = INDEX_INT(2,(1 + (M-1)*365))
		SOILFERT(M)    = SOIL_INT((1 + (M-1)*365))
         ENDDO

         ! Close file
         CLOSE( IUNIT )

         
!      !=================================================================
!      ! Read in climate type ###Rynda
!      !=================================================================
         ! Define arid/semi-arid 1-0 dataset location
  	 !FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_201004/clim.dat' !Maasatemp
		FILENAME = '/as/home/m/maasakkers/data/soil_NOx_201004/clim.dat'

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Open file
         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:4' )

         ! Read data -- save (I,J) pairs into INDEXSOIL array
         DO M = 1, NLAND
	    READ( IUNIT,*, IOSTAT=IOS ) 
     &           ( INDEXSOIL(I,M), I=1,2 ), CLIM(M)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:5' )
         ENDDO

         ! Close file
         CLOSE( IUNIT )

      ENDIF  
!      !=================================================================
!      ! Read in monthly soil precipitation data
!      !=================================================================

!      ! Only read data when we have entered a new month...
!      IF ( MONSAVE /= MONTH ) THEN
      
!         ! Save the current month 
!         MONSAVE = MONTH
      
!         ! M1 is the previous month
!         IF ( MONTH == 1 ) THEN
!            M1 = 12
!         ELSE
!            M1 = MONTH - 1
!         END IF

!        ! Define soil precip file name
!         FILENAME = TRIM( DATA_DIR ) // 'soil_NOx_200203/climatprep' //
!     &              GET_RES_EXT()    // '.dat'

!         ! Echo filename
!         WRITE( 6, 100 ) TRIM( FILENAME )

!         ! Open soil precip file
!         OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
!         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:6' )

!         ! Loop over Olson land types
!         DO M = 1, NLAND
!
!            ! Read monthly soil precip data for each (I,J) box
!            READ( IUNIT, *, IOSTAT=IOS ) I, J, ( TMP(K), K=1,12 )
!            IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'rdsoil:7' )
!
!            ! Error check -- make sure that each (I,J) pair has both
!            ! soil precip data and fertilizer data defined
!            IF ( INDEXSOIL(1,M) /= I .OR. INDEXSOIL(2,M) /= J ) THEN
!               WRITE(6,*) 'CORRUPTED TEMPCORR OR CLMATPRECIP DATA'
!               WRITE(6,*) 'CHECK (I,J)',I,J
!               CALL GEOS_CHEM_STOP
!            ELSE
!               SOILPREP(1,M) = TMP(M1)
!               SOILPREP(2,M) = TMP(MONTH)
!            ENDIF
!         ENDDO
!
!         ! Close file
!         CLOSE( IUNIT )
!
!      ENDIF

      ! Return to calling program
      END SUBROUTINE RDSOIL






