!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: readlai
!
! !DESCRIPTION: Subroutine READLAI reads the leaf area indices from disk 
!  for two months. (yhw, gmg, djj, 1994; bmy, 12/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READLAI( MM, YYYY )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE
      USE LOGICAL_MOD,   ONLY : LAVHRRLAI
      USE LOGICAL_MOD,   ONLY : LMODISLAI ! (mpb,2009)

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_VEL_MOD    ! XLAI, XLAI2
      USE CMN_DEP_MOD    ! IREG, ILAND, IUSE

      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
#     include "define.h"

!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MM
      INTEGER, INTENT(IN) :: YYYY ! (mpb,2009)
!
! !REVISION HISTORY:
!  06 Oct 1999 - R. Yantosca - Be sure to force double precision with the 
!                              DBLE function and the "D" exponent, wherever 
!                              necessary             
!  05 Jul 2001 - R. Yantosca - Now reads the LAI files directly from the data 
!                              directory, so you don't have to create symbolic 
!                              links anymore           
!  27 Feb 2002 - R. Yantosca - Deleted obsolete code                        
!  25 Jun 2002 - R. Yantosca - Replaced IMX with IIPAR and JMX with JJPAR
!  31 Jul 2002 - R. Yantosca - Now reference IU_FILE from "file_mod.f"   
!  13 Nov 2002 - R. Yantosca - Now define FILENAME and echo FILENAME to stdout.
!                              Now use F90 style declaration statements.
!                              Cleaned up old code.
!  20 Jul 2004 - R. Yantosca - Now references DATA_DIR from "directory_mod.f"
!  20 Dec 2004 - M. Fu       - Now use AVHRR LAI derived leaf-area index data 
!                              (stored in the leaf_area_index_200412 subdir of 
!                              DATA_DIR) if the logical switch LAVHRRLAI=T. 
!                              Otherwise use the old LAI data.
!EOP
!------------------------------------------------------------------------------
!BOC
!     
! ! LOCAL VARIABLES
!
      INTEGER             :: I, INDEX, J, K, MMM
      CHARACTER(LEN=2)    :: CMONTH(12) = (/ '01','02','03','04',
     &                                       '05','06','07','08',
     &                                       '09','10','11','12'/)
      CHARACTER(LEN=255)  :: FILENAME

      ! New for MODIS & AVHRR (mpb,2009)
      CHARACTER(LEN=4)    :: YEAR 
      INTEGER             :: IYYYY, NYYYY 

      !=================================================================
      ! READLAI begins here!
      !=================================================================

      ! Zero XLAI, XLAI2
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO K = 1, IREG(I,J)
         XLAI(I,J,K)  = 0.D0
         XLAI2(I,J,K) = 0.D0
      ENDDO
      ENDDO
      ENDDO

      ! Initialize (mpb,2009)
      IYYYY = 0
      NYYYY = 0

      !=================================================================
      ! Read current month's lai (XLAI) at (I,J) and for landtype K
      !=================================================================

      ! Pick proper filename for the old Yuhang Wang LAI, or
      ! for AVHRR satellite-derived LAI (tmf, bmy, 12/20/04)
      IF ( LAVHRRLAI ) THEN
         FILENAME = TRIM( DATA_DIR ) // 'leaf_area_index_200412/lai' //
     &        CMONTH(MM)       // '.global'
      ELSE IF ( LMODISLAI ) THEN ! Now include MODIS LAI (mpb,2009)
         
         IYYYY = YYYY
         
         IF ( IYYYY >= 2000 .AND. IYYYY <= 2008 ) THEN 
            
            ! Filename 
            WRITE( YEAR , '(I4)' ) IYYYY
            
            FILENAME = TRIM( DATA_DIR ) // 'MODIS_LAI_200911/' // 
     &                 YEAR // '/lai' // CMONTH(MM) // '.global'
         ELSE 
            FILENAME = TRIM( DATA_DIR ) // 'MODIS_LAI_200911/1985/' //
     &                 'lai' // CMONTH(MM) // '.global'
         END IF
         
      ELSE
         FILENAME = TRIM( DATA_DIR ) // 'leaf_area_index_200202/lai' //
     &              CMONTH(MM)       // '.global'
      ENDIF

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READLAI: Reading ', a )
      
      ! Open file
      IU_FILE = findFreeLUN()
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD' )

      ! Read until EOF
 10   READ( IU_FILE, '(3i3,20f5.1)', END=20 ) 
     &     I, J, INDEX, ( XLAI(I,J,K), K=1,INDEX )
      GOTO 10

      ! Close file
 20   CLOSE( IU_FILE )
      
      ! Save for next month
      MMM = MM
      IF(MMM .EQ. 12) THEN
         MMM = 0
         IYYYY = YYYY + 1 ! Increment year by 1 (mpb,2009)
      ENDIF

      !=================================================================
      ! Read following month's lai (XLAI2) at (I,J) and for landtype K 
      !=================================================================

      ! Pick proper filename for the old Yuhang Wang LAI, or
      ! for AVHRR satellite-derived LAI (tmf, bmy, 12/20/04)
      IF ( LAVHRRLAI ) THEN
         FILENAME = TRIM( DATA_DIR ) // 'leaf_area_index_200412/lai' //
     &              CMONTH(MMM+1)    // '.global'

      ELSE IF ( LMODISLAI ) THEN ! Now include MODIS LAI (mpb,2009)

         IF ( IYYYY >= 2000 .AND. IYYYY <= 2008 ) THEN

            WRITE( YEAR , '(I4)' ) IYYYY

            !FILENAME = TRIM( DATA_DIR ) // 'MODIS_LAIv_v5/' 
            FILENAME = TRIM( DATA_DIR ) // 'MODIS_LAI_200911/' // 
     &                 YEAR // '/lai' // CMONTH(MMM+1) // '.global'
         ELSE 
            !FILENAME = TRIM( DATA_DIR ) // 
            FILENAME = TRIM( DATA_DIR ) // 'MODIS_LAI_200911/1985/' // 
     &                 'lai' // CMONTH(MMM+1) // '.global'
         END IF

      ELSE
         FILENAME = TRIM( DATA_DIR ) // 'leaf_area_index_200202/lai' //
     &              CMONTH(MMM+1)    // '.global'
      ENDIF

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Open file
      IU_FILE = findFreeLUN()
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD' )

      ! Read until EOF
 30   READ( IU_FILE, '(3i3,20f5.1)', END=40 )
     &     I, J, INDEX, ( XLAI2(I,J,K), K=1,INDEX )
      GOTO 30

      ! Close file
 40   CLOSE( IU_FILE )
      
      END SUBROUTINE READLAI
!EOC
