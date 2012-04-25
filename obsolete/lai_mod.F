!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lai_mod
!
! !DESCRIPTION: Module LAI\_MOD reads and stores either AVHRR or MODIS LAI 
!  for calculating MEGAN biogenic VOC emissions.  These LAI values can also
!  be used to replace the data read from the obsolete lai*.global files 
!  (which is used for backwards compatibility w/ legacy code).
!\\
!\\
! !INTERFACE: 
!
      MODULE LAI_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
#     include "define.h"
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: READISOLAI_AVHRR
      PUBLIC :: READISOLAI_MODIS
      PUBLIC :: READISOLAI_05x0666
      PUBLIC :: RDISOLAI
      PUBLIC :: INIT_LAI
      PUBLIC :: CLEANUP_LAI
!
! !PUBLIC DATA MEMBERS:
!
      ! Scalars
      INTEGER, PUBLIC              :: DAYS_BTW_M     ! Days btw the current
                                                     !  and previous months 

      ! Arrays
      REAL*8,  PUBLIC, ALLOCATABLE :: ISOLAI(:,:)    ! LAI for current day
      REAL*8,  PUBLIC, ALLOCATABLE :: MISOLAI(:,:)   ! LAI for current month
      REAL*8,  PUBLIC, ALLOCATABLE :: NMISOLAI(:,:)  ! LAI for next month
      REAL*8,  PUBLIC, ALLOCATABLE :: PMISOLAI(:,:)  ! LAI for previous month
!
! !REVISION HISTORY:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!  (2 ) MEGAN is currently locked to use AVHRR LAI data.  
!        The LAVHRRLAI logical switch controls whether the AVHRR LAI data 
!        is used for the GEIA inventory and dry deposition.
!  (3 ) Modifications for 0.5 x 0.667 nested grid.  Added routine 
!        READISOLAI_05x0666 to read finer-resolution data for GEOS-5 nested
!        grids. (yxw, dan, bmy, 11/6/08)
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!  29 Mar 2012 - R. Yantosca - Updated out-of-date comments
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readisolai_avhrr
!
! !DESCRIPTION: Subroutine READISOLAI_AVHRR reads AVHRR LAI data from bpch 
!  file for the current month, the previous month, and the next month.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READISOLAI_AVHRR( MM )
!
! !USES:
!
      USE CMN_SIZE_MOD                           ! Size parameters
      USE BPCH2_MOD,      ONLY : GET_TAU0        ! Date to TAU0 conversion
      USE BPCH2_MOD,      ONLY : READ_BPCH2      ! Read from bpch format
      USE BPCH2_MOD,      ONLY : GET_RES_EXT     ! Get resolution string
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1    ! Data directory path
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A   ! Regrdidding routine
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D     ! Cast to REAL*8
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: MM                  ! Current month (1-12)
! 
! !REVISION HISTORY: 
!  18 Oct 2005 - D. Abbot, T-M Fu - Initial version
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!  (2 ) Renamed to READISOLAI_AVHRR to differentiate from the newer routine for
!       reading in the MODIS LAI (mpb,2009).
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, K, INDEX, MMM, PMM, IJLOOP
      REAL*4             :: ARRAY(I1x1,J1x1,1)
      REAL*8             :: TAU0
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255)     :: LLFILENAME
      REAL*8, ALLOCATABLE    :: OUTGRID(:,:)
      REAL*8, ALLOCATABLE    :: INGRID(:,:)

      !=================================================================
      ! READISOLAI begins here!
      !=================================================================     

      ! Zero arrays
      MISOLAI  = 0.d0
      NMISOLAI = 0.d0
      ARRAY    = 0.d0

      !------------------------------------
      ! Read current month's lai at (I,J) 
      !------------------------------------
      
      ! Filename
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'leaf_area_index_200412/avhrrlai.global.geos.1x1.2000'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - READISOLAI: Reading ', a )

      ! Get TAU0 value
      TAU0 = GET_TAU0( MM, 1, 2000 )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_200312/latlon_geos1x1.txt'

      ! Allocate arrays for regridding
      ALLOCATE(INGRID(I1x1,J1x1))
      ALLOCATE(OUTGRID(IIPAR,JJPAR))

      ! Regrid from 1x1 to current grid resolution
      INGRID  = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      MISOLAI = OUTGRID

      !------------------------------------
      ! Read next month's lai at (I,J) 
      !------------------------------------

      ! MMM is next month
      MMM = MM + 1
      IF ( MMM == 13 ) MMM = 1

      ! TAU0 for 1st day of next month
      TAU0 = GET_TAU0( MMM, 1, 2000 ) 
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      INGRID = ARRAY(:,:,1)        
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      NMISOLAI = OUTGRID

      !------------------------------------
      ! Read previous month's lai at (I,J) 
      !------------------------------------

      ! PMM is previous month
      PMM = MM - 1
      IF ( PMM == 0 ) PMM = 12

      ! TAU0 for 1st day of previous month
      TAU0 = GET_TAU0( PMM, 1, 2000 ) 
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      INGRID = ARRAY(:,:,1)   
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      PMISOLAI = OUTGRID

      ! Deallocate arrays for regridding
      DEALLOCATE(INGRID)
      DEALLOCATE(OUTGRID)

      END SUBROUTINE READISOLAI_AVHRR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readisolai_modis
!
! !DESCRIPTION: Subroutine READISOLAI\_MODIS reads MODIS LAI data from bpch 
!  file for the current  month, the previous month, and the next month.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READISOLAI_MODIS( MM, YYYY )
!
! !USES:
!
      USE CMN_SIZE_MOD                           ! Size parameters
      USE BPCH2_MOD,      ONLY : GET_TAU0        ! Date to TAU0 conversion
      USE BPCH2_MOD,      ONLY : READ_BPCH2      ! Read from bpch format
      USE BPCH2_MOD,      ONLY : GET_RES_EXT     ! Get resolution string
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1    ! Data directory path
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A   ! Regrdidding routine
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D     ! Cast to REAL*8
!                                                
! !INPUT PARAMETERS:                             
!                                                
      INTEGER, INTENT(IN)    :: MM               ! Current month (1-12)
      INTEGER, INTENT(IN)    :: YYYY             ! Current year
! 
! !REVISION HISTORY: 
!  03 Jan 2009 - M. Barkley  - Initial version
!  (1 ) I'm sure there is a much easier way of coding this (mpb,2009).
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, K, INDEX, MMM, PMM, IJLOOP
      REAL*4                 :: ARRAY(I1x1,J1x1,1)
      REAL*8                 :: TAU0
      CHARACTER(LEN=255)     :: FILENAME 
      CHARACTER(LEN=255)     :: LLFILENAME
      REAL*8, ALLOCATABLE    :: OUTGRID(:,:)
      REAL*8, ALLOCATABLE    :: INGRID(:,:)

      ! For MODIS add past, current and next years (mpb,2009)
      CHARACTER(LEN=4)       :: YEAR 
      INTEGER                :: IYYYY !Current Year
      INTEGER                :: NYYYY !Next Year
      INTEGER                :: PYYYY !Past Year

      CHARACTER(LEN=255)     :: MODIS_CLIM

      !=================================================================
      ! READISOLAI_MODIS begins here!
      !=================================================================     

      ! Zero arrays
      MISOLAI  = 0.d0
      NMISOLAI = 0.d0
      ARRAY    = 0.d0

      ! Short name for MODIS climatology
      !MODIS_CLIM =  TRIM( DATA_DIR_1x1 ) // 
      MODIS_CLIM = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' //
     &             '1985/MODIS.LAIv.V5.geos.1x1.1985'

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_200312/latlon_geos1x1.txt'

      IYYYY = 0
      PYYYY = 0
      NYYYY = 0

      ! ++++++++++++++++++++++++++++++++++
      ! Get LAI data for CURRENT month
      ! ++++++++++++++++++++++++++++++++++

      ! - Set current year -
      IYYYY = YYYY

      ! If valid (2000 to 2008) use data for that year 
      IF ( IYYYY >= 2000 .AND. IYYYY <= 2008 ) THEN      

         ! Filename 
         WRITE( YEAR , '(I4)' ) IYYYY   

         !FILENAME = TRIM( DATA_DIR_1x1 ) //  
         FILENAME = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' // 
     &              YEAR // '/MODIS.LAIv.V5.geos.1x1.' // YEAR

        PRINT*, '    - current Month = ' , MM , 'Year = ' , IYYYY

         ! Get TAU0 value
         TAU0 = GET_TAU0( MM, 1, IYYYY )

       ELSE
          
         ! Filename
         FILENAME = TRIM( MODIS_CLIM )

          PRINT*, '    - current : Month = ' , MM , 'Year = ' , 1985

         ! Get TAU0 value
         TAU0 = GET_TAU0( MM, 1, 1985 )

      END IF 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'MODIS', 2,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Allocate arrays for regridding
      ALLOCATE(INGRID(I1x1,J1x1))
      ALLOCATE(OUTGRID(IIPAR,JJPAR))

      ! Regrid from 1x1 to current grid resolution
      INGRID = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      MISOLAI = OUTGRID

      ! ++++++++++++++++++++++++++++++++++
      ! Get LAI data for NEXT month
      ! ++++++++++++++++++++++++++++++++++

      ! MMM is next month
      MMM = MM + 1

      ! Is this the start of a new year... 
      IF ( MMM == 13 ) THEN

         MMM = 1

         ! Increment the year by 1
         NYYYY = IYYYY + 1 

         ! ...does this new year fall into the valid MODIS range... 
         IF ( NYYYY >= 2000 .AND. NYYYY <= 2008 ) THEN 

            ! Filename 
            WRITE( YEAR , '(I4)' ) NYYYY   

            !FILENAME = TRIM( DATA_DIR_1x1 ) // 
            FILENAME = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' // 
     &                 YEAR // '/MODIS.LAIv.V5.geos.1x1.'// YEAR

            PRINT*, '    - next : Month = ' , MMM , 'Year = ' , NYYYY

            ! Get TAU0 value
            TAU0 = GET_TAU0( MMM, 1, NYYYY )

         ELSE ! ...if not use climatology...
  
            ! Filename
            FILENAME = TRIM( MODIS_CLIM )
            PRINT*, '    - next : Month = ' , MMM , 'Year = ' , 1985

            ! Get TAU0 value
            TAU0 = GET_TAU0( MMM, 1, 1985 )

         END IF 

      ELSE ! If the year is the same when MM --> MM + 1  

         ! If valid (2000 to 2008) use data for that year 
         IF ( IYYYY >= 2000 .AND. IYYYY <= 2008 ) THEN      
 
            ! Filename 
            WRITE( YEAR , '(I4)' ) IYYYY   

            !FILENAME = TRIM( DATA_DIR_1x1 ) // 
            FILENAME = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' // 
     &                 YEAR // '/MODIS.LAIv.V5.geos.1x1.' //YEAR

           PRINT*, '    - next : Month = ' , MMM , 'Year = ' , IYYYY

            ! Get TAU0 value
            TAU0 = GET_TAU0( MMM, 1, IYYYY )
         
         ELSE ! ...if not use climatology...
          
            ! Filename
            FILENAME = TRIM( MODIS_CLIM )
            PRINT*, '    - next : Month = ' , MMM , 'Year = ' , 1985

            ! Get TAU0 value
            TAU0 = GET_TAU0( MMM, 1, 1985 )

         END IF 
      END IF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'MODIS', 2,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      INGRID = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      NMISOLAI = OUTGRID

      ! ++++++++++++++++++++++++++++++++++
      ! Get LAI data for PAST month
      ! ++++++++++++++++++++++++++++++++++

      ! PMM is previous month
      PMM = MM - 1

      ! If we dip back in to the previous year 
      IF ( PMM == 0 ) THEN

         PMM   = 12 
         PYYYY = IYYYY - 1

         ! Does this past year fall into the valid MODIS range...
         IF ( PYYYY >= 2000 .AND. PYYYY <= 2008 ) THEN 

            ! Filename 
            WRITE( YEAR , '(I4)' ) PYYYY   

            !FILENAME = TRIM( DATA_DIR_1x1 ) // 
            FILENAME = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' //
     &                 YEAR // '/MODIS.LAIv.V5.geos.1x1.'// YEAR

            PRINT*, '    - past : Month = ' , PMM , 'Year = ' , PYYYY

            ! Get TAU0 value
            TAU0 = GET_TAU0( PMM, 1, PYYYY )

         ELSE ! ...if not use climatology...

            ! Filename
            FILENAME = TRIM( MODIS_CLIM )
            PRINT*, '    - past : Month = ' , PMM , 'Year = ' , 1985

            ! Get TAU0 value
            TAU0 = GET_TAU0( PMM, 1, 1985 )

         END IF 

      ELSE ! If the year is the same when MM -->  MM - 1  

         ! If valid (2000 to 2006) use data for that year 
         IF ( IYYYY >= 2000 .AND. IYYYY <= 2008 ) THEN      
 
            ! Filename 
            WRITE( YEAR , '(I4)' ) IYYYY   

            !FILENAME = TRIM( DATA_DIR_1x1 ) // 
            FILENAME = TRIM( DATA_DIR_1x1 ) // 'MODIS_LAI_200911/' // 
     &                 YEAR // '/MODIS.LAIv.V5.geos.1x1.' // YEAR

            PRINT*, '   -  past : Month = ' , PMM , 'Year = ' , IYYYY

            ! Get TAU0 value
            TAU0 = GET_TAU0( PMM, 1, IYYYY )
         
         ELSE ! ...if not use climatology...
          
            ! Filename
            FILENAME = TRIM( MODIS_CLIM )
            PRINT*, '    - past : Month = ' , PMM , 'Year = ' , 1985

            ! Get TAU0 value
            TAU0 = GET_TAU0( PMM, 1, 1985 )

         END IF 

      END IF 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'MODIS', 2,  
     &                 TAU0,      I1x1,   J1x1,    
     &                 1,         ARRAY,  QUIET=.TRUE. )

      ! Regrid from 1x1 to current grid resolution
      INGRID = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, 360, 182, INGRID, OUTGRID, 0 )
      PMISOLAI = OUTGRID

      ! Deallocate arrays for regridding
      DEALLOCATE(INGRID)
      DEALLOCATE(OUTGRID)

100   FORMAT( '     - READISOLAI: Reading ', a )
     
      END SUBROUTINE READISOLAI_MODIS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdisolai_05x0666
!
! !DESCRIPTION: Subroutine READISOLAI reads AVHRR LAI data from bpch file 
!  for the current month, the previous month, and the next month.  
!  Specially constructed to read hi-res data for the GEOS-5 0.5 x 0.666 
!  nested grid simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READISOLAI_05x0666( MM )
!
! !USES:
!
      USE CMN_SIZE_MOD                           ! Size parameters
      USE BPCH2_MOD,      ONLY : GET_TAU0        ! Date to TAU0 conversion
      USE BPCH2_MOD,      ONLY : READ_BPCH2      ! Read from bpch format
      USE BPCH2_MOD,      ONLY : GET_RES_EXT     ! Get resolution string
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1    ! Data directory path
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D     ! Cast to REAL*8
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: MM                  ! Current month (1-12)
!
! !REVISION HISTORY: 
!  06 Nov 2008 - Y. Wang, D. Chen - Initial version
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, K, INDEX, MMM, PMM, IJLOOP
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*8             :: TAU0
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READISOLAI begins here!
      !=================================================================     

      ! Zero arrays
      MISOLAI  = 0.d0
      NMISOLAI = 0.d0
      ARRAY    = 0.d0

      !------------------------------------
      ! Read current month's lai at (I,J) 
      !------------------------------------

      ! Filename
      FILENAME = TRIM( DATA_DIR ) //
     &      'leaf_area_index_200412/avhrrlai.global.geos.05x0666.2000'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READISOLAI: Reading ', a )

      ! Get TAU0 value
      TAU0 = GET_TAU0( MM, 1, 2000 )

      ! Read 1x1 LAI data [cm2/cm2]
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IIPAR,   JJPAR,
     &                 1,         ARRAY,  QUIET=.TRUE. )


      CALL TRANSFER_2D( ARRAY(:,:,1), MISOLAI )


      !------------------------------------
      ! Read next month's lai at (I,J) 
      !------------------------------------

      ! MMM is next month
      MMM = MM + 1
      IF ( MMM == 13 ) MMM = 1

      ! TAU0 for 1st day of next month
      TAU0 = GET_TAU0( MMM, 1, 2000 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IIPAR,   JJPAR,
     &                 1,         ARRAY,  QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), NMISOLAI )


      !------------------------------------
      ! Read previous month's lai at (I,J) 
      !------------------------------------

      ! PMM is previous month
      PMM = MM - 1
      IF ( PMM == 0 ) PMM = 12

      ! TAU0 for 1st day of previous month
      TAU0 = GET_TAU0( PMM, 1, 2000 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'AVHRR', 1,
     &                 TAU0,      IIPAR,   JJPAR,
     &                 1,         ARRAY,  QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), PMISOLAI )

      END SUBROUTINE READISOLAI_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdisolai
!
! !DESCRIPTION: Subroutine RDISOLAI sets ISOLAI daily.  The stored monthly 
!  LAI are used for the middle day in the month and LAIs are interpolated 
!  for other days.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RDISOLAI( JDAY, MONTH, YEAR )
!
! !USES:
!
      USE CMN_SIZE_MOD                         ! Size parameters
      USE TIME_MOD,    ONLY : ITS_A_LEAPYEAR   ! Leap year determination
      USE LOGICAL_MOD, ONLY : LMODISLAI        ! Use MODIS or AVHRR LAI
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: JDAY              ! Day of year (0-365 or 366)
      INTEGER, INTENT(IN) :: MONTH             ! Current month (1-12)
      INTEGER, INTENT(IN) :: YEAR              ! Current year
!
! !REVISION HISTORY: 
!  20 Oct 2005 - D. Abbot, T-M Fu - Initial version
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (7/8/03).  Updated  
!        and modified for the standard code by May Fu (11/2004).
!  (2 ) Now call READISOLAI_05x0666 to read hi-res LAI data if we are doing a 
!        GEOS-5 0.5 x 0.666 nested grid simulation. (yxw, dan, bmy, 11/6/08)
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, IMUL, ITD, IJLOOP, MM
      INTEGER             :: YYYY              ! Include year (mpb,2009)
      INTEGER, SAVE       :: LAST_MM = -1 
      REAL*8              :: FRACTION

      ! specify midmonth day for year 2000
      INTEGER, PARAMETER  :: STARTDAY(13) = 
     &                         (/  15,  45,  74, 105, 135, 166,
     &                            196, 227, 258, 288, 319, 349, 380/)

      !=================================================================
      ! RDISOLAI begins here!
      !=================================================================

      ! Find the month if we index by midmonth
      ! & now include the year (mpb,2009)
      CALL FINDMON( JDAY, MONTH, YEAR, MM, YYYY, STARTDAY )

      ! Read new data if it's a new LAI month
      IF ( MM /= LAST_MM ) THEN

#if   defined( GRID05x0666 )
         CALL READISOLAI_05x0666( MM )   ! GEOS-5 nested grid simulation
#else 

         ! Global simulations
         IF ( LMODISLAI ) THEN
            CALL READISOLAI_MODIS( MM , YYYY )
         ELSE 
            CALL READISOLAI_AVHRR( MM )
         END IF 
         ! CALL READISOLAI( MM )           ! Global simulations
#endif

         ! Save for next month
         LAST_MM = MM
      ENDIF
      
      ! IMUL is days since midmonth
      ! ITD  is days between midmonths
      IF ( JDAY < STARTDAY(1) ) THEN
         IMUL = 365 + JDAY - STARTDAY(12) 
         ITD  = 31
      ELSE
         IMUL = JDAY           - STARTDAY(MM)
         ITD  = STARTDAY(MM+1) - STARTDAY(MM)
      END IF

      ! Archive the days between midmonths in the LAI data
      DAYS_BTW_M     = ITD

      ! Fraction of the LAI month that we are in
      FRACTION       = DBLE( IMUL ) / DBLE( ITD ) 
       
      ! Interpolate to daily LAI value
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ISOLAI(I,J) = MISOLAI(I,J) + 
     &                ( FRACTION * ( NMISOLAI(I,J) - MISOLAI(I,J) ) )
      ENDDO
      ENDDO

      END SUBROUTINE RDISOLAI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lai
!
! !DESCRIPTION: Subroutine INIT_LAI allocates and initializes arrays for 
!  either AVHRR or MODIS LAI.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_LAI
!
! !USES:
!
      USE CMN_SIZE_MOD                  ! Size parameters
      USE ERROR_MOD, ONLY : ALLOC_ERR   ! Allocation error trapping
! 
! !REVISION HISTORY: 
!  08 Jul 2003 - D. Abbot, T-M. Fu - Initial version
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!  29 Mar 2012 - R. Yantosca - Exit if we have already called INIT_LAI before
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: AS
  
      !=================================================================
      ! INIT_LAI begins here!
      !=================================================================

      ! Skip if we have already called this routine before
      IF ( .not. FIRST ) THEN
         FIRST = .FALSE.
         RETURN
      ENDIF

      ALLOCATE( ISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ISOLAI' )
      ISOLAI = 0d0

      ALLOCATE( MISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MISOLAI' )
      MISOLAI = 0d0

      ALLOCATE( NMISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NMISOLAI' )
      NMISOLAI = 0d0

      ALLOCATE( PMISOLAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PMISOLAI' )
      PMISOLAI = 0d0        

      END SUBROUTINE INIT_LAI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_lai
!
! !DESCRIPTION: Subroutine CLEANUP\_ISOLAI deallocates all allocated arrays 
!  at the end of a GEOS-CHEM model run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_LAI
! 
! !REVISION HISTORY: 
!  08 Jul 2003 - D. Abbot, T-M. Fu - Initial version
!  29 Mar 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! CLEANUP_ISOLAI begins here!
      !=================================================================
      IF ( ALLOCATED( ISOLAI   ) ) DEALLOCATE( ISOLAI   )
      IF ( ALLOCATED( MISOLAI  ) ) DEALLOCATE( MISOLAI  )
      IF ( ALLOCATED( NMISOLAI ) ) DEALLOCATE( NMISOLAI )
      IF ( ALLOCATED( PMISOLAI ) ) DEALLOCATE( PMISOLAI )

      END SUBROUTINE CLEANUP_LAI
!EOC
      END MODULE LAI_MOD
