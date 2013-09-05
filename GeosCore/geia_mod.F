! $Id: geia_mod.f,v 1.3 2010/02/02 16:57:53 bmy Exp $
      MODULE GEIA_MOD
!
!******************************************************************************
!  Module GEIA_MOD contains routines used to read and scale the GEIA fossil 
!  fuel emissions for NOx, CO, and hydrocarbons (bmy, 7/28/00, 1/13/10)
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_TOTCO2         : reads total CO2 scale factors (for FF NOx)
!  (2 ) READ_LIQCO2         : reads liquid CO2 scale factors (for CO, HC's)
!  (3 ) READ_TODX           : reads "time-of-day" scale factors for GEIA em's
!  (4 ) READ_GEIA_ASCII     : reads GEIA fossil fuel emissions from ASCII file
!  (5 ) READ_GEIA           : reads GEIA fossil fuel emissions from binary file
!  (6 ) READ_C3H8_C2H6_NGAS : reads C2H6 and C3H8 based on Natural Gas em's
!  (7 ) GET_DAY_INDEX       : Determines if today is a weekday, Sat., or Sun. 
!  (8 ) GET_IHOUR           : Returns index for the "time-of-day" scale factor
!  (9 ) TOTAL_FOSSIL_TG     : Computes total fossil fuel emissions in Tg
!
!  GEOS-CHEM modules referenced by geia_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs  
!  (3 ) file_mod.f      : Module containing file unit numbers and error checks
!  (4 ) grid_mod.f      : Module containing horizontal grid information
!
!  NOTES:
!  (1 ) Renamed original READ_GEIA to READ_GEIA_ASCII.  The new READ_GEIA
!        now reads fossil fuel emissions from the newer binary punch
!        file format. (bmy, 4/23/01)
!  (2 ) Added new routine TOTAL_FOSSIL_TG (bmy, 4/27/01)
!  (3 ) Bug fix: now read C2H6 from punch file correctly (bmy, 7/2/01)
!  (4 ) Added new routine: READ_C3H8_C2H6_NGAS.  Also updated comments. 
!        (bmy, 9/4/01)
!  (5 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (6 ) Now read scalefoss* files directly from the DATA_DIR filesystem,
!        instead of relying on symbolic links.  Also updated comments.
!        (bmy, 1/25/02)
!  (7 ) Eliminated obsolete code (bmy, 2/27/02)
!  (8 ) Routine READ_TODX now reads files from DATA_DIR/fossil_200104/, and
!        routine READ_GEIA_ASCII now reads files from DATA_DIR/fossil_obsolete.
!        This eliminates the need for symbolic file links. (bmy, 4/3/02)
!  (9 ) Updated comments (bmy, 5/28/02)
!  (10) Now references "file_mod.f" (bmy, 6/27/02)
!  (11) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/10/03)
!  (12) Now references "directory_mod.f" (bmy, 7/20/04)
!  (13) Now can read data from both GEOS and GCAP grids (bmy, 8/16/05)
!  (14) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (15) Modifications for 0.5 x 0.666 nested grids (yxw, dan, bmy, 11/6/08)
!  (16) Modifications for reading data on GEOS-4 1 x 1.25 grids (bmy, 1/13/10)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE READ_TOTCO2( SCALEYEAR, TOTCO2 )
!
!******************************************************************************
!  Subroutine READ_TOTCO2 reads in the scale factors (SCALEYEAR/1985) based 
!  on total CO2 emissions.  These are used to scale anthropogenic NOx 
!  emissions from 1985 to the present. (bmy, 9/13/00, 7/20/04)
!
!  NOTES:
!  (1 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (2 ) Make sure parentheses end before column 73 (bmy, 9/25/00)
!  (3 ) Now read the "scalefoss.tot*" files directly from the 
!        scalefoss_200202/ subdirectory of DATA_DIR, w/o relying on 
!        symbolic links.  Also updated comments. (bmy, 1/24/02)
!  (4 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Now write file name to stdout (bmy, 4/3/02)
!  (6 ) Now use IU_FILE instead of IUNIT.  Also reference IU_FILE and IOERROR
!        from "file_mod.f" (bmy, 6/27/02)
!  (7 ) Now use ENCODE to define CYEAR string for PGI/Linux (bmy, 9/29/03)
!  (8 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: SCALEYEAR
      REAL*4,  INTENT(OUT) :: TOTCO2(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: IU_FILE
      INTEGER              :: I, J, IX, JX, IOS
      CHARACTER(LEN=4  )   :: CYEAR
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_TOTCO2 begins here!
      !=================================================================
      WRITE( 6, 10 ) SCALEYEAR
 10   FORMAT( 'READ_TOTCO2: Year for Total CO2 scale factor: ', i4 )
      
      ! Define the file name and the file unit
      ! Now use ENCODE for PGI/F90 on Linux (bmy, 9/29/03)
#if   defined( LINUX ) 
      ENCODE( 4, '(i4)', CYEAR ) SCALEYEAR
#else 
      WRITE( CYEAR, '(i4)' ) SCALEYEAR
#endif

      FILENAME = TRIM( DATA_DIR )                  // 
     &           'scalefoss_200202/scalefoss.tot.' // 
     &           GET_RES_EXT() // '.' // CYEAR

      ! 1985 is the base year -- TOTCO2 = 1 in 1985!
      IF ( SCALEYEAR == 1985 ) THEN
         TOTCO2 = 1d0

      ELSE

         ! Echo filename to stdout
         WRITE( 6, 20 ) TRIM( FILENAME )
 20      FORMAT( 'READ_TOTCO2: Reading ', a )

         ! Find a free file LUN         
         IU_FILE = findFreeLUN()

         ! Open the file containing liquid CO2 scale factors
         OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &                  FORM='UNFORMATTED',    IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_totco2_file:1')

         ! Read the array dimensions IX, JX
         READ( IU_FILE, IOSTAT=IOS ) IX, JX
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_totco2_file:2')

         ! Read the data block of Liquid CO2 scale factors
         READ( IU_FILE, IOSTAT=IOS ) ( ( TOTCO2(I,J), I=1,IX ), J=1,JX )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_totco2_file:3')

         ! Close the file
         CLOSE( IU_FILE )

      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_TOTCO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_LIQCO2( SCALEYEAR, LIQCO2 )
!
!******************************************************************************
!  Subroutine READ_LIQCO2 reads in the scale factors (SCALEYEAR/1985) based 
!  on liquid CO2 emissions.  These are used to scale anthropogenic CO and
!  hydrocarbon emissions from 1985 to the present. (bmy, 9/13/00, 7/20/04)
!
!  NOTES:
!  (1 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (2 ) Make sure parentheses end before column 73 (bmy, 9/25/00)
!  (3 ) Now read the "scalefoss.liq*" files directly from the 
!        scalefoss_200202/ subdirectory of DATA_DIR, w/o relying on 
!        symbolic links.  Also updated comments. (bmy, 1/24/02)
!  (4 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Now write file name to stdout (bmy, 4/3/02)
!  (6 ) Now use IU_FILE instead of IUNIT.  Also reference IU_FILE and IOERROR
!        from "file_mod.f" (bmy, 6/27/02)
!  (7 ) Now use ENCODE to define CYEAR string for PGI/Linux (bmy, 9/29/03)
!  (8 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: SCALEYEAR
      REAL*4,  INTENT(OUT) :: LIQCO2(IIPAR,JJPAR)

      ! Local variables
      INTEGER              :: IU_FILE
      INTEGER              :: I, J, IX, JX, IOS
      CHARACTER(LEN=4  )   :: CYEAR
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_LIQCO2 begins here!
      !=================================================================
      WRITE( 6, 10 ) SCALEYEAR
 10   FORMAT( 'READ_LIQCO2: Year for Liquid CO2 scale factor: ', i4 )

      ! Define the file name and the file unit
      ! Now use ENCODE to define CYEAR string for Linux (bmy, 9/29/03)
#if   defined( LINUX ) 
      ENCODE( 4, '(i4)', CYEAR ) SCALEYEAR
#else
      WRITE( CYEAR, '(i4)' ) SCALEYEAR
#endif
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'scalefoss_200202/scalefoss.liq.' // 
     &           GET_RES_EXT() // '.' // CYEAR

      ! 1985 is the base year -- LIQCO2 = 1 in 1985!
      IF ( SCALEYEAR == 1985 ) THEN
         LIQCO2 = 1d0

      ELSE

         ! Echo filename to stdout
         WRITE( 6, 20 ) TRIM( FILENAME )
 20      FORMAT( 'READ_TOTCO2: Reading ', a )

         ! Find a free file LUN
         IU_FILE = findFreeLUN()

         ! Open the file containing liquid CO2 scale factors
         OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &                  FORM='UNFORMATTED',    IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_liqco2_file:1')

         ! Read the array dimensions IX, JX
         READ( IU_FILE, IOSTAT=IOS ) IX, JX
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_liqco2_file:2')

         ! Read the data block of Liquid CO2 scale factors
         READ( IU_FILE, IOSTAT=IOS ) ( ( LIQCO2(I,J), I=1,IX ), J=1,JX )
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_FILE,'read_liqco2_file:3')

         ! Close the file
         CLOSE( IU_FILE )
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_LIQCO2

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_TODX( TODN, TODH, TODB, SCNR89 )
!
!******************************************************************************
!  Subroutine READ_TODX reads the time-of-day emission scale factors and 
!  weekday/weekend scale factors for GEIA emissions. (bmy, 7/18/00, 1/13/09)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TODN   (REAL*8) : Time-of-day scale factor for NOx
!  (2 ) TODH   (REAL*8) : Time-of-day scale factor for hydrocarbons
!  (3 ) TODB   (REAL*8) : Time-of-day scale factor for biogenic species
!  (4 ) SCNR89 (REAL*8) : Weekday/Saturday/Sunday emission scale factors
!
!  NOTES:
!  (1 ) Copied from routine "anthroems.f" (bmy, 7/18/00)
!  (2 ) Added code for 1 x 1 GEOS grid (bmy, 8/7/00)
!  (3 ) Now use IOS /= 0 to trap both I/O errors and EOF (bmy, 9/13/00)
!  (4 ) Now read files directly from DATA_DIR/fossil_200104 subdirectory.
!        Also echo the file name to stdout.  Now reference DATA_DIR from
!        the "CMN_SETUP" header file. (bmy, 4/3/02)
!  (5 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Also deleted
!        obsolete code from April 2002 (bmy, 6/27/02)
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (7 ) Added space in the #ifdef block for 1 x 1.25 grid (bmy, 12/1/04)
!  (8 ) Now reads appropriate file for 0.5 x 0.666 nested grid simulations
!        (yxw, dan, bmy, 11/6/08)
!  (9) BUG: Jintai Lin reported issues with some of the data read here. The
!       sum of TODB is not 6, which points to a problem. However this is not
!       used in the code. (phs, 4/16/09)
!  (10) Now read data for GEOS-4 1 x 1.25 grid (bmy, 1/13/09)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      REAL*8, INTENT(OUT) :: TODH(6), TODN(6), TODB(6), SCNR89(3,3)

      ! Local variables
      INTEGER             :: I, J, K, IOS, IUNIT
      REAL*4              :: DUMMY(IIPAR,JJPAR)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_TODX begins here!
      !=================================================================

      ! Find a free file LUN
      IUNIT = findFreeLUN()
      
#if   defined( GRID4x5  )

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'fossil_200104/MELD_N96_70_HC'

      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_TODX: Reading ', a ) 

      ! Open the 4 x 5 file
      OPEN( IUNIT,  FILE=TRIM( FILENAME ), STATUS='OLD', 
     &              FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:1' )

      ! For 4 x 5 grid, skip over plume data
      DO K = 1, 2 
         READ( IUNIT, '(6e12.4)', IOSTAT=IOS ) 
     &        ( ( DUMMY(I+3,J+23), I=1,29 ), J=1,15 )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:2' )
      ENDDO

      ! Read time-of-day emission scale factors
      READ( IUNIT, '(6e12.4)', IOSTAT=IOS ) TODN, TODH, TODB
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:3' )
      
      ! Read weekday/saturday/sunday emission scale factors
      READ( IUNIT, '(2f6.3, f7.4)', IOSTAT=IOS ) 
     &     ( ( SCNR89(I,J), J=1,3 ), I=1,3 )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:4' )

      ! Close the 4 x 5 file
      CLOSE( IUNIT )

#elif defined( GRID2x25 )

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'fossil_200104/MELD2x25'

      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_TODX: Reading: ', a ) 

      ! Open the 2 x 2.5 file
      OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &             FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:5' )

      ! Read time-of-day scale factors
      READ( IUNIT, '(6E12.4)', IOSTAT=IOS ) TODN, TODH, TODB
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:6' )
      
      ! Read Weekday/Saturday/Sunday emission scale factors
      READ( IUNIT, '(3F7.4)', IOSTAT=IOS ) 
     &     ( ( SCNR89(I,J), J=1,3 ), I=1,3 )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:7' )      

      ! Close the file
      CLOSE( IUNIT )

#elif defined( GRID1x125 )

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'fossil_200104/MELD1x125'

      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_TODX: Reading: ', a ) 

      ! Open the 1 x 1.25 file
      OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &             FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:5' )

      ! Read time-of-day scale factors
      READ( IUNIT, '(6E12.4)', IOSTAT=IOS ) TODN, TODH, TODB
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:6' )
      
      ! Read Weekday/Saturday/Sunday emission scale factors
      READ( IUNIT, '(3F7.4)', IOSTAT=IOS ) 
     &     ( ( SCNR89(I,J), J=1,3 ), I=1,3 )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:7' )      

      ! Close the file
      CLOSE( IUNIT )

#elif defined( GRID1x1 )

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'fossil_200104/MELD1x1'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_TODX: Reading: ', a ) 

      ! Open the 1 x 1 file
      OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD', 
     &             FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:5' )

      ! Read time-of-day scale factors
      READ( IUNIT, '(6E12.4)', IOSTAT=IOS ) TODN, TODH, TODB
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:6' )
      
      ! Read Weekday/Saturday/Sunday emission scale factors
      READ( IUNIT, '(3F7.4)', IOSTAT=IOS ) 
     &     ( ( SCNR89(I,J), J=1,3 ), I=1,3 )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:7' )      

      ! Close the file
      CLOSE( IUNIT )

#elif defined( GRID05x0666 )

      ! Define the file name
      FILENAME = TRIM( DATA_DIR ) // 'fossil_200104/MELD05x0666'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_TODX: Reading: ', a )

      ! Open the 05x0666 file
      OPEN( IUNIT, FILE=TRIM( FILENAME ), STATUS='OLD',
     &             FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:5' )

      ! Read time-of-day scale factors
      READ( IUNIT, '(6E12.4)', IOSTAT=IOS ) TODN, TODH, TODB
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:6' )


      ! Read Weekday/Saturday/Sunday emission scale factors
      READ( IUNIT, '(3F7.4)', IOSTAT=IOS )
     &     ( ( SCNR89(I,J), J=1,3 ), I=1,3 )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_todx:7' )

      ! Close the file
      CLOSE( IUNIT )


#endif      

      ! Return to calling program
      END SUBROUTINE READ_TODX

!------------------------------------------------------------------------------

      SUBROUTINE READ_GEIA_ASCII( E_NOX,  E_CO,   E_ETHE, E_PRPE, 
     &                            E_C2H6, E_C3H8, E_ALK4, E_ACET, 
     &                            E_MEK,  E_SOX )
!
!******************************************************************************
!  Subroutine READ_GEIA_ASCII reads the anthropogenic GEIA emissions from 
!  from the old-style ASCII "merge file". (bmy, 7/18/00, 7/20/04)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_NOX  (REAL*4) : GEIA anthro NOx  (4 seasons,      2 levels)
!  (2 ) E_CO   (REAL*4) : GEIA anthro CO   (no seasonality, 1 level )
!  (3 ) E_ETHE (REAL*4) : GEIA anthro ETHE (no seasonality, 1 level )
!  (4 ) E_PRPE (REAL*4) : GEIA anthro PRPE (no seasonality, 1 level )
!  (5 ) E_C2H6 (REAL*4) : GEIA anthro C2H6 (no seasonality, 1 level )
!  (6 ) E_C3H8 (REAL*4) : GEIA anthro C3H8 (no seasonality, 1 level ) 
!  (7 ) E_ALK4 (REAL*4) : GEIA anthro ALK4 (no seasonality, 1 level )
!  (8 ) E_ACET (REAL*4) : GEIA anthro ACET (no seasonality, 1 level )
!  (9 ) E_MEK  (REAL*4) : GEIA anthro MEK  (no seasonality, 1 level )
!  (10) E_SOX  (REAL*4) : GEIA anthro SOx  (4 seasons,      2 levels)
!
!  NOTES:
!  (1 ) All arguments are optional.
!  (2 ) Copied from routine "anthroems.f" (bmy, 7/18/00)
!  (3 ) Now use IOS /= 0 to trap both I/O errors and EOF (bmy, 9/13/00)
!  (4 ) Renamed to READ_GEIA_ASCII, since this only reads the old-style
!        ASCII merge file (bmy, 4/23/01)
!  (5 ) Now read files directly from DATA_DIR/fossil_obsolete subdirectory.
!        Also echo the file name to stdout.  Now reference DATA_DIR from
!        the "CMN_SETUP" header file. (bmy, 4/3/02)
!  (6 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Also deleted
!        obsolete code from April 2002 (bmy, 6/27/02)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT), OPTIONAL :: E_NOX (IIPAR,JJPAR,4,2)
      REAL*4, INTENT(OUT), OPTIONAL :: E_CO  (IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_ETHE(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_PRPE(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_C2H6(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_C3H8(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_ALK4(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_ACET(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_MEK (IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_SOX (IIPAR,JJPAR,4,2)

      ! Local variables
      INTEGER                       :: IOS, IUNIT
     
      REAL*4                        :: T_NOX (IIPAR,JJPAR,4,2)
      REAL*4                        :: T_CO  (IIPAR,JJPAR    ) 
      REAL*4                        :: T_ETHE(IIPAR,JJPAR    )
      REAL*4                        :: T_PRPE(IIPAR,JJPAR    )
      REAL*4                        :: T_C2H6(IIPAR,JJPAR    ) 
      REAL*4                        :: T_C3H8(IIPAR,JJPAR    ) 
      REAL*4                        :: T_ALK4(IIPAR,JJPAR    )
      REAL*4                        :: T_ACET(IIPAR,JJPAR    ) 
      REAL*4                        :: T_MEK (IIPAR,JJPAR    )
      REAL*4                        :: T_SOX (IIPAR,JJPAR,4,2)

      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_GEIA_ASCII begins here!
      !=================================================================

      ! Define the file name 
      FILENAME = TRIM( DATA_DIR ) // 'fossil_obsolete/merge.' // 
     &           GET_RES_EXT()    // '_CTM'

#if   defined( GRID4x5 )
      ! For 4x5, the old ASCII file had the "_SASS" extension.
      ! This is historical baggage (bmy, 4/3/02)
      FILENAME = TRIM( FILENAME ) // '_SASS' 
#endif

      ! Find a free file LUN
      IUNIT = findFreeLUN()
      
      ! Echo the filename to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_GEIA_ASCII: Reading ', a )

      ! Open the GEIA emissions merge file
      OPEN( IUNIT,  FILE=TRIM( FILENAME ), STATUS='OLD', 
     &              FORM='FORMATTED',      IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_geia:1' )

      ! Read data into temporary arrays
      READ( IUNIT, '(7e10.3)', IOSTAT=IOS )
     &     T_NOX,  T_CO,   T_ETHE, T_PRPE, T_C2H6,
     &     T_C3H8, T_ALK4, T_ACET, T_MEK,  T_SOX
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IUNIT, 'read_geia:2' )

      ! Close the merge file
      CLOSE( IUNIT )

      ! Assign data from temporary arrays into optional arguments
      IF ( PRESENT( E_NOX  ) ) E_NOX  = T_NOX
      IF ( PRESENT( E_CO   ) ) E_CO   = T_CO
      IF ( PRESENT( E_ETHE ) ) E_ETHE = T_ETHE
      IF ( PRESENT( E_PRPE ) ) E_PRPE = T_PRPE
      IF ( PRESENT( E_C2H6 ) ) E_C2H6 = T_C2H6
      IF ( PRESENT( E_C3H8 ) ) E_C3H8 = T_C3H8
      IF ( PRESENT( E_ALK4 ) ) E_ALK4 = T_ALK4
      IF ( PRESENT( E_ACET ) ) E_ACET = T_ACET
      IF ( PRESENT( E_MEK  ) ) E_MEK  = T_MEK
      IF ( PRESENT( E_SOX  ) ) E_SOX  = T_SOX

      ! Return to calling program
      END SUBROUTINE READ_GEIA_ASCII

!------------------------------------------------------------------------------

      SUBROUTINE READ_GEIA( E_NOX,  E_CO,   E_ALK4, E_ACET, E_MEK,
     &                      E_PRPE, E_C3H8, E_C2H6, E_ETHE, E_SOX )
!
!******************************************************************************
!  Subroutine READ_GEIA reads the anthropogenic GEIA emissions 
!  from a binary punch file. (bmy, 4/23/01, 8/16/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_NOX  (REAL*4) : GEIA anthro NOx  (4 seasons,      2 levels)
!  (2 ) E_CO   (REAL*4) : GEIA anthro CO   (no seasonality, 1 level )
!  (3 ) E_ALK4 (REAL*4) : GEIA anthro ALK4 (no seasonality, 1 level )
!  (4 ) E_ACET (REAL*4) : GEIA anthro ACET (no seasonality, 1 level )
!  (5 ) E_MEK  (REAL*4) : GEIA anthro MEK  (no seasonality, 1 level )
!  (6 ) E_PRPE (REAL*4) : GEIA anthro PRPE (no seasonality, 1 level )
!  (7 ) E_C3H8 (REAL*4) : GEIA anthro C3H8 (no seasonality, 1 level ) 
!  (8 ) E_C2H6 (REAL*4) : GEIA anthro C2H6 (no seasonality, 1 level )
!  (9 ) E_ETHE (REAL*4) : GEIA anthro ETHE (no seasonality, 1 level )
!  (10) E_SOX  (REAL*4) : GEIA anthro SOx  (4 seasons,      2 levels)
!
!  NOTES:
!  (1 ) Now reads from binary punch file format.  This is more convenient,
!        and is readable directly into GAMAP.  Read directly from the
!        DATA_DIR/fossil_200104/ subdirectory. (bmy, 4/23/01)
!  (2 ) Bug fix: T_C2H6 was being overwritten with T_ETHE.  This has now
!        been corrected. (bmy, 7/2/01)
!  (3 ) Now only read emissions for tracers whose keywords have been passed
!        (bmy, 9/6/01)
!  (4 ) Now write file name to stdout (bmy, 4/3/02)
!  (5 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/14/03)
!  (6 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (7 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT, READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT), OPTIONAL :: E_NOX (IIPAR,JJPAR,4,2)
      REAL*4, INTENT(OUT), OPTIONAL :: E_CO  (IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_ALK4(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_ACET(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_MEK (IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_PRPE(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_C3H8(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_C2H6(IIPAR,JJPAR    ) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_ETHE(IIPAR,JJPAR    )
      REAL*4, INTENT(OUT), OPTIONAL :: E_SOX (IIPAR,JJPAR,4,2)

      ! Local variables
      INTEGER                       :: L
      REAL*4                        :: ARRAY(IIPAR,JJPAR,2)
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_GEIA begins here!
      !=================================================================

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR )                          //
     &           'fossil_200104/merge_nobiofuels.'         //
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
      
      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_GEIA: Reading ', a )

      !=================================================================
      ! Read NOx (tracer #1): 4 seasons, 2 levels
      !=================================================================
      IF ( PRESENT( E_NOX ) ) THEN

         ! Read winter NOx (DJF)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,  
     &                    -744d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )
         
         DO L = 1, 2
            E_NOX(:,:,1,L) = ARRAY(:,:,L)
         ENDDO
      
         ! Read spring NOx (MAM)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,  
     &                    1416d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2
            E_NOX(:,:,2,L) = ARRAY(:,:,L)
         ENDDO

         ! Read summer NOx (JJA)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,  
     &                    3624d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2
            E_NOX(:,:,3,L) = ARRAY(:,:,L)
         ENDDO

         ! Read autumn NOx (SON)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,  
     &                    5832d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2
            E_NOX(:,:,4,L) = ARRAY(:,:,L) 
         ENDDO
      ENDIF
        
      !=================================================================
      ! Read CO (tracer #4): aseasonal
      !=================================================================
      IF ( PRESENT( E_CO ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    4,  
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_CO(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read ALK4 (tracer #5): aseasonal
      !=================================================================
      IF ( PRESENT( E_ALK4 ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    5,  
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_ALK4(:,:) = ARRAY(:,:,1)
      ENDIF
 
      !=================================================================
      ! Read ACET (tracer #9): aseasonal
      !=================================================================
      IF ( PRESENT( E_ACET ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    9,  
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_ACET(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read MEK (tracer #10): aseasonal
      !=================================================================
      IF ( PRESENT( E_MEK ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    10, 
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_MEK(:,:)  = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read PRPE (tracer #18): aseasonal
      !=================================================================
      IF ( PRESENT( E_PRPE ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    18, 
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_PRPE(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read C3H8 (tracer #19): aseasonal
      !=================================================================
      IF ( PRESENT( E_C3H8 ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    19, 
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         E_C3H8(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read C2H6 (tracer #20): aseasonal
      !=================================================================
      IF ( PRESENT( E_C2H6 ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    21, 
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )
         
         E_C2H6(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read ETHE (tracer #26): aseasonal
      !=================================================================
      IF ( PRESENT( E_ETHE ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    26, 
     &                    0d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )
         
         E_ETHE(:,:) = ARRAY(:,:,1)
      ENDIF

      !=================================================================
      ! Read SOx (tracer #27): 4 seasons, 2 levels
      !=================================================================
      IF ( PRESENT( E_SOX ) ) THEN 

         ! Read winter SOx (DJF)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
     &                    -744d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2 
            E_SOX(:,:,1,L) = ARRAY(:,:,L)
         ENDDO

         ! Read spring SOx (MAM)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
     &                    1416d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )
      
         DO L = 1, 2 
            E_SOX(:,:,2,L) = ARRAY(:,:,L)
         ENDDO

         ! Read summer SOx (JJA)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
     &                    3624d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2 
            E_SOX(:,:,3,L) = ARRAY(:,:,L)
         ENDDO

         ! Read autumn SOx (SON)
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 27, 
     &                    5832d0,    IIPAR,     JJPAR,     
     &                    2,         ARRAY,     QUIET=.TRUE. )

         DO L = 1, 2 
            E_SOX(:,:,4,L) = ARRAY(:,:,L)
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_GEIA

!------------------------------------------------------------------------------

      SUBROUTINE READ_C3H8_C2H6_NGAS( E_C3H8, E_C2H6 )
!
!******************************************************************************
!  Subroutine READ_C3H8_C2H6_NGAS reads the anthropogenic C3H8 and C2H6
!  emissions, which are scaled from Natural Gas (CH4) (bmy, 9/6/01, 8/16/05)
!
!  Emissions files are from Yaping Xiao (9/01)  Their path names are:
!     /data/ctm/GEOS_2x2.5/C3H8_C2H6_200109/C3H8_C2H6_ngas.geos.2x25
!     /data/ctm/GEOS_4x5/C3H8/C2H6_200109/C3H8_C2H6_ngas.geos.4x5
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_C3H8 (REAL*4) : Anthro C3H8 scaled from CH4 (aseasonal, 1 level) 
!  (2 ) E_C2H6 (REAL*4) : Anthro C2H6 scaled from CH4 (aseasonal, 1 level)
!
!  NOTES:
!  (1 ) Adapted from READ_GEIA (bmy, 9/6/01)
!  (2 ) Now echo filename to standard output (bmy, 1/25/02)
!  (3 ) Now call READ_BPCH2 with QUIET=.TRUE. (bmy, 3/11/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (5 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT), OPTIONAL :: E_C3H8(IIPAR,JJPAR) 
      REAL*4, INTENT(OUT), OPTIONAL :: E_C2H6(IIPAR,JJPAR) 

      ! Local variables
      INTEGER                       :: I, J, L, S
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1) 
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_GEIA begins here!
      !=================================================================

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR )                         //
     &           'C3H8_C2H6_200109/C3H8_C2H6_ngas.'       // 
     &           GET_NAME_EXT_2D() // '.' // GET_RES_EXT() 
      
      ! Echo filename to std output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'READ_C3H8_C2H6_NGAS: Reading ', a )

      ! Read C3H8 (tracer #19): aseasonal
      IF ( PRESENT( E_C3H8 ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 19, 
     &                    0d0,       IIPAR,     JJPAR,     
     &                    1,         ARRAY,     QUIET=.TRUE. )
         
         E_C3H8(:,:) = ARRAY(:,:,1)
      ENDIF

      ! Read C2H6 (tracer #21): aseasonal
      IF ( PRESENT( E_C2H6 ) ) THEN
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 21, 
     &                    0d0,       IIPAR,     JJPAR,     
     &                    1,         ARRAY,     QUIET=.TRUE. )

         E_C2H6(:,:) = ARRAY(:,:,1)
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_C3H8_C2H6_NGAS

!------------------------------------------------------------------------------

      FUNCTION GET_DAY_INDEX( NTAU ) RESULT( JSCEN )
!
!*****************************************************************************
!  Function GET_DAY_INDEX returns the day index (Saturday/Sunday/Weekday)
!  for the given time.  This is used to scale GEIA emissions. (bmy, 7/28/00)
!
!  Arguments as Input:
!  ===========================================================================
!  (1) NTAU  (INTEGER) : Integral hours since 0h, 1 Jan 1985
!
!  Return value:
!  ===========================================================================
!  (1) JSCEN (INTEGER) : Flag for Saturday (1), Sunday (2), or Weekday (3)
!
!  NOTES:
!  (1) Scale factors for Saturday/Sunday/Weekday must average out to 1!
!*****************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: NTAU

      ! Local variables
      INTEGER             :: NDAY

      ! Return value
      INTEGER             :: JSCEN

      !=================================================================
      ! GET_DAY_INDEX begins here!
      !=================================================================

      ! NDAY is the day of the week
      NDAY = NTAU / 24

      ! 1 Jan 1980 and 1 Jan 1985 were both Tuesdays, so NDAY mod 7 = 4 is a 
      ! Saturday and NDAY mod 7 = 5 is a Sunday (bmy, 3/23/98)
      SELECT CASE ( MOD( NDAY, 7 ) ) 

         ! Saturday
         CASE ( 4 )
            JSCEN = 1
            
         ! Sunday
         CASE ( 5 )
            JSCEN = 2

         ! Weekday
         CASE DEFAULT
            JSCEN = 3

      END SELECT

      ! Return to calling program
      END FUNCTION GET_DAY_INDEX

!------------------------------------------------------------------------------

      FUNCTION GET_IHOUR( I ) RESULT( IHOUR )
!
!******************************************************************************
!  Function GET_IHOUR returns the index for the TODH, TODN, TODB scale 
!  factors which are read by subroutine GET_TODX. (bmy, 7/28/00, 4/23/01)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : Grid box longitude index
!
!  Return value
!  ============================================================================
!  (1 ) IHOUR  (INTEGER) : Index for scale factor arrays
!
! !REMARKS:
!  For now use GET_LOCALTIME( I, 1, 1 ) which will be independent of J and L
!  for a pure cartesian grid.  This may need to be revisited once G-C is
!  interfaced into a GCM.
!
!  NOTES:
!  (1 ) COMPUTE IHOUR TO DETERMINE TIME OF DAY FACTOR FOR EMISSIONS
!        THIS TOFDAY IS IN GREENWICH TIME(GMT); WE NEED TO CHANGE IT TO
!        (I,J)BOX TIME=XLOCTM
!        TOFDAY is GMT at the BEGINNING of the time step.
!        (TOFDAY and NTAU do not refer to the end of time step until they
!        are updated at line 300 of the main driver, just before the
!        diagnostics are written.)
!        The 0.001 is added to remove roundoff amibuity when timestep
!        is exactly on boundary for emissions change.
!        1hr changed NCHEM->NDYN in following line
!  (2 ) For GEOS-CTM, NDYN is in minutes, NDYN/60 is in hours (bmy, 2/26/98)
!  (3 ) Make sure 0 <= XLOCTM < 24, to avoid subscript errors (bmy, 6/11/98)
!  (4 ) Middle of time step is between 10pm-2am when IHOUR = 1
!  (5 ) Updated comments (bmy, 4/23/01)
!  (6 ) Now use function GET_LOCALTIME from the new "time_mod.f".  Remove
!        IREF, NDYN, TOFDAY, DISIZE from arg list.  Add I to the arg list.
!        Removed XLOCTM variable. (bmy, 2/10/03)
!  (7 ) Modified to use NINT instead of INT to calculate the local time
!        (ccc, 4/15/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_LOCALTIME

      ! Arguments
      INTEGER, INTENT(IN) :: I

      ! Return value
      INTEGER             :: IHOUR
      
      !=================================================================
      ! GET_IHOUR begins here!
      !=================================================================

      ! IHOUR ranges from 1-6 
      ! Modified to use NINT instead of INT (ccc, 4/15/09)
      IHOUR = NINT( ( GET_LOCALTIME( I, 1, 1 ) ) / 4 ) + 1
      IF ( IHOUR == 7 ) IHOUR = 1

      ! Return to calling program
      END FUNCTION GET_IHOUR

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_FOSSIL_TG( FFARRAY, IX,   JX, LX, 
     &                            MOLWT,   NAME, NSEASON )
!
!******************************************************************************
!  Subroutine TOTAL_FOSSIL_TG prints the amount of biomass burning
!  emissions that are emitted each month in Tg or Tg C. (bmy, 4/27/01, 2/4/03)
!  
!  Arguments as Input:
!  ============================================================================
!  (1  ) FFARRAY  (REAL*8 ) : Fossil Fuel CO emissions [molec (C)/cm2/month]
!  (2-4) IX,JX,LX (INTEGER) : Dimensions of FFARRAY 
!  (5  ) MOLWT    (REAL*8 ) : Molecular wt [kg/mole] for the given tracer
!  (6  ) NAME     (REAL*8 ) : Tracer name
!  (7  ) NSEASON  (INTEGER) : Number of the season, for seasonal NOx/SOX
!
!  NOTES:
!  (1) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!      Bryan Duncan (bnd@io.harvard.edu), and Daniel Jacob (djj@io.harvard.edu)
!  (2) Now replace DXYP(J)*1d4 with routine GET_AREA_CM2 from "grid_mod.f".
!       (bmy, 2/4/03)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER,           INTENT(IN) :: IX, JX, LX
      INTEGER, OPTIONAL, INTENT(IN) :: NSEASON
      REAL*8,            INTENT(IN) :: FFARRAY(IX,JX,LX) 
      REAL*8,            INTENT(IN) :: MOLWT
      CHARACTER(LEN=*),  INTENT(IN) :: NAME

      ! Local variables
      INTEGER                       :: I, J, L
      REAL*8                        :: TOTAL, A_CM2
      CHARACTER(LEN=6)              :: UNIT

      !=================================================================
      ! TOTAL_FOSSIL_TG begins here!
      !=================================================================

      ! Initialize summing variable
      TOTAL = 0d0

      DO L = 1, LX
      DO J = 1, JX
      DO I = 1, IX

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( I, J, L )

         TOTAL = TOTAL + FFARRAY(I,J,L) * A_CM2 * ( MOLWT/ 6.023d23 )
      ENDDO
      ENDDO
      ENDDO

      IF ( PRESENT( NSEASON ) ) THEN
         
         ! Total for each season
         SELECT CASE( NSEASON )

            ! DJF is 90 days long
            CASE ( 1 )
               TOTAL = TOTAL * 1d-9 * 90d0 * 86400d0

            ! MAM, JJA are 92 days long
            CASE ( 2, 3 )
               TOTAL = TOTAL * 1d-9 * 92d0 * 86400d0

            ! SON is 91 days long
            CASE ( 4 )
               TOTAL = TOTAL * 1d-9 * 91d0 * 86400d0
         END SELECT

      ELSE

         ! Convert from kg --> Tg for aseasonal emissions
         TOTAL = TOTAL * 1d-9 * 365.25d0 * 86400d0

      ENDIF

      ! Define unit string
      SELECT CASE( TRIM( NAME ) ) 
         CASE ( 'NO' ) 
            UNIT = '[Tg N]'
         CASE ( 'CO', 'CH2O'  )
            UNIT = '[Tg  ]'
         CASE DEFAULT
            UNIT = '[Tg C]'
      END SELECT

      ! Write totals
      IF ( PRESENT( NSEASON ) ) THEN

         ! Seasonal
         WRITE( 6, 100 ) NAME, TOTAL, UNIT, NSEASON
 100     FORMAT( 'Total Anthropogenic ', a4, ': ', f9.3, 1x, a,
     &           ' Season =', i3 )

      ELSE

         ! Aseasonal
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Total Anthropogenic ', a4, ': ', f9.3, 1x, a )

      ENDIF

      ! Return to calling program
      END SUBROUTINE TOTAL_FOSSIL_TG

!------------------------------------------------------------------------------

      END MODULE GEIA_MOD
