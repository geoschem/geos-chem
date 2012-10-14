! $Id: soilnox_restart_mod.f,v 1 2009/10/15 14:14:45 rch Exp $
      MODULE SOILNOX_RESTART_MOD
!
!******************************************************************************
!  Module RESTART_MOD contains variables and routines which are used to read
!  and write GEOS-CHEM Soil NOx restart files, which contain the following:
!  DRYPERIOD - time since soil moisture increased by 0.01 (hours), 
!  PFACTOR   - If soil pulsing, pulse factor from prev. timestep (unitless)
!  This code was modified from restart_mod.f
!
!  Module Variables:
!  ============================================================================
!  (1 ) INPUT_RESTART_FILE   : Full path name of the restart file to be read
!  (2 ) OUTPUT_RESTART_FILE  : Full path name (w/ tokens!) of output file
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_SOILNOX_RESTART : Writes restart file to disk 
!  (2 ) READ_RESTART_FILE    : Reads restart file from disk 
!  (4 ) CHECK_DIMENSIONS     : Ensures that restart file contains global data
!  (6 ) CHECK_DATA_BLOCKS    : Makes sure we have read in data for each tracer
!
!  GEOS-CHEM modules referenced by restart_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f          : Module w/ routines for binary punch file I/O
!  (2 ) error_mod.f          : Module w/ NaN and other error check routines
!  (3 ) file_mod.f           : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f           : Module w/ horizontal grid information
!  (5 ) logical_mod.f        : Module w/ GEOS-CHEM logical switches
!  (6 ) time_mod.f           : Module w/ routines for computing time & date
!  (7 ) tracer_mod.f         : Module w/ GEOS-CHEM tracer array STT etc.
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "restart_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC  :: MAKE_SOILNOX_RESTART
      PUBLIC  :: READ_SOIL_RESTART
  

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_SOILNOX_RESTART( YYYYMMDD, HHMMSS, TAU )
!
!******************************************************************************
!  Subroutine MAKE_SOILNOX_RESTART creates GEOS-CHEM restart files of soilnox  
!  variables in binary punch file format. (rch, 10/15/2009)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!  (3 ) TAU      : GEOS-CHEM TAU value corresponding to YYYYMMDD, HHMMSS
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : BPCH2,         GET_MODELNAME
      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR, OPEN_BPCH2_FOR_WRITE
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IOERROR
      USE GRID_MOD,    ONLY : GET_XOFFSET,   GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE inquireMod,  ONLY : findfreeLUN
      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
      REAL*8,  INTENT(IN)  :: TAU

      ! Local Variables      
      INTEGER              :: I0, IOS,  J0, N, IU_SOIL
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      REAL*4               :: DRES(IGLOB, JGLOB) 
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1
      
      CHARACTER(LEN=255)   :: OUTPUT_RESTART_FILE 
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! MAKE_RESTART_FILE begins here!
      !=================================================================

      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Call GET_HALFPOLAR to return the proper value
      ! for either GCAP or GEOS grids (bmy, 6/28/05)
      HALFPOLAR = GET_HALFPOLAR()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      OUTPUT_RESTART_FILE = 'restart.soilnox.YYYYMMDDhh'

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_SOILNOX_RESTART_FILE: Writing: ', a )

      CALL FLUSH(6)

      IU_SOIL = findFreeLUN()

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_SOIL, FILENAME, TITLE )

      !=================================================================
      ! Write each variable to the restart file
      !=================================================================   

      ! DRY PERIOD 
      UNIT     = 'hours'
      CATEGORY = 'RST-SOIL'
      N = 1

         ! SAVE DRYPERIOD (HOURS)
         CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      TAU,       TAU,       RESERVED,   
     &               IIPAR,     JJPAR,     1,     I0+1,            
     &               J0+1,      1,         DRYPERIOD )


      ! PFACTOR 
      UNIT     = 'unitless'
      CATEGORY = 'RST-SOIL'
      N = 2


         ! SAVE PFACTOR
         CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      TAU,       TAU,       RESERVED,   
     &               IIPAR,     JJPAR,     1,     I0+1,            
     &               J0+1,      1,         PFACTOR )
      

      !GWET PREV
      UNIT     = 'unitless'
      CATEGORY = 'RST-SOIL'
      N = 3


         ! SAVE SOIL MOSITURE FROM PREVIOUS TIMESTEP
         CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      TAU,       TAU,       RESERVED,   
     &               IIPAR,     JJPAR,     1,     I0+1,            
     &               J0+1,      1,         GWET_PREV )

       WRITE(6,*) GWET_PREV(45,24)
       CALL FLUSH(6)



      ! N_RESERVOIR
      DRES = DEP_RESERVOIR
      UNIT     = 'ng N/m2'
      CATEGORY = 'RST-SOIL'
      N = 4

         ! SAVE N_Reservior
         CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      TAU,       TAU,       RESERVED,   
     &               IIPAR,     JJPAR,     1,     I0+1,            
     &               J0+1,      1,         DRES )



      ! Close file
      CLOSE( IU_SOIL )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_SOILNOX_RESTART: wrote file')

      ! Return to calling program
      END SUBROUTINE MAKE_SOILNOX_RESTART

!------------------------------------------------------------------------------

      SUBROUTINE READ_SOIL_RESTART( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_SOIL_RESTART initializes GEOS-CHEM Soil NOx 
!  parameters (binary punch file format) 
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  NOTES:
!  (1 )   
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IOERROR
      USE LOGICAL_MOD, ONLY : LSPLIT,      LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE inquireMod,  ONLY : findfreeLUN
      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ, NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      INTEGER             :: IU_SOIL
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      REAL*4               :: DRES(IGLOB, JGLOB) 
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=255)  :: INPUT_RESTART_FILE  
      !=================================================================
      ! READ_RESTART_FILE begins here!
      !=================================================================

  
      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================
      INPUT_RESTART_FILE = 'restart.soilnox.YYYYMMDDhh'
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_RESTART_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '-', 79 )
      WRITE( 6, '(a,/)' ) 'SOIL NOx Restart'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RESTART_FILE: Reading ', a )

      IU_SOIL = findFreeLUN()

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_SOIL, FILENAME )
      
      !=================================================================
      ! Read soil nox info 
      !=================================================================

      ! Read Dry Period ----------------------------------

      READ( IU_SOIL, IOSTAT=IOS )       
     &      MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:1')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NI,       NJ,    NL, IFIRST, JFIRST, LFIRST,
     &      NSKIP

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:2')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &    ( (   DRYPERIOD(I,J), I=1,NI ), J=1,NJ  )

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:3')

      ! Make sure array dimensions are of global size
      ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
      CALL CHECK_DIMENSIONS( NI, NJ )

      WRITE(6,*) 'DONE DRY PERIOD'

      ! Read PFACTOR    ----------------------------------

      READ( IU_SOIL, IOSTAT=IOS )       
     &      MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:4')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NI,       NJ,    NL, IFIRST, JFIRST, LFIRST,
     &      NSKIP

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:5')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &    ( (   PFACTOR(I,J), I=1,NI ), J=1,NJ  )

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:6')

      ! Make sure array dimensions are of global size
      ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
      CALL CHECK_DIMENSIONS( NI, NJ )

      WRITE(6,*) 'DONE PFACT'

      ! Read GWET_PREV  ----------------------------------

      READ( IU_SOIL, IOSTAT=IOS )       
     &      MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:7')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NI,       NJ,    NL, IFIRST, JFIRST, LFIRST,
     &      NSKIP

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:8')


      READ( IU_SOIL, IOSTAT=IOS ) 
     &    ( (   GWET_PREV(I,J), I=1,NI ), J=1,NJ  )


      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:9')

      ! Make sure array dimensions are of global size
      ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
      CALL CHECK_DIMENSIONS( NI, NJ )

      WRITE(6,*) 'DONE GWET PREV'


      ! Read N_RESERVOIR from prev. ts ----------------------------------
 
      READ( IU_SOIL, IOSTAT=IOS )       
     &      MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:10')

      READ( IU_SOIL, IOSTAT=IOS ) 
     &      CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &      NI,       NJ,    NL, IFIRST, JFIRST, LFIRST,
     &      NSKIP

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:11')


      READ( IU_SOIL, IOSTAT=IOS ) 
     &    ( (   DRES(I,J), I=1,NI ), J=1,NJ  )

      DEP_RESERVOIR = DRES

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_SOIL,'read_soil_restart:12')

      ! Make sure array dimensions are of global size
      ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
      CALL CHECK_DIMENSIONS( NI, NJ )


      WRITE(6,*) 'DONE NRES'
        
      ! Return to calling program
      END SUBROUTINE READ_SOIL_RESTART


!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIMENSIONS( NI, NJ ) 
!
!******************************************************************************
!  Subroutine CHECK_DIMENSIONS makes sure that the dimensions of the
!  restart file extend to cover the entire grid. (bmy, 6/25/02, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NI (INTEGER) : Number of longitudes read from restart file
!  (2 ) NJ (INTEGER) : Number of latitudes  read from restart file
!
!  NOTES:
!  (1 ) Added to "restart_mod.f".  Now no longer allow initialization with 
!        less than a globally-sized data block. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER, INTENT(IN) :: NI, NJ


      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check longitude dimension: NI must equal IIPAR
      IF ( NI /= IIPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
         WRITE( 6, '(a)' ) 'Wrong number of longitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Error check latitude dimension: NJ must equal JJPAR
      IF ( NJ /= JJPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in restart file!'
         WRITE( 6, '(a)' ) 'Wrong number of latitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      

      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS

!------------------------------------------------------------------------------
     

 
      ! End of module
      END MODULE SOILNOX_RESTART_MOD
