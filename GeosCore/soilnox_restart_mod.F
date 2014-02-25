!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: soilnox_restart_mod
!
! !DESCRIPTION: Module SOILNOX\_RESTART\_MOD contains variables and routines
!  used to read and write GEOS-CHEM Soil NOx restart files, which contain the
!  following:
!    DRYPERIOD - time since soil moisture increased by 0.01 (hours), 
!    PFACTOR   - If soil pulsing, pulse factor from prev. timestep (unitless)
!    This code was modified from restart\_mod.F
!\\
!\\
! !INTERFACE: 
!
      MODULE SOILNOX_RESTART_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: MAKE_SOILNOX_RESTART
      PUBLIC  :: READ_SOIL_RESTART
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_soilnox_restart
!
! !DESCRIPTION: Subroutine MAKE\_SOILNOX\_RESTART creates GEOS-CHEM restart
!  files of soilnox variables in binary punch file format. (rch, 10/15/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_SOILNOX_RESTART( am_I_Root, Input_Opt, YYYYMMDD,  
     &                                 HHMMSS,    TAU,       RC        )
! 
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR
      USE BPCH2_MOD,          ONLY : GET_MODELNAME
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE FILE_MOD,           ONLY : IOERROR
      USE GRID_MOD,           ONLY : GET_XOFFSET
      USE GRID_MOD,           ONLY : GET_YOFFSET
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE inquireMod,         ONLY : findfreeLUN
      USE TIME_MOD,           ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! Year-Month-Day
      INTEGER,        INTENT(IN)  :: HHMMSS      ! Hour-Min-Sec
      REAL*8,         INTENT(IN)  :: TAU         ! TAU value corresponding
                                                 !  to YYYYMMDD, HHMMSS
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  31 Oct 2012 - R. Yantosca - Now define a title string
!  22 Aug 2013 - R. Yantosca - Add am_I_Root, Input_Opt, RC arguments
!  22 Aug 2013 - R. Yantosca - Now get soil NOx file path from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: LPRT
      INTEGER              :: I0, IOS,  J0, N, IU_SOIL
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      REAL*4               :: DRES(IGLOB, JGLOB) 
      INTEGER              :: HALFPOLAR
      INTEGER, PARAMETER   :: CENTER180 = 1
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! MAKE_SOILNOX_RESTART begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Lon & lat sizes
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Call GET_HALFPOLAR to return the proper value
      ! for either GCAP or GEOS grids (bmy, 6/28/05)
      HALFPOLAR = GET_HALFPOLAR()

      ! Get the nested-grid offsets
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Top-of-file title string
      TITLE     = 'GEOS-CHEM Restart File: Soil NOx quantities'

      ! Print debug output?
      LPRT      = ( am_I_Root .and. Input_Opt%LPRT )

      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      ! Copy the soil NOx restart filename template into a local variable
      FILENAME = TRIM( Input_Opt%SOIL_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_SOILNOX_RESTART: Writing: ', a )

      CALL FLUSH(6)

      ! Find next free file unit
      IU_SOIL = findFreeLUN()

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_SOIL, FILENAME, TITLE )

      !=================================================================
      ! Write each variable to the restart file
      !=================================================================   

      !---------------------------
      ! DRY PERIOD [hours]
      !---------------------------
      UNIT     = 'hours'
      CATEGORY = 'RST-SOIL'
      N        = 1

      ! Save to disk
      CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  N,
     &            UNIT,      TAU,       TAU,       RESERVED,   
     &            IIPAR,     JJPAR,     1,         I0+1,            
     &            J0+1,      1,         DRYPERIOD            )

      !---------------------------
      ! PFACTOR [unitless]
      !---------------------------
      UNIT     = 'unitless'
      CATEGORY = 'RST-SOIL'
      N        = 2

      ! Write to disk
      CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  N,
     &            UNIT,      TAU,       TAU,       RESERVED,   
     &            IIPAR,     JJPAR,     1,         I0+1,            
     &            J0+1,      1,         PFACTOR              )
      
      !---------------------------
      ! Previous GWET [unitless]
      !---------------------------
      UNIT     = 'unitless'
      CATEGORY = 'RST-SOIL'
      N        = 3

      ! Save to disk
      CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  N,
     &            UNIT,      TAU,       TAU,       RESERVED,   
     &            IIPAR,     JJPAR,     1,         I0+1,            
     &            J0+1,      1,         GWET_PREV            )

      !---------------------------
      ! N_RESERVOIR [ng N/m2]
      !---------------------------
      DRES     = DEP_RESERVOIR
      UNIT     = 'ng N/m2'
      CATEGORY = 'RST-SOIL'
      N        = 4

      ! Save to disk
      CALL BPCH2( IU_SOIL,   MODELNAME, LONRES,    LATRES,    
     &            HALFPOLAR, CENTER180, CATEGORY,  N,
     &            UNIT,      TAU,       TAU,       RESERVED,   
     &            IIPAR,     JJPAR,     1,         I0+1,            
     &            J0+1,      1,         DRES                 )

      ! Close file
      CLOSE( IU_SOIL )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG('### MAKE_SOILNOX_RESTART: wrote file')

      ! Return to calling program
      END SUBROUTINE MAKE_SOILNOX_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_soil_restart
!
! !DESCRIPTION: Subroutine READ\_SOIL\_RESTART initializes GEOS-CHEM Soil NOx 
!  parameters (binary punch file format)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SOIL_RESTART( am_I_Root, Input_Opt, 
     &                              YYYYMMDD,  HHMMSS,    RC ) 
! 
! !USES:
!
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE CMN_SIZE_MOD
      USE COMMSOIL_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE inquireMod,         ONLY : findfreeLUN
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: YYYYMMDD    ! Year-Month-Day
      INTEGER,        INTENT(IN)  :: HHMMSS      ! Hour-Min-Sec
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  22 Aug 2013 - R. Yantosca - Add am_I_Root, Input_Opt, RC arguments
!  22 Aug 2013 - R. Yantosca - Now get soil NOx file path from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: LSPLIT, LPRT
      INTEGER             :: I, IOS, J, L, N
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,        NJ,       NL
      INTEGER             :: IFIRST,    JFIRST,   LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      INTEGER             :: IU_SOIL
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      REAL*4              :: DRES(IGLOB, JGLOB) 
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_SOIL_RESTART begins here!
      !
      ! Open restart file and read top-of-file header
      !=================================================================

      ! Assume success
      RC       = GIGC_SUCCESS

      ! Get logicals from Input_Opt
      LSPLIT   = Input_Opt%LSPLIT
      LPRT     = ( am_I_Root .and. Input_Opt%LPRT )

      ! Copy soil NOx restart filename template to a local variable
      FILENAME = TRIM( Input_Opt%SOIL_RST_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '-', 79 )
      WRITE( 6, '(a,/)' ) 'SOIL NOx Restart'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_SOIL_RESTART: Reading ', a )

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

      END SUBROUTINE READ_SOIL_RESTART
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS makes sure that the dimensions of
!  the restart file extend to cover the entire grid. (bmy, 6/25/02, 10/15/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_DIMENSIONS( NI, NJ ) 
! 
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NI ! Number of longitudes read from restart file
      INTEGER, INTENT(IN) :: NJ ! Number of latitudes  read from restart file
!
! !REVISION HISTORY:
!  (1 ) Added to "restart_mod.f".  Now no longer allow initialization with 
!        less than a globally-sized data block. (bmy, 6/25/02)
!  (2 ) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
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
      
      END SUBROUTINE CHECK_DIMENSIONS
!EOC
      END MODULE SOILNOX_RESTART_MOD
