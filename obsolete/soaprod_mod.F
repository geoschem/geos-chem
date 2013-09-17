! $Id: soaprod_mod.f,v 1.1 2010/03/15 19:33:20 ccarouge Exp $
      MODULE SOAPROD_MOD
!
!******************************************************************************
!  Module SOAPROD_MOD contains variables and routines which are used to read
!  and write APROD and GPROD restart files for SOA simulations. (dkh, 11/09/06)  
!
!  Module Variables:
!  ============================================================================
!  (1 ) INPUT_SOAPROD_FILE   : Full path name of the restart file to be read
!  (2 ) OUTPUT_SOAPROD_FILE  : Full path name (w/ tokens!) of output file
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_SOAPROD_FILE    : Writes restart file to disk 
!  (2 ) READ_SOAPROD_FILE    : Reads restart file from disk 
!  (4 ) CHECK_DIMENSIONS     : Ensures that restart file contains global data
!  (9 ) SET_SOAPROD          : Gets restart filenames from "input_mod.f"
!
!  GEOS-CHEM modules referenced by soaprod_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f   : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f   : Module containing NaN and other error check routines
!  (3 ) file_mod.f    : Module containing file unit numbers and error checks
!  (4 ) grid_mod.f    : Module containing horizontal grid information
!  (5 ) logical_mod.f : Module containing GEOS-CHEM logical switches
!  (6 ) time_mod.f    : Module containing routines for computing time & date
!  (7 ) tracer_mod.f  : Module containing GEOS-CHEM tracer array STT etc.
!  (8 ) carbon_mod.f  : Module containing routines for carbons simulation
!
!  NOTES:
!  (1 ) Based on restart_mod.f (dkh, 11/09/06) 
!  (2 ) Removed debug output (bmy, 4/21/10)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "restart_mod.f"
      !=================================================================

      ! Private module routines
      PRIVATE :: CHECK_DIMENSIONS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================    
      CHARACTER(LEN=255)   :: INPUT_SOAPROD_FILE  
      CHARACTER(LEN=255)   :: OUTPUT_SOAPROD_FILE 
      INTEGER, ALLOCATABLE :: JN(:),IN(:),NN(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_SOAPROD_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_SOAPROD_FILE creates GEOS-CHEM restart files of APROD
!  and GPROD. (dkh, 11/09/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a restart file       
!
!  NOTES:
!  (1 ) Based on MAKE_RESTART_FILE 
!  (2 ) Removed debug output (bmy, 4/21/10)
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD
      USE CARBON_MOD,  ONLY : MHC, NPROD, NNOX, APROD, GPROD
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE FILE_MOD,    ONLY : IOERROR
      USE inquireMod,  ONLY : findFreeLUN
      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD, ONLY : LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE, GET_TAU

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! TRCOFFSET

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER              :: IU_OAP
      INTEGER              :: I,    I0, IOS, J,  J0, L, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      INTEGER              :: JHC, IPR, NOX 
      REAL*4               :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      ! fixed after 1st runs (hotp 2/27/09)
      !INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER              :: HALFPOLAR 
      INTEGER, PARAMETER   :: CENTER180 = 1
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! MAKE_SOAPROD_FILE begins here!
      !=================================================================

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Restart File: ' // 
     &           'APROD and GPROD'
      UNIT     = 'unitless'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the restart file for output -- binary punch format
      !=================================================================

      ! Find a free file LUN
      IU_OAP   = findFreeLUN()

      ! Copy the output restart file name into a local variable
      FILENAME = TRIM( OUTPUT_SOAPROD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      WRITE( 6, 100 ) TRIM( FILENAME ), IU_OAP
 100  FORMAT( '     - MAKE_SOAPROD_FILE: Writing ', a, ' on unit ', i4 )

      ! Open restart file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_OAP, FILENAME, TITLE )

      N = 0 

      !=================================================================
      ! Write APROD and GRPOD to the restart file
      !=================================================================
      ! APROD
      CATEGORY = 'SOAAPROD'
      DO JHC = 1, MHC
      DO IPR = 1, NPROD(JHC) 
      DO NOX = 1, NNOX(JHC)
         
         N = N + 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = APROD(I,J,L,NOX,IPR,JHC)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_OAP,    MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
     &               J0+1,      1,         TRACER )
      ENDDO  
      ENDDO  
      ENDDO  

      ! Restart tracer # because we change category. (ccc, 3/17/10)
      N = 0


      ! GPROD
      CATEGORY = 'SOAGPROD'
      DO JHC = 1, MHC
      DO IPR = 1, NPROD(JHC) 
      DO NOX = 1, NNOX(JHC)
         
         N = N + 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = GPROD(I,J,L,NOX,IPR,JHC)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_OAP,    MODELNAME, LONRES,    LATRES,    
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
     &               J0+1,      1,         TRACER )
      ENDDO  
      ENDDO  
      ENDDO  
      ! Close file
      CLOSE( IU_OAP )
   
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SOAPROD_FILE: wrote file' )
      
      CALL CHECK_APRODGPROD()

      ! Return to calling program
      END SUBROUTINE MAKE_SOAPROD_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_SOAPROD_FILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_SOAPROD_FILE initializes GEOS-CHEM tracer concentrations 
!  from a restart file (binary punch file format) (bmy, 5/27/99, 7/20/04)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via F90 module dao_mod.f:
!  ============================================================================
!  (1 ) AD     : Air mass array (kg)       dim=(IIPAR,JJPAR,NNPAR) 
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!  (2 ) STT    : Tracer array (kg)         dim=(IIPAR,JJPAR,LLPAR,NNPAR)
!  (3 ) TCVV   : Air MW / Tracer MW array, dim=(NNPAR)
!
!  Passed via CMN_DIAG_mod:
!  ============================================================================
!  (1 ) TRCOFFSET : Offset for special tracers (like Rn, HCN, etc)
! 
!  NOTES:
!  (1 ) Now check that N = NTRACER - TRCOFFSET is valid.  
!        Also reorganize some print statements  (bmy, 10/25/99)
!  (2 ) Now pass LFORCE, LSPLIT via CMN_SETUP. (bmy, 11/4/99)
!  (3 ) Cosmetic changes, added comments (bmy, 3/17/00)
!  (4 ) Now use function NYMD_STRING from "time_mod.f" to generate a
!        Y2K compliant string for all data sets. (bmy, 6/22/00)
!  (5 ) Broke up sections of code into internal subroutines.  Also updated
!        comments & cleaned up a few things. (bmy, 7/17/00)
!  (6 ) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00)
!  (7 ) Print max & min of tracer regardless of the units (bmy, 10/5/00)
!  (8 ) Removed obsolete code from 10/00 (bmy, 12/21/00)
!  (9 ) Removed obsolete commented out code (bmy, 4/23/01)
!  (10) Added updates from amf for tagged Ox run.  Also updated comments
!        and made some cosmetic changes (bmy, 7/3/01)
!  (11) Bug fix: if starting from multiox restart file, then NTRACER 
!        will be greater than 40  but less than 60.  Adjust COPY_STT_FOR_OX
!        accordingly. (amf, bmy, 9/6/01)
!  (12) Now reference TRANUC from "charpak_mod.f" (bmy, 11/15/01)
!  (13) Updated comments (bmy, 1/25/02)
!  (14) Now reference AD from "dao_mod.f" (bmy, 9/18/02)
!  (15) Now added a call to DEBUG_MSG from "error_mod.f" (bmy, 2/11/03)
!  (16) Remove call to COPY_STT_FOR_OX, it's obsolete. (bmy, 8/18/03)
!  (17) Add fancy output string (bmy, 4/26/04)
!  (18) No longer use hardwired filename.  Also now reference "logical_mod.f"
!        and "tracer_mod.f" (bmy, 7/20/04)
!  (19) Removed debug output (bmy, 4/21/10)
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,   ONLY : OPEN_BPCH2_FOR_READ
      USE CARBON_MOD,  ONLY : APROD, GPROD
      USE ERROR_MOD,   ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,    ONLY : IOERROR
      USE inquireMod,  ONLY : findFreeLUN
      USE LOGICAL_MOD, ONLY : LSPLIT, LPRT
      USE TIME_MOD,    ONLY : EXPAND_DATE
      USE TRACER_MOD

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! TRCOFFSET

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: IU_OAP
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: JHC, IPR, NOX
      !INTEGER             :: NCOUNT(NNPAR) 
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_SOAPROD_FILE begins here!
      !=================================================================

      ! Initialize some variables
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open restart file and read top-of-file header
      !=================================================================
      
      ! Find a free file LUN
      IU_OAP   = findFreeLUN()

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_SOAPROD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'S O A P R O D   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME ), IU_OAP
 100  FORMAT( 'READ_SOAPROD_FILE: Reading ', a, ' on unit ', i4 )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_OAP, FILENAME )
      
      ! Echo more output
      WRITE( 6, 110 )
 110  FORMAT( /, 'Min and Max of APROD, GRPOD, as read from the file:',
     &        /, '(unitless)' )
      
      !=================================================================
      ! Read ratios -- store directly in APROD, GPROD
      !=================================================================
      DO 
         READ( IU_OAP, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_OAP,'read_soaprod_file:4' )

         READ( IU_OAP, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_OAP,'read_soaprod_file:5')

         READ( IU_OAP, IOSTAT=IOS ) 
     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_OAP,'read_soaprod_file:6')

         !==============================================================
         ! Assign data from the TRACER array to either APROD or GPROD
         !==============================================================
  
         ! APROD
         IF ( CATEGORY(1:8) == 'SOAAPROD' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

            ! Convert TRACER from its native units to [v/v] mixing ratio
            !CALL CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )

            JHC = JN(NTRACER)
            IPR = IN(NTRACER)
            NOX = NN(NTRACER)
 
            APROD(:,:,:,NOX,IPR,JHC) = TRACER(:,:,:)

         ! GPROD
         ELSEIF ( CATEGORY(1:8) == 'SOAGPROD' ) THEN 

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

            ! Convert TRACER from its native units to [v/v] mixing ratio
            !CALL CONVERT_TRACER_TO_VV( NTRACER, TRACER, UNIT )

            JHC = JN(NTRACER)
            IPR = IN(NTRACER)
            NOX = NN(NTRACER)
 
            GPROD(:,:,:,NOX,IPR,JHC) = TRACER(:,:,:)

         ELSE
            CALL ERROR_STOP('houston, we have a problem', 'soaprod_mod')
         ENDIF

      ENDDO

      ! Close file
      CLOSE( IU_OAP )      

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_SOAPROD_FILE: read file' )

      CALL CHECK_APRODGPROD()

      ! Return to calling program
      END SUBROUTINE READ_SOAPROD_FILE

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIMENSIONS( NI, NJ, NL ) 
!
!******************************************************************************
!  Subroutine CHECK_DIMENSIONS makes sure that the dimensions of the
!  restart file extend to cover the entire grid. (bmy, 6/25/02, 10/15/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NI (INTEGER) : Number of longitudes read from restart file
!  (2 ) NJ (INTEGER) : Number of latitudes  read from restart file
!  (3 ) NL (INTEGER) : Numbef of levels     read from restart file
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
      INTEGER, INTENT(IN) :: NI, NJ, NL

      !=================================================================
      ! CHECK_DIMENSIONS begins here!
      !=================================================================

      ! Error check longitude dimension: NI must equal IIPAR
      IF ( NI /= IIPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in soaprod file!'
         WRITE( 6, '(a)' ) 'Wrong number of longitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Error check latitude dimension: NJ must equal JJPAR
      IF ( NJ /= JJPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in soaprod file!'
         WRITE( 6, '(a)' ) 'Wrong number of latitudes encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
      
      ! Error check vertical dimension: NL must equal LLPAR
      IF ( NL /= LLPAR ) THEN
         WRITE( 6, '(a)' ) 'ERROR reading in soaprod file!'
         WRITE( 6, '(a)' ) 'Wrong number of levels encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS (restart_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS

!------------------------------------------------------------------------------
      SUBROUTINE CHECK_APRODGPROD()

!******************************************************************************
!     CHECK_APRODGPROD determines if there are bad APROD, GPROD values
!     no arguments
!
!******************************************************************************

      ! References to F90 modules
      USE CARBON_MOD, ONLY : APROD, GPROD, NNOX,
     &                       PRODPERCOFSTT, NPROD

      USE CMN_SIZE_MOD

      ! Local variables
      REAL*8           :: ATEMP, GTEMP
      INTEGER          :: IPR, NOX, JHC
      INTEGER          :: I, J, L

      !=================================================================
      ! CHECK_APRODGPROD begins here!
      !=================================================================
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR


      ! error checks
      ATEMP = 0.d0
      GTEMP = 0.d0
      DO JHC = 1, 3
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)
         ATEMP = ATEMP + APROD(I,J,L,NOX,IPR,JHC)
         GTEMP = GTEMP + GPROD(I,J,L,NOX,IPR,JHC)
      ENDDO
      ENDDO
      ENDDO
      ATEMP = (ATEMP - 1d0)**2.0d0
      GTEMP = (GTEMP - 1d0)**2.0d0
      ! update if more than 5% off
      IF( ATEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD APROD'
         print*, ATEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           APROD(I,J,L,NOX,IPR,1:3) = PRODPERCOFSTT(1:3)
         ENDDO
         ENDDO
      ENDIF
      IF( GTEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD GPROD'
         print*, GTEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           GPROD(I,J,L,NOX,IPR,1:3) = PRODPERCOFSTT(1:3)
         ENDDO
         ENDDO
      ENDIF

      ! error checks
      ATEMP = 0.d0
      GTEMP = 0.d0
      JHC = 4
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)
         ATEMP = ATEMP + APROD(I,J,L,NOX,IPR,JHC)
         GTEMP = GTEMP + GPROD(I,J,L,NOX,IPR,JHC)
      ENDDO
      ENDDO
      ATEMP = (ATEMP - 1d0)**2.0d0
      GTEMP = (GTEMP - 1d0)**2.0d0
      ! update if more than 5% off
      IF( ATEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD APROD'
         print*,ATEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           APROD(I,J,L,NOX,IPR,4) = PRODPERCOFSTT(4)
         ENDDO
         ENDDO
      ENDIF
      IF( GTEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD GPROD'
         print*,GTEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           GPROD(I,J,L,NOX,IPR,4) = PRODPERCOFSTT(4)
         ENDDO
         ENDDO
      ENDIF 


      ! error checks
      ATEMP = 0.d0
      GTEMP = 0.d0
      JHC = 5
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)
         ATEMP = ATEMP + APROD(I,J,L,NOX,IPR,JHC)
         GTEMP = GTEMP + GPROD(I,J,L,NOX,IPR,JHC)
      ENDDO
      ENDDO
      ATEMP = (ATEMP - 1d0)**2.0d0
      GTEMP = (GTEMP - 1d0)**2.0d0
      ! update if more than 5% off
      IF( ATEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD APROD'
         print*,ATEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           APROD(I,J,L,NOX,IPR,5) = PRODPERCOFSTT(5)
         ENDDO
         ENDDO
      ENDIF
      IF( GTEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD GPROD'
         print*,GTEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           GPROD(I,J,L,NOX,IPR,5) = PRODPERCOFSTT(5)
         ENDDO
         ENDDO
      ENDIF 

      ! error checks
      ATEMP = 0.d0
      GTEMP = 0.d0
      JHC = 6
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)
         ATEMP = ATEMP + APROD(I,J,L,NOX,IPR,JHC)
         GTEMP = GTEMP + GPROD(I,J,L,NOX,IPR,JHC)
      ENDDO
      ENDDO
      ATEMP = (ATEMP - 1d0)**2.0d0
      GTEMP = (GTEMP - 1d0)**2.0d0
      ! update if more than 5% off
      IF( ATEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD APROD'
         print*,ATEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           APROD(I,J,L,NOX,IPR,6) = PRODPERCOFSTT(6)
         ENDDO
         ENDDO
      ENDIF
      IF( GTEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD GPROD'
         print*,GTEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
           GPROD(I,J,L,NOX,IPR,6) = PRODPERCOFSTT(6)
         ENDDO
         ENDDO
      ENDIF 


      ! error checks
      ATEMP = 0.d0
      GTEMP = 0.d0
      DO JHC = 7, 9
      DO NOX = 1, NNOX(JHC)
      DO IPR = 1, NPROD(JHC)
         ATEMP = ATEMP + APROD(I,J,L,NOX,IPR,JHC)
         GTEMP = GTEMP + GPROD(I,J,L,NOX,IPR,JHC)
      ENDDO
      ENDDO
      ENDDO
      ATEMP = (ATEMP - 1d0)**2.0d0
      GTEMP = (GTEMP - 1d0)**2.0d0
      ! update if more than 5% off
      IF( ATEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD APROD'
         print*,ATEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
            APROD(I,J,L,NOX,IPR,7:9) = PRODPERCOFSTT(7:9)
         ENDDO
         ENDDO
      ENDIF
      IF( GTEMP>=1.0d-12 ) THEN
         WRITE( 6, '(a)' ) 'BAD GPROD'
         print*,GTEMP
         DO IPR = 1, NPROD(JHC)
         DO NOX = 1, NNOX(JHC)
            GPROD(I,J,L,NOX,IPR,7:9) = PRODPERCOFSTT(7:9)
         ENDDO
         ENDDO
      ENDIF 

      ENDDO                     ! loops over boxes
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CHECK_APRODGPROD

!------------------------------------------------------------------------------


      SUBROUTINE FIRST_APRODGPROD()
!******************************************************************************
!     Use this subroutine if you need to create values for APROD and
!     GPROD
!     for purposes of creating an APROD,GPROD restart file from scratch
!     (note that zeros may never be removed so we need this)
!
!******************************************************************************

      ! References to F90 modules
      USE CARBON_MOD, ONLY : APROD, GPROD, NNOX, NPROD,
     &                       PRODPERCOFSTT, MHC

      USE CMN_SIZE_MOD

      ! Local variables
      INTEGER          :: IPR, NOX, JHC
      INTEGER          :: I, J, L

      !=================================================================
      ! FIRST_APRODGPROD begins here!
      !=================================================================
      DO I = 1, IIPAR
      DO J = 1, JJPAR
      DO L = 1, LLPAR

      DO JHC = 1, 9
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)
        APROD(I,J,L,NOX,IPR,JHC) = PRODPERCOFSTT(JHC)
        GPROD(I,J,L,NOX,IPR,JHC) = PRODPERCOFSTT(JHC)
      ENDDO
      ENDDO
      ENDDO

      ENDDO
      ENDDO
      ENDDO


      END SUBROUTINE FIRST_APRODGPROD
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

      SUBROUTINE SET_SOAPROD( )
!
!******************************************************************************
!  Subroutine SET_SOAPROD initializes the variables INPUT_SOAPROD_FILE and
!  OUTPUT_SOAPROD_FILE with hardwired file  names and allocates 
!  ratios to be read in from restart file for SOA simulation. (dkh, 11/09/06)  
!
!  NOTES:
!******************************************************************************
!
      ! Reference to f90 modules
      USE CARBON_MOD, ONLY : MNOX, MPROD, MHC, APROD, GPROD, NNOX,
     &                       PRODPERCOFSTT, NPROD
      USE ERROR_MOD,    ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD   ! Size params

      ! Arguments
      CHARACTER(LEN=255) :: INFILE, OUTFILE
      
      ! Local variables
      INTEGER :: STAT, AS, MAXPROD, N, NOX, JHC, IPR

      !=================================================================
      ! SET_SOAPROD begins here
      !=================================================================
      INPUT_SOAPROD_FILE  = 'soaprod.YYYYMMDDhh' 
      OUTPUT_SOAPROD_FILE = 'soaprod.YYYYMMDDhh'


      ! Calculate max number of SOA products
!      MAXPROD = MHC * MPROD * MNOX * 2 
      MAXPROD = MHC * MPROD * MNOX
   
      ALLOCATE( GPROD( IIPAR, JJPAR, LLPAR, MNOX, MPROD, MHC ),
     &          STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GPROD' )
      GPROD = 0D0
      
      ALLOCATE( APROD( IIPAR, JJPAR, LLPAR, MNOX, MPROD, MHC ),
     &         STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'APROD' )
      APROD = 0D0
    
      ALLOCATE( JN( MAXPROD ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JN' )
      JN = 0

      ALLOCATE( IN( MAXPROD ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IN' )
      IN = 0

      ALLOCATE( NN( MAXPROD ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NN' )
      NN = 0

      ! Set size of NOX dimension for HCs (dkh, 11/05/06)  
      NNOX(1:6) = 1  ! Biogenic terps
      NNOX(7:9) = 2  ! aromatics

      ! Set size of IPR dimension for HCs 
      NPROD(1:4) = 3  
      NPROD(5)   = 2  
      NPROD(6)   = 2  
      NPROD(7:9) = 2  

    
      ! The number of products (GM0 or AM0) lumped into a
      ! particular tracer is calculated as:
      ! N =  ( sum of NOX groups ) * ( # of products per hydrocarbons )
      ! IDTSOA1 = HC 1-3
      ! IDTSOA2 = HC 4
      ! IDTSOA3 = HC 5
      ! IDTSOA4 = HC 6
      ! IDTSOA5 = HC 7-9
      ! likewise for the SOGs
      PRODPERCOFSTT      = 0.d0

      PRODPERCOFSTT(1:3) = 1.d0 / ( SUM(NNOX(1:3)) * 3.d0 )
      PRODPERCOFSTT(4)   = 1.d0 / (     NNOX(4)    * 3.d0 )
      PRODPERCOFSTT(5)   = 1.d0 / (     NNOX(5)    * 2.d0 )
      PRODPERCOFSTT(6)   = 1.d0 / (     NNOX(6)    * 2.d0 )

      PRODPERCOFSTT(7:9) = 1.d0 / ( SUM(NNOX(7:9)) * 2.d0 )

      ! Initialize mapping arrays for reading in products from file
      N = 0 
      ! for APROD
      DO JHC = 1, MHC
      DO IPR = 1, NPROD(JHC)
      DO NOX = 1, NNOX(JHC)

         N = N + 1
         JN(N) = JHC
         IN(N) = IPR
         NN(N) = NOX

      ENDDO
      ENDDO
      ENDDO
!      ! for GPROD
!      DO JHC = 1, MHC
!      DO IPR = 1, NPROD(JHC)
!      DO NOX = 1, NNOX(JHC)
!
!         N = N + 1
!         JN(N) = JHC
!         IN(N) = IPR
!         NN(N) = NOX
!
!      ENDDO
!      ENDDO
!      ENDDO
      
      ! Return to calling program
      END SUBROUTINE SET_SOAPROD

!------------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_SOAPROD
!
!******************************************************************************
!  Subroutine CLEANUP_SOAPROD deallocates module variables. (dkh, 11/09/06)  
!
!  NOTES:
!******************************************************************************
!     
      USE CARBON_MOD, ONLY : APROD, GPROD 

      !=================================================================
      ! CLEANUP_SOAPROD begins here!
      !=================================================================
      IF ( ALLOCATED( JN    ) ) DEALLOCATE( JN )
      IF ( ALLOCATED( IN    ) ) DEALLOCATE( IN )
      IF ( ALLOCATED( NN    ) ) DEALLOCATE( NN )
      IF ( ALLOCATED( APROD ) ) DEALLOCATE( APROD )
      IF ( ALLOCATED( GPROD ) ) DEALLOCATE( GPROD )

      ! Return to calling program 
      END SUBROUTINE CLEANUP_SOAPROD

!------------------------------------------------------------------------------

      ! End of module
      END MODULE SOAPROD_MOD
