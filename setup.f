! $Id: setup.f,v 1.9 2004/07/15 18:17:47 bmy Exp $
      SUBROUTINE SETUP( NYMDb, NYMDe,     NHMSb,   NHMSe,   NDT,      
     &                  NTDT,  NDIAGTIME, ALPHA_d, ALPHA_n, IORD,     
     &                  JORD,  KORD,      J1,      Umax,    IGD,      
     &                  KS,    PSTP )
!
!******************************************************************************
!  Subroutine SETUP reads in parameters specific to the GEOS-CHEM model.
!  (bmy, bey, bdf, 5/27/99, 7/13/04)
!                                                                           
!  NOTES:
!  (1 ) SETUP is written in Fixed-Form Fortran 90.
!  (2 ) Most of the arguments to SETUP are passed through the argument list.
!        Logical flags and file I/O units are passed in CMN_SETUP, since it
!        is easier to incorporate changes that way.
!  (3 ) Subroutines READ_PHIS, READ_I6, READ_A6, READ_A3, and READ_KZZ
!        are used to read in the GEOS-1 or GEOS-STRAT met fields.  
!  (4 ) All DAO met fields are now opened and closed in subroutine
!        OPEN_MET_FIELDS.  Subroutine CLOSE_FILES closes all open files
!        (including met fields files) at the end of the run.
!  (5 ) Use IMPLICIT NONE to explicitly declare all variables.                
!  (6 ) Read input information from unit 99, 'input.geos'.
!  (7 ) Add variables for reading KZZ fields (bdf, bmy, 3/18/99)
!  (8 ) Use subroutine IOERROR to trap I/O error conditions (bmy, 5/27/99)
!  (9 ) Bug fix -- now use proper 26-layer GEOS-STRAT sigma centers
!        and edges (bmy, 8/11/99)
!  (10) Add internal subroutines CHECK_STRING and CHECK_DIRECTORY for
!        error checking.  Also use F90 declaration syntax. (bmy, 2/1/00)
!  (11) Now flush the std output buffer before halting with error
!        conditions. (bmy, 2/7/00)
!  (12) Now only check either the GEOS-1 or GEOS-STRAT met fields directory
!        depending on the settings in "define.h" (bmy, 4/17/00)
!  (13) If LMFCT = T, set LFILL = F for TPCORE (bmy, 4/24/00)
!  (14) Make NYMDb, NYMDe Y2K compliant (i.e. 8 digits -- YYYYMMDD) for use
!        with new routines from "time_mod.f".  Also set LUNZIP = T, since
!        we need to keep the met field files zipped. (bmy, 6/27/00)
!  (15) Added additional internal routines for extracting character
!        data from the "input.geos" file, and for setting sigma level
!        variables.  Also check IJLN against NTRACE.   Reference
!        function NYMD6_2_NYMD8 from "time_mod.f".  (bmy, 7/17/00)
!  (16) Added internal routine CHECK_NTDT to make sure NTDT is not too 
!        small for the given grid (bmy, 9/8/00)
!  (17) Rename LGEOS1 to LBBSEA, for selecting seasonal or interannual
!        variability biomass burning emission (bmy, 9/11/00) 
!  (18) Now use IOS /= 0 to trap both I/O errors and EOF. (bmy, 9/13/00) 
!  (19) Rename LRWND to LTOMSAI, for computing interannual variability in
!        biomass burning emissions from TOMS Aerosol Index data (bmy, 9/28/00) 
!  (20) Echo NYMD1, NHMS1, NYMD2, NHMS2 using I6.6 format (bmy, 10/11/00)  
!  (21) Rename LUNZIP to LSTDRUN.  Updated comments. (bmy, 12/6/00)
!  (22) Removed obsolete code from 12/6/00 (bmy, 12/21/00)
!  (23) Added internal subroutine CHECK_NYMDb_NYMDe to make sure that the
!        end date (NYMDe) is not less than the start date (NYMDb).
!        (bmy, 3/15/01)
!  (24) Renamed LNOOPS to LWAIT -- this will cause GEOS-CHEM to wait 
!        until the met fields have been unzipped before proceeding.  This
!        is sometimes necessary for fast simulations. (bmy, 4/17/01)
!  (25) Removed obsolete commented out code (bmy, 4/23/01)
!  (26) Now test IJLN against LLPAR instead of LGLOB (bmy, 9/26/01)
!  (27) Now define sigma edges, centers for 30-level GEOS-3 grid (bmy, 9/27/01)
!  (28) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (29) For ALPHA, use ACCESS instead of INQUIRE.  Overwrite GRID_SUFFIX w/ 
!        the correct resolution string (bmy, 11/15/01) 
!  (30) Now don't read in the file unit #'s from "input.geos"; these are now
!        contained in "file_mod.f". Now reference IOERROR from "file_mod.f"
!        (bmy, 6/27/02)
!  (31) Removed all sigma center & edge variables from the arg list, since
!        we now use the hybrid grid formulation in "pressure_mod.f".  Also
!        eliminated obsolete, commented-out code. (dsa, bdf, bmy, 8/22/02)
!  (32) Now reference ERROR_STOP and GEOS_CHEM_STOP from "error_mod.f" 
!        Eliminated reference to obsolete IJLN.  Also renamed LFORCE to
!        LSULF for coupled fullchem/aerosol run. Bug fix: now stop if NTDT is 
!        greater than 1800s (4x5), 900s (2x2.5), or 600s (4x5).(bmy, 11/22/02)
!  (33) Now use external ACCESS function for Alpha compiler (bmy, 3/23/03)
!  (34) Hardwire LRERD = .FALSE. -- we now have many years of met fields that
!        we can use for testing.  Renamed LKZZ to LTPFV to select GEOS-4/fvDAS
!        transport or the TPCORE 7.1m for GEOS-3 winds (bmy, 4/28/03) 
!  (35) Removed GEOS_2_DIR.  Added GEOS_4_DIR (bmy, 6/18/03)
!  (36) Renamed LRERD to LUNZIP, this toggles met field unzipping on/off.
!        Disable checking of GEOS_1_DIR etc directory paths. (bmy, 12/9/03)
!  (37) Added flags for carbon & dust aerosols. Also do not read extra
!        lines from the top of the file. (rjp, tdf, bmy, 4/6/04) 
!  (38) Now reads LSHIPSO2 instead of LSPLIT.  LSPLIT is now automatically
!        set in "input.f" (bec, bmy, 5/20/04)
!  (39) Renamed LATEQ to LSOA (bmy, 7/13/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_RES_EXT
      USE FILE_MOD,  ONLY : IOERROR
      USE ERROR_MOD, ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE TIME_MOD,  ONLY : NYMD6_2_NYMD8, EXPAND_DATE

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_GCTM"   ! Physical constants
#     include "CMN_SETUP"  ! Logical switches other variables

      ! Arguments 
      INTEGER, INTENT(OUT) :: NYMDb, NHMSb, NYMDe, NHMSe
      INTEGER, INTENT(OUT) :: NDT,   NTDT,  NDIAGTIME  
      INTEGER, INTENT(OUT) :: IORD,  JORD,  KORD 
      INTEGER, INTENT(OUT) :: J1,    IGD,   KS  

      REAL*8,  INTENT(OUT) :: ALPHA_d, ALPHA_n, Umax, PSTP

      ! Local variables
      CHARACTER(LEN=79)    :: TITLE 
      CHARACTER(LEN=46)    :: STR1, STR2, STR3, STR4, STR5
      CHARACTER(LEN=30)    :: TMPSTR

      INTEGER              :: I, I1, I2, IS, J, JS, K ! Loop and 
      INTEGER              :: L, L1, L2, LS, N, NS    !  counter variables
      INTEGER              :: IDUM1, IOS 
!
!******************************************************************************
!  SETUP begins here!                                                       
!
!  Open the 'input.geos' file as unit 99 and read timing parameters
!  Echo NYMD1, NHMS1, NYMD2, NHMS2 using I6.6 format (bmy, 10/11/00)
!******************************************************************************
!
      OPEN ( 99, FILE='input.geos', STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:0' )

      READ ( 99, 10, IOSTAT=IOS ) TITLE  
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:1' )

      READ ( 99, 15, IOSTAT=IOS ) NYMDb, NHMSb, NYMDe,     NHMSe, STR1
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:2' )

      READ ( 99, 16, IOSTAT=IOS ) NDT,   NTDT,  NDIAGTIME,        STR2
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:3' )

      WRITE ( 6, 10 ) TITLE
      WRITE ( 6, 17 ) NYMDb, NHMSb, NYMDe,     NHMSe, STR1
      WRITE ( 6, 16 ) NDT,   NTDT,  NDIAGTIME,        STR2

      ! Error check NYMDb, NYMDe
      CALL CHECK_NYMDb_NYMDe

      ! Error check NTDT
      CALL CHECK_NTDT

      ! Change NYMDb, NYMDe to YYYYMMDD format
      NYMDb = NYMD6_2_NYMD8( NYMDb )
      NYMDe = NYMD6_2_NYMD8( NYMDe )
!
!******************************************************************************
!  Read in GEOS-CTM logical flags from INPUT.GEOS                            
!******************************************************************************
!
      READ ( 99, 10, IOSTAT=IOS ) TITLE  
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:8' )

      READ ( 99, 25, IOSTAT=IOS ) 
     &     LEMIS,  LDRYD,  LCHEM,   LTRAN,  LTPFV,  STR1
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:9' )
 
      READ ( 99, 25, IOSTAT=IOS ) 
     &     LTURB,  LCONV,  LWETD,   LDBUG,  LMONOT,  STR2 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:10' )

      READ ( 99, 25, IOSTAT=IOS )
     &     LWAIT,  LBBSEA, LUNZIP,   LSVGLB, LTOMSAI, STR3 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:11' )

      READ ( 99, 25, IOSTAT=IOS ) 
     &     LMFCT,  LFILL,  LSTDRUN, LDEAD,  LSHIPSO2, STR4   
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:12' )

      ! Added for various aerosol chemistry (rjp, tdf, bmy, 4/2/04)
      READ ( 99, 25, IOSTAT=IOS ) 
!----------------------------------------------------------------------
! Prior to 7/13/04:
! Renamed LATEQ to LSOA (bmy, 7/13/04)
!     &     LSULF,  LCARB,  LDUST,   LSSALT,  LATEQ,  STR5   
!----------------------------------------------------------------------
     &     LSULF,  LCARB,  LDUST,   LSSALT,  LSOA,    STR5   
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:12a' )

      WRITE ( 6, 10 ) TITLE
      WRITE ( 6, 25 ) LEMIS,  LDRYD,  LCHEM,   LTRAN,  LTPFV,   STR1 
      WRITE ( 6, 25 ) LTURB,  LCONV,  LWETD,   LDBUG,  LMONOT,  STR2 
      WRITE ( 6, 25 ) LWAIT,  LBBSEA, LUNZIP,  LSVGLB, LTOMSAI, STR3 
      WRITE ( 6, 25 ) LMFCT,  LFILL,  LSTDRUN, LDEAD,  LSHIPSO2,STR4
!-----------------------------------------------------------------------
! Prior to 7/13/04:
! Renamed LATEQ to LSOA (bmy, 7/13/04)
!      WRITE ( 6, 25 ) LSULF,  LCARB,  LDUST,   LSSALT, LATEQ,   STR5   
!-----------------------------------------------------------------------
      WRITE ( 6, 25 ) LSULF,  LCARB,  LDUST,   LSSALT, LSOA,    STR5   
!
!******************************************************************************
!  Read in other GEOS-CTM parameters                                        
!******************************************************************************
!
      READ ( 99, 10, IOSTAT=IOS ) TITLE 
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:13' )

      READ ( 99, 30, IOSTAT=IOS ) ALPHA_d, ALPHA_n,              STR1
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:14') 

      READ ( 99, 30, IOSTAT=IOS ) Umax,    PSTP,                 STR2
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:15' )

      READ ( 99, 35, IOSTAT=IOS ) IORD,    JORD,    KORD, IDUM1, STR3
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:16') 

      READ ( 99, 35, IOSTAT=IOS ) J1,      IGD,     KS,   IDUM1, STR4
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:17' )

      WRITE ( 6, 10 ) TITLE   
      WRITE ( 6, 30 ) ALPHA_d, ALPHA_n,              STR1
      WRITE ( 6, 30 ) Umax,    PSTP,                 STR2
      WRITE ( 6, 35 ) IORD,    JORD,    KORD, IDUM1, STR3
      WRITE ( 6, 35 ) J1,      IGD,     KS,   IDUM1, STR4
!
!******************************************************************************
!  Read in time limits for NO, J-Value and 2-PM diagnostics
!******************************************************************************
!
      READ ( 99, 10, IOSTAT=IOS ) TITLE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:21' )

      READ ( 99, 30, IOSTAT=IOS ) HR1_NO,  HR2_NO,  STR1
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:22' )

      READ ( 99, 30, IOSTAT=IOS ) HR1_JV,  HR2_JV,  STR2
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:23' )

      READ ( 99, 30, IOSTAT=IOS ) HR1_OH,  HR2_OH,  STR3
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:24' )

      READ ( 99, 30, IOSTAT=IOS ) HR1_OTH, HR2_OTH, STR4
      IF ( IOS /= 0 ) CALL IOERROR( IOS, 99, 'setup:25' )

      WRITE ( 6, 10 ) TITLE
      WRITE ( 6, 30 ) HR1_NO,  HR2_NO,  STR1
      WRITE ( 6, 30 ) HR1_JV,  HR2_JV,  STR2
      WRITE ( 6, 30 ) HR1_OH,  HR2_OH,  STR3
      WRITE ( 6, 30 ) HR1_OTH, HR2_OTH, STR4
! 
!******************************************************************************
!  Read OPERATING SYSTEM COMMAND STRINGS from 'input.geos'
!
!  Halt execution if "COMMAND STRINGS" is not found the input line, since 
!  this will likely denote a mis-indexed 'input.geos' file (bmy, 5/25/99)
!******************************************************************************
!
      ! Read separator line
      TITLE = READ_LINE()

      ! Error check
      IF ( TITLE(1:15) /= 'COMMAND STRINGS' ) THEN
         WRITE( 6, 100   )
         WRITE( 6, '(a)' ) '"COMMAND STRINGS" header not found!' 
         WRITE( 6, '(a)' ) 'The "input.geos" file is mis-indexed!' 
         WRITE( 6, '(a)' ) 'STOP in setup.f' 
         WRITE( 6, 100   )
         CALL GEOS_CHEM_STOP
      ENDIF
        
      ! Assign string values from data in the "input.geos" file
      BACKGROUND = ASSIGN_STRING( 'BACKGROUND'            )
      REDIRECT   = ASSIGN_STRING( 'REDIRECT'              )
      REMOVE_CMD = ASSIGN_STRING( 'REMOVE_CMD'            )
      SEPARATOR  = ASSIGN_STRING( 'SEPARATOR'             )
      SPACE      = ASSIGN_STRING( 'SPACE', NOCHECK=.TRUE. )
      UNZIP_CMD  = ASSIGN_STRING( 'UNZIP_CMD'             )
      WILD_CARD  = ASSIGN_STRING( 'WILD_CARD'             )
! 
!******************************************************************************
!  Read the FILE SUFFIXES from the 'input.geos' file
!******************************************************************************
!
      ! Read separator line
      TITLE = READ_LINE()

      ! Error check
      IF ( TITLE(1:13) /= 'FILE SUFFIXES' ) THEN
         write( 6, 100   )
         WRITE( 6, '(a)' ) '"FILE_SUFFIXES" header not found!' 
         WRITE( 6, '(a)' ) 'The "input.geos" file is mis-indexed!' 
         WRITE( 6, '(a)' ) 'STOP in setup.f' 
         WRITE( 6, 100   )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Assign string values from data in the "input.geos" file
      A3_SUFFIX   = ASSIGN_STRING( 'A3_SUFFIX'   )
      A6_SUFFIX   = ASSIGN_STRING( 'A6_SUFFIX'   )
      I6_SUFFIX   = ASSIGN_STRING( 'I6_SUFFIX'   )
      PH_SUFFIX   = ASSIGN_STRING( 'PH_SUFFIX'   )
      KZZ_SUFFIX  = ASSIGN_STRING( 'KZZ_SUFFIX'  )
      GRID_SUFFIX = ASSIGN_STRING( 'GRID_SUFFIX' )
      ZIP_SUFFIX  = ASSIGN_STRING( 'ZIP_SUFFIX'  )

      ! Overwrite GRID_SUFFIX w/ the correct res string (bmy, 11/15/01)
      GRID_SUFFIX = '.' // GET_RES_EXT()
! 
!******************************************************************************
!  Read the DATA DIRECTORIES from the 'input.geos' file
!
!  Halt execution if "COMMAND STRINGS" is not found the input line, since 
!  this will likely denote a mis-indexed 'input.geos' file (bmy, 5/25/99)
!******************************************************************************
!
      ! Read separator line
      TITLE = READ_LINE()

      ! Error check
      IF ( TITLE(1:16) /= 'DATA DIRECTORIES' ) THEN
         WRITE( 6, 100   )
         WRITE( 6, '(a)' ) '"DATA DIRECTORIES" header not found!' 
         WRITE( 6, '(a)' ) 'The "input.geos" file is mis-indexed!' 
         WRITE( 6, '(a)' ) 'STOP in setup.f' 
         WRITE( 6, 100   )
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Assign string values from data in the "input.geos" file
      DATA_DIR   = ASSIGN_STRING( 'DATA_DIR'   )
      GEOS_1_DIR = ASSIGN_STRING( 'GEOS_1_DIR' )
      GEOS_S_DIR = ASSIGN_STRING( 'GEOS_S_DIR' )
      GEOS_3_DIR = ASSIGN_STRING( 'GEOS_3_DIR' )
      GEOS_4_DIR = ASSIGN_STRING( 'GEOS_4_DIR' )
      TEMP_DIR   = ASSIGN_STRING( 'TEMP_DIR'   )

      ! Read last line of file
      TITLE = READ_LINE()

      ! Close the "input.geos" file -- we're done
      CLOSE ( 99 )
! 
!******************************************************************************
!  Call the internal subroutine CHECK_DIRECTORY to make sure that
!  each of the directories listed in "input.geos" are valid.
!******************************************************************************
!
      CALL CHECK_DIRECTORY( TRIM( DATA_DIR ) )

      ! Only check TEMP_DIR if we are unzipping met fields (bmy, 12/11/03)
      IF ( LUNZIP ) THEN
         CALL CHECK_DIRECTORY( TRIM( TEMP_DIR ) ) 
      ENDIF
!
!******************************************************************************
!  If LMFCT = T, then set LFILL = F (bmy, 4/24/00)
!  TPCORE does not need LFILL turned on when LMFCT is also turned on.
!******************************************************************************
!
      IF ( LMFCT ) THEN 
         LFILL = .false.
         WRITE( 6, '(a)' ) 'LMFCT = T; Setting LFILL = F for TPCORE!'
      ENDIF
! 
!******************************************************************************
!  Open the following GEOS-CTM file units:
!  =====================================================
!  (1) File unit #29  --> opened as 'debug' (reserved for debug!!!)
!******************************************************************************
!
      IF ( LDBUG ) THEN
         OPEN (29, FILE='debug', STATUS='unknown')
      ENDIF
!
!******************************************************************************
!  RETURN to calling program
!******************************************************************************
!
      RETURN
!
!******************************************************************************
!  FORMAT statements for reading in data from 'input.geos'
!******************************************************************************
!
 10   FORMAT( A79                  )
 15   FORMAT( 4I7,         5X, A46 )
 16   FORMAT( 3I7,        12X, A46 )
 17   FORMAT( 4(1x,I6.6),  5X, A46 )
 20   FORMAT( 5I4,        13X, A46 )
 21   FORMAT( 4I4,        17X, A46 )
 25   FORMAT( 5L4,        13X, A46 )
 26   FORMAT( 4L4,        17X, A46 )
 30   FORMAT( 2F8.3,      17X, A46 )
 35   FORMAT( 4I4,        17X, A46 )
 100  FORMAT( '=======================================',
     &        '=======================================' ) 
!
!******************************************************************************
!  Internal procedures -- Use the F90 CONTAINS command to inline 
!  subroutines that only can be called from subroutine SETUP.
!
!  All variables referenced in subroutine SETUP (local variables, F90 
!  module variables, or common block variables) also have scope within 
!  internal subroutines. 
!
!  List of Internal Procedures:
!  ============================================================================
!  (1) READ_LINE         - Reads one line at a time from the "input.geos" file
!  (2) CHECK_NYMDb_NYMDe - Makes sure that NYMDe is not less than NYMDb
!  (3) CHECK_NTDT        - Checks that NTDT is valid for the given grid
!  (4) CHECK_DIRECTORY  -- Ensures that a given directory exists
!  (5) CHECK_STRING     -- Ensures that a given string is not empty
!  (6) ASSIGN_STRING    -- Assigns data to character variables
!  (7) GET_SIGMA_LEVELS -- Initializes all of the sigma level variables 
!******************************************************************************
!
      CONTAINS

      FUNCTION READ_LINE() RESULT( LINE )

      !=================================================================
      ! Internal function READ_LINE reads one line from the
      ! "input.geos" file and prints it to the standard output.
      !=================================================================

      ! Local variable
      CHARACTER(LEN=79) :: LINE

      ! Read a one line from the "input.geos" file
      READ( 99, '(a79)', IOSTAT=IOS )  LINE
      IF ( IOS > 0 ) CALL IOERROR( IOS, 99, 'setup:read_line:1' )

      ! Echo to standard output
      WRITE( 6, '(a79)' )  LINE       

      ! Return to SETUP
      END FUNCTION READ_LINE

!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_NYMDb_NYMDe

      !=================================================================
      ! Internal Subroutine CHECK_NYMDb_NYMDe makes sure that NYMDe
      ! (the ending date) is not smaller than NYMDb (the start date)
      !=================================================================

      IF ( NYMDe < NYMDb ) THEN
         WRITE( 6, 100 )
         WRITE( 6, '(a   )' ) 'NYMDe cannot be less than NYMDb!'
         WRITE( 6, '(a,i8)' ) 'NYMDb = ', NYMDb 
         WRITE( 6, '(a,i8)' ) 'NYMDe = ', NYMDe 
         WRITE( 6, '(a,i8)' ) 'STOP in setup.f'
         WRITE( 6, 100 )
         CALL GEOS_CHEM_STOP
      ENDIF         

      ! Format statements
 100  FORMAT( '=======================================',
     &        '=======================================' ) 

      END SUBROUTINE CHECK_NYMDb_NYMDe

!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_NTDT

      !=================================================================
      ! Internal Subroutine CHECK_NTDT makes sure that the dynamic
      ! time step value NTDT is valid for the proper grid resolution
      !
      ! NOTES:
      ! (1) Bug fix: now stop if NTDT is greater than 1800s (4x5),
      !      900s (2x2.5), or 600s (4x5).  (bmy, 11/22/02)
      !=================================================================

#if   defined( GRID4x5 ) 
      IF ( NTDT > 1800 ) THEN
         CALL ERROR_STOP( 'NTDT > 1800s for 4x5 grid!', 'setup.f' )
      ENDIF
      
#elif defined( GRID2x25 )
      IF ( NTDT > 900 ) THEN
         CALL ERROR_STOP( 'NTDT > 900s for 2 x 2.5 grid!', 'setup.f' )
      ENDIF

#elif defined( GRID1x1 )
      IF ( NTDT > 600 ) THEN
         CALL ERROR_STOP( 'NTDT > 600s for 4x5 grid!', 'setup.f' )
      ENDIF

#endif

      ! Return to SETUP
      END SUBROUTINE CHECK_NTDT

!-----------------------------------------------------------------------------

      SUBROUTINE CHECK_DIRECTORY( DIR )
      
      !=================================================================
      ! Internal Subroutine CHECK_DIRECTORY makes sure that the 
      ! given directory is valid.
      !
      ! NOTES:
      ! (1 ) Now use external ACCESS function for Alpha compiler 
      !       (lyj, bmy, 3/20/03)
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: DIR

      ! Local variables
      LOGICAL                      :: IT_EXISTS
      INTEGER                      :: I

#if   defined( COMPAQ )
      
      ! For COMPAQ/Alpha use external ACCESS function
      INTEGER*4, EXTERNAL          :: ACCESS

      ! Test whether directory exists (lyj, bmy, 3/23/03)
      IT_EXISTS = ( ACCESS( TRIM( DIR ), ' ' ) == 0 )

#else

      ! Test whether directory exists w/ F90 INQUIRE function
      INQUIRE( FILE=TRIM( DIR ), EXIST=IT_EXISTS )

#endif     
 
      ! If the directory doesn't exist, write an error
      ! message, flush the std output buffer, and stop
      IF ( .not. IT_EXISTS ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) TRIM( DIR ) // ' does not exist!'
         WRITE( 6, '(a)' ) 'STOP in setup.f'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL FLUSH( 6 )
         CALL GEOS_CHEM_STOP
      ENDIF
         
      ! Also make sure that the last character of DIR is
      ! the separator character ('/' for Unix environment)
      I = LEN_TRIM( DIR )

      ! If the last character is not a separator, write an error
      ! message, flush the std output buffer, and stop
      IF ( DIR(I:I) /= TRIM( SEPARATOR ) ) THEN 
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, 110   ) TRIM( DIR ), TRIM( SEPARATOR )
         WRITE( 6, '(a)' ) 'STOP in setup.f'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL FLUSH( 6 )
         CALL GEOS_CHEM_STOP
      ENDIF

 110  FORMAT( a, ' must end with the "', a, '" character!' )

      ! Return to SETUP
      END SUBROUTINE CHECK_DIRECTORY

!------------------------------------------------------------------------------
      
      SUBROUTINE CHECK_STRING( STR, NAME )

      !=================================================================
      ! Internal subroutine CHECK_STRING makes sure that the 
      ! given string is not empty
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: STR, NAME

      ! Use the F90 LEN_TRIM command to make sure that the string has 
      ! more than zero characters.  Otherwise write an error message, 
      ! flush the std output buffer, and then stop.
      IF ( LEN_TRIM( STR ) == 0 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) TRIM( NAME ) // ' does not exist!'
         WRITE( 6, '(a)' ) 'STOP in setup.f'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL FLUSH( 6 ) 
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Return to SETUP
      END SUBROUTINE CHECK_STRING  

!------------------------------------------------------------------------------

      FUNCTION ASSIGN_STRING( NAME, NOCHECK ) RESULT( STR )

      !=================================================================
      ! Internal function ASSIGN_STRING reads a line from the
      ! "input.geos" file and extracts a 30-character string.
      !
      ! NAME is a character variable for the name of the string, 
      ! which is only used to display an error message.
      !=================================================================

      ! Arguments
      CHARACTER(LEN=*),  INTENT(IN)  :: NAME
      LOGICAL, OPTIONAL, INTENT(IN)  :: NOCHECK
      
      ! Return value of the function
      CHARACTER(LEN=30)              :: STR
      CHARACTER(LEN=79)              :: LINE
      LOGICAL                        :: NOCHECK_TMP
      

      ! Save the value of NOCHECK in a temp variable
      IF ( PRESENT( NOCHECK ) ) THEN 
         NOCHECK_TMP = NOCHECK
      ELSE
         NOCHECK_TMP = .FALSE.
      ENDIF
          
      ! Read a line from the "input.geos" file
      ! Save the first 30 characters in STR
      LINE = READ_LINE()
      STR  = LINE( 1:30 )

      ! If necessary, make sure STR is not a null string
      IF ( .not. NOCHECK_TMP ) CALL CHECK_STRING( STR, NAME )

      ! Return to SETUP
      END FUNCTION ASSIGN_STRING

!------------------------------------------------------------------------------

      END SUBROUTINE SETUP










