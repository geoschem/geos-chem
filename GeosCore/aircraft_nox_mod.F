! $Id: aircraft_nox_mod.f,v 1.1 2009/09/16 14:06:43 bmy Exp $
      MODULE AIRCRAFT_NOX_MOD
!
!******************************************************************************
!  Module AIRCRAFT_NOX_MOD contains variables and routines for emission
!  of aircraft NOx fields into arrays for SMVGEAR. (bmy, 2/15/02, 10/3/07)
!
!  The aircraft NOx fields are stored on grid with 1-km vertical resolution.
!  These fields will be interpolated onto the GEOS-CHEM vertical grid.
!
!  Module Variables:
!  ============================================================================
!  (1 ) NAIR     (INTEGER) : Max number of layers of the 1-km native grid
!  (2 ) LAIREMS  (INTEGER) : Highest GEOS-CHEM level we will emit NOx into
!  (3 ) AIR      (REAL*8 ) : Array for NOx emissions on the 1-km grid
!  (4 ) AIREMIS  (REAL*8 ) : Array for NOx emissions on the GEOS-CHEM grid
!  (5 ) AIRPRESS (REAL*8 ) : Approx. pressure edges of the 1-km native grid
!  (6 ) EMIS_AC_NOx (REAL*8) : Array to pass aircraft NOx to SMVGEAR
!
!  Module Routines:
!  ============================================================================
!  (1 ) READAIR              : Routine to read NOx emissions from disk
!  (2 ) AIREMISS             : Routine to emit aircraft NOx into GEOS-CHEM
!  (3 ) INIT_AIRCRAFT_NOX    : Routine to allocate/initialize module variables
!  (4 ) CLEANUP_AIRCRAFT_NOX : Routine to deallocate module variables
! 
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) file_mod.f     : Module containing file unit numbers and error checks
!  (5 ) grid_mod.f     : Module containing horizontal grid information
!  (6 ) pressure_mod.f : Module containing routines to compute P(I,J,L) 
! 
!  NOTES:
!  (1 ) Routines READAIR and AIREMISS were originally written by Yuhang Wang,
!        1993.  These have been bundled into "aircraft_nox_mod.f" for easier
!        bookkeeping.  They have been kept mostly as-is, save for some
!        cosmetic changes and improved I/O error trapping. (bmy, 2/14/02)
!  (2 ) Updated comments, deleted some obsolete code (bmy, 3/8/02)
!  (3 ) Bug fix: only allocate arrays the first call to READAIR.
!        (yxw, bmy, 4/2/02)
!  (4 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and 
!        MODULE ROUTINES sections.  Updated comments.  Also deleted obsolete 
!        code. (bmy, 5/28/02)
!  (5 ) Now references "file_mod.f" (bmy, 6/26/02)
!  (6 ) Now references "pressure_mod.f".  Also deleted obsolete, commented-out
!        code from 6/02. (bmy, 8/20/02)
!  (7 ) Now references BXHEIGHT from "dao_mod.f". Now references "error_mod.f".
!        Also deleted obsolete code from various routines (bmy, 10/15/02)
!  (8 ) Now references "grid_mod.f" and "time_mod.f" (bmy, 2/11/03)
!  (9 ) Now references "directory_mod.f" (bmy, 7/19/04)
!  (10) Replace GEMISNOX (from CMN_NOX) with a module variable (bmy, 10/3/07)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "aircraft_nox_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: NAIR 
      PRIVATE :: LAIREMS
      PRIVATE :: AIR
      PRIVATE :: AIREMIS
      PRIVATE :: AIRPRESS

      ! PRIVATE module routines
      PRIVATE INIT_AIRCRAFT_NOX

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! NAIR     - Maximum number (km) for aircraft NOx emissions   
      INTEGER, PARAMETER   :: NAIR = 20
               
      ! LAIREMS GEOS-CHEM level where we will put emissions into
      INTEGER              :: LAIREMS
               
      ! AIR    aft NOx emissions on native 1-km grid 
      REAL*8,  ALLOCATABLE :: AIR(:,:,:)
               
      ! Aircraft NOx emissions on GEOS-CHEM grid
      REAL*8,  ALLOCATABLE :: AIREMIS(:,:,:)

      ! AIRPRESS - Approx. pressure edges of the 1-km native aircraft NOx grid
      REAL*8               :: AIRPRESS(NAIR+1) = (/  
     &   1013.25, 954.61, 845.56, 746.83, 657.64, 577.28, 505.07, 
     &    440.35, 382.52, 330.99, 285.24, 244.75, 209.04, 178.65, 
     &    152.59, 130.34, 111.33,  95.09,  81.22,  69.37,  59.26 /)

      ! Array to pass aircraft NOx emissions to SMVGEAR
      ! (replacement for GEMISNOX array)
      REAL*8,  ALLOCATABLE :: EMIS_AC_NOx(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READAIR
!
!******************************************************************************
!  Subroutine READAIR reads the aircraft NOx emissions from disk.
!  (yhw, bmy, 7/5/01, 2/11/03)
!
!  NOTES: 
!  (1 ) Now reference DATA_DIR from include file "CMN_SETUP".  Also reference
!        routine GET_RES_EXT from BPCH2_MOD. (bmy, 7/5/01)
!  (2 ) Now also use F90 syntax for declarations.  Also use IOERROR to trap
!        I/O errors.  Use D exponent to force REAL*8 precision. Also updated 
!        comments. (bmy, 7/5/01)
!  (3 ) Removed obsolete code from 7/01 (bmy, 9/4/01)
!  (4 ) Now read aircraft NOx files from the aircraft_NOx_200202/ subdirectory
!        of DATA_DIR.  Also updated comments. (bmy, 1/24/02)
!  (5 ) Now bundled into "aircraft_nox_mod.f" (bmy, 2/14/02)
!  (6 ) Updated comments (bmy, 3/8/02)
!  (7 ) Bug fix: only call INIT_AIRCRAFT_NOX to allocate arrays on the
!        first call to READAIR (yxw, bmy, 4/2/02)
!  (8 ) Deleted obsolete code (bmy, 5/28/02)
!  (9 ) Now use IU_FILE instead of IUNIT as the file unit number.  Also 
!        reference IU_FILE and IOERROR from "file_mod.f" (bmy, 6/26/02)
!  (10) Deleted obsolete code from 6/02. (bmy, 8/20/02)
!  (11) Now use function GET_MONTH from "time_mod.f".  Renamed INIT to
!        FIRST and MONTHSAVE to LASTMONTH. (bmy, 2/11/03)
!  (12) Removed reference to "CMN", it's obsolete.  Now references DATA_DIR 
!        from "directory_mod.f" (bmy, 7/19/04)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IOERROR
      USE TIME_MOD,      ONLY : GET_MONTH
      USE CMN_SIZE_MOD        ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: FIRST     = .TRUE.
      INTEGER, SAVE      :: LASTMONTH = -99
      INTEGER            :: I, J, K, IOS
      INTEGER            :: IU_FILE
      REAL*8             :: FUEL

      ! Conversion factor from [kg NO2/4h] to [molec NO2/s]
      REAL*8             :: CONV=9.0927d+20

      ! Month array
      CHARACTER(LEN=7)   :: MONTHSTR(12) = (/
     &                      'airjan.', 'airfeb.', 'airmar.', 'airapr.',
     &                      'airmay.', 'airjun.', 'airjul.', 'airaug.',
     &                      'airsep.', 'airoct.', 'airnov.', 'airdec.'/)

      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READAIR begins here!
      !
      ! NOTE: Aircraft NOx emissions have been stored in the input 
      ! files in units of [kg NO2/4h] instead of the more logical
      ! [kg NO2/day].  This is primarily historical.  Yuhang Wang
      ! did it this way since the old GISS-II model had a timestep of
      ! 4 hours.  These routines were then ported into the GEOS-CHEM
      ! as-is, thus leaving the units in [kg NO2/4h].  However, the
      ! emissions are converted below into [molec NO2/s], by applying
      ! the conversion factor CONV, and are then passed to SMVGEAR. 
      ! (bmy, 3/8/02)
      !=================================================================
      IF ( FIRST .or. GET_MONTH() /= LASTMONTH ) THEN

         ! Allocate and initialize arrays
         IF ( FIRST ) THEN 
            CALL INIT_AIRCRAFT_NOX
            FIRST = .FALSE.
         ENDIF

         ! Save month in LASTMONTH
         LASTMONTH = GET_MONTH()
         
         ! Zero emissions
         AIR(:,:,:) = 0d0

         ! Construct filename
         FILENAME = TRIM( DATA_DIR )        // 'aircraft_NOx_200202/' //
     &              MONTHSTR( GET_MONTH() ) // GET_RES_EXT()

         ! Find a free file LUN
         IU_FILE = findFreeLUN()

         ! Open file
         OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'readair:1' )

         ! Read aircraft NOx
         DO 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% KLUDGE FOR NORTH AMERICAN NESTED GRID (bmy, 7/2/12)
!%%% For some reason there seems to be bad data within the N. American
!%%% nested grid data files.  Some of the longitude indices are out of bounds
!%%% for the grid.  Put in a simple kludge to just skip over these boxes.
!%%% We are going to be reinventing how emissions get done in GEOS-Chem,
!%%% so it's probably not worth recreating the the file at this time.
!%%%
!%%%            READ( IU_FILE, *, IOSTAT=IOS ) I, J, K, AIR(I,J,K)
!%%%          
!%%%            ! IOS < 0 is end of file
!%%%            ! IOS > 0 is an I/O error
!%%%            IF ( IOS < 0 ) EXIT
!%%%            IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'readair:2' )
  
            ! Now read emissions into a scalar
            READ( IU_FILE, *, IOSTAT=IOS ) I, J, K, FUEL

            ! IOS < 0 is end of file
            ! IOS > 0 is an I/O error
            IF ( IOS < 0 ) EXIT
            IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'readair:2' )

#if   defined( NESTED_NA )        
            ! Error check, cycle if out of bounds
            IF ( I > IIPAR ) CYCLE
            IF ( J > JJPAR ) CYCLE
#endif   

            ! Save emissions into the AIR array
            AIR(I,J,K) = FUEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ENDDO

         ! Close file
         CLOSE( IU_FILE )

         ! Convert the units from [kg NO2/4hr] -> [molecules/s]
         DO K = 1, NAIR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( AIR(I,J,K) > 0.0 ) AIR(I,J,K) = AIR(I,J,K) * CONV
         ENDDO
         ENDDO
         ENDDO
      ENDIF
     
      ! Return to calling program
      END SUBROUTINE READAIR

!------------------------------------------------------------------------------

      SUBROUTINE AIREMISS( State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine AIREMISS interpolates the aircraft NOx emissions from the 1-km
!  native grid onto the given GEOS-CHEM grid. (bmy, 2/14/02, 10/3/07)
!
!  Original code from Yuhang Wang (1993).
!
!  NOTES: 
!  (1 ) Now bundled into "aircraft_nox_mod.f" (bmy, 2/14/02)
!  (2 ) Replace P(IOFF,JOFF) with P(I,J), since P is now declared to be
!        of size (IIPAR,JJPAR) instead of (IIPAR,JJPAR) (bmy, 2/14/02)
!  (3 ) AIR has to be dimensioned (IIPAR,JJPAR,LGLOB), since it contains
!        global emissions.  AIREMIS can be declared (IIPAR,JJPAR,LAIREMS),
!        since that way it will have the same horizontal dimensions as
!        the GEMISNOX array. (bmy, 2/14/02)
!  (4 ) Removed obsolete code (bmy, 3/8/02)
!  (5 ) Now reference GET_PEDGE from "pressure_mod.f", which returns the
!        correct "floating" pressure (bmy, 8/20/02)
!  (6 ) Now reference BXHEIGHT from "dao_mod.f". (bmy, 9/18/02)
!  (7 ) I0 and J0 are now local variables.  Now use functions GET_XOFFSET
!        and GET_YOFFSET from "grid_mod.f" (bmy, 2/11/03)
!  (8 ) Replace GEMISNOX (from "CMN_NOX") with module variable EMIS_AC_NOx
!        in order to avoid common block errors. (ltm, bmy, 10/3/07)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept State_Chm argument
!  31 May 2013 - R. Yantosca - Replace IDTNOx with IDTNO for DEVEL tag
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,           ONLY : AD32_ac
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT
#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNO
#endif
      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches
      USE CMN_MOD                  ! PTOP, SIGE, AVP

      INTEGER             :: I,    J,     IREF, JREF,  L,     K
      INTEGER             :: I0,   J0
      REAL*8              :: PLOW, PHIGH, XSUM, PAIR1, PAIR2, TMP
      
      ! External functions
      REAL*8, EXTERNAL    ::  BOXVL
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      
      !=================================================================
      ! AIREMISS begins here!
      !=================================================================

      ! Read aircraft NOx emissions
      CALL READAIR

      ! Zero accumulator array
      EMIS_AC_NOx = 0d0

      ! Get nested-grid offsets
      I0 = GET_XOFFSET()
      J0 = GET_YOFFSET()
      
      ! Loop over surface grid boxes
      DO J = 1, JJPAR
         JREF = J + J0
      DO I = 1, IIPAR
         IREF = I + I0
         
         !==============================================================
         ! Loop over tropospheric GEOS-CHEM levels 
         !==============================================================
         DO L = 1, LAIREMS

            ! PLOW  is the pressure at the bottom edge of sigma level L
            ! PHIGH is the pressure at the top    edge of sigma level L 
            PLOW  = GET_PEDGE(I,J,L)
            PHIGH = GET_PEDGE(I,J,L+1)

            ! Make sure PLOW is not smaller than AIRPRESS(1)
            IF ( L == 1 .AND. PLOW < AIRPRESS(1) ) PLOW = AIRPRESS(1)
            
            ! Initialize the summing variable
            XSUM = 0.0

            !===========================================================
            ! Loop over the native 1-km aircraft NOx grid layers
            !===========================================================
            DO K = 1, NAIR

               ! PAIR1 is the pressure at the bottom of 1-km grid layer K
               ! PAIR2 is the pressure at the top    of 1-km grid layer K
               PAIR1 = AIRPRESS(K)
               PAIR2 = AIRPRESS(K+1)
            
               ! Compute the fraction of each 1-km layer K that
               ! lies within the given GEOS-CHEM layer L
               IF ( PHIGH >= PAIR1 ) THEN
                  GOTO 10

               ELSE IF ( PLOW < PAIR1 .AND. PLOW > PAIR2 ) THEN
                  IF ( PHIGH < PAIR2 ) THEN
                     XSUM = XSUM + 
     &                  AIR(IREF,JREF,K) * (PLOW-PAIR2) / (PAIR1-PAIR2)
                  ELSE 
                     XSUM = XSUM + 
     &                  AIR(IREF,JREF,K) * (PLOW-PHIGH) / (PAIR1-PAIR2)
                  ENDIF

               ELSE IF ( PHIGH <  PAIR1 .AND. 
     &                   PHIGH >  PAIR2 .AND.
     *                   PLOW  >= PAIR1) THEN
                  XSUM = XSUM + 
     &               AIR(IREF,JREF,K) * (PAIR1-PHIGH) / (PAIR1-PAIR2)

               ELSE IF ( PHIGH <= PAIR2 .AND. PLOW >= PAIR1 ) THEN
                  XSUM = XSUM + AIR(IREF,JREF,K)

               ENDIF
            ENDDO

            ! Store XSUM into AIREMIS array
 10         CONTINUE


            AIREMIS(I,J,L) = XSUM

            !===========================================================
            ! Store nonzero AIREMIS into GEMISNOX array,
            ! which will then be passed into SMVGEAR
            !===========================================================
            IF ( AIREMIS(I,J,L) > 0.0 ) THEN

               ! Convert from [molec/s] to [molec/cm3/s]
               TMP = AIREMIS(I,J,L)  / BOXVL(I,J,L,State_Met)

               ! No aircraft emissions in the stratosphere (cdh, 4/25/2013)
               IF ( ITS_IN_THE_STRAT(I, J, L, State_Met) ) TMP = 0d0

               ! Store in GEMISNOX in [molec/cm3/s]
               EMIS_AC_NOx(I,J,L) = TMP

               ! ND32 -- save NOx in [molec/cm2], will convert to
               ! [molec/cm2/s] in subroutine "diag3.f" (bmy, 3/16/00)
               IF ( ND32 > 0 ) THEN
                  AD32_ac(I,J,L) = AD32_ac(I,J,L) +
     &                     ( TMP * State_Met%BXHEIGHT(I,J,L) * 1d2 )
               ENDIF
            ENDIF

         ENDDO  ! L
      ENDDO     ! I 
      ENDDO     ! J

#if defined( DEVEL )
      ! Superimpose aircraft NOx onto exiting tendency
      ! Assumes dimension conformity
      State_Chm%TRAC_TEND(:,:,:,IDTNO) = 
     &State_Chm%TRAC_TEND(:,:,:,IDTNO) + EMIS_AC_NOx 
#endif

      ! Return to calling program
      END SUBROUTINE AIREMISS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_AIRCRAFT_NOX
!
!******************************************************************************
!  Subroutine INIT_AIRCRAFT_NOX allocates and initializes module variables.
!  (bmy, 2/14/02, 10/3/07)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Now allocate EMIS_AC_NOx array (bmy, 10/3/07)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD
      
      ! Local variables
      INTEGER :: AS

      ! LAIREMS is only defined up to the GEOS-CHEM tropopause 
      LAIREMS = LLTROP

      ! AIR holds the aircraft NOx on the native 1-km grid (NAIR levels)
      ALLOCATE( AIR( IIPAR, JJPAR, NAIR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIR' )
      AIR = 0d0

      ! AIREMIS holds the aircraft NOx on the GEOS grid (LAIREMS levels)
      ALLOCATE( AIREMIS( IIPAR, JJPAR, LAIREMS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIREMS' )
      AIREMIS = 0d0

      ! AIREMIS holds the aircraft NOx on the GEOS grid (LAIREMS levels)
      ALLOCATE( EMIS_AC_NOx( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIREMS' )
      EMIS_AC_NOx = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_AIRCRAFT_NOX

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_AIRCRAFT_NOX
!
!******************************************************************************
!  Subroutine CLEANUP_AIRCRAFT_NOX deallocates module variables. (bmy, 2/14/02)
!
!  NOTES:  
!  (1 ) Now deallocate EMIS_AC_NOx array (bmy, 10/3/07)
!******************************************************************************
!
      IF ( ALLOCATED( AIR         ) ) DEALLOCATE( AIR         )
      IF ( ALLOCATED( AIREMIS     ) ) DEALLOCATE( AIREMIS     )
      IF ( ALLOCATED( EMIS_AC_NOX ) ) DEALLOCATE( EMIS_AC_NOX ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_AIRCRAFT_NOX

!------------------------------------------------------------------------------

      ! End of module
      END MODULE AIRCRAFT_NOX_MOD
