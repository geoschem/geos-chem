! $Id: tpcore_bc_mod.f,v 1.2 2010/02/02 16:57:51 bmy Exp $
      MODULE TPCORE_BC_MOD
!
!******************************************************************************
!  Module TPCORE_BC_MOD contains modules and variables which are needed to
!  save and read TPCORE nested-grid boundary conditions to/from disk.
!  (yxw, bmy, 3/4/03, 12/18/09)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) CLEAN_BC (LOGICAL )    : Flag which denotes if we will zero BC's
!  (2 ) I0_W     (INTEGER )    : Lon offset of 1x1 TPCORE REGION [# boxes]
!  (3 ) J0_W     (INTEGER )    : Lat offset of 1x1 TPCORE REGION [# boxes]
!  (4 ) IM_W     (INTEGER )    : Lon extent of 1x1 TPCORE REGION [# boxes]
!  (5 ) JM_W     (INTEGER )    : Lat extent of 1x1 TPCORE REGION [# boxes]
!  (6 ) I1_W     (INTEGER )    : Lower left-hand  lon index of 1x1 WINDOW 
!  (7 ) J1_W     (INTEGER )    : Lower left-hand  lat index of 1x1 WINDOW
!  (8 ) I2_W     (INTEGER )    : Upper right-hand lon index of 1x1 WINDOW
!  (9 ) J2_W     (INTEGER )    : Upper right-hand lat index of 1x1 WINDOW
!  (10) IGZD     (INTEGER )    : ???
!  (110 TS_BC    (INTEGER )    : Timestep for reading in BC's from disk [min]
!  (11) I0_BC    (INTEGER )    : Lon offset of 4x5 WINDOW REGION [# boxes]
!  (12) J0_BC    (INTEGER )    : Lat offset of 4x5 WINDOW REGION [# boxes]
!  (13) IM_BC    (INTEGER )    : Lon extent of 4x5 WINDOW REGION [# boxes]
!  (14) JM_BC    (INTEGER )    : Lat extent of 4x5 WINDOW REGION [# boxes] 
!  (15) I1_BC    (INTEGER )    : Lower left-hand  lon index of 4x5 WINDOW 
!  (16) J1_BC    (INTEGER )    : Lower left-hand  lat index of 4x5 WINDOW 
!  (17) I2_BC    (INTEGER )    : Upper right-hand lon index of 4x5 WINDOW 
!  (18) J2_BC    (INTEGER )    : Upper right-hand lat index of 4x5 WINDOW 
!  (19) BC       (REAL*4  )    : Array containing CU tracers on coarse grid
!                                / any input region coarse grid
!  (19b) BC_NA   (REAL*4  )    : Array containing NA tracers on coarse grid
!  (19c) BC_EU   (REAL*4  )    : Array containing EU tracers on coarse grid
!  (19d) BC_CH   (REAL*4  )    : Array containing CH tracers on coarse grid
!  (20) MAP1x1   (INTEGER )    : Mapping array from 1x1 -> 4x5 grid
!
!  Module Routines:
!  ============================================================================
!  (1 ) SET_CLEAN_BC           : Initializes the CLEAN_BC module variable
!  (2 ) OPEN_BC_FILE           : Opens a new boundary conditions file (R or W)
!  (3 ) SAVE_GLOBAL_TPCORE_BC  : Saves boundary conditions from GLOBAL run
!  (4 ) DO_WINDOW_TPCORE_BC    : Calls routine to clean or read window BC's 
!  (5 ) CLEAN_WINDOW_TPCORE_BC : Zeroes the nested-grid TPCORE boundary cond's
!  (6 ) READ_WINDOW_TPCORE_BC  : Reads nested-grid TPCORE BC's from disk
!  (7 ) GET_4x5_BC             : Returns the 4x5 BC at a given 1x1 location
!  (7b) GET_2x25_BC             : Returns the 2x2.5 BC at a given 1x1 location
!  (8 ) ITS_TIME_FOR_BC        : Returns T if it is time to read BC's
!  (9 ) INIT_TPCORE_BC         : Initalizes and allocates module variables
!  (10) CLEANUP_TPCORE_BC      : Deallocates module vaf
!
!  GEOS-CHEM modules referenced by tpcore_call_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module containing routines for bpch file I/O
!  (2 ) directory_mod.f        : Module containing GEOS-CHEM data & metfld dirs
!  (3 ) error_mod.f            : Module containing I/O error/NaN check routines
!  (4 ) file_mod.f             : Contains file unit numbers and error checks
!  (5 ) grid_mod.f             : Module containing horizontal grid information
!  (6 ) logical_mod.f          : Module containing GEOS-CHEM logical switches
!  (7 ) time_mod.f             : Module containing routines for date & time
!  (8 ) tracer_mod.f           : Module containing GEOS-CHEM tracer array STT
!
!  Reference Diagram: 
!  ============================================================================
!
!  <-------------------------------------- IGLOB ---------------------->
!
!  +-------------------------------------------------------------------+   ^
!  | GLOBAL REGION                                                     |   |
!  |                                                                   |   |
!  |                       <-------------- IIPAR ------------->        |   |
!  |                                                                   |   |
!  |                       +=================================[Y]  ^    |   |
!  |                       |  WINDOW REGION (met field size)  |   |    |   |
!  |                       |                                  |   |    |   |
!  |                       |      <------- IM_W ------->      |   |    |   |
!  |                       |      +--------------------+  ^   |   |    |   |
!  |                       |      |  TPCORE REGION     |  |   |   |    |   |
!  |                       |      |  (transport is     |  |   |   |    |   |
!  |<------- I0 ---------->|<---->|   done in this     | JM_W | JJPAR  | JGLOB
!  |                       | I0_W |   window!!!)       |  |   |   |    |   |
!  |                       |      |                    |  |   |   |    |   |
!  |                       |      +--------------------+  V   |   |    |   |
!  |                       |        ^                         |   |    |   |
!  |                       |        | J0_W                    |   |    |   |
!  |                       |        V                         |   |    |   |
!  |                      [X]=================================+   V    |   |
!  |                                ^                                  |   |
!  |                                | J0                               |   |
!  |                                V                                  |   |
! [1]------------------------------------------------------------------+   V
!
!  DIAGRAM NOTES:
!  (a) The outermost box ("Global Region") is the global grid size.  This 
!      region has IGLOB boxes in longitude and JGLOB boxes in latitude.  
!      The origin of the "Global Region" is at the south pole, at the 
!      lower left-hand corner (point [1]). 
!
!  (b) The next innermost box ("Window Region") is the nested-grid window.
!      This region has IIPAR boxes in longitude and JJPAR boxes in latitude.
!      This is the size of the trimmed met fields that will be used for
!      a 1 x 1 "nested-grid" simulation.  
!          
!  (c) The innermost region ("TPCORE Region") is the actual area in which
!      TPCORE transport will be performed.  Note that this region is smaller
!      than the "Window Region".  It is set up this way since a cushion of 
!      grid boxes is needed TPCORE Region for boundary conditions.
!
!  (d) I0 is the longitude offset (# of boxes) and J0 is the latitude offset
!      (# of boxes) which translate between the "Global Region" and the
!      "Window Region". 
!
!  (e) I0_W is the longitude offset (# of boxes), and J0_W is the latitude
!      offset (# of boxes) which translate between the "Window Region"
!      and the "TPCORE Region".  
!
!  (f) The lower left-hand corner of the "Window Region" (point [X]) has
!      longitude and latitude indices (I1_W, J1_W).  Similarly, the upper
!      right-hand corner (point [Y]) has longitude and latitude indices 
!      (I2_W, J2_W).
!
!  (g) Note that if I0=0, J0=0, I0_W=0, J0_W=0, IIPAR=IGLOB, JJPAR=JGLOB
!      specifies a global simulation.  In this case the "Window Region"
!      totally coincides with the "Global Region".  
!
!  (h) In order for the nested-grid to work we must save out concentrations
!      over the WINDOW REGION from a coarse model (e.g. 4x5) corresponding to
!      the same WINDOW REGION at 1x1.  These concentrations are copied along
!      the edges of the 1x1 WINDOW REGION and are thus used as boundary
!      conditions for TPCORE.  We assume that we will save out concentrations
!      from the 4x5 model since presently it takes too long to run at 2x25.
!
!  NOTES:
!  (1 ) Bug fix for LINUX w/ TIMESTAMP_STRING (bmy, 9/29/03)
!  (2 ) Now references "tracer_mod.f", "directory_mod.f", and
!        "logical_mod.f" (bmy, 7/20/04)
!  (3 ) Now get HALFPOLAR for GEOS or GCAP grids (bmy, 6/28/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Rename arguments in GET_4x5_BC to avoid name conflict (bmy, 10/24/05)
!  (6 ) Now use EXPAND_DATE instead of obsolete DATE_STRING (bmy, 3/15/06)
!  (7 ) Added 2x2.5 boundary condition output (created GET_2x25_BC).
!        Added multi-boundary condition output (NA, EU, CH and Custom region).
!        Unternally defined boundary condition regions for NA, EU and CH.
!        (amv, bmy, 12/18/09)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tpcore_bc_mod.f"
      !=================================================================
      
      ! Make everything PRIVATE ...
      PRIVATE
      
      ! ... except these variables ...
      PUBLIC :: I0_W, J0_W 
      PUBLIC :: I1_W, J1_W
      PUBLIC :: I2_W, J2_W
      PUBLIC :: IM_W, JM_W
      PUBLIC :: IGZD

      ! ... and these routines
      PUBLIC :: INIT_TPCORE_BC
      PUBLIC :: DO_WINDOW_TPCORE_BC
      PUBLIC :: SET_CLEAN_BC
      PUBLIC :: SAVE_GLOBAL_TPCORE_BC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL              :: CLEAN_BC
      INTEGER              :: I0_W,  J0_W,  IM_W,  JM_W   
      INTEGER              :: I1_W,  J1_W,  I2_W,  J2_W
      INTEGER              :: I0_BC, J0_BC, I1_BC, J1_BC
      INTEGER              :: I2_BC, J2_BC, IM_BC, JM_BC
      INTEGER              :: I0_BC_NA, J0_BC_NA, I1_BC_NA, J1_BC_NA
      INTEGER              :: I2_BC_NA, J2_BC_NA, IM_BC_NA, JM_BC_NA
      INTEGER              :: I0_BC_EU, J0_BC_EU, I1_BC_EU, J1_BC_EU
      INTEGER              :: I2_BC_EU, J2_BC_EU, IM_BC_EU, JM_BC_EU
      INTEGER              :: I0_BC_CH, J0_BC_CH, I1_BC_CH, J1_BC_CH
      INTEGER              :: I2_BC_CH, J2_BC_CH, IM_BC_CH, JM_BC_CH
      INTEGER              :: IGZD,  TS_BC
      INTEGER, ALLOCATABLE :: MAP1x1(:,:,:)
      REAL*4,  ALLOCATABLE :: BC(:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_NA(:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_EU(:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_CH(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
     
      SUBROUTINE SET_CLEAN_BC( THIS_CLEAN_BC )
!
!******************************************************************************
!  Subroutine SET_CLEAN_BC initializes the CLEAN_BC logical flag.  CLEAN_BC
!  decides whether or not we will zero the nested-grid tpcore boundary.
!  (bmy, 3/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THIS_CLEAN_BC (LOGICAL) : Logical value (T/F) to assign to CLEAN_BC
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      LOGICAL, INTENT(IN) :: THIS_CLEAN_BC

      !=================================================================
      ! SET_CLEAN_BC begins here!
      !=================================================================
      CLEAN_BC = THIS_CLEAN_BC

      ! Return to calling program
      END SUBROUTINE SET_CLEAN_BC

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_BC_FILE( FOR_READ, FOR_WRITE, WINDOW ) 
!
!******************************************************************************
!  Subroutine OPEN_BC_FILE opens the file which contains boundary conditions
!  saved from the coarse-grid WINDOW REGION for either reading or writing.
!  (bmy, 3/7/03, 12/18/09)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FOR_READ  (LOGICAL) : If passed, opens binary punch file for reading
!  (2 ) FOR_WRITE (LOGICAL) : If passed, opens binary punch file for writing  
!
!  NOTES:
!  (1 ) Now use ITS_A_NEW_DAY from "time_mod.f".  Now references TPBC_DIR
!        from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) DATE_STRING is now obsolete; use EXPAND_DATE instead (bmy, 3/15/06)
!  (4 ) Can now read files from different directories (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_WRITE
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_MOD, ONLY : TPBC_DIR,    TPBC_DIR_NA
      USE DIRECTORY_MOD, ONLY : TPBC_DIR_CH, TPBC_DIR_EU
      USE FILE_MOD,      ONLY : IU_BC, IU_BC_NA, IU_BC_EU, IU_BC_CH
      USE TIME_MOD,      ONLY : EXPAND_DATE, GET_NYMD, ITS_A_NEW_DAY

      IMPLICIT NONE

#     include "define.h"     ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN), OPTIONAL :: FOR_READ, FOR_WRITE
      INTEGER, INTENT(IN)           :: WINDOW

      ! Local variables
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! OPEN_BC_FILE begins here!
      !=================================================================

      ! Only open file if it's a new day
      IF ( ITS_A_NEW_DAY() ) THEN

         ! File name for BC's
         IF ( WINDOW .eq. 5 ) THEN
#if   defined( NESTED_NA )
            FILENAME = TRIM( TPBC_DIR_NA ) // 'BC.YYYYMMDD'
#elif defined( NESTED_EU )
            FILENAME = TRIM( TPBC_DIR_EU ) // 'BC.YYYYMMDD'
#elif defined( NESTED_CH )
            FILENAME = TRIM( TPBC_DIR_CH ) // 'BC.YYYYMMDD'
#endif
         ELSEIF ( WINDOW .eq. 1 ) THEN
            FILENAME = TRIM( TPBC_DIR ) // 'BC.YYYYMMDD'
         ELSEIF ( WINDOW .eq. 2 ) THEN
            FILENAME = TRIM( TPBC_DIR_NA ) // 'BC.YYYYMMDD'
         ELSEIF ( WINDOW .eq. 3 ) THEN
            FILENAME = TRIM( TPBC_DIR_EU ) // 'BC.YYYYMMDD'
         ELSEIF ( WINDOW .eq. 4 ) THEN
            FILENAME = TRIM( TPBC_DIR_CH ) // 'BC.YYYYMMDD'
         ENDIF
         
         ! Replace YYYYMMDD with the actual date
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), 000000 )

         ! Echo file name to stdout
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - OPEN_BC_FILE: Opening ', a )

         IF ( WINDOW .eq. 5 ) THEN
            IF ( PRESENT( FOR_READ ) )
     &         CALL OPEN_BPCH2_FOR_READ( IU_BC, 
     &               FILENAME )
         ELSEIF ( WINDOW .eq. 1 ) THEN
            IF ( PRESENT( FOR_WRITE ) )
     &        CALL OPEN_BPCH2_FOR_WRITE( IU_BC, 
     &               FILENAME )
         ELSEIF ( WINDOW .eq. 2 ) THEN
            IF ( PRESENT( FOR_WRITE ) )
     &        CALL OPEN_BPCH2_FOR_WRITE( IU_BC_NA, 
     &               FILENAME )
         ELSEIF ( WINDOW .eq. 3 ) THEN
            IF ( PRESENT( FOR_WRITE ) )
     &        CALL OPEN_BPCH2_FOR_WRITE( IU_BC_EU, 
     &               FILENAME )
         ELSEIF ( WINDOW .eq. 4 ) THEN
            IF ( PRESENT( FOR_WRITE ) )
     &        CALL OPEN_BPCH2_FOR_WRITE( IU_BC_CH, 
     &               FILENAME )
         ENDIF

      ENDIF

      ! Return to calling program
      END SUBROUTINE OPEN_BC_FILE

!------------------------------------------------------------------------------
      
      SUBROUTINE SAVE_GLOBAL_TPCORE_BC
!
!******************************************************************************
!  Subroutine SAVE_GLOBAL_TPCORE_BC saves concentrations from the WINDOW
!  REGION of a coarse-resolution model run (e.g. 4x5) to a bpch file.  
!  A new boundary conditions file is created for each day. 
!  (yxw, bmy, 3/4/03, 12/18/09)
!
!  NOTES:
!  (1 ) Now references N_TRACERS and STT from "tracer_mod.f".  Also now 
!        references TIMESTAMP_STRING from "time_mod.f".  (bmy, 7/20/04) 
!  (2 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids (bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Can now save files to different directories (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : BPCH2,     GET_HALFPOLAR,   GET_MODELNAME
      USE FILE_MOD,   ONLY : IU_BC, IU_BC_NA, IU_BC_EU, IU_BC_CH
      USE TIME_MOD,   ONLY : GET_NYMD,  GET_NHMS, 
     &                       GET_TAU,   TIMESTAMP_STRING
      USE TRACER_MOD, ONLY : N_TRACERS, STT
      USE LOGICAL_MOD, ONLY : LWINDO_CU, LWINDO_NA
      USE LOGICAL_MOD, ONLY : LWINDO_CH, LWINDO_EU

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      LOGICAL, SAVE      :: FIRST     = .TRUE.
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180 =  1
      INTEGER            :: I, IOS, J, L, N, IC
      REAL*4             :: LONRES, LATRES
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*8             :: TAU
      CHARACTER(LEN=16)  :: STAMP
      CHARACTER(LEN=20)  :: MODELNAME 
      CHARACTER(LEN=40)  :: CATEGORY  = 'IJ-AVG-$'
      CHARACTER(LEN=40)  :: UNIT      = 'v/v'
      CHARACTER(LEN=40)  :: RESERVED  = ''

      !=================================================================
      ! SAVE_GLOBAL_TPCORE_BC begins here!
      !=================================================================

      ! Return if it's not time to write data to disk
      IF ( .not. ITS_TIME_FOR_BC() ) RETURN

      !=================================================================
      ! Save boundary conditions from coarse grid to a BPCH file
      !=================================================================

      ! Define variables for BPCH output
      LONRES    = DISIZE
      LATRES    = DJSIZE
      HALFPOLAR = GET_HALFPOLAR()
      MODELNAME = GET_MODELNAME()
      TAU       = GET_TAU()

      DO IC = 1, 4
         IF ((IC .eq. 1) .and. LWINDO_CU )THEN

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TRACERS

               ! Save concentrations in WINDOW REGION to disk
               DO L = 1, LLPAR
                  BC(1:IM_BC,1:JM_BC,L,N) = 
     &               STT(I1_BC:I2_BC,J1_BC:J2_BC,L,N) 
               ENDDO

               ! Write boundary conditions to binary punch file
               CALL BPCH2( IU_BC,     MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      TAU,       TAU,      RESERVED,
     &               IM_BC,     JM_BC,     LLPAR,    I1_BC,
     &               J1_BC,     1,   BC(1:IM_BC, 1:JM_BC, 1:LLPAR, N) )

               ENDDO

               ! Echo info
               STAMP = TIMESTAMP_STRING()
               WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote BC''s at ', a )

         ELSEIF ((IC .eq. 2) .and. LWINDO_NA )THEN

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TRACERS

               ! Save concentrations in WINDOW REGION to disk
               DO L = 1, LLPAR
                  BC_NA(1:IM_BC_NA,1:JM_BC_NA,L,N) = 
     &               STT(I1_BC_NA:I2_BC_NA,J1_BC_NA:J2_BC_NA,L,N)
               ENDDO

               ! Write boundary conditions to binary punch file
               CALL BPCH2( IU_BC_NA,  MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,   TAU,       TAU,      RESERVED,
     &               IM_BC_NA,  JM_BC_NA,   LLPAR,  I1_BC_NA,
     &               J1_BC_NA, 1, BC_NA(1:IM_BC_NA, 1:JM_BC_NA, 
     &               1:LLPAR, N) )

            ENDDO

            ! Echo info
            STAMP = TIMESTAMP_STRING()
            WRITE( 6, 111 ) STAMP
 111  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote NA BC''s at ', a )

         ELSEIF ((IC .eq. 3) .and. LWINDO_EU )THEN

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TRACERS

               ! Save concentrations in WINDOW REGION to disk
               DO L = 1, LLPAR
                  BC_EU(1:IM_BC_EU,1:JM_BC_EU,L,N) =
     &               STT(I1_BC_EU:I2_BC_EU,J1_BC_EU:J2_BC_EU,L,N)
               ENDDO
 
               ! Write boundary conditions to binary punch file
               CALL BPCH2( IU_BC_EU,  MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,   TAU,       TAU,      RESERVED,
     &               IM_BC_EU,  JM_BC_EU,   LLPAR,  I1_BC_EU,
     &               J1_BC_EU, 1, BC_EU(1:IM_BC_EU, 1:JM_BC_EU, 
     &               1:LLPAR, N) )

            ENDDO

            ! Echo info
            STAMP = TIMESTAMP_STRING()
            WRITE( 6, 112 ) STAMP
 112  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote EU BC''s at ', a )

         ELSEIF ((IC .eq. 4) .and. LWINDO_CH )THEN

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TRACERS

               ! Save concentrations in WINDOW REGION to disk
               DO L = 1, LLPAR
                  BC_CH(1:IM_BC_CH,1:JM_BC_CH,L,N) =
     &               STT(I1_BC_CH:I2_BC_CH,J1_BC_CH:J2_BC_CH,L,N)
               ENDDO
 
               ! Write boundary conditions to binary punch file
               CALL BPCH2( IU_BC_CH,  MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,   TAU,       TAU,      RESERVED,
     &               IM_BC_CH,  JM_BC_CH,   LLPAR,  I1_BC_CH,
     &               J1_BC_CH, 1, BC_CH(1:IM_BC_CH, 1:JM_BC_CH,
     &               1:LLPAR, N) )

            ENDDO

            ! Echo info
            STAMP = TIMESTAMP_STRING()
            WRITE( 6, 113 ) STAMP
 113  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote CH BC''s at ', a )

         ENDIF
      ENDDO

      ! Return to calling program  
      END SUBROUTINE SAVE_GLOBAL_TPCORE_BC

!------------------------------------------------------------------------------

      SUBROUTINE DO_WINDOW_TPCORE_BC
!
!******************************************************************************
!  Subroutine DO_WINDOW_TPCORE_BC is a driver routine for assigning TPCORE
!  boundary conditions to the tracer array STT.  (bmy, 3/7/03, 12/18/09)
!
!  At present, assume that we have saved 
!
!  NOTES:
!  (1 ) Now references N_TRACERS and STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now can use 2 x 2.5 BC's (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : N_TRACERS, STT
      USE LOGICAL_MOD, ONLY : LWINDO2x25

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, J, L, N

      !=================================================================
      ! DO_WINDOW_TPCORE_BC begins here!
      !=================================================================

      ! Either zero BC's or read them from disk
      IF ( CLEAN_BC ) THEN
         CALL CLEAN_WINDOW_TPCORE_BC
      ELSE
         CALL READ_WINDOW_TPCORE_BC
      ENDIF

      !=================================================================
      ! Copy the BC data into the proper elements of STT
      !
      ! NOTE: We assume that we have saved 4x5 BC's from a global 
      ! GEOS-CHEM run.  It takes too long to save the 2x25 BC's.  
      ! One may always another subroutine for 2x25 BC's later.
      !=================================================================

      IF ( .not. LWINDO2x25 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR

           ! First loop over all latitudes (WINDOW REGION)
           DO J = 1, JJPAR

              ! West BC
              DO I = 1, I0_W
                 STT(I,J,L,N) = GET_4x5_BC(I,J,L,N)
              ENDDO

              ! East BC
              DO I = IIPAR-I0_W+1, IIPAR
                 STT(I,J,L,N) = GET_4x5_BC(I,J,L,N)
              ENDDO
           ENDDO

           ! Now loop over the longitudes of the TPCORE REGION
           DO I = 1+I0_W, IM_W+I0_W
              
              ! South BC
              DO J = 1, J0_W
                 STT(I,J,L,N) = GET_4x5_BC(I,J,L,N)
              ENDDO

              ! North BC
              DO J = JJPAR-J0_W+1, JJPAR
                 STT(I,J,L,N) = GET_4x5_BC(I,J,L,N)
              ENDDO
           ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR

           ! First loop over all latitudes (WINDOW REGION)
           DO J = 1, JJPAR

              ! West BC
              DO I = 1, I0_W
                 STT(I,J,L,N) = GET_2x25_BC(I,J,L,N)
              ENDDO

              ! East BC
              DO I = IIPAR-I0_W+1, IIPAR
                 STT(I,J,L,N) = GET_2x25_BC(I,J,L,N)
              ENDDO
           ENDDO

           ! Now loop over the longitudes of the TPCORE REGION
           DO I = 1+I0_W, IM_W+I0_W
   
              ! South BC
              DO J = 1, J0_W
                 STT(I,J,L,N) = GET_2x25_BC(I,J,L,N)
              ENDDO

              ! North BC
              DO J = JJPAR-J0_W+1, JJPAR
                 STT(I,J,L,N) = GET_2x25_BC(I,J,L,N)
              ENDDO
           ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
      ! Return to calling program
      END SUBROUTINE DO_WINDOW_TPCORE_BC
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEAN_WINDOW_TPCORE_BC
!
!******************************************************************************
!  Subroutine CLEAN_WINDOW_TPCORE_BC zeroes the boundary conditions array
!  BC at each timestep.  (bmy, 3/7/03, 12/18/09) 
!
!  NOTES:
!  (1 ) Now references N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now zeroes the arrays for the different regions (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : N_TRACERS
      USE LOGICAL_MOD,   ONLY : LWINDO_NA, LWINDO_EU
      USE LOGICAL_MOD,   ONLY : LWINDO_CH, LWINDO_CU

#     include "CMN_SIZE"  ! Size parameters
#     include "define.h"

      ! Local variables
      INTEGER :: I, J, L, N

      !=================================================================
      ! CLEAN_WINDOW_TPCORE_BC begins here!
      !=================================================================
#if   defined(NESTED_CH) || defined(NESTED_NA) || defined(NESTED_EU)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JM_BC
         DO I = 1, IM_BC
            BC(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
#endif

      IF ( LWINDO_CU ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JM_BC
         DO I = 1, IM_BC
            BC(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF ( LWINDO_NA ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JM_BC_NA
         DO I = 1, IM_BC_NA
            BC_NA(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF ( LWINDO_EU ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JM_BC_EU
         DO I = 1, IM_BC_EU
            BC_EU(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF 
 
      IF ( LWINDO_CH ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JM_BC_CH
         DO I = 1, IM_BC_CH
            BC_CH(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CLEAN_WINDOW_TPCORE_BC

!------------------------------------------------------------------------------

      SUBROUTINE READ_WINDOW_TPCORE_BC
!
!******************************************************************************
!  Subroutine READ_WINDOW_TPCORE_BC reads tracer concentrations saved on the
!  WINDOW REGION of a coarse-grid simulation (e.g. 4x5).  These concentrations
!  will be used as boundary conditions for TPCORE transport. 
!  (bmy, 3/7/03, 12/18/09)
!
!  NOTES:
!  (1 ) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (2 ) Now references N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Rewritten to be more generic (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD,   ONLY : IOERROR, IU_BC
      USE TIME_MOD,   ONLY : GET_TAU, TIMESTAMP_STRING
      USE TRACER_MOD, ONLY : N_TRACERS

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      INTEGER           :: I,      IOS,     J,         L      
      INTEGER           :: NI,     NJ,      NL,        NT
      INTEGER           :: NFOUND, IFIRST,  JFIRST,    LFIRST
      INTEGER           :: NSKIP,  NTRACER, HALFPOLAR, CENTER180
      REAL*4            :: LONRES, LATRES
      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*8            :: ZTAU0,  ZTAU1
      CHARACTER(LEN=20) :: S = 'read_window_boundary'  
      CHARACTER(LEN=20) :: MODELNAME 
      CHARACTER(LEN=40) :: CATEGORY  
      CHARACTER(LEN=40) :: UNIT      
      CHARACTER(LEN=40) :: RESERVED  

      ! For LINUX fix
      CHARACTER(LEN=16) :: STAMP

      !=================================================================
      ! READ_WINDOW_TPCORE_BC begins here!
      !=================================================================
      IF ( ITS_TIME_FOR_BC() ) THEN

         ! Open boundary conditions file (if necessary) 
         CALL OPEN_BC_FILE( FOR_READ=.TRUE., WINDOW=5 )

         ! Initialize # of tracers found
         NFOUND = 0

         ! Loop
         DO 

            !===========================================================
            ! Read each data block one at a time
            !===========================================================

            ! Read first header line
            READ( IU_BC, IOSTAT=IOS ) 
     &           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT
         
            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS, IU_BC, S//':1' )

            ! Read second header line
            READ( IU_BC, IOSTAT=IOS ) 
     &           CATEGORY, NT, UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ, NL,   IFIRST, JFIRST, LFIRST, NSKIP
            
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_BC, S//':2' )

            ! Read data array
            READ( IU_BC, IOSTAT=IOS ) 
     &           ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_BC, S//':3' )
            
            !===========================================================
            ! If this is the right time, then save into BC array
            !===========================================================
            IF ( GET_TAU() == ZTAU0 ) THEN

               ! Copy into the BC array
               DO L = 1, NL
                  BC( 1:IM_BC, 1:JM_BC, L, NT ) = ARRAY(1:NI, 1:NJ, L)
               ENDDO

               ! Increment count of tracers found
               NFOUND = NFOUND + 1
            ENDIF

            !===========================================================
            ! Exit if we've found all tracers for this TAU value
            !===========================================================
            IF ( NFOUND == N_TRACERS ) THEN

               ! Echo output
               STAMP = TIMESTAMP_STRING()
               WRITE( 6, 100 ) N_TRACERS, STAMP
 100           FORMAT( '     - READ_WINDOW_TPCORE_BC: Found all ',
     &                       i3, ' BC''s at ', a )
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE READ_WINDOW_TPCORE_BC

!------------------------------------------------------------------------------

      FUNCTION GET_4x5_BC( I_1x1, J_1x1, L_1x1, N_1x1 ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_4x5_BC returns a value from the 4x5 BC boundary conditions
!  array at the location of a 1x1 grid box. (yxw, bmy, 3/7/03)
!
!  For now we will assume that we have saved tracer concentrations from a 
!  4x5 window whioh overlays the corresponding 1x1 WINDOW REGION.  These 4x5
!  tracer concentrations are used as boundary conditions for TPCORE.  
!  
!  We assume that we won't be saving 2x2.5 boundary conditions anytime in the 
!  near future since it currently takes too long to run a 2x2.5 full-chemistry
!  simulation for a whole year.  However, if we ever do save 2x25 boundary
!  conditions, one may always write a corresponding subroutine GET_2x25_BC to 
!  do the mapping from 1x1 -> 2x25.  (yxw, bmy, 3/7/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I_1x1 (INTEGER) : Longitude index of 1x1 grid box
!  (2 ) J_1x1 (INTEGER) : Latitude  index of 1x1 grid box
!  (3 ) L_1x1 (INTEGER) : Altitude  index of 1x1 grid box
!  (4 ) N_1x1 (INTEGER) : Tracer    index of 1x1 grid box
!
!  NOTES:
!  (1 ) Rename arguments to avoid conflict w/ I1x1, J1x1 parameters in 
!        CMN_SIZE. (bmy, 10/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! NTRACE

      ! Arguments
      INTEGER, INTENT(IN)  :: I_1x1, J_1x1, L_1x1, N_1x1

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, II, J, JJ, JJ1
      REAL*8               :: X,         Y
      REAL*8,  SAVE        :: XE4x5(73), YE4x5(47)

      ! Function value
      REAL*8               :: VALUE

      !=================================================================
      ! GET_4x5_BC begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         !==============================================================
         ! Define lon & lat edges of the 4x5 grid -- we need to 
         ! use these to do the mapping from 1x1 --> 4x5 
         !==============================================================

         ! Lon edges
         DO I = 1, 73
            XE4x5(I) = -182.5d0 + 5 * ( I - 1 )
         ENDDO

         ! Lat edges
         DO J = 2, 46
            YE4x5(J) = -92d0 + 4 * ( J - 1 ) 
         ENDDO

         ! Polar lats
         YE4x5(1)  = -90d0
         YE4x5(47) = +90d0

         !==============================================================
         ! Locate the 4x5 box(es) to which each 1x1 box belongs
         ! X, Y are lon & lat centers of the 1x1 boxes in degrees
         ! Save in the MAP1x1 array for future reference
         !==============================================================
         DO J = 1, JJPAR
            Y = GET_YMID( J )  

            DO I = 1, IIPAR
               X = GET_XMID( I )

               ! Loop over 4x5 longitudes
               DO II = 1, 72

                  ! If the 1x1 center lon lies w/in the 4x5 lon edges
                  ! then we have found the proper 4x5 box!
                  IF ( X > XE4x5(II) .and. X < XE4x5(II+1) ) THEN
                     MAP1x1(I,J,1) = II
                     EXIT
                  ENDIF
               ENDDO

               ! Loop over 4x5 latitudes
               DO JJ = 1, 46
            
                  ! If the 1x1 lat center lies between the 4x5 lat
                  ! edges, we have found the proper 4x5 box!
                  IF ( Y > YE4x5(JJ) .and. Y < YE4x5(JJ+1) ) THEN
                     MAP1x1(I,J,2) = JJ
                     MAP1x1(I,J,3) = 0
                     EXIT

                  ! If the 1x1 lat center equals the 4x5 lower lat
                  ! edge, then we need to average this 4x5 box and
                  ! the box just south of it
                  ELSE IF ( Y == YE4x5(JJ) ) THEN
                     MAP1x1(I,J,2) = JJ-1
                     MAP1x1(I,J,3) = JJ
                     EXIT

                  ! If the 1x1 lat center equals the 4x5 lower lat
                  ! edge, then we need to average this 4x5 box and
                  ! the box just north of it                 
                  ELSE IF ( Y == YE4x5(JJ+1) ) THEN
                     MAP1x1(I,J,2) = JJ
                     MAP1x1(I,J,3) = JJ+1
                     EXIT

                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=============================================================
      ! Locate the tracer concentration at the 4x5 box which 
      ! corresponds to the 1x1 box (I_1x1, J_1x1, L_1x1, N_1x1)
      !=============================================================

      ! Get lon indices
      II  = MAP1x1( I_1x1, J_1x1, 1 ) - I0_BC

      ! Get lat indices
      JJ  = MAP1x1( I_1x1, J_1x1, 2 ) - J0_BC
      JJ1 = MAP1x1( I_1x1, J_1x1, 3 ) - J0_BC

      ! Locate the 4x5 box(es) corresponding to the 1x1 box
      ! If our 1x1 box straddles 2 4x5 boxes, average the 4x5 values
      IF ( MAP1x1( I_1x1, J_1x1, 3 ) > 0 ) THEN
         VALUE = 0.5 * ( BC(II,JJ,L_1x1,N_1x1) + BC(II,JJ1,L_1x1,N_1x1))
      ELSE
         VALUE = BC(II,JJ,L_1x1,N_1x1)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_4x5_BC

!------------------------------------------------------------------------------

      FUNCTION GET_2x25_BC( I_1x1, J_1x1, L_1x1, N_1x1 ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_2x25_BC returns a value from the 2x2.5 BC boundary conditions
!  array at the location of a 1x1 grid box. (amv, bmy, 12/18/09)
!
!  For now we will assume that we have saved tracer concentrations from a
!  2x25 window which overlays the corresponding 1x1 WINDOW REGION.  These 2x2.5
!  tracer concentrations are used as boundary conditions for TPCORE.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I_1x1 (INTEGER) : Longitude index of 1x1 grid box
!  (2 ) J_1x1 (INTEGER) : Latitude  index of 1x1 grid box
!  (3 ) L_1x1 (INTEGER) : Altitude  index of 1x1 grid box
!  (4 ) N_1x1 (INTEGER) : Tracer    index of 1x1 grid box
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! NTRACE

      ! Arguments
      INTEGER, INTENT(IN)  :: I_1x1, J_1x1, L_1x1, N_1x1

      ! Local variables
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, II, J, JJ, JJ1
      REAL*8               :: X,         Y
      REAL*8,  SAVE        :: XE2x25(145), YE2x25(92)

      ! Function value
      REAL*8               :: VALUE

      !=================================================================
      ! GET_2x25_BC begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         !==============================================================
         ! Define lon & lat edges of the 2x25 grid -- we need to
         ! use these to do the mapping from 1x1 --> 2x25
         !==============================================================

         ! Lon edges
         DO I = 1, 145
            XE2x25(I) = -181.25d0 + 2.5 * ( I - 1 )
         ENDDO

         ! Lat edges
         DO J = 2, 91
            YE2x25(J) = -91d0 + 2 * ( J - 1 )
         ENDDO

         ! Polar lats
         YE2x25(1)  = -90d0
         YE2x25(92) = +90d0

         !==============================================================
         ! Locate the 2x25 box(es) to which each 1x1 box belongs
         ! X, Y are lon & lat centers of the 1x1 boxes in degrees
         ! Save in the MAP1x1 array for future reference
         !==============================================================
         DO J = 1, JJPAR
            Y = GET_YMID( J )

            DO I = 1, IIPAR
               X = GET_XMID( I )

               ! Loop over 2x25 longitudes
               DO II = 1, 144

                  ! If the 1x1 center lon lies w/in the 2x25 lon edges
                  ! then we have found the proper 2x25 box!
                  IF ( X > XE2x25(II) .and. X < XE2x25(II+1) ) THEN
                     MAP1x1(I,J,1) = II
                     EXIT
                  ENDIF
               ENDDO

               ! Loop over 2x25 latitudes
               DO JJ = 1, 91

                  ! If the 1x1 lat center lies between the 2x25 lat
                  ! edges, we have found the proper 2x25 box!
                  IF ( Y > YE2x25(JJ) .and. Y < YE2x25(JJ+1) ) THEN
                     MAP1x1(I,J,2) = JJ
                     MAP1x1(I,J,3) = 0
                     EXIT

                  ! If the 1x1 lat center equals the 2x25 lower lat
                  ! edge, then we need to average this 2x25 box and
                  ! the box just south of it
                  ELSE IF ( Y == YE2x25(JJ) ) THEN
                     MAP1x1(I,J,2) = JJ-1
                     MAP1x1(I,J,3) = JJ
                     EXIT

                  ! If the 1x1 lat center equals the 2x25 lower lat
                  ! edge, then we need to average this 2x25 box and
                  ! the box just north of it
                  ELSE IF ( Y == YE2x25(JJ+1) ) THEN
                     MAP1x1(I,J,2) = JJ
                     MAP1x1(I,J,3) = JJ+1
                     EXIT

                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=============================================================
      ! Locate the tracer concentration at the 4x5 box which
      ! corresponds to the 1x1 box (I_1x1, J_1x1, L_1x1, N_1x1)
      !=============================================================

      ! Get lon indices
      II  = MAP1x1( I_1x1, J_1x1, 1 ) - I0_BC

      ! Get lat indices
      JJ  = MAP1x1( I_1x1, J_1x1, 2 ) - J0_BC
      JJ1 = MAP1x1( I_1x1, J_1x1, 3 ) - J0_BC

      ! Locate the 2x25 box(es) corresponding to the 1x1 box
      ! If our 1x1 box straddles 2 2x25 boxes, average the 2x25 values
      IF ( MAP1x1( I_1x1, J_1x1, 3 ) > 0 ) THEN
         VALUE = 0.5 * ( BC(II,JJ,L_1x1,N_1x1) + BC(II,JJ1,L_1x1,N_1x1))
      ELSE
         VALUE = BC(II,JJ,L_1x1,N_1x1)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_2x25_BC

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_BC() RESULT( FLAG )
!
!******************************************************************************
!  Subroutine ITS_TIME_FOR_BC returns TRUE if it is time to read in the next
!  set of boundary conditions for TPCORE, or FALSE otherwise.  (bmy, 3/5/03)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_ELAPSED_MIN
      
      ! Function value
      LOGICAL :: FLAG

      !=================================================================
      ! ITS_TIME_FOR_BC begins here!
      !=================================================================
      FLAG = ( MOD( GET_ELAPSED_MIN(), TS_BC ) == 0 ) 
     
      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_BC

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TPCORE_BC( TS, I0W, J0W, I1, J1, I2, J2 )
!
!******************************************************************************
!  Subroutine INIT_TPCORE_BC initializes module variables and arrays 
!  (bmy, 2/10/03, 7/20/04)
!
!  NOTES:
!  (1 ) Now references N_TRACERS from "tracer_mod.f".  Now references LWINDO
!        from "logical_mod.f".  Now references TPBC_DIR from "directory_mod.f".
!        Now references ITS_A_NESTED_GRID from "grid_mod.f".  Also added 
!        arguments to take values from "input_mod.f". (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : TPBC_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE GRID_MOD,      ONLY : GET_XOFFSET, GET_YOFFSET 
      USE GRID_MOD,      ONLY : ITS_A_NESTED_GRID
      USE LOGICAL_MOD,   ONLY : LWINDO, LWINDO_NA, LWINDO_CU
      USE LOGICAL_MOD,   ONLY : LWINDO_EU, LWINDO_CH
      USE LOGICAL_MOD,   ONLY : LWINDO2x25
      USE TRACER_MOD,    ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN) :: TS, I0W, J0W, I1, J1, I2, J2

      ! Local variables
      INTEGER             :: AS, I, J

      !=================================================================
      ! INIT_TPCORE_BC begins here!
      !=================================================================

      ! Timestep for BC's [min]
      TS_BC = TS

      !---------------------
      ! *** NESTED GRID ***
      !---------------------

      ! TPCORE transport region offsets
      I0_W = I0W
      J0_W = J0W

      ! Extent of TPCORE transport region
      IM_W = IIPAR - ( 2 * I0_W )
      JM_W = JJPAR - ( 2 * J0_W )

      ! I and J at lower-left corner of TPCORE transport region
      I1_W = GET_XOFFSET( GLOBAL=.TRUE. ) + I0_W + 1
      J1_W = GET_YOFFSET( GLOBAL=.TRUE. ) + J0_W + 1

      ! I and J at upper-right corner of TPCORE transport region
      I2_W = GET_XOFFSET( GLOBAL=.TRUE. ) + IIPAR - I0_W
      J2_W = GET_YOFFSET( GLOBAL=.TRUE. ) + JJPAR - J0_W

      ! IGZD = ?
      IGZD = I0_W - 1

      !---------------------
      ! *** GLOBAL GRID ***
      !---------------------

      ! Lower-left corner of coarse-grid BC WINDOW region 
      I1_BC = I1
      J1_BC = J1
  
      ! Upper-right corner of coarse-grid BC WINDOW region 
      I2_BC = I2
      J2_BC = J2

      ! Extent of coarse-grid BC REGION
      IM_BC = I2_BC - I1_BC + 1
      JM_BC = J2_BC - J1_BC + 1

      ! TPCORE internal transport window offset
      I0_BC = I1_BC - 1
      J0_BC = J1_BC - 1

      IF (.not. LWINDO2x25) THEN
#if defined(GRID4x5) || defined(NESTED_NA)
         ! Lower-left corner of coarse-grid NA BC WINDOW region
         I1_BC_NA = 9
         J1_BC_NA = 26
   
         ! Upper-right corner of coarse-grid NA BC WINDOW region
         I2_BC_NA = 29
         J2_BC_NA = 41
#endif
      ELSE
#if defined(GRID2x25) || defined(NESTED_NA)
         ! Lower-left corner of coarse-grid NA BC WINDOW region
         I1_BC_NA = 17
         J1_BC_NA = 51
  
         ! Upper-right corner of coarse-grid NA BC WINDOW region
         I2_BC_NA = 57
         J2_BC_NA = 81
#endif
      ENDIF
 
      ! Extent of coarse-grid NA BC REGION
      IM_BC_NA = I2_BC_NA - I1_BC_NA + 1
      JM_BC_NA = J2_BC_NA - J1_BC_NA + 1
 
      ! TPCORE NA internal transport window offset
      I0_BC_NA = I1_BC_NA - 1
      J0_BC_NA = J1_BC_NA - 1

      IF (.not. LWINDO2x25) THEN
#if   defined(GRID4x5) || defined(NESTED_EU)
         ! Lower-left corner of coarse-grid EU BC WINDOW region
         I1_BC_EU = 31
         J1_BC_EU = 31
   
         ! Upper-right corner of coarse-grid EU BC WINDOW region
         I2_BC_EU = 47
         J2_BC_EU = 41
#endif
      ELSE
#if defined(GRID2x25) || defined(NESTED_EU)
         ! Lower-left corner of coarse-grid EU BC WINDOW region
         I1_BC_EU = 61
         J1_BC_EU = 61

         ! Upper-right corner of coarse-grid EU BC WINDOW region
         I2_BC_EU = 93
         J2_BC_EU = 81
#endif
      ENDIF
 
      ! Extent of coarse-grid EU BC REGION
      IM_BC_EU = I2_BC_EU - I1_BC_EU + 1
      JM_BC_EU = J2_BC_EU - J1_BC_EU + 1
 
      ! TPCORE NA internal transport window offset
      I0_BC_EU = I1_BC_EU - 1
      J0_BC_EU = J1_BC_EU - 1

      IF (.not. LWINDO2x25) THEN 
#if   defined(GRID4x5) || defined(NESTED_CH)
         ! Lower-left corner of coarse-grid CH BC WINDOW region
         I1_BC_CH = 51
         J1_BC_CH = 21
   
         ! Upper-right corner of coarse-grid CH BC WINDOW region
         I2_BC_CH = 67
         J2_BC_CH = 37
#endif
      ELSE
#if defined(GRID2x25) || defined(NESTED_CH)
      ! Lower-left corner of coarse-grid CH BC WINDOW region
      I1_BC_CH = 101
      J1_BC_CH = 40

      ! Upper-right corner of coarse-grid CH BC WINDOW region
      I2_BC_CH = 133
      J2_BC_CH = 74
#endif
      ENDIF
       
      ! Extent of coarse-grid EU BC REGION
      IM_BC_CH = I2_BC_CH - I1_BC_CH + 1
      JM_BC_CH = J2_BC_CH - J1_BC_CH + 1
 
      ! TPCORE NA internal transport window offset
      I0_BC_CH = I1_BC_CH - 1
      J0_BC_CH = J1_BC_CH - 1

      ! Return if we are not saving 4x5 BC's
      ! or if it's not a nested grid simulation
      IF ( .not. ITS_A_NESTED_GRID() ) THEN
         IF ( .not. LWINDO ) RETURN
      ENDIF

      !=================================================================
      ! Allocate and initialize arrays 
      !=================================================================
      IF ( LWINDO_CU ) THEN

         ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
         BC = 0e0

      ENDIF

      IF ( LWINDO_CH ) THEN

         ALLOCATE( BC_CH( IM_BC_CH, JM_BC_CH, LLPAR, N_TRACERS )
     &                   , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_CH' )
         BC_CH = 0e0

      ENDIF

      IF ( LWINDO_NA ) THEN

         ALLOCATE( BC_NA( IM_BC_NA, JM_BC_NA, LLPAR, N_TRACERS )
     &                   , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_NA' )
         BC_NA = 0e0

      ENDIF

      IF ( LWINDO_EU ) THEN

         ALLOCATE( BC_EU( IM_BC_EU, JM_BC_EU, LLPAR, N_TRACERS )
     &            , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_EU' )
         BC_EU = 0e0

      ENDIF

      ! If running a nested simulation, just allocate the custom BC
      ! array with nested dimensions
#if   defined( NESTED_CH )
      I1_BC = I1_BC_CH
      J1_BC = J1_BC_CH
      I2_BC = I2_BC_CH
      J2_BC = J2_BC_CH
      IM_BC = IM_BC_CH
      JM_BC = JM_BC_CH
      I0_BC = I0_BC_CH
      J0_BC = J0_BC_CH

      ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS )
     &             , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
      BC = 0e0
#elif defined( NESTED_NA )
      I1_BC = I1_BC_NA
      J1_BC = J1_BC_NA
      I2_BC = I2_BC_NA
      J2_BC = J2_BC_NA
      IM_BC = IM_BC_NA
      JM_BC = JM_BC_NA
      I0_BC = I0_BC_NA
      J0_BC = J0_BC_NA

      ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS )
     &             , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
      BC = 0e0
#elif defined( NESTED_EU )
      I1_BC = I1_BC_EU
      J1_BC = J1_BC_EU
      I2_BC = I2_BC_EU
      J2_BC = J2_BC_EU
      IM_BC = IM_BC_EU
      JM_BC = JM_BC_EU
      I0_BC = I0_BC_EU
      J0_BC = J0_BC_EU

      ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS )
     &             , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
      BC = 0e0
#endif

      ALLOCATE( MAP1x1( IIPAR, JJPAR, 3 ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAP1x1' )
      MAP1x1 = 0

      ! Return to calling program
      END SUBROUTINE INIT_TPCORE_BC

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TPCORE_BC
!
!******************************************************************************
!  Subroutine CLEANUP_TPCORE_BC deallocates all module arrays. (3/4/03)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TPCORE_BC begins here!
      !=================================================================
      IF ( ALLOCATED( BC     ) ) DEALLOCATE( BC     )
      IF ( ALLOCATED( BC_NA  ) ) DEALLOCATE( BC_NA  )
      IF ( ALLOCATED( BC_EU  ) ) DEALLOCATE( BC_EU  )
      IF ( ALLOCATED( BC_CH  ) ) DEALLOCATE( BC_CH  )
      IF ( ALLOCATED( MAP1x1 ) ) DEALLOCATE( MAP1x1 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TPCORE_BC

!------------------------------------------------------------------------------

      END MODULE TPCORE_BC_MOD
