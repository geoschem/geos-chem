!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tpcore_bc_mod
!
! !DESCRIPTION: Module TPCORE\_BC\_MOD contains modules and variables 
!  which are needed to save and read TPCORE nested-grid boundary conditions 
!  to/from disk.
!\\
!\\
! !INTERFACE: 
!
      MODULE TPCORE_BC_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      !-----------------------------------------------------------
      ! I0_W : Lon offset of TPCORE REGION [# boxes]
      ! J0_W : Lat offset of TPCORE REGION [# boxes]
      ! IM_W : Lon extent of TPCORE REGION [# boxes]
      ! JM_W : Lat extent of TPCORE REGION [# boxes]
      ! I1_W : Lower left-hand  (LL) lon index of NESTED WINDOW 
      ! J1_W : Lower left-hand  (LL) lat index of NESTED WINDOW
      ! I2_W : Upper right-hand (UR) lon index of NESTED WINDOW
      ! J2_W : Upper right-hand (UR) lat index of NESTED WINDOW
      ! IGZD : ???
      !
      ! Please also see the diagram in the REMARKS section.
      !------------------------------------------------------------
      INTEGER, PUBLIC :: I0_W, J0_W  
      INTEGER, PUBLIC :: I1_W, J1_W   
      INTEGER, PUBLIC :: I2_W, J2_W
      INTEGER, PUBLIC :: IM_W, JM_W
      INTEGER, PUBLIC :: IGZD
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: INIT_TPCORE_BC
      PUBLIC :: DO_WINDOW_TPCORE_BC
      PUBLIC :: SET_CLEAN_BC
      PUBLIC :: SAVE_GLOBAL_TPCORE_BC
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: OPEN_BC_FILE
      PRIVATE :: CLEAN_WINDOW_TPCORE_BC
      PRIVATE :: READ_WINDOW_TPCORE_BC
      PRIVATE :: GET_4x5_BC
      PRIVATE :: GET_2x25_BC
      PRIVATE :: ITS_TIME_FOR_BC
      PRIVATE :: CLEANUP_TPCORE_BC
!
! !REMARKS:
!  Reference Diagram: 
!  ===========================================================================
!                                                                             .
!  <-------------------------------------- IIPAR ---------------------->
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
!  |<------- I0 ---------->|<---->|   done in this     | JM_W | JJPAR  | JJPAR
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
!                                                                             .
!  DIAGRAM NOTES:
!  (a) The outermost box ("Global Region") is the global grid size.  This 
!      region has IIPAR boxes in longitude and JJPAR boxes in latitude.  
!      The origin of the "Global Region" is at the south pole, at the 
!      lower left-hand corner (point [1]). 
!                                                                             .
!  (b) The next innermost box ("Window Region") is the nested-grid window.
!      This region has IIPAR boxes in longitude and JJPAR boxes in latitude.
!      This is the size of the trimmed met fields that will be used for
!      a 1 x 1 "nested-grid" simulation.  
!                                                                             .
!  (c) The innermost region ("TPCORE Region") is the actual area in which
!      TPCORE transport will be performed.  Note that this region is smaller
!      than the "Window Region".  It is set up this way since a cushion of 
!      grid boxes is needed TPCORE Region for boundary conditions.
!                                                                             .
!  (d) I0 is the longitude offset (# of boxes) and J0 is the latitude offset
!      (# of boxes) which translate between the "Global Region" and the
!      "Window Region". 
!                                                                             .
!  (e) I0_W is the longitude offset (# of boxes), and J0_W is the latitude
!      offset (# of boxes) which translate between the "Window Region"
!      and the "TPCORE Region".  
!                                                                             .
!  (f) The lower left-hand corner of the "Window Region" (point [X]) has
!      longitude and latitude indices (I1_W, J1_W).  Similarly, the upper
!      right-hand corner (point [Y]) has longitude and latitude indices 
!      (I2_W, J2_W).
!                                                                             .
!  (g) Note that if I0=0, J0=0, I0_W=0, J0_W=0, IIPAR=IIPAR, JJPAR=JJPAR
!      specifies a global simulation.  In this case the "Window Region"
!      totally coincides with the "Global Region".  
!                                                                             .
!  (h) In order for the nested-grid to work we must save out concentrations
!      over the WINDOW REGION from a coarse model (e.g. 4x5) corresponding to
!      the same WINDOW REGION at 1x1.  These concentrations are copied along
!      the edges of the 1x1 WINDOW REGION and are thus used as boundary
!      conditions for TPCORE.
!
! !REVISION HISTORY:
!  04 Mar 2003 - R. Yantosca - Initial version
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
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - Now make IU_BC, IU_BC_NA, IU_BC_EU, IU_BC_CH,
!                              and IU_BC_SE local variables for findFreeLUN
!  07 Sep 2012 - R. Yantosca - Minor fixes for numerical stability
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! CLEAN_BC : Flag which denotes if we will zero BC's
      ! TS_BC    : Timestep for reading in BC's from disk [min]
      ! I0_BC    : Lon offset of 4x5 WINDOW REGION [# boxes]
      ! J0_BC    : Lat offset of 4x5 WINDOW REGION [# boxes]
      ! IM_BC    : Lon extent of 4x5 WINDOW REGION [# boxes]
      ! JM_BC    : Lat extent of 4x5 WINDOW REGION [# boxes] 
      ! I1_BC    : Lower left-hand  lon index of 4x5 WINDOW 
      ! J1_BC    : Lower left-hand  lat index of 4x5 WINDOW 
      ! I2_BC    : Upper right-hand lon index of 4x5 WINDOW 
      ! J2_BC    : Upper right-hand lat index of 4x5 WINDOW
      ! BC       : Array containing CU tracers on coarse grid
      !            / any input region coarse grid
      ! BC_NA    : Array containing NA tracers on coarse grid
      ! BC_EU    : Array containing EU tracers on coarse grid
      ! BC_CH    : Array containing CH tracers on coarse grid
      ! MAP1x1   : Mapping array from 1x1 -> 4x5 grid
      !=================================================================

      ! Scalars
      LOGICAL              :: CLEAN_BC
      INTEGER              :: I0_BC,    J0_BC,    I1_BC,    J1_BC
      INTEGER              :: I2_BC,    J2_BC,    IM_BC,    JM_BC
      INTEGER              :: I0_BC_NA, J0_BC_NA, I1_BC_NA, J1_BC_NA
      INTEGER              :: I2_BC_NA, J2_BC_NA, IM_BC_NA, JM_BC_NA
      INTEGER              :: I0_BC_EU, J0_BC_EU, I1_BC_EU, J1_BC_EU
      INTEGER              :: I2_BC_EU, J2_BC_EU, IM_BC_EU, JM_BC_EU
      INTEGER              :: I0_BC_CH, J0_BC_CH, I1_BC_CH, J1_BC_CH
      INTEGER              :: I2_BC_CH, J2_BC_CH, IM_BC_CH, JM_BC_CH
      INTEGER              :: I0_BC_SE, J0_BC_SE, I1_BC_SE, J1_BC_SE
      INTEGER              :: I2_BC_SE, J2_BC_SE, IM_BC_SE, JM_BC_SE
      INTEGER              :: TS_BC

      ! Logical unit numbers
      INTEGER              :: IU_BC      ! LUN for input  BC's
      INTEGER              :: IU_BC_NA   ! LUN for output BC's, N. Amer. grid
      INTEGER              :: IU_BC_EU   ! LUN for output BC's, Europe grid
      INTEGER              :: IU_BC_CH   ! LUN for output BC's, China grid
      INTEGER              :: IU_BC_SE   ! LUN for output BC's, SEAC4RS grid
 
      ! Arrays
      INTEGER, ALLOCATABLE :: MAP1x1(:,:,:  )
      REAL*4,  ALLOCATABLE :: BC    (:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_NA (:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_EU (:,:,:,:)
      REAL*4,  ALLOCATABLE :: BC_CH (:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine SET\_CLEAN\_BC initializes the CLEAN\_BC logical 
!  flag.  CLEAN\_BC decides whether or not we will zero the nested-grid 
!  tpcore boundary conditions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CLEAN_BC( THIS_CLEAN_BC )
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: THIS_CLEAN_BC
! 
! !REVISION HISTORY: 
!  04 Mar 2003 - R. Yantosca - Initial versioni
!  15 May 2012 - R. Yantosca - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! SET_CLEAN_BC begins here!
      !=================================================================
      CLEAN_BC = THIS_CLEAN_BC

      END SUBROUTINE SET_CLEAN_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bc_file
!
! !DESCRIPTION: Subroutine OPEN\_BC\_FILE opens the file which contains 
!  boundary conditions saved from the coarse-grid WINDOW REGION for 
!  either reading or writing.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPEN_BC_FILE( FOR_READ, FOR_WRITE, WINDOW ) 
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_WRITE
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_MOD, ONLY : TPBC_DIR,     TPBC_DIR_NA
      USE DIRECTORY_MOD, ONLY : TPBC_DIR_CH,  TPBC_DIR_EU
      USE DIRECTORY_MOD, ONLY : TPBC_DIR_SE
      USE inquireMod,    ONLY : findFreeLUN
      USE TIME_MOD,      ONLY : EXPAND_DATE,  GET_NYMD
      USE TIME_MOD,      ONLY : ITS_A_NEW_DAY
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN), OPTIONAL :: FOR_READ   ! Open file for read?
      LOGICAL, INTENT(IN), OPTIONAL :: FOR_WRITE  ! Open file for write?
      INTEGER, INTENT(IN)           :: WINDOW     ! Specifies nested region 
! 
! !REMARKS:
!   Pass these values via the WINDOW argument to do the following actions;
!     WINDOW = 1 : Save BC's to file (Custom window,        aka "CU")
!     WINDOW = 2 : Save BC's to file (North America window, aka "NA")
!     WINDOW = 3 : Save BC's to file (Europe window,        aka "EU")
!     WINDOW = 4 : Save BC's to file (China/SE Asia window, aka "CH")
!     WINDOW = 5 : Read BC's from file
!
! !REVISION HISTORY: 
!  07 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now use ITS_A_NEW_DAY from "time_mod.f".  Now references TPBC_DIR
!        from "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (3 ) DATE_STRING is now obsolete; use EXPAND_DATE instead (bmy, 3/15/06)
!  (4 ) Can now read files from different directories (amv, bmy, 12/18/09)
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  06 Aug 2012 - R. Yantosca - Cleaned up IF statement, added comments
!  06 Aug 2012 - R. Yantosca - Close existing files before opening new files
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Removed SEAC4RS C-preprocessor switch
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL            :: IS_OPEN    ! Is the file already open?
      CHARACTER(LEN=255) :: FILENAME   ! Name of file to read

      !=================================================================
      ! OPEN_BC_FILE begins here!
      !=================================================================

      ! Only open file if it's a new day
      IF ( ITS_A_NEW_DAY(.TRUE.) ) THEN

         ! File name for BC's
         IF ( WINDOW .eq. 6 ) THEN
#if   defined( NESTED_NA )
            FILENAME = TRIM( TPBC_DIR_NA ) // 'BC.YYYYMMDD' ! Input file, NA
#elif defined( NESTED_EU )
            FILENAME = TRIM( TPBC_DIR_EU ) // 'BC.YYYYMMDD' ! Input file, EU
#elif defined( NESTED_CH )
            FILENAME = TRIM( TPBC_DIR_CH ) // 'BC.YYYYMMDD' ! Input file, CH
#endif
         ELSEIF ( WINDOW .eq. 1 ) THEN
            FILENAME = TRIM( TPBC_DIR )    // 'BC.YYYYMMDD' ! Output file, CU
         ELSEIF ( WINDOW .eq. 2 ) THEN
            FILENAME = TRIM( TPBC_DIR_NA ) // 'BC.YYYYMMDD' ! Output file, NA
         ELSEIF ( WINDOW .eq. 3 ) THEN
            FILENAME = TRIM( TPBC_DIR_EU ) // 'BC.YYYYMMDD' ! Output file, EU
         ELSEIF ( WINDOW .eq. 4 ) THEN
            FILENAME = TRIM( TPBC_DIR_CH ) // 'BC.YYYYMMDD' ! Output file, CH
         ENDIF
         
         ! Replace YYYYMMDD with the actual date
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), 000000 )

         ! FORMAT string for use below
 100     FORMAT( '     - OPEN_BC_FILE: Opening ', a, ' on unit ', i4 )

         IF ( WINDOW .eq. 6 ) THEN

            !----------------------------------------------------------
            ! WINDOW == 6 : Read BC's from file
            !----------------------------------------------------------
            IF ( PRESENT( FOR_READ ) ) THEN

               ! Close previously opened file
               INQUIRE( IU_BC, OPENED=IS_OPEN )
               IF ( IS_OPEN ) CLOSE( IU_BC )

               ! Find a free LUN
               IU_BC = findFreeLun()

               ! Open file for input
               WRITE( 6, 100 ) TRIM( FILENAME ), IU_BC
               CALL OPEN_BPCH2_FOR_READ( IU_BC, FILENAME )
            ENDIF

         ELSE IF ( WINDOW .eq. 1 ) THEN

            !----------------------------------------------------------
            ! WINDOW == 1 : Save BC's to file (custom window)
            !----------------------------------------------------------

            IF ( PRESENT( FOR_WRITE ) ) THEN

               ! Close previously opened file
               INQUIRE( IU_BC, OPENED=IS_OPEN )
               IF ( IS_OPEN ) CLOSE( IU_BC )

               ! Find a free LUN
               IU_BC = findFreeLun()

               ! Open file for output
               WRITE( 6, 100 ) TRIM( FILENAME ), IU_BC
               CALL OPEN_BPCH2_FOR_WRITE( IU_BC,    FILENAME )
            ENDIF

         ELSE IF ( WINDOW .eq. 2 ) THEN

            !----------------------------------------------------------
            ! WINDOW == 2 : Save BC's to file (North America window)
            !----------------------------------------------------------
            IF ( PRESENT( FOR_WRITE ) ) THEN

               ! Close previously opened file
               INQUIRE( IU_BC_NA, OPENED=IS_OPEN )
               IF ( IS_OPEN ) CLOSE( IU_BC_NA )

               ! Find a free LUN
               IU_BC_NA = findFreeLun()

               ! Open file for output
               WRITE( 6, 100 ) TRIM( FILENAME ), IU_BC_NA
               CALL OPEN_BPCH2_FOR_WRITE( IU_BC_NA, FILENAME )
            ENDIF

         ELSE IF ( WINDOW .eq. 3 ) THEN

            !----------------------------------------------------------
            ! WINDOW == 3 : Save BC's to file (Europe window)
            !----------------------------------------------------------
            IF ( PRESENT( FOR_WRITE ) ) THEN

               ! Close previously opened file
               INQUIRE( IU_BC_EU, OPENED=IS_OPEN )
               IF ( IS_OPEN ) CLOSE( IU_BC_EU )

               ! Find a free LUN
               IU_BC_EU = findFreeLun()

               ! Open file for output
               WRITE( 6, 100 ) TRIM( FILENAME ), IU_BC_EU
               CALL OPEN_BPCH2_FOR_WRITE( IU_BC_EU, FILENAME )
            ENDIF

         ELSE IF ( WINDOW .eq. 4 ) THEN

            !----------------------------------------------------------
            ! WINDOW == 4 : Save BC's to file (China/SE Asia window)
            !----------------------------------------------------------
            IF ( PRESENT( FOR_WRITE ) ) THEN

               ! Close previously opened file
               INQUIRE( IU_BC_CH, OPENED=IS_OPEN )
               IF ( IS_OPEN ) CLOSE( IU_BC_CH )

               ! Find a free LUN
               IU_BC_CH = findFreeLun()

               ! Open file for output
               WRITE( 6, 100 ) TRIM( FILENAME ), IU_BC_CH
               CALL OPEN_BPCH2_FOR_WRITE( IU_BC_CH, FILENAME )
            ENDIF

         ENDIF

      ENDIF

      END SUBROUTINE OPEN_BC_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: save_global_tpcore_bc
!
! !DESCRIPTION: Subroutine SAVE\_GLOBAL\_TPCORE\_BC saves concentrations from 
!  the WINDOW REGION of a coarse-resolution model run to a bpch file.  
!  A new boundary conditions file is created for each day. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SAVE_GLOBAL_TPCORE_BC( am_I_Root, Input_Opt, 
     &                                  State_Chm, RC         )
!
! !USES:
!
      USE BPCH2_MOD
      USE CMN_SIZE_MOD 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_NYMD,  GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU,   TIMESTAMP_STRING
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  04 Mar 2003 - Y. Wang     - Initial version
!  (1 ) Now references N_TRACERS and STT from "tracer_mod.f".  Also now 
!        references TIMESTAMP_STRING from "time_mod.f".  (bmy, 7/20/04) 
!  (2 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids (bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Can now save files to different directories (amv, bmy, 12/18/09)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  06 Aug 2012 - R. Yantosca - Added comments & cosmetic changes
!  06 Aug 2012 - R. Yantosca - Now use file unit variables from this module
!                              instead of from GeosUtil/file_mod.F
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST     = .TRUE.
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180 =  1
      INTEGER            :: I, IOS, J, L, N, IC, N_TOT_TRC
      REAL*4             :: LONRES, LATRES
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*8             :: TAU
      CHARACTER(LEN=16)  :: STAMP
      CHARACTER(LEN=20)  :: MODELNAME 
      CHARACTER(LEN=40)  :: CATEGORY  = 'IJ-AVG-$'
      CHARACTER(LEN=40)  :: UNIT      = 'v/v'
      CHARACTER(LEN=40)  :: RESERVED  = ''

      ! For fields from Input_Opt
      LOGICAL            :: LWINDO_CU
      LOGICAL            :: LWINDO_NA 
      LOGICAL            :: LWINDO_CH
      LOGICAL            :: LWINDO_EU
      LOGICAL            :: LWINDO_SE
      INTEGER            :: N_TRACERS

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER    :: STT(:,:,:,:)

      !=================================================================
      ! SAVE_GLOBAL_TPCORE_BC begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LWINDO_CU = Input_Opt%LWINDO_CU
      LWINDO_NA = Input_Opt%LWINDO_NA
      LWINDO_CH = Input_Opt%LWINDO_CH
      LWINDO_EU = Input_Opt%LWINDO_EU
      LWINDO_SE = Input_Opt%LWINDO_SE
      N_TRACERS = Input_Opt%N_TRACERS

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over the possible combinations of nested grids
      DO IC = 1, 4
         IF ((IC .eq. 1) .and. LWINDO_CU )THEN

            !----------------------------------------------------------
            ! IC == 1 : Save BC's to file (custom window)
            !----------------------------------------------------------

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TOT_TRC

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

         ELSEIF ((IC .eq. 2) .and. LWINDO_NA )THEN

            !----------------------------------------------------------
            ! IC == 2 : Save BC's to file (North America window)
            !----------------------------------------------------------

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TOT_TRC

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

         ELSEIF ((IC .eq. 3) .and. LWINDO_EU )THEN

            !----------------------------------------------------------
            ! IC == 3 : Save BC's to file (Europe window)
            !----------------------------------------------------------

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TOT_TRC

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

         ELSEIF ((IC .eq. 4) .and. LWINDO_CH )THEN

            !----------------------------------------------------------
            ! IC == 4 : Save BC's to file (China/SE Asia window)
            !----------------------------------------------------------

            ! Open file for writing, if necessary
            CALL OPEN_BC_FILE( FOR_WRITE=.TRUE., WINDOW=IC)

            ! Loop over each tracer
            DO N = 1, N_TOT_TRC

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
            
         ENDIF
      ENDDO

      ! Free pointer
      NULLIFY( STT )

      ! FORMAT statements
 110  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote BC''s at ',    a )
 111  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote NA BC''s at ', a )
 112  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote EU BC''s at ', a )
 113  FORMAT( '     - SAVE_GLOBAL_TPCORE_BC: Wrote CH BC''s at ', a )

      END SUBROUTINE SAVE_GLOBAL_TPCORE_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_window_tpcore_bc
!
! !DESCRIPTION: Subroutine DO\_WINDOW\_TPCORE\_BC is a driver routine for 
!  assigning TPCORE boundary conditions to the tracer array STT. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, 
     &                                State_Chm, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  07 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references N_TRACERS and STT from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now can use 2 x 2.5 BC's (amv, bmy, 12/18/09)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE   :: FIRST = .TRUE.
      INTEGER         :: I, J, L, N, N_TOT_TRC

      ! For fields from Input_Opt
      LOGICAL         :: LWINDO2x25
      INTEGER         :: N_TRACERS

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! DO_WINDOW_TPCORE_BC begins here!
      !=================================================================

      ! Assume success
      RC         =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LWINDO2x25 = Input_Opt%LWINDO2x25
      N_TRACERS  = Input_Opt%N_TRACERS
      
      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Either zero BC's or read them from disk
      IF ( CLEAN_BC ) THEN
         CALL CLEAN_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, RC )
      ELSE
         CALL READ_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, RC )
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
         DO N = 1, N_TOT_TRC
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
         DO N = 1, N_TOT_TRC
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

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE DO_WINDOW_TPCORE_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clean_window_tpcore_bc
!
! !DESCRIPTION: Subroutine CLEAN\_WINDOW\_TPCORE\_BC zeroes the boundary 
!  conditions array BC at each timestep.  (bmy, 3/7/03, 12/18/09) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEAN_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD  
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! !REVISION HISTORY: 
!  07 Mar 2003 - M. Prather  - Initial version
!  (1 ) Now references N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now zeroes the arrays for the different regions (amv, bmy, 12/18/09)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  07 Sep 2012 - R. Yantosca - Simplify coding, remove parallel loops
!  07 Sep 2012 - R. Yantosca - Now use 0d0 instead of 0e0 to zero BC arrays
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Removed SEAC4RS C-preprocessor switch
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N, N_TOT_TRC

      ! For fields from Input_Opt
      LOGICAL :: LWINDO_NA, LWINDO_EU
      LOGICAL :: LWINDO_CH, LWINDO_CU
      LOGICAL :: LWINDO_SE
      INTEGER :: N_TRACERS

      !=================================================================
      ! CLEAN_WINDOW_TPCORE_BC begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LWINDO_NA = Input_Opt%LWINDO_NA
      LWINDO_EU = Input_Opt%LWINDO_EU
      LWINDO_CH = Input_Opt%LWINDO_CH
      LWINDO_CU = Input_Opt%LWINDO_CU
      LWINDO_SE = Input_Opt%LWINDO_SE
      N_TRACERS = Input_Opt%N_TRACERS

#if   defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_EU )
      BC = 0d0
#endif

      IF ( LWINDO_CU ) THEN
         BC = 0d0
      ENDIF

      IF ( LWINDO_NA ) THEN
         BC_NA = 0d0
      ENDIF

      IF ( LWINDO_EU ) THEN
         BC_EU = 0d0
      ENDIF 
 
      IF ( LWINDO_CH ) THEN
         BC_CH = 0d0
      ENDIF

      END SUBROUTINE CLEAN_WINDOW_TPCORE_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_window_tpcore_bc
!
! !DESCRIPTION: Subroutine READ\_WINDOW\_TPCORE\_BC reads tracer 
!  concentrations saved on the WINDOW REGION of a coarse-grid simulation 
!  (e.g. 4x5, 2x2.5).  These concentrations will be used as boundary 
!  conditions for TPCORE transport. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD 
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : GET_TAU, TIMESTAMP_STRING
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  07 Mar 2003 - R. Yantosca - Initial version
!  (1 ) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (2 ) Now references N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Rewritten to be more generic (amv, bmy, 12/18/09)
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      ! For fields from Input_Opt
      INTEGER           :: N_TRACERS

      !=================================================================
      ! READ_WINDOW_TPCORE_BC begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS

      IF ( ITS_TIME_FOR_BC() ) THEN

         ! Open boundary conditions file (if necessary) 
         CALL OPEN_BC_FILE( FOR_READ=.TRUE., WINDOW=6 )

         ! Initialize # of tracers found
         NFOUND = 0

         ! Loop
         DO 

            !===========================================================
            ! Read each data block one at a time
            !===========================================================

            ! Zero array, for safety's sake
            ARRAY = 0e0

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
               IF ( am_I_Root ) THEN
                  WRITE( 6, 100 ) N_TRACERS, STAMP
 100              FORMAT( '     - READ_WINDOW_TPCORE_BC: Found all ',
     &                       i3, ' BC''s at ', a )
               ENDIF
               EXIT
            ENDIF
         ENDDO
      ENDIF

      END SUBROUTINE READ_WINDOW_TPCORE_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_4x5_bc
!
! !DESCRIPTION: Function GET\_4x5\_BC returns a value from the 4x5 BC 
!  boundary conditions array at the location of a nested grid box.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_4x5_BC( I_1x1, J_1x1, L_1x1, N_1x1 ) RESULT( VALUE )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GRID_MOD,      ONLY : GET_XMID, GET_YMID
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I_1x1   ! Nested-grid lon    index
      INTEGER, INTENT(IN) :: J_1x1   ! Nested-grid lat    index
      INTEGER, INTENT(IN) :: L_1x1   ! Nested-grid level  index
      INTEGER, INTENT(IN) :: N_1x1   ! Nested-grid tracer index
!
! !RETURN VALUE:
!
      REAL*8              :: VALUE   ! 4 x 5 BC @ location of nested grid box
!
! !REMARKS:
!  NOTE: This routine was originally written for the 1 x 1 nested grid,
!  but this now works for the GEOS-5 0.5 x 0.666 nested grid data.  Keep
!  variable names I_1x1, J_1x1, etc. unchanged for the present.
!                                                                             .
!  For now we will assume that we have saved tracer concentrations from 
!  a  4x5 window which overlays the corresponding 1x1 WINDOW REGION.  
!  These 4x5 tracer concentrations are used as boundary conditions for 
!  TPCORE.  
! 
! !REVISION HISTORY: 
!  07 Mar 2003 - Y. Wang, R. Yantosca - Initial version
!  (1 ) Rename arguments to avoid conflict w/ I1x1, J1x1 parameters in 
!        CMN_SIZE. (bmy, 10/24/05)
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      !clb
      INTEGER       :: I, II, II1, J, JJ, JJ1
      REAL*8        :: X,         Y
      REAL*8,  SAVE :: XE4x5(73), YE4x5(47)

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
         DO I = 1, IIPAR

            X = GET_XMID( I, J, 1 )
            Y = GET_YMID( I, J, 1 )  

            ! Loop over 4x5 longitudes
            DO II = 1, 72

               ! If the 1x1 center lon lies w/in the 4x5 lon edges
               ! then we have found the proper 4x5 box!
               IF ( X > XE4x5(II) .and. X < XE4x5(II+1) ) THEN
                  MAP1x1(I,J,1) = II
                  MAP1x1(I,J,2) = 0
                  EXIT
               ! Bug fix for 4x5 BCs (clb, 5/29/12)
               ELSE IF ( X == XE4x5(II) .and. I == 1) THEN
                  MAP1x1(I,J,1) = II
                  MAP1x1(I,J,2) = 0
                  EXIT
               ELSE IF ( X == XE4x5(II) .and. I .ne. 1) THEN
                  MAP1x1(I,J,1) = II - 1
                  MAP1x1(I,J,2) = II
                  EXIT
               ELSE IF ( X == XE4x5(II+1) ) THEN
                  MAP1x1(I,J,1) = II
                  MAP1x1(I,J,2) = II + 1
                  EXIT
               ENDIF

            ENDDO

            ! Loop over 4x5 latitudes
            DO JJ = 1, 46
            
               ! If the 1x1 lat center lies between the 4x5 lat
               ! edges, we have found the proper 4x5 box!
               IF ( Y > YE4x5(JJ) .and. Y < YE4x5(JJ+1) ) THEN
                  MAP1x1(I,J,3) = JJ
                  MAP1x1(I,J,4) = 0
                  EXIT
                  
               ! If the 1x1 lat center equals the 4x5 lower lat
               ! edge, then we need to average this 4x5 box and
               ! the box just south of it
               ! Bug fix for 4x5 BCs (clb, 5/29/12)
               ELSE IF ( Y == YE4x5(JJ) .and. J == 1) THEN
                  MAP1x1(I,J,3) = JJ
                  MAP1x1(I,J,4) = 0
                  EXIT

               ELSE IF ( Y == YE4x5(JJ) .and. J .ne. 1) THEN
                  MAP1x1(I,J,3) = JJ-1
                  MAP1x1(I,J,4) = JJ
                  EXIT

               ! If the 1x1 lat center equals the 4x5 lower lat
               ! edge, then we need to average this 4x5 box and
               ! the box just north of it              
   
               ELSE IF ( Y == YE4x5(JJ+1) ) THEN
                  MAP1x1(I,J,3) = JJ
                  MAP1x1(I,J,4) = JJ+1
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
      II1 = MAP1x1( I_1x1, J_1x1, 2 ) - I0_BC

      ! Get lat indices
      JJ  = MAP1x1( I_1x1, J_1x1, 3 ) - J0_BC
      JJ1 = MAP1x1( I_1x1, J_1x1, 4 ) - J0_BC

      ! Locate the 4x5 box(es) corresponding to the 1x1 box
      ! If our 1x1 box straddles 2 4x5 boxes, average the 4x5 values
      IF (MAP1x1(I_1x1,J_1x1,2)>0 .and. MAP1x1(I_1x1,J_1x1,4)>0) THEN
         VALUE = 0.25 * (BC(II,JJ,L_1x1,N_1x1) + BC(II,JJ1,L_1x1,N_1x1)
     &               + BC(II1,JJ,L_1x1,N_1x1)+ BC(II1,JJ1,L_1x1,N_1x1))
      ELSE IF (MAP1x1(I_1x1,J_1x1,2)>0 ) THEN
         VALUE = 0.5 *( BC(II,JJ,L_1x1,N_1x1) + BC(II1,JJ,L_1x1,N_1x1))
      ELSE IF (MAP1x1(I_1x1,J_1x1,4)>0 ) THEN
         VALUE = 0.5 *( BC(II,JJ,L_1x1,N_1x1) + BC(II,JJ1,L_1x1,N_1x1))
      ELSE
         VALUE = BC(II,JJ,L_1x1,N_1x1)
      ENDIF

      END FUNCTION GET_4x5_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_2x25_bc
!
! !DESCRIPTION: Function GET\_2x25\_BC returns a value from the 2 x 2.5 
!  BC boundary conditions array at the location of a nested grid box. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_2x25_BC( I_1x1, J_1x1, L_1x1, N_1x1 ) RESULT( VALUE )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GRID_MOD,     ONLY : GET_XMID, GET_YMID
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I_1x1   ! Nested-grid lon    index
      INTEGER, INTENT(IN) :: J_1x1   ! Nested-grid lat    index
      INTEGER, INTENT(IN) :: L_1x1   ! Nested-grid level  index
      INTEGER, INTENT(IN) :: N_1x1   ! Nested-grid tracer index
!
! !RETURN VALUE:
!
      REAL*8              :: VALUE   ! 2 x 2.5 BC @ location of nested grid box
!
! !REMARKS:
!  NOTE: This routine was originally written for the 1 x 1 nested grid,
!  but this now works for the GEOS-5 0.5 x 0.666 nested grid data.  Keep
!  variable names I_1x1, J_1x1, etc. unchanged for the present.
!                                                                             .
!  For now we will assume that we have saved tracer concentrations from 
!  a 2 x 2.5 window which overlays the corresponding NESTED WINDOW REGION.  
!  These 2 x 2.5 tracer concentrations are used as boundary conditions
!  for TPCORE.
! 
! !REVISION HISTORY: 
!  18 Dec 2009 - A. van Donkeelar - Initial version
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, II, J, JJ, JJ1
      REAL*8        :: X,         Y
      REAL*8,  SAVE :: XE2x25(145), YE2x25(92)

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
         DO I = 1, IIPAR

            ! Grid box longitude and latitude [degrees]
            X = GET_XMID( I, J, 1 )
            Y = GET_YMID( I, J, 1 )

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

      END FUNCTION GET_2x25_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_bc
!
! !DESCRIPTION: Subroutine ITS\_TIME\_FOR\_BC returns TRUE if it is time to 
!  read in the next set of boundary conditions for TPCORE, or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_BC() RESULT( FLAG )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_ELAPSED_MIN
!
! !RETURN VALUE:
!
      LOGICAL :: FLAG   ! =T if it's time to read BC's from disk
!
! !REVISION HISTORY: 
!  05 Mar 2003 - R. Yantosca - Initial version
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! ITS_TIME_FOR_BC begins here!
      !=================================================================
      FLAG = ( MOD( GET_ELAPSED_MIN(), TS_BC ) == 0 ) 
     
      END FUNCTION ITS_TIME_FOR_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tpcore_bc
!
! !DESCRIPTION: Subroutine INIT\_TPCORE\_BC initializes module variables and 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TPCORE_BC( TS, I0W, J0W, I1, 
     &                           J1, I2,  J2,  Input_Opt )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : TPBC_DIR
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET 
      USE GRID_MOD,           ONLY : ITS_A_NESTED_GRID
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: TS          ! Timestep for BC data
      INTEGER,        INTENT(IN) :: I0W, J0W    ! Transport region offsets
      INTEGER,        INTENT(IN) :: I1,  J1     ! Lon, lat indices @ LL corner 
      INTEGER,        INTENT(IN) :: I2,  J2     ! Lon, lat indices @ UR corner
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
! 
! !REVISION HISTORY: 
!  10 Feb 2003 - R. Yantosca - Initial version
!  (1 ) Now references N_TRACERS from "tracer_mod.f".  Now references LWINDO
!        from "logical_mod.f".  Now references TPBC_DIR from "directory_mod.f".
!        Now references ITS_A_NESTED_GRID from "grid_mod.f".  Also added 
!        arguments to take values from "input_mod.f". (bmy, 7/20/04)
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!  07 Sep 2012 - R. Yantosca - Now use 0d0 instead of 0e0 to zero BC arrays
!  17 Oct 2012 - M. Payer    - Bug fix: Define extent of coarse grid BC region
!                              even if LWINDO_CH or LWINDO_SE = False to
!                              avoid out-of-bounds errors in array BC
!  05 Jun 2013 - K. Yu       - Define BC boundaries for GEOS_57 0.25 NA grid
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Removed SEAC4RS C-preprocessor switch
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS, I, J

      ! For values from Input_Opt
      LOGICAL :: LWINDO,    LWINDO_NA,  LWINDO_CU
      LOGICAL :: LWINDO_EU, LWINDO_CH,  LWINDO2x25
      LOGICAL :: LWINDO_SE
      INTEGER :: N_TRACERS

      !=================================================================
      ! INIT_TPCORE_BC begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LWINDO     = Input_Opt%LWINDO
      LWINDO_NA  = Input_Opt%LWINDO_NA
      LWINDO_CU  = Input_Opt%LWINDO_CU
      LWINDO_EU  = Input_Opt%LWINDO_EU
      LWINDO_CH  = Input_Opt%LWINDO_CH
      LWINDO2x25 = Input_Opt%LWINDO2x25
      N_TRACERS  = Input_Opt%N_TRACERS

      ! Timestep for BC's [min]
      TS_BC = TS

      !-----------------------------------------------------------------
      !       ****** NESTED GRID TPCORE TRANSPORT REGION ******
      !-----------------------------------------------------------------

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

      !-----------------------------------------------------------------
      !                   ****** GLOBAL GRID ******
      !-----------------------------------------------------------------

      ! BC window region on the coarse grid
      I1_BC = I1                             ! Lon index, LL corner
      J1_BC = J1                             ! Lat index, LL corner
      I2_BC = I2                             ! Lon index, UR corner
      J2_BC = J2                             ! Lat index, UR corner

      ! Extent of coarse-grid BC REGION
      IM_BC = I2_BC - I1_BC + 1
      JM_BC = J2_BC - J1_BC + 1

      ! TPCORE internal transport window offset
      I0_BC = I1_BC - 1
      J0_BC = J1_BC - 1

      !-----------------------------------------------------------------
      !        ****** BC REGION FOR N. AMERICA NESTED GRID ******
      !-----------------------------------------------------------------
      IF ( .not. LWINDO2x25 ) THEN

#if   defined( GRID4x5 ) || defined( NESTED_NA )

#if   defined( GEOS_FP )
         !%%%%% 4x5 BC REGION FOR NORTH AMERICA NESTED GRID (GEOS_FP)
         I1_BC_NA = 11                       ! 4x5  lon index, LL corner
         J1_BC_NA = 26                       ! 4x5  lat index, LL corner
         I2_BC_NA = 25                       ! 4x5  lon index, UR corner
         J2_BC_NA = 39                       ! 4x5  lat index, UR corner
#else
         !%%%%% 4x5 BC REGION FOR NORTH AMERICA NESTED GRID %%%%%
         I1_BC_NA = 9                        ! 4x5  lon index, LL corner
         J1_BC_NA = 26                       ! 4x5  lat index, LL corner
         I2_BC_NA = 29                       ! 4x5  lon index, UR corner
         J2_BC_NA = 41                       ! 4x5  lat index, UR corner
#endif

#endif

      ELSE

#if   defined( GRID2x25 ) || defined( NESTED_NA )

#if   defined( GEOS_FP )
         !%%%%% 2x25 BC REGION FOR NORTH AMERICA NESTED GRID (GEOS_FP)
         I1_BC_NA = 21                       ! 2x25 lon index, LL corner
         J1_BC_NA = 51                       ! 2x25 lat index, LL corner
         I2_BC_NA = 49                       ! 2x25 lon index, UR corner
         J2_BC_NA = 76                       ! 2x25 lat index, UR corner
#else
         !%%%%% 2x25 BC REGION FOR NORTH AMERICA NESTED GRID %%%%%
         I1_BC_NA = 17                       ! 2x25 lon index, LL corner
         J1_BC_NA = 51                       ! 2x25 lat index, LL corner
         I2_BC_NA = 57                       ! 2x25 lon index, UR corner
         J2_BC_NA = 81                       ! 2x25 lat index, UR corner
#endif

#endif

      ENDIF
 
      ! Extent of coarse-grid NA BC REGION
      IM_BC_NA = I2_BC_NA - I1_BC_NA + 1
      JM_BC_NA = J2_BC_NA - J1_BC_NA + 1
 
      ! TPCORE NA internal transport window offset
      I0_BC_NA = I1_BC_NA - 1
      J0_BC_NA = J1_BC_NA - 1

      !-----------------------------------------------------------------
      !         ****** BC REGION FOR EUROPE NESTED GRID ******
      !-----------------------------------------------------------------
      IF ( .not. LWINDO2x25 ) THEN

#if   defined( GRID4x5 ) || defined( NESTED_EU )
         !%%%%% 4x5 BC REGION FOR EUROPE NESTED GRID %%%%%
         I1_BC_EU = 31                       ! 4x5  lon index, LL corner
         J1_BC_EU = 31                       ! 4x5  lat index, LL corner
         I2_BC_EU = 47                       ! 4x5  lon index, UR corner
         J2_BC_EU = 41                       ! 4x5  lat index, UR corner
#endif

      ELSE

#if   defined( GRID2x25 ) || defined( NESTED_EU )
         !%%%%% 2x25 BC REGION FOR EUROPE NESTED GRID %%%%%
         I1_BC_EU = 61                       ! 2x25 lon index, LL corner
         J1_BC_EU = 61                       ! 2x25 lat index, LL corner
         I2_BC_EU = 93                       ! 2x25 lon index, UR corner
         J2_BC_EU = 81                       ! 2x25 lat index, UR corner
#endif

      ENDIF
 
      ! Extent of coarse-grid EU BC REGION
      IM_BC_EU = I2_BC_EU - I1_BC_EU + 1
      JM_BC_EU = J2_BC_EU - J1_BC_EU + 1
 
      ! TPCORE NA internal transport window offset
      I0_BC_EU = I1_BC_EU - 1
      J0_BC_EU = J1_BC_EU - 1

      !-----------------------------------------------------------------
      !      ****** BC REGION FOR CHINA/SE ASIA NESTED GRID ******
      !-----------------------------------------------------------------
      IF (.not. LWINDO2x25 ) THEN 

#if   defined( GRID4x5 ) || defined( NESTED_CH ) 

#if   defined( GEOS_FP )

         !%%%%% 4x5 BC REGION FOR NESTED CHINA GRID (GEOS-FP only) %%%%%
         I1_BC_CH = 51                       ! 4x5 lon index, LL corner
         J1_BC_CH = 27                       ! 4x5 lat index, LL corner
         I2_BC_CH = 65                       ! 4x5 lon index, UR corner
         J2_BC_CH = 37                       ! 4x5 lat index, UR corner

#else

         !%%%%% 4x5 BC REGION FOR TRADITIONAL GEOS-5 CHINA NESTED GRID %%%%%
         I1_BC_CH = 51                       ! 4x5 lon index, LL corner
         J1_BC_CH = 21                       ! 4x5 lat index, LL corner
         I2_BC_CH = 67                       ! 4x5 lon index, UR corner
         J2_BC_CH = 37                       ! 4x5 lat index, UR corner

#endif

#endif

      ELSE

#if   defined( GRID2x25 ) || defined( NESTED_CH )

#if   defined( GEOS_FP )

         !%%%%% 2x25 BC REGION FOR NESTED CHINA GRID (GEOS-FP only) %%%%%
         I1_BC_CH = 101                      ! 2x25 lon index, LL corner
         J1_BC_CH = 53                       ! 2x25 lat index, LL corner
         I2_BC_CH = 129                      ! 2x25 lon index, UR corner
         J2_BC_CH = 74                       ! 2x25 lat index, UR corner

#else

         !%%%%% 2x25 BC REGION FOR TRADITIONAL GEOS-5 CHINA NESTED GRID %%%%%
         I1_BC_CH = 101                      ! 2x25 lon index, LL corner
         J1_BC_CH = 40                       ! 2x25 lat index, LL corner
         I2_BC_CH = 133                      ! 2x25 lon index, UR corner
         J2_BC_CH = 74                       ! 2x25 lat index, UR corner

#endif

#endif

      ENDIF
      
      ! Extent of coarse-grid CH BC REGION
      IM_BC_CH = I2_BC_CH - I1_BC_CH + 1
      JM_BC_CH = J2_BC_CH - J1_BC_CH + 1
 
      ! TPCORE CH internal transport window offset
      I0_BC_CH = I1_BC_CH - 1
      J0_BC_CH = J1_BC_CH - 1

      ! Extent of coarse-grid SEAC4RS BC REGION
      IM_BC_SE = I2_BC_SE - I1_BC_SE + 1
      JM_BC_SE = J2_BC_SE - J1_BC_SE + 1
 
      ! TPCORE SEAC4RS internal transport window offset
      I0_BC_SE = I1_BC_SE - 1
      J0_BC_SE = J1_BC_SE - 1

      ! Return if we are not saving 4x5 BC's
      ! or if it's not a nested grid simulation
      IF ( .not. ITS_A_NESTED_GRID() ) THEN
         IF ( .not. LWINDO ) RETURN
      ENDIF

      !=================================================================
      ! Allocate and initialize arrays 
      !=================================================================
      
      ! Boundary condition array for CUSTOM NESTED REGION
      IF ( LWINDO_CU ) THEN
         ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
         BC = 0d0
      ENDIF

      ! Boundary condition array for NESTED CHINA/SE ASIA REGION
      IF ( LWINDO_CH ) THEN
         ALLOCATE( BC_CH( IM_BC_CH, JM_BC_CH, 
     &                    LLPAR,    N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_CH' )
         BC_CH = 0d0
      ENDIF

      ! Boundary condition array for NESTED N. AMERICA REGION
      IF ( LWINDO_NA ) THEN
         ALLOCATE( BC_NA( IM_BC_NA, JM_BC_NA, 
     &                    LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_NA' )
         BC_NA = 0d0
      ENDIF

      ! Boundary condition array for NESTED EUROPE REGION
      IF ( LWINDO_EU ) THEN
         ALLOCATE( BC_EU( IM_BC_EU, JM_BC_EU, 
     &                    LLPAR,    N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_EU' )
         BC_EU = 0d0
      ENDIF

      ! If running a nested simulation, just allocate the custom BC
      ! array with nested dimensions
#if   defined( NESTED_CH )

      !%%%%% Define boundaries for NESTED CHINA GRID %%%%%
      I1_BC = I1_BC_CH
      J1_BC = J1_BC_CH
      I2_BC = I2_BC_CH
      J2_BC = J2_BC_CH
      IM_BC = IM_BC_CH
      JM_BC = JM_BC_CH
      I0_BC = I0_BC_CH
      J0_BC = J0_BC_CH

#elif defined( NESTED_NA )

      !%%%%% Define boundaries for NESTED NORTH AMERICA GRID %%%%%
      I1_BC = I1_BC_NA
      J1_BC = J1_BC_NA
      I2_BC = I2_BC_NA
      J2_BC = J2_BC_NA
      IM_BC = IM_BC_NA
      JM_BC = JM_BC_NA
      I0_BC = I0_BC_NA
      J0_BC = J0_BC_NA

#elif defined( NESTED_EU )

      !%%%%% Define boundaries for NESTED EUROPE GRID %%%%%
      I1_BC = I1_BC_EU
      J1_BC = J1_BC_EU
      I2_BC = I2_BC_EU
      J2_BC = J2_BC_EU
      IM_BC = IM_BC_EU
      JM_BC = JM_BC_EU
      I0_BC = I0_BC_EU
      J0_BC = J0_BC_EU

#endif

      ALLOCATE( BC( IM_BC, JM_BC, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
      BC = 0d0

      ALLOCATE( MAP1x1( IIPAR, JJPAR, 4 ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAP1x1' )
      MAP1x1 = 0

      END SUBROUTINE INIT_TPCORE_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tpcore_bc
!
! !DESCRIPTION: Subroutine CLEANUP\_TPCORE\_BC deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TPCORE_BC
! 
! !REVISION HISTORY: 
!  04 Mar 2003 - R. Yantosca - Initial version
!  15 May 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TPCORE_BC begins here!
      !=================================================================
      IF ( ALLOCATED( BC     ) ) DEALLOCATE( BC     )
      IF ( ALLOCATED( BC_NA  ) ) DEALLOCATE( BC_NA  )
      IF ( ALLOCATED( BC_EU  ) ) DEALLOCATE( BC_EU  )
      IF ( ALLOCATED( BC_CH  ) ) DEALLOCATE( BC_CH  )
      IF ( ALLOCATED( MAP1x1 ) ) DEALLOCATE( MAP1x1 )

      END SUBROUTINE CLEANUP_TPCORE_BC
!EOC
      END MODULE TPCORE_BC_MOD
