!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_mod
!
! !DESCRIPTION: Module INPUT\_MOD contains routines that read the GEOS-Chem 
!  input file at the start of the run and pass the information to several 
!  other GEOS-Chem F90 modules.
!\\
!\\
! !INTERFACE:
!
      MODULE INPUT_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: INITIALIZE_GEOS_GRID
      PUBLIC  :: READ_INPUT_FILE
      PUBLIC  :: GIGC_Init_Extra
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: READ_ONE_LINE
      PRIVATE :: SPLIT_ONE_LINE        
      PRIVATE :: READ_SIMULATION_MENU
      PRIVATE :: READ_TRACER_MENU  
      PRIVATE :: READ_AEROSOL_MENU     
      PRIVATE :: READ_EMISSIONS_MENU
      PRIVATE :: READ_FUTURE_MENU
      PRIVATE :: READ_CHEMISTRY_MENU
      PRIVATE :: READ_TRANSPORT_MENU
      PRIVATE :: READ_CONVECTION_MENU
      PRIVATE :: READ_DEPOSITION_MENU
      PRIVATE :: READ_OUTPUT_MENU
      PRIVATE :: READ_DIAGNOSTIC_MENU
      PRIVATE :: SET_TINDEX
      PRIVATE :: READ_ND49_MENU      
      PRIVATE :: READ_ND50_MENU  
      PRIVATE :: READ_ND51_MENU  
      PRIVATE :: READ_ND51b_MENU  
      PRIVATE :: READ_ND63_MENU
      PRIVATE :: READ_PROD_LOSS_MENU 
      PRIVATE :: READ_UNIX_CMDS_MENU
      PRIVATE :: READ_NESTED_GRID_MENU
      PRIVATE :: READ_ARCHIVED_OH_MENU
      PRIVATE :: READ_O3PL_MENU
      PRIVATE :: READ_BENCHMARK_MENU  
      PRIVATE :: READ_CH4_MENU
      PRIVATE :: VALIDATE_DIRECTORIES  
      PRIVATE :: CHECK_DIRECTORY
      PRIVATE :: CHECK_TIME_STEPS 
      PRIVATE :: IS_LAST_DAY_GOOD
      PRIVATE :: INIT_INPUT
#if defined( TOMAS )
      PRIVATE :: INIT_TOMAS_MICROPHYSICS
#endif
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now references LSOA in READ_AEROSOL_MENU (bmy, 9/28/04)
!  (2 ) Fixed error checks and assign LSPLIT for tagged Hg.  Also now 
!        refernces LAVHRRLAI from "logical_mod.f" (eck, bmy, 12/20/04)
!  (3 ) Updated for crystalline/aqueous aerosol tracers.  Also moved routine
!        IS_LAST_DAY_GOOD here from "main.f".  Also now references 
!        "ocean_mercury_mod.f".  Also now open the bpch file for output in
!        READ_DIAGNOSTIC_MENU instead of in "main.f".  (cas, sas, bmy, 2/3/05)
!  (4 ) Now references "diag03_mod.f" and "diag41_mod.f".  Fixed minor
!        bugs.  Now references FILE_EXISTS from "file_mod.f".  Updated
!        comments. (bmy, 3/28/05)
!  (5 ) Now modified for GEOS-5 and GCAP met fields.  Also now set LSPLIT
!        correctly for HCN/CH3CN simulation. (swu, xyp, bmy, 6/30/05)
!  (6 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (7 ) Now read LMEGAN switch for MEGAN biogenics.  Now read variable
!        DATA_DIR_1x1 for 1x1 emissions files, etc.  Now reference XNUMOL and
!        XNUMOLAIR from "tracer_mod.f" (tmf, bmy, 10/25/05) 
!  (8 ) Now read LEMEP switch for EMEP emissions (bdf, bmy, 11/1/05)
!  (9 ) Now added MERCURY MENU section.  Also fixed bug in READ_ND48_MENU.
!        (eck, cdh, bmy, 3/6/06)
!  (10) Now read LGFED2BB switch for GFED2 biomass emissions (bmy, 4/5/06)
!  (11) Bug fix for GCAP in IS_LAST_DAY_GOOD.  Also now read LCTH, LMFLUX,
!        LPRECON in READ_EMISSIONS_MENU. (bmy, 5/10/06)
!  (12) Updated for ND42 SOA concentration diagnostic (dkh, bmy, 5/22/06)
!  (13) Modified for future emissions (swu, bmy, 6/1/06)
!  (14) Modified for BRAVO emissions (rjp, kfb, bmy, 6/26/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields.  Also modified 
!        for David Streets' emissions. (bmy, 8/17/06)
!  (16) Modified for variable tropopause.  Also set dimension of ND28 diag
!        for GFED2 or default biomass burning.  Now read if Time Spent in 
!        Troposphere is wanted (phs, bmy, 10/17/06)
!  (17) Now modified for OTD-LIS local redistribution.  Remove references
!        to GEOS-1 and GEOS-STRAT run dirs. (bmy, 11/5/07)
!  (18) New error traps for OTD-LIS scaling, dependent on met field type.
!        Bug fix, create string variables for ERROR_STOP.  Bug fix: use ND52
!        in call to SET_TINDEX in READ_DIAGNOSTIC_MENU. (ltm, bmy, 2/11/08)
!  (19) Bug fix: use (0,0) in call to INIT_TRANSFER (phs, 6/17/08)
!  (20) Minor fix in READ_TRANSPORT_MENU (cdh, bmy, 7/7/08)
!  (21) Fixed typo READ_EMISSIONS_MENU for GEOS-3 (bmy, 10/30/08)
!  (22) Set upper limit on dynamic timestep for 0.5 x 0.666 nested
!        grids (yxw, bmy, dan, 11/6/08)
!  (23) Now read LCAC switch for CAC emissions (amv, 1/09/2008)
!  (24) Move the call to NDXX_SETUP (phs, 11/18/08)
!  (25) Minor bug fix in READ_DIAGNOSTIC_MENU (tmf, 2/10/09)
!  (26) Add LMEGANMONO switch in emission menu (ccc, 3/2/09)
!  (27) Add LDICARB switch in aerosol menu (ccc, tmf, 3/10/09)
!  (28) Now read LCOOKE in aerosol menu (phs, 5/18/09)
!  (29) Add CH4_MENU in input.geos (kjw, 8/18/09)
!  (30) Corrected typos in CHECK_TIME_STEPS (bmy, 8/21/09)
!  (31) Now read LLINOZ in READ_SIMULATION_MENU (dbm, bmy, 10/16/09)
!  (32) Remove reference to obsolete embedded chemistry stuff (bmy, 2/25/10)
!  (33) Remove depreciated lightning options (ltm, bmy, 1/24/11)
!  25 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  29 Jul 2011 - R. Yantosca - Bug fix in READ_EMISSIONS_MENU for nested NA
!  07 Sep 2011 - P. Kasibhatla - Modified to include monthly GFED3
!  17 Jan 2012 - P. Kasibhatla - Modified to include daily and 3-hourly GFED3
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met data
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  30 Jul 2012 - R. Yantosca - READ_INPUT_FILE now accepts am_I_Root from
!                              both the ESMF interface and main.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Aug 2012 - R. Yantosca - Now make IU_GEOS a global module variable
!                              so that we can define it with findFreeLun
!  02 Nov 2012 - R. Yantosca - Now pass the Input Options object to routines;
!                              this will eventually replace logical_mod, etc.
!  26 Feb 2013 - M. Long     - Now make INITIALIZE_GEOS_GRID a public routine
!  04 Mar 2013 - R. Yantosca - Add routine GIGC_Init_Extra to move some init
!                              calls out of the run stage when using ESMF
!  23 Apr 2013 - R. Yantosca - For TOMAS, rename READ_MICROPHYSICS_MENU to 
!                              INIT_TOMAS_MICROPHYSICS
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  17 Sep 2013 - R. Yantosca - Increase MAXDIM from 255 to 500 for more tracers
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP everywhere
!  26 Sep 2013 - R. Yantosca - Now read GEOS_FP_DIR from Input_Opt everywhere
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      LOGICAL            :: VERBOSE  = .FALSE.
      INTEGER            :: IU_GEOS
      INTEGER, PARAMETER :: FIRSTCOL = 26
      INTEGER, PARAMETER :: MAXDIM   = 500
      INTEGER            :: TS_CHEM
      INTEGER            :: TS_DYN
      INTEGER            :: TS_CONV
      INTEGER            :: TS_EMIS
      INTEGER            :: TS_UNIT
      INTEGER            :: CT1, CT2, CT3
      CHARACTER(LEN=255) :: FILENAME = 'input.geos'
      CHARACTER(LEN=255) :: TOPTITLE
      CHARACTER(LEN=255) :: BPCH_FILE
      CHARACTER(LEN=255) :: DIAGINFO  
      CHARACTER(LEN=255) :: TRACERINFO

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_file
!
! !DESCRIPTION: Subroutine READ\_INPUT\_FILE is the driver program for 
!  reading the GEOS-Chem input file "input.geos" from disk. 
!\\
!\\
! !INTERFACE:
! 
      SUBROUTINE READ_INPUT_FILE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CHARPAK_MOD,        ONLY : STRREPL
      USE FILE_MOD,           ONLY : IOERROR
      USE GAMAP_MOD,          ONLY : DO_GAMAP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE SEASALT_MOD,        ONLY : INIT_SEASALT
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now call DO_GAMAP from "gamap_mod.f" to create "diaginfo.dat" and
!        "tracerinfo.dat" files after all diagnostic menus have been read in
!  (2 ) Now call NDXX_setup from this routine (phs, 11/18/08)
!  (3 ) Now call READ_ND51b_MENU (amv, bmy, 12/18/09)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root from main.F, so that we
!                              can get rid of duplicate code in DEVEL blocks
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  03 Aug 2012 - R. Yantosca - Now make IU_GEOS a global module variable
!                              so that we can define it with findFreeLun
!  15 Oct 2012 - R. Yantosca - Add EXTERNAL_GRID, XTERNAL_FORCING to #ifdef
!  16 Oct 2012 - R. Yantosca - Don't call CHECK_TIME_STEPS if we are calling
!                              READ_INPUT_FILE from the ESMF interface
!  09 Nov 2012 - R. Yantosca - Now pass Input_Opt to lower-level routines
!  06 Dec 2012 - R. Yantosca - Now call CHECK_TIME_STEPS when we are connecting
!                              to the GEOS-5 GCM via the ESMF environment,
!  19 Mar 2013 - R. Yantosca - When using ESMF interface to GEOS-5, append
!                              ".rc" to input.geos (instead of ___.rc)
!  04 Apr 2013 - R. Yantosca - Now pass objects to DO_GAMAP routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: EOF
      INTEGER            :: IOS
      CHARACTER(LEN=1)   :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)   :: SPACE = ' '
      CHARACTER(LEN=MAXDIM) :: LINE

      !=================================================================
      ! READ_INPUT_FILE begins here!
      !=================================================================  

      ! Echo output
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)' ) 'G E O S - C H E M   U S E R   I N P U T'
         WRITE( 6, 100   ) TRIM( FILENAME )
 100     FORMAT( 'READ_INPUT_FILE: Reading ', a )
      ENDIF

      ! Initialize directory & logical variables
      CALL INIT_INPUT

      ! Find a free file LUN
      IU_GEOS  = findFreeLUN()

      ! Assume success
      RC       = GIGC_SUCCESS

#if defined( ESMF_ )
      ! We need to add the ".rc" extension to the "input.geos" file
      ! so that GEOS-5 will copy it to the proper run directory.
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Open file
      OPEN( IU_GEOS, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:1' )

      ! Read TOPTITLE for binary punch file
      TOPTITLE = READ_ONE_LINE( EOF  )
      IF ( EOF ) RETURN

      ! Loop until EOF
      DO 
         
         ! Read a line from the file, exit if EOF
         LINE = READ_ONE_LINE( EOF ) 
         IF ( EOF ) EXIT
         
         ! Replace tab characters in LINE (if any) w/ spaces
         CALL STRREPL( LINE, TAB, SPACE )

         !=============================================================
         ! Call individual subroutines to read sections of the file
         ! 
         ! NOTE: You are pretty flexible in setting the order of the
         ! menus in the input file; however, a few guidelines apply:
         !
         ! (1) SIMULATION MENU should be listed first.
         ! (2) TRACER MENU should be listed second.
         ! (3) EMISSIONS, AEROSOL, CHEMISTRY, TRANSPORT, CONVECTION, 
         !      and DEPOSITION menus (in any order) should follow.
         ! (4) Diagnostic menus, including OUTPUT, DIAGNOSTIC,
         !      PLANEFLIGHT, ND48, ND49, ND50, ND51, and PROD_LOSS
         !      menus (in any order) should follow next.
         ! (5) The following menus have no other restriction and
         !      can be placed anywhere (but by convention we will
         !      place them after the diagnostic menu): NESTED GRID
         !      UNIX CMDS, ARCHIVED OH, and O3PL menus.
         !=============================================================
         IF      ( INDEX( LINE, 'SIMULATION MENU'  ) > 0 ) THEN
            CALL READ_SIMULATION_MENU( am_I_Root, Input_Opt, RC )             

         ELSE IF ( INDEX( LINE, 'TRACER MENU'      ) > 0 ) THEN
            CALL READ_TRACER_MENU( am_I_Root, Input_Opt, RC )                 

         ELSE IF ( INDEX( LINE, 'AEROSOL MENU'     ) > 0 ) THEN
            CALL READ_AEROSOL_MENU( am_I_Root, Input_Opt, RC )              

#if defined( TOMAS )
         ! Initialize the TOMAS microphysics package, if necessary
         CALL INIT_TOMAS_MICROPHYSICS( am_I_Root, Input_Opt, RC )    
#endif

#if defined( APM )
         ! Initialize APM microphysics package, if necessary
         ELSE IF ( INDEX( LINE, 'APM MENU'         ) > 0 ) THEN
            CALL READ_APM_MENU( am_I_Root, Input_Opt, RC )
#endif

         ELSE IF ( INDEX( LINE, 'EMISSIONS MENU'   ) > 0 ) THEN
            CALL READ_EMISSIONS_MENU( am_I_Root, Input_Opt, RC )             

         ELSE IF ( INDEX( LINE, 'FUTURE MENU'      ) > 0 ) THEN
            CALL READ_FUTURE_MENU( am_I_Root, Input_Opt, RC ) 
                                                   
         ELSE IF ( INDEX( LINE, 'CHEMISTRY MENU'   ) > 0 ) THEN
            CALL READ_CHEMISTRY_MENU( am_I_Root, Input_Opt, RC )              
                                                 
         ELSE IF ( INDEX( LINE, 'TRANSPORT MENU'   ) > 0 ) THEN
            CALL READ_TRANSPORT_MENU( am_I_Root, Input_Opt, RC )              
                                             
         ELSE IF ( INDEX( LINE, 'CONVECTION MENU'  ) > 0 ) THEN
            CALL READ_CONVECTION_MENU( am_I_Root, Input_Opt, RC )            
                                             
         ELSE IF ( INDEX( LINE, 'DEPOSITION MENU'  ) > 0 ) THEN
            CALL READ_DEPOSITION_MENU( am_I_Root, Input_Opt, RC )             

         ELSE IF ( INDEX( LINE, 'CO2 SIM MENU'     ) > 0 ) THEN
            CALL READ_CO2_SIM_MENU( am_I_Root, Input_Opt, RC )
                                    
         ELSE IF ( INDEX( LINE, 'GAMAP MENU'       ) > 0 ) THEN
            CALL READ_GAMAP_MENU( am_I_Root, Input_Opt, RC )                 
                                                  
         ELSE IF ( INDEX( LINE, 'OUTPUT MENU'      ) > 0 ) THEN
            CALL READ_OUTPUT_MENU( am_I_Root, Input_Opt, RC )                 
                                                  
         ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU'  ) > 0 ) THEN
            CALL READ_DIAGNOSTIC_MENU( am_I_Root, Input_Opt, RC )             

         ELSE IF ( INDEX( LINE, 'PLANEFLIGHT MENU' ) > 0 ) THEN
            CALL READ_PLANEFLIGHT_MENU( am_I_Root, Input_Opt, RC )             
                                                  
         ELSE IF ( INDEX( LINE, 'ND48 MENU'        ) > 0 ) THEN
            CALL READ_ND48_MENU( am_I_Root, Input_Opt, RC )                  

         ELSE IF ( INDEX( LINE, 'ND49 MENU'        ) > 0 ) THEN
            CALL READ_ND49_MENU( am_I_Root, Input_Opt, RC )                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND50 MENU'        ) > 0 ) THEN
            CALL READ_ND50_MENU( am_I_Root, Input_Opt, RC )                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND51 MENU'        ) > 0 ) THEN
            CALL READ_ND51_MENU( am_I_Root, Input_Opt, RC )                   
                                                  
         ELSE IF ( INDEX( LINE, 'ND51b MENU'       ) > 0 ) THEN
            CALL READ_ND51b_MENU( am_I_Root, Input_Opt, RC )  

         ELSE IF ( INDEX( LINE, 'ND63 MENU'        ) > 0 ) THEN
            CALL READ_ND63_MENU( am_I_Root, Input_Opt, RC )                   

         ELSE IF ( INDEX( LINE, 'PROD & LOSS MENU' ) > 0 ) THEN
            CALL READ_PROD_LOSS_MENU( am_I_Root, Input_Opt, RC )               
                                             
         ELSE IF ( INDEX( LINE, 'UNIX CMDS MENU'   ) > 0 ) THEN 
            CALL READ_UNIX_CMDS_MENU( am_I_Root, Input_Opt, RC )              

         ELSE IF ( INDEX( LINE, 'NESTED GRID MENU' ) > 0 ) THEN 
            CALL READ_NESTED_GRID_MENU( am_I_Root, Input_Opt, RC )

         ELSE IF ( INDEX( LINE, 'ARCHIVED OH MENU' ) > 0 ) THEN 
            CALL READ_ARCHIVED_OH_MENU( am_I_Root, Input_Opt, RC  )

         ELSE IF ( INDEX( LINE, 'O3 P/L MENU'      ) > 0 ) THEN 
            CALL READ_O3PL_MENU( am_I_Root, Input_Opt, RC  )

         ELSE IF ( INDEX( LINE, 'BENCHMARK MENU'   ) > 0 ) THEN 
            CALL READ_BENCHMARK_MENU( am_I_Root, Input_Opt, RC  )
         
         ELSE IF ( INDEX( LINE, 'POPS MENU'        ) > 0 ) THEN
            CALL READ_POPS_MENU( am_I_Root, Input_Opt, RC )
           
         ELSE IF ( INDEX( LINE, 'MERCURY MENU'     ) > 0 ) THEN 
            CALL READ_MERCURY_MENU( am_I_Root, Input_Opt, RC  )
                                   
         ELSE IF ( INDEX( LINE, 'CH4 MENU'         ) > 0 ) THEN 
            CALL READ_CH4_MENU( am_I_Root, Input_Opt, RC )
         ELSE IF ( INDEX( LINE, 'END OF FILE'      ) > 0 ) THEN 
            EXIT

         ENDIF  
      ENDDO

      ! Close input file
      CLOSE( IU_GEOS )

      !=================================================================
      ! Further error-checking and initialization
      !=================================================================

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do not initialize diagnostics when we are connecting to
      ! the GEOS-5 GCM via ESMF. (bmy, 11/2/12)
      !-----------------------------------------------------------------
#else

      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! We only need to initialize GEOS-5 diagnostics if we are not
      ! using ESMF.  The GEOS-5 GCM has its own History component
      ! for diagnostic output (bmy, 3/4/13)
      !-----------------------------------------------------------------

      ! Define GEOS-Chem diagnostic info
      CALL NDXX_SETUP( am_I_Root, Input_Opt, RC )

      ! Make sure all directories are valid
      CALL VALIDATE_DIRECTORIES

#endif

      ! Check GEOS-CHEM timesteps
      CALL CHECK_TIME_STEPS( am_I_Root )

      ! Create "diaginfo.dat" and "tracerinfo.dat" files for GAMAP
      CALL DO_GAMAP( am_I_Root, Input_Opt, RC )

      ! Echo output
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF

      END SUBROUTINE READ_INPUT_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_one_line
!
! !DESCRIPTION: Subroutine READ\_ONE\_LINE reads a line from the input file.  
!  If the global variable VERBOSE is set, the line will be printed to stdout.  
!  READ\_ONE\_LINE can trap an unexpected EOF if LOCATION is passed.  
!  Otherwise, it will pass a logical flag back to the calling routine, 
!  where the error trapping will be done.
!\\
!\\
! !INTERFACE:
!
      FUNCTION READ_ONE_LINE( EOF, LOCATION ) RESULT( LINE )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION    ! Msg to display
!
! !OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(OUT)          :: EOF         ! Denotes EOF 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Now make IU_GEOS a global module variable
!                              so that we can define it with findFreeLun
!  17 Sep 2013 - R. Yantosca - Extend line length to read in more tracers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
!------------------------------------------------------------------------------
! Prior to 9/17/13:
! Need to extend the line length for many more tracers (bmy, 9/17/13)
!      CHARACTER(LEN=255) :: LINE, MSG
!------------------------------------------------------------------------------
      CHARACTER(LEN=500) :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_one_line:1' )

      ! Print the line (if necessary)
      IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      END FUNCTION READ_ONE_LINE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: split_one_line
!
! !DESCRIPTION: Subroutine SPLIT\_ONE\_LINE reads a line from the input file 
!  (via routine READ\_ONE\_LINE), and separates it into substrings.
!\\
!\\
!  SPLIT\_ONE\_LINE also checks to see if the number of substrings found is 
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori, 
!  you can skip the error check.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION ) 
!
! !USES:
!
      USE CHARPAK_MOD, ONLY: STRSPLIT
!
! !INPUT PARAMETERS: 
!
      ! Number of substrings we expect to find
      INTEGER,            INTENT(IN)  :: N_EXP

      ! Name of routine that called SPLIT_ONE_LINE
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION 
!
! !OUTPUT PARAMETERS:
!
      ! Array of substrings (separated by " ")
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM)

      ! Number of substrings actually found
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  17 Sep 2013 - R. Yantosca - Extend LINE to 500 chars to allow more tracers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                         :: EOF
!----------------------------------------------------------------------------
! Prior to 9/17/13:
! Extend LINE to 500 chars to allow more tracers (bmy, 9/17/13)
!      CHARACTER(LEN=255)              :: LINE, MSG
!----------------------------------------------------------------------------
      CHARACTER(LEN=MAXDIM)           :: LINE
      CHARACTER(LEN=255)              :: MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!' 
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong 
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF
       
      END SUBROUTINE SPLIT_ONE_LINE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_simulation_menu
!
! !DESCRIPTION: Subroutine READ\_SIMULATION\_MENU reads the SIMULATION MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_SIMULATION_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!     
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR,    DATA_DIR_1x1, GCAP_DIR
      USE DIRECTORY_MOD,      ONLY : GEOS_4_DIR,  GEOS_5_DIR
      USE DIRECTORY_MOD,      ONLY : MERRA_DIR,   GEOS_FP_DIR
      USE DIRECTORY_MOD,      ONLY : RUN_DIR
      USE DIRECTORY_MOD,      ONLY : TEMP_DIR   
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : SET_XOFFSET, SET_YOFFSET
      USE LOGICAL_MOD,        ONLY : LSVGLB,      LUNZIP,       LWAIT
      USE LOGICAL_MOD,        ONLY : LVARTROP
      USE RESTART_MOD,        ONLY : SET_RESTART
      USE TIME_MOD,           ONLY : SET_BEGIN_TIME,   SET_END_TIME 
      USE TIME_MOD,           ONLY : SET_CURRENT_TIME, SET_DIAGb
      USE TIME_MOD,           ONLY : SET_NDIAGTIME,    GET_TAU
      USE TRANSFER_MOD,       ONLY : INIT_TRANSFER
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      USE TIME_MOD,           ONLY : Accept_External_Date_Time
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Bug fix: Read LSVGLB w/ the * format and not w/ '(a)'. (bmy, 2/23/05)
!  (2 ) Now read GEOS_5_DIR and GCAP_DIR from input.geos (swu, bmy, 5/25/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now references DATA_DIR_1x1 for 1x1 emissions files (bmy, 10/24/05)
!  (5 ) Now read switch for using variable tropopause or not (phs, 9/14/06)
!  (6 ) Remove references to GEOS-1 and GEOS-STRAT run dirs.  Now calls 
!        INIT_TRANSFER (bmy, 11/5/07)
!  (7 ) Fix typo in "print to screen" section  (phs, 6/1/08)
!  (8 ) Call INIT_TRANSFER w/ (0,0) instead of (I0,J0) (phs, 6/17/08)
!  (10) Now read LLINOZ switch from input.geos file (dbm, bmy, 10/16/09)
!  13 Aug 2010 - R. Yantosca - Now read MERRA_DIR
!  19 Aug 2010 - R. Yantosca - Set LUNZIP=F for MERRA met fields.
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Feb 2012 - R. Yantosca - Now read GEOS_57_DIR for GEOS-5.7.x met
!  08 Feb 2012 - R. Yantosca - Set LUNZIP=F for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now call routine INITIALIZE_GEOS_GRID to
!                              initialize horizontal grid parameters
!  10 Jun 2012 - L. Murray   - Move Linoz to chemistry menu
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  06 Dec 2012 - R. Yantosca - Now get NYMDb, NHMSb, NYMDe, NHMSe from the ESMF
!                              environment when connecting to the GEOS-5 GCM
!  11 Dec 2012 - R. Yantosca - ACCEPT_DATE_TIME_FROM_ESMF has now been renamed
!                              to ACCEPT_EXTERNAL_DATE_TIME
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I0,    J0
      INTEGER            :: N,     NDIAGTIME
      INTEGER            :: NYMDb, NHMSb 
      INTEGER            :: NYMDe, NHMSe
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      CHARACTER(LEN=255) :: IN_RST_FILE
      CHARACTER(LEN=255) :: OUT_RST_FILE

      !=================================================================
      ! READ_SIMULATION_MENU begins here!
      !=================================================================

      ! Start YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_simulation_menu:1' )
      READ( SUBSTRS(1:N), * ) NYMDb, NHMSb

      ! End YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_simulation_menu:2' )
      READ( SUBSTRS(1:N), * ) NYMDe, NHMSe

      ! Run directory
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) RUN_DIR

      ! Input restart file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) IN_RST_FILE

      ! Make new restart file?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:5' )
      READ( SUBSTRS(1:N), * ) LSVGLB

      ! Output restart file(s)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) OUT_RST_FILE

      ! Root data dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:7' )
      READ( SUBSTRS(1:N), '(a)' ) DATA_DIR

      ! GCAP subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:8' )
      READ( SUBSTRS(1:N), '(a)' ) GCAP_DIR

      ! GEOS-4 subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:9' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_4_DIR

      ! GEOS-5 subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:10' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_5_DIR

      ! GEOS-5.7.x subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:11' )
      READ( SUBSTRS(1:N), '(a)' ) GEOS_FP_DIR

      ! MERRA subdir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:12' )
      READ( SUBSTRS(1:N), '(a)' ) MERRA_DIR

      ! 1x1  dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:13' )
      READ( SUBSTRS(1:N), '(a)' ) DATA_DIR_1x1

      ! Temp dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:14' )
      READ( SUBSTRS(1:N), '(a)' ) TEMP_DIR

      ! Unzip met fields
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:15' )
      READ( SUBSTRS(1:N), *     ) LUNZIP

      ! Wait for met fields?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:16' )
      READ( SUBSTRS(1:N), *     ) LWAIT

      ! Variable Tropopause
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:17' )
      READ( SUBSTRS(1:N), *     ) LVARTROP

      ! I0, J0
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_simulation_menu:18' )
      READ( SUBSTRS(1:N), *     ) I0, J0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_simulation_menu:19' )

      !=================================================================
      ! Add safety checks for logical switches
      !=================================================================
#if   defined( MERRA ) || defined( GEOS_FP )
      ! Turn unzipping off for MERRA met fields
      LUNZIP = .FALSE. 
#endif

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
#if defined( ESMF_ ) 
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! If we are connecting to the ESMF interface, we need to take
      ! the start & end dates as defined in the ESMF resource file.
      ! (i.e. GEOSCHEMchem_GridComp_mod.rc) instead of those in 
      ! input.geos.  This is because the ESMF Clock object needs to be
      ! defined at the highest level (in the driver routine), before 
      ! input.geos is ever read.
      !
      ! Therefore, we will assign the start & end date fields (i.e.
      ! Input_Opt%NYMDb, Input_Opt%NYMDe, Input_Opt%NHMSb, and 
      ! Input_Opt%NHMSe) in the Gridded Component module file
      ! GEOSCHEMchem_GridComp_Mod.F90 (i.e. two levels higher
      ! in the code).  We don'need to define those fields here, so
      ! we have bracketed this with an #ifdef. 
      !
      ! Ideallly it would be good to have one single version of 
      ! "input.geos" that can double as the GEOS-Chem input file as 
      ! well as the ESMF resource.  We may be a way off from that.
      ! One potential solution is to use a Perl script to generate
      ! both the input.geos file and the ESMF resource file, so that
      ! both files are consistent with each other. (bmy, 12/6/12)
      !-----------------------------------------------------------------
      
      ! Save to local variables for printing below
      NYMDb                   = Input_Opt%NYMDb
      NHMSb                   = Input_Opt%NHMSb
      NYMDe                   = Input_Opt%NYMDe
      NHMSe                   = Input_Opt%NHMSe
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the standard GEOS-Chem is to initialize 
      ! the fields Input_Opt%NYMDb, Input_Opt%NYMDe, Input_Opt%NHMSb, 
      ! and Input_Opt%NHMSe from the values read in from input.geos.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------

      Input_Opt%NYMDb         = NYMDb 
      Input_Opt%NHMSb         = NHMSb
      Input_Opt%NYMDe         = NYMDe
      Input_Opt%NHMSe         = NHMSe
#endif
      Input_Opt%RUN_DIR       = RUN_DIR
      Input_Opt%IN_RST_FILE   = IN_RST_FILE
      Input_Opt%LSVGLB        = LSVGLB
      Input_Opt%OUT_RST_FILE  = OUT_RST_FILE
      Input_Opt%DATA_DIR      = DATA_DIR
      Input_Opt%GCAP_DIR      = GCAP_DIR
      Input_Opt%GEOS_4_DIR    = GEOS_4_DIR
      Input_Opt%GEOS_5_DIR    = GEOS_5_DIR
      Input_Opt%GEOS_FP_DIR   = GEOS_FP_DIR
      Input_Opt%MERRA_DIR     = MERRA_DIR
      Input_Opt%DATA_DIR_1x1  = DATA_DIR_1x1
      Input_Opt%TEMP_DIR      = TEMP_DIR
      Input_Opt%LUNZIP        = LUNZIP
      Input_Opt%LWAIT         = LWAIT
      Input_Opt%LVARTROP      = LVARTROP
      Input_Opt%NESTED_I0     = I0
      Input_Opt%NESTED_J0     = J0 

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'SIMULATION MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, 100 ) 'Start time of run           : ', NYMDb, NHMSb
         WRITE( 6, 100 ) 'End time of run             : ', NYMDe, NHMSe
         WRITE( 6, 110 ) 'Run directory               : ',
     &                    TRIM( RUN_DIR )
         WRITE( 6, 110 ) 'Data Directory              : ',
     &                    TRIM( DATA_DIR )
         WRITE( 6, 110 ) 'GCAP       sub-directory    : ', 
     &                    TRIM( GCAP_DIR )
         WRITE( 6, 110 ) 'GEOS-4     sub-directory    : ', 
     &                    TRIM( GEOS_4_DIR )
         WRITE( 6, 110 ) 'GEOS-5     sub-directory    : ', 
     &                    TRIM( GEOS_5_DIR )
         WRITE( 6, 110 ) 'GEOS-5.7.x sub-directory    : ', 
     &                    TRIM( GEOS_FP_DIR )
         WRITE( 6, 110 ) 'MERRA      sub-directory    : ', 
     &                    TRIM( MERRA_DIR )
         WRITE( 6, 110 ) '1x1 Emissions etc Data Dir  : ',
     &                    TRIM( DATA_DIR_1x1 )
         WRITE( 6, 110 ) 'Temporary Directory         : ', 
     &                    TRIM( TEMP_DIR )
         WRITE( 6, 110 ) 'Input restart file          : ', 
     &                    TRIM( IN_RST_FILE )
         WRITE( 6, 120 ) 'Create restart file?        : ', LSVGLB
         WRITE( 6, 110 ) 'Output restart file(s)      : ', 
     &                    TRIM( OUT_RST_FILE )
         WRITE( 6, 120 ) 'Unzip met fields?           : ', LUNZIP
         WRITE( 6, 120 ) 'Wait for met fields?        : ', LWAIT
         WRITE( 6, 120 ) 'Use variable tropopause?    : ', LVARTROP
         WRITE( 6, 130 ) 'Global offsets I0, J0       : ', I0, J0
      ENDIF

      ! Format statements
 100  FORMAT( A, I8.8, 1X, I6.6 )
 110  FORMAT( A, A              )
 120  FORMAT( A, L5             )
 130  FORMAT( A, 2I5            )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! If we are connecting to the GEOS-5 GCM via ESMF, then pass the
      ! values for the start & end times of the simulation directly to
      ! GeosUtil/time_mod.F via subroutine ACCEPT_EXTERNAL_DATE_TIME.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------
      CALL Accept_External_Date_Time( am_I_Root   = am_I_Root, 
     &                                value_NYMDb = Input_Opt%NYMDb,
     &                                value_NHMSb = Input_Opt%NHMSb,
     &                                value_NYMDe = Input_Opt%NYMDe,
     &                                value_NHMSe = Input_Opt%NHMSe,
     &                                value_NYMD  = Input_Opt%NYMDb,
     &                                value_NHMS  = Input_Opt%NHMSb,
     &                                RC          = RC               )
     &                                 

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! If we are not using ESMF, then call the traditional GEOS-Chem
      ! timing routines (from GeosUtil/time_mod.F) to set the start &
      ! end times of the simulation, as well as the current time.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------

      ! Set start time of run in "time_mod.f"
      CALL SET_BEGIN_TIME( NYMDb, NHMSb )

      ! Set end time of run in "time_mod.f"
      CALL SET_END_TIME( NYMDe, NHMSe )

      ! Set the current time
      CALL SET_CURRENT_TIME()

#endif

      ! Set the time of day for writing bpch files
      NDIAGTIME = NHMSe !### test
      CALL SET_NDIAGTIME( NDIAGTIME )

      ! Set the start of the 1st diagnostic interval
      CALL SET_DIAGb( GET_TAU() )
     
      ! Set input & output restart file names
      CALL SET_RESTART( IN_RST_FILE, OUT_RST_FILE )

      ! Set global offsets for the horizontal grid
      CALL SET_XOFFSET( I0 )
      CALL SET_YOFFSET( J0 )

      ! Initialize GEOS-Chem horizontal grid structure
      CALL INITIALIZE_GEOS_GRID( am_I_Root, RC )

      ! Initialze quantities for "transfer_mod.f"
      CALL INIT_TRANSFER( 0, 0 )

      ! Set counter
      CT1 = CT1 + 1

      END SUBROUTINE READ_SIMULATION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initialize_geos_grid
!
! !DESCRIPTION: Subroutine INITIALIZE\_GEOS\_GRID calls routines from
!  grid\_mod.F90 to initialize the horizontal grid parameters.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INITIALIZE_GEOS_GRID( am_I_Root, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GLOBAL_GRID_MOD, ONLY : COMPUTE_GLOBAL_GRID
      USE GRID_MOD,        ONLY : COMPUTE_GRID
      USE GRID_MOD,        ONLY : INIT_GRID
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root    ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  The module grid_mod.F90 has been modified to save grid parameters in 3D
!  format, which will facilitate interfacing GEOS-Chem to a GCM.
!                                                                             .
!  The module global_grid_mod.F90 contains several of the global grid arrays
!  (*_g) originally in grid_mod.F. These arrays are used in regridding GFED3
!  biomass emissions, which are available on a 0.5x0.5 global grid. The global
!  arrays may need to be used in the future for regridding other emissions for
!  nested grids.
!                                                                             .
! !REVISION HISTORY: 
!  01 Mar 2012 - R. Yantosca - Initial version
!  01 May 2012 - M. Payer    - Add call to COMPUTE_GLOBAL_GRID for nested grids
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  30 Nov 2012 - R. Yantosca - Accept external DLON, DLAT from ESMF interface
!  26 Feb 2013 - R. Yantosca - Now pass I_LO, J_LO to COMPUTE_GRID
!  28 Feb 2013 - R. Yantosca - Bug fix for GEOS-5 interface: Now call 
!                              Compute_Grid with 1..IIPAR, I..JJPAR
!  01 Jul 2013 - R. Yantosca - Don't use 1/2 sized polar boxes for GCAP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: JSP, JNP

      !=================================================================
      ! Define inputs depending on the grid that is selected
      !=================================================================
#if   defined( GRID4x5 )

      !-----------------------------
      ! Global 4 x 5 grid
      !-----------------------------
      JSP           = 1                  ! Lat index of S Pole
      JNP           = JM_WORLD           ! Lat index of N pole
      DLON          = 5.0d0              ! Delta-longitude array [degrees]
      DLAT          = 4.0d0              ! Delta-latitude array [degrees]
# if !defined( GCAP )
      DLAT(:,JSP,:) = 2.0d0              !  Half-sized S Pole boxes
      DLAT(:,JNP,:) = 2.0d0              !  Half-sized N Pole boxes
# endif

#elif defined( GRID2x25 )

      !-----------------------------
      ! Global 2 x 2.5 grid
      !-----------------------------
      JSP           = 1                  ! Lat index of S Pole
      JNP           = JM_WORLD           ! Lat index of N pole
      DLON          = 2.5d0              ! Delta-longitude array [degrees]
      DLAT          = 2.0d0              ! Delta-latitude array [degrees]
# if !defined( GCAP )
      DLAT(:,JSP,:) = 1.0d0              !  Half-sized S Pole boxes
      DLAT(:,JNP,:) = 1.0d0              !  Half-sized N Pole boxes
# endif

#elif defined( GRID1x125 )

      !-----------------------------
      ! Global 1 x 1.25 grid
      !-----------------------------
      JSP           = 1                  ! Lat index of S Pole
      JNP           = JM_WORLD           ! Lat index of N pole
      DLON          = 1.25d0             ! Delta-longitude array [degrees]
      DLAT          = 1.0d0              ! Delta-latitude array [degrees]
# if !defined( GCAP )
      DLAT(:,JSP,:) = 0.5d0              !  Half-sized S Pole boxes
      DLAT(:,JNP,:) = 0.5d0              !  Half-sized N Pole boxes
# endif

#elif defined( GRID1x1 )

      !-----------------------------
      ! Global 1 x 1 
      !-----------------------------
      JSP           = 1                  ! Lat index of S Pole
      JNP           = JM_WORLD           ! Lat index of N pole
      DLON          = 1.0d0              ! Delta-longitude array [degrees]
      DLAT          = 1.0d0              ! Delta-latitude array [degrees]
# if !defined( GCAP )
      DLAT(:,JSP,:) = 0.5d0              !  Half-sized S Pole boxes
      DLAT(:,JNP,:) = 0.5d0              !  Half-sized N Pole boxes
# endif

#elif defined( GRID05x0666 )

      CALL COMPUTE_GLOBAL_GRID( am_I_Root )

      !-------------------=---------
      ! Nested 0.5 x 0.666 grids
      !-----------------------------
      JSP           = 0                  ! Lat index of S Pole
      JNP           = 0                  ! Lat index of N pole
      DLON          = 2d0/3d0            ! Delta-longitude array [degrees]
      DLAT          = 0.5d0              ! Delta-latitude array [degrees]

#elif defined( GRID025x03125 )

      CALL COMPUTE_GLOBAL_GRID( am_I_Root )

      !-----------------------------
      ! Nested 0.25 x 0.3125 grids
      !-----------------------------
      JSP           = 0                  ! Lat index of S Pole
      JNP           = 0                  ! Lat index of N pole
      DLON          = 0.3125d0           ! Delta-longitude array [degrees]
      DLAT          = 0.25d0             ! Delta-latitude array [degrees]

#elif defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )

      !-----------------------------
      ! Connecting to GCM via ESMF
      !-----------------------------
      JSP           = 1
      JNP           = JM_WORLD

      ! NOTE: DLON, DLAT are defined in routine GIGC_Get_Options in
      ! ESMF/gigc_initialization_mod.F90, so no need to define them
      ! here.  (bmy, 12/3/12)

#endif

      !=================================================================
      ! Initialize the horizontal grid
      !=================================================================

      ! Allocate module arrays in grid_mod.F90
      CALL Init_Grid   ( am_I_Root = am_I_Root,     ! Are we on the root CPU?
     &                   IM        = IIPAR,         ! # of lons on this CPU
     &                   JM        = JJPAR,         ! # of lats on this CPU
     &                   LM        = LLPAR,         ! # of levs on this CPU
     &                   RC        = RC         )   ! Success or failure?

      ! Compute the horiziontal grid properties
      CALL Compute_Grid( am_I_Root = am_I_Root,     ! Are we on the root CPU?
     &                   I1        = 1,             ! Min lon index, this CPU
     &                   J1        = 1,             ! Min lat index, this CPU
     &                   L1        = 1,             ! Min lev index, this CPU
     &                   I2        = IIPAR,         ! Max lon index, this CPU
     &                   J2        = JJPAR,         ! Max lat index, this CPU
     &                   L2        = LLPAR,         ! Max lev index, this CPU
     &                   JSP       = JSP,           ! Lat index of South Pole
     &                   JNP       = JNP,           ! Lat index of North Pole
     &                   DLON      = DLON,          ! Delta-longitudes [deg]
     &                   DLAT      = DLAT,          ! Delta-latitudes  [deg]
     &                   I_LO      = I_LO,          ! Min global lon index
     &                   J_LO      = J_LO,          ! Min global lat index
     &                   RC        = RC         )   ! Success or failure?


      END SUBROUTINE INITIALIZE_GEOS_GRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_tracer_menu
!
! !DESCRIPTION: Subroutine READ\_TRACER\_MENU reads the TRACER MENU section 
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TRACER_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CHARPAK_MOD,  ONLY : ISDIGIT
      USE BIOFUEL_MOD,  ONLY : SET_BFTRACE
      USE BIOMASS_MOD,  ONLY : SET_BIOTRCE
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TRACER_MOD,   ONLY : ID_EMITTED,     ID_TRACER
      USE TRACER_MOD,   ONLY : SIM_TYPE,       N_TRACERS
      USE TRACER_MOD,   ONLY : TCVV,           TRACER_COEFF
      USE TRACER_MOD,   ONLY : TRACER_CONST,   TRACER_MW_G
      USE TRACER_MOD,   ONLY : TRACER_MW_KG,   TRACER_N_CONST
      USE TRACER_MOD,   ONLY : TRACER_NAME,    INIT_TRACER
      USE TRACER_MOD,   ONLY : XNUMOL,         XNUMOLAIR
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,   ONLY : ITS_A_HCN_SIM
      USE TRACER_MOD,   ONLY : ITS_A_MERCURY_SIM
      USE TRACERID_MOD, ONLY : TRACERID

      USE CMN_SIZE_MOD                    ! Size parameters
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now set LSPLIT correctly for Tagged Hg simulation (eck, bmy, 12/13/04)
!  (2 ) Now initialize ocean mercury module (sas, bmy, 1/20/05)
!  (3 ) Now set LSPLIT correctly for Tagged HCN/CH3CN sim (xyp, bmy, 6/30/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now reference XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Now move call to INIT_OCEAN_MERCURY to READ_MERCURY_MENU (bmy, 2/24/06)
!  (7 ) Now do not call SET_BIOTRCE anymore; it's obsolete (bmy, 4/5/06)
!  (8 ) Add SET_BIOTRCE to initialize IDBxxxs. (fp, 2/26/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2012 - R. Yantosca - Now pass CHEM_STATE as an argument (DEVEL only)
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  07 Nov 2012 - R. Yantosca - Now define Input_Opt%ITS_A_*_SIM fields
!  28 Oct 2013 - R. Yantosca - Set Input_Opt%ITS_A_SPECIALTY_SIM = .FALSE.
!                              when running GEOS-Chem in an ESMF environment
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N,  T,  C
      CHARACTER(LEN=255) :: C1, C2, LINE, NAME
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_TRACER_MENU begins here!
      !=================================================================

      ! NSRCX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:1' )
      READ( SUBSTRS(1:N), * ) SIM_TYPE

      ! NTRACE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:2' )
      READ( SUBSTRS(1:N), * ) N_TRACERS

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_tracer_menu:3' )
      ! NTRACE cannot exceed NNPAR

      IF ( N_TRACERS > NNPAR ) THEN
         CALL ERROR_STOP( 'Error: N_TRACERS > NNPAR!', 
     &                    'read_tracer_menu (input_mod.f)' )
      ENDIF
     
      ! Save N_TRACERS to Input_Opt
      Input_Opt%N_TRACERS = N_TRACERS

      ! Initialize tracer arrays in tracer_mod.f
      CALL INIT_TRACER( am_I_Root, Input_Opt, RC )

      !=================================================================
      ! Read tracer ID, emission ID, biomass ID, biofuel ID, & mol wt
      !=================================================================
      DO T = 1, N_TRACERS

         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_tracer_menu:4' )

         ! Save tracer number
         READ( SUBSTRS(1), * ) ID_TRACER(T)

         ! Save tracer name
         TRACER_NAME(T)  = TRIM( SUBSTRS(2) )

         ! Save tracer molecular wt [g/mole]
         READ( SUBSTRS(3), * ) TRACER_MW_G(T) 

         ! Save tracer molecular wt [kg/mole]
         TRACER_MW_KG(T) = TRACER_MW_G(T) * 1d-3

         ! Ratio of MW air / MW tracer 
         TCVV(T)         = 28.97d0 / TRACER_MW_G(T)

         ! Molecules tracer / kg tracer
         XNUMOL(T)       = 6.022d+23 / TRACER_MW_KG(T)
         
         !----------------------------------------------------
         ! If this is a family tracer, get member species
         !---------------------------------------------------- 
         IF ( N > 3 ) THEN
  
            ! Number of species in this family
            TRACER_N_CONST(T) = N-3
               
            ! Loop over substrings
            DO C = 1, N-3
         
               ! Trim the name
               NAME = TRIM( SUBSTRS(C+3) )
         
               !--------------------------------
               ! Emitted species denoted by ()
               !--------------------------------
               IF ( NAME(1:1) == '(' ) THEN 

                  ! Pull out the text between parentheses
                  NAME = NAME( 2:LEN_TRIM( NAME )-1 )

                  ! Set coeff of emitted tracer to 1 for now
                  TRACER_COEFF(T,C) = 1d0

                  ! Index of which family member is the emitted species
                  ID_EMITTED(T) = C

                  ! 2nd character of name is a number
                  IF ( ISDIGIT( NAME(1:1) ) ) THEN

                     ! Allow for tracers that contain >= 10C (tmf, 12/09/04) 
                     ! If 3rd character of name is a number
                     IF ( ISDIGIT( NAME(2:2) ) ) THEN

                        ! Coefficient of emitted tracer
                        READ( NAME(1:2), * ) TRACER_COEFF(T,C)

                        ! Get the rest of the name
                        NAME = NAME( 3:LEN_TRIM( NAME ) )

                     ELSE

                        ! Coefficient of emitted tracer
                        READ( NAME(1:1), * ) TRACER_COEFF(T,C)

                        ! Get the rest of the name
                        NAME = NAME( 2:LEN_TRIM( NAME ) )

                     ENDIF
                  ENDIF

                  ! Tracer constituent name
                  TRACER_CONST(T,C) = NAME

               !--------------------------------  
               ! Non-emitted species 
               !--------------------------------
               ELSE 

                  ! Set coeff of emitted tracer to 1 for now
                  TRACER_COEFF(T,C) = 1d0

                  ! 2nd character of name is a number
                  IF ( ISDIGIT( NAME(1:1) ) ) THEN

                     ! Allow for tracers that contain >= 10C (tmf, 12/09/04) 
                     ! If 3rd character of name is a number
                     IF ( ISDIGIT( NAME(2:2) ) ) THEN

                        ! Coefficient of emitted tracer
                        READ( NAME(1:2), * ) TRACER_COEFF(T,C)

                        ! Get the rest of the name
                        NAME = NAME( 3:LEN_TRIM( NAME ) )

                     ELSE

                        ! Coefficient of emitted tracer
                        READ( NAME(1:1), * ) TRACER_COEFF(T,C)

                        ! Get the rest of the name
                        NAME = NAME( 2:LEN_TRIM( NAME ) )

                     ENDIF

                  ENDIF

                  ! Tracer constituent name
                  TRACER_CONST(T,C) = NAME

               ENDIF

            ENDDO

         !----------------------------------------------------
         ! If not a family tracer, there is only one member
         !---------------------------------------------------- 
         ELSE

            ! Number of species in this family
            TRACER_N_CONST(T) = 1

            ! Index of which member is the emitted species
            ID_EMITTED(T)     = 0

            ! Coefficient of each tracer constituent
            TRACER_COEFF(T,1) = 1d0

            ! Tracer constituent name
            TRACER_CONST(T,1) = TRACER_NAME(T)

         ENDIF

      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_tracer_menu:5' )

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 
     &       'TRACER MENU (==> denotes SMVGEAR emitted species)' 
         WRITE( 6, '(  a)' ) REPEAT( '-', 48 )
         WRITE( 6, '(  a)' ) '  # Tracer          g/mole'

         ! Print info about each tracer
         DO T = 1, N_TRACERS

            ! Write tracer number, name, & mol wt
            WRITE( 6, 100 ) ID_TRACER(T), TRACER_NAME(T), TRACER_MW_G(T)
              
            ! If a family tracer (or just a tracer w/ emission)
            ! then also print info about species
            IF ( TRACER_N_CONST(T) > 1 .or. ID_EMITTED(T) > 0 ) THEN
            
               ! Loop over member species
               DO C = 1, TRACER_N_CONST(T)

                   ! Also flag which is the emitted tracer
                  IF ( ID_EMITTED(T) == C ) THEN
                     WRITE( 6,110 ) TRACER_COEFF(T,C), TRACER_CONST(T,C) 
                  ELSE
                     WRITE( 6,120 ) TRACER_COEFF(T,C), TRACER_CONST(T,C)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      ! Format statement
 100  FORMAT( I3, 1x, A10, 6x, F6.1 )
 110  FORMAT( 5x, '===> ', f4.1, 1x, A6  )
 120  FORMAT( 5x, '---> ', f4.1, 1x, A4  )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Split into tagged tracers (turn off for full-chemistry)
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! There are no tagged tracers for fullchem
         LSPLIT = .FALSE.

      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

         ! Need Hg0, Hg2, HgP for tagged Mercury
         LSPLIT = ( N_TRACERS > 3 )
         
      ELSE IF ( ITS_A_HCN_SIM() ) THEN

         ! Need HCN, CH3CN for HCN simulation
         LSPLIT = ( N_TRACERS > 2 )

      ELSE
         LSPLIT = ( N_TRACERS > 1 )

      ENDIF

      ! Set up tracer flags
      CALL TRACERID( am_I_Root )

      ! Set NBFTRACE in "biofuel_mod.f"
      CALL SET_BFTRACE( am_I_Root )

      ! Set NBIOTRACE in "biomass_mod.f"
      CALL SET_BIOTRCE( am_I_Root )

      ! Set counter
      CT1 = CT1 + 1

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%N_TRACERS              = N_TRACERS
      Input_Opt%LSPLIT                 = LSPLIT
      Input_Opt%ITS_A_RnPbBe_SIM       = ( SIM_TYPE == 1  )
      Input_Opt%ITS_A_CH3I_SIM         = ( SIM_TYPE == 2  )
      Input_Opt%ITS_A_FULLCHEM_SIM     = ( SIM_TYPE == 3  )
      Input_Opt%ITS_A_HCN_SIM          = ( SIM_TYPE == 4  )
      Input_Opt%ITS_A_TAGOX_SIM        = ( SIM_TYPE == 6  )
      Input_Opt%ITS_A_TAGCO_SIM        = ( SIM_TYPE == 7  )
      Input_Opt%ITS_A_C2H6_SIM         = ( SIM_TYPE == 8  )
      Input_Opt%ITS_A_CH4_SIM          = ( SIM_TYPE == 9  )
      Input_Opt%ITS_NOT_COPARAM_OR_CH4 = ( SIM_TYPE /= 9  )
      Input_Opt%ITS_AN_AEROSOL_SIM     = ( SIM_TYPE == 10 )
      Input_Opt%ITS_A_MERCURY_SIM      = ( SIM_TYPE == 11 )
      Input_Opt%ITS_A_CO2_SIM          = ( SIM_TYPE == 12 )
      Input_Opt%ITS_A_H2HD_SIM         = ( SIM_TYPE == 13 )
      Input_Opt%ITS_A_POPS_SIM         = ( SIM_TYPE == 14 )

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Set Input_Opt%ITS_A_SPECIALTY_SIM to FALSE, which will skip
      ! reading/allocating T_FULLGRID and AIRDEN_FULLGRID arrays.
      ! These are only needed for the offline "specialty" simulations
      ! so that the OH fields can be regridded properly. (bmy, 10/28/13)
      !-----------------------------------------------------------------
      Input_Opt%ITS_A_SPECIALTY_SIM    = .FALSE.
#else
      !-----------------------------------------------------------------
      !                %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Set Input_Opt%ITS_A_SPECIALTY_SIM to TRUE if we are performing
      ! any of the offline (non-full-chemistry) simulations that use
      ! offline OH fields.  This will cause both the T_FULLGRID and 
      ! AIRDEN_FULLGRID arrays to be allocated. (bmy, 10/28/13)
      !-----------------------------------------------------------------
      Input_Opt%ITS_A_SPECIALTY_SIM    =
     &     ( Input_Opt%ITS_A_CH3I_SIM      .or.
     &       Input_Opt%ITS_A_HCN_SIM       .or.
     &       Input_Opt%ITS_A_TAGOX_SIM     .or.
     &       Input_Opt%ITS_A_TAGCO_SIM     .or.
     &       Input_Opt%ITS_A_C2H6_SIM      .or.
     &       Input_Opt%ITS_A_CH4_SIM       .or.
     &       Input_Opt%ITS_AN_AEROSOL_SIM  .or.
     &       Input_Opt%ITS_A_MERCURY_SIM   .or.
     &       Input_Opt%ITS_A_CO2_SIM       .or.
     &       Input_Opt%ITS_A_H2HD_SIM      .or.
     &       Input_Opt%ITS_A_POPS_SIM           )
#endif

      ! Print info about each tracer
      DO T = 1, N_TRACERS
         Input_Opt%ID_EMITTED(T)       = ID_EMITTED(T)
         Input_Opt%ID_TRACER(T)        = ID_TRACER(T)
         Input_Opt%TRACER_NAME(T)      = TRIM( TRACER_NAME(T) )
         Input_Opt%TRACER_N_CONST(T)   = TRACER_N_CONST(T)
         Input_Opt%TRACER_MW_G(T)      = TRACER_MW_G(T)
         Input_Opt%TRACER_MW_KG(T)     = TRACER_MW_KG(T)
         Input_OPt%TCVV(T)             = TCVV(T)  
         Input_Opt%XNUMOL(T)           = XNUMOL(T)

         ! Loop over member species
         DO C = 1, TRACER_N_CONST(T)
            Input_Opt%TRACER_COEFF(T,C) = TRACER_COEFF(T,C)
            Input_Opt%TRACER_CONST(T,C) = TRIM( TRACER_CONST(T,C) )
         ENDDO
      ENDDO

      ! Return success
      RC = GIGC_SUCCESS

      END SUBROUTINE READ_TRACER_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_aerosol_menu
!
! !DESCRIPTION: Subroutine READ\_AEROSOL\_MENU reads the AEROSOL MENU 
!  section of the GEOS-Chem input file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_AEROSOL_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSULF, LCARB, LSOA,   LSVPOA
      USE LOGICAL_MOD,  ONLY : LDUST, LDEAD, LSSALT, LCRYST
      USE LOGICAL_MOD,  ONLY : LDICARB
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE TRACER_MOD,   ONLY : SALA_REDGE_um,      SALC_REDGE_um
      USE TRACER_MOD,   ONLY : ITS_AN_AEROSOL_SIM, ITS_A_FULLCHEM_SIM
      USE TRACERID_MOD, ONLY : IDTDMS,   IDTSO2,   IDTSO4,  IDTSO4s 
      USE TRACERID_MOD, ONLY : IDTMSA,   IDTNH3,   IDTNH4,  IDTNITs 
      USE TRACERID_MOD, ONLY : IDTAS,    IDTAHS,   IDTLET,  IDTNH4aq 
      USE TRACERID_MOD, ONLY : IDTSO4aq, IDTBCPO,  IDTBCPI, IDTOCPO 
      USE TRACERID_MOD, ONLY : IDTOCPI
      USE TRACERID_MOD, ONLY : IDTDST1,  IDTDST2,  IDTDST3, IDTDST4
      USE TRACERID_MOD, ONLY : IDTSALA,  IDTSALC 
      USE TRACERID_MOD, ONLY : IDTSOAG,  IDTSOAM
      ! SOAupdate: remove old tracers and add new mtp (hotp 5/24/10)
      USE TRACERID_MOD, ONLY : IDTMTPA,  IDTLIMO,  IDTMTPO
      USE TRACERID_MOD, ONLY : IDTBENZ,  IDTTOLU,  IDTXYLE
      USE TRACERID_MOD, ONLY : IDTTSOA1, IDTTSOA2, IDTTSOA3
      USE TRACERID_MOD, ONLY : IDTISOA1, IDTISOA2, IDTISOA3
      USE TRACERID_MOD, ONLY : IDTASOA1, IDTASOA2, IDTASOA3, IDTASOAN
      ! semivolpoa : add POA (hotp 8/23/09)
      ! POA now 1 and 2 (hotp 10/11/09), and NAP (hotp 7/28/10)
      USE TRACERID_MOD, ONLY : IDTPOA1,  IDTPOA2,  IDTPOG1,  IDTPOG2
      USE TRACERID_MOD, ONLY : IDTOPOA1, IDTOPOA2, IDTOPOG1, IDTOPOG2
      USE TRACERID_MOD, ONLY : IDTNAP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now reference LSOA (bmy, 9/28/04)
!  (2 ) Now stop run if LSOA=T and SOA tracers are undefined (bmy, 11/19/04)
!  (3 ) Now reference LCRYST from "logical_mod.f".  Also now check to make
!        prevent aerosol tracers from being undefined if the corresponding
!        logical switch is set. (cas, bmy, 1/14/05)
!  (4 ) Now also require LSSALT=T when LSULF=T, since we now compute the 
!        production of SO4 and NIT w/in the seasalt aerosol (bec, bmy, 4/13/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now update error check for SOG4, SOA4 (dkh, bmy, 6/1/06)
!  (7 ) Add LDICARB switch to cancel SOG condensation onto OC aerosols.
!      (ccc, tmf, 3/10/09)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N, T, I
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_AEROSOL_MENU begins here!
      !=================================================================

      ! Location string for ERROR_STOP
      LOCATION = 'READ_AEROSOL_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Use online sulfate aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:1' )
      READ( SUBSTRS(1:N), * ) LSULF

      ! Use crystalline sulfate aerosols  
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:2' )
      READ( SUBSTRS(1:N), * ) LCRYST

      ! Use online carbon aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:3' )
      READ( SUBSTRS(1:N), * ) LCARB

      ! Use secondary organic aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:4' )
      READ( SUBSTRS(1:N), * ) LSOA

      ! SOAupdate: Add Semi-volatile POA switch (hotp 8/9/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:5' )
      READ( SUBSTRS(1:N), * ) LSVPOA

      ! SOAupdate: Add NAP/IVOC surrogate emissions (hotp 8/9/09)
      ! SOAupdate: set NAP emissions (in Tg/yr)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:6' )
      READ( SUBSTRS(1:N), * ) Input_Opt%NAPEMISS

      ! SOAupdate: Add scaling for POA emiss (hotp 9/21/09)
      ! Fractional scaling, recommended value: 1.27
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:7' )
      READ( SUBSTRS(1:N), * ) Input_Opt%POAEMISSSCALE

      ! Use online dust aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:8' )
      READ( SUBSTRS(1:N), * ) LDUST

      ! Use DEAD dust mobilization (=T) or GINOUX (=F)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:9' )
      READ( SUBSTRS(1:N), * ) LDEAD      

      ! Use online sea-salt aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:10' )
      READ( SUBSTRS(1:N), * ) LSSALT      

      ! Accum mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aerosol_menu:11' )
      DO T = 1, N
         READ( SUBSTRS(T), * ) SALA_REDGE_um(T)
      ENDDO

      ! Coarse mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_aerosol_menu:12' )
      DO T = 1, N
         READ( SUBSTRS(T), * ) SALC_REDGE_um(T)
      ENDDO

      ! Switch to comment the SOG condensation in carbon_mod.f (ccc, 3/10/09)
      ! Use online dicarbonyl chemistry 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_aerosol_menu:13' )
      READ( SUBSTRS(1:N), * ) LDICARB

      !=================================================================
      ! Error checks
      !=================================================================

      ! Make sure that SALA, SALC bins are contiguous
      IF ( SALA_REDGE_um(2) /= SALC_REDGE_um(1) ) THEN
         MSG = 'SALA and SALC bin edges are not contiguous!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Turn off switches for simulations that don't use aerosols
      IF ( ( .not. ITS_A_FULLCHEM_SIM() )  .and. 
     &     ( .not. ITS_AN_AEROSOL_SIM() ) ) THEN
         LSULF  = .FALSE.
         LCRYST = .FALSE.
         LCARB  = .FALSE.
         LSOA   = .FALSE.
         LDUST  = .FALSE.
         LDEAD  = .FALSE.
         LSSALT = .FALSE.
         LSVPOA = .FALSE.
      ENDIF

      !%%% The cryst/aq code is currently under development so make
      !%%% sure that LCRYST = FALSE for now until further notice
      !%%% (rjp, bmy, 3/15/05)
      LCRYST = .FALSE.

      !---------------------------------
      ! Error check SEASALT AEROSOLS
      !---------------------------------
      I = IDTSALA + IDTSALC

      IF ( LSSALT ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check SULFUR AEROSOLS
      !---------------------------------
      I = IDTDMS + IDTSO2 + IDTSO4 + IDTSO4s + 
     &    IDTMSA + IDTNH3 + IDTNH4 + IDTNITs

      IF ( LSULF ) THEN

         ! We now compute the production of SO4s and NITs, so when 
         ! LSULF=T, then we must also have LSSALT=T (bec, bmy, 4/13/05)
         IF ( .not. LSSALT ) THEN 
            MSG = 'LSULF=T now also requires LSSALT=T!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF

         ! Stop w/ error if everything is undefined
         IF ( I == 0 ) THEN
            MSG = 'LSULF=T but ONLINE SULFUR AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF

      ELSE

         ! If LSULF=F but we have defined tracers, stop w/ error
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SULFUR AEROSOLS if LSULF=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF

      ENDIF

      !---------------------------------
      ! Error check CRYST /AQ AEROSOLS
      !---------------------------------
      I = IDTAS + IDTAHS + IDTLET + IDTNH4aq + IDTSO4aq

      IF ( LCRYST ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LCRYST=T but ONLINE CRYST/AQ AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE CRYST/AQ AEROSOLS if LCRYST=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check CARBON AEROSOLS
      !---------------------------------
      ! SOAupdate: Add POA (hotp 10/11/09)
      I = IDTBCPO + IDTBCPI + IDTOCPO + IDTOCPI + IDTPOA1

      IF ( LCARB ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LCARB=T but ONLINE CARBON AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE CARBON AEROSOLS if LCARB=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      IF ( LSVPOA .and. ( .NOT. LSOA ) ) THEN
         MSG = 'Semivolatile POA requires LSOA=T'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! SOAupdate: Error check (hotp 8/24/09)
      ! OCPI and OCPO are the non-volatile POA tracers
      ! POA (along w/ POG, OPOA, and OPOG) are the semivol POA tracers
      ! You can't have both!
      I = IDTOCPI + IDTOCPO
      IF ( IDTPOA1 > 0 ) THEN
         IF ( I > 0 ) THEN
            MSG = 'Semivolatile POA tracer is defined in addition to
     &             Nonvolatile POA'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
         IF ( ( .NOT. LSOA ) .or. (.NOT. LSVPOA ) ) THEN
            MSG = 'Semivolatile POA requires LSOA=T and LSVPOA=T'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      ! SOAupdate
      ! Options for organic aerosol tracers:
      ! IF LSOA = F: only OCPI and OCPO
      ! IF LSOA = T:
      !   OCPI OCPO SOA (non-vol + original traditional)
      !   POA POG OPOA OPOG SOA BTX NAP (semivol + orig trad + IVOC )
      ! NAP emissions are set in input.geos
      ! LSVPOA is just a check (doesn't do anything hotp 7/21/10)
      I = IDTPOA1  + IDTPOA2  + IDTPOG1  + IDTPOG2  + 
     &    IDTOPOA1 + IDTOPOA2 + IDTOPOG1 + IDTOPOG2
      IF ( LSVPOA ) THEN
         IF ( I < 8 ) THEN
            MSG = 'Not enough semivolatile POA tracers!'
                     CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
         IF ( IDTNAP == 0 ) THEN
            MSG = 'Semivolatile POA requires IVOCs/NAP!'
                     CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check 2dy ORG AEROSOLS
      !---------------------------------
      ! SOAupdate: update for new mtp (hotp 5/24/10)
      I = IDTMTPA  + IDTLIMO  + IDTMTPO  + 
     &    IDTBENZ  + IDTTOLU  + IDTXYLE  +
     &    IDTTSOA1 + IDTTSOA2 + IDTTSOA3 +
     &    IDTISOA1 + IDTISOA2 + IDTISOA3 +
     &    IDTASOA1 + IDTASOA2 + IDTASOA3 + IDTASOAN +
     &    IDTSOAG  + IDTSOAM

      IF ( LSOA ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSOA=T but ONLINE 2dy ORG AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE 2dy ORG AEROSOLS if LSOA=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check DUST AEROSOLS
      !---------------------------------
      I = IDTDST1 + IDTDST2 + IDTDST3 + IDTDST4

      IF ( LDUST ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LDUST=T but ONLINE DUST AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE DUST AEROSOLS if LDUST=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !---------------------------------
      ! Error check SEASALT AEROSOLS
      !---------------------------------
      I = IDTSALA + IDTSALC

      IF ( LSSALT ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LSULF              = LSULF  
      Input_Opt%LCRYST             = LCRYST
      Input_Opt%LCARB              = LCARB
      Input_Opt%LSOA               = LSOA 
      Input_Opt%LSVPOA             = LSVPOA
      Input_Opt%LDUST              = LDUST 
      Input_Opt%LDEAD              = LDEAD
      Input_Opt%LSSALT             = LSSALT
      Input_Opt%LDICARB            = LDICARB  
      Input_Opt%SALA_REDGE_um(1:2) = SALA_REDGE_um(1:2)
      Input_Opt%SALC_REDGE_um(1:2) = SALC_REDGE_um(1:2)

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'AEROSOL MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 100     ) 'Online SULFATE AEROSOLS?    : ', LSULF
         WRITE( 6, 100     ) 'Online CRYST & AQ AEROSOLS? : ', LCRYST
         WRITE( 6, 100     ) 'Online CARBON AEROSOLS?     : ', LCARB
         WRITE( 6, 100     ) 'Online 2dy ORGANIC AEROSOLS?: ', LSOA
         WRITE( 6, 100     ) 'Semivolatile POA?           : ', LSVPOA
         WRITE( 6, 105     ) 'Annual NAP emiss [TgC/yr]   : ',
     &                        Input_Opt%NAPEMISS
         WRITE( 6, 105     ) 'POA emission scale [dimless]: ', 
     &                        Input_Opt%POAEMISSSCALE
         WRITE( 6, 100     ) 'Online DUST AEROSOLS?       : ', LDUST
         WRITE( 6, 100     ) 'Use DEAD dust emissions?    : ', LDEAD
         WRITE( 6, 100     ) 'Online SEA SALT AEROSOLS?   : ', LSSALT
         WRITE( 6, 110     ) 'Accum  SEA SALT radii [um]  : ', 
     &                        SALA_REDGE_um(1), SALA_REDGE_um(2)
         WRITE( 6, 110     ) 'Coarse SEA SALT radii [um]  : ',
     &                        SALC_REDGE_um(1), SALC_REDGE_um(2)
      ENDIF

 100  FORMAT( A, L5     )
 110  FORMAT( A, f6.2, ' - ', f6.2 )
 105  FORMAT( A, f6.2 )

      END SUBROUTINE READ_AEROSOL_MENU
!EOC
#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tomas_microphys
!
! !DESCRIPTION: Subroutine INIT\_TOMAS\_MICROPHYS will initialize the 
!  TOMAS microphysics package.  This replaces the former subroutine
!  READ\_MICROPHYSICS\_MENU.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TOMAS_MICROPHYSICS( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TOMAS_MOD,          ONLY : INIT_TOMAS
      USE TRACERID_MOD,       ONLY : IDTNK1
      USE TRACERID_MOD,       ONLY : IDTSF1
      USE TRACERID_MOD,       ONLY : IDTSS1
      USE TRACERID_MOD,       ONLY : IDTECOB1, IDTECIL1, IDTOCOB1
      USE TRACERID_MOD,       ONLY : IDTOCIL1, IDTH2SO4, IDTDUST1
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REMARKS:
!  We now invoke TOMAS by compiling GEOS-Chem and setting either the TOMAS=yes
!  (30 bins, default) or TOMAS40=yes (40 bins, optional) switches.  The old
!  LTOMAS logical switch is now obsolete because all of the TOMAS code is
!  segregated from the rest of GEOS-Chem with #ifdef blocks.  Therefore,
!  we no longer need to read the microphysics menu, but we still need to
!  apply some error checks and then call INIT_TOMAS. (bmy, 4/23/13)
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now read LNEI99 -- switch for EPA/NEI99 emissions (bmy, 11/5/04)
!  (2 ) Now read LAVHRRLAI-switch for using AVHRR-derived LAI (bmy, 12/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read LMEGAN -- switch for MEGAN biogenics (tmf, bmy, 10/20/05)
!  (5 ) Now read LEMEP -- switch for EMEP emissions (bdf, bmy, 11/1/05)
!  (6 ) Now read LGFED2BB -- switch for GFED2 biomass emissions (bmy, 4/5/06)
!  (7 ) Now read LOTDLIS, LCTH, LMFLUX, LPRECON for lightning options 
!        (bmy, 5/10/06)
!  (8 ) Now read LBRAVO for BRAVO Mexican emissions (rjp, kfb, bmy, 6/26/06)
!  (9 ) Now read LEDGAR for EDGAR emissions (avd, bmy, 7/11/06)
!  (10) Now read LSTREETS for David Streets' emissions (bmy, 8/17/06)
!  (11) Kludge: Reset LMFLUX or LPRECON to LCTH, as the MFLUX and PRECON
!        lightning schemes have not yet been implemented.  Rename LOTDLIS
!        to LOTDREG.  Also read LOTDLOC for the OTD-LIS local redistribution
!        of lightning flashes (cf B. Sauvage).  Make sure LOTDREG and 
!        LOTDLOC are not both turned on at the same time. (bmy, 1/31/07)
!  (12) Add LOTDSCALE to the list of LNOx options (ltm, bmy, 9/24/07)
!  (13) Add new error traps for OTD-LIS options, dependent on met field type
!        (ltm, bmy, 11/29/07)
!  (14) Bug fix, create string variables for ERROR_STOP (bmy, 1/24/08)
!  (15) Now read LCAC for CAC emissions (amv, 1/09/2008)
!  (16) Now read LEDGARSHIP, LARCSHIP and LEMEPSHIP (phs, 12/5/08)
!  (17) Fixed typo in message for GEOS-3 (bmy, 10/30/08)
!  (18) Now read LVISTAS (amv, 12/2/08)
!  (19) Now read L8DAYBB, L3HRBB and LSYNOPBB for GFED2 8-days and 3hr
!        emissions, and LICARTT for corrected EPA (phs, yc, 12/17/08)
!  (20) Add a specific switch for MEGAN emissions for monoterpenes and MBO
!       (ccc, 2/2/09)
!  (21) Now read LICOADSSHIP (cklee, 6/30/09)
!  (22) Bug fix: for now, if LEMEPSHIP is turned on but LEMEP is turned off,
!        just turn off LEMEPSHIP and print a warning msg. (mak, bmy, 10/18/09)
!  (23) Now accounts for NEI2005 (amv, phs, 10/9/09)
!  (24) Included optional flag for using MODIS LAI data (mpb,2009).
!  (25) Included optional flag for using PCEEA model (mpb, 2009)
!  (26) Now force settings for EU, NA, CC nested grids (amv, bmy, 12/18/09)
!  (27) Now force MEGAN to use MODIS LAI (ccarouge, bmy, 2/24/10)
!  (28) Add separate switch for NOx fertilizer. (fp, 2/29/10)
!  (29) Add scaling for isoprene and NOx emissions. (fp, 2/29/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  23 Apr 2013 - R. Yantosca - Renamed to INIT_TOMAS_MICROPHYS
!  30 Jan 2014 - R. Yantosca - INIT_TOMAS accepts am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N, I
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! INIT_TOMAS_MICROPHYSICS begins here!
      !=================================================================

      ! Location string for ERROR_STOP
      LOCATION = 'INIT_TOMAS_MICROPHYS ("input_mod.f")'

      ! Halt with error if we are trying to run TOMAS in a simulation
      ! that does not have any defined aerosols
      ! Turn off switches for simulations that don't use microphysics
      IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and.
     &     ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
         MSG = 'TOMAS needs to run with either a full-chemistry ' //
     &         'or offline aerosol simulation!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Halt with error if none of the TOMAS aerosol species are defined
      I =  IDTNK1 + IDTSF1 + IDTSS1 + IDTECOB1 + IDTECIL1 + 
     &     IDTOCOB1 + IDTOCIL1  + IDTDUST1
      IF ( I == 0 ) THEN 
         MSG = 'None of the TOMAS aerosols are defined!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF
      
      ! Halt with error if sulfate aerosols are not defined
      IF( IDTSF1 > 0 .and. ( .not. Input_Opt%LSULF ) ) THEN
         MSG = 'Need LSULF on for TOMAS-Sulfate to work (for now)'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Halt with error if carbonaceous aerosols are not defined
      I = IDTECOB1 + IDTECIL1 + IDTOCOB1 + IDTOCIL1
      IF ( I > 0 .and. (.not. Input_Opt%LCARB ) ) THEN
         MSG = 'Need LCARB on for TOMAS-carb to work (for now)'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Halt with error if dust aerosols are turned on.
      ! TOMAS defines its own dust aerosol tracers.
!      IF ( IDTDUST1 > 0 .AND. Input_Opt%LDUST ) THEN
!         MSG = 'Need to turn off LDUST for TOMAS Dust to work'
!         CALL ERROR_STOP( MSG, LOCATION )
!      ENDIF

      !=================================================================
      ! All error checks are satisfied, so initialize TOMAS
      !=================================================================
      CALL INIT_TOMAS( am_I_Root, Input_Opt, RC )
 
      END SUBROUTINE INIT_TOMAS_MICROPHYSICS
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_emissions_menu
!
! !DESCRIPTION: Subroutine READ\_EMISSIONS\_MENU reads the EMISSIONS MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EMISSIONS_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE BROMOCARB_MOD,      ONLY : Br_SCALING
      USE CMN_O3_MOD               ! FSCALYR
      USE CMN_SIZE_MOD             ! Size parameters
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE EMISSIONS_MOD,      ONLY : ISOP_SCALING 
      USE EMISSIONS_MOD,      ONLY : NOx_SCALING
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE LOGICAL_MOD
      USE MODIS_LAI_MOD,      ONLY : INIT_MODIS_LAI
      USE RCP_MOD,            ONLY : RCPNAME, RCPYEAR
      USE TIME_MOD,           ONLY : SET_HISTYR
      USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now read LNEI99 -- switch for EPA/NEI99 emissions (bmy, 11/5/04)
!  (2 ) Now read LAVHRRLAI-switch for using AVHRR-derived LAI (bmy, 12/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read LMEGAN -- switch for MEGAN biogenics (tmf, bmy, 10/20/05)
!  (5 ) Now read LEMEP -- switch for EMEP emissions (bdf, bmy, 11/1/05)
!  (6 ) Now read LGFED2BB -- switch for GFED2 biomass emissions (bmy, 4/5/06)
!  (7 ) Now read LOTDLIS, LCTH, LMFLUX, LPRECON for lightning options 
!        (bmy, 5/10/06)
!  (8 ) Now read LBRAVO for BRAVO Mexican emissions (rjp, kfb, bmy, 6/26/06)
!  (9 ) Now read LEDGAR for EDGAR emissions (avd, bmy, 7/11/06)
!  (10) Now read LSTREETS for David Streets' emissions (bmy, 8/17/06)
!  (11) Kludge: Reset LMFLUX or LPRECON to LCTH, as the MFLUX and PRECON
!        lightning schemes have not yet been implemented.  Rename LOTDLIS
!        to LOTDREG.  Also read LOTDLOC for the OTD-LIS local redistribution
!        of lightning flashes (cf B. Sauvage).  Make sure LOTDREG and 
!        LOTDLOC are not both turned on at the same time. (bmy, 1/31/07)
!  (12) Add LOTDSCALE to the list of LNOx options (ltm, bmy, 9/24/07)
!  (13) Add new error traps for OTD-LIS options, dependent on met field type
!        (ltm, bmy, 11/29/07)
!  (14) Bug fix, create string variables for ERROR_STOP (bmy, 1/24/08)
!  (15) Now read LCAC for CAC emissions (amv, 1/09/2008)
!  (16) Now read LEDGARSHIP, LARCSHIP and LEMEPSHIP (phs, 12/5/08)
!  (17) Fixed typo in message for GEOS-3 (bmy, 10/30/08)
!  (18) Now read LVISTAS (amv, 12/2/08)
!  (19) Now read L8DAYBB, L3HRBB and LSYNOPBB for GFED2 8-days and 3hr
!        emissions, and LICARTT for corrected EPA (phs, yc, 12/17/08)
!  (20) Add a specific switch for MEGAN emissions for monoterpenes and MBO
!       (ccc, 2/2/09)
!  (21) Now read LICOADSSHIP (cklee, 6/30/09)
!  (22) Bug fix: for now, if LEMEPSHIP is turned on but LEMEP is turned off,
!        just turn off LEMEPSHIP and print a warning msg. (mak, bmy, 10/18/09)
!  (23) Now accounts for NEI2005 (amv, phs, 10/9/09)
!  (24) Included optional flag for using MODIS LAI data (mpb,2009).
!  (25) Included optional flag for using PCEEA model (mpb, 2009)
!  (26) Now force settings for EU, NA, CC nested grids (amv, bmy, 12/18/09)
!  (27) Now force MEGAN to use MODIS LAI (ccarouge, bmy, 2/24/10)
!  (28) Add separate switch for NOx fertilizer. (fp, 2/29/10)
!  (29) Add scaling for isoprene and NOx emissions. (fp, 2/29/10)
!  (30) Remove depreciated lightning options. (ltm, 1/25,11)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  27 Aug 2010 - R. Yantosca - Added warning msg for MERRA
!  29 Jul 2011 - L. Zhang    - Fix bug that turns off CAC/BRAVO emissions
!                             inadvertently during nested NA simulations
!  11 Aug 2011 - E. Leibensperger - Added flag for historical emissions and
!                                   base year
!  07 Sep 2011 - P. Kasibhatla - Add modifications for GFED3
!  14 Feb 2012 - R. Yantosca - Reorganize error checks for logical switches
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  05 Apr 2012 - R. Yantosca - Now call INIT_MODIS_LAI
!  05 Apr 2012 - R. Yantosca - Reorganized USE statements for clarity
!  10 Apr 2012 - R. Yantosca - Bug fix: do not turn off LAVHRRLAI or LMODISLAI
!                              when emissions are turned off.  LAI is used
!                              in other areas of the code.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP, LRCPSHIP, LRCPAIR to Input_Opt
!  31 Jul 2013 - M. Sulprizio- Now copy LAEIC to Input_Opt; Add check to
!                              make sure LAEIC and LRCPAIR are not both on
!  22 Aug 2013 - R. Yantosca - Now read path for soil NOx restart file
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocessor switch
!  03 Oct 2013 - M. Sulprizio- Removed obsolete options for LAVHRRLAI and
!                              LMODISLAI. MODIS LAI data are now read from
!                              netCDF files.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOC, SOIL_RST_FILE
![eml
      INTEGER            :: HISTYR
!eml]
      !=================================================================
      ! READ_EMISSIONS_MENU begins here!
      !=================================================================

      ! Location for error messages
      LOC = 'READ_EMISSIONS_MENU ("input_mod.F")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOC )
      ENDIF

      ! Turn on emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:1' )
      READ( SUBSTRS(1:N), * ) LEMIS

      ! Emissions timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:2' )
      READ( SUBSTRS(1:N), * ) TS_EMIS

      ! Include anthropogenic emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:3' )
      READ( SUBSTRS(1:N), * ) LANTHRO

      ! Scale 1985 to year ____ ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:4' )
      READ( SUBSTRS(1:N), * ) FSCALYR

      ! Include EMEP (Europe) anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:5' )
      READ( SUBSTRS(1:N), * ) LEMEP

      ! Include BRAVO (Mexico) anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:6' )
      READ( SUBSTRS(1:N), * ) LBRAVO

      ! Include EDGAR anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:7' )
      READ( SUBSTRS(1:N), * ) LEDGAR

      ! Include David Streets anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:8' )
      READ( SUBSTRS(1:N), * ) LSTREETS

      ! Include CAC anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:9' )
      READ( SUBSTRS(1:N), * ) LCAC

      ! Include NEI2005 emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:10' )
      READ( SUBSTRS(1:N), * ) LNEI05

      ! Use RETRO emissions? (wfr, 3/8/11)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:11' )
      READ( SUBSTRS(1:N), * ) LRETRO

      ! AEIC aircraft emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12' )
      READ( SUBSTRS(1:N), * ) LAEIC

      ! Use RCP emissions? (cdh, 10/18/11)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12b' )
      READ( SUBSTRS(1:N), * ) LRCP

      ! Which RCP scenario? (cdh, 10/18/11)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12c' )
      READ( SUBSTRS(1:N), * ) RCPNAME

      ! Which RCP year? (cdh, 10/18/11)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:12d' )
      READ( SUBSTRS(1:N), * ) RCPYEAR

      ! Use EPA/NEI99 emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:13' )
      READ( SUBSTRS(1:N), * ) LNEI99
      
      ! Include ICARTT-based corrections ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:14' )
      READ( SUBSTRS(1:N), * ) LICARTT

      ! Include VISTAS anthro emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:15' )
      READ( SUBSTRS(1:N), * ) LVISTAS

      ! Include biofuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:16' )
      READ( SUBSTRS(1:N), * ) LBIOFUEL

      ! Include biogenic emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:17' )
      READ( SUBSTRS(1:N), * ) LBIOGENIC

      ! Use MEGAN biogenic emissions for ISOP?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:18' )
      READ( SUBSTRS(1:N), * ) LMEGAN

      ! Use PCEEA model for biogenic emissions for? (mpb,2009)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:19' )
      READ( SUBSTRS(1:N), * ) LPECCA

      ! Use MEGAN biogenic emissions for MONOT and MBO ? (ccc, 2/2/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:20' )
      READ( SUBSTRS(1:N), * ) LMEGANMONO

      ! Isoprene Emissions scaling (fp, 06/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:21' )
      READ( SUBSTRS(1:N), * ) ISOP_SCALING

      ! Include biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:22' )
      READ( SUBSTRS(1:N), * ) LBIOMASS

      ! Seasonal biomass?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:23' )
      READ( SUBSTRS(1:N), * ) LBBSEA

      ! Scaled to TOMSAI?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:24' )
      READ( SUBSTRS(1:N), * ) LTOMSAI

      ! Separator line (start of GFED2 biomass emissions)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:25' )

      ! Use monthly GFED2 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:26' )
      READ( SUBSTRS(1:N), * ) LGFED2BB

      ! Use 8-day GFED2 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:27' )
      READ( SUBSTRS(1:N), * ) L8DAYBB

      ! Use 3-hr GFED2 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:28' )
      READ( SUBSTRS(1:N), * ) L3HRBB

      ! Use 3-hr synoptic GFED2 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:29' )
      READ( SUBSTRS(1:N), * ) LSYNOPBB

      ! Separator line (start of GFED3 biomass emissions)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:30' )

      ! Use monthly GFED3 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:31' )
      READ( SUBSTRS(1:N), * ) LGFED3BB

      ! Use daily GFED3 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:32' )
      READ( SUBSTRS(1:N), * ) LDAYBB3

      ! Use 3-hr GFED3 biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:33' )
      READ( SUBSTRS(1:N), * ) L3HRBB3

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:34' )

      ! Use aircraft NOx from RCP
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:34b' )
      READ( SUBSTRS(1:N), * ) LRCPAIR

      ! Use lightning NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:35' )
      READ( SUBSTRS(1:N), * ) LLIGHTNOX

      ! Apply OTD-LIS local redistribution for lightning spatial and
      ! seasonal constraint
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:36' )
      READ( SUBSTRS(1:N), * ) LOTDLOC

      ! Use soil NOx
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:37' )
      READ( SUBSTRS(1:N), * ) LSOILNOX

      ! Soil NOx restart file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:37a' )
      READ( SUBSTRS(1:N), '(a)' ) SOIL_RST_FILE

      ! separate use fertilizer and soil NOx (fp, 06/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:38' )
      READ( SUBSTRS(1:N), * ) LFERTILIZERNOX
      
      !(FP, 15/12/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:39' )
      READ( SUBSTRS(1:N), * ) NOx_SCALING

      ! Separator line (start of ship emissions)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:40' )

      ! Use ship EDGAR ship emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:41' )
      READ( SUBSTRS(1:N), * ) LEDGARSHIP

      ! Use ICOADS (NOx, SO2, CO) ship  emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:42' )
      READ( SUBSTRS(1:N), * ) LICOADSSHIP

      ! Use RCP (NOx, SO2, CO) ship  emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:42b' )
      READ( SUBSTRS(1:N), * ) LRCPSHIP

      ! Use ship EMEP emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:43' )
      READ( SUBSTRS(1:N), * ) LEMEPSHIP

      ! Use ship SO2 Corbett emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:44' )
      READ( SUBSTRS(1:N), * ) LSHIPSO2

      ! Use ship ARCTAS (SO2, CO2) emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:45' )
      READ( SUBSTRS(1:N), * ) LARCSHIP

      ! Use COOKE over North AMerica for BC/OC?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:46' )
      READ( SUBSTRS(1:N), * ) LCOOKE

      ! Use historical emissions (EDGAR NOx/SO2 and Bond BC/OC) for 
      ! 1890-2000?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:49' )
      READ( SUBSTRS(1:N), * ) LHIST
      
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:50' )
      READ( SUBSTRS(1:N), * ) HISTYR
      CALL SET_HISTYR( HISTYR )

      ! Separator line (start of bromine emissions)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:51' )

      ! Use the Warwick(2006), scenario 3? (400Gg/yr of CHBr3)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:52' )
      READ( SUBSTRS(1:N), * ) LWARWICK_VSLS

      ! Use sea-salt Br2 emissions? Currently Monahan w/ Yang et al. 2005
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:53' )
      READ( SUBSTRS(1:N), * ) LSSABr2

      ! Set PBL BrO concentrations equal to 1ppt during the day?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:54' )
      READ( SUBSTRS(1:N), * ) LFIX_PBL_BRO

      ! Bromine Emissions scaling
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:55' )
      READ( SUBSTRS(1:N), * ) Br_SCALING

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_emissions_menu:56' )

      !=================================================================
      ! Error check logical flags
      !=================================================================

      !-----------------------------------------------------------------
      ! Turn off individual emissions flags if the master emissions flag
      ! LEMIS is turned off.  This will prevent inadvertent calls to the
      ! emissions routines. (bmy, 2/2/12)
      !-----------------------------------------------------------------

      ! Master switches
      LANTHRO           = ( LANTHRO        .and. LEMIS     )
      LBIOFUEL          = ( LBIOFUEL       .and. LEMIS     )
      LBIOGENIC         = ( LBIOGENIC      .and. LEMIS     )
      LBIOMASS          = ( LBIOMASS       .and. LEMIS     )

      ! Anthro
      LAEIC             = ( LAEIC          .and. LANTHRO   )
      LEMEP             = ( LEMEP          .and. LANTHRO   )
      LBRAVO            = ( LBRAVO         .and. LANTHRO   )
      LEDGAR            = ( LEDGAR         .and. LANTHRO   )
      LSTREETS          = ( LSTREETS       .and. LANTHRO   )
      LCAC              = ( LCAC           .and. LANTHRO   )
      LNEI05            = ( LNEI05         .and. LANTHRO   )
      LRETRO            = ( LRETRO         .and. LANTHRO   )
      LNEI99            = ( LNEI99         .and. LANTHRO   )
      LICARTT           = ( LICARTT        .and. LANTHRO   )
      LVISTAS           = ( LVISTAS        .and. LANTHRO   )
      LEDGARSHIP        = ( LEDGARSHIP     .and. LANTHRO   )
      LICOADSSHIP       = ( LICOADSSHIP    .and. LANTHRO   )
      LEMEPSHIP         = ( LEMEPSHIP      .and. LANTHRO   )
      LSHIPSO2          = ( LSHIPSO2       .and. LANTHRO   )
      LARCSHIP          = ( LARCSHIP       .and. LANTHRO   )
      LCOOKE            = ( LCOOKE         .and. LANTHRO   )

      ! Biogenic
      LMEGAN            = ( LMEGAN         .and. LBIOGENIC )
      LPECCA            = ( LPECCA         .and. LBIOGENIC )
      LMEGANMONO        = ( LMEGANMONO     .and. LBIOGENIC )
      LWARWICK_VSLS     = ( LWARWICK_VSLS  .and. LBIOGENIC )
      LSSABr2           = ( LSSABr2        .and. LBIOGENIC )

      ! Biomass
      LTOMSAI           = ( LTOMSAI        .and. LBIOMASS  )
      LGFED2BB          = ( LGFED2BB       .and. LBIOMASS  )
      L8DAYBB           = ( L8DAYBB        .and. LBIOMASS  )
      L3HRBB            = ( L3HRBB         .and. LBIOMASS  )
      LSYNOPBB          = ( LSYNOPBB       .and. LBIOMASS  )
      LGFED3BB          = ( LGFED3BB       .and. LBIOMASS  )
      LDAYBB3           = ( LDAYBB3        .and. LBIOMASS  )
      L3HRBB3           = ( L3HRBB3        .and. LBIOMASS  )

      ! Other
      LLIGHTNOX         = ( LLIGHTNOX      .and. LEMIS     )
      LOTDLOC           = ( LOTDLOC        .and. LEMIS     )
      LSOILNOX          = ( LSOILNOX       .and. LEMIS     )
      LFERTILIZERNOX    = ( LFERTILIZERNOX .and. LEMIS     )

      ! Define these flags for backwards compatibility
      LFOSSIL           = LANTHRO
      LWOODCO           = LBIOFUEL
      LBIONOX           = LBIOMASS
      
      ! Turn off full-chem only switches 
      IF ( .not. ITS_A_FULLCHEM_SIM() ) THEN
         LLIGHTNOX      = .FALSE.
         LAEIC          = .FALSE.
         LSOILNOX       = .FALSE.
         LFERTILIZERNOX = .FALSE.  
      ENDIF
      
      ! Set other EDGAR switches (for now set all together)
      IF ( LEDGAR ) THEN
         LEDGARNOx      = .TRUE.
         LEDGARCO       = .TRUE.
         LEDGARSOx      = .TRUE.
      ENDIF

      ! If using historical emissions, make sure EDGAR NOx/SO2 are used
      IF ( LHIST ) THEN
         LEDGARNOx = .TRUE.
         LEDGARSOx = .TRUE.
      ENDIF
      
      ! RCP: Test for conflicts with other switches (cdh, 6/14/12) 
      IF ( LRCP ) THEN
         
         ! Force user to choose just 1 global inventory
         IF ( LEDGAR ) 
     &        CALL ERROR_STOP( 'INPUT.GEOS WARNING!! '//
     &           'Choose either RCP or EDGAR for global emissions', 
     &           'READ_EMISSIONS_MENU' )

         ! Force user to choose just 1 global inventory
         IF ( LRETRO ) 
     &        CALL ERROR_STOP( 'INPUT.GEOS WARNING!! '//
     &           'Choose either RCP or RETRO for global emissions', 
     &           'READ_EMISSIONS_MENU' )

         ! Force user to choose just 1 global ship inventory
         IF ( LEDGARSHIP .OR. LICOADSSHIP ) 
     &        CALL ERROR_STOP( 'INPUT.GEOS WARNING!! '//
     &           'Choose only 1 global ship emissions: '//
     &           'RCP, EDGAR, or ICOADS',
     &           'READ_EMISSIONS_MENU' )

         ! Force user to choose just 1 global aircraft inventory
         IF ( LAEIC )
     &        CALL ERROR_STOP( 'INPUT.GEOS WARNING!! '//
     &           'Choose either RCP or AEIC for global aircraft' //
     &           ' emissions', 'READ_EMISSIONS_MENU' )

         ! Avoid double fertilizer NOx emissions
         IF ( LFERTILIZERNOX ) THEN
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            WRITE( 6, '(a)' ) 'WARNING! RCP inventory already includes,'
            WRITE( 6, '(a)' ) '~0.7 Tg(N) from agriculture for y2010.'
            WRITE( 6, '(a)' ) 'Turning off "Fertilizer NOx"'
            WRITE( 6, '(a)' ) 'to avoid duplicate emissions!'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            LFERTILIZERNOX = .FALSE.
         ENDIF

      ENDIF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% Bug fix!  If LEMEPSHIP is turned on but LEMEP is turned %%%
      !%%% off, this will cause an error (because arrays are not   %%%
      !%%% allocated, etc.).  For now, just turn off LEMEPSHIP     %%%
      !%%% and print a warning message.  Whoever wants to fix this %%%
      !%%% in a more robust way is welcome to do so.               %%%
      !%%% (mak, bmy, 10/19/09)                                    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( LEMEPSHIP .and. ( .not. LEMEP ) ) THEN
         LEMEPSHIP = .FALSE.
         IF( am_I_Root ) THEN
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            WRITE( 6, '(a)' ) 'WARNING! EMEP emissions are turned off,'
            WRITE( 6, '(a)' ) 'so also turn off EMEP ship emissions'
            WRITE( 6, '(a)' ) 'in order to avoid crashing GEOS-Chem!'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
	 ENDIF
      ENDIF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% Bug fix!  If LMEGANMONO is turned on but LMEGAN is turned %%%
      !%%% off, this will cause an error (because arrays are not     %%%
      !%%% allocated, etc.).  For now, just turn off LMEGANMONO      %%%
      !%%% and print a warning message.  Whoever wants to fix this   %%%
      !%%% in a more robust way is welcome to do so.                 %%%
      !%%% (bmy, 10/19/09)                                           %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( LMEGANMONO .and. ( .not. LMEGAN ) ) THEN
         LMEGANMONO = .FALSE.
         IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'WARNING! MEGAN emissions are turned off, so'
         WRITE( 6, '(a)' ) 'we also turn off LMEGANMONO ship emissions'
         WRITE( 6, '(a)' ) 'in order to avoid crashing GEOS-Chem!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      !=================================================================
      ! Enforce emissions options for Nested China on native GEOS-5 grid
      !=================================================================
#if   defined( NESTED_CH )

      ! Turn off emissions options that are
      ! not appropriate for the Nested CH grid
      LEMEP    = .FALSE.
      LBRAVO   = .FALSE.
      LCAC     = .FALSE.
      LNEI99   = .FALSE.
      LVISTAS  = .FALSE.
      LICARTT  = .FALSE.

#endif

#if   defined( NESTED_NA )

      ! Turn off emissions options that are
      ! not appropriate for the Nested NA grid
      LSTREETS = .FALSE.
      LEMEP    = .FALSE.

#endif

#if   defined( NESTED_EU )

      ! Turn off emissions options that are
      ! not appropriate for the Nested EU grid
      LSTREETS = .FALSE.
      LEDGAR   = .FALSE.
      LBRAVO   = .FALSE.
      LCAC     = .FALSE.
      LNEI99   = .FALSE.
      LVISTAS  = .FALSE.
      LICARTT  = .FALSE.

#endif

      
      !=================================================================
      ! Check SO2 ship emissions options
      !=================================================================
      IF ( LARCSHIP ) LSHIPSO2 = .FALSE. ! SO2-Arctas precedes SO2-Corbett
   
      !=================================================================
      ! Check global ship emissions options
      !=================================================================
      ! Add an ship emissions options (cklee, 6/30/09)
      ! Replace with ICOADS ship emissions
      ! User can still use regional EMEP, and ARCTAS-SO2 if he wants (phs) 
      IF ( LICOADSSHIP ) THEN
         LEDGARSHIP = .FALSE. ! ICOADS precedes EDGAR
         LSHIPSO2   = .FALSE. ! SO2-ICOADS precedes SO2-Corbett
      ENDIF      

      !=================================================================
      ! Check EPA options
      !=================================================================
      
      ! VISTAS and ICARTT assumes that EPA/NEI 99 is ON
      IF ( ( ( LVISTAS ) .OR. ( LICARTT )) .AND. .NOT.( LNEI99 ) ) THEN
         LNEI99 = .TRUE.
         IF ( am_I_Root ) THEN
          WRITE( 6, '(/,a,/)' ) 
     &       ' EPA/NEI99 has been automatically ' //
     &       'switched ON to use V.I.S.T.A.S. or/and ICARTT-based' //
     &       ' modifications'
       ENDIF
      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     ICARTT corrections are NOT available at high resolution
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
#if !defined(GRID2x25) && !defined(GRID4x5)
      
      IF ( LICARTT ) THEN
         LICARTT = .FALSE.
         IF ( am_I_Root ) THEN
            WRITE( 6, '(/,a,/)' ) 
     &         ' ICARTT-based corrections to  ' //
     &         'EPA-NEI99 are not available at high resolution.'
         ENDIF
      ENDIF
         
#endif

      ! NEI 2005 precedes NEI 1999
      IF ( LNEI05 .AND. LNEI99 ) THEN
         LNEI99  = .FALSE.
         LVISTAS = .FALSE.
         LICARTT = .FALSE.
         IF ( am_I_Root ) THEN
          WRITE( 6, '(/,a,/)' ) 
     &       ' EPA/NEI99 (incl. VISTAS and ICARTT' //
     &        'options) has been automatically ' //
     &        'switched OFF to use NEI 2005.'
         ENDIF
      ENDIF
         
      
      !=================================================================
      ! Error check GFED switches (mpayer, 11/7/11)
      !=================================================================
      IF ( ( LGFED2BB .or. L8DAYBB  .or. LSYNOPBB  .or. L3HRBB  )  .AND.
     &     ( LGFED3BB .or. LDAYBB3  .or. L3HRBB3   ) )  THEN

         MSG = 'GFED2 and GFED3 emissions cannot both be turned on.'
         CALL ERROR_STOP( MSG, LOC )
      ENDIF

      IF ( LDAYBB3 ) THEN
         IF ( .not. LGFED3BB) THEN
            MSG = 'To use daily GFED3, also select monthly GFED3.'
            CALL ERROR_STOP( MSG, LOC )
         ENDIF
      ENDIF

      IF ( L3HRBB3 ) THEN
         IF ( .not. LGFED3BB) THEN
            MSG = 'To use 3-hourly GFED3, also select monthly GFED3.'
            CALL ERROR_STOP( MSG, LOC )
         ENDIF
         IF ( .not. LDAYBB3) THEN
            MSG = 'To use 3-hourly GFED3, also select daily GFED3.'
            CALL ERROR_STOP( MSG, LOC )
         ENDIF
      ENDIF

      !=================================================================
      ! Prioritize GFED2
      !=================================================================
      IF ( L3HRBB ) THEN
         LGFED2BB = .FALSE.
         L8DAYBB  = .FALSE.
         LSYNOPBB = .FALSE.
      ELSE IF ( LSYNOPBB ) THEN
         LGFED2BB = .FALSE.
         L8DAYBB  = .FALSE.
      ELSE IF ( L8DAYBB ) THEN
         LGFED2BB = .FALSE.
      ENDIF

      !=================================================================
      ! Error checks
      !=================================================================

#if    defined( GCAP )

      !--------------------------------
      ! GCAP error message 
      !--------------------------------

      IF ( LOTDLOC ) THEN
            MSG = 'Lightning rescaling not available nor very ' //
     &            'appropriate for GCAP sim because of window ' // 
     &            'of OTD/LIS satellite observations.  Select ' //
     &            'one of the raw params without redist.'
         CALL ERROR_STOP( MSG, LOC )
       ENDIF

#endif

      !=================================================================
      ! Error check Bromine switches
      !=================================================================
      ! jpp, 5/10/10: if we're using the 1ppt daytime PBL BrO option,
      !              then shut off the LSSABr2 emissions if they were
      !              accidentally left on.
      IF ( LFIX_PBL_BRO .and. LSSABr2 ) THEN
         LSSABr2 = .FALSE. ! turn off sea-salt emissions to avoid conflict.
         IF ( am_I_Root ) THEN
            write(6,'(/,a,/)') 'Warning in input_mod.f: LSSABr2 was'
     &           // ' left ON when LFIX_PBL_BRO was also selected.'
            write(6,'(a)') 'Turning off sea-salt emissions and'
     &           // ' selecting the 1pptv PBL daytime BrO simulation.'
         ENDIF

      ENDIF

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LEMIS          = LEMIS
      Input_Opt%TS_EMIS        = TS_EMIS
      Input_Opt%LANTHRO        = LANTHRO
      Input_Opt%FSCALYR        = FSCALYR
      Input_Opt%LEMEP          = LEMEP 
      Input_Opt%LBRAVO         = LBRAVO
      Input_Opt%LEDGAR         = LEDGAR
      Input_Opt%LSTREETS       = LSTREETS
      Input_Opt%LCAC           = LCAC
      Input_Opt%LNEI05         = LNEI05
      Input_Opt%LRETRO         = LRETRO
      Input_Opt%LNEI99         = LNEI99
      Input_Opt%LICARTT        = LICARTT
      Input_Opt%LVISTAS        = LVISTAS
      Input_Opt%LBIOFUEL       = LBIOFUEL
      Input_Opt%LBIOGENIC      = LBIOGENIC
      Input_Opt%LMEGAN         = LMEGAN
      Input_Opt%LPECCA         = LPECCA 
      Input_Opt%LMEGANMONO     = LMEGANMONO
      Input_Opt%ISOP_SCALING   = ISOP_SCALING 
      Input_Opt%LBIOMASS       = LBIOMASS 
      Input_Opt%LBBSEA         = LBBSEA
      Input_Opt%LTOMSAI        = LTOMSAI
      Input_Opt%LGFED2BB       = LGFED2BB
      Input_Opt%L8DAYBB        = L8DAYBB
      Input_Opt%L3HRBB         = L3HRBB
      Input_Opt%LSYNOPBB       = LSYNOPBB 
      Input_Opt%LGFED3BB       = LGFED3BB
      Input_Opt%LDAYBB3        = LDAYBB3
      Input_Opt%L3HRBB3        = L3HRBB3
      Input_Opt%LAEIC          = LAEIC
      Input_Opt%LLIGHTNOX      = LLIGHTNOX
      Input_Opt%LOTDLOC        = LOTDLOC
      Input_Opt%LSOILNOX       = LSOILNOX
      Input_Opt%SOIL_RST_FILE  = SOIL_RST_FILE
      Input_Opt%LFERTILIZERNOX = LFERTILIZERNOX
      Input_Opt%NOx_SCALING    = NOx_SCALING
      Input_Opt%LEDGARSHIP     = LEDGARSHIP
      Input_Opt%LICOADSSHIP    = LICOADSSHIP
      Input_Opt%LEMEPSHIP      = LEMEPSHIP
      Input_Opt%LSHIPSO2       = LSHIPSO2
      Input_Opt%LARCSHIP       = LARCSHIP
      Input_Opt%LCOOKE         = LCOOKE
      Input_Opt%LHIST          = LHIST
      Input_Opt%HISTYR         = HISTYR
      Input_Opt%LWARWICK_VSLS  = LWARWICK_VSLS
      Input_Opt%LSSABr2        = LSSABr2
      Input_Opt%LFIX_PBL_BRO   = LFIX_PBL_BRO
      Input_Opt%Br_SCALING     = Br_SCALING
      Input_Opt%LEDGARNOx      = LEDGARNOx
      Input_Opt%LEDGARCO       = LEDGARCO
      Input_Opt%LEDGARSOx      = LEDGARSOx
      Input_Opt%LRCP           = LRCP
      Input_Opt%LRCPSHIP       = LRCPSHIP
      Input_Opt%LRCPAIR        = LRCPAIR

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Initialization
      !=================================================================

      CALL INIT_MODIS_LAI()

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
        WRITE( 6, '(/,a)' ) 'EMISSIONS MENU'
        WRITE( 6, '(  a)' ) '--------------'
        WRITE( 6, 100 ) 'Turn on emissions?          : ', LEMIS
        WRITE( 6, 110 ) 'Emissions timestep [min]    : ', TS_EMIS
        WRITE( 6, 100 ) 'Turn on ANTHRO emissions    : ', LANTHRO
        WRITE( 6, 110 ) '  ANTHRO scale factor year  : ', FSCALYR
        WRITE( 6, 100 ) '  Use EDGAR anthro emissions: ', LEDGAR
        WRITE( 6, 100 ) '  Use EMEP anthro emissions : ', LEMEP
        WRITE( 6, 100 ) '  Use BRAVO anthro emissions: ', LBRAVO
        WRITE( 6, 100 ) '  Use CAC anthro emissions  : ', LCAC
        WRITE( 6, 100 ) '  Use STREETS anthro emiss. : ', LSTREETS
        WRITE( 6, 100 ) '  Use EPA/NEI 2005          : ', LNEI05
        WRITE( 6, 100 ) '  Use RETRO emissions       : ', LRETRO
        WRITE( 6, 100 ) '  Use AEIC emissions        : ', LAEIC
        WRITE( 6, 100 ) '  Use RCP emissions         : ', LRCP
        WRITE( 6, 130 ) '      RCP scenario?         : ', TRIM(RCPNAME)
        WRITE( 6, 110 ) '      RCP year?             : ', RCPYEAR
        WRITE( 6, 100 ) 'Use EPA/NEI99 (ANTH + BF)   : ', LNEI99
        WRITE( 6, 100 ) '      --> with NOx VISTAS?  : ', LVISTAS
        WRITE( 6, 100 ) '      --> with ICARTT modif?: ', LICARTT
        WRITE( 6, 100 ) 'Turn on BIOFUEL emissions?  : ', LFOSSIL
        WRITE( 6, 100 ) 'Turn on BIOGENIC emissions? : ', LBIOGENIC
        WRITE( 6, 100 ) 'Use MEGAN biogenic emissions: ', LMEGAN
        WRITE( 6, 100 ) 'Use PCEEA BVOC model        : ', LPECCA
        WRITE( 6, 100 ) 'Use MEGAN bio emissions MONO: ', LMEGANMONO
        WRITE( 6, 120 ) ' Use isop scaling factor    : ', ISOP_SCALING
        WRITE( 6, 100 ) 'Turn on BIOMASS EMISSIONS   : ', LBIOMASS
        WRITE( 6, 100 ) 'Use seasonal BIOMASS emiss? : ', LBBSEA
        WRITE( 6, 100 ) 'Scale BIOMASS to TOMS-AI?   : ', LTOMSAI
        WRITE( 6, 100 ) 'Use GFED2 BIOMASS emissions?: ',
     &       LGFED2BB .or.  L8DAYBB .or. LSYNOPBB .or. L3HRBB 
        WRITE( 6, 100 ) '    monthly GFED2 emissions?: ', LGFED2BB
        WRITE( 6, 100 ) '    8-day GFED2 emission?   : ', L8DAYBB
        WRITE( 6, 100 ) '    3hr GFED2 emission?     : ', L3HRBB 
        WRITE( 6, 100 ) '    synoptic GFED2 ?        : ', LSYNOPBB
        WRITE( 6, 100 ) 'Use GFED3 BIOMASS emissions?: ',
     &       LGFED3BB .or.  LDAYBB3 .or. L3HRBB3
        WRITE( 6, 100 ) '    monthly GFED3 emissions?: ', LGFED3BB
        WRITE( 6, 100 ) '    daily GFED3 emission?   : ', LDAYBB3
        WRITE( 6, 100 ) '    3hr GFED3 emission?     : ', L3HRBB3
        WRITE( 6, 100 ) 'Use AIRCRAFT NOx from RCP?  : ', LRCPAIR
        WRITE( 6, 100 ) 'Turn on LIGHTNING NOx?      : ', LLIGHTNOX
        WRITE( 6, 100 ) 'Use OTD-LIS local redist?   : ', LOTDLOC
        WRITE( 6, 100 ) 'Turn on SOIL NOx?           : ', LSOILNOX
        WRITE( 6, 130 ) 'Soil NOx restart file name? : ', 
     &       TRIM( SOIL_RST_FILE )
        WRITE( 6, 100 ) 'Turn on Fertilizer NOx?     : ', 
     &                                                  LFERTILIZERNOX
        WRITE( 6, 120 ) 'NOx scaling                 : ', NOx_SCALING
        WRITE( 6, 100 ) 'Turn on EDGAR   SHIP emiss.?: ', LEDGARSHIP
        WRITE( 6, 100 ) 'Turn on ICOADS  SHIP emiss.?: ', LICOADSSHIP
        WRITE( 6, 100 ) 'Turn on RCP     SHIP emiss.?: ', LRCPSHIP
        WRITE( 6, 100 ) 'Turn on  EMEP   SHIP emiss.?: ', LEMEPSHIP
        WRITE( 6, 100 ) 'Turn on Corbett SHIP SO2 ?  : ', LSHIPSO2
        WRITE( 6, 100 ) '     or ARCTAS  SHIP SO2 ?  : ', LARCSHIP
        WRITE( 6, 100 ) 'Use COOKE for OC/BC N.-Amer.: ', LCOOKE
        WRITE( 6, 100 ) 'Turn on historical emiss?   : ', LHIST
        WRITE( 6, 100 ) 'Turn on the Warwick VSLs?   : ', LWARWICK_VSLS
        WRITE( 6, 100 ) 'Use seasalt Br2 emissions?  : ', LSSABr2
        WRITE( 6, 100 ) 'Set day PBL BrO conc to 1ppt: ', LFIX_PBL_BRO
        WRITE( 6, 100 ) 'Bromine emiss scale factor  : ', Br_SCALING
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
      ! add formatting for ISOP_SCALING and NOx_scaling (hotp 8/4/09)
 120  FORMAT( A, f6.2 )
 130  FORMAT( A, A   )

      END SUBROUTINE READ_EMISSIONS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_co2_sim_menu
!
! !DESCRIPTION: Subroutine READ\_CO2\_SIM\_MENU reads the CO2 SIM MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CO2_SIM_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE LOGICAL_MOD
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE TRACER_MOD,  ONLY : ITS_A_CO2_SIM

      USE CMN_SIZE_MOD                   ! Size parameters
      USE CMN_O3_MOD                     ! FSCALYR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  02 Mar 2009 - R. Nassar   - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  07 Sep 2011 - P. Kasibhatla - Modified for GFED3
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOC

      !=================================================================
      ! READ_CO2_SIM_MENU begins here!
      !=================================================================

      ! Location for error messages
      LOC = 'READ_CO2_SIM_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOC )
      ENDIF

      !---------------------------------------------------------------
      ! Fossil Fuel Emissions
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:1' )

      ! Use Generic Fossil Fuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:2' )
      READ( SUBSTRS(1:N), * ) LGENFF

      ! Use annual Fossil Fuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:3' )
      READ( SUBSTRS(1:N), * ) LANNFF

      ! Use monthly Fossil Fuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:4' )
      READ( SUBSTRS(1:N), * ) LMONFF

      !---------------------------------------------------------------
      ! Chemical source of CO2 from oxidation of fossil fuel CO, CH4 etc.
      !---------------------------------------------------------------
      ! Chemical source of CO2 from oxidation of CO, CH4 and organics?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:5' )
      READ( SUBSTRS(1:N), * ) LCHEMCO2

      !---------------------------------------------------------------
      ! Biomass and Biofuel Emissions
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:6' )

      ! Use seasonal only biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:7' )
      READ( SUBSTRS(1:N), * ) LSEASBB

      ! Use GFED2 monthly biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:8' )
      READ( SUBSTRS(1:N), * ) LGFED2BB

      ! Use GFED2 8-day biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:9' )
      READ( SUBSTRS(1:N), * ) L8DAYBB

      ! Use GFED3 monthly biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:10' )
      READ( SUBSTRS(1:N), * ) LGFED3BB

      ! Use GFED3 daily biomass emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:11' )
      READ( SUBSTRS(1:N), * ) LDAYBB3
   
      ! Include biofuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:12' )
      READ( SUBSTRS(1:N), * ) LBIOFUEL

      !---------------------------------------------------------------
      ! Terrestrial Exchange
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:13' )

      ! Turn on CASA biosphere daily average?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:14' )
      READ( SUBSTRS(1:N), * ) LBIODAILY

      ! Turn on CASA biosphere with diurnal cycle?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:15' )
      READ( SUBSTRS(1:N), * ) LBIODIURNAL

      ! Use original Net Terrestrial Exchange?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:16' )
      READ( SUBSTRS(1:N), * ) LBIONETORIG

      ! Use Net Terrestrial Exchange Climatology?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:17' )
      READ( SUBSTRS(1:N), * ) LBIONETCLIM

      !---------------------------------------------------------------
      ! Ocean Exchange
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:18' )

      ! Use Ocean Exchange from Takahashi 1997?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:19' )
      READ( SUBSTRS(1:N), * ) LOCN1997

      ! Use annual Ocean Exchange from Takahashi 2009?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:20' )
      READ( SUBSTRS(1:N), * ) LOCN2009ANN

      ! Use monthly Ocean Exchange from Takahashi 2009?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:21' )
      READ( SUBSTRS(1:N), * ) LOCN2009MON

      !---------------------------------------------------------------
      ! Ship and Aircraft CO2 Emissions
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:22' )

      ! Turn on EDGAR ship emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:23' )
      READ( SUBSTRS(1:N), * ) LSHIPEDG

      ! Turn on ICOADS ship emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:24' )
      READ( SUBSTRS(1:N), * ) LSHIPICO

      ! Turn on Aircraft emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:25' )
      READ( SUBSTRS(1:N), * ) LPLANE

      !---------------------------------------------------------------
      ! Tagged CO2
      !---------------------------------------------------------------
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:26' )

      !Background CO2 (no emissions or exchange) for Tagged-CO2 runs
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:27' )
      READ( SUBSTRS(1:N), * ) LFFBKGRD

      ! Turn on biosphere and ocean exchange region tagged tracers?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:28' )
      READ( SUBSTRS(1:N), * ) LBIOSPHTAG

      ! Turn on fossil fuel emission region tagged tracers?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:27' )
      READ( SUBSTRS(1:N), * ) LFOSSILTAG

      ! Turn on global ship emissions tagged tracer?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:29' )
      READ( SUBSTRS(1:N), * ) LSHIPTAG

      ! Turn on global aircraft emissions tagged tracer?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CO2_sim_menu:30' )
      READ( SUBSTRS(1:N), * ) LPLANETAG

      !=================================================================
      ! Resolve conflicting selections by assuming the 'best' choice
      !=================================================================

	!--- Fossil Fuel ----------------
	IF (LMONFF) THEN
           LGENFF = .FALSE.
           LANNFF = .FALSE.
	ENDIF
	IF (LANNFF) THEN
           LGENFF = .FALSE.
	ENDIF

	!--- Biomass Burning ------------
	IF (L8DAYBB) THEN
           LGFED2BB = .FALSE.
           LSEASBB  = .FALSE.
	ENDIF
	IF (LGFED2BB) THEN
           LSEASBB = .FALSE.
	ENDIF

	!--- CASA Biosphere -------------
	IF (LBIODIURNAL) THEN
           LBIODAILY = .FALSE.
	ENDIF

	!--- Net Terrestrial Exchange ---
	IF (LBIONETCLIM) THEN
           LBIONETORIG = .FALSE.
	ENDIF
	
	!--- Ocean Exchange -------------
	IF (LOCN2009MON) THEN
           LOCN2009ANN = .FALSE.
           LOCN1997 = .FALSE.
	ENDIF
	IF (LOCN2009ANN) THEN
		LOCN1997 = .FALSE.
	ENDIF

	!--- Ship Emissions -------------
	IF (LSHIPICO) THEN
           LSHIPEDG = .FALSE.
	ENDIF

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LGENFF      = LGENFF
      Input_Opt%LANNFF      = LANNFF
      Input_Opt%LMONFF      = LMONFF
      Input_Opt%LCHEMCO2    = LCHEMCO2
      Input_Opt%LSEASBB     = LSEASBB
      Input_Opt%LGFED2BB    = LGFED2BB
      Input_Opt%L8DAYBB     = L8DAYBB 
      Input_Opt%LGFED3BB    = LGFED3BB
      Input_Opt%LDAYBB3     = LDAYBB3
      Input_Opt%LBIOFUEL    = LBIOFUEL
      Input_Opt%LBIODAILY   = LBIODAILY
      Input_Opt%LBIODIURNAL = LBIODIURNAL
      Input_Opt%LBIONETORIG = LBIONETORIG
      Input_Opt%LBIONETCLIM = LBIONETCLIM
      Input_Opt%LOCN1997    = LOCN1997
      Input_Opt%LOCN2009ANN = LOCN2009ANN
      Input_Opt%LOCN2009MON = LOCN2009MON
      Input_Opt%LSHIPEDG    = LSHIPEDG
      Input_Opt%LSHIPICO    = LSHIPICO
      Input_Opt%LPLANE      = LPLANE
      Input_Opt%LFFBKGRD    = LFFBKGRD
      Input_Opt%LBIOSPHTAG  = LBIOSPHTAG
      Input_Opt%LFOSSILTAG  = LFOSSILTAG
      Input_Opt%LSHIPTAG    = LSHIPTAG
      Input_Opt%LPLANETAG   = LPLANETAG

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF (ITS_A_CO2_SIM() ) THEN
        IF ( am_I_Root ) THEN
          WRITE(6,'(/,a)') 'CO2 SIMULATION MENU ' //
     &         '(overwrites any other settings related to CO2)'
          WRITE(6,'(  a)') '-------------------------------------'
	  WRITE(6, 90    ) 'Fossil Fuel emissions'
	  WRITE(6,100    ) '  Generic FF emissions        :', LGENFF
	  WRITE(6,100    ) '  Annual FF emissions         :', LANNFF
	  WRITE(6,100    ) '  Monthly FF emissions        :', LMONFF
	  WRITE(6,100    ) 'CO2 from oxidation (CO,CH4,..):', LCHEMCO2
	  WRITE(6, 90    ) 'Biomass Burning emissions'
	  WRITE(6,100    ) '  Seasonal only biomass       :', LSEASBB
	  WRITE(6,100    ) '  GFED2 monthly biomass       :', LGFED2BB
	  WRITE(6,100    ) '  GFED2 8-day biomass         :', L8DAYBB
	  WRITE(6,100    ) '  GFED3 monthly biomass       :', LGFED3BB
	  WRITE(6,100    ) '  GFED3 daily biomass         :', LDAYBB3
	  WRITE(6,100    ) 'Biofuel emissions             :', LBIOFUEL
	  WRITE(6, 90    ) 'Terrestrial Exchange'
	  WRITE(6,100    ) '  CASA Biosphere daily average:', LBIODAILY
	  WRITE(6,100    ) '  CASA Biosphere diurnal cycle:', LBIODIURNAL
	  WRITE(6,100    ) '  Net Terr Exch - Original    :', LBIONETORIG
	  WRITE(6,100    ) '  Net Terr Exch - Climatology :', LBIONETCLIM
	  WRITE(6,90     ) 'Ocean Exchange' 
	  WRITE(6,100    ) '  Takahashi 1997              :', LOCN1997 
	  WRITE(6,100    ) '  Takahashi 2009 annual       :', LOCN2009ANN 
	  WRITE(6,100    ) '  Takahashi 2009 monthly      :', LOCN2009MON
	  WRITE(6,90     ) 'Ship and Plane Emissions' 
	  WRITE(6,100    ) '  EDGAR ship emissions        :', LSHIPEDG 
	  WRITE(6,100    ) '  ICOADS ship emissions       :', LSHIPICO 
	  WRITE(6,100    ) '  Aviation emissions          :', LPLANE 
	  WRITE(6, 90    ) 'Tagged CO2 settings'
	  WRITE(6,100    ) '  Save Fossil CO2 in Bckgrnd  :', LFFBKGRD 
	  WRITE(6,100    ) '  Tag Biosphere/Ocean CO2     :', LBIOSPHTAG 
	  WRITE(6,100    ) '  Tag Fossil Fuel CO2         :', LFOSSILTAG 
	  WRITE(6,100    ) '  Tag Global Ship CO2         :', LSHIPTAG
	  WRITE(6,100    ) '  Tag Global Aviation CO2     :', LPLANETAG
          WRITE(6,'(  a)') '-------------------------------------'
	ENDIF
      ENDIF      

      ! FORMAT statements
  90  FORMAT( A )
 100  FORMAT( A, L5 )
 110  FORMAT( A, L5, A )

      END SUBROUTINE READ_CO2_SIM_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_future_menu
!
! !DESCRIPTION: Subroutine READ\_FUTURE\_MENU reads the FUTURE MENU section 
!  of the GEOS-Chem input file; this defines IPCC future emissions options.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_FUTURE_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD, ONLY : DO_FUTURE_EMISSIONS
      USE LOGICAL_MOD,          ONLY : LFUTURE
 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  01 Jun 2006 - S. Wu       - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      INTEGER            :: FUTURE_YEAR 
      CHARACTER(LEN=255) :: FUTURE_SCEN
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_FUTURE_MENU begins here!
      !=================================================================

      ! Use IPCC future emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_future_menu:1' )
      READ( SUBSTRS(1:N), * ) LFUTURE

      ! Future emission year
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_future_menu:2' )
      READ( SUBSTRS(1:N), * ) FUTURE_YEAR

      ! Future emission scenario
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_future_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) FUTURE_SCEN

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_future_menu:4' )

#if   !defined( GCAP )
      !### TEMPORARY KLUDGE: right now, future emissions are only defined
      !### for the GCAP met fields.  Set LFUTURE=F for other met fields
      !### for the time being.  We will implement the future emissions for
      !### other met fields at a later date. (swu, bmy, 6/1/06)
      LFUTURE = .FALSE. 
#endif

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LFUTURE     = LFUTURE
      Input_Opt%FUTURE_YEAR = FUTURE_YEAR
      Input_Opt%FUTURE_SCEN = FUTURE_SCEN

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'FUTURE MENU'
         WRITE( 6, '(  a)' ) '-----------'
         WRITE( 6, 100 ) 'Use IPCC future emissions   : ', LFUTURE
         WRITE( 6, 110 ) 'Future emissions for year   : ', FUTURE_YEAR 
         WRITE( 6, 120 ) 'Future emissions scenario   : ',  
     &      TRIM( FUTURE_SCEN )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, I4  )
 120  FORMAT( A, A   )
    
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize
      IF ( LFUTURE ) THEN
         CALL DO_FUTURE_EMISSIONS( FUTURE_YEAR, TRIM( FUTURE_SCEN ) )
      ENDIF

      END SUBROUTINE READ_FUTURE_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_chemistry_menu
!
! !DESCRIPTION: Subroutine READ\_CHEMISTRY\_MENU reads the CHEMISTRY MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CHEMISTRY_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ERROR_STOP 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE LOGICAL_MOD,        ONLY : LCHEM, LSCHEM, LLINOZ
      USE LOGICAL_MOD,        ONLY : LSVCSPEC, LKPP
      USE TIME_MOD,           ONLY : SET_CT_CHEM
      USE TRACER_MOD,         ONLY : N_TRACERS
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      USE TIME_MOD,           ONLY : Accept_External_Date_Time
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1) added optional test on KPPTRACER (phs, 6/17/09)
!  (2) Remove reference to obsolete embedded chemistry stuff in "CMN" 
!      (bmy, 2/25/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  10 Jun 2012 - L. Murray   - Move all strat chemistry switches here
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  06 Dec 2012 - R. Yantosca - Now get TS_CHEM from the ESMF environment
!                              when we are connecting to the GEOS-5 GCM
!  11 Dec 2012 - R. Yantosca - ACCEPT_DATE_TIME_FROM_ESMF has now been renamed
!                              to ACCEPT_EXTERNAL_DATE_TIME
!  22 May 2013 - M. Payer    - Now read in GAMMA_HO2. Recommended value is 0.2
!                              based on Jacon et al (2000) and Mao et al (2013).
!  22 Aug 2013 - R. Yantosca - Now read in path for species restart file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, SPEC_RST_FILE

      !=================================================================
      ! READ_CHEMISTRY_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_CHEMISTRY_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:1' )
      READ( SUBSTRS(1:N), * ) LCHEM

      ! Turn on stratospheric chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:2' )
      READ( SUBSTRS(1:N), * ) LSCHEM

      ! Use Linoz for stratospheric ozone? (Otherwise, Synoz is used)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:3' )
      READ( SUBSTRS(1:N), * ) LLINOZ

      ! Chemistry timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:4' )
      READ( SUBSTRS(1:N), * ) TS_CHEM

      ! Read and save CSPEC ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:5' )
      READ( SUBSTRS(1:N), * ) LSVCSPEC

      ! Read and save CSPEC ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) SPEC_RST_FILE

      ! Use KPP solver ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:7' )
      READ( SUBSTRS(1:N), * ) LKPP

      ! GAMMA HO2 ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:8' )
      READ( SUBSTRS(1:N), * ) Input_Opt%GAMMA_HO2

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_chemistry_menu:9' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
#if defined( ESMF_ ) && defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! If we are connecting to the ESMF interface, we need to take
      ! the chemistry timestep as defined in the ESMF resource file.
      ! (i.e. GEOSCHEMchem_GridComp_mod.rc) instead of those in 
      ! input.geos.  This is because the ESMF Clock object needs to be
      ! defined at the highest level (in the driver routine), before 
      ! input.geos is ever read.
      !
      ! Therefore, we will assign the chemistry timestep field (i.e.
      ! Input_Opt%TS_CHEM) in the Gridded Component module file
      ! GEOSCHEMchem_GridComp_Mod.F90 (i.e. two levels higher
      ! in the code).  We don'need to define those fields here, so
      ! we have bracketed this with an #ifdef. 
      !
      ! Also pass the chemistry timestep value directly to module
      ! GeosUtil/time_mod.F via subroutine ACCEPT_EXTERNAL_DATE_TIME.
      !
      ! Ideallly it would be good to have one single version of 
      ! "input.geos" that can double as the GEOS-Chem input file as 
      ! well as the ESMF resource.  We may be a way off from that.
      ! One potential solution is to use a Perl script to generate
      ! both the input.geos file and the ESMF resource file, so that
      ! both files are consistent with each other.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------
      
      ! Save to local variables for printing below
      TS_CHEM                 = Input_Opt%TS_CHEM

      ! Pass chemistry timestep to time_mod.F
      CALL Accept_External_Date_Time( am_I_Root     = am_I_Root, 
     &                                value_TS_CHEM = TS_CHEM,   
     &                                RC            = RC         )
     &                                 

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the standard GEOS-Chem is to initialize 
      ! the field Input_OPt%TS_CHEM with the value read in from 
      ! input.geos. (bmy, 12/6/12)
      !-----------------------------------------------------------------
      Input_Opt%TS_CHEM       = TS_CHEM

#endif
      Input_Opt%LCHEM         = LCHEM
      Input_Opt%LSCHEM        = LSCHEM
      Input_Opt%LLINOZ        = LLINOZ
      Input_Opt%LSVCSPEC      = LSVCSPEC
      Input_Opt%SPEC_RST_FILE = TRIM( SPEC_RST_FILE )
      Input_Opt%LKPP          = LKPP

      ! Return success
      RC                      = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CHEMISTRY MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Turn on chemistry?          : ', LCHEM
         WRITE( 6, 100     ) 'Use linear. strat. chem?    : ', LSCHEM
         WRITE( 6, 100     ) ' => Use Linoz for O3?       : ', LLINOZ
         WRITE( 6, 110     ) 'Chemistry timestep [min]    : ', TS_CHEM
         WRITE( 6, 100     ) 'Use CSPEC restart?          : ', LSVCSPEC
         WRITE( 6, 130     ) ' => CSPEC restart file?     : ', 
     &        TRIM( SPEC_RST_FILE )
         WRITE( 6, 100     ) 'Use solver coded by KPP?    : ', LKPP
         WRITE( 6, 100     ) 'GAMMA HO2                   : ',
     &      Input_Opt%GAMMA_HO2
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, I5  )
 120  FORMAT( A, 2I5 )
 130  FORMAT( A, A   )

      END SUBROUTINE READ_CHEMISTRY_MENU  
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_transport_menu
!
! !DESCRIPTION: Subroutine READ\_TRANSPORT\_MENU reads the TRANSPORT MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TRANSPORT_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE LOGICAL_MOD,        ONLY : LFILL
      USE LOGICAL_MOD,        ONLY : LTRAN
      USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,         ONLY : ITS_A_TAGOX_SIM
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      USE TIME_MOD,           ONLY : Accept_External_Date_Time
#endif
#if !defined( ESMF_ )
      USE TRANSPORT_MOD,      ONLY : SET_TRANSPORT
#endif

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now define MAX_DYN for 1 x 1.25 grid (bmy, 12/1/04)
!  (2 ) Update text in error message (bmy, 2/23/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Don't stop run if TS_DYN > MAX_DYN but transport is turned off
!        (cdh, bmy, 7/7/08)
!  (5 ) Set MAX_DYN for the 0.5 x 0.666 nested grid (yxw, dan, bmy, 11/6/08)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  10 Jun 2012 - L. Murray - Move strat to chemistry menu
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  06 Dec 2012 - R. Yantosca - Now get TS_DYN from the ESMF environment, if
!                              we are connecting to the GEOS-5 GCM
!  11 Dec 2012 - R. Yantosca - ACCEPT_DATE_TIME_FROM_ESMF has now been renamed
!                              to ACCEPT_EXTERNAL_DATE_TIME
!  03 Oct 2013 - M. Sulprizio- Removed obsolete option for flux correction. This
!                              was used for GEOS-3, which has been retired.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N, IORD, JORD, KORD, J1, KS, MAX_DYN
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_TRANSPORT_MENU begins here!
      !=================================================================

      ! Location for err msg
      LOCATION = 'READ_TRANSPORT_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Turn on transport?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:1' )
      READ( SUBSTRS(1:N), * ) LTRAN

      ! Fill negative values
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:3' )
      READ( SUBSTRS(1:N), * ) LFILL

      ! IORD, JORD, KORD
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'read_transport_menu:4' )
      READ( SUBSTRS(1:N), * ) IORD, JORD, KORD

      ! Transport timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:5' )
      READ( SUBSTRS(1:N), * ) TS_DYN

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_transport_menu:6' )

      !=================================================================
      ! Error check settings
      !=================================================================      
      
      ! Define maximum timestep for transport
#if   defined( GRID4x5   ) 
      MAX_DYN = 30
#elif defined( GRID2x25  )
      MAX_DYN = 15
#elif defined( GRID1x125 )
      MAX_DYN = 10
#elif defined( GRID1x1   ) 
      MAX_DYN = 10
#elif defined( GRID05x0666   )
      MAX_DYN = 10 
#elif defined( GRID025x03125 )
      MAX_DYN = 5
#elif defined( EXTERNAL_FORCING )
      MAX_DYN = 60
#endif

      ! If TS_DYN is greater than MAX_DYN, then stop w/ error
      IF ( TS_DYN > MAX_DYN .and. LTRAN ) THEN
         MSG = 'Transport timestep is too big!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LTRAN       = LTRAN
      Input_Opt%LFILL       = LFILL
      Input_Opt%TPCORE_IORD = IORD
      Input_Opt%TPCORE_JORD = JORD
      Input_Opt%TPCORE_KORD = KORD

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! If we are connecting to the ESMF interface, we need to take
      ! the dynamic timestep as defined in the ESMF resource file.
      ! (i.e. GEOSCHEMchem_GridComp_mod.rc) instead of that in 
      ! input.geos.  This is because the ESMF Clock object needs to be
      ! defined at the highest level (in the driver routine), before 
      ! input.geos is ever read.
      !
      ! Therefore, we will assign the chemistry timestep field (i.e.
      ! Input_Opt%TS_CHEM) in the Gridded Component module file
      ! GEOSCHEMchem_GridComp_Mod.F90 (i.e. two levels higher
      ! in the code).  We don'need to define those fields here, so
      ! we have bracketed this with an #ifdef. 
      !
      ! Also pass the chemistry timestep value directly to module
      ! GeosUtil/time_mod.F via subroutine ACCEPT_EXTERNAL_DATE_TIME.
      !
      ! Ideallly it would be good to have one single version of 
      ! "input.geos" that can double as the GEOS-Chem input file as 
      ! well as the ESMF resource.  We may be a way off from that.
      ! One potential solution is to use a Perl script to generate
      ! both the input.geos file and the ESMF resource file, so that
      ! both files are consistent with each other.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------
      
      ! Save to local variables for printing below
      TS_DYN                = Input_Opt%TS_DYN

      ! Pass chemistry timestep to time_mod.F
      CALL Accept_External_Date_Time( am_I_Root     = am_I_Root, 
     &                                value_TS_DYN  = TS_DYN,
     &                                RC            = RC         )
     &                                 

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the standard GEOS-Chem is to initialize 
      ! the field Input_Optt%TS_DYN with the value read in from 
      ! input.geos. (bmy, 12/6/12)
      !-----------------------------------------------------------------
      Input_Opt%TS_DYN      = TS_DYN

#endif

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'TRANSPORT MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Turn on transport?          : ', LTRAN
         WRITE( 6, 100     ) 'Let TPCORE Fill negatives?  : ', LFILL
         WRITE( 6, 110     ) 'IORD, JORD, KORD for TPCORE?: ', IORD, 
     &                                                   JORD, KORD
         WRITE( 6, 120     ) 'Transport timestep [min]    : ', TS_DYN
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, 5I5 )
 120  FORMAT( A, I5  )
      
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================
#if ! defined( ESMF_ )
      ! Pass parameters to "transport_mod.f"
      CALL SET_TRANSPORT( IORD, JORD, KORD )
#endif

      END SUBROUTINE READ_TRANSPORT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_convection_menu 
!
! !DESCRIPTION: Subroutine READ\_CONVECTION\_MENU reads the CONVECTION MENU 
!  section of the GEOS-Chem input file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CONVECTION_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE LOGICAL_MOD,        ONLY : LCONV, LTURB
      USE LOGICAL_MOD,        ONLY : LNLPBL        ! (Lin, 03/31/09)
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      USE TIME_MOD,           ONLY : Accept_External_Date_Time
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Add option for new non-local PBL scheme. And a check on GEOS-5, 
!        LNLPBL turned to false if GEOS-5 is not used (lin, ccc 5/13/09)
!  27 Aug 2010 - R. Yantosca - Now allow non-local PBL for MERRA met data
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  02 Feb 2012 - R. Yantosca - Added modifications for MERRA met data
!  13 Apr 2012 - R. Yantosca - Fixed typo ( defined( GEOS_FP ) should have 
!                              been !defined( GEOS_FP ) )
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  01 Mar 2013 - R. Yantosca - Now set TS_CONV to the same value as TS_DYN
!                              when connecting to the GEOS-5 GCM.         
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_CONVECTION_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_TRANSPORT_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on convection?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:1' )
      READ( SUBSTRS(1:N), * ) LCONV

      ! Turn on BL mixing
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:2' )
      READ( SUBSTRS(1:N), * ) LTURB

      ! Turn on non-local PBL scheme (Lin, 03/31/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:3' )
      READ( SUBSTRS(1:N), * ) LNLPBL

      ! Convection timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:4' )
      READ( SUBSTRS(1:N), * ) TS_CONV

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_convection_menu:5' )

      ! The non-local PBL scheme is valid only with GEOS-5/MERRA/GEOS-FP
#if   !defined( GEOS_5 ) && !defined( MERRA ) && !defined( GEOS_FP )
      IF ( LNLPBL ) THEN
         IF( am_I_Root ) THEN
            WRITE(6,*) REPEAT( '=', 79 )
            WRITE(6,*) 'The non-local PBL scheme is only valid for'
            WRITE(6,*) 'GEOS-5, GEOS-FP, and MERRA met fields!'
            WRITE(6,*) 'LNLPBL is automatically turned to false !'
            WRITE(6,*) REPEAT( '=', 79 )
         ENDIF
         LNLPBL = .FALSE.
      ENDIF
#endif

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LCONV   = LCONV
      Input_Opt%LTURB   = LTURB
      Input_Opt%LNLPBL  = LNLPBL

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! If we are connecting to the ESMF interface, we need to set
      ! the convection timestep equal to the dynamic timestep in
      ! order to avoid tripping an error message in the routine
      ! CHECK_TIMESTEPS. (bmy, 3/1/13)
      !-----------------------------------------------------------------
      
      ! Set TS_CONV from Input_Opt%TS_DYN
      TS_CONV           = Input_Opt%TS_DYN

      ! Pass chemistry timestep to time_mod.F
      CALL Accept_External_Date_Time( am_I_Root     = am_I_Root, 
     &                                value_TS_CONV = TS_CONV,
     &                                RC            = RC         )
     &                                 
#endif

      Input_Opt%TS_CONV = TS_CONV

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CONVECTION MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100     ) 'Turn on cloud convection?   : ', LCONV
         WRITE( 6, 100     ) 'Turn on PBL mixing?         : ', LTURB
         WRITE( 6, 100     ) 'Turn on non-local PBL?      : ', LNLPBL
         WRITE( 6, 110     ) 'Convection timestep [min]   : ', TS_CONV
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
      
      END SUBROUTINE READ_CONVECTION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_deposition_menu
!
! !DESCRIPTION: Subroutine READ\_DEPOSITION\_MENU reads the DEPOSITION MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_DEPOSITION_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE DRYDEP_MOD,  ONLY : INIT_DRYDEP
      USE LOGICAL_MOD, ONLY : LCONV,             LDRYD
      USE LOGICAL_MOD, ONLY : LWETD,             LSPLIT
      USE LOGICAL_MOD, ONLY : USE_OLSON_2001
      USE TRACER_MOD,  ONLY : ITS_A_C2H6_SIM,    ITS_A_CH3I_SIM
      USE TRACER_MOD,  ONLY : ITS_A_CH4_SIM,     ITS_A_HCN_SIM
      USE TRACER_MOD,  ONLY : ITS_A_MERCURY_SIM, ITS_A_TAGCO_SIM
      USE TRACER_MOD,  ONLY : ITS_A_TAGOX_SIM
      USE WETSCAV_MOD, ONLY : WETDEPID
#if   defined( APM )
      USE APM_WETS_MOD, ONLY : WETDEPBINID
#endif
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now print an informational message for tagged Hg (bmy, 12/15/04)
!  (2 ) We need to call WETDEPID for both wetdep and cloud convection
!        since this sets up the list of soluble tracers (bmy, 3/1/05)
!  (3 ) Remove references to obsolete CO_OH simulation (bmy, 6/24/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  31 Jul 2012 - R. Yantosca - Now pass am_I_Root to INIT_DRYDEP
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  26 Feb 2013 - R. Yantosca - Now call INIT_DUST jere to facilitate 
!                              connecting to the GEOS-5 GCM
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_DEPOSITION_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_DEPOSITION_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on drydep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:1' )
      READ( SUBSTRS(1:N), * ) LDRYD

      ! Turn on wetdep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:2' )
      READ( SUBSTRS(1:N), * ) LWETD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_deposition_menu:3' )

      !=================================================================
      ! Error check settings
      !=================================================================

      ! Turn off drydep for simulations that don't need it
      IF ( ITS_A_CH3I_SIM()    ) LDRYD = .FALSE.
      IF ( ITS_A_TAGCO_SIM()   ) LDRYD = .FALSE.

      ! Turn off wetdep for simulations that don't need it
      IF ( ITS_A_CH3I_SIM()    ) LWETD = .FALSE.
      IF ( ITS_A_HCN_SIM()     ) LWETD = .FALSE.
      IF ( ITS_A_TAGOX_SIM()   ) LWETD = .FALSE.
      IF ( ITS_A_TAGCO_SIM()   ) LWETD = .FALSE.
      IF ( ITS_A_C2H6_SIM()    ) LWETD = .FALSE.
      IF ( ITS_A_CH4_SIM()     ) LWETD = .FALSE.

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LDRYD = LDRYD
      Input_Opt%LWETD = LWETD

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'DEPOSITION MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, 100     ) 'Turn on dry deposition?     : ', LDRYD
         WRITE( 6, 100     ) 'Turn on wet deposition?     : ', LWETD
         WRITE( 6, 100     ) 'Use Olson 2001 land map?    : ', 
     &                        USE_OLSON_2001
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize dry deposition arrays
      IF ( Input_Opt%LDRYD ) THEN

         ! Setup for dry deposition
         CALL INIT_DRYDEP( am_I_Root, Input_Opt, RC )
      
         ! Print extra info message for Hg simulation
         IF ( ITS_A_MERCURY_SIM() .and. LSPLIT ) THEN
            WRITE ( 6, 120 )
            WRITE ( 6, 121 )
         ENDIF
      ENDIF

      ! Initialize wet deposition tracers
      IF ( LWETD .or. LCONV ) THEN
         CALL WETDEPID( am_I_Root, Input_Opt, RC )
      ENDIF

#if   defined( APM )
      ! Initialize wet deposition tracers for APM
      IF ( LWETD .or. LCONV ) CALL WETDEPBINID
#endif

      ! FORMAT strings
 120  FORMAT( /, 'All tagged Hg2 tracers have the same dep velocity '
     &           'as the total Hg2 tracer.' )
 121  FORMAT( 'All tagged HgP tracers have the same dep velocity '
     &        'as the total HgP tracer.' )

      END SUBROUTINE READ_DEPOSITION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gamap_menu
!
! !DESCRIPTION: Subroutine READ\_GAMAP\_MENU reads the GAMAP MENU section 
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_GAMAP_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  25 Apr 2005 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: EOF
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_GAMAP_MENU begins here!
      !=================================================================

      ! Background
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_gamap_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) DIAGINFO

      ! Redirect
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_gamap_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) TRACERINFO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_gamap_menu:3' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%GAMAP_DIAGINFO   = DIAGINFO
      Input_Opt%GAMAP_TRACERINFO = TRACERINFO

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'GAMAP MENU'
         WRITE( 6, '(  a)' ) '---------------'            
         WRITE( 6, '(a,a)' ) 'GAMAP "diaginfo.dat"   file : ', 
     &                        TRIM( DIAGINFO   )
         WRITE( 6, '(a,a)' ) 'GAMAP "tracerinfo.dat" file : ',
     &                        TRIM( TRACERINFO )
      ENDIF

      END SUBROUTINE READ_GAMAP_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_output_menu
!
! !DESCRIPTION: Subroutine READ\_OUTPUT\_MENU reads the OUTPUT MENU section of 
!  the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_OUTPUT_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
      
      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! NJDAY
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  03 Aug 2012 - R. Yantosca - IU_GEOS is now a global module variable
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  28 Feb 2013 - R. Yantosca - Don't call IS_LAST_DAY_GOOD when using ESMF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IOS

      !=================================================================
      ! READ_OUTPUT_MENU begins here!
      !=================================================================
      READ( IU_GEOS, 100, IOSTAT=IOS ) NJDAY
 100  FORMAT( 26x, 31i1, /  26x, 29i1, /, 26x, 31i1, /, 26x, 30i1, /, 
     &        26x, 31i1, /, 26x, 30i1, /, 26x, 31i1, /, 26x, 31i1, /,
     &        26x, 30i1, /  26x, 31i1, /, 26x, 30i1, /, 26x, 31i1 )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_output_menu:1' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%NJDAY = NJDAY
      
      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'OUTPUT MENU'
         WRITE( 6, '(  a)' ) '-----------'
         WRITE( 6, 110     )
         WRITE( 6, 120     )
         WRITE( 6, 130     )
         WRITE( 6, 140     ) NJDAY
      ENDIF

      ! FORMAT statements
 110  FORMAT( '              1111111111222222222233' )
 120  FORMAT( '     1234567890123456789012345678901' )
 130  FORMAT( '     -------------------------------' )
 140  FORMAT( 'JAN--', 31i1, /, 'FEB--', 29i1, /, 'MAR--', 31i1, /, 
     &        'APR--', 30i1, /, 'MAY--', 31i1, /, 'JUN--', 30i1, /, 
     &        'JUL--', 31i1, /, 'AUG--', 31i1, /, 'SEP--', 30i1, /,
     &        'OCT--', 31i1, /, 'NOV--', 30i1, /, 'DEC--', 31i1 )

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! When we connect to the GEOS-5 GCM via ESMF, we let the GCM
      ! handle all timing and diagnostic output via the GEOS-5 History
      ! Component.  Therefore we do not need to make sure that we have 
      ! scheduled output for the last day of the run.  This can cause
      ! a premature exit. (bmy, 2/28/13)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the std GEOS-Chem is to not let the run
      ! proceed unless the user has scheduled output for the last day.
      ! (bmy, 2/28/13)
      !-----------------------------------------------------------------

      ! Make sure we have output at end of run
      CALL IS_LAST_DAY_GOOD
#endif

      END SUBROUTINE READ_OUTPUT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_diagnostic_menu
!
! !DESCRIPTION: Subroutine READ\_DIAGNOSTIC\_MENU reads the DIAGNOSTIC MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_DIAGNOSTIC_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE BIOFUEL_MOD,  ONLY : NBFTRACE
      USE BPCH2_MOD,    ONLY : OPEN_BPCH2_FOR_WRITE
      USE DIAG03_MOD,   ONLY : ND03,      PD03,      INIT_DIAG03
      USE DIAG03_MOD,   ONLY : PD03_PL    !eds 9/9/10
      USE DIAG04_MOD,   ONLY : ND04,      PD04,      INIT_DIAG04
      USE DIAG41_MOD,   ONLY : ND41,      PD41,      INIT_DIAG41
      USE DIAG42_MOD,   ONLY : ND42,      PD42,      INIT_DIAG42
      USE DIAG53_MOD,   ONLY : ND53,      PD53,      INIT_DIAG53
      USE DIAG56_MOD,   ONLY : ND56,      PD56,      INIT_DIAG56
      USE DIAG_OH_MOD,  ONLY : INIT_DIAG_OH
      USE DRYDEP_MOD,   ONLY : NUMDEP
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_BPCH
      USE LOGICAL_MOD,  ONLY : LBIOMASS,  LBIOFUEL,  LCARB, LCONV    
      USE LOGICAL_MOD,  ONLY : LDRYD,     LDUST,     LPRT,  LSULF    
      USE LOGICAL_MOD,  ONLY : LSSALT,    LTURB,     LWETD, LGFED2BB  
      USE LOGICAL_MOD,  ONLY : LGFED3BB  
      USE TIME_MOD,     ONLY : GET_NYMDb, GET_NHMSb, EXPAND_DATE
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE TRACER_MOD,   ONLY : ITS_A_CO2_SIM,        ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,   ONLY : ITS_A_MERCURY_SIM,    ITS_A_RnPbBe_SIM
      USE TRACER_MOD,   ONLY : ITS_A_TAGOX_SIM,      ITS_A_CH3I_SIM
      USE TRACER_MOD,   ONLY : SALA_REDGE_um,        ITS_A_CH4_SIM
      USE TRACER_MOD,   ONLY : ITS_A_POPS_SIM
      USE TRACERID_MOD, ONLY : NEMANTHRO
      USE WETSCAV_MOD,  ONLY : GET_WETDEP_NMAX

      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_DIAG_MOD       ! NDxx flags
!-----------------------------------------------------------------------------
! Prior to 2/3/14:
! NOTE: These variables are in Headers/CMN_DIAG_mod.F, so we don't need
! to reference these.  Lee Murray says that this causes the code to choke
! when compiling GEOS-Chem on Mac. (bmy, 2/3/14)
!      USE DIAG_MOD,           ONLY: TINDEX
!      USE DIAG_MOD,           ONLY: TCOUNT
!      USE DIAG_MOD,           ONLY: TMAX
!-----------------------------------------------------------------------------
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now reference IU_BPCH from "file_mod.f" and OPEN_BPCH2_FOR_WRITE
!        from "bpch2_mod.f".  Now opens the bpch file for output here
!        instead of w/in "main.f" (bmy, 2/3/05)
!  (2 ) Now references "diag03_mod.f" and "diag41_mod.f".  Now turn off ND38
!        when both LWETD=F and LCONV=F.  Now calls EXPAND_DATE to replace
!        YYYYMMDD and HHMMSS tokens in the bpch file name with the actual
!        starting date & time of the run. (bmy, 3/25/05)
!  (3 ) Now get diag info for ND09 for HCN/CH3CN sim (bmy, 6/27/05)
!  (4 ) Now references "diag04_mod.f" (bmy, 7/26/05)
!  (5 ) Now make sure all USE statements are USE, ONLY.  Also remove reference
!        to DIAG_MOD, it's not needed. (bmy, 10/3/05)
!  (6 ) Now remove reference to NBIOTRCE; Replace w/ NBIOMAX. (bmy, 4/5/06)
!  (7 ) Now reference ND56, PD56, INIT_DIAG56 from "diag56_mod.f" 
!        (bmy, 5/10/06)
!  (8 ) Now reference ND42, PD42, INIT_DIAG42 from "diag42_mod.f"
!        (dkh, bmy, 5/22/06)
!  (9 ) Now set max dimension for GFED2 or default biomass (bmy, 9/22/06)
!  (10) Bug fix: Should use ND52 in call to SET_TINDEX (cdh, bmy, 2/11/08)
!  (11) Remove call to NDXX_SETUP; this is now called in READ_INPUT_FILE.
!        (phs, 11/18/08)
!  (12) Now set TINDEX with PD45=NNPAR+1 tracers instead of N_TRACERS.
!        (tmf, 2/10/09)
!  (13) NBIOMAX now in CMN_SIZE (fp, 6/2009)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  26 May 2011 - R. Yantosca - For ND17, ND18, ND37, ND38, ND39, we need to
!                              set N_TMP = N_TRACERS, or else wetdep tracers 
!                              with indices higher than #32 won't print out.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  08 Nov 2013 - M. Sulprizio- Remove HR1_NO, and HR2_NO from ND43 diagnostic.
!  03 Feb 2014 - R. Yantosca - Remove references to TINDEX, TCOUNT, and TMAX
!                              from diag_mod.  They are in CMN_SIZE_mod.F.
!                              Lee Murray reports this causes the compilation
!                              to choke on Macintosh platforms.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: M, N, N_MAX, N_TMP
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG, LOCATION

      !=================================================================
      ! READ_DIAGNOSTIC_MENU begins here!
      !=================================================================

      ! Location for ERROR_STOP
      LOCATION = 'READ_DIAGNOSTIC_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Get max number of tracers
      N_MAX = MIN( N_TRACERS, NNPAR )

      ! Binary punch file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:2' )

      !--------------------------
      ! ND01: Rn-Pb-Be source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:3' )
      READ( SUBSTRS(1), * ) ND01
      IF ( .not. ITS_A_RnPbBe_SIM() ) ND01 = 0 
      CALL SET_TINDEX( am_I_Root, 01, ND01, SUBSTRS(2:N), N-1, N_MAX )
      
      !--------------------------
      ! ND02: Rn-Pb-Be decay
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:4' )
      READ( SUBSTRS(1), * ) ND02
      IF ( .not. ITS_A_RnPbBe_SIM() ) ND02 = 0 
      CALL SET_TINDEX( am_I_Root, 02, ND02, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND03: Hg diagnostics
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:5' )
      READ( SUBSTRS(1), * ) ND03
      IF ( .not. ITS_A_MERCURY_SIM() ) ND03 = 0
      !CALL SET_TINDEX( am_I_Root, 03, ND03, SUBSTRS(2:N), N-1, PD03_PL ) !eds 9/9/10
      CALL SET_TINDEX( am_I_Root, 03, ND03, SUBSTRS(2:N), N-1, PD03 ) !yzh 12/15/12
	  

      !--------------------------
      ! ND04: CO2 emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:6' )
      READ( SUBSTRS(1), * ) ND04
      IF ( .not. ITS_A_CO2_SIM() ) ND04 = 0
      CALL SET_TINDEX( am_I_Root, 04, ND04, SUBSTRS(2:N), N-1, PD04 )

      !--------------------------
      ! ND05: Sulfate prod/loss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:7' )
      READ( SUBSTRS(1), * ) ND05
      IF ( .not. LSULF ) ND05 = 0
      CALL SET_TINDEX( am_I_Root, 05, ND05, SUBSTRS(2:N), N-1, PD05 )

      !--------------------------
      ! ND06: Dust emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:8' )
      READ( SUBSTRS(1), * ) ND06
      IF ( .not. LDUST ) ND06 = 0
      CALL SET_TINDEX( am_I_Root, 06, ND06, SUBSTRS(2:N), N-1, NDSTBIN )

      !--------------------------
      ! ND07: Carbon/SOA source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:9' )
      READ( SUBSTRS(1), * ) ND07
      IF ( .not. LCARB ) ND07 = 0
      CALL SET_TINDEX( am_I_Root, 07, ND07, SUBSTRS(2:N), N-1, PD07 )

      !--------------------------
      ! ND08: Sea salt source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:10' )
      READ( SUBSTRS(1), * ) ND08
      IF ( .not. LSSALT ) ND08 = 0
      CALL SET_TINDEX( am_I_Root, 08, ND08, SUBSTRS(2:N), N-1, PD08 )

      !--------------------------
      ! ND09: HCN/CH3CN
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:11' )
      READ( SUBSTRS(1), * ) ND09
      CALL SET_TINDEX( am_I_Root, 
     &     09, ND09, SUBSTRS(2:N), N-1, N_TRACERS+PD09 )

      !--------------------------
      ! ND10: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:12' )
      READ( SUBSTRS(1), * ) ND10
      CALL SET_TINDEX( am_I_Root, 10, ND10, SUBSTRS(2:N), N-1, PD10 )

      !--------------------------
      ! ND11: Acetone source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:13' )
      READ( SUBSTRS(1), * ) ND11
      CALL SET_TINDEX( am_I_Root, 11, ND11, SUBSTRS(2:N), N-1, PD11 )

      !--------------------------
      ! ND12: PBL distribution
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:14' )
      READ( SUBSTRS(1), * ) ND12
      CALL SET_TINDEX( am_I_Root, 12, ND12, SUBSTRS(2:N), N-1, PD12 )

      !--------------------------
      ! ND13: Sulfur sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:15' )
      READ( SUBSTRS(1), * ) ND13
      IF ( .not. LSULF ) ND13 = 0
      CALL SET_TINDEX( am_I_Root, 13, ND13, SUBSTRS(2:N), N-1, PD13 )

      !--------------------------
      ! ND14: Wet conv up flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:16' )
      READ( SUBSTRS(1), * ) ND14
      IF ( .not. LCONV ) ND14 = 0
      CALL SET_TINDEX( am_I_Root, 14, ND14, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND15: Mass change in PBL
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:17' )
      READ( SUBSTRS(1), * ) ND15
      IF ( .not. LTURB ) ND15 = 0
      CALL SET_TINDEX( am_I_Root, 15, ND15, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND16: Precip fractions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:18' )
      READ( SUBSTRS(1), * ) ND16
      CALL SET_TINDEX( am_I_Root, 16, ND16, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND17: Rainout losses
      !--------------------------
      N_TMP = N_TRACERS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:19' )
      READ( SUBSTRS(1), * ) ND17
      IF ( .not. LWETD ) ND17 = 0
      CALL SET_TINDEX( am_I_Root, 17, ND17, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND18: Washout losses
      !--------------------------
      N_TMP = N_TRACERS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:20' )
      READ( SUBSTRS(1), * ) ND18
      IF ( .not. LWETD ) ND18 = 0
      CALL SET_TINDEX( am_I_Root, 18, ND18, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND19: CH4 loss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:21' )
      READ( SUBSTRS(1), * ) ND19
      CALL SET_TINDEX( am_I_Root, 19, ND19, SUBSTRS(2:N), N-1, PD19 )

      !--------------------------
      ! ND21: Opt depths etc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:22' )
      READ( SUBSTRS(1), * ) ND21
      CALL SET_TINDEX( am_I_Root, 21, ND21, SUBSTRS(2:N), N-1, PD21 )

      ! Error check 
      IF ( ND21 > 0 .and. SALA_REDGE_um(2) /= 0.5 ) THEN
         MSG = 'Cannot output seasalt AOD''s when radius bin' //
     &         ' is not split at 0.5 um!!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      !--------------------------
      ! ND22: J-values
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:23' )
      READ( SUBSTRS(1), * ) ND22
      CALL SET_TINDEX( am_I_Root, 22, ND22, SUBSTRS(2:N), N-1, PD22 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:24' ) 
      READ( SUBSTRS(1:N), * ) HR1_JV, HR2_JV

      !--------------------------
      ! ND24: E/W transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:25' )
      READ( SUBSTRS(1), * ) ND24
      CALL SET_TINDEX( am_I_Root, 24, ND24, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND25: N/S transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:26' )
      READ( SUBSTRS(1), * ) ND25
      CALL SET_TINDEX( am_I_Root, 25, ND25, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND26: U/D transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:27' )
      READ( SUBSTRS(1), * ) ND26
      CALL SET_TINDEX( am_I_Root, 26, ND26, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND27: Strat fluxes
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:28' )
      READ( SUBSTRS(1), * ) ND27
      CALL SET_TINDEX( am_I_Root, 27, ND27, SUBSTRS(2:N), N-1, PD27 )

      !--------------------------
      ! ND28: Biomass emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:29' )
      READ( SUBSTRS(1), * ) ND28
      IF ( .not. LBIOMASS ) ND28 = 0
      CALL SET_TINDEX( am_I_Root, 28, ND28, SUBSTRS(2:N), N-1, NBIOMAX )

      !--------------------------
      ! ND29: CO sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:30' )
      READ( SUBSTRS(1), * ) ND29
      CALL SET_TINDEX( am_I_Root, 29, ND29, SUBSTRS(2:N), N-1, PD29 )

      !--------------------------
      ! ND30: Land map
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:31' )
      READ( SUBSTRS(1), * ) ND30
      CALL SET_TINDEX( am_I_Root, 30, ND30, SUBSTRS(2:N), N-1, PD30 )

      !--------------------------
      ! ND31: Surface pressure
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:32' )
      READ( SUBSTRS(1), * ) ND31
      CALL SET_TINDEX( am_I_Root, 31, ND31, SUBSTRS(2:N), N-1, PD31 )

      !--------------------------
      ! ND32: NOx sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:33' )
      READ( SUBSTRS(1), * ) ND32
      CALL SET_TINDEX( am_I_Root, 32, ND32, SUBSTRS(2:N), N-1, PD32 )

      !--------------------------
      ! ND33: Column tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:34' )
      READ( SUBSTRS(1), * ) ND33
      CALL SET_TINDEX( am_I_Root, 33, ND33, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND34: Biofuel sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:35' )
      READ( SUBSTRS(1), * ) ND34
      IF ( .not. LBIOFUEL ) ND34 = 0
      CALL SET_TINDEX( am_I_Root, 
     &                 34, ND34, SUBSTRS(2:N), N-1, NBFTRACE )

      !--------------------------
      ! ND35: 500 hPa tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:36' )
      READ( SUBSTRS(1), * ) ND35
      CALL SET_TINDEX( am_I_Root, 35, ND35, SUBSTRS(2:N), N-1, N_MAX )

      !--------------------------
      ! ND36: Anthro emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:37' )
      READ( SUBSTRS(1), * ) ND36
      CALL SET_TINDEX( am_I_Root, 
     &                 36, ND36, SUBSTRS(2:N), N-1, NEMANTHRO )      

      !--------------------------
      ! ND37: Updraft scav frac
      !--------------------------
      N_TMP = N_TRACERS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:38' )
      READ( SUBSTRS(1), * ) ND37
      CALL SET_TINDEX( am_I_Root, 37, ND37, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND38: Cld conv losses
      !--------------------------
      N_TMP = N_TRACERS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:39' )
      READ( SUBSTRS(1), * ) ND38
      IF ( .not. LWETD .and. .not. LCONV ) ND38 = 0
      CALL SET_TINDEX( am_I_Root, 38, ND38, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND39: Wet scav losses
      !--------------------------
      N_TMP = N_TRACERS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:40' )
      READ( SUBSTRS(1), * ) ND39
      IF ( .not. LWETD ) ND39 = 0
      CALL SET_TINDEX( am_I_Root, 39, ND39, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND41: Afternoon PBL ht
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:41' )
      READ( SUBSTRS(1), * ) ND41
      CALL SET_TINDEX( am_I_Root, 41, ND41, SUBSTRS(2:N), N-1, PD41 )

      !--------------------------
      ! ND42: SOA concentrations
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:42' )
      READ( SUBSTRS(1), * ) ND42
      CALL SET_TINDEX( am_I_Root, 42, ND42, SUBSTRS(2:N), N-1, PD42 )

      !--------------------------
      ! ND43: OH, HO2
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:43' )
      READ( SUBSTRS(1), * ) ND43
      IF ( .not. ( ITS_A_FULLCHEM_SIM().or.ITS_A_CH4_SIM() )) ND43 = 0
      CALL SET_TINDEX( am_I_Root, 43, ND43, SUBSTRS(2:N), N-1, PD43 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:44' )
      READ( SUBSTRS(1:N), * ) HR1_OH, HR2_OH

      !--------------------------
      ! ND44 drydep vel & flux
      !--------------------------

      ! Number of tracers depends on simulation type
      IF ( ITS_A_TAGOX_SIM() ) THEN
         N_TMP = N_TRACERS 
      ELSE
         N_TMP = NUMDEP
      ENDIF

      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:46' )
      READ( SUBSTRS(1), * ) ND44
      IF ( .not. LDRYD ) ND44 = 0
      CALL SET_TINDEX( am_I_Root, 44, ND44, SUBSTRS(2:N), N-1, N_TMP )

      !--------------------------
      ! ND45: Tracer conc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:47' )
      READ( SUBSTRS(1), * ) ND45
      CALL SET_TINDEX( am_I_Root, 45, ND45, SUBSTRS(2:N), N-1, PD45 )

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_diagnostic_menu:48' ) 
      READ( SUBSTRS(1:N), * ) HR1_OTH, HR2_OTH

      !--------------------------
      ! ND46: Biogenic sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:49' )
      READ( SUBSTRS(1), * ) ND46
      CALL SET_TINDEX( am_I_Root, 46, ND46, SUBSTRS(2:N), N-1, PD46 )

      !--------------------------
      ! ND47: 24h avg tracer
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:50' )
      READ( SUBSTRS(1), * ) ND47
      CALL SET_TINDEX( am_I_Root, 
     &                 47, ND47, SUBSTRS(2:N), N-1, N_TRACERS )

      !--------------------------
      ! ND52: GAMMA HO2
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:51' )
      READ( SUBSTRS(1), * ) ND52
      CALL SET_TINDEX( am_I_Root, 52, ND52, SUBSTRS(2:N), N-1, PD52 )

      !--------------------------
      ! ND53: POPS
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:52' )
      READ( SUBSTRS(1), * ) ND53
      IF ( .not. ITS_A_POPS_SIM() ) ND53 = 0
      CALL SET_TINDEX( am_I_Root, 53, ND53, SUBSTRS(2:N), N-1, PD53 )

      !--------------------------
      ! ND54: Time in troposphere
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:53' )
      READ( SUBSTRS(1), * ) ND54
      CALL SET_TINDEX( am_I_Root, 54, ND54, SUBSTRS(2:N), N-1, 1 )

      !--------------------------
      ! ND55: Tropopause diags
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:54' )
      READ( SUBSTRS(1), * ) ND55
      CALL SET_TINDEX( am_I_Root, 55, ND55, SUBSTRS(2:N), N-1, PD55 )

      !--------------------------
      ! ND56: Lightning flashes
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:55' )
      READ( SUBSTRS(1), * ) ND56
      CALL SET_TINDEX( am_I_Root, 56, ND56, SUBSTRS(2:N), N-1, PD56 )

      !--------------------------
      ! ND57: Potential temperature
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:56' )
      READ( SUBSTRS(1), * ) ND57
      CALL SET_TINDEX( am_I_Root, 57, ND57, SUBSTRS(2:N), N-1, PD57 )

      !--------------------------
      ! ND58: CH4 Emissions 
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:57' )
      READ( SUBSTRS(1), * ) ND58
      CALL SET_TINDEX( am_I_Root, 58, ND58, SUBSTRS(2:N), N-1, PD58 )

      !--------------------------
      ! ND59: TOMAS aerosol emiss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:58' )
      READ( SUBSTRS(1), * ) ND59
      CALL SET_TINDEX( am_I_Root, 59, ND59, SUBSTRS(2:N), N-1, PD59 )

      !--------------------------
      ! ND60: Wetland Fraction
      ! ND60: TOMAS rate
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:59' )
      READ( SUBSTRS(1), * ) ND60
      CALL SET_TINDEX( am_I_Root, 60, ND60, SUBSTRS(2:N), N-1, PD60 )

      !--------------------------
      ! ND61: 3-D TOMAS rate
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:60' )
      READ( SUBSTRS(1), * ) ND61
      CALL SET_TINDEX( am_I_Root, 61, ND61, SUBSTRS(2:N), N-1, PD61 )

      !--------------------------
      ! ND62: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:61' )
      READ( SUBSTRS(1), * ) ND62
      CALL SET_TINDEX( am_I_Root, 62, ND62, SUBSTRS(2:N), N-1, PD62 )

      !--------------------------
      ! ND64: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:63' )
      READ( SUBSTRS(1), * ) ND64
      CALL SET_TINDEX( am_I_Root, 64, ND64, SUBSTRS(2:N), N-1, PD64 )

      !--------------------------
      ! ND66: DAO 3-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:64' )
      READ( SUBSTRS(1), * ) ND66
      CALL SET_TINDEX( am_I_Root, 66, ND66, SUBSTRS(2:N), N-1, PD66 )

      !--------------------------
      ! ND67: DAO 2-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:65' )
      READ( SUBSTRS(1), * ) ND67
      CALL SET_TINDEX( am_I_Root, 67, ND67, SUBSTRS(2:N), N-1, PD67 )

      !--------------------------
      ! ND68: Air masses etc
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:66' )
      READ( SUBSTRS(1), * ) ND68
      CALL SET_TINDEX( am_I_Root, 68, ND68, SUBSTRS(2:N), N-1, PD68 )

      !--------------------------
      ! ND69: Surface areas
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:67' )
      READ( SUBSTRS(1), * ) ND69
      CALL SET_TINDEX( am_I_Root, 69, ND69, SUBSTRS(2:N), N-1, PD69 )

      !--------------------------
      ! ND70: Debug info
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_diagnostic_menu:68' )
      READ( SUBSTRS(1), * ) ND70
      LPRT = ( ND70 > 0 )
      CALL SET_TINDEX( am_I_Root, 70, ND70, SUBSTRS(2:N), N-1, PD70 )
     
      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_diagnostic_menu:69' )

      !=================================================================
      ! Call other setup routines
      !=================================================================
      Input_Opt%ND01               = ND01
      Input_Opt%ND02               = ND02
      Input_Opt%ND03               = ND03
      Input_Opt%ND04               = ND04
      Input_Opt%ND05               = ND05
      Input_Opt%ND06               = ND06
      Input_Opt%ND07               = ND07
      Input_Opt%ND08               = ND08
      Input_Opt%ND09               = ND09
      Input_Opt%ND10               = ND10
      Input_Opt%ND11               = ND11
      Input_Opt%ND12               = ND12
      Input_Opt%ND13               = ND13
      Input_Opt%ND14               = ND14
      Input_Opt%ND15               = ND15
      Input_Opt%ND16               = ND16
      Input_Opt%ND17               = ND17
      Input_Opt%ND18               = ND18
      Input_Opt%ND19               = ND19
      Input_Opt%ND20               = ND20
      Input_Opt%ND21               = ND21
      Input_Opt%ND22               = ND22
      Input_Opt%ND23               = ND23
      Input_Opt%ND24               = ND24
      Input_Opt%ND25               = ND25
      Input_Opt%ND26               = ND26
      Input_Opt%ND27               = ND27
      Input_Opt%ND28               = ND28
      Input_Opt%ND29               = ND29
      Input_Opt%ND30               = ND30
      Input_Opt%ND31               = ND31
      Input_Opt%ND32               = ND32
      Input_Opt%ND33               = ND33
      Input_Opt%ND34               = ND34
      Input_Opt%ND35               = ND35
      Input_Opt%ND36               = ND36
      Input_Opt%ND37               = ND37
      Input_Opt%ND38               = ND38
      Input_Opt%ND39               = ND39
      Input_Opt%ND40               = ND40
      Input_Opt%ND41               = ND41
      Input_Opt%ND42               = ND42
      Input_Opt%ND43               = ND43
      Input_Opt%ND44               = ND44
      Input_Opt%ND45               = ND45
      Input_Opt%ND46               = ND46
      Input_Opt%ND47               = ND47
      Input_Opt%ND48               = ND48
      Input_Opt%ND49               = ND49
      Input_Opt%ND50               = ND50
      Input_Opt%ND51               = ND51
      Input_Opt%ND52               = ND52
      Input_Opt%ND53               = ND53
      Input_Opt%ND54               = ND54
      Input_Opt%ND55               = ND55
      Input_Opt%ND56               = ND56
      Input_Opt%ND57               = ND57
      Input_Opt%ND58               = ND58
      Input_Opt%ND59               = ND59
      Input_Opt%ND60               = ND60
      Input_Opt%ND61               = ND61
      Input_Opt%ND62               = ND62
      Input_Opt%ND63               = ND63
      Input_Opt%ND64               = ND64
      Input_Opt%ND66               = ND66
      Input_Opt%ND67               = ND67
      Input_Opt%ND68               = ND68
      Input_Opt%ND69               = ND69
      Input_Opt%ND70               = ND70
      Input_Opt%LPRT               = ( ND70 > 0 )

      ! Loop over # of diagnostics
      DO M = 1, Input_Opt%MAX_DIAG
         Input_Opt%TCOUNT(M)       = TCOUNT(M)
         Input_Opt%TMAX(M)         = TMAX(M)

         ! Loop over tracers per diagnostic
         DO N = 1, Input_Opt%MAX_TRCS
            Input_Opt%TINDEX(M,N)  = TINDEX(M,N)
         ENDDO
      ENDDO

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Call other setup routines
      !================================================================

      ! Allocate diagnostic arrays
      CALL INIT_DIAG03
      CALL INIT_DIAG04
      CALL INIT_DIAG41
      CALL INIT_DIAG42
      CALL INIT_DIAG53
      CALL INIT_DIAG56

      ! Enable Mean OH (or CH3CCl3) diag for runs which need it
      CALL INIT_DIAG_OH( am_I_Root )

      ! Expand YYYYMMDD tokens in the bpch file name
      CALL EXPAND_DATE( BPCH_FILE, GET_NYMDb(), GET_NHMSb() )

      ! Find a free file LUN
      IU_BPCH = findFreeLUN()

      ! Open the binary punch file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, BPCH_FILE )

      END SUBROUTINE READ_DIAGNOSTIC_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_tindex
!
! !DESCRIPTION: Subroutine SET\_TINDEX sets the TINDEX and TMAX arrays, 
!  which determine how many tracers to print to the punch file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TINDEX( am_I_Root, 
     &                       N_DIAG, L_DIAG, SUBSTRS, N, NMAX )
!
! !USES:
!
#if   defined( TOMAS )
      USE CHARPAK_MOD, ONLY : TXTEXT   ! (win, 7/14/09)
#endif

      USE CMN_SIZE_MOD                              ! Size parameters
      USE CMN_DIAG_MOD                              ! TMAX, TINDEX
!
! !INPUT PARAMETERS: 
!
      INTEGER,            INTENT(IN) :: N_DIAG      ! GEOS-Chem diagnostic #
      INTEGER,            INTENT(IN) :: N           ! # of valid substrs passed
      INTEGER,            INTENT(IN) :: NMAX        ! Max # of tracers allowed
      INTEGER,            INTENT(IN) :: L_DIAG      ! # of levels to save
      CHARACTER(LEN=255), INTENT(IN) :: SUBSTRS(N)  ! Substrs passed from
                                                    !  READ_DIAGNOSTIC_MENU
      LOGICAL,            INTENT(IN) :: am_I_Root   ! Is this the root CPU?

! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Bug fix: now do not drop the last tracer number if "all" is not
!        explicitly specified (tmf, bmy, 11/15/04)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: IS_ALL 
      INTEGER             :: M
#if   defined( TOMAS )
      INTEGER             :: NN,     COL,     IFLAG, TC     ! (win, 7/14/09)
      CHARACTER (LEN=255) :: WORD,   SUBWORD, TMP1,  TMP2   ! (win, 7/14/09)
      INTEGER             :: MINTMP, MAXTMP                 ! (win, 7/14/09)
#endif

      !=================================================================
      ! SET_TINDEX begins here!
      !=================================================================     

      ! Error check
      IF ( N < 1 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) 'ERROR: N must be 1 or greater!'
            WRITE( 6, '(a)' ) 'STOP in SET_TINDEX (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         ENDIF
         STOP
      ENDIF

      !=================================================================
      ! If the word "all" is present, then set TMAX, TINDEX to all
      ! available tracers for the given diagnostic.  Otherwise, just
      ! use the tracers that were read in from the line
      !=================================================================
      IF ( TRIM( SUBSTRS(1) ) == 'all'  .or. 
     &     TRIM( SUBSTRS(1) ) == 'ALL' ) THEN 

         ! TMAX is the max # of tracers to print out
         TMAX(N_DIAG) = NMAX 

         ! Fill TINDEX with all possible diagnostic tracer numbers
         DO M = 1, TMAX(N_DIAG)
            TINDEX(N_DIAG,M) = M
         ENDDO

         ! Set flag
         IS_ALL = .TRUE. 

      ELSE 

#if   defined( TOMAS )
!(win, 7/14/09)  use TXTEXT and split the read in characters by -

         COL   = 1
         NN    = 0
         SUBWORD  = ''
         IFLAG = 0

         ! Use explicit DO-loop
         DO M = 1, N
            WORD = SUBSTRS(M)

            ! Check if the characters are a range with - in the middle
            CALL TXTEXT ( '-', WORD, COL, SUBWORD, IFLAG )

            ! Found a dash!  Get the numbers on both sides of the dash
            ! since these the min and max of the tracer range
            IF ( IFLAG == 0 ) THEN
               TMP1 = TRIM( WORD(      1:COL-1      ) )
               TMP2 = TRIM( WORD( COL+1:LEN_TRIM( WORD ) ) )

               READ( TMP1, * ) MINTMP
               READ( TMP2, * ) MAXTMP
               
               DO TC = MINTMP, MAXTMP
                  NN = NN + 1
                  TINDEX( N_DIAG, NN ) = TC
               ENDDO

            ! If we haven't found a dash, then there is only one number,
            ! so that number is both the min and max of the tracer range
            ELSE IF ( IFLAG == -1 ) THEN
               NN = NN + 1
               TMP1 = TRIM( WORD )
               READ( TMP1, * ) TINDEX( N_DIAG, NN )
            ENDIF
               
         ENDDO

         ! Set TMAX to the counted # of tracers
         TMAX( N_DIAG ) = NN
#else
         ! Otherwise, set TMAX, TINDEX to the # of tracers
         ! listed in "input.ctm" -- need some error checks too
         TMAX(N_DIAG) = N

         ! Use explicit DO-loop
         DO M = 1, N
            READ( SUBSTRS(M:M), * ) TINDEX(N_DIAG,M)
         ENDDO

#endif
         ! Set flag
         IS_ALL = .FALSE.

      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================

      ! First-time printing only
      IF ( FIRST ) THEN
         IF( am_I_Root ) THEN 
            WRITE( 6,'(/,a)' ) 'DIAGNOSTIC MENU'
            WRITE( 6,'(  a)' ) '---------------'
            WRITE( 6,'(  a)' ) 'Diag    L   Tracers being saved to disk'
         ENDIF
         FIRST = .FALSE.
      ENDIF

      ! Test if all tracers are being printed out
      IF ( IS_ALL ) THEN

         ! Print abbreviated output string
         IF ( L_DIAG > 0 ) THEN
           IF ( am_I_Root ) THEN
              WRITE( 6, 100 ) N_DIAG, L_DIAG, 1, TMAX(N_DIAG)
           ENDIF
 100       FORMAT( 'ND', i2.2, 2x, i3, 1x, i3, ' -', i3 ) 
         ENDIF

      ELSE

         ! Or just list each tracer
         ! Print each diagnostic and # of tracers that will print out
         IF ( L_DIAG > 0 ) THEN 
            IF ( am_I_Root ) THEN
               WRITE( 6, 110 ) N_DIAG, L_DIAG, 
     &                         ( TINDEX(N_DIAG,M), M=1,TMAX(N_DIAG) )
            ENDIF
 110        FORMAT( 'ND', i2, 2x, i3, 1x, 100i3 ) 
         ENDIF

      ENDIF

      END SUBROUTINE SET_TINDEX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_planeflight_menu
!
! !DESCRIPTION: Subroutine READ\_PLANEFLIGHT\_MENU reads the PLANEFLIGHT MENU 
!  section of the GEOS-Chem input file.  This turns on the ND40 flight track 
!  diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PLANEFLIGHT_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE PLANEFLIGHT_MOD, ONLY : SET_PLANEFLIGHT

      USE CMN_SIZE_MOD                   ! MAXFAM
      USE CMN_DIAG_MOD                   ! ND40
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: DO_PF
      INTEGER            :: N
      CHARACTER(LEN=255) :: IFILE
      CHARACTER(LEN=255) :: OFILE
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_PLANEFLIGHT_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_PLANEFLIGHT_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      ND40  = 0
      DO_PF = .FALSE.
      IFILE = ''
      OFILE = ''

      ! Turn on planeflight diagnostic?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_planeflight_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_PF

      ! Set ND40 flag from DO_PF
      IF ( DO_PF ) ND40 = 1

      ! Input file name (w/ flight track data points)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_planeflight_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) IFILE

      ! Output file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_planeflight_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) OFILE

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_PF    = DO_PF
      Input_Opt%PF_IFILE = IFILE
      Input_Opt%PF_OFILE = OFILE

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'PLANEFLIGHT MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100 ) 'Turn on planeflight diag?   : ', DO_PF
         WRITE( 6, 110 ) 'Flight track input file     : ', TRIM(IFILE)
         WRITE( 6, 110 ) 'Output file name            : ', TRIM(OFILE)
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass variables to "planeflight_mod.f"
      CALL SET_PLANEFLIGHT( DO_PF, IFILE, OFILE )
   
      END SUBROUTINE READ_PLANEFLIGHT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd48_menu
!
! !DESCRIPTION: Subroutine READ\_ND48\_MENU reads the ND48 MENU section of the 
!  GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND48_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG48_MOD, ONLY : INIT_DIAG48, ND48_MAX_STATIONS
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Bug fix: ND48 stations should now be read correctly. (bmy, 3/6/06)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: DO_ND48
      INTEGER            :: N, S
      INTEGER            :: FREQ
      INTEGER            :: N_STA
      INTEGER            :: IARR(ND48_MAX_STATIONS)
      INTEGER            :: JARR(ND48_MAX_STATIONS)
      INTEGER            :: LARR(ND48_MAX_STATIONS)
      INTEGER            :: NARR(ND48_MAX_STATIONS)
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255) :: FILE
      CHARACTER(LEN=10)  :: C

      !=================================================================
      ! READ_ND48_MENU begins here!
      !=================================================================
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND48_MENU ("input_mod.f")' )
      ENDIF

      ! Turn on ND48 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND48

      ! Timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Frequency
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:3' )
      READ( SUBSTRS(1:N), * ) FREQ

      ! Number of stations 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:4' )
      READ( SUBSTRS(1:N), * ) N_STA

      ! Initialize
      IARR(:) = 0
      JARR(:) = 0
      LARR(:) = 0
      NARR(:) = 0
      
      ! Read individual stations
      DO S = 1, N_STA
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'read_nd48_menu:5' )
         READ( SUBSTRS(1), * ) IARR(S) 
         READ( SUBSTRS(2), * ) JARR(S)
         READ( SUBSTRS(3), * ) LARR(S) 
         READ( SUBSTRS(4), * ) NARR(S) 
      ENDDO

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND48            = DO_ND48
      Input_Opt%ND48_FILE          = FILE
      Input_Opt%ND48_FREQ          = FREQ
      Input_Opt%ND48_N_STA         = N_STA
      Input_Opt%ND48_IARR(1:N_STA) = IARR(1:N_STA)
      Input_Opt%ND48_JARR(1:N_STA) = JARR(1:N_STA)
      Input_Opt%ND48_LARR(1:N_STA) = LARR(1:N_STA)
      Input_Opt%ND48_NARR(1:N_STA) = NARR(1:N_STA)

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND48 STATION TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '----------------------------'
         WRITE( 6, 100 ) 'Turn on ND48 timeseries?    : ', DO_ND48
         WRITE( 6, 110 ) 'ND48 timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 120 ) 'ND48 save frequency [min]   : ', FREQ
      ENDIF

      DO S = 1, N_STA
         WRITE( 6, 130 ) S, IARR(S), JARR(S), LARR(S), NARR(S) 
      ENDDO

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, I5    )                   
 130  FORMAT( 'ND48 timeseries station', i4, ' : ', 4i5 )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize for ND48 timeseries
      CALL INIT_DIAG48( DO_ND48, FREQ, N_STA, IARR, 
     &                  JARR,    LARR, NARR,  FILE )

      END SUBROUTINE READ_ND48_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd49_menu
!
! !DESCRIPTION: Subroutine READ\_ND49\_MENU reads the ND49 MENU section of 
!  the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND49_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG49_MOD, ONLY : INIT_DIAG49
      USE ERROR_MOD,  ONLY : ERROR_STOP

      USE CMN_SIZE_MOD     ! Size parameters
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: DO_ND49
      INTEGER             :: N,    I,         AS
      ! Increased to 121 from 100 (mpb,2009)
      INTEGER             :: ND49, N_TRACERS, TRACERS(121)
      INTEGER             :: IMIN, IMAX,      FREQ
      INTEGER             :: JMIN, JMAX,      N_ND49
      INTEGER             :: LMIN, LMAX
      CHARACTER(LEN=255)  :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255)  :: FILE

      !=================================================================
      ! READ_ND49_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND49_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      ND49       = 0
      TRACERS(:) = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND49

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND49, -1, 'read_nd49_menu:3' )
      DO N = 1, N_ND49
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! FREQ
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:4' )
      READ( SUBSTRS(1:N), * ) FREQ

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:5' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:6' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:7' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:8' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND49                = DO_ND49
      Input_Opt%ND49_FILE              = FILE
      Input_Opt%ND49_TRACERS(1:N_ND49) = TRACERS(1:N_ND49)
      Input_Opt%ND49_FREQ              = FREQ
      Input_Opt%ND49_IMIN              = IMIN
      Input_Opt%ND49_IMAX              = IMAX
      Input_Opt%ND49_JMIN              = JMIN
      Input_Opt%ND49_JMAX              = JMAX
      Input_Opt%ND49_LMIN              = LMIN
      Input_Opt%ND49_LMAX              = LMAX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND49 3-D INSTANTANEOUS TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '--------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND49 timeseries?    : ', DO_ND49
         WRITE( 6, 110 ) 'ND49 timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 120 ) 'ND49 timeseries tracers     : ', 
     &                    ( TRACERS(N), N=1, N_ND49 )
         WRITE( 6, 130 ) 'ND49 save frequency [min]   : ', FREQ
         WRITE( 6, 130 ) 'ND49 longitude limits       : ', IMIN, IMAX
         WRITE( 6, 130 ) 'ND49 latitude  limits       : ', JMIN, JMAX
         WRITE( 6, 130 ) 'ND49 level     limits       : ', LMIN, LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize for ND49 timeseries
      CALL INIT_DIAG49( DO_ND49, N_ND49, TRACERS, IMIN, 
     &                  IMAX,    JMIN,   JMAX,    LMIN,    
     &                  LMAX,    FREQ,   FILE )

      END SUBROUTINE READ_ND49_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd50_menu
!
! !DESCRIPTION: Subroutine READ\_ND50\_MENU reads the ND50 MENU section of 
!  the GEOS-Chem input file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND50_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG50_MOD,  ONLY : INIT_DIAG50
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE LOGICAL_MOD, ONLY : LND50_HDF
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now include option to save ND51 diagnostic to HDF5 file format
!        (amv, bmy, 12/21/09)
!  (2 ) Increase tracer number to 121. (ccc, 4/20/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: DO_ND50
      INTEGER              :: N,      I,       AS  
      ! Increased # of tracers to 121. (ccc, 4/20/10)
      INTEGER              :: N_ND50, IMIN,    IMAX, TRACERS(121)   
      INTEGER              :: JMIN,   JMAX,    LMIN, LMAX
      CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255)   :: FILE

      !=================================================================
      ! READ_ND50_MENU begins here!
      !=================================================================
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND50_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      TRACERS(:) = 0

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND50

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:3' )
      READ( SUBSTRS(1:N), * ) LND50_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      LND50_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND50, -1, 'read_nd50_menu:4' )
      DO N = 1, N_ND50
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:5' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:6' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:7' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:8' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND50                = DO_ND50
      Input_Opt%ND50_FILE              = FILE
      Input_Opt%LND50_HDF              = LND50_HDF
      Input_Opt%ND50_TRACERS(1:N_ND50) = TRACERS(1:N_ND50)
      Input_Opt%ND50_IMIN              = IMIN
      Input_Opt%ND50_IMAX              = IMAX
      Input_Opt%ND50_JMIN              = JMIN
      Input_Opt%ND50_JMAX              = JMAX
      Input_Opt%ND50_LMIN              = LMIN
      Input_Opt%ND50_LMAX              = LMAX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND50 3-D 24hr AVG TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '---------------------------------'
         WRITE( 6, 100 ) 'Turn on ND50 timeseries?    : ', DO_ND50
         WRITE( 6, 110 ) 'ND50 timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 100 ) 'Output as HDF?              : ', LND50_HDF
         WRITE( 6, 120 ) 'ND50 timeseries tracers     : ', 
     &                    ( TRACERS(N), N=1, N_ND50 )
         WRITE( 6, 130 ) 'ND50 longitude limits       : ', IMIN, IMAX
         WRITE( 6, 130 ) 'ND50 latitude  limits       : ', JMIN, JMAX
         WRITE( 6, 130 ) 'ND50 level     limits       : ', LMIN, LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize parameters for ND49 timeseries
      CALL INIT_DIAG50( DO_ND50, N_ND50, TRACERS, IMIN, IMAX, 
     &                  JMIN,    JMAX,   LMIN,    LMAX, FILE )

      END SUBROUTINE READ_ND50_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd51_menu
!
! !DESCRIPTION: Subroutine READ\_ND51\_MENU reads the ND51 MENU section of 
!  the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND51_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG51_MOD,  ONLY : INIT_DIAG51
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE LOGICAL_MOD, ONLY : LND51_HDF

      USE CMN_SIZE_MOD      ! Size parameters
      USE CMN_DIAG_MOD      ! NDxx flags
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now include option to save ND51 diagnostic to HDF5 file format
!        (amv, bmy, 12/21/09)
!  (2 ) Increase # of tracers to 121 (ccc, 4/20/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: DO_ND51
      INTEGER              :: N,      I,       AS
      ! Increased to 121 from 100 (mpb,2009)
      INTEGER              :: N_ND51, FREQ,    TRACERS(121)
      INTEGER              :: IMIN,   IMAX,    JMIN
      INTEGER              :: JMAX,   LMIN,    LMAX
      REAL*8               :: HR1,    HR2,     HR_WRITE
      CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255)   :: FILE

      !=================================================================
      ! READ_ND51_MENU begins here!
      !=================================================================
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND51_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      TRACERS(:) = 0

      ! Turn on ND51 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND51

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:3' )
      READ( SUBSTRS(1:N), * ) LND51_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      LND51_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND51, -1, 'read_nd51_menu:4' )
      DO N = 1, N_ND51
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! NHMS_WRITE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:6' )
      READ( SUBSTRS(1:N), * ) HR_WRITE

      ! HR1, HR2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7' )
      READ( SUBSTRS(1:N), * ) HR1, HR2

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:8' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:9' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:10' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:11' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND51                = DO_ND51
      Input_Opt%ND51_FILE              = FILE
      Input_Opt%LND51_HDF              = LND51_HDF
      Input_Opt%ND51_TRACERS(1:N_ND51) = TRACERS(1:N_ND51)
      Input_Opt%ND51_HR_WRITE          = HR_WRITE
      Input_Opt%ND51_HR1               = HR1 
      Input_Opt%ND51_HR2               = HR2
      Input_Opt%ND51_IMIN              = IMIN
      Input_Opt%ND51_IMAX              = IMAX
      Input_Opt%ND51_JMIN              = JMIN
      Input_Opt%ND51_JMAX              = JMAX
      Input_Opt%ND51_LMIN              = LMIN
      Input_Opt%ND51_LMAX              = LMAX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND51 MORNING OR AFTERNOON TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '-----------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND51 timeseries?    : ', DO_ND51
         WRITE( 6, 110 ) 'ND51 timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 100 ) 'Output as HDF?              : ', LND51_HDF
         WRITE( 6, 120 ) 'ND41 timeseries tracers     : ', 
     &                    ( TRACERS(N), N=1, N_ND51 )
         WRITE( 6, 140 ) 'ND51 hour to write to disk  : ', HR_WRITE
         WRITE( 6, 140 ) 'ND51 averaging period [GMT] : ', HR1,  HR2
         WRITE( 6, 130 ) 'ND51 longitude limits       : ', IMIN, IMAX
         WRITE( 6, 130 ) 'ND51 latitude  limits       : ', JMIN, JMAX
         WRITE( 6, 130 ) 'ND51 altitude  limits       : ', LMIN, LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )
 140  FORMAT( A, 2F5.1 )

      !=================================================================
      ! Call setup routine from other F90 modules
      !=================================================================

      ! Initialize parameters for ND51 timeseries
      CALL INIT_DIAG51( DO_ND51, N_ND51, TRACERS, HR_WRITE, 
     &                  HR1,     HR2,    IMIN,    IMAX,   
     &                  JMIN,    JMAX,   LMIN,    LMAX,  FILE )

      END SUBROUTINE READ_ND51_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd51b_menu
!
! !DESCRIPTION: Subroutine READ\_ND51b\_MENU reads the ND51 MENU section 
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND51b_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG51b_MOD, ONLY : INIT_DIAG51b
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE LOGICAL_MOD, ONLY : LND51b_HDF

      USE CMN_SIZE_MOD      ! Size parameters
      USE CMN_DIAG_MOD      ! NDxx flags
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  21 Dec 2009 - Aaron van D - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: DO_ND51
      INTEGER              :: N,      I,       AS
      ! Increase to 121 from 100 (ccc, 4/20/10)
      INTEGER              :: N_ND51, FREQ,    TRACERS(121)
      INTEGER              :: IMIN,   IMAX,    JMIN
      INTEGER              :: JMAX,   LMIN,    LMAX
      REAL*8               :: HR1,    HR2,     HR_WRITE
      CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255)   :: FILE

      !=================================================================
      ! READ_ND51_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND51b_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      TRACERS(:) = 0

      ! Turn on ND51 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND51

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:3' )
      READ( SUBSTRS(1:N), * ) LND51b_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      LND51b_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND51, -1, 'read_nd51b_menu:4' )
      DO N = 1, N_ND51
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! NHMS_WRITE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:5' )
      READ( SUBSTRS(1:N), * ) HR_WRITE

      ! HR1, HR2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:6' )
      READ( SUBSTRS(1:N), * ) HR1, HR2

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:7' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:8' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:9' )
      READ( SUBSTRS(1:N), * ) LMIN, LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:10' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND51b                 = DO_ND51
      Input_Opt%ND51b_FILE               = FILE
      Input_Opt%LND51b_HDF               = LND51b_HDF
      Input_Opt%ND51b_TRACERS(1:N_ND51)  = TRACERS(1:N_ND51)
      Input_Opt%ND51b_HR_WRITE           = HR_WRITE
      Input_Opt%ND51b_HR1                = HR1 
      Input_Opt%ND51b_HR2                = HR2
      Input_Opt%ND51b_IMIN               = IMIN
      Input_Opt%ND51b_IMAX               = IMAX
      Input_Opt%ND51b_JMIN               = JMIN
      Input_Opt%ND51b_JMAX               = JMAX
      Input_Opt%ND51b_LMIN               = LMIN
      Input_Opt%ND51b_LMAX               = LMAX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' )'ND51b MORNING OR AFTERNOON TIMESERIES MENU'
         WRITE( 6, '(  a)' )'-----------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND51b timeseries?    : ', DO_ND51
         WRITE( 6, 110 ) 'ND51b timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 100 ) 'Output as HDF?               : ', LND51b_HDF
         WRITE( 6, 120 ) 'ND41 timeseries tracers      : ',
     &                    ( TRACERS(N), N=1, N_ND51 )
         WRITE( 6, 140 ) 'ND51b hour to write to disk  : ', HR_WRITE
         WRITE( 6, 140 ) 'ND51b averaging period [GMT] : ', HR1,  HR2
         WRITE( 6, 130 ) 'ND51b longitude limits       : ', IMIN, IMAX
         WRITE( 6, 130 ) 'ND51b latitude  limits       : ', JMIN, JMAX
         WRITE( 6, 130 ) 'ND51b altitude  limits       : ', LMIN, LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )
 140  FORMAT( A, 2F5.1 )

      !=================================================================
      ! Call setup routine from other F90 modules
      !=================================================================

      ! Initialize parameters for ND51b timeseries
      CALL INIT_DIAG51b( DO_ND51, N_ND51, TRACERS, HR_WRITE,
     &                   HR1,     HR2,    IMIN,    IMAX,
     &                   JMIN,    JMAX,   LMIN,    LMAX,  FILE )

      END SUBROUTINE READ_ND51b_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nd63_menu
!
! !DESCRIPTION: Subroutine READ\_ND63\_MENU reads the ND63 MENU section
!  of the GEOS-Chem input file. (gvinken, 02/25/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ND63_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIAG63_MOD, ONLY : INIT_DIAG63
      USE ERROR_MOD,  ONLY : ERROR_STOP

      USE CMN_SIZE_MOD     ! Size parameters
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version
!  07 Feb 2012 - M. Payer    - Added ProTeX headers
!  24 Feb 2012 - M. Payer    - Renamed routine from READ_ND59_MENU to
!                              READ_ND63 MENU. ND59 is used by TOMAS.
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: DO_ND63
      INTEGER             :: N,    I,         AS
      INTEGER             :: ND63, N_TRACERS, TRACERS(121)
      INTEGER             :: IMIN, IMAX,      FREQ
      INTEGER             :: JMIN, JMAX,      N_ND63
      CHARACTER(LEN=255)  :: SUBSTRS(MAXDIM), MSG
      CHARACTER(LEN=255)  :: FILE

      !=================================================================
      ! READ_ND63_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_ND63_MENU ("input_mod.f")' )
      ENDIF

      ! Initialize
      ND63       = 0
      TRACERS(:) = 0

      ! Turn on ND63 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_ND63

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, N_ND63, -1, 'read_nd63_menu:3' )
      DO N = 1, N_ND63
         READ( SUBSTRS(N), * ) TRACERS(N)
      ENDDO

      ! FREQ
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:4' )
      READ( SUBSTRS(1:N), * ) FREQ

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd63_menu:5' )
      READ( SUBSTRS(1:N), * ) IMIN, IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd63_menu:6' )
      READ( SUBSTRS(1:N), * ) JMIN, JMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:8' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_ND63              = DO_ND63
      Input_Opt%ND63_FILE            = FILE
      Input_Opt%ND63_TRACERS(1:ND63) = TRACERS(1:N_ND63)
      Input_Opt%ND63_FREQ            = FREQ
      Input_Opt%ND63_IMIN            = IMIN
      Input_Opt%ND63_IMAX            = IMAX
      Input_Opt%ND63_JMIN            = JMIN
      Input_Opt%ND63_JMAX            = JMAX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND63 SHIP TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '--------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND63 timeseries?    : ', DO_ND63
         WRITE( 6, 110 ) 'ND63 timeseries file name   : ', TRIM( FILE )
         WRITE( 6, 120 ) 'ND63 timeseries tracers     : ', 
     &                    ( TRACERS(N), N=1, N_ND63 )
         WRITE( 6, 130 ) 'ND63 save frequency [min]   : ', FREQ
         WRITE( 6, 130 ) 'ND63 longitude limits       : ', IMIN, IMAX
         WRITE( 6, 130 ) 'ND63 latitude  limits       : ', JMIN, JMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Initialize for ND63 timeseries
      CALL INIT_DIAG63( DO_ND63, N_ND63, TRACERS, IMIN, 
     &                  IMAX,    JMIN,   JMAX,  FREQ,   FILE )

      END SUBROUTINE READ_ND63_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_prod_loss_menu 
!
! !DESCRIPTION: Subroutine READ\_PROD\_LOSS\_MENU reads the PROD AND LOSS MENU 
!  section of the GEOS-Chem input file 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PROD_LOSS_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CHARPAK_MOD, ONLY : ISDIGIT,         STRSPLIT
      USE DIAG_PL_MOD, ONLY : INIT_DIAG_PL
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE TRACER_MOD,  ONLY : N_TRACERS,       ITS_A_TAGCO_SIM
      USE TRACER_MOD,  ONLY : ITS_A_TAGOX_SIM, ITS_AN_AEROSOL_SIM
      USE LOGICAL_MOD, ONLY : LKPP

      USE CMN_SIZE_MOD      ! MAXFAM
      USE CMN_DIAG_MOD      ! ND65
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Bug fixes.  Only error check # of prod/loss families for TagOx and 
!        TagCO runs if DO_SAVE_PL=T.  Also turn off this diagnostic for
!        the offline aerosol run. (bmy, 10/29/04)
!  (2 ) Add error trap is P/L families are asked when using KPP. (ccc, 3/10/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  08 Nov 2012 - R. Yantosca - Now save fields to the Input_Opt object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: EOF, DO_SAVE_PL, DO_SAVE_O3
      !FP_ISOP moved to CMN_SIZE (6/2009)
      !INTEGER, PARAMETER   :: MAXMEM=10
      INTEGER              :: F, M, N, NFAM
      INTEGER              :: FAM_NMEM(MAXFAM)
      REAL*8               :: FAM_COEF(MAXMEM,MAXFAM)
      CHARACTER(LEN=14 )   :: FAM_NAME(MAXFAM)
      CHARACTER(LEN=14 )   :: FAM_TYPE(MAXFAM)
      CHARACTER(LEN=14 )   :: FAM_MEMB(MAXMEM,MAXFAM)
      CHARACTER(LEN=255)   :: LOCATION,        NAME         
      CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM), MSG 

      !=================================================================
      ! READ_PROD_LOSS_MENU begins here!
      !=================================================================

      ! Location string
      LOCATION = 'READ_PROD_LOSS_MENU ("input_mod.f")'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      ! Initialize
      FAM_NAME(:)   = ''
      FAM_TYPE(:)   = ''
      FAM_NMEM(:)   = 0
      FAM_MEMB(:,:) = ''
      FAM_COEF(:,:) = 0d0

      !=================================================================
      ! Read info about prod & loss families
      !=================================================================

      ! Turn on production & loss diagnostic (e.g. ND65 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:1' )
      READ( SUBSTRS(1:N), * ) DO_SAVE_PL

      ! Read number of levels for ND65 diagnostic 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:2' )
      READ( SUBSTRS(1:N), * ) ND65

      ! Save P(O3) & L(O3) for tagged Ox run? (i.e. ND20 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:3' )
      READ( SUBSTRS(1:N), * ) DO_SAVE_O3

      ! Read number of families 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:4' )
      READ( SUBSTRS(1:N), * ) NFAM

      ! Loop over families
      DO F = 1, NFAM

         ! Get family members
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,  'read_prod_loss_menu:5' )

         ! The first entry is the family name (strip colon at end)
         NAME        = SUBSTRS(1)
         FAM_NAME(F) = NAME( 1:LEN_TRIM( NAME )-1 )

         ! Get family type as prod or loss
         IF ( FAM_NAME(F)(1:1) == 'P'   .or. 
     &        FAM_NAME(F)(1:1) == 'p' ) THEN
            FAM_TYPE(F) = 'prod'
         ELSE
            FAM_TYPE(F) = 'loss'
         ENDIF

         ! Number of member species in this prodloss family
         FAM_NMEM(F) = N - 1

         ! Loop over substrings
         DO M = 1, N-1

            ! Family member name
            NAME          =  TRIM( SUBSTRS(M+1) )

            ! Family member coefficient (set to 1 for now)
            FAM_COEF(M,F) = 1d0

            ! If first char is a digit ...
            IF ( ISDIGIT( NAME(1:1) ) ) THEN
               
               ! Save new family coefficient 
               READ( NAME(1:1), * ) FAM_COEF(M,F)

               ! Get the rest of the member name (skip digit)
               NAME = NAME( 2:LEN_TRIM(NAME) )
            ENDIF

            ! Family member name
            FAM_MEMB(M,F) = NAME
         ENDDO
      ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_prod_loss_menu:6' )

      !=================================================================
      ! Error check: ND65 not available for KPP
      !=================================================================

      IF ( DO_SAVE_PL .AND. LKPP ) THEN
         MSG = 'P/L families are incompatible with KPP. ' //
     &         'Please change one option'
         CALL ERROR_STOP( MSG, LOCATION )
      ENDIF

      !=================================================================
      ! Error check families for certain types of simulations
      !=================================================================

      ! Tagged Ox
      IF ( DO_SAVE_PL .and. ITS_A_TAGOX_SIM() ) THEN
         IF ( NFAM /= 2*N_TRACERS ) THEN
            MSG = 'Wrong number of P/L families for Tagged Ox!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      ! Tagged CO
      IF ( DO_SAVE_PL .and. ITS_A_TAGCO_SIM() ) THEN
         IF ( NFAM /= N_TRACERS+5 ) THEN
            MSG = 'Wrong number of P/L families for Tagged CO!'
            CALL ERROR_STOP( MSG, LOCATION )
         ENDIF
      ENDIF

      ! Offline aerosol -- turn off DO_SAVE_PL, since we use ND05,
      ! ND06, ND07, ND08, ND13 etc diagnostics instead of ND65
      IF ( ITS_AN_AEROSOL_SIM() ) THEN 
         DO_SAVE_PL    = .FALSE.
         DO_SAVE_O3    = .FALSE.
         ND65          = 0
         NFAM          = 0
         FAM_NAME(:)   = ''
         FAM_TYPE(:)   = ''
         FAM_NMEM(:)   = 0
         FAM_MEMB(:,:) = ''
         FAM_COEF(:,:) = 0d0
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC MENU'
         WRITE( 6, '(  a)' ) '---------------------------'      
         WRITE( 6, 100 ) 'Turn on prod & loss diag?   : ', DO_SAVE_PL
         WRITE( 6, 110 ) '# of levels for P/L diag    : ', ND65
         WRITE( 6, 100 ) 'Save P(Ox), L(Ox) for TagOx?: ', DO_SAVE_O3
      
         ! Loop over families
         DO F = 1, NFAM

            ! Write family name, type and # of members
            WRITE( 6, 120 ) FAM_NAME(F), FAM_TYPE(F)
         
            ! Write info about each constituent member
            DO M = 1, FAM_NMEM(F)
               WRITE( 6, 130 ) M, FAM_COEF(M,F), FAM_MEMB(M,F)
            ENDDO
         ENDDO
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
 120  FORMAT( /, 'Family=', A5, '  Type=', a4 )
 130  FORMAT( I3, 1X, F4.1, 1X, A5 )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%DO_SAVE_PL          = DO_SAVE_PL
      Input_Opt%LFAMILY             = ( DO_SAVE_PL .and. NFAM > 0 )
      Input_Opt%ND65                = ND65
      Input_Opt%DO_SAVE_O3          = DO_SAVE_O3
      Input_Opt%NFAM                = NFAM

      ! Loop over # families
      DO F = 1, NFAM
         Input_Opt%FAM_NAME(F)      = FAM_NAME(F) 
         Input_Opt%FAM_NMEM(F)      = FAM_NMEM(F)
         Input_Opt%FAM_TYPE(F)      = FAM_TYPE(F)

         ! Loop over # members/family
         DO M = 1, FAM_NMEM(F)   
            Input_Opt%FAM_COEF(M,F) = FAM_COEF(M,F)
            Input_Opt%FAM_MEMB(M,F) = FAM_MEMB(M,F)
         ENDDO
      ENDDO

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass variables to "diag_pl_mod.f" 
      CALL INIT_DIAG_PL( DO_SAVE_PL, DO_SAVE_O3, NFAM,     FAM_NAME,   
     &                   FAM_TYPE,   FAM_NMEM,   FAM_MEMB, FAM_COEF )

      END SUBROUTINE READ_PROD_LOSS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_unix_cmds_menu
!
! !DESCRIPTION: Subroutine READ\_UNIX\_CMDS\_MENU reads the UNIX CMDS MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_UNIX_CMDS_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CHARPAK_MOD,   ONLY : STRSQUEEZE
      USE UNIX_CMDS_MOD, ONLY : BACKGROUND, REDIRECT,  REMOVE_CMD 
      USE UNIX_CMDS_MOD, ONLY : SEPARATOR,  SPACE,     UNZIP_CMD
      USE UNIX_CMDS_MOD, ONLY : WILD_CARD,  ZIP_SUFFIX 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: EOF
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_UNIX_CMDS_MENU begins here!
      !=================================================================

      ! Background
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) BACKGROUND

      ! Redirect
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) REDIRECT

      ! Remove command
      REMOVE_CMD = READ_ONE_LINE( EOF,    'read_unix_cmds_menu:3' ) 
      REMOVE_CMD = REMOVE_CMD(FIRSTCOL:)
      CALL STRSQUEEZE( REMOVE_CMD )

      ! Separator
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) SEPARATOR

      ! Wild Card
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:5' )
      READ( SUBSTRS(1:N), '(a)' ) WILD_CARD

      ! Unzip command
      UNZIP_CMD = READ_ONE_LINE( EOF,     'read_unix_cmds_menu:6' ) 
      UNZIP_CMD = UNZIP_CMD(FIRSTCOL:)
      CALL STRSQUEEZE( UNZIP_CMD )

      ! Zip suffix
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:7' )
      READ( SUBSTRS(1:N), '(a)' ) ZIP_SUFFIX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_unix_cmds_menu:8' )

      ! Just hardwire the SPACE character
      SPACE = ' '

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%BACKGROUND = BACKGROUND
      Input_Opt%REDIRECT   = REDIRECT
      Input_Opt%REMOVE_CMD = REMOVE_CMD
      Input_Opt%SEPARATOR  = SEPARATOR
      Input_Opt%WILD_CARD  = WILD_CARD
      Input_Opt%UNZIP_CMD  = UNZIP_CMD
      Input_Opt%ZIP_SUFFIX = ZIP_SUFFIX

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'UNIX CMDS MENU'
         WRITE( 6, '(  a)' ) '---------------'            
         WRITE( 6, 100     ) 'Unix BACKGROUND  command    : ', 
     &                        TRIM( BACKGROUND )
         WRITE( 6, 100     ) 'Unix REDIRECT    command    : ', 
     &                        TRIM( REDIRECT   )
         WRITE( 6, 100     ) 'Unix REMOVE      command    : ',
     &                        TRIM( REMOVE_CMD )
         WRITE( 6, 100     ) 'Unix SEPARATOR   command    : ',
     &                        TRIM( SEPARATOR  )
         WRITE( 6, 100     ) 'Unix WHITE SPACE command    : ',
     &                        TRIM( SPACE      )
         WRITE( 6, 100     ) 'Unix WILD CARD   command    : ',
     &                        TRIM( WILD_CARD  )
         WRITE( 6, 100     ) 'Unix UNZIP       command    : ',
     &                        TRIM( UNZIP_CMD  )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, A )

      END SUBROUTINE READ_UNIX_CMDS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nested_grid_menu
!
! !DESCRIPTION: Subroutine READ\_NESTED\_GRID\_MENU reads the NESTED GRID MENU 
!  section of the GEOS-CHEM input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_NESTED_GRID_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : TPBC_DIR
      USE DIRECTORY_MOD, ONLY : TPBC_DIR_NA, TPBC_DIR_EU, TPBC_DIR_CH
      USE DIRECTORY_MOD, ONLY : TPBC_DIR_SE
      USE LOGICAL_MOD,   ONLY : LWINDO,      LWINDO2x25,  LWINDO_CU
      USE LOGICAL_MOD,   ONLY : LWINDO_NA,   LWINDO_EU,   LWINDO_CH
      USE LOGICAL_MOD,   ONLY : LWINDO_SE
      USE TPCORE_BC_MOD, ONLY : INIT_TPCORE_BC
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now give user the option of saving out nested grid boundary conditions 
!        at 2 x 2.5 resolution for the EU, CH, or NA grids (amv, bmy, 12/18/09)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I0W, J0W, I1, I2, J1, J2, N, TS
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Save out TPCORE BC's at 4x5
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:1' )
      READ( SUBSTRS(1:N), * ) LWINDO

      ! Save out TPCORE BC's at 2x2.5?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:2' )
      READ( SUBSTRS(1:N), * ) LWINDO2x25

      ! Save out TPCORE BC's over NA
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:3' )
      READ( SUBSTRS(1:N), * ) LWINDO_NA

      ! Directory where NA TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:4' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR_NA

      ! Save out TPCORE BC's over EU
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:5' )
      READ( SUBSTRS(1:N), * ) LWINDO_EU

      ! Directory where EU TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR_EU

      ! Save out TPCORE BC's over CH
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:7' )
      READ( SUBSTRS(1:N), * ) LWINDO_CH

      ! Directory where CH TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:8' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR_CH

      ! Save out TPCORE BC's over the SEAC4RS region
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:9' )
      READ( SUBSTRS(1:N), * ) LWINDO_SE

      ! Directory where SEAC4RS TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:10' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR_SE

      ! Save out TPCORE BC's over Custom Region
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:11' )
      READ( SUBSTRS(1:N), * ) LWINDO_CU

      ! Directory where TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:12' )
      READ( SUBSTRS(1:N), '(a)' ) TPBC_DIR

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:13' )
      READ( SUBSTRS(1:N), * ) TS

      ! Lower left box
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:14' )
      READ( SUBSTRS(1:N), * ) I1, J1

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:15' )
      READ( SUBSTRS(1:N), * ) I2, J2

      ! 1x1 transport region offsets 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_nested_grid_menu:16' )
      READ( SUBSTRS(1:N), * ) I0W, J0W     

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nested_grid_menu:15' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LWINDO      = LWINDO
      Input_Opt%LWINDO2x25  = LWINDO2x25
      Input_Opt%LWINDO_NA   = LWINDO_NA
      Input_Opt%TPBC_DIR_NA = TPBC_DIR_NA
      Input_Opt%LWINDO_EU   = LWINDO_EU
      Input_Opt%TPBC_DIR_EU = TPBC_DIR_EU
      Input_Opt%LWINDO_CH   = LWINDO_CH
      Input_Opt%TPBC_DIR_CH = TPBC_DIR_CH
      Input_Opt%LWINDO_SE   = LWINDO_SE
      Input_Opt%TPBC_DIR_SE = TPBC_DIR_SE
      Input_Opt%LWINDO_CU   = LWINDO_CU
      Input_Opt%TPBC_DIR    = TPBC_DIR
      Input_Opt%NESTED_TS   = TS
      Input_Opt%NESTED_I1   = I1
      Input_Opt%NESTED_J1   = J1
      Input_Opt%NESTED_I2   = I2
      Input_Opt%NESTED_J2   = J2
      Input_Opt%NESTED_I0W  = I0W
      Input_Opt%NESTED_J0W  = J0W

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'NESTED GRID MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100 ) 'Save TPCORE BC''s at 4x5     : ', LWINDO
         WRITE( 6, 110 ) 'Dir w/ archived OH files     : ', 
     &                    TRIM( TPBC_DIR )
         WRITE( 6, 120 ) 'Timestep for 4x5 BC''s [min] : ', TS
         WRITE( 6, 120 ) 'Bot left  box of 4x5 BC area : ', I1,  J1
         WRITE( 6, 120 ) 'Top right box of 4x5 BC area : ', I2,  J2
         WRITE( 6, 120 ) '1 x 1 window offsets         : ', I0W, J0W
      END IF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
 120  FORMAT( A, 2I5 )
    
      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass values to "tpcore_bc_mod.f"
      CALL INIT_TPCORE_BC( TS, I0W, J0W, I1, J1, I2, J2, Input_Opt )

      END SUBROUTINE READ_NESTED_GRID_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_benchmark_menu
!
! !DESCRIPTION: Subroutine READ\_BENCHMARK\_MENU reads the BENCHMARK MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BENCHMARK_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE BENCHMARK_MOD, ONLY : INITIAL_FILE, FINAL_FILE
      USE LOGICAL_MOD,   ONLY : LSTDRUN
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N 
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Save benchmark diagnostic output?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:1' )
      READ( SUBSTRS(1:N), * ) LSTDRUN

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) INITIAL_FILE

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) FINAL_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_benchmark_menu:4' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LSTDRUN           = LSTDRUN
      Input_Opt%STDRUN_INIT_FILE  = INITIAL_FILE
      Input_Opt%STDRUN_FINAL_FILE = FINAL_FILE

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'BENCHMARK MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Save benchmark diag output? : ', LSTDRUN
         WRITE( 6, 110     ) 'File for initial tracer mass: ',  
     &                        TRIM( INITIAL_FILE )
         WRITE( 6, 110     ) 'File for final tracer mass  : ',  
     &                        TRIM( FINAL_FILE )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
    
      END SUBROUTINE READ_BENCHMARK_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_archived_oh_menu( am_I_Root )
!
! !DESCRIPTION: Subroutine READ\_ARCHIVED\_OH\_MENU reads the ARCHIVED OH MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ARCHIVED_OH_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : OH_DIR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ARCHIVED_OH_MENU begins here!
      !=================================================================

      ! Directory where offline OH files are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) OH_DIR

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%OH_DIR = OH_DIR

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ARCHIVED OH MENU'
         WRITE( 6, '(  a)' ) '-----------------'          
         WRITE( 6, 100 ) 'Dir w/ archived OH files    : ', TRIM(OH_DIR)
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, A )
    
      END SUBROUTINE READ_ARCHIVED_OH_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_o3pl_menu
!
! !DESCRIPTION: Subroutine READ\_O3PL\_MENU reads the O3 P/L MENU section 
!  of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_O3PL_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : O3PL_DIR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_O3PL_MENU begins here!
      !=================================================================

      ! Directory where archived P(O3) and L(O3) files are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) O3PL_DIR

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_archived_oh_menu:2' )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%O3PL_DIR = O3PL_DIR

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ARCHIVED O3PL MENU'
         WRITE( 6, '(  a)' ) '-------------------'          
         WRITE( 6, 100     ) 'Dir w/ archived O3 P/L files: ', 
     &                       TRIM( O3PL_DIR )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, A )
    
      END SUBROUTINE READ_O3PL_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_mercury_menu
!
! !DESCRIPTION: Subroutine READ\_MERCURY\_MENU reads the BENCHMARK MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_MERCURY_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      ! References to F90 modules
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN, LPREINDHG, LGTMM
      USE MERCURY_MOD,       ONLY : INIT_MERCURY
      USE OCEAN_MERCURY_MOD, ONLY : INIT_OCEAN_MERCURY
      USE DEPO_MERCURY_MOD,  ONLY : INIT_DEPO_MERCURY
      USE LAND_MERCURY_MOD,  ONLY : INIT_LAND_MERCURY
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  24 Feb 2006 - R. Yantosca - Initial version
!  ( 1) Update for Chris Holmes's mercury version. (ccc, 5/6/10)
!  ( 2) Add options to use GTMM for mercury soil emissions (ccc, 9/16/09)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: USE_CHECKS
      INTEGER            :: ANTHRO_Hg_YEAR,  N 
      CHARACTER(LEN=255) :: HG_SCENARIO  !eds 5/17/12
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), Hg_RST_FILE
      CHARACTER(LEN=255) :: GTMM_RST_FILE

      !=================================================================
      ! READ_MERCURY_MENU begins here!
      !=================================================================

      ! Year for anthro Hg emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:1' )
      READ( SUBSTRS(1:N), * ) ANTHRO_Hg_YEAR

      ! Use future emissions scenario?   !eds 5/17/12
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) HG_SCENARIO

      ! Use error check for tag/tot Hg?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:3' )
      READ( SUBSTRS(1:N), * ) USE_CHECKS

      ! Use dynamic ocean model?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:4' )
      READ( SUBSTRS(1:N), * ) LDYNOCEAN

      ! Use preindustrial simulation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read mercury_menu:5' )
      READ( SUBSTRS(1:N), * ) LPREINDHG

      ! Name of ocean restart file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:6' )
      READ( SUBSTRS(1:N), '(a)' ) Hg_RST_FILE

      ! Use GTMM ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:7' )
      READ( SUBSTRS(1:N), * ) LGTMM

      ! Name of GTMM restart file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:8' )
      READ( SUBSTRS(1:N), * ) GTMM_RST_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_mercury_menu:9' )

      ! Check on logical
      IF (.NOT.( ITS_A_MERCURY_SIM() ) ) THEN
         LGTMM     = .FALSE.
         LDYNOCEAN = .FALSE.
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'MERCURY MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 100 ) 'Anthro Hg emissions year    : ', 
     &                    ANTHRO_Hg_YEAR
         WRITE( 6, 120 ) 'Future emissions scenario   : ',  !eds 5/17/12
     &                    TRIM( HG_SCENARIO )
         WRITE( 6, 110 ) 'Error check tag & total Hg? : ', USE_CHECKS
         WRITE( 6, 110 ) 'Use dynamic ocean Hg model? : ', LDYNOCEAN
         WRITE( 6, 110 ) 'Preindustrial simulation?   : ', LPREINDHG
         WRITE( 6, 120 ) 'Ocean Hg restart file       : ',
     &                    TRIM( Hg_RST_FILE )
         WRITE( 6, 110 ) 'Use GTMM ?                  : ', LGTMM
         WRITE( 6, 120 ) '=> GTMM restart file        : ',
     &                        TRIM( GTMM_RST_FILE )
      END IF

      ! FORMAT statements
 100  FORMAT( A, I4  )
 110  FORMAT( A, L5  )
 120  FORMAT( A, A   )
    
      !=================================================================
      ! Copy to Input Options object
      !=================================================================
      Input_Opt%ANTHRO_Hg_YEAR = ANTHRO_Hg_YEAR
      Input_Opt%HG_SCENARIO    = HG_SCENARIO
      Input_Opt%USE_CHECKS     = USE_CHECKS
      Input_Opt%LDYNOCEAN      = LDYNOCEAN
      Input_Opt%LPREINDHG      = LPREINDHG
      Input_Opt%Hg_RST_FILE    = Hg_RST_FILE
      Input_Opt%LGTMM          = LGTMM
      Input_Opt%GTMM_RST_FILE  = GTMM_RST_FILE

      ! Return success
      RC = GIGC_SUCCESS

      ! If we are performing a Hg simulation ...
      IF ( ITS_A_MERCURY_SIM() ) THEN 

         ! Initialize "mercury_mod.f"
         CALL INIT_MERCURY( am_I_Root, Input_Opt, RC )
         CALL INIT_LAND_MERCURY()
         CALL INIT_DEPO_MERCURY(TRIM( GTMM_RST_FILE ))

         ! Initialize "ocean_mercury_mod.f"
         IF ( LDYNOCEAN ) THEN
            CALL INIT_OCEAN_MERCURY( Hg_RST_FILE, USE_CHECKS )
         ENDIF

      ENDIF

      END SUBROUTINE READ_MERCURY_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_ch4_menu
!
! !DESCRIPTION: Subroutine READ\_CH4\_MENU reads the CH4 MENU section of 
!  the GEOS-Chem input file; this defines emissions options for CH4 tagged 
!  simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CH4_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      ! References to F90 modules
      USE LOGICAL_MOD, ONLY : LGAO,    LCOL,    LLIV,   LWAST
      USE LOGICAL_MOD, ONLY : LBFCH4,  LBMCH4,  LWETL,  LRICE
      USE LOGICAL_MOD, ONLY : LOTANT,  LSOABS,  LOTNAT
      USE LOGICAL_MOD, ONLY : LCH4BUD
      ! kjw
      USE LOGICAL_MOD, ONLY : LBIOMASS
      USE LOGICAL_MOD, ONLY : LGFED3BB, LDAYBB3
      USE LOGICAL_MOD, ONLY : LGFED2BB, L8DAYBB
      USE ERROR_MOD,   ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  03 Aug 2009 - K. Wecht, C. Pickett-Heaps - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 Feb 2014 - R. Yantosca - Add warning for CH4 budget (which is now
!                              controlled by an #ifdef in global_ch4_mod.F)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      CHARACTER(LEN=255) :: MSG

      !=================================================================
      ! READ_CH4_MENU begins here!
      !=================================================================

      ! Compute CH4 budget
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:1' )
      READ( SUBSTRS(1:N), * ) LCH4BUD

      ! Use Gas & Oil emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:2' )
      READ( SUBSTRS(1:N), * ) LGAO

      ! Use Coal emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:3' )
      READ( SUBSTRS(1:N), * ) LCOL              
                                                
      ! Use Livestock emissions?                
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:4' )
      READ( SUBSTRS(1:N), * ) LLIV              
                                                
      ! Use Waste emissions?                    
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:5' )
      READ( SUBSTRS(1:N), * ) LWAST             
                                                
      ! Use Biofuel emissions?                  
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:6' )
      READ( SUBSTRS(1:N), * ) LBFCH4          
                                                
      ! Use Rice emissions?                     
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:7' )
      READ( SUBSTRS(1:N), * ) LRICE             
                                                
      ! Use Other Anthropogenic emissions?      
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:8' )
      READ( SUBSTRS(1:N), * ) LOTANT            
                                                
      ! Use Biomass emissions?                  
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:9' )
      READ( SUBSTRS(1:N), * ) LBMCH4          
                                                
      ! Use Wetlands emissions?                 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:10' )
      READ( SUBSTRS(1:N), * ) LWETL             
                                                
      ! Use Soil Absorption?          
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:11' )
      READ( SUBSTRS(1:N), * ) LSOABS            
                                                
      ! Use Other Natural emissions?            
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:12' )
      READ( SUBSTRS(1:N), * ) LOTNAT            
                                                
      ! Separator line                          
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_CH4_menu:13' )

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CH4 MENU'
         WRITE( 6, '(  a)' ) '-----------'
         WRITE( 6, 100     ) 'Compute CH4 budget?   : ', LCH4BUD
         WRITE( 6, 100     ) 'Use Gas & Oil emis?   : ', LGAO
         WRITE( 6, 100     ) 'Use Coal Mine emis?   : ', LCOL
         WRITE( 6, 100     ) 'Use Livestock emis?   : ', LLIV
         WRITE( 6, 100     ) 'Use Waste emis?       : ', LWAST
         WRITE( 6, 100     ) 'Use Biofuel emis?     : ', LBFCH4
         WRITE( 6, 100     ) 'Use Rice emis?        : ', LRICE
         WRITE( 6, 100     ) 'Use Ot. Anthro emis?  : ', LOTANT
         WRITE( 6, 100     ) 'Use Biomass emis?     : ', LBMCH4
         WRITE( 6, 100     ) 'Use Wetlands emis?    : ', LWETL
         WRITE( 6, 100     ) 'Use Soil Absorption?  : ', LSOABS
         WRITE( 6, 100     ) 'Use Ot. Natural emis? : ', LOTNAT
      END IF

      ! Error check (kjw)
      IF ( LBMCH4 .EQ. .TRUE. ) THEN
         IF ( ( LBIOMASS .EQ. .FALSE. ) .OR. 
     &        ( LGFED2BB .EQ. .FALSE. .AND. LGFED3BB .EQ. .FALSE.  .AND.  
     &          L8DAYBB  .EQ. .FALSE. .AND. LDAYBB3  .EQ. .FALSE.)) THEN
            MSG = 'Specify which biomass burning emissions to use.'  //
     &            '  If necessary, read CH4_MENU after EMISSIONS_MENU.'
            IF( am_I_Root ) THEN
               print*,'LBIOMASS = ',LBIOMASS
               print*,'LGFED3BB = ',LGFED3BB
               print*,'LGFED2BB = ',LGFED2BB
               print*,'L8DAYBB  = ',L8DAYBB
!              print*,'L8DAYBB3 = ',L8DAYBB3
	    ENDIF
            CALL ERROR_STOP( MSG, 'READ_CH4_MENU ("input_mod.F")' )
         ENDIF
      ENDIF

      ! Error check (kjw)
      IF ( LCH4BUD .EQ. .TRUE. ) THEN
         MSG = 'To print CH4 budget information, be sure that ' //
     &         'USE_CH4_BUDGET_DIAG is defined at the beginning' //
     &         ' of global_ch4_mod.F'
         IF ( am_I_Root ) THEN
            print*,'CH4 Budget warning:'
            print*,MSG
         ENDIF
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================
      Input_Opt%LCH4BUD = LCH4BUD
      Input_Opt%LGAO    = LGAO
      Input_Opt%LCOL    = LCOL         
      Input_Opt%LLIV    = LLIV         
      Input_Opt%LWAST   = LWAST         
      Input_Opt%LBFCH4  = LBFCH4       
      Input_Opt%LRICE   = LRICE         
      Input_Opt%LOTANT  = LOTANT         
      Input_Opt%LBMCH4  = LBMCH4       
      Input_Opt%LWETL   = LWETL         
      Input_Opt%LSOABS  = LSOABS          
      Input_Opt%LOTNAT  = LOTNAT 

      ! Return success
      RC = GIGC_SUCCESS         

      END SUBROUTINE READ_CH4_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_pops_menu
!
! !DESCRIPTION: Subroutine READ\_POPS\_MENU reads the POPS MENU section of 
!  the GEOS-Chem input file; this defines emissions options for POPs
!  simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_POPS_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE GET_POPSINFO_MOD,     ONLY : GET_POP_TYPE
      USE POPS_MOD,             ONLY : INIT_POPS
      USE DIRECTORY_MOD,        ONLY : POP_EMISDIR
      USE TRACER_MOD,           ONLY : ITS_A_POPS_SIM
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)           :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  01 Oct 2012 - C. Friedman - Initial version
!  26 Nov 2012 - M. Payer    - Added ProTeX headers
!  29 Nov 2012 - M. Payer    - Now pass Input_Opt, RC as arguments
!  26 Mar 2013 - R. Yantosca - Now pass Input_Opt to INIT_POPS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: N
      LOGICAL                    :: CHEM_PROCESS
      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)
      CHARACTER(LEN=3)           :: POP_TYPE
      REAL*8                     :: POP_XMW, POP_KOA, POP_KBC
      REAL*8                     :: POP_K_POPG_OH
      REAL*8                     :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL*8                     :: POP_HSTAR, POP_DEL_H, POP_DEL_Hw

      !=================================================================
      ! READ_POPS_MENU begins here!
      !=================================================================

      ! POP species
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:1' )
      READ( SUBSTRS(1:N), * ) POP_TYPE

      ! Dummy for future process logical switches
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:2' )
      READ( SUBSTRS(1:N), * ) CHEM_PROCESS

      ! Name of pops emissions file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) POP_EMISDIR

      ! Molecular weight
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:4' )
      READ( SUBSTRS(1:N), * ) POP_XMW

      ! KOA
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:5' )
      READ( SUBSTRS(1:N), * ) POP_KOA

      ! KBC
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:6' )
      READ( SUBSTRS(1:N), * ) POP_KBC

      ! OH oxidation
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:7' )
      READ( SUBSTRS(1:N), * ) POP_K_POPG_OH

      ! O3 oxidation 1
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:8' )
      READ( SUBSTRS(1:N), * ) POP_K_POPP_O3A

      ! O3 oxidation 2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:9' )
      READ( SUBSTRS(1:N), * ) POP_K_POPP_O3B

      ! H*
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:10' )
      READ( SUBSTRS(1:N), * ) POP_HSTAR

      ! DEL_H
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:11' )
      READ( SUBSTRS(1:N), * ) POP_DEL_H

      ! DEL_Hw
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:12' )
      READ( SUBSTRS(1:N), * ) POP_DEL_Hw

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_pops_menu:13' )

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'POPS MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 120     ) 'Species of POP    : ', 
     &                        POP_TYPE
         WRITE( 6, 110     ) 'Chemistry on? : ', CHEM_PROCESS
         WRITE( 6, 120     ) '=> POP Emissions file       : ',
     &                        TRIM( POP_EMISDIR )
         WRITE( 6, 130     ) 'Param 1:', POP_XMW
      ENDIF

      ! FORMAT statements
 110  FORMAT( A, L5  )
 120  FORMAT( A, A   )
 130  FORMAT( A, f6.5 )

      !=================================================================
      ! Copy to Input Options object
      !=================================================================
      Input_Opt%POP_TYPE       = POP_TYPE
      Input_Opt%CHEM_PROCESS   = CHEM_PROCESS
      Input_Opt%POP_EMISFILE   = POP_EMISDIR
      Input_Opt%POP_XMW        = POP_XMW
      Input_Opt%POP_KOA        = POP_KOA
      Input_Opt%POP_KBC        = POP_KBC
      Input_Opt%POP_K_POPG_OH  = POP_K_POPG_OH
      Input_Opt%POP_K_POPP_O3A = POP_K_POPP_O3A
      Input_Opt%POP_K_POPP_O3B = POP_K_POPP_O3B
      Input_Opt%POP_HSTAR      = POP_HSTAR
      Input_Opt%POP_DEL_H      = POP_DEL_H 
      Input_Opt%POP_DEL_Hw     = POP_DEL_Hw

      ! Return success
      RC = GIGC_SUCCESS

      ! If we are performing a POPS simulation ...
      IF ( ITS_A_POPS_SIM() ) THEN 

         ! Initially set pop type
         POP_TYPE = GET_POP_TYPE( POP_TYPE )

         ! Set emissions file to use
!         POP_EMISSFILE = GET_EMISSFILE( TRIM( POP_EMISSFILE ) )

         ! Initialize "pops_mod.f"
         CALL INIT_POPS(POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH,
     &                  POP_K_POPP_O3A, POP_K_POPP_O3B, 
     &                  POP_HSTAR, POP_DEL_H, POP_DEL_Hw, Input_Opt )

      ENDIF

      END SUBROUTINE READ_POPS_MENU
!EOC
#if   defined( APM ) 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_apm_menu
!
! !DESCRIPTION: Subroutine READ\_APM\_MENU reads the APM MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_APM_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE APM_INIT_MOD,       ONLY : APMTRACER_MW_G
      USE APM_INIT_MOD,       ONLY : APMTRACER_MW_Kg
      USE APM_INIT_MOD,       ONLY : IFNUCL
      USE APM_INIT_MOD,       ONLY : FE0
      USE APM_INIT_MOD,       ONLY : LAPM
      USE CHARPAK_MOD,        ONLY : STRREPL
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE FILE_MOD,           ONLY : IOERROR
      USE GAMAP_MOD,          ONLY : DO_GAMAP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACER_MOD,         ONLY : N_APMTRA
      USE TRACER_MOD,         ONLY : N_TRACERS
      USE TRACER_MOD,         ONLY : TCVV
      USE TRACER_MOD,         ONLY : XNUMOL
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  This subroutine is only compiled when you build GEOS-Chem with the
!  APM=yes makefile option.
!
! !REVISION HISTORY: 
!  30 Sep 2008 - G. Luo, F. Yu - Initial version
!  16 Feb 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N,               N_APMTRAC1, T
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM), MSG

      !=================================================================
      ! READ_APM_MENU begins here!
      !=================================================================

      ! Error check
      IF ( CT1 /= 2 ) THEN
         MSG = 'SIMULATION MENU & TRACER MENU must be read in first!'
         CALL ERROR_STOP( MSG, 'READ_APM_MENU ("input_mod.f")' )
      ENDIF

      ! Nucleation scheme
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_apm_menu:1' )
      READ( SUBSTRS(1:N), * ) IFNUCL

      ! Primary sulfate fraction
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_apm_menu:2' )
      READ( SUBSTRS(1:N), * ) FE0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_apm_menu:3' )

      ! Increment the TCVV and XNUMOL arrays for APM tracers
      DO T = 1, N_APMTRA
         TCVV(N_TRACERS+T)   = 28.97d0   / APMTRACER_MW_G(T)
         XNUMOL(N_TRACERS+T) = 6.022d+23 / APMTRACER_MW_KG(T) 
      ENDDO

      !=================================================================
      ! Set fields of the Input Options object
      !=================================================================
      Input_Opt%IFNUCL = IFNUCL
      Input_Opt%FE0    = FE0

      DO T = 1, N_APMTRA
         Input_Opt%TCVV  (N_TRACERS+T) = TCVV(T)
         Input_Opt%XNUMOL(N_TRACERS+T) = XNUMOL(T)
      ENDDO

      ! Return success
      RC = GIGC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'APM MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100     ) 'Turn on APM and add tracers?: ', LAPM
         WRITE( 6, 110     ) '# of additional APM tracers : ', N_APMTRA
         WRITE( 6, 110     ) 'Nucleation Scheme           : ', IFNUCL
         WRITE( 6, 120     ) 'Primary Sulfate Fraction    : ', FE0
      ENDIF

      ! FORMATs
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
 120  FORMAT( A, F7.3 )

      END SUBROUTINE READ_APM_MENU

!EOC
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: validate_directories
!
! !DESCRIPTION: Subroutine VALIDATE\_DIRECTORIES makes sure that each of the 
!  directories that we have read from the GEOS-Chem input file are valid.  
!  Also, trailing separator characters will be added.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE VALIDATE_DIRECTORIES
!
! !USES:
!
      USE DIRECTORY_MOD
      USE GRID_MOD,      ONLY : ITS_A_NESTED_GRID 
      USE LOGICAL_MOD,   ONLY : LWINDO_CU,   LUNZIP
      USE LOGICAL_MOD,   ONLY : LWINDO_NA,   LWINDO_EU,    LWINDO_CH
      USE LOGICAL_MOD,   ONLY : LWINDO_SE
      USE TIME_MOD,      ONLY : EXPAND_DATE, GET_NYMDb,    GET_NYMDe
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY.  Now also validate
!        GCAP and GEOS-5 directories. (bmy, 10/3/05)
!  (2 ) Now references DATA_DIR_1x1 from directory_mod.f (bmy, 10/24/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (4 ) Now check TPBC_DIR_NA, TPBC_DIR_CH, TPBC_DIR_EU (amv, bmy, 12/18/09)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  27 Aug 2010 - R. Yantosca -  Now check MERRA directory
!  08 Feb 2012 - R. Yantosca - Now check GEOS-5.7.x directory
!  09 Feb 2012 - R. Yantosca - Rewrote USE statements for clarity
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: NYMDb, NYMDe
      CHARACTER(LEN=255) :: DIR

      !=================================================================
      ! VALIDATE_DIRECTORIES begins here!
      !=================================================================

      ! Get starting & ending dates
      NYMDb = GET_NYMDb()
      NYMDe = GET_NYMDe()

      ! Check directories
      CALL CHECK_DIRECTORY( DATA_DIR     )
      CALL CHECK_DIRECTORY( DATA_DIR_1x1 )
      CALL CHECK_DIRECTORY( RUN_DIR      )
      CALL CHECK_DIRECTORY( OH_DIR       )
      CALL CHECK_DIRECTORY( O3PL_DIR     )
!      CALL CHECK_DIRECTORY( OC_DIR       )  !(clf, 1/20/2011)
!      CALL CHECK_DIRECTORY( BC_DIR       )  !(clf, 1/20/2011)

      ! Only validate the TEMP_DIR if we are unzipping met fields
      IF ( LUNZIP ) CALL CHECK_DIRECTORY( TEMP_DIR )

      ! Only validate the TPCORE BC directory if we need it
      IF ( ITS_A_NESTED_GRID() ) THEN
#if   defined( NESTED_NA)
         CALL CHECK_DIRECTORY( TPBC_DIR_NA )
#elif defined( NESTED_EU)
         CALL CHECK_DIRECTORY( TPBC_DIR_EU )
#elif defined( NESTED_CH)
         CALL CHECK_DIRECTORY( TPBC_DIR_CH )
#endif

      ENDIF
      IF ( LWINDO_CU ) THEN
         CALL CHECK_DIRECTORY( TPBC_DIR )
      ENDIF
      IF ( LWINDO_NA ) THEN
         CALL CHECK_DIRECTORY( TPBC_DIR_NA )
      ENDIF
      IF ( LWINDO_EU ) THEN
         CALL CHECK_DIRECTORY( TPBC_DIR_EU )
      ENDIF
      IF ( LWINDO_CH ) THEN
         CALL CHECK_DIRECTORY( TPBC_DIR_CH )
      ENDIF

#if   defined( GEOS_4 )

      ! Check GEOS-4 met field directory (starting date)
      DIR = GEOS_4_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )


      ! Check GEOS-4 met field directory (ending date)
      DIR = GEOS_4_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GEOS_5 )

      ! Check GEOS-5 met field directory (starting date)
      DIR = GEOS_5_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-5 met field directory (ending date)
      DIR = GEOS_5_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GEOS_FP )

      ! Check GEOS-FP met field directory (starting date)
      DIR = GEOS_FP_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-FP met field directory (ending date)
      DIR = GEOS_FP_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( MERRA )

      ! Check GEOS-5 met field directory (starting date)
      DIR = MERRA_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

      ! Check GEOS-5 met field directory (ending date)
      DIR = MERRA_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#elif defined( GCAP )

      ! Check GEOS-5 met field directory (starting date)
      DIR = GCAP_DIR
      CALL EXPAND_DATE( DIR, NYMDb, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )


      ! Check GEOS-5 met field directory (ending date)
      DIR = GCAP_DIR
      CALL EXPAND_DATE( DIR, NYMDe, 000000 )
      DIR = TRIM( DATA_DIR ) // TRIM( DIR )
      CALL CHECK_DIRECTORY( DIR )

#endif

      END SUBROUTINE VALIDATE_DIRECTORIES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_directory
!
! !DESCRIPTION: Subroutine CHECK\_DIRECTORY makes sure that the given 
!  directory is valid.  Also a trailing slash character will be added if 
!  necessary. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_DIRECTORY( DIR )
!
! !USES:
!
      ! References to F90 modules 
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE FILE_MOD,      ONLY : FILE_EXISTS
      USE UNIX_CMDS_MOD, ONLY : SEPARATOR
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(INOUT) :: DIR    ! Directory to be checked
! 
! !REVISION HISTORY: 
!  20 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: C
      CHARACTER(LEN=255) :: MSG
      
      !=================================================================
      ! CHECK_DIRECTORY begins here!
      !=================================================================

      ! Locate the last non-white-space character of NEWDIR
      C = LEN_TRIM( DIR )

      ! Add the trailing directory separator if it is not present
      IF ( DIR(C:C) /= TRIM( SEPARATOR ) ) THEN 
         DIR(C+1:C+1) = TRIM( SEPARATOR )
      ENDIF
     
      !=================================================================
      ! Test if the directory actually exists
      !=================================================================

      ! If the directory does not exist then stop w/ an error message
      IF ( .not. FILE_EXISTS( DIR ) ) THEN 
         MSG = 'Invalid directory: ' // TRIM( DIR ) 
         CALL ERROR_STOP( MSG, 'CHECK_DIRECTORY ("input_mod.f")' )
      ENDIF

      END SUBROUTINE CHECK_DIRECTORY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_time_steps
!
! !DESCRIPTION: Subroutine CHECK\_TIME\_STEPS computes the smallest dynamic 
!  time step for the model, based on which operation are turned on.  This 
!  is called from routine READ\_INPUT\_FILE, after all of the timesteps and 
!  logical flags have been read from "input.geos".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_TIME_STEPS( am_I_Root )
!
! !USES:
!
      USE LOGICAL_MOD, ONLY : LCONV, LCHEM, LDRYD 
      USE LOGICAL_MOD, ONLY : LEMIS, LTRAN, LTURB 
      USE TIME_MOD,    ONLY : SET_TIMESTEPS
      USE ERROR_MOD,   ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD,  ONLY : ITS_A_CH4_SIM
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Add TS_DIAG, the largest time steps used for diagnostics.
!        And test that all time steps are multiple of the smallest one.
!        (ccc, 5/13/09)
!  (3 ) Corrected typos -99999 instead of -999999 (phs, bmy, 8/21/09)
!  (4 ) Now compute TS_SUN_2 which is 1/2 of the chemistry timestep (or
!        smallest timestep if LCHEM=LEMIS=LDRYD=F).  This is used to compute
!        SUNCOS at the midpoint of the timestep instead of the beginning.
!        (bmy, 4/27/10)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  07 Oct 2011 - R. Yantosca - Add extra error checks for centralizing
!                              chemistry timestep algorithm
!  07 Oct 2011 - R. Yantosca - Remove TS_SUN_2 from call to SET_TIMESTEPS
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, K, L, TS_SMALLEST, TS_DIAG
      
      !=================================================================
      ! CHECK_TIME_STEPS begins here!
      !=================================================================

      ! NUNIT is time step in minutes for unit conversion
      TS_UNIT = -1

      ! For the centralizing chemistry timestep modification, the emissions 
      ! should be done at the same time as chemistry.  Stop the run if the
      ! emissions timestep differs from the chemistry timestep. (bmy, 10/7/11)
      IF ( TS_EMIS /= TS_CHEM ) THEN 
         IF ( am_I_Root ) THEN
          WRITE(6,*) 'The emission time step should equal the chemistry'
          WRITE(6,*) 'timestep.  This is required for the new central'
          WRITE(6,*) 'chemistry timestep algorithm.'
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      ! For the centralizing chemistry timestep modification, the convection 
      ! should be done at the same time as dynamics.  Stop the run if the
      ! convection timestep differs from the dynamic timestep. (bmy, 10/7/11)
      IF ( TS_CONV /= TS_DYN ) THEN 
         IF ( am_I_Root ) THEN
          WRITE(6,*) 'The convection time step should equal the dynamic'
          WRITE(6,*) 'timestep.  This is required for the new central'
          WRITE(6,*) 'chemistry timestep algorithm.'
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Only do unit conversion if necessary
      IF ( LTRAN .or. LCONV .or. LTURB ) THEN
         TS_UNIT = MAX( TS_DYN, TS_CONV )
      ENDIF

      ! Compute NSMALLEST as the minimum of NDYN, NCONV, NSRCE, NCHEM
      I = TS_DYN
      J = TS_CONV
      K = TS_EMIS
      L = TS_CHEM

      IF ( .not. LTRAN                  ) I = 999999 
      IF ( .not. LCONV .and..not. LTURB ) J = 999999
      IF ( .not. LDRYD .and..not. LEMIS ) K = 999999
      IF ( .not. LCHEM                  ) L = 999999

      ! Get the smallest of all of the above
      TS_SMALLEST = MIN( I, J, K, L )

      ! If all of the operators above are turned off, 
      ! then set TS_SMALLEST to TS_DYN.
      IF ( TS_SMALLEST == 999999 ) THEN 
         TS_SMALLEST = TS_DYN
      ENDIF
       
      IF ( LTRAN .and. TS_DYN /= TS_SMALLEST ) THEN
         IF ( am_I_Root ) THEN
         WRITE(6,*) 'The transport time step should be the smallest one'
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF
     
      ! If TS_DYN is smaller than TS_SMALLEST, reset TS_DYN
      ! to TS_SMALLEST.
      ! This is useful for runs where transport is turned off,
      ! but where chemistry is turned on. 
      IF ( TS_DYN < TS_SMALLEST ) THEN
         TS_DYN = TS_SMALLEST
      ENDIF

      ! Define the largest time step, TS_DIAG, for diagnostics.
      ! Diagnostics should be incremented at the end of multiples of
      ! TS_DIAG, so that the system is at a physical state.
      ! (ccc, 5/13/09)
      IF ( .not. LTRAN                  ) I = -999999 
      IF ( .not. LCONV .and..not. LTURB ) J = -999999
      IF ( .not. LDRYD .and..not. LEMIS ) K = -999999
      IF ( .not. LCHEM                  ) L = -999999

      TS_DIAG = MAX( I, J, K, L )

      ! If all the operators are turned off, then set TS_DIAG to TS_CHEM
      ! Usually the chemistry time step is large. (ccc, 5/13/09)
      IF ( TS_DIAG == -999999 ) THEN
         TS_DIAG = TS_CHEM
      ENDIF

      ! Change TS_DIAG to TS_DYN for methane simulation. (ccc, 8/27/09)
      IF ( ITS_A_CH4_SIM() ) TS_DIAG = TS_DYN

      ! Check if all time steps are multiples of the smallest.
      ! (ccc, 5/13/09)
      IF ( L /= -999999 .and. MOD( TS_CHEM, TS_SMALLEST ) /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) 'Chemistry', TS_CHEM, TS_SMALLEST
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF
      
      IF ( K /= -999999 .and. MOD( TS_EMIS, TS_SMALLEST ) /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) 'Emission', TS_EMIS, TS_SMALLEST
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      IF ( J /= -999999 .and. MOD( TS_CONV, TS_SMALLEST ) /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) 'Convection', TS_CONV, TS_SMALLEST
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      IF ( I /= -999999 .and. MOD( TS_DYN, TS_SMALLEST ) /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE( 6, 100 ) 'Transport', TS_DYN, TS_SMALLEST
         ENDIF
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Initialize timesteps in "time_mod.f"
      CALL SET_TIMESTEPS( am_I_Root,
     &                    CHEMISTRY  = TS_CHEM, EMISSION  = TS_EMIS, 
     &                    DYNAMICS   = TS_DYN,  UNIT_CONV = TS_UNIT,
     &                    CONVECTION = TS_CONV, DIAGNOS   = TS_DIAG)

 100  FORMAT( A, ' time step must be a multiple of the smallest one:',
     &        i5, i5 )

      END SUBROUTINE CHECK_TIME_STEPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_last_day_good
!
! !DESCRIPTION: Suborutine IS\_LAST\_DAY\_GOOD tests to see if there is 
!  output scheduled on the last day of the run. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE IS_LAST_DAY_GOOD
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE JULDAY_MOD, ONLY : JULDAY
      USE TIME_MOD,   ONLY : GET_NYMDe, ITS_A_LEAPYEAR, YMD_EXTRACT

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! NJDAY
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Moved to "input_mod.f" from "main.f" (bmy, 1/11/05)
!  (2 ) Now call ITS_A_LEAPYEAR with FORCE=.TRUE. to always return whether
!        the year Y would be a leap year, regardless of met field type.
!        (swu, bmy, 4/24/06)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_LEAPYEAR
      INTEGER :: NYMDe, Y, M, D, LASTDAY
      REAL*8  :: JD, JD0

      !=================================================================
      ! IS_LAST_DAY_GOOD begins here!
      !=================================================================

      ! Astronomical Julian Day corresponding to NYMDe
      NYMDe = GET_NYMDe()
      CALL YMD_EXTRACT( NYMDe, Y, M, D )
      JD = JULDAY( Y, M, DBLE( D ) )

      ! Astronomical Julian Day corresponding to the 1st of the year
      JD0 = JULDAY( Y, 1, 0d0 )

      ! LASTDAY is the day of year corresponding to NYMDe      
      LASTDAY = JD - JD0

      ! Skip past the element of NJDAY for Feb 29, if necessary
      IF ( .not. ITS_A_LEAPYEAR( Y, .TRUE. ) .and. LASTDAY > 59 ) THEN
         LASTDAY = LASTDAY + 1
      ENDIF

      ! Stop w/ error if THIS_NJDAY = 0 
      IF ( NJDAY(LASTDAY) == 0 ) THEN
         CALL ERROR_STOP( 'No output scheduled on last day of run!',
     &                    'IS_LAST_DAY_GOOD ("input_mod.f")' )
      ENDIF
     
      END SUBROUTINE IS_LAST_DAY_GOOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_extra
!
! !DESCRIPTION: Suborutine GIGC\_INIT\_EXTRA initializes other GEOS-Chem 
!  modules that have not been initialized in either GIGC\_Allocate\_All or
!  GIGC\_Init\_all.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GIGC_Init_Extra( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE Aerosol_Mod,        ONLY : Init_Aerosol
      USE Carbon_Mod,         ONLY : Init_Carbon
      USE Drydep_Mod,         ONLY : Init_Drydep
      USE Dust_Mod,           ONLY : Init_Dust
      USE Error_Mod,          ONLY : Debug_Msg
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE Linoz_Mod,          ONLY : Init_Linoz
      USE Seasalt_Mod,        ONLY : Init_SeaSalt
      USE Sulfate_Mod,        ONLY : Init_Sulfate
      USE Tropopause_Mod,     ONLY : Init_Tropopause
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Several of the INIT routines now called within GIGC_Init_Extra had 
!  originally been called from the Run method.  We now gather these INIT
!  routines here so that they may be called from the Initialization method.
!  This is necessary when connecting GEOS-Chem to the GEOS-5 GCM via ESMF.
!                                                                             .
!  GIGC_Init_Extra should be called after the call to DO_DRYDEP, since
!  these routines depend on dry deposition parameters being set up first.
!
! !REVISION HISTORY: 
!  04 Mar 2013 - R. Yantosca - Initial revision
!  05 Mar 2013 - R. Yantosca - Now call INIT_AEROSOL (GeosCore/aerosol_mod.F)
!  15 Mar 2013 - R. Yantosca - Now call INIT_LINOZ (GeosCore/linoz_mod.F)
!  29 Mar 2013 - R. Yantosca - Now call INIT_TROPOPAUSE (so that we can pass
!                              a LVARTROP from Input_Opt and not logical_mod.F)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: prtDebug

      !=================================================================
      ! GIGC_Init_Extra begins here!
      !=================================================================

      ! Assume success
      RC       = GIGC_SUCCESS

      ! Set a flag for debugging
      prtDebug = ( am_I_Root .and. Input_Opt%LPRT )

      ! Initialize "carbon_mod.F"
      IF ( Input_Opt%LCARB ) THEN
         CALL Init_Carbon( am_I_Root, Input_Opt, RC )
         IF ( RC /= GIGC_SUCCESS ) RETURN
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_CARBON' )
      ENDIF
    
      ! Initialize "dust_mod.F"
      IF ( Input_Opt%LDUST ) then
         CALL Init_Dust( am_I_Root, Input_Opt, RC )
         IF ( RC /= GIGC_SUCCESS ) RETURN
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_DUST' )
      ENDIF

      ! Initialize "seasalt_mod.F
      IF ( Input_Opt%LSSALT ) THEN
         CALL Init_Seasalt( am_I_Root, Input_Opt, RC )
         IF ( RC /= GIGC_SUCCESS ) RETURN
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_SEASALT' )
      ENDIF

      ! Initialize "sulfate_mod.F"
      IF ( Input_Opt%LSULF ) THEN
         CALL Init_Sulfate( am_I_Root, Input_Opt, RC )
         IF ( RC /= GIGC_SUCCESS ) RETURN
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_SULFATE' )
      ENDIF

      ! Initialize "aerosol_mod.F"
      IF ( Input_Opt%LSULF .or. Input_Opt%LCARB   .or. 
     &     Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN
         CALL Init_Aerosol( am_I_Root, Input_Opt, RC )
         IF ( RC /= GIGC_SUCCESS ) RETURN
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_AEROSOL' )
      ENDIF

      ! Initialize "linoz_mod.F"
      IF ( Input_Opt%LLINOZ ) THEN
         CALL INIT_LINOZ( am_I_Root, Input_Opt, RC )
         IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_LINOZ' )
      ENDIF

      ! Initialize "tropopause_mod.F"
      CALL INIT_TROPOPAUSE( am_I_Root, Input_Opt, RC )
      IF ( prtDebug ) CALL DEBUG_MSG( '### a INIT_TROPOPAUSE' )

      END SUBROUTINE GIGC_Init_Extra
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_input
!
! !DESCRIPTION: Subroutine INIT\_INPUT initializes all variables from 
!  "directory\_mod.f" and "logical\_mod.f" for safety's sake.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_INPUT
!
! !USES:
!
      USE DIRECTORY_MOD
      USE LOGICAL_MOD
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now also initialize LNEI99 from "logical_mod.f" (bmy, 11/5/04)
!  (2 ) Now also initialize LAVHRRLAI from "logical_mod.f" (bmy, 12/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now also initialize LMEGAN switch (tmf, bmy, 10/20/05)
!  (5 ) Now also initialize LEMEP, LGFED2BB switches and DATA_DIR_1x1
!        directory (bmy, 4/5/06)
!  (6 ) Now also intitialize LFUTURE (swu, bmy, 6/1/06)
!  (7 ) Now reference the EDGAR logical switches from "logical_mod.f"
!        (avd, bmy, 7/11/06)
!  (8 ) Now initialize the LVARTROP switch (phs, 9/14/06)
!  (9 ) Now initialize LOTDREG, LOTDLOC, LCTH, LMFLUX, LPRECON (bmy, 1/31/07)
!  (10) Now initialize LOTDSCALE (ltm, bmy, 9/24/07)
!  (11) Add MEGAN Monoterpenes switch (ccc, 2/2/09)
!  16 Oct 2009 - R. Yantosca - Now initialize LLINOZ
!  19 Nov 2009 - C. Carouge  - Initialize LMODISLAI and LPECCA
!  01 Dec 2009 - C. Carouge  - Initialize LNEI05 
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  07 Sep 2011 - P. Kasibhatla - Modified for GFED3
!  17 Jan 2012 - P. Kasibhatla - Modified for GFED3
!   8 Feb 2012 - R. Yantosca - Rewrote USE statements for clarity
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  26 Mar 2012 - R. Yantosca - Now turn off switch USE_OLSON_2001 by default
!                              (except for GEOS-5.7.2 met)
!  27 Mar 2012 - R. Yantosca - Cosmetic changes
!  13 Aug 2013 - M. Sulprizio- Add modifications for SOA + semivolatile POA
!                              simulation (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! INIT_INPUT begins here!
      !=================================================================

      ! Initialize directories
      DATA_DIR       = ''
      DATA_DIR_1x1   = ''
      GEOS_1_DIR     = ''
      GEOS_S_DIR     = ''
      GEOS_4_DIR     = ''
      TEMP_DIR       = ''
      RUN_DIR        = ''
      OH_DIR         = ''
      O3PL_DIR       = ''
      TPBC_DIR       = ''

      ! Initialize logicals
      LATEQ          = .FALSE.
      LCARB          = .FALSE.
      LDEAD          = .FALSE.
      LDUST          = .FALSE.
      LSULF          = .FALSE.
      LSOA           = .FALSE.
      LSSALT         = .FALSE.
      LCHEM          = .FALSE.
      LCONV          = .FALSE.
      LDBUG          = .FALSE.
      LDIAG          = .FALSE.
      LPRT           = .FALSE.
      LSTDRUN        = .FALSE.
      LDRYD          = .FALSE.
      LAEIC          = .FALSE.
      LANTHRO        = .FALSE.
      LBIONOX        = .FALSE.
      LBIOMASS       = .FALSE.
      LBIOFUEL       = .FALSE.
      LBIOGENIC      = .FALSE.
      LBBSEA         = .FALSE.
      LDYNOCEAN      = .FALSE.
      LEMEP          = .FALSE.
      LEMIS          = .FALSE.
      LEDGAR         = .FALSE.
      LEDGARNOx      = .FALSE. 
      LEDGARCO       = .FALSE. 
      LEDGARSHIP     = .FALSE. 
      LEDGARSOx      = .FALSE. 
      LFFNOX         = .FALSE.
      LFOSSIL        = .FALSE.
      LFUTURE        = .FALSE.
      LGFED2BB       = .FALSE.
      LGFED3BB       = .FALSE.
      LLIGHTNOX      = .FALSE.
      LMEGAN         = .FALSE.
      LMEGANMONO     = .FALSE.
      LMONOT         = .FALSE.
      LNEI99         = .FALSE.
      LNEI05         = .FALSE.
      LOTDLOC        = .FALSE.
      LSHIPSO2       = .FALSE.
      LSOILNOX       = .FALSE.
      LTOMSAI        = .FALSE.
      LWOODCO        = .FALSE.
      LFILL          = .FALSE.
      LTRAN          = .FALSE.
      LTPFV          = .FALSE.
      LWINDO         = .FALSE.
      LUNZIP         = .FALSE.
      LWAIT          = .FALSE.
      LTURB          = .FALSE.
      LSCHEM         = .FALSE.
      LSVGLB         = .FALSE.
      LSVCSPEC       = .FALSE.
      LSPLIT         = .FALSE.
      LWETD          = .FALSE.
      LVARTROP       = .FALSE.
      LLINOZ         = .FALSE.
      LPREINDHG      = .FALSE.
      USE_OLSON_2001 = .FALSE.
      LRETRO         = .FALSE.  !(wfr, 3/8/11)
      LRCP           = .FALSE.  !(cdh, 10/18/11)
      LRCPSHIP       = .FALSE.  !(cdh, 10/18/11)
      LRCPAIR        = .FALSE.  !(cdh, 10/18/11)
      LWARWICK_VSLS  = .FALSE.
      LSSABr2        = .FALSE.
      LFIX_PBL_BRO   = .FALSE.
      ! SOAupdate stuff for semivol POA (hotp 10/1/09)
      LSVPOA         = .FALSE.
      
      !Specifically for CO2 simulation (R Nassar, 2009-03-02)
      LGENFF         = .FALSE.
      LANNFF         = .FALSE.
      LMONFF         = .FALSE.
      LSEASBB        = .FALSE.
      LBIONETORIG    = .FALSE.
      LBIONETCLIM    = .FALSE.
      LBIODAILY      = .FALSE.
      LBIODIURNAL    = .FALSE.
      LOCN1997       = .FALSE.
      LOCN2009ANN    = .FALSE.
      LOCN2009MON    = .FALSE.
      LFFBKGRD       = .FALSE.
      LBIOSPHTAG     = .FALSE.
      LFOSSILTAG     = .FALSE.
      LSHIPEDG       = .FALSE.
      LSHIPICO       = .FALSE.
      LSHIPSCALE     = .FALSE.
      LSHIPTAG       = .FALSE.
      LPLANE         = .FALSE.
      LPLANESCALE    = .FALSE.
      LPLANETAG      = .FALSE.
      LCHEMCO2       = .FALSE.
                     
      ! Initialize counters
      CT1            = 0
      CT2            = 0
      CT3            = 0

#if   defined( GEOS_FP ) 
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% KLUDGE FOR TESTING GEOS-FP (bmy, 3/26/12)
      !%%% Use the Olson 2001 land map, which is at 0.25 x 0.25 resolution
      !%%% in order to facilitate the SEAC4RS simulations.  Eventually
      !%%% this will replace Olson 1992 for all other met field types,
      !%%% once more rigorous testing & validation has been done.
      USE_OLSON_2001 = .TRUE.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif

      END SUBROUTINE INIT_INPUT
!EOC
      END MODULE INPUT_MOD
