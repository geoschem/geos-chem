!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: input_mod.F
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
      USE CharPak_Mod,   ONLY : MAXDIM  => MAXSTRLEN
      USE inquireMod,    ONLY : findFreeLUN
      USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: Read_Input_File
      PUBLIC  :: Do_Error_Checks
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: SPLIT_ONE_LINE        
      PRIVATE :: READ_SIMULATION_MENU
      PRIVATE :: READ_ADVECTED_SPECIES_MENU
      PRIVATE :: READ_AEROSOL_MENU     
      PRIVATE :: READ_EMISSIONS_MENU
      PRIVATE :: READ_CHEMISTRY_MENU
      PRIVATE :: READ_RADIATION_MENU
      PRIVATE :: READ_TRANSPORT_MENU
      PRIVATE :: READ_CONVECTION_MENU
      PRIVATE :: READ_DEPOSITION_MENU
      PRIVATE :: READ_OUTPUT_MENU
      PRIVATE :: READ_DIAGNOSTIC_MENU
      PRIVATE :: READ_PLANEFLIGHT_MENU
      PRIVATE :: READ_OBSPACK_MENU
      PRIVATE :: SET_TINDEX
      PRIVATE :: READ_ND49_MENU      
      PRIVATE :: READ_ND50_MENU  
      PRIVATE :: READ_ND51_MENU  
      PRIVATE :: READ_ND51b_MENU  
      PRIVATE :: READ_ND63_MENU
      PRIVATE :: READ_PROD_LOSS_MENU 
      PRIVATE :: READ_NESTED_GRID_MENU
      PRIVATE :: READ_BENCHMARK_MENU  
      PRIVATE :: READ_CH4_MENU
      PRIVATE :: READ_PASSIVE_SPECIES_MENU
      PRIVATE :: VALIDATE_DIRECTORIES  
      PRIVATE :: CHECK_DIRECTORY
      PRIVATE :: CHECK_TIME_STEPS 
      PRIVATE :: IS_LAST_DAY_GOOD
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
!  23 Jun 2014 - R. Yantosca - Remove references to logical_mod.F
!  23 Jun 2014 - R. Yantosca - Removed INIT_INPUT routine
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  24 Nov 2014 - C. Keller   - Updates on timestep handling in ESMF environment
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  28 Jul 2016 - M. Sulprizio- Rename subroutine READ_TRACER_MENU to
!                              READ_ADVECTED_SPECIES_MENU
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and
!                              MERRA; Also remove LUNZIP, LWAIT, and TEMP_DIR
!                              options since they were only used for BPCH met
!  01 Nov 2017 - R. Yantosca - Now define MAXDIM from charpak_mod.F90
!  01 Nov 2017 - R. Yantosca - Now get ReadOneLine from charpak_mod.F90
!  08 Nov 2017 - R. Yantosca - Retire usage of ERROR_STOP; we need to send
!                              errors all the way up to the top level
!  26 Jan 2018 - M. Sulprizio- Moved Initialize_Geos_Grid and GC_Init_Extra
!                              routines to gc_environment_mod.F90
!  06 Feb 2018 - E. Lundgren - Change input timestep/freq units from min to sec
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!      LOGICAL            :: VERBOSE  = .FALSE.
      INTEGER            :: IU_GEOS
      INTEGER, PARAMETER :: FIRSTCOL = 26
      INTEGER            :: CT1      = 0
      CHARACTER(LEN=255) :: FILENAME = 'input.geos'

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
! In an ESMF environment, all time steps (chemistry, convection, emissions,
! dynamics) must be specified externally before calling this routine. This is
! done in routine GIGC\_Init\_Simulation (gigc\_initialization\_mod.F90).
! The time steps specified in input.geos are ignored.
!\\
!\\
! !INTERFACE:
! 
      SUBROUTINE READ_INPUT_FILE( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CHARPAK_MOD,        ONLY : ReadOneLine, STRREPL
      USE ErrCode_Mod
      USE FILE_MOD,           ONLY : IOERROR
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Now do not call INIT_INPUT, this is replaced
!                              by the INIT_GIGC_INPUT_OPT routine
!  15 Apr 2015 - R. Yantosca - Also define Input_Opt%ITS_A_NESTED_GRID here
!                              so that we can pass it to INITIALIZE_GEOS_GRID
!  04 Aug 2015 - M. Long     - Removed ".rc" file specifier. Not necessary.
!  16 Dec 2015 - R. Yantosca - Now pass State_Chm to READ_TRACER_MENU
!  04 Sep 2015 - C. Keller   - Added passive tracer menu
!  16 Jun 2016 - E. Lundgren - Move TOMAS init to GIGC_INIT_EXTRA
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL               :: EOF
      INTEGER               :: IOS
      CHARACTER(LEN=1)      :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)      :: SPACE = ' '
      CHARACTER(LEN=MAXDIM) :: LINE
      CHARACTER(LEN=255)    :: TOPTITLE, ErrMsg, ThisLoc

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

      ! Find a free file LUN
      IU_GEOS  = findFreeLUN()

      ! Assume success
      RC       = GC_SUCCESS

      ! For error handling
      ErrMsg   = ''
      ThisLoc  = 
     &   ' -> at Read_Input_File (in module GeosCore/input_mod.F)'

      ! Open file
      OPEN( IU_GEOS, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:1' )

      ! Read TOPTITLE for binary punch file
      TOPTITLE = ReadOneLine( IU_GEOS, EOF, IOS )
      IF ( EOF     ) RETURN
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:2' )

      ! Loop until EOF
      DO 
         
         ! Read a line from the file, exit if EOF
         LINE = ReadOneLine( IU_GEOS, EOF, IOS )
         IF ( EOF     ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS, 'read_input_file:3')
         
         ! Replace tab characters in LINE (if any) w/ spaces
         CALL STRREPL( LINE, TAB, SPACE )

         !=============================================================
         ! Call individual subroutines to read sections of the file
         ! 
         ! NOTE: You are pretty flexible in setting the order of the
         ! menus in the input file; however, a few guidelines apply:
         !
         ! (1) SIMULATION MENU should be listed first.
         ! (2) TIMESTEP MENU should be listed second.
         ! (3) ADVECTED SPECIES MENU should be listed third.
         ! (4) EMISSIONS, AEROSOL, CHEMISTRY, TRANSPORT, CONVECTION, 
         !      and DEPOSITION menus (in any order) should follow.
         ! (5) Diagnostic menus, including OUTPUT, DIAGNOSTIC,
         !      PLANEFLIGHT, ND48, ND49, ND50, ND51, and PROD_LOSS
         !      menus (in any order) should follow next.
         ! (6) The following menus have no other restriction and
         !      can be placed anywhere (but by convention we will
         !      place them after the diagnostic menu): NESTED GRID
         !      UNIX CMDS, ARCHIVED OH, and O3PL menus.
         !=============================================================
         IF      ( INDEX( LINE, 'SIMULATION MENU'  ) > 0 ) THEN
            CALL READ_SIMULATION_MENU( am_I_Root, Input_Opt, RC )             
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Simulation_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'TIMESTEP MENU' ) > 0 ) THEN
            CALL READ_TIMESTEP_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Timestep_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'ADVECTED SPECIES MENU' ) > 0 ) THEN
            CALL READ_ADVECTED_SPECIES_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Advected_Species_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'AEROSOL MENU'     ) > 0 ) THEN
            CALL READ_AEROSOL_MENU( am_I_Root, Input_Opt, RC )              
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Aerosol_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'EMISSIONS MENU'   ) > 0 ) THEN
            CALL READ_EMISSIONS_MENU( am_I_Root, Input_Opt, RC )             
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Emissions_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'CHEMISTRY MENU'   ) > 0 ) THEN
            CALL READ_CHEMISTRY_MENU( am_I_Root, Input_Opt, RC )              
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Chemistry_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                                 
         ELSE IF ( INDEX( LINE, 'TRANSPORT MENU'   ) > 0 ) THEN
            CALL READ_TRANSPORT_MENU( am_I_Root, Input_Opt, RC )              
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Transport_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                             
         ELSE IF ( INDEX( LINE, 'CONVECTION MENU'  ) > 0 ) THEN
            CALL READ_CONVECTION_MENU( am_I_Root, Input_Opt, RC )            
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Convection_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                             
         ELSE IF ( INDEX( LINE, 'DEPOSITION MENU'  ) > 0 ) THEN
            CALL READ_DEPOSITION_MENU( am_I_Root, Input_Opt, RC )             
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Depositon_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'RADIATION MENU'   ) > 0 ) THEN
            CALL READ_RADIATION_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Radiation_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'CO SIM MENU'     ) > 0 ) THEN
            CALL READ_CO_SIM_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_CO_Sim_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                    
         ELSE IF ( INDEX( LINE, 'CO2 SIM MENU'     ) > 0 ) THEN
            CALL READ_CO2_SIM_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_CO2_Sim_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                    
         ELSE IF ( INDEX( LINE, 'GAMAP MENU'       ) > 0 ) THEN
            CALL READ_GAMAP_MENU( am_I_Root, Input_Opt, RC )                 
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Gamap_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                                  
         ELSE IF ( INDEX( LINE, 'OUTPUT MENU'      ) > 0 ) THEN
            CALL READ_OUTPUT_MENU( am_I_Root, Input_Opt, RC )                 
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Output_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ! NOTE: This will get always get read, but the various bpch
         ! diagnostics will only be initialized if BPCH_DIAG=y (bmy, 1/16/18)
         ELSE IF ( INDEX( LINE, 'DIAGNOSTIC MENU'  ) > 0 ) THEN
            CALL READ_DIAGNOSTIC_MENU( am_I_Root, Input_Opt, RC )             
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Diagnostic_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'PLANEFLIGHT MENU' ) > 0 ) THEN
            CALL READ_PLANEFLIGHT_MENU( am_I_Root, Input_Opt, RC )             
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Planeflight_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'OBSPACK MENU' ) > 0 ) THEN
            CALL READ_OBSPACK_MENU( am_I_Root, Input_Opt, RC )   
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ObsPack_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ! Read even if BPCH_DIAG=n so that time series diags may be
         ! output if other diagnostics are output as netcdf (ewl, 4/5/18)
         ELSE IF ( INDEX( LINE, 'ND48 MENU'        ) > 0 ) THEN
            CALL READ_ND48_MENU( am_I_Root, Input_Opt, RC )                  
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND48_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'ND49 MENU'        ) > 0 ) THEN
            CALL READ_ND49_MENU( am_I_Root, Input_Opt, RC )                   
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND49_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                                  
         ELSE IF ( INDEX( LINE, 'ND50 MENU'        ) > 0 ) THEN
            CALL READ_ND50_MENU( am_I_Root, Input_Opt, RC )                   
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND50_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                                  
         ELSE IF ( INDEX( LINE, 'ND51 MENU'        ) > 0 ) THEN
            CALL READ_ND51_MENU( am_I_Root, Input_Opt, RC )                   
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND51_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                                  
         ELSE IF ( INDEX( LINE, 'ND51b MENU'       ) > 0 ) THEN
            CALL READ_ND51b_MENU( am_I_Root, Input_Opt, RC )  
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND51b_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'ND63 MENU'        ) > 0 ) THEN
            CALL READ_ND63_MENU( am_I_Root, Input_Opt, RC )                   
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_ND63_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

#if defined( BPCH_DIAG )
         ELSE IF ( INDEX( LINE, 'PROD & LOSS MENU' ) > 0 ) THEN
            CALL READ_PROD_LOSS_MENU( am_I_Root, Input_Opt, RC )               
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Prod_Loss_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
#endif

         ELSE IF ( INDEX( LINE, 'NESTED GRID MENU' ) > 0 ) THEN 
            CALL READ_NESTED_GRID_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Nested_Grid_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'BENCHMARK MENU'   ) > 0 ) THEN 
            CALL READ_BENCHMARK_MENU( am_I_Root, Input_Opt, RC  )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error encountered in ""!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
         
         ELSE IF ( INDEX( LINE, 'POPS MENU'        ) > 0 ) THEN
            CALL READ_POPS_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_POPS_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
           
         ELSE IF ( INDEX( LINE, 'MERCURY MENU'     ) > 0 ) THEN 
            CALL READ_MERCURY_MENU( am_I_Root, Input_Opt, RC  )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Mercury_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF
                                   
         ELSE IF ( INDEX( LINE, 'CH4 MENU'         ) > 0 ) THEN 
            CALL READ_CH4_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_CH4_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ELSE IF ( INDEX( LINE, 'PASSIVE SPECIES' ) > 0 ) THEN 
            CALL READ_PASSIVE_SPECIES_MENU( am_I_Root, Input_Opt, RC )
            IF ( RC /= GC_SUCCESS ) THEN
               ErrMsg = 'Error in "Read_Passive_Species_Menu"!'
               CALL GC_Error( ErrMsg, RC, ThisLoc )
               RETURN
            ENDIF

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
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! Do not initialize diagnostics when we are connecting to
      ! the GEOS-5 GCM via ESMF. (bmy, 11/2/12)
      !-----------------------------------------------------------------
#else

      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! We only need to initialize GEOS-5 diagnostics if we are not
      ! using ESMF.  The GEOS-5 GCM has its own History component
      ! for diagnostic output (bmy, 3/4/13)
      !-----------------------------------------------------------------

      ! Make sure all directories are valid
      CALL VALIDATE_DIRECTORIES( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error in "Validate_Directories"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      ! Check GEOS-CHEM timesteps
      CALL CHECK_TIME_STEPS( am_I_Root, Input_Opt, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error in "Check_Time_Steps"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      END SUBROUTINE READ_INPUT_FILE
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
      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, 
     &                           N_EXP,   LOCATION,  RC ) 
!
! !USES:
!
      USE Charpak_Mod, ONLY: ReadOneLine, StrSplit
      Use ErrCode_Mod
      USE File_Mod,    ONLY: IoError
!
! !INPUT PARAMETERS: 
!
      INTEGER,            INTENT(IN)  :: N_EXP      ! Expected # of substrs
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION   ! Name of calling routine
!
! !OUTPUT PARAMETERS:
!
      ! Array of substrings (separated by " ")
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM) ! Substrings
      INTEGER,            INTENT(OUT) :: N_SUBSTRS       ! # of substrings
      INTEGER,            INTENT(OUT) :: RC              ! Success/failure?  
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  17 Sep 2013 - R. Yantosca - Extend LINE to 500 chars to allow more tracers
!   8 Nov 2017 - R. Yantosca - Now return error conditon to calling routine 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER               :: IOS
      LOGICAL               :: EOF

      ! Strings
      CHARACTER(LEN=10)     :: IOS_Str
      CHARACTER(LEN=255)    :: ErrMsg, ThisLoc
      CHARACTER(LEN=MAXDIM) :: LINE

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================      

      ! Initialize
      RC      = GC_SUCCESS
      IOS     = 0
      EOF     = .FALSE.
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = ' -> at Split_One_Line (in module GeosCore/input_mod.F)'
      
      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = ReadOneLine( IU_GEOS, EOF, IOS )
      
      ! Trap premature end-of-file error
      IF ( EOF ) THEN
         ErrMsg = 'Unexpected end-of-file reading input.geos!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      
      ! Trap other I/O error conditions
      IF ( IOS > 0 ) THEN
         WRITE( IOS_Str, '(i10)' ) IOS
         ErrMsg = 'I/O error number: ' // TRIM( IOS_STR ) // 
     &            'encountered when readiang "input.geos"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong # of substrings
      IF ( N_EXP /= N_SUBSTRS ) THEN
         ErrMsg= 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )
         WRITE( 6, '(a)' ) TRIM( ErrMSg )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE TIME_MOD
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
!  19 May 2014 - C. Keller   - Now read HEMCO configuration file.
!  23 Jun 2014 - R. Yantosca - Remove references to logical_mod.F
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2
!  12 Aug 2015 - R. Yantosca - Add support for 05x0625 grids
!  13 Aug 2015 - R. Yantosca - Bug fix: prefix RES_DIR in front of MERRA2_DIR
!  27 Jul 2016 - M. Sulprizio- Remove LSVGLB and OUT_RST_FILE options. Restart
!                              files are now always saved out and the output
!                              restart file name is hardcoded in restart_mod.F.
!  09 Aug 2016 - E. Lundgren - Remove call to routine set_restart; use input
!                              rst filename directly in restart_mod.F.
!  07 Nov 2017 - R. Yantosca - Don't assign to Input_Opt%VARTROP, and will
!                              eventually remove the line from input.geos
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  26 Jan 2018 - M. Sulprizio- Remove obsolete options LVARTROP and reading of
!                              DATA_DIR_1x1
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, C
      REAL(fp)           :: JulianDateStart,  JulianDateEnd
#if defined( ESMF_ ) || defined( MODEL_ )
      INTEGER            :: H,       M,       S
      REAL(f4)           :: init_UTC
#endif

      ! Strings
      CHARACTER(LEN=6)   :: TimeStr
      CHARACTER(LEN=8)   :: DateStr
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_SIMULATION_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = ' -> at Read_Simulation_Menu (in GeosCore/input_mod.F)'

      !-----------------------------------------------------------------
      ! Simulation start and end time
      !-----------------------------------------------------------------

      ! Start YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'NYMDb, NHMSb', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NYMDb, Input_Opt%NHMSb

      ! Make sure the starting date NYMDb is valid
      IF ( .not. Valid_Date( Input_Opt%NYMDb ) ) THEN
         WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDb
         ErrMsg = 'Input%Opt%NYMDb = ' // DateStr        // 
     &            ' is not a valid calendar date!'       //
     &            ' Please check your "input.geos" file.'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Make sure the starting time NHMSb is valid
      IF ( .not. Valid_Time( Input_Opt%NHMSb ) ) THEN
         WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSb
         ErrMsg = 'Input%Opt%NHMSb = ' // TimeStr        // 
     &            ' is not a valid clock time!'          //
     &            ' Please check your "input.geos" file.'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! End YYYYMMDD, HHMMSS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'NYMDe, NHMSe', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NYMDe, Input_Opt%NHMSe

      ! Make sure the starting date NYMDb is valid
      IF ( .not. Valid_Date( Input_Opt%NYMDe ) ) THEN
         WRITE( DateStr, '(i8.8)' ) Input_Opt%NYMDe
         ErrMsg = 'Input%Opt%NYMDe = ' // DateStr        // 
     &            ' is not a valid calendar date!'       //
     &            ' Please check your "input.geos" file.'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Make sure the ending time NHMSe is valid
      IF ( .not. Valid_Time( Input_Opt%NHMSe ) ) THEN
         WRITE( TimeStr, '(i6.6)' ) Input_Opt%NHMSe
         ErrMsg = 'Input%Opt%NHMSe = ' // TimeStr        // 
     &            ' is not a valid clock time!'          //
     &            ' Please check your "input.geos" file.'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Compute the length of the simulation, in elapsed seconds
      JulianDateStart = GET_JD( Input_Opt%NymdB, Input_Opt%NhmsB )
      JulianDateEnd   = GET_JD( Input_Opt%NymdE, Input_Opt%NhmsE )
      Input_Opt%SimLengthSec = NINT( ( JulianDateEnd - JulianDateStart )
     &                         * 86400_f8)

      !-----------------------------------------------------------------
      ! Data directories
      !-----------------------------------------------------------------

      ! Run directory
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'RUN_DIR', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%RUN_DIR

      ! Root data dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DATA_DIR', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%DATA_DIR

      ! Make sure DATA-DIR ends with a "/" character
      C = LEN_TRIM( Input_Opt%DATA_DIR )
      IF ( Input_Opt%DATA_DIR(C:C) /= '/' ) THEN
         Input_Opt%DATA_DIR = TRIM( Input_Opt%DATA_DIR ) // '/'
      ENDIF

      ! Create CHEM_INPUTS directory
      Input_Opt%CHEM_INPUTS_DIR = TRIM( Input_Opt%DATA_DIR ) // 
     &                            'CHEM_INPUTS/'

      ! Define resolution-specific directory name RES_DIR
#if   defined( GRID4x5 )
      Input_Opt%RES_DIR = 'GEOS_4x5/'
#elif defined( GRID2x25 ) 
      Input_Opt%RES_DIR = 'GEOS_2x2.5/'
#elif defined( GRID05x0625 ) && defined( NESTED_AS )
      Input_Opt%RES_DIR = 'GEOS_0.5x0.625_AS/'
#elif defined( GRID05x0625 ) && defined( NESTED_EU )
      Input_Opt%RES_DIR = 'GEOS_0.5x0.625_EU/'
#elif defined( GRID05x0625 ) && defined( NESTED_NA )
      Input_Opt%RES_DIR = 'GEOS_0.5x0.625_NA/'
#elif defined( GRID025x03125 ) && defined( NESTED_CH ) 
      Input_Opt%RES_DIR = 'GEOS_0.25x0.3125_CH/'
#elif defined( GRID025x03125 ) && defined( NESTED_EU ) 
      Input_Opt%RES_DIR = 'GEOS_0.25x0.3125_EU/'
#elif defined( GRID025x03125 ) && defined( NESTED_NA ) 
      Input_Opt%RES_DIR = 'GEOS_0.25x0.3125_NA/'
#endif

      ! Determine if we are running a nested-grid simulation based on
      ! the settings of the C-preprocessor switches.
#if defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_NA ) || defined( NESTED_AS )
      Input_Opt%ITS_A_NESTED_GRID = .TRUE.
#endif

      ! I0, J0 for nested-grid simulations
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'Nested I0, J0', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NESTED_I0, 
     &                        Input_Opt%NESTED_J0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'SIMULATION MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, 100 ) 'Start time of run           : ', 
     &                    Input_Opt%NYMDb, Input_Opt%NHMSb
         WRITE( 6, 100 ) 'End time of run             : ',
     &                    Input_Opt%NYMDe, Input_Opt%NHMSe
         WRITE( 6, 110 ) 'Run directory               : ',
     &                    TRIM( Input_Opt%RUN_DIR )
         WRITE( 6, 110 ) 'Data Directory              : ',
     &                    TRIM( Input_Opt%DATA_DIR )
         WRITE( 6, 110 ) 'CHEM_INPUTS directory       : ',
     &                    TRIM( Input_Opt%CHEM_INPUTS_DIR )
         WRITE( 6, 110 ) 'Resolution-specific dir     : ',
     &                    TRIM( Input_Opt%RES_DIR )
         WRITE( 6, 120 ) 'Is this a nested-grid sim?  : ', 
     &                    Input_Opt%ITS_A_NESTED_GRID
         WRITE( 6, 130 ) 'Global offsets I0, J0       : ', 
     &                    Input_Opt%NESTED_I0, 
     &                    Input_Opt%NESTED_J0
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
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
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
      !-----------------------------------------------------------------

      ! Get the starting UTC time from Input_Opt%NHMSb for use below
      CALL YMD_Extract( Input_Opt%NHMSb, H, M, S )
      init_UTC = ( H + ( M / 60 ) + ( S / 3600 ) )

      ! Pass the values for the start & end times of the simulation directly
      ! to GeosUtil/time_mod.F via subroutine ACCEPT_EXTERNAL_DATE_TIME.
      ! (bmy, 12/6/12)
      CALL Accept_External_Date_Time( am_I_Root   = am_I_Root, 
     &                                value_NYMDb = Input_Opt%NYMDb,
     &                                value_NHMSb = Input_Opt%NHMSb,
     &                                value_NYMDe = Input_Opt%NYMDe,
     &                                value_NHMSe = Input_Opt%NHMSe,
     &                                value_NYMD  = Input_Opt%NYMDb,
     &                                value_NHMS  = Input_Opt%NHMSb,
     &                                value_UTC   = init_UTC,
     &                                RC          = RC               )
#else
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! If we are not using ESMF, then call the traditional GEOS-Chem
      ! timing routines (from GeosUtil/time_mod.F) to set the start &
      ! end times of the simulation, as well as the current time.
      ! (bmy, 12/6/12)
      !-----------------------------------------------------------------

      ! Set start time of run in "time_mod.f"
      CALL SET_BEGIN_TIME( Input_Opt%NYMDb, Input_Opt%NHMSb, RC  )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Set_Begin_Time"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Set end time of run in "time_mod.f"
      ErrMsg = 'Error encountered in "Set_Begin_Time"!'
      CALL SET_END_TIME( Input_Opt%NYMDe, Input_Opt%NHMSe, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Set the current time
      CALL SET_CURRENT_TIME()
#endif

#if defined( BPCH_DIAG )
      ! Set the start of the 1st diagnostic interval
      CALL SET_DIAGb( GET_TAU() )
#endif     

      ! Set menu counter
      CT1 = CT1 + 1

      END SUBROUTINE READ_SIMULATION_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_timestep_menu
!
! !DESCRIPTION: Subroutine READ\_TIMESTEP\_MENU reads the TIMESTEP MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TIMESTEP_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  09 Aug 2017 - M. Sulprizio- Initial version
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
     
      !=================================================================
      ! READ_TIMESTEP_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Timestep_Menu (in module GeosCore/input_mod.F)'

      ! Transport/convection timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_DYN', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%TS_DYN
      Input_Opt%TS_CONV = Input_Opt%TS_DYN

      ! Chemistry/emissions timestep
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_CHEM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%TS_CHEM
      Input_Opt%TS_EMIS = Input_Opt%TS_CHEM

      ! Separator line                          
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF


      !=================================================================
      ! Error checks
      !=================================================================

      IF ( Input_Opt%SimLengthSec < Input_Opt%TS_DYN .or.
     &     Input_Opt%SimLengthSec < Input_Opt%TS_CHEM ) THEN
         IF ( am_I_Root ) THEN

            WRITE( 6, *   ) '' 
            WRITE( 6, *   ) 'The length of the simulation is shorter '
            WRITE( 6, *   ) 'than the transport and/or chemistry '
            WRITE( 6, *   ) 'timesteps. Check the settings in '
            WRITE( 6, *   ) 'the "input.geos" file.'
            WRITE( 6, *   ) '' 
            WRITE( 6, 100 ) 'Transport/Convection [sec]  : ', 
     &                       Input_Opt%TS_DYN
            WRITE( 6, 100 ) 'Chemistry/Emissions  [sec]  : ', 
     &                       Input_Opt%TS_CHEM
            WRITE( 6, 100 ) 'Simulation duration  [sec]  : ',
     &                       Input_Opt%SimLengthSec
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN

         ENDIF
      ENDIF

#if   defined( MERRA2 ) && defined( GRID05x0625 )
      IF ( Input_Opt%ITS_A_CH4_SIM .or. Input_Opt%ITS_A_CO2_SIM ) THEN
         IF ( Input_Opt%TS_DYN > 300 .or. Input_Opt%TS_CHEM > 600 ) THEN
            IF ( am_I_Root ) THEN

            WRITE( 6,* ) '' 
            WRITE( 6,* ) 'It has been noted that MERRA-2 nested grid'
            WRITE( 6,* ) ' simulations can have very high species'
            WRITE( 6,* ) ' concentrations in the stratosphere caused'
            WRITE( 6,* ) ' by a violation of the CFL condition due to'
            WRITE( 6,* ) ' strong stratospheric winds. This is'
            WRITE( 6,* ) ' especially problematic when using total'
            WRITE( 6,* ) ' column concentrations. To avoid the issue,'
            WRITE( 6,* ) ' a timestep of 5/10 instead of 10/20 is'
            WRITE( 6,* ) ' recommended for CH4 and CO2 simulations.'
            WRITE( 6,* ) ''
            WRITE( 6,* ) 'You may remove this trap at your own peril,'
            WRITE( 6,* ) ' by commenting out the call to GC_ERROR in'
            WRITE( 6,* ) ' GeosCore/input_mod.F. '
            WRITE( 6,* ) ''
            WRITE( 6,* ) 'See the MERRA-2 implementation details page'
            WRITE( 6,* ) ' on the GEOS-Chem wiki for details'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN

            ENDIF
         ENDIF
      ENDIF
#endif

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'TIMESTEP MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, 100     ) 'Transport/Convection [sec]  : ', 
     &                        Input_Opt%TS_DYN
         WRITE( 6, 100     ) 'Chemistry/Emissions  [sec]  : ', 
     &                        Input_Opt%TS_CHEM
      ENDIF

      ! Format statements
 100  FORMAT( A, I5  )

      END SUBROUTINE READ_TIMESTEP_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_advected_species_menu
!
! !DESCRIPTION: Subroutine READ\_ADVECTED\_SPECIES\_MENU reads the ADVECTED
!  SPECIES MENU section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ADVECTED_SPECIES_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput, Set_Input_Opt_Advect
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
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  21 Aug 2014 - R. Yantosca - Bug fix: add ITS_A_RnPbBe_SIM to the test
!                              that defines the ITS_A_SPECIALTY_SIM flag
!  24 Mar 2015 - E. Lundgren - Move init_tracer to within APM ifdef block 
!  23 Jun 2016 - R. Yantosca - Remove call to TRACERID, it's obsolete
!  28 Jul 2016 - M. Sulprizio- Rename TRACER MENU to ADVECTED SPECIES MENU;
!                              Remove ID, MW, and member columns from menu
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  26 Jan 2018 - M. Sulprizio- Remove line for setting N_ADVECT and instead
!                              count the number of advected species listed in
!                              this menu and set N_ADVECT accordingly
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, T

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ADVECTED_SPECIES_MENU begins here!
      !
      ! Get the simulation type
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Advected_Species_Menu (in GeosCore/input_mod.F)'

      ! Simulation type
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'SIM_TYPE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%SIM_TYPE

      ! Set fields of Input Options object
      Input_Opt%ITS_A_RnPbBe_SIM       = ( Input_Opt%SIM_TYPE == 1  )
      Input_Opt%ITS_A_CH3I_SIM         = ( Input_Opt%SIM_TYPE == 2  )
      Input_Opt%ITS_A_FULLCHEM_SIM     = ( Input_Opt%SIM_TYPE == 3  )
      Input_Opt%ITS_A_HCN_SIM          = ( Input_Opt%SIM_TYPE == 4  )
      Input_Opt%ITS_A_TAGO3_SIM        = ( Input_Opt%SIM_TYPE == 6  )
      Input_Opt%ITS_A_TAGCO_SIM        = ( Input_Opt%SIM_TYPE == 7  )
      Input_Opt%ITS_A_C2H6_SIM         = ( Input_Opt%SIM_TYPE == 8  )
      Input_Opt%ITS_A_CH4_SIM          = ( Input_Opt%SIM_TYPE == 9  )
      Input_Opt%ITS_AN_AEROSOL_SIM     = ( Input_Opt%SIM_TYPE == 10 )
      Input_Opt%ITS_A_MERCURY_SIM      = ( Input_Opt%SIM_TYPE == 11 )
      Input_Opt%ITS_A_CO2_SIM          = ( Input_Opt%SIM_TYPE == 12 )
      Input_Opt%ITS_A_H2HD_SIM         = ( Input_Opt%SIM_TYPE == 13 )
      Input_Opt%ITS_A_POPS_SIM         = ( Input_Opt%SIM_TYPE == 14 )

      ! Also pick a name corresponding to each simulation type
      ! This replaces function GET_SIM_NAME from tracer_mod.F
      SELECT CASE( Input_Opt%SIM_TYPE )
         CASE( 1 ) 
            Input_Opt%SIM_NAME = 'Rn-Pb-Be'
         CASE( 2 ) 
            Input_Opt%SIM_NAME = 'CH3I'
         CASE( 3 ) 
            Input_Opt%SIM_NAME = 'NOx-Ox-Hydrocarbon-Aerosol'
         CASE( 4 )
            Input_Opt%SIM_NAME = 'HCN'
         CASE( 5 )
            Input_Opt%SIM_NAME = ''
         CASE( 6 )
            Input_Opt%SIM_NAME = 'Tagged O3'
         CASE( 7 )
            Input_Opt%SIM_NAME = 'Tagged CO'
         CASE( 8 ) 
            Input_Opt%SIM_NAME = 'Tagged C2H6'
         CASE( 9 )
            Input_Opt%SIM_NAME = 'CH4'
         CASE( 10 ) 
            Input_Opt%SIM_NAME = 'Offline Aerosol'
         CASE( 11 ) 
            Input_Opt%SIM_NAME = 'Mercury'
         CASE( 12 )
            Input_Opt%SIM_NAME = 'CO2'
         CASE( 13 )
            Input_Opt%SIM_NAME = 'H2 and HD'
         CASE( 14 )
            Input_Opt%SIM_NAME = 'POPs'
         CASE DEFAULT
            Input_Opt%SIM_NAME = 'UNKNOWN'
       END SELECT

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Read advected species name
      !=================================================================

      ! Initialize counter
      T = 0

      ! Do until exit
      DO

         ! Read species name
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'species', RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Exit here if separator line is encountered
         IF ( INDEX ( TRIM(SUBSTRS(1)), '-----' ) > 0 ) EXIT

         ! Increase number of advected species by one
         T = T + 1

         ! Stop simulation and print warning if we exceed maximum number
         ! of advected species hardcoded in input_opt_mod.F90
         IF ( T > Input_Opt%Max_AdvectSpc ) THEN
            ErrMsg = 'Number of advected species exceeds maximum. ' //
     &             'This value can be modified in input_opt_mod.F90.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Save advected species name
         Input_Opt%AdvectSpc_Name(T)  = TRIM( SUBSTRS(1) )

      ENDDO

      ! Total number of advected species
      Input_Opt%N_ADVECT = T

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Split into tagged species (turn off for full-chemistry)
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

         ! There are no tagged species for fullchem
         Input_Opt%LSPLIT = .FALSE.

      ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

         ! Need Hg0, Hg2, HgP for tagged Mercury
         Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 3 )
         
      ELSE IF ( Input_Opt%ITS_A_HCN_SIM ) THEN

         ! Need HCN, CH3CN for HCN simulation
         Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 2 )

      ELSE

         ! All other simulations: tagged if more than 1 species listed
         Input_Opt%LSPLIT = ( Input_Opt%N_ADVECT > 1 )

      ENDIF

      ! Initialize arrays in Input_Opt that depend on N_ADVECT
      CALL Set_Input_Opt_Advect( am_I_Root, Input_Opt, RC )

      ! Set menu counter
      CT1 = CT1 + 1

      END SUBROUTINE READ_ADVECTED_SPECIES_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  Move error checks that depend on species indices to the subroutine
!  DO_ERROR_CHECKS.  This is now called from GC_INIT_EXTRA, after the
!  initialization of the species database.
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
!  11 Apr 2013 - S.D. Eastham- Added gravitational settling flag
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA +
!                              semivolatile POA simulations (H. Pye)
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  16 Apr 2014 - M. Sulprizio- Now read path for PSC restart file
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  09 Apr 2015 - M. Sulprizio- Removed options for DEAD dust emissions, NAP
!                              emissions, POA emissions scale, and the PSC
!                              restart file. These options are now handled
!                              by HEMCO.
!  08 Jul 2015 - E. Lundgren - Added LMPOA option for marine organic aerosols
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function
!  22 Jun 2016 - R. Yantosca - Move error checks to DO_ERROR_CHECKS routine
!  16 Jan 2018 - M. Sulprizio- Add option for using spatially and seasonally
!                              varying OM/OC from Philip et al. (2014)
!  26 Jan 2018 - M. Sulprizio- Remove obsolete options for LCRYST and LDICARB
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, T

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_AEROSOL_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Aerosol_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &             'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Use online sulfate aerosols? 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSULF', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSULF

      ! Use metal catalyzed oxidation of SO2?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LMETALCATSO2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LMETALCATSO2

      ! Use online carbon aerosols? 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCARB', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCARB
      
      ! Use brown carbon aerosols? 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBRC', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBRC

      ! Use secondary organic aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSOA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSOA

      ! SOAupdate: Add Semi-volatile POA switch (hotp 8/9/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSVPOA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSVPOA

      ! Use spatially and seasonally varying OM/OC?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LOMOC', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LOMOC

      ! Use online dust aerosols ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDUST', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LDUST

      ! Use SO2 and HNO3 uptake on dust aerosols 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDSTUP', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LDSTUP

      ! Use online sea-salt aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSSALT', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSSALT      

      ! Accum mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'SALA_REDGE_um', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO T = 1, N
         READ( SUBSTRS(T), * ) Input_Opt%SALA_REDGE_um(T)
      ENDDO

      ! Coarse mode seasalt radii bin edges [um]
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'SALC_REDGE_um', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO T = 1, N
         READ( SUBSTRS(T), * ) Input_Opt%SALC_REDGE_um(T)
      ENDDO

      ! Use marine organic aerosols?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LMPOA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LMPOA

      ! Apply gravitational settling in stratosphere?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LGRAVSTRAT', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LGRAVSTRAT  

      ! Use solid polar stratospheric clouds (PSCs)?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSOLIDPSC', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSOLIDPSC   

      ! Allow homogeneous nucleation of NAT?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LHOMNUCNAT', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LHOMNUCNAT 

      ! NAT supercooling requirement (K)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'T_NAT_SUPERCOOL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%T_NAT_SUPERCOOL

      ! Ice supersaturation ratio requirement
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'P_ICE_SUPERSAT', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%P_ICE_SUPERSAT

      ! Perform PSC-related heterogeneous chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPSCCHEM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LPSCCHEM 

      ! Include stratospheric aerosols optical depths?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSTRATOD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSTRATOD

      ! Include BC absorption enhancement due to coating? (xnw, 8/24/15)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBCAE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBCAE

      ! Define BC absorption enhancement (xnw, 8/24/15)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'BCAE_1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%BCAE_1

      ! Define BC absorption enhancement (xnw, 8/24/15)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'BCAE_2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%BCAE_2

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error checks
      !=================================================================

      ! Make sure that SALA, SALC bins are contiguous
      IF ( Input_Opt%SALA_REDGE_um(2) /= 
     &     Input_Opt%SALC_REDGE_um(1)     ) THEN
         ErrMsg = 'SALA and SALC bin edges are not contiguous!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN         
      ENDIF

      ! Turn off switches for simulations that don't use aerosols
      IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. 
     &     ( .not. Input_OPt%ITS_AN_AEROSOL_SIM ) ) THEN
         Input_Opt%LSULF        = .FALSE.
         Input_Opt%LMETALCATSO2 = .FALSE.
         Input_Opt%LCARB        = .FALSE.
         Input_Opt%LBRC         = .FALSE.
         Input_Opt%LSOA         = .FALSE.
         Input_Opt%LDUST        = .FALSE.
         Input_Opt%LSSALT       = .FALSE.
         Input_Opt%LMPOA        = .FALSE.
         Input_Opt%LSVPOA       = .FALSE.
         Input_Opt%LOMOC        = .FALSE.
         Input_Opt%LBCAE        = .FALSE.
      ENDIF

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'AEROSOL MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 100     ) 'Online SULFATE AEROSOLS?    : ', 
     &                        Input_Opt%LSULF
         WRITE( 6, 100     ) 'Metal catalyzed SO2 ox.?    : ',
     &                        Input_Opt%LMETALCATSO2
         WRITE( 6, 100     ) 'Online CARBON AEROSOLS?     : ',
     &                        Input_Opt%LCARB
         WRITE( 6, 100     ) 'Brown Carbon Aerosol?       : ',
     &                        Input_Opt%LBRC
         WRITE( 6, 100     ) 'Online COMPLEX SOA?         : ',
     &                        Input_Opt%LSOA
         WRITE( 6, 100     ) 'Semivolatile POA?           : ',
     &                        Input_Opt%LSVPOA
         WRITE( 6, 100     ) 'Spatial & seasonal OM/OC?   : ',
     &                        Input_Opt%LOMOC
         WRITE( 6, 100     ) 'Online DUST AEROSOLS?       : ',
     &                        Input_Opt%LDUST
         WRITE( 6, 100     ) 'Acid uptake on dust?        : ',
     &                        Input_Opt%LDSTUP
         WRITE( 6, 100     ) 'Online SEA SALT AEROSOLS?   : ',
     &                        Input_Opt%LSSALT
         WRITE( 6, 110     ) 'Accum  SEA SALT radii [um]  : ', 
     &                        Input_Opt%SALA_REDGE_um(1), 
     &                        Input_Opt%SALA_REDGE_um(2)
         WRITE( 6, 110     ) 'Coarse SEA SALT radii [um]  : ',
     &                        Input_Opt%SALC_REDGE_um(1), 
     &                        Input_Opt%SALC_REDGE_um(2)
         WRITE( 6, 100     ) 'MARINE ORGANIC AEROSOLS?    : ',
     &                        Input_Opt%LMPOA
         WRITE( 6, 100     ) 'Settle strat. aerosols?     : ',
     &                        Input_Opt%LGRAVSTRAT
         WRITE( 6, 100     ) 'Online SOLID PSC aerosols?  : ',
     &                        Input_Opt%LSOLIDPSC
         WRITE( 6, 100     ) 'Allow hom. NAT nucleation?  : ',
     &                        Input_Opt%LHOMNUCNAT
         WRITE( 6, 120     ) 'NAT supercooling requirement: ',
     &                        Input_Opt%T_NAT_SUPERCOOL
         WRITE( 6, 120     ) 'Ice supersaturation req.    : ',
     &                        ((Input_Opt%P_ICE_SUPERSAT-1)*1.e+2_fp)
         WRITE( 6, 100     ) 'Perform PSC het. chemistry? : ',
     &                        Input_Opt%LPSCCHEM
         WRITE( 6, 100     ) 'Use strat. aerosol OD?      : ',
     &                        Input_Opt%LSTRATOD
         WRITE( 6, 100     ) 'BC Absorption Enhancement?  : ',
     &                        Input_Opt%LBCAE
         !for BC absorption enhancement, (xnw, 8/24/15)
         WRITE( 6, 105     ) 'Hydrophilic BC AE factor    : ',
     &                        Input_Opt%BCAE_1
         WRITE( 6, 105     ) 'Hydrophobic BC AE factor    : ',
     &                        Input_Opt%BCAE_2
      ENDIF

 100  FORMAT( A, L5     )
 105  FORMAT( A, f6.2 )
 110  FORMAT( A, f6.2, ' - ', f6.2 )
 120  FORMAT( A, f6.2, 'K' )

      END SUBROUTINE READ_AEROSOL_MENU
!EOC
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,      ONLY : GET_YEAR
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
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  Therefore now test for Ind_() > 0  for a
!  valid species.
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
!  03 Jun 2014 - R. Yantosca - Radically simplify this routine since most
!                              emissions options are now specified in HEMCO
!  23 Jun 2014 - R. Yantosca - Remove reference to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Move call to INIT_MODIS_LAIT to GIGC_INIT_EXTRA
!  04 Sep 2014 - R. Yantosca - Bug fix: Now use Input_Opt%OTDLOC since we have
!                              now removed logical_mod.F
!  24 Jun 2015 - R. Yantosca - Now always turn off emissions for mass cons test
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function 
!  22 Jun 2016 - R. Yantosca - Move some error checks to DO_ERROR_CHECKS
!  20 Sep 2016 - R. Yantosca - Use "I4" format to read in CFCYEAR
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_EMISSIONS_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Emissions_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLEMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LEMIS

      ! HEMCO Input file 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'HcoConfigFile', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%HcoConfigFile

      ! Set PBL BrO concentrations equal to 1ppt during the day?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFIX_PBL_BRO', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LFIX_PBL_BRO

      ! Separator line (start of UCX options)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Use variable methane emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCH4EMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCH4EMIS

      ! Separator line (start of surface BC options)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Use fixed methane surface mixing ratio?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCH4SBC', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCH4SBC

      ! Use fixed OCS surface mixing ratio?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LOCSEMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LOCSEMIS

      ! Use Montreal CFC emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCFCEMIS', RC )
      IF  ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCFCEMIS

      ! Use chlorine emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCLEMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCLEMIS

      ! Use bromine emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBREMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBREMIS

      ! Use fixed N2O surface mixing ratio?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LN2OEMIS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LN2OEMIS

      ! Separator line (start of initial global mixing ratio options)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Initialize strat H2O to GEOS-Chem baseline?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSETH2O', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSETH2O

      ! Starting year for CFC emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CFCYEAR', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(i4)') Input_Opt%CFCYEAR

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error check logical flags
      !=================================================================

      ! Fixed surface mixing ratio (basic)
      Input_Opt%LBASICEMIS = ( Input_Opt%LOCSEMIS .or.
     &                         Input_Opt%LCFCEMIS .or.
     &                         Input_Opt%LBREMIS  .or.
     &                         Input_Opt%LN2OEMIS )
      
      ! Turn off full-chem only switches 
      IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
         Input_Opt%LSETH2O = .FALSE.
      ENDIF

#if defined( MASSCONS )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%  MASS CONSERVATION TEST:                                  %%%
      !%%%  Always turn off emissions for the mass conservation      %%%
      !%%%  test.  (ewl, 6/24/15)                                    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

      ! Reset quantities
      Input_Opt%LEMIS = .FALSE.

      ! Print info to stdout stating that rint results
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) REPEAT( '%', 79 )
      WRITE( 6, 200   ) 
      WRITE( 6, 205   )
      WRITE( 6, 210   )
      WRITE( 6, '(a)' ) REPEAT( '%', 79 )
 200  FORMAT( '%%% MASS CONSERVATION TEST'                             )
 205  FORMAT( '%%% Automatically reset these EMISSIONS MENU settings:' )
 210  FORMAT( '%%% LEMIS is now FALSE'                                 )
#endif

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
        WRITE( 6, '(/,a)' ) 'EMISSIONS MENU'
        WRITE( 6, '(  a)' ) '--------------'
        WRITE( 6, 100 ) 'Turn on emissions?          : ',
     &                   Input_Opt%LEMIS
        WRITE( 6, 130 ) 'HEMCO Configuration file    : ',
     &                   TRIM( Input_Opt%HcoConfigFile )
        WRITE( 6, 100 ) 'Set day PBL BrO conc to 1ppt: ',
     &                   Input_Opt%LFIX_PBL_BRO
        WRITE( 6, 100 ) 'Use CH4 emissions inventory?: ',
     &                   Input_Opt%LCH4EMIS
        WRITE( 6, 100 ) 'Turn on surface BC for      : '
        WRITE( 6, 100 ) '      --> methane?          : ',
     &                         Input_Opt%LCH4SBC
        WRITE( 6, 100 ) '      --> OCS?              : ',
     &                         Input_Opt%LOCSEMIS
        WRITE( 6, 100 ) '      --> CFCs/HCFCs/halon? : ',
     &                         Input_Opt%LCFCEMIS
        WRITE( 6, 100 ) '      --> chlorine species? : ',
     &                         Input_Opt%LCLEMIS
        WRITE( 6, 100 ) '      --> bromine species?  : ',
     &                         Input_Opt%LBREMIS
        WRITE( 6, 100 ) '      --> N2O?              : ',
     &                         Input_Opt%LN2OEMIS
        WRITE( 6, 100 ) 'Set initial global MRs for  : '
        WRITE( 6, 100 ) '      --> strat. H2O?       : ',
     &                         Input_Opt%LSETH2O
        WRITE( 6, 110 ) 'CFC emissions read for year : ',
     &                   Input_Opt%CFCYEAR
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
 130  FORMAT( A, A   )

      END SUBROUTINE READ_EMISSIONS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_co_sim_menu
!
! !DESCRIPTION: Subroutine READ\_CO\_SIM\_MENU reads the CO SIM MENU 
!  section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CO_SIM_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  24 Mar 2017 - J. Fisher   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N
      CHARACTER(LEN=255) :: MSG, LOC

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CO_SIM_MENU begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Location for error messages
      LOC = ' -> at READ_CO_SIM_MENU (in GeosCore/input_mod.f)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         MSG = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &         'must be read in first!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

      ! Use P(CO) from CH4 from full chemistry simulation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_CH4', RC )
      READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_CH4

      ! Use P(CO) from CH4 from full chemistry simulation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPCO_NMVOC', RC )
      READ( SUBSTRS(1:N), * ) Input_Opt%LPCO_NMVOC

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( Input_Opt%ITS_A_TAGCO_SIM .and. am_I_Root ) THEN
         WRITE(6,'(/,a)') 'CO SIMULATION MENU ' //
     &        '(overwrites any other settings related to CO)'
         WRITE(6,'(  a)') '-------------------------------------'
         WRITE(6,100    ) 'Use full chem. P(CO) from CH4?:', 
     &                    Input_Opt%LPCO_CH4
         WRITE(6,100    ) '                   from NMVOC?:', 
     &                    Input_Opt%LPCO_NMVOC
         WRITE(6,'(  a)') '-------------------------------------'
      ENDIF      

      ! FORMAT statements
  90  FORMAT( A )
 100  FORMAT( A, L5 )
 110  FORMAT( A, L5, A )

      END SUBROUTINE READ_CO_SIM_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  03 Jun 2014 - R. Yantosca - Now specify biomass, biofuel options in HEMCO
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  13 Apr 2015 - R. Nassar - Simplified CO2 menu since options are now in HEMCO
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CO2_SIM_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_CO2_SIM_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Use Fossil Fuel emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSIL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSIL

      ! Use Ocean Exchange?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LOCEAN', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LOCEAN

       ! Turn on (balanced) biosphere with diurnal cycle?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIODIURNAL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBIODIURNAL

      ! Use Net Terrestrial Exchange Climatology?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIONETCLIM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBIONETCLIM

      ! Turn on Ship emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIP', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSHIP

      ! Turn on Aviation emissions?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LPLANE

      ! Turn on CO2 3D chemical source and surface correction?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCHEMCO2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCHEMCO2

      ! Tagged CO2 Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Background CO2 (no emissions or exchange) for Tagged-CO2 runs
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFFBKGRD',RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LFFBKGRD

      ! Turn on biosphere and ocean exchange region tagged species?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LBIOSPHTAG', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LBIOSPHTAG

      ! Turn on fossil fuel emission region tagged species?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFOSSILTAG', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LFOSSILTAG

      ! Turn on global ship emissions tagged species?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSHIPTAG', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSHIPTAG

      ! Turn on global aircraft emissions tagged species?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPLANETAG', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LPLANETAG

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( Input_Opt%ITS_A_CO2_SIM .and. am_I_Root ) THEN
         WRITE(6,'(/,a)') 'CO2 SIMULATION MENU ' //
     &        '(overwrites any other settings related to CO2)'
         WRITE(6,'(  a)') '-------------------------------------'
         WRITE(6,100    ) 'National Fossil Fuel Emission :',
     &                     Input_Opt%LFOSSIL
         WRITE(6,100    ) 'Ocean CO2 Uptake/Emission     :',
     &                     Input_Opt%LOCEAN
         WRITE(6,100    ) 'Biosphere seas/diurnal cycle  :',
     &                     Input_Opt%LBIODIURNAL
         WRITE(6,100    ) 'Net Terr Exch - Climatology   :',
     &                     Input_Opt%LBIONETCLIM
         WRITE(6,100    ) 'Intl/Domestic Ship emissions  :',
     &                     Input_Opt%LSHIP 
         WRITE(6,100    ) 'Intl/Domestic Aviation emiss  :',
     &                     Input_Opt%LPLANE 
         WRITE(6,100    ) 'CO2 from oxidation (CO,CH4,..):',
     &                     Input_Opt%LCHEMCO2
         WRITE(6, 90    ) 'Tagged CO2 settings'
         WRITE(6,100    ) '  Save Fossil CO2 in Bckgrnd  :',
     &                       Input_Opt%LFFBKGRD 
         WRITE(6,100    ) '  Tag Biosphere/Ocean CO2     :',
     &                       Input_Opt%LBIOSPHTAG 
         WRITE(6,100    ) '  Tag Fossil Fuel CO2         :',
     &                       Input_Opt%LFOSSILTAG 
         WRITE(6,100    ) '  Tag Global Ship CO2         :',
     &                       Input_Opt%LSHIPTAG
         WRITE(6,100    ) '  Tag Global Aviation CO2     :',
     &                       Input_Opt%LPLANETAG
         WRITE(6,'(  a)') '-------------------------------------'
      ENDIF      

      ! FORMAT statements
  90  FORMAT( A )
 100  FORMAT( A, L5 )

      END SUBROUTINE READ_CO2_SIM_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  27 Jul 2016 - M. Sulprizio- Remove LSVCSPEC and SPEC_RST_FILE. Restart files
!                              are now always saved out and the output restart
!                              file name is hardcoded in restart_mod.F.
!  14 Sep 2017 - M. Sulprizio- Add all options for overhead O3 used in FAST-JX
!                              (USE_ONLINE_O3, USE_O3_FROM_MET, & USE_TOMS_O3)
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CHEMISTRY_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Chemistry_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCHEM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCHEM

      ! Turn on stratospheric chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSCHEM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSCHEM

      ! Use Linoz for stratospheric ozone? (Otherwise, Synoz is used)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLINOZ', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LLINOZ

      ! Use Synoz if Linoz is turned off
      IF ( .not. Input_Opt%LLINOZ ) Input_Opt%LSYNOZ = .TRUE.

      ! Turn on unified strat-trop chemistry?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LUCX', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LUCX

      ! Turn on online stratospheric H2O?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LACTIVEH2O', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LACTIVEH2O

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Use online ozone in extinction calculations for FAST-JX?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_ONLINE_O3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%USE_ONLINE_O3

      ! Use ozone columns from met fields?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_O3_FROM_MET', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%USE_O3_FROM_MET

      ! Use ozone columns from TOMS?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_TOMS_O3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%USE_TOMS_O3

      ! GAMMA HO2 ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GAMMA_HO2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%GAMMA_HO2

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error check settings
      !=================================================================       

      ! Cannot use Synoz with linearized mesospheric chemistry
      IF ( Input_Opt%LUCX .and. Input_Opt%LSCHEM ) THEN
         IF (.not.Input_Opt%LLINOZ) THEN
            ErrMsg = 'Cannot use Synoz with linearized meso. chem.!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      ! Cannot have active H2O without stratospheric chemistry
      IF ( (.not.Input_Opt%LUCX) .and. Input_Opt%LACTIVEH2O ) THEN
         ErrMsg = 'Cannot have active H2O without full strat chem!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! FAST-JX is only used for fullchem, offline aerosol, and CH3I
      IF ( Input_Opt%ITS_A_FULLCHEM_SIM  .or.
     &     Input_Opt%ITS_AN_AEROSOL_SIM  .or.
     &     Input_Opt%ITS_A_CH3I_SIM      ) THEN

         ! Make sure either O3 from met or TOMS is selected
         IF ( .not. Input_Opt%USE_O3_FROM_MET .and.
     &        .not. Input_Opt%USE_TOMS_O3 ) THEN
            ErrMsg = 'Must select either O3 from met or TOMS/SBUV O3'
     &            // 'for O3 values above the chemistry grid!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
         IF ( Input_Opt%USE_O3_FROM_MET .and.
     &        Input_Opt%USE_TOMS_O3 ) THEN
            ErrMsg = 'Must select either O3 from met or TOMS/SBUV O3'
     &            // 'for O3 values above the chemistry grid!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Make sure specialty simulations select O3 from met or TOMS
         IF ( Input_Opt%ITS_AN_AEROSOL_SIM  .or.
     &        Input_Opt%ITS_A_CH3I_SIM      ) THEN
            IF ( Input_Opt%USE_ONLINE_O3 ) THEN
              ErrMsg= 'Cannot use online O3 for specialty simulations! '
     &             // 'Select O3 from met or TOMS O3 instead.'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
            ENDIF
         ENDIF

      ELSE

         Input_Opt%USE_ONLINE_O3   = .FALSE.
         Input_Opt%USE_O3_FROM_MET = .FALSE.
         Input_Opt%USE_TOMS_O3     = .FALSE.

      ENDIF

      ! Return success
      RC                      = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CHEMISTRY MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Turn on chemistry?          : ',
     &                        Input_Opt%LCHEM
         WRITE( 6, 100     ) 'Use linear. strat. chem?    : ',
     &                        Input_Opt%LSCHEM
         WRITE( 6, 100     ) ' => Use Linoz for O3?       : ',
     &                        Input_Opt%LLINOZ
         WRITE( 6, 100     ) 'Enable UCX?                 : ',
     &                        Input_Opt%LUCX
         WRITE( 6, 100     ) 'Online strat. H2O?          : ',
     &                        Input_Opt%LACTIVEH2O
         WRITE( 6, 100     ) 'Online ozone for FAST-JX?   : ',
     &                        Input_Opt%USE_ONLINE_O3
         WRITE( 6, 100     ) 'Ozone from met for FAST-JX? : ',
     &                        Input_Opt%USE_O3_FROM_MET
         WRITE( 6, 100     ) 'TOMS/SBUV ozone for FAST-JX?: ',
     &                        Input_Opt%USE_TOMS_O3
         WRITE( 6, 110     ) 'GAMMA HO2                   : ',
     &                        Input_Opt%GAMMA_HO2
         IF ( Input_Opt%USE_ONLINE_O3 ) THEN
         WRITE( 6, '(a)' ) ''
         WRITE( 6, '(a)' ) 'NOTE ABOUT OVERHEAD O3 FOR FAST-JX:'
         WRITE( 6, '(a)' ) ' Online O3 from GEOS-Chem will be used'
         WRITE( 6, '(a)' ) ' to weight the O3 column within the'
         WRITE( 6, '(a)' ) ' chemistry grid and O3 from met or TOMS'
         WRITE( 6, '(a)' ) ' will be used outside the chemistry grid.'
         ENDIF
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, F4.2 )

      END SUBROUTINE READ_CHEMISTRY_MENU  
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_radiation_menu
!
! !DESCRIPTION: Subroutine READ\_RADIATION\_MENU reads the RADIATION
! MENU section of the GEOS-Chem input file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_RADIATION_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  Flux outputs are now scheduled in the HISTORY.rc file, and the relevant
!  fields of Input_Opt will be populated in the RRTMG module routine
!  Init_RRTMG_Indices (called at startup).
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  03 Dec 2014 - M. Sulprizio- Now save fields to the Input_Opt object
!  10 Dec 2014 - M. Sulprizio- Add error checks for RRTMG switches
!  10 Mar 2017 - C. Keller   - Wrapped print out in am_I_Root statement.
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  23 Mar 2018 - E. Lundgren - Store wavelength data in Input_Opt instead of 
!                              variables formerly declared in CMN_FJX_MOD.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_RADIATION_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Radiation_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! AOD wavelength selection?
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%NWVSELECT, -1, 
     &                     'Wavelengths', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%NWVSELECT
         READ( SUBSTRS(N), * ) Input_Opt%WVSELECT(N)
         ! save the string version also
         Input_Opt%STRWVSELECT(N) = TRIM(SUBSTRS(N))
      ENDDO

      ! Turn on RRTMG?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LRAD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LRAD

      ! Turn on LW radiation calculation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LLWRAD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LLWRAD

      ! Turn on SW radiation calculation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSWRAD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSWRAD

      ! Calculate for clear-sky?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Clear sky', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSKYRAD(1)

      ! Calculate for all-sky?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'All Sky', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSKYRAD(2)

      ! Radiation timestep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS_RAD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%TS_RAD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error check settings
      !=================================================================       

      ! Use of RRTMG necessitates recompilation 
#if !defined( RRTMG )
      IF ( Input_Opt%LRAD ) THEN
         ErrMsg = 'LRAD=T but RRTMG not defined at compile time!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#else
      IF ( .not. Input_Opt%LRAD ) THEN
         ErrMsg = 'LRAD=F but RRTMG defined at compile time!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      ! Make sure radiation switches are turned off if RRTMG is off
      IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LLWRAD ) THEN
         ErrMsg = 'Cannot have LW fluxes turned on without RRTMG'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSWRAD ) THEN
         ErrMsg = 'Cannot have SW fluxes turned on without RRTMG'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(1) ) THEN
         ErrMsg = 'Cannot have clear-sky flux turned on without RRTMG'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      IF ( ( .not. Input_Opt%LRAD ) .and. Input_Opt%LSKYRAD(2) ) THEN
         ErrMsg = 'Cannot have all-sky flux turned on without RRTMG'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'RADIATION MENU'
         WRITE( 6, '(  a)' ) '--------------'
         DO N=1, Input_Opt%NWVSELECT
            WRITE( 6, 115     ) 'AOD output wavelength (nm)  : ', 
     &                           Input_Opt%WVSELECT(N)
         ENDDO
         WRITE( 6, 100 ) 'Turn on radiation?          : ',
     &                    Input_Opt%LRAD
         WRITE( 6, 100 ) 'Consider LW                 : ',
     &                    Input_Opt%LLWRAD
         WRITE( 6, 100 ) 'Consider SW                 : ',
     &                    Input_Opt%LSWRAD
         WRITE( 6, 125 ) 'Clear-sky/All-sky           : ',
     &                    Input_Opt%LSKYRAD(1), '/',
     &                    Input_Opt%LSKYRAD(2)
         WRITE( 6, 110 ) 'Radiation timestep [sec]    : ',
     &                    Input_Opt%TS_RAD
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5           )
 110  FORMAT( A, I5           )
 115  FORMAT( A, F7.1         ) 
 120  FORMAT( A, 11I1         )
 125  FORMAT( A, L5, A, L5    )
 130  FORMAT( A, 12( A2, 1x ) )

      END SUBROUTINE READ_RADIATION_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  20 Sep 2016 - R. Yantosca - Use "I8" format for write statement
!  18 Apr 2017 - C. Holmes   - Fix incorrect format statement in error msg
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_TRANSPORT_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Transport_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN
         ErrMSg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on transport?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LTRAN', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LTRAN

      ! Fill negative values
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LFILL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LFILL

      ! IORD, JORD, KORD
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 3, 'IORD, JORD, KORD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%TPCORE_IORD,
     &                        Input_Opt%TPCORE_JORD,
     &                        Input_Opt%TPCORE_KORD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'TRANSPORT MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Turn on transport?          : ',
     &                        Input_Opt%LTRAN
         WRITE( 6, 100     ) 'Let TPCORE Fill negatives?  : ',
     &                        Input_Opt%LFILL
         WRITE( 6, 110     ) 'IORD, JORD, KORD for TPCORE?: ',
     &                        Input_Opt%TPCORE_IORD, 
     &                        Input_Opt%TPCORE_JORD,
     &                        Input_Opt%TPCORE_KORD
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, 5I5 )
      
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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  12 Aug 2015 - R. Yantosca - Now allow non-local PBL mixing for MERRA2 met
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CONVECTION_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ThisLoc = 
     &   ' -> at Read_Convection_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ErrMsg  = 'Error reading the "input.geos" file!'

      ! Turn on convection?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LCONV', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LCONV

      ! Turn on BL mixing
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LTURB', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LTURB

      ! Turn on non-local PBL scheme (Lin, 03/31/09)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LNLPBL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LNLPBL

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CONVECTION MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100     ) 'Turn on cloud convection?   : ',
     &                        Input_Opt%LCONV
         WRITE( 6, 100     ) 'Turn on PBL mixing?         : ',
     &                        Input_Opt%LTURB
         WRITE( 6, 100     ) 'Turn on non-local PBL?      : ',
     &                        Input_Opt%LNLPBL
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
      
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  25 Jun 2014 - R. Yantosca - Now set USE_OLSON_2001 = T for GEOS-FP met;
!                              this was originally done in INIT_INPUT
!  05 Mar 2015 - C. Keller   - Added PBL_DRYDEP to define if dry 
!                              deposition rates are calculated over
!                              full PBL or 1st model level only.
!                              Set to .TRUE. (full PBL) if full PBL
!                              mixing is used, to false otherwise.
!  13 Sep 2017 - M. Sulprizio- Remove Input_Opt%USE_OLSON_2001. Olson 2001 is
!                              now the default.
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DEPOSITION_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Deposition_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &          'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on drydep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDRYD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LDRYD

      ! Turn on wetdep?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWETD', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_OPt%LWETD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error check settings
      !=================================================================

      ! Turn off drydep for simulations that don't need it
      IF ( Input_Opt%ITS_A_CH3I_SIM    ) Input_Opt%LDRYD = .FALSE.
      IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LDRYD = .FALSE.

      ! Turn off wetdep for simulations that don't need it
      IF ( Input_Opt%ITS_A_CH3I_SIM    ) Input_Opt%LWETD = .FALSE.
      IF ( Input_Opt%ITS_A_HCN_SIM     ) Input_Opt%LWETD = .FALSE.
      IF ( Input_Opt%ITS_A_TAGO3_SIM   ) Input_Opt%LWETD = .FALSE.
      IF ( Input_Opt%ITS_A_TAGCO_SIM   ) Input_Opt%LWETD = .FALSE.
      IF ( Input_Opt%ITS_A_C2H6_SIM    ) Input_Opt%LWETD = .FALSE.
      IF ( Input_Opt%ITS_A_CH4_SIM     ) Input_Opt%LWETD = .FALSE.

      ! Set the PBL drydep flag. This determines if dry deposition is 
      ! applied (and drydep frequencies are calculated) over the entire 
      ! PBL or the first model layer only. For now, set this value 
      ! automatically based upon the selected PBL scheme: 1st model layer 
      ! for the non-local PBL scheme, full PBL for the full-mixing scheme.
      IF ( Input_Opt%LNLPBL ) THEN 
         Input_Opt%PBL_DRYDEP = .FALSE.
      ELSE
         Input_Opt%PBL_DRYDEP = .TRUE.
      ENDIF 

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'DEPOSITION MENU'
         WRITE( 6, '(  a)' ) '---------------'
         WRITE( 6, 100     ) 'Turn on dry deposition?     : ', 
     &                        Input_Opt%LDRYD
         WRITE( 6, 100     ) 'Dry dep over full PBL?      : ', 
     &                        Input_Opt%PBL_DRYDEP
         WRITE( 6, 100     ) 'Turn on wet deposition?     : ', 
     &                        Input_Opt%LWETD
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_GAMAP_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Gamap_Menu (in module GeosCore/input_mod.F)'

      ! Background
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DIAGINFO', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%GAMAP_DIAGINFO

      ! Redirect
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TRACERINFO', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%GAMAP_TRACERINFO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'GAMAP MENU'
         WRITE( 6, '(  a)' ) '---------------'            
         WRITE( 6, '(a,a)' ) 'GAMAP "diaginfo.dat"   file : ', 
     &                        TRIM( Input_Opt%GAMAP_DIAGINFO   )
         WRITE( 6, '(a,a)' ) 'GAMAP "tracerinfo.dat" file : ',
     &                        TRIM( Input_Opt%GAMAP_TRACERINFO )
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
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
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
!  NJDAY is now located in CMN_SIZE_mod.F, so that we can eventually
!  retire the obsolete CMN_DIAG_mod.F. (bmy, 1/16/18)
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  03 Aug 2012 - R. Yantosca - IU_GEOS is now a global module variable
!  01 Nov 2012 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  28 Feb 2013 - R. Yantosca - Don't call IS_LAST_DAY_GOOD when using ESMF
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: IOS

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! READ_OUTPUT_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Output_Menu (in module GeosCore/input_mod.F)'

      ! Read info
      READ( IU_GEOS, 100, IOSTAT=IOS ) Input_Opt%NJDAY
 100  FORMAT( 26x, 31i1, /  26x, 29i1, /, 26x, 31i1, /, 26x, 30i1, /, 
     &        26x, 31i1, /, 26x, 30i1, /, 26x, 31i1, /, 26x, 31i1, /,
     &        26x, 30i1, /  26x, 31i1, /, 26x, 30i1, /, 26x, 31i1 )

      ! Trap potential errors
      IF ( IOS /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      
      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'OUTPUT MENU'
         WRITE( 6, '(  a)' ) '-----------'
         WRITE( 6, 110     )
         WRITE( 6, 120     )
         WRITE( 6, 130     )
         WRITE( 6, 140     ) Input_Opt%NJDAY
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
      !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
      !
      ! When we connect to the GEOS-5 GCM via ESMF, we let the GCM
      ! handle all timing and diagnostic output via the GEOS-5 History
      ! Component.  Therefore we do not need to make sure that we have 
      ! scheduled output for the last day of the run.  This can cause
      ! a premature exit. (bmy, 2/28/13)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
      !
      ! Current practice in the std GEOS-Chem is to not let the run
      ! proceed unless the user has scheduled output for the last day.
      ! (bmy, 2/28/13)
      !-----------------------------------------------------------------

#if defined ( BPCH_DIAG )
      ! Make sure we have output at end of run
      CALL IS_LAST_DAY_GOOD( Input_Opt, RC )
#endif

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
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
      USE CMN_SIZE_MOD        ! Size parameters
      USE CMN_DIAG_MOD        ! Needed for timeseries diags (binary only)
#if defined( BPCH_DIAG )
      USE BPCH2_MOD,     ONLY : OPEN_BPCH2_FOR_WRITE
      USE DIAG03_MOD,    ONLY : ND03,      PD03,      PD03_PL
      USE DIAG04_MOD,    ONLY : ND04,      PD04
      USE DIAG41_MOD,    ONLY : ND41,      PD41
      USE DIAG42_MOD,    ONLY : ND42,      PD42
      USE DIAG53_MOD,    ONLY : ND53,      PD53
      USE DIAG56_MOD,    ONLY : ND56,      PD56
#endif
      USE DRYDEP_MOD,    ONLY : NUMDEP
      USE ErrCode_Mod
      USE FILE_MOD,      ONLY : IU_BPCH
      USE Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,      ONLY : GET_NYMDb, GET_NHMSb, EXPAND_DATE
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
!  06 Mar 2013 - H. Amos     - merge C. Friedman's POP code
!  08 Nov 2013 - M. Sulprizio- Remove HR1_NO, and HR2_NO from ND43 diagnostic.
!  03 Feb 2014 - R. Yantosca - Remove references to TINDEX, TCOUNT, and TMAX
!                              from diag_mod.  They are in CMN_SIZE_mod.F.
!                              Lee Murray reports this causes the compilation
!                              to choke on Macintosh platforms.
!  23 Jun 2014 - R. Yantosca - Now pass Input_Opt, RC to INIT_DIAG_OH
!  24 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  15 Dec 2014 - M. Sulprizio- Moved radiation diagnostic from ND71 to ND72 to
!                              avoid conflicts with hourly max ppbv diagnostic.
!  15 Jan 2015 - R. Yantosca - Now define Input_Opt%DIAG_COLLECTION 
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
!  23 Jun 2016 - R. Yantosca - Now call INIT_DIAG_OH from GIGC_Init_Extra
!  20 Jul 2016 - R. Yantosca - Remove references to NNPAR
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: M, N, N_MAX, N_TMP, N_ADVECT

      ! Define shadow variables so that we can still use certain
      ! flags (e.g. ND70) even when compiling with BPCH_DIAG=n
      INTEGER            :: ND01x, ND02x, ND03x, ND04x, ND05x, ND06x
      INTEGER            :: ND07x, ND08x, ND09x, ND10x, ND11x, ND12x
      INTEGER            :: ND13x, ND14x, ND15x, ND16x, ND17x, ND18x
      INTEGER            :: ND19x, ND21x, ND22x, ND24x
      INTEGER            :: ND25x, ND26x, ND27x, ND28x, ND29x, ND30x
      INTEGER            :: ND31x, ND32x, ND33x, ND34x, ND35x, ND36x
      INTEGER            :: ND37x, ND38x, ND39x, ND41x, ND42x
      INTEGER            :: ND43x, ND44x, ND45x, ND46x, ND47x
      INTEGER            :: ND52x, ND53x, ND54x
      INTEGER            :: ND55x, ND56x, ND57x, ND58x, ND59x, ND60x
      INTEGER            :: ND61x, ND62x, ND64x, ND66x
      INTEGER            :: ND67x, ND68x, ND69x, ND70x, ND71x, ND72x
      INTEGER            :: ND73x, HR1_JVx, HR2_JVx, HR1_OHx, HR2_OHx
      INTEGER            :: HR1_OTHx, HR2_OTHx

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc, BPCH_FILE

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DIAGNOSTIC_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Diagnostic_Menu (in module GeosCore/input_mod.F)'

      ! Get fields from Input_Opt
      N_ADVECT = Input_Opt%N_ADVECT
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg= 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &           'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Get max number of species
      N_MAX = N_ADVECT

      IF ( .not. Input_Opt%LEMIS ) THEN
         WRITE( 6, '(a)' ) 'WARNING: Emissions are turned off. The'
         WRITE( 6, '(a)' ) ' following diagnostics will also be turned'
         WRITE( 6, '(a)' ) ' off:  ND01, ND04, ND06, ND07, ND08, ND11,'
         WRITE( 6, '(a)' ) ' ND13, ND28, ND29, ND32, ND34, ND36, ND46,'
         WRITE( 6, '(a)' ) ' ND53, ND56, ND58'
      ENDIF

      ! Binary punch file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'BPCH_FILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) BPCH_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !--------------------------
      ! ND01: Rn-Pb-Be source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND01', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND01x
      IF ( .not. Input_Opt%ITS_A_RnPbBe_SIM .or.
     &     .not. Input_Opt%LEMIS ) ND01x = 0 
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 01, ND01x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND02: Rn-Pb-Be decay
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND02', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND02x
      IF ( .not. Input_Opt%ITS_A_RnPbBe_SIM ) ND02x = 0 
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 02, ND02x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND03: Hg diagnostics
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND03', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND03x
      IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) ND03x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 03, ND03x, SUBSTRS(2:N), N-1, PD03 )
#endif

      !--------------------------
      ! ND04: CO2 emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND04', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND04x
      IF ( .not. Input_Opt%ITS_A_CO2_SIM .or.
     &     .not. Input_Opt%LEMIS ) ND04x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 04, ND04x, SUBSTRS(2:N), N-1, PD04 )
#endif

      !--------------------------
      ! ND05: Sulfate prod/loss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND05', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND05x
      IF ( .not. Input_Opt%LSULF ) ND05x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 05, ND05x, SUBSTRS(2:N), N-1, PD05 )
#endif

      !--------------------------
      ! ND06: Dust emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND06', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND06x
      IF ( .not. Input_Opt%LDUST .or.
     &     .not. Input_Opt%LEMIS ) ND06x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 06, ND06x, SUBSTRS(2:N), N-1, NDSTBIN)
#endif
 
      !--------------------------
      ! ND07: Carbon/SOA source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND07', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND07x
      IF ( .not. Input_Opt%LCARB .or.
     &     .not. Input_Opt%LEMIS ) ND07x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 07, ND07x, SUBSTRS(2:N), N-1, PD07 )
#endif

      !--------------------------
      ! ND08: Sea salt source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND08', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND08x
      IF ( .not. Input_Opt%LSSALT .or.
     &     .not. Input_Opt%LEMIS ) ND08x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 08, ND08x, SUBSTRS(2:N), N-1, PD08 )
#endif

      ! Included in HEMCO diagnostics for netCDF

      !--------------------------
      ! ND09: HCN/CH3CN (obsolete)
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND09', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND09x
      ND09x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 09, ND09x, SUBSTRS(2:N), N-1, PD09 )
#endif

      !--------------------------
      ! ND10: Free
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND10', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND10x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 10, ND10x, SUBSTRS(2:N), N-1, PD10 )
#endif

      !--------------------------
      ! ND11: Acetone source
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND11', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND11x
      IF ( .not. Input_Opt%LEMIS ) ND11x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 11, ND11x, SUBSTRS(2:N), N-1, PD11 )
#endif

      !--------------------------
      ! ND12: PBL distribution
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND12', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND12x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 12, ND12x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND13: Sulfur sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND13', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND13x
      IF ( .not. Input_Opt%LSULF .or.
     &     .not. Input_Opt%LEMIS ) ND13x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 13, ND13x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND14: Wet conv up flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND14', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND14x
      IF ( .not. Input_Opt%LCONV ) ND14x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 14, ND14x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND15: Mass change in PBL
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND15', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND15x
      IF ( .not. Input_Opt%LTURB ) ND15x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 15, ND15x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND16: Precip fractions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND16', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND16x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 16, ND16x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND17: Rainout losses
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND17', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND17x
      IF ( .not. Input_Opt%LWETD ) ND17x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 17, ND17x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND18: Washout losses
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND18', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND18x
      IF ( .not. Input_Opt%LWETD ) ND18x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 18, ND18x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND19: CH4 loss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND19', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND19x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 19, ND19x, SUBSTRS(2:N), N-1, 3 )
#endif

      !--------------------------
      ! ND21: Opt depths etc.
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND21', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND21x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 21, ND21x, SUBSTRS(2:N), N-1, PD21 )
#endif

      ! Error check 
      IF ( ND21x > 0 .and. Input_Opt%SALA_REDGE_um(2) /= 0.5 ) THEN
         ErrMsg = 'Cannot output seasalt AOD''s when radius bin' //
     &            ' is not split at 0.5 um!!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !--------------------------
      ! ND22: J-values
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND22', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND22x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 22, ND22x, SUBSTRS(2:N), N-1, PD22 )
#endif

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'separator 1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) HR1_JVx, HR2_JVx

      !--------------------------
      ! ND24: E/W transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND24', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND24x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 24, ND24x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND25: N/S transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND25', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND25x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 25, ND25x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND26: U/D transport flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND26', RC )
      READ( SUBSTRS(1), * ) ND26x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 26, ND26x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND27: Strat fluxes
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND27', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND27x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 27, ND27x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND28: Biomass emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND28', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND28x
      IF ( .not. Input_Opt%LEMIS ) ND28x = 0
 
      !--------------------------
      ! ND29: CO sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND29', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND29x
      IF ( .not. Input_Opt%LEMIS ) ND29x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 29, ND29, SUBSTRS(2:N), N-1, 5 )
#endif

      !--------------------------
      ! ND30: Land map
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND30', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND30x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 30, ND30x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND31: Surface pressure
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND31', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND31x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 31, ND31x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND32: NOx sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND32', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND32x
      IF ( .not. Input_Opt%LEMIS ) ND32x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 32, ND32x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND33: Column species
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND33', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND33x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 33, ND33x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND34: Biofuel sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND34', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND34x
      IF ( .not. Input_Opt%LEMIS ) ND34x = 0

      !--------------------------
      ! ND35: 500 hPa species
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND35', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND35x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 35, ND35x, SUBSTRS(2:N), N-1, N_MAX )
#endif

      !--------------------------
      ! ND36: Anthro emissions
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND36', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND36x
      IF ( .not. Input_Opt%LEMIS ) ND36x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 
     &                 36, ND36, SUBSTRS(2:N), N-1, N_ADVECT )
#endif

      !--------------------------
      ! ND37: Updraft scav frac
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND37', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND37x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 37, ND37x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND38: Cld conv losses
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND38', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND38x
      IF ( .not. Input_Opt%LWETD .and. .not. Input_Opt%LCONV ) THEN
         ND38x = 0
      ENDIF
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 38, ND38x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND39: Wet scav losses
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND39', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND39x
      IF ( .not. Input_Opt%LWETD ) ND39x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 39, ND39x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND41: Afternoon PBL ht
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND41', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND41x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 41, ND41x, SUBSTRS(2:N), N-1, PD41 )
#endif

      !--------------------------
      ! ND42: SOA concentrations
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND42', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND42x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 42, ND42x, SUBSTRS(2:N), N-1, PD42 )
#endif

      !--------------------------
      ! ND43: OH, HO2
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND43', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND43x
      IF ( .not. ( Input_Opt%ITS_A_FULLCHEM_SIM    .or.
     &             Input_Opt%ITS_A_CH4_SIM      ) ) ND43x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 43, ND43x, SUBSTRS(2:N), N-1, PD43 )
#endif

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'HR1_OH, HR2_OH', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) HR1_OHx, HR2_OHx

      !--------------------------
      ! ND44 drydep vel & flux
      !--------------------------

      ! Number of species depends on simulation type
      IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
         N_TMP = Input_Opt%N_ADVECT
      ELSE
         N_TMP = NUMDEP
      ENDIF

      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND44', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND44x
      IF ( .not. Input_Opt%LDRYD ) ND44x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 44, ND44, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND45: Species conc.
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND45', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND45x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 45, ND45x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'HR1_OTH, HR2_OTH', RC ) 
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) HR1_OTHx, HR2_OTHx

      !--------------------------
      ! ND46: Biogenic sources
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND46', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND46x
      IF ( .not. Input_Opt%LEMIS ) ND46x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 46, ND46, SUBSTRS(2:N), N-1, 26 )
#endif

      !--------------------------
      ! ND47: 24h avg species
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND47', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND47x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 
     &                 47, ND47, SUBSTRS(2:N), N-1, N_ADVECT )
#endif

      !--------------------------
      ! ND52: GAMMA values
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND52', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND52x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 52, ND52x, SUBSTRS(2:N), N-1, PD52 )
#endif

      !--------------------------
      ! ND53: POPS
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND53', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND53x
      IF ( .not. Input_Opt%ITS_A_POPS_SIM .or.
     &     .not. Input_Opt%LEMIS) ND53x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 53, ND53x, SUBSTRS(2:N), N-1, PD53 )
#endif

      !--------------------------
      ! ND54: Time in troposphere
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND54', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND54x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 54, ND54x, SUBSTRS(2:N), N-1, 1 )
#endif

      !--------------------------
      ! ND55: Tropopause diags
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND55', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND55x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 55, ND55x, SUBSTRS(2:N), N-1, PD55 )
#endif

      !--------------------------
      ! ND56: Lightning flashes
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND56', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND56x
      IF ( .not. Input_Opt%LEMIS ) ND56x = 0
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 56, ND56x, SUBSTRS(2:N), N-1, PD56 )
#endif

      !--------------------------
      ! ND57: Potential temperature
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND57', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND57x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 57, ND57x, SUBSTRS(2:N), N-1, PD57 )
#endif

      !--------------------------
      ! ND58: CH4 Emissions 
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND58', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND58x

      !--------------------------
      ! ND59: TOMAS aerosol emiss
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND59', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND59x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 59, ND59x, SUBSTRS(2:N), N-1, PD59 )
#endif

      !--------------------------
      ! ND60: Wetland Fraction
      ! ND60: TOMAS rate
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND60', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND60x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 60, ND60x, SUBSTRS(2:N), N-1, PD60 )
#endif

      !--------------------------
      ! ND61: 3-D TOMAS rate
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND61', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND61x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 61, ND61x, SUBSTRS(2:N), N-1, PD61 )
#endif

      !--------------------------
      ! ND62: Free
      !--------------------------
      N_TMP = N_ADVECT
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND62', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND62x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 62, ND62x, SUBSTRS(2:N), N-1, N_TMP )
#endif

      !--------------------------
      ! ND64: Radiative flux
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND64', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND64x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 64, ND64x, SUBSTRS(2:N), N-1, PD64*3 )
#endif

      !--------------------------
      ! ND66: DAO 3-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND66', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND66x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 66, ND66x, SUBSTRS(2:N), N-1, PD66 )
#endif

      !--------------------------
      ! ND67: DAO 2-D fields
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND67', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND67x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 67, ND67x, SUBSTRS(2:N), N-1, PD67 )
#endif

      !--------------------------
      ! ND68: Air masses etc
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND68', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND68x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 68, ND68x, SUBSTRS(2:N), N-1, PD68 )
#endif

      !--------------------------
      ! ND69: Surface areas
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND69', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND69x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 69, ND69x, SUBSTRS(2:N), N-1, PD69 )
#endif

      !--------------------------
      ! ND70: Debug info
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND70', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND70x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 70, ND70x, SUBSTRS(2:N), N-1, 0 )
#endif
     
      !--------------------------
      ! ND71: Maximum species
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND71', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND71x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 
     &                 71, ND71x, SUBSTRS(2:N), N-1, N_ADVECT )
#endif

      !--------------------------
      ! ND72: Radiation output
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND72', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND72x

#if defined( BPCH_DIAG )
      !output fields are nspecies*nradfields but user can only specify
      !rad fields (e.g. SW TOA ALL-SKY) so we set the max to the total
      !divided by number of allowed species (PD72R)
      CALL SET_TINDEX( am_I_Root, 72, ND72, SUBSTRS(2:N), N-1, PD72R )
#endif

      !If LRAD is on then ND72 must be on (so the diagnostic is
      !available to write into). Check for this
      IF ( (Input_Opt%LRAD) .AND. (ND72x.EQ.0) ) THEN
         ErrMsg = 'If LRAD is true then ' //
     &            'ND72 diagnostic must be switched on'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !--------------------------
      ! ND73: ISORROPIA pH and H+
      !--------------------------
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ND73', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1), * ) ND73x
#if defined( BPCH_DIAG )
      CALL SET_TINDEX( am_I_Root, 73, ND73x, SUBSTRS(2:N), N-1, PD73 )
#endif

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

#if defined( BPCH_DIAG )
      !=================================================================
      ! %%%%% IF BPCH DIAGNOSTICS ARE ACTIVATED (BPCH_DIAG=y) %%%%%
      !
      ! Copy shadow variables to diagnostic variables and 
      ! call various bpch diagnostic setup routines
      !=================================================================
      ND01                         = ND01x
      ND02                         = ND02x
      ND03                         = ND03x
      ND04                         = ND04x
      ND05                         = ND05x
      ND06                         = ND06x
      ND07                         = ND07x
      ND08                         = ND08x
      ND09                         = ND09x
      ND10                         = ND10x
      ND11                         = ND11x
      ND12                         = ND12x
      ND13                         = ND13x
      ND14                         = ND14x
      ND15                         = ND15x
      ND16                         = ND16x
      ND17                         = ND17x
      ND18                         = ND18x
      ND19                         = ND19x
      ND21                         = ND21x
      ND22                         = ND22x
      ND24                         = ND24x
      ND25                         = ND25x
      ND26                         = ND26x
      ND27                         = ND27x
      ND28                         = ND28x
      ND29                         = ND29x
      ND30                         = ND30x
      ND31                         = ND31x
      ND32                         = ND32x
      ND33                         = ND33x
      ND34                         = ND34x
      ND35                         = ND35x
      ND36                         = ND36x
      ND37                         = ND37x
      ND38                         = ND38x
      ND39                         = ND39x
      ND41                         = ND41x
      ND42                         = ND42x
      ND43                         = ND43x
      ND44                         = ND44x
      ND45                         = ND45x
      ND46                         = ND46x
      ND47                         = ND47x
      ND52                         = ND52x
      ND53                         = ND53x
      ND54                         = ND54x
      ND55                         = ND55x
      ND56                         = ND56x
      ND57                         = ND57x
      ND59                         = ND59x
      ND60                         = ND60x
      ND61                         = ND61x
      ND62                         = ND62x
      ND64                         = ND64x
      ND66                         = ND66x
      ND67                         = ND67x
      ND68                         = ND68x
      ND69                         = ND69x
      ND70                         = ND70x
      ND71                         = ND71x
      ND72                         = ND72x
      ND73                         = ND73x
      Input_Opt%ND01               = ND01x
      Input_Opt%ND02               = ND02x
      Input_Opt%ND03               = ND03x
      Input_Opt%ND04               = ND04x
      Input_Opt%ND05               = ND05x
      Input_Opt%ND06               = ND06x
      Input_Opt%ND07               = ND07x
      Input_Opt%ND08               = ND08x
      Input_Opt%ND09               = ND09x
      Input_Opt%ND10               = ND10x
      Input_Opt%ND11               = ND11x
      Input_Opt%ND12               = ND12x
      Input_Opt%ND13               = ND13x
      Input_Opt%ND14               = ND14x
      Input_Opt%ND15               = ND15x
      Input_Opt%ND16               = ND16x
      Input_Opt%ND17               = ND17x
      Input_Opt%ND18               = ND18x
      Input_Opt%ND19               = ND19x
      Input_Opt%ND21               = ND21x
      Input_Opt%ND22               = ND22x
      Input_Opt%ND24               = ND24x
      Input_Opt%ND25               = ND25x
      Input_Opt%ND26               = ND26x
      Input_Opt%ND27               = ND27x
      Input_Opt%ND28               = ND28x
      Input_Opt%ND29               = ND29x
      Input_Opt%ND30               = ND30x
      Input_Opt%ND31               = ND31x
      Input_Opt%ND32               = ND32x
      Input_Opt%ND33               = ND33x
      Input_Opt%ND34               = ND34x
      Input_Opt%ND35               = ND35x
      Input_Opt%ND36               = ND36x
      Input_Opt%ND37               = ND37x
      Input_Opt%ND38               = ND38x
      Input_Opt%ND39               = ND39x
      Input_Opt%ND41               = ND41x
      Input_Opt%ND42               = ND42x
      Input_Opt%ND43               = ND43x
      Input_Opt%ND44               = ND44x
      Input_Opt%ND45               = ND45x
      Input_Opt%ND46               = ND46x
      Input_Opt%ND47               = ND47x
      Input_Opt%ND52               = ND52x
      Input_Opt%ND53               = ND53x
      Input_Opt%ND54               = ND54x
      Input_Opt%ND55               = ND55x
      Input_Opt%ND56               = ND56x
      Input_Opt%ND57               = ND57x
      Input_Opt%ND59               = ND59x
      Input_Opt%ND60               = ND60x
      Input_Opt%ND61               = ND61x
      Input_Opt%ND62               = ND62x
      Input_Opt%ND64               = ND64x
      Input_Opt%ND66               = ND66x
      Input_Opt%ND67               = ND67x
      Input_Opt%ND68               = ND68x
      Input_Opt%ND69               = ND69x
      Input_Opt%ND70               = ND70x
      Input_Opt%ND71               = ND71x
      Input_Opt%ND72               = ND72x
      Input_Opt%ND73               = ND73x
      Input_Opt%LPRT               = ( ND70x > 0 )
      HR1_JV                       = HR1_JVx 
      HR2_JV                       = HR2_JVx
      HR1_OH                       = HR1_OHx
      HR2_OH                       = HR2_OHx
      HR1_OTH                      = HR1_OTHx
      HR2_OTH                      = HR2_OTHx

      ! Loop over # of diagnostics
      DO M = 1, Input_Opt%Max_BPCH_Diag
         Input_Opt%TCOUNT(M)       = TCOUNT(M)
         Input_Opt%TMAX(M)         = TMAX(M)

         ! Loop over tracers per diagnostic
         DO N = 1, Input_Opt%N_ADVECT
            Input_Opt%TINDEX(M,N)  = TINDEX(M,N)
         ENDDO
      ENDDO

      !=================================================================
      ! Call other setup routines
      !================================================================

      ! Expand YYYYMMDD tokens in the bpch file name
      CALL EXPAND_DATE( BPCH_FILE, GET_NYMDb(), GET_NHMSb() )

      ! Find a free file LUN
      IU_BPCH = findFreeLUN()

      ! Open the binary punch file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_BPCH, BPCH_FILE )

#else

      !=================================================================
      ! %%%%% IF BPCH DIAGNOSTICS ARE DEACTIVATED (BPCH_DIAG=n) %%%%%
      !
      ! When BPCH_DIAG=n, we still need to define the Input_Opt%LPRT 
      ! flag for debugging.  We will eventually move this elsewhere.
      ! Also set values for time series diagnostics which can still be
      ! output as binary if other diangostics are netcdf.
      !=================================================================
      Input_Opt%ND70 = ND70x
      Input_Opt%LPRT = ( ND70x > 0 )

#endif
#if defined( MODEL_GEOS )
      Input_Opt%ND21 = ND21x
#endif

      ! Return success
      RC = GC_SUCCESS

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
#if defined( BPCH_DIAG )
#if   defined( TOMAS )
      USE CHARPAK_MOD, ONLY : TXTEXT   ! (win, 7/14/09)
#endif
      USE CMN_DIAG_MOD                              ! TMAX, TINDEX
      USE CMN_SIZE_MOD                              ! Size parameters
#endif
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
! !REMARKS:
!  NOTE: This routine converts to a stub when BPCH_DIAG=n (bmy, 1/16/18)
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
#if defined( BPCH_DIAG )
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
#endif

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
      USE CMN_DIAG_MOD          ! ND40
      USE CMN_SIZE_MOD          ! MAXFAM
      USE ErrCode_Mod
      USE Input_Opt_Mod,   ONLY : OptInput
      USE PLANEFLIGHT_MOD, ONLY : SET_PLANEFLIGHT

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
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_PLANEFLIGHT_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ND40  = 0
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &  ' -> at Read_Planeflight_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on planeflight diagnostic?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DO_PF', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_PF

      ! Set ND40 flag from DO_PF
      IF ( Input_Opt%DO_PF ) ND40 = 1

      ! Input file name (w/ flight track data points)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'IFILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%PF_IFILE

      ! Output file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'OFILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%PF_OFILE

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'PLANEFLIGHT MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100 ) 'Turn on planeflight diag?   : ',
     &                    Input_Opt%DO_PF
         WRITE( 6, 110 ) 'Flight track input file     : ',
     &                    TRIM( Input_Opt%PF_IFILE )
         WRITE( 6, 110 ) 'Output file name            : ',
     &                    TRIM( Input_Opt%PF_OFILE )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )

      !=================================================================
      ! Call setup routines from other F90 modules
      !=================================================================

      ! Pass variables to "planeflight_mod.f"
      CALL SET_PLANEFLIGHT( Input_Opt%DO_PF,
     &                      Input_Opt%PF_IFILE,
     &                      Input_Opt%PF_OFILE )
   
      END SUBROUTINE READ_PLANEFLIGHT_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_obspack_menu
!
! !DESCRIPTION: Subroutine READ\_OBSPACK\_MENU reads the OBSPACK MENU
!  section of the GEOS-Chem input file.  This turns on the ObsPack diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_OBSPACK_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  04 Jan 2019 - R. Yantosca - Add extra ObsPack fields from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, S

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_OBSPACK_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc =
     &  ' -> at Read_ObsPack_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ObsPack diagnostic?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Do_ObsPack', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%Do_ObsPack

      ! ObsPack quiet output?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'Do_ObsPack', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ObsPack_Quiet

      ! Input file name (w/ coordinates and sampling strategy)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ObsPack_InputFile', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ObsPack_InputFile

      ! Output file name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ObsPack_OutputFile', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ObsPack_OutputFile

      ! Species names
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'ObsPack_SpecName', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      
      ! Populate the ObsPack species name list
      Input_Opt%ObsPack_nSpc = N
      DO S = 1, Input_Opt%ObsPack_nSpc
         Input_Opt%ObsPack_SpcName(S) = TRIM( SUBSTRS(S) )
      ENDDO

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'OBSPACK_MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100 ) 'Turn on ObsPack diagnostic? : ',
     &                    Input_Opt%Do_ObsPack
         WRITE( 6, 100 ) 'Suppress logfile output?    : ',
     &                    Input_Opt%ObsPack_Quiet
         WRITE( 6, 110 ) 'ObsPack input file          : ',
     &                    TRIM( Input_Opt%ObsPack_InputFile  )
         WRITE( 6, 110 ) 'ObsPack output file         : ',
     &                    TRIM( Input_Opt%ObsPack_OutputFile )
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )

      END SUBROUTINE READ_OBSPACK_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: N, S

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND48_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND48_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &             'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND48 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND48

      ! Timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND48_FILE

      ! Frequency
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND48_FREQ

      ! Number of stations 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_nd48_menu:4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND48_N_STA
      
      ! Read individual stations
      DO S = 1, Input_Opt%ND48_N_STA
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'read_nd48_menu:5', RC )
         READ( SUBSTRS(1), * ) Input_Opt%ND48_IARR(S) 
         READ( SUBSTRS(2), * ) Input_Opt%ND48_JARR(S)
         READ( SUBSTRS(3), * ) Input_Opt%ND48_LARR(S) 
         READ( SUBSTRS(4), * ) Input_Opt%ND48_NARR(S) 
      ENDDO

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND48 STATION TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '----------------------------'
         WRITE( 6, 100 ) 'Turn on ND48 timeseries?    : ',
     &                    Input_Opt%DO_ND48
         WRITE( 6, 110 ) 'ND48 timeseries file name   : ',
     &                    TRIM( Input_Opt%ND48_FILE )
         WRITE( 6, 120 ) 'ND48 save frequency [sec]   : ',
     &                    Input_Opt%ND48_FREQ

         DO S = 1, Input_Opt%ND48_N_STA
            WRITE( 6, 130 ) S, Input_Opt%ND48_IARR(S),
     &                         Input_Opt%ND48_JARR(S),
     &                         Input_Opt%ND48_LARR(S),
     &                         Input_Opt%ND48_NARR(S) 
         ENDDO
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, I5    )                   
 130  FORMAT( 'ND48 timeseries station', i4, ' : ', 4i5 )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: N
 
      ! Strings
      CHARACTER(LEN=255)  :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255)  :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND49_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND49_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg= 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND49

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND49_FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND49, -1,
     &                     'read_nd49_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%N_ND49
         READ( SUBSTRS(N), * ) Input_Opt%ND49_TRACERS(N)
      ENDDO

      ! FREQ
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND49_FREQ

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:5', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND49_IMIN, Input_Opt%ND49_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:6', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND49_JMIN, Input_Opt%ND49_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd49_menu:7', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND49_LMIN, Input_Opt%ND49_LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd49_menu:8', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

#if defined( MASSCONS )
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%  MASS CONSERVATION TEST:                                  %%%
      !%%%  Always turn on ND49 for the mass conservation test.      %%%
      !%%%  Also specify a single tracer for ND49 (in this case,     %%%
      !%%%  it will be total CO2). (ewl, 6/24/15)                    %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

      ! Reset quantities: Turn on ND49 but with 1 tracer
      Input_Opt%DO_ND49         = .TRUE.
      Input_Opt%N_ND49          = 1
      Input_Opt%ND49_TRACERS    = 0
      Input_Opt%ND49_TRACERS(1) = 1

      ! Print info to stdout stating that we have changed some quantities
      WRITE( 6, '(a)' ) REPEAT( '%', 79 )
      WRITE( 6, 200   ) 
      WRITE( 6, 205   )
      WRITE( 6, 210   )
      WRITE( 6, 215   )
      WRITE( 6, '(a)' ) REPEAT( '%', 79 )
 200  FORMAT( '%%% MASS CONSERVATION TEST'                         )
 205  FORMAT( '%%% Automatically reset these ND49 MENU settings: ' )
 210  FORMAT( '%%% DO_ND49 is now TRUE'                            )
 215  FORMAT( '%%% TRACERS is now 1 (i.e. Total CO2 only)'         )
#endif

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND49 3-D INSTANTANEOUS TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '--------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND49 timeseries?    : ',
     &                    Input_Opt%DO_ND49
         WRITE( 6, 110 ) 'ND49 timeseries file name   : ',
     &                    TRIM( Input_Opt%ND49_FILE )
         WRITE( 6, 120 ) 'ND49 timeseries tracers     : ', 
     &                    ( Input_Opt%ND49_TRACERS(N), N=1,
     &                      Input_Opt%N_ND49 )
         WRITE( 6, 130 ) 'ND49 save frequency [sec]   : ',
     &                    Input_Opt%ND49_FREQ
         WRITE( 6, 130 ) 'ND49 longitude limits       : ',
     &                    Input_Opt%ND49_IMIN, Input_Opt%ND49_IMAX
         WRITE( 6, 130 ) 'ND49 latitude  limits       : ',
     &                    Input_Opt%ND49_JMIN, Input_Opt%ND49_JMAX
         WRITE( 6, 130 ) 'ND49 level     limits       : ',
     &                    Input_Opt%ND49_LMIN, Input_Opt%ND49_LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I4 )
 130  FORMAT( A, 2I5   )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND50_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND50_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND49 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND50

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND50_FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LND50_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      Input_Opt%LND50_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND50, -1,
     &                     'read_nd50_menu:4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%N_ND50
         READ( SUBSTRS(N), * ) Input_Opt%ND50_TRACERS(N)
      ENDDO

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:5', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND50_IMIN, Input_Opt%ND50_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:6', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND50_JMIN, Input_Opt%ND50_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd50_menu:7', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND50_LMIN, Input_Opt%ND50_LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd50_menu:8', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND50 3-D 24hr AVG TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '---------------------------------'
         WRITE( 6, 100 ) 'Turn on ND50 timeseries?    : ',
     &                    Input_Opt%DO_ND50
         WRITE( 6, 110 ) 'ND50 timeseries file name   : ',
     &                    TRIM( Input_Opt%ND50_FILE )
         WRITE( 6, 100 ) 'Output as HDF?              : ',
     &                    Input_Opt%LND50_HDF
         WRITE( 6, 120 ) 'ND50 timeseries tracers     : ', 
     &                    ( Input_Opt%ND50_TRACERS(N), N=1,
     &                      Input_Opt%N_ND50 )
         WRITE( 6, 130 ) 'ND50 longitude limits       : ',
     &                    Input_Opt%ND50_IMIN, Input_Opt%ND50_IMAX
         WRITE( 6, 130 ) 'ND50 latitude  limits       : ',
     &                    Input_Opt%ND50_JMIN, Input_Opt%ND50_JMAX
         WRITE( 6, 130 ) 'ND50 level     limits       : ',
     &                    Input_Opt%ND50_LMIN, Input_Opt%ND50_LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I4 )
 130  FORMAT( A, 2I5   )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND51_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND51_Menu (in module GeosCore/input_mod.F)'
      
      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &           'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND51 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51_FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LND51_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      Input_Opt%LND51_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51, -1,
     &                     'read_nd51_menu:4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%N_ND51
         READ( SUBSTRS(N), * ) Input_Opt%ND51_TRACERS(N)
      ENDDO

      ! NHMS_WRITE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:6', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR_WRITE

      ! Make sure ND51_HR_WRITE is in the range 0-23.999 hrs
      Input_Opt%ND51_HR_WRITE = MOD( Input_Opt%ND51_HR_WRITE, 24e+0_fp )

      ! HR1, HR2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:7', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51_HR1, Input_Opt%ND51_HR2

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:8', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:9', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51_menu:10', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51_menu:11', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND51 MORNING OR AFTERNOON TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '-----------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND51 timeseries?    : ',
     &                    Input_Opt%DO_ND51
         WRITE( 6, 110 ) 'ND51 timeseries file name   : ',
     &                    TRIM( Input_Opt%ND51_FILE )
         WRITE( 6, 100 ) 'Output as HDF?              : ',
     &                    Input_Opt%LND51_HDF
         WRITE( 6, 120 ) 'ND41 timeseries tracers     : ', 
     &                    ( Input_Opt%ND51_TRACERS(N), N=1,
     &                      Input_Opt%N_ND51 )
         WRITE( 6, 140 ) 'ND51 hour to write to disk  : ',
     &                    Input_Opt%ND51_HR_WRITE
         WRITE( 6, 140 ) 'ND51 averaging period [GMT] : ',
     &                    Input_Opt%ND51_HR1,  Input_Opt%ND51_HR2
         WRITE( 6, 130 ) 'ND51 longitude limits       : ',
     &                    Input_Opt%ND51_IMIN, Input_Opt%ND51_IMAX
         WRITE( 6, 130 ) 'ND51 latitude  limits       : ',
     &                    Input_Opt%ND51_JMIN, Input_Opt%ND51_JMAX
         WRITE( 6, 130 ) 'ND51 altitude  limits       : ',
     &                    Input_Opt%ND51_LMIN, Input_Opt%ND51_LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I4 )
 130  FORMAT( A, 2I5   )
 140  FORMAT( A, 2F5.1 )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER              :: N

      ! Strings
      CHARACTER(LEN=255)   :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255)   :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND51b_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND51b_MENU  (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND51b diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND51b

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND51b_FILE

      ! Output as hdf
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LND51b_HDF

#if   !defined( USE_HDF5 )
      ! Reset HDF5 option to FALSE if we are not linking to HDF library
      ! (bmy, 12/21/09)
      Input_Opt%LND51b_HDF = .FALSE.
#endif

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND51b, -1,
     &                     'read_nd51b_menu:4',RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%N_ND51b
         READ( SUBSTRS(N), * ) Input_Opt%ND51b_TRACERS(N)
      ENDDO

      ! NHMS_WRITE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:5', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR_WRITE

      ! Make sure ND51b_HR_WRITE is in the range 0-23.999 hrs
      Input_Opt%ND51b_HR_WRITE = MOD(Input_Opt%ND51b_HR_WRITE,24e+0_fp)

      ! HR1, HR2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:6', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_HR1, Input_Opt%ND51b_HR2

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:7', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:8', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX

      ! LMIN, LMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd51b_menu:9', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd51b_menu:10', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' )'ND51b MORNING OR AFTERNOON TIMESERIES MENU'
         WRITE( 6, '(  a)' )'-----------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND51b timeseries?    : ',
     &                    Input_Opt%DO_ND51
         WRITE( 6, 110 ) 'ND51b timeseries file name   : ',
     &                    TRIM( Input_Opt%ND51b_FILE )
         WRITE( 6, 100 ) 'Output as HDF?               : ',
     &                    Input_Opt%LND51b_HDF
         WRITE( 6, 120 ) 'ND41 timeseries tracers      : ',
     &                    ( Input_Opt%ND51b_TRACERS(N), N=1,
     &                      Input_Opt%N_ND51b )
         WRITE( 6, 140 ) 'ND51b hour to write to disk  : ',
     &                    Input_Opt%ND51b_HR_WRITE
         WRITE( 6, 140 ) 'ND51b averaging period [GMT] : ',
     &                    Input_Opt%ND51b_HR1,  Input_Opt%ND51b_HR2
         WRITE( 6, 130 ) 'ND51b longitude limits       : ',
     &                    Input_Opt%ND51b_IMIN, Input_Opt%ND51b_IMAX
         WRITE( 6, 130 ) 'ND51b latitude  limits       : ',
     &                    Input_Opt%ND51b_JMIN, Input_Opt%ND51b_JMAX
         WRITE( 6, 130 ) 'ND51b altitude  limits       : ',
     &                    Input_Opt%ND51b_LMIN, Input_Opt%ND51b_LMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )
 140  FORMAT( A, 2F5.1 )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  20 Sep 2016 - R. Yantosca - Fixed typo (N_ND63 instead of ND63) when
!                              saving to Input_Opt%ND63_TRACERS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: N

      ! Strings
      CHARACTER(LEN=255)  :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255)  :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ND63_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_ND63_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &           'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Turn on ND63 diagnostic
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_ND63
      IF ( Input_Opt%DO_ND63 .and. .not. Input_Opt%LEMIS ) THEN
         WRITE( 6, '(a)' ) 'WARNING: Emissions are turned off. ND63'
         WRITE( 6, '(a)' ) 'will also be turned off.'
         Input_Opt%DO_ND63 = .FALSE.
      ENDIF

      ! Instantaneous 3-D timeseries file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%ND63_FILE

      ! Tracers to include
      CALL SPLIT_ONE_LINE( SUBSTRS, Input_Opt%N_ND63, -1,
     &                     'read_nd63_menu:3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      DO N = 1, Input_Opt%N_ND63
         READ( SUBSTRS(N), * ) Input_Opt%ND63_TRACERS(N)
      ENDDO

      ! FREQ
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:4', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND63_FREQ

      ! IMIN, IMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd63_menu:5', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND63_IMIN, Input_Opt%ND63_IMAX

      ! JMIN, JMAX
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2,  'read_nd63_menu:6', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND63_JMIN, Input_Opt%ND63_JMAX

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'read_nd63_menu:8', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Return success
      RC = GC_SUCCESS

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'ND63 SHIP TIMESERIES MENU'
         WRITE( 6, '(  a)' ) '--------------------------------------'
         WRITE( 6, 100 ) 'Turn on ND63 timeseries?    : ',
     &                    Input_Opt%DO_ND63
         WRITE( 6, 110 ) 'ND63 timeseries file name   : ', 
     &                    TRIM( Input_Opt%ND63_FILE )
         WRITE( 6, 120 ) 'ND63 timeseries tracers     : ', 
     &                    ( Input_Opt%ND63_TRACERS(N), N=1,
     &                      Input_Opt%N_ND63 )
         WRITE( 6, 130 ) 'ND63 save frequency [sec]   : ', 
     &                    Input_Opt%ND63_FREQ
         WRITE( 6, 130 ) 'ND63 longitude limits       : ', 
     &                    Input_Opt%ND63_IMIN, Input_Opt%ND63_IMAX
         WRITE( 6, 130 ) 'ND63 latitude  limits       : ', 
     &                    Input_Opt%ND63_JMIN, Input_Opt%ND63_JMAX
      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5    )
 110  FORMAT( A, A     )
 120  FORMAT( A, 100I3 )
 130  FORMAT( A, 2I5   )

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
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
#endif
      USE ErrCode_Mod
      USE gckpp_Parameters,   ONLY : NFAM
      USE gckpp_Monitor,      ONLY : FAM_NAMES
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  14 Jan 2016 - C. Keller   - Now refer to Input_Opt%LKPP instead of LKPP
!  27 Mar 2017 - R. Yantosca - For tagCO, ND65 slots < nAdvect are L(CO)
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( BPCH_DIAG )
      INTEGER            :: F,        N,          N_ADVECT

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
#endif

      !=================================================================
      ! READ_PROD_LOSS_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS

#if defined( BPCH_DIAG )
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Prod_Loss_Menu (in module GeosCore/input_mod.F)'

      ! Error check
      IF ( CT1 /= 2 ) THEN 
         ErrMsg = 'SIMULATION MENU & ADVECTED SPECIES MENU ' //
     &            'must be read in first!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Read info about prod & loss families
      !=================================================================

      ! Turn on production & loss diagnostic (e.g. ND65 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DO_SAVE_PL', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_SAVE_PL

      ! Read number of levels for ND65 diagnostic 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'ND65', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%ND65

      ! Copy field to variable in CMN_DIAG
      ND65 = Input_Opt%ND65

      ! Save P(O3) & L(O3) for tagged O3 run? (i.e. ND20 diagnostic)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'DO_SAVE_O3', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_SAVE_O3

      ! Save P(CO) from CH4 & NMVOC for tagged CO run?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,  'DO_SAVE_PCO', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%DO_SAVE_PCO

      !=================================================================
      ! Error check families for certain types of simulations
      !=================================================================

      ! Offline aerosol -- turn off DO_SAVE_PL, since we use ND05,
      ! ND06, ND07, ND08, ND13 etc diagnostics instead of ND65
      IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN 
         Input_Opt%DO_SAVE_PL    = .FALSE.
         Input_Opt%DO_SAVE_O3    = .FALSE.
         Input_Opt%DO_SAVE_PCO   = .FALSE.
         Input_Opt%ND65          = 0
      ENDIF

      !=================================================================
      ! Set fields of Input Options object
      !=================================================================

      ! Number of advected species
      N_ADVECT = Input_Opt%N_ADVECT 

      IF ( Input_Opt%DO_SAVE_PL ) THEN
         IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
            ! Fullchem - Obtain NFAM from KPP
            Input_Opt%NFAM = NFAM
         ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
            ! Tagged O3
            Input_Opt%NFAM = 2*N_ADVECT
         ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
            ! Tagged CO
            IF ( Input_Opt%LPCO_NMVOC ) THEN
               Input_Opt%NFAM = N_ADVECT+2
            ELSE
               Input_Opt%NFAM = N_ADVECT+6
            ENDIF
         ENDIF
      ENDIF

      ! Return if there are no prod/loss families
      ! or if we have turned off this diagnostic
      IF ( .not. ( Input_Opt%DO_SAVE_PL .and. Input_Opt%NFAM > 0 )) THEN
         Input_Opt%DO_SAVE_PL         = .FALSE.
         Input_Opt%DO_SAVE_O3         = .FALSE.
         Input_Opt%DO_SAVE_PCO        = .FALSE.
         Input_Opt%ND65               = 0
      ENDIF

      ! Loop over families
      DO F = 1, Input_Opt%NFAM

         IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

            ! Fullchem - Obtain FAM_NAME from KPP
            Input_Opt%FAM_NAME(F) = FAM_NAMES(F) 

         ELSEIF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

            ! Tagged O3
            IF ( F <= N_ADVECT ) THEN
               Input_Opt%FAM_NAME(F) =
     &            'P' // TRIM(Input_Opt%AdvectSpc_Name(F))
            ELSE
               Input_Opt%FAM_NAME(F) =
     &            'L' // TRIM(Input_Opt%AdvectSpc_Name(F-N_ADVECT))
            ENDIF

         ELSEIF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

            ! Tagged CO
            IF ( F <= N_ADVECT ) THEN
               Input_Opt%FAM_NAME(F) = 'L'//Input_Opt%AdvectSpc_Name(F)
            ELSEIF ( F == N_ADVECT+1 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_CH4'
            ELSEIF ( F == N_ADVECT+2 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_NMVOC'
            ELSEIF ( F == N_ADVECT+3 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_ISOP'
            ELSEIF ( F == N_ADVECT+4 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_CH3OH'
            ELSEIF ( F == N_ADVECT+5 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_MONO'
            ELSEIF ( F == N_ADVECT+6 ) THEN
               Input_Opt%FAM_NAME(F) = 'PCO_ACET'
            ENDIF

         ENDIF

         ! Get family type as prod or loss
         IF ( Input_Opt%FAM_NAME(F)(1:1) == 'P'   .or.
     &        Input_Opt%FAM_NAME(F)(1:1) == 'p' ) THEN
            Input_Opt%FAM_TYPE(F) = 'prod'
         ELSE
            Input_Opt%FAM_TYPE(F) = 'loss'
         ENDIF

      ENDDO

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'PROD & LOSS DIAGNOSTIC MENU'
         WRITE( 6, '(  a)' ) '---------------------------'      
         WRITE( 6, 100 ) 'Turn on prod & loss diag?   : ',
     &                    Input_Opt%DO_SAVE_PL
         WRITE( 6, 110 ) '# of levels for P/L diag    : ',
     &                    Input_Opt%ND65
         WRITE( 6, 100 ) 'Save P(O3), L(O3) for TagO3?: ',
     &                    Input_Opt%DO_SAVE_O3
         WRITE( 6, 100 ) 'Save CO Production?: ',
     &                    Input_Opt%DO_SAVE_PCO
      
         ! Loop over families
         DO F = 1, Input_Opt%NFAM

            ! Write family name and type
            WRITE( 6, 120 ) TRIM(Input_Opt%FAM_NAME(F)),
     &                      TRIM(Input_Opt%FAM_TYPE(F))

         ENDDO

      ENDIF

      ! FORMAT statements
 100  FORMAT( A, L5 )
 110  FORMAT( A, I5 )
 120  FORMAT( /, 'Family=', A10, '  Type=', A4 )
#endif

      END SUBROUTINE READ_PROD_LOSS_MENU
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg,      ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Nested_Grid_Menu (in module GeosCore/input_mod.F)'

      ! Save out TPCORE BC's at 4x5
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO

      ! Save out TPCORE BC's at 2x2.5?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO2x25', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO2x25

      ! Save out TPCORE BC's over NA
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO_NA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO_NA

      ! Directory where NA TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TPBC_DIR_NA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%TPBC_DIR_NA

      ! Save out TPCORE BC's over EU
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO_EU', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO_EU

      ! Directory where EU TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TPBC_DIR_EU', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%TPBC_DIR_EU

      ! Save out TPCORE BC's over CH
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' LWINDO_CH', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO_CH

      ! Directory where CH TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TPBC_DIR_CH', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%TPBC_DIR_CH

      ! Save out TPCORE BC's over AS (0.5x0.625 only)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO_AS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO_AS

      ! Directory where AS TPCORE BC's are stored (0.5x0.625 only)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TPBC_DIR_AS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%TPBC_DIR_AS

      ! Save out TPCORE BC's over Custom Region
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LWINDO_CU', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LWINDO_CU

      ! Directory where TPCORE BC's are stored
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TPBC_DIR', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%TPBC_DIR

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NESTED_TS

      ! Lower left box
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'I1, J1', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NESTED_I1, Input_Opt%NESTED_J1

      ! Timestep for 4x5 TPCORE BC's
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'I2, J2', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NESTED_I2, Input_Opt%NESTED_J2

      ! Nested grid transport offsets 
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 4, 'I0W, J0W, I0E, J0E', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%NESTED_I0W,Input_Opt%NESTED_J0W,
     &                        Input_Opt%NESTED_I0E,Input_Opt%NESTED_J0E

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Error checks
      !=================================================================
#if defined( MERRA2 )
      IF ( Input_Opt%LWINDO_CH ) THEN
         ErrMsg = 'MERRA-2 is not set up for the nested CH grid, use '//
     &            'nested AS instead.'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'NESTED GRID MENU'
         WRITE( 6, '(  a)' ) '----------------'
         WRITE( 6, 100 ) 'Save TPCORE BC''s at 4x5     : ',
     &                    Input_Opt%LWINDO
         WRITE( 6, 100 ) 'Save TPCORE BC''s at 2x25    : ',
     &                    Input_Opt%LWINDO2x25
         WRITE( 6, 100 ) 'Save BC''s over North America: ',
     &                    Input_Opt%LWINDO_NA
         WRITE( 6, 110 ) 'Directory with NA BC''s      : ', 
     &                    TRIM( Input_Opt%TPBC_DIR_NA )
         WRITE( 6, 100 ) 'Save BC''s over Europe       : ',
     &                    Input_Opt%LWINDO_EU
         WRITE( 6, 110 ) 'Directory with EU BC''s      : ', 
     &                    TRIM( Input_Opt%TPBC_DIR_EU )
         WRITE( 6, 100 ) 'Save BC''s over China/SE Asia: ',
     &                    Input_Opt%LWINDO_CH
         WRITE( 6, 110 ) 'Directory with CH BC''s      : ', 
     &                    TRIM( Input_Opt%TPBC_DIR_CH )
         WRITE( 6, 100 ) 'Save BC''s over Asia         : ',
     &                    Input_Opt%LWINDO_AS
         WRITE( 6, 110 ) 'Directory with AS BC''s      : ', 
     &                    TRIM( Input_Opt%TPBC_DIR_AS )
         WRITE( 6, 100 ) 'Save BC''s over custom region: ',
     &                    Input_Opt%LWINDO_CU
         WRITE( 6, 110 ) 'Directory with custom BC''s  : ', 
     &                    TRIM( Input_Opt%TPBC_DIR )
         WRITE( 6, 120 ) 'Timestep for BC''s [sec]     : ',
     &                    Input_Opt%NESTED_TS
         WRITE( 6, 120 ) 'LL box of BC region         : ',
     &                    Input_Opt%NESTED_I1,  Input_Opt%NESTED_J1
         WRITE( 6, 120 ) 'UR box of BC region         : ',
     &                    Input_Opt%NESTED_I2,  Input_Opt%NESTED_J2
         WRITE( 6, 130 ) 'Nested grid window offsets  : ',
     &                    Input_Opt%NESTED_I0W, Input_Opt%NESTED_J0W,
     &                    Input_Opt%NESTED_I0E, Input_Opt%NESTED_J0E
      END IF

      ! FORMAT statements
 100  FORMAT( A, L5  )
 110  FORMAT( A, A   )
 120  FORMAT( A, 2I5 )
 130  FORMAT( A, 4I5 )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N 

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_NESTED_GRID_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Benchmark_Menu (in module GeosCore/input_mod.F)'


      ! Save benchmark diagnostic output?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LSTDRUN', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LSTDRUN

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'INITIAL_FILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%STDRUN_INIT_FILE

      ! Filename for initial tracer mass
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'FINAL_FILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), '(a)' ) Input_Opt%STDRUN_FINAL_FILE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'BENCHMARK MENU'
         WRITE( 6, '(  a)' ) '--------------'
         WRITE( 6, 100     ) 'Save benchmark diag output?  : ',
     &                        Input_Opt%LSTDRUN
         WRITE( 6, 110     ) 'File for initial species mass: ',  
     &                        TRIM( Input_Opt%STDRUN_INIT_FILE )
         WRITE( 6, 110     ) 'File for final species mass  : ',  
     &                        TRIM( Input_Opt%STDRUN_FINAL_FILE )
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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!  25 Apr 2016 - R. Yantosca - Now init mercury modules from GIGC_INIT_EXTRA
!  27 Jul 2016 - M. Sulprizio- Remove IN_HG_RST_FILE and OUT_HG_RST_FILE.
!                              Hg restart fields are now stored in the GEOS-Chem
!                              netCDF restart file.
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LDYNOCEAN,      LPREINDHG
      LOGICAL            :: LGTMM,          USE_CHECKS
      LOGICAL            :: LARCTICRIV,     LKRedUV
      INTEGER            :: N 
      CHARACTER(LEN=255) :: GTMM_RST_FILE

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg,         ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_MERCURY_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_Mercury_Menu (in module GeosCore/input_mod.F)'

      ! Use error check for tag/tot Hg?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'USE_CHECKS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%USE_CHECKS

      ! Use dynamic ocean model?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LDYNOCEAN', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LDYNOCEAN

      ! Use preindustrial simulation?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LPREINDHG', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LPREINDHG

      ! Use GTMM ?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LGTMM', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LGTMM

      ! Name of GTMM restart file
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GTMM_RST_FILE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%GTMM_RST_FILE

      ! Use Arctic river Hg?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LARCTICRIV', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LARCTICRIV

      ! Tie reducible HgII(aq) to UV-B?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'LKREDUV', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%LKRedUV

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Check on logical
      IF (.NOT.( Input_Opt%ITS_A_MERCURY_SIM ) ) THEN
         Input_Opt%LGTMM      = .FALSE.
         Input_Opt%LDYNOCEAN  = .FALSE.
         Input_Opt%LARCTICRIV = .FALSE.
         Input_Opt%LKRedUV    = .FALSE.
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'MERCURY MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 110 ) 'Error check tag & total Hg? : ', 
     &                    Input_Opt%USE_CHECKS
         WRITE( 6, 110 ) 'Use dynamic ocean Hg model? : ', 
     &                    Input_Opt%LDYNOCEAN
         WRITE( 6, 110 ) 'Preindustrial simulation?   : ', 
     &                    Input_Opt%LPREINDHG
         WRITE( 6, 110 ) 'Use GTMM ?                  : ', 
     &                    Input_Opt%LGTMM
         WRITE( 6, 120 ) '=> GTMM restart file        : ',
     &                    TRIM( Input_Opt%GTMM_RST_FILE )
         WRITE( 6, 110 ) 'Use Arctic river Hg ?       : ', 
     &                    Input_Opt%LARCTICRIV
         WRITE( 6, 110 ) 'Tie HgII(aq) red. to UV-B?  : ', 
     &                    Input_Opt%LKRedUV
      END IF

      ! FORMAT statements
 100  FORMAT( A, I4  )
 110  FORMAT( A, L5  )
 120  FORMAT( A, A   )

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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  20 Sep 2016 - R. Yantosca - Rewrote IF statement to avoid Gfortran error
!  16 Jun 2017 - M. Sulprizio- Remove emission options since they are now
!                              handled by HEMCO
!  16 Jun 2017 - M. Sulprizio- Add option for GOSAT CH4 observation operator
!  17 Aug 2017 - M. Sulprizio- Add option for TCCON CH4 observation operator
!  18 Oct 2017 - M. Sulprizio- Remove CH4 budget option
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_CH4_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_CH4_Menu (in module GeosCore/input_mod.F)'
                                                
      ! Use GOSAT CH4 observation operator?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'GOSAT_CH4_OBS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%GOSAT_CH4_OBS

      ! Use TCCON CH4 observation operator?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'TCCON_CH4_OBS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%TCCON_CH4_OBS
                                 
      ! Separator line   
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF  

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'CH4 MENU'
         WRITE( 6, '(  a)' ) '-----------'
         WRITE( 6, 100     ) 'Use GOSAT obs operator: ',
     &                        Input_Opt%GOSAT_CH4_OBS
         WRITE( 6, 100     ) 'Use TCCON obs operator: ',
     &                        Input_Opt%TCCON_CH4_OBS
      END IF

      ! FORMAT statements
 100  FORMAT( A, L5  )

      ! Return success
      RC = GC_SUCCESS         

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
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  01 Oct 2012 - C. Friedman - Initial version
!  26 Nov 2012 - M. Payer    - Added ProTeX headers
!  29 Nov 2012 - M. Payer    - Now pass Input_Opt, RC as arguments
!  26 Mar 2013 - R. Yantosca - Now pass Input_Opt to INIT_POPS
!  14 Apr 2014 - R. Yantosca - Now echo all POPs parameters to log file
!  23 Jun 2014 - R. Yantosca - Now use Input_Opt%POP_EMISDIR
!  25 Aug 2014 - M. Sulprizio- Move call to INIT_POPs to GIGC_Init_Extra
!  27 Aug 2014 - M. Sulprizio- Remove POP_EMISDIR, emissions are now handled
!                              by HEMCO
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_POPS_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     &   ' -> at Read_POPs_Menu (in module GeosCore/input_mod.F)'

      ! POP species
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_TYPE', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_TYPE

      ! Dummy for future process logical switches
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'CHEM_PROCESS', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%CHEM_PROCESS

      ! Molecular weight
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_XMW', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_XMW

      ! KOA
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KOA', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_KOA

      ! KBC
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_KBC', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_KBC

      ! OH oxidation
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPG_OH', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPG_OH

      ! O3 oxidation 1
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, ' POP_K_POPP_O3A', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3A

      ! O3 oxidation 2
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_K_POPP_O3B', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_K_POPP_O3B

      ! H*
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_HSTAR', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_HSTAR

      ! DEL_H
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_H', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_H

      ! DEL_Hw
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'POP_DEL_Hw', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      READ( SUBSTRS(1:N), * ) Input_Opt%POP_DEL_Hw

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'separator', RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      IF( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'POPS MENU'
         WRITE( 6, '(  a)' ) '------------'
         WRITE( 6, 120     ) 'Species of POP        : ',
     &                        Input_Opt%POP_TYPE
         WRITE( 6, 110     ) 'Chemistry on?         : ',
     &                        Input_Opt%CHEM_PROCESS
         WRITE( 6, 130     ) 'POP_XMW               : ',
     &                        Input_Opt%POP_XMW
         WRITE( 6, 130     ) 'POP_KOA               : ',
     &                        Input_Opt%POP_KOA
         WRITE( 6, 130     ) 'POP_KBC               : ',
     &                        Input_Opt%POP_KBC
         WRITE( 6, 130     ) 'POP_K_POPG_OH         : ',
     &                        Input_Opt%POP_K_POPG_OH
         WRITE( 6, 130     ) 'POP_K_POPP_O3A        : ',
     &                        Input_Opt%POP_K_POPP_O3A
         WRITE( 6, 130     ) 'POP_K_POPP_O3B        : ',
     &                        Input_Opt%POP_K_POPP_O3B
         WRITE( 6, 130     ) 'POP_HSTAR             : ',
     &                        Input_Opt%POP_HSTAR
         WRITE( 6, 130     ) 'POP_DEL_H             : ',
     &                        Input_Opt%POP_DEL_H
         WRITE( 6, 130     ) 'POP_DEL_Hw            : ',
     &                        Input_Opt%POP_DEL_Hw 
      ENDIF

      ! FORMAT statements
 110  FORMAT( A, L5  )
 120  FORMAT( A, A   )
 130  FORMAT( A, ES10.2 )

      ! Return success
      RC = GC_SUCCESS  

      END SUBROUTINE READ_POPS_MENU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_passive_species_menu
!
! !DESCRIPTION: Subroutine READ\_PASSIVE\_SPECIES\_MENU reads the passive 
!  species menu section of the GEOS-Chem input file; this defines passive
!  species to be used for this simulation. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_PASSIVE_SPECIES_MENU( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller   - Initial version
!  14 Jul 2017 - E. Lundgren - Store passive species values in Input_Opt and
!                              remove dependency on passive_species_mod.F90
!  02 Aug 2017 - R. Yantosca - Bug fix: Now use NP as the loop index over
!                              passive species so as to not conflict with
!                              the returned # of substrings
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!  04 Apr 2018 - E. Lundgren - Call routine to init Input_Opt passive spc arrays
!  21 Dec 2018 - M. Sulprizio- Update passive species menu to include long name
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, P, D

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      ! Arrays
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_PASSIVE_SPECIES_MENU begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error reading the "input.geos" file!'
      ThisLoc = 
     & ' -> at Read_Passive_Species_Menu (in mod. GeosCore/input_mod.F)'

      ! Initialize
      P = 0
      D = 0

      IF ( am_I_Root ) THEN
         WRITE( 6, '(/,a)' ) 'PASSIVE SPECIES MENU'
         WRITE( 6, '(  a)' ) '---------------------'
      ENDIF

      ! Do until exit
      DO

         ! Read passive species information for each passive species
         ! Every passive species line is expected to have 4 entries:
         ! - Species name
         ! - Species molecular weight 
         ! - Atmospheric lifetime (s)
         ! - Initial atmospheric concentration (v/v)
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'Passive Species', RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
     
         ! Exit here if separator line is encountered
         IF ( INDEX ( TRIM(SUBSTRS(1)), '-----' ) > 0 ) EXIT

         ! Make sure there are at least 4 entries
         IF ( N < 4 ) THEN
            ErrMsg = 'Each passive species is expected to have ' //
     &             'at least four entries: Name, MW, TAU, initial ' //
     &             'concentration, and long name (optional).'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Increase number of passive species by one
         P = P + 1

         ! Stop simulation and print warning if we exceed maximum number
         ! of passive species hardcoded in input_opt_mod.F90.
         IF ( P > Input_Opt%Max_PassiveSpc ) THEN
            ErrMsg = 'Number of passive species exceeds maximum. ' //
     &             'This value can be modified in input_opt_mod.F90.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! Read and store species information
         Input_Opt%PASSIVE_ID(P)   = P
         Input_Opt%PASSIVE_NAME(P) = TRIM( SUBSTRS(1) )
         READ( SUBSTRS(2), * ) Input_Opt%PASSIVE_MW(P)
         READ( SUBSTRS(3), * ) Input_Opt%PASSIVE_TAU(P)
         READ( SUBSTRS(4), * ) Input_Opt%PASSIVE_INITCONC(P)

         ! Check if optional entry for long name is included
         IF ( N > 4 ) THEN
            READ( SUBSTRS(5), * ) Input_Opt%PASSIVE_LONGNAME(P)
         ENDIF

         ! Determine if the passive species decays (i.e. if it has 
         ! an atmospheric lifetime that is not -1).  This will allow
         ! us to skip those passive species that do not decay in
         ! routine CHEM_PASSIVE_SPECIES, to speed up execution.
         IF ( Input_Opt%PASSIVE_TAU(P) > 0.0_fp ) THEN
            D                            = D + 1
            Input_Opt%NPASSIVE_DECAY     = D
            Input_Opt%PASSIVE_DECAYID(D) = P
         ENDIF


         ! Verbose
         IF ( am_I_Root ) THEN
            WRITE( 6, '(a)' ) 'Added passive species: '
            WRITE( 6, 110   ) ' - Species name                : ', 
     &            TRIM( Input_Opt%PASSIVE_NAME(P) ) 
            WRITE( 6, 120   ) ' - Molec. weight [g/mol]       : ', 
     &            Input_Opt%PASSIVE_MW(P)
            WRITE( 6, 130   ) ' - Lifetime [s]                : ', 
     &            Input_Opt%PASSIVE_TAU(P)
            WRITE( 6, 130   ) ' - Initial concentration [v/v] : ', 
     &            Input_Opt%PASSIVE_INITCONC(P)
            WRITE( 6, 110   ) ' - Species long name           : ',
     &            TRIM( Input_Opt%PASSIVE_LONGNAME(P) )
         ENDIF

      ENDDO

      IF ( P < 0 ) THEN
         ErrMsg = 'Cannot add passive species ' //
     &             TRIM(Input_Opt%PASSIVE_NAME(P) ) //
     &            ': # of passive species is smaller than 1!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN 
      ENDIF

      ! Total number of passive species
      Input_Opt%NPASSIVE = P

 110  FORMAT( A, A )
 120  FORMAT( A, F10.2  )
 130  FORMAT( A, ES10.2 )

      ! Return success
      RC = GC_SUCCESS

      END SUBROUTINE READ_PASSIVE_SPECIES_MENU
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
      SUBROUTINE VALIDATE_DIRECTORIES( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE GC_GRID_MOD,   ONLY : ITS_A_NESTED_GRID 
      USE TIME_MOD,      ONLY : EXPAND_DATE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  20 Jun 2014 - R. Yantosca - Now use fields from Input_Opt
!  05 Mar 2015 - R. Yantosca - Now also check CHEM_INPUTS_DIR
!  12 Aug 2015 - R. Yantosca - Now validate Input_Opt%MERRA2_DIR
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Dir

      !=================================================================
      ! VALIDATE_DIRECTORIES begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Invalid directory encountered!'
      ThisLoc = 
     &   ' -> at Validate_Directories (in module GeosCore/input_mod.F)'

      ! Check directories
      CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%DATA_DIR, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%CHEM_INPUTS_DIR, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%RUN_DIR, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Only validate the TPCORE BC directories if we need it
      IF ( ITS_A_NESTED_GRID() ) THEN
#if   defined( NESTED_NA)
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_NA, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#elif defined( NESTED_EU)
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_EU, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#elif defined( NESTED_CH)
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_CH, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#elif defined( NESTED_AS)
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_AS, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#endif
      ENDIF

      IF ( Input_Opt%LWINDO_CU ) THEN
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      IF ( Input_Opt%LWINDO_NA ) THEN
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_NA, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      IF ( Input_Opt%LWINDO_EU ) THEN
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_EU, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      IF ( Input_Opt%LWINDO_CH ) THEN
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_CH, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      IF ( Input_Opt%LWINDO_AS ) THEN
         CALL CHECK_DIRECTORY( Input_Opt, Input_Opt%TPBC_DIR_AS, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

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
      SUBROUTINE CHECK_DIRECTORY( Input_Opt, DIR, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE FILE_MOD,      ONLY : File_Exists
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(OptInput),   INTENT(INOUT) :: Input_Opt   ! Input Options object
      CHARACTER(LEN=*), INTENT(INOUT) :: DIR         ! Dir to be checked
!
! !OUTPUT PARAMETERS: 
!
      INTEGER,          INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  20 Jun 2014 - R. Yantosca - Now added Input_Opt object via arg list
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: C

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! CHECK_DIRECTORY begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Check_Directory (in module GeosCore/input_mod.F)'


      ! Locate the last non-white-space character of NEWDIR
      C = LEN_TRIM( DIR )

      ! Add the trailing directory separator if it is not present
      IF ( DIR(C:C) /= '/' ) THEN 
         DIR(C+1:C+1) = TRIM( '/' )
      ENDIF
     
      !=================================================================
      ! Test if the directory actually exists
      !=================================================================

      ! If the directory does not exist then stop w/ an error message
      IF ( .not. FILE_EXISTS( DIR ) ) THEN 
         ErrMsg = 'Invalid directory: ' // TRIM( DIR ) 
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
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
      SUBROUTINE CHECK_TIME_STEPS( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,      ONLY : SET_TIMESTEPS
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
!  23 Jun 2014 - R. Yantosca - Removed references to logical_mod.F
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: LCONV, LCHEM,       LDRYD 
      LOGICAL            :: LEMIS, LTRAN,       LTURB 
      INTEGER            :: I,     J,           K 
      INTEGER            :: L,     TS_SMALLEST, TS_DIAG
      INTEGER            :: TS_CHEM, TS_EMIS, TS_CONV, TS_DYN
      INTEGER            :: TS_UNIT, TS_RAD,  MAX_DYN

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! CHECK_TIME_STEPS begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Check_Time_Steps (in module GeosCore/input_mod.F)'

      ! Copy fields from Input_Opt
      LCONV = Input_Opt%LCONV
      LCHEM = Input_Opt%LCHEM
      LDRYD = Input_Opt%LDRYD 
      LEMIS = Input_Opt%LEMIS
      LTRAN = Input_Opt%LTRAN
      LTURB = Input_Opt%LTURB 
     
      TS_CHEM = Input_Opt%TS_CHEM
      TS_EMIS = Input_Opt%TS_EMIS
      TS_CONV = Input_Opt%TS_CONV
      TS_DYN  = Input_Opt%TS_DYN 
      TS_RAD  = Input_Opt%TS_RAD
 
      ! NUNIT is time step in minutes for unit conversion
      TS_UNIT = -1

      ! Define maximum timestep for transport
#if   defined( GRID4x5   ) 
      MAX_DYN = 1800
#elif defined( GRID2x25  )
      MAX_DYN = 900
#elif defined( GRID05x0625 )
      MAX_DYN = 600 
#elif defined( GRID025x03125 )
      MAX_DYN = 300
#else
      MAX_DYN = 3600
#endif

      ! If TS_DYN is greater than MAX_DYN, then stop w/ error
      IF ( .not. Input_Opt%HPC ) THEN
         IF ( Input_Opt%TS_DYN > MAX_DYN .and. LTRAN ) THEN
            WRITE( ErrMsg, 300 ) 'Transport timestep exceeds max:',
     &                        Input_Opt%TS_DYN, MAX_DYN
 300        FORMAT( a, i8, ' >', i8 ) 
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
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

      ! SDE 2017-02-24: Always use LTRAN on the assumption that it will
      ! be used as a "heartbeat". This ensures that chemistry always
      ! takes place at the same time, regardless of whether or not
      ! transport is enabled.
      !IF ( .not. LTRAN                  ) I = 999999 
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
         ErrMsg = 'The transport time step should be the smallest one'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
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

      ! Check if all time steps are multiples of the smallest.
      ! (ccc, 5/13/09)
      IF ( L /= -999999 .and. MOD( TS_CHEM, TS_SMALLEST ) /= 0 ) THEN
         WRITE( ErrMsg, 100 ) 'Chemistry', TS_CHEM, TS_SMALLEST
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      
      IF ( K /= -999999 .and. MOD( TS_EMIS, TS_SMALLEST ) /= 0 ) THEN
         WRITE( ErrMSg, 100 ) 'Emission', TS_EMIS, TS_SMALLEST
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( J /= -999999 .and. MOD( TS_CONV, TS_SMALLEST ) /= 0 ) THEN
         WRITE( ErrMsg, 100 ) 'Convection', TS_CONV, TS_SMALLEST
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( I /= -999999 .and. MOD( TS_DYN, TS_SMALLEST ) /= 0 ) THEN
         WRITE( ErrMsg, 100 ) 'Transport', TS_DYN, TS_SMALLEST
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Initialize timesteps in "time_mod.f"
      CALL SET_TIMESTEPS( am_I_Root,
     &                    CHEMISTRY  = TS_CHEM, EMISSION  = TS_EMIS, 
     &                    DYNAMICS   = TS_DYN,  UNIT_CONV = TS_UNIT,
     &                    CONVECTION = TS_CONV, DIAGNOS   = TS_DIAG,
     &                    RADIATION  = TS_RAD )

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
      SUBROUTINE IS_LAST_DAY_GOOD( Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE JULDAY_MOD,    ONLY : JULDAY
      USE TIME_MOD,      ONLY : GET_NYMDe, ITS_A_LEAPYEAR, YMD_EXTRACT
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt  ! Input options
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Moved to "input_mod.f" from "main.f" (bmy, 1/11/05)
!  (2 ) Now call ITS_A_LEAPYEAR with FORCE=.TRUE. to always return whether
!        the year Y would be a leap year, regardless of met field type.
!        (swu, bmy, 4/24/06)
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  08 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: IS_LEAPYEAR
      INTEGER            :: NYMDe, Y, M, D, LASTDAY
      REAL(fp)           :: JD, JD0

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! IS_LAST_DAY_GOOD begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Is_Last_Day_Good (in module GeosCore/input_mod.F)'

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

      ! Exit w/ error if THIS_NJDAY = 0 
      IF ( Input_Opt%NJDAY(LASTDAY) == 0 ) THEN
         ErrMsg = 'No output scheduled on last day of run!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
     
      END SUBROUTINE IS_LAST_DAY_GOOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_error_checks
!
! !DESCRIPTION: Makes sure that certain species are defined in order to
!  proceed with a certain option.  Halts the simulation with an error message
!  if incorrect inputs  would have caused  a simulation to crash.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Do_Error_Checks( am_I_root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE State_Chm_Mod, ONLY : Ind_
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt
!
! !OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(OUT)   :: RC
!
! !REMARKS:
!  These error checks were originally called when the various menus were
!  read in from disk.  However, in order to use the Ind_() function to look
!  up given species indices, we need to call this after the Species Database
!  (which is in State_Chm) is initialized.  Therefore, we have now moved these
!  error checks to this routine, which is now called from GC_Init_Extra.
!                                                                             .
!  The Ind_() function now defines all species ID's.  It returns -1 if
!  a species cannot be found.  The prior behavior was to return 0 if a 
!  species wasn't found.  Therefore, in order to preserve the logic of the
!  error checks, we must force any -1's returned by Ind_() to 0's in
!  this subroutine.
!
! !REVISION HISTORY:
!  22 Jun 2016 - R. Yantosca - Initial version
!  08 Nov 2017 - R. Yantosca - Now return error condition to calling routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: I

      ! Strings
      CHARACTER(LEN=255) :: MSG, LOCATION

      !=================================================================
      ! Initialization
      !=================================================================

      ! Assume success
      RC       = GC_SUCCESS

      ! Define location string
      LOCATION = '-> at Do_Error_Checks (in GeosCore/input_mod.F)'

      !=================================================================
      ! Error check SEASALT AEROSOLS
      !=================================================================
      I = MAX( Ind_('SALA','A'), 0 ) 
     &  + MAX( Ind_('SALC','A'), 0 )

      IF ( Input_Opt%LSSALT ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check MARINE ORGANIC AEROSOLS
      !=================================================================
      I = MAX( Ind_('MOPO','A'), 0 )  
     &  + MAX( Ind_('MOPI','A'), 0 )

      IF ( Input_Opt%LMPOA ) THEN
         IF ( .not. Input_Opt%LSSALT ) THEN
            MSG = 'LMPOA=T but LSSALT=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
         IF ( I == 0 ) THEN
            MSG = 'LMPOA=T but MARINE ORGANIC AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use MARINE ORGANIC AEROSOLS if LMPOA=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check SULFUR AEROSOLS
      !=================================================================
      I = MAX( Ind_('DMS' ,'A'), 0 )
     &  + MAX( Ind_('SO2' ,'A'), 0 )
     &  + MAX( Ind_('SO4' ,'A'), 0 )
     &  + MAX( Ind_('SO4s','A'), 0 )
     &  + MAX( Ind_('MSA' ,'A'), 0 )
     &  + MAX( Ind_('NH3' ,'A'), 0 )
     &  + MAX( Ind_('NH4' ,'A'), 0 )
     &  + MAX( Ind_('NITs','A'), 0 ) 

      IF ( Input_Opt%LSULF ) THEN

         ! We now compute the production of SO4s and NITs, so when 
         ! LSULF=T, then we must also have LSSALT=T (bec, bmy, 4/13/05)
         IF ( .not. Input_Opt%LSSALT ) THEN 
            MSG = 'LSULF=T now also requires LSSALT=T!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF

         ! Stop w/ error if everything is undefined
         IF ( I == 0 ) THEN
            MSG = 'LSULF=T but ONLINE SULFUR AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF

      ELSE

         ! If LSULF=F but we have defined species, stop w/ error
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SULFUR AEROSOLS if LSULF=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF

      ENDIF

      !=================================================================
      ! Error check CARBON AEROSOLS
      !=================================================================

      ! SOAupdate: Add POA (hotp 10/11/09)
      I = MAX( Ind_('BCPO','A'), 0 )
     &  + MAX( Ind_('BCPI','A'), 0 )
     &  + MAX( Ind_('OCPO','A'), 0 )
     &  + MAX( Ind_('OCPI','A'), 0 )
     &  + MAX( Ind_('POA1','A'), 0 )

      IF ( Input_Opt%LCARB ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LCARB=T but ONLINE CARBON AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE CARBON AEROSOLS if LCARB=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      IF ( Input_Opt%LSVPOA .and. ( .NOT. Input_Opt%LSOA ) ) THEN
         MSG = 'Semivolatile POA requires COMPLEX SOA (LSOA=T)'
         CALL GC_Error( Msg, RC, Location )
         RETURN
      ENDIF

      ! SOAupdate: Error check (hotp 8/24/09)
      ! OCPI and OCPO are the non-volatile POA species
      ! POA (along w/ POG, OPOA, and OPOG) are the semivol POA species
      ! You can't have both!
      I = MAX( Ind_('OCPI','A'), 0 )  
     &  + MAX( Ind_('OCPO','A'), 0 )

      IF ( Ind_('POA1') > 0 ) THEN
         IF ( I > 0 ) THEN
            MSG = 'Semivolatile POA species is defined in addition to
     &             Nonvolatile POA'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
         IF ( ( .NOT. Input_Opt%LSOA   ) .or. 
     &        ( .NOT. Input_Opt%LSVPOA ) ) THEN
            MSG = 'Semivolatile POA requires LSOA=T and LSVPOA=T'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      ! SOAupdate
      ! Options for organic aerosol species:
      ! IF LSOA = F: only OCPI and OCPO
      ! IF LSOA = T:
      !   OCPI OCPO SOA (non-vol + original traditional)
      !   POA POG OPOA OPOG SOA BTX NAP (semivol + orig trad + IVOC )
      ! NAP emissions are set in input.geos
      ! LSVPOA is just a check (doesn't do anything hotp 7/21/10)
      I = MAX( Ind_('POA1' ,'A'), 0 )
     &  + MAX( Ind_('POA2' ,'A'), 0 )
     &  + MAX( Ind_('POG1' ,'A'), 0 )
     &  + MAX( Ind_('POG2' ,'A'), 0 )
     &  + MAX( Ind_('OPOA1','A'), 0 )
     &  + MAX( Ind_('OPOA2','A'), 0 )
     &  + MAX( Ind_('OPOG1','A'), 0 )
     &  + MAX( Ind_('OPOG2','A'), 0 )

      IF ( Input_Opt%LSVPOA ) THEN
         IF ( I < 8 ) THEN
            MSG = 'Not enough semivolatile POA species!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
         IF ( Ind_('NAP','A') < 0 ) THEN
            MSG = 'Semivolatile POA requires IVOCs/NAP!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check SECONDARY ORGANIC AEROSOLS
      !=================================================================

      ! Check for complex SOA species
      I = MAX( Ind_('TSOA1','A'), 0 )
     &  + MAX( Ind_('TSOA2','A'), 0 )
     &  + MAX( Ind_('TSOA3','A'), 0 )
     &  + MAX( Ind_('ISOA1','A'), 0 )
     &  + MAX( Ind_('ISOA2','A'), 0 )
     &  + MAX( Ind_('ISOA3','A'), 0 )
     &  + MAX( Ind_('ASOA1','A'), 0 )
     &  + MAX( Ind_('ASOA2','A'), 0 )
     &  + MAX( Ind_('ASOA3','A'), 0 )
     &  + MAX( Ind_('ASOAN','A'), 0 )

      IF ( Input_Opt%LSOA ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSOA=T but COMPLEX SOA species are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use COMPLEX SOA species if LSOA=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check DUST AEROSOLS
      !=================================================================

! JKODROS - Uncomment this if the dust fixes are back working. 
#if defined( TOMAS )
      ! For TOMAS only: If DUST1 is present, the other dust species are too
      I = MAX( Ind_('DUST1','A'), 0 )
#else
      ! Non-TOMAS simulations: Need all DST1-DST4 species
      I = MAX( Ind_('DST1','A'), 0 )
     &  + MAX( Ind_('DST2','A'), 0 )
     &  + MAX( Ind_('DST3','A'), 0 )
     &  + MAX( Ind_('DST4','A'), 0 )
#endif

      IF ( Input_Opt%LDUST ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LDUST=T but ONLINE DUST AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE DUST AEROSOLS if LDUST=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check DUST NITRATE    AEROSOLS
      !             DUST SULFATE    AEROSOLS
      !             DUST ALKALINITY AEROSOLS
      !=================================================================
      I = MAX( Ind_('NITd1'  ,'A'), 0 ) 
     &  + MAX( Ind_('NITd2'  ,'A'), 0 ) 
     &  + MAX( Ind_('NITd3'  ,'A'), 0 )
     &  + MAX( Ind_('NITd4'  ,'A'), 0 ) 
     &  + MAX( Ind_('SO4d1'  ,'A'), 0 ) 
     &  + MAX( Ind_('SO4d2'  ,'A'), 0 )
     &  + MAX( Ind_('SO4d3'  ,'A'), 0 ) 
     &  + MAX( Ind_('SO4d4'  ,'A'), 0 ) 
     &  + MAX( Ind_('DSTAL1' ,'A'), 0 )
     &  + MAX( Ind_('DSTAL2' ,'A'), 0 ) 
     &  + MAX( Ind_('DSTAL3' ,'A'), 0 ) 
     &  + MAX( Ind_('DSTAL4' ,'A'), 0 )

      IF ( Input_Opt%LDSTUP ) THEN
         IF ( I < 12 ) THEN
            MSG = 'LDSTUP=T but COATED DUST AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use COATED DUST AEROSOLS if LDSTUP=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Error check SEASALT AEROSOLS
      !=================================================================
      I = MAX( Ind_('SALA','A'), 0 ) 
     &  + MAX( Ind_('SALC','A'), 0 )

      IF ( Input_Opt%LSSALT ) THEN
         IF ( I == 0 ) THEN
            MSG = 'LSSALT=T but ONLINE SEASALT AEROSOLS are undefined!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ELSE
         IF ( I > 0 ) THEN
            MSG = 'Cannot use ONLINE SEASALT AEROSOLS if LSSALT=F!'
            CALL GC_Error( Msg, RC, Location )
            RETURN
         ENDIF
      ENDIF

      
      !=================================================================
      ! Error check CFC emissions
      !=================================================================
      IF ( Input_Opt%LCFCEMIS ) THEN
         IF ( Ind_('CFC113') < 0 .or.
     &        Ind_('CFC114') < 0 .or. 
     &        Ind_('CFC115') < 0) THEN
            Input_Opt%LCFCEMIS = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'CFCs are emitted but CFC species are undefined.'
            WRITE( 6, '(a)' ) 'Turning off CFC emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      !=================================================================
      ! Error check N2O emissions
      !=================================================================
      IF ( Input_Opt%LN2OEMIS ) THEN
         IF ( Ind_('N2O') < 0 ) THEN
            Input_Opt%LN2OEMIS = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'N2O is emitted but N2O species is undefined.'
            WRITE( 6, '(a)' ) 'Turning off nitrous oxide emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      !=================================================================
      ! Error check stratospheric H2O
      !=================================================================
      IF ( Input_Opt%LUCX .and. Input_Opt%LSETH2O ) THEN
         IF (Ind_('H2O') < 0 ) THEN
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'H2O is set but H2O species is undefined.'
            Input_Opt%LSETH2O = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ELSE
         Input_Opt%LSETH2O = .FALSE.
      ENDIF

      !=================================================================
      ! Error check halogen emissions
      !=================================================================
      IF ( Input_Opt%LCLEMIS ) THEN
         IF ( Ind_('Cl') < 0 ) THEN
            Input_Opt%LCLEMIS = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'Cl species emitted but Cl species is undefined.'
            WRITE( 6, '(a)' ) 'Turning off chlorine emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      IF ( Input_Opt%LBREMIS ) THEN
         IF ( Ind_('Br2') < 0 ) THEN
            Input_Opt%LBREMIS     = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'Br species emitted but Br2 species is undefined.'
            WRITE( 6, '(a)' ) 'Turning off bromine emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      !=================================================================
      ! Error check OCS emissions
      !=================================================================
      IF ( Input_Opt%LOCSEMIS ) THEN
         IF ( Ind_('OCS') < 0 ) THEN
            Input_Opt%LOCSEMIS = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'OCS is emitted but OCS species is undefined.'
            WRITE( 6, '(a)' ) 'Turning off OCS emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      !=================================================================
      ! Error check methane emissions
      !=================================================================
      IF ( Input_Opt%LCH4EMIS .AND. Input_Opt%LCH4SBC ) THEN
         MSG =  
     &  ' CH4 surface BC & emissions cannot both be turned on!'
         CALL GC_Error( Msg, RC, Location )
         RETURN
      ENDIF

      IF ( Input_Opt%LCH4EMIS   .or. 
     &     Input_Opt%LCH4SBC  ) THEN
         IF ( Ind_('CH4') < 0 ) THEN
            Input_Opt%LCH4EMIS = .FALSE.
            Input_Opt%LCH4SBC  = .FALSE.
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
            WRITE( 6, '(/,a,/)' ) 'Warning in input_mod.F: '
     &            // 'CH4 is emitted but CH4 species is undefined.'
            WRITE( 6, '(a)' ) 'Turning off methane emissions.'
            WRITE( 6, '(a)'     ) REPEAT( '=', 79 )
         ENDIF
      ENDIF

      END SUBROUTINE DO_ERROR_CHECKS
!EOC
      END MODULE INPUT_MOD
