!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: chemdr
!
! !DESCRIPTION: Subroutine CHEMDR is the driver subroutine for full chemistry 
!  with SMVGEAR or KPP.
!\\
!\\}
! !INTERFACE:
!
      SUBROUTINE CHEMDR( am_I_Root, Input_Opt, 
     &                   State_Met, State_Chm, RC  )
!
! !USES:
!
      USE AEROSOL_MOD,          ONLY : AEROSOL_CONC
      USE AEROSOL_MOD,          ONLY : RDAER
      USE AEROSOL_MOD,          ONLY : SOILDUST
      USE CMN_DIAG_MOD
      USE CMN_NOX_MOD
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE COMODE_MOD,           ONLY : ABSHUM
      USE COMODE_MOD,           ONLY : CSPEC
      USE COMODE_MOD,           ONLY : ERADIUS
      USE COMODE_MOD,           ONLY : JLOP
      USE COMODE_MOD,           ONLY : JLOP_PREVIOUS
      USE COMODE_MOD,           ONLY : TAREA
      USE COMODE_LOOP_MOD
      USE DAO_MOD,              ONLY : MAKE_AVGW
      USE DIAG_OH_MOD,          ONLY : DO_DIAG_OH
      USE DIAG_PL_MOD,          ONLY : DO_DIAG_PL
      USE DUST_MOD,             ONLY : RDUST_ONLINE
      USE DUST_MOD,             ONLY : RDUST_OFFLINE
      USE ERROR_MOD,            ONLY : DEBUG_MSG
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GIGC_State_Chm_Mod,   ONLY : ChmState
      USE GIGC_State_Met_Mod,   ONLY : MetState
      USE PLANEFLIGHT_MOD,      ONLY : SETUP_PLANEFLIGHT
      USE RESTART_MOD,          ONLY : READ_CSPEC_FILE 
      USE TIME_MOD,             ONLY : GET_MONTH
      USE TIME_MOD,             ONLY : GET_NHMS
      USE TIME_MOD,             ONLY : GET_NYMD
      USE TIME_MOD,             ONLY : GET_YEAR
      USE TIME_MOD,             ONLY : ITS_A_NEW_DAY
      USE TRACERID_MOD,         ONLY : IDTNO
      USE TRACERID_MOD,         ONLY : IDTO3
      USE TRACERID_MOD,         ONLY : SETTRACE
      USE TROPOPAUSE_MOD,       ONLY : SAVE_FULL_TROP

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Important input variables from "dao_mod.f" and "uvalbedo_mod.f":
!  ============================================================================
!  ALBD        : DAO visible albedo                         [unitless]
!  AVGW        : Mixing ratio of water vapor                [v/v] 
!  OPTD        : DAO grid-box optical depths (for FAST-J)   [unitless]
!  SUNCOS      : Cosine of solar zenith angle               [unitless]
!  SUNCOSB     : Cosine of solar zenith angle 1 hr from now [unitless]
!  UVALBEDO    : TOMS UV albedo 340-380 nm (for FAST-J)     [unitless]
!                                                                             .
!  Important input variables from "comode.h" or "comode_mod.f":
!  ============================================================================
!  NPTS        : Number of points (grid-boxes) to calculate
!  REMIS       : Emission rates                             [molec/cm3/s-1]
!  RAERSOL     : Frequency of gas-aerosol collisions        [s-1]
!  PRESS       : Pressure                                   [Pa]
!  TMPK        : Temperature                                [K]
!  ABSHUM      : Absolute humidity                          [molec/cm3]
!  CSPEC       : Initial species concentrations             [molec/cm3]
!                                                                             .
!  Important output variables in "comode.h" etc.
!  ============================================================================
!  NAMESPEC    : Character array of species names
!  NNSPEC      : # of ACTIVE + INACTIVE (not DEAD) species
!  CSPEC       : Final species concentrations               [molec/cm3]
!                                                                             .
!  Other Important Variables
!  ============================================================================
!  MAXPTS      : Maximum number of points or grid-boxes (in "comsol.h")
!                (NPTS must be <= MAXPTS, for SLOW-J only)
!  MAXDEP      : Maximum number of deposition species (note # of
!                depositing species listed in tracer.dat must be <= MAXDEP)
!  IGAS        : Maximum number of gases, ACTIVE + INACTIVE
!  IO93        : I/O unit for output for "ctm.chem" file
!                                                                             .
!  Input files for SMVGEAR II:
!  ============================================================================
!   mglob.dat  : control switches                       (read in "reader.f")
!  tracer.dat  : list of tracers, emitting species      (read in "reader.f")
!                and depositing species
! globchem.dat : species list, reaction list,           (read in "chemset.f")
!                photolysis reaction list
!                                                                             .
!  Input files for FAST-J photolysis:
!  ============================================================================
!     ratj.d   : Lists photo species, branching ratios  (read in "rd_js.f")
! jv_atms.dat  : Climatology of T and O3                (read in "rd_prof.f")
! jv_spec.dat  : Cross-sections for each species        (read in "RD_TJPL.f")
!
! !REVISION HISTORY:
!  15 Nov 2001 - Adapted from original code by lwh, jyl, gmg, djj. 
!  (1 ) Cleaned up a lot of stuff.  SUNCOS, OPTD, ALBD, and AVGW are now 
!        referenced from dao_mod.f.  IREF and JREF are obsolete.  Also 
!        updated comments. (bmy, 9/27/01)
!  (2 ) Do not declare LPRT or set LPRT = .FALSE. in "chemdr.f".  LPRT is 
!        included via "CMN" and is defined in "main.f". (bmy, 10/9/01)
!  (3 ) Removed obsolete data from 9/01 (bmy, 10/23/01)
!  (4 ) ERADIUS(JLOOP) is now ERADIUS(JLOOP,1) and TAREA(JLOOP) is now
!        TAREA(JLOOP,1) for sulfate aerosol.  Updated comments. (bmy, 11/15/01)
!  (5 ) Renamed routine PAFTOP to DEBUG_MSG.  Also deleted obsolete code
!        from 11/01.  Enhanced debug output via DEBUG_MSG.  Also reference
!        the UVALBEDO array directly from "uvalbedo_mod.f".  Remove UVALBEDO
!        from the argument list.  Removed obsolete comments. (bmy, 1/15/02)
!  (6 ) Now pass LPAUSE to "initgas.f" via the arg list (bmy, 2/14/02)
!  (7 ) Now call "rdaer.f" instead of RDAEROSOL to read the aerosol and dust 
!        fields from disk.  Also, ignore hygroscopic growth for dust.  Now
!        pass SAVEHO2 and FRACNO2 arrays in the arg list to "ohsave.f"; these 
!        return HO2 conc.'s and NO2 fraction.  Delete NTRACE from call
!        to "ohsave.f", it's obsolete.  Delete reference to DARSFCA from
!        "comode_mod.f", it's obsolete. (rvm, bmy, 2/27/02)
!  (8 ) Removed obsolete code from 2/02. (bmy, 4/15/02)
!  (9 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (10) Now reference IU_CTMCHEM from "file_mod.f".  Assign the value of
!        IU_CTMCHEM (which =93) to IO93 for SMVGEAR routines.  Also open 
!        "ctm.chem" file on the first call as file unit # IO93.  Add
!        references to "CMN_DIAG" and "planeflight_mod.f".  Call routine
!        SETUP_PLANEFLIGHT to initialize the plane track diagnostic
!        after reading the "chem.dat" file.  (bmy, 7/2/02)
!  (11) Now reference AD, T and BXHEIGHT from "dao_mod.f".  Also removed 
!        obsolete commented out code in various sections below.  Now also
!        references "tracerid_mod.f".  Also remove reference to BIOTRCE, since
!        this is now obsolete.  Now make FIRSTCHEM a local SAVED variable
!        instead of an argument.  Now calls MAKE_AVGW, which was formerly
!        called in "main.f". (bmy, 11/15/02)
!  (12) Now reference "AIRVOL" from "dao_mod.f".  Now declare local array
!        SO4_NH4_NIT, which will contain lumped SO4, NH3, NIT aerosol.  Now
!        pass SO4_NH4_NIT to "rdaer.f" via the arg list if sulfate chemistry
!        is turned on.  Now also references CMN_SETUP. (rjp, bmy, 3/23/03)
!  (13) Removed ITAU from the arg list.  Removed reference to IHOUR.  Now use
!        functions GET_MONTH, GET_YEAR from "time_mod.f" (bmy, 3/27/03)
!  (14) Remove KYEAR and TWO_PI, these are now obsolete for SMVGEAR II.  Now 
!        open unit #93 and call READER in the same FIRSTCHEM if-block.  Now
!        Replace call to CHEMSET with call to READCHEM.  JPARSE is now called 
!        from w/in READCHEM.  Replace call to INITGAS w/ call to GASCONC.
!        Removed reference to "file_mod.f".  Remove call to SETPL, we now must
!        call this in "readchem.f" before the call to JSPARSE. 
!        (bdf, ljm, bmy, 5/8/03)
!  (15) Now reference routine GET_GLOBAL_CH4 from "global_ch4_mod.f".  Also
!        added CH4_YEAR as a SAVEd variable. (bnd, bmy, 7/1/03)
!  (16) Remove references to MONTHP, IMIN, ISEC; they are obsolete and not 
!        defined anywhere. (bmy, 7/16/03)
!  (17) Now reference SUNCOSB from "dao_mod.f".  Now pass SUNCOSB to "chem.f". 
!        Also remove LSAMERAD from call to CHEM, since it's obsolete. 
!        (gcc, bmy, 7/30/03)
!  (18) Added BCPO, BCPI, OCPO, OCPI, and SOILDUST arrays.  Now pass SOILDUST
!       to RDUST_ONLINE (in "dust_mod.f").  Now pass PIEC, POEC, PIOC, POOC to
!       "rdaer.f".  Now references "dust_mod.f". (rjp, tdf, bmy, 4/1/04)
!  (19) Added SALA and SALC arrays for passing seasalt to rdaer.f.  Now
!        rearranged the DO loop for computational efficiency. (bmy, 4/20/04)
!  (20) Added OCF parameter to account for the other chemical components that
!        are attached to OC.  Also now handle hydrophilic OC differently for
!        online & offline SOA. (rjp, bmy, 7/15/04)
!  (21) Now reference "logical_mod.f".  Now reference STT and N_TRACERS from
!        "tracer_mod.f".  Now references DO_DIAG_PL from "diag_pl_mod.f".
!        Now references DO_DIAG_OH from "diag_oh_mod.f".  Now references
!        AEROSOL_CONC, RDAER, & SOILDUST from "aerosol_mod.f" (bmy, 7/20/04)
!  (22) Now references ITS_A_NEW_DAY from "time_mod.f".  Now calls routine
!        SETUP_PLANEFLIGHT at the start of each new day. (bmy, 3/24/05)
!  (23) FAST-J is now the default, so we don't need the LFASTJ C-preprocessor 
!        switch any more (bmy, 6/23/05)
!  (24) Now remove LPAUSE from the arg list to "ruralbox.f" and "gasconc.f".
!        (bmy, 8/22/05)
!  (25) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (26) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (27) Remove more obsolete SLOW-J code references.  Also now move function
!        calls from subroutine "chem.f" into "chemdr.f".  Remove obsolete
!        arguments from call to RURALBOX. (bmy, 4/10/06) 
!  (28) Remove reference to "global_ch4_mod.f".  Add error check for LISOPOH
!        when using the online SOA tracers. (dkh, bmy, 6/1/06)
!  (29) Now support variable tropopause (bdf, phs, bmy, 10/3/06)
!  (30) Now get CH4 concentrations for FUTURE_YEAR when using the future
!        emissions scale factors (swu, havala, bmy, 1/28/04)
!  (31) Now call "save_full_trop" at the end to account for "do_diag_pl" 
!        resetting some of CSPEC elements (phs, 6/3/08)
!  (32) Reading the CSPEC_FULL restart file if asked.(dkh, hotp, ccc 2/26/09)
!  (33) Added optional call to gckpp_driver (phs,ks,dhk, 09/15/09)
!  (34) CSPEC_FOR_KPP not used anymore (use CSPEC instead) (ccc, 12/3/09)
!  (35) Move the KPP interface in physproc.f to save memory (ccc, 12/3/09)
!  (36) Now bracket out dust emissions if TOMAS is invoked (bmy, 1/25/10)
!  (37) Now remove obsolete embedded chemistry stuff.  Modify arg list to
!        RURALBOX accordingly.   Removed obsolete LEMBED switch. (bmy, 2/26/10)
!  (38) Remove obsolete SUNCOSB (bmy, 4/28/10)
!  26 Aug 2010 - R. Yantosca - Added ProTeX headers
!  05 Oct 2011 - R. Yantosca - Now use SUNCOS30 (cosine of SZA at 30m after
!                              the GMT hour) instead of SUNCOS.  This fixes
!                              the photolysis issue.
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!   6 Mar 2012 - R. Yantosca - Now reference TO3 from dao_mod.F
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed Chemistry State object to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed Meteorology State object to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  24 Oct 2012 - R. Yantosca - Add fixes for connecting to GEOS-5 GCM via ESMF
!  25 Oct 2012 - R. Yantosca - Document #ifdefs for the GIGC implementation
!  25 Oct 2012 - R. Yantosca - Now add RC (success/failure) as an argument
!  08 Nov 2012 - R. Yantosca - Now add Input_Opt (Input Options) as an argument
!  11 Dec 2012 - R. Yantosca - Remove NI, NJ, NL arguments, they are not used
!  14 Mar 2013 - M. Payer    - Replace NOx and Ox with NO and O3 as part
!                              of removal of NOx-Ox partitioning
!  26 Feb 2013 - M. Long     - Now pass Input_Opt to INPHOT as an argument
!  19 Mar 2013 - R. Yantosca - Now copy Input_Opt%XNUMOL(1:N_TRACERS)
!  31 May 2013 - R. Yantosca - Now pass State_Met to DO_DIAG_PL, for TOMAS
!  06 Jun 2013 - M. Payer    - Move SETEMIS to after FAST_J (C. Holmes)
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER       :: I, J, JLOOP, L, NPTS, N, MONTH, YEAR
      INTEGER       :: WAVELENGTH
      INTEGER       :: DATE  

      ! SAVEd variables
      LOGICAL, SAVE :: FIRSTCHEM   = .TRUE.
      INTEGER, SAVE :: CH4_YEAR  = -1
      INTEGER, SAVE :: DATE_PREV = -1

      ! To use CSPEC_FULL restart (dkh, 02/12/09) 
      LOGICAL       :: IT_EXISTS
      REAL          :: START,FINISH

      ! For printing debug output
      LOGICAL       :: prtDebug

      ! Local variables for quantities from Input_Opt
      INTEGER       :: N_TRACERS
      REAL*8        :: XNUMOL(Input_Opt%N_TRACERS)
      LOGICAL       :: LCARB
      LOGICAL       :: LDUST
      LOGICAL       :: LPRT
      LOGICAL       :: LSSALT
      LOGICAL       :: LSULF  
      LOGICAL       :: LSOA
      LOGICAL       :: LVARTROP
      LOGICAL       :: LFUTURE
      LOGICAL       :: LSVCSPEC

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEMDR begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy logical fields from INPUT_OPT to local variables for use below
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)
      LCARB     = Input_Opt%LCARB
      LDUST     = Input_Opt%LDUST
      LPRT      = Input_Opt%LPRT
      LSSALT    = Input_Opt%LSSALT
      LSULF     = Input_Opt%LSULF
      LSOA      = Input_Opt%LSOA
      LVARTROP  = Input_Opt%LVARTROP
      LFUTURE   = Input_Opt%LFUTURE
      LSVCSPEC  = Input_Opt%LSVCSPEC

      ! Do we need to print debug output?
      prtDebug  = ( LPRT .and. am_I_Root )

      ! Set some size variables
      NLAT      = JJPAR
      NLONG     = IIPAR
      NVERT     = IVERT 
      NPVERT    = NVERT
      NPVERT    = NVERT + IPLUME

      ! Get month and year
      MONTH     = GET_MONTH()
      YEAR      = GET_YEAR()

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT       => State_Chm%Tracers

      !=================================================================
      ! Compute AVGW, the mixing ratio of water vapor
      !=================================================================
      CALL MAKE_AVGW( State_Met )

      !=================================================================
      ! Open "smv2.log" output file and read chem mechanism switches
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do nothing, since we will call READER etc. from the init 
      ! method of the ESMF interface (bmy, 10/24/12)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Call READER to read chemistry options from mglob.dat.  This is
      ! current practice in the std GEOS-Chem (bmy, 10/24/12)
      !-----------------------------------------------------------------
      IF ( FIRSTCHEM ) THEN

         ! Read from data file mglob.dat
         CALL READER( FIRSTCHEM, am_I_Root )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after READER' )
         ENDIF

         ! Set NCS for urban chemistry only (since that is where we
         ! have defined the GEOS-CHEM mechanism) (bdf, bmy, 4/21/03)
         NCS = NCSURBAN
      ENDIF
#endif

      !=================================================================      
      ! Call RURALBOX, which defines tropospheric boxes to be sent to
      ! the SMVGEAR solver, as well as setting up some SMVGEAR arrays.
      !=================================================================      

      ! Redefine NTLOOP since READER defines it initially (bmy, 9/28/04)
      NLOOP  = NLAT  * NLONG
      NTLOOP = NLOOP * NVERT

      !vartrop fix (dkh, 05/08/11)
      !Save a copy of JLOP from the previous time step before re-ordering
      !the grid mappings.
      JLOP_PREVIOUS(:,:,:) = JLOP(:,:,:)

      ! Locate the tropopsheric boxes
      CALL RURALBOX( State_Met )

      ! For the first chemistry timestep only, save JLOP into 
      ! JLOP_PREVIOUS, because there is no previous timestep.
      IF ( FIRSTCHEM ) THEN
         JLOP_PREVIOUS = JLOP
      ENDIF

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after RURALBOX' ) 
      ENDIF

      ! Reset NTTLOOP, the # of tropospheric grid boxes
      NTTLOOP = NTLOOP

      !=================================================================
      ! Call SETMODEL which defines number of grid-blocks in calculation,
      ! and copies meteorological parameters into local variables 
      !=================================================================
      CALL SETMODEL( am_I_Root )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after SETMODEL' )
      ENDIF

      !=================================================================
      ! Do the following only on the first call ...
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do nothing, since we will call these setup routines from the 
      ! init method of the ESMF interface (bmy, 10/24/12)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Call the following chemistry setup routines only on the very
      ! first chemistry timestep.  This is current practice in the
      ! std GEOS-Chem. (bmy, 10/24/12)
      !-----------------------------------------------------------------
      IF ( FIRSTCHEM ) THEN

         !---------------------------------
         ! Initialize chemistry mechanism
         !---------------------------------

         ! Read "globchem.dat" chemistry mechanism
         CALL READCHEM( am_I_Root, Input_Opt, RC )

         ! Set NCS=NCSURBAN here since we have defined our tropospheric
         ! chemistry mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
         NCS = NCSURBAN

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after READCHEM' )
         ENDIF

         !---------------------------------
         ! Check for LISOPOH for SOA
         !---------------------------------
         IF ( LSOA .and. ILISOPOH == 0 ) THEN
            CALL ERROR_STOP( 'LISOPOH needs to be defined for SOA!',
     &                       'chemdr.f' )
         ENDIF

         !---------------------------------------------------
         ! SOAupdate: Check ILISOPNO3 for SOA (hotp 6/1/10)
         !---------------------------------------------------
         IF ( LSOA .and. ILISOPNO3 == 0 ) THEN
            CALL ERROR_STOP( 'IDLISOPNO3 needs to be defined for SOA!',
     &                       'chemdr.f' )
         ENDIF

         !---------------------------------
         ! Set global concentration of CH4
         !---------------------------------
         IF ( ICH4 > 0 .and. ( CH4_YEAR /= GET_YEAR() ) ) THEN

            ! If CH4 is a SMVGEAR II species, then call GET_GLOBAL_CH4
            ! to return the globally-varying CH4 conc. as a function of
            ! year and latitude bin.  (ICH4 is defined in READCHEM.)
            ! (bnd, bmy, 7/1/03)
            !
            ! If we are using the future emissions, then get the CH4
            ! concentrations for FUTURE_YEAR.  Otherwise get the CH4
            ! concentrations for the current met field year. 
            ! (swu, havala, bmy, 1/24/08)
            IF ( LFUTURE ) THEN
               CH4_YEAR = GET_FUTURE_YEAR()
            ELSE
               CH4_YEAR = GET_YEAR()
            ENDIF

            ! Get CH4 [ppbv] in 4 latitude bins for each year
            CALL GET_GLOBAL_CH4( CH4_YEAR,  .TRUE.,  C3090S,
     &                           C0030S,    C0030N,  C3090N, 
     &                           am_I_Root, Input_Opt        )
         ENDIF

         !-------------------------------
         ! Initialize FAST-J photolysis
         !-------------------------------
         CALL INPHOT( LLTROP, NPHOT, Input_Opt, am_I_Root ) 
         
         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after INPHOT' )        
         ENDIF

         !-------------------------------
         ! Flag certain chemical species
         !-------------------------------
         CALL SETTRACE( am_I_Root, Input_Opt, State_Chm, RC )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after SETTRACE' )
         ENDIF

         !-------------------------------
         ! Flag emission & drydep rxns
         !-------------------------------
         CALL SETEMDEP( am_I_Root, Input_Opt, RC )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after SETEMDEP' )
         ENDIF

      ENDIF
#endif

      !=================================================================
      ! At the beginning of each new day, call SETUP_PLANEFLIGHT
      ! to see if there are any plane flight points to be processed
      !=================================================================

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do not call SETUP_PLANEFLIGHT, since we will not use this
      ! diagnostic with the ESMF interface. (bmy, 10/24/12)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! If this is the first chem timestep of a new day, then we need 
      ! to call SETUP_PLANEFLIGHT.  If chemistry is turned on, then 
      ! we need to place this call here, so as to make sure that the 
      ! chemical mechanism files (read by READER and READCHEM) have 
      ! been loaded.  This is current practice in the std GEOS-Chem.
      ! (bmy, 10/24/12)
      !-----------------------------------------------------------------

      ! Get todays' date
      DATE = GET_NYMD()

      IF ( ND40 .and. DATE /= DATE_PREV ) THEN
         CALL SETUP_PLANEFLIGHT( am_I_Root )
         DATE_PREV = DATE
      ENDIF

#endif

      !================================================================
      ! Get concentrations of aerosols in [kg/m3] 
      ! for FAST-J and optical depth diagnostics
      !=================================================================
      IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

         ! Skip this section if all these are turned off
         CALL AEROSOL_CONC
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ENDIF

      !=================================================================
      ! Call GASCONC which initializes gas concentrations and sets 
      ! miscellaneous parameters.  GASCONC also calls PARTITION, which
      ! splits up family tracers like NOx and Ox into individual
      ! chemical species for SMVGEAR.
      ! NOTE:
      ! (1) The call to GASCONC is modified to use CSPEC_FULL restart 
      !     file (dkh, hotp, ccc,2/26/09)
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Call GASCONC with FIRSTCHEM = .FALSE. and IS_EXISTS = .FALSE.
      ! Within GASCONC, we will call COPY_FULL_TROP to make sure that
      ! the CSPEC_FULL array gets copied to the CSPEC array.
      ! (bmy, 10/24/12)
      !-----------------------------------------------------------------
      CALL GASCONC(  am_I_Root, .FALSE.,   .FALSE.,   
     &               Input_Opt, State_Met, State_Chm, RC )
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! On the first chemistry timestep:
      !
      !   (1) If we are reading initial species values from the CSPEC
      !       restart file, then have GASCONC initialize the CSPEC 
      !       array with the values read in from the restart file.
      !
      !   (2) If we are not reading initial species values from the
      !       CSPEC restart file, then have GASCONC initialize the
      !       CSPEC array with the initial values listed in the
      !       "globchem.dat" chemistry mechanism file.
      !
      ! On all other chemistry timesteps:
      !
      !   (3) Call GASCONC with the current species concentrations in
      !       the CSPEC array.  Do not reset to initial values.
      !
      ! NOTE: In all 3 cases, GASCONC will call PARTITION to separate
      ! family tracers into individual species.
      !
      ! This is current practice in the std GEOS-Chem. (bmy, 10/24/12)
      !-----------------------------------------------------------------

      ! Assume there
      IT_EXISTS = .FALSE.

      ! Test if 
      IF ( FIRSTCHEM .and. LSVCSPEC ) THEN

         ! Read the CSPEC restart file.  If not found, then 
         ! return IT_EXISTS = .FALSE.
         CALL READ_CSPEC_FILE( am_I_Root,  Input_Opt, GET_NYMD(), 
     &                         GET_NHMS(), IT_EXISTS, State_Chm,  RC  )
   
         IF ( .not. IT_EXISTS ) THEN 

            !-----------------------------------------------------------
            ! (2) First chemistry timestep, Call GASCONC to initialize 
            !     CSPEC with default values from globchem.dat.
            !-----------------------------------------------------------
            IF ( am_I_Root ) THEN
               WRITE(6,*) 
     &  '    - CHEMDR: CSPEC restart not found, use background values'
            ENDIF

            CALL GASCONC( am_I_Root, FIRSTCHEM, IT_EXISTS, 
     &                    Input_Opt, State_Met, State_Chm, RC )
         ELSE 

            !-----------------------------------------------------------
            ! (1) First chemistry OK timestep, Call GASCONC to initialize 
            !     CSPEC with values from the CSPEC restart file
            !-----------------------------------------------------------
            IF ( am_I_Root ) THEN
               WRITE(6,*) 
     &        '    - CHEMDR: using CSPEC values from restart file'
            ENDIF

            CALL GASCONC( am_I_Root, .FALSE.,    IT_EXISTS, 
     &                    Input_Opt,  State_Met, State_Chm, RC )

         ENDIF 

      ELSE 

         !--------------------------------------------------------------
         ! (3) All other chemistry timesteps, call GASCONC with
         !     the current CSPEC values.  Do not re-initialize.
         !--------------------------------------------------------------
         CALL GASCONC(  am_I_Root, FIRSTCHEM, IT_EXISTS,
     &                  Input_Opt, State_Met, State_Chm, RC )

      ENDIF 

      IT_EXISTS = .FALSE.

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after READ_CSPEC_FILE' )
      ENDIF
#endif

      !=================================================================
      ! Call RDAER -- computes aerosol optical depths
      !=================================================================

      ! Call RDAER to compute AOD for FAST-J (skim, 02/03/11)
      WAVELENGTH = 0
      CALL RDAER( am_I_Root, Input_Opt, State_Met,  RC,
     &            MONTH,     YEAR,      WAVELENGTH      )

      !### Debug
      IF ( prtDebug ) THEN 
         CALL DEBUG_MSG( '### CHEMDR: after RDAER' )
      ENDIF

      !=================================================================
      ! If LDUST is turned on, then we have online dust aerosol in
      ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
      ! compute aerosol optical depth for FAST-J, etc.
      !
      ! If LDUST is turned off, then we do not have online dust aerosol
      ! in GEOS-CHEM...so read monthly-mean dust files from disk.
      ! (rjp, tdf, bmy, 4/1/04)
      !=================================================================
      IF ( LDUST ) THEN
         CALL RDUST_ONLINE(  SOILDUST,    WAVELENGTH, State_Met )
      ELSE

#if !defined( TOMAS )
            ! Don't read dust emissions from disk when using TOMAS,
            ! because TOMAS uses a different set of dust tracers than the 
            ! std code (win, bmy, 1/25/10)
            CALL RDUST_OFFLINE( MONTH, YEAR, WAVELENGTH, am_I_Root,
     &                          State_Met )
#endif
      ENDIF

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after RDUST' )
      ENDIF

      NPTS = NTTLOOP

      ! At present, we are only doing tropospheric chemistry, which 
      ! for the moment we are storing in SMVGEAR II's "urban" slot
      NCS = NCSURBAN

      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================
      CALL FAST_J( am_I_Root, State_Met, RC )              

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after FAST-J' )
      ENDIF

      !=================================================================  
      ! Call SETEMIS which sets emission rates REMIS 
      !=================================================================
      CALL SETEMIS( EMISRR,    EMISRRN,   am_I_Root, 
     &              Input_Opt, State_Met, RC         )
      
      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after SETEMIS' )
      ENDIf

      !================================================================
      ! Call chemistry routines
      !================================================================

      ! PHYSPROC calls both CALCRATE, which computes rxn rates 
      ! and SMVGEAR (if we do not use the solver coded by kpp), which
      ! is the chemistry solver
      CALL PHYSPROC( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after PHYSPROC' )
      ENDIF

      !=================================================================
      ! Call LUMP which lumps the species together after chemistry
      !=================================================================
      CALL LUMP( am_I_Root, Input_Opt, State_Chm, RC )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after LUMP' )
      ENDIF

      !=================================================================
      ! Call OHSAVE which saves info on OH AND HO2 concentrations 
      !=================================================================
      CALL OHSAVE( SAVEOH, SAVEHO2 )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after OHSAVE' )
      ENDIF

      !=================================================================
      ! Save quantities for computing mean OH lifetime
      !=================================================================
      CALL DO_DIAG_OH
      
      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after DO_DIAG_OH' )
      ENDIF

      !=================================================================
      ! Save production and loss for chemical families.  Also save
      ! P(Ox) and L(Ox) for a future tagged Ox run (if necessary).
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do not call DO_DIAG_PL, since we will not use this diagnostic
      ! with the ESMF interface. (bmy, 10/24/12)
      !-----------------------------------------------------------------
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Call DO_DIAG_PL to compute the family prod/loss diagnostics.
      ! This is current practice in the std GEOS-Chem. (bmy, 10/24/12)
      !-----------------------------------------------------------------
      IF ( Input_Opt%LFAMILY ) THEN
         CALL DO_DIAG_PL
     &        ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after DO_DIAG_PL' )
         ENDIF
      ENDIF
#endif
      
      !=================================================================
      ! Copy the modified values from CSPEC back to CSPEC_FULL
      !=================================================================

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! When we connect to the ESMF interface, the CSPEC_FULL array
      ! is defined from the values in the Chemistry State.  We need
      ! to copy these values to the CSPEC array for the rest of the
      ! chemistry routines to function properly.  Therefore, call
      ! COPY_FULL_TROP from tropopause_mod.F every time.
      !-----------------------------------------------------------------
      CALL SAVE_FULL_TROP( State_Chm )

      !###  Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### CHEMDR: after SAVE_FULL_TROP' )
      ENDIF

#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Preserve current practice in the std GEOS-Chem: (bmy, 10/24/12)
      ! 
      !    Copy the chemical species from CSPEC (actual troposphere) 
      !    to CSPEC_FULL (potential troposphere) array.  We only need 
      !    to do this if the variable troposphere is turned on. 
      !    (bdf, phs, bmy, 10/3/06)
      !
      !    NOTE: 
      !    (1) This has to be placed at the end of CHEMDR, after the
      !         call to the ND65 diagnostic DO_DIAG_PL. (phs, 6/3/08)
      !    (2) We also copy CSPEC to CSPEC_FULL if we want to write a 
      !        CSPEC_FULL restart file. (ccc, 2/26/09) 
      !-----------------------------------------------------------------
      IF ( LVARTROP .or. LSVCSPEC ) THEN
         CALL SAVE_FULL_TROP( State_Chm )

         !### Debug
         IF ( prtDebug ) THEN
            CALL DEBUG_MSG( '### CHEMDR: after SAVE_FULL_TROP' )
         ENDIF
      ENDIF
#endif

      ! Free pointer
      NULLIFY( STT )

      !=================================================================
      ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
      !=================================================================
      FIRSTCHEM = .FALSE.

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### Now exiting CHEMDR !' )
      ENDIF

      END SUBROUTINE CHEMDR
!EOC



