! $Id: ch3i_mod.f,v 1.2 2010/03/15 19:33:24 ccarouge Exp $
      MODULE CH3I_MOD
!
!******************************************************************************
!  Module CH3I_MOD contains emissions and chemistry routines for the CH3I
!  (Methyl Iodide) simulation. (bmy, 1/23/02, 9/27/06)
!
!  Module Routines:
!  ============================================================================
!  (1 ) OPEN_CH3I_FILES : Opens CH3I emissions files and reads data
!  (2 ) EMISSCH3I       : Emits CH3I from various sources into the STT array
!  (3 ) CHEMCH3I        : Performs CH3I chemistry on the STT tracer array
! 
!  GEOS-CHEM modules referenced by ch3i_mod.f
!  ============================================================================
!  (1 ) biofuel_mod.f   : Module w/ routines to read biofuel emissions
!  (2 ) biomass_mod.f   : Module w/ routines to read biomass emissions
!  (3 ) bpch2_mod.f     : Module w/ routines for binary punch file I/O
!  (4 ) dao_mod.f       : Module w/ arrays for DAO met fields  
!  (5 ) diag_mod.f      : Module w/ GEOS-CHEM diagnostic arrays
!  (6 ) diag_pl_mod.f   : Module w/ routines for prod & logs diag's
!  (7 ) error_mod.f     : Module w/ NaN and other error check routines
!  (8 ) file_mod.f      : Module w/ file unit numbers and error checks
!  (9 ) tracerid_mod.f  : Module w/ pointers to tracers & emissions
!  (10) transfer_mod.f  : Module w/ routines to cast & resize arrays
!  (11) uvalbedo_mod.f  : Module w/ routines to read UV albedo data
!
!  References
!  ============================================================================
!  (1 ) Bell, N. et al, "Methyl Iodide: Atmospheric budget and use as a tracer
!        of marine convection in global models", J. Geophys. Res, 107(D17), 
!        4340, 2002.
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!  (4 ) Wanninkhof, R., "Relation between wind speed and gas exchange over
!        the ocean", J. Geophys. Res, 97, 7373-7382, 1992.
!
!  NOTES:
!  (1 ) Removed obsolete code from 1/15/02 (bmy, 4/15/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (3 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (4 ) Now reference "file_mod.f" (bmy, 8/2/02)
!  (5 ) Updated call to INPHOT (bmy, 8/23/02)
!  (6 ) Now references BXHEIGHT from "dao_mod.f".  Now also references F90
!        modules "error_mod.f" and "tracerid_mod.f". (bmy, 11/6/02)
!  (7 ) Now references "grid_mod.f" and the new "time_mod.f" (bmy, 2/10/03)
!  (8 ) Added modifications for SMVGEAR II.  Removed reference to "file_mod.f".
!        (bdf, bmy, 4/21/03)
!  (9 ) Now references "directory_mod.f".  Now references "diag_pl_mod.f".
!        (bmy, 7/20/04)
!  (10) Now can read data for both GEOS and GCAP grids.  Now use Nightingale
!        et al formulation for piston velocity Kw. (bmy, 8/16/05)
!  (11) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (12) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE CMN_SIZE_MOD  ! Size parameters

      IMPLICIT NONE

      REAL*8, ALLOCATABLE :: OCDATA(:,:) ! ocean field for emiss. flux
      REAL*8, ALLOCATABLE :: EFCH4R(:,:) ! emission flux from rice
      REAL*8, ALLOCATABLE :: EFCH4W(:,:) ! emission flux from wetlands
      REAL*8              :: CH3ISUM(5)          ! sum of emissions in kg/yr

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE OPEN_CH3I_FILES( THISMONTH )
!
!******************************************************************************
!  Subroutine OPEN_CH3I_FILES loads surface emission fields for CH3I
!  (mgs, 3/15/99; bmy, hsu, 3/24/00,. bmy, 6/19/01, 10/3/05)
!
!  As of 16 June 1999, scale factors are applied in emissch3i.f (mgs)
!  and we use monthly RADSWG fields instead of NPP.
!
!  This routine is called at the first emission time step and on the
!  first of each month
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Shortwave radiation at the ground ...
!        *** used to be:
!        Ocean net primary productivity is used to estimate CH3I surface 
!        water concentration: parametrization derived from bilinear fit 
!        of ship cruise data with NPP from Rutgers University and RADSWG.
!        (surface water concentration should not exceed 8 ng/L)
!  (2 ) CH3I emissions from rice and wetlands use Fung's CH4 emission
!        inventory scaled with a constant factor from BIBLE observations.
!  (3 ) Added "CMN_SETUP" so that the proper path name to the /data/ctm
!        directories can be supplied. (bmy, 3/18/99)
!  (4 ) Trap I/O errors with subroutine IOERROR (bmy, 5/27/99)
!  (5 ) OCDATA now holds the aqueous CH3I in [ng/L], as read from disk.
!        No further unit conversion is necessary (hsu, bmy, 3/24/00)
!  (7 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (7 ) Now use function GET_TAU0 (from "bpch2_mod.f") to return the TAU0 
!        value used to index the binary punch file. (bmy, 7/20/00)
!  (8 ) Convert all input files to binary punch file format -- now call
!        READ_BPCH2 to read all binary punch files.  Also use GET_RES_EXT()
!        from BPCH2_MOD to get the proper extension string. (bmy, 8/8/00)
!  (9 ) Now read all CH3I files from the DATA_DIR/CH3I subdirectory.
!        Also updated comments & made cosmetic changes.  Also removed
!        reference to "CMN", which is not needed here. (bmy, 6/19/01)
!  (10) Now use routine TRANSFER_2D from "transfer_mod.f" to cast from REAL*4
!        to REAL*8 and to copy 2-D data to an array of size (IIPAR,JJPAR).
!        Also use 3 arguments (M/D/Y) in call to GET_TAU0.(bmy, 9/27/01)
!  (11) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (12) Now bundled into "ch3i_mod.f" Updated comments, cosmetic changes.
!        (bmy, 1/23/02)
!  (13) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (14) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (15) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) ::  THISMONTH    ! month of the year (1, 2, .., 12)

      ! Local variables
      REAL*4              :: Q1(IIPAR,JJPAR,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME 

      !=================================================================
      ! OPEN_CH3I_FILES begins here!
      !
      ! Get the TAU0 value for this month (use "generic" year 1985)
      !=================================================================
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      !=================================================================
      ! Read ocean field 
      !=================================================================

!      ! Uncomment this to read Ocean NPP
!      FILENAME = TRIM( DATA_DIR )       // 
!     &            'CH3I/ocean_npp.geos.' // GET_RES_EXT()

      ! Uncomment this to read aqueous CH3I
      FILENAME = TRIM( DATA_DIR )   // 
     &           'CH3I/ocean_ch3i.' // GET_NAME_EXT_2D() //
     &           '.'                // GET_RES_EXT()


      ! Read Caq in [ng/L] from the binary punch file
      CALL READ_BPCH2( TRIM( FILENAME ), 'IJ-AVG-$', 71, 
     &                 XTAU,              IIPAR,     JJPAR,     
     &                 1,                 Q1,        QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( Q1(:,:,1), OCDATA )

      !=================================================================
      ! Read rice paddy emissions
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) // 
     &           'CH3I/ch4_rice.' // GET_NAME_EXT_2D() //
     &           '.'              // GET_RES_EXT()

      ! Read CH4 rice paddy emissions in [kg/m2/s] 
      CALL READ_BPCH2( TRIM( FILENAME ), 'CH4-SRCE', 1, 
     &                 XTAU,              IIPAR,     JJPAR,     
     &                 1,                 Q1,        QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( Q1(:,:,1), EFCH4R )

      !=================================================================
      ! Read ocean DOC emissions and mixed-layer temp 
      !=================================================================
!      FILENAME = TRIM( DATA_DIR )       // 
!     &           'CH3I/ocean_DOC.geos.' // GET_RES_EXT()
! 
!      ! Read Short-Lived DOC emissions in [ng/L] 
!      CALL READ_BPCH2( FILENAME, 'DOC-SRCE', 1, XTAU,    
!     &                 IIPAR,     JJPAR,     1, Q1 )
! 
!      ! extract window
!      DOC_S(1:IIPAR,1:JJPAR) = Q1(1+I0:IIPAR+I0,1+J0:JJPAR+J0,1)
! 
!      ! Read Long-Lived DOC emissions in [ng/L] 
!      CALL READ_BPCH2( FILENAME, 'DOC-SRCE', 2, XTAU,    
!     &                 IIPAR,     JJPAR,     1, Q1 )
! 
!      ! extract window
!      DOC_L(1:IIPAR,1:JJPAR) = Q1(1+I0:IIPAR+I0,1+J0:JJPAR+J0,1)
! 
!      ! Read Mixed-Layer temp [T] 
!      CALL READ_BPCH2( FILENAME, 'DOC-SRCE', 3, XTAU,    
!     &                 IIPAR,     JJPAR,     1, Q1 )
! 
!      ! extract window
!      MLT(1:IIPAR,1:JJPAR) = Q1(1+I0:IIPAR+I0,1+J0:JJPAR+J0,1)

      !=================================================================
      ! Read CH4 wetland emissions
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) // 
     &           'CH3I/ch4_wetl.' // GET_NAME_EXT_2D() //
     &           '.'              // GET_RES_EXT()

      ! Read CH4 rice paddy emissions in [kg/m2/s] 
      CALL READ_BPCH2( TRIM( FILENAME ), 'CH4-SRCE', 2, 
     &                 XTAU,              IIPAR,     JJPAR,     
     &                 1,                 Q1,        QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( Q1(:,:,1), EFCH4W )

      ! Return to calling program
      END SUBROUTINE OPEN_CH3I_FILES

!------------------------------------------------------------------------------

      SUBROUTINE EMISSCH3I( am_I_Root, Input_Opt,
     &                      State_Met, State_Chm, RC )
!
!******************************************************************************
!  Subroutine EMISSCH3I (mgs, bmy, 11/23/98, 9/27/06) specifies methyl 
!  iodide (CH3I) emissions from the following sources:
!
!    Ocean: use correlation of surface water CH3I with net ocean
!           primary productivity and short wave radiation to get 
!           global fields of CH3I surface water concentrations.
!           Then compute sea-air exchange according to Liss&Slater, 1974
!           Use ND36 to get detailed output of ocean emissions!
!           (only shortwave radiation used, mgs, 06/21/99)
!   
!    Biomass burning: use CO emission database (J.A. Logan) and
!           scale with 0.4x10-6 v/v [Ferek et al., 1998]
!
!    Biofuel burning: same as biomass burning
!
!    Rice paddies and wetlands: use CH4 emission inventory from
!           Fung et al. [1991] and scale with 7.4x10-5 g/g
!           (BIBLE data over Japan, Blake pers. com., 1999)
!
!    Soil fumigation: CH3I could become the major replacement for
!           CH3Br in the future. Right now, we have no soil emissions,
!           but depending on availability, we might put the CH3Br
!           inventory from Brasseur et al., 1998 into the model
!           at some point.
!
!  Emissions from rice paddies, wetlands and biofuels are rather 
!  speculative at this point. There appears to be another terrestrial 
!  source from higher plants.  However, this is practically un-quantifiable.
!
!  NOTES:
!  (1 ) Starting point: cleaned version of EMISSRN (1.5 from bmy)
!  (2 ) initial version: simply specify surface layer concentrations 
!        for all ocean grid boxes (mgs, 11/20/98)
!  (3 ) As of 11/20/98, the following sources of CH3I are now used:
!        (a) Tracer #1 : CH3I from oceans
!        (b) Tracer #2 : CH3I from biomass burning (scaled from CO values)
!        (c) Tracer #3 : CH3I from wood burning (scaled from CO values)
!  (4 ) Added FIRSTEMISS as an argument...useful for later reference
!        (bmy, 11/23/98)
!  (5 ) Added ND29 diagnostics for CO woodburning and biomass burning 
!        (bmy, 11/23/98)
!  (6 ) Add FRCLND as an argument (bmy, 1/11/99)
!  (7 ) Replace constant surface concentrations for ocean source with
!        flux parametrization and add rice and wetland tracers: (mgs, 03/12/99)
!        (d) Tracer #4 : CH3I from rice paddies
!        (e) Tracer #5 : CH3I from wetlands
!  (8 ) DIAG36 is used for emission fluxes in ng/m2/s and surface water
!        concentrations in ng/L and a log of the ocean atmosphere exchange 
!        coefficient in cm/h. DIAG29 traces biomass burning and woodburning
!        CO emissions in ???.
!  (9 ) Added LOGMONTH for logging CH3I monthly mean output (mgs, 3/24/99)
!  (10) Now use F90 syntax for declarations.  Also added the OUTLOG
!       flag for sending monthly sums to a log file (bmy, 3/24/99)
!  (11) Fixed bugs in the expressions for H and FLUX (mgs, bmy, 5/15/99)
!  (12) Now uses bilinear correlation with NPP and RADSWG for ocean source
!       (before only NPP)   (mgs, 16 Jun 1999
!  (13) FRCLND removed as argument, because CMN_DEP now included 
!       (mgs, 06/16/99)
!  (14) Ocean emissions now differetn parametrizations for 3 latitude regions.
!        Emissions protocolled in more detail.
!  (15) added LASTEMISS flag for final summary output (mgs, 06/28/99)
!  (16) Replaced AIJ with AD36 allocatable array (bmy, 3/28/00)
!  (17) Removed obsolete code (bmy, 4/14/00)
!  (18) Now reference AIRVOL and TS from "dao_mod.f" instead of from
!        common block header files (bmy, 6/23/00)
!  (19) Eliminate obsolete code from 6/26/00 (bmy, 8/31/00)
!  (20) Added references to F90 modules "biomass_mod.f" and "biofuel_mod.f".
!       Also, TWOODIJ is now called BIOFUEL.  Finally, BURNEMIS is now 
!        referenced with IREF = I + I0 and JREF = J + J0. (bmy, 9/11/00)
!  (21) Removed obsolete code from 9/12/00 (bmy, 12/21/00)
!  (22) Now use IDBFCO to reference the biofuel CO emissions.  Also make
!        sure that IDBCO and IDBFCO are not zero. (bmy, 3/20/01) 
!  (23) Eliminated obsolete commented-out code (bmy, 4/20/01)
!  (24) Now prompt user to check IDBCO and IDBFCO in "tracer.dat" if
!        these switches are turned off.  Also now read all data files
!        from the CH3I subdirctory of DATA_DIR. (bmy, 6/19/01)
!  (25) BIOFUEL (N,IREF,JREF) is now BIOFUEL(N,I,J).  BURNEMIS(N,IREF,JREF)
!        is now BURNEMIS(N,I,J). (bmy, 9/28/01)
!  (26) Removed obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (27) Now bundled into "ch3i_mod.f".  Updated comments, cosmetic 
!        changes.  Removed LASTEMISS as an argument. (bmy, 1/23/02)
!  (28) Now reference file units from "file_mod.f" (bmy, 8/2/02)
!  (29) Now reference BXHEIGHT from "dao_mod.f".  Also references IDBCO and
!        IDBFCO from "tracerid_mod.f".  Now make FIRSTEMISS a local SAVEd
!        variable. (bmy, 11/15/02)
!  (30) Now use GET_AREA_M2 from "grid_mod.f" to compute grid box surface
!        areas.  Removed references to DXYP.  Now use functions GET_DAY,
!        GET_GMT, GET_TS_EMIS from the new "time_mod.f". (bmy, 2/10/03)
!  (31) Now reference STT & N_TRACERS from "tracer_mod.f".  Now reference
!        LEMIS from "logical_mod.f". (bmy, 7/20/04)
!  (32) Now modified for new "biomass_mod.f" (bmy, 4/5/06)
!  (33) BIOMASS(:,:,IDBCO) from "biomass_mod.f" is now in units of 
!        [molec CO/cm2/s].  Adjust unit conversion accordingly. (bmy, 9/27/06)
!  (34) IDBCO is now defined in TRACERID_MOD. (fp, hotp 7/31/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!******************************************************************************
!
      ! Reference to F90 modules
      USE BIOFUEL_MOD,        ONLY : BIOFUEL,   BIOFUEL_BURN
      USE BIOMASS_MOD,        ONLY : BIOMASS
      USE DIAG_MOD,           ONLY : AD29,      AD36
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE TIME_MOD,           ONLY : GET_DAY,   GET_GMT,
     &                               GET_MONTH, GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : IDBCO
      USE TRACERID_MOD,       ONLY : IDBFCO

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic switches
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

      ! Local common blocks (shared with OPEN_CH3I_FILES)
      REAL*8 :: OCDATA(IIPAR,JJPAR)   ! scaled ocean net primary prod.
      REAL*8 :: EFCH4R(IIPAR,JJPAR)  ! scaled emission flux from rice
      REAL*8 :: EFCH4W(IIPAR,JJPAR)  ! scaled emission flux from wetlands
      REAL*8 :: CH3ISUM(5)           ! sum of emissions in kg/yr

      ! Local variables:
      LOGICAL, SAVE      :: FIRSTEMISS = .TRUE.
      INTEGER            :: I, IJLOOP, J,    L,    N

      REAL*8             :: DTSRCE,  BIO_CO,  WOOD_CO
      REAL*8             :: SCALE,   MOLRAT,  BXHEIGHT_CM, AREA_M2

      ! local surface air temp. in K or degC
      REAL*8             :: TK, TC    

      ! Henry, Schmidt number and exchange coeff.
      REAL*8             :: H, Sc, KW 

      ! emission flux in kg m^-2 s^-1
      REAL*8             :: FLUX      
     
      ! surface layer gas concentration (ocean)
      REAL*8             :: CGAS      

      ! For Nightingale sea-air transfer formulation
      REAL*8             ::  W10
      REAL*8, PARAMETER  ::  ScCO2 = 600.0d0

      ! molar volume ratio CH3I/CH3Br for computation
      ! of exchange coefficient ocean atmosphere = (62.9/52.9)**0.6
      REAL*8, PARAMETER  :: MVR__ = 1.1094736d0    

      ! molar gas constant
      REAL*8, PARAMETER  :: R = 8.314d0

      ! (molecules/mole)^-1
      REAL*8, PARAMETER  :: XMOL = 1.d0 /6.0225d+23 

      ! molar weights (g/mole)
      REAL*8, PARAMETER  :: FMOL_CO   = 28.0d0
      REAL*8, PARAMETER  :: FMOL_CH3I = 141.939d0

      ! month number used for logging of emissions
      INTEGER            :: LOGMONTH    

      ! ECH3I: Emission ratio mole CH3I / mole CO for biomass burning
      REAL*8, PARAMETER  :: ECH3I = 4.0d-6 

      ! constant CH3I surface water conc. for tropical latitudes (0-20 deg
      ! ** use mean of 0.721D0 for a maximum estimate, median is 0.480D0
      REAL*8, PARAMETER  :: OCTROP   = 0.480D0     ! median

      ! scale factor for rice paddy emissions ( g CH3I / g CH4 )
      REAL*8, PARAMETER  :: SCCH4R = 7.4d-5

      ! scale factor for wetland emissions
      REAL*8, PARAMETER  :: SCCH4W = SCCH4R

      ! External functions
      REAL*8, EXTERNAL   :: BOXVL       ! grid box volume in cm^3
      REAL*8, EXTERNAL   :: SFCWINDSQR  ! square of surface wind speed (10m)

      ! output flag
      LOGICAL, PARAMETER :: OUTLOG = .FALSE.

      ! For fields from Input_Opt
      LOGICAL         :: LEMIS
      INTEGER         :: N_TRACERS

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! EMISSCH3I begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LEMIS     = Input_Opt%LEMIS
      N_TRACERS = Input_Opt%N_TRACERS

      !=================================================================
      ! If this is the first emission step, do the following...
      !=================================================================
      IF ( FIRSTEMISS ) THEN 

         ! Make sure that NTRACE = 5 since we have 5 CH3I tracers.
         !IF ( NTRACE /= 5 ) NTRACE = 5
         IF ( N_TRACERS /= 5 ) N_TRACERS = 5
         

         ! Make sure that emissions are turned on. 
         ! CH3I simulation doesn't make sense w/o emissions
         IF ( .not. LEMIS ) THEN
            PRINT*,'**** LEMIS=.FALSE.! I turn emissions on now!'
            LEMIS = .TRUE.
         ENDIF

         ! Output of accumulated emissions
         WRITE(97,*) ' CH3I emission log'
         WRITE(97,*) ' Total emissions in kg'
         WRITE(97,*)

         CH3ISUM = 0.0D0
      ENDIF
 
!### Debug -- comment out for now (bmy, 2/10/03)
!      !### Debug output in unit 97 ... set OUTLOG flag if wanted
!      !### print the first 24 hours, then every 24 hours
!      IF ( OUTLOG ) THEN
!         IF ( TAU-TAUI < 24 .OR. MOD(FLOOR(TAU-TAUI),24) == 0 ) THEN
!            write(97,*) 'before emiss, TAU=',tau
!            DO N = 1, NTRACE
!               write(97,'(A,I4,1p,E12.4)') ' N, SUM STT(N) = ', 
!     &                                   N, SUM( STT(:,:,:,N) )
!            ENDDO
!         ENDIF
!      ENDIF

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! DTSRCE = the number of seconds between emissions
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! On the first of each month and at the first emission time step:
      !
      ! (1) Read in the following fields
      !     Aqueous CH3I concentrations (ocean_ch3i.geosX.4x5)
      !     methane emissions from rice     (ch4_rice.4x5.MM)
      !     methane emissions from wetlands (ch4_wetl.4x5.MM)
      !
      ! (2) Log global total of emissions and reset sums
      !=================================================================
      IF ( FIRSTEMISS .OR. 
     &     ( GET_DAY() == 1 .and. GET_GMT() < 1d-3 ) ) THEN 

         CALL OPEN_CH3I_FILES( GET_MONTH() )

         ! apply scaling factors
         EFCH4R = EFCH4R * SCCH4R
         EFCH4W = EFCH4W * SCCH4W

         ! for surface water concentrations, use results from correlation
         ! analysis with NPP and solar radiation (i.e. only RADSWG):
         ! use latitude centers (YLMID) for different regions
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! set all values over land to zero
            IF ( State_Met%FRCLND(I,J) >= 0.8 ) OCDATA(I,J) = 0.0D0
         
            ! make sure we have no negative concentrations
            IF ( OCDATA(I,J) < 0. ) OCDATA(I,J) = 0.0D0
         ENDDO
         ENDDO

      ENDIF

!### Debug -- comment out for now (bmy, 2/10/03)
!      !### Debug output to fort.97
!      IF ( JDATE == 1 .AND. TOFDAY < 1.D-3 ) THEN
!         LOGMONTH = MONTH-1
!         IF ( LOGMONTH <= 0 ) LOGMONTH = 12
!
!         WRITE(97,*)
!         WRITE(97,*) 
!     &        ' ---------------------------------------------'
!         WRITE(97,*) ' Accumulated global emissions for month ',
!     &        LOGMONTH,': '
!         WRITE(97,'(3X,A20,1p,E12.4)') 
!     &        'ocean : ',CH3ISUM(1)
!         WRITE(97,'(3X,A20,1p,E12.4)') 
!     &        'biomass burning : ',CH3ISUM(2)
!         WRITE(97,'(3X,A20,1p,E12.4)') 
!     &        'wood burning : ',CH3ISUM(3)
!         WRITE(97,'(3X,A20,1p,E12.4)') 
!     &        'rice paddies: ',CH3ISUM(4)
!         WRITE(97,'(3X,A20,1p,E12.4)') 
!     &        'wetlands : ',CH3ISUM(5)
!         WRITE(97,*) 
!     &        ' ---------------------------------------------'
!         WRITE(97,*)
!
!         ! flush the output buffer
!         CALL FLUSH( 97 )
!
!         CH3ISUM = 0.0D0
!      ENDIF

      !=================================================================
      !          --------------------------------------------
      !                  Tracer #1: CH3I from Oceans
      !          --------------------------------------------
      !
      ! OCDATA contains estimated surface water CH3I concentrations 
      ! derived from ocean net primary productivity (in ng/L).  
      ! The net emission flux is given by
      !      F = KW ( Caq - Cg*H )
      !
      ! KW is exchange parameter (piston velocity) and given by
      !      KW = 0.31 u^2 ( Sc/660 )^(-1/2)     (cm/h) 
      !      [Wanninkhof et al., 1992]
      ! NOTE: As of 8/16/05, we now use the Nightingale et al [2000b]
      ! formulation for piston velocity which is:
      !      Kw   = ( 0.24 * u^2 + 0.061d0*u ) * SQRT( 600/Sc )  
      !
      ! u^2 is square of surface wind speed (10m above ground) in m/s
      !
      ! Sc is Schmidt number:
      !      Sc = (62.9/52.9)^0.6* (2004.-93.5*T+1.39*T^2)
      !
      ! with T in degC  [Moore and Groszko, 1999] 
      ! 660.0 is Schmidt number for CO2 in seawater (normalization)
      !
      ! Caq is the surface water concentration, 
      ! Cg  is the gas-phase concentration, and 
      ! H   is the (dimensionless!) Henry coefficient:
      !      H^-1 = 0.14 exp(-4300 * (T-298)/(T*298)) * R * T / 101.325  
      !      [R. Sander]
      !
      ! here T is in K (!)
      !
      ! To convert Cg from kg/gridbox into ng/L: * 1.d12*1.d-3/AIRVOL
      !
      ! To convert cm/h*ng/L to kg/m^2/s : *1.d-11/3600.
      !
      ! Since CH3I exhibits a pretty strong gradient near the surface,
      ! we may have to adjust the "surface" gas-phase concentration in 
      ! the future??
      !
      ! Apply emission flux only for grid boxes that contain at least
      ! 20% non-land (FRCLND) and where the surface temperature is above
      ! -2 degC (a little arbitrary).
      !
      ! NOTES: 
      ! (1 ) grid box surface area in m^2 is given by DXYP(JREF). 
      !       Attention: this is not window size!
      ! (2 ) Fixed bug with Henry's definition. Old code erroneously 
      !       used two different definitions of H and did not correctly 
      !       convert from H in mol/atm to dimensionless H. 
      !       (mgs, 05/14/1999)
      ! (3 ) Now we read in the aqueous CH3I (Caq) from disk into the 
      !       OCDATA array.  OCDATA now has units of [ng/L].  
      !       (hsu, bmy, 3/24/00)
      ! (4 ) DXYP(JREF) is now replaced by GET_AREA_M2(J) (bmy, 2/4/03)
      !=================================================================
      N = 1
      L = 1
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( I, J, 1 )

         IF ( State_Met%TS(I,J)     >= 271.15 .and.
     &        State_Met%FRCLND(I,J) <  0.8    ) THEN

            ! surface air temp [K]
            TK = State_Met%TS(I,J)        

            ! sea surface temp [C] 
            ! (use TS as surrogate for SST) 
            TC = TK - 273.15      

            ! Henry's law constant [unitless]
            H = 0.14*exp(-4300.*(TK-298.)/(TK*298.) ) * R / 101.325 * TK

            ! Schmidt # [unitless]
            Sc = MVR__ * ( 2004. - 93.5*TC + 1.39*TC**2 )
            
            ! 10-m wind speed 
            W10  = SQRT( SFCWINDSQR( I, J,
     &                               State_Met%U10M,
     &                               State_Met%V10M ) )

            ! Piston velocity [cm/h], cf Nightingale et al [2000b] 
            ! (swu, bmy, 8/16/05)
            Kw   = ( 0.24d0*W10*W10 + 0.061d0*W10 ) * SQRT( ScCO2/Sc )  

            ! convert gas-phase tracer mass to concentration in ng/L
            CGAS = STT(I,J,L,N) * 1.0D9 / State_Met%AIRVOL(I,J,L)

            ! Emission of CH3I from the ocean 
            FLUX = KW * ( OCDATA(I,J) - CGAS*H )

!###            !### debug output
!###            WRITE(97,9777) 'I,J,H,Sc,KW,CGAS,OCDATA(I,J),FLUX:',
!###  &         I,J,H,Sc,KW,CGAS,OCDATA(I,J),FLUX
!### 9777       FORMAT(A,2I3,1p,6E12.3)

            !###! make sure, flux is positive (really ??)
            !###IF (FLUX.LT.0.0D0) FLUX = 0.0D0

            ! convert flux to kg/m^2/time step
            FLUX = FLUX * ( 1.0D-11 / 3600.D0 ) * DTSRCE

            ! Add to diagnostic array
            IF (ND36.GT.0) THEN
               AD36(I,J,N) = AD36(I,J,N) + FLUX * 1.D+12
               AD36(I,J,6) = AD36(I,J,6) + KW * OCDATA(I,J) * DTSRCE
               AD36(I,J,7) = AD36(I,J,7) + KW * CGAS * H * DTSRCE
               !### debug output: store terms seperately !!
               !###AD36(I,J,N) = AD36(I,J,N) + FLUX * 1.D+12
               !###AD36(I,J,6) = AD36(I,J,6) + OCDATA(I,J)*DTSRCE
               !###AD36(I,J,7) = AD36(I,J,7) + KW*DTSRCE
            ENDIF

!###            !### debug output
!###            IF (I.eq.49 .AND. J.eq.26) THEN
!###               WRITE(97,'(A,2I5,1p,4e12.4)')
!###     &             'I,J, Caq, Cgas, KW, Ocean flux = ',
!###     &              I,J,OCDATA(I,J),CGAS,KW,FLUX
!###
!###               WRITE(97,'(A,1p,3e12.4)') 
!###     &              'STT,AIRVOL(I,J,L), H = ',
!###     &              STT(I,J,L,N),AIRVOL(I,J,L),H
!###            ENDIF

            ! convert flux to kg/gridbox/time step
            FLUX = FLUX * AREA_M2

            !### Debug
            !###write(97,'(A,1p,E12.4)') 'rescaled flux = ',FLUX

            ! add to tracer mass and to global sum
            STT(I,J,L,N) = STT(I,J,L,N) + FLUX

            !### make sure we get no negative concentrations
            IF (STT(I,J,L,N).LT.0.) THEN
               STT(I,J,L,N) = 0.0D0
               FLUX = 0.0D0
            ENDIF

            CH3ISUM(N) = CH3ISUM(N) + FLUX

            IF (ND36.GT.0) THEN
               AD36(I,J,8) = AD36(I,J,8) + FLUX
            ENDIF
         ENDIF
      ENDDO
      ENDDO

      !=================================================================
      !          --------------------------------------------
      !             Tracer #2: CH3I from Biomass burning
      !          --------------------------------------------
      !
      !  Biomass burning CO is stored in BURNEMIS(IDBCO,:,:) 
      !  in [molec/cm3/s].  Convert to kg CH3I as follows:
      !
      !     FLUXKG = flux * molar emission factor 
      !                   * mole weight CH3I / molec/mole
      !                   * grid box volume      
      !=================================================================
 
      ! Convert biomass CO into biomass CH3I
      N = 2
      L = 1
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box height [cm]
         BXHEIGHT_CM = State_Met%BXHEIGHT(I,J,L) * 100d0

         !-----------------------------------------------------------------
         ! Get emission flux in kg/cm3/time step
         !FLUX = ECH3I * BIOMASS(I,J,IDBCO)
         !-----------------------------------------------------------------

         ! Convert [molec/cm2/s] to [kg/cm3/timestep]
         FLUX = ECH3I * BIOMASS(I,J,IDBCO) / BXHEIGHT_CM
         FLUX = FLUX * 1.0D-3 * FMOL_CH3I * XMOL * DTSRCE

         ! Add to diagnostic array as kg/m2/time step
         IF ( ND36 > 0 ) THEN
            AD36(I,J,N) = AD36(I,J,N) + FLUX*BXHEIGHT_CM*1.0D4*1.0D+12
         ENDIF

         ! Convert to kg/grid box/time step
         FLUX = FLUX * BOXVL(I,J,L,State_Met)

         ! add to tracer mass and to global sum
         STT(I,J,L,N) = STT(I,J,L,N) + FLUX
         CH3ISUM(N) = CH3ISUM(N) + FLUX
      ENDDO
      ENDDO
      
      !=================================================================
      !          --------------------------------------------
      !                   Tracer #3: Wood burning
      !          --------------------------------------------
      !
      !  Wood burning CO is stored in BIOFUEL(IDBFCO,:,:) 
      !  in [molec/cm^3/s].  Proceed as in biomass burning emissions
      !=================================================================

      ! Make sure CO is a biofuel tracer
      IF ( IDBFCO == 0 ) THEN
         CALL ERROR_STOP( 'IDBFCO=0, check "tracer.dat"', 'EMISSCH3I' )
      ENDIF

      ! Now reference routine BIOFUEL_BURN from "biofuel_mod.f" (bmy, 9/12/00)
      CALL BIOFUEL_BURN
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      N = 3
      L = 1
      DO J = 1, JJPAR
      DO I = 1, IIPAR
 
         ! Get flux as kg/cm3/time step
         ! Now use IDBFCO to index biofuel burning CO (bmy, 3/20/01)
         FLUX = ECH3I * BIOFUEL(IDBFCO,I,J)
         FLUX = FLUX * 1.0D-3 * FMOL_CH3I * XMOL * DTSRCE

         ! Add to diagnostic array as kg/m2/time step
         IF ( ND36 > 0 ) THEN
            BXHEIGHT_CM = State_Met%BXHEIGHT(I,J,L) * 100d0
            AD36(I,J,N) = AD36(I,J,N)  + FLUX * BXHEIGHT_CM *
     &                    1.0D4 * 1.0D + 12
         ENDIF

         FLUX = FLUX * BOXVL(I,J,L,State_Met)
         
         ! add to tracer mass and to global sum
         STT(I,J,L,N) = STT(I,J,L,N) + FLUX
         CH3ISUM(N) = CH3ISUM(N) + FLUX
      ENDDO
      ENDDO
   
      !=================================================================
      !          --------------------------------------------
      !                Tracer #4: Rice paddy emissions
      !          --------------------------------------------
      !
      !  EFCH4R contains emission flux in kg(CH3I)/m^2/s.
      !  Simply convert to kg/grid box/time step and add to tracer mass
      !
      !  NOTES:
      !  (1) everything should be in window size except DXYP(J)
      !  (2) DXYP(JREF) is now replaced by GET_AREA_M2(J). (bmy, 2/4/03)
      !=================================================================
      N = 4
      L = 1
      DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid box surface area [m2]
            AREA_M2 = GET_AREA_M2( I, J, 1 )

            ! First compute flux as kg/m2/time step
            FLUX = EFCH4R(I,J)
            FLUX = FLUX * DTSRCE

            ! Add to diagnostic array as kg/m2/time step
            IF ( ND36 > 0 ) THEN
               AD36(I,J,N) = AD36(I,J,N) + FLUX*1.0D+12
            ENDIF

            ! Now convert to kg/grid box/time step
            FLUX = FLUX * AREA_M2

            ! add to tracer mass and to global sum
            STT(I,J,L,N) = STT(I,J,L,N) + FLUX
            CH3ISUM(N) = CH3ISUM(N) + FLUX
         ENDDO
      ENDDO
         
      !=================================================================
      !          --------------------------------------------
      !                 Tracer #5: Wetland emissions
      !          --------------------------------------------
      !
      !  EFCH4W contains emission flux in kg(CH3I)/m^2/s.
      !  Simply convert to kg/grid box/time step and add to tracer mass
      !
      !  NOTES:
      !  (1) everything should be in window size except DXYP(J)
      !  (2) DXYP(J) is now replaced by GET_AREA_M2(J) (bmy, 2/4/03)
      !=================================================================
      N = 5
      L = 1
      DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid box surface area [m2]
            AREA_M2 = GET_AREA_M2( I, J, 1 )

            ! First compute flux as kg/m2/time step
            FLUX = EFCH4W(I,J)
            FLUX = FLUX * DTSRCE

            ! Add to diagnostic array as kg/m2/time step
            IF ( ND36 > 0 ) THEN
               AD36(I,J,N) = AD36(I,J,N) + FLUX*1.0D+12
            ENDIF

            ! Now convert to kg/grid box/time step
            FLUX = FLUX * AREA_M2

            ! add to tracer mass and to global sum
            STT(I,J,L,N) = STT(I,J,L,N) + FLUX
            CH3ISUM(N) = CH3ISUM(N) + FLUX
         ENDDO
      ENDDO

      !=================================================================
      ! ** Future : 
      !  (1) add soil fumigation (CH3I as replacement for CH3Br)
      !      [may require change in ND36!]
      !=================================================================

!### Debug -- comment out for now
!      !### Debug output in unit 97 ... set OUTLOG flag if wanted
!      !### print the first 24 hours
!      IF ( OUTLOG ) THEN
!          IF ( TAU-TAUI < 24 .OR. MOD(FLOOR(TAU-TAUI),24) == 0 ) THEN
!             WRITE( 97, * ) 'after emiss'
! 
!             DO N = 1, NTRACE
!                WRITE( 97, '(''N, SUM STT(N) = '',i4,es12.4)' )
!     &               N, SUM( STT(:,:,:,N) )
!             ENDDO
!          ENDIF
!
!          WRITE( 97, * )
!       ENDIF

      ! Free pointer
      NULLIFY( STT )

      ! Make sure the next time is not the first emission time step ;-)
      FIRSTEMISS = .FALSE.

      ! Return to calling program
      END SUBROUTINE EMISSCH3I 

!------------------------------------------------------------------------------

      SUBROUTINE CHEMCH3I( am_I_Root, Input_Opt, 
     &                     State_Met, State_Chm, STATUS )
!
!******************************************************************************
!  Subroutine CHEMCH3I performs loss chemistry for methyl iodide (CH3I).  
!  (mgs, bey, bmy, 11/20/98, 8/16/05)
!
!  If the LFASTJ C-preprocessor switch is set, then CHEMCH3I will invokes 
!  the FAST-J subroutines to compute local photolysis rates, which in 
!  turn determine local CH3I loss rates. Otherwise, a constant loss rate 
!  of 1/4 day is applied.
!
!  NOTES:
!  (1 ) Based on subroutine CHEMRN.F (bey, bmy, 1998)
!  (2 ) Edited comments and changed constant lifetime from 3 to 4 days. 
!        (mgs, 3/12/99)
!  (3 ) Now call INPHOT directly, rather than via FJ_INIT. (bmy, 10/4/99)
!  (4 ) Make sure fast-J files "ratj.d" and "jv_spec.dat" include
!        the information for CH3I branching ratios & quantum yields.
!  (5 ) CHEMCH3I calls READER.F and CHEMSET.F to read in the "m.dat" and
!        "chem.dat" files for CH3I.  These is necessary to ensure that the 
!        J-Value mapping from Harvard indices to UCI indices will be done 
!        correctly.
!  (6 ) CH3I loss will be computed from the surface to layer NSKIPL-1,
!        which is specified in "input.ctm".
!  (7 ) CH3I loss is now computed only for places where it is daylight
!        (i.e. where SUNCOS > 0).  This will prevent computing the
!        exponential where the J-Values would be zero.  (bmy, 11/23/98)
!  (8 ) Add J-Value diagnostic for ND22  (bmy, 11/23/98)
!  (9 ) Now use F90 syntax for declarations (bmy, 3/24/99)
!  (10) Now "comsol.h" only contains variables relevant to SLOW-J, so
!        we don't have to #include it here.  
!  (11) AD22 is now declared allocatable in "diag_mod.f". (bmy, 11/29/99)
!  (12) LTJV is now declared allocatable in "diag_mod.f". 
!        Also made cosmetic changes, and updated comments. (bmy, 3/17/00)
!  (13) Added ND65 diagnostic for CH3I loss (nad, bmy, 3/27/01)
!  (14) Now reference the UVALBEDO array directly from "uvalbedo_mod.f".
!        Remove ALBD from the argument list.  Updated comments, cosmetic
!        changes. (bmy, 1/15/02)
!  (15) Now bundled into "ch3i_mod.f". (bmy, 1/23/02)
!  (16) Removed obsolete code from 1/15/02 (bmy, 4/15/02)
!  (17) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (18) Now reference file unit IU_CTMCHEM from "file_mod.f" (bmy, 8/2/02)
!  (19) Now reference SUNCOS, OPTD from "dao_mod.f".  Now make FIRSTCHEM
!        a local SAVEd variable. (bmy, 11/15/02)
!  (20) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 2/11/03)
!  (21) Replace call to CHEMSET with call to READCHEM for SMVGEAR II.  Replace
!        NAMESPEC array with NAMEGAS array.   Removed reference to "file_mod.f"
!        since now the "smv2.log" file is opened in READER. (bdf, bmy, 4/21/03)
!  (22) Now reference STT and N_TRACERS from "ch3i_mod.f".  Also replace
!        NSKIPL-1 with LLTROP for now.  Now references AD65 from 
!        "diag_pl_mod.f". (bmy, 7/20/04)
!  (23) FAST-J is now the default, so we don't need the LFASTJ C-preprocessor
!        switch any more (bmy, 6/23/05)
!  (24) Now use Nightingale et al [2000b] formulation for piston velocity
!        (swu, bmy, 8/16/05)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!******************************************************************************
!
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD22, LTJV
      USE DIAG_PL_MOD,        ONLY : AD65
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! SPECNAME
      USE CMN_DIAG_MOD             ! ND36
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: STATUS      ! Success or failure?
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      REAL*8                 :: DTCHEM, RLRAD, RDLOSS, T1L
      REAL*8                 :: TCHEMA, JVALUE

      ! Hardwired flag for CH3I species name
      CHARACTER (LEN=4)      :: SPECNAME

      ! Local loop and index variables
      INTEGER                :: I, IFNC, IBRCH, J, L, LN22, N, NK, NR
      INTEGER                :: RC

      ! Hardwired Logical flag for FAST-J
      LOGICAL, PARAMETER     :: DO_FASTJ = .TRUE. 

      ! External functions
      REAL*8, EXTERNAL       :: FJFUNC 

      ! For fields from Input_Opt
      INTEGER               :: N_TRACERS

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHEMCH3I begins here!
      !=================================================================

      ! Assume success
      STATUS = GIGC_SUCCESS
      
      ! Get N_TRACERS field from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! convert DTCHEM from mn to sec.
      DTCHEM = GET_TS_CHEM() * 60d0

!-----------------------------------------------------------------------------
!### Debug output in unit 97 ...comment out if necessary (bmy, 11/23/98)
!###  print the first 24 hours
!       IF (TAU-TAUI.LT.24) THEN
!          write(97,*) 'before chem'
!          DO N = 1, NTRACE
!             write(97,'(A,I4,1p,E12.4)') ' N, SUM STT(N) = ',
!     &                                   N, SUM( STT(:,:,:,N) )
!          ENDDO
!      ENDIF
!-----------------------------------------------------------------------------

         !==============================================================
         ! If LFASTJ is defined in "define.h", then invoke FAST-J 
         ! subroutines to compute photolysis rates.  The loss rate
         ! of CH3I is dependent on the local photolysis rates.
         !
         ! Initialize FAST-J quantities on the first timestep
         !==============================================================
         IF ( FIRSTCHEM ) THEN

            ! Call READER and READCHEM to read "mglob.dat" and 
            ! "globchem.dat" (these are needed for the J-value mapping). 
            CALL READER( FIRSTCHEM )
            CALL READCHEM
            
            ! Call INPHOT to initialize the fast-J variables. 
            CALL INPHOT( LLTROP, NPHOT )
            
            ! Echo output
            WRITE( 6,'(a)' ) 'Using U.C.I Fast-J Photolysis'

            ! Reset FIRSTCHEM
            FIRSTCHEM = .FALSE.
         ENDIF

         !==============================================================
         ! For each chemistry time step, compute J-values and store 
         ! them in an internal array.  SUNCOS, OPTD, and UVALBEDO are 
         ! needed for FAST-J.
         !==============================================================
         CALL FAST_J( am_I_Root, State_Met, RC )
        
         !==============================================================
         ! NR is the loop over the number of reactions (NR=1 for now!)
         ! Compute the proper branch number for each reaction, using 
         ! the same algorithm from CALCRATE.F.
         ! 
         ! For each photo reaction, loop over the grid boxes (I-J-L) 
         ! and test whether the grid box is in sunlight or not.  If the 
         ! grid box is a daytime box, then extract the proper photo 
         ! rate for that box.  
         ! 
         ! The photo rate for each grid box is in (s^-1) so multiply 
         ! this by the number of seconds in the chemistry interval and 
         ! use that as the loss rate (i.e. the arg of the exponential).
         ! 
         ! You must specify NTRACE in "input.ctm".  The number of CH3I 
         ! tracers from different sources ranges from N=1 to N=NTRACE:
         !     N = 1: CH3I from oceans
         !     N = 2: CH3I from biomass burning
         !     N = 3: CH3I from wood burning
         !     N = 4: CH3I from rice paddies
         !     N = 5: CH3I from wetlands
         ! 
         ! Also redefine RDLOSS so that it is just the exponential term, 
         ! which can then be multiplied by the tracer STT in one step 
         ! (bmy, 1/11/99)
         !==============================================================
         DO NR = 1, NPHOT
            NK       = NRATES(NCS) + NR
            IFNC     = DEFPRAT(NK,NCS) + 0.01D0
            IBRCH    = 10.D0*(DEFPRAT(NK,NCS)-IFNC) + 0.5D0   
            SPECNAME = NAMEGAS(IRM(1,NK,NCS))            

            ! Maybe later can replace this w/ the ann mean tropopause...
            DO L = 1, LLTROP
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! SUNCOS > 0 means daytime!
               IF ( State_Met%SUNCOS(I,J) > 0d0 ) THEN
                  JVALUE = FJFUNC( I, J, L, NR, IBRCH, SPECNAME )
                  RLRAD  = JVALUE * DTCHEM
                  RDLOSS = EXP( -RLRAD )

                  ! Loop over all individual CH3I tracers
                  ! (which have the same loss rate)
                  DO N = 1, N_TRACERS
                     STT(I,J,L,N) = STT(I,J,L,N) * RDLOSS
                  ENDDO

                  ! ND22: J-value diagnostic
                  IF ( ND22 > 0 ) THEN
                     IF ( LTJV(I,J) > 0 .and. L <= LD22 ) THEN
                        AD22(I,J,L,1) = AD22(I,J,L,1) + JVALUE
                     ENDIF
                  ENDIF
                  
                  ! ND65: Loss rates for each tracer 
                  IF ( ND65 > 0 ) THEN
                     IF ( L <= LD65 ) THEN
                        DO N = 1, N_TRACERS
                           AD65(I,J,L,N) = AD65(I,J,L,N) + 
     &                          ( STT(I,J,L,N) * JVALUE * DTCHEM )
                        ENDDO
                     ENDIF
                  ENDIF      
               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDDO

!-----------------------------------------------------------------------------
!### Debug output in unit 97 ...comment out if necessary (bmy, 11/23/98)
!### print the first 24 hours
!       IF (TAU-TAUI.LT.24) THEN
!          write(97,*) 'after chem'
!          DO N = 1, NTRACE
!             write(97,'(A,I4,1p,E12.4)') ' N, SUM STT(N) = ',
!     &                                   N, SUM( STT(:,:,:,N) )
!          ENDDO
!       ENDIF
!-----------------------------------------------------------------------------

      ! Free pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE CHEMCH3I

#if defined( DEVEL )
!------------------------------------------------------------------------------
      SUBROUTINE SET_CH3I_MOD

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(
     &  OCDATA(IIPAR,JJPAR),      ! ocean field for emiss. flux
     &  EFCH4R(IIPAR,JJPAR),      ! emission flux from rice
     &  EFCH4W(IIPAR,JJPAR),      ! emission flux from wetlands
     &  STAT=AS)
 
      END SUBROUTINE SET_CH3I_MOD
#endif
!------------------------------------------------------------------------------

      END MODULE CH3I_MOD
