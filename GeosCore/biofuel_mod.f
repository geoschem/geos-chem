! $Id: biofuel_mod.f,v 1.1 2009/09/16 14:06:40 bmy Exp $
      MODULE BIOFUEL_MOD
!
!******************************************************************************
!  Module BIOFUEL_MOD contains arrays and routines to compute yearly
!  biofuel emissions for NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, 
!  CH2O, and C2H6. (bmy, 9/12/00, 9/18/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NBFMAX             : Maximum # of biofuel burning species
!  (2 ) NBFTRACE           : # of emitted biofuel burning species (<= NBFMAX)
!  (3 ) BFTRACE            : Array of tracer #'s for emitted biofuel species
!  (4 ) BIOFUEL            : array containing biofuel emissions
!
!  Module Routines:
!  ============================================================================
!  (1 ) BIOFUEL_BURN       : Reads data from disk & computes biofuel emissions
!  (2 ) SCALE_BIOFUEL_CO   : Scales biofuel CO to account for VOC oxidation
!  (3 ) SCALE_BIOFUEL_ACET : Scales biofuel ACET to match a posteriori source
!  (4 ) SCALE_FUTURE       : Applies future scale factors to biofuel emissions
!  (5 ) SET_BFTRACE        : Initializes NBFTRACE counter and BFTRACE array
!  (6 ) INIT_BIOFUEL       : Initializes the BIOFUEL array
!  (7 ) CLEANUP_BIOFUEL    : Deallocates the BIOFUEL array
!
!  GEOS-CHEM modules referenced by biofuel_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f          : Module w/ DAO met field arrays
!  (3 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f    : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) epa_nei_mod.f      : Module w/ routines to read EPA/NEI99 data
!  (6 ) error_mod.f        : Module w/ NaN and other error check routines
!  (7 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (8 ) tracer_mod.f       : Module w/ GEOS-CHEM tracer array etc.
!  (9 ) tracerid_mod.f     : Module w/ pointers to tracers & emissions 
!  (10) transfer_mod.f     : Module w/ routines to cast & resize arrays
!
!  References:
!  ============================================================================
!
!  (1 ) Andreae, M.O., and P. Merlet, "Emissions of trace gases and aerosols
!        from biomass burning", Global Biogeochemical Cycles, Vol 15, pp
!        955-966, 2001.
!  (2 ) Hays, M.D., C.D. Geron, K.J. Linna, N.D. Smith, and J.J. Schauer, 
!        "Speciation of gas-phase and fine particle emissions from burning of
!        foliar fuels", Environ. Sci. Technol., Vol 36, pp 2281-2295, 2002.
!
!  NOTES:
!  (1 ) Now account for extra production of CO from VOC's for Tagged CO
!        and CO-OH simulations (bmy, 1/3/01)
!  (2 ) Now read NBIOFUEL=10 biofuel species.  Also archive biofuel emissions
!        in the ND34 diagnostic. (bmy, 4/17/01)
!  (3 ) Now dimension BFTRACE arrays to be of size NBFMAX instead of having 
!        them be made allocatable.  Also updated comments. (bmy, 4/17/01)
!  (4 ) Bug fix: now make sure to index biofuel tracers w/ the correct tracer
!        number, even when there are less than the maximum species being
!        requested (bmy, 8/24/01)
!  (5 ) Bug fix: now index biofuel CH2O correctly (bmy, 8/28/01)
!  (6 ) Now scale biofuel ACET by 0.82, in order to match the a posteriori
!        acetone source from Jacob et al 2001.  Also updated comments.
!        (bdf, bmy, 9/10/01)
!  (7 ) BIOFUEL is now declared (NBFTRACE,IIPAR,JJPAR).  Now use TRANSFER_2D
!        from "transfer_mod.f" to copy data into BIOFUEL. (bmy, 9/28/01)
!  (8 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (9 ) Now do unit conversion every time step.  Also added private
!        array BIOFUEL_KG to hold emissions in kg over the entire
!        month. (bmy, 5/9/02)
!  (10) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments.  BIOMASS_KG is now
!        an allocatable module array instead of a local array in routine
!        "biofuel_burn.f". (bmy, 5/28/02)  
!  (11) Now reference BXHEIGHT from "dao_mod.f". Now references "error_mod.f".
!        Also deleted obsolete code from various routines.  Also references
!        "tracerid_mod.f"  Added routine SET_NBFTRACE. (bmy, 11/6/02)
!  (12) Now call READ_BPCH2 with QUIET=.TRUE. to suppress output (bmy, 3/14/03)
!  (13) Now references "directory_mod.f" (bmy, 7/19/04)
!  (14) Now references "time_mod.f" and "epa_nei_mod.f" (bmy, 11/5/04)
!  (15) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (16) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (17) Rewrite IF statements to avoid seg fault errors when LNEI99 is turned 
!        off (bmy, 2/1/06)
!  (18) Modified for IPCC future emissions scale factors.  Added private
!        routine SCALE_FUTURE. (swu, bmy, 5/30/06)
!  (19) Modified for VOC-scaling of CO emissions for H2/HD sim (phs, 5/16/07)
!  (20) Added 9 gaseous biofuel emissions: GLYX, MGLY, BENZ, 
!        TOLU, XYLE, C2H4, C2H2, GLYC, HAC. (tmf, 1/7/09)
!  (21) Emissions for these 9 tracers are scaled from CO emissions. (tmf, 1/7/09)
!  (22) Updated for Havala's SOA + semivol POA code. Added for gas phase NAP
!       chemistry and NAP biofuel emissions. (mpayer, 7/6/11)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biofuel_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      PUBLIC :: NBFMAX
      PUBLIC :: NBFTRACE
      PUBLIC :: BFTRACE
      PUBLIC :: BIOFUEL

      ! ... and these routines
      PUBLIC :: BIOFUEL_BURN
      PUBLIC :: CLEANUP_BIOFUEL
      PUBLIC :: INIT_BIOFUEL
      PUBLIC :: SET_BFTRACE

      !=================================================================     
      ! MODULE VARIABLES
      !=================================================================

      !-----------------------------------------------------------------
      ! Prior to 7/6/11:
      ! Increase NMFMAX from 19 to 20 for NAP (hotp, mpayer, 7/6/11)
      ! INTEGER, PARAMETER  :: NBFMAX = 19
      !-----------------------------------------------------------------
      INTEGER, PARAMETER  :: NBFMAX = 20

      INTEGER             :: NBFTRACE
      INTEGER             :: BFTRACE(NBFMAX) 

      REAL*8, ALLOCATABLE :: BIOFUEL(:,:,:)
      REAL*8, ALLOCATABLE :: BIOFUEL_KG(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE BIOFUEL_BURN
!
!******************************************************************************
!  Subroutine BIOFUEL_BURN computes the yearly biofuel burning emissions
!  and also archives them into GEOS-CHEM diagnostics. 
!  (rvm, acs, bnd, bmy, 9/12/00, 5/30/06)
!
!  Biofuel emissions are based on estimates by Rose Yevich and Jennifer
!  Logan (reference TBA).
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Hayes et al, 2002
!  (2 ) Andreae & Merlet, 2001
!
!  NOTES:
!  (1 ) Renamed array that held biofuel emissions from TWOODIJ to BIOFUEL,
!        and also expaneded to hold NOx. (rvm, acs, bmy, 9/12/00)
!  (2 ) BIOFUEL is a true global array -- use offsets IREF, JREF to index
!        it as you would T, P, or other global CTM variables (bmy, 9/12/00)
!  (3 ) ND29 and ND32 diagnostics are now computed within BIOFUEL_BURN
!        instead of in "emissdr.f", as was done previously. (bmy, 9/12/00)
!  (4 ) Enhance CO from biofuel burning by 10% for Tagged CO and CO-OH
!        simulations, to account for extra production of CO from VOC's.
!        (bnd, bmy, 1/3/01)  
!  (5 ) Now read NBIOFUEL=10 biofuel species.  Also archive biofuel emissions
!        in the ND34 diagnostic.  Updated output information. (bmy, 4/17/01)
!  (6 ) Now read new biofuel emissions (Apr 2001) from the "biofuel_200104"
!        subdirectory of DATA_DIR (bmy, 4/18/01)
!  (7 ) Bug fix: now make sure to index biofuel tracers w/ the correct tracer
!        number, even when there are less than the maximum species being
!        requested (bmy, 8/24/01)
!  (8 ) Bug fix: now use tracer #20 to read biofuel CH2O (bmy, 8/28/01)
!  (9 ) Now scale biofuel ACET by 0.5, in order to match the a posteriori
!        acetone source from Jacob et al 2001 (bdf, bmy, 9/5/01)
!  (10) Remove IREF, JREF -- they are obsolete.  BIOFUEL(:,IREF,JREF) is now 
!        BIOFUEL(:,I,J).  Make sure to use IDBFCO and IDBFNOX in ND29 and
!        ND32 diagnostics.  Now use TRANSFER_2D from "transfer_mod.f" to
!        cast data from REAL*4 to REAL*8 and copy to BIOFUEL (bmy, 9/28/01)
!  (11) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (12) Bug fix -- need to convert from kg --> molec/cm3/s on every time
!        step since the box volumes change w/ the surface pressure over
!        the course of the year.  Add parallel DO-loop for unit conversion.
!        Also archive diagnostics w/in the parallel DO-loop.  MOLWT needs to 
!        be an array of size (NBFMAX).  Now read biofuel file w/ the correct 
!        amt of Tg C for ACET, C2H6, C3H8. (bmy, 6/11/02)
!  (13) Now reference BXHEIGHT from "dao_mod.f".  Also reference IDTNOX, 
!        IDBFNOX, etc. from "tracerid_mod.f". (bmy, 11/6/02)
!  (14) Now call READ_BPCH2 with QUIET=.TRUE. flag to suppress extra info 
!        from being printed (bmy, 3/14/03)
!  (15) Added fancy output (bmy, 4/26/04)
!  (16) Now references "tracer_mod.f" and "directory_mod.f" (bmy, 7/19/04)
!  (17) Now can overwrite USA with EPA/NEI biofuel emissions.  Now references 
!        function GET_DAY_OF_WEEK from "time_mod.f".  Now references LNEI99
!        from "logical_mod.f".  Now reference GET_EPA_BIOFUEL and
!        GET_USA_MASK from "epa_nei_mod.f". (rch, rjp, bmy, 11/5/04)
!  (18) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (19) Now make sure all USE statements are USE, ONLY.  Eliminate reference 
!        to TRACER_MOD, it's obsolete (bmy, 10/3/05)
!  (20) Rewrite IF statements to avoid seg fault errors when LNEI99 is turned 
!        off (bmy, 2/1/06)
!  (21) Now references LFUTURE from "logical_mod.f".  
!  (22) Now reference ITS_A_H2HD_SIM from "tracer_mod.f" for ND29.
!        (phs, 9/18/07)
!  (23) Switch off biofuel in S.E.-Asia if Streets 2006 inventory is used,
!        accounting for FSCLYR from CMN_O3 (phs,3/17/08)
!  (24) Add scaling of aromatic emissions over the US. (hotp, 11/23/09)
!  (25) Add IDTNAP, NAMEMISS, and IDBFNAP (hotp, mpayer, 7/6/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D,  GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,         READ_BPCH2
      USE DAO_MOD,              ONLY : BXHEIGHT
      USE DIAG_MOD,             ONLY : AD29,    AD32_bf, AD34
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE EPA_NEI_MOD,          ONLY : GET_EPA_BIOFUEL,  GET_USA_MASK
      USE LOGICAL_MOD,          ONLY : LFUTURE, LNEI99,  LSTREETS
      USE STREETS_ANTHRO_MOD,   ONLY : GET_SE_ASIA_MASK
      USE TIME_MOD,             ONLY : GET_DAY_OF_WEEK,  GET_YEAR
      USE TRACER_MOD,           ONLY : ITS_A_H2HD_SIM
      USE TRACERID_MOD,         ONLY : IDBFCO,  IDBFNOX, IDTACET
      USE TRACERID_MOD,         ONLY : IDTALD2, IDTALK4, IDTC2H6
      USE TRACERID_MOD,         ONLY : IDTC3H8, IDTCH2O, IDTCO
      USE TRACERID_MOD,         ONLY : IDTMEK,  IDTNOX,  IDTPRPE 
      USE TRACERID_MOD,         ONLY : IDTGLYX, IDTMGLY, IDTBENZ
      USE TRACERID_MOD,         ONLY : IDTTOLU, IDTXYLE, IDTC2H4
      USE TRACERID_MOD,         ONLY : IDTC2H2, IDTGLYC, IDTHAC
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      ! for US emission fix (hotp 11/20/09)
      USE TRACERID_MOD,         ONLY : IDBFBENZ,IDBFTOLU,IDBFXYLE
      USE TRACERID_MOD,         ONLY : IDBFGLYX,IDBFMGLY,IDBFC2H4
      USE TRACERID_MOD,         ONLY : IDBFC2H2,IDBFGLYC,IDBFHAC
      ! for gas phase NAP chemistry, NAP biofuel emiss (hotp, mpayer, 7/6/11)
      USE TRACERID_MOD,         ONLY : IDTNAP
      USE TRACER_MOD,           ONLY : ITS_A_FULLCHEM_SIM 
      USE TRACER_MOD,           ONLY : ITS_A_TAGCO_SIM
      USE LOGICAL_MOD,          ONLY : NAPEMISS
      USE TRACERID_MOD,         ONLY : IDBFNAP

      IMPLICIT NONE
      
#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_DIAG"             ! ND29, ND32, ND34
#     include "CMN_O3"               ! FSCLYR
      
      ! Local variables
      LOGICAL                       :: WEEKDAY
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      LOGICAL, SAVE                 :: DO_ND29, DO_ND32, DO_ND34
      INTEGER                       :: AS, I, J, N, NN, DAY_NUM
      INTEGER                       :: SIM_YEAR
      REAL*4                        :: ARRAY(IGLOB,JGLOB,1)
      REAL*8,  SAVE                 :: MOLWT(NBFMAX)
      REAL*8                        :: TOTAL, BXHEIGHT_CM, EPA_NEI
      REAL*8                        :: COSCALEDOWN  ! (hotp, mpayer, 7/6/11)
      CHARACTER(LEN=255)            :: FILENAME 
      
      ! External functions
      REAL*8,  EXTERNAL             :: BOXVL

      REAL*8                        :: BF_CO( IIPAR, JJPAR )  ! Biofuel emission of CO [molec/cm2/s]
      
      ! Scale up NAP emiss (hotp, mpayer, 7/6/11)
      REAL*8, PARAMETER             :: NAPTOTALSCALE = 66.09027d0

      !=================================================================
      !   B i o f u e l   B u r n i n g   B e g i n s   H e r e !!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate and zero the BIOFUEL array 
         CALL INIT_BIOFUEL

         ! Flags for whether or not to archive diagnostics
         DO_ND29 = ( ND29 > 0 .and. IDBFCO  /= 0 ) 
         DO_ND32 = ( ND32 > 0 .and. IDBFNOX /= 0 ) 
         DO_ND34 = ( ND34 > 0                    ) 

         ! Fancy output..
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)' ) 'B I O F U E L   E M I S S I O N S'

         !==============================================================
         ! GEOS-CHEM has the following biofuel burning species:
         !
         !    Species   Index   CTM Tracer #   Units as read from file
         !    ---------------------------------------------------------
         !      NOX       1          1         [kg NOx /box/year]
         !      CO        2          4         [kg CO  /box/year]
         !      ALK4      3          5         [kg C   /box/year]
         !      ACET      4          9         [kg C   /box/year]
         !      MEK       5          10        [kg C   /box/year]
         !      ALD2      6          11        [kg C   /box/year]
         !      PRPE      7          18        [kg C   /box/year]
         !      C3H8      8          19        [kg C   /box/year]
         !      CH2O      9          20        [kg CH2O/box/year]
         !      C2H6      10         21        [kg C   /box/year]
         !      GLYX                           [kg     /box/year]
         !      MGLY                           [kg     /box/year]
         !      BENZ                           [kg C   /box/year]
         !      TOLU                           [kg C   /box/year]
         !      XYLE                           [kg C   /box/year]
         !      C2H4                           [kg C   /box/year]
         !      C2H2                           [kg C   /box/year]
         !      GLYC                           [kg     /box/year]
         !      HAC                            [kg     /box/year]
         !      ! (hotp, mpayer, 7/6/11)
         !      NAP                            [kg C   /box/year]
         !
         ! These emissions are converted to [molec/cm3/s] (or 
         ! [molec C/cm3/s] for hydrocarbons), since the chemistry
         ! requires these units.
         !
         ! There are NBFMAX=20 maximum allowed biofuel species, but 
         ! only NBFTRACE of these are actually emitted.  Species are 
         ! turned off/on with the switches in the "tracer.dat" file.
         !
         ! The BIOFUEL array is only of size NBFTRACE, to save memory.
         ! We only read in the NBFTRACE species that are emitted.
         !
         ! Biofuel burning emissions are aseasonal, so we only have
         ! to read from disk on the very first model timestep. 
         ! However, we have to convert from kg --> molec/cm3/s every
         ! timestep to ensure that we use the box heights and box
         ! volumes throughout the year, instead of only at the
         ! first timestep (bmy, 5/3/02)
         !==============================================================
         FILENAME = TRIM( DATA_DIR )          // 
     &              'biofuel_200202/biofuel.' // GET_NAME_EXT_2D() // 
     &              '.'                       // GET_RES_EXT()

         ! Echo filename to log file
         WRITE( 6, 110 ) TRIM( FILENAME )

         ! Loop over the emitted biofuel burning tracers only
         DO N = 1, NBFTRACE

            ! NN is the actual CTM tracer number of species N
            NN = BFTRACE(N)

            ! Test for each tracer
            IF ( NN == IDTNOX ) THEN

               !----------------
               ! Biofuel NOx
               !----------------

               ! Read biofuel NOx emissions in [kg/box/yr] -- tracer #1
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 1 , 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 
                         
               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )
               
               ! Compute future NOx emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'NOxbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total of biofuel NOx
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'NOx ', TOTAL, '[Tg  /yr]'
               
               ! Define MOLWT for use below
               MOLWT(N) = 46d-3

            ELSE IF ( NN == IDTCO ) THEN

               !----------------
               ! Biofuel CO
               !----------------

               ! Read biofuel CO emissions in [kg/box/yr] -- tracer #4
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 4, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Save BF_CO before scaling (tmf, 6/15/07) 
               CALL TRANSFER_2D( ARRAY(:,:,1), BF_CO(:,:) )         

               ! Scale CO to account for oxidation of extra VOC's
               CALL SCALE_BIOFUEL_CO( BIOFUEL_KG(N,:,:) )
               
               ! Compute future CO emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'CObf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Print total CO
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'CO  ', TOTAL, '[Tg  /yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 28d-3

            ELSE IF ( NN == IDTALK4 ) THEN

               !----------------
               ! Biofuel ALK4
               !----------------

               ! Read biofuel ALK4 emissions in [kg/box/yr] -- tracer #5
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 5, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future ALK4 emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total ALK4
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ALK4', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTACET ) THEN

               !----------------
               ! Biofuel ACET
               !----------------

               ! Read biofuel ACET emissions in [kg/box/yr] -- tracer #9
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 9, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Scale to match a posteriori source (bdf, bmy, 9/10/01)
               CALL SCALE_BIOFUEL_ACET( BIOFUEL_KG(N,:,:) )

               ! Compute future ACET emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total ACET
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ACET', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTMEK ) THEN

               !----------------
               ! Biofuel MEK
               !----------------

               ! Read biofuel MEK emissions in [kg/box/yr] -- tracer #10
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 10, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future MEK emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total MEK
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'MEK ', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTALD2 ) THEN

               !----------------
               ! Biofuel ALD2
               !----------------

               ! Read biofuel ALD2 emissions in [kg/box/yr] -- tracer #11
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 11, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future ALD2 emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total ALD2
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ALD2', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTPRPE ) THEN

               !----------------
               ! Biofuel PRPE
               !----------------

               ! Read biofuel PRPE emissions in [kg/box/yr] -- tracer #18
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 18, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future PRPE emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total PRPE
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'PRPE', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTC3H8 ) THEN

               !----------------
               ! Biofuel C3H8
               !----------------

               ! Read biofuel C3H8 emissions in [kg/box/yr] -- tracer #19
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 19, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future C3H8 emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total C3H8
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C3H8', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTCH2O ) THEN

               !----------------
               ! Biofuel CH2O
               !----------------

               ! Read biofuel CH2O emissions in [kg/box/yr] -- tracer #20
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 20, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future CH2O emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF

               ! Compute total CH2O
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'CH2O', TOTAL, '[Tg  /yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 30d-3

            ELSE IF ( NN == IDTC2H6 ) THEN

               !----------------
               ! Biofuel C2H6
               !----------------

               ! Read biofuel C2H6 emissions in [kg/box/yr] -- tracer #21
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 21, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute future C2H6 emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'VOCbf', BIOFUEL_KG(N,:,:) )
               ENDIF 

               ! Compute total C2H6
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C2H6', TOTAL, '[Tg C/yr]'
               
               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTGLYX ) THEN

               !----------------
               ! Biofuel GLYX
               !----------------

               ! Emission ratio GLYX/CO = 6.62d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 58d-3 * 6.62d-3      ! [kg/box/yr]

               ! Compute total GLYX
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'GLYX', TOTAL, '[Tg/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 58d-3

 
            ELSE IF ( NN == IDTMGLY ) THEN

               !----------------
               ! Biofuel MGLY
               !----------------

               ! Emission ratio MGLY/CO = 3.47d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 72d-3 * 3.47d-3      ! [kg/box/yr]

               ! Compute total MGLY
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'MGLY', TOTAL, '[Tg/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 72d-3

            ELSE IF ( NN == IDTBENZ ) THEN

               !----------------
               ! Biofuel BENZ
               !----------------

               ! Emission ratio BENZ/CO = 4.06d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 12d-3 * 6d0 * 4.06d-3      ! [kg C/box/yr]

               ! Compute total BENZ
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'BENZ', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3


            ELSE IF ( NN == IDTTOLU ) THEN

               !----------------
               ! Biofuel TOLU
               !----------------

               ! Emission ratio TOLU/CO = 2.01d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 12d-3 * 7d0 * 2.01d-3      ! [kg C/box/yr]

               ! Compute total TOLU
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'TOLU', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTXYLE ) THEN

               !----------------
               ! Biofuel XYLE
               !----------------

               ! Emission ratio XYLE/CO = 0.82d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 12d-3 * 8d0 * 0.82d-3      ! [kg C/box/yr]

               ! Compute total XYLE
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'XYLE', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTC2H4 ) THEN

               !----------------
               ! Biofuel C2H4
               !----------------

               ! Emission ratio C2H4/CO = 15.7d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 12d-3 * 2d0 * 15.7d-3      ! [kg C/box/yr]

               ! Compute total C2H4
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C2H4', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTC2H2 ) THEN

               !----------------
               ! Biofuel C2H2
               !----------------

               ! Emission ratio C2H2/CO = 19d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 12d-3 * 2d0 * 19d-3      ! [kg C/box/yr]

               ! Compute total C2H2
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C2H2', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTGLYC ) THEN

               !----------------
               ! Biofuel GLYC
               !----------------

               ! Emission ratio GLYC/CO = 3.66d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 60d-3 * 3.66d-3      ! [kg/box/yr]

               ! Compute total GLYC
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'GLYC', TOTAL, '[Tg/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 60d-3

            ELSE IF ( NN == IDTHAC  ) THEN

               !----------------
               ! Biofuel HAC
               !----------------

               ! Emission ratio HAC/CO = 3.31d-3 [mole/mole]
               BIOFUEL_KG(N,:,:) = 
     &          BF_CO(:,:) / 28d-3 * 74d-3 * 3.31d-3      ! [kg/box/yr]

               ! Compute total HAC
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'HAC', TOTAL, '[Tg/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 74d-3

            ELSE IF ( NN == IDTNAP ) THEN

               !----------------
               ! Biofuel NAP (hotp, mpayer, 7/6/11)
               !----------------

               ! Scale down CO to account for VOC oxidation
               COSCALEDOWN = 1d0   ! default value

               ! Scale based on simulation type
               IF ( ITS_A_FULLCHEM_SIM() ) THEN
                   COSCALEDOWN = 1d0/1.086d0
               ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
                   COSCALEDOWN = 1d0/1.189d0
               ENDIF

               ! Emmision ratio NAP/CO = 0.0701d-3 [mole/mole]
               ! NAP emiss =  0.025 g NAP/kg DM (Table 4, Hays et al, 2002)
               ! CO  emiss =     78 g CO /kg DM (Table 1,Andreae & Merlet,2001)
               ! Scale emissions down if appropriate using COSCALEDOWN
               BIOFUEL_KG(N,:,:) = BIOFUEL_KG(IDBFCO,:,:) * 0.0701d-3 
     &                             * 120d0 / 28d0 * COSCALEDOWN ! [kg C/box/yr]

               ! Scale up total
               BIOFUEL_KG(N,:,:) = BIOFUEL_KG(N,:,:) * NAPTOTALSCALE

               ! Set NAP emissions according to input.geos
               BIOFUEL_KG(N,:,:) = BIOFUEL_KG(N,:,:) * NAPEMISS

               ! Compute future NAP emissions (if necessary)
               IF ( LFUTURE ) THEN
                  CALL SCALE_FUTURE( 'CO', BIOFUEL_KG(N,:,:) )
               ENDIF 

               ! Compute total NAP
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'NAP ', TOTAL, '[Tg C/yr]'
               
               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ENDIF
         ENDDO

         ! Reset first time flag
         FIRST = .FALSE.

         ! Fancy output
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF   

      !=================================================================
      ! Do the following on each emission timestep...
      ! 
      ! Convert from [kg/box/yr] (or [kg C/box/yr]) to [molec/cm3/s] 
      ! (or [molec C/cm3/s]), since the emissions need to be in these 
      ! units for the chemistry.  Now use parallel DO loops. 
      !
      ! NOTE: Need to do the unit conversion outside the IF (FIRST)
      ! block, so that we use the same airmass quantities as are used
      ! for the diagnostics. (bmy, 5/30/02)
      !
      ! Also archive diagnostics w/in parallel loop (bmy, 5/30/02)
      !=================================================================

      ! Get current day of the week
      DAY_NUM = GET_DAY_OF_WEEK()

      ! Is it a weekday?
      WEEKDAY = ( DAY_NUM > 0 .and. DAY_NUM < 6 )

      ! get emissions year to test Streets
      IF ( FSCALYR < 0 ) THEN
         SIM_YEAR = GET_YEAR()
      ELSE
         SIM_YEAR = FSCALYR
      ENDIF



!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, BXHEIGHT_CM, N, NN, EPA_NEI )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! BXHEIGHT_CM = the surface grid box height in cm
         BXHEIGHT_CM = BXHEIGHT(I,J,1) * 1d2

         ! Loop over all biofuel tracers
         DO N = 1, NBFTRACE

            ! Biofuel emissions in [molec/cm3/s]
            BIOFUEL(N,I,J) = BIOFUEL_KG(N,I,J)                  * 
     &                       ( 6.023d23 / MOLWT(N)            ) /
     &                       ( 365d0 * 86400d0 * BOXVL(I,J,1) )


            !-----------------------------------------------------------
            ! Overwrite biofuels w/ EPA/NEI emissions over the USA
            !-----------------------------------------------------------
            
            ! If EPA/NEI99 emissions are turned on....
            IF ( LNEI99 ) THEN

               ! If we are over the USA ...
               IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
  
                  ! Get GEOS-CHEM tracer number
                  NN      = BFTRACE(N)

                  ! We do not have EPA/NEI biofuel emission.  
                  ! Use default emission for the newly added species. 
                  ! (tmf, 1/8/08) 
                  ! Add NAP (hotp, mpayer, 7/6/11)
                  IF ( (NN /= IDTGLYX) .and. (NN /= IDTMGLY) .and. 
     &                 (NN /= IDTBENZ) .and. (NN /= IDTTOLU) .and.
     &                 (NN /= IDTXYLE) .and. (NN /= IDTC2H4) .and.
     &                 (NN /= IDTC2H2) .and. (NN /= IDTGLYC) .and.             
     &                 (NN /= IDTHAC ) .and. (NN /= IDTNAP ) ) THEN

                     ! Get EPA/NEI biofuel [molec/cm2/s or atoms C/cm2/s]
                     EPA_NEI = GET_EPA_BIOFUEL( I, J, NN, WEEKDAY )

                     ! Convert [molec/cm2/s] to [molec/cm3/s]
                     BIOFUEL(N,I,J) = EPA_NEI / BXHEIGHT_CM
                  
                  ENDIF
               ENDIF
            ENDIF

            !-----------------------------------------------------------
            ! If we are over SE ASIA and are using Streets 2006 (that is
            ! emission year is GE 2006), set BIOFUEL to zero since they
            ! are already accounted for (phs, 3/17/08)
            !-----------------------------------------------------------
            IF ( LSTREETS .and. ( SIM_YEAR >= 2006 ) ) THEN

               ! If we are over the SE Asia region
               IF ( GET_SE_ASIA_MASK( I, J ) > 0d0 ) THEN

                  BIOFUEL(N,I,J) = 0.d0

               ENDIF
            ENDIF

            ! move below so that BIOFUEL array is complete before
            ! archiving diagnostic info (hotp 11/23/09)
!            ! ND34 -- archive biofuel burning species [molec/cm2/s]
!            IF ( DO_ND34 ) THEN
!               AD34(I,J,N) = AD34(I,J,N) + ( BIOFUEL(N,I,J) * 
!     &                                       BXHEIGHT_CM )
!            ENDIF
         ENDDO  

         ! ND29 -- CO source diagnostics [molec/cm2/s]
         IF ( DO_ND29 ) THEN

            IF ( ITS_A_H2HD_SIM() ) THEN
               AD29(I,J,3) = AD29(I,J,3) + ( BIOFUEL(IDBFCO,I,J) * 
     &                                       BXHEIGHT_CM ) * 1.189d0
            ELSE
               AD29(I,J,3) = AD29(I,J,3) + ( BIOFUEL(IDBFCO,I,J) * 
     &                                       BXHEIGHT_CM )
            ENDIF

         ENDIF

         ! ND32 -- NOx source diagnostics [molec/cm2/s]
         IF ( DO_ND32 ) THEN 
            AD32_bf(I,J) = AD32_bf(I,J) + ( BIOFUEL(IDBFNOX,I,J) * 
     &                                      BXHEIGHT_CM )
         ENDIF
      ENDDO  
      ENDDO
!$OMP END PARALLEL DO

      ! update aromatics based on CO over US (hotp 11/23/09)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            ! If EPA/NEI99 emissions are turned on....
            IF ( LNEI99 ) THEN

               ! If we are over the USA ...
               IF ( GET_USA_MASK( I, J ) > 0d0 ) THEN
                  ! update aromatics if necessary (hotp 11/20/09)
                  ! molecC/cm3/s XXXX = molec/cm3/s CO * 
                  !                     mol XXXX/mol CO *
                  !                     molec C/molec XXXX

                  ! BENZ (6 carbon/molec)
                  IF ( IDBFBENZ > 0 ) THEN
                     BIOFUEL(IDBFBENZ,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 4.06d-3 * 6.d0
                  ENDIF

                  ! TOLU (7 carbon/molec)
                  IF ( IDBFTOLU > 0 ) THEN
                     BIOFUEL(IDBFTOLU,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 2.01d-3 * 7.d0
                  ENDIF

                  ! XYLE (8 carbon/molec)
                  IF ( IDBFXYLE > 0 ) THEN
                     BIOFUEL(IDBFXYLE,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 0.82d-3 * 8.d0
                  ENDIF

                  ! GLYX (molec/cm3/s)
                  IF ( IDBFGLYX > 0 ) THEN
                     BIOFUEL(IDBFGLYX,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 6.62d-3 
                  ENDIF

                  ! MGLY
                  IF ( IDBFMGLY > 0 ) THEN
                     BIOFUEL(IDBFMGLY,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 3.47d-3 
                  ENDIF

                  ! C2H4 (2 carbons/molec)
                  IF ( IDBFC2H4 > 0 ) THEN
                     BIOFUEL(IDBFC2H4,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 15.7d-3 * 2.d0
                  ENDIF

                  ! C2H2 (2 carbons/molec)
                  IF ( IDBFC2H2 > 0 ) THEN
                     BIOFUEL(IDBFC2H2,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 19.0d-3 * 2.d0
                  ENDIF

                  ! GLYC
                  IF ( IDBFGLYC > 0 ) THEN
                     BIOFUEL(IDBFGLYC,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 3.66d-3
                  ENDIF

                  ! HAC
                  IF ( IDBFHAC > 0 ) THEN
                     BIOFUEL(IDBFHAC,I,J) = 
     &                  BIOFUEL(IDBFCO,I,J) * 3.31d-3             
                  ENDIF

                  ! NAP (hotp, mpayer, 7/6/11)
                  IF ( IDBFNAP > 0 ) THEN
                  BIOFUEL(IDBFNAP,I,J) = 
     &                 BIOFUEL(IDBFCO,I,J) * 0.0701d-3 * 10.d0 *
     &                  NAPTOTALSCALE * NAPEMISS
                  ENDIF

               ENDIF ! USA_MASK

            ENDIF ! LNEI99
      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO


      ! Save diagnostic information (biofuel burning emissions) (hotp 11/23/09)
      ! ND34 -- archive biofuel burning species [molec/cm2/s] or
      ! [molecC/cm2/s] for species transported as carbon
      IF ( DO_ND34 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, BXHEIGHT_CM, N )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

            ! BXHEIGHT_CM = the surface grid box height in cm
            BXHEIGHT_CM = BXHEIGHT(I,J,1) * 1d2

            ! loop over biofuel species
            DO N = 1, NBFTRACE
               ! ND34 -- archive biofuel burning species [molec/cm2/s]
               AD34(I,J,N) = AD34(I,J,N) + ( BIOFUEL(N,I,J) * 
     &                                       BXHEIGHT_CM )
            ENDDO ! NBFTRACE

      ENDDO ! I
      ENDDO ! J
!$OMP END PARALLEL DO
      ENDIF ! ND34

      !=================================================================
      ! FORMAT statements
      !=================================================================
 110  FORMAT( 'BIOFUEL_BURN:  Reading ', a )
 120  FORMAT( 'Sum Biofuel ', a4, 1x, ': ', f9.3, 1x, a9  )

      ! Return to calling program
      END SUBROUTINE BIOFUEL_BURN

!------------------------------------------------------------------------------

      SUBROUTINE SCALE_BIOFUEL_CO( BFARRAY )
!
!******************************************************************************
!  Subroutine SCALE_BIOFUEL_CO multiplies the CO biofuel emissions by scale
!  factors to account for CO production from VOC's that are not explicitly 
!  carried in the chemistry mechanisms. (bnd, bmy, 3/19/01, 7/19/04)
!  
!  Arguments as Input:
!  ============================================================================
!  (1) BFARRAY (REAL*8) : Array containing biofuel burning CO emissions
!
!  NOTES:
!  (1 ) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!       Bryan Duncan (bnd@io.harvard.edu) and Daniel Jacob (djj@io.harvard.edu)
!  (2 ) BFARRAY is now of size (IIPAR,JJPAR) (bmy, 9/28/01)
!  (3 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (4 ) Now use inquiry functions in "tracer_mod.f" instead of the variable
!        NSRCX (bmy, 7/19/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGCO_SIM

#     include "CMN_SIZE"    ! Size parameters
      
      ! Arguments
      REAL*8, INTENT(INOUT) :: BFARRAY(IIPAR,JJPAR) 

      !=================================================================
      ! SCALE_BIOFUEL_CO begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN 

         ! Full chemistry w/ SMVGEAR  -- enhance by 8.6%
         BFARRAY = BFARRAY * 1.086d0

      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN 

         ! Tagged CO -- enhance by 18.9%
         BFARRAY = BFARRAY * 1.189d0

      ENDIF

      ! Return to calling program  
      END SUBROUTINE SCALE_BIOFUEL_CO

!------------------------------------------------------------------------------

      SUBROUTINE SCALE_BIOFUEL_ACET( BFARRAY )
!
!******************************************************************************
!  Subroutine SCALE_BIOFUEL_ACET multiplies the ACET biofuel emissions by a
!  scale factor in order to match the source from the Jacob et al 2001 paper.
!  (bdf, bmy, 9/10/01, 11/15/01)
!  
!  Arguments as Input:
!  ============================================================================
!  (1) BFARRAY (REAL*8) : Array containing biofuel burning ACET emissions
!
!  Reference:
!  ============================================================================
!  Jacob, D.J., B.D. Field, E. Jin, I. Bey, Q. Li, J.A. Logan, and 
!    R.M. Yantosca, Atmospheric budget of acetone, submitted to 
!    Geophys. Res. Lett., 2001. 
!
!  NOTES:
!  (1 ) Adapted from SCALE_BIOMASS_CO (bdf, bmy, 9/10/01)
!  (2 ) BFARRAY is now of size (IIPAR,JJPAR) (bmy, 9/28/01)
!  (3 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!******************************************************************************
!
#     include "CMN_SIZE"    ! Size parameters
      
      ! Arguments
      REAL*8, INTENT(INOUT) :: BFARRAY(IIPAR,JJPAR) 

      !=================================================================
      ! SCALE_BIOFUEL_ACET begins here!
      !=================================================================

      ! Scale by 0.82 to match the a posteriori source
      BFARRAY = BFARRAY * 0.82d0
        
      ! Return to calling program  
      END SUBROUTINE SCALE_BIOFUEL_ACET

!------------------------------------------------------------------------------

      SUBROUTINE SCALE_FUTURE( NAME, BF )
!
!******************************************************************************
!  Subroutine SCALE_FUTURE applies the IPCC future emissions scale factors
!  to the biofuel emisisons in order to compute the future biofuel emissions
!  for NOx, CO, and VOC's (swu, bmy, 5/30/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NAME (CHARACTER) : Denotes type of scale factor to use (e.g. NOx)
!  (2 ) BF   (REAL*8   ) : Array w/ biomass burning emisisons [molec/cm2]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObf
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbf
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbf

#     include "CMN_SIZE"               ! Size parameters

      ! Arguments
      REAL*8,           INTENT(INOUT) :: BF(IIPAR,JJPAR)
      CHARACTER(LEN=*), INTENT(IN)    :: NAME

      ! Local variables
      INTEGER                         :: I, J
      
      !=================================================================
      ! SCALE_FUTURE begins here!
      !=================================================================

      IF ( NAME == 'NOxbf' ) THEN

         ! Compute future NOx emissions
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            BF(I,J) = BF(I,J) * GET_FUTURE_SCALE_NOxbf( I, J )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE IF ( NAME == 'CObf' ) THEN

         ! Compute future CO emissions 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            BF(I,J) = BF(I,J) * GET_FUTURE_SCALE_CObf( I, J )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
         
      ELSE

         ! Compute future hydrocarbon emissions
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            BF(I,J) = BF(I,J) * GET_FUTURE_SCALE_VOCbf( I, J )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE SET_BFTRACE
!
!******************************************************************************
!  Subroutine SET_NBFTRACE sets the NBFTRACE variable with the number of
!  biofuel tracers that are turned on.  This was split off from "tracerid.f"
!  in order to prevent circular module references. (bmy, 11/6/02, 10/3/05)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Add IDTNAP, IDBFNAP (mpayer, 7/6/11)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDBFACET, IDBFALD2, IDBFALK4, IDBFC2H6 
      USE TRACERID_MOD, ONLY : IDBFC3H8, IDBFCH2O, IDBFCO,   IDBFMEK  
      USE TRACERID_MOD, ONLY : IDBFNOX,  IDBFPRPE, IDTACET,  IDTALD2 
      USE TRACERID_MOD, ONLY : IDTALK4,  IDTC2H6,  IDTC3H8,  IDTCH2O
      USE TRACERID_MOD, ONLY : IDTCO,    IDTMEK,   IDTNOX,   IDTPRPE 

      USE TRACERID_MOD, ONLY : IDBFGLYX, IDBFMGLY, IDBFBENZ, IDBFTOLU
      USE TRACERID_MOD, ONLY : IDBFXYLE, IDBFC2H4, IDBFC2H2, IDBFGLYC
      USE TRACERID_MOD, ONLY : IDBFHAC
      USE TRACERID_MOD, ONLY : IDTGLYX,  IDTMGLY,  IDTBENZ,  IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE,  IDTC2H4,  IDTC2H2,  IDTGLYC
      USE TRACERID_MOD, ONLY : IDTHAC
      ! for gas phase NAP chemistry, NAP biofuel emiss (hotp, mpayer, 7/6/11)
      USE TRACERID_MOD, ONLY : IDTNAP
      USE TRACERID_MOD, ONLY : IDBFNAP

      !=================================================================
      ! SET_BFTRACE begins here!
      !=================================================================

      ! Initialize
      NBFTRACE = 0
      
      ! Increment NBFTRACE for each turned on biofuel tracer
      IF ( IDBFNOX  /= 0 ) NBFTRACE = NBFTRACE + 1
      IF ( IDBFCO   /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFALK4 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFACET /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFMEK  /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFALD2 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFPRPE /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC3H8 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFCH2O /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H6 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFGLYX /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFMGLY /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFBENZ /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFTOLU /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFXYLE /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H4 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H2 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFGLYC /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFHAC  /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFNAP  /= 0 ) NBFTRACE = NBFTRACE + 1  ! (hotp, mpayer, 7/6/11)

      ! Fill BFTRACE w/ appropriate TRACER ID #'s
      IF ( IDBFNOX  /= 0 ) BFTRACE(IDBFNOX ) = IDTNOX
      IF ( IDBFCO   /= 0 ) BFTRACE(IDBFCO  ) = IDTCO
      IF ( IDBFALK4 /= 0 ) BFTRACE(IDBFALK4) = IDTALK4
      IF ( IDBFACET /= 0 ) BFTRACE(IDBFACET) = IDTACET
      IF ( IDBFMEK  /= 0 ) BFTRACE(IDBFMEK ) = IDTMEK
      IF ( IDBFALD2 /= 0 ) BFTRACE(IDBFALD2) = IDTALD2
      IF ( IDBFPRPE /= 0 ) BFTRACE(IDBFPRPE) = IDTPRPE
      IF ( IDBFC3H8 /= 0 ) BFTRACE(IDBFC3H8) = IDTC3H8
      IF ( IDBFCH2O /= 0 ) BFTRACE(IDBFCH2O) = IDTCH2O  
      IF ( IDBFC2H6 /= 0 ) BFTRACE(IDBFC2H6) = IDTC2H6  
      IF ( IDBFGLYX /= 0 ) BFTRACE(IDBFGLYX) = IDTGLYX
      IF ( IDBFMGLY /= 0 ) BFTRACE(IDBFMGLY) = IDTMGLY
      IF ( IDBFBENZ /= 0 ) BFTRACE(IDBFBENZ) = IDTBENZ
      IF ( IDBFTOLU /= 0 ) BFTRACE(IDBFTOLU) = IDTTOLU
      IF ( IDBFXYLE /= 0 ) BFTRACE(IDBFXYLE) = IDTXYLE
      IF ( IDBFC2H4 /= 0 ) BFTRACE(IDBFC2H4) = IDTC2H4
      IF ( IDBFC2H2 /= 0 ) BFTRACE(IDBFC2H2) = IDTC2H2
      IF ( IDBFGLYC /= 0 ) BFTRACE(IDBFGLYC) = IDTGLYC
      IF ( IDBFHAC  /= 0 ) BFTRACE(IDBFHAC ) = IDTHAC
      IF ( IDBFNAP  /= 0 ) BFTRACE(IDBFNAP ) = IDTNAP  ! (hotp, mpayer, 7/6/11)

      ! Echo biofuel tracer information
      WRITE( 6, 100 ) BFTRACE( 1:NBFTRACE )
 100  FORMAT( 'TRACERID: Biofuel burning tracers        :', 20i3 )
      
      ! Return to calling program
      END SUBROUTINE SET_BFTRACE

!------------------------------------------------------------------------------

      SUBROUTINE INIT_BIOFUEL
!
!******************************************************************************
!  Subroutine INIT_BIOFUEL allocates and zeroes the BIOFUEL array. 
!  (bmy, 9/12/00, 10/15/02)
!
!  NOTES:
!  (1 ) Increase BIOFUEL array from 2 to NBIOFUEL=10 elements (bmy, 3/15/01)
!  (2 ) Make sure NBFTRACE > 0 before allocating BIOFUEL (bmy, 4/17/01)
!  (3 ) BIOFUEL is now declared (NBFTRACE,IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (5 ) Now references ALLOC_ERR from "error_mod.f".  Also references IDBFNOX,
!        IDBFCO, etc from "tracerid_mod.f" (bmy, 11/6/02)
!  (6 ) Replace LWOODCO w/ LBIOFUEL from "logical_mod.f" (bmy, 7/19/04)
!  (7 ) Remove reference to TRACERID_MOD, it's obsolete (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE LOGICAL_MOD,  ONLY : LBIOFUEL

#     include "CMN_SIZE"  ! Size parameters, etc

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_BIOFUEL begins here!
      !=================================================================
      IF ( LBIOFUEL .and. NBFTRACE > 0 ) THEN
         ALLOCATE( BIOFUEL( NBFTRACE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL' )
         BIOFUEL = 0d0

         ! This is a local array to hold biofuel in kg
         ALLOCATE( BIOFUEL_KG( NBFTRACE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL_KG' )
         BIOFUEL_KG = 0d0
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_BIOFUEL

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BIOFUEL
!
!******************************************************************************
!  Subroutine CLEANUP_BIOFUEL deallocates the BIOFUEL array (bmy, 9/11/00)
!******************************************************************************
!      
      ! CLEANUP_BIOFUEL begins here!
      IF ( ALLOCATED( BIOFUEL    ) ) DEALLOCATE( BIOFUEL    )
      IF ( ALLOCATED( BIOFUEL_KG ) ) DEALLOCATE( BIOFUEL_KG )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOFUEL

!------------------------------------------------------------------------------
      
      END MODULE BIOFUEL_MOD
