! $Id: biofuel_mod.f,v 1.6 2005/10/20 14:03:13 bmy Exp $
      MODULE BIOFUEL_MOD
!
!******************************************************************************
!  Module BIOFUEL_MOD contains arrays and routines to compute yearly
!  biofuel emissions for NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, 
!  CH2O, and C2H6 (bmy, 9/12/00, 10/3/05)
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
!  (1 ) BIOFUEL_BURN       : reads data from disk & computes biofuel emissions
!  (2 ) SCALE_BIOFUEL_CO   : scales biofuel CO to account for VOC oxidation
!  (3 ) SCALE_BIOFUEL_ACET : scales biofuel ACET to match a posteriori source
!  (4 ) SET_BFTRACE        : Initializes NBFTRACE counter and BFTRACE array
!  (5 ) INIT_BIOFUEL       : initializes the BIOFUEL array
!  (6 ) CLEANUP_BIOFUEL    : deallocates the BIOFUEL array
!
!  GEOS-CHEM modules referenced by biofuel_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing DAO met field arrays
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) epa_nei_mod.f     : Module containing routines to read EPA/NEI99 data
!  (6 ) error_mod.f     : Module containing NaN and other error check routines
!  (7 ) logical_mod.f   : Module containing GEOS-CHEM logical switches
!  (8 ) tracer_mod.f    : Module containing GEOS-CHEM tracer array etc.
!  (9 ) tracerid_mod.f  : Module containing pointers to tracers & emissions 
!  (10) transfer_mod.f  : Module containing routines to cast & resize arrays
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biofuel_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE BIOFUEL_KG

      ! PRIVATE module routines
      PRIVATE SCALE_BIOFUEL_CO   
      PRIVATE SCALE_BIOFUEL_ACET 
      PRIVATE INIT_BIOFUEL      

      !=================================================================     
      ! MODULE VARIABLES
      !=================================================================
      INTEGER, PARAMETER  :: NBFMAX = 10

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
!  (rvm, acs, bnd, bmy, 9/12/00, 10/3/05)
!
!  Biofuel emissions are based on estimates by Rose Yevich and Jennifer
!  Logan (reference TBA).
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
!******************************************************************************
!
      ! References to F90 modules
      !USE BPCH2_MOD
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD29, AD32_bf, AD34
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE EPA_NEI_MOD,   ONLY : GET_EPA_BIOFUEL, GET_USA_MASK
      USE LOGICAL_MOD,   ONLY : LNEI99
      USE TIME_MOD,      ONLY : GET_DAY_OF_WEEK
      !----------------------------------------------------------------
      ! Prior to 10/3/05:
      !USE TRACER_MOD
      !----------------------------------------------------------------
      USE TRACERID_MOD,  ONLY : IDBFCO,  IDBFNOX, IDTACET, IDTALD2
      USE TRACERID_MOD,  ONLY : IDTALK4, IDTC2H6, IDTC3H8, IDTCH2O 
      USE TRACERID_MOD,  ONLY : IDTCO,   IDTMEK,  IDTNOX,  IDTPRPE 
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      IMPLICIT NONE
      
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND29, ND32, ND34
      
      ! Local variables
      LOGICAL            :: WEEKDAY
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL, SAVE      :: DO_ND29, DO_ND32, DO_ND34
      INTEGER            :: AS, I, J, N, NN, DAY_NUM
      REAL*4             :: ARRAY(IGLOB,JGLOB,1)
      REAL*8,  SAVE      :: MOLWT(NBFMAX)
      REAL*8             :: TOTAL, BXHEIGHT_CM, EPA_NEI
      CHARACTER(LEN=255) :: FILENAME 
      
      ! External functions
      REAL*8, EXTERNAL   :: BOXVL

      !=================================================================
      !   B i o f u e l   B u r n i n g   B e g i n s   H e r e !!
      !
      ! Do the following on the first call to BIOFUEL_BURN... 
      !=================================================================
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
         !
         ! These emissions are converted to [molec/cm3/s] (or 
         ! [molec C/cm3/s] for hydrocarbons), since the chemistry
         ! requires these units.
         !
         ! There are NBFMAX=10 maximum allowed biofuel species, but 
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

               ! Read biofuel NOx emissions in [kg/box/yr] -- tracer #1
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 1 , 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 
                         
               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )
               
               ! Compute total of biofuel NOx
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'NOx ', TOTAL, '[Tg  /yr]'
               
               ! Define MOLWT for use below
               MOLWT(N) = 46d-3

            ELSE IF ( NN == IDTCO ) THEN

               ! Read biofuel CO emissions in [kg/box/yr] -- tracer #4
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 4, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Scale CO to account for oxidation of extra VOC's
               CALL SCALE_BIOFUEL_CO( BIOFUEL_KG(N,:,:) )
               
               ! Print total CO
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'CO  ', TOTAL, '[Tg  /yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 28d-3

            ELSE IF ( NN == IDTALK4 ) THEN

               ! Read biofuel ALK4 emissions in [kg/box/yr] -- tracer #5
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 5, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total ALK4
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ALK4', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTACET ) THEN

               ! Read biofuel ACET emissions in [kg/box/yr] -- tracer #9
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 9, 
     &                          0d0,       IGLOB,     JGLOB,     
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Scale to match a posteriori source (bdf, bmy, 9/10/01)
               CALL SCALE_BIOFUEL_ACET( BIOFUEL_KG(N,:,:) )

               ! Compute total ACET
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ACET', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTMEK ) THEN

               ! Read biofuel MEK emissions in [kg/box/yr] -- tracer #10
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 10, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total MEK
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'MEK ', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTALD2 ) THEN

               ! Read biofuel ALD2 emissions in [kg/box/yr] -- tracer #11
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 11, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total ALD2
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'ALD2', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTPRPE ) THEN

               ! Read biofuel PRPE emissions in [kg/box/yr] -- tracer #18
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 18, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total PRPE
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'PRPE', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTC3H8 ) THEN

               ! Read biofuel C3H8 emissions in [kg/box/yr] -- tracer #19
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 19, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total C3H8
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C3H8', TOTAL, '[Tg C/yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 12d-3

            ELSE IF ( NN == IDTCH2O ) THEN

               ! Read biofuel CH2O emissions in [kg/box/yr] -- tracer #20
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 20, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total CH2O
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'CH2O', TOTAL, '[Tg  /yr]'

               ! Define MOLWT for use below
               MOLWT(N) = 30d-3

            ELSE IF ( NN == IDTC2H6 ) THEN

               ! Read biofuel C2H6 emissions in [kg/box/yr] -- tracer #21
               CALL READ_BPCH2( FILENAME, 'BIOFSRCE', 21, 
     &                          0d0,       IGLOB,     JGLOB,      
     &                          1,         ARRAY,     QUIET=.TRUE. ) 

               ! Cast from REAL*4 to REAL*8, resize to (IIPAR,JJPAR)
               CALL TRANSFER_2D( ARRAY(:,:,1), BIOFUEL_KG(N,:,:) )

               ! Compute total C2H6
               TOTAL = SUM( BIOFUEL_KG(N,:,:) ) * 1d-9
               WRITE( 6, 120 ) 'C2H6', TOTAL, '[Tg C/yr]'
               
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
            
            ! If we are over the USA ...
            IF ( LNEI99 .and. GET_USA_MASK( I, J ) > 0d0 ) THEN
               
               ! Get GEOS-CHEM tracer number
               NN      = BFTRACE(N)

               ! Get EPA/NEI biofuel [molec/cm2/s or atoms C/cm2/s]
               EPA_NEI = GET_EPA_BIOFUEL( I, J, NN, WEEKDAY )

               ! Convert [molec/cm2/s] to [molec/cm3/s]
               BIOFUEL(N,I,J) = EPA_NEI / BXHEIGHT_CM

            ENDIF

            ! ND34 -- archive biofuel burning species [molec/cm2/s]
            IF ( DO_ND34 ) THEN
               AD34(I,J,N) = AD34(I,J,N) + ( BIOFUEL(N,I,J) * 
     &                                       BXHEIGHT_CM )
            ENDIF
         ENDDO  

         ! ND29 -- CO source diagnostics [molec/cm2/s]
         IF ( DO_ND29 ) THEN
            AD29(I,J,3) = AD29(I,J,3) + ( BIOFUEL(IDBFCO,I,J) * 
     &                                    BXHEIGHT_CM )
         ENDIF

         ! ND32 -- NOx source diagnostics [molec/cm2/s]
         IF ( DO_ND32 ) THEN 
            AD32_bf(I,J) = AD32_bf(I,J) + ( BIOFUEL(IDBFNOX,I,J) * 
     &                                      BXHEIGHT_CM )
         ENDIF
      ENDDO  
      ENDDO
!$OMP END PARALLEL DO

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

      SUBROUTINE SET_BFTRACE
!
!******************************************************************************
!  Subroutine SET_NBFTRACE sets the NBFTRACE variable with the number of
!  biofuel tracers that are turned on.  This was split off from "tracerid.f"
!  in order to prevent circular module references. (bmy, 11/6/02, 10/3/05)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDBFACET, IDBFALD2, IDBFALK4, IDBFC2H6 
      USE TRACERID_MOD, ONLY : IDBFC3H8, IDBFCH2O, IDBFCO,   IDBFMEK  
      USE TRACERID_MOD, ONLY : IDBFNOX,  IDBFPRPE, IDTACET,  IDTALD2 
      USE TRACERID_MOD, ONLY : IDTALK4,  IDTC2H6,  IDTC3H8,  IDTCH2O
      USE TRACERID_MOD, ONLY : IDTCO,    IDTMEK,   IDTNOX,   IDTPRPE 

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
      !------------------------------------------
      ! Prior to 10/3/05
      !USE TRACERID_MOD
      !------------------------------------------

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
