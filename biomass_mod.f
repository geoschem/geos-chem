! $Id: biomass_mod.f,v 1.4 2004/09/21 18:04:08 bmy Exp $
      MODULE BIOMASS_MOD
!
!******************************************************************************
!  Module BIOMASS_MOD contains arrays and routines to compute monthly
!  biomass burning emissions for NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, 
!  C3H8, CH2O, C2H6, CH4, and CH3I. (bmy, 9/11/00, 7/20/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NBIOMAX            : maximum # of biomass burning tracers
!  (2 ) BIOTRCE            : index array of biomass burning tracers
!  (3 ) NBIOTRCE           : number of biomass burning tracers
!  (4 ) BIOMASS            : array of biomass burning emissions [molec/cm2/s]
!  (5 ) BURNEMIS           : array of biomass burning emissions [molec/cm3/s]
!  (6 ) TOMSAISCALE        : array for TOMS aerosol index values
!
!  Module Routines:
!  ============================================================================
!  (1 ) BIOBURN            : reads data from disk & computes biomass emissions
!  (2 ) READ_BIOMASS       : reads biomass burning data from binary punch file
!  (3 ) SCALE_BIOMASS_CO   : applies scale factors to CO for VOC oxidation
!  (4 ) SCALE_BIOMASS_ACET : applies scale factors to ACET 
!  (5 ) TOTAL_BIOMASS_TG   : prints monthly biomass emission totals in [Tg (C)]
!  (6 ) ADJUST_TO_TOMSAI   : wrapper for subroutine TOMSAI
!  (7 ) TOMSAI             : adjusts BB for interannual var'bilty w/ TOMS data
!  (8 ) SET_BIOTRCE        : Initializes NBIOTRCE counter and BIOTRCE array
!  (9 ) INIT_BIOMASS       : initializes the BURNEMIS and BIOTRCE arrays 
!  (10) CLEANUP_BIOMASS    : deallocates BURNEMIS, BIOTRCE
!
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f       : Module containing arrays for DAO met fields
!  (3 ) diag_mod.f      : Module containing GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f     : Module containing I/O error and NaN check routines
!  (6 ) grid_mod.f      : Module containing horizontal grid information
!  (7 ) logical_mod.f   : Module containing GEOS-CHEM logical switches
!  (8 ) time_mod.f      : Module containing routines for computing time & date
!  (9 ) tracerid_mod.f  : Module containing pointers to tracers & emissions
!  (10) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  Decision Tree for Biomass Burning Emissions:
!  ============================================================================
!
!  The cases below are described in "Interannual and Seasonal Variability of 
!  Biomass Burning Emissions Constrained by Remote-Sensed Observations" 
!  by Duncan et al.
!
!  Case  LBBSEA LTOMSAI  (LBBSEA and LTOMSAI are flags in "input.geos")
!
!  (1)     T       F     Average monthly BB emissions
!                        ------------------------------------------------------
!                        Read average monthly BB emissions. The mean monthly 
!                        emissions from biomass burning were estimated from 
!                        about four years of ATSR & AVHRR data. See 
!                        Sections 3.1 and 4 of Duncan et al.
!
!
!  (2)     T       T     Interannual varying monthly BB emissions
!                        ------------------------------------------------------
!                        (a) Read annual BB emissions (i.e., inventory of 
!                        Jennifer Logan & Rose Yevich) and impose 
!                        time-dependency by scaling to TOMS AI data for 
!                        those regions where TOMS AI was used. This option 
!                        allows the user to account for the interannual 
!                        variability of BB. 
!
!                        (b) Read average monthly BB emissions for Africa and
!                        areas where TOMS AI is not used.
!
!                        See Sections 3.2 and 5 of Duncan et al.
!
!
!  (3)     F       T     Same as Case 2, except with higher spatial resolution
!                        ------------------------------------------------------
!                        (a) Same as Case 2 prior to 8/1/1996.
!
!                        (b) After 8/1/1996, read monthly BB emissions from 
!                        disk.  The emissions are time-dependent as in Case 2
!                        and account for interannual variation. The spatial 
!                        resolution of emissions is greater than in Case 2 
!                        due to ATSR fire-counts.
!
!                        See Section 3.3 of Duncan et al.
!
!
!  (4)     F       F     Same as Case 3b
!                        ------------------------------------------------------
!                        Read interannual variability BB emissions from disk.
!
!                        See Section 3.3 of Duncan et al.
!
!  NOTES:
!  (1 ) Now treat BURNEMIS as a true global array of size (IGLOB,JGLOB);
!        use offsets IREF = I + I0 and JREF = J + J0 to index it (bmy, 9/12/00)
!  (2 ) Added subroutines READ_BIOMASS and TOMSAI (bmy, 9/25/00)
!  (3 ) Bug fixes in routines BIOBURN and READ_BIOMASS.  Added new decision
!        tree in BIOBURN.  Added routine ADJUST_TO_TOMSAI. (bmy, 10/12/00)
!  (4 ) Updated boundaries of geographic regions in TOMSAI (bnd, bmy, 10/16/00)
!  (5 ) Bug fix for CTM_LAT in TOMSAI (bnd, bmy, 11/28/00)
!  (6 ) Removed obsolete code in BIOBURN (bmy, 12/21/00)
!  (7 ) Now account for extra production of CO from VOC's for Tagged CO
!        and CO-OH simulations (bmy, 1/3/01)
!  (8 ) Now use routines from "error_mod.f" for trapping NaN's (bmy, 3/8/01)
!  (9 ) Moved NBIOMAX here from "CMN_SIZE" (bmy, 3/16/01)
!  (10) Now dimension BIOTRCE and to be of size NBIOMAX, instead of having 
!        them be allocatable.  Also change NBIOMAX from 9 to 10, since we 
!        will be adding ALK4 soon.  Elminate LDOBIOEMIT, since that is now
!        confusing and unnecessary. (bmy, 4/17/01)
!  (11) Bug fix: For option 2 in the decision tree above, scale annual
!        BB emissions to the TOMS aerosol index instead of seasonal.  This
!        will give the correct results.  Updated routines ADJUST_TO_TOMSAI
!        and TOMSAI accordingly. (bnd, bmy, 6/6/01)
!  (12) PRPE is already in molec C, so don't multiply it by 3 as we have
!        been doing before. (bmy, 6/29/01)
!  (13) Update comments for BB decision tree (bnd, bmy, 7/2/01)
!  (14) Now use correct scale factors for CO (bnd, bmy, 8/21/01)
!  (15) Bug fix: Make sure to read data from the biomass burning punch file
!        with the correct index for runs that have less than NBIOMAX species
!        turned on. (bmy, 8/24/01)
!  (16) Add new routine: SCALE_BIOMASS_ACET.  Also updated comments. 
!        (bmy, 9/6/01)
!  (17) Removed obsolete code (bmy, 9/18/01)
!  (18) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (19) Removed duplicate variable definitions.  Also now can specify 
!        biomass burning subdirectory via a variable in BIOBURN (bmy, 11/15/01)
!  (20) Now point to new biomass burning files from 10/2001 (bmy, 12/4/01)
!  (21) Updated comments (1/15/02)
!  (22) Fixed incorrect value for IPICK in "adjust_to_tomsai" (bmy, 2/27/02)
!  (23) Bug fix: convert from [molec/cm2/s] to [molec/cm3/s] every timestep.
!        Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments.  Renamed INIT_BURNEMIS
!        to INIT_BIOMASS.  BIOMASS is now an allocatable module array
!        instead of a SAVEd array within routine BIOBURN. (bmy, 5/30/02)
!  (24) Now reference BXHEIGHT from "dao_mod.f".  Now references "error_mod.f".
!        Also deleted obsolete code from various routines.  Now references
!        "tracerid_mod.f". (bmy, 11/6/02)
!  (25) Now references "grid_mod.f" and the new "time_mod.f".  Also suppresses
!        printing when calling routine READ_BPCH2.  Bug fix in routine TOMSAI.
!        Fixed bug in BIOBURN when passing arrays BIOMASS_SEA and BIOMASS_ANN
!        to routine READ_BIOMASS. (bmy, 4/28/03)
!  (26) Now references "directory_mod.f" & "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biomass_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE NAIREGIONS
      PRIVATE NAIYEARS
      PRIVATE NMONTHSAI
      PRIVATE TOMSAISCALE
      PRIVATE BIOMASS

      ! PRIVATE module routines
      PRIVATE INIT_BIOMASS
      PRIVATE READ_BIOMASS
      PRIVATE SCALE_BIOMASS_CO
      PRIVATE SCALE_BIOMASS_ACET
      PRIVATE TOTAL_BIOMASS_TG
      PRIVATE ADJUST_TO_TOMSAI

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER, PARAMETER   :: NBIOMAX = 10

      INTEGER              :: NBIOTRCE
      INTEGER              :: BIOTRCE(NBIOMAX)
      
      REAL*8,  ALLOCATABLE :: BURNEMIS(:,:,:)
      REAL*8,  ALLOCATABLE :: BIOMASS(:,:,:)

      ! TOMS AI interannual variability in biomass burning emissions
      INTEGER, PARAMETER   :: NAIREGIONS = 8
      INTEGER, PARAMETER   :: NAIYEARS   = 21
      INTEGER, PARAMETER   :: NMONTHSAI  = NAIYEARS * 12 

      REAL*8,  ALLOCATABLE :: TOMSAISCALE(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE BIOBURN
!
!******************************************************************************
!  Subroutine BIOBURN computes the biomass burning emissions for several
!  species for the given month (jal, acs, rvm, bmy, 9/11/00, 5/16/03)
!
!  NOTES:
!  (1 ) Incorporated original functionality of "bioburn.f" and "biomass.h"
!        in F90 module "bioburn_mod.f".  Biomass burning arrays now are only
!        allocated if biomass burning is turned on.  (bmy, 9/11/00)
!  (2 ) Split off calls to READ_BPCH2 into separate subroutine READ_BIOMASS 
!        for clarity.  Also now use logical switches LBBSEA and LTOMSAI to
!        switch between seasonal or interannual variability. (bmy, 9/28/00)
!  (3 ) Bug fixes: (a) Acetone is BIOMASS(5,:,:), not BIOMASS(9,:,:). 
!        (b) Make sure to read in all biomass burning tracers from the 
!        binary punch file, regardless of which tracers are actually emitted. 
!        (bmy, 10/11/00)
!  (4 ) Added new decision tree (see comments above) (bmy, 10/12/00)
!  (5 ) Removed obsolete code from 10/12/00 (bmy, 12/21/00)
!  (6 ) Enhance CO from biomass burning by 10% for Tagged CO and CO-OH
!        simulations, to account for extra production of CO from VOC's.
!        (bnd, bmy, 1/3/01)
!  (7 ) Now use interface IT_IS_NAN (from "error_mod.f") to trap NaN's.
!        This will work on DEC/Compaq and SGI platforms. (bmy, 3/8/01)
!  (8 ) Now call INIT_BURNEMIS on the very first call to BIOBURN.  Also
!        read biomass burning species w/o using LDOBIOEMIT, which is 
!        now unnecessary.  Call SCALE_BIOMASS_CO to multiply CO biomass
!        burning emissions by jal/bnd scale factors, to account for
!        oxidation of VOC's not carried (bmy, 4/17/01)
!  (9 ) Now read new biomass burning files (Apr 2001) from the 
!        "biomass_200104/" subdirectory of DATA_DIR.  (bmy, 4/18/01)
!  (10) Added BIOMASS_SEA and BIOMASS_ANN arrays for the scaling for Case #2
!        in the decision tree above.  This will scale the annual BB emissions
!        using TOMSAI in selected regions, but use the seasonal emissions
!        elsewhere. (bnd, bmy, 6/6/01)
!  (11) Now call SCALE_BIOMASS_ACET in order to enhance biomass burning ACET 
!        by 77%, to match results from Jacob et al 2001. (bdf, bmy, 9/4/01)
!  (12) BURNEMIS, BIOMASS, BIOMASS_SEA, and BIOMASS_ANN are now dimensioned
!        (NBIOTRCE,IIPAR,JJPAR).  BURNEMIS(:,IREF,JREF) is now 
!        BURNEMIS(:,I,J) and BIOMASS(:,IREF,JREF) is now BIOMASS(:,I,J).
!        Remove IREF, JREF, IOFF, JOFF -- these are obsolete. (bmy, 9/28/01)
!  (13) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (14) Removed duplicate definition of BOXVL.  Also added BIOMASS_DIR
!        string to specify the sub-directory of DATA_DIR where biomass
!        emissions are kept. (bmy, 11/15/01) 
!  (15) Now set BIOMASS_MOD = 'biomass_200110/' as the default.  This points
!        to newer biomass burning emissions from Randall Martin (bmy, 11/30/01
!  (16) Now set BIOMASS_DIR = 'biomass_200010/' in order to take advantage of 
!        new biomass burning files from Randall Martin (w/ firecounts thru 
!        2000).  
!  (17) Bug fix: convert from [molec/cm2/s] to [molec/cm3/s] every timestep.
!        Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections. Now call INIT_BIOMASS instead of 
!        INIT_BURNEMIS.  Added parallel DO-loops for unit conversion.  Now
!        archive diagnostics w/ in the parallel loop section. (bmy, 5/31/02)
!  (18) Now reference BXHEIGHT from "dao_mod.f".  Also call GEOS_CHEM_STOP
!        to free memory when stopping with an error.  Now call GET_TAU0 with
!        3 arguments instead of 2.  Now references IDTNOX, IDBNOX, etc. from
!        "tracerid_mod.f". (bmy, 11/6/02)
!  (19) Now remove IMONTH from the arg list.  Now use functions GET_MONTH,
!         GET_TAU, GET_YEAR, and ITS_A_LEAPYEAR from "time_mod.f".  
!         (bmy, 2/10/03)
!  (20) Bug fix: make sure only to pass BIOMASS_SEA(1:NBIOTRCE,:,:) and 
!         BIOMASS_ANN(1:NBIOTRCE,:,:) to READ_BIOMASS. (bnd, bmy, 5/16/03)
!  (21) Added fancy output (bmy, 4/26/04)
!  (22) Removed reference to CMN, it's obsolete.  Now reference DATA_DIR from
!        "directory_mod.f".  Now references LBBSEA and LTOMSAI from
!        "logical_mod.f". (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DAO_MOD,       ONLY : BXHEIGHT
      USE DIAG_MOD,      ONLY : AD28, AD29, AD32_bb
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : IT_IS_NAN,      GEOS_CHEM_STOP
      USE LOGICAL_MOD,   ONLY : LBBSEA,         LTOMSAI
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR, GET_MONTH, 
     &                          GET_TAU,        GET_YEAR       
      USE TRACERID_MOD

#     include "CMN_SIZE"   ! Size parameters
!------------------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN"        ! NSRCX
!------------------------------------------------------------
#     include "CMN_DIAG"   ! Diagnostic arrays & switches
!------------------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP"  ! DATA_DIR, LBBSEA, LTOMSAI
!------------------------------------------------------------

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL, SAVE       :: DO_ND28, DO_ND29, DO_ND32

      INTEGER             :: I, J, N, L
      INTEGER, SAVE       :: MONTHSAVE = -99

      REAL*4              :: TEMP(IGLOB,JGLOB,1)
      
      REAL*8              :: TIME, XTAU, BXHEIGHT_CM
      REAL*8              :: BIOMASS_SEA(NBIOMAX,IIPAR,JJPAR) 
      REAL*8              :: BIOMASS_ANN(NBIOMAX,IIPAR,JJPAR)

      CHARACTER(LEN=4  )  :: CYEAR
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: BIOMASS_DIR

      ! MONTHDATES = number of days per month
      INTEGER             :: MONTHDATES(12) = (/ 31, 28, 31, 30, 
     &                                           31, 30, 31, 31, 
     &                                           30, 31, 30, 31 /)

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL
      
      !=================================================================
      !    B i o m a s s   B u r n i n g   B e g i n s   H e r e !!
      ! 
      ! GEOS-CHEM has the following biomass burning species:
      !
      !    Species   Index   CTM Tracer #   Units as read from file
      !    ---------------------------------------------------------
      !      NOX       1          1         [molec NOx /cm2/month]
      !      CO        2          4         [molec CO  /cm2/month]
      !      ALK4      3          5         [molec C   /cm2/month]
      !      ACET      4          9         [molec ACET/cm2/month]
      !      MEK       5          10        [molec C   /cm2/month]
      !      ALD2      6          11        [molec C   /cm2/month]
      !      PRPE      7          18        [molec C   /cm2/month]
      !      C3H8      8          19        [molec C3H8/cm2/month]
      !      CH2O      9          20        [molec CH2O/cm2/month]
      !      C2H6      10         21        [molec C2H6/cm2/month]
      !
      ! Subsequent unit conversion is done on the following species:
      !      [molec ACET/cm2/month]  -->  [molec C/cm2/month]
      !      [molec C3H8/cm2/month]  -->  [molec C/cm2/month]
      !      [molec C2H6/cm2/month]  -->  [molec C/cm2/month]
      ! 
      ! There are NBIOMAX=10 maximum allowed biomass burning species, 
      ! but only NBIOTRCE of these are actually emitted.  Species are 
      ! turned off/on with the switches in the "tracer.dat" file.
      !
      ! There are two arrays in routine BIOBURN of that contain
      ! biomass burning emissions data:
      !   (1) BIOMASS  -- [molec/cm2/month] ( or [molec C/cm2/month] )
      !   (2) BURNEMIS -- [molec/cm3/s    ] ( or [molec C/cm3/s    ] )
      !
      ! Biomass burning emissions are first read from disk into the
      ! BIOMASS array.  After unit conversion to [molec/cm3/s] ( or
      ! [molec C/cm3/s]), the emissions are stored in BURNEMIS.
      !
      ! Biomass burning data is monthly, so we only have to read 
      ! emissions from disk once each month.
      !=================================================================
      IF ( FIRST ) THEN

         ! Initialize BURNEMIS and BIOMASS arrays on first call only
         CALL INIT_BIOMASS

         ! Set flags for diagnostics
         DO_ND28 = ( ND28 > 0 .and. NBIOTRCE > 0 ) 
         DO_ND29 = ( ND29 > 0 .and. IDBCO    > 0 )
         DO_ND32 = ( ND32 > 0 .and. IDBNOX   > 0 )

         ! Reset first-time flag
         FIRST   = .FALSE.
      ENDIF

      !=================================================================
      ! Do the following on the first day of a new month...
      !=================================================================
      IF ( GET_MONTH() /= MONTHSAVE ) THEN   
         
         ! Save the current month
         MONTHSAVE = GET_MONTH()

         ! Set MONTHDATES(2) = 29 for leapyears, = 28 otherwise (bmy, 4/19/99)
         IF ( GET_MONTH() == 2 ) THEN
            IF( ITS_A_LEAPYEAR() ) THEN
               MONTHDATES(2) = 29
            ELSE
               MONTHDATES(2) = 28
            ENDIF
         ENDIF

         ! TIME = conversion from [molec/cm2/month] to [molec/cm2/s]
         TIME = ( DBLE( MONTHDATES( GET_MONTH() ) ) * 86400d0 )

         ! Create a string for the 4-digit year
         WRITE( CYEAR, '(i4)' ) GET_YEAR()

         ! Fancy output...
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)' ) 
     &        'B I O M A S S   B U R N I N G   E M I S S I O N S'

         !==============================================================
         ! Set BIOMASS_DIR to the subdirectory where the current
         ! biomass burning files are stored 
         !==============================================================
         !BIOMASS_DIR = 'biomass_200104/'
         BIOMASS_DIR = 'biomass_200110/'

         !==============================================================
         ! Case 1: LBBSEA = T and LTOMSAI = F
         !
         ! Read seasonal biomass burning emissions from disk.
         !==============================================================
         IF ( LBBSEA .and. ( .not. LTOMSAI ) ) THEN

            ! Get TAU0 value to index the punch file -- use generic year 1985
            XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

            ! Filename for seasonal biomass burning emissions
            FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) // 
     &                 'bioburn.seasonal.geos.' // GET_RES_EXT()

            ! Read the seasonal biomass burning emissions from disk
            CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS )

         !==============================================================
         ! Case 2: LBBSEA = T and LTOMSAI = T
         !
         !
         ! Read annual biomass burning emissions from disk, but use
         ! TOMS aerosol index data to impose interannual variability.
         ! Read in seasonal biomass buring emissions for Africa and
         ! regions outside the regions adjusted by TOMS AI.
         !==============================================================
         ELSE IF ( LBBSEA .and. LTOMSAI ) THEN

            ! Get TAU0 value to index the punch file -- use generic year 1985
            XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )            

            ! Filename for seasonal biomass burning emissions
            FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //
     &                 'bioburn.seasonal.geos.' // GET_RES_EXT()

            ! Read the seasonal biomass burning emissions into BIOMASS_SEA
            CALL READ_BIOMASS( FILENAME, XTAU, 
     &                         BIOMASS_SEA(1:NBIOTRCE,:,:) )

            ! Get TAU0 value to index the punch file -- use generic year 1985
            XTAU = GET_TAU0( 1, 1, 1985 )

            ! Filename for annual biomass burning emissions 
            FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //
     &                 'bioburn.annual.geos.' // GET_RES_EXT()

            ! Read the annual biomass burning emissions into BIOMASS_ANN
            CALL READ_BIOMASS( FILENAME, XTAU, 
     &                         BIOMASS_ANN(1:NBIOTRCE,:,:) )

            ! Adjust the annual biomass burning to the TOMS Aerosol 
            ! Index data where necessary.  Otherwise, overwrite
            ! with seasonal data. Save result in the BIOMASS array.
            BIOMASS = 0d0
            CALL ADJUST_TO_TOMSAI( BIOMASS_ANN, BIOMASS_SEA, BIOMASS )

         !==============================================================
         ! Case 3: LBBSEA = F and LTOMSAI = T
         !
         ! (1) Prior to 8/1/1996, read seasonal biomass burning 
         !     emissions, and use TOMS AI data to impose int. var.
         !
         ! (2) On or after 8/1/1996, read the interannual variability
         !     biomass burning emissions (computed by Randall Martin:
         !     rvm@io.harvard.edu) directly from disk.
         !==============================================================
         ELSE IF ( ( .not. LBBSEA ) .and. LTOMSAI ) THEN

            ! 8/1/1996 is TAU value 101520
            IF ( GET_TAU() < 101520d0 ) THEN

               ! Get TAU0 value to index the punch file -- 
               ! use generic year 1985
               XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

               ! Filename for seasonal biomass burning emissions
               FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //
     &                    'bioburn.seasonal.geos.' // GET_RES_EXT()

               ! Read the seasonal biomass burning emissions into BIOMASS_SEA
               CALL READ_BIOMASS( FILENAME, XTAU, 
     &                            BIOMASS_SEA(1:NBIOTRCE,:,:) )

               ! Get TAU0 value to index the punch file -- 
               ! use generic year 1985
               XTAU = GET_TAU0( 1, 1, 1985 )

               ! Filename for annual biomass burning emissions
               FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //
     &                    'bioburn.annual.geos.' // GET_RES_EXT()

               ! Read the annual biomass burning emissions from disk
               CALL READ_BIOMASS( FILENAME, XTAU, 
     &                            BIOMASS_ANN(1:NBIOTRCE,:,:) )

               ! Adjust the annual biomass burning to the TOMS Aerosol 
               ! Index data where necessary.  Otherwise, overwrite
               ! with seasonal data. Save result in the BIOMASS array.
               BIOMASS = 0d0
               CALL ADJUST_TO_TOMSAI( BIOMASS_ANN, BIOMASS_SEA, BIOMASS)

            ELSE

               ! Use actual TAU0 value to index punch file
               XTAU = GET_TAU()

               ! Filename for interannual variability biomass burning emissions
               FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //
     &                    'bioburn.interannual.geos.' // 
     &                    GET_RES_EXT() // '.' // CYEAR           

               ! Read interannual variability biomass burning
               CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS )
            
         ENDIF

         !==============================================================
         ! Case 4: LBBSEA = F and LTOMSAI = F
         !
         ! Read the interannual variability biomass burning emissions 
         ! (computed by Randall Martin: rvm@io.harvard.edu) from disk.
         !==============================================================
         ELSE IF ( ( .not. LBBSEA ) .and. ( .not. LTOMSAI ) ) THEN

            ! Use actual TAU0 value to index punch file
            XTAU = GET_TAU()

            ! Filename for interannual variability biomass burning emissions
            FILENAME = TRIM( DATA_DIR ) // TRIM( BIOMASS_DIR ) //      
     &                'bioburn.interannual.geos.' // 
     &                 GET_RES_EXT() // '.' // CYEAR
               
            ! Read interannual variability biomass burning
            CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS )

         ENDIF

         !==============================================================
         ! Convert BIOMASS from:
         !
         !  [  molec/cm2/month] to [  molec/cm2/s] -> for NOx, CO, CH2O
         !  [atoms C/cm2/month] to [atoms C/cm2/s] -> for hydrocarbons
         !==============================================================
         BIOMASS = BIOMASS / TIME

         ! Fancy output...
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF

      !=================================================================
      ! Do the following on each emission timestep:
      !
      ! For NOx, CO, CH2O, compute:
      !    BURNEMIS [molec/cm3/s] from BIOMASS [molec/cm2/s].  
      !
      ! For ALK4, ACET, MEK, ALD2, PRPE, C3H8, C2HO, C2H6, compute:
      !    BURNEMIS [molec C/cm3/s] from BIOMASS [molec C/cm2/s].
      !
      ! NOTE: We must convert from [molec/cm2/s] to [molec/cm3/s] on
      ! every timestep, since the box heights change with the surface 
      ! pressure over the course of a year. (bmy, 5/28/02)
      !=================================================================

      ! Loop over latitudes and longitudes
      IF ( NBIOTRCE > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, BXHEIGHT_CM, N )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! BXHEIGHT_CM is boxheight in cm
            BXHEIGHT_CM = BXHEIGHT(I,J,1) * 100d0

            ! Loop over the emitted biomass burning tracers
            DO N = 1, NBIOTRCE

               ! Convert from [molec/cm2/s] to [molec/cm3/s]
               BURNEMIS(N,I,J) = BIOMASS(N,I,J) / BXHEIGHT_CM 
                  
               ! Make sure BURNEMIS is not negative
               IF ( BURNEMIS(N,I,J) < 0d0 ) THEN
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  PRINT*, 'BURNEMIS < 0'
                  PRINT*, 'I, J, N     : ', I, J, N
                  PRINT*, 'BIOMASS     : ', BIOMASS(N,I,J)
                  PRINT*, 'BXHEIGHT_CM : ', BXHEIGHT_CM
                  PRINT*, 'STOP in bioburn (biomass_mod.f)'
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  CALL GEOS_CHEM_STOP
               ENDIF

               ! Make sure BURNEMIS is not NaN
               IF ( IT_IS_NAN( BURNEMIS(N,I,J) ) ) THEN
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  PRINT*, 'BURNEMIS is NAN!'
                  PRINT*, 'I, J, N     : ', I, J, N
                  PRINT*, 'BIOMASS     : ', BIOMASS(N,I,J)
                  PRINT*, 'BXHEIGHT_CM : ', BXHEIGHT_CM
                  PRINT*, 'STOP in bioburn (biomass_mod.f)'
                  WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                  CALL GEOS_CHEM_STOP
               ENDIF
            ENDDO

            !===========================================================
            ! Archive the following GEOS-CHEM diagnostics:
            !                                                  
            ! ND28 = diagnostic for all biomass burning species
            ! ND29 = CO-source diagnostic                      
            ! ND32 = NOx-source diagnostic                     
            !===========================================================

            ! ND28 -- store biomass burning emissions in [molec/cm2/s]
            IF ( DO_ND28 ) THEN 
               DO N = 1, NBIOTRCE
                  AD28(I,J,N) = AD28(I,J,N) + ( BURNEMIS(N,I,J) * 
     &                                          BXHEIGHT_CM )
               ENDDO
            ENDIF

            ! ND29 -- store biomass burning of CO in [molec/cm2/s]
            IF ( DO_ND29 ) THEN 
               AD29(I,J,2) = AD29(I,J,2) + ( BURNEMIS(IDBCO,I,J) * 
     &                                       BXHEIGHT_CM )
            ENDIF

            ! ND32 -- store biomass burning of NOx in [molec/cm2/s]
            IF ( DO_ND32 ) THEN
               AD32_bb(I,J) = AD32_bb(I,J) + ( BURNEMIS(IDBNOX,I,J) *
     &                                         BXHEIGHT_CM )
            ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE BIOBURN

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOMASS( FILENAME, TAU0, BIOMASS )
!
!******************************************************************************
!  Subroutine READ_BIOMASS reads the biomass burning emissions from disk
!  in units of [molec (C)/cm2/month].  (bmy, 9/25/00, 3/14/03)
!      
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of the biomass burning file to read
!  (2 ) TAU0     (REAL*8   ) : TAU0 value used to index the BB data 
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) BIOMASS  (REAL*8   ) : Biomass burning emissions (for NBIOMAX tracers)
!
!  NOTES:
!  (1 ) Split off from "bioburn.f" to reduce code duplication (bmy, 9/25/00)
!  (2 ) Now read in all biomass burning tracers from the punch file, 
!        regardless of whether or not they are actually emitted. 
!        (bmy, 10/11/00)
!  (3 ) Now only read in the NBIOTRCE biomass burning tracers that
!        are actually emitted (bmy, 4/17/01)
!  (4 ) PRPE is already in molec C, so don't multiply it by 3 as we have
!        been doing before. (bmy, 6/29/01)
!  (5 ) Bug fix: make sure that tracers get read from the biomass burning
!        file w/ the right index number.  This was a bug for runs that had
!        less than NBIOMAX species specified.  (bmy, 8/24/01)
!  (6 ) Removed obsolete code from 8/24/01 (bmy, 9/18/01)
!  (7 ) BIOMASS is now of size (NBIOMAX,IIPAR,JJPAR).  Now call TRANSFER_2D
!        to copy data from REAL*4 to REAL*8 and also to resize from
!        (IGLOB,JGLOB) to (IIPAR,JJPAR).  (bmy, 9/28/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) BIOMASS needs to be of size (NBIOTRCE,IIPAR,JJPAR) (bmy, 5/31/02)
!  (10) Now references IDTNOX, etc. from "tracerid_mod.f" (bmy, 11/6/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. flag to suppress extra info 
!        from being printed (bmy, 3/14/03)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TRACERID_MOD
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

#     include "CMN_SIZE"

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)  :: FILENAME
      REAL*8,           INTENT(IN)  :: TAU0
      REAL*8,           INTENT(OUT) :: BIOMASS(NBIOTRCE,IIPAR,JJPAR)

      ! Local variables
      INTEGER                       :: N, NN, NT
      REAL*4                        :: ARRAY(IGLOB,JGLOB,1)

      !=================================================================
      ! READ_BIOMASS begins here!
      !
      ! NOTE: BIOMASS is dimensioned for NBIOTRCE species (bmy, 5/31/02)
      !=================================================================
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'BIOBURN:  Reading ', a )

      ! Initialize the BIOMASS array
      BIOMASS = 0d0

      ! Loop over only the emitted biomass tracers
      DO N = 1, NBIOTRCE
           
         ! NN is the actual CTM tracer number of species N
         NN = BIOTRCE(N)
         
         ! Do scaling if necessary and print totals in Tg
         IF ( NN == IDTNOX ) THEN

            ! NOx is stored in the biomass file as tracer #1
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 1, 
     &                       TAU0,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. )  

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! NOX -- print totals in [Tg/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 46d-3, 'NOx' )

         ELSE IF ( NN == IDTCO ) THEN
 
            ! CO is stored in the biomass file as tracer #4
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 4, 
     &                       TAU0,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! CO -- scale to account for oxidation of extra VOC's
            ! also print totals in [Tg/month]
            CALL SCALE_BIOMASS_CO( BIOMASS(N,:,:)              )
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 28d-3, 'CO' )

         ELSE IF ( NN == IDTALK4 ) THEN

            ! ALK4 is stored in the biomass file as tracer #5
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 5, 
     &                       TAU0,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! ALK4 -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'ALK4' )

         ELSE IF ( NN == IDTACET ) THEN

            ! ACET is stored in the biomass file as tracer #9
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 9, 
     &                       TAU0,      IGLOB,     JGLOB,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! ACET -- Convert from [molec/cm2/month] to [molec C/cm2/month] 
            BIOMASS(N,:,:) = BIOMASS(N,:,:) * 3d0  

            ! Scale to yearly value for biogenic acetone (bdf, bmy, 7/23/01)
            CALL SCALE_BIOMASS_ACET( BIOMASS(N,:,:) )

            ! Print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'ACET' )

         ELSE IF ( NN == IDTMEK ) THEN

            ! MEK is stored in the biomass file as tracer #10
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 10, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! MEK -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'MEK' )

         ELSE IF ( NN == IDTALD2 ) THEN

            ! ALD2 is stored in the biomass file as tracer #11
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 11, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! ALD2 -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'ALD2' )

         ELSE IF ( NN == IDTPRPE ) THEN

            ! PRPE is stored in the biomass file as tracer #18
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 18, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! PRPE -- convert from [molec/cm2/month] to [molec C/cm2/month]
            ! Print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'PRPE' )
            
         ELSE IF ( NN == IDTC3H8 ) THEN
               
            ! C3H8 is stored in the biomass file as tracer #19
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 19, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! C3H8 -- convert from [molec/cm2/month] to [molec C/cm2/month] 
            ! Print totals in [Tg C]
            BIOMASS(N,:,:) = BIOMASS(N,:,:) * 3d0 
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'C3H8' )

         ELSE IF ( NN == IDTCH2O ) THEN

            ! CH2O is stored in the biomass file as tracer #20
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 20, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! CH2O -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 30d-3, 'CH2O' )

         ELSE IF ( NN == IDTC2H6 ) THEN

            ! C2H6 is stored in the biomass file as tracer #21
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 21, 
     &                       TAU0,      IGLOB,     JGLOB,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(N,:,:) )

            ! C2H6 --convert from [molec/cm2/month] to [molec C/cm2/month]
            ! Print totals in [Tg C]
            BIOMASS(N,:,:) = BIOMASS(N,:,:) * 2d0 
            CALL TOTAL_BIOMASS_TG( BIOMASS(N,:,:), 12d-3, 'C2H6' )

         ENDIF
      ENDDO
     
      ! Return to calling program
      END SUBROUTINE READ_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE SCALE_BIOMASS_CO( BBARRAY )
!
!******************************************************************************
!  Subroutine SCALE_BIOMASS_CO multiplies the CO biomass emissions by scale 
!  factors to account for CO production from VOC's that are not explicitly 
!  carried in the chemistry mechanisms. (bnd, bmy, 8/21/01, 7/20/04)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Array containing biomass burning CO emissions
!
!  NOTES:
!  (1 ) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!       Bryan Duncan (bnd@io.harvard.edu) and Daniel Jacob (djj@io.harvard.edu)
!  (2 ) Scale factors have been corrected to 5% and 11% (bnd, bmy, 8/21/01)
!  (3 ) BBARRAY is now dimensioned (IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (5 ) Now references ITS_A_FULLCHEM_SIM, ITS_A_TAGCO_SIM from "tracer_mod.f"
!        (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGCO_SIM

#     include "CMN_SIZE"    ! Size parameters
!------------------------------------------------
!#     include "CMN"         ! NSRCX 
!------------------------------------------------
      
      ! Arguments
      REAL*8, INTENT(INOUT) :: BBARRAY(IIPAR,JJPAR) 

      !=================================================================
      ! SCALE_BIOMASS_CO begins here!
      !=================================================================
!------------------------------------------------------------------------------
! Prior to 7/20/04:
!      SELECT CASE ( NSRCX )
!
!         ! Full chemistry w/ SMVGEAR  -- enhance by 5%
!         CASE ( 3 ) 
!            BBARRAY = BBARRAY * 1.05d0
!         
!         ! CO-OH parameterization -- implement later
!         !CASE ( 5 )
!         !   BFARRAY = BFARRAY * 
!
!         ! Tagged CO -- enhance by 11%
!         CASE ( 7 )
!            BBARRAY = BBARRAY * 1.11d0
!         
!         CASE DEFAULT
!            ! Nothing
!
!      END SELECT 
!-----------------------------------------------------------------------------

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Full chemistry w/ SMVGEAR  -- enhance by 5%
         BBARRAY = BBARRAY * 1.05d0
         
      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN

         ! Tagged CO -- enhance by 11%
         BBARRAY = BBARRAY * 1.11d0

      !ELSE IF ( ITS_A_COPARAM_SIM()
      !   BBARRAY = BBARRAY * 

      ENDIF

      ! Return to calling program  
      END SUBROUTINE SCALE_BIOMASS_CO

!------------------------------------------------------------------------------

      SUBROUTINE SCALE_BIOMASS_ACET( BBARRAY )
!
!******************************************************************************
!  Subroutine SCALE_BIOMASS_ACET scales the seasonal acetone biomass
!  burning emissions (Case 1 in the decision tree above) to a given 
!  yearly value.  This is needed for the new biogenic emission fluxes.
!  (bdf, bmy, 9/4/01, 7/20/04)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Array containing biomass burning CO emissions
!
!  Reference:
!  ============================================================================
!  Jacob, D.J., B.D. Field, E. Jin, I. Bey, Q. Li, J.A. Logan, and 
!    R.M. Yantosca, Atmospheric budget of acetone, submitted to 
!    Geophys. Res. Lett., 2001. 
!
!  NOTES:
!  (1 ) Scale factors determined by Brendan Field, in order to match that
!        of the acetone paper: Jacob et al, 2001. (bdf, bmy, 9/4/01)
!  (2 ) BBARRAY is now dimensioned (IIPAR,JJPAR) (bmy, 9/28/01)
!  (3 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (4 ) Now reference LBBSEA, LTOMSAI, from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD, ONLY : LBBSEA, LTOMSAI

#     include "CMN_SIZE"    ! Size parameters
!----------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP"   ! LBBSEA, LTOMSAI
!----------------------------------------------------

      ! Arguments
      REAL*8, INTENT(INOUT) :: BBARRAY(IIPAR,JJPAR)

      !=================================================================
      ! SCALE_BIOMASS_ACET begins here!
      !
      ! Apply scale factor from Jacob et al 2001 (bdf)
      !=================================================================
      IF ( LBBSEA .and. .not. LTOMSAI ) THEN
         BBARRAY = BBARRAY * 1.77d0  
      ENDIF

      ! Return to calling program
      END SUBROUTINE SCALE_BIOMASS_ACET

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_BIOMASS_TG( BBARRAY, MOLWT, NAME )
!
!******************************************************************************
!  Subroutine TOTAL_BIOMASS_TG prints the amount of biomass burning emissions 
!  that are emitted each month in Tg or Tg C. (bmy, 3/20/01, 3/14/03)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Biomass burning CO emissions [molec/cm2/month]
!
!  NOTES:
!  (1 ) BBARRAY is now dimensioned (IIPAR,JJPAR).  Also, DXYP is dimensioned
!        as JGLOB, so use J+J0 to reference it. (bmy, 9/28/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (3 ) Now use function GET_AREA_CM2 from "grid_mod.f" to compute grid
!        box surface area in cm2.  Removed reference to CMN header file.
!        Cosmetic changes. (bmy, 3/14/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_CM2

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      REAL*8,           INTENT(IN) :: BBARRAY(IIPAR,JJPAR) 
      REAL*8,           INTENT(IN) :: MOLWT
      CHARACTER(LEN=*), INTENT(IN) :: NAME

      ! Local variables
      INTEGER                      :: I, J
      REAL*8                       :: TOTAL, A_CM2
      CHARACTER(LEN=6)             :: UNIT

      !=================================================================
      ! TOTAL_BIOMASS_TG begins here!
      !=================================================================

      ! Initialize summing variable
      TOTAL = 0d0

      ! Convert from [molec  /cm2/month] to [kg  /month]
      ! or      from [molec C/cm2/month] to [kg C/month]
      DO J = 1, JJPAR
         
         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR
            TOTAL = TOTAL + BBARRAY(I,J) * A_CM2 * ( MOLWT / 6.023d23 )
         ENDDO
      ENDDO
     
      ! Convert from kg --> Tg
      TOTAL = TOTAL * 1d-9

      ! Define unit string
      IF ( NAME == 'NOx' .or. NAME == 'CO' .or. NAME == 'CH2O' ) THEN
         UNIT = '[Tg  ]'
      ELSE
         UNIT = '[Tg C]'
      ENDIF

      ! Write totals
      WRITE( 6, 100 ) NAME, TOTAL, UNIT
 100  FORMAT( 'Sum Biomass ', a4, 1x, ': ', f9.3, 1x, a9  )
      ! Return to calling program
      END SUBROUTINE TOTAL_BIOMASS_TG

!------------------------------------------------------------------------------

      SUBROUTINE ADJUST_TO_TOMSAI( BIOMASS_ANN, BIOMASS_SEA, BIOMASS )
!
!******************************************************************************
!  Subroutine ADJUST_TO_TOMSAI is a wrapper for subroutine TOMSAI.
!  (bmy, 10/12/00, 2/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) BIOMASS_ANN (REAL*8 ) : Annual   biomass emissions  [molec/cm2/month]
!  (2 ) BIOMASS_SEA (REAL*8 ) : Seasonal biomass emissions  [molec/cm2/month]
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) BIOMASS     (REAL*8 ) : Adjusted biomass emisssions [molec/cm2/month]
!
!  NOTES:
!  (1 ) Bug fix: Now scale annual BB emissions to TOMS for selected
!        regions, and overwrite w/ seasonal BB emissions elsewhere. 
!        (bnd, bmy, 6/6/01)
!  (2 ) BIOMASS_ANN, BIOMASS_SEA, and BIOMASS are now all of size
!        (NBIOMAX,IIPAR,JJPAR). (bmy, 9/28/01)
!  (3 ) Removed obsolete code from 9/01 (bmy, 10/23/01) 
!  (4 ) Remove IMONTH from arg list.  Remove IMONTH from call to TOMSAI 
!        (bmy, 2/11/03)
!******************************************************************************
!
#     include "CMN_SIZE"

      ! Arguments 
      REAL*8,  INTENT(INOUT) :: BIOMASS_ANN(NBIOMAX,IIPAR,JJPAR)
      REAL*8,  INTENT(INOUT) :: BIOMASS_SEA(NBIOMAX,IIPAR,JJPAR)
      REAL*8,  INTENT(INOUT) :: BIOMASS(NBIOMAX,IIPAR,JJPAR)

      ! Local variables
      INTEGER                :: I, J, N

      ! ADJUST_TO_TOMSAI begins here!
      WRITE( 6, '(a)' ) 'BIOBURN: Adjusting to TOMS AI data...'

      ! Loop over all tracers & boxes -- adjust to TOMS Aerosol index
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO N = 1, NBIOTRCE
         CALL TOMSAI( I, J, BIOMASS_ANN(N,I,J), 
     &                      BIOMASS_SEA(N,I,J), BIOMASS(N,I,J) )
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ADJUST_TO_TOMSAI

!------------------------------------------------------------------------------

      SUBROUTINE TOMSAI( I, J, VAL_ANN, VAL_SEAS, ADJUSTED_VALUE )
!
!******************************************************************************
!  Subroutine TOMSAI uses TOMS aerosol index for the last two decades as a 
!  surrogate for biomass burning.  The biomass burning emission climatology 
!  is adjusted for each month and year.  For months without information,
!  the climatology is used.  There is no TOMS AI data for July-August 1990 
!  and May 1993 - August 1996. 
!
!  Written by Bryan Duncan 8/2000.
!  Inserted into F90 module "biomass_mod.f" (bmy, 9/25/00, 7/20/04)
!
!  Subroutine TOMSAI is called from routine BIOBURN of "biomass_mod.f".
!
!  Arguments as Input:
!  ===========================================================================
!  (1-2) I, J           (INTEGER) : indices of box
!  (3  ) VAL_SEAS       (REAL*4 ) : Seasonal biomass value
!  (4  ) VAL_ANN        (REAL*4 ) : Annual biomass value
!
!  Arguments as Output:
!  ===========================================================================
!  (5  ) ADJUSTED_VALUE (REAL*4)  : CO emission for box(I,J) after adjustment.
!
!
!  Other variables:
!  ===========================================================================
!  TOMSAISCALE    = scaling factor by region for a specific month and year.
!  NAIREGIONS     = number of regions for which there is data.
!  NAIYEARS       = number of years for which there is data.
!  NAIMONTHS      = 12*NAIYEARS; number of months for which there is data.
!
!  NOTES:
!  (1 ) Remove references to "CMN_CO", "CMN_OH", and "CMN". (bmy, 9/25/00)
!  (2 ) Updated lat/lon boundaries of geographic regions (bnd, bmy, 10/16/00)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (4 ) Now use functions GET_MONTH, GET_TAU, GET_YEAR from "time_mod.f" 
!        Removed IMONTH from the arg list.  IMONTH, JYEAR, and TAU are now
!        local variables. (bmy, 2/11/03)
!  (5 ) Change VAL_ANN and VAL_SEAS to INTENT(IN). (bmy, 4/28/03)
!  (6 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TAU, GET_YEAR

#     include "CMN_SIZE"   ! Size parameters
!------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN_SETUP"  ! DATA_DIR
!------------------------------------------------
  
      ! Arguments
      INTEGER, INTENT(IN)    :: I, J
      REAL*8,  INTENT(INOUT) :: ADJUSTED_VALUE
      REAL*8,  INTENT(IN)    :: VAL_ANN
      REAL*8,  INTENT(IN)    :: VAL_SEAS

      ! Local variables
      INTEGER                :: THEYEAR, IPICK, CTM_lat, CTM_lon
      INTEGER                :: II, JJ, KK, LL, AS, IMONTH, JYEAR
      INTEGER, SAVE          :: IFIRSTCALL = 1 
      REAL*8                 :: CONVERT_lon, READINTOMS(NMONTHSAI), TAU

      !=================================================================
      ! TOMSAI begins here!  
      !=================================================================

      ! Get time quantities
      IMONTH = GET_MONTH()
      JYEAR  = GET_YEAR()
      TAU    = GET_TAU()
      
      !=================================================================
      ! Read in scaling factors on first call to SR.
      ! The scaling factors are stored in TOMSAI.  
      ! They run from Jan 1979 to Dec 1999 = 252 months total.
      !=================================================================
      IF(IFIRSTCALL.EQ.1) THEN
         IFIRSTCALL = 0
         
         ! Allocate TOMSAISCALE array
         ALLOCATE( TOMSAISCALE( NAIREGIONS, NAIYEARS, 12 ), STAT=AS )
         IF ( AS / = 0 ) CALL ALLOC_ERR( 'TOMSAISCALE' )

         ! Read TOMS Aerosol index data
         OPEN( 199, FILE = TRIM( DATA_DIR ) // 'TOMSAI', STATUS='OLD' )
         DO JJ=1,NAIREGIONS
            READ( 199, * ) readinTOMS
            
            II=0
            DO KK=1,NAIYEARS
               DO LL=1,12
                  II=II+1
                  TOMSAISCALE(JJ,KK,LL)=readinTOMS(II)
               ENDDO
            ENDDO

         ENDDO
         CLOSE(199)
      ENDIF

      !=================================================================
      ! The AI data is on a 1.25 x 1 degree grid (lon,lat).  Therefore,
      ! convert the box number from the code to the corresponding box
      ! number of the AI data.
      !=================================================================
#if   defined( GRID4x5 )
      CONVERT_lon = ( DBLE(I) * 5.d0 ) * 1.d0 / 1.25d0
      CTM_lon     = INT( CONVERT_lon )
      CTM_lat     = ( J * 4 ) - 2

      IF (J == 1     ) CTM_LAT = 2
      IF (J == JJPAR ) CTM_LAT = 88 + 90

#elif defined( GRID2x25 )
      CONVERT_lon = ( DBLE(I) * 2.5d0 ) * 1.d0 / 1.25d0
      CTM_LON     = INT( CONVERT_LON )
      CTM_LAT     = ( J * 2 ) - 1

      IF (J == 1     ) CTM_LAT = 1
      IF (J == JJPAR ) CTM_LAT = 89 + 90

#elif defined( GRID1x1 )
      PRINT*, 'Need to compute CONVERT_LON for 1 x 1 grid!'
      PRINT*, 'STOP in TOMSAI (biomass_mod.f)'
      STOP

#endif

      !=================================================================
      ! See what region the box falls in to pick the appropriate 
      ! regional scaling factor.
      !=================================================================
      IPICK=0

      ! Indonesia
      IF(CTM_lat.GE.83.and.CTM_lat.LE.99) THEN
         IF(CTM_lon.GE.221.and.CTM_lon.LE.269) THEN
            IPICK=1
         ENDIF
      ENDIF
      
      ! Brazil
      IF(CTM_lat.GE.59.and.CTM_lat.LE.91) THEN
        IF(CTM_lon.GE.96.and.CTM_lon.LE.116) THEN
           IF(IMONTH.GE.6.AND.IMONTH.LE.12) THEN
               IPICK=2
            ELSE
               IPICK=20
            ENDIF
         ENDIF
      ENDIF

      ! Southern Africa
      IF(CTM_lat.GE.50.and.CTM_lat.LE.90) THEN
         IF(CTM_lon.GE.128.and.CTM_lon.LE.184) THEN
            IPICK=3
         ENDIF
      ENDIF

      ! Northern Africa
      IF(CTM_lat.GE.91.and.CTM_lat.LE.110) THEN
         IF(CTM_lon.GE.128.and.CTM_lon.LE.184) THEN
            IPICK=4
         ENDIF
      ENDIF

      ! Central America and Mexico
      IF(CTM_lat.GE.96.and.CTM_lat.LE.115) THEN
         IF(CTM_lon.GE.61.and.CTM_lon.LE.85) THEN
            IF(IMONTH.GE.2.AND.IMONTH.LE.5) THEN
               IPICK=5
            ELSE
               IPICK=20
            ENDIF
         ENDIF
      ENDIF
      
      ! Canada and Alaska
      ! We have fire burn estimates for Canada, so we can use
      ! this data to fill in the TOMS data gap.
      IF(CTM_lat.GE.141.and.CTM_lat.LE.161) THEN
        IF(CTM_lon.GE.16.and.CTM_lon.LE.96) THEN
           IF(IMONTH.GE.5.AND.IMONTH.LE.9) THEN
              IPICK=6
           ELSE
              IPICK=20
           ENDIF
        ENDIF
      ENDIF
      
      ! Asiatic Russia
      IF(CTM_lat.GE.136.and.CTM_lat.LE.161) THEN
         IF(CTM_lon.GE.211.and.CTM_lon.LE.291) THEN
            IF(IMONTH.GE.5.AND.IMONTH.LE.9) THEN
               IPICK=7
            ELSE
               IPICK=20
            ENDIF
         ENDIF
      ENDIF

      ! Southeast Asia
      IF(CTM_lat.GE.99.and.CTM_lat.LE.119) THEN
         IF(CTM_lon.GE.221.and.CTM_lon.LE.239) THEN
            IF(IMONTH.GE.1.AND.IMONTH.LE.5) THEN
               IPICK=8
            ELSE
               IPICK=20
            ENDIF
         ENDIF
      ENDIF
      
      ! Error Check.
      IF(IMONTH.LT.1.OR.IMONTH.GT.12) THEN
         PRINT*,'Error in SR TOMSAI:  Value of IMONTH is wrong.'
         PRINT*,'IMONTH = ',IMONTH
         STOP
      ENDIF

      !=================================================================
      ! During the TOMS data gaps, set IPICK = 0; emissions are
      ! not rescaled for the box and the seasonal variation is used,
      ! except for Indonesia and Canada & Alaska.
      !=================================================================

      ! July - August 1990
      IF ( IPICK /= 6 ) THEN
         IF ( TAU == 48168d0 ) IPICK = 0
         IF ( TAU == 48912d0 ) IPICK = 0
      ENDIF

      ! May 1993 - July 1996
      IF ( IPICK /= 6 .AND. IPICK /= 1 ) THEN
         IF ( TAU >= 73008d0 .AND. TAU <= 100776d0 ) IPICK = 0
      ENDIF

      !=================================================================
      ! Rescale CO emission.  If IPICK = 0 then emissions are
      ! not rescaled for the box and the seasonal variation is used.
      !=================================================================

      ! Adjust with TOMS AI
      IF( IPICK > 0 .AND. IPICK /= 3 .AND. IPICK /= 4 ) THEN

         THEYEAR = JYEAR - 1978
     
         IF ( THEYEAR > NAIYEARS .OR. THEYEAR < 0 ) THEN
            PRINT*,'Error in SR TOMSAI:  You have picked a year less'
            PRINT*,'than 1979 or greater than 1999.  The data in this'
            PRINT*,'SR used to scale biomass burning emissions is only'
            PRINT*,'good for 1979-1999. You may need to comment out'
            PRINT*,'the call to this SR in SR bioburn.'
            PRINT*,'Your year is ',JYEAR,'.'
            STOP
         ENDIF

         ! Do not adjust Africa with TOMSAI!!!!
         IF ( IPICK /= 20 ) THEN
            ADJUSTED_VALUE = VAL_ANN * 
     &                       TOMSAISCALE(IPICK,THEYEAR,IMONTH)
         ELSE

            ! Zero out IPICK regions when the biomass burning 
            ! season is not occuring.
            ADJUSTED_VALUE = 0d0
         ENDIF

      ELSE  ! IPICK=0; IPICK=3; IPICK=4 

         ! Use seasonal emissions instead
         ADJUSTED_VALUE = VAL_SEAS
      ENDIF

      ! Return to calling program
      END SUBROUTINE TOMSAI

!------------------------------------------------------------------------------
      
      SUBROUTINE SET_BIOTRCE
!
!******************************************************************************
!  Subroutine SET_BIOTRCE intializes the BIOTRCE array and the NBIOTRCE number
!  of biomass tracers.  This was split off from "tracerid.f" in order to 
!  prevent circular module references. (bmy, 11/6/02)
!
!  NOTES:
!  (1 ) SET_BIOTRCE is called from "input.f", after the call to TRACERID.
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD

      !=================================================================
      ! SET_BIOTRCE begins here!
      !=================================================================

      ! Initialize
      NBIOTRCE = 0
      
      ! Increment NBIOTRCE for each turned on biomass tracer
      IF ( IDBNOX  /= 0 ) NBIOTRCE = NBIOTRCE + 1
      IF ( IDBCO   /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBALK4 /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBACET /= 0 ) NBIOTRCE = NBIOTRCE + 1
      IF ( IDBMEK  /= 0 ) NBIOTRCE = NBIOTRCE + 1
      IF ( IDBALD2 /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBPRPE /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC3H8 /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBCH2O /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC2H6 /= 0 ) NBIOTRCE = NBIOTRCE + 1 

      ! Fill BIOTRCE w/ appropriate TRACER ID #'s
      IF ( IDBNOX  /= 0 ) BIOTRCE(IDBNOX ) = IDTNOX
      IF ( IDBCO   /= 0 ) BIOTRCE(IDBCO  ) = IDTCO
      IF ( IDBALK4 /= 0 ) BIOTRCE(IDBALK4) = IDTALK4
      IF ( IDBACET /= 0 ) BIOTRCE(IDBACET) = IDTACET
      IF ( IDBMEK  /= 0 ) BIOTRCE(IDBMEK ) = IDTMEK
      IF ( IDBALD2 /= 0 ) BIOTRCE(IDBALD2) = IDTALD2
      IF ( IDBPRPE /= 0 ) BIOTRCE(IDBPRPE) = IDTPRPE
      IF ( IDBC3H8 /= 0 ) BIOTRCE(IDBC3H8) = IDTC3H8
      IF ( IDBCH2O /= 0 ) BIOTRCE(IDBCH2O) = IDTCH2O  
      IF ( IDBC2H6 /= 0 ) BIOTRCE(IDBC2H6) = IDTC2H6  
      
      ! Echo biomass burning tracers 
      WRITE( 6, 100 ) BIOTRCE( 1:NBIOTRCE )
 100  FORMAT( 'TRACERID: Biomass burning tracers        :', 20i3 )

      ! Return to calling program
      END SUBROUTINE SET_BIOTRCE

!------------------------------------------------------------------------------

      SUBROUTINE INIT_BIOMASS
!
!******************************************************************************
!  Subroutine INIT_BIOMASS allocates and zeroes the BURNEMIS array.
!  (bmy, 9/11/00, 7/20/04)
!
!  NOTES:
!  (1 ) BURNEMIS does not need to be allocated if all of the elements of 
!        LDOBIOEMIT = .FALSE. or if biomass burning is not turned on.  
!        This will save memory from being allocated needlessly. (bmy, 9/11/00)
!  (2 ) Now use NBIOTRCE > 0 as an extra criterion for allocating 
!        the BURNEMIS array (bmy, 4/16/01)
!  (3 ) BURNEMIS is now allocated to be of size (NBIOTRCE,IIPAR,JJPAR).
!        (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (5 ) Renamed to "INIT_BIOMASS".  Also allocate the BIOMASS array here.
!        (bmy, 5/28/02)
!  (6 ) Now references ALLOC_ERR from "error_mod.f".  Now reference IDBNOX,
!        IDBCO, etc from "tracerid_mod.f". (bmy, 11/6/02)
!  (7 ) Now references LBIOMASS from "logical_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LBIOMASS
      USE TRACERID_MOD

#     include "CMN_SIZE"  ! Size parameters, etc
!----------------------------------------------------
! Prior to 7/20/04:
!#     include "CMN"       ! LBIONOX
!----------------------------------------------------

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_BIOMASS begins here!
      !=================================================================
      !------------------------------------------
      ! Prior to 7/20/04:
      !IF ( LBIONOX .and. NBIOTRCE > 0 ) THEN 
      !------------------------------------------
      IF ( LBIOMASS .and. NBIOTRCE > 0 ) THEN 

         ! Allocate and zero BURNEMIS array [molec/cm3/s]
         ALLOCATE( BIOMASS( NBIOTRCE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS' )
         BIOMASS = 0d0

         ! Allocate and zero BURNEMIS array [molec/cm3/s]
         ALLOCATE( BURNEMIS( NBIOTRCE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BURNEMIS' )
         BURNEMIS = 0d0

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE INIT_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_BIOMASS deallocates the BURNEMIS and
!  TOMSAISCALE arrays (bmy, 9/11/00, 3/19/01)
!
!  NOTES:
!  (1) Also deallocate TOMSAISCALE (bmy, 9/25/00)
!  (2) LDOBIOEMIT, BIOTRCE are no longer dynamically allocatable (bmy, 3/19/01)
!******************************************************************************
!      
      IF ( ALLOCATED( BURNEMIS    ) ) DEALLOCATE( BURNEMIS    )
      IF ( ALLOCATED( TOMSAISCALE ) ) DEALLOCATE( TOMSAISCALE )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOMASS

!------------------------------------------------------------------------------
      
      END MODULE BIOMASS_MOD
