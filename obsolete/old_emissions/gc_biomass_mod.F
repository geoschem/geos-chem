! $Id: gc_biomass_mod.f,v 1.4 2010/03/15 19:33:23 ccarouge Exp $
      MODULE GC_BIOMASS_MOD
!
!******************************************************************************
!  Module GC_BIOMASS_MOD contains arrays and routines to compute monthly
!  biomass burning emissions for NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, 
!  C3H8, CH2O, C2H6, CH4, and CH3I. (bmy, 9/11/00, 1/13/10)
!
!  NOTE: These biomass emissions are based on Bryan Duncan (Duncan et al 2001)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NBIOMAX              : maximum # of biomass burning tracers
!  (2 ) BIOTRCE              : index array of biomass burning tracers
!  (4 ) BIOMASS              : array of biomass burning emissions [molec/cm3/s]
!  (5 ) BIOMASS_SAVE         : array of biomass burning emissions [molec/cm2/s]
!  (6 ) TOMSAISCALE          : array for TOMS aerosol index values
!
!  Module Routines:
!  ============================================================================
!  (1 ) GC_COMPUTE_BIOMASS   : reads data, computes gas-phase biomass emissions
!  (2 ) GC_READ_BIOMASS_BCOC : reads biomass emissions of BC & OC
!  (3 ) GC_READ_BIOMASS_CO2  : reads biomass emissions of CO2
!  (4 ) GC_READ_BIOMASS_NH3  : reads biomass emissions of NH3
!  (5 ) GC_READ_BIOMASS_SO2  : reads biomass emissions of SO2
!  (6 ) READ_BIOMASS         : reads gas-phase biomass burning data from disk
!  (7 ) SCALE_BIOMASS_ACET   : applies scale factors to ACET 
!  (8 ) SCALE_FUTURE         : applies future scale factors to emissions
!  (9 ) TOTAL_BIOMASS_TG     : prints monthly emission totals in [Tg (C)]
!  (10) ADJUST_TO_TOMSAI     : wrapper for subroutine TOMSAI
!  (11) TOMSAI               : adjusts BB for int'annual var'bilty w/ TOMS data
!  (12) CLEANUP_BIOMASS      : deallocates BURNEMIS, BIOTRCE
!
!  GEOS-CHEM modules referenced by "gc_biomass_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f          : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f            : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f           : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f      : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f          : Module w/ I/O error and NaN check routines
!  (6 ) grid_mod.f           : Module w/ horizontal grid information
!  (7 ) logical_mod.f        : Module w/ GEOS-CHEM logical switches
!  (8 ) time_mod.f           : Module w/ routines for computing time & date
!  (9 ) transfer_mod.f       : Module w/ routines to cast & resize arrays
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
!  (1 ) Now treat BURNEMIS as a true global array of size (IIPAR,JJPAR);
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
!  (27) Bug fix in BIOBURN for TAU w/ interannual emissions (bmy, 3/18/05)
!  (28) Now can read data for both GEOS and GCAP grids (bmy, 3/18/05)
!  (29) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (30) Renamed to "gc_biomass_mod.f", so that we can use either these
!        "default" biomass emissions or GFED2 biomass emissions.  Cleaned up
!        a lot of obsolete stuff. (bmy, 4/5/06)
!  (31) Modified for IPCC future emissions scale factors.  Added private 
!        routine SCALE_FUTURE. (swu, bmy, 5/30/06)
!  (32) Added routines for reading BC, OC, SO2, NH3, CO2 biomass emissions.
!        (bmy, 9/28/06)
!  (33) Add 9 gaseous biomass burning emissions using emission ratios
!       w.r.t. CO.  Details in Fu et al. [2008] (tmf, 1/7/09)
!  (34) CO scaling for VOC production is transfered to biomass_mod.f.
!        (jaf, 2/6/09)
!  (35) Modification for 1 x 1.25 grid (lok, bmy, 1/13/10)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gc_biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except for these routines
      PUBLIC :: CLEANUP_GC_BIOMASS
      PUBLIC :: GC_COMPUTE_BIOMASS
      PUBLIC :: GC_READ_BIOMASS_BCOC
      PUBLIC :: GC_READ_BIOMASS_CO2
      PUBLIC :: GC_READ_BIOMASS_NH3
      PUBLIC :: GC_READ_BIOMASS_SO2
      PUBLIC :: TOTAL_BIOMASS_TG
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      ! NOTE: This is an internal declaration for the 
      ! gas-phase species only (bmy, 9/28/06)
      ! use NBIOMAX from CMN_SIZE_mod now (hotp 8/3/09)
      !INTEGER, PARAMETER   :: NBIOMAX = 19

      ! TOMS AI interannual variability in biomass burning emissions
      INTEGER, PARAMETER   :: NAIREGIONS = 8
      INTEGER, PARAMETER   :: NAIYEARS   = 21
      INTEGER, PARAMETER   :: NMONTHSAI  = NAIYEARS * 12 

      ! Arrays
      REAL*8,  ALLOCATABLE :: TOMSAISCALE(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GC_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )
!
!******************************************************************************
!  Subroutine GC_COMPUTE_BIOMASS computes the biomass burning emissions for 
!  several species for the given month (jal, acs, rvm, bmy, 9/11/00, 4/5/06)
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
!  (23) Bug fix: if using interannual biomass emissions then get the TAU value
!         for the first of the current month & year.  This will make sure that
!         runs which start mid-month will access the biomass data correctly.
!         (bmy, 3/18/05)
!  (24) Now can read data from both GEOS and GCAP grids (bmy, 8/16/05)
!  (25) Renamed to GC_COMPUTE_BIOMASS.  Now takes YEAR, MONTH, BIOMASS 
!        arguments. (bmy, 4/5/06)
!  (26) Add 9 biomass burning species (ccc, 1/7/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LBBSEA,          LTOMSAI
      USE LOGICAL_MOD,   ONLY : LFUTURE
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR,  GET_TAU

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments 
      INTEGER, INTENT(IN)    :: YEAR, MONTH
      REAL*8,  INTENT(OUT)   :: BIOMASS(IIPAR,JJPAR,NBIOMAX)

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, N
      INTEGER, SAVE          :: MONTHSAVE = -99
      REAL*8                 :: TIME, XTAU
      REAL*8                 :: BIOMASS_SEA(IIPAR,JJPAR,NBIOMAX) 
      REAL*8                 :: BIOMASS_ANN(IIPAR,JJPAR,NBIOMAX)
      CHARACTER(LEN=4  )     :: CYEAR
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=255)     :: BIOMASS_DIR

      ! MONTHDATES = number of days per month
      INTEGER                :: MONTHDATES(12) = (/ 31, 28, 31, 30, 
     &                                              31, 30, 31, 31, 
     &                                              30, 31, 30, 31 /)

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
      !      GLYX      15         55        [molec GLYX/cm2/month]
      !      MGLY      16         56        [molec MGLY/cm2/month]
      !      BENZ      17         57        [molec C   /cm2/month]
      !      TOLU      18         58        [molec C   /cm2/month]
      !      XYLE      19         59        [molec C   /cm2/month]
      !      C2H4      20         63        [molec C   /cm2/month]
      !      C2H2      21         64        [molec C   /cm2/month]
      !      GLYC      22         66        [molec GLYC/cm2/month]
      !      HAC       23         67        [molec HAC /cm2/month]
      !
      ! Subsequent unit conversion is done on the following species:
      !      [molec ACET/cm2/month]  -->  [molec C/cm2/month]
      !      [molec C3H8/cm2/month]  -->  [molec C/cm2/month]
      !      [molec C2H6/cm2/month]  -->  [molec C/cm2/month]
      ! 
      ! There are 19 biomass burning species in this module.
      !
      ! Biomass burning emissions are first read from disk into the
      ! BIOMASS array.  After unit conversion to [molec/cm3/s] ( or
      ! [atoms C/cm3/s] for hydrocarbons), the emissions are stored 
      ! in BIOMASS and passed back to the calling program.
      !
      ! Biomass burning data is monthly, so we only have to read 
      ! emissions from disk once each month.
      !=================================================================

      ! Do the following on the first day of a new month...
      IF ( MONTH /= MONTHSAVE ) THEN   
         
         ! Save the current month
         MONTHSAVE = MONTH

         ! Set MONTHDATES(2) = 29 for leapyears, = 28 otherwise (bmy, 4/19/99)
         IF ( MONTH == 2 ) THEN
            IF( ITS_A_LEAPYEAR() ) THEN
               MONTHDATES(2) = 29
            ELSE
               MONTHDATES(2) = 28
            ENDIF
         ENDIF

         ! TIME = conversion from [molec/cm2/month] to [molec/cm2/s]
         TIME = ( DBLE( MONTHDATES( MONTH ) ) * 86400d0 )

         ! Create a string for the 4-digit year
         WRITE( CYEAR, '(i4)' ) YEAR

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
            XTAU = GET_TAU0( MONTH, 1, 1985 )

            ! Filename for seasonal biomass burning emissions
            FILENAME = TRIM( DATA_DIR )    // TRIM( BIOMASS_DIR ) // 
     &                 'bioburn.seasonal.' // GET_NAME_EXT_2D()   //
     &                 '.'                 // GET_RES_EXT()

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
            XTAU = GET_TAU0( MONTH, 1, 1985 )            

            ! Filename for seasonal biomass burning emissions
            FILENAME = TRIM( DATA_DIR )    // TRIM( BIOMASS_DIR ) //
     &                 'bioburn.seasonal.' // GET_NAME_EXT_2D()   //
     &                 '.'                 // GET_RES_EXT()

            ! Read the seasonal biomass burning emissions into BIOMASS_SEA
            CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS_SEA )

            ! Get TAU0 value to index the punch file -- use generic year 1985
            XTAU = GET_TAU0( 1, 1, 1985 )

            ! Filename for annual biomass burning emissions 
            FILENAME = TRIM( DATA_DIR )  // TRIM( BIOMASS_DIR ) //
     &                 'bioburn.annual.' // GET_NAME_EXT_2D()   //
     &                 '.'               // GET_RES_EXT()

            ! Read the annual biomass burning emissions into BIOMASS_ANN
            CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS_ANN )

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
               XTAU = GET_TAU0( MONTH, 1, 1985 )

               ! Filename for seasonal biomass burning emissions
               FILENAME = TRIM( DATA_DIR )    // TRIM( BIOMASS_DIR ) //
     &                    'bioburn.seasonal.' // GET_NAME_EXT_2D()   //
     &                    '.'                 // GET_RES_EXT() 

               ! Read the seasonal biomass burning emissions into BIOMASS_SEA
               CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS_SEA )

               ! Get TAU0 value to index the punch file -- 
               ! use generic year 1985
               XTAU = GET_TAU0( 1, 1, 1985 )

               ! Filename for annual biomass burning emissions
               FILENAME = TRIM( DATA_DIR )  // TRIM( BIOMASS_DIR ) //
     &                    'bioburn.annual.' // GET_NAME_EXT_2D()   //
     &                    '.'               // GET_RES_EXT()

               ! Read the annual biomass burning emissions from disk
               CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS_ANN )

               ! Adjust the annual biomass burning to the TOMS Aerosol 
               ! Index data where necessary.  Otherwise, overwrite
               ! with seasonal data. Save result in the BIOMASS array.
               BIOMASS = 0d0
               CALL ADJUST_TO_TOMSAI( BIOMASS_ANN, BIOMASS_SEA, BIOMASS)

            ELSE

               ! Use actual TAU0 value to index punch file
               XTAU = GET_TAU()

               ! Filename for interannual variability biomass burning emissions
               FILENAME = TRIM( DATA_DIR )       // 
     &                    TRIM( BIOMASS_DIR )    //
     &                    'bioburn.interannual.' // GET_NAME_EXT_2D() // 
     &                    '.'                    // GET_RES_EXT()     // 
     &                    '.'                    // CYEAR           

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

            ! TAU0 value for 0 GMT on the first day of this month & year
            XTAU = GET_TAU0( MONTH, 1, YEAR )

            ! Filename for interannual variability biomass burning emissions
            FILENAME = TRIM( DATA_DIR )       // 
     &                 TRIM( BIOMASS_DIR )    //      
     &                 'bioburn.interannual.' // GET_NAME_EXT_2D() // 
     &                 '.'                    // GET_RES_EXT()     // 
     &                 '.'                    // CYEAR
               
            ! Read interannual variability biomass burning
            CALL READ_BIOMASS( FILENAME, XTAU, BIOMASS )

         ENDIF


         ! Convert to [molec/cm2/s] or [atoms C/cm2/s]
         BIOMASS = BIOMASS / TIME
 
         ! Fancy output...
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF

      ! Return to calling program
      END SUBROUTINE GC_COMPUTE_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE READ_BIOMASS( FILENAME, TAU0, BIOMASS )
!
!******************************************************************************
!  Subroutine READ_BIOMASS reads the biomass burning emissions from disk
!  in units of [molec/cm2/month] (or [atoms C/cm2/month] for hydrocarbons). 
!  (bmy, 9/25/00, 5/30/06)
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
!        (IIPAR,JJPAR) to (IIPAR,JJPAR).  (bmy, 9/28/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) BIOMASS needs to be of size (NBIOTRCE,IIPAR,JJPAR) (bmy, 5/31/02)
!  (10) Now references IDTNOX, etc. from "tracerid_mod.f" (bmy, 11/6/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. flag to suppress extra info 
!        from being printed (bmy, 3/14/03)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now make BIOMASS array argument (I,J,N) ordered instead of (N,I,J).
!        Also now read all NBIOMAX species.  (bmy, 4/5/06)
!  (14) Now refrerences LFUTURE from "logical_mod.f".  Also now calls private
!        routine SCALE_FUTURE to compute the future biomass emissions.
!        (swu, bmy, 5/30/06)
!  (15) Use IDBs from TRACERID_MOD. Replace hard-wired biomass # by
!       IDBs. (hotp 8/3/09)
!  (16) Treat all species with emissions dependant on CO emissions after
!       the main loop on biomass species. To make sure CO emissions are already
!       read. (ccc, 01/29/10)
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : READ_BPCH2
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      ! Add IDBs of all species handled here (hotp 8/3/09)
      USE TRACERID_MOD, ONLY : IDBNO,   IDBCO,   IDBALK4, IDBACET 
      USE TRACERID_MOD, ONLY : IDBMEK,  IDBALD2, IDBPRPE, IDBC3H8
      USE TRACERID_MOD, ONLY : IDBCH2O, IDBC2H6
      USE TRACERID_MOD, ONLY : IDBGLYX, IDBMGLY, IDBBENZ, IDBTOLU
      USE TRACERID_MOD, ONLY : IDBXYLE, IDBC2H4, IDBC2H2, IDBGLYC
      USE TRACERID_MOD, ONLY : IDBHAC
     
      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN)  :: FILENAME
      REAL*8,           INTENT(IN)  :: TAU0
      REAL*8,           INTENT(OUT) :: BIOMASS(IIPAR,JJPAR,NBIOMAX)

      ! Local variables
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)

      ! Add storage of CO emissions to calculate emissions 
      ! of the 9 new species (tmf, 12/18/06) 
      REAL*8   :: COEMIS(IIPAR, JJPAR)  ! CO emissions before scaling
      REAL*8   :: TRCEMIS(IIPAR, JJPAR) ! Tracer emissions scaled from CO

      !=================================================================
      ! READ_BIOMASS begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( 'GC_COMPUTE_BIOMASS: Reading ', a )

      ! Initialize the BIOMASS array
      BIOMASS = 0d0

      ! We don't need the loop any more. We calculate emissions if IDBxxx > 0
      ! (ccc, 01/29/10)
          
         ! Do scaling if necessary and print totals in Tg
         ! Replaced NOx emissions with NO emissions (mpayer, 3/14/13)
         IF  ( IDBNO > 0 ) THEN

            !----------
            ! NOx
            !----------

            ! NOx is stored in the biomass file as tracer #1
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 1, 
     &                       TAU0,      IIPAR,     JJPAR,     
     &                       1,         ARRAY,     QUIET=.TRUE. )  

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBNO) )

            ! Compute future NOx emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'NObb', BIOMASS(:,:,IDBNO) )
            ENDIF

            ! NO -- print totals in [Tg/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBNO), 14d-3, 'NO' )
         ENDIF

         !FP_ISOP (6/2009) removed hardwired IDBs
         IF  ( IDBCO > 0 ) THEN
 
            !----------
            ! CO
            !----------

            ! CO is stored in the biomass file as tracer #4
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 4, 
     &                       TAU0,      IIPAR,     JJPAR,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBCO) )

            ! Store CO emissions before scaling 
            ! for new gaseous emissions (tmf, 1/7/09)
            CALL TRANSFER_2D( ARRAY(:,:,1), COEMIS(:,:) )

            ! Compute future CO emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'CObb', BIOMASS(:,:,IDBCO) )
            ENDIF

            ! Print totals in [Tg/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBCO), 28d-3, 'CO' )
         ENDIF

         !FP_ISOP
         IF  ( IDBALK4 > 0 ) THEN

            !----------
            ! ALK4
            !----------

            ! ALK4 is stored in the biomass file as tracer #5
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 5, 
     &                       TAU0,      IIPAR,     JJPAR,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBALK4) )

            ! Compute future ALK4 emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBALK4) )
            ENDIF

            ! ALK4 -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBALK4), 12d-3, 'ALK4' )
         ENDIF

         !FP_ISOP
         IF  ( IDBACET > 0 ) THEN

            !----------
            ! ACET
            !----------

            ! ACET is stored in the biomass file as tracer #9
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 9, 
     &                       TAU0,      IIPAR,     JJPAR,     
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBACET) )

            ! ACET -- Convert from [molec/cm2/month] to [molec C/cm2/month] 
            BIOMASS(:,:,IDBACET) = BIOMASS(:,:,IDBACET) * 3d0  

            ! Scale to yearly value for biogenic acetone (bdf, bmy, 7/23/01)
            CALL SCALE_BIOMASS_ACET( BIOMASS(:,:,IDBACET) )

            ! Compute future ACET emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBACET) )
            ENDIF

            ! Print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBACET), 12d-3, 'ACET' )
         ENDIF

         !FP_ISOP
         IF  ( IDBMEK > 0 ) THEN

            !----------
            ! MEK
            !----------

            ! MEK is stored in the biomass file as tracer #10
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 10, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBMEK) )

            ! Compute future MEK emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBMEK) )
            ENDIF

            ! MEK -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBMEK), 12d-3, 'MEK' )
         ENDIF

         !FP_ISOP
         IF  ( IDBALD2 > 0 ) THEN

            !----------
            ! ALD2
            !----------

            ! ALD2 is stored in the biomass file as tracer #11
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 11, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBALD2) )

            ! Compute future ALD2 emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBALD2) )
            ENDIF

            ! ALD2 -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBALD2), 12d-3, 'ALD2' )
         ENDIF

         !FP_ISOP
         IF  ( IDBPRPE > 0 ) THEN

            !----------
            ! PRPE
            !----------

            ! PRPE is stored in the biomass file as tracer #18
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 18, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBPRPE) )

            ! Compute future PRPE emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBPRPE) )
            ENDIF

            ! PRPE -- convert from [molec/cm2/month] to [molec C/cm2/month]
            ! Print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBPRPE), 12d-3, 'PRPE' )
         ENDIF

         !FP_ISOP
         IF  ( IDBC3H8 > 0 ) THEN
               
            !----------
            ! C3H8
            !----------

            ! C3H8 is stored in the biomass file as tracer #19
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 19, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBC3H8) )

            ! C3H8 -- convert from [molec/cm2/month] to [molec C/cm2/month] 
            BIOMASS(:,:,IDBC3H8) = BIOMASS(:,:,IDBC3H8) * 3d0 

            ! Compute future C3H8 emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBC3H8) )
            ENDIF

            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBC3H8), 12d-3, 'C3H8' )
         ENDIF

         !FP_ISOP
         IF  ( IDBCH2O > 0 ) THEN

            !----------
            ! CH2O
            !----------

            ! CH2O is stored in the biomass file as tracer #20
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 20, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBCH2O) )

            ! Compute future CH2O emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBCH2O) )
            ENDIF

            ! CH2O -- print totals in [Tg C/month]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBCH2O), 30d-3, 'CH2O' )
         ENDIF

         !FP_ISOP
         IF  ( IDBC2H6 > 0 ) THEN

            !----------
            ! C2H6
            !----------

            ! C2H6 is stored in the biomass file as tracer #21
            CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 21, 
     &                       TAU0,      IIPAR,     JJPAR,      
     &                       1,         ARRAY,     QUIET=.TRUE. ) 

            ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
            CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS(:,:,IDBC2H6) )

            ! C2H6 --convert from [molec/cm2/month] to [molec C/cm2/month]
            BIOMASS(:,:,IDBC2H6) = BIOMASS(:,:,IDBC2H6) * 2d0 

            ! Compute future C2H6 emissions (if necessary)
            IF ( LFUTURE ) THEN
               CALL SCALE_FUTURE( 'VOCbb', BIOMASS(:,:,IDBC2H6) )
            ENDIF

            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBC2H6), 12d-3, 'C2H6' )

         ENDIF
!      ENDDO

!-------------------------------------------------------------------------
!  Add 9 gaseous BB emissions (tmf, 1/7/09)
!  These emissions depend on CO emissions. We have to make sure CO
!  emissions are compiled, so they are done at the end.
!  (ccc, 01/27/10)
!-------------------------------------------------------------------------
         ! Do scaling if necessary and print totals in Tg
         IF  ( IDBGLYX > 0 ) THEN

            !----------
            ! GLYX
            !----------

            ! Estimate GLYC emission by scaling CO emission
            ! GLYX [mole] / CO [mole]  = 0.00662  (from Andreae 2005 update)
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00662d0     ! [molecule GLYX/cm2/month]
            
            BIOMASS(:,:,IDBGLYX) = TRCEMIS(:,:)

            ! GLYX -- [molecule GLYX/cm2/month]
            ! Print totals in [Tg GLYX]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBGLYX), 58d-3, 'GLYX' )
         ENDIF

         IF  ( IDBMGLY > 0 ) THEN

            !----------
            ! MGLY
            !----------

            ! Estimate MGLY emission by scaling CO emission
            ! MGLY [mole] / CO [mole]  = 0.00347  (from Andreae 2005 update)
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00347d0     ! [molecule MGLY/cm2/month]
            
            BIOMASS(:,:,IDBMGLY) = TRCEMIS(:,:)

            ! MGLY -- [molecule MGLY/cm2/month]
            ! Print totals in [Tg MGLY]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBMGLY), 72d-3, 'MGLY' )
         ENDIF
         
         IF  ( IDBBENZ > 0 ) THEN

            !----------
            ! BENZ
            !----------
            ! Estimate BENZ emission by scaling CO emission
            ! BENZ [mole] / CO [mole]  =  0.00233
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00233d0 * 6.0d0   ! [molec C/cm2/month]
            
            BIOMASS(:,:,IDBBENZ) = TRCEMIS(:,:)

            ! BENZ -- [molec C/cm2/month]
            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBBENZ), 12d-3, 'BENZ' )

         ENDIF

         IF  ( IDBTOLU > 0 ) THEN

            !----------
            ! TOLU
            !----------
            ! Estimate TOLU emission by scaling CO emission
            ! TOLU [mole] / CO [mole]  =  0.00124
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00124d0 * 7.0d0   ! [molec C/cm2/month]
            
            BIOMASS(:,:,IDBTOLU) = TRCEMIS(:,:)

            ! TOLU -- [molec C/cm2/month]
            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBTOLU), 12d-3, 'TOLU' )
         ENDIF

         IF  ( IDBXYLE > 0) THEN

            !----------
            ! XYLE
            !----------
            ! Estimate XYLE emission by scaling CO emission
            ! XYLE [mole] / CO [mole]  =  0.00048
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00048d0 * 8.0d0   ! [molec C/cm2/month]
            
            BIOMASS(:,:,IDBXYLE) = TRCEMIS(:,:)

            ! XYLE -- [molec C/cm2/month]
            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBXYLE), 12d-3, 'XYLE' )
         ENDIF

         IF  ( IDBC2H4 > 0 ) THEN

            !----------
            ! C2H4
            !----------
            ! Estimate C2H4 emission by scaling CO emission
            ! C2H4 [mole] / CO [mole]  =  0.01381
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.01381d0 * 2.0d0   ! [molec C/cm2/month]
            
            BIOMASS(:,:,IDBC2H4) = TRCEMIS(:,:)

            ! C2H4 -- [molec C/cm2/month]
            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBC2H4), 12d-3, 'C2H4' )
         ENDIF

         IF  ( IDBC2H2 > 0) THEN

            !----------
            ! C2H2
            !----------
            ! Estimate C2H2 emission by scaling CO emission
            ! C2H2 [mole] / CO [mole]  =  0.004d0  from Xiao et al. [2007]
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.004d0 * 2.0d0   ! [molec C/cm2/month]
            
            BIOMASS(:,:,IDBC2H2) = TRCEMIS(:,:)

            ! C2H2 -- [molec C/cm2/month]
            ! Print totals in [Tg C]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBC2H2), 12d-3, 'C2H2' )
         ENDIF

         IF  ( IDBGLYC > 0 ) THEN

            !----------
            ! GLYC
            !----------
            ! Estimate GLYC emission by scaling CO emission
            ! GLYC [mole] / CO [mole]  = 0.00477  (from Andreae 2005 update)
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00477d0     ! [molecule GLYC/cm2/month]
            
            BIOMASS(:,:,IDBGLYC) = TRCEMIS(:,:)

            ! GLYC -- [molecule GLYC/cm2/month]
            ! Print totals in [Tg GLYC]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBGLYC), 60d-3, 'GLYC' )
         ENDIF

         IF  ( IDBHAC > 0 ) THEN

            !----------
            ! HAC
            !----------
            ! Estimate HAC emission by scaling CO emission
            ! HAC [mole] / CO [mole]  = 0.00331d0 (from Christian et al. [2003] for African biomass)
            TRCEMIS(:,:) = 
     &         COEMIS(:,:) * 0.00331d0     ! [molecule HAC/cm2/month]

            BIOMASS(:,:,IDBHAC) = TRCEMIS(:,:)

            ! HAC -- [molecule HAC/cm2/month]
            ! Print totals in [Tg HAC]
            CALL TOTAL_BIOMASS_TG( BIOMASS(:,:,IDBHAC), 74d-3, 'HAC' )
         ENDIF

     
      ! Return to calling program
      END SUBROUTINE READ_BIOMASS

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

      USE CMN_SIZE_MOD    ! Size parameters

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

      SUBROUTINE SCALE_FUTURE( NAME, BB )
!
!******************************************************************************
!  Subroutine SCALE_FUTURE applies the IPCC future emissions scale factors
!  to the biomass burning emisisons to compute the future emissions of biomass
!  burning for NOx, CO, and VOC's.  (swu, bmy, 5/30/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NAME (CHARACTER) : Denotes type of scale factor to use (e.g. NOx)
!  (2 ) BB   (REAL*8   ) : Array w/ biomass burning emissions [molec/cm2]
!
!  NOTES:
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbb

      USE CMN_SIZE_MOD               ! Size parameters

      ! Arguments
      REAL*8,           INTENT(INOUT) :: BB(IIPAR,JJPAR)
      CHARACTER(LEN=*), INTENT(IN)    :: NAME

      ! Local variables
      INTEGER                         :: I, J
      
      !=================================================================
      ! SCALE_FUTURE begins here!
      !=================================================================

      IF ( NAME == 'NObb' ) THEN

         ! Compute future NOx emissions
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            BB(I,J) = BB(I,J) * GET_FUTURE_SCALE_NOxbb( I, J )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE IF ( NAME == 'CObb' ) THEN

         ! Compute future CO emissions 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            BB(I,J) = BB(I,J) * GET_FUTURE_SCALE_CObb( I, J )
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
            BB(I,J) = BB(I,J) * GET_FUTURE_SCALE_VOCbb( I, J )
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_BIOMASS_TG( BBARRAY, MOLWT, NAME )
!
!******************************************************************************
!  Subroutine TOTAL_BIOMASS_TG prints the amount of biomass burning emissions 
!  that are emitted each month in Tg or Tg C. (bmy, 3/20/01, 4/5/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Biomass burning CO emissions [molec/cm2/month]
!
!  NOTES:
!  (1 ) BBARRAY is now dimensioned (IIPAR,JJPAR).  Also, DXYP is dimensioned
!        as JJPAR, so use J+J0 to reference it. (bmy, 9/28/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (3 ) Now use function GET_AREA_CM2 from "grid_mod.f" to compute grid
!        box surface area in cm2.  Removed reference to CMN header file.
!        Cosmetic changes. (bmy, 3/14/03)
!  (4 ) Now report sums of NOx as Tg N instead of Tg NOx (bmy, 4/5/06)
!  (5 ) Add unit choice for GLYX, MGLY, GLYC, HAC (tmf, 1/7/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,            ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*8,           INTENT(IN) :: BBARRAY(IIPAR,JJPAR) 
      REAL*8,           INTENT(IN) :: MOLWT
      CHARACTER(LEN=*), INTENT(IN) :: NAME

      ! Local variables
      INTEGER                      :: I, J
      REAL*8                       :: TOTAL, CONV
      CHARACTER(LEN=6)             :: UNIT

      !=================================================================
      ! TOTAL_BIOMASS_TG begins here!
      !=================================================================

      ! Initialize summing variable
      TOTAL = 0d0

      ! Convert to [Tg/month] (or [Tg C/month] for hydrocarbons)
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor to [Tg/month] (or [Tg C/month] for HC's)
         CONV  = GET_AREA_CM2( I, J, 1 ) * ( MOLWT / 6.023d23 ) * 1d-9
         
         ! Sum the emissions
         TOTAL = TOTAL + ( BBARRAY(I,J) * CONV )
      ENDDO
      ENDDO
     
      ! Define unit string
      SELECT CASE( NAME )
         CASE( 'NO' )
            UNIT = '[Tg N]'
         CASE( 'CO', 'CH2O', 'GLYX', 'MGLY', 'GLYC', 'HAC' )
            UNIT = '[Tg  ]'
         CASE DEFAULT
            UNIT = '[Tg C]'
      END SELECT

      ! Write totals
      ! increased length of NAME printed (hotp 7/31/09)
      WRITE( 6, 100 ) NAME, TOTAL, UNIT
 100  FORMAT( 'Sum Biomass ', a5, 1x, ': ', f9.3, 1x, a9  )

      ! Return to calling program
      END SUBROUTINE TOTAL_BIOMASS_TG

!------------------------------------------------------------------------------

      SUBROUTINE GC_READ_BIOMASS_BCOC( YEAR,       MONTH,
     &                                 BIOMASS_BC, BIOMASS_OC )
!
!******************************************************************************
!  Subroutine GC_READ_BIOMASS_BC_OC reads the GEOS-Chem default biomass 
!  emissions for black carbon and organic carbon.  (bmy, 9/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR       (INTEGER) : Current year
!  (2 ) MONTH      (INTEGER) : Current month
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) BIOMASS_BC (REAL*8 ) : Array for biomass BC emissions [atoms C/cm2/s]
!  (4 ) BIOMASS_OC (REAL*8 ) : Array for biomass OC emissions [atoms C/cm2/s]
! 
!  NOTES:
!  (1 ) Took the code that reads the emissions from disk from 
!        BIOMASS_CARB_GEOS in "carbon_mod.f". (bmy, 9/28/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  13 Aug 2013 - M. Sulprizio- Add POA1 for SOA + semivolatile POA (H. Pye) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCbb
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCbb
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LBBSEA,          LFUTURE
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR 
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTBCPO,         IDTOCPO
      ! SOAupdate (hotp 6/12/09)
      ! semivolpoa: replace OCPO with POA (hotp 2/17/09)
      USE TRACERID_MOD,         ONLY : IDTPOA1
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      INTEGER                       :: YEAR, MONTH
      REAL*8                        :: BIOMASS_BC(IIPAR,JJPAR)
      REAL*8                        :: BIOMASS_OC(IIPAR,JJPAR)

      ! Local variables
      INTEGER                       :: I,       J
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: CONV,    XTAU,   SEC_PER_MONTH
      CHARACTER(LEN=4)              :: CYEAR
      CHARACTER(LEN=255)            :: BC_FILE, OC_FILE

      ! Days per month (based on 1998)
      REAL*8 :: NDAYS(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0 /)

      !=================================================================
      ! GC_READ_BIOMASS_BC_OC begins here!
      !=================================================================

      ! Make sure BCPO, OCPO tracers are defined
      ! SOAupdate: add IDTPOA for svpoa (hotp 8/23/09)
      IF ( IDTBCPO == 0 .and. ( IDTOCPO == 0 .or. IDTPOA1 == 0 )) THEN
         BIOMASS_BC = 0d0
         BIOMASS_OC = 0d0
         RETURN
      ENDIF

      ! Test for leap year
      IF ( MONTH == 2 ) THEN
         IF( ITS_A_LEAPYEAR( YEAR ) ) THEN
            NDAYS(2) = 29d0
         ELSE
            NDAYS(2) = 28d0
         ENDIF
      ENDIF

      ! Number of seconds in this month
      SEC_PER_MONTH = 86400d0 * NDAYS(MONTH)

      ! Year string
      WRITE( CYEAR, '(i4)' ) YEAR

      !=================================================================
      ! Read BC/OC biomass burning [kg C/month] as tracers #34, 35
      !=================================================================

      ! Use seasonal or interannual emissions?
      IF ( LBBSEA ) THEN

         !------------------------------------
         ! Use seasonal biomass emissions
         !------------------------------------

         ! File name for seasonal BCPO biomass emissions
         BC_FILE = TRIM( DATA_DIR )                          //
     &             'biomass_200110/BCPO.bioburn.seasonal.'   //
     &             GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         ! File name for seasonal OCPO biomass emissions
         OC_FILE = TRIM( DATA_DIR )                         //
     &             'biomass_200110/OCPO.bioburn.seasonal.'  //
     &             GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( MONTH, 1, 1985 )      

      ELSE

         !------------------------------------
         ! Use interannual biomass emissions 
         ! for years between 1996 and 2002
         !------------------------------------

         ! File name for interannual BCPO biomass burning emissions`
         BC_FILE = TRIM( DATA_DIR )                           //
     &             'biomass_200110/BCPO.bioburn.interannual.' // 
     &             GET_NAME_EXT_2D() // '.'                   //
     &             GET_RES_EXT()     // '.'                   // CYEAR

         ! File name for interannual BCPO biomass burning emissions
         OC_FILE = TRIM( DATA_DIR )                           //
     &             'biomass_200110/OCPO.bioburn.interannual.' // 
     &             GET_NAME_EXT_2D() // '.'                   //
     &             GET_RES_EXT()     // '.'                   // CYEAR

         ! Use TAU0 value on the 1st of this month to index bpch file
         XTAU = GET_TAU0( MONTH, 1, YEAR )
  
      ENDIF

      !------------------
      ! Read BC biomass
      !------------------

      ! Echo info
      WRITE( 6, 100 ) TRIM( BC_FILE )
 100  FORMAT( '     - GC_READ_BIOMASS_BC_OC: Reading ', a )

      ! Read BC emission data [kg/mon]
      CALL READ_BPCH2( BC_FILE, 'BIOBSRCE', 34, 
     &                 XTAU,     IIPAR,     JJPAR,     
     &                 1,        ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), BIOMASS_BC )

      !------------------
      ! Read OC biomass
      !------------------

      ! Echo info
      WRITE( 6, 100 ) TRIM( OC_FILE )

      ! Read OC emission data [kg/mon]
      CALL READ_BPCH2( OC_FILE, 'BIOBSRCE', 35, 
     &                 XTAU,     IIPAR,     JJPAR,     
     &                 1,        ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), BIOMASS_OC )

      !=================================================================
      ! Convert from [kg C/mon] to [atoms C/cm2/s]
      !=================================================================
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Conversion factor for [1/cm2/s]
            CONV = 1d0 / ( GET_AREA_CM2( I, J, 1 ) * SEC_PER_MONTH )
          
            ! Convert [kg C/month] -> [atoms C/cm2/s]
            BIOMASS_BC(I,J) = BIOMASS_BC(I,J) * XNUMOL(IDTBCPO) * CONV

            ! SOAupdate: add POA fix (hotp 8/23/09)
            IF ( IDTOCPO > 0 ) THEN
              BIOMASS_OC(I,J) = BIOMASS_OC(I,J) * XNUMOL(IDTOCPO) * CONV
            ELSEIF ( IDTPOA1 > 0 ) THEN
              BIOMASS_OC(I,J) = BIOMASS_OC(I,J) * XNUMOL(IDTPOA1) * CONV
            ENDIF

            ! Scale to IPCC future scenario (if necessary)
            IF ( LFUTURE ) THEN

               ! Future scale BC biomass
               BIOMASS_BC(I,J) = BIOMASS_BC(I,J) * 
     &                           GET_FUTURE_SCALE_BCbb( I, J )

               ! Future scale OC biomass
               BIOMASS_OC(I,J) = BIOMASS_OC(I,J) * 
     &                           GET_FUTURE_SCALE_OCbb( I, J )
            ENDIF

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program 
      END SUBROUTINE GC_READ_BIOMASS_BCOC

!------------------------------------------------------------------------------
 
      SUBROUTINE GC_READ_BIOMASS_SO2( YEAR, MONTH, BIOMASS_SO2 )
!
!******************************************************************************
!  Subroutine GC_READ_BIOMASS_SO2 reads monthly mean biomass burning SO2
!  emissions.  (bmy, 9/28/06)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) YEAR        (INTEGER) : Current year
!  (2 ) MONTH       (INTEGER) : Current month
!
!  Arguments as Input:
!  ===========================================================================
!  (3 ) BIOMASS_SO2 (REAL*8 ) : Array for biomass SO2 [molec SO2/cm2/s]
!
!  NOTES:
!  (1 ) Took file reading code out of READ_BIOMASS_SO2 of "sulfate_mod.f"
!        and inserted here (bmy, 9/28/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2bb
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,          ONLY : LBBSEA,          LFUTURE
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTSO2
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D
      
      USE CMN_SIZE_MOD             ! Size parameters
                                    
      ! Arguments                   
      INTEGER, INTENT(IN)           :: YEAR, MONTH
      REAL*8,  INTENT(OUT)          :: BIOMASS_SO2(IIPAR,JJPAR)
                                    
      ! Local variables             
      INTEGER                       :: I, J, THISYEAR
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: CONV, XTAU, SEC_PER_MONTH
      CHARACTER(LEN=4  )            :: CYEAR
      CHARACTER(LEN=255)            :: FILENAME
                                    
      ! Days per month              
      REAL*8 :: NDAYS(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0 /)

      !=================================================================
      ! GC_READ_BIOMASS_SO2 begins here!
      !=================================================================

      ! Make sure BCPO, OCPO tracers are defined
      IF ( IDTSO2 == 0 ) THEN
         BIOMASS_SO2 = 0d0
         RETURN
      ENDIF

      ! Test for leap year
      IF ( MONTH == 2 ) THEN
         IF( ITS_A_LEAPYEAR( YEAR ) ) THEN
            NDAYS(2) = 29d0
         ELSE
            NDAYS(2) = 28d0
         ENDIF
      ENDIF

      ! Seconds in this month
      SEC_PER_MONTH = ( 86400d0 * NDAYS(MONTH) )

      ! Create a string for the 4-digit year
      WRITE( CYEAR, '(i4)' ) YEAR

      !=================================================================
      ! Read SO2 biomass emissions [kg SO2/month]
      !=================================================================

      ! Use seasonal or interannual emisisons?
      IF ( LBBSEA ) THEN

         !------------------------------------
         ! Use seasonal biomass emissions
         !------------------------------------

         ! File name for seasonal BB emissions
         FILENAME = TRIM( DATA_DIR )                         //
     &              'biomass_200110/SO2.bioburn.seasonal.'   //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( MONTH, 1, 1985 )      

      ELSE

         !------------------------------------
         ! Use interannual biomass emissions 
         ! for years between 1996 and 2002
         !------------------------------------

         ! File name for interannual biomass burning emissions
         FILENAME = TRIM( DATA_DIR )                          // 
     &              'biomass_200110/SO2.bioburn.interannual.' // 
     &              GET_NAME_EXT_2D() // '.'                  //
     &              GET_RES_EXT()     // '.'                  // CYEAR

         ! Use TAU0 value at start of this month to index punch file
         XTAU = GET_TAU0( MONTH, 1, YEAR )
    
      ENDIF

      !---------------------------------------
      ! Read biomass SO2 [kg SO2/month]
      !---------------------------------------

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GC_READ_BIOMASS_SO2: Reading ', a )

      ! Read SO2 emission data [kg/month]
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE',     26, 
     &                 XTAU,      IIPAR,         JJPAR,     
     &                 1,         ARRAY(:,:,1),  QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D ( ARRAY(:,:,1), BIOMASS_SO2 )

      !=================================================================
      ! Convert units [kg SO2/month] -> [molec SO2/cm2/s]
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor for [kg/month] -> [molec/cm2/s]
         CONV = XNUMOL(IDTSO2) 
     &        / ( GET_AREA_CM2( I, J, 1 ) * SEC_PER_MONTH )

         ! Convert [kg SO2/month] -> [molec SO2/cm2/s]
         BIOMASS_SO2(I,J) = BIOMASS_SO2(I,J) * CONV

         ! Scale to IPCC future scenario (if necessary)
         IF ( LFUTURE ) THEN
            BIOMASS_SO2(I,J) = BIOMASS_SO2(I,J) * 
     &                         GET_FUTURE_SCALE_SO2bb( I, J )
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GC_READ_BIOMASS_SO2     

!------------------------------------------------------------------------------

      SUBROUTINE GC_READ_BIOMASS_NH3( YEAR, MONTH, BIOMASS_NH3 ) 
!
!******************************************************************************
!  Subroutine GC_READ_BIOMASS_NH3 reads the monthly mean biomass NH3 
!  and biofuel emissions from disk and converts to [molec NH3/cm2/s]. 
!  (bmy, 9/28/06)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) YEAR        (INTEGER) : Current year
!  (2 ) MONTH       (INTEGER) : Current month
!
!  Arguments as Input:
!  ===========================================================================
!  (3 ) BIOMASS_NH3 (REAL*8 ) : Array for biomass NH3 [molec SO2/cm2/s]
!
!  NOTES:
!  (1 ) Took file reading code out of READ_BIOMASS_NH3 of "sulfate_mod.f"
!        and inserted here (bmy, 9/28/06)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3bb
      USE LOGICAL_MOD,          ONLY : LBBSEA,          LFUTURE
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE TIME_MOD,             ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNH3
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)           :: YEAR, MONTH
      REAL*8,  INTENT(OUT)          :: BIOMASS_NH3(IIPAR,JJPAR)

      ! Local variables
      INTEGER                       :: I, J
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: CONV, SEC_PER_MONTH, XTAU
      CHARACTER(LEN=4  )            :: CYEAR
      CHARACTER(LEN=255)            :: FILENAME

      ! Number of days in the month
      REAL*8 :: NDAYS(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0 /)

      !=================================================================
      ! READ_BIOMASS_NH3 begins here!
      !=================================================================

      ! Make sure BCPO, OCPO tracers are defined
      IF ( IDTNH3 == 0 ) THEN
         BIOMASS_NH3 = 0d0
         RETURN
      ENDIF

      ! Test for leap year
      IF ( MONTH == 2 ) THEN
         IF( ITS_A_LEAPYEAR( YEAR ) ) THEN
            NDAYS(2) = 29d0
         ELSE
            NDAYS(2) = 28d0
         ENDIF
      ENDIF

      ! Number of seconds in this month
      SEC_PER_MONTH = 86400d0 * NDAYS(MONTH)

      ! Create a string for the 4-digit year
      WRITE( CYEAR, '(i4)' ) YEAR

      !=================================================================
      ! Read biomass NH3 emissions [kg NH3/month]
      !=================================================================

      ! Use seasonal or interannual emisisons?
      IF ( LBBSEA ) THEN

         !------------------------------------
         ! Use seasonal biomass emissions
         !------------------------------------

         ! File name for seasonal BB emissions
         FILENAME = TRIM( DATA_DIR )                         // 
     &              'biomass_200110/NH3.bioburn.seasonal.'   //
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

         ! Get TAU0 value (use generic year 1985)
         XTAU = GET_TAU0( MONTH, 1, 1985 )      

      ELSE

         !------------------------------------
         ! Use interannual biomass emissions 
         ! for years between 1996 and 2002 
         !------------------------------------
    
         ! File name for interannual biomass burning emissions
         FILENAME = TRIM( DATA_DIR )                          // 
     &              'biomass_200110/NH3.bioburn.interannual.' // 
     &              GET_NAME_EXT_2D() // '.'                  //
     &              GET_RES_EXT()     // '.'                  // CYEAR

         ! Use TAU0 value on 1st day of this month to index bpch file
         XTAU = GET_TAU0( MONTH, 1, YEAR )
   
      ENDIF

      !---------------------------------------
      ! Read NH3 biomass [kg NH3/month]
      !---------------------------------------

      ! Echo filename
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BIOMASS_NH3: Reading ', a )

      ! Read NH3 emission data [kg/mon] as tracer 29
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE',    29, 
     &                 XTAU,      IIPAR,        JJPAR,      
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS_NH3 )

      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, CONV )

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Conversion factor for [kg/month] -> [molec/cm2/s]
         CONV = XNUMOL(IDTNH3) 
     &        / ( GET_AREA_CM2( I, J, 1 ) * SEC_PER_MONTH ) 


         ! Convert [kg NH3/month] -> [molec NH3/cm2/s]
         BIOMASS_NH3(I,J) = BIOMASS_NH3(I,J) * CONV

         ! Scale to IPCC future scenario (if necessary)
         IF ( LFUTURE ) THEN
            BIOMASS_NH3(I,J) = BIOMASS_NH3(I,J) * 
     &                         GET_FUTURE_SCALE_NH3bb( I, J )
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE GC_READ_BIOMASS_NH3

!------------------------------------------------------------------------------

       SUBROUTINE GC_READ_BIOMASS_CO2( YEAR, MONTH, BIOMASS_CO2 ) 
!
!******************************************************************************
!  Subroutine GC_READ_BIOMASS_CO2 reads in monthly values of CO for 
!  biomass burning from a binary punch file. (pns, bmy, 8/16/05, 9/28/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) YEAR  (INTEGER) : Current year (e.g. 1990)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Moved here from "co2_mod.f" (bmy, 9/28/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LBBSEA
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR
      USE TRACER_MOD,    ONLY : ITS_A_CO2_SIM
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: YEAR, MONTH
      REAL*8,  INTENT(OUT)   :: BIOMASS_CO2(IIPAR,JJPAR)

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                 :: EMFACTCO2CO, TAU, SEC_PER_MONTH
      CHARACTER(LEN=4)       :: SYEAR
      CHARACTER(LEN=255)     :: FILENAME

      ! Number of days in the month
      REAL*8 :: NDAYS(12) = (/ 31d0, 28d0, 31d0, 30d0, 31d0, 30d0, 
     &                         31d0, 31d0, 30d0, 31d0, 30d0, 31d0 /)

      !=================================================================
      ! READ_BIOMASS_CO2 begins here!
      !=================================================================

      ! Make sure it's a CO2 simulation
      IF ( .not. ITS_A_CO2_SIM() ) THEN
         BIOMASS_CO2 = 0d0
         RETURN
      ENDIF

      ! Test for leap year
      IF ( MONTH == 2 ) THEN
         IF( ITS_A_LEAPYEAR( YEAR ) ) THEN
            NDAYS(2) = 29d0
         ELSE
            NDAYS(2) = 28d0
         ENDIF
      ENDIF

      ! Seconds per month
      SEC_PER_MONTH = 86400d0 * NDAYS(MONTH)

      ! Currently calculate CO2 emissions as a function of CO emissions
      ! from biomass burning.   Set Emission factor (CO2 to CO)
      ! Calculation based on global totals of 
      !   5524.7  Tg dry matter (of which 45% is carbon)
      !    438.08 Tg CO (of which 187.75 Tg is carbon), and 
      !     32.6  Tg C of other species 
      !Refs : Staudt et al., Rose Yevich tables
      !Check with Rose Yevich for most recent estimates on emission factors 
      EMFACTCO2CO = 12.068d0

      ! Test for climatological or interannual emissions
      IF ( LBBSEA ) THEN

         !--------------------------------------
         ! Climatological biomass emissions
         !--------------------------------------
         
         ! TAU value for this month of "generic" year 1985
         TAU      = GET_TAU0( MONTH, 1, 1985 )

         ! Name of climatological biomass burning file
         FILENAME = TRIM( DATA_DIR )                   //
     &              'biomass_200110/bioburn.seasonal.' // 
     &              GET_NAME_EXT_2D() // '.' // GET_RES_EXT()

      ELSE

         !--------------------------------------
         ! Interannual biomass emissions
         !--------------------------------------

         ! Make a string for YEAR
         WRITE( SYEAR, '(i4)' ) YEAR 

         ! TAU value for the given month of this year
         TAU      = GET_TAU0( MONTH, 1, YEAR )

         ! Name of interannual biomass burning file
         FILENAME = TRIM( DATA_DIR )                      //
     &              'biomass_200110/bioburn.interannual.' // 
     &              GET_NAME_EXT_2D() // '.'              //
     &              GET_RES_EXT()     // '.'              // SYEAR 

      ENDIF

      !-----------------------------------------
      ! Read data from disk
      !-----------------------------------------

      ! Initialize ARRAY
      ARRAY = 0e0

      ! Read CO biomass emissions [molec CO/cm2/month]
      CALL READ_BPCH2( FILENAME, 'BIOBSRCE', 4, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), BIOMASS_CO2 )

      ! Convert from [molec CO/cm2/month] to [molec CO2/cm2/month]
      BIOMASS_CO2 = BIOMASS_CO2 * EMFACTCO2CO

      ! Print total CO2 biomass in Tg
      CALL TOTAL_BIOMASS_TG( BIOMASS_CO2, 44d-3, 'CO2' )

      ! Convert from [molec CO2/cm2/month] to [molec CO2/cm2/s]
      BIOMASS_CO2 = BIOMASS_CO2 / SEC_PER_MONTH

      ! Return to calling program
      END SUBROUTINE GC_READ_BIOMASS_CO2

!------------------------------------------------------------------------------

      SUBROUTINE ADJUST_TO_TOMSAI( BIOMASS_ANN, BIOMASS_SEA, BIOMASS )
!
!******************************************************************************
!  Subroutine ADJUST_TO_TOMSAI is a wrapper for subroutine TOMSAI.
!  (bmy, 10/12/00, 4/5/06)
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
!  (5 ) Now dimension arrays (I,J,N) instead of (N,I,J) (bmy, 4/5/06)
!******************************************************************************
!
      USE CMN_SIZE_MOD

      ! Arguments 
      REAL*8,  INTENT(INOUT) :: BIOMASS_ANN(IIPAR,JJPAR,NBIOMAX)
      REAL*8,  INTENT(INOUT) :: BIOMASS_SEA(IIPAR,JJPAR,NBIOMAX)
      REAL*8,  INTENT(INOUT) :: BIOMASS(IIPAR,JJPAR,NBIOMAX)

      ! Local variables
      INTEGER                :: I, J, N

      ! ADJUST_TO_TOMSAI begins here!
      WRITE( 6, '(a)' ) 'BIOBURN: Adjusting to TOMS AI data...'

      ! Loop over all tracers & boxes -- adjust to TOMS Aerosol index
      DO N = 1, NBIOMAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         CALL TOMSAI( I, J, BIOMASS_ANN(I,J,N), 
     &                      BIOMASS_SEA(I,J,N), BIOMASS(I,J,N) )
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
!  Inserted into F90 module "biomass_mod.f" (bmy, 9/25/00, 1/13/10)
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
!  (7 ) Added space in #ifdef block for 1 x 1.25 grid (bmy, 12/1/04)
!  (8 ) Modification for 1 x 1.25 grid (lok, bmy, 1/13/10)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE TIME_MOD,      ONLY : GET_MONTH, GET_TAU, GET_YEAR
      USE CMN_SIZE_MOD        ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J
      REAL*8,  INTENT(INOUT) :: ADJUSTED_VALUE
      REAL*8,  INTENT(IN)    :: VAL_ANN
      REAL*8,  INTENT(IN)    :: VAL_SEAS

      ! Local variables
      INTEGER                :: THEYEAR, IPICK, CTM_lat, CTM_lon
      INTEGER                :: II, JJ, KK, LL, N, AS, IMONTH, JYEAR
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

         ! Find a free file LUN
         N = findFreeLUN()

         ! Read TOMS Aerosol index data
         OPEN( N, FILE = TRIM( DATA_DIR ) // 'TOMSAI', STATUS='OLD' )
         DO JJ=1,NAIREGIONS
            READ( N, * ) readinTOMS
            
            II=0
            DO KK=1,NAIYEARS
               DO LL=1,12
                  II=II+1
                  TOMSAISCALE(JJ,KK,LL)=readinTOMS(II)
               ENDDO
            ENDDO

         ENDDO
         CLOSE(N)
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

#elif defined( GRID1x125 )
      CTM_LON     = I
      CTM_LAT     = J

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

      SUBROUTINE CLEANUP_GC_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_BIOMASS deallocates the BURNEMIS and
!  TOMSAISCALE arrays (bmy, 4/5/06)
!
!  NOTES:
!******************************************************************************
!      
      !=================================================================
      ! CLEANUP_GC_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( TOMSAISCALE ) ) DEALLOCATE( TOMSAISCALE )

      ! Return to calling program
      END SUBROUTINE CLEANUP_GC_BIOMASS

!------------------------------------------------------------------------------
      
      ! End of module
      END MODULE GC_BIOMASS_MOD
