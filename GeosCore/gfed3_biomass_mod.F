!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: gfed3_biomass_mod
!     
! !DESCRIPTION: Module GFED3\_BIOMASS\_MOD contains routines and variables 
!  used to incorporate GFED3 emissions into GEOS-Chem
!\\   
!\\   
! !INTERFACE: 
!
      MODULE GFED3_BIOMASS_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GFED3_COMPUTE_BIOMASS
      PUBLIC  :: CLEANUP_GFED3_BIOMASS
      PUBLIC  :: GFED3_IS_NEW
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_GFED3
      PRIVATE :: GFED3_AVAILABLE
      PRIVATE :: GFED3_SCALE_FUTURE
      PRIVATE :: GFED3_TOTAL_Tg
      PRIVATE :: INIT_GFED3_BIOMASS
      PRIVATE :: REARRANGE_BIOM
      PRIVATE :: READ_BPCH2_GFED3
!
! !REMARKS:
!  Monthly emissions of DM are read from disk,
!  multiplied by daily and 3hourly fractions (if necessary), and then
!  multiplied by the appropriate emission factors to produce biomass
!  burning emissions on the GFED3 0.5x0.5 degree grid  The emissions are 
!  then regridded to the current GEOS-Chem or GCAP grid (1x1, 2x25, or 4x5).
!                                                                             .
!  GFED3 biomass burning emissions are computed for the following gas-phase 
!  and aerosol-phase species:
!                                                                             .
!     (1 ) NOx  [  molec/cm2/s]     (13) BC   [atoms C/cm2/s]
!     (2 ) CO   [  molec/cm2/s]     (14) OC   [atoms C/cm2/s]                  
!     (3 ) ALK4 [atoms C/cm2/s]     (15) GLYX [  molec/cm2/s]    
!     (4 ) ACET [atoms C/cm2/s]     (16) MGLY [  molec/cm2/s]    
!     (5 ) MEK  [atoms C/cm2/s]     (17) BENZ [atoms C/cm2/s]  
!     (6 ) ALD2 [atoms C/cm2/s]     (18) TOLU [atoms C/cm2/s]     
!     (7 ) PRPE [atoms C/cm2/s]     (19) XYLE [atoms C/cm2/s]
!     (8 ) C3H8 [atoms C/cm2/s]     (20) C2H4 [atoms C/cm2/s]
!     (9 ) CH2O [  molec/cm2/s]     (21) C2H2 [atoms C/cm2/s]
!     (10) C2H6 [atoms C/cm2/s]     (22) GLYC [  molec/cm2/s]
!     (11) SO2  [  molec/cm2/s]     (23) HAC  [  molec/cm2/s]
!     (12) NH3  [  molec/cm2/s]     (24) CO2  [  molec/cm2/s]
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Original GFED3 database from Guido van der Werf 
!        http://www.falw.vu/~gwerf/GFED/GFED3/emissions/
!  (2 ) Giglio, L., Randerson, J. T., van der Werf, G. R., Kasibhatla, P. S.,
!       Collatz, G. J., Morton, D. C., and DeFries, R. S.: Assessing
!       variability and long-term trends in burned area by merging multiple 
!       satellite fire products, Biogeosciences, 7, 1171-1186, 
!       doi:10.5194/bg-7-1171-2010, 2010.
!  (3 ) van der Werf, G. R., Randerson, J. T., Giglio, L., Collatz, G. J.,
!       Mu, M., Kasibhatla, P. S., Morton, D. C., DeFries, R. S., Jin, Y., 
!       and van Leeuwen, T. T.: Global fire emissions and the contribution of 
!       deforestation, savanna, forest, agricultural, and peat fires 
!       (1997–2009), Atmos. Chem. Phys., 10, 11707-11735, 
!       doi:10.5194/acp-10-11707-2010, 2010. 
!
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  14 Feb 2012 - M. Payer      - Add modifications for CH4 (K. Wecht)
!  01 Mar 2012 - R. Yantosca   - Now reference new grid_mod.F90
!  06 Mar 2012 - P. Kasibhatla - Final version
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! N_EMFAC : Number of emission factors per species
      ! N_SPEC  : Number of species
      !=================================================================
      INTEGER,          PARAMETER   :: N_EMFAC = 6
      INTEGER,          PARAMETER   :: N_SPEC  = 25  !kjw

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! T3HR            : HH at start of the current 3-hr period.
      ! UPDATED         : flag to indicate if GFED3 emissions are updated
      ! UPDATED_MON     : flag to indicate if new month
      ! UPDATED_DAY     : flag to indicate if new day
      !                   - only set to true if daily emissions are used
      ! UPDATED_3HR     : flag to indicate if new 3-hour period
      !                   - only set to true if 3-hourly emissions are used
      ! SECONDS         : Number of seconds in the current month
      !
      ! Arrays
      !
      ! GFED3_SPEC_NAME : Array for GFED3 biomass species names
      ! GFED3_SPEC_MOLWT: Array for GFED3 biomass species molecular wts
      ! GFED3_SPEC_UNIT : Array for GFED3 biomass species emissions units
      ! GFED3_EMFAC     : Array for user-defined emission factors
      ! DM_GFED3_MON    : Array for monthly GFED3 DM burnt GFED3 grid
      ! DM_GFED3_DAY    : Array for daily GFED3 DM burnt on GFED3 grid
      ! FR_GFED3_3HR    : Array for 3hourly fractions on GFED3 grid
      ! HUMTROP_GFED3   : Array for GFED3 0.5x0.5 humid trop forest map
      ! BIOMASS_MODEL   : Array for GFED3 species emissions on model grid
      ! BIO_SAVE        : Index array to store IDBxxx values
      ! XEDGE_GFED3     : Array for lon edges of GFED3 grid
      ! YEDGE_GFED3     : Array for sin of at edges of GFED3 grid
      ! XEDGE_MODELG    : Array for lat edges of global grid at model res
      ! YEDGE_MODELG    : Array for sin of lat edges of global grid at model res
      !=================================================================

      ! Scalars
      LOGICAL                       :: UPDATED
      LOGICAL                       :: UPDATED_MON
      LOGICAL                       :: UPDATED_DAY
      LOGICAL                       :: UPDATED_3HR
      INTEGER                       :: T3HR
      REAL*8                        :: SECONDS
      INTEGER                       :: IIIPAR0
      INTEGER                       :: JJJPAR0

      ! Arrays
      CHARACTER(LEN=4), ALLOCATABLE :: GFED3_SPEC_NAME(:)
      REAL*8,           ALLOCATABLE :: GFED3_SPEC_MOLWT(:)
      CHARACTER(LEN=6), ALLOCATABLE :: GFED3_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: GFED3_EMFAC(:,:)
      REAL*8,           ALLOCATABLE :: DM_GFED3_MON(:,:,:)
      REAL*8,           ALLOCATABLE :: DM_GFED3_DAY(:,:,:)
      REAL*4,           ALLOCATABLE :: FR_GFED3_3HR(:,:,:)
      INTEGER,          ALLOCATABLE :: HUMTROP_GFED3(:,:)
      REAL*8,           ALLOCATABLE :: BIOMASS_MODEL(:,:,:)
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)
      REAL*8,           ALLOCATABLE :: XEDGE_GFED3(:)
      REAL*8,           ALLOCATABLE :: YEDGE_GFED3(:)
      REAL*8,           ALLOCATABLE :: XEDGE_MODELG(:)
      REAL*8,           ALLOCATABLE :: YEDGE_MODELG(:)

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
! !IROUTINE: gfed3_is_new
!
! !DESCRIPTION: Function GFED3\_IS\_NEW returns TRUE if GFED3 emissions 
!  have been updated.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GFED3_IS_NEW( ) RESULT( IS_UPDATED )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_UPDATED    ! =T if GFED3 is updated; =F otherwise
!
! !REMARKS:
!  Called from carbon_mod.f and sulfate_mod.f
!
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
      IS_UPDATED = UPDATED      

      END FUNCTION GFED3_IS_NEW
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_gfed3
!
! !DESCRIPTION: Subroutine CHECK\_GFED3 checks if we entered a new GFED period
!  since last emission timestep (ie, last call). The result depends
!  on the emissions time step, and the GFED time period used, as well
!  as MMDDHH at beginning of the GEOS-Chem run
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_GFED3( DOY, HH )
!
! !USES:
!
      USE LOGICAL_MOD, ONLY : LDAYBB3
      USE LOGICAL_MOD, ONLY : L3HRBB3
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,    ONLY : ITS_A_NEW_DAY
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: DOY   ! Day of year (0-365 or 0-366 leap years)
      INTEGER, INTENT(IN) :: HH    ! Hour of day (0-23)
!
! !REMARKS:
!  The routine computes the DOY (resp. HOUR) at start of the 1-day (resp.
!  3-hour) period we are in, if the 1-day (resp. 3-hr ) GFED3
!  option is on. Result is compared to previous value to indicate if new
!  data should be read.
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  06 Mar 2012 - P. Kasibhatla - final GFED3 version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: NEW_T3HR
      
      ! Reset to default
      UPDATED     = .FALSE.
      UPDATED_MON = .FALSE.
      UPDATED_DAY = .FALSE.
      UPDATED_3HR = .FALSE.

      ! Check if it is a new month
      IF ( ITS_A_NEW_MONTH() ) THEN
            UPDATED     = .TRUE.
            UPDATED_MON = .TRUE.
      ENDIF
      
      ! Check if it is a new day
      IF ( LDAYBB3 ) THEN
         IF ( ITS_A_NEW_DAY() ) THEN
            UPDATED     = .TRUE.
            UPDATED_DAY = .TRUE.
         ENDIF
      ENDIF

      ! Check if it is a new 3-hr period
      IF ( L3HRBB3 ) THEN

         NEW_T3HR = INT( HH / 3 ) * 3

         IF ( NEW_T3HR .NE. T3HR ) THEN
            UPDATED     = .TRUE.
            UPDATED_3HR = .TRUE.
            T3HR        = NEW_T3HR
         ENDIF     
    
      ENDIF
      
      END SUBROUTINE CHECK_GFED3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gfed3_available
!
! !DESCRIPTION: Function GFED3\_AVAILABLE checks an input YYYY year and MM
!  month against the available data dates.  If the requested YYYY and MM
!  lie outside of the valid range of dates, then GFED3\_AVAILABLE will return 
!  the last valid YYYY and MM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GFED3_AVAILABLE( YYYY, YMIN, YMAX, MM, MMIN, MMAX )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)              :: YMIN, YMAX   ! Min & max years 
      INTEGER, INTENT(IN),    OPTIONAL :: MMIN, MMAX   ! Min & max months
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT)           :: YYYY         ! Year of GFED3 data
      INTEGER, INTENT(INOUT), OPTIONAL :: MM           ! Month of GFED3 data
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Check year
      IF ( YYYY > YMAX .OR. YYYY < YMIN ) THEN
         
         YYYY = MAX( YMIN, MIN( YYYY, YMAX) )
         
         WRITE( 6, 100 ) YMAX, YMIN, YYYY
 100     FORMAT( 'YEAR > ', i4, ' or YEAR < ', i4, 
     $        '. Using GFED3 biomass for ', i4)
      ENDIF
      

      ! Check month
      IF ( PRESENT( MM ) ) THEN 
         IF ( MM > MMAX .OR. MM < MMIN ) THEN

            MM = MAX( MMIN, MIN( MM, MMAX) )
            
            WRITE( 6, 200 ) MMIN, MMAX, MM
 200        FORMAT( ' ** WARNING ** : MONTH is not within ', i2,'-',
     $              i2, '. Using GFED3 biomass for month #', i2)
         ENDIF
      ENDIF

      END SUBROUTINE GFED3_AVAILABLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gfed3_compute_biomass
!
! !DESCRIPTION: Subroutine GFED3\_COMPUTE\_BIOMASS computes the monthly 
!  GFED3 biomass burning emissions for a given year and month. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GFED3_COMPUTE_BIOMASS( THIS_YYYY, THIS_MM, BIOM_OUT )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_TAU0
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE JULDAY_MOD,      ONLY : JULDAY
      USE JULDAY_MOD,      ONLY : CALDATE
      USE LOGICAL_MOD,     ONLY : LFUTURE
      USE LOGICAL_MOD,     ONLY : LDAYBB3
      USE LOGICAL_MOD,     ONLY : L3HRBB3
      USE LOGICAL_MOD,     ONLY : LGFED3BB
      USE TIME_MOD,        ONLY : EXPAND_DATE
      USE TIME_MOD,        ONLY : TIMESTAMP_STRING
      USE TIME_MOD,        ONLY : GET_DAY
      USE TIME_MOD,        ONLY : GET_HOUR
      USE TIME_MOD,        ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,        ONLY : ITS_A_LEAPYEAR
      USE GLOBAL_GRID_MOD, ONLY : GET_XEDGE_G
      USE GLOBAL_GRID_MOD, ONLY : GET_YEDGE_G
      USE GRID_MOD,        ONLY : GET_XEDGE
      USE GRID_MOD,        ONLY : GET_YEDGE
      USE GRID_MOD,        ONLY : GET_XOFFSET
      USE GRID_MOD,        ONLY : GET_YOFFSET
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE REGRID_A2A_MOD,  ONLY : MAP_A2A

      USE CMN_SIZE_MOD                                       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THIS_YYYY                      ! Current year   
      INTEGER, INTENT(IN)  :: THIS_MM                        ! Current month
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: BIOM_OUT(IIPAR,JJPAR,NBIOMAX)  ! BB emissions
                                                             ! [molec/cm2/s]
!
! !REMARKS:
!  This routine has to be called on EVERY emissions-timestep if you use one
!  of the GFED3 options.
! 
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  03 May 2012 - M. Payer      - Now use 3D indices for grid_mod.F90 routines
!                                GET_XEDGE and GET_YEDGE. Nested grids still 
!                                use 1D indices for routines GET_XEDGE_G and
!                                GET_YEDGE_G found in global_grid_mod.90. This
!                                will maintain grid-independent functionality.
!  06 Jun 2013 - M. Payer      - Set TAU0 to first of month for 3-hr emissions
!                                to match time stamp of data file (C. Keller)
!  02 Jul 2013 - M. Payer      - Update GFED3 emissions to 2011 (P. Kasibhatla)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I,    J,  N,   NF , IT3
      INTEGER             :: AS
      INTEGER             :: II,   JJ
      INTEGER             :: I0,   J0
      INTEGER             :: YYYY, MM, MM1, YYYY1
      INTEGER             :: YYYYMMDD, HHMMSS
      REAL*8              :: GFED3_EMFACX
      REAL*4              :: ARRAY_GFED3(IGFED3, JGFED3, 1)
      REAL*4              :: FR_GFED3_DAY(IGFED3, JGFED3)
      REAL*8              :: DM_GFED3(IGFED3, JGFED3, N_EMFAC)
      REAL*8              :: BIOMASS_GFED3(IGFED3, JGFED3, N_SPEC)
      REAL*8              :: TAU0, TAU1
      REAL*4              :: TMP
      CHARACTER(LEN=255)  :: FILENAME1
      CHARACTER(LEN=255)  :: FILENAME2
      CHARACTER(LEN=255)  :: FILENAME3
      CHARACTER(LEN=255)  :: FILENAME4
      CHARACTER(LEN=255)  :: FILENAME5
      CHARACTER(LEN=255)  :: FILENAME6
      CHARACTER(LEN=255)  :: FILENAME7
      CHARACTER(LEN=255)  :: FILENAME8
      CHARACTER(LEN=16 )  :: TIME_STR
      INTEGER             :: DD, HH, DOY
      INTEGER             :: IT3HR, IT3H
      REAL*4              :: DEG2RAD
      REAL*8, ALLOCATABLE :: bq1(:,:)
      REAL*8, ALLOCATABLE :: bq2(:,:)

      !=================================================================
      ! GFED3_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_GFED3_BIOMASS

         DEG2RAD = (4. * ATAN(1.) ) /180.

         ! Define GFED3 grid box lat and lon edges
         XEDGE_GFED3( 1 ) = -180.d0
         DO I = 2,IGFED3+1
            XEDGE_GFED3( I ) = XEDGE_GFED3( I-1 ) + 5.d-1
         END DO

         YEDGE_GFED3( 1 ) = -90.d0
         DO J = 2, JGFED3+1
            YEDGE_GFED3( J ) = YEDGE_GFED3( J-1 ) + 5.d-1
         END DO

         DO J = 1,JGFED3+1
            YEDGE_GFED3( J ) = SIN( YEDGE_GFED3( J ) * DEG2RAD)
         END DO

         ! Define global grid box lat and lon edges at model resolution
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE_G ( I )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE_G ( J )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#else

         DO I = 1,IIIPAR0+1
            XEDGE_MODELG( I ) = GET_XEDGE( I, 1, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = GET_YEDGE( 1, J, 1 )
         END DO

         DO J = 1,JJJPAR0+1
            YEDGE_MODELG( J ) = SIN( YEDGE_MODELG( J ) * DEG2RAD)
         END DO

#endif

         FIRST = .FALSE.

      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM
      DD   = GET_DAY()
      HH   = GET_HOUR()
      DOY  = GET_DAY_OF_YEAR()
      
      ! Check if we need to update GFED3 
      CALL CHECK_GFED3( DOY, HH )
      
      ! If no updating is needed, module variable BIOMASS_MODEL
      ! from last update can be used
      IF ( .not. UPDATED ) THEN
         CALL REARRANGE_BIOM(BIOMASS_MODEL,BIOM_OUT)
         RETURN
      ENDIF

      ! If within same month, no nead to reread emission file.
      ! But go to statement 999 to check if daily fractiions
      ! need to be read.
      IF ( .not. UPDATED_MON ) THEN
            go to 999
      ENDIF
      
      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     &  'G F E D 3   B I O M A S S   B U R N I N G   E M I S S I O N S'

      
      !=================================================================
      ! Check GFED3 availability & get YYYYMMDD of data to read.
      !=================================================================
         
      ! Availability of MONTHLY data
      !-------------------------------
         
      CALL GFED3_AVAILABLE( YYYY, 1997, 2011 )

      ! Create YYYYMMDD integer values
      YYYYMMDD = YYYY*10000 + MM*100 + 01

      !=================================================================
      ! Filename, TAU0 and number of seconds
      !=================================================================
      
      ! for monthly data
      !-------------------------------
      ! TAU value at start of YYYY/MM
      TAU0     = GET_TAU0( MM, 1, YYYY )

      ! Get YYYY/MM value for next month
      MM1      = MM + 1
      YYYY1    = YYYY

      ! Increment year if necessary
      IF ( MM1 == 13 ) THEN
         MM1   = 1
         YYYY1 = YYYY + 1
      ENDIF

      ! TAU value at start of next month
      TAU1     = GET_TAU0( MM1, 1, YYYY1 )

      ! Number of seconds in this month 
      ! (NOTE: its value will be saved until the next month)
      SECONDS  = ( TAU1 - TAU0 ) * 3600d0

      ! File name with GFED3 DM emissions
      FILENAME1 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_AGW_YYYYMM'
      FILENAME2 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_DEF_YYYYMM'
      FILENAME3 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_FOR_YYYYMM'
      FILENAME4 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_PET_YYYYMM'
      FILENAME5 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_SAV_YYYYMM'
      FILENAME6 = TRIM( DATA_DIR_NATIVE ) //
     &            'GFED3_201212/YYYY/Monthly/GFED3_DM_WDL_YYYYMM'

      !=================================================================
      ! Read GFED3 DM burnt [g/m2/month]
      !=================================================================
      
      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME1, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME2, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME3, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME4, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME5, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME6, YYYYMMDD, 000000 )

      ! Read GFED3 DM emissions [g DM/m2/month] in the following order
      ! AGW, DEF, FOR, PET, SAV, WDL
      CALL READ_BPCH2_GFED3( FILENAME1, 'GFED3-BB',   91, 
     &                       TAU0,      IGFED3,       JGFED3,     
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. )
 
      DM_GFED3_MON(:,:,1) = ARRAY_GFED3(:,:,1)

      CALL READ_BPCH2_GFED3( FILENAME2, 'GFED3-BB',   92, 
     &                       TAU0,      IGFED3,       JGFED3,   
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. ) 

      DM_GFED3_MON(:,:,2) = ARRAY_GFED3(:,:,1)

      CALL READ_BPCH2_GFED3( FILENAME3, 'GFED3-BB',   93, 
     &                       TAU0,      IGFED3,       JGFED3,   
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. ) 

      DM_GFED3_MON(:,:,3) = ARRAY_GFED3(:,:,1)

      CALL READ_BPCH2_GFED3( FILENAME4, 'GFED3-BB',   94, 
     &                       TAU0,      IGFED3,       JGFED3,   
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. ) 

      DM_GFED3_MON(:,:,4) = ARRAY_GFED3(:,:,1)

      CALL READ_BPCH2_GFED3( FILENAME5, 'GFED3-BB',   95, 
     &                       TAU0,      IGFED3,       JGFED3,   
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. ) 

      DM_GFED3_MON(:,:,5) = ARRAY_GFED3(:,:,1)

      CALL READ_BPCH2_GFED3( FILENAME6, 'GFED3-BB',   96, 
     &                       TAU0,      IGFED3,       JGFED3,   
     &                       1,         ARRAY_GFED3,  QUIET=.TRUE. ) 

      DM_GFED3_MON(:,:,6) = ARRAY_GFED3(:,:,1)

      !=================================================================
      ! Convert [g DM/m2/month] to [kg DM/cm2/month]
      !
      ! Unit Conversions:
      ! (1) g    to kg    --> Divide by 1000
      ! (2) 1/m2 to 1/cm2 --> Divide by 10000
      !=================================================================

      ! Loop over GFED3 GRID
      DO J = 1, JGFED3
      DO I = 1, IGFED3
      DO NF = 1, N_EMFAC

         ! Set negatives to zero 
         DM_GFED3_MON(I,J,NF) = MAX( DM_GFED3_MON(I,J,NF), 0e0 )

         ! Convert [g DM/m2/month] to [kg DM/cm2/month]
         DM_GFED3_MON(I,J,NF) = DM_GFED3_MON(I,J,NF) * 1d-3 * 1d-4 

      ENDDO
      ENDDO
      ENDDO

      ! If 3-hourly emissions are used, read 3-hourly fractions
      ! at the same time that monthly emissions are read because
      ! these fractions are constant throughout the month
      ! - note that these should be applied after daily fractions
      !   are applied

      IF ( L3HRBB3 ) THEN
         DO IT3 = 1,8

            IT3H = (IT3-1)*3
            HHMMSS = IT3H*10000
            TAU0 = GET_TAU0( MM, 1, YYYY, IT3H )

            FILENAME7 = TRIM( DATA_DIR_NATIVE ) //
     &         'GFED3_201212/YYYY/3hourly/GFED3_FR_3HR_YYYYMMDDhh'

            ! Replace YYYY/MM/HH in the file name
            CALL EXPAND_DATE( FILENAME7, YYYYMMDD, HHMMSS )

            CALL READ_BPCH2_GFED3( FILENAME7, 'GFED3-BB',  89,
     &                             TAU0,      IGFED3,      JGFED3,
     &                             1,         ARRAY_GFED3, QUIET=.TRUE.)

            FR_GFED3_3HR(:,:,IT3) = ARRAY_GFED3(:,:,1)

         END DO

      ENDIF

999   CONTINUE

      !=================================================================
      ! At this point in the code, the following cases are possible:
      !
      ! UPDATED_MON=T, UPDATED_DAY=T, UPDATED_3HR=T
      ! UPDATED_MON=T, UPDATED_DAY=T, UODATED_3HR=F 
      ! UPDATED_MON=F, UPDATED_DAY=T, UPDATED_3HR=T
      ! UPDATED_MON=F, UPDATED_DAY=T, UPDATED_3HR=F
      ! UPDATED_MON=F, UPDATED_DAY=F, UPDATED_3HR=T
      !
      ! Note that the combination 
      ! UPDATED_MON=F, UPDATED_DAY=F, UPDATED_3HR=F
      ! is not possible at this point in code because of the
      ! of the RETURN statement when UPDATED=F near the
      ! start of this subroutine
      ! 
      ! Also note that the combinations
      ! UPDATED_MON=T, UPDATED_DAY=F, UPDATED_3HR=T
      ! UPDATED_MON=T, UPDATED_DAY=F, UODATED_3HR=F 
      ! are not possible because UPDATED_DAY=T
      ! when UPDATED_MON=T 
      !
      ! In the following code, the module variables
      ! DM_GFED3_MON and DM_GFED3_DAY contain the 
      ! latest updated monthly and daily (if applicable)
      ! emissions, respectively. Note that
      ! by making them module variables, it is ensured
      ! that the correct values are always available
      ! during repeated calls to this code.
      !
      ! The local variable DM_GFED3 is updated appropriately
      ! based on the 5 allowable cases describe above and passed
      ! to the DO_REGRID_G2G_1x1 subroutine.
      !
      !================================================================= 

      ! Read daily fractions
      IF ( UPDATED_DAY) THEN
         TAU0 = GET_TAU0( MM, DD, YYYY )

         ! Create YYYYMMDD integer value
         YYYYMMDD = YYYY*10000 + MM*100 + DD

         FILENAME8 = TRIM( DATA_DIR_NATIVE ) //
     &               'GFED3_201212/YYYY/Daily/GFED3_FR_DAY_YYYYMMDD'

         ! Replace YYYY/MM in the file name
         CALL EXPAND_DATE( FILENAME8, YYYYMMDD, 000000 )

         CALL READ_BPCH2_GFED3( FILENAME8, 'GFED3-BB',   88,
     &                          TAU0,      IGFED3,       JGFED3,
     &                          1,         ARRAY_GFED3,  QUIET=.TRUE. )

         FR_GFED3_DAY(:,:) = ARRAY_GFED3(:,:,1)

      ENDIF

! Convert DM burnt from kg/cm2/month to kg/cm2/day or
! kg/cm2/3hour if needed, and grid to model grid
      DO NF = 1, N_EMFAC

         DO J = 1, JGFED3
         DO I = 1, IGFED3
            DM_GFED3(I,J,NF)=DM_GFED3_MON(I,J,NF)

            IF ( UPDATED_DAY ) THEN
               DM_GFED3_DAY(I,J,NF)=DM_GFED3_MON(I,J,NF)
     &                                 *FR_GFED3_DAY(I,J)
               DM_GFED3(I,J,NF)=DM_GFED3_DAY(I,J,NF)
               SECONDS = 24 * 3600d0
            END IF

            IF ( UPDATED_3HR ) THEN
               IT3HR=T3HR/3+1
               DM_GFED3(I,J,NF)=DM_GFED3_DAY(I,J,NF)
     &                                *FR_GFED3_3HR(I,J,IT3HR)
               SECONDS =  3 * 3600d0
            END IF

         ENDDO
         ENDDO
 
      END DO

      !=================================================================
      ! Calculate biomass species emissions on GFED3 grid
      ! and regrid to model grid
      !
      ! Emission factors convert from [kg DM/cm2/timeperiod] to either
      ! [molec/cm2/timeperiod] or [atoms C/cm2/timeperiod]
      !
      ! Units:
      !  [  molec/cm2/month] : NOx,  CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/month] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8,
      !                        C2H6, BC,   OC
      !=================================================================

      ! Loop over biomass species 
      DO N = 1, N_SPEC

         DO J = 1, JGFED3
         DO I = 1, IGFED3
            BIOMASS_GFED3(I,J,N) = 0.0
            DO NF = 1, N_EMFAC
               GFED3_EMFACX=GFED3_EMFAC(N,NF)

               ! Use woodland emission factors for 'deforestation' outside
               ! humid tropical forest 
               IF(NF.EQ.2.AND.HUMTROP_GFED3(I,J).EQ.0)
     &            GFED3_EMFACX=GFED3_EMFAC(N,6)
                  BIOMASS_GFED3(I,J,N) =  BIOMASS_GFED3(I,J,N) +
     &                                    DM_GFED3(I,J,NF) *
     &                                    GFED3_EMFACX

            ENDDO 
         ENDDO
         ENDDO

         ! Regrid emissions from GFED3 grid to model grid
         ALLOCATE( bq1( IGFED3, JGFED3 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'bq1' )
         ALLOCATE( bq2( IIIPAR0, JJJPAR0 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'bq2' )

         bq1(:,:)=BIOMASS_GFED3(:,:,N)

         CALL MAP_A2A( IGFED3,      JGFED3,       XEDGE_GFED3,
     &                 YEDGE_GFED3, bq1,          IIIPAR0,
     &                 JJJPAR0,     XEDGE_MODELG, YEDGE_MODELG,
     &                 bq2,         0,            0 )

         I0    = GET_XOFFSET( GLOBAL=.TRUE. )
         J0    = GET_YOFFSET( GLOBAL=.TRUE. )

         DO JJ=1,JJPAR
         DO II=1,IIPAR
            BIOMASS_MODEL(II,JJ,N)=bq2(II+I0,JJ+J0)
         END DO
         END DO

         DEALLOCATE( bq1 )
         DEALLOCATE( bq2 )

      END DO

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL GFED3_SCALE_FUTURE( BIOMASS_MODEL )
      ENDIF

      ! Print totals in Tg/time period
      IF ( UPDATED_3HR ) THEN
         WRITE( 6, 412 ) YYYY, MM, DD, T3HR
 412     FORMAT( 'GFED3 3hourly emissions for year, month, day, 3hr: ',
     &            i4, '/', 3i2.2, / )
         go to 998
      ENDIF
      IF ( UPDATED_DAY ) THEN
         WRITE( 6, 411 ) YYYY, MM, DD
 411     FORMAT( 'GFED3 daily emissions for year, month, day: ',
     &            i4, '/', 2i2.2, / )
         go to 998
      ENDIF
      WRITE( 6, 410 ) YYYY, MM
 410     FORMAT( 'GFED3 monthly emissions for year, month: ',
     &            i4, '/', i2.2, / )
 998  CONTINUE
      CALL GFED3_TOTAL_Tg

      ! Convert from [molec/cm2/month], [molec/cm2/day] or
      ! [molec/cm2/3hr] to [molec/cm2/s]
      BIOMASS_MODEL = BIOMASS_MODEL / SECONDS

      ! Rearrange the species to the same order as in the IDBxxx (fp, 6/09)
      ! BIOMASS_MODEL is indexed as GFED3
      ! BIOM_OUT      is indexed as IDBs
      CALL REARRANGE_BIOM( BIOMASS_MODEL, BIOM_OUT )

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE GFED3_COMPUTE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gfed3_scale_future
!
! !DESCRIPTION: Subroutine GFED3\_SCALE\_FUTURE applies the IPCC future 
!  emissions scale factors to the GFED3 biomass burning emisisons in order 
!  to compute the future emissions of biomass burning for NOx, CO, and VOC's.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GFED3_SCALE_FUTURE( BB )
!
! !USES:
!
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_BCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NH3bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_OCbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_SO2bb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbb
      USE TRACER_MOD,             ONLY : ITS_A_CO2_SIM 
      USE TRACER_MOD,             ONLY : ITS_A_CH4_SIM
      USE TRACERID_MOD,           ONLY : IDBNO,   IDBCO,   IDBSO2 
      USE TRACERID_MOD,           ONLY : IDBNH3,  IDBBC,   IDBOC 

      USE CMN_SIZE_MOD                 ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      !  Array w/ biomass burning emisisons [molec/cm2]
      REAL*8, INTENT(INOUT) :: BB(IIPAR,JJPAR,N_SPEC)
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  14 Mar 2013 - M. Payer      - Replace NOx emissions with NO emissions as part
!                                of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_CO2
      LOGICAL :: ITS_CH4
      INTEGER :: I, J, N
      
      !=================================================================
      ! GFED3_SCALE_FUTURE begins here!
      !=================================================================

      ! Test if it's a CO2 simulation outside of the loop
      ITS_CO2 = ITS_A_CO2_SIM()
      ITS_CH4 = ITS_A_CH4_SIM()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )

      ! Loop over species and grid boxes
      DO N = 1, N_SPEC
      DO J = 1, JJPAR
      DO I = 1, IIPAR 

         ! Scale each species to IPCC future scenario
         IF ( BIO_SAVE(N) == IDBNO ) THEN

            ! Future biomass NOx [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NOxbb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBCO ) THEN

            ! Future biomass CO [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_CObb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBSO2 ) THEN

            ! Future biomass SO2 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_SO2bb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBNH3 ) THEN

            ! Future biomass NH3 [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NH3bb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBBC ) THEN

            ! Future biomass BC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_BCbb( I, J )

         ELSE IF ( BIO_SAVE(N) == IDBOC ) THEN

            ! Future biomass OC [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_OCbb( I, J )

         ! Don't scale future emissions if CO2 or CH4
         ELSE IF ( ITS_CO2 .OR. ITS_CH4 ) THEN

            ! Nothing

         ELSE

            ! Future biomass Hydrocarbons [atoms C/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_VOCbb( I, J )

         ENDIF
         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE GFED3_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gfed3_total_Tg
!
! !DESCRIPTION: Subroutine GFED3\_TOTAL\_Tg prints the amount of biomass 
!  burning emissions that are emitted each month/day/3-hr in Tg or Tg C. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GFED3_TOTAL_Tg
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD      ! Size parameters
!
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I,    J,     N
      REAL*8           :: CONV, MOLWT, TOTAL
      CHARACTER(LEN=4) :: NAME
      CHARACTER(LEN=6) :: UNIT

      !=================================================================
      ! GFED3_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         ! Initialize
         NAME  = GFED3_SPEC_NAME(N)
         MOLWT = GFED3_SPEC_MOLWT(N)
         UNIT  = GFED3_SPEC_UNIT(N)
         TOTAL = 0d0

         ! Loop over surface boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert to [Tg/gfed-period] (or [Tg C/gfed-period] for HC's)
            CONV  = GET_AREA_CM2( I, J, 1 ) 
     &            * ( MOLWT / 6.023d23 ) * 1d-9

            TOTAL = TOTAL + ( BIOMASS_MODEL(I,J,N) * CONV )
         ENDDO
         ENDDO
     
         ! Write totals
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', e12.5, 1x, a6 )
      ENDDO

      END SUBROUTINE GFED3_TOTAL_Tg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gfed3_biomass
!
! !DESCRIPTION: Subroutine INIT\_GFED3\_BIOMASS allocates all module arrays.  
!  It also reads the emission factors at the start of a GEOS-Chem
!  simulation. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GFED3_BIOMASS
!
! !USES:
!
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE FILE_MOD,        ONLY : IOERROR
      USE LOGICAL_MOD,     ONLY : LDICARB
      USE LOGICAL_MOD,     ONLY : LDAYBB3
      USE LOGICAL_MOD,     ONLY : L3HRBB3
      USE TRACERID_MOD,    ONLY : IDBNO,  IDBCO,   IDBALK4
      USE TRACERID_MOD,    ONLY : IDBACET, IDBMEK,  IDBALD2
      USE TRACERID_MOD,    ONLY : IDBPRPE, IDBC3H8, IDBCH2O
      USE TRACERID_MOD,    ONLY : IDBC2H6, IDBBC,   IDBOC
      USE TRACERID_MOD,    ONLY : IDBSO2,  IDBNH3,  IDBCO2
      USE TRACERID_MOD,    ONLY : IDBGLYX, IDBMGLY, IDBBENZ
      USE TRACERID_MOD,    ONLY : IDBTOLU, IDBXYLE, IDBC2H4
      USE TRACERID_MOD,    ONLY : IDBC2H2, IDBGLYC, IDBHAC
      USE TRACERID_MOD,    ONLY : IDBCH4   !kjw
      USE GLOBAL_GRID_MOD, ONLY : GET_IIIPAR
      USE GLOBAL_GRID_MOD, ONLY : GET_JJJPAR
      USE CMN_SIZE_MOD          ! Size parameters
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IU_FILE
      INTEGER            :: AS, IOS, M, N, NDUM
      REAL*4             :: ARRAY_LANDMAP(IGFED3,JGFED3,1)
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! INIT_GFED3_BIOMASS begins here!
      !=================================================================

      ! Get global longitude/latitude extent [# of boxes]
#if   defined( GRID05x0666 ) || defined( GRID025x03125 )

      ! Nested grids utilize global longitude and latidude extent
      ! parameters IIIPAR and JJJPAR (from global_grid_mod)
      IIIPAR0 = GET_IIIPAR()
      JJJPAR0 = GET_JJJPAR()

#else

      ! Global grids utilize window longitude and latidude extent
      ! parameters IIPAR and JJPAR (from CMN_SIZE_mod)
      IIIPAR0 = IIPAR
      JJJPAR0 = JJPAR

#endif

      ! Allocate array to hold GFED3 grid box lon edges
      ALLOCATE( XEDGE_GFED3( IGFED3+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XEDGE_GFED3' )
      XEDGE_GFED3 = 0.d0

      ! Allocate array to hold GFED3 grid box lat edges
      ALLOCATE( YEDGE_GFED3( JGFED3+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_GFED3' )
      YEDGE_GFED3 = 0.d0

      ! Allocate array to hold GEOS-Chem grid box lon edges
      ALLOCATE( XEDGE_MODELG( IIIPAR0+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XEDGE_MODELG' )
      XEDGE_MODELG = 0.d0

      ! Allocate array to hold GEOS-Chem grid box lat edges
      ALLOCATE( YEDGE_MODELG( JJJPAR0+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_MODELG' )
      YEDGE_MODELG = 0.d0

      ! Allocate array to hold GFED3 species emissions on model grid
      ALLOCATE( BIOMASS_MODEL( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_MODEL' )
      BIOMASS_MODEL = 0d0

      ! Allocate array to hold monthly GFED3 DM burnt GFED3 grid
      ALLOCATE( DM_GFED3_MON( IGFED3, JGFED3, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DM_GFED3_MON' )
      DM_GFED3_MON = 0d0

      ! Allocate array to hold daily GFED3 DM burnt GFED3 grid
      IF ( LDAYBB3 ) THEN
         ALLOCATE( DM_GFED3_DAY( IGFED3, JGFED3, N_SPEC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'DM_GFED3_DAY' )
         DM_GFED3_DAY = 0d0
      ENDIF

      ! Allocate array to hold 3hourly fractions on GFED3 grid
      IF ( L3HRBB3 ) THEN
         ALLOCATE( FR_GFED3_3HR( IGFED3,JGFED3, 8 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'FR_GFED3_3HR' )
         FR_GFED3_3HR = 0d0
      ENDIF

      ! Allocate array for emission factors
      ALLOCATE( GFED3_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_EMFAC' )
      GFED3_EMFAC = 0d0
      
      ! Allocate array for species molecular weight
      ALLOCATE( GFED3_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_SPEC_MOLWT' )
      GFED3_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( GFED3_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_SPEC_NAME' )
      GFED3_SPEC_NAME = ''

      ! Allocate array for GFED3 biomass buning species mass units
      ALLOCATE( GFED3_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_SPEC_UNIT' )
      GFED3_SPEC_UNIT = ''

      ! Allocate array for vegetation map
      ALLOCATE( HUMTROP_GFED3( IGFED3, JGFED3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HUMTROP_GFED3' )

      !IDBs are now the same as the ones in TRACERID AND BIOMASS_MOD
      !BIOSAVE INDEX IS THE LOCATION OF THE EMISSION IN THE GFED FILE
      !(fp)
      ALLOCATE( BIO_SAVE( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIO_SAVE' )
      BIO_SAVE = 0

      ! Set default values for module variables
      T3HR    = -1

      !=================================================================
      ! Read emission factors (which convert from kg DM to 
      ! either [molec species] or [atoms C]) from bpch file
      !=================================================================
     
      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE) // 
     &           'GFED3_201212/GFED3_emission_factors.txt'

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed3:1' )

      ! Skip header lines
      DO N = 1, 9 
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed3:2' )
      ENDDO

      ! Read emission factors for each species
      DO N = 1, N_SPEC
         READ( IU_FILE, 100, IOSTAT=IOS ) 
     &       NDUM, GFED3_SPEC_NAME(N), ( GFED3_EMFAC(N,M), M=1,N_EMFAC )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed3:3' )
      WRITE(6,100)NDUM,GFED3_SPEC_NAME(N),(GFED3_EMFAC(N,M),M=1,N_EMFAC)
      ENDDO
      
      ! FORMAT string
 100  FORMAT( 1x, i2, 1x, a4, 6(3x,es14.6) )

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Read GFED humid tropical forest map from bpch file
      ! This is used to assign emission factors for 'deforestation'
      ! 'Deforestation' occur outside of humid tropical forest
      ! is assigned a 'woodlands' emission factor'
      !
      ! Values:  1 = humid tropical forest
      !          0 = other
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE ) //
     &           'GFED3_201212/GFED3_humtropmap'

      ! Read GFED3 veg map
      CALL READ_BPCH2_GFED3( FILENAME, 'LANDMAP',     1,
     &                       0d0,      IGFED3,        JGFED3,     
     &                       1,        ARRAY_LANDMAP, QUIET=.TRUE. )

      ! Cast from REAL*4 to INTEGER
      HUMTROP_GFED3(:,:) = ARRAY_LANDMAP(:,:,1)
      
      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the GFED3 biomass species
      !=================================================================
      
      ! Save correspondance between GFED3 species order (N) and 
      ! species order of the simulation (IDBxxxs).(ccc, 2/4/10)
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( TRIM( GFED3_SPEC_NAME(N) ) ) 
            CASE( 'NOx'  )
               ! Replaced NOx with NO (mpayer, 3/14/13)
               BIO_SAVE(N)         = IDBNO
               GFED3_SPEC_MOLWT(N) = 14d-3
               GFED3_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO'   )
               BIO_SAVE(N)         = IDBCO
               GFED3_SPEC_MOLWT(N) = 28d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ALK4' )
               BIO_SAVE(N)         = IDBALK4
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ACET' )
               BIO_SAVE(N)         = IDBACET
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MEK'  )
               BIO_SAVE(N)         = IDBMEK
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ALD2' )
               BIO_SAVE(N)         = IDBALD2
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'PRPE' )
               BIO_SAVE(N)         = IDBPRPE
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C3H8' )
               BIO_SAVE(N)         = IDBC3H8
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'CH2O' )
               BIO_SAVE(N)         = IDBCH2O
               GFED3_SPEC_MOLWT(N) = 30d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'C2H6' )
               BIO_SAVE(N)         = IDBC2H6
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'SO2'  )
               BIO_SAVE(N)         = IDBSO2
               GFED3_SPEC_MOLWT(N) = 64d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'NH3'  )
               BIO_SAVE(N)         = IDBNH3
               GFED3_SPEC_MOLWT(N) = 17d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BC'   )
               !IDBBC = N
               BIO_SAVE(N)         = IDBBC
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'OC'   )
               BIO_SAVE(N)         = IDBOC
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYX' )
               BIO_SAVE(N)         = IDBGLYX
               GFED3_SPEC_MOLWT(N) = 58d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'MGLY' )
               BIO_SAVE(N)         = IDBMGLY
               GFED3_SPEC_MOLWT(N) = 72d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BENZ' )
               BIO_SAVE(N)         = IDBBENZ
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'TOLU' )
               BIO_SAVE(N)         = IDBTOLU
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'XYLE' )
               BIO_SAVE(N)         = IDBXYLE
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H4' )
               BIO_SAVE(N)         = IDBC2H4
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H2' )
               BIO_SAVE(N)         = IDBC2H2
               GFED3_SPEC_MOLWT(N) = 12d-3
               GFED3_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYC' )
               BIO_SAVE(N)         = IDBGLYC
               GFED3_SPEC_MOLWT(N) = 60d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'HAC' )
               BIO_SAVE(N)         = IDBHAC
               GFED3_SPEC_MOLWT(N) = 74d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CO2'  )
               BIO_SAVE(N)         = IDBCO2
               GFED3_SPEC_MOLWT(N) = 44d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CH4' )
               BIO_SAVE(N)         = IDBCH4
               GFED3_SPEC_MOLWT(N) = 16d-3
               GFED3_SPEC_UNIT(N)  = '[Tg  ]'
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0

            WRITE(*,*) 'NAME',TRIM( GFED3_SPEC_NAME(N) )
         END SELECT
      ENDDO

      END SUBROUTINE INIT_GFED3_BIOMASS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rearrange_biom
!
! !DESCRIPTION: Subroutine REARRANGE\_BIOM takes GFED3 emissions (which have 
!  their own, unique ID numbers and associates them with the IDBxxxs of 
!  tracerid\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE REARRANGE_BIOM( BIOM_OUT, BIOM_OUTM )

!
! !USES:
!
      USE CMN_SIZE_MOD     ! Size parameters
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)  :: BIOM_OUT (IIPAR,JJPAR,N_SPEC)
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT) :: BIOM_OUTM(IIPAR,JJPAR,NBIOMAX) !+1 from CO2
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      ! Loop over GFED3 species
      DO N = 1, N_SPEC

         ! Save into array w/ proper ordering for GEOS-Chem
         IF ( BIO_SAVE(N) .GT. 0 ) THEN
            BIOM_OUTM(:,:,BIO_SAVE(N)) = BIOM_OUT(:,:,N)
         ENDIF

      ENDDO

      END SUBROUTINE REARRANGE_BIOM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gfed3_biomass
!
! !DESCRIPTION: Subroutine CLEANUP\_GFED3\_BIOMASS deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GFED3_BIOMASS
! 
! !REVISION HISTORY: 
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GFED3_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( GFED3_EMFAC      ) ) DEALLOCATE( GFED3_EMFAC     )
      IF ( ALLOCATED( GFED3_SPEC_MOLWT ) ) DEALLOCATE( GFED3_SPEC_MOLWT)
      IF ( ALLOCATED( GFED3_SPEC_NAME  ) ) DEALLOCATE( GFED3_SPEC_NAME )
      IF ( ALLOCATED( HUMTROP_GFED3    ) ) DEALLOCATE( HUMTROP_GFED3   )
      IF ( ALLOCATED( BIOMASS_MODEL    ) ) DEALLOCATE( BIOMASS_MODEL   )
      
      END SUBROUTINE CLEANUP_GFED3_BIOMASS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_bpch2_gfed3
!
! !DESCRIPTION: Subroutine READ\_BPCH2\_GFED3 reads GFED3 DM burnt and
! and humid tropical forest map files 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BPCH2_GFED3( FILENAME, CATEGORY_IN, TRACER_IN,
     &                             TAU0_IN,  IX,          JX,
     &                             LX,       ARRAY,       QUIET )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR
      USE BPCH2_MOD, ONLY : OPEN_BPCH2_FOR_READ
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*),  INTENT(IN)  :: FILENAME         ! Bpch file to read
      CHARACTER(LEN=*),  INTENT(IN)  :: CATEGORY_IN      ! Diag. category name
      INTEGER,           INTENT(IN)  :: TRACER_IN        ! Tracer index #
      REAL*8,            INTENT(IN)  :: TAU0_IN          ! TAU timestamp
      INTEGER,           INTENT(IN)  :: IX, JX, LX       ! Dimensions of ARRAY
      LOGICAL, OPTIONAL, INTENT(IN)  :: QUIET            ! Don't print output
! 
! !OUTPUT PARAMETERS:
!
      REAL*4,            INTENT(OUT) :: ARRAY(IX,JX,LX)  ! Data array from file
! 
! !REVISION HISTORY:
!  (1 ) Adapted from READ_BPCH2 to facilitate reading of 0.5x0.5 GFED3 files (psk, 2/7/12)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: FOUND, TMP_QUIET
      INTEGER            :: IU_FILE
      INTEGER            :: I,  J,  L,  N,  IOS, M
      INTEGER            :: I1, I2, J1, J2, L1,  L2
      CHARACTER(LEN=255) :: MSG

      ! Make TEMPARRAY big enough for GFED3 files - 0.5x0.5 lat-lon grid
      REAL*4             :: TEMPARRAY_GFED3(720,360,1)

      ! For binary punch file, version 2.0
      INTEGER            :: NTRACER,   NSKIP
      INTEGER            :: HALFPOLAR, CENTER180
      INTEGER            :: NI,        NJ,        NL
      INTEGER            :: IFIRST,    JFIRST,    LFIRST
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: CATEGORY
      CHARACTER(LEN=40)  :: UNIT
      CHARACTER(LEN=40)  :: RESERVED

      !=================================================================
      ! READ_BPCH2_GFED3 begins here!
      !
      ! Initialize some variables
      !=================================================================
      FOUND            = .FALSE.
      ARRAY(:,:,:)     = 0e0
      TEMPARRAY_GFED3(:,:,:) = 0e0

      ! Define a temporary variable for QUIET
      IF ( PRESENT( QUIET ) ) THEN
         TMP_QUIET = QUIET
      ELSE
         TMP_QUIET = .FALSE.
      ENDIF

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      !=================================================================
      ! Open binary punch file and read top-of-file header.
      ! Do some error checking to make sure the file is the right format.
      !=================================================================
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      !=================================================================
      ! Read data from the binary punch file
      !
      ! NOTE: IOS < 0 is end-of-file, IOS > 0 is error condition
      !=================================================================
      DO
         READ( IU_FILE, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:4' )

         READ( IU_FILE, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:5' )

         READ( IU_FILE, IOSTAT=IOS )
     &        ( ( ( TEMPARRAY_GFED3(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_bpch2:6' )

         ! Test for a match
         IF ( TRIM( CATEGORY_IN ) == TRIM( CATEGORY ) .and.
     &        TRACER_IN           == NTRACER          .and.
     &        TAU0_IN             == ZTAU0 ) THEN
            FOUND = .TRUE.
            EXIT
         ENDIF

      ENDDO

      IF ( FOUND ) THEN

         I1 = IFIRST
         J1 = JFIRST
         L1 = LFIRST


         I2 = NI + I1 - 1
         J2 = NJ + J1 - 1
         L2 = NL + L1 - 1

         ARRAY( I1:I2, J1:J2, L1:L2 ) 
     &   = TEMPARRAY_GFED3( 1:NI, 1:NJ, 1:NL )

         ! Flag to decide whether or not we will echo info (bmy, 3/14/03)
         IF ( .not. TMP_QUIET ) THEN
            WRITE( 6, 100 ) ZTAU0, NTRACER
 100        FORMAT( 'READ_BPCH2_GFED3: Found data for TAU = ', f10.2,
     &              ' and tracer # ', i6 )
         ENDIF

      ELSE
         MSG = 'No matches found for file ' // TRIM( FILENAME ) // '!'
         CALL ERROR_STOP( MSG, 'READ_BPCH2_GFED3 (bpch2_mod.f)!' )
      ENDIF

      !=================================================================
      ! Close file and quit
      !=================================================================
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_BPCH2_GFED3
!EOC

      END MODULE GFED3_BIOMASS_MOD
