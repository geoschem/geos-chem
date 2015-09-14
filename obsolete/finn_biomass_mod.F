!------------------------------------------------------------------------------
!          J. Mao (NOAA/GFDL) & J.A. Fisher (U. Wollongong)                   !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: finn_biomass_mod
!     
! !DESCRIPTION: Module FINN\_BIOMASS\_MOD contains routines and variables 
!  used to incorporate FINN emissions into GEOS-Chem
!\\   
!\\   
! !INTERFACE: 
!
      MODULE FINN_BIOMASS_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: FINN_COMPUTE_BIOMASS
      PUBLIC  :: CLEANUP_FINN_BIOMASS
      PUBLIC  :: FINN_IS_NEW
!
! PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CHECK_FINN
      PRIVATE :: FINN_AVAILABLE
      PRIVATE :: READ_FINN
      PRIVATE :: FINN_SCALE_FUTURE
      PRIVATE :: FINN_TOTAL_Tg
      PRIVATE :: INIT_FINN_BIOMASS
      PRIVATE :: REARRANGE_BIOM
!
! !REMARKS:
!  Emissions of biomass burning species are read at monthly or daily
!  resolution. Note: no emission factors are used here - emissions of
!  individual species are given in input files. Emissions on the FINN 0.5x0.5 
!  degree grid are regridded to the current GEOS-Chem or GCAP grid.
!                                                                             .
!  FINN biomass burning emissions are computed for the following gas-phase 
!  and aerosol-phase species:
!                                                                             .
!     (1 ) NOx  [  molec/cm2/s]     (13) BC   [atoms C/cm2/s]
!     (2 ) CO   [  molec/cm2/s]     (14) OC   [atoms C/cm2/s]                  
!     (3 ) ALK4 [atoms C/cm2/s]     (15) MGLY [  molec/cm2/s]    
!     (4 ) ACET [atoms C/cm2/s]     (16) BENZ [atoms C/cm2/s]  
!     (5 ) MEK  [atoms C/cm2/s]     (17) TOLU [atoms C/cm2/s]     
!     (6 ) ALD2 [atoms C/cm2/s]     (18) C2H4 [atoms C/cm2/s]
!     (7 ) PRPE [atoms C/cm2/s]     (19) C2H2 [atoms C/cm2/s]
!     (8 ) C3H8 [atoms C/cm2/s]     (20) GLYC [  molec/cm2/s]
!     (9 ) CH2O [  molec/cm2/s]     (21) HAC  [  molec/cm2/s]
!     (10) C2H6 [atoms C/cm2/s]     (22) CO2  [  molec/cm2/s]
!     (11) SO2  [  molec/cm2/s]     (23) CH4  [  molec/cm2/s]
!     (12) NH3  [  molec/cm2/s]     (24) 
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Original FINN database from Christine Wiedinmyer
!        http://bai.acd.ucar.edu/Data/fire/
!  (2 ) Wiedinmyer, C., Akagi, S.K., Yokelson, R.J., Emmons, L.K.,
!       Al-Saadi, J.A., Orlando, J.J., and Soja, A.J.: The Fire
!       INventory from NCAR (FINN): a high resolution global model to
!       estimate the emissions from open burning, Geoscientific Model
!       Development, 4, 625-641, doi:10.5194/gmd-4-625-2011, 2011.
!
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J.A. Fisher - Initial version, based on GFED3
!  01 Oct 2013 - J.A. Fisher - Update to only use one input file
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
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
      ! I_NATIVE: Longitude dimension of input files
      ! J_NATIVE: Latitude dimension of input files
      !=================================================================
      INTEGER, PARAMETER   :: N_EMFAC = 6
      INTEGER, PARAMETER   :: N_SPEC  = 39
      INTEGER, PARAMETER   :: I_NATIVE = 720 ! 0.5 x 0.5 resolution
      INTEGER, PARAMETER   :: J_NATIVE = 360 ! 0.5 x 0.5 resolution

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! UPDATED         : flag to indicate if FINN emissions are updated
      ! SECONDS         : Number of seconds in the current month
      !
      ! Arrays
      !
      ! FINN_SPEC_NAME  : Array for FINN biomass species names
      ! FINN_SPEC_MOLWT : Array for FINN biomass species molecular wts
      ! FINN_SPEC_UNIT  : Array for FINN biomass species emissions units
      ! FINN_EMFAC      : Array for emission factors
      ! BIOMASS_MODEL   : Array for FINN species emissions on model grid
      ! BIO_SAVE        : Index array to store IDBxxx values
      ! LUMP_ID         : Index array to store IDBxxx for lumped species
      !=================================================================

      ! Scalars
      LOGICAL                       :: UPDATED
      REAL*8                        :: SECONDS

      ! Arrays
      REAL*8,           ALLOCATABLE :: FINN_SPEC_MOLWT(:)
      CHARACTER(LEN=6), ALLOCATABLE :: FINN_SPEC_NAME(:)
      CHARACTER(LEN=6), ALLOCATABLE :: FINN_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: FINN_EMFAC(:,:)
      REAL*8,           ALLOCATABLE :: BIOMASS_MODEL(:,:,:)
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)
      INTEGER,          ALLOCATABLE :: LUMP_ID(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finn_is_new
!
! !DESCRIPTION: Function FINN\_IS\_NEW returns TRUE if FINN emissions 
!  have been updated.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FINN_IS_NEW( ) RESULT( IS_UPDATED )
!
! !RETURN VALUE:
!
      LOGICAL :: IS_UPDATED    ! =T if FINN is updated; =F otherwise
!
! !REMARKS:
!  Called from carbon_mod.f and sulfate_mod.f
!
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
      IS_UPDATED = UPDATED      

      END FUNCTION FINN_IS_NEW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_finn
!
! !DESCRIPTION: Subroutine CHECK\_FINN checks if we entered a new FINN period
!  since last emission timestep (ie, last call). The result depends
!  on the emissions time step, and the FINN time period used, as well
!  as MMDD at beginning of the GEOS-Chem run
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_FINN( DOY )
!
! !USES:
!
      USE LOGICAL_MOD, ONLY : L1DAYBBF
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,    ONLY : ITS_A_NEW_DAY
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: DOY   ! Day of year (0-365 or 0-366 leap years)
!
! !REMARKS:
!  The routine computes the DOY at start of the 1-day period we are in, 
!  if the 1-day FINN. option is on. Result is compared to previous value
!  to indicate if new data should be read.
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! Reset to default
      UPDATED     = .FALSE.

      ! Check if it is a new month
      IF ( ITS_A_NEW_MONTH() ) UPDATED = .TRUE.
      
      ! Check if it is a new day
      IF ( L1DAYBBF .AND. ITS_A_NEW_DAY() ) UPDATED = .TRUE.

      END SUBROUTINE CHECK_FINN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finn_available
!
! !DESCRIPTION: Function FINN\_AVAILABLE checks an input YYYY year and MM
!  month against the available data dates.  If the requested YYYY and MM
!  lie outside of the valid range of dates, then FINN\_AVAILABLE will return 
!  the last valid YYYY and MM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FINN_AVAILABLE( YYYY, YMIN, YMAX, MM, MMIN, MMAX )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)              :: YMIN, YMAX   ! Min & max years 
      INTEGER, INTENT(IN),    OPTIONAL :: MMIN, MMAX   ! Min & max months
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT)           :: YYYY         ! Year of FINN data
      INTEGER, INTENT(INOUT), OPTIONAL :: MM           ! Month of FINN data
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Check year
      IF ( YYYY > YMAX .OR. YYYY < YMIN ) THEN
         
         YYYY = MAX( YMIN, MIN( YYYY, YMAX) )
         
         WRITE( 6, 100 ) YMAX, YMIN, YYYY
 100     FORMAT( 'YEAR > ', i4, ' or YEAR < ', i4, 
     $        '. Using FINN biomass for ', i4)
      ENDIF
      

      ! Check month
      IF ( PRESENT( MM ) ) THEN 
         IF ( MM > MMAX .OR. MM < MMIN ) THEN

            MM = MAX( MMIN, MIN( MM, MMAX) )
            
            WRITE( 6, 200 ) MMIN, MMAX, MM
 200        FORMAT( ' ** WARNING ** : MONTH is not within ', i2,'-',
     $              i2, '. Using FINN biomass for month #', i2)
         ENDIF
      ENDIF

      END SUBROUTINE FINN_AVAILABLE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finn_compute_biomass
!
! !DESCRIPTION: Subroutine FINN\_COMPUTE\_BIOMASS computes the monthly 
!  FINN biomass burning emissions for a given year and month. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FINN_COMPUTE_BIOMASS( THIS_YYYY, THIS_MM, BIOM_OUT )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_TAU0
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE JULDAY_MOD,      ONLY : JULDAY
      USE JULDAY_MOD,      ONLY : CALDATE
      USE LOGICAL_MOD,     ONLY : LFUTURE
      USE LOGICAL_MOD,     ONLY : L1DAYBBF
      USE REGRID_A2A_MOD,  ONLY : DO_REGRID_A2A
      USE TIME_MOD,        ONLY : EXPAND_DATE
      USE TIME_MOD,        ONLY : GET_DAY
      USE TIME_MOD,        ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,        ONLY : YMD_EXTRACT

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
! !REVISION HISTORY: 
!  02 Jan 2012 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.

      ! For reading
      INTEGER             :: N, NF
      INTEGER             :: YYYY, MM, MM1, YYYY1
      INTEGER             :: YYYYMMDD, YYYYMMDD1, HHMMSS
      REAL*8              :: ARRAY_FINN(IIPAR, JJPAR)
      REAL*8              :: ARRAY_FINN_NATIVE(I_NATIVE, J_NATIVE,1)
      REAL*8              :: CO2_FINN_IN(I_NATIVE, J_NATIVE, N_EMFAC)
      REAL*8              :: TAU0, TAU1
      REAL*8              :: JD,   JD1
      INTEGER             :: DD,   DD1,  DOY
      CHARACTER(LEN=255)  :: FINN_DIR
      CHARACTER(LEN=255)  :: FILENAME1
      CHARACTER(LEN=255)  :: FILENAME2
      CHARACTER(LEN=255)  :: FILENAME3
      CHARACTER(LEN=255)  :: FILENAME4
      CHARACTER(LEN=255)  :: FILENAME5
      CHARACTER(LEN=255)  :: FILENAME6
      CHARACTER(LEN=7  )  :: TIME_STR
      CHARACTER(LEN=4  )  :: SPEC_NAME

      ! For regridding
      CHARACTER(LEN=255)  :: LLFILENAME
      REAL*8              :: INGRID(I_NATIVE,J_NATIVE)

      !=================================================================
      ! FINN_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_FINN_BIOMASS
         FIRST = .FALSE.
      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM
      DD   = GET_DAY()
      DOY  = GET_DAY_OF_YEAR()
      
      ! Check if we need to update FINN 
      CALL CHECK_FINN( DOY )

      ! If no updating is needed, module variable BIOMASS_MODEL
      ! from last update can be used
      IF ( .not. UPDATED ) THEN
         CALL REARRANGE_BIOM(BIOMASS_MODEL,BIOM_OUT)
         RETURN
      ENDIF

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     &  'F I N N   B I O M A S S   B U R N I N G   E M I S S I O N S'

      !=================================================================
      ! Set-up for reading FINN daily/monthly emissions
      !=================================================================
      
      ! Check availability
      CALL FINN_AVAILABLE( YYYY, 2007, 2013 )

      ! DAILY data
      !-------------------------------
      IF ( L1DAYBBF ) THEN

         ! Create YYYYMMDD integer values
         YYYYMMDD = YYYY*10000 + MM*100 + DD

         ! Use Julian dates to automatically increase year/month
         JD  = JULDAY( YYYY, MM, DBLE(DD) )
         JD1 = JD + 1d0

         CALL CALDATE( JD1, YYYYMMDD1, HHMMSS )
         CALL YMD_EXTRACT( YYYYMMDD1, YYYY1, MM1, DD1 )

         ! TAU value at start of this day
         TAU0     = GET_TAU0( MM, DD, YYYY )
         ! TAU value at start of next day
         TAU1     = GET_TAU0( MM1, DD1, YYYY1 )

         ! Directory with FINN emissions
         FINN_DIR = 'FINN_201310/daily/'
         TIME_STR = 'daily'

      ! MONTHLY data
      !-------------------------------
      ELSE 

         ! Create YYYYMMDD integer values
         YYYYMMDD = YYYY*10000 + MM*100 + 01

         ! Time value is set to 1st of YYYY/MM.
         DD       = 1

         ! Get YYYY/MM value for next month
         MM1      = MM + 1
         YYYY1    = YYYY

         ! Increment year if necessary
         IF ( MM1 == 13 ) THEN
            MM1   = 1
            YYYY1 = YYYY + 1
         ENDIF

         ! TAU value at start of this month
         TAU0     = GET_TAU0( MM, 1, YYYY )
         ! TAU value at start of next month
         TAU1     = GET_TAU0( MM1, 1, YYYY1 )

         ! Directory with FINN emissions
         FINN_DIR = 'FINN_201310/monthly/'
         TIME_STR = 'monthly'

      ENDIF

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_NATIVE) //
     &             'MAP_A2A_Regrid_201203/' //
     &             'MAP_A2A_latlon_generic05x05.nc'

      !=================================================================
      ! Calculate biomass species emissions on FINN grid
      ! and regrid to model grid
      !
      ! Emission are in [molec/cm2/s] and converted here to
      ! [molec/cm2/timeperiod] or [atoms C/cm2/timeperiod]
      !
      ! Units:
      !  [  molec/cm2/month] : NOx,  CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/month] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8,
      !                        C2H6, BC,   OC
      !=================================================================

      ! Zero CO2 array
      CO2_FINN_IN = 0d0

      ! File names with FINN CO2 emissions
      FILENAME1 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype1_YYYY_05x05.nc'
      FILENAME2 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype2_YYYY_05x05.nc'
      FILENAME3 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype3_YYYY_05x05.nc'
      FILENAME4 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype4_YYYY_05x05.nc'
      FILENAME5 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype5_YYYY_05x05.nc'
      FILENAME6 = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
     &            TRIM( TIME_STR ) // '_CO2_vegtype9_YYYY_05x05.nc'

      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME1, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME2, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME3, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME4, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME5, YYYYMMDD, 000000 )
      CALL EXPAND_DATE( FILENAME6, YYYYMMDD, 000000 )

      ! Read FINN CO2 emissions [molec CO2/cm2/s] for each land type
      CALL READ_FINN( FILENAME1, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,1) = ARRAY_FINN_NATIVE(:,:,1)

      CALL READ_FINN( FILENAME2, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,2) = ARRAY_FINN_NATIVE(:,:,1)

      CALL READ_FINN( FILENAME3, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,3) = ARRAY_FINN_NATIVE(:,:,1)

      CALL READ_FINN( FILENAME4, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,4) = ARRAY_FINN_NATIVE(:,:,1)

      CALL READ_FINN( FILENAME5, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,5) = ARRAY_FINN_NATIVE(:,:,1)

      CALL READ_FINN( FILENAME6, ARRAY_FINN_NATIVE, YYYY, MM, DD )
      CO2_FINN_IN(:,:,6) = ARRAY_FINN_NATIVE(:,:,1)

      !=================================================================
      ! Calculate biomass species emissions on native FINN grid
      !
      ! Emission factors are given in [mole X]/[mole CO2]
      ! Some species need to be converted from [molec X] to [atoms C X]
      !=================================================================

      ! Loop over biomass species 
      DO N = 1, N_SPEC

         ! Initialize temporary arrays
         ARRAY_FINN        = 0d0
         ARRAY_FINN_NATIVE = 0d0

         ! Species with no emission factor have FINN_EMFAC=0
         ! Species that are being lumped with another will have
         ! BIO_SAVE=0 but still need to be calculated here, so use
         ! LUMP_ID as optional criterion.
         IF ( MAXVAL(FINN_EMFAC(N,:)) > 0 .AND.
     &        ( BIO_SAVE(N) > 0 .OR. LUMP_ID(N) > 0 ) ) THEN

         ! Loop over land types
         DO NF = 1, N_EMFAC

            ! Multiply CO2 emissions by appropriate ratio for each land
            ! type and sum to get total emissions for the species on the
            ! native grid - emissions are in molec/cm2/s

            ARRAY_FINN_NATIVE(:,:,1) = ARRAY_FINN_NATIVE(:,:,1) + 
     &                                 CO2_FINN_IN(:,:,NF)      *
     &                                 FINN_EMFAC(N,NF)

         ENDDO

         ! Regrid to current model resolution
         INGRID = ARRAY_FINN_NATIVE(:,:,1)
         CALL DO_REGRID_A2A( LLFILENAME, I_NATIVE,   J_NATIVE,
     &                       INGRID,     ARRAY_FINN, IS_MASS=0,
     &                       netCDF=.TRUE.                      )

         ! Add to biomass array
         BIOMASS_MODEL(:,:,N) = ARRAY_FINN

         ENDIF ! EMFAC > 0

      ENDDO

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL FINN_SCALE_FUTURE( BIOMASS_MODEL )
      ENDIF

      ! Print totals in Tg/time period
      IF ( L1DAYBBF ) THEN
         WRITE( 6, 411 ) YYYY, MM, DD
 411     FORMAT( 'FINN daily emissions for year, month, day: ',
     &            i4, '/', i2.2, '/', i2.2, / )
      ELSE
         WRITE( 6, 410 ) YYYY, MM
 410     FORMAT( 'FINN monthly emissions for year, month: ',
     &            i4, '/', i2.2, / )
      ENDIF
      ! Number of seconds in this period (day or month)
      ! (NOTE: its value will be saved until the next call)
      SECONDS  = ( TAU1 - TAU0 ) * 3600d0

      CALL FINN_TOTAL_Tg( SECONDS )

      ! Rearrange the species to the same order as in the IDBxxx (fp, 6/09)
      ! BIOMASS_MODEL is indexed as FINN
      ! BIOM_OUT      is indexed as IDBs
      CALL REARRANGE_BIOM( BIOMASS_MODEL, BIOM_OUT )

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE FINN_COMPUTE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!         J. Fisher - U. Wollongong                                           !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_finn
!
! !DESCRIPTION: Subroutine READ\_FINN\_NATIVE reads the FINN BB emissions at
!  0.5x0.5 degree resolution
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_FINN(FILENAME, ARR, YYYY, MM, DD)
!
! !USES:
!
      USE DIRECTORY_MOD,  ONLY: DATA_DIR_NATIVE => DATA_DIR_1x1
      USE NCDF_MOD,       ONLY: NC_READ

      USE CMN_SIZE_MOD          ! Size parameters

!
! !INPUT PARAMETERS:
!
      ! File and parameter to read
      CHARACTER(LEN=255),  INTENT(IN) :: FILENAME
      INTEGER,             INTENT(IN) :: YYYY
      INTEGER,             INTENT(IN) :: MM
      INTEGER,             INTENT(IN) :: DD
! 
! !OUTPUT PARAMETERS:
!
      ! Emissions output
      REAL*8,            INTENT(OUT) :: ARR(I_NATIVE,J_NATIVE,1)
! 
! !REVISION HISTORY:
!  02 Jan 2013 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! for ncdf-reading
      REAL*4, POINTER       :: NCARR(:,:,:) => NULL()

      !=================================================================
      ! READ_FINN begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_FINN : Reading ', a )

      ARR   = 0D0 

      CALL NC_READ ( NC_PATH = TRIM(FILENAME) ,
     &               PARA    = 'fire'         ,
     &               ARRAY   = NCARR          ,
     &               YEAR    = YYYY           ,
     &               MONTH   = MM             ,
     &               DAY     = DD             ,
     &               VERBOSE = .FALSE.        )

      ! Cast to REAL*8
      ARR(:,:,1) = NCARR(:,:,1)

      ! Return to calling program
      END SUBROUTINE READ_FINN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finn_scale_future
!
! !DESCRIPTION: Subroutine FINN\_SCALE\_FUTURE applies the IPCC future 
!  emissions scale factors to the FINN biomass burning emisisons in order 
!  to compute the future emissions of biomass burning for NOx, CO, and VOC's.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FINN_SCALE_FUTURE( BB )
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
      USE TRACERID_MOD,           ONLY : IDBNO,  IDBCO,   IDBSO2 
      USE TRACERID_MOD,           ONLY : IDBNH3,  IDBBC,   IDBOC 

      USE CMN_SIZE_MOD                 ! Size parameters
!
! !OUTPUT PARAMETERS:
!
      !  Array w/ biomass burning emisisons [molec/cm2]
      REAL*8, INTENT(INOUT) :: BB(IIPAR,JJPAR,N_SPEC)
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Fisher - Initial version, based on GFED3
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
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
      ! FINN_SCALE_FUTURE begins here!
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

      END SUBROUTINE FINN_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: finn_total_Tg
!
! !DESCRIPTION: Subroutine FINN\_TOTAL\_Tg prints the amount of biomass 
!  burning emissions that are emitted each month/day in Tg or Tg C. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FINN_TOTAL_Tg( SECONDS )
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD      ! Size parameters
!
! !INPUT PARAMETERS:
!
      ! Seconds in this time period
      REAL*8,  INTENT(IN)  :: SECONDS
!
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
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
      ! FINN_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         IF ( BIO_SAVE(N) > 0 .OR. LUMP_ID(N) > 0 ) THEN

         ! Initialize
         NAME  = FINN_SPEC_NAME(N)
         MOLWT = FINN_SPEC_MOLWT(N)
         UNIT  = FINN_SPEC_UNIT(N)
         TOTAL = 0d0

         ! Loop over surface boxes
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Convert to [Tg/finn-period] (or [Tg C/finn-period] for HC's)
            CONV  = GET_AREA_CM2( I, J, 1 ) 
     &            * ( MOLWT / 6.023d23 ) * 1d-9 * SECONDS

            TOTAL = TOTAL + ( BIOMASS_MODEL(I,J,N) * CONV )

         ENDDO
         ENDDO
     
         ! Write totals
         IF ( LUMP_ID(N) == 0 ) THEN
            WRITE( 6, 110 ) NAME, TOTAL, UNIT
         ELSE
            WRITE( 6, 111 ) NAME, TOTAL, UNIT, LUMP_ID(N)
         ENDIF

         ENDIF

      ENDDO

 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', f11.4, 1x, a6 )
 111     FORMAT( 'Sum Biomass ', a4, 1x, ': ', f11.4, 1x, a6,
     &           ' (lumped with BB tracer # ', i2,')')

      END SUBROUTINE FINN_TOTAL_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_finn_biomass
!
! !DESCRIPTION: Subroutine INIT\_FINN\_BIOMASS allocates all module arrays.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_FINN_BIOMASS
!
! !USES:
!
      USE CHARPAK_MOD,     ONLY : STRSPLIT
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE FILE_MOD,        ONLY : IOERROR
      USE inquireMod,      ONLY : findFreeLUN
      USE TRACER_MOD,      ONLY : TRACER_MW_KG
      USE TRACER_MOD,      ONLY : ITS_A_CO2_SIM, ITS_A_CH4_SIM
      USE TRACERID_MOD
      USE CMN_SIZE_MOD          ! Size parameters
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, PARAMETER :: N_SPEC_EMFAC = 12
      INTEGER, PARAMETER :: N_NMOC  = 28
      INTEGER            :: IU_FILE
      INTEGER            :: AS, IOS, M, N, NDUM, N_SPECSTRS, N_NMOCSTRS
      LOGICAL            :: IS_NMOC
      CHARACTER(LEN=255) :: FILENAME, ADUM
      CHARACTER(LEN=255) :: IN_SPEC_NAME(255)
      CHARACTER(LEN=255) :: IN_NMOC_NAME(255)
      CHARACTER(LEN=255) :: TMPNAME
      REAL*8             :: C_MOLEC
      REAL*8             :: EMFAC_IN(N_SPEC_EMFAC, N_EMFAC)
      REAL*8             :: NMOC_EMFAC(N_EMFAC), NMOC_RATIO(N_EMFAC)
      REAL*8             :: NMOC_RATIO_IN(N_NMOC, N_EMFAC)
      
      !=================================================================
      ! INIT_FINN_BIOMASS begins here!
      !=================================================================

      ! Allocate array to hold FINN species emissions on model grid
      ALLOCATE( BIOMASS_MODEL( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_MODEL' )
      BIOMASS_MODEL = 0d0

      ! Allocate array for emission factors (default = 0)
      ALLOCATE( FINN_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_EMFAC' )
      FINN_EMFAC = 0d0

      ! Allocate array for species molecular weight
      ALLOCATE( FINN_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_MOLWT' )
      FINN_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( FINN_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_NAME' )
      FINN_SPEC_NAME = ''

      ! Allocate array for FINN biomass burning species mass units
      ALLOCATE( FINN_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_UNIT' )
      FINN_SPEC_UNIT = '[Tg  ]'

      !IDBs are now the same as the ones in TRACERID AND BIOMASS_MOD
      !BIOSAVE INDEX IS THE LOCATION OF THE EMISSION IN THE FINN FILE
      ALLOCATE( BIO_SAVE( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIO_SAVE' )
      BIO_SAVE = 0

      !Some species are currently lumped with others. LUMP_ID saves
      !the ID of the primary species that secondary species will be
      !lumped into.
      ALLOCATE( LUMP_ID( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LUMP_ID' )
      LUMP_ID = 0

      !=================================================================
      ! Define FINN species names
      !=================================================================
      FINN_SPEC_NAME(1)  = 'CO2'
      FINN_SPEC_NAME(2)  = 'CO'
      FINN_SPEC_NAME(3)  = 'CH4'
      FINN_SPEC_NAME(4)  = 'NOx'
      FINN_SPEC_NAME(5)  = 'SO2'
      FINN_SPEC_NAME(6)  = 'OC'
      FINN_SPEC_NAME(7)  = 'BC'
      FINN_SPEC_NAME(8)  = 'NH3'
      FINN_SPEC_NAME(9)  = 'ACET'
      FINN_SPEC_NAME(10) = 'ACTA'   ! Not currently emitted by BB in GC
      FINN_SPEC_NAME(11) = 'ALD2'
      FINN_SPEC_NAME(12) = 'ALK4'
      FINN_SPEC_NAME(13) = 'BENZ'
      FINN_SPEC_NAME(14) = 'C2H2'
      FINN_SPEC_NAME(15) = 'C2H4'
      FINN_SPEC_NAME(16) = 'C2H6'
      FINN_SPEC_NAME(17) = 'C3H8'
      FINN_SPEC_NAME(18) = 'CH2O'
      FINN_SPEC_NAME(19) = 'EOH'    ! Not currently emitted in GC
      FINN_SPEC_NAME(20) = 'GLYC'
      FINN_SPEC_NAME(21) = 'GLYX'
      FINN_SPEC_NAME(22) = 'HAC'
      FINN_SPEC_NAME(23) = 'HCN'    ! Not currently emitted in GC
      FINN_SPEC_NAME(24) = 'HCOOH'  ! Not currently emitted by BB in GC
      FINN_SPEC_NAME(25) = 'ISOP'   ! Not currently emitted by BB in GC
      FINN_SPEC_NAME(26) = 'MACR'   ! Not currently emitted in GC
      FINN_SPEC_NAME(27) = 'MEK'
      FINN_SPEC_NAME(28) = 'MGLY'
      FINN_SPEC_NAME(29) = 'MOH'    ! Not currently emitted in GC
      FINN_SPEC_NAME(30) = 'MVK'    ! Not currently emitted in GC
      FINN_SPEC_NAME(31) = 'PRPE'
      FINN_SPEC_NAME(32) = 'R4N2'   ! Not currently emitted in GC
      FINN_SPEC_NAME(33) = 'RCHO'   ! Not currently emitted by BB in GC
      FINN_SPEC_NAME(34) = 'TOLU'
      FINN_SPEC_NAME(35) = 'XYLE'
      FINN_SPEC_NAME(36) = 'HNO2'   ! Not currently emitted in GC
      FINN_SPEC_NAME(37) = 'TMB'    ! Currently lumped with XYLE
      FINN_SPEC_NAME(38) = 'ETBENZ' ! Currently lumped with TOLU
      FINN_SPEC_NAME(39) = 'STYR'   ! Currently lumped with TOLU

      !=================================================================
      ! Read emission factors ([mole CO2]/[mole X])
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE) //
     &           'FINN_201310/FINN_EFratios_CO2.csv'

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:1' )

      ! Skip unnecessary header lines
      DO N = 1, 2
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:2' )
      ENDDO

      ! Read species names for emission ratio file
      READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
      CALL STRSPLIT(ADUM,',',IN_SPEC_NAME,N_SPECSTRS)

      ! Read emission factors for each species and land type
      DO N = 1, N_EMFAC
         READ( IU_FILE, *, IOSTAT=IOS ) NDUM, ADUM, EMFAC_IN(:,N)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:3' )
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Read NMOC factors ([mole X]/[kg NMOC])
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE) //
     &           'FINN_201310/FINN_VOC_speciation.csv'

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:4' )

      ! Skip unnecessary header lines
      DO N = 1, 2
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:5' )
      ENDDO

      ! Read species names for emission ratio file
      READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
      CALL STRSPLIT(ADUM,',',IN_NMOC_NAME,N_NMOCSTRS)

      ! Read emission factors for each species and land type
      DO N = 1, N_EMFAC
         READ( IU_FILE, *, IOSTAT=IOS ) NDUM, ADUM, NMOC_RATIO_IN(:,N)
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_finn:6' )
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the FINN biomass species
      !=================================================================
      ! Save correspondance between FINN species order (N) and 
      ! species order of the simulation (IDBxxxs).(ccc, 2/4/10)
      ! and also initialize arrays for mol wts and units
      ! Also get appropriate emission ratios to CO2 (jaf, 10/2/13)
      DO N = 1, N_SPEC

         ! Reset variables
         IS_NMOC    = .FALSE.
         C_MOLEC    = 1d0
         NMOC_RATIO = 0d0

         DO M = 1, N_SPECSTRS 
            TMPNAME = IN_SPEC_NAME(M)
            IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME(5:8)) ) THEN
               ! First two entries are not species. Also, EMFAC is
               ! stored as [mole CO2]/[mole X], but we want the inverse.
               FINN_EMFAC(N,:) = 1d0 / EMFAC_IN(M-2,:)
               WRITE( 6, 200 ) TRIM( FINN_SPEC_NAME(N) )
               EXIT
            ELSE IF ( TRIM(TMPNAME(5:8)) == 'NMOC' ) THEN
               NMOC_EMFAC = 1d0 / EMFAC_IN(M-2,:)
            ENDIF
         ENDDO

 200     FORMAT( 'Found FINN emission ratio for tracer ',a5 )

         DO M = 1, N_NMOCSTRS 
            TMPNAME = IN_NMOC_NAME(M)
            IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME) ) THEN
               ! First two entries are not species
               NMOC_RATIO = NMOC_RATIO_IN(M-2,:)
               IS_NMOC = .TRUE.
               WRITE( 6, 201 ) TRIM( FINN_SPEC_NAME(N) )
               EXIT
            ENDIF
         ENDDO

 201     FORMAT( 'Found FINN NMOC factor for tracer ',a5 )

         SELECT CASE ( TRIM( FINN_SPEC_NAME(N) ) ) 
            CASE( 'NOx'  )
               ! Replaced NOx with NO (jaf, 5/3/14)
               BIO_SAVE(N)        = IDBNO
               FINN_SPEC_MOLWT(N) = 14d-3 ! Tg N
               FINN_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO'   )
               BIO_SAVE(N)        = IDBCO
               IF ( IDTCO > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTCO)
            CASE( 'ALK4' )
               BIO_SAVE(N)        = IDBALK4
               IF ( IDTALK4 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTALK4)
               C_MOLEC            = 4d0
            CASE( 'ACET' )
               BIO_SAVE(N)        = IDBACET
               IF ( IDTACET > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTACET)
               C_MOLEC            = 3d0
            CASE( 'MEK'  )
               BIO_SAVE(N)        = IDBMEK
               IF ( IDTMEK > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTMEK)
               C_MOLEC            = 4d0
            CASE( 'ALD2' )
               BIO_SAVE(N)        = IDBALD2
               IF ( IDTALD2 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTALD2)
               C_MOLEC            = 2d0
            CASE( 'PRPE' )
               BIO_SAVE(N)        = IDBPRPE
               IF ( IDTPRPE > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTPRPE)
               C_MOLEC            = 3d0
            CASE( 'C3H8' )
               BIO_SAVE(N)        = IDBC3H8
               IF ( IDTC3H8 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTC3H8)
               C_MOLEC            = 3d0
            CASE( 'CH2O' )
               BIO_SAVE(N)        = IDBCH2O
               IF ( IDTCH2O > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTCH2O)
            CASE( 'C2H6' )
               BIO_SAVE(N)        = IDBC2H6
               IF ( IDTC2H6 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTC2H6)
               C_MOLEC            = 2d0
            CASE( 'SO2'  )
               BIO_SAVE(N)        = IDBSO2
               IF ( IDTSO2 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTSO2)
            CASE( 'NH3'  )
               BIO_SAVE(N)        = IDBNH3
               IF ( IDTNH3 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTNH3)
            CASE( 'MGLY' )
               BIO_SAVE(N)        = IDBMGLY
               IF ( IDTMGLY > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTMGLY)
            CASE( 'BENZ' )
               BIO_SAVE(N)        = IDBBENZ
               IF ( IDTBENZ > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTBENZ)
               C_MOLEC            = 6d0
            CASE( 'TOLU' )
               BIO_SAVE(N)        = IDBTOLU
               IF ( IDTTOLU > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTTOLU)
               C_MOLEC            = 7d0
            CASE( 'XYLE' )
               BIO_SAVE(N)        = IDBXYLE
               IF ( IDTXYLE > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTXYLE)
               C_MOLEC            = 8d0
            CASE( 'C2H4' )
               BIO_SAVE(N)        = IDBC2H4
               IF ( IDTC2H4 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTC2H4)
               C_MOLEC            = 2d0
            CASE( 'C2H2' )
               BIO_SAVE(N)        = IDBC2H2
               IF ( IDTC2H2 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTC2H2)
               C_MOLEC            = 2d0
            CASE( 'GLYC' )
               BIO_SAVE(N)        = IDBGLYC
               IF ( IDTGLYC > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTGLYC)
            CASE( 'GLYX' )
               BIO_SAVE(N)        = IDBGLYX
               IF ( IDTGLYX > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTGLYX)
            CASE( 'HAC' )
               BIO_SAVE(N)        = IDBHAC
               IF ( IDTHAC > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTHAC)

            ! Species for which FINN emissions are available but not yet
            ! used in any GEOS-Chem simulations. Use of these emissions
            ! will require user to create additional code for backbone
            ! of emissions. If these species are carried as atoms C,
            ! C_MOLEC should also be added here. (jaf, 5/12/14)
            CASE( 'ACTA' )
               BIO_SAVE(N)        = IDBACTA
               IF ( IDTACTA > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTACTA)
            CASE( 'HCOOH' )
               BIO_SAVE(N)        = IDBHCOOH
               IF ( IDTHCOOH > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTHCOOH)
            CASE( 'ISOP' )
               BIO_SAVE(N)        = IDBISOP
               IF ( IDTISOP > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTISOP)
               C_MOLEC            = 5d0
            CASE( 'RCHO' )
               BIO_SAVE(N)        = IDBRCHO
               IF ( IDTRCHO > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTRCHO)
            CASE( 'R4N2' )
               BIO_SAVE(N)        = IDBR4N2
               IF ( IDTR4N2 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTR4N2)
            CASE( 'MVK' )
               BIO_SAVE(N)        = IDBMVK
               IF ( IDTMVK > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTMVK)
            CASE( 'EOH' )
               BIO_SAVE(N)        = IDBEOH
               IF ( IDTEOH > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTEOH)
            CASE( 'HCN' )
               BIO_SAVE(N)        = IDBHCN
               IF ( IDTHCN > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTHCN)
            CASE( 'MOH' )
               BIO_SAVE(N)        = IDBMOH
               IF ( IDTMOH > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTMOH)
            CASE( 'HNO2' )
               BIO_SAVE(N)        = IDBHNO2
               IF ( IDTHNO2 > 0 )
     &            FINN_SPEC_MOLWT(N) = TRACER_MW_KG(IDTHNO2)

            ! Species that we are currently lumping with other species
            CASE( 'TMB' ) ! lumped with XYLE which is in atoms C
               LUMP_ID(N)         = IDBXYLE
               FINN_SPEC_MOLWT(N) = 12d-3
               C_MOLEC            = 9d0
            CASE( 'ETBENZ' ) ! lumped with TOLU which is in atoms C
               LUMP_ID(N)         = IDBTOLU
               FINN_SPEC_MOLWT(N) = 12d-3
               C_MOLEC            = 8d0
            CASE( 'STYR' ) ! lumped with TOLU which is in atoms C
               LUMP_ID(N)         = IDBTOLU
               FINN_SPEC_MOLWT(N) = 12d-3
               C_MOLEC            = 8d0

            ! Special case for offline species (don't have IDTs)
            CASE( 'CO2'  )
               IF ( ITS_A_CO2_SIM() ) THEN
                  BIO_SAVE(N)        = IDBCO2
                  FINN_SPEC_MOLWT(N) = 44d-3
               ENDIF
            CASE( 'CH4' )
               IF ( ITS_A_CH4_SIM() ) THEN
                  BIO_SAVE(N)        = IDBCH4
                  FINN_SPEC_MOLWT(N) = 16d-3
               ENDIF

            ! Special case for aerosols
            ! At this stage, EMFAC is [g BC or OC]/[mole CO2]
            CASE( 'BC'   )
               BIO_SAVE(N)        = IDBBC
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
               FINN_EMFAC(N,:)    = FINN_EMFAC(N,:) / 12d0
            CASE( 'OC'   )
               BIO_SAVE(N)        = IDBOC
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
               FINN_EMFAC(N,:)    = FINN_EMFAC(N,:) / 12d0
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0
         END SELECT

         ! Create emission factor for NMOC species
         ! NMOC_EMFAC is [mole NMOC] / [mole CO2]
         ! NMOC_RATIO is [mole X] / [kg NMOC]
         ! Assumed MW for NMOC is 68 g/mol
         ! Multiplying these gives [mole X] / [mole CO2].
         ! Most (not all) of these species are carried as atoms C, so
         ! we also multiply here by the number of carbon atoms/molec.
         IF ( IS_NMOC ) THEN
           DO M = 1, N_EMFAC
              FINN_EMFAC(N,M) = NMOC_EMFAC(M) * ( NMOC_RATIO(M)*68d-3 )
     &                                        * C_MOLEC
           ENDDO
           IF (FINN_SPEC_MOLWT(N) == 12d-3) FINN_SPEC_UNIT(N) = '[Tg C]'
         ENDIF

      ENDDO

      END SUBROUTINE INIT_FINN_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rearrange_biom
!
! !DESCRIPTION: Subroutine REARRANGE\_BIOM takes FINN emissions (which have 
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
!  02 Jan 2013 - J. Mao and J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      ! First loop over FINN species:
      !   Save into array w/ proper ordering for GEOS-Chem
      DO N = 1, N_SPEC

         IF ( BIO_SAVE(N) > 0 )
     &      BIOM_OUTM(:,:,BIO_SAVE(N)) = BIOM_OUT(:,:,N)

      ENDDO

      ! Second loop over FINN species:
      !   If necessary, lump into another species instead
      !   (separate loop to avoid previous loop overwriting addition of
      !   a secondary species)
      DO N = 1, N_SPEC

         IF ( LUMP_ID(N) > 0 ) BIOM_OUTM(:,:,LUMP_ID(N)) =
     &      BIOM_OUTM(:,:,LUMP_ID(N)) + BIOM_OUT(:,:,N)

      ENDDO

      END SUBROUTINE REARRANGE_BIOM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_finn_biomass
!
! !DESCRIPTION: Subroutine CLEANUP\_FINN\_BIOMASS deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_FINN_BIOMASS
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_FINN_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( FINN_SPEC_MOLWT ) ) DEALLOCATE( FINN_SPEC_MOLWT )
      IF ( ALLOCATED( FINN_SPEC_NAME  ) ) DEALLOCATE( FINN_SPEC_NAME  )
      IF ( ALLOCATED( FINN_SPEC_UNIT  ) ) DEALLOCATE( FINN_SPEC_UNIT  )
      IF ( ALLOCATED( FINN_EMFAC      ) ) DEALLOCATE( FINN_EMFAC      )
      IF ( ALLOCATED( BIOMASS_MODEL   ) ) DEALLOCATE( BIOMASS_MODEL   )
      IF ( ALLOCATED( BIO_SAVE        ) ) DEALLOCATE( BIO_SAVE        )
      IF ( ALLOCATED( LUMP_ID         ) ) DEALLOCATE( LUMP_ID         )
      
      END SUBROUTINE CLEANUP_FINN_BIOMASS
!EOC

      END MODULE FINN_BIOMASS_MOD
