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
#     include "define.h"
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! MODULE PARAMETERS
      !
      ! N_SPEC  : Number of species
      !=================================================================
      ! Consider increasing this to 25 (skim, 7/27/13)
      INTEGER,          PARAMETER   :: N_SPEC  = 23

!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! Scalars
      !
      ! UPDATED         : flag to indicate if FINN emissions are updated
      ! UPDATED_MON     : flag to indicate if new month
      ! UPDATED_DAY     : flag to indicate if new day
      !                   - only set to true if daily emissions are used
      ! SECONDS         : Number of seconds in the current month
      !
      ! Arrays
      !
      ! FINN_SPEC_NAME  : Array for FINN biomass species names
      ! FINN_SPEC_MOLWT : Array for FINN biomass species molecular wts
      ! FINN_SPEC_UNIT  : Array for FINN biomass species emissions units
      ! BIOMASS_MODEL   : Array for FINN species emissions on model grid
      ! BIO_SAVE        : Index array to store IDBxxx values
      !=================================================================

      ! Scalars
      LOGICAL                       :: UPDATED
      LOGICAL                       :: UPDATED_MON
      LOGICAL                       :: UPDATED_DAY
      REAL*8                        :: SECONDS

      ! Arrays
      CHARACTER(LEN=4), ALLOCATABLE :: FINN_SPEC_NAME(:)
      REAL*8,           ALLOCATABLE :: FINN_SPEC_MOLWT(:)
      CHARACTER(LEN=6), ALLOCATABLE :: FINN_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: BIOMASS_MODEL(:,:,:)
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)

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
!     Hardcode this for the moment (skim, 7/27/13)
!     USE LOGICAL_MOD, ONLY : L1DAYBBF
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
      UPDATED_MON = .FALSE.
      UPDATED_DAY = .FALSE.

      ! Check if it is a new month
      IF ( ITS_A_NEW_MONTH() ) THEN
            UPDATED     = .TRUE.
            UPDATED_MON = .TRUE.
      ENDIF
      
      ! Check if it is a new day
      ! Use the daily emissions, hardcoded (skim, 7/27/13)
      !IF ( L1DAYBBF ) THEN
         IF ( ITS_A_NEW_DAY() ) THEN
            UPDATED     = .TRUE.
            UPDATED_DAY = .TRUE.
         ENDIF
      !ENDIF

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

      !!!!! Kludge - also have 2012 emissions at the moment (jaf)
      IF ( YYYY == 2012 ) RETURN

      !!!!! Kludge - also have 2013 emissions for the SEAC4RS mission
      ! (skim, 7/27/13)
      IF ( YYYY == 2013 ) RETURN

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
!      HARDCODE IN FOR NOW (skim, 7/27/13)
!      USE LOGICAL_MOD,     ONLY : L1DAYBBF
!      USE LOGICAL_MOD,     ONLY : LFINNBB
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
      INTEGER             :: N
      INTEGER             :: YYYY, MM, MM1, YYYY1
      INTEGER             :: YYYYMMDD, YYYYMMDD1, HHMMSS
      REAL*8              :: ARRAY_FINN(IIPAR, JJPAR)
      REAL*8              :: TAU0, TAU1
      REAL*8              :: JD,   JD1
      INTEGER             :: DD,   DD1,  DOY
      REAL*4              :: TMP
      CHARACTER(LEN=255)  :: FILENAME,  FINN_DIR
      CHARACTER(LEN=7  )  :: TIME_STR
      CHARACTER(LEN=4  )  :: SPEC_NAME

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
      
      ! DAILY data
      !-------------------------------
      ! Hardcode for the SEAC4RS campaign
      IF ( .TRUE. ) THEN
      !IF ( L1DAYBBF ) THEN

         ! Check availability
         CALL FINN_AVAILABLE( YYYY, 2006, 2007 )

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
         ! FINN_DIR = 'FINN_201301/daily/'
         FINN_DIR = '/home/skim/FINN/FINN_files/daily/'
         TIME_STR = 'daily'

      ! MONTHLY data
      !-------------------------------
      ELSE 

         CALL FINN_AVAILABLE( YYYY, 2006, 2007 )

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
         FINN_DIR = 'FINN_201301/monthly/'
         TIME_STR = 'monthly'

      ENDIF

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

      ! Loop over biomass species 
      DO N = 1, N_SPEC

         ! Initialize temporary array
         ARRAY_FINN = 0d0

         ! Set filename for this species
         SPEC_NAME = FINN_SPEC_NAME(N)
!         FILENAME = TRIM( DATA_DIR_NATIVE ) // TRIM( FINN_DIR ) //
!     &              TRIM( TIME_STR ) // '_' // TRIM( SPEC_NAME ) //
!     &              '_YYYY_05x05.nc'
         ! This will probably need to be updated once we get the
         ! filenames (skim, 7/27/13)
         FILENAME = TRIM( FINN_DIR ) //
     &              TRIM( TIME_STR ) // '_' // TRIM( SPEC_NAME ) //
     &              '_YYYY_05x05.nc'
         CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

         ! Read FINN emissions for this species in [molec/cm2/s]
         CALL READ_FINN( FILENAME, 'fire', ARRAY_FINN, YYYY, MM, DD )

         ! Add to biomass array
         BIOMASS_MODEL(:,:,N) = ARRAY_FINN

      ENDDO

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL FINN_SCALE_FUTURE( BIOMASS_MODEL )
      ENDIF

      ! Print totals in Tg/time period
      IF ( UPDATED_DAY ) THEN
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
! !DESCRIPTION: Subroutine READ\_FINN reads the FINN BB emissions at 0.5x0.5
!  degree resolution and regrids them to the current model resolution
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_FINN( FILENAME, PARNAME, ARR, YYYY, MM, DD )
!
! !USES:
!
      USE DIRECTORY_MOD,  ONLY: DATA_DIR_NATIVE => DATA_DIR_1x1
      USE NCDF_MOD,       ONLY: NC_READ
      USE REGRID_A2A_MOD, ONLY: DO_REGRID_A2A

      USE CMN_SIZE_MOD          ! Size parameters

!
! !INPUT PARAMETERS:
!
      ! File and parameter to read
      CHARACTER(LEN=255),  INTENT(IN) :: FILENAME
      CHARACTER(LEN=4),    INTENT(IN) :: PARNAME
      INTEGER,             INTENT(IN) :: YYYY
      INTEGER,             INTENT(IN) :: MM
      INTEGER,             INTENT(IN) :: DD
! 
! !OUTPUT PARAMETERS:
!
      ! Emissions output
      REAL*8,            INTENT(OUT) :: ARR(IIPAR,JJPAR)
! 
! !REVISION HISTORY:
!  02 Jan 2013 - J. Fisher - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, PARAMETER    :: I_NATIVE = 720 ! 0.5 x 0.5 resolution
      INTEGER, PARAMETER    :: J_NATIVE = 360 ! 0.5 x 0.5 resolution
      REAL*8                :: GENERIC_NATIVE(I_NATIVE,J_NATIVE,1)
      CHARACTER(LEN=255)    :: LLFILENAME
      REAL*8                :: INGRID(I_NATIVE,J_NATIVE)
      ! for ncdf-reading
      REAL*4, POINTER       :: NCARR(:,:,:) => NULL()

      !=================================================================
      ! READ_FINN begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( PARNAME ), TRIM( FILENAME )
 100  FORMAT( '     - READ_FINN : Reading ', a , ' from ', a ) 

      ARR   = 0D0 

      CALL NC_READ ( NC_PATH = TRIM(FILENAME) ,
     &               PARA    = TRIM(PARNAME)  ,
     &               ARRAY   = NCARR          ,
     &               YEAR    = YYYY           ,
     &               MONTH   = MM             ,
     &               DAY     = DD              )

      ! Cast to REAL*8 before regridding
      GENERIC_NATIVE(:,:,1) = NCARR(:,:,1)

      ! Deallocate ncdf-array
      IF ( ASSOCIATED ( NCARR ) ) DEALLOCATE ( NCARR )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_NATIVE) //
     &             'MAP_A2A_Regrid_201203/' //
     &             'MAP_A2A_latlon_generic05x05.nc'

      ! Regrid to current model resolution
      INGRID = GENERIC_NATIVE(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I_NATIVE, J_NATIVE,
     &                    INGRID,     ARR,      IS_MASS=0,
     &                    netCDF=.TRUE.                    )

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
!  02 Jan 2013 - J. Fisher - Initial version, based on GFED3
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
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', e12.5, 1x, a6 )
      ENDDO

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
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE FILE_MOD,        ONLY : IOERROR
      USE LOGICAL_MOD,     ONLY : LDICARB
!     Hardcode for the SEAC4RS campaign (skim, 7/27/13) 
!     USE LOGICAL_MOD,     ONLY : L1DAYBBF
      USE TRACERID_MOD,    ONLY : IDBNO,   IDBCO,   IDBALK4
      USE TRACERID_MOD,    ONLY : IDBACET, IDBMEK,  IDBALD2
      USE TRACERID_MOD,    ONLY : IDBPRPE, IDBC3H8, IDBCH2O
      USE TRACERID_MOD,    ONLY : IDBC2H6, IDBBC,   IDBOC
      USE TRACERID_MOD,    ONLY : IDBSO2,  IDBNH3,  IDBCO2
      USE TRACERID_MOD,    ONLY : IDBGLYX, IDBMGLY, IDBBENZ
      USE TRACERID_MOD,    ONLY : IDBTOLU, IDBXYLE, IDBC2H4
      USE TRACERID_MOD,    ONLY : IDBC2H2, IDBGLYC, IDBHAC
      USE TRACERID_MOD,    ONLY : IDBCH4
      USE CMN_SIZE_MOD          ! Size parameters
! 
! !REVISION HISTORY: 
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IU_FILE
      INTEGER            :: AS, IOS, M, N, NDUM
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! INIT_FINN_BIOMASS begins here!
      !=================================================================

      ! Allocate array to hold FINN species emissions on model grid
      ALLOCATE( BIOMASS_MODEL( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_MODEL' )
      BIOMASS_MODEL = 0d0

      ! Allocate array for species molecular weight
      ALLOCATE( FINN_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_MOLWT' )
      FINN_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( FINN_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_NAME' )
      FINN_SPEC_NAME = ''

      ! Allocate array for FINN biomass buning species mass units
      ALLOCATE( FINN_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FINN_SPEC_UNIT' )
      FINN_SPEC_UNIT = ''

      !IDBs are now the same as the ones in TRACERID AND BIOMASS_MOD
      !BIOSAVE INDEX IS THE LOCATION OF THE EMISSION IN THE FINN FILE
      !(fp)
      ALLOCATE( BIO_SAVE( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIO_SAVE' )
      BIO_SAVE = 0

      !=================================================================
      ! Define FINN species names
      !=================================================================
      FINN_SPEC_NAME(1)  = 'NOx'
      FINN_SPEC_NAME(2)  = 'CO'
      FINN_SPEC_NAME(3)  = 'ALK4'
      FINN_SPEC_NAME(4)  = 'ACET'
      FINN_SPEC_NAME(5)  = 'MEK'
      FINN_SPEC_NAME(6)  = 'ALD2'
      FINN_SPEC_NAME(7)  = 'PRPE'
      FINN_SPEC_NAME(8)  = 'C3H8'
      FINN_SPEC_NAME(9)  = 'CH2O'
      FINN_SPEC_NAME(10) = 'C2H6'
      FINN_SPEC_NAME(11) = 'BC'
      FINN_SPEC_NAME(12) = 'OC'
      FINN_SPEC_NAME(13) = 'SO2'
      FINN_SPEC_NAME(14) = 'NH3'
      FINN_SPEC_NAME(15) = 'CO2'
      FINN_SPEC_NAME(16) = 'TOLU'
      FINN_SPEC_NAME(17) = 'BENZ'
      FINN_SPEC_NAME(18) = 'C2H2'
      FINN_SPEC_NAME(19) = 'C2H4'
      FINN_SPEC_NAME(20) = 'MGLY'
      FINN_SPEC_NAME(21) = 'GLYC'
      FINN_SPEC_NAME(22) = 'HAC'
      FINN_SPEC_NAME(23) = 'CH4'

      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the FINN biomass species
      !=================================================================
      ! Save correspondance between FINN species order (N) and 
      ! species order of the simulation (IDBxxxs).(ccc, 2/4/10)
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( TRIM( FINN_SPEC_NAME(N) ) ) 
            CASE( 'NOx'  )
               BIO_SAVE(N)         = IDBNO
               FINN_SPEC_MOLWT(N) = 14d-3
               FINN_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO'   )
               BIO_SAVE(N)         = IDBCO
               FINN_SPEC_MOLWT(N) = 28d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ALK4' )
               BIO_SAVE(N)         = IDBALK4
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ACET' )
               BIO_SAVE(N)         = IDBACET
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MEK'  )
               BIO_SAVE(N)         = IDBMEK
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ALD2' )
               BIO_SAVE(N)         = IDBALD2
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'PRPE' )
               BIO_SAVE(N)         = IDBPRPE
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C3H8' )
               BIO_SAVE(N)         = IDBC3H8
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'CH2O' )
               BIO_SAVE(N)         = IDBCH2O
               FINN_SPEC_MOLWT(N) = 30d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'C2H6' )
               BIO_SAVE(N)         = IDBC2H6
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'SO2'  )
               BIO_SAVE(N)         = IDBSO2
               FINN_SPEC_MOLWT(N) = 64d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'NH3'  )
               BIO_SAVE(N)         = IDBNH3
               FINN_SPEC_MOLWT(N) = 17d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BC'   )
               BIO_SAVE(N)         = IDBBC
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'OC'   )
               BIO_SAVE(N)         = IDBOC
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MGLY' )
               BIO_SAVE(N)         = IDBMGLY
               FINN_SPEC_MOLWT(N) = 72d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BENZ' )
               BIO_SAVE(N)         = IDBBENZ
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'TOLU' )
               BIO_SAVE(N)         = IDBTOLU
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H4' )
               BIO_SAVE(N)         = IDBC2H4
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H2' )
               BIO_SAVE(N)         = IDBC2H2
               FINN_SPEC_MOLWT(N) = 12d-3
               FINN_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYC' )
               BIO_SAVE(N)         = IDBGLYC
               FINN_SPEC_MOLWT(N) = 60d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'HAC' )
               BIO_SAVE(N)         = IDBHAC
               FINN_SPEC_MOLWT(N) = 74d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CO2'  )
               BIO_SAVE(N)         = IDBCO2
               FINN_SPEC_MOLWT(N) = 44d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CH4' )
               BIO_SAVE(N)         = IDBCH4
               FINN_SPEC_MOLWT(N) = 16d-3
               FINN_SPEC_UNIT(N)  = '[Tg  ]'
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0

            WRITE(*,*) 'NAME',TRIM( FINN_SPEC_NAME(N) )
         END SELECT
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

      ! Loop over FINN species
      DO N = 1, N_SPEC

         ! Save into array w/ proper ordering for GEOS-Chem
         IF ( BIO_SAVE(N) .GT. 0 ) THEN
            BIOM_OUTM(:,:,BIO_SAVE(N)) = BIOM_OUT(:,:,N)
         ENDIF

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
      IF ( ALLOCATED( BIOMASS_MODEL   ) ) DEALLOCATE( BIOMASS_MODEL   )
      
      END SUBROUTINE CLEANUP_FINN_BIOMASS
!EOC

      END MODULE FINN_BIOMASS_MOD
