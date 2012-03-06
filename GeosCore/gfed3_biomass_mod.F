!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: gfed3_biomass_mod
!     
! !DESCRIPTION: Module GFED3\_BIOMASS\_MOD contains routines and variables 
!  used in the wet scavenging of tracer in cloud updrafts, rainout, and 
!  washout. 
!\\   
!\\   
! !INTERFACE: 
!
      MODULE GFED3_BIOMASS_MOD
!
! !USES:
!
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
!
! !REMARKS:
!  Monthly/8-day/3-hr emissions of DM are read from disk and then
!  multiplied by the appropriate emission factors to produce biomass
!  burning emissions on a "generic" 1x1 grid.  The emissions are then
!  regridded to the current GEOS-Chem or GCAP grid (1x1, 2x25, or 4x5).
!                                                                             .
!  If several gfed3 options are switched on, the smaller period
!  product is used: 3-hr before 8-day before monthly.
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
!  01 Mar 2012 - R. Yantosca   - Now reference new grid_mod.F90
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
      INTEGER,          PARAMETER   :: N_SPEC  = 24
!
! PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! BIO_SAVE        : Index array to store IDBxxx values
      ! SECONDS         : Number of seconds in the current month
      ! VEG_GEN_1x1     : Array for GFED3 1x1 humid trop forest map
      ! GFED3_SPEC_NAME : Array for GFED3 biomass species names
      ! GFED3_EMFAC     : Array for user-defined emission factors
      ! BIOM_OUT        : Array for biomass emissions on model grid
      ! DOY8DAY         : Day Of the Year at start of the current
      !                   8-day period. 
      ! T3HR            : HH at start of the current 3-hr period.
      ! UPDATED         : flag to indicate if new data are read at
      !=================================================================

      ! Scalars
      LOGICAL                       :: UPDATED
      INTEGER                       :: DOY8DAY
      INTEGER                       :: T3HR
      REAL*8                        :: SECONDS

      ! Arrays
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)
      INTEGER,          ALLOCATABLE :: VEG_GEN_1x1(:,:)
      REAL*8,           ALLOCATABLE :: GFED3_EMFAC(:,:)
      REAL*8,           ALLOCATABLE :: GFED3_SPEC_MOLWT(:)
      CHARACTER(LEN=4), ALLOCATABLE :: GFED3_SPEC_NAME(:)
      CHARACTER(LEN=6), ALLOCATABLE :: GFED3_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: GFED3_BIOMASS(:,:,:)

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE LOGICAL_MOD, ONLY : LGFED3BB
      USE LOGICAL_MOD, ONLY : L8DAYBB3
      USE LOGICAL_MOD, ONLY : L3HRBB3
      USE LOGICAL_MOD, ONLY : LSYNOPBB3
      USE TIME_MOD,    ONLY : ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: DOY   ! Day of year (0-365 or 0-366 leap years)
      INTEGER, INTENT(IN) :: HH    ! Hour of day (0-23)
!
! !REMARKS:
!  The routine computes the DOY (resp. HOUR) at start of the 8-day (resp.
!  3-hour) period we are in, if the 8-day (resp. 3-hr or synoptic) GFED3
!  option is on. Result is compared to previous value to indicate if new
!  data should be read.
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
      INTEGER :: NEW_T3HR, NEW_DOY8DAY
      
      ! Reset to default
      UPDATED = .FALSE.

      ! Check if we enter a new 3hr GFED period (we assume that
      ! emissions time step is less than a day)
      IF ( L3HRBB3 .OR. LSYNOPBB3 ) THEN

         NEW_T3HR = INT( HH / 3 ) * 3

         IF ( NEW_T3HR .NE. T3HR ) THEN
            UPDATED = .TRUE.
            T3HR    = NEW_T3HR
         ENDIF         

      ! or a new 8-day GFED period
      ELSE IF ( L8DAYBB3 ) THEN

         NEW_DOY8DAY = DOY - MOD( DOY - 1, 8 )

         IF ( NEW_DOY8DAY .NE. DOY8DAY ) THEN
            UPDATED = .TRUE.
            DOY8DAY = NEW_DOY8DAY
         ENDIF

      ! or a new month (we assume that we always do emissions on
      ! 1st day 00 GMT of each month - except for the month the
      ! run starts, for which it is not required)
      ELSE IF ( LGFED3BB ) THEN 

         IF ( ITS_A_NEW_MONTH() ) UPDATED = .TRUE.
      
      ENDIF
      
      
      END SUBROUTINE CHECK_GFED3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE BPCH2_MOD,      ONLY : GET_TAU0
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE JULDAY_MOD,     ONLY : JULDAY
      USE JULDAY_MOD,     ONLY : CALDATE
      USE LOGICAL_MOD,    ONLY : LFUTURE
      USE LOGICAL_MOD,    ONLY : L8DAYBB3
      USE LOGICAL_MOD,    ONLY : L3HRBB3
      USE LOGICAL_MOD,    ONLY : LSYNOPBB3
      USE LOGICAL_MOD,    ONLY : LGFED3BB
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE TIME_MOD,       ONLY : TIMESTAMP_STRING
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1
      USE TIME_MOD,       ONLY : GET_DAY
      USE TIME_MOD,       ONLY : GET_HOUR
      USE TIME_MOD,       ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,       ONLY : ITS_A_LEAPYEAR

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: I,    J,  N,   NF 
      INTEGER            :: YYYY, MM, MM1, YYYY1
      INTEGER            :: YYYYMMDD, HHMMSS
      REAL*8             :: GFED3_EMFACX
      REAL*4             :: ARRAY_1x1(I1x1,J1x1-1)
      REAL*4             :: DM_GEN_1x1(I1x1,J1x1-1,6)
      REAL*8             :: BIOM_GEN_1x1(I1x1,J1x1-1,N_SPEC)
      REAL*8             :: BIOM_GEOS_1x1(I1x1,J1x1,N_SPEC)
      REAL*8             :: TAU0, TAU1, JD8DAY
      REAL*4             :: TMP
      CHARACTER(LEN=255) :: FILENAME1
      CHARACTER(LEN=255) :: FILENAME2
      CHARACTER(LEN=255) :: FILENAME3
      CHARACTER(LEN=255) :: FILENAME4
      CHARACTER(LEN=255) :: FILENAME5
      CHARACTER(LEN=255) :: FILENAME6
      CHARACTER(LEN=16 ) :: TIME_STR
      INTEGER            :: DD, HH, DOY

      !=================================================================
      ! GFED3_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_GFED3_BIOMASS
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
      
      IF ( UPDATED ) THEN
         GFED3_BIOMASS  = 0D0
      ELSE
         CALL REARRANGE_BIOM(GFED3_BIOMASS,BIOM_OUT)
         RETURN
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
      IF ( LGFED3BB ) THEN
         
         CALL GFED3_AVAILABLE( YYYY, 1997, 2010 )

         WRITE( 6, 410 ) YYYY, MM
 410     FORMAT( 'for year and month: ', i4, '/', i2.2, / )

         ! Create YYYYMMDD integer value
         YYYYMMDD = YYYY*10000 + MM*100 + 01

      ENDIF

     
      !=================================================================
      ! Filename, TAU0 and number of seconds
      !=================================================================
      
      ! for monthly data
      !-------------------------------
      IF ( LGFED3BB ) THEN 
      
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
         FILENAME1 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_AGW_YYYYMM.generic.1x1'
         FILENAME2 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_DEF_YYYYMM.generic.1x1'
         FILENAME3 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_FOR_YYYYMM.generic.1x1'
         FILENAME4 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_PET_YYYYMM.generic.1x1'
         FILENAME5 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_SAV_YYYYMM.generic.1x1'
         FILENAME6 = TRIM( DATA_DIR_1x1 ) //
     &       'GFED3_201110/YYYY/GFED3_DM_WDL_YYYYMM.generic.1x1'

      ENDIF
      
      !=================================================================
      ! Read GFED3 DM emissions [g/m2/month]
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
      CALL READ_BPCH2( FILENAME1, 'GFED3-BB',   91, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. )
 
      DM_GEN_1x1(:,:,1)=ARRAY_1x1(:,:)

      CALL READ_BPCH2( FILENAME2, 'GFED3-BB',   92, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. ) 

      DM_GEN_1x1(:,:,2)=ARRAY_1x1(:,:)

      CALL READ_BPCH2( FILENAME3, 'GFED3-BB',   93, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. ) 

      DM_GEN_1x1(:,:,3)=ARRAY_1x1(:,:)

      CALL READ_BPCH2( FILENAME4, 'GFED3-BB',   94, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. ) 

      DM_GEN_1x1(:,:,4)=ARRAY_1x1(:,:)

      CALL READ_BPCH2( FILENAME5, 'GFED3-BB',   95, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. ) 

      DM_GEN_1x1(:,:,5)=ARRAY_1x1(:,:)

      CALL READ_BPCH2( FILENAME6, 'GFED3-BB',   96, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         ARRAY_1x1,  QUIET=.TRUE. ) 

      DM_GEN_1x1(:,:,6)=ARRAY_1x1(:,:)

      !=================================================================
      ! Convert [g DM/m2/month] to [kg DM/cm2/month]
      !
      ! Unit Conversions:
      ! (1) g    to kg    --> Divide by 1000
      ! (2) 1/m2 to 1/cm2 --> Divide by 10000
      !=================================================================

      ! Loop over GENERIC 1x1 GRID
      DO J = 1, J1x1-1
      DO I = 1, I1x1
      DO NF = 1, N_EMFAC

         ! Set negatives to zero 
         DM_GEN_1x1(I,J,NF) = MAX( DM_GEN_1x1(I,J,NF), 0e0 )

         ! Convert [g DM/m2/month] to [kg DM/cm2/month]
         DM_GEN_1x1(I,J,NF) = DM_GEN_1x1(I,J,NF) * 1d-3 * 1d-4 

      ENDDO
      ENDDO
      ENDDO


      !=================================================================
      ! Calculate biomass species emissions on 1x1 emissions grid
      !
      ! Emission factors convert from [g DM/m2/month] to either
      ! [molec/cm2/month] or [atoms C/cm2/month]
      !
      ! Units:
      !  [  molec/cm2/month] : NOx,  CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/month] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8, 
      !                        C2H6, BC,   OC
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, NF, GFED3_EMFACX )
         DO J = 1, J1x1-1
         DO I = 1, I1x1
            BIOM_GEN_1x1(I,J,N) = 0.0
            DO NF = 1, N_EMFAC
               GFED3_EMFACX=GFED3_EMFAC(N,NF)

               ! Use woodland emission factors for 'deforestation' outside 
               ! humid tropical forest
               IF(NF.EQ.2.AND.VEG_GEN_1x1(I,J).EQ.0)
     &            GFED3_EMFACX=GFED3_EMFAC(N,6)
               BIOM_GEN_1x1(I,J,N) =  BIOM_GEN_1x1(I,J,N) +
     &                                DM_GEN_1x1(I,J,NF) * 
     &                                GFED3_EMFACX

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Regrid each species from GENERIC 1x1 GRID to GEOS-Chem 1x1 GRID
         CALL DO_REGRID_G2G_1x1( 'molec/cm2',
     &                            BIOM_GEN_1x1(:,:,N), 
     &                            BIOM_GEOS_1x1(:,:,N) )
      ENDDO

      ! Regrid from GEOS 1x1 grid to current grid.  (The unit 'molec/cm2' 
      ! is just used to denote that the quantity is per unit area.)
      CALL DO_REGRID_1x1( N_SPEC,       'molec/cm2', 
     &                    BIOM_GEOS_1x1, GFED3_BIOMASS ) 

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL GFED3_SCALE_FUTURE( GFED3_BIOMASS )
      ENDIF

      ! Print totals in Tg/month
      CALL GFED3_TOTAL_Tg( THIS_YYYY, THIS_MM )

      ! Convert from [molec/cm2/month], [molec/cm2/8day] or
      ! [molec/cm2/3hr] to [molec/cm2/s]
      GFED3_BIOMASS = GFED3_BIOMASS / SECONDS

      ! Rearrange the species to the same order as in the IDBxxx (fp, 6/09)
      ! GFED3_BIOMASS is indexed as GFED3
      ! BIOM_OUT      is indexed as IDBs
      CALL REARRANGE_BIOM( GFED3_BIOMASS, BIOM_OUT )

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE GFED3_COMPUTE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE TRACERID_MOD,           ONLY : IDBNOx,  IDBCO,   IDBSO2 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_CO2
      INTEGER :: I, J, N
      
      !=================================================================
      ! GFED3_SCALE_FUTURE begins here!
      !=================================================================

      ! Test if it's a CO2 simulation outside of the loop
      ITS_CO2 = ITS_A_CO2_SIM()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )

      ! Loop over species and grid boxes
      DO N = 1, N_SPEC
      DO J = 1, JJPAR
      DO I = 1, IIPAR 

         ! Scale each species to IPCC future scenario
         IF ( BIO_SAVE(N) == IDBNOx ) THEN

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

         ELSE IF ( ITS_CO2 ) THEN

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gfed3_total_Tg
!
! !DESCRIPTION: Subroutine GFED3\_TOTAL\_Tg prints the amount of biomass 
!  burning emissions that are emitted each month/8-day/3-hr in Tg or Tg C. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GFED3_TOTAL_Tg( YYYY, MM )
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD              ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: YYYY   ! Current year
      INTEGER, INTENT(IN) :: MM     ! Current month
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

            TOTAL = TOTAL + ( GFED3_BIOMASS(I,J,N) * CONV )
         ENDDO
         ENDDO
     
         ! Write totals
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', f9.4, 1x, a6 )
      ENDDO

      END SUBROUTINE GFED3_TOTAL_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD,      ONLY : IOERROR
      USE FILE_MOD,      ONLY : IU_FILE
      USE LOGICAL_MOD,   ONLY : LDICARB
      USE TRACERID_MOD,  ONLY : IDBNOx,  IDBCO,   IDBALK4
      USE TRACERID_MOD,  ONLY : IDBACET, IDBMEK,  IDBALD2
      USE TRACERID_MOD,  ONLY : IDBPRPE, IDBC3H8, IDBCH2O
      USE TRACERID_MOD,  ONLY : IDBC2H6, IDBBC,   IDBOC
      USE TRACERID_MOD,  ONLY : IDBSO2,  IDBNH3,  IDBCO2
      USE TRACERID_MOD,  ONLY : IDBGLYX, IDBMGLY, IDBBENZ
      USE TRACERID_MOD,  ONLY : IDBTOLU, IDBXYLE, IDBC2H4
      USE TRACERID_MOD,  ONLY : IDBC2H2, IDBGLYC, IDBHAC

      USE CMN_SIZE_MOD        ! Size parameters
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
      INTEGER            :: AS, IOS, M, N, NDUM
      REAL*4             :: ARRAY(I1x1,J1x1-1,1)
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! INIT_GFED3_BIOMASS begins here!
      !=================================================================

      ! Allocate array to hold emissions
      ALLOCATE( GFED3_BIOMASS( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_BIOMASS' )
      GFED3_BIOMASS = 0d0

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

      ! Allocate array for species molecular weight
      ALLOCATE( GFED3_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED3_SPEC_UNIT' )
      GFED3_SPEC_UNIT = ''

      ! Allocate array for vegetation map
      ALLOCATE( VEG_GEN_1x1( I1x1, J1x1-1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VEG_GEN_1x1' )

      !IDBs are now the same as the ones in TRACERID AND BIOMASS_MOD
      !BIOSAVE INDEX IS THE LOCATION OF THE EMISSION IN THE GFED FILE
      !(fp)
      ALLOCATE( BIO_SAVE( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIO_SAVE' )
      BIO_SAVE = 0


      ! Set default values for module variables
      T3HR    = -1
      DOY8DAY = -1

      !=================================================================
      ! Read emission factors (which convert from kg DM to 
      ! either [molec species] or [atoms C]) from bpch file
      !=================================================================
     
      ! File name
      FILENAME = TRIM( DATA_DIR_1x1) // 
     &           'GFED3_201110/GFED3_emission_factors.txt'

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
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'GFED3_201110/GFED3_vegmap.generic.1x1'

      ! Read GFED3 veg map
      CALL READ_BPCH2( FILENAME, 'LANDMAP',  1,
     &                 0d0,       I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to INTEGER
      VEG_GEN_1x1(:,:) = ARRAY(:,:,1)


      
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
               BIO_SAVE(N)         = IDBNOX
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
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0

              WRITE(*,*) 'NAME',TRIM( GFED3_SPEC_NAME(N) )
         END SELECT
      ENDDO

      END SUBROUTINE INIT_GFED3_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
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

      ! Loop over GFED2 species
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
      IF ( ALLOCATED( VEG_GEN_1x1      ) ) DEALLOCATE( VEG_GEN_1x1     )
      IF ( ALLOCATED( GFED3_BIOMASS    ) ) DEALLOCATE( GFED3_BIOMASS   )
      
      END SUBROUTINE CLEANUP_GFED3_BIOMASS
!EOC
      END MODULE GFED3_BIOMASS_MOD
