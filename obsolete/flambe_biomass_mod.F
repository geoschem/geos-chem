!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: flambe_biomass_mod
!     
! !DESCRIPTION: Module FLAMBE\_BIOMASS\_MOD containts routines to
! compute biomass burning emissions from the FLAMBE inventory for the
! SEAC4RS campaign. This code is based off of gfed3_biomass_mod.F90
!\\   
!\\   
! !INTERFACE: 
!
      MODULE FLAMBE_BIOMASS_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
#     include "define.h"
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: FLAMBE_COMPUTE_BIOMASS
      PUBLIC  :: CLEANUP_FLAMBE_BIOMASS
!
! !REVISION HISTORY: 
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
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

      ! Assume that we are using the GFED2 species and EF's for now
      ! Also ignore CH4 for now for SEAC4RS
      ! (skim, 6/27/13)
      INTEGER,          PARAMETER   :: N_EMFAC = 3   !6
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
      !LOGICAL                       :: UPDATED_MON
      !LOGICAL                       :: UPDATED_DAY
      !LOGICAL                       :: UPDATED_3HR
      INTEGER                       :: DOY8DAY
      INTEGER                       :: T3HR
      REAL*8                        :: SECONDS
      INTEGER                       :: IIIPAR0
      INTEGER                       :: JJJPAR0

      ! Arrays
      INTEGER,          ALLOCATABLE :: VEG_GEN_1x1(:,:)
      CHARACTER(LEN=4), ALLOCATABLE :: FLAMBE_SPEC_NAME(:)
      REAL*8,           ALLOCATABLE :: FLAMBE_SPEC_MOLWT(:)
      CHARACTER(LEN=6), ALLOCATABLE :: FLAMBE_SPEC_UNIT(:)
      REAL*8,           ALLOCATABLE :: FLAMBE_EMFAC(:,:)
      !REAL*8,           ALLOCATABLE :: DM_GFED3_MON(:,:,:)
      !REAL*8,           ALLOCATABLE :: DM_GFED3_DAY(:,:,:)
      !REAL*4,           ALLOCATABLE :: FR_GFED3_3HR(:,:,:)
      !INTEGER,          ALLOCATABLE :: HUMTROP_GFED3(:,:)
      REAL*8,           ALLOCATABLE :: BIOMASS_MODEL(:,:,:)
      INTEGER,          ALLOCATABLE :: BIO_SAVE(:)
      !REAL*8,           ALLOCATABLE :: XEDGE_GFED3(:)
      !REAL*8,           ALLOCATABLE :: YEDGE_GFED3(:)
      !REAL*8,           ALLOCATABLE :: XEDGE_MODELG(:)
      !REAL*8,           ALLOCATABLE :: YEDGE_MODELG(:)

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
! !IROUTINE: flambe_compute_biomass
!
! !DESCRIPTION: Subroutine FLAMBE\_COMPUTE\_BIOMASS computes the monthly 
!  FLAMBE biomass burning emissions for a given year, month, day, and
!  hour
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FLAMBE_COMPUTE_BIOMASS( BIOM_OUT, 
     &                                   THIS_YYYY, THIS_MM, 
     &                                   THIS_DD,   THIS_HH  )
!
! !USES:
!
      USE BPCH2_MOD
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,        ONLY : FILE_EXISTS
      USE JULDAY_MOD,      ONLY : JULDAY
      USE JULDAY_MOD,      ONLY : CALDATE
      USE LOGICAL_MOD,     ONLY : LFUTURE
      !USE LOGICAL_MOD,     ONLY : LDAYBB3
      !USE LOGICAL_MOD,     ONLY : L3HRBB3
      !USE LOGICAL_MOD,     ONLY : LGFED3BB
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
      USE REGRID_A2A_MOD,  ONLY : DO_REGRID_A2A
      ! Needed for SMOKE emissions (skim, 6/27/13)
      USE TRACERID_MOD,    ONLY : IDBBC,   IDBOC

      USE CMN_SIZE_MOD                                       ! Size parameters
      USE CMN_GCTM_MOD
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THIS_YYYY                      ! Current year   
      INTEGER, INTENT(IN)  :: THIS_MM                        ! Current month
      INTEGER, INTENT(IN)  :: THIS_DD                        ! Current day
      INTEGER, INTENT(IN)  :: THIS_HH                        ! Current hour
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: BIOM_OUT(IIPAR,JJPAR,NBIOMAX)  ! BB emissions
                                                             ! [molec/cm2/s]
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, J, N, N_VEG, NF, IT3, L, L1, L2
      INTEGER             :: AS
      INTEGER             :: II,   JJ
      INTEGER             :: I0,   J0
      INTEGER             :: YYYY, MM, MM1, YYYY1
      INTEGER             :: YYYYMMDD, HHMMSS
      INTEGER             :: IU_ADJ
      REAL*8              :: GFED3_EMFACX
      REAL*4              :: ARRAY_GFED3(IGFED3, JGFED3, 1)
      REAL*4              :: FR_GFED3_DAY(IGFED3, JGFED3)
      REAL*8              :: DM_GFED3(IGFED3, JGFED3, N_EMFAC)
      REAL*8              :: BIOMASS_GFED3(IGFED3, JGFED3, N_SPEC)
      REAL*8              :: TAU0, TAU1
      REAL*4              :: TMP
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: LLFILENAME
      CHARACTER(LEN=16 )  :: TIME_STR, STAMP
      CHARACTER(LEN=80 )  :: TITLE
      CHARACTER(LEN=30 )  :: OUTPUT_FILE
      CHARACTER(LEN=40 )  :: CATEGORY
      CHARACTER(LEN=40 )  :: UNIT
      CHARACTER(LEN=20 )  :: MODELNAME
      CHARACTER(LEN=40 )  :: RESERVED = ''
      CHARACTER(LEN=255)  :: FILENAME_OUT
      INTEGER             :: DD, HH, DOY
      INTEGER             :: IT3HR, IT3H
      REAL*4              :: DEG2RAD, LONRES, LATRES
      REAL*4, ALLOCATABLE :: temp(:,:)
      REAL*8, ALLOCATABLE :: bq1(:,:)
      REAL*8, ALLOCATABLE :: bq2(:,:)

      REAL*4              :: ARRAY(I1x1,J1x1-1)
      REAL*4              :: DM_GEN_1x1(I1x1,J1x1-1)
      REAL*8              :: SM_GEN_1x1(I1x1,J1x1-1)
      REAL*8              :: BIOM_GEN_1x1(I1x1,J1x1-1,N_SPEC)
      REAL*8              :: BIOM_GEOS_1x1(I1x1,J1x1,N_SPEC)
      REAL*8              :: JD8DAY
      REAL*8              :: Z, Z_MIN, Z_MED, Z_MAX, N_LEV

      !=================================================================
      ! FLAMBE_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_FLAMBE_BIOMASS
         FIRST = .FALSE.
      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM
      DD   = THIS_DD
      HH   = THIS_HH
      DOY  = GET_DAY_OF_YEAR()
      
      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     & 'F L A M B E   B I O M A S S   B U R N I N G   E M I S S I O N S'

      ! Zero local arrays
      DM_GEN_1x1    = 0d0
      SM_GEN_1x1    = 0d0
      BIOM_GEN_1x1  = 0d0
      BIOM_GEOS_1x1 = 0d0

      !=================================================================
      ! Filename, TAU0 and number of seconds
      !=================================================================
      
      TAU0     = GET_TAU0( MM, DD, YYYY, HH   )
      TAU1     = GET_TAU0( MM, DD, YYYY, HH+1 )

      ! Number of seconds in this hour 
      SECONDS  = ( TAU1 - TAU0 ) * 3600d0

      ! KLUDGE FOR SEAC4RS
      ! Hardwire the path for now (skim, 6/21/13)
      FILENAME = '/home/skim/Code_Development/FLAMBE/' //
     &           'arctas_hourly/YYYY/YYYYMM/'          // 
     &           'flambe_YYYYMMDDhhmm.generic.1x1'

      !=================================================================
      ! Read FLAMBE files C emissions [g/m2]
      !=================================================================
      
      ! Create YYYYMMDD and HHMMSS for file name timestamp
      YYYYMMDD = YYYY*10000 + MM*100 + DD
      HHMMSS   = HH*10000   + 0000

      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Check to make sure the file exists
      IF ( FILE_EXISTS( FILENAME ) ) THEN

         !--------------------------------------------
         ! Read FLAMBE C emissions [g C/cm2]
         !--------------------------------------------
         CALL READ_BPCH2( FILENAME, 'GFED2-BB', 99, 
     &                    TAU0,      I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Copy data into DM_GEN_1x1
         DM_GEN_1x1 = ARRAY

         !--------------------------------------------
         ! Read FLAMBE smoke emissions [g C/cm2]
         !--------------------------------------------
         CALL READ_BPCH2( FILENAME, 'OC-BIOB', 33, 
     &                    TAU0,      I1x1,     J1x1-1,     
     &                    1,         ARRAY,    QUIET=.TRUE. ) 

         ! Copy data into SM_GEN_1x1
         SM_GEN_1x1 = ARRAY

      ELSE

         PRINT*, 'FLAMBE emissions file '// TRIM( FILENAME ) //
     &           ' not found!'
         CALL FLUSH(6)
         CALL GEOS_CHEM_STOP

      ENDIF

      !=================================================================
      ! Convert C [g/cm2] to dry matter burned [kg/cm2]
      !
      ! Unit Conversions:
      ! (1) C    to DM    --> Divide by 0.45  
      ! (2) g    to kg    --> Divide by 1000  
      !=================================================================

      ! Loop over GENERIC 1x1 GRID
      DO J = 1, J1x1-1
      DO I = 1, I1x1

         ! Set negatives to zero
         DM_GEN_1x1(I,J) = MAX( DM_GEN_1x1(I,J), 0e0 )

         ! Convert [g C/cm2] to [kg DM/cm2]
         DM_GEN_1x1(I,J) = DM_GEN_1x1(I,J) / ( 0.45d0 * 1d3 )

         ! Convert [g smoke/cm2] to [molec/cm2] (assumes molwt of C)
         SM_GEN_1x1(I,J) = SM_GEN_1x1(I,J) * 6.0225d+23 / 12d0

      ENDDO
      ENDDO

      !=================================================================
      ! Calculate biomass species emissions on FLAMBE grid
      ! and regrid to model grid
      !
      ! Emission factors convert from [kg DM/cm2/hour] to either
      ! [molec/cm2/hour] or [atoms C/cm2/hour]
      !
      ! Units:
      !  [  molec/cm2/hour] : NOx,  CO,   CH2O, SO2,  NH3,  CO2
      !  [atoms C/cm2/hour] : ALK4, ACET, MEK,  ALD2, PRPE, C3H8,
      !                        C2H6, BC,   OC
      !
      ! NOTE: For smoke we use the ratio of emission factors to
      ! partition between BC and OC. This ratio is unitless, so that
      ! is why we had to convert SM_GEN_1x1 to molec/cm2 above
      ! (bmy, phs, 6/11/08)
      !=================================================================

      ! Loop over biomass species 
      DO N = 1, N_SPEC

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N_VEG )
         DO J = 1, J1x1-1
         DO I = 1, I1x1
 
            ! Vegetation type index
            N_VEG = VEG_GEN_1x1(I,J)

            IF ( BIO_SAVE(N) == IDBBC ) THEN 

               !-------------------------------------------------
               ! BC: derived from smoke aerosol
               ! NOTE: SM_GEN_1x1 is already in [molec/cm2]
               !-------------------------------------------------

               IF ( N_VEG == 0 ) THEN

                  ! No emissions over the ocean
                  BIOM_GEN_1x1(I,J,N) = 0d0

               ELSE

                  ! Temporary kludge
                  ! N = 13 is BC, N = 14 is OC
                  ! Hardwire N+1 for OC here

                  ! BC (multiply smoke by emission ratio BC/(BC+OC))
                  BIOM_GEN_1x1(I,J,N) = 
     &             ( SM_GEN_1x1(I,J) * FLAMBE_EMFAC(N,N_VEG) ) /  
     &             ( FLAMBE_EMFAC(N,N_VEG) + 
     &               FLAMBE_EMFAC(N+1,N_VEG) )

               ENDIF

            ELSE IF ( BIO_SAVE(N) == IDBOC ) THEN 

               !-------------------------------------------------
               ! OC: derived from smoke aerosol
               ! NOTE: SM_GEN_1x1 is already in [molec/cm2]
               !-------------------------------------------------

               IF ( N_VEG == 0 ) THEN

                  ! No emissions over the ocean
                  BIOM_GEN_1x1(I,J,N) = 0d0

               ELSE

                  ! Temporary kludge
                  ! N = 13 is BC, N = 14 is OC
                  ! Hardwire N-1 for BC here

                  ! OC (multiply smoke by emission ratio OC/(BC+OC))
                  BIOM_GEN_1x1(I,J,N) = 
     &             ( SM_GEN_1x1(I,J) * FLAMBE_EMFAC(N,N_VEG) ) /  
     &             ( FLAMBE_EMFAC(N-1,N_VEG) + 
     &               FLAMBE_EMFAC(N,N_VEG) ) 

               ENDIF

            ELSE

               !-------------------------------------------------
               ! Everything except BC & OC: derived from Carbon
               ! NOTE: Em factors change [kg/cm2] to [molec/cm2]
               !-------------------------------------------------            
            
            ! Multiply DM * EMISSION FACTOR to get biomass emissions
            ! for each species on the GENERIC 1x1 GRID 
            SELECT CASE( N_VEG )             

               ! Ocean 
               CASE( 0 ) 
                  BIOM_GEN_1x1(I,J,N) = 0d0

               ! Land
               CASE( 1:3 )
                  BIOM_GEN_1x1(I,J,N) = DM_GEN_1x1(I,J) * 
     &                                  FLAMBE_EMFAC(N,N_VEG)

               ! Otherwise do nothing
               CASE DEFAULT
                 
            END SELECT

            ENDIF
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Regrid emissions from GFED3 grid to model grid
         ALLOCATE( bq1( I1x1, J1x1-1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'bq1' )
         ALLOCATE( bq2( IIPAR, JJPAR ), STAT=AS )
         !ALLOCATE( bq2( IIIPAR0, JJJPAR0 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'bq2' )

         bq1(:,:)=BIOM_GEN_1x1(:,:,N)

!         CALL MAP_A2A( I1x1,         J1x1-1,       XEDGE_FLAMBE,
!     &                 YEDGE_FLAMBE, bq1,          IIIPAR0,
!     &                 JJJPAR0,      XEDGE_MODELG, YEDGE_MODELG,
!     &                 bq2,          0,            0 )

         ! File with lat/lon edges for regridding
         LLFILENAME = TRIM( DATA_DIR_NATIVE )        //
     &                'MAP_A2A_Regrid_201203/'       //
     &                'MAP_A2A_latlon_generic1x1.nc'

         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1,
     &                       bq1,        bq2,     IS_MASS=0,
     &                       netCDF=.TRUE.                   )

! THIS SEGMENT OF CODE IS FOR TESTING PURPOSES ONLY
!         I0    = GET_XOFFSET( GLOBAL=.TRUE. )
!         J0    = GET_YOFFSET( GLOBAL=.TRUE. )

!         IF ( N .EQ. 2 ) THEN
!
!            ALLOCATE( temp( IIPAR, JJPAR ), STAT=AS )
!            IF ( AS /= 0 ) CALL ALLOC_ERR('temp')

!            LONRES       = DISIZE
!            LATRES       = DJSIZE
!            MODELNAME    = GET_MODELNAME()
!            FILENAME_OUT = 'test.bpch'
!            CATEGORY     = 'IJ-AVG-$'
!            UNIT         = 'kg'
!            temp         = bq2

!            IU_ADJ = 32767

!            CALL OPEN_BPCH2_FOR_WRITE( IU_ADJ, FILENAME_OUT, TITLE )

!            CALL BPCH2( IU_ADJ, MODELNAME, LONRES,   LATRES,
!     &                  1,      1,         CATEGORY, 1,
!     &                  UNIT,   TAU0,      TAU0,     RESERVED,
!     &                  IIPAR,  JJPAR,     1,        160,
!     &                  399,    1,         TEMP              )

!            CLOSE( IU_ADJ )

!         ENDIF

         BIOMASS_MODEL(:,:,N) = bq2

         !DO JJ=1,JJPAR
         !DO II=1,IIPAR
         !   BIOMASS_MODEL(II,JJ,N)=bq2(II+I0,JJ+J0)
         !   IF ( bq2( II+I0, JJ+J0 ) .GT. 0 ) THEN
         !      PRINT*, bq2(II+I0, JJ+J0)
         !   ENDIF
         !END DO
         !END DO

         DEALLOCATE( bq1 )
         DEALLOCATE( bq2 )

      END DO

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL FLAMBE_SCALE_FUTURE( BIOMASS_MODEL )
      ENDIF

      WRITE( 6, 410 ) YYYY, MM, DD, HH
 410     FORMAT( 'FLAMBE BIOMASS BURNING EMISSIONS FOR YYYY/MM/DD/HH: ',
     &            i4, '/', i2.2, '/', i2.2, '/', i2.2, '/' )
 998  CONTINUE
      CALL FLAMBE_TOTAL_Tg

      ! Convert from [molec/cm2/hour] to [molec/cm2/s]
      BIOMASS_MODEL = BIOMASS_MODEL / SECONDS

      ! Rearrange the species to the same order as in the IDBxxx (fp, 6/09)
      ! BIOMASS_MODEL is indexed as FLAMBE
      ! BIOM_OUT      is indexed as IDBs
      CALL REARRANGE_BIOM( BIOMASS_MODEL, BIOM_OUT )

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE FLAMBE_COMPUTE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FLAMBE_scale_future
!
! !DESCRIPTION: Subroutine FLAMBE\_SCALE\_FUTURE applies the IPCC future 
!  emissions scale factors to the FLAMBE biomass burning emisisons in order 
!  to compute the future emissions of biomass burning for NOx, CO, and VOC's.  
!\\
!\\
! !INTERFACE:

      SUBROUTINE FLAMBE_SCALE_FUTURE( BB )
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
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
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
      ! FLAMBE_SCALE_FUTURE begins here!
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

      END SUBROUTINE FLAMBE_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FLAMBE_TOTAL_Tg
!
! !DESCRIPTION: Subroutine FLAMBE\_TOTAL\_Tg prints the amount of biomass 
!  burning emissions that are emitted each hour
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FLAMBE_TOTAL_Tg
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD      ! Size parameters
!
! !REVISION HISTORY: 
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
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
      ! FLAMBE_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         ! Initialize
         NAME  = FLAMBE_SPEC_NAME(N)
         MOLWT = FLAMBE_SPEC_MOLWT(N)
         UNIT  = FLAMBE_SPEC_UNIT(N)
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

      END SUBROUTINE FLAMBE_TOTAL_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_flambe_biomass
!
! !DESCRIPTION: Subroutine INIT\_FLAMBE\_BIOMASS allocates all module arrays.  
!  It also reads the emission factors at the start of a GEOS-Chem
!  simulation. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_FLAMBE_BIOMASS
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : READ_BPCH2
      USE DIRECTORY_MOD,   ONLY : DATA_DIR_NATIVE => DATA_DIR_1x1
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE FILE_MOD,        ONLY : IOERROR
      USE LOGICAL_MOD,     ONLY : LDICARB
      USE LOGICAL_MOD,     ONLY : LDAYBB3
      USE LOGICAL_MOD,     ONLY : L3HRBB3
      USE TRACERID_MOD,    ONLY : IDBNO,   IDBCO,   IDBALK4
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
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IU_FILE
      INTEGER            :: AS, IOS, M, N, NDUM
      REAL*4             :: ARRAY(I1x1, J1x1-1, 1)
      REAL*4             :: ARRAY_LANDMAP(IGFED3,JGFED3,1)
      CHARACTER(LEN=255) :: FILENAME
      
      !=================================================================
      ! INIT_FLAMBE_BIOMASS begins here!
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

      ! Allocate array to hold GFED3 species emissions on model grid
      ALLOCATE( BIOMASS_MODEL( IIPAR, JJPAR, N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_MODEL' )
      BIOMASS_MODEL = 0d0

      ! Allocate array to hold monthly GFED3 DM burnt GFED3 grid
      ! ALLOCATE( DM_GFED3_MON( IGFED3, JGFED3, N_SPEC ), STAT=AS )
      ! IF ( AS /= 0 ) CALL ALLOC_ERR( 'DM_GFED3_MON' )
      ! DM_GFED3_MON = 0d0

      ! Allocate array for emission factors
      ALLOCATE( FLAMBE_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FLAMBE_EMFAC' )
      FLAMBE_EMFAC = 0d0
      
      ! Allocate array for species molecular weight
      ALLOCATE( FLAMBE_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FLAMBE_SPEC_MOLWT' )
      FLAMBE_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( FLAMBE_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FLAMBE_SPEC_NAME' )
      FLAMBE_SPEC_NAME = ''

      ! Allocate array for FLAMBE biomass buning species mass units
      ALLOCATE( FLAMBE_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FLAMBE_SPEC_UNIT' )
      FLAMBE_SPEC_UNIT = ''

      ! Allocate array for vegetation map
      ! ALLOCATE( HUMTROP_GFED3( IGFED3, JGFED3 ), STAT=AS )
      ! IF ( AS /= 0 ) CALL ALLOC_ERR( 'HUMTROP_GFED3' )

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

      !=================================================================
      ! Read emission factors (which convert from kg DM to 
      ! either [molec species] or [atoms C]) from bpch file
      !=================================================================
     
      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE ) // !&'GFED3_201203/GFED3_emission_factors.txt'
     &           'GFED2_200601/GFED2_emission_factors_73t.txt'

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_flambe:1' )

      ! Skip header lines
      DO N = 1, 6
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_flambe:2' )
      ENDDO

      ! Read emission factors for each species
      DO N = 1, N_SPEC
         READ( IU_FILE, 100, IOSTAT=IOS ) 
     &     NDUM, FLAMBE_SPEC_NAME(N), ( FLAMBE_EMFAC(N,M), M=1,N_EMFAC )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_flambe:3' )
      WRITE (6,100) NDUM, FLAMBE_SPEC_NAME(N), 
     &                   (FLAMBE_EMFAC(N,M),M=1,N_EMFAC)
      ENDDO
      
      ! FORMAT string
 100  FORMAT( 1x, i2, 1x, a4, 3(3x,es14.6) )

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Read FLAMBE vegetation map from bpch file
      ! 
      ! Values:  3 = boreal forest 
      !          2 = tropical forest; 
      !          1 = savanna / herb / other land
      !          0 = water
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_NATIVE ) //
     &           'GFED2_200601/GFED2_vegmap.generic.1x1'

      ! Read FLAMBE veg map
      !  CALL READ_BPCH2_GFED3( FILENAME, 'LANDMAP',     1,
      ! &                       0d0,      I1x1,          J1x1,     
      ! &                       1,        ARRAY,         QUIET=.TRUE. )
      CALL READ_BPCH2( FILENAME, 'LANDMAP',  1,
     &                 0d0,       I1x1,      J1x1-1,
     &                 1,         ARRAY,     QUIET=.TRUE. )      

      ! Cast from REAL*4 to INTEGER
      VEG_GEN_1x1(:,:) = ARRAY(:,:,1)
      
      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the FLAMBE biomass species
      !=================================================================
      
      ! Save correspondance between FLAMBE species order (N) and 
      ! species order of the simulation (IDBxxxs).(ccc, 2/4/10)
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( TRIM( FLAMBE_SPEC_NAME(N) ) ) 
            CASE( 'NOx'  )
               ! Replaced NOx with NO (mpayer, 3/14/13)
               BIO_SAVE(N)          = IDBNO
               FLAMBE_SPEC_MOLWT(N) = 14d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO'   )
               BIO_SAVE(N)          = IDBCO
               FLAMBE_SPEC_MOLWT(N) = 28d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ALK4' )
               BIO_SAVE(N)          = IDBALK4
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ACET' )
               BIO_SAVE(N)          = IDBACET
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MEK'  )
               BIO_SAVE(N)          = IDBMEK
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ALD2' )
               BIO_SAVE(N)          = IDBALD2
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'PRPE' )
               BIO_SAVE(N)          = IDBPRPE
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C3H8' )
               BIO_SAVE(N)          = IDBC3H8
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'CH2O' )
               BIO_SAVE(N)          = IDBCH2O
               FLAMBE_SPEC_MOLWT(N) = 30d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'C2H6' )
               BIO_SAVE(N)          = IDBC2H6
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'SO2'  )
               BIO_SAVE(N)          = IDBSO2
               FLAMBE_SPEC_MOLWT(N) = 64d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'NH3'  )
               BIO_SAVE(N)          = IDBNH3
               FLAMBE_SPEC_MOLWT(N) = 17d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BC'   )
               BIO_SAVE(N)          = IDBBC
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'OC'   )
               BIO_SAVE(N)          = IDBOC
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYX' )
               BIO_SAVE(N)          = IDBGLYX
               FLAMBE_SPEC_MOLWT(N) = 58d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'MGLY' )
               BIO_SAVE(N)          = IDBMGLY
               FLAMBE_SPEC_MOLWT(N) = 72d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'BENZ' )
               BIO_SAVE(N)          = IDBBENZ
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'TOLU' )
               BIO_SAVE(N)          = IDBTOLU
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'XYLE' )
               BIO_SAVE(N)          = IDBXYLE
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H4' )
               BIO_SAVE(N)          = IDBC2H4
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C2H2' )
               BIO_SAVE(N)          = IDBC2H2
               FLAMBE_SPEC_MOLWT(N) = 12d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'GLYC' )
               BIO_SAVE(N)          = IDBGLYC
               FLAMBE_SPEC_MOLWT(N) = 60d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'HAC' )
               BIO_SAVE(N)          = IDBHAC
               FLAMBE_SPEC_MOLWT(N) = 74d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CO2'  )
               BIO_SAVE(N)          = IDBCO2
               FLAMBE_SPEC_MOLWT(N) = 44d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'CH4' )
               BIO_SAVE(N)          = IDBCH4
               FLAMBE_SPEC_MOLWT(N) = 16d-3
               FLAMBE_SPEC_UNIT(N)  = '[Tg  ]'
            CASE DEFAULT
               ! Nothing
               BIO_SAVE(N)         = 0

            WRITE(*,*) 'NAME',TRIM( FLAMBE_SPEC_NAME(N) )
         END SELECT
      ENDDO

      END SUBROUTINE INIT_FLAMBE_BIOMASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rearrange_biom
!
! !DESCRIPTION: Subroutine REARRANGE\_BIOM takes FLAMBE emissions (which have 
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
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: N

      ! Loop over FLAMBE species
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
! !IROUTINE: cleanup_flambe_biomass
!
! !DESCRIPTION: Subroutine CLEANUP\_FLAMBE\_BIOMASS deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_FLAMBE_BIOMASS
! 
! !REVISION HISTORY: 
!  21 Jun 2013 - S. Kim - Initial version, based on GFED3
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_FLAMBE_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED(FLAMBE_EMFAC     ) ) DEALLOCATE( FLAMBE_EMFAC     )
      IF ( ALLOCATED(FLAMBE_SPEC_MOLWT) ) DEALLOCATE( FLAMBE_SPEC_MOLWT)
      IF ( ALLOCATED(FLAMBE_SPEC_NAME ) ) DEALLOCATE( FLAMBE_SPEC_NAME )
      IF ( ALLOCATED(VEG_GEN_1x1      ) ) DEALLOCATE( VEG_GEN_1x1      )
      IF ( ALLOCATED(BIOMASS_MODEL    ) ) DEALLOCATE( BIOMASS_MODEL    )
     
      END SUBROUTINE CLEANUP_FLAMBE_BIOMASS
!EOC

      END MODULE FLAMBE_BIOMASS_MOD
