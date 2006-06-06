! $Id: gfed2_biomass_mod.f,v 1.2 2006/06/06 14:26:04 bmy Exp $
      MODULE GFED2_BIOMASS_MOD
!
!******************************************************************************
!  Module GFED2_BIOMASS_MOD contains variables and routines to compute the
!  GFED2 biomass burning emissions. (psk, bmy, 4/20/06, 5/30/06)
!
!  Monthly emissions of C are read from disk and then multiplied by the 
!  appropriate emission factors to produce biomass burning emissions on a 
!  "generic" 1x1 grid.  The emissions are then regridded to the current 
!  GEOS-Chem or GCAP grid (1x1, 2x25, or 4x5).
!
!  GFED2 biomass burning emissions are computed for the following 
!  gas-phase species:
!
!     (1 ) NOx  [  molec/cm2/s] 
!     (2 ) CO   [  molec/cm2/s]                
!     (3 ) ALK4 [atoms C/cm2/s]   
!     (4 ) ACET [atoms C/cm2/s]   
!     (5 ) MEK  [atoms C/cm2/s]   
!     (6 ) ALD2 [atoms C/cm2/s]   
!     (7 ) PRPE [atoms C/cm2/s]   
!     (8 ) C3H8 [atoms C/cm2/s]   
!     (9 ) CH2O [  molec/cm2/s]   
!     (10) C2H6 [atoms C/cm2/s]   
!
!  Module Variables:
!  ============================================================================
!  (1 ) IDBNOx          (INTEGER) : Local index for NOx  in BIOM_OUT array
!  (2 ) IDBCO           (INTEGER) : Local index for CO   in BIOM_OUT array
!  (3 ) IDBALK4         (INTEGER) : Local index for ALK4 in BIOM_OUT array
!  (4 ) IDBACET         (INTEGER) : Local index for ACET in BIOM_OUT array
!  (5 ) IDBMEK          (INTEGER) : Local index for MEK  in BIOM_OUT array
!  (6 ) IDBALD2         (INTEGER) : Local index for ALD2 in BIOM_OUT array
!  (7 ) IDBPRPE         (INTEGER) : Local index for PRPE in BIOM_OUT array
!  (8 ) IDBC3H8         (INTEGER) : Local index for C3H8 in BIOM_OUT array
!  (9 ) IDBCH2O         (INTEGER) : Local index for CH2O in BIOM_OUT array
!  (10) IDBC2H6         (INTEGER) : Local index for C2H6 in BIOM_OUT array
!  (11) SECONDS         (REAL*8 ) : Number of seconds in the current month
!  (12) N_EMFAC         (INTEGER) : Number of emission factors per species
!  (13) N_SPEC          (INTEGER) : Number of species
!  (14) VEG_GEN_1x1     (REAL*8 ) : Array for GFED2 1x1 vegetation ma
!  (15) GFED2_SPEC_NAME (CHAR*4 ) : Array for GFED2 biomass species names
!  (16) GFED2_EMFAC     (REAL*8 ) : Array for user-defined emission factors
!  (17) BIOM_OUT        (REAL*8 ) : Array for biomass emissions on model grid
!
!  Module Routines:
!  ============================================================================
!  (1 ) GFED2_COMPUTE_BIOMASS     : Computes biomass emissions once per month
!  (2 ) GFED2_SCALE_FUTURE        : Applies IPCC future scale factors to GFED2
!  (3 ) GFED2_TOTAL_Tg            : Totals GFED2 biomass emissions [Tg/month]
!  (4 ) INIT_GFED2_BIOMASS        : Initializes arrays and reads startup data
!  (5 ) CLEANUP_GFED2_BIOMASS     : Deallocates all module arrays
!
!  GEOS-Chem modules referenced by "gfed2_biomass_mod.f":
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f            : Module w/ error and NaN check routines
!  (4 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (5 ) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (6 ) grid_mod.f             : Module w/ horizontal grid information
!  (7 ) time_mod.f             : Module w/ routines for computing time & date
!  (8 ) regrid_1x1_mod.f       : Module w/ routines for regrid 1x1 data
!
!  References:
!  ============================================================================
!  (1 ) Original GFED2 database from Jim Randerson:
!        http://ess1.ess.uci.edu/~jranders/data/GFED2/
!  (2 ) Giglio, L., G.R. van der Werf, J.T. Randerson, G.J. Collatz, and
!        P. Kasibhatla, "Global estimation of burned area using MODIS active
!        fire observations", Atm. Chem. Phys. Discuss, Vol 5, 11091, 2005.
!        http://www.copernicus.org/EGU/acp/acpd/5/11091/acpd-5-11091.pdf
!  (3 ) G.R. van der Werf, J.T. Randerson, L. Giglio, G.J. Collatz, 
!        P.S. Kasibhatla, and A.F. Arellano, Jr., "Interannual variability
!        in global biomass burning emissions from 1997 to 2004", Atm. Chem.
!        Phys. Discuss., submitted, 2005, 
!        http://sheba.geo.vu.nl/~gwerf/pubs/VanderWerfEA2005ACPD.pdf
!
!  NOTES:
!  (1 ) Added private routine GFED2_SCALE_FUTURE (swu, bmy, 5/30/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gfed2_biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: GFED2_COMPUTE_BIOMASS
      PUBLIC :: CLEANUP_GFED2_BIOMASS

      !==================================================================
      ! MODULE VARIABLES
      !==================================================================

      ! Scalars
      INTEGER                       :: IDBNOx,  IDBCO,   IDBALK4
      INTEGER                       :: IDBACET, IDBMEK,  IDBALD2
      INTEGER                       :: IDBPRPE, IDBC3H8, IDBCH2O
      INTEGER                       :: IDBC2H6
      REAL*8                        :: SECONDS

      ! Parameters
      INTEGER,          PARAMETER   :: N_EMFAC = 3
      INTEGER,          PARAMETER   :: N_SPEC  = 10 

      ! Arrays
      INTEGER,          ALLOCATABLE :: VEG_GEN_1x1(:,:)
      REAL*8,           ALLOCATABLE :: GFED2_EMFAC(:,:)
      REAL*8,           ALLOCATABLE :: GFED2_SPEC_MOLWT(:)
      CHARACTER(LEN=4), ALLOCATABLE :: GFED2_SPEC_NAME(:)
      CHARACTER(LEN=6), ALLOCATABLE :: GFED2_SPEC_UNIT(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_COMPUTE_BIOMASS( THIS_YYYY, THIS_MM, BIOM_OUT )
!
!******************************************************************************
!  Subroutine GFED2_COMPUTE_BIOMASS computes the monthly GFED2 biomass burning
!  emissions for a given year and month. (psk, bmy, 4/20/06, 5/30/06)
!
!  This routine only has to be called on the first day of each month.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THIS_YYYY (INTEGER) : Current year 
!  (2 ) THIS_MM   (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Now references LFUTURE from "logical_mod.f".  Now call private routine
!        GFED2_SCALE_FUTURE to compute future biomass emissions, if necessary. 
!        (swu, bmy, 5/30/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2,    GET_TAU0
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE LOGICAL_MOD,    ONLY : LFUTURE
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      
#     include "CMN_SIZE"       ! Size parameters

      ! Arguments 
      INTEGER, INTENT(IN)     :: THIS_YYYY
      INTEGER, INTENT(IN)     :: THIS_MM
      REAL*8,  INTENT(OUT)    :: BIOM_OUT(IIPAR,JJPAR,N_SPEC)

      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: I,    J,  N,   N_VEG 
      INTEGER                 :: YYYY, MM, MM1, YYYY1, YYYYMMDD
      REAL*4                  :: DM_GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: BIOM_GEN_1x1(I1x1,J1x1-1,N_SPEC)
      REAL*8                  :: BIOM_GEOS_1x1(I1x1,J1x1,N_SPEC)
      REAL*8                  :: TAU0, TAU1 
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! GFED2_COMPUTE_BIOMASS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_GFED2_BIOMASS
         FIRST = .FALSE.
      ENDIF

      ! Save in local variables
      YYYY = THIS_YYYY
      MM   = THIS_MM

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 
     &  'G F E D 2   B I O M A S S   B U R N I N G   E M I S S I O N S'

      WRITE( 6, 100 ) YYYY, MM
 100  FORMAT( 'for year and month: ', i4, '/', i2.2, / )

      ! 1997 is the 1st year of available data
      IF ( YYYY < 1997 ) THEN
         YYYY = 1997
         WRITE( 6, 110 )
 110     FORMAT( 'YEAR < 1997; Using GFED2 biomass for 1997!' )
      ENDIF

      ! 2004 is currently the last year of available data
      IF ( YYYY > 2004 ) THEN
         YYYY = 2004
         WRITE( 6, 120 ) 
 120     FORMAT( 'YEAR > 2004; Using GFED2 biomass for 2004!' )
      ENDIF

      !=================================================================
      ! Read monthly GFED2 C emissions [g/m2/month]
      !=================================================================

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

      ! File name with GFED2 C emissions
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'GFED2_200601/YYYY/GFED2_C_YYYYMM.generic.1x1'

      ! Create YYYYMMDD integer value
      YYYYMMDD = YYYY*10000 + MM*100 + 01

      ! Replace YYYY/MM in the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read GFED2 C emissions [g C/m2/month]
      CALL READ_BPCH2( FILENAME, 'GFED2-BB',   99, 
     &                 TAU0,      I1x1,        J1x1-1,     
     &                 1,         DM_GEN_1x1,  QUIET=.TRUE. ) 

      !=================================================================
      ! Convert C [g/m2/month] to dry matter burned [kg/cm2/month]
      !
      ! Unit Conversions:
      ! (1) C    to DM    --> Divide by 0.45  
      ! (2) g    to kg    --> Divide by 1000  
      ! (3) 1/m2 to 1/cm2 --> Divide by 10000 
      !=================================================================

      ! Loop over GENERIC 1x1 GRID
      DO J = 1, J1x1-1
      DO I = 1, I1x1

         ! Set negatives to zero
         DM_GEN_1x1(I,J) = MAX( DM_GEN_1x1(I,J), 0e0 )

         ! Convert [g C/m2/month] to [kg DM/cm2/month]
         DM_GEN_1x1(I,J) = DM_GEN_1x1(I,J) / ( 0.45d0 * 1d3 * 1d4 )

      ENDDO
      ENDDO

      !=================================================================
      ! Calculate biomass species emissions on 1x1 emissions grid
      !
      ! Emission factors convert from [kg/cm2/month] to either
      ! [molec/cm2/month] or [atoms C/cm2/month]
      !
      ! Units:
      !  [  molec/cm2/month] : NOx,  CO,   CH2O 
      !  [atoms C/cm2/month] : ALK4, ACET, MEK, ALD2, PRPE  C3H8, C2H6
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
            
            ! Multiply DM * EMISSION FACTOR to get biomass emissions
            ! for each species on the GENERIC 1x1 GRID 
            SELECT CASE( N_VEG )

               ! Ocean 
               CASE( 0 ) 
                  BIOM_GEN_1x1(I,J,N) = 0d0

               ! Land
               CASE( 1:3 )
                  BIOM_GEN_1x1(I,J,N) = DM_GEN_1x1(I,J) * 
     &                                  GFED2_EMFAC(N,N_VEG)

               ! Otherwise
               CASE DEFAULT
                  ! Nothing

            END SELECT
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Regrid each species from GENERIC 1x1 GRID to GEOS-Chem 1x1 GRID
         CALL DO_REGRID_G2G_1x1( BIOM_GEN_1x1(:,:,N), 
     &                           BIOM_GEOS_1x1(:,:,N),
     &                           PER_UNIT_AREA=.TRUE. )
      ENDDO

      ! Regrid from GEOS 1x1 grid to current grid.  (The unit 'molec/cm2' 
      ! is just used to denote that the quantity is per unit area.)
      CALL DO_REGRID_1x1( N_SPEC,       'molec/cm2', 
     &                    BIOM_GEOS_1x1, BIOM_OUT ) 

      !### Debug
      !###! Print totals on 1x1 generic, 1x1 geos, & output grids
      !###CALL GFED2_DEBUG_PRINT( BIOM_GEN_1x1, BIOM_GEOS_1x1, BIOM_OUT )

      ! Compute future biomass emissions (if necessary)
      IF ( LFUTURE ) THEN
         CALL GFED2_SCALE_FUTURE( BIOM_OUT )
      ENDIF

      ! Print totals in Tg/month
      CALL GFED2_TOTAL_Tg( THIS_YYYY, THIS_MM, BIOM_OUT )

      ! Convert from [molec/cm2/month] to [molec/cm2/s]
      BIOM_OUT = BIOM_OUT / SECONDS

      ! Echo info
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE GFED2_COMPUTE_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_SCALE_FUTURE( BB )
!
!******************************************************************************
!  Subroutine GFED2_SCALE_FUTURE applies the IPCC future emissions scale 
!  factors to the GFED2 biomass burning emisisons in order to compute the 
!  future emissions of biomass burning for NOx, CO, and VOC's.  
!  (swu, bmy, 5/30/06)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) BB (REAL*8   ) : Array w/ biomass burning emisisons [molec/cm2]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_CObb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_NOxbb
      USE FUTURE_EMISSIONS_MOD,   ONLY : GET_FUTURE_SCALE_VOCbb

#     include "CMN_SIZE"               ! Size parameters

      ! Arguments
      REAL*8,           INTENT(INOUT) :: BB(IIPAR,JJPAR,N_SPEC)

      ! Local variables
      INTEGER                         :: I, J, N
      
      !=================================================================
      ! GFED2_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )
      DO N = 1, N_SPEC
      DO J = 1, JJPAR
      DO I = 1, IIPAR 

         IF ( N == IDBNOx ) THEN

            ! Future biomass NOx [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_NOxbb( I, J )

         ELSE IF ( N == IDBCO ) THEN

            ! Future biomass CO [molec/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_CObb( I, J )

         ELSE

            ! Future biomass Hydrocarbons [atoms C/cm2]
            BB(I,J,N) = BB(I,J,N) * GET_FUTURE_SCALE_VOCbb( I, J )

         ENDIF
         
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE GFED2_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE GFED2_TOTAL_Tg( YYYY, MM, BIOMASS )
!
!******************************************************************************
!  Subroutine TOTAL_BIOMASS_TG prints the amount of biomass burning emissions 
!  that are emitted each month in Tg or Tg C. (bmy, 3/20/01, 4/20/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYY    (INTEGER) : Current year
!  (2 ) MM      (INTEGER) : Currrent month
!  (3 ) BIOMASS (REAL*8)  : Biomass burning emissions [molec/cm2/month]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,   ONLY : GET_AREA_CM2

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYY, MM
      REAL*8,  INTENT(IN) :: BIOMASS(IIPAR,JJPAR,N_SPEC) 

      ! Local variables
      INTEGER             :: I,    J,     N
      REAL*8              :: CONV, MOLWT, TOTAL
      CHARACTER(LEN=4)    :: NAME
      CHARACTER(LEN=6)    :: UNIT

      !=================================================================
      ! GFED2_TOTAL_Tg begins here!
      !=================================================================

      ! Loop over biomass species
      DO N = 1, N_SPEC

         ! Initialize
         NAME  = GFED2_SPEC_NAME(N)
         MOLWT = GFED2_SPEC_MOLWT(N)
         UNIT  = GFED2_SPEC_UNIT(N)
         TOTAL = 0d0

         ! Loop over latitudes
         DO J = 1, JJPAR
         
            ! Convert to [Tg/month] (or [Tg C/month] for HC's)
            CONV = GET_AREA_CM2( J ) * ( MOLWT / 6.023d23 ) * 1d-9

            ! Loop over longitudes
            DO I = 1, IIPAR
               TOTAL = TOTAL + ( BIOMASS(I,J,N) * CONV )
            ENDDO
         ENDDO
     
         ! Write totals
         WRITE( 6, 110 ) NAME, TOTAL, UNIT
 110     FORMAT( 'Sum Biomass ', a4, 1x, ': ', f8.3, 1x, a6 )
      ENDDO

      ! Return to calling program
      END SUBROUTINE GFED2_TOTAL_Tg

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GFED2_BIOMASS
!
!******************************************************************************
!  Subroutine INIT_GFED2_BIOMASS allocates all module arrays.  It also reads
!  the emission factors and vegetation map files at the start of a GEOS-Chem
!  simulation. (psk, bmy, 4/20/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD,      ONLY : IOERROR, IU_FILE

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: AS, IOS, M, N, NDUM
      REAL*4                 :: ARRAY(I1x1,J1x1-1,1)
      CHARACTER(LEN=255)     :: FILENAME
      
      !=================================================================
      ! INIT_GFED2_BIOMASS begins here!
      !=================================================================

      ! Allocate array for emission factors
      ALLOCATE( GFED2_EMFAC( N_SPEC, N_EMFAC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_EMFAC' )
      GFED2_EMFAC = 0d0
      
      ! Allocate array for species molecular weight
      ALLOCATE( GFED2_SPEC_MOLWT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_MOLWT' )
      GFED2_SPEC_MOLWT = 0d0

      ! Allocate array for species name
      ALLOCATE( GFED2_SPEC_NAME( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_NAME' )
      GFED2_SPEC_NAME = ''

      ! Allocate array for species molecular weight
      ALLOCATE( GFED2_SPEC_UNIT( N_SPEC ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GFED2_SPEC_UNIT' )
      GFED2_SPEC_UNIT = ''

      ! Allocate array for vegetation map
      ALLOCATE( VEG_GEN_1x1( I1x1, J1x1-1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VEG_GEN_1x1' )

      !=================================================================
      ! Read emission factors (which convert from kg DM to 
      ! either [molec species] or [atoms C]) from bpch file
      !=================================================================
     
      ! File name
      FILENAME = TRIM( DATA_DIR_1x1) // 
     &           'GFED2_200601/GFED2_emission_factors.txt'

      ! Open emission factor file (ASCII format)
      OPEN( IU_FILE, file=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:1' )

      ! Skip header lines
      DO N = 1, 6 
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:2' )
      ENDDO

      ! Read emission factors for each species
      DO N = 1, N_SPEC
         READ( IU_FILE, 100, IOSTAT=IOS ) 
     &       NDUM, GFED2_SPEC_NAME(N), ( GFED2_EMFAC(N,M), M=1,N_EMFAC )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'init_gfed2:3' )
      ENDDO
      
      ! FORMAT string
 100  FORMAT( 1x, i2, 1x, a4, 3(3x,es14.6) )

      ! Close file
      CLOSE( IU_FILE )
      
      !=================================================================
      ! Read GFED2 vegetation map from bpch file
      ! 
      ! Values:  3 = boreal forest 
      !          2 = tropical forest; 
      !          1 = savanna / herb / other land
      !          0 = water
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'GFED2_200601/GFED2_vegmap.generic.1x1'

      ! Read GFED2 veg map 
      CALL READ_BPCH2( FILENAME, 'LANDMAP',  1, 
     &                 0d0,       I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to INTEGER
      VEG_GEN_1x1(:,:) = ARRAY(:,:,1)

      !=================================================================
      ! Define local ID flags and arrays for the names, units, 
      ! and molecular weights of the GFED2 biomass species
      !=================================================================
      
      ! Initialize 
      IDBNOx  = 0  
      IDBCO   = 0
      IDBALK4 = 0
      IDBACET = 0 
      IDBMEK  = 0 
      IDBALD2 = 0
      IDBPRPE = 0
      IDBC3H8 = 0
      IDBCH2O = 0
      IDBC2H6 = 0
 
      ! Save species # in IDBxxxx flags for future reference
      ! and also initialize arrays for mol wts and units
      DO N = 1, N_SPEC
         SELECT CASE ( GFED2_SPEC_NAME(N) ) 
            CASE( 'NOx ' )
               IDBNOx              = N
               GFED2_SPEC_MOLWT(N) = 14d-3
               GFED2_SPEC_UNIT(N)  = '[Tg N]'
            CASE( 'CO  ' )
               IDBCO               = N
               GFED2_SPEC_MOLWT(N) = 28d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'ALK4' )
               IDBALK4             = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ACET' )
               IDBACET = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'MEK ' )
               IDBMEK  = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'ALD2' )
               IDBALD2 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'PRPE' )
               IDBPRPE = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'C3H8' )
               IDBC3H8 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE( 'CH2O' )
               IDBCH2O = N
               GFED2_SPEC_MOLWT(N) = 30d-3
               GFED2_SPEC_UNIT(N)  = '[Tg  ]'
            CASE( 'C2H6' )
               IDBC2H6 = N
               GFED2_SPEC_MOLWT(N) = 12d-3
               GFED2_SPEC_UNIT(N)  = '[Tg C]'
            CASE DEFAULT
               ! Nothing
         END SELECT
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_GFED2_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GFED2_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_GFED2_BIOMASS deallocates all module arrays.
!  (psk, bmy, 4/20/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_GFED2_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( GFED2_EMFAC      ) ) DEALLOCATE( GFED2_EMFAC     )
      IF ( ALLOCATED( GFED2_SPEC_MOLWT ) ) DEALLOCATE( GFED2_SPEC_MOLWT)
      IF ( ALLOCATED( GFED2_SPEC_NAME  ) ) DEALLOCATE( GFED2_SPEC_NAME )
      IF ( ALLOCATED( VEG_GEN_1x1      ) ) DEALLOCATE( VEG_GEN_1x1     )
      
      ! Return to calling program
      END SUBROUTINE CLEANUP_GFED2_BIOMASS

!------------------------------------------------------------------------------
! This routine is for debugging only: uncomment it & recompile (bmy, 4/20/06)
!
!      SUBROUTINE GFED2_DEBUG_PRINT( GEN_1x1, GEOS_1x1, OUT )
!!
!!*****************************************************************************
!!  Subroutine PRINT_TOTALS prints monthly and cumulative totals of the GFED2 
!!  biomass burning emissions on the GENERIC 1x1 grid, the GEOS 1x1 grid, and 
!!  the current output grid (e.g. 2x25, 4x5). (bmy, 4/20/06)
!!
!!  Arguments as Input:
!!  ===========================================================================
!!  (1 ) GEN_1x1  (REAL*8) : Array of biomass burning on GENERIC   1x1    grid
!!  (2 ) GEOS_1x1 (REAL*8) : Array of biomass burning on GEOS-Chem 1x1    grid
!!  (3 ) OUT      (REAL*8) : Array of biomass burning on GEOS-Chem output grid
!!  
!!  NOTES:
!!*****************************************************************************
!!
!      ! References to F90 modules
!      USE GRID_MOD, ONLY : GET_AREA_CM2 
!      USE TIME_MOD, ONLY : ITS_A_NEW_YEAR
!
!#     include "CMN_SIZE" ! Size parameters
!#     include "CMN_GCTM" ! Re 
!
!      ! Local variables
!      INTEGER           :: I,        J,       N,      AS
!      REAL*8            :: SS,       NN,      Re_cm,  RLAT
!      REAL*8            :: SUM_GEOS, SUM_GEN, SUM_OUT
!      REAL*8            :: A_GEN(J1x1-1)
!      REAL*8            :: A_GEOS(J1x1)
!      REAL*8            :: GEN_1x1(I1x1,J1x1-1,N_SPEC)
!      REAL*8            :: GEOS_1x1(I1x1,J1x1,N_SPEC)
!      REAL*8            :: OUT(IIPAR,JJPAR,N_SPEC)
!      REAL*8            :: XNUMOL(N_SPEC)
!      REAL*8            :: YEDGE(J1x1+1)
!      REAL*8, SAVE      :: CUM_GEOS(N_SPEC)
!      REAL*8, SAVE      :: CUM_GEN(N_SPEC)
!      REAL*8, SAVE      :: CUM_OUT(N_SPEC)
!
!      !=================================================================
!      ! GFED2_DEBUG_PRINT begins here!      
!      !=================================================================
!
!      ! Zero cumulative totals once per year
!      IF ( ITS_A_NEW_YEAR() ) THEN
!         CUM_GEN  = 0d0
!         CUM_GEOS = 0d0
!         CUM_OUT  = 0d0
!      ENDIF
!
!      ! Radius of the earth [cm]
!      Re_cm       = Re * 100d0
!
!      ! XNUMOL: factor to convert from molecules to Tg
!      XNUMOL(1)   = 14d-3 / 6.0225d23 * 1d-9  ! NOx
!      XNUMOL(2)   = 28d-3 / 6.0225d23 * 1d-9  ! CO
!      XNUMOL(3)   = 12d-3 / 6.0225d23 * 1d-9  ! ALK4
!      XNUMOL(4)   = 12d-3 / 6.0225d23 * 1d-9  ! ACET
!      XNUMOL(5)   = 12d-3 / 6.0225d23 * 1d-9  ! MEK
!      XNUMOL(6)   = 12d-3 / 6.0225d23 * 1d-9  ! ALD2
!      XNUMOL(7)   = 12d-3 / 6.0225d23 * 1d-9  ! PRPE
!      XNUMOL(8)   = 12d-3 / 6.0225d23 * 1d-9  ! C3H8
!      XNUMOL(9)   = 30d-3 / 6.0225d23 * 1d-9  ! C2HO
!      XNUMOL(10)  = 12d-3 / 6.0225d23 * 1d-9  ! C2H6
! 
!      !---------------------------------------
!      ! Surface area on GEOS-Chem 1x1 grid
!      ! Uses same algorithm from "grid_mod.f"
!      !---------------------------------------
!
!      ! Initialize
!      YEDGE(:)      = 0d0
!
!      ! 1x1 latitude edges
!      DO J = 2, J1x1
!         YEDGE(J)   = -90.5d0 + ( J - 1 )
!      ENDDO
!
!      ! Special cases at poles
!      YEDGE(1)      = -90.0d0
!      YEDGE(2)      = -89.5d0
!      YEDGE(J1x1+1) =  90.0d0
!
!      ! Compute 1x1 surface area
!      DO J = 1, J1x1
!
!         ! Lat at S and N edges of 1x1 box [radians]
!         SS         = PI_180 * YEDGE(J  )
!         NN         = PI_180 * YEDGE(J+1) 
!
!         ! S to N extent of grid box [unitless]
!         RLAT       = SIN( NN ) - SIN( SS )
!
!         ! 1x1 surface area [m2] (see "grid_mod.f" for algorithm)
!         A_GEOS(J)  = 2d0 * PI * Re_cm**2 / DBLE( I1x1 ) * RLAT
!      ENDDO
!
!      !---------------------------------------
!      ! Surface area on GENERIC 1x1 grid
!      ! Uses same algorithm from "grid_mod.f"
!      !---------------------------------------
!
!      ! Initialize
!      YEDGE(:)    = 0d0
!
!      ! 1x1 latitude edges
!      DO J = 1, J1x1
!         YEDGE(J) = -90d0 + ( J - 1 )
!      ENDDO
!
!      ! Compute 1x1 surface area
!      DO J = 1, J1x1-1
!
!         ! Lat at S and N edges of 1x1 box [radians]
!         SS       = PI_180 * YEDGE(J  )
!         NN       = PI_180 * YEDGE(J+1) 
!
!         ! S to N extent of grid box [unitless]
!         RLAT     = SIN( NN ) - SIN( SS )
!
!         ! GENERIC GRID 1x1 surface area [cm2]
!         A_GEN(J) = 2d0 * PI * Re_cm**2 / DBLE( I1x1 ) * RLAT
!      ENDDO  
!
!      !---------------------------------------
!      ! Compute sums in molec (or atoms C)
!      ! as well as in Tg (or Tg C)
!      !---------------------------------------
!
!      ! echo inf
!      PRINT*, 'MONTHLY TOTALS'
!      PRINT*, 'SPECIES, GEN1x1, GEOS1x1, GEOS_OUT [molec and Tg]'
!
!      ! Loop over # of species
!      DO N = 1, N_SPEC
!
!         ! Sum on generic grid
!         SUM_GEN = 0d0
!         DO J = 1, J1x1-1
!         DO I = 1, I1x1
!            SUM_GEN = SUM_GEN + ( GEN_1x1(I,J,N) * A_GEN(J)  )
!         ENDDO
!         ENDDO
!
!         ! Sum on GEOS grid
!         SUM_GEOS = 0d0
!         DO J = 1, J1x1
!         DO I = 1, I1x1
!            SUM_GEOS = SUM_GEOS + ( GEOS_1x1(I,J,N) * A_GEOS(J) )
!         ENDDO
!         ENDDO
!
!         ! Sum on current grid (e.g. GEOS-4x5)
!         SUM_OUT = 0d0
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            SUM_OUT = SUM_OUT + ( OUT(I,J,N) * GET_AREA_CM2(J)  )
!         ENDDO
!         ENDDO
!
!         ! Update cumulative totals 
!         CUM_GEN(N)  = CUM_GEN(N)  + SUM_GEN
!         CUM_GEOS(N) = CUM_GEOS(N) + SUM_GEOS
!         CUM_OUT(N)  = CUM_OUT(N)  + SUM_OUT
!
!         ! Print sums in molecules and in Tg
!         WRITE( 6, '(a4, 2x, 3(es13.6,1x,f10.3,3x))' ) 
!     &     GFED2_SPEC_NAME(N), SUM_GEN,  SUM_GEN  * XNUMOL(N),
!     &                         SUM_GEOS, SUM_GEOS * XNUMOL(N),
!     &                         SUM_OUT,  SUM_OUT  * XNUMOL(N)
!      ENDDO
!
!      ! Echo info
!      PRINT*    
!      PRINT*, 'CUMULATIVE TOTALS!'
!      PRINT*, 'SPECIES, GEN1x1, GEOS1x1, GEOS_OUT [molec and Tg]'
!
!      ! Also print cumulative totals [molec and Tg]
!      DO N = 1, N_SPEC
!         WRITE( 6, '(a4, 2x, 3(es13.6,1x,f10.3,3x))' ) 
!     &     GFED2_SPEC_NAME(N), CUM_GEN(N),  CUM_GEN(N)  * XNUMOL(N),
!     &                         CUM_GEOS(N), CUM_GEOS(N) * XNUMOL(N),
!     &                         CUM_OUT(N),  CUM_OUT(N)  * XNUMOL(N)
!      ENDDO
!
!      ! Echo info
!      PRINT*
!
!      ! Return to calling program
!      END SUBROUTINE GFED2_DEBUG_PRINT
!
!------------------------------------------------------------------------------

      ! End of module 
      END MODULE GFED2_BIOMASS_MOD
