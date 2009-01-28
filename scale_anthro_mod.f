      MODULE SCALE_ANTHRO_MOD
!
!******************************************************************************
!  Module SCALE_ANTHRO_MOD contains routines to scale anthropogenic 
!  emissions from a base year to a simulation year (avd, phs, 3/10/08)
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_ANNUAL_SCALAR     : Returns int'annual scale factor file
!  (2 ) GET_ANNUAL_SCALAR_1x1 : Returns int'annual scale factor file at 1x1
!  (3 ) INIT_SCALE_ANTHRO     : Allocates and zeroes module arrays
!  (4 ) CLEANUP_SCALE_ANTHRO  : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "scale_anthro_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f       : Module w/ GEOS-CHEM data & met field dirs
!  (3 ) regrid_1x1_mod.f      : Module w/ routines to regrid 1x1 data
!
!  REFERENCE:
!  (1 ) van Donkelaar et al., ACPD, 8, 4017-4057, 2008
!
!  NOTES:
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "scale_anthro_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: GET_ANNUAL_SCALAR
      PUBLIC :: GET_ANNUAL_SCALAR_1x1

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      ! None yet

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_ANNUAL_SCALAR( TRACER, D_YEAR, N_YEAR, AS )
!
!******************************************************************************
!  Subroutine GET_ANNUAL_SCALAR returns annual scale factors to convert
!  D_YEAR (base year) to N_YEAR (simulation year), on the current model 
!  resolution (avd, bmy, phs, 7/14/06, 3/10/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TRACER   (INTEGER  ) : Tracer number
!  (2 ) DEN_YEAR (INTEGER  ) : Inventory base year
!  (3 ) NUM_YEAR (INTEGER  ) : Year of simulation
!  (4 ) AS       (REAL*4   ) : Array to hold scale factors
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE REGRID_1x1_MOD,       ONLY : DO_REGRID_1x1

#     include "CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER,           INTENT(IN) :: TRACER, D_YEAR, N_YEAR
      REAL*4,         INTENT(INOUT) :: AS(IIPAR,JJPAR)

      ! Local variables
      REAL*8                        :: AS_1x1(I1x1,J1x1)
      REAL*8                        :: AS_1x1x1(I1x1,J1x1,1)
      REAL*8                        :: AS_R8(IIPAR,JJPAR)

      CALL GET_ANNUAL_SCALAR_1x1( TRACER, D_YEAR, N_YEAR, AS_1x1 )

      AS_1x1x1(:,:,1) = AS_1x1(:,:)

      ! Regrid emissionS factors to current model resolution
      CALL DO_REGRID_1x1( 'unitless', AS_1x1x1, AS_R8 )

      AS(:,:) = AS_R8(:,:)

      ! Return to calling program
      END SUBROUTINE GET_ANNUAL_SCALAR

!------------------------------------------------------------------------------

      SUBROUTINE GET_ANNUAL_SCALAR_1x1( TRACER, D_YEAR, N_YEAR, AS_1x1 )
!
!******************************************************************************
!  Subroutine GET_ANNUAL_SCALAR_1x1 returns annual scale factors to convert
!  DEN_YEAR to NUM_YEAR, on the 1X1 GEOS-Chem grid.
!  (avd, bmy, phs, 3/10/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TRACER   (INTEGER  ) : Tracer number
!  (2 ) DEN_YEAR (INTEGER  ) : Inventory base year
!  (3 ) NUM_YEAR (INTEGER  ) : Year of simulation
!  (4 ) AS_1x1   (REAL*8   ) : Array to hold scale factors
!
!  NOTES:
!  (1 ) Scaling factors are for years between 1985 and 2005, on the GEOS-Chem
!        1x1 grid (phs, 3/10/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,        ONLY : DATA_DIR_1x1
      USE BPCH2_MOD,            ONLY : GET_TAU0, READ_BPCH2

#     include "CMN_SIZE"             ! Size parameters

      ! Arguments
      INTEGER,          INTENT(IN)  :: TRACER, D_YEAR, N_YEAR
      REAL*8,           INTENT(OUT) :: AS_1x1(I1x1,J1x1)

      ! Local variables
      REAL*4                        :: D_1x1(I1x1,J1x1)
      REAL*4                        :: N_1x1(I1x1,J1x1)
      REAL*8                        :: TAU2000
      CHARACTER(LEN=255)            :: FILENAME,     SCALE_DIR
      CHARACTER(LEN=4)              :: DEN_YYYY_STR, NUM_YYYY_STR
      INTEGER                       :: DEN_YEAR,     NUM_YEAR
      INTEGER                       :: I, J

      !=================================================================
      ! GET_ANNUAL_SCALAR_1x1 begins here!
      !=================================================================

      SCALE_DIR = TRIM( DATA_DIR_1x1 ) // 'anth_scale_factors_200811/'

      ! limit scaling between available years
      DEN_YEAR = MAX( MIN( D_YEAR, 2005 ), 1985 )
      NUM_YEAR = MAX( MIN( N_YEAR, 2005 ), 1985 )

      WRITE( DEN_YYYY_STR, '(i4.4)' ) DEN_YEAR
      WRITE( NUM_YYYY_STR, '(i4.4)' ) NUM_YEAR

      IF ( DEN_YEAR == 2000 ) THEN

         N_1x1(:,:) = 1.d0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) // 'NOxScalar-' // 
     &                 DEN_YYYY_STR // '-' // '2000.geos.1x1'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) // 'COScalar-' // 
     &                 DEN_YYYY_STR // '-' // '2000.geos.1x1'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) // 'SOxScalar-' //
     &                 DEN_YYYY_STR // '-' // '2000.geos.1x1'

         ENDIF

         ! Get Tau    
         TAU2000 = GET_TAU0(1,1,2000)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - GET_ANNUAL_SCALAR_1x1: Reading ', a )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU2000,  I1x1,       J1x1,
     &                    1,        N_1x1,      QUIET=.TRUE. )

      ENDIF

      IF ( NUM_YEAR == 2000 ) THEN

         D_1x1(:,:) = 1.d0

      ELSE

         ! Filename
         IF ( TRACER == 71 ) THEN

            ! NOx
            FILENAME = TRIM( SCALE_DIR ) // 'NOxScalar-' //
     &                 NUM_YYYY_STR // '-' // '2000.geos.1x1'

         ELSE IF ( TRACER == 72 ) THEN

            ! CO
            FILENAME = TRIM( SCALE_DIR ) // 'COScalar-' // 
     &                 NUM_YYYY_STR // '-' // '2000.geos.1x1'

         ELSE IF ( TRACER == 73 ) THEN

            ! SOx
            FILENAME = TRIM( SCALE_DIR ) // 'SOxScalar-' // 
     &                 NUM_YYYY_STR // '-' // '2000.geos.1x1'

         ENDIF

         ! Calc Tau
         TAU2000 = GET_TAU0(1,1,2000)

         ! Echo filename
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', TRACER,
     &                    TAU2000,  I1x1,       J1x1,
     &                    1,        D_1x1,      QUIET=.TRUE. )

      ENDIF

      ! Get scaling and cast as real*8
      AS_1x1(:,:) = D_1x1(:,:) / N_1x1(:,:)

      ! Return to calling program
      END SUBROUTINE GET_ANNUAL_SCALAR_1x1

!------------------------------------------------------------------------------

      ! End of module
      END MODULE SCALE_ANTHRO_MOD


