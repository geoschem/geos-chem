! $Id: scale_anthro_mod.f,v 1.3 2009/01/29 15:57:20 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: SCALE_ANTHRO_MOD
!
! !DESCRIPTION: Module SCALE\_ANTHRO\_MOD contains routines to scale 
!  anthropogenic emissions from a base year to a simulation year 
!  (avm, phs, 1/29/08)
!\\
!\\
! !INTERFACE: 
!
      MODULE SCALE_ANTHRO_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GET_ANNUAL_SCALAR
      PUBLIC  :: GET_ANNUAL_SCALAR_1x1 
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!EOP
!------------------------------------------------------------------------------

      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_ANNUAL_SCALAR
!
! !DESCRIPTION: Subroutine GET\_ANNUAL\_SCALAR returns annual scale 
!  factors to convert D\_YEAR (base year) to N\_YEAR (simulation year), 
!  on the current model resolution (avd, bmy, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR( TRACER, D_YEAR, N_YEAR, AS )
!
! !USES:
!
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1

#     include "CMN_SIZE"                         ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: TRACER           ! Tracer number
      INTEGER, INTENT(IN)    :: D_YEAR           ! Base year of emissions
      INTEGER, INTENT(IN)    :: N_YEAR           ! Target year of emissions
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*4,  INTENT(INOUT) :: AS(IIPAR,JJPAR)  ! Scale factor array
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                        :: AS_1x1(I1x1,J1x1)
      REAL*8                        :: AS_1x1x1(I1x1,J1x1,1)
      REAL*8                        :: AS_R8(IIPAR,JJPAR)

      ! Read 1x1 scale factors
      CALL GET_ANNUAL_SCALAR_1x1( TRACER, D_YEAR, N_YEAR, AS_1x1 )

      ! Cast to REAL*8
      AS_1x1x1(:,:,1) = AS_1x1(:,:)

      ! Regrid emissions factors to current model resolution
      CALL DO_REGRID_1x1( 'unitless', AS_1x1x1, AS_R8 )

      ! Cast to REAL*4
      AS(:,:) = AS_R8(:,:)

      ! Return to calling program
      END SUBROUTINE GET_ANNUAL_SCALAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_ANNUAL_SCALAR_1x1
!
! !DESCRIPTION: Subroutine GET\_ANNUAL\_SCALAR\_1x1 returns annual scale 
!  factors to convert D\_YEAR (base year) to N\_YEAR (target year), on the 1x1 
!  GEOS-Chem grid. (avd, bmy, phs, 1/28/09) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ANNUAL_SCALAR_1x1( TRACER, D_YEAR, N_YEAR, AS_1x1 )
!
! !USES:
!
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE BPCH2_MOD,     ONLY : GET_TAU0, READ_BPCH2

#     include "CMN_SIZE"                           ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: TRACER             ! Tracer number
      INTEGER, INTENT(IN)    :: D_YEAR             ! Base year of emissions
      INTEGER, INTENT(IN)    :: N_YEAR             ! Target year of emissions
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8,   INTENT(OUT)  :: AS_1x1(I1x1,J1x1)  ! Scale factor array
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!  (1) Scaling factors are for years between 1985 and 2005, on the GEOS-Chem
!       1x1 grid (phs, 3/10/08)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
!EOC


