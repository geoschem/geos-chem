! $Id: vistas_anthro_mod.f,v 1.1 2009/01/28 19:59:00 bmy Exp $
      MODULE VISTAS_ANTHRO_MOD
!
!******************************************************************************
!  Module VISTAS_ANTHRO_MOD contains variables and routines to read the 
!  VISTAS anthropogenic emissions 
!  (amv, 11/24/2008)
!
!  Module Variables:
!  ============================================================================
!  (1 ) VISTAS_WD_NOx (REAL*8 ) : Anthopogenic Weekday NOx emission [molec/cm2/s]
!  (1 ) VISTAS_WE_NOx (REAL*8 ) : Anthopogenic Weekend NOx emission [molec/cm2/s]
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_VISTAS_ANTHRO     : Gets emissions at (I,J) for emissions species 
!  (2 ) EMISS_VISTAS_ANTHRO   : Reads VISTAS emissions from disk
!  (3 ) VISTAS_SCALE_FUTURE   : Applies IPCC future scale factors to emissions
!  (4 ) INIT_VISTAS_ANTHRO    : Allocates and zeroes module arrays
!  (5 ) CLEANUP_VISTAS_ANTHRO : Dealocates module arrays
!  (6 ) TOTAL_ANTHRO_TG       : Computes the annual total emissions
!
!  GEOS-Chem modules referenced by "vistas_anthro_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f        : Module w/ GEOS-Chem data & met field dirs
!  (3 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (4 ) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (5 ) grid_mod.f             : Module w/ horizontal grid information
!  (6 ) logical_mod.f          : Module w/ GEOS-Chem logical switches
!  (7 ) regrid_1x1_mod.f       : Module w/ routines to regrid 1x1 data  
!  (8 ) time_mod.f             : Module w/ routines for computing time & date
!  (9 ) tracerid_mod.f         : Module w/ pointers to tracers & emissions  
!
!  References:
!  ============================================================================
!  (1 ) Visibility Improvment State and Tribal Association of the Southeast.
!     See: http://www.vistas-sesarm.org/index.asp
!
!
!  NOTES: 
!******************************************************************************
!

      USE EPA_NEI_MOD,          ONLY : GET_USA_MASK

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "cac_anthro_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_VISTAS_ANTHRO
      PUBLIC :: EMISS_VISTAS_ANTHRO
      PUBLIC :: GET_VISTAS_ANTHRO

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Arrays
      REAL*8,  ALLOCATABLE :: VISTAS_WD_NOx(:,:)
      REAL*8,  ALLOCATABLE :: VISTAS_WE_NOx(:,:)
      REAL*8,  ALLOCATABLE :: A_CM2(:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_VISTAS_ANTHRO( I,    J,     N,  WEEKDAY,
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_VISTAS_ANTHRO returns the VISTAS emission for 
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].  (amv, 1/09/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude index
!  (3 ) N           (INTEGER) : GEOS-Chem tracer number
!  (4 ) WEEKDAY     (LOGICAL) : is it a weekday?
!  (5 ) KG_S        (LOGICAL) : OPTIONAL -- return emissions in [kg/s]
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNOx

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J, N
      LOGICAL, INTENT(IN)           :: WEEKDAY
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
      
      ! Local variables
      LOGICAL                       :: DO_KGS
      REAL*8                        :: VALUE

      !=================================================================
      ! GET_VISTA_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S

      IF ( N == IDTNOx ) THEN

         ! NOx [molec/cm2/s]
         IF ( WEEKDAY ) THEN
            VALUE = VISTAS_WD_NOx(I,J)
         ELSE
            VALUE = VISTAS_WE_NOx(I,J)
         ENDIF

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no VISTAS emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [molec/c,2/s] to [kg/s]
         VALUE = VALUE * A_CM2(J) / XNUMOL(N)

      ENDIF

      ! Return to calling program
      END FUNCTION GET_VISTAS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_VISTAS_ANTHRO
!
!******************************************************************************
!  Subroutine EMISS_VISTAS_ANTHRO reads the VISTAS emission 
!  fields at 1x1 resolution and regrids them to the current model resolution.
!  (amv, 11/24/2008)
!  
!  NOTES:
!  (1 )
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_O3"            ! FSCALYR

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR
      INTEGER                    :: MN, SNo, ScNo
      REAL*4                     :: ARRAY(I1x1,J1x1,1)
      REAL*8                     :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                     :: SC_1x1(I1x1,J1x1)
      REAL*8                     :: TAU2002, TAU
      CHARACTER(LEN=255)         :: FILENAME, VISTAS_DIR
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER(LEN=2)           :: SMN
      CHARACTER(LEN=1)           :: SSMN

      !=================================================================
      ! EMISS_VISTAS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_VISTAS_ANTHRO
         FIRST = .FALSE.
      ENDIF

      VISTAS_DIR = TRIM( DATA_DIR_1x1 ) // 'VISTAS_200811/'

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! cap maximum scaling year
      IF ( THISYEAR .gt. 2007 ) THEN
         THISYEAR = 2007
      ENDIF

      SNAME = 'NOx'
      SNo = 1
      ScNo = 71
            
      TAU2002 = GET_TAU0( 1, 1, 2002)
      MN = GET_MONTH()

      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SMN // '.1x1'
      ENDIF

      WRITE( SYEAR, '(i4)' ) THISYEAR


      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_VISTAS_ANTHRO: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Load ozone season regulation factors
      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SMN // '.1x1'
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'RATIO-2D', ScNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Apply Ozone Season Scalars
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * ARRAY(:,:,1)

      ! Apply Annual Scalar
      IF ( THISYEAR .ne. 2002 ) THEN
         CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2002,
     &                             THISYEAR, SC_1x1 )

         GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)
      ENDIF

      ! Load/Apply weekend/weekday factors
      TAU = GET_TAU0( MN, 1, 1999)
      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkend_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 2,
     &                    TAU,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1 * ARRAY, 
     &                        VISTAS_WE_NOx )

      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkday_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 1,
     &                    TAU,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1 * ARRAY, 
     &                        VISTAS_WD_NOx )

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL VISTAS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR, MN )

      ! Return to calling program
      END SUBROUTINE EMISS_VISTAS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE VISTAS_SCALE_FUTURE
!
!******************************************************************************
!  Subroutine VISTAS_SCALE_FUTURE applies the IPCC future scale factors to 
!  the VISTAS anthropogenic emissions.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff

#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      INTEGER                       :: I, J

      !=================================================================
      ! STREETS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         VISTAS_WE_NOx(I,J)  = VISTAS_WE_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )
         VISTAS_WD_NOx(I,J)  = VISTAS_WD_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE VISTAS_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG ( YEAR, THISMONTH )
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the totals for the anthropogenic
!  emissions of NOx. (bmy, phs, 3/11/08)
!
!  NOTES:
!******************************************************************************
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTNOX

#     include "CMN_SIZE"   ! Size parameters

      ! argument
      INTEGER, INTENT(IN) :: YEAR, THISMONTH


      ! Local variables
      INTEGER             :: I,     J
      REAL*8              :: WD_NOX, WE_NOX, F_NOX, A
      CHARACTER(LEN=3)    :: UNIT

     ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      WD_NOX  = 0d0
      WE_NOX  = 0d0
      F_NOX   = TRACER_MW_KG(IDTNOX )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Weekday avg emissions
            WD_NOX  = WD_NOX  + VISTAS_WD_NOX (I,J) * A * F_NOX

            ! Weekend avg emissions
            WE_NOX  = WE_NOX  + VISTAS_WE_NOX (I,J) * A * F_NOX

         ENDDO
      ENDDO

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'VISTAS   U S A   E M I S S I O N S', / )

      ! Weekday avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NOx ', THISMONTH, WD_NOX,  '  '
 110  FORMAT( 'Total weekday avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Weekend avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 120   ) 'NOx ', THISMONTH, WE_NOX,  '  '
 120  FORMAT( 'Total weekend avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg

!------------------------------------------------------------------------------

      SUBROUTINE INIT_VISTAS_ANTHRO
!
!******************************************************************************
!  Subroutine INIT_VISTAS_ANTHRO allocates and zeroes all module arrays.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LVISTAS

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_VISTAS_ANTHRO begins here!
      !=================================================================

      ! Return if LVISTAS is false
      IF ( .not. LVISTAS ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( VISTAS_WD_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VISTAS_WD_NOx' )
      VISTAS_WD_NOx = 0d0

      ALLOCATE( VISTAS_WE_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VISTAS_WE_NOx' )
      VISTAS_WE_NOx = 0d0

      !---------------------------------------------------
      ! Pre-store array for grid box surface area in cm2
      !---------------------------------------------------

      ! Allocate array
      ALLOCATE( A_CM2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_VISTAS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_VISTAS_ANTHRO
!
!******************************************************************************
!  Subroutine CLEANUP_VISTAS_ANTHRO deallocates all module arrays 
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_STREETS begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2          ) ) DEALLOCATE( A_CM2          )
      IF ( ALLOCATED( VISTAS_WD_NOx  ) ) DEALLOCATE( VISTAS_WD_NOx  )
      IF ( ALLOCATED( VISTAS_WE_NOx  ) ) DEALLOCATE( VISTAS_WE_NOx  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_VISTAS_ANTHRO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE VISTAS_ANTHRO_MOD
