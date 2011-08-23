!------------------------------------------------------------------------------
!     Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: icoads_ship_mod
!     
! !DESCRIPTION: Module ICOADS\_SHIP\_MOD contains variables and routines to 
!  read the International Comprehensive Ocean-Atmosphere Data Set (ICOADS)
!  ship emissions.  Base year is 2002.
!\\   
!\\   
! !INTERFACE: 
!     
      MODULE ICOADS_SHIP_MOD
! 
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_ICOADS_SHIP
      PUBLIC  :: EMISS_ICOADS_SHIP
      PUBLIC  :: GET_ICOADS_SHIP
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: ICOADS_SCALE_FUTURE
      PRIVATE :: INIT_ICOADS_SHIP
      PRIVATE :: TOTAL_ICOADS_SHIP_TG
!
! !REMARKS:
!  Source: ICOADS Emissions data for NOx, SOx, and CO were downloaded from 
!  http://coast.cms.udel.edu/GlobalShipEmissions/Inventories/
!
!  Reference: Wang, C., J. J. Corbett, and J. Firestone, \emph{Improving
!  Spatial representation of Global Ship Emissions Inventories},
!  Environ. Sci. Technol., 42, (1), 193-199, 2008.
!
! !REVISION HISTORY:
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
      ! Array for surface area
      REAL*8,  ALLOCATABLE :: A_CM2(:)

      ! Arrays for emissions
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_icoads_ship
!
! !DESCRIPTION: Function GET\_ICOADS\_SHIP returns the ICOADS ship emissions 
!  for GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in 
!  units of [kg/s] or [molec/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ICOADS_SHIP( I,    J,     N, 
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNOx, IDTCO, IDTSO2, IDTNH3
      USE TIME_MOD,     ONLY : GET_YEAR, GET_MONTH
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N

      ! OPTIONAL -- return emissions in [molec/cm2/s]
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S  

      ! OPTIONAL -- return emissions in [kg/s]
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
!
! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE     
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS, DO_MCS
      INTEGER                       :: YEAR, MONTH
      REAL*8                        :: SEC_IN_MONTH

      !=================================================================
      ! GET_ICOADS_SHIP begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      IF ( N == IDTNOx ) THEN

         ! NOx [kg/month]
         VALUE = NOx(I,J)

      ELSE IF ( N == IDTCO ) THEN

         ! CO [kg/month]
         VALUE = CO(I,J)

      ELSE IF ( N == IDTSO2 ) THEN

         ! SO2 [kg/month]
         VALUE = SO2(I,J)

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no CAC emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      ! Get emissions year
      YEAR = GET_YEAR()

      ! Get emissions month      
      MONTH = GET_MONTH()

      IF ( (MONTH == 4) .OR. (MONTH == 6) .OR.
     &   (MONTH == 9) .OR. (MONTH == 11) ) THEN

         SEC_IN_MONTH = 86400D0*30.0D0

      ELSE IF (MONTH == 2) THEN

         ! ICOADS ship emissions for 2002
         IF (MOD(YEAR,4) == 0) THEN
            SEC_IN_MONTH = 86400D0*29.0D0
         ELSE
            SEC_IN_MONTH = 86400D0*28.0D0
         ENDIF

      ELSE

         SEC_IN_MONTH = 86400D0*31.0D0

      ENDIF

      IF ( DO_KGS ) THEN
            
         ! Convert from [kg/box/month] to [kg/box/s]
         VALUE = VALUE / SEC_IN_MONTH

      ELSE IF ( DO_MCS ) THEN

         ! Convert NOx from [kg/month] to [molec/cm2/s]
         VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_MONTH )

      ENDIF

      END FUNCTION GET_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_icoads_ship
!
! !DESCRIPTION: Subroutine EMISS\_ICOADS\_SHIP reads the ICOADS emission fields
!  at 1x1 resolution and regrids them to the current model resolution. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_ICOADS_SHIP
!
! !USES:
! 
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR,      GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1

      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR, SPECIES, SNo, ScNo
      INTEGER                    :: THISMONTH
      REAL*4                     :: ARRAY(I1x1,J1x1,1)
      REAL*8                     :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                     :: SC_1x1(I1x1,J1x1)
      REAL*8                     :: TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER (LEN=2)          :: SMONTH

      !=================================================================
      ! EMISS_ICOADS_SHIP begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_ICOADS_SHIP
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! Get emissions month      
      THISMONTH = GET_MONTH()

      WRITE( SMONTH, '(i2.2)' ) THISMONTH


      DO SPECIES = 1,3

         IF ( SPECIES .eq. 1 ) THEN
            SNAME = 'NOx'
            SNo = 1
            ScNo = 71
         ELSEIF ( SPECIES .eq. 2 ) THEN
            SNAME = 'CO'
            SNo = 4
            ScNo = 72
         ELSEIF ( SPECIES .eq. 3 ) THEN
            SNAME = 'SOx'
            SNo = 26
            ScNo = 73
         ENDIF
            

         ! TAU values for 2002
         TAU = GET_TAU0( 1, 1, 2002 )

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) //'ICOADS_200907/' //
     &               TRIM( SNAME ) // '_' // SMONTH // '.geos.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - EMISS_ICOADS_SHIP: Reading ', a )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'ICOADS-$', SNo, 
     &                    TAU,      I1x1,       J1x1,     
     &                    1,        ARRAY,      QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEOS_1x1(:,:,1) = ARRAY(:,:,1)

         ! Convert [kg S/month] to [kg SO2/month]
         IF ( SPECIES .eq. 3 ) THEN
            GEOS_1X1 = GEOS_1x1*64.0D0/32.0D0
         ENDIF

         ! Apply annual scalar factor
         CALL GET_ANNUAL_SCALAR_1x1( ScNo, 2002, THISYEAR, SC_1x1 )

         GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)



         ! Regrid from GEOS 1x1 --> current model resolution
         IF ( SPECIES .eq. 1 ) THEN

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, NOx )

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, CO )

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * 0.986

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, SO2 )

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL ICOADS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ICOADS_SHIP_TG( THISYEAR )

      END SUBROUTINE EMISS_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: icoads_scale_future
!
! !DESCRIPTION: applies the IPCC future scale factors
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE ICOADS_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

      USE CMN_SIZE_MOD             ! Size parameters
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! ICOADS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/month]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /month]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/month] 
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE ICOADS_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_icoads_ship_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ICOADS\_SHIP\_TG prints the totals for  
!   ship emissions of NOx, CO, and SO2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ICOADS_SHIP_TG( MONTH )
!
! !USES:
! 

      USE CMN_SIZE_MOD            ! Size parameters

!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MONTH   ! Month of data to compute totals
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      REAL*8              :: T_NOX, T_CO,  T_SO2
      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ICOADS_SHIP_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'I. C. O. A. D. S.   S H I P   E M I S S I O N S', / )


      ! Total NOx [Tg N]
      T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )

      ! Total CO  [Tg CO]
      T_CO  = SUM( CO  ) * 1d-9

      ! Total SO2 [Tg S]
      T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

      ! Print totals in [kg]
      WRITE( 6, 110 ) 'NOx ', MONTH, T_NOx,  '[Tg N  ]'
      WRITE( 6, 110 ) 'CO  ', MONTH, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'SO2 ', MONTH, T_SO2,  '[Tg S  ]'

      ! Format statement
 110  FORMAT( 'ICOADS ship ', a5, 
     &        'for month ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      END SUBROUTINE TOTAL_ICOADS_SHIP_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_icoads_ship
!
! !DESCRIPTION: Subroutine INIT\_ICOADS\_SHIP allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ICOADS_SHIP
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LICOADSSHIP

      USE CMN_SIZE_MOD    ! Size parameters
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_ICOADS_SHIP begins here!
      !=================================================================

      ! Return if LICOADSSHIP is false
      IF ( .not. LICOADSSHIP ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

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

      END SUBROUTINE INIT_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_icoads_ship
!
! !DESCRIPTION:  Subroutine CLEANUP\_ICOADS\_SHIP deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ICOADS_SHIP
!
! !REVISION HISTORY: 
!  21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( A_CM2          ) ) DEALLOCATE( A_CM2          )
      IF ( ALLOCATED( NOx            ) ) DEALLOCATE( NOx            )
      IF ( ALLOCATED( CO             ) ) DEALLOCATE( CO             )
      IF ( ALLOCATED( SO2            ) ) DEALLOCATE( SO2            )

      END SUBROUTINE CLEANUP_ICOADS_SHIP
!EOC
      END MODULE ICOADS_SHIP_MOD
