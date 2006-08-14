! $Id: streets_anthro_mod.f,v 1.1 2006/08/14 17:58:14 bmy Exp $
      MODULE STREETS_ANTHRO_MOD
!
!******************************************************************************
!  Module STREETS_ANTHRO_MOD contains variables and routines to read the 
!  David Streets et al Asian anthropogenic emissions for NOx and CO. 
!  (bmy, 8/8/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) A_CM2 (REAL*8)         : Array for grid box surface area [cm2]
!  (2 ) MASK  (REAL*8)         : Mask for the China/SE Asia region
!  (3 ) NOx   (REAL*8)         : Streets anthro NOx emissions [kg/yr]
!  (4 ) CO    (REAL*8)         : Streets anthro CO  emissions [kg/yr]
! 
!  Module Routines:
!  ============================================================================
!  (1 ) GET_STREETS_MASK       : Gets the China/SE Asia mask value at (I,J) 
!  (2 ) GET_STREETS_ANTHRO     : Gets emissions at (I,J) for emissions species 
!  (3 ) EMISS_STREETS_ANTHRO   : Reads STREETS emissions from disk
!  (4 ) STREETS_SCALE_FUTURE   : Applies IPCC future scale factors to emissions
!  (5 ) INIT_STREETS           : Allocates and zeroes module arrays
!  (6 ) CLEANUP_STREETS        : Dealocates module arrays
!
!  GEOS-Chem modules referenced by "streets_anthro_mod.f"
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
!  (1 ) Streets, D.G, Q. Zhang, L. Wang, K. He, J. Hao, Y. Wu, Y. Tang,
!        and G.C. Carmichael, "Revisiting China's CO emissions after the
!        Transport and Chemical Evolution over the Pacific (TRACE-P) mission:
!        Synthesis of inventories, atmospheric modeling, and observations",
!        J. Geophys. Res, 111, D14306, doi:10.1029/2006JD007118, 2006.
!  (2 ) Streets, D.G., T.C. Bond, G.R. Carmichael, S.D. Fernandes, Q. Fu,
!        Z. Klimont, S.M. Nelson, N.Y. Tsai, M.Q. Wang, J-H. Woo, and
!        K.F. Yarber, "An inventory of gaseous and primary aerosol emissions
!        in Asia in the year 2000", J. Geophys. Res, 108, D21, 
!        doi:10.1029/2002JD003093, 2003.
!  
!  NOTES: 
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "streets_anthro_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_STREETS_ANTHRO
      PUBLIC :: EMISS_STREETS_ANTHRO
      PUBLIC :: GET_STREETS_MASK
      PUBLIC :: GET_STREETS_ANTHRO

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Arrays
      REAL*8,  ALLOCATABLE :: A_CM2(:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: MASK(:,:)
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)

      ! Parameters
      INTEGER, PARAMETER   :: STREETS_YEAR = 2001
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_STREETS_MASK( I, J ) RESULT( THISMASK )
!
!******************************************************************************
!  Function GET_STREETS_MASK returns the value of the China mask for the David
!  Streets et al emissions at grid box (I,J).  MASK=1 if (I,J) is China, or 
!  MASK=0 otherwise. (bmy, 8/8/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-Chem longitude index 
!  (2 ) J (INTEGER) : GEOS-Chem latitude  index 
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      REAL*8              :: THISMASK

      !=================================================================
      ! GET_STREETS_MASK begins here!
      !=================================================================
      THISMASK = MASK(I,J)

      ! Return to calling program
      END FUNCTION GET_STREETS_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_STREETS_ANTHRO( I,    J,     N, 
     &                             MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_STREETS_ANTHRO returns the David Streets et al emission for 
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].  (bmy, 8/8/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude index
!  (3 ) N           (INTEGER) : GEOS-Chem tracer number
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL -- return emissions in [molec/cm2/s]
!  (5 ) KG_S        (LOGICAL) : OPTIONAL -- return emissions in [kg/s]
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNOx, IDTCO

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J, N
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
      
      ! Local variables
      LOGICAL                       :: DO_KGS, DO_MCS
      REAL*8                        :: VALUE

      !=================================================================
      ! GET_STREETS_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      !-----------------
      ! Get emissions
      !-----------------
      IF ( N == IDTNOx ) THEN

         ! NOx [kg/yr]
         VALUE = NOx(I,J)

      ELSE IF ( N == IDTCO ) THEN

         ! CO [kg/yr]
         VALUE = CO(I,J)

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no STREETS emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !-----------------
      ! Convert units
      !-----------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [kg/yr] to [kg/s]
         VALUE = VALUE / SEC_IN_YEAR

      ELSE IF ( DO_MCS ) THEN

         ! Convert NOx from [kg/yr] to [molec/cm2/s]
         VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_YEAR )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine EMISS_STREETS_ANTHRO reads the David Streets et al emission 
!  fields at 1x1 resolution and regrids them to the current model resolution.
!  (bmy, 8/8/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,    ONLY : LFUTURE
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: TAU0
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1,1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! EMISS_STREETS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_STREETS_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! TAU0 value
      TAU0 = GET_TAU0( 1, 1, STREETS_YEAR )
         
!------------------------------------------------------------------------------
! Prior to 8/8/06:
! Uncomment when we install NOx emissions (bmy, 8/8/06)
!      !--------------------------
!      ! Read NOx and regrid
!      !--------------------------
!
!      ! File name
!      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
!     &            'Streets_200607/Streets_NOx.generic.1x1'
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )
!
!      ! Read data [unitless]
!      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, 
!     &                 TAU0,      I1x1,      J1x1-1,     
!     &                 1,         ARRAY,     QUIET=.TRUE. ) 
!
!      ! Cast to REAL*8 before regridding
!      GEN_1x1(:,:,1) = ARRAY(:,:,1)
!
!      ! Regrid from GENERIC 1x1 --> GEOS 1x1 
!      CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1(:,:,1), GEOS_1x1(:,:,1) )
!
!      ! Regrid from GEOS 1x1 --> current model resolution
!      CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, STREETS_NOx )
!------------------------------------------------------------------------------
   
      !--------------------------
      ! Read CO and regrid
      !--------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'Streets_200607/Streets_CO_FF_2001.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:,1) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1(:,:,1), GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, CO )

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL STREETS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( STREETS_YEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE STREETS_SCALE_FUTURE
!
!******************************************************************************
!  Subroutine STREETS_SCALE_FUTURE applies the IPCC future scale factors to 
!  the STREETS anthropogenic emissions. (swu, bmy, 8/8/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
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

         ! Future NOx [molec/cm2/s]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO [molec/cm2/s]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE STREETS_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG( YEAR )
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the totals for the anthropogenic
!  emissions of NOx and CO. (bmy, 8/8/06)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR (INTGER) : Year of the emissions data 
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)   :: YEAR

      ! Local variables
      INTEGER               :: I,     J
      REAL*8                :: T_NOX, T_CO
      CHARACTER(LEN=3)      :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'S T R E E T S   A S I A N   E M I S S I O N S', / )

      !----------------
      ! Sum emissions
      !----------------
      
      ! Total NOx [Tg]
      T_NOX = SUM( NOx ) * 1d-9

      ! Total CO [Tg]
      T_CO  = SUM( CO  ) * 1d-9

      !----------------
      ! Print sums
      !----------------

      ! Print totals in [kg]
      WRITE( 6, 110   ) 'NOx ', YEAR, T_NOx,  ' N '
      WRITE( 6, 110   ) 'CO  ', YEAR, T_CO,   ' CO'
 110  FORMAT( 'David Streets anthropogenic ', a4, ' for ', i4, 
     &        ': ', f11.4, ' Tg', a3 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg

!------------------------------------------------------------------------------

      SUBROUTINE READ_CHINA_MASK
!
!******************************************************************************
!  Subroutine READ_CHINA_MASK allocates and initialize the array that 
      ! contains the mask for the China region
!******************************************************************************
!
      !--------------------------------------------------
      ! Allocate and initialize the array that 
      ! contains the mask for the China region
      !--------------------------------------------------


      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'Streets_200607/China_mask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - INIT_STREETS_ANTHRO: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,       I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:,1) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( GEN_1x1(:,:,1), 
     &                        GEOS_1x1(:,:,1),
     &                        PER_UNIT_AREA=.TRUE. )

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK )


!------------------------------------------------------------------------------

      SUBROUTINE INIT_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine INIT_STREETS_ANTHRO allocates and zeroes all module arrays.
!  (bmy, 8/8/06) 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE ERROR_MOD,      ONLY : ALLOC_ERR
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LSTREETS
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER                 :: AS, J
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1,1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! INIT_STREETS begins here!
      !=================================================================

      ! Return if LSTREETS is false
      IF ( .not. LSTREETS ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

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
      END SUBROUTINE INIT_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine CLEANUP_STREETS deallocates all module arrays (bmy, 8/8/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_STREETS begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2 ) ) DEALLOCATE( A_CM2 )
      IF ( ALLOCATED( CO    ) ) DEALLOCATE( CO    )
      IF ( ALLOCATED( MASK  ) ) DEALLOCATE( MASK  )
      IF ( ALLOCATED( NOx   ) ) DEALLOCATE( NOx   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_STREETS_ANTHRO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE STREETS_ANTHRO_MOD
