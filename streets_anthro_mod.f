! $Id: streets_anthro_mod.f,v 1.2 2006/08/18 20:32:39 bmy Exp $
      MODULE STREETS_ANTHRO_MOD
!
!******************************************************************************
!  Module STREETS_ANTHRO_MOD contains variables and routines to read the 
!  David Streets et al Asian anthropogenic emissions for NOx and CO. 
!  (yxw, bmy, 8/16/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) A_CM2        (REAL*8)  : Array for grid box surface area [cm2]
!  (2 ) MASK_CHINA   (REAL*8)  : Mask for the China region (for 2001 CO)
!  (3 ) MASK_SE_ASIA (REAL*8)  : Mask for the SE Asia region (for 2000 emiss.)
!  (4 ) NOx          (REAL*8)  : Streets anthro NOx emissions [kg/yr]
!  (5 ) CO           (REAL*8)  : Streets anthro CO  emissions [kg/yr]
!  (6 ) SO2          (REAL*8)  : Streets anthro SO2 emissions [kg/yr]
!  (7 ) NH3          (REAL*8)  : Streets anthro NH3 emissions [kg/yr]          
!  (8 ) CO2          (REAL*8)  : Streets anthro CO2 emissions [kg/yr]
!  (9 ) CH4          (REAL*8)  : Streets anthro CH4 emissions [kg/yr]
! 
!  Module Routines:
!  ============================================================================
!  (1 ) GET_CHINA_MASK         : Gets the China mask value at (I,J) 
!  (2 ) GET_SE_ASIA_MASK       : Gets the SE Asia mask value at (I,J) 
!  (3 ) GET_STREETS_ANTHRO     : Gets emissions at (I,J) for emissions species 
!  (4 ) EMISS_STREETS_ANTHRO   : Reads Streets' emissions from disk
!  (5 ) STREETS_SCALE_FUTURE   : Applies IPCC future scale factors to emissions
!  (6 ) READ_STREETS_MASKS     : Reads mask info from disk
!  (7 ) INIT_STREETS_ANTHRO    : Allocates and zeroes module arrays
!  (8 ) CLEANUP_STREETS_ANTHRO : Dealocates module arrays
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
      PUBLIC :: GET_CHINA_MASK
      PUBLIC :: GET_SE_ASIA_MASK
      PUBLIC :: GET_STREETS_ANTHRO

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Arrays
      REAL*8,  ALLOCATABLE :: A_CM2(:)
      REAL*8,  ALLOCATABLE :: MASK_CHINA(:,:)
      REAL*8,  ALLOCATABLE :: MASK_SE_ASIA(:,:)
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)
      REAL*8,  ALLOCATABLE :: NH3(:,:)
      REAL*8,  ALLOCATABLE :: CO2(:,:)
      REAL*8,  ALLOCATABLE :: CH4(:,:)

      ! Parameters
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_CHINA_MASK( I, J ) RESULT( THISMASK )
!
!******************************************************************************
!  Function GET_STREETS_MASK returns the value of the China mask for the David
!  Streets et al emissions at grid box (I,J).  MASK=1 if (I,J) is China, or 
!  MASK=0 otherwise. (bmy, 8/16/06)
!
!  NOTE: The China Mask is used with the 2001 CO emissions.
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
      ! GET_CHINA_MASK begins here!
      !=================================================================
      THISMASK = MASK_CHINA(I,J)

      ! Return to calling program
      END FUNCTION GET_CHINA_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_SE_ASIA_MASK( I, J ) RESULT( THISMASK )
!
!******************************************************************************
!  Function GET_SE_ASIA_MASK returns the value of the China mask for the David
!  Streets et al emissions at grid box (I,J).  MASK=1 if (I,J) is China, or 
!  MASK=0 otherwise. (bmy, 8/16/06)
!
!  NOTE: The SE Asia Mask is used with the 2000 emissions for 
!        NOx, CO, CO2, SO2, NH3, and CH4. 
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
      ! GET_SE_ASIA_MASK begins here!
      !=================================================================
      THISMASK = MASK_SE_ASIA(I,J)

      ! Return to calling program
      END FUNCTION GET_SE_ASIA_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_STREETS_ANTHRO( I,    J,     N, 
     &                             MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
!******************************************************************************
!  Function GET_STREETS_ANTHRO returns the David Streets et al emission for 
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].  (bmy, 8/16/06)
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
      USE TRACER_MOD,           ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,           ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,           ONLY : XNUMOL
      USE TRACERID_MOD,         ONLY : IDTNOx, IDTCO, IDTSO2, IDTNH3

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

      ! Test for simulation type
      IF ( ITS_A_CH4_SIM() ) THEN

         !-------------------
         ! CH4 simulation
         !-------------------
         VALUE = CH4(I,J)

      ELSE IF ( ITS_A_CO2_SIM() ) THEN
         
         !-------------------
         ! CH4 simulation
         !-------------------
         VALUE = CO2(I,J)

      ELSE

         !-------------------
         ! Other simulations
         !-------------------
         IF ( N == IDTNOx ) THEN

            ! NOx [kg/yr]
            VALUE = NOx(I,J)

         ELSE IF ( N == IDTCO ) THEN

            ! CO [kg/yr]
            VALUE = CO(I,J)

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2(I,J)

         ELSE IF ( N == IDTNH3 ) THEN

            ! NH3 [kg/yr]
            VALUE = NH3(I,J)

         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no STREETS emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
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
!  (bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,    ONLY : LFUTURE
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1, DO_REGRID_G2G_1x1
      USE TRACER_MOD,     ONLY : ITS_A_CO2_SIM, ITS_A_CH4_SIM

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                  :: TAU2000, TAU2001
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! EMISS_STREETS_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_STREETS_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! TAU0 values for 2000 and 2001
      TAU2000 = GET_TAU0( 1, 1, 2000 )
      TAU2001 = GET_TAU0( 1, 1, 2001 )

      ! Test for simulation type
      IF ( ITS_A_CH4_SIM() ) THEN

         !--------------------------
         ! Read CH4 and regrid
         ! (CH4 simulations only)
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_CH4_FF_2000.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - EMISS_STREETS_ANTHRO: Reading ', a )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 1, 
     &                    TAU2000,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 --> GEOS 1x1 
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 --> current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, CH4 )

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !--------------------------
         ! Read CO2 and regrid
         ! (CO2 simulations only)
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_CO2_FF_2000.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                    TAU2000,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 --> GEOS 1x1 
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 --> current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, CO2 )

      ELSE

         !--------------------------
         ! Read NOx and regrid
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_NOx_FF_2000.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, 
     &                    TAU2000,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 --> GEOS 1x1 
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 --> current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, NOx )

         !--------------------------
         ! Read CO and regrid
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_CO_FF_2001.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, 
     &                    TAU2001,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 GRID to current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, CO )

         !--------------------------
         ! Read SO2 and regrid
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_SO2_FF_2000.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 26, 
     &                    TAU2000,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 GRID to current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, SO2 )

         !--------------------------
         ! Read NH3 and regrid
         !--------------------------

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &               'Streets_200607/Streets_NH3_FF_2000.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data [unitless]
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 30, 
     &                    TAU2000,   I1x1,      J1x1-1,     
     &                    1,         ARRAY,     QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEN_1x1(:,:) = ARRAY(:,:,1)

         ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
         CALL DO_REGRID_G2G_1x1( 'kg/yr', GEN_1x1, GEOS_1x1(:,:,1) )

         ! Regrid from GEOS 1x1 GRID to current model resolution
         CALL DO_REGRID_1x1( 'kg/yr', GEOS_1x1, NH3 )

      ENDIF

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL STREETS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg

      ! Return to calling program
      END SUBROUTINE EMISS_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE STREETS_SCALE_FUTURE
!
!******************************************************************************
!  Subroutine STREETS_SCALE_FUTURE applies the IPCC future scale factors to 
!  the David Streets' anthropogenic emissions. (swu, bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

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
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /yr]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/yr] 
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future SO2 [kg SO2/yr] 
         NH3(I,J)  = NH3(I,J) * GET_FUTURE_SCALE_NH3an( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE STREETS_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the totals for the anthropogenic
!  emissions of NOx and CO. (bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : ITS_A_CH4_SIM, ITS_A_CO2_SIM 

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER             :: I,     J
      REAL*8              :: T_NOX, T_CO,  T_SO2
      REAL*8              :: T_NH3, T_CH4, T_CO2
      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'S T R E E T S   A S I A N   E M I S S I O N S', / )

      ! Test for simulation type
      IF ( ITS_A_CH4_SIM() ) THEN

         !-----------------------
         ! CH4 simulation
         !-----------------------

         ! Total CH4 [Tg CH4]
         T_CH4 = SUM( CH4 ) * 1d-9

         ! Print totals
         WRITE( 6, 110 ) 'CH4 ', 2000, T_NOx,  ' CH4'

      ELSE IF ( ITS_A_CO2_SIM() ) THEN

         !-----------------------
         ! CO2 simulation
         !-----------------------

         ! Total CO2 [Tg CO2]
         T_CH4 = SUM( CO2 ) * 1d-9

         ! Print totals
         WRITE( 6, 110 ) 'CO2 ', 2000, T_NOx,  ' CO2'

      ELSE

         !-----------------------
         ! Other simulations
         !-----------------------

         ! Total NOx [Tg N]
         T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )
 
         ! Total CO  [Tg CO]
         T_CO  = SUM( CO  ) * 1d-9

         ! Total SO2 [Tg S]
         T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

         ! Total NH3 [Tg NH3]
         T_NH3 = SUM( NH3 ) * 1d-9

         ! Print totals in [kg]
         WRITE( 6, 110 ) 'NOx ', 2000, T_NOx, '[Tg N  ]'
         WRITE( 6, 110 ) 'CO  ', 2001, T_CO,  '[Tg CO ]'
         WRITE( 6, 110 ) 'SO2 ', 2000, T_SO2, '[Tg S  ]'
         WRITE( 6, 110 ) 'NH3 ', 2000, T_NH3, '[Tg NH3]'

      ENDIF

      ! Format statement
 110  FORMAT( 'David Streets anthro ', a4, 
     &        'for base year ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg

!------------------------------------------------------------------------------

      SUBROUTINE READ_STREETS_MASKS
!
!******************************************************************************
!  Subroutine READ_STREETS_MASKS reads and regrids the China and SE Asia
!  masks that define the David Streets' emission regions (bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_G2G_1x1, DO_REGRID_1x1

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_STREETS_MASKS begins here!
      !=================================================================

      !------------------------------------
      ! China Mask (for 2001 CO emisisons)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'Streets_200607/China_mask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_STREETS_MASKS: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,       I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK_CHINA )

      !------------------------------------
      ! SE Asia Mask (for 2000 emissions)
      !------------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'Streets_200607/SE_Asia_mask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,       I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, MASK_SE_ASIA )

      ! Return to calling program
      END SUBROUTINE READ_STREETS_MASKS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine INIT_STREETS_ANTHRO allocates and zeroes all module arrays.
!  (bmy, 8/16/06) 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LSTREETS

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER              :: AS, J

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

      ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0

      ALLOCATE( CO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO2' )
      CO2 = 0d0

      ALLOCATE( CH4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4' )
      CH4 = 0d0

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

      !---------------------------------------------------
      ! Read & Regrid masks for Streets' emissions
      !---------------------------------------------------

      ALLOCATE( MASK_CHINA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CHINA' )
      MASK_CHINA = 0d0

      ALLOCATE( MASK_SE_ASIA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_SE_ASIA' )
      MASK_SE_ASIA = 0d0

      ! Read China & SE Asia masks from disk
      CALL READ_STREETS_MASKS

      ! Return to calling program
      END SUBROUTINE INIT_STREETS_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_STREETS_ANTHRO
!
!******************************************************************************
!  Subroutine CLEANUP_STREETS deallocates all module arrays (bmy, 8/16/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_STREETS begins here!
      !=================================================================
      IF ( ALLOCATED( A_CM2        ) ) DEALLOCATE( A_CM2        )
      IF ( ALLOCATED( MASK_CHINA   ) ) DEALLOCATE( MASK_CHINA   ) 
      IF ( ALLOCATED( MASK_SE_ASIA ) ) DEALLOCATE( MASK_SE_ASIA )
      IF ( ALLOCATED( NOx          ) ) DEALLOCATE( NOx          )
      IF ( ALLOCATED( CO           ) ) DEALLOCATE( CO           )
      IF ( ALLOCATED( SO2          ) ) DEALLOCATE( SO2          )
      IF ( ALLOCATED( NH3          ) ) DEALLOCATE( NH3          )
      IF ( ALLOCATED( CH4          ) ) DEALLOCATE( CH4          )
      IF ( ALLOCATED( CO2          ) ) DEALLOCATE( CO2          )

      ! Return to calling program
      END SUBROUTINE CLEANUP_STREETS_ANTHRO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE STREETS_ANTHRO_MOD
