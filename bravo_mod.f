! $Id: bravo_mod.f,v 1.3 2006/09/08 19:20:51 bmy Exp $
      MODULE BRAVO_MOD
!
!******************************************************************************
!  Module BRAVO_MOD contains variables and routines to read the BRAVO 
!  Mexican anthropogenic emission inventory for NOx, CO, and SO2. 
!  (rjp, kfb, bmy, 6/22/06, 8/9/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BRAVO_CO   (REAL*8)    : BRAVO anthro CO  emissions [molec/cm2/s]
!  (2 ) BRAVO_MASK (REAL*8)    : Array used to mask out the Europe region
!  (3 ) BRAVO_NOx  (REAL*8)    : BRAVO anthro NOx emissions [molec/cm2/s]
!  (4 ) BRAVO_SO2  (REAL*8)    : BRAVO anthro SO2 emissions [molec/cm2/s]
! 
!  Module Routines:
!  ============================================================================
!  (1 ) GET_BRAVO_MASK         : Gets the value of the Mexico mask at (I,J) 
!  (2 ) GET_BRAVO_ANTHRO       : Gets emissions at (I,J) for BRAVO species 
!  (3 ) EMISS_BRAVO            : Reads BRAVO emissions from disk once per year
!  (4 ) BRAVO_SCALE_FUTURE     : Applies IPCC future scale factors to BRAVO
!  (5 ) INIT_BRAVO             : Allocates and zeroes module arrays
!  (6 ) CLEANUP_BRAVO          : Dealocates module arrays
!
!  GEOS-CHEM modules referenced by "bravo_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f        : Module w/ GEOS-Chem data & met field dirs
!  (3 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (5 ) future_emissions_mod.f : Module w/ routines for IPCC future emissions
!  (6 ) grid_mod.f             : Module w/ horizontal grid information
!  (7 ) logical_mod.f          : Module w/ GEOS-Chem logical switches
!  (8 ) regrid_1x1_mod.f       : Module w/ routines to regrid 1x1 data  
!  (9 ) time_mod.f             : Module w/ routines for computing time & date
!  (10) tracerid_mod.f         : Module w/ pointers to tracers & emissions  
!
!  References:
!  ============================================================================
!  (1 ) Kuhns, H., M. Green, and Etyemezian, V, "Big Bend Regional Aerosol and
!        Visibility Observational (BRAVO) Study Emissions Inventory", Desert
!        Research Institute, 2003.
!
!  NOTES: 
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "bravo_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_BRAVO
      PUBLIC :: EMISS_BRAVO
      PUBLIC :: GET_BRAVO_MASK
      PUBLIC :: GET_BRAVO_ANTHRO

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Arrays
      REAL*8,  ALLOCATABLE :: BRAVO_MASK(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_NOx(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_CO(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_SO2(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION GET_BRAVO_MASK( I, J ) RESULT( MASK )
!
!******************************************************************************
!  Function GET_BRAVO_MASK returns the value of the Mexico mask for BRAVO
!  emissions at grid box (I,J).  MASK=1 if (I,J) is in the BRAVO Mexican
!  region, or MASK=0 otherwise. (rjp, kfb, bmy, 6/22/06)
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

      ! Function return value
      REAL*8              :: MASK

      !=================================================================
      ! GET_BRAVO_MASK begins here!
      !=================================================================
      MASK = BRAVO_MASK(I,J)

      ! Return to calling program
      END FUNCTION GET_BRAVO_MASK

!------------------------------------------------------------------------------

      FUNCTION GET_BRAVO_ANTHRO( I, J, N ) RESULT( BRAVO )
!
!******************************************************************************
!  Function GET_BRAVO_ANTHRO returns the BRAVO emission for GEOS-Chem grid box
!  (I,J) and tracer N.  Units are [molec/cm2/s]. (rjp, kfb, bmy, 6/22/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-Chem longitude index
!  (2 ) J (INTEGER) : GEOS-Chem latitude  index
!  (3 ) N (INTEGER) : GEOS-Chem tracer    number
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTNOX, IDTCO, IDTSO2

      ! Arguments
      INTEGER, INTENT(IN)   :: I, J, N

      ! Function return value
      REAL*8                :: BRAVO
      
      !=================================================================
      ! GET_BRAVO_ANTHRO begins here!
      !=================================================================

      ! NOx
      IF ( N  == IDTNOX ) THEN
         BRAVO = BRAVO_NOx(I,J)

      ! CO
      ELSE IF ( N == IDTCO ) THEN
         BRAVO = BRAVO_CO(I,J)

      ! SO2 
      ELSE IF ( N == IDTSO2 ) THEN
         BRAVO = BRAVO_SO2(I,J)

      ! Otherwise return a negative value to indicate
      ! that there are no BRAVO emissions for tracer N
      ELSE
         BRAVO = -1d0

      ENDIF

      ! Return to calling program
      END FUNCTION GET_BRAVO_ANTHRO

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_BRAVO
!
!******************************************************************************
!  Subroutine EMISS_BRAVO reads the BRAVO emission fields at 1x1 
!  resolution and regrids them to the current model resolution. 
!  (rjp, kfb, bmy, 6/22/06, 8/9/06)
!
!  NOTES:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
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
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                  :: TAU0
      CHARACTER(LEN=255)      :: FILENAME
      
      !=================================================================
      ! EMISS_BRAVO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_BRAVO
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read data from disk
      !=================================================================

      ! Use 1999 for BRAVO emission files
      TAU0  = GET_TAU0( 1, 1, 1999 )
        
      !---------------------
      ! Read and regrid NOx
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.NOx.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_BRAVO: Reading ', a )
      
      ! Read NOx [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid NOx [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid NOx [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_NOx )

      !---------------------
      ! Read and regrid CO
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.CO.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read CO [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid CO [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid CO [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_CO )

      !---------------------
      ! Read and regrid SO2
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.SO2.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read SO2 [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 26, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid SO2 [molec/cm2/s] to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'molec/cm2/s', GEN_1x1, GEOS_1x1(:,:,1) )

      ! Regrid SO2 [molec/cm2/s] to current model resolution
      CALL DO_REGRID_1x1( 'molec/cm2/s', GEOS_1x1, BRAVO_SO2 )

      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================
      IF ( LFUTURE ) THEN 
         CALL BRAVO_SCALE_FUTURE
      ENDIF

      !=================================================================
      ! Print emission totals
      !=================================================================
      CALL TOTAL_ANTHRO_TG

      ! Return to calling program
      END SUBROUTINE EMISS_BRAVO

!------------------------------------------------------------------------------

      SUBROUTINE BRAVO_SCALE_FUTURE
!
!******************************************************************************
!  Subroutine BRAVO_SCALE_FUTURE applies the IPCC future scale factors to 
!  the BRAVO anthropogenic emissions. (swu, bmy, 5/30/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      INTEGER                       :: I, J

      !=================================================================
      ! BRAVO_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [molec/cm2/s]
         BRAVO_NOx(I,J) = BRAVO_NOx(I,J)                * 
     &                    GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO [molec/cm2/s]
         BRAVO_CO(I,J)  = BRAVO_CO(I,J)                 *
     &                    GET_FUTURE_SCALE_COff( I, J )

         ! Future ALK4 [atoms C/cm2/s]
         BRAVO_SO2(I,J) = BRAVO_SO2(I,J)                *
     &                    GET_FUTURE_SCALE_SO2ff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE BRAVO_SCALE_FUTURE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG
!
!******************************************************************************
!  Subroutine TOTAL_ANTHRO_TG prints the amount of BRAVO anthropogenic 
!  emissions that are emitted each month,(rjp, kfb, bmy, 6/26/06)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTNOX, IDTCO, IDTSO2

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER               :: I, J
      REAL*8                :: A, B(3), NOx, CO, SO2
      CHARACTER(LEN=3)      :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'B R A V O   M E X I C A N   E M I S S I O N S', / )
      
      !----------------
      ! Sum emissions
      !----------------
      
      ! Define conversion factors for kg/molec
      ! (Undefined tracers will be zero)
      B(:) = 0d0
      IF ( IDTNOx > 0 ) B(1) = 1d0 / ( 6.0225d23 / 14d-3 )  ! Tg N
      IF ( IDTCO  > 0 ) B(2) = 1d0 / ( 6.0225d23 / 28d-3 )  ! Tg CO
      IF ( IDTSO2 > 0 ) B(3) = 1d0 / ( 6.0225d23 / 32d-3 )  ! Tg S

      ! Summing variables
      NOX = 0d0   
      CO  = 0d0 
      SO2 = 0d0 

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Convert [molec/cm2/s] to [Tg]
         ! (Multiply by 1d-9 to convert from [kg] to [Tg])
         A = GET_AREA_CM2( J ) * 365.25d0 * 86400d0 * 1d-9 
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Sum emissions (list NOx as Tg N)
            NOX = NOX + ( BRAVO_NOX(I,J) * A * B(1) )
            CO  = CO  + ( BRAVO_CO (I,J) * A * B(2) )
            SO2 = SO2 + ( BRAVO_SO2(I,J) * A * B(3) )
         ENDDO
      ENDDO
 
      !----------------
      ! Print sums
      !----------------

      ! Print totals in [kg]
      WRITE( 6, 110   ) 'NOx', NOx, ' N'
      WRITE( 6, 110   ) 'CO ', CO,  '  '
      WRITE( 6, 110   ) 'SO2', SO2, ' S'
 110  FORMAT( 'BRAVO anthropogenic ', a3, ': ', f9.4, ' Tg', a2 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_TG

!------------------------------------------------------------------------------

      SUBROUTINE READ_BRAVO_MASK
!
!******************************************************************************
!  Subroutine READ_BRAVO_MASK reads the Mexico mask from disk.  The Mexico
!  mask is the fraction of the grid box (I,J) which lies w/in the BRAVO
!  Mexican emissions region. (rjp, kfb, bmy, 6/22/06, 8/9/06)
!
!  NOTES:
!  (1 ) Now pass UNIT to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1,   DO_REGRID_G2G_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1)
      REAL*8                  :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                  :: XTAU
      CHARACTER(LEN=255)      :: FILENAME 

      !=================================================================
      ! READ_BRAVO_MASK begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'BRAVO_200607/BRAVO.MexicoMask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_MEXICO_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1999 )

      ! Mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 XTAU,      I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1) 

      ! Regrid from GENERIC 1x1 GRID to GEOS 1x1 GRID
      CALL DO_REGRID_G2G_1x1( 'unitless', GEN_1x1, GEOS_1x1(:,:,1) )
      
      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, BRAVO_MASK )

      ! Return to calling program
      END SUBROUTINE READ_BRAVO_MASK

!------------------------------------------------------------------------------

      SUBROUTINE INIT_BRAVO
!
!******************************************************************************
!  Subroutine INIT_BRAVO allocates and zeroes BRAVO module arrays, and also
!  creates the mask which defines the Mexico region (rjp, kfb, bmy, 6/26/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_XMID, GET_YMID
      USE LOGICAL_MOD, ONLY : LBRAVO

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER              :: AS

      !=================================================================
      ! INIT_BRAVO begins here!
      !=================================================================

      ! Return if LBRAVO is false
      IF ( .not. LBRAVO ) RETURN
      
      !--------------------------
      ! Allocate and zero arrays
      !--------------------------

      ALLOCATE( BRAVO_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_NOx' )
      BRAVO_NOx = 0d0

      ALLOCATE( BRAVO_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_CO' )
      BRAVO_CO = 0d0

      ALLOCATE( BRAVO_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_SO2' )
      BRAVO_SO2 = 0d0

      !--------------------------
      ! Read Mexico mask
      !--------------------------
     
      ALLOCATE( BRAVO_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_MASK' )
      BRAVO_MASK = 0d0
      
      ! Read the mask
      CALL READ_BRAVO_MASK

      ! Return to calling program
      END SUBROUTINE INIT_BRAVO

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BRAVO
!
!******************************************************************************
!  Subroutine CLEANUP_BRAVO deallocates all BRAVO module arrays.
!  (rjp, kfb, bmy, 6/26/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_BRAVO begins here!
      !=================================================================
      IF ( ALLOCATED( BRAVO_NOx  ) ) DEALLOCATE( BRAVO_NOx  )
      IF ( ALLOCATED( BRAVO_CO   ) ) DEALLOCATE( BRAVO_CO   )
      IF ( ALLOCATED( BRAVO_SO2  ) ) DEALLOCATE( BRAVO_SO2  )
      IF ( ALLOCATED( BRAVO_MASK ) ) DEALLOCATE( BRAVO_MASK )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BRAVO

!------------------------------------------------------------------------------

      ! End of module
      END MODULE BRAVO_MOD
