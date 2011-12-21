! $Id: Kr85_mod.f,v 1.1 2009/11/20 21:42:59 bmy Exp $
      MODULE Kr85_MOD
!
!******************************************************************************
!  Module Kr85_MOD contains routines and variables for the Kr85 radionuclide
!  simulation. (jsw, bmy, 8/21/03, 11/6/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) N_SOURCES (INTEGER) : Maximum number of Kr85 point sources
!  (2 ) N_YEARS   (INTEGER) : Maximum number of years for Kr85 emissions
!  (3 ) SMALLNUM  (REAL*8 ) : A small number, used to prevent underflow
!
!  Module Routines:
!  ============================================================================
!  (1 ) GET_SOURCE_IJ       : Returns (I,J) location of each Kr85 point source
!  (2 ) GET_EMITTED_Kr85    : Returns Kr85 emission from a point src in [kg]
!  (3 ) EMISSKr85           : Adds Kr85 emissions into the tracer array
!  (4 ) CHEMKr85            : Performs radioactive (1st-order) loss for Kr85 
!
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (2 ) error_mod.f    : Module containing I/O error and NaN check routines
!  (3 ) time_mod.f     : Module containing routines for computing time & date
!  (4 ) tracer_mod.f   : Module containing GEOS-CHEM tracer array STT etc.
!
!  References:
!  ============================================================================
!  (1 ) Jacob, D.J., M.J. Prather, S.C. Wofsy, M.B. McElroy, "Atmospheric
!        distribution of 85Kr simulated with a general circulation model",
!        JGR, 92(D6), pp. 6614-6626, June 20, 1987.
!
!  NOTES:
!  (1 ) Now references "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Modifications for GEOS-5 nested grids (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!     
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "Kr85_mod.f"
      !=================================================================

      ! PRIVATE module routines
      PRIVATE :: GET_SOURCE_IJ
      PRIVATE :: GET_EMITTED_KR85
      
      ! PRIVATE module variables
      PRIVATE :: N_SOURCES, N_YEARS, SMALLNUM

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER, PARAMETER :: N_SOURCES = 8
      INTEGER, PARAMETER :: N_YEARS   = 6
      REAL*8,  PARAMETER :: SMALLNUM  = 1d-20

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================      
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_SOURCE_IJ( N_SOURCE, I, J )
!
!******************************************************************************
!  Subroutine GET_SOURCE_IJ returns the (I,J) grid box indices for each
!  Kr85 point source.  For now we have hardwired this, since there are only
!  a few stations.  Worry about making this more general at some future time.
!  (bmy, 8/21/03, 11/6/08) 
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) N_SOURCE (INTEGER) : Number of Kr85 point source (1-8)
!
!  Arguments as Output
!  ============================================================================
!  (1-2) I, J     (INTEGER) : Lon & lat indices for the N_SOURCEth Kr85 source
!
!  NOTES:
!  (1 ) Updated for 0.5 x 0.666 nested grids (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

#     include "define.h"   ! Switches

      ! Arguments
      INTEGER, INTENT(IN)  :: N_SOURCE
      INTEGER, INTENT(OUT) :: I, J
     
      ! Local variables
      INTEGER              :: GRID

      ! Station (I,J) arrays             4x5 2x25  1x1
      INTEGER              :: I1(3) = (/  21,  41, 101 /)
      INTEGER              :: J1(3) = (/  32,  63, 105 /)
      !---
      INTEGER              :: I2(3) = (/  14,  27,  66 /)
      INTEGER              :: J2(3) = (/  35,  69, 137 /)
      !---
      INTEGER              :: I3(3) = (/  36,  71, 176 /)
      INTEGER              :: J3(3) = (/  37,  73, 145 /)
      !---
      INTEGER              :: I4(3) = (/  38,  75, 186 /)
      INTEGER              :: J4(3) = (/  35,  69, 137 /)
      !---
      INTEGER              :: I5(3) = (/  37,  73, 181 /)
      INTEGER              :: J5(3) = (/  36,  71, 141 /)
      !---
      INTEGER              :: I6(3) = (/  39,  77, 191 /)
      INTEGER              :: J6(3) = (/  36,  71, 141 /)
      !---
      INTEGER              :: I7(3) = (/  65, 129, 321 /)
      INTEGER              :: J7(3) = (/  32,  63, 125 /)
      !---
      INTEGER              :: I8(3) = (/  49,  97, 241 /)
      INTEGER              :: J8(3) = (/  37,  73, 145 /)

      !=================================================================
      ! GET_SOURCE_IJ begins here!
      !=================================================================

      ! Select flag for grid type
#if   defined( GRID4x5  )
      GRID = 1
#elif defined( GRID2x25 )
      GRID = 2 
#elif defined( GRID1x1  )
      GRID = 3
#elif defined( GRID05x0666 )
      GRID = 3   !(dan )
#endif

      ! Select proper (I,J) for each station
      SELECT CASE( N_SOURCE )
         CASE( 1 )
            I = I1(GRID)
            J = J1(GRID)
         CASE( 2 )
            I = I2(GRID)
            J = J2(GRID)
         CASE( 3 )
            I = I3(GRID)
            J = J3(GRID)
         CASE( 4 )
            I = I4(GRID)
            J = J4(GRID)
         CASE( 5 )
            I = I5(GRID)
            J = J5(GRID)
         CASE( 6 )
            I = I6(GRID)
            J = J6(GRID)
         CASE( 7 )
            I = I7(GRID)
            J = J7(GRID)
         CASE( 8 )
            I = I8(GRID)
            J = J8(GRID)
         CASE DEFAULT
            CALL ERROR_STOP( 'N_SOURCE must be between 1-8!',
     &                       'GET_SOURCE_IJ (Kr85_mod.f)' )
      END SELECT
     
      ! Return to calling program
      END SUBROUTINE GET_SOURCE_IJ

!------------------------------------------------------------------------------

      FUNCTION GET_EMITTED_Kr85( N_SOURCE, YEARCOUNT ) RESULT( Kr85 )
!
!******************************************************************************
!  Subroutine GET_EMITTED_Kr85 returns the amount of Kr85 emitted from a
!  given point source 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N_SOURCE  (INTEGER) : Kr85 point source index (1-N_SOURCES)
!  (2 ) YEARCOUNT (INTEGER) : Year of Kr85 emissions to use (1-N_YEARS)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE TIME_MOD,  ONLY : GET_TS_EMIS 

      ! Arguments
      INTEGER, INTENT(IN)  :: N_SOURCE, YEARCOUNT

      ! Local variables
      REAL*8               :: Kr85, DTSRCE
      REAL*8, PARAMETER    :: SEC_PER_YR = 86400d0 * 365.25d0

      ! Kr85 point src emissions
      ! Units: MCi/year     Year1  Year2  Year3  Year4  Year5  Year6
      REAL*8 :: S1(6) = (/ 0.48d0,0.46d0,0.47d0,0.54d0,0.56d0,0.69d0 /)
      REAL*8 :: S2(6) = (/ 0.10d0,0.00d0,0.09d0,0.06d0,0.01d0,0.00d0 /)
      REAL*8 :: S3(6) = (/ 0.70d0,0.94d0,0.84d0,1.40d0,1.19d0,1.13d0 /)
      REAL*8 :: S4(6) = (/ 0.31d0,0.28d0,0.54d0,0.31d0,0.31d0,0.31d0 /)
      REAL*8 :: S5(6) = (/ 0.79d0,0.64d0,0.83d0,0.91d0,1.27d0,1.95d0 /)
      REAL*8 :: S6(6) = (/ 0.00d0,0.05d0,0.03d0,0.00d0,0.00d0,0.08d0 /)
      REAL*8 :: S7(6) = (/ 0.06d0,0.00d0,0.28d0,0.11d0,0.19d0,0.09d0 /)
      REAL*8 :: S8(6) = (/ 3.56d0,3.77d0,3.19d0,3.07d0,3.00d0,2.40d0 /)

      !=================================================================
      ! GET_EMITTED_Kr85 begins here!
      !=================================================================

      ! Error check year
      IF ( YEARCOUNT < 1 .or. YEARCOUNT > 6 ) THEN
         CALL ERROR_STOP( 'YEARCOUNT must be between 1-6!',
     &                    'GET_EMITTED_KR85 (Kr85_mod.f)' )
      ENDIF

      ! Return Kr85 for the given point source & year
      SELECT CASE( N_SOURCE )
         CASE( 1 )
            Kr85 = S1(YEARCOUNT)
         CASE( 2 )
            Kr85 = S2(YEARCOUNT)
         CASE( 3 )
            Kr85 = S3(YEARCOUNT)
         CASE( 4 )
            Kr85 = S4(YEARCOUNT)
         CASE( 5 )
            Kr85 = S5(YEARCOUNT)
         CASE( 6 )
            Kr85 = S6(YEARCOUNT)
         CASE( 7 )
            Kr85 = S7(YEARCOUNT)
         CASE( 8 )
            Kr85 = S8(YEARCOUNT)
         CASE DEFAULT
            CALL ERROR_STOP( 'N_SOURCE must be between 1-8!',
     &                       'GET_SOURCE_IJ (Kr85_mod.f)' )
      END SELECT

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Convert from [MCi/yr] to [kg/emission timestep] 
      ! 1 kg of Kr85 is equivalent to 2.55 MCi (cf Jacob et al 1987)
      Kr85   = Kr85 * 2.55d0 * ( DTSRCE / SEC_PER_YR )

      ! Return to calling program
      END FUNCTION GET_EMITTED_Kr85

!------------------------------------------------------------------------------

      SUBROUTINE EMISSKr85
!
!******************************************************************************
!  Subroutine EMISSKr85 places Kr85 emissions from point sources (e.g. nuclear
!  reprocessing plants) into the tracer array. (jsw, bmy, 8/21/03, 7/20/04)
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      !--------------------------------------------------------------
      ! Prior to 12/7/04:
      ! Need to reassign the diagnostic number
      !USE DIAG_MOD,   ONLY : AD03
      !--------------------------------------------------------------      
      USE TIME_MOD,   ONLY : GET_TS_EMIS, GET_YEAR
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! Diagnostics
#     include "CMN_O3"        ! FSCALYR

      ! Local Variables
      LOGICAL, SAVE          :: FIRSTEMISS = .TRUE. 
      INTEGER, SAVE          :: YEARCOUNT, LASTYEAR
      INTEGER                :: I, J, N
      REAL*8                 :: Kr85_KG, DTSRCE

      !=================================================================
      ! EMISSKr85 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      ! First-time initialization
      IF ( FIRSTEMISS ) THEN
         YEARCOUNT  = FSCALYR
         LASTYEAR   = GET_YEAR()
         FIRSTEMISS = .FALSE.

         WRITE( 6, 100 ) YEARCOUNT
 100     FORMAT( '     - EMISSKr85: Using Kr85 emissions from year ',i3)
      ENDIF

      ! If it's a new year, increment YEARCOUNT
      IF ( GET_YEAR() /= LASTYEAR ) THEN
         YEARCOUNT = YEARCOUNT + 1
         LASTYEAR  = GET_YEAR()
         
         WRITE( 6, 100 ) YEARCOUNT
      ENDIF

      !=================================================================
      ! Add Kr85 emissions [kg] to the STT tracer array
      ! NOTE: Assumes a global (not a window!) simulation
      !=================================================================
      DO N = 1, N_SOURCES
       
         ! Get (I,J) for each Kr85 point source
         CALL GET_SOURCE_IJ( N, I, J )

         ! Get emitted Kr85 from each point source [kg]
         Kr85_KG = GET_EMITTED_Kr85( N, YEARCOUNT )

         ! Add Kr85 into STT array
         STT(I,J,1,1) = STT(I,J,1,1) + Kr85_KG

         !--------------------------------------------------------------
         ! Prior to 12/7/04:
         ! Need to reassign the diagnostic number (bmy, 12/7/04)
         !! Archive emitted Kr85 for ND04 diagnostic [kg]
         !IF ( ND03 > 0 ) THEN
         !   AD03(I,J,1,1) = AD03(I,J,1,1) + Kr85_KG
         !ENDIF
         !--------------------------------------------------------------
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSKr85

!------------------------------------------------------------------------------

      SUBROUTINE CHEMKr85
!
!******************************************************************************
!  Subroutine CHEMKr85 applies first-order loss to the Kr85 tracer.
!  (jsw, bmy, 8/21/03, 7/20/04)
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      !------------------------------------------
      ! Prior to 12/7/04:
      ! Need to reassign diagnostic number
      !USE DIAG_MOD,   ONLY : AD03
      !------------------------------------------
      USE TIME_MOD,   ONLY : GET_TS_CHEM
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! ND03

      ! Local variables
      INTEGER :: I,      J,     L
      REAL*8  :: DTCHEM, KRATE, LOSS_FACTOR, Kr_LOST

      !=================================================================
      ! CHEMKr85 begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM      = GET_TS_CHEM() * 60d0 

      ! The decay for 85Kr is calculated by: dC/dt = -kC
      ! where k = 1/15.52yr = 2.042E-9 s^-1
      KRATE       = 2.042d-9 

      ! Multiplication factor to compute tracer lost 
      LOSS_FACTOR = 1d0 - EXP( -2.042d-9 * DTCHEM )

      !=================================================================
      ! Apply radioactive decay to Kr85 tracer
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, Kr_LOST )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Amount of Kr85 lost to radioactive decay [kg]
         Kr_LOST = STT(I,J,L,1) * LOSS_FACTOR

         ! Prevent underflow
         IF ( Kr_LOST < SMALLNUM ) Kr_LOST = 0d0
             
         ! Subtract Kr85 lost from the tracer array 
         STT(I,J,L,1) = STT(I,J,L,1) - Kr_LOST
         
         !-------------------------------------------------------------
         ! Prior to 12/7/04:
         ! Need to reassign the diagnostic number
         !! Archive Kr85 lost by decay [kg] in ND04 diagnostic
         !IF ( ND03 > 0 ) THEN
         !   AD03(I,J,L,2) = AD03(I,J,L,2) + Kr_LOST
         !ENDIF
         !-------------------------------------------------------------
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEMKr85

!------------------------------------------------------------------------------

      END MODULE Kr85_MOD
