! $Id: diag_oh_mod.f,v 1.1 2004/09/21 18:04:11 bmy Exp $
      MODULE DIAG_OH_MOD
!
!******************************************************************************
!  Module DIAG_OH_MOD contains routines and variables to archive OH mass
!  and air mass concentrations.  These are then used to print out the mass-
!  weighted mean OH concentration in 1e5 molec/cm3.  This is a metric of
!  how certain chemisry simulations are performing. (bmy, 7/20/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AIR_MASS (REAL*8) : Array used to sum mean air mass [molec/cm3]
!  (2 ) OH_MASS  (REAL*8) : Array used to sum mean OH mass  [molec/cm3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_DIAG_OH        : Driver routine for mean OH diagnostic (fullchem)
!  (2 ) DO_DIAG_OH_CH4    : Driver routine for mean OH diagnostic (CH4 sim)
!  (3 ) PRINT_DIAG_OH     : Prints the mean OH concentration [1e5 molec/cm3]
!  (4 ) INIT_DIAG_OH      : Allocates and zeroes module arrays
!  (5 ) CLEANUP_DIAG_OH   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "diag_oh_mod.f":
!  ============================================================================
!  (1 ) comode_mod.f      : Module containing SMVGEAR allocatable arrays
!  (2 ) error_mod.f       : Module containing I/O error and NaN check routines 
!  (3 ) logical_mod.f     : Module containing GEOS-CHEM logical switches
!  (4 ) tracer_mod.f      : Module containing GEOS-CHEM tracer array STT etc.  
!  (5 ) tracerid_mod.f    : Module containing pointers to tracers & emissions 
!
!  NOTES:
!******************************************************************************
! 
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag_oh_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_DIAG_OH
      PUBLIC :: DO_DIAG_OH
      PUBLIC :: DO_DIAG_OH_CH4
      PUBLIC :: INIT_DIAG_OH
      PUBLIC :: PRINT_DIAG_OH
      
      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================
      
      ! Scalars
      LOGICAL             :: DO_SAVE_OH

      ! Arrays 
      REAL*8, ALLOCATABLE :: OH_MASS(:,:,:)
      REAL*8, ALLOCATABLE :: AIR_MASS(:,:,:)
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_DIAG_OH
!
!******************************************************************************
!  Subroutine DO_DIAG_OH sums the OH and air mass (from SMVGEAR arrays) for
!  the mean OH concentration diagnostic. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! Reference to F90 modules
      USE COMODE_MOD,   ONLY : AIRDENS, CSPEC, JLOP, T3, VOLUME
      USE TRACERID_MOD, ONLY : IDOH

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NPVERT, NLAT, NLONG

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I,     J,       L,       JLOOP
      REAL*8        :: XLOSS, XOHMASS, XAIRMASS

      !=================================================================
      ! DO_DIAG_OH begins here!
      !=================================================================

      ! Safety valve -- avoid seg faults
      IF ( .not. DO_SAVE_OH ) RETURN

      ! Loop over boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, XAIRMASS, XOHMASS, XLOSS )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, NPVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! 1-D grid box index
         JLOOP = JLOP(I,J,L)

         ! Cycle if this isn't a valid SMVGEAR gridbox
         IF ( JLOOP > 0 ) THEN
         
            ! Sum air mass term into AIR_MASS array
            XAIRMASS        = AIRDENS(JLOOP)    * VOLUME(JLOOP)
            AIR_MASS(I,J,L) = AIR_MASS(I,J,L)   + XAIRMASS

            ! Sum OH mass term into OH_MASS array
            XOHMASS         = CSPEC(JLOOP,IDOH) * XAIRMASS
            OH_MASS(I,J,L)  = OH_MASS(I,J,L)    + XOHMASS

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
         
      ! Return to calling program
      END SUBROUTINE DO_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE DO_DIAG_OH_CH4( I, J, L, XOHMASS, XAIRMASS )
! 
!******************************************************************************
!  Subroutine DO_DIAG_OH_CH4 passes the OH loss, OH mass, and air mass terms
!  from "global_ch4_mod.f" to "diag_oh_mod.f" (bmy, 7/20/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, & altitude indices
!  (4  ) XOHMASS  (REAL*8 ) : OH mass term from "global_ch4_mod.f"
!  (5  ) XAIRMASS (REAL*8 ) : air mass term from "global_ch4_mod.f"
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: XOHMASS, XAIRMASS

      !=================================================================
      ! DO_DIAG_OH_CH4 begins here!
      !=================================================================

      ! Sum air mass & OH mass into arrays
      AIR_MASS(I,J,L) = AIR_MASS(I,J,L) + XAIRMASS
      OH_MASS(I,J,L)  = OH_MASS(I,J,L)  + XOHMASS
      
      ! Return to calling program
      END SUBROUTINE DO_DIAG_OH_CH4

!------------------------------------------------------------------------------

      SUBROUTINE PRINT_DIAG_OH
!
!******************************************************************************
!  Subroutine PRINT_DIAG_OH prints the mass-weighted OH concentration at
!  the end of a simulation. (bmy, 10/21/03, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables 
      REAL*8 :: SUM_OHMASS, SUM_MASS, OHCONC
     
      !=================================================================
      ! PRINT_DIAG_OH begins here!
      !=================================================================

      ! Return if this diagnostic is turned off
      IF ( .not. DO_SAVE_OH ) RETURN

      ! Total Mass-weighted OH [molec OH/cm3] * [molec air]
      SUM_OHMASS = SUM( OH_MASS )

      ! Atmospheric air mass [molec air]
      SUM_MASS   = SUM( AIR_MASS ) 
         
      ! Avoid divide-by-zero errors 
      IF ( SUM_MASS > 0d0 ) THEN 
            
         ! Divide OH by [molec air] and report as [1e5 molec/cm3]
         OHCONC = ( SUM_OHMASS / SUM_MASS ) / 1d5
         
         ! Write value to log file
         WRITE( 6, '(/,a)' ) REPEAT( '=', 79 ) 
         WRITE( 6, *       ) 'ND23: Mass-Weighted OH Concentration'
         WRITE( 6, *       ) 'Mean OH = ', OHCONC, ' [1e5 molec/cm3]' 
         WRITE( 6, '(  a)' ) REPEAT( '=', 79 ) 

      ELSE

         ! Write error msg if SUM_MASS is zero
         WRITE( 6, '(/,a)' ) REPEAT( '=', 79 ) 
         WRITE( 6, '(  a)' ) 'Could not print mass-weighted OH!'
         WRITE( 6, '(  a)' ) 'Atmospheric air mass is zero!'
         WRITE( 6, '(  a)' ) REPEAT( '=', 79 ) 
            
      ENDIF

      ! Return to MAIN program
      END SUBROUTINE PRINT_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG_OH
!
!******************************************************************************
!  Subroutine INIT_DIAG_OH initializes all module arrays. (bmy, 7/20/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LCHEM
      USE TRACER_MOD,  ONLY : ITS_A_FULLCHEM_SIM,  
     &                        ITS_A_COPARAM_SIM, 
     &                        ITS_A_CH4_SIM

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND23

      ! Local variables
      INTEGER :: AS, LMAX

      !=================================================================
      ! INIT_DIAG_OH begins here!
      !=================================================================

      ! Initialize
      DO_SAVE_OH = .FALSE.

      ! Return if we are not doing chemistry
      IF ( .not. LCHEM ) RETURN

      ! Set vertical levels and decide whether to print CH3CCl3
      ! lifetime or just mean mass-weighted OH concentration
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Fullchem: tropopshere only
         LMAX       = LLTROP
         DO_SAVE_OH = .TRUE.

      ELSE IF ( ITS_A_COPARAM_SIM() ) THEN

         ! CO w/ OH param: all levels
         LMAX       = LLPAR
         DO_SAVE_OH = .TRUE.

      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         ! CH4: all levels
         LMAX       = LLPAR
         DO_SAVE_OH = .TRUE.

      ENDIF

      ! Echo info
      WRITE( 6, 100 ) DO_SAVE_OH
 100  FORMAT( /, 'Turn on Mean OH diagnostic (ND23)? :', L5 )
  
      ! Return if we aren't saving mean OH
      IF ( .not. DO_SAVE_OH ) RETURN

      !=================================================================
      ! Allocate arrays
      !=================================================================

      ! Air mass array
      ALLOCATE( AIR_MASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH_LOSS' )
      AIR_MASS = 0d0

      ! OH mass array
      ALLOCATE( OH_MASS( IIPAR, JJPAR, LMAX ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH_LOSS' )
      OH_MASS = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_DIAG_OH

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG_OH
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG_OH deallocates all module arrays. (bmy, 7/20/04)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG_OH begins here!
      !=================================================================
      IF ( ALLOCATED( OH_MASS  ) ) DEALLOCATE( OH_MASS  )
      IF ( ALLOCATED( AIR_MASS ) ) DEALLOCATE( AIR_MASS )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG_OH

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG_OH_MOD
