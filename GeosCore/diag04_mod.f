! $Id: diag04_mod.f,v 1.1 2009/09/16 14:06:36 bmy Exp $
      MODULE DIAG04_MOD
!
!******************************************************************************
!  Module DIAG04_MOD contains arrays and routines for archiving the ND04
!  diagnostic -- CO2 emissions and fluxes (bmy, 7/26/05, 9/5/06) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD04         (REAL*4) : Array for Hg emissions & ocean masses 
!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG04           : Sets all module arrays to zero
!  (2 ) WRITE_DIAG04          : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG04           : Allocates all module arrays
!  (4 ) CLEANUP_DIAG04        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag04_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f           : Module w/ NaN and other error check routines
!  (3 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f            : Module w/ horizontal grid information
!  (5 ) time_mod.f            : Module w/ routines to compute date & time
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag04_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND04
      INTEGER, PARAMETER   :: PD04 = 6

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD04(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG04
!
!******************************************************************************
!  Subroutine ZERO_DIAG04 zeroes the ND04 diagnostic array (bmy, 7/26/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! ZERO_DIAG04 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Zero array
      AD04(:,:,:) = 0e0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG04
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 7/26/05, 9/3/06)
!
!   # : Field     : Description                  : Units       : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) CO2-SRCE  : CO2 fossil fuel emissions    : molec/cm2/s : SCALE
!  (2 ) CO2-SRCE  : CO2 ocean emissions          : molec/cm2/s : SCALE
!  (3 ) CO2-SRCE  : CO2 balanced biosphere       : molec/cm2/s : SCALE
!  (4 ) CO2-SRCE  : CO2 biomass emissions        : molec/cm2/s : SCALE
!  (5 ) CO2-SRCE  : CO2 biofuel emissions        : molec/cm2/s : SCALE
!  (6 ) CO2-SRCE  : CO2 net terrestrial exchange : molec/cm2/s : SCALE
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,  ONLY : IU_BPCH
      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,  ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! TINDEX

      ! Local variables
      INTEGER            :: CENTER180, HALFPOLAR, IFIRST, JFIRST 
      INTEGER            :: LFIRST,    LMAX,      M,      N       
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*4             :: LONRES,    LATRES
      REAL*8             :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)  :: MODELNAME 
      CHARACTER(LEN=40)  :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG04 begins here!
      !=================================================================

      ! Exit if ND04 is turned off
      IF ( ND04 == 0 ) RETURN

      ! Initialize
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = GET_HALFPOLAR()
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      SCALE     = DBLE( GET_CT_EMIS() ) + 1d-32

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND04 diagnostic tracers
      DO M = 1, TMAX(4)

         ! Get quantities
         N            = TINDEX(4,M)
         CATEGORY     = 'CO2-SRCE'
         !UNIT         = 'molec/cm2/s'
         UNIT         = ''                     ! Let GAMAP pick the unit
         ARRAY(:,:,1) = AD04(:,:,N) / SCALE

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG04
!
!******************************************************************************
!  Subroutine INIT_DIAG04 allocates all module arrays (bmy, 7/26/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
   
#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG04 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND04 == 0 ) RETURN

      ! 2-D array ("CO2-SRCE")
      ALLOCATE( AD04( IIPAR, JJPAR, PD04 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD04' )

      ! Zero arrays
      CALL ZERO_DIAG04

      ! Return to calling program
      END SUBROUTINE INIT_DIAG04

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG04
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG04 deallocates all module arrays (bmy, 7/26/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG04 begins here!
      !=================================================================
      IF ( ALLOCATED( AD04 ) ) DEALLOCATE( AD04 ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG04

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG04_MOD
