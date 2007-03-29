! $Id: diag56_mod.f,v 1.5 2007/03/29 20:31:14 bmy Exp $
      MODULE DIAG56_MOD
!
!******************************************************************************
!  Module DIAG56_MOD contains arrays and routines for archiving the ND56
!  diagnostic -- lightning flash rates. (bmy, 5/11/06, 3/7/07) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD56 (REAL*4)  : Diagnostic array for lightning flash rates
!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG56    : Sets all module arrays to zero
!  (2 ) WRITE_DIAG56   : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG56    : Allocates all module arrays
!  (4 ) CLEANUP_DIAG56 : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag03_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f    : Module w/ NaN and other error check routines
!  (3 ) file_mod.f     : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f     : Module w/ horizontal grid information
!  (5 ) time_mod.f     : Module w/ routines to compute date & time
!
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Now divide AD56 by the # of A-6 timesteps (ltm, bmy, 3/7/07)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag56_mod.f"
      !=================================================================

      ! Make everything PUBLIC
      PUBLIC 

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND56
      INTEGER, PARAMETER   :: PD56 = 3

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD56(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG56
!
!******************************************************************************
!  Subroutine ZERO_DIAG03 zeroes the ND03 diagnostic arrays. 
!  (bmy, 5/11/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! ZERO_DIAG56 begins here!
      !=================================================================

      ! Exit if ND56 is turned off
      IF ( ND56 == 0 ) RETURN

      ! Zero arrays
      AD56(:,:,:) = 0e0

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG56

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG56
!
!******************************************************************************
!  Subroutine WRITE_DIAG56 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 5/11/06, 3/7/06)
!
!   # : Field    : Description              : Units          : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) LFLASH-$ : Lightning flash rate     : flashes/min/km2 : SCALE_A6
!  (2 ) LFLASH-$ : Intra-cloud flash rate   : flashes/min/km2 : SCALE_A6
!  (3 ) LFLASH-$ : Cloud-ground flash rate  : flashes/min/km2 : SCALE_A6
!
!  NOTES:
!  (1 ) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (2 ) Now scale AD56 by the # of A-6 timesteps (ltm, bmy, 3/7/07)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      !--------------------------------------------------------------------
      ! Prior to 3/7/07:
      !USE TIME_MOD,     ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
      !--------------------------------------------------------------------
      USE TIME_MOD,     ONLY : GET_CT_A6,   GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! TINDEX

      ! Local variables
      INTEGER               :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER               :: JFIRST,    LFIRST,      M,      N         
      REAL*4                :: ARRAY(IIPAR,JJPAR,1)
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG56 begins here!
      !=================================================================

      ! Exit if ND03 is turned off
      IF ( ND56 == 0 ) RETURN

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
      !------------------------------------------------
      ! Prior to 3/7/07
      !SCALE     = DBLE( GET_CT_EMIS() ) + 1d-32
      !------------------------------------------------
      SCALE     = DBLE( GET_CT_A6() ) + 1d-32
        
      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND03 diagnostic tracers
      DO M = 1, TMAX(56)

         ! Define quantities
         N            = TINDEX(56,M)
         CATEGORY     = 'LFLASH-$'
         UNIT         = 'flashes/min/km2'
         !---------------------------------------------
         ! Prior to 3/7/07:
         !ARRAY(:,:,1) = AD56(:,:,N)
         !---------------------------------------------
         ARRAY(:,:,1) = AD56(:,:,N) / SCALE

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG56

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG56
!
!******************************************************************************
!  Subroutine INIT_DIAG56 allocates all module arrays (bmy, 5/11/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
   
#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG03 begins here!
      !=================================================================

      ! Exit if ND56 is turned off
      IF ( ND56 == 0 ) RETURN

      ! 2-D array ("LFLASH-$")
      ALLOCATE( AD56( IIPAR, JJPAR, PD56 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD56' )

      ! Zero arrays
      CALL ZERO_DIAG56

      ! Return to calling program
      END SUBROUTINE INIT_DIAG56

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG56
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG56 deallocates all module arrays (bmy, 5/11/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG56 begins here!
      !=================================================================
      IF ( ALLOCATED( AD56 ) ) DEALLOCATE( AD56 ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG56

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG56_MOD
