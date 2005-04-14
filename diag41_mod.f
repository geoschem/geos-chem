! $Id: diag41_mod.f,v 1.1 2005/04/14 14:34:30 bmy Exp $
      MODULE DIAG41_MOD
!
!******************************************************************************
!  Module DIAG41_MOD contains arrays and routines for archiving the ND41
!  diagnostic -- Afternoon PBL heights. (bmy, 2/17/05) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) AD41     (REAL*4 ) : Array for afternoon PBL height
!  (2 ) GOOD_CT  (INTEGER) : Counter of grid boxes where it's afternoon
!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG41        : Sets all module arrays to zero
!  (2 ) WRITE_DIAG41       : Writes data in module arrays to bpch file
!  (3 ) DIAG41             : Archives afternoon PBL heights
!  (4 ) INIT_DIAG41        : Allocates all module arrays
!  (4 ) CLEANUP_DIAG41     : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag41_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f        : Module w/ NaN and other error check routines
!  (3 ) file_mod.f         : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f         : Module w/ horizontal grid information
!  (5 ) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (6 ) time_mod.f         : Module w/ routines to compute date & time!      
!
!  NOTES:
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag41_mod.f"
      !=================================================================

      ! Make everything PUBLIC ...
      PUBLIC

      ! ... except these routines
      PRIVATE :: AD41
      PRIVATE :: GOOD_CT

      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND41
      INTEGER, PARAMETER   :: PD41 = 2

      ! Arrays
      INTEGER, ALLOCATABLE :: GOOD_CT(:)
      REAL*4,  ALLOCATABLE :: AD41(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG41
!
!******************************************************************************
!  Subroutine ZERO_DIAG41 zeroes the ND41 diagnostic arrays (bmy, 2/17/05)
!
!  NOTES:
!******************************************************************************
!
#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: I, J, N

      !=================================================================
      ! ZERO_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Zero GOOD_CT
      DO I = 1, IIPAR
         GOOD_CT(I) = 0
      ENDDO

      ! Zero AD41
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )
      DO N = 1, PD41
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD41(I,J,N) = 0e0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ZERO_DIAG41

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG41
!
!******************************************************************************
!  Subroutine WRITE_DIAG41 writes the ND41 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 2/17/05)
!
!  ND41: Afternoon PBL depth (between 1200 and 1600 Local Time)
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) PBLDEPTH : Afternoon PBL heights       : m         : GOOD_CT
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE FILE_MOD, ONLY : IU_BPCH
      USE GRID_MOD, ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD, ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! TINDEX

      ! Local variables
      INTEGER           :: I,         J,           M,      N
      INTEGER           :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER           :: JFIRST,    LFIRST,      LMAX
      REAL*4            :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4            :: LONRES,    LATRES,      EPS
      REAL*8            :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20) :: MODELNAME 
      CHARACTER(LEN=40) :: CATEGORY,  RESERVED,    UNIT

      !=================================================================
      ! WRITE_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Initialize
      CATEGORY  = 'PBLDEPTH'
      CENTER180 = 1
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      HALFPOLAR = 1
      IFIRST    = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST    = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LATRES    = DJSIZE
      LFIRST    = 1
      LONRES    = DISIZE
      MODELNAME = GET_MODELNAME()
      RESERVED  = ''
      EPS       = TINY( 1d0 )
         
      !=================================================================
      ! Write data to the bpch file
      !=================================================================
      
      ! Loop over ND41 diagnostic tracers
      DO M = 1, TMAX(41)
         N = TINDEX(41,M)
         IF ( N > PD41 ) CYCLE

         ! Select proper unit string
         IF ( N == 1 ) UNIT = 'm' 
         IF ( N == 2 ) UNIT = 'level'
                     
         ! Divide by # of afternoon boxes at each longitude
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            SCALE        = DBLE( GOOD_CT(I) ) + EPS
            ARRAY(I,J,1) = AD41(I,J,N)        / SCALE
         ENDDO
         ENDDO

         ! Write to bpch file
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG41

!------------------------------------------------------------------------------

      SUBROUTINE DIAG41 
!
!******************************************************************************
!  Subroutine DIAG41 produces monthly mean boundary layer height in meters 
!  between 1200-1600 local time for the U.S. geographical domain. 
!  (amf, swu, bmy, 11/18/99, 11/6/03)
!
!  Input via "CMN" header file:
!  ===========================================================================
!  (1 ) XTRA2 : Height of PBL in boxes
!
!  NOTES:
!  (1 ) DIAG41 is written in Fixed-Format F90. 
!  (2 ) XTRA2 must be computed by turning TURBDAY on first.  Also,
!        XTRA2 is a global-size array, so use window offsets IREF, JREF
!        to index it correctly. (bmy, 11/18/99)
!  (3 ) Do a little rewriting so that the DO-loops get executed
!        in the correct order (J first, then I). (bmy, 11/18/99)
!  (4 ) AD41 is now declared allocatable in "diag_mod.f". (bmy, 12/6/99)
!  (5 ) AFTTOT is now declared allocatable in "diag_mod.f". (bmy, 3/17/00)
!  (6 ) Remove NYMD from the argument list -- it wasn't used (bmy, 6/22/00) 
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        (bmy, 9/25/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) Now reference BXHEIGHT from "dao_mod.f".  Also removed obsolete
!        code. (bmy, 9/18/02)
!  (10) Now use function GET_LOCALTIME from "dao_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in DO-loop for calculating local time (bmy, 7/9/03)
!  (12) For GEOS-4, PBL depth is already in meters, so we only have to
!        multiply that by the GOOD array.  Also now references PBL array
!        from "dao_mod.f".  Bug fix: now use barometric law to compute PBL 
!        height in meters for GEOS-1, GEOS-STRAT, GEOS-3.  This eliminates an 
!        overprediction of the PBL height. (swu, bmy, 11/6/03)
!******************************************************************************
!
      ! References to F90 modules
      USE PBL_MIX_MOD, ONLY : GET_PBL_TOP_L, GET_PBL_TOP_m
      USE TIME_MOD,    ONLY : GET_LOCALTIME

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      INTEGER              :: I, J, N, GOOD(IIPAR)
      REAL*8               :: LT, PBLTOP

      !=================================================================
      ! DIAG41 begins here!
      !=================================================================

      !-----------------------------------
      ! Find boxes where it is afternoon
      !-----------------------------------
      DO I = 1, IIPAR

         ! Local time
         LT = GET_LOCALTIME( I )
   
         ! Find points between 12 and 16 GMT
         IF ( LT >= 12d0 .and. LT <= 16d0 ) THEN
            GOOD(I) = 1
         ELSE
            GOOD(I) = 0
         ENDIF

         ! Increment counter of afternoon boxes
         GOOD_CT(I)  = GOOD_CT(I) + GOOD(I)
      ENDDO

      !-----------------------------------
      ! Archive afternoon PBL heights
      !-----------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, PBLTOP )
      DO N = 1, PD41
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         IF ( N == 1 ) THEN

            ! Afternoon PBL top [m]
            PBLTOP = GET_PBL_TOP_m( I, J ) * GOOD(I)

         ELSE IF ( N == 2 ) THEN

            ! Afternoon PBL top [model layers]
            PBLTOP = GET_PBL_TOP_L( I, J ) * GOOD(I)

         ENDIF
           
         ! Store in AD41 array
         AD41(I,J,N) = AD41(I,J,N) + PBLTOP

      ENDDO    
      ENDDO    
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE DIAG41

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG41
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG41 allocates all module arrays (bmy, 2/17/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
   
#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER            :: AS
      
      !=================================================================
      ! INIT_DIAG41 begins here!
      !=================================================================

      ! Exit if ND41 is turned off
      IF ( ND41 == 0 ) RETURN

      ! Counter of afternoon pts
      ALLOCATE( GOOD_CT( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD_CT' )

      ! Diagnostic array
      ALLOCATE( AD41( IIPAR, JJPAR, PD41 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD41' )

      ! Zero arrays
      CALL ZERO_DIAG41

      ! Return to calling program
      END SUBROUTINE INIT_DIAG41

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG41
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG41 deallocates all module arrays (bmy, 2/17/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG41 begins here!
      !=================================================================
      IF ( ALLOCATED( AD41    ) ) DEALLOCATE( AD41    ) 
      IF ( ALLOCATED( GOOD_CT ) ) DEALLOCATE( GOOD_CT )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG41

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG41_MOD
