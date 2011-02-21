!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag53_mod.f
!
! !DESCRIPTION: !  Module DIAG53_MOD contains arrays and routines 
!  for archiving the ND53
!  diagnostic --POPS emissions, mass, and production. 
!  ( eck 9/20/10 based on DIAG03_MOD, bmy, 1/21/05, 9/5/06) 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG53_MOD
! 
! !USES:
!
      IMPLICIT NONE
! Make everything Public ...
      PUBLIC
!
! !PUBLIC TYPES:
!

! !PUBLIC MEMBER FUNCTIONS:
!
!
! !PUBLIC DATA MEMBERS:
!
! !REVISION HISTORY:
!  20 September 2010 N.E. Selin - Initial Version
!
! !REMARKS:
!EOP


!  
!  Module Variables:
!  ============================================================================
!  (1 ) AD53         (REAL*4) : Array for POPS emissions
!  (2 ) AD53_PG_PP   (REAL*4) : Array for partitioning of gas phase POP to particles 
!  (3 ) AD53_POPG_OH (REAL*4) : Array for POPG oxidized by OH

!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG53           : Sets all module arrays to zero
!  (2 ) WRITE_DIAG53          : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG53           : Allocates all module arrays
!  (4 ) CLEANUP_DIAG53        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag03_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary pch file I/O
!  (2 ) error_mod.f           : Module w/ NaN and other error check routines
!  (3 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f            : Module w/ horizontal grid information
!  (5 ) time_mod.f            : Module w/ routines to compute date & time
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) POPG                  : Gas phase POP
!  (2 ) POPP                  : PARTICULATE PHASE POP


      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Scalars
      INTEGER              :: ND53, LD53
      INTEGER, PARAMETER   :: PD53 = 4  ! Number of diagnostics for AD03

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD53(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_PG_PP(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPG_OH(:,:,:)
!      REAL*4,  ALLOCATABLE :: AD53_POPG_NO3(:,:,:)
!      REAL*4,  ALLOCATABLE :: AD53_POPG_OX(:,:,:)
 
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
     
!------------------------------------------------------------------------------

      SUBROUTINE ZERO_DIAG53
!
!******************************************************************************
!  Subroutine ZERO_DIAG53 zeroes the ND53 diagnostic arrays. 
!  (bmy, 1/21/05, 4/6/06) (now diag 53, eck 9/20/10)
!
!  NOTES:
!  (1 ) Now references N_Hg_CATS from "tracerid_mod.f".  Now zero AD03_Hg2_SS
!        array. (bmy, 4/6/06)
!  (2 ) Now use broadcast assignment and double precision 0D0 to zero arrays,
!        rather than nested DO loops and single precision 0E0. (cdh, 8/14/08)
!******************************************************************************
!
      ! References to F90 modules
  

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: I, J, L, N

      !=================================================================
      ! ZERO_DIAG53 begins here!
      !=================================================================

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) RETURN

      ! Zero arrays
      AD53         = 0D0
      AD53_PG_PP   = 0D0
      AD53_POPG_OH = 0D0
  
      ! Return to calling program
      END SUBROUTINE ZERO_DIAG53

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG53
!
!******************************************************************************
!  Subroutine WRITE_DIAG03 writes the ND03 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 1/21/05, 2/24/06)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) PG-SRCE  : Total POP emission                  : kg       : 1
!  (2 ) PG-SRCE  : Total gas phase POP reacted with OH : kg       : 1
!  (21) PL-PP-$  : Amount of POPG CONVERTED TO PARTICLE: kg       : 1 
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
      USE TIME_MOD,     ONLY : GET_CT_CHEM ! CDH for sea salt loss rate
!      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! TINDEX

      ! Local variables
      INTEGER               :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER               :: JFIRST,    LFIRST,      LMAX
      INTEGER               :: M,         N,           NN
      REAL*4                :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: DIAGb,     DIAGe,       SCALE
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY,  RESERVED,    UNIT
      REAL*8                :: NCHEMSTEP !CDH for sea salt loss rate

      !=================================================================
      ! WRITE_DIAG53 begins here!
      !=================================================================

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) RETURN

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
      NCHEMSTEP = DBLE( GET_CT_CHEM() ) + TINY( 1d0 ) !CDH for sea salt loss rat         
      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND53 diagnostic tracers
      DO M = 1, TMAX(53)

         ! Get ND53 tracer #
         N = TINDEX(53,M)

         ! Pick the proper array & dimensions
         IF ( N == 1 .or. N == 2 .or. N == 3 .or. N == 4  ) THEN
               
            !--------------------------------
            ! #1 POP emissions (total, OC-sorbed, BC-sorbed, and gas phase)
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)     = AD53(:,:,N)

         ELSE IF ( N == 5  ) THEN

            !--------------------------------
            ! #2 Amount converted to OC particle
            !--------------------------------
            CATEGORY          = 'PL-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_PP(:,:,1:LMAX)

         ELSE IF ( N == 6  ) THEN

            !--------------------------------
            ! #2 Amount converted to BC particle
            !--------------------------------
            CATEGORY          = 'PL-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_PP(:,:,1:LMAX)
        
         ELSE IF ( N == 7  ) THEN

            !--------------------------------
            ! #3 Production of oxidized POPG from rxn with OH (clf, 1/27/11)
            !--------------------------------
            CATEGORY          = 'PL-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_POPG_OH(:,:,1:LMAX)               
  

            !--------------------------------
            ! Otherwise skip to next N
            !--------------------------------
            CYCLE

         ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LMAX,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )
      ENDDO

      ! Return to calling program
      END SUBROUTINE WRITE_DIAG53

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG53
!
!******************************************************************************
!  Subroutine INIT_DIAG53 allocates all module arrays (bmy, 1/21/05, 4/6/06)
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR

   
#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS
      
      !=================================================================
      ! INIT_DIAG53 begins here!
      !=================================================================

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) THEN
         LD53 = 0
         RETURN
      ENDIF

      ! Get number of levels for 3-D arrays
      LD53 = MIN( ND53, LLPAR )

      ! 2-D array ("PP-SRCE")
      ALLOCATE( AD53( IIPAR, JJPAR, PD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53' )

      ! 3-D arrays ("PL-PG-$")
      ALLOCATE( AD53_PG_PP( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PP_PG' )

      ALLOCATE( AD53_POPG_OH( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPG_OH' )

    
      ! Zero arrays
      CALL ZERO_DIAG53

      ! Return to calling program
      END SUBROUTINE INIT_DIAG53

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG53
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG53 deallocates all module arrays 
!  (bmy, 1/21/05, 4/6/06)
!
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG03 begins here!
      !=================================================================
      IF ( ALLOCATED( AD53         ) ) DEALLOCATE( AD53         ) 
      IF ( ALLOCATED( AD53_PG_PP ) ) DEALLOCATE( AD53_PG_PP)
      IF ( ALLOCATED( AD53_POPG_OH ) ) DEALLOCATE( AD53_POPG_OH )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG53

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG53_MOD
