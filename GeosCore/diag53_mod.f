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
!  (2 ) AD53_PG_OC   (REAL*4) : Array for archiving POP in OC particles 
!  (3 ) AD53_PG_BC   (REAL*4) : Array for archiving POP on BC particles 
!  (4 ) AD53_PG_G    (REAL*4) : Array for archiving gas phase POP
!  (4 ) AD53_POPG_OH (REAL*4) : Array for POPG oxidized by OH

!
!  Module Routines:
!  ============================================================================
!  (1 ) ZERO_DIAG53           : Sets all module arrays to zero
!  (2 ) WRITE_DIAG53          : Writes data in module arrays to bpch file
!  (3 ) INIT_DIAG53           : Allocates all module arrays
!  (4 ) CLEANUP_DIAG53        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag53_mod.f
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
      INTEGER, PARAMETER   :: PD53 = 29  ! Number of diagnostics for AD53

      ! Arrays
      REAL*4,  ALLOCATABLE :: AD53(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_PG_OC_NEG(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_PG_OC_POS(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_PG_BC_NEG(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_PG_BC_POS(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPG_OH(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_OCPO_O3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_OCPI_O3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_BCPO_O3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_BCPI_O3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_OCPO_NO3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_OCPI_NO3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_BCPO_NO3(:,:,:)
      REAL*4,  ALLOCATABLE :: AD53_POPP_BCPI_NO3(:,:,:)

 
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
!  (1 ) 
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
      AD53               = 0D0
      AD53_PG_OC_NEG     = 0D0
      AD53_PG_OC_POS     = 0D0
      AD53_PG_BC_NEG     = 0D0
      AD53_PG_BC_POS     = 0D0
      AD53_POPG_OH       = 0D0
      AD53_POPP_OCPO_O3  = 0D0
      AD53_POPP_OCPI_O3  = 0D0
      AD53_POPP_BCPO_O3  = 0D0
      AD53_POPP_BCPI_O3  = 0D0
      AD53_POPP_OCPO_NO3 = 0D0
      AD53_POPP_OCPI_NO3 = 0D0
      AD53_POPP_BCPO_NO3 = 0D0
      AD53_POPP_BCPI_NO3 = 0D0
  
      ! Return to calling program
      END SUBROUTINE ZERO_DIAG53

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG53
!
!******************************************************************************
!  Subroutine WRITE_DIAG53 writes the ND53 diagnostic arrays to the binary
!  punch file at the proper time. (bmy, 1/21/05, 2/24/06)
!
!   # : Field    : Description                                   : Units    : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) PG-SRCE  : POP emissions                                 : kg       : 1
!  (2 ) PG-PP-$  : Gas phase POP reacted with OH or partitioned  : kg       : 1
!   
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
         IF ( N == 1 .or. N == 2 .or. N == 3 .or. N == 4 
     &       .or. N == 5 .or. N == 6 .or. N == 7) THEN
               
            !--------------------------------
            ! #1 POP emissions 
            ! (N = 1-5 is Total (including secondary), OC-sorbed, BC-sorbed,
            ! primary gas phase, secondary soil emissions, secondary 
            ! lake emissions, and secondary vegatation emissions respectively [kg])
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'kg'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N)

         ELSE IF ( N == 8) THEN
               
            !--------------------------------
            ! #2 Secondary gas phase positive soil POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE


         ELSE IF ( N == 9  ) THEN

            !--------------------------------
            ! #2 Secondary gas phase negative soil POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 10 ) THEN
               
            !--------------------------------
            ! #2 Secondary gas phase positive lake POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE


         ELSE IF ( N == 11  ) THEN

            !--------------------------------
            ! #2 Secondary gas phase negative lake POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 12 ) THEN
               
            !--------------------------------
            ! #2 Secondary gas phase positive leaf surface POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE


         ELSE IF ( N == 13  ) THEN

            !--------------------------------
            ! #2 Secondary gas phase negative leaf surface POP fluxes
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'ng/m2/day'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 14  ) THEN

            !--------------------------------
            ! #2 fugacity ratios for soils (soil/air)
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'unitless'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 15 ) THEN

            !--------------------------------
            ! #2 fugacity ratios for lakes (water/air)
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'unitless'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 16 ) THEN

            !--------------------------------
            ! #2 fugacity ratios for leaves (leaf surface/air)
            !--------------------------------
            CATEGORY          = 'PG-SRCE'
            UNIT              = 'unitless'
            LMAX              = 1
            NN                = N
            ARRAY(:,:,1)      = AD53(:,:,N) / SCALE

         ELSE IF ( N == 17  ) THEN

            !--------------------------------
            ! #3 New gas phase from OC (negative formation of OC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_NEG(:,:,1:LMAX)

         ELSE IF ( N == 18  ) THEN

            !--------------------------------
            ! #4 New OC phase from gas (positive formation of OC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_POS(:,:,1:LMAX)

         ELSE IF ( N == 19  ) THEN

            !--------------------------------
            ! #5 New gas phase from BC (negative formation of BC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_NEG(:,:,1:LMAX)

         ELSE IF ( N == 20 ) THEN

            !--------------------------------
            ! #6 New BC phase from gas (positive formation of BC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_POS(:,:,1:LMAX)


         ELSE IF ( N == 21  ) THEN

            !--------------------------------
            ! #7 Production of oxidized POPG from rxn with OH (clf, 1/27/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPG_OH(:,:,1:LMAX)     

         ELSE IF ( N == 22  ) THEN

            !--------------------------------
            ! #8 Production of oxidized POPOCPO from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPO_O3(:,:,1:LMAX) 

         ELSE IF ( N == 23  ) THEN

            !--------------------------------
            ! #9 Production of oxidized POPOCPI from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPI_O3(:,:,1:LMAX) 

         ELSE IF ( N == 24  ) THEN

            !--------------------------------
            ! #10 Production of oxidized POPBCPO from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPO_O3(:,:,1:LMAX)    

         ELSE IF ( N == 25  ) THEN

            !--------------------------------
            ! #11 Production of oxidized POPBCPI from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPI_O3(:,:,1:LMAX)  

         ELSE IF ( N == 26  ) THEN

            !--------------------------------
            ! #12 Production of oxidized POPOCPO from rxn with NO3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPO_NO3(:,:,1:LMAX) 

         ELSE IF ( N == 27  ) THEN

            !--------------------------------
            ! #9 Production of oxidized POPOCPI from rxn with NO3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPI_NO3(:,:,1:LMAX) 

         ELSE IF ( N == 28  ) THEN

            !--------------------------------
            ! #10 Production of oxidized POPBCPO from rxn with NO3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPO_NO3(:,:,1:LMAX)    

         ELSE IF ( N == 29  ) THEN

            !--------------------------------
            ! #11 Production of oxidized POPBCPI from rxn with NO3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPI_NO3(:,:,1:LMAX)       
  
         ELSE

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

      ! 2-D array ("PG-SRCE")
      ALLOCATE( AD53( IIPAR, JJPAR, PD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53' )

      ! 3-D arrays ("PP-PG-$")
      ALLOCATE( AD53_PG_OC_NEG( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_OC_NEG' )

      ALLOCATE( AD53_PG_OC_POS( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_OC_POS' )

      ALLOCATE( AD53_PG_BC_NEG( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_BC_NEG' )

      ALLOCATE( AD53_PG_BC_POS( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_BC_POS' )

      ALLOCATE( AD53_POPG_OH( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPG_OH' )

      ALLOCATE( AD53_POPP_OCPO_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPO_O3' )

      ALLOCATE( AD53_POPP_OCPI_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPI_O3' )

      ALLOCATE( AD53_POPP_BCPO_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPO_O3' )

      ALLOCATE( AD53_POPP_BCPI_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPI_O3' )

      ALLOCATE( AD53_POPP_OCPO_NO3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPO_NO3' )

      ALLOCATE( AD53_POPP_OCPI_NO3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPI_NO3' )

      ALLOCATE( AD53_POPP_BCPO_NO3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPO_NO3' )

      ALLOCATE( AD53_POPP_BCPI_NO3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPI_NO3' )

    
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
      ! CLEANUP_DIAG53 begins here!
      !=================================================================
      IF ( ALLOCATED( AD53         ) ) DEALLOCATE( AD53         )
      IF ( ALLOCATED( AD53_PG_OC_NEG ) ) DEALLOCATE( AD53_PG_OC_NEG ) 
      IF ( ALLOCATED( AD53_PG_OC_POS ) ) DEALLOCATE( AD53_PG_OC_POS )
      IF ( ALLOCATED( AD53_PG_BC_NEG ) ) DEALLOCATE( AD53_PG_BC_NEG )
      IF ( ALLOCATED( AD53_PG_BC_POS ) ) DEALLOCATE( AD53_PG_BC_POS )
      IF ( ALLOCATED( AD53_POPG_OH ) ) DEALLOCATE( AD53_POPG_OH )
      IF ( ALLOCATED( AD53_POPP_OCPO_O3)) DEALLOCATE(AD53_POPP_OCPO_O3)
      IF ( ALLOCATED( AD53_POPP_OCPI_O3)) DEALLOCATE(AD53_POPP_OCPI_O3)
      IF ( ALLOCATED( AD53_POPP_BCPO_O3)) DEALLOCATE(AD53_POPP_BCPO_O3)
      IF ( ALLOCATED( AD53_POPP_BCPI_O3)) DEALLOCATE(AD53_POPP_BCPI_O3)
      IF (ALLOCATED( AD53_POPP_OCPO_NO3)) DEALLOCATE(AD53_POPP_OCPO_NO3)
      IF (ALLOCATED( AD53_POPP_OCPI_NO3)) DEALLOCATE(AD53_POPP_OCPI_NO3)
      IF (ALLOCATED( AD53_POPP_BCPO_NO3)) DEALLOCATE(AD53_POPP_BCPO_NO3)
      IF (ALLOCATED( AD53_POPP_BCPI_NO3)) DEALLOCATE(AD53_POPP_BCPI_NO3)

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG53

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG53_MOD
