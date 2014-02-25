!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag53_mod
!
! !DESCRIPTION: Module DIAG53\_MOD contains arrays and routines for archiving
!  the ND53 diagnostic -- POPS emissions, mass, and production. (eck 9/20/10)
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG53_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PUBLIC, PARAMETER   :: PD53 = 11  ! # of AD53 diags
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: ZERO_DIAG53
      PUBLIC :: WRITE_DIAG53
      PUBLIC :: INIT_DIAG53
      PUBLIC :: CLEANUP_DIAG53
!
! !PUBLIC DATA MEMBERS:
!
      ! Scalars
      INTEGER, PUBLIC              :: ND53  ! ND53 on/off flag
      INTEGER, PUBLIC              :: LD53  ! # of levels

      ! Arrays
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_OC_NEG(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_OC_POS(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_BC_NEG(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_BC_POS(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPG_OH(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_OC_O3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_BC_O3(:,:,:)

!      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPG_NO3(:,:,:)
!      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPG_OX(:,:,:)
!
! !REMARKS:
!  Nomenclature: 
!  ============================================================================
!  (1 ) POPG                  : Gas phase POP
!  (2 ) POPP                  : PARTICULATE PHASE POP
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on DIAG03_MOD
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: zero_diag53
!
! !DESCRIPTION: Subroutine ZERO\_DIAG53 zeroes all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ZERO_DIAG53
!
! !USES:
!
      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N

      !=================================================================
      ! ZERO_DIAG53 begins here!
      !=================================================================

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) RETURN

      ! Zero arrays
      AD53              = 0D0
      AD53_PG_OC_NEG    = 0D0
      AD53_PG_OC_POS    = 0D0
      AD53_PG_BC_NEG    = 0D0
      AD53_PG_BC_POS    = 0D0
      AD53_POPG_OH      = 0D0

      END SUBROUTINE ZERO_DIAG53
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diag53
!
! !DESCRIPTION: Subroutine WRITE\_DIAG53 writes the ND53 diagnostic arrays to 
!  the binary punch file at the proper time.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WRITE_DIAG53
!
! !USES:
!
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
      USE TIME_MOD,     ONLY : GET_CT_CHEM ! CDH for sea salt loss rate

      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_DIAG_MOD       ! TINDEX
!
! !REMARKS:
!   # : Field    : Description                     : Units    : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) PG-SRCE  : POP emissions                   : kg       : 1
!  (2 ) PG-PP-$  : Gas phase POP reacted with OH   : kg       : 1
!                   or partitioned
! 
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
            ARRAY(:,:,1)      = AD53(:,:,N)

         ELSE IF ( N == 5  ) THEN

            !--------------------------------
            ! #2 New gas phase from OC (negative formation of OC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_NEG(:,:,1:LMAX)

         ELSE IF ( N == 6  ) THEN

            !--------------------------------
            ! #3 New OC phase from gas (positive formation of OC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_POS(:,:,1:LMAX)

         ELSE IF ( N == 7  ) THEN

            !--------------------------------
            ! #4 New gas phase from BC (negative formation of BC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_NEG(:,:,1:LMAX)

         ELSE IF ( N == 8 ) THEN

            !--------------------------------
            ! #5 New BC phase from gas (positive formation of BC)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_POS(:,:,1:LMAX)


         ELSE IF ( N == 9  ) THEN

            !--------------------------------
            ! #6 Production of oxidized POPG from rxn with OH (clf, 1/27/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_POPG_OH(:,:,1:LMAX)     

         ELSE IF ( N == 10  ) THEN

            !--------------------------------
            ! #7 Production of oxidized POPOC from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_POPP_OC_O3(:,:,1:LMAX) 

         ELSE IF ( N == 11  ) THEN

            !--------------------------------
            ! #8 Production of oxidized POPBC from rxn with O3 (clf, 6/28/11)
            !--------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg'
            LMAX              = LD53 
            NN                = N
            ARRAY(:,:,1:LMAX) = AD53_POPP_BC_O3(:,:,1:LMAX)       
  
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

      END SUBROUTINE WRITE_DIAG53
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag53
!
! !DESCRIPTION: Subroutine INIT\_DIAG53 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG53
!
! !USES:
!
      USE ERROR_MOD,    ONLY : ALLOC_ERR
   
      USE CMN_SIZE_MOD
! 
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      ALLOCATE( AD53_POPP_OC_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OC_O3' )

      ALLOCATE( AD53_POPP_BC_O3( IIPAR, JJPAR, LD53 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BC_O3' )
    
      ! Zero arrays
      CALL ZERO_DIAG53

      END SUBROUTINE INIT_DIAG53
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_diag53
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG53 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DIAG53
! 
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_DIAG53 begins here!
      !=================================================================
      IF ( ALLOCATED( AD53         ) ) DEALLOCATE( AD53         )
      IF ( ALLOCATED( AD53_PG_OC_NEG ) ) DEALLOCATE( AD53_PG_OC_NEG ) 
      IF ( ALLOCATED( AD53_PG_OC_POS ) ) DEALLOCATE( AD53_PG_OC_POS )
      IF ( ALLOCATED( AD53_PG_BC_NEG ) ) DEALLOCATE( AD53_PG_BC_NEG )
      IF ( ALLOCATED( AD53_PG_BC_POS ) ) DEALLOCATE( AD53_PG_BC_POS )
      IF ( ALLOCATED( AD53_POPG_OH ) ) DEALLOCATE( AD53_POPG_OH )
      IF ( ALLOCATED( AD53_POPP_OC_O3 ) ) DEALLOCATE( AD53_POPP_OC_O3 )
      IF ( ALLOCATED( AD53_POPP_BC_O3 ) ) DEALLOCATE( AD53_POPP_BC_O3 )

      END SUBROUTINE CLEANUP_DIAG53
!EOC
      END MODULE DIAG53_MOD
