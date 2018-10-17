!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag53_mod.F
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
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
!
      INTEGER, PUBLIC, PARAMETER   :: PD53 = 29  ! # of AD53 diags
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
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_OC_NEG(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_OC_POS(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_BC_NEG(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_PG_BC_POS(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPG_OH(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_OCPO_O3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_OCPI_O3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_BCPO_O3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_BCPI_O3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_OCPO_NO3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_OCPI_NO3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_BCPO_NO3(:,:,:)
      REAL*4,  PUBLIC, ALLOCATABLE :: AD53_POPP_BCPI_NO3(:,:,:)
!
! !REMARKS:
!  Nomenclature: 
!  ============================================================================
!  (1 ) POPG                  : Gas phase POP
!  (2 ) POPP                  : Particulate phase POP
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on DIAG03_MOD
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!  26 Aug 2014 - M. Sulprizio- Remove AD53 array, POPs emissions are now
!                              tracked by HEMCO
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
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
      !=================================================================
      ! ZERO_DIAG53 begins here!
      !=================================================================

#if defined( BPCH_DIAG )

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) RETURN

      ! Zero arrays
      AD53_PG_OC_NEG     = 0e+0_fp
      AD53_PG_OC_POS     = 0e+0_fp
      AD53_PG_BC_NEG     = 0e+0_fp
      AD53_PG_BC_POS     = 0e+0_fp
      AD53_POPG_OH       = 0e+0_fp
      AD53_POPP_OCPO_O3  = 0e+0_fp
      AD53_POPP_OCPI_O3  = 0e+0_fp
      AD53_POPP_BCPO_O3  = 0e+0_fp
      AD53_POPP_BCPI_O3  = 0e+0_fp
      AD53_POPP_OCPO_NO3 = 0e+0_fp
      AD53_POPP_OCPI_NO3 = 0e+0_fp
      AD53_POPP_BCPO_NO3 = 0e+0_fp
      AD53_POPP_BCPI_NO3 = 0e+0_fp
#endif

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
#if defined( BPCH_DIAG )
      USE BPCH2_MOD,    ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
#endif
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GC_GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,     ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
      USE TIME_MOD,     ONLY : GET_CT_CHEM ! CDH for sea salt loss rate

      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_DIAG_MOD       ! TINDEX
!
! !REMARKS:
!   # : Field    : Description                     : Units    : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) PG-SRCE  : POP emissions                   : kg/s     : 1
!  (2 ) PG-PP-$  : Gas phase POP reacted with OH   : kg/s     : SCALE
!                   or partitioned
! 
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  27 Nov 2012 - M. Payer    - Added ProTeX headers
!  26 Aug 2014 - M. Sulprizio- Now track POPs emissions through HEMCO and
!                              write to bpch in diag3.F
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  15 Oct 2018 - R. Yantosca - Data for #17-#29 are now divided by seconds
!                              in pops_mod.F. Divide by the # of chemistry
!                              timesteps. This is done to better match the
!                              netCDF output.  BPCH will be removed soon.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: CENTER180, HALFPOLAR,   IFIRST
      INTEGER               :: JFIRST,    LFIRST,      LMAX
      INTEGER               :: M,         N,           NN
      REAL(f4)              :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL(f4)              :: LONRES,    LATRES
      REAL(fp)              :: SCALE,     SECONDS
      REAL(f8)              :: DIAGb,     DIAGe       
      CHARACTER(LEN=20)     :: MODELNAME 
      CHARACTER(LEN=40)     :: CATEGORY,  RESERVED,    UNIT
      REAL(fp)              :: NCHEMSTEP !CDH for sea salt loss rate

      !=================================================================
      ! WRITE_DIAG53 begins here!
      !=================================================================
#if defined( BPCH_DIAG )

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
      SCALE     = DBLE( GET_CT_EMIS() ) + 1e-32_fp
      NCHEMSTEP = DBLE( GET_CT_CHEM() ) + TINY( 1e+0_fp )
      SECONDS   = ( DIAGe - DIAGb ) * 3600e+0_fp

      !=================================================================
      ! Write data to the bpch file
      !=================================================================

      ! Loop over ND53 diagnostic tracers
      DO M = 1, TMAX(53)

         ! Get ND53 tracer #
         N = TINDEX(53,M)

         !--------------------------------------------------------------
         ! POPs emissions
         !
         ! N = 1-16 are now tracked by HEMCO and saved out in diag3.F
         !--------------------------------------------------------------

         IF ( N == 17  ) THEN

            !-----------------------------------------------------------
            ! New gas phase from OC (negative formation of OC)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_NEG(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 18  ) THEN

            !-----------------------------------------------------------
            ! New OC phase from gas (positive formation of OC)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_OC_POS(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 19  ) THEN

            !-----------------------------------------------------------
            ! New gas phase from BC (negative formation of BC)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_NEG(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 20 ) THEN

            !-----------------------------------------------------------
            ! New BC phase from gas (positive formation of BC)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_PG_BC_POS(:,:,1:LMAX) / SCALE


         ELSE IF ( N == 21  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPG from rxn with OH (clf, 1/27/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPG_OH(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 22  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPOCPO from rxn with O3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPO_O3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 23  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPOCPI from rxn with O3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPI_O3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 24  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPBCPO from rxn with O3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPO_O3(:,:,1:LMAX) / SCALE 

         ELSE IF ( N == 25  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPBCPI from rxn with O3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPI_O3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 26  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPOCPO from rxn with NO3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPO_NO3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 27  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPOCPI from rxn with NO3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_OCPI_NO3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 28  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPBCPO from rxn with NO3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPO_NO3(:,:,1:LMAX) / SCALE

         ELSE IF ( N == 29  ) THEN

            !-----------------------------------------------------------
            ! Production of oxidized POPBCPI from rxn with NO3 (clf, 6/28/11)
            !-----------------------------------------------------------
            CATEGORY          = 'PG-PP-$'
            UNIT              = 'kg/s'
            LMAX              = LD53 
            NN                = N-16
            ARRAY(:,:,1:LMAX) = AD53_POPP_BCPI_NO3(:,:,1:LMAX) / SCALE
  
         ELSE

            ! Otherwise skip to next N
            CYCLE

         ENDIF

         ! Write data to disk
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LMAX,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )
      ENDDO
#endif

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
#if defined( BPCH_DIAG )

      ! Exit if ND53 is turned off
      IF ( ND53 == 0 ) THEN
         LD53 = 0
         RETURN
      ENDIF

      ! Get number of levels for 3-D arrays
      LD53 = MIN( ND53, LLPAR )

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
#endif

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
#if defined( BPCH_DIAG )
      IF (ALLOCATED(AD53_PG_OC_NEG    )) DEALLOCATE(AD53_PG_OC_NEG    ) 
      IF (ALLOCATED(AD53_PG_OC_POS    )) DEALLOCATE(AD53_PG_OC_POS    )
      IF (ALLOCATED(AD53_PG_BC_NEG    )) DEALLOCATE(AD53_PG_BC_NEG    )
      IF (ALLOCATED(AD53_PG_BC_POS    )) DEALLOCATE(AD53_PG_BC_POS    )
      IF (ALLOCATED(AD53_POPG_OH      )) DEALLOCATE(AD53_POPG_OH      )
      IF (ALLOCATED(AD53_POPP_OCPO_O3 )) DEALLOCATE(AD53_POPP_OCPO_O3 )
      IF (ALLOCATED(AD53_POPP_OCPI_O3 )) DEALLOCATE(AD53_POPP_OCPI_O3 )
      IF (ALLOCATED(AD53_POPP_BCPO_O3 )) DEALLOCATE(AD53_POPP_BCPO_O3 )
      IF (ALLOCATED(AD53_POPP_BCPI_O3 )) DEALLOCATE(AD53_POPP_BCPI_O3 )
      IF (ALLOCATED(AD53_POPP_OCPO_NO3)) DEALLOCATE(AD53_POPP_OCPO_NO3)
      IF (ALLOCATED(AD53_POPP_OCPI_NO3)) DEALLOCATE(AD53_POPP_OCPI_NO3)
      IF (ALLOCATED(AD53_POPP_BCPO_NO3)) DEALLOCATE(AD53_POPP_BCPO_NO3)
      IF (ALLOCATED(AD53_POPP_BCPI_NO3)) DEALLOCATE(AD53_POPP_BCPI_NO3)
#endif

      END SUBROUTINE CLEANUP_DIAG53
!EOC
      END MODULE DIAG53_MOD
