!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag53_mod.F90
!
! !DESCRIPTION: Module DIAG53\_MOD contains arrays and routines for archiving
!  the ND53 diagnostic -- POPS emissions, mass, and production. (eck 9/20/10)
!\\
!\\
! !INTERFACE:
!
MODULE DIAG53_MOD
#ifdef BPCH_DIAG
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
!  See https://github.com/geoschem/geos-chem for complete history
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
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! ZERO_DIAG53 begins here!
    !=================================================================

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
  SUBROUTINE WRITE_DIAG53( Input_Opt, State_Grid )
!
! !USES:
!
    USE BPCH2_MOD,      ONLY : BPCH2, GET_MODELNAME, GET_HALFPOLAR
    USE FILE_MOD,       ONLY : IU_BPCH
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_CT_EMIS, GET_DIAGb,  GET_DIAGe
    USE TIME_MOD,       ONLY : GET_CT_CHEM ! CDH for sea salt loss rate
    USE CMN_DIAG_MOD         ! TINDEX
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: CENTER180, HALFPOLAR,   IFIRST
    INTEGER               :: JFIRST,    LFIRST,      LMAX
    INTEGER               :: M,         N,           NN
    REAL(f4)              :: ARRAY(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(f4)              :: LONRES,    LATRES
    REAL(fp)              :: SCALE,     SECONDS
    REAL(f8)              :: DIAGb,     DIAGe
    CHARACTER(LEN=20)     :: MODELNAME
    CHARACTER(LEN=40)     :: CATEGORY,  RESERVED,    UNIT
    REAL(fp)              :: NCHEMSTEP !CDH for sea salt loss rate

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
    IFIRST    = State_Grid%XMinOffset + 1
    JFIRST    = State_Grid%YMinOffset + 1
    LATRES    = State_Grid%DY
    LFIRST    = 1
    LONRES    = State_Grid%DX
    MODELNAME = GET_MODELNAME( Input_Opt, State_Grid )
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
       CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES, &
                   HALFPOLAR, CENTER180, CATEGORY, NN, &
                   UNIT,      DIAGb,     DIAGe,    RESERVED, &
                   State_Grid%NX, State_Grid%NY, LMAX, IFIRST, &
                   JFIRST,    LFIRST,    ARRAY(:,:,1:LMAX) )
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
  SUBROUTINE INIT_DIAG53( State_Grid )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
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
    LD53 = MIN( ND53, State_Grid%NZ )
    
    ! 3-D arrays ("PP-PG-$")
    ALLOCATE( AD53_PG_OC_NEG( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_OC_NEG' )

    ALLOCATE( AD53_PG_OC_POS( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_OC_POS' )

    ALLOCATE( AD53_PG_BC_NEG( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_BC_NEG' )

    ALLOCATE( AD53_PG_BC_POS( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_PG_BC_POS' )

    ALLOCATE( AD53_POPG_OH( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPG_OH' )

    ALLOCATE( AD53_POPP_OCPO_O3( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPO_O3' )

    ALLOCATE( AD53_POPP_OCPI_O3( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPI_O3' )

    ALLOCATE( AD53_POPP_BCPO_O3( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPO_O3' )

    ALLOCATE( AD53_POPP_BCPI_O3( State_Grid%NX, State_Grid%NY, LD53 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPI_O3' )

    ALLOCATE( AD53_POPP_OCPO_NO3( State_Grid%NX, State_Grid%NY, LD53), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPO_NO3' )
    
    ALLOCATE( AD53_POPP_OCPI_NO3( State_Grid%NX, State_Grid%NY, LD53), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_OCPI_NO3' )

    ALLOCATE( AD53_POPP_BCPO_NO3( State_Grid%NX, State_Grid%NY, LD53), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPO_NO3' )

    ALLOCATE( AD53_POPP_BCPI_NO3( State_Grid%NX, State_Grid%NY, LD53), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD53_POPP_BCPI_NO3' )

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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_DIAG53 begins here!
    !=================================================================
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

  END SUBROUTINE CLEANUP_DIAG53
!EOC
#endif
END MODULE DIAG53_MOD
