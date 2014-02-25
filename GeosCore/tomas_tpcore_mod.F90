#if   defined ( TOMAS )
! $Id: tomas_tpcore_mod.f90,v 1.1 2010/02/02 16:57:48 bmy Exp $
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tomas_tpcore_mod
!
! !DESCRIPTION: Module TOMAS\_tpcore\_mod contains routine for GEOS-4 and 
!  GEOS-5 tpcore module (tpcore_fvdas_mod.f90) used for TOMAS aerosol 
!  microphysics.
!
! !INTERFACE: 
!
MODULE tomas_tpcore_mod
!
! !USES:
! 
  IMPLICIT NONE
  PUBLIC
!
! !PRIVATE MEMBER FUNCTIONS:
! 
  PRIVATE :: MASSFLUX
!
! !REVISION HISTORY:
!  22 Jan 2010 - Win T.      - Modified for TOMAS
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
! !IROUTINE: IDFLUX
!
! !DESCRIPTION: \subsection*{Overview}
!  Function IDFLUX return an ID of 0, 1, or 2 which defines if the tracer 
!  should be advected with the TPCORE's original flux computation 
!  or use flux calculated from subroutine MASSFLUX.  This works with the 
!  size-resolved aerosol number and mass. (win, 7/24/07)
!
!\subsection*{More detail}
!  TOMAS aerosol algorithm requires the ratio of aerosol mass over number 
!  (mass per particle) for each size bin to be in size bin boundaries.
!  With PPM advection scheme, if aerosol mass and number tracers are advected 
!  individually, the mass per particle will not be restrained in the 
!  boundaries.
!  Therefore, the aerosol number is chosen to use TPCORE advection fluxes, 
!  the fluxing of aerosol mass tracers will then be computed according to the 
!  number fluxes.
!     IDFLUX = 0 is a non-aerosol tracers --> normal TPCORE calculation
!     IDFLUX = 1 is an aerosol number     --> need TPCORE flux calculation 
!                                             and save flux weight
!     IDFLUX = 2 is an aerosol mass       --> skip flux calculation
!
!\\
!\\
! !INTERFACE:
!
  FUNCTION IDFLUX ( IC ) RESULT(IDF) 
!
! !USES:
!
    USE TOMAS_MOD,    ONLY : IBINS, ICOMP, IDIAG
    USE TRACERID_MOD, ONLY : IDTNK1
!
! !INPUT PARAMETERS: 
!
    INTEGER                :: IC       ! Tracer number
!
! !RETURN VALUE: 
!
    INTEGER                :: IDF      ! Tracer type
!
! !REVISION HISTORY: 
!  24 Jul 2007 - Win T. - Initial version
!  27 Sep 2008 - Win T. - Use IDIAG. Just skip normal flux calculation for
!                         aerosol mass tracer excluding the diagnostic species,
!                         that is aerosol water
!EOP
!------------------------------------------------------------------------------
!BOC
    IDF = 0

    ! Exit if aerosol number is not at work
    IF ( IDTNK1 > 0 ) THEN 
       IF ( IC >= IDTNK1 .AND. IC < IDTNK1 + IBINS ) IDF = 1
       IF ( IC >= IDTNK1 + IBINS .AND.                   &
            IC < IDTNK1 + (IBINS * (ICOMP-IDIAG+1)) )    &
            IDF = 2
    ENDIF

  END FUNCTION IDFLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aermassdr
!
! !DESCRIPTION: Subroutine AERMASSDR is a driver to calculate the 
!  advection of aerosols in X and Y directions for GEOS-4 and GEOS-5
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AERMASSDR( IDFLUX, IM,   JM,    J1, J2, NL, &
                        FLUXNK, DQNK, FLUXM, DQ, DIR     )
!
! !INPUT PARAMETERS: 
!
    INTEGER  :: IM,JM,J1,J2, NL, IDFLUX
    INTEGER  :: DIR
    REAL*8   :: DQ(IM,JM,NL)  ! DQ of aerosol mass tracers, i.e., SF1..SF30.
!
! !INPUT/OUTPUT PARAMETERS
!  
    REAL*8   :: FLUXNK(IM,JM,NL), FLUXM(IM,JM,NL)
    REAL*8   :: DQNK(IM,JM,NL)! the original DQ of aerosol number
!
! !REVISION HISTORY: 
!  24 Jul 2007 - Win T. - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, K

    !=======================================================================
    ! AERMASSDR begins here!
    !=======================================================================

    ! FOR AEROSOL NUMBER:
    ! Before the high-order fluxes are added to DQ array, save the 
    ! current DQ and the fluxes arrays in all three direction
    IF ( IDFLUX == 1 ) THEN
       DO K = 1, NL
       DO J = J1, J2
       DO I = 1, IM
          DQNK(I,J,K)   = DQ(I,J,K)
          FLUXNK(I,J,K) = FLUXM(I,J,K)
       ENDDO
       ENDDO
       ENDDO
    ENDIF

    ! FOR AEROSOL MASS:
    ! Compute the high-order flux 'fx' based on aerosol number
    IF ( IDFLUX == 2 ) THEN
       CALL MASSFLUX(IM,JM,J1,J2,NL,FLUXNK,        &
            DQNK,FLUXM,DQ,DIR )
    ENDIF

  END SUBROUTINE AERMASSDR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: massflux
!
! !DESCRIPTION: Subroutine MASSFLUX computes fluxes of aerosol mass 
!  based on the flux of aerosol number in X and Y directions and 
!  the ratio of DQ_mass/DQ_num.
!  Flux of mass = Flux of number * DQ(aerosol mass)/DQ(aerosol number)
!  Have to choose the DQ for this multiplication according to the upstream 
!  direction of a flux.  (win, 6/8/05, 7/24/07, 7/17/09), (ccc, 8/24/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MASSFLUX ( IM, JM, J1, J2, NL,FLUXNK, DQNK,      &
                             FLUXM, DQ, DIR )
!
! !USES:
!
    USE ERROR_MOD,  ONLY: SAFE_DIV
!
! !INPUT PARAMETERS:
!
    INTEGER  :: IM            ! Total number of E-W grid boxes
    INTEGER  :: JM            ! Total number of N-S grid boxes
    INTEGER  :: J1            ! North edge of the South polar cap
    INTEGER  :: J2            ! South edge of the North polar cap
    INTEGER  :: NL               ! Total number of vertical grid levels
    INTEGER  :: DIR           ! Direction chosen to calculate mass flux
    REAL*8   :: DQ(IM,J1:J2,NL)  ! DQ of aerosol mass tracers, i.e., SF1..SF30.
    REAL*8   :: DQNK(IM,J1:J2,NL)     ! Saved DQ of aerosol number
    REAL*8   :: FLUXNK(IM,J1:J2,NL)   ! Saved flux of aerosol number
!
! !OUTPUT PARAMETERS:
!
    REAL*8   :: FLUXM(IM,J1:J2,NL)    ! Flux of aerosol mass to be calculated
!
!  !REVISION HISTORY: 
!  (1 ) About the size of fluxes array, the high-order fluxes eventually I want
!       fx(im+1,jm,nl), fy(im,jm,nl) and fz(im,jm,nl+1). But here, let's just 
!       use all fluxes array of(im,jm,nl).  I can add the extra strip later 
!       fx(im+1,:,:)=fx(1,:,:)
!       fz(:,:,nl+1)=0d0
!       also fz(:,:,1)=0d0
!  (2 ) dir = 1 for X-direction flux
!           = 2 for Y-direction flux
!  (3 ) Introduce the courant number to massflux for accurate reference of 
!       the upstream box. (win, 3/12/06)
!  (4 ) Modify the subroutine to calculate only in one selected direction
!       with only one set of input argument for the chosen direction
!       (win, 3/16/06)
!  (5 ) Import to new version of G-C and remove courant from the subroutine
!       (win,7/24/07)
!  (6 ) Import to G-C v.8-02-02 (win, 7/17/09)
!  (7 ) Import from tpcore_mod.f to tpcore_fvdas_mod.f90 (ccc, 8/24/09)
!  (8 ) Move the Z-direction to a different routine because arrays must be
!       2D for X and Y directions and 3D for Z direction. (ccc, 8/24/09) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER  :: I, J, K

    !=======================================================================
    ! MASSFLUX begins here!
    !=======================================================================
    
    FLUXM = 0d0          ! Initialize the output flux arrays as zeroes

    ! ## PART 1 -- x-direction fluxes ##
    IF ( DIR == 1 ) THEN
#if defined( GEOS_3 )
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( K, J, I ) 
#endif
       DO K = 1,NL
       DO J = J1, J2

          !Start calculating flux for aerosol mass
          IF ( FLUXNK(1,J,K) >= 0d0 ) THEN
             FLUXM(1,J,K) = SAFE_DIV(FLUXNK(1,J,K)* DQ(IM,J,K),  &
                                     DQNK(IM,J,K), FLUXM(1,J,K))
          ELSE
             FLUXM(1,J,K) = SAFE_DIV(FLUXNK(1,J,K)* DQ(1,J,K),  &
                                     DQNK(1,J,K), FLUXM(1,J,K))
          ENDIF

          DO I = 2, IM
             IF ( FLUXNK(I,J,K) >= 0d0 ) THEN
                FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I-1,J,K),  &
                                        DQNK(I-1,J,K), FLUXM(I,J,K))
             ELSE                  
                FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I,J,K),  &
                                        DQNK(I,J,K), FLUXM(I,J,K))
             ENDIF
          ENDDO               !I-loop

       ENDDO                  !j-loop
       ENDDO                  !k-loop
#if defined( GEOS_3 )
!$OMP END PARALLEL DO 
#endif
    ENDIF                     !DIR = 1

    ! ## PART 2 -- y-direction fluxes ##
    IF ( DIR == 2 ) THEN
#if defined( GEOS_3 )
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( K, J, I ) 
#endif
       DO K = 1,NL
       DO J = J1, J2+1
       DO I = 1, IM
          IF ( FLUXNK(I,J,K) >= 0d0 ) THEN
             FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I,J-1,K),  &
                                     DQNK(I,J-1,K), FLUXM(I,J,K))
          ELSE
             FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I,J,K),  &
                                     DQNK(I,J,K), FLUXM(I,J,K))
          ENDIF
       ENDDO                  !i-loop
       ENDDO                  !j-loop
       ENDDO                  !k-loop
#if defined( GEOS_3 )
!$OMP END PARALLEL DO 
#endif
    ENDIF                     !DIR = 2

    ! ## PART 3 -- z-direction fluxes ##
    IF ( DIR == 3 ) THEN
#if defined( GEOS_3 )
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( K, J, I ) 
#endif
       DO K = 2, NL
       DO J = 1, JM
       DO I = 1, IM
          IF ( FLUXNK(I,J,K) >= 0d0 ) THEN
             FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I,J,K-1),  &
                                     DQNK(I,J,K-1), FLUXM(I,J,K))
          ELSE
             FLUXM(I,J,K) = SAFE_DIV(FLUXNK(I,J,K)* DQ(I,J,K),  &
                                     DQNK(I,J,K), FLUXM(I,J,K))
          ENDIF
       ENDDO                  !i-loop
       ENDDO                  !j-loop
       ENDDO                  !k-loop
#if defined( GEOS_3 )
!$OMP END PARALLEL DO 
#endif
    ENDIF                     !DIR = 2

    ! Return to calling subroutine
  END SUBROUTINE MASSFLUX
!EOC      
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nbin
!
! !DESCRIPTION: Function NBIN returns the bin number of the size-resolved 
!  aerosol according to the current tracers list.
!\\
!\\
! !INTERFACE:
!
  FUNCTION NBIN ( N ) RESULT( NB )
!
! !USES:
! 
    USE TOMAS_MOD,    ONLY : IBINS
    USE TRACERID_MOD, ONLY : IDTNK1
!
! !REVISION HISTORY:
!  24 Jul 2007 - Win T.   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
! 
    INTEGER               :: NB, N, M

    !=======================================================================
    ! NBIN begins here!
    !=======================================================================

    M = N - IDTNK1 + 1
      
    IF ( M > IBINS ) THEN 
       NB = MOD( M, IBINS )
       IF ( NB == 0 ) NB = IBINS
    ELSE
       NB = M
    ENDIF

  END FUNCTION NBIN
!EOC
END MODULE tomas_tpcore_mod
#endif
