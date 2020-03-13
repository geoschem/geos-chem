#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndxx_setup.F90
!
! !DESCRIPTION: Subroutine NDXX\_SETUP dynamically allocates memory for
!  certain diagnostic arrays that  are declared allocatable in "diag\_mod.F90".
!\\
!\\
!  This allows us to reduce the amount of memory that needs to be declared
!  globally.  We only allocate memory for arrays if the corresponding
!  diagnostic is turned on.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE NDXX_SETUP( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
  USE CMN_DIAG_MOD
  USE DIAG_MOD
  USE ErrCode_Mod
  USE ERROR_MOD
  USE Input_Opt_Mod,      ONLY : OptInput
  USE State_Chm_Mod,      ONLY : ChmState
  USE State_Chm_Mod,      ONLY : Ind_
  USE State_Grid_Mod,     ONLY : GrdState
  USE State_Met_Mod,      ONLY : MetState
#ifdef TOMAS
  USE TOMAS_MOD,          ONLY : IBINS, ICOMP, IDIAG   !(win, 7/9/09)
#endif

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
  TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
  TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
  INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  16 Jun 1998 - I. Bey, R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER :: NMAX, AS

  !=================================================================
  ! NDXX_SETUP begins here!
  !=================================================================

#ifdef TOMAS
  !=================================================================
  ! ND44: Drydep fluxes [s-1] and drydep velocities [cm/s]
  !       --> uses AD44 arrays (allocatable)
  !=================================================================
  IF ( .not. Input_Opt%LDRYD ) ND44 = 0

  IF ( ND44 > 0 ) THEN
     ! add space in diag array for TOMAS aerosol mass (win, 7/14/09)
     ! Now use State_Chm%nDryDep for # dry depositing species
     ! (ewl, 10/14/15)
     IF ( Ind_('NK1') > 1 ) THEN
        NMAX = State_Chm%nDryDep + ( ICOMP - IDIAG )* IBINS
     ELSE
        NMAX = State_Chm%nDryDep
     ENDIF
     ! Allocate AD44 array
     ALLOCATE( AD44( State_Grid%NX, State_Grid%NY, NMAX, 2 ), STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD44' )
  ENDIF

  !=================================================================
  ! ND59: Size-resolved primary aerosol emissions      !(win, 7/9/09)
  !         Emissions to number, sulfate, sea-salt, carb, dust
  !      ----> save 3-D (I,J,1) or up to (I,J,2)
  !=================================================================
  IF ( ND59 > 0 ) THEN
     LD59 = MIN( ND59, State_Grid%NZ )

     ! Number emission
     ALLOCATE( AD59_NUMB( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_NUMB' )

     ALLOCATE( AD59_SULF( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SULF' )

     ALLOCATE( AD59_SALT( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SALT' )

     ALLOCATE( AD59_ECIL( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECIL' )

     ALLOCATE( AD59_ECOB( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECOB' )

     ALLOCATE( AD59_OCIL( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCIL' )

     ALLOCATE( AD59_OCOB( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCOB' )

     ALLOCATE( AD59_DUST( State_Grid%NX, State_Grid%NY, 2, IBINS ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_DUST' )
  ENDIF

  !=================================================================
  ! ND60: TOMAS microphysical process rates (condensation,
  !           coagulation, nucleation, aqueous oxidation,
  !           error-fudging)
  !       ---> save 2-D (J,L) for 30-bin of each aerosol species
  !=================================================================
  IF ( ND60 > 0 ) THEN
     LD60 = MIN( ND60, State_Grid%NZ )

     ! Now the array dimension is IBINS*(ICOMP-IDIAG+1) because
     ! we need it for all prognostic mass species + 1 for number
     ! IDIAG = # of diagnostic species.  (win, 9/27/08)
     !Condensation rate
     ALLOCATE( AD60_COND(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COND' )

     !Coagulation rate
     ALLOCATE( AD60_COAG(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COAG' )

     !Nucleation rate
     ALLOCATE( AD60_NUCL(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_NUCL' )

     !Aqueous oxidation rate
     ALLOCATE( AD60_AQOX(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_AQOX' )

     !Accumulated error-fudging
     ALLOCATE( AD60_ERROR(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_ERROR' )

     !SOA Condensation rate
     ALLOCATE( AD60_SOA(1,State_Grid%NY,LD60,IBINS*(ICOMP-IDIAG+1)), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_SOA' )
  ENDIF

  !=================================================================
  ! ND61: 3-D TOMAS process rate diagnostic
  !       --> Uses AD61 array (allocatable)
  !  NOTE: ND61 is used for 10-nm particle formation
  !        rate and cluster-size nucleation rate.  So the array
  !        is declared for (NX,NY,LD61,2) (win, 10/6/08)
  !=================================================================
  IF ( ND61 > 0 ) THEN
     LD61 = MIN( ND61, State_Grid%NZ )

     ALLOCATE( AD61( State_Grid%NX, State_Grid%NY, LD61, PD61 ), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61' )
     ALLOCATE( AD61_INST( State_Grid%NX, State_Grid%NY, LD61, PD61), &
               STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61_INST' )
  ENDIF
#endif

#ifdef RRTMG
  !=================================================================
  ! ND72: Radiative output diagnostics (TOASW, SRFSW, TOALW, SRFLW,
  !       GTOASW, GSRFSW, GTOALW, GSRFLW, ATOASW, ASRFSW, ATOALW,
  !       ASRFLW) for both clear and all sky (prefix CLR and ALL)
  !       --> uses AD72 array (allocatable)
  !=================================================================
  IF ( ND72 > 0 ) THEN
     ALLOCATE( AD72( State_Grid%NX, State_Grid%NY, PD72 ), STAT=AS )
     IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD72' )
  ENDIF
#endif

END SUBROUTINE NDXX_SETUP
!EOC
#endif
