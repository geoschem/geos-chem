!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_Br_mod.F90
!
! !DESCRIPTION: Module GLOBAL\_Br\_MOD contains variables and routines for
!  storing the global monthly mean Br concentration evaluated from HEMCO.
!\\
!\\
! !INTERFACE:
!
MODULE GLOBAL_Br_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
  ! Array to store global monthly mean BR field
  REAL(fp), PUBLIC, ALLOCATABLE :: BR_TROP(:,:,:)
  REAL(fp), PUBLIC, ALLOCATABLE :: BR_STRAT(:,:,:)
  REAL(fp), PUBLIC, ALLOCATABLE :: BR_MERGE(:,:,:)

  ! Array to store global monthly mean BrO field
  REAL(fp), PUBLIC, ALLOCATABLE :: BRO_TROP(:,:,:)
  REAL(fp), PUBLIC, ALLOCATABLE :: BRO_STRAT(:,:,:)
  REAL(fp), PUBLIC, ALLOCATABLE :: BRO_MERGE(:,:,:)

  ! Array to store global monthly J-BrO field
  REAL(fp), PUBLIC, ALLOCATABLE :: J_BRO(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GET_GLOBAL_Br
  PUBLIC :: INIT_GLOBAL_Br
  PUBLIC :: CLEANUP_GLOBAL_Br
!
! !REMARKS:
!  References:
!  (1 ) Holmes, C. D., et al. (2006), Global lifetime of elemental mercury
!       against oxidation by atomic bromine in the free troposphere, Geophys.
!       Res. Lett., 33(20).
!  (2 ) Holmes, C.D., et al. (2010) Global atmospheric model for mercury
!       including oxidation by bromine atoms, AC&P, 10, 12,037-12,057.
!  (3 ) Parrella, J. et al. (2012), Tropospheric bromine chemistry:
!       implications for present and pre-industrial ozone and mercury, ACP.
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes - Initial version
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
! !IROUTINE: get_global_Br
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_Br evaluates global Br from
!  data stored in HEMCO, including application of any masks or scaling
!  set in the HEMCO configuration file. This Br data is needed as oxidant
!  for mercury chemistry.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_GLOBAL_Br( Input_Opt, State_Grid, State_Met, THISMONTH, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,      ONLY : OptInput
    USE OCEAN_MERCURY_MOD,  ONLY : LGCBROMINE     !eds 4/19/12
    USE State_Met_Mod,      ONLY : MetState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L
    INTEGER            :: TPL
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! GET_GLOBAL_BR begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at GET_GLOBAL_Br (in GeosCore/global_br_mod.F)'

    ! Allocate BR array, if this is the first call
    IF ( FIRST ) THEN

       ! Initialize arrays
       CALL INIT_GLOBAL_BR( State_Grid, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in call to "INIT_GLOBAL_BR"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    IF ( LGCBROMINE ) THEN

       !-----------------------------------------------------------------
       ! Evaluate Br_GC from GEOS-Chem to set Br_MERGE (trop+strat)
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'Br_GC', Br_MERGE, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find Br_GC in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Convert v/v -> pptv
       Br_MERGE = Br_MERGE * 1e+12_fp

       !-----------------------------------------------------------------
       ! Evaluate BrO_GC from GEOS-Chem to set BrO_MERGE (trop+strat)
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'BrO_GC', BrO_MERGE, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find BrO_GC in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Convert v/v -> pptv
       BrO_MERGE = BrO_MERGE * 1e+12_fp

    ELSE

       !-----------------------------------------------------------------
       ! Evaluate Br from pTOMCAT biogenic bromocarbons [pptv]
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'Br_TOMCAT', Br_TROP, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find Br_TOMCAT in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! Evaluate BrO from pTOMCAT biogenic bromocarbons [pptv]
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'BrO_TOMCAT', Br_TROP, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find BrO_TOMCAT in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! Evaluate Br from GMI for stratosphere [pptv]
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'Br_GMI', Br_STRAT, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find Br_GMI in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! Evaluate BrO from GMI for stratosphere [pptv]
       !-----------------------------------------------------------------
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'BrO_GMI', BrO_STRAT, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find BrO_GMI in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! Use pTOMCAT exclusively in the troposphere.
       ! In the stratosphere, use the greater value from either COMBO or
       ! the tropospheric model. COMBO source gases include CH3Br and
       ! halons, while pTOMCAT includes CH3Br and shorter-lived gases.
       !-----------------------------------------------------------------
       BR_MERGE  = BR_TROP
       BRO_MERGE = BRO_TROP

       !$OMP PARALLEL DO        &
       !$OMP DEFAULT( SHARED )  &
       !$OMP PRIVATE( I, J, TPL )
       DO J=1, State_Grid%NY
       DO I=1, State_Grid%NX

          ! First layer in the stratosphere
          TPL = State_Met%TropLev(I,J)

          BR_MERGE(I,J,TPL:State_Grid%NZ) = MERGE(      &
               BR_STRAT(I,J,TPL:State_Grid%NZ),         &
               BR_TROP(I,J,TPL:State_Grid%NZ),          &
               MASK = BR_STRAT(I,J,TPL:State_Grid%NZ) > &
                      BR_TROP(I,J,TPL:State_Grid%NZ) )

          BRO_MERGE(I,J,TPL:State_Grid%NZ) = MERGE(     &
               BRO_STRAT(I,J,TPL:State_Grid%NZ),        &
               BRO_TROP(I,J,TPL:State_Grid%NZ),         &
               MASK = BR_STRAT(I,J,TPL:State_Grid%NZ) > &
                      BR_TROP(I,J,TPL:State_Grid%NZ) )

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    !-----------------------------------------------------------------
    ! Evaluate J_BrO
    !-----------------------------------------------------------------
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'JBrO', J_BrO, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find JBrO in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Only set values up to max chemistry level
    IF ( State_Grid%MaxChemLev < State_Grid%NZ ) THEN
       J_BrO(:,:,State_Grid%MaxChemLev+1:State_Grid%NZ) = 0e+0_fp 
    ENDIF

  END SUBROUTINE GET_GLOBAL_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_Br
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_Br allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_GLOBAL_Br( State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_GLOBAL_BR begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Global_Br (in module GeosCore/global_br_mod.F90)'

    !-------------------------------------
    ! Br Arrays
    !-------------------------------------

    ! Allocate BR_TROP array
    ALLOCATE( BR_TROP( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BR_TROP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BR_TROP = 0e+0_fp

    ! Allocate BR_STRAT array
    ALLOCATE( BR_STRAT( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BR_STRAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BR_STRAT = 0e+0_fp

    ! Allocate BR_MERGE array
    ALLOCATE( BR_MERGE( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BR_MERGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BR_MERGE = 0e+0_fp

    !-------------------------------------
    ! BrO Arrays
    !-------------------------------------

    ! Allocate J_BrO array
    ALLOCATE( J_BrO( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:J_BrO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    J_BrO = 0e+0_fp

    ! Allocate BrO_TROP array
    ALLOCATE( BrO_TROP( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BRO_TROP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BrO_TROP = 0e+0_fp

    ! Allocate BrO_STRAT array
    ALLOCATE( BrO_STRAT( State_Grid%NX, State_Grid%NY, State_Grid%NZ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BRO_STRAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BrO_STRAT = 0e+0_fp

    ! Allocate BrO_MERGE array
    ALLOCATE( BrO_MERGE( State_Grid%NX, State_Grid%NY, State_Grid%NZ), &
              STAT=RC )
    CALL GC_CheckVar( 'global_br_mod.F:BrO_MERGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BrO_MERGE = 0e+0_fp

  END SUBROUTINE INIT_GLOBAL_BR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_Br
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_Br deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_GLOBAL_Br( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate variables
    IF ( ALLOCATED( Br_TROP ) ) THEN
       DEALLOCATE( Br_TROP, STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:Br_TROP', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Br_STRAT ) ) THEN
       DEALLOCATE( Br_STRAT, STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:Br_STRAT', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Br_MERGE ) ) THEN
       DEALLOCATE( Br_MERGE , STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:Br_MERGE ', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( J_BrO ) ) THEN
       DEALLOCATE( J_BrO, STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:J_BrO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( BrO_TROP ) ) THEN
       DEALLOCATE( BrO_TROP, STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:BrO_TROP', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( BrO_STRAT ) ) THEN
       DEALLOCATE( BrO_STRAT, STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:BrO_STRAT', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( BrO_MERGE ) ) THEN
       DEALLOCATE( BrO_MERGE , STAT=RC )
       CALL GC_CheckVar( 'global_br_mod.F90:BrO_MERGE ', 2, RC )
       RETURN
    ENDIF

  END SUBROUTINE CLEANUP_GLOBAL_Br
!EOC
END MODULE GLOBAL_Br_MOD
