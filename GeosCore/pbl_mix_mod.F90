!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pbl_mix_mod.F90
!
! !DESCRIPTION: Module PBL\_MIX\_MOD contains routines and variables used to
!  compute the planetary boundary layer (PBL) height and to mix tracers
!  underneath the PBL top.
!\\
!\\
! !INTERFACE:
!
MODULE Pbl_Mix_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_Full_Pbl_Mixing
  PUBLIC  :: Compute_Pbl_Height
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: TurbDay
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
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
! !IROUTINE: do_pbl_mix
!
! !DESCRIPTION: Subroutine DO\_PBL\_MIX is the driver routine for planetary
!  boundary layer mixing.  The PBL layer height and related quantities are
!  always computed.  Complete mixing of tracers underneath the PBL top is
!  toggled by the DO\_TURBDAY switch.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Full_Pbl_Mixing( Input_Opt,  State_Chm, State_Diag,          &
                                 State_Grid, State_Met, RC                  )
!
! !USES:
!
    USE ErrCode_Mod
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE Time_Mod,        ONLY : Get_Ts_Conv
    USE Time_Mod,        ONLY : Get_Ts_Dyn
    USE Timers_Mod,      ONLY : Timer_End, Timer_Start
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    INTEGER,        INTENT(INOUT) :: RC          ! Return code
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N
    INTEGER            :: NA
    INTEGER            :: TS_Dyn
    INTEGER            :: previous_units
    REAL(f8)           :: DT_Dyn

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=======================================================================
    ! Do_Full_Pbl_Mixing begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_PBL_MIX (in module GeosCore/pbl_mix_mod.F90)'

    !========================================================================
    ! Mixing budget diagnostics - Part 1 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get initial column masses (full, trop, PBL)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetMixingFull,               &
            diagFull    = NULL(),                                            &
            mapDataFull = State_Diag%Map_BudgetMixingFull,                   &
            isTrop      = State_Diag%Archive_BudgetMixingTrop,               &
            diagTrop    = NULL(),                                            &
            mapDataTrop = State_Diag%Map_BudgetMixingTrop,                   &
            isPBL       = State_Diag%Archive_BudgetMixingPBL,                &
            diagPBL     = NULL(),                                            &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            isLevs      = State_Diag%Archive_BudgetMixingLevs,               &
            diagLevs    = NULL(),                                            &
            mapDataLevs = State_Diag%Map_BudgetMixingLevs,                   &
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Proceed to do full PBL mixing only if it has been selected in
    ! geoschem_config.yml
    IF ( Input_Opt%LTURB .and. ( .not. Input_Opt%LNLPBL ) ) THEN

       !=====================================================================
       ! Unit conversion #1
       !=====================================================================

       ! Halt mixing timer (so that unit conv can be timed separately)
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "Boundary layer mixing", RC )
       ENDIF

       ! Convert species to [v/v dry] aka [mol/mol dry]
       CALL Convert_Spc_Units(                                               &
            Input_Opt      = Input_Opt,                                      &
            State_Chm      = State_Chm,                                      &
            State_Grid     = State_Grid,                                     &
            State_Met      = State_Met,                                      &
            mapping        = State_Chm%Map_Advect,                           &
            new_units      = MOLES_SPECIES_PER_MOLES_DRY_AIR,                &
            previous_units = previous_units,                                 &
            RC             = RC                                             )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (to mol/mol dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Start mixing timer again
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "Boundary layer mixing", RC )
       ENDIF

       !=====================================================================
       ! Do full PBL mixing
       !=====================================================================

       ! Do complete mixing of tracers in the PBL
       CALL TurbDay( Input_Opt,  State_Chm, State_Diag,                      &
                     State_Grid, State_Met, RC                              )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "TURBDAY"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !=====================================================================
       ! Unit conversion #2
       !=====================================================================

       ! Halt mixing timer (so that unit conv can be timed separately)
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "Boundary layer mixing", RC )
       ENDIF

       ! Convert species back to original units
       CALL Convert_Spc_Units(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            State_Met  = State_Met,                                          &
            mapping    = State_Chm%Map_Advect,                               &
            new_units  = previous_units,                                     &
            RC         = RC                                                 )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (from mol/mol dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Start mixing timer again
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "Boundary layer mixing", RC )
       ENDIF
    ENDIF

    !========================================================================
    ! Full PBL mixing budget diagnostics - Part 2 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get dynamics timestep [s]
       TS_Dyn = Get_Ts_Dyn()
       DT_Dyn = DBLE( TS_Dyn )

       ! Compute change in column masses (after mixing - before mixing)
       ! and store in diagnostic arrays.  Units are [kg/s].
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetMixingFull,               &
            diagFull    = State_Diag%BudgetMixingFull,                       &
            mapDataFull = State_Diag%Map_BudgetMixingFull,                   &
            isTrop      = State_Diag%Archive_BudgetMixingTrop,               &
            diagTrop    = State_Diag%BudgetMixingTrop,                       &
            mapDataTrop = State_Diag%Map_BudgetMixingTrop,                   &
            isPBL       = State_Diag%Archive_BudgetMixingPBL,                &
            diagPBL     = State_Diag%BudgetMixingPBL,                        &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            isLevs      = State_Diag%Archive_BudgetMixingLevs,               &
            diagLevs    = State_Diag%BudgetMixingLevs,                       &
            mapDataLevs = State_Diag%Map_BudgetMixingLevs,                   &
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = DT_Dyn,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Do_Full_Pbl_Mixing
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_pbl_height
!
! !DESCRIPTION: Subroutine COMPUTE\_PBL\_HEIGHT computes the PBL height and
!  other related quantities.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Pbl_Height( Input_Opt, State_Grid, State_Chm,  &
                                 State_Met, State_Diag, RC )
!
! !USES:
!
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE PhysConstants,  ONLY : Scale_Height, Rd, g0
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Met_Mod,   ONLY : MetState
    USE State_Diag_Mod,  ONLY : DgnState
    USE Time_Mod,        ONLY : Get_TS_Dyn
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
<<<<<<< HEAD
    LOGICAL  :: Bad_Sum
    INTEGER  :: I,      J,      L,    LTOP
    REAL(fp) :: Lower_Edge_Height
=======
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    LOGICAL            :: Bad_Sum
    INTEGER            :: I,      J,      L,    LTOP, TS_Dyn
    REAL(fp)           :: BLTOP,  BLTHIK, DELP
    REAL(f8)           :: DT_Dyn
>>>>>>> cc64cd801 (Compute mixing PBL budget diag for PBL ht change in compute_pbl_height)

    ! Arrays
    REAL(fp) :: P(0:State_Grid%NZ)

    !=================================================================
    ! COMPUTE_PBL_HEIGHT begins here!
    !=================================================================

    ! Initialize
<<<<<<< HEAD
    RC                       = GC_SUCCESS
    Bad_Sum                  = .FALSE.
    State_Met%InPbl          = .FALSE.
    State_Met%F_of_PBL       = 0.0_fp
    State_Met%F_Under_PBLTop = 0.0_fp

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, LTOP, Lower_Edge_Height                         )&
    !$OMP COLLAPSE( 2                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 8                                              )
=======
    RC              = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Compute_PBL_Height (in module GeosCore/pbl_mix_mod.F90)'

    Bad_Sum         = .FALSE.
    State_Met%InPbl = .FALSE.

    !------------------------------------------------------------------------
    ! Change in PBL mass for use with budget mixing PBL diagnostic - 1 of 2
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_BudgetMixingPBL ) THEN

       ! Get initial column masses (full, trop, PBL)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = .FALSE.,                                           &
            diagFull    = NULL(),                                            &
            mapDataFull = NULL(),                                            &
            isTrop      = .FALSE.,                                           &
            diagTrop    = NULL(),                                            &
            mapDataTrop = NULL(),                                            &
            isPBL       = .TRUE.,                                            &
            diagPBL     = State_Diag%BudgetMixingPBLHeight,                  &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing PBL height budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( I, J, L, P, BLTOP, BLTHIK, LTOP, DELP )
>>>>>>> cc64cd801 (Compute mixing PBL budget diag for PBL ht change in compute_pbl_height)
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize variables
       State_Met%PBL_Top_m(I,J) = State_Met%PBLH(I,J) ! Pbl ht [m] above sfc
       LTOP                     = 0                   ! Level w/ PBL top
       Lower_Edge_Height        = 0.0_fp              ! Lower edge height
                                                      !  above sfc [m]

       ! Find PBL top level (L) and pressure (hPa)
       DO L = 1, State_Grid%NZ

         !-------------------------------------------------------------------
         ! The PBL top occurs in this level if the condition is true
         !-------------------------------------------------------------------
         IF ( Lower_Edge_Height + State_Met%BXHEIGHT(I,J,L) >=               &
              State_Met%PBLH(I,J)                                ) THEN

            ! PBL top is in this level
            LTOP = L

            ! Pressure at the PBL top altitude, hPa
            ! Use pressure lapse equation:
            !    p(PBLH) = p(z1) * exp( -(PBLH-z1) / Scale_Height )
            !
            ! p(z1) = State_Met%PEDGE(I,J,L) = Pressure @ lower level edge
            !
            ! PBLH - z1 = (State_Met%PBLH(I,J,L) - Lower_Edge_Height) =
            !   Height above the lower level edge
            !
            ! Scale_Height = Rd * Tv / g0
            State_Met%PBL_Top_hPa(I,J) = State_Met%PEdge(I,J,L)     *        &
                  EXP( -( State_Met%PBLH(I,J) - Lower_Edge_Height ) *        &
                        g0 / ( Rd * State_Met%TV(I,J,L) )            )

            ! Fraction of PBL mass in layer L, will be normalized below
            State_Met%F_of_PBL(I,J,L) = State_Met%PEdge(I,J,L)               &
                                      - State_Met%PBL_Top_hPa(I,J)
      
            ! Fraction of the grid cell mass under PBL top
            State_Met%F_Under_PBLTop(I,J,L) = State_Met%F_of_PBL(I,J,L) /    &
                 ( State_Met%PEdge(I,J,L) - State_Met%PEdge(I,J,L+1) )

            ! Model level of PBL top (integer+fraction).
            ! The top is within level CEILING(PBL_Top_L)
            State_Met%PBL_Top_L(I,J) = ( LTOP - 1 )                          &
                                     + State_Met%F_Under_PBLTop(I,J,L)

            ! PBL Thickness from surface to top, hPa
            State_Met%PBL_Thick(I,J) = State_Met%PEdge(I,J,1)                &
                                     - State_Met%PBL_Top_hPa(I,J)

            ! Exit Do loop after we found PBL top level
            EXIT
         ENDIF

         !-------------------------------------------------------------------
         ! The PBL top does not occur in this level.
         ! Update variables and go to next level.
         !-------------------------------------------------------------------

         ! Grid cell fully within PBL
         State_Met%inPBL(I,J,L) = .True.

         ! Fraction of the grid cell mass under PBL top
         State_Met%F_Under_PBLTop(I,J,L) = 1.0_fp

         ! Fraction of PBL mass in layer L, will be normalized below
         State_Met%F_of_PBL(I,J,L) = State_Met%PEdge(I,J,L)                  &
                                   - State_Met%PEdge(I,J,L+1)

         ! Update lower edge height, m
         Lower_Edge_Height = Lower_Edge_Height + State_Met%BXHeight(I,J,L)

       ENDDO

       ! Fraction of PBL mass in layer L, now normalize to sum of 1
       State_Met%F_of_PBL(I,J,:) = State_Met%F_of_PBL(I,J,:)                 &
                                 / State_Met%PBL_Thick(I,J)

       ! Error check
       IF ( ABS( SUM( State_Met%F_OF_PBL(I,J,:) ) - 1.0_fp) > 1.0e-3_fp) THEN
          !$OMP CRITICAL
          PRINT*, 'bad sum at: ', I, J, SUM( State_Met%F_OF_PBL(I,J,:) )
          Bad_Sum = .TRUE.
          !$OMP END CRITICAL
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Exit to main program level if bad sum was encountered
    IF ( Bad_Sum ) THEN
       CALL GC_Error( 'Error in computing F_OF_PBL !', RC, &
                      'COMPUTE_PBL_HEIGHT ("pbl_mix_mod.F90")' )
       RETURN
    ENDIF

    ! Model level where PBL top occurs
    State_Met%PBL_MAX_L = MAXVAL( CEILING( State_Met%PBL_Top_L ) )

    !------------------------------------------------------------------------
    ! Change in PBL mass for use with budget mixing PBL diagnostic - 2 of 2
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_BudgetMixingPBL ) THEN

       ! Dynamic timestep [s]
       TS_Dyn = GET_TS_DYN()
       DT_Dyn = DBLE( TS_Dyn )

       ! Compute change in column masses (after PBL ht change minus before)
       ! and store in diagnostic arrays.  Units are [kg/s].
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = .FALSE.,                                           &
            diagFull    = NULL(),                                            &
            mapDataFull = NULL(),                                            &
            isTrop      = .FALSE.,                                           &
            diagTrop    = NULL(),                                            &
            mapDataTrop = NULL(),                                            &
            isPBL       = .TRUE.,                                            &
            diagPBL     = State_Diag%BudgetMixingPBLHeight,                  &
            mapDataPBL  = State_Diag%Map_BudgetMixingPBL,                    &
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = DT_Dyn,                                            &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing PBL height budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Compute_Pbl_Height
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: turbday
!
! !DESCRIPTION: !  Subroutine TURBDAY executes the GEOS-Chem boundary layer
!  mixing algorithm (full PBL mixing).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TurbDay( Input_Opt,  State_Chm, State_Diag,                     &
                      State_Grid, State_Met, RC                             )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AIRMW
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Get_Ts_Conv
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options Object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Metoerology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object

!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Original subroutine by Dale Allen, Univ of MD.
!
! !REVISION HISTORY:
!  30 Jan 1998 - I. Bey, R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: I,    J,  L
    INTEGER            :: LTOP, N,  NA
    REAL(fp)           :: AA,   CC, CC_AA, DTCONV

    ! Arrays
    REAL(fp)           :: DTC
    REAL(fp)           :: A(State_Grid%NX,State_Grid%NY)
    REAL(fp)           :: FPBL(State_Grid%NX,State_Grid%NY)
    INTEGER            :: IMIX(State_Grid%NX,State_Grid%NY)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    REAL(fp),      POINTER  :: AD  (:,:,:  )
    TYPE(SpcConc), POINTER  :: TC  (:      )

    !========================================================================
    ! TURBDAY begins here!
    !========================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! First-time initialization
    IF ( FIRST .and. Input_Opt%amIRoot ) THEN

       ! Echo info
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) 'T U R B D A Y  -- by Dale Allen, U. Md.'
       WRITE( 6, '(a)' ) 'Adapted for GEOS-Chem by the GCST'
       WRITE( 6, '(a)' ) 'Last Modification Date: 15 May 2020'
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )

       ! Reset first time flag
       FIRST = .FALSE.
    ENDIF

    !========================================================================
    ! Do the boundary layer mixing
    !========================================================================

    ! Initalize
    AD      => State_Met%AD         ! Dry air mass
    TC      => State_Chm%Species    ! Chemical species [v/v]
    IMIX    = ceiling( State_Met%PBL_Top_L ) ! Integer level where PBL top occurs
    FPBL    = State_Met%PBL_Top_L - (IMIX-1) ! Fractional level above IMIX to PBL top

    ! Convection timestep [s]
    DTCONV = GET_TS_CONV()

    ! Loop over Lat/Long grid boxes (I,J)
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, NA, N, AA, CC, CC_AA, DTC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! We assume full mixing in the boundary layer, so the A
       ! coefficients are 1 everywhere, day & night (bmy, 2/11/03)
       A(I,J) = 1e+0_fp

       ! Calculate air mass within PBL at grid box (I,J,L)
       AA = 0.e+0_fp
       DO L = 1, IMIX(I,J)-1
          AA = AA + AD(I,J,L)
       ENDDO

       L  = IMIX(I,J)
       AA = AA + AD(I,J,L) * FPBL(I,J)

       ! Loop over only the advected species
       DO NA = 1, State_Chm%nAdvect

          ! Species ID
          N = State_Chm%Map_Advect(NA)

          !===========================================================
          ! Calculate tracer mass within PBL at grid box (I,J,L)
          !===========================================================

          ! Sum mass from (I,J,L) below the PBL top
          CC = 0.e+0_fp
          DO L = 1, IMIX(I,J)-1
             CC = CC + AD(I,J,L) * TC(N)%Conc(I,J,L)
          ENDDO

          ! Then also sum mass from (I,J,L) which straddle the PBL top
          L     = IMIX(I,J)
          CC    = CC + AD(I,J,L) * TC(N)%Conc(I,J,L) * FPBL(I,J)

          ! CC/AA is the mean mixing ratio of tracer at
          ! (I,J) from L=1 to L=LTOP
          CC_AA = CC / AA

          !========================================================
          ! TC(N)%Conc(I,J,L) new  = TC(N)%Conc(I,J,L) old +
          !                          ( DTC / AD(I,J,L) )
          !
          ! where
          !
          ! DTC = [ alpha * (mean MR below PBL) *
          !                  Airmass at (I,J,L) ] -
          !                [ alpha * TC(N)%Conc(I,J,L) old     *
          !                  Airmass at (I,J,L) ]
          !
          ! DTC is thus the change in mass (kg) due to BL mixing,
          ! so DTC/AD is the change in (V/V) mixing ratio units.
          !========================================================

          ! For grid boxes (I,J,L) which lie below the PBL top
          DO L = 1, IMIX(I,J)-1
             DTC = ( A(I,J) * CC_AA       * AD(I,J,L)  - &
                     A(I,J) * TC(N)%Conc(I,J,L) * AD(I,J,L) )

             TC(N)%Conc(I,J,L) = TC(N)%Conc(I,J,L) + DTC / AD(I,J,L)
          ENDDO

          ! For grid boxes (I,J,L) which straddle the PBL top
          L = IMIX(I,J)

          DTC = ( A(I,J) * FPBL(I,J)  * CC_AA       * AD(I,J,L) - &
                  A(I,J) * FPBL(I,J)  * TC(N)%Conc(I,J,L) * AD(I,J,L) )

          TC(N)%Conc(I,J,L) = TC(N)%Conc(I,J,L) + DTC / AD(I,J,L)

       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    AD   => NULL()
    TC   => NULL()

  END SUBROUTINE TurbDay
!EOC
END MODULE Pbl_Mix_Mod
