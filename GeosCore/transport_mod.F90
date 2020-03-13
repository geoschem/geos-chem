!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transport_mod.F90
!
! !DESCRIPTION: Module TRANSPORT\_MOD is used to call the proper version of
!  the TPCORE advection scheme for different met field datasets and their
!  nested or global grids.
!\\
!\\
! !INTERFACE:
!
MODULE TRANSPORT_MOD
!
! !USES:
!
  USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)
  USE PRESSURE_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CLEANUP_TRANSPORT
  PUBLIC  :: DO_TRANSPORT
  PUBLIC  :: INIT_TRANSPORT
  PUBLIC  :: INIT_WINDOW_TRANSPORT
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DO_GLOBAL_ADV
  PRIVATE :: DO_WINDOW_TRANSPORT
!
! !REVISION HISTORY:
!  10 Mar 2003 - Y. Wang, R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  !=================================================================
  ! MODULE VARIABLES:
  !
  ! (1 ) Ap     (REAL(fp) ) : Vertical coordinate array for TPCORE
  ! (2 ) A_M2   (REAL(fp) ) : Grid box surface areas [m2]
  ! (3 ) Bp     (REAL(fp) ) : Vertical coordinate array for TPCORE
  ! (7 ) JLAST  (INTEGER)   : For fvDAS TPCORE
  ! (8 ) MG     (INTEGER)   : For fvDAS TPCORE
  ! (9 ) NG     (INTEGER)   : For fvDAS TPCORE
  ! (10) N_ADJ  (INTEGER)   : For fvDAS TPCORE
  !=================================================================
  INTEGER                       :: JFIRST
  INTEGER                       :: JLAST, NG,   MG,   N_ADJ
  REAL(fp), ALLOCATABLE         :: Ap(:)
  REAL(fp), ALLOCATABLE         :: Bp(:)
  REAL(fp), ALLOCATABLE, TARGET :: A_M2(:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_transport
!
! !DESCRIPTION: Subroutine DO\_TRANSPORT is the driver routine for the proper
!  TPCORE program for GEOS-3, GEOS-4/GEOS-5, or window simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_TRANSPORT( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_TS_DYN
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  10 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    REAL(fp)           :: DT_Dyn

    !=================================================================
    ! DO_TRANSPORT begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Transport  (in module GeosCore/transport_mod.F90)'

    !----------------------------------------------------------
    ! Transport (advection) budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetTransport ) THEN
       ! Get initial column masses
       CALL Compute_Column_Mass( Input_Opt, &
                                 State_Chm, State_Grid, State_Met, &
                                 State_Chm%Map_Advect, &
                                 State_Diag%Archive_BudgetTransportFull, &
                                 State_Diag%Archive_BudgetTransportTrop, &
                                 State_Diag%Archive_BudgetTransportPBL,  &
                                 State_Diag%BudgetMass1, &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Transport budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! First-time initialization
    IF ( FIRST ) THEN

       IF ( State_Grid%NestedGrid ) THEN

          ! All nested grid simulations
          CALL INIT_WINDOW_TRANSPORT( Input_Opt, State_Diag, State_Grid, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Init_Window_Transport"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ELSE

          ! All global simulations
          CALL INIT_TRANSPORT( Input_Opt, State_Diag, State_Grid, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Init_Transport"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDIF

       FIRST = .FALSE.

    ENDIF

    !=================================================================
    ! Choose the proper version of TPCORE for the nested-grid window
    ! region (usually 1x1 grids) or for the entire globe
    !=================================================================
    IF ( State_Grid%NestedGrid ) THEN

       ! Nested-grid simulation with GEOS-FP/MERRA2 met
       CALL DO_WINDOW_TRANSPORT( Input_Opt,  State_Chm, State_Diag, &
                                 State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Window_Transport"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !=================================================================
    ! Choose the proper version of TPCORE for global simulations
    !=================================================================
    ELSE

       ! Call TPCORE w/ proper settings for the GEOS-FP/MERRA2 met fields
       CALL DO_GLOBAL_ADV( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Global_Adv"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Dynamic timestep [s]
    DT_Dyn = GET_TS_DYN()

    !----------------------------------------------------------
    ! Transport (advection) budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetTransport ) THEN
       ! Get final column masses and compute diagnostics
       CALL Compute_Column_Mass( Input_Opt, &
                                 State_Chm, State_Grid, State_Met, &
                                 State_Chm%Map_Advect, &
                                 State_Diag%Archive_BudgetTransportFull, &
                                 State_Diag%Archive_BudgetTransportTrop, &
                                 State_Diag%Archive_BudgetTransportPBL, &
                                 State_Diag%BudgetMass2, &
                                 RC )
       CALL Compute_Budget_Diagnostics( State_Grid, &
                                 State_Chm%Map_Advect, &
                                 DT_Dyn, &
                                 State_Diag%Archive_BudgetTransportFull, &
                                 State_Diag%Archive_BudgetTransportTrop, &
                                 State_Diag%Archive_BudgetTransportPBL, &
                                 State_Diag%BudgetTransportFull, &
                                 State_Diag%BudgetTransportTrop, &
                                 State_Diag%BudgetTransportPBL, &
                                 State_Diag%BudgetMass1, &
                                 State_Diag%BudgetMass2, &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Transport budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE DO_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_global_adv
!
! !DESCRIPTION: Subroutine DO\_GLOBAL\_ADV is the driver routine
!  for TPCORE with the GMAO GEOS-FP or MERRA-2 met fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_GLOBAL_ADV( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : AIRQNT
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE PhysConstants            ! Physical constants
    USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
    USE TIME_MOD,           ONLY : GET_TS_DYN
    USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  As of July 2016, we assume that all of the advected species are listed
!  first in the species database.  This is the easiest way to pass a slab
!  to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)
!                                                                             .
!  Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
!  are incremented upside-down (level 1 = top of the atmosphere).
!  The levels order is reversed only when written out to diagnostic output.
!
! !REVISION HISTORY:
!  30 Oct 2007 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: LFILL
    LOGICAL            :: prtDebug
    INTEGER            :: IORD, JORD, KORD
    INTEGER            :: I, J, L, L2, N, N_DYN, NA, nAdvect
    INTEGER            :: ND24x, ND25x, ND26x
    REAL(fp)           :: D_DYN

    ! Arrays
    REAL(fp)           :: P_TP1 (State_Grid%NX,State_Grid%NY)
    REAL(fp)           :: P_TP2 (State_Grid%NX,State_Grid%NY)
    REAL(fp)           :: P_TEMP(State_Grid%NX,State_Grid%NY)
    REAL(fp),  TARGET  :: XMASS (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp),  TARGET  :: YMASS (State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    REAL(fp),  POINTER :: p_UWND (:,:,:  )
    REAL(fp),  POINTER :: p_VWND (:,:,:  )
    REAL(fp),  POINTER :: p_XMASS(:,:,:  )
    REAL(fp),  POINTER :: p_YMASS(:,:,:  )
    REAL(fp),  POINTER :: p_MFLEW(:,:,:,:)
    REAL(fp),  POINTER :: p_MFLNS(:,:,:,:)
    REAL(fp),  POINTER :: p_MFLUP(:,:,:,:)
    REAL(fp),  POINTER :: p_Spc  (:,:,:,:)

    !=================================================================
    ! DO_GLOBAL_ADV begins here!
    !=================================================================

    ! Assume success
    RC          =  GC_SUCCESS

    ! Initialize
    LFILL       =  Input_Opt%LFILL
    prtDebug    =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    IORD        =  Input_Opt%TPCORE_IORD
    JORD        =  Input_Opt%TPCORE_JORD
    KORD        =  Input_Opt%TPCORE_KORD
    nAdvect     =  State_Chm%nAdvect
    p_MFLEW     => NULL()
    p_MFLNS     => NULL()
    p_MFLUP     => NULL()

    ! Dynamic timestep [s]
    N_DYN       =  GET_TS_DYN()
    D_DYN       =  DBLE( N_DYN )

    ! Define shadow variables for ND24, ND25, ND26
    ND24x       = 0
    ND25x       = 0
    ND26x       = 0

    !=================================================================
    ! Prepare variables for calls to PJC pressure-fixer and TPCORE
    !
    ! For hybrid grids, the pressure at the
    ! bottom edge of grid box (I,J,L) is given by:
    !
    !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    !
    ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
    ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
    !=================================================================

    !!### DEBUG: Print a few global species sums
    !IF ( prtDebug ) THEN
    !   CALL Print_Global_Species_Kg( 20, 20, 1, 'O3',       &
    !                                 Input_Opt, State_Chm,  &
    !                                 State_Grid, State_Met, &
    !                                 "do_global_adv: pre-advection", &
    !                                 RC )
    !ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
       P_TP1(I,J) = GET_PEDGE_DRY(I,J,1)

       ! Set true dry sfc pressure at end of dynamic timestep [hPa]
       P_TP2(I,J) = State_Met%PSC2_DRY(I,J)

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=================================================================
    ! Call the PJC/LLNL pressure fixer to get the adjusted air
    ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
    ! TPCORE_FVDAS in order to ensure mass conservation.
    !=================================================================

    ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
    CALL DO_PJC_PFIX( State_Grid,  D_DYN, P_TP1, P_TP2, &
                      State_Met%U, State_Met%V, XMASS, YMASS )

    !=================================================================
    ! Call TPCORE_FVDAS to perform the advection
    !=================================================================

    ! Flip array indices in the vertical using pointer storage
    p_UWND    => State_Met%U      (:,:,State_Grid%NZ:1:-1)
    p_VWND    => State_Met%V      (:,:,State_Grid%NZ:1:-1)
    p_XMASS   => XMASS            (:,:,State_Grid%NZ:1:-1)
    p_YMASS   => YMASS            (:,:,State_Grid%NZ:1:-1)

    ! NOTE: For now, so as to avoid having to rewrite the internals
    ! of the TPCORE routines, just point to 1:nAdvect entries of
    ! State_Chm%Species.  This is OK for now, as of July 2016, all of
    ! the advected species are listed first.  This may change in the
    ! future, but we'll worry about that later.  The units of p_SPC
    ! will be converted to [kg/kg moist air] below. (bmy, 7/13/16)
    p_Spc     => State_Chm%Species(:,:,State_Grid%NZ:1:-1,1:nAdvect)

    ! Do the advection
    CALL TPCORE_FVDAS( D_DYN,     Re, &
                       State_Grid%NX, &
                       State_Grid%NY, &
                       State_Grid%NZ, &
                       JFIRST,    JLAST,    NG, &
                       MG,        nAdvect,  Ap,       Bp,    &
                       p_UWND,    p_VWND,   P_TP1,    P_TP2, &
                       P_TEMP,    p_Spc,    IORD,     JORD,  &
                       KORD,      N_ADJ,    p_XMASS,  p_YMASS, &
                       LFILL,     &
                       A_M2,      ND24x,    ND25x,    ND26x, &
                       State_Diag                            )

    ! Free pointer memory
    p_UWND  => NULL()
    p_VWND  => NULL()
    p_XMASS => NULL()
    p_YMASS => NULL()
    p_Spc   => NULL()
    p_MFLEW => NULL()
    p_MFLNS => NULL()
    p_MFLUP => NULL()

    !=================================================================
    ! Reset surface pressure and ensure mass conservation
    !=================================================================

    ! Update dry and wet floating pressures to the most recently
    ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
    ! (ewl, 7/6/16)
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC)

    ! Update State_Met air quantities with new pressures.
    ! Do not update tracer mixing ratio because after advection
    ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
    CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, &
                 RC, update_mixing_ratio=.FALSE. )

    !!### DEBUG: Print a few global species sums
    !IF ( prtDebug ) THEN
    !   CALL Print_Global_Species_Kg( 20, 20, 1, 'O3',       &
    !                                 Input_Opt, State_Chm,  &
    !                                 State_Grid, State_Met, &
    !                                 "do_global_adv: post-airqnt", &
    !                                 RC )
    !   CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' )
    !ENDIF

  END SUBROUTINE DO_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_window_transport
!
! !DESCRIPTION: Subroutine DO\_WINDOW\_TRANSPORT is the driver program
!  for the proper TPCORE program for the GEOS-FP/MERRA2 nested-grid
!  simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_WINDOW_TRANSPORT( Input_Opt,  State_Chm, State_Diag, &
                                  State_Grid, State_Met, RC )
!
! !USES:
!
    USE Calc_Met_Mod,         ONLY : AIRQNT
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE PhysConstants              ! Physical constants
    USE PJC_PFIX_WINDOW_MOD,  ONLY : DO_PJC_PFIX_WINDOW
    USE TIME_MOD,             ONLY : GET_TS_DYN
    USE TPCORE_WINDOW_MOD,    ONLY : TPCORE_WINDOW
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  As of July 2016, we assume that all of the advected species are listed
!  first in the species database.  This is the easiest way to pass a slab
!  to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)
!                                                                             .
!  Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
!  are incremented upside-down (level 1 = top of the atmosphere).
!  The levels order is reversed only when written out to diagnostic output.
!
! !REVISION HISTORY:
!  10 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: LFILL
    LOGICAL            :: prtDebug
    INTEGER            :: IORD,  JORD,  KORD
    INTEGER            :: I0,    J0,    NA,    nAdvect, N_DYN
    INTEGER            :: IM_W1, JM_W1, I0_W1, J0_W1,   BUFF_SIZE
    INTEGER            :: I,     J,     L,     L2,      N
    INTEGER            :: ND24x, ND25x, ND26x
    REAL(fp)           :: D_DYN

    ! Arrays
    REAL(fp),  TARGET  :: P_TP1 (State_Grid%NX,State_Grid%NY)
    REAL(fp),  TARGET  :: P_TP2 (State_Grid%NX,State_Grid%NY)
    REAL(fp),  TARGET  :: P_TEMP(State_Grid%NX,State_Grid%NY)
    REAL(fp),  TARGET  :: XMASS (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp),  TARGET  :: YMASS (State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    REAL(fp),  POINTER :: p_A_M2  (  :    )
    REAL(fp),  POINTER :: p_P_TP1 (:,:    )
    REAL(fp),  POINTER :: p_P_TP2 (:,:    )
    REAL(fp),  POINTER :: p_P_TEMP(:,:    )
    REAL(fp),  POINTER :: p_UWND  (:,:,:  )
    REAL(fp),  POINTER :: p_VWND  (:,:,:  )
    REAL(fp),  POINTER :: p_XMASS (:,:,:  )
    REAL(fp),  POINTER :: p_YMASS (:,:,:  )
    REAL(fp),  POINTER :: p_MFLEW (:,:,:,:)
    REAL(fp),  POINTER :: p_MFLNS (:,:,:,:)
    REAL(fp),  POINTER :: p_MFLUP (:,:,:,:)
    REAL(fp),  POINTER :: p_Spc   (:,:,:,:)

    !=================================================================
    ! DO_FVDAS_WINDOW_TRANSPORT begins here!
    !=================================================================

    ! Assume success
    RC          =  GC_SUCCESS

    ! Copy values from Input_Opt
    LFILL       =  Input_Opt%LFILL
    IORD        =  Input_Opt%TPCORE_IORD
    JORD        =  Input_Opt%TPCORE_JORD
    KORD        =  Input_Opt%TPCORE_KORD
    nAdvect     =  State_Chm%nAdvect
    prtDebug    =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Initialize pointers
    p_A_M2      => NULL()
    p_P_TP1     => NULL()
    p_P_TP2     => NULL()
    p_P_TEMP    => NULL()
    p_UWND      => NULL()
    p_VWND      => NULL()
    p_XMASS     => NULL()
    p_YMASS     => NULL()
    p_MFLEW     => NULL()
    p_MFLNS     => NULL()
    p_MFLUP     => NULL()
    p_Spc       => NULL()

    ! Get nested-grid lon/lat offsets [# boxes]
    I0          =  State_Grid%XMinOffset
    J0          =  State_Grid%YMinOffset

    ! Dynamic timestep [s]
    N_DYN       =  GET_TS_DYN()
    D_DYN       =  DBLE( N_DYN )

    ! (lzh, 09/01/2014)
    BUFF_SIZE   =  2
    IM_W1       =  ( State_Grid%NX - State_Grid%WestBuffer - &
                     State_Grid%EastBuffer  ) + 2 * BUFF_SIZE
    JM_W1       =  ( State_Grid%NY - State_Grid%SouthBuffer - &
                     State_Grid%NorthBuffer ) + 2 * BUFF_SIZE
    I0_W1       =  State_Grid%WestBuffer  - BUFF_SIZE
    J0_W1       =  State_Grid%SouthBuffer - BUFF_SIZE

    ! Define shadow variables for ND24, ND25, ND26
    ND24x       = 0
    ND25x       = 0
    ND26x       = 0

    !=================================================================
    ! Prepare variables for calls to PJC pressure-fixer and TPCORE
    !
    ! For hybrid grids, the pressure at the
    ! bottom edge of grid box (I,J,L) is given by:
    !
    !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    !
    ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
    ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
    !=================================================================

    !IF ( prtDebug ) THEN
    !   CALL Print_Global_Species_Kg( 20, 20, 1, 'SPC_O3',   &
    !                                 Input_Opt,  State_Chm, &
    !                                 State_Grid, State_Met, &
    !                                 "do_window_transport: pre-advection", &
    !                                 RC )
    !ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
       P_TP1(I,J) = GET_PEDGE_DRY(I,J,1)

       ! Set true dry sfc pressure at end of dynamic timestep [hPa]
       P_TP2(I,J) = State_Met%PSC2_DRY(I,J)

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=================================================================
    ! Call the PJC/LLNL pressure fixer to get the adjusted air
    ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
    ! TPCORE_FVDAS in order to ensure mass conservation.
    !=================================================================
    XMASS = 0e+0_fp !(dan)
    YMASS = 0e+0_fp
    ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
    CALL DO_PJC_PFIX_WINDOW( State_Grid,  D_DYN, &
                             P_TP1,       P_TP2, &
                             State_Met%U, State_Met%V, &
                             XMASS,       YMASS )

    IF ( prtDebug ) CALL DEBUG_MSG( '### FVDAS_WINDOW: a PJC_PFIX_WINDOW')

    ! Flip array indices in the vertical using pointer storage

    ! Exclude the buffer zone (lzh, 4/1/2015)
    p_UWND  => State_Met%U( I0_W1+1 : I0_W1+IM_W1, &
                            J0_W1+1 : J0_W1+JM_W1, &
                            State_Grid%NZ:1:-1 )

    p_VWND  => State_Met%V( I0_W1+1 : I0_W1+IM_W1, &
                            J0_W1+1 : J0_W1+JM_W1, &
                            State_Grid%NZ:1:-1  )

    ! NOTE: For now, so as to avoid having to rewrite the internals
    ! of the TPCORE routines, just point to 1:nAdvect entries of
    ! State_Chm%Species.  This is OK for now, as of July 2016, all of
    ! the advected species are listed first.  This may change in the
    ! future, but we'll worry about that later. (bmy, 7/13/16)
    p_Spc => State_Chm%Species( I0_W1+1 : I0_W1+IM_W1, &
                                J0_W1+1 : J0_W1+JM_W1, &
                                State_Grid%NZ:1:-1,    &
                                1:nAdvect )

    p_XMASS  => XMASS( I0_W1+1 : I0_W1+IM_W1, &
                       J0_W1+1 : J0_W1+JM_W1, &
                       State_Grid%NZ:1:-1 )

    p_YMASS  => YMASS( I0_W1+1 : I0_W1+IM_W1, &
                       J0_W1+1 : J0_W1+JM_W1, &
                       State_Grid%NZ:1:-1 )

    p_P_TP1  => P_TP1( I0_W1+1 : I0_W1+IM_W1, &
                       J0_W1+1 : J0_W1+JM_W1 )

    p_P_TP2  => P_TP2( I0_W1+1 : I0_W1+IM_W1, &
                       J0_W1+1 : J0_W1+JM_W1)

    p_P_TEMP => P_TEMP( I0_W1+1 : I0_W1+IM_W1, &
                        J0_W1+1 : J0_W1+JM_W1 )

    p_A_M2   => A_M2( J0_W1+1 : J0_W1+JM_W1 )

    ! Do the advection
    CALL TPCORE_WINDOW(D_DYN,   Re,       IM_W1,   JM_W1,   State_Grid%NZ, &
                       JFIRST,  JLAST,    NG,      MG,      nAdvect,       &
                       Ap,      Bp,       p_UWND,  p_VWND,  p_P_TP1,       &
                       p_P_TP2, p_P_TEMP, p_Spc,   IORD,    JORD,          &
                       KORD,    N_ADJ,    p_XMASS, p_YMASS,                &
                       p_A_M2,  ND24x,    ND25x,   ND26x )

    ! Free pointer memory
    p_UWND   => NULL()
    p_VWND   => NULL()
    p_Spc    => NULL()
    p_XMASS  => NULL()
    p_YMASS  => NULL()
    p_P_TP1  => NULL()
    p_P_TP2  => NULL()
    p_P_TEMP => NULL()
    p_A_M2   => NULL()
    p_MFLEW  => NULL()
    p_MFLNS  => NULL()
    p_MFLUP  => NULL()

    !=================================================================
    ! Reset surface pressure and ensure mass conservation
    !=================================================================

    ! Update dry and wet floating pressures to the most recently
    ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
    ! (ewl, 7/6/16)
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC)

    ! Update State_Met air quantities with new pressures.
    ! Do not update tracer mixing ratio because after advection
    ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
    CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
                 Update_Mixing_Ratio=.FALSE. )

    !!### Debug
    !IF ( prtDebug ) THEN
    !   CALL Print_Global_Species_Kg( 20, 20, 1, 'SPC_O3',   &
    !                                 Input_Opt, State_Chm,  &
    !                                 State_Grid, State_Met, &
    !                                 "do_window_transport: post-airqnt", &
    !                                 RC )
    !   CALL DEBUG_MSG( '### NESTED_ADV: a TPCORE' )
    !ENDIF

  END SUBROUTINE DO_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_transport
!
! !DESCRIPTION: Subroutine INIT\_TRANSPORT initializes all module variables
!  and arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TRANSPORT( Input_Opt, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,        ONLY : ALLOC_ERR
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE PhysConstants          ! Re
    USE TIME_MOD,         ONLY : GET_TS_DYN
    USE TPCORE_FVDAS_MOD, ONLY : INIT_TPCORE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  10 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: LTRAN
    INTEGER            :: J, K, L, N_DYN
    REAL(fp)           :: YMID_R(State_Grid%NY)
    REAL(fp)           :: REAL_N_DYN

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! Initialize
    !=================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Transport (in module GeosCore/transport_mod.F90)'

    !=================================================================
    ! Allocate arrays for TPCORE vertical coordinates
    !
    ! For fvDAS TPCORE with for GEOS-FP or MERRA-2 met fields:
    !
    !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    !
    ! Also here Ap, Bp will be flipped since both TPCORE versions
    ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    !=================================================================
    ALLOCATE( Ap( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'transport_mod.F:Ap', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ALLOCATE( Bp( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'transport_mod.F:Bp', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Flip Ap and Bp for TPCORE
    DO L = 1, State_Grid%NZ+1

       ! As L runs from the surface up,
       ! K runs from the top down
       K = ( State_Grid%NZ + 1 ) - L + 1

       Ap(L) = GET_AP(K)          ! Ap(L) is in [hPa]
       Bp(L) = GET_BP(K)
    ENDDO

    !=================================================================
    ! Allocate arrays for surface area and layer thickness
    !=================================================================
    ALLOCATE( A_M2( State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'transport_mod.F:A_m2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Surface area [m2]
    DO J = 1, State_Grid%NY
       A_M2(J) = State_Grid%Area_M2(1,J)
    ENDDO

    !=================================================================
    ! Additional setup for the fvDAS version of TPCORE
    !=================================================================

    ! Initialize
    N_DYN = GET_TS_DYN()
    N_ADJ = 0
    NG    = 0
    MG    = 0

    ! YMID_R is latitude of grid box center [radian]
    DO J = 1,State_Grid%NY
       YMID_R(J) = State_Grid%YMid_R(1,J)
    ENDDO

    REAL_N_DYN = N_DYN

    ! Call INIT routine from "tpcore_fvdas_mod.f"
    CALL INIT_TPCORE( State_Grid%NX, State_Grid%NY,  State_Grid%NZ, &
                      JFIRST, JLAST, NG, MG,         REAL_N_DYN,    &
                      Re,    YMID_R, State_Diag, RC          )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Init_Tpcore"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE INIT_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_window_transport
!
! !DESCRIPTION: Subroutine INIT\_WINDOW\_TRANSPORT initializes all
!  module variables and arrays for the GEOS-FP/MERRA2 nested grid
!  simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_WINDOW_TRANSPORT( Input_Opt, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants            ! Re
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_TS_DYN
    USE TPCORE_FVDAS_MOD,   ONLY : INIT_TPCORE
    USE TPCORE_WINDOW_MOD,  ONLY : INIT_WINDOW
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  06 Jun 2008 - D. Chen & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL  :: LTRAN
    INTEGER  :: BUFF_SIZE
    INTEGER  :: J,     K,     L,     N_DYN
    INTEGER  :: IM_W1, JM_W1, I0_W1, J0_W1
    REAL(fp) :: REAL_N_DYN

    ! Arrays
    REAL(fp) :: YMID_R_W(0:State_Grid%NY+1)

    !=================================================================
    ! Initialize
    !=================================================================

    ! Assume success
    RC        =  GC_SUCCESS

    ! Copy values from Input_Opt
    LTRAN     = Input_Opt%LTRAN

    ! Cast N_DYN to flexible precision
#if defined( USE_REAL8 )
    REAL_N_DYN = DBLE( N_DYN )
#else
    REAL_N_DYN = FLOAT( N_DYN )
#endif

    !=================================================================
    ! Allocate arrays for TPCORE vertical coordinates
    ! GEOS-FP/MERRA2 nested grid simulation only!!!
    !
    ! For fvDAS TPCORE with for GEOS-FP/MERRA2 met fields:
    !
    !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    !
    ! Also here Ap, Bp will be flipped since both TPCORE versions
    ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    !=================================================================
    ALLOCATE( Ap( State_Grid%NZ+1 ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )

    ALLOCATE( Bp( State_Grid%NZ+1 ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

    ! Flip Ap and Bp for TPCORE
    DO L = 1, State_Grid%NZ+1

       ! As L runs from the surface up,
       ! K runs from the top down
       K = ( State_Grid%NZ + 1 ) - L + 1

       Ap(L) = GET_AP(K)
       Bp(L) = GET_BP(K)
    ENDDO

    !=================================================================
    ! Allocate arrays for surface area and layer thickness
    !=================================================================
    ALLOCATE( A_M2( State_Grid%NY ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

    ! Surface area [m2]
    DO J = 1, State_Grid%NY
       A_M2(J) = State_Grid%Area_M2(1,J)
    ENDDO

    !=================================================================
    ! Additional setup for the fvDAS version of TPCORE
    !=================================================================

    ! Initialize
    N_DYN = GET_TS_DYN()
    N_ADJ = 0
    NG    = 0
    MG    = 0

    ! (lzh, 4/1/2015)
    BUFF_SIZE = 2
    IM_W1       =  ( State_Grid%NX - State_Grid%WestBuffer - &
                     State_Grid%EastBuffer  ) + 2 * BUFF_SIZE
    JM_W1       =  ( State_Grid%NY - State_Grid%SouthBuffer - &
                     State_Grid%NorthBuffer ) + 2 * BUFF_SIZE
    I0_W1     = State_Grid%WestBuffer  - BUFF_SIZE
    J0_W1     = State_Grid%SouthBuffer - BUFF_SIZE

    ! YMID_R is latitude of grid box center [radians]
    DO J = 1, State_Grid%NY
       YMID_R_W(J) = State_Grid%YMid_R(1,J)
    ENDDO

    ! Compute YMID_R_W at southern edge of nested region
    J = 0
    YMID_R_W(J) = State_Grid%YMid_R(1,J+1) - (State_Grid%DY * PI_180)

    ! Compute YMID_R_W at northern edge of nested region
    J = State_Grid%NY+1
    YMID_R_W(J) = State_Grid%YMid_R(1,J-1) + (State_Grid%DY * PI_180)

    REAL_N_DYN = N_DYN

    ! Call INIT routine from "tpcore_window_mod.F90"
    CALL INIT_WINDOW( State_Grid,    &
                      IM_W1,         &
                      JM_W1,         &
                      State_Grid%NZ, &
                      JFIRST,        &
                      JLAST,         &
                      NG,            &
                      MG,            &
                      REAL_N_DYN,    &
                      Re,            &
                      YMID_R_W( J0_W1:(J0_W1+JM_W1+1) ) )

  END SUBROUTINE INIT_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_transport
!
! !DESCRIPTION: Subroutine CLEANUP\_TRANSPORT deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_TRANSPORT
!
! !REVISION HISTORY:
!  10 Mar 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
    IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
    IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )

  END SUBROUTINE CLEANUP_TRANSPORT
!EOC
END MODULE TRANSPORT_MOD
