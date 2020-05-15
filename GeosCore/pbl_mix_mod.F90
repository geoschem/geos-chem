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
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE Time_Mod,        ONLY : Get_Ts_Conv
    USE Time_Mod,        ONLY : Get_Ts_Dyn
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units
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

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    REAL(fp)           :: DT_Dyn

    !=======================================================================
    ! Do_Full_Pbl_Mixing begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_PBL_MIX (in module GeosCore/pbl_mix_mod.F90)'

    !========================================================================
    ! Full PBL mixing budget diagnostics - Part 1 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get initial column masses
       CALL Compute_Column_Mass( Input_Opt,                                  &
                                 State_Chm,                                  &
                                 State_Grid,                                 &
                                 State_Met,                                  &
                                 State_Chm%Map_Advect,                       &
                                 State_Diag%Archive_BudgetMixingFull,        &
                                 State_Diag%Archive_BudgetMixingTrop,        &
                                 State_Diag%Archive_BudgetMixingPBL,         &
                                 State_Diag%BudgetMass1,                     &
                                 RC                                         )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing budget diagnostics error 1 (full PBL mixing)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Unit conversion #1
    !========================================================================

    ! Convert species to v/v dry
    CALL Convert_Spc_Units( Input_Opt,         State_Chm, State_Grid,        &
                            State_Met,        'v/v dry',  RC,                &
                            OrigUnit=OrigUnit                               )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Do full PBL mixing
    !========================================================================

    ! Do complete mixing of tracers in the PBL
    CALL TurbDay( Input_Opt,  State_Chm, State_Diag,                         &
                  State_Grid, State_Met, RC                                 )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "TURBDAY"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Unit conversion #2
    !========================================================================

    ! Convert species back to original units
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (from v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Full PBL mixing budget diagnostics - Part 1 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetMixing ) THEN

       ! Get dynamics timestep [s]
       DT_Dyn = Get_Ts_Dyn()

       ! Get final column masses and compute diagnostics
       CALL Compute_Column_Mass( Input_Opt,                                  &
                                 State_Chm,                                  &
                                 State_Grid,                                 &
                                 State_Met,                                  &
                                 State_Chm%Map_Advect,                       &
                                 State_Diag%Archive_BudgetMixingFull,        &
                                 State_Diag%Archive_BudgetMixingTrop,        &
                                 State_Diag%Archive_BudgetMixingPBL,         &
                                 State_Diag%BudgetMass2,                     &
                                 RC                                         )


       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Mixing budget diagnostics error 2 (full PBL mixing)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Archive budget diagnostic
       CALL Compute_Budget_Diagnostics( State_Grid,                          &
                                        State_Chm%Map_Advect,                &
                                        DT_Dyn,                              &
                                        State_Diag%Archive_BudgetMixingFull, &
                                        State_Diag%Archive_BudgetMixingTrop, &
                                        State_Diag%Archive_BudgetMixingPBL,  &
                                        State_Diag%BudgetMixingFull,         &
                                        State_Diag%BudgetMixingTrop,         &
                                        State_Diag%BudgetMixingPBL,          &
                                        State_Diag%BudgetMass1,              &
                                        State_Diag%BudgetMass2,              &
                                        RC                                  )

       ! Trap potential errors 
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Compute_Budget_Diagnostics!"'
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
  SUBROUTINE Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : Scale_Height
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
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
    LOGICAL  :: Bad_Sum
    INTEGER  :: I,      J,      L,    LTOP
    REAL(fp) :: BLTOP,  BLTHIK, DELP

    ! Arrays
    REAL(fp) :: P(0:State_Grid%NZ)

    !=================================================================
    ! COMPUTE_PBL_HEIGHT begins here!
    !=================================================================

    ! Initialize
    RC              = GC_SUCCESS
    Bad_Sum         = .FALSE.
    State_Met%InPbl = .FALSE.

    !$OMP PARALLEL DO                                      &
    !$OMP DEFAULT( SHARED                                ) &
    !$OMP PRIVATE( I, J, L, P, BLTOP, BLTHIK, LTOP, DELP )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !----------------------------------------------
       ! Define pressure edges:
       ! P(L-1) = P at bottom edge of box (I,J,L)
       ! P(L  ) = P at top    edge of box (I,J,L)
       !----------------------------------------------

       ! Pressure at level edges [hPa]
       DO L = 0, State_Grid%NZ
          P(L) = State_Met%PEDGE(I,J,L+1)
       ENDDO

       !----------------------------------------------
       ! Find PBL top and thickness [hPa]
       !----------------------------------------------

       ! BLTOP = pressure at PBL top [hPa]
       ! Use barometric law since PBL is in [m]
       BLTOP  = P(0) * EXP( -State_Met%PBLH(I,J) / SCALE_HEIGHT )

       ! BLTHIK is PBL thickness [hPa]
       BLTHIK = P(0) - BLTOP

       !----------------------------------------------
       ! Find model level where BLTOP occurs
       !----------------------------------------------

       ! Initialize
       LTOP = 0

       ! Loop over levels
       DO L = 1, State_Grid%NZ

          ! Exit when we get to the PBL top level
          IF ( BLTOP > P(L) ) THEN
             LTOP = L
             EXIT
          ELSE
             State_Met%InPbl(I,J,L) = .TRUE.
          ENDIF

       ENDDO

       !----------------------------------------------
       ! Define various related quantities
       !----------------------------------------------

       ! IMIX(I,J)   is the level where the PBL top occurs at (I,J)
       ! IMIX(I,J)-1 is the number of whole levels below the PBL top
       State_Met%IMIX(I,J) = LTOP

       ! Fraction of the IMIXth level underneath the PBL top
       State_Met%FPBL(I,J) = 1.0_fp                                          &
                           - ( BLTOP - P(LTOP) ) / ( P(LTOP-1) - P(LTOP) )

       ! PBL top [model layers]
       State_Met%PBL_TOP_L(I,J) = FLOAT( State_Met%IMIX(I,J) - 1 )           &
                                + State_Met%FPBL(I,J)

       ! PBL top [hPa]
       State_Met%PBL_TOP_hPa(I,J) = BLTOP

       ! Zero PBL top [m] -- compute below
       State_Met%PBL_TOP_m(I,J) = 0.0_fp

       ! PBL thickness [hPa]
       State_Met%PBL_THICK(I,J) = BLTHIK

       !==============================================================
       ! Loop up to edge of chemically-active grid
       !==============================================================
       DO L = 1, State_Grid%MaxChemLev

          ! Thickness of grid box (I,J,L) [hPa]
          DELP = P(L-1) - P(L)

          IF ( L < State_Met%IMIX(I,J) ) THEN

             !--------------------------------------------
             ! (I,J,L) lies completely below the PBL top
             !--------------------------------------------

             ! Fraction of grid box (I,J,L) w/in the PBL
             State_Met%F_OF_PBL(I,J,L) = DELP / BLTHIK

             ! Fraction of grid box (I,J,L) underneath PBL top
             State_Met%F_UNDER_PBLTOP(I,J,L) = 1.0_fp

             ! PBL height [m]
             State_Met%PBL_TOP_m(I,J) = State_Met%PBL_TOP_m(I,J) + &
                                        State_Met%BXHEIGHT(I,J,L)

          ELSE IF ( L == State_Met%IMIX(I,J) ) THEN

             !--------------------------------------------
             ! (I,J,L) straddles the PBL top
             !--------------------------------------------

             ! Fraction of grid box (I,J,L) w/in the PBL
             State_Met%F_OF_PBL(I,J,L) = ( P(L-1) - BLTOP ) / BLTHIK

             ! Fraction of grid box (I,J,L) underneath PBL top
             State_Met%F_UNDER_PBLTOP(I,J,L) = State_Met%FPBL(I,J)

             ! PBL height [m]
             State_Met%PBL_TOP_m(I,J) = State_Met%PBL_TOP_m(I,J)             &
                                      + ( State_Met%BXHEIGHT(I,J,L)          &
                                      *   State_Met%FPBL(I,J)               )

          ELSE

             !--------------------------------------------
             ! (I,J,L) lies completely above the PBL top
             !--------------------------------------------

             ! Fraction of grid box (I,J,L) w/in the PBL
             State_Met%F_OF_PBL(I,J,L)       = 0.0_fp

             ! Fraction of grid box (I,J,L) underneath PBL top
             State_Met%F_UNDER_PBLTOP(I,J,L) = 0.0_fp

          ENDIF

          !### Debug
          !IF ( I==23 .and. J==34 .and. L < 6 ) THEN
          !   PRINT*, '###--------------------------------------'
          !   PRINT*, '### COMPUTE_PBL_HEIGHT'
          !   PRINT*, '### I, J, L     : ', I, J, L
          !   PRINT*, '### P(L-1)      : ', P(L-1)
          !   PRINT*, '### P(L)        : ', P(L)
          !   PRINT*, '### F_OF_PBL    : ', State_Met%F_OF_PBL(I,J,L)
          !   PRINT*, '### F_UNDER_TOP : ', &
          !            State_Met%F_UNDER_PBLTOP(I,J,L)
          !   PRINT*, '### IMIX        : ', IMIX(I,J)
          !   PRINT*, '### FPBL        : ', FPBL(I,J)
          !   PRINT*, '### PBL_TOP_hPa : ', State_Met%PBL_TOP_hPa(I,J)
          !   PRINT*, '### PBL_TOP_L   : ', State_Met%PBL_TOP_L(I,J)
          !   PRINT*, '### DELP        : ', DELP
          !   PRINT*, '### BLTHIK      : ', BLTHIK
          !   PRINT*, '### BLTOP       : ', BLTOP
          !   PRINT*, '### BXHEIGHT    : ', State_Met%BXHEIGHT(I,J,L)
          !   PRINT*, '### PBL_TOP_m   : ', State_Met%PBL_TOP_m(I,J)
          !   PRINT*, '### other way m : ', &
          !    P(0) * EXP( -State_Met%PBL_TOP_hPa(I,J) / SCALE_HEIGHT )
          !ENDIF

       ENDDO

       ! Error check
       IF ( ABS( SUM( State_Met%F_OF_PBL(I,J,:) ) - 1.0_fp) > 1.0e-3_fp) THEN
          !$OMP CRITICAL
          PRINT*, 'bad sum at: ', I, J
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
    State_Met%PBL_MAX_L = MAXVAL( State_Met%IMIX )

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
    REAL(fp)           :: A(State_Grid%NX,State_Grid%NY)
    REAL(fp)           :: DTC(State_Grid%NX,State_Grid%NY,                   &
                              State_Grid%NZ,State_Chm%nAdvect)
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    INTEGER,  POINTER  :: IMIX(:,:    )
    REAL(fp), POINTER  :: FPBL(:,:    )
    REAL(fp), POINTER  :: AD  (:,:,:  )
    REAL(fp), POINTER  :: TC  (:,:,:,:)

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
    AD      => State_Met%AD        ! Dry air mass
    IMIX    => State_Met%IMIX      ! Integer level where PBL top occurs
    FPBL    => State_Met%FPBL      ! Fractional level above IMIX to PBL top
    TC      => State_Chm%Species   ! Chemical species [v/v]

    ! Convection timestep [s]
    DTCONV = GET_TS_CONV()

    ! Loop over Lat/Long grid boxes (I,J)
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, NA, N, AA, CC, CC_AA )
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
             CC = CC + AD(I,J,L) * TC(I,J,L,N)
          ENDDO

          ! Then also sum mass from (I,J,L) which straddle the PBL top
          L     = IMIX(I,J)
          CC    = CC + AD(I,J,L) * TC(I,J,L,N) * FPBL(I,J)

          ! CC/AA is the mean mixing ratio of tracer at
          ! (I,J) from L=1 to L=LTOP
          CC_AA = CC / AA

          !========================================================
          ! TC(I,J,L,N) new  = TC(I,J,L,N) old +
          !                    ( DTC(I,J,L,N) / AD(I,J,L) )
          !
          ! where
          !
          ! DTC(I,J,L,N) = [ alpha * (mean MR below PBL) *
          !                  Airmass at (I,J,L) ] -
          !                [ alpha * TC(I,J,L,N) old     *
          !                  Airmass at (I,J,L) ]
          !
          ! DTC is thus the change in mass (kg) due to BL mixing,
          ! so DTC/AD is the change in (V/V) mixing ratio units.
          !========================================================

          ! For grid boxes (I,J,L) which lie below the PBL top
          DO L = 1, IMIX(I,J)-1
             DTC(I,J,L,N) = ( A(I,J) * CC_AA       * AD(I,J,L)  - &
                              A(I,J) * TC(I,J,L,N) * AD(I,J,L) )

             TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L)
          ENDDO

          ! For grid boxes (I,J,L) which straddle the PBL top
          L = IMIX(I,J)

          DTC(I,J,L,N)  = ( A(I,J) * FPBL(I,J)  * CC_AA       * AD(I,J,L) - &
                            A(I,J) * FPBL(I,J)  * TC(I,J,L,N) * AD(I,J,L) )

          TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L)

       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    AD   => NULL()
    IMIX => NULL()
    FPBL => NULL()
    TC   => NULL()

  END SUBROUTINE TurbDay
!EOC
END MODULE Pbl_Mix_Mod
