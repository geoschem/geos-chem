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
MODULE PBL_MIX_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CLEANUP_PBL_MIX
  PUBLIC  :: DO_PBL_MIX
  PUBLIC  :: INIT_PBL_MIX
  PUBLIC  :: COMPUTE_PBL_HEIGHT
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: TURBDAY
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Arrays
  INTEGER,   ALLOCATABLE :: IMIX(:,:)  ! # of levels under PBL top
  REAL(fp),  ALLOCATABLE :: FPBL(:,:)  ! frac # of levels under PBL top

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
  SUBROUTINE DO_PBL_MIX( DO_TURBDAY, Input_Opt, State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Diagnostics_Mod,    ONLY : Compute_Column_Mass
    USE Diagnostics_Mod,    ONLY : Compute_Budget_Diagnostics
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CONV
    USE Time_Mod,           ONLY : Get_Ts_Dyn
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: DO_TURBDAY  ! =T means call TURBDAY
                                                 !    for full PBL mixing
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
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: N, NA

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    REAL(fp)           :: DT_Dyn

    !=================================================================
    ! DO_PBL_MIX begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_PBL_MIX (in module GeosCore/pbl_mix_mod.F90)'

    !-----------------------------
    ! First-time initialization
    !-----------------------------
    IF ( FIRST ) THEN

       ! Allocate arrays etc.
       CALL INIT_PBL_MIX( Input_Opt, State_Grid, RC )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "INIT_PBL_MIX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! The following only needs be done if full PBL mixing is on...
    !=================================================================
    IF ( DO_TURBDAY ) THEN

       !-------------------------------------------------
       ! Full PBL mixing budget diagnostics - Part 1 of 2
       !-------------------------------------------------
       IF ( State_Diag%Archive_BudgetMixing ) THEN
          ! Get initial column masses
          CALL Compute_Column_Mass( Input_Opt,                           &
                                    State_Chm,                           &
                                    State_Grid,                          &
                                    State_Met,                           &
                                    State_Chm%Map_Advect,                &
                                    State_Diag%Archive_BudgetMixingFull, &
                                    State_Diag%Archive_BudgetMixingTrop, &
                                    State_Diag%Archive_BudgetMixingPBL,  &
                                    State_Diag%BudgetMass1,              &
                                    RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Mixing budget diagnostics error 1 (full PBL mixing)'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------
       ! Unit conversion #1
       !--------------------------

       ! Convert species to v/v dry
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                               'v/v dry', RC, OrigUnit=OrigUnit )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !--------------------------
       ! Do full PBL mixing
       !--------------------------

       ! Do complete mixing of tracers in the PBL
       CALL TURBDAY( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "TURBDAY"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !--------------------------
       ! Unit conversion #2
       !--------------------------

       ! Convert species back to original units
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                               OrigUnit,  RC )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "Convert_Spc_Units" (from v/v dry)!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !-------------------------------------------------
       ! Full PBL mixing budget diagnostics - Part 2 of 2
       !-------------------------------------------------
       IF ( State_Diag%Archive_BudgetMixing ) THEN

          ! Get dynamics timestep [s]
          DT_Dyn = Get_Ts_Dyn()

          ! Get final column masses and compute diagnostics
          CALL Compute_Column_Mass( Input_Opt,                           &
                                    State_Chm,                           &
                                    State_Grid,                          &
                                    State_Met,                           &
                                    State_Chm%Map_Advect,                &
                                    State_Diag%Archive_BudgetMixingFull, &
                                    State_Diag%Archive_BudgetMixingTrop, &
                                    State_Diag%Archive_BudgetMixingPBL,  &
                                    State_Diag%BudgetMass2,              &
                                    RC )
          CALL Compute_Budget_Diagnostics( State_Grid,                         &
                                           State_Chm%Map_Advect,               &
                                           DT_Dyn,                             &
                                           State_Diag%Archive_BudgetMixingFull,&
                                           State_Diag%Archive_BudgetMixingTrop,&
                                           State_Diag%Archive_BudgetMixingPBL, &
                                           State_Diag%BudgetMixingFull,        &
                                           State_Diag%BudgetMixingTrop,        &
                                           State_Diag%BudgetMixingPBL,         &
                                           State_Diag%BudgetMass1,             &
                                           State_Diag%BudgetMass2,             &
                                           RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Mixing budget diagnostics error 2 (full PBL mixing)'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

  END SUBROUTINE DO_PBL_MIX
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
  SUBROUTINE COMPUTE_PBL_HEIGHT( State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE PhysConstants            ! g0, Rd
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
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
    LOGICAL  :: Bad_Sum
    INTEGER  :: I,     J,      L,    LTOP
    REAL(fp) :: BLTOP, BLTHIK, DELP, Lower_Edge_Height
    REAL(fp) :: P(0:State_Grid%NZ)

    !=================================================================
    ! COMPUTE_PBL_HEIGHT begins here!
    !=================================================================

    ! Initialize
    RC              = GC_SUCCESS
    Bad_Sum         = .FALSE.
    State_Met%InPbl = .FALSE.
    State_Met%F_Under_PBLTop = 0d0
    State_Met%F_of_PBL       = 0d0

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, P, BLTOP, BLTHIK, LTOP, DELP ) &
    !$OMP Private( Lower_Edge_Height )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! PBL height above surface, m 
       State_Met%PBL_Top_m(I,J) = State_Met%PBLH(I,J)

       ! Height of lower edge above surface, m
       Lower_Edge_Height = 0d0

       ! Find PBL top level (L) and pressure (hPa)
       Do L=1, State_Grid%NZ

         If ( Lower_Edge_Height + State_Met%BXHEIGHT(I,J,L) > State_Met%PBLH(I,J) ) then

            ! PBL top is in this level
            LTOP = L

            ! Pressure at the PBL top altitude, hPa
            ! Use pressure lapse equation: p(PBLH) = p(z1) * exp( -(PBLH-z1) / Scale_Height )
            ! p(z1) = State_Met%PEDGE(I,J,L) = Pressure at the lower level edge
            ! PBLH - z1 = (State_Met%PBLH(I,J,L) - Lower_Edge_Height) = Height above the lower level edge
            ! Scale_Height = Rd * Tv / g0
            State_Met%PBL_Top_hPa(I,J) = State_Met%PEdge(I,J,L) * &
                  EXP( -(State_Met%PBLH(I,J) - Lower_Edge_Height) * g0 / ( Rd * State_Met%Tv(I,J,L) ) )

            ! Fraction of PBL mass in layer L, will be normalized below
            State_Met%F_of_PBL(I,J,L) = State_Met%PEdge(I,J,L) - State_Met%PBL_Top_hPa(I,J)
      
            ! Fraction of the grid cell mass under PBL top
            State_Met%F_Under_PBLTop(I,J,L) =   State_Met%F_of_PBL(I,J,L) / &
                                              ( State_Met%PEdge(I,J,L) - State_Met%PEdge(I,J,L+1) )

            ! Model level of PBL top (integer+fraction). The top is within level CEILING(PBL_Top_L) 
            State_Met%PBL_Top_L(I,J) = (LTOP-1) + State_Met%F_Under_PBLTop(I,J,L)

            ! PBL Thickness from surface to top, hPa
            State_Met%PBL_Thick(I,J) = State_Met%PEdge(I,J,1) - State_Met%PBL_Top_hPa(I,J)

            !! Exit Do loop after we found PBL top level
            Exit

         Else 

            ! Grid cell fully within PBL
            State_Met%inPBL(I,J,L) = .True.

            ! Fraction of the grid cell mass under PBL top
            State_Met%F_Under_PBLTop(I,J,L) = 1.0d0

            ! Fraction of PBL mass in layer L, will be normalized below
            State_Met%F_of_PBL(I,J,L) = State_Met%PEdge(I,J,L) - State_Met%PEdge(I,J,L+1)

            ! Update lower edge height, m
            Lower_Edge_Height = Lower_Edge_Height + State_Met%BXHeight(I,J,L)

         EndIf
         
       EndDo

       ! Fraction of PBL mass in layer L, now normalize to sum of 1
       State_Met%F_of_PBL(I,J,:) = State_Met%F_of_PBL(I,J,:) / State_Met%PBL_Thick(I,J)

       ! Error check
       IF ( ABS( SUM( State_Met%F_OF_PBL(I,J,:) ) - 1.e+0_fp) > 1.e-3_fp) THEN
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
    State_Met%PBL_MAX_L = MAXVAL( CEILING( State_Met%PBL_Top_L ) )

  END SUBROUTINE COMPUTE_PBL_HEIGHT
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
  SUBROUTINE TURBDAY( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
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
    USE TIME_MOD,       ONLY : GET_TS_CONV
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
    LOGICAL, SAVE :: FIRST = .TRUE.

    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    INTEGER            :: I,    J,  L
    INTEGER            :: LTOP, N,  NA,    nAdvect
    REAL(fp)           :: AA,   CC, CC_AA, AREA_M2, DTCONV

    ! Arrays
    REAL(fp)           :: A(State_Grid%NX,State_Grid%NY)
    REAL(fp)           :: DTC(State_Grid%NX,State_Grid%NY, &
                              State_Grid%NZ,State_Chm%nAdvect)

    ! Pointers
    REAL(fp), POINTER  :: AD(:,:,:)
    REAL(fp), POINTER  :: TC(:,:,:,:)

    !=================================================================
    ! TURBDAY begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! First-time initialization
    IF ( FIRST .and. Input_Opt%amIRoot ) THEN

       ! Echo info
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) 'T U R B D A Y  -- by Dale Allen, U. Md.'
       WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
       WRITE( 6, '(a)' ) 'Last Modification Date: 2/4/03'
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )

       ! Reset first time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Do the boundary layer mixing
    !=================================================================

    ! Initalize
    AD      => State_Met%AD        ! Dry air mass
    nAdvect =  State_Chm%nAdvect   ! # of advected species
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
       DO NA = 1, nAdvect

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
    AD => NULL()
    TC => NULL()

  END SUBROUTINE TURBDAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pbl_mix
!
! !DESCRIPTION: Subroutine INIT\_PBL\_MIX allocates and zeroes module arrays
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PBL_MIX( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
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
    !=================================================================
    ! INIT_PBL_MIX begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Error trap: in an ESMF environment, it is possible that this
    ! routine is called twice. No need to allocate arrays if already
    ! done so.
    IF ( ALLOCATED( IMIX ) ) RETURN

    ! Arrays
    ALLOCATE( IMIX( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pbl_mix_mod:IMIX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IMIX = 0

    ALLOCATE( FPBL( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pbl_mix_mod:FPBL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    FPBL = 0.0_fp

  END SUBROUTINE INIT_PBL_MIX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_pbl_mix
!
! !DESCRIPTION: Subroutine CLEANUP\_PBL\_MIX allocates and zeroes
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PBL_MIX
!
! !REVISION HISTORY:
!  11 Feb 2005 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_PBL_MIX begins here!
    !=================================================================
    IF ( ALLOCATED( IMIX        ) ) DEALLOCATE( IMIX        )
    IF ( ALLOCATED( FPBL        ) ) DEALLOCATE( FPBL        )

  END SUBROUTINE CLEANUP_PBL_MIX
!EOC
END MODULE PBL_MIX_MOD
