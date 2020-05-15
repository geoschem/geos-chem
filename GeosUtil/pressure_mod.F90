!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pressure_mod.F90
!
! !DESCRIPTION: Module PRESSURE\_MOD contains variables and routines which
!  specify the grid box pressures for both hybrid or pure-sigma models.
!  This is necessary for running GEOS-Chem with the hybrid grids.
!\\
!\\
! !INTERFACE:
!
MODULE PRESSURE_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GET_AP
  PUBLIC  :: GET_BP
  PUBLIC  :: GET_PEDGE                ! wet air P at lower grid edge
  PUBLIC  :: GET_PCENTER              ! wet air P at grid center
  PUBLIC  :: GET_PEDGE_FULLGRID
  PUBLIC  :: GET_PEDGE_DRY
  PUBLIC  :: GET_DELP_DRY
  PUBLIC  :: INIT_PRESSURE
  PUBLIC  :: SET_FLOATING_PRESSURES
  PUBLIC  :: CLEANUP_PRESSURE
#if defined( ESMF_ ) || defined( MODEL_ )
  PUBLIC  :: Accept_External_Pedge
#endif
#if defined( MODEL_WRF )
  PUBLIC  :: Accept_External_ApBp
#endif
!
! !REMARKS:
!
!  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 2/2/12)
!  ============================================================================
!
!  The pressure at the bottom edge of grid box (I,J,L) is defined as follows:
!                                                                             .
!     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
!                                                                             .
!  where
!                                                                             .
!     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
!     Ap(L)         has the same units as surface pressure [hPa]
!     Bp(L)         is  a unitless constant given at level edges
!                                                                             .
!  Ap(L) and Bp(L) are given to us by GMAO.
!                                                                             .
!  The following are true:
!  ----------------------------------------------------------------------------
!  (1) Bp(NZ+1) = 0.0       (L=NZ+1 is the atmosphere top)
!  (2) Bp(1)    = 1.0       (L=1    is the surface       )
!  (3) PTOP     = Ap(NZ+1)  (L=NZ+1 is the atmosphere top)
!
! !REVISION HISTORY:
!  27 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  REAL(fp), ALLOCATABLE :: AP(:)                  ! "A" term for hybrid grid
  REAL(fp), ALLOCATABLE :: BP(:)                  ! "B" term for hybrid grid
  REAL(fp), ALLOCATABLE :: PFLT_DRY(:,:)          ! "Floating" dry sfc pres
  REAL(fp), ALLOCATABLE :: PFLT_WET(:,:)          ! "Floating" wet sfc pres
  REAL(fp), ALLOCATABLE :: AP_FULLGRID(:)         ! "A" term for full grid
  REAL(fp), ALLOCATABLE :: BP_FULLGRID(:)         ! "B" term for full grid
#if defined( ESMF_ ) || defined( MODEL_ )
  REAL(fp), ALLOCATABLE :: EXTERNAL_PEDGE(:,:,:)  ! Pressure edges from
                                                  !  external grid
#endif

  REAL(fp)              :: a132_loc(133), b132_loc(133)

! GEOS-5 132 levels (from /src/GMAO_Shared/GMAO_hermes/m_set_eta.F90
! --------------
  data a132_loc /  1.000000,     1.996276,     3.093648,     4.651099,    &
                   6.804155,     9.711212,     13.553898,    18.536953,   &
                  24.887674,    32.854966,    42.708057,    54.734916,    &
                  69.240493,    86.544776,    106.980758,   130.892382,   &
                  158.632424,   190.560538,   227.041195,   268.441904,   &
                  315.131439,   367.478204,   425.848769,   490.606509,   &
                  562.110455,   640.714290,   726.765342,   820.603888,   &
                  922.562490,   1032.965616,  1152.128995,  1280.359406,  &
                  1417.954457,  1565.202880,  1722.383803,  1889.767115,  &
                  2067.613829,  2256.175598,  2455.695564,  2666.408361,  &
                  2888.539866,  3122.308425,  3367.924596,  3625.591648,  &
                  3895.506041,  4177.787642,  4472.464900,  4779.536600,  &
                  5098.971133,  5430.705281,  5774.647623,  6130.914868,  &
                  6500.271455,  6883.621876,  7281.985387,  7695.829790,  &
                  8126.006088,  8573.341452,  9039.303976,  9523.598485,  &
                  10024.837122, 10541.370406, 11071.225963, 11612.410025, &
                  12161.636274, 12714.691534, 13270.207397, 13824.594107, &
                  14373.151226, 14914.405313, 15444.869700, 15960.611311, &
                  16459.769620, 16939.268383, 17396.217121, 17828.450893, &
                  18233.600515, 18609.343488, 18953.501254, 19264.447677, &
                  19539.848583, 19778.217887, 19977.939176, 20137.018678, &
                  20254.734748, 20328.875760, 20358.523606, 20342.231101, &
                  20278.589963, 20166.744330, 20004.982477, 19792.792832, &
                  19528.424768, 19211.380327, 18840.138412, 18414.132983, &
                  17933.325139, 17400.426408, 16819.657745, 16195.578563, &
                  15532.946677, 14837.558610, 14115.393726, 13372.886551, &
                  12616.479397, 11852.696266, 11087.800514, 10327.790957, &
                  9578.207359,  8844.157660,  8129.832058,  7440.098773,  &
                  6777.003948,  6143.217998,  5541.186971,  4972.725810,  &
                  4438.905073,  3940.077056,  3475.984433,  3045.886238,  &
                  2648.697264,  2283.946319,  1951.862407,  1652.526827,  &
                  1385.902714,  1151.874101,  950.288155,   780.991556,   &
                  643.875906,   538.919476,   466.225293,   426.071190,   &
                  0.000000      /

  data b132_loc / 0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000000, &
                  0.000000, 0.000000, 0.000000, 0.000007, &
                  0.000024, 0.000059, 0.000112, 0.000198, &
                  0.000339, 0.000560, 0.000886, 0.001347, &
                  0.001984, 0.002845, 0.003955, 0.005356, &
                  0.007104, 0.009223, 0.011758, 0.014755, &
                  0.018243, 0.022264, 0.026854, 0.032044, &
                  0.037871, 0.044366, 0.051561, 0.059484, &
                  0.068168, 0.077639, 0.087925, 0.099055, &
                  0.111049, 0.123939, 0.137748, 0.152499, &
                  0.168220, 0.184930, 0.202659, 0.221424, &
                  0.241254, 0.262166, 0.284188, 0.307337, &
                  0.331578, 0.356790, 0.382792, 0.409444, &
                  0.436599, 0.464098, 0.491782, 0.519487, &
                  0.547056, 0.574335, 0.601181, 0.627461, &
                  0.653056, 0.677861, 0.701765, 0.724759, &
                  0.746767, 0.767710, 0.787535, 0.806224, &
                  0.823790, 0.840276, 0.855742, 0.870260, &
                  0.883905, 0.896733, 0.908781, 0.920085, &
                  0.930681, 0.940600, 0.949868, 0.958500, &
                  0.966498, 0.973850, 0.980526, 0.986474, 1.000000 /

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Ap
!
! !DESCRIPTION: Function GET\_AP returns the "A" term [hPa] for the
!  hybrid ETA coordinate.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_AP( L ) RESULT( AP_TEMP )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(fp)            :: AP_TEMP  ! Corresponding "A" value [hPa]
                                    !  at bottom edge of level L
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    AP_TEMP = AP(L)

  END FUNCTION GET_AP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Bp
!
! !DESCRIPTION: Function GET\_BP returns the "B" term [unitless] for the
!  hybrid ETA coordinate.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_BP( L ) RESULT( BP_TEMP )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(fp)            :: BP_TEMP  ! Corresponding "B" value [unitless]
                                    !  at bottom edge of level L
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    BP_TEMP = BP(L)

  END FUNCTION GET_BP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Floating_Pressures
!
! !DESCRIPTION: Subroutine SET\_FLOATING\_PRESSURES initializes the
!  dry and wet floating pressure fields PFLT\_DRY and PFLT\_WET with the 
!  "true" surface pressures PSC2\_DRY and PSC2\_WET, stored in State\_Met.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
!
! !USES:
!
    USE ERROR_MOD, ONLY : CHECK_VALUE
    USE ErrCode_Mod
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology state object
!
! !OUTPUT ARGUMENTS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!   The surface pressures PSC2_DRY and PSC2_WET represent the most recently
!   interpolated values derived from GMAO instantaneous atmospheric pressure
!   at the surface (including moisture).
!
! !REVISION HISTORY:
!  21 Jun 2016 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L
    INTEGER            :: ERR_LOC(4)
    CHARACTER(LEN=255) :: ERR_VAR
    CHARACTER(LEN=255) :: ERR_MSG
    REAL(fp)           :: PEDGE1, PEDGE2, SPHU_KGKG

    !=================================================================
    ! SET_FLOATING_PRESSURES begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    !! DEBUGGING (ewl)
    !PRINT *, " "
    !PRINT *, "In SET_FLOATING_PRESSURES"
    !PRINT *, "   Old PFLT_DRY(56,20): ", PFLT_DRY(56,20)
    !PRINT *, "   Old PFLT_WET(56,20): ", PFLT_WET(56,20)
    !! END DEBUGGING

    ! Set PFLT_DRY equal to input value PS
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, ERR_LOC, ERR_VAR, ERR_MSG )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Set the floating pressures to the most recently interpolated
       ! instantaneous pressures
       PFLT_DRY(I,J) = State_Met%PSC2_DRY(I,J)
       PFLT_WET(I,J) = State_Met%PSC2_WET(I,J)

       ! Check for NaN or Infinities in PFLT_DRY and PFLT_WET
       ERR_LOC = (/ I, J, 0, 0 /)
       ERR_VAR = 'PFLT_DRY'
       ERR_MSG = 'set_floating_pressures:1'
       CALL CHECK_VALUE( PFLT_DRY(I,J), ERR_LOC, ERR_VAR, ERR_MSG )
       ERR_VAR = 'PFLT_WET'
       ERR_MSG = 'set_floating_pressures:2'
       CALL CHECK_VALUE( PFLT_WET(I,J), ERR_LOC, ERR_VAR, ERR_MSG )

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !! DEBUGGING (ewl)
    !PRINT *, "   New PFLT_DRY(56,20): ", PFLT_DRY(56,20)
    !PRINT *, "   New PFLT_WET(56,20): ", PFLT_WET(56,20)
    !! END DEBUGGING

  END SUBROUTINE SET_FLOATING_PRESSURES
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Pedge
!
! !DESCRIPTION: Function GET\_PEDGE returns the pressure at the bottom edge
!  of level L.  L=1 is the surface, L=State\_Grid%NZ+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_PEDGE( I, J, L ) RESULT( PEDGE )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)   :: I        ! GEOS-Chem lon   index
    INTEGER, INTENT(IN)   :: J        ! GEOS-Chem lat   index
    INTEGER, INTENT(IN)   :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(f8)              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( ESMF_ ) || defined( MODEL_ )
    ! Pressure [hPa] at bottom edge of level L (see documentation header)
    ! Taken from the GCM fields
    PEDGE = EXTERNAL_PEDGE(I,J,L)
#else
    ! Pressure [hPa] at bottom edge of level L (see documentation header)
    ! Computed for use w/in GEOS-Chem
    PEDGE = AP(L) + ( BP(L) * PFLT_WET(I,J) )
#endif

  END FUNCTION GET_PEDGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Pcenter
!
! !DESCRIPTION: Function GET\_PCENTER returns the pressure at the vertical
!  midpoint of level L.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_PCENTER( I, J, L ) RESULT( PCENTER )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I        ! GEOS-Chem lon   index
    INTEGER, INTENT(IN) :: J        ! GEOS-Chem lat   index
    INTEGER, INTENT(IN) :: L        ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: PCENTER  ! Pressure @ center of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! The pressure at the center of a grid-box is found
    ! by averaging the pressures at the box's two edges
    PCENTER = 0.5e+0_fp * ( GET_PEDGE(I,J,L) + GET_PEDGE(I,J,L+1) )

  END FUNCTION GET_PCENTER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Pedge_Fullgrid
!
! !DESCRIPTION: Function GET\_PEDGE\_FULLGRID returns the pressure at the
!  bottom edge of level L of the unreduced vertical grid.  L=1 is the surface,
!  L=LState\_Grid%NZ+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_PEDGE_FULLGRID( I, J, L ) RESULT( PEDGE )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
    INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
    INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(fp)              :: PEDGE  ! Pressure @ bottom edge of (I,J,L) [hPa]
!
! !REVISION HISTORY:
!  (1 ) Modified from GET_PEDGE (cdh, 1/22/09)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! GET_PEDGE_FULLGRID begins here!
    !=================================================================

    ! Pressure [hPa] at bottom edge of level L (see documentation header)
    PEDGE = AP_FULLGRID(L) + ( BP_FULLGRID(L) * PFLT_WET(I,J) )

  END FUNCTION GET_PEDGE_FULLGRID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Pedge_Dry
!
! !DESCRIPTION: Function GET\_PEDGE\_DRY returns the pressure at the
!  bottom edge of level L, reconstructed using the dry surface pressure. 
!  L=1 is the surface, L=State\_Grid%NZ+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_PEDGE_DRY( I, J, L ) RESULT( PEDGE_DRY )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
    INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
    INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(f8) :: PEDGE_DRY  ! Dry prssr @ bottom edge of (I,J,L) [hPa]
!
! !REMARKS:
!  Dry pressures at the edges calculated within this routine should not
!  be used as height proxies. Wet pressure edge should be used instead.
!
! !REVISION HISTORY:
!  16 Jun 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    PEDGE_DRY = AP(L) + ( BP(L) * PFLT_DRY(I,J) )

  END FUNCTION GET_PEDGE_DRY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Delp_Dry
!
! !DESCRIPTION: Function GET\_DELP\_DRY returns the delta dry pressure
!  between the bottom edge of level L and top edge of level L+1,
!  constructed using the dry surface pressure and A and B parameters.
!  L=1 is the surface, L=State\_Grid%NZ+1 is the atm top.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_DELP_DRY( I, J, L ) RESULT( DELP_DRY )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I      ! GEOS-Chem lon   index
    INTEGER, INTENT(IN) :: J      ! GEOS-Chem lat   index
    INTEGER, INTENT(IN) :: L      ! GEOS-Chem level index
!
! !RETURN VALUE:
!
    REAL(f8) :: DELP_DRY          ! Prssr difference [hPa] between
                                  ! bottom edge of (I,J,L) and
                                  ! bottom edge of (I,J,L+1)
!
! !REVISION HISTORY:
!  06 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)           :: PEDGE_DRY_BOT, PEDGE_DRY_TOP

    PEDGE_DRY_BOT = AP(L)   + ( BP(L)   * PFLT_DRY(I,J) )
    PEDGE_DRY_TOP = AP(L+1) + ( BP(L+1) * PFLT_DRY(I,J) )

    DELP_DRY = PEDGE_DRY_BOT - PEDGE_DRY_TOP

  END FUNCTION GET_DELP_DRY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Pressure
!
! !DESCRIPTION: Subroutine INIT\_PRESSURE allocates and initializes the AP
!  and BP arrays.  It must be called in "main.f", after SIGE is defined.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PRESSURE( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  27 Aug 2002 - D. Abbot, S. Wu, & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L

    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, nLev

    !=================================================================
    ! INIT_PRESSURE begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Init_Pressure (in GeosUtil/pressure_mod.F90)'

    ALLOCATE( PFLT_DRY( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:PFLT_DRY', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PFLT_DRY = 0e+0_fp

    ALLOCATE( PFLT_WET( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:PFLT_WET', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PFLT_WET = 0e+0_fp

    ALLOCATE( AP( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:AP', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AP = 1e+0_fp

    ALLOCATE( BP( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:BP', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BP = 0e+0_fp

    ALLOCATE( AP_FULLGRID( State_Grid%NativeNZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:AP_FULLGRID', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AP = 1e+0_fp

    ALLOCATE( BP_FULLGRID( State_Grid%NativeNZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:BP_FULLGRID', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BP = 0e+0_fp

#if defined( ESMF_ ) || defined( MODEL_ )
    ALLOCATE( EXTERNAL_PEDGE( State_Grid%NX, State_Grid%NY, State_Grid%NZ+1 ), &
              STAT=RC )
    CALL GC_CheckVar( 'vdiff_mod.F90:EXTERNAL_PEDGE', 2, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'EXTERNAL_PEDGE' )
    EXTERNAL_PEDGE = 0e+0_fp
#endif

    IF ( State_Grid%NZ == 47 ) THEN

       !-----------------------------------------------------------------
       ! 47-level reduced vertical grid
       !
       !  Bottom   Bottom    # levels
       !  edge of  edge prs  lumped
       !  level    (hPa)     together
       !
       !   PTOP       0.010
       !    47        0.066     4
       !    46        0.211     4
       !    45        0.617     4
       !    44        1.651     4
       !    43        4.077     4
       !    42        9.293     4
       !    41       19.792     4
       !    40       28.368     2
       !    39       40.175     2
       !    38       56.388     2
       !    37       78.512     2
       ! %%%% START LUMPING LEVELS ABOVE HERE %%%%%
       !    36       92.366
       !    35      108.663
       !    34      127.837
       !    33      150.393
       !    32      176.930
       ! %%%% FIXED-PRESSURE LEVELS BEGIN HERE %%%%
       !-----------------------------------------------------------------

       ! Ap [hPa] for 47 levels (48 edges)
       AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
               1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
               4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
               7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
               1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
               1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
               2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
               2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02, &
               1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
               7.851231d+01, 5.638791d+01, 4.017541d+01, 2.836781d+01, &
               1.979160d+01, 9.292942d+00, 4.076571d+00, 1.650790d+00, &
               6.167791d-01, 2.113490d-01, 6.600001d-02, 1.000000d-02 /)

       ! Bp [unitless] for 47 levels (48 edges)
       BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01, &
               9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01, &
               8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01, &
               7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01, &
               6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01, &
               4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01, &
               2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01, &
               6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00 /)

       !-----------------------------------------------------------------
       ! 72 level grid
       !-----------------------------------------------------------------

       ! Ap [hPa] for 72 levels (73 edges)
       AP_FULLGRID = &
            (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
               1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
               4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
               7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
               1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
               1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
               2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
               2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02, &
               1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
               7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01, &
               4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01, &
               1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01, &
               9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00, &
               4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00, &
               1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01, &
               6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01, &
               2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02, &
               6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02, &
               1.000000d-02 /)

       ! Bp [unitless] for 72 levels (73 edges)
       BP_FULLGRID = &
            (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01, &
               9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01, &
               8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01, &
               7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01, &
               6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01, &
               4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01, &
               2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01, &
               6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00 /)

    ELSE IF (State_Grid%NZ == 72 ) THEN

       !-----------------------------------------------------------------
       ! 72 level grid
       !-----------------------------------------------------------------

       ! Ap [hPa] for 72 levels (73 edges)
       AP = (/ 0.000000d+00, 4.804826d-02, 6.593752d+00, 1.313480d+01, &
               1.961311d+01, 2.609201d+01, 3.257081d+01, 3.898201d+01, &
               4.533901d+01, 5.169611d+01, 5.805321d+01, 6.436264d+01, &
               7.062198d+01, 7.883422d+01, 8.909992d+01, 9.936521d+01, &
               1.091817d+02, 1.189586d+02, 1.286959d+02, 1.429100d+02, &
               1.562600d+02, 1.696090d+02, 1.816190d+02, 1.930970d+02, &
               2.032590d+02, 2.121500d+02, 2.187760d+02, 2.238980d+02, &
               2.243630d+02, 2.168650d+02, 2.011920d+02, 1.769300d+02, &
               1.503930d+02, 1.278370d+02, 1.086630d+02, 9.236572d+01, &
               7.851231d+01, 6.660341d+01, 5.638791d+01, 4.764391d+01, &
               4.017541d+01, 3.381001d+01, 2.836781d+01, 2.373041d+01, &
               1.979160d+01, 1.645710d+01, 1.364340d+01, 1.127690d+01, &
               9.292942d+00, 7.619842d+00, 6.216801d+00, 5.046801d+00, &
               4.076571d+00, 3.276431d+00, 2.620211d+00, 2.084970d+00, &
               1.650790d+00, 1.300510d+00, 1.019440d+00, 7.951341d-01, &
               6.167791d-01, 4.758061d-01, 3.650411d-01, 2.785261d-01, &
               2.113490d-01, 1.594950d-01, 1.197030d-01, 8.934502d-02, &
               6.600001d-02, 4.758501d-02, 3.270000d-02, 2.000000d-02, &
               1.000000d-02 /)

       ! Bp [unitless] for 72 levels (73 edges)
       BP = (/ 1.000000d+00, 9.849520d-01, 9.634060d-01, 9.418650d-01, &
               9.203870d-01, 8.989080d-01, 8.774290d-01, 8.560180d-01, &
               8.346609d-01, 8.133039d-01, 7.919469d-01, 7.706375d-01, &
               7.493782d-01, 7.211660d-01, 6.858999d-01, 6.506349d-01, &
               6.158184d-01, 5.810415d-01, 5.463042d-01, 4.945902d-01, &
               4.437402d-01, 3.928911d-01, 3.433811d-01, 2.944031d-01, &
               2.467411d-01, 2.003501d-01, 1.562241d-01, 1.136021d-01, &
               6.372006d-02, 2.801004d-02, 6.960025d-03, 8.175413d-09, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00, 0.000000d+00, 0.000000d+00, 0.000000d+00, &
               0.000000d+00 /)

    ELSE IF ( State_Grid%NZ == 132 ) THEN

       !-----------------------------------------------------------------
       ! 132 level grid
       !-----------------------------------------------------------------

       DO L=1,State_Grid%NZ+1
          ! Ap [hPa] for 132 levels (133 edges)
          AP(L) = a132_loc(State_Grid%NZ+2-L)/100.0

          ! Bp [unitless] for 132 levels (133 edges)
          BP(L) = b132_loc(State_Grid%NZ+2-L)
       ENDDO

    ELSE

#if !defined( MODEL_WRF )
       WRITE( nLev, * ) State_Grid%NZ
       ErrMSg = 'Ap and Bp not defined for ' // TRIM( nLev ) // &
                ' levels. Please add these defintions in pressure_mod.F90.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
#endif

    ENDIF

#if ( !defined( ESMF_ ) && !defined( MODEL_ ) ) || defined( MODEL_GEOS )
    ! Echo info to std output (skip if interfacing with external models)
    IF ( Input_Opt%amIRoot .and. ( .not. Input_Opt%DryRun ) ) THEN
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
       WRITE( 6, '(a,/)' ) 'INIT_PRESSURE: Vertical coordinates!'
       WRITE( 6, '( ''Ap '', /, 6(f11.6,1x) )' ) AP(1:State_Grid%NZ+1)
       WRITE( 6, '(a)'   )
       WRITE( 6, '( ''Bp '', /, 6(f11.6,1x) )' ) BP(1:State_Grid%NZ+1)
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
    ENDIF
#endif

  END SUBROUTINE INIT_PRESSURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Pressure
!
! !DESCRIPTION: Subroutine CLEANUP\_PRESSURE deallocates all allocated arrays
!  at the end of a GEOS-Chem model run.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PRESSURE
!
! !REVISION HISTORY:
!  20 Aug 2002 - D. Abbot & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( AP          ) ) DEALLOCATE( AP          )
    IF ( ALLOCATED( BP          ) ) DEALLOCATE( BP          )
    IF ( ALLOCATED( AP_FULLGRID ) ) DEALLOCATE( AP_FULLGRID )
    IF ( ALLOCATED( BP_FULLGRID ) ) DEALLOCATE( BP_FULLGRID )
    IF ( ALLOCATED( PFLT_DRY    ) ) DEALLOCATE( PFLT_DRY    )
    IF ( ALLOCATED( PFLT_WET    ) ) DEALLOCATE( PFLT_WET    )
#if defined( ESMF_ ) || defined( MODEL_ )
    IF ( ALLOCATED( EXTERNAL_PEDGE ) ) DEALLOCATE( EXTERNAL_PEDGE )
#endif

  END SUBROUTINE CLEANUP_PRESSURE
!EOC
#if defined( ESMF_ ) || defined( MODEL_ )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Accept_External_Pedge
!
! !DESCRIPTION: Subroutine ACCEPT\_EXTERNAL\_PEDGE sets the GEOS-Chem
!  pressure edge variable with the values obtained from an external GCM
!  (such as the NASA GEOS-5 GCM).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Accept_External_Pedge( State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology state object
!
! !OUTPUT ARGUMENTS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is a setter for EXTERNAL_PEDGE.  It allows us to keep the
!  EXTERNAL_PEDGE array PRIVATE to this module, which is good programming
!  practice.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Set EXTERNAL_PEDGE to the pressure edges [hPa] carried in the
    ! State_Met object, which were obtained from the external GCM
    EXTERNAL_PEDGE = State_Met%PEDGE

    ! Return successfully
    RC             = GC_SUCCESS

  END SUBROUTINE Accept_External_Pedge
!EOC
#endif
#if defined ( MODEL_WRF )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Accept_External_ApBp
!
! !DESCRIPTION: Subroutine ACCEPT\_EXTERNAL\_ApBp sets the GEOS-Chem
!  hybrid grid AP, BP values with values obtained from an external model,
!  such as the WRF model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Accept_External_ApBp( State_Grid, ApIn, BpIn, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid            ! Grid State object
    REAL(fp),       INTENT(IN)  :: ApIn(State_Grid%NZ+1) ! "A" term for hybrid grid
    REAL(fp),       INTENT(IN)  :: BpIn(State_Grid%NZ+1) ! "B" term for hybrid grid
!
! !OUTPUT ARGUMENTS:
!
    INTEGER,        INTENT(OUT) :: RC              ! Success or failure?
!
! !REMARKS:
!  This routine is a setter for AP, BP.  It allows us to keep the
!  AP, BP array PRIVATE to this module, which is good programming
!  practice.
!
!  For WRF-GC, you need to enable the v3.9+ hybrid-sigma vertical coordinate system
!  in the WRF model. Set hybrid_opt = 2 in &dynamics, and ./configure -hyb.
!  Like WRF-GC itself, hybrid-sigma grids are experimental in WRF 3.9+.
!
! !REVISION HISTORY:
!  17 Aug 2018 - H.P. Lin    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    AP             = ApIn
    BP             = BpIn

    ! Return successfully
    RC             = GC_SUCCESS

  END SUBROUTINE Accept_External_ApBp
!EOC
#endif
END MODULE PRESSURE_MOD
