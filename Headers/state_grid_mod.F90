!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: state_grid_mod.F90
!
! !DESCRIPTION: Module STATE\_GRID\_MOD contains the derived type
!  used to define the Grid State object for GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE State_Grid_Mod
!
! USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Grid
  PUBLIC :: Allocate_State_Grid
  PUBLIC :: Cleanup_State_Grid
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Grid State
  !=========================================================================
  TYPE, PUBLIC :: GrdState

     !----------------------------------------
     ! User-defined grid fields
     !----------------------------------------
     CHARACTER(LEN=255) :: GridRes     ! Grid resolution
     REAL(fp)           :: DX          ! Delta X         [degrees longitude]
     REAL(fp)           :: DY          ! Delta Y         [degrees latitude]
     REAL(fp)           :: XMin        ! Minimum X value [degrees longitude]
     REAL(fp)           :: XMax        ! Maximum X value [degrees longitude]
     REAL(fp)           :: YMin        ! Minimum Y value [degrees latitude]
     REAL(fp)           :: YMax        ! Maximum Y value [degrees latitude]
     INTEGER            :: NX          ! # of grid boxes in X-direction
     INTEGER            :: NY          ! # of grid boxes in Y-direction
     INTEGER            :: NZ          ! # of grid boxes in Z-direction
     LOGICAL            :: HalfPolar   ! Use half-sized polar boxes?
     LOGICAL            :: NestedGrid  ! Is it a nested grid sim?
     INTEGER            :: NorthBuffer ! # buffer grid boxes on North edge
     INTEGER            :: SouthBuffer ! # buffer grid boxes on South edge
     INTEGER            :: EastBuffer  ! # buffer grid boxes on East edge
     INTEGER            :: WestBuffer  ! # buffer grid boxes on West edge

     !----------------------------------------
     ! Grid fields computed in gc_grid_mod.F90
     !----------------------------------------
     INTEGER            :: GlobalNX    ! NX on the global grid
     INTEGER            :: GlobalNY    ! NY on the global grid
     INTEGER            :: NativeNZ    ! NZ on the native-resolution grid
     INTEGER            :: MaxChemLev  ! Max # levels in chemistry grid
     INTEGER            :: MaxStratLev ! Max # levels below strat
     INTEGER            :: MaxTropLev  ! Max # levels below trop
     INTEGER            :: XMinOffset  ! X offset from global grid
     INTEGER            :: XMaxOffset  ! X offset from global grid
     INTEGER            :: YMinOffset  ! Y offset from global grid
     INTEGER            :: YMaxOffset  ! Y offset from global grid

     ! Arrays
     REAL(fp),  POINTER :: GlobalXMid(:,:) ! Lon centers on global grid [deg]
     REAL(fp),  POINTER :: GlobalYMid(:,:) ! Lat centers on global grid [deg]
     REAL(fp),  POINTER :: XMid      (:,:) ! Lon centers [degrees]
     REAL(fp),  POINTER :: XEdge     (:,:) ! Lon edges   [degrees]
     REAL(fp),  POINTER :: YMid      (:,:) ! Lat centers [degrees]
     REAL(fp),  POINTER :: YEdge     (:,:) ! Lat edges   [degrees]
     REAL(fp),  POINTER :: YMid_R    (:,:) ! Lat centers [radians]
     REAL(fp),  POINTER :: YEdge_R   (:,:) ! Lat edges   [radians]
     REAL(fp),  POINTER :: YSIN      (:,:) ! SIN( lat edges )
     REAL(fp),  POINTER :: Area_M2   (:,:) ! Grid box area [m2]

  END TYPE GrdState
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 Nov 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_State_Grid
!
! !DESCRIPTION: Subroutine INIT\_STATE\_GRID initializes all fields of
!  the Grid State derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Grid( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 Nov 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Assume success
    RC            =  GC_SUCCESS

    !----------------------------------------
    ! User-defined grid fields
    !----------------------------------------
    State_Grid%GridRes      = ''
    State_Grid%DX           = 0e+0_fp
    State_Grid%DY           = 0e+0_fp
    State_Grid%XMin         = 0e+0_fp
    State_Grid%XMax         = 0e+0_fp
    State_Grid%YMin         = 0e+0_fp
    State_Grid%YMax         = 0e+0_fp
    State_Grid%NX           = 0
    State_Grid%NY           = 0
    State_Grid%NZ           = 0
    State_Grid%HalfPolar    = .FALSE.
    State_Grid%NestedGrid   = .FALSE.
    State_Grid%NorthBuffer  = 0
    State_Grid%SouthBuffer  = 0
    State_Grid%EastBuffer   = 0
    State_Grid%WestBuffer   = 0

    !----------------------------------------
    ! Grid fields computed in gc_grid_mod.F90
    !----------------------------------------
    State_Grid%GlobalNX     = 0
    State_Grid%GlobalNY     = 0
    State_Grid%NativeNZ     = 0
    State_Grid%MaxChemLev   = 0
    State_Grid%MaxStratLev  = 0
    State_Grid%MaxTropLev   = 0
    State_Grid%XMinOffset   = 0
    State_Grid%XMaxOffset   = 0
    State_Grid%YMinOffset   = 0
    State_Grid%YMaxOffset   = 0

    !---------------------------------------------------------------
    ! Nullify all fields for safety's sake before allocating them
    !---------------------------------------------------------------
    State_Grid%GlobalXMid   => NULL()
    State_Grid%GlobalYMid   => NULL()
    State_Grid%XMid         => NULL()
    State_Grid%XEdge        => NULL()
    State_Grid%YMid         => NULL()
    State_Grid%YEdge        => NULL()
    State_Grid%YMid_R       => NULL()
    State_Grid%YEdge_R      => NULL()
    State_Grid%YSIN         => NULL()
    State_Grid%Area_M2      => NULL()

   END SUBROUTINE Init_State_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_State_Grid
!
! !DESCRIPTION: Subroutine ALLOCATE\_STATE\_GRID initializes variables and!
!  allocates module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Allocate_State_Grid( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  10 Mar 2019 - M. Sulprizio- Initial version, based on Init_Grid formerly in
!                              gc_grid_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: AS

    !======================================================================
    ! Initialize module variables
    !======================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! NOTE: State_Grid%GlobalXMid and State_Grid%GlobalYMid are allocated
    ! in gc_grid_mod.F90 after computing State_Grid%GlobalNX and
    ! State_Grid%GlobalNY

    ALLOCATE( State_Grid%XMid( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%XMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%XMid = 0e+0_fp

    ALLOCATE( State_Grid%XEdge( State_Grid%NX+1, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%XEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%XEdge = 0e+0_fp

    ALLOCATE( State_Grid%YMid( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YMid = 0e+0_fp

    ALLOCATE( State_Grid%YEdge( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YEdge = 0e+0_fp

    ALLOCATE( State_Grid%YMid_R( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YMid_R', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YMid_R = 0e+0_fp

    ALLOCATE( State_Grid%YEdge_R( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YEdge_R', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YEdge_R = 0e+0_fp

    ALLOCATE( State_Grid%YSIN( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YSIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YSIN = 0e+0_fp

    ALLOCATE( State_Grid%Area_M2( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Area_M2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Area_M2 = 0e+0_fp

  END SUBROUTINE Allocate_State_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Grid
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_GRID deallocates all fields
!  of the grid state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Grid( State_Grid, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Return code
!
! !REVISION HISTORY:
!  11 Nov 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC      = GC_SUCCESS

    !=======================================================================
    ! Deallocate arrays
    !=======================================================================
    IF ( ASSOCIATED( State_Grid%GlobalXMid ) ) THEN
       DEALLOCATE( State_Grid%GlobalXMid, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%GlobalXMid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%GlobalXMid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%GlobalYMid ) ) THEN
       DEALLOCATE( State_Grid%GlobalYMid, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%GlobalYMid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%GlobalYMid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%XMid ) ) THEN
       DEALLOCATE( State_Grid%XMid, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%XMid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%XMid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%XEdge ) ) THEN
       DEALLOCATE( State_Grid%XEdge, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%XEdge', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%XEdge => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%YMid ) ) THEN
       DEALLOCATE( State_Grid%YMid, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%XMid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%XMid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%YEdge ) ) THEN
       DEALLOCATE( State_Grid%YEdge, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%YEdge', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%YEdge => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%YMid_R ) ) THEN
       DEALLOCATE( State_Grid%Ymid_R, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%YMid_R', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%YMid_R => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%YEdge_R ) ) THEN
       DEALLOCATE( State_Grid%YEdge_R, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%YEdge_R', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%YEdge_R => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%YSIN ) ) THEN
       DEALLOCATE( State_Grid%YSIN, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%YSIN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%YSIN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%Area_M2 ) ) THEN
       DEALLOCATE( State_Grid%Area_M2, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Area_M2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Area_M2 => NULL()
    ENDIF

  END SUBROUTINE Cleanup_State_Grid
!EOC
END MODULE State_Grid_Mod
