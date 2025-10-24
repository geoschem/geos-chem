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
  USE Dictionary_M, ONLY : dictionary_t
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Grid
  PUBLIC :: Allocate_State_Grid
  PUBLIC :: Register_State_Grid
  PUBLIC :: Cleanup_State_Grid
  PUBLIC :: Lookup_Grid
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Grid State
  !=========================================================================
  TYPE, PUBLIC :: GrdState

     !----------------------------------------------------------------------
     ! General grid fields
     !----------------------------------------------------------------------
     CHARACTER(LEN=255) :: GridRes          ! Grid resolution
     REAL(fp)           :: DX               ! Delta X         [deg longitude]
     REAL(fp)           :: DY               ! Delta Y         [degs latitude]
     REAL(fp)           :: XMin             ! Minimum X value [deg longitude]
     REAL(fp)           :: XMax             ! Maximum X value [deg longitude]
     REAL(fp)           :: YMin             ! Minimum Y value [degs latitude]
     REAL(fp)           :: YMax             ! Maximum Y value [degs latitude]
     INTEGER            :: NX               ! # of grid boxes in X-direction
     INTEGER            :: NY               ! # of grid boxes in Y-direction
     INTEGER            :: NZ               ! # of grid boxes in Z-direction
     LOGICAL            :: HalfPolar        ! Use half-sized polar boxes?
     LOGICAL            :: Center180        ! Is the Int'l Date Line a model
                                            !  midpoint (T) or edge (F)?
     LOGICAL            :: NestedGrid       ! Is it a nested grid sim?
     INTEGER            :: NorthBuffer      ! # buffer grid boxes on N edge
     INTEGER            :: SouthBuffer      ! # buffer grid boxes on S edge
     INTEGER            :: EastBuffer       ! # buffer grid boxes on E edge
     INTEGER            :: WestBuffer       ! # buffer grid boxes on W edge
     INTEGER            :: GlobalNX         ! NX on the global grid
     INTEGER            :: GlobalNY         ! NY on the global grid
     INTEGER            :: NativeNZ         ! NZ on the native-resolution grid
     INTEGER            :: MaxChemLev       ! Max # levels in chemistry grid
     INTEGER            :: MaxStratLev      ! Max # levels below strat
     INTEGER            :: MaxTropLev       ! Max # levels below trop
     INTEGER            :: XMinOffset       ! Min X offset from global grid
     INTEGER            :: XMaxOffset       ! Max X offset from global grid
     INTEGER            :: YMinOffset       ! Min Y offset from global grid
     INTEGER            :: YMaxOffset       ! Max Y offset from global grid
     REAL(f8),  POINTER :: GlobalXMid (:,:) ! Lon centers on global grid [deg]
     REAL(f8),  POINTER :: GlobalYMid (:,:) ! Lat centers on global grid [deg]
     REAL(f8),  POINTER :: GlobalXEdge(:,:) ! Lon centers on global grid [deg]
     REAL(f8),  POINTER :: GlobalYEdge(:,:) ! Lat centers on global grid [deg]
     REAL(f8),  POINTER :: XMid       (:,:) ! Lon centers [degrees]
     REAL(f8),  POINTER :: XEdge      (:,:) ! Lon edges   [degrees]
     REAL(f8),  POINTER :: YMid       (:,:) ! Lat centers [degrees]
     REAL(f8),  POINTER :: YEdge      (:,:) ! Lat edges   [degrees]
     REAL(f8),  POINTER :: YMid_R     (:,:) ! Lat centers [radians]
     REAL(f8),  POINTER :: YEdge_R    (:,:) ! Lat edges   [radians]
     REAL(f8),  POINTER :: YSIN       (:,:) ! SIN( lat edges )
     REAL(f8),  POINTER :: Area_M2    (:,:) ! Grid box area [m2]

     !----------------------------------------------------------------------
     ! Coordinate variables for GC-Classic History netCDF files
     !----------------------------------------------------------------------
     REAL(f4), POINTER  :: Area       (:,:) ! Surface area (REAL*4)
     REAL(f8), POINTER  :: HyAi       (:  ) ! Hybrid Ap at level interface
     REAL(f8), POINTER  :: HyAm       (:  ) ! Hybrid Ap at level midpoint
     REAL(f8), POINTER  :: HyBi       (:  ) ! Hybrid B  at level interface
     REAL(f8), POINTER  :: HyBm       (:  ) ! Hybrid B  at level midpoint
     REAL(f8), POINTER  :: ILev       (:  ) ! Level interface coordinate
     REAL(f8), POINTER  :: Lat        (:  ) ! Latitude centers
     REAL(f8), POINTER  :: LatBnd     (:,:) ! CF-compliant lat bounds
     REAL(f8), POINTER  :: LatE       (:  ) ! Latitude edges
     REAL(f8), POINTER  :: Lev        (:  ) ! Level midpoint coordinate
     REAL(f8), POINTER  :: Lon        (:  ) ! Longitude centers
     REAL(f8), POINTER  :: LonBnd     (:,:) ! Cf-compliant lon bounds
     REAL(f8), POINTER  :: LonE       (:  ) ! Longitude edges
     REAL(f8)           :: P0               ! Reference pressure (hPa)
     REAL(f8), POINTER  :: Time       (:  ) ! Time

#ifdef LUO_WETDEP
     !----------------------------------------------------------------------
     ! Fields needed for the Luo et al wet deposition scheme
     !----------------------------------------------------------------------
     REAL(fp),  POINTER :: DXSN_M     (:,:) ! Averange grid box width [m]
                                            !  at the S and N edges
     REAL(fp),  POINTER :: DYWE_M     (:,:) ! Averange grid box width [m]
                                            !  at the W and E edges
#endif

#if defined( MODEL_WRF ) || defined( MODEL_CESM )
     !----------------------------------------------------------------------
     ! Grid numbers for WRF and CESM, for each CPU to run multiple
     ! instances  of GEOS-Chem. These numbers are unique-per-core (local).
     ! A pair of (Input_Opt%thisCPU, State_Grid%CPU_Subdomain_ID) is needed
     ! to uniquely identify a geographical region.
     !----------------------------------------------------------------------

     ! Grid identifier number (local)
     ! WRF  : domain number
     ! CESM : chunk number/lchnk
     INTEGER            :: CPU_Subdomain_ID

     ! First grid identifier number (local) in this CPU
     INTEGER            :: CPU_Subdomain_FirstID
#endif

#ifdef MODEL_GEOS
     !----------------------------------------------------------------------
     ! NASA GEOS ESM-specific fields
     !----------------------------------------------------------------------
     LOGICAL            :: PredictorIsActive   ! Are we in the predictor step?
#endif

     !-----------------------------------------------------------------------
     ! Registry of variables contained within thje State_Grid object
     !-----------------------------------------------------------------------
     CHARACTER(LEN=4)            :: State     = 'GRID'   ! Name of this state
     TYPE(MetaRegItem),  POINTER :: Registry  => NULL()  ! Registry object
     TYPE(dictionary_t)          :: RegDict              ! Lookup table

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
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Mod,       ONLY : Registry_AddField
    USE Registry_Params_Mod
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
    ! Scalars
    INTEGER            :: I,      J,       L

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Variable, Desc, Units

    !========================================================================
    ! Init State_Grid begins here!
    !========================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
       " -> at Init_State_Grid (located in Headers/state_grid_mod.F90)"

    !========================================================================
    ! Zero scalar fields
    !========================================================================
    State_Grid%GridRes               =  ''
    State_Grid%DX                    =  0e+0_fp
    State_Grid%DY                    =  0e+0_fp
    State_Grid%XMin                  =  0e+0_fp
    State_Grid%XMax                  =  0e+0_fp
    State_Grid%YMin                  =  0e+0_fp
    State_Grid%YMax                  =  0e+0_fp
    State_Grid%NX                    =  0
    State_Grid%NY                    =  0
    State_Grid%NZ                    =  0
    State_Grid%HalfPolar             =  .FALSE.
    State_Grid%Center180             =  .FALSE.
    State_Grid%NestedGrid            =  .FALSE.
    State_Grid%NorthBuffer           =  0
    State_Grid%SouthBuffer           =  0
    State_Grid%EastBuffer            =  0
    State_Grid%WestBuffer            =  0
    State_Grid%GlobalNX              =  0
    State_Grid%GlobalNY              =  0
    State_Grid%NativeNZ              =  0
    State_Grid%MaxChemLev            =  0
    State_Grid%MaxStratLev           =  0
    State_Grid%MaxTropLev            =  0
    State_Grid%XMinOffset            =  0
    State_Grid%XMaxOffset            =  0
    State_Grid%YMinOffset            =  0
    State_Grid%YMaxOffset            =  0
    State_Grid%P0                    =  1000.0_f8 ! reference pressure, hPa
#if defined( MODEL_WRF ) || defined( MODEL_CESM )
    State_Grid%CPU_Subdomain_ID      =  -1
    State_Grid%CPU_Subdomain_FirstID =  -1
#endif
#if defined( MODEL_GEOS )
    State_Grid%PredictorIsActive     =  .FALSE.
#endif

    !========================================================================
    ! Nullify pointer array fields
    !========================================================================
    State_Grid%GlobalXMid            => NULL()
    State_Grid%GlobalYMid            => NULL()
    State_Grid%GlobalXEdge           => NULL()
    State_Grid%GlobalYEdge           => NULL()
    State_Grid%XMid                  => NULL()
    State_Grid%XEdge                 => NULL()
    State_Grid%YMid                  => NULL()
    State_Grid%YEdge                 => NULL()
    State_Grid%YMid_R                => NULL()
    State_Grid%YEdge_R               => NULL()
    State_Grid%YSIN                  => NULL()
    State_Grid%Area_M2               => NULL()
    State_Grid%Area                  => NULL()
    State_Grid%Time                  => NULL()
    State_Grid%HyAm                  => NULL()
    State_Grid%HyBm                  => NULL()
    State_Grid%Lev                   => NULL()
    State_Grid%HyAi                  => NULL()
    State_Grid%HyBi                  => NULL()
    State_Grid%ILev                  => NULL()
    State_Grid%Lat                   => NULL()
    State_Grid%LatE                  => NULL()
    State_Grid%LatBnd                => NULL()
    State_Grid%Lon                   => NULL()
    State_Grid%LonE                  => NULL()
    State_Grid%LonBnd                => NULL()
#ifdef LUO_WETDEP
    State_Grid%DXSN_M                => NULL()
    State_Grid%DYWE_M                => NULL()
#endif



   END SUBROUTINE Init_State_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate_State_Grid
!
! !DESCRIPTION: Subroutine ALLOCATE\_STATE\_GRID initializes variables and
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
! !REMARKS:
!  State_Grid fields are allocated here.  They will be registered
!  separately after the call to COMPUTE_GRID.
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
    ! Allocate_State_Grid begins here!
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS

    !======================================================================
    ! Allocate general grid fields
    !
    ! NOTE: State_Grid%GlobalXMid and State_Grid%GlobalYMid are allocated
    ! in gc_grid_mod.F90 after computing State_Grid%GlobalNX and
    ! State_Grid%GlobalNY
    !======================================================================
    ALLOCATE( State_Grid%Area_M2( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Area_M2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Area_M2 = 0.0_f8

    ALLOCATE( State_Grid%XEdge( State_Grid%NX+1, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%XEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%XEdge = 0.0_f8

    ALLOCATE( State_Grid%XMid( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%XMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%XMid = 0.0_f8

    ALLOCATE( State_Grid%YEdge( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YEdge', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YEdge = 0.0_f8

    ALLOCATE( State_Grid%YEdge_R( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YEdge_R', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YEdge_R = 0.0_f8

    ALLOCATE( State_Grid%YMid( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YMid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YMid = 0.0_f8

    ALLOCATE( State_Grid%YMid_R( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YMid_R', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YMid_R = 0.0_f8

    ALLOCATE( State_Grid%YSIN( State_Grid%NX, State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%YSIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%YSIN = 0.0_f8

#if !defined( MODEL_GCHPCTM ) && !defined( MODEL_GEOS )
    !========================================================================
    ! Allocate coordinate variables for GC-Classic History diagnostics
    !========================================================================
    ALLOCATE( State_Grid%Area( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Area_M2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Area_M2 = 0.0_f4

    ALLOCATE( State_Grid%HyAi( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%HyAi', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%HyAi = 0.0_f8

    ALLOCATE( State_Grid%HyAm( State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%HyAm', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%HyAm = 0.0_f8

    ALLOCATE( State_Grid%HyBi( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%HyBi', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%HyBi = 0.0_f8

    ALLOCATE( State_Grid%HyBm( State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%HyBm', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%HyBm = 0.0_f8

    ALLOCATE( State_Grid%ILev( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%ILev', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%ILev = 0.0_f8

    ALLOCATE( State_Grid%Lat( State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Lat', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Lat = 0.0_f8

    ALLOCATE( State_Grid%LatBnd( 2, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%LatBnd', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%LatBnd = 0.0_f8

    ALLOCATE( State_Grid%LatE( State_Grid%NY+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%LatE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%LatE = 0.0_f8

    ALLOCATE( State_Grid%Lev( State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Lev', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Lev = 0.0_f8

    ALLOCATE( State_Grid%Lon( State_Grid%NX ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Lon', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Lon = 0.0_f8

    ALLOCATE( State_Grid%LonBnd( 2, State_Grid%NX ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%LonBnd', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%LonBnd = 0.0_f8

    ALLOCATE( State_Grid%LonE( State_Grid%NX+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%LonE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%LonE = 0.0_f8

    ALLOCATE( State_Grid%Time( 1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%Time', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%Area_M2 = 0.0_f4
#endif

#ifdef LUO_WETDEP
    !========================================================================
    ! Allocate grid arrays for Luo et al wetdep
    !========================================================================
    ALLOCATE( State_Grid%DXSN_M( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%DXSN_M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%DXSN_M = 0e+0_fp

    ALLOCATE( State_Grid%DYWE_M( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'State_Grid%DYWE_M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Grid%DYWE_M = 0e+0_fp
#endif

  END SUBROUTINE Allocate_State_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_State_Grid
!
! !DESCRIPTION: Allocates and registers all module variables, which hold
!  horizontal and vertical grid information.  This will be used for netCDF
!  metadata in the History component
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_State_Grid( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Mod
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Register_State_Grid begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
       ' -> at Register_State_Grid (located in Headers/state_grid_mod.F90)'

    !========================================================================
    ! Register general grid fields
    !========================================================================

    !---------------------------
    ! State_Grid%Area_M2
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'AREAM2',                                          &
         Description    = 'Surface area',                                    &
         Units          = 'm2',                                              &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%Area_M2,                                &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%Area_M2', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%GlobalXEdge
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'GLOBALXEDGE',                                     &
         Description    = 'Global longitude edges',                          &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%GlobalXEdge,                            &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%GlobalXEdge', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%GlobalXMid
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'GLOBALXMID',                                      &
         Description    = 'Global longitude centers',                        &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%GlobalXMid,                             &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%GlobalXMid', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%GlobalYEdge
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'GLOBALYEDGE',                                     &
         Description    = 'Global latitude edges',                           &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%GlobalYEdge,                            &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%GlobalYEdge', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%GlobalYMid
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'GLOBALYMID',                                      &
         Description    = 'Global latitude centers',                         &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%GlobalYMid,                             &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%GlobalYMid', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%XEdge
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'XEDGE',                                           &
         Description    = 'Longitude edges',                                 &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%XEdge,                                  &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%XEdge', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%XMid
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'XMID',                                            &
         Description    = 'Longitude centers',                               &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%XMid,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%XMid', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%YEdge
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'YEDGE',                                           &
         Description    = 'Latitude edges',                                  &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%YEdge,                                  &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%YEdge', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%YEdge_R
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'YEDGER',                                          &
         Description    = 'Latitude edges',                                  &
         Units          = 'radians',                                         &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%YEdge_R,                                &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%YEdge_R', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%YMid
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'YMID',                                            &
         Description    = 'Latitude centers',                                &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%YMid,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%YMid', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%YMid_R
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'YMIDR',                                           &
         Description    = 'Latitude centers',                                &
         Units          = 'radians',                                         &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%YMid_R,                                 &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%YMid_R', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%YSIN
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'YSIN',                                            &
         Description    = 'Sine of latitude edges',                          &
         Units          = '1',                                               &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'xy',                                              &
         Data2d_8       = State_Grid%YSin,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%YSIN', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

#if !defined( MODEL_GCHPCTM ) && !defined( MODEL_GEOS )
    !========================================================================
    ! Register coordinate variables for GC-Classic History diagnostics
    ! (these may also be needed for WRF-GC)
    !========================================================================

    !---------------------------
    ! State_Grid%Area
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'AREA',                                            &
         Description    = 'Surface area',                                    &
         Units          = 'm2',                                              &
         Output_KindVal = KINDVAL_F4,                                        &
         DimNames       = 'xy',                                              &
         Data2d_4       = State_Grid%Area,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%Area', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%HyAi
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'HYAI',                                            &
         Description    = 'hybrid A coefficient at layer interfaces',        &
         Units          = 'hPa',                                             &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         Data1d_8       = State_Grid%HyAi,                                   &
         OnLevelEdges   = .TRUE.,                                            &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%HyAi', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%HyAm
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'HYAM',                                            &
         Description    = 'hybrid A coefficient at layer midpoints',         &
         Units          = 'hPa',                                             &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         Data1d_8       = State_Grid%HyAm,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%HyAm', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%HyBi
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'HYBI',                                            &
         Description    = 'hybrid B coefficient at layer interfaces',        &
         Units          = '1',                                               &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         Data1d_8       = State_Grid%HyBi,                                   &
         OnLevelEdges   = .TRUE.,                                            &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%HyBi', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%HyBm
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'HYBM',                                            &
         Description    = 'hybrid B coefficient at layer midpoints',         &
         Units          = '1',                                               &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         Data1d_8       = State_Grid%HyBm,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%HyBm', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%ILev
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'ILEV',                                            &
         Description    = 'hybrid level at interfaces ((A/P0)+B)',           &
         Units          = 'level',                                           &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         OnLevelEdges   = .TRUE.,                                            &
         Data1d_8       = State_Grid%ILev,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%ILev', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%Lat
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LAT',                                             &
         Description    = 'Latitude',                                        &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'y',                                               &
         Data1d_8       = State_Grid%Lat,                                    &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%Lat', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%LatBnd
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LATBND',                                          &
         Description    = 'Latitude bounds (CF-compliant)',                  &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'by',                                              &
         Data2d_8       = State_Grid%LatBnd,                                 &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%LatBnd', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%LatE
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LATE',                                            &
         Description    = 'Latitude edges',                                  &
         Units          = 'degrees_north',                                   &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'y',                                               &
         Data1d_8       = State_Grid%LatE,                                   &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%LatE', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%Lev
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LEV',                                             &
         Description    = 'hybrid level at midpoints ((A/P0)+B)',            &
         Units          = 'level',                                           &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'z',                                               &
         Data1d_8       = State_Grid%Lev,                                    &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%Lev', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%Lon
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LON',                                             &
         Description    = 'Longitude',                                       &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'x',                                               &
         Data1d_8       = State_Grid%Lon,                                    &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%Lat', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%LonBnd
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LONBND',                                          &
         Description    = 'Longitude bounds (CF-compliant)',                 &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'bx',                                              &
         Data2d_8       = State_Grid%LonBnd,                                 &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%LonBnd', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%LonE
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'LONE',                                            &
         Description    = 'Longitude edges',                                 &
         Units          = 'degrees_east',                                    &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 'x',                                               &
         Data1d_8       = State_Grid%Lon,                                    &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%LatE', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%P0
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'P0',                                              &
         Description    = 'reference pressure',                              &
         Units          = 'hPa',                                             &
         Output_KindVal = KINDVAL_F8,                                        &
         Data0d_8       = State_Grid%P0,                                     &
         RC             = RC                                                )

    CALL GC_CheckVar( 'State_Grid%P0', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------
    ! State_Grid%Time
    !---------------------------
    CALL Registry_AddField(                                                  &
         Input_Opt      = Input_Opt,                                         &
         Registry       = State_Grid%Registry,                               &
         State          = State_Grid%State,                                  &
         Variable       = 'TIME',                                            &
         Description    = 'Time',                                            &
         Units          = 'minutes since YYYY-MM-DD hh:mm:ss',               &
         Output_KindVal = KINDVAL_F8,                                        &
         DimNames       = 't',                                               &
         Data1d_8       = State_Grid%Time,                                   &
         RC             = RC                                                )

    ! Allocate
    CALL GC_CheckVar( 'State_Grid%Time', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    !========================================================================
    ! Once we are done registering all fields, we need to define the
    ! registry lookup table.  This algorithm will avoid hash collisions.
    !========================================================================
    CALL Registry_Set_LookupTable(                                           &
         Registry = State_Grid%Registry,                                     &
         RegDict  = State_Grid%RegDict,                                      &
         RC       = RC                                                      )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in routine "Registry_Set_LookupTable"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Print list of fields
    !========================================================================
    CALL Print_Grid( Input_Opt, State_Grid, RC, ShortFormat=.TRUE. )

    ! Write spacer line for log file
    WRITE( 6, '(a)' )

  END SUBROUTINE Register_State_Grid
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
! USES:
!
    USE Registry_Mod, ONLY : Registry_Destroy
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
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Cleanup_State_Grid begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
       ' -> at Cleanup_State_Grid (located in Headers/state_grid_mod.F90)'

    !========================================================================
    ! Deallocate general grid fields
    !========================================================================
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

    IF ( ASSOCIATED( State_Grid%GlobalXEdge ) ) THEN
       DEALLOCATE( State_Grid%GlobalXEdge, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%GlobalXEdge', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%GlobalXEdge => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%GlobalYEdge ) ) THEN
       DEALLOCATE( State_Grid%GlobalYEdge, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%GlobalYEdge', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%GlobalYEdge => NULL()
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
       CALL GC_CheckVar( 'State_Grid %Area_M2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Area_M2 => NULL()
    ENDIF

#if !defined( MODEL_GCHPCTM ) && !defined( MODEL_GEOS )
    !========================================================================
    ! Deallocate coordinate variables for GC-Classic History diagnostics
    ! (These fields may also be needed for WRF-GC)
    !========================================================================
    IF ( ASSOCIATED( State_Grid%Area ) ) THEN
       DEALLOCATE( State_Grid%Area, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Area_M2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Area_M2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%HyAi ) ) THEN
       DEALLOCATE( State_Grid%HyAi, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%HyAi', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%HyAi => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%HyAm ) ) THEN
       DEALLOCATE( State_Grid%HyAm, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%HyAm', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%HyAm => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%HyBi ) ) THEN
       DEALLOCATE( State_Grid%HyBi, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%HyBi', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%HyBi => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%HyBm ) ) THEN
       DEALLOCATE( State_Grid%HyBm, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%HyBm', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%HyBm => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%ILev ) ) THEN
       DEALLOCATE( State_Grid%ILev, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%ILev', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%ILev => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%Lat ) ) THEN
       DEALLOCATE( State_Grid%Lat, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Lat', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Lat => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%LatBnd ) ) THEN
       DEALLOCATE( State_Grid%LatBnd, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%LatBnd', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%LatBnd => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%LatE ) ) THEN
       DEALLOCATE( State_Grid%LatE, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%LatE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%LatE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%Lev ) ) THEN
       DEALLOCATE( State_Grid%Lev, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Lev', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Lev => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%Lon ) ) THEN
       DEALLOCATE( State_Grid%Lon, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Lon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Lon => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%LonBnd ) ) THEN
       DEALLOCATE( State_Grid%LonBnd, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%LonBnd', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%LonBnd => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%LonE ) ) THEN
       DEALLOCATE( State_Grid%LonE, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%LonE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%LonE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Grid%Time ) ) THEN
       DEALLOCATE( State_Grid%Time, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%Time', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%Area_M2 => NULL()
    ENDIF
#endif

#ifdef LUO_WETDEP
    !========================================================================
    ! Deallocate grid fields for Luo et al wetdep
    !========================================================================
    IF ( ASSOCIATED( State_Grid%DXSN_M ) ) THEN
       DEALLOCATE( State_Grid%DXSN_M, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%DXSN_M', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%DXSN_M => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Grid%DYWE_M ) ) THEN
       DEALLOCATE( State_Grid%DYWE_M, STAT=RC )
       CALL GC_CheckVar( 'State_Grid%DYWE_M', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Grid%DYWE_M => NULL()
    ENDIF
#endif

    !========================================================================
    ! Destroy the registry of fields for this module
    !========================================================================
    CALL Registry_Destroy( State_Grid%Registry, State_Grid%RegDict, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object "Registry"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       State_Grid%Registry => NULL()
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_State_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Grid
!
! !DESCRIPTION: Print information about all the registered variables
!  contained within the gc\_grid\_mod.F90 module. This is basically a wrapper
!  for routine REGISTRY\_PRINT in registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Grid( Input_Opt, State_Grid, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Registry_Mod,  ONLY : Registry_Print
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    LOGICAL,        OPTIONAL    :: ShortFormat ! Print truncated info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success/failure?
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Print_GC_Grid (in GeosUtil/grid_registry_mod.F90)'

    ! Only print information on the root CPU
    IF ( .not. Input_Opt%amIRoot ) RETURN

    !========================================================================
    ! Print info about registered variables
    !========================================================================
    IF ( Input_Opt%amIRoot .and. Input_Opt%Verbose ) THEN

       ! Header line
       WRITE( 6, 10 )
 10    FORMAT( /, &
         'Registered variables contained within the State_Grid object:'     )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )

       ! Print registry info in truncated format
       CALL Registry_Print( Input_Opt   = Input_Opt,                         &
                            Registry    = State_Grid%Registry,               &
                            ShortFormat = ShortFormat,                       &
                            RC          = RC                                )

       ! Trap error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "Registry_Print"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Print_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lookup_Grid
!
! !DESCRIPTION: Return metadata and/or a pointer to the data for any
!  variable contained within the GRID registry by searching for its name.
!  This is basically a wrapper for routine REGISTRY\_LOOKUP in
!  registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Lookup_Grid( Input_Opt,       State_Grid,      Variable,        &
                          RC,              Description,     Dimensions,      &
                          Source_KindVal,  Output_KindVal,  MemoryInKb,      &
                          Rank,            Units,           OnLevelEdges,    &
                          Ptr0d_8,         Ptr1d_8,         Ptr2d_4,         &
                          Ptr2d_8                                           )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE Registry_Mod,  ONLY : Registry_Lookup
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt       ! Input Options object
    TYPE(GrdState),      INTENT(IN)  :: State_Grid      ! Grid State object
    CHARACTER(LEN=*),    INTENT(IN)  :: Variable        ! Variable name
!
! !OUTPUT PARAMETERS:
!
    ! Required outputs
    INTEGER,             INTENT(OUT) :: RC              ! Success or failure?

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL    :: Description     ! Description of data
    INTEGER,             OPTIONAL    :: Dimensions(3)   ! Dimensions of data
    INTEGER,             OPTIONAL    :: Source_KindVal  ! KIND value of data
    INTEGER,             OPTIONAL    :: Output_KindVal  ! KIND value of data
    REAL(fp),            OPTIONAL    :: MemoryInKb      ! Memory usage
    INTEGER,             OPTIONAL    :: Rank            ! Size of data
    CHARACTER(LEN=255),  OPTIONAL    :: Units           ! Units of data
    LOGICAL,             OPTIONAL    :: OnLevelEdges    ! =T if data is defined
                                                        !  on vertical grid
                                                        !  edges; =F if center

    ! Pointers to data
    REAL(f8),   POINTER, OPTIONAL    :: Ptr0d_8         ! 0D 8-byte data
    REAL(f8),   POINTER, OPTIONAL    :: Ptr1d_8(:  )    ! 1D 8-byte data
    REAL(f4),   POINTER, OPTIONAL    :: Ptr2d_4(:,:)    ! 2D 4-byte data
    REAL(f8),   POINTER, OPTIONAL    :: Ptr2d_8(:,:)    ! 2D 8-byte data
!
! !REMARKS:
!  We keep the StateName variable private to this module. Users only have
!  to supply the name of each module variable.
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Lookup_Grid (in GeosUtil/grid_registry_mod.F90)'

    !========================================================================
    ! Look up a variable; Return metadata and/or a pointer to the data
    !========================================================================
    CALL Registry_Lookup( am_I_Root      = Input_Opt%amIRoot,                &
                          Registry       = State_Grid%Registry,              &
                          RegDict        = State_Grid%RegDict,               &
                          State          = State_Grid%State,                 &
                          Variable       = Variable,                         &
                          Description    = Description,                      &
                          Dimensions     = Dimensions,                       &
                          Source_KindVal = Source_KindVal,                   &
                          Output_KindVal = Output_KindVal,                   &
                          MemoryInKb     = MemoryInKb,                       &
                          Rank           = Rank,                             &
                          Units          = Units,                            &
                          OnLevelEdges   = OnLevelEdges,                     &
                          Ptr0d_8        = Ptr0d_8,                          &
                          Ptr1d_8        = Ptr1d_8,                          &
                          Ptr2d_8        = Ptr2d_8,                          &
                          Ptr2d_4        = Ptr2d_4,                          &
                          RC             = RC                               )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable "' // TRIM( Variable ) //           &
                '" in the State_Grid registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_Grid
!EOC
END MODULE State_Grid_Mod
