!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_registry_mod.F90
!
! !DESCRIPTION: Contains a registry of the GEOS-Chem horizontal and vertical
!  grid parameters.  The History Component uses this information to create
!  the metadata for the netCDF diagnostic files.
!\\
!\\
! !INTERFACE:
!
MODULE Grid_Registry_Mod
!
! !USES:
!
  USE Precision_Mod
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_Grid_Registry
  PUBLIC :: Cleanup_Grid_Registry
  PUBLIC :: Lookup_Grid
  PUBLIC :: Print_Grid
!
! !REMARKS:
!  Init_Grid_Registry has to be called AFTER both Init_Grid and Init_Pressure.
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Module variables
  REAL(f8), ALLOCATABLE, TARGET :: Area        (:,:)
  REAL(f8), ALLOCATABLE, TARGET :: Ap_Center   (:  )
  REAL(f8), ALLOCATABLE, TARGET :: B_Center    (:  )
  REAL(f8), ALLOCATABLE, TARGET :: Ap_Edge     (:  )
  REAL(f8), ALLOCATABLE, TARGET :: B_Edge      (:  )
  REAL(f8), ALLOCATABLE, TARGET :: Latitude    (:  )
  REAL(f8), ALLOCATABLE, TARGET :: Level_Center(:  )
  REAL(f8), ALLOCATABLE, TARGET :: Level_Edge  (:  )
  REAL(f8), ALLOCATABLE, TARGET :: Longitude   (:  )
  REAL(f8), ALLOCATABLE, TARGET :: Time        (:  )

  ! Registry of variables contained within gc_grid_mod.F90
  CHARACTER(LEN=4)              :: State     = 'GRID'   ! Name of this state
  TYPE(MetaRegItem),    POINTER :: Registry  => NULL()  ! Registry object

CONtAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Grid_Registry
!
! !DESCRIPTION: Allocates and registers all module variables, which hold
!  horizontal and vertical grid information.  This will be used for netCDF
!  metadata in the History component
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Grid_Registry( am_I_Root, RC )
!
! !USES:
!

    USE CMN_Size_Mod, ONLY: IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE GC_Grid_Mod
    USE Pressure_Mod   
    USE Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I,      J,       L

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Variable, Desc, Units

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Grid_Registry (in GeosUtil/grid_registry_mod.F90)'

    !======================================================================
    ! Allocate and register AREA
    !======================================================================
    Variable = 'AREA'
    Desc     = 'Surface area'
    Units    = 'm2'

    ! Allocate
    ALLOCATE( Area( IIPAR, JJPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_AREA', 0, RC )  
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       Area(I,J) = Get_Area_M2( I, J, 1 )
    ENDDO
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data2d      = Area,                              &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_AREA', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register AP_EDGE
    !======================================================================
    Variable = 'AP_EDGE'
    Desc     = 'Hybrid Ap at level edges'
    Units    = 'hPa'

    ! Allocate
    ALLOCATE( Ap_Edge( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_AP_EDGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR+1
       Ap_Edge(L) = Get_Ap( L )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Ap_Edge,                           &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_AP_EDGE', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register B_EDGE
    !======================================================================
    Variable = 'B_EDGE'
    Desc     = 'Hybrid B at level edges'
    Units    = '1'

    ! Allocate
    ALLOCATE( B_Edge( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_B_EDGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR+1
       B_Edge(L) = Get_Bp( L )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = B_Edge,                            &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_B_EDGE', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register AP_CENTER
    !======================================================================
    Variable = 'AP_CNTR'
    Desc     = 'Hybrid Ap at level centers'
    Units    = 'hPa'

    ! Allocate
    ALLOCATE( Ap_Center( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_AP_CNTR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR
       Ap_Center(L) = ( Ap_Edge(L) + Ap_Edge(L+1) ) * 0.5_f8
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Ap_Center,                         &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_AP_CNTR', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register B_EDGE
    !======================================================================
    Variable = 'B_CNTR'
    Desc     = 'Hybrid B at level centers'
    Units    = '1'

    ! Allocate
    ALLOCATE( B_Center( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_B_CNTR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR
       B_Center(L) = ( B_Edge(L) + B_Edge(L+1) ) * 0.5_f8
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = B_Center,                          &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_B_CNTR', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register TIME
    !======================================================================
    Variable = 'TIME'
    Desc     = 'Time'
    Units    = 'minutes since YYYY-MM-DD hh:mm:ss UTC'

    ! Allocate 
    ALLOCATE( Time( 1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_TIME', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    Time = 0.0_f8

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Time,                              &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_TIME', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LEV_CNTR
    !======================================================================
    Variable = 'LEV_CNTR'
    Desc     = 'atmosphere_hybrid_sigma_pressure_coordinate'
    Units    = 'a: ap b: b" ps: ps'

    ! Allocate
    ALLOCATE( Level_Center( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LEV_CNTR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    Level_Center = 0.0_fp

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Level_Center,                      &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_LEV_CNTR', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LEV_EDGE
    !======================================================================
    Variable = 'LEV_EDGE'
    Desc     = 'atmosphere_hybrid_sigma_pressure_coordinate'
    Units    = 'a: ap b: b" ps: ps'

    ! Allocate
    ALLOCATE( Level_Edge( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LEV_EDGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    Level_Edge = 0.0_fp

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Level_Edge,                        &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_LEV_EDGE', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LAT
    !======================================================================
    Variable = 'LAT'
    Desc     = 'Latitude'
    Units    = 'degrees_north'

    ! Allocate
    ALLOCATE( Latitude( JJPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO J = 1, JJPAR
       Latitude(J) = Get_YMid( 1, J, 1 )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Latitude,                          &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_LAT', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LON
    !======================================================================
    Variable = 'LON'
    Desc     = 'Longitude'
    Units    = 'degrees_east'

    ! Allocate
    ALLOCATE( Longitude( IIPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LON', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO I = 1, IIPAR
       Longitude(I) = Get_XMid( I, 1, 1 )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = Registry,                          &
                            State       = State,                             &
                            Variable    = Variable,                          &
                            Description = Desc,                              &
                            Units       = Units,                             &
                            Data1d      = Longitude,                         &
                            RC          = RC                                )

    CALL GC_CheckVar( 'GRID_LON', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Print list of fields
    !======================================================================
    CALL Print_Grid( am_I_Root, RC, ShortFormat=.TRUE. )

    ! Write spacer line for log file
    WRITE( 6, '(a)' )

  END SUBROUTINE Init_Grid_Registry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Grid_Registry 
!
! !DESCRIPTION: Deallocates all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Grid_Registry( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Destroy
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS: 
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
         ' -> at Cleanup_Grid_Registry (in GeosUtil/grid_registry_mod.F90)'

    !=======================================================================
    ! Deallocate fields
    !=======================================================================
    IF ( ALLOCATED( Area ) ) THEN
       DEALLOCATE( Area ) 
       CALL GC_CheckVar( 'GRID_AREA', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Ap_Center ) ) THEN
       DEALLOCATE( Ap_Center ) 
       CALL GC_CheckVar( 'GRID_AP_CNTR', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Ap_Edge ) ) THEN
       DEALLOCATE( Ap_Edge ) 
       CALL GC_CheckVar( 'GRID_AP_EDGE', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( B_Center ) ) THEN
       DEALLOCATE( B_Center ) 
       CALL GC_CheckVar( 'GRID_B_CNTR', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( B_Edge ) ) THEN
       DEALLOCATE( B_Edge ) 
       CALL GC_CheckVar( 'GRID_B_EDGE', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Latitude ) ) THEN
       DEALLOCATE( Latitude ) 
       CALL GC_CheckVar( 'GRID_LAT', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Level_Center ) ) THEN
       DEALLOCATE( Level_Center ) 
       CALL GC_CheckVar( 'GRID_LEV_CNTR', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Level_Edge ) ) THEN
       DEALLOCATE( Level_Edge ) 
       CALL GC_CheckVar( 'GRID_LEV_EDGE', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Longitude ) ) THEN
       DEALLOCATE( Longitude ) 
       CALL GC_CheckVar( 'GRID_LON', 3, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( Time ) ) THEN
       DEALLOCATE( Time ) 
       CALL GC_CheckVar( 'GRID_TIME', 3, RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object "Registry"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_Grid_Registry
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
  SUBROUTINE Print_Grid( am_I_Root, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Print
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Root CPU?  
    LOGICAL,        OPTIONAL    :: ShortFormat ! Print truncated info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success/failure?
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Print_GC_Grid (in GeosUtil/grid_registry_mod.F90)'

    !=======================================================================
    ! Print info about registered variables
    !=======================================================================

    ! Header line
    WRITE( 6, 10 )
 10 FORMAT( /, 'Registered variables contained within grid_registry_mod.F90:' )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    ! Print registry info in truncated format
    CALL Registry_Print( am_I_Root   = am_I_Root,           &
                         Registry    = Registry,            &
                         ShortFormat = ShortFormat,         &
                         RC          = RC                  )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Print"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
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
  SUBROUTINE Lookup_Grid( am_I_Root,   Variable,   RC,                       &
                          Description, Dimensions, KindVal,                  &
                          MemoryInKb,  Rank,       Units,                    &
                          Ptr1d,       Ptr2d                                )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Lookup
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root       ! Is this the root CPU? 
    CHARACTER(LEN=*),    INTENT(IN)  :: Variable        ! Variable name
!
! !OUTPUT PARAMETERS:
!
    ! Required outputs
    INTEGER,             INTENT(OUT) :: RC              ! Success or failure?

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL    :: Description     ! Description of data
    INTEGER,             OPTIONAL    :: Dimensions(3)   ! Dimensions of data
    INTEGER,             OPTIONAL    :: KindVal         ! Numerical KIND value
    REAL(fp),            OPTIONAL    :: MemoryInKb      ! Memory usage
    INTEGER,             OPTIONAL    :: Rank            ! Size of data
    CHARACTER(LEN=255),  OPTIONAL    :: Units           ! Units of data

    ! Pointers to data
    REAL(fp),   POINTER, OPTIONAL    :: Ptr1d(:  )      ! 1D flex-prec data
    REAL(fp),   POINTER, OPTIONAL    :: Ptr2d(:,:)      ! 2D flex-prec data
!
! !REMARKS:
!  We keep the StateName variable private to this module. Users only have
!  to supply the name of each module variable.
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Lookup_Grid (in GeosUtil/grid_registry_mod.F90)'

    !=======================================================================
    ! Look up a variable; Return metadata and/or a pointer to the data
    !=======================================================================
    CALL Registry_Lookup( am_I_Root   = am_I_Root,           &
                          Registry    = Registry,            &
                          State       = State,               &
                          Variable    = Variable,            &
                          Description = Description,         &
                          Dimensions  = Dimensions,          &
                          KindVal     = KindVal,             &
                          MemoryInKb  = MemoryInKb,          &
                          Rank        = Rank,                &
                          Units       = Units,               &
                          Ptr1d       = Ptr1d,               &
                          Ptr2d       = Ptr2d,               &
                          RC          = RC                  )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable "' // TRIM( Variable ) // &
               '" in the gc_grid_mod.F90 registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_Grid
!EOC
END MODULE Grid_Registry_Mod
