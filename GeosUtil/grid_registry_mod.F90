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
!  24 Aug 2017 - R. Yantosca - Rename variables to be consistent w/ convention
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Module variables
  REAL(f4), ALLOCATABLE, TARGET :: Area(:,:)  ! Surface area 
  REAL(f8), ALLOCATABLE, TARGET :: Time(:  )  ! Time
  REAL(f8), ALLOCATABLE, TARGET :: HyAm(:  )  ! Hybrid Ap at level midpoint
  REAL(f8), ALLOCATABLE, TARGET :: HyBm(:  )  ! Hybrid B  at level midpoint
  REAL(f8), ALLOCATABLE, TARGET :: Lev (:  )  ! Level midpoint coordinate
  REAL(f8), ALLOCATABLE, TARGET :: HyAi(:  )  ! Hybrid Ap at level interface
  REAL(f8), ALLOCATABLE, TARGET :: HyBi(:  )  ! Hybrid B  at level interface
  REAL(f8), ALLOCATABLE, TARGET :: ILev(:  )  ! Level interface coordinate
  REAL(f8), ALLOCATABLE, TARGET :: Lat (:  )  ! Latitude
  REAL(f8), ALLOCATABLE, TARGET :: Lon (:  )  ! Longitude
  REAL(f8),              TARGET :: P0         ! Reference pressure

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
!  25 Aug 2017 - R. Yantosca - Now register the P0 as 0-D REAL*8; AREA as
!                              2-D REAL*4, and all others as 1-D REAL*8
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
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'xy',                             &
                            Data2d_4     = Area,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_AREA', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register P0
    !======================================================================
    Variable = 'P0'
    Desc     = 'reference pressure'
    Units    = 'hPa'

    ! Initialize
    P0       = 1000.0_f8

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = '-',                              &
                            Data0d_8     = P0,                               &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_P0', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register Time
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
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 't',                              &
                            Data1d_8     = Time,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_TIME', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register HyAm
    !======================================================================
    Variable = 'HYAM'
    Desc     = 'hybrid A coefficient at layer midpoints'
    Units    = 'hPa'

    ! Allocate
    ALLOCATE( HyAm( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_HYAM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR
       HyAm(L) = ( Get_Ap( L ) + Get_Ap( L+1 ) ) * 0.5_f8
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            Data1d_8     = HyAm,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_HYAM', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register HyBm
    !======================================================================
    Variable = 'HYBM'
    Desc     = 'hybrid B coefficient at layer midpoints'
    Units    = '1'

    ! Allocate
    ALLOCATE( HyBm( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_HYBM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR
       HyBm(L) = ( Get_Bp( L ) + Get_Bp( L+1 ) ) * 0.5_f8
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            Data1d_8     = HyBm,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_HYBM', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register Lev
    !======================================================================
    Variable = 'LEV'
    Desc     = 'hybrid level at midpoints ((A/P0)+B)'
    Units    = 'level'

    ! Allocate
    ALLOCATE( Lev( LLPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LEV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR
       Lev(L) = ( HyAm(L) / P0 ) + HyBm(L)
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            Data1d_8     = Lev,                              &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_LEV', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register HyAi
    !======================================================================
    Variable = 'HYAI'
    Desc     = 'hybrid A coefficient at layer interfaces'
    Units    = 'hPa'

    ! Allocate
    ALLOCATE( HyAi( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_HYAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR+1
       HyAi(L) = Get_Ap(L)
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         & 
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            OnLevelEdges = .TRUE.,                           &
                            Data1d_8     = HyAi,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_HYAI', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register HyBi
    !======================================================================
    Variable = 'HYBI'
    Desc     = 'hybrid B coefficient at layer interfaces'
    Units    = '1'

    ! Allocate
    ALLOCATE( HyBi( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_HYBI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR+1
       HyBi(L) = Get_Bp(L)
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            OnLevelEdges = .TRUE.,                           &
                            Data1d_8     = HyBi,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_HYBI', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register ILev
    !======================================================================
    Variable = 'ILEV'
    Desc     = 'hybrid level at interfaces ((A/P0)+B)'
    Units    = 'level'

    ! Allocate
    ALLOCATE( ILev( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'GRID_ILEV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO L = 1, LLPAR+1
       ILev(L) = ( HyAi(L) / P0 ) + HyBi(L)
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'z',                              &
                            OnLevelEdges = .TRUE.,                           &
                            Data1d_8     = ILev,                             &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_ILEV', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LAT
    !======================================================================
    Variable = 'LAT'
    Desc     = 'Latitude'
    Units    = 'degrees_north'

    ! Allocate
    ALLOCATE( Lat( JJPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO J = 1, JJPAR
       Lat(J) = Get_YMid( 1, J, 1 )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'y',                              &
                            Data1d_8     = Lat,                              &
                            RC           = RC                               )

    CALL GC_CheckVar( 'GRID_LAT', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !======================================================================
    ! Allocate and register LON
    !======================================================================
    Variable = 'LON'
    Desc     = 'Longitude'
    Units    = 'degrees_east'

    ! Allocate
    ALLOCATE( Lon( IIPAR ), STAT=RC )
    CALL GC_CheckVar( 'GRID_LON', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    DO I = 1, IIPAR
       Lon(I) = Get_XMid( I, 1, 1 )
    ENDDO

    ! Register
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = Registry,                         &
                            State        = State,                            &
                            Variable     = Variable,                         &
                            Description  = Desc,                             &
                            Units        = Units,                            &
                            DimNames     = 'x',                              &
                            Data1d_8     = Lon,                              &
                            RC           = RC                               )

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
!  05 Nov 2018 - R. Yantosca - Now deallocate AND nullify the Registry object
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
       DEALLOCATE( Area, STAT=RC )
       CALL GC_CheckVar( 'GRID_AREA', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( Time ) ) THEN
       DEALLOCATE( Time, STAT=RC )
       CALL GC_CheckVar( 'GRID_TIME', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( Lev ) ) THEN
       DEALLOCATE( Lev, STAT=RC )
       CALL GC_CheckVar( 'GRID_LEV', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( HyAm ) ) THEN
       DEALLOCATE( HyAm, STAT=RC )
       CALL GC_CheckVar( 'GRID_HYAM', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( HyBm ) ) THEN
       DEALLOCATE( HyBm, STAT=RC )
       CALL GC_CheckVar( 'GRID_HYBM', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( ILev ) ) THEN
       DEALLOCATE( ILev, STAT=RC )
       CALL GC_CheckVar( 'GRID_ILEV', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( HyAi ) ) THEN
       DEALLOCATE( HyAi, STAT=RC )
       CALL GC_CheckVar( 'GRID_HYAI', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( HyBi ) ) THEN
       DEALLOCATE( HyBi, STAT=RC )
       CALL GC_CheckVar( 'GRID_HYBI', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( Lat ) ) THEN
       DEALLOCATE( Lat, STAT=RC )
       CALL GC_CheckVar( 'GRID_LAT', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( Lon ) ) THEN
       DEALLOCATE( Lon, STAT=RC )
       CALL GC_CheckVar( 'GRID_LON', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object "Registry"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       Registry => NULL()
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
  SUBROUTINE Lookup_Grid( am_I_Root,  Variable,     RC,         Description, &
                          Dimensions, KindVal,      MemoryInKb, Rank,        &
                          Units,      OnLevelEdges, Ptr0d_8,    Ptr1d_8,     &
                          Ptr2d_4                                           )
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
    LOGICAL,             OPTIONAL    :: OnLevelEdges    ! =T if data is defined
                                                        !  on vertical grid
                                                        !  edges; =F if center

    ! Pointers to data
    REAL(f8),   POINTER, OPTIONAL    :: Ptr0d_8         ! 0D 8-byte data
    REAL(f8),   POINTER, OPTIONAL    :: Ptr1d_8(:  )    ! 1D 8-byte data
    REAL(f4),   POINTER, OPTIONAL    :: Ptr2d_4(:,:)    ! 2D 4-byte data
!
! !REMARKS:
!  We keep the StateName variable private to this module. Users only have
!  to supply the name of each module variable.
!
! !REVISION HISTORY:
!  23 Aug 2017 - R. Yantosca - Initial version
!  24 Aug 2017 - R. Yantosca - Added optional OnLevelEdges argument
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
    CALL Registry_Lookup( am_I_Root    = am_I_Root,                          &
                          Registry     = Registry,                           &
                          State        = State,                              & 
                          Variable     = Variable,                           &
                          Description  = Description,                        &
                          Dimensions   = Dimensions,                         &
                          KindVal      = KindVal,                            &
                          MemoryInKb   = MemoryInKb,                         &
                          Rank         = Rank,                               &
                          Units        = Units,                              &
                          OnLevelEdges = OnLevelEdges,                       &
                          Ptr0d_8      = Ptr0d_8,                            &
                          Ptr1d_8      = Ptr1d_8,                            &
                          Ptr2d_4      = Ptr2d_4,                            &
                          RC           = RC                                 )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable "' // TRIM( Variable ) // &
               '" in the grid_registry_mod.F90 registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_Grid
!EOC
END MODULE Grid_Registry_Mod
