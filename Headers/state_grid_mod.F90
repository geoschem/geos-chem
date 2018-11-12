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
  PUBLIC :: Cleanup_State_Grid
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Grid State
  !=========================================================================
  TYPE, PUBLIC :: GrdState

     !----------------------------------------
     ! GRID MENU fields
     !----------------------------------------
     CHARACTER(LEN=255)          :: GridRes
     REAL(fp)                    :: dx
     REAL(fp)                    :: dy
     REAL(fp)                    :: MinLon
     REAL(fp)                    :: MaxLon
     REAL(fp)                    :: MinLat
     REAL(fp)                    :: MaxLat
     INTEGER                     :: nx
     INTEGER                     :: ny
     INTEGER                     :: nz
     LOGICAL                     :: Its_A_Nested_Grid
     CHARACTER(LEN=255)          :: Res_Dir
     INTEGER                     :: Nested_I0
     INTEGER                     :: Nested_J0
 
     !----------------------------------------
     ! NESTED GRID MENU fields
     !----------------------------------------
     LOGICAL                     :: LWINDO
     LOGICAL                     :: LWINDO2x25
     LOGICAL                     :: LWINDO_NA
     CHARACTER(LEN=255)          :: TPBC_DIR_NA
     LOGICAL                     :: LWINDO_EU
     CHARACTER(LEN=255)          :: TPBC_DIR_EU
     LOGICAL                     :: LWINDO_CH
     CHARACTER(LEN=255)          :: TPBC_DIR_CH
     LOGICAL                     :: LWINDO_AS
     CHARACTER(LEN=255)          :: TPBC_DIR_AS
     LOGICAL                     :: LWINDO_CU
     CHARACTER(LEN=255)          :: TPBC_DIR
     INTEGER                     :: NESTED_TS
     INTEGER                     :: NESTED_I1
     INTEGER                     :: NESTED_J1
     INTEGER                     :: NESTED_I2
     INTEGER                     :: NESTED_J2
     INTEGER                     :: NESTED_I0W
     INTEGER                     :: NESTED_J0W
     INTEGER                     :: NESTED_I0E
     INTEGER                     :: NESTED_J0E

  END TYPE GrdState
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  11 Nov 2018 - M. Sulprizio- Initial version
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
  SUBROUTINE Init_State_Grid( am_I_Root, State_Grid, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  11 Nov 2018 - M. Sulprizio- Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Assume success
    RC            =  GC_SUCCESS

    !----------------------------------------
    ! GRID MENU fields
    !----------------------------------------
    Input_Opt%GridRes                = ''
    Input_Opt%dx                     = 0e+0_fp
    Input_Opt%dy                     = 0e+0_fp
    Input_Opt%MinLon                 = 0e+0_fp
    Input_Opt%MaxLon                 = 0e+0_fp
    Input_Opt%MinLat                 = 0e+0_fp
    Input_Opt%MaxLat                 = 0e+0_fp
    Input_Opt%nx                     = 0
    Input_Opt%ny                     = 0
    Input_Opt%nz                     = 0
    Input_Opt%Its_A_Nested_Grid      = .FALSE.
    Input_Opt%Res_Dir                = './'

    Input_Opt%NESTED_I0              = 0
    Input_Opt%NESTED_J0              = 0

    !----------------------------------------
    ! NESTED GRID MENU fields
    !----------------------------------------
    Input_Opt%LWINDO                 = .FALSE.
    Input_Opt%LWINDO2x25             = .FALSE.
    Input_Opt%LWINDO_NA              = .FALSE.
    Input_Opt%TPBC_DIR_NA            = './'
    Input_Opt%LWINDO_EU              = .FALSE.
    Input_Opt%TPBC_DIR_EU            = './'
    Input_Opt%LWINDO_CH              = .FALSE.
    Input_Opt%TPBC_DIR_CH            = './'
    Input_Opt%LWINDO_AS              = .FALSE.
    Input_Opt%TPBC_DIR_AS            = './'
    Input_Opt%LWINDO_CU              = .FALSE.
    Input_Opt%TPBC_DIR               = './'
    Input_Opt%NESTED_TS              = 0
    Input_Opt%NESTED_I1              = 0
    Input_Opt%NESTED_J1              = 0
    Input_Opt%NESTED_I2              = 0
    Input_Opt%NESTED_J2              = 0
    Input_Opt%NESTED_I0W             = 0
    Input_Opt%NESTED_J0W             = 0 
    Input_Opt%NESTED_I0E             = 0
    Input_Opt%NESTED_J0E             = 0 

   END SUBROUTINE Init_State_Grid
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
  SUBROUTINE Cleanup_State_Grid( am_I_Root, State_Grid, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(INOUT) :: State_Grid   ! Obj for grid state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY: 
!  11 Nov 2018 - M. Sulprizio- Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC      = GC_SUCCESS

  END SUBROUTINE Cleanup_State_Grid
!EOC
END MODULE State_Grid_Mod
