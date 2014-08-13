!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_arr_mod.F90
!
! !DESCRIPTION: Module HCO\_Arr\_Mod contains routines and variables to 
! initialize, validate, and cleanup HEMCO data arrays. HEMCO data arrays
! can be 2D or 3D. They can be organized as single arrays or as vector
! of arrays to represent an additional dimension (time).
!\\
!\\
! The public data types Arr2D\_Hp and Arr3D\_Hp represent the 2D/3D arrays 
! used by HEMCO. The HEMCO precision HP is defined in HCO\_Error\_Mod.
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_Arr_Mod 
!
! !USES:
!
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_ArrInit
  PUBLIC  :: HCO_ValInit
  PUBLIC  :: HCO_ArrAssert
  PUBLIC  :: HCO_ArrCleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: HCO_ArrInit_3D_Hp
  PRIVATE :: HCO_ArrInit_2D_Hp
  PRIVATE :: HCO_ArrInit_2D_I
  PRIVATE :: HCO_ArrVecInit_3D_Hp
  PRIVATE :: HCO_ArrVecInit_2D_Hp
  PRIVATE :: HCO_ValInit_3D_Sp
  PRIVATE :: HCO_ValInit_3D_Dp
  PRIVATE :: HCO_ValInit_2D_Sp
  PRIVATE :: HCO_ValInit_2D_Dp
  PRIVATE :: HCO_ValInit_2D_I
  PRIVATE :: HCO_ArrAssert_2D_Hp
  PRIVATE :: HCO_ArrAssert_3D_Hp
  PRIVATE :: HCO_ArrCleanup_3D_Hp
  PRIVATE :: HCO_ArrCleanup_2D_Hp
  PRIVATE :: HCO_ArrCleanup_2D_I
  PRIVATE :: HCO_ArrVecCleanup_3D_Hp
  PRIVATE :: HCO_ArrVecCleanup_2D_Hp
  PRIVATE :: HCO_ValCleanup_3D_Sp
  PRIVATE :: HCO_ValCleanup_3D_Dp
  PRIVATE :: HCO_ValCleanup_2D_Sp
  PRIVATE :: HCO_ValCleanup_2D_Dp
  PRIVATE :: HCO_ValCleanup_2D_I
!
! !PUBLIC DATA MEMBERS:
!
  ! 2D arrays
  TYPE, PUBLIC :: Arr2D_Hp
     REAL(hp), POINTER :: Val(:,:)    ! x,y
  END TYPE Arr2D_Hp

  TYPE, PUBLIC :: Arr2D_I
     INTEGER,  POINTER :: Val(:,:)    ! x,y
  END TYPE Arr2D_I

  ! 3D arrays
  TYPE, PUBLIC :: Arr3D_Hp
     REAL(hp), POINTER :: Val(:,:,:)  ! x,y,z
  END TYPE Arr3D_Hp
!
! !PRIVATE TYPES:
!
  INTERFACE HCO_ArrInit 
     MODULE PROCEDURE HCO_ArrInit_3D_Hp
     MODULE PROCEDURE HCO_ArrInit_2D_Hp
     MODULE PROCEDURE HCO_ArrInit_2D_I
     MODULE PROCEDURE HCO_ArrVecInit_3D_Hp
     MODULE PROCEDURE HCO_ArrVecInit_2D_Hp
  END INTERFACE HCO_ArrInit

  INTERFACE HCO_ValInit 
     MODULE PROCEDURE HCO_ValInit_3D_Sp
     MODULE PROCEDURE HCO_ValInit_3D_Dp
     MODULE PROCEDURE HCO_ValInit_2D_Sp
     MODULE PROCEDURE HCO_ValInit_2D_Dp
     MODULE PROCEDURE HCO_ValInit_2D_I
  END INTERFACE HCO_ValInit

  INTERFACE HCO_ArrAssert
     MODULE PROCEDURE HCO_ArrAssert_2D_Hp
     MODULE PROCEDURE HCO_ArrAssert_3D_Hp
  END INTERFACE HCO_ArrAssert

  INTERFACE HCO_ArrCleanup
     MODULE PROCEDURE HCO_ArrCleanup_3D_Hp
     MODULE PROCEDURE HCO_ArrCleanup_2D_Hp
     MODULE PROCEDURE HCO_ArrCleanup_2D_I
     MODULE PROCEDURE HCO_ArrVecCleanup_3D_Hp
     MODULE PROCEDURE HCO_ArrVecCleanup_2D_Hp
  END INTERFACE HCO_ArrCleanup

  INTERFACE HCO_ValCleanup
     MODULE PROCEDURE HCO_ValCleanup_3D_Sp
     MODULE PROCEDURE HCO_ValCleanup_3D_Dp
     MODULE PROCEDURE HCO_ValCleanup_2D_Sp
     MODULE PROCEDURE HCO_ValCleanup_2D_Dp
     MODULE PROCEDURE HCO_ValCleanup_2D_I
  END INTERFACE HCO_ValCleanup
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  01 Jul 2014 - R. Yantosca - Corrected errors in ProTeX headers
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrInit_2D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_2D\_Hp initializes the given data
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrInit_2D_Hp( Arr, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Hp), POINTER       :: Arr       ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! HCO_ArrInit_2D_Hp begins here
    ! ================================================================

    IF ( .NOT. ASSOCIATED( Arr) ) ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrInit_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrInit_2D_I
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_2D\_I initializes the given data
! container integer 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrInit_2D_I( Arr, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_I),  POINTER       :: Arr   ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! HCO_ArrInit_2D_I begins here
    ! ================================================================

    IF ( .NOT. ASSOCIATED(Arr) ) ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrInit_2D_I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrInit_3D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_3D\_Hp initializes the given data
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrInit_3D_Hp( Arr, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Hp), POINTER       :: Arr   ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nz        ! z-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! HCO_ArrInit_3D_Hp begins here
    ! ================================================================

    IF ( .NOT. ASSOCIATED(Arr) ) ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, nz, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrInit_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecInit_2D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_2D\_Hp initializes the given data
! container 2D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecInit_2D_Hp( ArrVec, nn, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Hp),   POINTER       :: ArrVec(:) ! Array vector
    INTEGER,          INTENT(IN)    :: nn        ! vector length 
    INTEGER,          INTENT(IN)    :: nx        ! x-dim
    INTEGER,          INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecInit_2D_Hp begins here
    ! ================================================================

    ! Init
    NULLIFY( ArrVec ) 
   
    IF ( nn > 0 ) THEN
       IF ( .NOT. ASSOCIATED(ArrVec) ) ALLOCATE(ArrVec(nn))
       DO I = 1, nn
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, RC )
          IF ( RC/=HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrVecInit_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecInit_3D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_3D\_Hp initializes the given data
! container 3D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecInit_3D_Hp( ArrVec, nn, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Hp),   POINTER       :: ArrVec(:) ! Array vector
    INTEGER,          INTENT(IN)    :: nn        ! vector length 
    INTEGER,          INTENT(IN)    :: nx        ! x-dim
    INTEGER,          INTENT(IN)    :: ny        ! y-dim
    INTEGER,          INTENT(IN)    :: nz        ! z-dim
!
! !INPUT/OUTPUT PARAMETERS:
!          
    INTEGER,          INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecInit_3D_Hp begins here
    ! ================================================================

    ! Init
    NULLIFY( ArrVec ) 
  
    IF ( nn > 0 ) THEN 
       IF ( .NOT. ASSOCIATED(ArrVec) ) ALLOCATE(ArrVec(nn))
       DO I = 1, nn
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, nz, RC )
          IF ( RC/=HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrVecInit_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_2D\_Sp initializes the given data
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_Sp( Val, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    REAL(sp),       POINTER       :: Val(:,:)  ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! ================================================================
    ! HCO_ValInit_2D_Sp begins here
    ! ================================================================

    Val => NULL()
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          CALL HCO_ERROR ( 'Arr2D value allocation error', RC )
          RETURN
       ENDIF
       Val(:,:) = 0.0_sp
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValInit_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_2D_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_2D\_Dp initializes the given data
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_Dp( Val, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    REAL(dp),       POINTER       :: Val(:,:)  ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! ================================================================
    ! HCO_ValInit_2D_Dp begins here
    ! ================================================================

    Val => NULL()
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          CALL HCO_ERROR ( 'Arr2D value allocation error', RC )
          RETURN
       ENDIF
       Val(:,:) = 0.0_dp
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValInit_2D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_2D_I
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_2D\_I initializes the given data
! container integer 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_I( Val, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,        POINTER       :: Val(:,:)  ! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! ================================================================
    ! HCO_ValInit_2D_I begins here
    ! ================================================================

    Val => NULL()
    IF ( nx > 0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          CALL HCO_ERROR ( 'Arr2D value allocation error', RC )
          RETURN
       ENDIF
       Val(:,:) = 0 
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValInit_2D_I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_3D_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_3D\_Dp initializes the given data
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_3D_Dp( Val, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    REAL(dp),       POINTER       :: Val(:,:,:)! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nz        ! z-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! ================================================================
    ! HCO_ValInit_3D_Dp begins here
    ! ================================================================

    Val => NULL()
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny,nz),STAT=AS)
       IF(AS/=0) THEN
          CALL HCO_ERROR ( 'Arr3D value allocation error', RC )
          RETURN
       ENDIF
       Val(:,:,:) = 0.0_dp
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValInit_3D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_3D\_Sp initializes the given data
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_3D_Sp( Val, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    REAL(sp),       POINTER       :: Val(:,:,:)! Array 
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
    INTEGER,        INTENT(IN)    :: nz        ! z-dim
!
! !INPUT/OUTPUT PARAMETERS:
!          
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! ================================================================
    ! HCO_ValInit_3D_Sp begins here
    ! ================================================================

    Val => NULL()
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny,nz),STAT=AS)
       IF(AS/=0) THEN
          CALL HCO_ERROR ( 'Arr3D value allocation error', RC )
          RETURN
       ENDIF
       Val(:,:,:) = 0.0_sp
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ValInit_3D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrAssert_3D_Hp 
!
! !DESCRIPTION: Routine HCO\_ArrAssert\_3D\_Hp makes sure that the passed 
! 3D array is allocated. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrAssert_3D_Hp( ThisArr3D, I, J, L, RC )
!
! !INPUT PARAMETERS:
! 
    TYPE(Arr3D_Hp),  POINTER         :: ThisArr3D ! 3D array
    INTEGER,         INTENT(IN   )   :: I, J, L   ! Array dims 
!
! !INPUT/OUTPUT PARAMETERS:
! 

    INTEGER,         INTENT(INOUT)   :: RC        ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_ArrAssert_3D_Hp begins here!
    !=====================================================================
  
    ! Check flux array
    IF ( .NOT. ASSOCIATED ( ThisArr3D ) ) THEN
       CALL HCO_ArrInit( ThisArr3D, I, J, L, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ELSEIF ( .NOT. ASSOCIATED ( ThisArr3D%Val ) ) THEN
       CALL HCO_ValInit ( ThisArr3D%Val, I, J, L, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrAssert_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrAssert_2D_Hp 
!
! !DESCRIPTION: Routine HCO\_ArrAssert\_2D\_Hp makes sure that the passed 
! 2D array is allocated. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrAssert_2D_Hp( ThisArr2D, I, J, RC )
!
! !INPUT PARAMETERS:
! 
    TYPE(Arr2D_Hp),  POINTER         :: ThisArr2D ! 2D array
    INTEGER,         INTENT(IN   )   :: I, J      ! Array dims 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)   :: RC        ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_ArrAssert_2D_Hp begins here!
    !=====================================================================
  
    ! Check flux array
    IF ( .NOT. ASSOCIATED ( ThisArr2D ) ) THEN
       CALL HCO_ArrInit( ThisArr2D, I, J, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ELSEIF ( .NOT. ASSOCIATED ( ThisArr2D%Val ) ) THEN
       CALL HCO_ValInit ( ThisArr2D%Val, I, J, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ENDIF
  
    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrAssert_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrCleanup_2D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrCleanup\_2D\_Hp cleans up the given 
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrCleanup_2D_Hp( Arr, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Hp),      POINTER  :: Arr       ! Array 
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC

    ! ================================================================
    ! HCO_ArrCleanup_2D_Hp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(Arr) ) THEN 
       CALL HCO_ValCleanup( Arr%Val, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrCleanup_2D_I
!
! !DESCRIPTION: Subroutine HCO\_ArrCleanup\_2D\_I cleans up the given 
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrCleanup_2D_I( Arr, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_I),       POINTER  :: Arr       ! Array 
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC

    ! ================================================================
    ! HCO_ArrCleanup_2D_I begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(Arr) ) THEN 
       CALL HCO_ValCleanup( Arr%Val, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_2D_I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrCleanup_3D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrCleanup\_3D\_Hp cleans up the given 
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrCleanup_3D_Hp( Arr, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Hp),      POINTER  :: Arr       ! Array 
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC

    ! ================================================================
    ! HCO_ArrCleanup_3D_Hp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(Arr) ) THEN 
       CALL HCO_ValCleanup( Arr%Val, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecCleanup_2D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecCleanup\_2D\_Hp cleans up the given 
! container 2D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecCleanup_2D_Hp( ArrVec, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Hp),      POINTER  :: ArrVec(:) ! Array 
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecCleanup_2D_Hp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(ArrVec) ) THEN 
       DO I = 1, SIZE(ArrVec,1)
          CALL HCO_ValCleanup( ArrVec(I)%Val, DC )
       ENDDO

       DEALLOCATE ( ArrVec )
    ENDIF

  END SUBROUTINE HCO_ArrVecCleanup_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecCleanup_3D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecCleanup\_3D\_Hp cleans up the given 
! container 3D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecCleanup_3D_Hp( ArrVec, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Hp),      POINTER  :: ArrVec(:) ! Array 
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecCleanup_3D_Hp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(ArrVec) ) THEN 
       DO I = 1, SIZE(ArrVec,1)
          CALL HCO_ValCleanup( ArrVec(I)%Val, DC )
       ENDDO

       DEALLOCATE ( ArrVec )
    ENDIF

  END SUBROUTINE HCO_ArrVecCleanup_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_2D_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_2D\_Dp cleans up the given 
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_Dp( Val, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    REAL(dp),            POINTER  :: Val(:,:)  ! Array 
    LOGICAL, INTENT(IN)           :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) ) THEN 
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_2D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_2D\_Sp cleans up the given 
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_Sp( Val, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    REAL(sp), POINTER    :: Val(:,:)  ! Array 
    LOGICAL,  INTENT(IN) :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) ) THEN 
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_2D_I
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_2D\_I cleans up the given 
! container 2D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_I( Val, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, POINTER    :: Val(:,:)  ! Array 
    LOGICAL, INTENT(IN) :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    IF ( DeepClean .AND. ASSOCIATED(Val) ) THEN 
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_2D_I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_3D_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_3D\_Dp cleans up the given 
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_3D_Dp( Val, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    REAL(dp), POINTER    :: Val(:,:,:) ! Array 
    LOGICAL,  INTENT(IN) :: DeepClean  ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) ) THEN 
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_3D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_3D\_Sp cleans up the given 
! container 3D array. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_3D_Sp( Val, DeepClean ) 
!
! !INPUT PARAMETERS:
!
    REAL(sp), POINTER    :: Val(:,:,:) ! Array 
    LOGICAL,  INTENT(IN) :: DeepClean  ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) ) THEN 
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_3D_Sp
!EOC
END MODULE HCO_Arr_Mod
