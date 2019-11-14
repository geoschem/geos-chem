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
! used by HEMCO. The HEMCO precision HP is defined in HCO\_Error\_Mod. There
! is also an integer 2D array (type Arr2D\_I) that can be used to store
! integer data. Additional data types can be added as needed.
!\\
!\\
! Data values are stored in array 'Val'. Val can be self-allocated or a
! pointer to existing data, as denoted by the Alloc flag.
!\\
!\\
! At the moment, the following HEMCO structures use HEMCO arrays:
!\begin{itemize}
!\item FileData: emission data (base emissions, scale factors, masks) stored
!  in the FileData derived type. These data are read from disk as specified in
!  the configuration file. See HCO\_FileData\_Mod.F90.
!\item FluxArr: the HEMCO flux arrays (emissions and deposition velocities)
!  stored in the HEMCO state object. See HCO\_State\_Mod.F90.
!\item Grid: all grid information arrays (x midpoints, y midpoints, etc.)
!  stored in the HEMCO state object.
!\item ExtDat: external data required by the extensions (primarily met fields).
!  See HCOX\_State\_Mod.F90.
!\end{itemize}
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
  PUBLIC  :: HCO_ArrAssert
  PUBLIC  :: HCO_ArrCleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: HCO_ArrInit_3D_Hp
  PRIVATE :: HCO_ArrInit_2D_Hp
  PRIVATE :: HCO_ArrInit_3D_Sp
  PRIVATE :: HCO_ArrInit_2D_Sp
  PRIVATE :: HCO_ArrInit_2D_I
  PRIVATE :: HCO_ArrVecInit_3D_Hp
  PRIVATE :: HCO_ArrVecInit_2D_Hp
  PRIVATE :: HCO_ArrVecInit_3D_Sp
  PRIVATE :: HCO_ArrVecInit_2D_Sp
  PRIVATE :: HCO_ValInit
  PRIVATE :: HCO_ValInit_3D_Sp
  PRIVATE :: HCO_ValInit_3D_Dp
  PRIVATE :: HCO_ValInit_2D_Sp
  PRIVATE :: HCO_ValInit_2D_Dp
  PRIVATE :: HCO_ValInit_2D_I
  PRIVATE :: HCO_ArrAssert_2D_Hp
  PRIVATE :: HCO_ArrAssert_3D_Hp
  PRIVATE :: HCO_ArrAssert_2D_Sp
  PRIVATE :: HCO_ArrAssert_3D_Sp
  PRIVATE :: HCO_ArrCleanup_3D_Hp
  PRIVATE :: HCO_ArrCleanup_2D_Hp
  PRIVATE :: HCO_ArrCleanup_3D_Sp
  PRIVATE :: HCO_ArrCleanup_2D_Sp
  PRIVATE :: HCO_ArrCleanup_2D_I
  PRIVATE :: HCO_ArrVecCleanup_3D_Hp
  PRIVATE :: HCO_ArrVecCleanup_2D_Hp
  PRIVATE :: HCO_ArrVecCleanup_3D_Sp
  PRIVATE :: HCO_ArrVecCleanup_2D_Sp
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
     LOGICAL           :: Alloc       ! Allocated?
  END TYPE Arr2D_Hp

  TYPE, PUBLIC :: Arr2D_I
     INTEGER,  POINTER :: Val(:,:)    ! x,y
     LOGICAL           :: Alloc       ! Allocated?
  END TYPE Arr2D_I

  TYPE, PUBLIC :: Arr2D_Sp
     REAL(sp), POINTER :: Val(:,:)    ! x,y
     LOGICAL           :: Alloc       ! Allocated?
  END TYPE Arr2D_Sp

  ! 3D arrays
  TYPE, PUBLIC :: Arr3D_Hp
     REAL(hp), POINTER :: Val(:,:,:)  ! x,y,z
     LOGICAL           :: Alloc       ! Allocated?
  END TYPE Arr3D_Hp

  TYPE, PUBLIC :: Arr3D_Sp
     REAL(sp), POINTER :: Val(:,:,:)  ! x,y,z
     LOGICAL           :: Alloc       ! Allocated?
  END TYPE Arr3D_Sp
!
! !PRIVATE TYPES:
!
  INTERFACE HCO_ArrInit
     MODULE PROCEDURE HCO_ArrInit_3D_Hp
     MODULE PROCEDURE HCO_ArrInit_2D_Hp
     MODULE PROCEDURE HCO_ArrInit_3D_Sp
     MODULE PROCEDURE HCO_ArrInit_2D_Sp
     MODULE PROCEDURE HCO_ArrInit_2D_I
     MODULE PROCEDURE HCO_ArrVecInit_3D_Hp
     MODULE PROCEDURE HCO_ArrVecInit_2D_Hp
     MODULE PROCEDURE HCO_ArrVecInit_3D_Sp
     MODULE PROCEDURE HCO_ArrVecInit_2D_Sp
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
     MODULE PROCEDURE HCO_ArrAssert_2D_Sp
     MODULE PROCEDURE HCO_ArrAssert_3D_Sp
     MODULE PROCEDURE HCO_ArrAssert_2D_I
  END INTERFACE HCO_ArrAssert

  INTERFACE HCO_ArrCleanup
     MODULE PROCEDURE HCO_ArrCleanup_3D_Hp
     MODULE PROCEDURE HCO_ArrCleanup_2D_Hp
     MODULE PROCEDURE HCO_ArrCleanup_2D_I
     MODULE PROCEDURE HCO_ArrCleanup_3D_Sp
     MODULE PROCEDURE HCO_ArrCleanup_2D_Sp
     MODULE PROCEDURE HCO_ArrVecCleanup_3D_Hp
     MODULE PROCEDURE HCO_ArrVecCleanup_2D_Hp
     MODULE PROCEDURE HCO_ArrVecCleanup_3D_Sp
     MODULE PROCEDURE HCO_ArrVecCleanup_2D_Sp
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
!  01 Oct 2014 - C. Keller   - Added Alloc flag
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
! container 2D array. nx and ny denote the array size dimensions. If nx is
! set to 0, no data is allocated but Val is set to a (nullified) pointer
! instead.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! HCO_ArrInit_2D_Hp begins here
    ! ================================================================

    NULLIFY (Arr)
    ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, Arr%Alloc, RC )
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
! !IROUTINE: HCO_ArrInit_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_2D\_Sp initializes the given data
! container 2D array. nx and ny denote the array size dimensions. If nx is
! set to 0, no data is allocated but Val is set to a (nullified) pointer
! instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrInit_2D_Sp( Arr, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Sp), POINTER       :: Arr       ! Array
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! HCO_ArrInit_2D_Sp begins here
    ! ================================================================

    NULLIFY (Arr)
    ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, Arr%Alloc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrInit_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrInit_2D_I
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_2D\_I initializes the given data
! container integer 2D array. nx and ny denote the array size dimensions.
! If nx is set to 0, no data is allocated but Val is set to a (nullified)
! pointer instead.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! HCO_ArrInit_2D_I begins here
    ! ================================================================

    NULLIFY (Arr)
    ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, Arr%Alloc, RC )
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
! container 3D array. nx and ny denote the array size dimensions. If nx is
! set to 0, no data is allocated but Val is set to a (nullified) pointer
! instead.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! HCO_ArrInit_3D_Hp begins here
    ! ================================================================

    NULLIFY (Arr)
    ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, nz, Arr%Alloc, RC )
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
! !IROUTINE: HCO_ArrInit_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrInit\_3D\_Sp initializes the given data
! container 3D array. nx and ny denote the array size dimensions. If nx is
! set to 0, no data is allocated but Val is set to a (nullified) pointer
! instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrInit_3D_Sp( Arr, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Sp), POINTER       :: Arr   ! Array
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! HCO_ArrInit_3D_Hp begins here
    ! ================================================================

    NULLIFY (Arr)
    ALLOCATE(Arr)
    CALL HCO_ValInit( Arr%Val, nx, ny, nz, Arr%Alloc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrInit_3D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecInit_2D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_2D\_Hp initializes the given data
! container 2D array vector. nn denotes the number of 2D arrays, and nx and ny
! denote the array size dimensions. If nx is set to 0, no data is allocated but
! Val is set to a (nullified) pointer instead.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
    NULLIFY(ArrVec)

    IF ( nn > 0 ) THEN
       IF ( .NOT. ASSOCIATED(ArrVec) ) ALLOCATE(ArrVec(nn))
       DO I = 1, nn
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, ArrVec(I)%Alloc, RC )
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
! !IROUTINE: HCO_ArrVecInit_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_2D\_Sp initializes the given data
! container 2D array vector. nn denotes the number of 2D arrays, and nx and ny
! denote the array size dimensions. If nx is set to 0, no data is allocated but
! Val is set to a (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecInit_2D_Sp( ArrVec, nn, nx, ny, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Sp),   POINTER       :: ArrVec(:) ! Array vector
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecInit_2D_Sp begins here
    ! ================================================================

    ! Init
    NULLIFY(ArrVec)

    IF ( nn > 0 ) THEN
       IF ( .NOT. ASSOCIATED(ArrVec) ) ALLOCATE(ArrVec(nn))
       DO I = 1, nn
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, ArrVec(I)%Alloc, RC )
          IF ( RC/=HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrVecInit_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrVecInit_3D_Hp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_3D\_Hp initializes the given data
! container 3D array vector. nn denotes the number of 2D arrays, and nx and ny
! denote the array size dimensions. If nx is set to 0, no data is allocated but
! Val is set to a (nullified) pointer instead.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, nz, ArrVec(I)%Alloc, RC )
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
! !IROUTINE: HCO_ArrVecInit_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecInit\_3D\_Sp initializes the given data
! container 3D array vector. nn denotes the number of 2D arrays, and nx and ny
! denote the array size dimensions. If nx is set to 0, no data is allocated but
! Val is set to a (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecInit_3D_Sp( ArrVec, nn, nx, ny, nz, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Sp),   POINTER       :: ArrVec(:) ! Array vector
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecInit_3D_Sp begins here
    ! ================================================================

    ! Init
    NULLIFY( ArrVec )

    IF ( nn > 0 ) THEN
       IF ( .NOT. ASSOCIATED(ArrVec) ) ALLOCATE(ArrVec(nn))
       DO I = 1, nn
          CALL HCO_ValInit( ArrVec(I)%Val, nx, ny, nz, ArrVec(I)%Alloc, RC )
          IF ( RC/=HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrVecInit_3D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValInit_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ValInit\_2D\_Sp initializes the given data
! container 2D single precision array. nx and ny denote the array size
! dimensions. If nx is set to 0, no data is allocated but Val is set to a
! (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_Sp( Val, nx, ny, Alloc, RC )
!
! !INPUT PARAMETERS:
!
    REAL(sp),       POINTER       :: Val(:,:)  ! Array
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(  OUT) :: Alloc     ! allocated?
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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

    Val   => NULL()
    ALLOC =  .FALSE.

    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          WRITE(*,*) 'Arr2D value allocation error'
          RC = HCO_FAIL
          RETURN
       ENDIF
       Val(:,:) = 0.0_sp
       ALLOC    = .TRUE.
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
! container 2D double precision array. nx and ny denote the array size
! dimensions. If nx is set to 0, no data is allocated but Val is set to a
! (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_Dp( Val, nx, ny, Alloc, RC )
!
! !INPUT PARAMETERS:
!
    REAL(dp),       POINTER       :: Val(:,:)  ! Array
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(  OUT) :: Alloc     ! allocated?
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
    Alloc = .FALSE.
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          WRITE(*,*) 'Arr2D value allocation error'
          RC = HCO_FAIL
          RETURN
       ENDIF
       Val(:,:) = 0.0_dp
       Alloc = .TRUE.
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
! container 2D integer array. nx and ny denote the array size
! dimensions. If nx is set to 0, no data is allocated but Val is set to a
! (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_2D_I( Val, nx, ny, Alloc, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,        POINTER       :: Val(:,:)  ! Array
    INTEGER,        INTENT(IN)    :: nx        ! x-dim
    INTEGER,        INTENT(IN)    :: ny        ! y-dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(  OUT) :: Alloc     ! allocated?
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
    Alloc = .FALSE.
    IF ( nx > 0 ) THEN
       ALLOCATE(Val(nx,ny),STAT=AS)
       IF(AS/=0) THEN
          WRITE(*,*) 'Arr2D value allocation error'
          RC = HCO_FAIL
          RETURN
       ENDIF
       Val(:,:) = 0
       Alloc = .TRUE.
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
! container 3D double precision array. nx and ny denote the array size
! dimensions. If nx is set to 0, no data is allocated but Val is set to a
! (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_3D_Dp( Val, nx, ny, nz, Alloc, RC )
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
    LOGICAL,        INTENT(  OUT) :: Alloc     ! allocated?
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
    Alloc = .FALSE.
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny,nz),STAT=AS)
       IF(AS/=0) THEN
          WRITE(*,*) 'Arr3D value allocation error'
          RC = HCO_FAIL
          RETURN
       ENDIF
       Val(:,:,:) = 0.0_dp
       Alloc = .TRUE.
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
! container 3D single precision array. nx and ny denote the array size
! dimensions. If nx is set to 0, no data is allocated but Val is set to a
! (nullified) pointer instead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValInit_3D_Sp( Val, nx, ny, nz, Alloc, RC )
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
    LOGICAL,        INTENT(  OUT) :: Alloc     ! allocated?
    INTEGER,        INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
    Alloc = .FALSE.
    IF ( nx>0 ) THEN
       ALLOCATE(Val(nx,ny,nz),STAT=AS)
       IF(AS/=0) THEN
          WRITE(*,*) 'Arr3D value allocation error'
          RC = HCO_FAIL
          RETURN
       ENDIF
       Val(:,:,:) = 0.0_sp
       Alloc = .TRUE.
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
       CALL HCO_ValInit ( ThisArr3D%Val, I, J, L, ThisArr3D%Alloc, RC )
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
! !IROUTINE: HCO_ArrAssert_3D_Sp
!
! !DESCRIPTION: Routine HCO\_ArrAssert\_3D\_Sp makes sure that the passed
! 3D array is allocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrAssert_3D_Sp( ThisArr3D, I, J, L, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Sp),  POINTER         :: ThisArr3D ! 3D array
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_ArrAssert_3D_Sp begins here!
    !=====================================================================

    ! Check flux array
    IF ( .NOT. ASSOCIATED ( ThisArr3D ) ) THEN
       CALL HCO_ArrInit( ThisArr3D, I, J, L, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ELSEIF ( .NOT. ASSOCIATED ( ThisArr3D%Val ) ) THEN
       CALL HCO_ValInit ( ThisArr3D%Val, I, J, L, ThisArr3D%Alloc, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrAssert_3D_Sp
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
       CALL HCO_ValInit ( ThisArr2D%Val, I, J, ThisArr2D%Alloc, RC )
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
! !IROUTINE: HCO_ArrAssert_2D_Sp
!
! !DESCRIPTION: Routine HCO\_ArrAssert\_2D\_Sp makes sure that the passed
! 2D array is allocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrAssert_2D_Sp( ThisArr2D, I, J, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Sp),  POINTER         :: ThisArr2D ! 2D array
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_ArrAssert_2D_Sp begins here!
    !=====================================================================

    ! Check flux array
    IF ( .NOT. ASSOCIATED ( ThisArr2D ) ) THEN
       CALL HCO_ArrInit( ThisArr2D, I, J, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ELSEIF ( .NOT. ASSOCIATED ( ThisArr2D%Val ) ) THEN
       CALL HCO_ValInit ( ThisArr2D%Val, I, J, ThisArr2D%Alloc, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrAssert_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrAssert_2D_I
!
! !DESCRIPTION: Routine HCO\_ArrAssert\_2D\_I makes sure that the passed
! 2D array is allocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrAssert_2D_I( ThisArr2D, I, J, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_I),   POINTER         :: ThisArr2D ! 2D array
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_ArrAssert_2D_I begins here!
    !=====================================================================

    ! Check flux array
    IF ( .NOT. ASSOCIATED ( ThisArr2D ) ) THEN
       CALL HCO_ArrInit( ThisArr2D, I, J, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ELSEIF ( .NOT. ASSOCIATED ( ThisArr2D%Val ) ) THEN
       CALL HCO_ValInit ( ThisArr2D%Val, I, J, ThisArr2D%Alloc, RC )
       IF ( RC/= HCO_SUCCESS ) RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ArrAssert_2D_I
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
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate allocated array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
       CALL HCO_ValCleanup( Arr%Val, Arr%Alloc, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_2D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrCleanup_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrCleanup\_2D\_Sp cleans up the given
! container 2D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrCleanup_2D_Sp( Arr, DeepClean )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Sp),      POINTER  :: Arr       ! Array
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate allocated array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC

    ! ================================================================
    ! HCO_ArrCleanup_2D_Sp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(Arr) ) THEN
       CALL HCO_ValCleanup( Arr%Val, Arr%Alloc, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_2D_Sp
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
       CALL HCO_ValCleanup( Arr%Val, Arr%Alloc, DC )
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
       CALL HCO_ValCleanup( Arr%Val, Arr%Alloc, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_3D_Hp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ArrCleanup_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrCleanup\_3D\_Sp cleans up the given
! container 3D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrCleanup_3D_Sp( Arr, DeepClean )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Sp),      POINTER  :: Arr       ! Array
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC

    ! ================================================================
    ! HCO_ArrCleanup_3D_Sp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(Arr) ) THEN
       CALL HCO_ValCleanup( Arr%Val, Arr%Alloc, DC )
       DEALLOCATE ( Arr )
    ENDIF

  END SUBROUTINE HCO_ArrCleanup_3D_Sp
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
          CALL HCO_ValCleanup( ArrVec(I)%Val, ArrVec(I)%Alloc, DC )
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
! !IROUTINE: HCO_ArrVecCleanup_2D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecCleanup\_2D\_Sp cleans up the given
! container 2D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecCleanup_2D_Sp( ArrVec, DeepClean )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr2D_Sp),      POINTER  :: ArrVec(:) ! Array
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecCleanup_2D_Sp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(ArrVec) ) THEN
       DO I = 1, SIZE(ArrVec,1)
          CALL HCO_ValCleanup( ArrVec(I)%Val, ArrVec(I)%Alloc, DC )
       ENDDO

       DEALLOCATE ( ArrVec )
    ENDIF

  END SUBROUTINE HCO_ArrVecCleanup_2D_Sp
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
!  01 Oct 2014 - C. Keller - Added Alloc flag
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
          CALL HCO_ValCleanup( ArrVec(I)%Val, ArrVec(I)%Alloc, DC )
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
! !IROUTINE: HCO_ArrVecCleanup_3D_Sp
!
! !DESCRIPTION: Subroutine HCO\_ArrVecCleanup\_3D\_Sp cleans up the given
! container 3D array vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ArrVecCleanup_3D_Sp( ArrVec, DeepClean )
!
! !INPUT PARAMETERS:
!
    TYPE(Arr3D_Sp),      POINTER  :: ArrVec(:) ! Array
    LOGICAL, INTENT(IN), OPTIONAL :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: DC
    INTEGER :: I

    ! ================================================================
    ! HCO_ArrVecCleanup_3D_Sp begins here
    ! ================================================================

    IF ( PRESENT(DeepClean) ) THEN
       DC = DeepClean
    ELSE
       DC = .TRUE.
    ENDIF

    IF ( ASSOCIATED(ArrVec) ) THEN
       DO I = 1, SIZE(ArrVec,1)
          CALL HCO_ValCleanup( ArrVec(I)%Val, ArrVec(I)%Alloc, DC )
       ENDDO

       DEALLOCATE ( ArrVec )
    ENDIF

  END SUBROUTINE HCO_ArrVecCleanup_3D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ValCleanup_2D_Dp
!
! !DESCRIPTION: Subroutine HCO\_ValCleanup\_2D\_Dp cleans up the given
! container 2D array. If DeepClean is set to TRUE and the array is
! indeed allocated (as determined by the Alloc flag), the array becomes
! deallocated. Otherwise, it is just nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_Dp( Val, Alloc, DeepClean )
!
! !INPUT PARAMETERS:
!
    REAL(dp),            POINTER  :: Val(:,:)  ! Array
    LOGICAL, INTENT(IN)           :: Alloc     ! Allocated?
    LOGICAL, INTENT(IN)           :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) .AND. Alloc ) THEN
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
! container 2D array. If DeepClean is set to TRUE and the array is
! indeed allocated (as determined by the Alloc flag), the array becomes
! deallocated. Otherwise, it is just nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_Sp( Val, Alloc, DeepClean )
!
! !INPUT PARAMETERS:
!
    REAL(sp), POINTER    :: Val(:,:)  ! Array
    LOGICAL,  INTENT(IN) :: Alloc     ! Allocated?
    LOGICAL,  INTENT(IN) :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) .AND. Alloc ) THEN
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
! container 2D array. If DeepClean is set to TRUE and the array is
! indeed allocated (as determined by the Alloc flag), the array becomes
! deallocated. Otherwise, it is just nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_2D_I( Val, Alloc, DeepClean )
!
! !INPUT PARAMETERS:
!
    INTEGER, POINTER    :: Val(:,:)  ! Array
    LOGICAL, INTENT(IN) :: Alloc     ! Allocated?
    LOGICAL, INTENT(IN) :: DeepClean ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    IF ( DeepClean .AND. ASSOCIATED(Val) .AND. Alloc ) THEN
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
! container 3D array. If DeepClean is set to TRUE and the array is
! indeed allocated (as determined by the Alloc flag), the array becomes
! deallocated. Otherwise, it is just nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_3D_Dp( Val, Alloc, DeepClean )
!
! !INPUT PARAMETERS:
!
    REAL(dp), POINTER    :: Val(:,:,:) ! Array
    LOGICAL,  INTENT(IN) :: Alloc      ! Allocated?
    LOGICAL,  INTENT(IN) :: DeepClean  ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) .AND. Alloc ) THEN
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
! container 3D array. If DeepClean is set to TRUE and the array is
! indeed allocated (as determined by the Alloc flag), the array becomes
! deallocated. Otherwise, it is just nullified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ValCleanup_3D_Sp( Val, Alloc, DeepClean )
!
! !INPUT PARAMETERS:
!
    REAL(sp), POINTER    :: Val(:,:,:) ! Array
    LOGICAL,  INTENT(IN) :: Alloc      ! Allocated?
    LOGICAL,  INTENT(IN) :: DeepClean  ! Deallocate array?
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  01 Oct 2014 - C. Keller - Added Alloc flag
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( DeepClean .AND. ASSOCIATED(Val) .AND. Alloc ) THEN
       DEALLOCATE( Val )
    ENDIF
    Val => NULL()

  END SUBROUTINE HCO_ValCleanup_3D_Sp
!EOC
END MODULE HCO_Arr_Mod
