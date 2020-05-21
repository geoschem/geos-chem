!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_scale_mod.F90
!
! !DESCRIPTION: Module hco\_scale\_mod contains a collection of
! routines to uniformly scale emissions by species-specific scale
! scale factors.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Scale_Mod
!
! !USES:
!
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_ScaleInit
  PUBLIC :: HCO_ScaleGet
  PUBLIC :: HCO_ScaleArr
  PUBLIC :: HCO_ScaleFinal
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: HCO_ScaleArr3D_sp
  PRIVATE :: HCO_ScaleArr3D_dp
  PRIVATE :: HCO_ScaleArr2D_sp
  PRIVATE :: HCO_ScaleArr2D_dp
  PRIVATE :: HCO_ScaleArr1D_sp
  PRIVATE :: HCO_ScaleArr1D_dp
!
! !PRIVATE VARIABLES:
!
  REAL(hp), ALLOCATABLE :: SpcScal(:)
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller   - Initial version
!EOP
!
! !INTERFACES:
!
  INTERFACE HCO_ScaleArr
     MODULE PROCEDURE HCO_ScaleArr3D_dp
     MODULE PROCEDURE HCO_ScaleArr3D_sp
     MODULE PROCEDURE HCO_ScaleArr2D_dp
     MODULE PROCEDURE HCO_ScaleArr2D_sp
     MODULE PROCEDURE HCO_ScaleArr1D_dp
     MODULE PROCEDURE HCO_ScaleArr1D_sp
  END INTERFACE
!
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleInit
!
! !DESCRIPTION: Function HCO\_ScaleInit initialized the uniform scale
!  factors for every HEMCO species.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleInit( HcoState, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    LOGICAL            :: FOUND
    CHARACTER(LEN=255) :: iName, MSG
    REAL(hp)           :: ScalFactor

    !--------------------------
    ! HCO_ScaleInit begins here
    !--------------------------

    ! Allocate scale factors and initialize all scale factors to 1.0
    IF ( ALLOCATED(SpcScal) ) DEALLOCATE(SpcScal)
    ALLOCATE( SpcScal( HcoState%nSpc ) )
    SpcScal(:) = 1.0_hp

    ! Search for scale factors in HEMCO configuration file
    DO N = 1, HcoState%nSpc
       iName = 'EmisScale_'//TRIM(HcoState%Spc(N)%SpcName)
       CALL GetExtOpt( HcoState%Config, -999, iName, OptValHp=ScalFactor, &
                       FOUND=FOUND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! If there is a scale factor, update scale vector accordingly
       IF ( FOUND ) THEN
          SpcScal(N) = ScalFactor

          ! Verbose mode
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE (MSG,*) 'Will use universal emission scale factor for ', &
                TRIM(HcoState%Spc(N)%SpcName),': ',SpcScal(N)
             CALL HCO_MSG ( HcoState%Config%Err, MSG )
          ENDIF
       ENDIF
    ENDDO

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleInit
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleGet
!
! !DESCRIPTION: Function HCO\_ScaleGet returns the scale factor for the given
!  species ID.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_ScaleGet ( HcoID ) RESULT ( ScalFact )
!
! !USES
!
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )  :: HcoID          ! HEMCO species ID
!
! !RESULT:
!
    REAL(hp)                        :: ScalFact
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------
    ! HCO_ScaleGet begins here
    !--------------------------
    IF ( HcoID > SIZE(SpcScal,1) .OR. HcoID <= 0 ) THEN
       ScalFact = 1.0_hp
    ELSE
       ScalFact = SpcScal(HcoID)
    ENDIF

  END FUNCTION HCO_ScaleGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr3D_sp
!
! !DESCRIPTION: Function HCO\_ScaleArr3D scales the 3D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr3D_sp( HcoState, HcoID, Arr3D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),        INTENT(INOUT)  :: Arr3D( HcoState%NX, &
                                              HcoState%NY, &
                                              HcoState%NZ )
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr3D_sp begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr3D = Arr3D * ScalFact
       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) '3D field scaled by factor of ',ScalFact
          CALL HCO_MSG ( HcoState%Config%Err, MSG )
       ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr3D_sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr3D_dp
!
! !DESCRIPTION: Function HCO\_ScaleArr3D scales the 3D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr3D_dp( HcoState, HcoID, Arr3D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp),        INTENT(INOUT)  :: Arr3D( HcoState%NX, &
                                              HcoState%NY, &
                                              HcoState%NZ )
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr3D_dp begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr3D = Arr3D * ScalFact
       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) '3D field scaled by factor of ',ScalFact
          CALL HCO_MSG ( HcoState%Config%Err, MSG )
       ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr3D_dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr2D_sp
!
! !DESCRIPTION: Function HCO\_ScaleArr2D scales the 2D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr2D_sp( HcoState, HcoID, Arr2D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),        INTENT(INOUT)  :: Arr2D( HcoState%NX, &
                                              HcoState%NY )
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr2D begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr2D = Arr2D * ScalFact
       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) '2D field scaled by factor of ',ScalFact
          CALL HCO_MSG ( HcoState%Config%Err, MSG )
       ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr2D_sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr2D_dp
!
! !DESCRIPTION: Function HCO\_ScaleArr2D scales the 2D array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr2D_dp( HcoState, HcoID, Arr2D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp),        INTENT(INOUT)  :: Arr2D( HcoState%NX, &
                                              HcoState%NY )
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr2D begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr2D = Arr2D * ScalFact
       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) '2D field scaled by factor of ',ScalFact
          CALL HCO_MSG ( HcoState%Config%Err, MSG )
       ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr2D_dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr1D_sp
!
! !DESCRIPTION: Function HCO\_ScaleArr1D scales a single value.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr1D_sp( HcoState, HcoID, Arr1D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp),        INTENT(INOUT)  :: Arr1D
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr2D begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr1D = Arr1D * ScalFact
       ! Verbose mode
       !IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       !   WRITE(MSG,*) '1D field scaled by factor of ',ScalFact
       !   CALL HCO_MSG ( HcoState%Config%Err, MSG )
       !ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr1D_sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleArr1D_dp
!
! !DESCRIPTION: Function HCO\_ScaleArr1D scales a single value.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleArr1D_dp( HcoState, HcoID, Arr1D, RC )
!
! !USES
!
    USE HCO_STATE_MOD,    ONLY : HCO_STATE
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState       ! HEMCO state object
    INTEGER,         INTENT(IN   )  :: HcoID          ! Species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(dp),        INTENT(INOUT)  :: Arr1D
    INTEGER ,        INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)           :: ScalFact
    CHARACTER(LEN=255) :: MSG

    !--------------------------
    ! HCO_ScaleArr2D begins here
    !--------------------------
    ScalFact = HCO_ScaleGet( HcoID )
    IF ( .NOT. HcoState%Options%ScaleEmis ) ScalFact = 1.0_hp
    IF ( ScalFact /= 1.0_hp ) THEN
       Arr1D = Arr1D * ScalFact
       ! Verbose mode
       !IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       !   WRITE(MSG,*) '1D field scaled by factor of ',ScalFact
       !   CALL HCO_MSG ( HcoState%Config%Err, MSG )
       !ENDIF
    ENDIF

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_ScaleArr1D_dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ScaleFinal
!
! !DESCRIPTION: Function HCO\_ScaleFinal finalizes the module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ScaleFinal()
!
! !USES
!
!
! !INPUT PARAMETERS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 May 2017 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------
    ! HCO_ScaleFinal begins here
    !--------------------------

    ! Allocate scale factors and initialize all scale factors to 1.0
    IF ( ALLOCATED(SpcScal) ) DEALLOCATE(SpcScal)

  END SUBROUTINE HCO_ScaleFinal
!EOC
END MODULE HCO_Scale_Mod
!EOM
