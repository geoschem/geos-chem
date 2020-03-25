!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_interface_mod.F90
!
! !DESCRIPTION: Module hco\_interface\_mod.F90 contains routines and
! variables to interface GEOS-Chem and HEMCO. It contains the HEMCO
! state object (HcoState) as well as some wrapper routines to exchange
! values and fields between HEMCO and GEOS-Chem.
!\\
!\\
! The HEMCO driver routines are located in hcoi\_gc\_main\_mod.F90.
! This module just contains some high level variables and routines
! that can be accessed from everywhere within GeosCore.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_INTERFACE_MOD
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCOX_State_Mod, ONLY : Ext_State
  USE HCO_State_Mod,  ONLY : HCO_State
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: GetHcoVal
  PUBLIC  :: GetHcoID
  PUBLIC  :: GetHcoDiagn
  PUBLIC  :: SetHcoTime
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Feb 2016 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MODULE VARIABLES:
!
  !--------------------------
  ! %%% Pointers %%%
  !--------------------------

  ! HEMCO state
  TYPE(HCO_State), POINTER, PUBLIC :: HcoState => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State), POINTER, PUBLIC :: ExtState => NULL()
!
! !DEFINED PARAMETERS:
!

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetHcoTime
!
! !DESCRIPTION: SUBROUTINE SetHcoTime sets the current simulation
! datetime in HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetHcoTime( TimeForEmis, RC )
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_Set
    USE TIME_MOD,      ONLY : GET_YEAR, GET_MONTH,  GET_DAY
    USE TIME_MOD,      ONLY : GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: TimeForEmis
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: cYr, cMt, cDy, cHr, cMin, cSec, cDOY

    !=================================================================
    ! SetHcoTime begins here
    !=================================================================

    cYr      = GET_YEAR()
    cMt      = GET_MONTH()
    cDy      = GET_DAY()
    cHr      = GET_HOUR()
    cMin     = GET_MINUTE()
    cSec     = GET_SECOND()
    cDOY     = GET_DAY_OF_YEAR()

    CALL HcoClock_Set ( HcoState, cYr, cMt, cDy, cHr, &
                        cMin, cSec, cDoy, IsEmisTime=TimeForEmis, RC=RC )

  END SUBROUTINE SetHcoTime
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoVal
!
! !DESCRIPTION: Subroutine GetHcoVal is a wrapper routine to return an
! emission (kg/m2/s) or deposition (1/s) value from the HEMCO state object
! for a given GEOS-Chem tracer at position I, J, L.
! A value of zero is returned if no HEMCO species is defined for the given
! tracer, and the output parameter Found is set to false.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoVal ( TrcID, I, J, L, Found, Emis, Dep )
!
! !USES
!
!
! !INPUT ARGUMENTS:
!
    INTEGER,            INTENT(IN   )  :: TrcID   ! GEOS-Chem tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L ! Position
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: FOUND   ! Was this tracer ID found?
    REAL(hp), OPTIONAL, INTENT(  OUT)  :: Emis    ! Emissions  [kg/m2/s]
    REAL(hp), OPTIONAL, INTENT(  OUT)  :: Dep     ! Deposition [1/s]
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER   :: HcoID, tID

    !=================================================================
    ! GetHcoVal begins here
    !=================================================================

    ! Init
    FOUND = .FALSE.
    IF ( PRESENT(Emis) ) Emis = 0.0_hp
    IF ( PRESENT(Dep ) ) Dep  = 0.0_hp

    ! Define tracer ID to be used.
    HcoID = TrcID

!    ! HEMCO species ID corresponding to this GEOS-Chem tracer
!    IF ( tID > 0 ) HcoID = M2HID(tID)%ID

    ! If HEMCO species exists, get value from HEMCO state
    IF ( HcoID > 0 ) THEN
       IF ( PRESENT(Emis) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(HcoID)%Emis%Val) ) THEN
             Emis  = HcoState%Spc(HcoID)%Emis%Val(I,J,L)
             FOUND = .TRUE.
          ENDIF
       ENDIF
       IF ( PRESENT(Dep) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(HcoID)%Depv%Val) ) THEN
             Dep   = HcoState%Spc(HcoID)%Depv%Val(I,J)
             FOUND = .TRUE.
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE GetHcoVal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoID
!
! !DESCRIPTION: Function GetHcoID is a convenience wrapper function to
! return the HEMCO ID by name or by GC tracer ID.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetHcoID( name, SpcID ) RESULT ( HcoID )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: Name  ! Tracer name
    INTEGER,          INTENT(IN   ), OPTIONAL :: SpcID ! Tracer ID
!
! !OUTPUT PARAMETERS:
!
    INTEGER                                   :: HcoID
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Init
    HcoID = -1

    ! To get HEMCO ID by tracer ID
    IF ( PRESENT(SpcID) ) THEN
!       IF ( TrcID > 0 ) HcoID = M2HID(TrcID)%ID
       IF ( SpcID > 0 ) HcoID = SpcID
    ENDIF
    IF ( PRESENT(name) ) THEN
       HcoID = HCO_GetHcoID( name, HcoState )
    ENDIF

  END FUNCTION GetHcoID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoDiagn
!
! !DESCRIPTION: Subroutine GetHcoDiagn is a convenience wrapper routine to
! get a HEMCO diagnostics from somewhere within GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoDiagn ( DiagnName, StopIfNotFound, RC, &
                           Ptr2D,     Ptr3D,     COL, AutoFill        )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_TYPES_MOD,      ONLY : DiagnCont
    USE HCO_DIAGN_MOD
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)           :: DiagnName      ! Name of diagnostics
    LOGICAL,          INTENT(IN)           :: StopIfNotFound ! Stop if diagnostics
                                                             ! does not exist?
    INTEGER,          INTENT(IN), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN), OPTIONAL :: AutoFill       ! AutoFill diagnostics only?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)        :: RC             ! Error return code
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         POINTER, OPTIONAL    :: Ptr2D(:,:)      ! Pointer to 2D data
    REAL(sp),         POINTER, OPTIONAL    :: Ptr3D(:,:,:)    ! Pointer to 3D data
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Sep 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: FLAG, LevIDx, PS, AF
    TYPE(DiagnCont), POINTER  :: DgnCont  => NULL()

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=======================================================================
    ! GetHcoDiagn begins here
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at GetHcoDiagn (in module GeosCore/hcoi_gc_diagn_mod.F90)'

    ! Set collection number
    PS = HcoState%Diagn%HcoDiagnIDManual
    IF ( PRESENT(COL) ) PS = COL

    ! Set AutoFill flag
    AF = -1
    IF ( PRESENT(AutoFill) ) AF = AutoFill

    ! Get diagnostics by name. Search all diagnostics, i.e. both AutoFill
    ! and manually filled diagnostics. Also include those with a manual
    ! output interval.
    CALL Diagn_Get( HcoState,    .FALSE., DgnCont,                 &
                    FLAG,        RC,      cName=TRIM(DiagnName),   &
                    AutoFill=AF, COL=PS                            )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error in getting diagnostics: ' // TRIM(DiagnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( (FLAG /= HCO_SUCCESS) .AND. StopIfNotFound ) THEN
       ErrMsg = 'Cannot get diagnostics for this time stamp: ' //    &
                 TRIM(DiagnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Pass data to output pointer (only if diagnostics defined):
    IF ( FLAG == HCO_SUCCESS ) THEN

       ! 2D pointer
       IF ( PRESENT(Ptr2D) ) THEN

          ! Pass 2D data
          IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
             Ptr2D => DgnCont%Arr2D%Val

          ! Pass 3D data. Get level index from diagnostics (if set)
          ELSEIF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             LevIDx = DgnCont%LevIdx
             IF ( LevIdx < 1 ) LevIdx = 1
             Ptr2D => DgnCont%Arr3D%Val(:,:,LevIDx)

          ! Error if no 2D or 3D data available
          ELSE
             ErrMsg = 'no data defined: '// TRIM(DiagnName)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ! 3D pointer: must point to 3D data
       ELSEIF ( PRESENT(Ptr3D) ) THEN
          IF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             Ptr3D => DgnCont%Arr3D%Val
          ELSE
             ErrMsg = 'no 3D data defined: '// TRIM(DiagnName)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ! Error otherwise
       ELSE
          ErrMsg = 'Please define output data pointer: ' // TRIM(DiagnName)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Free pointer
    DgnCont  => NULL()

    ! Leave with success
    RC = HCO_SUCCESS

  END SUBROUTINE GetHcoDiagn
!EOC
END MODULE HCO_INTERFACE_MOD
