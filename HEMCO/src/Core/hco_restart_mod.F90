!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_restart_mod.F90
!
! !DESCRIPTION: Module HCO\_RESTART\_MOD contains wrapper routines to define,
! get and write restart fields.
!\\
!\\
! Restart variables are required by some of the HEMCO extensions. The
! HEMCO restart variables can be organized through the HEMCO restart
! diagnostics collection. At the end of a simulation, all diagnostic fields
! ('containers') of the restart collection are written to the HEMCO restart
! file.
!\\
!\\
! All fields from the HEMCO restart file can be easily read back into HEMCO
! via the HEMCO I/O infrastructure, i.e. by listing them in the HEMCO
! configuration file.
!\\
!\\
! In an ESMF/MAPL environment, restart variables should be organized through
! the ESMF internal state object. This is particularly important for
! simulation that rely on checkpoint files (e.g. replay simulations).
! In this cases, the restart variables are obtained from / written to the
! ESMF internal state object, rather than the HEMCO restart file.
!\\
!\\
! This module contains wrapper routines to define, obtain and write restart
! fields for the two aforementioned restart types. The routines work both for
! 'traditional' and ESMF restart variables. In an ESMF application, the first
! check is always performed within the internal state, e.g. it is first
! checked if the given field variable exists in the internal state of the
! gridded component that HEMCO sits in. If so, the field is obtained from /
! written to the internal state. If no internal state object exist, an
! attempt is made to obtain the field through the HEMCO data list, i.e. it
! is checked if the restart variable is specified in the HEMCO configuration
! file.
! A HEMCO diagnostics container is created in the diagnostics restart
! collection in both the traditional and the ESMF environment.
!\\
!\\
! Routine HCO\_RestartDefine should be called during the initialization stage.
! Routine HCO\_RestartGet should be called on the first run call and after
! each rewinding of the clock. HEMCO routines HcoClock\_First and
! HcoClock\_Rewind can be used to determine if it's time to read/update the
! restart variable.
! HCO\_RestartWrite should be called *on every time step*. This is important
! in ESMF applications that rely on checkpoint files (e.g. replay simulations)
! that are written out by ESMF/MAPL throughout the simulation. In a non-ESMF
! environment, the HCO\_RestartWrite call is basically void but it should be
! called nevertheless.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_RESTART_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  ! defined in all environment
  PUBLIC :: HCO_RestartDefine
  PUBLIC :: HCO_RestartGet
  PUBLIC :: HCO_RestartWrite
!
! !PRIVATE MEMBER FUNCTIONS:
!
#if defined(ESMF_)
  PRIVATE :: HCO_CopyFromIntnal_ESMF
#endif

  INTERFACE HCO_RestartDefine
     MODULE PROCEDURE HCO_RestartDefine_3D
     MODULE PROCEDURE HCO_RestartDefine_2D
  END INTERFACE HCO_RestartDefine

  INTERFACE HCO_RestartGet
     MODULE PROCEDURE HCO_RestartGet_3D
     MODULE PROCEDURE HCO_RestartGet_2D
  END INTERFACE HCO_RestartGet

  INTERFACE HCO_RestartWrite
     MODULE PROCEDURE HCO_RestartWrite_3D
     MODULE PROCEDURE HCO_RestartWrite_2D
  END INTERFACE HCO_RestartWrite
!
! !REVISION HISTORY:
!  10 Mar 2015 - C. Keller   - Initial version
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
! !ROUTINE: HCO_RestartDefine_3D
!
! !DESCRIPTION: Subroutine HCO\_RestartDefine\_3D defines a restart diagnostics.
! This adds a diagnostics with output frequency 'End' to the HEMCO diagnostics
! list. Arr3D is the 3D field of interest. The diagnostics will not copy the
! current content of Arr3D but establish a 'link' (e.g. pointer) to it. This
! way, any updates to Arr3D will automatically be seen by the diagnostics
! and there is no need to explicitly update the content of the diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartDefine_3D( HcoState, Name, Arr3D, Unit, RC )
!
! !USES:
!
    USE HCO_DIAGN_MOD,    ONLY : Diagn_Create
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER               :: HcoState  ! HEMCO state obj.
    CHARACTER(LEN=*),    INTENT(IN   )         :: Name      ! Name of restart variable
    ! Array with data of interest
    REAL(sp),            INTENT(IN   ), TARGET :: Arr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
    CHARACTER(LEN=*),    INTENT(IN   )         :: Unit      ! Units of Arr3D
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,             INTENT(INOUT)         :: RC        ! Return code
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! ================================================================
    ! HCO_RestartDefine_3D begins here
    ! ================================================================

    ! Define diagnostics array
    CALL Diagn_Create ( HcoState,                               &
                        cName      = TRIM(Name),                &
                        ExtNr      = -1,                        &
                        Cat        = -1,                        &
                        Hier       = -1,                        &
                        HcoID      = -1,                        &
                        SpaceDim   =  3,                        &
                        OutUnit    = TRIM(Unit),                &
                        COL = HcoState%Diagn%HcoDiagnIDRestart, &
                        AutoFill   = 0,                         &
                        Trgt3D     = Arr3D,                     &
                        RC         = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartDefine_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartDefine_2D
!
! !DESCRIPTION: Subroutine HCO\_RestartDefine\_2D defines a restart diagnostics.
! This adds a diagnostics with output frequency 'End' to the HEMCO diagnostics
! list. Arr2D is the 2D field of interest. The diagnostics will not copy the
! current content of Arr2D but establish a 'link' (e.g. pointer) to it. This
! way, any updates to Arr2D will automatically be seen by the diagnostics
! and there is no need to explicitly update the content of the diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartDefine_2D( HcoState, Name, Arr2D, Unit, RC )
!
! !USES:
!
    USE HCO_DIAGN_MOD,    ONLY : Diagn_Create
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER               :: HcoState  ! HEMCO state obj.
    CHARACTER(LEN=*),    INTENT(IN   )         :: Name      ! Name of restart variable
    ! Array with data of interest
    REAL(sp),            INTENT(IN   ), TARGET :: Arr2D(HcoState%NX,HcoState%NY)
    CHARACTER(LEN=*),    INTENT(IN   )         :: Unit      ! Units of Arr2D
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,             INTENT(INOUT)         :: RC        ! Return code
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! ================================================================
    ! HCO_RestartDefine_2D begins here
    ! ================================================================

    ! Define diagnostics array
    CALL Diagn_Create ( HcoState,                               &
                        cName      = TRIM(Name),                &
                        ExtNr      = -1,                        &
                        Cat        = -1,                        &
                        Hier       = -1,                        &
                        HcoID      = -1,                        &
                        SpaceDim   =  2,                        &
                        OutUnit    = TRIM(Unit),                &
                        COL = HcoState%Diagn%HcoDiagnIDRestart, &
                        AutoFill   = 0,                         &
                        Trgt2D     = Arr2D,                     &
                        RC         = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartDefine_2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartGet_3D
!
! !DESCRIPTION: Subroutine HCO\_RestartGet\_3D attempts to read a restart field.
! In an ESMF environment, it first checks if the given field (name) is included
! in the internal state object, in which case the data object is filled with
! these values. If not found or if not in an ESMF environment, the HEMCO data
! list (specified in the HEMCO configuration file) is searched. A default value
! can be specified in case that no field could be imported via ESMF and/or the
! HEMCO interface.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartGet_3D( HcoState, Name,   Arr3D, &
                                RC,       FILLED, Def3D, DefVal  )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER                 :: HcoState  ! HEMCO state object
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name      ! Name of restart variable
    ! Default value to be used if restart variable could not be found
    REAL(sp),            INTENT(IN   ), OPTIONAL :: Def3D(HcoState%NX,HcoState%NY,HcoState%NZ)
    ! Default uniform value to be used if restart variable could not be found and
    ! Def2D is not defined.
    REAL(sp),            INTENT(IN   ), OPTIONAL :: DefVal
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FILLED    ! Was the restart variable found?
!
! !INPUT/OUTPUT ARGUMENTS:
!
    ! Data field with restart variable
    REAL(sp),            INTENT(INOUT)           :: Arr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
    INTEGER,             INTENT(INOUT)           :: RC       ! Return code
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp), POINTER    :: Ptr3D(:,:,:)
    LOGICAL              :: FLD
    CHARACTER(LEN=255)   :: MSG

    ! ================================================================
    ! HCO_RestartGet begins here
    ! ================================================================

    ! Init
    Ptr3D => NULL()

    ! Is the output array filled yet?
    FLD = .FALSE.

    ! ------------------------------------------------------------------
    ! Try to get from ESMF internal state
    ! ------------------------------------------------------------------
#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( HcoState, TRIM(Name),   &
                                  1,        FLD,      RC, Arr3D=Arr3D )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If field is all zero assume it to be not filled
    IF ( FLD ) THEN
       IF ( SUM(Arr3D) == 0.0 ) FLD = .FALSE.
    ENDIF

    ! Log output
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       IF ( HcoState%amIRoot .AND. FLD ) THEN
          MSG = 'Obtained restart variable from ESMF internal state: '//TRIM(Name)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF
#endif

    ! ------------------------------------------------------------------
    ! If not yet filled, try to get from HEMCO configuration
    ! ------------------------------------------------------------------
    IF ( .NOT. FLD ) THEN

       ! Try to get pointer from HEMCO configuration
       CALL HCO_GetPtr( HcoState, TRIM(Name), Ptr3D, RC, FILLED=FLD )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually pass data
       IF ( FLD ) THEN
          Arr3D = Ptr3D

          ! Log output
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Obtained restart variable from HEMCO config: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
          ENDIF
       ENDIF

       ! Cleanup
       Ptr3D => NULL()
    ENDIF

    ! ------------------------------------------------------------------
    ! If still not filled, assign default values
    ! ------------------------------------------------------------------
    IF ( .NOT. FLD ) THEN
       IF ( PRESENT(Def3D) ) THEN
          Arr3D = Def3D
          FLD   = .TRUE.
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Filled restart variable with default 3D field: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
          ENDIF
       ELSEIF( PRESENT(DefVal) ) THEN
          Arr3D = DefVal
          FLD   = .TRUE.
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Filled restart variable with default scalar: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! Leave
    ! ------------------------------------------------------------------
    IF ( PRESENT(FILLED) ) THEN
       FILLED = FLD
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartGet_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartGet_2D
!
! !DESCRIPTION: Subroutine HCO\_RestartGet\_2D attempts to read a restart field.
! In an ESMF environment, it first checks if the given field (name) is included
! in the internal state object, in which case the data object is filled with
! these values. If not found or if not in an ESMF environment, the HEMCO data
! list (specified in the HEMCO configuration file) is searched. A default value
! can be specified in case that no field could be imported via ESMF and/or the
! HEMCO interface.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartGet_2D( HcoState, Name,   Arr2D, &
                                RC,       FILLED, Def2D, DefVal  )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER                 :: HcoState  ! HEMCO state object
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name      ! Name of restart variable
    ! Default value to be used if restart variable could not be found
    REAL(sp),            INTENT(IN   ), OPTIONAL :: Def2D(HcoState%NX,HcoState%NY)
    ! Default uniform value to be used if restart variable could not be found and
    ! Def2D is not defined.
    REAL(sp),            INTENT(IN   ), OPTIONAL :: DefVal
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FILLED    ! Was the restart variable found?
!
! !INPUT/OUTPUT ARGUMENTS:
!
    ! Data field with restart variable
    REAL(sp),            INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,             INTENT(INOUT)           :: RC       ! Return code
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp), POINTER    :: Ptr2D(:,:)
    LOGICAL              :: FLD
    CHARACTER(LEN=255)   :: MSG

    ! ================================================================
    ! HCO_RestartGet_2D begins here
    ! ================================================================

    ! Init
    Ptr2D => NULL()

    ! Is the output array filled yet?
    FLD = .FALSE.

    ! ------------------------------------------------------------------
    ! Try to get from ESMF internal state
    ! ------------------------------------------------------------------
#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( HcoState, TRIM(Name),   &
                                  1,        FLD,      RC, Arr2D=Arr2D )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If field is all zero assume it to be not filled
    IF ( FLD ) THEN
       IF ( SUM(Arr2D) == 0.0 ) FLD = .FALSE.
    ENDIF

    ! Log output
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       IF ( HcoState%amIRoot .AND. FLD ) THEN
          MSG = 'Obtained restart variable from ESMF internal state: '//TRIM(Name)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(*,*) TRIM(MSG)
       ENDIF
    ENDIF
#endif

    ! ------------------------------------------------------------------
    ! If not yet filled, try to get from HEMCO configuration
    ! ------------------------------------------------------------------
    IF ( .NOT. FLD ) THEN

       ! Try to get pointer from HEMCO configuration
       CALL HCO_GetPtr( HcoState, TRIM(Name), Ptr2D, RC, FOUND=FLD )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually pass data
       IF ( FLD ) THEN
          Arr2D = Ptr2D

          ! Log output
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Obtained restart variable from HEMCO config: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
                WRITE(*,*) TRIM(MSG)
             ENDIF
          ENDIF
       ENDIF

       ! Cleanup
       Ptr2D => NULL()
    ENDIF

    ! ------------------------------------------------------------------
    ! If still not filled, assign default values
    ! ------------------------------------------------------------------
    IF ( .NOT. FLD ) THEN
       IF ( PRESENT(Def2D) ) THEN
          Arr2D = Def2D
          FLD   = .TRUE.
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Filled restart variable with default 2D field: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
                WRITE(*,*) TRIM(MSG)
             ENDIF
          ENDIF
       ELSEIF( PRESENT(DefVal) ) THEN
          Arr2D = DefVal
          FLD   = .TRUE.
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             IF ( HcoState%amIRoot ) THEN
                MSG = 'Filled restart variable with default scalar: '//TRIM(Name)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
                WRITE(*,*) TRIM(MSG)
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! Leave
    ! ------------------------------------------------------------------
    IF ( PRESENT(FILLED) ) FILLED = FLD

    ! Verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       IF ( HcoState%amIRoot .AND. .NOT. FLD ) THEN
          MSG = 'No restart field found: '//TRIM(Name)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(*,*) TRIM(MSG)
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartGet_2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartWrite_3D
!
! !DESCRIPTION: Subroutine HCO\_RestartWrite\_3D writes a restart variable to
! the ESMF internal state. This is only of relevance in an ESMF environment.
! The 'regular' HEMCO diagnostics created in HCO\_RestartDefine becomes
! automatically written to disk.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartWrite_3D( HcoState, Name, Arr3D, RC, FOUND )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER                 :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FOUND
!
! !INPUT/OUTPUT ARGUMENTS:
!
    REAL(sp),            INTENT(INOUT)           :: Arr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
    INTEGER,             INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL   :: WRITTEN

    ! ================================================================
    ! HCO_RestartWrite_3D begins here
    ! ================================================================

    ! Data written to internal state?
    WRITTEN = .FALSE.

#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( HcoState, TRIM(Name), &
                                  -1,       WRITTEN,    RC, Arr3D=Arr3D )
#endif

    ! Pass to output
    IF ( PRESENT(FOUND) ) THEN
       FOUND = WRITTEN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartWrite_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartWrite_2D
!
! !DESCRIPTION: Subroutine HCO\_RestartWrite\_2D writes a restart variable to
! the ESMF internal state. This is only of relevance in an ESMF environment.
! The 'regular' HEMCO diagnostics created in HCO\_RestartDefine becomes
! automatically written to disk.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartWrite_2D( HcoState, Name, Arr2D, RC, FOUND )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_State),     POINTER                 :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FOUND
!
! !INPUT/OUTPUT ARGUMENTS:
!
    REAL(sp),            INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,             INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  11 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL   :: WRITTEN

    ! ================================================================
    ! HCO_RestartWrite_2D begins here
    ! ================================================================

    ! Data written to internal state?
    WRITTEN = .FALSE.

#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( HcoState, TRIM(Name), &
                                  -1,       WRITTEN,    RC, Arr2D=Arr2D )
#endif

    ! Pass to output
    IF ( PRESENT(FOUND) ) THEN
       FOUND = WRITTEN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartWrite_2D
!EOC
#if defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_CopyFromIntnal_ESMF
!
! !DESCRIPTION: Subroutine HCO\_CopyFromIntnal\_ESMF attempts to transfer
! data to and from the ESMF/MAPL internal state.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CopyFromIntnal_ESMF ( HcoState,  Name,    &
                                       Direction, Found, RC, Arr2D, Arr3D )
!
! !USES:
!
#include "MAPL_Generic.h"
    USE ESMF
    USE MAPL_Mod
    USE HCO_STATE_MOD,   ONLY : Hco_State
!
! !ARGUMENTS:
!
    TYPE(HCO_State),     POINTER                 :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name
    INTEGER,             INTENT(IN   )           :: Direction    ! 1: internal to Arr2D; -1: Arr2D to internal
    LOGICAL,             INTENT(  OUT)           :: Found
    INTEGER,             INTENT(INOUT)           :: RC
    REAL(sp),            INTENT(INOUT), OPTIONAL :: Arr2D(HcoState%NX,HcoState%NY)
    REAL(sp),            INTENT(INOUT), OPTIONAL :: Arr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
!
! !REVISION HISTORY:
!  10 Mar 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                      :: STAT
    TYPE(MAPL_MetaComp), POINTER :: STATE
    TYPE(ESMF_STATE)             :: INTERNAL
    REAL,                POINTER :: Ptr2D(:,:)
    REAL,                POINTER :: Ptr3D(:,:,:)

    ! ================================================================
    ! HCO_CopyFromIntnal_ESMF begins here
    ! ================================================================

    ! For MAPL/ESMF error handling (defines Iam and STATUS)
    __Iam__('HCO_CopyFromIntnal_ESMF (HCOI_ESMF_MOD.F90)')

    ! Init
    Ptr2D => NULL()
    Ptr3D => NULL()

    ! Get internal state
    CALL MAPL_GetObjectFromGC( HcoState%GridComp, STATE, __RC__ )
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

    ! Try to import field
    IF ( PRESENT(Arr2D) ) THEN
       CALL MAPL_GetPointer( INTERNAL, Ptr2D, TRIM(Name), &
                             NotFoundOk=.TRUE., __RC__ )

       ! Eventually copy data to or from output array
       IF ( ASSOCIATED(Ptr2D) ) THEN

          ! Make sure we can copy the data
          ASSERT_(SIZE(Arr2D,1)==SIZE(Ptr2D,1))
          ASSERT_(SIZE(Arr2D,2)==SIZE(Ptr2D,2))

          ! transfer direction must be 1 or -1
          ASSERT_(Direction==1 .OR. Direction==-1)

          ! Transfer data
          IF ( Direction == 1 ) THEN
             Arr2D = Ptr2D
          ELSEIF ( Direction == -1 ) THEN
             Ptr2D = Arr2D
          ENDIF
          Found = .TRUE.
       ELSE
          Found = .FALSE.
       ENDIF

       ! Cleanup
       Ptr2D => NULL()
    ENDIF

    ! Try to import field
    IF ( PRESENT(Arr3D) ) THEN
       CALL MAPL_GetPointer( INTERNAL, Ptr3D, TRIM(Name), &
                             NotFoundOk=.TRUE., __RC__ )

       ! Eventually copy data to or from output array
       IF ( ASSOCIATED(Ptr3D) ) THEN

          ! Make sure we can copy the data
          ASSERT_(SIZE(Arr3D,1)==SIZE(Ptr3D,1))
          ASSERT_(SIZE(Arr3D,2)==SIZE(Ptr3D,2))
          ASSERT_(SIZE(Arr3D,3)==SIZE(Ptr3D,3))

          ! transfer direction must be 1 or -1
          ASSERT_(Direction==1 .OR. Direction==-1)

          ! Transfer data
          IF ( Direction == 1 ) THEN
             Arr3D = Ptr3D
          ELSEIF ( Direction == -1 ) THEN
             Ptr3D = Arr3D
          ENDIF
          Found = .TRUE.
       ELSE
          Found = .FALSE.
       ENDIF

       ! Cleanup
       Ptr3D => NULL()
    ENDIF

    ! Return success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_CopyFromIntnal_ESMF
!EOC
#endif
END MODULE HCO_RESTART_MOD
