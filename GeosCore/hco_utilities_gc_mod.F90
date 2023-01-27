!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_utilities_gc_mod.F90
!
! !DESCRIPTION: Module hco\_utilities\_gc\_mod.F90 is the high-level component
!  used by GEOS-Chem routines to receive emissions and generally data from
!  HEMCO. It interacts with HEMCO routines directly but also GEOS-Chem data
!  structures.
!\\
!\\
!  Please note the following guidelines for inclusion of routines in this module
!  to prevent feature creep like the former HCOI\_GC\_Main\_Mod:
!  - The routines must use both HEMCO and GEOS-Chem data structures. If not,
!    they may be better suited for EMISSIONS\_MOD or FLEXGRID\_READ\_MOD,
!    or implemented as common tools in HEMCO/Interfaces/hco\_interface\_common.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Utilities_GC_Mod
!
! !USES:
!
  USE ErrCode_Mod
  USE Precision_Mod
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: InquireHco         ! Inquire availability of emis/drydep field
  PUBLIC   :: LoadHcoValEmis
  PUBLIC   :: LoadHcoValDep
  PUBLIC   :: GetHcoValEmis
  PUBLIC   :: GetHcoValDep
  PUBLIC   :: HCO_GC_EvalFld     ! Shim interface for HCO_EvalFld
  PUBLIC   :: HCO_GC_GetPtr      ! Shim interface for HCO_GetPtr
  PUBLIC   :: HCO_GC_GetDiagn    ! Shim interface for GetHcoDiagn
  PUBLIC   :: HCO_GC_GetOption   ! Shim interface for HCO_GetOpt
  PUBLIC   :: HCO_GC_HcoStateOK  ! Test if HCO_State is allocated

#if defined( MODEL_CLASSIC )
  !=========================================================================
  ! These are only needed for GEOS-Chem "Classic"
  ! Intermediate grid (IMGrid) functionality
  !=========================================================================
  PUBLIC   :: Regrid_MDL2HCO
  PUBLIC   :: Regrid_HCO2MDL
  PUBLIC   :: Init_IMGrid

  !=========================================================================
  ! These are only needed for GEOS-Chem "Classic"
  !=========================================================================
  PUBLIC   :: Get_GC_Restart
  PUBLIC   :: Get_Boundary_Conditions
#endif

!
! !REMARKS:
!  Mostly wrapper functions migrated from the former HCO_Interface_Mod
!
! !REVISION HISTORY:
!  12 Mar 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
  INTERFACE HCO_GC_EvalFld
    MODULE PROCEDURE HCO_GC_EvalFld_2D
    MODULE PROCEDURE HCO_GC_EvalFld_3D
  END INTERFACE HCO_GC_EvalFld

  INTERFACE HCO_GC_GetPtr
    MODULE PROCEDURE HCO_GC_GetPtr_2D
    MODULE PROCEDURE HCO_GC_GetPtr_3D
  END INTERFACE HCO_GC_GetPtr

  INTERFACE HCO_GC_GetDiagn
    MODULE PROCEDURE HCO_GC_GetDiagn_2D
    MODULE PROCEDURE HCO_GC_GetDiagn_3D
  END INTERFACE HCO_GC_GetDiagn
!
! !PRIVATE TYPES:
!
#if defined( MODEL_CLASSIC )
  !------------------------------------------------------
  ! HEMCO Intermediate Grid functionality.
  ! Array buffers for storing regridded data. These are essentially
  ! temporary pointer targets.
  !
  ! They are refreshed every Regrid_x2y so do not point data to here
  ! beyond one subroutine call.
  !
  ! For performance optimization, there are TWO model buffers (H2M, H2Mb)
  ! which allow for some edge cases where GetHcoValDep is interspersed
  ! with GetHcoValEmis, or there may be two species. By default,
  ! GetHcoValDep uses H2Mb and GetHcoValEmis can choose the default or b,
  ! using an optional argument. This allows shimming as much of the
  ! operating specifics from GEOS-Chem core code as possible.
  !------------------------------------------------------
  REAL(hp), POINTER, PUBLIC            :: TMP_MDL (:,:,:)
  REAL(hp), POINTER                    :: TMP_MDLb(:,:,:)
  REAL(hp), POINTER                    :: TMP_HCO (:,:,:)

  ! f4 variant temporaries.
  ! not directly used for regrid, used for pointing and downgrading data
  REAL(f4), POINTER                    :: TMP_MDL_r4 (:,:,:)
  REAL(f4), POINTER                    :: TMP_MDL_r4b(:,:,:)
  REAL(f4), POINTER                    :: TMP_HCO_r4 (:,:,:)

  CHARACTER(LEN=90)                    :: LAST_TMP_REGRID_M2H       ! Last regridded container name
  CHARACTER(LEN=90)                    :: LAST_TMP_REGRID_H2M       ! ... HEMCO to Model
  CHARACTER(LEN=90)                    :: LAST_TMP_REGRID_H2Mb      ! ... HEMCO to Model (alt bfr)
  INTEGER                              :: LAST_TMP_MDL_ZBND         ! Last z-boundary for TMP_MDL_r4 ptr

  ! Temporaries for Map_A2A shadow input variables.
  ! Only need to be allocated once.
  REAL(hp), POINTER                    :: LonEdgeH(:)               ! HEMCO lon, lat edges (NX+1, NY+1)
  REAL(hp), POINTER                    :: LatEdgeH(:)

  REAL(hp), POINTER                    :: LonEdgeM(:)               ! Model lon, lat edges (NX+1, NY+1)
  REAL(hp), POINTER                    :: LatEdgeM(:)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InquireHco
!
! !DESCRIPTION: Subroutine InquireHco INQUIRES to the HEMCO emissions list whether
!  the given tracer ID has emissions or dry deposition spec.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InquireHco ( TrcID, Emis, Dep )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoVal
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT ARGUMENTS:
!
    INTEGER,            INTENT(IN   )  :: TrcID      ! GEOS-Chem tracer ID
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL, OPTIONAL,  INTENT(  OUT)  :: Dep        ! Dep?
    LOGICAL, OPTIONAL,  INTENT(  OUT)  :: Emis       ! Emis?
!
! !REMARKS:
!  Note this assumes TrcID == HcoID.
!
! !REVISION HISTORY:
!  13 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
  IF ( PRESENT( Dep ) ) THEN
    Dep = ASSOCIATED( HcoState%Spc(TrcID)%Depv%Val )
  ENDIF

  IF ( PRESENT(Emis) ) THEN
    Emis = ASSOCIATED( HcoState%Spc(TrcID)%Emis%Val )
  ENDIF

  END SUBROUTINE InquireHco
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LoadHcoValEmis
!
! !DESCRIPTION: For GC-Classic intermediate grid: Load emissions regridded onto
!  model grid into the regridding buffer. Does nothing in other models.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LoadHcoValEmis ( Input_Opt, State_Grid, TrcID, AltBuffer )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    INTEGER,            INTENT(IN   )  :: TrcID      ! GEOS-Chem tracer ID
    LOGICAL, OPTIONAL,  INTENT(IN   )  :: AltBuffer  ! Alternate buffer? (Use B)
!
! !REMARKS:
!  This is achieved through a OMP CRITICAL failsafe and calls GetHcoValEmis/Dep to
!  trigger the regridding.
!
! !REVISION HISTORY:
!  27 Sep 2020 - H.P. Lin  - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
!
    LOGICAL            :: TMP_Found
    REAL(hp)           :: TMP_Value                   ! Dummy values

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN

      ! Check if we have to load the data.
      IF ( TrcID > 0 .and. (.not. ASSOCIATED(HcoState%Spc(TrcID)%Emis%Val)) ) RETURN

      ! The below section must be OMP CRITICAL because it is stateful.
      ! The first call to the critical section will update the container!!
      !$OMP CRITICAL

      ! due to a compiler bug in ifort 19
      ! we have to use PRESENT and not copy the value, as sometimes it becomes
      ! flipped! (hplin, 9/29/20)
      IF ( PRESENT(AltBuffer) ) THEN
        CALL GetHcoValEmis ( Input_Opt, State_Grid, TrcID, 1, 1, 1, Found=TMP_Found, &
                           Emis=TMP_Value, AltBuffer=.true. )
      ELSE
        CALL GetHcoValEmis ( Input_Opt, State_Grid, TrcID, 1, 1, 1, Found=TMP_Found, &
                           Emis=TMP_Value )
      ENDIF

      IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# LoadHcoValEmis/ImGrid: Loading", TrcID, PRESENT( AltBuffer )

      !$OMP END CRITICAL
      ! End of LIMGRID OMP Critical section
    ENDIF
#endif
  END SUBROUTINE LoadHcoValEmis
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LoadHcoValDep
!
! !DESCRIPTION: For GC-Classic intermediate grid: Load deposition value regridded onto
!  model grid into the regridding buffer. Does nothing in other models.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LoadHcoValDep ( Input_Opt, State_Grid, TrcID )
!
! !USES:
!
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    INTEGER,            INTENT(IN   )  :: TrcID      ! GEOS-Chem tracer ID
!
! !REMARKS:
!  This is achieved through a OMP CRITICAL failsafe and calls GetHcoValEmis/Dep to
!  trigger the regridding.
!
! !REVISION HISTORY:
!  27 Sep 2020 - H.P. Lin  - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
!
    LOGICAL            :: TMP_Found
    REAL(hp)           :: TMP_Value                   ! Dummy values

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN

      ! Check if we have to load the data.
      IF ( TrcID > 0 .and. (.not. ASSOCIATED(HcoState%Spc(TrcID)%Emis%Val)) ) RETURN

      ! The below section must be OMP CRITICAL because it is stateful.
      ! The first call to the critical section will update the container!!
      !$OMP CRITICAL

      CALL GetHcoValDep ( Input_Opt, State_Grid, TrcID, 1, 1, 1, Found=TMP_Found, &
                          Dep=TMP_Value )
      !$OMP END CRITICAL
      ! End of LIMGRID OMP Critical section

      IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# LoadHcoValDep/ImGrid: Loading", TrcID
    ENDIF
#endif
  END SUBROUTINE LoadHcoValDep
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoValEmis
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
  SUBROUTINE GetHcoValEmis ( Input_Opt, State_Grid, TrcID, I, J, L, Found, Emis, AltBuffer, SkipCheck )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoVal
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState

    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    INTEGER,            INTENT(IN   )  :: TrcID      ! GEOS-Chem tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L    ! Position
    LOGICAL, OPTIONAL,  INTENT(IN   )  :: AltBuffer  ! Alternate buffer? (Use B)
    LOGICAL, OPTIONAL,  INTENT(IN   )  :: SkipCheck  ! Skip buffer validity check - Dangerous, use in tight loops
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: Found      ! Was this tracer ID found?
    REAL(hp),           INTENT(  OUT)  :: Emis       ! Emissions  [kg/m2/s]
!
! !REMARKS:
!  This subroutine is currently just a stub to call the equivalent in HEMCO
!  utilities. This assumes that HEMCO and GEOS-Chem grids match.
!  When the intermediate grid is implemented, a regridding routine will live
!  here and regrid data on-demand.
!
!  This also assumes that TrcID matches HEMCO tracer ID. If not, a mapping
!  needs to be performed here.
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  12 Mar 2020 - H.P. Lin  - Now wrapper around common utilities
!  05 Jun 2020 - H.P. Lin  - Add GC-Classic on-demand regridding
!EOP
!------------------------------------------------------------------------------
!BOC
#ifdef MODEL_CLASSIC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=32)  :: TMP_TrcIDFldName            ! Temporary tracer field name _HCO_Trc_<id>
    INTEGER            :: ZBND
    REAL(hp), POINTER  :: TMP_MDL_target(:,:,:)       ! Pointer to ease switcheroo of the model target buffer

    TMP_MDL_target => NULL()

    IF ( Input_Opt%LIMGRID ) THEN
      ! Check if we have to load the data.
      IF ( TrcID > 0 .and. (.not. ASSOCIATED(HcoState%Spc(TrcID)%Emis%Val)) ) RETURN

      ! due to a compiler bug in ifort 19
      ! we have to use PRESENT and not copy the value, as sometimes it becomes
      ! flipped! (hplin, 9/29/20)
      IF ( PRESENT( AltBuffer ) ) THEN
        TMP_MDL_target => TMP_MDLb
      ELSE
        TMP_MDL_target => TMP_MDL
      ENDIF

      ! ... on-demand intermediate regridding. Check if we already have this field
      Found = .false.

      WRITE(TMP_TrcIDFldName, '(a,i4)') "_HCO_Trc_", TrcID
      IF( .not. PRESENT( SkipCheck ) .and. &
          .not. ( (.not. PRESENT( AltBuffer ) .and. TMP_TrcIDFldName == LAST_TMP_REGRID_H2M) .or. &
                  (      PRESENT( AltBuffer ) .and. TMP_TrcIDFldName == LAST_TMP_REGRID_H2Mb) ) ) THEN   ! Not already in buffer
        ! Do not use GetHcoVal: load the entire chunk into memory
        ! Note: TrcID matches HcoID here. If not, remap the tracer ID to HEMCO ID below.

        IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValEmis/ImGrid: Attempting to load", TMP_TrcIDFldName, "was", LAST_TMP_REGRID_H2M, LAST_TMP_REGRID_H2Mb

        IF ( TrcID > 0 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(TrcID)%Emis%Val) ) THEN  ! Present! Read in the data
            ! Retrieve data into the HEMCO temporary!
            ! First, clear the buffer name. We are not sure if this was found yet.
            IF (.not. PRESENT( AltBuffer )) THEN
              LAST_TMP_REGRID_H2M  = "_HCO_Trc_TBD"
            ELSE
              LAST_TMP_REGRID_H2Mb = "_HCO_Trc_TBD"
            ENDIF

            ZBND = SIZE( HcoState%Spc(TrcID)%Emis%Val, 3 )
            TMP_HCO = 0.0_fp ! Clear the output first
            TMP_HCO(:,:,1:ZBND) = HcoState%Spc(TrcID)%Emis%Val(:,:,1:ZBND)

            ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValEmis/ImGrid: Read from HEMCO"

            ! Now perform the on-demand regrid
            CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL_target, ZBND )

            ! Now in TMP_MDL_target, read the pointer data
            ! FIXME: Could use a little DRY here (hplin, 6/6/20)
            Found = .true.
            Emis = TMP_MDL_target(I,J,L)

            IF (.not. PRESENT( AltBuffer )) THEN
              LAST_TMP_REGRID_H2M  = TMP_TrcIDFldName
            ELSE
              LAST_TMP_REGRID_H2Mb = TMP_TrcIDFldName
            ENDIF

            IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValEmis/ImGrid: Regrid OK return", TMP_TrcIDFldName, PRESENT( AltBuffer )
          ENDIF ! Associated
        ENDIF ! TrcID > 0
      ELSE ! Already in buffer! Just read the pointer data
        Found = .true.

        IF ( PRESENT(AltBuffer) ) THEN
          Emis  = TMP_MDLb(I,J,L)
        ELSE
          Emis  = TMP_MDL(I,J,L)
        ENDIF
      ENDIF
    ELSE
#endif
    ! Not GC-Classic or not on-demand intermediate grid, just shim around calls
      CALL GetHcoVal( HcoState, ExtState, TrcID, I, J, L, Found, Emis=Emis )
#ifdef MODEL_CLASSIC
    ENDIF

    TMP_MDL_target => NULL()
#endif

  END SUBROUTINE GetHcoValEmis
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoValDep
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
  SUBROUTINE GetHcoValDep ( Input_Opt, State_Grid, TrcID, I, J, L, Found, Dep )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoVal
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState

    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    INTEGER,            INTENT(IN   )  :: TrcID   ! GEOS-Chem tracer ID
    INTEGER,            INTENT(IN   )  :: I, J, L ! Position
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,            INTENT(  OUT)  :: Found   ! Was this tracer ID found?
    REAL(hp),           INTENT(  OUT)  :: Dep     ! Deposition [1/s]
!
! !REMARKS:
!  This subroutine is currently just a stub to call the equivalent in HEMCO
!  utilities. This assumes that HEMCO and GEOS-Chem grids match.
!  When the intermediate grid is implemented, a regridding routine will live
!  here and regrid data on-demand.
!
!  This also assumes that TrcID matches HEMCO tracer ID. If not, a mapping
!  needs to be performed here.
!
! !REVISION HISTORY:
!  20 Oct 2014 - C. Keller - Initial Version
!  12 Mar 2020 - H.P. Lin  - Now wrapper around common utilities
!  08 Jun 2020 - H.P. Lin  - Add GC-Classic on-demand regridding
!EOP
!------------------------------------------------------------------------------
!BOC
#ifdef MODEL_CLASSIC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=32)  :: TMP_TrcIDFldName            ! Temporary tracer field name _HCO_TrcD_<id>

    IF ( Input_Opt%LIMGRID ) THEN
      ! Check if we have to load the data.
      IF ( TrcID > 0 .and. (.not. ASSOCIATED(HcoState%Spc(TrcID)%Depv%Val)) ) RETURN

      ! ... on-demand intermediate regridding. Check if we already have this field
      Found = .false.
      WRITE(TMP_TrcIDFldName, '(a,i4)') "_HCO_TrcD_", TrcID

      IF( TMP_TrcIDFldName /= LAST_TMP_REGRID_H2Mb ) THEN   ! Not already in buffer
        ! Do not use GetHcoVal: load the entire chunk into memory

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValDep/ImGrid: Attempting to load", TMP_TrcIDFldName, "was", LAST_TMP_REGRID_H2Mb
        ! Note: TrcID matches HcoID here. If not, remap the tracer ID to HEMCO ID below.
        IF ( TrcID > 0 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(TrcID)%Depv%Val) ) THEN  ! Present! Read in the data

            ! Retrieve data into the HEMCO temporary!
            ! First, clear the buffer name. We are not sure if this was found yet.
            LAST_TMP_REGRID_H2Mb = "_HCO_TrcD_TBD"

            TMP_HCO = 0.0_fp ! Clear the output first
            TMP_HCO(:,:,1) = HcoState%Spc(TrcID)%Depv%Val(:,:)

            ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValDep/ImGrid: Read from HEMCO"

            ! Now perform the on-demand regrid
            ! Note deposition is surface layer, so L is always 1. Read ZBND = 1
            CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDLb, 1 )

            ! Now in TMP_MDLb, read the pointer data
            ! FIXME: Could use a little DRY here (hplin, 6/6/20)
            Found = .true.
            Dep = TMP_MDLb(I,J,1)
            LAST_TMP_REGRID_H2Mb = TMP_TrcIDFldName

            IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# GetHcoValDep/ImGrid: Regrid OK return", TMP_TrcIDFldName
          ENDIF
        ENDIF
      ELSE ! Already in buffer! Just read the pointer data
        Found = .true.
        Dep  = TMP_MDLb(I,J,1)
      ENDIF
    ELSE
#endif
    ! Not GC-Classic or not on-demand intermediate grid, just shim around calls
      CALL GetHcoVal( HcoState, ExtState, TrcID, I, J, L, Found, Dep=Dep )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

  END SUBROUTINE GetHcoValDep
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_EvalFld_3D
!
! !DESCRIPTION: Subroutine HCO\_GC\_EvalFld\_3D is a wrapper routine to obtain
!  the 3D data field belonging to the emissions list data.
!  It is a stub to simply map the call and route it to HEMCO in most cases,
!  and for GC-Classic with a different HEMCO grid, perform a transparent
!  HEMCO to model regrid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_EvalFld_3D ( Input_Opt, State_Grid, cName, Arr3D, RC, FOUND )
!
! !USES:
!
    USE HCO_Calc_Mod,         ONLY : HCO_EvalFld
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*), INTENT(IN   )  :: cName
!
! !OUTPUT ARGUMENTS:
!
    REAL(hp),         INTENT(INOUT)  :: Arr3D(:,:,:) ! 3D array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FOUND
!
! !REVISION HISTORY:
!  04 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                          :: FND
    CHARACTER(LEN=255)               :: ThisLoc
    CHARACTER(LEN=512)               :: ErrMsg

    INTEGER                          :: ZBND

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCO_GC_EvalFld_3D (in module GeosCore/hco_utilities_gc_mod.F90)'

    ! Empty the target first
    Arr3D   = 0.0_hp

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
      !=====================================================================
      ! We ARE USING the HEMCO intermediate grid
      !=====================================================================

      ! Sanity check - output array must be sized correctly for MODEL grid
      IF ( SIZE(Arr3D, 1) /= State_Grid%NX   .or.                            &
           SIZE(Arr3D, 2) /= State_Grid%NY ) THEN
        RC = GC_FAILURE
        ErrMsg = 'Input array dimensions are incorrect!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

      ! Z-level boundary for 3-D data and sanity check, 1st pass
      ZBND = MAX(1, SIZE(Arr3D, 3))
      IF( ZBND .ge. State_Grid%NZ ) THEN
        RC = GC_FAILURE
        ErrMsg = 'Input array Z-dimension higher than model maximum!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

      ! Check if cName is existing in the regrid buffer.
      !If not, regrid on-the-fly
      IF ( cName /= LAST_TMP_REGRID_H2M ) THEN
        !IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_3D: Last regrid not equal, looking up field ", cName

        ! Now retrieve data into the HEMCO temporary!
        ! The bdy is a slice to ensure safety
        CALL HCO_EvalFld( HcoState, cName, TMP_HCO(:,:,1:ZBND), RC, FND )

        ! If failure, return up the chain. The calls to this function will
        ! be able to propagate the error above.
        IF ( RC /= GC_SUCCESS ) RETURN

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_3D: Lookup complete", cName, FND

        IF ( FND ) THEN

          ! For safety, overwrite the temporary
          LAST_TMP_REGRID_H2M = "_HCO_Eval3D_TBD"

          ! Regrid the buffer appropriately.
          ! We do not use TMP_MDL here in EvalFld,
          ! because the field target is given.
          ! ( Input_Opt, State_Grid, State_Grid_HCO, PtrIn, PtrOut, ZBND )
          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, Arr3D, ZBND )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_3D: Regrid complete, ", cName, " z-boundary", ZBND

          ! The output should be in Arr3D and ready to go.
          LAST_TMP_REGRID_H2M = cName

        ENDIF

      ELSE
        ! Already existing in the buffer. Simply copy the data
        Arr3D(:,:,1:ZBND) = TMP_MDL(:,:,1:ZBND)
        FND = .true.

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_3D: Last regrid equal, reading from buffer"
      ENDIF
    ELSE
      !=====================================================================
      ! We ARE NOT USING the HEMCO intermediate grid
      !=====================================================================
#endif
      ! In which case, we just pass the call through
      CALL HCO_EvalFld( HcoState, cName, Arr3D, RC, FND )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

    IF( PRESENT(FOUND) ) THEN
      FOUND = FND
    ENDIF

  END SUBROUTINE HCO_GC_EvalFld_3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_EvalFld_2D
!
! !DESCRIPTION: Subroutine HCO\_GC\_EvalFld\_2D is a wrapper routine to obtain
!  the 3D data field belonging to the emissions list data.
!  It is a stub to simply map the call and route it to HEMCO in most cases,
!  and for GC-Classic with a different HEMCO grid, perform a transparent
!  HEMCO to model regrid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_EvalFld_2D ( Input_Opt, State_Grid, cName, Arr2D, RC, FOUND )
!
! !USES:
!
    USE HCO_Calc_Mod,         ONLY : HCO_EvalFld
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*), INTENT(IN   )  :: cName
!
! !OUTPUT ARGUMENTS:
!
    REAL(hp),         INTENT(INOUT)  :: Arr2D(:,:)   ! 2D array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FOUND
!
! !REVISION HISTORY:
!  04 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                          :: FND
    CHARACTER(LEN=255)               :: ThisLoc
    CHARACTER(LEN=512)               :: ErrMsg

    ! debug
    INTEGER                          :: II

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCO_GC_EvalFld_2D (in module GeosCore/hco_utilities_gc_mod.F90)'

    ! Empty the target first
    Arr2D   = 0.0_hp

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
      !=====================================================================
      ! We ARE USING the HEMCO intermediate grid
      !=====================================================================

      ! Sanity check - output array must be sized correctly for MODEL grid
      IF ( SIZE(Arr2D, 1) /= State_Grid%NX .or. SIZE(Arr2D, 2) /= State_Grid%NY ) THEN
        RC = GC_FAILURE
        ErrMsg = 'Input array dimensions are incorrect!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
      ENDIF

      ! Check if cName is existing in the regrid buffer. If not, regrid on-the-fly
      IF ( cName /= LAST_TMP_REGRID_H2M ) THEN
        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_2D: Last regrid not equal, looking up field ", cName

        ! Now retrieve data into the HEMCO temporary!
        ! The bdy is a slice to ensure safety
        CALL HCO_EvalFld( HcoState, cName, TMP_HCO(:,:,1), RC, FND )

        ! If failure, return up the chain. The calls to this function will
        ! be able to propagate the error above.
        IF ( RC /= GC_SUCCESS ) RETURN

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_2D: Lookup complete", cName, FND

        ! If not found, return
        IF ( FND ) THEN

          ! For safety, overwrite the temporary
          LAST_TMP_REGRID_H2M = "_HCO_Eval2D_TBD"

          ! Regrid the buffer appropriately. We do not use TMP_MDL here in EvalFld,
          ! because the field target is given.
          ! Z-boundary is 1 because 2-D field.
          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL, 1 )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_2D: Regrid", cName, " complete"

          ! Note that we cannot pass Arr2D to call above directly because it accepts a
          ! 3-D argument. So we regrid to the target dummy and copy
          Arr2D(:,:) = TMP_MDL(:,:,1)

          ! The output should be in Arr2D and ready to go.
          LAST_TMP_REGRID_H2M = cName

        ENDIF
      ELSE
        ! Already existing in the buffer. Simply copy the data
        Arr2D(:,:) = TMP_MDL(:,:,1)
        FND = .true.

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_EvalFld_2D: Last regrid equal, reading from buffer"
      ENDIF
    ELSE
      !=====================================================================
      ! We ARE NOT USING the HEMCO intermediate grid
      !=====================================================================
#endif
      ! In which case, we just pass the call through
      CALL HCO_EvalFld( HcoState, cName, Arr2D, RC, FND )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

    IF( PRESENT(FOUND) ) THEN
      FOUND = FND
    ENDIF

  END SUBROUTINE HCO_GC_EvalFld_2D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_GetPtr_2D
!
! !DESCRIPTION: Subroutine HCO\_GC\_GetPtr_2D is a wrapper routine to obtain
!  the 2D data pointer "directly" to HEMCO.
!  It is a stub to simply map the call and route it to HEMCO in most cases,
!  and for GC-Classic with a different HEMCO grid, perform a transparent
!  HEMCO to model regrid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_GetPtr_2D ( Input_Opt, State_Grid, DctName, Ptr2D, RC, &
                                TIDX, FOUND, FILLED )
!
! !USES:
!
    USE HCO_EmisList_Mod,     ONLY : HCO_GetPtr
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*), INTENT(IN   )  :: DctName    ! Container name
    INTEGER, OPTIONAL,INTENT(IN   )  :: TIDX       ! Time index (default = 1)
!
! !OUTPUT ARGUMENTS:
!
    REAL(sp), POINTER                :: Ptr2D(:,:)   ! Output array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FOUND
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FILLED
!
! !REMARKS:
!  Note that there is some code duplication here, because we have to handle
!  the optional arguments. Ideally we want to refactor away all calls to GetPtr,
!  but the met field reading requires a time index, so unfortunately we have to
!  replicate this here.
!
!  Note that for GC-Classic IMGrid, only ONE pointer may be kept at a time.
!  This is an underlying assumption that may bite and needs to be taken care of.
!  Once the data is received in the form of a pointer from HCO_GC_GetPtr,
!  you must COPY it and operate on it, otherwise you will operate in the shimmed
!  temporary. BE WARNED! (hplin, 6/8/20)
!
! !REVISION HISTORY:
!  08 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Internal copies of logicals to pass to the GetPtr call, because
    ! we need to replicate the previously optional arguments...
    LOGICAL                          :: iFOUND, iFILLED
    INTEGER                          :: iTIDX

    ! Debug
    INTEGER                          :: II, JJ

    CHARACTER(LEN=255)               :: ThisLoc
    CHARACTER(LEN=512)               :: ErrMsg
    CHARACTER(LEN=80)                :: TMP_GetPtrFldName

    ! Pointers
    REAL(sp), POINTER                :: TMP_Ptr2D(:,:)

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCO_GC_GetPtr_2D (in module GeosCore/hco_utilities_gc_mod.F90)'

    iTIDX = 1
    IF ( PRESENT(TIDX ) ) iTIDX = TIDX

    ! Nullify pointer
    TMP_Ptr2D => NULL()

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
      ! Build the name which requires a unique recognition of the time index,
      ! in case they are read sequentially
      write(TMP_GetPtrFldName, '(a,i4)') DctName, iTIDX

      ! Check if is existing in the regrid buffer. If not, regrid on-the-fly
      IF ( TMP_GetPtrFldName /= LAST_TMP_REGRID_H2M ) THEN
        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Last regrid not equal, looking up field ", TMP_GetPtrFldName

        ! For safety, overwrite the temporary
        LAST_TMP_REGRID_H2M = "_HCO_Ptr2D_TBD"

        ! Now retrieve data into the HEMCO temporary!
        CALL HCO_GetPtr( HcoState, DctName, TMP_Ptr2D, RC, iTIDX, iFOUND, iFILLED )

        ! If failure, return up the chain. The calls to this function will
        ! be able to propagate the error above.
        IF ( RC /= GC_SUCCESS ) RETURN

        ! If not found, return
        IF ( iFOUND .and. iFILLED ) THEN

          ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Lookup complete", TMP_GetPtrFldName, iFOUND

          ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Dim Debug", SIZE(TMP_HCO, 1), SIZE(TMP_HCO, 2), SIZE(TMP_Ptr2D, 1), SIZE(TMP_Ptr2D, 2)

          ! Copy data to the temporary
          TMP_HCO(:,:,1) = TMP_Ptr2D(:,:)

          ! Regrid the buffer appropriately. We do not use TMP_MDL here in EvalFld,
          ! because the field target is given.
          ! Z-boundary is 1 because 2-D field.
          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL, 1 )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Regrid", TMP_GetPtrFldName, "complete"

          ! Free the pointer
          TMP_Ptr2D => NULL()

          ! Copy to the r4 so you can downgrade the precision for GetPtr calls
          ! for backwards compatibility
          TMP_MDL_r4(:,:,1) = TMP_MDL(:,:,1)

          ! Point to target dummy
          Ptr2D => TMP_MDL_r4(:,:,1)

          ! The output should be pointing to Ptr2D (in TMP_MDL:,:,1) and ready to go.
          LAST_TMP_REGRID_H2M = TMP_GetPtrFldName

        ENDIF
      ELSE
        ! Already existing in the buffer. Simply point to the data
        Ptr2D => TMP_MDL_r4(:,:,1)
        iFOUND = .true.
        iFILLED = .true.

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Last regrid equal, pointing to buffer"
      ENDIF
    ELSE
      ! Not within the intermediate grid code path
#endif
      ! In which case, we just pass the call through
      CALL HCO_GetPtr( HcoState, DctName, Ptr2D, RC, iTIDX, iFOUND, iFILLED )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

    IF( PRESENT(FOUND) ) THEN
      FOUND = iFOUND
    ENDIF

    IF( PRESENT(FILLED) ) THEN
      FILLED = iFILLED
    ELSEIF ( .not. iFILLED ) THEN
      ! RC = GC_FAILURE
      ! ErrMsg = 'Could not fill last GetPtr container!'
      ! CALL GC_Error( ErrMsg, RC, ThisLoc )

      Ptr2D => NULL()
      ! FIXME: Maybe need to throw a HEMCO error from here. See behavior
      ! in HCO_EmisList_Mod
    ENDIF

  END SUBROUTINE HCO_GC_GetPtr_2D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_GetPtr_3D
!
! !DESCRIPTION: Subroutine HCO\_GC\_GetPtr_3D is a wrapper routine to obtain
!  the 3D data pointer "directly" to HEMCO.
!  It is a stub to simply map the call and route it to HEMCO in most cases,
!  and for GC-Classic with a different HEMCO grid, perform a transparent
!  HEMCO to model regrid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_GetPtr_3D ( Input_Opt, State_Grid, DctName, Ptr3D, RC, &
                                TIDX, FOUND, FILLED )
!
! !USES:
!
    USE HCO_EmisList_Mod,     ONLY : HCO_GetPtr
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*), INTENT(IN   )  :: DctName    ! Container name
    INTEGER, OPTIONAL,INTENT(IN   )  :: TIDX       ! Time index (default = 1)
!
! !OUTPUT ARGUMENTS:
!
    REAL(sp), POINTER                :: Ptr3D(:,:,:) ! Output array
    INTEGER,          INTENT(INOUT)  :: RC           ! Return code
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FOUND
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FILLED
!
! !REMARKS:
!  See 2D version.
!
! !REVISION HISTORY:
!  08 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Internal copies of logicals to pass to the GetPtr call, because
    ! we need to replicate the previously optional arguments...
    LOGICAL                          :: iFOUND, iFILLED
    INTEGER                          :: iTIDX
    INTEGER                          :: ZBND                     ! Maximum z-boundary

    CHARACTER(LEN=255)               :: ThisLoc
    CHARACTER(LEN=512)               :: ErrMsg
    CHARACTER(LEN=80)                :: TMP_GetPtrFldName

    ! Pointers
    REAL(sp), POINTER                :: TMP_Ptr3D(:,:,:)

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCO_GC_GetPtr_3D (in module GeosCore/hco_utilities_gc_mod.F90)'

    iTIDX = 1
    IF ( PRESENT(TIDX ) ) iTIDX = TIDX

    ! Nullify pointer
    TMP_Ptr3D => NULL()

#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN

      ! Build the name which requires a unique recognition of the time index,
      ! in case they are read sequentially
      write(TMP_GetPtrFldName, '(a,i4)') DctName, iTIDX

      ! Check if is existing in the regrid buffer. If not, regrid on-the-fly
      IF ( TMP_GetPtrFldName /= LAST_TMP_REGRID_H2M ) THEN
        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_3D: Last regrid not equal, looking up field ", TMP_GetPtrFldName

        ! For safety, overwrite the temporary
        LAST_TMP_REGRID_H2M = "_HCO_Ptr3D_TBD"

        ! Now retrieve data into the HEMCO temporary!
        CALL HCO_GetPtr( HcoState, DctName, TMP_Ptr3D, RC, iTIDX, iFOUND, iFILLED )

        ! If failure, return up the chain. The calls to this function will
        ! be able to propagate the error above.
        IF ( RC /= GC_SUCCESS ) RETURN

        ! If not found, return
        IF ( iFOUND .and. iFILLED ) THEN

          ! Get z-boundary
          ZBND = MAX(1, SIZE(TMP_Ptr3D, 3))

          ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_3D: Lookup complete", TMP_GetPtrFldName, iFOUND, ZBND

          ! Copy data to the temporary
          TMP_HCO(:,:,1:ZBND) = TMP_Ptr3D(:,:,1:ZBND)

          ! Regrid the buffer appropriately. We do not use TMP_MDL here in EvalFld,
          ! because the field target is given.
          ! Z-boundary is 1 because 2-D field.
          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL, ZBND )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_3D: Regrid", TMP_GetPtrFldName, "complete"

          ! Free the pointer
          TMP_Ptr3D => NULL()

          ! Copy to the r4 so you can downgrade the precision for GetPtr calls
          ! for backwards compatibility
          TMP_MDL_r4(:,:,1:ZBND) = TMP_MDL(:,:,1:ZBND)

          ! Point to target dummy
          Ptr3D => TMP_MDL_r4(:,:,1:ZBND)

          ! The output should be pointing to Ptr2D (in TMP_MDL:,:,1) and ready to go.
          LAST_TMP_REGRID_H2M = TMP_GetPtrFldName
          LAST_TMP_MDL_ZBND   = ZBND ! Remember the z-boundary ... will be used later

        ENDIF
      ELSE
        ! Already existing in the buffer. Simply point to the data
        Ptr3D => TMP_MDL_r4(:,:,1:LAST_TMP_MDL_ZBND)
        iFOUND = .true.
        iFILLED = .true.

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetPtr_2D: Last regrid equal, pointing to buffer"
      ENDIF
    ELSE
      ! Not within the intermediate grid code path
#endif
      ! In which case, we just pass the call through
      CALL HCO_GetPtr( HcoState, DctName, Ptr3D, RC, iTIDX, iFOUND, iFILLED )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

    IF( PRESENT(FOUND) ) THEN
      FOUND = iFOUND
    ENDIF

    IF( PRESENT(FILLED) ) THEN
      FILLED = iFILLED
    ELSEIF ( .not. iFILLED ) THEN
      ! RC = GC_FAILURE
      ! ErrMsg = 'Could not fill last GetPtr_3D container!'
      ! CALL GC_Error( ErrMsg, RC, ThisLoc )

      Ptr3D => NULL()
    ENDIF

  END SUBROUTINE HCO_GC_GetPtr_3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_GetDiagn_3D
!
! !DESCRIPTION: Subroutine HCO_GC_GetDiagn is a shim for the original GetHcoDiagn.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_GetDiagn_3D ( Input_Opt,      State_Grid, DiagnName, &
                                  StopIfNotFound, RC,         Ptr3D,     &
                                  COL,            AutoFill,   AltBuffer )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoDiagn
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState

    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*),   INTENT(IN   )  :: DiagnName  ! Name of diagnostics
    LOGICAL,            INTENT(IN   )  :: StopIfNotFound

    INTEGER, OPTIONAL,  INTENT(IN   )  :: COL        ! Collection Nr.
    INTEGER, OPTIONAL,  INTENT(IN   )  :: AutoFill   ! Autofill diagnostics only?

    LOGICAL, OPTIONAL,  INTENT(IN   )  :: AltBuffer  ! Alternate buffer? (Use B)
!
! !OUTPUT ARGUMENTS:
!
    REAL(sp),           POINTER        :: Ptr3D(:,:,:)
    INTEGER,            INTENT(INOUT)  :: RC
!
! !REMARKS:
!  See GetPtr. Note that this allows an alternative buffer to be used to
!  avoid conflict in reading the cached pointer.
!
!  FOR SAFETY SAKE, DESTROY THE POINTER READ FROM THIS ROUTINE BEFORE THE NEXT
!  CALL TO HCO_GC_GetDiagn OR YOU WILL HAVE WRONG DATA FED TO YOU!
!
! !REVISION HISTORY:
!  21 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: iCOL, iAF

#ifdef MODEL_CLASSIC
    CHARACTER(LEN=90)  :: TMP_DiagnFldName            ! Temporary diagnostic fld name
    INTEGER            :: ZBND
    REAL(sp), POINTER  :: TMP_Ptr3D(:,:,:)
    REAL(hp), POINTER  :: TMP_MDL_target(:,:,:)       ! Pointer to ease switcheroo of the model target buffer
    REAL(sp), POINTER  :: TMP_MDL_target4(:,:,:)      ! Pointer to ease switcheroo of the model target buffer
#endif

    ! Initialize local variables
    iCOL = HcoState%Diagn%HcoDiagnIDManual
    IF ( PRESENT(COL) ) THEN
      iCOL = COL
    ENDIF

    iAF = -1
    IF ( PRESENT(AutoFill) ) THEN
      iAF = AutoFill
    ENDIF

#ifdef MODEL_CLASSIC
    TMP_MDL_target => NULL()
    TMP_Ptr3D => NULL()

    IF ( Input_Opt%LIMGRID ) THEN
      ! The below section must be OMP CRITICAL because it is stateful.
      ! The first call to the critical section will update the container!!
      !$OMP CRITICAL

      IF ( PRESENT( AltBuffer ) ) THEN
        TMP_MDL_target => TMP_MDLb
        TMP_MDL_target4 => TMP_MDL_r4b
      ELSE
        TMP_MDL_target => TMP_MDL
        TMP_MDL_target4 => TMP_MDL_r4
      ENDIF

      ! ... on-demand intermediate regridding. Check if we already have this field
      WRITE(TMP_DiagnFldName, *) "_Dgn_", DiagnName
      IF( .not. ( (.not. PRESENT( AltBuffer ) .and. TMP_DiagnFldName == LAST_TMP_REGRID_H2M) .or. &
                  (      PRESENT( AltBuffer ) .and. TMP_DiagnFldName == LAST_TMP_REGRID_H2Mb) ) ) THEN   ! Not already in buffer

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_3D: Last regrid not equal, looking up field ", TMP_DiagnFldName

        ! Now retrieve data into the HEMCO temporary!
        ! Note that the data is in sp and has to be promoted for regridding,
        ! then demoted again for output.
        CALL GetHcoDiagn( HcoState, ExtState, DiagnName, StopIfNotFound, RC, &
                          Ptr3D=TMP_Ptr3D,    COL=iCOL,  AutoFill=iAF )

        ! If not found, return
        IF ( ASSOCIATED( TMP_Ptr3D ) ) THEN

          ! For safety, overwrite the temporary
          LAST_TMP_REGRID_H2M = "_HCO_Dgn3D_TBD"

          ! Get z-boundary
          ZBND = MAX(1, SIZE(TMP_Ptr3D, 3))

          ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_3D: Lookup complete", TMP_DiagnFldName, ZBND

          ! Copy data to the temporary
          TMP_HCO(:,:,1:ZBND) = TMP_Ptr3D(:,:,1:ZBND)

          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL, ZBND )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_3D: Regrid", TMP_DiagnFldName, "complete", PRESENT(AltBuffer)

          ! Free the pointer
          TMP_Ptr3D => NULL()

          ! Copy to the r4 so you can downgrade the precision for GetPtr calls
          ! for backwards compatibility
          TMP_MDL_target4(:,:,1:ZBND) = TMP_MDL(:,:,1:ZBND)

          ! Point to target dummy
          Ptr3D => TMP_MDL_target4(:,:,1:ZBND)

          ! The output should be pointing to Ptr2D (in TMP_MDL:,:,1) and ready to go.
          LAST_TMP_REGRID_H2M = TMP_DiagnFldName
          LAST_TMP_MDL_ZBND   = ZBND ! Remember the z-boundary ... will be used later

        ENDIF

      ELSE ! Already in buffer! Just read the pointer data
        Ptr3D => TMP_MDL_target4(:,:,1:LAST_TMP_MDL_ZBND)
        ! ... fill the ptr
      ENDIF
      !$OMP END CRITICAL
      ! End of LIMGRID OMP Critical section

    ELSE
#endif
      ! Not GC-Classic or not on-demand intermediate grid, just shim around calls
      CALL GetHcoDiagn( HcoState, ExtState, DiagnName, StopIfNotFound, RC, &
                        Ptr3D=Ptr3D,        COL=iCOL,  AutoFill=iAF )
#ifdef MODEL_CLASSIC
    ENDIF

    TMP_MDL_target => NULL()
    TMP_MDL_target4 => NULL()
#endif

  END SUBROUTINE HCO_GC_GetDiagn_3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_GetDiagn_2D
!
! !DESCRIPTION: Subroutine HCO_GC_GetDiagn is a shim for the original GetHcoDiagn.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GC_GetDiagn_2D ( Input_Opt,      State_Grid, DiagnName, &
                                  StopIfNotFound, RC,         Ptr2D,     &
                                  COL,            AutoFill,   AltBuffer )
!
! !USES:
!
    USE HCO_Interface_Common, ONLY : GetHcoDiagn
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState

    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
#ifdef MODEL_CLASSIC
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO
#endif
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),     INTENT(IN   )  :: State_Grid ! Grid State
    CHARACTER(LEN=*),   INTENT(IN   )  :: DiagnName  ! Name of diagnostics
    LOGICAL,            INTENT(IN   )  :: StopIfNotFound

    INTEGER, OPTIONAL,  INTENT(IN   )  :: COL        ! Collection Nr.
    INTEGER, OPTIONAL,  INTENT(IN   )  :: AutoFill   ! Autofill diagnostics only?

    LOGICAL, OPTIONAL,  INTENT(IN   )  :: AltBuffer  ! Alternate buffer? (Use B)
!
! !OUTPUT ARGUMENTS:
!
    REAL(sp),           POINTER        :: Ptr2D(:,:)
    INTEGER,            INTENT(INOUT)  :: RC
!
! !REMARKS:
!  See GetPtr. Note that this allows an alternative buffer to be used to
!  avoid conflict in reading the cached pointer.
!
!  FOR SAFETY SAKE, DESTROY THE POINTER READ FROM THIS ROUTINE BEFORE THE NEXT
!  CALL TO HCO_GC_GetDiagn OR YOU WILL HAVE WRONG DATA FED TO YOU!
!
! !REVISION HISTORY:
!  21 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: iCOL, iAF

#ifdef MODEL_CLASSIC
    CHARACTER(LEN=90)  :: TMP_DiagnFldName            ! Temporary diagnostic fld name
    REAL(sp), POINTER  :: TMP_Ptr2D(:,:)
    REAL(hp), POINTER  :: TMP_MDL_target(:,:,:)       ! Pointer to ease switcheroo of the model target buffer
    REAL(sp), POINTER  :: TMP_MDL_target4(:,:,:)      ! Pointer to ease switcheroo of the model target buffer
#endif

    ! initialize local variables
    iCOL = HcoState%Diagn%HcoDiagnIDManual
    IF ( PRESENT(COL) ) THEN
      iCOL = COL
    ENDIF

    iAF = -1
    IF ( PRESENT(AutoFill) ) THEN
      iAF = AutoFill
    ENDIF

#ifdef MODEL_CLASSIC
    TMP_MDL_target => NULL()
    TMP_Ptr2D => NULL()

    IF ( Input_Opt%LIMGRID ) THEN
      ! The below section must be OMP CRITICAL because it is stateful.
      ! The first call to the critical section will update the container!!
      !$OMP CRITICAL

      IF ( PRESENT( AltBuffer ) ) THEN
        TMP_MDL_target => TMP_MDLb
        TMP_MDL_target4 => TMP_MDL_r4b
      ELSE
        TMP_MDL_target => TMP_MDL
        TMP_MDL_target4 => TMP_MDL_r4
      ENDIF

      ! ... on-demand intermediate regridding. Check if we already have this field
      WRITE(TMP_DiagnFldName, *) "_Dgn_", DiagnName
      IF( .not. ( (.not. PRESENT( AltBuffer ) .and. TMP_DiagnFldName == LAST_TMP_REGRID_H2M) .or. &
                  (      PRESENT( AltBuffer ) .and. TMP_DiagnFldName == LAST_TMP_REGRID_H2Mb) ) ) THEN   ! Not already in buffer

        ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_2D: Last regrid not equal, looking up field ", TMP_DiagnFldName

        ! Now retrieve data into the HEMCO temporary!
        ! Note that the data is in sp and has to be promoted for regridding,
        ! then demoted again for output.
        CALL GetHcoDiagn( HcoState, ExtState, DiagnName, StopIfNotFound, RC, &
                          Ptr2D=TMP_Ptr2D,    COL=iCOL,  AutoFill=iAF )

        ! If not found, return
        IF ( ASSOCIATED( TMP_Ptr2D ) ) THEN

          ! For safety, overwrite the temporary
          LAST_TMP_REGRID_H2M = "_HCO_Dgn2D_TBD"

          ! IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_2D: Lookup complete", TMP_DiagnFldName

          ! Copy data to the temporary
          TMP_HCO(:,:,1) = TMP_Ptr2D(:,:)

          CALL Regrid_HCO2MDL( Input_Opt, State_Grid, State_Grid_HCO, TMP_HCO, TMP_MDL, 1 )

          IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) WRITE(6,*) "# HCO_GC_GetDiagn_2D: Regrid", TMP_DiagnFldName, "complete"

          ! Free the pointer
          TMP_Ptr2D => NULL()

          ! Copy to the r4 so you can downgrade the precision for GetDgn calls
          ! for backwards compatibility
          TMP_MDL_target4(:,:,1) = TMP_MDL(:,:,1)

          ! Point to target dummy
          Ptr2D => TMP_MDL_target4(:,:,1)

          ! The output should be pointing to Ptr2D (in TMP_MDL:,:,1) and ready to go.
          LAST_TMP_REGRID_H2M = TMP_DiagnFldName
          LAST_TMP_MDL_ZBND   = 1 ! Remember the z-boundary ... will be used later

        ENDIF

      ELSE ! Already in buffer! Just read the pointer data
        Ptr2D => TMP_MDL_target4(:,:,1)
        ! ... fill the ptr
      ENDIF
      !$OMP END CRITICAL
      ! End of LIMGRID OMP Critical section

    ELSE
#endif
      ! Not GC-Classic or not on-demand intermediate grid, just shim around calls
      CALL GetHcoDiagn( HcoState, ExtState, DiagnName, StopIfNotFound, RC, &
                        Ptr2D=Ptr2D,        COL=iCOL,  AutoFill=iAF )
#ifdef MODEL_CLASSIC
    ENDIF

    TMP_MDL_target => NULL()
    TMP_MDL_target4 => NULL()
#endif

  END SUBROUTINE HCO_GC_GetDiagn_2D
!EOC
#if defined ( MODEL_CLASSIC )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gc_restart
!
! !DESCRIPTION: Subroutine GET\_GC\_RESTART reads species concentrations
!  [mol/mol] from the GEOS-Chem restart file and uses them to initialize
!  species concentrations in [kg/kg dry]. If species data are missing from
!  the restart file, pre-configured background values are used. If using the
!  mercury simulation, additional restart data are read from file.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_GC_Restart( Input_Opt, State_Chm, State_Grid, &
                            State_Met, RC )
!
! !USES:
!
   USE CMN_SIZE_Mod,     ONLY  : NDUST
   USE ErrCode_Mod
   USE Error_Mod
   USE HCO_State_GC_Mod,  ONLY : HcoState
   USE PhysConstants,     ONLY : AIRMW
   USE Input_Opt_Mod,     ONLY : OptInput
   USE Species_Mod,       ONLY : Species, SpcConc
   USE State_Chm_Mod,     ONLY : ChmState
   USE State_Grid_Mod,    ONLY : GrdState
   USE State_Met_Mod,     ONLY : MetState
   USE Time_Mod,          ONLY : Expand_Date
   USE UnitConv_Mod,      ONLY : Convert_Spc_Units
#ifdef APM
   USE APM_Init_Mod,      ONLY : APMIDS
#endif
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
   TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!
!  09 Feb 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER                   :: I, J, L, M, N      ! lon, lat, lev, indexes
   LOGICAL                   :: FOUND              ! Found in restart file?
   CHARACTER(LEN=60)         :: Prefix             ! utility string
   CHARACTER(LEN=255)        :: LOC                ! routine location
   CHARACTER(LEN=255)        :: MSG                ! message
   CHARACTER(LEN=255)        :: v_name             ! variable name
   REAL(fp)                  :: MW_g               ! species molecular weight
   REAL(fp)                  :: SMALL_NUM          ! small number threshold
   CHARACTER(LEN=63)         :: OrigUnit

   ! Temporary arrays and pointers
   REAL*4,  TARGET           :: Temp2D(State_Grid%NX,State_Grid%NY)
   REAL*4,  TARGET           :: Temp3D(State_Grid%NX,State_Grid%NY, &
                                       State_Grid%NZ)
   REAL*4,  POINTER          :: Ptr2D(:,:  )
   REAL*4,  POINTER          :: Ptr3D(:,:,:)

   ! For Hg simulation
   INTEGER                   :: Num_Hg_Categories
   INTEGER                   :: Total_Hg_Id
   CHARACTER(LEN=60)         :: HgSpc
   CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

   ! Default background concentration
   REAL(fp)                  :: Background_VV

   ! Objects
   TYPE(SpcConc),    POINTER :: Spc(:)
   TYPE(Species),    POINTER :: SpcInfo

   !=================================================================
   ! READ_GC_RESTART begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()
   SpcInfo     => NULL()
   Hg_Cat_Name => NULL()

   ! Name of this routine
   LOC = ' -> at Get_GC_Restart (in GeosCore/hco_utilities_gc_mod.F90)'

   ! Set minimum value threshold for [mol/mol]
   SMALL_NUM = 1.0e-30_fp

   ! Set pointer to species concentrations
   Spc => State_Chm%Species

   !=================================================================
   ! Open GEOS-Chem restart file
   !=================================================================

   ! Write read message to log
   WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
   WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'

   !=================================================================
   ! Read species concentrations from NetCDF or use default
   ! background [mol/mol]; store in State_Chm%Species%Conc in [kg/kg dry]
   !=================================================================

   ! Print header for min/max concentration to log
   WRITE( 6, 110 )
110 FORMAT( 'Min and Max of each species in restart file [mol/mol]:' )

   ! Loop over species
   DO N = 1, State_Chm%nSpecies

      ! Initialize species concentration to all zeroes
      Spc(N)%Conc = 0.e+0_fp

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%MW_g

      ! Define variable name
      v_name = 'SPC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM(v_name), &
                          Ptr3D,     RC,         FOUND=FOUND )

      ! Check if species data is in file
      IF ( FOUND ) THEN
         SpcInfo%Is_InRestart = .TRUE.
      ELSE
         SpcInfo%Is_InRestart = .FALSE.
      ENDIF

      ! If data is in file, read in as [mol/mol] and convert to
      ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
      ! either stored in species database (advected species all levels and
      ! non-advected species levels in the chemistry grid) or a small number
      ! (non-advected species levels above the chemistry grid) converted to
      ! [kg/kg dry]
      IF ( SpcInfo%Is_InRestart ) THEN

         ! Print the min & max of each species as it is read from
         ! the restart file in mol/mol
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                            MINVAL( Ptr3D ), MAXVAL( Ptr3D ), SUM ( Ptr3D(:,:,1:State_Grid%NZ) )
120         FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                    '  Max = ',es15.9, '  Sum = ',es15.9)
         ENDIF

         ! Convert file value [mol/mol] to [kg/kg dry] for storage
         !$OMP PARALLEL DO       &
         !$OMP DEFAULT( SHARED ) &
         !$OMP PRIVATE( I, J, L )
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX
            Spc(N)%Conc(I,J,L) = Ptr3D(I,J,L) * MW_g / AIRMW
         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

      ELSE

         ! Set species to the background value converted to [kg/kg dry]
         !$OMP PARALLEL DO       &
         !$OMP DEFAULT( SHARED ) &
         !$OMP PRIVATE( I, J, L )
         ! Loop over all grid boxes
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX

            ! For non-advected species at levels above chemistry grid,
            ! use a small number for background
            IF ( L > State_Grid%MaxChemLev .and. &
                     .NOT. SpcInfo%Is_Advected ) THEN

               Spc(N)%Conc(I,J,L) = SMALL_NUM * MW_g / AIRMW

            ! For all other cases, use the background value
            ! stored in the species database
            ELSE

               Spc(N)%Conc(I,J,L) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

               ! Print to log if debugging is on
               IF ( Input_Opt%amIRoot .AND. &
                    I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                  WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
140               FORMAT('Species ', i3, ', ', a9, &
                         ': Use background = ', es15.9)
               ENDIF


            ENDIF

         ENDDO
         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

#ifdef APM
         !================================================================
         ! APM MICROPHYSICS
         !================================================================
         WRITE(*,*)'APM run does not find '// TRIM( SpcInfo%Name ),N
         IF(SpcInfo%Name(1:9)=='APMSPBIN2')THEN
            !$OMP PARALLEL DO        &
            !$OMP DEFAULT( SHARED  ) &
            !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               Spc(N)%Conc(I,J,L) = &
               Spc(APMIDS%id_SO4)%Conc(I,J,L)/20.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:9)=='APMSPBIN3')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               Spc(N)%Conc(I,J,L) = &
               Spc(APMIDS%id_SO4)%Conc(I,J,L)/20.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           !GanLuotest
           IF(SpcInfo%Name(1:10)=='APMSEABIN0')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               Spc(N)%Conc(I,J,L) = &
               Spc(APMIDS%id_SALA)%Conc(I,J,L)/9.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMSEABIN1')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               Spc(N)%Conc(I,J,L) = &
               Spc(APMIDS%id_SALC)%Conc(I,J,L)/10.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMDSTBIN1')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               Spc(N)%Conc(I,J,L) = &
                  ( Spc(APMIDS%id_DST1)%Conc(I,J,L)    &
                  + Spc(APMIDS%id_DST2)%Conc(I,J,L)    &
                  + Spc(APMIDS%id_DST3)%Conc(I,J,L)    &
                  + Spc(APMIDS%id_DST4)%Conc(I,J,L) )/6.D0
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMBCBIN')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               Spc(N)%Conc(I,J,L) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMOCBIN')THEN
              !$OMP PARALLEL DO        &
              !$OMP DEFAULT( SHARED  ) &
              !$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               Spc(N)%Conc(I,J,L) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
           ENDIF
#endif
      ENDIF

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Set species units
   State_Chm%Spc_Units = 'kg/kg dry'

   ! If in debug mode, print out species min and max in [molec/cm3]
   IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN

      ! Convert units
      PRINT *, " "
      PRINT *, "Species min and max in molec/cm3"
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              'molec/cm3', RC, OrigUnit=OrigUnit )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #1!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

      ! Print values
      DO N = 1, State_Chm%nSpecies
         SpcInfo => State_Chm%SpcData(N)%Info
         WRITE(6,150) N, TRIM( SpcInfo%Name ),         &
                         MINVAL( Spc(N)%Conc(:,:,:) ), &
                         MAXVAL( Spc(N)%Conc(:,:,:) )
150      FORMAT( 'Species ', i3, ', ', a9,             &
                 ': Min = ', es15.9, ', Max = ', es15.9 )
         SpcInfo => NULL()
      ENDDO

      ! Convert units back
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              OrigUnit,  RC )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #2!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

   ENDIF

   !=========================================================================
   ! Get variables for KPP mechanisms (right now just fullchem and Hg)
   !=========================================================================
   IF ( ( Input_Opt%ITS_A_FULLCHEM_SIM .or.                                  &
          Input_Opt%ITS_A_MERCURY_SIM        ) .and. Input_Opt%LCHEM ) THEN

      !----------------------------------------------------------------------
      ! KPP_HVALUE (count of internal timesteps at each grid box)
      !----------------------------------------------------------------------
      v_name = 'KPP_HVALUE'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr3D,     RC,         FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%KPPHvalue = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name              ),                  &
                            MINVAL(  State_Chm%KPPHvalue ),                  &
                            MAXVAL(  State_Chm%KPPHvalue ),                  &
                            SUM(     State_Chm%KPPHvalue )
         ENDIF
      ELSE
         State_Chm%KPPHvalue = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! FORMAT strings
500   FORMAT( a                                                              )
510   FORMAT( a21, ': Min = ', es15.9, '  Max = ', es15.9, '  Sum = ',es15.9 )
520   FORMAT( a21, ': not found in restart, set to zero'                      )

   ENDIF

   !=========================================================================
   ! Get variables for Soil NOx emissions
   !=========================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !----------------------------------------------------------------------
      ! WETDEP_N (wet-deposited nitrogen)
      !----------------------------------------------------------------------
      v_name = 'WETDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr2D,     RC,         FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%WetDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name                   ),             &
                            MINVAL(  State_Chm%WetDepNitrogen ),             &
                            MAXVAL(  State_Chm%WetDepNitrogen ),             &
                            SUM(     State_Chm%WetDepNitrogen )
         ENDIF
      ELSE
         State_Chm%WetDepNitrogen = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) TRIM( v_name )
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !----------------------------------------------------------------------
      ! DRYDEP_N (dry-deposited nitrogen)
      !----------------------------------------------------------------------
      v_name = 'DRYDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr2D,     RC,         FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%DryDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name                   ),             &
                            MINVAL(  State_Chm%DryDepNitrogen ),             &
                            MAXVAL(  State_Chm%DryDepNitrogen ),             &
                            SUM(     State_Chm%DryDepNitrogen )
         ENDIF
      ELSE
         State_Chm%DryDepNitrogen = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

   ENDIF

   !=========================================================================
   ! Read variables for sulfate chemistry and aerosols
   !=========================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
        Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      !----------------------------------------------------------------------
      ! H2O2_AFTERCHEM
      !----------------------------------------------------------------------
      v_name = 'H2O2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr3D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%H2O2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name                  ),              &
                            MINVAL(  State_Chm%H2O2AfterChem ),              &
                            MAXVAL(  State_Chm%H2O2AfterChem ),              &
                            SUM(     State_Chm%H2O2AfterChem )
        ENDIF
      ELSE
         State_Chm%H2O2AfterChem = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !----------------------------------------------------------------------
      ! SO2_AFTERCHEM
      !----------------------------------------------------------------------
      v_name = 'SO2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr3D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%SO2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name                 ),               &
                            MINVAL(  State_Chm%SO2AfterChem ),               &
                            MAXVAL(  State_Chm%SO2AfterChem ),               &
                            SUM(     State_Chm%SO2AfterChem )
         ENDIF
      ELSE
         State_Chm%SO2AfterChem = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !----------------------------------------------------------------------
      ! AeroH2O_SNA
      !----------------------------------------------------------------------
      v_name = 'AEROH2O_SNA'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr3D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%AeroH2O(:,:,:,NDUST+1) = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name                           ),     &
                            MINVAL(  State_Chm%AeroH2O(:,:,:,NDUST+1) ),     &
                            MAXVAL(  State_Chm%AeroH2O(:,:,:,NDUST+1) ),     &
                            SUM(     State_Chm%AeroH2O(:,:,:,NDUST+1) )
         ENDIF
      ELSE
         State_Chm%AeroH2O(:,:,:,NDUST+1) = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !----------------------------------------------------------------------
      ! ORVCsesq
      !----------------------------------------------------------------------
      IF ( Input_Opt%LCARB .AND. Input_Opt%LSOA ) THEN

         v_name = 'ORVCSESQ'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                             Ptr3D,     RC,         FOUND=FOUND                )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Chm%ORVCsesq(:,:,:) = Ptr3D
            IF ( Input_Opt%amIRoot ) THEN
               WRITE( 6, 510 ) ADJUSTL( v_name                           ),     &
                               MINVAL(  State_Chm%ORVCsesq(:,:,:) ),     &
                               MAXVAL(  State_Chm%ORVCsesq(:,:,:) ),     &
                               SUM(     State_Chm%ORVCsesq(:,:,:) )
            ENDIF
         ELSE
            State_Chm%ORVCsesq(:,:,:) = 0.0_fp
            IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
         ENDIF

         ! Nullify pointer
         Ptr3D => NULL()

      ENDIF

   ENDIF

   !=========================================================================
   ! Read variables for UCX and the HEMCO PARANOx extension
   !=========================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      !----------------------------------------------------------------------
      ! STATE_PSC (needed to initialize UCX)
      !----------------------------------------------------------------------
      v_name = 'STATE_PSC'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr3D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%STATE_PSC = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name              ),                  &
                            MINVAL(  State_Chm%STATE_PSC ),                  &
                            MAXVAL(  State_Chm%STATE_PSC ),                  &
                            SUM(     State_Chm%STATE_PSC )

         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
#ifdef ESMF_
            ! ExtData and HEMCO behave ambiguously - if the file was found
            ! but was full of zeros throughout the domain of interest, it
            ! will result in the same output from ExtData as if the field
            ! was missing from the file. As such, HEMCO cannot distinguish
            ! between a missing file and a field of zeros
            WRITE(6,*) 'PSC restart either all zeros in the '
            WRITE(6,*) 'root domain, or the restart file did '
            WRITE(6,*) 'not contain STATE_PSC. Root domain '
            WRITE(6,*) 'will be initialized PSC-free'
         ENDIF
#else
            WRITE( 6, 500 ) &
               'STATE_PSC not found in restart, initialize PSC-free'
         ENDIF
#endif
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !----------------------------------------------------------------------
      ! JOH (needed to initialize PARANOx)
      !----------------------------------------------------------------------
      v_name = 'JOH'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr2D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%JOH = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name       ),                         &
                            MINVAL(  State_Chm%JOH ),                        &
                            MAXVAL(  State_Chm%JOH ),                        &
                            SUM(     State_Chm%JOH )
         ENDIF
      ELSE
         State_Chm%JOH = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !----------------------------------------------------------------------
      ! JNO2 (needed to initialize PARANOx)
      !----------------------------------------------------------------------
      v_name = 'JNO2'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),             &
                          Ptr2D,     RC,         FOUND=FOUND                )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%JNO2 = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 510 ) ADJUSTL( v_name         ),                       &
                            MINVAL(  State_Chm%JNO2 ),                       &
                            MAXVAL(  State_Chm%JNO2 ),                       &
                            SUM(     State_Chm%JNO2 )
         ENDIF
      ELSE
         State_Chm%JNO2 = 0.0_fp
         IF ( Input_Opt%amIRoot ) WRITE( 6, 520 ) ADJUSTL( v_name )
      ENDIF
      ! Nullify pointer
      Ptr2D => NULL()

   ENDIF

   !=========================================================================
   ! Read ocean mercury variables
   !=========================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Print total mass to log
      WRITE( 6, 220 )
220   FORMAT(/, 'Total mass of each ocean and snow Hg species:')

      !----------------------------------------------------------------------
      ! Total Hg in ocean
      !----------------------------------------------------------------------
      DO M = 1, 3

         ! Define variable name
         SELECT CASE( M )
           CASE ( 1 )
              HgSpc    = 'Hg0'
           CASE ( 2 )
              HgSpc    = 'Hg2'
           CASE ( 3 )
              HgSpc    = 'HgP'
         END SELECT
         v_name = 'OCEAN_' // TRIM( HgSpc )

         ! Get variable from HEMCO and store in local array
         CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),          &
                             Ptr2D,     RC,         FOUND=FOUND             )

         ! Check if variable is in file
         IF ( FOUND ) THEN

            ! Check for negative concentrations (jaf, 7/6/11)
            DO I = 1, State_Grid%NX
            DO J = 1, State_Grid%NY
               IF ( Ptr2D(I,J) < 0.0d4 ) THEN
                  Ptr2D(I,J) = 0.0d4
               ENDIF
            ENDDO
            ENDDO

            ! Assign ocean mercury data and write total mass to log file
            SELECT CASE( M )
               CASE ( 1 )
                  State_Chm%OceanHg0 = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name            ),                 &
                                  SUM( State_Chm%OceanHg0 ), 'kg'
               CASE ( 2 )
                  State_Chm%OceanHg2 = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name            ),                 &
                                  SUM( State_Chm%OceanHg2 ), 'kg'
               CASE ( 3 )
                  State_Chm%OceanHgP = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name            ),                 &
                                  SUM( State_Chm%OceanHgP ), 'kg'
            END SELECT

         ELSE
            WRITE( 6, 230 ) TRIM( v_name )
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()

      ENDDO

      !--------------------------------------------------------------
      ! Hg snowpack on land and ocean
      !--------------------------------------------------------------
      DO M = 1, 4

         ! Define variable name prefix
         SELECT CASE( M )
            CASE ( 1 )
               Prefix = 'SNOW_HG_OCEAN'        ! Reducible on ocean
            CASE ( 2 )
               Prefix = 'SNOW_HG_OCEAN_STORED' ! Non-reducible on ocean
            CASE ( 3 )
               Prefix = 'SNOW_HG_LAND'         ! Reducible on land
            CASE ( 4 )
               Prefix = 'SNOW_HG_LAND_STORED'  ! Non-reducible on land
         END SELECT

         v_name = TRIM( Prefix )

         ! Get variable from HEMCO and store in local array
         CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM( v_name ),          &
                             Ptr2D,     RC,         FOUND=FOUND             )

         ! Check if variable is in file
         IF ( FOUND ) THEN

            ! Assign ocean mercury data and write total mass to file
            SELECT CASE( M )
               CASE ( 1 )
                  State_Chm%SnowHgOcean = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name                     ),        &
                                  SUM( State_Chm%SnowHgOcean       ), 'kg'
               CASE ( 2 )
                  State_Chm%SnowHgOceanStored = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name                     ),        &
                                  SUM( State_Chm%SnowHgOceanStored ),'kg'
               CASE ( 3 )
                  State_Chm%SnowHgLand = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name                     ),        &
                                  SUM( State_Chm%SnowHgLand        ), 'kg'
               CASE ( 4 )
                  State_Chm%SnowHgLandStored = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name                     ),        &
                                  SUM( State_Chm%SnowHgLandStored  ), 'kg'
               END SELECT

         ELSE
            WRITE( 6, 230 ) TRIM( v_name )
         ENDIF

      ENDDO

      ! Format strings
230   FORMAT( a24, ' not found in restart file, set to zero')
240   FORMAT( a24, ':   ', es15.9, 1x, a4)

      ! Print note that variables are initialized to zero if not
      ! found (currently only happens in tagged Hg simulation)
      IF ( Input_Opt%LSPLIT ) THEN
         WRITE( 6, 250 )
250      FORMAT( /, 'NOTE: all variables not found in restart ', &
                    'are initialized to zero')
      ENDIF

      ! Free pointers for Hg indexing
      Hg_Cat_Name => NULL()

   ENDIF

   !=================================================================
   ! Clean up
   !=================================================================

   ! Free pointer
   Spc => NULL()

   ! Mark end of section in log
   IF ( Input_Opt%LPRT .AND. Input_Opt%amIRoot ) THEN
      CALL DEBUG_MSG('### DONE GET_GC_RESTART')
   ENDIF
   WRITE( 6, '(a)' ) REPEAT( '=', 79 )

 END SUBROUTINE Get_GC_Restart
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_boundary_conditions
!
! !DESCRIPTION: Subroutine GET\_BOUNDARY\_CONDITIONS calls the various routines
! to get boundary conditions from HEMCO for nested grid simulations.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Boundary_Conditions( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, YYYYMMDD,  HHMMSS, RC )
!
! ! USES:
!
   USE ErrCode_Mod
   USE HCO_State_GC_Mod, ONLY : HcoState
   USE Input_Opt_Mod,    ONLY : OptInput
   USE PhysConstants,    ONLY : AIRMW
   USE Species_Mod,      ONLY : Species, SpcConc
   USE State_Chm_Mod,    ONLY : ChmState
   USE State_Grid_Mod,   ONLY : GrdState
   USE State_Met_Mod,    ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid State
   INTEGER,          INTENT(IN   )          :: YYYYMMDD   ! GMT date
   INTEGER,          INTENT(IN   )          :: HHMMSS     ! GMT time
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Meteorology State
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry State
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Apr 2019 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, N, NA     ! lon, lat, lev, spc indexes
   INTEGER              :: t_index            ! Time index
   LOGICAL              :: FOUND              ! Found in restart file?
   LOGICAL, SAVE        :: FIRST = .TRUE.     ! Is this the first routine call?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message
   CHARACTER(LEN=255)   :: v_name             ! variable name
   REAL(fp)             :: MW_g               ! species molecular weight
   CHARACTER(LEN=16)    :: STAMP

   ! Temporary arrays and pointers
   REAL*4,  TARGET        :: Temp3D(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
   REAL*4,  POINTER       :: Ptr3D(:,:,:)
   TYPE(SpcConc), POINTER :: Spc(:)

   ! Objects
   TYPE(Species), POINTER :: SpcInfo

   !=================================================================
   ! GET_BOUNDARY_CONDITIONS begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! We only need to get boundary conditions if this is a nested-grid
   ! simulation.  Otherwise the BoundaryCond field won't be allocated.
   IF ( .not. State_Grid%NestedGrid ) RETURN

   ! Initialize pointers
   Ptr3D     => NULL()
   SpcInfo   => NULL()

   ! Point to species array [kg/kg]
   Spc       => State_Chm%Species

   ! Name of this routine
   LOC = ' -> at Get_Boundary_Conditions (in GeosCore/hco_utilities_gc_mod.F90)'

   ! Find the proper time-slice to read from disk
   t_index = ( HHMMSS / 030000 ) + 1

   ! Stop w/ error if the time index is invalid
   IF ( t_index < 1 .or. t_index > 8 ) THEN
      WRITE( MSG, 100 ) t_index
100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
      CALL GC_Error( MSG, RC, LOC)
      RETURN
   ENDIF

   !=================================================================
   ! Read species concentrations from NetCDF [mol/mol] and
   ! store in State_Chm%BoundaryCond in [kg/kg dry]
   !=================================================================

   ! Print header for min/max concentration to log
   IF ( Input_Opt%amIRoot ) THEN
      WRITE( 6, 110 )
110   FORMAT( 'Min and Max of each species in BC file [mol/mol]:' )
   ENDIF

   ! Initialize BCs to all zeroes
   State_Chm%BoundaryCond = 0.e+0_fp

   ! Loop over advected species
   DO NA = 1, State_Chm%nAdvect

      ! Get the species ID from the advected species ID
      N = State_Chm%Map_Advect(NA)

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%MW_g

      ! Define variable name
      v_name = 'BC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM(v_name), Ptr3D, RC, &
                       TIDX=t_index, FOUND=FOUND )

      ! Check if BCs are found
      IF ( FOUND ) THEN

         ! Print the min & max of each species as it is read from
         ! the BC file in mol/mol if debug is turned on in geoschem_config.yml
         IF ( Input_Opt%amIRoot ) THEN
            IF ( FIRST .or. Input_Opt%LPRT ) THEN
               WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                               MINVAL( Ptr3D ), MAXVAL( Ptr3D )
120            FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                       '  Max = ',es15.9)
            ENDIF
         ENDIF

         ! Copy data from file to State_Chm%BoundaryCond
         ! and convert from [mol/mol] to [kg/kg dry]
         State_Chm%BoundaryCond(:,:,:,N) = Ptr3D(:,:,:) * MW_g / AIRMW

         ! Debug
         ! Print*, 'BCs found for ', TRIM( SpcInfo%Name ), &
         !         MINVAL(State_Chm%BoundaryCond(:,:,:,N)), &
         !         MAXVAL(State_Chm%BoundaryCond(:,:,:,N)), &
         !         SUM(State_Chm%BoundaryCond(:,:,:,N))

      ELSE

         ! Print to log if debug is turned on in geoschem_config.yml
         IF ( Input_Opt%amIRoot ) THEN
            IF ( FIRST .or. Input_Opt%LPRT ) THEN
               WRITE( 6, 130 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
130            FORMAT('Species ', i3, ', ', a9, ': Use background = ', es15.9)
            ENDIF
         ENDIF

         ! Use the background value stored in the species database
         State_Chm%BoundaryCond(:,:,:,N) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

      ENDIF


      ! Loop over grid boxes and apply BCs to the specified buffer zone
      !$OMP PARALLEL DO       &
      !$OMP DEFAULT( SHARED ) &
      !$OMP PRIVATE( I, J, L )
      DO L = 1, State_Grid%NZ

         ! First loop over all latitudes of the nested domain
         DO J = 1, State_Grid%NY

            ! West BC
            DO I = 1, State_Grid%WestBuffer
               Spc(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

            ! East BC
            DO I = (State_Grid%NX-State_Grid%EastBuffer)+1, State_Grid%NX
               Spc(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

         ENDDO

         ! Then loop over the longitudes of the nested domain
         DO I = 1+State_Grid%WestBuffer,(State_Grid%NX-State_Grid%EastBuffer)

            ! South BC
            DO J = 1, State_Grid%SouthBuffer
               Spc(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO

            ! North BC
            DO J = (State_Grid%NY-State_Grid%NorthBuffer)+1, State_Grid%NY
               Spc(N)%Conc(I,J,L) = State_Chm%BoundaryCond(I,J,L,N)
            ENDDO
         ENDDO

      ENDDO
      !OMP END PARALLEL DO

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Free pointer
   Spc => NULL()

   ! Reset FIRST flag
   FIRST = .FALSE.

   ! Echo output
   IF ( Input_Opt%amIRoot ) THEN
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 140 ) STAMP
140   FORMAT( 'GET_BOUNDARY_CONDITIONS: Done applying BCs at ', a )
   ENDIF

 END SUBROUTINE Get_Boundary_Conditions
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Regrid_HCO2MDL
!
! !DESCRIPTION: Subroutine Regrid\_HCO2MDL is a buffer function to regrid a given
!  field described on the HEMCO intermediate grid ("IMGrid") to the GEOS-Chem model
!  grid. Only horizontal interpolation is performed via Regrid\_A2A.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Regrid_HCO2MDL ( Input_Opt, State_Grid, State_Grid_HCO, &
                              PtrIn,     PtrOut,     ZBND, Debug )
!
! !USES:
!
    USE State_Grid_Mod,    ONLY : GrdState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE Regrid_A2A_Mod,    ONLY : Map_A2A
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )     :: Input_Opt        ! Input opts
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid       ! Grid state
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid_HCO   ! Optional, HEMCO intermediate
    REAL(hp),         INTENT(IN   )     :: PtrIn (:,:,:)    ! 3-D input data
    INTEGER, OPTIONAL,INTENT(IN   )     :: ZBND             ! z-level bounds of data
    LOGICAL, OPTIONAL,INTENT(IN   )     :: Debug            ! For debugging
!
! !OUTPUT ARGUMENTS:
!
    REAL(hp),         INTENT(  OUT)     :: PtrOut(:,:,:)    ! 3-D output data
!
! !REMARKS:
!  Usually, the regridded quantities are stored in the array temporaries in this module,
!  so PtrOut should be a pointer array temporary here.
!  This module is NOT multiple domain safe and should be kept to MODEL\_CLASSIC only.
!
!  For some cases, the PtrOut may point to a pointer array allocated somewhere else.
!
!  If ZBND is optionally specified, the regridding operation is capped up to the
!  given bound. This is usually used for 2-D data but may also limit useless regridding
!  beyond the pointer size. The array temporaries in the module are usually always
!  allocated to the model top to save memory thrashing, but we also do not want to
!  sacrifice compute.
!
! !REVISION HISTORY:
!  03 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                             :: NZ, L
    LOGICAL                             :: DBG

    IF ( PRESENT( ZBND ) ) THEN
      NZ = MIN(ZBND, State_Grid%NZ+1)                ! May be equal to NZ+1 maximum
    ELSE
      NZ = State_Grid%NZ
    ENDIF

    IF ( PRESENT( Debug ) ) THEN
      DBG = Debug
    ELSE
      DBG = .false.
    ENDIF

    ! Intermediate grid functionality?
    IF ( .not. Input_Opt%LIMGRID ) RETURN

    ! Initialize shadow regridding handles if they are not ready
    IF ( .not. ASSOCIATED( LonEdgeM ) ) THEN
      CALL Init_IMGrid ( Input_Opt, State_Grid, State_Grid_HCO )
    ENDIF

    ! Empty the target first
    PtrOut(:,:,1:NZ) = 0.0_hp

    ! Do the regridding layer by layer
    DO L = 1, NZ
      call Map_A2A (                             &
        State_Grid_HCO%NX, State_Grid_HCO%NY,    &   ! Input grid dimensions
        LonEdgeH,          LatEdgeH,             &   ! Lons and lat sines at input edges
        PtrIn(:,:,L),                            &   ! Input data, 2-D slice
        State_Grid%NX,     State_Grid%NY,        &   ! Output grid dimensions
        LonEdgeM,          LatEdgeM,             &
        PtrOut(:,:,L),                           &
        0, 0,                                    &   ! Pole treatment, scalar quantity
        missval=-1.e31_hp)                           ! Important to prevent crash
    ENDDO

    ! Return happily ever after

  END SUBROUTINE Regrid_HCO2MDL
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Regrid_MDL2HCO
!
! !DESCRIPTION: Subroutine Regrid\_MDL2HCO is a buffer function to regrid a given
!  field described on the GEOS-Chem model grid to the HEMCO intermediate grid ("IMGrid").
!  Only horizontal interpolation is performed via Regrid\_A2A.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Regrid_MDL2HCO ( Input_Opt, State_Grid, State_Grid_HCO, PtrIn, PtrOut, ZBND, &
                              ResetRegrName )
!
! !USES:
!
    USE State_Grid_Mod,    ONLY : GrdState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE Regrid_A2A_Mod,    ONLY : Map_A2A
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )     :: Input_Opt        ! Input opts
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid       ! Grid state
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid_HCO   ! Optional, HEMCO intermediate
    REAL(hp),         INTENT(IN   )     :: PtrIn (:,:,:)    ! 3-D input data
    INTEGER, OPTIONAL,INTENT(IN   )     :: ZBND             ! z-level bounds of data

    LOGICAL, OPTIONAL,INTENT(IN   )     :: ResetRegrName    ! Reset regridding name?
                                                            ! Used when called by outside routine
!
! !OUTPUT ARGUMENTS:
!
    REAL(hp),         INTENT(  OUT)     :: PtrOut(:,:,:)    ! 3-D output data
!
! !REMARKS:
!  Usually, the regridded quantities are stored in the array temporaries in this module,
!  so PtrOut should be a pointer array temporary here.
!  This module is NOT multiple domain safe and should be kept to MODEL\_CLASSIC only.
!
!  For some cases, the PtrOut may point to a pointer array allocated somewhere else.
!
!  If ZBND is optionally specified, the regridding operation is capped up to the
!  given bound. This is usually used for 2-D data but may also limit useless regridding
!  beyond the pointer size. The array temporaries in the module are usually always
!  allocated to the model top to save memory thrashing, but we also do not want to
!  sacrifice compute.
!
! !REVISION HISTORY:
!  03 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                             :: NZ, L

    IF ( PRESENT( ZBND ) ) THEN
      NZ = MIN(ZBND, State_Grid%NZ+1)                ! May be equal to NZ+1 maximum
    ELSE
      NZ = State_Grid%NZ
    ENDIF

    ! Intermediate grid functionality?
    IF ( .not. Input_Opt%LIMGRID ) RETURN

    ! Initialize shadow regridding handles if they are not ready
    IF ( .not. ASSOCIATED( LonEdgeM ) ) THEN
      CALL Init_IMGrid ( Input_Opt, State_Grid, State_Grid_HCO )
    ENDIF

    ! Reset as instructed
    IF ( PRESENT(ResetRegrName) ) THEN
      IF ( ResetRegrName ) THEN
        LAST_TMP_REGRID_M2H = "_HCO_Dummy_Reset_"
      ENDIF
    ENDIF

    ! Empty the target first
    PtrOut(:,:,1:NZ) = 0.0_hp

    ! Do the regridding layer by layer
    DO L = 1, NZ
      call Map_A2A (                             &
        State_Grid%NX,     State_Grid%NY,        &   ! Input grid dimensions
        LonEdgeM,          LatEdgeM,             &   ! Lons and lat sines at input edges
        PtrIn(:,:,L),                            &   ! Input data, 2-D slice
        State_Grid_HCO%NX, State_Grid_HCO%NY,    &   ! Output grid dimensions
        LonEdgeH,          LatEdgeH,             &
        PtrOut(:,:,L),                           &
        0, 0, missval=-1.e31_hp )                    ! Pole treatment, scalar quantity
    ENDDO

    ! Return happily ever after

  END SUBROUTINE Regrid_MDL2HCO
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_IMGrid
!
! !DESCRIPTION: Subroutine Init\_IMGrid initializes regridding shadow variables
!  used for regridding from GEOS-Chem Classic model grid to the HEMCO grid,
!  the InterMediate Grid ("IMGrid")
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_IMGrid ( Input_Opt, State_Grid, State_Grid_HCO )
!
! !USES:
!
    USE State_Grid_Mod,    ONLY : GrdState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE Regrid_A2A_Mod,    ONLY : Map_A2A
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),   INTENT(IN   )     :: Input_Opt        ! Input opts
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid       ! Grid state
    TYPE(GrdState),   INTENT(IN   )     :: State_Grid_HCO   ! Optional, HEMCO intermediate
!
! !REMARKS:
!  Note that at this point we are committed to IMGrid being turned on and the grids
!  are different. Otherwise the memory is not even allocated.
!
! !REVISION HISTORY:
!  04 Jun 2020 - H.P. Lin  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                             :: RC
    CHARACTER(LEN=255)                  :: ThisLoc
    CHARACTER(LEN=512)                  :: ErrMsg

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_IMGrid (in module GeosCore/hco_utilities_gc_mod.F90)'

    ! Intermediate grid functionality?
    IF ( .not. Input_Opt%LIMGRID ) RETURN

    ! Are we allocated?
    IF ( ASSOCIATED(LonEdgeM) ) RETURN

    ! TMP_MDL, TMP_HCO 3D IJK, TMP_MDLb (alternate model buffer)
    ! LAST_TMP_REGRID_M2H, LAST_TMP_REGRID_H2M, H2Mb
    ! LonEdgeH, LatEdgeH (actually SIN), LonEdgeM, LatEdgeM (actually SIN) -- NX+1, NY+1

    ! Allocate arrays
    ALLOCATE(LonEdgeM(State_Grid%NX+1    ), LatEdgeM(State_Grid%NY+1    ), STAT=RC)
    ALLOCATE(LonEdgeH(State_Grid_HCO%NX+1), LatEdgeH(State_Grid_HCO%NY+1), STAT=RC)

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in allocating model and HEMCO coords!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF

    ! Describe the model grid first
    LonEdgeM(:) = State_Grid%XEDGE(:,1)
    LatEdgeM(:) = State_Grid%YSIN (1,:)

    ! The HEMCO grid
    LonEdgeH(:) = State_Grid_HCO%XEDGE(:,1)
    LatEdgeH(:) = State_Grid_HCO%YSIN (1,:)

    ! Init the temporary targets
    ! These may be conservatively larger than necessary and ZBND is capped
    ! in the regridding routine depending on the given array dims.
    !
    ! Using NZ+1 in case the data is on vertical grid edges.
    ! Probably used for some met field stuff.
    ALLOCATE(TMP_MDL (State_Grid    %NX, State_Grid    %NY, State_Grid%NZ+1), STAT=RC)
    ALLOCATE(TMP_MDLb(State_Grid    %NX, State_Grid    %NY, State_Grid%NZ+1), STAT=RC)
    ALLOCATE(TMP_HCO (State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid%NZ+1), STAT=RC)

    ALLOCATE(TMP_MDL_r4 (State_Grid    %NX, State_Grid    %NY, State_Grid%NZ+1), STAT=RC)
    ALLOCATE(TMP_MDL_r4b(State_Grid    %NX, State_Grid    %NY, State_Grid%NZ+1), STAT=RC)
    ALLOCATE(TMP_HCO_r4 (State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid%NZ+1), STAT=RC)

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in allocating model and HEMCO temporaries!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF

    ! Initialize defaults
    LAST_TMP_REGRID_M2H = ''
    LAST_TMP_REGRID_H2M = ''
    LAST_TMP_REGRID_H2Mb = ''
    LAST_TMP_MDL_ZBND = 1

  END SUBROUTINE Init_IMGrid
!EOC
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_GetOption
!
! !DESCRIPTION: Shim interface for HCO_GetOpt.  Used to return a value
!  from the HEMCO config file into GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GC_GetOption( optName, extNr ) RESULT( optVal )
!
! !USES:
!
    USE HCO_State_GC_Mod, ONLY : HcoState
    USE HCO_ExtList_Mod,  ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: optName
    INTEGER,          OPTIONAL   :: extNr
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255)           :: optVal
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER :: extension

    ! Extension number
    extension = 0
    IF ( PRESENT( extNr ) ) extension = extNr

    ! Return character variable
    optVal = HCO_GetOpt( HcoState%Config%ExtList, optName, extNr=extension )

  END FUNCTION HCO_GC_GetOption
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_HcoStateOK
!
! !DESCRIPTION: Returns TRUE if the HcoState object is allocated, or FALSE
!  otherwise.  This allows us to abstract the HcoState object out of
!  GEOS-Chem source code.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GC_HcoStateOk() RESULT( HcoStateIsOk )
!
! !USES:
!
    USE HCO_State_GC_Mod, ONLY : HcoState
!
! !RETURN VALUE:
!
    LOGICAL :: HcoStateIsOK
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Return status of HcoState
    HcoStateIsOK = ASSOCIATED( HcoState )

  END FUNCTION HCO_GC_HcoStateOk
!EOC

END MODULE HCO_Utilities_GC_Mod
