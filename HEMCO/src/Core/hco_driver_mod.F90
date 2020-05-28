!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_driver_mod.F90
!
! !DESCRIPTION: Module HCO\_Driver\_Mod contains the driver routines
! (INIT, RUN, FINAL) for the HEMCO core module. It calls all the
! subroutines to initialize, execute and finalize the HEMCO core
! emissions calculations, i.e. all emissions not calculated in a HEMCO
! extension (See module hcox\_driver\_mod.F90 for the extensions).
!\\
!\\
! Call this module at the HEMCO - model interface level to execute the
! HEMCO core operations.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Driver_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_State_Mod, ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_Init
  PUBLIC :: HCO_Run
  PUBLIC :: HCO_Final
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller   - Initialization
!  11 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
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
! !IROUTINE: HCO_Run
!
! !DESCRIPTION: Subroutine HCO\_Run is the HEMCO core run routine. It
! calculates the HEMCO emissions as specified in the HEMCO configuration
! file. All calculation settings, such as the extension number, the
! lowest and highest species ID, and the lowest and highest emission
! category, are passed through the HEMCO options object (HcoState%Opt).
! The time stamp is taken from the HEMCO clock object. Subroutine
! HcoClock\_Set should be used to update the HEMCO clock (see module
! hco\_clock\_mod.F90). This should be done at the HEMCO - model
! interface level.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Run( HcoState, Phase, RC, IsEndStep )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_CalcEmis
    USE HCO_ReadList_Mod, ONLY : ReadList_Read
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_Clock_Mod,    ONLY : HcoClock_InitTzPtr
    USE HCOIO_DIAGN_MOD,  ONLY : HcoDiagn_Write
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ) :: Phase       ! Run phase (1 or 2)
    LOGICAL,         INTENT(IN   ), OPTIONAL :: IsEndStep ! Last timestep of simulation?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state object
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller   - Initialization
!  16 Jul 2014 - R. Yantosca - Cosmetic changes
!  23 Dec 2014 - C. Keller   - ReadList_to_EmisList is now obsolete.
!                              Containers are added to EmisList within
!                              routine ReadList_Read.
!  23 Feb 2015 - R. Yantosca - Now call HcoClock_InitTzPtr on the first
!                              emissions timestep to initialize the pointer
!                              to the timezones data (i.e. hours from UTC)
!  01 Apr 2015 - C. Keller   - Added run phases
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IsEmisTime
    LOGICAL            :: notDryRun
    LOGICAL            :: ItIsEndStep

    ! Strings
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! HCO_RUN begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCO_RUN (hco_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Define a local convenience variable to negate HcoState%Options%isDryRun
    notDryRun = ( .not. HcoState%Options%isDryRun )

    ! Define a shadow variable for optional argument IsEndStep
    IF ( PRESENT( IsEndStep ) ) THEN
       ItIsEndStep = IsEndStep
    ELSE
       ItIsEndStep = .FALSE.
    ENDIF

    !--------------------------------------------------------------
    ! 1. Check if it's time for emissions
    !--------------------------------------------------------------
    CALL HcoClock_Get ( HcoState%Clock, IsEmisTime=IsEmisTime, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !--------------------------------------------------------------
    ! 2. Write HEMCO diagnostics. Do this only if the corresponding
    ! option is enabled. Otherwise, let the user decide when to
    ! call HcoDiagn_Write.  Skip if it is a GEOS-Chem "dry-run".
    !--------------------------------------------------------------
    IF ( HcoState%Options%HcoWritesDiagn .and. notDryRun ) THEN
       CALL HcoDiagn_Write( HcoState, .FALSE., RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Check if this is the last timestep of simulation. If so, return
    ! and not read in update data.
    IF ( ItIsEndStep .and. notDryRun ) RETURN

    !--------------------------------------------------------------
    ! 3. Read/update data
    !
    ! Check if there are any data files that need to be read or
    ! updated, e.g. on the first call of HEMCO or if we enter a new
    ! month, year, etc.
    !
    ! NOTE: If this is a GEOS-Chem "dry-run", then HEMCO will
    ! print the files that will be read to either the stdout
    ! (log file) and HEMCO log file, but will not read them.
    !--------------------------------------------------------------

    ! Update data, as specified in ReadList.
    IF ( Phase /= 2 ) THEN
       CALL ReadList_Read( HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *, "Error in ReadList_Read called from hco_run"
          RETURN
       ENDIF

       ! If we are reading timezone data (i.e. offsets from UTC in hours)
       ! from a file, then we need to initialize the TIMEZONES pointer
       ! variable in hco_clock_mod.F90.  This has to be done only on the
       ! very first emissions timestep, after the call to READLIST_READ.
       ! We must leave this call here (instead of in the more customary
       ! initialization routine HCO_INIT) because the HEMCO configuration
       ! file has to be read in its entirety before the timezone data
       ! is loaded into a data container. (bmy, 2/23/15)
       IF ( HcoClock_First(HcoState%Clock,.FALSE.) ) THEN
          CALL HcoClock_InitTzPtr( HcoState, RC )
       ENDIF

    ENDIF


    !-----------------------------------------------------------------
    ! 4. Calculate the emissions for current time stamp based on the
    ! content of EmisList. Emissions become written into HcoState.
    ! Do this only if it's time for emissions and NOT a dry-run.
    !-----------------------------------------------------------------
    IF ( IsEmisTime .AND. Phase == 2 .and. notDryRun ) THEN

       ! Use emission data only
       CALL HCO_CalcEmis( HcoState, .FALSE., RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Use concentration data only
       ! This is currently not being used. Concentrations can be read
       ! through HEMCO but should be assembled manually.
       !CALL HCO_CalcEmis( HcoState, .TRUE., RC )
       !IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCO_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Init
!
! !DESCRIPTION: Subroutine HCO\_INIT initializes the HEMCO core modules.
! This routine assumes that the HEMCO configuration file has already been
! read to buffer (subroutine Config\_ReadFile in HCO\_CONFIG\_MOD.F90)
! and that the HEMCO state object has already been initialized. This has
! to be done at the HEMCO - model interface level.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Init( HcoState, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,    ONLY : HcoDiagn_Init
    USE HCO_tIdx_Mod,     ONLY : tIDx_Init
    USE HCO_Clock_Mod,    ONLY : HcoClock_Init
    USE HCO_Config_Mod,   ONLY : SetReadList
    USE HCO_Scale_Mod,    ONLY : HCO_ScaleInit
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! HcoState object
    INTEGER,          INTENT(INOUT) :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller   - Initialization
!  11 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  16 Jul 2014 - R. Yantosca - Remove reference to gigc_errcode_mdo.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! !LOCAL VARIABLES:
!
    !=================================================================
    ! HCO_INIT begins here!
    !=================================================================

    ! Enter
    CALL HCO_Enter( HcoState%Config%Err, 'HCO_INIT (hco_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HCO_Enter called from HCO_Init"
       RETURN
    ENDIF

    ! Initialize time slice pointers
    CALL tIDx_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in tIDx_Init called from HCO_Init"
       RETURN
    ENDIF

    ! Initialize HEMCO Clock
    HcoState%Clock => NULL()
    CALL HcoClock_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HcoClock_Init called from HCO_Init"
       RETURN
    ENDIF

    ! Initialize the HEMCO diagnostics
    CALL HcoDiagn_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in HcoDiagn_Init called from HCO_Init"
       RETURN
    ENDIF

    ! Set ReadList based upon the content of the configuration file.
    CALL SetReadList( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in SetReadList called from HCO_Init"
       RETURN
    ENDIF

    ! Define universal scale factor for each HEMCO species
    CALL Hco_ScaleInit( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *, "Error in Hco_ScaleInit called from HCO_Init"
       RETURN
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE HCO_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Final
!
! !DESCRIPTION: Subroutine HCO\_Final finalizes HEMCO core.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Final( HcoState, ERROR, RC )
!
! !USES:
!
    USE HCO_Clock_Mod,     ONLY : HcoClock_Cleanup
    USE HCO_tIdx_Mod,      ONLY : tIDx_Cleanup
    USE HCO_DataCont_Mod,  ONLY : cIDList_Cleanup
    USE HCO_ReadList_Mod,  ONLY : ReadList_Cleanup
    USE HCO_DataCont_Mod,  ONLY : ListCont_Cleanup
    USE HCO_ExtList_Mod,   ONLY : ExtFinal
    USE HCOIO_DIAGN_MOD,   ONLY : HcoDiagn_Write
    USE HCO_Scale_Mod,     ONLY : HCO_ScaleFinal
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: ERROR      ! Cleanup because of crash?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! HcoState object
    INTEGER,          INTENT(INOUT) :: RC         ! Failure or success
!
! !REMARKS:
!  (1) ConfigFile_Cleanup also cleans up the data containers, while routine
!       EmisList_Cleanup and ReadList_Cleanup only removes the pointers to
!       them. Hence, we have to call these routines before ConfigFile_Cleanup!
!  (2) HcoState is cleaned up in the HEMCO-module interface.
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller   - Initialization
!  11 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  16 Jul 2014 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCO_FINAL begins here
    !=================================================================

    ! Write diagnostics if needed
    IF ( HcoState%Options%HcoWritesDiagn .AND. .NOT. ERROR ) THEN
       CALL  HcoDiagn_Write( HcoState, .FALSE., RC )
       CALL  HcoDiagn_Write( HcoState, .TRUE.,  RC )
    ENDIF

    CALL cIDList_Cleanup  ( HcoState                           )
    CALL HcoClock_Cleanup ( HcoState%Clock                     )
    CALL tIDx_Cleanup     ( HcoState%AlltIDx                   )
    CALL ReadList_Cleanup ( HcoState%ReadLists,        .FALSE. )
    CALL ListCont_Cleanup ( HcoState%EmisList,         .FALSE. )
    CALL ListCont_Cleanup ( HcoState%Config%ConfigList, .TRUE. )
    HcoState%nnEmisCont = 0
    HcoState%SetReadListCalled = .FALSE.

    ! Cleanup scaling factors
    CALL HCO_ScaleFinal()

    ! Cleanup the extension list object
    CALL ExtFinal         ( HcoState%Config%ExtList )

    ! Close the logfile and cleanup error object.
    CALL HCO_Error_Final  ( HcoState%Config%Err )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Final
!EOC
END MODULE HCO_Driver_Mod
