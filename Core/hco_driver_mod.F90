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
  SUBROUTINE HCO_Run( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_CalcEmis
    USE HCO_tIdx_Mod,     ONLY : tIDx_Update
    USE HCO_ReadList_Mod, ONLY : ReadList_Read 
    USE HCO_ReadList_Mod, ONLY : ReadList_to_EmisList 
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root   ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state object
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller   - Initialization
!  16 Jul 2014 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCO_RUN begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( 'HCO_RUN (hco_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! 1. Update the time slice indeces
    ! This is to make sure that the correct time slices will be used
    ! for all emission fields. See also hco\_tidx\_mod.F90. 
    !=================================================================
    CALL tIDx_Update( am_I_Root, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! 2. Read/update data
    ! Check if there are any data files that need to be read or 
    ! updated, e.g. on the first call of HEMCO or if we enter a new
    ! month, year, etc. 
    !=================================================================

    ! Read/update data, as specified in ReadList.
    CALL ReadList_Read( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Make sure that content of EmisList is up-to-date. This call
    ! primarily ensures that pre-calculation optimizations are
    ! correctly applied (e.g. to unify arrays with same properties). 
    CALL ReadList_to_EmisList ( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    !=================================================================
    ! 3. Calculate the emissions for current time stamp based on the
    ! content of EmisList. Emissions become written into HcoState. 
    !=================================================================
    CALL HCO_CalcEmis( am_I_Root, HcoState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    CALL HCO_LEAVE( RC ) 

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
  SUBROUTINE HCO_Init( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_tIdx_Mod,     ONLY : tIDx_Init
    USE HCO_Clock_Mod,    ONLY : HcoClock_Init
    USE HCO_ReadList_Mod, ONLY : ReadList_Init
    USE HCO_Config_Mod,   ONLY : SetReadList 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! root CPU?
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
    !=================================================================
    ! HCO_INIT begins here!
    !=================================================================

    ! Enter
    CALL HCO_Enter( 'HCO_INIT (hco_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize time slice pointers 
    CALL tIDx_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize HEMCO Clock
    CALL HcoClock_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize the HEMCO ReadList. This has to be done before
    ! the call to SetReadList below. 
    CALL ReadList_Init() 

    ! Set ReadList based upon the content of the configuration file. 
    CALL SetReadList ( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

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
  SUBROUTINE HCO_Final() 
!
! !USES:
!
    USE HCO_EmisList_Mod,  ONLY : EmisList_Cleanup
    USE HCO_Clock_Mod,     ONLY : HcoClock_Cleanup
    USE HCO_tIdx_Mod,      ONLY : tIDx_Cleanup
    USE HCO_ReadList_Mod,  ONLY : ReadList_Cleanup
    USE HCO_Config_Mod,    ONLY : Config_Cleanup
    USE HCO_DataCont_Mod,  ONLY : cIDList_Cleanup
    USE HCO_DataCont_Mod,  ONLY : Reset_nnDataCont
    USE HCO_ExtList_Mod,   ONLY : ExtFinal
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

    CALL cIDList_Cleanup  (         ) 
    CALL HcoClock_Cleanup (         )
    CALL tIDx_Cleanup     (         )
    CALL EmisList_Cleanup ( .FALSE. )
    CALL ReadList_Cleanup ( .FALSE. )
    CALL Config_Cleanup   ( .TRUE.  )
    CALL Reset_nnDataCont

    ! Cleanup the extension list object
    CALL ExtFinal         (         )

    ! Close the logfile and cleanup error object. 
    CALL HCO_Error_Final  (         )

  END SUBROUTINE HCO_Final
!EOC
END MODULE HCO_Driver_Mod
