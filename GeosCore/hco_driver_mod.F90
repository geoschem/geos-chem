!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_driver_mod
!
! !DESCRIPTION: Module HCO\_DRIVER\_MOD contains the driver routines 
! (INIT, RUN, FINAL) for the HEMCO core module. It calls all the 
! subroutines to initialize, execute and finalize the HEMCO core 
! emissions calculations, i.e. all emissions not calculated in a HEMCO
! extension (See module HCOX\_DRIVER\_MOD.F90 for the extensions).\\
! Call this module at the HEMCO - model interface level to execute the
! HEMCO core operations.
!\\
!\\
! !INTERFACE: 
!
      MODULE HCO_DRIVER_MOD 
! 
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_STATE_MOD,       ONLY : HCO_State 

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCO_INIT 
      PUBLIC :: HCO_RUN 
      PUBLIC :: HCO_FINAL 
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL TYPES:
!
      CONTAINS

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_run 
!
! !DESCRIPTION: Subroutine HCO\_RUN is the HEMCO core run routine. It
! calculates the HEMCO emissions as specified in the HEMCO configuration 
! file. All calculation settings, such as the extension number, the 
! lowest and highest species ID, and the lowest and highest emission 
! category, are passed through the HEMCO options object (HcoState%Opt).
! The time stamp is taken from the HEMCO clock object. Subroutine
! HcoClock\_Set should be used to update the HEMCO clock (see module 
! HCO\_CLOCK\_MOD.F90). This should be done at the HEMCO - model 
! interface level. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_RUN ( am_I_Root, HcoState, RC ) 
!
! !USES:
!
      USE HCO_CALC_MOD,        ONLY : HCO_CalcEmis
      USE HCO_TIDX_MOD,        ONLY : tIDx_Update
      USE HCO_READLIST_MOD,    ONLY : ReadList_Read 
      USE HCO_READLIST_MOD,    ONLY : ReadList_to_EmisList 
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root   ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState    ! HEMCO state object
      INTEGER,         INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCO_RUN begins here!
      !=================================================================

      ! Enter
      CALL HCO_ENTER ( 'HCO_RUN (HCO_MOD.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! 1. Update the time slice indeces
      ! This is to make sure that the correct time slices will be used
      ! for all emission fields. See also hco\_tidx\_mod.F90. 
      !=================================================================
      CALL tIDx_Update ( am_I_Root, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! 2. Read/update data
      ! Check if there are any data files that need to be read or 
      ! updated, e.g. on the first call of HEMCO or if we enter a new
      ! month, year, etc. 
      !=================================================================

      ! Read/update data, as specified in ReadList.
      CALL ReadList_Read ( am_I_Root, HcoState, RC )
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
      CALL HCO_CalcEmis ( am_I_Root, HcoState, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Leave w/ success
      CALL HCO_LEAVE ( RC ) 

      END SUBROUTINE HCO_RUN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_INIT
!
! !DESCRIPTION: Subroutine HCO\_INIT initializes the HEMCO core modules.
! This routine assumes that the HEMCO configuration file has been read
! to buffer already (subroutine Config\_ReadFile in HCO\_CONFIG\_MOD.F90)
! and that the HEMCO state object has already been initialized. This has
! to be done at the HEMCO - model interface level.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_INIT( am_I_Root, HcoState, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod

      USE HCO_TIDX_MOD,        ONLY : tIDx_Init
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_Init
      USE HCO_READLIST_MOD,    ONLY : ReadList_Init
      USE HCO_CONFIG_MOD,      ONLY : SetReadList 
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(HCO_State),  POINTER        :: HcoState   ! HcoState object
      INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! HCO_INIT begins here!
      !=================================================================

      ! Enter
      CALL HCO_ENTER ( 'HCO_INIT (HCO_DRIVER_MOD.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize time slice pointers 
      CALL tIDx_Init ( HcoState, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize HEMCO Clock
      CALL HcoClock_Init ( HcoState, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize the HEMCO ReadList. This has to be done before
      ! the call to SetReadList below. 
      CALL ReadList_Init 

      ! Set ReadList based upon the content of the configuration file. 
      CALL SetReadList ( am_I_Root, HcoState, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Leave w/ success
      CALL HCO_LEAVE ( RC ) 

      END SUBROUTINE HCO_INIT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_final 
!
! !DESCRIPTION: Subroutine HCO\_FINAL finalizes HEMCO core. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_FINAL 
!
! !USES:
!
      USE HCO_EMISLIST_Mod,  ONLY : EmisList_Cleanup
      USE HCO_CLOCK_MOD,     ONLY : HcoClock_Cleanup
      USE HCO_TIDX_MOD,      ONLY : tIDx_Cleanup
      USE HCO_READLIST_MOD,  ONLY : ReadList_Cleanup
      USE HCO_CONFIG_MOD,    ONLY : Config_Cleanup
      USE HCO_DATACONT_MOD,  ONLY : cIDList_Cleanup
!
! !ARGUMENTS:
!
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCO_FINAL begins here 
      !=================================================================

      ! Notes: 
      ! - ConfigFile_Cleanup also cleans up the data containers, 
      !   while EmisList_Cleanup and ReadList_Cleanup only removes the 
      !   pointers to them. Hence, we have to call these routines
      !   before ConfigFile_Cleanup! 
      ! - HcoState is cleaned up in the HEMCO-module interface. 
      CALL cIDList_Cleanup
      CALL HcoClock_Cleanup
      CALL tIDx_Cleanup
      CALL EmisList_Cleanup ( .FALSE. )
      CALL ReadList_Cleanup ( .FALSE. )
      CALL Config_Cleanup   ( .TRUE.  )

      ! Close the logfile and cleanup error object. 
      CALL HCO_ERROR_FINAL

      END SUBROUTINE HCO_FINAL 
!EOC
      END MODULE HCO_DRIVER_MOD 
