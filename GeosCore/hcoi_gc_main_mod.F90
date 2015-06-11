!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_main_mod.F90 
!
! !DESCRIPTION: Module hcoi\_gc\_main\_mod.F90 is the HEMCO-to-GEOS-Chem
! interface module, providing the link between GEOS-Chem and HEMCO. It 
! contains wrapper routines to initialize, execute and finalize HEMCO from
! within GEOS-Chem. These routines are called from emissions\_mod.F90.
!\\
!\\
! Notes:
! \begin{itemize}
! \item HEMCO is used to calculate all emission fields. The emission tendencies
!  are passed to GEOS-Chem in module mixing\_mod.F90. 
! \item Most meteorological fields needed by the HEMCO extensions are provided
!  through the GEOS-Chem meteorological state object Met\_State. Few fields 
!  such as the pressure edges or J-values are defined and updated explicitly 
!  within this module.
! \end{itemize}
! !INTERFACE:
!
MODULE HCOI_GC_Main_Mod
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
  PUBLIC  :: HCOI_GC_Init
  PUBLIC  :: HCOI_GC_Run
  PUBLIC  :: HCOI_GC_Final
  PUBLIC  :: HCOI_GC_WriteDiagn
  PUBLIC  :: GetHcoState
  PUBLIC  :: GetHcoVal
  PUBLIC  :: GetHcoID
  PUBLIC  :: GetHcoDiagn
  PUBLIC  :: SetHcoTime
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Calc_SumCosZa
  PRIVATE :: ExtState_InitTargets
  PRIVATE :: ExtState_SetFields
  PRIVATE :: ExtState_UpdateFields
  PRIVATE :: Get_SzaFact
  PRIVATE :: GridEdge_Set
  PRIVATE :: Set_Grid
  PRIVATE :: CheckSettings
  PRIVATE :: SetHcoSpecies 
!
! !REMARKS:
!  This module is ignored if you are using HEMCO in an ESMF environment.
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version. 
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  30 Jul 2014 - C. Keller   - Added GetHcoState 
!  20 Aug 2014 - M. Sulprizio- Modify for POPs simulation
!  21 Aug 2014 - R. Yantosca - Added routine EmissRnPbBe; cosmetic changes
!  06 Oct 2014 - C. Keller   - Removed PCENTER. Now calculate from pressure edges
!  21 Oct 2014 - C. Keller   - Removed obsolete routines MAP_HCO2GC and 
!                              Regrid_Emis2Sim. Added wrapper routine GetHcoID
!  18 Feb 2015 - C. Keller   - Added routine CheckSettings.
!  04 Mar 2015 - C. Keller   - Now register all GEOS-Chem species as HEMCO
!                              species. 
!  11 Mar 2015 - R. Yantosca - Now move computation of SUMCOSZA here from 
!                              the obsolete global_oh_mod.F.  Add routines
!                              GET_SZAFACT and CALC_SUMCOSA.
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
  TYPE(HCO_State),      POINTER :: HcoState               => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State),      POINTER :: ExtState               => NULL()

  !--------------------------
  ! %%% Arrays %%%
  !--------------------------

  ! Internal met fields (will be used by some extensions)
  INTEGER,               TARGET :: HCO_PBL_MAX            ! level
  REAL(hp), ALLOCATABLE, TARGET :: HCO_FRAC_OF_PBL(:,:,:) ! unitless
  REAL(hp), ALLOCATABLE, TARGET :: HCO_SZAFACT(:,:)       ! -

  ! Arrays to store J-values (used by Paranox extension)
  REAL(hp), ALLOCATABLE, TARGET :: JNO2(:,:)
  REAL(hp), ALLOCATABLE, TARGET :: JO1D(:,:)

  ! Sigma coordinate (temporary)
  REAL(hp), ALLOCATABLE, TARGET :: ZSIGMA(:,:,:)

  ! Sum of cosine of the solar zenith angle. Used to impose a
  ! diurnal variability on OH concentrations
  REAL(fp), ALLOCATABLE         :: SUMCOSZA(:,:)
!
! !DEFINED PARAMETERS:
!
  ! Temporary toggle for diagnostics
  LOGICAL,  PARAMETER           :: DoDiagn = .TRUE.

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_INIT initializes the HEMCO derived
! types and arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE TIME_MOD,           ONLY : GET_TS_EMIS, GET_TS_DYN
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE ERROR_MOD,          ONLY : ERROR_STOP
#if defined( TOMAS ) 
    USE TOMAS_MOD,          ONLY : IBINS
    USE TOMAS_MOD,          ONLY : Xk
#endif

    ! HEMCO routines 
    USE HCO_Config_Mod,     ONLY : Config_ReadFile
    USE HCO_State_Mod,      ONLY : HcoState_Init
    USE HCO_Driver_Mod,     ONLY : HCO_Init
    USE HCOI_GC_Diagn_Mod,  ONLY : HCOI_GC_Diagn_Init
    USE HCOX_Driver_Mod,    ONLY : HCOX_Init
    USE HCOX_State_Mod,     ONLY : ExtStateInit
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!  07 Jul 2014 - C. Keller    - Now match species and set species properties
!                               via module variables.
!  30 Sep 2014 - R. Yantosca  - Now pass fields for aerosol and microphysics
!                               options to extensions via HcoState
!  13 Feb 2015 - C. Keller    - Now read configuration file in two steps.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                         :: LSTRAT,  FOUND
    INTEGER                         :: nHcoSpc, HMRC
    CHARACTER(LEN=255)              :: OptName, LOC

    !=================================================================
    ! HCOI_GC_INIT begins here!
    !=================================================================

    ! Error handling 
    LOC = 'HCOI_GC_Init (hcoi_gc_main_mod.F90)'

    ! Set return code flag to HCO success. This value should be
    ! preserved throughout all HCO calls, otherwise an error
    ! will be returned!
    HMRC = HCO_SUCCESS

    !=================================================================
    ! Read HEMCO configuration file and save into buffer. This also
    ! sets the HEMCO error properties (verbose mode? log file name, 
    ! etc.) based upon the specifications in the configuration file.
    ! The log file is now read in two phases: phase 1 reads only the
    ! settings and extensions; phase 2 reads all data fields. This 
    ! way, settings and extension options can be updated before 
    ! reading all the associated fields. For instance, if the LEMIS
    ! toggle is set to false (=no emissions), all extensions can be
    ! deactivated. Similarly, certain brackets can be set explicitly
    ! to make sure that these data is only read by HEMCO if the 
    ! corresponding GEOS-Chem switches are turned on.
    ! (ckeller, 2/13/15).
    !=================================================================

    ! Phase 1: read settings and switches
    CALL Config_ReadFile( am_I_Root, Input_Opt%HcoConfigFile, 1, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Config_ReadFile', LOC )

    ! Check settings
    CALL CheckSettings( am_I_Root, Input_Opt, State_Met, State_Chm, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'CheckSettings', LOC )

    ! Phase 2: read fields
    CALL Config_ReadFile( am_I_Root, Input_Opt%HcoConfigFile, 2, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Config_ReadFile', LOC )

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LOGFILE_OPEN( RC=HMRC ) 
       IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Open Logfile', LOC )
    ENDIF

    !=================================================================
    ! Initialize HEMCO state object and populate it 
    !=================================================================

    !-----------------------------------------------------------------
    ! Extract species to use in HEMCO. nHcoSpc denotes the number of
    ! species that shall be used in HEMCO. The species properties are
    ! defined in the Register_Species call below.
    ! Typically, nHcoSpc is just the number of species defined in both 
    ! the HEMCO configuration file and GEOS-Chem. However, additional
    ! species can be defined, e.g. those not transported in GEOS-Chem
    ! (e.g. SESQ) or tagged species (e.g. specialty simulations).
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, HcoState, & 
                         nHcoSpc,   1,         HMRC       )
!    CALL Get_nHcoSpc( am_I_Root, Input_Opt, nHcoSpc, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SetHcoSpecies-1', LOC )

    !-----------------------------------------------------------------
    ! Now that number of HEMCO species are known, initialize HEMCO
    ! state object.
    CALL HcoState_Init( am_I_Root, HcoState, nHcoSpc, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HcoState_Init', LOC )

    !-----------------------------------------------------------------
    ! Register species. This will define all species properties
    ! (names, molecular weights, etc.) of the HEMCO species.
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, HcoState, & 
                         nHcoSpc,   2,         HMRC       )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SetHcoSpecies-2', LOC )

    !-----------------------------------------------------------------
    ! Set grid. 
    CALL Set_Grid( am_I_Root, State_Met, HcoState, RC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Set_Grid', LOC )

    !=================================================================
    ! Set misc. parameter
    !=================================================================

    ! Emission, chemistry and dynamics timestep in seconds
    HcoState%TS_EMIS = GET_TS_EMIS() * 60.0
    HcoState%TS_CHEM = GET_TS_CHEM() * 60.0
    HcoState%TS_DYN  = GET_TS_DYN()  * 60.0

    ! Is this an ESMF simulation or not?
    ! The ESMF flag must be set before calling HCO_Init because the
    ! source file name is set differently in an ESMF environment
    ! compared to a stand-alone version: in ESMF, the source file name
    ! is set to the container name since this is the identifying name
    ! used by ExtData.
#if defined(ESMF_)
    HcoState%isESMF = .TRUE.
#else 
    HcoState%isESMF = .FALSE.
#endif

    ! HEMCO configuration file
    HcoState%ConfigFile = Input_Opt%HcoConfigFile

    ! Set deposition length scale. This determines if dry deposition
    ! frequencies are calculated over the entire PBL or the first
    ! model layer only.
    HcoState%Options%PBL_DRYDEP = Input_Opt%PBL_DRYDEP

    !=================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and 
    ! the HEMCO configuration file is removed from buffer in this
    ! step. This also initializes the HEMCO clock as well as the
    ! HEMCO emissions diagnostics collection.
    !=================================================================
    CALL HCO_Init( am_I_Root, HcoState, HMRC )
    IF( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC )

    ! Save # of defined dust species in HcoState
    HcoState%nDust                     =  Input_Opt%N_DUST_BINS

#if defined( TOMAS )

    ! Save # of TOMAS size bins in HcoState
    HcoState%MicroPhys%nBins           =  IBINS

    ! Point to TOMAS bin boundaries array (Xk) in HcoState
    HcoState%MicroPhys%BinBound        => Xk

    ! Save # of TOMAS active mode bins in HcoState
# if defined( TOMAS40 )
    HcoState%MicroPhys%nActiveModeBins =  10
# elif defined( TOMAS15 )
    HcoState%MicroPhys%nActiveModeBins =  3
# else 
    HcoState%MicroPhys%nActiveModeBins =  0
# endif
#endif

    !=================================================================
    ! Initialize all HEMCO extensions. This also selects the required 
    ! met fields used by each extension.
    !=================================================================
    CALL HCOX_Init( am_I_Root, HcoState, ExtState, HMRC )
    IF( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC )

    !-----------------------------------------------------------------
    ! Update logical switches in Input_Opt 
    !-----------------------------------------------------------------

    ! Soil NOx
    Input_Opt%LSOILNOX      = ExtState%SoilNOx

    ! Ginoux dust emissions
    IF ( ExtState%DustGinoux ) THEN
       Input_Opt%LDUST      = .TRUE.
       Input_Opt%LDEAD      = .FALSE.
    ENDIF

    ! DEAD dust emissions
    IF ( ExtState%DustDead ) THEN
       Input_Opt%LDUST      = .TRUE.
       Input_Opt%LDEAD      = .TRUE.
    ENDIF

    !-----------------------------------------------------------------
    ! Set constants for POPs simulation
    !-----------------------------------------------------------------
    IF ( ExtState%GC_POPs ) THEN
       ExtState%POP_DEL_H   = Input_Opt%POP_DEL_H
       ExtState%POP_KOA     = Input_Opt%POP_KOA
       ExtState%POP_KBC     = Input_Opt%POP_KBC
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize ExtState target arrays. 
    ! Extensions typically depend on environmental dependent met. 
    ! variables such as wind speed, surface temp., etc. Pointers 
    ! to these (2D or 3D) fields are defined in the extension object. 
    ! Here, we need to make sure that these pointers are properly 
    ! connected.
    !-----------------------------------------------------------------
    CALL ExtState_InitTargets( am_I_Root, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Define diagnostics
    !-----------------------------------------------------------------
    IF ( DoDiagn ) THEN 

       ! Set up traditional GEOS-Chem NDxx diagnostics for emissions
       CALL HCOI_GC_DIAGN_INIT                                &
            ( am_I_Root, Input_Opt, HcoState, ExtState, HMRC )

       ! Exit if any of the diagnostics could not be initialized
       IF ( HMRC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'HCOI_GC_DIAGN_INIT', LOC )
       ENDIF
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Leave w/ success
    RC = GIGC_SUCCESS

    END SUBROUTINE HCOI_GC_INIT 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Run
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Run runs HEMCO from GEOS-Chem. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Run( am_I_Root, Input_Opt, State_Met, State_Chm, & 
                          EmisTime,  Phase,     RC                     )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE GIGC_State_Chm_Mod,    ONLY : ChmState

    ! HEMCO routines 
    USE HCO_CLOCK_MOD,         ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,         ONLY : HcoClock_EmissionsDone
    USE HCO_DIAGN_MOD,         ONLY : HcoDiagn_AutoUpdate
    USE HCO_FLUXARR_MOD,       ONLY : HCO_FluxarrReset 
    USE HCO_DRIVER_MOD,        ONLY : HCO_RUN
    USE HCOX_DRIVER_MOD,       ONLY : HCOX_RUN

    ! For soilnox
    USE GET_NDEP_MOD,          ONLY : RESET_DEP_N
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    LOGICAL,          INTENT(IN   )  :: EmisTime   ! Is this an emission time step? 
    INTEGER,          INTENT(IN   )  :: Phase      ! Run phase: 1, 2, -1 (all) 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input options
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Meteo state 
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Modifi
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  22 Aug 2014 - R. Yantosca - Now pass State_Met to MAP_HCO2GC
!  02 Oct 2014 - C. Keller   - PEDGE is now in HcoState%Grid
!  13 Jan 2015 - C. Keller   - Now check if it's time for emissions. Added
!                              call to HcoClock_EmissionsDone.
!  06 Mar 2015 - R. Yantosca - Now create splash page for HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE                  :: FIRST = .TRUE.
    INTEGER                        :: HMRC 
    LOGICAL                        :: IsEmisTime
    CHARACTER(LEN=255), PARAMETER  :: LOC='HCOI_GC_RUN (hcoi_gc_main_mod.F90)'

    !=======================================================================
    ! HCOI_GC_RUN begins here!
    !=======================================================================

    ! Set return code flag to HCO success. This value should be
    ! preserved throughout all HCO calls, otherwise an error
    ! will be returned!
    HMRC = HCO_SUCCESS

    ! Create a splash page
    IF ( am_I_Root .and. FIRST ) THEN 
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
       WRITE( 6, 100   ) 'HEMCO: Harvard-NASA Emissions Component'
       WRITE( 6, 101   ) 'You are using HEMCO version ', ADJUSTL(HCO_VERSION)
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
 100   FORMAT( '%%%%%', 15x, a,      15x, '%%%%%' )
 101   FORMAT( '%%%%%', 15x, a, a12, 14x  '%%%%%' )
       FIRST = .FALSE.
    ENDIF

    !=======================================================================
    ! Make sure HEMCO time is in sync with simulation time
    ! This is now done in main.F 
    !=======================================================================
    CALL SetHcoTime ( am_I_Root, EmisTime, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SetHcoTime', LOC )

    !=======================================================================
    ! See if it's time for emissions. Don't just use the EmisTime flag in
    ! case that we call this routine multiple times. IsEmisTime will only
    ! be true if this is an emission time step AND emissions have not yet
    ! been calculated for that time step.
    !=======================================================================
    CALL HcoClock_Get( IsEmisTime=IsEmisTime, RC=HMRC )

    ! ======================================================================
    ! Reset all emission and deposition values. Do this only if it is time
    ! for emissions, i.e. if those values will be refilled.
    ! ======================================================================
    IF ( IsEmisTime .AND. Phase /= 1 ) THEN
       CALL HCO_FluxarrReset ( HcoState, HMRC )
       IF ( HMRC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP('ResetArrays', LOC )
          RETURN 
       ENDIF
    ENDIF

    !=======================================================================
    ! Define pressure edges [Pa] on HEMCO grid.
    !=======================================================================
    CALL GridEdge_Set ( State_Met, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('GridEdge_Update', LOC )
       RETURN 
    ENDIF
 
    !=======================================================================
    ! Set HCO options 
    !=======================================================================

    ! Range of tracers and emission categories.
    ! Set Extension number ExtNr to 0, indicating that the core
    ! module shall be executed. 
    HcoState%Options%SpcMin = 1 
    HcoState%Options%SpcMax = -1
    HcoState%Options%CatMin = 1 
    HcoState%Options%CatMax = -1 
    HcoState%Options%ExtNr  = 0

    ! Use temporary array?
    HcoState%Options%FillBuffer = .FALSE. 

    !=======================================================================
    ! Run HCO core module
    ! Pass phase as argument. Phase 1 will update the emissions list,
    ! phase 2 will calculate the emissions. Emissions will be written into 
    ! the corresponding flux arrays in HcoState. 
    !=======================================================================
    CALL HCO_RUN ( am_I_Root, HcoState, Phase, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('HCO_RUN', LOC )
       RETURN 
    ENDIF

    !=======================================================================
    ! Do the following only if it's time to calculate emissions 
    !=======================================================================
    IF ( Phase /= 1 .AND. IsEmisTime ) THEN 

       !-----------------------------------------------------------------
       ! Set / update ExtState fields.
       ! Extensions typically depend on environmental dependent met. 
       ! variables such as wind speed, surface temp., etc. Pointers 
       ! to these (2D or 3D) fields are defined in the extension object. 
       ! Here, we need to make sure that these pointers are properly 
       ! connected.
       !-----------------------------------------------------------------
       CALL ExtState_SetFields( am_I_Root, State_Met, State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP('ExtState_SetFields', LOC )
          RETURN 
       ENDIF 
   
       CALL ExtState_UpdateFields( am_I_Root, State_Met, State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) RETURN
   
       !=======================================================================
       ! Run HCO extensions. Emissions will be added to corresponding
       ! flux arrays in HcoState.
       !=======================================================================
       CALL HCOX_RUN ( am_I_Root, HcoState, ExtState, HMRC )
       IF ( HMRC/= HCO_SUCCESS ) THEN
          CALL ERROR_STOP('HCOX_RUN', LOC )
          RETURN
       ENDIF
   
       !=======================================================================
       ! Update all 'AutoFill' diagnostics. This makes sure that all 
       ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
       ! AutoFill flag is specified when creating a diagnostics container
       ! (Diagn_Create).
       !=======================================================================
       IF ( DoDiagn ) THEN
          CALL HcoDiagn_AutoUpdate ( am_I_Root, HcoState, HMRC )
          IF( HMRC /= HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_UPDATE', LOC )
       ENDIF
   
       !=======================================================================
       ! Reset the accumulated nitrogen dry and wet deposition to zero. Will
       ! be re-filled in drydep and wetdep.
       !=======================================================================
       CALL RESET_DEP_N()
   
       !=======================================================================
       ! Emissions are now done for this time step
       !=======================================================================
       CALL HcoClock_EmissionsDone( am_I_Root, RC )

    ENDIF  
 
    ! We are now back in GEOS-Chem environment, hence set 
    ! return flag accordingly! 
    RC = GIGC_SUCCESS

  END SUBROUTINE HCOI_GC_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Final
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Final cleans up HEMCO. This routine
! should be called before the finalize routines of State\_Chm in order to 
! make sure that the emissions flux pointers are properly removed!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Final( am_I_Root, ERROR )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE Error_Mod,           ONLY : Error_Stop
    USE CMN_SIZE_Mod,        ONLY : IIPAR, JJPAR, LLPAR
    USE HCO_Driver_Mod,      ONLY : HCO_Final
    USE HCO_Diagn_Mod,       ONLY : DiagnCollection_Cleanup
    USE HCO_State_Mod,       ONLY : HcoState_Final
    USE HCOX_Driver_Mod,     ONLY : HCOX_Final
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)              :: am_I_Root
    LOGICAL, INTENT(IN)              :: ERROR
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  19 Feb 2015 - R. Yantosca - Change restart file name back to HEMCO_restart
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: HMRC
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! HCOI_GC_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_GC_Final (hcoi_gc_main_mod.F90)'

    ! Cleanup HCO core. 
    CALL HCO_FINAL( am_I_Root, HcoState, ERROR, HMRC )

    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    CALL HCOX_FINAL ( am_I_Root, HcoState, ExtState, HMRC ) 

    ! Cleanup HcoState object 
    CALL HcoState_Final ( HcoState ) 

    ! Cleanup all diagnostics
    CALL DiagnCollection_Cleanup

    ! Module variables
    IF ( ALLOCATED ( ZSIGMA          ) ) DEALLOCATE( ZSIGMA          )
    IF ( ALLOCATED ( HCO_FRAC_OF_PBL ) ) DEALLOCATE( HCO_FRAC_OF_PBL )
    IF ( ALLOCATED ( HCO_SZAFACT     ) ) DEALLOCATE( HCO_SZAFACT     )
    IF ( ALLOCATED ( JNO2            ) ) DEALLOCATE( JNO2            )
    IF ( ALLOCATED ( JO1D            ) ) DEALLOCATE( JO1D            )
    IF ( ALLOCATED ( SUMCOSZA        ) ) DEALLOCATE( SUMCOSZA        ) 

  END SUBROUTINE HCOI_GC_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_WriteDiagn
!
! !DESCRIPTION: Subroutine HCOI\_GC\_WriteDiagn is the wrapper routine to
! write the HEMCO diagnostics. This will only write the diagnostics of 
! diagnostics collection 1.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, Restart, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE Error_Mod,           ONLY : Error_Stop
    USE GIGC_Input_Opt_Mod,  ONLY : OptInput
    USE HCOIO_Diagn_Mod,     ONLY : HcoDiagn_Write 
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   )    :: am_I_Root    ! Root CPU?
    TYPE(OptInput), INTENT(IN   )    :: Input_Opt    ! Input options
    LOGICAL,        INTENT(IN   )    :: Restart      ! write restart (enforced)?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT)    :: RC         ! Return code
!
! !REVISION HISTORY: 
!  01 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: HMRC
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! HCOI_GC_WriteDiagn begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_GC_WriteDiagn (hcoi_gc_main_mod.F90)'

    ! Make sure HEMCO time is in sync 
    CALL SetHcoTime ( am_I_Root, .FALSE., HMRC )
    IF( HMRC /= HCO_SUCCESS) CALL ERROR_STOP ( 'SetHcoTime', LOC )

    ! Write diagnostics
    CALL HcoDiagn_Write( am_I_Root, HcoState, RESTART, HMRC )
    IF(HMRC/=HCO_SUCCESS) THEN
       WRITE(MSG,*) 'Error writing HEMCO diagnostics' 
       CALL ERROR_STOP ( MSG, LOC )
    ENDIF

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE HCOI_GC_WriteDiagn
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
  SUBROUTINE SetHcoTime( am_I_Root, TimeForEmis, RC ) 
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_Set
    USE TIME_MOD,      ONLY : GET_YEAR, GET_MONTH,  GET_DAY
    USE TIME_MOD,      ONLY : GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR, GET_DAY_OF_WEEK
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root
    LOGICAL,         INTENT(IN   ) :: TimeForEmis 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 Jan 2013 - C. Keller - Now call MAP_A2A instead of DO_REGRID_A2A
!  12 Jan 2015 - C. Keller - Added argument TimeForEmis 
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

    CALL HcoClock_Set ( am_I_Root,  HcoState, cYr, cMt, cDy, cHr, &
                        cMin, cSec, cDoy, IsEmisTime=TimeForEmis, RC=RC )

  END SUBROUTINE SetHcoTime
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_InitTargets
!
! !DESCRIPTION: SUBROUTINE ExtState\_InitTargets allocates some local arrays
! that act as targets for the ExtState object.
!\\
! Note that for now, this explicitly assumes that the HEMCO emissions grid is 
! the same as the GEOS-Chem simulation grid. To support other grids, the met 
! data has to be regridded explicitly at every time step!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_InitTargets( am_I_Root, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller    - Initial Version
!  20 Aug 2014 - M. Sulprizio - Add PBL_MAX and FRAC_OF_PBL for POPs simulation
!  02 Oct 2014 - C. Keller    - PEDGE is now in HcoState%Grid
!  16 Oct 2014 - C. Keller    - Removed SUNCOSmid5. This is now calculated
!                               internally in Paranox.
!  12 Mar 2015 - R. Yantosca  - Allocate SUMCOSZA array for SZAFACT
!  12 Mar 2015 - R. Yantosca  - Use 0.0e0_hp when zeroing REAL(hp) variables
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: AS
    CHARACTER(LEN=255) :: LOC = 'ExtState_InitTargets (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! ExtState_InitTargets begins here
    !=================================================================

    ! Init
    RC = GIGC_SUCCESS

    ! ----------------------------------------------------------------
    ! HCO_SZAFACT is not defined in Met_State.  Hence need to 
    ! define here so that we can point to them.
    !
    ! Now include HCO_FRAC_OF_PBL and HCO_PBL_MAX for POPs specialty
    ! simulation (mps, 8/20/14)
    ! ----------------------------------------------------------------
    IF ( ExtState%SZAFACT%DoUse ) THEN 

       ALLOCATE( SUMCOSZA( IIPAR, JJPAR ), STAT=AS )
       IF ( AS/=0 ) CALL ERROR_STOP ( 'SUMCOSZA', LOC )
       SUMCOSZA = 0.0e0_fp

       ALLOCATE( HCO_SZAFACT( IIPAR, JJPAR      ),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_SZAFACT', LOC )
       HCO_SZAFACT = 0e0_hp
    ENDIF

    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN 
       ALLOCATE(HCO_FRAC_OF_PBL(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_FRAC_OF_PBL', LOC )
       HCO_FRAC_OF_PBL = 0.0_hp
    ENDIF

    ! Initialize max. PBL
    HCO_PBL_MAX = 0

    ! ----------------------------------------------------------------
    ! The J-Values for NO2 and O3 are not defined in Met_State. We
    ! need to compute them separately. 
    ! ----------------------------------------------------------------
    IF ( ExtState%JNO2%DoUse ) THEN 
       ALLOCATE( JNO2(IIPAR,JJPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JNO2', LOC )
       JNO2 = 0.0e0_hp
    ENDIF

    IF ( ExtState%JO1D%DoUse ) THEN 
       ALLOCATE( JO1D(IIPAR,JJPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JO1D', LOC )
       JO1D = 0.0e0_hp
    ENDIF

    ! ----------------------------------------------------------------
    ! Arrays to be copied physically because HEMCO units are not the
    ! same as in GEOS-Chem
    ! ----------------------------------------------------------------

    ! TROPP: GEOS-Chem TROPP is in hPa, while HEMCO uses Pa.
    IF ( ExtState%TROPP%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%TROPP%Arr, IIPAR, JJPAR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%TROPP', LOC )
       ENDIF
    ENDIF

    ! SPHU: GEOS-Chem SPHU is in g/kg, while HEMCO uses kg/kg.
    ! NOTE: HEMCO only uses SPHU surface values.
    IF ( ExtState%SPHU%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SPHU%Arr, IIPAR, JJPAR, 1, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%SPHU', LOC )
       ENDIF
    ENDIF

    ! SUNCOS: HEMCO now calculates SUNCOS values based on its own
    ! subroutine 
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, IIPAR, JJPAR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%SUNCOS', LOC )
       ENDIF
    ENDIF

    ! Leave with success
    RC = GIGC_SUCCESS

  END SUBROUTINE ExtState_InitTargets
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_SetFields
!
! !DESCRIPTION: SUBROUTINE ExtState\_SetFields connects the ExtState fields 
! of the HEMCO ExtState object to its target data. This can be a field in
! State\_Met, State\_Chm, or any other 2D/3D field defined within GEOS-Chem or
! even explicitly calculated in this module. All these fields are expected 
! to be of the same type as the corresponding ExtState object, and a pointer 
! link is established between the two fields on the first call.
!\\
!\\
! The ExtState object fields can also be linked to data fields read through
! the HEMCO configuration file. In this case, the data fields will be copied
! from the HEMCO data list into the ExtState object every time this routine
! is called. The field name of the HEMCO field must match the field name 
! passed to ExtState\_Set. 
!\\
!\\
! Fields from the HEMCO data list are given priority over the target fields from
! Met\_State, Chm\_State, etc. For example, if the HEMCO data list contains 
! a field named 'U10M', this field will be used in ExtState%U10M in lieu of
! State\_Met%U10M. 
!\\  
!\\
! Note that for now, this explicitly assumes that the HEMCO emissions grid is 
! the same as the GEOS-Chem simulation grid. To support other grids, the met 
! data has to be regridded explicitly at every time step!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_SetFields( am_I_Root, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE HCOX_STATE_MOD,        ONLY : ExtDat_Set

    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE GIGC_State_Chm_Mod,    ONLY : ChmState

    ! For ParaNOx
    USE TRACERID_MOD,          ONLY : IDTO3, IDTNO, IDTNO2, IDTHNO3

    ! For SoilNox
    USE Drydep_Mod,            ONLY : DRYCOEFF
    USE Get_Ndep_Mod,          ONLY : DRY_TOTN
    USE Get_Ndep_Mod,          ONLY : WET_TOTN

#if !defined(ESMF_)
    USE MODIS_LAI_MOD,         ONLY : GC_LAI
#endif

#if defined(ESMF_)
    USE HCOI_ESMF_MOD,      ONLY : HCO_SetExtState_ESMF
#endif

!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    INTEGER,          INTENT(INOUT)  :: RC         ! Return code
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller    - Initial Version
!  20 Aug 2014 - M. Sulprizio - Add PBL_MAX and FRAC_OF_PBL for POPs simulation
!  02 Oct 2014 - C. Keller    - PEDGE is now in HcoState%Grid
!  16 Oct 2014 - C. Keller    - Removed SUNCOSmid5. This is now calculated
!                               internally in Paranox.
!  12 Mar 2015 - R. Yantosca  - Allocate SUMCOSZA array for SZAFACT
!  12 Mar 2015 - R. Yantosca  - Use 0.0e0_hp when zeroing REAL(hp) variables
!  03 Apr 2015 - C. Keller    - Now call down to ExtDat_Set for all fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRST = .TRUE.
    INTEGER            :: HCRC
    CHARACTER(LEN=255) :: LOC = 'ExtState_SetFields (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! ExtState_SetFields begins here
    !=================================================================

    ! Init
    RC = GIGC_FAILURE

    ! ----------------------------------------------------------------
    ! Pointers to local module arrays 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SZAFACT, & 
          'SZAFACT', HCRC,      FIRST,    HCO_SZAFACT        )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%JNO2, &
            'JNO2',  HCRC,      FIRST,    JNO2            )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%JO1D, &
            'JO1D',  HCRC,      FIRST,    JO1D            )  
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRAC_OF_PBL, &
     'FRAC_OF_PBL',  HCRC,      FIRST,    HCO_FRAC_OF_PBL )  
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! 2D fields 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%U10M, &
            'U10M',  HCRC,      FIRST,    State_Met%U10M )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%V10M, &
            'V10M',  HCRC,      FIRST,    State_Met%V10M )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%ALBD, &
            'ALBD',  HCRC,      FIRST,    State_Met%ALBD )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%WLI, &
             'WLI',  HCRC,      FIRST,    State_Met%LWI  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%T2M, &
             'T2M',  HCRC,      FIRST,    State_Met%TS   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%TSKIN, &
           'TSKIN',  HCRC,      FIRST,    State_Met%TSKIN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%GWETROOT, &
         'GWETROOT', HCRC,      FIRST,    State_Met%GWETROOT  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%GWETTOP, &
          'GWETTOP', HCRC,      FIRST,    State_Met%GWETTOP  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%USTAR, &
            'USTAR', HCRC,      FIRST,    State_Met%USTAR  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%Z0, &
               'Z0', HCRC,      FIRST,    State_Met%Z0  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PARDR, &
            'PARDR', HCRC,      FIRST,    State_Met%PARDR  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PARDF, &
            'PARDF', HCRC,      FIRST,    State_Met%PARDF  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%RADSWG, &
           'RADSWG', HCRC,      FIRST,    State_Met%RADSWG  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRCLND, &
           'FRCLND', HCRC,      FIRST,    State_Met%FRCLND  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CLDFRC, &
            'CLDFRC', HCRC,      FIRST,    State_Met%CLDFRC  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

#if   defined( GEOS_4 ) 
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNOWHGT, &
          'SNOWHGT', HCRC,      FIRST,    State_Met%SNOW     )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNODP, &
            'SNODP', HCRC,      FIRST,    State_Met%SNOW   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else

    ! SNOWHGT is is mm H2O, which is the same as kg H2O/m2.
    ! This is the unit of SNOMAS.
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNOWHGT, &
          'SNOWHGT', HCRC,      FIRST,    State_Met%SNOMAS   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! SNOWDP is in m
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNODP, &
            'SNODP', HCRC,      FIRST,    State_Met%SNODP  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLAND, &
           'FRLAND', HCRC,      FIRST,    State_Met%FRLAND  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FROCEAN, &
          'FROCEAN', HCRC,      FIRST,    State_Met%FROCEAN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLAKE, &
           'FRLAKE', HCRC,      FIRST,    State_Met%FRLAKE  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLANDIC, &
         'FRLANDIC', HCRC,      FIRST,    State_Met%FRLANDIC  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! Use 'offline' LAI in standard GEOS-Chem
#if defined(ESMF_)
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%LAI, &
              'LAI', HCRC,      FIRST,    State_Met%LAI  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%LAI, &
              'LAI', HCRC,      FIRST,    GC_LAI         )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    ! ----------------------------------------------------------------
    ! 3D fields 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CNV_MFC, &
          'CNV_MFC', HCRC,      FIRST,    State_Met%CMFMC,  &
          OnLevEdge=.TRUE. )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%TK, &
               'TK', HCRC,      FIRST,    State_Met%T   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIR, &
              'AIR', HCRC,      FIRST,    State_Met%AD   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIRVOL, &
           'AIRVOL', HCRC,      FIRST,    State_Met%AIRVOL  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
 
    ! ----------------------------------------------------------------
    ! Tracer fields
    ! ----------------------------------------------------------------
    IF ( IDTO3 > 0 ) THEN
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%O3, &
            'HEMCO_O3', HCRC,      FIRST,    State_Chm%Tracers(:,:,:,IDTO3))
       IF ( HCRC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( IDTNO2 > 0 ) THEN
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%NO2, &
           'HEMCO_NO2', HCRC,      FIRST,    State_Chm%Tracers(:,:,:,IDTNO2))
       IF ( HCRC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( IDTNO > 0 ) THEN
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%NO, &
            'HEMCO_NO', HCRC,      FIRST,    State_Chm%Tracers(:,:,:,IDTNO))
       IF ( HCRC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( IDTHNO3 > 0 ) THEN
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%HNO3, &
          'HEMCO_HNO3', HCRC,      FIRST,    State_Chm%Tracers(:,:,:,IDTHNO3))
       IF ( HCRC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Deposition parameter
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%DRY_TOTN, &
         'DRY_TOTN', HCRC,      FIRST,    DRY_TOTN            ) 
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%WET_TOTN, &
         'WET_TOTN', HCRC,      FIRST,    WET_TOTN            ) 
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Other pointers to be set on first call
    ! ----------------------------------------------------------------
    IF ( FIRST ) THEN
       IF ( ExtState%WET_TOTN%DoUse .OR. ExtState%DRY_TOTN%DoUse ) THEN
          ExtState%DRYCOEFF => DRYCOEFF
       ENDIF

       ExtState%PBL_MAX => HCO_PBL_MAX
    ENDIF

    ! ----------------------------------------------------------------
    ! ESMF environment: add some additional variables to ExtState.
    ! These values must be defined here and not in the initialization
    ! because it seems like the IMPORT state is not yet properly
    ! defined during initialization. 
    ! ----------------------------------------------------------------
#if defined( ESMF_ )
    IF ( FIRST ) THEN
       CALL HCO_SetExtState_ESMF ( am_I_Root, HcoState, ExtState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP ( 'Error in HCO_SetExtState_ESMF!', LOC )
       ENDIF
    ENDIF
#endif

    ! Not first call any more
    FIRST = .FALSE.

    ! Leave with success
    RC = GIGC_SUCCESS

  END SUBROUTINE ExtState_SetFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_UpdateFields
!
! !DESCRIPTION: SUBROUTINE ExtState\_UpdateFields updates the extension 
! object data pointers. Updates are only required for the shadow arrays
! defined in this module, such as J-values, SZAFACT, etc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_UpdateFields( am_I_Root, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR

    USE PBL_MIX_MOD,           ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L

    USE FAST_JX_MOD,           ONLY : FJXFUNC
    USE COMODE_LOOP_MOD,       ONLY : NCS, JPHOTRAT, NRATES
    USE COMODE_LOOP_MOD,       ONLY : NAMEGAS, IRM

    USE HCO_GeoTools_Mod,      ONLY : HCO_GetSUNCOS
#if defined(ESMF_) 
    USE HCOI_ESMF_MOD,      ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller   - Initial Version
!  20 Aug 2014 - M. Sulprizio- Add PBL_MAX and FRAC_OF_PBL for POPs simulation
!  02 Oct 2014 - C. Keller   - PEDGE is now in HcoState%Grid
!  11 Mar 2015 - R. Yantosca - Now call GET_SZAFACT in this module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                       :: I, J, L

    ! To calculate J-values
    INTEGER                       :: K, NK, KMAX
    CHARACTER(LEN=  8)            :: SPECNAME

    CHARACTER(LEN=255), PARAMETER :: &
       LOC = 'ExtState_UpdateFields (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! ExtState_UpdateFields begins here
    !=================================================================

    ! Init
    RC = GIGC_SUCCESS

    ! TROPP: convert from hPa to Pa
    IF ( ExtState%TROPP%DoUse ) THEN
       ExtState%TROPP%Arr%Val = State_Met%TROPP * 100.0_hp
    ENDIF

    ! SPHU: convert from g/kg to kg/kg. Only need surface value. 
    IF ( ExtState%SPHU%DoUse ) THEN
       ExtState%SPHU%Arr%Val(:,:,1) = State_Met%SPHU(:,:,1) / 1000.0_hp
    ENDIF

    ! SUNCOS
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_GetSUNCOS( am_I_Root, HcoState, ExtState%SUNCOS%Arr%Val, 0, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in HCO_GetSUNCOS', LOC )
       ENDIF
    ENDIF

    ! If we need to use the SZAFACT scale factor (i.e. to put a diurnal
    ! variation on monthly mean OH concentrations), then call CALC_SUMCOSZA
    ! here.  CALC_SUMCOSZA computes the sum of cosine of the solar zenith
    ! angle over a 24 hour day, as well as the total length of daylight.
    ! This information is required by GET_SZAFACT. (bmy, 3/11/15)
    IF ( ExtState%SZAFACT%DoUse ) THEN
       CALL Calc_SumCosZa()
    ENDIF

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, K, NK, KMAX, SPECNAME )
    ! Loop over all grid boxes
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Current SZA divided by total daily SZA (2D field only)
       ! (This is mostly needed for offline simulations where a diurnal
       ! scale factor has to be imposed on monthly mean OH concentrations.)
       IF ( ExtState%SZAFACT%DoUse .AND. L==1 ) THEN
          HCO_SZAFACT(I,J) = GET_SZAFACT(I,J,State_Met)
       ENDIF

       ! Fraction of PBL for each box [unitless]
       IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
          HCO_FRAC_OF_PBL(I,J,L) = GET_FRAC_OF_PBL(I,J,L)
       ENDIF

       ! Maximum extent of the PBL [model level]
       HCO_PBL_MAX = GET_PBL_MAX_L()

       ! J-values for NO2 and O3 (2D field only)
       ! This code was moved from hcox_paranox_mod.F90 to break
       ! dependencies to GC specific code (ckeller, 07/28/14).
       IF ( L==1 .AND.                                    &
            (ExtState%JNO2%DoUse .OR. ExtState%JO1D%DoUse) ) THEN

          ! Check if sun is up
          IF ( State_Met%SUNCOS(I,J) == 0d0 ) THEN
             IF ( ExtState%JNO2%DoUse ) JNO2 = 0.0_hp
             IF ( ExtState%JO1D%DoUse ) JO1D = 0.0_hp
          ELSE
             ! Loop over photolysis reactions to find NO2, O3
             KMAX = JPHOTRAT(NCS)
             DO K = 1, KMAX 
                ! Reaction number
                NK = NRATES(NCS) + K

                ! Nae of species being photolyzed
                SPECNAME = NAMEGAS(IRM(1,NK,NCS))

                ! Check if this is NO2 or O3, store values, 1/s
                SELECT CASE ( TRIM( SPECNAME ) )
                   CASE ( 'NO2' )
                      IF ( ExtState%JNO2%DoUse ) &
                         JNO2(I,J) = FJXFUNC(I,J,1,K,1,SPECNAME)
                   CASE ( 'O3'  )
                      IF ( ExtState%JO1D%DoUse ) &
#if defined( UCX )
                      ! IMPORTANT: Need branck *2* for O1D
                      ! Branch 1 is O3P!
                      JO1D(I,J) = FJXFUNC(I,J,1,L,2,SPECNAME)
#else
                      JO1D(I,J) = FJXFUNC(I,J,1,L,1,SPECNAME)
#endif
                   CASE DEFAULT
                END SELECT
             ENDDO !K
          ENDIF

       ENDIF
    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE ExtState_UpdateFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GridEdge_Set 
!
! !DESCRIPTION: SUBROUTINE GridEdge\_Set sets the grid edge pressure values
! on the HEMCO grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GridEdge_Set ( State_Met, RC )
!
! !USES:
!
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR

! !INPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  08 Oct 2014 - C. Keller   - Initial version
!  03 Mar 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: I, J, L

    !=================================================================
    ! GridEdge_Set begins here
    !=================================================================

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Get pressure edges [Pa] and pass to HEMCO grid. 
       HcoState%Grid%PEDGE%Val(I,J,L) = State_Met%PEDGE(I,J,L) * 100.0_hp
       IF ( L==LLPAR ) THEN
          HcoState%Grid%PEDGE%Val(I,J,L+1) = State_Met%PEDGE(I,J,L+1)  & 
                                             * 100.0_hp
       ENDIF
    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE GridEdge_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetHcoSpecies 
!
! !DESCRIPTION: Subroutine SetHcoSpecies defines the HEMCO species. These 
! are typically just the GEOS-Chem tracers. Some additional species may be 
! manually added, e.g. SESQ (which is not a tracer) or individual CO2 tracers
! per emission source (for CO2 specialty sim). 
!\\
!\\
! This routine has two phases: phase 1 simply returns the number of species
! to be used by HEMCO. This is useful as this number needs to be passed to
! the HEMCO initialization call.
! Phase 2 sets the HEMCO species information in the HEMCO state object. This
! needs to be done after initialization of the HEMCO state object.
! !INTERFACE:
!
  SUBROUTINE SetHcoSpecies( am_I_Root, Input_Opt, HcoState, &
                            nSpec,     Phase,     RC         )
!
! !USES:
!
    USE GIGC_State_Chm_Mod,    ONLY : Get_Indx
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE HENRY_COEFFS,          ONLY : Get_Henry_Constant
    USE HCO_LogFile_Mod,       ONLY : HCO_SPEC2LOG
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )   :: am_I_Root
    INTEGER,          INTENT(IN   )   :: Phase 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)   :: Input_Opt  ! Input Options object
    TYPE(Hco_State),  POINTER         :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT)   :: nSpec
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  06 Mar 2015 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: nSpc
    INTEGER            :: N,  IDTLIMO
    REAL(dp)           :: K0, CR,  pKa
    CHARACTER(LEN= 31) :: ThisName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'SetHcoSpecies (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! SetHcoSpecies begins here
    !=================================================================

    !-----------------------------------------------------------------
    ! For most simulations (e.g. full-chem simulation, most of the
    ! specialty sims), just use the GEOS-Chem species definitions.
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM   .or. &
         Input_OPt%ITS_A_CH4_SIM        .or. &
         Input_Opt%ITS_A_MERCURY_SIM    .or. &
         Input_Opt%ITS_A_POPS_SIM       .or. &
         Input_Opt%ITS_A_RnPbBe_SIM     .or. &
         Input_Opt%ITS_A_TAGOX_SIM      .or. &
         Input_Opt%ITS_A_TAGCO_SIM    ) THEN

       ! Get number of model species
       nSpc = Input_Opt%N_TRACERS
  
       ! Check for SESQ: SESQ is not transported due to its short lifetime,
       ! but emissions are still calculated (in MEGAN). SESQ is only used
       ! in the SOA simulation, i.e. if LIMO is defined. Thus, add one more
       ! species here if LIMO is a model species and calculate SESQ emissions
       ! along with LIMO!
       IDTLIMO = Get_Indx('LIMO', Input_Opt%ID_TRACER, Input_Opt%TRACER_NAME )
       IF ( IDTLIMO > 0 ) THEN
          nSpc = nSpc + 1
       ENDIF
 
       ! Assign species variables
       IF ( PHASE == 2 ) THEN

          ! Verbose
          IF ( am_I_Root ) THEN
             MSG = 'Registering HEMCO species:'
             CALL HCO_MSG(MSG)
          ENDIF

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(MSG,*) 'Input species /= expected species: ', nSpec, nSpc 
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          DO N = 1, Input_Opt%N_TRACERS

             ! Model ID and species name 
             HcoState%Spc(N)%ModID      = Input_Opt%ID_TRACER(N)
             HcoState%Spc(N)%SpcName    = Input_Opt%TRACER_NAME(N)
  
             ! Get species molecular weight, emitted molecular weight,
             ! and molecule emission ratio (molecules of emitted 
             ! compound / molecules of species). For now, hardcode
             ! the species MW for some species (VOCs), since GEOS-Chem
             ! doesn't carry them.
             CALL SetSpcMW ( am_I_Root, Input_Opt, HcoState, N, N, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
 
             ! Set Henry coefficients
             CALL GET_HENRY_CONSTANT ( TRIM(HcoState%Spc(N)%SpcName), K0, CR, pKa, RC )
             HcoState%Spc(N)%HenryK0    = K0 
             HcoState%Spc(N)%HenryCR    = CR 
             HcoState%Spc(N)%HenryPKA   = PKA

             ! Write to logfile
             IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, N )
          ENDDO      
      
          ! Eventually add SESQ. This is the last entry
          IF ( IDTLIMO > 0 ) THEN
             N                           = nSpec
             HcoState%Spc(N)%ModID       = N
             HcoState%Spc(N)%SpcName     = 'SESQ'
             HcoState%Spc(N)%MW_g        = 150.0_hp
             HcoState%Spc(N)%EmMW_g      = 150.0_hp
             HcoState%Spc(N)%MolecRatio  = 1.0_hp
             HcoState%Spc(N)%HenryK0     = 0.0_hp
             HcoState%Spc(N)%HenryCR     = 0.0_hp
             HcoState%Spc(N)%HenryPKa    = 0.0_hp

             ! Write to logfile
             IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, N )
          ENDIF


          ! Add line to log-file
          IF ( am_I_Root ) CALL HCO_MSG(SEP1='-')
       ENDIF ! Phase = 2   

    !-----------------------------------------------------------------
    ! CO2 specialty simulation 
    ! For the CO2 specialty simulation, define here all tagged 
    ! tracer. This will let HEMCO calculate emissions for each
    ! tagged tracer individually. The emissions will be passed
    ! to the CO2 arrays in co2_mod.F 
    !-----------------------------------------------------------------
    ELSEIF ( Input_Opt%ITS_A_CO2_SIM ) THEN

       ! There are up to 11 tracers
       nSpc = 11 
   
       ! Set species
       IF ( PHASE == 2 ) THEN

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(MSG,*) 'Input species /= expected species: ', nSpec, nSpc 
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Henry constants are the same for all tracers
          CALL GET_HENRY_CONSTANT ( 'CO2', K0, CR, pKa, RC )
   
          ! Assign variables
          DO N = 1, nSpec 
      
             ! Define species names. These are the names that must also be 
             ! used in the HEMCO configuration file!
             SELECT CASE ( N )
   
                CASE ( 1  )
                   ThisName = 'CO2'
                CASE ( 2  ) 
                   ThisName = 'CO2ff'
                CASE ( 3  ) 
                   ThisName = 'CO2oc'
                CASE ( 4  ) 
                   ThisName = 'CO2bal'
                CASE ( 5  ) 
                   ThisName = 'CO2bb'
                CASE ( 6  ) 
                   ThisName = 'CO2bf'
                CASE ( 7  ) 
                   ThisName = 'CO2nte'
                CASE ( 8  ) 
                   ThisName = 'CO2se'
                CASE ( 9  ) 
                   ThisName = 'CO2av'
                CASE ( 10 ) 
                   ThisName = 'CO2ch'
                CASE ( 11 ) 
                   ThisName = 'CO2corr'
   
                CASE DEFAULT
                   MSG = 'Only 11 species defined for CO2 simulation!'
                   CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                   RETURN
   
             END SELECT
 
             ! Model ID and species name 
             HcoState%Spc(N)%ModID      = N 
             HcoState%Spc(N)%SpcName    = ThisName
   
             ! Molecular weights of species & emitted species.
             HcoState%Spc(N)%MW_g       = Input_Opt%Tracer_MW_G(N) 
             HcoState%Spc(N)%EmMW_g     = Input_Opt%Tracer_MW_G(N) 
   
             ! Emitted molecules per molecule of species.
             HcoState%Spc(N)%MolecRatio = 1.0_hp
   
             ! Set Henry coefficients
             HcoState%Spc(N)%HenryK0    = K0 
             HcoState%Spc(N)%HenryCR    = CR 
             HcoState%Spc(N)%HenryPKA   = PKA

             ! Write to logfile
             IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, N )
          ENDDO
          IF ( am_I_Root ) CALL HCO_MSG(SEP1='-')

       ENDIF ! Phase = 2

    !-----------------------------------------------------------------
    ! DEFAULT (RETURN W/ ERROR) 
    !-----------------------------------------------------------------
    ELSE
       MSG = 'Invalid simulation type - cannot define model species' 
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! For phase 1, pass species to output
    nSpec = nSpc

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE SetHcoSpecies 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetSpcMW
!
! !DESCRIPTION: Subroutine SetSpcMW sets the species molecular weight, the
! emitted molecular weight, and the molecular ratio of a HEMCO species. 
! The species molecular weight is the real molecular weight of the species
! of interest, whereas the emitted molecular weight is the MW of the emitted
! quantity (e.g. 12 gC/mol if species is emitted as kgC). The molecular ratio
! is the ratio of emitted molecules per molecules of species (e.g. three 
! molecules of carbon per molecule of species).
!\\
!\\
! The quantities defined here are used for unit conversion. HEMCO will 
! attempt to convert input data onto kg emitted species/m2/s (kg emitted 
! species/m3 for concentrations) based on the species properties set here,
! the data units found in the data file, and the srcUnit attribute set in
! the HEMCO configuration file. More details on input data unit conversion 
! are given in hcoio\_dataread\_mod.F90.
!\\
!\\
! Currently, GEOS-Chem does not carry the 'real' species molecular weight,
! only the emitted molecular weight. The species molecular weight therefore
! gets hardcoded here unit there is a more flexible species structure in 
! place in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetSpcMW( am_I_Root, Input_Opt, HcoState, N, M, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )   :: am_I_Root  ! Root CPU
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)   :: Input_Opt  ! Input Options object
    TYPE(Hco_State),  POINTER         :: HcoState   ! HEMCO state
    INTEGER,          INTENT(IN   )   :: N          ! species index in Input_Opt
    INTEGER,          INTENT(IN   )   :: M          ! species index in HcoState
    INTEGER,          INTENT(INOUT)   :: RC         ! Return code
!
! !REVISION HISTORY:
!  30 Mar 2015 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: ID_EMIT
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! SetSpcMW begins here
    !=================================================================

    ! Molecular weights of species & emitted species.
    HcoState%Spc(M)%EmMW_g = Input_Opt%Tracer_MW_G(N) 
   
    ! Emitted molecules per molecule of species.
    ID_EMIT = Input_Opt%ID_EMITTED(N)
    IF ( ID_EMIT <= 0 ) THEN
       HcoState%Spc(M)%MolecRatio = 1.0_hp
    ELSE
       HcoState%Spc(M)%MolecRatio = Input_Opt%TRACER_COEFF(N,ID_EMIT)
    ENDIF

    ! Set species MW for GEOS-Chem tracers transported (and emitted) 
    ! in units other than molecules of species. 
    SELECT CASE ( TRIM(Input_Opt%TRACER_NAME(N)) )
       CASE ( 'ALK4' )
          HcoState%Spc(M)%MW_g = 58.12_hp ! Butane
       CASE ( 'ISOP' )
          HcoState%Spc(M)%MW_g = 68.12_hp
       CASE ( 'ACET' )
          HcoState%Spc(M)%MW_g = 58.08_hp 
       CASE ( 'MEK' )
          HcoState%Spc(M)%MW_g = 72.11_hp
       CASE ( 'ALD2' )
          HcoState%Spc(M)%MW_g = 44.05_hp
       CASE ( 'PRPE' )
          HcoState%Spc(M)%MW_g = 42.08_hp ! Propene
       CASE ( 'C3H8' )
          HcoState%Spc(M)%MW_g = 44.10_hp
       CASE ( 'C2H6' )
          HcoState%Spc(M)%MW_g = 30.07_hp
       CASE ( 'BCPI' )
          HcoState%Spc(M)%MW_g = 12.01_hp
       CASE ( 'OCPI' )
          HcoState%Spc(M)%MW_g = 12.01_hp
       CASE ( 'BCPO' )
          HcoState%Spc(M)%MW_g = 12.01_hp
       CASE ( 'OCPO' )
          HcoState%Spc(M)%MW_g = 12.01_hp
       CASE ( 'POA1' )
          HcoState%Spc(M)%MW_g = 12.01_hp ! ??
       CASE ( 'POG1' )
          HcoState%Spc(M)%MW_g = 12.01_hp ! ??
       CASE ( 'BENZ' )
          HcoState%Spc(M)%MW_g = 78.11_hp
       CASE ( 'TOLU' )
          HcoState%Spc(M)%MW_g = 92.14_hp
       CASE ( 'XYLE' )
          HcoState%Spc(M)%MW_g = 106.16_hp
       CASE ( 'NAP' )
          HcoState%Spc(M)%MW_g = 128.27_hp ! Naphtalene
       CASE ( 'POG2' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'POA2' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'OPOG1' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'OPOG2' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'OPOA1' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'OPOA2' )
          HcoState%Spc(M)%MW_g = 12.01_hp  ! ??
       CASE ( 'MONX' )
          HcoState%Spc(M)%MW_g = 136.0_hp  ! C10H16 
       CASE ( 'C2H4' )
          HcoState%Spc(M)%MW_g = 28.05_hp  
       CASE ( 'C2H2' )
          HcoState%Spc(M)%MW_g = 26.04_hp  
       CASE ( 'MBO' )
          HcoState%Spc(M)%MW_g = 86.13_hp 
       CASE DEFAULT
          HcoState%Spc(M)%MW_g = Input_Opt%Tracer_MW_G(N) 
    END SELECT

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE SetSpcMW
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Grid 
!
! !DESCRIPTION: Subroutine Set\_Grid tells HEMCO about the grid that is being
!  used by the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Grid( am_I_Root, State_Met, HcoState, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE GRID_MOD,           ONLY : XMID,  YMID
    USE GRID_MOD,           ONLY : XEDGE, YEDGE, YSIN
    USE GRID_MOD,           ONLY : AREA_M2
    USE HCO_ARR_MOD,        ONLY : HCO_ArrInit
!
! !INPUT ARGUMENTS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT ARGUMENTS:
!
    TYPE(Hco_State),  POINTER        :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!

    !=================================================================
    ! SET_GRID begins here
    !=================================================================

    ! NOTE: for now, just copy GEOS-Chem grid, i.e. HEMCO calculations 
    ! are performed on the GEOS-Chem simulation grid. 
    ! It is possible to define a different emissions grid below. 
    ! In this case, all arrays have to be regridded when passing 
    ! them between HEMCO and GEOS-Chem (this is also true for the 
    ! met-fields used by the extensions)! 

    ! Grid dimensions
    HcoState%NX = IIPAR
    HcoState%NY = JJPAR
    HcoState%NZ = LLPAR

    ! Set pointers to grid variables
    HcoState%Grid%XMID%Val       => XMID   (:,:,1)
    HcoState%Grid%YMID%Val       => YMID   (:,:,1)
    HcoState%Grid%XEDGE%Val      => XEDGE  (:,:,1)
    HcoState%Grid%YEDGE%Val      => YEDGE  (:,:,1)
    HcoState%Grid%YSIN%Val       => YSIN   (:,:,1)
    HcoState%Grid%AREA_M2%Val    => AREA_M2(:,:,1)
    HcoState%Grid%ZSFC%Val       => State_Met%PHIS      ! Surface geopotential height
    HcoState%Grid%BXHEIGHT_M%Val => State_Met%BXHEIGHT  ! Grid box heights

    ! Allocate PEDGE. Will be updated every time step!
    CALL HCO_ArrInit( HcoState%Grid%PEDGE, HcoState%NX, HcoState%NY, HcoState%NZ+1, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE Set_Grid
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHcoState 
!
! !DESCRIPTION: Subroutine GetHcoState is a wrapper routine to connect the 
! passed pointer to the internal HcoState object. This routine can be called
! from outside of HEMCO to obtain the HcoState object (e.g. for diagnostics).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoState ( HcoStatePtr ) 
!
! !INPUT/OUTPUT ARGUMENTS:
!
    TYPE(Hco_State),    POINTER        :: HcoStatePtr  ! HEMCO state pointer
!
! !REVISION HISTORY:
!  01 Aug 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! GetHcoState begins here
    !=================================================================

    HcoStatePtr => HcoState

  END SUBROUTINE GetHcoState
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
    USE TRACERID_MOD
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
!  12 Dec 2014 - M. Yannetti - Changed real(dp) to real(hp)
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
  FUNCTION GetHcoID( name, TrcID ) RESULT ( HcoID )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: Name  ! Tracer name 
    INTEGER,          INTENT(IN   ), OPTIONAL :: TrcID ! Tracer ID 
!
! !OUTPUT PARAMETERS:
!
    INTEGER                                   :: HcoID 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  21 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Init
    HcoID = -1

    ! To get HEMCO ID by tracer ID
    IF ( PRESENT(TrcID) ) THEN
!       IF ( TrcID > 0 ) HcoID = M2HID(TrcID)%ID
       IF ( TrcID > 0 ) HcoID = TrcID 
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
  SUBROUTINE GetHcoDiagn ( am_I_Root, DiagnName, StopIfNotFound, RC, &
                           Ptr2D,     Ptr3D,     COL                  )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_DIAGN_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)           :: am_I_Root      ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN)           :: DiagnName      ! Name of diagnostics
    LOGICAL,          INTENT(IN)           :: StopIfNotFound ! Stop if diagnostics 
                                                             ! does not exist?
    INTEGER,          INTENT(IN), OPTIONAL :: COL            ! Collection Nr. 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: FLAG, ERR, LevIDx, PS
    TYPE(DiagnCont), POINTER  :: DgnCont  => NULL()

    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetHcoDiagn (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! GetHcoDiagn begins here 
    !=======================================================================

    ! Set collection number
    PS = HcoDiagnIDManual
    IF ( PRESENT(COL) ) PS = COL

    ! Get diagnostics by name. Search all diagnostics, i.e. both AutoFill
    ! and manually filled diagnostics. Also include those with a manual
    ! output interval.
    CALL Diagn_Get( am_I_Root,   .FALSE.,  DgnCont,               &
                    FLAG,        ERR,      cName=TRIM(DiagnName), &
                    AutoFill=-1, COL=PS                            )     

    ! Error checks
    IF ( ERR /= HCO_SUCCESS ) THEN
       MSG = 'Error in getting diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP ( MSG, LOC )
    ENDIF
    IF ( (FLAG /= HCO_SUCCESS) .AND. StopIfNotFound ) THEN
       MSG = 'Cannot get diagnostics for this time stamp: ' // TRIM(DiagnName)
       CALL ERROR_STOP ( MSG, LOC )
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
             MSG = 'no data defined: '// TRIM(DiagnName)
             CALL ERROR_STOP ( MSG, LOC )
          ENDIF 
  
       ! 3D pointer: must point to 3D data
       ELSEIF ( PRESENT(Ptr3D) ) THEN
          IF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             Ptr3D => DgnCont%Arr3D%Val
          ELSE
             MSG = 'no 3D data defined: '// TRIM(DiagnName)
             CALL ERROR_STOP ( MSG, LOC )
          ENDIF 

       ! Error otherwise 
       ELSE
          MSG = 'Please define output data pointer: ' // TRIM(DiagnName)
          CALL ERROR_STOP ( MSG, LOC )
       ENDIF
    ENDIF

    ! Free pointer
    DgnCont  => NULL()

    ! Leave with success 
    RC = HCO_SUCCESS

  END SUBROUTINE GetHcoDiagn 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckSettings
!
! !DESCRIPTION: Subroutine CheckSettings performs some sanity checks of the
! switches provided in the HEMCO configuration file (in combination with the
! settings specified in input.geos). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CheckSettings( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP

    USE HCO_ExtList_Mod,    ONLY : GetExtNr,  SetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt, AddExtOpt 
    USE HCO_ExtList_Mod,    ONLY : CoreNr 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  18 Feb 2015 - C. Keller   - Initial Version
!  04 Mar 2015 - R. Yantosca - Now determine if we need to read UV albedo
!                              data from the settings in input.geos
!  16 Mar 2015 - R. Yantosca - Now also toggle TOMS_SBUV_O3 based on
!                              met field type and input.geos settings
!  25 Mar 2015 - C. Keller   - Added switch for STATE_PSC (for UCX)
!EOP
!------------------------------------------------------------------------------

    ! Local variables
    INTEGER                       :: ExtNr
    LOGICAL                       :: LTMP
    LOGICAL                       :: FOUND
    CHARACTER(LEN= 31)            :: OptName
    CHARACTER(LEN=255)            :: MSG
    CHARACTER(LEN=255), PARAMETER :: LOC = 'CheckSettings (hcoi_gc_main_mod.F90'

    !=======================================================================
    ! CheckSettings begins here
    !=======================================================================

    !-----------------------------------------------------------------------
    ! If emissions shall not be used, reset all extension numbers to -999. 
    ! This will make sure that none of the extensions will be initialized 
    ! and none of the input data related to any of the extensions will be 
    ! used.  The only exception is the NON-EMISSIONS DATA.
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LEMIS ) THEN
       CALL SetExtNr( am_I_Root, -999, RC=RC )
       IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( 'SetExtNr', LOC )
    ENDIF

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #1: UV Albedoes
    !
    ! Set the UV albedo toggle according to options in input.geos.  This 
    ! will enable/disable all fields in input.geos that are  bracketed by 
    ! '+UValbedo+'.  Check first if this bracket values has been set 
    ! explicitly in the HEMCO configuration file, in which case it will
    ! not be changed.
    !
    ! UV albedoes are needed for photolysis.  Photolysis is only used in 
    ! fullchem and aerosol-only simulations that have chemistry switched on.
    ! Now search through full list of extensions (ExtNr = -999).
    !-----------------------------------------------------------------------
    CALL GetExtOpt( -999, '+UValbedo+',  OptValBool=LTMP, &
                            FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +UValbedo+', LOC )
    ENDIF

    IF ( FOUND ) THEN

       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LCHEM .AND. ( .NOT. LTMP ) ) THEN
          MSG = 'Setting +UValbedo+ in the HEMCO configuration file ' // &
                'must not be disabled if chemistry is turned on. '    // &
                'If you don`t set that setting explicitly, it will '  // &
                'be set automatically during run-time (recommended)'
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF

    ELSE

       ! If this collection is not found in the HEMCO config file, then
       ! activate it for those simulations requiring photolysis (i.e. 
       ! fullchem or aerosols), and only if chemistry is turned on.
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%LCHEM ) THEN
             OptName = '+UValbedo+ : true'
          ELSE
             OptName = '+UValbedo+ : false'
          ENDIF
       ELSE
          OptName = '+UValbedo+ : false'
       ENDIF
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +UValbedo+', LOC )
       ENDIF

    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #2: PSC STATE (for UCX) 
    !-----------------------------------------------------------------------
    CALL GetExtOpt( -999, '+STATE_PSC+', OptValBool=LTMP, &
                           FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +STATE_PSC+', LOC )
    ENDIF
    IF ( FOUND ) THEN
       IF ( Input_Opt%LUCX /= LTMP ) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Setting +STATE_PSC+ in the HEMCO configuration'
          WRITE(*,*) 'file does not agree with stratospheric chemistry'
          WRITE(*,*) 'settings in input.geos. This may be inefficient' 
          WRITE(*,*) 'and/or yield to wrong results!' 
       ENDIF
    ELSE
       IF ( Input_Opt%LUCX ) THEN
          OptName = '+STATE_PSC+ : true'
       ELSE
          OptName = '+STATE_PSC+ : false'
       ENDIF
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +STATE_PSC+', LOC )
       ENDIF
    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #3: GMI linear stratospheric chemistry
    !
    ! Set stratospheric chemistry toggle according to options in the
    ! input.geos file.  This will enable/disable all fields in the HEMCO 
    ! configuration file that are bracketed by '+LinStratChem+'.  Check 
    ! first if +LinStratChem+  has been set explicitly in the HEMCO 
    ! configuration file, in which case it will not be changed. Search
    ! through all extensions (--> ExtNr = -999).
    !-----------------------------------------------------------------------
    CALL GetExtOpt( -999, '+LinStratChem+', OptValBool=LTMP, &
                           FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +LinStratChem+', LOC )
    ENDIF

    IF ( FOUND ) THEN

       ! Print a warning if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LSCHEM /= LTMP ) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Setting +LinStratChem+ in the HEMCO configuration'
          WRITE(*,*) 'file does not agree with stratospheric chemistry'
          WRITE(*,*) 'settings in input.geos. This may be inefficient' 
          WRITE(*,*) 'and/or may yield wrong results!' 
       ENDIF

    ELSE

       ! If this collection is not found in the HEMCO config file, then
       ! activate it only if stratospheric chemistry is turned on in
       ! the input.geos file.
       IF ( Input_Opt%LSCHEM ) THEN
          OptName = '+LinStratChem+ : true'
       ELSE
          OptName = '+LinStratChem+ : false'
       ENDIF
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +LinStratChem+', LOC )
       ENDIF

    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #4: TOMS/SBUV overhead O3 columns
    !
    ! If we are using the GEOS-FP met fields, then we will not read in 
    ! the TOMS/SBUV O3 columns.  We will instead use the O3 columns from
    ! the GEOS-FP met fields.  In this case, we will toggle the 
    ! +TOMS_SBUV_O3+ collection OFF.
    !
    ! All other met fields use the TOMS/SBUV data in one way or another,
    ! so we will have to read these data from netCDF files.  In this
    ! case, toggle the +TOMS_SBUV_O3+ collection ON if wphotolysis is
    ! required (i.e. for fullchem/aerosol simulations w/ chemistry on).
    !-----------------------------------------------------------------------
    CALL GetExtOpt( -999, '+TOMS_SBUV_O3+', OptValBool=LTMP, &
                           FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +TOMS_SBUV_O3+', LOC )
    ENDIF

#if defined( GEOS_FP )
    
    ! Disable for GEOS-FP met fields, no matter what it is set to
    ! in the HEMCO configuration file.
    IF ( FOUND ) THEN
       OptName = '+TOMS_SBUV_O3+ : false'
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt GEOS-FP +TOMS_SBUV_O3+', LOC )
       ENDIF
    ENDIF

#else

    IF ( FOUND ) THEN

       ! Print a warning if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LCHEM /= LTMP ) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Setting +TOMS_SBUV_O3+ in the HEMCO configuration'
          WRITE(*,*) 'file does not agree with the chemistry settings'
          WRITE(*,*) 'in input.geos. This may be inefficient and/or' 
          WRITE(*,*) 'may yield wrong results!' 
       ENDIF

    ELSE

       ! If this collection is not found in the HEMCO config file, then
       ! activate it only for those simulations that use photolysis 
       ! (e.g. fullchem or aerosol) and only when chemistry is turned on.
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%LCHEM ) THEN
             OptName = '+TOMS_SBUV_O3+ : true'
          ELSE
             OptName = '+TOMS_SBUV_O3+ : false'
          ENDIF
       ELSE
          OptName = '+TOMS_SBUV_O3+ : false'
       ENDIF
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +Uvalbedo+', LOC )
       ENDIF
    ENDIF 

#endif

    !-----------------------------------------------------------------
    ! NON-EMISSIONS DATA #5: Ocean Hg input data (for Hg sims only)
    !
    ! If we have turned on the Ocean Mercury simulation in the
    ! input.geos file, then we will also toggle the +OCEAN_Hg+ 
    ! collection so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------
    CALL GetExtOpt( -999, '+OCEAN_Hg+', OptValBool=LTMP, &
                            FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +OCEAN_Hg+', LOC )
    ENDIF

    IF ( FOUND ) THEN
       
       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LDYNOCEAN .AND. ( .NOT. LTMP ) ) THEN
          MSG = 'Setting +UValbedo+ in the HEMCO configuration file ' // &
                'must not be disabled if chemistry is turned on. '    // &
                'If you don`t set that setting explicitly, it will '  // &
                'be set automatically during run-time (recommended)'
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
       
    ELSE

       ! If this collection is not found in the HEMCO config file, then
       ! activate it for those simulations requiring photolysis (i.e. 
       ! fullchem or aerosols), and only if chemistry is turned on.
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
          IF ( Input_Opt%LDYNOCEAN ) THEN
             OptName = '+OCEAN_Hg+ : true'
          ELSE
             OptName = '+OCEAN_Hg+ : false'
          ENDIF
       ELSE
          OptName = '+OCEAN_Hg+ : false'
       ENDIF
       CALL AddExtOpt( TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +OCEAN_Hg+', LOC )
       ENDIF

    ENDIF 

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE CheckSettings 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_szafact
!
! !DESCRIPTION:
!  Subroutine GET\_SZAFACT returns diurnal scale factors from dividing
!  the sza by the sum of the total sza per day. These factors are mainly
!  imposed to the monthly OH climatology. 
!  However, the same scale factors are dimensionless and can hence be 
!  applied to other compounds too (e.g. O3).
!\\
! !INTERFACE:
!
  FUNCTION Get_SzaFact( I, J, State_Met ) RESULT( FACT )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE Time_Mod,           ONLY : Get_TS_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J
    TYPE(MetState), INTENT(IN) :: State_Met
!
! !RETURN VALUE:
!
    REAL(fp)                   :: FACT
!
! !REMARKS:
!  Moved here from the obsolete global_oh_mod.F.
!
! !REVISION HISTORY: 
!  01 Mar 2013 - C. Keller - Imported from carbon_mod.F, where these
!                            calculations are done w/in GET_OH
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! GET_SZAFACT begins here!
    !=======================================================================

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e+0_fp  .AND. & 
         SUMCOSZA(I,J)         > 0e+0_fp ) THEN

       ! Impose a diurnal variation on OH during the day
       FACT = ( State_Met%SUNCOS(I,J) / SUMCOSZA(I,J) ) &
            *  ( 1440e+0_fp           / GET_TS_CHEM() )

    ELSE

       ! At night, OH goes to zero
       FACT = 0e+0_fp

    ENDIF

  END FUNCTION Get_SzaFact
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_sumcosza
!
! !DESCRIPTION:
!  Subroutine CALC\_SUMCOSZA computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_SumCosZa
!
! !USES:
!
    USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
    USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

    USE CMN_SIZE_MOD  ! Size parameters
    USE CMN_GCTM_MOD
!
! !REMARKS:
!  Moved here from the obsolete global_oh_mod.F.
!
! !REVISION HISTORY: 
! 01 Mar 2013 - C. Keller - Imported from carbon_mod.F, where it's
! called OHNO3TIME
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE :: SAVEDOY = -1
    INTEGER       :: I, IJLOOP, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(MAXIJ)

    !=======================================================================
    ! CALC_SUMCOSZA begins here!
    !=======================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) & 
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( SAVEDOY /= GET_DAY_OF_YEAR() ) THEN 

       ! Zero arrays
       SUMCOSZA(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()      

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day 
       DO N = 1, NDYSTEP
            
          ! Zero SUNTMP array
          SUNTMP(:) = 0e+0_fp

          ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
          DO J = 1, JJPAR
          DO I = 1, IIPAR

             ! Grid box latitude center [radians]
             YMID_R = GET_YMID_R( I, J, 1 )

             ! Increment IJLOOP
             IJLOOP = ( (J-1) * IIPAR ) + I
             TIMLOC = real(LHR0) + real(NT)/3600.0 + &
                      GET_XMID( I, J, 1 ) / 15.0
         
             DO WHILE (TIMLOC .lt. 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC .gt. 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !     
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
             !                   
             ! where LAT = the latitude angle, 
             !       DEC = the solar declination angle,  
             !       AHR = the hour angle, all in radians. 
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and 
             ! therefore does not contribute to any solar heating.  
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +          &
                              cos(YMID_R) * cos(DEC) * cos(AHR)

             ! SUMCOSZA is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             SUMCOSZA(I,J) = SUMCOSZA(I,J) +             &
                             MAX(SUNTMP(IJLOOP),0e+0_fp)

         ENDDO
         ENDDO
!!$OMP END PARALLEL DO

         ! Increment elapsed time [sec]
         NT = NT + ( GET_TS_CHEM() * 60 )             
      ENDDO

      ! Set saved day of year to current day of year 
      SAVEDOY = GET_DAY_OF_YEAR()

   ENDIF

   ! Return to calling program
 END SUBROUTINE Calc_SumCosZa
!EOC
END MODULE HCOI_GC_MAIN_MOD
