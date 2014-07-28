# if !defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_main_mod.F90 
!
! !DESCRIPTION: Module hcoi\_gc\_main\_mod.F90 is the HEMCO-to-GEOS-Chem
! interface module, providing the link between GEOS-Chem and HEMCO.
!\\
!\\
! This module contains wrapper routines to initialize, execute and finalize
! HEMCO from within GEOS-Chem. Typically, these routines are called from
! main.F.
!\\
!\\
! !INTERFACE:
!
MODULE HCOI_GC_Main_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCOX_State_Mod, ONLY : Ext_State 
  USE HCO_State_Mod,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOI_GC_Init
  PUBLIC :: HCOI_GC_Run
  PUBLIC :: HCOI_GC_Final
!
! !REMARKS:
!  This module is ignored if you are using HEMCO in an ESMF environment.
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version. 
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MODULE VARIABLES:
!
  ! HEMCO state 
  TYPE(HCO_State), POINTER        :: HcoState  => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State), POINTER        :: ExtState  => NULL()

  ! Internal met fields (will be used by some extensions)
  REAL*8, ALLOCATABLE, TARGET     :: HCO_PCENTER(:,:,:)
  REAL*8, ALLOCATABLE, TARGET     :: HCO_PEDGE  (:,:,:)
  REAL*8, ALLOCATABLE, TARGET     :: HCO_SZAFACT(:,:)

  ! Sigma coordinate (temporary)
  REAL(df), ALLOCATABLE, TARGET   :: ZSIGMA(:,:,:)

  ! Pointers used during initialization (for species matching)
  INTEGER                     :: nHcoSpec
  CHARACTER(LEN= 31), POINTER :: HcoSpecNames(:) => NULL()
  INTEGER                     :: nModelSpec
  CHARACTER(LEN= 31), POINTER :: ModelSpecNames(:) => NULL()
  INTEGER,            POINTER :: ModelSpecIDs  (:) => NULL()
  REAL(hp),           POINTER :: ModelSpecMW   (:) => NULL()
  REAL(hp),           POINTER :: ModelSpecEmMW (:) => NULL()
  REAL(hp),           POINTER :: ModelSpecMolecRatio(:) => NULL()
  REAL(hp),           POINTER :: ModelSpecK0   (:) => NULL()
  REAL(hp),           POINTER :: ModelSpecCR   (:) => NULL()
  REAL(hp),           POINTER :: ModelSpecPKA  (:) => NULL()
  INTEGER,            POINTER :: matchidx(:) => NULL()

  ! Temporary toggle for diagnostics
  LOGICAL, PARAMETER :: DoDiagn = .TRUE.

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
  SUBROUTINE HCOI_GC_Init( am_I_Root, Input_Opt, State_Chm, State_Met, RC ) 
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

    ! HEMCO routines 
    USE HCO_Config_Mod,     ONLY : Config_ReadFile
    USE HCO_State_Mod,      ONLY : HcoState_Init
    USE HCO_Driver_Mod,     ONLY : HCO_INIT
    USE HCOX_Driver_Mod,    ONLY : HCOX_INIT
    USE HCOI_GC_diagn_Mod,  ONLY : HCOI_DIAGN_INIT
    USE HCO_LogFile_Mod,    ONLY : HCO_SPEC2LOG
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: nnMatch, HMRC
    CHARACTER(LEN=255)              :: LOC

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
    !=================================================================
    CALL Config_ReadFile ( am_I_Root, Input_Opt%HcoConfigFile, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Config_ReadFile', LOC )

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LOGFILE_OPEN ( RC=HMRC ) 
       IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Open Logfile', LOC )

    ! If this is not the root CPU, always disable verbose mode.
    ELSE
       CALL HCO_VERBOSE_SET ( .FALSE. )
    ENDIF

    !=================================================================
    ! Initialize HEMCO state object and populate it 
    !=================================================================

    !-----------------------------------------------------------------
    ! Extract species to use in HEMCO 
    CALL Get_nnMatch( Input_Opt, nnMatch, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Get_nnMatch', LOC )

    !-----------------------------------------------------------------
    ! Initialize HCO state. Use only species that are used
    ! in GEOS-Chem and are also found in the HEMCO config. file.
    CALL HcoState_Init( am_I_Root, HcoState, nnMatch, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HcoState_Init', LOC )

    !-----------------------------------------------------------------
    ! Set grid
    CALL Set_Grid( am_I_Root, State_Met, RC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Set_Grid', LOC )

    !-----------------------------------------------------------------
    ! Register species
    CALL Register_Species( am_I_Root, State_Chm, RC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Register_Species', LOC )

    !=================================================================
    ! Set misc. parameter
    !=================================================================

    ! Emission, chemistry and dynamics timestep in seconds
    HcoState%TS_EMIS = GET_TS_EMIS() * 60.0
    HcoState%TS_CHEM = GET_TS_CHEM() * 60.0
    HcoState%TS_DYN  = GET_TS_DYN()  * 60.0

    ! Set ESMF flag 
    HcoState%isESMF = .FALSE.  

    ! HEMCO configuration file
    HcoState%ConfigFile = Input_Opt%HcoConfigFile

    !=================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and 
    ! the HEMCO configuration file is removed from buffer in this
    ! step. Also initializes the HEMCO clock
    !=================================================================
    CALL HCO_INIT ( am_I_Root, HcoState, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCO_INIT', LOC )

    !=================================================================
    ! Initialize extensions.
    ! This initializes all (enabled) extensions and selects all met.
    ! fields needed by them. 
    !=================================================================
    CALL HCOX_INIT ( am_I_Root, HcoState, ExtState, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCO_INIT', LOC )

    !-----------------------------------------------------------------
    ! Update logical switches in Input_Opt 
    !-----------------------------------------------------------------
    Input_Opt%LSOILNOX = ExtState%SoilNOx

    !-----------------------------------------------------------------
    ! Set pointers to met fields.
    ! Extensions typically depend on environmental dependent met. 
    ! variables such as wind speed, surface temp., etc. Pointers 
    ! to these (2D or 3D) fields are defined in the extension object. 
    ! Here, we need to make sure that these pointers are properly 
    ! connected.
    !-----------------------------------------------------------------
    CALL SET_EXTOPT_FIELDS ( State_Met, State_Chm, RC )
    IF ( RC/=GIGC_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Define diagnostics
    !-----------------------------------------------------------------
    IF ( DoDiagn ) THEN
    CALL HCOI_DIAGN_INIT ( am_I_Root, HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP('HCOI_DIAGN_INIT',LOC)
    ENDIF

    !-----------------------------------------------------------------
    ! Leave 
    !-----------------------------------------------------------------

    ! Deallocate local variables
    IF (ASSOCIATED(ModelSpecNames)) DEALLOCATE(ModelSpecNames)
    IF (ASSOCIATED(ModelSpecIDs)) DEALLOCATE(ModelSpecIDs)
    IF (ASSOCIATED(ModelSpecMW)) DEALLOCATE(ModelSpecMW)
    IF (ASSOCIATED(ModelSpecEmMW)) DEALLOCATE(ModelSpecEmMW)
    IF (ASSOCIATED(ModelSpecMolecRatio)) DEALLOCATE(ModelSpecMolecRatio)
    IF (ASSOCIATED(ModelSpecK0)) DEALLOCATE(ModelSpecK0)
    IF (ASSOCIATED(ModelSpecCR)) DEALLOCATE(ModelSpecCR)
    IF (ASSOCIATED(ModelSpecPKA)) DEALLOCATE(ModelSpecPKA)

    IF (ASSOCIATED(matchIDx    )) DEALLOCATE(matchIDx    )
    IF (ASSOCIATED(HcoSpecNames)) DEALLOCATE(HcoSpecNames)

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
  SUBROUTINE HCOI_GC_Run( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE GIGC_State_Chm_Mod,    ONLY : ChmState

    ! HEMCO routines 
    USE HCO_DIAGN_MOD,         ONLY : HCO_DIAGN_UPDATE
    USE HCO_FLUXARR_MOD,       ONLY : HCO_FluxarrReset 
    USE HCO_DRIVER_MOD,        ONLY : HCO_RUN
    USE HCOX_DRIVER_MOD,       ONLY : HCOX_RUN
    USE HCOI_GC_DIAGN_MOD,     ONLY : HCOI_DIAGN_WRITEOUT

    USE PRESSURE_MOD,          ONLY : GET_PEDGE, GET_PCENTER
    USE GLOBAL_OH_MOD,         ONLY : GET_SZAFACT

    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR

    ! For dust PBL mixing
    USE TRACERID_MOD,          ONLY : IDTDST1, IDTDST2
    USE TRACERID_MOD,          ONLY : IDTDST3, IDTDST4

    ! For soilnox
    USE GET_NDEP_MOD,          ONLY : RESET_DEP_N

    ! for temporary dust workaround 
    USE HCO_STATE_MOD,         ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Meteo state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: I, J, L
    INTEGER                   :: HMRC 
    CHARACTER(LEN=255)        :: MSG, LOC

    ! temporary only (PBL mixing of dust emissions)
    INTEGER :: HcoDST1, HcoDST2, HcoDST3, HcoDST4
    REAL*8  :: TMP

    !=================================================================
    ! HCOI_GC_RUN begins here!
    !=================================================================

    ! For error handling
    LOC = 'HCOI_GC_RUN (hcoi_gc_main_mod.F90)'

    ! Set return code flag to HCO success. This value should be
    ! preserved throughout all HCO calls, otherwise an error
    ! will be returned!
    HMRC = HCO_SUCCESS

    !=================================================================
    ! Set HcoClock 
    !=================================================================
    CALL SET_CURRENT_TIME ( am_I_Root, HcoState, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SET_CURRENT_TIME', LOC )

    !=================================================================
    ! Output diagnostics 
    !=================================================================
    IF ( DoDiagn ) THEN
    CALL HCOI_DIAGN_WRITEOUT ( am_I_Root, HcoState, .FALSE., HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_WRITEOUT', LOC )
    ENDIF

    ! ================================================================
    ! Reset all emission and deposition values
    ! ================================================================
    CALL HCO_FluxarrReset ( HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('ResetArrays', LOC )
       RETURN 
    ENDIF
 
    ! ================================================================
    ! Set HCO options and define all arrays needed by core module 
    ! and the extensions 
    ! ================================================================

    ! Range of tracers and emission categories.
    ! Set Extension number ExtNr to 0, indicating that the core
    ! module shall be executed. 
    HcoState%Options%SpcMin = 1 
    HcoState%Options%SpcMax = Input_Opt%N_Tracers 
    HcoState%Options%CatMin = 1 
    HcoState%Options%CatMax = -1 
    HcoState%Options%ExtNr  = 0

    ! Use temporary array?
    HcoState%Options%FillBuffer = .FALSE. 

    ! ================================================================
    ! Run HCO core module
    ! Emissions will be written into the corresponding flux arrays 
    ! in HcoState. 
    ! ================================================================
    CALL HCO_RUN ( am_I_Root, HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('HCO_RUN', LOC )
       RETURN 
    ENDIF

    ! ================================================================
    ! Run HCO extensions
    ! ================================================================

    ! Update HCO_PEDGE, HCO_PCENTER and HCO_SZAFACT
    IF ( ExtState%PEDGE%DoUse .OR. ExtState%PCENTER%DoUse .OR. &
         ExtState%SZAFACT%DoUse                               ) THEN

       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          HCO_PEDGE(I,J,L)   = GET_PEDGE(  I,J,L)
          HCO_PCENTER(I,J,L) = GET_PCENTER(I,J,L)
          IF ( L==1 ) HCO_SZAFACT(I,J) = GET_SZAFACT(I,J,State_Met)
          IF ( L==LLPAR ) &
                  HCO_PEDGE(I,J,LLPAR+1) = GET_PEDGE(I,J,LLPAR+1)
       ENDDO
       ENDDO
       ENDDO
    ENDIF

    ! Execute all enabled emission extensions. Emissions will be 
    ! added to corresponding flux arrays in HcoState.
    CALL HCOX_RUN ( am_I_Root, HcoState, ExtState, HMRC )
    IF ( HMRC/= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('HCOX_RUN', LOC )
       RETURN
    ENDIF 

    !=================================================================
    ! Update diagnostics 
    !=================================================================

    IF ( DoDiagn ) THEN
    CALL HCO_DIAGN_UPDATE ( am_I_Root, HcoState, HMRC )
    IF( HMRC /= HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_UPDATE', LOC )
    ENDIF

!-----------------------------------------------------------------------
    ! For now, make sure that dust emissions are just emitted
    ! into lowest layer. This is just a workaround until we
    ! have a cleaner PBL mixing implementation in place!
    ! TODO: needs nicer implementation
!-----------------------------------------------------------------------

    ! toggle (don't use in classic run)
    IF ( ( ExtState%DustDead .OR. ExtState%DustGinoux) ) THEN 

       ! Get HEMCO IDs
       HcoDST1 = HCO_GetHcoID( 'DST1', HcoState ) 
       HcoDST2 = HCO_GetHcoID( 'DST2', HcoState ) 
       HcoDST3 = HCO_GetHcoID( 'DST3', HcoState ) 
       HcoDST4 = HCO_GetHcoID( 'DST4', HcoState ) 

       DO J = 1, JJPAR
       DO I = 1, IIPAR
          IF ( IDTDST1 > 0 ) THEN
             TMP = HcoState%Spc(HcoDST1)%Emis%Val(I,J,1) * &
                   HcoState%Grid%AREA_M2(I,J) * HcoState%TS_EMIS
             State_Chm%Tracers(I,J,1,IDTDST1) =    &
             State_Chm%Tracers(I,J,1,IDTDST1) + TMP
             HcoState%Spc(HcoDST1)%Emis%Val(I,J,1) = 0d0
          ENDIF
          IF ( IDTDST2 > 0 ) THEN
             TMP = HcoState%Spc(HcoDST2)%Emis%Val(I,J,1) * &
                   HcoState%Grid%AREA_M2(I,J) * HcoState%TS_EMIS
             State_Chm%Tracers(I,J,1,IDTDST2) =    &
             State_Chm%Tracers(I,J,1,IDTDST2) + TMP
             HcoState%Spc(HcoDST2)%Emis%Val(I,J,1) = 0d0
          ENDIF
          IF ( IDTDST3 > 0 ) THEN
             TMP = HcoState%Spc(HcoDST3)%Emis%Val(I,J,1) * &
                   HcoState%Grid%AREA_M2(I,J) * HcoState%TS_EMIS
             State_Chm%Tracers(I,J,1,IDTDST3) =    &
              State_Chm%Tracers(I,J,1,IDTDST3) + TMP
             HcoState%Spc(HcoDST3)%Emis%Val(I,J,1) = 0d0
          ENDIF
          IF ( IDTDST4 > 0 ) THEN
             TMP = HcoState%Spc(HcoDST4)%Emis%Val(I,J,1) * &
                   HcoState%Grid%AREA_M2(I,J) * HcoState%TS_EMIS
             State_Chm%Tracers(I,J,1,IDTDST4) =    &
             State_Chm%Tracers(I,J,1,IDTDST4) + TMP
             HcoState%Spc(HcoDST4)%Emis%Val(I,J,1) = 0d0
          ENDIF
       ENDDO
       ENDDO
       ! end dust mixing 
    ENDIF

    ! ================================================================
    ! Translate emissions array from HCO state onto GC arrays
    ! This step also converts emissions from kg/m2/s to molec/cm2/s!
    ! ================================================================
    CALL MAP_HCO2GC ( HcoState, Input_Opt, State_Chm, RC )

    ! Reset deposition arrays  
    ! TODO: Do somewhere else? e.g. in drydep/wetdep routines?
    CALL RESET_DEP_N

    ! We are now back in GEOS-Chem environment, hence set 
    ! return flag accordingly! 
    !      first = .false.
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
  SUBROUTINE HCOI_GC_Final( am_I_Root )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE Error_Mod,           ONLY : Error_Stop
    USE CMN_SIZE_Mod,        ONLY : IIPAR, JJPAR, LLPAR
    USE HCO_Driver_Mod,      ONLY : HCO_Final
    USE HCO_Diagn_Mod,       ONLY : Diagn_Cleanup
    USE HCO_State_Mod,       ONLY : HcoState_Final
    USE HCOI_GC_Diagn_Mod,   ONLY : HCOI_Diagn_WriteOut
    USE HCOX_Driver_Mod,     ONLY : HCOX_Final
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)     :: am_I_Root
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, HMRC
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! HCOI_GC_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_GC_Final (hcoi_gc_main_mod.F90)'

    ! Write out 'standard' diagnostics
    CALL HCOI_DIAGN_WRITEOUT ( am_I_Root, HcoState, .FALSE., HMRC, &
                               UsePrevTime=.FALSE. )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCOI_DIAGN_FINAL', LOC )
 
    ! Also write all other diagnostics into netCDF file
    CALL HCOI_DIAGN_WRITEOUT ( am_I_Root, HcoState, .TRUE., HMRC, &
                               PREFIX='HEMCO_restart', UsePrevTime=.FALSE.)
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCOI_DIAGN_FINAL', LOC )
 
    ! Cleanup diagnostics
    CALL Diagn_Cleanup()
 
    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    CALL HCOX_FINAL ( ExtState ) 

    ! Cleanup HCO core
    CALL HCO_FINAL()
  
    ! Remove emission array pointers.
    ! Note: Only need to do this if HEMCO grid is equal to GC 
    ! simulation grid. Otherwise, nothing to do here (arrays
    ! will be explicitly deallocated in HCO\_FINAL call below).
    IF ( ( HcoState%NX == IIPAR ) .AND.      &
         ( HcoState%NY == JJPAR ) .AND.      &
         ( HcoState%NZ == LLPAR )      ) THEN 
       DO I = 1, HcoState%nSpc
          HcoState%Spc(I)%Emis%Val => NULL()
          HcoState%Spc(I)%Depv%Val => NULL()
       ENDDO
    ENDIF

    ! Make sure HcoState variables pointing to shared data are 
    ! nullified (otherwise, the target arrays become deallocated)
    HcoState%Grid%XMID       => NULL()
    HcoState%Grid%YMID       => NULL()
    HcoState%Grid%ZSIGMA     => NULL()
    HcoState%Grid%XEDGE      => NULL()
    HcoState%Grid%YEDGE      => NULL()
    HcoState%Grid%YSIN       => NULL()
    HcoState%Grid%AREA_M2    => NULL()
    HcoState%Grid%BXHEIGHT_M => NULL()

    ! Cleanup HcoState object 
    CALL HcoState_Final ( HcoState ) 

    ! Module variables
    IF ( ALLOCATED  ( ZSIGMA      ) ) DEALLOCATE ( ZSIGMA      )
    IF ( ALLOCATED  ( HCO_PEDGE   ) ) DEALLOCATE ( HCO_PEDGE   )
    IF ( ALLOCATED  ( HCO_PCENTER ) ) DEALLOCATE ( HCO_PCENTER )
    IF ( ALLOCATED  ( HCO_SZAFACT ) ) DEALLOCATE ( HCO_SZAFACT )

  END SUBROUTINE HCOI_GC_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Map_Hco2Gc
!
! !DESCRIPTION: Function MAP\_Hco2Gc fills the GEOS-Chem tracer
! tendency array (in chemistry state) with the corresponding emission
! values of the HCO state object. Regridding is performed if
! necessary. Emissions are kept in units of kg/m2/s. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Map_Hco2Gc( HcoState, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE HCO_STATE_MOD,         ONLY : HCO_State

    USE CMN_SIZE_MOD,          ONLY : IIPAR,  JJPAR,  LLPAR
    USE ERROR_MOD,             ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState    ! HCO state
    TYPE(OptInput),  INTENT(IN   )  :: Input_Opt   ! Input options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(INOUT)  :: State_Chm   ! Chemistry state 
    INTEGER,         INTENT(INOUT)  :: RC          ! Failure?
!
! !REVISION HISTORY:
!  01 May 2012 - C. Keller - Initial Version
!  20 Aug 2013 - C. Keller - Now pass from HEMOC state to chemistry state
!  14 Sep 2013 - C. Keller - Now keep in units of kg/m2/s.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                  :: N, nSpc, trcID
    REAL*8                   :: COEFF 
    REAL*8, PARAMETER        :: N_0 = 6.022d+23
    CHARACTER(LEN=255)       :: MSG, LOC

    !=================================================================
    ! MAP_HCO2GC begins here
    !=================================================================

    ! For error handling
    LOC = 'Map_HCO2GC (hcoi_gc_main_mod.F90)'

    ! Assume success until otherwise
    RC = GIGC_SUCCESS

    ! Number of HEMCO species
    nSpc = HcoState%nSpc

    ! Loop over all specified tracers 
    DO N = 1, nSpc 

       ! GC tracer ID
       trcID = HcoState%Spc(N)%ModID

       ! Skip if tracer ID is not defined
       IF ( trcID <= 0 ) CYCLE

       ! Get conversion coefficient from kg/m2/s to molec/cm2/s
       ! --> kg to molec: [kg] * [g/kg] / [g/mol] * [molec/mol]
       !                : [kg] * 1000   / MW      * N_0
       ! --> m2 to cm2  : [m2] * [cm2/m2]
       !                : [m2] * 10000
       ! --> kg/m2/s to molec/cm2/s: 
       !     [kg/m2/s] * 1000 / MW * N_0 / 10000
       !     [kg/m2/s] * 0.1 / MW * N_0
!       COEFF = 0.1d0 / HcoState%Spc(N)%EmMW_g * N_0

       ! If simulation grid and emission grid are equal, the 
       ! HEMCO state emission array already points to
       ! State_Chm%Trac_Tend and there is nothing to do here. 
       IF ( (HcoState%NX == IIPAR) .AND.       &
            (HcoState%NY == JJPAR) .AND.       &
            (HcoState%NZ == LLPAR)       ) THEN

          ! ... nothing to do here        
 
       ! For different grids, regrid emissions onto simulation
       ! grid
       ! TODO: needs testing
       ELSE
          ! Do regridding
          CALL REGRID_EMIS2SIM ( HcoState, N, State_Chm, trcID, Input_Opt ) 
       ENDIF

!      ! Convert to molec/cm2/s
!      State_Chm%Trac_Tend(:,:,:,trcID) = &
!         State_Chm%Trac_Tend(:,:,:,trcID) * COEFF

    ENDDO !N 

  END SUBROUTINE Map_HCO2GC
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Regrid_Emis2Sim
!
! !DESCRIPTION: Function REGRID\_Emis2Sim regrids the original emission
!  field onto the simulation grid. Emissions are in kg/m2/s, i.e. no
!  area adjustment is applied.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Regrid_Emis2Sim( HcoState, hcoID, State_Chm, trcID, Input_Opt ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE HCO_STATE_MOD,         ONLY : HCO_State 
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE REGRID_A2A_MOD,        ONLY : MAP_A2A
    USE GRID_MOD,              ONLY : XEDGE, YSIN
    USE DRYDEP_MOD,            ONLY : NUMDEP, NTRAIND
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState
    INTEGER,         INTENT(IN   )  :: hcoID 
    INTEGER,         INTENT(IN   )  :: trcID 
    TYPE(OptInput),  INTENT(IN   )  :: Input_Opt  ! Input options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(INOUT)  :: State_Chm        ! Chemistry state 
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 Jan 2013 - C. Keller - Now call MAP_A2A instead of DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: II, JJ, LL
    INTEGER            :: NX, NY
    INTEGER            :: NN, DepIdx
    INTEGER            :: nLevIn, nLevOut
    REAL*8, POINTER    :: ORIG_2D(:,:) => NULL()
    REAL*8, POINTER    :: REGR_2D(:,:) => NULL()
    REAL*8, POINTER    :: InArr(:,:,:) => NULL()
    REAL*8, POINTER    :: OutArr(:,:,:) => NULL()
    CHARACTER(LEN=255) :: LOC
    REAL*8                   :: XEDGE_IN ( HcoState%NX+1)
    REAL*8                   :: YSIN_IN  ( HcoState%NY+1)
    REAL*8                   :: XEDGE_OUT( IIPAR+1)
    REAL*8                   :: YSIN_OUT ( JJPAR+1)

    !=================================================================
    ! REGRID_EMIS2SIM begins here
    !=================================================================

    ! For error handling
    LOC = 'Regrid_Emis2Sim ("hcoi_gc_main_mod.F90")'

    ! Pass horizontal emission grid dimensions
    NX = HcoState%NX
    NY = HcoState%NY

    ! Get grid edges on input grid
    XEDGE_IN(:) = HcoState%Grid%XEDGE(:,1)
    YSIN_IN (:) = HcoState%Grid%YSIN (1,:)

    ! Get grid edges on output grid
    XEDGE_OUT(:) = XEDGE(:,1,1)
    YSIN_OUT (:) = YSIN (1,:,1)

    ! Point to tracer tendency slice
    InArr  => HcoState%Spc(hcoID)%Emis%Val
    OutArr => State_Chm%Trac_Tend(:,:,:,trcID) 

    ! Only do for existing (i.e. associated) emission arrays
    IF ( ASSOCIATED( InArr ) ) THEN
   
       ! Get vertical levels
       nLevIn  = SIZE(InArr,3)
       nLevOut = SIZE(OutArr,3)
   
       ! Output array must not have less vertical levels than input
       ! array.
       IF ( nLevOut < nLevIn ) THEN
          CALL ERROR_STOP ( 'nLevIn < nLevOut', LOC )
       ENDIF
   
       ! Loop over all vertical levels of the input array
       DO LL = 1, nLevIn
   
          ! No need to regrid if all emissions in this
          ! layer and for this tracer are zero!
          IF ( SUM(InArr(:,:,LL)) == 0d0 ) THEN
             CYCLE
          ENDIF
          
          ! Point to 2D slice to be regridded
          ORIG_2D => InArr (:,:,LL)
          REGR_2D => OutArr(:,:,LL)

          ! Do the regridding
          CALL MAP_A2A( NX, NY, XEDGE_IN, YSIN_IN,  ORIG_2D, &
                        IIPAR, JJPAR, XEDGE_OUT,YSIN_OUT, REGR_2D, &
                        0, 0 )

          ! Clear pointers 
          ORIG_2D => NULL()
          REGR_2D => NULL()

       ENDDO !LL
    ENDIF

    ! Free pointer
    InArr  => NULL()               
    OutArr => NULL()

    ! Eventually also regrid deposition array
    ! --> deposition is in m/s, i.e. no mass-conservative regridding is required.
    IF ( Input_Opt%LDRYD ) THEN
       ORIG_2D => HcoState%Spc(hcoID)%Depv%Val
       IF ( ASSOCIATED(ORIG_2D) ) THEN
   
          ! Find drydep index
          DepIdx = -1
          DO NN = 1, NUMDEP
             IF ( NTRAIND(NN) == trcID ) THEN
                DepIdx = NN
                EXIT
             ENDIF
          ENDDO
   
          IF ( DepIdx > 0 ) THEN
             REGR_2D => State_Chm%DepSav(:,:,DepIdx)
   
             ! Do the regridding
             CALL MAP_A2A( NX, NY, XEDGE_IN, YSIN_IN,  ORIG_2D, &
                           IIPAR, JJPAR, XEDGE_OUT,YSIN_OUT, REGR_2D, &
                           0, 0 )
          ENDIF
   
          ! Clear pointers 
          ORIG_2D => NULL()
          REGR_2D => NULL()
       ENDIF
    ENDIF

  END SUBROUTINE  Regrid_Emis2Sim
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Current_Time 
!
! !DESCRIPTION: SUBROUTINE Set\_Current\_Time sets the current simulation 
! datetime in HcoState. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Current_Time( am_I_Root, HcoState, RC ) 
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
    TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 Jan 2013 - C. Keller - Now call MAP_A2A instead of DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: cYr, cMt, cDy, cHr, cMin, cSec, cDOY 

    !=================================================================
    ! SET_CURRENT_TIME begins here
    !=================================================================

    cYr      = GET_YEAR()
    cMt      = GET_MONTH()
    cDy      = GET_DAY()
    cHr      = GET_HOUR()
    cMin     = GET_MINUTE()
    cSec     = GET_SECOND()
    cDOY     = GET_DAY_OF_YEAR()

    CALL HcoClock_Set ( am_I_Root,  HcoState, cYr, cMt, cDy, cHr, &
                        cMin, cSec, cDoy, RC=RC )

  END SUBROUTINE Set_Current_Time
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_ExtOpt_Fields
!
! !DESCRIPTION: SUBROUTINE Set\_ExtOpt\_Fields sets the extension object data
! pointers. This routine must be called after HCOX\_Init. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_ExtOpt_Fields( State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState

    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR

    ! For ParaNOx
    USE TRACERID_MOD,       ONLY : IDTO3, IDTNO, IDTNO2, IDTHNO3

    ! For SoilNox
    USE Drydep_Mod,         ONLY : DRYCOEFF
    USE Get_Ndep_Mod,       ONLY : DRY_TOTN
    USE Get_Ndep_Mod,       ONLY : WET_TOTN

    ! MODIS variables (used by MEGAN & Soil NOx) 
    USE Modis_LAI_Mod,      ONLY : DAYS_BTW_MON
    USE Modis_LAI_Mod,      ONLY : GC_LAI
    USE Modis_LAI_Mod,      ONLY : GC_LAI_PM 
    USE Modis_LAI_Mod,      ONLY : GC_LAI_CM
    USE Modis_LAI_Mod,      ONLY : GC_LAI_NM
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller - Initial Version
!  23 Jan 2013 - C. Keller - Now call MAP_A2A instead of DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL :: DoDryCoeff
    INTEGER :: AS
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! SET_EXTOPT_FIELDS begins here
    !=================================================================

    ! Init
    RC = GIGC_SUCCESS

    LOC = 'Set_ExtOpt_Fields (hcoi_gc_main_mod.F90)'

    ! ----------------------------------------------------------------
    ! HCO_PEDGE, HCO_PCENTER and HCO_SZAFACT aren't defined as 3D 
    ! arrays in Met_State. Hence need to construct here so that we 
    ! can point to them.
    ! ----------------------------------------------------------------
    IF ( ExtState%PEDGE%DoUse .OR. ExtState%PCENTER%DoUse .OR. &
         ExtState%SZAFACT%DoUse                               ) THEN

       ALLOCATE(HCO_PEDGE  (IIPAR,JJPAR,LLPAR+1),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_PEDGE', LOC )
       HCO_PEDGE = 0d0

       ALLOCATE(HCO_PCENTER(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_PCENTER', LOC )
       HCO_PCENTER = 0d0

       ALLOCATE(HCO_SZAFACT(IIPAR,JJPAR      ),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_SZAFACT', LOC )
       HCO_SZAFACT = 0d0

       ExtState%PEDGE%Arr%Val => HCO_PEDGE
       ExtState%PCENTER%Arr%Val => HCO_PCENTER
       ExtState%SZAFACT%Arr%Val => HCO_SZAFACT
    ENDIF

    ! ----------------------------------------------------------------
    ! All other met fields: point to corresponding state-met field
    ! ----------------------------------------------------------------

    ! ----------------------------------------------------------------
    ! 2D fields
    IF ( ExtState%U10M%DoUse ) THEN
       ExtState%U10M%Arr%Val => State_Met%U10M
    ENDIF
    IF ( ExtState%V10M%DoUse ) THEN
       ExtState%V10M%Arr%Val => State_Met%V10M
    ENDIF
    IF ( ExtState%ALBD%DoUse ) THEN
       ExtState%ALBD%Arr%Val => State_Met%ALBD
    ENDIF
    IF ( ExtState%WLI%DoUse ) THEN
       ExtState%WLI%Arr%Val => State_Met%LWI
    ENDIF
    IF ( ExtState%TSURFK%DoUse ) THEN
       ExtState%TSURFK%Arr%Val => State_Met%TS
    ENDIF
    IF ( ExtState%TSKIN%DoUse ) THEN
       ExtState%TSKIN%Arr%Val => State_Met%TSKIN
    ENDIF
    IF ( ExtState%GWETTOP%DoUse ) THEN
       ExtState%GWETTOP%Arr%Val => State_Met%GWETTOP
    ENDIF
    IF ( ExtState%USTAR%DoUse ) THEN
       ExtState%USTAR%Arr%Val => State_Met%USTAR
    ENDIF
    IF ( ExtState%Z0%DoUse ) THEN
       ExtState%Z0%Arr%Val => State_Met%Z0
    ENDIF
    IF ( ExtState%TROPP%DoUse ) THEN
       ExtState%TROPP%Arr%Val => State_Met%TROPP
    ENDIF
    IF ( ExtState%SUNCOSmid%DoUse ) THEN
       ExtState%SUNCOSmid%Arr%Val => State_Met%SUNCOSmid
    ENDIF
    IF ( ExtState%SUNCOSmid5%DoUse ) THEN
       ExtState%SUNCOSmid5%Arr%Val => State_Met%SUNCOSmid5
    ENDIF
    IF ( ExtState%PARDR%DoUse ) THEN
       ExtState%PARDR%Arr%Val => State_Met%PARDR
    ENDIF
    IF ( ExtState%PARDF%DoUse ) THEN
       ExtState%PARDF%Arr%Val => State_Met%PARDF
    ENDIF
    IF ( ExtState%RADSWG%DoUse ) THEN
       ExtState%RADSWG%Arr%Val => State_Met%RADSWG
    ENDIF
    IF ( ExtState%PSURF%DoUse ) THEN
       ExtState%PSURF%Arr%Val => State_Met%PSC2
    ENDIF
    IF ( ExtState%FRCLND%DoUse ) THEN
       ExtState%FRCLND%Arr%Val => State_Met%FRCLND
    ENDIF
    IF ( ExtState%CLDFRC%DoUse ) THEN
       ExtState%CLDFRC%Arr%Val => State_Met%CLDFRC
    ENDIF
    IF ( ExtState%CLDTOPS%DoUse ) THEN
       ExtState%CLDTOPS%Arr%Val => State_Met%CLDTOPS
    ENDIF
    IF ( ExtState%SNOWHGT%DoUse ) THEN

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
       ExtState%SNOWHGT%Arr%Val => State_Met%SNOMAS
#else
       ExtState%SNOWHGT%Arr%Val => State_Met%SNOW
#endif
    ENDIF

    ! ----------------------------------------------------------------
    ! Modis LAI parameter
    IF ( ExtState%GC_LAI%DoUse ) THEN
       ExtState%GC_LAI%Arr%Val => GC_LAI
    ENDIF
    IF ( ExtState%GC_LAI_PM%DoUse ) THEN
       ExtState%GC_LAI_PM%Arr%Val => GC_LAI_PM
    ENDIF
    IF ( ExtState%GC_LAI_CM%DoUse ) THEN
       ExtState%GC_LAI_CM%Arr%Val => GC_LAI_CM
    ENDIF
    IF ( ExtState%GC_LAI_NM%DoUse ) THEN
       ExtState%GC_LAI_NM%Arr%Val => GC_LAI_NM
    ENDIF
    ExtState%DAYS_BTW_M => DAYS_BTW_MON
    
    ! 3D fields
    IF ( ExtState%SPHU%DoUse ) THEN
       ExtState%SPHU%Arr%Val => State_Met%SPHU
    ENDIF
    IF ( ExtState%TK%DoUse ) THEN
       ExtState%TK%Arr%Val => State_Met%T
    ENDIF
    IF ( ExtState%AIR%DoUse ) THEN
       ExtState%AIR%Arr%Val => State_Met%AD
    ENDIF
    IF ( ExtState%AIRVOL%DoUse ) THEN
       ExtState%AIRVOL%Arr%Val => State_Met%AIRVOL
    ENDIF
    
    ! ----------------------------------------------------------------
    ! Tracer fields
    IF ( ExtState%O3%DoUse ) THEN
       IF ( IDTO3 <= 0 ) THEN
          CALL ERROR_STOP ( 'IDTO3 need to be positive!', LOC )
       ELSE
          ExtState%O3%Arr%Val => State_Chm%Tracers(:,:,:,IDTO3)
       ENDIF
    ENDIF
    IF ( ExtState%NO2%DoUse ) THEN
       IF ( IDTNO2 <= 0 ) THEN
          CALL ERROR_STOP ( 'IDTNO2 need to be positive!', LOC )
       ELSE
          ExtState%NO2%Arr%Val => State_Chm%Tracers(:,:,:,IDTNO2)
       ENDIF
    ENDIF
    IF ( ExtState%NO%DoUse ) THEN
       IF ( IDTNO <= 0 ) THEN
          CALL ERROR_STOP ( 'IDTNO need to be positive!', LOC )
       ELSE
          ExtState%NO%Arr%Val => State_Chm%Tracers(:,:,:,IDTNO)
       ENDIF
    ENDIF
    IF ( ExtState%HNO3%DoUse ) THEN
       IF ( IDTHNO3 <= 0 ) THEN
          CALL ERROR_STOP ( 'IDTHNO3 need to be positive!', LOC )
       ELSE
          ExtState%HNO3%Arr%Val => State_Chm%Tracers(:,:,:,IDTHNO3)
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Deposition parameter
    DoDryCoeff = .FALSE.
    IF ( ExtState%DRY_TOTN%DoUse ) THEN
       ExtState%DRY_TOTN%Arr%Val => DRY_TOTN 
       DoDryCoeff = .TRUE.
    ENDIF
    IF ( ExtState%WET_TOTN%DoUse ) THEN
       ExtState%WET_TOTN%Arr%Val => WET_TOTN
       DoDryCoeff = .TRUE.
    ENDIF
    IF ( DoDryCoeff ) ExtState%DRYCOEFF => DRYCOEFF

  END SUBROUTINE Set_ExtOpt_Fields
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Model_GetSpecies 
!
! !DESCRIPTION: Subroutine Model\_GetSpecies returns 'model' species 
! information from the HEMCO standalone input file. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Model_GetSpecies( Input_Opt,                          &
                               nModelSpec,     ModelSpecNames,     &
                               ModelSpecIDs,   ModelSpecMW,        &
                               ModelSpecEmMW,  ModelSpecMolecRatio,&
                               ModelSpecK0,    ModelSpecCR,        &
                               ModelSpecPKA,   RC                   )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE Henry_Mod,             ONLY : Get_Henry_Constants
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT) :: Input_Opt  ! Input Options object
!
! !OUPTUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: nModelSpec
    CHARACTER(LEN= 31), POINTER     :: ModelSpecNames(:)
    INTEGER,            POINTER     :: ModelSpecIDs  (:)
    REAL(hp),           POINTER     :: ModelSpecMW   (:)
    REAL(hp),           POINTER     :: ModelSpecEmMW (:)
    REAL(hp),           POINTER     :: ModelSpecMolecRatio(:)
    REAL(hp),           POINTER     :: ModelSpecK0   (:)
    REAL(hp),           POINTER     :: ModelSpecCR   (:)
    REAL(hp),           POINTER     :: ModelSpecPKA  (:)
    INTEGER,            INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller   - Initial Version
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER   :: N, ID_EMIT
    REAL(dp)  :: K0, CR, pKa

    !=================================================================
    ! Model_GetSpecies begins here
    !=================================================================

    ! # of model species
    nModelSpec = Input_Opt%N_TRACERS

    ! Allocate
    ALLOCATE(ModelSpecNames     (nModelSpec))
    ALLOCATE(ModelSpecIDs       (nModelSpec))
    ALLOCATE(ModelSpecMW        (nModelSpec))
    ALLOCATE(ModelSpecEmMW      (nModelSpec))
    ALLOCATE(ModelSpecMolecRatio(nModelSpec))
    ALLOCATE(ModelSpecK0        (nModelSpec))
    ALLOCATE(ModelSpecCR        (nModelSpec))
    ALLOCATE(ModelSpecPKA       (nModelSpec))

    ! Assign variables
    DO N = 1, nModelSpec

       ! Species names
       ModelSpecNames(N) = Input_Opt%TRACER_NAME(N)
       
       ! Species ID
       ModelSpecIDs(N)   = Input_Opt%ID_TRACER(N)

       ! Molecular weights
       ModelSpecMW(N)    = Input_Opt%Tracer_MW_G(N)
       ModelSpecEmMW(N)  = Input_Opt%Tracer_MW_G(N)

       ! Emitted molecules per molecule of species
       ID_EMIT = Input_Opt%ID_EMITTED(N)
       IF ( ID_EMIT <= 0 ) THEN
          ModelSpecMolecRatio(N) = 1.0_hp
       ELSE
          ModelSpecMolecRatio(N) = Input_Opt%TRACER_COEFF(N,ID_EMIT)
       ENDIF

       ! Henry coefficients
       CALL GET_HENRY_CONSTANTS ( N, K0, CR, pKa, RC )
       ModelSpecK0(N)  = K0
       ModelSpecCR(N)  = CR
       ModelSpecPKA(N) = PKA

    ENDDO      

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE Model_GetSpecies 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Set_Grid( am_I_Root, State_Met, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE GRID_MOD,           ONLY : XMID,  YMID
    USE GRID_MOD,           ONLY : XEDGE, YEDGE, YSIN
    USE GRID_MOD,           ONLY : AREA_M2
!
! !INPUT ARGUMENTS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT ARGUMENTS:
!
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
    INTEGER :: I

    !=================================================================
    ! SET_GRID begins here
    !=================================================================

    ! NOTE: for now, just copy GEOS-Chem grid, i.e. HEMCO calculations 
    ! are performed on the GEOS-Chem simulation grid. 
    ! It is possible to define a different emissions grid below. 
    ! In this case, all arrays have to be regridded when passing 
    ! them between HEMCO and GEOS-Chem (this is also true for the 
    ! met-fields used by the extensions)! 

    ! ZSIGMA (temporary, for testing)
    ALLOCATE(ZSIGMA(1,1,LLPAR+1))
    DO I = 1,LLPAR+1
       ZSIGMA(:,:,I) = I
    ENDDO

    ! Grid dimensions
    HcoState%NX = IIPAR
    HcoState%NY = JJPAR
    HcoState%NZ = LLPAR

    ! Set pointers to grid variables
    HcoState%Grid%XMID       => XMID   (:,:,1)
    HcoState%Grid%YMID       => YMID   (:,:,1)
    HcoState%Grid%ZSIGMA     => ZSIGMA (:,:,:)
    HcoState%Grid%XEDGE      => XEDGE  (:,:,1)
    HcoState%Grid%YEDGE      => YEDGE  (:,:,1)
    HcoState%Grid%YSIN       => YSIN   (:,:,1)
    HcoState%Grid%AREA_M2    => AREA_M2(:,:,1)
    HcoState%Grid%BXHEIGHT_M => State_Met%BXHEIGHT

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE Set_Grid
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nnMatch 
!
! !DESCRIPTION: Subroutine Get\_nnMatch returns the number of HEMCO species
! that are also used in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_nnMatch( Input_Opt, nnMatch, RC ) 
!
! !USES:
!
    USE HCO_CharTools_Mod,  ONLY : HCO_CharMatch
    USE GIGC_State_Chm_Mod, ONLY : Get_Indx
    USE HCO_Config_MOD,     ONLY : Config_GetnSpecies
    USE HCO_Config_MOD,     ONLY : Config_GetSpecNames
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS
!
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt  ! Input Options object
    INTEGER,        INTENT(INOUT)  :: RC         ! Success or fialure
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(  OUT)  :: nnMatch    ! # of HEMCO species that
                                                 ! are also GEOS-Chem species
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller   - Initial Version
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: AS, IDX
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! Get_nnMatch begins here
    !=================================================================

    ! For error handling
    LOC = 'Get_nnMatch (hcoi_gc_main_mod.F90)'

    ! Extract number of HEMCO species and corresponding species names 
    ! as read from the HEMCO config. file.
    nHcoSpec = Config_GetnSpecies ( )
 
    CALL Config_GetSpecNames( HcoSpecNames, nHcoSpec, RC )
    IF( RC /= HCO_SUCCESS) RETURN 

    ! Extract species to be used from input file
    CALL Model_GetSpecies( Input_Opt,                           &
                           nModelSpec,     ModelSpecNames,      &
                           ModelSpecIDs,   ModelSpecMW,         &
                           ModelSpecEmMW,  ModelSpecMolecRatio, &
                           ModelSpecK0,    ModelSpecCR,         &
                           ModelSpecPKA,   RC                    )
    IF ( RC /= HCO_SUCCESS) RETURN

    ! See how many species are also used in GEOS-Chem
    ALLOCATE(matchIDx(nHcoSpec),STAT=AS)
    IF ( AS/=0 ) THEN 
       CALL HCO_ERROR ('Allocation error matchIDx', RC, THISLOC=LOC )
       RETURN
    ENDIF
    matchIDx(:) = -1
    CALL HCO_CharMatch( HcoSpecNames,   nHcoSpec,      &
                        ModelSpecNames, nModelSpec,    &
                        matchIDx,       nnMatch         )
    IF ( nnMatch == 0 ) THEN
       CALL HCO_ERROR ('No matching species!', RC, THISLOC=LOC )
       RETURN
    ENDIF

!=============================================================================
    ! KLUDGE for SESQ: SESQ is not transported due to its short lifetime,
    ! but emissions are still calculated (in MEGAN). SESQ is only used
    ! in the SOA simulation, i.e. if LIMO is defined. Thus, add one more
    ! species here if LIMO is a model species and calculate SESQ emissions
    ! along with LIMO!
    IDX = Get_Indx('LIMO', ModelSpecIDs, ModelSpecNames)
    IF ( IDX > 0 ) THEN
       nnMatch = nnMatch + 1 
    ENDIF
!=============================================================================

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Get_nnMatch 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Species 
!
! !DESCRIPTION: Subroutine Register\_Species registers all emissions
!  species in the HEMCO state object. 
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE Register_Species( am_I_Root, State_Chm, RC )
!
! !USES:
!
    USE HCO_LogFile_Mod,    ONLY : HCO_SPEC2LOG
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR

    ! For SOA mechanism
    USE CARBON_MOD,         ONLY : BIOG_SESQ
!
! !INPUT ARGUMENTS:
!
    LOGICAL,            INTENT(IN   )  :: am_I_Root
    TYPE(ChmState),     INTENT(IN   )  :: State_Chm  ! Chem state
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,            INTENT(INOUT)  :: RC
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
    INTEGER     :: CNT, I, IDX

    !=================================================================
    ! REGISTER_SPECIES begins here
    !=================================================================

    ! Loop over all possible HEMCO species
    cnt = 0 
    DO I = 1, nHcoSpec

       ! Skip if this HEMCO species is not used in GEOS-Chem
       IF ( MatchIDx(I) < 0 ) CYCLE

       ! increase counter: this is the index in HcoState%Spc!
       cnt = cnt + 1

       ! Set species name and GEOS-Chem tracer ID 
       IDX = ModelSpecIDs(MatchIDx(I))
       HcoState%Spc(cnt)%SpcName  = HcoSpecNames(I) 
       HcoState%Spc(cnt)%ModID    = IDX

       ! Molecular weights of species & emitted species.
       HcoState%Spc(cnt)%MW_g   = ModelSpecMW(IDX)
       HcoState%Spc(cnt)%EmMW_g = ModelSpecEmMW(IDX)

       ! Emitted molecules per molecule of species.
       HcoState%Spc(cnt)%MolecRatio = ModelSpecMolecRatio(IDX)

       ! Set Henry coefficients
       HcoState%Spc(cnt)%HenryK0  = ModelSpecK0(IDX)
       HcoState%Spc(cnt)%HenryCR  = ModelSpecCR(IDX)
       HcoState%Spc(cnt)%HenryPKA = ModelSpecPKA(IDX)

       ! Set pointer to trac_tend
       ! NOTE: As long as the HEMCO grid is equal to the GEOS-Chem
       ! grid, the HEMCO emission arrays can directly point to the 
       ! tracer tendency arrays (trac_tend) in the GEOS-Chem chemistry
       ! state.
       IF ( ( HcoState%NX == IIPAR ) .AND.      &
            ( HcoState%NY == JJPAR ) .AND.      &
            ( HcoState%NZ == LLPAR )      ) THEN 
          HcoState%Spc(cnt)%Emis%Val => State_Chm%Trac_Tend(:,:,:,IDX)

          ! Also check for drydep array
          HcoState%Spc(cnt)%Depv%Val => State_Chm%DepSav(:,:,IDX)
       ENDIF

       ! Logfile I/O
       CALL HCO_SPEC2LOG( am_I_Root, HcoState, Cnt )

!=============================================================================
       ! KLUDGE for SESQ
       IF ( TRIM(HcoState%Spc(cnt)%SpcName) == 'LIMO' ) THEN 

          cnt = cnt + 1
          HcoState%Spc(cnt)%ModID    = -1
          HcoState%Spc(cnt)%SpcName  = 'SESQ'
          HcoState%Spc(cnt)%Emis%Val => BIOG_SESQ 
          HcoState%Spc(cnt)%MW_g       = 150.0d0 
          HcoState%Spc(cnt)%EmMW_g     = 150.0d0
          HcoState%Spc(cnt)%MolecRatio = 1.0d0
          HcoState%Spc(cnt)%HenryK0    = 0.0d0
          HcoState%Spc(cnt)%HenryCR    = 0.0d0
          HcoState%Spc(cnt)%HenryPKA   = 0.0d0

          ! Logfile I/O
          CALL HCO_SPEC2LOG( am_I_Root, HcoState, Cnt )
       ENDIF
!=============================================================================

    ENDDO !I
    CALL HCO_MSG(SEP1='-')

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Register_Species
!EOC
END MODULE HCOI_GC_MAIN_MOD
#endif
