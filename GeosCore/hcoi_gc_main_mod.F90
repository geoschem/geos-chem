# if !defined(ESMF_)
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
! within GEOS-Chem. Typically, these routines are called from main.F.
!\\
!\\
! Notes:
! \begin{itemize}
! \item HEMCO is used to calculate all emission fields. The emission tendencies
!  are passed to GEOS-Chem through the Trac\_Tend array of State\_Chm.
! \item Carbon aerosols are treated in HEMCO as black carbon (BC) and organic 
!  carbon (OC). The speciation into hydrophobic and hydrophilic carbon is 
!  done when passing HEMCO emissions to GEOS-Chem (MAP\_HCO2GC). Speciation
!  factors are defined below.
! \item Dust aerosol emissions become directly added to the Tracers array 
!  instead of Trac\_Tend. This is to avoid unrealistic vertical mixing of dust
!  particles.
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

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_GC_Init
  PUBLIC  :: HCOI_GC_Run
  PUBLIC  :: HCOI_GC_Final
  PUBLIC  :: GetHcoState
  PUBLIC  :: GetHcoDiagn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Map_HCO2GC
  PRIVATE :: Regrid_Emis2Sim
  PRIVATE :: Set_Current_Time
  PRIVATE :: ExtState_SetPointers
  PRIVATE :: ExtState_UpdtPointers
  PRIVATE :: ModelSpec_Allocate
  PRIVATE :: Model_SetSpecies
  PRIVATE :: Set_Grid
  PRIVATE :: Get_nHcoSpc
  PRIVATE :: Register_Species
!
! !REMARKS:
!  This module is ignored if you are using HEMCO in an ESMF environment.
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version. 
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  30 Jul 2014 - C. Keller   - Added GetHcoState 
!  01 Aug 2014 - C. Keller   - Now use only OC and BC within HEMCO. 
!  20 Aug 2014 - M. Sulprizio- Modify for POPs simulation
!  21 Aug 2014 - R. Yantosca - Added routine EmissRnPbBe; cosmetic changes
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

  ! Pointers used during initialization (for species matching)
  INTEGER                       :: nHcoSpec
  CHARACTER(LEN= 31),   POINTER :: HcoSpecNames       (:) => NULL()
  INTEGER                       :: nModelSpec
  CHARACTER(LEN= 31),   POINTER :: ModelSpecNames     (:) => NULL()
  INTEGER,              POINTER :: ModelSpecIDs       (:) => NULL()
  REAL(hp),             POINTER :: ModelSpecMW        (:) => NULL()
  REAL(hp),             POINTER :: ModelSpecEmMW      (:) => NULL()
  REAL(hp),             POINTER :: ModelSpecMolecRatio(:) => NULL()
  REAL(hp),             POINTER :: ModelSpecK0        (:) => NULL()
  REAL(hp),             POINTER :: ModelSpecCR        (:) => NULL()
  REAL(hp),             POINTER :: ModelSpecPKA       (:) => NULL()
  INTEGER,              POINTER :: MatchIDx           (:) => NULL()

  !--------------------------
  ! %%% Arrays %%%
  !--------------------------

  ! Internal met fields (will be used by some extensions)
  INTEGER,               TARGET :: HCO_PBL_MAX            ! level
  REAL(hp), ALLOCATABLE, TARGET :: HCO_FRAC_OF_PBL(:,:,:) ! unitless
  REAL(hp), ALLOCATABLE, TARGET :: HCO_PCENTER(:,:,:)     ! Pa
  REAL(hp), ALLOCATABLE, TARGET :: HCO_PEDGE  (:,:,:)     ! Pa
  REAL(hp), ALLOCATABLE, TARGET :: HCO_SZAFACT(:,:)       ! -

  ! Arrays to store J-values (used by Paranox extension)
  REAL(hp), ALLOCATABLE, TARGET :: JNO2(:,:)
  REAL(hp), ALLOCATABLE, TARGET :: JO1D(:,:)

  ! Sigma coordinate (temporary)
  REAL(hp), ALLOCATABLE, TARGET :: ZSIGMA(:,:,:)
!
! !DEFINED PARAMETERS:
!
  ! Hydrophilic and hydrophobic fraction of black carbon
  REAL(dp),           PARAMETER :: BC2BCPI = 0.2_dp  ! hydrophilic
  REAL(dp),           PARAMETER :: BC2BCPO = 0.8_dp  ! hydrophobic

  ! Hydrophilic and hydrophobic fraction of organic carbon
  REAL(dp),           PARAMETER :: OC2OCPI = 0.5_dp  ! hydrophilic
  REAL(dp),           PARAMETER :: OC2OCPO = 0.5_dp  ! hydrophobic

  ! Logical switch that determines if species arrays in HEMCO state
  ! can be pointing to corresponding arrays in State_Chm.
  LOGICAL                       :: UsePtrs2GC = .TRUE.

  ! Temporary toggle for diagnostics
  LOGICAL,            PARAMETER :: DoDiagn = .TRUE.

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
    USE HCO_LogFile_Mod,    ONLY : HCO_SPEC2LOG
    USE HCOI_GC_Diagn_Mod,  ONLY : HCOI_GC_DIagn_Init
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: nHcoSpc, HMRC
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
    CALL Config_ReadFile( am_I_Root, Input_Opt%HcoConfigFile, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Config_ReadFile', LOC )

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LOGFILE_OPEN( RC=HMRC ) 
       IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Open Logfile', LOC )
    ELSE
       ! If this is not the root CPU, always disable verbose mode.
       CALL HCO_VERBOSE_SET ( .FALSE. )
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
    CALL Get_nHcoSpc( Input_Opt, nHcoSpc, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Get_nHcoSpc', LOC )

    !-----------------------------------------------------------------
    ! Now that number of HEMCO species are known, initialize HEMCO
    ! state object.
    CALL HcoState_Init( am_I_Root, HcoState, nHcoSpc, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HcoState_Init', LOC )

    !-----------------------------------------------------------------
    ! Set grid. This has to be done before register the species.
    CALL Set_Grid( am_I_Root, State_Met, HcoState, RC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Set_Grid', LOC )

    !-----------------------------------------------------------------
    ! Register species. This will define all species properties
    ! (names, molecular weights, etc.) of the HEMCO species.
    ! If the HEMCO grid is the same as the GEOS-Chem grid, each HEMCO
    ! species is connected to the corresponding Trac_Tend array of the
    ! GEOS-Chem chemistry state object, so that HEMCO directly writes
    ! emissions into these arrays.
    CALL Register_Species( am_I_Root, Input_Opt, State_Chm, HcoState, RC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'Register_Species', LOC )

    !=================================================================
    ! Set misc. parameter
    !=================================================================

    ! Emission, chemistry and dynamics timestep in seconds
    HcoState%TS_EMIS = GET_TS_EMIS() * 60.0
    HcoState%TS_CHEM = GET_TS_CHEM() * 60.0
    HcoState%TS_DYN  = GET_TS_DYN()  * 60.0

    ! This is not an ESMF simulation
    HcoState%isESMF = .FALSE.  

    ! HEMCO configuration file
    HcoState%ConfigFile = Input_Opt%HcoConfigFile

    !=================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and 
    ! the HEMCO configuration file is removed from buffer in this
    ! step. Also initializes the HEMCO clock
    !=================================================================
    CALL HCO_Init( am_I_Root, HcoState, HMRC )
    IF( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC )

    !=================================================================
    ! Initialize the ExtState object
    !
    ! NOTE: This used to be done in routine HCOX_INIT.  Moved this
    ! call here so that we pass some additional quantities to HEMCO
    ! via scalar or logical fields of ExtState. (bmy, 9/29/14)
    !=================================================================

    ! Initialize extension object
    CALL ExtStateInit( ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Define fields of ExtState to be passed to HEMCO
    ExtState%N_DUST_BINS    =  Input_Opt%N_DUST_BINS   ! # of dust bins

#if defined( TOMAS )
    ExtState%IBINS          =  IBINS                   ! # of TOMAS size bins
    ExtState%Xk             => Xk                      ! Size bin edges

# if defined( TOMAS40 )
    ExtState%ACTMODEBINS    =  10                      ! # of activation
# elif defined( TOMAS15 )                              ! mode bins for TOMAS
    ExtState%ACTMODEBINS    =  3
# else 
    ExtState%ACTMODEBINS    =  0
# endif
#endif

    !=================================================================
    ! Initialize all HEMCO extensions.  
    ! Also selects the required met fields used by each extension.
    !=================================================================
    CALL HCOX_Init( am_I_Root, HcoState, ExtState, HMRC, NoExtStateInit=.TRUE. )
    IF( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC )

    !-----------------------------------------------------------------
    ! Update logical switches in Input_Opt 
    !-----------------------------------------------------------------
    Input_Opt%LSOILNOX      = ExtState%SoilNOx

    !-----------------------------------------------------------------
    ! Set constants for POPs simulation
    !-----------------------------------------------------------------
    IF ( ExtState%GC_POPs ) THEN
       ExtState%POP_DEL_H   = Input_Opt%POP_DEL_H
       ExtState%POP_KOA     = Input_Opt%POP_KOA
       ExtState%POP_KBC     = Input_Opt%POP_KBC
    ENDIF

    !-----------------------------------------------------------------
    ! Set pointers to met fields.
    ! Extensions typically depend on environmental dependent met. 
    ! variables such as wind speed, surface temp., etc. Pointers 
    ! to these (2D or 3D) fields are defined in the extension object. 
    ! Here, we need to make sure that these pointers are properly 
    ! connected.
    !-----------------------------------------------------------------
    CALL ExtState_SetPointers( State_Met, State_Chm, RC )
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

    ! Deallocate local variables
    IF ( ASSOCIATED( ModelSpecNames      ) ) DEALLOCATE( ModelSpecNames      )
    IF ( ASSOCIATED( ModelSpecIDs        ) ) DEALLOCATE( ModelSpecIDs        )
    IF ( ASSOCIATED( ModelSpecMW         ) ) DEALLOCATE( ModelSpecMW         )
    IF ( ASSOCIATED( ModelSpecEmMW       ) ) DEALLOCATE( ModelSpecEmMW       )
    IF ( ASSOCIATED( ModelSpecMolecRatio ) ) DEALLOCATE( ModelSpecMolecRatio )
    IF ( ASSOCIATED( ModelSpecK0         ) ) DEALLOCATE( ModelSpecK0         )
    IF ( ASSOCIATED( ModelSpecCR         ) ) DEALLOCATE( ModelSpecCR         )
    IF ( ASSOCIATED( ModelSpecPKA        ) ) DEALLOCATE( ModelSpecPKA        )
    IF ( ASSOCIATED( matchIDx            ) ) DEALLOCATE( matchIDx            )
    IF ( ASSOCIATED( HcoSpecNames        ) ) DEALLOCATE( HcoSpecNames        )

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
    USE HCO_DIAGN_MOD,         ONLY : HCO_DIAGN_AUTOUPDATE
    USE HCO_FLUXARR_MOD,       ONLY : HCO_FluxarrReset 
    USE HCO_DRIVER_MOD,        ONLY : HCO_RUN
    USE HCOX_DRIVER_MOD,       ONLY : HCOX_RUN
    USE HCOIO_DIAGN_MOD,       ONLY : HCOIO_DIAGN_WRITEOUT

    ! For soilnox
    USE GET_NDEP_MOD,          ONLY : RESET_DEP_N
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
! !REMARKS:
!  Modifi
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  22 Aug 2014 - R. Yantosca - Now pass State_Met to MAP_HCO2GC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: HMRC 
    CHARACTER(LEN=255), PARAMETER  :: LOC='HCOI_GC_RUN (hcoi_gc_main_mod.F90)'

    !=======================================================================
    ! HCOI_GC_RUN begins here!
    !=======================================================================

    ! Set return code flag to HCO success. This value should be
    ! preserved throughout all HCO calls, otherwise an error
    ! will be returned!
    HMRC = HCO_SUCCESS

    !=======================================================================
    ! Set HcoClock 
    !=======================================================================
    CALL SET_CURRENT_TIME ( am_I_Root, HcoState, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SET_CURRENT_TIME', LOC )

    !=======================================================================
    ! Output diagnostics 
    !=======================================================================
    IF ( DoDiagn ) THEN
    CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root, HcoState, .FALSE., HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_WRITEOUT', LOC )
    ENDIF

    ! ======================================================================
    ! Reset all emission and deposition values
    ! ======================================================================
    CALL HCO_FluxarrReset ( HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('ResetArrays', LOC )
       RETURN 
    ENDIF
 
    !=======================================================================
    ! Set HCO options and define all arrays needed by core module 
    ! and the extensions 
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
    ! Emissions will be written into the corresponding flux arrays 
    ! in HcoState. 
    !=======================================================================
    CALL HCO_RUN ( am_I_Root, HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('HCO_RUN', LOC )
       RETURN 
    ENDIF

    !=======================================================================
    ! Eventually update variables in ExtState 
    !=======================================================================
    CALL ExtState_UpdtPointers ( State_Met, State_Chm, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('ExtState_UpdtPointers', LOC )
       RETURN 
    ENDIF

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
    ! Update diagnostics 
    !=======================================================================
    IF ( DoDiagn ) THEN
       CALL HCO_DIAGN_AUTOUPDATE ( am_I_Root, HcoState, HMRC )
       IF( HMRC /= HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_UPDATE', LOC )
    ENDIF

    !=======================================================================
    ! Translate emissions array from HCO state onto GC arrays
    !=======================================================================
    CALL MAP_HCO2GC( HcoState, Input_Opt, State_Met, State_Chm, RC )

    !=======================================================================
    ! Reset deposition arrays
    ! TODO: Do somewhere else? e.g. in drydep/wetdep routines?
    !=======================================================================
    CALL RESET_DEP_N()

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
    USE HCOIO_Diagn_Mod,     ONLY : HCOIO_Diagn_WriteOut
    USE HCOX_Driver_Mod,     ONLY : HCOX_Final
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)              :: am_I_Root
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

    ! File prefix for restart file
    CHARACTER(LEN= 31), PARAMETER :: RST = 'HEMCO_restart'

    !=================================================================
    ! HCOI_GC_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_GC_Final (hcoi_gc_main_mod.F90)'

    ! Set HcoClock to current time. This is to make sure that the 
    ! diagnostics are properly written.
    CALL SET_CURRENT_TIME ( am_I_Root, HcoState, HMRC )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SET_CURRENT_TIME', LOC )

    ! Write out 'standard' diagnostics. Use previous time.
    CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root, HcoState, .FALSE., HMRC, &
                                UsePrevTime=.TRUE. )
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCOI_DIAGN_FINAL A', LOC )
 
    ! Also write all other diagnostics into restart file. Use current time.
    CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root, HcoState, .TRUE., HMRC, &
                                UsePrevTime=.FALSE., PREFIX=RST ) 
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'HCOI_DIAGN_FINAL B', LOC )
 
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
    IF ( UsePtrs2GC ) THEN 
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
    IF ( ALLOCATED  ( ZSIGMA          ) ) DEALLOCATE ( ZSIGMA          )
    IF ( ALLOCATED  ( HCO_FRAC_OF_PBL ) ) DEALLOCATE ( HCO_FRAC_OF_PBL )
    IF ( ALLOCATED  ( HCO_PEDGE       ) ) DEALLOCATE ( HCO_PEDGE       )
    IF ( ALLOCATED  ( HCO_PCENTER     ) ) DEALLOCATE ( HCO_PCENTER     )
    IF ( ALLOCATED  ( HCO_SZAFACT     ) ) DEALLOCATE ( HCO_SZAFACT     )
    IF ( ALLOCATED  ( JNO2            ) ) DEALLOCATE ( JNO2            )
    IF ( ALLOCATED  ( JO1D            ) ) DEALLOCATE ( JO1D            )

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
  SUBROUTINE Map_Hco2Gc( HcoState, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,          ONLY : IIPAR,  JJPAR,  LLPAR
    USE Error_Mod,             ONLY : ERROR_STOP
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE Tracerid_Mod 
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState    ! HCO state
    TYPE(OptInput),  INTENT(IN   )  :: Input_Opt   ! Input Options object
    TYPE(MetState),  INTENT(IN   )  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(INOUT)  :: State_Chm   ! Chemistry State object
    INTEGER,         INTENT(INOUT)  :: RC          ! Failure?
!
! !REMARKS:
!  For a detailed discussion of how emissions get added from HEMCO into
!  GEOS-Chem (and how emissions are distributed throughout the boundary
!  layer), please see this wiki page:
!  http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
!  
! !REVISION HISTORY:
!  01 May 2012 - C. Keller   - Initial Version
!  20 Aug 2013 - C. Keller   - Now pass from HEMOC state to chemistry state
!  14 Sep 2013 - C. Keller   - Now keep in units of kg/m2/s.
!  22 Aug 2014 - R. Yantosca - Now get surface area from State_Met%AREA_M2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                  :: N, nSpc, trcID
    CHARACTER(LEN=255)       :: MSG, LOC

    !=======================================================================
    ! MAP_HCO2GC begins here
    !=======================================================================

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

       ! testing only
       IF ( ASSOCIATED(HcoState%Spc(N)%Conc%Val) ) THEN
          write(*,*) 'HEMCO concentrations found for GC tracer ', trcID
          write(*,*) 'total conc: ', SUM(HcoState%Spc(N)%Conc%Val)
       ENDIF

       ! Skip if no emissions defined
       IF ( .NOT. ASSOCIATED(HcoState%Spc(N)%Emis%Val) ) CYCLE

       ! If simulation grid and emission grid are equal, the 
       ! HEMCO state emission array already points to
       ! State_Chm%Trac_Tend and there is nothing to do here. 
       IF ( ( HcoState%NX == IIPAR ) .AND.       &
            ( HcoState%NY == JJPAR ) .AND.       &
            ( HcoState%NZ == LLPAR )       ) THEN

          ! ... nothing to do here        
 
       ! For different grids, regrid emissions onto simulation
       ! grid
       ! TODO: needs testing
       ELSE
          ! Do regridding
          CALL REGRID_EMIS2SIM( HcoState, N, State_Chm, trcID, Input_Opt ) 
       ENDIF

       !----------------------------------------------------------------------
       ! HEMCO holds total OC and BC. Those have GEOS-Chem species IDs of 
       ! OCPI and BCPI (see Model_SetSpecies). Split here into hydrophobic
       ! and hydrophilic fractions
       !----------------------------------------------------------------------
       IF ( trcID == IDTBCPI ) THEN
          State_Chm%Trac_Tend(:,:,:,IDTBCPO) = &
             State_Chm%Trac_Tend(:,:,:,IDTBCPI) * BC2BCPO
          State_Chm%Trac_Tend(:,:,:,IDTBCPI) = &
             State_Chm%Trac_Tend(:,:,:,IDTBCPI) * BC2BCPI

       ELSEIF ( trcID == IDTOCPI ) THEN
          State_Chm%Trac_Tend(:,:,:,IDTOCPO) = &
             State_Chm%Trac_Tend(:,:,:,IDTOCPI) * OC2OCPO
          State_Chm%Trac_Tend(:,:,:,IDTOCPI) = &
             State_Chm%Trac_Tend(:,:,:,IDTOCPI) * OC2OCPI
       ENDIF

       !----------------------------------------------------------------------
       ! Dust emissions shall be directly added to the tracer array instead 
       ! of Trac_Tend. This is to avoid unrealistic vertical mixing of dust
       ! particles. 
       !----------------------------------------------------------------------
       IF ( trcID == IDTDST1 .OR. trcID == IDTDST2 .OR. &
            trcID == IDTDST3 .OR. trcID == IDTDST4       ) THEN
          State_Chm%Tracers(:,:,1,trcID) = State_Chm%Tracers  (:,:,1,trcID) &
                                         + State_Chm%Trac_Tend(:,:,1,trcID) &
                                         * State_Met%AREA_M2(:,:,1)         &
                                         * HcoState%TS_EMIS
          State_Chm%Trac_Tend(:,:,:,trcID) = 0.0d0
       ENDIF

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
! !IROUTINE: ExtState_SetPointers
!
! !DESCRIPTION: SUBROUTINE ExtState\_SetPointers sets the extension object data
! pointers. 
!\\
! Note that for now, this explicitly assumes that the HEMCO emissions grid is 
! the same as the GEOS-Chem simulation grid. To support other grids, the met 
! data has to be regridded explicitly at every time step!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_SetPointers( State_Met, State_Chm, RC ) 
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
!  23 Oct 2012 - C. Keller   - Initial Version
!  20 Aug 2014 - M. Sulprizio- Add PBL_MAX and FRAC_OF_PBL for POPs simulation
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
    ! ExtState_SetPointers begins here
    !=================================================================

    ! Init
    RC = GIGC_SUCCESS

    LOC = 'ExtState_SetPointers (hcoi_gc_main_mod.F90)'

    ! ----------------------------------------------------------------
    ! HCO_PEDGE, HCO_PCENTER and HCO_SZAFACT aren't defined as 3D 
    ! arrays in Met_State. Hence need to construct here so that we 
    ! can point to them.
    !
    ! Now include HCO_FRAC_OF_PBL and HCO_PBL_MAX for POPs specialty
    ! simulation (mps, 8/20/14)
    ! ----------------------------------------------------------------
    IF ( ExtState%PEDGE%DoUse ) THEN 
       ALLOCATE(HCO_PEDGE  (IIPAR,JJPAR,LLPAR+1),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_PEDGE', LOC )
       HCO_PEDGE = 0d0
       ExtState%PEDGE%Arr%Val   => HCO_PEDGE
    ENDIF

    IF ( ExtState%PCENTER%DoUse ) THEN 
       ALLOCATE(HCO_PCENTER(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_PCENTER', LOC )
       HCO_PCENTER = 0d0
       ExtState%PCENTER%Arr%Val => HCO_PCENTER
    ENDIF

    IF ( ExtState%SZAFACT%DoUse ) THEN 
       ALLOCATE(HCO_SZAFACT(IIPAR,JJPAR      ),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_SZAFACT', LOC )
       HCO_SZAFACT = 0d0
       ExtState%SZAFACT%Arr%Val => HCO_SZAFACT
    ENDIF

    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN 
       ALLOCATE(HCO_FRAC_OF_PBL(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_FRAC_OF_PBL', LOC )
       HCO_FRAC_OF_PBL = 0d0
       ExtState%FRAC_OF_PBL%Arr%Val => HCO_FRAC_OF_PBL
    ENDIF

    HCO_PBL_MAX = 0d0
    ExtState%PBL_MAX => HCO_PBL_MAX

    ! ----------------------------------------------------------------
    ! The J-Values for NO2 and O3 are not defined in Met_State. We
    ! need to compute them separately. 
    ! ----------------------------------------------------------------
    IF ( ExtState%JNO2%DoUse ) THEN 
       ALLOCATE( JNO2(IIPAR,JJPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JNO2', LOC )
       JNO2 = 0d0
       ExtState%JNO2%Arr%Val => JNO2 
    ENDIF

    IF ( ExtState%JO1D%DoUse ) THEN 
       ALLOCATE( JO1D(IIPAR,JJPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JO1D', LOC )
       JO1D = 0d0
       ExtState%JO1D%Arr%Val => JO1D 
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
    IF ( ExtState%SNODP%DoUse ) THEN
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
       ExtState%SNODP%Arr%Val => State_Met%SNODP
#else
       ExtState%SNODP%Arr%Val => State_Met%SNOW
#endif
    ENDIF
    IF ( ExtState%FRLAND%DoUse ) THEN
       ExtState%FRLAND%Arr%Val => State_Met%FRLAND
    ENDIF
    IF ( ExtState%FROCEAN%DoUse ) THEN
       ExtState%FROCEAN%Arr%Val => State_Met%FROCEAN
    ENDIF
    IF ( ExtState%FRLAKE%DoUse ) THEN
       ExtState%FRLAKE%Arr%Val => State_Met%FRLAKE
    ENDIF
    IF ( ExtState%FRLANDIC%DoUse ) THEN
       ExtState%FRLANDIC%Arr%Val => State_Met%FRLANDIC
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

  END SUBROUTINE ExtState_SetPointers
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_UpdtPointers
!
! !DESCRIPTION: SUBROUTINE ExtState\_UpdtPointers updates the extension 
! object data pointers. Updates are only required for the shadow arrays
! defined in this module, such as pressure edges, J-values, etc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_UpdtPointers( State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Met_Mod,    ONLY : MetState
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR

    USE PRESSURE_MOD,          ONLY : GET_PEDGE, GET_PCENTER
    USE GLOBAL_OH_MOD,         ONLY : GET_SZAFACT
    USE PBL_MIX_MOD,           ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L

    USE FAST_JX_MOD,           ONLY : FJXFUNC
    USE COMODE_LOOP_MOD,       ONLY : NCS, JPHOTRAT, NRATES
    USE COMODE_LOOP_MOD,       ONLY : NAMEGAS, IRM
!
! !INPUT PARAMETERS:
!
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
       LOC = 'ExtState_UpdtPointers (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! ExtState_UpdtPointers begins here
    !=================================================================

    ! Init
    RC = GIGC_SUCCESS

    ! If necessary, calculate internal met fields
    IF ( ExtState%PEDGE%DoUse   .OR. ExtState%PCENTER%DoUse .OR. &
         ExtState%SZAFACT%DoUse .OR. ExtState%JNO2%DoUse    .OR. &
         ExtState%JO1D%DoUse    .OR. ExtState%FRAC_OF_PBL%DoUse ) THEN

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, K, NK, KMAX, SPECNAME )
       ! Loop over all grid boxes
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR

          ! pressure edges [Pa]
          IF ( ExtState%PEDGE%DoUse ) THEN
             HCO_PEDGE(I,J,L) = GET_PEDGE( I,J,L) * 100.0_hp
             IF ( L==LLPAR ) HCO_PEDGE(I,J,L+1) = GET_PEDGE(I,J,L+1) * 100.0_hp
          ENDIF

          ! pressure centers [Pa]
          IF ( ExtState%PCENTER%DoUse ) THEN
             HCO_PCENTER(I,J,L) = GET_PCENTER(I,J,L) * 100.0_hp
          ENDIF

          ! current SZA divided by total daily SZA (2D field only)
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
             IF ( State_Met%SUNCOSmid(I,J) == 0d0 ) THEN
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
    ENDIF

  END SUBROUTINE ExtState_UpdtPointers
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Model_SetSpecies 
!
! !DESCRIPTION: Subroutine Model\_SetSpecies defines information on the 
! GEOS-Chem species. The names of these species will then be compared against 
! the species names found in the HEMCO configuration file, and only matching
! species will be used by HEMCO (with the properties defined here being copied
! to the HEMCO state object).
!\\
!\\
! Typically, the GEOS-Chem species defined here are nothing else than the
! GEOS-Chem tracers. However, species may be manually added to that list (e.g. 
! SESQ) or be specified entirely manually for specialty simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Model_SetSpecies( Input_Opt, nModelSpec, RC )
!
! !USES:
!
    USE GIGC_State_Chm_Mod,    ONLY : Get_Indx
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE HENRY_COEFFS,          ONLY : Get_Henry_Constant
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT) :: Input_Opt  ! Input Options object
!
! !OUPTUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: nModelSpec
    INTEGER,            INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller   - Initial Version
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  04 Sep 2014 - R. Yantosca - Include more specialty sims in IF statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: N,  IDTLIMO, ID_EMIT
    REAL(dp)           :: K0, CR,  pKa
    CHARACTER(LEN= 31) :: ThisName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'Model_SetSpecies (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! Model_SetSpecies begins here
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

       ! # of model species
       nModelSpec = Input_Opt%N_TRACERS
   
       ! Check for SESQ: SESQ is not transported due to its short lifetime,
       ! but emissions are still calculated (in MEGAN). SESQ is only used
       ! in the SOA simulation, i.e. if LIMO is defined. Thus, add one more
       ! species here if LIMO is a model species and calculate SESQ emissions
       ! along with LIMO!
       IDTLIMO = Get_Indx('LIMO', Input_Opt%ID_TRACER, Input_Opt%TRACER_NAME )
       IF ( IDTLIMO > 0 ) THEN
          nModelSpec = nModelSpec + 1
       ENDIF
 
       ! Allocate model species variables
       CALL ModelSpec_Allocate ( nModelSpec, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
 
       ! Assign variables
       DO N = 1, Input_Opt%N_TRACERS 
   
          ! Get species names
          ! ==> Treat BCPI as BC and OCPI as OC. Will be split into
          !     hydrophobic and hydrophilic fraction lateron!
          IF ( TRIM(Input_Opt%TRACER_NAME(N)) == 'BCPI' ) THEN
             ModelSpecNames(N) = 'BC' 
          ELSEIF ( TRIM(Input_Opt%TRACER_NAME(N)) == 'OCPI' ) THEN
             ModelSpecNames(N) = 'OC' 
          ELSE 
             ModelSpecNames(N) = Input_Opt%TRACER_NAME(N)
          ENDIF
          
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
          CALL GET_HENRY_CONSTANT ( TRIM(ModelSpecNames(N)), K0, CR, pKa, RC )
          ModelSpecK0(N)  = K0
          ModelSpecCR(N)  = CR
          ModelSpecPKA(N) = PKA
       ENDDO      
   
       ! Eventually add SESQ
       IF ( IDTLIMO > 0 ) THEN
          N                      = nModelSpec
          ModelSpecIDs(N)        = -1
          ModelSpecNames(N)      = 'SESQ'
          ModelSpecEmMW(N)       = 150.0_hp
          ModelSpecMW(N)         = 150.0_hp
          ModelSpecMolecRatio(N) = 1.0_hp
          ModelSpecK0(N)         = 0.0_hp
          ModelSpecCR(N)         = 0.0_hp
          ModelSpecPKA(N)        = 0.0_hp
       ENDIF

    !-----------------------------------------------------------------
    ! CO2 specialty simulation 
    !-----------------------------------------------------------------
    ELSEIF ( Input_Opt%ITS_A_CO2_SIM ) THEN

       !--------------------------------------------------------------
       ! For the CO2 specialty simulation, define here all tagged 
       ! tracer. This will let HEMCO calculate emissions for each
       ! tagged tracer individually. The emissions will be passed
       ! to the CO2 arrays in co2_mod.F 
       !--------------------------------------------------------------

       ! There are up to 10 tracers
       nModelSpec = 10 
   
       ! Allocate model species variables
       CALL ModelSpec_Allocate ( nModelSpec, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
 
       ! Henry constants are the same for all tracers
       CALL GET_HENRY_CONSTANT ( 'CO2', K0, CR, pKa, RC )

       ! Assign variables
       DO N = 1, nModelSpec
   
          ! Define species names. These are the names that must also be 
          ! used in the HEMCO configuration file!
          SELECT CASE ( N )

             CASE ( 1  ) 
                ThisName = 'CO2ff'
             CASE ( 2  ) 
                ThisName = 'CO2oc'
             CASE ( 3  ) 
                ThisName = 'CO2bal'
             CASE ( 4  ) 
                ThisName = 'CO2bb'
             CASE ( 5  ) 
                ThisName = 'CO2bf'
             CASE ( 6  ) 
                ThisName = 'CO2nte'
             CASE ( 7  ) 
                ThisName = 'CO2se'
             CASE ( 8  ) 
                ThisName = 'CO2av'
             CASE ( 9  ) 
                ThisName = 'CO2ch'
             CASE ( 10 ) 
                ThisName = 'CO2corr'

             CASE DEFAULT
                MSG = 'Only 10 species defined for CO2 simulation!'
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN

          END SELECT

          ! Species name 
          ModelSpecNames(N)      = ThisName 
          
          ! Species ID. Tracer 1 is total CO2, thus add one to get 
          ! emission tracer ID.
          ModelSpecIDs(N)        = N + 1
   
          ! Molecular weights and molecule ratio
          ModelSpecMW(N)         = Input_Opt%Tracer_MW_G(N)
          ModelSpecEmMW(N)       = Input_Opt%Tracer_MW_G(N)
          ModelSpecMolecRatio(N) = 1.0_hp
   
          ! Henry coefficients
          ModelSpecK0(N)         = K0
          ModelSpecCR(N)         = CR
          ModelSpecPKA(N)        = PKA

          ! For CO2 sim, never set pointers to Trac_Tend arrays
          UsePtrs2GC = .FALSE.

       ENDDO

    !-----------------------------------------------------------------
    ! DEFAULT (RETURN W/ ERROR) 
    !-----------------------------------------------------------------
    ELSE
       MSG = 'Invalid simulation type - cannot define model species' 
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE Model_SetSpecies 
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
  SUBROUTINE Set_Grid( am_I_Root, State_Met, HcoState, RC ) 
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
! !IROUTINE: Get_nHcoSpc 
!
! !DESCRIPTION: Subroutine Get\_nHcoSpc returns the number of species that
! shall be used by HEMCO. This number depends on the definitions of the HEMCO
! configuration file (i.e. how many species are defined in there) and the 
! GEOS-Chem species definitions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_nHcoSpc( Input_Opt, nHcoSpec, RC ) 
!
! !USES:
!
    USE HCO_CharTools_Mod,  ONLY : HCO_CharMatch
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
    INTEGER,        INTENT(  OUT)  :: nHcoSpec   ! # of species to be 
                                                 ! used by HEMCO 
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
    INTEGER            :: nConfigSpec, nModelSpec, AS
    CHARACTER(LEN=255) :: LOC = 'Get_nHcoSpc (hcoi_gc_main_mod.F90)'
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! Get_nHcoSpc begins here
    !=================================================================

    ! Extract number of species found in the HEMCO config. file.
    nConfigSpec = Config_GetnSpecies ( )

    ! If at least one species is set in the configuration file, try
    ! to match those species against the GEOS-Chem species.
    IF ( nConfigSpec == 0 ) THEN
       nHcoSpec = 0
    ELSE

       ! Get list of all species names found in the HEMCO config file.
       CALL Config_GetSpecNames( HcoSpecNames, nConfigSpec, RC )
       IF( RC /= HCO_SUCCESS) RETURN 

       ! Extract GC species names and properties. Those will be written
       ! into the module arrays ModelSpec*.
       CALL Model_SetSpecies( Input_Opt, nModelSpec, RC )
       IF ( RC /= HCO_SUCCESS) RETURN
   
       ! This returns the matching indeces of the HEMCO species (HcoSpecNames)
       ! in ModelSpecNames. A value of -1 is returned if no matching species
       ! is found.
       ALLOCATE(MatchIDx(nConfigSpec),STAT=AS)
       IF ( AS/=0 ) THEN 
          CALL HCO_ERROR ('Allocation error matchIDx', RC, THISLOC=LOC )
          RETURN
       ENDIF
       MatchIDx(:) = -1
       CALL HCO_CharMatch( HcoSpecNames,   nConfigSpec,   &
                           ModelSpecNames, nModelSpec,    &
                           MatchIDx,       nHcoSpec        )
    ENDIF

    IF ( nHcoSpec == 0 ) THEN
       MSG = 'No matching species between HEMCO and the model!'
       CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Get_nHcoSpc 
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
  SUBROUTINE Register_Species( am_I_Root, Input_Opt, State_Chm, HcoState, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_LogFile_Mod,    ONLY : HCO_SPEC2LOG
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR

    ! For SOA mechanism
    USE CARBON_MOD,         ONLY : BIOG_SESQ
!
! !INPUT ARGUMENTS:
!
    TYPE(OptInput),     INTENT(IN   )  :: Input_Opt  ! Input Options object
    LOGICAL,            INTENT(IN   )  :: am_I_Root
    TYPE(ChmState),     INTENT(IN   )  :: State_Chm  ! Chem state
!
! !INPUT/OUTPUT ARGUMENTS:
!
    TYPE(Hco_State),    POINTER        :: HcoState   ! HEMCO state
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

    ! Only if # of HEMCO species is not zero
    IF ( HcoState%nSpc > 0 ) THEN

       ! Loop over all possible HEMCO species
       cnt = 0 
       DO I = 1, SIZE(MatchIDx)
   
          ! Skip if this HEMCO species is not used in GEOS-Chem
          IF ( MatchIDx(I) < 0 ) CYCLE
   
          ! increase counter: this is the index in HcoState%Spc!
          cnt = cnt + 1
   
          ! Set species name and GEOS-Chem tracer ID 
          IDX                          = MatchIDx(I)
          HcoState%Spc(cnt)%ModID      = ModelSpecIDs(IDX)
          HcoState%Spc(cnt)%SpcName    = HcoSpecNames(I) 
   
          ! Molecular weights of species & emitted species.
          HcoState%Spc(cnt)%MW_g       = ModelSpecMW(IDX)
          HcoState%Spc(cnt)%EmMW_g     = ModelSpecEmMW(IDX)
   
          ! Emitted molecules per molecule of species.
          HcoState%Spc(cnt)%MolecRatio = ModelSpecMolecRatio(IDX)
   
          ! Set Henry coefficients
          HcoState%Spc(cnt)%HenryK0    = ModelSpecK0(IDX)
          HcoState%Spc(cnt)%HenryCR    = ModelSpecCR(IDX)
          HcoState%Spc(cnt)%HenryPKA   = ModelSpecPKA(IDX)
   
          !----------------------------------------------------------------------
          ! Set pointer to trac_tend
          !
          ! As long as the HEMCO grid is equal to the GEOS-Chem grid, the 
          ! HEMCO emission arrays can directly point to the tracer tendency 
          ! arrays (trac_tend) in the GEOS-Chem chemistry state. The same 
          ! applies to the deposition array.
          ! 
          ! SESQ emissions are not transported and hence don't have a Trac_Tend 
          ! array. Point them to the BIOG_SESQ array of CARBON_MOD instead.
          ! 
          ! For the CO2 simulation, don't link the HEMCO emissions to Trac_Tend
          ! because PBL mixing is not applied. Instead, we will add emissions to
          ! the surface layer in co2_mod.F. To enable PBL mixing for the CO2
          ! simulation, just remove the logical switch below as well as the
          ! corresponding code in co2_mod.F.
          !----------------------------------------------------------------------
          IF ( UsePtrs2GC ) THEN
   
             ! Only if HEMCO and GEOS-Chem are on the same grid...
             IF ( ( HcoState%NX == IIPAR ) .AND.      &
                  ( HcoState%NY == JJPAR ) .AND.      &
                  ( HcoState%NZ == LLPAR )      ) THEN 
      
                IF ( TRIM(HcoState%Spc(cnt)%SpcName) == 'SESQ' ) THEN
                   HcoState%Spc(cnt)%Emis%Val => BIOG_SESQ 
                   HcoState%Spc(cnt)%Depv%Val => NULL()
                ELSE
                   HcoState%Spc(cnt)%Emis%Val => State_Chm%Trac_Tend(:,:,:,IDX)
                   HcoState%Spc(cnt)%Depv%Val => State_Chm%DepSav(:,:,IDX)
                ENDIF
   
             ELSE
                UsePtrs2GC = .FALSE.
             ENDIF ! Same grid
          ENDIF 
   
          ! Write to logfile
          CALL HCO_SPEC2LOG( am_I_Root, HcoState, Cnt )
   
       ENDDO !I
       CALL HCO_MSG(SEP1='-')

    ENDIF 

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Register_Species
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
! !IROUTINE: GetHcoDiagn 
!
! !DESCRIPTION: Subroutine GetHcoDiagn is a convenience wrapper routine to 
! get a HEMCO diagnostics from somewhere within GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetHcoDiagn ( am_I_Root, DiagnName, Force, RC, Ptr2D, Ptr3D )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_DIAGN_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )      :: am_I_Root  ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN   )      :: DiagnName  ! Name of diagnostics
    LOGICAL,          INTENT(IN   )      :: Force      ! Force error if diagn. not found?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)      :: RC         ! Error return code
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         POINTER, OPTIONAL  :: Ptr2D(:,:)   ! Pointer to 2D data
    REAL(hp),         POINTER, OPTIONAL  :: Ptr3D(:,:,:) ! Pointer to 3D data
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
    INTEGER                   :: FLAG, ERR, LevIDx
    TYPE(DiagnCont), POINTER  :: DgnCont  => NULL()

    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'GetHcoDiagn (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! GetHcoDiagn begins here 
    !=======================================================================

    ! Check HEMCO state object
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       CALL ERROR_STOP ( 'HcoState not defined', LOC )
    ENDIF

    ! Get diagnostics by name. Search all diagnostics, i.e. both AutoFill
    ! and manually filled diagnostics. Also include those with a manual
    ! output interval.
    CALL Diagn_Get( am_I_Root,   HcoState, .FALSE., DgnCont,       &
                    FLAG,        ERR,      cName=TRIM(DiagnName),  &
                    AutoFill=-1, InclManual=.TRUE. )     

    ! Error checks
    IF ( ERR /= HCO_SUCCESS ) THEN
       MSG = 'Error in getting diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP ( MSG, LOC )
    ENDIF
    IF ( (FLAG /= HCO_SUCCESS) .AND. Force ) THEN
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
             MSG = 'no data defined: ' // TRIM(DiagnName)
             CALL ERROR_STOP ( MSG, LOC )
          ENDIF 
   
       ! 3D pointer: must point to 3D data
       ELSEIF ( PRESENT(Ptr3D) ) THEN
          IF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
             Ptr3D => DgnCont%Arr3D%Val
          ELSE
             MSG = 'no 3D data defined: ' // TRIM(DiagnName)
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ModelSpec_Allocate 
!
! !DESCRIPTION: Subroutine ModelSpec\_Allocate allocates the model species
! arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ModelSpec_Allocate ( N, RC )
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER, INTENT(IN   ) :: N     ! Array size
    INTEGER, INTENT(INOUT) :: RC    ! Return code 
!
! !REVISION HISTORY:
!  01 Aug 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER            :: AS
    CHARACTER(LEN=255) :: LOC = 'ModelSpec_Allocate (hcoi_gc_main_mod.F90)'

    !=================================================================
    ! ModelSpec_Allocate begins here
    !=================================================================

    ALLOCATE(ModelSpecNames     (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecNames', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecIDs       (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecIDs', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecMW        (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecMW', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecEmMW      (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecEmMW', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecMolecRatio(N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecMolecRatio', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecK0        (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecK0', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecCR        (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecCR', RC, THISLOC=LOC )
    ENDIF

    ALLOCATE(ModelSpecPKA       (N), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR ( 'Allocation error: ModelSpecPKA', RC, THISLOC=LOC )
    ENDIF

  END SUBROUTINE ModelSpec_Allocate
!EOC
END MODULE HCOI_GC_MAIN_MOD
#endif
