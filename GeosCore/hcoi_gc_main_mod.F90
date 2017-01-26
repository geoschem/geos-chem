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
  USE Precision_Mod
  USE HCO_Error_Mod
  USE HCO_Interface_Mod
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
  PUBLIC  :: HCOI_GC_WriteDiagn
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
!  01 Sep 2015 - R. Yantosca - Remove routine SetSpcMw; we now get parameters
!                              for species from the species database object.
!  27 Feb 2016 - C. Keller   - Update to HEMCO v2.0
!  02 May 2016 - R. Yantosca - Now define IDTPOPG as a module variable
!  16 Jun 2016 - J. Sheng    - Add species index retriever
!  20 Jun 2016 - R. Yantosca - Now define species ID's as module variables
!                              so that we can define them in HCOI_GC_INIT
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !--------------------------
  ! %%% Species ID's %%%
  !--------------------------
  INTEGER                       :: id_HNO3
  INTEGER                       :: id_LIMO
  INTEGER                       :: id_NO
  INTEGER                       :: id_NO2
  INTEGER                       :: id_O3
  INTEGER                       :: id_POPG

  !--------------------------
  ! %%% Pointers %%%
  !--------------------------

  !--------------------------
  ! %%% Arrays %%%
  !--------------------------

  ! Internal met fields (will be used by some extensions)
  INTEGER,               TARGET :: HCO_PBL_MAX                      ! level
  REAL(hp), POINTER             :: HCO_FRAC_OF_PBL(:,:,:)
  REAL(hp), POINTER             :: HCO_SZAFACT(:,:)

  ! Arrays to store J-values (used by Paranox extension)
  REAL(hp), POINTER             :: JNO2(:,:)
  REAL(hp), POINTER             :: JOH(:,:)

  ! Sigma coordinate (temporary)
  REAL(hp), POINTER             :: ZSIGMA(:,:,:)

  ! Sum of cosine of the solar zenith angle. Used to impose a
  ! diurnal variability on OH concentrations
  REAL(fp), POINTER             :: SUMCOSZA(:,:)
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
! types and arrays. The HEMCO configuration is read from the HEMCO 
! configuration file (as listed in Input_Opt%HcoConfigFile) and stored in 
! the HEMCO configuration object. The entire HEMCO setup is based upon the
! entries in the HEMCO configuration object. It is possible to explicitly 
! provide a (previously read) HEMCO configuration object via input argument 
! `HcoConfig`. In this case the HEMCO configuration file will not be read
! any more. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Init( am_I_Root, Input_Opt, State_Met, State_Chm, &
                           RC,        HcoConfig ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE TIME_MOD,           ONLY : GET_TS_EMIS, GET_TS_DYN
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#if defined( TOMAS ) 
    USE TOMAS_MOD,          ONLY : IBINS
    USE TOMAS_MOD,          ONLY : Xk
#endif

    ! HEMCO routines 
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE HCO_Config_Mod,     ONLY : Config_ReadFile
    USE HCO_State_Mod,      ONLY : HcoState_Init
    USE HCO_Driver_Mod,     ONLY : HCO_Init
    USE HCOI_GC_Diagn_Mod,  ONLY : HCOI_GC_Diagn_Init
    USE HCOX_Driver_Mod,    ONLY : HCOX_Init
    USE HCOX_State_Mod,     ONLY : ExtStateInit
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )          :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )          :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )          :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)          :: Input_Opt  ! Input opts
    TYPE(ConfigObj),  POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
    INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  07 Jul 2014 - C. Keller   - Now match species and set species properties
!                              via module variables.
!  30 Sep 2014 - R. Yantosca - Now pass fields for aerosol and microphysics
!                              options to extensions via HcoState
!  13 Feb 2015 - C. Keller   - Now read configuration file in two steps.
!  04 Apr 2016 - C. Keller   - Now accept optional input argument HcoConfig.
!  16 Jun 2016 - J. Sheng    - Add species index retriever
!  20 Jun 2016 - R. Yantosca - Now initialize all species ID's here
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                         :: LSTRAT,  FOUND
    INTEGER                         :: nHcoSpc, HMRC

    ! Strings
    CHARACTER(LEN=255)              :: OptName, LOC, MSG, INS

    ! Pointers
    TYPE(ConfigObj), POINTER        :: iHcoConfig => NULL()

    !=================================================================
    ! HCOI_GC_INIT begins here!
    !=================================================================

    ! Error handling 
    LOC  = 'HCOI_GC_Init (GeosCore/hcoi_gc_main_mod.F90)'
    INS  = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    ! Set return code flag to HCO success. This value should be
    ! preserved throughout all HCO calls, otherwise an error
    ! will be returned!
    HMRC = HCO_SUCCESS

    !=================================================================
    ! Define all species ID's here, for use in module routines below
    !=================================================================
    id_HNO3 = Ind_('HNO3')
    id_LIMO = Ind_('LIMO')
    id_NO   = Ind_('NO'  )
    id_NO2  = Ind_('NO2' )
    id_O3   = Ind_('O3'  )
    id_POPG = Ind_('POPG')

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

    ! If HcoConfig is provided
    IF ( PRESENT(HcoConfig) ) iHcoConfig => HcoConfig

    ! Phase 1: read settings and switches
    CALL Config_ReadFile( am_I_Root, iHcoConfig, Input_Opt%HcoConfigFile, 1, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Config_ReadFile', LOC, INS )

    ! Check settings
    CALL CheckSettings( am_I_Root, iHcoConfig, Input_Opt, &
                        State_Met, State_Chm, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'CheckSettings', LOC, INS )

    ! Phase 2: read fields
    CALL Config_ReadFile( am_I_Root, iHcoConfig, Input_Opt%HcoConfigFile, 2, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Config_ReadFile', LOC, INS )

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LOGFILE_OPEN( iHcoConfig%Err, RC=HMRC ) 
       IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Open Logfile', LOC, INS )
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
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, State_Chm,  &
                         HcoState,  nHcoSpc, 1, HMRC ) 
!    CALL Get_nHcoSpc( am_I_Root, Input_Opt, nHcoSpc, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP ( 'SetHcoSpecies-1', LOC, INS )

    !-----------------------------------------------------------------
    ! Now that number of HEMCO species are known, initialize HEMCO
    ! state object. Links the HEMCO configuration file object 
    ! iHcoConfig to HcoState%Config.
    CALL HcoState_Init( am_I_Root, HcoState, iHcoConfig, nHcoSpc, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP ( 'HcoState_Init', LOC )

    !-----------------------------------------------------------------
    ! Register species. This will define all species properties
    ! (names, molecular weights, etc.) of the HEMCO species.
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, State_Chm,  &
                         HcoState,  nHcoSpc, 2, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP ( 'SetHcoSpecies-2', LOC, INS )

    !-----------------------------------------------------------------
    ! Set grid. 
    CALL Set_Grid( am_I_Root, State_Met, HcoState, RC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Set_Grid', LOC, INS )

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
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC, INS )

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
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP( 'HCO_INIT', LOC, INS )

    !-----------------------------------------------------------------
    ! Update and check logical switches in Input_Opt 
    !-----------------------------------------------------------------

    ! Soil NOx
    Input_Opt%LSOILNOX      = ( ExtState%SoilNOx > 0 )

    ! Ginoux dust emissions
    IF ( ExtState%DustGinoux ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          MSG = 'DustGinoux is on in HEMCO but LDUST=F in input.geos'
          CALL ERROR_STOP( MSG, LOC, INS )
       ENDIF
       Input_Opt%LDEAD      = .FALSE.
    ENDIF

    ! DEAD dust emissions
    IF ( ExtState%DustDead > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          MSG = 'DustDead is on in HEMCO but LDUST=F in input.geos'
          CALL ERROR_STOP( MSG, LOC, INS )
       ENDIF
       Input_Opt%LDEAD      = .TRUE.
    ENDIF

    ! Dust alkalinity
    IF ( ExtState%DustAlk ) THEN
       IF ( .not. Input_Opt%LDSTUP ) THEN
          MSG = 'DustAlk is on in HEMCO but LDSTUP=F in input.geos'
          CALL ERROR_STOP( MSG, LOC, INS )
       ENDIF
    ENDIF

    ! Marine organic aerosols
    IF ( ExtState%MarinePOA ) THEN
       IF ( .not. Input_Opt%LMPOA ) THEN
          MSG = 'MarinePOA is on in HEMCO but LMPOA=F in input.geos'
          CALL ERROR_STOP( MSG, LOC, INS )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Set constants for POPs simulation
    !-----------------------------------------------------------------
    IF ( ExtState%GC_POPs ) THEN
       ExtState%POP_DEL_H   = Input_Opt%POP_DEL_H
       ExtState%POP_KOA     = Input_Opt%POP_KOA
       ExtState%POP_KBC     = Input_Opt%POP_KBC
       ExtState%POP_DEL_Hw  = Input_Opt%POP_DEL_Hw
       ExtState%POP_XMW     = Input_Opt%POP_XMW
       ExtState%POP_HSTAR   = Input_Opt%POP_HSTAR
    ENDIF

    !-----------------------------------------------------------------
    ! Initialize ExtState target arrays. 
    ! Extensions typically depend on environmental dependent met. 
    ! variables such as wind speed, surface temp., etc. Pointers 
    ! to these (2D or 3D) fields are defined in the extension object. 
    ! Here, we need to make sure that these pointers are properly 
    ! connected.
    !-----------------------------------------------------------------
    CALL ExtState_InitTargets( am_I_Root, HcoState, ExtState, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Define diagnostics
    !-----------------------------------------------------------------
    IF ( DoDiagn ) THEN 

       ! Set up traditional GEOS-Chem NDxx diagnostics for emissions
       CALL HCOI_GC_DIAGN_INIT                                &
            ( am_I_Root, Input_Opt, HcoState, ExtState, HMRC )

       ! Exit if any of the diagnostics could not be initialized
       IF ( HMRC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'HCOI_GC_DIAGN_INIT', LOC, INS )
       ENDIF
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Eventually remove pointer
    IF ( PRESENT(HcoConfig) ) iHcoConfig => NULL()

    ! Leave w/ success
    RC = GC_SUCCESS

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
    USE ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE GET_NDEP_MOD,          ONLY : RESET_DEP_N ! For soilnox
    USE Input_Opt_Mod,         ONLY : OptInput
    USE State_Met_Mod,         ONLY : MetState
    USE State_Chm_Mod,         ONLY : ChmState

    ! HEMCO routines 
    USE HCO_CLOCK_MOD,         ONLY : HcoClock_Get
    USE HCO_CLOCK_MOD,         ONLY : HcoClock_EmissionsDone
    USE HCO_DIAGN_MOD,         ONLY : HcoDiagn_AutoUpdate
    USE HCO_FLUXARR_MOD,       ONLY : HCO_FluxarrReset 
    USE HCO_DRIVER_MOD,        ONLY : HCO_RUN
    USE HCOX_DRIVER_MOD,       ONLY : HCOX_RUN
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
!  06 Mar 2015 - R. Yantosca - Now create splash page for HEMC
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: HMRC 
    LOGICAL            :: IsEmisTime

    ! Strings
    CHARACTER(LEN=255) :: LOC, INS

    !=======================================================================
    ! HCOI_GC_RUN begins here!
    !=======================================================================

    ! Error handling
    LOC  = 'HCOI_GC_RUN (GeosCore/hcoi_gc_main_mod.F90)'
    INS  = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

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
    IF(HMRC/=HCO_SUCCESS) CALL ERROR_STOP ( 'SetHcoTime', LOC, INS )

    !=======================================================================
    ! See if it's time for emissions. Don't just use the EmisTime flag in
    ! case that we call this routine multiple times. IsEmisTime will only
    ! be true if this is an emission time step AND emissions have not yet
    ! been calculated for that time step.
    !=======================================================================
    CALL HcoClock_Get( am_I_Root, HcoState%Clock, &
                       IsEmisTime=IsEmisTime, RC=HMRC )

    ! ======================================================================
    ! Reset all emission and deposition values. Do this only if it is time
    ! for emissions, i.e. if those values will be refilled.
    ! ======================================================================
    IF ( IsEmisTime .AND. Phase /= 1 ) THEN
       CALL HCO_FluxarrReset ( HcoState, HMRC )
       IF ( HMRC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP('ResetArrays', LOC, INS )
          RETURN 
       ENDIF
    ENDIF

    !=======================================================================
    ! Define pressure edges [Pa] on HEMCO grid.
    !=======================================================================
    CALL GridEdge_Set ( am_I_Root, State_Met, HcoState, HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP('GridEdge_Update', LOC, INS )
       RETURN 
    ENDIF
 
    !=======================================================================
    ! Set HCO options 
    !=======================================================================

    ! Range of species and emission categories.
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
       CALL ERROR_STOP('HCO_RUN', LOC, INS )
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
       CALL ExtState_SetFields( am_I_Root, State_Met, State_Chm, &
                                HcoState,  ExtState,  RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP('ExtState_SetFields', LOC, INS )
          RETURN 
       ENDIF 
   
       CALL ExtState_UpdateFields( am_I_Root, State_Met, State_Chm, &
                                   HcoState,  ExtState,  RC )
       IF ( RC /= GC_SUCCESS ) RETURN
   
       !=======================================================================
       ! Run HCO extensions. Emissions will be added to corresponding
       ! flux arrays in HcoState.
       !=======================================================================
       CALL HCOX_RUN ( am_I_Root, HcoState, ExtState, HMRC )
       IF ( HMRC/= HCO_SUCCESS ) THEN
          CALL ERROR_STOP('HCOX_RUN', LOC, INS )
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
          IF( HMRC /= HCO_SUCCESS) CALL ERROR_STOP ( 'DIAGN_UPDATE', LOC, INS )
       ENDIF
   
       !=======================================================================
       ! Reset the accumulated nitrogen dry and wet deposition to zero. Will
       ! be re-filled in drydep and wetdep.
       !=======================================================================
       CALL RESET_DEP_N()
   
       !=======================================================================
       ! Emissions are now done for this time step
       !=======================================================================
       CALL HcoClock_EmissionsDone( am_I_Root, HcoState%Clock, RC )

    ENDIF  
 
    ! We are now back in GEOS-Chem environment, hence set 
    ! return flag accordingly! 
    RC = GC_SUCCESS

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
    USE CMN_SIZE_Mod,        ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE Error_Mod,           ONLY : Error_Stop
    USE HCO_Driver_Mod,      ONLY : HCO_Final
    USE HCO_Diagn_Mod,       ONLY : DiagnBundle_Cleanup
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
    ! Scalars
    INTEGER             :: HMRC

    ! Strings
    CHARACTER(LEN=255)  :: LOC, INS

    !=================================================================
    ! HCOI_GC_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_GC_Final (GeosCore/hcoi_gc_main_mod.F90)'

    ! Cleanup HCO core. 
    CALL HCO_FINAL( am_I_Root, HcoState, ERROR, HMRC )

    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    CALL HCOX_FINAL ( am_I_Root, HcoState, ExtState, HMRC ) 

    ! Cleanup diagnostics
    CALL DiagnBundle_Cleanup ( HcoState%Diagn )

    ! Cleanup HcoState object 
    CALL HcoState_Final ( HcoState ) 

    ! Module variables
    IF ( ASSOCIATED ( ZSIGMA          ) ) DEALLOCATE( ZSIGMA          )
    IF ( ASSOCIATED ( HCO_FRAC_OF_PBL ) ) DEALLOCATE( HCO_FRAC_OF_PBL )
    IF ( ASSOCIATED ( HCO_SZAFACT     ) ) DEALLOCATE( HCO_SZAFACT     )
    IF ( ASSOCIATED ( JNO2            ) ) DEALLOCATE( JNO2            )
    IF ( ASSOCIATED ( JOH             ) ) DEALLOCATE( JOH             )
    IF ( ASSOCIATED ( SUMCOSZA        ) ) DEALLOCATE( SUMCOSZA        ) 

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
    USE ErrCode_Mod
    USE Error_Mod,           ONLY : Error_Stop
    USE HCOIO_Diagn_Mod,     ONLY : HcoDiagn_Write 
    USE Input_Opt_Mod,       ONLY : OptInput
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
!  06 Jan 2017 - R. Yantosca - Now tell user to check HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: HMRC
    CHARACTER(LEN=255)              :: MSG, LOC, INS

    !=================================================================
    ! HCOI_GC_WriteDiagn begins here!
    !=================================================================

    ! Error handling
    LOC = 'HCOI_GC_WriteDiagn (GeosCore/hcoi_gc_main_mod.F90)'
    INS = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    ! Make sure HEMCO time is in sync 
    CALL SetHcoTime ( am_I_Root, .FALSE., HMRC )
    IF ( HMRC /= HCO_SUCCESS ) CALL ERROR_STOP ( 'SetHcoTime', LOC, INS )

    ! Write diagnostics
    CALL HcoDiagn_Write( am_I_Root, HcoState, RESTART, HMRC )
    IF ( HMRC/=HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Error writing HEMCO diagnostics' 
       CALL ERROR_STOP ( MSG, LOC, INS )
    ENDIF

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE HCOI_GC_WriteDiagn
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
  SUBROUTINE ExtState_InitTargets( am_I_Root, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER        :: HcoState
    TYPE(EXT_STATE),  POINTER        :: ExtState
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
    ! Scalars
    INTEGER            :: AS

    ! Strings
    CHARACTER(LEN=255) :: LOC, INS

    !=================================================================
    ! ExtState_InitTargets begins here
    !=================================================================

    ! Init
    RC  = GC_SUCCESS

    ! Error handling
    LOC = 'ExtState_InitTargets (GeosCore/hcoi_gc_main_mod.F90)'
    INS = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    ! ----------------------------------------------------------------
    ! HCO_SZAFACT is not defined in Met_State.  Hence need to 
    ! define here so that we can point to them.
    !
    ! Now include HCO_FRAC_OF_PBL and HCO_PBL_MAX for POPs specialty
    ! simulation (mps, 8/20/14)
    ! ----------------------------------------------------------------
    IF ( ExtState%SZAFACT%DoUse ) THEN 

       ALLOCATE( SUMCOSZA( IIPAR, JJPAR ), STAT=AS )
       IF ( AS/=0 ) CALL ERROR_STOP ( 'SUMCOSZA', LOC, INS )
       SUMCOSZA = 0.0e0_fp

       ALLOCATE( HCO_SZAFACT( IIPAR, JJPAR      ),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_SZAFACT', LOC, INS )
       HCO_SZAFACT = 0e0_hp
    ENDIF

    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN 
       ALLOCATE(HCO_FRAC_OF_PBL(IIPAR,JJPAR,LLPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'HCO_FRAC_OF_PBL', LOC, INS )
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
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JNO2', LOC, INS )
       JNO2 = 0.0e0_hp
    ENDIF

    IF ( ExtState%JOH%DoUse ) THEN 
       ALLOCATE( JOH(IIPAR,JJPAR),STAT=AS)
       IF ( AS/=0 ) CALL ERROR_STOP ( 'JOH', LOC, INS )
       JOH = 0.0e0_hp
    ENDIF

    ! ----------------------------------------------------------------
    ! Arrays to be copied physically because HEMCO units are not the
    ! same as in GEOS-Chem
    ! ----------------------------------------------------------------

    ! TROPP: GEOS-Chem TROPP is in hPa, while HEMCO uses Pa.
    IF ( ExtState%TROPP%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%TROPP%Arr, IIPAR, JJPAR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%TROPP', LOC, INS )
       ENDIF
    ENDIF

    ! SPHU: GEOS-Chem SPHU is in g/kg, while HEMCO uses kg/kg.
    ! NOTE: HEMCO only uses SPHU surface values.
    IF ( ExtState%SPHU%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SPHU%Arr, IIPAR, JJPAR, 1, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%SPHU', LOC, INS )
       ENDIF
    ENDIF

    ! SUNCOS: HEMCO now calculates SUNCOS values based on its own
    ! subroutine 
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, IIPAR, JJPAR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Allocate ExtState%SUNCOS', LOC, INS )
       ENDIF
    ENDIF

    ! Leave with success
    RC = GC_SUCCESS

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
  SUBROUTINE ExtState_SetFields( am_I_Root, State_Met, State_Chm, &
                                 HcoState,  ExtState,  RC ) 
!
! !USES:
!
    USE HCOX_STATE_MOD,        ONLY : ExtDat_Set

    USE ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE State_Met_Mod,         ONLY : MetState
    USE State_Chm_Mod,         ONLY : ChmState

    ! For SoilNox
    USE Drydep_Mod,            ONLY : DRYCOEFF
    USE Get_Ndep_Mod,          ONLY : DRY_TOTN
    USE Get_Ndep_Mod,          ONLY : WET_TOTN

#if !defined(ESMF_)
    USE MODIS_LAI_MOD,         ONLY : GC_LAI
#endif
    USE MODIS_LAI_MOD,         ONLY : GC_CHLR

#if defined(ESMF_)
    USE HCOI_ESMF_MOD,         ONLY : HCO_SetExtState_ESMF
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
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
    TYPE(EXT_STATE),  POINTER        :: ExtState   ! HEMCO ext. state
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
!  14 Mar 2016 - C. Keller    - Append '_FOR_EMIS' to all HEMCO met field names
!                               to avoid conflict if met-fields are read via
!                               HEMCO.
!  02 May 2016 - R. Yantosca  - Now define IDTPOPG locally
!  30 Jun 2016 - R. Yantosca  - Remove instances of STT.  Now get the advected
!                               species ID from State_Chm%Map_Advect.
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Pointers
    REAL(hp), POINTER  :: Trgt3D(:,:,:) => NULL()

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: HCRC

    ! Strings
    CHARACTER(LEN=255) :: LOC, INS

    !=================================================================
    ! ExtState_SetFields begins here
    !=================================================================

    ! Init
    RC = GC_FAILURE

    ! Error handling
    LOC = 'ExtState_SetFields (GeosCore/hcoi_gc_main_mod.F90)'
    INS = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    ! ----------------------------------------------------------------
    ! Pointers to local module arrays 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SZAFACT, & 
          'SZAFACT_FOR_EMIS',   HCRC, FIRST, HCO_SZAFACT        )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%JNO2, &
            'JNO2_FOR_EMIS',    HCRC, FIRST, JNO2            )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%JOH, &
            'JOH_FOR_EMIS',     HCRC,     FIRST,    JOH     )  
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRAC_OF_PBL, &
     'FRAC_OF_PBL_FOR_EMIS',    HCRC, FIRST, HCO_FRAC_OF_PBL )  
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! 2D fields 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%U10M, &
            'U10M_FOR_EMIS',    HCRC, FIRST, State_Met%U10M )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%V10M, &
            'V10M_FOR_EMIS',    HCRC, FIRST, State_Met%V10M )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%ALBD, &
            'ALBD_FOR_EMIS',    HCRC, FIRST, State_Met%ALBD )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

#if defined ( GCAP )
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%WLI, &
             'WLI_FOR_EMIS',    HCRC, FIRST, State_Met%LWI_GISS  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%WLI, &
             'WLI_FOR_EMIS',    HCRC, FIRST, State_Met%LWI  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%T2M, &
             'T2M_FOR_EMIS',    HCRC, FIRST, State_Met%TS   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%TSKIN, &
           'TSKIN_FOR_EMIS',    HCRC, FIRST, State_Met%TSKIN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%GWETROOT, &
         'GWETROOT_FOR_EMIS',   HCRC, FIRST, State_Met%GWETROOT  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%GWETTOP, &
          'GWETTOP_FOR_EMIS',   HCRC, FIRST, State_Met%GWETTOP  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%USTAR, &
            'USTAR_FOR_EMIS',   HCRC, FIRST, State_Met%USTAR  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%Z0, &
               'Z0_FOR_EMIS',   HCRC, FIRST, State_Met%Z0  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PARDR, &
            'PARDR_FOR_EMIS',   HCRC, FIRST, State_Met%PARDR  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PARDF, &
            'PARDF_FOR_EMIS',   HCRC, FIRST, State_Met%PARDF  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

!    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PSC2, &
!             'PSC2_FOR_EMIS',   HCRC, FIRST, State_Met%PSC2  )
!    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%PSC2_WET, &
            'PSC2_WET_FOR_EMIS', HCRC, FIRST, State_Met%PSC2_WET  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! NOTE: State_Met%RADSWG is net radiation at ground for all
    ! MET fields except GEOS-FP and MERRA2. For those data sets,
    ! net radiation was not processed and so we use incident radiation 
    ! at ground (SWGDN). For simplicity we store radiation as a single 
    ! variable in HEMCO. SWGDN is also available for MERRA but is not 
    ! used in HEMCO to preserve legacy usage of net radiation (ewl, 9/23/15)
#if defined ( GEOS_FP ) || ( MERRA2 )
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%RADSWG, &
           'RADSWG_FOR_EMIS',   HCRC, FIRST, State_Met%SWGDN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%RADSWG, &
           'RADSWG_FOR_EMIS',   HCRC, FIRST, State_Met%RADSWG  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRCLND, &
           'FRCLND_FOR_EMIS',   HCRC, FIRST, State_Met%FRCLND  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CLDFRC, &
            'CLDFRC_FOR_EMIS',   HCRC, FIRST, State_Met%CLDFRC  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

#if defined( GEOS_4 ) 
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNOWHGT, &
          'SNOWHGT_FOR_EMIS',   HCRC, FIRST, State_Met%SNOW     )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNODP, &
            'SNODP_FOR_EMIS',   HCRC, FIRST, State_Met%SNOW   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#elif defined ( GCAP )
   CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNOWHGT, &
         'SNOWHGT_FOR_EMIS',   HCRC, FIRST, State_Met%SNOW     )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNODP, &
            'SNODP_FOR_EMIS',   HCRC, FIRST, State_Met%SNOW   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNICE, &
            'SNICE_FOR_EMIS',   HCRC, FIRST, State_Met%SNICE  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else
    ! SNOWHGT is is mm H2O, which is the same as kg H2O/m2.
    ! This is the unit of SNOMAS.
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNOWHGT, &
          'SNOWHGT_FOR_EMIS',   HCRC, FIRST, State_Met%SNOMAS   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! SNOWDP is in m
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%SNODP, &
            'SNODP_FOR_EMIS',   HCRC, FIRST, State_Met%SNODP  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLAND, &
           'FRLAND_FOR_EMIS',   HCRC, FIRST, State_Met%FRLAND  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FROCEAN, &
          'FROCEAN_FOR_EMIS',   HCRC, FIRST, State_Met%FROCEAN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLAKE, &
           'FRLAKE_FOR_EMIS',   HCRC, FIRST, State_Met%FRLAKE  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%FRLANDIC, &
         'FRLANDIC_FOR_EMIS',   HCRC, FIRST, State_Met%FRLANDIC  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! Use 'offline' MODIS LAI in standard GEOS-Chem
#if defined(ESMF_)
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%LAI, &
              'LAI_FOR_EMIS',   HCRC, FIRST, State_Met%LAI  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#else
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%LAI, &
              'LAI_FOR_EMIS',   HCRC, FIRST, GC_LAI         )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
#endif

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CHLR, &
             'CHLR', HCRC,      FIRST,   GC_CHLR          )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! Convective fractions
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CNV_FRC,  &
          'CNV_FRC_FOR_EMIS',   HCRC, FIRST, State_Met%CNV_FRC, &
          NotFillOk=.TRUE.)
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! 3D fields 
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%CNV_MFC, &
          'CNV_MFC_FOR_EMIS',   HCRC, FIRST, State_Met%CMFMC,  &
          OnLevEdge=.TRUE. )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%TK, &
               'TK_FOR_EMIS',   HCRC, FIRST, State_Met%T   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! Air mass [kg/grid box]
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIR, &
              'AIR_FOR_EMIS',   HCRC, FIRST, State_Met%AD   )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIRVOL, &
           'AIRVOL_FOR_EMIS',   HCRC, FIRST, State_Met%AIRVOL  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    ! Dry air density [kg/m3]
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIRDEN, &
           'AIRDEN', HCRC,      FIRST,    State_Met%AIRDEN  )
    IF ( HCRC /= HCO_SUCCESS ) RETURN
 
    ! ----------------------------------------------------------------
    ! Species concentrations
    ! ----------------------------------------------------------------
    IF ( id_O3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_O3)
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%O3, &
            'HEMCO_O3_FOR_EMIS', HCRC,      FIRST,    Trgt3D )
       IF ( HCRC /= HCO_SUCCESS ) RETURN
       Trgt3D => NULL()
    ENDIF
    IF ( id_NO2 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO2)
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%NO2, &
           'HEMCO_NO2_FOR_EMIS', HCRC,      FIRST,    Trgt3D ) 
       IF ( HCRC /= HCO_SUCCESS ) RETURN
       Trgt3D => NULL()
    ENDIF
    IF ( id_NO > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO)
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%NO, &
            'HEMCO_NO_FOR_EMIS', HCRC,      FIRST,    Trgt3D ) 
       IF ( HCRC /= HCO_SUCCESS ) RETURN
       Trgt3D => NULL()
    ENDIF
    IF ( id_HNO3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_HNO3)
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%HNO3, &
          'HEMCO_HNO3_FOR_EMIS', HCRC,      FIRST,    Trgt3D ) 
       IF ( HCRC /= HCO_SUCCESS ) RETURN
       Trgt3D => NULL()
    ENDIF
    IF ( id_POPG > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_POPG)
       CALL ExtDat_Set( am_I_Root, HcoState, ExtState%POPG, &
          'HEMCO_POPG_FOR_EMIS', HCRC,      FIRST,    Trgt3D ) 
       IF ( HCRC /= HCO_SUCCESS ) RETURN
       Trgt3D => NULL()
    ENDIF

    ! ----------------------------------------------------------------
    ! Deposition parameter
    ! ----------------------------------------------------------------
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%DRY_TOTN, &
         'DRY_TOTN_FOR_EMIS',   HCRC, FIRST, DRY_TOTN            ) 
    IF ( HCRC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%WET_TOTN, &
         'WET_TOTN_FOR_EMIS',   HCRC, FIRST, WET_TOTN            ) 
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
          CALL ERROR_STOP ( 'Error in HCO_SetExtState_ESMF!', LOC, INS )
       ENDIF
    ENDIF
#endif

    ! Not first call any more
    FIRST = .FALSE.

    ! Leave with success
    RC = GC_SUCCESS

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
  SUBROUTINE ExtState_UpdateFields( am_I_Root, State_Met, State_Chm, &
                                    HcoState,  ExtState,  RC          ) 
!
! !USES:
!
    USE CMN_FJX_MOD,           ONLY : ZPJ
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE FAST_JX_MOD,           ONLY : RXN_NO2, RXN_O3_1, RXN_O3_2a
    USE HCO_GeoTools_Mod,      ONLY : HCO_GetSUNCOS
    USE PBL_MIX_MOD,           ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
    USE State_Met_Mod,         ONLY : MetState
    USE State_Chm_Mod,         ONLY : ChmState
#if defined(ESMF_) 
    USE HCOI_ESMF_MOD,         ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chm state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
    TYPE(EXT_STATE),  POINTER        :: ExtState   ! HEMCO ext. state
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller   - Initial Version
!  20 Aug 2014 - M. Sulprizio- Add PBL_MAX and FRAC_OF_PBL for POPs simulation
!  02 Oct 2014 - C. Keller   - PEDGE is now in HcoState%Grid
!  11 Mar 2015 - R. Yantosca - Now call GET_SZAFACT in this module
!  11 Sep 2015 - E. Lundgren - Remove State_Chm from passed args since not used
!  20 Apr 2016 - M. Sulprizio- Change JO1D to JOH to reflect that the array now
!                              holds the effective O3 + hv -> 2OH rates
!  27 Jun 2016 - M. Sulprizio- Obtain photolysis rate directly from ZPJ array
!                              and remove reference to FJXFUNC and obsolete
!                              SMVGEAR variables like NKSO4PHOT, NAMEGAS, etc.
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L

    ! Strings
    CHARACTER(LEN=255) :: LOC, INS

    !=================================================================
    ! ExtState_UpdateFields begins here
    !=================================================================

    ! Init
    RC   = GC_SUCCESS

    ! Error handling
    LOC  = 'ExtState_UpdateFields (GeosCore/hcoi_gc_main_mod.F90)'
    INS  = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

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
          CALL ERROR_STOP( 'Error in HCO_GetSUNCOS', LOC, INS )
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
!$OMP PRIVATE( I, J, L )
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
            (ExtState%JNO2%DoUse .OR. ExtState%JOH%DoUse) ) THEN

          ! Check if sun is up
          IF ( State_Met%SUNCOS(I,J) == 0d0 ) THEN
             IF ( ExtState%JNO2%DoUse ) JNO2 = 0.0_hp
             IF ( ExtState%JOH%DoUse  ) JOH  = 0.0_hp
          ELSE
             IF ( ExtState%JNO2%DoUse ) THEN
                ! RXN_NO2: NO2 + hv --> NO  + O
                JNO2(I,J) = ZPJ(L,RXN_NO2,I,J)
             ENDIF
             IF ( ExtState%JOH%DoUse ) THEN
#if defined( UCX )
                ! RXN_O3_1: O3  + hv --> O2  + O
                JOH(I,J) = ZPJ(L,RXN_O3_1,I,J)
#else
                ! RXN_O3_2a: O3 + hv --> 2OH
                JOH(I,J) = ZPJ(L,RXN_O3_2a,I,J)
#endif
             ENDIF
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
  SUBROUTINE GridEdge_Set ( am_I_Root, State_Met, HcoState, RC )
!
! !USES:
!
    USE ERROR_MOD,             ONLY : ERROR_STOP
    USE State_Met_Mod,         ONLY : MetState
    USE HCOX_STATE_MOD,        ONLY : ExtDat_Set
    USE HCO_GeoTools_MOD,      ONLY : HCO_CalcVertGrid
    USE HCO_GeoTools_MOD,      ONLY : HCO_SetPBLm
    USE CMN_SIZE_MOD,          ONLY : IIPAR, JJPAR, LLPAR
    USE PBL_MIX_MOD,           ONLY : GET_PBL_TOP_M
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root 
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REMARKS:
!  GridEdge_Set defines the HEMCO vertical grid used in GEOS-Chem "classic"
!  simulations.  (GCHP uses its own interface to HEMCO.)
!
! !REVISION HISTORY:
!  08 Oct 2014 - C. Keller   - Initial version
!  28 Sep 2015 - C. Keller   - Now call HCO_CalcVertGrid
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  06 Jun 2016 - R. Yantosca - Now declar PEDGE array for edge pressures [Pa]
!  06 Jun 2016 - R. Yantosca - PSFC now points to PEDGE(:,:,1)
!  06 Jun 2016 - R. Yantosca - Now add error traps
!  26 Oct 2016 - R. Yantosca - Now improve error traps for PBLM
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Integer
    INTEGER           :: I, J

    ! Strings
    CHARACTER(LEN=80) :: MSG, LOC, INS

    ! Pointers
    REAL(hp), POINTER :: PBLM    (:,:  )    ! PBL height           [m ] 
    REAL(hp), POINTER :: BXHEIGHT(:,:,:)    ! Grid box height      [m ]
    REAL(hp), POINTER :: PEDGE   (:,:,:)    ! Pressure @ lvl edges [Pa]
    REAL(hp), POINTER :: PSFC    (:,:  )    ! Surface pressure     [Pa]
    REAL(hp), POINTER :: TK      (:,:,:)    ! Temperature          [K ]
    REAL(hp), POINTER :: ZSFC    (:,:  )    ! Surface geopotential [m ]

    !=======================================================================
    ! GridEdge_Set begins here
    !=======================================================================

    ! Assume success
    RC   =  HCO_SUCCESS

    ! Error handling
    LOC  = 'CheckSettings (GeosCore/hcoi_gc_main_mod.F90)'
    INS  = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    ! Allocate the PEDGE array, which holds level edge pressures [Pa]
    ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
    ! make PEDGE be a pointer and allocate/deallocate it on each call.
    ALLOCATE( PEDGE( IIPAR, JJPAR, LLPAR+1 ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'ERROR allocating the PEDGE pointer-based array!'
       LOC = 'GridEdge_Set (GeosCore/hcoi_gc_main_mod.F90)'
       CALL ERROR_STOP( MSG, LOC, INS )
    ENDIF

    ! Edge and surface pressures [Pa]
    PEDGE    =  State_Met%PEDGE * 100.0_hp  ! Convert hPa -> Pa
    PSFC     => PEDGE(:,:,1)

    ! Point to other fields of State_Met
    ZSFC     => State_Met%PHIS               
    BXHEIGHT => State_Met%BXHEIGHT           
    TK       => State_Met%T                  

    ! Calculate missing quantities
    CALL HCO_CalcVertGrid( am_I_Root, HcoState, PSFC,    &
                           ZSFC,  TK, BXHEIGHT, PEDGE, RC )
    ! Stop with an error if the vertical grid was not computed properly
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'ERROR returning from HCO_CalcVertGrid!'
       LOC = 'GridEdge_Set (GeosCore/hcoi_gc_main_mod.F90)'
       CALL ERROR_STOP( MSG, LOC, INS )
    ENDIF

    ! Set PBL heights
    ALLOCATE( PBLM( IIPAR, JJPAR ), STAT=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'ERROR allocating the PBLM pointer-based array!'
       LOC = 'GridEdge_Set (GeosCore/hcoi_gc_main_mod.F90)'
       CALL ERROR_STOP( MSG, LOC, INS )
    ENDIF

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J )
    DO J=1,JJPAR
    DO I=1,IIPAR
       PBLM(I,J) = GET_PBL_TOP_m(I,J)
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Use the met field PBL field to initialize HEMCO
    CALL HCO_SetPBLm ( am_I_Root, HcoState, FldName='PBL_HEIGHT', &
                       PBLM=PBLM, DefVal=1000.0_hp, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'ERROR returning from HCO_SetPBLm!'
       LOC = 'GridEdge_Set (GeosCore/hcoi_gc_main_mod.F90)'
       CALL ERROR_STOP( MSG, LOC, INS )
    ENDIF

    ! Nullify local pointers
    ZSFC     => NULL()
    BXHEIGHT => NULL()
    TK       => NULL()
    PSFC     => NULL()

    ! Deallocate and nullify PEDGE
    IF ( ASSOCIATED( PEDGE ) ) DEALLOCATE( PEDGE )
    PEDGE    => NULL()

    ! Deallocate the PBLM array
    IF ( ASSOCIATED( PBLM  ) ) DEALLOCATE( PBLM  )
    PBLM     => NULL()

    ! Return w/ success
    RC       = HCO_SUCCESS

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
! are typically just the GEOS-Chem species. Some additional species may be 
! manually added, e.g. SESQ (which is not a species) or individual CO2 species
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
  SUBROUTINE SetHcoSpecies( am_I_Root, Input_Opt, State_Chm, &
                            HcoState,  nSpec,     Phase, RC   )
!
! !USES:
!
    USE HCO_LogFile_Mod,       ONLY : HCO_SPEC2LOG
    USE Input_Opt_Mod,         ONLY : OptInput
    USE Species_Mod,           ONLY : Species
    USE HCO_Types_Mod,         ONLY : ConfigObj
    USE State_Chm_Mod,         ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )   :: am_I_Root  ! Are we on the root CPU
    INTEGER,          INTENT(IN   )   :: Phase      ! 1=Init, 2=Run 
    TYPE(ChmState),   INTENT(IN   )   :: State_Chm  ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)   :: Input_Opt  ! Input Options object
    TYPE(Hco_State),  POINTER         :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT)   :: nSpec      ! # of species for HEMCO
    INTEGER,          INTENT(INOUT)   :: RC         ! Success or failure?
!
! !REMARKS:
!  (1) We now get physical parameters for species from the species database,
!       which is part of the State_Chm object.  
!  (2) In the future, it will be easier to specify non-advected species
!       like SESQ and the CO2 regional species from the species database.
!       The species database flags if a species is advected or not.
!
! !REVISION HISTORY:
!  06 Mar 2015 - C. Keller   - Initial Version
!  01 Sep 2015 - R. Yantosca - Remove reference to GET_HENRY_CONSTANT; we now
!                              get Henry constants from the species database
!  02 May 2016 - R. Yantosca - Now initialize IDTPOPG here
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  20 Jun 2016 - R. Yantosca - All species IDs are now set in the init phase
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: nSpc
    INTEGER                :: N,  L,  M
    REAL(dp)               :: K0, CR, pKa

    ! Strings
    CHARACTER(LEN= 31)     :: ThisName
    CHARACTER(LEN=255)     :: MSG, LOC

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! SetHcoSpecies begins here
    !=================================================================

    ! Error handling
    LOC = 'SetHcoSpecies (GeosCore/hcoi_gc_main_mod.F90)'

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
         Input_Opt%ITS_A_TAGO3_SIM      .or. &
         Input_Opt%ITS_A_TAGCO_SIM    ) THEN

       ! Get number of model species
       nSpc = State_Chm%nAdvect
  
       !%%%%% FOR SOA SIMULATIONS %%%%%
       ! Check for SESQ: SESQ is not transported due to its short lifetime,
       ! but emissions are still calculated (in MEGAN). SESQ is only used
       ! in the SOA simulation, i.e. if LIMO is defined. Thus, add one more
       ! species here if LIMO is a model species and calculate SESQ emissions
       ! along with LIMO!
       IF ( id_LIMO > 0 ) THEN
          nSpc = nSpc + 1
       ENDIF

       !%%%%% FOR THE TAGGED CO SIMULATION %%%%%
       ! Add 3 extra species (ISOP, ACET, MONX) for tagged CO 
       IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
          nSpc = nSpc + 3 
       ENDIF

       ! Assign species variables
       IF ( PHASE == 2 ) THEN

          ! Verbose
          IF ( am_I_Root ) THEN
             MSG = 'Registering HEMCO species:'
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(MSG,*) 'Input species /= expected species: ', nSpec, nSpc 
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          DO N = 1, State_Chm%nAdvect

             ! Get info for this species from the species database
             SpcInfo => State_Chm%SpcData(N)%Info

             ! Model ID and species name 
             HcoState%Spc(N)%ModID      = SpcInfo%ModelID
             HcoState%Spc(N)%SpcName    = TRIM( SpcInfo%Name )
             
             ! Actual molecular weight of species [g/mol]
             HcoState%Spc(N)%MW_g       = SpcInfo%MW_g
             
             ! Emitted molecular weight of species [g/mol].  
             ! Some hydrocarbon species (like ISOP) are emitted and 
             ! transported as a number of equivalent carbon atoms.
             ! For these species, the emitted molecular weight will 
             ! be 12.0 (the weight of 1 carbon atom).
             HcoState%Spc(N)%EmMW_g     = SpcInfo%EmMw_g

             ! Emitted molecules per molecules of species [1].  
             ! For most species, this will be 1.0.  For hydrocarbon 
             ! species (like ISOP) that are emitted and transported
             ! as equivalent carbon atoms, this will be be the number
             ! of moles carbon per mole species.
             HcoState%Spc(N)%MolecRatio = SpcInfo%MolecRatio
   
             ! Set Henry's law coefficients
             HcoState%Spc(N)%HenryK0    = SpcInfo%Henry_K0   ! [M/atm]
             HcoState%Spc(N)%HenryCR    = SpcInfo%Henry_CR   ! [K    ]
             HcoState%Spc(N)%HenryPKA   = SpcInfo%Henry_pKa  ! [1    ]

             ! Write to logfile
             IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, N )

             ! Free pointer memory
             SpcInfo => NULL()
          ENDDO      

          !------------------------------------------------------------------
          ! %%%%% FOR SOA SIMULATIONS %%%%%
          !
          ! Add the non-advected species SESQ in the last species slot
          !------------------------------------------------------------------
          IF ( id_LIMO > 0 ) THEN
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

          !------------------------------------------------------------------
          ! %%%%% FOR THE TAGGED CO SIMULATION %%%%%
          !
          ! Add the non-advected species ISOP, ACET, and MONX
          ! in the last 3 species slots (bmy, ckeller, 6/1/16)
          !------------------------------------------------------------------
          IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
       
             ! Add 3 additional species
             DO L = 1, 3
                
                ! ISOP, ACET, MONX follow the regular tagged CO species
                M = State_Chm%nAdvect + L

                ! Get the species name
                SELECT CASE( L )
                   CASE( 1 ) 
                      ThisName = 'ISOP'
                   CASE( 2 )
                      ThisName = 'ACET'
                   CASE( 3 )
                      ThisName = 'MONX'
                END SELECT

                ! Add physical properties to the HEMCO state
                HcoState%Spc(M)%ModID      = M
                HcoState%Spc(M)%SpcName    = TRIM( ThisName )
                HcoState%Spc(M)%MW_g       = 12.0_hp
                HcoState%Spc(M)%EmMW_g     = 12.0_hp
                HcoState%Spc(M)%MolecRatio = 1.0_hp
                HcoState%Spc(M)%HenryK0    = 0.0_hp
                HcoState%Spc(M)%HenryCR    = 0.0_hp
                HcoState%Spc(M)%HenryPKa   = 0.0_hp

                ! Write to log file
                IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, M )
             ENDDO
          ENDIF

          ! Add line to log-file
          IF ( am_I_Root ) CALL HCO_MSG(HcoState%Config%Err,SEP1='-')
       ENDIF ! Phase = 2   

    !-----------------------------------------------------------------
    ! CO2 specialty simulation 
    ! For the CO2 specialty simulation, define here all tagged 
    ! species. This will let HEMCO calculate emissions for each
    ! tagged species individually. The emissions will be passed
    ! to the CO2 arrays in co2_mod.F 
    !-----------------------------------------------------------------
    ELSEIF ( Input_Opt%ITS_A_CO2_SIM ) THEN

       ! There are up to 11 species
       nSpc = 11 
   
       ! Set species
       IF ( PHASE == 2 ) THEN

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(MSG,*) 'Input species /= expected species: ', nSpec, nSpc 
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Get info about the total CO2 species (i.e. N=1) from the 
          ! species database object.  All tagged CO2 species will
          ! have the same properties as the total CO2 species.
          SpcInfo => State_Chm%SpcData(1)%Info
 
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
                   CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
                   RETURN
   
             END SELECT

             ! Model ID and species name 
             HcoState%Spc(N)%ModID      = N
             HcoState%Spc(N)%SpcName    = TRIM( ThisName )
   
             ! Molecular weights of species & emitted species.
             ! NOTE: Use the species database emMW_g to replicate the 
             ! prior behavior.  The MW's of the tagged species in 
             ! the prior code all have MW_g = 0 and EmMW_g = 0.
             ! Ask Christoph about this. (bmy, 9/1/15)
             HcoState%Spc(N)%MW_g       = SpcInfo%emMW_g           ! [g/mol]
             HcoState%Spc(N)%EmMW_g     = SpcInfo%emMW_g           ! [g/mol]
             HcoState%Spc(N)%MolecRatio = SpcInfo%MolecRatio       ! [1    ]
 
             ! Set Henry coefficients
             HcoState%Spc(N)%HenryK0    = SpcInfo%Henry_K0         ! [M/atm]
             HcoState%Spc(N)%HenryCR    = SpcInfo%Henry_CR         ! [K    ]
             HcoState%Spc(N)%HenryPKA   = SpcInfo%Henry_PKA        ! [1    ]

             ! Write to logfile
             IF ( am_I_Root ) CALL HCO_SPEC2LOG( am_I_Root, HcoState, N )

          ENDDO
          IF ( am_I_Root ) CALL HCO_MSG(HcoState%Config%Err,SEP1='-')

          ! Free pointer
          SpcInfo => NULL()

       ENDIF ! Phase = 2

    !-----------------------------------------------------------------
    ! DEFAULT (RETURN W/ ERROR) 
    !-----------------------------------------------------------------
    ELSE
       MSG = 'Invalid simulation type - cannot define model species' 
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
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
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE GC_GRID_MOD,        ONLY : XMID,  YMID
    USE GC_GRID_MOD,        ONLY : XEDGE, YEDGE, YSIN
    USE GC_GRID_MOD,        ONLY : AREA_M2
    USE HCO_ARR_MOD,        ONLY : HCO_ArrInit
    USE HCO_VERTGRID_MOD,   ONLY : HCO_VertGrid_Define
    USE PRESSURE_MOD,       ONLY : GET_AP, GET_BP
    USE State_Met_Mod,      ONLY : MetState
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
!  13 Sep 2013 - C. Keller   - Initial Version
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  28 Sep 2015 - C. Keller   - Now use HCO_VertGrid_Mod for vertical grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                :: L
    REAL(hp), ALLOCATABLE  :: Ap(:), Bp(:)

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

    ! Initialize the vertical grid. Pass Ap, Bp values from GEOS-Chem grid.
    ALLOCATE(Ap(LLPAR+1),Bp(LLPAR+1))
    DO L = 1, LLPAR+1
       Ap(L) = GET_AP(L) * 100_hp ! hPa to Pa 
       Bp(L) = GET_BP(L)          ! unitless
    ENDDO

    CALL HCO_VertGrid_Define( am_I_Root, HcoState%Config,       &
                              zGrid      = HcoState%Grid%zGrid, &
                              nz         = LLPAR,               &
                              Ap         = Ap,                  & 
                              Bp         = Bp,                  & 
                              RC         = RC                    )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set pointers to grid variables
    HcoState%Grid%XMID%Val       => XMID   (:,:,1)
    HcoState%Grid%YMID%Val       => YMID   (:,:,1)
    HcoState%Grid%XEDGE%Val      => XEDGE  (:,:,1)
    HcoState%Grid%YEDGE%Val      => YEDGE  (:,:,1)
    HcoState%Grid%YSIN%Val       => YSIN   (:,:,1)
    HcoState%Grid%AREA_M2%Val    => AREA_M2(:,:,1)
!    HcoState%Grid%ZSFC%Val       => State_Met%PHIS      ! Surface geopotential height
!    HcoState%Grid%BXHEIGHT_M%Val => State_Met%BXHEIGHT  ! Grid box heights

!    ! Allocate PEDGE. Will be updated every time step!
!    CALL HCO_ArrInit( HcoState%Grid%PEDGE, HcoState%NX, HcoState%NY, HcoState%NZ+1, RC )
!    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE Set_Grid
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
  SUBROUTINE CheckSettings( am_I_Root, HcoConfig, Input_Opt, &
                            State_Met, State_Chm, RC )
!
! !USES:
!
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_ExtList_Mod,    ONLY : GetExtNr,  SetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt, AddExtOpt 
    USE HCO_ExtList_Mod,    ONLY : CoreNr 
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(ConfigObj),  POINTER        :: HcoConfig  ! HEMCO config obj 
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Gfortran will choke unless we use the .eqv. operator to compare LOGICAL 
!  variables for equality (or .neqv. for inequality).

! !REVISION HISTORY:
!  18 Feb 2015 - C. Keller   - Initial Version
!  04 Mar 2015 - R. Yantosca - Now determine if we need to read UV albedo
!                              data from the settings in input.geos
!  16 Mar 2015 - R. Yantosca - Now also toggle TOMS_SBUV_O3 based on
!                              met field type and input.geos settings
!  25 Mar 2015 - C. Keller   - Added switch for STATE_PSC (for UCX)
!  27 Aug 2015 - E. Lundgren - Now always read TOMS for mercury simulation when
!                              photo-reducible HgII(aq) to UV-B radiation is on
!  11 Sep 2015 - E. Lundgren - Remove State_Met and State_Chm from passed args
!  03 Dec 2015 - R. Yantosca - Bug fix: pass am_I_Root to AddExtOpt
!  19 Sep 2016 - R. Yantosca - Now compare LOGICALs with .eqv. and .neqv.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ExtNr
    LOGICAL            :: LTMP
    LOGICAL            :: FOUND

    ! Strings
    CHARACTER(LEN= 31) :: OptName
    CHARACTER(LEN=255) :: MSG, LOC, INS

    !=======================================================================
    ! CheckSettings begins here
    !=======================================================================

    ! Error handling
    LOC  = 'CheckSettings (GeosCore/hcoi_gc_main_mod.F90)'
    INS  = 'HEMCO ERROR: Please check the HEMCO log file for error messages!'

    !-----------------------------------------------------------------------
    ! If emissions shall not be used, reset all extension numbers to -999. 
    ! This will make sure that none of the extensions will be initialized 
    ! and none of the input data related to any of the extensions will be 
    ! used.  The only exception is the NON-EMISSIONS DATA.
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LEMIS ) THEN
       CALL SetExtNr( am_I_Root, HcoConfig, -999, RC=RC )
       IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( 'SetExtNr', LOC, INS )
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
    CALL GetExtOpt( HcoConfig, -999, '+UValbedo+',  OptValBool=LTMP, &
                    FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +UValbedo+', LOC, INS )
    ENDIF

    IF ( FOUND ) THEN

       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LCHEM .AND. ( .NOT. LTMP ) ) THEN
          MSG = 'Setting +UValbedo+ in the HEMCO configuration file ' // &
                'must not be disabled if chemistry is turned on. '    // &
                'If you don`t set that setting explicitly, it will '  // &
                'be set automatically during run-time (recommended)'
          CALL ERROR_STOP( MSG, LOC, INS ) 
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +UValbedo+', LOC, INS )
       ENDIF

    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #2: PSC STATE (for UCX) 
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig, -999, '+STATE_PSC+', OptValBool=LTMP, &
                    FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +STATE_PSC+', LOC, INS )
    ENDIF
    IF ( FOUND ) THEN
       IF ( Input_Opt%LUCX .neqv. LTMP ) THEN
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +STATE_PSC+', LOC, INS )
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
    CALL GetExtOpt( HcoConfig, -999, '+LinStratChem+', OptValBool=LTMP, &
                    FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +LinStratChem+', LOC, INS )
    ENDIF

    IF ( FOUND ) THEN

       ! Print a warning if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LSCHEM .neqv. LTMP ) THEN
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +LinStratChem+', LOC, INS )
       ENDIF

    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #4: TOMS/SBUV overhead O3 columns
    !
    ! If we are using the GEOS-FP met fields, then we will not read in 
    ! the TOMS/SBUV O3 columns unless running a mercury simulation.
    ! We will instead use the O3 columns from the GEOS-FP met fields.  
    ! In this case, we will toggle the +TOMS_SBUV_O3+ collection OFF.
    !
    ! All other met fields use the TOMS/SBUV data in one way or another,
    ! so we will have to read these data from netCDF files.  In this
    ! case, toggle the +TOMS_SBUV_O3+ collection ON if photolysis is
    ! required (i.e. for fullchem/aerosol simulations w/ chemistry on).
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig, -999, '+TOMS_SBUV_O3+', OptValBool=LTMP, &
                    FOUND=FOUND,     RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +TOMS_SBUV_O3+', LOC, INS )
    ENDIF

#if defined( GEOS_FP ) || defined( MERRA2 )
    
    ! Disable for GEOS-FP met fields no matter what it is set to in the 
    ! HEMCO configuration file unless it is a mercury simulation done 
    ! with photo-reducible HgII(aq) to UV-B radiation turned on (jaf)
    IF ( Input_Opt%ITS_A_MERCURY_SIM .and.   &
         Input_Opt%LKRedUV ) THEN
       OptName = '+TOMS_SBUV_O3+ : true'          
    ELSE
       OptName = '+TOMS_SBUV_O3+ : false'          
    ENDIF
    CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC ) 
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'AddExtOpt GEOS-FP +TOMS_SBUV_O3+', LOC, INS )
    ENDIF

#else

    IF ( FOUND ) THEN

       ! Print a warning if this collection is defined in the HEMCO config
       ! file, but is set to a value inconsistent with input.geos file.
       IF ( Input_Opt%LCHEM .neqv. LTMP .and.                       &
            .not. Input_Opt%ITS_A_MERCURY_SIM ) THEN
          WRITE(*,*) ' '
          WRITE(*,*) 'Setting +TOMS_SBUV_O3+ in the HEMCO configuration'
          WRITE(*,*) 'file does not agree with the chemistry settings'
          WRITE(*,*) 'in input.geos. This may be inefficient and/or' 
          WRITE(*,*) 'may yield wrong results!' 
       ENDIF

    ELSE

       ! If this collection is not found in the HEMCO config file, then
       ! activate it only for those simulations that use photolysis 
       ! (e.g. fullchem or aerosol) and only when chemistry is turned on,
       ! or when a mercury simulation is done with photo-reducible
       ! HgII(aq) to UV-B radiation turned on.
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          IF ( Input_Opt%LCHEM ) THEN
             OptName = '+TOMS_SBUV_O3+ : true'
          ENDIF
       ELSE IF ( Input_Opt%ITS_A_MERCURY_SIM .and. &
                 Input_Opt%LKRedUV ) THEN
          OptName = '+TOMS_SBUV_O3+ : true'
       ELSE
          OptName = '+TOMS_SBUV_O3+ : false'
       ENDIF
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +TOMS_SBUV_O3+', LOC, INS )
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
    CALL GetExtOpt( HcoConfig, -999, '+OCEAN_Hg+', OptValBool=LTMP, &
                    FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'GetExtOpt +OCEAN_Hg+', LOC, INS )
    ENDIF

    IF ( FOUND ) THEN
       
       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LDYNOCEAN .AND. ( .NOT. LTMP ) ) THEN
          MSG = 'Setting +OCEAN_Hg+ in the HEMCO configuration file ' // &
                'must not be disabled if chemistry is turned on. '    // &
                'If you don`t set that setting explicitly, it will '  // &
                'be set automatically during run-time (recommended)'
          CALL ERROR_STOP( MSG, LOC, INS ) 
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'AddExtOpt +OCEAN_Hg+', LOC, INS )
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
    USE State_Met_Mod,      ONLY : MetState
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
    USE CMN_SIZE_MOD        ! Size parameters
    USE GC_GRID_MOD,   ONLY : GET_XMID,    GET_YMID_R
    USE PhysConstants
    USE TIME_MOD,      ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,      ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !REMARKS:
!  Moved here from the obsolete global_oh_mod.F.
!
! !REVISION HISTORY: 
!  01 Mar 2013 - C. Keller   - Imported from carbon_mod.F, where it's
!                              called OHNO3TIME
!  16 May 2016 - M. Sulprizio- Remove IJLOOP and change SUNTMP array dimensions
!                              from (MAXIJ) to (IIPAR,JJPAR)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE :: SAVEDOY = -1
    INTEGER       :: I, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(IIPAR,JJPAR)

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
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, YMID_R, TIMLOC, AHR )
          DO J = 1, JJPAR
          DO I = 1, IIPAR

             ! Grid box latitude center [radians]
             YMID_R = GET_YMID_R( I, J, 1 )

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
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) +          &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! SUMCOSZA is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             SUMCOSZA(I,J) = SUMCOSZA(I,J) +             &
                             MAX(SUNTMP(I,J),0e+0_fp)

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
