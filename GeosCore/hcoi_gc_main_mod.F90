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
! \End{itemize}
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
#if !defined(ESMF_) && !defined( MODEL_WRF )
  PRIVATE :: Get_GC_Restart
  PRIVATE :: Get_Met_Fields
#endif
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
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !------------------------------------
  ! %%% Species ID's %%%
  !------------------------------------
  INTEGER             :: id_HNO3
  INTEGER             :: id_LIMO
  INTEGER             :: id_NO
  INTEGER             :: id_NO2
  INTEGER             :: id_O3
  INTEGER             :: id_POPG

  !------------------------------------
  ! %%% Arrays, Pointers, Targets %%%
  !------------------------------------

  ! Internal met fields (will be used by some extensions)
  INTEGER,  TARGET    :: HCO_PBL_MAX                      ! level
  REAL(hp), POINTER   :: HCO_FRAC_OF_PBL(:,:,:)
  REAL(hp), POINTER   :: HCO_SZAFACT(:,:)

  ! Arrays to store J-values (used by Paranox extension)
  REAL(hp), POINTER   :: JNO2(:,:)
  REAL(hp), POINTER   :: JOH(:,:)

  ! Sigma coordinate (temporary)
  REAL(hp), POINTER   :: ZSIGMA(:,:,:)

  ! Sum of cosine of the solar zenith angle. Used to impose a
  ! diurnal variability on OH concentrations
  REAL(fp), POINTER   :: SUMCOSZA(:,:)
!
! !DEFINED PARAMETERS:
!
  ! Temporary toggle for diagnostics
  LOGICAL,  PARAMETER :: DoDiagn = .TRUE.

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
! configuration file (as listed in Input\_Opt%HcoConfigFile) and stored in 
! the HEMCO configuration object. The entire HEMCO setup is based upon the
! entries in the HEMCO configuration object. It is possible to explicitly 
! provide a (previously read) HEMCO configuration object via input argument 
! `HcoConfig`. In this case the HEMCO configuration file will not be read
! any more. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Init( am_I_Root, Input_Opt, State_Met,                  &
                           State_Chm, RC,        HcoConfig                  ) 
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
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
    USE HCO_Config_Mod,     ONLY : Config_ReadFile, ConfigInit
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
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                   :: LSTRAT,  FOUND
    INTEGER                   :: nHcoSpc, HMRC
    INTEGER                   :: N

    ! Strings
    CHARACTER(LEN=255)        :: OptName, ThisLoc, Instr
    CHARACTER(LEN=512)        :: ErrMSg

    ! Pointers
    TYPE(ConfigObj), POINTER  :: iHcoConfig => NULL()
    TYPE(Species),   POINTER  :: SpcInfo
    
    !=======================================================================
    ! HCOI_GC_INIT begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at HCOI_GC_Init (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    ! Define all species ID's here, for use in module routines below
    id_HNO3  = Ind_('HNO3')
    id_LIMO  = Ind_('LIMO')
    id_NO    = Ind_('NO'  )
    id_NO2   = Ind_('NO2' )
    id_O3    = Ind_('O3'  )
    id_POPG  = Ind_('POPG')

    !=======================================================================
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
    !=======================================================================

    ! If HcoConfig is provided
    IF ( PRESENT( HcoConfig ) ) iHcoConfig => HcoConfig

    !---------------------------------------
    ! Initialize HEMCO config object
    !---------------------------------------
    CALL ConfigInit ( iHcoConfig, HMRC, State_Chm%nSpecies )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ConfigInit"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
    ENDIF
    
    ! Pass GEOS-Chem species information to HEMCO config object to
    ! facilitate reading GEOS-Chem restart file via HEMCO
    iHcoConfig%nModelSpc = State_Chm%nSpecies
    DO N = 1, State_Chm%nSpecies
       ! Get info for this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Model ID and species name 
       iHcoConfig%ModelSpc(N)%ModID      = SpcInfo%ModelID
       iHcoConfig%ModelSpc(N)%SpcName    = TRIM( SpcInfo%Name )
    ENDDO

    !---------------------------------------
    ! Phase 1: read settings and switches
    !---------------------------------------
    CALL Config_ReadFile( am_I_Root,               iHcoConfig,               &
                          Input_Opt%HcoConfigFile, 1,          HMRC         )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Config_Readfile" (Phase 1)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
    ENDIF

    ! Check settings
    CALL CheckSettings( am_I_Root, iHcoConfig, Input_Opt,                    &
                        State_Met, State_Chm,  HMRC                         )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "CheckSettings"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !---------------------------------------
    ! Phase 2: read fields
    !---------------------------------------
    CALL Config_ReadFile( am_I_Root,               iHcoConfig,               &
                          Input_Opt%HcoConfigFile, 2,            HMRC       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Config_Readfile" (Phase 2)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !=======================================================================
    ! Open logfile
    !=======================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LOGFILE_OPEN( iHcoConfig%Err, RC=HMRC ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_LogFile_Open"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Initialize HEMCO state object and populate it 
    !=======================================================================

    !-----------------------------------------------------------------------
    ! Extract species to use in HEMCO. nHcoSpc denotes the number of
    ! species that shall be used in HEMCO. The species properties are
    ! defined in the Register_Species call below.
    ! Typically, nHcoSpc is just the number of species defined in both 
    ! the HEMCO configuration file and GEOS-Chem. However, additional
    ! species can be defined, e.g. those not transported in GEOS-Chem
    ! (e.g. SESQ) or tagged species (e.g. specialty simulations).
    !-----------------------------------------------------------------------
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, State_Chm,                    &
                         HcoState,  nHcoSpc,   1,         HMRC              ) 

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = &
         'Error encountered in ""SetHcoSpecies" (first call, to get species)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Now that number of HEMCO species are known, initialize HEMCO
    ! state object.  Links the HEMCO configuration file object 
    ! iHcoConfig to HcoState%Config.
    !-----------------------------------------------------------------------
    CALL HcoState_Init( am_I_Root, HcoState, iHcoConfig, nHcoSpc, HMRC      )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOState_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Register species. This will define all species properties
    ! (names, molecular weights, etc.) of the HEMCO species.
    !-----------------------------------------------------------------------
    CALL SetHcoSpecies ( am_I_Root, Input_Opt, State_Chm,                    &
                         HcoState,  nHcoSpc,   2,         HMRC              )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = &
     'Error encountered in "SetHcoSpecies" (second call, to register species)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set the grid
    !-----------------------------------------------------------------------
    CALL Set_Grid( am_I_Root, State_Met, HcoState, RC                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Set_Grid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Set misc. parameter
    !=======================================================================

    ! Emission, chemistry and dynamics timestep in seconds
    HcoState%TS_EMIS = GET_TS_EMIS()
    HcoState%TS_CHEM = GET_TS_CHEM()
    HcoState%TS_DYN  = GET_TS_DYN()

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

    !=======================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and 
    ! the HEMCO configuration file is removed from buffer in this
    ! step. This also initializes the HEMCO clock as well as the
    ! HEMCO emissions diagnostics collection.
    !=======================================================================
    CALL HCO_Init( am_I_Root, HcoState, HMRC                                )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! Save # of defined dust species in HcoState
    HcoState%nDust                     =  NDSTBIN

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

    !=======================================================================
    ! Initialize all HEMCO extensions. This also selects the required 
    ! met fields used by each extension.
    !=======================================================================
    CALL HCOX_Init( am_I_Root, HcoState, ExtState, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOX_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Update and check logical switches in Input_Opt 
    !-----------------------------------------------------------------------

    ! Soil NOx
    Input_Opt%LSOILNOX      = ( ExtState%SoilNOx > 0 )

    ! Ginoux dust emissions
    IF ( ExtState%DustGinoux > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          ErrMsg = 'DustGinoux is on in HEMCO but LDUST=F in input.geos!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .FALSE.
    ENDIF

    ! DEAD dust emissions
    IF ( ExtState%DustDead > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          ErrMsg = 'DustDead is on in HEMCO but LDUST=F in input.geos!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .TRUE.
    ENDIF

    ! Dust alkalinity
    IF ( ExtState%DustAlk > 0 ) THEN
       IF ( .not. Input_Opt%LDSTUP ) THEN
          ErrMsg = 'DustAlk is on in HEMCO but LDSTUP=F in input.geos'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! Marine organic aerosols
    IF ( ExtState%MarinePOA > 0 ) THEN
       IF ( .not. Input_Opt%LMPOA ) THEN
          ErrMsg = 'MarinePOA is on in HEMCO but LMPOA=F in input.geos'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Set constants for POPs simulation
    !-----------------------------------------------------------------------
    IF ( ExtState%GC_POPs > 0 ) THEN
       ExtState%POP_DEL_H   = Input_Opt%POP_DEL_H
       ExtState%POP_KOA     = Input_Opt%POP_KOA
       ExtState%POP_KBC     = Input_Opt%POP_KBC
       ExtState%POP_DEL_Hw  = Input_Opt%POP_DEL_Hw
       ExtState%POP_XMW     = Input_Opt%POP_XMW
       ExtState%POP_HSTAR   = Input_Opt%POP_HSTAR
    ENDIF

    !-----------------------------------------------------------------------
    ! Initialize ExtState target arrays. 
    ! Extensions typically depend on environmental dependent met. 
    ! variables such as wind speed, surface temp., etc. Pointers 
    ! to these (2D or 3D) fields are defined in the extension object. 
    ! Here, we need to make sure that these pointers are properly 
    ! connected.
    !-----------------------------------------------------------------------
    CALL ExtState_InitTargets( am_I_Root, HcoState, ExtState, HMRC          )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtState_InitTargets"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Define diagnostics
    !-----------------------------------------------------------------------
    IF ( DoDiagn ) THEN 

       ! Set up traditional GEOS-Chem NDxx diagnostics for emissions
       CALL HCOI_GC_Diagn_Init( am_I_Root, Input_Opt,                        &
                                HcoState,  ExtState,  HMRC                  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCOI_GC_Diagn_Init"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Eventually remove pointer
    IF ( PRESENT( HcoConfig ) ) iHcoConfig => NULL()

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
  SUBROUTINE HCOI_GC_Run( am_I_Root, Input_Opt, State_Met,                   &
                          State_Chm, EmisTime,  Phase,     RC               )
!
! !USES:
!
    USE ErrCode_Mod
    USE Get_Ndep_Mod,    ONLY : Reset_Dep_N   ! For soilnox
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState

    ! HEMCO routines 
    USE HCO_Clock_Mod,   ONLY : HcoClock_Get
    USE HCO_Clock_Mod,   ONLY : HcoClock_EmissionsDone
    USE HCO_Diagn_Mod,   ONLY : HcoDiagn_AutoUpdate
    USE HCO_FluxArr_Mod, ONLY : HCO_FluxarrReset 
    USE HCO_Driver_Mod,  ONLY : HCO_Run
    USE HCOX_Driver_Mod, ONLY : HCOX_Run
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    LOGICAL,          INTENT(IN   )  :: EmisTime   ! Is this an emission time step? 
    INTEGER,          INTENT(IN   )  :: Phase      ! Run phase (see remarks)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input options
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Meteo state 
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Phase -1 : Used for GCHP
!  Phase  0 : Simplified Phase 1 for reading initial met fields
!  Phase  1 : Update HEMCO clock and HEMCO data list
!  Phase  2 : Perform emissions calculation
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  22 Aug 2014 - R. Yantosca - Now pass State_Met to MAP_HCO2GC
!  02 Oct 2014 - C. Keller   - PEDGE is now in HcoState%Grid
!  13 Jan 2015 - C. Keller   - Now check if it's time for emissions. Added
!                              call to HcoClock_EmissionsDone.
!  06 Mar 2015 - R. Yantosca - Now create splash page for HEMC
!  06 Jan 2017 - R. Yantosca - Now tell user to look at HEMCO log for err msgs
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
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
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! HCOI_GC_RUN begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at HCOI_GC_Run (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

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
    CALL SetHcoTime( am_I_Root, EmisTime, HMRC                             )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "SetHcoTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! See if it's time for emissions. Don't just use the EmisTime flag in
    ! case that we call this routine multiple times. IsEmisTime will only
    ! be true if this is an emission time step AND emissions have not yet
    ! been calculated for that time step.
    !=======================================================================
    CALL HcoClock_Get( am_I_Root,             HcoState%Clock,                &
                       IsEmisTime=IsEmisTime, RC=HMRC                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HcoClock_Get"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Reset all emission and deposition values. Do this only if it is time
    ! for emissions, i.e. if those values will be refilled.
    !=======================================================================
    IF ( IsEmisTime .AND. Phase == 2 ) THEN
       CALL HCO_FluxArrReset( HcoState, HMRC                                )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_FluxArrReset"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Define pressure edges [Pa] on HEMCO grid.
    ! At Phase 0, the pressure field is not known yet.
    !=======================================================================
    IF ( Phase /= 0 ) THEN
       CALL GridEdge_Set( am_I_Root, State_Met, HcoState, HMRC              )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GridEdge_Set"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
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
    CALL HCO_Run( am_I_Root, HcoState, Phase, HMRC                          )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Run"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

#if !defined(ESMF_) && !defined( MODEL_WRF )
    !=======================================================================
    ! Get met fields from HEMCO
    !=======================================================================
    CALL Get_Met_Fields( am_I_Root, Input_Opt, State_Met, State_Chm, &
                         Phase, RC )

    !=======================================================================
    ! Get fields from GEOS-Chem restart file
    !=======================================================================
    IF ( Phase == 0 ) THEN
       CALL Get_GC_Restart( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    ENDIF
#endif

    !=======================================================================
    ! Do the following only if it's time to calculate emissions 
    !=======================================================================
    IF ( Phase == 2 .AND. IsEmisTime ) THEN 

       !--------------------------------------------------------------------
       ! Set / update ExtState fields.
       ! Extensions typically depend on environmental dependent met. 
       ! variables such as wind speed, surface temp., etc. Pointers 
       ! to these (2D or 3D) fields are defined in the extension object. 
       ! Here, we need to make sure that these pointers are properly 
       ! connected.
       !--------------------------------------------------------------------
       CALL ExtState_SetFields( am_I_Root, State_Met, State_Chm,             &
                                HcoState,  ExtState,  HMRC                  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtState_SetFields"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
   
       CALL ExtState_UpdateFields( am_I_Root, Input_Opt, State_Met, State_Chm,&
                                   HcoState,  ExtState,  HMRC               )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtState_UpdateFields"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
   
       !====================================================================
       ! Run HCO extensions. Emissions will be added to corresponding
       ! flux arrays in HcoState.
       !====================================================================
       CALL HCOX_Run( am_I_Root, HcoState, ExtState, HMRC                   )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCOX_Run"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
   
       !====================================================================
       ! Update all 'AutoFill' diagnostics. This makes sure that all 
       ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
       ! AutoFill flag is specified when creating a diagnostics container
       ! (Diagn_Create).
       !====================================================================
       IF ( DoDiagn ) THEN
          CALL HcoDiagn_AutoUpdate( am_I_Root, HcoState, HMRC               )

          ! Trap potential errors
          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "HcoDiagn_AutoUpdate"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             CALL Flush( HcoState%Config%Err%Lun )
             RETURN
          ENDIF
       ENDIF
   
       !====================================================================
       ! Reset the accumulated nitrogen dry and wet deposition to zero. 
       ! Will be re-filled in drydep and wetdep.
       !====================================================================
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          CALL RESET_DEP_N( State_Chm )
       ENDIF

       !====================================================================
       ! Emissions are now done for this time step
       !====================================================================
       CALL HcoClock_EmissionsDone( am_I_Root, HcoState%Clock, HMRC         )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HcoClock_EmissionsDone"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       
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
  SUBROUTINE HCOI_GC_Final( am_I_Root, Error, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,     ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE HCO_Driver_Mod,   ONLY : HCO_Final
    USE HCO_Diagn_Mod,    ONLY : DiagnBundle_Cleanup
    USE HCO_State_Mod,    ONLY : HcoState_Final
    USE HCOX_Driver_Mod,  ONLY : HCOX_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU
    LOGICAL, INTENT(IN)  :: Error       ! Cleanup after exit?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  19 Feb 2015 - R. Yantosca - Change restart file name back to HEMCO_restart
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! HCOI_GC_Final begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at HCOI_GC_Final (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Cleanup HEMCO core 
    !-----------------------------------------------------------------------
    CALL HCO_Final( am_I_Root, HcoState, Error, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    !-----------------------------------------------------------------------
    CALL HCOX_Final( am_I_Root, HcoState, ExtState, HMRC ) 

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOX_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup diagnostics
    !-----------------------------------------------------------------------
    CALL DiagnBundle_Cleanup( HcoState%Diagn )

    !-----------------------------------------------------------------------
    ! Cleanup HcoState object 
    !-----------------------------------------------------------------------
    CALL HcoState_Final( HcoState ) 

    !-----------------------------------------------------------------------
    ! Deallocate module variables
    !-----------------------------------------------------------------------
    IF ( ASSOCIATED( ZSIGMA ) ) THEN
       DEALLOCATE( ZSIGMA, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:ZSIGMA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( HCO_FRAC_OF_PBL ) ) THEN
       DEALLOCATE( HCO_FRAC_OF_PBL, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_FRAC_OF_PBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( HCO_SZAFACT ) ) THEN
       DEALLOCATE( HCO_SZAFACT, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_SZAFACT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( JNO2 ) ) THEN
       DEALLOCATE( JNO2, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JNO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( JOH ) ) THEN
       DEALLOCATE( JOH, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( SUMCOSZA ) ) THEN
       DEALLOCATE( SUMCOSZA, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:SUMCOSZA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
 
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
    USE HCOIO_Diagn_Mod, ONLY : HcoDiagn_Write 
    USE Input_Opt_Mod,   ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   ) :: am_I_Root    ! Root CPU?
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt    ! Input options
    LOGICAL,        INTENT(IN   ) :: Restart      ! write restart (enforced)?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY: 
!  01 Apr 2015 - C. Keller   - Initial version 
!  06 Jan 2017 - R. Yantosca - Now tell user to check HEMCO log for err msgs
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: HMRC

    ! Strings
    CHARACTER(LEN=255)  :: ThisLoc, Instr
    CHARACTER(LEN=512)  :: ErrMsg

    !=======================================================================
    ! HCOI_GC_WriteDiagn begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at HCOI_GC_WriteDiagn (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Make sure HEMCO time is in sync 
    !-----------------------------------------------------------------------
    CALL SetHcoTime ( am_I_Root, .FALSE., HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "SetHcoTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Write diagnostics
    !-----------------------------------------------------------------------
    CALL HcoDiagn_Write( am_I_Root, HcoState, RESTART, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HcoDiagn_Write"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

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
    USE CMN_SIZE_Mod, ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE HCO_Arr_Mod,  ONLY : HCO_ArrAssert
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
!  19 Jan 2018 - R. Yantosca  - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_InitTargets begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_InitTargets (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! HCO_SZAFACT is not defined in Met_State.  Hence need to 
    ! define here so that we can point to them.
    !
    ! Now include HCO_FRAC_OF_PBL and HCO_PBL_MAX for POPs specialty
    ! simulation (mps, 8/20/14)
    ! ----------------------------------------------------------------------
    IF ( ExtState%SZAFACT%DoUse ) THEN 

       ALLOCATE( SUMCOSZA( IIPAR, JJPAR ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:SUMCOSZA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMCOSZA = 0.0_fp

       ALLOCATE( HCO_SZAFACT( IIPAR, JJPAR ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_SZAFACT', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       HCO_SZAFACT = 0e0_hp

    ENDIF

    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN 
       ALLOCATE( HCO_FRAC_OF_PBL( IIPAR, JJPAR, LLPAR ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_FRAC_OF_PBL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       HCO_FRAC_OF_PBL = 0.0_hp
    ENDIF

    ! Initialize max. PBL
    HCO_PBL_MAX = 0

    ! ----------------------------------------------------------------------
    ! The J-Values for NO2 and O3 are not defined in Met_State. We
    ! need to compute them separately. 
    ! ----------------------------------------------------------------------
    IF ( ExtState%JNO2%DoUse ) THEN 
       ALLOCATE( JNO2( IIPAR, JJPAR ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JNO2', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       JNO2 = 0.0e0_hp
    ENDIF

    IF ( ExtState%JOH%DoUse ) THEN 
       ALLOCATE( JOH( IIPAR, JJPAR ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JOH', 0, RC )
       JOH = 0.0e0_hp
    ENDIF

    ! ----------------------------------------------------------------------
    ! Arrays to be copied physically because HEMCO units are not the
    ! same as in GEOS-Chem
    ! ----------------------------------------------------------------------

    ! TROPP: GEOS-Chem TROPP is in hPa, while HEMCO uses Pa.
    IF ( ExtState%TROPP%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%TROPP%Arr, IIPAR, JJPAR, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( TROPP )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! SPHU: GEOS-Chem SPHU is in g/kg, while HEMCO uses kg/kg.
    ! NOTE: HEMCO only uses SPHU surface values.
    IF ( ExtState%SPHU%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SPHU%Arr, IIPAR, JJPAR, 1, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( SPHU )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! SUNCOS: HEMCO now calculates SUNCOS values based on its own
    ! subroutine 
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, IIPAR, JJPAR, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( SUNCOS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
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
! State\_Met, State\_Chm, etc. For example, if the HEMCO data list contains 
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
  SUBROUTINE ExtState_SetFields( am_I_Root, State_Met, State_Chm,            &
                                 HcoState,  ExtState,  RC                   ) 
!
! !USES:
!
    USE Hcox_State_Mod, ONLY : ExtDat_Set
    USE ErrCode_Mod
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState

    ! For SoilNox
    USE Drydep_Mod,     ONLY : DryCoeff

#if defined(ESMF_)
    USE HCOI_Esmf_Mod,  ONLY : HCO_SetExtState_ESMF
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
!  06 Jan 2017 - R. Yantosca  - Now tell user to look at HEMCO log for err msgs
!  23 May 2017 - R. Yantosca  - Fixed typo for MERRA2 in #ifdef at line 1193
!  02 Jun 2017 - C. Keller    - Call HCO_SetExtState_ESMF every time.
!  19 Jan 2018 - R. Yantosca  - Now return error code to calling program
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
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_SetFields begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_SetFields (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Pointers to local module arrays 
    !-----------------------------------------------------------------------

    ! SZAFACT
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%SZAFACT,         & 
                    'SZAFACT_FOR_EMIS', HMRC,     FIRST,                    &
                     HCO_SZAFACT                                           )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SZAFACT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JNO2
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%JNO2,             &
                    'JNO2_FOR_EMIS',    HMRC,     FIRST,                     &
                     JNO2                                                   )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JNO2_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JOH
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%JOH,              &
                    'JOH_FOR_EMIS',     HMRC,     FIRST,                     &
                     JOH                                                    )  

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JOH_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Frac of PBL
    CALL ExtDat_Set( am_I_Root,             HcoState,                        &
                     ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS',           &
                     HMRC,                  FIRST,                           &
                     HCO_FRAC_OF_PBL                                        )  

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRAC_OF_PBL_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! 2D fields 
    ! ----------------------------------------------------------------------

    ! U10M
    CALL ExtDat_Set( am_I_Root,      HcoState, ExtState%U10M,                &
                    'U10M_FOR_EMIS',  HMRC   , FIRST,                        &
                     State_Met%U10M                                         )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( U10M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! V10M
    CALL ExtDat_Set( am_I_Root,      HcoState, ExtState%V10M,                &
                    'V10M_FOR_EMIS', HMRC,     FIRST,                        & 
                     State_Met%V10M                                         )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( V10M_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! ALBD
    CALL ExtDat_Set( am_I_Root,      HcoState, ExtState%ALBD,                &
                    'ALBD_FOR_EMIS', HMRC,     FIRST,                        &
                     State_Met%ALBD                                         )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( ALBD_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! WLI
    CALL ExtDat_Set( am_I_Root,     HcoState, ExtState%WLI,                  &
                    'WLI_FOR_EMIS', HMRC,     FIRST,                         &
                     State_Met%LWI                                          )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( WLI_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( am_I_Root,     HcoState, ExtState%T2M,                  &
                    'T2M_FOR_EMIS', HMRC,     FIRST,                         &
                     State_Met%TS                                           )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( T2M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! TSKIN
    CALL ExtDat_Set( am_I_Root,       HcoState, ExtState%TSKIN,              &
                    'TSKIN_FOR_EMIS', HMRC,     FIRST,                       &
                     State_Met%TSKIN                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TSKIN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETROOT
    CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%GWETROOT,        &
                    'GWETROOT_FOR_EMIS', HMRC,     FIRST,                    &
                     State_Met%GWETROOT                                     )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETROOT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETTOP
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%GWETTOP,          &
                    'GWETTOP_FOR_EMIS', HMRC,     FIRST,                     &
                    State_Met%GWETTOP                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETTOP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( am_I_Root,       HcoState, ExtState%USTAR,              &
                    'USTAR_FOR_EMIS', HMRC,     FIRST,                       &
                     State_Met%USTAR                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( USTAR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Z0
    CALL ExtDat_Set( am_I_Root,    HcoState, ExtState%Z0,                    &
                    'Z0_FOR_EMIS', HMRC,     FIRST,                          &
                     State_Met%Z0                                           )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( Z0_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( am_I_Root,       HcoState, ExtState%PARDR,              &
                    'PARDR_FOR_EMIS', HMRC,     FIRST,                       &
                     State_Met%PARDR                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PARDF
    CALL ExtDat_Set( am_I_Root,        HcoState, ExtState%PARDF,             &
                    'PARDF_FOR_EMIS',  HMRC, FIRST,                          &
                     State_Met%PARDF                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDF_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PSC2_WET
    CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%PSC2_WET,        &
                    'PSC2_WET_FOR_EMIS', HMRC,     FIRST,                    &
                     State_Met%PSC2_WET                                     )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PSC2_WET_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! RADSWG
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%RADSWG,           &
                    'RADSWG_FOR_EMIS',  HMRC,     FIRST,                     &
                     State_Met%SWGDN                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( RADSWG_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRCLND
    CALL ExtDat_Set( am_I_Root,        HcoState, ExtState%FRCLND,            &
                    'FRCLND_FOR_EMIS', HMRC,     FIRST,                      &
                     State_Met%FRCLND                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRCLND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! CLDFRC
    CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%CLDFRC,          &
                     'CLDFRC_FOR_EMIS',  HMRC,     FIRST,                    &
                     State_Met%CLDFRC                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CLDFRC_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWHGT is is mm H2O, which is the same as kg H2O/m2.
    ! This is the unit of SNOMAS.
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%SNOWHGT,          &
                    'SNOWHGT_FOR_EMIS', HMRC,     FIRST,                     &
                     State_Met%SNOMAS                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWHGT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWDP is in m
    CALL ExtDat_Set( am_I_Root,        HcoState, ExtState%SNODP,             &
                    'SNODP_FOR_EMIS',  HMRC,    FIRST,                       & 
                     State_Met%SNODP                                        )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWDP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLAND
    CALL ExtDat_Set( am_I_Root,        HcoState, ExtState%FRLAND,            &
                    'FRLAND_FOR_EMIS', HMRC,     FIRST,                      &
                     State_Met%FRLAND                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FROCEAN
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%FROCEAN,          &
                    'FROCEAN_FOR_EMIS', HMRC,     FIRST,                     &
                     State_Met%FROCEAN                                      )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FROCEAN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLAKE
    CALL ExtDat_Set( am_I_Root,        HcoState, ExtState%FRLAKE,            &
                    'FRLAKE_FOR_EMIS', HMRC,     FIRST,                      &
                     State_Met%FRLAKE                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAKE_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLANDIC
    CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%FRLANDIC,        &
                    'FRLANDIC_FOR_EMIS', HMRC,     FIRST,                    &
                     State_Met%FRLANDIC                                     )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLANDIC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! LAI
    CALL ExtDat_Set( am_I_Root,      HcoState, ExtState%LAI,                 &
                    'LAI_FOR_EMIS',  HMRC,     FIRST,                        &
                     State_Met%MODISLAI  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( LAI_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! CHLR
    CALL ExtDat_Set( am_I_Root,           HcoState, ExtState%CHLR,           &
                    'CHLR',               HMRC,     FIRST,                   &
                     State_Met%MODISCHLR                                    )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CHLR )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Convective fractions
    CALL ExtDat_Set( am_I_Root,         HcoState,        ExtState%CNV_FRC,   &
                    'CNV_FRC_FOR_EMIS', HMRC,            FIRST,              &
                     State_Met%CNV_FRC, NotFillOk=.TRUE.                    )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_FRC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! 3D fields 
    !-----------------------------------------------------------------------

    ! CNV_MFC
    CALL ExtDat_Set( am_I_Root,         HcoState,        ExtState%CNV_MFC,   &
                    'CNV_MFC_FOR_EMIS', HMRC,            FIRST,              &
                     State_Met%CMFMC,   OnLevEdge=.TRUE.                    )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_MFC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( am_I_Root,    HcoState, ExtState%TK,                    &
                    'TK_FOR_EMIS', HMRC,     FIRST,                          &
                     State_Met%T                                            )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TK_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Air mass [kg/grid box]
    CALL ExtDat_Set( am_I_Root,     HcoState, ExtState%AIR,                  &
                    'AIR_FOR_EMIS', HMRC,     FIRST,                         &
                     State_Met%AD                                           )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! AIRVOL_FOR_EMIS
    CALL ExtDat_Set( am_I_Root,         HcoState, ExtState%AIRVOL,           &
                    'AIRVOL_FOR_EMIS',  HMRC,     FIRST,                     &
                     State_Met%AIRVOL                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRVOL_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Dry air density [kg/m3]
    CALL ExtDat_Set( am_I_Root, HcoState, ExtState%AIRDEN,                   &
                    'AIRDEN',   HMRC,     FIRST,                             &
                     State_Met%AIRDEN                                       )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRDEN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
 
    ! ----------------------------------------------------------------
    ! Species concentrations
    ! ----------------------------------------------------------------
    IF ( id_O3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_O3)
       CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%O3,           &
                       'HEMCO_O3_FOR_EMIS', HMRC,     FIRST,         Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_O3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_NO2 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO2)
       CALL ExtDat_Set( am_I_Root,           HcoState, ExtState%NO2,         &
                       'HEMCO_NO2_FOR_EMIS', HMRC,     FIRST,        Trgt3D ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO2_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_NO > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO)
       CALL ExtDat_Set( am_I_Root,          HcoState, ExtState%NO,           &
                       'HEMCO_NO_FOR_EMIS', HMRC,     FIRST,         Trgt3D ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_HNO3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_HNO3)
       CALL ExtDat_Set( am_I_Root,            HcoState, ExtState%HNO3,       &
                       'HEMCO_HNO3_FOR_EMIS', HMRC,     FIRST,       Trgt3D ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_HNO3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_POPG > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_POPG)
       CALL ExtDat_Set( am_I_Root,            HcoState, ExtState%POPG,       &
                       'HEMCO_POPG_FOR_EMIS', HMRC,     FIRST,       Trgt3D ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_POPG_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    ! ----------------------------------------------------------------------
    ! Deposition parameter
    ! ----------------------------------------------------------------------
    
    ! DRY_TOTN
    CALL ExtDat_Set( am_I_Root,           HcoState, ExtState%DRY_TOTN,       &
                    'DRY_TOTN_FOR_EMIS',  HMRC,     FIRST,                   &
                     State_Chm%DryDepNitrogen                               ) 

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( DRY_TOTN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! WET_TOTN
    CALL ExtDat_Set( am_I_Root,           HcoState, ExtState%WET_TOTN,       &
                    'WET_TOTN_FOR_EMIS',  HMRC,     FIRST,                   &
                     State_Chm%WetDepNitrogen                               ) 

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( WET_TOTN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

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
    ! ckeller, 06/02/17: now call this on every time step. Routine
    ! HCO_SetExtState_ESMF copies the fields to ExtState.
    ! ----------------------------------------------------------------
#if defined( ESMF_ )
    ! IF ( FIRST ) THEN
    CALL HCO_SetExtState_ESMF ( am_I_Root, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_SetExtState_ESMF"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
    !ENDIF
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
  SUBROUTINE ExtState_UpdateFields( am_I_Root, Input_Opt, State_Met,  &
                                    State_Chm, HcoState,  ExtState,  RC ) 
!
! !USES:
!
    USE CMN_FJX_MOD,      ONLY : ZPJ
    USE CMN_SIZE_MOD,     ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE FAST_JX_MOD,      ONLY : RXN_NO2, RXN_O3_1, RXN_O3_2a
    USE HCO_GeoTools_Mod, ONLY : HCO_GetSUNCOS
    USE Input_Opt_Mod,    ONLY : OptInput
    USE PBL_MIX_MOD,      ONLY : GET_FRAC_OF_PBL
    USE PBL_MIX_MOD,      ONLY : GET_PBL_MAX_L
    USE State_Met_Mod,    ONLY : MetState
    USE State_Chm_Mod,    ONLY : ChmState
#if defined(ESMF_) 
    USE HCOI_ESMF_MOD,    ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input opts
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
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_UpdateFields begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_UpdateFields (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !=======================================================================
    ! Update fields in the HEMCO Extension state
    !=======================================================================

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
       CALL HCO_GetSUNCOS( am_I_Root,               HcoState,                &
                           ExtState%SUNCOS%Arr%Val, 0,        HMRC          ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_GetSuncos"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
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
                IF ( Input_Opt%LUCX ) THEN
                   ! RXN_O3_1: O3  + hv --> O2  + O
                   JOH(I,J) = ZPJ(L,RXN_O3_1,I,J)
                ELSE
                   ! RXN_O3_2a: O3 + hv --> 2OH
                   JOH(I,J) = ZPJ(L,RXN_O3_2a,I,J)
                ENDIF
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
  SUBROUTINE GridEdge_Set( am_I_Root, State_Met, HcoState, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Met_Mod,    ONLY : MetState
    USE HCOX_STATE_MOD,   ONLY : ExtDat_Set
    USE HCO_GeoTools_Mod, ONLY : HCO_CalcVertGrid
    USE HCO_GeoTools_Mod, ONLY : HCO_SetPBLm
    USE CMN_SIZE_MOD,     ONLY : IIPAR, JJPAR, LLPAR
    USE PBL_MIX_MOD,      ONLY : GET_PBL_TOP_M
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
    ! Scalars
    INTEGER            :: I, J
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    ! Pointers
    REAL(hp), POINTER  :: PBLM    (:,:  )    ! PBL height           [m ] 
    REAL(hp), POINTER  :: BXHEIGHT(:,:,:)    ! Grid box height      [m ]
    REAL(hp), POINTER  :: PEDGE   (:,:,:)    ! Pressure @ lvl edges [Pa]
    REAL(hp), POINTER  :: PSFC    (:,:  )    ! Surface pressure     [Pa]
    REAL(hp), POINTER  :: TK      (:,:,:)    ! Temperature          [K ]
    REAL(hp), POINTER  :: ZSFC    (:,:  )    ! Surface geopotential [m ]

    !=======================================================================
    ! GridEdge_Set begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at GridEdge_Set (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Allocate the PEDGE array, which holds level edge pressures [Pa]
    ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
    ! make PEDGE be a pointer and allocate/deallocate it on each call.
    !-----------------------------------------------------------------------
    ALLOCATE( PEDGE( IIPAR, JJPAR, LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PEDGE', 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Edge and surface pressures [Pa]
    PEDGE    =  State_Met%PEDGE * 100.0_hp  ! Convert hPa -> Pa
    PSFC     => PEDGE(:,:,1)

    ! Point to other fields of State_Met
    ZSFC     => State_Met%PHIS               
    BXHEIGHT => State_Met%BXHEIGHT           
    TK       => State_Met%T                  

    !-----------------------------------------------------------------------
    ! Calculate vertical grid properties
    !-----------------------------------------------------------------------
    CALL HCO_CalcVertGrid( am_I_Root, HcoState, PSFC,  ZSFC,                &
                           TK,        BXHEIGHT, PEDGE, HMRC                )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_CalcVertGrid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set PBL heights
    !-----------------------------------------------------------------------
    ALLOCATE( PBLM( IIPAR, JJPAR ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PBLM', 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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
    CALL HCO_SetPBLm( am_I_Root, HcoState,         FldName='PBL_HEIGHT',     &
                      PBLM=PBLM, DefVal=1000.0_hp, RC=HMRC                  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_SetPblM"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup and quit
    !-----------------------------------------------------------------------

    ! Deallocate and PEDGE
    IF ( ASSOCIATED( PEDGE ) ) THEN
       DEALLOCATE( PEDGE, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PEDGE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Deallocate PBLM
    IF ( ASSOCIATED( PBLM ) ) THEN
       DEALLOCATE( PBLM  )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PBLM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Free pointers
    ZSFC     => NULL()
    BXHEIGHT => NULL()
    TK       => NULL()
    PSFC     => NULL()
    PEDGE    => NULL()
    PBLM     => NULL()

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
  SUBROUTINE SetHcoSpecies( am_I_Root, Input_Opt, State_Chm,                 &
                            HcoState,  nSpec,     Phase,     RC             )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_LogFile_Mod, ONLY : HCO_SPEC2LOG
    USE Input_Opt_Mod,   ONLY : OptInput
    USE Species_Mod,     ONLY : Species
    USE HCO_Types_Mod,   ONLY : ConfigObj
    USE State_Chm_Mod,   ONLY : ChmState
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
    INTEGER                :: nSpc, HMRC
    INTEGER                :: N,    L,    M
    REAL(dp)               :: K0,   CR,   pKa

    ! Strings
    CHARACTER(LEN= 31)     :: ThisName
    CHARACTER(LEN=255)     :: ThisLoc, Instr
    CHARACTER(LEN=512)     :: ErrMsg,  Msg

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! SetHcoSpecies begins here
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at SetHcoSpecies (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

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
         Input_Opt%ITS_A_TAGCO_SIM      .or. &
         Input_Opt%ITS_A_CO2_SIM              ) THEN

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
       ! Add 5 extra species (ISOP, ACET, MTPA, LIMO, MTPO) for tagged CO 
       IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
          nSpc = nSpc + 5
       ENDIF

       ! Assign species variables
       IF ( PHASE == 2 ) THEN

          ! Verbose
          IF ( am_I_Root ) THEN
             Msg = 'Registering HEMCO species:'
             CALL HCO_MSG( HcoState%Config%Err, Msg )
          ENDIF

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(ErrMsg,*) 'Input species /= expected species: ', nSpec, nSpc 
             CALL HCO_ERROR( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
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
          ! Add the non-advected species ISOP, ACET, MTPA, LIMO, MTPO
          ! in the last 5 species slots (bmy, ckeller, 6/1/16)
          !------------------------------------------------------------------
          IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
       
             ! Add 5 additional species
             DO L = 1, 5
                
                ! ISOP, ACET, MONX follow the regular tagged CO species
                M = State_Chm%nAdvect + L

                ! Get the species name
                SELECT CASE( L )
                   CASE( 1 ) 
                      ThisName = 'ISOP'
                   CASE( 2 )
                      ThisName = 'ACET'
                   CASE( 3 )
                      ThisName = 'MTPA'
                   CASE( 4 )
                      ThisName = 'LIMO'
                   CASE( 5 )
                      ThisName = 'MTPO'
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
          IF ( am_I_Root ) CALL HCO_MSG( HcoState%Config%Err, SEP1='-' )
       ENDIF ! Phase = 2   

    !-----------------------------------------------------------------
    ! DEFAULT (RETURN W/ ERROR) 
    !-----------------------------------------------------------------
    ELSE
       ErrMsg = 'Invalid simulation type - cannot define model species' 
       CALL HCO_ERROR( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
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
    USE ErrCode_Mod
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
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: L
    INTEGER               :: HMRC

    ! Arrays
    REAL(hp), ALLOCATABLE :: Ap(:),   Bp(:)

    ! Strings
    CHARACTER(LEN=255)    :: ThisLoc, Instr
    CHARACTER(LEN=512)    :: ErrMsg

    !=======================================================================
    ! SET_Grid begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Set_Grid (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !=======================================================================
    ! NOTE: for now, just copy GEOS-Chem grid, i.e. HEMCO calculations 
    ! are performed on the GEOS-Chem simulation grid. 
    ! It is possible to define a different emissions grid below. 
    ! In this case, all arrays have to be regridded when passing 
    ! them between HEMCO and GEOS-Chem (this is also true for the 
    ! met-fields used by the extensions)! 
    !=======================================================================

    ! Grid dimensions
    HcoState%NX = IIPAR
    HcoState%NY = JJPAR
    HcoState%NZ = LLPAR

    ! Allocate Ap array
    ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod:Set_Grid:Ap', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Allocate Bp array
    ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod:Set_Grid:Bp', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Get Ap and Bp values from GEOS-Chem
    DO L = 1, LLPAR+1
       Ap(L) = GET_AP(L) * 100_hp ! hPa to Pa 
       Bp(L) = GET_BP(L)          ! unitless
    ENDDO

    ! Define the vertical grid
    CALL HCO_VertGrid_Define( am_I_Root, HcoState%Config,                    &
                              zGrid      = HcoState%Grid%zGrid,              &
                              nz         = LLPAR,                            &
                              Ap         = Ap,                               & 
                              Bp         = Bp,                               & 
                              RC         = HMRC                             )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_VertGrid_Define"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

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
  SUBROUTINE CheckSettings( am_I_Root, HcoConfig, Input_Opt,                 &
                            State_Met, State_Chm, RC                        )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Types_Mod,      ONLY : ConfigObj
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
!  19 Jan 2018 - R. Yantosca - Now return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ExtNr
    INTEGER            :: HMRC
    LOGICAL            :: LTMP
    LOGICAL            :: FOUND

    ! Strings
    CHARACTER(LEN=31 ) :: OptName
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! CheckSettings begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at CheckSettings (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! If emissions shall not be used, reset all extension numbers to -999. 
    ! This will make sure that none of the extensions will be initialized 
    ! and none of the input data related to any of the extensions will be 
    ! used.  The only exception is the NON-EMISSIONS DATA.
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LEMIS ) THEN
       CALL SetExtNr( am_I_Root, HcoConfig, -999, RC=HMRC                   )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "SetExtNr"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
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
    CALL GetExtOpt( HcoConfig,       -999,       '+UValbedo+',               &
                    OptValBool=LTMP, FOUND=FOUND, RC=HMRC                   )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( +UValbedo+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    IF ( FOUND ) THEN

       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LCHEM .AND. ( .NOT. LTMP ) ) THEN
          ErrMsg= 'Setting +UValbedo+ in the HEMCO configuration file '   // &
                  'must not be disabled if chemistry is turned on. '      // &
                  'If you don`t set that setting explicitly, it will '    // &
                  'be set automatically during run-time (recommended)'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( +UValbedo+ )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

    ENDIF 

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #2: PSC STATE (for UCX) 
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig,       -999,        '+STATE_PSC+',             &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC                  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( +STATE_PSC+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    IF ( FOUND ) THEN
#if defined( ESMF_ )
       ! If this is in an ESMF environment, then we should not get STATE_PSC
       ! through HEMCO - instead it is an internal restart variable
       If (LTMP) Then
          ErrMsg = 'Error encountered in "GetExtOpt( +STATE_PSC+ in ESMF )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
#else
       IF ( Input_Opt%LUCX .neqv. LTMP ) THEN
          ErrMsg = 'Setting +STATE_PSC+ in the HEMCO configuration'       // &
                   'file does not agree with stratospheric chemistry'     // &
                   'settings in input.geos. This may be inefficient'      // &
                   'and/or yield to wrong results!'
          CALL GC_Warning( ErrMsg, RC, ThisLoc )
       ENDIF
#endif
    ELSE
#if defined( ESMF_ )
       OptName = '+STATE_PSC+ : false'
#else
       IF ( Input_Opt%LUCX ) THEN
          OptName = '+STATE_PSC+ : true'
       ELSE
          OptName = '+STATE_PSC+ : false'
       ENDIF
#endif
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt(  +STATE_PSC+ )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF 

#if defined( ESMF_ )
    ! Also check that HEMCO_RESTART is not set
    CALL GetExtOpt( HcoConfig,       -999,        'HEMCO_RESTART',           &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC                  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( HEMCO_RESTART in ESMF)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    If ( FOUND .and. LTMP ) Then
       ErrMsg = 'Error encountered in "ESMF HEMCO_RESTART"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    End If
#endif

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
    CALL GetExtOpt( HcoConfig,      -999,         '+LinStratChem+',          &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC                  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( +LinStratChem+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    IF ( FOUND ) THEN

       ! Print a warning if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LSCHEM .neqv. LTMP ) THEN
          ErrMsg = 'Setting +LinStratChem+ in the HEMCO configuration'    // & 
                   'file does not agree with stratospheric chemistry'     // &
                   'settings in input.geos. This may be inefficient'      // &
                   'and/or may yield wrong results!' 
          CALL GC_Warning( ErrMsg, RC, ThisLoc )
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC ) 

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( +LinStratChem+ )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
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
    CALL GetExtOpt( HcoConfig,       -999,       '+TOMS_SBUV_O3+',           &
                    OptValBool=LTMP, FOUND=FOUND, RC=HMRC                   )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( +TOMS_SBUV_O3+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Disable TOMS/SBUV O3 no matter what it is set to in the 
    ! HEMCO configuration file unless it is a mercury simulation done 
    ! with photo-reducible HgII(aq) to UV-B radiation turned on (jaf)
    IF ( Input_Opt%ITS_A_MERCURY_SIM .and.   &
         Input_Opt%LKRedUV ) THEN
       OptName = '+TOMS_SBUV_O3+ : true'          
    ELSE
       OptName = '+TOMS_SBUV_O3+ : false'          
    ENDIF
    CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC    ) 

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "AddExtOpt( +TOMS_SBUV_O3+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #5: Ocean Hg input data (for Hg sims only)
    !
    ! If we have turned on the Ocean Mercury simulation in the
    ! input.geos file, then we will also toggle the +OCEAN_Hg+ 
    ! collection so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig,       -999,       '+OCEAN_Hg+',               &
                    OptValBool=LTMP, FOUND=FOUND, RC=HMRC                   )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( +OCEAN_Hg+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    IF ( FOUND ) THEN
       
       ! Stop the run if this collection is defined in the HEMCO config
       ! file, but is set to an value inconsistent with input.geos file.
       IF ( Input_Opt%LDYNOCEAN .AND. ( .NOT. LTMP ) ) THEN
          ErrMsg = 'Setting +OCEAN_Hg+ in the HEMCO configuration file '  // &
                   'must not be disabled if chemistry is turned on. '     // &
                   'If you don`t set that setting explicitly, it will '   // &
                   'be set automatically during run-time (recommended)'
          CALL GC_Warning( ErrMsg, RC, ThisLoc )
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
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( +OCEAN_Hg+ )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! NON-EMISSIONS DATA #6: RRTMG input data
    !
    ! If we have turned on the Ocean Mercury simulation in the
    ! input.geos file, then we will also toggle the +OCEAN_Hg+ 
    ! collection so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig,       -999,       '+RRTMG+',                  &
                    OptValBool=LTMP, FOUND=FOUND, RC=HMRC                   )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       RETURN
       ErrMsg = 'Error encountered in "GetExtOpt( +RRTMG+ )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )

    ENDIF
    
    IF ( FOUND ) THEN

       ! If this collection is explicitly found in the HEMCO_Config file,
       ! but RRTMG is turned off, then throw an error and stop the run
       IF ( ( .not. Input_Opt%LRAD               )   .and.                   &
            ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) ) THEN 

          ErrMsg = 'Setting +RRTMG+ explicitly in the HEMCO '             // &
                   'configuration file must only be done if the '         // &
                   'RRTMG radiative transfer model is turned on, and '    // &
                   'GEOS-Chem is using one of the full-chemistry '        // &
                   'simulations.'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

    ELSE
          
       ! If this collection is not explicitly found in the HEMCO Config file,
       ! then turn it on if (1) RRTMG is turned on, and (2) the current
       ! simulation is one of the full-chemistry simulations (bmy, 10/31/18)
       IF ( Input_Opt%LRAD .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          OptName = '+RRTMG+ : true'
       ELSE
          OptName = '+RRTMG+ : false'
       ENDIF
       CALL AddExtOpt( am_I_Root, HcoConfig, TRIM(OptName), CoreNr, RC=HMRC  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( +RRTMG+ )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
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
!  01 Mar 2013 - C. Keller   - Imported from carbon_mod.F, where these
!                              calculations are done w/in GET_OH
!  06 Feb 2018 - E. Lundgren - Update unit conversion factor for timestep
!                              unit changed from min to sec
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
              * ( 86400e+0_fp / GET_TS_CHEM() )

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
  SUBROUTINE Calc_SumCosZa()
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
!  06 Feb 2018 - E. Lundgren - Update time conversion factors for timestep
!                              unit change from min to sec
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
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()      

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day 
       DO N = 1, NDYSTEP
            
          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP DEFAULT( SHARED )
!!$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
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
         NT = NT + GET_TS_CHEM()             
      ENDDO

      ! Set saved day of year to current day of year 
      SAVEDOY = GET_DAY_OF_YEAR()

   ENDIF

   ! Return to calling program
 END SUBROUTINE Calc_SumCosZa
!EOC
#if !defined(ESMF_) && !defined( MODEL_WRF )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_met_fields
!
! !DESCRIPTION: Subroutine GET\_MET\_FIELDS calls the various routines to get
! met fields from HEMCO.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Met_Fields( am_I_Root, Input_Opt, State_Met, State_Chm, &
                            Phase,     RC )
!
! ! USES:
!
   USE DAO_Mod
   USE ErrCode_Mod
   USE FlexGrid_Read_Mod
   USE HCO_INTERFACE_MOD,      ONLY : HcoState
   USE HCO_EMISLIST_MOD,       ONLY : HCO_GetPtr 
   USE Input_Opt_Mod,          ONLY : OptInput
   USE Pressure_Mod,           ONLY : Set_Floating_Pressures
   USE State_Chm_Mod,          ONLY : ChmState
   USE State_Met_Mod,          ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   LOGICAL,          INTENT(IN   )          :: am_I_Root  ! root CPU?
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input opts
   INTEGER,          INTENT(IN   )          :: Phase      ! Run phase
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Met state
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry state
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
! 
! !REMARKS:
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  23 Oct 2013 - R. Yantosca - Now pass Input_Opt to GET_A6_FIELDS
!  23 Oct 2013 - R. Yantosca - Now pass Input_Opt to GET_MERRA_A3_FIELDS
!  24 Jun 2014 - R. Yantosca - Now pass Input_Opt to other routines
!  24 Jun 2014 - R. Yantosca - Cosmetic changes, line up arguments
!  12 Aug 2015 - R. Yantosca - Call routines for reading MERRA2 fields
!  25 Oct 2018 - M. Sulprizio- Move READ_INITIAL_MET_FIELDS to hcoi_gc_main_mod
!                              and rename Get_Met_Fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: N_DYN              ! Dynamic timestep in seconds
   INTEGER              :: D(2)               ! Variable for date and time
   LOGICAL              :: FOUND              ! Found in restart file?
   CHARACTER(LEN=255)   :: v_name             ! Variable name 

   ! Pointers
   REAL*4,  POINTER     :: Ptr2D(:,:)
   REAL*4,  POINTER     :: Ptr3D(:,:,:)
   
   !=================================================================
   !    *****  R E A D   M E T   F I E L D S    *****
   !    *****  At the start of the GEOS-Chem simulation  *****
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS
   
   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()
   
   !----------------------------------
   ! Read time-invariant data (Phase 0 only)
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      CALL FlexGrid_Read_CN( Input_Opt, State_Met )
   ENDIF
      
   !----------------------------------
   ! Read 1-hr time-averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A1_TIME()
   ELSE
      D = GET_A1_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A1() ) THEN
      CALL FlexGrid_Read_A1  ( D(1), D(2), Input_Opt, State_Met )
   ENDIF
   
   !----------------------------------
   ! Read 3-hr time averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A3_TIME()
   ELSE
      D = GET_A3_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A3() ) THEN
      CALL FlexGrid_Read_A3  ( D(1), D(2), Input_Opt, State_Met )
   ENDIF

   !----------------------------------
   ! Read 3-hr instantanous data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_I3_TIME()
   ELSE
      D = GET_I3_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_I3() ) THEN
      CALL FlexGrid_Read_I3( D(1), D(2), Input_Opt, State_Met )

      ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
      ! and compute avg dry pressure near polar caps
      CALL Set_Dry_Surface_Pressure( State_Met, 2 )
      CALL AvgPole( State_Met%PS2_DRY )

      ! Compute avg moist pressure near polar caps
      CALL AvgPole( State_Met%PS2_WET ) 

      ! On first call, attempt to get instantaneous met fields for prior
      ! timestep from the GEOS-Chem restart file. Otherwise, initialize
      ! to met fields for this timestep.
      IF ( PHASE == 0 ) THEN

         !-------------
         ! TMPU
         !-------------

         ! Define variable name
         v_name = 'TMPU1'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr3D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Met%TMPU1 = Ptr3D
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Initialize TMPU1    from restart file'
            ENDIF
         ELSE
            State_Met%TMPU1 = State_Met%TMPU2
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'TMPU1    not found in restart, set to TMPU2'
            ENDIF
         ENDIF

         ! Nullify pointer
         Ptr3D => NULL()

         !-------------
         ! SPHU
         !-------------

         ! Define variable name
         v_name = 'SPHU1'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr3D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Met%SPHU1 = Ptr3D
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Initialize SPHU1    from restart file'
            ENDIF
         ELSE
            State_Met%SPHU1 = State_Met%SPHU2
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'SPHU1    not found in restart, set to SPHU2'
            ENDIF
         ENDIF

         ! Nullify pointer
         Ptr3D => NULL()

         !-------------
         ! PS1_WET
         !-------------

         ! Define variable name
         v_name = 'PS1WET'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr2D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Met%PS1_WET = Ptr2D
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Initialize PS1_WET  from restart file'
            ENDIF
         ELSE
            State_Met%PS1_WET = State_Met%PS2_WET
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'PS1_WET  not found in restart, set to PS2_WET'
            ENDIF
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()
         
         !-------------
         ! PS1_DRY
         !-------------

         ! Define variable name
         v_name = 'PS1DRY'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr2D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Met%PS1_DRY = Ptr2D
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Initialize PS1_DRY  from restart file'
            ENDIF
         ELSE
            State_Met%PS1_DRY = State_Met%PS2_DRY
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'PS1_DRY  not found in restart, set to PS2_DRY'
            ENDIF
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()

         !-------------
         ! DELP_DRY
         !-------------

         ! Define variable name
         v_name = 'DELPDRY'

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr3D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN
            State_Met%DELP_DRY = Ptr3D
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'Initialize DELP_DRY from restart file'
            ENDIF
         ELSE
            IF ( am_I_Root ) THEN
               WRITE(6,*) 'DELP_DRY not found in restart, set to zero'
            ENDIF
         ENDIF

         ! Nullify pointer
         Ptr3D => NULL()

         ! Interpolate I-3 fields to current dynamic timestep
         N_DYN = GET_TS_DYN()
         CALL Interp( 0, 0, N_DYN, Input_Opt, State_Met )
         
         ! Initialize surface pressures prior to interpolation
         ! to allow initialization of floating pressures
         State_Met%PSC2_WET = State_Met%PS2_WET
         State_Met%PSC2_DRY = State_Met%PS2_DRY
         CALL Set_Floating_Pressures( am_I_Root, State_Met, RC )

         !=================================================================
         ! Call AIRQNT to compute initial air mass quantities
         !=================================================================
         ! Do not update initial tracer concentrations since not read 
         ! from restart file yet (ewl, 10/28/15)
         CALL AirQnt( am_I_Root, Input_Opt, State_Met, State_Chm, RC, &
                update_mixing_ratio=.FALSE. )

      ENDIF
         
   ENDIF
   
 END SUBROUTINE Get_Met_Fields
!EOC
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
 SUBROUTINE Get_GC_Restart( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!     
   USE CMN_SIZE_Mod
   USE DAO_Mod,            ONLY : AIRQNT
   USE ErrCode_Mod
   USE Error_Mod
   USE HCO_INTERFACE_MOD,  ONLY : HcoState
   USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr 
   USE OCEAN_MERCURY_MOD,  ONLY : CHECK_OCEAN_MERCURY
   USE PHYSCONSTANTS,      ONLY : BOLTZ, AIRMW
   USE Input_Opt_Mod,      ONLY : OptInput
   USE Species_Mod,        ONLY : Species
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Met_Mod,      ONLY : MetState
   USE TIME_MOD,           ONLY : EXPAND_DATE
   USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
   LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
   TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
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
!  20 Apr 2016 - E. Lundgren - Implement ocean and snow Hg variables
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  31 May 2016 - E. Lundgren - Replace Input_Opt%TRACER_MW_G with species
!                              database field emMW_g (emitted species g/mol)
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  22 Jun 2016 - R. Yantosca - Now refer to Hg0_Id_List, Hg2_Id_List, and
!                              HgP_Id_List fields of State_Chm
!  11 Jul 2016 - E. Lundgren - Remove tracers and read only species
!  12 Jul 2016 - E. Lundgren - Rename from read_gc_restart_nc
!  18 Jul 2016 - M. Sulprizio- Remove special handling of ISOPN, MMN, CFCX, and
!                              HCFCX. Family tracers have been eliminated.
!  25 Jul 2016 - E. Lundgren - Store whether species in rst file in species db
!                              rather than module-level variable
!  03 Aug 2016 - E. Lundgren - Remove tracers; now only use species
!  11 Aug 2016 - E. Lundgren - Move source of background values to spc database
!  28 Nov 2017 - R. Yantosca - Now only replace tokens in filename but not
!                              in the rest of the path
!  24 Oct 2018 - M. Sulprizio- Move READ_GC_RESTART to hcoi_gc_main_mod.F90
!                              and rename GET_GC_RESTART
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, M, N      ! lon, lat, lev, indexes
   LOGICAL              :: FOUND              ! Found in restart file?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message 
   CHARACTER(LEN=255)   :: v_name             ! variable name 
   REAL(fp)             :: MW_g               ! species molecular weight
   REAL(fp)             :: SMALL_NUM          ! small number threshold
   CHARACTER(LEN=63)    :: OrigUnit

   ! Temporary arrays and pointers
   REAL*4,  TARGET           :: Temp2D(IIPAR,JJPAR) 
   REAL*4,  TARGET           :: Temp3D(IIPAR,JJPAR,LLPAR)
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
   LOC = ' -> at Get_GC_Restart (in GeosCore/hcoi_gc_main_mod.F)'

   ! Set minimum value threshold for [mol/mol]
   SMALL_NUM = 1.0e-30_fp
      
   !=================================================================
   ! If running Hg simulation, set Hg-specific local variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Set the # of tagHg categories from State_Chm
      Num_Hg_Categories   =  State_Chm%N_Hg_CATS

      ! Set variable storing names for each of the Hg categories
      Hg_Cat_Name => State_Chm%Hg_Cat_Name
    
      ! Set Hg species index corresponding to a given Hg category number;
      ! total is always the first category
      Total_Hg_Id   =  State_Chm%Hg0_Id_List(1)

   ENDIF

   !=================================================================
   ! Open GEOS-Chem restart file
   !=================================================================

   ! Write read message to log
   WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
   WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'

   !=================================================================
   ! Read species concentrations from NetCDF or use default 
   ! background [mol/mol]; store in State_Chm%Species in [kg/kg dry]
   !=================================================================

   ! IMPORTANT NOTE: the unit conversion from mol/mol to kg/kg uses
   ! the molecular weight stored in the species database which is
   ! a meaningful value for advected species but is a bad value (-1)
   ! for all others. Non-advected species should NOT be used when 
   ! State_Chm%Species units are in mass mixing ratio. Current
   ! units can be determined at any point by looking at 
   ! State_Chm%Spc_Units. (ewl, 8/11/16)

   ! Print header for min/max concentration to log
   WRITE( 6, 110 )
110 FORMAT( 'Min and Max of each species in restart file [mol/mol]:' )

   ! Initialize species to all zeroes
   State_Chm%Species = 0.e+0_fp

   ! Loop over species
   DO N = 1, State_Chm%nSpecies

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%emMW_g

      ! Define variable name
      v_name = 'SPC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr3D,     RC,       FOUND=FOUND )

      ! Check if species data is in file
      IF ( FOUND ) THEN
         SpcInfo%Is_InRestart = .TRUE.
      ELSE
         SpcInfo%Is_InRestart = .FALSE.
      ENDIF
      
      ! If data is in file, read in as [mol/mol] and convert to 
      ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
      ! either stored in species database (advected species all levels and
      ! non-advected species levels up to LLCHEM) or a small number
      ! (non-advected species levels above LLCHEM) converted to 
      ! [kg/kg dry]
      IF ( SpcInfo%Is_InRestart ) THEN

         ! Print the min & max of each species as it is read from 
         ! the restart file in mol/mol
         IF ( am_I_Root ) THEN
            WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                            MINVAL( Ptr3D ), MAXVAL( Ptr3D )
120         FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                    '  Max = ',es15.9)
         ENDIF

         ! Convert file value [mol/mol] to [kg/kg dry] for storage
!$OMP PARALLEL DO                                                       &
!$OMP DEFAULT( SHARED )                                                 &
!$OMP PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ! Apply minimum value threshold where input conc is very low
            IF ( Ptr3D(I,J,L) < SMALL_NUM ) THEN
                 Ptr3D(I,J,L) = SMALL_NUM
            ENDIF
            State_Chm%Species(I,J,L,N) = Ptr3D(I,J,L) * MW_g / AIRMW
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         ! Set species to the background value converted to [kg/kg dry] 
!$OMP PARALLEL DO                                                       &
!$OMP DEFAULT( SHARED )                                                 &
!$OMP PRIVATE( I, J, L )
         ! Loop over all grid boxes
         DO L = 1, LLPAR 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
               
            ! Special handling for MOH
            IF ( TRIM( SpcInfo%Name ) == 'MOH' ) THEN

               !----------------------------------------------------
               ! For methanol (MOH), use different initial
               ! background concentrations for different regions of
               ! the atmosphere:
               !
               ! (a) 2.0 ppbv MOH -- continental boundary layer
               ! (b) 0.9 ppbv MOH -- marine boundary layer
               ! (c) 0.6 ppbv MOH -- free troposphere
               !
               ! The concentrations listed above are from Heikes et
               ! al, "Atmospheric methanol budget and ocean
               ! implication", _Global Biogeochem. Cycles_, 2002.
               ! These represent the best estimates for the methanol
               ! conc.'s in the troposphere based on various
               ! measurements.
               !
               ! MOH is an inactive chemical species in GEOS-CHEM,
               ! so these initial concentrations will never change.
               ! However, MOH acts as a sink for OH, and therefore
               ! will affect both the OH concentration and the
               ! methylchloroform lifetime.
               !
               ! We specify the MOH concentration as ppbv, but then
               ! we need to multiply by CONV_FACTOR in order to
               ! convert to [molec/cm3].  (bdf, bmy, 2/22/02)
               !----------------------------------------------------
                  
               ! Test for altitude (L < 9 is always in the trop)
               IF ( L <= 9 ) THEN
                  ! Test for ocean/land boxes
                  IF ( State_Met%FRCLND(I,J) >= 0.5 ) THEN
                     ! Continental boundary layer: 2 ppbv MOH
                     State_Chm%Species(I,J,L,N) = 2.000e-9_fp &
                                                  * MW_g / AIRMW
                  ELSE
                     ! Marine boundary layer: 0.9 ppbv MOH
                     State_Chm%Species(I,J,L,N) = 0.900e-9_fp &
                                                  * MW_g / AIRMW
                  ENDIF
               ELSE
                  ! Test for troposphere
                  IF ( State_Met%InTroposphere(I,J,L) ) THEN
                     ! Free troposphere: 0.6 ppbv MOH
                     State_Chm%Species(I,J,L,N) = 0.600e-9_fp &
                                                  * MW_g / AIRMW
                  ELSE
                     ! Strat/mesosphere:
                     State_Chm%Species(I,J,L,N) =             &
                                        SMALL_NUM * MW_g / AIRMW
                  ENDIF
               ENDIF

               ! Print to log if debugging is on
               IF ( am_I_Root .AND. I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                  WRITE( 6, 130 ) N, TRIM( SpcInfo%Name )
130               FORMAT('Species ', i3, ', ', a9, &
                         ': see READ_GC_RESTART for special MOH values')
               ENDIF

            ! For non-advected species at levels above LLCHEM, use a 
            ! small number for background
            ELSEIF ( L > LLCHEM .AND. &
                   ( .NOT. SpcInfo%Is_Advected ) ) THEN

               State_Chm%Species(I,J,L,N) = SMALL_NUM * MW_g / AIRMW 

            ! For all other cases except MOH, use the background value  
            ! stored in the species database
            ELSE

               State_Chm%Species(I,J,L,N) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

               ! Print to log if debugging is on
               IF ( am_I_Root .AND. I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                  WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
140               FORMAT('Species ', i3, ', ', a9, &
                         ': Use background = ', es15.9)
               ENDIF


            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Set species units
   State_Chm%Spc_Units = 'kg/kg dry'

   ! If in debug mode, print out species min and max in [molec/cm3]
   IF ( am_I_Root .and. Input_Opt%LPRT ) THEN

      ! Convert units
      PRINT *, " "
      PRINT *, "Species min and max in molec/cm3"
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met,  &
                              State_Chm, 'molec/cm3', RC,       &
                              OrigUnit=OrigUnit )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #1!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

      ! Print values
      DO N = 1, State_Chm%nSpecies
         SpcInfo => State_Chm%SpcData(N)%Info
         WRITE(6,150) N, TRIM( SpcInfo%Name ),                 &
                         MINVAL( State_Chm%Species(:,:,:,N) ), & 
                         MAXVAL( State_Chm%Species(:,:,:,N) ) 
150      FORMAT( 'Species ', i3, ', ', a9,                     &
                 ': Min = ', es15.9, ', Max = ', es15.9 )
         SpcInfo => NULL()
      ENDDO

      ! Convert units back
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                              State_Chm, OrigUnit,  RC )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #2!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

   ENDIF

   !=================================================================
   ! Get variables for FlexChem
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LCHEM ) THEN

      ! Define variable name
      v_name = 'KPP_HVALUE'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%KPPHvalue = Ptr3D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize KPP H-value from restart file'
            WRITE(6,160) MINVAL( State_Chm%KPPHvalue(:,:,:) ), & 
                         MAXVAL( State_Chm%KPPHvalue(:,:,:) ) 
160         FORMAT( 'KPP_HVALUE: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%KPPHvalue = 0e+0_fp
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'KPP_HVALUE     not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()
      
   ENDIF

   !=================================================================
   ! Get variables for Soil NOx emissions
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      ! Define variable name
      v_name = 'WETDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr2D,     RC,       FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%WetDepNitrogen = Ptr2D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize wet deposited nitrogen from restart file'
            WRITE(6,170) MINVAL( State_Chm%WetDepNitrogen(:,:) ), & 
                         MAXVAL( State_Chm%WetDepNitrogen(:,:) ) 
170         FORMAT( 12x, '  WETDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%WetDepNitrogen = 0e+0_fp
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'WETDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      ! Define variable name
      v_name = 'DRYDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%DryDepNitrogen = Ptr2D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,180) MINVAL( State_Chm%DryDepNitrogen(:,:) ), & 
                         MAXVAL( State_Chm%DryDepNitrogen(:,:) ) 
180         FORMAT( 12x, '  DRYDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%DryDepNitrogen = 0e+0_fp
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'DRYDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

   ENDIF

   !=================================================================
   ! Read variables for sulfate chemistry
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
        Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      ! Define variable name
      v_name = 'H2O2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%H2O2AfterChem = Ptr3D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize H2O2 from restart file'
            WRITE(6,190) MINVAL( State_Chm%H2O2AfterChem(:,:,:) ), & 
                         MAXVAL( State_Chm%H2O2AfterChem(:,:,:) ) 
190         FORMAT( 12x, 'H2O2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%H2O2AfterChem = 0e+0_fp
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'H2O2_AFTERCHEM not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! Define variable name
      v_name = 'SO2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%SO2AfterChem = Ptr3D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,200) MINVAL( State_Chm%SO2AfterChem(:,:,:) ), & 
                         MAXVAL( State_Chm%SO2AfterChem(:,:,:) ) 
200         FORMAT( 12x, ' SO2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%SO2AfterChem = 0e+0_fp
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'SO2_AFTERCHEM  not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF
   
   !=================================================================
   ! Read variables for UCX
   !=================================================================
   IF ( Input_Opt%LUCX ) THEN

      ! Define variable name
      v_name = 'STATE_PSC'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                       Ptr3D,     RC,       FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%STATE_PSC = Ptr3D
         IF ( am_I_Root ) THEN
            WRITE(6,*) 'Initialize PSC from restart for UCX'
            WRITE(6,210) MINVAL( State_Chm%STATE_PSC(:,:,:) ), & 
                         MAXVAL( State_Chm%STATE_PSC(:,:,:) ) 
210         FORMAT( 12x, ' STATE_PSC: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         IF ( am_I_Root ) THEN
#if defined(ESMF_)
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
            WRITE(6,*) 'STATE_PSC      not found in restart, initialize PSC-free'
         ENDIF
#endif
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF
   
   !=================================================================
   ! Read ocean mercury variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Print total mass to log
      WRITE( 6, 220 )
220   FORMAT(/, 'Total mass of each ocean and snow Hg species:')

      !--------------------------------------------------------------
      ! Total Hg in ocean
      !--------------------------------------------------------------
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
         CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                          Ptr2D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN

            ! Check for negative concentrations (jaf, 7/6/11)
            DO I = 1, IIPAR
            DO J = 1, JJPAR
               IF ( Ptr2D(I,J) < 0.0d4 ) THEN
                  Ptr2D(I,J) = 0.0d4
               ENDIF
            ENDDO
            ENDDO

            ! Assign ocean mercury data and write total mass to log file
            SELECT CASE( M )
            CASE ( 1 )
               State_Chm%OceanHg0(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ), &
                            SUM( State_Chm%OceanHg0(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 2 )
               State_Chm%OceanHg2(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHg2(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 3 )
               State_Chm%OceanHgP(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHgP(:,:,Total_Hg_Id) ), 'kg'
            END SELECT

         ELSE
            WRITE( 6, 230 ) TRIM( v_name )
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()

      ENDDO

      !--------------------------------------------------------------
      ! Additional tagged ocean Hg species
      !--------------------------------------------------------------
      IF ( Input_Opt%LSPLIT ) THEN
         DO M = 1, 3
            DO N = 2, Num_Hg_Categories

               ! Define variable name. Include appended region.
               SELECT CASE( M )
               CASE ( 1 )
                  HgSpc = 'Hg0' 
               CASE ( 2 )
                  HgSpc = 'Hg2'
               CASE ( 3 )
                  HgSpc = 'HgP' 
               END SELECT
               v_name = 'OCEAN_' // TRIM( HgSpc ) //  &
                        '_'      // TRIM( Hg_Cat_Name(N) )

               ! Get variable from HEMCO and store in local array
               CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                                Ptr2D, RC, FOUND=FOUND )

               ! Check if variable is in file
               IF ( FOUND ) THEN

                  ! Assign ocean mercury data and write total mass to log
                  SELECT CASE( M )
                  CASE ( 1 )
                     State_Chm%OceanHg0(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg0(:,:,N) ), 'kg'
                  CASE ( 2 )
                     State_Chm%OceanHg2(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg2(:,:,N) ), 'kg'
                  CASE ( 3 )
                     State_Chm%OceanHgP(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHgP(:,:,N) ), 'kg'
                  END SELECT

               ELSE
                  WRITE( 6, 230 ) TRIM( v_name )
               ENDIF

               ! Nullify pointer
               Ptr2D => NULL()

            ENDDO
         ENDDO

         ! Make sure tagged & total species sum up
         IF ( Input_Opt%USE_CHECKS ) THEN
            CALL CHECK_OCEAN_MERCURY( State_Chm, 'end of READ_GC_RESTART' )
         ENDIF
      ENDIF

      !--------------------------------------------------------------
      ! Hg snowpack on land and ocean
      !--------------------------------------------------------------
      DO M = 1, 4               
         DO N = 1, Num_Hg_Categories               

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

            IF ( N == 1 ) THEN
               v_name = TRIM( Prefix )
            ELSE
               ! Append category name if tagged
               v_name = TRIM( Prefix         ) // '_' // &
                        TRIM( Hg_Cat_Name(N) ) 
            ENDIF

            ! Get variable from HEMCO and store in local array
            CALL HCO_GetPtr( am_I_Root, HcoState, TRIM(v_name), &
                             Ptr2D, RC, FOUND=FOUND )

            ! Check if variable is in file
            IF ( FOUND ) THEN

               ! Assign ocean mercury data and write total mass to file
               SELECT CASE( M )
               CASE ( 1 )
                  State_Chm%SnowHgOcean(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOcean(:,:,N) ), 'kg'
               CASE ( 2 )
                  State_Chm%SnowHgOceanStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOceanStored(:,:,N) ),'kg'
               CASE ( 3 )
                  State_Chm%SnowHgLand(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLand(:,:,N) ), 'kg'
               CASE ( 4 )
                  State_Chm%SnowHgLandStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLandStored(:,:,N) ), 'kg'
               END SELECT

            ELSE
               WRITE( 6, 230 ) TRIM( v_name )
            ENDIF

            ! Nullify pointer
            Ptr2D => NULL()

         ENDDO
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
      
   ! Mark end of section in log
   IF ( Input_Opt%LPRT .AND. am_I_Root ) THEN
      CALL DEBUG_MSG('### DONE GET_GC_RESTART')
   ENDIF
   WRITE( 6, '(a)' ) REPEAT( '=', 79 )

 END SUBROUTINE Get_GC_Restart
#endif
!EOC
END MODULE Hcoi_GC_Main_Mod
