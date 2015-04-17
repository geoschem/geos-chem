!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_state_mod.F90
!
! !DESCRIPTION: Module HCOX\_State\_Mod contains routines and variables
! to organize the extensions state type ExtState. ExtState contains the
! logical switches for each extension (denoting whether or not it is 
! enabled) as well as pointers to all met fields used by the extensions. 
! ExtState is passed to all extension modules, and the met fields
! defined in here are thus available to all extensions. Additional met 
! fields (and extension switches) can be added as required.
!\\
! This module contains the routines to initialize and finalize the
! ExtState object, but doesn't link the met field pointers to the 
! corresponding fields. This is done in the HEMCO-model interface
! routines (e.g. hcoi\_standalone\_mod.F90, hcoi\_gc\_main\_mod.F90).
! Newly added met fields will only work if the corresponding pointer
! assignments are added to these interface routines!
!\\
!\\
! !INTERFACE: 
!
MODULE HCOX_STATE_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_ARR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: ExtStateInit
  PUBLIC :: ExtStateFinal
  PUBLIC :: ExtDat_Set
!
! !DERIVED TYPES:
!
  !=========================================================================
  ! ExtDat_*: Derived types containing pointers to the met field arrays 
  ! (Arr) and a logical flag whether or not the field is used by any of 
  ! the extensions (DoUse).  Arrays can be 3D reals or 2D reals or integer 
  ! All real values are of default precision! (df), as specified in 
  ! HCO\_ERROR\_MOD.  You can add more types if necessary.
  !=========================================================================
 
  ! 2D real, default precision
  TYPE, PUBLIC :: ExtDat_2R
     TYPE(Arr2D_HP), POINTER :: Arr
     LOGICAL                 :: DoUse
     LOGICAL                 :: FromList 
  END TYPE ExtDat_2R

  ! 2D real, single precision
  TYPE, PUBLIC :: ExtDat_2S
     TYPE(Arr2D_SP), POINTER :: Arr
     LOGICAL                 :: DoUse
     LOGICAL                 :: FromList 
  END TYPE ExtDat_2S

  ! 2D integer
  TYPE, PUBLIC :: ExtDat_2I
     TYPE(Arr2D_I),  POINTER :: Arr
     LOGICAL                 :: DoUse
     LOGICAL                 :: FromList 
  END TYPE ExtDat_2I

  ! 3D real, default precision
  TYPE, PUBLIC :: ExtDat_3R
     TYPE(Arr3D_HP), POINTER :: Arr
     LOGICAL                 :: DoUse
     LOGICAL                 :: FromList 
  END TYPE ExtDat_3R

  ! 3D real, single precision
  TYPE, PUBLIC :: ExtDat_3S
     TYPE(Arr3D_SP), POINTER :: Arr
     LOGICAL                 :: DoUse
     LOGICAL                 :: FromList 
  END TYPE ExtDat_3S
  !=========================================================================
  ! Ext_State: Derived type declaration for the State object containing 
  ! pointers to all met fields and related quantities used by the HEMCO 
  ! extensions. An 'Ext_State' type called ExtState is defined at the 
  ! beginning of a HEMCO run and populated according to the specifications
  ! set in the configuration file.  You can add more fields if necessary.
  !=========================================================================
  TYPE, PUBLIC :: Ext_State
 
     !----------------------------------------------------------------------
     ! Extension switches (enabled?)
     ! NOTE: When adding a new extension, don't forget to initialize this
     ! switch in subroutine ExtStateInit below!
     !----------------------------------------------------------------------
     LOGICAL                   :: Custom         ! Customizable ext.
     LOGICAL                   :: DustDead       ! DEAD dust model
     LOGICAL                   :: DustGinoux     ! Ginoux dust emissions
     LOGICAL                   :: LightNOx       ! Lightning NOx
     LOGICAL                   :: ParaNOx        ! PARANOX ship emissions
     LOGICAL                   :: SoilNOx        ! Soil NOx emissions
     LOGICAL                   :: Megan          ! MEGAN biogenic emissions
     LOGICAL                   :: SeaFlux        ! air-sea exchange
     LOGICAL                   :: SeaSalt        ! Seasalt emissions
     LOGICAL                   :: GFED           ! GFED biomass burning
     LOGICAL                   :: FINN           ! FINN biomass burning
     LOGICAL                   :: GC_RnPbBe      ! GEOS-Chem Rn-Pb-Be simulation
     LOGICAL                   :: GC_POPs        ! GEOS-Chem POPs simulation
     LOGICAL                   :: Wetland_CH4    ! Methane emissions from wetlands
     LOGICAL                   :: TOMAS_SeaSalt  ! TOMAS sectional sea salt

     !----------------------------------------------------------------------
     ! Data directory
     !----------------------------------------------------------------------
     CHARACTER(LEN=255)        :: DATA_DIR    ! Directory for data 

     !----------------------------------------------------------------------
     ! Met fields
     !----------------------------------------------------------------------
     TYPE(ExtDat_2R),  POINTER :: U10M        ! E/W 10m wind speed [m/s]
     TYPE(ExtDat_2R),  POINTER :: V10M        ! N/S 10m wind speed [m/s]
     TYPE(ExtDat_2R),  POINTER :: ALBD        ! Surface albedo [-] 
     TYPE(ExtDat_2R),  POINTER :: WLI         ! 0=water, 1=land, 2=ice
     TYPE(ExtDat_2R),  POINTER :: T2M         ! 2m Sfce temperature [K] 
     TYPE(ExtDat_2R),  POINTER :: TSKIN       ! Surface skin temperature [K]
     TYPE(ExtDat_2R),  POINTER :: GWETROOT    ! Root soil wetness [1]
     TYPE(ExtDat_2R),  POINTER :: GWETTOP     ! Top soil moisture [-]
     TYPE(ExtDat_2R),  POINTER :: SNOWHGT     ! Snow height [mm H2O = kg H2O/m2]
     TYPE(ExtDat_2R),  POINTER :: SNODP       ! Snow depth [m ] 
     TYPE(ExtDat_2R),  POINTER :: USTAR       ! Friction velocity [m/s] 
     TYPE(ExtDat_2R),  POINTER :: Z0          ! Sfc roughness height [m]
     TYPE(ExtDat_2R),  POINTER :: TROPP       ! Tropopause pressure [hPa] 
     TYPE(ExtDat_2R),  POINTER :: SUNCOSmid   ! COS (SZA) 
     TYPE(ExtDat_2R),  POINTER :: SZAFACT     ! current SZA/total daily SZA
     TYPE(ExtDat_2R),  POINTER :: PARDR       ! direct photsyn radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: PARDF       ! diffuse photsyn radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: RADSWG      ! surface radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: FRCLND      ! Olson land fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLAND      ! land fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FROCEAN     ! ocean fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLAKE      ! lake fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLANDIC    ! land ice fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: CLDFRC      ! cloud fraction [-]
     TYPE(ExtDat_2R),  POINTER :: JNO2        ! J-Value for NO2 [1/s] 
     TYPE(ExtDat_2R),  POINTER :: JO1D        ! J-Value for O3  [1/s]
     TYPE(ExtDat_2R),  POINTER :: GC_LAI      ! daily leaf area index [cm2/cm2] 
     INTEGER,          POINTER :: PBL_MAX     ! Max height of PBL [level]
     TYPE(ExtDat_3R),  POINTER :: CNV_MFC     ! Convective cloud mass flux [kg/m2/s] 
     TYPE(ExtDat_3R),  POINTER :: FRAC_OF_PBL ! Fraction of grid box in PBL
     TYPE(ExtDat_3R),  POINTER :: SPHU        ! Spec. humidity [g H2O/kg air] 
     TYPE(ExtDat_3R),  POINTER :: TK          ! Air temperature [K]
     TYPE(ExtDat_3R),  POINTER :: AIR         ! Air mass [kg]
     TYPE(ExtDat_3R),  POINTER :: AIRVOL      ! Air volume [m3] 
     TYPE(ExtDat_3R),  POINTER :: O3          ! O3 mass [kg]
     TYPE(ExtDat_3R),  POINTER :: NO          ! NO mass [kg]
     TYPE(ExtDat_3R),  POINTER :: NO2         ! NO2 mass [kg]
     TYPE(ExtDat_3R),  POINTER :: HNO3        ! HNO3 mass [kg]

     !----------------------------------------------------------------------
     ! Deposition parameter
     ! DRY_TOTN and WET_TOTN are the total (dry/wet) deposited N since the
     ! last emission timestep. Even though these numbers are per second,
     ! they may represent accumulated deposition velocities if chemistry
     ! and/or dynamic timestep are not equal to the emission timestep.
     ! These values are used by the soil NOx module. Note that it is assumed
     ! that DRY_TOTN and WET_TOTN are summed over chemistry and transport 
     ! timesteps, respectively!
     !----------------------------------------------------------------------
     TYPE(ExtDat_2R),  POINTER :: DRY_TOTN    ! Dry deposited N   [molec/cm2/s] 
     TYPE(ExtDat_2R),  POINTER :: WET_TOTN    ! Wet deposited N   [kg N/s] 
     REAL(hp),         POINTER :: DRYCOEFF(:) ! Baldocci drydep coeff.
     
     !----------------------------------------------------------------------
     ! Constants for POPs emissions module
     !----------------------------------------------------------------------
     REAL(dp)                  :: POP_DEL_H   ! Delta H [J/mol]
     REAL(dp)                  :: POP_KOA     ! POP octanol-water partition coef
     REAL(dp)                  :: POP_KBC     ! POP BC-air partition coeff.

     !----------------------------------------------------------------------
     ! Fields used in ESMF environment only. These arrays won't be used
     ! in a classic environment. They become filled in HCO_SetExtState_ESMF
     ! in hcoi_esmf_mod.F90 (called from within hcoi_gc_main_mod.F90). 
     ! Note: CNV_TOPP is currently not being used. 
     !----------------------------------------------------------------------
     TYPE(ExtDat_2S),  POINTER :: CNV_TOPP    ! Convective cloud top height 
     TYPE(ExtDat_3S),  POINTER :: RCCODE      ! Convection return code
     TYPE(ExtDat_3S),  POINTER :: BYNCY       ! Buoyancy 

  END TYPE Ext_State
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!  02 Oct 2013 - C. Keller   - Initial version
!  23 Jun 2014 - R. Yantosca - Now add DATA_DIR to Ext_State declaration
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  27 Jun 2014 - C. Keller   - Added FINN biomass burning extension
!  07 Jul 2014 - R. Yantosca - Modified for GEOS-Chem Rn-Pb-Be simulation
!  28 Jul 2014 - C. Keller   - Added J-Values for NO2 and O3 to state obj. 
!  20 Aug 2014 - M. Sulprizio- Modified for GEOS-Chem POPs emissions module
!  01 Oct 2014 - R. Yantosca - Modified for TOMAS sea salt emissions module
!  11 Dec 2014 - M. Yannetti - Updated DRYCOEFF to REAL(hp)
!  10 Mar 2015 - C. Keller   - Fields can now be in HEMCO precision or single
!                              precision. Single precision is useful for 
!                              fields used in ESMF setting. 
!  03 Apr 2015 - C. Keller   - Added ExtDat_Set.
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES: 
!
  INTERFACE ExtDat_Init
     MODULE PROCEDURE ExtDat_Init_2R
     MODULE PROCEDURE ExtDat_Init_2S
     MODULE PROCEDURE ExtDat_Init_2I
     MODULE PROCEDURE ExtDat_Init_3R
     MODULE PROCEDURE ExtDat_Init_3S
  END INTERFACE ExtDat_Init
 
  INTERFACE ExtDat_Set
     MODULE PROCEDURE ExtDat_Set_2R
     MODULE PROCEDURE ExtDat_Set_2S
     MODULE PROCEDURE ExtDat_Set_2I
     MODULE PROCEDURE ExtDat_Set_3R
     MODULE PROCEDURE ExtDat_Set_3S
  END INTERFACE ExtDat_Set
 
  INTERFACE ExtDat_Cleanup
     MODULE PROCEDURE ExtDat_Cleanup_2R
     MODULE PROCEDURE ExtDat_Cleanup_2S
     MODULE PROCEDURE ExtDat_Cleanup_2I
     MODULE PROCEDURE ExtDat_Cleanup_3R
     MODULE PROCEDURE ExtDat_Cleanup_3S
  END INTERFACE ExtDat_Cleanup

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ExtStateInit
!
! !DESCRIPTION: Initializes all fields of the ExtState object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtStateInit( ExtState, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState   ! ExtState object
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REMARKS:
!  You can add more initialization statements as is necessary.
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !======================================================================
    ! ExtStateInit begins here
    !======================================================================

    ! Allocate object 
    IF ( .NOT. ASSOCIATED ( ExtState ) ) ALLOCATE ( ExtState )

    !-----------------------------------------------------------------------
    ! Set all switches to FALSE
    !-----------------------------------------------------------------------
    ExtState%Custom        = .FALSE.
    ExtState%DustDead      = .FALSE.
    ExtState%DustGinoux    = .FALSE.
    ExtState%LightNOx      = .FALSE.
    ExtState%ParaNOx       = .FALSE.
    ExtState%SoilNOx       = .FALSE.
    ExtState%Megan         = .FALSE.
    ExtState%SeaFlux       = .FALSE.
    ExtState%SeaSalt       = .FALSE.
    ExtState%GFED          = .FALSE.
    ExtState%FINN          = .FALSE.
    ExtState%GC_RnPbBe     = .FALSE.
    ExtState%GC_POPs       = .FALSE.
    ExtState%Wetland_CH4   = .FALSE.
    ExtState%TOMAS_SeaSalt = .FALSE.

    !-----------------------------------------------------------------------
    ! Initialize constants for POPs emissions module
    !-----------------------------------------------------------------------
    ExtState%POP_DEL_H   = 0d0
    ExtState%POP_KOA     = 0d0
    ExtState%POP_KBC     = 0d0

    !-----------------------------------------------------------------------
    ! Initialize all met arrays.
    ! This defines a nullified pointer for every met field and sets the
    ! corresponding DoUse flag to FALSE. The pointers to the met fields 
    ! need to be defined in the HEMCO-model interface routine.
    !-----------------------------------------------------------------------
    CALL ExtDat_Init( ExtState%U10M, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%V10M, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%ALBD, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%WLI , RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%T2M, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%TSKIN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%GWETROOT, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%GWETTOP, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SNOWHGT, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SNODP, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%USTAR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%Z0, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%TROPP, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SUNCOSmid, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SZAFACT, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%PARDR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%PARDF, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%RADSWG, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%FRCLND, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%FRLAND, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%FROCEAN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%FRLAKE, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%FRLANDIC, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%CLDFRC, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%GC_LAI, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%JNO2, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%JO1D, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%CNV_MFC, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtState%PBL_MAX    => NULL()

    CALL ExtDat_Init ( ExtState%FRAC_OF_PBL, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SPHU, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%TK, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%AIR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%AIRVOL, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%O3, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%NO, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%NO2, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%HNO3, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%DRY_TOTN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%WET_TOTN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%CNV_TOPP, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%RCCODE, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%BYNCY, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ExtStateInit
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtStateFinal
!
! !DESCRIPTION: Finalizes the ExtState object. This removes all defined 
!  pointer links (i.e. nullifies ExtDat\%Arr), but does not deallocate 
!  the target array!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtStateFinal( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER  :: ExtState
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    !======================================================================
    ! ExtStateFinal begins here
    !======================================================================

    IF ( ASSOCIATED(ExtState) ) THEN

       ! Cleanup arrays. Don't do deepclean, i.e. only nullify pointers!
       CALL ExtDat_Cleanup( ExtState%U10M       )
       CALL ExtDat_Cleanup( ExtState%V10M       )
       CALL ExtDat_Cleanup( ExtState%ALBD       )
       CALL ExtDat_Cleanup( ExtState%WLI        )
       CALL ExtDat_Cleanup( ExtState%T2M        )
       CALL ExtDat_Cleanup( ExtState%TSKIN      )
       CALL ExtDat_Cleanup( ExtState%GWETROOT   )
       CALL ExtDat_Cleanup( ExtState%GWETTOP    )
       CALL ExtDat_Cleanup( ExtState%SNOWHGT    )
       CALL ExtDat_Cleanup( ExtState%SNODP      )
       CALL ExtDat_Cleanup( ExtState%USTAR      )
       CALL ExtDat_Cleanup( ExtState%Z0         )
       CALL ExtDat_Cleanup( ExtState%TROPP      )
       CALL ExtDat_Cleanup( ExtState%SUNCOSmid  )
       CALL ExtDat_Cleanup( ExtState%SZAFACT    )
       CALL ExtDat_Cleanup( ExtState%PARDR      )
       CALL ExtDat_Cleanup( ExtState%PARDF      )
       CALL ExtDat_Cleanup( ExtState%RADSWG     )
       CALL ExtDat_Cleanup( ExtState%FRCLND     )
       CALL ExtDat_Cleanup( ExtState%FRLAND     )
       CALL ExtDat_Cleanup( ExtState%FROCEAN    )
       CALL ExtDat_Cleanup( ExtState%FRLAKE     )
       CALL ExtDat_Cleanup( ExtState%FRLANDIC   )
       CALL ExtDat_Cleanup( ExtState%CLDFRC     )
       CALL ExtDat_Cleanup( ExtState%GC_LAI     )
       CALL ExtDat_Cleanup( ExtState%JNO2       )
       CALL ExtDat_Cleanup( ExtState%JO1D       )
       CALL ExtDat_Cleanup( ExtState%CNV_MFC    )
       CALL ExtDat_Cleanup( ExtState%FRAC_OF_PBL)
       CALL ExtDat_Cleanup( ExtState%SPHU       )
       CALL ExtDat_Cleanup( ExtState%TK         )
       CALL ExtDat_Cleanup( ExtState%AIR        )
       CALL ExtDat_Cleanup( ExtState%AIRVOL     )
       CALL ExtDat_Cleanup( ExtState%O3         )
       CALL ExtDat_Cleanup( ExtState%NO         )
       CALL ExtDat_Cleanup( ExtState%NO2        )
       CALL ExtDat_Cleanup( ExtState%HNO3       )
       CALL ExtDat_Cleanup( ExtState%DRY_TOTN   )
       CALL ExtDat_Cleanup( ExtState%WET_TOTN   )
       CALL ExtDat_Cleanup( ExtState%CNV_TOPP   )
       CALL ExtDat_Cleanup( ExtState%RCCODE     )
       CALL ExtDat_Cleanup( ExtState%BYNCY      )

       ExtState%DRYCOEFF   => NULL()
       ExtState%PBL_MAX    => NULL()

    ENDIF

  END SUBROUTINE ExtStateFinal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Init_2R
!
! !DESCRIPTION: Subroutine ExtDat\_Init\_2R initializes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Init_2R ( ExtDat, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2R), POINTER       :: ExtDat
    INTEGER,         INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ExtDat_Init_2R begins here
    ! ================================================================

    ExtDat     => NULL()
    ALLOCATE(ExtDat)
    ExtDat%Arr => NULL()

    ! Establish pointer to ExtDat%Arr%Val
    CALL HCO_ArrInit( ExtDat%Arr, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtDat%DoUse = .FALSE.
    ExtDat%FromList = .FALSE.

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE ExtDat_Init_2R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Init_2S
!
! !DESCRIPTION: Subroutine ExtDat\_Init\_2S initializes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Init_2S ( ExtDat, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2S), POINTER       :: ExtDat
    INTEGER,         INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ExtDat_Init_2S begins here
    ! ================================================================

    ExtDat     => NULL()
    ALLOCATE(ExtDat)
    ExtDat%Arr => NULL()

    ! Establish pointer to ExtDat%Arr%Val
    CALL HCO_ArrInit( ExtDat%Arr, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtDat%DoUse = .FALSE.
    ExtDat%FromList = .FALSE.

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE ExtDat_Init_2S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Init_2I
!
! !DESCRIPTION: Subroutine ExtDat\_Init\_2I initializes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Init_2I ( ExtDat, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2I), POINTER       :: ExtDat
    INTEGER,         INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

    ! ================================================================
    ! ExtDat_Init_2I begins here
    ! ================================================================

    ExtDat => NULL()
    ALLOCATE(ExtDat)
    ExtDat%Arr => NULL()

    ! Establish pointer to ExtDat%Arr%Val
    CALL HCO_ArrInit( ExtDat%Arr, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtDat%DoUse = .FALSE.
    ExtDat%FromList = .FALSE.

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE ExtDat_Init_2I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Init_3R
!
! !DESCRIPTION: Subroutine ExtDat\_Init\_3R initializes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Init_3R ( ExtDat, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_3R), POINTER       :: ExtDat
    INTEGER,         INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Init_3R begins here
    ! ================================================================

    ExtDat => NULL()
    ALLOCATE(ExtDat)
    ExtDat%Arr => NULL()

    ! Establish pointer to ExtDat%Arr%Val
    CALL HCO_ArrInit( ExtDat%Arr, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtDat%DoUse = .FALSE.
    ExtDat%FromList = .FALSE.

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE ExtDat_Init_3R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Init_3S
!
! !DESCRIPTION: Subroutine ExtDat\_Init\_3S initializes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Init_3S ( ExtDat, RC ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_3S), POINTER       :: ExtDat
    INTEGER,         INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Init_3S begins here
    ! ================================================================

    ExtDat => NULL()
    ALLOCATE(ExtDat)
    ExtDat%Arr => NULL()

    ! Establish pointer to ExtDat%Arr%Val
    CALL HCO_ArrInit( ExtDat%Arr, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ExtDat%DoUse = .FALSE.
    ExtDat%FromList = .FALSE.

    ! Leave
    RC = HCO_SUCCESS

  END SUBROUTINE ExtDat_Init_3S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Cleanup_2R
!
! !DESCRIPTION: Subroutine ExtDat\_Cleanup\_2R removes the given ExtDat type.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Cleanup_2R ( ExtDat ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2R), POINTER       :: ExtDat
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Cleanup_2R begins here
    ! ================================================================

    IF ( ASSOCIATED( ExtDat) ) THEN 
       CALL HCO_ArrCleanup( ExtDat%Arr, DeepClean=.TRUE. ) 
       DEALLOCATE ( ExtDat )
    ENDIF

  END SUBROUTINE ExtDat_Cleanup_2R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Cleanup_2S
!
! !DESCRIPTION: Subroutine ExtDat\_Cleanup\_2S removes the given ExtDat type.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Cleanup_2S ( ExtDat ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2S), POINTER       :: ExtDat
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Cleanup_2S begins here
    ! ================================================================

    IF ( ASSOCIATED( ExtDat) ) THEN 
       CALL HCO_ArrCleanup( ExtDat%Arr, DeepClean=.TRUE. ) 
       DEALLOCATE ( ExtDat )
    ENDIF

  END SUBROUTINE ExtDat_Cleanup_2S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Cleanup_2I
!
! !DESCRIPTION: Subroutine ExtDat\_Cleanup\_2I removes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Cleanup_2I ( ExtDat ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_2I), POINTER       :: ExtDat
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Cleanup_2I begins here
    ! ================================================================

    IF ( ASSOCIATED( ExtDat) ) THEN 
       CALL HCO_ArrCleanup( ExtDat%Arr, DeepClean=.TRUE. ) 
       DEALLOCATE ( ExtDat )
    ENDIF

  END SUBROUTINE ExtDat_Cleanup_2I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Cleanup_3R
!
! !DESCRIPTION: Subroutine ExtDat\_Cleanup\_3R removes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Cleanup_3R( ExtDat ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_3R), POINTER       :: ExtDat
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Cleanup_3R begins here
    ! ================================================================

    IF ( ASSOCIATED( ExtDat) ) THEN 
       CALL HCO_ArrCleanup( ExtDat%Arr, DeepClean=.TRUE. ) 
       DEALLOCATE ( ExtDat )
    ENDIF

  END SUBROUTINE ExtDat_Cleanup_3R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Cleanup_3S
!
! !DESCRIPTION: Subroutine ExtDat\_Cleanup\_3S removes the given ExtDat type. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Cleanup_3S( ExtDat ) 
!
! !INPUT PARAMETERS:
!
    TYPE(ExtDat_3S), POINTER       :: ExtDat
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  23 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ================================================================
    ! ExtDat_Cleanup_3S begins here
    ! ================================================================

    IF ( ASSOCIATED( ExtDat) ) THEN 
       CALL HCO_ArrCleanup( ExtDat%Arr, DeepClean=.TRUE. ) 
       DEALLOCATE ( ExtDat )
    ENDIF

  END SUBROUTINE ExtDat_Cleanup_3S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Set_2R
!
! !DESCRIPTION: Subroutine ExtDat\_Set\_2R sets/updates the data array of an
! ExtDat object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Set_2R ( am_I_Root, HcoState, ExtDat, &
                             FldName,   RC,       First,  Trgt ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2R),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(hp),         INTENT(INOUT), OPTIONAL, TARGET :: Trgt(:,:)
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: NX, NY
    REAL(sp), POINTER  :: Ptr2D(:,:) => NULL() 
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtDat_Set_2R (hcox_state_mod.F90)'
    LOGICAL            :: FRST
    LOGICAL            :: FOUND 

    ! ================================================================
    ! ExtDat_Set_2R begins here
    ! ================================================================

    ! Leave
    RC = HCO_SUCCESS

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Try to get data from list
       CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.
  
             ! Pointer dimensions 
             NX = SIZE(Ptr2D,1)
             NY = SIZE(Ptr2D,2)

             ! Must cover the horizontal grid 
             IF ( ( (NX /= 1) .AND. (NX /= HcoState%NX) ) .OR. &
                  ( (NY /= 1) .AND. (NY /= HcoState%NY) )       ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of input data do not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
   
          ! Target to data
          ELSE
   
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
  
             ! Make sure dimensions agree
             NX = SIZE(Trgt,1)
             NY = SIZE(Trgt,2)
 
             ! Must cover the horizontal grid 
             IF ( (NX/=HcoState%NX) .OR. (NY/=HcoState%NY) ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of target data do not ', &
                   'correspond to simulation grid: ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             ! Link data to target
             ExtDat%Arr%Val => Trgt
   
             ! Make sure it's not from list
             ExtDat%FromList = .FALSE.
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
          ENDIF
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( .NOT. FOUND ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          ! Copy values
          IF ( SIZE(Ptr2D,1) == 1 ) THEN
             ExtDat%Arr%Val(:,:) = Ptr2D(1,1)
          ELSE
             ExtDat%Arr%Val(:,:) = Ptr2D(:,:)
          ENDIF
       ENDIF ! FromList
    ENDIF  

    ! Make sure array exists
    IF ( .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_2R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Set_2S
!
! !DESCRIPTION: Subroutine ExtDat\_Set\_2S sets/updates the data array of an
! ExtDat object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Set_2S ( am_I_Root, HcoState, ExtDat, &
                             FldName,   RC,       First,  Trgt ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2S),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(sp),         INTENT(INOUT), OPTIONAL, TARGET :: Trgt(:,:)
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: NX, NY
    REAL(sp), POINTER  :: Ptr2D(:,:) => NULL() 
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtDat_Set_2S (hcox_state_mod.F90)'
    LOGICAL            :: FRST
    LOGICAL            :: FOUND 

    ! ================================================================
    ! ExtDat_Set_2S begins here
    ! ================================================================

    ! Leave
    RC = HCO_SUCCESS

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Try to get data from list
       CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
 
             ! Pointer dimensions 
             NX = SIZE(Ptr2D,1)
             NY = SIZE(Ptr2D,2)
     
             ! Must cover the horizontal grid 
             IF ( ( (NX /= 1) .AND. (NX /= HcoState%NX) ) .OR. &
                  ( (NY /= 1) .AND. (NY /= HcoState%NY) )       ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of input data do not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Check if input data is of same size. In this case, we can set
             ! a pointer
             IF ( (NX == HcoState%NX) .AND. (NY == HcoState%NY) ) THEN
   
                ExtDat%Arr%Val  => Ptr2D
                ExtDat%FromList = .FALSE.
   
                ! Verbose
                IF ( HCO_IsVerb(2) ) THEN
                   MSG = 'Set extension field pointer to HEMCO data list field ' // TRIM(FldName)
                   CALL HCO_MSG(MSG)
                ENDIF
   
             ELSE
                ExtDat%FromList = .TRUE.
   
                ! Make sure array is allocated
                CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
   
                ! Verbose
                IF ( HCO_IsVerb(2) ) THEN
                   MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                   CALL HCO_MSG(MSG)
                ENDIF
             ENDIF
   
          ! Target to data
          ELSE
   
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
  
             ! Make sure dimensions agree
             NX = SIZE(Trgt,1)
             NY = SIZE(Trgt,2)
 
             ! Must cover the horizontal grid 
             IF ( (NX/=HcoState%NX) .OR. (NY/=HcoState%NY) ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of target data do not ', &
                   'correspond to simulation grid: ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Link data to target
             ExtDat%Arr%Val => Trgt
   
             ! Make sure it's not from list
             ExtDat%FromList = .FALSE.
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( .NOT. FOUND ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          ! Copy values
          IF ( SIZE(Ptr2D,1) == 1 ) THEN
             ExtDat%Arr%Val(:,:) = Ptr2D(1,1)
          ELSE
             ExtDat%Arr%Val(:,:) = Ptr2D(:,:)
          ENDIF
       ENDIF !FromList
    ENDIF

    ! Make sure array exists
    IF ( .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_2S
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Set_2I
!
! !DESCRIPTION: Subroutine ExtDat\_Set\_2I sets/updates the data array of an
! ExtDat object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Set_2I ( am_I_Root, HcoState, ExtDat, &
                             FldName,   RC,       First,  Trgt ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2I),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    INTEGER,          INTENT(INOUT), OPTIONAL, TARGET :: Trgt(:,:)
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: NX, NY
    REAL(sp), POINTER  :: Ptr2D(:,:) => NULL() 
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtDat_Set_2I (hcox_state_mod.F90)'
    LOGICAL            :: FRST
    LOGICAL            :: FOUND 

    ! ================================================================
    ! ExtDat_Set_2I begins here
    ! ================================================================

    ! Leave
    RC = HCO_SUCCESS

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Try to get data from list
       CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.
  
             ! Pointer dimensions 
             NX = SIZE(Ptr2D,1)
             NY = SIZE(Ptr2D,2)

             ! Must cover the horizontal grid 
             IF ( ( (NX /= 1) .AND. (NX /= HcoState%NX) ) .OR. &
                  ( (NY /= 1) .AND. (NY /= HcoState%NY) )       ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of input data do not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
   
          ! Target to data
          ELSE
   
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
   
             ! Make sure dimensions agree
             NX = SIZE(Trgt,1)
             NY = SIZE(Trgt,2)
 
             ! Must cover the horizontal grid 
             IF ( (NX /= HcoState%NX) .OR. (NY /= HcoState%NY) ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of target data do not ', &
                   'correspond to simulation grid: ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Link data to target
             ExtDat%Arr%Val => Trgt
   
             ! Make sure it's not from list
             ExtDat%FromList = .FALSE.
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( .NOT. FOUND ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          ! Copy values
          IF ( SIZE(Ptr2D,1) == 1 ) THEN
             ExtDat%Arr%Val(:,:) = Ptr2D(1,1)
          ELSE
             ExtDat%Arr%Val(:,:) = Ptr2D(:,:)
          ENDIF
       ENDIF !FromList
    ENDIF 
   
    ! Make sure array exists
    IF ( .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_2I
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Set_3R
!
! !DESCRIPTION: Subroutine ExtDat\_Set\_3R sets/updates the data array of an
! ExtDat object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Set_3R ( am_I_Root, HcoState, ExtDat, FldName, & 
                             RC,        First,    Trgt,   OnLevEdge ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_3R),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(hp),         INTENT(INOUT), OPTIONAL, TARGET :: Trgt(:,:,:)
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: OnLevEdge 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: NX, NY, NZ, NZ_EXPECTED
    INTEGER            :: L
    LOGICAL            :: FRST
    LOGICAL            :: FOUND 
    REAL(sp), POINTER  :: Ptr2D(:,:)   => NULL() 
    REAL(sp), POINTER  :: Ptr3D(:,:,:) => NULL() 
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtDat_Set_3R (hcox_state_mod.F90)'

    ! ================================================================
    ! ExtDat_Set_3R begins here
    ! ================================================================

    ! Leave
    RC = HCO_SUCCESS

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! Expected number of vertical levels: NZ if not on edge, NZ+1 if on edge
    NZ_EXPECTED = HcoState%NZ
    IF ( PRESENT(OnLevEdge) ) THEN
       IF ( OnLevEdge ) THEN
          NZ_EXPECTED = HcoState%NZ + 1
       ENDIF
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Try to get data from list
       CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr3D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! Also try to read 2D if not found 
       IF ( .NOT. FOUND ) THEN
          CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr2D, RC, FOUND=FOUND )
          IF ( RC /= HCO_SUCCESS ) RETURN     
       ENDIF 

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.

             ! Dimensions of data pointer
             IF ( ASSOCIATED(Ptr3D) ) THEN
                NX = SIZE(Ptr3D,1)
                NY = SIZE(Ptr3D,2)
                NZ = SIZE(Ptr3D,3)
             ELSEIF ( ASSOCIATED(Ptr2D) ) THEN 
                NX = SIZE(Ptr2D,1)
                NY = SIZE(Ptr2D,2)
                NZ = 1
             ! will cause error:
             ELSE
                NX = 0
                NY = 0
                NZ = 0
             ENDIF

             ! Must cover the horizontal grid 
             IF ( ( (NX /= 1) .AND. (NX /= HcoState%NX) ) .OR. &
                  ( (NY /= 1) .AND. (NY /= HcoState%NY) )       ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of input data do not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
   
             ! Must also cover vertical extension
             IF ( (NZ /= 1) .AND. (NZ /= NZ_EXPECTED) ) THEN
                WRITE(MSG,*) 'Vertical dimension of input data does not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimension: ', NZ_EXPECTED, &
                   '; encountered dimensions: ', NZ, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ_EXPECTED, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
   
          ! Target to data
          ELSE
   
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
  
             ! Make sure dimensions agree
             NX = SIZE(Trgt,1)
             NY = SIZE(Trgt,2)
             NZ = SIZE(Trgt,3)
 
             ! Must cover the horizontal grid 
             IF ( (NX/=HcoState%NX) .OR. (NY/=HcoState%NY) .OR. (NZ/=NZ_EXPECTED) ) THEN
                WRITE(MSG,*) 'Dimensions of target data do not ', &
                   'correspond to simulation grid: ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, NZ_EXPECTED, &
                   '; encountered dimensions: ', NX, NY, NZ, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Link data to target
             ExtDat%Arr%Val => Trgt
   
             ! Make sure it's not from list
             ExtDat%FromList = .FALSE.
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( .NOT. FOUND ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          ! If it's a 3D array...
          IF ( ASSOCIATED(Ptr3D) ) THEN
             ExtDat%Arr%Val(:,:,:) = Ptr3D(:,:,:)
   
          ! If it's a 2D array...
          ELSEIF ( ASSOCIATED(Ptr2D) ) THEN
             IF ( SIZE(Ptr2D,1) == 1 ) THEN
                ExtDat%Arr%Val(:,:,:) = Ptr2D(1,1)
             ELSE
                DO L = 1, NZ_EXPECTED 
                   ExtDat%Arr%Val(:,:,L) = Ptr2D(:,:)
                ENDDO
             ENDIF
          ENDIF
       ENDIF !FromList
    ENDIF 

    ! Make sure array exists
    IF ( .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_3R
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtDat_Set_3S
!
! !DESCRIPTION: Subroutine ExtDat\_Set\_3S sets/updates the data array of an
! ExtDat object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtDat_Set_3S ( am_I_Root, HcoState, ExtDat, FldName, &
                             RC,        First,    Trgt,   OnLevEdge ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_3S),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(sp),         INTENT(INOUT), OPTIONAL, TARGET :: Trgt(:,:,:)
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: OnLevEdge 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: NX, NY, NZ, NZ_EXPECTED
    INTEGER            :: L
    LOGICAL            :: FRST
    LOGICAL            :: FOUND 
    REAL(sp), POINTER  :: Ptr2D(:,:)   => NULL() 
    REAL(sp), POINTER  :: Ptr3D(:,:,:) => NULL() 
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtDat_Set_3S (hcox_state_mod.F90)'

    ! ================================================================
    ! ExtDat_Set_3S begins here
    ! ================================================================

    ! Leave
    RC = HCO_SUCCESS

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! Expected number of vertical levels: NZ if not on edge, NZ+1 if on edge
    NZ_EXPECTED = HcoState%NZ
    IF ( PRESENT(OnLevEdge) ) THEN
       IF ( OnLevEdge ) THEN
          NZ_EXPECTED = HcoState%NZ + 1
       ENDIF
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Try to get data from list
       CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr3D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! Also try to read 2D if not found 
       IF ( .NOT. FOUND ) THEN
          CALL HCO_GetPtr( am_I_Root, TRIM(FldName), Ptr2D, RC, FOUND=FOUND )
          IF ( RC /= HCO_SUCCESS ) RETURN     
       ENDIF 

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
  
             ! Dimensions of data pointer
             IF ( ASSOCIATED(Ptr3D) ) THEN
                NX = SIZE(Ptr3D,1)
                NY = SIZE(Ptr3D,2)
                NZ = SIZE(Ptr3D,3)
             ELSEIF ( ASSOCIATED(Ptr2D) ) THEN 
                NX = SIZE(Ptr2D,1)
                NY = SIZE(Ptr2D,2)
                NZ = 1
             ! will cause error:
             ELSE
                NX = 0
                NY = 0
                NZ = 0
             ENDIF

             ! Must cover the horizontal grid 
             IF ( ( (NX /= 1) .AND. (NX /= HcoState%NX) ) .OR. &
                  ( (NY /= 1) .AND. (NY /= HcoState%NY) )       ) THEN
                WRITE(MSG,*) 'Horizontal dimensions of input data do not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, &
                   '; encountered dimensions: ', NX, NY, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
   
             ! Must also cover vertical extension
             IF ( (NZ /= 1) .AND. (NZ /= NZ_EXPECTED) ) THEN
                WRITE(MSG,*) 'Vertical dimension of input data does not ', &
                   'correspond to simulation grid (and is not 1, either): ', &
                   'Expected dimension: ', NZ_EXPECTED, &
                   '; encountered dimensions: ', NZ, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Set pointer to 3D field
             IF ( ASSOCIATED(Ptr3D) .AND. (NX==HcoState%NX) ) THEN
                ExtDat%Arr%Val  => Ptr3D
                ExtDat%FromList = .FALSE.
   
                ! Verbose
                IF ( HCO_IsVerb(2) ) THEN
                   MSG = 'Set extension field pointer to HEMCO data list field ' // TRIM(FldName)
                   CALL HCO_MSG(MSG)
                ENDIF
   
             ELSE
                ExtDat%FromList = .TRUE.
   
                ! Make sure array is allocated
                CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ_EXPECTED, RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
   
                ! Verbose
                IF ( HCO_IsVerb(2) ) THEN
                   MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                   CALL HCO_MSG(MSG)
                ENDIF
             ENDIF
   
          ! Target to data
          ELSE
   
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
  
             ! Make sure dimensions agree
             NX = SIZE(Trgt,1)
             NY = SIZE(Trgt,2)
             NZ = SIZE(Trgt,3)
 
             ! Must cover the horizontal grid 
             IF ( (NX/=HcoState%NX) .OR. (NY/=HcoState%NY) .OR. (NZ/=NZ_EXPECTED) ) THEN
                WRITE(MSG,*) 'Dimensions of target data do not ', &
                   'correspond to simulation grid: ', &
                   'Expected dimensions: ', HcoState%NX, HcoState%NY, NZ_EXPECTED, &
                   '; encountered dimensions: ', NX, NY, NZ, '. Error occured ', &
                   'for field ', TRIM(FldName)
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
 
             ! Link data to target
             ExtDat%Arr%Val => Trgt
   
             ! Make sure it's not from list
             ExtDat%FromList = .FALSE.
   
             ! Verbose
             IF ( HCO_IsVerb(2) ) THEN
                MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                CALL HCO_MSG(MSG)
             ENDIF
          ENDIF 
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( .NOT. FOUND ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
          ! If it's a 3D array...
          IF ( ASSOCIATED(Ptr3D) ) THEN
             ExtDat%Arr%Val(:,:,:) = Ptr3D(:,:,:)
   
          ! If it's a 2D array...
          ELSEIF ( ASSOCIATED(Ptr2D) ) THEN
             IF ( SIZE(Ptr2D,1) == 1 ) THEN
                ExtDat%Arr%Val(:,:,:) = Ptr2D(1,1)
             ELSE
                DO L = 1, NZ_EXPECTED
                   ExtDat%Arr%Val(:,:,L) = Ptr2D(:,:)
                ENDDO
             ENDIF
          ENDIF
       ENDIF ! FromList
    ENDIF  
 
    ! Make sure array exists
    IF ( .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_3S
!EOC
END MODULE HCOX_STATE_MOD
