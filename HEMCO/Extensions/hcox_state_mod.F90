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
     INTEGER                   :: Custom         ! Customizable ext.
     INTEGER                   :: DustDead       ! DEAD dust model
     INTEGER                   :: DustGinoux     ! Ginoux dust emissions
     INTEGER                   :: DustAlk        ! Dust alkalinity
     INTEGER                   :: LightNOx       ! Lightning NOx
     INTEGER                   :: ParaNOx        ! PARANOX ship emissions
     INTEGER                   :: SoilNOx        ! Soil NOx emissions
     INTEGER                   :: Megan          ! MEGAN biogenic emissions
     INTEGER                   :: SeaFlux        ! air-sea exchange
     INTEGER                   :: SeaSalt        ! Seasalt emissions
     INTEGER                   :: MarinePOA      ! Marine organic aerosols
     INTEGER                   :: GFED           ! GFED biomass burning
     INTEGER                   :: FINN           ! FINN biomass burning
     INTEGER                   :: GC_RnPbBe      ! GEOS-Chem Rn-Pb-Be simulation
     INTEGER                   :: GC_POPs        ! GEOS-Chem POPs simulation
     INTEGER                   :: Wetland_CH4    ! Methane emiss from wetlands
     INTEGER                   :: TOMAS_Jeagle   ! TOMAS Jeagle sea salt
     INTEGER                   :: TOMAS_DustDead ! TOMAS sectional Dead Dust
     INTEGER                   :: AeroCom        ! AeroCom volcano 
     INTEGER                   :: Inorg_Iodine   ! Oceanic inorganic iodine emissions

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
     TYPE(ExtDat_2R),  POINTER :: SNICE       ! Fraction of snow/ice [1]
     TYPE(ExtDat_2R),  POINTER :: USTAR       ! Friction velocity [m/s] 
     TYPE(ExtDat_2R),  POINTER :: Z0          ! Sfc roughness height [m]
     TYPE(ExtDat_2R),  POINTER :: TROPP       ! Tropopause pressure [Pa] 
     TYPE(ExtDat_2R),  POINTER :: SUNCOS      ! COS (SZA) 
     TYPE(ExtDat_2R),  POINTER :: SZAFACT     ! current SZA/total daily SZA
     TYPE(ExtDat_2R),  POINTER :: PARDR       ! direct photsyn radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: PARDF       ! diffuse photsyn radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: PSC2_WET    ! Interpolated sfc pressure [hPa]
     TYPE(ExtDat_2R),  POINTER :: RADSWG      ! surface radiation [W/m2]
     TYPE(ExtDat_2R),  POINTER :: FRCLND      ! Olson land fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLAND      ! land fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FROCEAN     ! ocean fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLAKE      ! lake fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: FRLANDIC    ! land ice fraction [-] 
     TYPE(ExtDat_2R),  POINTER :: CLDFRC      ! cloud fraction [-]
     TYPE(ExtDat_2R),  POINTER :: JNO2        ! J-Value for NO2 [1/s] 
     TYPE(ExtDat_2R),  POINTER :: JOH         ! J-Value for O3->OH  [1/s]
     TYPE(ExtDat_2R),  POINTER :: LAI         ! daily leaf area index [cm2/cm2]
     TYPE(ExtDat_2R),  POINTER :: CHLR        ! daily chlorophyll-a [mg/m3]
     INTEGER,          POINTER :: PBL_MAX     ! Max height of PBL [level]
     TYPE(ExtDat_3R),  POINTER :: CNV_MFC     ! Convective cloud mass flux [kg/m2/s] 
     TYPE(ExtDat_3R),  POINTER :: FRAC_OF_PBL ! Fraction of grid box in PBL
     TYPE(ExtDat_3R),  POINTER :: SPHU        ! Spec. humidity [kg H2O/kg total air] 
     TYPE(ExtDat_3R),  POINTER :: TK          ! Air temperature [K]
     TYPE(ExtDat_3R),  POINTER :: AIR         ! Dry air mass [kg]
     TYPE(ExtDat_3R),  POINTER :: AIRVOL      ! Air volume [m3] 
     TYPE(ExtDat_3R),  POINTER :: AIRDEN      ! Dry air density [kg/m3] 
     TYPE(ExtDat_3R),  POINTER :: O3          ! O3 mass [kg/kg dry air]
     TYPE(ExtDat_3R),  POINTER :: NO          ! NO mass [kg/kg dry air]
     TYPE(ExtDat_3R),  POINTER :: NO2         ! NO2 mass [kg/kg dry air]
     TYPE(ExtDat_3R),  POINTER :: HNO3        ! HNO3 mass [kg/kg dry air]
     TYPE(ExtDat_3R),  POINTER :: POPG        ! POPG mass [kg/kg dry air]

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
     REAL(dp)                  :: POP_DEL_Hw  ! Delta Hw [J/mol]
     REAL(dp)                  :: POP_HSTAR   ! Henry's law constant [atm/M/L]
     REAL(dp)                  :: POP_KOA     ! POP octanol-water partition coef
     REAL(dp)                  :: POP_KBC     ! POP BC-air partition coeff.
     REAL(dp)                  :: POP_XMW     ! POP molecular weight [kg/mol]

     !----------------------------------------------------------------------
     ! Fields used in ESMF environment only. These arrays won't be used
     ! in a classic environment. They become filled in HCO_SetExtState_ESMF
     ! in hcoi_esmf_mod.F90 (called from within hcoi_gc_main_mod.F90). 
     !----------------------------------------------------------------------
     TYPE(ExtDat_3S),  POINTER :: BYNCY       ! Buoyancy 
     TYPE(ExtDat_2S),  POINTER :: LFR         ! Lightning flash rate 
     TYPE(ExtDat_2R),  POINTER :: CNV_FRC     ! convective fraction (filled
                                              ! from State_Met) 
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
!  21 Feb 2016 - C. Keller   - Update to HEMCO v2.0
!  03 Mar 2016 - C. Keller   - Added CNV_FRC
!  20 Apr 2016 - M. Sulprizio- Change JO1D pointer to JOH to reflect that it now
!                              points to the effective O3 + hv -> 2OH rates
!  01 Nov 2016 - M. Sulprizio- Rename TOMAS sea salt to TOMAS Jeagle (J. Kodros)
!  17 Oct 2017 - C. Keller   - Add lightning flash rate 
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
    ! Set all switches to -1
    !-----------------------------------------------------------------------
    ExtState%Custom         = -1
    ExtState%DustDead       = -1
    ExtState%DustGinoux     = -1
    ExtState%DustAlk        = -1
    ExtState%LightNOx       = -1
    ExtState%ParaNOx        = -1
    ExtState%SoilNOx        = -1
    ExtState%Megan          = -1
    ExtState%SeaFlux        = -1
    ExtState%SeaSalt        = -1
    ExtState%MarinePOA      = -1
    ExtState%GFED           = -1
    ExtState%FINN           = -1
    ExtState%GC_RnPbBe      = -1
    ExtState%GC_POPs        = -1
    ExtState%Wetland_CH4    = -1 
    ExtState%TOMAS_Jeagle   = -1
    ExtState%TOMAS_DustDead = -1
    ExtState%AeroCom        = -1
    ExtState%Inorg_Iodine   = -1

    !-----------------------------------------------------------------------
    ! Initialize constants for POPs emissions module
    !-----------------------------------------------------------------------
    ExtState%POP_DEL_H      = 0d0
    ExtState%POP_DEL_Hw     = 0d0
    ExtState%POP_HSTAR      = 0d0
    ExtState%POP_KOA        = 0d0
    ExtState%POP_KBC        = 0d0
    ExtState%POP_XMW        = 0d0

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

    CALL ExtDat_Init ( ExtState%SNICE, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%USTAR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%Z0, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%TROPP, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SUNCOS, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%SZAFACT, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%PARDR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%PARDF, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%PSC2_WET, RC ) 
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

    CALL ExtDat_Init ( ExtState%LAI, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%CHLR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%JNO2, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%JOH, RC ) 
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

    CALL ExtDat_Init ( ExtState%AIRDEN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%O3, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%NO, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%NO2, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%HNO3, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%POPG, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%DRY_TOTN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%WET_TOTN, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%BYNCY, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%LFR, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Init ( ExtState%CNV_FRC, RC ) 
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
!  09 Jul 2015 - E. Lundgren - Add chlorophyll-a (CHLR)
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
       CALL ExtDat_Cleanup( ExtState%SNICE      )
       CALL ExtDat_Cleanup( ExtState%USTAR      )
       CALL ExtDat_Cleanup( ExtState%Z0         )
       CALL ExtDat_Cleanup( ExtState%TROPP      )
       CALL ExtDat_Cleanup( ExtState%SUNCOS     )
       CALL ExtDat_Cleanup( ExtState%SZAFACT    )
       CALL ExtDat_Cleanup( ExtState%PARDR      )
       CALL ExtDat_Cleanup( ExtState%PARDF      )
       CALL ExtDat_Cleanup( ExtState%PSC2_WET   )
       CALL ExtDat_Cleanup( ExtState%RADSWG     )
       CALL ExtDat_Cleanup( ExtState%FRCLND     )
       CALL ExtDat_Cleanup( ExtState%FRLAND     )
       CALL ExtDat_Cleanup( ExtState%FROCEAN    )
       CALL ExtDat_Cleanup( ExtState%FRLAKE     )
       CALL ExtDat_Cleanup( ExtState%FRLANDIC   )
       CALL ExtDat_Cleanup( ExtState%CLDFRC     )
       CALL ExtDat_Cleanup( ExtState%LAI        )
       CALL ExtDat_Cleanup( ExtState%CHLR       )
       CALL ExtDat_Cleanup( ExtState%JNO2       )
       CALL ExtDat_Cleanup( ExtState%JOH        )
       CALL ExtDat_Cleanup( ExtState%CNV_MFC    )
       CALL ExtDat_Cleanup( ExtState%FRAC_OF_PBL)
       CALL ExtDat_Cleanup( ExtState%SPHU       )
       CALL ExtDat_Cleanup( ExtState%TK         )
       CALL ExtDat_Cleanup( ExtState%AIR        )
       CALL ExtDat_Cleanup( ExtState%AIRVOL     )
       CALL ExtDat_Cleanup( ExtState%AIRDEN     )
       CALL ExtDat_Cleanup( ExtState%O3         )
       CALL ExtDat_Cleanup( ExtState%NO         )
       CALL ExtDat_Cleanup( ExtState%NO2        )
       CALL ExtDat_Cleanup( ExtState%HNO3       )
       CALL ExtDat_Cleanup( ExtState%POPG       )
       CALL ExtDat_Cleanup( ExtState%DRY_TOTN   )
       CALL ExtDat_Cleanup( ExtState%WET_TOTN   )
       CALL ExtDat_Cleanup( ExtState%CNV_FRC    )
       CALL ExtDat_Cleanup( ExtState%BYNCY      )
       CALL ExtDat_Cleanup( ExtState%LFR        )

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
  SUBROUTINE ExtDat_Set_2R ( am_I_Root, HcoState, ExtDat,  &
                             FldName,   RC,       First,   &
                             Trgt,      Filled,   NotFillOk ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_CALC_MOD,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2R),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(hp),         POINTER      , OPTIONAL         :: Trgt(:,:)
    LOGICAL,          INTENT(  OUT), OPTIONAL         :: Filled
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: NotFillOk 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!  11 May 2015 - C. Keller - Now use HCO_EvalFld instead of HCO_GetPtr. This
!                            allows the application of scale factors to
!                            ExtState fields read through the HEMCO interface. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NX, NY
    REAL(hp), ALLOCATABLE  :: Arr2D(:,:)
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'ExtDat_Set_2R (hcox_state_mod.F90)'
    LOGICAL                :: FRST
    LOGICAL                :: FOUND 
    LOGICAL                :: FailIfNotFilled 

    ! ================================================================
    ! ExtDat_Set_2R begins here
    ! ================================================================

    ! Initialize 
    RC = HCO_SUCCESS
    IF ( PRESENT(Filled) ) Filled = .FALSE.

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! Check for fill requirement
    IF ( PRESENT(NotFillOk) ) THEN
       FailIfNotFilled = .NOT. NotFillOk
    ELSE
       FailIfNotFilled = .TRUE.
    ENDIF

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Allocate temporary array
       ALLOCATE(Arr2D(HcoState%NX,HcoState%NY),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, "Arr2D allocation error", RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Try to get data from list
       CALL HCO_EvalFld( am_I_Root, HcoState, TRIM(FldName), Arr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.
  
             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
   
          ! Target to data
          ELSEIF ( PRESENT(Trgt) ) THEN
  
             ! If target is not associated: 
             IF ( .NOT. ASSOCIATED(Trgt) ) THEN
                IF ( FailIfNotFilled ) THEN
                   MSG = 'Cannot fill extension field ' // TRIM(FldName) // &
                         ' because target field is not associated.'
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
          
             ! If target is associated:
             ELSE
 
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
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
   
                ! Link data to target
                ExtDat%Arr%Val => Trgt
      
                ! Make sure it's not from list
                ExtDat%FromList = .FALSE.
      
                ! This array is now filled
                IF ( PRESENT(Filled) ) Filled = .TRUE.
   
                ! Verbose
                IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                   MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                   CALL HCO_MSG(HcoState%Config%Err,MSG)
                ENDIF
             ENDIF

          ! Field not found and no target defined 
          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot fill extension field ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( FOUND ) THEN
             ! Copy values and mark array as filled
             ExtDat%Arr%Val(:,:) = Arr2D(:,:)
             IF ( PRESENT(Filled) ) Filled = .TRUE.
          ELSEIF ( FailIfNotFilled ) Then
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF ! FromList
    ENDIF  

    ! Make sure array exists
    IF ( FailIfNotFilled .AND. .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
    ENDIF

    ! Cleanup
    IF ( ALLOCATED(Arr2D) ) DEALLOCATE(Arr2D)
 
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
  SUBROUTINE ExtDat_Set_2S ( am_I_Root, HcoState, ExtDat,  &
                             FldName,   RC,       First,   &
                             Trgt,      Filled,   NotFillOk ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_CALC_MOD,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2S),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(sp),         POINTER      , OPTIONAL         :: Trgt(:,:)
    LOGICAL,          INTENT(  OUT), OPTIONAL         :: Filled
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: NotFillOk 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!  11 May 2015 - C. Keller - Now use HCO_EvalFld instead of HCO_GetPtr. This
!                            allows the application of scale factors to
!                            ExtState fields read through the HEMCO interface. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NX, NY
    REAL(hp), ALLOCATABLE  :: Arr2D(:,:)
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'ExtDat_Set_2S (hcox_state_mod.F90)'
    LOGICAL                :: FRST
    LOGICAL                :: FOUND 
    LOGICAL                :: FailIfNotFilled 

    ! ================================================================
    ! ExtDat_Set_2S begins here
    ! ================================================================

    ! Init 
    RC = HCO_SUCCESS
    IF ( PRESENT(Filled) ) Filled = .FALSE.

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! Check for fill requirement
    IF ( PRESENT(NotFillOk) ) THEN
       FailIfNotFilled = .NOT. NotFillOk
    ELSE
       FailIfNotFilled = .TRUE.
    ENDIF

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Allocate temporary array
       ALLOCATE(Arr2D(HcoState%NX,HcoState%NY),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, "Arr2D allocation error", RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Try to get data from list
       CALL HCO_EvalFld( am_I_Root, HcoState, TRIM(FldName), Arr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.
  
             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
   
          ! Target to data
          ELSEIF ( PRESENT(Trgt) ) THEN
 
             ! If target is not associated: 
             IF ( .NOT. ASSOCIATED(Trgt) ) THEN
                IF ( FailIfNotFilled ) THEN
                   MSG = 'Cannot fill extension field ' // TRIM(FldName) // &
                         ' because target field is not associated.'
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
          
             ! If target is associated:
             ELSE
 
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
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
   
                ! Link data to target
                ExtDat%Arr%Val => Trgt
      
                ! Make sure it's not from list
                ExtDat%FromList = .FALSE.
     
                ! Mark as filled 
                IF ( PRESENT(Filled) ) Filled = .TRUE.
   
                ! Verbose
                IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                   MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                   CALL HCO_MSG(HcoState%Config%Err,MSG)
                ENDIF
             ENDIF

          ! Field not found and no target defined 
          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot fill extension field ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( FOUND ) THEN
             ! Copy values and mark as filled
             ExtDat%Arr%Val(:,:) = Arr2D(:,:)
             IF ( PRESENT(Filled) ) Filled = .TRUE.
          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF 
       ENDIF ! FromList
    ENDIF  

    ! Make sure array exists
    IF ( FailIfNotFilled .AND. .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
    ENDIF

    ! Cleanup
    IF ( ALLOCATED(Arr2D) ) DEALLOCATE(Arr2D)
 
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
                             FldName,   RC,       First,  &
                             Trgt,      Filled,   NotFillOk ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_CALC_MOD,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_2I),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    INTEGER,          POINTER,       OPTIONAL         :: Trgt(:,:)
    LOGICAL,          INTENT(  OUT), OPTIONAL         :: Filled
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: NotFillOk 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!  11 May 2015 - C. Keller - Now use HCO_EvalFld instead of HCO_GetPtr. This
!                            allows the application of scale factors to
!                            ExtState fields read through the HEMCO interface. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NX, NY
    REAL(hp), ALLOCATABLE  :: Arr2D(:,:)
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'ExtDat_Set_2I (hcox_state_mod.F90)'
    LOGICAL                :: FRST
    LOGICAL                :: FOUND 
    LOGICAL                :: FailIfNotFilled 

    ! ================================================================
    ! ExtDat_Set_2I begins here
    ! ================================================================

    ! Init 
    RC = HCO_SUCCESS
    IF ( PRESENT(Filled) ) Filled = .FALSE.

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! Check for fill requirement
    IF ( PRESENT(NotFillOk) ) THEN
       FailIfNotFilled = .NOT. NotFillOk
    ELSE
       FailIfNotFilled = .TRUE.
    ENDIF

    ! On first call or if data is flagged as being read from list, get data
    ! from emissions list 
    IF ( FRST .OR. ExtDat%FromList ) THEN

       ! Allocate temporary array
       ALLOCATE(Arr2D(HcoState%NX,HcoState%NY),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, "Arr2D allocation error", RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Try to get data from list
       CALL HCO_EvalFld( am_I_Root, HcoState, TRIM(FldName), Arr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.
  
             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
   
          ! Target to data
          ELSEIF ( PRESENT(Trgt) ) THEN
   
             ! If target is not associated: 
             IF ( .NOT. ASSOCIATED(Trgt) ) THEN
                IF ( FailIfNotFilled ) THEN
                   MSG = 'Cannot fill extension field ' // TRIM(FldName) // &
                         ' because target field is not associated.'
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
          
             ! If target is associated:
             ELSE
 
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
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
    
                ! Link data to target
                ExtDat%Arr%Val => Trgt
      
                ! Make sure it's not from list
                ExtDat%FromList = .FALSE.
      
                ! Mark as filled
                IF ( PRESENT(Filled) ) Filled = .TRUE.
   
                ! Verbose
                IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                   MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                   CALL HCO_MSG(HcoState%Config%Err,MSG)
                ENDIF
             ENDIF

          ! Not found in list and no target defined 
          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot fill extension field ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( FOUND ) THEN

             ! Copy values and mark as filled
             ExtDat%Arr%Val(:,:) = Arr2D(:,:)
             IF ( PRESENT(Filled) ) Filled = .TRUE.

          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
   
       ENDIF !FromList
    ENDIF 
   
    ! Make sure array exists
    IF ( FailIfNotFilled .AND. .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Cleanup
    IF ( ALLOCATED(Arr2D) ) DEALLOCATE(Arr2D)
 
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
  SUBROUTINE ExtDat_Set_3R ( am_I_Root, HcoState, ExtDat, FldName,   & 
                             RC,        First,    Trgt,   OnLevEdge, & 
                             Filled,    NotFillOk                     ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_CALC_MOD,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_3R),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(hp),         POINTER      , OPTIONAL         :: Trgt(:,:,:)
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: OnLevEdge 
    LOGICAL,          INTENT(  OUT), OPTIONAL         :: Filled
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: NotFillOk 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!  11 May 2015 - C. Keller - Now use HCO_EvalFld instead of HCO_GetPtr. This
!                            allows the application of scale factors to
!                            ExtState fields read through the HEMCO interface. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NX, NY, NZ, NZ_EXPECTED
    INTEGER                :: L
    LOGICAL                :: FRST
    LOGICAL                :: FOUND 
    LOGICAL                :: FailIfNotFilled 
    REAL(hp), ALLOCATABLE  :: Arr3D(:,:,:) 
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'ExtDat_Set_3R (hcox_state_mod.F90)'

    ! ================================================================
    ! ExtDat_Set_3R begins here
    ! ================================================================

    ! Init 
    RC = HCO_SUCCESS
    IF ( PRESENT(Filled) ) Filled = .FALSE.

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! Check for fill requirement
    IF ( PRESENT(NotFillOk) ) THEN
       FailIfNotFilled = .NOT. NotFillOk
    ELSE
       FailIfNotFilled = .TRUE.
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

       ! Allocate temporary array
       ALLOCATE(Arr3D(HcoState%NX,HcoState%NY,NZ_EXPECTED),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, "Arr3D allocation error", RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Try to get data from list
       CALL HCO_EvalFld( am_I_Root, HcoState, TRIM(FldName), Arr3D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.

             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ_EXPECTED, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
   
          ! Target to data
          ELSEIF ( PRESENT(Trgt) ) THEN
   
             ! If target is not associated: 
             IF ( .NOT. ASSOCIATED(Trgt) ) THEN
                IF ( FailIfNotFilled ) THEN
                   MSG = 'Cannot fill extension field ' // TRIM(FldName) // &
                         ' because target field is not associated.'
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
          
             ! If target is associated:
             ELSE
    
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
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
    
                ! Link data to target
                ExtDat%Arr%Val => Trgt
      
                ! Make sure it's not from list
                ExtDat%FromList = .FALSE.
      
                ! Mark as filled
                IF ( PRESENT(Filled) ) Filled = .TRUE.
   
                ! Verbose
                IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                   MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                   CALL HCO_MSG(HcoState%Config%Err,MSG)
                ENDIF
             ENDIF
          
          ! Not found in list and no target defined 
          ELSEIF ( FailIfNotFilled ) THEN
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( FOUND ) THEN

             ! Copy data and mark as filled 
             ExtDat%Arr%Val(:,:,:) = Arr3D(:,:,:)
             IF ( PRESENT(Filled) ) Filled = .TRUE.

          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN

          ENDIF 
       ENDIF !FromList
    ENDIF 

    ! Make sure array exists
    IF ( FailIfNotFilled .AND. .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Cleanup
    IF ( ALLOCATED(Arr3D) ) DEALLOCATE(Arr3D)
 
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
  SUBROUTINE ExtDat_Set_3S ( am_I_Root, HcoState, ExtDat, FldName,   & 
                             RC,        First,    Trgt,   OnLevEdge, &
                             Filled,    NotFillOk                     ) 
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_STATE_MOD,      ONLY : HCO_State
    USE HCO_CALC_MOD,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root
    TYPE(HCO_State),  POINTER                         :: HcoState
    TYPE(ExtDat_3S),  POINTER                         :: ExtDat
    CHARACTER(LEN=*), INTENT(IN   )                   :: FldName
    INTEGER,          INTENT(INOUT)                   :: RC     
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: First
    REAL(sp),         POINTER      , OPTIONAL         :: Trgt(:,:,:)
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: OnLevEdge 
    LOGICAL,          INTENT(  OUT), OPTIONAL         :: Filled
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: NotFillOk 
!
! !REVISION HISTORY:
!  03 Apr 2015 - C. Keller - Initial version
!  11 May 2015 - C. Keller - Now use HCO_EvalFld instead of HCO_GetPtr. This
!                            allows the application of scale factors to
!                            ExtState fields read through the HEMCO interface. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: AS, NX, NY, NZ, NZ_EXPECTED
    INTEGER                :: L
    LOGICAL                :: FRST
    LOGICAL                :: FOUND 
    LOGICAL                :: FailIfNotFilled 
    REAL(hp), ALLOCATABLE  :: Arr3D(:,:,:) 
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'ExtDat_Set_3S (hcox_state_mod.F90)'

    ! ================================================================
    ! ExtDat_Set_3S begins here
    ! ================================================================

    ! Init 
    RC = HCO_SUCCESS
    IF ( PRESENT(Filled) ) Filled = .FALSE.

    ! Nothing to do if this ExtDat field is not in use
    IF ( .NOT. ExtDat%DoUse ) RETURN

    ! First time
    IF ( PRESENT(FIRST) ) THEN
       FRST = FIRST
    ELSE
       FRST = .FALSE.
    ENDIF

    ! Check for fill requirement
    IF ( PRESENT(NotFillOk) ) THEN
       FailIfNotFilled = .NOT. NotFillOk
    ELSE
       FailIfNotFilled = .TRUE.
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

       ! Allocate temporary array
       ALLOCATE(Arr3D(HcoState%NX,HcoState%NY,NZ_EXPECTED),STAT=AS)
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, "Arr3D allocation error", RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Try to get data from list
       CALL HCO_EvalFld( am_I_Root, HcoState, TRIM(FldName), Arr3D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) RETURN     

       ! On first call, need to make additional checks
       IF ( FRST ) THEN
   
          ! If read from list
          IF ( FOUND ) THEN
             ExtDat%FromList = .TRUE.

             ! Make sure array is allocated
             CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ_EXPECTED, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
   
             ! Verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                MSG = 'Will fill extension field from HEMCO data list field ' // TRIM(FldName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF
   
          ! Target to data
          ELSEIF ( PRESENT(Trgt) ) THEN
  
             ! If target is not associated: 
             IF ( .NOT. ASSOCIATED(Trgt) ) THEN
                IF ( FailIfNotFilled ) THEN
                   MSG = 'Cannot fill extension field ' // TRIM(FldName) // &
                         ' because target field is not associated.'
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
          
             ! If target is associated:
             ELSE
    
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
                   CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
    
                ! Link data to target
                ExtDat%Arr%Val => Trgt
      
                ! Make sure it's not from list
                ExtDat%FromList = .FALSE.
      
                ! Mark as filled
                IF ( PRESENT(Filled) ) Filled = .TRUE.
   
                ! Verbose
                IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                   MSG = 'Set extension field pointer to external data: ' // TRIM(FldName)
                   CALL HCO_MSG(HcoState%Config%Err,MSG)
                ENDIF
             ENDIF
   
          ! Not found in list and no target defined 
          ELSEIF ( FailIfNotFilled ) THEN
             ! Target array must be present
             IF ( .NOT. PRESENT(Trgt) ) THEN
                MSG = 'Cannot fill extension field ' // TRIM(FldName)
                CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
          ENDIF
    
       ENDIF ! FIRST
   
       ! Eventually copy field from HEMCO list to ExtState. We need to
       ! make a copy and cannot just set a pointer because ExtState fields
       ! are in HEMCO precision but the EmisList fields are in single 
       ! precisions.
       IF ( ExtDat%FromList ) THEN
          IF ( FOUND ) THEN
             ! Copy data and mark as filled 
             ExtDat%Arr%Val(:,:,:) = Arr3D(:,:,:)
             IF ( PRESENT(Filled) ) Filled = .TRUE.
          ELSEIF ( FailIfNotFilled ) THEN
             MSG = 'Cannot find extension field in HEMCO data list: ' // TRIM(FldName)
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF 
       ENDIF !FromList
    ENDIF 

    ! Make sure array exists
    IF ( FailIfNotFilled .AND. .NOT. ASSOCIATED(ExtDat%Arr%Val) ) THEN
       MSG = 'ExtState array not filled: ' // TRIM(FldName)
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC=LOC )
    ENDIF
 
    ! Cleanup
    IF ( ALLOCATED(Arr3D) ) DEALLOCATE(Arr3D)
 
    ! Return w/ success
    RC = HCO_SUCCESS  

  END SUBROUTINE ExtDat_Set_3S
!EOC
END MODULE HCOX_STATE_MOD
