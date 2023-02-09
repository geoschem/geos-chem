!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: state_met_mod.F90
!
! !DESCRIPTION: Module STATE\_MET\_MOD contains the derived type
!  used to define the Meteorology State object for GEOS-Chem.
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory
!  to the Meteorology State object.  The Meteorology State object is not
!  defined in this module.  It must be be declared as variable in the top-level
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE:
!
MODULE State_Met_Mod
!
! USES:
!
  USE Cmn_Size_Mod, ONLY : NSURFTYPE
  USE Dictionary_M, ONLY : dictionary_t
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod


  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Met
  PUBLIC :: Cleanup_State_Met
  PUBLIC :: Get_Metadata_State_Met
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Init_and_Register
  PRIVATE :: Register_MetField
  PRIVATE :: Zero_State_Met
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Meteorology State
  !=========================================================================
  TYPE, PUBLIC :: MetState

     !----------------------------------------------------------------------
     ! Surface fields
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: ALBD          (:,:  ) ! Visible surface albedo [1]
     REAL(fp), POINTER :: AREA_M2       (:,:  ) ! Grid box surface area [m2]
     INTEGER,  POINTER :: ChemGridLev   (:,:  ) ! Chemistry grid level
     REAL(fp), POINTER :: CLDFRC        (:,:  ) ! Column cloud fraction [1]
     INTEGER,  POINTER :: CLDTOPS       (:,:  ) ! Max cloud top height [levels]
     REAL(fp), POINTER :: CONV_DEPTH    (:,:  ) ! Convective cloud depth [m]
     REAL(fp), POINTER :: EFLUX         (:,:  ) ! Latent heat flux [W/m2]
     REAL(fp), POINTER :: FLASH_DENS    (:,:  ) ! Lightning flash density [#/km2/s]
     REAL(fp), POINTER :: FRCLND        (:,:  ) ! Olson land fraction [1]
     REAL(fp), POINTER :: FRLAKE        (:,:  ) ! Fraction of lake [1]
     REAL(fp), POINTER :: FRLAND        (:,:  ) ! Fraction of land [1]
     REAL(fp), POINTER :: FRLANDIC      (:,:  ) ! Fraction of land ice [1]
     REAL(fp), POINTER :: FROCEAN       (:,:  ) ! Fraction of ocean [1]
     REAL(fp), POINTER :: FRSEAICE      (:,:  ) ! Sfc sea ice fraction
     REAL(fp), POINTER :: FRSNO         (:,:  ) ! Sfc snow fraction
     REAL(fp), POINTER :: GWETROOT      (:,:  ) ! Root soil wetness [1]
     REAL(fp), POINTER :: GWETTOP       (:,:  ) ! Top soil moisture [1]
     REAL(fp), POINTER :: HFLUX         (:,:  ) ! Sensible heat flux [W/m2]
     LOGICAL,  POINTER :: IsLand        (:,:  ) ! Is this a land  grid box?
     LOGICAL,  POINTER :: IsWater       (:,:  ) ! Is this a water grid box?
     LOGICAL,  POINTER :: IsIce         (:,:  ) ! Is this a ice   grid box?
     LOGICAL,  POINTER :: IsSnow        (:,:  ) ! Is this a snow  grid box?
     REAL(fp), POINTER :: LAI           (:,:  ) ! Leaf area index [m2/m2]
                                                !  (online)
     REAL(fp), POINTER :: PARDR         (:,:  ) ! Direct photsynthetically
                                                !  active radiation [W/m2]
     REAL(fp), POINTER :: PARDF         (:,:  ) ! Diffuse photsynthetically
                                                !  active radiation [W/m2]
     REAL(fp), POINTER :: PBLH          (:,:  ) ! PBL height [m]
     REAL(fp), POINTER :: PBL_TOP_hPa   (:,:  ) ! PBL top [hPa]
     REAL(fp), POINTER :: PBL_TOP_L     (:,:  ) ! PBL top [level]
     REAL(fp), POINTER :: PBL_TOP_m     (:,:  ) ! PBL top [m]
     REAL(fp), POINTER :: PBL_THICK     (:,:  ) ! PBL thickness [hPa]
     REAL(fp), POINTER :: PHIS          (:,:  ) ! Surface geopotential height
                                                !  [m2/s2]
     REAL(fp), POINTER :: PRECANV       (:,:  ) ! Anvil previp @ ground
                                                !  [kg/m2/s] -> [mm/day]
     REAL(fp), POINTER :: PRECCON       (:,:  ) ! Conv  precip @ ground
                                                !  [kg/m2/s] -> [mm/day]
     REAL(fp), POINTER :: PRECLSC       (:,:  ) ! Large-scale precip @ ground
                                                !  [kg/m2/s] -> [mm/day]
     REAL(fp), POINTER :: PRECTOT       (:,:  ) ! Total precip @ ground
                                                !  [kg/m2/s] -> [mm/day]
     REAL(fp), POINTER :: PS1_WET       (:,:  ) ! Wet surface pressure at
                                                !  start of timestep [hPa]
     REAL(fp), POINTER :: PS2_WET       (:,:  ) ! Wet surface pressure at
                                                !  end of timestep [hPa]
     REAL(fp), POINTER :: PSC2_WET      (:,:  ) ! Wet interpolated surface
                                                !  pressure [hPa]
     REAL(fp), POINTER :: PS1_DRY       (:,:  ) ! Dry surface pressure at
                                                !  start of timestep [hPa]
     REAL(fp), POINTER :: PS2_DRY       (:,:  ) ! Dry surface pressure at
                                                !  end of timestep [hPa]
     REAL(fp), POINTER :: PSC2_DRY      (:,:  ) ! Dry interpolated surface
                                                !  pressure [hPa]
     REAL(fp), POINTER :: QV2M          (:,:  ) ! Specific Humidity at 2m [kg/kg]
     REAL(fp), POINTER :: SEAICE00      (:,:  ) ! Sea ice coverage 00-10%
     REAL(fp), POINTER :: SEAICE10      (:,:  ) ! Sea ice coverage 10-20%
     REAL(fp), POINTER :: SEAICE20      (:,:  ) ! Sea ice coverage 20-30%
     REAL(fp), POINTER :: SEAICE30      (:,:  ) ! Sea ice coverage 30-40%
     REAL(fp), POINTER :: SEAICE40      (:,:  ) ! Sea ice coverage 40-50%
     REAL(fp), POINTER :: SEAICE50      (:,:  ) ! Sea ice coverage 50-60%
     REAL(fp), POINTER :: SEAICE60      (:,:  ) ! Sea ice coverage 60-70%
     REAL(fp), POINTER :: SEAICE70      (:,:  ) ! Sea ice coverage 70-80%
     REAL(fp), POINTER :: SEAICE80      (:,:  ) ! Sea ice coverage 80-90%
     REAL(fp), POINTER :: SEAICE90      (:,:  ) ! Sea ice coverage 90-100%
     REAL(fp), POINTER :: SLP           (:,:  ) ! Sea level pressure [hPa]
     REAL(fp), POINTER :: SNODP         (:,:  ) ! Snow depth [m]
     REAL(fp), POINTER :: SNOMAS        (:,:  ) ! Snow mass [kg/m2]
     REAL(fp), POINTER :: SUNCOS        (:,:  ) ! COS(solar zenith angle) at
                                                !   current time
     REAL(fp), POINTER :: SUNCOSmid     (:,:  ) ! COS(solar zenith angle) at
                                                !  midpoint of chem timestep
     REAL(fp), POINTER :: SUNCOSsum     (:,:  ) ! Sum of COS(SZA) for HEMCO OH
                                                !  diurnal variability
     REAL(fp), POINTER :: SZAFACT       (:,:  ) ! Diurnal scale factor for HEMCO OH
                                                !  diurnal variability (computed) [1]
     REAL(fp), POINTER :: SWGDN         (:,:  ) ! Incident radiation @ ground
                                                !  [W/m2]
     REAL(fp), POINTER :: TO3           (:,:  ) ! Total overhead O3 column [DU]
     REAL(fp), POINTER :: TROPP         (:,:  ) ! Tropopause pressure [hPa]
     INTEGER,  POINTER :: TropLev       (:,:  ) ! Tropopause level [1]
     REAL(fp), POINTER :: TropHt        (:,:  ) ! Tropopause height [km]
     REAL(fp), POINTER :: TS            (:,:  ) ! Surface temperature [K]
     REAL(fp), POINTER :: TSKIN         (:,:  ) ! Surface skin temperature [K]
     REAL(fp), POINTER :: U10M          (:,:  ) ! E/W wind speed @ 10m ht [m/s]
     REAL(fp), POINTER :: USTAR         (:,:  ) ! Friction velocity [m/s]
     REAL(fp), POINTER :: UVALBEDO      (:,:  ) ! UV surface albedo [1]
     REAL(fp), POINTER :: V10M          (:,:  ) ! N/S wind speed @ 10m ht [m/s]
     REAL(fp), POINTER :: Z0            (:,:  ) ! Surface roughness height [m]
     REAL(fp), POINTER :: CNV_FRC       (:,:  ) ! Convective fraction [1]

     !----------------------------------------------------------------------
     ! 3-D Fields
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: CLDF          (:,:,:) ! 3-D cloud fraction [1]
     REAL(fp), POINTER :: CMFMC         (:,:,:) ! Cloud mass flux [kg/m2/s]
     REAL(fp), POINTER :: DQRCU         (:,:,:) ! Conv precip production rate
                                                !  [kg/kg/s] (assume per
                                                !  dry air)
     REAL(fp), POINTER :: DQRLSAN       (:,:,:) ! LS precip prod rate [kg/kg/s]
                                                !  (assume per dry air)
     REAL(fp), POINTER :: DTRAIN        (:,:,:) ! Detrainment flux [kg/m2/s]
     REAL(fp), POINTER :: F_OF_PBL      (:,:,:) ! Fraction of box within PBL [1]
     REAL(fp), POINTER :: F_UNDER_PBLTOP(:,:,:) ! Fraction of box under PBL top
     REAL(fp), POINTER :: OMEGA         (:,:,:) ! Updraft velocity [Pa/s]
     REAL(fp), POINTER :: OPTD          (:,:,:) ! Visible optical depth [1]
     REAL(fp), POINTER :: PEDGE         (:,:,:) ! Wet air press @ level
                                                !  edges [hPa]
     REAL(fp), POINTER :: PFICU         (:,:,:) ! Dwn flux ice prec:conv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PFILSAN       (:,:,:) ! Dwn flux ice prec:LS+anv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PFLCU         (:,:,:) ! Dwn flux liq prec:conv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PFLLSAN       (:,:,:) ! Dwn flux ice prec:LS+anv
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: QI            (:,:,:) ! Ice mixing ratio
                                                !  [kg/kg dry air]
     REAL(fp), POINTER :: QL            (:,:,:) ! Water mixing ratio
                                                !  [kg/kg dry air]
     REAL(fp), POINTER :: REEVAPCN      (:,:,:) ! Evap of precip conv [kg/kg/s]
                                                !  (assume per dry air)
     REAL(fp), POINTER :: REEVAPLS      (:,:,:) ! Evap of precip LS+anvil
                                                !  [kg/kg/s] (assume per
                                                !  dry air)
     REAL(fp), POINTER :: RH            (:,:,:) ! Relative humidity [%]
     REAL(fp), POINTER :: SPHU          (:,:,:) ! Specific humidity
                                                !  [g H2O/kg tot air]
     REAL(fp), POINTER :: SPHU1         (:,:,:) ! Specific humidity at start
                                                !  of timestep [g/kg]
     REAL(fp), POINTER :: SPHU2         (:,:,:) ! Specific humidity at end
                                                !  of timestep [g/kg]
     REAL(fp), POINTER :: T             (:,:,:) ! Temperature [K]
     REAL(fp), POINTER :: TAUCLI        (:,:,:) ! Opt depth of ice clouds [1]
     REAL(fp), POINTER :: TAUCLW        (:,:,:) ! Opt depth of H2O clouds [1]
     REAL(fp), POINTER :: TMPU1         (:,:,:) ! Temperature at start of
                                                !  timestep [K]
     REAL(fp), POINTER :: TMPU2         (:,:,:) ! Temperature at end of
                                                !  timestep [K]
     REAL(fp), POINTER :: U             (:,:,:) ! E/W component of wind [m s-1]
     REAL(fp), POINTER :: UPDVVEL       (:,:,:) ! Updraft vertical velocity
                                                !  [hPa/s]
     REAL(fp), POINTER :: V             (:,:,:) ! N/S component of wind [m s-1]

     !----------------------------------------------------------------------
     ! Air quantities assigned in AIRQNT
     !----------------------------------------------------------------------
     ! Note on pressures: PMID is calculated from PEDGE,
     ! and dry air pressures assume constant RH and T across grid box
     REAL(fp), POINTER :: PEDGE_DRY     (:,:,:) ! Dry air partial pressure
                                                !  @ level edges [hPa]
     REAL(fp), POINTER :: PMID          (:,:,:) ! Average wet air pressure [hPa]
                                                !  defined as arithmetic
                                                !  average of edge pressures
     REAL(fp), POINTER :: PMID_DRY      (:,:,:) ! Dry air partial pressure [hPa]
                                                !  defined as arithmetic avg
                                                !  of edge pressures
     REAL(fp), POINTER :: THETA         (:,:,:) ! Potential temperature [K]
     REAL(fp), POINTER :: TV            (:,:,:) ! Virtual temperature [K]
     REAL(fp), POINTER :: MAIRDEN       (:,:,:) ! Moist air density [kg/m3]
     REAL(fp), POINTER :: AIRDEN        (:,:,:) ! Dry air density [kg/m3]
     REAL(fp), POINTER :: AIRNUMDEN     (:,:,:) ! Dry air density [molec/cm3]
     REAL(fp), POINTER :: AVGW          (:,:,:) ! Water vapor volume mixing
                                                !  ratio [vol H2O/vol dry air]
     REAL(fp), POINTER :: BXHEIGHT      (:,:,:) ! Grid box height [m] (dry air)
     REAL(fp), POINTER :: DELP          (:,:,:) ! Delta-P (wet) across box [hPa]
     REAL(fp), POINTER :: DELP_DRY      (:,:,:) ! Delta-P (dry) across box [hPa]
     REAL(fp), POINTER :: AD            (:,:,:) ! Dry air mass [kg] in grid box
     REAL(fp), POINTER :: AIRVOL        (:,:,:) ! Grid box volume [m3] (dry air)
     REAL(fp), POINTER :: DP_DRY_PREV   (:,:,:) ! Previous State_Met%DELP_DRY
     REAL(fp), POINTER :: SPHU_PREV     (:,:,:) ! Previous State_Met%SPHU

     !----------------------------------------------------------------------
     ! Offline land type, leaf area index, and chlorophyll fields
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: IREG          (:,:  ) ! # of landtypes in box (I,J)
     INTEGER,  POINTER :: ILAND         (:,:,:) ! Land type at (I,J);
                                                !  1..IREG(I,J)
     INTEGER,  POINTER :: IUSE          (:,:,:) ! Fraction (per mil) of box
                                                !  (I,J) occupied by each land
                                                !  type
     REAL(fp), POINTER :: MODISLAI      (:,:  ) ! Daily LAI computed from
                                                !  monthly offline MODIS [m2/m2]
     REAL(fp), POINTER :: XLAI          (:,:,:) ! MODIS LAI per land type,
                                                !  for this month
     REAL(fp), POINTER :: LandTypeFrac  (:,:,:) ! Olson frac per type (I,J,type)
     REAL(fp), POINTER :: XLAI_NATIVE   (:,:,:) ! avg LAI per type (I,J,type)
     REAL(fp), POINTER :: XLAI2         (:,:,:) ! MODIS LAI per land type,
                                                !  for next month

     !----------------------------------------------------------------------
     ! Fields for querying in which vertical regime a grid box is in
     ! or if a grid box is near local noon solar time
     !----------------------------------------------------------------------
     LOGICAL,  POINTER :: InChemGrid    (:,:,:) ! Are we in the chemistry grid?
     LOGICAL,  POINTER :: InPbl         (:,:,:) ! Are we in the PBL?
     LOGICAL,  POINTER :: InStratMeso   (:,:,:) ! Are we in the stratosphere
                                                !            or mesosphere?
     LOGICAL,  POINTER :: InStratosphere(:,:,:) ! Are we in the stratosphere?
     LOGICAL,  POINTER :: InTroposphere (:,:,:) ! Are we in the troposphere?
     REAL(fp), POINTER :: LocalSolarTime(:,:  ) ! Local solar time
     LOGICAL,  POINTER :: IsLocalNoon   (:,:  ) ! Is it local noon (between 11
                                                !  and 13 local solar time?

     !----------------------------------------------------------------------
     ! Fields for wet scavenging module
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: C_H2O         (:,:,:) ! Mix ratio of H2O [v/v]
     REAL(fp), POINTER :: CLDICE        (:,:,:) ! Cloud ice mixing ratio [cm3 ice/cm3 air]
     REAL(fp), POINTER :: CLDLIQ        (:,:,:) ! Cloud liquid water mixing ratio [cm3 H2O/cm3 air]
     REAL(fp), POINTER :: PDOWN         (:,:,:) ! Precipitation thru the bottom of the grid box
                                                ! [cm3 H2O/cm2 area/s]
     REAL(fp), POINTER :: QQ            (:,:,:) ! Rate of new precip formation [cm3 H2O/cm3 air/s]
     REAL(fp), POINTER :: REEVAP        (:,:,:) ! Rate of precip reevaporation
     REAL(fp), POINTER :: PSO4_SO2APM2  (:,:,:)

     !----------------------------------------------------------------------
     ! Fields for boundary layer mixing
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: IMIX          (:,:  ) ! Integer and fractional level
     REAL(fp), POINTER :: FPBL          (:,:  ) !  where PBL top occurs
     INTEGER           :: PBL_MAX_L             ! Max level where PBL top occurs

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Met
     !----------------------------------------------------------------------
     CHARACTER(LEN=3)             :: State     = 'MET'    ! Name of this state
     TYPE(MetaRegItem), POINTER   :: Registry  => NULL()  ! Registry object
     TYPE(dictionary_t)           :: RegDict              ! Reg. lookup table

  END TYPE MetState
!
! !REMARKS:
!  In MERRA2, PS and SLP are kept in Pa (not converted to hPa).
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, split off from gc_type_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Init_and_Register
     MODULE PROCEDURE Init_and_Register_R4_2D
     MODULE PROCEDURE Init_and_Register_R4_3D
     MODULE PROCEDURE Init_and_Register_R8_2D
     MODULE PROCEDURE Init_and_Register_R8_3D
     MODULE PROCEDURE Init_and_Register_Log_2D
     MODULE PROCEDURE Init_and_Register_Log_3D
     MODULE PROCEDURE Init_and_Register_Int_2D
     MODULE PROCEDURE Init_and_Register_Int_3D
  END INTERFACE Init_and_Register

  INTERFACE Register_MetField
     MODULE PROCEDURE Register_MetField_R4_2D
     MODULE PROCEDURE Register_MetField_R4_3D
     MODULE PROCEDURE Register_MetField_R8_2D
     MODULE PROCEDURE Register_MetField_R8_3D
     MODULE PROCEDURE Register_MetField_Int_2D
     MODULE PROCEDURE Register_MetField_Int_3D
  END INTERFACE Register_MetField

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Zero_State_Met
!
! !DESCRIPTION: Nullifies and/or zeroes all fields of State\_Met.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Zero_State_Met( State_Met, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Initialize
    RC = GC_SUCCESS

    !=======================================================================
    ! Nullify all fields for safety's sake before allocating them
    ! This can prevent compilation errors caused by uninitialized values
    !=======================================================================
    State_Met%ALBD           => NULL()
    State_Met%AREA_M2        => NULL()
    State_Met%ChemGridLev    => NULL()
    State_Met%CLDFRC         => NULL()
    State_Met%CLDTOPS        => NULL()
    State_Met%CONV_DEPTH     => NULL()
    State_Met%EFLUX          => NULL()
    State_Met%FLASH_DENS     => NULL()
    State_Met%FRCLND         => NULL()
    State_Met%FRLAKE         => NULL()
    State_Met%FRLAND         => NULL()
    State_Met%FRLANDIC       => NULL()
    State_Met%FROCEAN        => NULL()
    State_Met%FRSEAICE       => NULL()
    State_Met%FRSNO          => NULL()
    State_Met%GWETROOT       => NULL()
    State_Met%GWETTOP        => NULL()
    State_Met%HFLUX          => NULL()
    State_Met%IsLand         => NULL()
    State_Met%IsWater        => NULL()
    State_Met%IsIce          => NULL()
    State_Met%IsSnow         => NULL()
    State_Met%LAI            => NULL()
    State_Met%PARDR          => NULL()
    State_Met%PARDF          => NULL()
    State_Met%PBLH           => NULL()
    State_Met%PBL_TOP_hPa    => NULL()
    State_Met%PBL_TOP_L      => NULL()
    State_Met%PBL_TOP_m      => NULL()
    State_Met%PBL_THICK      => NULL()
    State_Met%PHIS           => NULL()
    State_Met%PRECANV        => NULL()
    State_Met%PRECCON        => NULL()
    State_Met%PRECLSC        => NULL()
    State_Met%PRECTOT        => NULL()
    State_Met%PS1_WET        => NULL()
    State_Met%PS2_WET        => NULL()
    State_Met%PSC2_WET       => NULL()
    State_Met%PS1_DRY        => NULL()
    State_Met%PS2_DRY        => NULL()
    State_Met%PSC2_DRY       => NULL()
    State_Met%QV2M           => NULL()
    State_Met%SEAICE00       => NULL()
    State_Met%SEAICE10       => NULL()
    State_Met%SEAICE20       => NULL()
    State_Met%SEAICE30       => NULL()
    State_Met%SEAICE40       => NULL()
    State_Met%SEAICE50       => NULL()
    State_Met%SEAICE60       => NULL()
    State_Met%SEAICE70       => NULL()
    State_Met%SEAICE80       => NULL()
    State_Met%SEAICE90       => NULL()
    State_Met%SLP            => NULL()
    State_Met%SNODP          => NULL()
    State_Met%SNOMAS         => NULL()
    State_Met%SUNCOS         => NULL()
    State_Met%SUNCOSmid      => NULL()
    State_Met%SUNCOSsum      => NULL()
    State_Met%SZAFACT        => NULL()
    State_Met%SWGDN          => NULL()
    State_Met%TO3            => NULL()
    State_Met%TROPP          => NULL()
    State_Met%TropLev        => NULL()
    State_Met%TropHt         => NULL()
    State_Met%TS             => NULL()
    State_Met%TSKIN          => NULL()
    State_Met%U10M           => NULL()
    State_Met%USTAR          => NULL()
    State_Met%UVALBEDO       => NULL()
    State_Met%V10M           => NULL()
    State_Met%Z0             => NULL()
    State_Met%CNV_FRC        => NULL()
    State_Met%CLDF           => NULL()
    State_Met%CMFMC          => NULL()
    State_Met%DQRCU          => NULL()
    State_Met%DQRLSAN        => NULL()
    State_Met%DTRAIN         => NULL()
    State_Met%F_OF_PBL       => NULL()
    State_Met%F_UNDER_PBLTOP => NULL()
    State_Met%OMEGA          => NULL()
    State_Met%OPTD           => NULL()
    State_Met%PEDGE          => NULL()
    State_Met%PFICU          => NULL()
    State_Met%PFILSAN        => NULL()
    State_Met%PFLCU          => NULL()
    State_Met%PFLLSAN        => NULL()
    State_Met%QI             => NULL()
    State_Met%QL             => NULL()
    State_Met%REEVAPCN       => NULL()
    State_Met%REEVAPLS       => NULL()
    State_Met%RH             => NULL()
    State_Met%SPHU           => NULL()
    State_Met%SPHU1          => NULL()
    State_Met%SPHU2          => NULL()
    State_Met%T              => NULL()
    State_Met%TAUCLI         => NULL()
    State_Met%TAUCLW         => NULL()
    State_Met%TMPU1          => NULL()
    State_Met%TMPU2          => NULL()
    State_Met%U              => NULL()
    State_Met%UPDVVEL        => NULL()
    State_Met%V              => NULL()
    State_Met%PEDGE_DRY      => NULL()
    State_Met%PMID           => NULL()
    State_Met%PMID_DRY       => NULL()
    State_Met%THETA          => NULL()
    State_Met%TV             => NULL()
    State_Met%MAIRDEN        => NULL()
    State_Met%AIRDEN         => NULL()
    State_Met%AIRNUMDEN      => NULL()
    State_Met%AVGW           => NULL()
    State_Met%BXHEIGHT       => NULL()
    State_Met%DELP           => NULL()
    State_Met%DELP_DRY       => NULL()
    State_Met%AD             => NULL()
    State_Met%AIRVOL         => NULL()
    State_Met%DP_DRY_PREV    => NULL()
    State_Met%SPHU_PREV      => NULL()
    State_Met%IREG           => NULL()
    State_Met%ILAND          => NULL()
    State_Met%IUSE           => NULL()
    State_Met%MODISLAI       => NULL()
    State_Met%XLAI           => NULL()
    State_Met%LandTypeFrac   => NULL()
    State_Met%XLAI_NATIVE    => NULL()
    State_Met%XLAI2          => NULL()
    State_Met%InChemGrid     => NULL()
    State_Met%InPbl          => NULL()
    State_Met%InStratMeso    => NULL()
    State_Met%InStratosphere => NULL()
    State_Met%InTroposphere  => NULL()
    State_Met%LocalSolarTime => NULL()
    State_Met%IsLocalNoon    => NULL()
    State_Met%C_H2O          => NULL()
    State_Met%CLDICE         => NULL()
    State_Met%CLDLIQ         => NULL()
    State_Met%PDOWN          => NULL()
    State_Met%QQ             => NULL()
    State_Met%REEVAP         => NULL()
    State_Met%IMIX           => NULL()
    State_Met%FPBL           => NULL()
    State_Met%REEVAP         => NULL()
    State_Met%PSO4_SO2APM2   => NULL()
    State_Met%PBL_MAX_L      = 0

  END SUBROUTINE Zero_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_State_Met
!
! !DESCRIPTION: Subroutine INIT\_STATE\_MET allocates all fields of
!  the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Met( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: errMsg_ir, thisLoc, metId
    CHARACTER(LEN=512) :: errMsg

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    errMsg     =  ''
    errMsg_ir  =  'Error encountered in "Init_and_Register", metId = '
    thisLoc    =  &
         ' -> at Init_State_Met (in module Headers/state_met_mod.F90)'

    ! Nullify or zero all State_Met variables
    CALL Zero_State_Met( State_Met, RC )

    !========================================================================
    ! Exit if this is a dry-run simulation
    !========================================================================
    IF ( Input_Opt%DryRun ) THEN
       RC = GC_SUCCESS
       RETURN
    ENDIF

    !========================================================================
    ! Allocate 2-D Fields
    !========================================================================

    !------------------------------------------------------------------------
    ! ALBD [1]
    !------------------------------------------------------------------------
    metId = 'ALBD'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%ALBD,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! AREA_M2 [m2]
    !------------------------------------------------------------------------
    metId = 'AREAM2'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AREA_M2,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! ChemGridLev [1]
    !------------------------------------------------------------------------
    metId = 'ChemGridLev'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%ChemGridLev,                                 &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! CLDFRC [1]
    !------------------------------------------------------------------------
    metId = 'CLDFRC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CLDFRC,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! CLDTOPS [level]
    !------------------------------------------------------------------------
    metId = 'CLDTOPS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CLDTOPS,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef MODEL_GEOS
    !------------------------------------------------------------------------
    ! CNV_FRC [1]
    ! Convective fractions are not yet a standard GEOS-FP
    ! field. Only available to online model (ckeller, 3/4/16)
    !------------------------------------------------------------------------
    metId = 'CNVFRC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CNV_FRC,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    !------------------------------------------------------------------------
    ! Convective Depth [m]
    !------------------------------------------------------------------------
    metId = 'CONVDEPTH'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CONV_DEPTH,                                  &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! EFLUX [W m-2]
    !------------------------------------------------------------------------
    metId = 'EFLUX'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%EFLUX,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !---------------------------------------------------------------------
    ! Lightning density [#/km2/s]
    !---------------------------------------------------------------------
    metId = 'FLASHDENS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FLASH_DENS,                                  &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FPBL [1] : Local variable for PBL mixing -- do not register this
    !------------------------------------------------------------------------
    metId = 'FPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FPBL,                                        &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRCLND [1]
    !------------------------------------------------------------------------
    metId = 'FRCLND'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRCLND,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRLAKE [1]
    !------------------------------------------------------------------------
    metId = 'FRLAKE'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRLAKE,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRLAND [1]
    !------------------------------------------------------------------------
    metId = 'FRLAND'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRLAND,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRLANDIC [1]
    !------------------------------------------------------------------------
    metId = 'FRLANDIC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRLANDIC,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FROCEAN [1]
    !------------------------------------------------------------------------
    metId = 'FROCEAN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FROCEAN,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRESEAICE [1]
    !------------------------------------------------------------------------
    metId = 'FRSEAICE'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRSEAICE,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! FRSNO [1]
    !------------------------------------------------------------------------
    metId = 'FRSNO'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%FRSNO,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! GWETROOT [1]
    !------------------------------------------------------------------------
    metId = 'GWETROOT'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%GWETROOT,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! GWETTOP [1]
    !------------------------------------------------------------------------
    metId = 'GWETTOP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%GWETTOP,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! HFLUX [W m-2]
    !------------------------------------------------------------------------
    metId = 'HFLUX'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%HFLUX,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IMIX [1]: Local variable for PBL mixing -- Do not register this
    !------------------------------------------------------------------------
    metId = 'IMIX'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IMIX,                                        &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IREG [1]
    !------------------------------------------------------------------------
    metId = 'IREG'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IREG,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IsWater (do not register for diagnostics)
    !------------------------------------------------------------------------
    metId = 'IsWater'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IsWater,                                     &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IsWater (do not register for diagnostics)
    !------------------------------------------------------------------------
    metId = 'IsLand'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IsLand,                                      &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IsIce (do not register for diagnostics)
    !------------------------------------------------------------------------
    metId = 'IsIce'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IsIce,                                       &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IsSnow (do not register for diagnostics)
    !------------------------------------------------------------------------
    metId = 'IsSnow'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IsSnow,                                      &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! LAI [1]
    !------------------------------------------------------------------------
    metId = 'LAI'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%LAI,                                         &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! MODISLAI [1]
    !------------------------------------------------------------------------
    metId = 'MODISLAI'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%MODISLAI,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PARDR [W m-2]
    !------------------------------------------------------------------------
    metId = 'PARDR'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PARDR,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PARDF [W m-2]
    !------------------------------------------------------------------------
    metId = 'PARDF'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PARDF,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PBLH [m]
    !------------------------------------------------------------------------
    metId = 'PBLH'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PBLH,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PBL top [hPa]
    !------------------------------------------------------------------------
    metId = 'PBLTOPhPa'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PBL_TOP_hPa,                                 &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PBL top [level]
    !------------------------------------------------------------------------
    metId = 'PBLTOPL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PBL_TOP_L,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PBL top [m]
    !------------------------------------------------------------------------
    metId = 'PBLTOPM'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PBL_TOP_m,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PBL thickness [hPa]
    !------------------------------------------------------------------------
    metId = 'PBLTHICK'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PBL_THICK,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PHIS [m2 s-2], converted to [m] after data read
    !------------------------------------------------------------------------
    metId = 'PHIS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PHIS,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PRECANV [kg m-2 s-1], converted to [mm day-1]
    !------------------------------------------------------------------------
    metId = 'PRECANV'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PRECANV,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PRECCON [kg m-2 s-1], converted to [mm day-1]
    !------------------------------------------------------------------------
    metId = 'PRECCON'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PRECCON,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PRECLSC [kg m-2 s-1], converted to [mm day-1]
    !------------------------------------------------------------------------
    metId = 'PRECLSC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PRECLSC,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PRECTOT [kg m-2 s-1], converted to [mm day-1]
    !------------------------------------------------------------------------
    metId = 'PRECTOT'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PRECTOT,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PS1_WET [hPa]
    !------------------------------------------------------------------------
    metId = 'PS1WET'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PS1_WET,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PS2_WET [hPa]
    !------------------------------------------------------------------------
    metId = 'PS2WET'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PS2_WET,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PSC2_WET [hPa]
    !------------------------------------------------------------------------
    metId = 'PSC2WET'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PSC2_WET,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PS1_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'PS1DRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PS1_DRY,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PS2_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'PS2DRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PS2_DRY,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PSC2_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'PSC2DRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PSC2_DRY,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !-------------------------
    ! QV2M [kg/kg]
    !-------------------------
    metId = 'QV2M'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%QV2M,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE00 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE00'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE00,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE10 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE10'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE10,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE20 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE20'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE20,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE30 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE30'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE30,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE40 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE40'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE40,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE50 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE50'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE50,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE60 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE60'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE60,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE70 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE70'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE70,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE80 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE80'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE80,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SEAICE90 [1]
    !------------------------------------------------------------------------
    metId = 'SEAICE90'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SEAICE90,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SLP [hPa]
    !------------------------------------------------------------------------
    metId = 'SLP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SLP,                                         &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SNODP [m]
    !------------------------------------------------------------------------
    metId = 'SNODP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SNODP,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SNOMAS [kg m-2]
    !------------------------------------------------------------------------
    metId = 'SNOMAS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SNOMAS,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SUNCOS [1]
    !------------------------------------------------------------------------
    metId = 'SUNCOS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SUNCOS,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SUNCOSmid [1]
    !------------------------------------------------------------------------
    metId = 'SUNCOSmid'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SUNCOSmid,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SUNCOSsum [1] (for HEMCO)
    !------------------------------------------------------------------------
    metId = 'SUNCOSsum'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SUNCOSsum,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SZAFACT [1] (for HEMCO)
    !------------------------------------------------------------------------
    metId = 'SZAFACT'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SZAFACT,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SWGDN [W m-2]
    !------------------------------------------------------------------------
    metId = 'SWGDN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SWGDN,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TO3 [dobsons]
    !------------------------------------------------------------------------
    metId = 'TO3'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TO3,                                         &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TropLev [1]
    !------------------------------------------------------------------------
    metId = 'TropLev'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TropLev,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TropHt [m]
    !------------------------------------------------------------------------
    metId = 'TropHt'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TropHt,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TropP [hPa]
    !------------------------------------------------------------------------
    metId = 'TropP '
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TropP,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TS [K]
    !------------------------------------------------------------------------
    metId = 'TS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TS,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TSKIN [K]
    !------------------------------------------------------------------------
    metId = 'TSKIN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TSKIN,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! U10M [m s-1]
    !------------------------------------------------------------------------
    metId = 'U10M'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%U10M,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! USTAR [m -s]
    !------------------------------------------------------------------------
    metId = 'USTAR'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%USTAR,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! UVALBEDO [1]
    !------------------------------------------------------------------------
    metId = 'UVALBEDO'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%UVALBEDO,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! V10M [m s-1]
    !------------------------------------------------------------------------
    metId = 'V10M'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%V10M,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Z0 [m]
    !------------------------------------------------------------------------
    metId = 'Z0'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%Z0,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================

    !------------------------------------------------------------------------
    ! AD [kg]
    !------------------------------------------------------------------------
    metId = 'AD'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AD,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! AIRDEN [kg m-3]
    !------------------------------------------------------------------------
    metId = 'AIRDEN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AIRDEN,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! AIRNUMDEN [1]
    !------------------------------------------------------------------------
    metId = 'AIRNUMDEN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AIRNUMDEN,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! AIRVOL [m3]
    !------------------------------------------------------------------------
    metId = 'AIRVOL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AIRVOL,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! AVGW [v/v]
    !------------------------------------------------------------------------
    metId = 'AVGW'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%AVGW,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! BXHEIGHT [m]
    !------------------------------------------------------------------------
    metId = 'BXHEIGHT'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%BXHEIGHT,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! CLDF [1]
    !------------------------------------------------------------------------
    metId = 'CLDF'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CLDF,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! CMFMC [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'CMFMC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%CMFMC,                                       &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DELP [hPa]
    !------------------------------------------------------------------------
    metId = 'DELP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DELP,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DELP_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'DELPDRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DELP_DRY,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DP_DRY_PREV [hPa]
    !------------------------------------------------------------------------
    metId = 'DPDRYPREV'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DP_DRY_PREV,                                 &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Fraction of PBL
    !------------------------------------------------------------------------
    metId = 'FOFPBL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%F_OF_PBL,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Fraction of box under PBL top
    !------------------------------------------------------------------------
    metId = 'FUNDERPBLTOP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%F_UNDER_PBLTOP,                              &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DQRCU [kg kg-1 s-1]
    !------------------------------------------------------------------------
    metId = 'DQRCU'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DQRCU,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DQRLSAN [kg kg-1 s-1]
    !------------------------------------------------------------------------
    metId = 'DQRLSAN '
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DQRLSAN,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! DTRAIN [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'DTRAIN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%DTRAIN,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! ILAND [1]
    !------------------------------------------------------------------------
    metId = 'ILAND'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%ILAND,                                       &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IUSE [1]
    !------------------------------------------------------------------------
    metId = 'IUSE'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IUSE,                                        &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! LANDTYPEFRAC [1]
    !------------------------------------------------------------------------
    metId = 'LANDTYPEFRAC'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%LANDTYPEFRAC,                                &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! MAIRDEN [kg m-3]
    !------------------------------------------------------------------------
    metId = 'MAIRDEN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%MAIRDEN,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! OMEGA [Pa s-1]
    !------------------------------------------------------------------------
    metId = 'OMEGA'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%OMEGA,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! OPTD [1]
    !------------------------------------------------------------------------
    metId = 'OPTD'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%OPTD,                                         &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PEDGE [hPa]
    !------------------------------------------------------------------------
    metId = 'PEDGE'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PEDGE,                                       &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PEDGE_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'PEDGEDRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PEDGE_DRY,                                   &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PFICU [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'PFICU'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PFICU,                                       &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PFILSAN [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'PFILSAN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PFILSAN,                                     &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PFLCU [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'PFLCU'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PFLCU,                                       &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PFLLSAN [kg m-2 s-1]
    !------------------------------------------------------------------------
    metId = 'PFLLSAN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PFLLSAN,                                     &
         onEdges    = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PMID [1]
    !------------------------------------------------------------------------
    metId = 'PMID'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PMID,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! PMID_DRY [hPa]
    !------------------------------------------------------------------------
    metId = 'PMIDDRY'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%PMID_DRY,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! QI [kg kg-1]
    !------------------------------------------------------------------------
    metId = 'QI'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%QI,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! QL [kg kg-1]
    !------------------------------------------------------------------------
    metId = 'QL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%QL,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! REEVAPCN [kg kg-1 s-1]
    !------------------------------------------------------------------------
    metId = 'REEVAPCN'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%REEVAPCN,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! REEVAPLS [kg kg-1 s-1]
    !------------------------------------------------------------------------
    metId = 'REEVAPLS'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%REEVAPLS,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! RH [%]
    !------------------------------------------------------------------------
    metId = 'RH'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%RH,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SPHU [g kg-1]
    !------------------------------------------------------------------------
    metId = 'SPHU'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SPHU,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SPHU1 [g kg-1]
    !------------------------------------------------------------------------
    metId = 'SPHU1'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SPHU1,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SPHU2 [g kg-1]
    !------------------------------------------------------------------------
    metId = 'SPHU2'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SPHU2,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SPHU_PREV [g/kg]
    !------------------------------------------------------------------------
    metId = 'SPHUPREV'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%SPHU_PREV,                                   &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! T [K]
    !------------------------------------------------------------------------
    metId = 'T'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%T,                                           &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TAUCLI [1]
    !------------------------------------------------------------------------
    metId = 'TAUCLI'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TAUCLI,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TAUCLW [1]
    !------------------------------------------------------------------------
    metId = 'TAUCLW'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TAUCLW,                                      &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! THETA [K]
    !------------------------------------------------------------------------
    metId = 'THETA'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%THETA,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TMPU1 [K]
    !------------------------------------------------------------------------
    metId = 'TMPU1'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TMPU1,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TMPU2 [K]
    !------------------------------------------------------------------------
    metId = 'TMPU2'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TMPU2,                                       &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! TV [K]
    !------------------------------------------------------------------------
    metId = 'TV'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%TV,                                          &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! U [m s-1]
    !------------------------------------------------------------------------
    metId = 'U'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%U,                                           &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef MODEL_GEOS
    !------------------------------------------------------------------------
    ! UPDVVEL [hPa s-1]
    ! Updraft vertical velocity is not yet a standard GEOS-FP
    ! field. Only available to online model (ckeller, 3/4/16)
    !------------------------------------------------------------------------
    metId = 'UPDVVEL'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%UPDVVEL,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    !------------------------------------------------------------------------
    ! V [m s-1]
    !------------------------------------------------------------------------
    metId = 'V'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%V,                                           &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! XLAI [1]
    !------------------------------------------------------------------------
    metId = 'XLAI'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%XLAI,                                        &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! XLAI2 [1]
    !------------------------------------------------------------------------
    metId = 'XLAI2'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%XLAI2,                                       &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! XLAI_NATIVE [1]
    !------------------------------------------------------------------------
    metId = 'XLAINATIVE'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%XLAI_NATIVE,                                 &
         nSlots     = NSURFTYPE,                                             &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Allocate fields used by wet scavenging and other GeosCore modules.
    ! Note some are memory-order ZXY and some are regular XYZ.
    ! Only allocate arrays if wetdep or convection is turned on
    !========================================================================
    IF ( Input_Opt%LWETD .or. Input_Opt%LCONV ) THEN

       !---------------------------------------------------------------------
       ! C_H2O
       !---------------------------------------------------------------------
       metId = 'C_H2O'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%C_H2O,                                    &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! CLDICE
       !-----------------------------------------------------------------
       metId = 'CLDICE'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%CLDICE,                                   &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !-----------------------------------------------------------------
       ! CLDLIQ
       !-----------------------------------------------------------------
       metId = 'CLDLIQ'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%CLDLIQ,                                   &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! PDOWN (ZXY order)
       !---------------------------------------------------------------------
       metId = 'PDOWN'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%PDOWN,                                    &
            noRegister = .TRUE.,                                             &
            zxyOrder   = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! QQ (ZXY order)
       !---------------------------------------------------------------------
       metId = 'QQ'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%QQ,                                       &
            noRegister = .TRUE.,                                             &
            zxyOrder   = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! REEVAP (ZXY order)
       !---------------------------------------------------------------------
       metId = 'REEVAP'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%REEVAP,                                   &
            noRegister = .TRUE.,                                             &
            zxyOrder   = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

#ifdef APM
       !---------------------------------------------------------------------
       ! PSO4_SO2APM2
       !---------------------------------------------------------------------
       metId = 'PSO4SO2APM2'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Met  = State_Met,                                          &
            State_Grid = State_Grid,                                         &
            metId      = metId,                                              &
            Ptr2Data   = State_Met%PSO4_SO2APM2,                             &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( metId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
#endif

    ENDIF

    !=======================================================================
    ! Allocate fields for querying which vertical regime a grid box is in
    ! or if a grid box is near local solar noontime.
    !
    ! %%%%% NOTE: Most of these are logical fields and thus %%%%%
    ! %%%%%   cannot be archived to HISTORY diagnostics.    %%%%%
    !=======================================================================

    !------------------------------------------------------------------------
    ! InChemGrid
    !------------------------------------------------------------------------
    metId = 'InChemGrid'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%InChemGrid,                                  &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! InPBL
    !------------------------------------------------------------------------
    metId = 'InPbl'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%InPbl,                                       &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! InStratosphere
    !------------------------------------------------------------------------
    metId = 'InStratosphere'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%InStratosphere,                              &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! InStratMeso
    !------------------------------------------------------------------------
    metId = 'InStratMeso'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%InStratMeso,                                 &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! InTroposphere
    !------------------------------------------------------------------------
    metId = 'InTroposphere'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%InTroposphere,                               &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! IsLocalNoon
    !------------------------------------------------------------------------
    metId = 'IsLocalNoon'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%IsLocalNoon,                                 &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! LocalSolarTime (register this for diagnostics)
    !------------------------------------------------------------------------
    metId = 'LocalSolarTime'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Met  = State_Met,                                             &
         State_Grid = State_Grid,                                            &
         metId      = metId,                                                 &
         Ptr2Data   = State_Met%LocalSolarTime,                              &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( metId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Once we are done registering all fields, we need to define the
    ! registry lookup table.  This algorithm will avoid hash collisions.
    !========================================================================
    CALL Registry_Set_LookupTable( Registry  = State_Met%Registry,           &
                                   RegDict   = State_Met%RegDict,            &
                                   RC        = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Set_LookupTable"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Print information about the registered fields (short format)
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 10 )
10     FORMAT(/, 'Registered variables contained within the State_Met object:')
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( Input_Opt   = Input_Opt,                            &
                         Registry    = State_Met%Registry,                   &
                         ShortFormat = .TRUE.,                               &
                         RC          = RC                                   )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

   END SUBROUTINE Init_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Met
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_MET deallocates all fields
!  of the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Met( State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Obj for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Cleanup_State_Met (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Deallocate 2-D fields
    !========================================================================
    IF ( ASSOCIATED( State_Met%ALBD ) ) THEN
       DEALLOCATE( State_Met%ALBD, STAT=RC )
       CALL GC_CheckVar( 'State_Met%ALBD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%ALBD => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%AREA_M2 ) ) THEN
       DEALLOCATE( State_Met%AREA_M2, STAT=RC )
       CALL GC_CheckVar( 'State_Met%AREA_M2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AREA_M2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%ChemGridLev ) ) THEN
       DEALLOCATE( State_Met%ChemGridLev, STAT=RC )
       CALL GC_CheckVar( 'State_Met%ChemGridLev', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%ChemGridLev => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CLDFRC ) ) THEN
       DEALLOCATE( State_Met%CLDFRC, STAT=RC )
       CALL GC_CheckVar( 'State_Met%CLDFRC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CLDFRC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CLDTOPS ) ) THEN
       DEALLOCATE( State_Met%CLDTOPS, STAT=RC )
       CALL GC_CheckVar( 'State_Met%CLDTOPS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CLDTOPS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CONV_DEPTH ) ) THEN
       DEALLOCATE( State_Met%CONV_DEPTH, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%CONV_DEPTH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CONV_DEPTH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%EFLUX ) ) THEN
       DEALLOCATE( State_Met%EFLUX, STAT=RC )
       CALL GC_CheckVar( 'State_Met%EFLUX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%EFLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FLASH_DENS ) ) THEN
       DEALLOCATE( State_Met%FLASH_DENS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%FLASH_DENS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FLASH_DENS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRCLND ) ) THEN
       DEALLOCATE( State_Met%FRCLND, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRCLND', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRCLND => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRLAKE ) ) THEN
       DEALLOCATE( State_Met%FRLAKE, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRLAKE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRLAKE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRLAND ) ) THEN
       DEALLOCATE( State_Met%FRLAND, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRLAND', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRLAND => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRLANDIC ) ) THEN
       DEALLOCATE( State_Met%FRLANDIC, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRLANDIC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRLANDIC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FROCEAN ) ) THEN
       DEALLOCATE( State_Met%FROCEAN, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FROCEAN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FROCEAN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRSEAICE ) ) THEN
       DEALLOCATE( State_Met%FRSEAICE, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRSEAICE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRSEAICE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FRSNO ) ) THEN
       DEALLOCATE( State_Met%FRSNO, STAT=RC )
       CALL GC_CheckVar( 'State_Met%FRSNO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FRSNO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%GWETROOT ) ) THEN
       DEALLOCATE( State_Met%GWETROOT, STAT=RC )
       CALL GC_CheckVar( 'State_Met%GWETROOT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%GWETROOT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%GWETTOP ) ) THEN
       DEALLOCATE( State_Met%GWETTOP, STAT=RC )
       CALL GC_CheckVar( 'State_Met%GWETTOP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%GWETTOP => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%HFLUX ) ) THEN
       DEALLOCATE( State_Met%HFLUX, STAT=RC )
       CALL GC_CheckVar( 'State_Met%HFLUX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%HFLUX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IsWater ) ) THEN
       DEALLOCATE( State_Met%IsWater, STAT=RC )
       CALL GC_CheckVar( 'State_Met%IsWater', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IsWater => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IsLand ) ) THEN
       DEALLOCATE( State_Met%IsLand, STAT=RC )
       CALL GC_CheckVar( 'State_Met%IsLand', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IsLand => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IsIce ) ) THEN
       DEALLOCATE( State_Met%IsIce, STAT=RC )
       CALL GC_CheckVar( 'State_Met%IsIce', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IsIce => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IsSnow ) ) THEN
       DEALLOCATE( State_Met%IsSnow, STAT=RC )
       CALL GC_CheckVar( 'State_Met%IsSnow', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IsSnow => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%LAI ) ) THEN
       DEALLOCATE( State_Met%LAI, STAT=RC )
       CALL GC_CheckVar( 'State_Met%LAI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%LAI => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PARDR ) ) THEN
       DEALLOCATE( State_Met%PARDR, STAT=RC )
       CALL GC_CheckVar( 'State_Met%PARDR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PARDR => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PARDF ) ) THEN
       DEALLOCATE( State_Met%PARDF, STAT=RC )
       CALL GC_CheckVar( 'State_Met%PARDF', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PARDF => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PBLH ) ) THEN
       DEALLOCATE( State_Met%PBLH, STAT=RC )
       CALL GC_CheckVar( 'State_Met%PBLH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBLH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PBL_TOP_hPa ) ) THEN
       DEALLOCATE( State_Met%PBL_TOP_hPa, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PBL_TOP_hPa', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBL_TOP_hPa => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PBL_TOP_L ) ) THEN
       DEALLOCATE( State_Met%PBL_TOP_L, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PBL_TOP_L', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBL_TOP_L => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PBL_TOP_m ) ) THEN
       DEALLOCATE( State_Met%PBL_TOP_m, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PBL_TOP_m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBL_TOP_m => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PBL_THICK ) ) THEN
       DEALLOCATE( State_Met%PBL_THICK, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PBL_THICK', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBL_THICK => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PHIS ) ) THEN
       DEALLOCATE( State_Met%PHIS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PHIS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PHIS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PRECANV ) ) THEN
       DEALLOCATE( State_Met%PRECANV, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PRECANV', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PRECANV => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PRECCON ) ) THEN
       DEALLOCATE( State_Met%PRECCON, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PRECCON', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PRECCON => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PRECLSC ) ) THEN
       DEALLOCATE( State_Met%PRECLSC, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PRECLSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PRECLSC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PRECTOT ) ) THEN
       DEALLOCATE( State_Met%PRECTOT, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PRECTOT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PRECTOT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PS1_WET ) ) THEN
       DEALLOCATE( State_Met%PS1_WET, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PS1_WET', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PS1_WET => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PS1_DRY ) ) THEN
       DEALLOCATE( State_Met%PS1_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PS1_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PS1_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PS2_WET ) ) THEN
       DEALLOCATE( State_Met%PS2_WET, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PS2_WET', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PS2_WET => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PS2_DRY ) ) THEN
       DEALLOCATE( State_Met%PS2_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PS2_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PS2_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PSC2_WET ) ) THEN
       DEALLOCATE( State_Met%PSC2_WET, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PSC2_WET', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PSC2_WET => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PSC2_DRY ) ) THEN
       DEALLOCATE( State_Met%PSC2_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PSC2_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PSC2_DRY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%QV2M ) ) THEN
       DEALLOCATE( State_Met%QV2M, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%QV2M', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%QV2M => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE00 ) ) THEN
       DEALLOCATE( State_Met%SEAICE00, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE00', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE00 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE10 ) ) THEN
       DEALLOCATE( State_Met%SEAICE10, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE10', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE10 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE20 ) ) THEN
       DEALLOCATE( State_Met%SEAICE20, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE20', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE20 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE30 ) ) THEN
       DEALLOCATE( State_Met%SEAICE30, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE30', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE30 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE40 ) ) THEN
       DEALLOCATE( State_Met%SEAICE40, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE40', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE40 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE50 ) ) THEN
       DEALLOCATE( State_Met%SEAICE50, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE50', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE50 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE60 ) ) THEN
       DEALLOCATE( State_Met%SEAICE60, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE60', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE60 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE70 ) ) THEN
       DEALLOCATE( State_Met%SEAICE70, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE70', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE70 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE80 ) ) THEN
       DEALLOCATE( State_Met%SEAICE80, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE80', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE80 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SEAICE90 ) ) THEN
       DEALLOCATE( State_Met%SEAICE90, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SEAICE90', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SEAICE90 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SLP ) ) THEN
       DEALLOCATE( State_Met%SLP, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SLP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SLP => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SNODP ) ) THEN
       DEALLOCATE( State_Met%SNODP, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SNODP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SNODP => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SNOMAS ) ) THEN
       DEALLOCATE( State_Met%SNOMAS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SNOMAS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SNOMAS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SUNCOS ) ) THEN
       DEALLOCATE( State_Met%SUNCOS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SUNCOS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SUNCOS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SUNCOSmid ) ) THEN
       DEALLOCATE( State_Met%SUNCOSmid, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SUNCOSmid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SUNCOSmid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%SWGDN ) ) THEN
       DEALLOCATE( State_Met%SWGDN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SWGDN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SWGDN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TropLev ) ) THEN
       DEALLOCATE( State_Met%TropLev, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TropLev', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TropLev => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TropHt ) ) THEN
       DEALLOCATE( State_Met%TropHt, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TropHt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TropHt => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TROPP ) ) THEN
       DEALLOCATE( State_Met%TROPP, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TROPP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TROPP => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TS ) ) THEN
       DEALLOCATE( State_Met%TS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TSKIN ) ) THEN
       DEALLOCATE( State_Met%TSKIN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TSKIN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TSKIN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%TO3 ) ) THEN
       DEALLOCATE( State_Met%TO3, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TO3 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%U10M ) ) THEN
       DEALLOCATE( State_Met%U10M, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%U10M', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%U10M => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%USTAR ) ) THEN
       DEALLOCATE( State_Met%USTAR, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%USTAR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%USTAR => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%UVALBEDO ) ) THEN
       DEALLOCATE( State_Met%UVALBEDO, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%UVALBEDO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%UVALBEDO => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%V10M ) ) THEN
       DEALLOCATE( State_Met%V10M, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%V10M', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%V10M => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%Z0 ) ) THEN
       DEALLOCATE( State_Met%Z0, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%Z0', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%Z0 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CNV_FRC ) ) THEN
       DEALLOCATE( State_Met%CNV_FRC, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%CNV_FRC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CNV_FRC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%ILAND ) ) THEN
       DEALLOCATE( State_Met%ILAND, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%ILAND', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%ILAND => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IREG ) ) THEN
       DEALLOCATE( State_Met%IREG, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%IREG', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IREG => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IMIX) ) THEN
       DEALLOCATE( State_Met%IMIX, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%IMIX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IMIX => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%FPBL ) ) THEN
       DEALLOCATE( State_Met%FPBL, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%FPBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%FPBL => NULL()
    ENDIF

    !========================================================================
    ! Deallocate 3-D fields
    !
    ! NOTE: If using GEOS-Chem as GCHP, or coupled to GMAO/GEOS, then just
    ! nullify the fields without deallocating.  This will prevent abnormal
    ! exits in MAPL.  This is probably due to the fact that the State_Met
    ! fields point to ESMF/MAPL Imports, and cannot be deallocated
    ! before the Import itself is finalized.
    !
    ! ALSO NOTE: If using GEOS-Chem coupled to WRF, then do the same
    ! as for GCHP or GEOS, as WRF does its own separate deallocation.
    !
    !  -- Lizzie Lundgren and Bob Yantosca, 05 Nov 2018
    !========================================================================
    IF ( ASSOCIATED( State_Met%IUSE ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%IUSE => NULL()
#else
       DEALLOCATE( State_Met%IUSE, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%IUSE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IUSE => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%LANDTYPEFRAC ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%LANDTYPEFRAC => NULL()
#else
       DEALLOCATE( State_Met%LANDTYPEFRAC, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%LANDTYPEFRAC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%LANDTYPEFRAC => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%MODISLAI ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%MODISLAI => NULL()
#else
       DEALLOCATE( State_Met%MODISLAI, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%MODISLAI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%MODISLAI => NULL()
#endif
    ENDIF

    !------------------------------------------------------------------------
    ! 3-D fields
    !------------------------------------------------------------------------
    IF ( ASSOCIATED( State_Met%AD ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%AD => NULL()
#else
       DEALLOCATE( State_Met%AD, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%AD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AD => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%AIRDEN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%AIRDEN => NULL()
#else
       DEALLOCATE( State_Met%AIRDEN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%AIRDEN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AIRDEN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%MAIRDEN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%MAIRDEN => NULL()
#else
       DEALLOCATE( State_Met%MAIRDEN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%MAIRDEN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%MAIRDEN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%AIRVOL ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%AIRVOL => NULL()
#else
       DEALLOCATE( State_Met%AIRVOL, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%AIRVOL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AIRVOL => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%BXHEIGHT ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%BXHEIGHT => NULL()
#else
       DEALLOCATE( State_Met%BXHEIGHT, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%BXHEIGHT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%BXHEIGHT => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%CLDF ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%CLDF => NULL()
#else
       DEALLOCATE( State_Met%CLDF, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%CLDF', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CLDF => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%CMFMC ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%CMFMC => NULL()
#else
       DEALLOCATE( State_Met%CMFMC, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%CMFMC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%CMFMC => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DELP ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DELP => NULL()
#else
       DEALLOCATE( State_Met%DELP, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DELP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DELP => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DELP_DRY ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DELP_DRY => NULL()
#else
       DEALLOCATE( State_Met%DELP_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DELP_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DELP_DRY => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DP_DRY_PREV ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DP_DRY_PREV => NULL()
#else
       DEALLOCATE( State_Met%DP_DRY_PREV, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DP_DRY_PREV', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DP_DRY_PREV => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DQRCU ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DQRCU => NULL()
#else
       DEALLOCATE( State_Met%DQRCU, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DQRCU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DQRCU => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DQRLSAN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DQRLSAN => NULL()
#else
       DEALLOCATE( State_Met%DQRLSAN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DQRLSAN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DQRLSAN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%DTRAIN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%DTRAIN => NULL()
#else
       DEALLOCATE( State_Met%DTRAIN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%DTRAIN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%DTRAIN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%F_OF_PBL ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%F_OF_PBL => NULL()
#else
       DEALLOCATE( State_Met%F_OF_PBL, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%F_OF_PBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%F_OF_PBL => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%F_UNDER_PBLTOP ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%F_UNDER_PBLTOP => NULL()
#else
       DEALLOCATE( State_Met%F_UNDER_PBLTOP, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%F_UNDER_PBLTOP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%F_UNDER_PBLTOP => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%OMEGA ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%OMEGA => NULL()
#else
       DEALLOCATE( State_Met%OMEGA, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%OMEGA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%OMEGA => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%OPTD ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%OPTD => NULL()
#else
       DEALLOCATE( State_Met%OPTD, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%OPTD', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%OPTD => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PEDGE ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PEDGE => NULL()
#else
       DEALLOCATE( State_Met%PEDGE, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PEDGE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PEDGE => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PEDGE_DRY ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PEDGE_DRY => NULL()
#else
       DEALLOCATE( State_Met%PEDGE_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PEDGE_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PEDGE_DRY => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PFICU ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PFICU => NULL()
#else
       DEALLOCATE( State_Met%PFICU, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PFICU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PFICU => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PFILSAN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PFILSAN => NULL()
#else
       DEALLOCATE( State_Met%PFILSAN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PFILSAN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PFILSAN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PFLCU ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PFLCU => NULL()
#else
       DEALLOCATE( State_Met%PFLCU, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PFLCU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PFLCU => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PFLLSAN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PFLLSAN => NULL()
#else
       DEALLOCATE( State_Met%PFLLSAN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PFLLSAN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PFLLSAN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PMID ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PMID => NULL()
#else
       DEALLOCATE( State_Met%PMID, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PMID', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PMID => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%PMID_DRY ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%PMID_DRY => NULL()
#else
       DEALLOCATE( State_Met%PMID_DRY, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PMID_DRY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PMID_DRY => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%QI ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%QI => NULL()
#else
       DEALLOCATE( State_Met%QI, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%QI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%QI => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%QL ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%QL => NULL()
#else
       DEALLOCATE( State_Met%QL, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%QL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%QL => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%REEVAPCN ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%REEVAPCN => NULL()
#else
       DEALLOCATE( State_Met%REEVAPCN, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%REEVAPCN', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%REEVAPCN => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%REEVAPLS ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%REEVAPLS => NULL()
#else
       DEALLOCATE( State_Met%REEVAPLS, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%REEVAPLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%REEVAPLS => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%RH ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%RH => NULL()
#else
       DEALLOCATE( State_Met%RH, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%RH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%RH => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%SPHU ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%SPHU => NULL()
#else
       DEALLOCATE( State_Met%SPHU, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SPHU', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SPHU => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%SPHU1 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%SPHU1 => NULL()
#else
       DEALLOCATE( State_Met%SPHU1, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SPHU1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SPHU1 => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%SPHU2 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%SPHU2 => NULL()
#else
       DEALLOCATE( State_Met%SPHU2, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%SPHU2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%SPHU2 => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%T ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%T => NULL()
#else
       DEALLOCATE( State_Met%T, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%T', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%T => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%TMPU1 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%TMPU1 => NULL()
#else
       DEALLOCATE( State_Met%TMPU1, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TMPU1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TMPU1 => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%TMPU2 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%TMPU2 => NULL()
#else
       DEALLOCATE( State_Met%TMPU2, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TMPU2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TMPU2 => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%TV ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%TV => NULL()
#else
       DEALLOCATE( State_Met%TV, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TV', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TV => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%TAUCLI ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%TAUCLI => NULL()
#else
       DEALLOCATE( State_Met%TAUCLI, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TAUCLI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TAUCLI => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%TAUCLW ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%TAUCLW => NULL()
#else
       DEALLOCATE( State_Met%TAUCLW, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%TAUCLW', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%TAUCLW => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%U ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%U => NULL()
#else
       DEALLOCATE( State_Met%U, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%U', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%U => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%UPDVVEL ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%UPDVVEL => NULL()
#else
       DEALLOCATE( State_Met%UPDVVEL, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%UPDVVEL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%UPDVVEL => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%V ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%V => NULL()
#else
       DEALLOCATE( State_Met%V, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%V', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%V => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%XLAI ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XLAI => NULL()
#else
       DEALLOCATE( State_Met%XLAI, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XLAI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XLAI => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%XLAI2 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XLAI2 => NULL()
#else
       DEALLOCATE( State_Met%XLAI2, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XLAI2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XLAI2 => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%XLAI_NATIVE ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XLAI_NATIVE => NULL()
#else
       DEALLOCATE( State_Met%XLAI_NATIVE, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XLAI_NATIVE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XLAI_NATIVE => NULL()
#endif
    ENDIF

    !=======================================================================
    ! Fields for querying which vertical regime a grid box is in
    ! or if it is near local solar noon at a grid box
    !=======================================================================
    IF ( ASSOCIATED( State_Met%InChemGrid ) ) THEN
       DEALLOCATE( State_Met%InChemGrid, STAT=RC )
       CALL GC_CheckVar( 'State_Met%InChemGrid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%InChemGrid => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%InPbl ) ) THEN
       DEALLOCATE( State_Met%InPbl, STAT=RC )
       CALL GC_CheckVar( 'State_Met%InPbl', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%InPbl => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%InStratMeso ) ) THEN
       DEALLOCATE( State_Met%InStratMeso, STAT=RC )
       CALL GC_CheckVar( 'State_Met%InStratMeso', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%InStratMeso => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%InStratosphere ) ) THEN
       DEALLOCATE( State_Met%InStratosphere, STAT=RC )
       CALL GC_CheckVar( 'State_Met%InStratosphere', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%InStratosphere => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%InTroposphere ) ) THEN
       DEALLOCATE( State_Met%InTroposphere, STAT=RC )
       CALL GC_CheckVar( 'State_Met%InTroposphere', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%InTroposphere => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%IsLocalNoon ) ) THEN
       DEALLOCATE( State_Met%IsLocalNoon, STAT=RC )
       CALL GC_CheckVar( 'State_Met%IsLocalNoon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%IsLocalNoon => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%LocalSolarTime ) ) THEN
       DEALLOCATE( State_Met%LocalSolarTime, STAT=RC )
       CALL GC_CheckVar( 'State_Met%LocalSolarTime', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%LocalSolarTime => NULL()
    ENDIF

    !========================================================================
    ! Fields temporaries used in other modules
    !========================================================================
    IF ( ASSOCIATED( State_Met%C_H2O ) ) THEN
      DEALLOCATE( State_Met%C_H2O, STAT=RC )
      CALL GC_CheckVar( 'State_Met%C_H2O', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%C_H2O => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CLDICE ) ) THEN
      DEALLOCATE( State_Met%CLDICE, STAT=RC )
      CALL GC_CheckVar( 'State_Met%CLDICE', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%CLDICE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%CLDLIQ ) ) THEN
      DEALLOCATE( State_Met%CLDLIQ, STAT=RC )
      CALL GC_CheckVar( 'State_Met%CLDLIQ', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%CLDLIQ => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%PDOWN ) ) THEN
      DEALLOCATE( State_Met%PDOWN, STAT=RC )
      CALL GC_CheckVar( 'State_Met%PDOWN', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%PDOWN => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%QQ ) ) THEN
      DEALLOCATE( State_Met%QQ, STAT=RC )
      CALL GC_CheckVar( 'State_Met%QQ', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%QQ => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%REEVAP ) ) THEN
      DEALLOCATE( State_Met%REEVAP, STAT=RC )
      CALL GC_CheckVar( 'State_Met%REEVAP', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%REEVAP => NULL()
    ENDIF

#ifdef APM
    IF ( ASSOCIATED( State_Met%PSO4_SO2APM2 ) ) THEN
      DEALLOCATE( State_Met%PSO4_SO2APM2, STAT=RC )
      CALL GC_CheckVar( 'State_Met%PSO4_SO2APM2', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      State_Met%PSO4_SO2APM2 => NULL()
    ENDIF
#endif

   !-------------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
   !-------------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Met%xxx ) ) THEN
    !   DEALLOCATE( State_Met%xxx, STAT=RC )
    !   CALL GC_CheckVar( 'State_Met%xxx', 2, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Met%xxx => NULL()
    !ENDIF

    !========================================================================
    ! Destroy the registry of fields for this module
    !========================================================================
    CALL Registry_Destroy( State_Met%Registry, State_Met%RegDict, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Met%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Nullify the registry object
    State_Met%Registry => NULL()

  END SUBROUTINE Cleanup_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Metadata_State_Met
!
! !DESCRIPTION: Subroutine GET\_METDATA\_STATE\_MET retrieves basic
!  information about each State\_Met field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Met( am_I_Root, metadataID, Found,   RC,     &
                                     Desc,      Units,      Rank,    Type,   &
                                     VLoc,      perQnt                      )
!
! !USES:
!
    USE Charpak_Mod,         ONLY: To_UpperCase
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)  :: metadataID ! State_Met field ID
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT) :: Found      ! Item found?
    INTEGER,             INTENT(OUT) :: RC         ! Return code
    CHARACTER(LEN=255),  OPTIONAL    :: Desc       ! Long name string
    CHARACTER(LEN=255),  OPTIONAL    :: Units      ! Units string
    INTEGER,             OPTIONAL    :: Rank       ! # of dimensions
    INTEGER,             OPTIONAL    :: Type       ! Desc of data type
    INTEGER,             OPTIONAL    :: VLoc       ! Vertical placement
    CHARACTER(LEN=255),  OPTIONAL    :: perQnt     ! "Quantity" dimension?
!
! !REVISION HISTORY:
!  28 Aug 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Name_AllCaps
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc, isQnt

    !========================================================================
    ! Initialize
    !========================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Met (in Headers/state_met_mod.F90)'
    Found = .TRUE.

    ! Optional arguments present?
    isDesc  = PRESENT( Desc   )
    isUnits = PRESENT( Units  )
    isRank  = PRESENT( Rank   )
    isType  = PRESENT( Type   )
    isVLoc  = PRESENT( VLoc   )
    isQnt   = PRESENT( perQnt )

    ! Set defaults for optional arguments. Assume type and vertical
    ! location are real (flexible precision) and center unless specified
    ! otherwise
    IF ( isUnits ) Units  = ''
    IF ( isDesc  ) Desc   = ''
    IF ( isRank  ) Rank   = -1              ! initialize as bad value
    IF ( isType  ) Type   = KINDVAL_FP      ! Assume real with flex precision
    IF ( isVLoc  ) VLoc   = VLocationNone   ! Assume no vertical location
    IF ( isQnt   ) perQnt = ''              ! Assume no "species" dimension

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !========================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !========================================================================
    SELECT CASE ( TRIM( Name_AllCaps) )

       !---------------------------------------------------------------------
       ! 2-D Fields
       !---------------------------------------------------------------------
       CASE ( 'ALBD' )
          IF ( isDesc  ) Desc  = 'Visible surface albedo'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'AREAM2' )
          IF ( isDesc  ) Desc  = 'Surface area of grid box'
          IF ( isUnits ) Units = 'm2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'CHEMGRIDLEV' )
          IF ( isDesc  ) Desc  = 'Highest level of the chemistry grid'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'CLDFRC' )
          IF ( isDesc  ) Desc  = 'Column cloud fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'CLDTOPS' )
          IF ( isDesc  ) Desc  = 'Maximum cloud top height'
          IF ( isUnits ) Units = 'level'
          IF ( isRank  ) Rank  = 2
          IF ( isType  ) Type  = KINDVAL_I4

#ifdef MODEL_GEOS
       CASE ( 'CNVFRC' )
          IF ( isDesc  ) Desc  = 'Convective fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2
#endif

       CASE ( 'CONVDEPTH' )
          IF ( isDesc  ) Desc  = 'Convective cloud depth'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 2

       CASE ( 'EFLUX' )
          IF ( isDesc  ) Desc  = 'Latent heat flux'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FLASHDENS' )
          IF ( isDesc  ) Desc  = 'Lightning flash density'
          IF ( isUnits ) Units = 'km-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRCLND' )
          IF ( isDesc  ) Desc  = 'Olson land fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRLAKE' )
          IF ( isDesc  ) Desc  = 'Fraction of lake'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRLAND' )
          IF ( isDesc  ) Desc  = 'Fraction of land'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRLANDIC' )
          IF ( isDesc  ) Desc  = 'Fraction of land ice'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FROCEAN' )
          IF ( isUnits ) Units = '1'
          IF ( isDesc  ) Desc  = 'Fraction of ocean'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRSEAICE' )
          IF ( isDesc  ) Desc  = 'Fraction of sea ice'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRSNO' )
          IF ( isDesc  ) Desc  = 'Fraction of snow on surface'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'GWETROOT' )
          IF ( isDesc  ) Desc  = 'Root soil wetness'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'GWETTOP' )
          IF ( isDesc  ) Desc  = 'Top soil moisture'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'HFLUX' )
          IF ( isDesc  ) Desc  = 'Sensible heat flux'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'LAI' )
          IF ( isDesc  ) Desc  = 'Leaf area index from GMAO'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PARDR' )
          IF ( isDesc  ) Desc  = 'Direct photosynthetically-active radiation'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PARDF' )
          IF ( isDesc  ) Desc  = 'Diffuse photosynthetically-active radiation'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PBLH' )
          IF ( isDesc  ) Desc  = 'Planetary boundary layer height'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PBLTOPHPA' )
          IF ( isDesc  ) Desc  = 'Planetary boundary layer top'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PBLTOPL' )
          IF ( isDesc  ) Desc  = 'Planetary boundary layer top'
          IF ( isUnits ) Units = 'layer'
          IF ( isRank  ) Rank  = 2
          IF ( isType  ) Type  = KINDVAL_I4

       CASE ( 'PBLTOPM' )
          IF ( isDesc  ) Desc  = 'Planetary boundary layer top'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PBLTHICK' )
          IF ( isDesc  ) Desc  = 'Planetary boundary layer thickness'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PHIS' )
          IF ( isDesc  ) Desc  = 'Surface geopotential height'
          IF ( isUnits ) Units = 'm2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECANV' )
          IF ( isDesc  ) Desc  = 'Anvil precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECCON' )
          IF ( isDesc  ) Desc  = 'Convective precipitation at the ground'
          IF ( isUnits ) Units = 'mm day-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECLSC' )
          IF ( isDesc  ) Desc  = 'Large-scale precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECTOT' )
          IF ( isDesc  ) Desc  = 'Total precipitation at the ground'
          IF ( isUnits ) Units = 'mm day-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS1WET' )
          IF ( isDesc  ) Desc  = 'Wet surface pressure at dt start'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS2WET' )
          IF ( isDesc  ) Desc  = 'Wet surface pressure at dt end'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PSC2WET' )
          IF ( isDesc  ) Desc  = 'Wet interpolated surface pressure'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS1DRY' )
          IF ( isDesc  ) Desc  = 'Dry surface pressure at dt start'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS2DRY' )
          IF ( isDesc  ) Desc  = 'Dry surface pressure at dt end'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PSC2DRY' )
          IF ( isDesc  ) Desc  = 'Dry interpolated surface pressure'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE00' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 00-10%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'QV2M' )
          IF ( isDesc  ) Desc  = 'Specific humidity at 2 m'
          IF ( isUnits ) Units = 'kg kg-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE10' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 10-20%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE20' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 20-30%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE30' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 30-40%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE40' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 40-50%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE50' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 50-60%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE60' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 60-70%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE70' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 70-80%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE80' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 80-90%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE90' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 90-100%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SLP' )
          IF ( isDesc  ) Desc  = 'Sea level pressure'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SNODP' )
          IF ( isDesc  ) Desc  = 'Snow depth'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SNOMAS' )
          IF ( isDesc  ) Desc  = 'Snow mass'
          IF ( isUnits ) Units = 'kg m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SST' )
          IF ( isDesc  ) Desc  = 'Sea surface temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SUNCOS' )
          IF ( isDesc  ) Desc  = 'Cosine of solar zenith angle, current time'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SUNCOSMID' )
          IF ( isDesc  ) Desc  = 'Cosine of solar zenith angle, at ' // &
                                 'midpoint of chemistry timestep'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SUNCOSSUM' )
          IF ( isDesc  ) Desc  = 'Sum of Cosine of solar zenith angle, current time (HEMCO)'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SZAFACT' )
          IF ( isDesc  ) Desc  = 'Diurnal scale factor from dividing the sza by the sum of the total sza per day (HEMCO)'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SWGDN' )
          IF ( isDesc  ) Desc  = 'Incident shortwave radiation at ground'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TO3' )
          IF ( isDesc  ) Desc  = 'Total overhead ozone column'
          IF ( isUnits ) Units = 'dobsons'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TROPLEV' )
          IF ( isDesc  ) Desc  = 'GEOS-Chem level where the tropopause occurs'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TROPHT' )
          IF ( isDesc  ) Desc  = 'Tropopause height'
          IF ( isUnits ) Units = 'km'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TROPP' )
          IF ( isDesc  ) Desc  = 'Tropopause pressure'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TS' )
          IF ( isDesc  ) Desc  = 'Surface temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 2

       CASE ( 'TSKIN' )
          IF ( isDesc  ) Desc  = 'Surface skin temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 2

       CASE ( 'U10M' )
          IF ( isDesc  ) Desc  = 'East-west wind at 10 meter height'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'USTAR' )
          IF ( isDesc  ) Desc  = 'Friction velocity'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'UVALBEDO' )
          IF ( isDesc  ) Desc  = 'Ultraviolet surface albedo'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'V10M' )
          IF ( isDesc  ) Desc  = 'North-south wind at 10 meter height'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'Z0' )
          IF ( isDesc  ) Desc  = 'Surface roughness height'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 2

       CASE ( 'LOCALSOLARTIME' )
          IF ( isDesc  ) Desc  = 'Local solar time'
          IF ( isUnits ) Units = 'hours'
          IF ( isRank  ) Rank  = 2

       !---------------------------------------------------------------------
       ! 3-D Fields
       !---------------------------------------------------------------------
       CASE ( 'AD' )
          IF ( isDesc  ) Desc  = 'Dry air mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AIRDEN' )
          IF ( isDesc  ) Desc  = 'Dry air density'
          IF ( isUnits ) Units = 'kg m-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'MAIRDEN' )
          IF ( isDesc  ) Desc  = 'Moist air density'
          IF ( isUnits ) Units = 'kg m-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AIRNUMDEN' )
          IF ( isDesc  ) Desc  = 'Dry air density'
          IF ( isUnits ) Units = 'molec cm-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AIRVOL' )
          IF ( isDesc  ) Desc  = 'Volume of grid box'
          IF ( isUnits ) Units = 'm3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AVGW' )
          IF ( isDesc  ) Desc  = 'Water vapor mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'vol vol-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'BXHEIGHT' )
          IF ( isDesc  ) Desc  = 'Grid box height'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'CLDF' )
          IF ( isDesc  ) Desc  = '3-D cloud fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'CMFMC' )
          IF ( isDesc  ) Desc  = 'Cloud mass flux'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'DELP' )
          IF ( isDesc  ) Desc  = 'Delta-pressure across grid box(wet air)'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'DELPDRY' )
          IF ( isDesc  ) Desc  = 'Delta-pressure across grid box (dry air)'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'DPDRYPREV' )
          IF ( isDesc  ) Desc  = 'Previous State_Met%DELP_DRY'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'SPHUPREV' )
          IF ( isDesc  ) Desc  = 'Previous State_Met%SPHU_PREV'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'DQRCU' )
          IF ( isDesc  ) Desc  = 'Production rate of convective ' // &
                                 'precipitation (per dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'DQRLSAN' )
          IF ( isDesc  ) Desc  = 'Production rate of large-scale ' // &
                                 'precipitation (per dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'DTRAIN' )
          IF ( isDesc  ) Desc  = 'Detrainment flux'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'FOFPBL' )
          IF ( isDesc  ) Desc  = 'Fraction of PBL'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'FUNDERPBLTOP' )
          IF ( isDesc  ) Desc  = 'Fraction of box under PBL top'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'OMEGA' )
          IF ( isDesc  ) Desc  = 'Updraft velocity'
          IF ( isUnits ) Units = 'Pa s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'OPTD' )
          IF ( isDesc  ) Desc  = 'Visible optical depth'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'PEDGE' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t moist air) at level edges'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PEDGEDRY' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t dry air) at level edges'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PFICU' )
          IF ( isDesc  ) Desc  = 'Downward flux of ice precipitation ' // &
                                 '(convective)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PFILSAN' )
          IF ( isDesc  ) Desc  = 'Downwared flux of ice precipitation ' // &
                                 '(large-scale + anvil)'
          IF ( isRank  ) Rank  = 3
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PFLCU' )
          IF ( isDesc  ) Desc  = 'Downward flux of liquid precipitation ' // &
                                 '(convective)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PFLLSAN' )
          IF ( isDesc  ) Desc  = 'Downward flux of liquid precipitation ' // &
                                 '(large-scale + anvil)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationEdge

       CASE ( 'PMID' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t moist air) at level centers'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'PMIDDRY' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t dry air) at level centers'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'QI' )
          IF ( isDesc  ) Desc  = 'Ice mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'QL' )
          IF ( isDesc  ) Desc  = 'Water mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'REEVAPCN' )
          IF ( isDesc  ) Desc  = 'Evaporation of convective ' // &
                                 'precipitation (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'REEVAPLS' )
          IF ( isDesc  ) Desc  = 'Evaporation of large-scale + anvil ' // &
                                 'precipitation (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg '
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'RH' )
          IF ( isDesc  ) Desc  = 'Relative humidity'
          IF ( isUnits ) Units = '%'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'SPHU' )
          IF ( isDesc  ) Desc  = 'Specific humidity (w/r/t moist air)'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'SPHU1' )
          IF ( isDesc  ) Desc  = 'Instantaneous specific humidity at time=T'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'SPHU2' )
          IF ( isDesc  ) Desc  = 'Instantaneous specific humidity at time=T+dt'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'T' )
          IF ( isDesc  ) Desc  = 'Temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'THETA' )
          IF ( isDesc  ) Desc  = 'Potential temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'TV' )
          IF ( isDesc  ) Desc  = 'Virtual temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'TAUCLI' )
          IF ( isDesc  ) Desc  = 'Optical depth of ice clouds'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'TAUCLW' )
          IF ( isDesc  ) Desc  = 'Optical depth of H2O clouds'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'TMPU1' )
          IF ( isDesc  ) Desc  = 'Instantaneous temperature at time=T'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'TMPU2' )
          IF ( isDesc  ) Desc  = 'Instantaneous temperature at time T+dt'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'U' )
          IF ( isDesc  ) Desc  = 'East-west component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

#ifdef MODEL_GEOS
       CASE ( 'UPDVVEL' )
          IF ( isDesc  ) Desc  = 'Updraft vertical velocity'
          IF ( isUnits ) Units = 'hPa s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter
#endif

       CASE ( 'V' )
          IF ( isDesc  ) Desc  = 'North-south component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       !---------------------------------------------------------------------
       ! Offline land type, leaf area index, and chlorophyll fields
       !---------------------------------------------------------------------
       CASE ( 'IREG' )
          IF ( isDesc  ) Desc   = 'Number of Olson land types in each grid box'
          IF ( isUnits ) Units  = '1'
          IF ( isRank  ) Rank   = 2
          IF ( isType  ) Type   = KINDVAL_I4

       CASE ( 'ILAND' )
          IF ( isDesc  ) Desc   = 'Olson land type indices in each grid box'
          IF ( isUnits ) Units  = '1'
          IF ( isRank  ) Rank   = 2
          IF ( isType  ) Type   = KINDVAL_I4
          IF ( isQnt   ) perQnt = 'OLSON'

       CASE ( 'IUSE' )
          IF ( isDesc  ) Desc   = 'Fraction (per mil) occupied by each '  // &
                                   'Olson land type in the grid box'
          IF ( isUnits ) Units  = 'o/oo'
          IF ( isRank  ) Rank   = 2
          IF ( isType  ) Type   = KINDVAL_I4
          IF ( isQnt   ) perQnt = 'OLSON'

       CASE ( 'XLAI' )
          IF ( isDesc  ) Desc   = 'MODIS LAI for each Olson land type, '  // &
                                  'current month'
          IF ( isUnits ) Units  = 'm2 m-2'
          IF ( isRank  ) Rank   = 2
          IF ( isQnt   ) perQnt = 'OLSON'

       CASE ( 'XLAI2' )
          IF ( isDesc  ) Desc   = 'MODIS LAI for each Olson land type, '  // &
                                   'next month'
          IF ( isUnits ) Units  = 'm2 m-2'
          IF ( isRank  ) Rank   = 2
          IF ( isQnt   ) perQnt = 'OLSON'

       CASE ( 'MODISLAI' )
          IF ( isDesc  ) Desc   = 'Daily LAI computed from monthly '      // &
                                 'offline MODIS values'
          IF ( isUnits ) Units  = 'm2 m-2'
          IF ( isRank  ) Rank   = 2

       CASE ( 'LANDTYPEFRAC' )
          IF ( isDesc  ) Desc   = 'Olson fraction per land type'
          IF ( isUnits ) Units  = '1'
          IF ( isRank  ) Rank   = 2
          IF ( isQnt   ) perQnt = 'OLSON'

       CASE ( 'XLAINATIVE' )
          IF ( isDesc  ) Desc   = 'Average LAI per Olson land type'
          IF ( isUnits ) Units  = 'm2 m-2'
          IF ( isRank  ) Rank   = 2
          IF ( isQnt   ) perQnt = 'OLSON'

!       CASE ( 'INCHEMGRID' )
!          IF ( isDesc  ) Desc  = 'Is each grid box in the chemistry grid?'
!          IF ( isUnits ) Units = 'boolean'
!          IF ( isRank  ) Rank  = 3
!
!       CASE ( 'INTROPOSPHERE' )
!          IF ( isDesc  ) Desc  = 'Is each grid box in the troposphere?'
!          IF ( isUnits ) Units = 'boolean'
!          IF ( isRank  ) Rank  = 3
!
!       CASE ( 'INPBL' )
!          IF ( isDesc  ) Desc  = 'Is each grid box in the planetary boundary layer?'
!          IF ( isUnits ) Units = 'boolean'
!          IF ( isRank  ) Rank  = 3

       CASE DEFAULT
          Found = .False.
          ErrMsg = 'Metadata not found for State_Met field ID: ' &
                   // TRIM( metadataID )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

    END SELECT

    ! Set VLoc to undefined if variable is 2d
    IF ( isVLoc .AND. Rank == 2 ) THEN
       VLoc = VLocationNone
    ENDIF

   END SUBROUTINE Get_Metadata_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_2D
!
! !DESCRIPTION: Allocates the data array for a State_Met field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 2-dimensional array fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_2D( Input_Opt, State_Met, State_Grid,      &
                                      Ptr2Data,  metId,     RC,              &
                                      noRegister                            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doRegister

    ! Arrays
    INTEGER            :: dims(2)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       dims(1) = State_Grid%NX
       dims(2) = State_Grid%NY

       ! Allocate the array
       ALLOCATE( Ptr2Data( dims(1), dims(2) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f4

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_3D
!
! !DESCRIPTION: Allocates the data array for a State_Met field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 3-dimensional array fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_3D( Input_Opt,  State_Met, State_Grid,     &
                                      Ptr2Data,   metId,     RC,             &
                                      noRegister, onEdges,   zxyOrder,       &
                                      nSlots                                )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
    LOGICAL,          OPTIONAL    :: onEdges             ! Data on vert edges?
    LOGICAL,          OPTIONAL    :: zxyOrder            ! Data array (Z,X,Y)?
    INTEGER,          OPTIONAL    :: nSlots              ! # slots for Z dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doEdges, doRegister, doSlots, doZxy

    ! Arrays
    INTEGER            :: dims(3)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_3D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( onEdges ) ) THEN
       doEdges = onEdges
    ELSE
       doEdges = .FALSE.
    ENDIF

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    doSlots = PRESENT( nSlots )

    IF ( PRESENT( zxyOrder ) ) THEN
       doZxy = zxyOrder
    ELSE
       doZxy = .FALSE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array ID and dimensions
       IF ( doZxy ) THEN

          ! ZXY order
          dims(1) = State_Grid%NZ
          dims(2) = State_Grid%NX
          dims(3) = State_Grid%NY

          ! If we have specified nSlots, use that for 1st dimension
          ! Otherwise, if data is on vertical edges, increment 1st dimension
          IF ( doSlots ) THEN
             dims(1) = nSlots
          ELSE
             IF ( doEdges ) dims(1) = dims(1) + 1
          ENDIF

       ELSE

          ! XYZ order
          dims(1) = State_Grid%NX
          dims(2) = State_Grid%NY
          dims(3) = State_Grid%NZ

          ! If we have specified nSlots, use that for 3rd dimension
          ! Otherwise, if data is on vertical edges, increment 3rd dimension
          IF ( doSlots ) THEN
             dims(3) = nSlots
          ELSE
             IF ( doEdges ) dims(3) = dims(3) + 1
          ENDIF

       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( dims(1), dims(2), dims(3) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f4

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_2D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 8-byte, 2-dimensional fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_2D( Input_Opt, State_Met, State_Grid,      &
                                      Ptr2Data,  metId,     RC,              &
                                      noRegister                            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),         POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doRegister

    ! Arrays
    INTEGER            :: dims(2)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !=======================================================================
    ! Init_and_Register_R8_2D begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !=======================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !=======================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       dims(1) = State_Grid%NX
       dims(2) = State_Grid%NY

       ! Allocate the data
       ALLOCATE( Ptr2Data( dims(1), dims(2) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f8

    ENDIF

    !=======================================================================
    ! Register the field
    !=======================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_3D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 8-byte, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_3D( Input_Opt,  State_Met, State_Grid,     &
                                      Ptr2Data,   metId,     RC,             &
                                      noRegister, onEdges,   zxyOrder,       &
                                      nSlots                                )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
    LOGICAL,          OPTIONAL    :: onEdges             ! Data on vert edges?
    LOGICAL,          OPTIONAL    :: zxyOrder            ! Data array (Z,X,Y)?
    INTEGER,          OPTIONAL    :: nSlots              ! # slots, Z dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doEdges, doRegister, doSlots, doZxy

    ! Arrays
    INTEGER            :: dims(3)


    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( onEdges ) ) THEN
       doEdges = onEdges
    ELSE
       doEdges = .FALSE.
    ENDIF

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    doSlots = PRESENT( nSlots )

    IF ( PRESENT( zxyOrder ) ) THEN
       doZxy = zxyOrder
    ELSE
       doZxy = .FALSE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array ID and dimensions
       IF ( doZxy ) THEN

          ! ZXY order
          dims(1) = State_Grid%NZ
          dims(2) = State_Grid%NX
          dims(3) = State_Grid%NY

          ! If we have specified nSlots, use that for 1st dimension
          ! Otherwise, if data is on vertical edges, increment 1st dimension
          IF ( doSlots ) THEN
             dims(1) = nSlots
          ELSE
             IF ( doEdges ) dims(1) = dims(1) + 1
          ENDIF

       ELSE

          ! XYZ order
          dims(1) = State_Grid%NX
          dims(2) = State_Grid%NY
          dims(3) = State_Grid%NZ

          ! If we have specified nSlots, use that for 3rd dimension
          ! Otherwise, if data is on vertical edges, increment 3rd dimension
          IF ( doSlots ) THEN
             dims(3) = nSlots
          ELSE
             IF ( doEdges ) dims(3) = dims(3) + 1
          ENDIF

       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( dims(1), dims(2), dims(3) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f8

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_Log_2D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field.
!  This particular routine is for logical, 2-dimensional arrays.
!  NOTE: At present, it is not possible to archive logical fields
!  to HISTORY diagnostics, so we will skip registering logical fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_Log_2D( Input_Opt, State_Met, State_Grid,     &
                                       Ptr2Data,  metId,     RC,             &
                                       noRegister                           )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_Log_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    !========================================================================
    ! Allocate the field array
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY

       ALLOCATE( Ptr2Data( NX, NY ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = .FALSE.
    ENDIF

  END SUBROUTINE Init_and_Register_Log_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_Log_3D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field.
!  This particular routine is for logical, 2-dimensional arrays.
!  NOTE: At present, it is not possible to archive logical fields
!  to HISTORY diagnostics, so we will skip registering logical fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_Log_3D( Input_Opt,  State_Met, State_Grid,    &
                                       Ptr2Data,   metId,     RC,            &
                                       noRegister, onEdges,   zxyOrder,      &
                                       nSlots                               )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
    LOGICAL,          OPTIONAL    :: onEdges             ! Data on vert edges?
    LOGICAL,          OPTIONAL    :: zxyOrder            ! Data array (Z,X,Y)?
    INTEGER,          OPTIONAL    :: nSlots              ! # slots for Z dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure!
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doEdges, doRegister, doSlots, doZxy

    ! Arrays
    INTEGER            :: dims(3)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_Log_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( onEdges ) ) THEN
       doEdges = onEdges
    ELSE
       doEdges = .FALSE.
    ENDIF

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    doSlots = PRESENT( nSlots )

    IF ( PRESENT( zxyOrder ) ) THEN
       doZxy = zxyOrder
    ELSE
       doZxy = .FALSE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array ID and dimensions
       IF ( doZxy ) THEN

          ! ZXY order
          dims(1) = State_Grid%NZ
          dims(2) = State_Grid%NX
          dims(3) = State_Grid%NY

          ! If we have specified nSlots, use that for 1st dimension
          ! Otherwise, if data is on vertical edges, increment 1st dimension
          IF ( doSlots ) THEN
             dims(1) = nSlots
          ELSE
             IF ( doEdges ) dims(1) = dims(1) + 1
          ENDIF

       ELSE

          ! XYZ order
          dims(1) = State_Grid%NX
          dims(2) = State_Grid%NY
          dims(3) = State_Grid%NZ

          ! If we have specified nSlots, use that for 3rd dimension
          ! Otherwise, if data is on vertical edges, increment 3rd dimension
          IF ( doSlots ) THEN
             dims(3) = nSlots
          ELSE
             IF ( doEdges ) dims(3) = dims(3) + 1
          ENDIF

       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( dims(1), dims(2), dims(3) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = .FALSE.

    ENDIF

  END SUBROUTINE Init_and_Register_Log_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_Int_2D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for integer, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_Int_2D( Input_Opt, State_Met, State_Grid,      &
                                       Ptr2Data,  metId,     RC,              &
                                       noRegister                            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name

    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doRegister

    ! Arrays
    INTEGER            :: dims(2)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_Int_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       dims(1) = State_Grid%NX
       dims(2) = State_Grid%NY

       ALLOCATE( Ptr2Data( dims(1), dims(2) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0
    ENDIF

    !========================================================================
    ! Register the field
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_Int_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_Int_3D
!
! !DESCRIPTION: Allocates the data array for a State_Met field,
!  and also adds the field into the State_Met registry.
!  This particular routine is for integer, 3-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_Int_3D( Input_Opt,  State_Met, State_Grid,    &
                                       Ptr2Data,   metId,     RC,            &
                                       noRegister, onEdges,   zxyOrder,      &
                                       nSlots                               )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: metId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
    LOGICAL,          OPTIONAL    :: onEdges             ! Data on vert edges?
    LOGICAL,          OPTIONAL    :: zxyOrder            ! Data array (Z,X,Y)?
    INTEGER,          OPTIONAL    :: nSlots              ! # slots for Z dim
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure!
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: doEdges, doRegister, doSlots, doZxy

    ! Arrays
    INTEGER            :: dims(3)

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_Int_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Met%' // TRIM( metId )

    IF ( PRESENT( onEdges ) ) THEN
       doEdges = onEdges
    ELSE
       doEdges = .FALSE.
    ENDIF

    doSlots = PRESENT( nSlots )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    doSlots = PRESENT( nSlots )

    IF ( PRESENT( zxyOrder ) ) THEN
       doZxy = zxyOrder
    ELSE
       doZxy = .FALSE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array ID and dimensions
       IF ( doZxy ) THEN

          ! ZXY order
          dims(1) = State_Grid%NZ
          dims(2) = State_Grid%NX
          dims(3) = State_Grid%NY

          ! If we have specified nSlots, use that for 1st dimension
          ! Otherwise, if data is on vertical edges, increment 1st dimension
          IF ( doSlots ) THEN
             dims(1) = nSlots
          ELSE
             IF ( doEdges ) dims(1) = dims(1) + 1
          ENDIF

       ELSE

          ! XYZ order
          dims(1) = State_Grid%NX
          dims(2) = State_Grid%NY
          dims(3) = State_Grid%NZ

          ! If we have specified nSlots, use that for 3rd dimension
          ! Otherwise, if data is on vertical edges, increment 3rd dimension
          IF ( doSlots ) THEN
             dims(3) = nSlots
          ELSE
             IF ( doEdges ) dims(3) = dims(3) + 1
          ENDIF

       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( dims(1), dims(2), dims(3) ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_MetField( Input_Opt, metId, Ptr2Data, State_Met, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_Int_3D
!
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_R4_2D
!
! !DESCRIPTION: Registers a 2-D State\_Met field (4-byte real).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_R4_2D( Input_Opt, metadataID, Ptr2Data,       &
                                      State_Met, RC                         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    REAL(f4),         POINTER     :: Ptr2Data(:,:)       ! Pointer to array
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,  type,       vloc
    LOGICAL            :: found

    ! Strings
    CHARACTER(LEN=255) :: desc,  ErrMsg_reg, units, ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    !========================================================================
    ! Register_MetField_R4_2D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    ErrMsg     = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'
    ThisLoc    = &
      ' -> at Register_MetField_R4_2D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Check dimensions
    !========================================================================
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Register the entire 2-D field
    !========================================================================
    CALL Registry_AddField(                                                  &
         Input_Opt   = Input_Opt,                                            &
         Registry    = State_Met%Registry,                                   &
         State       = State_Met%State,                                      &
         Variable    = TRIM( MetadataID ),                                   &
         Units       = TRIM( Units      ),                                   &
         Description = TRIM( Desc       ),                                   &
         Data2d_4    = Ptr2Data,                                             &
         RC          = RC                                                   )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
          '; Abnormal exit from routine "Registry_AddField"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_R4_3D
!
! !DESCRIPTION: Registers a 3-D State\_Met field (4-byte real).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_R4_3D( Input_Opt, metadataID, Ptr2Data,       &
                                      State_Met, RC                         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,     type,     vloc,       N
    LOGICAL            :: found,    onEdges

    ! Strings
    CHARACTER(LEN=2  ) :: numStr
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: desc,     units,    errMsg_reg
    CHARACTER(LEN=255) :: thisDesc, thisName, ThisLoc,    perQnt

    !========================================================================
    ! Register_MetField_R4_3D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    thisName   = ''
    thisDesc   = ''
    ErrMsg     = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'
    ThisLoc    = &
         ' -> at Register_MetField_R4_3D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perQnt     = perQnt,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == vLocationEdge )

    !========================================================================
    ! If there is an Olson landtype dimension,
    ! then register each land type as a 2-D field
    !========================================================================
    IF ( TRIM( perQnt ) == 'OLSON' ) THEN

       ! Check dimensions
       IF ( rank /= 2 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Loop over # of Olson types
       DO N = 1, NSURFTYPE

          ! Append the Olson land type index (zero based)
          ! to the name & description from the metadata
          WRITE( numStr, '(I2.2)') N - 1
          thisName = TRIM( metaDataId )        // numStr
          thisDesc = TRIM( desc       ) // " " // numStr

          ! Register each 2-D field per Olson landtype separately
          CALL Registry_AddField(                                            &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Met%Registry,                            &
               State        = State_Met%State,                               &
               Variable     = TRIM( thisName   ),                            &
               Description  = TRIM( thisDesc   ),                            &
               Units        = TRIM( Units      ),                            &
               Data2d_4     = Ptr2Data(:,:,N),                               &
               RC           = RC                                            )
       ENDDO

    !========================================================================
    ! Otherwise, register as a single 3-D field
    !========================================================================
    ELSE

       ! Check dimensions
       IF ( rank /= 3 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register the entire field
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Met%Registry,                               &
            State        = State_Met%State,                                  &
            Variable     = TRIM( MetadataID ),                               &
            Description  = TRIM( Desc       ),                               &
            Units        = TRIM( Units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data3d_4     = Ptr2Data,                                         &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
             '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_MetField_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_R8_2D
!
! !DESCRIPTION: Registers a 2-D State\_Met field (8-byte real).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_R8_2D( Input_Opt, metadataID, Ptr2Data,      &
                                      State_Met, RC                        )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    REAL(f8),         POINTER     :: Ptr2Data(:,:)       ! Pointer to array
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,  type,       vloc
    LOGICAL            :: found

    ! Strings
    CHARACTER(LEN=255) :: desc,  ErrMsg_reg, units, ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    !========================================================================
    ! Register_MetField_R8_2D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    ErrMsg     = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'
    ThisLoc    = &
      ' -> at Register_MetField_R8_2D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Check dimensions
    !========================================================================
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Register the 2-D field
    !========================================================================
    CALL Registry_AddField(                                                  &
         Input_Opt   = Input_Opt,                                            &
         Registry    = State_Met%Registry,                                   &
         State       = State_Met%State,                                      &
         Variable    = TRIM( MetadataID ),                                   &
         Description = TRIM( Desc       ),                                   &
         Units       = TRIM( Units      ),                                   &
         Data2d_8    = Ptr2Data,                                             &
         RC          = RC                                                   )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
            '; Abnormal exit from routine "Registry_AddField"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_R8_3D
!
! !DESCRIPTION: Registers a 3-D State\_Met field (8-byte real).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_R8_3D( Input_Opt, metadataID, Ptr2Data,        &
                                      State_Met, RC                          )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    REAL(f8),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,     type,     vloc,       N
    LOGICAL            :: found,    onEdges

    ! Strings
    CHARACTER(LEN=2  ) :: numStr
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: desc,     units,    errMsg_reg
    CHARACTER(LEN=255) :: thisDesc, thisName, ThisLoc,    perQnt

    !========================================================================
    ! Register_MetField_R8_3D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    thisDesc   = ''
    thisName   = ''
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Met%'
    thisLoc    = &
       ' -> at Register_MetField_R8_3D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perQnt     = perQnt,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == vLocationEdge )

    !========================================================================
    ! If there is an Olson landtype dimension,
    ! then register each land type as a 2-D field
    !========================================================================
    IF ( TRIM( perQnt ) == 'OLSON' ) THEN

       ! Check dimensions
       IF ( rank /= 2 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Loop over # of Olson types
       DO N = 1, NSURFTYPE

          ! Append the Olson land type index (zero-based)
          ! to the name & description from the metadata
          WRITE( numStr, '(I2.2)') N - 1
          thisName = TRIM( metaDataId )        // numStr
          thisDesc = TRIM( desc       ) // " " // numStr

          ! Register each 2-D field per Olson landtype separately
          CALL Registry_AddField(                                            &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Met%Registry,                            &
               State        = State_Met%State,                               &
               Variable     = TRIM( thisName   ),                            &
               Description  = TRIM( thisDesc   ),                            &
               Units        = TRIM( units      ),                            &
               Data2d_8     = Ptr2Data(:,:,N),                               &
               RC           = RC                                            )

       ENDDO

    !========================================================================
    ! Otherwise, register as a single 3-D field
    !========================================================================
    ELSE

       ! Check dimensions
       IF ( rank /= 3 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register the entire 3-D field
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Met%Registry,                               &
            State        = State_Met%State,                                  &
            Variable     = TRIM( MetadataID ),                               &
            Description  = TRIM( Desc       ),                               &
            Units        = TRIM( Units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data3d_8     = Ptr2Data,                                         &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
             '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_MetField_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_Int_2D
!
! !DESCRIPTION: Registers a 2-D State\_Met field (4-byte integer).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_2D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    INTEGER,          POINTER     :: Ptr2Data(:,:)       ! Pointer to array
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,  type,       vloc
    LOGICAL            :: found, doRegister

    ! Strings
    CHARACTER(LEN=255) :: desc,  errMsg_reg, units, thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Register_MetField_Int_2D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Met%'
    thisLoc    = &
      ' -> at Register_MetField_Int_2D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Check dimensions
    !========================================================================
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Register the 2-D field
    !========================================================================
    CALL Registry_AddField(                                                  &
         Input_Opt   = Input_Opt,                                            &
         Registry    = State_Met%Registry,                                   &
         State       = State_Met%State,                                      &
         Variable    = TRIM( MetadataID ),                                   &
         Description = TRIM( Desc       ),                                   &
         Units       = TRIM( Units      ),                                   &
         Data2d_I    = Ptr2Data,                                             &
         RC          = RC                                                   )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_reg ) // TRIM( metadataID ) //                  &
          '; Abnormal exit from routine "Registry_AddField"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_Int_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_Int_3D
!
! !DESCRIPTION: Registers a 3-D State\_Met field (4-byte integer).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_3D( Input_Opt, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID          ! Field name
    INTEGER,          POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
    TYPE(MetState),   INTENT(IN)  :: State_Met           ! Meteorology State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: rank,     type,     vloc,       N
    LOGICAL            :: found,    onEdges

    ! Strings
    CHARACTER(LEN=2  ) :: numStr
    CHARACTER(LEN=512) :: errMsg
    CHARACTER(LEN=255) :: desc,     units,    errMsg_reg
    CHARACTER(LEN=255) :: thisDesc, thisName, thisLoc,    perQnt

    !========================================================================
    ! Register_MetField_Int_3D begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    thisDesc   = ''
    thisName   = ''
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Met%'
    thisLoc    = &
         ' -> at Register_MetField_Int_3D (in Headers/state_met_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Met(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perQnt     = perQnt,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! If there is an Olson landtype dimension,
    ! then register each land type as a 2-D field
    !========================================================================
    IF ( TRIM( perQnt ) == 'OLSON' ) THEN

       ! Check dimensions
       IF ( rank /= 2 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Loop over # of Olson types
       DO N = 1, NSURFTYPE

          ! Append the Olson land type index to the name & description
          WRITE( numStr, '(I2.2)') N - 1
          thisName = TRIM( metaDataId )        // numStr
          thisDesc = TRIM( desc       ) // " " // numStr

          ! Register each 2-D array per Olson landtype individually
          CALL Registry_AddField(                                            &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Met%Registry,                            &
               State        = State_Met%State,                               &
               Variable     = TRIM( thisName   ),                            &
               Description  = TRIM( thisDesc   ),                            &
               Units        = TRIM( Units      ),                            &
               Data2d_I     = Ptr2Data(:,:,N),                               &
               RC           = RC                                            )

       ENDDO

    !========================================================================
    ! Otherwise, register as a single 3-D field
    !========================================================================
    ELSE

       ! Check dimensions
       IF ( rank /= 3 ) THEN
          ErrMsg = &
             'Data and metadata rank do not match for ' // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register the entire 3-D field
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Met%Registry,                               &
            State        = State_Met%State,                                  &
            Variable     = TRIM( MetadataID ),                               &
            Description  = TRIM( Desc       ),                               &
            Units        = TRIM( Units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data3d_I     = Ptr2Data,                                         &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_reg ) // TRIM( metadataID ) //               &
             '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( errMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_MetField_Int_3D
!EOC
END MODULE State_Met_Mod
