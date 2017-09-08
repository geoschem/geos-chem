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
  USE Precision_Mod
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Met
  PUBLIC :: Cleanup_State_Met
  PUBLIC :: Lookup_State_Met
  PUBLIC :: Print_State_Met
  PUBLIC :: Get_MetField_Metadata
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Register_MetField
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
     REAL(fp), POINTER :: ALBD      (:,:  ) ! Visible surface albedo [1]
     REAL(fp), POINTER :: CLDFRC    (:,:  ) ! Column cloud fraction [1]
     INTEGER,  POINTER :: CLDTOPS   (:,:  ) ! Max cloud top height [levels]
     REAL(fp), POINTER :: EFLUX     (:,:  ) ! Latent heat flux [W/m2]
     REAL(fp), POINTER :: EVAP      (:,:  ) ! Surface evap [kg/m2/s]
     REAL(fp), POINTER :: FRCLND    (:,:  ) ! Olson land fraction [1]
     REAL(fp), POINTER :: FRLAKE    (:,:  ) ! Fraction of lake [1]
     REAL(fp), POINTER :: FRLAND    (:,:  ) ! Fraction of land [1]
     REAL(fp), POINTER :: FRLANDIC  (:,:  ) ! Fraction of land ice [1]
     REAL(fp), POINTER :: FROCEAN   (:,:  ) ! Fraction of ocean [1]
     REAL(fp), POINTER :: FRSEAICE  (:,:  ) ! Sfc sea ice fraction
     REAL(fp), POINTER :: FRSNO     (:,:  ) ! Sfc snow fraction
     REAL(fp), POINTER :: GRN       (:,:  ) ! Greenness fraction
     REAL(fp), POINTER :: GWETROOT  (:,:  ) ! Root soil wetness [1]
     REAL(fp), POINTER :: GWETTOP   (:,:  ) ! Top soil moisture [1]
     REAL(fp), POINTER :: HFLUX     (:,:  ) ! Sensible heat flux [W/m2]
     REAL(fp), POINTER :: LAI       (:,:  ) ! Leaf area index [m2/m2] (online)
     REAL(fp), POINTER :: ITY       (:,:  ) ! Land surface type index
     REAL(fp), POINTER :: LWI       (:,:  ) ! Land/water indices [1]
     REAL(fp), POINTER :: LWI_GISS  (:,:  ) ! Land fraction [1]
     REAL(fp), POINTER :: MOLENGTH  (:,:  ) ! Monin-Obhukov length [m]
     REAL(fp), POINTER :: OICE      (:,:  ) ! Fraction of ocean ice [1]
     REAL(fp), POINTER :: PARDR     (:,:  ) ! Direct  photsyn active rad [W/m2]
     REAL(fp), POINTER :: PARDF     (:,:  ) ! Diffuse photsyn active rad [W/m2]
     REAL(fp), POINTER :: PBLH      (:,:  ) ! PBL height [m]
     INTEGER,  POINTER :: PBL_TOP_L (:,:  ) ! PBL top layer [1]
     REAL(fp), POINTER :: PHIS      (:,:  ) ! Sfc geopotential height [m2/s2]
     REAL(fp), POINTER :: PRECANV   (:,:  ) ! Anvil previp @ ground [kg/m2/s]
     REAL(fp), POINTER :: PRECCON   (:,:  ) ! Conv  precip @ ground [kg/m2/s]
     REAL(fp), POINTER :: PRECTOT   (:,:  ) ! Total precip @ ground [kg/m2/s]
     REAL(fp), POINTER :: PRECLSC   (:,:  ) ! LS precip @ ground [kg/m2/s]
     REAL(fp), POINTER :: PRECSNO   (:,:  ) ! Snow precip [kg/m2/s]
     REAL(fp), POINTER :: PS1_WET   (:,:  ) ! Wet sfc press at dt start[hPa]
     REAL(fp), POINTER :: PS2_WET   (:,:  ) ! Wet sfc press at dt end [hPa]
     REAL(fp), POINTER :: PSC2_WET  (:,:  ) ! Wet interpolated sfc press [hPa]
     REAL(fp), POINTER :: PS1_DRY   (:,:  ) ! Dry sfc press at dt start[hPa]
     REAL(fp), POINTER :: PS2_DRY   (:,:  ) ! Dry sfc press at dt end [hPa]
     REAL(fp), POINTER :: PSC2_DRY  (:,:  ) ! Dry interpolated sfc press [hPa]
     REAL(fp), POINTER :: RADLWG    (:,:  ) ! Net LW radiation @ ground [W/m2]
     REAL(fp), POINTER :: RADSWG    (:,:  ) ! Solar radiation @ ground [W/m2]
     REAL(fp), POINTER :: SEAICE00  (:,:  ) ! Sea ice coverage 00-10%
     REAL(fp), POINTER :: SEAICE10  (:,:  ) ! Sea ice coverage 10-20%
     REAL(fp), POINTER :: SEAICE20  (:,:  ) !  Sea ice coverage 20-30%
     REAL(fp), POINTER :: SEAICE30  (:,:  ) ! Sea ice coverage 30-40%
     REAL(fp), POINTER :: SEAICE40  (:,:  ) ! Sea ice coverage 40-50%
     REAL(fp), POINTER :: SEAICE50  (:,:  ) ! Sea ice coverage 50-60%
     REAL(fp), POINTER :: SEAICE60  (:,:  ) ! Sea ice coverage 60-70%
     REAL(fp), POINTER :: SEAICE70  (:,:  ) ! Sea ice coverage 70-80%
     REAL(fp), POINTER :: SEAICE80  (:,:  ) ! Sea ice coverage 80-90%
     REAL(fp), POINTER :: SEAICE90  (:,:  ) ! Sea ice coverage 90-100%
     REAL(fp), POINTER :: SLP       (:,:  ) ! Sea level pressure [hPa]
     REAL(fp), POINTER :: SNICE     (:,:  ) ! Fraction of snow/ice [1]
     REAL(fp), POINTER :: SNODP     (:,:  ) ! Snow depth [m]
     REAL(fp), POINTER :: SNOMAS    (:,:  ) ! Snow mass [kg/m2]
     REAL(fp), POINTER :: SNOW      (:,:  ) ! Snow depth (H2O equiv) [mm H2O]
     REAL(fp), POINTER :: SST       (:,:  ) ! Sea surface temperature [K]
     REAL(fp), POINTER :: SUNCOS    (:,:  ) ! COS(SZA), current time
     REAL(fp), POINTER :: SUNCOSmid (:,:  ) ! COS(SZA), midpt of chem timestep
!     REAL(fp), POINTER :: SUNCOSmid5(:,:  ) ! COS(SZA), midpt of chem timestep
!                                            !  5 hrs ago (for PARANOX)
     REAL(fp), POINTER :: SWGDN     (:,:  ) ! Incident radiation @ grnd [W/m2]
     REAL(fp), POINTER :: TO3       (:,:  ) ! Total overhead O3 column [DU]
     REAL(fp), POINTER :: TO31      (:,:  ) ! Total O3 at timestep start [DU]
     REAL(fp), POINTER :: TO32      (:,:  ) ! Total O3 at timestep end [DU]
     REAL(fp), POINTER :: TROPP     (:,:  ) ! Tropopause pressure [hPa]
     REAL(fp), POINTER :: TROPP1    (:,:  ) ! Trop P at timestep start [hPa]
     REAL(fp), POINTER :: TROPP2    (:,:  ) ! Trop P at timestep end [hPa]
     REAL(fp), POINTER :: TS        (:,:  ) ! Surface temperature [K]
     REAL(fp), POINTER :: TSKIN     (:,:  ) ! Surface skin temperature [K]
     REAL(fp), POINTER :: TTO3      (:,:  ) ! Tropospheric ozone column [DU]
     REAL(fp), POINTER :: U10M      (:,:  ) ! E/W wind speed @ 10m height [m/s]
     REAL(fp), POINTER :: USTAR     (:,:  ) ! Friction velocity [m/s]
     REAL(fp), POINTER :: UVALBEDO  (:,:  ) ! UV surface albedo [1]
     REAL(fp), POINTER :: V10M      (:,:  ) ! N/S wind speed @ 10m height [m/s]
     REAL(fp), POINTER :: Z0        (:,:  ) ! Surface roughness height [m]
     REAL(fp), POINTER :: CNV_FRC   (:,:  ) ! Convective fraction [1] 
            
     !----------------------------------------------------------------------
     ! 3-D Fields                  
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: AREA_M2   (:,:,:) ! Grid box surface area [cm2]
     REAL(fp), POINTER :: CLDF      (:,:,:) ! 3-D cloud fraction [1]
     REAL(fp), POINTER :: CMFMC     (:,:,:) ! Cloud mass flux [kg/m2/s]
     REAL(fp), POINTER :: DETRAINE  (:,:,:) ! Detrainment (entrain plume)[Pa/s]
     REAL(fp), POINTER :: DETRAINN  (:,:,:) ! Detrainment (non-entr plume)[Pa/s]
     REAL(fp), POINTER :: DNDE      (:,:,:) ! Downdraft (entr plume) [Pa/s]
     REAL(fp), POINTER :: DNDN      (:,:,:) ! Downdraft (non-entr plume) [Pa/s]
     REAL(fp), POINTER :: DQRCU     (:,:,:) ! Conv precip prod rate [kg/kg/s]
                                            ! (assume per dry air)
     REAL(fp), POINTER :: DQRLSAN   (:,:,:) ! LS precip prod rate [kg/kg/s]
                                            ! (assume per dry air)
     REAL(fp), POINTER :: DQIDTMST  (:,:,:) ! Ice tendency, mst proc [kg/kg/s]
     REAL(fp), POINTER :: DQLDTMST  (:,:,:) ! H2O tendency, mst proc [kg/kg/s]
     REAL(fp), POINTER :: DQVDTMST  (:,:,:) ! Vapor tendency, mst proc [kg/kg/s]
     REAL(fp), POINTER :: DTRAIN    (:,:,:) ! Detrainment flux [kg/m2/s]
     REAL(fp), POINTER :: ENTRAIN   (:,:,:) ! GCAP entrainment [Pa/s]
     REAL(fp), POINTER :: HKBETA    (:,:,:) ! Hack overshoot parameter [1]
     REAL(fp), POINTER :: HKETA     (:,:,:) ! Hack conv mass flux [kg/m2/s]
     REAL(fp), POINTER :: MOISTQ    (:,:,:) ! Tendency in sp. humidity 
                                            ! [kg/kg tot air/s]
     REAL(fp), POINTER :: OMEGA     (:,:,:) ! Updraft velocity [Pa/s]
     REAL(fp), POINTER :: OPTD      (:,:,:) ! Visible optical depth [1]
     REAL(fp), POINTER :: PEDGE     (:,:,:) ! Wet air press @ level edges [hPa]
     REAL(fp), POINTER :: PFICU     (:,:,:) ! Dwn flux ice prec:conv [kg/m2/s]
     REAL(fp), POINTER :: PFILSAN   (:,:,:) ! Dwn flux ice prec:LS+anv [kg/m2/s]
     REAL(fp), POINTER :: PFLCU     (:,:,:) ! Dwn flux liq prec:conv [kg/m2/s]
     REAL(fp), POINTER :: PFLLSAN   (:,:,:) ! Dwn flux ice prec:LS+anv [kg/m2/s]
     REAL(fp), POINTER :: PV        (:,:,:) ! Potential vort [kg*m2/kg/s]
     REAL(fp), POINTER :: QI        (:,:,:) ! Ice mixing ratio [kg/kg dry air]
     REAL(fp), POINTER :: QL        (:,:,:) ! Water mixing ratio [kg/kg dry air]
     REAL(fp), POINTER :: REEVAPCN  (:,:,:) ! Evap of precip conv [kg/kg/s]
                                            ! (assume per dry air)
     REAL(fp), POINTER :: REEVAPLS  (:,:,:) ! Evap of precip LS+anvil [kg/kg/s]
                                            ! (assume per dry air)
     REAL(fp), POINTER :: RH        (:,:,:) ! Relative humidity [%]
     REAL(fp), POINTER :: RH1       (:,:,:) ! RH at timestep start [%]
     REAL(fp), POINTER :: RH2       (:,:,:) ! RH at timestep end [%]
     REAL(fp), POINTER :: SPHU      (:,:,:) ! Spcfc humidity [g H2O/kg tot air]
     REAL(fp), POINTER :: SPHU1     (:,:,:) ! Spec hum at timestep start [g/kg]
     REAL(fp), POINTER :: SPHU2     (:,:,:) ! Spec hum at timestep end [g/kg] 
     REAL(fp), POINTER :: T         (:,:,:) ! Temperature [K]
     REAL(fp), POINTER :: TAUCLI    (:,:,:) ! Opt depth of ice clouds [1]
     REAL(fp), POINTER :: TAUCLW    (:,:,:) ! Opt depth of H2O clouds [1]
     REAL(fp), POINTER :: TMPU1     (:,:,:) ! Temperature at timestep start [K]
     REAL(fp), POINTER :: TMPU2     (:,:,:) ! Temperature at timestep end [K]
     REAL(fp), POINTER :: U         (:,:,:) ! E/W component of wind [m s-1]
     REAL(fp), POINTER :: UPDE      (:,:,:) ! Updraft (entraining plume) [Pa/s]
     REAL(fp), POINTER :: UPDN      (:,:,:) ! Updraft (non-entr'n plume) [Pa/s]
     REAL(fp), POINTER :: UPDVVEL   (:,:,:) ! Updraft vertical velocity [hPa/s]
     REAL(fp), POINTER :: V         (:,:,:) ! N/S component of wind [m s-1]
     REAL(fp), POINTER :: ZMEU      (:,:,:) ! Z/M updraft entrainment [Pa/s]
     REAL(fp), POINTER :: ZMMD      (:,:,:) ! Z/M downdraft mass flux [Pa/s]
     REAL(fp), POINTER :: ZMMU      (:,:,:) ! Z/M updraft   mass flux [Pa/s]

     !----------------------------------------------------------------------
     ! Air quantities assigned in AIRQNT
     !----------------------------------------------------------------------
     ! Note on pressures: PMID is calculated from PEDGE, 
     ! and dry air pressures assume constant RH and T across grid box
     REAL(fp), POINTER :: PEDGE_DRY (:,:,:) ! Dry air partial pressure [hPa] 
                                            ! @ level edges [hPa]
     REAL(fp), POINTER :: PMID      (:,:,:) ! Average wet air pressure [hPa]
                                            ! defined as arithmetic
                                            ! average of edge pressures
     REAL(fp), POINTER :: PMID_DRY  (:,:,:) ! Dry air partial pressure [hPa]
                                            ! defined as arithmetic average
                                            ! of edge pressures 
     REAL(fp), POINTER :: TV        (:,:,:) ! Virtual temperature [K]
     REAL(fp), POINTER :: MAIRDEN   (:,:,:) ! Moist air density [kg/m3]
     REAL(fp), POINTER :: AIRDEN    (:,:,:) ! Dry air density [kg/m3]
     REAL(fp), POINTER :: AIRNUMDEN (:,:,:) ! Dry air density [molec/cm3]
     REAL(fp), POINTER :: AVGW      (:,:,:) ! Water vapor volume mixing ratio
                                            ! [vol H2O / vol dry air]
     REAL(fp), POINTER :: BXHEIGHT  (:,:,:) ! Grid box height [m] (dry air)
     REAL(fp), POINTER :: DELP      (:,:,:) ! Delta-P (wet) across box [hPa]
     REAL(fp), POINTER :: DELP_DRY  (:,:,:) ! Delta-P (dry) across box [hPa]
     REAL(fp), POINTER :: AD        (:,:,:) ! Dry air mass [kg] in grid box
     REAL(fp), POINTER :: AIRVOL    (:,:,:) ! Grid box volume [m3] (dry air)
     REAL(fp), POINTER :: DELP_PREV (:,:,:) ! Previous State_Met%DELP
     REAL(fp), POINTER :: DP_DRY_PREV (:,:,:) ! Previous State_Met%DELP_DRY
     REAL(fp), POINTER :: SPHU_PREV (:,:,:) ! Previous State_Met%SPHU

     !----------------------------------------------------------------------
     ! Offline land type, leaf area index, and chlorophyll fields
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: IREG      (:,:  ) ! # of landtypes in grid box (I,J) 
     INTEGER,  POINTER :: ILAND     (:,:,:) ! Land type at (I,J); 1..IREG(I,J)
     INTEGER,  POINTER :: IUSE      (:,:,:) ! Fraction (per mil) of grid box
                                            ! (I,J) occupied by each land type
     REAL(fp), POINTER :: MODISLAI  (:,:  ) ! Daily LAI computed from monthly
                                            ! offline MODIS values [m2/m2]
     REAL(fp), POINTER :: MODISCHLR (:,:  ) ! Daily chlorophyll-a computed from
                                            ! offline MODIS monthly values
     REAL(fp), POINTER :: XLAI      (:,:,:) ! MODIS LAI per land type, this mo
     REAL(fp), POINTER :: XLAI2     (:,:,:) ! MODIS LAI per land type, next mo
     REAL(fp), POINTER :: XCHLR     (:,:,:) ! MODIS CHLR per land type, this mo
     REAL(fp), POINTER :: XCHLR2    (:,:,:) ! MODIS CHLR per land type, next mo
     REAL(fp), POINTER :: LandTypeFrac(:,:,:) ! Olson frac per type (I,J,type)
     REAL(fp), POINTER :: XLAI_NATIVE(:,:,:)  ! avg LAI per type (I,J,type)
     REAL(fp), POINTER :: XCHLR_NATIVE(:,:,:) ! avg CHLR per type (I,J,type)

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Met
     !----------------------------------------------------------------------
     CHARACTER(LEN=3)             :: State     = 'MET'    ! Name of this state
     TYPE(MetaRegItem), POINTER   :: Registry  => NULL()  ! Registry object  

  END TYPE MetState
!
! !REMARKS:
!  In MERRA2, PS and SLP are kept in Pa (not converted to hPa).
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, split off from gc_type_mod.F90
!  23 Oct 2012 - R. Yantosca - Added QI, QL met fields to the derived type
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  12 Dec 2012 - R. Yantosca - Add IREG, ILAND, IUSE fields for dry deposition
!  13 Dec 2012 - R. Yantosca - Add XLAI, XLAI2 fields for dry deposition
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  15 Nov 2013 - R. Yantosca - Now denote that RH fields have units of [%]
!  10 Oct 2014 - C. Keller   - Added ITY (needed for GIGC). For now, this value
!                              is just initialized to 1.0 and not modified any
!                              more.
!  05 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp)
!  12 Feb 2015 - C. Keller   - Added UPDVVEL (for use in wet scavenging).
!  24 Feb 2015 - E. Lundgren - Add PEDGE_DRY, PMID_DRY, and MAIRDEN
!  03 Mar 2015 - E. Lundgren - Add TV (virtual temperature)
!  16 Apr 2015 - E. Lundgren - Add mean pressures PMEAN and PMEAN_DRY. Clarify
!                              definition of PMID as arithmetic average P. 
!                              Add MOISTMW to use TCVV with moist mixing ratio. 
!  25 May 2015 - C. Keller   - Removed SUNCOSmid5 (now calculated by HEMCO).
!  08 Jul 2015 - E. Lundgren - Add XCHLR and XCHLR2 for organic marine aerosols
!  11 Aug 2015 - R. Yantosca - Extend #ifdefs for MERRA2 met fields
!  22 Sep 2015 - E. Lundgren - Add SWGDN for incident radiation at ground
!  28 Oct 2015 - E. Lundgren - Add previous delta-P and specific humidity for
!                              tracer mass conservation in mixing ratio update
!  04 Mar 2016 - C. Keller   - Add CNV_FRC for convective fraction. Currently
!                              not a standard GEOS-FP output, only used in 
!                              online model (ESMF). 
!  21 Dec 2015 - M. Sulprizio- Add AIRNUMDEN, which is the same as AIRDEN but
!                              has units molec/cm3 for the chemistry routines.
!  17 Mar 2016 - M. Sulprizio- Remove OPTDEP. Instead, we now solely use OPTD.
!  03 May 2016 - E. Lundgren - Add PSC2_DRY, PS1_DRY, and PS2_DRY
!  06 Jul 2016 - E. Lundgren - Rename PS1, PS2, and PSC1: add '_WET' suffix
!  06 Jul 2016 - E. Lundgren - Add DELP_DRY and DP_DRY_PREV
!  19 Jul 2016 - E. Lundgren - Remove PMEAN, PMEAN_DRY, MOISTMW, and ADMOIST  
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_state_chm_mod.F90 to
!                              state_chm_mod.F90. The "gigc" nomenclature is
!                              no longer used.
!  18 Oct 2016 - E. Lundgren - Remove XLAI2, CHLR2; add MODISLAI, MODISCHLR to
!                              replace modis_lai_mod-level GC_LAI and GC_CHLR
!  19 Oct 2016 - E. Lundgren - Use NSURFTYPE as the # of land types
!  03 Feb 2017 - M. Sulprizio- Add OMEGA for use in sulfate_mod.F (Q. Chen)
!  26 Jun 2017 - R. Yantosca - Added StateName and Registry to type MetState
!  27 Jun 2017 - R. Yantosca - Add fields of State_Met to the registry
!  07 Sep 2017 - E. Lundgren - Add Register_MetField interface for init
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE REGISTER_METFIELD
     MODULE PROCEDURE REGISTER_METFIELD_Rfp_2D
     MODULE PROCEDURE REGISTER_METFIELD_Rfp_3D
     MODULE PROCEDURE REGISTER_METFIELD_Int_2D
     MODULE PROCEDURE REGISTER_METFIELD_Int_3D
  END INTERFACE REGISTER_METFIELD

CONTAINS
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
  SUBROUTINE Init_State_Met( am_I_Root, IM, JM, LM, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE CMN_SIZE_MOD, ONLY : NSURFTYPE
    USe Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
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
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!  23 Oct 2012 - R. Yantosca - Now allocate QI, QL fields
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  16 Nov 2012 - R. Yantosca - Now zero all fields after allocating
!  27 Nov 2012 - R. Yantosca - Now allocate SUNCOS fields (IM,JM)
!  12 Dec 2012 - R. Yantosca - Now allocate the IREG, ILAND, IUSE fields
!  13 Dec 2012 - R. Yantosca - Now allocate the XLAI, XLAI2 fields
!  07 Mar 2013 - R. Yantosca - Now allocate PF*LSAN, PF*CU fields properly
!                              for GEOS-5.7.x met (they are edged)
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  22 Aug 2014 - R. Yantosca - Allocate PBL_TOP_L field
!  05 Nov 2014 - R. Yantosca - Now use 0.0_fp instead of 0d0
!  06 Nov 2014 - R. Yantosca - Now make all fields (IM,JM,LM) instead of 
!                              (LM,JM,IM), to facilitate use w/in GEOS-5 GCM
!  05 Oct 2016 - R. Yantosca - Swapped order of HKETA and HKBETA allocation
!  28 Nov 2016 - R. Yantosca - Nullify fields that may or may not be allocated
!  01 Jun 2017 - C. Keller   - Initialize UPDVVEL to -999.0 to ensure that 
!                              GET_VUD (wetscav_mod.F) works properly.
!  26 Jun 2017 - R. Yantosca - Now register each variable after it's allocated
!  24 Aug 2017 - R. Yantosca - Now register level-edged variables appropriately
!  07 Sep 2017 - E. Lundgren - Abstract the metadata and adding to registry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: LX
    
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC    =  GC_SUCCESS

    !=======================================================================
    ! The following fields of State_Met may or may not get allocated
    ! depending on the met field being used, or if we are using GEOS-Chem
    ! in the ESMF/HPC configuration.  Make sure to nullify these fields
    ! in order to prevent issues with unintialized fields.  In particular,
    ! the GNU Fortran compiler may cause simulations to die with an error 
    ! when encountering uninitialized fields of State_Met.
    ! 
    ! We do not have to nullify the fields that always get allocated,
    ! since they will be defined for each GEOS-Chem simulation. 
    ! (bmy, 11/28/16) 
    !=======================================================================
    State_Met%CNV_FRC  => NULL()
    State_Met%DETRAINE => NULL()
    State_Met%DETRAINN => NULL()
    State_Met%DNDE     => NULL()
    State_Met%DNDN     => NULL()
    State_Met%ENTRAIN  => NULL()
    State_Met%FRSEAICE => NULL()
    State_Met%FRSNO    => NULL()
    State_Met%HKETA    => NULL()
    State_Met%HKBETA   => NULL()
    State_Met%LWI_GISS => NULL()
    State_Met%MOLENGTH => NULL()
    State_Met%OICE     => NULL()
    State_Met%PFICU    => NULL()
    State_Met%PFILSAN  => NULL()
    State_Met%PFLCU    => NULL()
    State_Met%PFLLSAN  => NULL()
    State_Met%PRECANV  => NULL()
    State_Met%PRECLSC  => NULL()
    State_Met%REEVAPCN => NULL()
    State_Met%REEVAPLS => NULL()
    State_Met%RH1      => NULL()
    State_Met%RH2      => NULL()
    State_Met%SEAICE00 => NULL()
    State_Met%SEAICE10 => NULL()
    State_Met%SEAICE20 => NULL()
    State_Met%SEAICE30 => NULL()
    State_Met%SEAICE40 => NULL()
    State_Met%SEAICE50 => NULL()
    State_Met%SEAICE60 => NULL()
    State_Met%SEAICE70 => NULL()
    State_Met%SEAICE80 => NULL()
    State_Met%SEAICE90 => NULL()
    State_Met%SNICE    => NULL()
    State_Met%SNOW     => NULL()
    State_Met%SNOW     => NULL()
    State_Met%SPHU     => NULL()
    State_Met%SPHU1    => NULL()
    State_Met%SPHU2    => NULL()
    State_Met%TMPU1    => NULL()
    State_Met%TMPU2    => NULL()
    State_Met%TO31     => NULL()
    State_Met%TO32     => NULL()
    State_Met%TROPP1   => NULL()
    State_Met%TROPP2   => NULL()
    State_Met%TTO3     => NULL()
    State_Met%UPDE     => NULL()
    State_Met%UPDN     => NULL()
    State_Met%ZMEU     => NULL()
    State_Met%ZMMD     => NULL()
    State_Met%ZMMU     => NULL()

    !=======================================================================
    ! Allocate 2-D Fields
    !=======================================================================

    !-------------------------
    ! ALBD [1]
    !-------------------------
    ALLOCATE( State_Met%ALBD( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%ALBD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ALBD = 0.0_fp
    CALL Register_MetField( am_I_Root, 'ALBD', State_Met%ALBD, &
                            State_Met, RC )

    !-------------------------
    ! CLDFRC [1]
    !-------------------------
    ALLOCATE( State_Met%CLDFRC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDFRC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CLDFRC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CLDFRC', State_Met%CLDFRC, &
                            State_Met, RC )

    !-------------------------
    ! CLDTOPS [level]
    !-------------------------
    ALLOCATE( State_Met%CLDTOPS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDTOPS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CLDTOPS = 0
    CALL Register_MetField( am_I_Root, 'CLDTOPS', State_Met%CLDTOPS, &
                            State_Met, RC )

    !-------------------------
    ! EFLUX [W m-2]
    !-------------------------
    ALLOCATE( State_Met%EFLUX( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%EFLUX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%EFLUX    = 0.0_fp
    CALL Register_MetField( am_I_Root, 'EFLUX', State_Met%EFLUX, &
                            State_Met, RC )

    !-------------------------
    ! EVAP [W m-2]
    !-------------------------
    ALLOCATE( State_Met%EVAP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%EVAP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%EVAP= 0.0_fp
    CALL Register_MetField( am_I_Root, 'EVAP', State_Met%EVAP, &
                            State_Met, RC )

    !-------------------------
    ! FRCLND [1]
    !-------------------------
    ALLOCATE( State_Met%FRCLND( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRCLND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRCLND = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRCLND', State_Met%FRCLND, &
                            State_Met, RC )

    !-------------------------
    ! FRLAKE [1]
    !-------------------------
    ALLOCATE( State_Met%FRLAKE( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLAKE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLAKE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRLAKE', State_Met%FRLAKE, &
                            State_Met, RC )

    !-------------------------
    ! FRLAND [1]
    !-------------------------
    ALLOCATE( State_Met%FRLAND( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLAND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLAND = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'FRLAND', State_Met%FRLAND, &
                            State_Met, RC )

    !-------------------------
    ! FRLANDIC [1]
    !-------------------------
    ALLOCATE( State_Met%FRLANDIC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLANDIC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLANDIC = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'FRLANDIC', State_Met%FRLANDIC, &
                            State_Met, RC )

    !-------------------------
    ! FROCEAN [1]
    !-------------------------
    ALLOCATE( State_Met%FROCEAN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FROCEAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FROCEAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FROCEAN', State_Met%FROCEAN, &
                            State_Met, RC )

    !-------------------------
    ! GRN [1]
    !-------------------------
    ALLOCATE( State_Met%GRN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%GRN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%GRN = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'GRN', State_Met%GRN, &
                            State_Met, RC )

    !-------------------------
    ! GWETROOT [1]
    !-------------------------
    ALLOCATE( State_Met%GWETROOT( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%GWETROOT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%GWETROOT = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'GWETROOT', State_Met%GWETROOT, &
                            State_Met, RC )

    !-------------------------
    ! GWETTOP [1]
    !-------------------------
    ALLOCATE( State_Met%GWETTOP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%GWETTOP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%GWETTOP = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'GWETTOP', State_Met%GWETTOP, &
                            State_Met, RC )

    !-------------------------
    ! HFLUX [W m-2]
    !-------------------------
    ALLOCATE( State_Met%HFLUX( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%HFLUX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%HFLUX = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'HFLUX', State_Met%HFLUX, &
                            State_Met, RC )

    !-------------------------
    ! LAI [1]
    !-------------------------  
    ALLOCATE( State_Met%LAI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%LAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LAI', State_Met%LAI, &
                            State_Met, RC )

    !-------------------------
    ! ITY [1]
    !-------------------------
    ALLOCATE( State_Met%ITY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%ITY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ITY = 1.0_fp
    CALL Register_MetField( am_I_Root, 'ITY', State_Met%ITY, &
                            State_Met, RC )

    !-------------------------
    ! LWI [1]
    !-------------------------
    ALLOCATE( State_Met%LWI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%LWI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LWI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LWI', State_Met%LWI, &
                            State_Met, RC )

    !-------------------------
    ! PARDR [W m-2]
    !-------------------------
    ALLOCATE( State_Met%PARDR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PARDR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PARDR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PARDR', State_Met%PARDR, &
                            State_Met, RC )

    !-------------------------
    ! PARDF [W m-2]
    !-------------------------
    ALLOCATE( State_Met%PARDF( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PARDF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PARDF= 0.0_fp
    CALL Register_MetField( am_I_Root, 'PARDF', State_Met%PARDF, &
                            State_Met, RC )

    !-------------------------
    ! PBLH [m]
    !-------------------------
    ALLOCATE( State_Met%PBLH( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PBLH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PBLH = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PBLH', State_Met%PBLH, &
                            State_Met, RC )

    !-------------------------
    ! PBL_TOP_L [1]
    !-------------------------
    ALLOCATE( State_Met%PBL_TOP_L( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PBL_TOP_L', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PBL_TOP_L = 0
    CALL Register_MetField( am_I_Root, 'PBL_TOP_L', State_Met%PBL_TOP_L, &
                            State_Met, RC )

    !-------------------------
    ! PHIS [m2 s-2]
    !-------------------------
    ALLOCATE( State_Met%PHIS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PHIS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PHIS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PHIS', State_Met%PHIS, &
                            State_Met, RC )

    !-------------------------
    ! PRECCON [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECCON( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECCON', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECCON = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECCON', State_Met%PRECCON, &
                            State_Met, RC )

    !-------------------------
    ! PRECSNO [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECSNO( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECSNO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECSNO = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECSNO', State_Met%PRECSNO, &
                            State_Met, RC )

    !-------------------------
    ! PRECTOT [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECTOT( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECTOT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECTOT = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECTOT', State_Met%PRECTOT, &
                            State_Met, RC )

    !-------------------------
    ! PS1_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS1_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS1_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS1_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS1_WET', State_Met%PS1_WET, &
                            State_Met, RC )

    !-------------------------
    ! PS2_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS2_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS2_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS2_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS2_WET', State_Met%PS2_WET, &
                            State_Met, RC )

    !-------------------------
    ! PSC2_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PSC2_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PSC2_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PSC2_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PSC2_WET', State_Met%PSC2_WET, &
                            State_Met, RC )

    !-------------------------
    ! PS1_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS1_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS1_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS1_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS1_DRY', State_Met%PS1_DRY, &
                            State_Met, RC )

    !-------------------------
    ! PS2_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS2_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS2_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS2_DRY   = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS2_DRY', State_Met%PS2_DRY, &
                            State_Met, RC )

    !-------------------------
    ! PSC2_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PSC2_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PSC2_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PSC2_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PSC2_DRY', State_Met%PSC2_DRY, &
                            State_Met, RC )

    !-------------------------
    ! RADLWG [W m-2]
    !-------------------------
    ALLOCATE( State_Met%RADLWG( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RADLWG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%RADLWG = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RADLWG', State_Met%RADLWG, &
                            State_Met, RC )

    !-------------------------
    ! RADSWG [W m-2]
    !-------------------------
    ALLOCATE( State_Met%RADSWG( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RADSWG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%RADSWG = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RADSWG', State_Met%RADSWG, &
                            State_Met, RC )

    !-------------------------
    ! SLP [hPa]
    !-------------------------
    ALLOCATE( State_Met%SLP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SLP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SLP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SLP', State_Met%SLP, &
                            State_Met, RC )

    !-------------------------
    ! SNODP [m]
    !-------------------------
    ALLOCATE( State_Met%SNODP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SNODP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNODP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SNODP', State_Met%SNODP, &
                            State_Met, RC )

    !-------------------------
    ! SNOMAS [kg m-2]
    !-------------------------
    ALLOCATE( State_Met%SNOMAS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SNOMAS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNOMAS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SNOMAS', State_Met%SNOMAS, &
                            State_Met, RC )

    !-------------------------
    ! SST [K]
    !-------------------------
    ALLOCATE( State_Met%SST( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SST = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SST', State_Met%SST, &
                            State_Met, RC )

    !-------------------------
    ! SUNCOS [1]
    !-------------------------
    ALLOCATE( State_Met%SUNCOS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SUNCOS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SUNCOS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SUNCOS', State_Met%SUNCOS, &
                            State_Met, RC )

    !-------------------------
    ! SUNCOSmid [1]
    !-------------------------
    ALLOCATE( State_Met%SUNCOSmid( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SUNCOSmid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SUNCOSmid = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SUNCOSmid', State_Met%SUNCOSmid, &
                            State_Met, RC )

    !-------------------------
    ! SWGDN [W m-2]
    !-------------------------
    ALLOCATE( State_Met%SWGDN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SWGDN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SWGDN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SWGDN', State_Met%SWGDN, &
                            State_Met, RC )

    !-------------------------
    ! TO3 [dobsons]
    !-------------------------
    ALLOCATE( State_Met%TO3( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TO3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TO3 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TO3', State_Met%TO3, &
                            State_Met, RC )

    !-------------------------
    ! TROPP [hPa]
    !-------------------------
    ALLOCATE( State_Met%TROPP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TROPP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TROPP', State_Met%TROPP, &
                            State_Met, RC )

    !-------------------------
    ! TS [K]
    !-------------------------
    ALLOCATE( State_Met%TS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TS', State_Met%TS, &
                            State_Met, RC )

    !-------------------------
    ! TSKIN [1]
    !-------------------------
    ALLOCATE( State_Met%TSKIN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TSKIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TSKIN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TSKIN', State_Met%TSKIN, &
                            State_Met, RC )

    !-------------------------
    ! U10M [m s-1]
    !-------------------------
    ALLOCATE( State_Met%U10M( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%U10M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%U10M = 0.0_fp
    CALL Register_MetField( am_I_Root, 'U10M', State_Met%U10M, &
                            State_Met, RC )

    !-------------------------
    ! USTAR [m -s]
    !-------------------------
    ALLOCATE( State_Met%USTAR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%USTAR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%USTAR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'USTAR', State_Met%USTAR, &
                            State_Met, RC )

    !-------------------------
    ! UVALBEDO [1]
    !-------------------------
    ALLOCATE( State_Met%UVALBEDO( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%UVALBEDO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UVALBEDO = 0.0_fp
    CALL Register_MetField( am_I_Root, 'UVALBEDO', State_Met%UVALBEDO, &
                            State_Met, RC )

    !-------------------------
    ! V10M [m s-1]
    !-------------------------
    ALLOCATE( State_Met%V10M( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%V10M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%V10M = 0.0_fp
    CALL Register_MetField( am_I_Root, 'V10M', State_Met%V10M, &
                            State_Met, RC )

    !-------------------------
    ! Z0 [m]
    !-------------------------
    ALLOCATE( State_Met%Z0( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Z0', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%Z0 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'Z0', State_Met%Z0, &
                            State_Met, RC )

    ! Convective fractions are not yet a standard GEOS-FP
    ! field. Only available to online model (ckeller, 3/4/16) 
#if defined( ESMF_ )
    !-------------------------
    ! CNV_FRC [1]
    !-------------------------
    ALLOCATE( State_Met%CNV_FRC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CNV_FRC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CNV_FRC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CNV_FRC', State_Met%CNV_FRC, &
                            State_Met, RC )
#endif

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !
    ! NOTE: Do not add the GCAP fields to the registry, since these
    !       may be deprecated (GCAP2 is probably used instead).
    !=======================================================================
    ALLOCATE( State_Met%LWI_GISS  ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LWI_GISS = 0.0_fp

    ALLOCATE( State_Met%MOLENGTH  ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%MOLENGTH = 0.0_fp

    ALLOCATE( State_Met%OICE      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%OICE     = 0.0_fp

    ALLOCATE( State_Met%SNICE     ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNICE    = 0.0_fp

    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNOW     = 0.0_fp

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0.0_fp

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0.0_fp

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !
    ! NOTE: Do not add the GEOS-4 fields to the registry, since these
    !       are deprecated and are slated to be removed soon.
    !=======================================================================
    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNOW     = 0.0_fp

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0.0_fp

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0.0_fp

#elif defined( GEOS_5 )

    !=======================================================================
    ! GEOS-5 met fields
    !
    ! NOTE: Do not add the GCAP fields to the registry, since these
    !       may be deprecated and are slated to be removed soon.
    !=======================================================================
    ALLOCATE( State_Met%TO31      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TO31     = 0.0_fp

    ALLOCATE( State_Met%TO32      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TO32     = 0.0_fp

    ALLOCATE( State_Met%TTO3      ( IM, JM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TTO3     = 0.0_fp

#elif defined( GEOS_FP ) || defined( MERRA ) || defined( MERRA2 )

    !=======================================================================
    ! MERRA, GEOS-FP, and MERRA2 met fields
    !=======================================================================

    !-------------------------
    ! FRESEAICE [1]
    !-------------------------
    ALLOCATE( State_Met%FRSEAICE( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRSEAICE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRSEAICE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRSEAICE', State_Met%FRSEAICE, &
                            State_Met, RC )

    !-------------------------
    ! FRSNO [1]
    !-------------------------
    ALLOCATE( State_Met%FRSNO( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRSNO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRSNO = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRSNO', State_Met%FRSNO, &
                            State_Met, RC )

    !-------------------------
    ! PRECANV [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECANV( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECANV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECANV = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECANV', State_Met%PRECANV, &
                            State_Met, RC )

    !-------------------------
    ! PRECLSC [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECLSC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECLSC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECLSC  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECLSC', State_Met%PRECLSC, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE00 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE00( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE00', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE00 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE00', State_Met%SEAICE00, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE10 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE10( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE10', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE10 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE10', State_Met%SEAICE10, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE20 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE20( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE20', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE20 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE20', State_Met%SEAICE20, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE30 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE30( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE30', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE30 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE30', State_Met%SEAICE30, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE40 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE40( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE40', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE40 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE40', State_Met%SEAICE40, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE50 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE50( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE50', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE50 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE50', State_Met%SEAICE50, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE60 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE60( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE60', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE60 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE60', State_Met%SEAICE60, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE70 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE70( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE70', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE70 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE70', State_Met%SEAICE70, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE80 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE80( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE80', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE80 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE80', State_Met%SEAICE80, &
                            State_Met, RC )

    !-------------------------
    ! SEAICE90 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE90( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE90', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE90 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE90', State_Met%SEAICE90, &
                            State_Met, RC )

#endif

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================

    !-------------------------
    ! AD [kg]
    !-------------------------
    ALLOCATE( State_Met%AD( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AD = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AD', State_Met%AD, &
                            State_Met, RC )

    !-------------------------
    ! AIRDEN [kg m-3]
    !-------------------------
    ALLOCATE( State_Met%AIRDEN( IM, JM, LM ), STAT=RC )  
    CALL GC_CheckVar( 'State_Met%AIRDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AIRDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRDEN', State_Met%AIRDEN, &
                            State_Met, RC )

    !-------------------------
    ! MAIRDEN [kg m-3]
    !-------------------------
    ALLOCATE( State_Met%MAIRDEN( IM, JM, LM   ), STAT=RC )  
    CALL GC_CheckVar( 'State_Met%MAIRDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%MAIRDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MAIRDEN', State_Met%MAIRDEN, &
                            State_Met, RC )

    !-------------------------
    ! AIRNUMDEN [1]
    !-------------------------        
    ALLOCATE( State_Met%AIRNUMDEN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AIRNUMDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%AIRNUMDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRNUMDEN', State_Met%AIRNUMDEN, &
                            State_Met, RC )

    !-------------------------
    ! AIRVOL [m3]
    !-------------------------
    ALLOCATE( State_Met%AIRVOL( IM, JM, LM  ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AIRVOL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AIRVOL = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRVOL', State_Met%AIRVOL, &
                            State_Met, RC )

    !-------------------------
    ! AREA_M2 [m2]
    !-------------------------
    ALLOCATE( State_Met%AREA_M2( IM, JM, 1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AREA_M2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AREA_M2  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AREA_M2', State_Met%AREA_M2, &
                            State_Met, RC )

    !-------------------------
    ! AVGW [v/v]
    !-------------------------
    ALLOCATE( State_Met%AVGW( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AVGW', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%AVGW = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AVGW', State_Met%AVGW, &
                            State_Met, RC )

    !-------------------------
    ! BXHEIGHT [m]
    !-------------------------
    ALLOCATE( State_Met%BXHEIGHT( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%BXHEIGHT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%BXHEIGHT = 0.0_fp
    CALL Register_MetField( am_I_Root, 'BXHEIGHT', State_Met%BXHEIGHT, &
                            State_Met, RC )

    !-------------------------
    ! CLDF [1]
    !-------------------------                 
    ALLOCATE( State_Met%CLDF( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%CLDF = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CLDF', State_Met%CLDF, &
                            State_Met, RC )

    !-------------------------
    ! CMFMC [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%CMFMC( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CMFMC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%CMFMC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CMFMC', State_Met%CMFMC, &
                            State_Met, RC )

    !-------------------------
    ! DELP [hPa]
    !-------------------------           
    ALLOCATE( State_Met%DELP( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DELP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DELP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DELP', State_Met%DELP, &
                            State_Met, RC )

    !-------------------------
    ! DELP_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%DELP_DRY( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DELP_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DELP_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DELP_DRY', State_Met%DELP_DRY, &
                            State_Met, RC )

    !-------------------------
    ! DELP_PREV [hPa]
    !-------------------------
    ALLOCATE( State_Met%DELP_PREV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DELP_PREV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DELP_PREV= 0.0_fp
    CALL Register_MetField( am_I_Root, 'DELP_PREV', State_Met%DELP_PREV, &
                            State_Met, RC )

    !-------------------------
    ! DP_DRY_PREV [hPa]
    !-------------------------
    ALLOCATE( State_Met%DP_DRY_PREV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DP_DRY_PREV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DP_DRY_PREV= 0.0_fp
    CALL Register_MetField( am_I_Root, 'DP_DRY_PREV', State_Met%DP_DRY_PREV, &
                            State_Met, RC )

    !-------------------------
    ! DQRCU [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQRCU( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQRCU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DQRCU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQRCU', State_Met%DQRCU, &
                            State_Met, RC )

    !-------------------------
    ! DQRLSAN [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQRLSAN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQRLSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DQRLSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQRLSAN', State_Met%DQRLSAN, &
                            State_Met, RC )

    !-------------------------
    ! DQIDTMST [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQIDTMST( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQIDTMST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%DQIDTMST = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQIDTMST', State_Met%DQIDTMST, &
                            State_Met, RC )

    !-------------------------
    ! DQLDTMST [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQLDTMST( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQLDTMST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%DQLDTMST = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQLDTMST', State_Met%DQLDTMST, &
                            State_Met, RC )

    !-------------------------
    ! DQVDTMST [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQVDTMST( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQVDTMST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%DQVDTMST = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQVDTMST', State_Met%DQVDTMST, &
                            State_Met, RC )

    !-------------------------
    ! DTRAIN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%DTRAIN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DTRAIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%DTRAIN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DTRAIN', State_Met%DTRAIN, &
                            State_Met, RC )

    !-------------------------
    ! MOISTQ [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%MOISTQ( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%MOISTQ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%MOISTQ = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MOISTQ', State_Met%MOISTQ, &
                            State_Met, RC )

    !-------------------------
    ! OMEGA [Pa s-1]
    !-------------------------
    ALLOCATE( State_Met%OMEGA( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%OMEGA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%OMEGA = 0.0_fp
    CALL Register_MetField( am_I_Root, 'OMEGA', State_Met%OMEGA, &
                            State_Met, RC )

    !-------------------------
    ! OPTD [1]
    !-------------------------
    ALLOCATE( State_Met%OPTD( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%OPTD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%OPTD = 0.0_fp
    CALL Register_MetField( am_I_Root, 'OPTD', State_Met%OPTD, &
                            State_Met, RC )
                            
    !-------------------------
    ! PEDGE [hPa]
    !-------------------------
    ALLOCATE( State_Met%PEDGE( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PEDGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PEDGE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PEDGE', State_Met%PEDGE, &
                            State_Met, RC )

    !-------------------------
    ! PEDGE_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PEDGE_DRY ( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PEDGE_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PEDGE_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PEDGE_DRY', State_Met%PEDGE_DRY, &
                            State_Met, RC )

    !-------------------------
    ! PMID [1]
    !-------------------------
    ALLOCATE( State_Met%PMID( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PMID', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PMID = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PMID', State_Met%PMID, &
                            State_Met, RC )

    !-------------------------
    ! PMID_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PMID_DRY( IM, JM, LM   ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PMID_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PMID_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PMID_DRY', State_Met%PMID_DRY, &
                            State_Met, RC )

    !-------------------------
    ! PV [kg m2 kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%PV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PV = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PV', State_Met%PV, &
                            State_Met, RC )

    !-------------------------
    ! QI [kg kg-1]
    !-------------------------
    ALLOCATE( State_Met%QI( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%QI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%QI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'QI', State_Met%QI, &
                            State_Met, RC )

    !-------------------------
    ! QL [kg kg-1]
    !-------------------------
    ALLOCATE( State_Met%QL( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%QL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%QL = 0.0_fp
    CALL Register_MetField( am_I_Root, 'QL', State_Met%QL, &
                            State_Met, RC )

    !-------------------------
    ! RH [%]
    !-------------------------
    ALLOCATE( State_Met%RH ( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN                               
    State_Met%RH = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RH', State_Met%RH, &
                            State_Met, RC )

    !-------------------------
    ! SPHU [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%SPHU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU', State_Met%SPHU, &
                            State_Met, RC )

    !-------------------------
    ! SPHU_PREV [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU_PREV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU_PREV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%SPHU_PREV= 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU_PREV', State_Met%SPHU_PREV, &
                            State_Met, RC )

    !-------------------------
    ! T [K]
    !-------------------------                                               
    ALLOCATE( State_Met%T( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%T', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%T = 0.0_fp
    CALL Register_MetField( am_I_Root, 'T', State_Met%T, &
                            State_Met, RC )

    !-------------------------
    ! TV [K]
    !-------------------------
    ALLOCATE( State_Met%TV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%TV = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TV', State_Met%TV, &
                            State_Met, RC )

    !-------------------------
    ! TAUCLI [1]
    !-------------------------
    ALLOCATE( State_Met%TAUCLI( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TAUCLI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%TAUCLI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TAUCLI', State_Met%TAUCLI, &
                            State_Met, RC )

    !-------------------------
    ! TAUCLW [1]
    !-------------------------
    ALLOCATE( State_Met%TAUCLW( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TAUCLW', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TAUCLW = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TAUCLW', State_Met%TAUCLW, &
                            State_Met, RC )

    !-------------------------
    ! U [m s-1]
    !-------------------------
    ALLOCATE( State_Met%U( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%U', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%U = 0.0_fp
    CALL Register_MetField( am_I_Root, 'U', State_Met%U, &
                            State_Met, RC )

    !-------------------------
    ! V [m s-1]
    !-------------------------
    ALLOCATE( State_Met%V( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%V', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%V = 0.0_fp
    CALL Register_MetField( am_I_Root, 'V', State_Met%V, &
                            State_Met, RC )

    !-------------------------
    ! UPDVVEL [hPa s-1]
    !-------------------------
    ALLOCATE( State_Met%UPDVVEL( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%UPDVVEL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UPDVVEL  = -999.0_fp
    CALL Register_MetField( am_I_Root, 'UPDVVEL', State_Met%UPDVVEL, &
                            State_Met, RC )

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !
    ! NOTE: Do not add the GCAP fields to the registry, since these
    !       may be deprecated (GCAP2 is probably used instead).
    !=======================================================================
    ALLOCATE( State_Met%DETRAINE  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DETRAINE = 0.0_fp

    ALLOCATE( State_Met%DETRAINN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DETRAINN = 0.0_fp

    ALLOCATE( State_Met%DNDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DNDE     = 0.0_fp

    ALLOCATE( State_Met%DNDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DNDN     = 0.0_fp

    ALLOCATE( State_Met%ENTRAIN   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ENTRAIN  = 0.0_fp

    ALLOCATE( State_Met%UPDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UPDE     = 0.0_fp

    ALLOCATE( State_Met%UPDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UPDN     = 0.0_fp

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !
    ! NOTE: Do not add the GEOS-4 fields to the registry, since these
    !       are deprecated and are slated to be removed soon.
    !=======================================================================

    ALLOCATE( State_Met%HKETA     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%HKETA    = 0.0_fp

    ALLOCATE( State_Met%HKBETA    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%HKBETA   = 0.0_fp

    ALLOCATE( State_Met%ZMEU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ZMEU     = 0.0_fp

    ALLOCATE( State_Met%ZMMD      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ZMMD     = 0.0_fp

    ALLOCATE( State_Met%ZMMU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ZMMU     = 0.0_fp

#elif defined( GEOS_FP ) || defined( MERRA ) || defined( MERRA2 )

    !=======================================================================
    ! MERRA, GEOS-FP, and MERRA2 met fields
    !=======================================================================

    ! Pick the proper vertical dimension
#if defined( GEOS_FP ) || defined( MERRA2 )
    LX = LM + 1           ! For fields that are on level edges
#else
    LX = LM               ! For fields that are on level centers
#endif
    
    !-------------------------
    ! PFICU [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFICU( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFICU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFICU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFICU', State_Met%PFICU, &
                            State_Met, RC )

    !-------------------------
    ! PFILSAN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFILSAN( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFILSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFILSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFILSAN', State_Met%PFILSAN, &
                            State_Met, RC )

    !-------------------------
    ! PFLCU [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFLCU( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFLCU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFLCU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFLCU', State_Met%PFLCU, &
                            State_Met, RC )

    !-------------------------
    ! PFLLSAN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFLLSAN( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFLLSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFLLSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFLLSAN', State_Met%PFLLSAN, &
                            State_Met, RC )

    !-------------------------
    ! REEVAPCN [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%REEVAPCN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%REEVAPCN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%REEVAPCN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'REEVAPCN', State_Met%REEVAPCN, &
                            State_Met, RC )

    !-------------------------
    ! REEVAPLS [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%REEVAPLS( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%REEVAPLS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%REEVAPLS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'REEVAPLS', State_Met%REEVAPLS, &
                            State_Met, RC )

    !-------------------------
    ! RH1 [%]
    !-------------------------
    ALLOCATE( State_Met%RH1( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RH1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%RH1 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RH1', State_Met%RH1, &
                            State_Met, RC )

    !-------------------------
    ! RH2 [%]
    !-------------------------
    ALLOCATE( State_Met%RH2( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RH2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%RH2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RH2', State_Met%RH2, &
                            State_Met, RC )

    !-------------------------
    ! SPHU1 [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU1( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SPHU1 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU1', State_Met%SPHU1, &
                            State_Met, RC )

    !-------------------------
    ! SPHU2 [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU2( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SPHU2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU2', State_Met%SPHU2, &
                            State_Met, RC )

    !-------------------------
    ! TMPU1 [K]
    !-------------------------
    ALLOCATE( State_Met%TMPU1( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TMPU1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TMPU1 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TMPU1', State_Met%TMPU1, &
                            State_Met, RC )

    !-------------------------
    ! TMPU2 [K]
    !-------------------------
    ALLOCATE( State_Met%TMPU2( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TMPU2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TMPU2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TMPU2', State_Met%TMPU2, &
                            State_Met, RC )

#endif

    !=======================================================================
    ! Allocate land type and leaf area index fields for dry deposition
    !=======================================================================

    !-------------------------
    ! IREG [1]
    !-------------------------
    ALLOCATE( State_Met%IREG( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%IREG', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%IREG = 0
    CALL Register_MetField( am_I_Root, 'IREG', State_Met%IREG, &
                            State_Met, RC )

    !-------------------------
    ! ILAND [1]
    !-------------------------
    ALLOCATE( State_Met%ILAND( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%ILAND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ILAND = 0
    CALL Register_MetField( am_I_Root, 'ILAND', State_Met%ILAND, &
                            State_Met, RC )

    !-------------------------
    ! IUSE [1]
    !-------------------------
    ALLOCATE( State_Met%IUSE( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%IUSE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%IUSE = 0
    CALL Register_MetField( am_I_Root, 'IUSE', State_Met%IUSE, &
                            State_Met, RC )

    !-------------------------
    ! XLAI [1]
    !-------------------------
    ALLOCATE( State_Met%XLAI( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XLAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAI', State_Met%XLAI, &
                            State_Met, RC )

    !-------------------------
    ! MODISLAI [1]
    !-------------------------
    ALLOCATE( State_Met%MODISLAI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%MODISLAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%MODISLAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MODISLAI', State_Met%MODISLAI, &
                            State_Met, RC )

    !-------------------------
    ! XCHLR [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%XCHLR( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XCHLR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLR', State_Met%XCHLR, &
                            State_Met, RC )

    !-------------------------
    ! MODISCHLR [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%MODISCHLR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%MODISCHLR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%MODISCHLR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MODISCHLR', State_Met%MODISCHLR, &
                            State_Met, RC )

    !-------------------------
    ! LANDTYPEFRAC [1]
    !-------------------------    
    ALLOCATE( State_Met%LANDTYPEFRAC( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%LANDTYPEFRAC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LANDTYPEFRAC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LANDTYPEFRAC', State_Met%LANDTYPEFRAC, &
                            State_Met, RC )

    !-------------------------
    ! XLAI_NATIVE [1]
    !-------------------------
    ALLOCATE( State_Met%XLAI_NATIVE( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XLAI_NATIVE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI_NATIVE  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAI_NATIVE', State_Met%XLAI_NATIVE, &
                            State_Met, RC )

    !-------------------------
    ! XCHLR_NATIVE [1]
    !-------------------------
    ALLOCATE( State_Met%XCHLR_NATIVE( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XCHLR_NATIVE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR_NATIVE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLR_NATIVE', State_Met%XCHLR_NATIVE, &
                            State_Met, RC )

    !-------------------------
    ! XLAI2 [m2 m-2]
    !-------------------------
    ALLOCATE( State_Met%XLAI2( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XLAI2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAI2', State_Met%XLAI2, &
                            State_Met, RC )

    !-------------------------
    ! XCHLR2 [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%XCHLR2( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XCHLR2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLR2', State_Met%XCHLR2, &
                            State_Met, RC )

    !=======================================================================
    ! Print information about the registered fields (short format)
    !=======================================================================
    CALL Print_State_Met( am_I_Root, State_Met, RC, ShortFormat=.TRUE.)
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Print_State_Met"'
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
  SUBROUTINE Cleanup_State_Met( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Destroy
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  23 Oct 2012 - R. Yantosca - Now deallocate QI, QL fields
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  19 Nov 2012 - R. Yantosca - Segregate DEALLOCATE statements w/ #ifdefs
!                              for each met field data product type
!  27 Nov 2012 - R. Yantosca - Now deallocate the SUNCOS fields
!  12 Dec 2012 - R. Yantosca - Now deallocate the IREG, ILAND, IUSE fields
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  22 Aug 2014 - R. Yantosca - Deallocate PBL_TOP_L field
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
    ! These met fields are used for all data products
    !========================================================================

    ! 2-D fields
    IF ( ASSOCIATED( State_Met%ALBD       )) DEALLOCATE( State_Met%ALBD       ) 
    IF ( ASSOCIATED( State_Met%CLDFRC     )) DEALLOCATE( State_Met%CLDFRC     ) 
    IF ( ASSOCIATED( State_Met%CLDTOPS    )) DEALLOCATE( State_Met%CLDTOPS    )
    IF ( ASSOCIATED( State_Met%EFLUX      )) DEALLOCATE( State_Met%EFLUX      ) 
    IF ( ASSOCIATED( State_Met%EVAP       )) DEALLOCATE( State_Met%EVAP       )
    IF ( ASSOCIATED( State_Met%FRCLND     )) DEALLOCATE( State_Met%FRCLND     )
    IF ( ASSOCIATED( State_Met%FRLAKE     )) DEALLOCATE( State_Met%FRLAKE     )
    IF ( ASSOCIATED( State_Met%FRLAND     )) DEALLOCATE( State_Met%FRLAND     )
    IF ( ASSOCIATED( State_Met%FRLANDIC   )) DEALLOCATE( State_Met%FRLANDIC   )
    IF ( ASSOCIATED( State_Met%FROCEAN    )) DEALLOCATE( State_Met%FROCEAN    )
    IF ( ASSOCIATED( State_Met%GRN        )) DEALLOCATE( State_Met%GRN        )
    IF ( ASSOCIATED( State_Met%GWETROOT   )) DEALLOCATE( State_Met%GWETROOT   )
    IF ( ASSOCIATED( State_Met%GWETTOP    )) DEALLOCATE( State_Met%GWETTOP    )
    IF ( ASSOCIATED( State_Met%HFLUX      )) DEALLOCATE( State_Met%HFLUX      )
    IF ( ASSOCIATED( State_Met%LAI        )) DEALLOCATE( State_Met%LAI        )
    IF ( ASSOCIATED( State_Met%LWI        )) DEALLOCATE( State_Met%LWI        )
    IF ( ASSOCIATED( State_Met%PARDR      )) DEALLOCATE( State_Met%PARDR      )
    IF ( ASSOCIATED( State_Met%PARDF      )) DEALLOCATE( State_Met%PARDF      )
    IF ( ASSOCIATED( State_Met%PBLH       )) DEALLOCATE( State_Met%PBLH       )
    IF ( ASSOCIATED( State_Met%PBL_TOP_L  )) DEALLOCATE( State_Met%PBL_TOP_L  )
    IF ( ASSOCIATED( State_Met%PHIS       )) DEALLOCATE( State_Met%PHIS       )
    IF ( ASSOCIATED( State_Met%PRECCON    )) DEALLOCATE( State_Met%PRECCON    )
    IF ( ASSOCIATED( State_Met%PRECTOT    )) DEALLOCATE( State_Met%PRECTOT    )
    IF ( ASSOCIATED( State_Met%PRECSNO    )) DEALLOCATE( State_Met%PRECSNO    )
    IF ( ASSOCIATED( State_Met%PS1_WET    )) DEALLOCATE( State_Met%PS1_WET    )
    IF ( ASSOCIATED( State_Met%PS2_WET    )) DEALLOCATE( State_Met%PS2_WET    )
    IF ( ASSOCIATED( State_Met%PSC2_WET   )) DEALLOCATE( State_Met%PSC2_WET   )
    IF ( ASSOCIATED( State_Met%PS1_DRY    )) DEALLOCATE( State_Met%PS1_DRY    )
    IF ( ASSOCIATED( State_Met%PS2_DRY    )) DEALLOCATE( State_Met%PS2_DRY    )
    IF ( ASSOCIATED( State_Met%PSC2_DRY   )) DEALLOCATE( State_Met%PSC2_DRY   )
    IF ( ASSOCIATED( State_Met%RADLWG     )) DEALLOCATE( State_Met%RADLWG     )
    IF ( ASSOCIATED( State_Met%RADSWG     )) DEALLOCATE( State_Met%RADSWG     )
    IF ( ASSOCIATED( State_Met%SLP        )) DEALLOCATE( State_Met%SLP        )
    IF ( ASSOCIATED( State_Met%SNODP      )) DEALLOCATE( State_Met%SNODP      )
    IF ( ASSOCIATED( State_Met%SNOMAS     )) DEALLOCATE( State_Met%SNOMAS     )
    IF ( ASSOCIATED( State_Met%SST        )) DEALLOCATE( State_Met%SST        )
    IF ( ASSOCIATED( State_Met%SUNCOS     )) DEALLOCATE( State_Met%SUNCOS     )
    IF ( ASSOCIATED( State_Met%SUNCOSmid  )) DEALLOCATE( State_Met%SUNCOSmid  )
!    IF ( ASSOCIATED( State_Met%SUNCOSmid5 )) DEALLOCATE( State_Met%SUNCOSmid5 )
    IF ( ASSOCIATED( State_Met%SWGDN      )) DEALLOCATE( State_Met%SWGDN      )
    IF ( ASSOCIATED( State_Met%TROPP      )) DEALLOCATE( State_Met%TROPP      )
    IF ( ASSOCIATED( State_Met%TS         )) DEALLOCATE( State_Met%TS         )
    IF ( ASSOCIATED( State_Met%TSKIN      )) DEALLOCATE( State_Met%TSKIN      )
    IF ( ASSOCIATED( State_Met%TO3        )) DEALLOCATE( State_Met%TO3        )
    IF ( ASSOCIATED( State_Met%U10M       )) DEALLOCATE( State_Met%U10M       )
    IF ( ASSOCIATED( State_Met%USTAR      )) DEALLOCATE( State_Met%USTAR      )
    IF ( ASSOCIATED( State_Met%UVALBEDO   )) DEALLOCATE( State_Met%UVALBEDO   )
    IF ( ASSOCIATED( State_Met%V10M       )) DEALLOCATE( State_Met%V10M       )
    IF ( ASSOCIATED( State_Met%Z0         )) DEALLOCATE( State_Met%Z0         )
    IF ( ASSOCIATED( State_Met%CNV_FRC    )) DEALLOCATE( State_Met%CNV_FRC    )

#if defined( ESMF_ ) 

    !=========================================================================
    ! SDE 2016-03-28: GCHP requires that these are nullified rather than being
    ! deallocated. Not yet sure why, but deallocating causes it to hang during
    ! cleanup
    !=========================================================================

    ! 3-D fields
    IF ( ASSOCIATED( State_Met%AD         )) NULLIFY( State_Met%AD         )
    IF ( ASSOCIATED( State_Met%AIRDEN     )) NULLIFY( State_Met%AIRDEN     )
    IF ( ASSOCIATED( State_Met%AIRVOL     )) NULLIFY( State_Met%AIRVOL     )
    IF ( ASSOCIATED( State_Met%AREA_M2    )) NULLIFY( State_Met%AREA_M2    )
    IF ( ASSOCIATED( State_Met%AVGW       )) NULLIFY( State_Met%AVGW       )
    IF ( ASSOCIATED( State_Met%BXHEIGHT   )) NULLIFY( State_Met%BXHEIGHT   )
    IF ( ASSOCIATED( State_Met%CLDF       )) NULLIFY( State_Met%CLDF       )
    IF ( ASSOCIATED( State_Met%CMFMC      )) NULLIFY( State_Met%CMFMC      )
    IF ( ASSOCIATED( State_Met%DELP       )) NULLIFY( State_Met%DELP       )
    IF ( ASSOCIATED( State_Met%DELP_PREV  )) NULLIFY( State_Met%DELP_PREV  )
    IF ( ASSOCIATED( State_Met%DQRCU      )) NULLIFY( State_Met%DQRCU      )
    IF ( ASSOCIATED( State_Met%DQRLSAN    )) NULLIFY( State_Met%DQRLSAN    )
    IF ( ASSOCIATED( State_Met%DQIDTMST   )) NULLIFY( State_Met%DQIDTMST   )
    IF ( ASSOCIATED( State_Met%DQLDTMST   )) NULLIFY( State_Met%DQLDTMST   )
    IF ( ASSOCIATED( State_Met%DQVDTMST   )) NULLIFY( State_Met%DQVDTMST   )
    IF ( ASSOCIATED( State_Met%DTRAIN     )) NULLIFY( State_Met%DTRAIN     )
    IF ( ASSOCIATED( State_Met%MOISTQ     )) NULLIFY( State_Met%MOISTQ     )
    IF ( ASSOCIATED( State_Met%OMEGA      )) NULLIFY( State_Met%OMEGA      )
    IF ( ASSOCIATED( State_Met%OPTD       )) NULLIFY( State_Met%OPTD       )
    IF ( ASSOCIATED( State_Met%PEDGE      )) NULLIFY( State_Met%PEDGE      )
    IF ( ASSOCIATED( State_Met%PEDGE_DRY  )) NULLIFY( State_Met%PEDGE_DRY  )
    IF ( ASSOCIATED( State_Met%PMID       )) NULLIFY( State_Met%PMID       )
    IF ( ASSOCIATED( State_Met%PMID_DRY   )) NULLIFY( State_Met%PMID_DRY   )
    IF ( ASSOCIATED( State_Met%PV         )) NULLIFY( State_Met%PV         )
    IF ( ASSOCIATED( State_Met%QI         )) NULLIFY( State_Met%QI         )
    IF ( ASSOCIATED( State_Met%QL         )) NULLIFY( State_Met%QL         )
    IF ( ASSOCIATED( State_Met%RH         )) NULLIFY( State_Met%RH         )
    IF ( ASSOCIATED( State_Met%SPHU       )) NULLIFY( State_Met%SPHU       )
    IF ( ASSOCIATED( State_Met%SPHU_PREV  )) NULLIFY( State_Met%SPHU_PREV  )
    IF ( ASSOCIATED( State_Met%T          )) NULLIFY( State_Met%T          )
    IF ( ASSOCIATED( State_Met%TV         )) NULLIFY( State_Met%TV         )
    IF ( ASSOCIATED( State_Met%TAUCLI     )) NULLIFY( State_Met%TAUCLI     )
    IF ( ASSOCIATED( State_Met%TAUCLW     )) NULLIFY( State_Met%TAUCLW     ) 
    IF ( ASSOCIATED( State_Met%U          )) NULLIFY( State_Met%U          )
    IF ( ASSOCIATED( State_Met%UPDVVEL    )) NULLIFY( State_Met%UPDVVEL    )
    IF ( ASSOCIATED( State_Met%V          )) NULLIFY( State_Met%V          )

#else

    !=========================================================================
    ! For GEOS-Chem "Classic" simulations, we should use DEALLOCATE instead
    ! of NULLIFY.  Typically if you allocate memory to a pointer-based
    ! variable, you also need to DEALLOCATE.  If a pointer points to a
    ! target variable, then you should use NULLIFY. (bmy, 6/27/16)
    !=========================================================================

    ! 3-D fields
    IF ( ASSOCIATED( State_Met%AD         )) DEALLOCATE( State_Met%AD         )
    IF ( ASSOCIATED( State_Met%AIRDEN     )) DEALLOCATE( State_Met%AIRDEN     )
    IF ( ASSOCIATED( State_Met%MAIRDEN    )) DEALLOCATE( State_Met%MAIRDEN    )
    IF ( ASSOCIATED( State_Met%AIRVOL     )) DEALLOCATE( State_Met%AIRVOL     )
    IF ( ASSOCIATED( State_Met%AREA_M2    )) DEALLOCATE( State_Met%AREA_M2    )
    IF ( ASSOCIATED( State_Met%AVGW       )) DEALLOCATE( State_Met%AVGW       )
    IF ( ASSOCIATED( State_Met%BXHEIGHT   )) DEALLOCATE( State_Met%BXHEIGHT   )
    IF ( ASSOCIATED( State_Met%CLDF       )) DEALLOCATE( State_Met%CLDF       )
    IF ( ASSOCIATED( State_Met%CMFMC      )) DEALLOCATE( State_Met%CMFMC      )
    IF ( ASSOCIATED( State_Met%DELP       )) DEALLOCATE( State_Met%DELP       )
    IF ( ASSOCIATED( State_Met%DELP_DRY   )) DEALLOCATE( State_Met%DELP_DRY   )
    IF ( ASSOCIATED( State_Met%DELP_PREV  )) DEALLOCATE( State_Met%DELP_PREV  )
    IF ( ASSOCIATED( State_Met%DP_DRY_PREV)) DEALLOCATE( State_Met%DP_DRY_PREV)
    IF ( ASSOCIATED( State_Met%DQRCU      )) DEALLOCATE( State_Met%DQRCU      )
    IF ( ASSOCIATED( State_Met%DQRLSAN    )) DEALLOCATE( State_Met%DQRLSAN    )
    IF ( ASSOCIATED( State_Met%DQIDTMST   )) DEALLOCATE( State_Met%DQIDTMST   )
    IF ( ASSOCIATED( State_Met%DQLDTMST   )) DEALLOCATE( State_Met%DQLDTMST   )
    IF ( ASSOCIATED( State_Met%DQVDTMST   )) DEALLOCATE( State_Met%DQVDTMST   )
    IF ( ASSOCIATED( State_Met%DTRAIN     )) DEALLOCATE( State_Met%DTRAIN     )
    IF ( ASSOCIATED( State_Met%MOISTQ     )) DEALLOCATE( State_Met%MOISTQ     )
    IF ( ASSOCIATED( State_Met%OMEGA      )) DEALLOCATE( State_Met%OMEGA      )
    IF ( ASSOCIATED( State_Met%OPTD       )) DEALLOCATE( State_Met%OPTD       )
    IF ( ASSOCIATED( State_Met%PEDGE      )) DEALLOCATE( State_Met%PEDGE      )
    IF ( ASSOCIATED( State_Met%PEDGE_DRY  )) DEALLOCATE( State_Met%PEDGE_DRY  )
    IF ( ASSOCIATED( State_Met%PMID       )) DEALLOCATE( State_Met%PMID       )
    IF ( ASSOCIATED( State_Met%PMID_DRY   )) DEALLOCATE( State_Met%PMID_DRY   )
    IF ( ASSOCIATED( State_Met%PV         )) DEALLOCATE( State_Met%PV         )
    IF ( ASSOCIATED( State_Met%QI         )) DEALLOCATE( State_Met%QI         )
    IF ( ASSOCIATED( State_Met%QL         )) DEALLOCATE( State_Met%QL         )
    IF ( ASSOCIATED( State_Met%RH         )) DEALLOCATE( State_Met%RH         )
    IF ( ASSOCIATED( State_Met%SPHU       )) DEALLOCATE( State_Met%SPHU       )
    IF ( ASSOCIATED( State_Met%SPHU_PREV  )) DEALLOCATE( State_Met%SPHU_PREV  )
    IF ( ASSOCIATED( State_Met%T          )) DEALLOCATE( State_Met%T          )
    IF ( ASSOCIATED( State_Met%TV         )) DEALLOCATE( State_Met%TV         )
    IF ( ASSOCIATED( State_Met%TAUCLI     )) DEALLOCATE( State_Met%TAUCLI     )
    IF ( ASSOCIATED( State_Met%TAUCLW     )) DEALLOCATE( State_Met%TAUCLW     ) 
    IF ( ASSOCIATED( State_Met%U          )) DEALLOCATE( State_Met%U          )
    IF ( ASSOCIATED( State_Met%UPDVVEL    )) DEALLOCATE( State_Met%UPDVVEL    )
    IF ( ASSOCIATED( State_Met%V          )) DEALLOCATE( State_Met%V          )

#endif

#if defined( GCAP )
    !========================================================================
    ! Fields specific to the GCAP met data product
    !========================================================================

    ! 2-D fields
    IF ( ASSOCIATED( State_Met%LWI_GISS   )) DEALLOCATE( State_Met%LWI_GISS   )
    IF ( ASSOCIATED( State_Met%MOLENGTH   )) DEALLOCATE( State_Met%MOLENGTH   )
    IF ( ASSOCIATED( State_Met%OICE       )) DEALLOCATE( State_Met%OICE       )
    IF ( ASSOCIATED( State_Met%SNICE      )) DEALLOCATE( State_Met%SNICE      )
    IF ( ASSOCIATED( State_Met%SNOW       )) DEALLOCATE( State_Met%SNOW       )
    IF ( ASSOCIATED( State_Met%TROPP1     )) DEALLOCATE( State_Met%TROPP1     )
    IF ( ASSOCIATED( State_Met%TROPP2     )) DEALLOCATE( State_Met%TROPP2     )

    ! 3-D fields
    IF ( ASSOCIATED( State_Met%DETRAINE   )) DEALLOCATE( State_Met%DETRAINE   )
    IF ( ASSOCIATED( State_Met%DETRAINN   )) DEALLOCATE( State_Met%DETRAINN   )
    IF ( ASSOCIATED( State_Met%DNDE       )) DEALLOCATE( State_Met%DNDE       )
    IF ( ASSOCIATED( State_Met%DNDN       )) DEALLOCATE( State_Met%DNDN       )
    IF ( ASSOCIATED( State_Met%ENTRAIN    )) DEALLOCATE( State_Met%ENTRAIN    )
    IF ( ASSOCIATED( State_Met%UPDE       )) DEALLOCATE( State_Met%UPDE       )
    IF ( ASSOCIATED( State_Met%UPDN       )) DEALLOCATE( State_Met%UPDN       )

#elif defined( GEOS_4 )
    !========================================================================
    ! Fields specific to the GMAO GEOS-4 met data product
    !========================================================================

    ! 2-D fields
    IF ( ASSOCIATED( State_Met%SNOW       )) DEALLOCATE( State_Met%SNOW       )
    IF ( ASSOCIATED( State_Met%TROPP1     )) DEALLOCATE( State_Met%TROPP1     )
    IF ( ASSOCIATED( State_Met%TROPP2     )) DEALLOCATE( State_Met%TROPP2     )
                                          
    ! 3-D fields                          
    IF ( ASSOCIATED( State_Met%HKBETA     )) DEALLOCATE( State_Met%HKBETA     )
    IF ( ASSOCIATED( State_Met%HKETA      )) DEALLOCATE( State_Met%HKETA      )
    IF ( ASSOCIATED( State_Met%ZMEU       )) DEALLOCATE( State_Met%ZMEU       )
    IF ( ASSOCIATED( State_Met%ZMMD       )) DEALLOCATE( State_Met%ZMMD       )
    IF ( ASSOCIATED( State_Met%ZMMU       )) DEALLOCATE( State_Met%ZMMU       )

#elif defined( GEOS_5 )
    !========================================================================
    ! Fields specific to the GMAO GEOS-5 met data product
    !========================================================================

    ! 2-D fields
    IF ( ASSOCIATED( State_Met%TO31       )) DEALLOCATE( State_Met%TO31       )
    IF ( ASSOCIATED( State_Met%TO32       )) DEALLOCATE( State_Met%TO32       )
    IF ( ASSOCIATED( State_Met%TTO3       )) DEALLOCATE( State_Met%TTO3       )

#elif defined( GEOS_FP ) || defined( MERRA ) || defined( MERRA2 )
    !========================================================================
    ! Fields specific to the GMAO GEOS-FP, MERRA, and MERRA-2 met products
    !========================================================================

    ! 2-D fields
    IF ( ASSOCIATED( State_Met%FRSEAICE   )) DEALLOCATE( State_Met%FRSEAICE   )
    IF ( ASSOCIATED( State_Met%FRSNO      )) DEALLOCATE( State_Met%FRSNO      )
    IF ( ASSOCIATED( State_Met%PRECANV    )) DEALLOCATE( State_Met%PRECANV    )
    IF ( ASSOCIATED( State_Met%PRECLSC    )) DEALLOCATE( State_Met%PRECLSC    )
    IF ( ASSOCIATED( State_Met%SEAICE00   )) DEALLOCATE( State_Met%SEAICE00   )
    IF ( ASSOCIATED( State_Met%SEAICE10   )) DEALLOCATE( State_Met%SEAICE10   )
    IF ( ASSOCIATED( State_Met%SEAICE20   )) DEALLOCATE( State_Met%SEAICE20   )
    IF ( ASSOCIATED( State_Met%SEAICE30   )) DEALLOCATE( State_Met%SEAICE30   )
    IF ( ASSOCIATED( State_Met%SEAICE40   )) DEALLOCATE( State_Met%SEAICE40   )
    IF ( ASSOCIATED( State_Met%SEAICE50   )) DEALLOCATE( State_Met%SEAICE50   )
    IF ( ASSOCIATED( State_Met%SEAICE60   )) DEALLOCATE( State_Met%SEAICE60   )
    IF ( ASSOCIATED( State_Met%SEAICE70   )) DEALLOCATE( State_Met%SEAICE70   )
    IF ( ASSOCIATED( State_Met%SEAICE80   )) DEALLOCATE( State_Met%SEAICE80   )
    IF ( ASSOCIATED( State_Met%SEAICE90   )) DEALLOCATE( State_Met%SEAICE90   )

    ! 3-D fields
    IF ( ASSOCIATED( State_Met%PFICU      )) DEALLOCATE( State_Met%PFICU      )
    IF ( ASSOCIATED( State_Met%PFILSAN    )) DEALLOCATE( State_Met%PFILSAN    )
    IF ( ASSOCIATED( State_Met%PFLCU      )) DEALLOCATE( State_Met%PFLCU      )
    IF ( ASSOCIATED( State_Met%PFLLSAN    )) DEALLOCATE( State_Met%PFLLSAN    )
    IF ( ASSOCIATED( State_Met%REEVAPCN   )) DEALLOCATE( State_Met%REEVAPCN   )
    IF ( ASSOCIATED( State_Met%REEVAPLS   )) DEALLOCATE( State_Met%REEVAPLS   )
    IF ( ASSOCIATED( State_Met%RH1        )) DEALLOCATE( State_Met%RH1        )
    IF ( ASSOCIATED( State_Met%RH2        )) DEALLOCATE( State_Met%RH2        )
    IF ( ASSOCIATED( State_Met%SPHU1      )) DEALLOCATE( State_Met%SPHU1      )
    IF ( ASSOCIATED( State_Met%SPHU2      )) DEALLOCATE( State_Met%SPHU2      )
    IF ( ASSOCIATED( State_Met%TMPU1      )) DEALLOCATE( State_Met%TMPU1      )
    IF ( ASSOCIATED( State_Met%TMPU2      )) DEALLOCATE( State_Met%TMPU2      )

#endif

    !========================================================================
    ! Land type and leaf area index (LAI) fields for dry deposition
    !========================================================================
    IF ( ASSOCIATED( State_Met%IREG       )) DEALLOCATE( State_Met%IREG       )
    IF ( ASSOCIATED( State_Met%ILAND      )) DEALLOCATE( State_Met%ILAND      )
    IF ( ASSOCIATED( State_Met%IUSE       )) DEALLOCATE( State_Met%IUSE       )
    IF ( ASSOCIATED( State_Met%XLAI       )) DEALLOCATE( State_Met%XLAI       )
    IF ( ASSOCIATED( State_Met%XLAI2      )) DEALLOCATE( State_Met%XLAI2      )
    IF ( ASSOCIATED( State_Met%MODISLAI   )) DEALLOCATE( State_Met%MODISLAI   )
    IF ( ASSOCIATED( State_Met%XCHLR      )) DEALLOCATE( State_Met%XCHLR      )
    IF ( ASSOCIATED( State_Met%XCHLR2     )) DEALLOCATE( State_Met%XCHLR2     )
    IF ( ASSOCIATED( State_Met%MODISCHLR  )) DEALLOCATE( State_Met%MODISCHLR  )
    IF (ASSOCIATED( State_Met%LANDTYPEFRAC)) DEALLOCATE( State_Met%LANDTYPEFRAC)
    IF (ASSOCIATED( State_Met%XLAI_NATIVE )) DEALLOCATE( State_Met%XLAI_NATIVE )
    IF (ASSOCIATED( State_Met%XCHLR_NATIVE)) DEALLOCATE( State_Met%XCHLR_NATIVE)

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, State_Met%Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Met%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_State_Met
!
! !DESCRIPTION: Print information about all the registered variables
!  contained within the State\_Met object.  This is basically a wrapper for
!  routine REGISTRY\_PRINT in registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_State_Met( am_I_Root, State_Met, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Print
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Root CPU?  
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
    LOGICAL,        OPTIONAL    :: ShortFormat ! Print truncated info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success/failure?
!
! !REVISION HISTORY:
!  29 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Print_State_Met (in Headers/state_met_mod.F90)'

    !=======================================================================
    ! Print info about registered variables
    !=======================================================================

    ! Header line
    if ( am_I_Root ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the ' // &
               'State_Met object:')
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    ! Print registry info in truncated format
    CALL Registry_Print( am_I_Root   = am_I_Root,           &
                         Registry    = State_Met%Registry,  &
                         ShortFormat = ShortFormat,         &
                         RC          = RC                  )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Print"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Print_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lookup_State_Met
!
! !DESCRIPTION: Return metadata and/or a pointer to the data for any
!  variable contained within the State\_Met object by searching for its name.
!  This is basically a wrapper for routine REGISTRY\_LOOKUP in 
!  registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Lookup_State_Met( am_I_Root, State_Met,    Variable,            &
                               RC,        Description,  Dimensions,          &
                               KindVal,   MemoryInKb,   Rank,                &
                               Units,     OnLevelEdges, Ptr2d,               &
                               Ptr3d,     Ptr2d_I,      Ptr3d_I             )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Lookup
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)  :: am_I_Root       ! Is this the root CPU? 
    TYPE(MetState),   INTENT(IN)  :: State_Met       ! Meteorology State
    CHARACTER(LEN=*), INTENT(IN)  :: Variable        ! Variable name
!
! !OUTPUT PARAMETERS:
!
    ! Required outputs
    INTEGER,          INTENT(OUT) :: RC              ! Success or failure?

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL :: Description     ! Description of data
    INTEGER,             OPTIONAL :: Dimensions(3)   ! Dimensions of data
    INTEGER,             OPTIONAL :: KindVal         ! Numerical KIND value
    REAL(fp),            OPTIONAL :: MemoryInKb      ! Memory usage
    INTEGER,             OPTIONAL :: Rank            ! Size of data
    CHARACTER(LEN=255),  OPTIONAL :: Units           ! Units of data
    LOGICAL,             OPTIONAL :: OnLevelEdges    ! =T if data is defined
                                                     !    on level edges
                                                     ! =F if on centers

    ! Pointers to data
    REAL(fp),   POINTER, OPTIONAL :: Ptr2d  (:,:  )  ! 2D flex-prec data
    REAL(fp),   POINTER, OPTIONAL :: Ptr3d  (:,:,:)  ! 3D flex-prec data
    INTEGER,    POINTER, OPTIONAL :: Ptr2d_I(:,:  )  ! 2D integer data
    INTEGER,    POINTER, OPTIONAL :: Ptr3d_I(:,:,:)  ! 3D integer data
!
! !REMARKS:
!  We keep the StateName variable private to this module. Users only have
!  to supply the name of each module variable.
!
! !REVISION HISTORY:
!  29 Jun 2017 - R. Yantosca - Initial version
!  30 Jun 2017 - R. Yantosca - Rename variables Ptr{2,3}*dI to Ptr{2,3}d_I
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Lookup_State_Met (in Headers/state_met_mod.F90)'

    !=======================================================================
    ! Look up a variable; Return metadata and/or a pointer to the data
    !=======================================================================
    CALL Registry_Lookup( am_I_Root    = am_I_Root,                         &
                          Registry     = State_Met%Registry,                &
                          State        = State_Met%State,                   &
                          Variable     = Variable,                          &
                          Description  = Description,                       &
                          Dimensions   = Dimensions,                        &
                          KindVal      = KindVal,                           &
                          MemoryInKb   = MemoryInKb,                        &
                          Rank         = Rank,                              &
                          Units        = Units,                             &
                          OnLevelEdges = OnLevelEdges,                      &
                          Ptr2d        = Ptr2d,                             &
                          Ptr3d        = Ptr3d,                             &
                          Ptr2d_I      = Ptr2d_I,                           &
                          Ptr3d_I      = Ptr3d_I,                           &
                          RC           = RC                                )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable "' // TRIM( Variable ) // &
               '" in the State_Met registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_MetField_Metadata
!
! !DESCRIPTION: Subroutine GET\_METFIELD\_METADATA retrieves basic 
!  information about each State_Met field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_MetField_Metadata( am_I_Root, Name,  Desc, Units, &
                                    Rank,      Type,  VLoc, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
    CHARACTER(LEN=*),    INTENT(IN)  :: Name       ! Sate_Met field name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Return code

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL    :: Desc       ! Long name string
    CHARACTER(LEN=255),  OPTIONAL    :: Units      ! Units string
    INTEGER,             OPTIONAL    :: Rank       ! # of dimensions
    INTEGER,             OPTIONAL    :: Type       ! Desc of data type
    INTEGER,             OPTIONAL    :: VLoc       ! Vertical placement
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  28 Aug 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc
    
    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_MetField_Metadata (in Headers/state_met_mod.F90)'

    ! Optional arguments present?
    isDesc  = PRESENT( Desc  )
    isUnits = PRESENT( Units )
    isRank  = PRESENT( Rank  )
    isType  = PRESENT( Type  )
    isVLoc  = PRESENT( VLoc  )

    ! Set defaults for optional arguments
    IF ( isUnits ) Units = ''
    IF ( isDesc  ) Desc  = ''
    IF ( isRank  ) Rank  = 0
    IF ( isType  ) Type  = KINDVAL_FP ! Is this always true if real?
    IF ( isVLoc  ) VLoc  = VLocCenter ! Assume centered

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM(Name) )

       CASE ( 'ALBD' )
          IF ( isDesc  ) Desc  = 'Visible surface albedo'
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

       CASE ( 'EFLUX' )
          IF ( isDesc  ) Desc  = 'Latent heat flux'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'EVAP' )
          IF ( isDesc  ) Desc  = 'Surface evaporation'
          IF ( isUnits ) Units = 'kg m-2 s-1'
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

       CASE ( 'GRN' )
          IF ( isDesc  ) Desc  = 'Greenness fraction'
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

       CASE ( 'ITY' )
          IF ( isDesc  ) Desc  = 'Land surface type index'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'LWI' )
          IF ( isDesc  ) Desc  = 'Land-water-ice indices'
          IF ( isUnits ) Units = '1'
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

       CASE ( 'PBL_TOP_L' )
          IF ( isDesc  ) Desc  = 'Model layer of the planetary boundary ' // &
                                 'layer top occurs'
          IF ( isUnits ) Units = 'layer'
          IF ( isRank  ) Rank  = 2
          IF ( isType  ) Type  = KINDVAL_I4

       CASE ( 'PHIS' )
          IF ( isDesc  ) Desc  = 'Surface geopotential height'
          IF ( isUnits ) Units = 'm2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECCON' )
          IF ( isDesc  ) Desc  = 'Convective precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECSNO' )
          IF ( isDesc  ) Desc  = 'Snow precipitation'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECTOT' )
          IF ( isDesc  ) Desc  = 'Total precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS1_WET' )
          IF ( isDesc  ) Desc  = 'Wet surface pressure at dt start'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS2_WET' )
          IF ( isDesc  ) Desc  = 'Wet surface pressure at dt end'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PSC2_WET' )
          IF ( isDesc  ) Desc  = 'Wet interpolated surface pressure'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS1_DRY' )
          IF ( isDesc  ) Desc  = 'Dry surface pressure at dt start'
          IF ( isUnits ) Units = ''
          IF ( isRank  ) Rank  = 2

       CASE ( 'PS2_DRY' )
          IF ( isDesc  ) Desc  = 'Dry surface pressure at dt end'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PSC2_DRY' )
          IF ( isDesc  ) Desc  = 'Dry interpolated surface pressure'
          IF ( isUnits ) Units = 'hPA'
          IF ( isRank  ) Rank  = 2

       CASE ( 'RADLWG' )
          IF ( isDesc  ) Desc  = 'Net longwave radiation at ground'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2
 
       CASE ( 'RADSWG' )
          IF ( isDesc  ) Desc  = 'Shortwave radiation at ground'
          IF ( isUnits ) Units = 'W m-2'
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

       CASE ( 'SUNCOSmid' )
          IF ( isDesc  ) Desc  = 'Cosine of solar zenith angle, at ' // &
                                 'midpoint of chemistry timestep'
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

#if defined( ESMF_ )
       CASE ( 'CNV_FRC' )
          IF ( isDesc  ) Desc  = 'Convective fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

#endif
       CASE ( 'FRESEAICE' )
          IF ( isDesc  ) Desc  = 'Fraction of sea ice'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'FRSNO' )
          IF ( isDesc  ) Desc  = 'Fraction of snow on surface'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECANV' )
          IF ( isDesc  ) Desc  = 'Anvil precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'PRECLSC' )
          IF ( isDesc  ) Desc  = 'Large-scale precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEASICE00' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 00-10%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEASICE10' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 10-20%'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SEAICE20' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 20-30%'
          IF ( isUnits ) Units = '1'

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

       CASE ( 'AD' )
          IF ( isDesc  ) Desc  = 'Dry air mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AIRDEN' )
          IF ( isDesc  ) Desc  = 'Dry air density'
          IF ( isUnits ) Units = 'kg m-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'MAIRDEN' )
          IF ( isDesc  ) Desc  = 'Moist air density'
          IF ( isUnits ) Units = 'kg m-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AIRNUMDEN' )
          IF ( isDesc  ) Desc  = 'Dry air density'
          IF ( isUnits ) Units = 'm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AIRVOL' )
          IF ( isDesc  ) Desc  = 'Volume of dry air in grid box'
          IF ( isUnits ) Units = 'm3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AREA_M2' )
          IF ( isDesc  ) Desc  = 'Surface area of grid box'
          IF ( isUnits ) Units = 'm2'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AVGW' )
          IF ( isDesc  ) Desc  = 'Water vapor mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'vol vol-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'BXHEIGHT' )
          IF ( isDesc  ) Desc  = 'Grid box height (w/r/t dry air)'
          IF ( isUnits ) Units = 'm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'CLDF' )
          IF ( isDesc  ) Desc  = '3-D cloud fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'CMFMC' )
          IF ( isDesc  ) Desc  = 'Cloud mass flux'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocEdge

       CASE ( 'DELP' )
          IF ( isDesc  ) Desc  = 'Delta-pressure across grid box(wet air)'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DELP_DRY' )
          IF ( isDesc  ) Desc  = 'Delta-pressure across grid box (dry air)'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DELP_PREV' )
          IF ( isDesc  ) Desc  = 'Previous State_Met%DELP'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DP_DRY_PREV' )
          IF ( isDesc  ) Desc  = 'Previous State_Met%DELP_DRY'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DQRCU' )
          IF ( isDesc  ) Desc  = 'Production rate of convective ' // &
                                 'precipitation (per dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DQRLSAN' )
          IF ( isDesc  ) Desc  = 'Production rate of large-scale ' // &
                                 'precipitation (per dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DQIDTMST' )
          IF ( isDesc  ) Desc  = 'Ice tendency from moist processes'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DQLDTMST' )
          IF ( isDesc  ) Desc  = 'H2O tendency from moist processes'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DQVDTMST' )
          IF ( isDesc  ) Desc  = 'Vapor tendency from moist processes'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'DTRAIN' )
          IF ( isDesc  ) Desc  = 'Detrainment flux'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'MOISTQ' )
          IF ( isDesc  ) Desc  = 'Tendency in specific humidity'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'OMEGA' )
          IF ( isDesc  ) Desc  = 'Updraft velocity'
          IF ( isUnits ) Units = 'Pa s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'OPTD' )
          IF ( isDesc  ) Desc  = 'Visible optical depth'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'PEDGE' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t moist air) at level edges'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'PEDGE_DRY' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t dry air) at level edges'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'PMID' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t moist air) at level centers'
          IF ( isUnits ) Units = 'hPa'
          IF ( isRank  ) Rank  = 3

       CASE ( 'PMID_DRY' )
          IF ( isDesc  ) Desc  = 'Pressure (w/r/t dry air) at level centers'
          IF ( isUnits ) Units = ''
          IF ( isRank  ) Rank  = 3

       CASE ( 'PV' )
          IF ( isDesc  ) Desc  = 'Ertel potential vorticity'
          IF ( isUnits ) Units = 'kg m2 kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'QI' )
          IF ( isDesc  ) Desc  = 'Ice mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'QL' )
          IF ( isDesc  ) Desc  = 'Water mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'RH' )
          IF ( isDesc  ) Desc  = 'Relative humidity'
          IF ( isUnits ) Units = '%'
          IF ( isRank  ) Rank  = 3

       CASE ( 'SPHU' )
          IF ( isDesc  ) Desc  = 'Specific humidity (w/r/t moist air)'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'SPHU_PREV' )
          IF ( isDesc  ) Desc  = 'Previous SPHU'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'T' )
          IF ( isDesc  ) Desc  = 'Temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3

       CASE ( 'TV' )
          IF ( isDesc  ) Desc  = 'Virtual temperature'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3

       CASE ( 'TAUCLI' )
          IF ( isDesc  ) Desc  = 'Optical depth of ice clouds'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'TAUCLW' )
          IF ( isDesc  ) Desc  = 'Optical depth of H2O clouds'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'U' )
          IF ( isDesc  ) Desc  = 'East-west component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'V' )
          IF ( isDesc  ) Desc  = 'North-south component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'UPDVVEL' )
          IF ( isDesc  ) Desc  = 'Updraft vertical velocity'
          IF ( isUnits ) Units = 'hPa s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'PFICU' )
          IF ( isDesc  ) Desc  = 'Downward flux of ice precipitation ' // &
                                 '(convective)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocEdge

       CASE ( 'PFILSAN' )
          IF ( isDesc  ) Desc  = 'Downwared flux of ice precipitation ' // &
                                 '(large-scale + anvil)'
          IF ( isRank  ) Rank  = 3
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isVLoc  ) VLoc  = VLocEdge

       CASE ( 'PFLCU' )
          IF ( isDesc  ) Desc  = 'Downward flux of liquid precipitation ' // &
                                 '(convective)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocEdge

       CASE ( 'PFLLSAN' )
          IF ( isDesc  ) Desc  = 'Downward flux of liquid precipitation ' // &
                                 '(large-scale + anvil)'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocEdge

       CASE ( 'REEVAPCN' )
          IF ( isDesc  ) Desc  = 'Evaporation of convective ' // &
                                 'precipitation (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg kg-1 s-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'REEVAPLS' )
          IF ( isDesc  ) Desc  = 'Evaporation of large-scale + anvil ' // &
                                 'precipitation (w/r/t dry air)'
          IF ( isUnits ) Units = 'kg '
          IF ( isRank  ) Rank  = 3

       CASE ( 'RH1' )
          IF ( isDesc  ) Desc  = 'Instantaneous relative humidity at time=T'
          IF ( isUnits ) Units = '%'
          IF ( isRank  ) Rank  = 3

       CASE ( 'RH2' )
          IF ( isDesc  ) Desc  = 'Instantaneous relative humidity at time=T+dt'
          IF ( isUnits ) Units = '%'
          IF ( isRank  ) Rank  = 3

       CASE ( 'SPHU1' )
          IF ( isDesc  ) Desc  = 'Instantaneous specific humidity at time=T'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocCenter

       CASE ( 'SPHU2' )
          IF ( isDesc  ) Desc  = 'Instantaneous specific humidity at time=T+dt'
          IF ( isUnits ) Units = 'g kg-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'TMPU1' )
          IF ( isDesc  ) Desc  = 'Instantaneous temperature at time=T'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3

       CASE ( 'TMPU2' )
          IF ( isDesc  ) Desc  = 'Instantaneous temperature at time T+dt'
          IF ( isUnits ) Units = 'K'
          IF ( isRank  ) Rank  = 3

       CASE ( 'IREG' )
          IF ( isDesc  ) Desc  = 'Number of Olson land types in each grid box'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2
          IF ( isType  ) Type  = KINDVAL_I4
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'ILAND' )
          IF ( isDesc  ) Desc  = 'Olson land type indices in each grid box'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isType  ) Type  = KINDVAL_I4
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'IUSE' )
          IF ( isDesc  ) Desc  = 'Fraction (per mil) occupied by each ' // &
                                 'Olson land type in the grid box'
          IF ( isUnits ) Units = 'o/oo'
          IF ( isRank  ) Rank  = 3
          IF ( isType  ) Type  = KINDVAL_I4
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'XLAI' )
          IF ( isDesc  ) Desc  = 'MODIS LAI for each Olson land type, ' // &
                                 'current month'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'MODISLAI' )
          IF ( isDesc  ) Desc  = 'Daily LAI computed from monthly ' // &
                                 'offline MODIS values'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 2

       CASE ( 'XCHLR' )
          IF ( isDesc  ) Desc  = 'MODIS chlorophyll-a per land type, ' // &
                                 'current month'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'MODISCHLR' )
          IF ( isDesc  ) Desc  = 'Daily chlorophyll-a computed ' // &
                                 'from offline MODIS monthly values'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 2

       CASE ( 'LANDTYPEFRAC' )
          IF ( isDesc  ) Desc  = 'Olson fraction per land type'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'XLAI_NATIVE' )
          IF ( isDesc  ) Desc  = 'Average LAI per Olson land type'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'XCHLR_NATIVE' )
          IF ( isDesc  ) Desc  = 'Average CHLR per Olson type'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'XLAI2' )
          IF ( isDesc  ) Desc  = 'MODIS chlorophyll-a per Olson land ' // &
                                 'type, next month'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE ( 'XCHLR2' )
          IF ( isDesc  ) Desc  = 'MODIS chlorophyll-a per Olson land ' // &
                                 'type, next month'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocUndefined

       CASE DEFAULT
          ErrMsg = 'No information available for field: ' // TRIM( Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

    END SELECT

    ! Set VLoc to undefined if variable is 2d
    IF ( isVLoc .AND. Rank == 2 ) THEN
       VLoc = VLocUndefined
    ENDIF

   END SUBROUTINE Get_MetField_Metadata
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_Rfp_2D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Rfp_2D( am_I_Root, varName,  Ptr2Data,         &
                                      State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: VarName         ! Variable name
    REAL(fp),          POINTER       :: Ptr2Data(:,:)   ! pointer to met
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),    INTENT(INOUT) :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC               ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: desc, units, ErrMsg,   ThisLoc
    INTEGER                :: rank, type,  vloc

    ! Initialize
    RC             = GC_SUCCESS
    units = ''
    desc = ''
    rank = -1
    type = -1
    vloc = -1
    ThisLoc        = ' -> at Register_MetField_Rfp_2D ' // &
                     '(in Headers/state_met_mod.F90)'

    CALL Get_MetField_Metadata( am_I_Root,   TRIM(VarName),  desc=desc, &
                                units=units, rank=rank,      type=type, &
                                vloc=vloc,   RC=RC                     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_MetField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
            TRIM(Varname)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL Registry_AddField( am_I_Root,         State_Met%Registry,    &
                            State_Met%State,   TRIM(VarName),         &
                            units=TRIM(units), Data2d=Ptr2Data,       &
                            Description=TRIM(desc),  RC=RC                 )
    CALL GC_CheckVar( 'State_Met%' // TRIM(VarName), 1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Registry_AddField'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_Rfp_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_Rfp_3D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Rfp_3D( am_I_Root, varName,  Ptr2Data,  &
                                      State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: VarName         ! Variable name
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),          POINTER       :: Ptr2Data(:,:,:) ! pointer to met
    TYPE(MetState),    INTENT(INOUT) :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC               ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: desc, units, ErrMsg,   ThisLoc
    INTEGER                :: rank, type,  vloc

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_MetField_Rfp_3D ' // &
                     '(in Headers/state_met_mod.F90)'

    CALL Get_MetField_Metadata( am_I_Root,   TRIM(VarName),  desc=desc, &
                                units=units, rank=rank,      type=type, &
                                vloc=vloc,   RC=RC                     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_MetField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
            TRIM(Varname)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL Registry_AddField( am_I_Root,              State_Met%Registry,    &
                            State_Met%State,        TRIM(VarName),         &
                            units=TRIM(units),      Data3d=Ptr2Data,       &
                            Description=TRIM(desc), RC=RC                 )
    CALL GC_CheckVar( 'State_Met%' // TRIM(VarName), 1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Registry_AddField'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_Rfp_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_MetField_Int_2D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_2D( am_I_Root, varName,  Ptr2Data, &
                                      State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: VarName         ! Variable name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           POINTER       :: Ptr2Data(:,:)   ! pointer to met
    TYPE(MetState),    INTENT(INOUT) :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC               ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: desc, units, ErrMsg,   ThisLoc
    INTEGER                :: rank, type,  vloc

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_MetField_Int_2D ' // &
                     '(in Headers/state_met_mod.F90)'

    CALL Get_MetField_Metadata( am_I_Root,   TRIM(VarName),  desc=desc, &
                                units=units, rank=rank,      type=type, &
                                vloc=vloc,   RC=RC                     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_MetField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
            TRIM(Varname)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL Registry_AddField( am_I_Root,              State_Met%Registry,    &
                            State_Met%State,        TRIM(VarName),         &
                            units=TRIM(units),      Data2d_I=Ptr2Data,     &
                            Description=TRIM(desc), RC=RC                 )
    CALL GC_CheckVar( 'State_Met%' // TRIM(VarName), 1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Registry_AddField'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
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
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_3D( am_I_Root, varName,  Ptr2Data,  &
                                       State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: VarName         ! Variable name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           POINTER       :: Ptr2Data(:,:,:) ! pointer to met
    TYPE(MetState),    INTENT(INOUT) :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC               ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: desc, units, ErrMsg,   ThisLoc
    INTEGER                :: rank, type,  vloc

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_MetField_Int_3D ' // &
                     '(in Headers/state_met_mod.F90)'

    CALL Get_MetField_Metadata( am_I_Root,   TRIM(VarName),  desc=desc, &
                                units=units, rank=rank,      type=type, &
                                vloc=vloc,   RC=RC                     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_MetField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
            TRIM(Varname)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL Registry_AddField( am_I_Root,              State_Met%Registry,    &
                            State_Met%State,        TRIM(VarName),         &
                            units=TRIM(units),      Data3d_I=Ptr2Data,     &
                            Description=TRIM(desc), RC=RC                 )
    CALL GC_CheckVar( 'State_Met%' // TRIM(VarName), 1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Registry_AddField'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_Int_3D
!EOC
END MODULE State_Met_Mod
