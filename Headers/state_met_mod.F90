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
     REAL(fp), POINTER :: ALBD          (:,:  ) ! Visible surface albedo [1]
     INTEGER,  POINTER :: ChemGridLev   (:,:  ) ! Chemistry grid level
     REAL(fp), POINTER :: CLDFRC        (:,:  ) ! Column cloud fraction [1]
     INTEGER,  POINTER :: CLDTOPS       (:,:  ) ! Max cloud top height [levels]
     REAL(fp), POINTER :: EFLUX         (:,:  ) ! Latent heat flux [W/m2]
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
     REAL(fp), POINTER :: LAI           (:,:  ) ! Leaf area index [m2/m2]  
                                                !  (online)
     REAL(fp), POINTER :: LWI           (:,:  ) ! Land/water indices [1]
     REAL(fp), POINTER :: PARDR         (:,:  ) ! Direct photsynthetically
                                                !  active radiation [W/m2]
     REAL(fp), POINTER :: PARDF         (:,:  ) ! Diffuse photsynthetically
                                                !  active radiation [W/m2]
     REAL(fp), POINTER :: PBLH          (:,:  ) ! PBL height [m]
     INTEGER,  POINTER :: PBL_TOP_L     (:,:  ) ! PBL top layer [1]
     REAL(fp), POINTER :: PHIS          (:,:  ) ! Surface geopotential height 
                                                !  [m2/s2]
     REAL(fp), POINTER :: PRECANV       (:,:  ) ! Anvil previp @ ground 
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PRECCON       (:,:  ) ! Conv  precip @ ground 
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PRECTOT       (:,:  ) ! Total precip @ ground 
                                                !  [kg/m2/s]
     REAL(fp), POINTER :: PRECLSC       (:,:  ) ! LS precip @ ground [kg/m2/s]
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
     REAL(fp), POINTER :: AREA_M2       (:,:,:) ! Grid box surface area [cm2]
     REAL(fp), POINTER :: CLDF          (:,:,:) ! 3-D cloud fraction [1]
     REAL(fp), POINTER :: CMFMC         (:,:,:) ! Cloud mass flux [kg/m2/s]
     REAL(fp), POINTER :: DQRCU         (:,:,:) ! Conv precip production rate 
                                                !  [kg/kg/s] (assume per 
                                                !  dry air)
     REAL(fp), POINTER :: DQRLSAN       (:,:,:) ! LS precip prod rate [kg/kg/s]
                                                !  (assume per dry air)
     REAL(fp), POINTER :: DTRAIN        (:,:,:) ! Detrainment flux [kg/m2/s]
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

     !----------------------------------------------------------------------
     ! Age of air for diagnosing transport
     !----------------------------------------------------------------------
     INTEGER,  POINTER :: AgeOfAir      (:,:,:) ! Age of air [s]
     
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
     REAL(fp), POINTER :: MODISCHLR     (:,:  ) ! Daily chlorophyll-a computed 
                                                !  from offline monthly MODIS
     REAL(fp), POINTER :: XLAI          (:,:,:) ! MODIS LAI per land type, 
                                                !  for this month
     REAL(fp), POINTER :: XCHLR         (:,:,:) ! MODIS CHLR per land type,
                                                !  for this month
     REAL(fp), POINTER :: LandTypeFrac  (:,:,:) ! Olson frac per type (I,J,type)
     REAL(fp), POINTER :: XLAI_NATIVE   (:,:,:) ! avg LAI per type (I,J,type)
     REAL(fp), POINTER :: XCHLR_NATIVE  (:,:,:) ! avg CHLR per type (I,J,type)

     REAL(fp), POINTER :: XLAI2         (:,:,:) ! MODIS LAI per land type, 
                                                !  for next month
     REAL(fp), POINTER :: XCHLR2        (:,:,:) ! MODIS CHLR per land type,
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
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!                              and remove obsolete met fields from State_Met
!  07 Sep 2017 - E. Lundgren - Add Register_MetField interface for init
!  13 Sep 2017 - M. Sulprizio- Remove DELP_PREV and SPHU_PREV; they're not used
!  14 Sep 2017 - M. Sulprizio- Comment out met fields that aren't actually used
!                              in GEOS-Chem (EVAP, GRN, PRECSNO, PV, RADLWG)
!  26 Sep 2017 - E. Lundgren - Remove Lookup_State_Met and Print_State_Met
!  07 Nov 2017 - R. Yantosca - Add tropht and troplev fields
!  08 Jan 2018 - R. Yantosca - Added logical query fields
!  31 Jan 2018 - E. Lundgren - Remove underscores from diagnostic names
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Register_MetField
     MODULE PROCEDURE Register_MetField_Rfp_2D
     MODULE PROCEDURE Register_MetField_Rfp_3D
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
! !IROUTINE: Init_State_Met
!
! !DESCRIPTION: Subroutine INIT\_STATE\_MET allocates all fields of 
!  the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Met( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR, NSURFTYPE
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
!  07 Sep 2017 - E. Lundgren - Abstract the metadata and method add to registry
!  16 Nov 2017 - E. Lundgren - Get grid params from CMN_Size_Mod not arguments
!  05 Nov 2018 - R. Yantosca - Now nullify all fields before allocating
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: LX, IM, JM, LM
    
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC    =  GC_SUCCESS
    ThisLoc = ' -> Init_State_Met (in Headers/state_met_mod.F90)'

    ! Shorten grid parameters for readability
    IM = IIPAR ! # latitudes
    JM = JJPAR ! # longitudes
    LM = LLPAR ! # levels

    !=======================================================================
    ! Nullify all fields for safety's sake before allocating them
    !=======================================================================
    State_Met%ALBD           => NULL()
    State_Met%ChemGridLev    => NULL()
    State_Met%CLDFRC         => NULL()
    State_Met%CLDTOPS        => NULL()
    State_Met%EFLUX          => NULL()
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
    State_Met%LAI            => NULL()
    State_Met%LWI            => NULL()
    State_Met%PARDR          => NULL()
    State_Met%PARDF          => NULL()
    State_Met%PBLH           => NULL()
    State_Met%PBL_TOP_L      => NULL()
    State_Met%PHIS           => NULL()
    State_Met%PRECANV        => NULL()
    State_Met%PRECCON        => NULL()
    State_Met%PRECLSC        => NULL()
    State_Met%PRECTOT        => NULL()
    State_Met%PS1_WET        => NULL()
    State_Met%PS1_DRY        => NULL()
    State_Met%PS2_WET        => NULL()
    State_Met%PS2_DRY        => NULL()
    State_Met%PSC2_WET       => NULL()
    State_Met%PSC2_DRY       => NULL()
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
    State_Met%SWGDN          => NULL()
    State_Met%TropLev        => NULL()
    State_Met%TropHt         => NULL()
    State_Met%TROPP          => NULL()
    State_Met%TS             => NULL()
    State_Met%TSKIN          => NULL()
    State_Met%TO3            => NULL()
    State_Met%U10M           => NULL()
    State_Met%USTAR          => NULL()
    State_Met%UVALBEDO       => NULL()
    State_Met%V10M           => NULL()
    State_Met%Z0             => NULL()
    State_Met%CNV_FRC        => NULL()
    State_Met%ILAND          => NULL()
    State_Met%IREG           => NULL()
    State_Met%IUSE           => NULL()
    State_Met%LANDTYPEFRAC   => NULL()
    State_Met%MODISLAI       => NULL()
    State_Met%MODISCHLR      => NULL()
    State_Met%AD             => NULL()
    State_Met%AIRDEN         => NULL()
    State_Met%MAIRDEN        => NULL()
    State_Met%AIRVOL         => NULL()
    State_Met%AREA_M2        => NULL()
    State_Met%BXHEIGHT       => NULL()
    State_Met%CLDF           => NULL()
    State_Met%CMFMC          => NULL()
    State_Met%DELP           => NULL()
    State_Met%DELP_DRY       => NULL()
    State_Met%DP_DRY_PREV    => NULL()
    State_Met%DQRCU          => NULL()
    State_Met%DQRLSAN        => NULL()
    State_Met%DTRAIN         => NULL()
    State_Met%OMEGA          => NULL()
    State_Met%OPTD           => NULL()
    State_Met%PEDGE          => NULL()
    State_Met%PEDGE_DRY      => NULL()
    State_Met%PFICU          => NULL()
    State_Met%PFILSAN        => NULL()
    State_Met%PFLCU          => NULL()
    State_Met%PFLLSAN        => NULL()
    State_Met%PMID           => NULL()
    State_Met%PMID_DRY       => NULL()
    State_Met%QI             => NULL()
    State_Met%QL             => NULL()
    State_Met%REEVAPCN       => NULL()
    State_Met%REEVAPLS       => NULL()
    State_Met%RH             => NULL()
    State_Met%SPHU           => NULL()
    State_Met%SPHU1          => NULL()
    State_Met%SPHU2          => NULL()
    State_Met%T              => NULL()
    State_Met%TMPU1          => NULL()
    State_Met%TMPU2          => NULL()
    State_Met%TV             => NULL()
    State_Met%TAUCLI         => NULL()
    State_Met%TAUCLW         => NULL()
    State_Met%U              => NULL()
    State_Met%UPDVVEL        => NULL()
    State_Met%V              => NULL()
    State_Met%XLAI           => NULL()
    State_Met%XCHLR          => NULL()
    State_Met%XLAI_NATIVE    => NULL()
    State_Met%XCHLR_NATIVE   => NULL()
    State_Met%InChemGrid     => NULL()
    State_Met%InPbl          => NULL()
    State_Met%InStratMeso    => NULL()
    State_Met%InStratosphere => NULL()
    State_Met%InTroposphere  => NULL()
    State_Met%IsLocalNoon    => NULL()
    State_Met%LocalSolarTime => NULL()
    State_Met%AgeOfAir       => NULL()

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
    IF ( RC /= GC_SUCCESS ) RETURN


    !-------------------------
    ! ChemGridLev [1]
    !-------------------------
    ALLOCATE( State_Met%ChemGridLev( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%ChemGridLev', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ChemGridLev = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CHEMGRIDLEV',  &
                            State_Met%ChemGridLev,     &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! CLDFRC [1]
    !-------------------------
    ALLOCATE( State_Met%CLDFRC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDFRC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CLDFRC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CLDFRC', State_Met%CLDFRC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! CLDTOPS [level]
    !-------------------------
    ALLOCATE( State_Met%CLDTOPS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDTOPS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CLDTOPS = 0
    CALL Register_MetField( am_I_Root, 'CLDTOPS', State_Met%CLDTOPS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! EFLUX [W m-2]
    !-------------------------
    ALLOCATE( State_Met%EFLUX( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%EFLUX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%EFLUX    = 0.0_fp
    CALL Register_MetField( am_I_Root, 'EFLUX', State_Met%EFLUX, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FRCLND [1]
    !-------------------------
    ALLOCATE( State_Met%FRCLND( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRCLND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRCLND = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRCLND', State_Met%FRCLND, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FRLAKE [1]
    !-------------------------
    ALLOCATE( State_Met%FRLAKE( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLAKE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLAKE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRLAKE', State_Met%FRLAKE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FRLAND [1]
    !-------------------------
    ALLOCATE( State_Met%FRLAND( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLAND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLAND = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'FRLAND', State_Met%FRLAND, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FRLANDIC [1]
    !-------------------------
    ALLOCATE( State_Met%FRLANDIC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRLANDIC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRLANDIC = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'FRLANDIC', State_Met%FRLANDIC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FROCEAN [1]
    !-------------------------
    ALLOCATE( State_Met%FROCEAN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FROCEAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FROCEAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FROCEAN', State_Met%FROCEAN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! GWETROOT [1]
    !-------------------------
    ALLOCATE( State_Met%GWETROOT( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%GWETROOT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%GWETROOT = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'GWETROOT', State_Met%GWETROOT, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! GWETTOP [1]
    !-------------------------
    ALLOCATE( State_Met%GWETTOP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%GWETTOP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%GWETTOP = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'GWETTOP', State_Met%GWETTOP, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! HFLUX [W m-2]
    !-------------------------
    ALLOCATE( State_Met%HFLUX( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%HFLUX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%HFLUX = 0.0_fp 
    CALL Register_MetField( am_I_Root, 'HFLUX', State_Met%HFLUX, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! LAI [1]
    !-------------------------  
    ALLOCATE( State_Met%LAI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%LAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LAI', State_Met%LAI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! LWI [1]
    !-------------------------
    ALLOCATE( State_Met%LWI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%LWI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LWI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LWI', State_Met%LWI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PARDR [W m-2]
    !-------------------------
    ALLOCATE( State_Met%PARDR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PARDR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PARDR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PARDR', State_Met%PARDR, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PARDF [W m-2]
    !-------------------------
    ALLOCATE( State_Met%PARDF( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PARDF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PARDF= 0.0_fp
    CALL Register_MetField( am_I_Root, 'PARDF', State_Met%PARDF, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PBLH [m]
    !-------------------------
    ALLOCATE( State_Met%PBLH( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PBLH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PBLH = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PBLH', State_Met%PBLH, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PBL_TOP_L [1]
    !-------------------------
    ALLOCATE( State_Met%PBL_TOP_L( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PBL_TOP_L', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PBL_TOP_L = 0
    CALL Register_MetField( am_I_Root, 'PBLTOPL', State_Met%PBL_TOP_L, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PHIS [m2 s-2]
    !-------------------------
    ALLOCATE( State_Met%PHIS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PHIS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PHIS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PHIS', State_Met%PHIS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PRECCON [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECCON( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECCON', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECCON = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECCON', State_Met%PRECCON, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PRECTOT [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECTOT( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECTOT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECTOT = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECTOT', State_Met%PRECTOT, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PS1_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS1_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS1_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS1_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS1WET', State_Met%PS1_WET, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PS2_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS2_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS2_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS2_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS2WET', State_Met%PS2_WET, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PSC2_WET [hPa]
    !-------------------------
    ALLOCATE( State_Met%PSC2_WET( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PSC2_WET', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PSC2_WET = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PSC2WET', State_Met%PSC2_WET, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PS1_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS1_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS1_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS1_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS1DRY', State_Met%PS1_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PS2_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PS2_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PS2_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PS2_DRY   = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PS2DRY', State_Met%PS2_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PSC2_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PSC2_DRY( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PSC2_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PSC2_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PSC2DRY', State_Met%PSC2_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SLP [hPa]
    !-------------------------
    ALLOCATE( State_Met%SLP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SLP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SLP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SLP', State_Met%SLP, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SNODP [m]
    !-------------------------
    ALLOCATE( State_Met%SNODP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SNODP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNODP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SNODP', State_Met%SNODP, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SNOMAS [kg m-2]
    !-------------------------
    ALLOCATE( State_Met%SNOMAS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SNOMAS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SNOMAS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SNOMAS', State_Met%SNOMAS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SUNCOS [1]
    !-------------------------
    ALLOCATE( State_Met%SUNCOS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SUNCOS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SUNCOS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SUNCOS', State_Met%SUNCOS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SUNCOSmid [1]
    !-------------------------
    ALLOCATE( State_Met%SUNCOSmid( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SUNCOSmid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SUNCOSmid = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SUNCOSmid', State_Met%SUNCOSmid, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SWGDN [W m-2]
    !-------------------------
    ALLOCATE( State_Met%SWGDN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SWGDN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SWGDN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SWGDN', State_Met%SWGDN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TO3 [dobsons]
    !-------------------------
    ALLOCATE( State_Met%TO3( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TO3', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TO3 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TO3', State_Met%TO3, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TropLev [1]
    !-------------------------
    ALLOCATE( State_Met%TropLev( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TropLev', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TropLev = 0
    CALL Register_MetField( am_I_Root, 'TROPLEV', State_Met%TropLev, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TropHt [hPa]
    !-------------------------
    ALLOCATE( State_Met%TropHt( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TropHt', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TropHt = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TROPHT', State_Met%TropHt, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TROPP [hPa]
    !-------------------------
    ALLOCATE( State_Met%TROPP( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TROPP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TROPP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TROPP', State_Met%TROPP, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TS [K]
    !-------------------------
    ALLOCATE( State_Met%TS( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TS', State_Met%TS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TSKIN [1]
    !-------------------------
    ALLOCATE( State_Met%TSKIN( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TSKIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TSKIN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TSKIN', State_Met%TSKIN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! U10M [m s-1]
    !-------------------------
    ALLOCATE( State_Met%U10M( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%U10M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%U10M = 0.0_fp
    CALL Register_MetField( am_I_Root, 'U10M', State_Met%U10M, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! USTAR [m -s]
    !-------------------------
    ALLOCATE( State_Met%USTAR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%USTAR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%USTAR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'USTAR', State_Met%USTAR, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! UVALBEDO [1]
    !-------------------------
    ALLOCATE( State_Met%UVALBEDO( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%UVALBEDO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UVALBEDO = 0.0_fp
    CALL Register_MetField( am_I_Root, 'UVALBEDO', State_Met%UVALBEDO, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! V10M [m s-1]
    !-------------------------
    ALLOCATE( State_Met%V10M( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%V10M', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%V10M = 0.0_fp
    CALL Register_MetField( am_I_Root, 'V10M', State_Met%V10M, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! Z0 [m]
    !-------------------------
    ALLOCATE( State_Met%Z0( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Z0', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%Z0 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'Z0', State_Met%Z0, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Convective fractions are not yet a standard GEOS-FP
    ! field. Only available to online model (ckeller, 3/4/16) 
#if defined( ESMF_ ) || defined( MODEL_ )
    !-------------------------
    ! CNV_FRC [1]
    !-------------------------
    ALLOCATE( State_Met%CNV_FRC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CNV_FRC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%CNV_FRC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CNVFRC', State_Met%CNV_FRC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    !-------------------------
    ! FRESEAICE [1]
    !-------------------------
    ALLOCATE( State_Met%FRSEAICE( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRSEAICE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRSEAICE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRSEAICE', State_Met%FRSEAICE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! FRSNO [1]
    !-------------------------
    ALLOCATE( State_Met%FRSNO( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%FRSNO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%FRSNO = 0.0_fp
    CALL Register_MetField( am_I_Root, 'FRSNO', State_Met%FRSNO, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PRECANV [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECANV( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECANV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECANV = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECANV', State_Met%PRECANV, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PRECLSC [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PRECLSC( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PRECLSC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PRECLSC  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PRECLSC', State_Met%PRECLSC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE00 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE00( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE00', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE00 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE00', State_Met%SEAICE00, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE10 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE10( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE10', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE10 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE10', State_Met%SEAICE10, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE20 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE20( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE20', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE20 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE20', State_Met%SEAICE20, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE30 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE30( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE30', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE30 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE30', State_Met%SEAICE30, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE40 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE40( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE40', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE40 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE40', State_Met%SEAICE40, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE50 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE50( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE50', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE50 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE50', State_Met%SEAICE50, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE60 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE60( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE60', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE60 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE60', State_Met%SEAICE60, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE70 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE70( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE70', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE70 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE70', State_Met%SEAICE70, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE80 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE80( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE80', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE80 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE80', State_Met%SEAICE80, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SEAICE90 [1]
    !-------------------------
    ALLOCATE( State_Met%SEAICE90( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SEAICE90', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SEAICE90 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SEAICE90', State_Met%SEAICE90, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

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
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! AIRDEN [kg m-3]
    !-------------------------
    ALLOCATE( State_Met%AIRDEN( IM, JM, LM ), STAT=RC )  
    CALL GC_CheckVar( 'State_Met%AIRDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AIRDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRDEN', State_Met%AIRDEN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! MAIRDEN [kg m-3]
    !-------------------------
    ALLOCATE( State_Met%MAIRDEN( IM, JM, LM   ), STAT=RC )  
    CALL GC_CheckVar( 'State_Met%MAIRDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%MAIRDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MAIRDEN', State_Met%MAIRDEN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! AIRNUMDEN [1]
    !-------------------------        
    ALLOCATE( State_Met%AIRNUMDEN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AIRNUMDEN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%AIRNUMDEN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRNUMDEN', State_Met%AIRNUMDEN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! AIRVOL [m3]
    !-------------------------
    ALLOCATE( State_Met%AIRVOL( IM, JM, LM  ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AIRVOL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AIRVOL = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AIRVOL', State_Met%AIRVOL, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! AREA_M2 [m2]
    !-------------------------
    ALLOCATE( State_Met%AREA_M2( IM, JM, 1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AREA_M2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%AREA_M2  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AREAM2', State_Met%AREA_M2, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! AVGW [v/v]
    !-------------------------
    ALLOCATE( State_Met%AVGW( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AVGW', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%AVGW = 0.0_fp
    CALL Register_MetField( am_I_Root, 'AVGW', State_Met%AVGW, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! BXHEIGHT [m]
    !-------------------------
    ALLOCATE( State_Met%BXHEIGHT( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%BXHEIGHT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%BXHEIGHT = 0.0_fp
    CALL Register_MetField( am_I_Root, 'BXHEIGHT', State_Met%BXHEIGHT, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! CLDF [1]
    !-------------------------                 
    ALLOCATE( State_Met%CLDF( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CLDF', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%CLDF = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CLDF', State_Met%CLDF, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! CMFMC [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%CMFMC( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%CMFMC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%CMFMC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'CMFMC', State_Met%CMFMC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DELP [hPa]
    !-------------------------           
    ALLOCATE( State_Met%DELP( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DELP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DELP = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DELP', State_Met%DELP, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DELP_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%DELP_DRY( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DELP_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DELP_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DELPDRY', State_Met%DELP_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DP_DRY_PREV [hPa]
    !-------------------------
    ALLOCATE( State_Met%DP_DRY_PREV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DP_DRY_PREV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DP_DRY_PREV= 0.0_fp
    CALL Register_MetField( am_I_Root, 'DPDRYPREV', State_Met%DP_DRY_PREV, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DQRCU [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQRCU( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQRCU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DQRCU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQRCU', State_Met%DQRCU, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DQRLSAN [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%DQRLSAN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DQRLSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%DQRLSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DQRLSAN', State_Met%DQRLSAN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! DTRAIN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%DTRAIN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%DTRAIN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%DTRAIN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'DTRAIN', State_Met%DTRAIN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! OMEGA [Pa s-1]
    !-------------------------
    ALLOCATE( State_Met%OMEGA( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%OMEGA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%OMEGA = 0.0_fp
    CALL Register_MetField( am_I_Root, 'OMEGA', State_Met%OMEGA, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! OPTD [1]
    !-------------------------
    ALLOCATE( State_Met%OPTD( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%OPTD', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%OPTD = 0.0_fp
    CALL Register_MetField( am_I_Root, 'OPTD', State_Met%OPTD, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
                            
    !-------------------------
    ! PEDGE [hPa]
    !-------------------------
    ALLOCATE( State_Met%PEDGE( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PEDGE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PEDGE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PEDGE', State_Met%PEDGE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PEDGE_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PEDGE_DRY ( IM, JM, LM+1 ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PEDGE_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PEDGE_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PEDGEDRY', State_Met%PEDGE_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PMID [1]
    !-------------------------
    ALLOCATE( State_Met%PMID( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PMID', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PMID = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PMID', State_Met%PMID, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PMID_DRY [hPa]
    !-------------------------
    ALLOCATE( State_Met%PMID_DRY( IM, JM, LM   ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PMID_DRY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%PMID_DRY = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PMIDDRY', State_Met%PMID_DRY, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! QI [kg kg-1]
    !-------------------------
    ALLOCATE( State_Met%QI( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%QI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%QI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'QI', State_Met%QI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! QL [kg kg-1]
    !-------------------------
    ALLOCATE( State_Met%QL( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%QL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%QL = 0.0_fp
    CALL Register_MetField( am_I_Root, 'QL', State_Met%QL, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! RH [%]
    !-------------------------
    ALLOCATE( State_Met%RH ( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%RH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN                               
    State_Met%RH = 0.0_fp
    CALL Register_MetField( am_I_Root, 'RH', State_Met%RH, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SPHU [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%SPHU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU', State_Met%SPHU, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! T [K]
    !-------------------------                                               
    ALLOCATE( State_Met%T( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%T', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%T = 0.0_fp
    CALL Register_MetField( am_I_Root, 'T', State_Met%T, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! THETA [K]
    !-------------------------
    ALLOCATE( State_Met%THETA( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%THETA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%THETA = 0.0_fp
    CALL Register_MetField( am_I_Root, 'THETA', State_Met%THETA, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TV [K]
    !-------------------------
    ALLOCATE( State_Met%TV( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TV', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%TV = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TV', State_Met%TV, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TAUCLI [1]
    !-------------------------
    ALLOCATE( State_Met%TAUCLI( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TAUCLI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN           
    State_Met%TAUCLI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TAUCLI', State_Met%TAUCLI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TAUCLW [1]
    !-------------------------
    ALLOCATE( State_Met%TAUCLW( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TAUCLW', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TAUCLW = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TAUCLW', State_Met%TAUCLW, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! U [m s-1]
    !-------------------------
    ALLOCATE( State_Met%U( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%U', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%U = 0.0_fp
    CALL Register_MetField( am_I_Root, 'U', State_Met%U, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! V [m s-1]
    !-------------------------
    ALLOCATE( State_Met%V( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%V', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%V = 0.0_fp
    CALL Register_MetField( am_I_Root, 'V', State_Met%V, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Updraft vertical velocity is not yet a standard GEOS-FP
    ! field. Only available to online model (ckeller, 3/4/16) 
#if defined( ESMF_ ) || defined( MODEL_ )
    !-------------------------
    ! UPDVVEL [hPa s-1]
    !-------------------------
    ALLOCATE( State_Met%UPDVVEL( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%UPDVVEL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%UPDVVEL  = -999.0_fp
    CALL Register_MetField( am_I_Root, 'UPDVVEL', State_Met%UPDVVEL, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
#endif

    ! Pick the proper vertical dimension
    LX = LM + 1           ! For fields that are on level edges
    
    !-------------------------
    ! PFICU [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFICU( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFICU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFICU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFICU', State_Met%PFICU, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PFILSAN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFILSAN( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFILSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFILSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFILSAN', State_Met%PFILSAN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PFLCU [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFLCU( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFLCU', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFLCU = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFLCU', State_Met%PFLCU, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! PFLLSAN [kg m-2 s-1]
    !-------------------------
    ALLOCATE( State_Met%PFLLSAN( IM, JM, LX ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%PFLLSAN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%PFLLSAN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'PFLLSAN', State_Met%PFLLSAN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! REEVAPCN [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%REEVAPCN( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%REEVAPCN', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%REEVAPCN = 0.0_fp
    CALL Register_MetField( am_I_Root, 'REEVAPCN', State_Met%REEVAPCN, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! REEVAPLS [kg kg-1 s-1]
    !-------------------------
    ALLOCATE( State_Met%REEVAPLS( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%REEVAPLS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%REEVAPLS = 0.0_fp
    CALL Register_MetField( am_I_Root, 'REEVAPLS', State_Met%REEVAPLS, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SPHU1 [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU1( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SPHU1 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU1', State_Met%SPHU1, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! SPHU2 [g kg-1]
    !-------------------------
    ALLOCATE( State_Met%SPHU2( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%SPHU2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%SPHU2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'SPHU2', State_Met%SPHU2, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TMPU1 [K]
    !-------------------------
    ALLOCATE( State_Met%TMPU1( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TMPU1', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TMPU1 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TMPU1', State_Met%TMPU1, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! TMPU2 [K]
    !-------------------------
    ALLOCATE( State_Met%TMPU2( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%TMPU2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%TMPU2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'TMPU2', State_Met%TMPU2, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! Age of Air [s]
    !-------------------------
    ALLOCATE( State_Met%AgeOfAir( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%AgeOfAir', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%AgeOfAir = 0
    CALL Register_MetField( am_I_Root, 'AgeOfAir', State_Met%AgeOfAir, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    
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
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! ILAND [1]
    !-------------------------
    ALLOCATE( State_Met%ILAND( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%ILAND', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%ILAND = 0
    CALL Register_MetField( am_I_Root, 'ILAND', State_Met%ILAND, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! IUSE [1]
    !-------------------------
    ALLOCATE( State_Met%IUSE( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%IUSE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%IUSE = 0
    CALL Register_MetField( am_I_Root, 'IUSE', State_Met%IUSE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XLAI [1]
    !-------------------------
    ALLOCATE( State_Met%XLAI( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XLAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAI', State_Met%XLAI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! MODISLAI [1]
    !-------------------------
    ALLOCATE( State_Met%MODISLAI( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%MODISLAI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%MODISLAI = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MODISLAI', State_Met%MODISLAI, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XCHLR [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%XCHLR( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XCHLR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLR', State_Met%XCHLR, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XLAI2 [1]
    !-------------------------
    ALLOCATE( State_Met%XLAI2( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XLAI2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAI2', State_Met%XLAI2, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XCHLR2 [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%XCHLR2( IM, JM, NSURFTYPE ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%XCHLR2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR2 = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLR2', State_Met%XCHLR2, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! MODISCHLR [mg m-3]
    !-------------------------
    ALLOCATE( State_Met%MODISCHLR( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%MODISCHLR', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%MODISCHLR = 0.0_fp
    CALL Register_MetField( am_I_Root, 'MODISCHLR', State_Met%MODISCHLR, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! LANDTYPEFRAC [1]
    !-------------------------    
    ALLOCATE( State_Met%LANDTYPEFRAC( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%LANDTYPEFRAC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LANDTYPEFRAC = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LANDTYPEFRAC', State_Met%LANDTYPEFRAC, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XLAI_NATIVE [1]
    !-------------------------
    ALLOCATE( State_Met%XLAI_NATIVE( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XLAI_NATIVE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XLAI_NATIVE  = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XLAINATIVE', State_Met%XLAI_NATIVE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !-------------------------
    ! XCHLR_NATIVE [1]
    !-------------------------
    ALLOCATE( State_Met%XCHLR_NATIVE( IM, JM, NSURFTYPE ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%XCHLR_NATIVE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%XCHLR_NATIVE = 0.0_fp
    CALL Register_MetField( am_I_Root, 'XCHLRNATIVE', State_Met%XCHLR_NATIVE, &
                            State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN


    !=======================================================================
    ! Allocate fields for querying which vertical regime a grid box is in
    ! or if a grid box is near local solar noontime.
    ! 
    ! %%%%% NOTE: Do not register these query fields %%%%%
    !=======================================================================

    !-------------------------
    ! InChemGrid
    !-------------------------
    ALLOCATE( State_Met%InChemGrid( IM, JM, LM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%IsChemGrid', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%InChemGrid = .FALSE.

    !-------------------------
    ! InPBL
    !-------------------------
    ALLOCATE( State_Met%InPbl( IM, JM, LM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%InPbl', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%InPbl = .FALSE.

    !-------------------------
    ! InStratosphere
    !-------------------------
    ALLOCATE( State_Met%InStratosphere( IM, JM, LM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%InStratosphere', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%InStratosphere = .FALSE.

    !-------------------------
    ! InStratMeso
    !-------------------------
    ALLOCATE( State_Met%InStratMeso( IM, JM, LM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%InStratMeso', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%InStratMeso = .FALSE.

    !-------------------------
    ! InTroposphere
    !-------------------------
    ALLOCATE( State_Met%InTroposphere( IM, JM, LM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%InTropoSphere', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%InTroposphere = .FALSE.

    !-------------------------
    ! IsLocalNoon
    !-------------------------
    ALLOCATE( State_Met%IsLocalNoon( IM, JM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%IsLocalNoon', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%IsLocalNoon = .FALSE.

    !-------------------------
    ! LocalSolarTime
    !-------------------------
    ALLOCATE( State_Met%LocalSolarTime( IM, JM ), STAT=RC )        
    CALL GC_CheckVar( 'State_Met%LocalSolarTime', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Met%LocalSolarTime = 0.0_fp
    CALL Register_MetField( am_I_Root, 'LOCALSOLARTIME',                     &
                            State_Met%LocalSolarTime,                        &
                            State_Met, RC                             )
    IF ( RC /= GC_SUCCESS ) RETURN

    !=======================================================================
    ! Print information about the registered fields (short format)
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the State_Met object:')
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( am_I_Root   = am_I_Root,                            &
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
  SUBROUTINE Cleanup_State_Met( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  05 Nov 2018 - R. Yantosca - Now deallocate AND nullify all pointer fields
!                              (except 3D fields for GCHP/GEOS/WRF interfaces)
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

    IF ( ASSOCIATED( State_Met%EFLUX ) ) THEN
       DEALLOCATE( State_Met%EFLUX, STAT=RC )
       CALL GC_CheckVar( 'State_Met%EFLUX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%EFLUX => NULL()
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

    IF ( ASSOCIATED( State_Met%LAI ) ) THEN
       DEALLOCATE( State_Met%LAI, STAT=RC )
       CALL GC_CheckVar( 'State_Met%LAI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%LAI => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Met%LWI ) ) THEN
       DEALLOCATE( State_Met%LWI, STAT=RC )
       CALL GC_CheckVar( 'State_Met%LWI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%LWI => NULL()
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

    IF ( ASSOCIATED( State_Met%PBL_TOP_L ) ) THEN
       DEALLOCATE( State_Met%PBL_TOP_L, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%PBL_TOP_L', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%PBL_TOP_L => NULL()
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

    IF ( ASSOCIATED( State_Met%MODISCHLR ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%MODISCHLR => NULL()
#else
       DEALLOCATE( State_Met%MODISCHLR, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%MODISCHLR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%MODISCHLR => NULL()
#endif
    ENDIF

    !---------------------------
    ! 3-D fields
    !---------------------------
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

    IF ( ASSOCIATED( State_Met%AREA_M2 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%AREA_M2 => NULL()
#else
       DEALLOCATE( State_Met%AREA_M2, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%AREA_M2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AREA_M2 => NULL()
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

    IF ( ASSOCIATED( State_Met%XCHLR ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XCHLR => NULL()
#else
       DEALLOCATE( State_Met%XCHLR, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XCHLR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XCHLR => NULL()
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

    IF ( ASSOCIATED( State_Met%XCHLR2 ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XCHLR2 => NULL()
#else
       DEALLOCATE( State_Met%XCHLR2, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XCHLR2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XCHLR2 => NULL()
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

    IF ( ASSOCIATED( State_Met%XCHLR_NATIVE ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%XCHLR_NATIVE => NULL()
#else
       DEALLOCATE( State_Met%XCHLR_NATIVE, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%XCHLR_NATIVE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%XCHLR_NATIVE => NULL()
#endif
    ENDIF

    IF ( ASSOCIATED( State_Met%AgeOfAir ) ) THEN
#if defined( ESMF_ ) || defined( MODEL_WRF )
       State_Met%AgeOfAir => NULL()
#else
       DEALLOCATE( State_Met%AgeOfAir, STAT=RC  )
       CALL GC_CheckVar( 'State_Met%AgeOfAir', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Met%AgeOfAir => NULL()
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

    !-----------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
    !-----------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Met%xxx ) ) THEN
    !   DEALLOCATE( State_Met%xxx, STAT=RC )
    !   CALL GC_CheckVar( 'State_Met%xxx', 2, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Met%xxx => NULL()
    !ENDIF

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, State_Met%Registry, RC )
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
  SUBROUTINE Get_Metadata_State_Met( am_I_Root, metadataID, Found,  &
                                     RC,        Desc,       Units,  &
                                     Rank,      Type,       VLoc )
!
! !USES:
!
    USE Charpak_Mod,         ONLY: To_UpperCase
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
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
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  28 Aug 2017 - E. Lundgren - Initial version
!  01 Nov 2017 - R. Yantosca - Now get To_UpperCase from charpak_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Name_AllCaps
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Met (in Headers/state_met_mod.F90)'
    Found = .TRUE.

    ! Optional arguments present?
    isDesc  = PRESENT( Desc  )
    isUnits = PRESENT( Units )
    isRank  = PRESENT( Rank  )
    isType  = PRESENT( Type  )
    isVLoc  = PRESENT( VLoc  )

    ! Set defaults for optional arguments. Assume type and vertical 
    ! location are real (flexible precision) and center unless specified 
    ! otherwise
    IF ( isUnits ) Units = ''
    IF ( isDesc  ) Desc  = ''              
    IF ( isRank  ) Rank  = -1              ! initialize as bad value 
    IF ( isType  ) Type  = KINDVAL_FP      ! Assume real with flex precision
    IF ( isVLoc  ) VLoc  = VLocationNone   ! Assume no vertical location

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM( Name_AllCaps) )

       CASE ( 'ALBD' )
          IF ( isDesc  ) Desc  = 'Visible surface albedo'
          IF ( isUnits ) Units = '1'
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

       CASE ( 'EFLUX' )
          IF ( isDesc  ) Desc  = 'Latent heat flux'
          IF ( isUnits ) Units = 'W m-2'
          IF ( isRank  ) Rank  = 2

!------------------------------------------------------------------------------
! Comment out for now. State_Met%EVAP is not used in the code. (mps, 9/14/17)
!       CASE ( 'EVAP' )
!          IF ( isDesc  ) Desc  = 'Surface evaporation'
!          IF ( isUnits ) Units = 'kg m-2 s-1'
!          IF ( isRank  ) Rank  = 2
!------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
! Comment out for now. State_Met%EVAP is not used in the code. (mps, 9/14/17)
!       CASE ( 'GRN' )
!          IF ( isDesc  ) Desc  = 'Greenness fraction'
!          IF ( isUnits ) Units = '1'
!          IF ( isRank  ) Rank  = 2
!------------------------------------------------------------------------------

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

       CASE ( 'PBLTOPL' )
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

!------------------------------------------------------------------------------
! Comment out for now. State_Met%EVAP is not used in the code. (mps, 9/14/17)
!       CASE ( 'PRECSNO' )
!          IF ( isDesc  ) Desc  = 'Snow precipitation'
!          IF ( isUnits ) Units = 'kg m-2 s-1'
!          IF ( isRank  ) Rank  = 2
!------------------------------------------------------------------------------

       CASE ( 'PRECTOT' )
          IF ( isDesc  ) Desc  = 'Total precipitation at the ground'
          IF ( isUnits ) Units = 'kg m-2 s-1'
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

!------------------------------------------------------------------------------
! Comment out for now. State_Met%EVAP is not used in the code. (mps, 9/14/17)
!       CASE ( 'RADLWG' )
!          IF ( isDesc  ) Desc  = 'Net longwave radiation at ground'
!          IF ( isUnits ) Units = 'W m-2'
!          IF ( isRank  ) Rank  = 2
!------------------------------------------------------------------------------
 
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

#if defined( ESMF_ ) || defined( MODEL_ )
       CASE ( 'CNVFRC' )
          IF ( isDesc  ) Desc  = 'Convective fraction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

#endif
       CASE ( 'FRSEAICE' )
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

       CASE ( 'SEAICE00' )
          IF ( isDesc  ) Desc  = 'Sea ice coverage 00-10%'
          IF ( isUnits ) Units = '1'
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

       CASE ( 'LOCALSOLARTIME' )
          IF ( isDesc  ) Desc  = 'Local solar time'
          IF ( isUnits ) Units = 'hours'
          IF ( isRank  ) Rank  = 2

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
          IF ( isUnits ) Units = 'm-3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AIRVOL' )
          IF ( isDesc  ) Desc  = 'Volume of dry air in grid box'
          IF ( isUnits ) Units = 'm3'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AREAM2' )
          IF ( isDesc  ) Desc  = 'Surface area of grid box'
          IF ( isUnits ) Units = 'm2'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'AVGW' )
          IF ( isDesc  ) Desc  = 'Water vapor mixing ratio (w/r/t dry air)'
          IF ( isUnits ) Units = 'vol vol-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'BXHEIGHT' )
          IF ( isDesc  ) Desc  = 'Grid box height (w/r/t dry air)'
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

!------------------------------------------------------------------------------
! Comment out for now. State_Met%EVAP is not used in the code. (mps, 9/14/17)
!       CASE ( 'PV' )
!          IF ( isDesc  ) Desc  = 'Ertel potential vorticity'
!          IF ( isUnits ) Units = 'kg m2 kg-1 s-1'
!          IF ( isRank  ) Rank  = 3
!          IF ( isVLoc  ) VLoc  = VLocationCenter
!------------------------------------------------------------------------------

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

       CASE ( 'U' )
          IF ( isDesc  ) Desc  = 'East-west component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'V' )
          IF ( isDesc  ) Desc  = 'North-south component of wind'
          IF ( isUnits ) Units = 'm s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

       CASE ( 'UPDVVEL' )
          IF ( isDesc  ) Desc  = 'Updraft vertical velocity'
          IF ( isUnits ) Units = 'hPa s-1'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter

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

       CASE ( 'IREG' )
          IF ( isDesc  ) Desc  = 'Number of Olson land types in each grid box'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2
          IF ( isType  ) Type  = KINDVAL_I4

       CASE ( 'ILAND' )
          IF ( isDesc  ) Desc  = 'Olson land type indices in each grid box'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3
          IF ( isType  ) Type  = KINDVAL_I4

       CASE ( 'IUSE' )
          IF ( isDesc  ) Desc  = 'Fraction (per mil) occupied by each ' // &
                                 'Olson land type in the grid box'
          IF ( isUnits ) Units = 'o/oo'
          IF ( isRank  ) Rank  = 3
          IF ( isType  ) Type  = KINDVAL_I4

       CASE ( 'XLAI' )
          IF ( isDesc  ) Desc  = 'MODIS LAI for each Olson land type, ' // &
                                 'current month'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3

       CASE ( 'XLAI2' )
          IF ( isDesc  ) Desc  = 'MODIS LAI for each Olson land type, ' // &
                                 'next month'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3

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

       CASE ( 'XCHLR2' )
          IF ( isDesc  ) Desc  = 'MODIS chlorophyll-a per land type, ' // &
                                 'next month'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'MODISCHLR' )
          IF ( isDesc  ) Desc  = 'Daily chlorophyll-a computed ' // &
                                 'from offline MODIS monthly values'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 2

       CASE ( 'LANDTYPEFRAC' )
          IF ( isDesc  ) Desc  = 'Olson fraction per land type'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'XLAINATIVE' )
          IF ( isDesc  ) Desc  = 'Average LAI per Olson land type'
          IF ( isUnits ) Units = 'm2 m-2'
          IF ( isRank  ) Rank  = 3

       CASE ( 'XCHLRNATIVE' )
          IF ( isDesc  ) Desc  = 'Average CHLR per Olson type'
          IF ( isUnits ) Units = 'mg m-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AGEOFAIR' )
          IF ( isDesc  ) Desc  = 'Age of air'
          IF ( isUnits ) Units = 's'
          IF ( isRank  ) Rank  = 3
          IF ( isVLoc  ) VLoc  = VLocationCenter
          
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
! !IROUTINE: Register_MetField_Rfp_2D
!
! !DESCRIPTION: Registers a 2-D State\_Met field (flexible precision).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Rfp_2D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(fp),          POINTER       :: Ptr2Data(:,:)   ! pointer to met
    TYPE(MetState),    INTENT(IN)    :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
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
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: desc,  units, ErrMsg_reg, ThisLoc
    INTEGER                :: rank,  type,  vloc
    LOGICAL                :: found

    !---------------------   
    ! Initialize
    !---------------------   
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_MetField_Rfp_2D (in Headers/state_met_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'

    !---------------------
    ! Get metadata
    !---------------------
    CALL Get_Metadata_State_Met( am_I_Root, metadataID,  found, RC,          &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc                       )
    
    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------
    ! Check dimensions
    !---------------------
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------
    ! Add to registry
    !---------------------
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = State_Met%Registry,                &
                            State       = State_Met%State,                   &
                            Variable    = TRIM( MetadataID ),                &
                            Units       = TRIM( Units      ),                &
                            Description = TRIM( Desc       ),                &
                            Data2d      = Ptr2Data,                          &
                            RC          = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Registry_AddField"!'
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
! !DESCRIPTION: Registers a 3-D State\_Met field (flexible precision).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Rfp_3D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(fp),          POINTER       :: Ptr2Data(:,:,:) ! pointer to met
    TYPE(MetState),    INTENT(IN)    :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
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
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: desc,  units,  ErrMsg_reg, ThisLoc
    INTEGER                :: rank,  type,   vloc
    LOGICAL                :: found, onEdges

    !---------------------   
    ! Initialize
    !---------------------   
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_MetField_Rfp_3D (in Headers/state_met_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'

    !---------------------   
    ! Get metadata
    !---------------------
    CALL Get_Metadata_State_Met( am_I_Root, metadataID,  found, RC,          &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------   
    ! Check dimensions
    !---------------------  
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == vLocationEdge )

    !---------------------   
    ! Add to registry
    !---------------------   
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = State_Met%Registry,               &
                            State        = State_Met%State,                  &
                            Variable     = TRIM( MetadataID ),               &
                            Units        = TRIM( Units      ),               &
                            Description  = TRIM( Desc       ),               &
                            OnLevelEdges = onEdges,                          &
                            Data3d       = Ptr2Data,                         &
                            RC           = RC                               )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Registry_AddField"!'
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
! !DESCRIPTION: Registers a 2-D State\_Met field (integer precision).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_2D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
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
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: desc, units, ErrMsg_reg, ThisLoc
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found

    !---------------------   
    ! Initialize
    !---------------------   
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_MetField_Int_2D (in Headers/state_met_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'

    !---------------------   
    ! Get metadata
    !---------------------   
    CALL Get_Metadata_State_Met( am_I_Root, metadataID,  found, RC,          &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc                       )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------
    ! Check dimensions
    !---------------------
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------   
    ! Add to registry
    !---------------------   
    CALL Registry_AddField( am_I_Root   = am_I_Root,                         &
                            Registry    = State_Met%Registry,                &
                            State       = State_Met%State,                   &
                            Variable    = TRIM( metadataID ),                &
                            units       = TRIM( units      ),                &
                            Description = TRIM( desc       ),                &
                            Data2d_I    = Ptr2Data,                          &
                            RC          = RC                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Registry_AddField"!'
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
! !DESCRIPTION: Registers a 3-D State\_Met field (integer precision).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_MetField_Int_3D( am_I_Root, metadataID, Ptr2Data,      &
                                       State_Met, RC                        )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    INTEGER,           POINTER       :: Ptr2Data(:,:,:) ! pointer to met
    TYPE(MetState),    INTENT(IN)    :: State_Met       ! Obj for met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
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
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: desc,   units,  ErrMsg_reg, ThisLoc
    INTEGER                :: rank,   type,   vloc
    LOGICAL                :: found,  onEdges

    !---------------------   
    ! Initialize
    !---------------------   
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_MetField_Int_3D (in Headers/state_met_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Met%'

    !---------------------   
    ! Get metadata
    !---------------------   
    CALL Get_Metadata_State_Met( am_I_Root, metadataID,  found, RC,          &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc                       )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Met"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !---------------------
    ! Check dimensions
    !---------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data and metadata rank do not match for ' // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == vLocationEdge )

    !---------------------   
    ! Add to registry
    !---------------------   
    CALL Registry_AddField( am_I_Root    = am_I_Root,                        &
                            Registry     = State_Met%Registry,               &
                            State        = State_Met%State,                  &
                            Variable     = TRIM( metadataID ),               &
                            units        = TRIM( units      ),               &
                            Description  = TRIM( desc       ),               &
                            OnLevelEdges = onEdges,                          &
                            Data3d_I     = Ptr2Data,                         &
                            RC           = RC                               )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Registry_AddField"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_MetField_Int_3D
!EOC
END MODULE State_Met_Mod
