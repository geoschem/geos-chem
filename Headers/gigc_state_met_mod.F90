!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_state_met_mod
!
! !DESCRIPTION: Module GIGC\_STATE\_MET\_MOD contains the derived type
!  used to define the Meteorology State object for the Grid-Independent 
!  GEOS-Chem implementation (abbreviated "GIGC").
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
MODULE GIGC_State_Met_Mod
!
! USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_GIGC_State_Met
  PUBLIC :: Cleanup_GIGC_State_Met
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
     REAL(fp), POINTER :: ALBD      (:,:  ) ! Visible surface albedo [1
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
     REAL(fp), POINTER :: PS1       (:,:  ) ! Sfc press at timestep start[hPa]
                                            ! (wet air pressure)
     REAL(fp), POINTER :: PS2       (:,:  ) ! Sfc press at timestep end [hPa]
                                            ! (wet air pressure)
     REAL(fp), POINTER :: PSC2      (:,:  ) ! Interpolated sfc pressure [hPa]
                                            ! (wet air pressure)
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
            
     !----------------------------------------------------------------------
     ! 3-D Fields                  
     !----------------------------------------------------------------------
     REAL(fp), POINTER :: AREA_M2   (:,:,:) ! Grid box surface area [cm2]
     REAL(fp), POINTER :: AIRDEN    (:,:,:) ! Dry air density [kg/m3]
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
                                            ! (assume per moist air)
     REAL(fp), POINTER :: DQLDTMST  (:,:,:) ! H2O tendency, mst proc [kg/kg/s]
                                            ! (assume per moist air)
     REAL(fp), POINTER :: DQVDTMST  (:,:,:) ! Vapor tendency, mst proc [kg/kg/s]
                                            ! (assume per moist air)
     REAL(fp), POINTER :: DTRAIN    (:,:,:) ! Detrainment flux [kg/m2/s]
     REAL(fp), POINTER :: ENTRAIN   (:,:,:) ! GCAP entrainment [Pa/s]
     REAL(fp), POINTER :: HKBETA    (:,:,:) ! Hack overshoot parameter [1]
     REAL(fp), POINTER :: HKETA     (:,:,:) ! Hack conv mass flux [kg/m2/s]
     REAL(fp), POINTER :: MOISTQ    (:,:,:) ! Tendency in sp. humidity 
                                            ! [kg/kg moist air/s]
     REAL(fp), POINTER :: OPTD      (:,:,:) ! Visible optical depth [1]
     REAL(fp), POINTER :: PEDGE     (:,:,:) ! Wet air pressure [hPa] @ level 
                                            ! edges [hPa]
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
     REAL(fp), POINTER :: SPHU      (:,:,:) ! Specific humidity 
                                            ! [g water vapor / kg moist air]
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
     ! Note on pressures: PMID and PMEAN are calculated from PEDGE, 
     ! and dry air partial pressures assume constant RH and T across grid box
     REAL(fp), POINTER :: PEDGE_DRY (:,:,:) ! Dry air partial pressure [hPa] 
                                            ! @ level edges [hPa]
     REAL(fp), POINTER :: PMID      (:,:,:) ! Average wet air pressure [hPa]
                                            ! defined as arithmetic
                                            ! average of edge pressures
     REAL(fp), POINTER :: PMID_DRY  (:,:,:) ! Dry air partial pressure [hPa]
                                            ! defined as arithmetic average
                                            ! of edge pressures 
     REAL(fp), POINTER :: PMEAN     (:,:,:) ! Grid box mean wet air pressure 
                                            ! [hPa] calculated as P
                                            ! integrated over z
     REAL(fp), POINTER :: PMEAN_DRY (:,:,:) ! Grid box mean dry air partial 
                                            ! pressure calculated using water 
                                            ! vapor partial pressure derived
                                            ! from PMEAN
     REAL(fp), POINTER :: MOISTMW   (:,:,:) ! Moist air molec weight [g/mol]
                                            ! at the grid box altitude-weighted
                                            ! pressure level (temporary)        
     REAL(fp), POINTER :: TV        (:,:,:) ! Virtual temperature [K]
     REAL(fp), POINTER :: MAIRDEN   (:,:,:) ! Moist air density [kg/m3]
     REAL(fp), POINTER :: AVGW      (:,:,:) ! Water vapor volume mixing ratio
                                            ! [vol H2O / vol dry air]
     REAL(fp), POINTER :: BXHEIGHT  (:,:,:) ! Grid box height [m] (dry air)
     REAL(fp), POINTER :: DELP      (:,:,:) ! Delta-P extent of grid box [hPa]
                                            ! (wet air)
     REAL(fp), POINTER :: AD        (:,:,:) ! Dry air mass [kg] in grid box
     REAL(fp), POINTER :: ADMOIST   (:,:,:) ! Moist air mass [kg] in grid box
     REAL(fp), POINTER :: AIRVOL    (:,:,:) ! Grid box volume [m3] (dry air)
     REAL(fp), POINTER :: DELP_PREV (:,:,:) ! Previous State_Met%DELP
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
     REAL(fp), POINTER :: XCHLR     (:,:,:) ! MODIS CHLR per land type, this mo
     REAL(fp), POINTER :: LandTypeFrac(:,:,:) ! Olson frac per type (I,J,type)
     REAL(fp), POINTER :: XLAI_NATIVE(:,:,:)  ! avg LAI per type (I,J,type)
     REAL(fp), POINTER :: XCHLR_NATIVE(:,:,:) ! avg CHLR per type (I,J,type)

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
!  17 Mar 2016 - M. Sulprizio- Remove OPTDEP. Instead, we now solely use OPTD.
!  18 Oct 2016 - E. Lundgren - Remove XLAI2, CHLR2; add MODISLAI, MODISCHLR to
!                              replace modis_lai_mod-level GC_LAI and GC_CHLR
!  19 Oct 2016 - E. Lundgren - Use NSURFTYPE as the # of land types
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gigc_state_met
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_STATE\_MET allocates all fields of 
!  the Grid-Indpendent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_GIGC_State_Met( am_I_Root, IM, JM, LM, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
    USE CMN_SIZE_MOD,    ONLY : NSURFTYPE         ! # of land types

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: LX

    ! Assume success
    RC = GIGC_SUCCESS

    !=======================================================================
    ! Allocate 2-D Fields
    !=======================================================================
    ALLOCATE( State_Met%ALBD      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ALBD     = 0.0_fp

    ALLOCATE( State_Met%CLDFRC    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%CLDFRC   = 0.0_fp

    ALLOCATE( State_Met%CLDTOPS   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%CLDTOPS  = 0.0_fp

    ALLOCATE( State_Met%EFLUX     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%EFLUX    = 0.0_fp

    ALLOCATE( State_Met%EVAP      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%EVAP     = 0.0_fp

    ALLOCATE( State_Met%FRCLND    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRCLND   = 0.0_fp

    ALLOCATE( State_Met%FRLAKE    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLAKE   = 0.0_fp

    ALLOCATE( State_Met%FRLAND    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLAND   = 0.0_fp 

    ALLOCATE( State_Met%FRLANDIC  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLANDIC = 0.0_fp 

    ALLOCATE( State_Met%FROCEAN   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FROCEAN  = 0.0_fp
 
    ALLOCATE( State_Met%GRN       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GRN      = 0.0_fp 

    ALLOCATE( State_Met%GWETROOT  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GWETROOT = 0.0_fp 

    ALLOCATE( State_Met%GWETTOP   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GWETTOP  = 0.0_fp 

    ALLOCATE( State_Met%HFLUX     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HFLUX    = 0.0_fp 

    ALLOCATE( State_Met%LAI       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LAI      = 0.0_fp

    ALLOCATE( State_Met%ITY       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ITY      = 1.d0

    ALLOCATE( State_Met%LWI       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LWI      = 0.0_fp

    ALLOCATE( State_Met%PARDR     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PARDR    = 0.0_fp

    ALLOCATE( State_Met%PARDF     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PARDF    = 0.0_fp

    ALLOCATE( State_Met%PBLH      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PBLH     = 0.0_fp

    ALLOCATE( State_Met%PBL_TOP_L ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PBL_TOP_L = 0

    ALLOCATE( State_Met%PHIS      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PHIS     = 0.0_fp

    ALLOCATE( State_Met%PRECCON   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECCON  = 0.0_fp

    ALLOCATE( State_Met%PRECSNO   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECSNO  = 0.0_fp

    ALLOCATE( State_Met%PRECTOT   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECTOT  = 0.0_fp

    ALLOCATE( State_Met%PS1       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PS1      = 0.0_fp

    ALLOCATE( State_Met%PS2       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PS2      = 0.0_fp

    ALLOCATE( State_Met%PSC2      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PSC2     = 0.0_fp

    ALLOCATE( State_Met%RADLWG    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RADLWG   = 0.0_fp

    ALLOCATE( State_Met%RADSWG    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RADSWG   = 0.0_fp

    ALLOCATE( State_Met%SLP       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SLP      = 0.0_fp

    ALLOCATE( State_Met%SNODP     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNODP    = 0.0_fp

    ALLOCATE( State_Met%SNOMAS    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOMAS   = 0.0_fp

    ALLOCATE( State_Met%SST       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SST      = 0.0_fp

    ALLOCATE( State_Met%SUNCOS    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SUNCOS   = 0.0_fp

    ALLOCATE( State_Met%SUNCOSmid ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SUNCOSmid = 0.0_fp

!    ALLOCATE( State_Met%SUNCOSmid5( IM, JM ), STAT=RC )
!    IF ( RC /= GIGC_SUCCESS ) RETURN
!    State_Met%SUNCOSmid5 = 0.0_fp

    ALLOCATE( State_Met%SWGDN     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SWGDN    = 0.0_fp

    ALLOCATE( State_Met%TO3       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO3      = 0.0_fp

    ALLOCATE( State_Met%TROPP     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP    = 0.0_fp

    ALLOCATE( State_Met%TS        ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TS       = 0.0_fp

    ALLOCATE( State_Met%TSKIN     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TSKIN    = 0.0_fp

    ALLOCATE( State_Met%U10M      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%U10M     = 0.0_fp

    ALLOCATE( State_Met%USTAR     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%USTAR    = 0.0_fp

    ALLOCATE( State_Met%UVALBEDO  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UVALBEDO = 0.0_fp

    ALLOCATE( State_Met%V10M      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%V10M     = 0.0_fp

    ALLOCATE( State_Met%Z0        ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%Z0       = 0.0_fp

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !=======================================================================
    ALLOCATE( State_Met%LWI_GISS  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LWI_GISS = 0.0_fp

    ALLOCATE( State_Met%MOLENGTH  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%MOLENGTH = 0.0_fp

    ALLOCATE( State_Met%OICE      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%OICE     = 0.0_fp

    ALLOCATE( State_Met%SNICE     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNICE    = 0.0_fp

    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOW     = 0.0_fp

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0.0_fp

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0.0_fp

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !=======================================================================
    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOW     = 0.0_fp

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0.0_fp

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0.0_fp

#elif defined( GEOS_5 )

    !=======================================================================
    ! GEOS-5 met fields
    !=======================================================================
    ALLOCATE( State_Met%TO31      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO31     = 0.0_fp

    ALLOCATE( State_Met%TO32      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO32     = 0.0_fp

    ALLOCATE( State_Met%TTO3      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TTO3     = 0.0_fp

#elif defined( GEOS_FP ) || defined( MERRA ) || defined( MERRA2 )

    !=======================================================================
    ! MERRA, GEOS-FP, and MERRA2 met fields
    !=======================================================================
    ALLOCATE( State_Met%FRSEAICE  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRSEAICE = 0.0_fp

    ALLOCATE( State_Met%FRSNO     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRSNO    = 0.0_fp

    ALLOCATE( State_Met%PRECANV   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECANV  = 0.0_fp

    ALLOCATE( State_Met%PRECLSC   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECLSC  = 0.0_fp

    ALLOCATE( State_Met%SEAICE00  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE00 = 0.0_fp

    ALLOCATE( State_Met%SEAICE10  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE10 = 0.0_fp

    ALLOCATE( State_Met%SEAICE20  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE20 = 0.0_fp

    ALLOCATE( State_Met%SEAICE30  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE30 = 0.0_fp

    ALLOCATE( State_Met%SEAICE40  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE40 = 0.0_fp

    ALLOCATE( State_Met%SEAICE50  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE50 = 0.0_fp

    ALLOCATE( State_Met%SEAICE60  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE60 = 0.0_fp

    ALLOCATE( State_Met%SEAICE70  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE70 = 0.0_fp

    ALLOCATE( State_Met%SEAICE80  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE80 = 0.0_fp

    ALLOCATE( State_Met%SEAICE90  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE90 = 0.0_fp

#endif

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================
    ALLOCATE( State_Met%AD        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AD       = 0.0_fp

    ALLOCATE( State_Met%ADMOIST   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%ADMOIST  = 0.0_fp
                                               
    ALLOCATE( State_Met%AIRDEN    ( IM, JM, LM   ), STAT=RC )  
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AIRDEN   = 0.0_fp
    
    ALLOCATE( State_Met%MAIRDEN    ( IM, JM, LM   ), STAT=RC )  
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%MAIRDEN   = 0.0_fp
                              
    ALLOCATE( State_Met%AIRVOL    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AIRVOL   = 0.0_fp
                                               
!    ALLOCATE( State_Met%AREA_M2   ( IM, JM, LM   ), STAT=RC )
    ALLOCATE( State_Met%AREA_M2   ( IM, JM, 1    ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AREA_M2  = 0.0_fp

    ALLOCATE( State_Met%AVGW      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%AVGW     = 0.0_fp
                                               
    ALLOCATE( State_Met%BXHEIGHT  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%BXHEIGHT = 0.0_fp
                                               
    ALLOCATE( State_Met%CLDF      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%CLDF     = 0.0_fp
                                               
    ALLOCATE( State_Met%CMFMC     ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%CMFMC    = 0.0_fp
                                               
    ALLOCATE( State_Met%DELP      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DELP     = 0.0_fp

    ALLOCATE( State_Met%DELP_PREV ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DELP_PREV= 0.0_fp

    ALLOCATE( State_Met%DQRCU     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DQRCU    = 0.0_fp

    ALLOCATE( State_Met%DQRLSAN   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DQRLSAN  = 0.0_fp

    ALLOCATE( State_Met%DQIDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQIDTMST = 0.0_fp
                                               
    ALLOCATE( State_Met%DQLDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQLDTMST = 0.0_fp
                                               
    ALLOCATE( State_Met%DQVDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQVDTMST = 0.0_fp
                                               
    ALLOCATE( State_Met%DTRAIN    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DTRAIN   = 0.0_fp

    ALLOCATE( State_Met%MOISTMW ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%MOISTMW = 0.0_fp
                                               
    ALLOCATE( State_Met%MOISTQ    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%MOISTQ   = 0.0_fp
                                               
    ALLOCATE( State_Met%OPTD      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%OPTD     = 0.0_fp
                                               
    ALLOCATE( State_Met%PEDGE     ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PEDGE    = 0.0_fp

    ALLOCATE( State_Met%PEDGE_DRY ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PEDGE_DRY = 0.0_fp

    ALLOCATE( State_Met%PMID      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PMID     = 0.0_fp

    ALLOCATE( State_Met%PMID_DRY   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PMID_DRY  = 0.0_fp

    ALLOCATE( State_Met%PMEAN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PMEAN     = 0.0_fp

    ALLOCATE( State_Met%PMEAN_DRY   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PMEAN_DRY  = 0.0_fp

    ALLOCATE( State_Met%PV        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PV       = 0.0_fp

    ALLOCATE( State_Met%QI        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%QI       = 0.0_fp
                                               
    ALLOCATE( State_Met%QL        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%QL       = 0.0_fp
                         
    ALLOCATE( State_Met%RH        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN                                 
    State_Met%RH       = 0.0_fp
                        
    ALLOCATE( State_Met%SPHU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%SPHU     = 0.0_fp

    ALLOCATE( State_Met%SPHU_PREV ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%SPHU_PREV= 0.0_fp
                                               
    ALLOCATE( State_Met%T         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%T        = 0.0_fp

    ALLOCATE( State_Met%TV        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%TV       = 0.0_fp
                                               
    ALLOCATE( State_Met%TAUCLI    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%TAUCLI   = 0.0_fp
                                               
    ALLOCATE( State_Met%TAUCLW    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TAUCLW   = 0.0_fp
    
    ALLOCATE( State_Met%U         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%U        = 0.0_fp

    ALLOCATE( State_Met%V         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%V        = 0.0_fp

    ALLOCATE( State_Met%UPDVVEL   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UPDVVEL  = 0.0_fp

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !=======================================================================
    ALLOCATE( State_Met%DETRAINE  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DETRAINE = 0.0_fp

    ALLOCATE( State_Met%DETRAINN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DETRAINN = 0.0_fp

    ALLOCATE( State_Met%DNDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DNDE     = 0.0_fp

    ALLOCATE( State_Met%DNDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DNDN     = 0.0_fp

    ALLOCATE( State_Met%ENTRAIN   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ENTRAIN  = 0.0_fp

    ALLOCATE( State_Met%UPDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UPDE     = 0.0_fp

    ALLOCATE( State_Met%UPDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UPDN     = 0.0_fp

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !=======================================================================

    ALLOCATE( State_Met%HKBETA    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HKETA    = 0.0_fp

    ALLOCATE( State_Met%HKETA     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HKBETA   = 0.0_fp

    ALLOCATE( State_Met%ZMEU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ZMEU     = 0.0_fp

    ALLOCATE( State_Met%ZMMD      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ZMMD     = 0.0_fp

    ALLOCATE( State_Met%ZMMU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
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
    
    ALLOCATE( State_Met%PFICU     ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFICU    = 0.0_fp

    ALLOCATE( State_Met%PFILSAN   ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFILSAN  = 0.0_fp

    ALLOCATE( State_Met%PFLCU     ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFLCU    = 0.0_fp

    ALLOCATE( State_Met%PFLLSAN   ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFLLSAN  = 0.0_fp

    ALLOCATE( State_Met%REEVAPCN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%REEVAPCN = 0.0_fp

    ALLOCATE( State_Met%REEVAPLS  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%REEVAPLS = 0.0_fp

    ALLOCATE( State_Met%RH1       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RH1      = 0.0_fp

    ALLOCATE( State_Met%RH2       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RH2      = 0.0_fp

    ALLOCATE( State_Met%SPHU1     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SPHU1    = 0.0_fp

    ALLOCATE( State_Met%SPHU2     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SPHU2    = 0.0_fp

    ALLOCATE( State_Met%TMPU1     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TMPU1    = 0.0_fp

    ALLOCATE( State_Met%TMPU2     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TMPU2    = 0.0_fp

#endif

    !=======================================================================
    ! Allocate land type and leaf area index fields for dry deposition
    !=======================================================================
    ALLOCATE( State_Met%IREG       ( IM, JM           ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%IREG = 0

    ALLOCATE( State_Met%ILAND      ( IM, JM, NSURFTYPE), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ILAND = 0

    ALLOCATE( State_Met%IUSE       ( IM, JM, NSURFTYPE), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%IUSE = 0

    ALLOCATE( State_Met%XLAI       ( IM, JM, NSURFTYPE), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XLAI = 0.0_fp

    ALLOCATE( State_Met%MODISLAI   ( IM, JM           ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%MODISLAI = 0.0_fp

    ALLOCATE( State_Met%XCHLR      ( IM, JM, NSURFTYPE), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XCHLR = 0.0_fp

    ALLOCATE( State_Met%MODISCHLR  ( IM, JM           ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%MODISCHLR = 0.0_fp

    ALLOCATE( State_Met%LANDTYPEFRAC( IM, JM, NSURFTYPE ), STAT=RC )        
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LANDTYPEFRAC = 0.0_fp

    ALLOCATE( State_Met%XLAI_NATIVE ( IM, JM, NSURFTYPE ), STAT=RC )        
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XLAI_NATIVE  = 0.0_fp

    ALLOCATE( State_Met%XCHLR_NATIVE( IM, JM, NSURFTYPE ), STAT=RC )        
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XCHLR_NATIVE = 0.0_fp

  END SUBROUTINE Init_GIGC_State_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gigc_state_met
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_MET allocates all fields 
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Meteorology State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_GIGC_State_Met( am_I_Root, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
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
    ! Return success
    RC = GIGC_SUCCESS

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
    IF ( ASSOCIATED( State_Met%PS1        )) DEALLOCATE( State_Met%PS1        )
    IF ( ASSOCIATED( State_Met%PS2        )) DEALLOCATE( State_Met%PS2        )
    IF ( ASSOCIATED( State_Met%PSC2       )) DEALLOCATE( State_Met%PSC2       )
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

    ! SDE 2016-03-28: GCHP requires that these are nullified rather than being
    ! deallocated. Not yet sure why, but deallocating causes it to hang during
    ! cleanup
    ! 3-D fields
    IF ( ASSOCIATED( State_Met%AD         )) NULLIFY( State_Met%AD         )
    IF ( ASSOCIATED( State_Met%ADMOIST    )) NULLIFY( State_Met%ADMOIST    )
    IF ( ASSOCIATED( State_Met%AIRDEN     )) NULLIFY( State_Met%AIRDEN     )
    IF ( ASSOCIATED( State_Met%MAIRDEN    )) NULLIFY( State_Met%MAIRDEN    )
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
    IF ( ASSOCIATED( State_Met%MOISTMW    )) NULLIFY( State_Met%MOISTMW    )
    IF ( ASSOCIATED( State_Met%DTRAIN     )) NULLIFY( State_Met%DTRAIN     )
    IF ( ASSOCIATED( State_Met%MOISTQ     )) NULLIFY( State_Met%MOISTQ     )
    IF ( ASSOCIATED( State_Met%OPTD       )) NULLIFY( State_Met%OPTD       )
    IF ( ASSOCIATED( State_Met%PEDGE      )) NULLIFY( State_Met%PEDGE      )
    IF ( ASSOCIATED( State_Met%PEDGE_DRY  )) NULLIFY( State_Met%PEDGE_DRY  )
    IF ( ASSOCIATED( State_Met%PMID       )) NULLIFY( State_Met%PMID       )
    IF ( ASSOCIATED( State_Met%PMID_DRY   )) NULLIFY( State_Met%PMID_DRY   )
    IF ( ASSOCIATED( State_Met%PMEAN      )) NULLIFY( State_Met%PMEAN      )
    IF ( ASSOCIATED( State_Met%PMEAN_DRY  )) NULLIFY( State_Met%PMEAN_DRY  )
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
    IF ( ASSOCIATED( State_Met%MODISLAI   )) DEALLOCATE( State_Met%MODISLAI   )
    IF ( ASSOCIATED( State_Met%XCHLR      )) DEALLOCATE( State_Met%XCHLR      )
    IF ( ASSOCIATED( State_Met%MODISCHLR  )) DEALLOCATE( State_Met%MODISCHLR  )
    IF (ASSOCIATED( State_Met%LANDTYPEFRAC)) DEALLOCATE( State_Met%LANDTYPEFRAC)
    IF (ASSOCIATED( State_Met%XLAI_NATIVE )) DEALLOCATE( State_Met%XLAI_NATIVE )
    IF (ASSOCIATED( State_Met%XCHLR_NATIVE)) DEALLOCATE( State_Met%XCHLR_NATIVE)

   END SUBROUTINE Cleanup_GIGC_State_Met
!EOC
END MODULE GIGC_State_Met_Mod
