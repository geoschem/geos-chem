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
     REAL*8,  POINTER :: ALBD      (:,:  )  ! Visible surface albedo [1]
     REAL*8,  POINTER :: CLDFRC    (:,:  )  ! Column cloud fraction [1]
     INTEGER, POINTER :: CLDTOPS   (:,:  )  ! Max cloud top height [levels]
     REAL*8,  POINTER :: EFLUX     (:,:  )  ! Latent heat flux [W/m2]
     REAL*8,  POINTER :: EVAP      (:,:  )  ! Surface evap [kg/m2/s]
     REAL*8,  POINTER :: FRCLND    (:,:  )  ! Olson land fraction [1]
     REAL*8,  POINTER :: FRLAKE    (:,:  )  ! Fraction of lake [1]
     REAL*8,  POINTER :: FRLAND    (:,:  )  ! Fraction of land [1]
     REAL*8,  POINTER :: FRLANDIC  (:,:  )  ! Fraction of land ice [1]
     REAL*8,  POINTER :: FROCEAN   (:,:  )  ! Fraction of ocean [1]
     REAL*8,  POINTER :: FRSEAICE  (:,:  )  ! Sfc sea ice fraction
     REAL*8,  POINTER :: FRSNO     (:,:  )  ! Sfc snow fraction
     REAL*8,  POINTER :: GRN       (:,:  )  ! Greenness fraction
     REAL*8,  POINTER :: GWETROOT  (:,:  )  ! Root soil wetness [1]
     REAL*8,  POINTER :: GWETTOP   (:,:  )  ! Top soil moisture [1]
     REAL*8,  POINTER :: HFLUX     (:,:  )  ! Sensible heat flux [W/m2]
     REAL*8,  POINTER :: LAI       (:,:  )  ! Leaf area index [m2/m2]
     REAL*8,  POINTER :: LWI       (:,:  )  ! Land/water indices [1]
     REAL*8,  POINTER :: LWI_GISS  (:,:  )  ! Land fraction [1]
     REAL*8,  POINTER :: MOLENGTH  (:,:  )  ! Monin-Obhukov length [m]
     REAL*8,  POINTER :: OICE      (:,:  )  ! Fraction of ocean ice [1]
     REAL*8,  POINTER :: PARDR     (:,:  )  ! Direct  photsyn active rad [W/m2]
     REAL*8,  POINTER :: PARDF     (:,:  )  ! Diffuse photsyn active rad [W/m2]
     REAL*8,  POINTER :: PBLH      (:,:  )  ! PBL height [m]
     REAL*8,  POINTER :: PHIS      (:,:  )  ! Sfc geopotential height [m2/s2]
     REAL*8,  POINTER :: PRECANV   (:,:  )  ! Anvil previp @ ground [kg/m2/s]
     REAL*8,  POINTER :: PRECCON   (:,:  )  ! Conv  precip @ ground [kg/m2/s]
     REAL*8,  POINTER :: PRECTOT   (:,:  )  ! Total precip @ ground [kg/m2/s]
     REAL*8,  POINTER :: PRECLSC   (:,:  )  ! LS precip @ ground [kg/m2/s]
     REAL*8,  POINTER :: PRECSNO   (:,:  )  ! Snow precip [kg/m2/s]
     REAL*8,  POINTER :: PS1       (:,:  )  ! Sfc press at timestep start[hPa]
     REAL*8,  POINTER :: PS2       (:,:  )  ! Sfc press at timestep end [hPa]
     REAL*8,  POINTER :: PSC2      (:,:  )  ! Interpolated sfc pressure [hPa]
     REAL*8,  POINTER :: RADLWG    (:,:  )  ! Net LW radiation @ ground [W/m2]
     REAL*8,  POINTER :: RADSWG    (:,:  )  ! Solar radiation @ ground [W/m2]
     REAL*8,  POINTER :: SEAICE00  (:,:  )  ! Sea ice coverage 00-10%
     REAL*8,  POINTER :: SEAICE10  (:,:  )  ! Sea ice coverage 10-20%
     REAL*8,  POINTER :: SEAICE20  (:,:  )  !  Sea ice coverage 20-30%
     REAL*8,  POINTER :: SEAICE30  (:,:  )  ! Sea ice coverage 30-40%
     REAL*8,  POINTER :: SEAICE40  (:,:  )  ! Sea ice coverage 40-50%
     REAL*8,  POINTER :: SEAICE50  (:,:  )  ! Sea ice coverage 50-60%
     REAL*8,  POINTER :: SEAICE60  (:,:  )  ! Sea ice coverage 60-70%
     REAL*8,  POINTER :: SEAICE70  (:,:  )  ! Sea ice coverage 70-80%
     REAL*8,  POINTER :: SEAICE80  (:,:  )  ! Sea ice coverage 80-90%
     REAL*8,  POINTER :: SEAICE90  (:,:  )  ! Sea ice coverage 90-100%
     REAL*8,  POINTER :: SLP       (:,:  )  ! Sea level pressure [hPa]
     REAL*8,  POINTER :: SNICE     (:,:  )  ! Fraction of snow/ice [1]
     REAL*8,  POINTER :: SNODP     (:,:  )  ! Snow depth [m]
     REAL*8,  POINTER :: SNOMAS    (:,:  )  ! Snow mass [kg/m2]
     REAL*8,  POINTER :: SNOW      (:,:  )  ! Snow depth (H2O equiv) [mm H2O]
     REAL*8,  POINTER :: SST       (:,:  )  ! Sea surface temperature [K]
     REAL*8,  POINTER :: SUNCOS    (:,:  )  ! COS(SZA), current time
     REAL*8,  POINTER :: SUNCOSmid (:,:  )  ! COS(SZA), midpt of chem timestep
     REAL*8,  POINTER :: SUNCOSmid5(:,:  )  ! COS(SZA), midpt of chem timestep
                                            !  5 hrs ago (for PARANOX)
     REAL*8,  POINTER :: TO3       (:,:  )  ! Total overhead O3 column [DU]
     REAL*8,  POINTER :: TO31      (:,:  )  ! Total O3 at timestep start [DU]
     REAL*8,  POINTER :: TO32      (:,:  )  ! Total O3 at timestep end [DU]
     REAL*8,  POINTER :: TROPP     (:,:  )  ! Tropopause pressure [hPa]
     REAL*8,  POINTER :: TROPP1    (:,:  )  ! Trop P at timestep start [hPa]
     REAL*8,  POINTER :: TROPP2    (:,:  )  ! Trop P at timestep end [hPa]
     REAL*8,  POINTER :: TS        (:,:  )  ! Surface temperature [K]
     REAL*8,  POINTER :: TSKIN     (:,:  )  ! Surface skin temperature [K]
     REAL*8,  POINTER :: TTO3      (:,:  )  ! Tropospheric ozone column [DU]
     REAL*8,  POINTER :: U10M      (:,:  )  ! E/W wind speed @ 10m height [m/s]
     REAL*8,  POINTER :: USTAR     (:,:  )  ! Friction velocity [m/s]
     REAL*8,  POINTER :: UVALBEDO  (:,:  )  ! UV surface albedo [1]
     REAL*8,  POINTER :: V10M      (:,:  )  ! N/S wind speed @ 10m height [m/s]
     REAL*8,  POINTER :: Z0        (:,:  )  ! Surface roughness height [m]
            
     !----------------------------------------------------------------------
     ! 3-D Fields                  
     !----------------------------------------------------------------------
     REAL*8,  POINTER :: AD        (:,:,:)  ! Air mass [kg]
     REAL*8,  POINTER :: AIRDEN    (:,:,:)  ! Air density [kg/m3]
     REAL*8,  POINTER :: AIRVOL    (:,:,:)  ! Grid box volume [m3]
     REAL*8,  POINTER :: AREA_M2   (:,:,:)  ! Grid box surface area [cm2]
     REAL*8,  POINTER :: AVGW      (:,:,:)  ! Mixing ratio of water vapor
     REAL*8,  POINTER :: BXHEIGHT  (:,:,:)  ! Grid box height [m]
     REAL*8,  POINTER :: CLDF      (:,:,:)  ! 3-D cloud fraction [1]
     REAL*8,  POINTER :: CMFMC     (:,:,:)  ! Cloud mass flux [kg/m2/s]
     REAL*8,  POINTER :: DELP      (:,:,:)  ! Delta-P extent  of a grid box [mb]
     REAL*8,  POINTER :: DETRAINE  (:,:,:)  ! Detrainment (entrain plume)[Pa/s]
     REAL*8,  POINTER :: DETRAINN  (:,:,:)  ! Detrainment (non-entr plume)[Pa/s]
     REAL*8,  POINTER :: DNDE      (:,:,:)  ! Downdraft (entr plume) [Pa/s]
     REAL*8,  POINTER :: DNDN      (:,:,:)  ! Downdraft (non-entr plume) [Pa/s]
     REAL*8,  POINTER :: DQRCU     (:,:,:)  ! Conv precip prod rate [kg/kg/s]
     REAL*8,  POINTER :: DQRLSAN   (:,:,:)  ! LS precip prod rate [kg/kg/s]
     REAL*8,  POINTER :: DQIDTMST  (:,:,:)  ! Ice tendency, mst proc [kg/kg/s]
     REAL*8,  POINTER :: DQLDTMST  (:,:,:)  ! H2O tendency, mst proc [kg/kg/s]
     REAL*8,  POINTER :: DQVDTMST  (:,:,:)  ! Vapor tendency, mst proc [kg/kg/s]
     REAL*8,  POINTER :: DTRAIN    (:,:,:)  ! Detrainment flux [kg/m2/s]
     REAL*8,  POINTER :: ENTRAIN   (:,:,:)  ! GCAP entrainment [Pa/s]
     REAL*8,  POINTER :: HKBETA    (:,:,:)  ! Hack overshoot parameter [1]
     REAL*8,  POINTER :: HKETA     (:,:,:)  ! Hack conv mass flux [kg/m2/s]
     REAL*8,  POINTER :: MOISTQ    (:,:,:)  ! Tendency in sp. humidity [kg/kg/s]
     REAL*8,  POINTER :: OPTD      (:,:,:)  ! Visible optical depth [1]
     REAL*8,  POINTER :: OPTDEP    (:,:,:)  ! Visible optical depth [1]
     REAL*8,  POINTER :: PEDGE     (:,:,:)  ! Pressure @ level edges [Pa]
     REAL*8,  POINTER :: PMID      (:,:,:)  ! Pressure @ level centers [Pa]
     REAL*8,  POINTER :: PFICU     (:,:,:)  ! Dwn flux ice prec:conv [kg/m2/s]
     REAL*8,  POINTER :: PFILSAN   (:,:,:)  ! Dwn flux ice prec:LS+anv [kg/m2/s]
     REAL*8,  POINTER :: PFLCU     (:,:,:)  ! Dwn flux liq prec:conv [kg/m2/s]
     REAL*8,  POINTER :: PFLLSAN   (:,:,:)  ! Dwn flux ice prec:LS+anv [kg/m2/s]
     REAL*8,  POINTER :: PV        (:,:,:)  ! Potential vort [kg*m2/kg/s]
     REAL*8,  POINTER :: QI        (:,:,:)  ! Ice mixing ratio [kg/kg]
     REAL*8,  POINTER :: QL        (:,:,:)  ! Water mixing ratio [kg/kg]
     REAL*8,  POINTER :: REEVAPCN  (:,:,:)  ! Evap of precip conv [kg/kg/s]
     REAL*8,  POINTER :: REEVAPLS  (:,:,:)  ! Evap of precip LS+anvil [kg/kg/s]
     REAL*8,  POINTER :: RH        (:,:,:)  ! Relative humidity [%]
     REAL*8,  POINTER :: RH1       (:,:,:)  ! RH at timestep start [%]
     REAL*8,  POINTER :: RH2       (:,:,:)  ! RH at timestep end [%]
     REAL*8,  POINTER :: SPHU      (:,:,:)  ! Specific humidity [kg/kg]
     REAL*8,  POINTER :: SPHU1     (:,:,:)  ! Spec hum at timestep start [kg/kg]
     REAL*8,  POINTER :: SPHU2     (:,:,:)  ! Spec hum at timestep end [kg/kg]
     REAL*8,  POINTER :: T         (:,:,:)  ! Temperature [K]
     REAL*8,  POINTER :: TAUCLI    (:,:,:)  ! Opt depth of ice clouds [1]
     REAL*8,  POINTER :: TAUCLW    (:,:,:)  ! Opt depth of H2O clouds [1]
     REAL*8,  POINTER :: TMPU1     (:,:,:)  ! Temperature at timestep start [K]
     REAL*8,  POINTER :: TMPU2     (:,:,:)  ! Temperature at timestep end [K]
     REAL*8,  POINTER :: U         (:,:,:)  ! E/W component of wind [m s-1]
     REAL*8,  POINTER :: UPDE      (:,:,:)  ! Updraft (entraining plume) [Pa/s]
     REAL*8,  POINTER :: UPDN      (:,:,:)  ! Updraft (non-entr'n plume) [Pa/s]
     REAL*8,  POINTER :: V         (:,:,:)  ! N/S component of wind [m s-1]
     REAL*8,  POINTER :: ZMEU      (:,:,:)  ! Z/M updraft entrainment [Pa/s]
     REAL*8,  POINTER :: ZMMD      (:,:,:)  ! Z/M downdraft mass flux [Pa/s]
     REAL*8,  POINTER :: ZMMU      (:,:,:)  ! Z/M updraft   mass flux [Pa/s]

     !----------------------------------------------------------------------
     ! Land type and leaf area index (LAI) fields for dry deposition
     !----------------------------------------------------------------------
     INTEGER, POINTER :: IREG      (:,:  )  ! # of landtypes in grid box (I,J) 
     INTEGER, POINTER :: ILAND     (:,:,:)  ! Land type at (I,J); 1..IREG(I,J)
     INTEGER, POINTER :: IUSE      (:,:,:)  ! Fraction (per mil) of grid box
                                            !  (I,J) occupied by each land type
     REAL*8,  POINTER :: XLAI      (:,:,:)  ! LAI per land type, this month
     REAL*8,  POINTER :: XLAI2     (:,:,:)  ! LAI per land type, next month

  END TYPE MetState
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, split off from gc_type_mod.F90
!  23 Oct 2012 - R. Yantosca - Added QI, QL met fields to the derived type
!  15 Nov 2012 - M. Payer    - Added all remaining met fields
!  12 Dec 2012 - R. Yantosca - Add IREG, ILAND, IUSE fields for dry deposition
!  13 Dec 2012 - R. Yantosca - Add XLAI, XLAI2 fields for dry deposition
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  15 Nov 2013 - R. Yantosca - Now denote that RH fields have units of [%]
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
    USE CMN_SIZE_MOD,    ONLY : NTYPE            ! # of land types

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
    State_Met%ALBD     = 0d0

    ALLOCATE( State_Met%CLDFRC    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%CLDFRC   = 0d0

    ALLOCATE( State_Met%CLDTOPS   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%CLDTOPS  = 0d0

    ALLOCATE( State_Met%EFLUX     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%EFLUX    = 0d0

    ALLOCATE( State_Met%EVAP      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%EVAP     = 0d0

    ALLOCATE( State_Met%FRCLND    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRCLND   = 0d0

    ALLOCATE( State_Met%FRLAKE    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLAKE   = 0d0

    ALLOCATE( State_Met%FRLAND    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLAND   = 0d0 

    ALLOCATE( State_Met%FRLANDIC  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRLANDIC = 0d0 

    ALLOCATE( State_Met%FROCEAN   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FROCEAN  = 0d0
 
    ALLOCATE( State_Met%GRN       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GRN      = 0d0 

    ALLOCATE( State_Met%GWETROOT  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GWETROOT = 0d0 

    ALLOCATE( State_Met%GWETTOP   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%GWETTOP  = 0d0 

    ALLOCATE( State_Met%HFLUX     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HFLUX    = 0d0 

    ALLOCATE( State_Met%LAI       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LAI      = 0d0

    ALLOCATE( State_Met%LWI       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LWI      = 0d0

    ALLOCATE( State_Met%PARDR     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PARDR    = 0d0

    ALLOCATE( State_Met%PARDF     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PARDF    = 0d0

    ALLOCATE( State_Met%PBLH      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PBLH     = 0d0

    ALLOCATE( State_Met%PHIS      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PHIS     = 0d0

    ALLOCATE( State_Met%PRECCON   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECCON  = 0d0

    ALLOCATE( State_Met%PRECSNO   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECSNO  = 0d0

    ALLOCATE( State_Met%PRECTOT   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECTOT  = 0d0

    ALLOCATE( State_Met%PS1       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PS1      = 0d0

    ALLOCATE( State_Met%PS2       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PS2      = 0d0

    ALLOCATE( State_Met%PSC2      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PSC2     = 0d0

    ALLOCATE( State_Met%RADLWG    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RADLWG   = 0d0

    ALLOCATE( State_Met%RADSWG    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RADSWG   = 0d0

    ALLOCATE( State_Met%SLP       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SLP      = 0d0

    ALLOCATE( State_Met%SNODP     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNODP    = 0d0

    ALLOCATE( State_Met%SNOMAS    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOMAS   = 0d0

    ALLOCATE( State_Met%SST       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SST      = 0d0

    ALLOCATE( State_Met%SUNCOS    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SUNCOS   = 0d0

    ALLOCATE( State_Met%SUNCOSmid ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SUNCOSmid = 0d0

    ALLOCATE( State_Met%SUNCOSmid5( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SUNCOSmid5 = 0d0

    ALLOCATE( State_Met%TO3       ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO3      = 0d0

    ALLOCATE( State_Met%TROPP     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP    = 0d0

    ALLOCATE( State_Met%TS        ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TS       = 0d0

    ALLOCATE( State_Met%TSKIN     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TSKIN    = 0d0

    ALLOCATE( State_Met%U10M      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%U10M     = 0d0

    ALLOCATE( State_Met%USTAR     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%USTAR    = 0d0

    ALLOCATE( State_Met%UVALBEDO  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UVALBEDO = 0d0

    ALLOCATE( State_Met%V10M      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%V10M     = 0d0

    ALLOCATE( State_Met%Z0        ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%Z0       = 0d0

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !=======================================================================
    ALLOCATE( State_Met%LWI_GISS  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%LWI_GISS = 0d0

    ALLOCATE( State_Met%MOLENGTH  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%MOLENGTH = 0d0

    ALLOCATE( State_Met%OICE      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%OICE     = 0d0

    ALLOCATE( State_Met%SNICE     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNICE    = 0d0

    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOW     = 0d0

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0d0

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0d0

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !=======================================================================
    ALLOCATE( State_Met%SNOW      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SNOW     = 0d0

    ALLOCATE( State_Met%TROPP1    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP1   = 0d0

    ALLOCATE( State_Met%TROPP2    ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TROPP2   = 0d0

#elif defined( GEOS_5 )

    !=======================================================================
    ! GEOS-5 met fields
    !=======================================================================
    ALLOCATE( State_Met%TO31      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO31     = 0d0

    ALLOCATE( State_Met%TO32      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TO32     = 0d0

    ALLOCATE( State_Met%TTO3      ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TTO3     = 0d0

#elif defined( GEOS_FP ) || defined( MERRA )

    !=======================================================================
    ! GEOS-5.7.x / MERRA met fields
    !=======================================================================
    ALLOCATE( State_Met%FRSEAICE  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRSEAICE = 0d0

    ALLOCATE( State_Met%FRSNO     ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%FRSNO    = 0d0

    ALLOCATE( State_Met%PRECANV   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECANV  = 0d0

    ALLOCATE( State_Met%PRECLSC   ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PRECLSC  = 0d0

    ALLOCATE( State_Met%SEAICE00  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE00 = 0d0

    ALLOCATE( State_Met%SEAICE10  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE10 = 0d0

    ALLOCATE( State_Met%SEAICE20  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE20 = 0d0

    ALLOCATE( State_Met%SEAICE30  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE30 = 0d0

    ALLOCATE( State_Met%SEAICE40  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE40 = 0d0

    ALLOCATE( State_Met%SEAICE50  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE50 = 0d0

    ALLOCATE( State_Met%SEAICE60  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE60 = 0d0

    ALLOCATE( State_Met%SEAICE70  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE70 = 0d0

    ALLOCATE( State_Met%SEAICE80  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE80 = 0d0

    ALLOCATE( State_Met%SEAICE90  ( IM, JM ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SEAICE90 = 0d0

#endif

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================
    ALLOCATE( State_Met%AD        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AD       = 0d0
                                               
    ALLOCATE( State_Met%AIRDEN    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AIRDEN   = 0d0
                              
    ALLOCATE( State_Met%AIRVOL    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AIRVOL   = 0d0
                                               
    ALLOCATE( State_Met%AREA_M2   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%AREA_M2  = 0d0

    ALLOCATE( State_Met%AVGW      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%AVGW     = 0d0
                                               
    ALLOCATE( State_Met%BXHEIGHT  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%BXHEIGHT = 0d0
                                               
    ALLOCATE( State_Met%CLDF      ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%CLDF     = 0d0
                                               
    ALLOCATE( State_Met%CMFMC     ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%CMFMC    = 0d0
                                               
    ALLOCATE( State_Met%DELP      ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DELP     = 0d0

    ALLOCATE( State_Met%DQRCU     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DQRCU    = 0d0

    ALLOCATE( State_Met%DQRLSAN   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DQRLSAN  = 0d0

    ALLOCATE( State_Met%DQIDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQIDTMST = 0d0
                                               
    ALLOCATE( State_Met%DQLDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQLDTMST = 0d0
                                               
    ALLOCATE( State_Met%DQVDTMST  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DQVDTMST = 0d0
                                               
    ALLOCATE( State_Met%DTRAIN    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%DTRAIN   = 0d0
                                               
    ALLOCATE( State_Met%MOISTQ    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%MOISTQ   = 0d0
                                               
    ALLOCATE( State_Met%OPTD      ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%OPTD     = 0d0
                                               
    ALLOCATE( State_Met%OPTDEP    ( LM, IM, JM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%OPTDEP   = 0d0

    ALLOCATE( State_Met%PEDGE      ( IM, JM, LM+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PEDGE    = 0d0
                                               
    ALLOCATE( State_Met%PMID       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%PMID     = 0d0

    ALLOCATE( State_Met%PV        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PV       = 0d0

    ALLOCATE( State_Met%QI        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%QI       = 0d0
                                               
    ALLOCATE( State_Met%QL        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%QL       = 0d0
                         
    ALLOCATE( State_Met%RH        ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN                                 
    State_Met%RH       = 0d0
                        
    ALLOCATE( State_Met%SPHU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%SPHU     = 0d0
                                               
    ALLOCATE( State_Met%T         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%T        = 0d0
                                               
    ALLOCATE( State_Met%TAUCLI    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN           
    State_Met%TAUCLI   = 0d0
                                               
    ALLOCATE( State_Met%TAUCLW    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TAUCLW   = 0d0
    
    ALLOCATE( State_Met%U         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%U        = 0d0

    ALLOCATE( State_Met%V         ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%V        = 0d0

#if defined( GCAP )

    !=======================================================================
    ! GCAP met fields
    !=======================================================================
    ALLOCATE( State_Met%DETRAINE  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DETRAINE = 0d0

    ALLOCATE( State_Met%DETRAINN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DETRAINN = 0d0

    ALLOCATE( State_Met%DNDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DNDE     = 0d0

    ALLOCATE( State_Met%DNDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%DNDN     = 0d0

    ALLOCATE( State_Met%ENTRAIN   ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ENTRAIN  = 0d0

    ALLOCATE( State_Met%UPDE      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UPDE     = 0d0

    ALLOCATE( State_Met%UPDN      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%UPDN     = 0d0

#elif defined( GEOS_4 )

    !=======================================================================
    ! GEOS-4 met fields
    !=======================================================================

    ALLOCATE( State_Met%HKBETA    ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HKETA    = 0d0

    ALLOCATE( State_Met%HKETA     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%HKBETA   = 0d0

    ALLOCATE( State_Met%ZMEU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ZMEU     = 0d0

    ALLOCATE( State_Met%ZMMD      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ZMMD     = 0d0

    ALLOCATE( State_Met%ZMMU      ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ZMMU     = 0d0

#elif defined( GEOS_FP ) || defined( MERRA )

    !=======================================================================
    ! GEOS-5.7.x / MERRA met fields
    !=======================================================================

    ! Pick the proper vertical dimension
#if defined( GEOS_FP )
    LX = LM + 1           ! For fields that are on level edges
#else
    LX = LM               ! For fields that are on level centers
#endif
    
    ALLOCATE( State_Met%PFICU     ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFICU    = 0d0

    ALLOCATE( State_Met%PFILSAN   ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFILSAN  = 0d0

    ALLOCATE( State_Met%PFLCU     ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFLCU    = 0d0

    ALLOCATE( State_Met%PFLLSAN   ( IM, JM, LX   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%PFLLSAN  = 0d0

    ALLOCATE( State_Met%REEVAPCN  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%REEVAPCN = 0d0

    ALLOCATE( State_Met%REEVAPLS  ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%REEVAPLS = 0d0

    ALLOCATE( State_Met%RH1       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RH1      = 0d0

    ALLOCATE( State_Met%RH2       ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%RH2      = 0d0

    ALLOCATE( State_Met%SPHU1     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SPHU1    = 0d0

    ALLOCATE( State_Met%SPHU2     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%SPHU2    = 0d0

    ALLOCATE( State_Met%TMPU1     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TMPU1    = 0d0

    ALLOCATE( State_Met%TMPU2     ( IM, JM, LM   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%TMPU2    = 0d0

#endif

    !=======================================================================
    ! Allocate land type and leaf area index fields for dry deposition
    !=======================================================================
    ALLOCATE( State_Met%IREG      ( IM, JM        ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%IREG     = 0

    ALLOCATE( State_Met%ILAND     ( IM, JM, NTYPE ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%ILAND    = 0

    ALLOCATE( State_Met%IUSE      ( IM, JM, NTYPE ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%IUSE     = 0

    ALLOCATE( State_Met%XLAI      ( IM, JM, NTYPE ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XLAI     = 0d0

    ALLOCATE( State_Met%XLAI2     ( IM, JM, NTYPE ), STAT=RC )        
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Met%XLAI2    = 0d0

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
    IF ( ASSOCIATED( State_Met%SUNCOSmid5 )) DEALLOCATE( State_Met%SUNCOSmid5 )
    IF ( ASSOCIATED( State_Met%TROPP      )) DEALLOCATE( State_Met%TROPP      )
    IF ( ASSOCIATED( State_Met%TS         )) DEALLOCATE( State_Met%TS         )
    IF ( ASSOCIATED( State_Met%TSKIN      )) DEALLOCATE( State_Met%TSKIN      )
    IF ( ASSOCIATED( State_Met%TO3        )) DEALLOCATE( State_Met%TO3        )
    IF ( ASSOCIATED( State_Met%U10M       )) DEALLOCATE( State_Met%U10M       )
    IF ( ASSOCIATED( State_Met%USTAR      )) DEALLOCATE( State_Met%USTAR      )
    IF ( ASSOCIATED( State_Met%UVALBEDO   )) DEALLOCATE( State_Met%V10M       )
    IF ( ASSOCIATED( State_Met%V10M       )) DEALLOCATE( State_Met%UVALBEDO   )
    IF ( ASSOCIATED( State_Met%Z0         )) DEALLOCATE( State_Met%Z0         )

    ! 3-D fields
    IF ( ASSOCIATED( State_Met%AD         )) DEALLOCATE( State_Met%AD         )
    IF ( ASSOCIATED( State_Met%AIRDEN     )) DEALLOCATE( State_Met%AIRDEN     )
    IF ( ASSOCIATED( State_Met%AIRVOL     )) DEALLOCATE( State_Met%AIRVOL     )
    IF ( ASSOCIATED( State_Met%AREA_M2    )) DEALLOCATE( State_Met%AREA_M2    )
    IF ( ASSOCIATED( State_Met%AVGW       )) DEALLOCATE( State_Met%AVGW       )
    IF ( ASSOCIATED( State_Met%BXHEIGHT   )) DEALLOCATE( State_Met%BXHEIGHT   )
    IF ( ASSOCIATED( State_Met%CLDF       )) DEALLOCATE( State_Met%CLDF       )
    IF ( ASSOCIATED( State_Met%CMFMC      )) DEALLOCATE( State_Met%CMFMC      )
    IF ( ASSOCIATED( State_Met%DELP       )) DEALLOCATE( State_Met%DELP       )
    IF ( ASSOCIATED( State_Met%DQRCU      )) DEALLOCATE( State_Met%DQRCU      )
    IF ( ASSOCIATED( State_Met%DQRLSAN    )) DEALLOCATE( State_Met%DQRLSAN    )
    IF ( ASSOCIATED( State_Met%DQIDTMST   )) DEALLOCATE( State_Met%DQIDTMST   )
    IF ( ASSOCIATED( State_Met%DQLDTMST   )) DEALLOCATE( State_Met%DQLDTMST   )
    IF ( ASSOCIATED( State_Met%DQVDTMST   )) DEALLOCATE( State_Met%DQVDTMST   )
    IF ( ASSOCIATED( State_Met%DTRAIN     )) DEALLOCATE( State_Met%DTRAIN     )
    IF ( ASSOCIATED( State_Met%MOISTQ     )) DEALLOCATE( State_Met%MOISTQ     )
    IF ( ASSOCIATED( State_Met%OPTD       )) DEALLOCATE( State_Met%OPTD       )
    IF ( ASSOCIATED( State_Met%OPTDEP     )) DEALLOCATE( State_Met%OPTDEP     )
    IF ( ASSOCIATED( State_Met%PEDGE      )) DEALLOCATE( State_Met%PEDGE      )
    IF ( ASSOCIATED( State_Met%PMID       )) DEALLOCATE( State_Met%PMID       )
    IF ( ASSOCIATED( State_Met%PV         )) DEALLOCATE( State_Met%PV         )
    IF ( ASSOCIATED( State_Met%QI         )) DEALLOCATE( State_Met%QI         )
    IF ( ASSOCIATED( State_Met%QL         )) DEALLOCATE( State_Met%QL         )
    IF ( ASSOCIATED( State_Met%RH         )) DEALLOCATE( State_Met%RH         )
    IF ( ASSOCIATED( State_Met%SPHU       )) DEALLOCATE( State_Met%SPHU       )
    IF ( ASSOCIATED( State_Met%T          )) DEALLOCATE( State_Met%T          )
    IF ( ASSOCIATED( State_Met%TAUCLI     )) DEALLOCATE( State_Met%TAUCLI     )
    IF ( ASSOCIATED( State_Met%TAUCLW     )) DEALLOCATE( State_Met%TAUCLW     ) 
    IF ( ASSOCIATED( State_Met%U          )) DEALLOCATE( State_Met%U          )
    IF ( ASSOCIATED( State_Met%V          )) DEALLOCATE( State_Met%V          )

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

#elif defined( GEOS_FP ) || defined( MERRA )
    !========================================================================
    ! Fields specific to the GMAO MERRA and GEOS-5.7.x met data products
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

   END SUBROUTINE Cleanup_GIGC_State_Met
!EOC
END MODULE GIGC_State_Met_Mod
