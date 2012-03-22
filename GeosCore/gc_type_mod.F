!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_type_mod.f
!
! !DESCRIPTION: Module GC\_TYPE\_MOD contains derived type definitions for 
!  GEOS-Chem.  These definitions are used to create objects for:
!  \begin{itemize}
!  \item ID flags for chemical species
!  \item ID flags for advected tracers
!  \item Coefficients \& other quantities that translate between chemical
!        species and advected tracers
!  \item ID flags for advected tracers that dry deposit
!  \item Logical flags for selecting the various GEOS-Chem options
!  \item GEOS-Chem columnized meteorological fields and related quantities
!  \item Error traceback output
!  \end{itemize}
!
! !INTERFACE:
!
      MODULE GC_TYPE_MOD
!
! USES:
!
      IMPLICIT NONE
#     include "define.h"
      PUBLIC
!
! !PUBLIC TYPES:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Leave these types commented out for now ... activate them later
!%%% (bmy, 9/29/10)
!
!      !-----------------------------------------------------------------
!      ! ID flags for SMVGEAR species
!      !-----------------------------------------------------------------
!      TYPE :: ID_SPEC
!         INTEGER          :: A3O2,    ACET,    ACTA,    ALD2,    ALK4
!         INTEGER          :: ATO2,    B3O2,    C2H6,    C3H8,    CH2O
!         INTEGER          :: CH4,     CO,      CO2,     DMS,     DRYCH2O
!         INTEGER          :: DRYDEP,  DRYH2O2, DRYHNO3, DRYN2O5, DRYNO2
!         INTEGER          :: DRYO3,   DRYPAN,  DRYPMN,  DRYPPN,  DRYR4N2
!         INTEGER          :: EMISSION,EOH,     ETO2,    ETP,     GCO3
!         INTEGER          :: GLCO3,   GLP,     GLPAN,   GLYC,    GLYX
!         INTEGER          :: GP,      GPAN,    H,       H2,      H2O
!         INTEGER          :: H2O2,    HAC,     HCOOH,   HNO2,    HNO3
!         INTEGER          :: HNO4,    HO2,     IALD,    IAO2,    IAP
!         INTEGER          :: INO2,    INPN,    ISN1,    ISNO3,   ISNP
!         INTEGER          :: ISOP,    KO2,     LISOPOH, M,       MACR
!         INTEGER          :: MAN2,    MAO3,    MAOP,    MAP,     MCO3
!         INTEGER          :: MEK,     MGLY,    MNO3,    MO2,     MOH
!         INTEGER          :: MP,      MRO2,    MRP,     MSA,     MVK
!         INTEGER          :: MVN2,    N2,      N2O,     N2O5,    NH2
!         INTEGER          :: NH3,     NO,      NO2,     NO3,     O       
!         INTEGER          :: O1D,     O2,      O2CH2OH, O3,      OH
!         INTEGER          :: PAN,     PMN,     PO2,     PP,      PPN
!         INTEGER          :: PRN1,    PRPE,    PRPN,    R4N1,    R4N2    
!         INTEGER          :: R4O2,    R4P,     RA3P,    RB3P,    RCHO
!         INTEGER          :: RCO3,    RCOOH,   RIO1,    RIO2,    RIP
!         INTEGER          :: ROH,     RP,      SO2,     SO4,     VRO2
!         INTEGER          :: VRP
!      END TYPE ID_SPEC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !-----------------------------------------------------------------
      ! ID flags for advected tracers
      !-----------------------------------------------------------------
      TYPE :: ID_TRAC
         INTEGER          :: NOx,     Ox,      PAN,     CO,      ALK4  
         INTEGER          :: ISOP,    HNO3,    H2O2,    ACET,    MEK   
         INTEGER          :: ALD2,    RCHO,    MVK,     MACR,    PMN
         INTEGER          :: PPN,     R4N2,    PRPE,    C3H8,    CH2O
         INTEGER          :: C2H6,    N2O5,    HNO4,    MP,      DMS
         INTEGER          :: SO2,     SO4,     SO4s,    MSA,     NH3
         INTEGER          :: NH4,     NIT,     NITs,    BCPI,    BCPO
         INTEGER          :: OCPI,    OCPO,    ALPH,    LIMO,    ALCO
         INTEGER          :: SOG1,    SOG2,    SOG3,    SOG4,    SOA1
         INTEGER          :: SOA2,    SOA3,    SOA4,    DST1,    DST2
         INTEGER          :: DST3,    DST4,    SALA,    SALC
         INTEGER          :: Hg0,     Hg2,     HgP
      END TYPE ID_TRAC

      !-----------------------------------------------------------------
      ! Coefficients & arrays that link species & tracers
      !-----------------------------------------------------------------
      TYPE :: SPEC_2_TRAC
         REAL*8,  POINTER :: SPEC_COEF(:,:)
         INTEGER, POINTER :: SPEC_ID(:,:)
         INTEGER, POINTER :: SPEC_EMITTED(:)
         INTEGER, POINTER :: SPEC_PER_TRAC(:)
         REAL*8,  POINTER :: TRAC_COEF(:)
         REAL*8,  POINTER :: MOLWT_KG(:)
         REAL*8,  POINTER :: XNUMOL(:)
      END TYPE SPEC_2_TRAC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Leave these types commented out for now ... activate them later
!%%% (bmy, 9/29/10)
!      !-----------------------------------------------------------------
!      ! ID flags for advected tracers that dry deposit
!      !-----------------------------------------------------------------
!      TYPE :: ID_DRYD
!         INTEGER          :: NOx,     Ox,      PAN,     HNO3,    H2O2
!         INTEGER          :: PMN,     PPN,     R4N2,    CH2O,    N2O5
!         INTEGER          :: SO2,     SO4,     SO4s,    MSA,     NH3
!         INTEGER          :: NH4,     NIT,     NITs,    BCPI,    BCPO
!         INTEGER          :: OCPI,    OCPO,    ALPH,    LIMO,    ALCO
!         INTEGER          :: SOG1,    SOG2,    SOG3,    SOG4,    SOA1
!         INTEGER          :: SOA2,    SOA3,    SOA4,    DST1,    DST2
!         INTEGER          :: DST3,    DST4,    SALA,    SALC
!      END TYPE ID_DRYD
!
!      !-----------------------------------------------------------------
!      ! ID flags for advected tracers that wet deposit
!      !-----------------------------------------------------------------
!      TYPE :: ID_WETD
!         INTEGER          :: HNO3,    H2O2,    CH2O,    MP  
!         INTEGER          :: SO2,     SO4,     SO4s,    MSA,     NH3
!         INTEGER          :: NH4,     NIT,     NITs,    BCPI,    OCPI
!         INTEGER          :: BCPO,    OCPO,    ALPH,    LIMO,    ALCO
!         INTEGER          :: SOG1,    SOG2,    SOG3,    SOG4,    SOA1
!         INTEGER          :: SOA2,    SOA3,    SOA4,    DST1,    DST2
!         INTEGER          :: DST3,    DST4,    SALA,    SALC
!      END TYPE ID_WETD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !-----------------------------------------------------------------
      ! Logical flags that turn various options on/off
      !-----------------------------------------------------------------
      TYPE :: GC_OPTIONS
         LOGICAL          :: USE_ANTHRO
         LOGICAL          :: USE_ANTHRO_BRAVO  
         LOGICAL          :: USE_ANTHRO_CAC       
         LOGICAL          :: USE_ANTHRO_EDGAR     
         LOGICAL          :: USE_ANTHRO_EPA 
         LOGICAL          :: USE_ANTHRO_VISTAS    
         LOGICAL          :: USE_ANTHRO_EMEP      
         LOGICAL          :: USE_BIOGENIC
         LOGICAL          :: USE_BIOMASS   
         LOGICAL          :: USE_BIOMASS_GFED2
         LOGICAL          :: USE_CARBON_AEROSOLS
         LOGICAL          :: USE_CHEMISTRY
         LOGICAL          :: USE_CONVECTION
         LOGICAL          :: USE_DEAD_DUST
         LOGICAL          :: USE_DEBUG_PRINT
         LOGICAL          :: USE_DRYDEP
         LOGICAL          :: USE_DUST_AEROSOLS
         LOGICAL          :: USE_EMISSIONS
         LOGICAL          :: USE_NOx_AIRCRAFT
         LOGICAL          :: USE_NOx_LIGHTNING
         LOGICAL          :: USE_NOx_SOIL
         LOGICAL          :: USE_PBL_MIXING
         LOGICAL          :: USE_SEC_ORG_AEROSOLS
         LOGICAL          :: USE_SHIP_ARCTAS   
         LOGICAL          :: USE_SEASALT_AEROSOLS
         LOGICAL          :: USE_SULFATE_AEROSOLS
         LOGICAL          :: USE_WETDEP
         LOGICAL          :: USE_Hg
         LOGICAL          :: USE_Hg_DYNOCEAN
         LOGICAL          :: USE_DIAG14
         LOGICAL          :: USE_DIAG38
      END TYPE GC_OPTIONS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Leave these types commented out for now ... activate them later
!%%% (bmy, 9/29/10)
!      !-----------------------------------------------------------------
!      ! Time & date values
!      !-----------------------------------------------------------------
!      TYPE :: GC_TIME
!         INTEGER          :: YEAR         ! Current year (YYYY)
!         INTEGER          :: MONTH        ! Current month (1-12)
!         INTEGER          :: DAY          ! Current day (1-31)
!         INTEGER          :: DOY          ! Day of year (0-365/366)
!         INTEGER          :: HOUR         ! Current hour (0-23)
!         INTEGER          :: MINUTE       ! Current minute (0-59)
!         LOGICAL          :: FIRST_TIME   ! Is it the first timestep?
!         REAL*8           :: T_ELAPSED    ! Elapsed simulation time [min]
!         REAL*8           :: TS_DYN       ! Dynamic timestep [min]
!         REAL*8           :: TS_CHEM      ! Chemistry timestep [min]
!      END TYPE GC_TIME
!
!      !-----------------------------------------------------------------
!      ! Geographic location
!      !-----------------------------------------------------------------
!      TYPE :: GC_GEOLOC
!         REAL*8           :: LON          ! Longitude [degrees]
!         REAL*8           :: LAT          ! Latitude  [degrees]
!         REAL*8           :: LOCALTIME    ! Local solar time [hrs]
!      END TYPE GC_GEOLOC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !-----------------------------------------------------------------
      ! Dimension information
      !-----------------------------------------------------------------
      TYPE :: GC_DIMS  
         INTEGER          :: L_COLUMN     ! # of boxes in the vertical column
         INTEGER          :: N_AER        ! # of aerosol tracers
         INTEGER          :: N_DRYDEP     ! # of dry deposited tracers
         INTEGER          :: N_DUST       ! # of dust tracers
         INTEGER          :: N_JV         ! # of J-value reactions
         INTEGER          :: N_MEMBERS    ! Max # of species per family tracer
         INTEGER          :: N_RH         ! # of RH bins for photolysis
         INTEGER          :: N_REACTIONS  ! # of chemical reactons
         INTEGER          :: N_SOA_HC     ! # of SOA HC classes and products
         INTEGER          :: N_SOA_PROD   !   (dimensions of GPROD/APROD)
         INTEGER          :: N_SPECIES    ! # of chemical species
         INTEGER          :: N_TRACERS    ! # of advected tracers
         INTEGER          :: N_WETDEP     ! # of wet deposited tracers
      END TYPE GC_DIMS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Leave these types commented out for now ... activate them later
!%%% (bmy, 9/29/10)
!      !-----------------------------------------------------------------
!      ! Column met fields
!      !-----------------------------------------------------------------
      TYPE :: GC_MET_LOCAL
! 2-D Arrays
         REAL*8, ALLOCATABLE :: ALBD(:,:)   ! Visible surface albedo [unitless]
         REAL*8, ALLOCATABLE :: AREA_M2(:,:)! Grid box surface area [cm2]
         REAL*8, ALLOCATABLE :: CLDFRC(:,:) ! Column cloud fraction [unitless]
         REAL*8, ALLOCATABLE :: FRCLND(:,:) ! Olson land fraction [unitless]
         REAL*8, ALLOCATABLE :: GWETTOP(:,:)! Top soil moisture [unitless]
         REAL*8, ALLOCATABLE :: HFLUX(:,:)  ! Sensible heat flux [W/m2]
         REAL*8, ALLOCATABLE :: LWI(:,:)    ! Land/water indices [unitless]
         REAL*8, ALLOCATABLE :: PARDR(:,:)  ! Direct  photsyn active rad [W/m2]
         REAL*8, ALLOCATABLE :: PARDF(:,:)  ! Diffuse photsyn active rad [W/m2]
         REAL*8, ALLOCATABLE :: PBLH(:,:)   ! PBL height [m]
         REAL*8, ALLOCATABLE :: PRECCON(:,:)! Conv  precip @ ground [kg/m2/s]
         REAL*8, ALLOCATABLE :: PRECTOT(:,:)! Total precip @ ground [kg/m2/s]
         REAL*8, ALLOCATABLE :: RADSWG(:,:) ! Solar radiation @ ground [W/m2]
         REAL*8, ALLOCATABLE :: SST(:,:)    ! Sea surface temperature [K]
         REAL*8, ALLOCATABLE :: SUNCOS(:,:) ! Cosine of solar zenith angle
         REAL*8, ALLOCATABLE :: TO3(:,:)    ! Total overhead O3 column [DU]
         REAL*8, ALLOCATABLE :: TROPP(:,:)  ! Tropopause pressure [hPa]
         REAL*8, ALLOCATABLE :: TS(:,:)     ! Surface temperature [K]
         REAL*8, ALLOCATABLE :: U10M(:,:)   ! E/W wind speed @ 10m height [m/s]
         REAL*8, ALLOCATABLE :: USTAR(:,:)  ! Friction velocity [m/s]
         REAL*8, ALLOCATABLE :: UVALBEDO(:,:)! UV surface albedo [unitless]
         REAL*8, ALLOCATABLE :: V10M(:,:)   ! N/S wind speed @ 10m height [m/s]
         REAL*8, ALLOCATABLE :: Z0(:,:)     ! Surface roughness height [m]
! 3-D Arrays
         REAL*8, ALLOCATABLE :: AD(:,:,:)        ! Air mass [kg]
         REAL*8, ALLOCATABLE :: AIRDENS(:,:,:)   ! Air density [kg/m3]
         REAL*8, ALLOCATABLE :: AIRVOL(:,:,:)    ! Grid box volume [m3]
         REAL*8, ALLOCATABLE :: BXHEIGHT(:,:,:)  ! Grid box height [m]
         REAL*8, ALLOCATABLE :: CLDF(:,:,:)      ! 3-D cloud fraction [unitless]
         REAL*8, ALLOCATABLE :: CMFMC(:,:,:)     ! Cloud mass flux [kg/m2/s]
         REAL*8, ALLOCATABLE :: DQIDTMST(:,:,:)  ! Ice tendency, mst proc [kg/kg/s]
         REAL*8, ALLOCATABLE :: DQLDTMST(:,:,:)  ! H2O tendency, mst proc [kg/kg/s]
         REAL*8, ALLOCATABLE :: DQVDTMST(:,:,:)  ! Vapor tendency, mst proc [kg/kg/s]
         REAL*8, ALLOCATABLE :: DTRAIN(:,:,:)    ! Detrainment flux [kg/m2/s]
         REAL*8, ALLOCATABLE :: MOISTQ(:,:,:)    ! Tendency in sp. humidity [kg/kg/s]
         REAL*8, ALLOCATABLE :: OPTD(:,:,:)      ! Visible optical depth [unitless]
         REAL*8, ALLOCATABLE :: PEDGE(:,:,:)     ! Pressure @ level edges [Pa]
         REAL*8, ALLOCATABLE :: PMID(:,:,:)      ! Pressure @ level centers [Pa]
         REAL*8, ALLOCATABLE :: DELP(:,:,:)      ! Delta-P extent  of a grid box [mb]
         REAL*8, ALLOCATABLE :: RH(:,:,:)        ! Relative humidity [unitless]
         REAL*8, ALLOCATABLE :: SPHU(:,:,:)      ! Specific humidity [kg/kg]
         REAL*8, ALLOCATABLE :: T(:,:,:)         ! Temperature [K]
         REAL*8, ALLOCATABLE :: TAUCLI(:,:,:)    ! Opt depth of ice clouds [unitless]
         REAL*8, ALLOCATABLE :: TAUCLW(:,:,:)    ! Opt depth of H2O clouds [unitless]
      END TYPE GC_MET_LOCAL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !-----------------------------------------------------------------
      ! Derived type to pass ID information to GEOS-Chem routines
      !-----------------------------------------------------------------
      TYPE :: GC_IDENT
         INTEGER            :: PET          ! # of the CPU we are on
         INTEGER            :: STDOUT_LUN   ! LUN for stdout redirect
         CHARACTER(LEN=255) :: STDOUT_FILE  ! Filename for stdout redirect
         INTEGER            :: LEV          ! Stack pointer
         CHARACTER(LEN=40)  :: I_AM(20)     ! Stack for routine names
         CHARACTER(LEN=999) :: ERRMSG       ! Error message to display
         LOGICAL            :: VERBOSE      ! Should we print out debug info?
      END TYPE GC_IDENT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Leave these types commented out for now ... activate them later
!%%% (bmy, 9/29/10)
!      !-----------------------------------------------------------------
!      ! Type for column oxidant fields used by SCHEM
!      !-----------------------------------------------------------------
!      TYPE :: SCOX_1d
!         REAL*8, POINTER    :: OH(:)        ! Monthly archived OH [molec/cm3]
!         REAL*8, POINTER    :: JVALUE(:,:)  ! Monthly archived J-values [1/s]
!         REAL*8, POINTER    :: PCO(:)       ! CO production rate [mol/mol/s]
!         REAL*8, POINTER    :: LCO(:)       ! CO loss rate [1/s]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY: 
!  24 Mar 2009 - R. Yantosca - Initial version
!  21 Apr 2009 - R. Yantosca - Renamed from "id_type_mod.f" to
!                              "gc_type_mod.f"; added type GC_OPTIONS
!  05 Jun 2009 - R. Yantosca - added LISOPOH to ID_SPEC
!                            - added MOLWT, XNUMOL to SPEC_2_TRAC
!  08 Jul 2009 - R. Yantosca - Added USE_CONVECTION, USE_PBL_MIXING,
!                              USE_WETDEP, USE_EMISSIONS flags to GC_OPTIONS
!  24 Aug 2009 - R. Yantosca - Added GC_TIME type
!  30 Oct 2008 - R. Yantosca - Added GC_GEOLOC type
!  05 Nov 2009 - R. Yantosca - Added GC_MET_LOCAL type
!  13 Apr 2010 - R. Yantosca - Added N_MEMBERS to GC_DIMS type
!  16 Apr 2010 - R. Yantosca - Added all chemical species to ID_SPEC
!  16 Apr 2010 - R. Yantosca - Added N_SPECIES, N_JV to GC_DIMS
!  23 Apr 2010 - R. Yantosca - Added GC_IDENT object
!  29 Apr 2010 - R. Yantosca - Added TO3, FRCLND to GC_MET_LOCAL object
!  30 Apr 2010 - R. Yantosca - Change character lengths in GC_IDENT
!  13 May 2010 - R. Yantosca - Updated comments  
!  17 May 2010 - R. Yantosca - Added L_COLUMN to DIMINFO
!  02 Jun 2010 - R. Yantosca - Added VERBOSE to GC_IDENT
!  02 Jun 2010 - R. Yantosca - Added TAUCLI, TAUCLW to GC_MET_LOCAL
!  02 Jun 2010 - R. Yantosca - Added DQ*DTMST fields to GC_MET_LOCAL
!  02 Jun 2010 - R. Yantosca - Moved type SCOX_1d here from schem_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE GC_TYPE_MOD
!EOC
