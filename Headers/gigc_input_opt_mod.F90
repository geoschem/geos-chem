#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_input_opt_mod
!
! !DESCRIPTION: Module GIGC\_INPUT\_OPT\_MOD contains the derived type
!  for GEOS-Chem options and logical switches.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Input_Opt_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_GIGC_Input_Opt
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Input Options 
  !=========================================================================
  TYPE, PUBLIC :: OptInput

     !%%%% Time fields %%%%
     INTEGER            :: NYMDb         ! YYYY/MM/DD @ beginning of run
     INTEGER            :: NHMSb         ! hh:mm:ss   @ beginning of run
     INTEGER            :: NYMDe         ! YYYY/MM/DD @ end of run
     INTEGER            :: NHMSe         ! hh:mm:ss   @ end of run
     INTEGER            :: TS_CHEM       ! Timestep for chemistry  [min]
     INTEGER            :: TS_DYN        ! Timestep for dynamics   [min]
     INTEGER            :: TS_EMIS       ! Timestep for emissions  [min]
     INTEGER            :: TS_CONV       ! Timestep for convection [min]

     !%%%% GEOS-Chem directories and files %%%%
     CHARACTER(LEN=255) :: DATA_DIR      ! Main DAO met field directory
     CHARACTER(LEN=255) :: DATA_DIR_1x1  ! Root data dir for 1x1 emissions
     CHARACTER(LEN=255) :: GCAP_DIR      ! Subdir for GCAP met data
     CHARACTER(LEN=255) :: GEOS_1_DIR    ! !%%% OBSOLETE %%%
     CHARACTER(LEN=255) :: GEOS_S_DIR    ! !%%% OBSOLETE 
     CHARACTER(LEN=255) :: GEOS_4_DIR    ! Subdir for GEOS-4 met data
     CHARACTER(LEN=255) :: GEOS_5_DIR    ! Subdir for GEOS-5 met data
     CHARACTER(LEN=255) :: GEOS_57_DIR   ! Subdir for GEOS-5.7.x met data
     CHARACTER(LEN=255) :: MERRA_DIR     ! Subdir for MERRA  met data
     CHARACTER(LEN=255) :: TEMP_DIR      ! Temp dir for unzipping met data
     CHARACTER(LEN=255) :: RUN_DIR       ! Run directory for GEOS-Chem
     CHARACTER(LEN=255) :: OH_DIR        ! Dir w/ mean OH files
     CHARACTER(LEN=255) :: O3PL_DIR      ! Dir w/ archived O3 P/L rate files 
     CHARACTER(LEN=255) :: TPBC_DIR      ! Dir w/ TPCORE boundary conditions
     CHARACTER(LEN=255) :: TPBC_DIR_NA   ! Dir w/ TPCORE BC's for NA nest grid
     CHARACTER(LEN=255) :: TPBC_DIR_EU   ! Dir w/ TPCORE BC's for EU nest grid
     CHARACTER(LEN=255) :: TPBC_DIR_CH   ! Dir w/ TPCORE BC's for CH nest grid
     CHARACTER(LEN=255) :: TPBC_DIR_SE   ! Dir w/ TPCORE BC's for SEAC4RS grid
     CHARACTER(LEN=255) :: IN_RST_FILE   ! Input restart file
     CHARACTER(LEN=255) :: OUT_RST_FILE  ! Input restart file

     !%%%% Nested grid fields %%%%
     INTEGER            :: I0
     INTEGER            :: J0

     !%%%% Aerosols %%%%
     LOGICAL            :: LATEQ           ! ??     
     LOGICAL            :: LCARB           ! Use carbon aerosol tracers?
     LOGICAL            :: LCRYST          ! Use Crystalline aerosols? (not implemented)
     LOGICAL            :: LCOOKE          ! Use the C
     LOGICAL            :: LDEAD           ! Use the DEAD/Zender dust algorithm?
     LOGICAL            :: LDUST           ! Use dust aerosol tracers?
     LOGICAL            :: LSULF           ! Use sulfate aerosol tracers?
     LOGICAL            :: LSOA            ! Use secondary organic aerosol tracers?
     LOGICAL            :: LSSALT          ! Use sea-salt aerosol tracers?
     LOGICAL            :: LDICARB         ! Use dicarbonyl chemistry mechanism?
     
     !%%%% Chemistry %%%%
     LOGICAL            :: LCHEM           ! Use chemistry?
     LOGICAL            :: LKPP            ! Use KPP solver instead of SMVGEAR?
     LOGICAL            :: LSCHEM          ! Use linear strat chem? (vs. no strat schem)
     LOGICAL            :: LLINOZ          ! T = Linoz O3; F = Synoz O3
     
     !%%%% Cloud convection %%%%
     LOGICAL            :: LCONV           ! Use cloud convection?
     
     !%%%% Diagnostics %%%%
     LOGICAL            :: LDBUG
     LOGICAL            :: LDIAG           ! Use diagnostics?
     LOGICAL            :: LPRT            ! Save debug output to log file?
     LOGICAL            :: LSTDRUN         ! Save out Ox.mass files for benchmarks
     
     !%%%% Deposition %%%%
     LOGICAL            :: LDRYD           ! Use dry deposition?
     LOGICAL            :: LWETD           ! Use wet deposition?
     LOGICAL            :: USE_OLSON_2001  ! Use the Olson 2001 land map & drydep info?

     !%%%% Emissions %%%%
     LOGICAL            :: LAIRNOX         ! Use aircraft NOx emissions?
     LOGICAL            :: LANTHRO         ! Turns all anthropogenic emissions on/off
     LOGICAL            :: LBBSEA          ! Use seasonal biomass burning emissions?
     LOGICAL            :: LBIONOX         ! <-- deprecated: replace w/ LBIOMASS soon
     LOGICAL            :: LBIOMASS        ! Use biomass emissions?
     LOGICAL            :: LBIOFUEL        ! Use biofuel emissions?
     LOGICAL            :: LBIOGENIC       ! Use biogenic emissions?
     LOGICAL            :: LCAC            ! Use CAC Canadian regional emissions
     LOGICAL            :: LBRAVO          ! Use BRAVO Mexican regional emissions?
     LOGICAL            :: LEDGAR          ! Use EDGAR emissions?
     LOGICAL            :: LEDGARNOx       ! Use EDGAR NOx emissions?
     LOGICAL            :: LEDGARCO        ! Use EDGAR CO emissions?
     LOGICAL            :: LEDGARSHIP      ! Use EDGAR ship emissions?
     LOGICAL            :: LEDGARSOx       ! Use EDGAR SOx emissions?
     LOGICAL            :: LEMEP           ! Use EMEP European regional emissions?
     LOGICAL            :: LEMIS           ! Use emissions in GEOS-Chem (main switch)
     LOGICAL            :: LFFNOX          ! Use anthropogenic NOx emissions
     LOGICAL            :: LFOSSIL         ! <-- deprecated: replace w/ LANTHRO soon
     LOGICAL            :: LSTREETS        ! Use David Streets SE Asian emissions?
     LOGICAL            :: LICARTT         ! Use ICARTT fix to EPA emissions?
     LOGICAL            :: LICOADSSHIP     ! Use ICOADS ship emissions inventory
     LOGICAL            :: LLIGHTNOX       ! Use lightning NOx emissions?
     LOGICAL            :: LOTDLOC         ! Use OTD-LIS local flash redistribution?
     LOGICAL            :: LMEGAN          ! Use MEGAN biogenic emissions?
     LOGICAL            :: LMEGANMONO      ! Use MEGAN monoterpenes?
     LOGICAL            :: LMONOT          ! Use old 
     LOGICAL            :: LNEI99          ! Use EPA 1999 regional emissions?
     LOGICAL            :: LNEI05          ! Use EPA 2005 regional emissions?
     LOGICAL            :: LSHIPSO2        ! Use SO2 from ship emissions?
     LOGICAL            :: LSOILNOX        ! Use soil NOx emissions
     LOGICAL            :: LFERTILIZERNOX  ! Use fertilizer NOx emissions
     LOGICAL            :: LTOMSAI         ! Scale biomass burning to TOMS AI index?
     LOGICAL            :: LWOODCO         ! <-- deprecated: replace w/ LBIOFUEL soon
     LOGICAL            :: LAVHRRLAI       ! Use AVHRR leaf-area-indices?
     LOGICAL            :: LGFED2BB        ! Use GFED2 biomass burning?
     LOGICAL            :: LGFED3BB        ! Use GFED3 biomass burning?
     LOGICAL            :: LFUTURE         ! Use future-years emission scaling (GCAP)?
     LOGICAL            :: LARCSHIP        ! Use ARCTAS ship emissions?
     LOGICAL            :: LEMEPSHIP       ! Use EMEP ship emissions?
     LOGICAL            :: LVISTAS         ! Use VISTAS NOx emissions?
     LOGICAL            :: L8DAYBB         ! Use GFED2 8-day biomass burning?
     LOGICAL            :: L3HRBB          ! Use GFED2 3-hr biomass burning?
     LOGICAL            :: LSYNOPBB        ! Use GFED2 synoptic biomass burning
     LOGICAL            :: LDAYBB3         ! Use GFED3 daily biomass burning?
     LOGICAL            :: L3HRBB3         ! Use GFED3 3-hr biomass burning?
     LOGICAL            :: LMODISLAI       ! MODIS LAI (mpb, 2009)
     LOGICAL            :: LPECCA          ! PCEEA BVOC emission model (mpb,2009)
     LOGICAL            :: LRETRO          ! RETRO anthropogenic emissions (wfr, 3/8/11)
     LOGICAL            :: LHIST           ! Use historical emissions?
     LOGICAL            :: LWARWICK_VSLS   ! Use Warwick_s3 R(Br)n? (jpp 8/23/07)
     LOGICAL            :: LSSABr2         ! Use sea salt Br2 emiss? (jpp, 3/28/2010)
     LOGICAL            :: LFIX_PBL_BrO    ! Run 1ppt MBL BRO Sim.? (jpp 5/10/10)
     REAL*8             :: Br_SCALING  
     REAL*8             :: ISOP_SCALING
     REAL*8             :: NOx_SCALING

     !%%%% Transport and strat BC's %%%%
     LOGICAL            :: LFILL           ! Fill negative values in TPCORE?
     LOGICAL            :: LMFCT           ! Turns TPCORE MFCT option on/off
     LOGICAL            :: LTRAN           ! Turns advection on/off
     LOGICAL            :: LTPFV           ! Are we using the GEOS-4 
     LOGICAL            :: LWINDO          ! Use nested-grid simulation?
     LOGICAL            :: LWINDO2x25      ! Use nested-grid BC's @ 2 x 2.5 resolution?
     LOGICAL            :: LWINDO_NA       ! Use 
     LOGICAL            :: LWINDO_EU
     LOGICAL            :: LWINDO_CH
     LOGICAL            :: LWINDO_CU
     LOGICAL            :: LWINDO_SE
     LOGICAL            :: LUPBD           ! Use stratospheric O3, NOY bdry conditions

     
     !%%%% Met fields %%%%
     LOGICAL            :: LUNZIP          ! Unzip met fields on-the-fly?
     LOGICAL            :: LWAIT           ! Wait for met fields to be unzipped?
     
     !%%%% PBL mixing %%%%
     LOGICAL            :: LTURB           ! Use PBL mixing?
     LOGICAL            :: LNLPBL          ! Use non-local PBL scheme instead of default?
     LOGICAL            :: LARPBLH         ! --> pblh_ar in vdiff_mod.f
     LOGICAL            :: LDEPBCK         ! --> drydep_back_cons in vdiff_mod.f
     
     !%%%% Restart files %%%%
     LOGICAL            :: LSVGLB                   ! Save tracer restart file?
     LOGICAL            :: LSVCSPEC                 ! Save CSPEC restart file?
     
     !%%%% Tagged simulations %%%%
     LOGICAL            :: LSPLIT                   ! Tagged tracers?
     
     !%%%% Variable Tropopause %%%%
     LOGICAL            :: LVARTROP                 ! Dynamic tropopause?
     
     !%%%% Dynamic ocean Hg model %%%%
     LOGICAL            :: LDYNOCEAN                ! Dynamic ocean Hg model?
     LOGICAL            :: LPREINDHG                ! Preindustrial Hg sim
     LOGICAL            :: LGTMM                    ! GTMM soil model
     
     !%%%% For the CH4 offline simulation only %%%%
     LOGICAL            :: LGAO                     ! Use gas & oil emissions?
     LOGICAL            :: LCOL                     ! Use coal emissions?
     LOGICAL            :: LLIV                     ! Use livestock emissions?
     LOGICAL            :: LWAST                    ! Use waste emissions?
     LOGICAL            :: LRICE                    ! Use rice emissions?
     LOGICAL            :: LOTANT                   ! Use other anthro emiss?
     LOGICAL            :: LWETL                    ! Use wetland emissions?
     LOGICAL            :: LSOABS                   ! Use soil absorption?
     LOGICAL            :: LOTNAT                   ! Use other natural emiss?
     LOGICAL            :: LBFCH4                   ! Use CH4 biofuel emissions?
     LOGICAL            :: LBMCH4                   ! Use CH4 biomass emissions?
     LOGICAL            :: LCH4BUD                  ! Use computing CH4 budget
     
     !%%%% For the CO2 offline simulation only %%%%
     LOGICAL            :: LGENFF      
     LOGICAL            :: LANNFF      
     LOGICAL            :: LMONFF      
     LOGICAL            :: LSEASBB
     LOGICAL            :: LBIONETORIG  
     LOGICAL            :: LBIONETCLIM  
     LOGICAL            :: LBIODAILY
     LOGICAL            :: LBIODIURNAL    
     LOGICAL            :: LOCEAN     
     LOGICAL            :: LFFBKGRD
     LOGICAL            :: LBIOSPHTAG     
     LOGICAL            :: LFOSSILTAG     
     LOGICAL            :: LOCN1997
     LOGICAL            :: LOCN2009ANN
     LOGICAL            :: LOCN2009MON
     LOGICAL            :: LSHIPEDG
     LOGICAL            :: LSHIPICO
     LOGICAL            :: LSHIPSCALE
     LOGICAL            :: LSHIPTAG
     LOGICAL            :: LPLANE
     LOGICAL            :: LPLANESCALE
     LOGICAL            :: LPLANETAG
     LOGICAL            :: LCHEMCO2
     
     !%%%% HDF5 output %%%%
     LOGICAL            :: LND50_HDF                ! Save ND50  in HDF5?
     LOGICAL            :: LND51_HDF                ! Save ND51  in HDF5?
     LOGICAL            :: LND51b_HDF               ! Save ND51b c in HDF5?
    
     !%%%% For interface with GCM's in MPI environments
     LOGICAL            :: DO_DIAG_WRITE = .TRUE.   ! TURN ON/OFF DIAG WRITING
     
     !%%%% Flags that denote simulation types %%%%
     LOGICAL            :: ITS_A_RnPbBe_SIM         ! Rn-Pb-Be  simulation?
     LOGICAL            :: ITS_A_CH3I_SIM           ! CH3I      simulation?
     LOGICAL            :: ITS_A_FULLCHEM_SIM       ! Fullchem  simulation?
     LOGICAL            :: ITS_A_HCN_SIM            ! HCN       simulation?
     LOGICAL            :: ITS_A_TAGOX_SIM          ! Tagged Ox simulation?
     LOGICAL            :: ITS_A_TAGCO_SIM          ! Tagged CO simulation?
     LOGICAL            :: ITS_A_C2H6_SIM           ! C2H6      simulation?
     LOGICAL            :: ITS_A_CH4_SIM            ! CH4       simulation?
     LOGICAL            :: ITS_AN_AEROSOL_SIM       ! Aerosol   simulation?
     LOGICAL            :: ITS_A_MERCURY_SIM        ! Mercury   simulation?
     LOGICAL            :: ITS_A_CO2_SIM            ! CO2       simulation
     LOGICAL            :: ITS_A_H2HD_SIM           ! H2/HD     simulation?
     LOGICAL            :: ITS_NOT_COPARAM_OR_CH4   ! Not COPARAM or CH4?

  END TYPE OptInput
!
! !REMARKS:
!  This will eventually replace the switches in logical_mod.F.  
!
! !REVISION HISTORY:
!  01 Nov 2012 - R. Yantosca - Initial version, based on logical_mod.F
!                              newer Olson 2001 land map & drydep inputs
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_gigc_input_opt
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_INPUT\_OPT intializes all GEOS-Chem
!  options carried in Input Options derived type object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  May need to add better error checking.
! 
! !REVISION HISTORY: 
!  01 Nov 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Assume success
    RC                               = GIGC_SUCCESS

    ! Initialize time fields
    Input_Opt%NYMDb                  = 0
    Input_Opt%NHMSb                  = 0
    Input_Opt%NYMDe                  = 0
    Input_Opt%NHMSe                  = 0
    Input_Opt%TS_CHEM                = 0
    Input_Opt%TS_DYN                 = 0
    Input_Opt%TS_EMIS                = 0
    Input_Opt%TS_CONV                = 0

    ! Initialize directories         = ''
    Input_Opt%DATA_DIR               = ''
    Input_Opt%DATA_DIR_1x1           = ''
    Input_Opt%GCAP_DIR               = ''
    Input_Opt%GEOS_1_DIR             = ''
    Input_Opt%GEOS_S_DIR             = ''
    Input_Opt%GEOS_4_DIR             = ''
    Input_Opt%GEOS_5_DIR             = ''
    Input_Opt%GEOS_57_DIR            = ''
    Input_Opt%MERRA_DIR              = ''
    Input_Opt%TEMP_DIR               = ''
    Input_Opt%RUN_DIR                = ''
    Input_Opt%OH_DIR                 = ''
    Input_Opt%O3PL_DIR               = ''
    Input_Opt%TPBC_DIR               = ''
    Input_Opt%TPBC_DIR_NA            = ''
    Input_Opt%TPBC_DIR_EU            = ''
    Input_Opt%TPBC_DIR_CH            = ''
    Input_Opt%TPBC_DIR_SE            = ''

    ! Initialize nested grid fields
    Input_Opt%I0                     = 0
    Input_Opt%J0                     = 0
    
    ! Initialize scaling fields
    Input_Opt%Br_SCALING             = 0d0 
    Input_Opt%ISOP_SCALING           = 0d0
    Input_Opt%NOx_SCALING            = 0d0

    ! Initialize logical fields
    Input_Opt%LATEQ                  = .FALSE.
    Input_Opt%LCARB                  = .FALSE.
    Input_Opt%LCRYST                 = .FALSE.
    Input_Opt%LCOOKE                 = .FALSE.
    Input_Opt%LDEAD                  = .FALSE.
    Input_Opt%LDUST                  = .FALSE.
    Input_Opt%LSULF                  = .FALSE.
    Input_Opt%LSOA                   = .FALSE.
    Input_Opt%LSSALT                 = .FALSE.
    Input_Opt%LDICARB                = .FALSE.
    Input_Opt%LCHEM                  = .FALSE.
    Input_Opt%LKPP                   = .FALSE.
    Input_Opt%LSCHEM                 = .FALSE.
    Input_Opt%LLINOZ                 = .FALSE.
    Input_Opt%LCONV                  = .FALSE.
    Input_Opt%LDBUG                  = .FALSE.
    Input_Opt%LDIAG                  = .FALSE.
    Input_Opt%LPRT                   = .FALSE.
    Input_Opt%LSTDRUN                = .FALSE.
    Input_Opt%LDRYD                  = .FALSE.
    Input_Opt%LWETD                  = .FALSE.
    Input_Opt%USE_OLSON_2001         = .FALSE.
    Input_Opt%LAIRNOX                = .FALSE.
    Input_Opt%LANTHRO                = .FALSE.
    Input_Opt%LBBSEA                 = .FALSE.
    Input_Opt%LBIONOX                = .FALSE.
    Input_Opt%LBIOMASS               = .FALSE.
    Input_Opt%LBIOFUEL               = .FALSE.
    Input_Opt%LBIOGENIC              = .FALSE.
    Input_Opt%LCAC                   = .FALSE.
    Input_Opt%LBRAVO                 = .FALSE.
    Input_Opt%LEDGAR                 = .FALSE.
    Input_Opt%LEDGARNOx              = .FALSE.
    Input_Opt%LEDGARCO               = .FALSE.
    Input_Opt%LEDGARSHIP             = .FALSE.
    Input_Opt%LEDGARSOx              = .FALSE.
    Input_Opt%LEMEP                  = .FALSE.
    Input_Opt%LEMIS                  = .FALSE.
    Input_Opt%LFFNOX                 = .FALSE.
    Input_Opt%LFOSSIL                = .FALSE.
    Input_Opt%LSTREETS               = .FALSE.
    Input_Opt%LICARTT                = .FALSE.
    Input_Opt%LICOADSSHIP            = .FALSE.
    Input_Opt%LLIGHTNOX              = .FALSE.
    Input_Opt%LOTDLOC                = .FALSE.
    Input_Opt%LMEGAN                 = .FALSE.
    Input_Opt%LMEGANMONO             = .FALSE.
    Input_Opt%LMONOT                 = .FALSE.
    Input_Opt%LNEI99                 = .FALSE.
    Input_Opt%LNEI05                 = .FALSE.
    Input_Opt%LSHIPSO2               = .FALSE.
    Input_Opt%LSOILNOX               = .FALSE.
    Input_Opt%LFERTILIZERNOX         = .FALSE. 
    Input_Opt%LTOMSAI                = .FALSE.
    Input_Opt%LWOODCO                = .FALSE.
    Input_Opt%LAVHRRLAI              = .FALSE.
    Input_Opt%LGFED2BB               = .FALSE.
    Input_Opt%LGFED3BB               = .FALSE.
    Input_Opt%LFUTURE                = .FALSE.
    Input_Opt%LARCSHIP               = .FALSE.
    Input_Opt%LEMEPSHIP              = .FALSE.
    Input_Opt%LVISTAS                = .FALSE.
    Input_Opt%L8DAYBB                = .FALSE.
    Input_Opt%L3HRBB                 = .FALSE.
    Input_Opt%LSYNOPBB               = .FALSE.
    Input_Opt%LDAYBB3                = .FALSE.
    Input_Opt%L3HRBB3                = .FALSE.
    Input_Opt%LMODISLAI              = .FALSE.
    Input_Opt%LPECCA                 = .FALSE.
    Input_Opt%LRETRO                 = .FALSE.
    Input_Opt%LHIST                  = .FALSE.
    Input_Opt%LWARWICK_VSLS          = .FALSE.
    Input_Opt%LSSABr2                = .FALSE.
    Input_Opt%LFIX_PBL_BrO           = .FALSE.
    Input_Opt%LFILL                  = .FALSE.
    Input_Opt%LMFCT                  = .FALSE.
    Input_Opt%LTRAN                  = .FALSE.
    Input_Opt%LTPFV                  = .FALSE.
    Input_Opt%LWINDO                 = .FALSE.
    Input_Opt%LWINDO2x25             = .FALSE.
    Input_Opt%LWINDO_NA              = .FALSE.
    Input_Opt%LWINDO_EU              = .FALSE.
    Input_Opt%LWINDO_CH              = .FALSE.
    Input_Opt%LWINDO_CU              = .FALSE.
    Input_Opt%LWINDO_SE              = .FALSE.
    Input_Opt%LUPBD                  = .FALSE.
    Input_Opt%LUNZIP                 = .FALSE.
    Input_Opt%LWAIT                  = .FALSE.
    Input_Opt%LTURB                  = .FALSE.
    Input_Opt%LNLPBL                 = .FALSE.
    Input_Opt%LARPBLH                = .FALSE.
    Input_Opt%LDEPBCK                = .FALSE.
    Input_Opt%LSVGLB                 = .FALSE.
    Input_Opt%LSVCSPEC               = .FALSE.
    Input_Opt%LSPLIT                 = .FALSE.
    Input_Opt%LVARTROP               = .FALSE.
    Input_Opt%LDYNOCEAN              = .FALSE.
    Input_Opt%LPREINDHG              = .FALSE.
    Input_Opt%LGTMM                  = .FALSE.
    Input_Opt%LGAO                   = .FALSE.
    Input_Opt%LCOL                   = .FALSE.
    Input_Opt%LLIV                   = .FALSE.
    Input_Opt%LWAST                  = .FALSE.
    Input_Opt%LRICE                  = .FALSE.
    Input_Opt%LOTANT                 = .FALSE.
    Input_Opt%LWETL                  = .FALSE.
    Input_Opt%LSOABS                 = .FALSE.
    Input_Opt%LOTNAT                 = .FALSE.
    Input_Opt%LBFCH4                 = .FALSE.
    Input_Opt%LBMCH4                 = .FALSE.
    Input_Opt%LCH4BUD                = .FALSE.
    Input_Opt%LGENFF                 = .FALSE.
    Input_Opt%LANNFF                 = .FALSE.
    Input_Opt%LMONFF                 = .FALSE.
    Input_Opt%LSEASBB                = .FALSE.
    Input_Opt%LBIONETORIG            = .FALSE.
    Input_Opt%LBIONETCLIM            = .FALSE.
    Input_Opt%LBIODAILY              = .FALSE.
    Input_Opt%LBIODIURNAL            = .FALSE.
    Input_Opt%LOCEAN                 = .FALSE.
    Input_Opt%LFFBKGRD               = .FALSE.
    Input_Opt%LBIOSPHTAG             = .FALSE.
    Input_Opt%LFOSSILTAG             = .FALSE.
    Input_Opt%LOCN1997               = .FALSE.
    Input_Opt%LOCN2009ANN            = .FALSE.
    Input_Opt%LOCN2009MON            = .FALSE.
    Input_Opt%LSHIPEDG               = .FALSE.
    Input_Opt%LSHIPICO               = .FALSE.
    Input_Opt%LSHIPSCALE             = .FALSE.
    Input_Opt%LSHIPTAG               = .FALSE.
    Input_Opt%LPLANE                 = .FALSE.
    Input_Opt%LPLANESCALE            = .FALSE.
    Input_Opt%LPLANETAG              = .FALSE.
    Input_Opt%LCHEMCO2               = .FALSE.
    Input_Opt%LND50_HDF              = .FALSE.
    Input_Opt%LND51_HDF              = .FALSE.
    Input_Opt%LND51b_HDF             = .FALSE.
    Input_Opt%DO_DIAG_WRITE          = .FALSE.
    Input_Opt%ITS_A_RnPbBe_SIM       = .FALSE.
    Input_Opt%ITS_A_CH3I_SIM         = .FALSE.
    Input_Opt%ITS_A_FULLCHEM_SIM     = .FALSE.
    Input_Opt%ITS_A_HCN_SIM          = .FALSE.
    Input_Opt%ITS_A_TAGOX_SIM        = .FALSE.
    Input_Opt%ITS_A_TAGCO_SIM        = .FALSE.
    Input_Opt%ITS_A_C2H6_SIM         = .FALSE.
    Input_Opt%ITS_A_CH4_SIM          = .FALSE.
    Input_Opt%ITS_AN_AEROSOL_SIM     = .FALSE.
    Input_Opt%ITS_A_MERCURY_SIM      = .FALSE.
    Input_Opt%ITS_A_CO2_SIM          = .FALSE.
    Input_Opt%ITS_A_H2HD_SIM         = .FALSE.
    Input_Opt%ITS_NOT_COPARAM_OR_CH4 = .FALSE.

  END SUBROUTINE Set_GIGC_Input_Opt
!EOC
END MODULE GIGC_Input_Opt_Mod
#endif
