!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gckpp_HetRates
!
! !DESCRIPTION: FlexChem module for heterogeneous chemistry, via KPP.
!\\
!\\
! !INTERFACE:
!
MODULE GCKPP_HETRATES
!
! !USES:
!
  USE CMN_FJX_MOD,        ONLY : NDUST
  USE CMN_FJX_MOD,        ONLY : NAER
  USE ERROR_MOD,          ONLY : ERROR_STOP
  USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
  USE ERROR_MOD,          ONLY : IS_SAFE_DIV, SAFE_DIV
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global,       ONLY : HET
  USE State_Chm_Mod,      ONLY : ChmState
  USE State_Chm_Mod,      ONLY : Ind_
  USE State_Met_Mod,      ONLY : MetState
  USE Input_Opt_Mod,      ONLY : OptInput
  USE PhysConstants,      ONLY : AVO, RGASLATM, CONSVAP, RSTARG, PI
  USE Precision_Mod,      ONLY : fp

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SET_HET
!
! !PRIVATE MEMBER FUNCTIONS:
!
  ! These functions are used for all mechanisms
  PRIVATE :: HetNO3
  PRIVATE :: HetNO2
  PRIVATE :: HetHO2
  PRIVATE :: HetGLYX
  PRIVATE :: HetMGLY
  PRIVATE :: HetIEPOX
  PRIVATE :: HetIMAE
  PRIVATE :: HetLVOC
  PRIVATE :: HetISOPND
  PRIVATE :: HetISOPNB
  PRIVATE :: HetMACRN
  PRIVATE :: HetMVKN
  PRIVATE :: HetR4N2
  PRIVATE :: HetISN1OG
  PRIVATE :: HetDHDN
  PRIVATE :: HetMONITS
  PRIVATE :: HetMONITU
  PRIVATE :: HetHONIT
  PRIVATE :: HetIONITA
  PRIVATE :: HetMONITA
  PRIVATE :: HetHBr
  PRIVATE :: HetN2O5
  PRIVATE :: N2O5
  PRIVATE :: HO2

  ! New iodine heterogeneous chemistry
  PRIVATE :: HETIUptake
  PRIVATE :: HETIXCycleSSA

  ! These are the new Br/Cl functions from J. Schmidt
  PRIVATE :: HETBrNO3_JS
  PRIVATE :: HETClNO3_JS
  PRIVATE :: HETHOBr_HBr_JS
  PRIVATE :: HETHOBr_HCl_JS
  PRIVATE :: HETClNO3_HBr_JS
  PRIVATE :: HETO3_HBr_JS
  PRIVATE :: HETHOBr_SS_JS
  PRIVATE :: HETClNO3_SS_JS
  PRIVATE :: HETO3_SS_JS
  PRIVATE :: HETHXUptake_JS
  PRIVATE :: HETN2O5_SS

  ! New subroutines required by the JS functions
  PRIVATE :: Gamma_ClNO3_Br
  PRIVATE :: Gamma_O3_Br
  PRIVATE :: Gamma_HOBr_X
  PRIVATE :: Gamma_HX_Uptake
  PRIVATE :: Coth
  PRIVATE :: ReactoDiff_Corr
  PRIVATE :: Gamma_HOBr_CLD     !qjc
  PRIVATE :: Gamma_HOBr_AER     !qjc

  ! These are formerly strat-only reactions extended to take place in the
  ! troposphere on sulfate aerosol
  PRIVATE :: HETClNO3_HCl
  PRIVATE :: HETHOCl_HBr
  PRIVATE :: HETHOCl_HCl
  PRIVATE :: HETBrNO3_HCl

  ! These are subfunctions to calculate rates on/in clouds and SSA
  PRIVATE :: CLD_PARAMS
  PRIVATE :: GET_HALIDE_CLDConc
  Private :: Get_Halide_SSAConc
  PRIVATE :: COMPUTE_L2G_LOCAL
  PRIVATE :: CLD1K_XNO3
  PRIVATE :: EPOXUPTK
  PRIVATE :: FCRO2HO2
  PRIVATE :: FYHORO
  PRIVATE :: FYRNO3
  PRIVATE :: ARSL1K
  PRIVATE :: kIIR1Ltd
  PRIVATE :: kIIR1R2L
!
! !PRIVATE DATA MEMBERS:
!
  ! Scalars
  INTEGER  :: NAEROTYPE
  LOGICAL  :: NATSURFACE,   PSCBOX,    STRATBOX
  REAL(fp) :: TEMPK,        RELHUM,    SUNCOS,  SPC_SO4
  REAL(fp) :: SPC_NIT,      GAMMA_HO2, XTEMP,   XDENA
  REAL(fp) :: QLIQ,         QICE,      SPC_SALA
  REAL(fp) :: H_PLUS,       MSO4,      MNO3,    MHSO4
  REAL(fp) :: MW_HO2,       MW_NO2,    MW_NO3
  REAL(fp) :: MW_N2O5,      MW_GLYX,   MW_MGLY
  REAL(fp) :: MW_IEPOXA,    MW_IEPOXB, MW_IEPOXD
  REAL(fp) :: MW_HMML,      MW_LVOC,   MW_ICHE
  REAL(fp) :: MW_ITHN,      MW_ITCN,   MW_IDN
  REAL(fp) :: MW_MVKN,      MW_MCRHN,  MW_MCRHNB
  REAL(fp) :: MW_R4N2,      MW_INPB,   MW_INPD
  REAL(fp) :: MW_IHN1,      MW_IHN2,   MW_IHN3
  REAL(fp) :: MW_IONITA,    MW_MONITA, MW_IHN4
  REAL(fp) :: MW_MONITS,    MW_MONITU, MW_HONIT
  REAL(fp) :: MW_HOBr,      MW_HBr,    MW_ClNO3
  REAL(fp) :: MW_HOCl,      MW_HI,     MW_HOI
  REAL(fp) :: MW_I2O2,      MW_I2O3,   MW_I2O4
  REAL(fp) :: MW_IONO,      MW_IONO2,  MW_HCl
  REAL(fp) :: MW_O3,        MW_BrNO3,  MW_PYAC
  REAL(fp) :: H_K0_O3,      H_CR_O3,   H_O3_T
  REAL(fp) :: H_K0_HOBr,    H_CR_HOBr, H_HOBr_T
  REAL(fp) :: H_K0_HBr,     H_CR_HBr
  REAL(fp) :: H_K0_HCl,     H_CR_HCl
  REAL(fp) :: HSO3conc_Cld, SO3conc_Cld, fupdateHOBr
  REAL(fp) :: OMOC_POA, OMOC_OPOA
  ! Arrays
  REAL(fp) :: XAREA(25), XRADI(25), XVOL(25), XH2O(25)
  REAL(fp) :: KHETI_SLA(11)

!$OMP THREADPRIVATE( NAEROTYPE,        NATSURFACE, PSCBOX,   STRATBOX )
!$OMP THREADPRIVATE( TEMPK,        RELHUM,     SPC_NIT,  SPC_SO4  )
!$OMP THREADPRIVATE( GAMMA_HO2,    XTEMP,      XDENA,    QLIQ     )
!$OMP THREADPRIVATE( QICE,         KHETI_SLA,  SUNCOS             )
!$OMP THREADPRIVATE( XAREA,        XRADI,      XVOL,     XH2O     )
!$OMP THREADPRIVATE( OMOC_POA,     OMOC_OPOA,  SPC_SALA           )
!$OMP THREADPRIVATE( H_PLUS,       MSO4,       MNO3,     MHSO4    )
!$OMP THREADPRIVATE( HSO3conc_Cld, SO3conc_Cld, fupdateHOBr       )
!
! !DEFINED PARAMETERS:
!
  REAL(fp), PARAMETER :: HetMinLife = 1.e-3_fp

  ! Critical RH for uptake of GLYX, MGLYX, and GLYC:
  REAL(fp), PARAMETER :: CRITRH = 35.0e+0_fp

  ! Effective Henry's Law constant of IEPOX for reactive
  ! uptake to aqueous aerosols (M/atm)
  !REAL(fp), PARAMETER :: HSTAR_EPOX = 5.0e+6_fp ! Prior to 3/2/18
  REAL(fp), PARAMETER :: HSTAR_EPOX = 1.7e+7_fp


  ! Conversion factor from atm to bar
  REAL(fp), PARAMETER :: con_atm_bar = 1.0e+0_fp/1.01325e+0_fp

  ! Universal gas consatant [bar/(mol/kg)/K]
  REAL(fp), PARAMETER :: con_R = RStarG*1.0e-2_fp
!
! !REMARKS:
!  Need
!  - TOTAREA (previously used for archiving N2O5 hydrolysis in the planeflight
!             diagnostic only)
!  - Air NUM. DENSITY
!  - TEMPERATURE
!  - Aerosol Surface Area
!  - Aerosol Type
!  - Gamma (XSTKCF; sticking factor)
!  - ARR
!  - Species num density (mcl cm-3)
!  - Continental PBL or no?
!  - In stratosphere or no?
!  - Reaction index (e.g. NK1HBr, NK2HBr)
!
!  According to S. Eastham, we should also include
!  cloud and ice area explicitly, in addition to
!  aerosol area
!
!  C.D. Holmes: Cloud uptake of HBr, HOBr, BrNO3, and Cl equivalents should be handled by CloudHet
!  so that they account for entrainment limits.
!
! !REFERENCES:
!  Eastham et al., Development and evaluation of the unified tropospheric-
!    stratospheric chemistry extension (UCX) for the global chemistry-transport
!    model GEOS-Chem, Atmos. Env., doi:10.1016/j.atmosenv.2014.02.001, 2014.
!  Fisher et al, Organic nitrate chemistry and its implications for nitrogen
!    budgets in an isoprene- and monoterpene-rich atmosphere: constraints from
!    aircraft (SEAC4RS) and ground-based (SOAS) observations in the Southeast
!    US. Atmos. Chem. Phys., 16, 2961-2990, 2016.
!  Marais et al., Aqueous-phase mechanism for secondary organic aerosol
!    formation from isoprene: application to the southeast United States and
!    co-benefit of SO2 emission controls, Atmos. Chem. Phys., 16, 1603-1618,
!    doi:10.5194/acp-16-1603-2016, 2016.
!  Parrella et al, Tropospheric bromine chemistry: implications for present and
!    pre-industrial ozone and mercury, Atmos. Chem. Phys., 12, 6,723-6,740,
!    doi:10.5194/acp-12-6723-2012, 2012.
!  Schmidt, J., et al., “Modelling the observed tropospheric BrO background:
!    Importance of multiphase chemistry & implications for ozone, OH, &
!    mercury”, J Geophys. Res-Atmos., 121, 024229,
!    https://doi.org/10.1002/2015JD024229, 2016.
!  Sherwen, T., et al., Global impacts of tropospheric halogens (Cl, Br, I) on
!    oxidants and composition in GEOS-Chem, Atmos. Chem. Phys., 16, 12239-12271,
!    https://doi.org/10.5194/acp-16-12239-2016, 2016.
!
! !REVISION HISTORY:
!  14 Dec 2015 - M. Long     - Initial version
!  29 Jan 2016 - M. Sulprizio- Update to include heterogeneous chemistry for
!                              UCX mechanism
!  29 Mar 2016 - R. Yantosca - NOTE: SPC_HBR and SPC_HOBR are defined
!                              for trop-only mechanisms
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!  29 Mar 2016 - R. Yantosca - Moved all the UCX-based functions to the
!                              end of the module, for clarity
!  01 Apr 2016 - R. Yantosca - Remove many global variables that can be
!                              declared locally from the THREADPRIVATEs
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  14 Jun 2016 - M. Sulprizio- Replace Spc_GetIndx with Ind_
!  15 Jun 2017 - M. Sulprizio- Add heterogeneous chemistry for isoprene SOA from
!                              E. Marais (Marais et al., 2016)
!  14 Jul 2017 - M. Sulprizio- Add heterogeneous chemistry for monoterpenes from
!                              J. Fisher (Fisher et al., 2017)
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  15 Nov 2017 - M. Sulprizio- Add modifications for HOBr + S(IV) based on work
!                              by Qianjie Chen
!  02 Mar 2018 - M. Sulprizio- Update HSTAR_EPOX following recommendation from
!                              E. Marais to address SOAIE being a factor of 2
!                              lower in v11-02d than in Marais et al. [2016]
!  17 Oct 2018 - C.D. Holmes - Added cloud heterogeneous chemistry (CloudHet);
!                              Gamma updates for NOx species
!  13 Dec 2018 - E. McDuffie - Report gammaN2O5 as a State Chem parameter
!  04 Nov 2019 - C.D. Holmes - Bug fixes for gammaN2O5 (RH dependence)
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
! !IROUTINE: Set_Het
!
! !DESCRIPTION: Main heterogenous chemistry driver routine.  Sets up the
!  vector of heterogeneous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE SET_HET( I, J, L, Input_Opt, State_Chm, State_Met )
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, level indices
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  01 Apr 2016 - R. Yantosca - Define many variables locally that don't
!                              need to be in the THREADPRIVATE statements
!  01 Apr 2016 - R. Yantosca - Remove KII_KI; we now declare that locally
!  31 May 2016 - E. Lundgren - Replace Input_Opt%XNUMOL with emMW_g from species
!                              database (emitted species g/mol)
!  26 Jul 2017 - M. Sulprizio- Remove hardcoded molecular weights from calls to
!                              Het* functions and use MW from species database
!                              instead
!  03 Jan 2018 - M. Sulprizio- Remove SCF argument. It was apparently added for
!                              diagnostic purposes and is no longer used. Also
!                              rename IO,SM,SC to Input_Opt,State_Met,State_Chm
!                              for consistency with other GEOS-Chem routines.
!  27 Feb 2018 - M. Sulprizio- Obtain Henry's law parameters from species
!                              database instead of hardcoding in halogens code
!  17 Oct 2018 - C.D. Holmes - Added cloud heterogeneous chemistry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL  :: SAFEDIV
      INTEGER  :: IND
      REAL(fp) :: ADJUSTEDRATE, CONSEXP,  DUMMY,    HBr_RTEMP
      REAL(fp) :: HOBr_RTEMP,   QICE,     QLIQ,     SPC_BrNO3
      REAL(fp) :: SPC_ClNO3,    SPC_H2O,  SPC_HBr,  SPC_HCl
      REAL(fp) :: SPC_HOBr,     SPC_HOCl, SPC_N2O5, VPRESH2O
      LOGICAL, SAVE :: FIRST = .TRUE.

      ! New treatment for educt removal
      Real(fp),Pointer :: spcVec(:)
      Real(fp)         :: kITemp, kIITemp
      Real(fp)         :: HetTemp(3)      !new temp array to hold gamma results

      ! Cloud parameters
      Real(fp)         :: rLiq, ALiq, VLiq, CLDFr
      Real(fp)         :: rIce, AIce, VIce

      ! Volume of air (cm3)
      Real(fp)         :: VAir

      ! New bromine/chlorine chemistry
      Logical, Parameter :: fixedSaltBr=.True.
      Logical            :: useSaltBr
      Real(fp)           :: hConc_Sul
      Real(fp)           :: hConc_LCl
      Real(fp)           :: hConc_ICl
      Real(fp)           :: hConc_SSA
      Real(fp)           :: hConc_SSC
      Real(fp)           :: brConc_Base
      Real(fp)           :: brConc_Cld, clConc_Cld
      Real(fp)           :: brConc_SSA, brConc_SSC
      Real(fp)           :: pHCloud
      Real(fp)           :: SSAlk(2)

      !====================================================================
      ! SET_HET begins here!
      !====================================================================

      ! Zero scalars and arrays
      ADJUSTEDRATE  = 0.0_fp
      CONSEXP       = 0.0_fp
      DUMMY         = 0.0_fp
      HBr_RTEMP     = 0.0_fp
      HOBr_RTEMP    = 0.0_fp
      KHETI_SLA     = 0.0_fp
      NAEROTYPE         = State_Chm%nAeroType
      QICE          = 0.0_fp
      QLIQ          = 0.0_fp
      SPC_BrNO3     = 0.0_fp
      SPC_ClNO3     = 0.0_fp
      SPC_H2O       = 0.0_fp
      SPC_HBr       = 0.0_fp
      SPC_HCl       = 0.0_fp
      SPC_HOBr      = 0.0_fp
      SPC_HOCl      = 0.0_fp
      SPC_N2O5      = 0.0_fp
      VPRESH2O      = 0.0_fp
      HetTemp       = 0.0_fp
      OMOC_POA      = State_Chm%OMOC_POA(I,J)
      OMOC_OPOA     = State_Chm%OMOC_OPOA(I,J)

      ! Initialize logicals
      SAFEDIV       = .FALSE.
      PSCBOX        = .FALSE.
      STRATBOX      = .FALSE.
      NATSURFACE    = .FALSE.

      ! KHETI_SLA = sticking coefficients for PSC reactions on SLA
      IF ( Input_Opt%LUCX ) THEN
         KHETI_SLA  = State_Chm%KHETI_SLA(I,J,L,:)
      ENDIF

      ! Point to the chemical species array [molec/cm3]
      spcVec          => State_Chm%Species(I,J,L,:)

      !--------------------------------------------------------------------
      ! Calculate RH [%]
      ! Not clear why this calc is slightly different than State_Met%RH
      !--------------------------------------------------------------------
      RELHUM        = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)
      CONSEXP       = 17.2693882e+0_fp * (State_Met%T(I,J,L) - 273.16e+0_fp) /&
                      (State_Met%T(I,J,L) - 35.86e+0_fp)
      VPRESH2O      = CONSVAP * EXP(CONSEXP) / State_Met%T(I,J,L)
      RELHUM        = RELHUM / VPRESH2O
      RELHUM        = RELHUM * 100e+0_fp

      !--------------------------------------------------------------------
      ! Get species molecular weights [g/mol]
      !--------------------------------------------------------------------
      IF ( FIRST) THEN
         ! Hardcode HO2 for now
         ! MW_g is not defined for HO2 in the species database but model
         ! output changes when it is added there (mps, 7/26/17)
         MW_HO2    = 33.0_fp

         IND = Ind_( 'NO2' )
         IF ( IND > 0 ) MW_NO2    = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'NO3' )
         IF ( IND > 0 ) MW_NO3    = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'N2O5' )
         IF ( IND > 0 ) MW_N2O5   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'GLYX' )
         IF ( IND > 0 ) MW_GLYX   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MGLY' )
         IF ( IND > 0 ) MW_MGLY   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IEPOXA' )
         IF ( IND > 0 ) MW_IEPOXA = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IEPOXB' )
         IF ( IND > 0 ) MW_IEPOXB = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IEPOXD' )
         IF ( IND > 0 ) MW_IEPOXD = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HMML' )
         IF ( IND > 0 ) MW_HMML   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'PYAC' )
         IF ( IND > 0 ) MW_PYAC   = State_Chm%SpcData(IND)%Info%MW_g
	 
         IND = Ind_( 'LVOC' )
         IF ( IND > 0 ) MW_LVOC   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MVKN' )
         IF ( IND > 0 ) MW_MVKN   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'R4N2' )
         IF ( IND > 0 ) MW_R4N2   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'ICHE' )
         IF ( IND > 0 ) MW_ICHE   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'ITHN' )
         IF ( IND > 0 ) MW_ITHN   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'ITCN' )
         IF ( IND > 0 ) MW_ITCN   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IDN' )
         IF ( IND > 0 ) MW_IDN    = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MCRHN' )
         IF ( IND > 0 ) MW_MCRHN  = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MCRHNB' )
         IF ( IND > 0 ) MW_MCRHNB = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'INPB' )
         IF ( IND > 0 ) MW_INPB   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'INPD' )
         IF ( IND > 0 ) MW_INPD   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IHN1' )
         IF ( IND > 0 ) MW_IHN1   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IHN2' )
         IF ( IND > 0 ) MW_IHN2   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IHN3' )
         IF ( IND > 0 ) MW_IHN3   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IHN4' )
         IF ( IND > 0 ) MW_IHN4   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MONITS' )
         IF ( IND > 0 ) MW_MONITS = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MONITU' )
         IF ( IND > 0 ) MW_MONITU = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HONIT' )
         IF ( IND > 0 ) MW_HONIT  = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IONITA' )
         IF ( IND > 0 ) MW_IONITA = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'MONITA' )
         IF ( IND > 0 ) MW_MONITA = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'BrNO3' )
         IF ( IND > 0 ) MW_BrNO3  = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HOBr' )
         IF ( IND > 0 ) THEN
            MW_HOBr   = State_Chm%SpcData(IND)%Info%MW_g

            ! Henry's law parameters
            H_K0_HOBr = State_Chm%SpcData(IND)%Info%Henry_K0 * con_atm_bar
            H_CR_HOBr = State_Chm%SpcData(IND)%Info%Henry_CR
            H_HOBr_T  = 298.15
         ENDIF

         IND = Ind_( 'O3' )
         IF ( IND > 0 ) THEN
            MW_O3     = State_Chm%SpcData(IND)%Info%MW_g

            ! Henry's law parameters
            H_K0_O3   = 1.1e-2_fp * con_atm_bar
            H_CR_O3   = 2300.0
            H_O3_T    = 298.15
         ENDIF

         IND = Ind_( 'HBr' )
         IF ( IND > 0 ) THEN
            MW_HBr    = State_Chm%SpcData(IND)%Info%MW_g

            ! Henry's law parameters
            H_K0_HBr  = State_Chm%SpcData(IND)%Info%Henry_K0
            H_CR_HBr  = State_Chm%SpcData(IND)%Info%Henry_CR
         ENDIF

         IND = Ind_( 'HCl' )
         IF ( IND > 0 ) THEN
            MW_HCl    = State_Chm%SpcData(IND)%Info%MW_g

            ! Henry's law parameters
            H_K0_HCl  = State_Chm%SpcData(IND)%Info%Henry_K0
            H_CR_HCl  = State_Chm%SpcData(IND)%Info%Henry_CR
         ENDIF

         IND = Ind_( 'ClNO3' )
         IF ( IND > 0 ) MW_ClNO3  = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HOCl' )
         IF ( IND > 0 ) MW_HOCl   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HI' )
         IF ( IND > 0 ) MW_HI     = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'HOI' )
         IF ( IND > 0 ) MW_HOI    = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'I2O2' )
         IF ( IND > 0 ) MW_I2O2   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'I2O3' )
         IF ( IND > 0 ) MW_I2O3   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'I2O4' )
         IF ( IND > 0 ) MW_I2O4   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IONO' )
         IF ( IND > 0 ) MW_IONO   = State_Chm%SpcData(IND)%Info%MW_g

         IND = Ind_( 'IONO2' )
         IF ( IND > 0 ) MW_IONO2  = State_Chm%SpcData(IND)%Info%MW_g

         ! Reset flag
         FIRST = .FALSE.

      ENDIF

      !--------------------------------------------------------------------
      ! Get species concentrations [molec/cm3]
      !--------------------------------------------------------------------
      IND = Ind_( 'NIT' )
      IF (IND .le. 0) THEN
         SPC_NIT    = 0.0e+0_fp
      ELSE
         SPC_NIT    = spcVec(IND)
      ENDIF

      !EEM: Add for ClNO2 yield calculation
      IND = Ind_( 'SALA' )
      IF (IND .le. 0) THEN
         SPC_SALA    = 0.0e+0_fp
      ELSE
         SPC_SALA    = spcVec(IND)
      ENDIF

      IND = Ind_('SO4')
      IF (IND .le. 0) THEN
         SPC_SO4    = 0.0e+0_fp
      ELSE
         SPC_SO4    = spcVec(IND)
      ENDIF

      IND = Ind_('HBr')
      IF (IND .le. 0) THEN
         SPC_HBr    = 0.0e+0_fp
      ELSE
         SPC_HBr    = spcVec(IND)
      ENDIF

      IND = Ind_('HOBr')
      IF (IND .le. 0) THEN
         SPC_HOBr   = 0.0e+0_fp
      ELSE
         SPC_HOBr   = spcVec(IND)
      ENDIF

      !--------------------------------------------------------------------
      ! Get fields required for epoxide uptake hydrolysis (EPOXUPTK)
      ! These values are saved in isorropiaII_mod.F
      !--------------------------------------------------------------------
      ! Proton activity [unitless] and H+ concentration [M]
      ! (assumed equivalent - for now):
      H_PLUS = State_Chm%HplusSav(I,J,L)

      ! Sulfate concentration [M]:
      MSO4   = State_Chm%SulRatSav(I,J,L)

      ! Nitrate concentration [M]:
      MNO3   = State_Chm%NaRatSav(I,J,L)

      ! Bisulfate (general acid) concentration [M]:
      MHSO4  = State_Chm%BisulSav(I,J,L)

      !--------------------------------------------------------------------
      ! Aerosol Physical Properties
      !--------------------------------------------------------------------

      ! Aerosol specific surface area, cm2(aerosol)/cm3(air)
      XAREA(1:State_Chm%nAeroType) = State_Chm%AeroArea(I,J,L,:)

      ! Aerosol effective radius, cm
      XRADI(1:State_Chm%nAeroType) = State_Chm%AeroRadi(I,J,L,:)

      ! Aerosol specific volume, cm3(aerosol)/cm3(air)
      XVOL(1:State_Chm%nAeroType)  = XAREA(1:State_Chm%nAeroType)              &
                                   * XRADI(1:State_Chm%nAeroType) / 3e+0_fp

      ! Aerosol water content, cm3(H2O)/cm3(air) [note: AeroH2O has units g/m3]
      XH2O(1:State_Chm%nAeroType)  = State_Chm%AeroH2O(I,J,L,:) * 1e-6_fp

      !--------------------------------------------------------------------
      ! Get fields from State_Met, State_Chm, and Input_Opt
      !--------------------------------------------------------------------

      TEMPK  = State_Met%T(I,J,L)              ! Temperature [K]
      XTEMP  = sqrt(State_Met%T(I,J,L))        ! Square root of temperature
      XDENA  = State_Met%AIRNUMDEN(I,J,L)      ! Dry air density [molec/cm3]
      SUNCOS = State_Met%SUNCOSmid(I,J)        ! COS(SZA),midpt of chem timestep
      VAir   = State_Met%AIRVOL(I,J,L)*1.0e6_fp! Volume of air (cm3)
      QICE   = State_Met%QI(I,J,L)             ! Ice   mix ratio [kg/kg dry air]
      QLIQ   = State_Met%QL(I,J,L)             ! Water mix ratio [kg/kg dry air]

      GAMMA_HO2 = Input_Opt%GAMMA_HO2

      !--------------------------------------------------------------------
      ! UCX-based mechanisms: Check surface type of PSCs (SDE 04/17/13)
      !--------------------------------------------------------------------
      IF ( Input_Opt%LUCX ) THEN
         CALL CHECK_NAT( I,  J,  L, NATSURFACE, PSCBOX, STRATBOX, &
                         Input_Opt, State_Met, State_Chm )
      ENDIF

      !--------------------------------------------------------------------
      !  Calculate parameters for cloud halogen chemistry
      !  under the new scheme (SDE 2016-12-21)
      !--------------------------------------------------------------------

      ! Get cloud physical parameters
      CALL Cld_Params( I, J, L, XDenA, VAir, TempK, QLiq, QIce, State_Met, &
                       rLiq,  ALiq,  VLiq, rIce,  AIce,  VIce, CLDFr )

      ! Retrieve cloud pH and alkalinity
      pHCloud    = State_Chm%pHCloud(I,J,L)
      SSAlk(1:2) = State_Chm%SSAlk(I,J,L,1:2)

      ! Estimate liquid phase pH (H+ concentration)
      hConc_Sul = 10.0**(-0.0e+0_fp)
      hConc_LCl = 10.0**(-1.0e+0_fp*pHCloud)
      hConc_ICl = 10.0**(-4.5e+0_fp)
      hConc_SSA = 10.0**(-5.0e+0_fp)
      hConc_SSC = 10.0**(-5.0e+0_fp)

      ! If not using BrSALA, manually set a depleted Br- concentration (mol/l)
      useSaltBr = ((.not.fixedSaltBr).and.(Ind_('BrSALA') > 0))
      IF (useSaltBr) THEN
         brConc_Base = 0.0e+0_fp
      ELSE
         brConc_Base = 1.0e+4_fp
      ENDIF

      ! Get the concentration of Br/Cl in clouds
      CALL Get_Halide_CldConc(spcVec(Ind_('HBr')),spcVec(Ind_('HCl')),&
                              VLiq, VIce, VAir, CLDFr, TempK, xArea(8),&
                              xRadi(8), brConc_Cld, clConc_Cld)

      ! Get the concentration of Br in sea-salt (in excess of any assumed
      ! baseline)
      CALL Get_Halide_SSAConc(spcVec(Ind_('BrSALA')),xArea(11),xRadi(11), &
                              brConc_SSA)
      CALL Get_Halide_SSAConc(spcVec(Ind_('BrSALC')),xArea(12),xRadi(12), &
                              brConc_SSC)

      !--------------------------------------------------------------------
      !  Get parameters for HOBr + S(IV)
      !--------------------------------------------------------------------

      ! Cloud bisulfite (HSO3-) concentration [mol/l] from sulfate_mod.F
      HSO3conc_Cld = State_Chm%HSO3_AQ(I,J,L)

      ! Cloud sulfite (SO3--) concentration [mol/l] from sulfate_mod.F
      SO3conc_Cld  = State_Chm%SO3_AQ(I,J,L)

      ! Avoid div-by-zero issues in GAMMA_HOBr_X
      !IF ( HSO3conc_Cld <= 0.0_fp) HSO3conc_Cld = 1e-20_fp
      !IF (  SO3conc_Cld  <= 0.0_fp)  SO3conc_Cld = 1e-20_fp

      ! Correction factor for HOBr removal by SO2 [unitless]
      fupdateHOBr  = State_Chm%fupdateHOBr(I,J,L)

      !--------------------------------------------------------------------
      ! Calculate and pass het rates to the KPP rate array
      !--------------------------------------------------------------------

      ! Zero the HET array
      HET = 0.0_dp

      ! Calculate genuine first-order uptake reactions first
      HET(ind_HO2,    1) = HetHO2(        MW_HO2,    2E-1_fp)
      HET(ind_NO2,    1) = HetNO2(        MW_NO2,    1E-4_fp)
      HET(ind_NO3,    1) = HetNO3(        MW_NO3,    1E-1_fp)
      HET(ind_GLYX,   1) = HetGLYX(       MW_GLYX,   1E-1_fp)
      HET(ind_MGLY,   1) = HetMGLY(       MW_MGLY,   1E-1_fp)
      HET(ind_IEPOXA, 1) = HetIEPOX(      MW_IEPOXA, 1E-1_fp)
      HET(ind_IEPOXB, 1) = HetIEPOX(      MW_IEPOXB, 1E-1_fp)
      HET(ind_IEPOXD, 1) = HetIEPOX(      MW_IEPOXD, 1E-1_fp)
      HET(ind_HMML,   1) = HetIMAE(       MW_HMML,   1E-1_fp)
      HET(ind_PYAC,   1) = HetMGLY(       MW_PYAC,   1E-1_fp)
      HET(ind_ICHE,   1) = HetIEPOX(      MW_ICHE,   1E-1_fp)
      HET(ind_LVOC,   1) = HetLVOC(       MW_LVOC,   1E+0_fp)
      HET(ind_IHN1,   1) = HetISOPND(     MW_IHN1,   5E-3_fp)
      HET(ind_IHN2,   1) = HetISOPNB(     MW_IHN2,   5E-2_fp)
      HET(ind_IHN3,   1) = HetISOPNB(     MW_IHN3,   5E-3_fp)
      HET(ind_IHN4,   1) = HetISOPND(     MW_IHN4,   5E-3_fp)
      HET(ind_INPB,   1) = HetISOPNB(     MW_INPB,   5E-3_fp)
      HET(ind_INPD,   1) = HetISOPND(     MW_INPD,   5E-3_fp)
      HET(ind_MCRHN,  1) = HetMACRN(      MW_MCRHN,  5E-3_fp)
      HET(ind_MCRHNB, 1) = HetMACRN(      MW_MCRHNB, 5E-3_fp)
      HET(ind_MVKN,   1) = HetMVKN(       MW_MVKN,   5E-3_fp)
      HET(ind_R4N2,   1) = HetR4N2(       MW_R4N2,   5E-3_fp)
      HET(ind_IDN,    1) = HetDHDN(       MW_IDN,    5E-3_fp)
      HET(ind_ITHN,   1) = HetDHDN(       MW_ITHN,   5E-3_fp)
      HET(ind_ITCN,   1) = HetDHDN(       MW_ITCN,   5E-3_fp)
      HET(ind_MONITS, 1) = HetMONITS(     MW_MONITS, 1E-2_fp)
      HET(ind_MONITU, 1) = HetMONITU(     MW_MONITU, 1E-2_fp)
      HET(ind_HONIT,  1) = HetHONIT(      MW_HONIT,  1E-2_fp)
      HET(ind_IONITA, 1) = HetIONITA(     MW_IONITA, 1E-1_fp)
      HET(ind_MONITA, 1) = HetMONITA(     MW_MONITA, 1E-1_fp)

      ! First-order loss in clouds
      HET(ind_NO3,    1) = HET(ind_NO3, 1) + &
                           CloudHet( 'NO3', CldFr, Aliq, Aice, rLiq, rIce, TempK, XDenA )
      HET(ind_N2O5,   4) = CloudHet( 'N2O5',CldFr, Aliq, Aice, rLiq, rIce, TempK, XDenA )

      ! Now calculate reaction rates where the educt can be consumed.
      ! kIIR1Ltd: Assume that the first reactant is limiting. Assume that the
      ! second reactant is "abundant" and calculate the overall rate based on
      ! the uptake rate of the first reactant only.
      HetTemp(1:2) = HETN2O5(1.08E2_fp, 1E-1_fp)
      HET(ind_N2O5,  1) = kIIR1Ltd( spcVec, Ind_('N2O5'), Ind_('H2O'), &
                                        HetTemp(1))
      State_Chm%GammaN2O5(I,J,L,1) = HetTemp(2)

      !--------------------------------------------------------------------
      ! Br/Cl heterogeneous chemistry
      !--------------------------------------------------------------------
      IF (Ind_('ClNO3') > 0) THEN

         !----------------------------------------------------------------
         ! ClNO3 and BrNO3 hydrolysis (SDE 2016-12-21)
         !----------------------------------------------------------------
         kITemp = HETBrNO3_JS( XDenA, rLiq, rIce, ALiq, AIce, TempK )
         HET(ind_BrNO3, 1) = kIIR1Ltd( spcVec, Ind_('BrNO3'), Ind_('H2O'), &
                                       kITemp, HetMinLife)
         kITemp = HETClNO3_JS( XDenA, rLiq, rIce, ALiq, AIce, TempK )
         HET(ind_ClNO3, 1) = kIIR1Ltd( spcVec, Ind_('ClNO3'), Ind_('H2O'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! HOBr + HBr (TMS index: hhc06)
         !----------------------------------------------------------------
         kITemp = HETHOBr_HBr_JS( XDenA, rLiq, rIce, ALiq, AIce, VAir, TempK, &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc_Cld, &
                           brConc_Cld, HSO3conc_Cld, SO3conc_Cld )
         HET(ind_HOBr,  1) = kIIR1Ltd( spcVec, Ind_('HOBr'),  Ind_('HBr'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! HOBr + HCl (TMS index: hhc03)
         !----------------------------------------------------------------
         kITemp = HETHOBr_HCl_JS( XDenA, rLiq, rIce, ALiq, AIce, VAir, TempK, &
                                  hConc_Sul, hConc_LCl, hConc_ICl, clConc_Cld, &
                                  brConc_Cld, HSO3conc_Cld, SO3conc_Cld )
         HET(ind_HOBr,  2) = kIIR1Ltd( spcVec, Ind_('HOBr'),  Ind_('HCl'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! HOBr + BrSalA/C (TMS index: hhc07/08)
         !----------------------------------------------------------------
         ! NOTE: This has not been fully tested, as the initial simulations had
         ! near-zero BrSALA and BrSALC
         kITemp = HETHOBr_SS_JS( XDenA, xRadi(11), xArea(11), SSAlk(1), TempK, &
                                 hConc_SSA, 0.5e+0_fp, brConc_SSA, 2 )
         HET(ind_HOBr,  4) = kIIR1Ltd( spcVec, Ind_('HOBr'),  Ind_('BrSALA'), &
                                       kITemp, HetMinLife)

         kITemp = HETHOBr_SS_JS( XDenA, xRadi(12), xArea(12), SSAlk(2), TempK, &
                                 hConc_SSC, 0.5e+0_fp, brConc_SSC, 2 )
         HET(ind_HOBr,  5) = kIIR1Ltd( spcVec, Ind_('HOBr'),  Ind_('BrSALC'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! HOBr + ClSALA/C (TMS index: hhc04/05)
         !----------------------------------------------------------------
         ! NOTE: Cl- in salt is assumed to always be in excess, so we assume a
         ! molarity of 0.5 mol/L. This reaction is also pseudo-first order, so
         ! conversion to a second-order rate constant is not necessary.
         kITemp = HETHOBr_SS_JS( XDenA, xRadi(11), xArea(11), SSAlk(1), &
                                 TempK, hConc_SSA, 0.5e+0_fp, brConc_SSA, 1 )
         kITemp = kITemp + &
                  HETHOBr_SS_JS( XDenA, xRadi(12), xArea(12), SSAlk(2), &
                                 TempK, hConc_SSC, 0.5e+0_fp, brConc_SSC, 1 )
         HET(ind_HOBr,  3) = kITemp

         !----------------------------------------------------------------
         ! HOBr + HSO3-(aq) (QJC index: EhcHSHOBCld)
         !----------------------------------------------------------------
         ! This reaction is first order, so no kII calculation is required
         kITemp = HETHOBr_HSO3( XDenA, rLiq, rIce, ALiq, AIce, VAir, TempK, &
                             hConc_Sul, hConc_LCl, hConc_ICl, clConc_Cld, &
                             brConc_Cld, HSO3conc_Cld, SO3conc_Cld )

         ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
         HET(ind_HOBr,  6) = kITemp * fupdateHOBr

         !----------------------------------------------------------------
         ! HOBr + SO3--(aq) (QJC index: EhcSOHOBCld)
         !----------------------------------------------------------------
         ! This reaction is first order, so no kII calculation is required
         kITemp = HETHOBr_SO3( XDenA, rLiq, rIce, ALiq, AIce, VAir, TempK, &
                             hConc_Sul, hConc_LCl, hConc_ICl, clConc_Cld, &
                             brConc_Cld, HSO3conc_Cld, SO3conc_Cld )

         ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
         HET(ind_HOBr,  7) = kITemp * fupdateHOBr

         !----------------------------------------------------------------
         ! ClNO3 + BrSALA/C (TMS index: hhc10/11)
         !----------------------------------------------------------------
         ! NOTE: This has not been fully tested, as the initial simulations had
         ! near-zero BrSALA and BrSALC
         kITemp = HETClNO3_SS_JS( XDenA, xRadi(11), xArea(11), SSAlk(1), &
                                  TempK, brConc_SSA)
         HET(ind_ClNO3, 4) = kIIR1Ltd( spcVec, Ind_('ClNO3'), Ind_('BrSALA'), &
                                       kITemp, HetMinLife)
         kITemp = HETClNO3_SS_JS( XDenA, xRadi(12), xArea(12), SSAlk(2), &
                                  TempK, brConc_SSC)
         HET(ind_ClNO3, 5) = kIIR1Ltd( spcVec, Ind_('ClNO3'), Ind_('BrSALC'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! ClNO3 + HCl
         !----------------------------------------------------------------
	 ! NOTE: the restriction of these reactions to the troposphere has been
         ! restored - TMS (2017/04/06 )
         HET(ind_ClNO3, 2) = kIIR1Ltd( spcVec, Ind_('ClNO3'), Ind_('HCl'), &
                             HETClNO3_HCl( 0.97E2_fp, 0E+0_fp), HetMinLife)

         !----------------------------------------------------------------
         ! ClNO3 + HBr (TMS index: hhc09)
         !----------------------------------------------------------------
         kITemp = HETClNO3_HBr_JS( xDenA, rLiq, rIce, ALiq, AIce, VAir, &
                                   TempK, brConc_Cld, Input_Opt )
         HET(ind_ClNO3, 3) = kIIR1Ltd( spcVec, Ind_('ClNO3'), Ind_('HBr'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! HOCl + HCl and HOCl + HBr to take place in the troposphere
         !----------------------------------------------------------------
	 ! NOTE: the restriction of these reactions to the troposphere has been
         ! restored - TMS (2017/04/06 )
         HET(ind_HOCl,  1) = kIIR1Ltd( spcVec, Ind_('HOCl'),  Ind_('HCl'), &
                             HETHOCl_HCl(  0.52E2_fp, 0E+0_fp, Input_Opt), &
                             HetMinLife)
         HET(ind_HOCl,  2) = kIIR1Ltd( spcVec, Ind_('HOCl'),  Ind_('HBr'), &
                             HETHOCl_HBr(  0.52E2_fp, 0E+0_fp, Input_Opt), &
                             HetMinLife)

         !----------------------------------------------------------------
         ! O3 + Br- calculation (TMS index: hhc12)
         !----------------------------------------------------------------
         kITemp = HETO3_HBr_JS( XDenA, rLiq, rIce, ALiq, AIce, VAir, &
                                TempK, brConc_Cld, spcVec(Ind_('O3')))
         HET(ind_O3,    1) = kIIR1Ltd( spcVec, Ind_('O3'), Ind_('HBr'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! O3 + BrSALA/C calculations (TMS index: hhc13/14)
         !----------------------------------------------------------------
         kITemp = HETO3_SS_JS( XDenA, xRadi(11), xArea(11), SSAlk(1), &
                               TempK, brConc_SSA, spcVec(Ind_('O3')))
         HET(ind_O3,    2) = kIIR1Ltd( spcVec, Ind_('O3'), Ind_('BrSALA'), &
                                       kITemp, HetMinLife)
         kITemp = HETO3_SS_JS( XDenA, xRadi(12), xArea(12), SSAlk(2), &
                               TempK, brConc_SSC, spcVec(Ind_('O3')))
         HET(ind_O3,    3) = kIIR1Ltd( spcVec, Ind_('O3'), Ind_('BrSALC'), &
                                       kITemp, HetMinLife)

         !----------------------------------------------------------------
         ! Cl uptake calculations (TMS index: hhc15/16)
         !----------------------------------------------------------------
         ! Cl is always assumed to be in excess in sea salt, so any HCl "taken
         ! up" is just removed. This may change in the future. This reaction is
         ! also first order, so no kII calculation is required
         kITemp = HETHXUptake_JS( XDenA, xRadi(11), xArea(11), TempK, 1)
         HET(ind_HCl,   1) = kITemp
         kITemp = HETHXUptake_JS( XDenA, xRadi(12), xArea(12), TempK, 1)
         HET(ind_HCl,   2) = kITemp

         !----------------------------------------------------------------
         ! Br uptake calculation - forms BrSALA/C (TMS index: hhc17/18)
         !----------------------------------------------------------------
         ! First-order reactions, no calculation of kII required
         kITemp = HETHXUptake_JS( XDenA, xRadi(11), xArea(11), TempK, 2)
         HET(ind_HBr,   1) = kITemp
         kITemp = HETHXUptake_JS( XDenA, xRadi(12), xArea(12), TempK, 2)
         HET(ind_HBr,   2) = kITemp

         !----------------------------------------------------------------
         ! BrNO3 + HCl into the troposphere
         !----------------------------------------------------------------
	 ! NOTE: the restriction of these reactions to the troposphere has been
         ! restored - TMS (2017/04/06 )
         HET(ind_BrNO3, 2) = kIIR1Ltd( spcVec, Ind_('BrNO3'), Ind_('HCl'), &
                             HETBrNO3_HCl(  1.42E2_fp, 0E+0_fp), HetMinLife)

         !----------------------------------------------------------------
         ! N2O5 + HCl on sulfate
         !----------------------------------------------------------------
	 ! NOTE: this extension of calculation in troposphere has been removed
         !  (TMS 17/04/10)
         ! 1st HetTemp() parameter is kN2O5 loss rate coefficient, 2nd is gamma
         HetTemp(1:2) = HETN2O5_HCl( 1.08E2_fp, 0.0e+0_fp, Input_Opt )
         HET(ind_N2O5,  2) = kIIR1Ltd( spcVec, Ind_('N2O5'), Ind_('HCl'), &
                                       HetTemp(1), HetMinLife)
         State_Chm%GammaN2O5(I,J,L,2) = HetTemp(2)

         !----------------------------------------------------------------
         ! Reaction of N2O5 with sea-salt Cl-
         ! (assumed to be in excess, so no kII calculation)
         !----------------------------------------------------------------
         ! 1st HetTemp() parameter is kN2O5 loss rate coefficient, 2nd is gamma
         ! 3rd is the ClNO2 yield (currently set to zero, add in future update)
         HetTemp(1:3) = HETN2O5_SS(1.08E2_fp, 1E-1_fp)
         HET(ind_N2O5,  3) = HetTemp(1)
         State_Chm%GammaN2O5(I,J,L,3) = HetTemp(2)
         State_Chm%GammaN2O5(I,J,L,4) = HetTemp(3)
      ENDIF

      ! Iodine chemistry
      IF (Ind_('I2').gt.0) THEN

         ! Uptake reactions (forming AERI, ISALA and ISALC)
         HET(ind_HI,   1) = HETIUptake( MW_HI,   0.10e+0_fp,  8, Input_Opt )
         HET(ind_HI,   2) = HETIUptake( MW_HI,   0.10e+0_fp, 11, Input_Opt )
         HET(ind_HI,   3) = HETIUptake( MW_HI,   0.10e+0_fp, 12, Input_Opt )
         HET(ind_I2O2, 1) = HETIUptake( MW_I2O2, 0.02e+0_fp,  8, Input_Opt )
         HET(ind_I2O2, 2) = HETIUptake( MW_I2O2, 0.02e+0_fp, 11, Input_Opt )
         HET(ind_I2O2, 3) = HETIUptake( MW_I2O2, 0.02e+0_fp, 12, Input_Opt )
         HET(ind_I2O3, 1) = HETIUptake( MW_I2O3, 0.02e+0_fp,  8, Input_Opt )
         HET(ind_I2O3, 2) = HETIUptake( MW_I2O3, 0.02e+0_fp, 11, Input_Opt )
         HET(ind_I2O3, 3) = HETIUptake( MW_I2O3, 0.02e+0_fp, 12, Input_Opt )
         HET(ind_I2O4, 1) = HETIUptake( MW_I2O4, 0.02e+0_fp,  8, Input_Opt )
         HET(ind_I2O4, 2) = HETIUptake( MW_I2O4, 0.02e+0_fp, 11, Input_Opt )
         HET(ind_I2O4, 3) = HETIUptake( MW_I2O4, 0.02e+0_fp, 12, Input_Opt )

         ! These uptake reactions require non-acidic aerosol
         ! Fine sea salt first
         IF (SSAlk(1).gt.0.05) THEN
            HET(ind_HOI,  1) = HETIUptake( MW_HOI,   0.01e+0_fp, 11, Input_Opt )
            HET(ind_IONO, 1) = HETIUptake( MW_IONO,  0.02e+0_fp, 11, Input_Opt )
            HET(ind_IONO2,1) = HETIUptake( MW_IONO2, 0.01e+0_fp, 11, Input_Opt )
         ENDIF

         ! Now coarse sea salt
         IF (SSAlk(2).gt.0.05) THEN
            HET(ind_HOI,  2) = HETIUptake( MW_HOI,   0.01e+0_fp, 12, Input_Opt )
            HET(ind_IONO, 2) = HETIUptake( MW_IONO,  0.02e+0_fp, 12, Input_Opt )
            HET(ind_IONO2,2) = HETIUptake( MW_IONO2, 0.01e+0_fp, 12, Input_Opt )
         ENDIF

         ! Breakdown of iodine compounds on sea-salt
         HET(ind_HOI,  3) = HETIXCycleSSA( MW_HOI,   0.01E+0_fp, SSAlk )
         HET(ind_IONO, 3) = HETIXCycleSSA( MW_IONO,  0.02E+0_fp, SSAlk )
         HET(ind_IONO2,3) = HETIXCycleSSA( MW_IONO2, 0.01E+0_fp, SSAlk )

      ENDIF

      ! Nullify pointers
      NULLIFY( spcVec )

      RETURN

    END SUBROUTINE SET_HET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CloudHet
!
! !DESCRIPTION: Function CloudHet calculates the loss frequency (1/s) of gas
!  species due to heterogeneous chemistry on clouds in a partially cloudy grid
!  cell. The function uses the "entrainment limited uptake" equations of
!  Holmes et al. (2018). Both liquid and ice water clouds are treated.
!\\
!\\
! !INTERFACE:
!
    function CloudHet( SpeciesName, fc, Aliq, Aice, rLiq, rIce, T, airNumDen ) result( kHet )
!
! !INPUT PARAMETERS:
!
      character(len=*),intent(in) :: SpeciesName
      real(fp),intent(in)         :: fc, &   ! Cloud Fraction [0-1]
                                     Aliq, & ! Surface area density of cloud liquid & ice, cm2/cm3
                                     Aice, & !  (grid average, not in-cloud)
                                     rLiq, & ! Effective radius for liquid and ice clouds, cm
                                     rIce, &
                                     T,    & ! Temperature, K
                                     airNumDen ! Air number density, molec/cm3
!
! !RETURN VALUE:
!
      real(fp)                    :: kHet ! Grid-average loss frequency, 1/s
!
! !REMARKS:
!
! !REVISION HISTORY:
!  23 Aug 2018 - C. D. Holmes - Initial version
!  17 Oct 2018 - C. D. Holmes - Re-implemented for v12.0.2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Residence time of air in clouds, s
      real(fp),parameter  :: tauc = 3600

      ! Molar mass of species, kg/mol
      real(fp), parameter :: mN2O5 = 0.108e+0_fp
      real(fp), parameter :: mNO2  = 0.046e+0_fp
      real(fp), parameter :: mNO3  = 0.062e+0_fp
      real(fp), parameter :: mHO2  = 0.033e+0_fp
!
! !LOCAL VARIABLES:
!
      real(fp) :: kI, gam, rd, gammaLiq, gammaIce, &
           area, alpha, beta, molmass
      real(fp) :: kk, ff, xx
      integer  :: K
!
!------------------------------------------------------------------------------
!
      ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface area,
      ! then return zero uptake
      if ( (fc < 0.0001) .or. ((ALiq + AIce) <= 0) ) then
         kHet = 0
         return
      endif

      !------------------------------------------------------------------------
      ! Select Gamma and molar mass for this species
      !------------------------------------------------------------------------

      select case (trim(speciesName))

      case ('HO2')

         gammaLiq = 0.1e+0_fp
         gammaIce = 0.025e+0_fp
         molmass  = mHO2

      case ('NO2')

         gammaLiq = 1e-8_fp
         gammaIce = 0.0e+0_fp
         molmass  = mNO2

      case ('NO3')

         gammaLiq = 0.002e+0_fp
         gammaIce = 0.001e+0_fp
         molmass  = mNO3

      case ('N2O5')

         ! Reactive uptake coefficient for N2O5 on liquid water cloud
         ! Value is 0.03 at 298 K (JPL, Burkholder et al., 2015)
         ! For temperature dependence, JPL recommends the same as
         ! sulfuric acid aerosol at zero percent H2SO4, which is 0.019 at 298 K.
         ! Then apply constant scale factor (0.03/0.019)
         gammaLiq = ( 0.03e+0_fp / 0.019e+0_fp ) * &
              exp( -25.5265e+0_fp + 9283.76e+0_fp / T - 851801e+0_fp / T**2 )

         ! Reactive uptake coefficient for N2O5 on water ice
         gammaIce = 0.02e+0_fp

         molmass  = mN2O5

      case default

         print*, speciesName // ' not found in CloudHet function'
         call GEOS_CHEM_STOP

      end select

      !------------------------------------------------------------------------
      ! Loss frequency inside cloud
      !
      ! Assume both water and ice phases are inside the same cloud, so mass
      ! transport to both phases works in parallel (additive)
      !------------------------------------------------------------------------

      ! initialize loss, 1/s
      kI = 0

      ! Loop over water (K=1) and ice (K=2)
      do K=1, 2

         ! Liquid water cloud
         if (K==1) then

            ! gamma, unitless
            gam = gammaLiq

            ! Liquid particle radius, cm
            rd = rLiq

            ! Liquid surface area density, cm2/cm3
            area = Aliq

         ! Ice water cloud
         elseif (K==2) then

            ! gamma, unitless
            gam = gammaIce

            ! Ice particle effective radius, cm
            rd = rIce

            ! Ice surface area density, cm2/cm3
            area = Aice

         else

            print*, 'CloudHet: index value exceeded'
            call GEOS_CHEM_STOP

         endif
!         print'(I4,6E10.3)',K,area,rd,airnumden,gam,T,molmass

         ! Skip calculation if there is no surface area
         if ( area <= 0.0e+0_fp ) cycle

         ! In-cloud loss frequency, combining ice and liquid in parallel, 1/s
         ! Pass radius in cm and mass in g.
         kI = kI + arsl1k( area, rd, airnumden, gam, sqrt(T), sqrt(molmass*1000) )

      end do

!      !------------------------------------------------------------------------
!      ! Grid-average loss frequency; Add in-cloud and entrainment rates in series
!      !
!      ! APPROXIMATE expression for entrainment-limited uptake
!      !   Approximation error in loss frequency is typically <2% and always <50%.
!      !------------------------------------------------------------------------
!
!      ! Requires declaring...
!      ! real(fp) :: kIinv, kEinv
!
!      ! Entrainment rate, inverse, s
!      kEinv = safe_div( tauc * ( 1e+0_fp - fc ), fc, 1e+30_fp )
!
!      ! In-cloud loss rate, inverse, s
!      kIinv = safe_div( 1e+0_fp, fc*kI, 1e+30_fp )
!
!      ! Overall heterogeneous loss rate, grid average, 1/s
!      kHet = safe_div( 1e+0_fp, ( kEinv + kIinv ), 0e+0_fp )
!

      !------------------------------------------------------------------------
      ! Grid-average loss frequency
      !
      ! EXACT expression for entrainment-limited uptake
      !------------------------------------------------------------------------

      ! Ratio (in cloud) of heterogeneous loss to detrainment, s/s
      kk = kI * tauc

      ! Ratio of volume inside to outside cloud
      ! ff has a range [0,+inf], so cap it at 1e30
      ff = safe_div( fc, (1e0_fp - fc), 1e30_fp )
      ff = max( ff, 1e30_fp )

      ! Ratio of mass inside to outside cloud
      ! xx has range [0,+inf], but ff is capped at 1e30, so this shouldn't overflow
      xx =     ( ff - kk - 1e0_fp ) / 2e0_fp + &
           sqrt( 1e0_fp + ff**2 + kk**2 + 2*ff**2 + 2*kk**2 - 2*ff*kk ) / 2e0_fp

      ! Overall heterogeneous loss rate, grid average, 1/s
      ! kHet = kI * xx / ( 1e+0_fp + xx )
      !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
      !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
      kHet = kI / ( 1e0_fp + safe_div( 1e0_fp, xx, 1e30_fp ) )

    end function CloudHet
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kiir1ltd
!
! !DESCRIPTION: Determine removal rates for both species in an uptake reaction.
!\\
!\\
! !INTERFACE:
!
    FUNCTION kIIR1Ltd( spcVec, indGas, indEduct, kISource, minLife ) &
       RESULT( kII )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN)           :: spcVec(:)
      INTEGER,  INTENT(IN)           :: indGas
      INTEGER,  INTENT(IN)           :: indEduct
      REAL(fp), INTENT(IN)           :: kISource
      REAL(fp), INTENT(IN), OPTIONAL :: minLife
!
! !RETURN VALUE:
!
      REAL(fp)                       :: kII
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: kIGas, kIEduct, concGas, concEduct
      REAL(fp) :: lifeA, lifeB, kIMult

      concGas = spcVec(indGas)
      concEduct = spcVec(indEduct)

      ! Copy kI as calculated assuming no limitation
      kIGas = kISource
      kIEduct = 0.0e+0_fp
      kII = 0.0e+0_fp

      IF (concEduct.lt.100.0e+0_fp) THEN
         kIGas = 0.0e+0_fp
         kIEduct = 0.0e+0_fp
         kII = 0.0e+0_fp
      ELSE
         ! Safe division here is probably overkill - may remove this
         IF (Is_Safe_Div(concGas*kIGas,concEduct)) THEN
            kIEduct = kIGas*concGas/concEduct
            kII = kIGas/concEduct
         ELSE
            kIGas = 0.0e+0_fp
            kIEduct = 0.0e+0_fp
            kII = 0.0e+0_fp
         ENDIF
      ENDIF

      ! Enforce a minimum lifetime?
      IF (PRESENT(minLife)) THEN
         IF ((kIGas.gt.0.0e+0_fp).and.(minLife.gt.0.0e+0_fp)) THEN
            ! Calculate lifetime of each reactant against removal
            lifeA = Safe_Div(1.0e+0_fp,kIGas,0.0e+0_fp)
            lifeB = Safe_Div(1.0e+0_fp,kIEduct,0.0e+0_fp)
            ! Check if either lifetime is "too short"
            IF ((lifeA.lt.lifeB).and.(lifeA.lt.minLife)) THEN
               IF (Is_Safe_Div(concGas*kIGas,concEduct)) THEN
                  kIGas = 1.0e+0_fp/minLife
                  kII = kIGas/concEduct
               ELSE
                  kIGas = 0.0e+0_fp
                  kII = 0.0e+0_fp
               ENDIF
            ELSEIF (lifeB.lt.minLife) THEN
               IF (Is_Safe_Div(concEduct*kIEduct,concGas)) THEN
                  kIEduct = 1.0e+0_fp/minLife
                  kII = kIEduct/concGas
               ELSE
                  kIEduct = 0.0e+0_fp
                  kII = 0.0e+0_fp
               ENDIF
            ENDIF
         ENDIF
      ENDIF

    END FUNCTION kIIR1Ltd
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kiir1r2l
!
! !DESCRIPTION: Determine removal rates for both species in an uptake reaction
! without assuming which reactant is limiting.
!\\
!\\
! !INTERFACE:
!
    FUNCTION kIIR1R2L( spcVec, indGasA, indGasB, kIASource, kIBSource ) &
       RESULT( kII )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN)    :: spcVec(:)
      INTEGER,  INTENT(IN)    :: indGasA
      INTEGER,  INTENT(IN)    :: indGasB
      REAL(fp), INTENT(IN)    :: kIASource
      REAL(fp), INTENT(IN)    :: kIBSource
!
! !RETURN VALUE:
!
      REAL(fp)                :: kII

!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: concGasA, concGasB, kIA, kIB
      REAL(fp) :: R_GasA, R_GasB
      LOGICAL  :: nonZeroRate

      ! Get the base concentrations
      concGasA = spcVec(indGasA)
      concGasB = spcVec(indGasB)

      ! Copy the first estimates of kI for each species
      kIA = kIASource
      kIB = kIBSource

      ! Assume for now that the reaction will not proceed
      nonZeroRate = .False.
      kII = 0.0e+0_fp

      ! Prevent reaction if either concentration is too low
      IF ((concGasA.gt.100.0e+0_fp).and.(concGasB.gt.100.0e+0_fp).and.&
          (kIA.gt.0.0e+0_fp).and.(kIB.gt.0.0e+0_fp)) THEN
         ! Calculate the overall rate based on each reactant
         R_GasA = kIA*concGasA
         R_GasB = kIB*concGasB
         IF (R_GasA > R_GasB) THEN
            ! Limited by uptake of B
            nonZeroRate = Is_Safe_Div( R_GasB, concGasA )
            IF (nonZeroRate) THEN
               kII = kIB / concGasA
            ENDIF
         ELSE
            ! Limited by uptake of A
            nonZeroRate = Is_Safe_Div( R_GasA, concGasB )
            IF (nonZeroRate) THEN
               kII = kIA / concGasB
            ENDIF
         ENDIF
      ENDIF

      ! If no tests were passed, zero out both rates
      IF (.not.nonZeroRate) THEN
         kII = 0.0e+0_fp
      ENDIF

    END FUNCTION kIIR1R2L
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIXCycleSSA
!
! !DESCRIPTION: Set the iodine reaction rate on sea salt, assuming a fixed ratio
! of ICl and IBr (85:15) is produced.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIXCycleSSA( A, B, SSAlk ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
      ! Sea salt alkalinity
      REAL(fp), INTENT(IN) :: SSAlk(2)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      INTEGER  :: N, NAer

      ! Initialize
      kISum        = 0.0_fp

      DO N=1,2
         NAer = N + 10
         ! Only allow reaction on acidic aerosol
         IF (SSAlk(N).le.0.05e+0_fp) THEN
            ! Reaction rate for surface of aerosol
            AdjustedRate = ARSL1K(XAREA(NAer),XRADI(NAer),XDENA,B,XTEMP,(A**0.5_fp))
            kISum = kISum + AdjustedRate
         ENDIF
      ENDDO

    END FUNCTION HETIXCycleSSA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIUptake
!
! !DESCRIPTION: Set the uptake rate for iodine species.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIUptake( A, B, N, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      INTEGER,        INTENT(IN) :: N          ! Which aerosol?
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version
!  03 Jan 2018 - M. Sulprizio  - Replace UCX CPP switch with Input_Opt%LUCX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum        = 0.0_fp

      ! Reaction rate for surface of aerosol
      kISum = ARSL1K(XAREA(N),XRADI(N),XDENA,B,XTEMP,(A**0.5_fp))

      IF ( Input_Opt%LUCX ) THEN
         ! For UCX-based mechanisms also allow reaction on stratospheric
         ! sulfate (N=13) if tropospheric sulfate is requested (N=8)
         IF (N.eq.8) THEN
            kISum = kISum + ARSL1K(XAREA(13),XRADI(13),XDENA,B,XTEMP, &
                    (A**0.5_fp))
         ENDIF
      ENDIF

    END FUNCTION HETIUptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetNO3
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO3.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO3( A, B ) RESULT( HET_NO3 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_NO3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca  - Added ProTeX header
!  01 Apr 2016 - R. Yantosca  - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca  - Replace KII_KI with DO_EDUCT local variable
!  23 Aug 2018 - C. D. Holmes - Updated Gamma values
!  06 Nov 2019 - R. Yantosca  - Force flexible precision with _fp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_NO3      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = B

         ! Set uptake coefficients
         select case (N)
            case (1:7)
               ! dust
               xstkcf = 0.01_fp
            case (8)
               ! sulfate
               if ( relhum < 40.0_fp ) then
                  xstkcf = 0.001_fp
               else
                  xstkcf = 0.002_fp
               endif
            case (9)
               ! BC
               if ( relhum < 50.0_fp ) then
                  xstkcf = 2e-4_fp
               else
                  xstkcf = 1e-3_fp
               endif
            case (10)
               ! OC
               xstkcf = 0.005_fp
            case (11:12)
               ! sea salt
               if ( relhum < 40.0_fp ) then
                  xstkcf = 0.05_fp
               elseif ( relhum > 70.0_fp ) then
                  xstkcf = 0.002_fp
               else
                  xstkcf = 0.05_fp + (0.002_fp - 0.05_fp)                    &
                         * (relhum - 40.0_fp) / (70.0_fp - 40.0_fp)
               endif
         end select

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_NO3 = HET_NO3 + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetNO2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO2( A, B ) RESULT( HET_NO2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_NO2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  23 Aug 2018 - C. D. Holmes - Updated Gamma values
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_NO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = B

         ! Set uptake coefficients
         select case (N)
            case (1:7)
               ! dust
               xstkcf = 1e-8_fp
            case (8)
               ! sulfate
               xstkcf = 5e-6_fp
            case (9)
               ! BC
               xstkcf = 1e-4_fp
            case (10)
               ! OC
               xstkcf = 1e-6_fp
            case (11:12)
               ! sea salt
               if ( relhum < 40 ) then
                  xstkcf = 1e-8_fp
               elseif ( relhum > 70 ) then
                  xstkcf = 1e-4_fp
               else
                  xstkcf = 1e-8_fp + (1e-4_fp-1e-8_fp) * (relhum-40)/(70-40)
               endif
         end select

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_fp))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_NO2 = HET_NO2 + ADJUSTEDRATE

      END DO

    END FUNCTION HETNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHO2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for HO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHO2( A, B ) RESULT( HET_HO2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HO2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         IF (N.gt.12) THEN
            XSTKCF = TINY(1e+0_fp)
         ELSE
            XSTKCF = GAMMA_HO2
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HO2 = HET_HO2 + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETHO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHBr
!
! !DESCRIPTION: Set the heterogeneous rate for HBr.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHBr( A, B ) RESULT( HET_HBr )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HBr
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Define local variable for educt adjustment
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HBr      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Only apply PSC rate adjustment if at high altitude
      DO_EDUCT     = STRATBOX

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! jpp, 3/22/11: set the sticking coefficient to
         !  ~0 for aerosol types we don't want reactions on
         !  for the HBr and HOBr surface reaction

         ! Select proper aerosol type
         IF ( (N == 8) .OR. (N == 11) .OR. (N == 12)) THEN
            ! sulfate, 2 modes of sea-salt
            XSTKCF = B
         ELSEIF ( N == 13 ) THEN
            XSTKCF = KHETI_SLA(11)
         ELSEIF ( N == 14 ) THEN
            XSTKCF = 0.1e+0_fp
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HBr = HET_HBr + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETHBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5( A, B ) RESULT( HET_N2O5 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      ! HET_N2O5(1) = rate coefficient; HET_N2O5(2) = SA-weighted gamma
      REAL(fp)             :: HET_N2O5(2)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  02 Jan 2019 - E. McDuffie - Update for gamma for SO4-NIT-NH4 & OC aerosol
!                              (following McDuffie et al., JGR 2018)
!  23 Jul 2019 - C.D. Holmes - Consolidated McDuffie gamma into one function
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      REAL(fp) :: output(4)
      REAL(fp) :: ClNO2_yield, kClNO2, Rp, SA, SA_sum

!------------------------------------------------------------------------------

      ! Initialize
      HET_N2O5      = 0.0_fp    !Output, holds ADJUSTEDRATE and XSTCKF results
      ADJUSTEDRATE  = 0.0_fp    !Total N2O5 loss rate coefficient for all aerosol types
      XSTKCF        = 0.0_fp    !Total N2O5 gamma for all aerosol types
      ClNO2_yield   = 0.0_fp    !ClNO2 production yield (for future update)
      kClNO2        = 0.0_fp
      output        = 0.0_fp    !temporary variable
      Rp            = 0.0_fp    !Total SNA+ORG radius
      SA            = 0.0_fp    !Total SNA+ORG surface area
      SA_sum        = 0.0_fp    ! Used for weighted mean

      ! Always apply PSC rate adjustment
      DO_EDUCT     = .TRUE.

      !NOTE: This calculations follows these general steps:
      !Loop over all aerosol types
      ! A) Calculate gamma for each aerosol type
      ! B) Use gamma to calculate N2O5 loss rate coefficient for each aerosol
      !    type
      ! C) Add results from all aerosol types to get running sum of ADJUSTEDRATE
      !    and surface-area-weighted gamma
      !
      !After looping through all aerosol types:
      ! D) Divide gamma by SA to get single SA-weighted gamma

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE
         ! A) Calculate gamma for each aerosol type

         ! Get GAMMA for N2O5 hydrolysis, which is
         ! a function of aerosol type, temp, and RH
         IF (N.eq.14) THEN
            IF (NATSURFACE) THEN
               XSTKCF = 4.0e-4_fp ! NAT
            ELSE
               XSTKCF = 0.02e+0_fp ! Ice
            ENDIF
         ELSEIF (N.eq.13) THEN
            ! Stratospheric aerosol
            XSTKCF = KHETI_SLA(1)
         ELSEIF ((N==8) .or. (N==10) .or. (N==11) .or. (N==12)) THEN
            ! Set all of these to zero because they are handled elsewhere
            ! For inorganic and organic aerosol (N=8,10), Calculation occurs
            ! below
            ! For Sea salt (N=11,12) - follows the N2O5 + Cl- channel
            XSTKCF = 0.0e+0_fp
         ELSE
            ! For UCX-based mechanisms ABSHUMK is set to Spc(I,J,L,id_H2O)
            XSTKCF = N2O5( N, TEMPK, RELHUM )
         ENDIF

         ! B) Use gamma to calculate N2O5 loss rate coefficient for each
         !    aerosol type
         IF (N.eq.8) THEN

            ! Properties of inorganic (sulfate-nitrate-ammonium-sea salt) coated
            ! with organics
            output = N2O5_InorgOrg( XVOL(8), XVOL(10), XH2O(8), XH2O(10), &
                 XRADI(8),    SPC_NIT,    0e+0_fp,   TEMPK,    RELHUM )

            ! Gamma
            XSTKCF = output(1)
            ! ClNO2 yield, fraction [0,1]
            ClNO2_yield = output(2)
            ! Particle radius with coating, cm
            Rp     = output(3)
            ! Surface area of coated particles, cm2/cm3
            SA     = output(4)

            !For SNA aerosol...
            ! Total loss rate of N2O5 (kN2O5) on SNA+ORG aerosol
            ADJUSTEDRATE=ARSL1K( SA, Rp, XDENA, XSTKCF, XTEMP, (A**0.5_FP) )

            !Calculate ClNO2 yield on SNA(+ORG) aerosol using kN2O5 from above
            ! phi = kClNO2/kN2O5 (from ClNO2 and HNO3 production pathways)
            kClNO2 = ClNO2_yield * ADJUSTEDRATE

            ! reduce the rate of this HNO3 pathway in accordance with the yield
            ADJUSTEDRATE = ADJUSTEDRATE - kClNO2
         ELSEIF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            !For all other aerosol types...
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         ! C) Add results from all aerosol types to get running sum of
         !    ADJUSTEDRATE and surface-area-weighted gamma
         HET_N2O5(1) = HET_N2O5(1) + ADJUSTEDRATE         !loss rate coefficient

         IF (N.eq.8) THEN
              HET_N2O5(2) = HET_N2O5(2) + ( XSTKCF * SA ) !SA-weighted gamma
              SA_sum = SA_sum + SA                        !total SA in cm2/cm3
         ELSEIF (N.eq.10) THEN
              !don't include ORG SA since it has already been included above
              SA_sum = SA_sum
         ELSE
              HET_N2O5(2) = HET_N2O5(2) + ( XSTKCF * XAREA(N))!SA-weighted gamma
              SA_sum = SA_sum + XAREA(N)                      !total SA
         ENDIF

      END DO

      ! D) Divide gamma by SA to get SA-weighted gamma
      HET_N2O5(2) = safe_div( HET_N2O5(2), SA_sum, 0.0e+0_fp )

    END FUNCTION HETN2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N2O5_InorgOrg
!
! !DESCRIPTION: Function N2O5\_inorg\_org computes the GAMMA reaction probability
!     for N2O5 loss in inorganic (sulfate-nitrate-ammonium-sea salt) aerosols
!     with organic coatings, based on the recommendation of McDuffie (2018) JGR.
!     The inorganic core is based on Bertram and Thornton ACP (2009).
!\\
!\\
! !INTERFACE:
!
    FUNCTION N2O5_InorgOrg( volInorg, volOrg, H2Oinorg, H2Oorg, Rcore, &
                            NIT,      Cl,     T,        RH ) &
         RESULT ( output )
!
! !USES:
!
  USE PhysConstants,      ONLY : AVO, RGASLATM, RSTARG, PI
!
! !INPUT PARAMETERS:
!
      real(fp), intent(in) :: volInorg ! volume of wet inorganic aerosol core
                                       !  [cm3(aerosol)/cm3(air)]
      real(fp), intent(in) :: volOrg   ! volume of wet organic aerosol coating
                                       !  [cm3(aerosol)/cm3(air)]
      real(fp), intent(in) :: H2Oinorg ! volume of H2O in inorganic core
                                       !  [cm3(H2O)/cm3(air)]
      real(fp), intent(in) :: H2Oorg   ! volume of H2O in organic coating
                                       !  [cm3(H2O)/cm3(air)]
      real(fp), intent(in) :: Rcore    ! radius of inorganic core [cm]
      real(fp), intent(in) :: NIT      ! aerosol nitrate concentration
                                       !  [molecule/cm3(air)
      real(fp), intent(in) :: Cl       ! aerosol chloride concentration
                                       !  [molecule/cm3(air)
      real(fp), intent(in) :: T        ! air temperature [K]
      real(fp), intent(in) :: RH       ! relative humidity [%]
!
! !RETURN VALUE:
!
      real(fp) :: output(4) ! output(1) = gamma;
                            ! output(2) = ClNO2_yield
                            ! output(2) = particle radius, cm
                            ! output(4) = surface area, cm2/cm3
!
! !REMARKS:
!
! !REVISION HISTORY:
!    13 Dec 2018 - E. McDuffie - Initial version
!    23 Jul 2019 - C.D. Holmes - Consolidated into one function
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters from Bertram and Thornton (2009) ACP and McDuffie 2018 JGR
      real(fp), parameter :: KH    = 5.1e+1_fp   !unitless
      real(fp), parameter :: k3k2b = 4.0e-2_fp   !unitless
      real(fp), parameter :: beta  = 1.15e+6_fp  ![s-1]
      real(fp), parameter :: delta = 1.3e-1_fp   ![M-1]

      ! Organic Parameters from Antilla et al., 2006 and Riemer et al., 2009
      ! aq. Henry's Law coef [mol m-3 atm-1] (Antilla 2006)
      real(fp), parameter :: Haq  = 5e+3_fp
      ! aq. diffusion coef [m2 s-1] (Riemer, 2009)
      real(fp), parameter :: Daq  = 1e-9_fp
      ! Gas constant [m3 atm K-1 mol-1]
      real(fp), parameter :: R    = RGASLATM * 1e-3_fp
!
! !LOCAL VARIABLES
!
      real(fp)            :: k2f   ! [s-1]
      real(fp)            :: A     ! [s]
      real(fp)            :: speed ! mean molecular speed of N2O5 [m/s]

      real(fp) :: gamma, gamma_core, gamma_coat
      real(fp) :: volTotal, H2Ototal, areaTotal, volRatioDry
      real(fp) :: Rp, l, eps, OCratio
      real(fp) :: M_H2O, M_NIT, M_Cl, ClNO2_yield

      !------------------------------------------------------------------------
      ! Particle size & volume, coating thickness, molar concentrations
      !------------------------------------------------------------------------

      ! Total volume (organic + inorganic), cm3(aerosol)/cm3(air)
      volTotal = volInorg + volOrg

      ! Total H2O (organic + inorganic), cm3(H2O)/cm3(air)
      H2Ototal = H2Oinorg + H2Oorg

      ! Ratio of inorganic to total (organic+inorganic) volumes when dry, unitless
      volRatioDry = safe_div((volInorg - H2Oinorg), (volTotal - H2Ototal), 0e+0_fp)

      ! Particle radius, cm
      ! Derived from spherical geometry
      ! [note: The radius and surface area of a wet particle are
      ! properly calculated from the wet volume volume ratio (including water).
      ! We use the dry volume ratio here because McDuffie et al. (2018) fitted
      ! the N2O5 gamma parameters to field data in a model using the
      ! dry ratio. cdholmes 7/22/2019]
      Rp = safe_div( Rcore, volRatioDry**(1e+0_fp/3e+0_fp), Rcore)

      ! Coating thickness, cm
      l = Rp - Rcore

      ! mean molecular speed [m s-1]
      ! sqrt( 8RT / (pi M) )
      speed = sqrt( 8e+0_fp * RSTARG * T / ( PI * 0.108+0_fp) )

      ! H2O molar concentration, mol/L
      M_H2O = H2Ototal / 18e+0_fp / volTotal * 1e+3_fp
      ! Nitrate molar concentration, mol/L
      M_NIT = NIT / volTotal / AVO * 1e+3_fp
      ! Chloride molar concentration, mol/L
      M_Cl = Cl   / volTotal / AVO * 1e+3_fp

      !------------------------------------------------------------------------
      ! Gamma for the organic shell
      ! Implements recommendations by McDuffie (2018) JGR,
      !------------------------------------------------------------------------

      !O:C ratio from Eq. 10 of Canagaratna et al., 2015 (ACP)
      ! Take average OM/OC ratio from /GeosCore/aerosol_mod.F90
      OCratio = ( ( ( OMOC_POA + OMOC_OPOA ) / 2 ) - 1.17e+0_fp ) / 1.29e+0_fp

      ! organic scaling factor (eps(Haq*Daq) = Horg*Dorg)
      ! from McDuffie (2018) JGR
      eps = 1.5e-1_fp * OCratio + 1.6e-3_fp * RH

      ! Gamma for coating
      ! [Rcore, Rp, and l converted cm -> m here]
      IF ( l <= 0.0e+0_fp ) THEN
         gamma_coat = 0.0e+0_fp
      ELSE
         gamma_coat = ( 4e+0_fp * R * T * eps * Haq * Daq * Rcore/1e+2_fp ) / &
                      ( speed * l/1e+2_fp * Rp/1e+2_fp )
      ENDIF

      ! Total particle surface area, cm2/cm3
      areaTotal = 3e+0_fp * volTotal / Rp

      !------------------------------------------------------------------------
      ! Gamma for the inorganic core
      ! Implements recommendations by McDuffie (2018) JGR,
      ! following the general approach from Bertram and Thornton ACP (2009).
      !------------------------------------------------------------------------

      ! Select dry or deliquesed aerosol based on molar concentration of H2O
      IF ( M_H2O < 0.1e+0_fp ) THEN

         ! When H2O is nearly zero, use dry aerosol value
         gamma_core = 0.005e+0_fp

      ELSE

         ! mean molecular speed [cm/s]
         speed = speed * 1e+2_fp

         ! A factor from Bertram and Thornton (2009), s
         ! Their paper suggested an approximated value of A = 3.2D-8
         A = ( ( 4 * volTotal ) / ( speed * areaTotal ) ) * KH

         ! k2f - reaction rate constant of N2O5 with H2O
         ! From McDuffie (2018): k2f = 2.14D5 * H2O
         ! This linear water dependence is not accurate at large
         ! (>20 M) aerosol water concentrations. Therefore, k2f is
         ! calculated following Bertram and Thornton ACP (2009).
         ! Eq 11 from Bertram and Thronton (2009):
         ! Modified to avoid underflow when exp(-delta*H2O) ~1
         IF ( delta * M_H2O < 1e-2_fp ) THEN
            k2f = beta * ( delta * M_H2O )
         ELSE
            k2f = beta * ( 1e+0_fp - exp( -delta * M_H2O ) )
         ENDIF

         ! Eq 12 from Bertram and Thornton (2009)
         ! Use safe_div to avoid overflow when NIT ~ 0
         gamma_core = A * k2f * ( 1e+0_fp - &
            1e+0_fp / ( 1e+0_fp + safe_div( k3k2b * M_H2O, M_NIT, 1e+30_fp ) ) )

      ENDIF

      !------------------------------------------------------------------------
      ! Gamma for overall uptake
      !------------------------------------------------------------------------

      IF (gamma_coat <= 0.0e+0_fp ) THEN
         gamma = gamma_core
      ELSEIF (gamma_core <= 0.0e+0_fp ) THEN
         gamma = 0.0e+0_fp
      ELSE
         gamma = 1e+0_fp / (( 1e+0_fp / gamma_core ) + ( 1e+0_fp / gamma_coat ))
      ENDIF

      !------------------------------------------------------------------------
      ! ClNO2 yield
      !------------------------------------------------------------------------

      ! Calculate the ClNO2 yield following Bertram and Thornton 2009 ACP
      ! Reduce by 74% following McDuffie (2018) JGR
      ClNO2_yield = ClNO2_BT( M_Cl, M_H2O ) * 0.25e+0_fp

      output(1) = gamma
      output(2) = ClNO2_yield
      output(3) = Rp
      output(4) = areaTotal

    END FUNCTION N2O5_InorgOrg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ClNO2_BT
!
! !DESCRIPTION: Function ClNO2\_BT computes the PHI production yield of ClNO2
!  from N2O5 loss in sulfate-nitrate-ammonium (SNA) aerosols based on the
!  recommendation of Bertram and Thornton (2009) ACP.
!  In this implementation, we assume that SNA and sea salt aerosol are
!  externally mixed, so [Cl-] = 10% SALA in SNA. We do this because the
!  surface area of sea salt (coarse and fine) is calculated separately from
!  the surface area of SNA.
!\\
!\\
! !INTERFACE:
!
    FUNCTION ClNO2_BT( Cl, H2O ) RESULT( PHI )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: Cl         !Aerosol chloride, mol/L
      REAL(fp), INTENT(IN) :: H2O        ! Aerosol H2O, mol/L

!
! !RETURN VALUE:
!
      REAL(fp)             :: PHI
!
! !REMARKS:
!
! !REVISION HISTORY:
!    13 Dec 2018 - E. McDuffie - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters from Bertram and Thornton (2009) ACP
      REAL(fp), parameter :: k2k3  = 1e+0_fp / 4.5e+2_fp

      ! Initialize
      PHI     = 0.0_fp

      !When H2O is nearly zero, assign phi accordingly
      IF ( H2O < 0.1 ) THEN
           IF ( Cl > 1e-3_fp ) THEN
               PHI = 1e+0_fp
           ELSE
               PHI = 0e+0_fp
           ENDIF
      ELSE
           ! Eq from Bertram and Thronton (2009)
           ! Use safe_div to avoid overflow when Cl ~ 0
           PHI = 1e+0_fp / ( 1e+0_fp + k2k3 * ( safe_div( H2O, Cl, 0e+0_fp ) ) )

     ENDIF

    END FUNCTION ClNO2_BT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetGLYX
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for GLYX.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETGLYX( A, B ) RESULT( HET_GLYX )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_GLYX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  02 Mar 2018 - M. Sulprizio- Change daytime gamma to 4.4e-3 and nighttime
!                              gamma to 8.0e-6 based on recommendation from E.
!                              Marais to address that SOAGX is a factor of 1.5
!                              lower in v11-02d than in Marais et al. [2016]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_GLYX     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         ! Uptake during the day (higher uptake than night)
         ! (Sumner et al., 2014):
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define gamma for GLYX:
            IF ( SUNCOS .gt. 0 ) THEN

               ! Uptake during the day (use Liggio et al., 2005):
               ! XSTKCF = 2.9e-3_fp ! Prior to 3/2/18
               XSTKCF = 4.4e-3_fp

            ELSE

               ! Uptake at night (lower uptake than day)
               ! Value is within the range 1d-5 to 1d-6
               ! (Faye McNeill personal communication, eam, 2015):
               ! XSTKCF = 5.0e-6_fp ! Prior to 3/2/18
               XSTKCF = 8.0e-6_fp

            ENDIF

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_GLYX = HET_GLYX + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETGLYX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMGLY
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MGLY.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMGLY( A, B ) RESULT( HET_MGLY )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MGLY
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MGLY     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define gamma for MGLY:
            ! Obtained by scaling gamma GLYX by the
            ! ratio of effective Henry's law constants
            ! for GLYX (3d7) and MGLY (3.7d3) (eam, 02/2015):
            XSTKCF = 3.6e-7_fp

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MGLY = HET_MGLY + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMGLY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIEPOX
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IEPOX.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIEPOX( A, B ) RESULT( HET_IEPOX )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IEPOX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      REAL(fp) :: HSTAR, K_HPLUS, K_NUC, K_HSO4, K_HYDRO

      ! Initialize
      HET_IEPOX    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define Henry's Law constant
            ! Changes H* for IEPOX again to accommodate
            ! reduction in yields of RIP, precursor
            ! of IEPOX (eam, 07/2015):
            HSTAR = HSTAR_EPOX    ! (Nguyen et al., 2014)

            ! Define first-order particle phase reaction rates
            ! specific to IEPOX (from Gaston et al., 2014):
            K_HPLUS = 3.6e-2_fp   ! Alternate: 1.2d-3 (Edding)
            K_NUC   = 2.e-4_fp    ! Alternate: 5.2d-1 (Piletic)
            K_HSO4  = 7.3e-4_fp
            K_HYDRO = 0.0e+0_fp

            ! Get GAMMA for IEPOX hydrolysis:
            XSTKCF = EPOXUPTK( XAREA(N), XRADI(N),            &
                               TEMPK,    (A**0.5_fp),         &
                               HSTAR,    K_HPLUS,    H_PLUS,  &
                               K_NUC,    MSO4,       MNO3,    &
                               K_HSO4,   MHSO4,      K_HYDRO )

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_IEPOX = HET_IEPOX + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETIEPOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIMAE
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IMAE.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIMAE( A, B ) RESULT( HET_IMAE )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IMAE
!
! !REMARKS:
! Here use the same values as are read in for IEPOX, but scale down gamma by
! a factor of 30 to get the value for IMAE. Gamma for the two species are
! similar under neutral conditions, so use IEPOX gamma when [H+] <= 8d-5.
! Implemented by (eam, 01/2015) using lab study findings from Riedel et al.,
! EST, 2015.
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      REAL(fp) :: HSTAR, K_HPLUS, K_NUC, K_HSO4, K_HYDRO

      ! Initialize
      HET_IMAE     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define Henry's Law constant.
            ! Changes H* for IEPOX again to accommodate
            ! reduction in yields of RIP, precursor
            ! of IEPOX (eam, 07/2015):
            HSTAR = HSTAR_EPOX   ! (Nguyen et al., 2014)

            ! Define first-order particle phase reaction rates
            ! specific to IEPOX (from Gaston et al., 2014):
            K_HPLUS = 3.6e-2_fp   ! Alternate: 1.2d-3 (Edding)
            K_NUC   = 2.6e-4_fp   ! Alternate: 5.2d-1 (Piletic)
            K_HSO4  = 7.3e-4_fp
            K_HYDRO = 0.e+0_fp

            ! Get GAMMA for IMAE hydrolysis:
            XSTKCF = EPOXUPTK( XAREA(N), XRADI(N),            &
                               TEMPK,    (A**0.5_fp),         &
                               HSTAR,    K_HPLUS,    H_PLUS,  &
                               K_NUC,    MSO4,       MNO3,    &
                               K_HSO4,   MHSO4,      K_HYDRO )

            ! Scale down gamma if H+ > 8d-5 (30x less than gamma for IEPOX)
            ! (Riedel et al., 2015)
            IF ( H_PLUS .gt. 8.e-5_fp ) THEN
               XSTKCF = XSTKCF / 30.e+0_fp
            ENDIF

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_IMAE = HET_IMAE + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETIMAE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetLVOC
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for LVOC: condensation of
! low-volatility ISOPOOH oxidation products.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETLVOC( A, B ) RESULT( HET_LVOC )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_LVOC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_LVOC     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_LVOC = HET_LVOC + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETLVOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISN1OG
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISN1OG: uptake of 2nd
! generation organic nitrates formed from ISOP+NO3 reaction (eam, 02/2015).
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISN1OG( A, B ) RESULT( HET_ISN1OG )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISN1OG
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISN1OG   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISN1OG = HET_ISN1OG + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISN1OG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISOPND
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISOPND.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISOPND( A, B ) RESULT( HET_ISOPND )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISOPND
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISOPND   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISOPND = HET_ISOPND + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISOPND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISOPNB
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISOPNB.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISOPNB( A, B ) RESULT( HET_ISOPNB )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISOPNB
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISOPNB   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISOPNB = HET_ISOPNB + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISOPNB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMACRN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MACRN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMACRN( A, B ) RESULT( HET_MACRN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MACRN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MACRN    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MACRN = HET_MACRN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMACRN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMVKN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MVKN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMVKN( A, B ) RESULT( HET_MVKN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MVKN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MVKN     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MVKN = HET_MVKN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMVKN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetR4N2
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for R4N2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETR4N2( A, B ) RESULT( HET_R4N2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_R4N2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_R4N2   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_R4N2 = HET_R4N2 + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETR4N2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetDHDN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for DHDN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETDHDN( A, B ) RESULT( HET_DHDN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_DHDN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  14 Jul 2017 - M. Sulprizio- Product has been changed to IONITA, which also
!                              has heterogeneous reaction. Remove call to
!                              EPOXUPTK here and use gamma value specified in
!                              SET_HET (Fisher et al., 2016).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_DHDN     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_DHDN = HET_DHDN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETDHDN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITS
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITS
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITS( A, B ) RESULT( HET_MONITS )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITS
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITS   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MONITS = HET_MONITS + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMONITS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITU
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITU
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITU( A, B ) RESULT( HET_MONITU )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITU
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITU   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MONITU = HET_MONITU + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMONITU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHONIT
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for HONIT
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHONIT( A, B ) RESULT( HET_HONIT )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HONIT
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HONIT    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         XSTKCF = B

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HONIT = HET_HONIT + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETHONIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIONITA
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IONITA: Aerosol-phase
!  organic nitrate formed from monoterpene precursors.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIONITA( A, B ) RESULT( HET_IONITA )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IONITA
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_IONITA   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma for IONITA
         ! Imposed lifetime = 1 hour (Fisher et al., 2016)
         XSTKCF = 2.78e-4_fp

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_IONITA = HET_IONITA + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETIONITA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITA
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITA: Aerosol-phase
!  organic nitrate formed from monoterpene precursors.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITA( A, B ) RESULT( HET_MONITA )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITA
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITA   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma for MONITA
         ! Imposed lifetime = 1 hour (Fisher et al., 2016)
         XSTKCF = 2.78e-4_fp

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MONITA = HET_MONITA + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMONITA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5_SS
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5 on sea salt. This
!  reaction follows the N2O5 + Cl- channel, and Cl- is assumed to be in excess.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5_SS( A, B ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum(3) !(1) = rate coefficient, (2) = gamma (SS), (3) = phi
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  13 Dec 2018 - E. McDuffie - Update to report ClNO2 yield & N2O5 gamma
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, SA_total, ClNO2_yield, kClNO2
!
! !DEFINED PARAMETERS:
!
      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      SA_total     = 0.0_fp
      ClNO2_yield  = 0.0_fp
      kClNO2       = 0.0_fp

      ! Directly calculate for sea salt only
      ! Get GAMMA for N2O5 hydrolysis, which is
      ! a function of aerosol type, temp, and RH
      DO N=8, NAEROTYPE
         IF ((N.eq.11).or.(N.eq.12)) THEN
            ! Sea salt - follows the N2O5 + Cl- channel
            XSTKCF = N2O5( N, TEMPK, RELHUM )

            ! Convert to first-order rate constant
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))

            ! Add to overall reaction rate
            kISum(1) = kISum(1) + ADJUSTEDRATE
            kISum(2) = kISum(2) + ( XSTKCF * XAREA(N) )
            SA_total = SA_total + XAREA(N)
         ELSEIF (N.eq.8) THEN
            ClNO2_yield = 0 ! (add in future update)
            ! phi = kClNO2/kN2O5 (kN2O5 = sum (ClNO2+HNO3) and 2*HNO3 production pathways)
            kClNO2 = ClNO2_yield * ADJUSTEDRATE
            ! rate of this HNO3 + ClNO2 pathway for this aerosol type is kClNO2
            !Calculate the total rate of this pathway from both SS and SNA aerosol
            ADJUSTEDRATE = kClNO2
            kISum(1) = kISum(1) + ADJUSTEDRATE
         ENDIF
      END DO

      ! Only report gamma for SS aerosol here
      ! SA-weighted gamma from other aerosol types calc'd by HETN2O5()
      ! This is the only place where ClNO2_yield is recorded
      kISum(2) = safe_div( kISum(2),  SA_total, 0.0e+0_fp )

    END FUNCTION HETN2O5_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5(g) + HCl(l,s)
!  in polar stratospheric clouds and on tropospheric sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum(2) !(1) = rate coefficient, (2) = gamma (SS)
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  04 May 2016 - M. Sulprizio- Add fixes for setting rate if not a STRATBOX
!  24 Dec 2016 - S. D. Eastham - Extended into the troposphere. Also now use the
!                              standard N2O5 calculation to establish gamma for
!                              sulfate, rather than relying on a fixed factor.
!  03 Jan 2018 - M. Sulprizio  - Replace UCX CPP switch with Input_Opt%LUCX
!  13 Dec 2018 - E. McDuffie - Now report kN2O5 and SA-weighted gamma
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, SA_total

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      SA_total     = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Assume zero
         XStkCf = 0.0e+0_fp

	 ! restore stratosphere only limitation - TMS 17/04/10
         IF ( STRATBOX ) THEN
            IF (N.eq.8) THEN
               ! Fixed gamma?
               !XSTKCF = 0.1e-4_fp ! Sulfate
               ! RH dependence
               ! Note gamma calculation on sulfate in troposphere uses
               ! the McDuffie parameterization (func: HetN2O5())
      	       XSTKCF = N2O5( N, TEMPK, RELHUM )
	    ENDIF
         ENDIF

         ! For UCX-based mechanisms only consider PSC reactions in strat
         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
            IF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(2)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.003e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.03e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            ! EEM: update to report surface area weighted gamma and rate
            kISum(1) = kISum(1) + ADJUSTEDRATE
            kISum(2) = kISum(2) + ( XSTKCF * XAREA(N) )
            SA_total = SA_total + XAREA(N)
         ENDIF

      END DO

      !Calculate SA-weighted gamma at end
      IF ( SA_total > 0.0e+0_fp ) THEN
         kISum(2) = kISum(2) / SA_total
      ELSE
         kISum(2) = 0.0e+0_fp
      ENDIF

    END FUNCTION HETN2O5_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHXUptake_js
!
! !DESCRIPTION: Sets the uptake rate of HCl and HBr on sea salt using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHXUptake_JS( denAir, rAer, AAer, TK, X ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      INTEGER,  INTENT(IN) :: X           ! 1: Cl-, 2: Br-
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, XSqM
      Real(fp), Parameter :: XMolWeightHCl=36.5e+0_fp
      Real(fp), Parameter :: XSqMHCl=SQRT(XMolWeightHCl)
      Real(fp), Parameter :: XMolWeightHBr=81.0e+0_fp
      Real(fp), Parameter :: XSqMHBr=SQRT(XMolWeightHBr)

      ! Initialize
      kISum        = 0.0_fp

      ! Select between halogens
      IF (X.eq.1) THEN
         XSqM = XSqMHCl
      ELSEIF (X.eq.2) THEN
         XSqM = XSqMHBr
      ENDIF

      XStkCf = Gamma_HX_Uptake( rAer, denAir, X, TK )

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,XTemp,XSqM)

    END FUNCTION HETHXUptake_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetO3_SS_JS
!
! !DESCRIPTION: Sets the O3 + Br- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETO3_SS_JS( denAir, rAer, AAer, alkAer, TK, halConc, O3Conc ) &
                             RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: alkAer      ! Aerosol alkalinity (?)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: halConc     ! Halide concentration (mol/L)
      REAL(fp), INTENT(IN) :: O3Conc      ! Ozone concentration (#/cm3)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=48.0e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      kISum        = 0.0_fp

      ! Reaction can only proceed on acidic aerosol
      IF (alkAer > 0.05e+0_fp) THEN
         XStkCf = 0.0e+0_fp
      ELSE
         XStkCf = Gamma_O3_Br( rAer, denAir, TK, halConc, O3Conc )
      ENDIF

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,XTemp,XSqM)

    END FUNCTION HETO3_SS_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_SS_JS
!
! !DESCRIPTION: Sets the ClNO3 + Br- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_SS_JS( denAir, rAer, AAer, alkAer, TK, halConc ) &
                             RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: alkAer      ! Aerosol alkalinity (?)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: halConc     ! Halide concentration (mol/L)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=97.5e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      kISum        = 0.0_fp

      ! Reaction can only proceed on acidic aerosol
      IF (alkAer > 0.05e+0_fp) THEN
         XStkCf = 0.0e+0_fp
      ELSE
         XStkCf = Gamma_ClNO3_Br( rAer, denAir, TK, halConc )
      ENDIF

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,XTemp,XSqM)

    END FUNCTION HETClNO3_SS_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_SS_JS
!
! !DESCRIPTION: Sets the HOBr + Br- or Cl- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_SS_JS( denAir, rAer, AAer, alkAer, TK, hConc, clConc, &
                            brConc, X ) &
                            RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: alkAer      ! Aerosol alkalinity (?)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: hConc       ! H+ concentration (mol/L)
      REAL(fp), INTENT(IN) :: clConc      ! Cloride concentration (mol/L)
      REAL(fp), INTENT(IN) :: brConc      ! Bromide concentration (mol/L)
      Integer,  INTENT(IN) :: X           ! 1: Cl-, 2: Br-
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!  01 Dec 2017 - Q.J. Chen     - Updated to account for Cl- and Br- separately;
!                                Now calls routine Gamma_HOBr_AER
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=96.9e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp

      ! Reaction can only proceed on acidic aerosol
      IF (alkAer > 0.05e+0_fp) THEN
         XStkCf = 0.0e+0_fp
         r_gp   = 1.0e+0_fp
      ELSE
!-----------------------------------------------------------------------------
! Prior to 7/17/18:
! Xuan Wang says to replace 2 with X in the call to GAMMA_HOBR_AER
! (bmy, 7/17/18)
!         CALL Gamma_HOBr_AER(rAer, denAir, 2, TK, clConc, brConc, &
!-----------------------------------------------------------------------------
         CALL Gamma_HOBr_AER(rAer, denAir, X, TK, clConc, brConc, &
                             hConc, GAM_HOBr, r_gp)
         XStkCf = GAM_HOBr
      ENDIF

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,XTemp,XSqM)*r_gp

    END FUNCTION HETHOBr_SS_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetO3_HBr_JS
!
! !DESCRIPTION: Sets the O3 + Br- rate using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETO3_HBr_JS( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, brConc, O3Conc ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: brConc      ! Bromide concentration (mol/L)
      REAL(fp), INTENT(IN) :: O3Conc      ! Ozone concentration (mol/L)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=48.0e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         XStkCf = 0.0e+0_fp
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            XSTKCF = Gamma_O3_Br( xRadi(8), denAir, TK, brConc, O3Conc )
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF
      END DO

    ! Reaction on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       IF (ALiq.gt.0.0e+0_fp) THEN
          XStkCf = Gamma_O3_Br( rLiq, denAir, TK, brConc, O3Conc )
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)
       ENDIF
       IF (AIce.gt.0.0e+0_fp) THEN
          XStkCf = Gamma_O3_Br( rIce, denAir, TK, brConc, O3Conc )
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)
       ENDIF
    ENDIF

    END FUNCTION HETO3_HBr_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_O3_Br
!
! !DESCRIPTION: Function GAMMA\_O3\_Br calculates reactive uptake coef.
!               for bromide oxidation by O3
!\\
!\\
! !INTERFACE:
!
      FUNCTION GAMMA_O3_Br( Radius, n_air, T, C_Y, C_X_g ) RESULT( GAM )
!
! !USES:
!
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                         :: GAM
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)             :: Radius, n_air
      REAL(fp), INTENT(IN)             :: T, C_Y, C_X_g
!
! !REVISION HISTORY:
!  24 Sep 2015 - J. Schmidt  - Initial version
!  27 Feb 2018 - M. Sulprizio- Obtain Henry's law parameters from species
!                              database in SET_HET instead of hardcoding here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)       :: ab, gb, gd, gs, M_X
      REAL(fp)       :: cavg, H_X
      REAL(fp)       :: KLangC, k_s, C_Y_surf, Nmax
      REAL(fp)       :: k_b, D_l, l_r

      H_X = H_K0_O3*dexp(H_CR_O3*(1.0e0_fp/T - 1.0e0_fp/H_O3_T))
      M_X = MW_O3 * 1e-3_fp

      cavg    = dsqrt(8*RStarG*T/(Pi*M_X)) *1.0e2_fp ! thermal velocity (cm/s)

      Nmax = 3.0e14_fp ! #/cm2
      KLangC = 1.0e-13_fp !cm3
      k_s = 1.0e-16_fp !cm2s-1, from ks*Nmax=0.03s-1
      C_Y_surf= min(3.41e14_fp*C_Y, Nmax) ! [Br-(surf)] = 3.41E14 cm-2/M * [Br-(bulk)], but not gt Nmax.
      gs = (4.0e0_fp * k_s * C_Y_surf * KLangC * Nmax) / &
                    (cavg * (1.0e0_fp + KLangC * C_X_g) )

      k_b = 6.3e8_fp *  dexp(-4.45e3_fp / T) !M-1 s-1
      D_l = 8.9e-6_fp !cm2 s-1.
      l_r = dsqrt( D_l / (k_b * C_Y ) )! cm
      gb  = 4.0e0_fp * H_X * con_R * T * l_r * k_b * C_Y / cavg
      gb  = gb * REACTODIFF_CORR( Radius, l_r)

      GAM = gb + gs

      END FUNCTION GAMMA_O3_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_HBr_JS
!
! !DESCRIPTION: Sets the ClNO3 + Br- rate using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_HBr_JS( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, &
                              brConc, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: brConc      ! Bromide concentration (mol/L)
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  22 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!  03 Jan 2018 - M. Sulprizio  - Replace UCX CPP switch with Input_Opt%LUCX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=97.5e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         XStkCf = 0.0e+0_fp

         IF ( N == 8 ) THEN

            ! sulfate aerosol
            XSTKCF = Gamma_ClNO3_Br( xRadi(8), denAir, TK, brConc )

         ELSEIF ( Input_Opt%LUCX .and. STRATBOX ) THEN

            ! For UCX-based mechanisms only consider PSC reactions in strat
            IF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(5)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ENDIF

         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF
      END DO

      ! Reaction on liquid and ice clouds (tropospheric only)
      IF (.not. StratBox) THEN
         IF (ALiq.gt.0.0e+0_fp) THEN
          XStkCf = Gamma_ClNO3_Br( rLiq, denAir, TK, brConc )
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)
         ENDIF
         IF (AIce.gt.0.0e+0_fp) THEN
          XStkCf = Gamma_ClNO3_Br( rIce, denAir, TK, brConc )
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)
         ENDIF
      ENDIF

    END FUNCTION HETClNO3_HBr_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_JS
!
! !DESCRIPTION: Sets the hydrolysis rate for ClNO3 using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_JS( denAir, rLiq, rIce, ALiq, AIce, TK ) RESULT( HET_ClNO3 )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ClNO3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  16 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=97.5e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      HET_ClNO3    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Only apply PSC rate adjustment if at high altitude
      DO_EDUCT     = STRATBOX

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            XSTKCF = 0.024e+0_fp
         ELSEIF ( (N == 11) .OR. ( N == 12) ) THEN
            ! 2 modes of sea-salt
            XSTKCF = 0.024e+0_fp
         ELSEIF (N.eq.13) THEN
            XSTKCF = KHETI_SLA(3)
         ELSEIF (N.eq.14) THEN
            IF (NATSURFACE) THEN
               XSTKCF = 0.004e+0_fp ! NAT
            ELSE
               XSTKCF = 0.3e+0_fp ! Ice
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         HET_ClNO3 = HET_ClNO3 + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       HET_ClNO3 = HET_ClNO3 + Cld1K_XNO3(denAir,TK,rLiq,rIce,ALiq,AIce,XMolWeight,2.4E-2_fp)
    ENDIF

    END FUNCTION HETClNO3_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_HCl_JS
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Cl- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_HCl_JS( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc, brConc, &
                           hso3Conc, so3Conc ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp), INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp), INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp), INTENT(IN) :: clConc      ! Chloride concentration (mol/L)
      REAL(fp), INTENT(IN) :: brConc      ! Bromide  concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L), qjc
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Generated code based on Johan Schmidt's work
!  01 Dec 2017 - Q.J. Chen     - Updated to account for Br-, HSO3-, and SO3--;
!                                Now calls routine Gamma_HOBr_AER
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=96.9e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            CALL Gamma_HOBr_AER(xRadi(8), denAir, 1, TK, clConc, brConc, &
                                hConc_Sul, GAM_HOBr, r_gp)
            XSTKCF = GAM_HOBr
         ELSEIF (N.eq.13) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(10)
         ELSEIF (N.eq.14) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.1e+0_fp ! NAT
            ELSE
               XSTKCF = 0.3e+0_fp ! Ice
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            CALL Gamma_HOBr_AER(xRadi(8), denAir, 1, TK, clConc, brConc, &
                                hConc_Sul, GAM_HOBr, r_gp)
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)*r_gp

         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       IF (ALiq.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rLiq, denAir, 1, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_LCl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
       IF (AIce.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rIce, denAir, 1, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_ICl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
    ENDIF

    END FUNCTION HETHOBr_HCl_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_HBr_JS
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Br- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_HBr_JS( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc, brConc, &
                           hso3Conc, so3Conc ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp), INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp), INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp), INTENT(IN) :: clConc      ! Chloride concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: brConc      ! Bromide  concentration (mol/L)
      REAL(fp), INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L), qjc
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Generated code based on Johan Schmidt's work
!  01 Dec 2017 - Q.J. Chen     - Updated to account for Cl-, HSO3-, and SO3--;
!                                Now calls routine Gamma_HOBr_AER
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=96.9e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
      Real(fp) :: SADen
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            CALL Gamma_HOBr_AER(xRadi(8), denAir, 2, TK, clConc, brConc, &
                                hConc_Sul, GAM_HOBr, r_gp)
            XSTKCF = GAM_HOBr
         ELSEIF ( N == 13 ) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(6)
         ELSEIF ( N == 14 ) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.001e+0_fp
            ELSE
               XSTKCF = 0.3e+0_fp
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            CALL Gamma_HOBr_AER(xRadi(8), denAir, 2, TK, clConc, brConc, &
                                hConc_Sul, GAM_HOBr, r_gp)
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTemp,XSQM)*r_gp
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       IF (ALiq.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rLiq, denAir, 2, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_LCl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
       IF (AIce.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rIce, denAir, 2, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_ICl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
    ENDIF

    END FUNCTION HETHOBr_HBr_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_HSO3
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + HSO3- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_HSO3( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc, brConc, &
                           hso3Conc, so3Conc ) &
                           RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp), INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp), INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp), INTENT(IN) :: clConc      ! Chloride concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: brConc      ! Bromide  concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L)
      REAL(fp), INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L), qjc
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Nov 2017 - M. Sulprizio  - Generated code based on Qianjie Chen's work
!  01 Dec 2017 - Q.J. Chen     - Updated to account for Cl-, Br-, and SO3--;
!                                Now calls routine Gamma_HOBr_CLD
!                                Remove HOBr+S(IV) on sulfate aerosols
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=96.9e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! No HOBr+S(IV) on sulfate aerosols
            XSTKCF = 0e+0_fp
         ELSEIF (N.eq.13) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(10)
         ELSEIF (N.eq.14) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.1e+0_fp ! NAT
            ELSE
               XSTKCF = 0.3e+0_fp ! Ice
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       IF (ALiq.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rLiq, denAir, 3, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_LCl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
       IF (AIce.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rIce, denAir, 3, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_ICl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
    ENDIF

    END FUNCTION HETHOBr_HSO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_SO3
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + SO3-- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_SO3( denAir, rLiq, rIce, ALiq, AIce, VAir, TK, &
                          hConc_Sul, hConc_LCl, hConc_ICl, clConc, brConc, &
                          hso3Conc, so3Conc ) &
                          RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp), INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp), INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp), INTENT(IN) :: clConc      ! Chloride concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: brConc      ! Bromide  concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L), qjc
      REAL(fp), INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Nov 2017 - M. Sulprizio  - Generated code based on Qianjie Chen's work
!  01 Dec 2017 - Q.J. Chen     - Updated to account for Cl-, Br-, and HSO3-;
!                                Now calls routine Gamma_HOBr_CLD;
!                                Remove HOBr+S(IV) on sulfate aerosols
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=96.9e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! No HOBr+S(IV) on sulfate aerosols
            XSTKCF = 0e+0_fp
         ELSEIF (N.eq.13) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(10)
         ELSEIF (N.eq.14) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.1e+0_fp ! NAT
            ELSE
               XSTKCF = 0.3e+0_fp ! Ice
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       IF (ALiq.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rLiq, denAir, 4, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_LCl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(ALiq, rLiq, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
       IF (AIce.gt.0.0e+0_fp) THEN
          CALL Gamma_HOBr_CLD(rIce, denAir, 4, TK, clConc, brConc, &
                              hso3Conc, so3Conc, hConc_ICl, GAM_HOBr, r_gp)
          XSTKCF = GAM_HOBr
          kISum = kISum + Arsl1K(AIce, rIce, denAir, XStkCf, XTemp, XSqM)*r_gp
       ENDIF
    ENDIF

    END FUNCTION HETHOBr_SO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_HX_Uptake
!
! !DESCRIPTION: Function GAMMA\_HX\_uptake calculates mass accomidation coef.
!               for uptake of HX (HCl or HBr)
!\\
!\\
! !INTERFACE:
!
      FUNCTION Gamma_HX_Uptake( Radius, n_air, X, T )  RESULT( GAM )
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                         :: GAM
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)             :: Radius, n_air
      INTEGER, INTENT(IN)              :: X
      REAL(fp), INTENT(IN)             :: T
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      REAL(fp)       :: ab

      ! 1: Cl-, 2: Br-
      IF (X==1) THEN
         ab = 4.4e-6_fp * dexp( 2898.0e0_fp / T ) ! ab(RT) = 0.069
      ELSE
         ab = 1.3e-8_fp * dexp( 4290.0e0_fp / T ) ! ab(RT) = 0.021
      ENDIF

      GAM = ab

      END FUNCTION Gamma_HX_Uptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_HOBr_X
!
! !DESCRIPTION: Function GAMMA\_HOBr\_X calculates reactive uptake coef.
!               for halide (Cl- and Br-) and S(IV) (HSO3- and SO3--)
!               oxidation by HOBr
!\\
!\\
! !INTERFACE:
!
      FUNCTION GAMMA_HOBr_X( Radius, n_air, X, T, C_Y, C_Hp ) RESULT( GAM )
!
! !OUTPUT PARAMETERS:
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                       :: GAM
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)           :: Radius   ! Radius (cm)
      REAL(fp), INTENT(IN)           :: n_air    ! n_air (#/cm)
      INTEGER,  INTENT(IN)           :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(fp), INTENT(IN)           :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)           :: C_Y      ! Concentration (mol/L)
      REAL(fp), INTENT(IN)           :: C_Hp     ! Concentration (mol/L)
!
! !REVISION HISTORY:
!  24 Sep 2015 - J. Schmidt  - Initial version
!  15 Nov 2017 - M. Sulprizio- Added options for HSO3- and SO3-- based on
!                              Qianjie Chen's work
!  27 Feb 2018 - M. Sulprizio- Obtain Henry's law parameters from species
!                              database in SET_HET instead of hardcoding here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!


      REAL(fp)       :: ab, gb, gd, M_X
      REAL(fp)       :: cavg, k_b, D_l, l_r, H_X

      IF ( X==1 ) THEN
         ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
         !k_b  = 5.9e+9_fp
         ! (Liu and Margerum, Environ. Sci. Tech., 2001)
         k_b  = 2.3e+10_fp ! (qjc, 12/28/16)
      ELSEIF ( X==2 ) THEN
         ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
         k_b  = 1.6e+10_fp
      ELSEIF ( X==3 ) THEN
         ! Reaction rate coefficient for HOBr + HSO3- [M-2 s-1]
         ! (Liu and Margerum, Environ. Sci. Tech., 2001)
         k_b  = 3.2e+9_fp
      ELSEIF ( X==4 ) THEN
         ! Reaction rate coefficient for HOBr + HSO3-- [M-2 s-1]
         ! (Troy and Margerum, Inorg. Chem., 1991)
         k_b  = 5.0e+9_fp
      ENDIF

      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l  = 1.4e-5_fp

      H_X = H_K0_HOBr*dexp(H_CR_HOBr*(1.0e0_fp/T - 1.0e0_fp/H_HOBr_T))
      M_X = MW_HOBr * 1e-3_fp

      ! Mass accommodation coefficient
      ab = 0.6e0_fp

      ! Thermal velocity [cm/s]
      cavg = dsqrt(8*RStarG*T/(pi*M_X)) *1.0e2_fp

      ! Diffusive length scale [cm]
      l_r  = dsqrt( D_l / (k_b * C_Y * C_Hp ) )

      ! Bulk reaction coefficient [unitless]
      gb = 4.0e0_fp * H_X * con_R * T * l_r * k_b * C_Y * C_Hp / cavg
      gb = gb * REACTODIFF_CORR( Radius, l_r)

      ! Reactive uptake coefficient [unitless]
      GAM = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb)

    END FUNCTION GAMMA_HOBr_X
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gamma_HOBr_cld
!
! !DESCRIPTION: Returns GAMMA for HOBr in clouds (need better description)
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_HOBr_CLD( Radius, n_air, X, T, C_Y1, C_Y2, &
                               C_Y3, C_Y4, C_Hp, GAM_HOBr, r_gp )
!
! !OUTPUT PARAMETERS:
      ! Reactive uptake coefficient (unitless)
      REAL(fp), INTENT(OUT)          :: GAM_HOBr, r_gp
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)           :: Radius   ! Radius (cm)
      REAL(fp), INTENT(IN)           :: n_air    ! n_air (#/cm)
      INTEGER,  INTENT(IN)           :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(fp), INTENT(IN)           :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)           :: C_Y1, C_Y2, C_Y3, C_Y4      ! Concentration (mol/L)
      REAL(fp), INTENT(IN)           :: C_Hp     ! Concentration (mol/L)
!
! !REVISION HISTORY:
!  30 Nov 2017 - Q.J. Chen   - Initial version
!  27 Feb 2018 - M. Sulprizio- Obtain Henry's law parameters from species
!                              database in SET_HET instead of hardcoding here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      REAL(fp)       :: ab, gd, M_X
      REAL(fp)       :: cavg, H_X
      REAL(fp)       :: gb1, gb2, gb3, gb4, gb_tot
      REAL(fp)       :: k_b1, k_b2, k_b3, k_b4
      REAL(fp)       :: D_l, ybr2
      REAL(fp)       :: l_r1, l_r2, l_r3, l_r4

!      IF ( X==1 ) THEN
         ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
         !k_b  = 5.9e+9_fp
         ! (Liu and Margerum, Environ. Sci. Tech., 2001)
         k_b1  = 2.3e+10_fp ! (qjc, 12/28/16)
!      ELSEIF ( X==2 ) THEN
         ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
         k_b2  = 1.6e+10_fp
!      ELSEIF ( X==3 ) THEN
         ! Reaction rate coefficient for HOBr + HSO3- [M-2 s-1]
         ! (Liu and Margerum, Environ. Sci. Tech., 2001)
         k_b3  = 3.2e+9_fp
!      ELSEIF ( X==4 ) THEN
         ! Reaction rate coefficient for HOBr + HSO3-- [M-2 s-1]
         ! (Troy and Margerum, Inorg. Chem., 1991)
         k_b4  = 5.0e+9_fp
!      ENDIF

      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l  = 1.4e-5_fp

      H_X = H_K0_HOBr*dexp(H_CR_HOBr*(1.0e0_fp/T - 1.0e0_fp/H_HOBr_T))
      M_X = MW_HOBr * 1e-3_fp

      ! Mass accommodation coefficient
      ab = 0.6e0_fp

      ! Thermal velocity [cm/s]
      cavg = dsqrt(8*RStarG*T/(pi*M_X)) *1.0e2_fp

      ! l_r is diffusive length scale [cm]; gb is Bulk reaction coefficient [unitless]
      l_r1  = dsqrt( D_l / (k_b1 * C_Y1 * C_Hp ) )
      gb1 = 4.0e0_fp * H_X * con_R * T * l_r1 * k_b1 * C_Y1 * C_Hp / cavg
      gb1 = gb1 * REACTODIFF_CORR( Radius, l_r1)

      l_r2  = dsqrt( D_l / (k_b2 * C_Y2 * C_Hp ) )
      gb2 = 4.0e0_fp * H_X * con_R * T * l_r2 * k_b2 * C_Y2 * C_Hp / cavg
      gb2 = gb2 * REACTODIFF_CORR( Radius, l_r2)

      l_r3  = dsqrt( D_l / (k_b3 * C_Y3 ) )
      gb3 = 4.0e0_fp * H_X * con_R * T * l_r3 * k_b3 * C_Y3 / cavg
      gb3 = gb3 * REACTODIFF_CORR( Radius, l_r3)

      l_r4  = dsqrt( D_l / (k_b4 * C_Y4 ) )
      gb4 = 4.0e0_fp * H_X * con_R * T * l_r4 * k_b4 * C_Y4 / cavg
      gb4 = gb4 * REACTODIFF_CORR( Radius, l_r4)

      gb_tot = gb1 + gb2 + gb3 + gb4

      ! Reactive uptake coefficient [unitless]
      GAM_HOBr = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb_tot)

      ybr2 = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
      ybr2 = MIN(ybr2, 0.9e0)
      ybr2 = MAX(ybr2, TINY(1.0e0))

      IF ( X==1 ) THEN

         r_gp = (gb1+gb2)/gb_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = TINY(1.0e0)
         ELSE
            r_gp = r_gp * (1.e0 - ybr2)
         ENDIF

      ELSEIF ( X==2 ) THEN

         r_gp = (gb1+gb2)/gb_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.9e0 * r_gp
         ELSE
            r_gp = r_gp * ybr2
         ENDIF

      ELSEIF ( X==3 ) THEN

         r_gp = gb3/gb_tot

      ELSEIF ( X==4 ) THEN

         r_gp = gb4/gb_tot

      ENDIF

    END SUBROUTINE GAMMA_HOBr_CLD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gamma_HOBr_aer
!
! !DESCRIPTION: Returns GAMMA for HOBr in aerosol? (need better description)
!\\
!\\!
! !INTERFACE:
!
    SUBROUTINE GAMMA_HOBr_AER( Radius, n_air, X, T, C_Y1, C_Y2, &
                               C_Hp, GAM_HOBr, r_gp )
!
! !USES:
!
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(fp), INTENT(OUT)          :: GAM_HOBr, r_gp
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)           :: Radius   ! Radius (cm)
      REAL(fp), INTENT(IN)           :: n_air    ! n_air (#/cm)
      INTEGER,  INTENT(IN)           :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(fp), INTENT(IN)           :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)           :: C_Y1, C_Y2      ! Concentration (mol/L)
      REAL(fp), INTENT(IN)           :: C_Hp     ! Concentration (mol/L)
!
! !REVISION HISTORY:
!  30 Nov 2017 - Q.J. Chen   - Initial version
!  27 Feb 2018 - M. Sulprizio- Obtain Henry's law parameters from species
!                              database in SET_HET instead of hardcoding here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)       :: ab, gd, M_X
      REAL(fp)       :: cavg, H_X
      REAL(fp)       :: gb1, gb2, gb_tot
      REAL(fp)       :: k_b1, k_b2
      REAL(fp)       :: D_l, ybr2
      REAL(fp)       :: l_r1, l_r2

      ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
      !k_b  = 5.9e+9_fp
      ! (Liu and Margerum, Environ. Sci. Tech., 2001)
      k_b1  = 2.3e+10_fp ! (qjc, 12/28/16)
      ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
      k_b2  = 1.6e+10_fp

      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l  = 1.4e-5_fp

      H_X = H_K0_HOBr*dexp(H_CR_HOBr*(1.0e0_fp/T - 1.0e0_fp/H_HOBr_T))
      M_X = MW_HOBr * 1e-3_fp

      ! Mass accommodation coefficient
      ab = 0.6e0_fp

      ! Thermal velocity [cm/s]
      cavg = dsqrt(8*RStarG*T/(pi*M_X)) *1.0e2_fp

      ! l_r is diffusive length scale [cm]; gb is Bulk reaction coefficient [unitless]
      l_r1  = dsqrt( D_l / (k_b1 * C_Y1 * C_Hp ) )
      gb1 = 4.0e0_fp * H_X * con_R * T * l_r1 * k_b1 * C_Y1 * C_Hp / cavg
      gb1 = gb1 * REACTODIFF_CORR( Radius, l_r1)

      l_r2  = dsqrt( D_l / (k_b2 * C_Y2 * C_Hp ) )
      gb2 = 4.0e0_fp * H_X * con_R * T * l_r2 * k_b2 * C_Y2 * C_Hp / cavg
      gb2 = gb2 * REACTODIFF_CORR( Radius, l_r2)

      gb_tot = gb1 + gb2

      ! Reactive uptake coefficient [unitless]
      GAM_HOBr = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb_tot)

      ybr2 = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
      ybr2 = MIN(ybr2, 0.9e0)
      ybr2 = MAX(ybr2, TINY(1.0e0))

      IF ( X==1 ) THEN

         r_gp = (gb1+gb2)/gb_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = TINY(1.0e0)
         ELSE
            r_gp = r_gp * (1.e0 - ybr2)
         ENDIF

      ELSEIF ( X==2 ) THEN

         r_gp = (gb1+gb2)/gb_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.9e0 * r_gp
         ELSE
            r_gp = r_gp * ybr2
         ENDIF

      ENDIF

    END SUBROUTINE GAMMA_HOBr_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO3_Br
!
! !DESCRIPTION: Function GAMMA\_ClNO3\_Br calculates reactive uptake coef.
!               for bromide oxidation by ClNO3
!\\
!\\
! !INTERFACE:
!
      FUNCTION GAMMA_ClNO3_Br( Radius, n_air, T, C_Y ) RESULT( GAM )
!
! !USES:
!
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                         :: GAM
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm)
      REAL(fp), INTENT(IN)             :: Radius, n_air
      REAL(fp), INTENT(IN)             :: T, C_Y
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)       :: ab, gb, gd, M_X
      REAL(fp)       :: cavg, D_l

      M_X = MW_ClNO3 * 1e-3_fp
      ab = 0.11e0_fp

      cavg = dsqrt(8.0e+0_fp*RStarG*T/(Pi*M_X)) *1.0e2_fp ! thermal velocity (cm/s)

      D_l  = 5.0e-6_fp !cm2 s-1.
      gb   = 4.0e0_fp * con_R * T * 1.0e6_fp * dsqrt(C_Y*D_l) / cavg ! H*sqrt(kb)=10^6 (M/s)^½ s-1

      GAM = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb)

      END FUNCTION GAMMA_ClNO3_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COTH
!
! !DESCRIPTION: COTH (Hyperbolic cotangent)
!               coth(x) = cosh(x)/sinh(x) = (1 + exp(-2x))/(1 - exp(-2x))
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION COTH( X)
!
! !INPUT PARAMETERS:
!
      REAL(fp),         INTENT(IN)  :: X           ! The argument
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL(fp)                 :: exp_temp

      exp_temp = dexp(-2.0e0_fp*X)
      COTH = (1.0e0_fp + exp_temp)/(1.0e0_fp - exp_temp)

      RETURN

      END FUNCTION COTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: REACTODIFF_CORR
!
! !DESCRIPTION: REACTODIFF\_CORR
!     Correction =  COTH( x ) - ( 1/x )
!              x = radius / l
!     Correction approaches 1 as x becomes large, corr(x>1000)~1
!     Correction approaches x/3 as x goes towards 0
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION REACTODIFF_CORR( radius, l)
!
! !INPUT PARAMETERS:
!
      REAL(fp),         INTENT(IN)  :: radius, l           ! [cm] and [cm]
!
! !REVISION HISTORY:
!  14 Oct 2013 - J. Schmidt - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL(fp)                 :: x

      x = radius / l

      IF (x<0.0e0_fp) THEN
         PRINT *, 'ERROR x<0, particle radius or C_Y is neg!'
      ELSEIF (x>1.0e3_fp) THEN
         REACTODIFF_CORR = 1.0e0_fp
      ELSEIF (x<1.0e-1_fp) THEN
         REACTODIFF_CORR = x/3.0e0_fp
      ELSE
         REACTODIFF_CORR = COTH(x) - (1.0e0_fp/x)
      ENDIF


      RETURN

      END FUNCTION REACTODIFF_CORR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetBrNO3_JS
!
! !DESCRIPTION: Sets the hydrolysis rate for BrNO3 using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3_JS( denAir, rLiq, rIce, ALiq, AIce, TK ) RESULT( HET_BrNO3 )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rLiq        ! Radius of liquid cloud droplets (cm)
      REAL(fp), INTENT(IN) :: rIce        ! Radius of ice cloud crystals (cm)
      REAL(fp), INTENT(IN) :: ALiq        ! Area of liquid cloud droplets (cm2/cm3)
      REAL(fp), INTENT(IN) :: AIce        ! Area of ice cloud crystals (cm2/cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_BrNO3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  16 Dec 2016 - S. D. Eastham - Updated code based on Johan Schmidt's work
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp), Parameter :: XMolWeight=142.0e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)

      ! Initialize
      HET_BrNO3    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Only apply PSC rate adjustment if at high altitude
      DO_EDUCT     = STRATBOX

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            XSTKCF = 0.02e+0_fp
         ELSEIF ( (N == 11) .OR. ( N == 12) ) THEN
            ! 2 modes of sea-salt
            XSTKCF = 0.02e+0_fp
         ELSEIF ( N == 13 ) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(6)
         ELSEIF ( N == 14 ) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.001e+0_fp
            ELSE
               XSTKCF = 0.3e+0_fp
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         HET_BrNO3 = HET_BrNO3 + ADJUSTEDRATE
      END DO

    ! Hydrolysis on liquid and ice clouds (tropospheric only)
    IF (.not. StratBox) THEN
       HET_BrNO3 = HET_BrNO3 + Cld1K_XNO3(denAir,TK,rLiq,rIce,ALiq,AIce,XMolWeight,2.0E-2_fp)
    ENDIF

    END FUNCTION HETBrNO3_JS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_HCl
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ClNO3(g) + HCl(l,s)
! in polar stratospheric clouds and on tropospheric sulfate.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_HCl( A, B ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  04 May 2016 - M. Sulprizio- Add fixes for setting rate if not a STRATBOX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum          = 0.0_fp
      ADJUSTEDRATE   = 0.0_fp
      XSTKCF         = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Assume zero unless proven otherwise
         XSTKCF = 0e+0_fp

!         IF (N.eq.8) THEN
!            XSTKCF = 0.1e-4_fp ! Sulfate
!         ELSEIF ( STRATBOX ) THEN
!            IF (N.eq.13) THEN
	 ! restore limitation to stratosphere - TMS 17/04/10
         IF  ( STRATBOX ) THEN
            IF (N.eq.8) THEN
               XSTKCF = 0.1e-4_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(4)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.2e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (STRATBOX.and.(N.eq.13)) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE

      END DO

    END FUNCTION HETClNO3_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetBrNO3_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for BrNO3(g) + HCl(l,s)
!  in polar stratospheric clouds and on tropospheric sulfate.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3_HCl( A, B ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  04 May 2016 - M. Sulprizio- Add fixes for setting rate if not a STRATBOX
!  24 Dec 2016 - S. D. Eastham - Extended into the troposphere
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default to zero
         XSTKCF = 0.0e+0_fp

	 ! restore limitation to stratosphere - TMS 17/04/10
	 IF ( STRATBOX ) THEN
	    IF (N.eq.8) THEN
               XSTKCF = 0.9e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(7)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF

      END DO

    END FUNCTION HETBrNO3_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HCl(l,s)
!  in polar stratospheric clouds and on sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  04 May 2016 - M. Sulprizio- Add fixes for setting rate if not a STRATBOX
!  22 Dec 2016 - S. D. Eastham - Now active for non-UCX mechanisms
!  03 Jan 2018 - M. Sulprizio  - Replace UCX CPP switch with Input_Opt%LUCX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF        = 0.0_fp

         ! For UCX-based mechanisms only consider PSC reactions in strat
         ! restore limitation to stratosphere - TMS 17/04/10
         IF ( Input_Opt%LUCX .and. STRATBOX) THEN
	    IF (N.eq.8) THEN
	       XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(8)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.1e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.2e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + AdjustedRate

         ENDIF

      END DO

    END FUNCTION HETHOCl_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_HBr
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HBr(l,s)
!  in polar stratospheric clouds and on sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_HBr( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  29 Mar 2016 - R. Yantosca - Added ProTeX header
!  01 Apr 2016 - R. Yantosca - Define N, XSTKCF, ADJUSTEDRATE locally
!  01 Apr 2016 - R. Yantosca - Replace KII_KI with DO_EDUCT local variable
!  04 May 2016 - M. Sulprizio- Add fixes for setting rate if not a STRATBOX
!  22 Dec 2016 - S. D. Eastham - Now active for non-UCX mechanisms
!  03 Jan 2018 - M. Sulprizio  - Replace UCX CPP switch with Input_Opt%LUCX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF        = 0.0_fp

         ! For UCX-based mechanisms only consider PSC reactions in strat
         ! restore limitation to stratosphere - TMS 17/04/10
         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
	    IF (N.eq.8) THEN
 	       XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(9)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + AdjustedRate

         ENDIF

      END DO

    END FUNCTION HETHOCl_HBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N2O5
!
! !DESCRIPTION: Internal function N2O5 computes the GAMMA sticking factor
!  for N2O5 hydrolysis. (mje, bmy, 8/7/03)
!\\
!\\
! !INTERFACE:
!
    FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )
!
! !INPUT PARAMETERS:
!
      INTEGER,   INTENT(IN) :: AEROTYPE  ! Denoting aerosol type (cf FAST_JX)
      REAL(fp),  INTENT(IN) :: TEMP      ! Temperature [K]
      REAL(fp),  INTENT(IN) :: RH        ! Relative humidity [%]
!
! !RETURN VALUE:
!
      REAL(fp)              :: GAMMA

!
! !REMARKS:
!  Taken from the old SMVGEAR function calcrate.F.
!
! !REVISION HISTORY:
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!  15 Jun 2017 - M. Sulprizio- Move conversion of RH from fraction to % to
!                              SET_HET above
!  23 Aug 2018 - C. D. Holmes - Updated Gamma values
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL(fp) :: RH_P, FACT, TTEMP

      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! RH percent max = 100%
      RH_P  = MIN( RH, 100e0_fp )

      ! Default value
      GAMMA = 0.01e+0_fp

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )

            GAMMA = 0.02e+0_fp

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )
            ! NOTE: The following calculation for sulfate aerosol is
            ! done only for the stratosphere following the UCX implementation
            ! For troposphere, we use the McDuffie parameterization (HetN2O5())

            !===========================================================
            ! RH dependence from Kane et al., Heterogenous uptake of
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470
            !===========================================================
            ! No RH dependence above 50.0% (lzh, 10/26/2011)
            ! According to Bertram and Thornton, ACP, 9, 8351-8363, 2009
            RH_P  = MIN( RH_P, 50e+0_fp )

            GAMMA = 2.79e-4_fp + RH_P*(  1.30e-4_fp +    &
                              RH_P*( -3.43e-6_fp +       &
                              RH_P*(  7.52e-8_fp ) ) )

            !===========================================================
            ! Temperature dependence factor (Cox et al, Cambridge UK)
            ! is of the form:
            !
            !          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
            ! FACT = -------------------------------------------------
            !                     10^( LOG10( G294 ) )
            !
            ! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
            !
            ! For computational speed, replace LOG10( 1e-2 ) with -2
            ! and replace 10^( LOG10( G294 ) ) with G294
            !===========================================================
            TTEMP = MAX( TEMP, 282e0_fp )
            FACT  = 10.e0_fp**( -2e+0_fp - 4e-2_fp       &
                  *( TTEMP - 294.e+0_fp ) ) / 1e-2_fp

            ! Apply temperature dependence
            GAMMA = GAMMA * FACT

         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )

             ! From IUPAC
             GAMMA = 0.005e+0_fp

         !----------------
         ! Organic Carbon
         !----------------
         CASE ( 10 )

            ! Escoria et al (2010)
            IF ( RH_P < 30e+0_fp ) THEN
               GAMMA = 0.6e-4_fp
            ELSE
               GAMMA = 1.5e-4_fp
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )

            ! IUPAC and Thornton and Abbatt (2005)
            if (RH_P < 40) then
               gamma = 0.005e0_fp
            elseif (RH_P > 70) then
               gamma = 0.02e0_fp
            else
               gamma = 0.005e0_fp + (0.02e0_fp - 0.005e0_fp) * &
                    ( RH_P - 40e0_fp ) / ( 70e0_fp - 40e0_fp )
            endif

         !----------------
         ! Strat. aerosols
         !----------------
         CASE ( 13, 14 )

            ! These are handled outside this routine - something
            ! is wrong if AEROTYPE=13 or 14 reaches this point
            WRITE (6,*) 'Stratospheric aerosols should not '
            WRITE (6,*) 'be passed to general N2O5 het. '
            WRITE (6,*) 'chem. subroutine'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for N2O5 hydrolysis'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT

    END FUNCTION N2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HO2
!
! !DESCRIPTION: Function HO2 computes the GAMMA reaction probability
!  for HO2 loss in aerosols based on the recommendation of
!  Thornton, Jaegle, and McNeill, "Assessing Known Pathways For HO2 Loss in
!  Aqueous Atmospheric Aerosols: Regional and Global Impacts on Tropospheric
!  Oxidants" J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008
!\\
!\\
! !INTERFACE:
!
      FUNCTION HO2( RADIUS,          TEMP,     DENAIR,       &
                    SQM,             HO2DENS,  AEROTYPE,     &
                    CONTINENTAL_PBL, Input_Opt            )  &
                    RESULT( GAMMA )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
      Use PhysConstants, ONLY : AVO, RGASLATM
!
! !INPUT PARAMETERS:
!
      ! Arguments
      REAL(fp),       INTENT(IN) :: RADIUS          ! Aerosol radius [cm]
      REAL(fp),       INTENT(IN) :: TEMP            ! Temperature [K]
      REAL(fp),       INTENT(IN) :: DENAIR          ! Air density [molec/cm3]
      REAL(fp),       INTENT(IN) :: HO2DENS         ! HO2 density [molec/cm3]
      REAL(fp),       INTENT(IN) :: SQM             ! Square root of MW [g/mol]
      INTEGER,        INTENT(IN) :: AEROTYPE        ! Aerosol type (cf FAST-JX)
      INTEGER,        INTENT(IN) :: CONTINENTAL_PBL ! Flag set to 1 if the box
                                                    !  box is located in the
                                                    !  continenal boundary
                                                    !  layer, otherwise 0.
                                                    !  Also check for ICE/SNOW
                                                    !  (to disable this at
                                                    !  high latitudes).
      TYPE(OptInput), INTENT(IN) :: Input_Opt       ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: GAMMA           ! Reaction probability

! !REMARKS:
!  Taken from the old SMVGEAR routine calcrate.F.
!  Gamma(HO2) is a function of aerosol type, radius, temperature.

!  References:
!  ---------------------------------------------------------------
!  (1) Jacob, D.J., Heterogeneous chemistry and tropospheric ozone,
!       Atmos. Environ., 34, 2131-2159, 2000. [full text (pdf)]
!  (2) J. Mao, Fan, S., Jacob, D. J., and Travis, K. R.: Radical
!       loss in the atmosphere from Cu-Fe redox coupling in aerosols,
!       Atmos. Chem. Phys., 13, 509-519, doi:10.5194/acp-13-509-2013,
!       2013.
!
! !REVISION HISTORY:
!  17 May 2013 - M. Payer    - Add improved HO2 uptake (J. Mao)
!  22 May 2013 - M. Payer    - Added option to read GAMMA_HO2 from
!                              input.geos. Recommended value is 0.2
!                              based on Jacob et al (2000) and Mao
!                              et al. (2013).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  29 Mar 2016 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: ALPHA
      REAL(fp)             :: delG, Keq, w, H_eff
      REAL(fp)             :: A1, B1, k1, k2, A, B, C
      REAL(fp)             :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL(fp)             :: TEST
!
! !DEFINED PARAMETERS:
!
      ! Ideal gas constant [atm cm3/mol/K]
      REAL(fp),  PARAMETER :: Raq = RGASLATM * 1e+3_fp

      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0e+0_fp

      ! Error check
      IF (RADIUS.le.1e-30_fp) THEN
         RETURN
      ENDIF

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )

            ! Assume default gamma=0.1 on dust aerosols
            ! This is tentative as no lab measurements presently exist
            ! for gamma(HO2) on dust aerosols. We assume the rate to
            ! be fast on dust aerosols as transition metal ion induced
            ! chemistry is likely to occur in a thin aqueous surface layer.
            GAMMA = 0.1e+0_fp

         !----------------
         ! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
         ! Sea-salt accum & coarse (11,12) calculate the
         ! reaction probability due to self reaction
         ! by using the algebraic expression in Thornton et al.  (2008)
         ! (equation 7) which is a function of temperature, aerosol radius,
         ! air density and HO2 concentration.
         !
         ! Transition metal ions (such as copper and iron) in sea-salt and
         ! carbonaceous aerosols are complexed to ligands and/or exist at
         ! a concentration too low to catalyze HO2 loss efficiently, so we
         ! apply the HO2 self reaction expression directly for these aerosols.
         !
         ! In the case of sulfate aerosol, the aerosols likely
         ! contain copper in the continental boundary layer and
         ! HO2 uptake proceeds rapidly. To account for the metal catalyzed
         ! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
         ! 0.04-0.1 by Thornton et al, based on observed copper concentrations
         ! in the US boundary layer). Outside the continental boundary layer, we
         ! use the HO2-only algebraic expression.
         !
         ! SDE 04/18/13: Added stratospheric sulfur aerosols
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12, 13 )

            ! Mean molecular speed [cm/s]
            w = 14550.5e+0_fp * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45E+17_fp/DENAIR * SQRT(TEMP) *      &
                    SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus)
            ! with Keq = 2.1d-5 [M], Equilibrium constant for
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.e+0_fp^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1e-5_fp / (10.e+0_fp**(-5) ) )
            B1 = 1.987e-3_fp * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298e+0_fp)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9e+0_fp - (TEMP-298.e+0_fp) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58e+10_fp * exp( -3. / B1 )
            k2  =   2.4e+9_fp   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.e+0_fp) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4e+0_fp/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed
            ! to be 1 (Thornton et al.)
            kmt = 1.e+0_fp/( RADIUS/DFKG + 4e+0_fp/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2e+0_fp * kaq
            B = -3e+0_fp * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3e+0_fp * kmt * HO2DENS * 1000e+0_fp / RADIUS / AVO

            ! Error check that B^2-(4e+0_fp*A*C) is not negative
            TEST= B**2-(4e+0_fp*A*C)
            IF ( TEST < 0e+0_fp ) THEN
                GAMMA = 0e+0_fp
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4e+0_fp*A*C)) )/(2e+0_fp*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*AVO*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0e0_fp ) THEN
                   GAMMA = 0e+0_fp
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) -              &
                                  ( RADIUS/DFKG ) ) * w / 4.e+0_fp )
                ENDIF
            ENDIF
            ! For sulfate aerosols, check whether we are in
            ! the continental boundary layer, in which case
            ! copper catalyzed HO2 uptake likely dominates and
            ! speeds up the reaction: we assume gamma=0.07,
            ! which is in the middle of the 0.04-0.1 range recommended
            ! by Thornton et al. (2008)
            !
            IF ( AEROTYPE == 8 .and. CONTINENTAL_PBL == 1) THEN
                GAMMA = 0.07
            ENDIF

         !----------------
         ! NAT/ice (SDE 04/18/13)
         !----------------
         CASE ( 14 )

            GAMMA = 0.e+0_fp

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for HO2 uptake'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT

      ! If negative value is calculated, set it to zero
      IF ( GAMMA  <= 0e+0_fp ) GAMMA = 0e+0_fp

      ! This is for the improved HO2 uptake (J. Mao)
      GAMMA = Input_Opt%GAMMA_HO2

    END FUNCTION HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EPOXUPTK
!
! !DESCRIPTION: Function EPOXUPTK computes the GAMMA sticking factor for
! EPOXUPTK hydrolysis to form 2-methyltetrols (AITET). (eam, 2014).
!\\
!\\
! !INTERFACE:
!
    FUNCTION EPOXUPTK( AERAREA, AERRAD,  TEMP,  SQMW,              &
                       HENRY,   KHPLUS,  HPLUS, KNUC,  SULF, NITR, &
                       KGACID,  BISULF,  KHYDRO ) &
             RESULT( GAMMA )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
      USE PhysConstants, ONLY : RGASLATM
      USE ERROR_MOD,     ONLY : IT_IS_NAN
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: AERRAD   ! Aerosol radius [cm]
      REAL(fp), INTENT(IN) :: AERAREA  ! Aerosol surf. area [cm2/cm3]
      REAL(fp), INTENT(IN) :: TEMP     ! Temperature [K]
      REAL(fp), INTENT(IN) :: SQMW     ! Square root of the molecular weight
      REAL(fp), INTENT(IN) :: HENRY    ! Henry's Law constant [M/atm]
      REAL(fp), INTENT(IN) :: KHPLUS   ! 1st order rxn rate (acid-catalyzed ring
                                       ! opening)
      REAL(fp), INTENT(IN) :: HPLUS    ! Proton activity [unitless] and [H+] [M]
      REAL(fp), INTENT(IN) :: KNUC     ! 1st order rxn rate due to specific
                                       ! nucleophiles (SO4, NO3)
      REAL(fp), INTENT(IN) :: SULF     ! Sulfate concentration [M]
      REAL(fp), INTENT(IN) :: NITR     ! Nitrate concentration [M]
      REAL(fp), INTENT(IN) :: KGACID   ! 1st order rxn rate due to general acids
                                       ! (bisulfate in this case)
      REAL(fp), INTENT(IN) :: BISULF   ! Bisulfate concentration [M]
      REAL(fp), INTENT(IN) :: KHYDRO   ! Hydrolysis rate of alkylnitrates [1/s]
!
! !RETURN VALUE:
!
      REAL(fp)             :: GAMMA    ! Reaction probability

! !REMARKS:
! Calculation is only done for inorganic aqueous phase aerosols
!                                                                             .
! This calculation uses the parameterization of Gaston et al., EST, 2014.
!                                                                             .
! Redistribution of products (e.g. AITET) to yield organosulfates and
! organonitrates is done in SOA_CHEMISTRY in carbon_mod.F.
! This is only done for IEPOX and IMAE if it's an SOA simulation
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL(fp)             :: AERVOL   ! Aerosol volume [cm3/cm3]
      REAL(fp)             :: KPART    ! Particle-phase reaction rate [1/s]
      REAL(fp)             :: XMMS     ! Mean molecular speed [cm/s]
      REAL(fp)             :: VAL1, VAL2, VAL3  ! Terms for calculating GAMMA
      REAL(fp)             :: VALTMP
!
! !DEFINED PARAMETERS:
!
      ! Gas-phase diffusion constant [cm2/s]:
      REAL(fp), PARAMETER  :: DIFF_N2O5_STD = 1.0e-1_fp

      ! Mass accommodation coefficient [unitless]:
      REAL(fp), PARAMETER  :: MACOEFF = 1.0e-1_fp

      !=================================================================
      ! EPOXUPTK begins here!
      !=================================================================

      ! Initialize
      GAMMA  = 0.0_fp
      AERVOL = 0.0_fp
      KPART  = 0.0_fp
      XMMS   = 0.0_fp
      VAL1   = 0.0_fp
      VAL2   = 0.0_fp
      VAL3   = 0.0_fp
      VALTMP = 0.0_fp

      ! Calculate aerosol volume (use formula in aerosol_mod.F):
      AERVOL = (AERAREA * AERRAD)/3.0e+0_fp

      ! Calculate mean molecular speed [cm/s]:
      XMMS = SQRT( (2.117e+8_fp * TEMP) / (SQMW * SQMW) )

      ! Calculate first-order particle-phase reaction rate:
      ! (assume [H+] = proton activity)
      ! KHYDRO is only important for alkylnitrates (not currently used).
      KPART = ( KHPLUS*HPLUS )               + &
              ( KNUC*HPLUS*( NITR + SULF ) ) + &
              ( KGACID*BISULF )              + &
              ( KHYDRO )

      ! Calculate the first uptake parameterization term:
      VAL1 = ( AERRAD * XMMS )/( 4.e+0_fp * DIFF_N2O5_STD )

      ! Calculate the second uptake parameterization term:
      VAL2 = ( 1.e+0_fp/MACOEFF )

      ! Calculate the third uptake parameterization term:
      IF ( AERAREA > 0.0_fp .and. XMMS > 0.0_fp ) THEN
         VALTMP = ( 4.e+0_fp * AERVOL * RGASLATM * TEMP * HENRY * KPART ) / &
                  ( AERAREA * XMMS )
      ENDIF
      IF ( VALTMP .GT. 0 ) THEN
         VAL3 = 1.e+0_fp / VALTMP
      ELSE
         VAL3 = 0.0e+0_fp
      ENDIF

      ! Account for small reaction rates:
      IF ( KPART .LT. 1.e-8_fp ) THEN

         GAMMA = TINY(1e+0_fp)

      ELSE

         ! Calculate the uptake coefficient:
         GAMMA = 1.e+0_fp/( VAL1 + VAL2 + VAL3 )

      ENDIF

      ! Fail safes for negative, very very small, and NAN GAMMA values:
      IF ( GAMMA  .lt. 0.0e+0_fp )    GAMMA = TINY(1e+0_fp)
      IF ( IT_IS_NAN( GAMMA ) )       GAMMA = TINY(1e+0_fp)
      IF ( GAMMA .lt. TINY(1e+0_fp) ) GAMMA = TINY(1e+0_fp)

      END FUNCTION EPOXUPTK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cld_Params
!
! !DESCRIPTION: Subroutine CLD\_PARAMS returns ice and liquid cloud
!  parameters based on State\_Met off of cloud particles.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE CLD_PARAMS( I,      J,      L,    DENAIR,            &
                           VAir,   T,      QL,   QI,     State_Met, &
                           rLiq,   ALiq,   VLiq, &
                           rIce,   AIce,   VIce, CLDFr )

!
! !USES:
!
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I         ! Longitude index
      INTEGER,        INTENT(IN)  :: J         ! Latitude  index
      INTEGER,        INTENT(IN)  :: L         ! Altitude  index
      REAL(fp),       INTENT(IN)  :: DENAIR    ! Density of air [#/cm3]
      REAL(fp),       INTENT(IN)  :: VAir      ! Volume of air [cm3]
      REAL(fp),       INTENT(IN)  :: T         ! Temperature [K]
      REAL(fp),       INTENT(IN)  :: QL, QI    ! Cloud water mixing ratio [kg/kg]
      TYPE(MetState), INTENT(IN)  :: State_Met ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: rLiq     ! Radius of liquid cloud droplets (cm)
      REAL(fp),       INTENT(OUT) :: rIce     ! Radius of ice cloud crystals (cm)
      REAL(fp),       INTENT(OUT) :: ALiq     ! Sfc area of liq. cloud (cm2/cm3)
      REAL(fp),       INTENT(OUT) :: AIce     ! Sfc area of ice cloud (cm2/cm3)
      REAL(fp),       INTENT(OUT) :: VLiq     ! Volume of liq. cloud (cm3/cm3)
      REAL(fp),       INTENT(OUT) :: VIce     ! Volume of ice cloud (cm3/cm3)
      REAL(fp),       INTENT(OUT) :: CLDFr     ! cloud fraction
!
! !REMARKS:
!  References:
!   Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships between ice water content and volume extinction coefficient from in situ observations for temperatures from 0° to –86°C: implications for spaceborne lidar retrievals. Journal of Applied Meteorology and Climatology, 53(2), 479–505. https://doi.org/10.1175/JAMC-D-13-087.1
!
!   Schmitt, C. G., & Heymsfield, A. J. (2005). Total Surface Area Estimates for Individual Ice Particles and Particle Populations. Journal of Applied Meteorology, 44(4), 467–474. https://doi.org/10.1175/JAM2209.1
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Adapted from CLD1K_BrNO3
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  15 Oct 2018 - C.D. Holmes - Corrections for ice radius, volume, surface area
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Cloud droplet radius in continental warm clouds [cm]
      REAL(fp), PARAMETER :: xCLDR_CONT =  6.e-4_fp

      ! Cloud droplet radius in marine warm clouds [cm]
      REAL(fp), PARAMETER :: xCLDR_MARI = 10.e-4_fp

      ! Ice cloud droplet radius [cm]
      REAL(fp), PARAMETER :: xCLDrIce = 75.e-4_fp

      ! Density of H2O liquid [kg/cm3]
      REAL(fp), PARAMETER :: dens_h2o = 1.0e-3_fp

      ! Density of H2O ice [kg/cm3]
      REAL(fp), PARAMETER :: dens_ice = 0.91e-3_fp

! !LOCAL VARIABLES:
!
      LOGICAL             :: IS_LAND, IS_ICE, Is_Warm
      REAL(fp)            :: alpha,   beta

      ! Pointers
      REAL(fp), POINTER   :: AD(:,:,:)
      REAL(fp), POINTER   :: CLDF(:,:,:)
      REAL(fp), POINTER   :: FRLAND(:,:)
      REAL(fp), POINTER   :: FROCEAN(:,:)

      !=================================================================
      ! CLD_PARAMS begins here!
      !=================================================================

      ! Initialize pointers
      AD      => State_Met%AD
      CLDF    => State_Met%CLDF
      FRLAND  => State_Met%FRLAND
      FROCEAN => State_Met%FROCEAN

      CLDFr = CLDF(I,J,L)
      IF ( CLDFR.le.0.0e+0_fp ) CLDFr = 1.0e-32_fp

      ! Quick test - is there any cloud?
      IF ( ( (QL+QI) <= 0.0e+0_fp) .or. (CLDF(I,J,L) <= 0.0e+0_fp) ) THEN
         rLiq = xCldR_Cont
         rIce = xCLDrIce
         ALiq = 0.0e+0_fp
         VLiq = 0.0e+0_fp
         AIce = 0.0e+0_fp
         VIce = 0.0e+0_fp
         Return
      ENDIF

      ! Volume of cloud condensate, water or ice [cm3]
      ! QL is [g/g]
      VLiq = QL * AD(I,J,L) / dens_h2o
      VIce = QI * AD(I,J,L) / dens_ice

      ! ----------------------------------------------
      ! In GC 12.0 and earlier, the liquid water volume was
      ! set to zero at temperatures colder than 258K and
      ! over land ice (Antarctica & Greenland). That
      ! was likely legacy code from GEOS-4, which provided
      ! no information on cloud phase. As of GC 12.0,
      ! all met data sources provide explicit liquid and
      ! ice condensate amounts, so we use those as provided.
      ! (C.D. Holmes)
      ! ----------------------------------------------

      !-----------------------------------------------
      ! Liquid water clouds
      !
      ! Droplets are spheres, so
      ! Surface area = 3 * Volume / Radius
      !
      ! Surface area density = Surface area / Grid volume
      !-----------------------------------------------

      IF ( FRLAND(I,J) > FROCEAN(I,J) ) THEN
         ! Continental cloud droplet radius [cm]
         rLiq = XCLDR_CONT
      ELSE
         ! Marine cloud droplet radius [cm]
         rLiq = XCLDR_MARI
      ENDIF

      ! Surface area density, cm2/cm3
      ALiq = 3.e+0_fp * (VLiq/VAir) / rLiq

      !-----------------------------------------------
      ! Ice water clouds
      !
      ! Surface area calculation requires information about
      ! ice crystal size and shape, which is a function of
      ! temperature. Use Heymsfield (2014) empirical relationships
      ! between temperature, effective radius, surface area
      ! and ice water content.
      !
      ! Schmitt and Heymsfield (2005) found that ice surface area
      ! is about 9 times its cross-sectional area.
      ! For any shape,
      ! Cross section area = pi * (Effective Radius)^2, so
      ! Cross section area = 3 * Volume / ( 4 * Effective Radius ).
      !
      ! Thus, for ice
      ! Surface area = 9 * Cross section area
      !              = 2.25 * 3 * Volume / Effective Radius
      ! (C.D. Holmes)
      !-----------------------------------------------

      ! Heymsfield (2014) ice size parameters
      if (T < 273e+0_fp - 71e+0_fp ) then
          alpha = 83.3e+0_fp
          beta  = 0.0184e+0_fp
      elseif ( T < 273e+0_fp - 56e+0_fp ) then
          alpha = 9.1744e+4_fp
          beta = 0.117e+0_fp
      else
          alpha = 308.4e+0_fp
          beta  = 0.0152e+0_fp
      endif

      ! Effective radius, cm
      rIce = 0.5e+0_fp * alpha * exp( beta * (T-273.15e+0_fp) ) / 1e+4_fp

      ! Ice surface area density, cm2/cm3
      AIce = 3.e+0_fp * (VIce/VAir) / rIce * 2.25e+0_fp

      ! Free Pointers
      NULLIFY( AD      )
      NULLIFY( CLDF    )
      NULLIFY( FRLAND  )
      NULLIFY( FROCEAN )

    END SUBROUTINE Cld_Params
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halide_CldConc
!
! !DESCRIPTION: Subroutine GET\_HALIDE\_CLDCONC returns the in-cloud
!  concentration of bromide and chloride (Br- and Cl-).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HALIDE_CLDCONC( HBr, HCl, VLiq, VIce, VAir, CLDFr, &
                 TK, SA_SULF, R_SULF, br_conc, cl_conc )

!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: HCl, HBr  ! Number density [#/cm3]
      REAL(fp),  INTENT(IN) :: VAir    ! Volume of air [cm3]
      REAL(fp),  INTENT(IN) :: SA_SULF, R_SULF! Sulfate aerosol surface area (cm2/cm3) and radius (cm)
      REAL(fp),  INTENT(IN) :: VLiq, VIce ! Volume of the cloud (liq and ice) [cm3]
      REAL(fp),  INTENT(IN) :: CLDFr ! cloud fraction
      REAL(fp),  INTENT(IN) :: TK      ! Air temperature [K]

!
! !RETURN VALUE:
!
      REAL(fp), INTENT(OUT) :: cl_conc, br_conc ! Liq. phase molar concentration [mol/kg-water]
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Initial version
!  27 Feb 2018 - M. Sulprizio  - Obtain Henry's law parameters from species
!                                database in SET_HET instead of hardcoding here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: n_br, n_cl ! dissolved bromide and chloride [#/cm3(air)]
      REAL(fp)            :: V_tot, dr_ratio, t2l !
      REAL(fp)            :: L2G, F_L

      !=================================================================
      ! GET_HALIDE_CLDCONC begins here!
      !=================================================================

      !---------------------------------------------------------------
      ! jas, 07/30/2014 (SETUP d/r ratio for ice cloud droplets)
      ! V_liq = 4pi/3 ( r^3 - (r - r*(d/r))^3 = (r^3 - r^3*(1 - d/r)^3) = r^3 (1 - (1 - d/r)^3
      ! V_tot / V_liq = 1 / (1 - (1 - d/r)^3))
      DR_RATIO = 2e-2_fp
      T2L = 1.0e0_fp / ( 1.0e0_fp - (1.0e0_fp - DR_RATIO)**3.0e0_fp )
      !---------------------------------------------------------------

      V_tot = (VLiq/CLDFr/VAir) + ((VIce/CLDFr/VAir) / T2L) + &
               SA_SULF * R_SULF / 3.0e0_fp  ! (cm3(liq)/cm3(air)

      IF (V_tot.lt.1.0e-20) THEN
         br_conc = 1.0e-20_fp
         cl_conc = 1.0e-20_fp
         Return
      ENDIF

      ! Bromide (Assuming ph=4.5)
      CALL COMPUTE_L2G_LOCAL( H_K0_HBr, H_CR_HBr, 0.0e+0_fp, TK, V_tot, L2G)
      F_L = L2G/(1.0e0_fp + L2G)
      br_conc = F_L * HBr / (V_tot * AVO * 1.0e-3_fp) ! [Br-] in (mol/L)

      br_conc = min(br_conc,5.0e0_fp)
      br_conc = max(br_conc,1.0e-20_fp)

      ! Chloride (Assuming ph=4.5)
      CALL COMPUTE_L2G_LOCAL( H_K0_HCl, H_CR_HCl, 0.0e+0_fp, TK, V_tot, L2G)
      F_L = L2G/(1.0e0_fp + L2G)
      cl_conc = F_L * HCl / (V_tot * AVO * 1.0e-3_fp) ! [Cl-] in (mol/L)
      cl_conc = min(cl_conc,5.0e0_fp)
      cl_conc = max(cl_conc,1.0e-20_fp)

      END SUBROUTINE GET_HALIDE_CLDCONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halide_SSAConc
!
! !DESCRIPTION: Function GET\_HALIDE\_SSACONC calculates concentration of a
!               halide in sea salt aerosol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HALIDE_SSACONC( n_x, surf_area, r_w, conc_x )
!
! !OUTPUT PARAMETER:
      ! concentration of X- in SALX (mol/L)
      REAL(fp)                         :: conc_x
! !INPUT PARAMETERS:
      ! n_x = X-(ssa) number density (#/cm3), surf_area = AERO surface area
      ! conc (cm2/cm3), r_w = AERO wet radius (cm)
      REAL(fp), INTENT(IN)             :: n_x, surf_area, r_w

!
! !REVISION HISTORY:
!  25 Jul 2014 - J. Schmidt - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      !REAL(fp),  PARAMETER :: con_NA      = 6.0221413e23_fp ! #/mol
      REAL(fp)             :: V_tot

      V_tot = surf_area * r_w * 0.3333333e0_fp * 1e-3_fp ! L(liq)/cm3(air)
      IF (V_tot .le. 1.0e-20) THEN
         conc_x = 1.0e-20_fp
         Return
      ELSE
         conc_x =  (n_x / AVO) / V_tot ! mol/L
         conc_x = MIN(conc_x,5.0e0_fp)
         conc_x = MAX(conc_x,1.0e-20_fp)
      ENDIF

      END SUBROUTINE GET_HALIDE_SSACONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_L2G_Local
!
! !DESCRIPTION: Subroutine COMPUTE\_L2G\_LOCAL is a local copy of the
!  liquid-gas partitioning routine in GEOS-Chem's wetscav\_mod.F file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPUTE_L2G_LOCAL( K0, CR, pKa, TK, H2OLIQ, L2G )
!
! !USES:
!
      USE Henry_Mod, ONLY : Calc_KH
      USE Henry_Mod, ONLY : Calc_Heff
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: K0     ! Henry's solubility constant [M/atm]
      REAL(fp), INTENT(IN)  :: CR     ! Henry's volatility constant [K]
      REAL(fp), INTENT(IN)  :: pKa    ! Henry's pH correction factor [1]
      REAL(fp), INTENT(IN)  :: TK     ! Temperature [K]
      REAL(fp), INTENT(IN)  :: H2OLIQ ! Liquid water content [cm3 H2O/cm3 air]
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: L2G    ! Cliq/Cgas ratio [1]
!
! !REMARKS:
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate
!  values of Kstar298 and H298_R must be supplied for each species.
!  (cf Jacob et al 2000, p. 3)
!
! !REVISION HISTORY:
!  23 Feb 2000 - R. Yantosca - Initial version
!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10-Jan-2011 - H. Amos - Corrected the units on KStar298 from moles/atm
!                          to M/atm
!  15-May-2013 - F. Paulot - Fix R constant
!  08 Dec 2015 - R. Yantosca - Now use functions from henry_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: RC
      REAL(fp) :: HEFF, KH, pH, TK_8

      !=================================================================
      ! COMPUTE_L2G_LOCAL begins here!
      !=================================================================

      ! Cast temperature to REAL*8
      TK_8 = TK

      ! For wetdep, we assume a pH of 4.5 for rainwater
      pH = 4.5_fp

      ! Calculate the Henry's law constant
      CALL CALC_KH( K0, CR, TK_8, KH, RC )

      ! Calculate effective Henry's law constant, corrected for pH
      ! (for those species that have a defined pKa value)
      CALL CALC_HEFF( pKa, pH, KH, HEFF, RC )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = HEFF * H2OLIQ

      END SUBROUTINE COMPUTE_L2G_LOCAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cld1k_XNO3
!
! !DESCRIPTION: Function CLD1K\_XNO3 calculates the rate constant for
!  heterogeneous cycling of XNO3 off of cloud particles.
!\\
!\\
! !INTERFACE:
!
    FUNCTION CLD1K_XNO3( denAir, TK, rLiq, rIce, ALiq, AIce, &
                         MX_gmol, AlphaX )    RESULT( cld1k )

!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL(fp),       Intent(IN) :: DENAIR   ! Density of air [#/cm3]
      REAL(fp),       Intent(In) :: TK       ! Air temperature [K]
      REAL(fp),       Intent(In) :: rLiq     ! Radius of liquid cloud drops [cm]
      REAL(fp),       Intent(In) :: rIce     ! Radius of ice cloud crystals [cm]
      REAL(fp),       Intent(In) :: ALiq     ! Surface area (liquid) [cm2/cm3]
      REAL(fp),       Intent(In) :: AIce     ! Surface area (ice) ) [cm2/cm3]
      REAL(fp),       Intent(IN) :: MX_gmol  ! Molecular mass of XNO3 [g/mol]
      REAL(fp),       Intent(IN) :: AlphaX   ! XNO3 accomodation coef [unitless]
!
! !RETURN VALUE:
!
      REAL(fp)                   :: cld1k    ! Rate constant for heterogeneous
                                             ! cycling of BrNO3 off of cloud
                                             ! particles
!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Feb 2011 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  06 Nov 2014 - R. Yantosca - Now use State_Met%CLDF(I,J,L)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: SQM        ! Square root of molec. weight [g/mol]
      REAL(fp)             :: STK        ! Square root of temperature   [K]
      REAL(fp)             :: DFKG       ! Gas diffusion coefficient    [cm2/s]

      !=================================================================
      ! CLD1K_XNO3 begins here!
      !=================================================================

      ! Quick test - is there any cloud?
      IF ((ALiq.le.0.0e+0_fp).and.(AIce.le.0.0e+0_fp)) THEN
         cld1k = 0.0e+0_fp
         Return
      ENDIF

      ! ------------------------------------------------------------
      !   Calculate the 1st order rate constant for XNO3 hydrolysis.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ------------------------------------------------------------
      SQM = sqrt(MX_gmol)    ! square root of molar mass [g/mole]
      STK = sqrt(TK) ! square root of temperature [K]

      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

      ! Compute ARSL1K according to the formula listed above
      ! Sum contribution from ice and liquid clouds
      cld1k = ALiq / ( rLiq/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) )
      cld1k = AIce / ( rIce/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) ) + cld1k

    END FUNCTION CLD1K_XNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FcrO2HO2
!
! !DESCRIPTION: !fgap, based on saunder 2003 k14
!\\
!\\
! !INTERFACE:
!
    FUNCTION FCRO2HO2( XCARBN ) RESULT( FC_RO2HO2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: XCARBN
!
! !RETURN VALUE:
!
      REAL(fp)             :: FC_RO2HO2
!
! !REVISION HISTORY:
!  24 Jul 2014 - R. Yantosca - Now inlined to calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC

      FC_RO2HO2 = 1E+0_fp - EXP( -0.245E+0_fp * XCARBN )

    END FUNCTION FCRO2HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FyHORO
!
! !DESCRIPTION: \subsection*{Overview}
!  Function FYHORO returns returns the branching ratio between
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O
!
!\subsection*{References}
!  \begin{enumerate}
!  \item Orlando et al., 1998: \emph{Laboratory and theoretical study of the
!         oxyradicals in the OH- and Cl-initiated oxidation of ethene},
!        \underline{J. Phys. Chem. A}, \textbf{102}, 8116-8123.
!  \item Orlando et al., 2003: \emph{The atmospheric chemistry of alkoxy
!         radicals}, \underline{Chem. Rev.}, \textbf{103}, 4657-4689.
!  \end{enumerate}
!
! !INTERFACE:
!
    FUNCTION FYHORO( ZDNUM, TT ) RESULT( FY_HORO )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: ZDNUM   ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT      ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FY_HORO
!
! !REVISION HISTORY:
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM  = ZDNUM * 0.21E+0_fp
      K1      = 6.0E-14_fp * EXP(-550.E+0_fp/TT) * O2DNUM
      K2      = 9.5E+13_fp * EXP(-5988.E+0_fp/TT)

      FY_HORO = K1 / (K1 + K2)

    END FUNCTION FYHORO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FyrNO3
!
! !DESCRIPTION: Function FYRNO3 returns organic nitrate yields
!  YN = RKA/(RKA+RKB) from RO2+NO reactions as a function of the number
!  N of carbon atoms.
!\\
!\\
! !INTERFACE:
!
    FUNCTION FYRNO3( XCARBN, ZDNUM, TT ) RESULT( FYR_NO3 )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: XCARBN   ! Number of C atoms in RO2
      REAL(fp), INTENT(IN) :: ZDNUM    ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT       ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FYR_NO3
!
! !REVISION HISTORY:
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  (2 ) Updated following Atkinson 1990.
!  (3 ) Change yield from Isoprene Nitrate (ISN2) from 0.44% to 12%,
!        according to Sprengnether et al., 2002. (amf, bmy, 1/7/02)
!  (4 ) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Updated comment description of XCARBN (bmy, 6/26/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  25 Jul 2014 - R. Yantosca - Now inlined into calcrate.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
      REAL(fp) :: XF,   ALPHA, Y300, BETA, XMINF, XM0

      ! Initialize static variables
      DATA   Y300,ALPHA,BETA,XM0,XMINF,XF/.826,1.94E-22,.97,0.,8.1,.411/

      !=================================================================
      ! FYRNO3 begins here!
      !=================================================================
      XXYN    = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300./TT)**XM0)
      YYYN    = Y300*((300./TT)**XMINF)
      AAA     = LOG10(XXYN/YYYN)
      ZZYN    = 1./(1.+ AAA*AAA )
      RARB    = (XXYN/(1.+ (XXYN/YYYN)))*(XF**ZZYN)
      FYR_NO3 = RARB/(1. + RARB)

    END FUNCTION FYRNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Arsl1k
!
! !DESCRIPTION: Function ARSL1K calculates the 1st-order loss rate of species
!  on wet aerosol surface.
!\\
!\\
! !INTERFACE:
!
    FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM ) &
         RESULT( ARS_L1K )
!
! !INPUT PARAMETERS:
!
      ! Surface  area of wet aerosols/volume of air [cm2/cm3]
      REAL(fp), INTENT(IN) :: AREA

      ! Radius of wet aerosol [cm], order of 0.01-10 um;
      ! Note that radius here is Rd, not Ro
      REAL(fp), INTENT(IN) :: RADIUS

      ! Density of air [#/cm3]
      REAL(fp), INTENT(IN) :: DENAIR

      ! Sticking coefficient [unitless], order of 0.1
      REAL(fp), INTENT(IN) :: STKCF

      ! Square root of temperature [K]
      REAL(fp), INTENT(IN) :: STK

      ! Square root of molecular weight [g/mole]
      REAL(fp), INTENT(IN) :: SQM
!
! !RETURN VALUE:
!
      REAL(fp)             :: ARS_L1K
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)

! !REVISION HISTORY:
!  01 Jul 1994 - lwh, jyl, gmg, djj - Initial version
!  04 Apr 2003 - R. Yantosca - Updated comments, cosmetic changes
!  07 Apr 2004 - R. Yantosca - Now return w/ default value if RADIUS is zero
!                              (i.e. is smaller than a very small number)
!  03 Dec 2009 - R. Yantosca - Prevent div-by-zero errors by returning the
!                              default value if any of the args are zero
!  03 Dec 2009 - R. Yantosca - Added ProTeX Header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA < 0e0_fp      .or. DENAIR < 1e-30_fp    .or.  &
           RADIUS < 1e-30_fp  .or. SQM  < 1e-30_fp      .or.  &
           STK    < 1e-30_fp  .or.  STKCF  < 1e-30_fp ) THEN

         ! Use default value if any of the above values are zero
         ! This will prevent div-by-zero errors in the eqns below
         ! Value changed from 1d-3 to 1d-30 (bhh, jmao, eam, 7/18/2011)
         ARS_L1K = 1.E-30_fp

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E0_fp/ &
          (SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARS_L1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )

      ENDIF

    END FUNCTION ARSL1K
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Nat
!
! !DESCRIPTION: Subroutine CHECK\_NAT determines whether the solid PSC is
!  composed of ice or NAT (nitric acid trihydrate) (needed for heterogeneous
!  chemistry), or indeed if there is any direct PSC calculation at all. This
!  is important for determining whether to use the JPP or Kirner scheme for
!  ice cloud radii.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_NAT( I, J, L, IS_NAT, IS_PSC, IS_STRAT, &
                            Input_Opt, State_Met, State_Chm )
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I,J,L      ! Grid indices
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
      LOGICAL,        INTENT(OUT) :: IS_NAT     ! Is surface NAT?
      LOGICAL,        INTENT(OUT) :: IS_PSC     ! Are there solid PSCs?
      LOGICAL,        INTENT(OUT) :: IS_STRAT   ! Are we in the strat?
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms
!
! !REVISION HISTORY:
!  17 Apr 2013 - S. D. Eastham - Initial version
!  21 Feb 2014 - M. Sulprizio  - Now pass Input_Opt, State_Met, and State_Chm
!                                objects via the arg list
!  08 Apr 2015 - R. Yantosca   - Remove call to READ_PSC_FILE, this is
!                                now done from DO_CHEMISTRY (chemistry_mod.F)
!  29 Jan 2016 - M. Sulprizio  - Moved this routine from ucx_mod.F to
!                                gckpp_HetRates.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_TROP

      !=================================================================
      ! CHECK_NAT begins here!
      !=================================================================

      ! Check if box is in the troposphere
      IS_TROP  = ( State_Met%InTroposphere(I,J,L) )

      ! Check if box is in the stratosphere
      IS_STRAT = ( State_Met%InStratosphere(I,J,L) .and. ( .not. IS_TROP) )

      ! Check if there are solid PSCs
      IS_PSC   = ( ( Input_Opt%LPSCCHEM ) .and. &
                 ( State_Chm%STATE_PSC(I,J,L) >= 2.0 ) .and. ( IS_STRAT ) )

      ! Check if there is surface NAT
      IS_NAT   = ( ( IS_PSC ) .and. ( SPC_NIT .gt. TINY(1e+0_fp) ) )

    END SUBROUTINE CHECK_NAT
!EOC
  END MODULE GCKPP_HETRATES
