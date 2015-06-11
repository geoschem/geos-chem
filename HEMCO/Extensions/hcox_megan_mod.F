!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_megan_mod.F90
!
! !DESCRIPTION: Module HCOX\_Megan\_Mod contains variables and routines 
!  specifying the algorithms that control the MEGAN inventory of biogenic 
!  emissions (as implemented into the GEOS-Chem model).
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! MEGAN calculates gamma activity factor based upon temperature and 
! radiation information from the past 10 days. In the original GEOS-Chem 
! code, the initial 10-d averages were explicitly calculated during 
! initialization of MEGAN. This is not feasible in an ESMF environment,
! and the following restart variables can now be provided through the 
! HEMCO configuration file:
! \begin{itemize}
! \item T\_DAVG: temperature (10 day average) 
! \item PARDR\_DAVG: direct radiation (10 day average) 
! \item PARDF\_DAVG: diffuse radiation (10 day average)
! \item T\_PREVDAY: temperature (1 day average) 
! \end{itemize}
! These variables are automatically searched for on the first run
! call. If not defined, the current tempeature/radiation value is 
! used instead. Especially for the radiation (PARDR and PARDF), this 
! may result in very different average values - which in turn may 
! result in local differences of >50%. HEMCO emissions are identical 
! to original GEOS-Chem emissions once we pass the first 10 days 
! (ckeller, 10/8/2013).
!\\
!\\
! A similar procedure is also applied to the leaf area index variables.
! The original GEOS-Chem MEGAN code used three LAI variables: current month
! LAI (LAI\_CM), previous month LAI (LAI\_PM), and instantaneous LAI (LAI), 
! which was a daily interpolation of LAI\_CM and next month' LAI, (LAI\_NM). 
! The HEMCO implementation uses only the instantaneous LAI, assuming it is
! updated every day. The previous day LAI is kept in memory and used to
! determine the LAI change over time (used to calculate the gamma leaf age).
! For the first simulation day, the previous' day LAI is taken from the 
! restart file (field LAI\_PREVDAY). If no restart variable is defined, a 
! LAI change of zero is assumed.
! (ckeller, 10/9/2014)
!\\
!\\
! !References:
!
!  \begin{itemize}
!  \item Guenther, A., et al., \emph{The Model of Emissions of Gases and
!        Aerosols from Nature version 2.1 (MEGAN2.1): an extended and updated
!        framework for modeling biogenic emissions}, \underline{Geosci. Model
!        Dev.}, \textbf{5}, 1471-1792, 2012.
!  \item Guenther, A., et al., \emph{A global model of natural volatile 
!        organic compound emissions}, \underline{J.Geophys. Res.}, 
!        \textbf{100}, 8873-8892, 1995.
!  \item Wang, Y., D. J. Jacob, and J. A. Logan, \emph{Global simulation of 
!        tropospheric O3-Nox-hydrocarbon chemistry: 1. Model formulation}, 
!        \underline{J. Geophys. Res.}, \textbf{103}, D9, 10713-10726, 1998.
!  \item Guenther, A., B. Baugh, G. Brasseur, J. Greenberg, P. Harley, L. 
!        Klinger, D. Serca, and L. Vierling, \emph{Isoprene emission estimates 
!        and uncertanties for the Central African EXPRESSO study domain}, 
!        \underline{J. Geophys. Res.}, \textbf{104}, 30,625-30,639, 1999.
!  \item Guenther, A. C., T. Pierce, B. Lamb, P. Harley, and R. Fall, 
!        \emph{Natural emissions of non-methane volatile organic compounds, 
!        carbon monoxide, and oxides of nitrogen from North America}, 
!        \underline{Atmos. Environ.}, \textbf{34}, 2205-2230, 2000.
!  \item Guenther, A., and C. Wiedinmyer, \emph{User's guide to Model of 
!        Emissions of Gases and Aerosols from Nature}. http://cdp.ucar.edu. 
!        (Nov. 3, 2004) 
!  \item Guenther, A., \emph{AEF for methyl butenol}, personal commucation. 
!        (Nov, 2004)
!  \item Sakulyanontvittaya, T., T. Duhl, C. Wiedinmyer, D. Helmig, S. 
!        Matsunaga, M. Potosnak, J. Milford, and A. Guenther, \emph{Monoterpene
!        and sesquiterpene emission estimates for the United States}, 
!        \underline{Environ. Sci. Technol}, \textbf{42}, 1623-1629, 2008.
!  \end{itemize}
!
! !INTERFACE:
!
      MODULE HCOX_MEGAN_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_DIAGN_MOD
      USE HCOX_State_MOD,    ONLY : Ext_State
      USE HCO_STATE_MOD,     ONLY : HCO_STATE

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC  :: HCOX_Megan_Init
      PUBLIC  :: HCOX_Megan_Run
      PUBLIC  :: HCOX_Megan_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_MEGAN_EMISSIONS  ! dbm, new MEGAN driver routine
                                      ! for all compounds (6/21/2012)
      PRIVATE :: UPDATE_HCOT_DAY      
      PRIVATE :: UPDATE_HCOT_15_AVG   
      PRIVATE :: GET_MEGAN_PARAMS
      PRIVATE :: GET_MEGAN_AEF
      PRIVATE :: GET_GAMMA_PAR_PCEEA
      PRIVATE :: GET_GAMMA_T_LI
      PRIVATE :: GET_GAMMA_T_LD
      PRIVATE :: GET_GAMMA_LAI
      PRIVATE :: GET_GAMMA_AGE
      PRIVATE :: GET_GAMMA_SM
      PRIVATE :: CALC_NORM_FAC
      PRIVATE :: SOLAR_ANGLE
      PRIVATE :: FILL_RESTART_VARS 
      PRIVATE :: CALC_AEF
! 
! !REVISION HISTORY: 
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (6/2003).  Updated to 
!        latest algorithm and modified for the standard code by May Fu 
!        (11/2004).
!  (2 ) All emission are currently calculated using TS from DAO met field.
!        TS is the surface air temperature, which should be carefully 
!        distinguished from TSKIN. (tmf, 11/20/2004)
!  (3 ) In GEOS4, the TS used here are the T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  (4 ) Bug fix: change #if block to also cover GCAP met fields (bmy, 12/6/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Bug fix: Skip Feb 29th if GCAP in INIT_MEGAN (phs, 9/18/07)
!  (7 ) Added routine GET_AEF_05x0666 to read hi-res AEF data for the GEOS-5
!        0.5 x 0.666 nested grid simulations (yxw, dan, bmy, 11/6/08)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  09 Mar 2010 - R. Yantosca - Minor bug fix in GET_EMMONOT_MEGAN
!  17 Mar 2010 - H. Pye      - AEF_SPARE must be a scalar local variable
!                              in GET_EMMONOT_MEGAN for parallelization.
!  20 Aug 2010 - R. Yantosca - Move CMN_SIZE to top of module
!  20 Aug 2010 - R. Yantosca - Now set DAY_DIM = 24 for MERRA, since the
!                              surface temperature is now an hourly field.
!  01 Sep 2010 - R. Yantosca - Bug fix in INIT_MEGAN: now only read in 
!                              NUM_DAYS (instead of 15) days of sfc temp data
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids 
!  06 Dec 2011 - E. Fischer  - Added Acetone emissions 
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Use updated GET_LOCALTIME from time_mod.F
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  13 Aug 2013 - M. Sulprizio- Modifications for updated SOA sim (H. Pye):
!                               Add sesquiterpenes to MEGAN group;
!                               Add plant functional types (PFT_xx);
!                               Rename GET_EMMONOG_MEGAN to GET_EMTERP_MEGAN;
!                               Add routines READ_PFT and GET_AEF_GEN
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  05 Oct 2013 - C. Keller   - Now a HEMCO extension 
!  04 Aug 2014 - C. Keller   - Added 'manual' diagnostics for Acetone.
!  09 Oct 2014 - C. Keller   - Now use only GC_LAI (keep prev. LAI in memory)
!  22 Dec 2014 - C. Keller   - Now use flexible precision (hp) everywhere.
!                              Option to read temperature/irradiation from
!                              restart.
!  26 Jan 2015 - M. Sulprizio- Update from D. Millet (19 Jan 2013): Streamlined
!                              computations into a single driver routine
!                              and updated emissions according to MEGAN 2.1 as
!                              described in:
!                              Guenther et al., The Model of Emissions of Gases
!                              and Aerosols from Nature version 2.1 (MEGAN2.1):
!                              an extended and updated framework for modeling 
!                              biogenic emissions, GMD, 5, 1471-1492, 2012.
!  12 Feb 2015 - M. Sulprizio- Remove GET_AEF_GEN routine. We now calculate AEFs
!                              for FARN, BCAR, and OSQT in CALC_AEF using
!                              parameters from Guenther et al., 2012.
!  18 Feb 2015 - M. Sulprizio- Remove LPECCA logical flag since we use this
!                              scheme exclusively now.
!                              Restore emissions of individual MEGAN species to
!                              diagnostics for consistency with pre-HEMCO code.
!  10 Jun 2015 - M. Sulprizio- Bug fix for SOA simulation: Now convert AEFs for
!                              sesquiterpenes to kg/m2/s.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
! 
      INTEGER                 :: ExtNr
      INTEGER                 :: ExtNrMono 
      INTEGER                 :: ExtNrSoa
      REAL(hp)                :: ISOP_SCALING

      ! Physical parameter
      REAL(hp)                :: D2RAD  ! Degrees to radians 
      REAL(hp)                :: RAD2D  ! Radians to degrees

      ! HEMCO species IDs
      INTEGER                 :: IDTISOP
      INTEGER                 :: IDTMONX
      INTEGER                 :: IDTACET
      INTEGER                 :: IDTCO 
      INTEGER                 :: IDTPRPE
      INTEGER                 :: IDTC2H4
      INTEGER                 :: IDTALD2
      INTEGER                 :: IDTOCPI
      INTEGER                 :: IDTLIMO
      INTEGER                 :: IDTMTPA
      INTEGER                 :: IDTMTPO
      INTEGER                 :: IDTSESQ

      ! Pointers to annual emission factor arrays.
      ! These fields are obtained from ext. data (from config file)
      ! These arrays represent the base emission values in kg(C)/m2/s
      ! that will be scaled based upon meteorological conditions.
      REAL(sp), POINTER       :: AEF_ISOP(:,:)  => NULL()   ! Isoprene 
      REAL(sp), POINTER       :: AEF_MBOX(:,:)  => NULL()   ! Methyl butenol
      REAL(sp), POINTER       :: AEF_BPIN(:,:)  => NULL()   ! Beta-pinene
      REAL(sp), POINTER       :: AEF_CARE(:,:)  => NULL()   ! 3-Carene
      REAL(sp), POINTER       :: AEF_LIMO(:,:)  => NULL()   ! Limonene
      REAL(sp), POINTER       :: AEF_OCIM(:,:)  => NULL()   ! Ocimene
      REAL(sp), POINTER       :: AEF_SABI(:,:)  => NULL()   ! Sabinene
      REAL(sp), POINTER       :: GEIA_ORVC(:,:) => NULL()   ! Organic carbon

      ! Annual emission factors arrays
      ! These fields are not read from file, but are computed in CALC_AEF
      REAL(hp), ALLOCATABLE   :: AEF_APIN(:,:)              ! Alpha-pinene
      REAL(hp), ALLOCATABLE   :: AEF_MYRC(:,:)              ! Myrcene
      REAL(hp), ALLOCATABLE   :: AEF_OMON(:,:)              ! Other monoterpenes
      REAL(hp), ALLOCATABLE   :: AEF_MOHX(:,:)              ! Methanol
      REAL(hp), ALLOCATABLE   :: AEF_ACET(:,:)              ! Acetone
      REAL(hp), ALLOCATABLE   :: AEF_ETOH(:,:)              ! Ethanol
      REAL(hp), ALLOCATABLE   :: AEF_CH2O(:,:)              ! Formaldehyde
      REAL(hp), ALLOCATABLE   :: AEF_ALD2(:,:)              ! Acetaldehyde
      REAL(hp), ALLOCATABLE   :: AEF_FAXX(:,:)              ! Formic acid
      REAL(hp), ALLOCATABLE   :: AEF_AAXX(:,:)              ! Acetic acid
      REAL(hp), ALLOCATABLE   :: AEF_C2H4(:,:)              ! Ethene
      REAL(hp), ALLOCATABLE   :: AEF_TOLU(:,:)              ! Toluene
      REAL(hp), ALLOCATABLE   :: AEF_HCNX(:,:)              ! HCN
      REAL(hp), ALLOCATABLE   :: AEF_PRPE(:,:)              ! >= C3 alkenes
      REAL(hp), ALLOCATABLE   :: AEF_FARN(:,:)              ! a-Farnesene
      REAL(hp), ALLOCATABLE   :: AEF_BCAR(:,:)              ! b-Caryophyllene
      REAL(hp), ALLOCATABLE   :: AEF_OSQT(:,:)              ! Other sesquiterp.

      ! Pointers to CLM4 plant functional types
      REAL(sp), POINTER       :: PFT_BARE(:,:)                => NULL()
      REAL(sp), POINTER       :: PFT_NDLF_EVGN_TMPT_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_NDLF_EVGN_BORL_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_NDLF_DECD_BORL_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_EVGN_TROP_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_EVGN_TMPT_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_DECD_TROP_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_DECD_TMPT_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_DECD_BORL_TREE(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_EVGN_SHRB(:,:)      => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_DECD_TMPT_SHRB(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_BDLF_DECD_BORL_SHRB(:,:) => NULL()
      REAL(sp), POINTER       :: PFT_C3_ARCT_GRSS(:,:)        => NULL()
      REAL(sp), POINTER       :: PFT_C3_NARC_GRSS(:,:)        => NULL()
      REAL(sp), POINTER       :: PFT_C4_GRSS(:,:)             => NULL()
      REAL(sp), POINTER       :: PFT_CROP(:,:)                => NULL()

      ! Normalization factor
      REAL(hp), ALLOCATABLE   :: NORM_FAC(:)

      ! Scalars
      INTEGER,  PARAMETER  :: DAY_DIM        = 24       ! # of 1-hr periods/day
!      INTEGER, PARAMETER  :: DAY_DIM        = 8        ! # of 3-hr periods/day
      INTEGER,  PARAMETER  :: NUM_DAYS       = 10       ! # of days to avg 
      REAL(hp), PARAMETER  :: WM2_TO_UMOLM2S = 4.766_hp ! W/m2 -> umol/m2/s

      ! Past light & temperature conditions (mpb,2009)
      ! (1) Temperature at 2m (TS):
!      REAL(hp), ALLOCATABLE, TARGET :: HCOT_DAILY(:,:)       ! Daily averaged sfc temp
      REAL(sp), POINTER             :: T_PREVDAY(:,:)        ! Pointer to previous day T.
      REAL(sp), ALLOCATABLE, TARGET :: HCOT_DAY(:,:,:)       ! Holds 1 day of sfc temp data
      REAL(sp), ALLOCATABLE, TARGET :: HCOT_15(:,:,:)        ! Holds 15 days of daily avg T
      REAL(sp), ALLOCATABLE, TARGET :: HCOT_15_AVG(:,:)      ! Sfc temp avg'd over NUM_DAYS

      ! (2) PAR Direct: 
!      REAL(hp), ALLOCATABLE, TARGET :: HCOPARDR_DAILY(:,:)   ! Average daily PARDR
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDR_DAY(:,:,:)   ! Holds 1 day of PARDR data
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDR_15(:,:,:)    ! 10 days of daily avg'd PARDR
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDR_15_AVG(:,:)  ! PARDR averaged over NUM_DAYS

      ! (3) PAR Diffuse: 
!      REAL(hp), ALLOCATABLE, TARGET :: HCOPARDF_DAILY(:,:)   ! Average daily PARDR
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDF_DAY(:,:,:)   ! Holds 1-day of PARDR data
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDF_15(:,:,:)    ! 10 days of daily avg'd PARDR 
      REAL(sp), ALLOCATABLE, TARGET :: HCOPARDF_15_AVG(:,:)  ! PARDF averaged over NUM_DAYS

      ! previous LAI values (ckeller, 10/09/2014)
      REAL(sp), ALLOCATABLE, TARGET :: LAI_PREVSTEP(:,:)     ! LAI of prev. time step
      REAL(sp), ALLOCATABLE, TARGET :: LAI_PREVDAY(:,:)      ! LAI of prev. day

      ! Days between mid-months (updated by HEMCO clock)
      INTEGER             :: DAYS_BTW_M

      ! testing only
      integer, parameter  :: ix = 20 !20 !25 !13 !20
      integer, parameter  :: iy = 43 !43 !22 !38 !31

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Megan_Run
!
! !DESCRIPTION: Subroutine HCOX\_MEGAN\_Run is the driver routine for
! the MEGAN model within the new emissions structure. 
! Note that all diagnostics are commented since those are still written
! as part of the old emission structure.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_Megan_Run( am_I_Root, ExtState, HcoState, RC ) 
!
! !USES:
!
      USE HCO_FLUXARR_MOD,     ONLY : HCO_EmisAdd 
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_First 
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_Rewind
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewHour 
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewDay
      USE HCO_Restart_Mod,     ONLY : HCO_RestartWrite
!
! !INPUT PARAMETERS: 
!
      LOGICAL,         INTENT(IN)    :: am_I_Root
      TYPE(Ext_State), POINTER       :: ExtState 
      TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY: 
!  05 Aug 2013 - C. Keller   - Initial version 
!  08 Oct 2014 - C. Keller   - Now use HEMCO clock to get days between months
!  09 Oct 2014 - C. Keller   - Now use only GC_LAI. This makes days between months
!                              obsolete
!  16 Oct 2014 - C. Keller   - Initialize flux arrays to avoid float invalid
!   9 Mar 2015 - R. Yantosca - Bug fix: add EMIS_ALD2 to $OMP PRIVATE clause
!  08 May 2015 - C. Keller   - Now read/write restart variables from here to
!                              accomodate replay runs in GEOS-5.
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J
      REAL(hp)            :: EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN
      REAL(hp)            :: EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE
      REAL(hp)            :: EMIS_OCIM, EMIS_OMON, EMIS_MONO, EMIS_ALD2
      REAL(hp)            :: EMIS_MOHX, EMIS_ETOH, EMIS_FAXX, EMIS_AAXX
      REAL(hp)            :: EMIS_ACET, EMIS_PRPE, EMIS_C2H4
      REAL(hp)            :: EMIS_FARN, EMIS_BCAR, EMIS_OSQT
      REAL(hp)            :: BIOG_ALPH, BIOG_LIMO
      REAL(hp)            :: BIOG_ALCO
      REAL(hp)            :: MONO_MOL
      REAL(hp)            :: MBOX_MOL
      REAL(hp)            :: X, Y 
      REAL(hp)            :: DIURFACT, DIURORVC
      REAL(hp)            :: TMP
      LOGICAL             :: ERR
      LOGICAL             :: FIRST, IsNewHour, IsNewDay 

      ! Emission arrays (kgC/m2/s)
      REAL(hp), TARGET    :: FLUXISOP  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMONO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETmo(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETmb(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXACETbg(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXCO    (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXPRPE  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXC2H4  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOCPI  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXLIMO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMTPA  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMTPO  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSESQ  (HcoState%NX,HcoState%NY)

      ! Emission arrays for use in diagnostics only (kgC/m2/s)
      REAL(hp), TARGET    :: FLUXMBOX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXAPIN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXBPIN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXSABI  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMYRC  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXCARE  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOCIM  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOMON  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXALD2  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXMOHX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXETOH  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXFAXX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXAAXX  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXFARN  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXBCAR  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXOSQT  (HcoState%NX,HcoState%NY)

      ! For diagnostics
      REAL(hp), POINTER   :: Arr2D(:,:) => NULL()
      CHARACTER(LEN=63)   :: DiagnName

      ! Molecules C / kg C
      REAL(hp), PARAMETER   :: XNUMOL_C = 6.022e+23_hp / 12e-3_hp
      REAL(hp), PARAMETER   :: UNITCONV = 1.0_hp

      ! Fraction of yield of OC (hydrophilic) from terpene emissions
      REAL(hp), PARAMETER   :: EMMT2OC = 1.0e-1_hp

      ! Conversion factors for OC/BC
      REAL(hp), PARAMETER   :: FC1 = 136.2364_hp / 120.11_hp
      REAL(hp), PARAMETER   :: FC2 = 154.2516_hp / 120.11_hp
      REAL(hp), PARAMETER   :: FC3 = 204.3546_hp / 180.165_hp
      REAL(hp), PARAMETER   :: FC4 = 152.0_hp    / 120.11_hp

      ! Conversion factors for acetone calculations
      REAL(hp), PARAMETER   :: YIELD_MO   = 0.116_hp
      REAL(hp), PARAMETER   :: MONO_SCALE = 0.89_hp
      REAL(hp), PARAMETER   :: MB_SCALE1  = 0.6_hp
      REAL(hp), PARAMETER   :: MB_SCALE2  = 0.76_hp

      !=================================================================
      ! HCOX_Megan_Run begins here! 
      !================================================================= 

      ! Enter 
      CALL HCO_ENTER ( 'HCOX_Megan_Run (hcox_megan_mod.F)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      ERR = .FALSE.

      ! Time information
      FIRST     = HcoClock_First  ( .TRUE. )
      IsNewHour = HcoClock_NewHour( .TRUE. )
      IsNewDay  = HcoClock_NewDay ( .TRUE. )

      ! Initialize flux arrays
      FLUXISOP   = 0.0_hp
      FLUXMBOX   = 0.0_hp
      FLUXAPIN   = 0.0_hp
      FLUXBPIN   = 0.0_hp
      FLUXSABI   = 0.0_hp
      FLUXMYRC   = 0.0_hp
      FLUXCARE   = 0.0_hp
      FLUXOCIM   = 0.0_hp
      FLUXOMON   = 0.0_hp
      FLUXMONO   = 0.0_hp
      FLUXALD2   = 0.0_hp
      FLUXMOHX   = 0.0_hp
      FLUXETOH   = 0.0_hp
      FLUXFAXX   = 0.0_hp
      FLUXAAXX   = 0.0_hp
      FLUXACETmo = 0.0_hp
      FLUXACETmb = 0.0_hp
      FLUXACETbg = 0.0_hp
      FLUXCO     = 0.0_hp
      FLUXPRPE   = 0.0_hp
      FLUXC2H4   = 0.0_hp
      FLUXOCPI   = 0.0_hp
      FLUXLIMO   = 0.0_hp
      FLUXMTPA   = 0.0_hp
      FLUXMTPO   = 0.0_hp
      FLUXFARN   = 0.0_hp
      FLUXBCAR   = 0.0_hp
      FLUXOSQT   = 0.0_hp
      FLUXSESQ   = 0.0_hp

      IF ( FIRST ) THEN
         ! Generate annual emission factors for MEGAN inventory
         CALL CALC_AEF( am_I_Root, HcoState, ExtState, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Calculate normalization factor (dbm, 11/2012)
         CALL CALC_NORM_FAC( am_I_Root, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Fill restart variables. This has to be done on the first call
      ! or after rewinding the clock (in replay simulations). 
      IF ( FIRST .OR. HcoClock_Rewind( .TRUE. ) ) THEN
         CALL FILL_RESTART_VARS( am_I_Root, HcoState, ExtState, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

      ! On a new day, copy previous LAI values to array that holds
      ! previous day LAI.
      ELSEIF ( IsNewDay ) THEN
         LAI_PREVDAY = LAI_PREVSTEP
      ENDIF

      ! Set to 1 day since we know that LAI is updated every day.
      DAYS_BTW_M = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J,      EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN   )
!$OMP+PRIVATE( EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE, EMIS_OCIM   )
!$OMP+PRIVATE( EMIS_OMON, EMIS_MONO, EMIS_MOHX, EMIS_ETOH, EMIS_FAXX   )
!$OMP+PRIVATE( EMIS_AAXX, EMIS_ACET, EMIS_PRPE, EMIS_C2H4              )
!$OMP+PRIVATE( EMIS_FARN, EMIS_BCAR, EMIS_OSQT, EMIS_ALD2              )
!$OMP+PRIVATE( MONO_MOL,  TMP,       MBOX_MOL                          )
!$OMP+PRIVATE( DIURFACT,  DIURORVC,  BIOG_ALPH                         ) 
!$OMP+PRIVATE( BIOG_LIMO, BIOG_ALCO, X, Y, RC                          )

      !-----------------------------------------------------------------
      ! Loop over all grid boxes 
      !-----------------------------------------------------------------
      DO J = 1, HcoState%NY 
      DO I = 1, HcoState%NX

         ! Zero biogenic species
         EMIS_ISOP  = 0.0_hp
         EMIS_MBOX  = 0.0_hp
         EMIS_APIN  = 0.0_hp
         EMIS_BPIN  = 0.0_hp
         EMIS_LIMO  = 0.0_hp
         EMIS_SABI  = 0.0_hp
         EMIS_MYRC  = 0.0_hp
         EMIS_CARE  = 0.0_hp
         EMIS_OCIM  = 0.0_hp
         EMIS_OMON  = 0.0_hp
         EMIS_MONO  = 0.0_hp
         EMIS_ALD2  = 0.0_hp
         EMIS_MOHX  = 0.0_hp
         EMIS_ETOH  = 0.0_hp
         EMIS_FAXX  = 0.0_hp
         EMIS_AAXX  = 0.0_hp
         EMIS_ACET  = 0.0_hp
         EMIS_PRPE  = 0.0_hp
         EMIS_C2H4  = 0.0_hp
         EMIS_FARN  = 0.0_hp
         EMIS_BCAR  = 0.0_hp
         EMIS_OSQT  = 0.0_hp

         !-----------------------------------------------------------------
         ! Update internal temperature arrays 
         !-----------------------------------------------------------------
   
         ! Arrays holding hourly temperature bins. Update every hour. 
         ! For now, all historical values are set to the current values
         ! on first call. This avoids reading all met-fields of the past
         ! 10 days, but introduces errors in the MEGAN emissions of the
         ! first 10 days.
         IF ( IsNewHour .AND. .NOT. FIRST ) THEN 
            CALL UPDATE_HCOT_DAY ( HcoState, ExtState, I, J, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( 'UPDATE_HCOT_DAY', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
         ENDIF
   
         ! Arrays holding the daily temperatures. Update every new day
         ! or on first call. 
         ! This routine relies on the arrays filled in UPDATE_HCOT_DAY, so
         ! make sure it's called afterwards!
         IF ( IsNewDay .AND. .NOT. FIRST ) THEN
            CALL UPDATE_HCOT_15_AVG ( HcoState, I, J, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( 'UPDATE_HCOT_15_AVG', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
         ENDIF

         !--------------------------------------------------------------
         ! Calculate VOC emissions 
         !
         ! The GET_EMIS*_MEGAN calls now use the annual scale factors
         ! imported through the HEMCO list, which are already in units 
         ! of kg(C)/m2/s. Thus, no further unit conversion is required 
         ! anymore (--> UNITCONV = 1). ckeller, 14/01/25.
         !
         ! Updated to new MEGAN routine (dbm, 12/2012)
         !--------------------------------------------------------------

         !--------------------------------------------------------------------
         ! MEGAN Isoprene
         !--------------------------------------------------------------------

         ! Isoprene [kg C/m2/s]
         CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                             'ISOP',   EMIS_ISOP, RC )
         IF ( RC /= HCO_SUCCESS ) THEN 
            CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_ISOP', RC )
            ERR = .TRUE.
            EXIT
         ENDIF 

         !FP_ISOP (6/2009)
         EMIS_ISOP = ISOP_SCALING * EMIS_ISOP

         ! Write isoprene emissions into tracer tendency array
         IF ( IDTISOP > 0 ) THEN
            FLUXISOP(I,J) = EMIS_ISOP 

!            ! testing only
!            if ( i==ix .and. j==iy ) then
!               write(*,*) 'HEMCO MEGAN @ ', ix, iy
!               write(*,*) 'ISOP_SCAL      : ', ISOP_SCALING
!               write(*,*) 'ISOP (kgC/m2/s): ', EMIS_ISOP
!            endif

         ENDIF

         !--------------------------------------------------------------------
         ! MEGAN monoterpenes
         !--------------------------------------------------------------------
         ! Updated the calls below for new MEGAN routine (dbm, 12/2012)

         IF ( ExtNrMono > 0 ) THEN

            ! ---------------------------------------------------
            ! Methyl butenol
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'MBOX',   EMIS_MBOX, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_MBOX', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXMBOX(I,J) = EMIS_MBOX

            ! ---------------------------------------------------
            ! Alpha Pinene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'APIN',   EMIS_APIN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_APIN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXAPIN(I,J) = EMIS_APIN

            ! ---------------------------------------------------
            ! Beta Pinene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'BPIN',   EMIS_BPIN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_BPIN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXBPIN(I,J) = EMIS_BPIN

            ! ---------------------------------------------------
            ! Limonene emissions 
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'LIMO',   EMIS_LIMO, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_LIMO', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            IF ( ExtNrSoa < 1 ) THEN
               ! Update FLUXLIMO here for diagnostic purposes if the
               ! SOA mechanism is turned off.
               ! Otherwise, FLUXLIMO is updated below.
               FLUXLIMO(I,J) = EMIS_LIMO
            ENDIF

            ! ---------------------------------------------------
            ! Sabinene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'SABI',   EMIS_SABI, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_SABI', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXSABI(I,J) = EMIS_SABI

            ! ---------------------------------------------------
            ! Mycrene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'MYRC',   EMIS_MYRC, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_MYRC', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXMYRC(I,J) = EMIS_MYRC

            ! ---------------------------------------------------
            ! 3-Carene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'CARE',   EMIS_CARE, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_CARE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXCARE(I,J) = EMIS_CARE

            ! ---------------------------------------------------
            ! Ocimene emissions
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'OCIM',   EMIS_OCIM, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_OCIM', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXOCIM(I,J) = EMIS_OCIM

            ! ---------------------------------------------------
            ! Other monoterpenes
            ! (added 12/2012; dbm)
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'OMON',   EMIS_OMON, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_OMON', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 
            FLUXOMON(I,J) = EMIS_OMON

            ! ---------------------------------------------------
            ! Total monoterpenes = sum of individual
            ! dbm, now add other lumped monoterpenes (11/2012)
            EMIS_MONO = EMIS_APIN + EMIS_BPIN + EMIS_LIMO + EMIS_SABI + 
     &                  EMIS_MYRC + EMIS_CARE + EMIS_OCIM + EMIS_OMON

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXMONO(I,J) = EMIS_MONO

         ELSE
            EMIS_MONO = 0.0_hp
            EMIS_MBOX = 0.0_hp
         ENDIF


!         ! testing only
!         if ( i == ix .AND. j == iy .and. .FALSE. ) then
!                 write(*,*) ' '
!                 write(*,*) 'HEMCO MEGAN (kg/m2/s) @ ', i, j 
!                 write(*,*) 'ISOP: ', EMIS_ISOP
!                 write(*,*) 'MBOX: ', EMIS_MBOX
!                 write(*,*) 'MONO: ', EMIS_MONO
!                 write(*,*) 'APIN: ', EMIS_APIN
!                 write(*,*) 'BPIN: ', EMIS_BPIN
!                 write(*,*) 'LIMO: ', EMIS_LIMO
!                 write(*,*) 'SABI: ', EMIS_SABI
!                 write(*,*) 'MYRC: ', EMIS_MYRC
!                 write(*,*) 'CARE: ', EMIS_CARE
!                 write(*,*) 'OCIM: ', EMIS_OCIM
!                 write(*,*) 'OMON: ', EMIS_OMON
!         endif

         !--------------------------------------------------------------------
         ! MEGAN Acetaldehyde
         !--------------------------------------------------------------------
         IF ( IDTALD2 > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'ALD2',   EMIS_ALD2, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_ALD2', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXALD2(I,J) = EMIS_ALD2
         ENDIF

         !--------------------------------------------------------------------
         ! Other MEGAN biogenics
         ! Calls included here for future incorporation or
         ! specialized simulations. (dbm, 12/2012)
         !--------------------------------------------------------------------

         ! ---------------------------------------------------
         ! MEGAN Methanol
         ! Uncomment this to compute biogenic methanol emissions
!         CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
!     &                             'MOHX',   EMIS_MOHX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN 
!            CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_MOHX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF 
!         FLUXMOHX(I,J) = EMIS_MOHX

         ! ---------------------------------------------------
         ! MEGAN Ethanol
         ! Uncomment this to compute biogenic ethanol emissions
!         CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
!     &                             'ETOH',   EMIS_ETOH, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN 
!            CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_ETOH', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXETOH(I,J) = EMIS_ETOH

         ! ---------------------------------------------------
         ! MEGAN Formic acid
         ! Uncomment this to compute biogenic formic acid emissions
!         CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
!     &                             'FAXX',   EMIS_FAXX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN 
!            CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_FAXX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXFAXX(I,J) = EMIS_FAXX

         ! ---------------------------------------------------
         ! MEGAN Acetic acid
         ! Uncomment this to compute biogenic acetic acid emissions
!         CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
!     &                             'AAXX',   EMIS_AAXX, RC )
!         IF ( RC /= HCO_SUCCESS ) THEN 
!            CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_AAXX', RC )
!            ERR = .TRUE.
!            EXIT
!         ENDIF
!         FLUXAAXX(I,J) = EMIS_AAXX

         !--------------------------------------------------------------------
         ! Biogenic acetone emissions
         !--------------------------------------------------------------------
         IF ( IDTACET > 0 ) THEN

            !-----------------------------------------------------------------
            ! (1) BIOGENIC EMISSIONS OF ACETONE FROM MONOTERPENES
            ! Monoterpenes has same # molecules/kg of carbon as isoprene
            ! The yield for monoterpenes is .12 mol/mol from Reisell et.al. 
            ! 1999 (this does not includes direct acetone emissions)

            ! Convert [kg C/m2/s] to [kg MONOTERPENE/m2/s]
            ! There are 10 C atoms per molecule of MONOTERPENE
            MONO_MOL = EMIS_MONO / 10.0_hp
      
            ! Apply yield from monoterpenes to get [kg ACET/m2/s]
            TMP      = MONO_MOL * YIELD_MO
      
            ! Convert acetone emissions back into [kg C/m2/s] 
            TMP      = TMP * 3.0_hp
      
            ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
            TMP      = TMP * MONO_SCALE
      
            ! Add to total biogenic acetone emissions
            FLUXACETmo(I,J) = TMP

            !-----------------------------------------------------------------
            ! (2) BIOGENIC ACETONE FROM METHYL BUTENOL -- NORTH AMERICA
            !
            ! Methyl Butenol (a.k.a. MBO) produces acetone with a molar yield 
            ! of 0.6 [Alvarado (1999)].  The biogenic source of MBO is thought 
            ! to be restricted to North America.  According to Guenther (1999) 
            ! North america emits 3.2Tg-C of MBO, producing 1.15 Tg-C of 
            ! Acetone in North America.
            !=================================================================
            TMP = 0.0_hp
      
            ! Lon and lat of grid box (I,J) in degrees
            X = HcoState%Grid%XMID%Val( I, J )
            IF ( X >= 180.0_hp ) X = X - 360.0_hp
            Y = HcoState%Grid%YMID%Val( I, J )
      
            ! Methyl butenol is emitted only in North America, where
            ! ( -167.5 <= lon <= -52.5 ) and ( 16.0 <= lat <= 72.0 ) 
            IF ( ( X >= -167.5_hp .and. X <= -52.5_hp ) .AND.
     &           ( Y >=   16.0_hp .and. Y <=  72.0_hp ) ) THEN
      
               ! Convert from [kg C/m2/s] to [kg MBO/m2/s] 
               ! There are 5 C atoms per molecule MBO
               MBOX_MOL = EMIS_MBOX / 5.0_hp
      
               ! Apply yield from MBO to get [kg ACET/m2/s]
               TMP      = MBOX_MOL * MB_SCALE1
      
               ! Convert from [kg ACET/m2/s] to [kg C/m2/s]
               ! There are 3 C atoms per acetone molecule
               TMP      = TMP * 3.0_hp
      
               ! Scale to a posteriori source from Jacob et al 2001 (bdf,
               ! 9/5/01)
               FLUXACETmb(I,J) = TMP * MB_SCALE2

            ENDIF

            !-----------------------------------------------------------------
            ! (3) BIOGENIC ACETONE -- DIRECT EMISSION 
            ! evf, removed obsolete code, replaced with MEGAN acetone
            ! emissions (5/25/2011) Direct Emission now includes emission
            ! from grasses and emission from dry leaf matter
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'ACET',   EMIS_ACET, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_ACET', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! Write to array
            FLUXACETbg(I,J) = EMIS_ACET
         ENDIF

         !--------------------------------------------------------------------
         !  Biogenic source of CO -- from oxidation of METHANOL and MONOTERPENES
         !
         !  CO from METHANOL oxidation -- scaled from ISOPRENE (bnd, 1/2/01)
         !
         !    We need to scale the Isoprene flux to get the CH3OH (methanol)
         !    flux. Currently, the annual isoprene flux in GEOS-CHEM is
         !    ~ 397 Tg C.
         !
         !    Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH
         !    oxidation based on Singh et al. 2000 [JGR 105, 3795-3805] who
         !    estimate a global methanol source of 122 Tg yr-1, of which most
         !    (75 Tg yr-1) is "primary biogenic". He also recommends for now
         !    that the CO flux from CH3OH oxidation be scaled to monthly mean
         !    isoprene flux.
         !
         !    To get CO from METHANOL oxidation, we must therefore multiply
         !    the ISOPRENE flux by the following scale factor:
         !    ( 100 Tg CO from CH3OH Oxidation / 397 Tg C from Isoprene Flux ) *
         !    (  12 g C/mole                   / 28 g CO/mole                )
         !
         !  CO from MONOTERPENE oxidation (bnd, bmy, 1/2/01)
         !
         !    Assume the production of CO from monoterpenes is instantaneous
         !    even though the lifetime of intermediate species may be on the
         !    order of hours or days. This assumption will likely cause CO from
         !    monoterpene oxidation to be too high in the box in which the
         !    monoterpene is emitted.
         !
         !    The CO yield here is taken from:
         !
         !    Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
         !      "The ultimate yield of CO from the tropospheric oxidation of
         !       terpenes (including both O3 and OH reactions) was estimated
         !       to be 20% on the carbon number basis." They studied
         !       ALPHA- & BETA-pinene.
         !
         !    Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
         !      "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
         !
         !    EMIS_MONO is now kgC, conversion factor then becomes: 
         !    28g CO / 60g C = 46.666%
         !    (ckeller, 4/15/14).
         !--------------------------------------------------------------------
 
         ! NOTE: Don't save into EMISRR for the tagged CO simulation.  
         ! Also for tagged CO we don't use monoterpenes  ??????
         ! (jaf, mak, bmy, 2/14/08)
         IF ( (IDTCO>0) .AND. (ExtNrMono>0) ) THEN
!            FLUX = EMMO * 0.2d0
            FLUXCO(I,J) = EMIS_MONO * 0.4666_hp

!            ! testing only
!            if ( i==ix .and. j==iy .and. .FALSE. ) then
!               write(*,*) ' '
!               write(*,*) 'HEMCO MEGAN @ ', ix, iy
!               write(*,*) 'CO (kg/m2/s): ', FLUXCO(I,J)
!            endif

         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of PRPE
         !
         ! Now uses MEGAN2.1 (dbm, 12/2012)
         !--------------------------------------------------------------------
         IF ( IDTPRPE > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'PRPE',   EMIS_PRPE, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_PRPE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXPRPE(I,J) = EMIS_PRPE

!           ! testing only
!           if ( i==ix .and. j==iy .and. .FALSE. ) then
!              write(*,*) ' '
!              write(*,*) 'HEMCO MEGAN @ ', ix, iy
!              write(*,*) 'PRPE (kg/m2/s): ', FLUXPRPE(I,J)
!           endif

         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of ethene (C2H4)
         !
         ! Now uses MEGAN2.1 (dbm, 12/2012)
         !--------------------------------------------------------------------
         IF ( IDTC2H4 > 0 ) THEN
            CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                'C2H4',   EMIS_C2H4, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_C2H4', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! Add to tracer tendency array [kg C/m2/s]
            FLUXC2H4(I,J) = EMIS_C2H4
         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of OC
         ! Assume a 10% yield of hydrophilic OC from terpene emissions
         ! if secondary organic aerosols are turned off. Otherwise,
         ! use CALTECH algorithm to calculate biogenic VOC emissions.
         !--------------------------------------------------------------------

         ! Get monoterpenes from MEGAN [kg C/m2/s]
         ! terpene emissions
         IF ( (ExtNrMono > 0 ) .AND. ( IDTOCPI > 0 ) ) THEN           

            ! Convert monoterpenes into hydrophilic OC and pass to
            ! OC emissions array 
            FLUXOCPI(I,J) = EMIS_MONO * EMMT2OC

!            ! testing only
!            if ( i==ix .and. j==iy ) then
!                  write(*,*) ' '
!                  write(*,*) 'HEMCO MEGAN @ ', ix, iy
!                  write(*,*) 'OCPI (NOSOA) (kg/m2/s): ', FLUXOCPI(I,J)
!                  write(*,*) 'EMIS_MONO             : ', EMIS_MONO
!                  write(*,*) 'EMMT2OC               : ', EMMT2OC
!            endif
         ENDIF

         ! Use CALTECH algorithm if SOA is turned on
         IF ( ExtNrSoa > 0 ) THEN
           
            ! Error trap: check if SZAFACT is indeed defined
            IF ( .NOT. ASSOCIATED(ExtState%SZAFACT%Arr%Val) ) THEN
                  CALL HCO_ERROR('SZAFACT not set - needed for SOA', RC)
                  ERR = .TRUE.
                  EXIT
            ENDIF
 
            ! Impose a diurnal variation on NVOC during the day
            DIURFACT = ExtState%SZAFACT%Arr%Val(I,J)
            DIURORVC = DIURFACT * GEIA_ORVC(I,J) 

            ! ----------------------------------------------------------------
            ! The new MEGAN implementation has speciated information
            ! (hotp 3/7/10)
            !! For GCAP Meteorology year 2000 in Tg/yr:
            !! ------------------------------
            !! HC Class  New MEGAN  Old MEGAN
            !! --------  ---------  ---------
            !! ALPH        84         92
            !! LIMO        10         27
            !! TERP         3.2        3.5
            !! ALCO        47         38
            !! SESQ        15         15    (no change for SESQ)
            !!           -----      -----
            !! TOTAL      159        176
            ! ------------------------------
            ! updated again to include MEGAN sesquiterpenes 
            ! and other monoterpenes (hotp 7/28/10)
            ! as of 7/28/10 for year 2000 GEOS4 2x2.5 in Tg/yr:
            ! ------------------------------
            ! HC Class  New MEGAN  Old MEGAN
            ! --------  ---------  ---------
            ! MTPA        73         86
            ! LIMO        10         25
            ! MTPO        33          3.2+38 
            ! SESQ        13         15
            !           -----      --------
            ! TOTAL      129        167.2
            ! see Pye et al. 2010 ACP
            ! ----------------------------------------------------------------

            IF ( ExtNrMono > 0 ) THEN

               ! MEGAN emission update for SESQ and other monoterps
               ! (hotp 3/10/10)
               ! terpenoid ketones no longer lumped into ALPH
               ! terpenoid alcohols no longer lumped in ALCO
               ! Other monoterpenes (OMON from MEGAN) all placed in TERP
               ! All sesq lumped together into SESQ

               ! new mtp lumping (as of 5/20/10) hotp
               ! GEOS-Chem   MEGAN
               ! =========   ==========================================
               ! MTPA        a-pinene   (APIN), b-pinene  (BPIN), 
               !             sabinene   (SABI), carene    (CARE)
               ! LIMO        limonene   (LIMO)
               ! MTPO        myrcene    (MYRC), ocimene   (OCIM),
               !             other mono (OMON)
               ! SESQ        farnesene  (FARN), b-caryoph (BCAR), 
               !             other sesq (OSQT)
               ! =========   ==========================================

               !--------------------------------------------------------------
               ! MEGAN MTPA 
               !--------------------------------------------------------------
               ! MTPA=a-,b-pinene,sabinene,carene (hotp 5/20/10)
               IF ( IDTMTPA > 0 ) THEN
                  FLUXMTPA(I,J) = ( EMIS_APIN + EMIS_BPIN +
     &                              EMIS_SABI + EMIS_CARE ) * FC1
               ENDIF

               !--------------------------------------------------------------
               ! MEGAN Limonene 
               !--------------------------------------------------------------
               ! [kg C/m2/s]
               IF ( IDTLIMO > 0 ) THEN
                  FLUXLIMO(I,J) = EMIS_LIMO * FC1 
               ENDIF

               !--------------------------------------------------------------
               ! MEGAN MTPO
               !--------------------------------------------------------------
               ! MTPO is all other monoterpenes (MEGAN categories:
               ! myrcene, ocimene, OMON) (hotp 5/20/10)
               ! All other monoterpenes (mostly camphene, linalool,
               ! terpinolene, terpinolene, phellandrene) (hotp 3/10/10)
               ! 14-18% of OMTP is terpinene and terpinolene
               IF ( IDTMTPO > 0 ) THEN
                  FLUXMTPO(I,J) = ( EMIS_MYRC + EMIS_OCIM +
     &                              EMIS_OMON ) * FC1
               ENDIF

               !--------------------------------------------------------------
               ! MEGAN sesquiterpenes
               !--------------------------------------------------------------

               ! ---------------------------------------------------
               ! a-Farnesene 
               CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                   'FARN',   EMIS_FARN, RC )
               IF ( RC /= HCO_SUCCESS ) THEN 
                  CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_FARN', RC )
                  ERR = .TRUE.
                  EXIT
               ENDIF 
               FLUXFARN(I,J) = EMIS_FARN

               ! ---------------------------------------------------
               ! b_Caryophyllene 
               CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                   'BCAR',   EMIS_BCAR, RC )
               IF ( RC /= HCO_SUCCESS ) THEN 
                  CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_BCAR', RC )
                  ERR = .TRUE.
                  EXIT
               ENDIF 
               FLUXBCAR(I,J) = EMIS_BCAR

               ! ---------------------------------------------------
               ! Other sesquiterpene
               CALL GET_MEGAN_EMISSIONS( HcoState, ExtState,  I, J,
     &                                   'OSQT',   EMIS_OSQT, RC )
               IF ( RC /= HCO_SUCCESS ) THEN 
                  CALL HCO_ERROR( 'GET_MEGAN_EMISSIONS_OSQT', RC )
                  ERR = .TRUE.
                  EXIT
               ENDIF 
               FLUXOSQT(I,J) = EMIS_OSQT

               ! ---------------------------------------------------
               ! Total sesquiterpenes from MEGAN (hotp 3/10/10)
               IF ( IDTSESQ > 0 ) THEN
                  FLUXSESQ(I,J) = ( EMIS_FARN + EMIS_BCAR +
     &                              EMIS_OSQT ) * FC3
               ENDIF

            ! If MEGAN Monoterpene is turned off:
            ! NOTE: In the original GEOS-Chem code, many of these fluxes
            ! were derived from the monoterpenese fluxes. These fluxes
            ! are only calculated if the MEGAN_MONO option is turned on
            ! and their contributions are hence skipped in the following!
            ! (ckeller, 05/19/14)
            ELSE

               ! MTPA
               ! old ALPH is now MTPA (hotp 5/20/10)
               ! ALPHA-PINENE (0.35)
               ! BETA-PINENE lumped with ALPHA-PINENE (0.23)
               ! SABINENE    lumped with ALPHA-PINENE (0.05)
               ! D3-CARENE   lumped with ALPHA-PINENE (0.04)

               ! TERPENOID KETONE is lumped with SABINENE
               ! Then SABINENE    is lumped with ALPHA-PINENE
               IF ( IDTMTPA > 0 ) THEN
                  FLUXMTPA(I,J) = DIURORVC * FC4 * 0.04_hp ! using campher
               ENDIF

               ! LIMONENE not defined (derived from TERP)
               IF ( IDTLIMO > 0 ) THEN
                  FLUXLIMO(I,J) = 0.0_hp
               ENDIF  

               ! MTPO
               IF ( IDTMTPO > 0 ) THEN
                  FLUXMTPO(I,J) = DIURORVC * FC2 * 0.09_hp ! using LINALOOL 
               ENDIF

               ! IDTSESQ 
               IF ( IDTSESQ > 0 ) THEN
                  FLUXSESQ(I,J) = DIURORVC * FC3 * 0.05_hp
               ENDIF

            ENDIF
         ENDIF

      ENDDO !I
      ENDDO !J
!$OMP END PARALLEL DO

      IF ( ERR ) THEN
         RC = HCO_FAIL
         RETURN
      ENDIF

      !=================================================================
      ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
      !=================================================================

      ! ----------------------------------------------------------------
      ! ISOPRENE
      IF ( IDTISOP > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXISOP, IDTISOP, RC,
     &                     ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXISOP', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MONOTERPENES 
      IF ( ExtNrMono > 0 ) THEN

         IF ( IDTMONX > 0 ) THEN
            ! Add flux to emission array
            CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXMONO, IDTMONX, 
     &                        RC,        ExtNr=ExtNr )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXMONO', RC )
               RETURN 
            ENDIF
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ACETALDEHYDE
      IF ( IDTALD2 > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXALD2, IDTALD2, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXALD2', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ACETONE 
      IF ( IDTACET > 0 ) THEN

         ! Eventually add individual diagnostics. These are assumed to
         ! have names MEGAN_ACET_MONO, MEGAN_ACET_MBO, and 
         ! MEGAN_ACET_DIRECT
         DiagnName =  'MEGAN_ACET_MONO'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &      cName=TRIM(DiagnName), Array2D=FLUXACETmo, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         DiagnName =  'MEGAN_ACET_MBO'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &      cName=TRIM(DiagnName), Array2D=FLUXACETmb, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         DiagnName =  'MEGAN_ACET_DIRECT'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &      cName=TRIM(DiagnName), Array2D=FLUXACETbg, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Add flux to emission array
         FLUXACETbg = FLUXACETbg + FLUXACETmo + FLUXACETmb
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXACETbg, IDTACET, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXACETbg', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! CO
      IF ( ( ExtNrMono > 0 ) .AND. ( IDTCO > 0 ) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXCO, IDTCO, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXCO', RC )
            RETURN 
         ENDIF

         ! Also archive the BIOGENIC_CO diagnostic (bmy, 4/29/15)
         DiagnName =  'BIOGENIC_CO'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono,
     &      cName=TRIM(DiagnName), Array2D=FLUXCO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN
         Arr2D => NULL()
      ENDIF

      ! ----------------------------------------------------------------
      ! ALKENES 
      IF ( IDTPRPE > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXPRPE, IDTPRPE, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXPRPE', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ETHENE 
      IF ( IDTC2H4 > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXC2H4, IDTC2H4, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXC2H4', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! OC
      ! Note: OC emissions will only be non-zero if OC species is defined
      ! and Monoterpene option is enabled.
      IF ( IDTOCPI > 0 ) THEN 

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXOCPI, IDTOCPI, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXOCPI', RC )
            RETURN 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MTPA 
      IF ( ( ExtNrSoa > 0 ) .AND. ( IDTMTPA > 0 ) ) THEN 

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXMTPA, IDTMTPA, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXMTPA', RC )
            RETURN 
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! MTPO 
      IF ( ( ExtNrSoa > 0 ) .AND. ( IDTMTPO > 0 ) ) THEN 

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXMTPO, IDTMTPO, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXMTPO', RC )
            RETURN 
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! LIMONENE 
      IF ( ( ExtNrSoa > 0 ) .AND. ( IDTLIMO > 0 ) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXLIMO, IDTLIMO, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXLIMO', RC )
            RETURN 
         ENDIF

      ENDIF

      ! ----------------------------------------------------------------
      ! SESQ 
      IF ( (ExtNrSoa>0) .AND. (IDTSESQ>0) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXSESQ, IDTSESQ, 
     &                     RC,        ExtNr=ExtNr )
         IF ( RC /= HCO_SUCCESS ) THEN
            CALL HCO_ERROR( 'HCO_EmisAdd error: FLUXSESQ', RC )
            RETURN 
         ENDIF

      ENDIF

      ! Restore MEGAN monoterpenes and sesquiterpenes to diagnostics
      ! Add individual species to diagnostics rather than lumped species,
      ! following method used in pre-HEMCO diagnostics (mps, 2/18/15)
      IF ( ExtNrMono > 0 ) THEN

         ! -------------------------------------------------------------
         ! Methyl butenol
         DiagnName =  'BIOGENIC_MBOX'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &     cName=TRIM(DiagnName), Array2D=FLUXMBOX, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Alpha Pinene
         DiagnName =  'BIOGENIC_APIN'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXAPIN, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Beta Pinene
         DiagnName =  'BIOGENIC_BPIN'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXBPIN, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Limonene
         DiagnName =  'BIOGENIC_LIMO'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXLIMO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Sabinene
         DiagnName =  'BIOGENIC_SABI'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXSABI, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Mycrene
         DiagnName =  'BIOGENIC_MYRC'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXMYRC, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! 3-Carene
         DiagnName =  'BIOGENIC_CARE'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXCARE, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Ocimene
         DiagnName =  'BIOGENIC_OCIM'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXOCIM, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Other monoterpenes
         DiagnName =  'BIOGENIC_OMON'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXOMON, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Total monoterpenes
         DiagnName =  'BIOGENIC_MONX'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXMONO, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! a-Farnesene
         DiagnName =  'BIOGENIC_FARN'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXFARN, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! b_Caryophyllene
         DiagnName =  'BIOGENIC_BCAR'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXBCAR, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! -------------------------------------------------------------
         ! Other sesquiterpenes
         DiagnName =  'BIOGENIC_OSQT'
         CALL Diagn_Update( am_I_Root, ExtNr=ExtNrMono, 
     &      cName=TRIM(DiagnName), Array2D=FLUXOSQT, RC=RC )
         IF ( RC /= HCO_SUCCESS ) RETURN 

      ENDIF

      ! ----------------------------------------------------------------
      ! Other MEGAN biogenics
      !
      ! Calls included here for future incorporation or
      ! specialized simulations. (dbm, 12/2012)
      !
      ! These diagnostics will be zero because FLUXMOHX,
      !  FLUXETOH, FLUXFAXX, and FLUXAAXX are commented out
      !  above (mps, 2/19/15)
      ! ----------------------------------------------------------------
      ! Methanol
      DiagnName =  'BIOGENIC_MOHX'
      CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &   cName=TRIM(DiagnName), Array2D=FLUXMOHX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! ----------------------------------------------------------------
      ! Ethanol
      DiagnName =  'BIOGENIC_ETOH'
      CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &   cName=TRIM(DiagnName), Array2D=FLUXETOH, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! ----------------------------------------------------------------
      ! Formic acid
      DiagnName =  'BIOGENIC_FAXX'
      CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &   cName=TRIM(DiagnName), Array2D=FLUXFAXX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! ----------------------------------------------------------------
      ! Acetic acid
      DiagnName =  'BIOGENIC_AAXX'
      CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, 
     &   cName=TRIM(DiagnName), Array2D=FLUXAAXX, RC=RC)
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! ----------------------------------------------------------------
      ! Eventually copy internal values to ESMF internal state object
      ! ----------------------------------------------------------------

      ! LAI_PREVDAY
      CALL HCO_RestartWrite( am_I_Root, HcoState, 
     &                       'LAI_PREVDAY', LAI_PREVDAY, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! T_DAVG
      CALL HCO_RestartWrite( am_I_Root, HcoState, 
     &                       'T_DAVG', HCOT_15_AVG, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! T_PREVDAY
      CALL HCO_RestartWrite( am_I_Root, HcoState, 
     &                       'T_PREVDAY', HCOT_15(:,:,1), RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! PARDR_DAVG
      CALL HCO_RestartWrite( am_I_Root, HcoState, 
     &                       'PARDR_DAVG', HCOPARDR_15_AVG, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! PARDF_DAVG
      CALL HCO_RestartWrite( am_I_Root, HcoState, 
     &                       'PARDF_DAVG', HCOPARDF_15_AVG, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! ALL DONE! 
      !=================================================================

!      ! testing only
!      write(*,*) ''
!      write(*,*) 'MEGAN done!'
!      write(*,*) 'total LAI : ', SUM(ExtState%GC_LAI%Arr%Val)
!      write(*,*) 'total ISOP: ', SUM(FLUXISOP)
!      write(*,*) 'total ACET: ', SUM(FLUXACETmo)+SUM(FLUXACETmb)
!     &                          +SUM(FLUXACETbg)
!      write(*,*) 'total CO  : ', SUM(FLUXCO  )
!      write(*,*) 'total PRPE: ', SUM(FLUXPRPE)
!      write(*,*) 'total OCPI: ', SUM(FLUXOCPI)
!      write(*,*) ''

      ! Keep LAI in memory (ckeller, 10/9/2014)
      LAI_PREVSTEP = ExtState%LAI%Arr%Val

      ! Return w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE HCOX_Megan_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_Emissions
!
! !DESCRIPTION: Subroutine Get\_Megan\Emissions computes biogenic emissions in 
!  units of [kgC/m2/s] or [kg/m2/s] using the MEGAN inventory. (dbm, 12/2012)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MEGAN_EMISSIONS( HcoState, ExtState, I, J, 
     &                                CMPD, MEGAN_EMIS, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE),  POINTER     :: HcoState
      TYPE(Ext_State),  POINTER     :: ExtState
      INTEGER,          INTENT(IN)  :: I, J      ! GEOS-Chem lon & lat indices
      CHARACTER(LEN=*), INTENT(IN)  :: CMPD      ! Compound name (dbm,6/21/2012)
!
! !OUTPUT PARAMETERS:
!
      REAL(hp),         INTENT(OUT) :: MEGAN_EMIS ! VOC emission in kgC/m2/s 
                                                  ! or kg/m2/s, depending on
                                                  ! units the compound is
                                                  ! carried in
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT) :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
!  (4 ) Guenther et al, 2012 GMD MEGANv2.1 description and associated code at
!                                http://acd.ucar.edu/~guenther/MEGAN/
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Updated to the latest 
!        algorithm and modified for the standard code by May Fu (11/20/04)
!  (2 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (3 ) Restructing of function & implementation of activity factors (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  11 Apr 2012 - R. Yantosca - Now use data from modis_lai_mod.F90
!  11 Apr 2012 - R. Yantosca - Cosmetic changes
!  26 Jan 2015 - M. Sulprizio- Update from D. Millet (21 Jun 2012): New driver
!                              routine for all MEGAN compounds
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)  :: GAMMA_LAI
      REAL(hp)  :: GAMMA_AGE
      REAL(hp)  :: GAMMA_PAR
      REAL(hp)  :: GAMMA_T_LD
      REAL(hp)  :: GAMMA_T_LI
      REAL(hp)  :: GAMMA_SM
      REAL(hp)  :: AEF
      REAL(hp)  :: D_BTW_M
      REAL(hp)  :: TS, SUNCOS
      REAL(hp)  :: Q_DIR_2, Q_DIFF_2
      REAL(hp)  :: BETA, LDF, CT1, CEO
      REAL(hp)  :: ANEW, AGRO, AMAT, AOLD
      REAL(hp)  :: ISOLAI, PMISOLAI, MISOLAI
      LOGICAL   :: BIDIR    

      !=================================================================
      ! GET_MEGAN_EMISSIONS begins here!
      !================================================================= 

      ! Initialize parameters, gamma values, and return value
      MEGAN_EMIS = 0.0_hp
      GAMMA_LAI  = 0.0_hp
      GAMMA_AGE  = 0.0_hp
      GAMMA_T_LD = 0.0_hp
      GAMMA_T_LI = 0.0_hp
      GAMMA_PAR  = 0.0_hp
      GAMMA_SM   = 0.0_hp
      BETA       = 0.0_hp
      AEF        = 0.0_hp
      LDF        = 0.0_hp
      CT1        = 0.0_hp
      CEO        = 0.0_hp
      ANEW       = 0.0_hp
      AGRO       = 0.0_hp
      AMAT       = 0.0_hp
      AOLD       = 0.0_hp
      BIDIR      = .FALSE.  

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Pass met variables. Now use only LAI (ckeller, 10/9/2014)
      ISOLAI   = ExtState%LAI%Arr%Val(I,J)
      PMISOLAI = LAI_PREVDAY(I,J)
      MISOLAI  = ISOLAI
      TS       = ExtState%T2M%Arr%Val(I,J)
      SUNCOS   = ExtState%SUNCOS%Arr%Val(I,J)

      ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
      Q_DIR_2  = ExtState%PARDR%Arr%Val(I,J) * WM2_TO_UMOLM2S
      Q_DIFF_2 = ExtState%PARDF%Arr%Val(I,J) * WM2_TO_UMOLM2S

      ! --------------------------------------------
      ! Get MEGAN parameters for this compound
      ! --------------------------------------------
      CALL GET_MEGAN_PARAMS ( CMPD, BETA, LDF,  CT1,  CEO, 
     &                        ANEW, AGRO, AMAT, AOLD, BIDIR, RC )

      ! --------------------------------------------
      ! Get base emission factor for this compound and grid square
      ! Units: kgC/m2/s or kg/m2/s
      ! --------------------------------------------
      CALL GET_MEGAN_AEF ( I, J, CMPD, AEF, RC )
     
      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( ISOLAI * AEF > 0.0_hp ) THEN

         ! --------------------------------------------------
         ! GAMMA_par (light activity factor)
         ! --------------------------------------------------

         ! Calculate GAMMA PAR only during day
         IF ( SUNCOS > 0.0_hp ) THEN

            GAMMA_PAR = GET_GAMMA_PAR_PCEEA(HcoState, ExtState, I, J, 
     &                                      Q_DIR_2,  Q_DIFF_2,
     &                                      HCOPARDR_15_AVG(I,J),
     &                                      HCOPARDF_15_AVG(I,J)  )

         ELSE             
 
            ! If night
            GAMMA_PAR = 0.0_hp  
            
         ENDIF 

         ! --------------------------------------------------
         ! GAMMA_T_LI (temperature activity factor for
         ! light-independent fraction)
         ! --------------------------------------------------
         GAMMA_T_LI = GET_GAMMA_T_LI( TS, BETA )

         ! --------------------------------------------------
         ! GAMMA_T_LD (temperature activity factor for
         ! light-dependent fraction)
         ! --------------------------------------------------
         GAMMA_T_LD = GET_GAMMA_T_LD( TS, HCOT_15_AVG(I,J), 
     &        HCOT_15(I,J,1), CT1, CEO )

         ! --------------------------------------------------
         ! GAMMA_LAI (leaf area index activity factor)
         ! --------------------------------------------------
         GAMMA_LAI = GET_GAMMA_LAI( MISOLAI, BIDIR )

         ! --------------------------------------------------
         ! GAMMA_AGE (leaf age activity factor)
         ! --------------------------------------------------
         GAMMA_AGE = GET_GAMMA_AGE( MISOLAI, PMISOLAI, 
     &        D_BTW_M, HCOT_15_AVG(I,J), ANEW, AGRO, AMAT, AOLD )

         ! --------------------------------------------------
         ! GAMMA_SM (soil moisture activity factor)
         ! --------------------------------------------------
         GAMMA_SM = GET_GAMMA_SM( ExtState, I, J, CMPD )

      ELSE

         ! set activity factors to zero
         GAMMA_PAR  = 0.0_hp
         GAMMA_T_LI = 0.0_hp
         GAMMA_T_LD = 0.0_hp
         GAMMA_LAI  = 0.0_hp
         GAMMA_AGE  = 0.0_hp
         GAMMA_SM   = 0.0_hp

      END IF

      ! Emission is the product of all of these.
      ! Units here are kgC/m2/s or kg/m2/s as appropriate for the compound.
      ! Normalization factor ensures product of GAMMA values is 1.0 under
      !  standard conditions. 
      MEGAN_EMIS = NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * GAMMA_LAI
     &           * ((1.0_hp - LDF) * GAMMA_T_LI +
     &             (LDF * GAMMA_PAR * GAMMA_T_LD))


!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' ' 
!         write(*,*) '--- GET_MEGAN_EMISSIONS --- '
!         write(*,*) 'HEMCO MEGAN @    ', i,j
!         write(*,*) 'Compound       : ', TRIM(CMPD)
!         write(*,*) 'MEGAN_EMIS     : ', MEGAN_EMIS
!         write(*,*) 'SUNCOS         : ', SUNCOS
!         write(*,*) 'AEF [kgC/m2/s] : ', AEF
!         write(*,*) 'GAMMA_LAI      : ', GAMMA_LAI
!         write(*,*) 'GAMMA_AGE      : ', GAMMA_AGE
!         write(*,*) 'GAMMA_T_LI     : ', GAMMA_T_LI
!         write(*,*) 'GAMMA_T_LD     : ', GAMMA_T_LD
!         write(*,*) 'GAMMA_PAR      : ', GAMMA_PAR
!         write(*,*) 'GAMMA_SM       : ', GAMMA_SM
!         write(*,*) 'TS             : ', TS
!         write(*,*) 'HCOT_15_AVG    : ', HCOT_15_AVG(I,J)
!         write(*,*) 'HCOT_DAILY     : ', HCOT_DAILY(I,J)
!         write(*,*) 'HCOPARDR_15_AVG: ', HCOPARDR_15_AVG(I,J)
!         write(*,*) 'HCOPARDF_15_AVG: ', HCOPARDF_15_AVG(I,J)
!         write(*,*) 'ISOLAI         : ', ISOLAI
!         write(*,*) 'MISOLAI        : ', MISOLAI
!         write(*,*) 'PMISOLAI       : ', PMISOLAI
!         write(*,*) 'D_BTW_M        : ', D_BTW_M
!         write(*,*) ' ' 
!      endif

      ! Leave w/ success
      RC = HCO_SUCCESS
 
      END SUBROUTINE GET_MEGAN_EMISSIONS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_Params
!
! !DESCRIPTION: Subroutine Get\_Megan\_Params returns the emission parameters
!  for each MEGAN compound needed to compute emissions. Called from
!  GET\_MEGAN\_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MEGAN_PARAMS( CPD,   BTA,   LIDF,  C_T1,  C_EO,
     &                             A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR,
     &                             RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: CPD    ! Compound name
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(hp), INTENT(INOUT) :: BTA    ! Beta coefficient for temperature activity 
                                        ! factor for light-independent fraction
      REAL(hp), INTENT(INOUT) :: LIDF   ! Light-dependent fraction of emissions
      REAL(hp), INTENT(INOUT) :: C_T1   ! CT1 parameter for temperature activity
                                        ! factor for light-dependent fraction
      REAL(hp), INTENT(INOUT) :: C_EO   ! Ceo parameter for temperature activity
                                        ! factor for light-dependent fraction
      REAL(hp), INTENT(INOUT) :: A_NEW  ! Relative emission factor (new leaves)
      REAL(hp), INTENT(INOUT) :: A_GRO  ! Relative emission factor (growing leaves)
      REAL(hp), INTENT(INOUT) :: A_MAT  ! Relative emission factor (mature leaves)
      REAL(hp), INTENT(INOUT) :: A_OLD  ! Relative emission factor (old leaves)
      LOGICAL,  INTENT(INOUT) :: BI_DIR ! Logical flag to indicate bidirectional exchange
      INTEGER,  INTENT(INOUT) :: RC
! 
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code
!
! !REVISION HISTORY: 
!  (1 ) Created by dbm 07/2012
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255):: MSG

      !=================================================================
      ! GET_MEGAN_PARAMS begins here!
      !================================================================= 

      ! Initialize values
      BTA    = 0.0_hp           
      LIDF   = 0.0_hp    
      C_T1   = 0.0_hp    
      C_EO   = 0.0_hp     
      A_NEW  = 0.0_hp  
      A_GRO  = 0.0_hp   
      A_MAT  = 0.0_hp    
      A_OLD  = 0.0_hp   
      BI_DIR = .FALSE. 

      ! ----------------------------------------------------------------
      ! Note that not all the above compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or 
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding CPD name and the appropriate paramaters.
      ! (dbm, 01/2013)
      !
      ! Values are from Table 4 in Guenther et al., 2012
      ! ----------------------------------------------------------------

      ! Isoprene, MBO
      IF ( TRIM(CPD) == 'ISOP' .OR.
     &     TRIM(CPD) == 'MBOX' ) THEN
         BTA    = 0.13_hp  ! Not actually used for ISOP, MBO
         LIDF   = 1.0_hp
         C_T1   = 95.0_hp
         C_EO   = 2.0_hp
         A_NEW  = 0.05_hp
         A_GRO  = 0.6_hp
         A_MAT  = 1.0_hp
         A_OLD  = 0.9_hp
         BI_DIR = .FALSE. 

      ! Myrcene, sabinene, alpha-pinene
      ELSE IF ( TRIM(CPD) == 'MYRC' .OR. 
     &          TRIM(CPD) == 'SABI' .OR.
     &          TRIM(CPD) == 'APIN' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.6_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE. 

      ! Limonene, 3-carene, beta-pinene
      ELSE IF ( TRIM(CPD) == 'LIMO' .OR. 
     &          TRIM(CPD) == 'CARE' .OR.
     &          TRIM(CPD) == 'BPIN' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE. 
       
      ! t-beta-ocimene
      ELSE IF ( TRIM(CPD) == 'OCIM' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.8_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE. 

      ! Other monoterpenes (lumped)
      ELSE IF ( TRIM(CPD) == 'OMON' ) THEN
         BTA    = 0.10_hp
         LIDF   = 0.4_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 2.0_hp
         A_GRO  = 1.8_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.05_hp
         BI_DIR = .FALSE. 
         
      ! Methanol
      ELSE IF ( TRIM(CPD) == 'MOHX' ) THEN
         BTA    = 0.08_hp
         LIDF   = 0.8_hp
         C_T1   = 60.0_hp
         C_EO   = 1.6_hp
         A_NEW  = 3.5_hp
         A_GRO  = 3.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.2_hp
         BI_DIR = .FALSE. 

      ! Acetone
      ELSE IF ( TRIM(CPD) == 'ACET' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE. 
         
      ! Bidirectional VOC: Ethanol, formaldehyde, acetaldehyde, formic acid,
      ! acetic acid
      ELSE IF ( TRIM(CPD) == 'ETOH' .OR.
     &          TRIM(CPD) == 'CH2O' .OR. 
     &          TRIM(CPD) == 'ALD2' .OR.
     &          TRIM(CPD) == 'FAXX' .OR.
     &          TRIM(CPD) == 'AAXX' ) THEN
         BTA    = 0.13_hp
         LIDF   = 0.8_hp
         C_T1   = 95.0_hp
         C_EO   = 2.0_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .TRUE.  
        
      ! Stress VOCs: ethene, toluene, HCN
      ! There are others species in this category but none are currently 
      ! used in GEOS-Chem
      ELSE IF ( TRIM(CPD) == 'C2H4' .OR.
     &          TRIM(CPD) == 'TOLU' .OR. 
     &          TRIM(CPD) == 'HCNX' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.8_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE. 
         
      ! Other VOCs: >C2 alkenes
      ! This includes propene, butene and very minor contribution from
      ! larger alkenes
      ELSE IF ( TRIM(CPD) == 'PRPE' ) THEN
         BTA    = 0.1_hp
         LIDF   = 0.2_hp
         C_T1   = 80.0_hp
         C_EO   = 1.83_hp
         A_NEW  = 1.0_hp
         A_GRO  = 1.0_hp
         A_MAT  = 1.0_hp
         A_OLD  = 1.0_hp
         BI_DIR = .FALSE. 
         
      ! SOAupdate: Sesquiterpenes hotp 3/2/10
      ! alpha-Farnesene, beta-Caryophyllene, other sesquiterpenes
      ELSE IF ( TRIM(CPD) == 'FARN' .OR.
     &          TRIM(CPD) == 'BCAR' .OR.
     &          TRIM(CPD) == 'OSQT' ) THEN
         BTA    = 0.17_hp
         LIDF   = 0.5_hp
         C_T1   = 130.0_hp
         C_EO   = 2.37_hp
         A_NEW  = 0.4_hp
         A_GRO  = 0.6_hp
         A_MAT  = 1.0_hp
         A_OLD  = 0.95_hp
         BI_DIR = .FALSE. 

      ! Calls for any other MEGAN compounds (e.g. sesquiterpenes, CO, etc.) can
      ! be added following the above format based on the parameters in
      ! Guenther 2012 or the MEGAN source code (dbm, 6/21/2012).
      ELSE

         MSG = 'Invalid compound name'
         CALL HCO_ERROR( MSG, RC, THISLOC='GET_MEGAN_PARAMS' )
         RETURN 

      ENDIF
     
      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_MEGAN_PARAMS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_AEF
!
! !DESCRIPTION: Function Get\_Megan\_AEF returns the appropriate AEF value
!  for a given compound and grid square.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_MEGAN_AEF( I, J, CPD, EMFAC, RC )
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN)  :: I, J        ! Lon & lat indices
      CHARACTER(LEN=*),  INTENT(IN)  :: CPD         ! Compound name
!
! !OUTPUT PARAMETERS: 
!
      REAL(hp),          INTENT(OUT) :: EMFAC       ! MEGAN base emission factor
                                                    ! (kgC/m2/s or kg/m2/s)
                                                    ! for grid cell (I,J)
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,           INTENT(INOUT) :: RC        ! Return code
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2012, MEGANv2.1 source code
! 
! !REVISION HISTORY: 
!  (1 ) Created 11/2012 by dbm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255):: MSG

      !=================================================================
      ! GET_MEGAN_AEF begins here!
      !================================================================= 

      ! Find appropriate tracer
      SELECT CASE ( TRIM( CPD ) )
         CASE ( 'ISOP' )
            EMFAC = AEF_ISOP(I,J)
         CASE ( 'MBOX' )
            EMFAC = AEF_MBOX(I,J)
         CASE ( 'MYRC' )
            EMFAC = AEF_MYRC(I,J)
         CASE ( 'SABI' )
            EMFAC = AEF_SABI(I,J)
         CASE ( 'APIN' )
            EMFAC = AEF_APIN(I,J)
         CASE ( 'LIMO' )
            EMFAC = AEF_LIMO(I,J)
         CASE ( 'CARE' )
            EMFAC = AEF_CARE(I,J)
         CASE ( 'BPIN' )
            EMFAC = AEF_BPIN(I,J)
         CASE ( 'OCIM' )
            EMFAC = AEF_OCIM(I,J)
         CASE ( 'OMON' )
            EMFAC = AEF_OMON(I,J)
         CASE ( 'MOHX' )
            EMFAC = AEF_MOHX(I,J)
         CASE ( 'ACET' )
            EMFAC = AEF_ACET(I,J)
         CASE ( 'ETOH' )
            EMFAC = AEF_ETOH(I,J)
         CASE ( 'CH2O' )
            EMFAC = AEF_CH2O(I,J)
         CASE ( 'ALD2' )
            EMFAC = AEF_ALD2(I,J)
         CASE ( 'FAXX' )
            EMFAC = AEF_FAXX(I,J)
         CASE ( 'AAXX' )
            EMFAC = AEF_AAXX(I,J)
         CASE ( 'C2H4' )
            EMFAC = AEF_C2H4(I,J)
         CASE ( 'TOLU' )
            EMFAC = AEF_TOLU(I,J)
         CASE ( 'HCNX' )
            EMFAC = AEF_HCNX(I,J)
         CASE ( 'PRPE' )
            EMFAC = AEF_PRPE(I,J)
         CASE ( 'FARN' )
            EMFAC = AEF_FARN(I,J)
         CASE ( 'BCAR' )
            EMFAC = AEF_BCAR(I,J)
         CASE ( 'OSQT' )
            EMFAC = AEF_OSQT(I,J)
         CASE DEFAULT     
            MSG = 'Invalid compound name'
            CALL HCO_ERROR( MSG, RC, THISLOC='GET_MEGAN_AEF' )
            RETURN 
      END SELECT

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_MEGAN_AEF
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_PAR_PCEEA
!
! !DESCRIPTION: Computes the PCEEA gamma activity factor with sensitivity 
!  to LIGHT.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_PAR_PCEEA( HcoState, ExtState, I, J, Q_DIR_2, 
     &                              Q_DIFF_2, HCOPARDR_AVG_SIM, 
     &                              HCOPARDF_AVG_SIM )
     &              RESULT( GAMMA_P_PCEEA )
!
! !USES:
!
      USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal 
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER    :: HcoState
      TYPE(Ext_State), POINTER    :: ExtState
      INTEGER,         INTENT(IN) :: I,  J             ! Lon & lat indices 
      REAL(sp),        INTENT(IN) :: HCOPARDR_AVG_SIM  ! Avg direct PAR [W/m2]
      REAL(sp),        INTENT(IN) :: HCOPARDF_AVG_SIM  ! Avg diffuse PAR [W/m2]
      REAL(hp),        INTENT(IN) :: Q_DIR_2           ! Direct PAR [umol/m2/s]
      REAL(hp),        INTENT(IN) :: Q_DIFF_2          ! Diffuse PAR [umol/m2/s]
!
! !RETURN VALUE:
!
      REAL(hp)                    :: GAMMA_P_PCEEA     ! GAMMA factor for light
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 user guide
! 
! !REVISION HISTORY: 
!  (1 ) Here PAR*_AVG_SIM is the average light conditions over the simulation 
!       period. I've set this = 10 days to be consistent with temperature & as 
!       outlined in Guenther et al, 2006. (mpb,2009)
!  (2 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!       (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)   :: mmHCOPARDR_DAILY
      REAL(hp)   :: mmHCOPARDF_DAILY
      REAL(hp)   :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL(hp)   :: PTOA, PHI
      REAL(hp)   :: BETA,   SINbeta 
      INTEGER    :: DOY, RC 
      REAL(hp)   :: AAA, BBB
      REAL(hp)   :: LocalHour
      REAL(hp)   :: LAT

      !-----------------------------------------------------------------
      ! Compute GAMMA_PAR_PCEEA
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD   = 0.0_hp
      PTOA     = 0.0_hp

      ! Convert past light conditions to micromol/m2/s 
      mmHCOPARDR_DAILY   = HCOPARDR_AVG_SIM  * WM2_TO_UMOLM2S
      mmHCOPARDF_DAILY   = HCOPARDF_AVG_SIM  * WM2_TO_UMOLM2S

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmHCOPARDR_DAILY + mmHCOPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      ! Get latitude
      LAT = HcoState%Grid%YMID%Val(I,J)

      ! Get day of year, local-time and latitude
      ! TODO: Evaluate RC?
      CALL HcoClock_Get ( cDOY = DOY, RC=RC )
      CALL HcoClock_GetLocal ( HcoState, I, J, cH = LocalHour, RC=RC )

      ! Get solar elevation angle
      SINbeta      =  SOLAR_ANGLE( HcoState, DOY, LocalHour, LAT )
      BETA         =  ASIN( SINbeta ) * RAD2D       

      IF ( SINbeta < 0.0_hp ) THEN

         GAMMA_P_PCEEA = 0.0_hp

      ELSEIF ( SINbeta > 0.0_hp ) THEN       

         ! PPFD at top of atmosphere
         PTOA    = 3000.0_hp + 99.0_hp * 
     &             COS( 2._hp * 3.14_hp *( DOY - 10.0_hp ) / 365.0_hp )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0_hp + 0.0005_hp *( PAC_DAILY - 400.0_hp ) 
         AAA     = ( 2.46_hp * BBB * PHI ) - ( 0.9_hp * PHI**2 )

         GAMMA_P_PCEEA = SINbeta * AAA

      ENDIF

       ! Screen unforced errors. IF solar elevation angle is 
       ! less than 1 THEN gamma_p can not be greater than 0.1.
       IF ( BETA < 1.0_hp .AND. GAMMA_P_PCEEA > 0.1_hp ) THEN
          GAMMA_P_PCEEA  = 0.0_hp
       ENDIF
   
      ! Prevent negative values
      GAMMA_P_PCEEA = MAX( GAMMA_P_PCEEA , 0.0_hp )

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) 'HEMCO GAMMA_PAR_PCEEA: ', GAMMA_P_PCEEA
!         write(*,*) 'DOY                  : ', DOY
!         write(*,*) 'LAT                  : ', LAT
!         write(*,*) 'SINbeta              : ', SINbeta
!         write(*,*) 'BETA                 : ', BETA
!         write(*,*) 'PTOA                 : ', PTOA
!         write(*,*) 'PHI                  : ', PHI
!         write(*,*) 'BBB                  : ', BBB
!         write(*,*) 'AAA                  : ', AAA
!         write(*,*) 'PAC_DAILY            : ', PAC_DAILY
!         write(*,*) 'PAC_INSTANT          : ', PAC_INSTANT
!      endif

      END FUNCTION GET_GAMMA_PAR_PCEEA
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Solar_Angle
!
! !DESCRIPTION: Function SOLAR\_ANGLE computes the local solar angle for a 
!  given day of year, latitude and longitude (or local time).  Called from  
!  routine Get\_Gamma\_P\_Pecca.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SOLAR_ANGLE( HcoState, DOY, SHOUR, LAT ) RESULT(SINbeta)
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      TYPE(HCO_State),   POINTER    :: HcoState
      INTEGER,           INTENT(IN) :: DOY       ! Day of year
      REAL(hp),          INTENT(IN) :: SHOUR     ! Local time 
      REAL(hp),          INTENT(IN) :: LAT       ! Latitude
!
! !RETURN VALUE:
!
      REAL(hp)                      :: SINbeta   ! Sin of the local solar angle
!
! !REMARKS:
!  References (see above for full citations):
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN v2.1 user mannual 2007-09
! 
! !REVISION HISTORY: 
!  (1 ) This code was taken directly from the MEGAN v2.1 source code.(mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp) :: BETA                        ! solar elevation angle
      REAL(hp) :: sindelta, cosdelta, A, B

      ! Calculation of sin beta 
      sindelta = -SIN( 0.40907_hp ) * 
     &            COS( 6.28_hp * ( DOY + 10_dp ) / 365_dp )

      cosdelta = (1-sindelta**2.0_hp)**0.5_hp

      A = SIN( LAT * D2RAD ) * sindelta
      B = COS( LAT * D2RAD ) * cosdelta

      SINbeta = A + B * 
     &          COS( 2.0_hp * HcoState%Phys%PI * ( SHOUR-12_dp )/24_dp )

      END FUNCTION SOLAR_ANGLE 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_LI
!
! !DESCRIPTION: Function Get\_Gamma\_T\_LI computes the temperature activity
!  factor (GAMMA\_T\_LI) for the light-independent fraction of emissions
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LI( T, BETA ) RESULT( GAMMA_T_LI )
!
! !INPUT PARAMETERS: 
!
      ! Current leaf temperature, the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL(hp),  INTENT(IN) :: T 

      ! Temperature factor per species
      REAL(hp),  INTENT(IN) :: BETA
!
! !RETURN VALUE:
!
      ! Activity factor for the light-independent fraction of emissions
      REAL(hp)              :: GAMMA_T_LI
!
! !REMARKS:
!  GAMMA_T =  exp[Beta*(T - T_Standard)]
!                                                                             .
!             where Beta   = temperature dependent parameter
!                   Ts     = standard temperature (normally 303K, 30C)
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Michael Barkley (2009).
!       Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to GAMMA_T_LI (dbm, 6/21/2012)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: T_STANDARD = 303.d0

      !=================================================================
      ! GET_GAMMAT_T_LI begins here!
      !================================================================= 

      GAMMA_T_LI = EXP( BETA * ( T - T_STANDARD ) )

      END FUNCTION GET_GAMMA_T_LI
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_LD
!
! !DESCRIPTION: Function Get\_Gamma\_T\_LD computes the temperature 
!  sensitivity for the light-dependent fraction of emissions.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_LD( T, PT_15, PT_1, CT1, CEO ) 
     &     RESULT( GAMMA_T_LD )
!
! !INPUT PARAMETERS: 
!
      ! Current leaf temperature [K], the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL(hp), INTENT(IN) :: T        

       ! Average leaf temperature over the past 15 days
      REAL(sp), INTENT(IN) :: PT_15

      ! Average leaf temperature over the past arbitray day(s).
      ! This is not used at present
      REAL(sp), INTENT(IN) :: PT_1

      ! Compound-specific parameters for light-dependent temperature activity
      ! factor (dbm, 6/21/2012)
      REAL(hp), INTENT(IN) :: CT1, CEO
!
! !RETURN VALUE:
!
      ! Temperature activity factor for the light-dependent fraction of 
      ! emissions
      REAL(hp)             :: GAMMA_T_LD
!
! !REMARKS:
!  References (see above for full citations):
!  (1 ) Guenther et al, 1995
!  (2 ) Guenther et al, 2006
!  (3 ) Guenther et al, MEGAN v2.1 user mannual 2007-08
!  (4 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY: 
!  (1 ) Includes the latest MEGAN v2.1 temperature algorithm (mpb, 2009).
!       Note, this temp-dependence is the same for the PCEEA & hybrid models.
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  (2 ) Modified to gamma_t_ld and to permit compound specific parameters
!       CT1 and Ceo (dbm, 6/21/2012)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)              :: C_T,   CT2
      REAL(hp)              :: E_OPT, T_OPT, X
!
! !DEFINED PARAMETERS:
!
      ! Ideal gas constant [J/mol/K]
      REAL(hp), PARAMETER   :: R   = 8.314e-3_hp

      !=================================================================
      ! GET_GAMMA_T_LD begins here!
      !================================================================= 
      E_OPT = CEO * EXP( 0.08_hp * ( PT_15  - 2.97e2_hp ) )
      T_OPT = 3.13e2_hp + ( 6.0e-1_hp * ( PT_15 - 2.97e2_hp ) )
      CT2   = 200.0_hp

      ! Variable related to temperature 
      X     = ( 1.0_hp/T_OPT - 1.0_hp/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / 
     &        ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      GAMMA_T_LD = MAX( C_T, 0.0_hp )

      END FUNCTION GET_GAMMA_T_LD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Lai
!
! !DESCRIPTION: Function Get\_Gamma\_Lai computes the gamma exchange activity 
!  factor which is sensitive to leaf area (= GAMMA\_LAI).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_LAI( CMLAI, BIDIREXCH ) RESULT( GAMMA_LAI )
!
! !INPUT PARAMETERS: 
!
      REAL(hp), INTENT(IN) :: CMLAI       ! Current month's LAI [cm2/cm2]
      LOGICAL,  INTENT(IN) :: BIDIREXCH   ! Logical flag indicating whether
                                          ! the compound undergoes bidirectional
                                          ! exchange
!
! !RETURN VALUE:
!
      REAL(hp)             :: GAMMA_LAI 
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest MEGAN v2.1 User's Guide (mpb,2009)
!  (4 ) Updated to treat bidirectional exchange compounds appropriately (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !-----------------------
      ! Compute GAMMA_LAI
      !-----------------------

      ! Formulation for birectional compounds is as described for 
      ! ALD2 in Millet et al., ACP 2010
      IF ( BIDIREXCH ) THEN
         
         IF ( CMLAI <= 6.0_hp) THEN
               
            ! if lai less than 2:
            IF ( CMLAI <= 2.0_hp ) THEN
               GAMMA_LAI = 0.5_hp * CMLAI
               
            ! if between 2 and 6:
            ELSE
               GAMMA_LAI = 1.0_hp - 0.0625_hp * ( CMLAI - 2.0_hp )
            END IF

         ELSE
            ! keep at 0.75 for LAI > 6
            GAMMA_LAI = 0.75_hp
         END IF

      ! For all other compounds use the standard gamma_lai formulation
      ELSE
         GAMMA_LAI = 0.49_hp * CMLAI / SQRT( 1.0_hp + 0.2_hp *
     &               CMLAI*CMLAI)
      ENDIF

      END FUNCTION GET_GAMMA_LAI 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Age
!
! !DESCRIPTION: Function Get\_Gamma\_Age computes the gamma exchange 
!  activity factor which is sensitive to leaf age (= Gamma\_Age).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_AGE( CMLAI, PMLAI, DBTWN, TT,
     &                        AN, AG, AM, AO )
     &         RESULT( GAMMA_AGE )
!
! !INPUT PARAMETERS: 
!
      REAL(hp), INTENT(IN) :: CMLAI     ! Current month's LAI [cm2/cm2]
      REAL(hp), INTENT(IN) :: PMLAI     ! Previous months LAI [cm2/cm2]
      REAL(hp), INTENT(IN) :: DBTWN     ! Number of days between 
      REAL(sp), INTENT(IN) :: TT        ! Daily average temperature [K]
      REAL(hp), INTENT(IN) :: AN        ! Relative emiss factor (new leaves)
      REAL(hp), INTENT(IN) :: AG        ! Relative emiss factor (growing leaves)
      REAL(hp), INTENT(IN) :: AM        ! Relative emiss factor (mature leaves)
      REAL(hp), INTENT(IN) :: AO        ! Relative emiss factor (old leaves)
!
! !RETURN VALUE:
!
      REAL(hp)             :: GAMMA_AGE ! Activity factor    
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003). Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest User's Guide (tmf, 11/19/04)
!  (4 ) Renamed & now includes specific relative emission activity factors for
!       each BVOC based on MEGAN v2.1 algorithm (mpb,2008)
!  (5 ) Now calculate TI (number of days after budbreak required to induce 
!       iso. em.) and TM (number of days after budbreak required to reach 
!       peak iso. em. rates) using the daily average temperature, instead 
!       of using fixed values (mpb,2008)
!       NOTE: Can create 20% increases in tropics (Guenther et al 2006)
!  (6 ) Implemented change for the calculation of FGRO if ( CMLAI > PMLAI ),
!       i.e. if LAI has increased with time, and used new values for 
!       all foilage fractions if ( CMLAI = PMLAI ). Also removed TG variable 
!       as not now needed. (mpb,2000)
!  (7 ) Updated to pass leaf age activity factors as arguments (dbm, 6/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Updated for sesquiterpenes (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)             :: FNEW  ! Fraction of new leaves in canopy
      REAL(hp)             :: FGRO  ! Fraction of growing leaves
      REAL(hp)             :: FMAT  ! Fraction of mature leaves
      REAL(hp)             :: FOLD  ! Fraction of old leaves

      ! TI: number of days after budbreak required to induce emissions
      REAL(hp)             :: TI   

      ! TM: number of days after budbreak required to reach peak emissions
      REAL(hp)             :: TM

      !=================================================================
      ! GET_GAMMA_AGE begins here!
      !================================================================= 
      
      !-----------------------
      ! Compute TI and TM 
      ! (mpb,2009)
      !-----------------------

      IF ( TT <= 303.0_hp ) THEN
         TI = 5.0_hp + 0.7_hp * ( 300.0_hp - TT )
      ELSEIF ( TT >  303.0_hp ) THEN
         TI = 2.9_hp
      ENDIF   
      TM = 2.3_hp * TI

      !-----------------------
      ! Compute GAMMA_AGE
      !-----------------------

      IF ( CMLAI == PMLAI ) THEN !(i.e. LAI stays the same) 

         FNEW = 0.0_hp        
         FGRO = 0.1_hp
         FMAT = 0.8_hp
         FOLD = 0.1_hp

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has INcreased) 

         ! Calculate Fnew
         IF ( DBTWN > TI ) THEN
            FNEW = ( TI / DBTWN ) * ( 1.0_hp -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.0_hp - ( PMLAI / CMLAI )
         ENDIF

         ! Calculate FMAT
         IF ( DBTWN > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +
     &            (( DBTWN - TM ) / DBTWN )*( 1.0_hp -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Calculate Fgro and Fold
         FGRO = 1.0_hp - FNEW - FMAT
         FOLD = 0.0_hp

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has DEcreased) 

         FNEW = 0.0_hp
         FGRO = 0.0_hp
         FOLD = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.0_hp - FOLD

      ENDIF

      ! Age factor
      GAMMA_AGE = FNEW * AN + FGRO * AG +
     &            FMAT * AM + FOLD * AO

      ! Prevent negative values
      GAMMA_AGE = MAX( GAMMA_AGE , 0.0_hp )

      END FUNCTION GET_GAMMA_AGE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gamma_sm
!
! !DESCRIPTION: Function GET\_GAMMA\_SM computes activity factor for soil
!  moisture
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_SM( ExtState, I, J, CMPD )
     &            RESULT( GAMMA_SM )
!
! !INPUT PARAMETERS: 
!
      TYPE(Ext_State),  POINTER     :: ExtState
      INTEGER,          INTENT(IN)  :: I, J     ! GEOS-Chem lon & lat indices
      CHARACTER(LEN=*), INTENT(IN)  :: CMPD     ! Compound name (dbm, 6/21/2012)
      
! !RETURN VALUE:
!
      REAL(hp)                      :: GAMMA_SM ! Activity factor    
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, ACP 2006
!  (2 ) Guenther et al., GMD 2012 and MEGANv2.1 source code
! 
! !REVISION HISTORY: 
!  (1 ) Created by dbm (6/2012). We are not currently using a soil moisture
!       effect for isoprene. For all compounds other than acetaldehyde and
!       ethanol, gamma_sm =1 presently.
!  16 Apr 2015 - C. Keller   - Now restrict GWETROOT to values between 0.0 and
!                              1.0. This only seems to be a problem within the
!                              GEOS-5 ESM, where GWETROOT values over the ocean
!                              become 1e+15 (= missing value). 
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp)  :: GWETROOT

      !=================================================================
      ! GET_GAMMA_SM begins here!
      !================================================================= 
      
      ! By default gamma_sm is 1.0
      GAMMA_SM = 1.0_hp

      ! GWETROOT is only defined for GEOS-5
#if defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA )

      !GWETROOT = ExtState%GWETROOT%Arr%Val(I,J)
      ! Error trap: GWETROOT must be between 0.0 and 1.0 (ckeller, 4/16/15)
      GWETROOT = MIN(MAX(ExtState%GWETROOT%Arr%Val(I,J),0.0_hp),1.0_hp)

      IF ( TRIM( CMPD ) == 'ALD2' .OR. TRIM ( CMPD ) == 'ETOH' ) THEN
         
         ! GWETROOT = degree of saturation or wetness in the root-zone
         ! (top meter of soil). This is defined as the ratio of the volumetric
         ! soil moisture to the porosity. We use a soil moisture activity factor
         ! for ALD2 to account for stimulation of emission by flooding.
         ! (Millet et al., ACP 2010)
         ! Constant value of 1.0 for GWETROOT = 0-0.9, increasing linearly to
         ! 3.0 at GWETROOT =1.
         GAMMA_SM = MAX( 20.0_hp * GWETROOT - 17.0_hp, 1.0_hp)

      ENDIF

#endif

      ! return to calling program
      END FUNCTION GET_GAMMA_SM

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_NORM_FAC
!
! !DESCRIPTION: Function CALC\_NORM\_FAC calculates the normalization factor 
!  needed to compute emissions. Called from GET_MEGAN_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_NORM_FAC( am_I_Root, RC )
!
! !INPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN)     :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS
!
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, (GMD 2012) and associated MEGANv2.1 source code 
! 
! !REVISION HISTORY: 
!  (1 ) Created by dbm 11/2012. We calculate only 1 normalization factor for all
!       compounds based on the isoprene gamma values. Formally there should be a
!       different normalization factor for each compound, but we are following
!       Alex Guenther's approach here and the MEGAN source code.
!       "Hi Dylan, sorry for being so slow to get back to you.
!        Since the change is only a few percent or less, I didn't
!        bother to assign a different normalization factor to each
!        compound.  Since the MEGAN canopy environment model also 
!        has 8 different canopy types (tropical broadleaf tree, 
!        conifer tree, etc.) then to be correct we should have a 
!        different CCE for each canopy type for each compound class
!        (which would be 160 slightly different values of CCE)."
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(hp) :: PAC_DAILY, PHI, BBB, AAA, GAMMA_P_STANDARD
      REAL(hp) :: GAMMA_T_LI_STANDARD
      REAL(hp) :: GAMMA_SM_STANDARD
      REAL(hp) :: CMLAI, GAMMA_LAI_STANDARD
      REAL(hp) :: GAMMA_AGE_STANDARD
      REAL(hp) :: PT_15, T, R, CEO, CT1, E_OPT, T_OPT, CT2, X
      REAL(hp) :: GAMMA_T_LD_STANDARD
      REAL(hp) :: LDF, GAMMA_STANDARD

      !-----------------------------------------------------------------
      ! CALC_NORM_FAC
      !-----------------------------------------------------------------  

      ! -----------------
      ! GAMMA_P for standard conditions
      ! -----------------
      ! Based on Eq. 11b from Guenther et al., 2006
      ! Using standard conditions of phi = 0.6, solar angle of 60 deg,
      ! and P_daily = 400
      ! Note corrigendum for Eq. 11b in that paper, should be a
      ! minus sign before the 0.9.
      PAC_DAILY = 400.0_hp
      PHI       = 0.6_hp
      BBB       = 1.0_hp + 0.0005_hp *( PAC_DAILY - 400.0_hp ) 
      AAA       = ( 2.46_hp * BBB * PHI ) - ( 0.9_hp * PHI**2 )
      ! sin(60) = 0.866
      GAMMA_P_STANDARD = 0.866_hp * AAA

      ! -----------------
      ! GAMMA_T_LI for standard conditions
      ! -----------------
      ! gamma_t_li = EXP( Beta * ( T - T_Standard ) )
      ! This is 1.0 for T = T_Standard
      GAMMA_T_LI_STANDARD = 1.0_hp

      ! -----------------
      ! GAMMA_SM for standard conditions
      ! -----------------
      ! Standard condition is soil moisture = 0.3 m^3/m^3
      ! GAMMA_SM = 1.0 for all compounds under this condition
      GAMMA_SM_STANDARD = 1.0_hp

      ! -----------------
      ! GAMMA_LAI for standard conditions 
      ! -----------------
      ! Standard condition is LAI = 5
      CMLAI = 5.0_hp
      GAMMA_LAI_STANDARD = 0.49_hp * 
     &     CMLAI / SQRT( 1.0_hp + 0.2_hp * CMLAI*CMLAI )
     
      ! -----------------
      ! GAMMA_AGE for standard conditions
      ! -----------------
      ! Standard condition is 0% new, 10% growing, 80% mature, 10% old foliage
      ! Isoprene uses A_NEW = 0.05d0, A_GRO = 0.6d0, A_MAT = 1.d0, A_OLD = 0.9d0
      GAMMA_AGE_STANDARD = 0.1_hp*0.6_hp + 0.8_hp*1.0_hp + 0.1_hp*0.9_hp

      ! -----------------
      ! GAMMA_T_LD for standard conditions
      ! -----------------
      ! Standard condition is
      ! PT_15 = average leaf temp over past 24-240 hours = 297K
      ! T = air temperature = 303K
      PT_15 = 297.0_hp
      T     = 303.0_hp
      R     = 8.314e-3_hp
      ! parameters for isoprene
      CEO = 2.0_hp
      CT1 = 95.0_hp

      E_OPT = CEO * EXP( 0.08_hp * ( PT_15  - 2.97e2_hp ) )           
      T_OPT = 3.13e2_hp + ( 6.0e-1_hp * ( PT_15 - 2.97e2_hp ) )
      CT2   = 200.0_hp

      ! Variable related to temperature 
      X     = ( 1.0_hp/T_OPT - 1.0_hp/T ) / R

      GAMMA_T_LD_STANDARD = E_OPT * CT2 * EXP( CT1 * X ) / 
     &     ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

      ! -----------------
      ! Overall GAMMA_STANDARD
      ! -----------------
      ! LDF = 1d0 for isoprene
      LDF = 1.0_hp
      GAMMA_STANDARD = 
     &     GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD * GAMMA_LAI_STANDARD
     &     * ((1.0_hp - LDF) * GAMMA_T_LI_STANDARD 
     &     + (LDF * GAMMA_P_STANDARD * GAMMA_T_LD_STANDARD))
      ! This ends up being 1.0101081.
      
      NORM_FAC = 1.0_hp / GAMMA_STANDARD

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE CALC_NORM_FAC
!EOC

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fill_Restart_Vars 
!
! !DESCRIPTION: Subroutine FILL\_RESTART\_VARS fills the megan restart
!  variables.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FILL_RESTART_VARS( am_I_Root, HcoState, ExtState, RC )
!
! !USES:
!
      USE HCO_Restart_Mod,    ONLY : HCO_RestartGet
!
! !INPUT PARAMETERS: 
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root 
      TYPE(HCO_State), POINTER        :: HcoState 
      TYPE(Ext_State), POINTER        :: ExtState
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,         INTENT(INOUT)  :: RC
! 
! !REVISION HISTORY: 
!  19 Dec 2014 - C. Keller    - Initial version 
!  09 Mar 2015 - C. Keller    - Now use HCO_RestartGet to fill restarts 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I
      LOGICAL           :: FND

      !=================================================================
      ! FILL_RESTART_VARS begins here!
      !=================================================================

      ! On first call, see if there are restart variables available for
      ! historic values of temperature and irradiation. If so, fill the
      ! daily arrays with those average values to start with. If not,
      ! simply use the current values.
      ! The same procedure is applied to the 15 day average values. For
      ! temperature, we also restart previous day temperature so that
      ! we can determine a trend in the daily average temperatures
      ! (between the long term daily averages and the most recent daily
      ! average).
      ! We now first call HCO_CopyFromInt_ESMF to check if the variable
      ! is found in the (ESMF) internal state object. In an non-ESMF
      ! environment, this call will always return FND=.FALSE.
      ! (ckeller, 3/9/15).

      !-----------------------------------------------------------------
      ! LAI
      !-----------------------------------------------------------------
      CALL HCO_RestartGet( am_I_Root,     HcoState, 
     &                     'LAI_PREVDAY', LAI_PREVDAY, RC, FILLED=FND )

      ! Default value
      IF ( .NOT. FND ) THEN
         LAI_PREVDAY = ExtState%LAI%Arr%Val
      ENDIF

      !-----------------------------------------------------------------
      ! TEMPERATURE
      !-----------------------------------------------------------------
      CALL HCO_RestartGet( am_I_Root, HcoState, 
     &                    'T_DAVG',  HCOT_15_AVG, RC, FILLED=FND )

      ! Default value
      IF ( .NOT. FND ) THEN
         HCOT_15_AVG = ExtState%T2M%Arr%Val
      ENDIF

      ! Set all daily averages to this value. 
      DO I = 1, NUM_DAYS
         HCOT_15(:,:,I) = HCOT_15_AVG 
      ENDDO 

      ! -------------------------------------------------------------
      ! Check for most recent daily average. This will overwrite the
      ! first slice of HCOT_15.
      CALL HCO_RestartGet( am_I_Root,   HcoState, 
     &                    'T_PREVDAY', HCOT_15(:,:,1), RC, FILLED=FND )

      ! Default value
      IF ( .NOT. FND ) THEN
         HCOT_15(:,:,1) = HCOT_15_AVG
      ENDIF

      ! Set hourly variables
      DO I = 1, DAY_DIM 
         HCOT_DAY(:,:,I) = HCOT_15(:,:,1)
      ENDDO

      !-----------------------------------------------------------------
      ! DIRECT RADIATION ( PARDR )
      !-----------------------------------------------------------------
      CALL HCO_RestartGet( am_I_Root,    HcoState, 
     &                    'PARDR_DAVG', HCOPARDR_15_AVG, RC, FILLED=FND)

      ! Default value
      IF ( .NOT. FND ) THEN
         HCOPARDR_15_AVG = ExtState%PARDR%Arr%Val
      ENDIF

      ! Set daily averages
      DO I = 1, NUM_DAYS
         HCOPARDR_15(:,:,I) = HCOPARDR_15_AVG
      ENDDO 

      ! Set hourly variables
      DO I = 1, DAY_DIM 
         HCOPARDR_DAY(:,:,I) = HCOPARDR_15(:,:,1) 
      ENDDO

      !-----------------------------------------------------------------
      ! DIFFUSE RADIATION ( PARDF )
      !-----------------------------------------------------------------
      CALL HCO_RestartGet( am_I_Root,    HcoState, 
     &                    'PARDF_DAVG', HCOPARDF_15_AVG, RC, FILLED=FND)

      ! Default value
      IF ( .NOT. FND ) THEN
         HCOPARDF_15_AVG = ExtState%PARDF%Arr%Val
      ENDIF

      ! Set daily averages
      DO I = 1, NUM_DAYS
         HCOPARDF_15(:,:,I) = HCOPARDF_15_AVG
      ENDDO 

      ! Set hourly variables
      DO I = 1, DAY_DIM 
         HCOPARDF_DAY(:,:,I) = HCOPARDF_15(:,:,1)
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS
       
      END SUBROUTINE FILL_RESTART_VARS 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_HCOT_Day
!
! !DESCRIPTION: Subroutine Update\_HCOT\_Day must be called every time the 
!  A-3 fields are updated. Each 3h TS value for each gridbox is moved up one 
!  spot in the matrix and the current value is put in the last spot.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE UPDATE_HCOT_DAY( HcoState, ExtState, I, J, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER        :: HcoState 
      TYPE(Ext_State), POINTER        :: ExtState
      INTEGER,         INTENT(IN   )  :: I 
      INTEGER,         INTENT(IN   )  :: J 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,         INTENT(INOUT)  :: RC
! 
! !REVISION HISTORY: 
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: D

      !=================================================================
      ! UPDATE_HCOT_DAY begins here!
      !=================================================================

      ! Move each day up
      DO D = DAY_DIM, 2, -1
         ! Need PAR as well as Temp  (mpb,2009)
         HCOT_DAY(I,J,D)     = HCOT_DAY(I,J,D-1)
         HCOPARDR_DAY(I,J,D) = HCOPARDR_DAY(I,J,D-1)
         HCOPARDF_DAY(I,J,D) = HCOPARDF_DAY(I,J,D-1)
      ENDDO
      
      ! Store 
      HCOT_DAY(I,J,1)     = ExtState%T2M%Arr%Val(I,J)
      HCOPARDF_DAY(I,J,1) = ExtState%PARDF%Arr%Val(I,J)
      HCOPARDR_DAY(I,J,1) = ExtState%PARDR%Arr%Val(I,J)

      ! Return w/ success
      RC = HCO_SUCCESS
       
      END SUBROUTINE UPDATE_HCOT_DAY
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Update_HCOT_15_Avg
!
! !DESCRIPTION: Subroutine Update\_HCOT\_15\_Avg should be called at the 
!  beginning of each day. It loops through the gridboxes doing the following:
!
!  \begin{enumerate}
!  \item Average T\_DAY over the 8 TS values during the day.  
!  \item Push the daily average TS values through T\_15, throwing out the 
!        oldest and putting the newest (the T\_DAY average) in the last spot 
!  \item Get T\_15\_AVG by averaging T\_15 over the 15 day period. 
!  \end{enumerate}
!
! !INTERFACE:
!
      SUBROUTINE UPDATE_HCOT_15_AVG ( HcoState, I, J, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER       :: HcoState
      INTEGER,         INTENT(IN   ) :: I 
      INTEGER,         INTENT(IN   ) :: J 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,         INTENT(INOUT) :: RC
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - M. Prather  - Initial version
!  (1 ) All MEGAN biogenic emission are currently calculated using TS from DAO 
!        met field. TS is the surface air temperature, which should be 
!        carefully distinguished from TSKIN. (tmf, 11/20/04)
!  (2 ) In GEOS4, TS are originally T2M in the A3 files, read in 
!        'a3_read_mod.f'. 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: D
      REAL(hp)          :: D_DIM, D_DAYS, TMP_T
      REAL(hp)          :: TMP_PARDR , TMP_PARDF ! (mpb,2009)

      !=================================================================
      ! UPDATE_HCOT_15_AVG begins here!
      !=================================================================

      ! Convert to REAL*8
      D_DIM  = REAL( DAY_DIM , KIND=hp )
      D_DAYS = REAL( NUM_DAYS, KIND=hp )

      ! Average HCOT_DAY over the TS values during the day.
      TMP_T = SUM( HCOT_DAY(I,J,:) ) / D_DIM
         
      ! Do the same for light (mpb,2009)
      TMP_PARDR = SUM( HCOPARDR_DAY(I,J,:) ) / D_DIM
      TMP_PARDF = SUM( HCOPARDF_DAY(I,J,:) ) / D_DIM

      ! Push the daily average TS values through HCOT_15,
      ! throwing out the oldest 
      DO D = NUM_DAYS, 2, -1
         HCOT_15(I,J,D)     = HCOT_15(I,J,D-1)
         HCOPARDR_15(I,J,D) = HCOPARDR_15(I,J,D-1)
         HCOPARDF_15(I,J,D) = HCOPARDF_15(I,J,D-1)
      ENDDO

      ! Put the newest daily average TS value in the first spot
      HCOT_15(I,J,1) = TMP_T

      ! Get HCOT_15_AVG by averaging HCOT_15 over the 15 day period.
      HCOT_15_AVG(I,J) = SUM( HCOT_15(I,J,:) ) / D_DAYS 

      ! Assign daily average temperature to HCOT_DAILY (mpb,2009)
      ! HCOT_DAILY(I,J)  = TMP_T

      ! Repeat for PAR diffuse & direct (mpb,2009)
      HCOPARDR_15(I,J,1)   = TMP_PARDR
      HCOPARDR_15_AVG(I,J) = SUM( HCOPARDR_15(I,J,:) ) / D_DAYS 
      ! HCOPARDR_DAILY(I,J)  = TMP_PARDR

      HCOPARDF_15(I,J,1)   = TMP_PARDF
      HCOPARDF_15_AVG(I,J) = SUM( HCOPARDF_15(I,J,:) ) / D_DAYS 
      ! HCOPARDF_DAILY(I,J)  = TMP_PARDF

      ! Return w/ success 
      RC = HCO_SUCCESS

      END SUBROUTINE UPDATE_HCOT_15_AVG
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_AEF
!
! !DESCRIPTION: Subroutine CALC\_AEF reads Emission Factors for all
!  biogenic VOC species from disk. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_AEF( am_I_Root, HcoState, ExtState, RC )
!
! !USES:
!
      ! References to F90 modules
      USE HCO_EMISLIST_MOD,    ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS: 
!
      LOGICAL,         INTENT(IN)    :: am_I_Root
      TYPE(Ext_State), POINTER       :: ExtState 
      TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,         INTENT(INOUT) :: RC
!
! !REMARKS:
!  Reference: (5 ) Guenther et al, 2004 
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) AEF detailed in the latest MEGAN User's Guide (tmf, 11/19/04)
!  (3 ) Bug fix (tmf, 11/30/04)
!  (4 ) Now reads 1x1 files and regrids to current resolution (bmy, 10/24/05)
!  (5 ) Uses new v2.1 emission factors maps for isoprene, MBO and 7 monoterpene
!     species, download in 2009. (mpb,2009)
!  (6 ) Now use 2.1 emission factors for isoprene, MBO, and 7 monoterpenes. EFs for
!       other compounds are computed by reading in the PFT fractions and multiplying
!       the fractions by corresponding EF values. (dbm, 11/2012)
!  (7 ) Also, now read in the EF values already gridded to model resolution.
!       (dbm, 11/2012)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  10 Jun 2015 - M. Sulprizio- Bug fix: Now convert sesquiterpenes (FARN, BCAR,
!                              OSQT) from ug compound/m2/h to kg C/m2/s
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                 :: I, J, P, ARR_IND
      REAL(hp)                :: FACTOR, SPECIES2CARBON
      REAL(hp)                :: ARRAY_16(HcoState%NX,HcoState%NY,16)
      REAL(hp)                :: PFT_EF_OMON(15), PFT_EF_MOHX(15)
      REAL(hp)                :: PFT_EF_ACET(15), PFT_EF_BIDR(15)
      REAL(hp)                :: PFT_EF_STRS(15), PFT_EF_OTHR(15)
! --->
! dbm, compute EF maps for a-pinene and myrcene as well since there seems to 
! be an issue with the EF maps for these species provided on the MEGAN
! data portal
      REAL(hp)                :: PFT_EF_APIN(15), PFT_EF_MYRC(15)
! <---
      REAL(hp)                :: PFT_EF_FARN(15), PFT_EF_BCAR(15)
      REAL(hp)                :: PFT_EF_OSQT(15)
      REAL(hp)                :: EM_FRAC_ALD2(15), EM_FRAC_ETOH(15)
      REAL(hp)                :: EM_FRAC_FAXX(15), EM_FRAC_AAXX(15)
      REAL(hp)                :: EM_FRAC_CH2O(15)
      
      !=================================================================
      ! CALC_AEF begins here!
      !=================================================================

      ! ----------------------------------------------------------------
      ! Note that not all these compounds are used in standard chemistry
      ! simulations, but they are provided here for future incorporation or 
      ! specialized applications. More compounds can be added as needed
      ! by adding the corresponding PFT-specific emission factors and 
      ! emission category fraction.
      ! (dbm, 01/2013)
      ! ----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Point to external data 
      !-----------------------------------------------------------------

      CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_ISOP', AEF_ISOP, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

! ---> Now compute EF maps for acetone (dbm, 12/2012)
!     CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_ACET', AEF_ACET, RC )
!     IF ( RC /= HCO_SUCCESS ) RETURN
! <---

      ! Arrays used if monoterpene switch is turned on
      IF ( ExtNrMono > 0 ) THEN
         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_MBOX', AEF_MBOX, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      
! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
!         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_APIN',AEF_APIN, RC )
!         IF ( RC /= HCO_SUCCESS ) RETURN
! <---

         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_BPIN',AEF_BPIN, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_CARE',AEF_CARE, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_LIMO',AEF_LIMO, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

! ---> Now compute EF maps for a-pinene and myrcene (dbm, 12/2012)
!         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_MYRC',AEF_MYRC, RC )
!         IF ( RC /= HCO_SUCCESS ) RETURN
! <---

         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_OCIM',AEF_OCIM, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_SABI',AEF_SABI, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Only used for SOA simulation:
      IF ( ExtNrSoa > 0 ) THEN
         CALL HCO_GetPtr( am_I_Root, 'MEGAN_ORVC', GEIA_ORVC, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Point to PFT fractions
      !-----------------------------------------------------------------

      ! CLM4 PFT coverage (unitless)
      ! From Table 3 in Guenther et al., 2012
      ! PFT_BARE                : Bare
      ! PFT_NDLF_EVGN_TMPT_TREE : Needleleaf evergreen temperate tree
      ! PFT_NDLF_EVGN_BORL_TREE : Needleleaf evergreen boreal tree
      ! PFT_NDLF_DECD_BORL_TREE : Needleleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_TROP_TREE : Broadleaf evergreen tropical tree
      ! PFT_BDLF_EVGN_TMPT_TREE : Broadleaf evergreen temperate tree
      ! PFT_BDLF_DECD_TROP_TREE : Broadleaf deciduous tropical tree
      ! PFT_BDLF_DECD_TMPT_TREE : Broadleaf deciduous temperate tree
      ! PFT_BDLF_DECD_BORL_TREE : Broadleaf deciduous boreal tree
      ! PFT_BDLF_EVGN_SHRB      : Broadleaf evergreen temperate shrub
      ! PFT_BDLF_DECD_TMPT_SHRB : Broadleaf deciduous temperate shrub
      ! PFT_BDLF_DECD_BORL_SHRB : Broadleaf deciduous boreal shrub
      ! PFT_C3_ARCT_GRSS        : Arctic C3 grass
      ! PFT_C3_NARC_GRSS        : Cool C3 grass
      ! PFT_C4_GRSS             : Warm C4 grass
      ! PFT_CROP                : Crop

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BARE', PFT_BARE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_NDLF_EVGN_TMPT_TREE',
     &                 PFT_NDLF_EVGN_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_NDLF_EVGN_BORL_TREE',
     &                 PFT_NDLF_EVGN_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_NDLF_DECD_BORL_TREE',
     &                 PFT_NDLF_DECD_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_EVGN_TROP_TREE',
     &                 PFT_BDLF_EVGN_TROP_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_EVGN_TMPT_TREE',
     &                 PFT_BDLF_EVGN_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_DECD_TROP_TREE',
     &                 PFT_BDLF_DECD_TROP_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_DECD_TMPT_TREE',
     &                 PFT_BDLF_DECD_TMPT_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_DECD_BORL_TREE',
     &                 PFT_BDLF_DECD_BORL_TREE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_EVGN_SHRB',
     &                 PFT_BDLF_EVGN_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_DECD_TMPT_SHRB',
     &                 PFT_BDLF_DECD_TMPT_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_BDLF_DECD_BORL_SHRB',
     &                 PFT_BDLF_DECD_BORL_SHRB, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_C3_ARCT_GRSS',
     &                 PFT_C3_ARCT_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_C3_NARC_GRSS',
     &                 PFT_C3_NARC_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_C4_GRSS',
     &                 PFT_C4_GRSS, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'CLM4_PFT_CROP', PFT_CROP, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Copy PFTs into ARRAY_16
      ARRAY_16(:,:, 1) = PFT_BARE
      ARRAY_16(:,:, 2) = PFT_NDLF_EVGN_TMPT_TREE
      ARRAY_16(:,:, 3) = PFT_NDLF_EVGN_BORL_TREE
      ARRAY_16(:,:, 4) = PFT_NDLF_DECD_BORL_TREE
      ARRAY_16(:,:, 5) = PFT_BDLF_EVGN_TROP_TREE
      ARRAY_16(:,:, 6) = PFT_BDLF_EVGN_TMPT_TREE
      ARRAY_16(:,:, 7) = PFT_BDLF_DECD_TROP_TREE
      ARRAY_16(:,:, 8) = PFT_BDLF_DECD_TMPT_TREE
      ARRAY_16(:,:, 9) = PFT_BDLF_DECD_BORL_TREE
      ARRAY_16(:,:,10) = PFT_BDLF_EVGN_SHRB
      ARRAY_16(:,:,11) = PFT_BDLF_DECD_TMPT_SHRB
      ARRAY_16(:,:,12) = PFT_BDLF_DECD_BORL_SHRB
      ARRAY_16(:,:,13) = PFT_C3_ARCT_GRSS
      ARRAY_16(:,:,14) = PFT_C3_NARC_GRSS
      ARRAY_16(:,:,15) = PFT_C4_GRSS
      ARRAY_16(:,:,16) = PFT_CROP

      ! --------------------------------------------------------------------------------
      ! PFT-specific EFs from Table 2 in Guenther et al., 2012
      ! in ug compound/m2/h
      ! PFTs 1-15 in the table correspond to #2-16 
      ! (i.e., excluding bare ground #1) in the above array.
      ! --------------------------------------------------------------------------------
      ! Compound Class EF1 EF2 EF3 EF4 EF5 EF6 EF7 EF8 EF9 EF10 EF11 EF12 EF13 EF14 EF15
      ! --------------------------------------------------------------------------------
      ! Other Monoterp 180 180 170 150 150 150 150 150 110 200  110  5    5    5    5
      ! Methanol       900 900 900 500 900 500 900 900 900 900  900  500  500  500  900
      ! Acetone        240 240 240 240 240 240 240 240 240 240  240  80   80   80   80
      ! Bidirect VOC   500 500 500 500 500 500 500 500 500 500  500  80   80   80   80
      ! Stress VOC     300 300 300 300 300 300 300 300 300 300  300  300  300  300  300
      ! Other VOC      140 140 140 140 140 140 140 140 140 140  140  140  140  140  140  
      ! a-Pinene       500 500 510 600 400 600 400 400 200 300  200    2    2    2    2
      ! Myrcene         70  70  60  80  30  80  30  30  30  50   30  0.3  0.3  0.3  0.3
      ! a-Farnesene     40  40  40  60  40  60  40  40  40  40   40    3    3    3    4
      ! b-Carophyllene  80  80  80  60  40  60  40  40  50  50   50    1    1    1    4
      ! Other sesqt.   120 120 120 120 100 120 100 100 100 100  100    2    2    2    2
      ! --------------------------------------------------------------------------------

      ! One thing to note is these are net emissions to the canopy atmosphere
      ! but not net emissions to the above canopy atmosphere since they don't
      !  account for within-canopy deposition. Only an issue for OVOCs.
      
!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_OMON = (/180.0_hp, 180.0_hp, 170.0_hp, 150.0_hp, 150.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                150.0_hp, 150.0_hp, 150.0_hp, 110.0_hp, 200.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                110.0_hp, 5.0_hp  , 5.0_hp  , 5.0_hp  , 5.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_MOHX = (/900.0_hp, 900.0_hp, 900.0_hp, 500.0_hp, 900.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                500.0_hp, 900.0_hp, 900.0_hp, 900.0_hp, 900.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                900.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 900.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_ACET = (/240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                240.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_BIDR = (/500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                500.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_STRS = (/300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_OTHR = (/140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp/)

! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_APIN = (/500.0_hp, 500.0_hp, 510.0_hp, 600.0_hp, 400.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                600.0_hp, 400.0_hp, 400.0_hp, 200.0_hp, 300.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                200.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_MYRC = (/70.0_hp,  70.0_hp,  60.0_hp,  80.0_hp,  30.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                80.0_hp,  30.0_hp,  30.0_hp,  30.0_hp,  50.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                30.0_hp,  0.3_hp,   0.3_hp,   0.3_hp,   0.3_hp/)
! <---

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_FARN = (/40.0_hp,  40.0_hp,  40.0_hp,  60.0_hp,  40.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                60.0_hp,  40.0_hp,  40.0_hp,  40.0_hp,  40.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                40.0_hp,  3.0_hp,   3.0_hp,   3.0_hp,   4.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_BCAR = (/80.0_hp,  80.0_hp,  80.0_hp,  60.0_hp,  40.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                60.0_hp,  40.0_hp,  40.0_hp,  50.0_hp,  50.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                50.0_hp,  1.0_hp,   1.0_hp,   1.0_hp,   4.0_hp/)

!                     EF1       EF2       EF3       EF4       EF5    
      PFT_EF_OSQT = (/120.0_hp, 120.0_hp, 120.0_hp, 120.0_hp, 100.0_hp, 
!                     EF6       EF7       EF8       EF9       EF10
     &                120.0_hp, 100.0_hp, 100.0_hp, 100.0_hp, 100.0_hp, 
!                     EF11      EF12      EF13      EF14      EF15
     &                100.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)


      ! Other monoterpenes, methanol, acetone, MBO are each 100% of thier
      ! respective categories. The VOCs within the stress category each 
      ! account for a specific fraction of emissions across all PFTs
      ! (ethene 58%, toluene 3%, HCN 1.5%). The VOCs within the 
      ! other category also account for a given fraction of emissions
      ! across all PFTs (propene 48%, butene 24%, other alkenes 0.2%). But
      ! VOCs in the bidirectional category account for a different amount of
      ! the total flux for the different PFTs. So in this case we define a
      ! vector containing these fractions.
      
      ! Acetaldehyde: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_ALD2 = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

      ! Ethanol: 40% of bidirectional category flux, except 25%
      ! for grasses and crops
      EM_FRAC_ETOH = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp,
     &                 0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

      ! Formic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_FAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

      ! Acetic acid: 6% of bidirectional category flux, except 15%
      ! for grasses and crops
      EM_FRAC_AAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp,
     &                 0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

      ! Formaldehyde: 8% of bidirectional category flux, except 20%
      ! for grasses and crops
      EM_FRAC_CH2O = (/0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp,
     &                 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp,
     &                 0.08_hp, 0.20_hp, 0.20_hp, 0.20_hp, 0.20_hp/)

      ! Loop through plant types
      DO P = 1, 15

         ! Add 1 to Array_16 index to skip bare ground
         ARR_IND = P + 1
         
         ! Don't need to divide ARRAY_16 by 100 since data from netCDF
         ! file has already been converted to fraction (mps, 2/12/15)

         IF ( ExtNrMono > 0 ) THEN

! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
            ! a-pinene: 100% of category
            AEF_APIN(:,:) = AEF_APIN(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_APIN(P)

            ! Myrcene: 100% of category
            AEF_MYRC(:,:) = AEF_MYRC(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_MYRC(P)
! <---
            ! Other monoterpenes: 100% of category
            AEF_OMON(:,:) = AEF_OMON(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_OMON(P)

            ! a-Farnesene: 100% of category
            AEF_FARN(:,:) = AEF_FARN(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_FARN(P)

            ! b-Caryophyllene: 100% of category
            AEF_BCAR(:,:) = AEF_BCAR(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_BCAR(P)

            ! Other sesquiterpenes: 100% of category
            AEF_OSQT(:,:) = AEF_OSQT(:,:) + 
     &           ARRAY_16(:,:,ARR_IND) * PFT_EF_OSQT(P)

         ENDIF

         ! Methanol: 100% of category
         AEF_MOHX(:,:) = AEF_MOHX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_MOHX(P)

         ! Acetone: 100% of category
         AEF_ACET(:,:) = AEF_ACET(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_ACET(P)

         ! Ethanol: variable fraction of category
         AEF_ETOH(:,:) = AEF_ETOH(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_BIDR(P) * EM_FRAC_ETOH(P)

         ! Formaldehyde: variable fraction of category
         AEF_CH2O(:,:) = AEF_CH2O(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_BIDR(P) * EM_FRAC_CH2O(P)

         ! Acetaldehyde: variable fraction of category
         AEF_ALD2(:,:) = AEF_ALD2(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_BIDR(P) * EM_FRAC_ALD2(P)

         ! Formic acid: variable fraction of category
         AEF_FAXX(:,:) = AEF_FAXX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_BIDR(P) * EM_FRAC_FAXX(P)

         ! Acetic acid: variable fraction of category
         AEF_AAXX(:,:) = AEF_AAXX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_BIDR(P) * EM_FRAC_AAXX(P)

         ! Ethene: 58% of "stress" category for all PFTs
         AEF_C2H4(:,:) = AEF_C2H4(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.58_hp

         ! Toluene: 3% of "stress" category for all PFTs
         AEF_TOLU(:,:) = AEF_TOLU(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.03_hp

         ! HCN: 1.5% of "stress" category for all PFTs
         AEF_HCNX(:,:) = AEF_HCNX(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.015_hp

         ! Propene: 48% of "other" category for all PFTs
         ! Butene:  24% of "other" category for all PFTs
         ! Larger alkenes: 0.2% of "other" category for all PFTs
         ! Total: 72.2% 
         AEF_PRPE(:,:) = AEF_PRPE(:,:) + 
     &        ARRAY_16(:,:,ARR_IND) * PFT_EF_OTHR(P) * 0.722_hp
         
      ENDDO

      ! Nullify pointers
      PFT_BARE                => NULL()
      PFT_NDLF_EVGN_TMPT_TREE => NULL()
      PFT_NDLF_EVGN_BORL_TREE => NULL()
      PFT_NDLF_DECD_BORL_TREE => NULL()
      PFT_BDLF_EVGN_TROP_TREE => NULL()
      PFT_BDLF_EVGN_TMPT_TREE => NULL()
      PFT_BDLF_DECD_TROP_TREE => NULL()
      PFT_BDLF_DECD_TMPT_TREE => NULL()
      PFT_BDLF_DECD_BORL_TREE => NULL()
      PFT_BDLF_EVGN_SHRB      => NULL()
      PFT_BDLF_DECD_TMPT_SHRB => NULL()
      PFT_BDLF_DECD_BORL_SHRB => NULL()
      PFT_C3_ARCT_GRSS        => NULL()
      PFT_C3_NARC_GRSS        => NULL()
      PFT_C4_GRSS             => NULL()
      PFT_CROP                => NULL()

      ! Conversion factor from [ug compound/m2/hr] to [kg compound/m2/s]
      FACTOR = 1.0e-9_hp / 3600.0_hp
          
      ! Loop over grid boxes
      DO J = 1, HcoState%NY 
      DO I = 1, HcoState%NX

         ! Convert AEF arrays to [kgC/m2/s]
         ! Multiply arrays by FACTOR and ratio [g C/g compound]
         ! NOTE: AEFs for ISOP, MBOX, BPIN, CARE, LIMO, OCIM, SABI 
         ! are read from file in [kgC/m2/s], so no need to convert here
         IF ( ExtNrMono > 0 ) THEN
            AEF_APIN(I,J) = AEF_APIN(I,J) 
     &                    * FACTOR * 120.0_hp / 136.234_hp
            AEF_MYRC(I,J) = AEF_MYRC(I,J) 
     &                    * FACTOR * 120.0_hp / 136.234_hp
            AEF_OMON(I,J) = AEF_OMON(I,J) 
     &                    * FACTOR * 120.0_hp / 136.234_hp

            ! Sesquiterpenes
            SPECIES2CARBON = 15.d0 * 12.01d0 / 
     &                     ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
            AEF_FARN(I,J) = AEF_FARN(I,J) * FACTOR * SPECIES2CARBON
            AEF_BCAR(I,J) = AEF_BCAR(I,J) * FACTOR * SPECIES2CARBON
            AEF_OSQT(I,J) = AEF_OSQT(I,J) * FACTOR * SPECIES2CARBON

         ENDIF

         AEF_ACET(I,J) = AEF_ACET(I,J) * FACTOR *  36.0_hp /  58.079_hp 
         AEF_ETOH(I,J) = AEF_ETOH(I,J) * FACTOR *  24.0_hp /  46.068_hp 
         AEF_ALD2(I,J) = AEF_ALD2(I,J) * FACTOR *  24.0_hp /  44.053_hp 
         AEF_C2H4(I,J) = AEF_C2H4(I,J) * FACTOR *  24.0_hp /  28.053_hp
         AEF_TOLU(I,J) = AEF_TOLU(I,J) * FACTOR *  84.0_hp /  92.138_hp 
         AEF_PRPE(I,J) = AEF_PRPE(I,J) * FACTOR *  36.0_hp /  42.080_hp

         ! Methanol, formaldehyde, formic acid, acetic acid, HCN are
         ! carried in kg, not kg C
         ! Convert AEF arrays to [kg/m2/s]
         AEF_MOHX(I,J) = AEF_MOHX(I,J) * FACTOR
         AEF_CH2O(I,J) = AEF_CH2O(I,J) * FACTOR
         AEF_FAXX(I,J) = AEF_FAXX(I,J) * FACTOR
         AEF_AAXX(I,J) = AEF_AAXX(I,J) * FACTOR
         AEF_HCNX(I,J) = AEF_HCNX(I,J) * FACTOR

      ENDDO
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE CALC_AEF
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Megan_Init 
!
! !DESCRIPTION: Subroutine HCOX\_Megan\_Init allocates and initializes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_Megan_Init( am_I_Root, HcoState, ExtName, 
     &                            ExtState,    RC                ) 
!
! !USES:
!
      USE HCO_STATE_MOD,    ONLY : Hco_GetHcoID
      USE HCO_STATE_MOD,    ONLY : HCO_GetExtHcoID
      USE HCO_ExtList_Mod,  ONLY : GetExtNr, GetExtOpt
      USE HCO_Restart_Mod,  ONLY : HCO_RestartDefine
!
! !INPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     
      TYPE(Ext_State),  POINTER        :: ExtState
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(INOUT)           :: RC
! 
! !REVISION HISTORY: 
!  (1 ) Change the logic in the #if block for G4AHEAD. (bmy, 12/6/05)
!  (2 ) Bug fix: skip Feb 29th if GCAP (phs, 9/18/07)
!  (3 ) Now call GET_AEF_05x0666 for GEOS-5 nested grids (yxw,dan,bmy, 11/6/08)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  26 Aug 2010 - R. Yantosca - Now reference merra_a1_mod.f
!  01 Sep 2010 - R. Yantosca - Now read in NUM_DAYS of sfc temp data (this had
!                              been hardwired to 15 days previously)
!  07 Feb 2011 - R. Yantosca - Fix typos: make sure to zero out the proper 
!                              HCOPARDF_* and HCOPARDR_* arrays after allocation
!  22 Nov 2011 - R. Yantosca - Do not use erroneous AEF's for nested grids
!  08 Feb 2012 - R. Yantosca - Now read surface temperature for GEOS-5.7.x
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  11 Apr 2012 - R. Yantosca - Now remove the call to INIT_LAI; we shall 
!                              initialize the LAI arrays from main.F
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  11 Apr 2013 - R. Yantosca - Now pass directory info with Input_Opt
!  18 Feb 2015 - M. Sulprizio- Now allocate AEF arrays for species not read
!                              from file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                        :: I, nSpc, AS, NX, NY
      INTEGER, ALLOCATABLE           :: HcoIDs(:)
      REAL*8                         :: PI_180
      REAL(hp), POINTER              :: Ptr2D(:,:) => NULL()
      CHARACTER(LEN=255)             :: MSG
      CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:) 
 
      !=================================================================
      ! HCOX_MEGAN_INIT begins here!
      !=================================================================

      ! Extension Nr.
      ExtNr = GetExtNr( TRIM(ExtName) )
      IF ( ExtNr <= 0 ) RETURN

      ! Enter
      CALL HCO_ENTER ( 'HCOX_Megan_Init (hcox_megan_mod.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !-----------------------------------------------------------------
      ! Set settings and species IDs
      !-----------------------------------------------------------------

      ! Read settings specified in configuration file
      ! Note: the specified strings have to match those in 
      !       the config. file!
      CALL GetExtOpt ( ExtNr, 'Isoprene scaling', 
     &                 OptValHp=ISOP_SCALING, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Verbose mode
      IF ( am_I_Root ) THEN
         MSG = 'Use MEGAN biogenic emissions (extension module)'
         CALL HCO_MSG( MSG, SEP1='-' )
         MSG = ' - Use the following species:'
         CALL HCO_MSG( MSG )
      ENDIF

      ! Get species IDs
      ! --> Assume that species are ordered ISOP, ACET, PRPE, C2H4 in
      !     config. file!
      CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
      IF ( RC /= HCO_SUCCESS ) RETURN
      
      ! Assign species IDs
      IDTISOP = -1
      IDTACET = -1
      IDTPRPE = -1
      IDTC2H4 = -1
      IDTALD2 = -1
      DO I = 1, nSpc
         SELECT CASE ( TRIM(SpcNames(I)) )
            CASE( 'ISOP' )
               IDTISOP = HcoIDs(I)
               WRITE(MSG,*) '   Isoprene   = ',TRIM(SpcNames(I)),IDTISOP
            CASE( 'ACET' )
               IDTACET = HcoIDs(I)
               WRITE(MSG,*) '   Acetone    = ',TRIM(SpcNames(I)),IDTACET
            CASE( 'PRPE' )
               IDTPRPE = HcoIDs(I)
               WRITE(MSG,*) '   C3 Alkenes = ',TRIM(SpcNames(I)),IDTPRPE
            CASE( 'ALD2' )
               IDTALD2 = HcoIDs(I)
               WRITE(MSG,*) '   ALD2       = ',TRIM(SpcNames(I)),IDTALD2
            CASE( 'C2H4' )
               IDTC2H4 = HcoIDs(I)
               WRITE(MSG,*) '   Ethene     = ',TRIM(SpcNames(I)),IDTC2H4
            CASE DEFAULT
               MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
               CALL HCO_ERROR( MSG, RC )
               RETURN
         END SELECT
    
         ! Verbose
         IF ( am_I_Root ) CALL HCO_MSG( MSG )
      ENDDO

      ! Verbose mode
      IF ( am_I_Root ) THEN
         WRITE(MSG,*) ' --> Isoprene scale factor is ', ISOP_SCALING
         CALL HCO_MSG( MSG )
      ENDIF

      ! Check for monoterpene option --> this one is required for 
      ! CO, OCPI, MONX
      IDTCO     = -1 
      IDTOCPI   = -1 
      IDTMONX   = -1 
      ExtNrMono = GetExtNr('MEGAN_Mono')
      IF ( ExtNrMono > 0 ) THEN

         ! Verbose
         IF ( am_I_Root ) THEN
            MSG = ' - MEGAN monoterpene option enabled' 
            CALL HCO_MSG( MSG )
         ENDIF

         ! Get HEMCO species IDs of CO, OCPI and MONX
         CALL HCO_GetExtHcoID( HcoState, ExtNrMono, HcoIDs,
     &                         SpcNames, nSpc,      RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN
    
         DO I = 1, nSpc
            SELECT CASE ( TRIM(SpcNames(I)) )
               CASE( 'CO' )
                  IDTCO = HcoIDs(I)
                  WRITE(MSG,*) '   CO         = ',
     &                         TRIM(SpcNames(I)),IDTCO
               CASE( 'OC', 'OCPI' )
                  IDTOCPI = HcoIDs(I)
                  WRITE(MSG,*) '   OC aerosol = ',
     &                         TRIM(SpcNames(I)),IDTOCPI
               CASE( 'MONX' )
                  IDTMONX = HcoIDs(I)
                  WRITE(MSG,*) '   Monoterp.  = ',
     &                         TRIM(SpcNames(I)),IDTMONX
               CASE DEFAULT
                  MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
                  CALL HCO_ERROR( MSG, RC )
                  RETURN
            END SELECT

            ! Verbose
            IF ( am_I_Root ) CALL HCO_MSG( MSG )
         ENDDO

         ! Check if at least one of the Mono species is defined 
         IF ( IDTCO   < 0 .AND. IDTOCPI < 0 .AND. IDTMONX < 0 ) THEN
            MSG = 'MEGAN_MONO option is anabled but none of the '// 
     &            'species is found (CO, OC/OCPI, MONX)!'
            CALL HCO_ERROR( MSG, RC )
            RETURN
         ENDIF

      !-----------------------------------------------------------------
      ELSE

         ! CO emissions and acetone emissions will be zero/incomplete!
         IF ( IDTACET > 0 ) THEN
            MSG = 'Cannot fully calculate MEGAN acetone ' //
     &            'emissions because monoterpene switch is off!'
            CALL HCO_WARNING ( MSG, RC )
         ENDIF
         MSG = 'Will not calculate MEGAN CO emissions ' //
     &         'because monoterpene switch is off!'
         CALL HCO_WARNING ( MSG, RC )
      ENDIF

      ! --> Check for SOA option. This option is only valid together
      !     w/ monoterpene!
      ! This further adds species ALPH/LIMO/ALCO
      IDTMTPA = -1 
      IDTMTPO = -1 
      IDTLIMO = -1 
      IDTSESQ = -1
      ExtNrSoa = GetExtNr('MEGAN_SOA')
      IF ( ExtNrSoa > 0 ) THEN
         IF ( ExtNrMono <= 0 ) THEN
            MSG = 'MEGAN SOA option only valid with monoterpene!'
            CALL HCO_ERROR ( MSG, RC )
            RETURN
         ENDIF

         ! Verbose
         IF ( am_I_Root ) THEN
            MSG = '   - SOA mechanism enabled:' 
            CALL HCO_MSG( MSG )
         ENDIF

         ! Get HEMCO species IDs of ALPH/LIMO/ALCO 
         CALL HCO_GetExtHcoID( HcoState, ExtNrSoa, HcoIDs,
     &                         SpcNames, nSpc,     RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN

         DO I = 1, nSpc
            SELECT CASE ( TRIM(SpcNames(I)) )
               CASE( 'MTPA' )
                  IDTMTPA = HcoIDs(I)
                  WRITE(MSG,*) '     a-,b-pinene    = ',
     &               TRIM(SpcNames(I)),IDTMTPA
               CASE( 'MTPO' )
                  IDTMTPO = HcoIDs(I)
                  WRITE(MSG,*) '     Other monoterp.= ',
     &               TRIM(SpcNames(I)),IDTMTPO
               CASE( 'LIMO' )
                  IDTLIMO = HcoIDs(I)
                  WRITE(MSG,*) '     Limonene       = ',
     &               TRIM(SpcNames(I)),IDTLIMO
               CASE( 'SESQ' )
                  IDTSESQ = HcoIDs(I)
                  WRITE(MSG,*) '     Sesquiterpenes = ',
     &               TRIM(SpcNames(I)),IDTSESQ
               CASE DEFAULT
                  MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
                  CALL HCO_ERROR( MSG, RC )
                  RETURN
            END SELECT

            ! Verbose
            IF ( am_I_Root ) CALL HCO_MSG( MSG )
         ENDDO

         ! Check if at least one of the SOA species is defined 
         IF ( IDTMTPA < 0 .AND. IDTMTPO < 0 .AND.  
     &        IDTLIMO < 0 .AND. IDTSESQ < 0        ) THEN
            MSG = 'MEGAN_SOA option is anabled but none of the '// 
     &            'SOA species is found (MTPA, MTPO, LIMO, SESQ)!'
            CALL HCO_ERROR( MSG, RC )
            RETURN
         ENDIF

         ! Never calculate OCPI if SOA is turned on!!
         IF ( IDTOCPI > 0 ) THEN
            IDTOCPI = -1
            IF ( am_I_Root ) THEN
               MSG = 'SOA option is turned on - '// 
     &               'no OCPI emissions will be calculated'
               CALL HCO_WARNING(MSG,RC)
            ENDIF
         ENDIF
      ENDIF

      !-----------------------------------------------------------------
      ! Allocate module variables 
      !-----------------------------------------------------------------

      ! Get horizontal grid extensions on this CPU
      NX = HcoState%NX
      NY = HcoState%NY

      ! Allocate arrays
      ALLOCATE( HCOT_DAY( NX, NY, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_DAY', RC )
        RETURN
      ENDIF
      HCOT_DAY = 0.0_sp

      ALLOCATE( HCOT_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_15', RC )
        RETURN
      ENDIF
      HCOT_15 = 0.0_sp
      T_PREVDAY => HCOT_15(:,:,1)

      ALLOCATE( HCOT_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_15_AVG', RC )
        RETURN
      ENDIF
      HCOT_15_AVG = 0.0_sp

!      ! Daily averaged temperature (mpb,2009)
!      ALLOCATE( HCOT_DAILY( NX, NY ), STAT=AS )
!      IF ( AS /= 0 ) THEN 
!        CALL HCO_ERROR( 'HCOT_DAILY', RC )
!        RETURN
!      ENDIF
!      HCOT_DAILY = 0.0_hp
 
      ! Allocate arrays for light (mpb,2009)
      ! -- Direct --
      ALLOCATE( HCOPARDR_DAY( NX, NY, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_DAY', RC )
        RETURN
      ENDIF
      HCOPARDR_DAY = 0.0_sp

      ALLOCATE( HCOPARDR_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_15', RC )
        RETURN
      ENDIF
      HCOPARDR_15 = 0.0_sp

      ALLOCATE( HCOPARDR_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_15_AVG', RC )
        RETURN
      ENDIF
      HCOPARDR_15_AVG = 0.0_sp

!      ALLOCATE( HCOPARDR_DAILY( NX, NY ), STAT=AS )
!      IF ( AS /= 0 ) THEN 
!        CALL HCO_ERROR( 'HCOPARDR_DAILY', RC )
!        RETURN
!      ENDIF
!      HCOPARDR_DAILY = 0.0_hp
 
      ! -- Diffuse --
      ALLOCATE( HCOPARDF_DAY( NX, NY, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_DAY', RC )
        RETURN
      ENDIF
      HCOPARDF_DAY = 0.0_sp

      ALLOCATE( HCOPARDF_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_15', RC )
        RETURN
      ENDIF
      HCOPARDF_15 = 0.0_sp

      ALLOCATE( HCOPARDF_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_15_AVG', RC )
        RETURN
      ENDIF
      HCOPARDF_15_AVG = 0.0_sp

!      ALLOCATE( HCOPARDF_DAILY( NX, NY ), STAT=AS )
!      IF ( AS /= 0 ) THEN 
!        CALL HCO_ERROR( 'HCOPARDF_DAILY', RC )
!        RETURN
!      ENDIF
!      HCOPARDF_DAILY = 0d0

      ! LAI from previous time step
      ALLOCATE( LAI_PREVSTEP( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'LAI_PREVSTEP', RC )
        RETURN
      ENDIF
      LAI_PREVSTEP = 0.0_sp

      ALLOCATE( LAI_PREVDAY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'LAI_PREVDAY', RC )
        RETURN
      ENDIF
      LAI_PREVDAY = 0.0_sp

      ! Normalization factor
      ! There should be a different normalization factor for each compound, but
      ! we calculate only 1 normalization factor for all compounds (dbm 11/2012)
      ALLOCATE( NORM_FAC( 1 ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'NORM_FAC', RC )
        RETURN
      ENDIF
      NORM_FAC = -99d0

      IF ( ExtNrMono > 0 ) THEN

         ALLOCATE( AEF_APIN( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_APIN', RC )
            RETURN
         ENDIF
         AEF_APIN = 0.0_hp

         ALLOCATE( AEF_MYRC( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_MYRC', RC )
            RETURN
         ENDIF
         AEF_MYRC = 0.0_hp

         ALLOCATE( AEF_OMON( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_OMON', RC )
            RETURN
         ENDIF
         AEF_OMON = 0.0_hp

         ALLOCATE( AEF_FARN( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_FARN', RC )
            RETURN
         ENDIF
         AEF_FARN = 0.0_hp

         ALLOCATE( AEF_BCAR( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_BCAR', RC )
            RETURN
         ENDIF
         AEF_BCAR = 0.0_hp

         ALLOCATE( AEF_OSQT( NX, NY ), STAT=AS )
         IF ( AS /= 0 ) THEN 
            CALL HCO_ERROR( 'AEF_OSQT', RC )
            RETURN
         ENDIF
         AEF_OSQT = 0.0_hp

      ENDIF

      ALLOCATE( AEF_MOHX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_MOHX', RC )
         RETURN
      ENDIF
      AEF_MOHX = 0.0_hp

      ALLOCATE( AEF_ACET( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_ACET', RC )
         RETURN
      ENDIF
      AEF_ACET = 0.0_hp

      ALLOCATE( AEF_ETOH( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_ETOH', RC )
         RETURN
      ENDIF
      AEF_ETOH = 0.0_hp

      ALLOCATE( AEF_CH2O( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_CH2O', RC )
         RETURN
      ENDIF
      AEF_CH2O = 0.0_hp

      ALLOCATE( AEF_ALD2( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_ALD2', RC )
         RETURN
      ENDIF
      AEF_ALD2 = 0.0_hp

      ALLOCATE( AEF_FAXX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_FAXX', RC )
         RETURN
      ENDIF
      AEF_FAXX = 0.0_hp

      ALLOCATE( AEF_AAXX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_AAXX', RC )
         RETURN
      ENDIF
      AEF_AAXX = 0.0_hp

      ALLOCATE( AEF_C2H4( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_C2H4', RC )
         RETURN
      ENDIF
      AEF_C2H4 = 0.0_hp

      ALLOCATE( AEF_TOLU( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_TOLU', RC )
         RETURN
      ENDIF
      AEF_TOLU = 0.0_hp

      ALLOCATE( AEF_HCNX( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_HCNX', RC )
         RETURN
      ENDIF
      AEF_HCNX = 0.0_hp

      ALLOCATE( AEF_PRPE( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
         CALL HCO_ERROR( 'AEF_PRPE', RC )
         RETURN
      ENDIF
      AEF_PRPE = 0.0_hp

      !=================================================================
      ! Initialize internal diagnostics. These are the restart variables
      ! that can be used for a 'warm' start of MEGAN.
      !=================================================================
      CALL HCO_RestartDefine ( am_I_Root, HcoState, 
     &                         'T_DAVG',  HCOT_15_AVG, 'K', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
   
      CALL HCO_RestartDefine ( am_I_Root,   HcoState, 
     &                         'T_PREVDAY', T_PREVDAY, 'K', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine ( am_I_Root,     HcoState, 
     &                         'LAI_PREVDAY', LAI_PREVDAY, '1', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_RestartDefine ( am_I_Root,    HcoState, 
     &                         'PARDR_DAVG', HCOPARDR_15_AVG, 
     &                         'W/m2',       RC )

      CALL HCO_RestartDefine ( am_I_Root,    HcoState, 
     &                         'PARDF_DAVG', HCOPARDF_15_AVG, 
     &                         'W/m2',       RC )

      !=================================================================
      ! The original MEGAN code used to read the emission factors here.
      ! We now get the emisson factors through HEMCO, hence no need to
      ! do this anymore (ckeller, 05/19/2014).
      !=================================================================

      ! Set physical constants
      PI_180 = HcoState%Phys%PI / 180.0_dp
      D2RAD  = PI_180                ! Degrees to radians
      RAD2D  = 1.0_hp / PI_180       ! Radians to degrees

      ! Enable met. fields
      ExtState%T2M%DoUse       = .TRUE.
      ExtState%SUNCOS%DoUse    = .TRUE.
      ExtState%PARDR%DoUse     = .TRUE.
      ExtState%PARDF%DoUse     = .TRUE.
      ExtState%LAI%DoUse       = .TRUE.
      ExtState%GWETROOT%DoUse  = .TRUE.

      ! SZAFACT field is only used for SOA simulation
      IF ( ExtNrSoa > 0 ) THEN
         ExtState%SZAFACT%DoUse = .TRUE.
      ENDIF

      ! Enable module 
      ExtState%Megan = .TRUE.

      ! Leave w/ success
      IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
      IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
      CALL HCO_LEAVE( RC )

      END SUBROUTINE HCOX_Megan_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Megan_Final
!
! !DESCRIPTION: Subroutine HCOX\_Megan\_Final deallocates all allocated arrays 
!  at the end of a GEOS-Chem model run.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_MEGAN_FINAL ( am_I_Root, HcoState, RC )
!
! !USES
!
!
! !INPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root     ! Root CPU?
      TYPE(HCO_State), POINTER        :: HcoState      ! HEMCO State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT)  :: RC 
! 
! !REVISION HISTORY: 
!  17 Dec 2009 - R. Yantosca  - Added ProTeX headers
!  16 Aug 2013 - C. Keller    - Now a HEMCO extension 
!  18 Feb 2015 - M. Sulprizio - Added AEF arrays
!  09 Mar 2015 - C. Keller    - Write variables to internal state in ESMF mode
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCOX_MEGAN_FINAL begins here!
      !=================================================================

      ! Nullify pointers
      GEIA_ORVC => NULL()
      AEF_ISOP  => NULL()
      AEF_MBOX  => NULL()
      AEF_BPIN  => NULL()
      AEF_CARE  => NULL()
      AEF_LIMO  => NULL()
      AEF_OCIM  => NULL()
      AEF_SABI  => NULL()

      ! Cleanup module arrays
      T_PREVDAY => NULL()
      IF ( ALLOCATED( HCOT_DAY        ) ) DEALLOCATE( HCOT_DAY        )
      IF ( ALLOCATED( HCOT_15         ) ) DEALLOCATE( HCOT_15         )
      IF ( ALLOCATED( HCOT_15_AVG     ) ) DEALLOCATE( HCOT_15_AVG     )
!      IF ( ALLOCATED( HCOT_DAILY      ) ) DEALLOCATE( HCOT_DAILY      )
      IF ( ALLOCATED( HCOPARDR_DAY    ) ) DEALLOCATE( HCOPARDR_DAY    )
      IF ( ALLOCATED( HCOPARDR_15     ) ) DEALLOCATE( HCOPARDR_15     )
      IF ( ALLOCATED( HCOPARDR_15_AVG ) ) DEALLOCATE( HCOPARDR_15_AVG )
!      IF ( ALLOCATED( HCOPARDR_DAILY  ) ) DEALLOCATE( HCOPARDR_DAILY  )
      IF ( ALLOCATED( HCOPARDF_DAY    ) ) DEALLOCATE( HCOPARDF_DAY    )
      IF ( ALLOCATED( HCOPARDF_15     ) ) DEALLOCATE( HCOPARDF_15     )
      IF ( ALLOCATED( HCOPARDF_15_AVG ) ) DEALLOCATE( HCOPARDF_15_AVG )
!      IF ( ALLOCATED( HCOPARDF_DAILY  ) ) DEALLOCATE( HCOPARDF_DAILY  )
      IF ( ALLOCATED( LAI_PREVSTEP    ) ) DEALLOCATE( LAI_PREVSTEP    )
      IF ( ALLOCATED( LAI_PREVDAY     ) ) DEALLOCATE( LAI_PREVDAY     )
      IF ( ALLOCATED( NORM_FAC        ) ) DEALLOCATE( NORM_FAC        )
      IF ( ALLOCATED( AEF_APIN        ) ) DEALLOCATE( AEF_APIN        )
      IF ( ALLOCATED( AEF_MYRC        ) ) DEALLOCATE( AEF_MYRC        )
      IF ( ALLOCATED( AEF_OMON        ) ) DEALLOCATE( AEF_OMON        )
      IF ( ALLOCATED( AEF_MOHX        ) ) DEALLOCATE( AEF_MOHX        )
      IF ( ALLOCATED( AEF_ACET        ) ) DEALLOCATE( AEF_ACET        )
      IF ( ALLOCATED( AEF_ETOH        ) ) DEALLOCATE( AEF_ETOH        )
      IF ( ALLOCATED( AEF_CH2O        ) ) DEALLOCATE( AEF_CH2O        )
      IF ( ALLOCATED( AEF_ALD2        ) ) DEALLOCATE( AEF_ALD2        )
      IF ( ALLOCATED( AEF_FAXX        ) ) DEALLOCATE( AEF_FAXX        )
      IF ( ALLOCATED( AEF_AAXX        ) ) DEALLOCATE( AEF_AAXX        )
      IF ( ALLOCATED( AEF_C2H4        ) ) DEALLOCATE( AEF_C2H4        )
      IF ( ALLOCATED( AEF_TOLU        ) ) DEALLOCATE( AEF_TOLU        )
      IF ( ALLOCATED( AEF_HCNX        ) ) DEALLOCATE( AEF_HCNX        )
      IF ( ALLOCATED( AEF_PRPE        ) ) DEALLOCATE( AEF_PRPE        )
      IF ( ALLOCATED( AEF_FARN        ) ) DEALLOCATE( AEF_FARN        )
      IF ( ALLOCATED( AEF_BCAR        ) ) DEALLOCATE( AEF_BCAR        )
      IF ( ALLOCATED( AEF_OSQT        ) ) DEALLOCATE( AEF_OSQT        )

      END SUBROUTINE HCOX_MEGAN_FINAL
!EOC
      END MODULE HCOX_MEGAN_MOD
