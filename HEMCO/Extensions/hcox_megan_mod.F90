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
! radiation information from the past. In the original GEOS-Chem
! code, the initial 10-d averages were explicitly calculated during
! initialization of MEGAN. This is not feasible in an ESMF environment,
! and the following restart variables can now be provided through the
! HEMCO configuration file:
! \begin{itemize}
! \item T\_DAVG: long-term historical temperature
! \item PARDR\_DAVG: long-term historical direct radiation
! \item PARDF\_DAVG: long-term historical diffuse radiation
! \item T\_PREVDAY: short-term historical temperature
! \end{itemize}
! These variables are automatically searched for on the first call of
! the run call. If not defined, default values will be used. The values
! of T\_DAVG, T\_PREVDAY, PARDR\_DAVG, and PARDF\_DAVG are continuously
! updated at the end of the run sequence, e.g. they represent the
! instantaneous running average. The e-folding times to be used when
! calculating the short=term and long-term running averages are defined
! as module parameter below (parameter TAU\_HOURS and TAU\_DAYS).
!\\
!\\
! A similar procedure is also applied to the leaf area index variables.
! The original GEOS-Chem MEGAN code used three LAI variables: current month
! LAI (LAI\_CM), previous month LAI (LAI\_PM), and instantaneous LAI (LAI),
! which was a daily interpolation of LAI\_CM and next month' LAI, (LAI\_NM).
! The HEMCO implementation uses only the instantaneous LAI, assuming it is
! updated every day. The short term historical LAI is kept in memory and
! used to determine the LAI change over time (used to calculate the gamma
! leaf age). It is also updated on every time step.
! For the first simulation day, the previous' day LAI is taken from the
! restart file (field LAI\_PREVDAY). If no restart variable is defined, a
! LAI change of zero is assumed (ckeller, 10/9/2014).
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
  PRIVATE :: GET_GAMMA_CO2  ! (Tai, Jan 2013)
!
! !REVISION HISTORY:
!  (1 ) Original code (biogen_em_mod.f) by Dorian Abbot (6/2003).  Updated to
!        latest algorithm and modified for the standard code by May Fu
!        (11/2004).
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  TYPE :: MyInst
     INTEGER                 :: Instance
     INTEGER                 :: ExtNr        ! Extension number for MEGAN
     REAL(hp)                :: ISOP_SCALING ! Isoprene emissions scaling
     REAL(hp)                :: GLOBCO2      ! Global CO2 conc (ppmv)
     REAL(hp)                :: ISOPTOSOAP   ! SOAP Emission factor
     REAL(hp)                :: MONOTOSOAP   ! SOAP Emission factor
     REAL(hp)                :: OTHRTOSOAP   ! SOAP Emission factor
     REAL(hp)                :: ISOPTOSOAS   ! Direct SOA Emission factor
     REAL(hp)                :: MONOTOSOAS   ! Direct SOA Emission factor
     REAL(hp)                :: OTHRTOSOAS   ! Direct SOA Emission factor
     LOGICAL                 :: LISOPCO2     ! Include CO2 inhibition of ISOP?
     LOGICAL                 :: NORMLAI      ! Normalize LAI by PFT?
     LOGICAL                 :: OFFLINE_BIOGENICVOC ! Use offline emiss?

     ! Physical parameter
     REAL(hp)                :: D2RAD  ! Degrees to radians
     REAL(hp)                :: RAD2D  ! Radians to degrees

     ! HEMCO species IDs
     INTEGER                 :: IDTISOP
     INTEGER                 :: IDTACET
     INTEGER                 :: IDTPRPE
     INTEGER                 :: IDTC2H4
     INTEGER                 :: IDTALD2
     INTEGER                 :: IDTMOH
     INTEGER                 :: IDTEOH
     INTEGER                 :: IDTLIMO
     INTEGER                 :: IDTMTPA
     INTEGER                 :: IDTMTPO
     INTEGER                 :: IDTSESQ
     INTEGER                 :: IDTSOAP
     INTEGER                 :: IDTSOAS

     ! Pointers to annual emission factor arrays.
     ! These fields are obtained from ext. data (from config file)
     ! These arrays represent the base emission values in kg(C)/m2/s
     ! that will be scaled based upon meteorological conditions.
     REAL(hp), POINTER       :: AEF_ISOP(:,:)
     REAL(hp), POINTER       :: AEF_MBOX(:,:)
     REAL(hp), POINTER       :: AEF_BPIN(:,:)
     REAL(hp), POINTER       :: AEF_CARE(:,:)
     REAL(hp), POINTER       :: AEF_LIMO(:,:)
     REAL(hp), POINTER       :: AEF_OCIM(:,:)
     REAL(hp), POINTER       :: AEF_SABI(:,:)
     REAL(hp), POINTER       :: GEIA_ORVC(:,:)

     ! Annual emission factors arrays
     ! These fields are not read from file, but are computed in CALC_AEF
     REAL(hp), POINTER       :: AEF_APIN(:,:)              ! Alpha-pinene
     REAL(hp), POINTER       :: AEF_MYRC(:,:)              ! Myrcene
     REAL(hp), POINTER       :: AEF_OMON(:,:)              ! Other monoterpenes
     REAL(hp), POINTER       :: AEF_ACET(:,:)              ! Acetone
     REAL(hp), POINTER       :: AEF_MOH(:,:)               ! Methanol
     REAL(hp), POINTER       :: AEF_EOH(:,:)               ! Ethanol
     REAL(hp), POINTER       :: AEF_CH2O(:,:)              ! Formaldehyde
     REAL(hp), POINTER       :: AEF_ALD2(:,:)              ! Acetaldehyde
     REAL(hp), POINTER       :: AEF_FAXX(:,:)              ! Formic acid
     REAL(hp), POINTER       :: AEF_AAXX(:,:)              ! Acetic acid
     REAL(hp), POINTER       :: AEF_C2H4(:,:)              ! Ethene
     REAL(hp), POINTER       :: AEF_TOLU(:,:)              ! Toluene
     REAL(hp), POINTER       :: AEF_HCNX(:,:)              ! HCN
     REAL(hp), POINTER       :: AEF_PRPE(:,:)              ! >= C3 alkenes
     REAL(hp), POINTER       :: AEF_FARN(:,:)              ! a-Farnesene
     REAL(hp), POINTER       :: AEF_BCAR(:,:)              ! b-Caryophyllene
     REAL(hp), POINTER       :: AEF_OSQT(:,:)              ! Other sesquiterp.

     ! Emission arrays (kgC/m2/s)
     REAL(hp), POINTER       :: FLUXISOP(:,:)
     REAL(hp), POINTER       :: FLUXMONO(:,:)
     REAL(hp), POINTER       :: FLUXACETmo(:,:)
     REAL(hp), POINTER       :: FLUXACETmb(:,:)
     REAL(hp), POINTER       :: FLUXACETbg(:,:)
     REAL(hp), POINTER       :: FLUXPRPE(:,:)
     REAL(hp), POINTER       :: FLUXC2H4(:,:)
     REAL(hp), POINTER       :: FLUXLIMO(:,:)
     REAL(hp), POINTER       :: FLUXMTPA(:,:)
     REAL(hp), POINTER       :: FLUXMTPO(:,:)
     REAL(hp), POINTER       :: FLUXSESQ(:,:)
     REAL(hp), POINTER       :: FLUXSOAP(:,:)
     REAL(hp), POINTER       :: FLUXSOAS(:,:)
     REAL(hp), POINTER       :: FLUXALD2(:,:)
     REAL(hp), POINTER       :: FLUXMOH(:,:)
     REAL(hp), POINTER       :: FLUXEOH(:,:)
         
     ! Emission arrays for use in diagnostics only (kgC/m2/s)
     REAL(hp), POINTER       :: FLUXAPIN(:,:)
     REAL(hp), POINTER       :: FLUXBPIN(:,:)
     REAL(hp), POINTER       :: FLUXSABI(:,:)
     REAL(hp), POINTER       :: FLUXMYRC(:,:)
     REAL(hp), POINTER       :: FLUXCARE(:,:)
     REAL(hp), POINTER       :: FLUXOCIM(:,:)
     REAL(hp), POINTER       :: FLUXOMON(:,:)
     REAL(hp), POINTER       :: FLUXFARN(:,:)
     REAL(hp), POINTER       :: FLUXBCAR(:,:)
     REAL(hp), POINTER       :: FLUXOSQT(:,:)
     REAL(hp), POINTER       :: FLUXMBOX(:,:)
     REAL(hp), POINTER       :: FLUXFAXX(:,:)
     REAL(hp), POINTER       :: FLUXAAXX(:,:)

     ! Normalization factor
     REAL(hp), POINTER       :: NORM_FAC(:)

     ! New restart variables (ckeller, 11/05/2015)
     REAL(sp), POINTER       :: T_LASTXDAYS    (:,:) ! Avg. temp  of last X days
     REAL(sp), POINTER       :: T_LAST24H      (:,:) ! Avg. temp  of last 24 hrs
     REAL(sp), POINTER       :: PARDF_LASTXDAYS(:,:) ! Avg. PARDF of last X days
     REAL(sp), POINTER       :: PARDR_LASTXDAYS(:,:) ! Avg. PARDR of last X days

     ! previous LAI values (ckeller, 10/09/2014)
     REAL(sp), POINTER       :: LAI_PREVDAY(:,:)     ! LAI of prev. day

     ! Array for PFT
     REAL(hp), POINTER       :: ARRAY_16(:,:,:)

     ! Days between mid-months (updated by HEMCO clock)
     INTEGER                 :: DAYS_BTW_M

     ! Input data suffix
     CHARACTER(LEN=255)      :: SUFFIX

     TYPE(MyInst), POINTER   :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER      :: AllInst => NULL()

  ! Scalars
  ! e-folding time to be applied to long-term past conditions (in days)
  REAL(hp), PARAMETER  :: TAU_DAYS  = 5.0_hp

  ! e-folding time to be applied to short-term past conditions (in hours)
  REAL(hp), PARAMETER  :: TAU_HOURS = 12.0_hp

  ! W/m2 -> umol/m2/s
  REAL(hp), PARAMETER  :: WM2_TO_UMOLM2S = 4.766_hp

  ! Maximum LAI value [cm2/cm2]
  REAL(hp), PARAMETER  :: LAI_MAX = 6.0_hp

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
  SUBROUTINE HCOX_Megan_Run( ExtState, HcoState, RC )
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
    TYPE(Ext_State), POINTER       :: ExtState
    TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  05 Aug 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J
    REAL(hp)            :: EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN
    REAL(hp)            :: EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE
    REAL(hp)            :: EMIS_OCIM, EMIS_OMON, EMIS_MONO, EMIS_ALD2
    REAL(hp)            :: EMIS_MOH,  EMIS_EOH,  EMIS_FAXX, EMIS_AAXX
    REAL(hp)            :: EMIS_ACET, EMIS_PRPE, EMIS_C2H4
    REAL(hp)            :: EMIS_FARN, EMIS_BCAR, EMIS_OSQT
    REAL(hp)            :: EMIS_OTHR
    REAL(hp)            :: BIOG_ALPH, BIOG_LIMO
    REAL(hp)            :: BIOG_ALCO
    REAL(hp)            :: MONO_MOL
    REAL(hp)            :: MBOX_MOL
    REAL(hp)            :: X, Y
    REAL(hp)            :: DIURFACT, DIURORVC
    REAL(hp)            :: TMP
    REAL(hp)            :: TS_EMIS
    REAL(hp)            :: DNEWFRAC
    REAL(hp)            :: DOLDFRAC
    REAL(hp)            :: HNEWFRAC
    REAL(hp)            :: HOLDFRAC
    LOGICAL             :: ERR
    LOGICAL             :: FIRST, IsNewHour, IsNewDay

    ! For diagnostics
    REAL(hp), POINTER   :: Arr2D(:,:)
    CHARACTER(LEN=63)   :: DiagnName
    CHARACTER(LEN=255)  :: MSG

    ! Molecules C / kg C
    REAL(hp), PARAMETER   :: XNUMOL_C = 6.022140857e+23_hp / 12e-3_hp
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

    TYPE(MyInst), POINTER :: Inst

    !=================================================================
    ! HCOX_Megan_Run begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, &
                    'HCOX_Megan_Run (hcox_megan_mod.F)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    ERR = .FALSE.

    ! Nullify
    Arr2d => NULL()
    Inst  => NULL()

    ! Get instance
    CALL InstGet ( ExtState%Megan, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find Megan instance Nr. ', ExtState%Megan
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! Time information
    FIRST     = HcoClock_First  ( HcoState%Clock, .TRUE. )
    IsNewHour = HcoClock_NewHour( HcoState%Clock, .TRUE. )
    IsNewDay  = HcoClock_NewDay ( HcoState%Clock, .TRUE. )

    ! Initialize flux arrays
    Inst%FLUXISOP   = 0.0_hp
    Inst%FLUXMBOX   = 0.0_hp
    Inst%FLUXAPIN   = 0.0_hp
    Inst%FLUXBPIN   = 0.0_hp
    Inst%FLUXSABI   = 0.0_hp
    Inst%FLUXMYRC   = 0.0_hp
    Inst%FLUXCARE   = 0.0_hp
    Inst%FLUXOCIM   = 0.0_hp
    Inst%FLUXOMON   = 0.0_hp
    Inst%FLUXMONO   = 0.0_hp
    Inst%FLUXALD2   = 0.0_hp
    Inst%FLUXMOH    = 0.0_hp
    Inst%FLUXEOH    = 0.0_hp
    Inst%FLUXFAXX   = 0.0_hp
    Inst%FLUXAAXX   = 0.0_hp
    Inst%FLUXACETmo = 0.0_hp
    Inst%FLUXACETmb = 0.0_hp
    Inst%FLUXACETbg = 0.0_hp
    Inst%FLUXPRPE   = 0.0_hp
    Inst%FLUXC2H4   = 0.0_hp
    Inst%FLUXLIMO   = 0.0_hp
    Inst%FLUXMTPA   = 0.0_hp
    Inst%FLUXMTPO   = 0.0_hp
    Inst%FLUXFARN   = 0.0_hp
    Inst%FLUXBCAR   = 0.0_hp
    Inst%FLUXOSQT   = 0.0_hp
    Inst%FLUXSESQ   = 0.0_hp
    Inst%FLUXSOAP   = 0.0_hp
    Inst%FLUXSOAS   = 0.0_hp

#if defined( ESMF_ ) || defined( MODEL_GEOS )

    !----------------------------------------------------------------
    ! %%%%% ESMF environment: execute these on every call  %%%%%
    ! %%%%% because this will fill from the External State %%%%%
    !----------------------------------------------------------------

    ! Generate annual emission factors for MEGAN inventory
    CALL CALC_AEF( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in MEGAN routine CALC_AEF!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Calculate normalization factor (dbm, 11/2012)
    CALL CALC_NORM_FAC( Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in MEGAN routine CALC_NORM_FAC!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! Fill restart variables.
    CALL FILL_RESTART_VARS( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Error encountered in MEGAN routine FILL_RESTART_VARS!'
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

#else

    !----------------------------------------------------------------
    ! %%%%% Standard environment: Execute these only once     %%%%%
    ! %%%%% to avoid constantly overwriting restart variables %%%%%
    !----------------------------------------------------------------
    IF ( FIRST ) THEN

       ! Generate annual emission factors for MEGAN inventory
       CALL CALC_AEF( HcoState, ExtState, Inst, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Error encountered in MEGAN routine CALC_AEF!'
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Calculate normalization factor (dbm, 11/2012)
       CALL CALC_NORM_FAC( Inst, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Error encountered in MEGAN routine CALC_NORM_FAC!'
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Fill restart variables
       CALL FILL_RESTART_VARS( HcoState, ExtState, Inst, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Error encountered in MEGAN routine FILL_RESTART_VARS!'
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF
    ENDIF

#endif

    ! Set to 1 day since we know that LAI is updated every day.
    Inst%DAYS_BTW_M = 1

    ! Calculate weights for running means of historic variables
    ! DNEWFRAC and DOLDFRAC are the weights given to the current
    ! and existing value, respectively, when updating running means
    ! over the last X days. HNEWFRAC and HOLDFRAC are the same but
    ! for the 24H means. (ckeller, 11/05/2015)
    TS_EMIS  = HcoState%TS_EMIS
    DNEWFRAC = TS_EMIS / ( TAU_DAYS * 24.0_hp * 3600.0_hp )
    DOLDFRAC = 1.0_hp - DNEWFRAC
    HNEWFRAC = TS_EMIS / ( TAU_HOURS * 3600.0_hp )
    HOLDFRAC = 1.0_hp - HNEWFRAC

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J,      EMIS_ISOP, EMIS_MBOX, EMIS_APIN, EMIS_BPIN   ) &
    !$OMP PRIVATE( EMIS_LIMO, EMIS_SABI, EMIS_MYRC, EMIS_CARE, EMIS_OCIM   ) &
    !$OMP PRIVATE( EMIS_OMON, EMIS_MONO, EMIS_MOH,  EMIS_EOH,  EMIS_FAXX   ) &
    !$OMP PRIVATE( EMIS_AAXX, EMIS_ACET, EMIS_PRPE, EMIS_C2H4              ) &
    !$OMP PRIVATE( EMIS_FARN, EMIS_BCAR, EMIS_OSQT, EMIS_ALD2, EMIS_OTHR   ) &
    !$OMP PRIVATE( MONO_MOL,  TMP,       MBOX_MOL                          ) &
    !$OMP PRIVATE( DIURFACT,  DIURORVC,  BIOG_ALPH                         ) &
    !$OMP PRIVATE( BIOG_LIMO, BIOG_ALCO, X, Y, RC                          )

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
       EMIS_OTHR  = 0.0_hp
       EMIS_ALD2  = 0.0_hp
       EMIS_MOH   = 0.0_hp
       EMIS_EOH   = 0.0_hp
       EMIS_FAXX  = 0.0_hp
       EMIS_AAXX  = 0.0_hp
       EMIS_ACET  = 0.0_hp
       EMIS_PRPE  = 0.0_hp
       EMIS_C2H4  = 0.0_hp
       EMIS_FARN  = 0.0_hp
       EMIS_BCAR  = 0.0_hp
       EMIS_OSQT  = 0.0_hp

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
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, Inst, &
                                 I, J, 'ISOP', EMIS_ISOP, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_ISOP', RC )
          ERR = .TRUE.
          EXIT
       ENDIF

       !FP_ISOP (6/2009)
       EMIS_ISOP = Inst%ISOP_SCALING * EMIS_ISOP

       ! Write isoprene emissions into tracer tendency array
       IF ( Inst%IDTISOP > 0 ) THEN
          Inst%FLUXISOP(I,J) = EMIS_ISOP
       ENDIF

       ! Biogenic emissions of SOA and SOA-Precursor from isopene
       IF ( (Inst%IDTISOP>0) .AND. (Inst%IDTSOAP>0) ) THEN
          Inst%FLUXSOAP(I,J) = Inst%FLUXSOAP(I,J) + EMIS_ISOP * Inst%ISOPTOSOAP
       ENDIF
       IF ( (Inst%IDTISOP>0) .AND. (Inst%IDTSOAS>0) ) THEN
          Inst%FLUXSOAS(I,J) = Inst%FLUXSOAS(I,J) + EMIS_ISOP * Inst%ISOPTOSOAS
       ENDIF

       !--------------------------------------------------------------------
       ! MEGAN monoterpenes
       !--------------------------------------------------------------------

       ! ---------------------------------------------------
       ! Alpha Pinene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'APIN', EMIS_APIN, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_APIN', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXAPIN(I,J) = EMIS_APIN

       ! ---------------------------------------------------
       ! Beta Pinene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'BPIN', EMIS_BPIN, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_BPIN', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXBPIN(I,J) = EMIS_BPIN

       ! ---------------------------------------------------
       ! Limonene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'LIMO', EMIS_LIMO, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_LIMO', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXLIMO(I,J) = EMIS_LIMO

       ! ---------------------------------------------------
       ! Sabinene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'SABI', EMIS_SABI, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_SABI', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXSABI(I,J) = EMIS_SABI

       ! ---------------------------------------------------
       ! Mycrene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'MYRC', EMIS_MYRC, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_MYRC', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXMYRC(I,J) = EMIS_MYRC

       ! ---------------------------------------------------
       ! 3-Carene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'CARE', EMIS_CARE, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_CARE', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXCARE(I,J) = EMIS_CARE

       ! ---------------------------------------------------
       ! Ocimene emissions
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'OCIM', EMIS_OCIM, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_OCIM', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXOCIM(I,J) = EMIS_OCIM

       ! ---------------------------------------------------
       ! Other monoterpenes
       ! (added 12/2012; dbm)
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'OMON', EMIS_OMON, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_OMON', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXOMON(I,J) = EMIS_OMON

       ! ---------------------------------------------------
       ! Total monoterpenes = sum of individual
       ! dbm, now add other lumped monoterpenes (11/2012)
       EMIS_MONO = EMIS_APIN + EMIS_BPIN + EMIS_LIMO + EMIS_SABI + &
                   EMIS_MYRC + EMIS_CARE + EMIS_OCIM + EMIS_OMON

       ! Add to tracer tendency array [kg C/m2/s]
       Inst%FLUXMONO(I,J) = EMIS_MONO

       !! testing only
       !if ( i == ix .AND. j == iy .and. .FALSE. ) then
       !        write(*,*) ' '
       !        write(*,*) 'HEMCO MEGAN (kg/m2/s) @ ', i, j
       !        write(*,*) 'ISOP: ', EMIS_ISOP
       !        write(*,*) 'MONO: ', EMIS_MONO
       !        write(*,*) 'APIN: ', EMIS_APIN
       !        write(*,*) 'BPIN: ', EMIS_BPIN
       !        write(*,*) 'LIMO: ', EMIS_LIMO
       !        write(*,*) 'SABI: ', EMIS_SABI
       !        write(*,*) 'MYRC: ', EMIS_MYRC
       !        write(*,*) 'CARE: ', EMIS_CARE
       !        write(*,*) 'OCIM: ', EMIS_OCIM
       !        write(*,*) 'OMON: ', EMIS_OMON
       !endif

       !--------------------------------------------------------------------
       ! MEGAN Acetaldehyde
       !--------------------------------------------------------------------
       IF ( Inst%IDTALD2 > 0 ) THEN
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'ALD2', EMIS_ALD2, RC)
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'GET_MEGAN_EMISSIONS_ALD2', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Add to tracer tendency array [kg C/m2/s]
          Inst%FLUXALD2(I,J) = EMIS_ALD2
       ENDIF

       ! ---------------------------------------------------
       ! MEGAN Methanol
       IF ( Inst%IDTMOH > 0 ) THEN
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'MOH', EMIS_MOH, RC)
          IF ( RC /= HCO_SUCCESS ) THEN 
             CALL HCO_ERROR( HcoState%Config%Err,  &
                             'GET_MEGAN_EMISSIONS_MOH', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Add to tracer tendency array [kg/m2/s]
          Inst%FLUXMOH(I,J) = EMIS_MOH
       ENDIF

       ! ---------------------------------------------------
       ! MEGAN Ethanol
       IF ( Inst%IDTEOH > 0 ) THEN
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'EOH', EMIS_EOH, RC)
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'GET_MEGAN_EMISSIONS_EOH', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Add to tracer tendency array [kg C/m2/s]
          Inst%FLUXEOH(I,J) = EMIS_EOH
       ENDIF

       !--------------------------------------------------------------------
       ! Other MEGAN biogenics
       ! Calls included here for future incorporation or
       ! specialized simulations. (dbm, 12/2012)
       !--------------------------------------------------------------------

       ! ---------------------------------------------------
       ! Methyl butenol
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'MBOX', EMIS_MBOX, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_MBOX', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXMBOX(I,J) = EMIS_MBOX

       ! ---------------------------------------------------
       ! MEGAN Formic acid
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'FAXX', EMIS_FAXX, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_FAXX', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXFAXX(I,J) = EMIS_FAXX

       ! ---------------------------------------------------
       ! MEGAN Acetic acid
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'AAXX', EMIS_AAXX, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_AAXX', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXAAXX(I,J) = EMIS_AAXX

       !--------------------------------------------------------------------
       ! Biogenic acetone emissions
       !--------------------------------------------------------------------
       IF ( Inst%IDTACET > 0 ) THEN

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
          Inst%FLUXACETmo(I,J) = TMP

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
          IF ( ( X >= -167.5_hp .and. X <= -52.5_hp ) .AND. &
               ( Y >=   16.0_hp .and. Y <=  72.0_hp ) ) THEN

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
             Inst%FLUXACETmb(I,J) = TMP * MB_SCALE2

          ENDIF

          !-----------------------------------------------------------------
          ! (3) BIOGENIC ACETONE -- DIRECT EMISSION
          ! evf, removed obsolete code, replaced with MEGAN acetone
          ! emissions (5/25/2011) Direct Emission now includes emission
          ! from grasses and emission from dry leaf matter
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'ACET', EMIS_ACET, RC)
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'GET_MEGAN_EMISSIONS_ACET', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Write to array
          Inst%FLUXACETbg(I,J) = EMIS_ACET

       ENDIF

       !--------------------------------------------------------------------
       ! Biogenic emissions of SOA and SOA-Precursor from monoterpenes
       !--------------------------------------------------------------------
       IF ( Inst%IDTSOAP>0 ) THEN
          Inst%FLUXSOAP(I,J) = Inst%FLUXSOAP(I,J) + EMIS_MONO * Inst%MONOTOSOAP
       ENDIF
       IF ( Inst%IDTSOAS>0 ) THEN
          Inst%FLUXSOAS(I,J) = Inst%FLUXSOAS(I,J) + EMIS_MONO * Inst%MONOTOSOAS
       ENDIF

       !--------------------------------------------------------------------
       ! Biogenic emissions of PRPE
       !
       ! Now uses MEGAN2.1 (dbm, 12/2012)
       !--------------------------------------------------------------------
       IF ( Inst%IDTPRPE > 0 ) THEN
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'PRPE', EMIS_PRPE, RC)
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'GET_MEGAN_EMISSIONS_PRPE', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Add to tracer tendency array [kg C/m2/s]
          Inst%FLUXPRPE(I,J) = EMIS_PRPE

          !! testing only
          !if ( i==ix .and. j==iy .and. .FALSE. ) then
          !   write(*,*) ' '
          !   write(*,*) 'HEMCO MEGAN @ ', ix, iy
          !   write(*,*) 'PRPE (kg/m2/s): ', Inst%FLUXPRPE(I,J)
          !endif

       ENDIF

       !--------------------------------------------------------------------
       ! Biogenic emissions of ethene (C2H4)
       !
       ! Now uses MEGAN2.1 (dbm, 12/2012)
       !--------------------------------------------------------------------
       IF ( Inst%IDTC2H4 > 0 ) THEN
          CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                    Inst, I, J, 'C2H4', EMIS_C2H4, RC)
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( HcoState%Config%Err, &
                             'GET_MEGAN_EMISSIONS_C2H4', RC )
             ERR = .TRUE.
             EXIT
          ENDIF

          ! Add to tracer tendency array [kg C/m2/s]
          Inst%FLUXC2H4(I,J) = EMIS_C2H4
       ENDIF

       ! ----------------------------------------------------------------
       ! The new MEGAN implementation has speciated information
       ! (hotp 3/7/10)
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
       !
       ! Monoterpene lumping:
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
       IF ( Inst%IDTMTPA > 0 ) THEN
          Inst%FLUXMTPA(I,J) = ( EMIS_APIN + EMIS_BPIN + &
                                 EMIS_SABI + EMIS_CARE ) * FC1
       ENDIF

       !--------------------------------------------------------------
       ! MEGAN Limonene
       !--------------------------------------------------------------
       ! [kg C/m2/s]
       IF ( Inst%IDTLIMO > 0 ) THEN
          Inst%FLUXLIMO(I,J) = EMIS_LIMO * FC1
       ENDIF

       !--------------------------------------------------------------
       ! MEGAN MTPO
       !--------------------------------------------------------------
       ! MTPO is all other monoterpenes (MEGAN categories:
       ! myrcene, ocimene, OMON) (hotp 5/20/10)
       ! All other monoterpenes (mostly camphene, linalool,
       ! terpinolene, terpinolene, phellandrene) (hotp 3/10/10)
       ! 14-18% of OMTP is terpinene and terpinolene
       IF ( Inst%IDTMTPO > 0 ) THEN
          Inst%FLUXMTPO(I,J) = ( EMIS_MYRC + EMIS_OCIM + EMIS_OMON ) * FC1
       ENDIF

       !--------------------------------------------------------------
       ! MEGAN sesquiterpenes
       !--------------------------------------------------------------

       ! ---------------------------------------------------
       ! a-Farnesene
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'FARN', EMIS_FARN, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_FARN', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXFARN(I,J) = EMIS_FARN

       ! ---------------------------------------------------
       ! b_Caryophyllene
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'BCAR', EMIS_BCAR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_BCAR', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXBCAR(I,J) = EMIS_BCAR

       ! ---------------------------------------------------
       ! Other sesquiterpene
       CALL GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                 Inst, I, J, 'OSQT', EMIS_OSQT, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'GET_MEGAN_EMISSIONS_OSQT', RC )
          ERR = .TRUE.
          EXIT
       ENDIF
       Inst%FLUXOSQT(I,J) = EMIS_OSQT

       ! ---------------------------------------------------
       ! Total sesquiterpenes from MEGAN (hotp 3/10/10)
       IF ( Inst%IDTSESQ > 0 ) THEN
          Inst%FLUXSESQ(I,J) = ( EMIS_FARN + EMIS_BCAR + EMIS_OSQT ) * FC3
       ENDIF

       ! Other terpenes
       EMIS_OTHR = EMIS_FARN + EMIS_BCAR + EMIS_OSQT

       !--------------------------------------------------------------------
       ! Biogenic emissions of SOA and SOA-Precursor from Other terpenes
       !--------------------------------------------------------------------
       IF ( Inst%IDTSOAP>0 ) THEN
          Inst%FLUXSOAP(I,J) = Inst%FLUXSOAP(I,J) + EMIS_OTHR * Inst%OTHRTOSOAP
       ENDIF
       IF ( Inst%IDTSOAS>0 ) THEN
          Inst%FLUXSOAS(I,J) = Inst%FLUXSOAS(I,J) + EMIS_OTHR * Inst%OTHRTOSOAS
       ENDIF

       !-----------------------------------------------------------------
       ! Update historical temperature / radiation values
       !-----------------------------------------------------------------
       ! Do this now on every time step. The arrays are simply the
       ! the running means over the intentend time window (24 hours
       ! and NUM_DAYS, respectively). This hugely faciliates warm
       ! model restarts, irrespective of simulation start/end dates
       ! and times. It also makes sure that all environmental
       ! variables are incorporated into the time averages (e.g. if
       ! emission time step is less than 60 minutes, all values will
       ! be used to calculate the daily mean).
       ! (ckeller, 11/05/2015)

       ! Updated LAI of last 24 hours
       Inst%LAI_PREVDAY(I,J) = ( HOLDFRAC * Inst%LAI_PREVDAY(I,J) ) + &
                               ( HNEWFRAC * ExtState%LAI%Arr%Val(I,J) )

       ! Updated temperature of last 24 hours
       Inst%T_LAST24H(I,J) = ( HOLDFRAC * Inst%T_LAST24H(I,J) ) + &
                             ( HNEWFRAC * ExtState%T2M%Arr%Val(I,J) )

       ! Updated temperature of last NUM_DAYS
       Inst%T_LASTXDAYS(I,J) = ( DOLDFRAC * Inst%T_LASTXDAYS(I,J) ) + &
                               ( DNEWFRAC * ExtState%T2M%Arr%Val(I,J) )

       ! Updated direct radiation of last NUM_DAYS
       Inst%PARDR_LASTXDAYS(I,J) = ( DOLDFRAC * Inst%PARDR_LASTXDAYS(I,J) ) + &
                                   ( DNEWFRAC * ExtState%PARDR%Arr%Val(I,J) )

       ! Updated diffuse radiation of last NUM_DAYS
       Inst%PARDF_LASTXDAYS(I,J) = ( DOLDFRAC * Inst%PARDF_LASTXDAYS(I,J) ) + &
                                   ( DNEWFRAC * ExtState%PARDF%Arr%Val(I,J) )

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

    !! testing only
    !if( HcoState%amIRoot) then
    ! write(*,*) 'ISOP emissions instance ',ExtState%Megan
    ! write(*,*) 'A MGN ',ExtState%Megan,Inst%IDTISOP,SUM(Inst%FLUXISOP)
    ! write(*,*) 'B MGN ',ExtState%Megan,SUM(Inst%LAI_PREVDAY)
    ! write(*,*) 'C MGN ',ExtState%Megan,SUM(Inst%T_LAST24H)
    ! write(*,*) 'D MGN ',ExtState%Megan,SUM(Inst%T_LASTXDAYS)
    ! write(*,*) 'E MGN ',ExtState%Megan,SUM(Inst%PARDR_LASTXDAYS)
    ! write(*,*) 'F MGN ',ExtState%Megan,SUM(Inst%PARDF_LASTXDAYS)
    ! write(*,*) 'G MGN ',ExtState%Megan,SUM(ExtState%T2M%Arr%Val)
    ! write(*,*) 'H MGN ',ExtState%Megan,SUM(ExtState%SUNCOS%Arr%Val)
    ! write(*,*) 'I MGN ',ExtState%Megan,SUM(ExtState%PARDR%Arr%Val)
    ! write(*,*) 'J MGN ',ExtState%Megan,SUM(ExtState%PARDF%Arr%Val)
    ! write(*,*) 'K MGN ',ExtState%Megan,SUM(ExtState%LAI%Arr%Val)
    ! write(*,*) 'L MGN ',ExtState%Megan,SUM(ExtState%GWETROOT%Arr%Val)
    !endif

    ! ----------------------------------------------------------------
    ! ISOPRENE
    IF ( Inst%IDTISOP > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXISOP, Inst%IDTISOP, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXISOP', RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! ACETALDEHYDE
    IF ( Inst%IDTALD2 > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXALD2, Inst%IDTALD2, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXALD2', RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! METHANOL
    IF ( Inst%IDTMOH > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXMOH, Inst%IDTMOH, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXMOH', RC )
          RETURN 
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! ETHANOL
    IF ( Inst%IDTEOH > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXEOH, Inst%IDTEOH, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXEOH', RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! ACETONE
    IF ( Inst%IDTACET > 0 ) THEN

       ! Eventually add individual diagnostics. These are assumed to
       ! have names MEGAN_ACET_MONO, MEGAN_ACET_MBO, and
       ! MEGAN_ACET_DIRECT
       DiagnName =  'InvMEGAN_ACET_MONO'
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Inst%FLUXACETmo, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN

       DiagnName =  'InvMEGAN_ACET_MBOX'
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Inst%FLUXACETmb, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN

       DiagnName =  'InvMEGAN_ACET_DIRECT'
       CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                          cName=TRIM(DiagnName), Array2D=Inst%FLUXACETbg, RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Add flux to emission array
       Inst%FLUXACETbg = Inst%FLUXACETbg + Inst%FLUXACETmo + Inst%FLUXACETmb
       CALL HCO_EmisAdd( HcoState, Inst%FLUXACETbg, Inst%IDTACET, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXACETbg', RC )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! SOA-Precursor (SOAP)
    IF (  Inst%IDTSOAP > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXSOAP, Inst%IDTSOAP, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXSOAP', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! SOA-Simplified (SOAS)
    IF (  Inst%IDTSOAS > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXSOAS, Inst%IDTSOAS, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXSOAS', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! ALKENES
    IF ( Inst%IDTPRPE > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXPRPE, Inst%IDTPRPE, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXPRPE', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! ETHENE
    IF ( Inst%IDTC2H4 > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXC2H4, Inst%IDTC2H4, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXC2H4', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! MTPA
    IF ( Inst%IDTMTPA > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXMTPA, Inst%IDTMTPA, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXMTPA', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! MTPO
    IF ( Inst%IDTMTPO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXMTPO, Inst%IDTMTPO, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXMTPO', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! LIMONENE
    IF ( Inst%IDTLIMO > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXLIMO, Inst%IDTLIMO, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXLIMO', RC )
          RETURN
       ENDIF

    ENDIF

    ! ----------------------------------------------------------------
    ! SESQ
    IF ( Inst%IDTSESQ > 0 ) THEN

       ! Add flux to emission array
       CALL HCO_EmisAdd( HcoState, Inst%FLUXSESQ, Inst%IDTSESQ, &
                         RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, &
                          'HCO_EmisAdd error: FLUXSESQ', RC )
          RETURN
       ENDIF

    ENDIF

    !=================================================================
    ! Manual diagnostics
    !=================================================================

    ! -------------------------------------------------------------
    ! Alpha Pinene
    DiagnName =  'InvMEGAN_APIN'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXAPIN, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Beta Pinene
    DiagnName =  'InvMEGAN_BPIN'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXBPIN, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Sabinene
    DiagnName =  'InvMEGAN_SABI'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXSABI, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Mycrene
    DiagnName =  'InvMEGAN_MYRC'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXMYRC, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! 3-Carene
    DiagnName =  'InvMEGAN_CARE'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXCARE, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Ocimene
    DiagnName =  'InvMEGAN_OCIM'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXOCIM, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Other monoterpenes
    DiagnName =  'InvMEGAN_OMON'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXOMON, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Total monoterpenes
    DiagnName =  'InvMEGAN_MONX'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXMONO, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! a-Farnesene
    DiagnName =  'InvMEGAN_FARN'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXFARN, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! b_Caryophyllene
    DiagnName =  'InvMEGAN_BCAR'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXBCAR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Other sesquiterpenes
    DiagnName =  'InvMEGAN_OSQT'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXOSQT, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Methyl butenol
    DiagnName =  'InvMEGAN_MBOX'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXMBOX, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! -------------------------------------------------------------
    ! Formic acid
    DiagnName =  'InvMEGAN_FAXX'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXFAXX, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Acetic acid
    DiagnName =  'InvMEGAN_AAXX'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(DiagnName), Array2D=Inst%FLUXAAXX, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Eventually copy internal values to ESMF internal state object
    ! ----------------------------------------------------------------

    ! LAI_PREVDAY
    CALL HCO_RestartWrite( HcoState, &
                           'LAI_PREVDAY', Inst%LAI_PREVDAY, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! T_LAST24H
    CALL HCO_RestartWrite( HcoState, &
                           'T_PREVDAY',  Inst%T_LAST24H, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! T_LASTXDAYS
    CALL HCO_RestartWrite( HcoState, &
                           'T_DAVG',     Inst%T_LASTXDAYS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! PARDR_LASTXDAYS
    CALL HCO_RestartWrite( HcoState, &
                           'PARDR_DAVG', Inst%PARDR_LASTXDAYS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! PARDF_LASTXDAYS
    CALL HCO_RestartWrite( HcoState, &
                           'PARDF_DAVG', Inst%PARDF_LASTXDAYS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! ALL DONE!
    !=================================================================

    !! testing only
    !write(*,*) ''
    !write(*,*) 'MEGAN done!'
    !write(*,*) 'total LAI : ', SUM(ExtState%GC_LAI%Arr%Val)
    !write(*,*) 'total ISOP: ', SUM(Inst%FLUXISOP)
    !write(*,*) 'total ACET: ', SUM(Inst%FLUXACETmo)+SUM(Inst%FLUXACETmb)+&
    !                           SUM(Inst%FLUXACETbg)
    !write(*,*) 'total PRPE: ', SUM(Inst%FLUXPRPE)
    !write(*,*) ''

    ! Cleanup
    Inst => NULL()

    ! Reset first-time flag
    FIRST = .FALSE.

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Megan_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Megan_Emissions
!
! !DESCRIPTION: Subroutine Get\_Megan\_Emissions computes biogenic emissions in
!  units of [kgC/m2/s] or [kg/m2/s] using the MEGAN inventory. (dbm, 12/2012)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_MEGAN_EMISSIONS( HcoState, ExtState, &
                                  Inst, I, J, CMPD, MEGAN_EMIS, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER     :: HcoState
    TYPE(Ext_State),  POINTER     :: ExtState
    TYPE(MyInst),     POINTER     :: Inst
    INTEGER,          INTENT(IN)  :: I, J      ! GEOS-Chem lon & lat indices
    CHARACTER(LEN=*), INTENT(IN)  :: CMPD      ! Compound name (dbm,6/21/2012)
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(OUT) :: MEGAN_EMIS ! VOC emission in kgC/m2/s or
                                                ! kg/m2/s, depending on units
                                                ! the compound is carried in
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL(hp)  :: GAMMA_CO2  ! (Tai, Jan 2013)
    REAL(hp)  :: AEF
    REAL(hp)  :: D_BTW_M
    REAL(hp)  :: TS, SUNCOS
    REAL(hp)  :: Q_DIR_2, Q_DIFF_2
    REAL(hp)  :: BETA, LDF, CT1, CEO
    REAL(hp)  :: ANEW, AGRO, AMAT, AOLD
    REAL(hp)  :: ISOLAI, PMISOLAI, MISOLAI
    REAL(hp)  :: PFTSUM
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
    GAMMA_CO2  = 0.0_hp  ! (Tai, Jan 2013)
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
    D_BTW_M  = DBLE( Inst%DAYS_BTW_M )

    ! Pass met variables. Now use only LAI (ckeller, 10/9/2014)
    ISOLAI   = ExtState%LAI%Arr%Val(I,J)
    PMISOLAI = Inst%LAI_PREVDAY(I,J)
    MISOLAI  = ISOLAI
    TS       = ExtState%T2M%Arr%Val(I,J)
    SUNCOS   = ExtState%SUNCOS%Arr%Val(I,J)

    ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
    Q_DIR_2  = ExtState%PARDR%Arr%Val(I,J) * WM2_TO_UMOLM2S
    Q_DIFF_2 = ExtState%PARDF%Arr%Val(I,J) * WM2_TO_UMOLM2S

    ! Eventually normalize LAI by PFT (if setting is set
    ! accordingly). ckeller, 7/17/17.
    IF ( Inst%NORMLAI ) THEN
       PFTSUM = SUM( Inst%ARRAY_16(I,J,2:16) )
       IF ( PFTSUM > 0.0_hp ) THEN
          MISOLAI  = MIN( MISOLAI  / PFTSUM, LAI_MAX )
          PMISOLAI = MIN( PMISOLAI / PFTSUM, LAI_MAX )
       ENDIF
    ENDIF

    ! --------------------------------------------
    ! Get MEGAN parameters for this compound
    ! --------------------------------------------
    CALL GET_MEGAN_PARAMS ( HcoState,                         &
                            CMPD, BETA, LDF,  CT1,  CEO,      &
                            ANEW, AGRO, AMAT, AOLD, BIDIR, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! --------------------------------------------
    ! Get base emission factor for this compound and grid square
    ! Units: kgC/m2/s or kg/m2/s
    ! --------------------------------------------
    CALL GET_MEGAN_AEF ( HcoState, Inst, I, J, CMPD, AEF, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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

          GAMMA_PAR = GET_GAMMA_PAR_PCEEA( HcoState,                  &
                                           ExtState, Inst, I, J,      &
                                           Q_DIR_2,  Q_DIFF_2,        &
                                           Inst%PARDR_LASTXDAYS(I,J), &
                                           Inst%PARDF_LASTXDAYS(I,J) )
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
       GAMMA_T_LD = GET_GAMMA_T_LD( TS, Inst%T_LASTXDAYS(I,J), &
                                    Inst%T_LAST24H(I,J), CT1, CEO )

       ! --------------------------------------------------
       ! GAMMA_LAI (leaf area index activity factor)
       ! --------------------------------------------------
       GAMMA_LAI = GET_GAMMA_LAI( MISOLAI, BIDIR )

       ! --------------------------------------------------
       ! GAMMA_AGE (leaf age activity factor)
       ! --------------------------------------------------
       GAMMA_AGE = GET_GAMMA_AGE( MISOLAI, PMISOLAI, D_BTW_M, &
                                  Inst%T_LASTXDAYS(I,J),      &
                                  ANEW, AGRO, AMAT, AOLD )

       ! --------------------------------------------------
       ! GAMMA_SM (soil moisture activity factor)
       ! --------------------------------------------------
       GAMMA_SM = GET_GAMMA_SM( ExtState, I, J, CMPD )

       ! CO2 inhibition of isoprene (Tai, Jan 2013)
       IF ( Inst%LISOPCO2 ) THEN
          GAMMA_CO2 = GET_GAMMA_CO2( Inst%GLOBCO2 )
       ELSE
          GAMMA_CO2 = 1.0_hp
       ENDIF

    ELSE

       ! set activity factors to zero
       GAMMA_PAR  = 0.0_hp
       GAMMA_T_LI = 0.0_hp
       GAMMA_T_LD = 0.0_hp
       GAMMA_LAI  = 0.0_hp
       GAMMA_AGE  = 0.0_hp
       GAMMA_SM   = 0.0_hp
       GAMMA_CO2  = 0.0_hp  ! (Tai, Jan 2013)

    END IF

    ! Emission is the product of all of these.
    ! Units here are kgC/m2/s or kg/m2/s as appropriate for the compound.
    ! Normalization factor ensures product of GAMMA values is 1.0 under
    !  standard conditions.
    IF ( CMPD == 'ISOP' ) THEN
       ! Only apply CO2 inhibition to isoprene (mps, 9/15/15)
       ! Amos Tai wrote:
       !  In my opinion, we should not apply the CO2 inhibition factor on
       !  other monoterpene species yet, because the empirical data I've used
       !  are for isoprene emissions only. But we generally agree that CO2
       !  inhibition should affect monoterpenes too, so we'll leave room for
       !  future incorporation when new data arise.
       MEGAN_EMIS = Inst%NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * &
                    GAMMA_LAI * ((1.0_hp - LDF) * GAMMA_T_LI +      &
                    (LDF * GAMMA_PAR * GAMMA_T_LD)) * GAMMA_CO2
    ELSE
       MEGAN_EMIS = Inst%NORM_FAC(1) * AEF * GAMMA_AGE * GAMMA_SM * &
                    GAMMA_LAI * ((1.0_hp - LDF) * GAMMA_T_LI +      &
                    (LDF * GAMMA_PAR * GAMMA_T_LD))
    ENDIF

    !! testing only
    !if ( i==ix .and. j==iy ) then
    !   write(*,*) ' '
    !   write(*,*) '--- GET_MEGAN_EMISSIONS --- '
    !   write(*,*) 'HEMCO MEGAN @    ', i,j
    !   write(*,*) 'Compound       : ', TRIM(CMPD)
    !   write(*,*) 'MEGAN_EMIS     : ', MEGAN_EMIS
    !   write(*,*) 'SUNCOS         : ', SUNCOS
    !   write(*,*) 'AEF [kgC/m2/s] : ', AEF
    !   write(*,*) 'GAMMA_LAI      : ', GAMMA_LAI
    !   write(*,*) 'GAMMA_AGE      : ', GAMMA_AGE
    !   write(*,*) 'GAMMA_T_LI     : ', GAMMA_T_LI
    !   write(*,*) 'GAMMA_T_LD     : ', GAMMA_T_LD
    !   write(*,*) 'GAMMA_PAR      : ', GAMMA_PAR
    !   write(*,*) 'GAMMA_SM       : ', GAMMA_SM
    !   write(*,*) 'TS             : ', TS
    !   write(*,*) 'HCOT_DAILY     : ', HCOT_DAILY(I,J)
    !   write(*,*) 'ISOLAI         : ', ISOLAI
    !   write(*,*) 'MISOLAI        : ', MISOLAI
    !   write(*,*) 'PMISOLAI       : ', PMISOLAI
    !   write(*,*) 'D_BTW_M        : ', D_BTW_M
    !   write(*,*) ' '
    !endif

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
  SUBROUTINE GET_MEGAN_PARAMS( HcoState,                           &
                               CPD,   BTA,   LIDF,  C_T1,  C_EO,   &
                               A_NEW, A_GRO, A_MAT, A_OLD, BI_DIR, &
                               RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER    :: HcoState
    CHARACTER(LEN=*), INTENT(IN) :: CPD       ! Compound name
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    IF ( TRIM(CPD) == 'ISOP' .OR. &
         TRIM(CPD) == 'MBOX' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'MYRC' .OR. &
              TRIM(CPD) == 'SABI' .OR. &
              TRIM(CPD) == 'APIN' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'LIMO' .OR. &
              TRIM(CPD) == 'CARE' .OR. &
              TRIM(CPD) == 'BPIN' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'MOH' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'EOH'  .OR. &
              TRIM(CPD) == 'CH2O' .OR. &
              TRIM(CPD) == 'ALD2' .OR. &
              TRIM(CPD) == 'FAXX' .OR. &
              TRIM(CPD) == 'AAXX' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'C2H4' .OR. &
              TRIM(CPD) == 'TOLU' .OR. &
              TRIM(CPD) == 'HCNX' ) THEN
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
    ELSE IF ( TRIM(CPD) == 'FARN' .OR. &
              TRIM(CPD) == 'BCAR' .OR. &
              TRIM(CPD) == 'OSQT' ) THEN
       BTA    = 0.17_hp
       LIDF   = 0.5_hp
       C_T1   = 130.0_hp
       C_EO   = 2.37_hp
       A_NEW  = 0.4_hp
       A_GRO  = 0.6_hp
       A_MAT  = 1.0_hp
       A_OLD  = 0.95_hp
       BI_DIR = .FALSE.

    ! Calls for any other MEGAN compounds (e.g. sesquiterpenes, etc.) can
    ! be added following the above format based on the parameters in
    ! Guenther 2012 or the MEGAN source code (dbm, 6/21/2012).
    ELSE

       MSG = 'Invalid compound name'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, &
                      THISLOC='GET_MEGAN_PARAMS' )
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
  SUBROUTINE GET_MEGAN_AEF( HcoState, Inst, I, J, CPD, EMFAC, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),   POINTER     :: HcoState
    TYPE(MyInst),      POINTER     :: Inst
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
!  See https://github.com/geoschem/geos-chem for complete history
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
       EMFAC = Inst%AEF_ISOP(I,J)
    CASE ( 'MBOX' )
       EMFAC = Inst%AEF_MBOX(I,J)
    CASE ( 'MYRC' )
       EMFAC = Inst%AEF_MYRC(I,J)
    CASE ( 'SABI' )
       EMFAC = Inst%AEF_SABI(I,J)
    CASE ( 'APIN' )
       EMFAC = Inst%AEF_APIN(I,J)
    CASE ( 'LIMO' )
       EMFAC = Inst%AEF_LIMO(I,J)
    CASE ( 'CARE' )
       EMFAC = Inst%AEF_CARE(I,J)
    CASE ( 'BPIN' )
       EMFAC = Inst%AEF_BPIN(I,J)
    CASE ( 'OCIM' )
       EMFAC = Inst%AEF_OCIM(I,J)
    CASE ( 'OMON' )
       EMFAC = Inst%AEF_OMON(I,J)
    CASE ( 'MOH' )
       EMFAC = Inst%AEF_MOH(I,J)
    CASE ( 'ACET' )
       EMFAC = Inst%AEF_ACET(I,J)
    CASE ( 'EOH' )
       EMFAC = Inst%AEF_EOH(I,J)
    CASE ( 'CH2O' )
       EMFAC = Inst%AEF_CH2O(I,J)
    CASE ( 'ALD2' )
       EMFAC = Inst%AEF_ALD2(I,J)
    CASE ( 'FAXX' )
       EMFAC = Inst%AEF_FAXX(I,J)
    CASE ( 'AAXX' )
       EMFAC = Inst%AEF_AAXX(I,J)
    CASE ( 'C2H4' )
       EMFAC = Inst%AEF_C2H4(I,J)
    CASE ( 'TOLU' )
       EMFAC = Inst%AEF_TOLU(I,J)
    CASE ( 'HCNX' )
       EMFAC = Inst%AEF_HCNX(I,J)
    CASE ( 'PRPE' )
       EMFAC = Inst%AEF_PRPE(I,J)
    CASE ( 'FARN' )
       EMFAC = Inst%AEF_FARN(I,J)
    CASE ( 'BCAR' )
       EMFAC = Inst%AEF_BCAR(I,J)
    CASE ( 'OSQT' )
       EMFAC = Inst%AEF_OSQT(I,J)
    CASE DEFAULT
       MSG = 'Invalid compound name'
       CALL HCO_ERROR(HcoState%Config%Err, MSG, &
                      RC, THISLOC='GET_MEGAN_AEF' )
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
  FUNCTION GET_GAMMA_PAR_PCEEA( HcoState, ExtState,            &
                                Inst, I, J, Q_DIR_2, Q_DIFF_2, &
                                PARDR_AVG_SIM, PARDF_AVG_SIM ) &
    RESULT( GAMMA_P_PCEEA )
!
! !USES:
!
    USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER    :: HcoState
    TYPE(Ext_State), POINTER    :: ExtState
    TYPE(MyInst),    POINTER    :: Inst
    INTEGER,         INTENT(IN) :: I,  J             ! Lon & lat indices
    REAL(sp),        INTENT(IN) :: PARDR_AVG_SIM     ! Avg direct PAR [W/m2]
    REAL(sp),        INTENT(IN) :: PARDF_AVG_SIM     ! Avg diffuse PAR [W/m2]
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
!  (1 ) Code was taken & adapted directly from the MEGAN v2.1 source code.
!      (mpb,2009)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp)   :: mmPARDR_DAILY
    REAL(hp)   :: mmPARDF_DAILY
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
    C_PPFD = 0.0_hp
    PTOA   = 0.0_hp

    ! Convert past light conditions to micromol/m2/s
    mmPARDR_DAILY = PARDR_AVG_SIM  * WM2_TO_UMOLM2S
    mmPARDF_DAILY = PARDF_AVG_SIM  * WM2_TO_UMOLM2S

    ! Work out the light at the top of the canopy.
    PAC_DAILY    = mmPARDR_DAILY + mmPARDF_DAILY
    PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

    ! Get latitude
    LAT = HcoState%Grid%YMID%Val(I,J)

    ! Get day of year, local-time and latitude
    ! TODO: Evaluate RC?
    CALL HcoClock_Get( HcoState%Clock, cDOY = DOY, RC=RC )
    CALL HcoClock_GetLocal( HcoState, I, J, cH = LocalHour, RC=RC )

    ! Get solar elevation angle
    SINbeta =  SOLAR_ANGLE( HcoState, Inst, DOY, LocalHour, LAT )
    BETA    =  ASIN( SINbeta ) * Inst%RAD2D

    IF ( SINbeta < 0.0_hp ) THEN

       GAMMA_P_PCEEA = 0.0_hp

    ELSEIF ( SINbeta > 0.0_hp ) THEN

       ! PPFD at top of atmosphere
       PTOA    = 3000.0_hp + 99.0_hp * &
                 COS( 2._hp * 3.14159265358979323_hp * &
                 ( DOY - 10.0_hp ) / 365.0_hp )

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

    !! testing only
    !if ( i==ix .and. j==iy ) then
    !   write(*,*) ' '
    !   write(*,*) 'HEMCO GAMMA_PAR_PCEEA: ', GAMMA_P_PCEEA
    !   write(*,*) 'DOY                  : ', DOY
    !   write(*,*) 'LAT                  : ', LAT
    !   write(*,*) 'SINbeta              : ', SINbeta
    !   write(*,*) 'BETA                 : ', BETA
    !   write(*,*) 'PTOA                 : ', PTOA
    !   write(*,*) 'PHI                  : ', PHI
    !   write(*,*) 'BBB                  : ', BBB
    !   write(*,*) 'AAA                  : ', AAA
    !   write(*,*) 'PAC_DAILY            : ', PAC_DAILY
    !   write(*,*) 'PAC_INSTANT          : ', PAC_INSTANT
    !endif

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
  FUNCTION SOLAR_ANGLE( HcoState, Inst, DOY, SHOUR, LAT ) &
       RESULT(SINbeta)
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),   POINTER    :: HcoState
    TYPE(MyInst),      POINTER    :: Inst
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(hp) :: BETA                        ! solar elevation angle
    REAL(hp) :: sindelta, cosdelta, A, B

    ! Calculation of sin beta
    sindelta = -SIN( 0.40907_hp ) * COS( 6.28_hp * ( DOY + 10_dp ) / 365_dp )

    cosdelta = (1-sindelta**2.0_hp)**0.5_hp

    A = SIN( LAT * Inst%D2RAD ) * sindelta
    B = COS( LAT * Inst%D2RAD ) * cosdelta

    SINbeta = A + B * COS( 2.0_hp * HcoState%Phys%PI * ( SHOUR-12_dp )/24_dp )

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
  FUNCTION GET_GAMMA_T_LI( T, BETA ) &
       RESULT( GAMMA_T_LI )
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
!  Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
!  (3 ) Guenther et al., GMD 2012 and MEGANv2.1 source code.
!
! !REVISION HISTORY:
!  (1 ) Original code by Michael Barkley (2009).
!  See https://github.com/geoschem/geos-chem for complete history
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
  FUNCTION GET_GAMMA_T_LD( T, PT_15, PT_1, CT1, CEO ) &
       RESULT( GAMMA_T_LD )
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    REAL(hp), PARAMETER   :: R   = 8.3144598e-3_hp

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
    C_T   = E_OPT * CT2 * EXP( CT1 * X ) / &
            ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

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
  FUNCTION GET_GAMMA_LAI( CMLAI, BIDIREXCH ) &
       RESULT( GAMMA_LAI )
!
! !INPUT PARAMETERS:
!
    REAL(hp),     INTENT(IN) :: CMLAI       ! Current month's LAI [cm2/cm2]
    LOGICAL,      INTENT(IN) :: BIDIREXCH   ! Logical flag indicating whether
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

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
       GAMMA_LAI = 0.49_hp * CMLAI / SQRT( 1.0_hp + 0.2_hp * CMLAI*CMLAI)
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
  FUNCTION GET_GAMMA_AGE( CMLAI, PMLAI, DBTWN, TT, AN, AG, AM, AO ) &
       RESULT( GAMMA_AGE )
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
!  See https://github.com/geoschem/geos-chem for complete history
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
          FMAT = ( PMLAI / CMLAI ) + &
                 (( DBTWN - TM ) / DBTWN )*( 1.0_hp -  PMLAI / CMLAI )
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
    GAMMA_AGE = FNEW * AN + FGRO * AG + FMAT * AM + FOLD * AO

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
  FUNCTION GET_GAMMA_SM( ExtState, I, J, CMPD ) &
       RESULT( GAMMA_SM )
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
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

    ! Error trap: GWETROOT must be between 0.0 and 1.0 (ckeller, 4/16/15)
    GWETROOT = MIN(MAX(ExtState%GWETROOT%Arr%Val(I,J),0.0_hp),1.0_hp)

    IF ( TRIM( CMPD ) == 'ALD2' .OR. TRIM ( CMPD ) == 'EOH' ) THEN

       ! GWETROOT = degree of saturation or wetness in the root-zone
       ! (top meter of soil). This is defined as the ratio of the volumetric
       ! soil moisture to the porosity. We use a soil moisture activity factor
       ! for ALD2 to account for stimulation of emission by flooding.
       ! (Millet et al., ACP 2010)
       ! Constant value of 1.0 for GWETROOT = 0-0.9, increasing linearly to
       ! 3.0 at GWETROOT =1.
       GAMMA_SM = MAX( 20.0_hp * GWETROOT - 17.0_hp, 1.0_hp)

    ENDIF

  END FUNCTION GET_GAMMA_SM
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gamma_co2
!
! !DESCRIPTION: Function GET\_GAMMA\_CO2 computes the CO2 activity factor
!  associated with CO2 inhibition of isoprene emission. Called from
!  GET\_MEGAN\_EMISSIONS only.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_GAMMA_CO2( CO2a ) &
       RESULT( GAMMA_CO2 )
!
! !INPUT PARAMETERS:
    REAL(hp), INTENT(IN) :: CO2a       ! Atmospheric CO2 conc [ppmv]
!
! !RETURN VALUE:
    REAL(hp)             :: GAMMA_CO2  ! CO2 activity factor [unitless]
!
! !LOCAL VARIABLES:
    REAL(hp)             :: CO2i       ! Intercellular CO2 conc [ppmv]
    REAL(hp)             :: ISMAXi     ! Asymptote for intercellular CO2
    REAL(hp)             :: HEXPi      ! Exponent for intercellular CO2
    REAL(hp)             :: CSTARi     ! Scaling coef for intercellular CO2
    REAL(hp)             :: ISMAXa     ! Asymptote for atmospheric CO2
    REAL(hp)             :: HEXPa      ! Exponent for atmospheric CO2
    REAL(hp)             :: CSTARa     ! Scaling coef for atmospheric CO2
    LOGICAL              :: LPOSSELL   ! Use Possell & Hewitt (2011)?
    LOGICAL              :: LWILKINSON ! Use Wilkinson et al. (2009)?
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Heald, C. L., Wilkinson, M. J., Monson, R. K., Alo, C. A.,
!       Wang, G. L., and Guenther, A.: Response of isoprene emission
!       to ambient co(2) changes and implications for global budgets,
!       Global Change Biology, 15, 1127-1140, 2009.
!  (2 ) Wilkinson, M. J., Monson, R. K., Trahan, N., Lee, S., Brown, E.,
!       Jackson, R. B., Polley, H. W., Fay, P. A., and Fall, R.: Leaf
!       isoprene emission rate as a function of atmospheric CO2
!       concentration, Global Change Biology, 15, 1189-1200, 2009.
!  (3 ) Possell, M., and Hewitt, C. N.: Isoprene emissions from plants
!       are mediated by atmospheric co2 concentrations, Global Change
!       Biology, 17, 1595-1610, 2011.
!
! !REVISION HISTORY:
!  (1 ) Implemented in the standard code by A. Tai (Jun 2012).
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !----------------------------------------------------------
    ! Choose between two alternative CO2 inhibition schemes
    !----------------------------------------------------------

    ! Empirical relationship of Possell & Hewitt (2011) based on nine
    ! experimental studies including Wilkinson et al. (2009). This is
    ! especially recommended for sub-ambient CO2 concentrations:
    LPOSSELL    = .TRUE.   ! Default option

    ! Semi-process-based parameterization of Wilkinson et al. (2009),
    ! taking into account of sensitivity to intercellular CO2
    ! fluctuation, which is here set as a constant fraction of
    ! atmospheric CO2:
    LWILKINSON  = .FALSE.   ! Set .TRUE. only if LPOSSELL = .FALSE.

    !-----------------------
    ! Compute GAMMA_CO2
    !-----------------------

    IF ( LPOSSELL ) THEN

       ! Use empirical relationship of Possell & Hewitt (2011):
       GAMMA_CO2 = 8.9406_hp / ( 1.0_hp + 8.9406_hp * 0.0024_hp * CO2a )

    ELSEIF ( LWILKINSON ) THEN

       ! Use parameterization of Wilkinson et al. (2009):

       ! Parameters for intercellular CO2 using linear interpolation:
       IF ( CO2a <= 600.0_hp ) THEN
          ISMAXi = 1.036_hp  - (1.036_hp - 1.072_hp) / &
                   (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
          HEXPi  = 2.0125_hp - (2.0125_hp - 1.7000_hp) / &
                   (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
          CSTARi = 1150.0_hp - (1150.0_hp - 1218.0_hp) / &
                   (600.0_hp - 400.0_hp) * (600.0_hp - CO2a)
       ELSEIF ( CO2a > 600.0_hp .AND. CO2a < 800.0_hp ) THEN
          ISMAXi = 1.046_hp  - (1.046_hp - 1.036_hp) / &
                   (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
          HEXPi  = 1.5380_hp - (1.5380_hp - 2.0125_hp) / &
                   (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
          CSTARi = 2025.0_hp - (2025.0_hp - 1150.0_hp) / &
                   (800.0_hp - 600.0_hp) * (800.0_hp - CO2a)
       ELSE
          ISMAXi = 1.014_hp - (1.014_hp - 1.046_hp) / &
                   (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
          HEXPi  = 2.8610_hp - (2.8610_hp - 1.5380_hp) / &
                   (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
          CSTARi = 1525.0_hp - (1525.0_hp - 2025.0_hp) / &
                   (1200.0_hp - 800.0_hp) * (1200.0_hp - CO2a)
       ENDIF

       ! Parameters for atmospheric CO2:
       ISMAXa    = 1.344_hp
       HEXPa     = 1.4614_hp
       CSTARa    = 585.0_hp

       ! For now, set CO2_Ci = 0.7d0 * CO2_Ca as recommended by Heald
       ! et al. (2009):
       CO2i      = 0.7_hp * CO2a

       ! Compute GAMMA_CO2:
       GAMMA_CO2 = ( ISMAXi -  ISMAXi * CO2i**HEXPi / &
                   ( CSTARi**HEXPi + CO2i**HEXPi ) )  &
                   * ( ISMAXa - ISMAXa * ( 0.7_hp * CO2a )**HEXPa / &
                   ( CSTARa**HEXPa + ( 0.7_hp * CO2a )**HEXPa ) )

    ELSE

       ! No CO2 inhibition scheme is used; GAMMA_CO2 set to unity:
       GAMMA_CO2 = 1.0_hp

    ENDIF

  END FUNCTION GET_GAMMA_CO2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_NORM_FAC
!
! !DESCRIPTION: Function CALC\_NORM\_FAC calculates the normalization factor
!  needed to compute emissions. Called from GET\_MEGAN\_EMISSIONS.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_NORM_FAC( Inst, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(MyInst),    POINTER        :: Inst
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    GAMMA_LAI_STANDARD = 0.49_hp * CMLAI / SQRT( 1.0_hp + 0.2_hp * CMLAI*CMLAI )

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
    R     = 8.3144598e-3_hp
    ! parameters for isoprene
    CEO = 2.0_hp
    CT1 = 95.0_hp

    E_OPT = CEO * EXP( 0.08_hp * ( PT_15  - 2.97e2_hp ) )
    T_OPT = 3.13e2_hp + ( 6.0e-1_hp * ( PT_15 - 2.97e2_hp ) )
    CT2   = 200.0_hp

    ! Variable related to temperature
    X     = ( 1.0_hp/T_OPT - 1.0_hp/T ) / R

    GAMMA_T_LD_STANDARD = E_OPT * CT2 * EXP( CT1 * X ) / &
                          ( CT2 - CT1 * ( 1.0_hp - EXP( CT2 * X ) ) )

    ! -----------------
    ! Overall GAMMA_STANDARD
    ! -----------------
    ! LDF = 1d0 for isoprene
    LDF = 1.0_hp
    GAMMA_STANDARD = &
         GAMMA_AGE_STANDARD * GAMMA_SM_STANDARD * GAMMA_LAI_STANDARD &
         * ((1.0_hp - LDF) * GAMMA_T_LI_STANDARD &
         + (LDF * GAMMA_P_STANDARD * GAMMA_T_LD_STANDARD))
    ! This ends up being 1.0101081.

    Inst%NORM_FAC = 1.0_hp / GAMMA_STANDARD

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
  SUBROUTINE FILL_RESTART_VARS( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_Restart_Mod,    ONLY : HCO_RestartGet
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState
    TYPE(Ext_State), POINTER        :: ExtState
    TYPE(MyInst),    POINTER        :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  19 Dec 2014 - C. Keller    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
    CALL HCO_RestartGet( HcoState, 'LAI_PREVDAY', &
                         Inst%LAI_PREVDAY, RC, FILLED=FND )

    ! Default value
    IF ( .NOT. FND ) THEN
       Inst%LAI_PREVDAY = ExtState%LAI%Arr%Val
    ENDIF

    !-----------------------------------------------------------------
    ! DIFFUSE RADIATION ( PARDF )
    !-----------------------------------------------------------------

    ! Temperature over last 24 hours
    CALL HCO_RestartGet( HcoState, 'T_PREVDAY', &
                         Inst%T_LAST24H, RC, DefVal = 288.15 )

    ! Temperature over last X days
    CALL HCO_RestartGet( HcoState, 'T_DAVG', &
                         Inst%T_LASTXDAYS, RC, DefVal = 288.15 )

    ! Direct radiation (PARDR) over last X days
    CALL HCO_RestartGet( HcoState, 'PARDR_DAVG', &
                         Inst%PARDR_LASTXDAYS, RC, DefVal = 30.0 )

    ! Diffuse radiation (PARDF) over last X days
    CALL HCO_RestartGet( HcoState, 'PARDF_DAVG', &
                         Inst%PARDF_LASTXDAYS, RC, DefVal = 48.0 )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE FILL_RESTART_VARS
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
  SUBROUTINE CALC_AEF( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_EMISLIST_MOD,    ONLY : HCO_GetPtr
    USE HCO_CALC_MOD,        ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState
    TYPE(HCO_State), POINTER       :: HcoState
    TYPE(MyInst),    POINTER       :: Inst
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: I, J, P, ARR_IND
    REAL(hp)                :: FACTOR, SPECIES2CARBON
    REAL(hp)                :: PFT_EF_OMON(15), PFT_EF_MOH(15)
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
    REAL(hp)                :: EM_FRAC_ALD2(15), EM_FRAC_EOH(15)
    REAL(hp)                :: EM_FRAC_FAXX(15), EM_FRAC_AAXX(15)
    REAL(hp)                :: EM_FRAC_CH2O(15)

    ! Pointers to CLM4 plant functional types
    REAL(hp)  :: PFT_BARE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_NDLF_EVGN_TMPT_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_NDLF_EVGN_BORL_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_NDLF_DECD_BORL_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_EVGN_TROP_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_EVGN_TMPT_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_DECD_TROP_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_DECD_TMPT_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_DECD_BORL_TREE(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_EVGN_SHRB(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_DECD_TMPT_SHRB(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_BDLF_DECD_BORL_SHRB(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_C3_ARCT_GRSS(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_C3_NARC_GRSS(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_C4_GRSS(HcoState%NX,HcoState%NY)
    REAL(hp)  :: PFT_CROP(HcoState%NX,HcoState%NY)

    ! Suffix
    CHARACTER(LEN=255)      :: SFX

    !=================================================================
    ! CALC_AEF begins here!
    !=================================================================

    ! Suffix
    SFX = Inst%SUFFIX

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
    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_ISOP'//TRIM(SFX), &
                      Inst%AEF_ISOP, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_MBOX'//TRIM(SFX), &
                      Inst%AEF_MBOX, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_BPIN'//TRIM(SFX), &
                      Inst%AEF_BPIN, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_CARE'//TRIM(SFX), &
                      Inst%AEF_CARE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_LIMO'//TRIM(SFX), &
                      Inst%AEF_LIMO, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_OCIM'//TRIM(SFX), &
                      Inst%AEF_OCIM, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_AEF_SABI'//TRIM(SFX), &
                      Inst%AEF_SABI, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, 'MEGAN_ORVC'//TRIM(SFX), &
                      Inst%GEIA_ORVC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BARE'//TRIM(SFX), &
                      PFT_BARE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_NDLF_EVGN_TMPT_TREE'//TRIM(SFX), &
                      PFT_NDLF_EVGN_TMPT_TREE, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_NDLF_EVGN_BORL_TREE'//TRIM(SFX), &
                      PFT_NDLF_EVGN_BORL_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_NDLF_DECD_BORL_TREE'//TRIM(SFX), &
                      PFT_NDLF_DECD_BORL_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_EVGN_TROP_TREE'//TRIM(SFX), &
                      PFT_BDLF_EVGN_TROP_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_EVGN_TMPT_TREE'//TRIM(SFX), &
                      PFT_BDLF_EVGN_TMPT_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_DECD_TROP_TREE'//TRIM(SFX), &
                      PFT_BDLF_DECD_TROP_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_DECD_TMPT_TREE'//TRIM(SFX), &
                      PFT_BDLF_DECD_TMPT_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_DECD_BORL_TREE'//TRIM(SFX), &
                      PFT_BDLF_DECD_BORL_TREE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_EVGN_SHRB'//TRIM(SFX), &
                      PFT_BDLF_EVGN_SHRB, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_DECD_TMPT_SHRB'//TRIM(SFX), &
                      PFT_BDLF_DECD_TMPT_SHRB, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_BDLF_DECD_BORL_SHRB'//TRIM(SFX), &
                      PFT_BDLF_DECD_BORL_SHRB, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_C3_ARCT_GRSS'//TRIM(SFX), &
                      PFT_C3_ARCT_GRSS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_C3_NARC_GRSS'//TRIM(SFX), &
                      PFT_C3_NARC_GRSS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_C4_GRSS'//TRIM(SFX), &
                      PFT_C4_GRSS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_EvalFld( HcoState, &
                      'CLM4_PFT_CROP'//TRIM(SFX), &
                      PFT_CROP, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Copy PFTs into ARRAY_16
    Inst%ARRAY_16(:,:, 1) = PFT_BARE
    Inst%ARRAY_16(:,:, 2) = PFT_NDLF_EVGN_TMPT_TREE
    Inst%ARRAY_16(:,:, 3) = PFT_NDLF_EVGN_BORL_TREE
    Inst%ARRAY_16(:,:, 4) = PFT_NDLF_DECD_BORL_TREE
    Inst%ARRAY_16(:,:, 5) = PFT_BDLF_EVGN_TROP_TREE
    Inst%ARRAY_16(:,:, 6) = PFT_BDLF_EVGN_TMPT_TREE
    Inst%ARRAY_16(:,:, 7) = PFT_BDLF_DECD_TROP_TREE
    Inst%ARRAY_16(:,:, 8) = PFT_BDLF_DECD_TMPT_TREE
    Inst%ARRAY_16(:,:, 9) = PFT_BDLF_DECD_BORL_TREE
    Inst%ARRAY_16(:,:,10) = PFT_BDLF_EVGN_SHRB
    Inst%ARRAY_16(:,:,11) = PFT_BDLF_DECD_TMPT_SHRB
    Inst%ARRAY_16(:,:,12) = PFT_BDLF_DECD_BORL_SHRB
    Inst%ARRAY_16(:,:,13) = PFT_C3_ARCT_GRSS
    Inst%ARRAY_16(:,:,14) = PFT_C3_NARC_GRSS
    Inst%ARRAY_16(:,:,15) = PFT_C4_GRSS
    Inst%ARRAY_16(:,:,16) = PFT_CROP

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

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_OMON = (/180.0_hp, 180.0_hp, 170.0_hp, 150.0_hp, 150.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    150.0_hp, 150.0_hp, 150.0_hp, 110.0_hp, 200.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    110.0_hp, 5.0_hp  , 5.0_hp  , 5.0_hp  , 5.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_MOH  = (/900.0_hp, 900.0_hp, 900.0_hp, 500.0_hp, 900.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    500.0_hp, 900.0_hp, 900.0_hp, 900.0_hp, 900.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    900.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 900.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_ACET = (/240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, 240.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    240.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_BIDR = (/500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, 500.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    500.0_hp, 80.0_hp , 80.0_hp , 80.0_hp , 80.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_STRS = (/300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp, 300.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_OTHR = (/140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp, 140.0_hp/)

    ! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_APIN = (/500.0_hp, 500.0_hp, 510.0_hp, 600.0_hp, 400.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    600.0_hp, 400.0_hp, 400.0_hp, 200.0_hp, 300.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    200.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_MYRC = (/70.0_hp,  70.0_hp,  60.0_hp,  80.0_hp,  30.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    80.0_hp,  30.0_hp,  30.0_hp,  30.0_hp,  50.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    30.0_hp,  0.3_hp,   0.3_hp,   0.3_hp,   0.3_hp/)
    ! <---

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_FARN = (/40.0_hp,  40.0_hp,  40.0_hp,  60.0_hp,  40.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    60.0_hp,  40.0_hp,  40.0_hp,  40.0_hp,  40.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    40.0_hp,  3.0_hp,   3.0_hp,   3.0_hp,   4.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_BCAR = (/80.0_hp,  80.0_hp,  80.0_hp,  60.0_hp,  40.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    60.0_hp,  40.0_hp,  40.0_hp,  50.0_hp,  50.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    50.0_hp,  1.0_hp,   1.0_hp,   1.0_hp,   4.0_hp/)

    !               EF1       EF2       EF3       EF4       EF5
    PFT_EF_OSQT = (/120.0_hp, 120.0_hp, 120.0_hp, 120.0_hp, 100.0_hp, &
    !               EF6       EF7       EF8       EF9       EF10
                    120.0_hp, 100.0_hp, 100.0_hp, 100.0_hp, 100.0_hp, &
    !               EF11      EF12      EF13      EF14      EF15
                    100.0_hp, 2.0_hp,   2.0_hp,   2.0_hp,   2.0_hp/)


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
    EM_FRAC_ALD2 = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, &
                     0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, &
                     0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

    ! Ethanol: 40% of bidirectional category flux, except 25%
    ! for grasses and crops
    EM_FRAC_EOH  = (/0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, &
                     0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, 0.40_hp, &
                     0.40_hp, 0.25_hp, 0.25_hp, 0.25_hp, 0.25_hp/)

    ! Formic acid: 6% of bidirectional category flux, except 15%
    ! for grasses and crops
    EM_FRAC_FAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, &
                     0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, &
                     0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

    ! Acetic acid: 6% of bidirectional category flux, except 15%
    ! for grasses and crops
    EM_FRAC_AAXX = (/0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, &
                     0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, 0.06_hp, &
                     0.06_hp, 0.15_hp, 0.15_hp, 0.15_hp, 0.15_hp/)

    ! Formaldehyde: 8% of bidirectional category flux, except 20%
    ! for grasses and crops
    EM_FRAC_CH2O = (/0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, &
                     0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, 0.08_hp, &
                     0.08_hp, 0.20_hp, 0.20_hp, 0.20_hp, 0.20_hp/)

    ! Loop through plant types
    DO P = 1, 15

       ! Add 1 to Array_16 index to skip bare ground
       ARR_IND = P + 1

       ! Don't need to divide ARRAY_16 by 100 since data from netCDF
       ! file has already been converted to fraction (mps, 2/12/15)

       ! ---> Now compute EFs for a-pinene and myrcene as well (dbm, 12/2012)
       ! a-pinene: 100% of category
       Inst%AEF_APIN(:,:) = Inst%AEF_APIN(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_APIN(P)

       ! Myrcene: 100% of category
       Inst%AEF_MYRC(:,:) = Inst%AEF_MYRC(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_MYRC(P)
       ! <---

       ! Other monoterpenes: 100% of category
       Inst%AEF_OMON(:,:) = Inst%AEF_OMON(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OMON(P)

       ! a-Farnesene: 100% of category
       Inst%AEF_FARN(:,:) = Inst%AEF_FARN(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_FARN(P)

       ! b-Caryophyllene: 100% of category
       Inst%AEF_BCAR(:,:) = Inst%AEF_BCAR(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_BCAR(P)

       ! Other sesquiterpenes: 100% of category
       Inst%AEF_OSQT(:,:) = Inst%AEF_OSQT(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OSQT(P)

       ! Methanol: 100% of category
       Inst%AEF_MOH(:,:) = Inst%AEF_MOH(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_MOH(P)

       ! Acetone: 100% of category
       Inst%AEF_ACET(:,:) = Inst%AEF_ACET(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_ACET(P)

       ! Ethanol: variable fraction of category
       Inst%AEF_EOH(:,:) = Inst%AEF_EOH(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_EOH(P)

       ! Formaldehyde: variable fraction of category
       Inst%AEF_CH2O(:,:) = Inst%AEF_CH2O(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_CH2O(P)

       ! Acetaldehyde: variable fraction of category
       Inst%AEF_ALD2(:,:) = Inst%AEF_ALD2(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_ALD2(P)

       ! Formic acid: variable fraction of category
       Inst%AEF_FAXX(:,:) = Inst%AEF_FAXX(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_FAXX(P)

       ! Acetic acid: variable fraction of category
       Inst%AEF_AAXX(:,:) = Inst%AEF_AAXX(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND)*PFT_EF_BIDR(P)*EM_FRAC_AAXX(P)

       ! Ethene: 58% of "stress" category for all PFTs
       Inst%AEF_C2H4(:,:) = Inst%AEF_C2H4(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.58_hp

       ! Toluene: 3% of "stress" category for all PFTs
       Inst%AEF_TOLU(:,:) = Inst%AEF_TOLU(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.03_hp

       ! HCN: 1.5% of "stress" category for all PFTs
       Inst%AEF_HCNX(:,:) = Inst%AEF_HCNX(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_STRS(P) * 0.015_hp

       ! Propene: 48% of "other" category for all PFTs
       ! Butene:  24% of "other" category for all PFTs
       ! Larger alkenes: 0.2% of "other" category for all PFTs
       ! Total: 72.2%
       Inst%AEF_PRPE(:,:) = Inst%AEF_PRPE(:,:) + &
            Inst%ARRAY_16(:,:,ARR_IND) * PFT_EF_OTHR(P) * 0.722_hp

    ENDDO

    !! Nullify pointers
    !PFT_BARE                => NULL()
    !PFT_NDLF_EVGN_TMPT_TREE => NULL()
    !PFT_NDLF_EVGN_BORL_TREE => NULL()
    !PFT_NDLF_DECD_BORL_TREE => NULL()
    !PFT_BDLF_EVGN_TROP_TREE => NULL()
    !PFT_BDLF_EVGN_TMPT_TREE => NULL()
    !PFT_BDLF_DECD_TROP_TREE => NULL()
    !PFT_BDLF_DECD_TMPT_TREE => NULL()
    !PFT_BDLF_DECD_BORL_TREE => NULL()
    !PFT_BDLF_EVGN_SHRB      => NULL()
    !PFT_BDLF_DECD_TMPT_SHRB => NULL()
    !PFT_BDLF_DECD_BORL_SHRB => NULL()
    !PFT_C3_ARCT_GRSS        => NULL()
    !PFT_C3_NARC_GRSS        => NULL()
    !PFT_C4_GRSS             => NULL()
    !PFT_CROP                => NULL()

    ! Conversion factor from [ug compound/m2/hr] to [kg compound/m2/s]
    FACTOR = 1.0e-9_hp / 3600.0_hp

    ! Loop over grid boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Convert AEF arrays to [kgC/m2/s]
       ! Multiply arrays by FACTOR and ratio [g C/g compound]
       ! NOTE: AEFs for ISOP, MBOX, BPIN, CARE, LIMO, OCIM, SABI
       ! are read from file in [kgC/m2/s], so no need to convert here
       Inst%AEF_APIN(I,J) = Inst%AEF_APIN(I,J) * FACTOR * 120.0_hp / 136.234_hp
       Inst%AEF_MYRC(I,J) = Inst%AEF_MYRC(I,J) * FACTOR * 120.0_hp / 136.234_hp
       Inst%AEF_OMON(I,J) = Inst%AEF_OMON(I,J) * FACTOR * 120.0_hp / 136.234_hp

       ! Sesquiterpenes
       SPECIES2CARBON = 15.d0 * 12.01d0 / ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
       Inst%AEF_FARN(I,J) = Inst%AEF_FARN(I,J) * FACTOR * SPECIES2CARBON
       Inst%AEF_BCAR(I,J) = Inst%AEF_BCAR(I,J) * FACTOR * SPECIES2CARBON
       Inst%AEF_OSQT(I,J) = Inst%AEF_OSQT(I,J) * FACTOR * SPECIES2CARBON

       Inst%AEF_ACET(I,J) = Inst%AEF_ACET(I,J) * FACTOR *  36.0_hp /  58.079_hp
       Inst%AEF_EOH(I,J)  = Inst%AEF_EOH(I,J)  * FACTOR *  24.0_hp /  46.068_hp
       Inst%AEF_ALD2(I,J) = Inst%AEF_ALD2(I,J) * FACTOR *  24.0_hp /  44.053_hp
       Inst%AEF_C2H4(I,J) = Inst%AEF_C2H4(I,J) * FACTOR *  24.0_hp /  28.053_hp
       Inst%AEF_TOLU(I,J) = Inst%AEF_TOLU(I,J) * FACTOR *  84.0_hp /  92.138_hp
       Inst%AEF_PRPE(I,J) = Inst%AEF_PRPE(I,J) * FACTOR *  36.0_hp /  42.080_hp

       ! Methanol, formaldehyde, formic acid, acetic acid, HCN are
       ! carried in kg, not kg C
       ! Convert AEF arrays to [kg/m2/s]
       Inst%AEF_MOH(I,J)  = Inst%AEF_MOH(I,J)  * FACTOR
       Inst%AEF_CH2O(I,J) = Inst%AEF_CH2O(I,J) * FACTOR
       Inst%AEF_FAXX(I,J) = Inst%AEF_FAXX(I,J) * FACTOR
       Inst%AEF_AAXX(I,J) = Inst%AEF_AAXX(I,J) * FACTOR
       Inst%AEF_HCNX(I,J) = Inst%AEF_HCNX(I,J) * FACTOR

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
  SUBROUTINE HCOX_Megan_Init( HcoState, ExtName, ExtState, RC )
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
    TYPE(HCO_State),  POINTER        :: HcoState
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName
    TYPE(Ext_State),  POINTER        :: ExtState
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  05 Aug 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: ExtNr, iExtNr
    INTEGER                        :: I, nSpc, AS, NX, NY
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    LOGICAL                        :: FOUND
    REAL*8                         :: PI_180
    REAL(hp), POINTER              :: Ptr2D(:,:)
    TYPE(MyInst), POINTER          :: Inst
    CHARACTER(LEN=255)             :: MSG
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    LOGICAL                        :: Optfound

    !=================================================================
    ! HCOX_MEGAN_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, &
                    'HCOX_Megan_Init (hcox_megan_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Nullify
    Ptr2D => NULL()
    Inst  => NULL()

    ! Create an instance for this extension
    CALL InstCreate ( ExtNr, ExtState%Megan, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, &
                       'Cannot create MEGAN instance', RC )
       RETURN
    ENDIF

    ! Check to see if offline biogenic VOCs are used. If so, those
    ! will be used in place of online MEGAN emissions (mps, 12/17/19).
    Inst%OFFLINE_BIOGENICVOC = .FALSE.
    CALL GetExtOpt( HcoState%Config, 0, 'OFFLINE_BIOGENICVOC', &
                    OptValBool=OptFound, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) Inst%OFFLINE_BIOGENICVOC = OptFound

    ! Verbose mode
    IF ( HcoState%amIRoot) THEN
       MSG = 'Use MEGAN biogenic emissions (extension module)'
       CALL HCO_MSG( HcoState%Config%Err, MSG, SEP1='-' )
       WRITE(MSG,*) '- Use offline biogenic VOCs? ', Inst%OFFLINE_BIOGENICVOC
       CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDIF

    !-----------------------------------------------------------------
    ! Read settings
    !-----------------------------------------------------------------

    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in
    !       the config. file!
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Isoprene scaling', &
                    OptValHp=Inst%ISOP_SCALING, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CO2 inhibition', &
                    OptValBool=Inst%LISOPCO2, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CO2 conc (ppmv)', &
                    OptValHp=Inst%GLOBCO2, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Normalize LAI by PFT? Default setting is 'yes'
    ! ckeller, 7/17/17.
    Inst%NORMLAI = .TRUE.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Normalize LAI', &
                    OptValBool=OptFound, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) Inst%NORMLAI = OptFound

    ! Check GLOBCO2 if CO2 inhibition is turned on (LISOPCO2 = .TRUE.)
    ! GLOBCO2 should be between 150-1250 ppmv. Isoprene response to
    ! CO2 outside this range has no empirical basis. (Tai, Jan 2013)
    IF ( Inst%LISOPCO2 ) THEN
       IF ( Inst%GLOBCO2 <  150.0_hp .OR. &
            Inst%GLOBCO2 > 1250.0_hp     ) THEN
          MSG = 'Global CO2 outside valid range of 150-1250 ppmv!'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
    ENDIF

    Optfound = .FALSE.
    !optional arguments for SOAP
    CALL GetExtOpt ( HcoState%Config, ExtNr, 'Isoprene to SOAP', &
                     OptValHp=Inst%ISOPTOSOAP, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%ISOPTOSOAP = Inst%ISOPTOSOAP * 1.134
    ELSE
       !set to zero if not specified
       Inst%ISOPTOSOAP = 0.0
    ENDIF

    Optfound = .FALSE.
    !optional arguments for SOAS
    CALL GetExtOpt ( HcoState%Config, ExtNr, 'Isoprene to SOAS', &
                     OptValHp=Inst%ISOPTOSOAS, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%ISOPTOSOAS = Inst%ISOPTOSOAS * 1.134
    ELSE
       !set to zero if not specified
       Inst%ISOPTOSOAS = 0.0
    ENDIF

    Optfound = .FALSE.
    CALL GetExtOpt ( HcoState%Config, Inst%ExtNr, 'Monoterp to SOAP', &
                     OptValHp=Inst%MONOTOSOAP, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%MONOTOSOAP = Inst%MONOTOSOAP ! * 1.134
    ELSE
       !set to zero if not specified
       Inst%MONOTOSOAP = 0.0
    ENDIF

    Optfound = .FALSE.
    CALL GetExtOpt ( HcoState%Config, Inst%ExtNr, 'Monoterp to SOAS', &
                     OptValHp=Inst%MONOTOSOAS, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%MONOTOSOAS = Inst%MONOTOSOAS ! * 1.134
    ELSE
       !set to zero if not specified
       Inst%MONOTOSOAS = 0.0
    ENDIF

    Optfound = .FALSE.
    CALL GetExtOpt ( HcoState%Config, Inst%ExtNr, 'Othrterp to SOAP', &
                     OptValHp=Inst%OTHRTOSOAP, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%OTHRTOSOAP = Inst%OTHRTOSOAP * 1.134
    ELSE
       !set to zero if not specified
       Inst%OTHRTOSOAP = 0.0
    ENDIF

    Optfound = .FALSE.
    CALL GetExtOpt ( HcoState%Config, Inst%ExtNr, 'Othrterp to SOAS', &
                     OptValHp=Inst%OTHRTOSOAS, FOUND=Optfound, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Optfound ) THEN
       !convert from carbon basis to mass basis
       Inst%OTHRTOSOAS = Inst%OTHRTOSOAS * 1.134
    ELSE
       !set to zero if not specified
       Inst%OTHRTOSOAS = 0.0
    ENDIF

    !-----------------------------------------------------------------
    ! Set species IDs
    !-----------------------------------------------------------------

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = ' - Use the following species:'
       CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDIF

    ! Get species IDs
    ! --> Assume that species are ordered ISOP, ACET, PRPE, C2H4 in
    !     config. file!
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Assign species IDs
    Inst%IDTISOP = -1
    Inst%IDTACET = -1
    Inst%IDTPRPE = -1
    Inst%IDTC2H4 = -1
    Inst%IDTALD2 = -1
    Inst%IDTMOH  = -1
    Inst%IDTEOH  = -1
    Inst%IDTSOAP = -1
    Inst%IDTSOAS = -1
    Inst%IDTMTPA = -1
    Inst%IDTMTPO = -1
    Inst%IDTLIMO = -1
    Inst%IDTSESQ = -1
    DO I = 1, nSpc
       SELECT CASE ( TRIM(SpcNames(I)) )
       CASE( 'ISOP' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTISOP = HcoIDs(I)
             WRITE(MSG,*) '   Isoprene   = ',TRIM(SpcNames(I)),Inst%IDTISOP
          ELSE
             WRITE(MSG,*) '   Isoprene will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'ACET' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTACET = HcoIDs(I)
             WRITE(MSG,*) '   Acetone    = ',TRIM(SpcNames(I)),Inst%IDTACET
          ELSE
             WRITE(MSG,*) '   Acetone will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'PRPE' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTPRPE = HcoIDs(I)
             WRITE(MSG,*) '   PRPE       = ',TRIM(SpcNames(I)),Inst%IDTPRPE
          ELSE
             WRITE(MSG,*) '   PRPE will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'ALD2' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTALD2 = HcoIDs(I)
             WRITE(MSG,*) '   ALD2       = ',TRIM(SpcNames(I)),Inst%IDTALD2
          ELSE
             WRITE(MSG,*) '   ALD2 will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'MOH' )
          ! No offline biogenic VOC emissions available for MOH as
          ! of GEOS-Chem 12.7.0 (mps, 12/17/19)
          Inst%IDTMOH  = HcoIDs(I)
          WRITE(MSG,*) '   MOH        = ',TRIM(SpcNames(I)),Inst%IDTMOH
       CASE( 'EOH' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTEOH  = HcoIDs(I)
             WRITE(MSG,*) '   EOH        = ',TRIM(SpcNames(I)),Inst%IDTEOH
          ELSE
             WRITE(MSG,*) '   EOH will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'C2H4' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTC2H4 = HcoIDs(I)
             WRITE(MSG,*) '   C2H4       = ',TRIM(SpcNames(I)),Inst%IDTC2H4
          ELSE
             WRITE(MSG,*) '   C2H4 will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'SOAP' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTSOAP = HcoIDs(I)
             WRITE(MSG,*) '   SOAP       = ',TRIM(SpcNames(I)),Inst%IDTSOAP
          ELSE
             WRITE(MSG,*) '   SOAP will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'SOAS' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTSOAS = HcoIDs(I)
             WRITE(MSG,*) '   SOAS       = ',TRIM(SpcNames(I)),Inst%IDTSOAS
          ELSE
             WRITE(MSG,*) '   SOAS will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'MTPA' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTMTPA = HcoIDs(I)
             WRITE(MSG,*) '   MTPA       = ',TRIM(SpcNames(I)),Inst%IDTMTPA
          ELSE
             WRITE(MSG,*) '   MTPA will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'MTPO' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTMTPO = HcoIDs(I)
             WRITE(MSG,*) '   MTPO       = ',TRIM(SpcNames(I)),Inst%IDTMTPO
          ELSE
             WRITE(MSG,*) '   MTPO will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'LIMO' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTLIMO = HcoIDs(I)
             WRITE(MSG,*) '   LIMO       = ',TRIM(SpcNames(I)),Inst%IDTLIMO
          ELSE
             WRITE(MSG,*) '   LIMO will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE( 'SESQ' )
          IF ( .not. Inst%OFFLINE_BIOGENICVOC ) THEN
             Inst%IDTSESQ = HcoIDs(I)
             WRITE(MSG,*) '   SESQ       = ',TRIM(SpcNames(I)),Inst%IDTSESQ
          ELSE
             WRITE(MSG,*) '   SESQ will be obtained from ' // &
                          'offline biogenic VOC emissions'
          ENDIF
       CASE DEFAULT
          MSG = 'Invalid species names: ' // TRIM(SpcNames(I))
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       END SELECT

       ! Verbose
       IF ( HcoState%amIRoot ) CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDDO

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       WRITE(MSG,*) ' --> Isoprene scale factor is ',Inst%ISOP_SCALING
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Use CO2 inhibition on isoprene option ', &
                    Inst%LISOPCO2
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Global atmospheric CO2 concentration : ', &
                    Inst%GLOBCO2, ' ppmv'
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Normalize LAI by PFT: ', &
                    Inst%NORMLAI
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Isoprene to SOA-Precursor', &
                    !convert back to direct mass basis just for show
                    Inst%ISOPTOSOAP / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Isoprene direct to SOA (Simple)', &
                    !convert back to direct mass basis just for show
                    Inst%ISOPTOSOAS / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Monoterpene to SOA-Precursor', &
                    !convert back to direct mass basis just for show
                    Inst%MONOTOSOAP / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Monoterpene direct to SOA (Simple)', &
                    !convert back to direct mass basis just for show
                    Inst%MONOTOSOAS / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Othrterpene to SOA-Precursor', &
                    !convert back to direct mass basis just for show
                    Inst%OTHRTOSOAP / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
       WRITE(MSG,*) ' --> Othrterpene direct to SOA (Simple)', &
                    !convert back to direct mass basis just for show
                    Inst%OTHRTOSOAS / 1.134
       CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDIF

    CALL GetExtOpt( HcoState%Config, ExtNr, 'MEGAN_SUFFIX', &
                    OptValChar=Inst%SUFFIX, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) Inst%SUFFIX = ''

    !-----------------------------------------------------------------
    ! Allocate module variables
    !-----------------------------------------------------------------

    ! Get horizontal grid extensions on this CPU
    NX = HcoState%NX
    NY = HcoState%NY

    ALLOCATE( Inst%T_LAST24H( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'T_LAST24H', RC )
       RETURN
    ENDIF
    Inst%T_LAST24H = 0.0_hp

    ALLOCATE( Inst%T_LASTXDAYS( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'T_LASTXDAYS', RC )
       RETURN
    ENDIF
    Inst%T_LASTXDAYS = 0.0_hp

    ALLOCATE( Inst%PARDR_LASTXDAYS( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'PARDR_LASTXDAYS', RC )
       RETURN
    ENDIF
    Inst%PARDR_LASTXDAYS = 0.0_hp

    ALLOCATE( Inst%PARDF_LASTXDAYS( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'PARDF_LASTXDAYS', RC )
       RETURN
    ENDIF
    Inst%PARDF_LASTXDAYS = 0.0_hp

    ALLOCATE( Inst%LAI_PREVDAY( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'LAI_PREVDAY', RC )
      RETURN
    ENDIF
    Inst%LAI_PREVDAY = 0.0_sp

    ALLOCATE( Inst%ARRAY_16( NX, NY, 16 ), STAT=AS )
    IF ( AS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'ARRAY_16', RC )
      RETURN
    ENDIF
    Inst%ARRAY_16 = 0.0_hp

    ! Normalization factor
    ! There should be a different normalization factor for each compound, but
    ! we calculate only 1 normalization factor for all compounds (dbm 11/2012)
    ALLOCATE( Inst%NORM_FAC( 1 ), STAT=AS )
    IF ( AS /= 0 ) THEN
      CALL HCO_ERROR( HcoState%Config%Err, 'NORM_FAC', RC )
      RETURN
    ENDIF
    Inst%NORM_FAC = -99d0

    ALLOCATE( Inst%AEF_APIN( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_APIN', RC )
       RETURN
    ENDIF
    Inst%AEF_APIN = 0.0_hp

    ALLOCATE( Inst%AEF_MYRC( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_MYRC', RC )
       RETURN
    ENDIF
    Inst%AEF_MYRC = 0.0_hp

    ALLOCATE( Inst%AEF_OMON( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_OMON', RC )
       RETURN
    ENDIF
    Inst%AEF_OMON = 0.0_hp

    ALLOCATE( Inst%AEF_FARN( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_FARN', RC )
       RETURN
    ENDIF
    Inst%AEF_FARN = 0.0_hp

    ALLOCATE( Inst%AEF_BCAR( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_BCAR', RC )
       RETURN
    ENDIF
    Inst%AEF_BCAR = 0.0_hp

    ALLOCATE( Inst%AEF_OSQT( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_OSQT', RC )
       RETURN
    ENDIF
    Inst%AEF_OSQT = 0.0_hp

    ALLOCATE( Inst%AEF_MOH( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_MOH', RC )
       RETURN
    ENDIF
    Inst%AEF_MOH = 0.0_hp

    ALLOCATE( Inst%AEF_ACET( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_ACET', RC )
       RETURN
    ENDIF
    Inst%AEF_ACET = 0.0_hp

    ALLOCATE( Inst%AEF_EOH( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_EOH', RC )
       RETURN
    ENDIF
    Inst%AEF_EOH = 0.0_hp

    ALLOCATE( Inst%AEF_CH2O( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_CH2O', RC )
       RETURN
    ENDIF
    Inst%AEF_CH2O = 0.0_hp

    ALLOCATE( Inst%AEF_ALD2( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_ALD2', RC )
       RETURN
    ENDIF
    Inst%AEF_ALD2 = 0.0_hp

    ALLOCATE( Inst%AEF_FAXX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_FAXX', RC )
       RETURN
    ENDIF
    Inst%AEF_FAXX = 0.0_hp

    ALLOCATE( Inst%AEF_AAXX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_AAXX', RC )
       RETURN
    ENDIF
    Inst%AEF_AAXX = 0.0_hp

    ALLOCATE( Inst%AEF_C2H4( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err,'AEF_C2H4', RC )
       RETURN
    ENDIF
    Inst%AEF_C2H4 = 0.0_hp

    ALLOCATE( Inst%AEF_TOLU( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_TOLU', RC )
       RETURN
    ENDIF
    Inst%AEF_TOLU = 0.0_hp

    ALLOCATE( Inst%AEF_HCNX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_HCNX', RC )
       RETURN
    ENDIF
    Inst%AEF_HCNX = 0.0_hp

    ALLOCATE( Inst%AEF_PRPE( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF_PRPE', RC )
       RETURN
    ENDIF
    Inst%AEF_PRPE = 0.0_hp

    ALLOCATE( Inst%FLUXISOP( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXISOP', RC )
       RETURN
    ENDIF
    Inst%FLUXISOP = 0.0_hp

    ALLOCATE( Inst%FLUXMONO( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMONO', RC )
       RETURN
    ENDIF
    Inst%FLUXMONO = 0.0_hp

    ALLOCATE( Inst%FLUXACETmo( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXACETmo', RC )
       RETURN
    ENDIF
    Inst%FLUXACETmo = 0.0_hp

    ALLOCATE( Inst%FLUXACETmb( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err,'FLUXACETmb', RC )
       RETURN
    ENDIF
    Inst%FLUXACETmb = 0.0_hp

    ALLOCATE( Inst%FLUXACETbg( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err,'FLUXACETbg', RC )
       RETURN
    ENDIF
    Inst%FLUXACETbg = 0.0_hp

    ALLOCATE( Inst%FLUXPRPE( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXPRPE', RC )
       RETURN
    ENDIF
    Inst%FLUXPRPE = 0.0_hp

    ALLOCATE( Inst%FLUXC2H4( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXC2H4', RC )
       RETURN
    ENDIF
    Inst%FLUXC2H4 = 0.0_hp

    ALLOCATE( Inst%FLUXLIMO( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXLIMO', RC )
       RETURN
    ENDIF
    Inst%FLUXLIMO = 0.0_hp

    ALLOCATE( Inst%FLUXMTPA( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMTPA', RC )
       RETURN
    ENDIF
    Inst%FLUXMTPA = 0.0_hp

    ALLOCATE( Inst%FLUXMTPO( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMTPO', RC )
       RETURN
    ENDIF
    Inst%FLUXMTPO = 0.0_hp

    ALLOCATE( Inst%FLUXSESQ( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXSESQ', RC )
       RETURN
    ENDIF
    Inst%FLUXSESQ = 0.0_hp

    ALLOCATE( Inst%FLUXSOAP( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXSOAP', RC )
       RETURN
    ENDIF
    Inst%FLUXSOAP = 0.0_hp

    ALLOCATE( Inst%FLUXSOAS( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXSOAS', RC )
       RETURN
    ENDIF
    Inst%FLUXSOAS = 0.0_hp

    ALLOCATE( Inst%FLUXALD2( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXALD2', RC )
       RETURN
    ENDIF
    Inst%FLUXALD2 = 0.0_hp

    ALLOCATE( Inst%FLUXMOH( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMOH', RC )
       RETURN
    ENDIF
    Inst%FLUXMOH = 0.0_hp

    ALLOCATE( Inst%FLUXEOH( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXEOH', RC )
       RETURN
    ENDIF
    Inst%FLUXEOH = 0.0_hp

    ALLOCATE( Inst%FLUXAPIN( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXAPIN', RC )
       RETURN
    ENDIF
    Inst%FLUXAPIN = 0.0_hp

    ALLOCATE( Inst%FLUXBPIN( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXBPIN', RC )
       RETURN
    ENDIF
    Inst%FLUXBPIN = 0.0_hp

    ALLOCATE( Inst%FLUXSABI( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXBPIN', RC )
       RETURN
    ENDIF
    Inst%FLUXBPIN = 0.0_hp

    ALLOCATE( Inst%FLUXSABI( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXSABI', RC )
       RETURN
    ENDIF
    Inst%FLUXSABI = 0.0_hp

    ALLOCATE( Inst%FLUXMYRC( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMYRC', RC )
       RETURN
    ENDIF
    Inst%FLUXMYRC = 0.0_hp

    ALLOCATE( Inst%FLUXCARE( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXCARE', RC )
       RETURN
    ENDIF
    Inst%FLUXCARE = 0.0_hp

    ALLOCATE( Inst%FLUXOCIM( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXOCIM', RC )
       RETURN
    ENDIF
    Inst%FLUXOCIM = 0.0_hp

    ALLOCATE( Inst%FLUXOMON( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err,  'FLUXOMON', RC )
       RETURN
    ENDIF
    Inst%FLUXOMON = 0.0_hp

    ALLOCATE( Inst%FLUXFARN( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXFARN', RC )
       RETURN
    ENDIF
    Inst%FLUXFARN = 0.0_hp

    ALLOCATE( Inst%FLUXBCAR( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXBCAR', RC )
       RETURN
    ENDIF
    Inst%FLUXBCAR = 0.0_hp

    ALLOCATE( Inst%FLUXOSQT( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXOSQT', RC )
       RETURN
    ENDIF
    Inst%FLUXOSQT = 0.0_hp

    ALLOCATE( Inst%FLUXMBOX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXMBOX', RC )
       RETURN
    ENDIF
    Inst%FLUXMBOX = 0.0_hp

    ALLOCATE( Inst%FLUXFAXX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXFAXX', RC )
       RETURN
    ENDIF
    Inst%FLUXFAXX = 0.0_hp

    ALLOCATE( Inst%FLUXAAXX( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'FLUXAAXX', RC )
       RETURN
    ENDIF
    Inst%FLUXAAXX = 0.0_hp

    ALLOCATE( Inst%ARRAY_16( NX, NY, 16 ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'ARRAY_16', RC )
       RETURN
    ENDIF
    Inst%ARRAY_16 = 0.0_hp

    ALLOCATE ( Inst%AEF_ISOP ( NX, NY ), &
               Inst%AEF_MBOX ( NX, NY ), &
               Inst%AEF_BPIN ( NX, NY ), &
               Inst%AEF_CARE ( NX, NY ), &
               Inst%AEF_LIMO ( NX, NY ), &
               Inst%AEF_OCIM ( NX, NY ), &
               Inst%AEF_SABI ( NX, NY ), &
               Inst%GEIA_ORVC( NX, NY ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'AEF allocation error', RC )
       RETURN
    ENDIF
    Inst%AEF_ISOP  = 0.0_hp
    Inst%AEF_MBOX  = 0.0_hp
    Inst%AEF_BPIN  = 0.0_hp
    Inst%AEF_CARE  = 0.0_hp
    Inst%AEF_LIMO  = 0.0_hp
    Inst%AEF_OCIM  = 0.0_hp
    Inst%AEF_SABI  = 0.0_hp
    Inst%GEIA_ORVC = 0.0_hp

    !=================================================================
    ! Create manual diagnostics
    !=================================================================
    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_ACET_MONO',  &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_ACET_MBOX',  &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_ACET_DIRECT',&
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_APIN',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_BPIN',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_SABI',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_MYRC',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_CARE',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_OCIM',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_OMON',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_MONX',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_FARN',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_BCAR',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_OSQT',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_MBOX',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_FAXX',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create( HcoState  = HcoState,              &
                       cName     = 'InvMEGAN_AAXX',       &
                       ExtNr     = ExtNr,                 &
                       Cat       = -1,                    &
                       Hier      = -1,                    &
                       HcoID     = -1,                    &
                       SpaceDim  = 2,                     &
                       OutUnit   = 'kg/m2/s',             &
                       OutOper   = 'Mean',                &
                       AutoFill  = 1,                     &
                       COL       = HcoState%Diagn%HcoDiagnIDManual, &
                       RC        = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Initialize internal diagnostics. These are the restart variables
    ! that can be used for a 'warm' start of MEGAN.
    !=================================================================
    CALL HCO_RestartDefine ( HcoState, 'LAI_PREVDAY', &
                             Inst%LAI_PREVDAY, '1', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine ( HcoState, 'T_PREVDAY', &
                             Inst%T_LAST24H,   'K', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine ( HcoState, 'T_DAVG', &
                             Inst%T_LASTXDAYS, 'K', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine (  HcoState, 'PARDR_DAVG', &
                             Inst%PARDR_LASTXDAYS, 'W/m2', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_RestartDefine ( HcoState, 'PARDF_DAVG', &
                             Inst%PARDF_LASTXDAYS, 'W/m2', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! The original MEGAN code used to read the emission factors here.
    ! We now get the emisson factors through HEMCO, hence no need to
    ! do this anymore (ckeller, 05/19/2014).
    !=================================================================

    ! Set physical constants
    PI_180     = HcoState%Phys%PI_180
    Inst%D2RAD = PI_180                ! Degrees to radians
    Inst%RAD2D = 1.0_hp / PI_180       ! Radians to degrees

    ! Enable met. fields
    ExtState%T2M%DoUse       = .TRUE.
    ExtState%SUNCOS%DoUse    = .TRUE.
    ExtState%PARDR%DoUse     = .TRUE.
    ExtState%PARDF%DoUse     = .TRUE.
    ExtState%LAI%DoUse       = .TRUE.
    ExtState%GWETROOT%DoUse  = .TRUE.

    ! Cleanup
    Inst => NULL()

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

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
  SUBROUTINE HCOX_MEGAN_FINAL ( HcoState, ExtState, RC )
!
! !USES
!
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState      ! HEMCO State obj
    TYPE(Ext_State), POINTER        :: ExtState      ! Extension State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  05 Aug 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CALL InstRemove ( ExtState%Megan )

  END SUBROUTINE HCOX_MEGAN_FINAL
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate adds a new instance to the list of
!  instances, assigns a unique instance number to this new instance, and
!  archives this instance number to output argument Instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove removes an instance from the list of
! instances.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                     :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrevInst => NULL()
    Inst     => NULL()
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       IF ( ASSOCIATED( Inst%GEIA_ORVC   ) ) DEALLOCATE( Inst%GEIA_ORVC   )
       IF ( ASSOCIATED( Inst%ARRAY_16    ) ) DEALLOCATE( Inst%ARRAY_16    )
       IF ( ASSOCIATED( Inst%NORM_FAC    ) ) DEALLOCATE( Inst%NORM_FAC    )
       IF ( ASSOCIATED( Inst%AEF_ISOP    ) ) DEALLOCATE( Inst%AEF_ISOP    )
       IF ( ASSOCIATED( Inst%AEF_MBOX    ) ) DEALLOCATE( Inst%AEF_MBOX    )
       IF ( ASSOCIATED( Inst%AEF_BPIN    ) ) DEALLOCATE( Inst%AEF_BPIN    )
       IF ( ASSOCIATED( Inst%AEF_CARE    ) ) DEALLOCATE( Inst%AEF_CARE    )
       IF ( ASSOCIATED( Inst%AEF_LIMO    ) ) DEALLOCATE( Inst%AEF_LIMO    )
       IF ( ASSOCIATED( Inst%AEF_OCIM    ) ) DEALLOCATE( Inst%AEF_OCIM    )
       IF ( ASSOCIATED( Inst%AEF_SABI    ) ) DEALLOCATE( Inst%AEF_SABI    )
       IF ( ASSOCIATED( Inst%AEF_APIN    ) ) DEALLOCATE( Inst%AEF_APIN    )
       IF ( ASSOCIATED( Inst%AEF_MYRC    ) ) DEALLOCATE( Inst%AEF_MYRC    )
       IF ( ASSOCIATED( Inst%AEF_OMON    ) ) DEALLOCATE( Inst%AEF_OMON    )
       IF ( ASSOCIATED( Inst%AEF_ACET    ) ) DEALLOCATE( Inst%AEF_ACET    )
       IF ( ASSOCIATED( Inst%AEF_MOH     ) ) DEALLOCATE( Inst%AEF_MOH     )
       IF ( ASSOCIATED( Inst%AEF_EOH     ) ) DEALLOCATE( Inst%AEF_EOH     )
       IF ( ASSOCIATED( Inst%AEF_CH2O    ) ) DEALLOCATE( Inst%AEF_CH2O    )
       IF ( ASSOCIATED( Inst%AEF_ALD2    ) ) DEALLOCATE( Inst%AEF_ALD2    )
       IF ( ASSOCIATED( Inst%AEF_FAXX    ) ) DEALLOCATE( Inst%AEF_FAXX    )
       IF ( ASSOCIATED( Inst%AEF_AAXX    ) ) DEALLOCATE( Inst%AEF_AAXX    )
       IF ( ASSOCIATED( Inst%AEF_C2H4    ) ) DEALLOCATE( Inst%AEF_C2H4    )
       IF ( ASSOCIATED( Inst%AEF_TOLU    ) ) DEALLOCATE( Inst%AEF_TOLU    )
       IF ( ASSOCIATED( Inst%AEF_HCNX    ) ) DEALLOCATE( Inst%AEF_HCNX    )
       IF ( ASSOCIATED( Inst%AEF_PRPE    ) ) DEALLOCATE( Inst%AEF_PRPE    )
       IF ( ASSOCIATED( Inst%AEF_FARN    ) ) DEALLOCATE( Inst%AEF_FARN    )
       IF ( ASSOCIATED( Inst%AEF_BCAR    ) ) DEALLOCATE( Inst%AEF_BCAR    )
       IF ( ASSOCIATED( Inst%AEF_OSQT    ) ) DEALLOCATE( Inst%AEF_OSQT    )
       IF ( ASSOCIATED( Inst%FLUXISOP    ) ) DEALLOCATE( Inst%FLUXISOP    )
       IF ( ASSOCIATED( Inst%FLUXMONO    ) ) DEALLOCATE( Inst%FLUXMONO    )
       IF ( ASSOCIATED( Inst%FLUXACETmo  ) ) DEALLOCATE( Inst%FLUXACETmo  )
       IF ( ASSOCIATED( Inst%FLUXACETmb  ) ) DEALLOCATE( Inst%FLUXACETmb  )
       IF ( ASSOCIATED( Inst%FLUXACETbg  ) ) DEALLOCATE( Inst%FLUXACETbg  )
       IF ( ASSOCIATED( Inst%FLUXPRPE    ) ) DEALLOCATE( Inst%FLUXPRPE    )
       IF ( ASSOCIATED( Inst%FLUXC2H4    ) ) DEALLOCATE( Inst%FLUXC2H4    )
       IF ( ASSOCIATED( Inst%FLUXLIMO    ) ) DEALLOCATE( Inst%FLUXLIMO    )
       IF ( ASSOCIATED( Inst%FLUXMTPA    ) ) DEALLOCATE( Inst%FLUXMTPA    )
       IF ( ASSOCIATED( Inst%FLUXMTPO    ) ) DEALLOCATE( Inst%FLUXMTPO    )
       IF ( ASSOCIATED( Inst%FLUXSESQ    ) ) DEALLOCATE( Inst%FLUXSESQ    )
       IF ( ASSOCIATED( Inst%FLUXSOAP    ) ) DEALLOCATE( Inst%FLUXSOAP    )
       IF ( ASSOCIATED( Inst%FLUXSOAS    ) ) DEALLOCATE( Inst%FLUXSOAS    )
       IF ( ASSOCIATED( Inst%FLUXALD2    ) ) DEALLOCATE( Inst%FLUXALD2    )
       IF ( ASSOCIATED( Inst%FLUXMOH     ) ) DEALLOCATE( Inst%FLUXMOH     )
       IF ( ASSOCIATED( Inst%FLUXEOH     ) ) DEALLOCATE( Inst%FLUXEOH     )
       IF ( ASSOCIATED( Inst%FLUXAPIN    ) ) DEALLOCATE( Inst%FLUXAPIN    )
       IF ( ASSOCIATED( Inst%FLUXBPIN    ) ) DEALLOCATE( Inst%FLUXBPIN    )
       IF ( ASSOCIATED( Inst%FLUXSABI    ) ) DEALLOCATE( Inst%FLUXSABI    )
       IF ( ASSOCIATED( Inst%FLUXMYRC    ) ) DEALLOCATE( Inst%FLUXMYRC    )
       IF ( ASSOCIATED( Inst%FLUXCARE    ) ) DEALLOCATE( Inst%FLUXCARE    )
       IF ( ASSOCIATED( Inst%FLUXOCIM    ) ) DEALLOCATE( Inst%FLUXOCIM    )
       IF ( ASSOCIATED( Inst%FLUXOMON    ) ) DEALLOCATE( Inst%FLUXOMON    )
       IF ( ASSOCIATED( Inst%FLUXMBOX    ) ) DEALLOCATE( Inst%FLUXMBOX    )
       IF ( ASSOCIATED( Inst%FLUXFAXX    ) ) DEALLOCATE( Inst%FLUXFAXX    )
       IF ( ASSOCIATED( Inst%FLUXAAXX    ) ) DEALLOCATE( Inst%FLUXAAXX    )
       IF ( ASSOCIATED( Inst%FLUXFARN    ) ) DEALLOCATE( Inst%FLUXFARN    )
       IF ( ASSOCIATED( Inst%FLUXBCAR    ) ) DEALLOCATE( Inst%FLUXBCAR    )
       IF ( ASSOCIATED( Inst%FLUXOSQT    ) ) DEALLOCATE( Inst%FLUXOSQT    )

       IF ( ASSOCIATED( Inst%LAI_PREVDAY ) ) DEALLOCATE( Inst%LAI_PREVDAY )
       IF ( ASSOCIATED( Inst%T_LASTXDAYS ) ) DEALLOCATE( Inst%T_LASTXDAYS )
       IF ( ASSOCIATED( Inst%T_LAST24H   ) ) DEALLOCATE( Inst%T_LAST24H   )

       IF ( ASSOCIATED( Inst%PARDF_LASTXDAYS ) ) &
            DEALLOCATE( Inst%PARDF_LASTXDAYS )
       IF ( ASSOCIATED( Inst%PARDR_LASTXDAYS ) ) &
            DEALLOCATE( Inst%PARDR_LASTXDAYS )

       ! ----------------------------------------------------------------
       ! Pop off instance from list
       ! ----------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()

    ENDIF

  END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_MEGAN_MOD
