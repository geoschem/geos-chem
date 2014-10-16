!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
! initialization of MEGAN. This is not feasible in an ESMF environment. 
! As a first approximation, the averages are now set to the initial 
! temperature/radiation value. Especially for the radiation (PARDR and 
! PARDF), this may result in very different average values - which in 
! turn may result in local differences of >50%. HEMCO emissions are 
! identical to original GEOS-Chem emissions once we pass the first 10 
! days (ckeller, 10/8/2013)! 
!\\
!\\
! A similar procedure is now also applied for the leaf area index variables.
! The original GEOS-Chem MEGAN code used three LAI variables: current month
! LAI (LAI\_CM), previous month LAI (LAI\_PM), and instantaneous LAI (LAI, 
! which was a daily interpolation of LAI\_CM and next month' LAI, LAI_NM). 
! The HEMCO implementation uses only the instantaneous LAI any more, assuming 
! it changes every day. The previous day LAI is kept in memory and used to
! determine the LAI change over time (used to calculate the gamma leaf age).
! For the first simulation day, a LAI change of zero is assumed.
! (ckeller, 10/9/2014)
!\\
!\\
! !References:
!
!  \begin{itemize}
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
      PRIVATE :: GET_EMACET_MEGAN    
      PRIVATE :: GET_EMISOP_MEGAN  
      PRIVATE :: GET_EMMBO_MEGAN    
!     SOAudpate: renamed, includes SESQ (hotp 3/2/10)
      PRIVATE :: GET_EMTERP_MEGAN
      PRIVATE :: GET_EMMONOT_MEGAN 
      PRIVATE :: UPDATE_HCOT_DAY      
      PRIVATE :: UPDATE_HCOT_15_AVG   
      PRIVATE :: GET_GAMMA_LAI
      PRIVATE :: GET_GAMMA_LEAF_AGE
      PRIVATE :: GET_GAMMA_P
      PRIVATE :: GET_GAMMA_T_ISOP
      PRIVATE :: GET_GAMMA_T_NISOP
      PRIVATE :: GET_GAMMA_P_PECCA
      PRIVATE :: SOLAR_ANGLE
      PRIVATE :: GET_AEF_GEN
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
! 
      INTEGER                 :: ExtNr
      INTEGER                 :: ExtNrMono 
      INTEGER                 :: ExtNrSoa
      LOGICAL                 :: LPECCA
      REAL*8                  :: ISOP_SCALING

      ! Physical parameter
      REAL*8                  :: D2RAD  ! Degrees to radians 
      REAL*8                  :: RAD2D  ! Radians to degrees

      ! HEMCO species IDs
      INTEGER                 :: IDTISOP
      INTEGER                 :: IDTMONX
      INTEGER                 :: IDTACET
      INTEGER                 :: IDTCO 
      INTEGER                 :: IDTPRPE
      INTEGER                 :: IDTC2H4
      INTEGER                 :: IDTOCPI
      INTEGER                 :: IDTLIMO
      INTEGER                 :: IDTMTPA
      INTEGER                 :: IDTMTPO
      INTEGER                 :: IDTSESQ

      ! Pointers to annual emission factor arrays.
      ! These fields are obtained from ext. data (from config file)
      ! These arrays represent the base emission values in kg(C)/m2/s
      ! that will be scaled based upon meteorological conditions.
      REAL(hp), POINTER       :: AEF_ISOP(:,:) => NULL()     ! Isoprene 
      REAL(hp), POINTER       :: AEF_MONOT(:,:) => NULL()    ! Total monoterpenes
      REAL(hp), POINTER       :: AEF_MBO(:,:) => NULL()      ! Methyl butenol
      REAL(hp), POINTER       :: AEF_APINE(:,:) => NULL()    ! Alpha-pinene
      REAL(hp), POINTER       :: AEF_BPINE(:,:) => NULL()    ! Beta-pinene
      REAL(hp), POINTER       :: AEF_LIMON(:,:) => NULL()    ! Limonene
      REAL(hp), POINTER       :: AEF_SABIN(:,:) => NULL()    ! Sabine
      REAL(hp), POINTER       :: AEF_MYRCN(:,:) => NULL()    ! Myrcene
      REAL(hp), POINTER       :: AEF_CAREN(:,:) => NULL()    ! 3-Carene
      REAL(hp), POINTER       :: AEF_OCIMN(:,:) => NULL()    ! Ocimene
      REAL(hp), POINTER       :: AEF_ACET(:,:) => NULL()     ! Acetone
      REAL(hp), POINTER       :: GEIA_ORVC(:,:) => NULL()

      ! Generic emission factors (all 20 MEGAN groups). Will be determined from 
      ! plant functional types.
      REAL(hp), ALLOCATABLE   :: AEF_GEN(:,:,:)

      ! Scalars
      INTEGER, PARAMETER  :: DAY_DIM        = 24       ! # of 1-hr periods/day
!      INTEGER, PARAMETER  :: DAY_DIM        = 8        ! # of 3-hr periods/day
      INTEGER, PARAMETER  :: NUM_DAYS       = 10       ! # of days to avg 
      REAL*8,  PARAMETER  :: WM2_TO_UMOLM2S = 4.766d0  ! W/m2 -> umol/m2/s

      ! SOAupdate: Megan group IDs (from MEGAN CDP) (hotp 3/1/10)
      ! Used to locate species within EF and AEF_GEN arrays
      INTEGER, PARAMETER  :: IDMGFARN = 10  ! farnesene
      INTEGER, PARAMETER  :: IDMGBCAR = 11  ! beta-caryophyllene
      INTEGER, PARAMETER  :: IDMGOSQT = 12  ! other sesquiterpenes
      INTEGER, PARAMETER  :: IDMGOMTP = 9   ! other monoterpenes

      ! Past light & temperature conditions (mpb,2009)
      ! (1) Temperature at 2m (TS):
      REAL*8, ALLOCATABLE :: HCOT_DAILY(:,:)       ! Daily averaged sfc temp
      REAL*8, ALLOCATABLE :: HCOT_DAY(:,:,:)       ! Holds 1 day of sfc temp data
      REAL*8, ALLOCATABLE :: HCOT_15(:,:,:)        ! Holds 15 days of daily avg T
      REAL*8, ALLOCATABLE :: HCOT_15_AVG(:,:)      ! Sfc temp avg'd over NUM_DAYS

      ! (2) PAR Direct: 
      REAL*8, ALLOCATABLE :: HCOPARDR_DAILY(:,:)   ! Average daily PARDR
      REAL*8, ALLOCATABLE :: HCOPARDR_DAY(:,:,:)   ! Holds 1 day of PARDR data
      REAL*8, ALLOCATABLE :: HCOPARDR_15(:,:,:)    ! 10 days of daily avg'd PARDR
      REAL*8, ALLOCATABLE :: HCOPARDR_15_AVG(:,:)  ! PARDR averaged over NUM_DAYS

      ! (3) PAR Diffuse: 
      REAL*8, ALLOCATABLE :: HCOPARDF_DAILY(:,:)   ! Average daily PARDR
      REAL*8, ALLOCATABLE :: HCOPARDF_DAY(:,:,:)   ! Holds 1-day of PARDR data
      REAL*8, ALLOCATABLE :: HCOPARDF_15(:,:,:)    ! 10 days of daily avg'd PARDR 
      REAL*8, ALLOCATABLE :: HCOPARDF_15_AVG(:,:)  ! PARDF averaged over NUM_DAYS

      ! previous LAI values (ckeller, 10/09/2014)
      REAL*8, ALLOCATABLE :: LAI_PREVSTEP(:,:)     ! LAI of prev. time step
      REAL*8, ALLOCATABLE :: LAI_PREVDAY(:,:)      ! LAI of prev. day

      ! Days between mid-months (updated by HEMCO clock)
      INTEGER             :: DAYS_BTW_M

      ! testing only
      integer, parameter  :: ix = -1 !25 !13 !20
      integer, parameter  :: iy = -1 !22 !38 !31

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewHour 
      USE HCO_CLOCK_MOD,       ONLY : HcoClock_NewDay
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
! !REVISION HISTORY: 
!  05 Aug 2013 - C. Keller - Initial version 
!  08 Oct 2014 - C. Keller - Now use HEMCO clock to get days between months
!  09 Oct 2014 - C. Keller - Now use only GC_LAI. This makes days between months
!                            obsolete
!  16 Oct 2014 - C. Keller - Initialize flux arrays to avoid float invalid
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J
      REAL*8              :: SUNCOS,   D_BTW_M,  Q_DIR_2
      REAL*8              :: Q_DIFF_2, ISOLAI,   MISOLAI,  PMISOLAI
      REAL*8              :: GAMMA_P,  GAMMA_SM, GAMMA_LAI
      REAL*8              :: EMISOP, EMMO,  EMMB,  APINE, EMMT
      REAL*8              :: BPINE,  LIMON, SABIN, MYRCN, CAREN, OCIMN 
      REAL*8              :: FARN,   BCAR,  OSQT,  OMTP
      REAL*8              :: TERP_ORGC, BIOG_ALPH, BIOG_LIMO
      REAL*8              :: BIOG_ALCO
      REAL*8              :: EMMO_MOL
      REAL*8              :: MB_MOL
      REAL*8              :: X, Y 
      REAL*8              :: DIURFACT, DIURORVC
      REAL*8              :: TMP
      LOGICAL             :: ERR
      LOGICAL             :: FIRST, IsNewHour, IsNewDay 

      ! Emission arrays (kgC/m2/s)
      REAL(hp), TARGET    :: FLUXISOP  (HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: FLUXEMMO  (HcoState%NX,HcoState%NY)
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

      ! For diagnostics
      REAL(hp), POINTER   :: Arr2D(:,:) => NULL()
      CHARACTER(LEN=63)   :: DiagnName

      ! Molecules C / kg C
      REAL*8, PARAMETER   :: XNUMOL_C = 6.022d+23 / 12d-3
      REAL*8, PARAMETER   :: UNITCONV = 1.0d0

      ! Fraction of yield of OC (hydrophilic) from terpene emissions
      REAL*8, PARAMETER   :: EMMT2OC = 1.0d-1

      ! Conversion factors for OC/BC
      REAL*8, PARAMETER   :: FC1 = 136.2364d0 / 120.11d0
      REAL*8, PARAMETER   :: FC2 = 154.2516d0 / 120.11d0
      REAL*8, PARAMETER   :: FC3 = 204.3546d0 / 180.165d0
      REAL*8, PARAMETER   :: FC4 = 152.d0     / 120.11d0

      ! Conversion factors for acetone calculations
      REAL*8, PARAMETER   :: YIELD_MO   = 0.116d0
      REAL*8, PARAMETER   :: MONO_SCALE = 0.89d0
      REAL*8, PARAMETER   :: MB_SCALE1  = 0.6d0
      REAL*8, PARAMETER   :: MB_SCALE2  = 0.76d0

      !--------------------------------------------------------------------
      !  Biogenic source of PRPE -- scaled to ISOPRENE
      !
      !  Also, add biogenic emissions of alkenes. We do this by scaling to
      !  isoprene emissions (probably OK for summertime conditions). The
      !  scaling factor is based on work by Allen Goldstein. His values
      !  indicate emission ratios of ethene:propene:butene=4:2:1 (on a
      !  per molecule basis), with total emissions approx. equal to
      !  10% of isoprene emissions (again, on molecule basis).
      !  BIOSCAL is in units of atoms C (alkenes) / atoms C (isoprene)
      !*******************************************************************
      ! Change this factor to exclude ethene (bey, ljm)
      !    (10 molec alkenes / 100 molec isop) * (1 molec isop / 5 atoms C isop)
      !    *(3 molec butene + propene / 7 molec total alkenes)
      !    *(3.3333 atoms C but+prop mix/ 1 molec but+prop mix)
      !    = 0.0286 atoms C butene+propene / atom C isop
      ! Note that 3.3333 atoms C/molecule is the weighted average for this mix.
      ! Scale factor is the same on a mass basis, i.e. kgC/kgC 
      ! (ckeller, 08/05/13) 
      !--------------------------------------------------------------------
      REAL*8, PARAMETER   :: BIOSCAL = 0.0286d0 ! new factor, (ljm, bey, 9/28/98)

      !--------------------------------------------------------------------
      ! Biogenic source of ethene (scaled to isoprene)
      ! Scale factor BIOSCALEC2H4 =
      !   ( 10 molec alkenes / 100 molec isop ) * ( 1 molec isop / 5 atoms C )
      !   * ( 4 molec ethene / 7 molec alkenes )
      !   * ( 2 atoms C / 1 molec ethene )
      !   = 0.022857d0  [atoms C / atoms C isop]
      ! (tmf, 1/13/06)
      ! Scale factor is the same on a mass basis, i.e. kgC/kgC 
      ! (ckeller, 08/05/13) 
      !--------------------------------------------------------------------
      REAL*8, PARAMETER   :: BIOSCALEC2H4 = 0.022857d0

      !=================================================================
      ! HCOX_Megan_Run begins here! 
      !================================================================= 

      ! Enter 
      CALL HCO_ENTER ( 'HCOX_Megan_Run (hcox_megan_mod.F)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      ERR = .FALSE.

      ! Time information
      IsNewHour = HcoClock_NewHour()
      IsNewDay  = HcoClock_NewDay()
      FIRST     = HcoClock_First()

      ! Initialize flux arrays
      FLUXISOP   = 0.0_hp
      FLUXEMMO   = 0.0_hp
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
      FLUXSESQ   = 0.0_hp

      !-----------------------------------------------------------------
      ! Point to external data 
      !-----------------------------------------------------------------

      ! Isoprene and acetone are always calculated
      IF ( FIRST ) THEN
         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_ISOP', 
     &                             AEF_ISOP,   RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_ACET', 
     &                             AEF_ACET,   RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         ! Arrays used if monoterpene switch is turned on
            IF ( ExtNrMono > 0 ) THEN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_MONOT', 
     &                                AEF_MONOT,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_MBO', 
     &                                AEF_MBO,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
      
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_APINE', 
     &                                AEF_APINE,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_BPINE', 
     &                                AEF_BPINE,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_LIMON', 
     &                                AEF_LIMON,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_SABIN', 
     &                                AEF_SABIN,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_MYRCN', 
     &                                AEF_MYRCN,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_CAREN', 
     &                                AEF_CAREN,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_AEF_OCIMN', 
     &                                AEF_OCIMN,   RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
         ENDIF
   
         ! Only used for SOA simulation:
         IF ( ExtNrSoa>0 ) THEN
            CALL HCO_GetPtr( am_I_Root, 'MEGAN_ORVC', 
     &                                GEIA_ORVC, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
         ENDIF
      ENDIF ! FIRST

      ! SET AEF_GEN (only on first call). These values are needed for
      ! subroutine GET_EMMONOT_MEGAN, which is only called for OCPI
      IF ( FIRST .AND. (ExtNrMono>0) .AND. (ExtNrSoa<=0) .AND. 
     &     (IDTOCPI>0) ) THEN
         CALL GET_AEF_GEN( am_I_Root, HcoState, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! On first call, fill previous day LAI with current LAI since we
      ! don't know the previous day values yet... 
      IF ( FIRST ) THEN
         LAI_PREVDAY = ExtState%GC_LAI%Arr%Val

      ! On a new day, copy previous LAI values to array that holds
      ! previous day LAI.
      ELSEIF ( IsNewDay ) THEN
         LAI_PREVDAY = LAI_PREVSTEP
      ENDIF

      ! For testing, set to 1 day since we know that LAI is updated every
      ! day.
      DAYS_BTW_M = 1

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, SUNCOS, Q_DIR_2, Q_DIFF_2, ISOLAI, PMISOLAI       )
!$OMP+PRIVATE( MISOLAI, GAMMA_P, GAMMA_LAI, GAMMA_SM, EMISOP           ) 
!$OMP+PRIVATE( EMMB, APINE, BPINE, LIMON, SABIN, MYRCN, CAREN, OCIMN   ) 
!$OMP+PRIVATE( EMMO, EMMO_MOL, TMP, MB_MOL, TERP_ORGC                  ) 
!$OMP+PRIVATE( DIURFACT, DIURORVC, BIOG_ALPH                           ) 
!$OMP+PRIVATE( BIOG_LIMO, BIOG_ALCO, X, Y, RC                          )

      !-----------------------------------------------------------------
      ! Loop over all grid boxes 
      !-----------------------------------------------------------------
      DO J = 1, HcoState%NY 
      DO I = 1, HcoState%NX

         !-----------------------------------------------------------------
         ! Update internal temperature arrays 
         !-----------------------------------------------------------------
   
         ! Arrays holding hourly temperature bins. Update every hour. 
         ! For now, all historical values are set to the current values
         ! on first call. This avoids reading all met-fields of the past
         ! 10 days, but introduces errors in the MEGAN emissions of the
         ! first 10 days.
         IF ( IsNewHour .OR. FIRST ) THEN 
            CALL UPDATE_HCOT_DAY ( HcoState, ExtState, I, J, FIRST, RC )
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
         IF ( IsNewDay ) THEN
            CALL UPDATE_HCOT_15_AVG ( HcoState, I, J, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
               CALL HCO_ERROR( 'UPDATE_HCOT_15_AVG', RC )
               ERR = .TRUE.
               EXIT
            ENDIF
         ENDIF

         !-----------------------------------------------------------------
         ! Extract meteorological variables 
         !-----------------------------------------------------------------

         ! Cosine of solar zenith angle
         SUNCOS   = ExtState%SUNCOSmid%Arr%Val(I,J)
      
         ! Convert Q_DIR and Q_DIFF from (W/m2) to (micromol/m2/s)
         Q_DIR_2  = ExtState%PARDR%Arr%Val(I,J) * WM2_TO_UMOLM2S
         Q_DIFF_2 = ExtState%PARDF%Arr%Val(I,J) * WM2_TO_UMOLM2S

         ! Pass met variables. Now use only GC_LAI (ckeller, 10/9/2014)
         ISOLAI   = ExtState%GC_LAI%Arr%Val(I,J)
         PMISOLAI = LAI_PREVDAY(I,J)
         MISOLAI  = ISOLAI
!         PMISOLAI = ExtState%GC_LAI_PM%Arr%Val(I,J)
!         MISOLAI  = ExtState%GC_LAI_CM%Arr%Val(I,J)

         !-----------------------------------------------------------------
         ! Get species independent scale factors for particular grid box 
         !-----------------------------------------------------------------

         !---------------------------------------------------
         ! Do only during day and over continent
         ! Only interested in terrestrial biosphere (pip)
         !---------------------------------------------------
         IF ( SUNCOS > 0d0 ) THEN

            IF ( LPECCA ) THEN 
                
               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P_PECCA( HcoState, ExtState, I, J, 
     &                   Q_DIR_2, Q_DIFF_2, HCOPARDR_15_AVG(I,J), 
     &                   HCOPARDF_15_AVG(I,J))

            ELSE 

               ! Activity factor for light 
               GAMMA_P = GET_GAMMA_P(  ISOLAI,  SUNCOS,
     &                                 Q_DIR_2, Q_DIFF_2 )

            ENDIF
      
         ELSE
            
            GAMMA_P = 0.0

         ENDIF

         ! Activity factor for leaf area
         GAMMA_LAI = GET_GAMMA_LAI( MISOLAI )

         ! Activity factor for soil moisture
         GAMMA_SM = 1.d0

         !--------------------------------------------------------------
         ! Calculate VOC emissions 
         ! The GET_EMIS*_MEGAN calls now use the annual scale factors
         ! imported through the HEMCO list, which are already in units 
         ! of kg(C)/m2/s. Thus, no further unit conversion is required 
         ! anymore (--> UNITCONV = 1). ckeller, 14/01/25.
         !--------------------------------------------------------------

!         ! testing only
!         if ( i==ix .and. j==iy ) then
!                  write(*,*) ' '
!                  write(*,*) 'do HEMCO MEGAN @i,j ', i,j
!                  write(*,*) 'SUNCOSmid:             ', SUNCOS
!                  write(*,*) 'TEMPERATURE:',ExtState%TSURFK%Arr%Val(I,J)
!                  write(*,*) 'Q_DIR: ', Q_DIR_2 / WM2_TO_UMOLM2S
!                  write(*,*) 'Q_DIFF:', Q_DIFF_2 / WM2_TO_UMOLM2S
!                  write(*,*) 'HCOT_15_AVG: ', HCOT_15_AVG(I,J)
!                  write(*,*) 'HCOPARDR_15_AVG: ', HCOPARDR_15_AVG(I,J)
!                  write(*,*) 'HCOPARDF_15_AVG: ', HCOPARDF_15_AVG(I,J)
!                  write(*,*) 'HCOT_DAILY : ', HCOT_DAILY (I,J)
!                  write(*,*) ' '
!         endif
!
         !--------------------------------------------------------------------
         ! MEGAN Isoprene
         !--------------------------------------------------------------------

         ! Isoprene [kg C/m2/s]
         CALL GET_EMISOP_MEGAN( HcoState, ExtState, I, J,
     &                          GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                          UNITCONV, EMISOP, RC )
         IF ( RC /= HCO_SUCCESS ) THEN 
            CALL HCO_ERROR( 'GET_EMISOP_MEGAN', RC )
            ERR = .TRUE.
            EXIT
         ENDIF 

         !FP_ISOP (6/2009)
         EMISOP = ISOP_SCALING * EMISOP

         ! Write isoprene emissions into tracer tendency array
         IF ( IDTISOP > 0 ) THEN
            FLUXISOP(I,J) = EMISOP 

            ! testing only
            if ( i==ix .and. j==iy ) then
               write(*,*) 'HEMCO MEGAN @ ', ix, iy
               write(*,*) 'ISOP_SCAL: ', ISOP_SCALING
               write(*,*) 'ISOP (kgC/m2/s): ', EMISOP
            endif

         ENDIF

         !--------------------------------------------------------------------
         ! MEGAN biogenics
         !--------------------------------------------------------------------

         IF ( ExtNrMono > 0 ) THEN

            ! ---------------------------------------------------
            ! Methyl butenol
            CALL GET_EMMBO_MEGAN( HcoState, ExtState, I, J, 
     &                            GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                            UNITCONV, EMMB, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMMBO_MEGAN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Alpha Pinene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'APINE', APINE, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_APINE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Beta Pinene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'BPINE', BPINE, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_BPINE', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Limonene emissions 
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'LIMON', LIMON, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_LIMON', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Sabinene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'SABIN', SABIN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_SABIN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Mycrene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'MYRCN', MYRCN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_MYRCN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! 3-Carene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'CAREN', CAREN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_CAREN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Ocimene emissions
            CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV , 'OCIMN', OCIMN, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMTERP_MEGAN_OCIMN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! ---------------------------------------------------
            ! Total monoterpenes = sum of individual
            EMMO = APINE + BPINE + LIMON + SABIN + 
     &             MYRCN + CAREN + OCIMN

            ! Add to tracer tendency array [kg C/m2/s]
            IF ( IDTMONX > 0 ) THEN
               FLUXEMMO(I,J) = EMMO
            ENDIF

         ELSE
            EMMO = 0d0
            EMMB = 0d0
         ENDIF

! -------------------------------------------------------------
!         EMISS_BVOC( I , J , 1 ) =  EMIS   / AREA_CM2 / DTSRCE     
!         EMISS_BVOC( I , J , 2 ) =  EMMO   / AREA_CM2 / DTSRCE
!         EMISS_BVOC( I , J , 3 ) =  EMMB   / AREA_CM2 / DTSRCE                 
!         EMISS_BVOC( I , J , 4 ) =  APINE  / AREA_CM2 / DTSRCE     
!         EMISS_BVOC( I , J , 5 ) =  BPINE  / AREA_CM2 / DTSRCE
!         EMISS_BVOC( I , J , 6 ) =  LIMON  / AREA_CM2 / DTSRCE  
!         EMISS_BVOC( I , J , 7 ) =  SABIN  / AREA_CM2 / DTSRCE     
!         EMISS_BVOC( I , J , 8 ) =  MYRCN  / AREA_CM2 / DTSRCE
!         EMISS_BVOC( I , J , 9 ) =  CAREN  / AREA_CM2 / DTSRCE  
!         EMISS_BVOC( I , J , 10) =  OCIMN  / AREA_CM2 / DTSRCE 
! -------------------------------------------------------------

          ! testing only
          if ( i == ix .AND. j == iy ) then
                  write(*,*) ' '
                  write(*,*) 'HEMCO MEGAN (kg/m2/s) @ ', i, j 
                  write(*,*) 'EMISOP: ', EMISOP
                  write(*,*) 'EMMO: ', EMMO
                  write(*,*) 'EMMB: ', EMMB
                  write(*,*) 'APINE: ', APINE
                  write(*,*) 'BPINE: ', BPINE
                  write(*,*) 'LIMON: ', LIMON
                  write(*,*) 'SABIN: ', SABIN
                  write(*,*) 'MYRCN: ', MYRCN
                  write(*,*) 'CAREN: ', CAREN
                  write(*,*) 'OCIMN: ', OCIMN
          endif

!-----------------------------------------------------------------------------
         

         !--------------------------------------------------------------------
         ! BIOGENIC ACETONE EMISSIONS
         !--------------------------------------------------------------------
         IF ( IDTACET > 0 ) THEN

            !-----------------------------------------------------------------
            ! (1) BIOGENIC EMISSIONS OF ACETONE FROM MONOTERPENES
            ! Monoterpenes has same # molecules/kg of carbon as isoprene
            ! The yield for monoterpenes is .12 mol/mol from Reisell et.al. 
            ! 1999 (this does not includes direct acetone emissions)

            ! Convert [kg C/m2/s] to [kg MONOTERPENE/m2/s]
            ! There are 10 C atoms per molecule of MONOTERPENE
            EMMO_MOL   = EMMO / 10d0
      
            ! Apply yield from monoterpenes to get [kg ACET/m2/s]
            TMP        = EMMO_MOL * YIELD_MO
      
            ! Convert acetone emissions back into [kg C/m2/s] 
            TMP       = TMP * 3.d0
      
            ! Scale to a posteriori source from Jacob et al 2001 (bdf, 9/5/01)
            TMP        = TMP * MONO_SCALE
      
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
            TMP = 0D0
      
            ! Lon and lat of grid box (I,J) in degrees
            X = HcoState%Grid%XMID%Val( I, J )
            Y = HcoState%Grid%YMID%Val( I, J )
      
            ! Methyl butenol is emitted only in North America, where
            ! ( -167.5 <= lon <= -52.5 ) and ( 16.0 <= lat <= 72.0 ) 
            IF ( ( X >= -167.5d0 .and. X <= -52.5d0 ) .AND.
     &           ( Y >=   16.0d0 .and. Y <=  72.0d0 ) ) THEN
      
               ! Convert from [kg C/m2/s] to [kg MBO/m2/s] 
               ! There are 5 C atoms per molecule MBO
               MB_MOL   = EMMB / 5.d0
      
               ! Apply yield from MBO to get [kg ACET/m2/s]
               TMP = MB_MOL * MB_SCALE1
      
               ! Convert from [kg ACET/m2/s] to [kg C/m2/s]
               ! There are 3 C atoms per acetone molecule
               TMP = TMP * 3.d0
      
               ! Scale to a posteriori source from Jacob et al 2001 (bdf,
               ! 9/5/01)
               FLUXACETmb(I,J) = TMP * MB_SCALE2

            ENDIF

            !-----------------------------------------------------------------
            ! (3) BIOGENIC ACETONE -- DIRECT EMISSION 
            ! evf, removed obsolete code, replaced with MEGAN acetone
            ! emissions (5/25/2011) Direct Emission now includes emission
            ! from grasses and emission from dry leaf matter
            CALL GET_EMACET_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             UNITCONV, TMP, RC )
            IF ( RC /= HCO_SUCCESS ) THEN 
               CALL HCO_ERROR( 'GET_EMACET_MEGAN', RC )
               ERR = .TRUE.
               EXIT
            ENDIF 

            ! Write to array
            FLUXACETbg(I,J) = TMP
         ENDIF

         !--------------------------------------------------------------------
         !  Biogenic source of CO -- from oxidation of METHANOL and MONOTERPENES
         !
         !  CO from METHANOL oxidation -- scaled from ISOPRENE (bnd, 1/2/01)
         !
         !    We need to scale the Isoprene flux to get the CH3OH (methanol) flux.
         !    Currently, the annual isoprene flux in GEOS-CHEM is ~ 397 Tg C.
         !
         !    Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH oxidation  
         !    based on Singh et al. 2000 [JGR 105, 3795-3805] who estimate a global 
         !    methanol source of 122 Tg yr-1, of which most (75 Tg yr-1) is 
         !    "primary biogenic".  He also recommends for now that the CO flux 
         !    from CH3OH oxidation be scaled to monthly mean isoprene flux.
         !
         !    To get CO from METHANOL oxidation, we must therefore multiply
         !    the ISOPRENE flux by the following scale factor:
         !      ( 100 Tg CO from CH3OH Oxidation  / 397 Tg C from Isoprene Flux ) *
         !      (  12 g C/mole                    / 28 g CO/mole                )
         !
         !  CO from MONOTERPENE oxidation (bnd, bmy, 1/2/01)
         !
         !    Assume the production of CO from monoterpenes is instantaneous even 
         !    though the lifetime of intermediate species may be on the order of hours 
         !    or days.  This assumption will likely cause CO from monoterpene oxidation
         !    to be too high in the box in which the monoterpene is emitted.
         !
         !    The CO yield here is taken from:
         !
         !    Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
         !      "The ultimate yield of CO from the tropospheric oxidation of terpenes 
         !       (including both O3 and OH reactions) was estimated to be 20% on the 
         !       carbon number basis."  They studied ALPHA- & BETA-pinene.
         !
         !    Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
         !      "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
         !
         !    EMMO is now kgC, conversion factor then becomes: 
         !    28g CO / 60g C = 46.666%
         !    (ckeller, 4/15/14).
         !--------------------------------------------------------------------
 
         ! NOTE: Don't save into EMISRR for the tagged CO simulation.  
         ! Also for tagged CO we don't use monoterpenes  ??????
         ! (jaf, mak, bmy, 2/14/08)
         IF ( (IDTCO>0) .AND. (ExtNrMono>0) ) THEN
!            FLUX = EMMO * 0.2d0
            FLUXCO(I,J) = EMMO * 0.4666d0

            ! testing only
            if ( i==ix .and. j==iy ) then
               write(*,*) ' '
               write(*,*) 'HEMCO MEGAN @ ', ix, iy
               write(*,*) 'CO (kg/m2/s): ', FLUXCO(I,J)
            endif

!            ! ND29: CO-source from monoterpenes [molec/cm2/s]
!            IF ( ND29 > 0 ) THEN
!               AD29(I,J,5) = AD29(I,J,5) + ( TMPVAL / AREA_CM2 )
!            ENDIF
         ENDIF

         ! NOTE: Don't save into EMISRR for the tagged 
         ! CO simulation. (jaf, mak, bmy, 2/14/08)
         IF ( IDTPRPE > 0 ) THEN
            FLUXPRPE(I,J) = EMISOP * BIOSCAL

            ! testing only
            if ( i==ix .and. j==iy ) then
               write(*,*) ' '
               write(*,*) 'HEMCO MEGAN @ ', ix, iy
               write(*,*) 'PRPE (kg/m2/s): ', FLUXPRPE(I,J)
            endif

         ENDIF

         !--------------------------------------------------------------------
         ! Add biogenic emission of ethene (C2H4) --> scaled to isoprene
         !--------------------------------------------------------------------

         IF ( IDTC2H4 > 0 ) THEN
            FLUXC2H4(I,J) = EMISOP * BIOSCALEC2H4
         ENDIF

         !--------------------------------------------------------------------
         ! Biogenic emissions of OC
         ! Assume a 10% yield of hydrophilic OC from terpene emissions
         ! if secondary organic aerosols are turned off. Otherwise,
         ! use CALTECH algorithm to calculate biogenic VOC emissions.
         !--------------------------------------------------------------------

         ! IF secondary organic aerosol is turned off...
         IF ( ExtNrSoa <= 0 ) THEN

            ! Get monoterpenes from MEGAN [kg C/m2/s]
            ! terpene emissions
            IF ( (ExtNrMono > 0 ) .AND. ( IDTOCPI > 0 ) ) THEN 

               CALL GET_EMMONOT_MEGAN( HcoState, ExtState, I, J, 
     &                                 GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                                 UNITCONV, TERP_ORGC, RC )
               IF ( RC /= HCO_SUCCESS ) THEN 
                  CALL HCO_ERROR( 'GET_EMMONOT_MEGAN', RC )
                  ERR = .TRUE.
                  EXIT
               ENDIF               

               ! Convert monoterpenes into hydrophilic OC and pass to
               ! OC emissions array 
               FLUXOCPI(I,J) = TERP_ORGC * EMMT2OC

               ! testing only
               if ( i==ix .and. j==iy ) then
                  write(*,*) ' '
                  write(*,*) 'HEMCO MEGAN @ ', ix, iy
                  write(*,*) 'OCPI (NOSOA) (kg/m2/s): ', FLUXOCPI(I,J)
               endif

            ELSE
               FLUXOCPI(I,J) = 0.0d0 
            ENDIF

         ! Use CALTECH algorithm if SOA is turned on
         ELSE
           
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
               ! Other monoterpenes (OMTP from MEGAN) all placed in TERP
               ! All sesq lumped together into SESQ

               ! new mtp lumping (as of 5/20/10) hotp
               ! GEOS-Chem   MEGAN
               ! =========   ==========================================
               ! MTPA        a-pinene, b-pinene, sabinene, carene
               ! LIMO        limonene
               ! MTPO        myrcene, ocimene, other monot (OMTP)
               ! SESQ        farnesene (FARN), b-caryoph (BCAR), OSQT
               ! =========   ==========================================

               !----------------------------------------
               ! MTPA 
               !----------------------------------------
               ! rename EMMONOG to EMTERP (hotp 7/29/10)
               ! New lumping: MTPA=a-,b-pinene,sabinene,carene (hotp 5/20/10)
               ! a-pinene
               IF ( IDTMTPA > 0 ) THEN
                  FLUXMTPA(I,J) = (APINE + BPINE + SABIN + CAREN) * FC1
               ENDIF

               !----------------------------------------
               ! LIMO 
               !----------------------------------------
               ! LIMONENE [kg C/m2/s]
               IF ( IDTLIMO > 0 ) THEN
                  FLUXLIMO(I,J) = LIMON * FC1 
               ENDIF

               !----------------------------------------
               ! MTPO
               !----------------------------------------
               ! MTPO is all other monoterpenes (MEGAN categories:
               ! myrcene, ocimene, OMTP) (hotp 5/20/10)
               ! All other monoterpenes (mostly camphene, linalool,
               ! terpinolene, terpinolene, phellandrene) (hotp 3/10/10)
               ! 14-18% of OMTP is terpinene and terpinolene
               IF ( IDTMTPO > 0 ) THEN

                  ! ---------------------------------------------------
                  ! Other monoterpene
                  CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                                   GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                                   UNITCONV , 'OMTPE', OMTP, RC )
                  IF ( RC /= HCO_SUCCESS ) THEN 
                     CALL HCO_ERROR( 'GET_EMTERP_MEGAN_FARNE', RC )
                     ERR = .TRUE.
                     EXIT
                  ENDIF 

                  FLUXMTPO(I,J) = (OMTP + MYRCN + OCIMN ) * FC1
               ENDIF

               !----------------------------------------
               ! SESQ updated
               !----------------------------------------
               ! Sesquiterpenes from MEGAN (hotp 3/10/10)
               IF ( IDTSESQ > 0 ) THEN

                  ! ---------------------------------------------------
                  ! a-Farnesene 
                  CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                                   GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                                   UNITCONV , 'FARNE', FARN, RC )
                  IF ( RC /= HCO_SUCCESS ) THEN 
                     CALL HCO_ERROR( 'GET_EMTERP_MEGAN_FARNE', RC )
                     ERR = .TRUE.
                     EXIT
                  ENDIF 
      
                  ! ---------------------------------------------------
                  ! b_Caryophyllene 
                  CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                                   GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                                   UNITCONV , 'BCARE', BCAR, RC )
                  IF ( RC /= HCO_SUCCESS ) THEN 
                     CALL HCO_ERROR( 'GET_EMTERP_MEGAN_BCARE', RC )
                     ERR = .TRUE.
                     EXIT
                  ENDIF 
      
                  ! ---------------------------------------------------
                  ! Other sesquiterpene
                  CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                                   GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                                   UNITCONV , 'QSQTE', OSQT, RC )
                  IF ( RC /= HCO_SUCCESS ) THEN 
                     CALL HCO_ERROR( 'GET_EMTERP_MEGAN_QSQTE', RC )
                     ERR = .TRUE.
                     EXIT
                  ENDIF 

                  FLUXSESQ(I,J) = ( FARN + BCAR + OSQT ) * FC3
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
                  FLUXMTPA(I,J) = DIURORVC * FC4 * 0.04d0 ! using campher
               ENDIF

               ! LIMONENE not defined (derived from TERP)
               IF ( IDTLIMO > 0 ) THEN
                  FLUXLIMO(I,J) = 0.0d0
               ENDIF  

               ! MTPO
               IF ( IDTMTPO > 0 ) THEN
                  FLUXMTPO(I,J) = DIURORVC * FC2 * 0.09d0 ! using LINALOOL 
               ENDIF

               ! IDTSESQ 
               IF ( IDTSESQ > 0 ) THEN
                  FLUXSESQ(I,J) = DIURORVC * FC3 * 0.05d0 
               ENDIF

            ENDIF
         ENDIF
 
         !--------------------------------------------------------------------
         !  ND46 diagnostic: Biogenic emissions 
         !
         !     AD46(:,:,1)  = Total biogenic ISOP     emissions [atoms C/cm2/s]
         !     AD46(:,:,2)  = Total biogenic ACET     emissions [atoms C/cm2/s]
         !     AD46(:,:,3)  = Total biogenic PRPE     emissions [atoms C/cm2/s]
         !     AD46(:,:,4)  = Total biogenic MONOT    emissions [atoms C/cm2/s]
         !     AD46(:,:,5)  = Total biogenic MBO      emissions [atoms C/cm2/s]
         !     AD46(:,:,6)  = Total biogenic C2H4     emissions [atoms C/cm2/s]
         !     AD46(:,:,7)  = Total biogenic a-pinene emissions [atoms C/cm2/s]
         !     AD46(:,:,8)  = Total biogenic b-pinene emissions [atoms C/cm2/s]
         !     AD46(:,:,9)  = Total biogenic limonene emissions [atoms C/cm2/s]
         !     AD46(:,:,10) = Total biogenic sabinene emissions [atoms C/cm2/s]
         !     AD46(:,:,11) = Total biogenic mycrene  emissions [atoms C/cm2/s]
         !     AD46(:,:,12) = Total biogenic 3-carene emissions [atoms C/cm2/s]
         !     AD46(:,:,13) = Total biogenic ocimene  emissions [atoms C/cm2/s]
         !     AD46(:,:,14) = Total biogenic CHBr3    emissions [kg/m2/s]
         !     AD46(:,:,15) = Total biogenic CH2Br2   emissions [kg/m2/s]
         !     AD46(:,:,16) = Total biogenic SSBr2    emissions [kg/m2/s]
         !
         !  NOTES: 
         !  (1 ) Now make ACET tracer #2 and PRPE tracer #3 (bmy, 9/13/01)
         !  (2 ) Now archive ND46 as [atoms C/cm2/s] here (bmy, 9/13/01)
         !  (3 ) Added MBO emission diagnostics [atoms C/cm2/s] (bmy, tmf, 10/20/05)
         !  22 Oct 2012 - R. Yantosca - Prevent out-of-bounds errors in ND46 if bromine
!                                       emissions are turned off (e.g. for tagged CO)
         !--------------------------------------------------------------------

!               IF ( ND46 > 0 ) THEN
!
!                  ! ISOP emissions [atoms C/cm2/s] -- tracer #1
!                  AD46(I,J,1) = AD46(I,J,1) + ( EMIS / AREA_CM2 /DTSRCE)
!
!                  ! ACET emissions [atoms C/cm2/s] -- tracer #2
!                  AD46(I,J,2) = AD46(I,J,2) + ( BIO_ACET / AREA_CM2 )
!
!                  ! PRPE emissions [atoms C/cm2/s] -- tracer #3
!                  AD46(I,J,3) = AD46(I,J,3) + 
!     &                          ( EMIS * BIOSCAL / AREA_CM2 / DTSRCE )
!
!                  ! Monoterpene emissions [atoms C/cm2/s] -- tracer #4
!                  AD46(I,J,4) = AD46(I,J,4) + ( EMMO / AREA_CM2 /DTSRCE) 
!
!                  ! MBO emissions [atoms C/cm2/s] -- tracer #5
!                  AD46(I,J,5) = AD46(I,J,5) + ( EMMB / AREA_CM2 /DTSRCE) 
!
!                  ! C2H4 emissions [atoms C/cm2/s] -- tracer #6 (tmf, 1/13/06)
!                  AD46(I,J,6) = AD46(I,J,6) + 
!     &            ( EMIS * BIOSCALEC2H4 / AREA_CM2 / DTSRCE )
!
!               ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !              MEGAN v2.1 - (mpb,2009)               !
!               ! ++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!               ! Aplha Pinene emissions [atoms C/cm2/s] -- tracer #6
!               AD46(I,J,7) = AD46(I,J,7) + ( APINE / AREA_CM2 / DTSRCE ) 
!
!               ! Beta Pinene emissions [atoms C/cm2/s] -- tracer #7
!               AD46(I,J,8) = AD46(I,J,8) + ( BPINE / AREA_CM2 / DTSRCE ) 
!
!               ! Limonene emissions [atoms C/cm2/s] -- tracer #8
!               AD46(I,J,9) = AD46(I,J,9) + ( LIMON / AREA_CM2 / DTSRCE ) 
! 
!               ! Sabinene emissions [atoms C/cm2/s] -- tracer #9
!              AD46(I,J,10) = AD46(I,J,10) + ( SABIN / AREA_CM2 / DTSRCE) 
!
!               ! Mycrene emissions [atoms C/cm2/s] -- tracer #10
!              AD46(I,J,11) = AD46(I,J,11) + ( MYRCN / AREA_CM2 / DTSRCE) 
!
!               ! 3-Carene emissions [atoms C/cm2/s] -- tracer #11
!              AD46(I,J,12) = AD46(I,J,12) + ( CAREN / AREA_CM2 / DTSRCE) 
!
!               ! Ocimene emissions [atoms C/cm2/s] -- tracer #12
!              AD46(I,J,13) = AD46(I,J,13) + ( OCIMN / AREA_CM2 / DTSRCE) 
!
!               
!               ENDIF !ND46 

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
         CALL HCO_EmisAdd( HcoState, FLUXISOP, IDTISOP, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXISOP
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr,
     &                         Cat=-1, Hier=-1, HcoID=IDTISOP,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MONOTERPENES 
      IF ( (ExtNrMono > 0 ) .AND. (IDTMONX > 0) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXEMMO, IDTMONX, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXEMMO
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrMono,
     &                         Cat=-1, Hier=-1, HcoID=IDTMONX,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ACETONE 
      IF ( IDTACET > 0 ) THEN

         ! Eventually add individual diagnostics. These are assumed to
         ! have names MEGAN_ACET_MONO, MEGAN_ACET_MBO, and 
         ! MEGAN_ACET_DIRECT
         Arr2D     => FLUXACETmo 
         DiagnName =  'MEGAN_ACET_MONO'
         CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, 
     &                      cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         Arr2D     => FLUXACETmb
         DiagnName =  'MEGAN_ACET_MBO'
         CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, 
     &                      cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         Arr2D     => FLUXACETbg 
         DiagnName =  'MEGAN_ACET_DIRECT'
         CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, 
     &                      cName=TRIM(DiagnName), Array2D=Arr2D, RC=RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Add flux to emission array
         FLUXACETbg = FLUXACETbg + FLUXACETmo + FLUXACETmb
         CALL HCO_EmisAdd( HcoState, FLUXACETbg, IDTACET, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXACETbg
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr,
     &                         Cat=-1, Hier=-1, HcoID=IDTACET,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! CO
      IF ( ( ExtNrMono > 0 ) .AND. ( IDTCO > 0 ) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXCO, IDTCO, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXCO
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrMono,
     &                         Cat=-1, Hier=-1, HcoID=IDTCO,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ALKENES 
      IF ( IDTPRPE > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXPRPE, IDTPRPE, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXPRPE
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr,
     &                         Cat=-1, Hier=-1, HcoID=IDTPRPE,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! ETHENE 
      IF ( IDTC2H4 > 0 ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXC2H4, IDTC2H4, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXC2H4
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr,
     &                         Cat=-1, Hier=-1, HcoID=IDTC2H4,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! OC 
      IF ( (ExtNrMono>0) .AND. (ExtNrSoa<=0) .AND. (IDTOCPI>0) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXOCPI, IDTOCPI, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXOCPI
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrMono,
     &                         Cat=-1, Hier=-1, HcoID=IDTOCPI,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MTPA 
      IF ( (ExtNrSoa>0) .AND. (IDTMTPA>0) ) THEN 

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXMTPA, IDTMTPA, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXMTPA
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrSoa,
     &                         Cat=-1, Hier=-1, HcoID=IDTMTPA,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! MTPO 
      IF ( (ExtNrSoa>0) .AND. (IDTMTPO>0) ) THEN 

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXMTPO, IDTMTPO, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXMTPO
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrSoa,
     &                         Cat=-1, Hier=-1, HcoID=IDTMTPO,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! LIMONENE 
      IF ( (ExtNrSoa>0) .AND. (IDTLIMO>0) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXLIMO, IDTLIMO, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXLIMO
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrSoa,
     &                         Cat=-1, Hier=-1, HcoID=IDTLIMO,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      ! ----------------------------------------------------------------
      ! SESQ 
      IF ( (ExtNrSoa>0) .AND. (IDTSESQ>0) ) THEN

         ! Add flux to emission array
         CALL HCO_EmisAdd( HcoState, FLUXSESQ, IDTSESQ, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => FLUXSESQ
            CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNrSoa,
     &                         Cat=-1, Hier=-1, HcoID=IDTSESQ,
     &                         AutoFill=1, Array2D=Arr2D, RC=RC   )
            IF ( RC /= HCO_SUCCESS ) RETURN 
            Arr2D => NULL() 
         ENDIF
      ENDIF

      !=================================================================
      ! ALL DONE! 
      !=================================================================

!      ! testing only
!      write(*,*) ''
!      write(*,*) 'MEGAN done!'
!      write(*,*) 'total LAI: ', SUM(ExtState%GC_LAI%Arr%Val)
!      write(*,*) 'total ISOP: ', SUM(FLUXISOP)
!      write(*,*) 'total ACET: ', SUM(FLUXACETmo)+SUM(FLUXACETmb)
!     &                          +SUM(FLUXACETbg)
!      write(*,*) 'total CO  : ', SUM(FLUXCO  )
!      write(*,*) 'total PRPE: ', SUM(FLUXPRPE)
!      write(*,*) 'total OCPI: ', SUM(FLUXOCPI)
!      write(*,*) ''

      ! Keep LAI in memory (ckeller, 10/9/2014)
      LAI_PREVSTEP = ExtState%GC_LAI%Arr%Val

      ! Return w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE HCOX_Megan_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_EmIsop_Megan
!
! !DESCRIPTION: Subroutine Get\_EmIsop\_Megan computes isoprene emissions in 
!  units of [atoms C/box] using the MEGAN inventory.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_EMISOP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM,
     &                             XNUMOL, EMISOP, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      INTEGER, INTENT(IN) :: I, J     ! GEOS-Chem lon & lat indices
      REAL*8,  INTENT(IN) :: GAMMA_LAI ! Leaf area correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_P   ! Light correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_SM  ! Soil moisture correction factor 
      REAL*8,  INTENT(IN) :: XNUMOL   ! Number of atoms C / kg C 
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: EMISOP   ! Isoprene emissions [atoms C/box]
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(INOUT) :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  :: TS, SUNCOS
      REAL*8  :: GAMMA_LEAF_AGE, GAMMA_T, GAMMA_LAI_ISOP
      REAL*8  :: D_BTW_M
      REAL*8  :: ISOLAI, PMISOLAI, MISOLAI

      !=================================================================
      ! GET_EMISOP_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMISOP          = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_LAI_ISOP  = 1.d0

      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Pass met variables. Now use only GC_LAI (ckeller, 10/9/2014)
      ISOLAI   = ExtState%GC_LAI%Arr%Val(I,J)
      PMISOLAI = LAI_PREVDAY(I,J)
      MISOLAI  = ISOLAI
!      PMISOLAI = ExtState%GC_LAI_PM%Arr%Val(I,J)
!      MISOLAI  = ExtState%GC_LAI_CM%Arr%Val(I,J)
      TS       = ExtState%TSURFK%Arr%Val(I,J)
      SUNCOS   = ExtState%SUNCOSmid%Arr%Val(I,J)


      !---------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !---------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         IF ( ISOLAI * AEF_ISOP(I,J) > 0d0 ) THEN

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI, 
     &                                 PMISOLAI, D_BTW_M, 
     &                                'ISOP' , HCOT_15_AVG(I,J) )


            ! Activity factor for temperature
            GAMMA_T  = GET_GAMMA_T_ISOP( TS, HCOT_15_AVG(I,J), 
     &                                       HCOT_DAILY(I,J)  )

            ! Activity factor for LAI. Use default (1.0) if PECCA model
            ! is not used!   
            IF ( LPECCA ) THEN
               GAMMA_LAI_ISOP = GAMMA_LAI
            ENDIF

            ! Isoprene emission is the product of all these
            EMISOP    = AEF_ISOP(I,J) * GAMMA_LAI_ISOP 
     &                * GAMMA_LEAF_AGE * GAMMA_T   * GAMMA_P * GAMMA_SM

            ! Convert from [kg/box] to [atoms C/box]
            EMISOP    = EMISOP * XNUMOL

         ENDIF
      ENDIF
   
!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' ' 
!         write(*,*) 'HEMCO MEGAN @ ', i,j
!         write(*,*) 'EMISOP: ', EMISOP
!         write(*,*) 'AEF_ISOP [kgC/m2/s]: ', AEF_ISOP(I,J,1,1)
!         write(*,*) 'GAMMA_LAI: ', GAMMA_LAI 
!         write(*,*) 'GAMMA_LEAF_AGE: ', GAMMA_LEAF_AGE
!         write(*,*) 'GAMMA_T: ', GAMMA_T
!         write(*,*) 'GAMMA_P: ', GAMMA_P 
!         write(*,*) 'GAMMA_SM: ', GAMMA_SM 
!         write(*,*) 'TS: ', TS 
!         write(*,*) 'HCOT_15_AVG: ', HCOT_15_AVG(I,J) 
!         write(*,*) 'HCOT_DAILY: ', HCOT_DAILY(I,J) 
!         write(*,*) 'HCOPARDR_15_AVG: ', HCOPARDR_15_AVG(I,J) 
!         write(*,*) 'HCOPARDF_15_AVG: ', HCOPARDF_15_AVG(I,J) 
!         write(*,*) 'ISOLAI  : ', ISOLAI
!         write(*,*) 'MISOLAI : ', MISOLAI
!         write(*,*) 'PMISOLAI: ', PMISOLAI
!         write(*,*) 'D_BTW_M : ', D_BTW_M 
!         write(*,*) ' ' 
!      endif

      ! Leave w/ success
      RC = HCO_SUCCESS
 
      END SUBROUTINE GET_EMISOP_MEGAN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_EmMbo_Megan
!
! !DESCRIPTION: Subroutine Get\_EmMbo\_Megan computes methylbutenol emissions
!  in units of [atoms C/box] using the MEGAN inventory.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_EMMBO_MEGAN( HcoState, ExtState, I, J, 
     &                            GAMMA_LAI, GAMMA_P, GAMMA_SM,
     &                            XNUMOL, EMMBO, RC ) 
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      INTEGER, INTENT(IN) :: I, J      ! GEOS-Chem lon & lat indices
      REAL*8,  INTENT(IN) :: GAMMA_LAI ! Leaf area correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_P   ! Light correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_SM  ! Soil moisture correction factor 
      REAL*8,  INTENT(IN) :: XNUMOL    ! Number of atoms C / kg C 
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT) :: EMMBO   ! Methylbutenol emissions [atoms C/box]
!
! !INPUT/OUTPUT PARAMETERS: 
!

      INTEGER, INTENT(INOUT) :: RC
! 
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2004, 2006
!  (2 ) Wang,    et al, 1998
!  (3 ) Guenther et al, 2007, MEGAN v2.1 User mannual 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8            :: GAMMA_LEAF_AGE, GAMMA_T, TS, SUNCOS
      REAL*8            :: GAMMA_LAI_EMMB, D_BTW_M
      REAL*8            :: ISOLAI, PMISOLAI, MISOLAI
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: BETA = 0.09

      !=================================================================
      ! GET_EMMBO_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMMBO           = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      GAMMA_LAI_EMMB  = 1.d0

      ! Pass met variables. Now use only GC_LAI (ckeller, 10/9/2014)
      ISOLAI   = ExtState%GC_LAI%Arr%Val(I,J)
      PMISOLAI = LAI_PREVDAY(I,J)
      MISOLAI = ISOLAI
!      PMISOLAI = ExtState%GC_LAI_PM%Arr%Val(I,J)
!      MISOLAI  = ExtState%GC_LAI_CM%Arr%Val(I,J)
      TS       = ExtState%TSURFK%Arr%Val(I,J)
      SUNCOS   = ExtState%SUNCOSmid%Arr%Val(I,J)


      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      !-------------------------------------------------
      ! Do only during day and over continent
      ! Only interested in terrestrial biosphere (pip)
      !-------------------------------------------------
      IF ( SUNCOS > 0d0 ) THEN

         IF ( ISOLAI * AEF_MBO(I,J) > 0d0 ) THEN

            ! Activity factor for leaf age 
            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI, 
     &                                 PMISOLAI, D_BTW_M, 
     &                                'MBOT' , HCOT_15_AVG(I,J) )
         
            ! Activity factor for temperature
            GAMMA_T = GET_GAMMA_T_NISOP( TS, BETA )

            ! Activity factor for LAI. Use default (1.0) if PECCA model
            ! is not used!   
            IF ( LPECCA ) THEN
               GAMMA_LAI_EMMB = GAMMA_LAI
            ENDIF

            ! MBO emissions in [kgC/m2/s]
            EMMBO = AEF_MBO(I,J)     * GAMMA_LAI_EMMB 
     &            * GAMMA_LEAF_AGE * GAMMA_T  * GAMMA_P   * GAMMA_SM

            ! Convert if necessary 
            EMMBO = EMMBO * XNUMOL

         ENDIF
      ENDIF

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) 'HEMCO GET_EMMBO_MEGAN @',i,j
!         write(*,*) 'GAMMA_LEAF_AGE: ', GAMMA_LEAF_AGE
!         write(*,*) 'GAMMA_T: ', GAMMA_T
!         write(*,*) 'GAMMA_SM: ', GAMMA_SM
!         write(*,*) 'GAMMA_P: ', GAMMA_P
!         write(*,*) 'GAMMA_LAI: ', GAMMA_LAI_EMMB
!      endif
!
      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_EMMBO_MEGAN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_EmTerp_Megan
!
! !DESCRIPTION: Function Get\_EmTerp\_Megan computes monoterpene and
!  sesquiterpene emissions for individual species in units of [atoms C/box]
!  using the MEGAN v2.1 inventory.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_EMTERP_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM, 
     &                             XNUMOL, TERP_SPECIES, EMTERP, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      INTEGER,          INTENT(IN) :: I, J          ! Lon & lat indices
      REAL*8,  INTENT(IN) :: GAMMA_LAI ! Leaf area correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_P   ! Light correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_SM  ! Soil moisture correction factor 
      REAL*8,           INTENT(IN) :: XNUMOL        ! Number of atoms C / kg C 
      CHARACTER(LEN=5), INTENT(IN) :: TERP_SPECIES  ! Monoterpene species name
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,          INTENT(OUT) :: EMTERP        ! Emissions [kg/m2/s]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER,       INTENT(INOUT) :: RC            ! Return code 
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2004, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 User Manual
!  (3 ) Sakulyanontvittaya et al, 2008
! 
! !REVISION HISTORY: 
!  (1 ) Written by Michael Barkley (2008), based on old monoterpene code by 
!       dsa,tmf.
!  (2 ) Uses gamma factors instead of exchange factors, this includes
!        calling of a new temperature algorithm which use a beta factor.
!        (mpb,2008)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  11 Apr 2012 - R. Yantosca - Now use data from modis_lai_mod.F90
!  13 Aug 2013 - M. Sulprizio- Renamed from GET_EMMONOG_MEGAN and added
!                              sesquiterpenes (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8            :: BETA
      REAL*8            :: GAMMA_LEAF_AGE, GAMMA_T, TS
      REAL*8            :: D_BTW_M
      REAL*8            :: LDF 
      REAL*8            :: AEF_SPARE
      REAL*8            :: ISOLAI, PMISOLAI, MISOLAI
      CHARACTER(LEN=4)  :: AGE_SPECIES ! Identify age param & pass to GAMMA_AGE
      CHARACTER(LEN=255):: MSG

      !=================================================================
      ! GET_EMMONOT_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMTERP          = 0.d0
      GAMMA_LEAF_AGE  = 0.d0
      AEF_SPARE       = 0.d0

      ! Pass met variables. Now use only GC_LAI (ckeller, 10/9/2014)
      ISOLAI   = ExtState%GC_LAI%Arr%Val(I,J)
      PMISOLAI = LAI_PREVDAY(I,J)
      MISOLAI = ISOLAI
!      PMISOLAI = ExtState%GC_LAI_PM%Arr%Val(I,J)
!      MISOLAI  = ExtState%GC_LAI_CM%Arr%Val(I,J)
      TS       = ExtState%TSURFK%Arr%Val(I,J)


      ! Number of days between MISOLAI and PMISOLAI
      D_BTW_M  = DBLE( DAYS_BTW_M )

      ! Need to determine which monoterpene AEFs 
      ! we need to use (mpb,2009)
      ! SOAudpate: add BETA, AGE_SPECIES (hotp 7/29/10)
      SELECT CASE( TERP_SPECIES )
      CASE( 'APINE' )
         AEF_SPARE      = AEF_APINE(I,J)
         LDF            = 0.1
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'BPINE' )
         AEF_SPARE      = AEF_BPINE(I,J)
         LDF            = 0.1
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'LIMON' )
         AEF_SPARE      = AEF_LIMON(I,J)
         LDF            = 0.05
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'SABIN' )
         AEF_SPARE      = AEF_SABIN(I,J)
         LDF            = 0.1
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'MYRCN' )
         AEF_SPARE      = AEF_MYRCN(I,J)
         LDF            = 0.05
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'CAREN' )
         AEF_SPARE      = AEF_CAREN(I,J)
         LDF            = 0.05
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      CASE( 'OCIMN' )
         AEF_SPARE      = AEF_OCIMN(I,J)
         LDF            = 0.8
         BETA           = 0.09
         AGE_SPECIES    = 'MONO'
      ! SOAupdate: add Sesquiterpenes (hotp 3/1/10)
      ! Sesquiterpene LDF based on Sakulyanontvittaya 2008
      ! beta from MEGANv2.04 (megan CDP)
      ! make name 5 characters
      CASE( 'FARNE' )
         AEF_SPARE      = AEF_GEN(I,J,IDMGFARN)
         LDF            = 0.5d0
         BETA           = 0.17d0
         AGE_SPECIES    = 'SESQ'
      CASE( 'BCARE' )
         AEF_SPARE      = AEF_GEN(I,J,IDMGBCAR)
         LDF            = 0.5d0
         BETA           = 0.17d0
         AGE_SPECIES    = 'SESQ'
      CASE( 'OSQTE' )
         AEF_SPARE      = AEF_GEN(I,J,IDMGOSQT)
         LDF            = 0.5d0
         BETA           = 0.17d0
         AGE_SPECIES    = 'SESQ'
      ! Other monoterpenes (hotp 3/5/10)
      CASE( 'OMTPE' )
         AEF_SPARE      = AEF_GEN(I,J,IDMGOMTP)
         LDF            = 0.1d0
         BETA           = 0.1d0
         AGE_SPECIES    = 'MONO'
      CASE DEFAULT
         MSG = 'Invalid TERPENE species'
         CALL HCO_ERROR( MSG, RC, THISLOC='GET_EMTERP_MEGAN' )
         RETURN 
      END SELECT   

      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      ! AEF_SPARE no longer I,J (hotp 3/10/10)
      IF ( ISOLAI * AEF_SPARE > 0d0 ) THEN

         ! Activity factor for leaf age 
         GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI, 
     &                              PMISOLAI, D_BTW_M, 
     &                              AGE_SPECIES, HCOT_15_AVG(I,J) )

         ! Activity factor for temperature
         GAMMA_T = GET_GAMMA_T_NISOP( TS, BETA )
    
         ! Monoterpene emission is the product of all these; must be 
         ! careful to distinguish between canopy & PECCA models.
         IF ( LPECCA ) THEN
            EMTERP = AEF_SPARE * GAMMA_LEAF_AGE * GAMMA_T 
     &                         * GAMMA_SM       * GAMMA_LAI
     &                         * ( (1.d0 - LDF) + (LDF * GAMMA_P) )

         ELSE 

            ! removed I,J from AEF_SPARE (hotp 3/10/10)
            EMTERP = AEF_SPARE * GAMMA_LEAF_AGE * GAMMA_T 
     &                         * GAMMA_SM       
     &             * ( GAMMA_LAI * (1.d0 - LDF) + (LDF * GAMMA_P) )

         END IF 

         ! Convert from [kg/box] to [atoms C/box]
         EMTERP  = EMTERP * XNUMOL

      END IF

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_EMTERP_MEGAN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_EmaAet_megan
!
! !DESCRIPTION: Subroutine Get\_EmAcet\_Megan computes acetone emissions in
!  units of [atomsC/box] using the MEGAN inventory.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_EMACET_MEGAN( HcoState, ExtState, I, J, 
     &                             GAMMA_LAI, GAMMA_P, GAMMA_SM,
     &                             XNUMOL, EMACET, RC ) 
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      INTEGER, INTENT(IN) :: I,  J   ! GEOS-Chem lon & lat indices
      REAL*8,  INTENT(IN) :: GAMMA_LAI ! Leaf area correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_P   ! Light correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_SM  ! Soil moisture correction factor 
      REAL*8,  INTENT(IN) :: XNUMOL  ! Number of atoms C / kg C 
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT) :: EMACET  ! Acetone emissions [atoms C/box]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT) :: RC
! 
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2004, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 User Manual
!
! !REVISION HISTORY:
!  (1 ) Written by Michael Barkley (2008), based on old monoterpene code by 
!       dsa,tmf.
!  (2 ) Uses gamma factors instead of exchange factors, this includes
!        calling of a new temperature algorithm which use a beta factor.
!        (mpb,2008) 
!  24 May 2011 - E. Fischer  - Modified for acetone. Function GET_EMACET_MEGAN
!                              is called from "acetone_mod.f"
!  06 Dec 2011 - M. Payer    - Added ProTeX headers
!  27 Mar 2012 - R. Yantosca - Avoid segfault errors if LMEGAN=.FALSE.
!  11 Apr 2012 - R. Yantosca - Now use data from modis_lai_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Fraction of emission that is light dependent (A. Guenther 5/25/2011)
      REAL*8, PARAMETER   :: LDF    = 0.20d0

      ! Temperature response factor for light-independent emissions
      ! (A. Guenther 5/25/2011)
      REAL*8, PARAMETER   :: BETA   = 0.10d0
!
! !LOCAL VARIABLES:
!
      REAL*8              :: GAMMA_LEAF_AGE
      REAL*8              :: GAMMA_T
      REAL*8              :: D_BTW_M, TS, SUNCOS
      REAL*8              :: ISOLAI, PMISOLAI, MISOLAI

      !=================================================================
      ! GET_EMACET_MEGAN begins here!
      !================================================================= 

      ! Initialize return value & activity factors
      EMACET          = 0.d0
      GAMMA_T         = 0.d0
      GAMMA_LEAF_AGE  = 1.d0

      ! Pass met variables. Now use only GC_LAI (ckeller, 10/9/2014)
      ISOLAI   = ExtState%GC_LAI%Arr%Val(I,J)
      PMISOLAI = LAI_PREVDAY(I,J)
      MISOLAI = ISOLAI
!      PMISOLAI = ExtState%GC_LAI_PM%Arr%Val(I,J)
!      MISOLAI  = ExtState%GC_LAI_CM%Arr%Val(I,J)
      TS       = ExtState%TSURFK%Arr%Val(I,J)
      SUNCOS   = ExtState%SUNCOSmid%Arr%Val(I,J)


      !-----------------------------------------------------
      ! Only interested in terrestrial biosphere (pip)
      ! If (local LAI != 0 .AND. baseline emission !=0 ) 
      !-----------------------------------------------------
      IF ( ISOLAI * AEF_ACET(I,J) > 0d0 ) THEN

         ! Calculate gamma PAR only if sunlight conditions
         IF ( SUNCOS > 0d0 ) THEN

            ! --------------------------------------------------
            ! GAMMA_age
            ! --------------------------------------------------
            ! Activity factor for leaf age
            ! there is no leaf age dependence for acetone so 
            ! GAMMA_LEAF_AGE should be 1
!            GAMMA_LEAF_AGE = GET_GAMMA_LEAF_AGE( MISOLAI(I,J), 
!     &                                 PMISOLAI(I,J), D_BTW_M, 
!     &                                 'MOHT', HCOT_15_AVG(I,J) )

         ENDIF
    
         ! --------------------------------------------------
         ! GAMMA_T
         ! --------------------------------------------------
         ! Activity factor for temperature
         GAMMA_T = GET_GAMMA_T_NISOP( TS, BETA )

         ! Acetone emission is the product of all these; must be 
         ! careful to distinguish between canopy & PECCA models.
         IF ( LPECCA ) THEN

            EMACET   = AEF_ACET(I,J) * GAMMA_LEAF_AGE 
     &               * GAMMA_T   * GAMMA_SM   * GAMMA_LAI
     &               * ( (1.d0 - LDF) + (LDF * GAMMA_P) )

         ELSE 

            EMACET   = AEF_ACET(I,J)     * GAMMA_LEAF_AGE 
     &               * GAMMA_T                    * GAMMA_SM       
     &               * (  GAMMA_LAI * (1.d0 - LDF) + (LDF * GAMMA_P) )

         END IF 

         ! Convert from [kg/box] to [atoms C/box]
         EMACET  = EMACET * XNUMOL
 
      END IF

      ! Return w/ success
      RC = HCO_SUCCESS

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) 'HEMCO GET_EMACET_MEGAN @',i,j
!         write(*,*) 'GAMMA_LEAF_AGE: ', GAMMA_LEAF_AGE
!         write(*,*) 'GAMMA_T: ', GAMMA_T
!         write(*,*) 'GAMMA_SM: ', GAMMA_SM
!         write(*,*) 'GAMMA_P: ', GAMMA_P
!         write(*,*) 'GAMMA_LAI: ', GAMMA_LAI
!      endif
!
      END SUBROUTINE GET_EMACET_MEGAN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_EmMonot_Megan
!
! !DESCRIPTION: Subroutine Get\_EmMonot\_Megan computes the total 
!  monoterpene emissions in units of [atoms C/box] using the MEGAN v2.1 
!  inventory.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_EMMONOT_MEGAN( HcoState, ExtState, I, J, 
     &                              GAMMA_LAI, GAMMA_P, GAMMA_SM,
     &                              XNUMOL, EMMONOT, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_STATE), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      INTEGER, INTENT(IN) :: I, J      ! Lon & lat indices
      REAL*8,  INTENT(IN) :: GAMMA_LAI ! Leaf area correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_P   ! Light correction factor 
      REAL*8,  INTENT(IN) :: GAMMA_SM  ! Soil moisture correction factor 
      REAL*8,  INTENT(IN) :: XNUMOL    ! Number of atoms C / kg C 
!
! !OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(OUT) :: EMMONOT   ! Monoterpene emissions [atoms C/box]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(INOUT) :: RC
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995, 1999, 2000, 2006
!  (2 ) Guenther et al, 2007, MEGAN v2.1 User Manual
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Michael Barkley (mpb,2009).
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  09 Mar 2010 - H.O.T. Pye  - Change order of arguments in call to 
!                              routine GET_EMMONOG_MEGAN
!  11 Apr 2012 - R. Yantosca - Now use data from modis_lai_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER :: N = 7
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL*8             :: MONO
      INTEGER            :: K 

      ! Arrays
      CHARACTER(LEN=5)   :: SPECIES(N) = (/ 'APINE', 'BPINE', 'LIMON', 
     &                                      'SABIN', 'MYRCN', 'CAREN', 
     &                                      'OCIMN'                   /)

      !=================================================================
      ! GET_EMMONOT_MEGAN begins here!
      !================================================================= 

      ! Initialize
      EMMONOT = 0.d0

      ! Loop over all species
      DO K = 1, N

         ! Get contribution from each individual species [kg C/box]
         CALL GET_EMTERP_MEGAN( HcoState, ExtState, I, J,  
     &                          GAMMA_LAI, GAMMA_P, GAMMA_SM,
     &                          XNUMOL, SPECIES(K), MONO, RC )

         ! Error trap
         IF ( RC /= HCO_SUCCESS ) THEN
            RC = HCO_FAIL 
            RETURN
         ENDIF

         ! Compute overall total
         EMMONOT = EMMONOT + MONO
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_EMMONOT_MEGAN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_P_Pecca
!
! !DESCRIPTION: Computes the PECCA gamma activity factor with sensitivity 
!  to LIGHT.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_P_PECCA( HcoState, ExtState, I, J, Q_DIR_2, 
     &                            Q_DIFF_2, HCOPARDR_AVG_SIM, 
     &                            HCOPARDF_AVG_SIM )
     &              RESULT( GAMMA_P_PECCA )
!
! !USES:
!
      USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal 
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER  :: HcoState
      TYPE(Ext_State), POINTER  :: ExtState
      INTEGER,       INTENT(IN) :: I,  J             ! Lon & lat indices 
      REAL*8,        INTENT(IN) :: HCOPARDR_AVG_SIM     ! Average direct PAR [W/m2]
      REAL*8,        INTENT(IN) :: HCOPARDF_AVG_SIM     ! Average diffuse PAR [W/m2]
      REAL*8,        INTENT(IN) :: Q_DIR_2           ! Direct PAR [umol/m2/s]
      REAL*8,        INTENT(IN) :: Q_DIFF_2          ! Diffuse PAR [umol/m2/s]
!
! !RETURN VALUE:
!
      REAL*8                    :: GAMMA_P_PECCA     ! GAMMA factor for light
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
      REAL*8  :: LAT 
      REAL*8  :: mmHCOPARDR_DAILY
      REAL*8  :: mmHCOPARDF_DAILY
      REAL*8  :: PAC_DAILY, PAC_INSTANT, C_PPFD
      REAL*8  :: PTOA, PHI
      REAL*8  :: BETA,   SINbeta 
      INTEGER :: DOY, RC 
      REAL*8  :: AAA, BBB
      REAL*4  :: LocalHour

      !-----------------------------------------------------------------
      ! Compute GAMMA_P_PECCA
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD   = 0.0d0
      PTOA     = 0.0d0

      ! Convert past light conditions to micromol/m2/s 
      mmHCOPARDR_DAILY   = HCOPARDR_AVG_SIM  * WM2_TO_UMOLM2S
      mmHCOPARDF_DAILY   = HCOPARDF_AVG_SIM  * WM2_TO_UMOLM2S

      ! Work out the light at the top of the canopy.
      PAC_DAILY    = mmHCOPARDR_DAILY + mmHCOPARDF_DAILY
      PAC_INSTANT  = Q_DIR_2       +  Q_DIFF_2

      ! Get day of year, local-time and latitude
      ! TODO: Evaluate RC?
      CALL HcoClock_Get ( cDOY = DOY, RC=RC )
      CALL HcoClock_GetLocal ( I, cH = LocalHour, RC=RC )
      LAT = HcoState%Grid%YMID%Val( I, J )

      ! Get solar elevation angle
      SINbeta      =  SOLAR_ANGLE( HcoState, DOY , LocalHour, LAT )
      BETA         =  ASIN( SINbeta ) * RAD2D       

      IF ( SINbeta .LE. 0.0d0 ) THEN

         GAMMA_P_PECCA = 0.0d0

      ELSEIF ( SINbeta .GT. 0.0d0 ) THEN       

         ! PPFD at top of atmosphere
         PTOA    = 3000.0d0 + 99.0d0 * 
     &             COS( 2.d0 * 3.14d0 *( DOY - 10 ) / 365 )

         ! Above canopy transmission
         PHI     = PAC_INSTANT / ( SINbeta * PTOA )

         ! Work out gamma P
         BBB     = 1.0d0 + 0.0005d0 *( PAC_DAILY - 400.0d0  ) 
         AAA     = ( 2.46d0 * BBB * PHI ) - ( 0.9d0 * PHI**2 )

         GAMMA_P_PECCA = SINbeta * AAA

      ENDIF

       ! Screen unforced errors. IF solar elevation angle is 
       ! less than 1 THEN gamma_p can not be greater than 0.1.
       IF ( BETA .LT. 1.0 .AND. GAMMA_P_PECCA .GT. 0.1) THEN
          GAMMA_P_PECCA  = 0.0
       ENDIF
   
      ! Prevent negative values
      GAMMA_P_PECCA = MAX( GAMMA_P_PECCA , 0d0 )

!      ! testing only
!      if ( i==ix .and. j==iy ) then
!         write(*,*) ' '
!         write(*,*) 'HEMCO GAMMA_P_PECCA: ', GAMMA_P_PECCA
!         write(*,*) 'DOY: ', DOY
!         write(*,*) 'LUT: ', LUT
!         write(*,*) 'LAT: ', LAT
!         write(*,*) 'SINbeta: ', SINbeta
!         write(*,*) 'BETA: ', BETA 
!         write(*,*) 'PTOA: ', PTOA 
!         write(*,*) 'PHI : ', PHI 
!         write(*,*) 'BBB: ', BBB 
!         write(*,*) 'AAA: ', AAA
!         write(*,*) 'PAC_DAILY: ', PAC_DAILY
!         write(*,*) 'PAC_INSTANT: ', PAC_INSTANT
!      endif

      END FUNCTION GET_GAMMA_P_PECCA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      TYPE(HCO_State), POINTER    :: HcoState
      INTEGER,         INTENT(IN) :: DOY       ! Day of year
      REAL*4,          INTENT(IN) :: SHOUR     ! Local time 
      REAL*8,          INTENT(IN) :: LAT       ! Latitude
!
! !RETURN VALUE:
!
      REAL*8                      :: SINbeta   ! Sin of the local solar angle
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
      REAL*8 :: BETA                        ! solar elevation angle
      REAL*8 :: sindelta, cosdelta, A, B

      ! Calculation of sin beta 
      sindelta = -SIN( 0.40907d0 ) * 
     &            COS( 6.28d0 * ( DOY + 10 ) / 365 )

      cosdelta = (1-sindelta**2.)**0.5

      A = SIN( LAT * D2RAD ) * sindelta
      B = COS( LAT * D2RAD ) * cosdelta

      SINbeta = A + B * 
     &          COS( 2.0d0 * HcoState%Phys%PI * ( SHOUR-12 ) / 24 ) 

      END FUNCTION SOLAR_ANGLE 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_Isop
!
! !DESCRIPTION: Function Get\_Gamma\_T\_Isop computes the temperature 
!  sensitivity for ISOPRENE ONLY.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_ISOP( T, PHCOT_15, PT_1 ) RESULT( GAMMA_T )
!
! !INPUT PARAMETERS: 
!

      ! Current leaf temperature, the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL*8,  INTENT(IN) :: T 

 
      ! Average leaf temperature over the past 15 days
      REAL*8,  INTENT(IN) :: PHCOT_15

      ! Average leaf temperature over the past arbitray day(s).
      ! This is not used at present (but might be soon!).
      REAL*8,  INTENT(IN) :: PT_1
!
! !RETURN VALUE:
!
      ! GAMMA factor for temperature (isoprene only)
      REAL*8              :: GAMMA_T
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Guenther et al, 2006
!  (3 ) Guenther et al, MEGAN v2.1 user mannual 2007-08
! 
! 
! !REVISION HISTORY: 
!  (1 ) Includes the latest MEGAN v2.1 temperature algorithm (mpb, 2009).
!       Note, this temp-dependence is the same for the PECCA & hybrid models.
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: C_T,   CT1,   CT2
      REAL*8              :: E_OPT, T_OPT, X

      ! Ideal gas constant [J/mol/K]
      REAL*8, PARAMETER   :: R   = 8.314d-3


      ! Normalization Factor 
      REAL*8, PARAMETER   ::  NORM = 0.99558166894501043d0    
                                   
      !=================================================================!
      !             ALWAYS CHECK THE ABOVE NORMALIZATION                !
      !=================================================================!
       E_OPT = 1.75d0 * EXP( 0.08d0 * ( PHCOT_15  - 2.97d2 ) )     
       T_OPT = 3.13d2 + ( 6.0d-1 * ( PHCOT_15 - 2.97d2 ) )
       CT1   = 80d0
       CT2   = 200d0

      ! Variable related to temperature 
      X     = ( 1.d0/T_OPT - 1.d0/T ) / R

      ! C_T: Effect of temperature on leaf BVOC emission, including 
      ! effect of average temperature over previous 15 days, based on 
      ! Eq 5a, 5b, 5c from Guenther et al, 1999.
      C_T   = E_OPT * CT2 * EXP( CT1 * X ) / 
     &        ( CT2 - CT1 * ( 1.d0 - EXP( CT2 * X ) ) )

      ! Hourly emission activity = C_T
      ! Prevent negative values
      GAMMA_T = MAX( C_T * NORM , 0d0 )

      END FUNCTION GET_GAMMA_T_ISOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_T_nIsop
!
! !DESCRIPTION: Function Get\_Gamma\_T\_nIsop computes the temperature 
!  activity factor (GAMMA\_T) for BVOCs OTHER than isoprene.  Called from 
!  routines Get\_EmMonog\_Megan and Get\_Emmbo\_Megan.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_T_NISOP( T, BETA ) RESULT( GAMMA_T )
!
! !INPUT PARAMETERS: 
!
      ! Current leaf temperature [K], the surface air temperature field (TS) 
      ! is assumed equivalent to the leaf temperature over forests.
      REAL*8, INTENT(IN) :: T        

      ! Temperature factor per species (from MEGAN user manual). 
      ! Beta = 0.09 for MBO and for monoterpene species (APINE, BPINE, LIMON, 
      ! SABIN, MYRCN, CAREN, OCIMN). Pass as an argument in case this changes.
      REAL*8, INTENT(IN) :: BETA     
!
! !RETURN VALUE:
!
      REAL*8             :: GAMMA_T  ! 
!
! !REMARKS:
!  GAMMA_T =  exp[BETA*(T-Ts)]
!                                                                             .
!             where BETA   = temperature dependent parameter
!                   Ts     = standard temperature (normally 303K, 30C)
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08

! !REVISION HISTORY: 
!  (1 ) Original code by Michael Barkley (2009).
!       Note: If T = Ts  (i.e. standard conditions) then GAMMA_T = 1 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Standard reference temperature [K]
      REAL*8, PARAMETER   :: Ts = 303.0

      !=================================================================
      ! GET_GAMMAT_NISOP begins here!
      !================================================================= 

      GAMMA_T = EXP( BETA * ( T - Ts ) )

      END FUNCTION GET_GAMMA_T_NISOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_P
!
! !DESCRIPTION: Function Get\_Gamma\_P computes the gamma activity factor with
!  sensitivity to LIGHT (aka 'PAR').  Called by the functions !
!  GEt\_EmIsop\_Megan, Get\_EmMbo\_Megan, and Get\_EmMonog\_Megan.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_P( LAI, SUNCOS1, Q_DIR_2, Q_DIFF_2 ) 
     &              RESULT( GAMMA_P )
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN) :: LAI        ! Cumulative leaf area index 
      REAL*8,  INTENT(IN) :: SUNCOS1    ! Cosine of solar zenith angle 
      REAL*8,  INTENT(IN) :: Q_DIR_2    ! Direct PAR above canopy [umol/m2/s]
      REAL*8,  INTENT(IN) :: Q_DIFF_2   ! Diffuse PAR above canopy [umol/m2/s]
!
! !RETURN VALUE:
!
      REAL*8              :: GAMMA_P    ! Gamma activity factor w/r/t light
!
! !REMARKS:
!  *** REVAMPED FUNCTION ***
!                                                                             .
!  C_PPFD: Effect of increasing PPFD up to a saturation point, where emission 
!          level off, based on Eq 4abc from Guenther et al. (1999)
!          In addition, a 5 layered canopy model based on Eqs 12-16 
!          from Guenther et al. (1995) is included to correct for light 
!          attenuation in the canopy.
!                                                                             .
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 1995
!  (2 ) Wang     et al, 1998
!  (3 ) Guenther et al, 1999
!  (5 ) Guenther et al, 2004
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot and by May Fu.
!  (2 ) This code was extracted from the previous GET_HEA_TL function. 
!       (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !-----------------------------------------------------------------
      ! Canopy model variables (Eqs 12-16 from Guenther et al, 1995)
      !-----------------------------------------------------------------


      ! C_PPFD: Effect of increasing PPFD up to a saturation point, where 
      ! emissions level off, based on Eq 4abc from Guenther et al. (1999)
      ! In addition, a 5 layered canopy model based on Eqs 12-16 
      ! from Guenther et al. (1995) is included to correct for light 
      ! attenuation in the canopy.
      REAL*8              :: C_PPFD

      ! LAYERS: number of layers in canopy model
      INTEGER, PARAMETER  :: LAYERS = 5

      ! WEIGHTGAUSS : weights for gaussian integration
      REAL*8,  PARAMETER  :: WEIGHTGAUSS(LAYERS) = (/ 1.1846d-1, 
     &                                                2.3931d-1, 
     &                                                2.8444d-1, 
     &                                                2.3931d-1, 
     &                                                1.1846d-1 /)

      ! DISTGAUSS: points to evaluate fcn for gaussian integration
      REAL*8,  PARAMETER  :: DISTGAUSS(LAYERS)   = (/ 4.6910d-2, 
     &                                                2.3075d-1, 
     &                                                5.0000d-1, 
     &                                                7.6924d-1, 
     &                                                9.5308d-1 /)

      ! SCAT: Scattering coefficient
      REAL*8, PARAMETER   :: SCAT = 2.0d-1

      ! REFLD: Reflection coefficient for diffuse light
      REAL*8, PARAMETER   :: REFLD = 5.7d-2

      ! CLUSTER: clustering coefficient (accounts for leaf clumping 
      ! influence on mean projected leaf area in the direction 
      ! of the sun's beam) use 0.85 for default
      REAL*8, PARAMETER   :: CLUSTER = 8.5d-1

      ! NORMAL_FACTOR : C_PPFD calculated with LAI = 5, and above canopy 
      ! total PPFD = 1500 umol/m2/s, and Q_DIFF = 0.2 * above canopy total 
      ! PPFD.  May Fu calculated this to be 1.8967d0
      ! Quickly checked by mpb in 2008 & found to be:
      REAL*8, PARAMETER   :: NORMAL_FACTOR = 1.80437817285271950d0

      ! F_SUN: Fraction leaves sunlit
      REAL*8              :: F_SUN

      ! F_SHADE: Fraction leaves shaded
      REAL*8              :: F_SHADE

      ! LAI_DEPTH: Cumulative LAI above current layer
      REAL*8              :: LAI_DEPTH

      ! Q_SUN: Flux density of PAR on sunny leaves [umol/m2/s]
      REAL*8              :: Q_SUN
 
      ! Q_SHADE: Flux density of PAR on shaded leaves [umol/m2/s]
      REAL*8              :: Q_SHADE, Q_SHADE_1, Q_SHADE_2

      ! KB: Extinction coefficient for black leaves for direct (beam)
      REAL*8              :: KB, KBP

      ! KD: Extinction coefficient for black leaves for diffuse light
      REAL*8              :: KD,   KDP
      REAL*8              :: REFLB

      ! C_P_SUN: C_PPFD at layer I for sunlit leaves
      REAL*8              :: C_P_SUN
      
      ! C_P_SHADE: C_PPFD at layer I for shaded leaves
      REAL*8              :: C_P_SHADE

      ! C_PPFD_I    : C_PPFD at layer I
      REAL*8              :: C_PPFD_I

      ! Empirical functions (Eq 4a, 4b, 4c from Guenther et al, 1999)
      REAL*8              :: ALPHA, CL

      ! ??
      REAL*8              :: P

      ! Index for layers
      INTEGER             :: I

      !-----------------------------------------------------------------
      ! Compute C_PPFD
      !-----------------------------------------------------------------  

      ! Initialize
      C_PPFD  = 0.d0

      ! 0.5 and 0.8 assume a spherical leaf-angle distribution
      KB      = 0.5d0 * CLUSTER / SUNCOS1
      KD      = 0.8d0 * CLUSTER
      P       = SQRT( 1.d0 - SCAT )
      REFLB   = 1.d0 - 
     &          EXP(-2.d0*KB * ( 1.d0-P ) / ( 1.d0+P ) / ( 1.d0+KB ) )
      KBP     = KB * P
      KDP     = KD * P

      ! 5-layer Gaussian integration over canopy
      DO I = 1, LAYERS 

         ! Cumulative LAI above layer I
         LAI_DEPTH  = LAI * DISTGAUSS( I )

         ! Fraction sun and shade leaves at layer I
         F_SUN      = EXP( -KB * LAI_DEPTH ) 
         F_SHADE    = 1.d0 - F_SUN

         ! For PAR on shaded leaves
         Q_SHADE_1  = Q_DIFF_2 * KDP * ( 1.d0 - REFLD ) *
     &                EXP( -KDP * LAI_DEPTH ) 

         ! For PAR on shaded leaves
         Q_SHADE_2  = Q_DIR_2 * ( KBP * ( 1.d0 - REFLB ) * 
     &                EXP( -KBP * LAI_DEPTH ) - KB * ( 1.d0 - SCAT ) * 
     &                EXP( -KB * LAI_DEPTH ) )

         ! PAR on shaded leaves
         Q_SHADE    = ( Q_SHADE_1 + Q_SHADE_2 ) / ( 1.d0 - SCAT )

         ! PAR on sunlit leaves
         Q_SUN      = Q_SHADE + KB * Q_DIR_2

         ! Update C_P_SUN and C_P_SHADE at layer I
         ! The following already accounts for canopy attenuation
         ! (Guenther et al, 1999)
         ALPHA      = 1.0d-3 + ( 8.5d-4 ) * LAI_DEPTH
         CL         = ( 1.42d0 ) * EXP( -( 3.0d-1 ) * LAI_DEPTH )

         C_P_SUN    = ALPHA * CL * Q_SUN / 
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SUN*Q_SUN )

         C_P_SHADE  = ALPHA * CL * Q_SHADE / 
     &                SQRT( 1.d0 + ALPHA*ALPHA * Q_SHADE*Q_SHADE )

         ! Update C_PPFD_I at layer I
         C_PPFD_I   = F_SUN * C_P_SUN + F_SHADE * C_P_SHADE

         ! Add on to total C_PPFD
         C_PPFD     = C_PPFD + WEIGHTGAUSS( I ) * C_PPFD_I * LAI

      ENDDO

      ! Prevent negative values.
      GAMMA_P = MAX( C_PPFD / NORMAL_FACTOR, 0d0 )

      END FUNCTION GET_GAMMA_P
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Leaf_Age
!
! !DESCRIPTION: Function Get\_Gamma\_Leaf\_Age computes the gamma exchange 
!  activity factor which is sensitive to leaf age (= Gamma\_Leaf\_Age).
!  Called from Get\_EmIsop\_Megan, Get\_EmMbo\_Megan, and Get\_EmMonog\_Megan.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_LEAF_AGE( CMLAI, PMLAI, T, SPECIES, TT ) 
     &         RESULT( GAMMA_LEAF_AGE )
!
! !INPUT PARAMETERS: 
!
      REAL*8,           INTENT(IN) :: T        ! Number of days between 
                                               !  current and previous LAI.
      REAL*8,           INTENT(IN) :: CMLAI    ! Current month's LAI [cm2/cm2]
      REAL*8,           INTENT(IN) :: PMLAI    ! Previous months LAI [cm2/cm2]
      CHARACTER(LEN=4), INTENT(IN) :: SPECIES  ! BVOC species name
      REAL*8,           INTENT(IN) :: TT       ! Daily average temperature [K]
!
! !RETURN VALUE:
!
      REAL*8                       :: GAMMA_LEAF_AGE   ! Activity factor    
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (3 ) Guenther et al, 2006
!  (5 ) Guenther et al, MEGAN user mannual 2007-08
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
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  13 Aug 2013 - M. Sulprizio- Updated for sesquiterpenes (H. Pye)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! M_AGE: leaf age factor
      REAL*8             :: M_AGE

      ! FNEW, ANEW: new foliage that emits negligible amounts of isoprene
      REAL*8             :: FNEW

      ! FGRO, AGRO: growing foliage that emits isoprene at low rates
      REAL*8             :: FGRO

      ! FMAT, AMAT: mature foliage that emits isoprene at peak rates
      REAL*8             :: FMAT

      ! FSEN, ASEN: senescing foliage that emits isoprene at reduced rates
      REAL*8             :: FSEN

      ! TI: number of days after budbreak required to induce iso. em.
      REAL*8             :: TI   

      ! TM: number of days after budbreak required to reach peak iso. em. rates
      REAL*8             :: TM

      ! Index for the relative emission activity (mpb,2008)
      INTEGER            :: AINDX
 
      ! Use species specific relative emission activity factors (mpb,2008)
      INTEGER, PARAMETER :: N_CAT = 5
      REAL*8             :: ANEW(N_CAT)
      REAL*8             :: AGRO(N_CAT)
      REAL*8             :: AMAT(N_CAT)
      REAL*8             :: ASEN(N_CAT)
      REAL*8             :: NORM_V(N_CAT)

      ! Constant Factors (not used)
      DATA    ANEW(  1),  AGRO(  1),  AMAT(  1),  ASEN(  1)
     &     /  1.0d0    ,  1.0d0    ,  1.0d0    ,  1.0d0      /
      ! Monoterpenes 
      DATA    ANEW(  2),  AGRO(  2),  AMAT(  2),  ASEN(  2)
     &     /  2.0d0   ,   1.8d0    ,  0.95d0   ,  1.0d0      /
      ! Isoprene and MBO 
      DATA    ANEW(  3),  AGRO(  3),  AMAT(  3),  ASEN(  3)
     &     /  0.05d0   ,  0.6d0    ,  1.125d0  ,  1.0d0      /
      ! Current set-up (v7-04-11) - used for other VOCs (OVOC)
      DATA    ANEW(  4),  AGRO(  4),  AMAT(  4),  ASEN(  4)
     &     /  0.01d0   ,  0.5d0    ,  1.00d0   ,  0.33d0     /
      ! SOAupdate: Sesquiterpenes hotp 3/2/10
      ! from MEGANv2.04/src/EMPROC/INCLDIR
      DATA    ANEW(  5),  AGRO(  5),  AMAT(  5),  ASEN(  5)
     &     /  0.4d0   ,  0.6d0    ,  1.075d0   ,  1.0d0     /

      ! Normalization factors
      ! SOAupdate: Sesq normalization factor added (hotp 3/16/2010)
      ! NORM_V(5) = 1/( 0.1*0.6 + 0.1*1 + 0.8*1.075 )
      DATA   NORM_V(1) , NORM_V(2) , NORM_V(3) , NORM_V(4), NORM_V(5)
     &     /                 1.0d0 , ! Constant
     &       0.96153846153846145d0 , ! Monoterpenes
     &       0.94339622641509424d0 , ! Isoprene & MBO
     &         1.132502831257078d0 , ! Other VOCs
     &            0.980392156863d0 / ! Sesquiterpenes (hotp 3/2/10,3/16/10)

      !=================================================================
      ! GET_GAMMA_LEAF_AGE begins here!
      !================================================================= 
      
      !-----------------------
      ! Compute TI and TM 
      ! (mpb,2009)
      !-----------------------

      IF ( TT <= 303.d0 ) THEN
         TI = 5.0d0 + 0.7 * ( 300.0d0 - TT )
      ELSEIF ( TT >  303.d0 ) THEN
         TI = 2.9d0
      ENDIF   
      TM = 2.3d0 * TI

      !-----------------------
      ! Compute M_AGE
      !-----------------------

      IF ( CMLAI == PMLAI ) THEN !(i.e. LAI stays the same) 

         !New values           Old Vlaues (mpb, 2009)
         FMAT = 0.8d0            !1.d0  
         FNEW = 0.d0             !0.d0
         FGRO = 0.1d0            !0.d0
         FSEN = 0.1d0            !0.d0

      ELSE IF ( CMLAI > PMLAI ) THEN !(i.e. LAI has INcreased) 

         IF ( T > TI ) THEN
            FNEW = ( TI / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE
            FNEW = 1.d0 - ( PMLAI / CMLAI )
         ENDIF

         IF ( T > TM ) THEN
            FMAT = ( PMLAI / CMLAI ) +
     &             ( ( T - TM ) / T ) * ( 1.d0 -  PMLAI / CMLAI )
         ELSE 
            FMAT = ( PMLAI / CMLAI )
         ENDIF

         ! Implement new condition for FGRO (mpb,2009)
         FSEN = 0.d0
         FGRO = 1.d0 - FNEW - FMAT

      ELSE ! This is the case if  PMLAI > CMLAI (i.e. LAI has DEcreased) 

         FSEN = ( PMLAI - CMLAI ) / PMLAI
         FMAT = 1.d0 - FSEN
         FGRO = 0.d0
         FNEW = 0.d0

      ENDIF

      ! Choose relative emission activity (mpb,2009)
      SELECT CASE ( TRIM(SPECIES) )
      CASE ('CONSTANT')
         AINDX = 1
      CASE ('MONO')
         AINDX = 2
      CASE ( 'ISOP','MBOT' )
         AINDX = 3
      CASE ('OVOC')
         AINDX = 4
      ! SOAupdate: Sesquiterpenes (hotp 3/2/10)
      CASE ('SESQ')
         AINDX = 5
      CASE DEFAULT
         GAMMA_LEAF_AGE = -999d0
         RETURN
      END SELECT

      ! Age factor
      M_AGE = FNEW*ANEW(AINDX) + FGRO*AGRO(AINDX) + 
     &        FMAT*AMAT(AINDX) + FSEN*ASEN(AINDX)  

      !Note: I think this should be normalized.
      !But, in the megan code I've consistently download 
      !from NCAR it never is...

      ! Normalize & prevent negative values
      GAMMA_LEAF_AGE = MAX( M_AGE * NORM_V(AINDX) , 0d0 )

      END FUNCTION GET_GAMMA_LEAF_AGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Gamma_Lai
!
! !DESCRIPTION: Function Get\_Gamma\_Lai computes the gamma exchange activity 
!  factor which is sensitive to leaf area (= GAMMA\_LAI).  Called from
!  Get\_EmIsop\_Megan, Get\_EmMbo\_Megan, and Get\_EmMonog\_Megan.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_GAMMA_LAI( CMLAI ) RESULT( GAMMA_LAI )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: CMLAI       ! Current month's LAI [cm2/cm2]
!
! !RETURN VALUE:
!
      REAL*8             :: GAMMA_LAI 
!
! !REMARKS:
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Guenther et al, 2006
!  (2 ) Guenther et al, MEGAN user mannual 2007-08
! 
! !REVISION HISTORY: 
!  (1 ) Original code by Dorian Abbot (9/2003).  Modified for the standard 
!        code by May Fu (11/2004)
!  (2 ) Update to publically released (as of 11/2004) MEGAN algorithm and 
!        modified for the standard code by May Fu (11/2004).
!  (3 ) Algorithm is based on the latest MEGAN v2.1 User's Guide (mpb,2009)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !-----------------------
      ! Compute GAMMA_LAI
      !-----------------------
      GAMMA_LAI = 0.49d0 * CMLAI / SQRT( 1.d0 + 0.2d0 * CMLAI*CMLAI )

      END FUNCTION GET_GAMMA_LAI 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE UPDATE_HCOT_DAY( HcoState, ExtState, I, J, FIRST, RC )
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER        :: HcoState 
      TYPE(Ext_State), POINTER        :: ExtState
      INTEGER,         INTENT(IN   )  :: I 
      INTEGER,         INTENT(IN   )  :: J 
      LOGICAL,         INTENT(IN   )  :: FIRST 
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

      ! On first call, set all memory values to current value
      ! Also fill the 15 day values with current values... 
      IF ( FIRST ) THEN
         HCOT_DAY(I,J,:)     = ExtState%TSURFK%Arr%Val(I,J)
         HCOPARDR_DAY(I,J,:) = ExtState%PARDR%Arr%Val(I,J)
         HCOPARDF_DAY(I,J,:) = ExtState%PARDF%Arr%Val(I,J)

         HCOT_15(I,J,:)      = ExtState%TSURFK%Arr%Val(I,J)
         HCOPARDR_15(I,J,:)  = ExtState%PARDR%Arr%Val(I,J) 
         HCOPARDF_15(I,J,:)  = ExtState%PARDF%Arr%Val(I,J)
      ELSE

         ! Move each day up
         DO D = DAY_DIM, 2, -1
            ! Need PAR as well as Temp  (mpb,2009)
            HCOT_DAY(I,J,D)     = HCOT_DAY(I,J,D-1)
            HCOPARDR_DAY(I,J,D) = HCOPARDR_DAY(I,J,D-1)
            HCOPARDF_DAY(I,J,D) = HCOPARDF_DAY(I,J,D-1)
         ENDDO
         
         ! Store 
         HCOT_DAY(I,J,1)     = ExtState%TSURFK%Arr%Val(I,J)
         HCOPARDF_DAY(I,J,1) = ExtState%PARDF%Arr%Val(I,J)
         HCOPARDR_DAY(I,J,1) = ExtState%PARDR%Arr%Val(I,J)
      ENDIF

      ! Return w/ success
      RC = HCO_SUCCESS
       
      END SUBROUTINE UPDATE_HCOT_DAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      REAL*8            :: D_DIM, D_DAYS, TMP_T
      REAL*8            :: TMP_PARDR , TMP_PARDF ! (mpb,2009)

      !=================================================================
      ! UPDATE_HCOT_15_AVG begins here!
      !=================================================================

      ! Convert to REAL*8
      D_DIM  = DBLE( DAY_DIM  )
      D_DAYS = DBLE( NUM_DAYS )

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
      HCOT_DAILY(I,J)  = TMP_T

      ! Repeat for PAR diffuse & direct (mpb,2009)
      HCOPARDR_15(I,J,1)   = TMP_PARDR
      HCOPARDR_15_AVG(I,J) = SUM( HCOPARDR_15(I,J,:) ) / D_DAYS 
      HCOPARDR_DAILY(I,J)  = TMP_PARDR

      HCOPARDF_15(I,J,1)   = TMP_PARDF
      HCOPARDF_15_AVG(I,J) = SUM( HCOPARDF_15(I,J,:) ) / D_DAYS 
      HCOPARDF_DAILY(I,J)  = TMP_PARDF

      ! Return w/ success 
      RC = HCO_SUCCESS

      END SUBROUTINE UPDATE_HCOT_15_AVG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Aef_Gen
!
! !DESCRIPTION: Subroutine Get\_Aef\_Gen creates AEFs for species without 
!  global explicit georeferenced emission factor maps. This includes 
!  sesquiterpenes (isoprene, MBO, and monoterpenes all have global AEF maps 
!  available).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_AEF_GEN( am_I_Root, HcoState, RC )
!
! !USES:
!
      USE HCO_EMISLIST_MOD,    ONLY : HCO_GetPtr
      USE EF_MGN20_MOD      ! EF_MGN20.EXT from MEGAN CDP, EFs
!
! !INPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN)     :: am_I_Root
      TYPE(HCO_State), POINTER        :: HcoState
!
! !INPUT/OUTPUT PARAMETERS
!
      INTEGER,         INTENT(INOUT)  :: RC
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) MEGAN CDP: MEGAN EFs originally from:
!       http://cdp.ucar.edu (search for MEGAN)
!       MEGAN_v2.04/src/MECHCONV/INCDIR/EF_MGN20.EXT
!       Files dated: 23 Aug 2007
!
!  Notes:
!  ============================================================================
!  (1 ) EFs are from EF_MGN20.EXT from the MEGAN CDP
!
! !REVISION HISTORY: 
!  25 Feb 2010 - H.O.T. Pye  - Wrote original code
!  13 Jul 2011 - M. Payer    - Added ProTeX headers
!  19 May 2014 - C. Keller   - Now a HEMCO extension: get plant functional
!                              types through HEMCO framework.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, SPECIES
      REAL*8               :: FACTOR
      REAL*8               :: SPECIES2CARBON
      CHARACTER(LEN=255)   :: MSG
      REAL*8               :: SUMFARN, SUMBCAR, SUMOSQT, SUMOMTP

      ! SOAupdate: Plant functional types (hotp 2/26/10)
      REAL(hp), POINTER       :: PFT_BT(:,:) => NULL()       ! broadleaf trees
      REAL(hp), POINTER       :: PFT_NT(:,:) => NULL()       ! needleleaf trees
      REAL(hp), POINTER       :: PFT_SH(:,:) => NULL()       ! shrubs
      REAL(hp), POINTER       :: PFT_GR(:,:) => NULL()       ! grasses
      REAL(hp), POINTER       :: PFT_CR(:,:) => NULL()       ! crops

      !=================================================================
      ! GET_AEF_GEN begins here!
      !=================================================================

      ! Init
      SUMFARN = 0.d0
      SUMBCAR = 0.d0
      SUMOSQT = 0.d0
      SUMOMTP = 0.d0

      ! Get pointers to data
      CALL HCO_GetPtr( am_I_Root, 'MEGAN_PFT_BT', PFT_BT, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'MEGAN_PFT_NT', PFT_NT, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'MEGAN_PFT_SH', PFT_SH, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'MEGAN_PFT_GR', PFT_GR, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetPtr( am_I_Root, 'MEGAN_PFT_CR', PFT_CR, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Calculation:
      ! AEF = PFT (m2/m2) * EF (ug/m2/h) * FACTOR
      ! FACOTR = kgC/kg * kg/ug * h/s
      ! AEF in kgC/m2/s

      ! Loop over the 20 MEGAN groups
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( SPECIES, SPECIES2CARBON )
!$OMP+PRIVATE( I, J, FACTOR )
      DO SPECIES = 1, 20

         ! Ratio of carbon mass to total compound mass [kgC/kg]
         IF ( SPECIES == IDMGFARN ) THEN
            SPECIES2CARBON = 15.d0 * 12.01d0 /
     &                       ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
         ELSEIF( SPECIES == IDMGBCAR ) THEN
            SPECIES2CARBON = 15.d0 * 12.01d0 /
     &                       ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
         ELSEIF( SPECIES == IDMGOSQT ) THEN
            SPECIES2CARBON = 15.d0 * 12.01d0 /
     &                       ( 15.d0 * 12.01d0 + 24.d0 * 1.01d0 )
         ! Other monoterpenes (hotp 3/5/10)
         ELSEIF( SPECIES == IDMGOMTP ) THEN
            SPECIES2CARBON = 10.d0 * 12.01d0 /
     &                       ( 10.d0 * 12.01d0 + 16.d0 * 1.01d0 )
         ELSE 
            ! species other than those listed above have not been
            ! checked so let's zero them
            SPECIES2CARBON = 0d0
         ENDIF

         ! Conversion factor from ug/m2/hr to 
         FACTOR = 1.d-9 / 3600.d0 * SPECIES2CARBON 

         ! Compute AEF and convert from [ug/m2/hr] to [kgC/m2/s]
         ! Loop over surface boxes
         DO J = 1, HcoState%NY 
         DO I = 1, HcoState%NX

            ! Calculate AEF in [ug/m2/h]
            AEF_GEN( I, J, SPECIES ) =
     &                        PFT_BT( I, J ) * EF_BT( SPECIES ) + 
     &                        PFT_NT( I, J ) * EF_NT( SPECIES ) +
     &                        PFT_SH( I, J ) * EF_SB( SPECIES ) +
     &                        PFT_GR( I, J ) * EF_HB( SPECIES ) +
     &                        PFT_CR( I, J ) * EF_HB( SPECIES )

            ! Convert AEF to [kgC/m2/s]
            AEF_GEN( I, J, SPECIES ) = AEF_GEN( I, J, SPECIES ) *
     &                                 FACTOR 

            ! Eventually write out total sum (TgC/yr)
            IF ( SPECIES == IDMGFARN ) THEN
               SUMFARN = SUMFARN +
     &                 ( AEF_GEN(I,J,SPECIES) * 1d-9 * 365d0 * 24d0 *
     &                   3600d0 * HcoState%Grid%AREA_M2%Val(I,J) )
            ELSEIF( SPECIES == IDMGBCAR ) THEN
               SUMBCAR = SUMBCAR +
     &                 ( AEF_GEN(I,J,SPECIES) * 1d-9 * 365d0 * 24d0 *
     &                   3600d0 * HcoState%Grid%AREA_M2%Val(I,J) )
            ELSEIF( SPECIES == IDMGOSQT ) THEN
               SUMOSQT = SUMOSQT +
     &                 ( AEF_GEN(I,J,SPECIES) * 1d-9 * 365d0 * 24d0 *
     &                   3600d0 * HcoState%Grid%AREA_M2%Val(I,J) )
            ELSEIF( SPECIES == IDMGOMTP ) THEN
               SUMOMTP = SUMOMTP +
     &                 ( AEF_GEN(I,J,SPECIES) * 1d-9 * 365d0 * 24d0 *
     &                   3600d0 * HcoState%Grid%AREA_M2%Val(I,J) )
            ENDIF
         ENDDO
         ENDDO

      ENDDO       ! species
!$OMP END PARALLEL DO

      ! Diagnostic info to log/screen
      write(MSG,*) 'Annual Sesquiterpene Base AEFs [TgC/yr]:'
      CALL HCO_MSG(MSG)
      write(MSG,*) '   FARN', SUMFARN
      CALL HCO_MSG(MSG)
      write(MSG,*) '   BCAR', SUMBCAR
      CALL HCO_MSG(MSG)
      write(MSG,*) '   OSQT', SUMOSQT
      CALL HCO_MSG(MSG)
      write(MSG,*) '   OMTP', SUMOMTP
      CALL HCO_MSG(MSG)

      ! Nullify pointers
      PFT_BT => NULL()
      PFT_NT => NULL()
      PFT_SH => NULL()
      PFT_GR => NULL()
      PFT_CR => NULL()

      ! Return to calling program
      RC = HCO_SUCCESS

      END SUBROUTINE GET_AEF_GEN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                        :: nSpc, AS, NX, NY
      INTEGER, ALLOCATABLE           :: HcoIDs(:)
      REAL*8                         :: PI_180
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
      CALL GetExtOpt ( ExtNr, 'PECCA model', OptValBool=LPECCA, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL GetExtOpt ( ExtNr, 'Isoprene scaling', 
     &                 OptValDp=ISOP_SCALING, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Get species IDs
      ! --> Assume that species are ordered ISOP, ACET, PRPE, C2H4 in
      !     config. file!
      CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
      IF ( RC /= HCO_SUCCESS ) RETURN
      IF ( nSpc /= 4 ) THEN
         MSG = 'Expected four species for extension ' // TRIM(ExtName)
         CALL HCO_ERROR ( MSG, RC )
         RETURN
      ENDIF
      IDTISOP = HcoIDs(1)
      IDTACET = HcoIDs(2)
      IDTPRPE = HcoIDs(3)
      IDTC2H4 = HcoIDs(4)

      ! Verbose mode
      MSG = 'Use MEGAN biogenic emissions (extension module)'
      CALL HCO_MSG( MSG, SEP1='-' )
      WRITE(MSG,*) '   - Use PECCA model: ', LPECCA
      CALL HCO_MSG( MSG )
      MSG = '   - Use the following species:'
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '     Isoprene   = ', TRIM(SpcNames(1)), IDTISOP 
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '     --> Isoprene scale factor is ', ISOP_SCALING
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '     Acetone    = ', TRIM(SpcNames(2)), IDTACET 
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '     C3 Alkenes = ', TRIM(SpcNames(3)), IDTPRPE 
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '     Ethene     = ', TRIM(SpcNames(4)), IDTC2H4
      CALL HCO_MSG( MSG )

      ! Check for monoterpene option --> this one is required for 
      ! CO, OCPI, MONX
      ExtNrMono = GetExtNr('MEGAN_Mono')
      IF ( ExtNrMono > 0 ) THEN

         ! Get HEMCO species IDs of CO, OCPI and MONX
         CALL HCO_GetExtHcoID( HcoState, ExtNrMono, HcoIDs,
     &                         SpcNames, nSpc,      RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN
         IF ( nSpc /= 3 ) THEN
            MSG = 'Expected three species for MEGAN_Mono extension' 
            CALL HCO_ERROR ( MSG, RC )
            RETURN
         ENDIF
         IDTCO   = HcoIDs(1)
         IDTOCPI = HcoIDs(2)
         IDTMONX = HcoIDs(3)

         ! Verbose
         MSG = '   - MEGAN monoterpene option enabled:' 
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     CO         = ', TRIM(SpcNames(1)), IDTCO
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     OC aerosol = ', TRIM(SpcNames(2)), IDTOCPI
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     Monoterp.  = ', TRIM(SpcNames(3)), IDTMONX
         CALL HCO_MSG( MSG )
      ELSE
         IDTCO      = -1 
         IDTOCPI    = -1 
         IDTMONX    = -1 

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
      ! This also adds species ALPH/LIMO/ALCO
      ExtNrSoa = GetExtNr('MEGAN_SOA')
      IF ( ExtNrSoa > 0 ) THEN
         IF ( ExtNrMono <= 0 ) THEN
            MSG = 'MEGAN SOA option only valid with monoterpene!'
            CALL HCO_ERROR ( MSG, RC )
            RETURN
         ENDIF

         ! Get HEMCO species IDs of ALPH/LIMO/ALCO 
         CALL HCO_GetExtHcoID( HcoState, ExtNrSoa, HcoIDs,
     &                         SpcNames, nSpc,     RC       )
         IF ( RC /= HCO_SUCCESS ) RETURN
         IF ( nSpc /= 4 ) THEN
            MSG = 'Expected four species for MEGAN_SOA extension' 
            CALL HCO_ERROR ( MSG, RC )
            RETURN
         ENDIF
         IDTMTPA = HcoIDs(1)
         IDTMTPO = HcoIDs(2)
         IDTLIMO = HcoIDs(3)
         IDTSESQ = HcoIDs(4)

         ! Verbose
         MSG = '   - SOA mechanism enabled:' 
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     a-,b-pinene    = ',TRIM(SpcNames(1)),IDTMTPA
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     Other monoterp.= ',TRIM(SpcNames(2)),IDTMTPO
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     Limonene       = ',TRIM(SpcNames(3)),IDTLIMO
         CALL HCO_MSG( MSG )
         WRITE(MSG,*) '     Sesquiterpenes = ',TRIM(SpcNames(4)),IDTSESQ
         CALL HCO_MSG( MSG )
      ELSE
         IDTMTPA = -1
         IDTMTPO = -1
         IDTLIMO = -1
         IDTSESQ = -1
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
      HCOT_DAY = 0d0

      ALLOCATE( HCOT_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_15', RC )
        RETURN
      ENDIF
      HCOT_15 = 0d0

      ALLOCATE( HCOT_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_15_AVG', RC )
        RETURN
      ENDIF
      HCOT_15_AVG = 0d0

      ! Daily averaged temperature (mpb,2009)
      ALLOCATE( HCOT_DAILY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOT_DAILY', RC )
        RETURN
      ENDIF
      HCOT_DAILY = 0d0

      ! Allocate arrays for light (mpb,2009)
      ! -- Direct --
      ALLOCATE( HCOPARDR_DAY( NX, NY, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_DAY', RC )
        RETURN
      ENDIF
      HCOPARDR_DAY = 0d0

      ALLOCATE( HCOPARDR_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_15', RC )
        RETURN
      ENDIF
      HCOPARDR_15 = 0d0

      ALLOCATE( HCOPARDR_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_15_AVG', RC )
        RETURN
      ENDIF
      HCOPARDR_15_AVG = 0d0

      ALLOCATE( HCOPARDR_DAILY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDR_DAILY', RC )
        RETURN
      ENDIF
      HCOPARDR_DAILY = 0d0

      ! -- Diffuse --
      ALLOCATE( HCOPARDF_DAY( NX, NY, DAY_DIM ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_DAY', RC )
        RETURN
      ENDIF
      HCOPARDF_DAY = 0d0

      ALLOCATE( HCOPARDF_15( NX, NY, NUM_DAYS ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_15', RC )
        RETURN
      ENDIF
      HCOPARDF_15 = 0d0

      ALLOCATE( HCOPARDF_15_AVG( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_15_AVG', RC )
        RETURN
      ENDIF
      HCOPARDF_15_AVG = 0d0

      ALLOCATE( HCOPARDF_DAILY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'HCOPARDF_DAILY', RC )
        RETURN
      ENDIF
      HCOPARDF_DAILY = 0d0

      ! -- Generic emission factors --
      ! SOAupdate: For all 20 MEGAN groups (hotp 2/27/10)
      ALLOCATE( AEF_GEN( NX, NY, 20 ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'AEF_GEN', RC )
        RETURN
      ENDIF
      AEF_GEN = 0.0_hp

      ! LAI from previous time step
      ALLOCATE( LAI_PREVSTEP( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'LAI_PREVSTEP', RC )
        RETURN
      ENDIF
      LAI_PREVSTEP = 0d0

      ALLOCATE( LAI_PREVDAY( NX, NY ), STAT=AS )
      IF ( AS /= 0 ) THEN 
        CALL HCO_ERROR( 'LAI_PREVDAY', RC )
        RETURN
      ENDIF
      LAI_PREVDAY = 0d0

      !=================================================================
      ! The original MEGAN code used to read the emission factors here.
      ! We now get the emisson factors through HEMCO, hence no need to
      ! do this anymore (ckeller, 05/19/2014).
      !=================================================================

      ! Set physical constants
      PI_180 = HcoState%Phys%PI / 180.0d0
      D2RAD  = PI_180               ! Degrees to radians
      RAD2D  = 1.0d0 / PI_180       ! Radians to degrees

      ! Enable met. fields
      ExtState%TSURFK%DoUse    = .TRUE.
      ExtState%SUNCOSmid%DoUse = .TRUE.
      ExtState%PARDR%DoUse     = .TRUE.
      ExtState%PARDF%DoUse     = .TRUE.
      ExtState%SZAFACT%DoUse   = .TRUE.
      ExtState%GC_LAI%DoUse    = .TRUE.

      ! Enable module 
      ExtState%Megan = .TRUE.

      ! Leave w/ success
      IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
      IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
      CALL HCO_LEAVE( RC )

      END SUBROUTINE HCOX_Megan_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE HCOX_MEGAN_FINAL
! 
! !REVISION HISTORY: 
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  16 Aug 2013 - C. Keller   - Now a HEMCO extension 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! HCOX_MEGAN_FINAL begins here!
      !=================================================================

      ! Nullify pointers
      GEIA_ORVC => NULL()
      AEF_ISOP  => NULL()
      AEF_MONOT => NULL()
      AEF_MBO   => NULL()
      AEF_APINE => NULL()
      AEF_BPINE => NULL()
      AEF_LIMON => NULL()
      AEF_SABIN => NULL()
      AEF_MYRCN => NULL()
      AEF_CAREN => NULL()
      AEF_OCIMN => NULL()
      AEF_ACET  => NULL()

      ! Cleanup module arrays
      IF ( ALLOCATED( HCOT_DAY        ) ) DEALLOCATE( HCOT_DAY        )
      IF ( ALLOCATED( HCOT_15         ) ) DEALLOCATE( HCOT_15         )
      IF ( ALLOCATED( HCOT_15_AVG     ) ) DEALLOCATE( HCOT_15_AVG     )
      IF ( ALLOCATED( HCOT_DAILY      ) ) DEALLOCATE( HCOT_DAILY      )
      IF ( ALLOCATED( HCOPARDR_DAY    ) ) DEALLOCATE( HCOPARDR_DAY    )
      IF ( ALLOCATED( HCOPARDR_15     ) ) DEALLOCATE( HCOPARDR_15     )
      IF ( ALLOCATED( HCOPARDR_15_AVG ) ) DEALLOCATE( HCOPARDR_15_AVG )
      IF ( ALLOCATED( HCOPARDR_DAILY  ) ) DEALLOCATE( HCOPARDR_DAILY  )
      IF ( ALLOCATED( HCOPARDF_DAY    ) ) DEALLOCATE( HCOPARDF_DAY    )
      IF ( ALLOCATED( HCOPARDF_15     ) ) DEALLOCATE( HCOPARDF_15     )
      IF ( ALLOCATED( HCOPARDF_15_AVG ) ) DEALLOCATE( HCOPARDF_15_AVG )
      IF ( ALLOCATED( HCOPARDF_DAILY  ) ) DEALLOCATE( HCOPARDF_DAILY  )
      IF ( ALLOCATED( LAI_PREVSTEP    ) ) DEALLOCATE( LAI_PREVSTEP    )
      IF ( ALLOCATED( LAI_PREVDAY     ) ) DEALLOCATE( LAI_PREVDAY     )

      END SUBROUTINE HCOX_MEGAN_FINAL
!EOC
      END MODULE HCOX_MEGAN_MOD
