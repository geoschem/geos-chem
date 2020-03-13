!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ucx_mod.F90
!
! !DESCRIPTION: Module UCX\_MOD contains routines and variables which
!  are associated with the addition of full stratospheric chemistry to
!  GEOS-Chem (based on the NASA GMI implementation, forming the Unified
!  Chemistry eXtension (UCX).
!\\
!\\
! !INTERFACE:
!
MODULE UCX_MOD
!
! !USES:
!
  USE ERROR_MOD,     ONLY : DEBUG_MSG
  USE inquireMod,    ONLY : findFreeLUN
  USE PhysConstants       ! Physical constants
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

  !#if !defined(ESMF_)
  ! NcdfUtil modules for netCDF I/O
  USE m_netcdf_io_open                    ! netCDF open
  USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
  USE m_netcdf_io_read                    ! netCDF data reads
  USE m_netcdf_io_close                   ! netCDF close
  !#endif

  IMPLICIT NONE
  !#if !defined(ESMF_)
# include "netcdf.inc"
  !#endif

  PRIVATE

!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC :: NOON_FILE_ROOT ! Directory for noontime data
  PUBLIC :: T_STS          ! Max temperature of STS formation (K)
  PUBLIC :: NDENS_AER      ! See below
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SET_H2O_TRAC
  PUBLIC  :: SETTLE_STRAT_AER
  PUBLIC  :: SO4_PHOTFRAC
  PUBLIC  :: UCX_NOX
  PUBLIC  :: UCX_H2SO4PHOT
  PUBLIC  :: CALC_STRAT_AER
  PUBLIC  :: GET_STRAT_OPT
  PUBLIC  :: KG_STRAT_AER
  PUBLIC  :: RHO_STRAT_AER
  PUBLIC  :: INIT_UCX
  PUBLIC  :: CLEANUP_UCX
!
! PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: TERNARY
  PRIVATE :: CARSLAW_DENSITY
  PRIVATE :: CALC_H2SO4_GAS
  PRIVATE :: CALC_SLA_GAMMA
  PRIVATE :: MOLEC_SPEED
  PRIVATE :: NOXCOEFF_INIT
  PRIVATE :: GET_JJNOX
!
! !REVISION HISTORY:
!  26 Mar 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  !=================================================================
  ! MODULE PARAMETERS
  !
  ! UCX_NLEVS       : Number of levels in AER data
  ! UCX_NLAT        : Number of latitudes in AER data
  ! T_STS           : Maximum temperature of STS formation (K)
  ! I_SLA           : Index of liquid aerosols
  ! I_SPA           : Index of particulate PSCs
  ! INITMR_BASIS    : Year for which the initializing mixing ratios
  !                   were calculated (needed for future-scaling)
  ! UCXNETCDF       : Read data from NetCDF
  !
  !=================================================================

  INTEGER,  PARAMETER           :: UCX_NLEVS=51
  INTEGER,  PARAMETER           :: UCX_NLAT=19
  REAL(fp), PARAMETER           :: T_STS=240.0e+0_fp
  INTEGER,  PARAMETER           :: I_SLA=1
  INTEGER,  PARAMETER           :: I_SPA=2
  INTEGER,  PARAMETER           :: INITMR_BASIS = 2005

#if defined( ESMF_ ) || defined( MODEL_WRF )
  ! Never use NETCDF in ESMF environment (ckeller, 12/05/14).
  ! Don't use UCXNETCDF for WRF-GC as NetCDF NOx coeffs not
  ! preprocessed for external grid (hplin, 8/15/18).
  LOGICAL, PARAMETER            :: UCXNETCDF = .FALSE.
#else
  LOGICAL, PARAMETER            :: UCXNETCDF = .TRUE.
#endif
!
! PRIVATE TYPES:
!
  !=================================================================
  ! MODULE VARIABLES:
  !
  ! Scalars
  !
  ! TRAC_IDX           : Species index for output
  ! SLA_VA             : SLA volume-area conversion
  ! SLA_RR             : SLA effective-liquid radius conversion
  ! SLA_VR             : SLA volume-effective radius conversion
  ! NATMW              : Molar mass of NAT (kg/kmol)
  ! ICEMEW             : Molar mass of ice (kg/kmol)
  ! DENSNAT            : Density of pure NAT (kg/m3)
  ! DENSICE            : Density of pure ice (kg/m3)
  ! ISR_ClNO3          : ClNO3 MW (inverse sqrt) (kg/kmol)^-0.5
  ! ISR_BrNO3          : BrNO3 MW (inverse sqrt) (kg/kmol)^-0.5
  ! ISR_N2O5           : N2O5 MW (inverse sqrt) (kg/kmol)^-0.5
  ! ISR_HOCl           : HOCl MW (inverse sqrt) (kg/kmol)^-0.5
  ! ISR_HOBr           : HOBr MW (inverse sqrt) (kg/kmol)^-0.5
  !
  ! Arrays
  !
  ! UCX_PLEVS          : Pressure levels of 2D data (hPa)
  ! UCX_LATS           : Latitude edges of 2D data (deg)
  ! RAD_AER            : Strat. aerosol radius (cm)
  ! KG_AER             : Aerosol mass (kg/box)
  ! SAD_AER            : Aerosol surface area density (cm2/cm3)
  ! NDENS_AER          : Aerosol number density (#/m3)
  ! RHO_AER            : Aerosol mass density (kg/m3 aerosol)
  ! AERFRAC            : Mass fraction of species in liquid aerosols
  ! AERFRACIND         : Indices of liquid aerosol species
  ! NOX_O              : Monthly mean noontime O3P/O1D for NOx calcs
  ! NOX_J              : Monthly mean noontime J-rates for NOx calcs
  ! SO4_TOPPHOT        : Photolysis rate at the top of the chemgrid (1/s)
  !
  !=================================================================

  ! Scalars
  CHARACTER(LEN=255)   :: NOON_FILE_ROOT
  REAL(fp)             :: SLA_VA
  REAL(fp)             :: SLA_RR
  REAL(fp)             :: SLA_VR
  REAL(fp), PARAMETER  :: NATMW   = 117.0
  REAL(fp), PARAMETER  :: ICEMW   = 18.0
  REAL(fp), PARAMETER  :: DENSNAT = 1626.e+0_fp
  REAL(fp), PARAMETER  :: DENSICE = 990.0e+0_fp
  REAL(fp), PARAMETER  :: ISR_ClNO3=1.e+0_fp/sqrt(97.46e+0_fp)
  REAL(fp), PARAMETER  :: ISR_BrNO3=1.e+0_fp/sqrt(141.9e+0_fp)
  REAL(fp), PARAMETER  :: ISR_N2O5 =1.e+0_fp/sqrt(108.0e+0_fp)
  REAL(fp), PARAMETER  :: ISR_HOCl =1.e+0_fp/sqrt(52.46e+0_fp)
  REAL(fp), PARAMETER  :: ISR_HOBr =1.e+0_fp/sqrt(96.91e+0_fp)

  ! Arrays
  REAL(fp),DIMENSION(:,:),ALLOCATABLE     :: UCX_REGRID
  REAL(fp),DIMENSION(:),ALLOCATABLE       :: UCX_PLEVS
  REAL(fp),DIMENSION(:),ALLOCATABLE       :: UCX_LATS
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: RAD_AER
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: KG_AER
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: SAD_AER
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NDENS_AER
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: RHO_AER
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: AERFRAC
  INTEGER,DIMENSION(:),ALLOCATABLE      :: AERFRACIND
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_O
  REAL(fp),DIMENSION(:,:,:,:),ALLOCATABLE :: NOX_J
  REAL(fp),DIMENSION(:,:),ALLOCATABLE     :: SO4_TOPPHOT

  !=================================================================
  ! Variables to use NOx coefficients in ESMF / grid-independent
  ! envionment. The NOx coefficients are climatological 2D
  ! (lat/lev/12 months) data that are currently available for
  ! horizontal (latitude) resolutions of 2 and 4 degrees. For other
  ! resolutions, the horizontal data becomes mapped onto the
  ! simulation grid (see GET_JJNOX).
  ! Similar to the surface mixing ratio boundary conditions, we now
  ! read all the NOx coefficients during initialization to avoid
  ! additional I/O calls during run time (only if UCXNETCDF = false)
  ! (ckeller, 05/12/2014).
  !=================================================================
  REAL(fp), ALLOCATABLE, TARGET :: NOXCOEFF(:,:,:,:)
  REAL(fp), ALLOCATABLE         :: NOXLAT(:)
  INTEGER                     :: JJNOXCOEFF

  !=================================================================
  ! Species ID flags
  ! These are now defined in INIT_UCX and used where needed.
  ! (sde, bmy, 6/21/16)
  !=================================================================
  INTEGER :: id_BCPI,    id_Br,       id_Br2,    id_BrCl,  id_BrNO2
  INTEGER :: id_BrNO3,   id_BrO,      id_CCl4
  INTEGER :: id_CH3Br,   id_CHBr3,    id_CH2Br2, id_CH3Cl
  INTEGER :: id_CH3CCl3, id_CH4,      id_Cl,     id_Cl2,   id_Cl2O2
  INTEGER :: id_ClNO2,   id_ClNO3,    id_ClO,    id_ClOO,  id_H1211
  INTEGER :: id_H1301,   id_H2,       id_H2402,  id_H2O,   id_HBr
  INTEGER ::             id_HCl,      id_HNO2,   id_HNO3,  id_N
  INTEGER :: id_HNO4,    id_HOBr,     id_HOCl,   id_N2O,   id_N2O5
  INTEGER :: id_NIT,     id_NO,       id_NO2,    id_NO3,   id_O3
  INTEGER :: id_OClO,    id_PAN,      id_SO2,   id_SO4

CONTAINS
!
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_nox
!
! !DESCRIPTION: Subroutine UCX\_NOX calculates NOx and N2O loss
!  rates above the chemistry grid, based on estimates of j-rates from
!  a 2D model and simple photochemical assumptiones.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE UCX_NOX( Input_Opt, State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : ITS_A_LEAPYEAR
    USE TIME_MOD,           ONLY : GET_HOUR
    USE TIME_MOD,           ONLY : GET_LOCALTIME
    USE TIME_MOD,           ONLY : GET_MINUTE
    USE CMN_FJX_MOD,        ONLY : ZPJ
    USE FAST_JX_MOD,        ONLY : RXN_NO, RXN_NO2, RXN_NO3, RXN_N2O
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  16 Jul 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! For indexing the O1D/O3P concentration array
    INTEGER, PARAMETER :: O3PIDX = 1
    INTEGER, PARAMETER :: O1DIDX = 2

    ! For indexing the J-rate coefficient array
    INTEGER, PARAMETER :: JNOIDX = 1
    INTEGER, PARAMETER :: JNO2IDX= 2
    INTEGER, PARAMETER :: JNO3IDX= 3
    INTEGER, PARAMETER :: JN2OIDX= 4

    ! Reaction rates and indices
    REAL(fp),DIMENSION(12):: RRATE
    INTEGER, PARAMETER :: k_JNO = 6
    INTEGER, PARAMETER :: k_JNO2 = 4
    INTEGER, PARAMETER :: k_JNO3 = 5
    INTEGER, PARAMETER :: k_JN2O = 12

    ! Intermediate variables
    REAL(fp) :: LOCALNOX, LOCALN2O, LOCALO3, LOCALO1D, LOCALO3P
    REAL(fp) :: OLD_NO3, OLD_NO2, OLD_NO, OLD_N, OLD_N2O
    REAL(fp) :: NEW_NO3, NEW_NO2, NEW_NO, NEW_N, NEW_N2O
    REAL(fp) :: FRACNO, FRACNO2, FRACNO3, FRACN
    REAL(fp) :: NOXRATE, N2ORATE, KGNOX, KGN2O
    REAL(fp) :: NO_ALPHA, NO_BETA, NO_GAMMA, NO_EPSILON
    REAL(fp) :: NO_QA, NO_QX, NO_QC, DNOX, DN2O
    REAL(fp) :: MESONOX_DELTA
    REAL(fp) :: MESON2O_DELTA

    ! Local air number density (molec/cm3) and box mass (kg)
    REAL(fp) :: NDAIR
    REAL(fp) :: XAIR

    ! Local temperature (K) and inverted T (1/K)
    REAL(fp) :: T3K, TINV

    ! Chemistry timestep in seconds
    REAL(f8) :: DTCHEM

    ! Timestep in hours
    REAL(f8) :: DTCHEM_HR

    ! Timing information
    INTEGER, SAVE :: LASTMONTH=0

    ! Grid indexing
    INTEGER :: I,J,L
    REAL(fp)  :: MIDLAT, ZDEL, ZBASE
    REAL(fp)  :: ZMID(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Local daylight fraction
    REAL(fp)             :: DAYFRAC
    LOGICAL              :: CYCLEBOX
    CHARACTER(LEN=255)   :: DBGMSG

    ! Local variables for quantities from Input_Opt
    LOGICAL              :: prtDebug

    ! Pointers
    REAL(fp), POINTER    :: Spc (:,:,:,:)

    ! Required for updated chemistry
    Integer            :: LMinPhot

    !=================================================================
    ! UCX_NOX begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Initialize GEOS-Chem species array [kg]
    Spc => State_Chm%Species

    ! Retrieve monthly mean data if necessary
    IF (LASTMONTH.ne.GET_MONTH()) THEN
       LASTMONTH = GET_MONTH()
       CALL GET_NOXCOEFF( LASTMONTH, Input_Opt, State_Grid, State_Met)
    ENDIF

    ! Get chemistry step length in seconds
    DTCHEM = GET_TS_CHEM()
    DTCHEM_HR = DTCHEM/3600.e+0_f8

    ! Reset NOx/N2O mass counters
    MESONOX_DELTA = 0e+0_fp
    MESON2O_DELTA = 0e+0_fp

    ! First compute ZMID outside of main parallel loop
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, ZBASE, ZDEL )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Zero base height
       ZBASE = 0e+0_fp

       ! Compute the array of midpoint heights [m]
       DO L = 1, State_Grid%NZ
          ZDEL        = State_Met%BXHEIGHT(I,J,L)
          ZMID(I,J,L) = ZBASE + (ZDEL/2.e+0_fp)
          ZBASE       = ZBASE + ZDEL
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Main parallel DO loop over lon, lat, alt
    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I,        J,        L                       ) &
    !$OMP PRIVATE( DAYFRAC,            DN2O,     CYCLEBOX      ) &
    !$OMP PRIVATE(           RRATE,    T3K,      TINV          ) &
    !$OMP PRIVATE( NDAIR,    XAIR,     KGNOX,    KGN2O         ) &
    !$OMP PRIVATE( LOCALNOX, LOCALN2O, LOCALO3,  LOCALO3P      ) &
    !$OMP PRIVATE( LOCALO1D, NO_ALPHA, NO_BETA,  NO_GAMMA      ) &
    !$OMP PRIVATE( FRACNO2,  FRACNO3,  FRACNO,   FRACN         ) &
    !$OMP PRIVATE( NOXRATE,  N2ORATE,  DNOX                    ) &
    !$OMP PRIVATE( OLD_NO3,  OLD_NO2,  OLD_NO,   OLD_N         ) &
    !$OMP PRIVATE( NEW_NO3,  NEW_NO2,  NEW_NO,   NEW_N         ) &
    !$OMP PRIVATE( OLD_N2O,  NEW_N2O,  LMinPhot                )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Calculate daylight fraction
       DAYFRAC = State_Met%SUNCOSmid(I,J)
       LMINPHOT  = State_Met%ChemGridLev(I,J)
       CycleBox = (L.le.LMinPhot)

       ! Are the pre-calculated J-rates zero (ie in darkness)?
       IF (CYCLEBOX) CYCLE

       ! Reset reaction rates
       RRATE = 0e+0_fp

       ! Retrieve air mass (kg) and local temperature (K)
       T3K   = State_Met%T(I,J,L)
       TINV  = 1.e+0_fp/T3K

       ! Calculate air number density (molec/cm3)
       ! Now using dry air partial pressure (ewl, 3/2/15)
       NDAIR = State_Met%PMID_DRY(I,J,L)*1.d-4*AVO/(RSTARG*T3K)
       XAIR  = NDAIR / State_Met%AD(I,J,L)

       ! Get total mass of N2O
       kgN2O = Spc(I,J,L,id_N2O)

       ! Get local concentrations in molec/cm3
       ! Ignore the data for atomic N - not well defined for
       ! non-chemgrid boxes, and in any case the molar mass
       ! may not be present or correctly defined
       OLD_N   = 0.0_fp
       OLD_NO  = Spc(I,J,L,id_NO) * &
                 (AIRMW / State_Chm%SpcData(id_NO )%Info%emMW_g) * XAIR
       OLD_NO2 = Spc(I,J,L,id_NO2) * &
                 (AIRMW / State_Chm%SpcData(id_NO2)%Info%emMW_g) * XAIR
       OLD_NO3 = Spc(I,J,L,id_NO3) * &
                 (AIRMW / State_Chm%SpcData(id_NO3)%Info%emMW_g) * XAIR
       OLD_N2O = kgN2O * &
                 (AIRMW / State_Chm%SpcData(id_N2O)%Info%emMW_g) * XAIR

       ! Total concentrations
       localNOx = OLD_N + OLD_NO + OLD_NO2 + OLD_NO3
       localN2O = OLD_N2O

       ! Get ozone from model
       LOCALO3  = Spc(I,J,L,id_O3) *  ( AIRMW / &
                  State_Chm%SpcData(id_O3)%Info%emMW_g ) * XAIR

       ! These reactions are relevant during both day and night-time
       ! chemistry
       ! 2:  NO + O3 -> NO2 + O2
       RRATE(2)  = 3.0e-12_fp*exp(-1500.e+0_fp*TINV)
       ! 3:  NO2 + O3 -> NO3 + O2
       RRATE(3)  = 1.2e-13_fp*exp(-2450.e+0_fp*TINV)

       ! During night time, assume everything is NO2 and
       ! that there is no significant production or loss.
       ! The only exception is NO3, which is assumed to
       ! remain unchanged
       If (DayFrac.lt.1.0e-10_fp) Then
          ! Simple assumption: all N and NO go to NO2
          ! NO3 remains unchanged
          DN2O = 0.0e+0_fp
          DNOx = 0.0e+0_fp
          ! Complex assumption: all N goes to NO
          ! Night-time chemistry is slow; assume zero
          ! O1D or O3P.
          ! Simplified reactions - allow slow conversion of NO -> NO2
          ! and NO2 -> NO3. Ignore all others. This can be solved
          ! analytically (decay chain)
          NEW_N   = 0.0e+0_fp
          NEW_NO  = (OLD_NO + OLD_N)*dexp(-dtChem*RRATE(2)*localO3)
          NEW_NO2 = OLD_NO2*dexp(-dtChem*RRATE(3)*localO3) + &
                    ((OLD_NO + OLD_N)/((RRATE(3)/RRATE(2))-1.0e+0_fp)) &
                    *(dexp(-dtChem*RRATE(2)*localO3) - &
                      dexp(-dtChem*RRATE(3)*localO3))
          NEW_NO3 = OLD_NO3 + (OLD_NO2 + OLD_NO + OLD_N) - &
                              (NEW_NO2 + NEW_NO + NEW_N)
          FRACN   = 0.0e+0_fp
          FRACNO  = NEW_NO/localNOx
          FRACNO2 = NEW_NO2/localNOx
          FRACNO3 = NEW_NO3/localNOx
       Else
          ! Calculate remaining rate constants
          ! 1:  NO2 + O -> NO + O2
          RRATE(1)  = 5.1e-12_fp*exp(210.e+0_fp*TINV)
          ! 4:  NO2 + hv -> NO + O1D
          !RRATE(k_JNO2) = NOX_J(I,J,L,JNO2IDX)*DAYFRAC
          RRATE(k_JNO2) = ZPJ(LMINPHOT,RXN_NO2,I,J)
          ! 5:  NO3 + hv -> NO2 + O
          !RRATE(k_JNO3) = NOX_J(I,J,L,JNO3IDX)*DAYFRAC
          RRATE(k_JNO3) = ZPJ(LMINPHOT,RXN_NO3,I,J)
          ! 6:  NO + hv -> N + O
          !RRATE(k_JNO ) = NOX_J(I,J,L,JNOIDX)*DAYFRAC
          RRATE(k_JNO) = ZPJ(LMINPHOT,RXN_NO,I,J)
          ! 7:  N + NO2 -> N2O + O
          RRATE(7)  = 5.8e-12_fp*exp(220.e+0_fp*TINV)
          ! 8:  N + NO -> N2 + O
          RRATE(8)  = 2.1e-11_fp*exp(100.e+0_fp*TINV)
          ! 9:  N + O2 -> NO + O
          RRATE(9)  = 1.5e-11_fp*exp(-3600.e+0_fp*TINV)
          ! 10:  N2O + O1D -> N2 + O2
          RRATE(10) = 4.63e-11_fp*exp(20.e+0_fp*TINV)
          ! 11:  N2O + O1D -> 2NO
          RRATE(11) = 7.25e-11_fp*exp(20.e+0_fp*TINV)
          ! 12:  N2O + hv -> N2 + O1D
          !RRATE(k_JN2O) = NOX_J(I,J,L,JN2OIDX)*DAYFRAC
          RRATE(k_JN2O) = ZPJ(LMINPHOT,RXN_N2O,I,J)

          ! Sanity check
          Where(RRate.lt.0.0e+0_fp) RRate = 0.0e+0_fp

          ! Retrieve local O3P/O1D mixing ratios and relevant
          ! j-rates from interpolated 2D arrays
          LOCALO3P = NOX_O(I,J,L,1)*NDAIR*DAYFRAC
          LOCALO1D = NOX_O(I,J,L,2)*NDAIR*DAYFRAC

          ! Partition NOx into N, NO, NO2 and NO3 based on PSSA
          ! Two cases: Daytime/nighttime
          NO_ALPHA = RRATE(k_JNO) / (RRATE(9)*0.21e+0_fp*NDAIR)
          NO_BETA  = (RRATE(k_JNO2)+(RRATE(1)*LOCALO3P)) / (RRATE(2)*LOCALO3)
          NO_GAMMA = (RRATE(3)*LOCALO3) / RRATE(k_JNO3)

          ! Calculate the partition fractions
          FRACNO2 = 1.e+0_fp/(1.e+0_fp+NO_GAMMA+(NO_BETA* &
                    (1.e+0_fp+NO_ALPHA)))
          FRACNO3 = NO_GAMMA * FRACNO2
          FRACNO  = NO_BETA  * FRACNO2
          FRACN   = NO_ALPHA * FRACNO

          ! Estimate net production rates for NOx and N2O
          NOXRATE = 2.e+0_fp*((RRATE(11)*LOCALN2O*LOCALO1D) - &
                    (((RRATE(7)*FRACN*FRACNO2)+(RRATE(8)*FRACN*FRACNO)) &
                    *LOCALNOX*LOCALNOX))
          N2ORATE = (RRATE(7)*FRACN*FRACNO2*LOCALNOX*LOCALNOX) - &
                    ((RRATE(k_JN2O)+((RRATE(10)+RRATE(11))*LOCALO1D))*LOCALN2O)

          ! Calculate NOx in kg NO at the beginning
          kgNOx = localNOx / (AIRMW / &
                  State_Chm%SpcData(id_NO )%Info%emMW_g*XAIR)

          ! Calculate total change in NOx and N2O
          ! Explicit Euler method (fast)
          DNOX = NOXRATE * DTCHEM / ( ( AIRMW / &
                 State_Chm%SpcData(id_NO)%Info%emMW_g) * XAIR )
          DN2O = N2ORATE * DTCHEM / ( AIRMW / &
                 State_Chm%SpcData(id_N2O)%Info%emMW_g * XAIR )

          ! Safety check - ensure NOx and N2O are positive
          IF ((DNOX*-1e+0_fp).gt.KGNOX) THEN
             DNOX = -1e+0_fp*KGNOX
             KGNOX = 0e+0_fp
             LOCALNOX = 0e+0_fp
          ELSE
             KGNOX = KGNOX + DNOX
             LOCALNOX = LOCALNOX + (NOXRATE*DTCHEM)
          ENDIF
          IF ((DN2O*-1e+0_fp).gt.KGN2O) THEN
             DN2O = -1e+0_fp*KGN2O
             KGN2O = 0e+0_fp
          ELSE
             KGN2O = KGN2O + DN2O
          ENDIF
       EndIf ! DayFrac < 0

       ! Convert to kg
       ! Now add all N to NO to get around problems caused by
       ! using negative molar mass for N (SDE 2018-03-19)
       NEW_N   = 0.0e+0_fp
       NEW_NO  = localNOx*(fracN+fracNO) / (AIRMW/ &
                 State_Chm%SpcData(id_NO )%Info%emMW_g * XAIR)
       NEW_NO2 = localNOx*fracNO2/(AIRMW/ &
                 State_Chm%SpcData(id_NO2)%Info%emMW_g * XAIR)
       NEW_NO3 = localNOx*fracNO3/(AIRMW/ &
                 State_Chm%SpcData(id_NO3)%Info%emMW_g * XAIR)
       NEW_N2O = kgN2O

       Spc(I,J,L,id_N)   = 0.d0
       Spc(I,J,L,id_NO)  = NEW_NO
       Spc(I,J,L,id_NO2) = NEW_NO2
       Spc(I,J,L,id_NO3) = NEW_NO3
       Spc(I,J,L,id_N2O) = NEW_N2O

       MESONOX_DELTA = MESONOX_DELTA + DNOX
       MESON2O_DELTA = MESON2O_DELTA + DN2O

    ENDDO ! J
    ENDDO ! I
    ENDDO ! L
    !$OMP END PARALLEL DO

    IF ( prtDebug ) THEN
       ! Print mean NOx tendency in mesosphere
       DBGMSG = ' ### UCX_NOX: Mesospheric NOx processed'
       CALL DEBUG_MSG(TRIM(DBGMSG))
       WRITE(DBGMSG,'(a,1(1x,F10.5),a)') ' ### Timestep:         ', &
            DTCHEM_HR, ' hours'
       CALL DEBUG_MSG(TRIM(DBGMSG))
       WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### NOx direct delta: ', &
            MESONOX_DELTA, ' kg'
       CALL DEBUG_MSG(TRIM(DBGMSG))
       WRITE(DBGMSG,'(a,1x,E10.4,1x,a)') ' ### N2O direct delta: ', &
            MESON2O_DELTA, ' kg'
       CALL DEBUG_MSG(TRIM(DBGMSG))
    ENDIF

    ! Free pointer
    NULLIFY( Spc )

  END SUBROUTINE UCX_NOX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_noxcoeff
!
! !DESCRIPTION: Subroutine GET\_NOXCOEFF reads in O1D and O3P mixing
!  ratios along with NO, NO2, NO3 and N2O J-rates from 2D data,
!  interpolating onto the 3D grid and storing in NOX\_O and NOX\_J.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_NOXCOEFF( TARG_MONTH, Input_Opt, State_Grid, State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE FILE_MOD,           ONLY : IOERROR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: TARG_MONTH
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  At some later point we should attempt to rewrite the parallel DO loop so
!  that the loop order is L-J-I.  Not sure how easy that is. (bmy, 2/14/14)
!
! !REVISION HISTORY:
!  26 Mar 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                               :: IOS
    INTEGER                               :: ILEV
    INTEGER                               :: I, J, L, ITRAC
    INTEGER                               :: JJNOX
    REAL(fp)                              :: PCENTER
    REAL(fp)                              :: CURRVAL
    INTEGER                               :: VERTCOUNT, AS
    LOGICAL                               :: FOUNDLEV,EXTRAP
    INTEGER                               :: IU_FILE, fId
    INTEGER                               :: st3d(3), ct3d(3) ! Start/count
    LOGICAL                               :: ISRATE
    CHARACTER(LEN=255)                    :: NOX_FILE
    CHARACTER(LEN=255)                    :: TARG_TRAC
    CHARACTER(LEN=255)                    :: DBGMSG
    REAL(fp), DIMENSION(:,:,:), POINTER   :: NOXDATA2D => NULL()
    REAL(fp), DIMENSION(:,:), ALLOCATABLE :: NOXD2D_IN
    INTEGER                               :: LSTART

    ! Local variables for quantities from Input_Opt
    LOGICAL                               :: prtDebug

    !=================================================================
    ! GET_NOXCOEFF begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Clear interpolated arrays
    NOX_O = 0e+0_fp
    NOX_J = 0e+0_fp

    IF (UCXNETCDF) THEN
       ! Allocate and zero 2D data array
       ALLOCATE( NOXDATA2D( State_Grid%NY, 51, 6 ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXDATA2D' )
       NOXDATA2D = 0e+0_fp

       ! Allocate and zero 2D input array
       ALLOCATE( NOXD2D_IN( UCX_NLAT, UCX_NLEVS ), STAT=AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXD2D_IN' )
       NOXD2D_IN = 0e+0_fp

       NOX_FILE = TRIM(NOON_FILE_ROOT)

       IF ( prtDebug ) THEN
          WRITE(DBGMSG,'(3a)') ' ### UCX: Reading mesospheric NOx coeff from ',&
               TRIM( NOX_FILE )
          CALL DEBUG_MSG( TRIM(DBGMSG) )
       ENDIF

       CALL NcOp_Rd (fId,TRIM(NOX_FILE))

       ! Start and count indices
       st3d = (/ 1,        1,          TARG_MONTH /)
       ct3d = (/ UCX_NLAT, UCX_NLEVS,  1          /)

       DO ITRAC = 1,6
          SELECT CASE (ITRAC)
          CASE ( 1 )
             TARG_TRAC = 'O'
          CASE ( 2 )
             TARG_TRAC = 'O1D'
          CASE ( 3 )
             TARG_TRAC = 'JNO'
          CASE ( 4 )
             TARG_TRAC = 'JNO2'
          CASE ( 5 )
             TARG_TRAC = 'JNO3'
          CASE ( 6 )
             TARG_TRAC = 'JN2O'
          END SELECT
          CALL NcRd( NOXD2D_IN, fId, TRIM(TARG_TRAC), st3d, ct3d )

          !! Debug
          !IF ( prtDebug ) THEN
          !   WRITE(DBGMSG,'(a,a,a)') ' ### UCX: Base ', &
          !      TRIM(TARG_TRAC), ', native: '
          !   CALL DEBUG_MSG( DBGMSG )
          !   DO J=1,19
          !      WRITE(DBGMSG,'(I02,x,E16.4)') J,NOXD2D_IN(J,1)
          !      CALL DEBUG_MSG( DBGMSG )
          !   ENDDO
          !ENDIF

          ! Regrid each level from 19 to State_Grid%NY latitudes
          ! Precalculated matrix for simple linear algebra
          ! Need to reverse levels - have layer 1 = TOA
          DO ILEV = 1,51
             NOXDATA2D(:,UCX_NLEVS+1-ILEV,ITRAC) = &
                  MATMUL(UCX_REGRID,NOXD2D_IN(:,ILEV))
          ENDDO

          !! Debug
          !IF ( prtDebug ) THEN
          !   WRITE(DBGMSG,'(a,a,a)') ' ### UCX: Base ', &
          !      TRIM(TARG_TRAC), ', regridded: '
          !   CALL DEBUG_MSG( DBGMSG )
          !   DO J=1,State_Grid%NY
          !      WRITE(DBGMSG,'(I02,x,E16.4)') J,NOXDATA2D(J,1,ITRAC)
          !      CALL DEBUG_MSG( DBGMSG )
          !   ENDDO
          !ENDIF
       ENDDO

       DEALLOCATE( NOXD2D_IN )
       CALL NcCl( fId )
    ELSE

       ! All the coefficients are now stored in NOXCOEFF.
       ! NOXDATA2D points to the desired month slice.
       ! (ckeller, 05/12/14)
       NOXDATA2D => NOXCOEFF(:,:,:,TARG_MONTH)

       !! Allocate and zero 2D data array
       !ALLOCATE( NOXDATA2D( State_Grid%NY, 51, 6 ), STAT=AS )
       !IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXDATA2D' )
       !NOXDATA2D = 0d0
       !DO ITRAC = 1,6
       !   SELECT CASE (ITRAC)
       !   CASE ( 1 )
       !      TARG_TRAC = 'O'
       !   CASE ( 2 )
       !      TARG_TRAC = 'O1D'
       !   CASE ( 3 )
       !      TARG_TRAC = 'JNO'
       !   CASE ( 4 )
       !      TARG_TRAC = 'JNO2'
       !   CASE ( 5 )
       !      TARG_TRAC = 'JNO3'
       !   CASE ( 6 )
       !      TARG_TRAC = 'JN2O'
       !   END SELECT
       !   WRITE(NOX_FILE,'(a,a,a,I0.2,a)') TRIM(NOON_FILE_ROOT), &
       !       TRIM(TARG_TRAC), '_', TARG_MONTH, '.dat'
       !
       !   ! Get a free LUN
       !   IU_FILE = findFreeLUN()
       !
       !   IOS = 1
       !   OPEN( IU_FILE,FILE=TRIM(NOX_FILE),STATUS='OLD',IOSTAT=IOS)
       !   IF ( IOS /= 0 ) THEN
       !      WRITE(6,*) 'UCX_MOD: Could not read ', TRIM(NOX_FILE)
       !      CALL IOERROR( IOS, IU_FILE,'UCX_MOD:GET_NOXCOEFF')
       !   ENDIF
       !
       !   IF ( prtDebug ) THEN
       !      WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', &
       !                TRIM( NOX_FILE )
       !      CALL DEBUG_MSG( TRIM(DBGMSG) )
       !   ENDIF
       !
       !   ! Read in data
       !   DO ILEV = 1,51
       !      READ(IU_FILE, 110, IOSTAT=IOS ) NOXDATA2D(:,ILEV,ITRAC)
       !      IF ( IOS /= 0 ) THEN
       !         WRITE(6,'(a,a,I4,a,1x,a)') 'UCX_MOD: Error reading ' &
       !            'line ', ILEV, ' in file ', TRIM( NOX_FILE )
       !         CALL IOERROR( IOS, IU_FILE,'UCX_MOD:GET_NOXCOEFF')
       !      ENDIF
       !   ENDDO
       !
       !   IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       !110   FORMAT(91E10.3)
       !   ELSE IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN
       !110   FORMAT(46E10.3)
       !   ELSE
       !      ! use 2x25 as default
       !110   FORMAT(91E10.3)
       !   ENDIF
       !
       !   CLOSE(IU_FILE)
       !
       !   ! Debug
       !   IF ( prtDebug ) THEN
       !      WRITE(DBGMSG,'(a,a,a)') ' ### UCX: Base ', &
       !         TRIM(TARG_TRAC), ', regridded: '
       !      CALL DEBUG_MSG( DBGMSG )
       !      DO J=1,State_Grid%NY
       !         WRITE(DBGMSG,'(I02,x,E16.4)') J,NOXDATA2D(J,1,ITRAC)
       !         CALL DEBUG_MSG( DBGMSG )
       !      ENDDO
       !   ENDIF
       !ENDDO

    ENDIF

    ! Scan through target array, element by element
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,      J,      L,       VERTCOUNT ) &
    !$OMP PRIVATE( EXTRAP, LSTART, PCENTER, FOUNDLEV  ) &
    !$OMP PRIVATE( ITRAC,  ISRATE, CURRVAL, JJNOX     )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Get corresponding J value in NOXDATA2D. This can be different
       ! from J if the simulation grid is not 2x25 or 4x5.
       JJNOX = GET_JJNOX( I, J, State_Grid )

       ! Vertcount is the layer count for the input, where layer 1
       ! is at the *top* of the atmosphere
       VERTCOUNT = 51
       EXTRAP = .TRUE.
       LSTART = State_Met%ChemGridLev(I,J)

       DO L = LSTART, State_Grid%NZ
          IF ( State_Met%InChemGrid(I,J,L) ) CYCLE

          ! Pressure at center of cell
          PCENTER = State_Met%PMID(I,J,L)
          FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
          DO WHILE (.not. FOUNDLEV)
             IF (VERTCOUNT.eq.1) THEN
                ! At top layer; use it anyway
                FOUNDLEV = .TRUE.
                EXTRAP = .TRUE.
             ELSE
                VERTCOUNT = VERTCOUNT - 1
                FOUNDLEV = (PCENTER.gt.UCX_PLEVS(VERTCOUNT))
                EXTRAP = .FALSE.
             ENDIF
          ENDDO

          DO ITRAC=1,6
             ISRATE = (ITRAC.gt.2)
             ! Interpolate data
             IF (EXTRAP) THEN
                ! Just take outside value if at edges
                CURRVAL = NOXDATA2D(JJNOX,VERTCOUNT,ITRAC)
             ELSE
                ! Interpolate by pressure
                CURRVAL = (UCX_PLEVS(VERTCOUNT+1)-PCENTER)
                CURRVAL = CURRVAL/(UCX_PLEVS(VERTCOUNT+1)-UCX_PLEVS(VERTCOUNT))
                CURRVAL = CURRVAL * &
                          (NOXDATA2D(JJNOX,VERTCOUNT+1,ITRAC)- &
                           NOXDATA2D(JJNOX,VERTCOUNT,ITRAC))
                CURRVAL = NOXDATA2D(JJNOX,VERTCOUNT+1,ITRAC) - CURRVAL
             ENDIF
             IF (.not.ISRATE) THEN
                ! Reading in mixing ratios (v/v)
                NOX_O(I,J,L,ITRAC) = CURRVAL
             ELSE
                ! J-rate (no conversion necessary)
                NOX_J(I,J,L,ITRAC-2) = CURRVAL
             ENDIF
          ENDDO ! ITRAC
       ENDDO ! L

    ENDDO ! I
    ENDDO ! J
    !$OMP END PARALLEL DO

    ! Cleanup
    IF (UCXNETCDF) THEN
       IF ( ASSOCIATED(NOXDATA2D) ) DEALLOCATE( NOXDATA2D )
    ELSE
       NOXDATA2D => NULL()
    ENDIF

  END SUBROUTINE GET_NOXCOEFF
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: settle_strat_aer
!
! !DESCRIPTION: Subroutine SETTLE\_STRAT\_AER performs gravitational settling
!  of stratospheric aerosols. It is copied largely from GRAV\_SETTLING in
!  sulfate\_mod.F90. All of this is ignored if APM is active.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SETTLE_STRAT_AER( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_FJX_MOD,        ONLY : RAA, IND999
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : IT_IS_NAN,ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC, GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    INTEGER,        INTENT(INOUT) :: RC          ! Return code
!
! !REVISION HISTORY:
!  11 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Common
    INTEGER                :: I,      J,     L,        DTCHEM
    REAL(fp)               :: DELZ,   DELZ1
    REAL(fp)               :: P,      TEMP,  DP,       PDP
    REAL(fp)               :: CONST,  SLIP,  VISC,     FAC1
    REAL(fp)               :: FAC2,   FLUX,  AREA_CM2, RHB
    REAL(fp)               :: RCM
    REAL(fp)               :: TOT1,   TOT2,  NIT_MW_G, HNO3_MW_G
    INTEGER                :: K
    LOGICAL                :: NATCOL

    ! Specific to each class
    REAL(fp)               :: RWET(2),CONST_V(2)
    REAL(fp)               :: RHO(2),RATIO_R(2),REFF(2)
    REAL(fp)               :: VTS(State_Grid%NZ,2)

    ! Used for old Seinfeld & Pandis slip factor calc
    REAL(fp)               :: sp_Lambda, sp_Num

    ! Parameters
    REAL(fp), PARAMETER    :: BCDEN = 1000.e+0_fp ! density (kg/m3)

    ! Indexing
    INTEGER, PARAMETER     :: IBC  = 1
    INTEGER, PARAMETER     :: ILIQ = 2
    INTEGER, PARAMETER     :: NSETTLE = 2
    INTEGER                :: IAERO
    LOGICAL                :: RUNCALC

    ! Partitioning
    REAL(fp)               :: PHASEMASS(3,2)
    REAL(fp)               :: SEDMASS
    INTEGER                :: IDTCURRENT

    ! NAT only
    REAL(fp)               :: VFALLMAX
    REAL(fp)               :: VNAT(State_Grid%NZ)
    INTEGER                :: MAXALT, MINALT
    REAL(fp)               :: BXMIN, SEDSTEP, TEMPREAL, P_0, P_ABOVE
    REAL(fp)               :: PSED0, PSEDABOVE, INVAIR_0, INVAIR_ABOVE
    REAL(fp)               :: XPSC_ABOVE, XPSC_0, XPSC_BELOW
    REAL(fp)               :: XNAT_ABOVE, XNAT_0, XNAT_BELOW
    REAL(fp)               :: XNO3_ABOVE, XNO3_0, XNO3_BELOW
    REAL(fp)               :: XICE_ABOVE, XICE_0, XICE_BELOW
    REAL(fp)               :: SEDPSC, SEDNAT, SEDICE, SEDQUANT
    REAL(fp)               :: SEDH2O, SEDNO3
    REAL(fp)               :: BELOWGRAD, ABOVEGRAD
    INTEGER                :: LOCALPROFILE, NUMSEDSTEPS, STARTPT, ISED

    ! Local variables for quantities from Input_Opt
    LOGICAL                :: LGRAVSTRAT

    ! Pointers
    REAL(fp), POINTER      :: Spc (:,:,:,:)
    REAL(f4), POINTER      :: STATE_PSC(:,:,:)
    REAL(fp), POINTER      :: WERADIUS(:,:,:,:)

    !=================================================================
    ! SETTLE_STRAT_AER begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Copy fields from INPUT_OPT
    LGRAVSTRAT  = Input_Opt%LGRAVSTRAT

    ! Copy fields from species database
    NIT_MW_G    = State_Chm%SpcData(id_NIT)%Info%emMW_g  ! g/mol
    HNO3_MW_G   = State_Chm%SpcData(id_HNO3)%Info%emMW_g ! g/mol

    ! Initialize pointers
    Spc       => State_Chm%Species     ! Chemical species [kg]
    STATE_PSC => State_Chm%STATE_PSC   ! PSC type (Kirner et al. 2011, GMD)
    WERADIUS  => State_Chm%WetAeroRadi ! Aerosol Radius [cm]

    ! Return if gravitational settling disabled
    IF (.not. LGRAVSTRAT) RETURN

    ! Return if it's the start of the run
    IF ( GET_ELAPSED_SEC() == 0 ) RETURN

    ! Chemistry timestep [s]
    DTCHEM = GET_TS_CHEM()

    ! First settle liquid aerosols (SLA) using scheme found
    ! elsewhere in GEOS-Chem
    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( J,            I,            L,          VTS        ) &
    !$OMP PRIVATE( VNAT,         NATCOL,       BXMIN,      MINALT     ) &
    !$OMP PRIVATE( MAXALT,       VFALLMAX,     TEMP,       P          ) &
    !$OMP PRIVATE( RUNCALC,      RWET,         RHO                    ) &
    !$OMP PRIVATE( SP_NUM,       SP_LAMBDA,    VISC,       IAERO      ) &
    !$OMP PRIVATE( DP,           PDP,          CONST,      SLIP       ) &
    !$OMP PRIVATE( DELZ,         CONST_V,      PHASEMASS,  K          ) &
    !$OMP PRIVATE( IDTCURRENT,   SEDMASS,      DELZ1,      SEDSTEP    ) &
    !$OMP PRIVATE( NUMSEDSTEPS,  TEMPREAL,     ISED,       STARTPT    ) &
    !$OMP PRIVATE( XNO3_0,       XNAT_0,       XICE_0,     XPSC_0     ) &
    !$OMP PRIVATE( INVAIR_0,     INVAIR_ABOVE, XNO3_ABOVE, XNAT_ABOVE ) &
    !$OMP PRIVATE( XICE_ABOVE,   XPSC_ABOVE,   P_ABOVE,    ABOVEGRAD  ) &
    !$OMP PRIVATE( XNO3_BELOW,   XNAT_BELOW,   XICE_BELOW, XPSC_BELOW ) &
    !$OMP PRIVATE( P_0,          PSED0,        PSEDABOVE,  BELOWGRAD  ) &
    !$OMP PRIVATE( LOCALPROFILE, SEDQUANT,     SEDNO3,     SEDICE     ) &
    !$OMP PRIVATE( SEDH2O,       SEDNAT,       SEDPSC                 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize
       DO L = 1, State_Grid%NZ
          VTS(L,:) = 0e+0_fp
          VNAT(L) = 0e+0_fp
       ENDDO
       NATCOL = .FALSE.

       ! Arbitrary limits
       BXMIN = 1.e+20_fp
       MINALT = State_Grid%NZ+1
       MAXALT = 0
       VFALLMAX = 0e+0_fp

       ! Loop over levels
       DO L = 1, State_Grid%NZ

          ! Temperature [K]
          TEMP    = State_Met%T(I,J,L)

          ! Pressure at center of the level [kPa]
          ! Use moist air pressure for mean free path (ewl, 3/2/15)
          P       = State_Met%PMID(I,J,L)  * 0.1e+0_fp

          RUNCALC = State_Met%InStratMeso(I,J,L)

          IF (RUNCALC) THEN
             ! Need to translate for BC radii
             IF ( State_Met%InChemGrid(I,J,L) ) THEN
                RWET(IBC) = WERADIUS(I,J,L,2)*1.e-2_fp
             ELSE
                ! Use defaults, assume dry (!)
                RWET(IBC) = RAA(IND999,29) * 1.0e-6_fp
             ENDIF

             ! Taken from aerosol_mod (MSDENS(2))
             RHO(IBC) = BCDEN

             ! Get aerosol properties
             RWET(ILIQ) = RAD_AER(I,J,L,I_SLA)*1.e-2_fp
             RHO(ILIQ)  = RHO_AER(I,J,L,I_SLA)

             ! Do we need to sediment NAT?
             IF (NDENS_AER(I,J,L,I_SPA).gt.TINY(0e+0_fp)) THEN
                NATCOL = .TRUE.
                BXMIN  = MIN(BXMIN, State_Met%BXHEIGHT(I,J,L))
                MINALT = MIN(L,MINALT)
                MAXALT = MAX(L,MAXALT)
             ENDIF
          ENDIF

          IF (.not.RUNCALC) THEN
             VTS(L,:) = 0e+0_fp
             VNAT(L) = 0e+0_fp
          ELSE
             ! Calculate common variables first
             sp_Num = P * 1.e+3_fp * AVO / (RSTARG * Temp)
             sp_Lambda = 1.e+6_fp / ( SQRT(2e+0_fp) * sp_Num * PI &
                         * (3.7e-10_fp)**2 )

             ! Viscosity [Pa*s] of air as a function of temperature
             VISC = 1.458e-6_fp * (Temp)**(1.5e+0_fp) &
                    / ( Temp + 110.4e+0_fp )

             DO IAERO=1,2
                IF (RWET(IAERO).le.TINY(0e+0_fp)) THEN
                   VTS(L,IAERO) = 0e+0_fp
                ELSE
                   ! Dp = particle diameter [um]
                   DP = 2.e+0_fp * RWET(IAERO) * 1.e+6_fp

                   ! PdP = P * dP [hPa * um]
                   PDp     = P * Dp

                   ! Constant
                   CONST = 2.e+0_fp * RHO(IAERO) * RWET(IAERO)**2 &
                           * g0 / 9.e+0_fp

                   !=========================================================
                   ! NOTE: Slip correction factor calculations following
                   ! Seinfeld, pp464 which is thought to be more accurate
                   ! but more computation required. (rjp, 1/24/02)
                   !
                   ! # air molecule number density
                   ! num = P * 1d3 * 6.023d23 / (8.314 * Temp)
                   !
                   ! # gas mean free path
                   ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
                   !
                   ! # Slip correction
                   ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp
                   !     &     / (2. * lamda))) / Dp
                   !
                   ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
                   ! which produces slip correction factore with small error
                   ! compared to the above with less computation.
                   !========================================================= 

                   ! Slip correction factor (as function of P*dp)
                   ! Slip = 1.e+0_fp+(15.60e+0_fp + 7.0e+0_fp * &
                   !           EXP(-0.059e+0_fp * PDp)) / PDp
                   ! Reverting to Seinfeld and Pandis
                   Slip = 1. + 2. * sp_Lambda * (1.257 + 0.4*exp(-1.1* &
                          Dp / (2. * sp_Lambda))) / Dp

                   ! Settling velocity [m/s]
                   VTS(L,IAERO) = CONST * Slip / VISC
                ENDIF ! RWET
             ENDDO ! IAERO

             ! Now solid PSC particles
             IF (NATCOL) THEN
                ! sp_Num: Air molecule#/m3
                VNAT(L) = CALC_FALLVEL(RHO_AER(I,J,L,I_SPA), &
                          RAD_AER(I,J,L,I_SPA),TEMP,P)
                IF (VNAT(L).gt.VFALLMAX) VFALLMAX = VNAT(L)
             ELSE
                VNAT(L) = 0.e+0_fp
             ENDIF
          ENDIF ! RUNCALC
       ENDDO

       ! First apply simpler SLA sedimentation scheme
       ! Handle model top condition
       L    = State_Grid%NZ
       DELZ = State_Met%BXHEIGHT(I,J,L)

       DO IAERO=1,2
          CONST_V(IAERO) = 1.e+0_fp / (1.e+0_fp + DTCHEM * VTS(L,IAERO) / DELZ)
       ENDDO

       ! Zero arrays
       PHASEMASS(:,:) = 0e+0_fp

       ! Only want to sediment fraction of species currently
       ! in the aerosol
       DO K = 1,7
          ! Only process transported species
          IDTCURRENT = AERFRACIND(K)
          IF (IDTCURRENT.ne.0) THEN
             ! Calculate local phase partitioning
             ! Total upper gridbox mass
             PHASEMASS(3,2) = Spc(I,J,L,IDTCURRENT)
             ! Aerosol-phase upper gridbox mass
             PHASEMASS(2,2) = AERFRAC(I,J,L,K)*PHASEMASS(3,2)
             ! Gas-phase upper gridbox mass
             PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)

             ! Calculate total sedimented mass
             SEDMASS = PHASEMASS(2,2) * (1.e+0_fp-CONST_V(ILIQ))

             ! Remove from upper gridbox
             PHASEMASS(2,2) = PHASEMASS(2,2) - SEDMASS
             PHASEMASS(3,2) = PHASEMASS(1,2) + PHASEMASS(2,2)

             ! Recalculate phase fractions
             IF (PHASEMASS(3,2).gt.TINY(1e+0_fp)) THEN
                AERFRAC(I,J,L,K) = PHASEMASS(2,2)/PHASEMASS(3,2)
             ELSE
                AERFRAC(I,J,L,K) = 0e+0_fp
             ENDIF

             ! Store result
             Spc(I,J,L,IDTCURRENT) = PHASEMASS(3,2)
          ENDIF
       ENDDO
       Spc(I,J,L,id_BCPI)= Spc(I,J,L,id_BCPI)* CONST_V(IBC)

       DO L = State_Grid%NZ-1,1,-1
          IF ( State_Met%InTroposphere(I,J,L+1) ) CYCLE
          DELZ  = State_Met%BXHEIGHT(I,J,L)
          DELZ1 = State_Met%BXHEIGHT(I,J,L+1)

          DO K=1,7
             IDTCURRENT = AERFRACIND(K)
             IF (IDTCURRENT.ne.0) THEN
                ! Total upper gridbox mass
                PHASEMASS(3,2) = Spc(I,J,L+1,IDTCURRENT)
                ! Aerosol-phase upper gridbox mass
                PHASEMASS(2,2) = AERFRAC(I,J,L+1,K)*PHASEMASS(3,2)
                ! Gas-phase upper gridbox mass
                PHASEMASS(1,2) = PHASEMASS(3,2) - PHASEMASS(2,2)

                ! Total lower gridbox mass
                PHASEMASS(3,1) = Spc(I,J,L,IDTCURRENT)
                ! Aerosol-phase lower gridbox mass
                PHASEMASS(2,1) = AERFRAC(I,J,L,K)*PHASEMASS(3,1)
                ! Gas-phase lower gridbox mass
                PHASEMASS(1,1) = PHASEMASS(3,1) - PHASEMASS(2,1)

                ! New lower gridbox mass
                PHASEMASS(2,1) = 1.e+0_fp/(1.e+0_fp+DTCHEM &
                                 * VTS(L,ILIQ) / DELZ) &
                                 * (PHASEMASS(2,1)+DTCHEM*VTS(L+1,ILIQ)/DELZ1 &
                                 * PHASEMASS(2,2))

                ! Calculate new total mass in lower gridbox
                PHASEMASS(3,1) = PHASEMASS(2,1) + PHASEMASS(1,1)

                ! Recalculate phase fraction
                IF (PHASEMASS(3,1).gt.TINY(1e+0_fp)) THEN
                   AERFRAC(I,J,L,K) = PHASEMASS(2,1)/PHASEMASS(3,1)
                ELSE
                   AERFRAC(I,J,L,K) = 0e+0_fp
                ENDIF

                ! Store result
                Spc(I,J,L,IDTCURRENT) = PHASEMASS(3,1)
             ENDIF
          ENDDO
          Spc(I,J,L,id_BCPI) = 1.e+0_fp/(1.e+0_fp+DTCHEM &
                               * VTS(L,IBC) / DELZ) &
                               * (Spc(I,J,L,id_BCPI)+DTCHEM*VTS(L+1,IBC)/DELZ1 &
                               * Spc(I,J,L+1,id_BCPI))
       ENDDO

       ! Now perform trapezoidal scheme for particulates
       ! Calculate maximum allowable timestep (seconds)
       IF (VFALLMAX.gt.TINY(1e+0_fp)) THEN
          SEDSTEP = BXMIN/VFALLMAX
          IF(DTCHEM.le.SEDSTEP)THEN
             NUMSEDSTEPS = 1
             SEDSTEP = DTCHEM
          ELSE
             ! Will need to run iteratively
             ! Calculate minimum necessary number of steps, limiting to
             ! 10 steps if excessive
             TEMPREAL = DTCHEM/SEDSTEP
             NUMSEDSTEPS = CEILING(TEMPREAL)
             NUMSEDSTEPS = MIN(10,NUMSEDSTEPS)
             SEDSTEP = DTCHEM/(NUMSEDSTEPS*1.e+0_fp)
             VFALLMAX = BXMIN/SEDSTEP
          ENDIF

          SEDSTEPLOOP: DO ISED=1,NUMSEDSTEPS
             STARTPT = MAX(1,MINALT-1)
             ! XPSC is the number of molecules tied up in solid particles
             ! per m3 in a grid box
             L = STARTPT
             IF (STATE_PSC(I,J,L) < 2.0e+0_f4 ) THEN
                XNO3_0 = 0e+0_fp
                XNAT_0 = 0e+0_fp
                XICE_0 = 0e+0_fp
             ELSE
                XNO3_0 = Spc(I,J,L,id_NIT) * AIRMW / &
                         ( NIT_MW_G * State_Met%AD(I,J,L) )
                XNAT_0 = XNO3_0! * 4.e+0_fp
                XICE_0 = (KG_AER(I,J,L,I_SPA)- &
                         ( ( NATMW / NIT_MW_G ) * &
                         Spc(I,J,L,id_NIT)))*AIRMW/(ICEMW* &
                         State_Met%AD(I,J,L))
             ENDIF
             XPSC_0   = XNAT_0 + XICE_0
             P_0      = 100.0e+0_fp * State_Met%PMID(I,J,L)
             INVAIR_0 = AIRMW/State_Met%AD(I,J,L)
             IF (L .lt. (State_Grid%NZ-1)) THEN
                INVAIR_ABOVE  = AIRMW/State_Met%AD(I,J,L+1)
                IF (STATE_PSC(I,J,L+1) < 2.0e+0_f4 ) THEN
                   XNO3_ABOVE = 0e+0_fp
                   XNAT_ABOVE = 0e+0_fp
                   XICE_ABOVE = 0e+0_fp
                ELSE
                   XNO3_ABOVE = Spc(I,J,L+1,id_NIT) * INVAIR_ABOVE / NIT_MW_G
                   ! NAT = HNO3.3H2O = 4 molecules
                   XNAT_ABOVE = XNO3_ABOVE! * 4.e+0_fp
                   ! XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)-Spc(I,J,L+1,id_NIT)) &
                   !                            *INVAIR_ABOVE/ICEMW
                   XICE_ABOVE = (KG_AER(I,J,L+1,I_SPA)- &
                                ( ( NATMW / NIT_MW_G ) * &
                                Spc(I,J,L+1,id_NIT)))*INVAIR_ABOVE/ICEMW
                ENDIF
                XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                P_ABOVE    = 100.0e+0_fp * State_Met%PMID(I,J,L+1)

                ! Replace PEDGE(I,J,L+1) - PEDGE(I,J,L+2)
                ! with equivalent State_Met%DELP(I,J,L+1)
                ! (ewl, 3/2/15)
                PSEDABOVE  = VNAT(L+1) * SEDSTEP * 1.e+2_fp * &
                             State_Met%DELP(I,J,L+1) &
                             / State_Met%BXHEIGHT(I,J,L+1)
                !PSEDABOVE=g0 * (1.0d3*AIRMW) * P_ABOVE * &
                !   VNAT(L+1) / (Rd * T(I,J,L+1))*SEDSTEP
             ELSE
                INVAIR_ABOVE = 1.0e+0_fp
                XNO3_ABOVE = 0.0e+0_fp
                XNAT_ABOVE = 0.0e+0_fp
                XICE_ABOVE = 0.0e+0_fp
                XPSC_ABOVE = 0.0e+0_fp
                P_ABOVE = State_Met%PEDGE(I,J,L+2) * 100.0e+0_fp
                PSEDABOVE = 0.0e+0_fp
             ENDIF
             ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)
             SED_LLOOP: DO L=STARTPT,MAXALT
                ! Actually calculating sedimentation for the box above this
                ! one. By the time we get around to processing a given
                ! box, it is actually the (i-1)th box
                ! We are therefore concerned with calculating
                ! sedimentation *out* of the box above us and *into* this
                ! one
                XNO3_BELOW = XNO3_0
                XNAT_BELOW = XNAT_0
                XICE_BELOW = XICE_0
                XPSC_BELOW = XPSC_0
                XNO3_0 = XNO3_ABOVE
                XNAT_0 = XNAT_ABOVE
                XICE_0 = XICE_ABOVE
                XPSC_0 = XPSC_ABOVE
                P_0 = P_ABOVE
                PSED0 = PSEDABOVE
                INVAIR_0 = INVAIR_ABOVE

                IF (L.lt.(State_Grid%NZ-1)) THEN
                   INVAIR_ABOVE  = AIRMW/State_Met%AD(I,J,L+2)
                   IF (STATE_PSC(I,J,L+2) >= 2.0e+0_f4 ) THEN
                      XNO3_ABOVE = Spc(I,J,L+2,id_NIT)*INVAIR_ABOVE/NIT_MW_G
                      XNAT_ABOVE = XNO3_ABOVE! * 4.e+0_fp
                      XICE_ABOVE = (KG_AER(I,J,L+2,I_SPA)- &
                                   ( ( NATMW / NIT_MW_G )* &
                                   Spc(I,J,L+2,id_NIT)))*INVAIR_ABOVE/ICEMW
                   ELSE
                      XNO3_ABOVE = 0e+0_fp
                      XNAT_ABOVE = 0e+0_fp
                      XICE_ABOVE = 0e+0_fp
                   ENDIF
                   XPSC_ABOVE = XNAT_ABOVE + XICE_ABOVE
                   P_ABOVE    = State_Met%PMID(I,J,L+2) * 100.0e+0_fp

                   ! Replace PEDGE(I,J,L+2) - PEDGE(I,J,L+3)
                   ! with equivalent State_Met%DELP(I,J,L+2)
                   ! (ewl, 3/2/15)
                   PSEDABOVE  = VNAT(L+2) * SEDSTEP * 1.e+2_fp * &
                                State_Met%DELP(I,J,L+2) &
                                / State_Met%BXHEIGHT(I,J,L+2)
                   !PSEDABOVE = g0 * (1.0d3*AIRMW) * P_ABOVE * &
                   !   VNAT(L+2)/(Rd * T(I,J,L+2))*SEDSTEP
                ELSE
                   INVAIR_ABOVE = 0.0e+0_fp
                   XNO3_ABOVE = 0.0e+0_fp
                   XNAT_ABOVE = 0.0e+0_fp
                   XICE_ABOVE = 0.0e+0_fp
                   XPSC_ABOVE = 0.0e+0_fp
                   P_ABOVE = State_Met%PEDGE(I,J,L+2) * 100.0e+0_fp
                   PSEDABOVE = 0.0e+0_fp
                ENDIF
                ! Note reversal of pressure values as pressure falls with
                ! height. If a value is positive, PSC particle substance
                ! mixing ratios are increasing with height, and vice versa
                BELOWGRAD = ABOVEGRAD
                ABOVEGRAD = (XPSC_ABOVE-XPSC_0)/(P_0-P_ABOVE)

                ! NB: Order of following conditionals is important! Could
                ! theoretically do without localprofile as a variable, but
                ! should be sure to check for the possible (but phenomenally
                ! unlikely) condition that the gradients are both zero
                IF ((ABOVEGRAD .eq. 0e+0_fp) .and. &
                    (BELOWGRAD .eq. 0e+0_fp)) THEN
                   ! Extremelely unlikely case, but use standard scheme
                   ! for flat vertical profile
                   LOCALPROFILE = 0
                ELSEIF ((ABOVEGRAD .le. 0e+0_fp) .and. &
                        (BELOWGRAD .le. 0e+0_fp)) THEN
                   ! Consistently decreasing with height - above a peak
                   LOCALPROFILE = -1
                ELSEIF ((ABOVEGRAD .ge. 0e+0_fp) .and. &
                        (BELOWGRAD .ge. 0e+0_fp)) THEN
                   ! Consistently increasing with height - below a peak
                   LOCALPROFILE = +1
                ELSE
                   ! Local minmum or maximum
                   LOCALPROFILE = 0
                ENDIF

                IF (LOCALPROFILE.eq.0) THEN
                   ! Standard sedimentation scheme
                   SEDQUANT = XPSC_0 * PSED0
                   SEDNO3   = XNO3_0 * PSED0
                   SEDICE   = XICE_0 * PSED0
                ELSEIF (LOCALPROFILE.eq.1) THEN
                   ! Currently below a peak
                   ! Use gradient of bottom and centre boxes
                   IF (XPSC_0 .le. (0.5e+0_fp * XPSC_ABOVE)) THEN
                      ! Interpret as nearing peak; use lower gradient
                      SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5e+0_fp*PSED0
                      SEDNO3   = (XNO3_BELOW+XNO3_0)*0.50*PSED0
                      SEDICE   = (XICE_BELOW+XICE_0)*0.50*PSED0
                   ELSE
                      SEDQUANT = (XPSC_ABOVE+XPSC_0)*0.50*PSED0
                      SEDNO3   = (XNO3_ABOVE+XNO3_0)*0.50*PSED0
                      SEDICE   = (XICE_ABOVE+XICE_0)*0.50*PSED0
                   ENDIF
                ELSE
                   ! Above a peak
                   SEDQUANT = (XPSC_0+XPSC_BELOW)*0.5e+0_fp*PSED0
                   SEDNO3   = (XNO3_0+XNO3_BELOW)*0.5e+0_fp*PSED0
                   SEDICE   = (XICE_0+XICE_BELOW)*0.5e+0_fp*PSED0
                ENDIF
                ! Divide sedimenting quantity by the pressure difference
                ! across the box being sedimented from
                ! Note conversion from hPa to Pa for denominator, ie
                ! multiply numerator by (1/100)
                IF (L .ne. State_Grid%NZ) THEN

                   ! Replace PEDGE(I,J,L+1) - PEDGE(I,J,L+2)
                   ! with equivalent State_Met%DELP(I,J,L+1)
                   ! (ewl, 3/2/15)
                   SEDQUANT = (1.e-2_fp) * SEDQUANT / State_Met%DELP(I,J,L+1)
                   SEDNO3   = (1.e-2_fp) * SEDNO3   / State_Met%DELP(I,J,L+1)
                   SEDICE   = (1.e-2_fp) * SEDICE   / State_Met%DELP(I,J,L+1)
                ELSE
                   ! This shouldn't be possible?
                   CALL ERROR_STOP('Unknown sedimentation error', 'UCX_mod.F90')
                ENDIF

                ! Apply limits so that sedimented quantity is:
                ! Greater than or equal to zero
                ! Less than or equal to the total available quantity

                ! Note that we are calculating sedimentation using
                ! the total ice and NAT, but are only actually
                ! transporting the local NO3
                SEDQUANT = MAX(0.0e+0_fp,SEDQUANT)
                SEDNO3 = MAX(0.0e+0_fp,SEDNO3)
                SEDICE = MAX(0.0e+0_fp,SEDICE)

                ! Convert v/v to kg/box
                SEDNO3 = SEDNO3 * NIT_MW_G / INVAIR_0
                SEDICE = SEDICE*ICEMW/INVAIR_0

                SEDNO3 = MIN(SEDNO3,Spc(I,J,L+1,id_NIT))
                Spc(I,J,L,id_NIT) = Spc(I,J,L,id_NIT) + SEDNO3
                Spc(I,J,L+1,id_NIT)=Spc(I,J,L+1,id_NIT)-SEDNO3

                ! Settle the ice out too
                SEDH2O = MIN(SEDICE,Spc(I,J,L+1,id_H2O))
                Spc(I,J,L,id_H2O) = Spc(I,J,L,id_H2O) + SEDH2O
                Spc(I,J,L+1,id_H2O)=Spc(I,J,L+1,id_H2O)-SEDH2O

                ! Now correct aerosol totals
                SEDNAT = SEDNO3 * NATMW / NIT_MW_G
                SEDPSC = SEDNAT + SEDICE
                SEDPSC = MIN(SEDPSC,KG_AER(I,J,L+1,I_SPA))
                KG_AER(I,J,L,I_SPA) = KG_AER(I,J,L,I_SPA) + SEDPSC
                KG_AER(I,J,L+1,I_SPA)=KG_AER(I,J,L+1,I_SPA)-SEDPSC

             ENDDO SED_LLOOP
          ENDDO SEDSTEPLOOP
       ENDIF ! VFALLMAX > 0

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    Spc       => NULL()
    STATE_PSC => NULL()
    WERADIUS  => NULL()

  END SUBROUTINE SETTLE_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_h2so4_gas
!
! !DESCRIPTION: Subroutine CALC\_H2SO4\_GAS calculates the fraction of strat.
!  SO4 aerosol which can be considered to be gaseous H2SO4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_H2SO4_GAS( Input_Opt, State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  11 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRST=.TRUE.
    REAL(fp),PARAMETER :: GF_THRESHOLD = 0.0e+0_fp
    REAL(fp),PARAMETER :: GF_RANGE     = 1.0e-8_fp
    REAL(fp),PARAMETER :: GF_DELTAHBYR = 10156.e+0_fp
    REAL(fp),PARAMETER :: GF_T0        = 360.e+0_fp
    REAL(fp),PARAMETER :: GF_TC        = 905.e+0_fp
    REAL(fp),SAVE      :: GF_LOGP0
    REAL(fp),SAVE      :: GF_BFACTOR
    REAL(fp),SAVE      :: GF_ATMCONV
    REAL(fp),SAVE      :: GF_INVT0
    REAL(fp)           :: GF_INVT,GF_LOGPSULFATE,GF_CFACTOR
    REAL(fp)           :: GF_AFACTOR
    REAL(fp)           :: GF_PP,GF_PVAP,GF_DIFF

    INTEGER            :: I, J, L
    REAL(fp)           :: PCENTER, PCENTER_P, TCENTER, H2SO4SUM
    REAL(fp)           :: INVAIR
    REAL(fp)           :: SO4_MW_G

    ! Local variables for quantities from Input_Opt
    LOGICAL            :: prtDebug

    ! Pointers
    REAL(fp), POINTER  :: Spc (:,:,:,:)

    !=================================================================
    ! CALC_H2SO4_GAS begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    prtDebug = Input_Opt%LPRT .and. Input_Opt%amIRoot

    ! Copy fields from species database
    SO4_MW_G = State_Chm%SpcData(id_SO4)%Info%emMW_g ! g/mol

    ! Initialize GEOS-Chem species array [kg]
    Spc => State_Chm%Species

    IF (FIRST) THEN
       FIRST = .FALSE.
       ! Calculate H2SO4 gas phase prefactors
       GF_INVT0 = 1.e+0_fp/GF_T0
       GF_LOGP0 = (-1.e+0_fp*GF_DELTAHBYR*GF_INVT0) + 16.259e+0_fp
       GF_BFACTOR = 0.38e+0_fp/(GF_TC - GF_T0)
       GF_ATMCONV = LOG(ATM)
    ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( J,          I,              L       ) &
    !$OMP PRIVATE( PCENTER,    PCENTER_P,      TCENTER ) &
    !$OMP PRIVATE( INVAIR )                              &
    !$OMP PRIVATE( H2SO4SUM,   GF_PP,          GF_INVT ) &
    !$OMP PRIVATE( GF_CFACTOR, GF_LOGPSULFATE, GF_PVAP ) &
    !$OMP PRIVATE( GF_DIFF                             ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       ! Only interested in low-pressure boxes
       TCENTER = State_Met%T(I,J,L)

       ! Use moist air pressure thresholding (ewl, 3/26/15)
       PCENTER = State_Met%PMID(I,J,L)
       ! Use dry air partial pressure for part. P calc (ewl, 3/26/15)
       PCENTER_P = State_Met%PMID_DRY(I,J,L)

       INVAIR  = AIRMW / State_Met%AD(I,J,L)
       IF (PCENTER.ge.1.e+2_fp) THEN
          AERFRAC(I,J,L,1) = 1e+0_fp
       ELSEIF ( State_Met%InTroposphere(I,J,L) ) THEN
          ! Don't want to interfere with tropospheric aerosols
          AERFRAC(I,J,L,1) = 1e+0_fp
       ELSE
          H2SO4SUM = Spc(I,J,L,id_SO4) * INVAIR / SO4_MW_G
          ! Use approximation from Kumala (1990)
          ! Use dry air partial pressure (ewl, 3/26/15)
          GF_PP = H2SO4SUM*PCENTER_P
          GF_INVT = 1./TCENTER
          GF_CFACTOR = 1.e+0_fp+(LOG(GF_T0*GF_INVT))-(GF_T0*GF_INVT)
          GF_LOGPSULFATE = GF_LOGP0 + (GF_DELTAHBYR*(GF_INVT0 - &
                           GF_INVT + (GF_BFACTOR*GF_CFACTOR)))
          GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
          GF_PVAP = 1.e-2_fp * EXP(GF_LOGPSULFATE)
          GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
          IF (GF_DIFF .lt. 0) THEN
             AERFRAC(I,J,L,1) = 1.e+0_fp
          ELSEIF (GF_DIFF .lt. GF_RANGE) THEN
             AERFRAC(I,J,L,1) = 1.e+0_fp-(GF_DIFF/GF_RANGE)
          ELSE
             AERFRAC(I,J,L,1) = 0e+0_fp
          ENDIF
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    NULLIFY( Spc )

    IF ( prtDebug ) CALL DEBUG_MSG( '### UCX: H2SO4 partitioned' )

  END SUBROUTINE CALC_H2SO4_GAS
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: so4_photfrac
!
! !DESCRIPTION: FUNCTION SO4\_PHOTFRAC returns the fraction of H2SO4 which
!  is available for photolysis.
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION SO4_PHOTFRAC(I,J,L)
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: I,J,L      ! Location indices
!
! !OUTPUT VARIABLES:
!
!   REAL(fp), INTENT(OUT) :: PHOTFRAC   ! Gaseous fraction of H2SO4
!
! !REVISION HISTORY:
!  11 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! SO4_PHOTFRAC begins here!
    !=================================================================

    SO4_PHOTFRAC = 1.e+0_fp - AERFRAC(I,J,L,1)

  END FUNCTION SO4_PHOTFRAC
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_strat_aer
!
! !DESCRIPTION: Subroutine CALC\_STRAT\_AER calculates aerosol properties
!  stratosphere using the thermodynamic parameterization described in
!  Kirner et al. (`Simulation of polar stratospheric clouds in the
!  chemistry-climate-model EMAC via the submodel PSC', Geosci. Mod. Dev.,
!  4, 169-182, 2011).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_STRAT_AER( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : IS_SAFE_DIV
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  13 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Limits on PSC formation
    REAL(fp),PARAMETER      :: PSC_MAXLAT =  45.0e+0_fp
    REAL(fp),PARAMETER      :: PSC_MINLAT = -55.0e+0_fp
    REAL(fp),PARAMETER      :: PSC_PMAX   =  18.0e+3_fp
    REAL(fp),PARAMETER      :: PSC_PMIN   =   5.0e+2_fp

    ! Allow PSC formation outsize Kirner limits?
    LOGICAL,PARAMETER       :: PSC_FULL   =.FALSE.

    ! Saturation and nucleation properties of HNO3
    REAL(fp), PARAMETER     :: TSATHNO3_A = -2.7836e+0_fp
    REAL(fp), PARAMETER     :: TSATHNO3_B = -0.00088e+0_fp
    REAL(fp), PARAMETER     :: TSATHNO3_C = 38.9855e+0_fp
    REAL(fp), PARAMETER     :: TSATHNO3_D = -11397.0e+0_fp
    REAL(fp), PARAMETER     :: TSATHNO3_E = 0.009179e+0_fp

    ! Saturation and nucleation properties of water
    REAL(fp), PARAMETER     :: TSATH2O_A  = -2663.5e+0_fp
    REAL(fp), PARAMETER     :: TSATH2O_B  = 12.537e+0_fp

    ! Peak pressure at which NAT can form homogeneously
    REAL(fp), PARAMETER     :: P_MAXNAT = 1.40e+4_fp ! Pa

    ! Maximum temperature for PSC formation (K)
    REAL(fp), PARAMETER     :: T_MAX = 215.0e+0_fp

    ! Limits on NAT/ice formation
    REAL(fp), PARAMETER     :: MIN_RAD = 1.0e-7_fp ! m
    REAL(fp), PARAMETER     :: MAX_NDENS=42.0e+3_fp ! #/m3

    ! Local conditions
    REAL(fp)                :: TCENTER, PCENTER, PCENTER_PA, DENAIR
    REAL(fp)                :: INVAIR, TINV, TOFFSET

    ! Gridbox mixing ratios and partial pressures
    REAL(fp)                :: HNO3SUM, H2OSUM
    REAL(fp)                :: HNO3PP,  H2OPP
    REAL(fp)                :: PSATHNO3_SUPERCOOL
    REAL(fp)                :: PSATH2O_SUPERSAT
    REAL(fp)                :: PSATHNO3, PSATH2O
    REAL(fp)                :: H2SO4SUM
    REAL(fp)                :: ClNO3SUM, HClSUM, HOClSUM
    REAL(fp)                :: BrNO3SUM, HBrSUM, HOBrSUM

    ! Gridbox aerosol and phase data
    REAL(fp)                :: HNO3_BOX_G, HNO3_BOX_L, HNO3_BOX_S
    REAL(fp)                :: H2O_BOX_G,  H2O_BOX_L,  H2O_BOX_S
    REAL(fp)                :: H2SO4_BOX_G,H2SO4_BOX_L
    REAL(fp)                :: HCl_BOX_G,  HCl_BOX_L
    REAL(fp)                :: HOCl_BOX_G, HOCl_BOX_L
    REAL(fp)                :: HBr_BOX_G,  HBr_BOX_L
    REAL(fp)                :: HOBr_BOX_G, HOBr_BOX_L
    REAL(fp)                :: HNO3GASFRAC, HClGASFRAC, HOClGASFRAC
    REAL(fp)                :: HBrGASFRAC, HOBrGASFRAC
    REAL(fp)                :: VOL_NAT, VOL_ICE, VOL_SLA, VOL_TOT
    REAL(fp)                :: RAD_AER_BOX,RHO_AER_BOX
    REAL(fp)                :: KG_AER_BOX,NDENS_AER_BOX,SAD_AER_BOX
    REAL(fp)                :: KG_NAT, KG_ICE, KG_NO3

    ! SLA weight fractions
    REAL(fp)                :: W_H2SO4, W_H2O, W_HNO3
    REAL(fp)                :: W_HCl, W_HOCl, W_HBr, W_HOBr

    ! Reaction prefactors
    REAL(fp)                :: KHET_COMMON
    REAL(fp)                :: KHET_SPECIFIC

    ! Grid box location
    REAL(fp)                :: BOX_LAT_S, BOX_LAT_N, BOX_LAT
    LOGICAL                 :: IS_VALID, IS_POLAR, IS_STRAT

    ! Local properties
    REAL(fp), DIMENSION(11) :: GAMMA_BOX
    INTEGER                 :: STATE_LOCAL

    ! Loop variables
    INTEGER                 :: I, J, L, K

    ! Local variables for quantities from Input_Opt
    LOGICAL                 :: prtDebug
    LOGICAL                 :: LHOMNUCNAT
    LOGICAL                 :: LSOLIDPSC
    LOGICAL                 :: LACTIVEH2O

    ! Local variables for quantities from species database
    REAL(fp)                :: NIT_MW_G, HNO3_MW_G, H2O_MW_G

    ! Pointers
    REAL(fp), POINTER       :: Spc (:,:,:,:)
    REAL(fp), POINTER       :: KHETI_SLA(:,:,:,:)
    REAL(f4), POINTER       :: STATE_PSC(:,:,:)

    !=================================================================
    ! CALC_STRAT_AER begins here!
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Copy fields from INPUT_OPT
    LHOMNUCNAT  = Input_Opt%LHOMNUCNAT
    LSOLIDPSC   = Input_Opt%LSOLIDPSC
    LACTIVEH2O  = Input_Opt%LACTIVEH2O

    ! Do we have to print debug output?
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Copy fields from species database
    NIT_MW_G  = State_Chm%SpcData(id_NIT)%Info%emMW_g   ! g/mol
    HNO3_MW_G = State_Chm%SpcData(id_HNO3)%Info%emMW_g  ! g/mol
    H2O_MW_G  = State_Chm%SpcData(id_H2O)%Info%emMW_g   ! g/mol

    ! Initialize GEOS-Chem species array [kg]
    Spc => State_Chm%Species

    ! Initialize sticking coefficients for PSC reactions on SLA
    KHETI_SLA => State_Chm%KHETI_SLA

    ! Initialize gridbox PSC type (see Kirner et al. 2011, GMD)
    STATE_PSC => State_Chm%STATE_PSC

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### UCX: start CALC_STRAT_AER' )
    ENDIF

    ! Partition H2SO4 before proceeding
    CALL CALC_H2SO4_GAS( Input_Opt, State_Chm, State_Grid, State_Met )

    ! Loop over latitude boxes first
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,            J,                  L             ) &
    !$OMP PRIVATE( K,            IS_POLAR,           IS_STRAT      ) &
    !$OMP PRIVATE( PCENTER,      PCENTER_PA,         DENAIR        ) &
    !$OMP PRIVATE( INVAIR,       PSATHNO3,           PSATH2O       ) &
    !$OMP PRIVATE( VOL_NAT,      VOL_ICE                           ) &
    !$OMP PRIVATE( VOL_SLA,      PSATHNO3_SUPERCOOL, TCENTER       ) &
    !$OMP PRIVATE( TINV,         IS_VALID                          ) &
    !$OMP PRIVATE( RAD_AER_BOX,  RHO_AER_BOX                       ) &
    !$OMP PRIVATE( KG_AER_BOX,   NDENS_AER_BOX,      SAD_AER_BOX   ) &
    !$OMP PRIVATE( KG_NAT,       KG_ICE,             KG_NO3        ) &
    !$OMP PRIVATE( GAMMA_BOX,    PSATH2O_SUPERSAT,   H2OSUM        ) &
    !$OMP PRIVATE( H2OPP,        H2O_BOX_S,          H2O_BOX_L     ) &
    !$OMP PRIVATE( H2O_BOX_G,    H2SO4SUM,           HNO3SUM       ) &
    !$OMP PRIVATE( HNO3PP,       HNO3_BOX_S,         HNO3_BOX_L    ) &
    !$OMP PRIVATE( HNO3_BOX_G,   BrNO3SUM,           HBrSUM        ) &
    !$OMP PRIVATE( HOBrSUM,      ClNO3SUM,           HClSUM        ) &
    !$OMP PRIVATE( HOClSUM,      STATE_LOCAL,        HBrGASFRAC    ) &
    !$OMP PRIVATE( HOBrGASFRAC,  HNO3GASFRAC,        HClGASFRAC    ) &
    !$OMP PRIVATE( HOClGASFRAC,  TOFFSET,            W_H2SO4       ) &
    !$OMP PRIVATE( W_H2O,        W_HCl,              W_HOCl        ) &
    !$OMP PRIVATE( W_HBr,        W_HOBr,             W_HNO3        ) &
    !$OMP PRIVATE( HCl_BOX_G,    HCl_BOX_L,          HOCl_BOX_G    ) &
    !$OMP PRIVATE( HOCl_BOX_L,   H2SO4_BOX_G,        HBr_BOX_G     ) &
    !$OMP PRIVATE( HBr_BOX_L,    HOBr_BOX_G,         HOBr_BOX_L    ) &
    !$OMP PRIVATE( H2SO4_BOX_L,  KHET_COMMON,        KHET_SPECIFIC ) &
    !$OMP PRIVATE( VOL_TOT,      BOX_LAT                           ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Now do IS_POLAR check for every grid box separately
       ! and using mid-point instead of edges. This is to be
       ! applicable to cube-sphere systems (ckeller, 1/16/15).
       BOX_LAT  = State_Grid%YMid(I,J)
       IS_POLAR = ( (BOX_LAT <= PSC_MINLAT) .OR. &
                    (BOX_LAT >= PSC_MAXLAT)     )

       ! Get local conditions
       ! Now using dry air partial pressure (ewl, 3/2/15)
       PCENTER    = State_Met%PMID_DRY(I,J,L)
       PCENTER_PA = PCENTER*1.e+2_fp
       TCENTER    = State_Met%T(I,J,L)
       TINV       = 1.e+0_fp/TCENTER

       ! Apply other limits from Kirner et al.
       IS_STRAT =  State_Met%InStratMeso(I,J,L)
       IS_VALID = (IS_POLAR).and.(IS_STRAT).and. &
                  (.not.((PCENTER_PA.lt.PSC_PMIN).or. &
                         (PCENTER_PA.gt.PSC_PMAX)))
       IS_VALID = (IS_VALID.or.PSC_FULL)

       ! Initialize variables
       VOL_NAT      = 0.0e+0_fp
       KG_NAT       = 0.0e+0_fp
       VOL_ICE      = 0.0e+0_fp
       KG_ICE       = 0.0e+0_fp

       RAD_AER_BOX  = 0.0e+0_fp
       RHO_AER_BOX  = 0.0e+0_fp
       KG_AER_BOX   = 0.0e+0_fp
       NDENS_AER_BOX= 0.0e+0_fp

       H2O_BOX_S    = 0.0e+0_fp
       H2O_BOX_L    = 0.0e+0_fp
       H2O_BOX_G    = 0.0e+0_fp
       HNO3_BOX_S   = 0.0e+0_fp
       HNO3_BOX_L   = 0.0e+0_fp
       HNO3_BOX_G   = 0.0e+0_fp
       H2SO4_BOX_L  = 0.0e+0_fp
       H2SO4_BOX_G  = 0.0e+0_fp
       HCl_BOX_L    = 0.0e+0_fp
       HCl_BOX_G    = 0.0e+0_fp
       HOCl_BOX_L   = 0.0e+0_fp
       HOCl_BOX_G   = 0.0e+0_fp
       HBr_BOX_L    = 0.0e+0_fp
       HBr_BOX_G    = 0.0e+0_fp
       HOBr_BOX_L   = 0.0e+0_fp
       HOBr_BOX_G   = 0.0e+0_fp

       STATE_LOCAL  = NINT(STATE_PSC(I,J,L))

       ! Calculate local air density
       DENAIR = AVO * PCENTER_PA / (TCENTER*RSTARG)
       INVAIR = AIRMW / State_Met%AD(I,J,L)

       ! Get available NO3 mass
       KG_NO3 = ( Spc(I,J,L,id_HNO3) * NIT_MW_G / HNO3_MW_g ) + &
                  Spc(I,J,L,id_NIT)

       ! Calculate HNO3 mixing ratio
       HNO3SUM = Spc(I,J,L,id_HNO3) * INVAIR / HNO3_MW_G
       HNO3SUM = HNO3SUM + Spc(I,J,L,id_NIT) * INVAIR / NIT_MW_G

       H2OSUM = Spc(I,J,L,id_H2O) * INVAIR / H2O_MW_G

       ! Calculate partial pressures (Pa)
       HNO3PP = PCENTER_PA * HNO3SUM
       H2OPP  = PCENTER_PA * H2OSUM

       IF (.not.IS_VALID) THEN
          ! No PSCs (SSA only)
          STATE_LOCAL = 0
       ELSEIF (TCENTER.gt.T_MAX) THEN
          ! STS/SSA only
          STATE_LOCAL = 1
       ELSE
          ! Calculate saturation pressures
          PSATHNO3 = (10.e+0_fp**(((TSATHNO3_A+(TCENTER* &
                     TSATHNO3_B))*LOG10(H2OPP*760.0e+0_fp/ATM))+ &
                     (TSATHNO3_C+(TSATHNO3_D*TINV)+(TSATHNO3_E* &
                     TCENTER))))*ATM/760.0e+0_fp

          PSATH2O = 10.0e+0_fp**((TSATH2O_A*TINV)+TSATH2O_B)

          ! Supersaturation requirement for ice
          PSATH2O_SUPERSAT = PSATH2O * Input_Opt%P_ICE_SUPERSAT

          ! If homogeneous NAT nucleation allowed, calculate
          ! threshold saturation pressure
          IF (LHOMNUCNAT) THEN
             ! Calculate as if temperature is (T+T_NAT_SUPERCOOL)
             TOFFSET = TCENTER + Input_Opt%T_NAT_SUPERCOOL
             PSATHNO3_SUPERCOOL = 10.e+0_fp**(((TSATHNO3_A+ &
                     (TOFFSET*TSATHNO3_B))*LOG10(H2OPP*760.0e+0_fp/ATM))+ &
                     (TSATHNO3_C+(TSATHNO3_D/TOFFSET)+(TSATHNO3_E* &
                     TOFFSET)))*ATM/760.0e+0_fp
          ELSE
             ! Make homogeneous nucleation impossible
             TOFFSET = 280.0e+0_fp
             PSATHNO3_SUPERCOOL = HNO3PP + 1e+0_fp
          ENDIF

          IF (LACTIVEH2O) THEN
             ! Only interested in sign
             IF (STATE_LOCAL .gt. 1) THEN
                H2O_BOX_S = H2OPP-PSATH2O
             ELSE
                H2O_BOX_S = H2OPP-PSATH2O_SUPERSAT
             ENDIF
          ELSE
             ! Use local ice mass ratio from GEOS-5 data
             ! Note that we are using this only for its sign!
             H2O_BOX_S = State_Met%QI(I,J,L)
          ENDIF

          ! Ice exists/possible?
          IF (H2O_BOX_S.gt.TINY(0e+0_fp)) THEN
             STATE_LOCAL = 3
             HNO3_BOX_S = HNO3PP-PSATHNO3
          ELSE
             ! If ice not possible could still have NAT
             H2O_BOX_S = 0e+0_fp
             ! 1. Homogeneous nucleation
             IF ((LHOMNUCNAT).and.(STATE_LOCAL.eq.1)) THEN
                HNO3_BOX_S = HNO3PP-PSATHNO3_SUPERCOOL
             ENDIF
             ! 2. Box formerly contained ice or NAT
             IF (STATE_LOCAL.eq.2) THEN
                HNO3_BOX_S = HNO3PP-PSATHNO3
             ENDIF
             IF (HNO3_BOX_S.gt.TINY(1e+0_fp)) THEN
                STATE_LOCAL = 2
             ELSE
                STATE_LOCAL = 1
             ENDIF
          ENDIF
       ENDIF

       ! Store state
       STATE_PSC(I,J,L) = REAL(STATE_LOCAL,f4)

       ! Only continue if we want online solid PSCs
       IF (LSOLIDPSC) THEN

          IF (STATE_LOCAL.eq.3) THEN
             ! Form ice PSCs
             IF (LACTIVEH2O) THEN
                H2O_BOX_S = (H2OPP-PSATH2O)/PCENTER_PA
                H2O_BOX_S = MAX(0e+0_fp,H2O_BOX_S)
                KG_ICE    = H2O_BOX_S*ICEMW*State_Met%AD(I,J,L)/AIRMW
             ELSE
                H2O_BOX_S = State_Met%QI(I,J,L)   * &
                            State_Met%CLDF(I,J,L) * AIRMW / ICEMW
                KG_ICE    = State_Met%QI(I,J,L)   * &
                            State_Met%CLDF(I,J,L) * State_Met%AD(I,J,L)
             ENDIF
             VOL_ICE = H2O_BOX_S * DENAIR * (1.e-3_fp) * ICEMW / &
                       (DENSICE*AVO) ! m3 ice/m3 air
          ELSE
             VOL_ICE = 0e+0_fp
             H2O_BOX_S = 0e+0_fp
             KG_ICE = 0e+0_fp
          ENDIF

          ! Calculate NAT if relevant
          IF ((HNO3_BOX_S.gt.TINY(1e+0_fp)).and.(STATE_LOCAL.ge.2)) THEN
             HNO3_BOX_S = (HNO3PP-PSATHNO3)/PCENTER_PA
             HNO3_BOX_S = MAX(0e+0_fp,HNO3_BOX_S)

             ! Calculate m3 NAT/m3 air
             ! HNO3_BOX_S is the number of moles of HNO3
             ! which will be frozen into HNO3.3H2O (NAT)
             ! Therefore volume calculation must be done
             ! with care!
             VOL_NAT = HNO3_BOX_S * DENAIR * (1.e-3_fp) * NATMW / (DENSNAT*AVO)
             KG_NAT  = HNO3_BOX_S * NATMW * State_Met%AD(I,J,L) / AIRMW

          ELSE
             HNO3_BOX_S = 0e+0_fp
             VOL_NAT = 0e+0_fp
             KG_NAT = 0e+0_fp
          ENDIF

          ! Calculate particle properties
          IF (STATE_LOCAL.lt.2) THEN
             ! Zero all!
             KG_AER_BOX   = 0e+0_fp
             RAD_AER_BOX  = 0e+0_fp
             RHO_AER_BOX  = DENSICE
             NDENS_AER_BOX= 0e+0_fp
          ELSE
             VOL_TOT = VOL_NAT + VOL_ICE
             KG_AER_BOX = KG_NAT + KG_ICE
             RAD_AER_BOX = MIN_RAD
             NDENS_AER_BOX = (3.0e+0_fp*(VOL_TOT)/ &
                             (4.0e+0_fp*PI*(RAD_AER_BOX**3.0e+0_fp)))
             IF (NDENS_AER_BOX.gt.MAX_NDENS) THEN
                NDENS_AER_BOX = MAX_NDENS
                RAD_AER_BOX = (3.0e+0_fp*(VOL_TOT)/ &
                              (4.0e+0_fp*PI*MAX_NDENS))**(1.e+0_fp/3.e+0_fp)
             ENDIF

             ! Prevent div-zero (ckeller, 1/16/15)
             IF ( VOL_TOT > 0.0_fp ) THEN
                RHO_AER_BOX = ((VOL_ICE*DENSICE)+(VOL_NAT*DENSNAT))/VOL_TOT
             ELSE
                RHO_AER_BOX = DENSICE ! Is that correct?
             ENDIF
          ENDIF

          ! Calculate SAD (cm2/cm3)
          SAD_AER_BOX = 4.0e-2_fp * RAD_AER_BOX * &
                        RAD_AER_BOX * NDENS_AER_BOX * PI
       ELSE
          ! Solid PSCs not active
          RAD_AER_BOX = 0e+0_fp
          RHO_AER_BOX = 1000e+0_fp
          KG_AER_BOX = 0e+0_fp
          NDENS_AER_BOX = 0e+0_fp
          SAD_AER_BOX = 0e+0_fp
          HNO3_BOX_S = 0e+0_fp
          H2O_BOX_S = 0e+0_fp
       ENDIF

       ! Store in outer arrays
       RAD_AER(I,J,L,I_SPA)   = RAD_AER_BOX*1.e+2_fp ! cm
       RHO_AER(I,J,L,I_SPA)   = RHO_AER_BOX          ! kg/m3
       KG_AER(I,J,L,I_SPA)    = KG_AER_BOX           ! kg
       NDENS_AER(I,J,L,I_SPA) = NDENS_AER_BOX        !#/m3
       SAD_AER(I,J,L,I_SPA)   = SAD_AER_BOX          ! cm2/cm3

       ! Repartition NIT and HNO3 in strat/meso
       IF (LSOLIDPSC.and.IS_STRAT) THEN

          ! Convert NAT from kg NAT to kg NO3
          Spc(I,J,L,id_NIT) = KG_NAT * NIT_MW_G / NATMW

          ! Remove (kg NO3 as NAT) from total kg NO3
          ! then convert to kg HNO3
          Spc(I,J,L,id_HNO3) = (KG_NO3-Spc(I,J,L,id_NIT)) &
                               * HNO3_MW_G / NIT_MW_G
       ENDIF

       ! Now start liquid aerosol consideration
       ! Start by assuming all non-solid H2O/HNO3 is gaseous
       HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
       HNO3_BOX_L = 0e+0_fp
       H2O_BOX_G  = H2OSUM - H2O_BOX_S
       H2O_BOX_L  = 0e+0_fp

       ! Calculate mixing ratios of other relevant species
       H2SO4SUM = Spc(I,J,L,id_SO4) * INVAIR / &
                  State_Chm%SpcData(id_SO4)%Info%emMW_g
       BrNO3SUM = Spc(I,J,L,id_BrNO3) * INVAIR / &
                  State_Chm%SpcData(id_BrNO3)%Info%emMW_g
       ClNO3SUM = Spc(I,J,L,id_ClNO3) * INVAIR / &
                  State_Chm%SpcData(id_ClNO3)%Info%emMW_g
       HOClSUM  = Spc(I,J,L,id_HOCl) * INVAIR / &
                  State_Chm%SpcData(id_HOCl)%Info%emMW_g
       HClSUM   = Spc(I,J,L,id_HCl) * INVAIR / &
                  State_Chm%SpcData(id_HCl)%Info%emMW_g
       HOBrSUM  = Spc(I,J,L,id_HOBr) * INVAIR / &
                  State_Chm%SpcData(id_HOBr)%Info%emMW_g
       HBrSUM   = Spc(I,J,L,id_HBr) * INVAIR / &
                  State_Chm%SpcData(id_HBr)%Info%emMW_g

       ! H2SO4 gas fraction calculated earlier throughout grid
       ! Consider gaseoues H2SO4 to be unavailable for SLA
       H2SO4_BOX_L = H2SO4SUM * AERFRAC(I,J,L,1)
       H2SO4_BOX_G = H2SO4SUM - H2SO4_BOX_L

       ! Zero local properties
       RHO_AER_BOX   = 1000e+0_fp
       RAD_AER_BOX   = 0e+0_fp
       KG_AER_BOX    = 0e+0_fp
       NDENS_AER_BOX = 0e+0_fp
       SAD_AER_BOX   = 0e+0_fp
       VOL_SLA       = 0e+0_fp
       W_H2O         = 0e+0_fp
       W_H2SO4       = 1e+0_fp

       IF (.not.IS_STRAT) THEN
          ! Use JPL 10-06/Oslo CTM data, where available,
          ! for conventional sulfates/H2SO4
          GAMMA_BOX(1)  = 0.1e+0_fp
          GAMMA_BOX(2)  = 0.0e+0_fp
          GAMMA_BOX(3)  = 0.0e+0_fp
          GAMMA_BOX(4)  = 0.0e+0_fp
          GAMMA_BOX(5)  = 0.3e+0_fp
          GAMMA_BOX(6)  = 0.4e+0_fp
          GAMMA_BOX(7)  = 0.9e+0_fp
          GAMMA_BOX(8)  = 0.0e+0_fp
          GAMMA_BOX(9)  = 0.0e+0_fp
          GAMMA_BOX(10) = 0.2e+0_fp
          GAMMA_BOX(11) = 0.0e+0_fp
       ELSEIF (H2SO4_BOX_L.lt.1e-15_fp) THEN
          ! No aerosol to speak of
          DO K=1,11
             GAMMA_BOX(K) = 0.0e+0_fp
          ENDDO
       ELSE
          IF (STATE_LOCAL.eq.0) THEN
             ! Allow binary H2SO4.nH2O only
             CALL TERNARY( PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L, &
                           0.e+0_fp   ,HClSUM,HOClSUM,HBrSUM,HOBrSUM, &
                           W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr, &
                           HNO3GASFRAC,HClGASFRAC,HOClGASFRAC, &
                           HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)

             ! For safety's sake, zero out HNO3 uptake
             HNO3GASFRAC = 1.e+0_fp
             W_H2O = W_H2O + W_HNO3
             W_HNO3 = 0.e+0_fp
             HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
             HNO3_BOX_L = 0.e+0_fp
          ELSE
             ! As per Buchholz, use only non-NAT HNO3 for STS
             HNO3_BOX_G = HNO3SUM - HNO3_BOX_S
             CALL TERNARY( PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L, &
                           HNO3_BOX_G,HClSUM,HOClSUM,HBrSUM,HOBrSUM, &
                           W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr, &
                           HNO3GASFRAC,HClGASFRAC,HOClGASFRAC, &
                           HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)

             ! Partition HNO3 here for safety
             HNO3_BOX_G = HNO3_BOX_G*HNO3GASFRAC
             HNO3_BOX_L = HNO3SUM - (HNO3_BOX_G+HNO3_BOX_S)
          ENDIF

          ! Partition minor species
          HCl_BOX_G  = HClSUM *HClGASFRAC
          HCl_BOX_L  = HClSUM -HCl_BOX_G
          HOCl_BOX_G = HOClSUM*HOClGASFRAC
          HOCl_BOX_L = HOClSUM-HOCl_BOX_G
          HBr_BOX_G  = HBrSUM *HBrGASFRAC
          HBr_BOX_L  = HBrSUM -HBr_BOX_G
          HOBr_BOX_G = HOBrSUM*HOBrGASFRAC
          HOBr_BOX_L = HOBrSUM-HOBr_BOX_G

          ! Calculate SLA parameters (Grainger 1995)
          SAD_AER_BOX = SLA_VA*(VOL_SLA**0.751e+0_fp)        ! cm2/cm3
          RAD_AER_BOX = SLA_VR*SLA_RR*(VOL_SLA**0.249e+0_fp) ! m
          KG_AER_BOX  = RHO_AER_BOX*VOL_SLA*State_Met%AIRVOL(I,J,L) ! kg

          IF (VOL_SLA.gt.1.e-30_fp) THEN
             ! Approximate particles as spherical for calculation
             ! of aerosol number density
             NDENS_AER_BOX = VOL_SLA*3.e+0_fp/  &
                             (4.e+0_fp*PI*(RAD_AER_BOX**3.e+0_fp))

             ! DENAIR in #/m3 - convert to #/cm3
             ! RHO_AER_BOX in kg/m3 - convert to g/cm3
             ! RAD_AER_BOX in m - convert to cm
             CALL CALC_SLA_GAMMA(DENAIR*1.e-6_fp,TCENTER,PCENTER, &
                                 W_H2SO4,H2OSUM,HClSUM,HBrSUM,HOBrSUM, &
                                 ClNO3SUM,BrNO3SUM,RHO_AER_BOX*1.e-3_fp, &
                                 RAD_AER_BOX*1.e+2_fp,GAMMA_BOX)
          ELSE
             ! Ignore SLA
             DO K=1,11
                GAMMA_BOX(K) = 0.0e+0_fp
             ENDDO
          ENDIF
       ENDIF

       ! Store liquid fractions
       ! Liquid H2O is removed from the sum, then it is assumed
       ! that the pre-calculated solid H2O is taken out of this
       ! liquid total
       H2O_BOX_L = (98.09e+0_fp/18.02e+0_fp)*H2SO4_BOX_L * (W_H2O/W_H2SO4)
       H2O_BOX_L = MAX(0e+0_fp,MIN(H2O_BOX_L-H2O_BOX_S,H2OSUM))
       H2O_BOX_G = MAX(0e+0_fp,H2O_BOX_G-(H2O_BOX_L+H2O_BOX_S))

       ! If very low number density, ignore settling
       AERFRAC(I,J,L,2) = 0e+0_fp
       AERFRAC(I,J,L,3) = 0e+0_fp
       AERFRAC(I,J,L,4) = 0e+0_fp
       AERFRAC(I,J,L,5) = 0e+0_fp
       AERFRAC(I,J,L,6) = 0e+0_fp
       AERFRAC(I,J,L,7) = 0e+0_fp

       IF ((HNO3SUM.gt.1e+0_fp).and.(IS_SAFE_DIV(HNO3_BOX_L,HNO3SUM))) THEN
          AERFRAC(I,J,L,2) = HNO3_BOX_L/HNO3SUM
       ENDIF

       IF ((HClSUM.gt.1e+0_fp).and.(IS_SAFE_DIV(HCl_BOX_L,HClSUM))) THEN
          AERFRAC(I,J,L,3) = HCl_BOX_L/HClSUM
       ENDIF

       IF ((HOClSUM.gt.1e+0_fp).and. (IS_SAFE_DIV(HOCl_BOX_L,HOClSUM))) THEN
          AERFRAC(I,J,L,4) = HOCl_BOX_L/HOClSUM
       ENDIF

       IF ((HBrSUM.gt.1e+0_fp).and.(IS_SAFE_DIV(HBr_BOX_L,HBrSUM))) THEN
          AERFRAC(I,J,L,5) = HBr_BOX_L/HBrSUM
       ENDIF

       IF ((HOBrSUM.gt.1e+0_fp).and.(IS_SAFE_DIV(HOBr_BOX_L,HOBrSUM))) THEN
          AERFRAC(I,J,L,6) = HOBr_BOX_L/HOBrSUM
       ENDIF

       IF ((H2OSUM.gt.1e+0_fp).and.(IS_SAFE_DIV(H2O_BOX_L,H2OSUM))) THEN
          AERFRAC(I,J,L,7) = H2O_BOX_L/H2OSUM
       ENDIF

       ! Send properties to larger array
       ! Convert sticking coefficients into
       ! premultiplying factors ((Kirner)
       KHET_COMMON = 0.25e+0_fp*MOLEC_SPEED(TCENTER,1e+0_fp)

       ! N2O5 + H2O/HCl
       KHET_SPECIFIC= KHET_COMMON*ISR_N2O5
       KHETI_SLA(I,J,L,1)  = GAMMA_BOX(1 )*KHET_SPECIFIC
       KHETI_SLA(I,J,L,2)  = GAMMA_BOX(2 )*KHET_SPECIFIC

       ! ClNO3 + H2O/HCl/HBr
       KHET_SPECIFIC= KHET_COMMON*ISR_ClNO3
       KHETI_SLA(I,J,L,3)  = GAMMA_BOX(3 )*KHET_SPECIFIC
       KHETI_SLA(I,J,L,4)  = GAMMA_BOX(4 )*KHET_SPECIFIC
       KHETI_SLA(I,J,L,5)  = GAMMA_BOX(5 )*KHET_SPECIFIC

       ! BrNO3 + H2O/HCl
       KHET_SPECIFIC= KHET_COMMON*ISR_BrNO3
       KHETI_SLA(I,J,L,6)  = GAMMA_BOX(6 )*KHET_SPECIFIC
       KHETI_SLA(I,J,L,7)  = GAMMA_BOX(7 )*KHET_SPECIFIC

       ! HOCl + HCl/HBr
       KHET_SPECIFIC= KHET_COMMON*ISR_HOCl
       KHETI_SLA(I,J,L,8)  = GAMMA_BOX(8 )*KHET_SPECIFIC
       KHETI_SLA(I,J,L,9)  = GAMMA_BOX(9 )*KHET_SPECIFIC

       ! HOBr + HBr/HCl
       KHET_SPECIFIC= KHET_COMMON*ISR_HOBr
       KHETI_SLA(I,J,L,10) = GAMMA_BOX(10)*KHET_SPECIFIC
       KHETI_SLA(I,J,L,11) = GAMMA_BOX(11)*KHET_SPECIFIC

       RAD_AER(I,J,L,I_SLA)  = RAD_AER_BOX*1.e+2_fp ! cm
       RHO_AER(I,J,L,I_SLA)  = RHO_AER_BOX      ! kg/m3
       KG_AER(I,J,L,I_SLA)   = KG_AER_BOX       ! kg
       NDENS_AER(I,J,L,I_SLA)= NDENS_AER_BOX    ! #/m3
       SAD_AER(I,J,L,I_SLA)  = SAD_AER_BOX      ! cm2/cm3

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    NULLIFY( Spc, STATE_PSC, KHETI_SLA )

  END SUBROUTINE CALC_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kg_strat_aer
!
! !DESCRIPTION: Function KG\_STRAT\_AER returns the calculated mass of a
!  stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays.
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION KG_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
    INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
    INTEGER,INTENT(IN)          :: IAER       ! Aerosol index
                                              ! 1 = SSA (pure H2SO4)
                                              ! 2 = STS
                                              ! 3 = Solid PSC
!
! !REVISION HISTORY:
!  18 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! KG_STRAT_AER begins here!
    !=================================================================

    KG_STRAT_AER = KG_AER(I,J,L,IAER)

  END FUNCTION KG_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rho_strat_aer
!
! !DESCRIPTION: Function RHO\_STRAT\_AER returns the calculated
!  stratospheric aerosol mass density.
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION RHO_STRAT_AER (I,J,L,IAER)
!
! !INPUT PARAMETERS:
!
    INTEGER,INTENT(IN)          :: I,J,L      ! Grid indices
    INTEGER,INTENT(IN)          :: IAER       ! Aerosol index:
                                              ! 1 = Liquid aerosol
                                              ! 2 = Solid aerosol
!
! !REVISION HISTORY:
!  18 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! RHO_STRAT_AER begins here!
    !=================================================================

    RHO_STRAT_AER = RHO_AER(I,J,L,IAER)

  END FUNCTION RHO_STRAT_AER
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_strat_opt
!
! !DESCRIPTION: Subroutine GET\_STRAT\_OPT returns local optical properties
!  for a given stratospheric aerosol. The routine is essentially just an
!  interface to allow external routines to "see" the arrays. However, local
!  aerosol radius is adjusted from liquid to effective radius for aerosol
!  optical depth calculations with liquid aerosols.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_STRAT_OPT (I,J,L,IAER,RAER,REFF,SAD,XSA)
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: I, J, L    ! Grid indices
    INTEGER,  INTENT(IN)  :: IAER       ! Aerosol index
                                        ! 1 = Liquid aerosols
                                        ! 2 = Solid PSC
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) :: REFF       ! Effective radius (cm)
    REAL(fp), INTENT(OUT) :: RAER       ! Physical radius (cm)
    REAL(fp), INTENT(OUT) :: SAD        ! Surface area density (cm2/cm3)
    REAL(fp), INTENT(OUT) :: XSA        ! X-S area density (m2/m3)
!
! !REMARKS:
!  Seb Eastham writes: "I would edit GET_STRAT_OPT so that, when SAD is less
!  than some small value (say 1 nm2/cm3,  which is a vanishingly small surface
!  area), it returns SADSTRAT=XSASTRAT=0.d0 and RAER=REFF=0.1d0 for safety's
!  sake. I think that will prevent code blow-up later on."
!
! !REVISION HISTORY:
!  17 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! GET_STRAT_OPT begins here!
    !=================================================================

    ! Surface area density [cm2/cm3]
    SAD  = SAD_AER(I,J,L,IAER)

    ! Add error check: threshold is 1e-14 cm2/cm3 (bmy, 4/7/15)
    IF ( SAD < 1e-14_fp ) THEN

       !--------------------------------------------------------------
       ! FOR SAFETY'S SAKE: Set outputs to "safe" values and
       ! exit if very small SAD is encountered. (bmy, 4/7/15)
       !--------------------------------------------------------------
       SAD  = 0.0e+0_fp
       XSA  = 0.0e+0_fp
       RAER = 0.1e+0_fp
       REFF = 0.1e+0_fp

    ELSE

       !--------------------------------------------------------------
       ! Otherwise, compute RAER, REFF, XSA normally
       !--------------------------------------------------------------

       ! For SLA, convert liquid radius to effective optical radius
       RAER = RAD_AER(I,J,L,IAER)
       IF (IAER.eq.I_SLA) THEN
          REFF = RAER/SLA_RR
       ELSE
          REFF = RAER
       ENDIF

       ! Standard log-normal distribution approach
       ! Probably OK for PSCs too? Note cm2/cm3 to m2/m3 = x100
       XSA = 0.25e+2_fp*SAD

    ENDIF

  END SUBROUTINE GET_STRAT_OPT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ternary
!
! !DESCRIPTION: Subroutine TERNARY calculates the composition of SSA/STS
!  aerosols using a paramaterization from Carslaw et al. "A Thermodynamic
!  Model of the System HCl-HNO3-H2SO4-H2O, Including Solubilities of HBr,
!  from <200 to 328 K". The bulk of this code was taken directly from the
!  Global Modeling Initiative implementation by David Considine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TERNARY( PCENTER_IN,TCENTER_IN,H2OSUM_IN,H2SO4SUM,       &
                      HNO3SUM,HClSUM,HOClSUM,HBRSum,HOBrSUM,          &
                      W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr, &
                      HNO3GASFRAC,HClGASFRAC,HOClGASFRAC,             &
                      HBrGASFRAC,HOBrGASFRAC,SLA_VOL,SLA_RHO)
!
! !USES:
!
    ! Temporary - for debug
    USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
    USE ERROR_MOD,     ONLY : SAFE_EXP,SAFE_DIV, DEBUG_MSG
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: PCENTER_IN   ! Pressure (hPa)
    REAL(fp), INTENT(IN)  :: TCENTER_IN   ! Temperature (K)
    REAL(fp), INTENT(IN)  :: H2OSUM_IN    ! Total H2O mixing ratio
    REAL(fp), INTENT(IN)  :: H2SO4SUM     ! Liquid H2SO4 mixing ratio
    REAL(fp), INTENT(IN)  :: HNO3SUM      ! Total HNO3 mixing ratio
    REAL(fp), INTENT(IN)  :: HClSUM       ! Total HCl mixing ratio
    REAL(fp), INTENT(IN)  :: HOClSUM      ! Total HOCl mixing ratio
    REAL(fp), INTENT(IN)  :: HBrSUM       ! Total HBr mixing ratio
    REAL(fp), INTENT(IN)  :: HOBrSUM      ! Total HOBr mixing ratio
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) :: W_H2SO4      ! kg H2SO4/kg SLA
    REAL(fp), INTENT(OUT) :: W_H2O        ! kg H2O  /kg SLA
    REAL(fp), INTENT(OUT) :: W_HNO3       ! kg HNO3 /kg SLA
    REAL(fp), INTENT(OUT) :: W_HCl        ! kg HCl  /kg SLA
    REAL(fp), INTENT(OUT) :: W_HOCl       ! kg HOCl /kg SLA
    REAL(fp), INTENT(OUT) :: W_HBr        ! kg HBr  /kg SLA
    REAL(fp), INTENT(OUT) :: W_HOBr       ! kg HOBr /kg SLA
    REAL(fp), INTENT(OUT) :: HNO3GASFRAC  ! Gas fraction HNO3
    REAL(fp), INTENT(OUT) :: HClGASFRAC   ! Gas fraction HCl
    REAL(fp), INTENT(OUT) :: HOClGASFRAC  ! Gas fraction HOCl
    REAL(fp), INTENT(OUT) :: HBrGASFRAC   ! Gas fraction HBr
    REAL(fp), INTENT(OUT) :: HOBrGASFRAC  ! Gas fraction HOBr
    REAL(fp), INTENT(OUT) :: SLA_VOL      ! Aerosol volume (m3/m3)
    REAL(fp), INTENT(OUT) :: SLA_RHO      ! Aer. mass density (kg/m3)
!
! !REVISION HISTORY:
!  19 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Derived inputs
    REAL(fp) :: H2OSUM
    REAL(fp) :: TCENTER
    REAL(fp) :: PCENTER

    ! Partial pressures
    REAL(fp) :: PATMH2O
    REAL(fp) :: PATMHNO3
    REAL(fp) :: PATMHCl
    REAL(fp) :: PATMHOCl
    REAL(fp) :: PATMHBr
    REAL(fp) :: PATMHOBr

    ! Molar densities (mol/m3)
    REAL(fp) :: MOLDENS_H2SO4

    ! Mass totals
    REAL(fp) :: M_H2SO4,M_HNO3
    REAL(fp) :: M_HCl,M_HOCl
    REAL(fp) :: M_HBr,M_HOBr

    ! Binary solutions denoted with BIN
    ! Mole fractions
    REAL(fp) :: X_H2SO4_BIN
    REAL(fp) :: X_HNO3_BIN

    ! Mass fractions
    REAL(fp) :: M_H2SO4_BIN
    REAL(fp) :: M_HNO3_BIN

    ! Effective Henry's Law coefficients
    REAL(fp) :: H_H2SO4_BIN
    REAL(fp) :: H_HNO3_BIN
    REAL(fp) :: H_HCL, H_HOCL
    REAL(fp) :: H_HBr, H_HOBr

    ! Frost point
    REAL(fp) :: T_ICE

    ! Equilibrium vapor pressure
    REAL(fp) :: PVAP_HNO3
    REAL(fp) :: PVAP_HCl
    REAL(fp) :: PVAP_HBr
    REAL(fp) :: PVAP_HOBr

    ! R in m3.atm/(mol K)
    REAL(fp), PARAMETER :: R_ATM = Rd * AIRMW * 1e-3_fp / ATM

    ! Transitional variables
    REAL(fp) :: DENSITY
    REAL(fp) :: TEMPERATURE
    REAL(fp) :: PR
    REAL(fp) :: TR
    REAL(fp) :: TT

    ! Coefficients (q, k) for calculation of H* for H2SO4 and HNO3
    REAL(fp),DIMENSION(10) :: QN,QS
    REAL(fp),DIMENSION(7)  :: KN,KS

    ! Derived parameters
    REAL(fp) :: A,B,C,PHI

    ! Error trapping
    REAL(fp) :: TMP1,TMP2,TMP3,TMP4,TMP5,TMP6

    REAL(fp), PARAMETER :: TNY = 1.0e-28_fp

    ! Debug message
    CHARACTER(LEN=255)   :: DBGMSG

    ! ======================================================================
    DATA QN/14.5734e+0_fp,0.0615994e+0_fp,-1.14895e+0_fp, &
            0.691693e+0_fp,-0.098863e+0_fp, &
            0.0051579e+0_fp,0.123472e+0_fp,-0.115574e+0_fp, &
            0.0110113e+0_fp,0.0097914e+0_fp/
    DATA QS/14.4700e+0_fp,0.0638795e+0_fp,-3.29597e+0_fp, &
            1.778224e+0_fp,-0.223244e+0_fp, &
            0.0086486e+0_fp,0.536695e+0_fp,-0.335164e+0_fp, &
            0.0265153e+0_fp,0.0157550e+0_fp/
    DATA KN/-39.136e+0_fp,6358.4e+0_fp,83.29e+0_fp, &
            -17650.0e+0_fp,198.53e+0_fp, &
            -11948.e+0_fp,-28.469e+0_fp/
    DATA KS/-21.661e+0_fp,2724.2e+0_fp,51.81e+0_fp, &
            -15732.0e+0_fp,47.004e+0_fp, &
            -6969.0e+0_fp,-4.6183e+0_fp/
    ! ======================================================================

    SAVE QN, QS, KN, KS

    !=================================================================
    ! TERNARY begins here!
    !=================================================================

    ! Routine only valid for certain limits
    H2OSUM = MAX(H2OSUM_IN,5.0e-7_fp)
    PCENTER = MAX(PCENTER_IN,5.0e+0_fp)
    TCENTER = TCENTER_IN

    ! Calculate partial pressure of H2O & HNO3
    ! PCENTER is in hPa, so need to convert ATM from Pa to hPa
    PATMH2O  = H2OSUM  * PCENTER / (ATM*1e-2_fp)

    ! Carslaw only valid for 2e-5 < PPH2O < 2e-3 (hPa)
    PATMH2O = MAX(PATMH2O,1.9738465e-8_fp)
    PATMH2O = MIN(PATMH2O,1.9738465e-6_fp)

    PATMHNO3 = HNO3SUM * PCENTER / (ATM*1e-2_fp)
    PATMHCl  = HClSUM  * PCENTER / (ATM*1e-2_fp)
    PATMHOCl = HOClSUM * PCENTER / (ATM*1e-2_fp)
    PATMHBr  = HBrSUM  * PCENTER / (ATM*1e-2_fp)
    PATMHOBr = HOBrSUM * PCENTER / (ATM*1e-2_fp)

    ! Moles of H2SO4 per m3 air
    MOLDENS_H2SO4 = 100.e+0_fp*PCENTER*H2SO4SUM/(RSTARG*TCENTER)

    ! Nucleation temperature of ice
    T_ICE = 2668.70e+0_fp/ &
            (10.4310e+0_fp-(LOG(PATMH2O)+LOG(760.0e+0_fp))/LOG(10.0e+0_fp))

    ! Pressure relation
    PR = LOG(PATMH2O)+18.4e+0_fp

    ! Therefore if temperature lower, set to T_ICE-3
    IF (TCENTER .lt. (T_ICE-3.0e+0_fp)) THEN
       TCENTER = (T_ICE-3.0e+0_fp)
    ENDIF

    IF (TCENTER .lt. 185.0e+0_fp) THEN
       TCENTER = 185.0e+0_fp
    ENDIF

    ! ??
    TT = TCENTER * R_ATM * MOLDENS_H2SO4

    ! Temperature relation
    TR = 1.0e+4_fp/TCENTER-43.4782608e+0_fp

    ! Determine H2SO4/H2O pure solution concentration
    ! Mole fraction of H2SO4 in binary solution
    TMP1 = (KS(1)+KS(2)/TCENTER)** &
           2.0e+0_fp-4.0e+0_fp*(KS(3)+KS(4)/TCENTER)*(KS(5)+KS(6)/ &
           TCENTER+KS(7)*LOG(TCENTER)-LOG(PATMH2O))
    IF ( TMP1 > 0.0_fp ) THEN
       X_H2SO4_BIN = 1.0e+0_fp/(2.0e+0_fp*(KS(3)+KS(4)/TCENTER))* &
           (-KS(1)-KS(2)/TCENTER-(TMP1)**0.5e+0_fp)
    ELSE
       !X_H2SO4_BIN = 1.0e+0_fp/(2.0e+0_fp*(KS(3)+KS(4)/TCENTER))* &
       !   (-KS(1)-KS(2)/TCENTER)
       X_H2SO4_BIN = 0.0_fp
    ENDIF
    !X_H2SO4_BIN = 1.0e+0_fp/(2.0e+0_fp*(KS(3)+KS(4)/TCENTER))* &
    !   (-KS(1)-KS(2)/TCENTER-((KS(1)+KS(2)/TCENTER)** &
    !   2.0e+0_fp-4.0e+0_fp*(KS(3)+KS(4)/TCENTER)*(KS(5)+KS(6)/ &
    !   TCENTER+KS(7)*LOG(TCENTER)-LOG(PATMH2O)))**0.5e+0_fp)

    ! Molality (mol H2SO4/kg H2O) in binary solution
    M_H2SO4_BIN = 55.51e+0_fp*X_H2SO4_BIN/(1.0e+0_fp-X_H2SO4_BIN)

    IF ((TCENTER.le.215.0e+0_fp).AND.(PATMHNO3.gt.TNY)) THEN
       ! Determine HNO3/H2SO4/H2O solution composition
       H_H2SO4_BIN = EXP(QS(1)+QS(2)*TR**2+(QS(3)+QS(4)*TR+ &
            QS(5)*TR**2+QS(6)*TR**3)*PR+(QS(7)+QS(8)*TR+ &
            QS(9)*TR**2)*PR**2+QS(10)*TR*PR**3)
       X_HNO3_BIN=1.0e+0_fp/(2.0e+0_fp*(KN(3)+KN(4)/TCENTER))* &
            (-KN(1)-KN(2)/TCENTER-((KN(1)+KN(2)/TCENTER)** &
            2-4.0e+0_fp*(KN(3)+KN(4)/TCENTER)*(KN(5)+ &
            KN(6)/TCENTER+KN(7)*LOG(TCENTER)-LOG(PATMH2O) &
            ))**0.5e+0_fp)
       M_HNO3_BIN=55.51e+0_fp*X_HNO3_BIN/(1.0e+0_fp-X_HNO3_BIN)
       H_HNO3_BIN=EXP(QN(1)+QN(2)*TR**2+(QN(3)+QN(4)*TR+QN(5)* &
            TR**2+QN(6)*TR**3)*PR+(QN(7)+QN(8)*TR+QN(9)*TR**2)* &
            PR**2+QN(10)*TR*PR**3)
       A=(TT*H_HNO3_BIN*M_HNO3_BIN**2-TT*H_H2SO4_BIN*M_HNO3_BIN* &
            M_H2SO4_BIN-2.0e+0_fp*M_HNO3_BIN**2*M_H2SO4_BIN+ &
            M_HNO3_BIN*M_H2SO4_BIN**2+H_HNO3_BIN*M_HNO3_BIN* &
            M_H2SO4_BIN*PATMHNO3-H_H2SO4_BIN*M_H2SO4_BIN**2* &
            PATMHNO3)/(M_HNO3_BIN**2-M_HNO3_BIN*M_H2SO4_BIN)
       B=M_H2SO4_BIN*(-2.0e+0_fp*TT*H_HNO3_BIN*M_HNO3_BIN+TT* &
            H_H2SO4_BIN*M_H2SO4_BIN+M_HNO3_BIN*M_H2SO4_BIN- &
            H_HNO3_BIN*M_H2SO4_BIN*PATMHNO3)/(M_HNO3_BIN- &
            M_H2SO4_BIN)
       C=(TT*H_HNO3_BIN*M_HNO3_BIN*M_H2SO4_BIN**2)/ &
            (M_HNO3_BIN-M_H2SO4_BIN)
       PHI=ATAN(SQRT(4.0e+0_fp*(A**2-3.0e+0_fp*B)**3-(-2.0e+0_fp*A**3+ &
            9.0e+0_fp*A*B-27.0e+0_fp*C)**2)/(-2.0e+0_fp*A**3+9.0e+0_fp*A &
            *B-27.0e+0_fp*C))
       IF (PHI.lt.0.e+0_fp) THEN
          PHI = PHI + PI
       ENDIF
       M_H2SO4=-1.0e+0_fp/3.0e+0_fp*(A+2.0e+0_fp* &
            SQRT(A**2-3.0e+0_fp*B)* &
            COS((PI+PHI)/3.0e+0_fp))
       M_HNO3=M_HNO3_BIN*(1.0e+0_fp-M_H2SO4/M_H2SO4_BIN)
       W_H2SO4 = M_H2SO4*0.098076e+0_fp/(1.0e+0_fp+M_H2SO4* &
            0.098076e+0_fp+M_HNO3*0.063012e+0_fp)

       ! Check for low H2SO4
       IF (M_H2SO4 .lt. TNY) THEN
          M_H2SO4 = 0.0e+0_fp
          M_HNO3 = M_HNO3_BIN
          W_H2SO4 = 0.0e+0_fp
       ENDIF

       PVAP_HNO3=M_HNO3/(H_HNO3_BIN*M_HNO3/(M_HNO3+ &
                 M_H2SO4)+H_H2SO4_BIN*M_H2SO4/(M_HNO3+M_H2SO4))
       W_HNO3 = (M_HNO3*0.063012e+0_fp)/(1.0e+0_fp+M_H2SO4* &
                 0.098076e+0_fp+M_HNO3*0.063012e+0_fp)

       HNO3GASFRAC=(1.0e+0_fp-(PATMHNO3-PVAP_HNO3)/PATMHNO3)

    ELSE
       ! Solution is pure H2SO4/H2O
       M_H2SO4 = M_H2SO4_BIN
       M_HNO3 = 0.0e+0_fp
       W_H2SO4 = M_H2SO4_BIN*0.098076e+0_fp/(1.0e+0_fp+M_H2SO4_BIN* &
                 0.098076e+0_fp)
       W_HNO3 = 0.0e+0_fp
       PVAP_HNO3 = 0.0e+0_fp
       HNO3GASFRAC=1.0e+0_fp
    ENDIF

    ! ckeller: restrict values to range that prevents float invalids
    IF ( M_H2SO4 < TNY ) M_H2SO4 = 0.0_fp
    IF ( M_HNO3  < TNY ) M_HNO3  = 0.0_fp
    IF ( W_HNO3  < TNY ) W_HNO3  = 0.0_fp

    ! Handle HCl (Luo et al., Vapor pressures of
    ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
    ! temperatures, 1995)
#ifdef MODEL_GEOS
    IF (PATMHCL.gt.0.0e+0_fp) THEN
       TMP1 = W_HNO3+0.610e+0_fp*W_H2SO4
       IF ( TMP1 <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = LOG(TMP1)
       ENDIF
       TMP3 = 36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
              M_H2SO4+63.012e+0_fp* M_HNO3)
       IF ( TMP3 <= 0.0_fp ) THEN
          TMP4 = 0.0_fp
       ELSE
          TMP4 = LOG(TMP3)
       ENDIF
       H_HCL = EXP(-(21.0e+0_fp+46.610e+0_fp*W_HNO3+ &
               4.0690e+0_fp*W_H2SO4- &
               4.8370e+0_fp*SQRT(W_HNO3)+2.1860e+0_fp* &
               SQRT(W_H2SO4)-63.0e+0_fp* &
               W_HNO3**2-40.170e+0_fp*W_HNO3*W_H2SO4- &
               1.5710e+0_fp*W_H2SO4**2)- &
               1.0e+0_fp/TCENTER*(-7437.0e+0_fp-8327.80e+0_fp* &
               W_HNO3+1300.90e+0_fp* &
               W_H2SO4+1087.20e+0_fp*SQRT(W_HNO3)-242.710e+0_fp* &
               SQRT(W_H2SO4)+ &
               18749.0e+0_fp*W_HNO3**2+18500.0e+0_fp*W_HNO3*W_H2SO4+ &
               5632.0e+0_fp*W_H2SO4**2)-TMP2-TMP4)*(ATM*1e-2_fp)
       !H_HCL = EXP(-(21.0e+0_fp+46.610e+0_fp*W_HNO3+ &
       !        4.0690e+0_fp*W_H2SO4- &
       !        4.8370e+0_fp*SQRT(W_HNO3)+2.1860e+0_fp* &
       !        SQRT(W_H2SO4)-63.0e+0_fp* &
       !        W_HNO3**2-40.170e+0_fp*W_HNO3*W_H2SO4- &
       !        1.5710e+0_fp*W_H2SO4**2)- &
       !        1.0e+0_fp/TCENTER*(-7437.0e+0_fp-8327.80e+0_fp* &
       !        W_HNO3+1300.90e+0_fp* &
       !        W_H2SO4+1087.20e+0_fp*SQRT(W_HNO3)-242.710e+0_fp* &
       !        SQRT(W_H2SO4)+ &
       !        18749.0e+0_fp*W_HNO3**2+18500.0e+0_fp*W_HNO3*W_H2SO4+ &
       !        5632.0e+0_fp*W_H2SO4**2)-LOG(W_HNO3+0.610e+0_fp*W_H2SO4)- &
       !        LOG(36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
       !        M_H2SO4+63.012e+0_fp* &
       !        M_HNO3)))*(ATM*1e-2_fp)
       IF ( M_H2SO4 <= 0.0_fp ) THEN
          TMP1 = 0.0_fp
       ELSE
          TMP1 = MOLDENS_H2SO4/M_H2SO4
       ENDIF
       IF ( H_HCL <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HCL
       ENDIF
       IF ( TMP1 == 0.0_fp .AND. TMP2 == 0.0_fp ) THEN
          M_HCl = 0.0_fp
       ELSE
          M_HCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHCL)/(TMP1 + TMP2)
       ENDIF
       !M_HCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHCL)/ &
       !        (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HCL)
       TMP1 = 1.0e+3_fp+98.076e+0_fp*M_H2SO4+63.012e+0_fp*M_HNO3
       IF ( TMP1 <= 0.0_fp ) THEN
          W_HCL=0.0_fp
       ELSE
          W_HCL=M_HCL*36.461e+0_fp/TMP1
       ENDIF
       !W_HCL=M_HCL*36.461e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+
       !      63.012e+0_fp*M_HNO3)
       IF ( H_HCl <= 0.0_fp ) THEN
          PVAP_HCl = 0.0_fp
       ELSE
          PVAP_HCl = M_HCl/H_HCl
       ENDIF
       HClGASFRAC=1.0e+0_fp-(PATMHCL-PVAP_HCL)/PATMHCL
    ELSE
       W_HCl=0.0e+0_fp
       HClGASFRAC=1.0e+0_fp
    ENDIF
#else
    W_HCl=0.0e+0_fp
    HClGASFRAC=1.0e+0_fp
    IF ( PATMHCL .gt. TNY .AND. &
         M_H2SO4 .gt. TNY .AND. &
         M_HNO3  .gt. TNY       ) THEN
       TMP1 = W_HNO3+0.610e+0_fp*W_H2SO4
       IF ( TMP1 > TNY ) THEN
          TMP2 = LOG(TMP1)
          TMP3 = 36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
                 M_H2SO4+63.012e+0_fp*M_HNO3)
          IF ( TMP3 > TNY ) THEN
             TMP4 = LOG(TMP3)
             H_HCL = EXP(-(21.0e+0_fp+46.610e+0_fp*W_HNO3+ &
                     4.0690e+0_fp*W_H2SO4- &
                     4.8370e+0_fp*SQRT(W_HNO3)+2.1860e+0_fp* &
                     SQRT(W_H2SO4)-63.0e+0_fp* &
                     W_HNO3**2-40.170e+0_fp*W_HNO3*W_H2SO4- &
                     1.5710e+0_fp*W_H2SO4**2)- &
                     1.0e+0_fp/TCENTER*(-7437.0e+0_fp-8327.80e+0_fp* &
                     W_HNO3+1300.90e+0_fp* &
                     W_H2SO4+1087.20e+0_fp*SQRT(W_HNO3)-242.710e+0_fp* &
                     SQRT(W_H2SO4)+ &
                     18749.0e+0_fp*W_HNO3**2+18500.0e+0_fp*W_HNO3*W_H2SO4+ &
                     5632.0e+0_fp*W_H2SO4**2)-TMP2-TMP4)*(ATM*1e-2_fp)
             IF ( H_HCL > TNY ) THEN
                TMP1 = MOLDENS_H2SO4/M_H2SO4
                TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HCL
                IF ( TMP1+TMP2 > TNY ) THEN
                   M_HCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHCL)/(TMP1 + TMP2)
                   TMP1 = 1.0e+3_fp+98.076e+0_fp*M_H2SO4+63.012e+0_fp*M_HNO3
                   W_HCL=M_HCL*36.461e+0_fp/TMP1
                   PVAP_HCl = M_HCl/H_HCl
                   HClGASFRAC=1.0e+0_fp-(PATMHCL-PVAP_HCL)/PATMHCL
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF
#endif

    ! Now HOCl
    W_HOCl=0.0e+0_fp
    HOClGASFRAC=1.0e+0_fp
    IF (PATMHOCl>TNY .AND. M_H2SO4>TNY ) THEN
       H_HOCl=EXP(6.49460e+0_fp-(-0.041070e+0_fp+54.56e+0_fp/TCENTER)* &
              (M_H2SO4+M_HNO3)-5862.0e+0_fp*(1.0e+0_fp/298.15e+0_fp- &
              1.0e+0_fp/TCENTER))

#ifdef MODEL_GEOS
       IF ( M_H2SO4 <= 0.0_fp ) THEN
          TMP1 = 0.0_fp
       ELSE
          TMP1 = MOLDENS_H2SO4/M_H2SO4
       ENDIF
       IF ( H_HOCL <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HOCL
       ENDIF
       IF ( TMP1 == 0.0_fp .AND. TMP2 == 0.0_fp ) THEN
          M_HOCl = 0.0_fp
       ELSE
          M_HOCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHOCl)/(TMP1+TMP2)
       ENDIF
       !M_HOCl=(1.0e+0_fp/R_ATM/TCENTER*PATMHOCl)/ &
       !   (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HOCL)
       W_HOCL=M_HOCL*52.46e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+ &
              63.012e+0_fp*M_HNO3)
       ! Realistically expect no gas phase removal
       HOCLGASFRAC=1.0e+0_fp
    ELSE
       W_HOCl=0.0e+0_fp
       HOClGASFRAC=1.0e+0_fp
#else
       IF ( H_HOCL > TNY ) THEN
          TMP1 = MOLDENS_H2SO4/M_H2SO4
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HOCL
          IF ( TMP1+TMP2 > TNY ) THEN
             M_HOCl = (1.0e+0_fp/R_ATM/TCENTER*PATMHOCl)/(TMP1+TMP2)
             W_HOCL = M_HOCL*52.46e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+ &
                      63.012e+0_fp*M_HNO3)
             ! Realistically expect no gas phase removal
             HOCLGASFRAC=1.0e+0_fp
          ENDIF
       ENDIF
#endif
    ENDIF

    ! Now HBr (Luo et al., Vapor pressures of
    ! H2SO4/HNO3/HCl/HBr/H2O solutions to low stratospheric
    ! temperatures, 1995)
#ifdef MODEL_GEOS
    IF (PATMHBr.gt.0.0e+0_fp) THEN
       TMP1 = W_HNO3+0.410e+0_fp*W_H2SO4
       IF ( TMP1 <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = LOG(TMP1)
       ENDIF
       TMP3 = 36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
              M_H2SO4+63.012e+0_fp*M_HNO3)
       IF ( TMP3 > 0.0_fp ) THEN
          TMP4 = LOG(TMP3)
       ELSE
          TMP4 = 0.0_fp
       ENDIF
       H_HBr = EXP(-(17.83e+0_fp+1.02e+0_fp*W_HNO3-1.08e+0_fp*W_H2SO4+ &
               3.9e+0_fp*SQRT(W_HNO3)+4.38e+0_fp*SQRT(W_H2SO4)-8.87e+0_fp* &
               W_HNO3**2-17.0e+0_fp*W_HNO3*W_H2SO4+3.73e+0_fp*W_H2SO4**2)- &
               1.0e+0_fp/TCENTER*(-8220.5e+0_fp-362.76e+0_fp* &
               W_HNO3+658.93e+0_fp* &
               W_H2SO4-914.0e+0_fp*SQRT(W_HNO3)-955.3e+0_fp*SQRT(W_H2SO4)+ &
               9976.6e+0_fp*W_HNO3**2+19778.5e+0_fp*W_HNO3*W_H2SO4+ &
               7680.0e+0_fp*W_H2SO4**2)-TMP2-TMP4)*(ATM*1e-2_fp)
       !H_HBr = EXP(-(17.83e+0_fp+1.02e+0_fp*W_HNO3-1.08e+0_fp*W_H2SO4+ &
       !        3.9e+0_fp*SQRT(W_HNO3)+4.38e+0_fp*SQRT(W_H2SO4)-8.87e+0_fp* &
       !        W_HNO3**2-17.0e+0_fp*W_HNO3*W_H2SO4+3.73e+0_fp*W_H2SO4**2)- &
       !        1.0e+0_fp/TCENTER*(-8220.5e+0_fp-362.76e+0_fp* &
       !        W_HNO3+658.93e+0_fp* &
       !        W_H2SO4-914.0e+0_fp*SQRT(W_HNO3)-955.3e+0_fp*SQRT(W_H2SO4)+ &
       !        9976.6e+0_fp*W_HNO3**2+19778.5e+0_fp*W_HNO3*W_H2SO4+ &
       !        7680.0e+0_fp*W_H2SO4**2)-LOG(W_HNO3+0.410e+0_fp*W_H2SO4)- &
       !        LOG(36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
       !        M_H2SO4+63.012e+0_fp* &
       !        M_HNO3)))*(ATM*1e-2_fp)
       IF ( M_H2SO4 <= 0.0_fp ) THEN
          TMP1 = 0.0_fp
       ELSE
          TMP1 = MOLDENS_H2SO4/M_H2SO4
       ENDIF
       IF ( H_HBr <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HBr
       ENDIF
       IF ( TMP1==0.0_fp .AND. TMP2==0.0_fp ) THEN
          M_HBr = 0.0_fp
       ELSE
          M_HBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHBr)/(TMP1+TMP2)
       ENDIF
       !M_HBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHBr)/ &
       !        (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HBr)
       W_HBr = M_HBr*80.91e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+ &
               63.012e+0_fp*M_HNO3)
       IF ( H_HBr <= 0.0_fp ) THEN
          PVAP_HBr = 0.0_fp
       ELSE
          PVAP_HBr = M_HBr/H_HBr
       ENDIF
       HBrGASFRAC=1.0e+0_fp-(PATMHBr-PVAP_HBr)/PATMHBr
    ELSE
       W_HBr=0.0e+0_fp
       HBrGASFRAC=1.0e+0_fp
    ENDIF
#else
    W_HBr=0.0e+0_fp
    HBrGASFRAC=1.0e+0_fp
    IF (PATMHBr > TNY .AND. M_H2SO4 > TNY ) THEN
       TMP1 = W_HNO3+0.410e+0_fp*W_H2SO4
       IF ( TMP1 > TNY ) THEN
          TMP2 = LOG(TMP1)
          TMP3 = 36.461e+0_fp/(1000.0e+0_fp+98.076e+0_fp* &
                 M_H2SO4+63.012e+0_fp*M_HNO3)
          IF ( TMP3 > TNY ) THEN
             TMP4 = LOG(TMP3)
             H_HBr = EXP(-(17.83_fp+1.02_fp*W_HNO3-1.08_fp*W_H2SO4+ &
                   3.9e+0_fp*SQRT(W_HNO3)+4.38e+0_fp*SQRT(W_H2SO4)-8.87e+0_fp* &
                   W_HNO3**2-17.0e+0_fp*W_HNO3*W_H2SO4+3.73e+0_fp*W_H2SO4**2)- &
                   1.0e+0_fp/TCENTER*(-8220.5e+0_fp-362.76e+0_fp* &
                   W_HNO3+658.93e+0_fp* &
                   W_H2SO4-914.0e+0_fp*SQRT(W_HNO3)-955.3e+0_fp*SQRT(W_H2SO4)+ &
                   9976.6e+0_fp*W_HNO3**2+19778.5e+0_fp*W_HNO3*W_H2SO4+ &
                   7680.0e+0_fp*W_H2SO4**2)-TMP2-TMP4)*(ATM*1e-2_fp)
             IF ( H_HBr > TNY ) THEN
                TMP1 = MOLDENS_H2SO4/M_H2SO4
                TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HBr
                IF ( TMP1+TMP2 > TNY ) THEN
                   M_HBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHBr)/(TMP1+TMP2)
                   W_HBr = M_HBr*80.91e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+ &
                           63.012e+0_fp*M_HNO3)
                   PVAP_HBr = M_HBr/H_HBr
                   HBrGASFRAC=1.0e+0_fp-(PATMHBr-PVAP_HBr)/PATMHBr
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF
#endif

    ! Finally HOBr (Hanson and Ravishankara, Heterogeneous
    ! chemistry of Bromine species in sulfuric acid under
    ! stratospheric conditions, 1995)
    W_HOBr=0e+0_fp
    HOBrGASFRAC=1.0e+0_fp
    IF (PATMHOBr > TNY .and. M_H2SO4 > TNY ) THEN
       ! Hanson and Ravishankara state that the volume-based
       ! Henry's Law coefficient for HOBr in H2SO4 is 10^6 M/atm.
       ! The molality-based Henry's law constant, H_HOBr, is
       ! therefore:
#ifdef MODEL_GEOS
       IF ( M_H2SO4 <= 0.0_fp ) THEN
          H_HOBr = 0.0_fp
       ELSE
          H_HOBr = (1.0e+6_fp) * MOLDENS_H2SO4 / M_H2SO4
       ENDIF
       IF ( M_H2SO4 <= 0.0_fp ) THEN
          TMP1 = 0.0_fp
       ELSE
          TMP1 = MOLDENS_H2SO4/M_H2SO4
       ENDIF
       IF ( H_HOBr <= 0.0_fp ) THEN
          TMP2 = 0.0_fp
       ELSE
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HOBr
       ENDIF
       IF ( TMP1 == 0.0_fp .AND. TMP2 == 0.0_fp ) THEN
          M_HOBr = 0.0_fp
       ELSE
          M_HOBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHOBr)/(TMP1 + TMP2)
       ENDIF
       M_HOBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHOBr)/ &
                (MOLDENS_H2SO4/M_H2SO4 + 1.0e+0_fp/R_ATM/TCENTER/H_HOBr)
       W_HOBr = M_HOBr*96.911e+0_fp/(1.0e+3_fp+98.076e+0_fp*M_H2SO4+ &
                63.012e+0_fp*M_HNO3)
       IF ( H_HOBr <= 0.0_fp ) THEN
          PVAP_HOBr = 0.0_fp
       ELSE
          PVAP_HOBr = M_HOBr/H_HOBr
       ENDIF
       HOBrGASFRAC=1.0e+0_fp-(PATMHOBr-PVAP_HOBr)/PATMHOBr
    ELSE
       W_HOBr=0e+0_fp
       HOBrGASFRAC=1.0e+0_fp
    ENDIF
#else
       H_HOBr = (1.0e+6_fp) * MOLDENS_H2SO4 / M_H2SO4
       IF ( H_HOBr > TNY ) THEN
          TMP1 = MOLDENS_H2SO4/M_H2SO4
          TMP2 = 1.0e+0_fp/R_ATM/TCENTER/H_HOBr
          IF ( TMP1+TMP2 > TNY ) THEN
             M_HOBr = (1.0e+0_fp/R_ATM/TCENTER*PATMHOBr)/(TMP1+TMP2)
             W_HOBr = M_HOBr*96.911_fp/(1.0e+3_fp+98.076_fp*M_H2SO4+ &
                      63.012e+0_fp*M_HNO3)
             PVAP_HOBr = M_HOBr/H_HOBr
             HOBrGASFRAC=1.0e+0_fp-(PATMHOBr-PVAP_HOBr)/PATMHOBr
          ENDIF
       ENDIF
    ENDIF
#endif

    ! Take W_H2O as remainder
    W_H2O = 1.e+0_fp-(W_H2SO4+W_HNO3+W_HCl+W_HOCl+W_HBr+W_HOBr)

#ifdef MODEL_GEOS
    ! restrict values to range that prevents invalids (12/29/17)
    IF ( M_H2SO4 < 1.0e-30_fp ) M_H2SO4 = 0.0_fp
    IF ( M_HNO3  < 1.0e-30_fp ) M_HNO3  = 0.0_fp
#endif

    ! Aerosol mass density in kg/m3 aerosol
    SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER)

    ! Aerosol volume in m3/m3 air
    IF ( W_H2SO4 < TNY .OR. SLA_RHO < TNY ) THEN
       SLA_VOL = 0.0_fp
    ELSE
       SLA_VOL = (MOLDENS_H2SO4*98.076e+0_fp/W_H2SO4/SLA_RHO)*1.e-3_fp
    ENDIF

  END SUBROUTINE TERNARY
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carslaw_density
!
! !DESCRIPTION: Function CARSLAW\_DENSITY determines the density of a
!  sol'n through a relationship from Carslaw et al.. Result is in kg/m3.
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION CARSLAW_DENSITY(CS,CN,T)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: CS         ! H2SO4 molality (mol H2SO4/kg solvent)
    REAL(fp), INTENT(IN) :: CN         ! HNO3 molality (mol HNO3/kg solvent)
    REAL(fp), INTENT(IN) :: T          ! Temperature (K)
!
! !REVISION HISTORY:
!  19 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)                :: DENSS,DENSN

    !=================================================================
    ! CARSLAW_DENSITY begins here!
    !=================================================================

    DENSS=1000.0e+0_fp+123.64e+0_fp*CS-5.6e-4_fp*CS*T**2 &
          -29.54e+0_fp*CS**1.5e+0_fp + 1.814e-4_fp*CS**1.5e+0_fp &
          *T**2+2.343e+0_fp*CS**2  -1.487e-3_fp*CS**2*T &
          -1.324e-5_fp*CS**2*T**2

    DENSN=1000.0e+0_fp+85.107e+0_fp*CN-5.043e-4_fp*CN*T**2 &
          -18.96e+0_fp*CN**1.5e+0_fp + 1.427e-4_fp*CN**1.5e+0_fp &
          *T**2+1.458e+0_fp*CN**2  -1.198e-3_fp*CN**2*T &
          -9.703e-6_fp*CN**2*T**2

    ! Error trap for zeros (ckeller, 12/29/17)
    IF ( CS == 0.0_fp .AND. CN == 0.0_fp ) THEN
       CARSLAW_DENSITY = 0.0_fp
    ELSE
       CARSLAW_DENSITY=1.0e+0_fp/((1.0e+0_fp/DENSS*CS/(CS+CN) &
                      +1.0e+0_fp/DENSN*CN/(CS+CN)))
    ENDIF
    RETURN

  END FUNCTION CARSLAW_DENSITY
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_fallvel
!
! !DESCRIPTION: Function CALC\_FALLVEL calculates the terminal velocity of a
!  solid particle.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CALC_FALLVEL(DENSITY,RADIUS,TCENTER,PCENTER) RESULT(VEL)
!
! !USES:
!
    USE ERROR_MOD,       ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    REAL(fp),INTENT(IN)   :: RADIUS  ! Particle radius (cm)
    REAL(fp),INTENT(IN)   :: DENSITY ! Particle density (kg/m3)
    REAL(fp),INTENT(IN)   :: TCENTER ! Local temperature (K)
    REAL(fp),INTENT(IN)   :: PCENTER ! Local pressure (kPa)
!
! !OUTPUT VARIABLES:
!
    REAL(fp)              :: VEL ! Fall velocity (m/s)
!
! !REMARKS:
! (1) A remark
!
! !REVISION HISTORY:
!  11 Aug 2012 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Vy ! Intermediate velocity (m/s)
    REAL(fp),PARAMETER    :: eta=6.45e-8_fp ! Constant (kg/(msK))
    REAL(fp)              :: val_x ! Dimensionless variable
    REAL(fp),DIMENSION(3) :: alpha ! Auxiliary variables
    REAL(fp)              :: PR ! Pressure times radius

    !=================================================================
    ! CALC_FALLVEL begins here!
    !=================================================================

    DATA ALPHA/1.49e-5_fp,5.02e-6_fp,2.64e-5_fp/

    ! Sanity check
    IF ((RADIUS.le.0.e+0_fp).or.(DENSITY.le.0.e+0_fp)) THEN
       VEL=0.e+0_fp
    ELSE
       ! PCENTER (kPa -> Pa) = *1.d3
       ! RADIUS  (cm  -> m ) = *1.d-2
       ! Therefore multiply PR by 10
       PR = PCENTER * RADIUS * 10e+0_fp
       VAL_X = -1.0e+0_fp*PR/(ALPHA(3)*TCENTER)
       VAL_X = ALPHA(2)*TCENTER*EXP(VAL_X)/PR
       VAL_X = 1.0e+0_fp + VAL_X + (ALPHA(1)*TCENTER/PR)
       Vy = g0*DENSITY*RADIUS*RADIUS*(1.e-4_fp)/(4.5*ETA*TCENTER)
       VEL = 0.893e+0_fp * Vy * VAL_X
    ENDIF

    ! Velocities should be of the order of 0.1 m/s
    IF (VEL.gt.10.e+0_fp) THEN
       CALL ERROR_STOP(' Excessive fall velocity? ', ' CALC_FALLVEL, UCX_mod')
    ENDIF

  END FUNCTION CALC_FALLVEL
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cacl_sla_gamma
!
! !DESCRIPTION: Subroutine CALC\_SLA\_GAMMA calculates 11 different sticking 
!  coefficients on the surface of local stratospheric liquid aerosols,
!  relevant to each of the 11 reactions listed in Kirner's paper.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_SLA_GAMMA( NDENS, T, P, WT_FRC, H2OSUM, HClSUM, &
                             HBrSUM, HOBrSUM, ClNO3SUM, BrNO3SUM, &
                             RHO, ARAD, RXNGAMMA )
!
! !USES:
!
    USE ERROR_MOD,     ONLY : IT_IS_NAN,ERROR_STOP     ! Test for NaN
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: NDENS    ! Air number density (molec/cm3)
    REAL(fp), INTENT(IN)  :: T        ! Temperature (K)
    REAL(fp), INTENT(IN)  :: P        ! Pressure (hPa)
    REAL(fp), INTENT(IN)  :: WT_FRC   ! Weight fraction of H2SO4 (kg/kg)
    REAL(fp), INTENT(IN)  :: H2OSUM   ! H2O mixing ratio
    REAL(fp), INTENT(IN)  :: HClSUM   ! HCl mixing ratio
    REAL(fp), INTENT(IN)  :: HBrSUM   ! HBr mixing ratio
    REAL(fp), INTENT(IN)  :: HOBrSUM  ! HOBr mixing ratio
    REAL(fp), INTENT(IN)  :: ClNO3SUM ! ClNO3 mixing ratio
    REAL(fp), INTENT(IN)  :: BrNO3SUM ! BrNO3 mixing ratio
    REAL(fp), INTENT(IN)  :: RHO      ! STS density (g/cm3)
    REAL(fp), INTENT(IN)  :: ARAD     ! SLA radius (cm)
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) :: RXNGAMMA(11) ! Premultiplying factors
!
! !REVISION HISTORY:
!  10 Oct 2012 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(f8)               :: WT      ! Weight percentage H2SO4 (100*kg/kg)
    REAL(f8)               :: H2OPP   ! Partial pressure of H2O (hPa)
    REAL(f8)               :: HClPP   ! Partial pressure of HCl (atm)
    REAL(f8)               :: HBrPP   ! Partial pressure of HBr (atm)
    REAL(f8)               :: HOBrPP  ! Partial pressure of HOBr (atm)
    REAL(f8)               :: ClNO3PP ! Partial p. of ClONO2 (atm)
    REAL(f8)               :: BrNO3PP ! Partial p. of BrONO2 (atm)
    REAL(f8)               :: PSATH2O ! Water vapor sat. pressure (hPa)
    REAL(f8)               :: ACTH2O  ! Activity of water
    REAL(f8)               :: MOLAL   ! Molality of H2SO4 (mol H2SO4/kg solvent)
    REAL(f8), DIMENSION(3) :: Z       ! Parameters for H2SO4 sol'n
    REAL(f8)               :: M_H2SO4 ! Mass of H2SO4

    ! HOBr parameters
    REAL(f8)               :: c_HOBr
    REAL(f8)               :: SHOBr
    REAL(f8)               :: HHOBr
    REAL(f8)               :: DHOBr
    REAL(f8)               :: kHOBr_HCl
    REAL(f8)               :: GHOBrrxn
    REAL(f8)               :: lHOBr
    REAL(f8)               :: fHOBr
    REAL(f8)               :: gHOBr_HCl

    ! HOCl parameters
    REAL(f8)               :: c_HOCl
    REAL(f8)               :: SHOCl
    REAL(f8)               :: HHOCl
    REAL(f8)               :: DHOCl
    REAL(f8)               :: kHOCl_HCl
    REAL(f8)               :: GHOClrxn
    REAL(f8)               :: lHOCl
    REAL(f8)               :: fHOCl
    REAL(f8)               :: gHOCl_HCl

    ! ClNO3 parameters
    REAL(f8)               :: c_ClNO3
    REAL(f8)               :: SClNO3
    REAL(f8)               :: HClNO3
    REAL(f8)               :: DClNO3
    REAL(f8)               :: GClNO3rxn
    REAL(f8)               :: lClNO3
    REAL(f8)               :: fClNO3
    REAL(f8)               :: gClNO3
    REAL(f8)               :: gClNO3_HCl
    REAL(f8)               :: gClNO3_H2O

    ! N2O5 parameters
    REAL(f8), DIMENSION(3) :: AK

    ! Other parameters
    REAL(f8) :: kH2O,kH,khdr,GbH2O,HHCl,MHCl,kHCl,GbHCl,Gs,FHCl,Gsp
    REAL(f8) :: GbHClp, Gb, khydr, kII, k_dl

    ! Interim variables
    REAL(f8)               :: X,A,H,T_THRESHOLD,aH
    REAL(f8), PARAMETER    :: MAX_T_DIFF = 6.0e+0_fp

    ! Control whether to run calculations
    LOGICAL              :: HClOK, HOBrOK

    ! Debug variables
    INTEGER              :: I
    CHARACTER(LEN=255)   :: DBGMSG

    !=================================================================
    ! CALC_SLA_GAMMA begins here!
    !=================================================================

    PSATH2O = EXP(18.452406985e+0_f8-3505.1578807e+0_f8 &
              /T-330918.55082e+0_f8/(T*T) &
              +12725068.262e+0_f8/(T*T*T))       ! Saturation pressure of H2O
    H2OPP = H2OSUM * P                  ! Partial pressure of H2O
    ACTH2O = MAX((H2OPP/PSATH2O),1.0e+0_f8) ! Water activity

    ! Calculate molality of solution
    !WT = MIN(100.0e+0_fp,100.0e+0_fp*WT_FRC) ! Convert from fraction to %
    WT = 100.0e+0_f8*WT_FRC ! Convert from fraction to %
    MOLAL = 1000.0e+0_f8 * (WT/98.0e+0_f8/(100.0-WT))

    ! Parameters for H2SO4 solution
    !----------------------------------------------------------
    ! The solution density is calculated earlier, including
    ! contributions from HNO3. This code treats it as a binary
    ! solution - so far this is just a kludge. Need to update
    ! all this code to acknowledge the presence of at least
    ! HNO3 (e.g. X is still calculated based on pure H2O
    ! solvent!)
    !----------------------------------------------------------
    !Z(1) =   0.12364e+0_fp-5.6d-7*T*T
    !Z(2) =   -0.02954e+0_fp+1.814d-7*T*T
    !Z(3) =   2.343d-3-1.487d-6*T-1.324d-8*T*T
    !RHO  =   1.0e+0_fp+Z(1)*MOLAL+Z(2)*MOLAL**1.5+Z(3)*MOLAL*MOLAL
    !----------------------------------------------------------
    M_H2SO4 = RHO*WT/9.8 ! Molality (mol H2SO4/kg solvent)
    X       = WT/(WT+(100.-WT)*98./18.)
    A       = 169.5+5.18*WT-0.0825*WT*WT+3.27e-3_f8*WT*WT*WT
    T_THRESHOLD = 144.11+0.166*WT-0.015*WT*WT+2.18e-4_f8*WT*WT*WT
    IF ((T-T_THRESHOLD).gt.MAX_T_DIFF) THEN
       H = A*T**(-1.43)*EXP(448./(T-T_THRESHOLD))
    ELSE
       H = A*T**(-1.43)*EXP(448./MAX_T_DIFF)
    ENDIF

    aH = EXP(60.51-0.095*WT+0.0077*WT*WT-1.61e-5_f8*WT*WT*WT &
         -(1.76+2.52e-4_f8*WT*WT)*SQRT(T) + &
         (-805.89+253.05*WT**0.076)/SQRT(T))

    HClPP   = HClSUM  *P/(ATM*1e-2_f8) ! Note atm, not hPa
    ClNO3PP = ClNO3SUM*P/(ATM*1e-2_f8)
    BrNO3PP = BrNO3SUM*P/(ATM*1e-2_f8)
    HOBrPP  = HOBrSUM *P/(ATM*1e-2_f8) ! Note atm, not hPa

    ! Should we bother running calculations?
    HClOK  = (HClPP  .gt. 1.e-30_f8)
    HOBrOK = (HOBrPP .gt. 1.e-30_f8)

    ! Reaction 1. N2O5 + H2O (hydrolysis of N2O5)
    AK(1)=-25.5265-0.133188*WT+0.0093084*WT**2-9.0194E-5_f8*WT**3
    AK(2)=9283.76+115.345*WT-5.19258*WT**2+0.0483464*WT**3
    AK(3)=-851801-22191.2*WT+766.916*WT**2-6.85427*WT**3
    RXNGAMMA(1)=exp(AK(1)+AK(2)/T+AK(3)/T**2)

    ! Reaction 2. N2O5 + HCl
    ! JPL 10-06 suggests near-zero gamma
    RXNGAMMA(2) = TINY(1e+0_f8)

    ! Reactions 3/4. ClNO3 + H2O/HCl
    ! Now od only if HCl concentrations are large enough
    ! (to avoid div-by-zero errors), ckeller, 2/10/15.
    IF (HClOK) THEN
       c_ClNO3     =  1474.e+0_f8*SQRT(T)
       SClNO3      =  0.306e+0_f8+24.e+0_f8/T
       HClNO3      =  1.6e-6_f8*EXP(4710.e+0_f8/T) * EXP(-SClNO3*M_H2SO4)
       DClNO3      =  5e-8_f8*T/h
       kH2O        =  1.95e+10_f8*EXP(-2800.e+0_f8/T)
       kH          =  1.22e+12_f8*EXP(-6200.e+0_f8/T)
       khydr       =  kH2O*ACTH2O + kH*aH*ACTH2O
       GbH2O       =  4*HClNO3*0.082*T*SQRT(DClNO3*khydr) / c_ClNO3
       HHCl        =  (0.094e+0_f8-0.61e+0_f8*X+1.2e+0_f8*X*X) &
                      * EXP(-8.68+(8515-10718*X**0.7)/T)
       MHCl        =  HHCl *HClPP
       kHCl        =  7.9e+11_f8*aH*DClNO3*MHCl
       lClNO3      =  SQRT(DClNO3/(khydr+kHCl))
       if (lClNO3.gt.(1.e+5_f8*arad)) then
          ! Limiting rate
          fClNO3   =  arad/(3.e+0_f8*lClNO3)
       else
          fClNO3   =  1.e+0_f8/tanh(arad/lClNO3)- lClNO3/arad
       endif
       GClNO3rxn   =  fClNO3*GbH2O *SQRT(1.e+0_f8+kHCl/khydr)
       GbHCl       =  GClNO3rxn* kHCl/(kHCl+ khydr)
       Gs          =  66.12e+0_f8*EXP(-1374.e+0_f8/T)*HClNO3*MHCl
       FHCl        =  1.e+0_f8/(1.e+0_f8+0.612e+0_f8*(Gs+GbHCl)* &
                      ClNO3PP/ HClPP)
       Gsp         =  FHCl*Gs
       GbHClp      =  FHCl*GbHCl
       Gb          =  GbHClp  + GClNO3rxn* khydr/( kHCl+ khydr)
       ! Catch for zero (ckeller, 12/29/17)
       IF ( Gsp == 0.0_fp .AND. Gb == 0.0_fp ) THEN
          gClNO3_HCl  =  TINY(1e+0_f8)
          gClNO3_H2O  =  TINY(1e+0_f8)
       ELSE
          gClNO3      =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/(Gsp + Gb))
          gClNO3_HCl  =  gClNO3 *(Gsp + GbHClp)/(Gsp + Gb)
          gClNO3_H2O  =  gClNO3 - gClNO3_HCl
       ENDIF

       !IF (HClOK) THEN
       RXNGAMMA(3) =  gClNO3_H2O
       RXNGAMMA(4) =  gClNO3_HCl
    ELSE
       RXNGAMMA(3) = TINY(1e+0_f8)
       RXNGAMMA(4) = TINY(1e+0_f8)
    ENDIF

    ! Reaction 5. ClNO3 + HBr
    ! Not present in JPL 10-06 for H2SO4
    RXNGAMMA(5) = TINY(1e+0_f8)

    ! Reaction 6. BrNO3 + H2O
    !RXNGAMMA(6) = 1.0/(1.0/0.88+exp(-17.832+0.245*WT))
    RXNGAMMA(6) = 1.e+0_f8/(1.e+0_f8/0.80e+0_f8+1.e+0_f8/ &
                  (exp(29.2e+0_f8-0.4e+0_f8*WT )+0.11))

    ! Reaction 7. BrNO3 + HCl
    RXNGAMMA(7) = 0.9e+0_f8 ! JPL 10-06

    ! Reaction 8. HOCl + HCl
    IF (HClOK) THEN
       c_HOCl    =  MOLEC_SPEED(T,52.46e+0_fp) ! input variable, fp
       SHOCl     =  0.0776e+0_f8+59.18e+0_f8/T
       HHOCl     =  1.91e-6_f8*EXP(5862.4e+0_f8/T)*EXP(-SHOCl*M_H2SO4)
       DHOCl     =  6.4e-8_f8*T/H
       kHOCl_HCl =  1.25e+9_f8*aH*DHOCl*MHCl
       lHOCl     =  SQRT(DHOCl/kHOCl_HCl)
       if (lHOCl.gt.(1.e+5_f8*arad)) then
          ! Limiting rate
          fHOCl  = arad/(3.e+0_f8*lHOCl)
       else
          fHOCl  =  1.e+0_f8/tanh(arad/lHOCl)- lHOCl/arad
       endif
       GHOClrxn  =  4.e+0_f8*HHOCl*0.082e+0_f8*T* &
                    sqrt(DHOCl*kHOCl_HCl)/c_HOCl
       !IF (fHOCl.eq.0.) THEN
       ! Catch for zero (ckeller, 12/29/17)
       IF ((fHOCl*GHOClrxn*FHCl).eq.0.) THEN
          gHOCl_HCl =  TINY(1e+0_f8)
       ELSE
          gHOCl_HCl =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/(fHOCl*GHOClrxn*FHCl))
       ENDIF
    ELSE
       gHOCl_HCl = TINY(1e+0_f8)
    ENDIF

    RXNGAMMA(8) = gHOCl_HCl

    ! Reaction 9. HOCl + HBr
    ! Not yet implemented for STS; JPL 10-06 suggests complex
    ! relationship, not yet sufficiently well understood or
    ! parameterized for the purposes of simulation. Ignore for now
    RXNGAMMA(9) = TINY(1e+0_f8)

    ! Reaction 10. HOBr + HCl
    IF ((HClOK).and.(HOBrOK)) THEN
       c_HOBr    =  MOLEC_SPEED(T,96.91e+0_fp) ! input variable, fp
       SHOBr     =  0.0776e+0_f8+59.18e+0_f8/T
       !HHOBr     =  30.D0
       HHOBR     = exp(-9.86e+0_f8+5427.e+0_f8/T)
       DHOBr     =  1.E-8_f8
       kII       = exp(154.e+0_f8-1.63e+0_f8*WT)*exp(-(3.85e+4_f8- &
                   478.e+0_f8*WT)/T)
       k_dl      = 7.5e+14_f8*(DHOBr*arad*1.e7_f8)
       IF (kII.gt.k_dl) kII=k_dl
       kHOBr_HCl =  kII*HHOBr*HOBrPP
       !IF (kHOBr_HCl.eq.0.) THEN   ! catch for zero (mdy, 04/15)
       !  kHOBr_HCl = TINY(1e+0_fp)
       !ENDIF
       GHOBrrxn  =  4.e+0_f8*HHCl*0.082e+0_f8*T*sqrt(DHOBr*kHOBr_HCl)/c_HOBr
       lHOBr     =  sqrt(DHOBr/kHOBr_HCl)
       if (lHOBr.gt.(1.e+3_f8*arad)) then
          ! Limiting rate
          fHOBr = arad/(3.e+0_f8*lHOBr)
       else
          fHOBr     =  1.e+0_f8/tanh(arad/lHOBr)- lHOBr/arad
       endif
       ! catch for zeros (ckeller, 12/29/17)
       IF ((fHOBr*GHOBrrxn).eq.0.) THEN
          gHOBr_HCl = TINY(1e+0_f8)
          !gHOBr_HCl =  1.e+0_fp/(1.e+0_fp+1.e+0_fp/TINY(1e+0_fp))
       ELSE
          gHOBr_HCl =  1.e+0_f8/(1.e+0_f8+1.e+0_f8/(fHOBr*GHOBrrxn))
       ENDIF
       RXNGAMMA(10) = gHOBr_HCl
    ELSE
       RXNGAMMA(10) = TINY(1e+0_f8)
    ENDIF

    ! Reaction 11. HOBr + HBr
    ! Data from JPL limited; ignore for now
    RXNGAMMA(11) = TINY(1e+0_f8)

    ! SDE 2013-10-18: DEBUG
    DO I=1,11
       IF (IT_IS_NAN(RXNGAMMA(I))) THEN
          WRITE(DBGMSG,'(a,I2)') 'RXNGAMMA NaN: ', I
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'NDENS: ', NDENS
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'T: ', T
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'P: ', P
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'WT_FRC: ', WT_FRC
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'H2OSUM: ', H2OSUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'HClSUM: ', HClSUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'HBrSUM: ', HBrSUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'HOBrSUM: ', HOBrSUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'ClNO3SUM: ', ClNO3SUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'BrNO3SUM: ', BrNO3SUM
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'RHO: ', RHO
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          WRITE(DBGMSG,'(a,E10.4)') 'ARAD: ', ARAD
          CALL DEBUG_MSG( TRIM(DBGMSG) )
          CALL ERROR_STOP('BAD GAMMA','UCX_mod')
       ENDIF
    ENDDO

  END SUBROUTINE CALC_SLA_GAMMA
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: molec_speed
!
! !DESCRIPTION: Function MOLEC\_SPEED calculates the mean velocity of gas
!  phase particles based on temperature and molecular mass.
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION MOLEC_SPEED(T,MOLMASS)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)          :: T       ! Temperature (K)
    REAL(fp), INTENT(IN)          :: MOLMASS ! Molecular mass (g/mol)
!
! !REVISION HISTORY:
!  10 Oct 2012 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! MOLEC_SPEED begins here!
    !=================================================================

    MOLEC_SPEED=SQRT(8.0e+0_fp*RSTARG*1e+7_fp*T/(PI*MOLMASS))

  END FUNCTION MOLEC_SPEED
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_h2o_trac
!
! !DESCRIPTION: Subroutine SET\_H2O\_TRAC sets the H2O species throughout
!  the selected domain (either troposphere only or the full grid).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_H2O_TRAC( SETSTRAT,  Input_Opt, State_Chm, State_Grid, &
                           State_Met, RC )
!
! !USES:
!
    USE Calc_Met_Mod,       ONLY : AIRQNT
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState, Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: SETSTRAT    ! Set strat H2O?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: I,J,L, LEVCPT
    INTEGER           :: TPLEV
    LOGICAL           :: READ_SPHU
    REAL(fp)          :: SPHU_kgkg, H2OVV_moist, Ev_mid
    REAL(fp)          :: Esat, EsatA, EsatB, EsatC, EsatD
    CHARACTER(LEN=255):: MSG, LOC

    ! Local variables for quantities from Input_Opt
    LOGICAL           :: LACTIVEH2O

    ! Empirical parameters for water vapor saturation pressure
    ! (Source: Nordquist, 1973. "Numerical Approximiations of
    !  Selected Meteorological Parameters Related to Cloud Physics"
    !  Text quality clarifications from Stipanuk, 1973. "Algorithms
    !  for Generating a Skew-T, Log P Diagram and Computing Selected
    !  Meteorological Quantities")
    REAL(fp), PARAMETER   :: ESATP1  = 2.3832241e+1_fp
    REAL(fp), PARAMETER   :: ESATP2  = -5.02808e+0_fp
    REAL(fp), PARAMETER   :: ESATP3  = 8.1328e-3_fp
    REAL(fp), PARAMETER   :: ESATP4  = 3.49149e+0_fp
    REAL(fp), PARAMETER   :: ESATP5  = -1.3028844e+3_fp
    REAL(fp), PARAMETER   :: ESATP6  = -1.3816e-7_fp
    REAL(fp), PARAMETER   :: ESATP7  = 1.1344e+1_fp
    REAL(fp), PARAMETER   :: ESATP8  = -3.03998e-2_fp
    REAL(fp), PARAMETER   :: ESATP9  = -2.949076e+3_fp

    ! Pointers
    REAL(fp), POINTER     :: Spc (:,:,:,:)

    ! If H2O in strat. is not prescribed, set H2O up to tropopause
    ! level plus a given offset (to avoid H2O leaking into
    ! stratosphere).
    INTEGER, PARAMETER    :: TROPPLEV_OFFSET = 0

    !=================================================================
    ! SET_H2O_TRAC begins here
    !=================================================================

    ! Assume success
    RC  = GC_SUCCESS

    ! Copy fields from INPUT_OPT
    LACTIVEH2O = Input_Opt%LACTIVEH2O

    ! Check that species concentration units are as expected
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'UCX_MOD: SET_H2O_TRAC'
       CALL GC_Error( TRIM(MSG), RC, TRIM(LOC) )
    ENDIF

    ! Initialize GEOS-Chem species array [kg/kg dry]
    Spc => State_Chm%Species

    ! Error trap: make sure id_H2O is defined. There are instances
    ! in GEOS-5 where we want to call this routine even with UCX
    ! turned off, in which case id_H2O might be undefined.
    ! (ckeller, 3/13/17)
    IF ( id_H2O <= 0 ) THEN
       id_H2O = Ind_('H2O')
    ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, TPLEV, READ_SPHU, LEVCPT) &
    !$OMP PRIVATE( SPHU_kgkg, H2OVV_moist, Ev_mid   ) &
    !$OMP PRIVATE( Esat, EsatA, EsatB, EsatC, EsatD )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Tropopause level for this grid box
       TPLEV = State_Met%TropLev(I,J)

       ! Between 30S-30N, find the level of cold point tropopause,
       ! the altitude with the coldest temperature below 10 hPa
       IF ( ABS(State_Grid%YMid(I,J)) <= 30 ) THEN
          ! Level with minimum temperature,
          ! use MASK to screen for >10 hPa
          LEVCPT = MINLOC( State_Met%T(I,J,:), DIM=1, &
                           MASK=(State_Met%PMID(I,J,:) >= 10) )
       ELSE
          LEVCPT = -1
       ENDIF

       DO L = 1, State_Grid%NZ

          ! Mark as tropospheric air everything that is at current
          ! tropopause or offset levels above. This is to avoid H2O
          ! leaks into the stratosphere due to jumps in the model
          ! tropopause. The offset is 0 by default (ckeller, 11/12/17)
          ! Logical flag determines if we are in a region where H2O from
          ! met files should overwrite the value from the chemical solver
          ! C. Holmes)
          READ_SPHU = ( ( (L-TPLEV) < TROPPLEV_OFFSET ) .or. &
                        SETSTRAT .or. ( .not. LACTIVEH2O ) .or. &
                        (L < LEVCPT)  )

          IF ( READ_SPHU ) THEN

             ! Calculate specific humidity [g H2O/kg total air] as
             ! [kg/kg]
             SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.e-3_fp

             ! Set H2O species concentration [kg/kg dry] using SPHU
             ! [kg/kg]
             ! without using box mass (ewl, 7/18/16)
             Spc(I,J,L,id_H2O) = SPHU_kgkg / ( 1.0e+0_fp - SPHU_kgkg )

          ELSE

             ! Calculate specific humidity in [kg H2O / kg total air]
             ! using transported H2O [kg/kg dry] (ewl, 7/18/16)
             SPHU_kgkg =  Spc(I,J,L,id_H2O) / &
                          ( 1.0e+0_fp + Spc(I,J,L,id_H2O) )

             ! Set State_Met specific humidity [g/kg]
             State_Met%SPHU(I,J,L) = SPHU_kgkg * 1.e+3_fp

             ! Calculate water vapor saturation pressure [hPa] from T
             ! (see local variables above for method reference)
             EsatA = ESATP1 + ESATP2 * log10( State_Met%T(I,J,L) )
             EsatB = ESATP3 * 10**( ESATP4+ESATP5/State_Met%T(I,J,L) )
             EsatC = ESATP6 * 10**( ESATP7+ESATP8*State_Met%T(I,J,L) )
             EsatD = ESATP9 / State_Met%T(I,J,L)
             Esat = 10**( EsatA + EsatB + EsatC + EsatD )

             ! Calculate mol water vapor per mol moist air from SPHU
             ! Note that SPHU must be converted to [kg/kg] to use 1-SPHU
             ! as mass dry air / mass moist air
             H2OVV_moist = SPHU_kgkg * AIRMW &
                       / ( SPHU_kgkg * AIRMW + (1 - SPHU_kgkg) * H2OMW )

             ! Calculate water vapor partial pressure at grid box center [hPa]
             ! Note that grid box center is defined as the arithmetic average
             ! of the grid box pressure edges not the vertical mid-point.
             Ev_mid = State_Met%PMID(I,J,L) * H2OVV_MOIST

             ! Set State_Met relative humidity [%]
             State_Met%RH(I,J,L) = ( Ev_mid / Esat ) * 100e+0_fp

          ENDIF
       ENDDO ! L
    ENDDO ! I
    ENDDO ! J
    !$OMP END PARALLEL DO

    ! If humidity was updated, update all moist-dependent air quantities
    ! and species mixing ratio with the new moisture content (ewl, 4/29/15)
    IF ( LActiveH2O .and. ( .not. SetStrat ) ) THEN
       CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, &
                    RC, Update_Mixing_Ratio=.TRUE. )
    ENDIF

    ! Free pointer
    NULLIFY( Spc )

  END SUBROUTINE SET_H2O_TRAC
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ucx_h2so4phot
!
! !DESCRIPTION: Subroutine UCX\_H2SO4PHOT propagates the calculated H2SO4
!  photolysis (J) rate at the top of the chemistry grid through to the top
!  of the transport grid, approximating H2SO4 photolysis in the mesosphere.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE UCX_H2SO4PHOT( Input_Opt, State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE CMN_FJX_MOD,        ONLY : ZPJ
    USE FAST_JX_MOD,        ONLY : RXN_H2SO4
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  17 Aug 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, J, L, IJWINDOW
    REAL(fp)              :: GMU,SO4_IN,PHOTDELTA,DTCHEM
    REAL(fp)              :: SO2_MW_G, SO4_MW_g, SO4_DELTA
    LOGICAL               :: DAYCOLUMN
    INTEGER               :: LMINPHOT
    REAL(fp)              :: RELWT
    LOGICAL,SAVE          :: FIRST=.TRUE.
    INTEGER               :: ICS

    ! Pointers
    REAL(fp), POINTER     :: Spc (:,:,:,:)

    !=================================================================
    ! UCX_H2SO4PHOT begins here!
    !=================================================================

    ! Copy fields from species database
    SO2_MW_G = State_Chm%SpcData(id_SO2)%Info%emMW_g ! g/mol
    SO4_MW_G = State_Chm%SpcData(id_SO4)%Info%emMW_g ! g/mol
    RELWT    = SO2_MW_G / SO4_MW_G

    ! Initialize GEOS-Chem species array [kg]
    Spc => State_Chm%Species

    ! Allow for the possibility of variable timestep
    DTCHEM = GET_TS_CHEM()

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, DAYCOLUMN, GMU ) &
    !$OMP PRIVATE( LMINPHOT, PHOTDELTA, SO4_IN, SO4_DELTA )
    DO J=1,State_Grid%NY
    DO I=1,State_Grid%NX

       GMU       = State_Met%SUNCOSmid(I,J)
       DAYCOLUMN = (GMU.gt.0e+0_fp)

       IF (DAYCOLUMN) THEN
          LMINPHOT  = State_Met%ChemGridLev(I,J)

          ! Retrieve photolysis rate as a fraction of gaseous SO4
          PHOTDELTA = ZPJ(LMINPHOT,RXN_H2SO4,I,J) * DTCHEM
          PHOTDELTA = MIN(1.e+0_fp,PHOTDELTA)

          DO L=LMINPHOT+1,State_Grid%NZ
             ! Apply photolysis to SO4
             ! First retrieve gaseous fraction
             SO4_IN = Spc(I,J,L,id_SO4)*SO4_PHOTFRAC(I,J,L)
             SO4_DELTA = PHOTDELTA*SO4_IN
             ! Remove from SO4
             Spc(I,J,L,id_SO4) = Spc(I,J,L,id_SO4) - SO4_DELTA
             ! Add to SO2. Note change in molar mass
             Spc(I,J,L,id_SO2) = Spc(I,J,L,id_SO2) + (SO4_DELTA*RELWT)
          ENDDO

       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    NULLIFY( Spc )

  END SUBROUTINE UCX_H2SO4PHOT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: noxcoeff_init
!
! !DESCRIPTION: Subroutine NOXCOEFF\_INIT initializes the NOX 2D interpolation
! values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NOXCOEFF_INIT( Input_Opt, State_Grid )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE FILE_MOD,           ONLY : IoError
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN) :: Input_Opt          ! Input options
    TYPE(GrdState),   INTENT(IN) :: State_Grid         ! Grid State object
!
! !REVISION HISTORY:
!  05 Dec 2014 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug, FileExists
    INTEGER            :: I, AS, IOS
    INTEGER            :: IMON, ITRAC, ILEV
    INTEGER            :: IU_FILE

    ! Strings
    CHARACTER(LEN=255) :: NOX_FILE
    CHARACTER(LEN=255) :: TARG_TRAC
    CHARACTER(LEN=255) :: DBGMSG
    CHARACTER(LEN=255) :: FileMsg

    !=================================================================
    ! NOXCOEFF_INIT begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    !=================================================================
    ! In dry-run mode, print file paths to dryrun log and exit.
    !=================================================================
    IF ( Input_Opt%DryRun ) THEN

       ! Loop over months and tracers
       DO IMON  = 1,12
          DO ITRAC = 1,6

             ! Pick tracer name
             SELECT CASE (ITRAC)
             CASE ( 1 )
                TARG_TRAC = 'O'
             CASE ( 2 )
                TARG_TRAC = 'O1D'
             CASE ( 3 )
                TARG_TRAC = 'JNO'
             CASE ( 4 )
                TARG_TRAC = 'JNO2'
             CASE ( 5 )
                TARG_TRAC = 'JNO3'
             CASE ( 6 )
                TARG_TRAC = 'JN2O'
             END SELECT

             ! Create file name for each tracer and month
             WRITE( NOx_File, '(a,a,a,I0.2,a)' ) &
                  TRIM(NOON_FILE_ROOT), TRIM(TARG_TRAC), '_', IMON, '.dat'

             ! Test if the file exists
             INQUIRE( FILE=TRIM( NOx_File ), EXIST=FileExists )

             ! Test if the file exists and define an output string
             IF ( FileExists ) THEN
                FileMsg = 'UCX (SFCMR_READ): Opening'
             ELSE
                FileMsg = 'UCX (SFCMR_READ): REQUIRED FILE NOT FOUND'
             ENDIF

             ! Write to stdout for both regular and dry-run simulations
             IF ( Input_Opt%amIRoot ) THEN
                WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NOx_File )
300             FORMAT( a, ' ', a )
             ENDIF
          ENDDO
       ENDDO

       ! exit without doing any computation (dry-run only)
       RETURN
    ENDIF

    !=================================================================
    ! If this is a regular simulation, then read data from files
    !=================================================================

    ! Number of latitude levels of NOXCOEFF array. Data is only
    ! available for 2x25 and 4x5. Use 2x25 for any different grid
    ! and map NOx coeffs onto simulation grid when calling
    ! GET_NOXCOEFF.
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' .or. &
         TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       JJNOXCOEFF = State_Grid%NY
    ELSE
       JJNOXCOEFF = 91
    ENDIF

    ! Fill NOx latitudes
    ALLOCATE(NOXLAT(JJNOXCOEFF+1), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXLAT' )

    ! Fill manually
    NOXLAT(2) = -89.0e+0_fp
    DO I = 3,JJNOXCOEFF
       NOXLAT(I) = NOXLAT(I-1) + 2e+0_fp
    ENDDO
    ! Overshoot to make sure that a latitude of 90.0 will be properl
    ! matched onto JJNOXCOEFF.
    NOXLAT(JJNOXCOEFF+1) = 90.5e+0_fp

    ! Initialize the NOXCOEFF array. This array holds monthly NOx
    ! coefficients on 51 levels and for 6 species.
    ALLOCATE(NOXCOEFF(JJNOXCOEFF,UCX_NLEVS,6,12), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOXCOEFF' )
    NOXCOEFF = 0.0e+0_fp

    ! Fill array
    DO IMON  = 1,12
    DO ITRAC = 1,6
       SELECT CASE (ITRAC)
       CASE ( 1 )
          TARG_TRAC = 'O'
       CASE ( 2 )
          TARG_TRAC = 'O1D'
       CASE ( 3 )
          TARG_TRAC = 'JNO'
       CASE ( 4 )
          TARG_TRAC = 'JNO2'
       CASE ( 5 )
          TARG_TRAC = 'JNO3'
       CASE ( 6 )
          TARG_TRAC = 'JN2O'
       END SELECT
       WRITE(NOX_FILE,'(a,a,a,I0.2,a)') TRIM(NOON_FILE_ROOT), &
            TRIM(TARG_TRAC), '_', IMON, '.dat'

       ! Get a free LUN
       IU_FILE = findFreeLUN()

       IOS = 1

       ! Test if the file exists
       INQUIRE( FILE=TRIM( NOx_File ), EXIST=FileExists )

       ! Define an output string
       IF ( FileExists ) THEN
          FileMsg = 'UCX (SFCMR_READ): Opening'
       ELSE
          FileMsg = 'UCX (SFCMR_READ): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Print the name of the file being read to stdout (bmy, 10/23/19)
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NOx_File )
       ENDIF

       OPEN( IU_FILE,FILE=TRIM(NOX_FILE),STATUS='OLD',IOSTAT=IOS)
       IF ( IOS /= 0 ) THEN
          WRITE(6,*) 'UCX_MOD: Could not read ', TRIM(NOX_FILE)
          CALL IOERROR( IOS, IU_FILE,'UCX_MOD:NOXCOEFF_INIT')
       ENDIF

       IF ( prtDebug ) THEN
          WRITE(DBGMSG,'(a,a)') ' ### UCX: Reading ', TRIM( NOX_FILE )
          CALL DEBUG_MSG( TRIM(DBGMSG) )
       ENDIF

       ! Read in data
       DO ILEV = 1,UCX_NLEVS

          IF ( TRIM(State_Grid%GridRes) =='4.0x5.0' ) THEN
             READ(IU_FILE, 110, IOSTAT=IOS ) NOXCOEFF(:,ILEV,ITRAC,IMON)
110          FORMAT(46E10.3)
          ELSE
             ! Use 2x25 as default
             READ(IU_FILE, 120, IOSTAT=IOS ) NOXCOEFF(:,ILEV,ITRAC,IMON)
120          FORMAT(91E10.3)
          ENDIF
          IF ( IOS /= 0 ) THEN
             WRITE(6,'(a,a,I4,a,1x,a)') 'UCX_MOD: Error reading ', &
                  'line ', ILEV, ' in file ', TRIM( NOX_FILE )
             CALL IOERROR( IOS, IU_FILE,'UCX_MOD:NOXCOEFF_INIT')
          ENDIF
       ENDDO

       CLOSE(IU_FILE)

    ENDDO !ITRAC
    ENDDO !IMON

  END SUBROUTINE NOXCOEFF_INIT
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jjnox
!
! !DESCRIPTION: Subroutine GET\_JJNOX maps grid box at location IISIM, JJSIM of
! the simulation grid onto the latitude grid of the NOXCOEFF array. JJNOX can
! differ from JJSIM if it's not a 4x5 or 2x25 simulation.
!\\
!\\
! This routine simply returns the index of the NOx latitude vector that covers
! the latitude value of interest. No grid box weighting, etc. is performed.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_JJNOX( IISIM, JJSIM, State_Grid ) RESULT ( JJNOX )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: IISIM    ! Latitude index on simulation grid
    INTEGER,          INTENT(IN) :: JJSIM    ! Latitude index on simulation grid
    TYPE(GrdState),   INTENT(IN) :: State_Grid ! Grid State object
!
! !OUTPUT VARIABLES:
!
    INTEGER                      :: JJNOX    ! Latitude index on NOXCOEFF grid
!
! !REVISION HISTORY:
!  05 Dec 2014 - C. Keller - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I
    REAL(fp)              :: LAT

    !=================================================================
    ! GET_JJNOX begins here!
    !=================================================================

    ! Nothing to do for 'standard' grids
    IF ( .not. State_Grid%NestedGrid ) THEN
       JJNOX = JJSIM
    ELSE

       ! Init
       JJNOX = -1

       ! Get latitude in degrees north on simulation grid
       LAT = State_Grid%YMid( IISIM, JJSIM )

       ! Loop over all latitudes of the NOx grid until we reach the grid
       ! box where the simulation latitude sits in.
       DO I = 1,JJNOXCOEFF
          IF ( LAT < NOXLAT(I+1) ) THEN
             JJNOX = I
             EXIT
          ENDIF
       ENDDO

    ENDIF

  END FUNCTION GET_JJNOX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ucx
!
! !DESCRIPTION: Subroutine INIT\_UCX initializes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NSTRATAER
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE ERROR_MOD,      ONLY : IS_SAFE_DIV
    USE ERROR_MOD,      ONLY : ERROR_STOP
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : PI_180
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN) :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  04 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: N, AS
    INTEGER            :: JIN, JOUT
    REAL(fp)           :: JMIN_IN, JMAX_IN, JDIF_IN
    REAL(fp)           :: JMIN_OUT,JMAX_OUT,JDIF_OUT
    REAL(fp)           :: JMIN_TMP,JMAX_TMP,JDIF_TMP
    REAL(fp)           :: JRATIO
    REAL(fp)           :: DEG_SUM

    ! Local variables for quantities from Input_Opt
    LOGICAL            :: prtDebug
    LOGICAL            :: LUCX

    ! Strings
    CHARACTER(LEN=255) :: DBGMSG, GRIDSPEC, FileMsg, FileName

    ! Return code
    INTEGER :: RC

    !=================================================================
    ! INIT_UCX begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    LUCX     = Input_Opt%LUCX

    ! Initialize species ID flags
    id_BCPI  = Ind_('BCPI'      )
    id_BrNO3 = Ind_('BrNO3'     )
    id_ClNO3 = Ind_('ClNO3'     )
    id_H2O   = Ind_('H2O'       )
    id_HBr   = Ind_('HBr'       )
    id_HCl   = Ind_('HCl'       )
    id_HNO3  = Ind_('HNO3'      )
    id_HOBr  = Ind_('HOBr'      )
    id_HOCl  = Ind_('HOCl'      )
    id_N2O   = Ind_('N2O'       )
    id_NIT   = Ind_('NIT'       )
    id_N     = Ind_('N'         )
    id_NO    = Ind_('NO'        )
    id_NO2   = Ind_('NO2'       )
    id_NO3   = Ind_('NO3'       )
    id_O3    = Ind_('O3'        )
    id_SO2   = Ind_('SO2'       )
    id_SO4   = Ind_('SO4'       )

    ! Print info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,'(a)') REPEAT( '=', 79 )
       WRITE( 6,'(a)') 'U N I F I E D   C H E M I S T R Y'
       WRITE( 6,'(a)') 'Routines written by SEBASTIAN D. EASTHAM'
       WRITE( 6,'(a)') REPEAT( '=', 79 )
    ENDIF

    ! --------------------------------------------------------------
    ! Input data sources
    ! --------------------------------------------------------------

    ! Determine folder paths from root folder
    ! Regridding of netCDF input data is currently not supported.
    ! NOTE: netCDF inputs not read by GCHP, GEOS-5, WRF-GC
    IF (UCXNETCDF) THEN

       IF ( TRIM(State_Grid%GridRes) /= '4.0x5.0'   .and. &
            TRIM(State_Grid%GridRes) /= '2.0x2.5' ) THEN
          DBGMSG = 'NetCDF NOx coeffs only preprocessed for 2x25' // &
                   ' and 4x5 - try setting the module variable'   // &
                   ' UCXNETCDF=.FALSE.'
          CALL ERROR_STOP( DBGMSG, 'INIT_UCX (UCX_mod.F90)!' )
          RETURN
       ENDIF

       WRITE( NOON_FILE_ROOT,'(a,a)') &
            TRIM(Input_Opt%CHEM_INPUTS_DIR), 'UCX_201403/Init2D/Noontime.nc'

    ! For ASCII input, use 2x25 grid for all other grids than 4x5.
    ! This is ok for the NOx coeffs which can be regridded on the fly
    ! from 2x25 onto any other grid. This won't work for the 2D
    ! boundary conditions, but those have been checked in the logical
    ! check above (USE2DDATA).
    ELSE

       IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN
          GRIDSPEC = 'Grid4x5/InitCFC_'
       ELSE
          GRIDSPEC = 'Grid2x25/InitCFC_'
       ENDIF
       WRITE(   NOON_FILE_ROOT,'(a,a,a)') TRIM(Input_Opt%CHEM_INPUTS_DIR), &
#ifdef MODEL_GEOS
            'UCX_201710/NoonTime/', TRIM(GRIDSPEC)
#else
            'UCX_201403/NoonTime/', TRIM(GRIDSPEC)
#endif

    ENDIF

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and return.
    !=================================================================
    IF ( Input_Opt%DryRun ) THEN
         
       ! Test if the file exists and define an output string
       FileName = Noon_File_Root
       INQUIRE( FILE=TRIM( FileName ), EXIST=FileExists )
       IF ( FileExists ) THEN
          FileMsg = 'UCX (SFCMR_READ): Opening'
       ELSE
          FileMsg = 'UCX (SFCMR_READ): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write to stdout for the dry-run simulation
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( FileName )
300       FORMAT( a, ' ', a )
       ENDIF

       ! Exit without doing any computation (dry-run only)
       RETURN
    ENDIF

    !=================================================================
    ! For regular simulations, read data from files
    !=================================================================

    ! Write the status of Noon_File_Root to stdout
    INQUIRE( FILE=TRIM( Noon_File_Root ), EXIST=FileExists )
    IF ( FileExists ) THEN
       FileMsg = 'UCX (SFCMR_READ): Opening'
    ELSE
       FileMsg = 'UCX (SFCMR_READ): REQUIRED FILE NOT FOUND'
    ENDIF
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( Noon_File_Root )
    ENDIF

    IF ( prtDebug ) THEN
       WRITE(DBGMSG,'(a,a)') '### UCX: Reading O1D/O3P from ', &
            TRIM(NOON_FILE_ROOT)
       CALL DEBUG_MSG( TRIM(DBGMSG) )
    ENDIF

    ! Allocate arrays of input pressure levels and lat edges
    ALLOCATE( UCX_PLEVS( UCX_NLEVS ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_PLEVS' )

    ALLOCATE( UCX_LATS( UCX_NLAT+1 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_LATS' )

    ! Set input pressure levels (hPa)
    UCX_PLEVS = (/ 0.2200e+00_fp, 0.2600e+00_fp, &
                   0.3100e+00_fp, 0.3600e+00_fp, &
                   0.4300e+00_fp, 0.5100e+00_fp, &
                   0.6000e+00_fp, 0.7100e+00_fp, &
                   0.8400e+00_fp, 0.9900e+00_fp, &
                   1.1700e+00_fp, 1.3800e+00_fp, &
                   1.6300e+00_fp, 1.9300e+00_fp, &
                   2.2800e+00_fp, 2.6900e+00_fp, &
                   3.1800e+00_fp, 3.7600e+00_fp, &
                   4.4400e+00_fp, 5.2500e+00_fp, &
                   6.2000e+00_fp, 7.3200e+00_fp, &
                   8.6500e+00_fp, 1.0220e+01_fp, &
                   1.2070e+01_fp, 1.4260e+01_fp, &
                   1.6850e+01_fp, 1.9910e+01_fp, &
                   2.3520e+01_fp, 2.7780e+01_fp, &
                   3.2820e+01_fp, 3.8770e+01_fp, &
                   4.5810e+01_fp, 5.4110e+01_fp, &
                   6.3930e+01_fp, 7.5220e+01_fp, &
                   8.9220e+01_fp, 1.0540e+02_fp, &
                   1.2451e+02_fp, 1.4710e+02_fp, &
                   1.7377e+02_fp, 2.0529e+02_fp, &
                   2.4252e+02_fp, 2.8650e+02_fp, &
                   3.3847e+02_fp, 3.9985e+02_fp, &
                   4.7237e+02_fp, 5.5804e+02_fp, &
                   6.5924e+02_fp, 7.7880e+02_fp, &
                   9.2004e+02_fp /)

    ! Set input latitude edges (degrees)
    UCX_LATS(1)  = -90.0e+0_fp
    UCX_LATS(2)  = (-90.0e+0_fp) + (9.5e+0_fp/2.0e+0_fp)
    UCX_LATS(20) = 90.0e+0_fp
    DO N=2,18
       UCX_LATS(N+1) = UCX_LATS(N) + 9.5e+0_fp
    ENDDO

    ! Calculate conversion factors for SLA
    ! Factor to convert volume (m3 SLA/m3 air) to
    ! surface area density (cm2 SLA/cm3 air)
    SLA_VA = (8.406e-8_fp)*(10.e+0_fp**(12.e+0_fp*0.751e+0_fp))

    ! Factor to convert effective radius to
    ! liquid radius (unitless)
    SLA_RR = EXP(-0.173e+0_fp)

    ! Factor to convert volume (m3/m3) to effective
    ! radius (m)
    SLA_VR = (0.357e-6_fp)*(10.e+0_fp**(12.e+0_fp*0.249))

    ! Initialize NOx coefficient arrays
    ALLOCATE( NOX_O( State_Grid%NX, State_Grid%NY, State_Grid%NZ, 2 ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_O' )
    NOX_O = 0e+0_fp

    ALLOCATE( NOX_J( State_Grid%NX, State_Grid%NY, State_Grid%NZ, 4 ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX_J' )
    NOX_J = 0e+0_fp

    ! Initialize PSC variables
    ALLOCATE( RAD_AER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                       NSTRATAER ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'RAD_AER' )
    RAD_AER = 0e+0_fp

    ALLOCATE( KG_AER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                      NSTRATAER ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'KG_AER' )
    KG_AER = 0e+0_fp

    ALLOCATE( SAD_AER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                       NSTRATAER ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SAD_AER' )
    SAD_AER = 0e+0_fp

    ALLOCATE( NDENS_AER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                         NSTRATAER ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NDENS_AER' )
    NDENS_AER = 0e+0_fp

    ALLOCATE( RHO_AER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                       NSTRATAER ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'RHO_AER' )
    RHO_AER = 0e+0_fp

    ! Mass fraction of species contained in liquid aerosol
    ! Indices: 1 - SO4
    !          2 - HNO3
    !          3 - HCl
    !          4 - HOCl
    !          5 - HBr
    !          6 - HOBr
    !          7 - H2O
    ALLOCATE( AERFRAC( State_Grid%NX, State_Grid%NY, State_Grid%NZ,7), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRAC' )
    AERFRAC = 0e+0_fp

    ALLOCATE( AERFRACIND( 7 ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'AERFRACIND' )
    AERFRACIND(1) = id_SO4
    AERFRACIND(2) = id_HNO3
    AERFRACIND(3) = id_HCl
    AERFRACIND(4) = id_HOCl
    AERFRACIND(5) = id_HBr
    AERFRACIND(6) = id_HOBr
    AERFRACIND(7) = id_H2O

    ! H2SO4 photolysis rate at the top of the chemgrid
    ALLOCATE( SO4_TOPPHOT( State_Grid%NX,State_Grid%NY ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_TOPPHOT' )
    SO4_TOPPHOT = 0.e+0_fp

    ALLOCATE( UCX_REGRID( State_Grid%NY, UCX_NLAT ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'UCX_REGRID' )
    UCX_REGRID = 0e+0_fp

    ! Calculate the scaling matrix
    ! Note cosine (area-weighted)
    JMAX_OUT = State_Grid%YEdge(1,1)
    DO JOUT=1,State_Grid%NY
       JMIN_OUT = JMAX_OUT
       JMAX_OUT = State_Grid%YEdge(1,JOUT+1)
       JDIF_OUT = SIN( JMAX_OUT * PI_180 ) &
                - SIN( JMIN_OUT * PI_180 )
       DEG_SUM = 0e+0_fp
       DO JIN=1,UCX_NLAT
          JMIN_IN = UCX_LATS(JIN)
          JMAX_IN = UCX_LATS(JIN+1)
          IF ((JMAX_OUT.ge.JMIN_IN).and.(JMIN_OUT.le.JMAX_IN)) THEN
             JMAX_TMP = MIN(JMAX_IN,JMAX_OUT)
             JMIN_TMP = MAX(JMIN_IN,JMIN_OUT)
             JDIF_TMP = SIN( JMAX_TMP * PI_180 ) &
                      - SIN( JMIN_TMP * PI_180 )
             IF (IS_SAFE_DIV(JDIF_TMP,JDIF_OUT)) THEN
                JRATIO = JDIF_TMP/JDIF_OUT
                UCX_REGRID(JOUT,JIN) = JRATIO
                DEG_SUM = DEG_SUM + JRATIO
             ENDIF
          ENDIF
       ENDDO
       ! Normalize
       IF (DEG_SUM.gt.0e+0_fp) THEN
          DO JIN=1,UCX_NLAT
             UCX_REGRID(JOUT,JIN) = UCX_REGRID(JOUT,JIN)/DEG_SUM
          ENDDO
       ELSE
          UCX_REGRID(JOUT,:) = 0e+0_fp
       ENDIF

       !! Debug
       !IF ( prtDebug ) THEN
       !   WRITE(DBGMSG,'(a,I03,a,3(F6.2,x))') '### UCX: Exgrid: J-', &
       !     JOUT, '->',JMIN_OUT,JMAX_OUT,JDIF_OUT
       !   CALL DEBUG_MSG( TRIM(DBGMSG) )
       !   WRITE(DBGMSG,'(a,I03,a,F6.2)') '### UCX: Regrid: J-', &
       !     JOUT, '->',DEG_SUM
       !   CALL DEBUG_MSG( TRIM(DBGMSG) )
       !ENDIF
    ENDDO

    ! Initialize NOXCOEFF arrays
    IF ( .NOT. UCXNETCDF ) THEN
       CALL NOXCOEFF_INIT( Input_Opt, State_Grid )
    ENDIF

  END SUBROUTINE INIT_UCX
!EOC
!------------------------------------------------------------------------------
!               MIT Laboratory for Aviation and the Environment               !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_ucx
!
! !DESCRIPTION: Subroutine CLEANUP\_UCX deallocates module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_UCX ( )
!
! !REVISION HISTORY:
!  04 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: N

    !=================================================================
    ! CLEANUP_UCX begins here!
    !=================================================================

    IF ( ALLOCATED( RAD_AER    ) ) DEALLOCATE( RAD_AER    )
    IF ( ALLOCATED( SAD_AER    ) ) DEALLOCATE( SAD_AER    )
    IF ( ALLOCATED( KG_AER     ) ) DEALLOCATE( KG_AER     )
    IF ( ALLOCATED( RHO_AER    ) ) DEALLOCATE( RHO_AER    )
    IF ( ALLOCATED( NDENS_AER  ) ) DEALLOCATE( NDENS_AER  )
    IF ( ALLOCATED( AERFRAC    ) ) DEALLOCATE( AERFRAC    )
    IF ( ALLOCATED( AERFRACIND ) ) DEALLOCATE( AERFRACIND )
    IF ( ALLOCATED( UCX_REGRID ) ) DEALLOCATE( UCX_REGRID )
    IF ( ALLOCATED( UCX_PLEVS  ) ) DEALLOCATE( UCX_PLEVS  )
    IF ( ALLOCATED( UCX_LATS   ) ) DEALLOCATE( UCX_LATS   )
    IF ( ALLOCATED( NOX_O      ) ) DEALLOCATE( NOX_O      )
    IF ( ALLOCATED( NOX_J      ) ) DEALLOCATE( NOX_J      )
    IF ( ALLOCATED( SO4_TOPPHOT) ) DEALLOCATE( SO4_TOPPHOT)

    ! Cleanup the NOx coeff arrays
    IF ( ALLOCATED( NOXCOEFF ) ) DEALLOCATE( NOXCOEFF )
    IF ( ALLOCATED( NOXLAT   ) ) DEALLOCATE( NOXLAT   )

  END SUBROUTINE CLEANUP_UCX
!EOC
END MODULE UCX_MOD
