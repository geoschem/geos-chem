#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tomas_mod.F90
!
! !DESCRIPTION: Module TOMAS\_MOD contains variable specific to the TOMAS
!  aerosol microphysics simulation, e.g. number of species, number of size-bins,
!  mass per particle bin boundaries and arrays used inside the microphysics
!  subroutine. (win, 7/9/07)
!\\
!\\
! !INTERFACE:
!
MODULE TOMAS_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !REMARKS:
!  This module also contains what used to be in sizecode.COM header file
!  containing common blocks for TOMAS aerosol microphysics algorithm
!  originally implemented in GISS GCM-II' by Peter Adams. Below is the original
!  comment from sizecode.COM
!
!     This header file includes all the variables used by the
!     size-resolved aerosol microphysics code incorporated into
!     the GISS GCM II' by Peter Adams.  The microphysics algorithm
!     conserves aerosol number and mass using the schemes developed
!     by Graham Feingold and others.
!
!  References:
!  ============================================================================
!  Tzivion, S., Feingold, G., and Levin, Z., An Efficient
!   Numerical Solution to the Stochastic Collection Equation,
!   J. Atmos. Sci., 44, 3139-3149, 1987.
!  Feingold, G., Tzivion, S., and Levin, Z., Evolution of
!   Raindrop Spectra. Part I: Solution to the Stochastic
!   Collection/Breakup Equation Using the Method of Moments,
!   J. Atmos. Sci., 45, 3387-3399, 1988.
!  Tzivion, S., Feingold, G., and Levin, Z., The Evolution of
!   Raindrop Spectra. Part II: Collisional Collection/Breakup
!   and Evaporation in a Rainshaft, J. Atmos. Sci., 46, 3312-
!   3327, 1989.
!  Feingold, G., Levin, Z., and Tzivion, S., The Evolution of
!   Raindrop Spectra. Part III: Downdraft Generation in an
!   Axisymmetrical Rainshaft Model, J. Atmos. Sci., 48, 315-
!   330, 1991.
!
!  The algorithms described in these papers have been extended
!  to include multicomponent aerosols and modified for a moving
!  sectional approach.  Using this approach, the boundaries
!  between size bins are defined in terms of dry aerosol mass
!  such that the actual sizes of the sections move as water
!  is added to or lost from the aerosol.
!
!  All of the subroutines needed for this aerosol microphysics
!  algorithm use only their own internal variables or the ones
!  listed here.  GISS GCM II' variables are not used (a driver
!  subroutine performs the necessary swapping between the GCM
!  and the microphysics code).  The microphysics code is,
!  therefore, completely modular.
!
! !REVISION HISTORY:
!  09 Jul 2006 - W. Trivitayanurak - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  !========================================================================
  ! Module Variables:
  ! IBINS   : Number of size bins
  ! ICOMP   : Number of aerosol mass species + 1 for number
  ! Nk      : Aerosol number internal array
  ! Mk      : Aerosol mass internal array
  ! Gc      : Condensing gas
  ! Nkd     : Aerosol number diagnostic array
  ! Mkd     : Aerosol mass diagnostic array
  ! Gcd     : Condensing gas diagnostic array
  ! Xk      : Size bin boundary in dry mass per particle
  ! MOLWT   : Aerosol molecular weight
  ! SRTSO4  : ID of sulfate
  ! SRTNACL : ID of sea-salt
  ! SRTH2O  : ID of aerosol water
  ! SRTECOB : ID of hydrophobic EC
  ! SRTECIL : ID of hydrophilic EC
  ! SRTOCOB : ID of hydrophobic OC
  ! SRTOCIL : ID of hydrophilic OC
  ! dust??  : ID of internally mixed dust (future work)
  ! dust??  : ID of externally mixed dust (future work)
  ! BOXVOL  : Grid box volume  (cm3)
  ! BOXMASS : Grid box air mass (kg)
  ! TEMPTMS : Temperature (K) of each grid box
  ! PRES    : Pressure (Pa) in grid box
  ! RHTOMAS : Relative humidity (0-1)
  ! BINACT1 : Activated bin as a function of composition
  ! FRACTION: Activated fraction as a fcn of composition
  ! IDIAG   : Number of diagnostic tracer (NH4 and H2O)
  !========================================================================
!
! !DEFINED PARAMETERS:
!
  ! Note that there is a parameter declared in CMN_SIZE called TOMASBIN -->
  ! which defines how many size bins (win, 7/7/09)
#if   defined ( TOMAS12 )
  INTEGER, PARAMETER   :: IBINS = 12
#elif defined ( TOMAS15 )
  INTEGER, PARAMETER   :: IBINS = 15
#elif defined ( TOMAS40 )
  INTEGER, PARAMETER   :: IBINS = 40
#else
  INTEGER, PARAMETER   :: IBINS = 30
#endif

  INTEGER, PARAMETER   :: SRTSO4  = 1
  INTEGER, PARAMETER   :: SRTNACL = 2
  INTEGER, PARAMETER   :: SRTECIL = 3
  INTEGER, PARAMETER   :: SRTECOB = 4
  INTEGER, PARAMETER   :: SRTOCIL = 5
  INTEGER, PARAMETER   :: SRTOCOB = 6
  INTEGER, PARAMETER   :: SRTDUST = 7
  INTEGER, PARAMETER   :: SRTNH4  = 8
  INTEGER, PARAMETER   :: SRTH2O  = 9

  INTEGER, PARAMETER   :: ICOMPHARD = 9
!
! !PUBLIC DATA MEMBERS:
!
  ! Scalars
  INTEGER                       :: ICOMP,   IDIAG

  ! Arrays
  REAL(fp), TARGET, SAVE        :: Xk(IBINS+1)
  REAL*4,   save,   ALLOCATABLE :: MOLWT(:)
  INTEGER,          SAVE        :: BINACT1(101,101,101)
  REAL(fp),         SAVE        :: FRACTION1(101,101,101)
  INTEGER,          SAVE        :: BINACT2(101,101,101)
  REAL(fp),         SAVE        :: FRACTION2(101,101,101)

  REAL(fp) :: AVGMASS(IBINS)       ! Average mass per particle
                                   ! mid-range of size bin
                                   ! [kg/no.]
  REAL(fp) :: cosmic_ions(72,46,9) !careful, this is not GCHP safe!
                                   ! [ion pairs kg^-1 s^-1]

  REAL(fp),         SAVE        :: OCSCALE30(IBINS)
  REAL(fp),         SAVE        :: OCSCALE100(IBINS)
  REAL(fp),         SAVE        :: ECSCALE30(IBINS)
  REAL(fp),         SAVE        :: ECSCALE100(IBINS)

  INTEGER  :: bin_nuc = 1, tern_nuc = 1  ! Switches for nucleation type.
  INTEGER  :: act_nuc = 0 ! in BL
  INTEGER  :: ion_nuc = 0 ! 1 for modgil, 2 for Yu
  INTEGER  :: absall  = 1 ! 1 for soa absorb to all specnapari
                          ! nucleation tuned by factor of 1.0D-5

  REAL(fp) :: soaareafrac=1.0e+0_fp ! fraction of SOA that goes
                                    ! to SA rather than into mass

  INTEGER :: xSOA = 1     !Switch for xSOA. If set to one, emit 100
                          ! Tg/yr SOA correlated with anrtho CO
                          ! (JKodros 6/3/15)
                          ! switch to 0 for anthro-free
                          ! runs (Pengfei Liu 4/18/2018)
  INTEGER :: lowRH = 1    !This is to match AW more with AERONET (JKODROS 6/15)

  REAL(fp), ALLOCATABLE :: H2SO4_RATE(:,:,:) ! H2SO4 prod rate [kg s-1]

#if  defined( TOMAS12 ) || defined( TOMAS15 )
  !tomas12 or tomas15
  data OCSCALE30/ &
# if  defined( TOMAS15)
       0.0e+0_fp     , 0.0e+0_fp     , 0.0e+0_fp     , &
# endif
       1.1291E-03, 4.9302E-03, 1.2714E-02, 3.6431E-02, &
       1.0846E-01, 2.1994E-01, 2.7402E-01, 2.0750E-01, &
       9.5304E-02, 2.6504E-02, 1.2925E-02, 1.6069E-05/! use for fossil fuel (bimodal)

  data OCSCALE100/ &
# if  defined( TOMAS15)
       0.0e+0_fp     , 0.0e+0_fp     , 0.0e+0_fp     , &
# endif
       1.9827E-06, 3.9249E-05, 5.0202E-04, 4.1538E-03, &
       2.2253E-02, 7.7269E-02, 1.7402E-01, 2.5432E-01, &
       2.4126E-01, 1.4856E-01, 7.6641E-02, 9.8120E-04/! use for biomass burning

  data ECSCALE30/ &
# if  defined( TOMAS15)
       0.0e+0_fp     , 0.0e+0_fp     , 0.0e+0_fp     , &
# endif
       1.1291E-03, 4.9302E-03, 1.2714E-02, 3.6431E-02, &
       1.0846E-01, 2.1994E-01, 2.7402E-01, 2.0750E-01, &
       9.5304E-02, 2.6504E-02, 1.2925E-02, 1.6069E-05/! use for fossil fuel (bimodal)

  data ECSCALE100/ &
# if  defined( TOMAS15)
       0.0e+0_fp     , 0.0e+0_fp     , 0.0e+0_fp     , &
# endif
       1.9827E-06, 3.9249E-05, 5.0202E-04, 4.1538E-03, &
       2.2253E-02, 7.7269E-02, 1.7402E-01, 2.5432E-01, &
       2.4126E-01, 1.4856E-01, 7.6641E-02, 9.8120E-04/  ! use for biomass burning

#else
  !tomas30 or tomas40
  DATA OCSCALE30/  &     ! use for fossil fuel
# if  defined( TOMAS40)
       0.0     , 0.0     , 0.0     , 0.0     , 0.0     , &
       0.0     , 0.0     , 0.0     , 0.0     , 0.0     , &
# endif
       1.04E-03, 2.77E-03, 6.60E-03, 1.41E-02, 2.69E-02, &
       4.60E-02, 7.06E-02, 9.69E-02, 1.19E-01, 1.31E-01, &
       1.30E-01, 1.15E-01, 9.07E-02, 6.44E-02, 4.09E-02, &
       2.33E-02, 1.19E-02, 5.42E-03, 2.22E-03, 8.12E-04, &
       2.66E-04, 7.83E-05, 2.06E-05, 4.86E-06, 1.03E-06, &
       1.94E-07, 3.29E-08, 4.99E-09, 6.79E-10, 8.26E-11/

  DATA OCSCALE100/  &    ! use for biomass burning
# if  defined( TOMAS40)
       0.0        , 0.0        , 0.0        , 0.0        , 0.0       , &
       0.0        , 0.0        , 0.0        , 0.0        , 0.0       , &
# endif
       3.2224e-07 , 1.6605e-06 , 7.6565e-06 , 3.1592e-05 , 0.00011664, &
       0.00038538 , 0.0011394  , 0.0030144  , 0.0071362  , 0.015117  , &
       0.028657   , 0.048612   , 0.073789   , 0.10023    , 0.12182   , &
       0.1325     , 0.12895    , 0.11231    , 0.087525   , 0.061037  , &
       0.038089   , 0.02127    , 0.010628   , 0.0047523  , 0.0019015 , &
       0.00068081 , 0.00021813 , 6.2536e-05 , 1.6044e-05 , 3.6831e-06/

  DATA ECSCALE30/  &     ! use for fossil fuel
# if  defined( TOMAS40)
       0.0     , 0.0     , 0.0     , 0.0     , 0.0     , &
       0.0     , 0.0     , 0.0     , 0.0     , 0.0     , &
# endif
       1.04E-03, 2.77E-03, 6.60E-03, 1.41E-02, 2.69E-02, &
       4.60E-02, 7.06E-02, 9.69E-02, 1.19E-01, 1.31E-01, &
       1.30E-01, 1.15E-01, 9.07E-02, 6.44E-02, 4.09E-02, &
       2.33E-02, 1.19E-02, 5.42E-03, 2.22E-03, 8.12E-04, &
       2.66E-04, 7.83E-05, 2.06E-05, 4.86E-06, 1.03E-06, &
       1.94E-07, 3.29E-08, 4.99E-09, 6.79E-10, 8.26E-11/

  DATA ECSCALE100/  &    ! use for biomass burning
# if  defined( TOMAS40)
       0.0        , 0.0        , 0.0        , 0.0        , 0.0       , &
       0.0        , 0.0        , 0.0        , 0.0        , 0.0       , &
# endif
       3.2224e-07 , 1.6605e-06 , 7.6565e-06 , 3.1592e-05 , 0.00011664, &
       0.00038538 , 0.0011394  , 0.0030144  , 0.0071362  , 0.015117  , &
       0.028657   , 0.048612   , 0.073789   , 0.10023    , 0.12182   , &
       0.1325     , 0.12895    , 0.11231    , 0.087525   , 0.061037  , &
       0.038089   , 0.02127    , 0.010628   , 0.0047523  , 0.0019015 , &
       0.00068081 , 0.00021813 , 6.2536e-05 , 1.6044e-05 , 3.6831e-06/
#endif

  ! Subgrid coagulation timescale (win, 10/28/08)
  REAL*4 :: SGCTSCALE
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER, PRIVATE :: id_AW1
  INTEGER, PRIVATE :: id_DUST1
  INTEGER, PRIVATE :: id_ECIL1
  INTEGER, PRIVATE :: id_ECOB1
  INTEGER, PRIVATE :: id_H2SO4
  INTEGER, PRIVATE :: id_NH3
  INTEGER, PRIVATE :: id_NH4
  INTEGER, PRIVATE :: id_NK1
  INTEGER, PRIVATE :: id_OCIL1
  INTEGER, PRIVATE :: id_OCOB1
  INTEGER, PRIVATE :: id_SF1
  INTEGER, PRIVATE :: id_SO4
  INTEGER, PRIVATE :: id_SS1

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_tomas
!
! !DESCRIPTION: Subroutine DO\_TOMAS is the driver subroutine that calls the
!  appropriate aerosol microphysics and dry deposition subroutines. (win,
!  7/23/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_TOMAS( Input_Opt,  State_Chm, State_Diag, State_Grid, &
                       State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
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
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! DO_TOMAS begins here
    !=================================================================

    ! Assume success
    RC    = GC_SUCCESS

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine DO_TOMAS in tomas_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
    ENDIF

    ! Do TOMAS aerosol microphysics
    CALL AEROPHYS( Input_Opt, State_Chm, State_Grid, State_Met, RC )

#if defined ( DEBUG )
    print *, 'call checkmn in tomas_mod:222'
#endif
    CALL CHECKMN( 0, 0, 0, Input_Opt, State_Chm, State_Grid, &
                  State_Met,'Before Aerodrydep', RC)

    ! in kg

    ! Do dry deposition
    IF ( Input_Opt%LDRYD ) THEN
       CALL AERO_DRYDEP( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
    ENDIF

    ! not in kg

#if defined ( DEBUG )
    print *, 'call checkmn in tomas_mod:229'
#endif
    CALL CHECKMN( 0, 0, 0, Input_Opt, State_Chm, State_Grid, &
                  State_Met, 'Before exiting DO_TOMAS', RC )

  END SUBROUTINE DO_TOMAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerophys
!
! !DESCRIPTION: Subroutine AEROPHYS does aerosol microphysics, including
!  nucleation, coagulation, and condensation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AEROPHYS( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG_MOD,           ONLY : AD61,  AD61_INST
#endif
    USE ErrCode_Mod
    USE ERROR_MOD
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
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER             :: I, J, L, N, JC, K !counters
    INTEGER             :: MPNUM    !microphysical process id #
    REAL*4              :: ADT      !aerosol microphysics time step (seconds)
    REAL(fp)            :: QSAT     !used in RH calculation
    INTEGER             :: TRACNUM
    REAL(fp)            :: FRAC
    CHARACTER(LEN=255)  :: MSG, LOC ! species unit check

    ! Arguments for CHECK_VALUE; avoids array temporaries (bmy, 1/28/14)
    CHARACTER(LEN=255) :: ERR_VAR
    CHARACTER(LEN=255) :: ERR_MSG
    INTEGER            :: ERR_IND(4)

    !---------
    !sfarina - move definitions of these from module common
    !to within each sub and pass them around for openmp
    REAL(fp)           :: Nk(IBINS), Nkd(IBINS)
    REAL(fp)           :: Mk(IBINS, ICOMP)
    REAL(fp)           :: Mkd(IBINS,ICOMP)
    REAL(fp)           :: Gc(ICOMP - 1)
    REAL(fp)           :: Gcd(ICOMP - 1)

    REAL*4             :: BOXVOL,  BOXMASS, TEMPTMS
    REAL*4             :: PRES,    RHTOMAS

    REAL(fp)           ::  surf_area     ! aerosol surface area [micon^2 cm^-3]
    REAL(fp)           ::  ionrate       ! ion pair formation
                                         ! rate [ion pairs cm^-3 s^-1]
    !---------

    REAL(fp)          :: Nkout(ibins), Mkout(ibins,icomp)
    REAL(fp)          :: Gcout(icomp-1)
    REAL(fp)          :: Nknuc(ibins), Mknuc(ibins,icomp)
    REAL(fp)          :: Nkcond(ibins), Mkcond(ibins,icomp)
    REAL(fp)          :: fn  ! nucleation rate of clusters cm-3 s-1
    REAL(fp)          :: fn1 ! formation rate of particles to first size bin cm-3 s-1
    REAL(fp)          :: nucrate(State_Grid%NY,State_Grid%NZ)
    REAL(fp)          :: nucrate1(State_Grid%NY,State_Grid%NZ)
    REAL(fp)          :: tot_n_1, tot_n_1a, tot_n_2, tot_n_i ! used for nitrogen mass checks
    REAL(fp)          :: tot_s_1, tot_s_1a, tot_s_2 ! used for sulfur mass checks
    REAL(fp)          :: h2so4rate_o ! H2SO4rate for the specific grid cell
    REAL(fp)          :: TOT_Mk, TOT_nk  ! for checking mass and number

    REAL(fp)          :: transfer(ibins)
    LOGICAL           :: PRINTNEG  !<step4.0-temp> (win, 3/24/05)
    logical           :: ERRORSWITCH  !<step4.2> To see where mnfix found negative value (win, 9/12/05)
    logical           :: errspot   !<step4.4> To see where so4cond found errors (win, 9/21/05)
    logical           :: printdebug !<step4.3> Print out for debugging (win, 9/16/05)
    logical           :: COND, COAG, NUCL !<step5.1> switch for each process (win 4/8/06)
    integer :: iob, job,lob !Just declare in case I want to debug (4/8/06)
    data       iob, job, lob /  1  ,       1    ,       1 /
    real(fp)           :: NH3_to_NH4, CEPS
    parameter ( CEPS=1.e-17_fp )

    real(fp) igR
    parameter (igR=8.314) !Ideal gas constant J/mol.K

    !The following are constants used in calculating rel. humidity
    real(fp) axcons, bxcons, bytf, tf  !for RH calculation
    parameter(axcons=1.8094520287589733, &
              bxcons=0.0021672473136556273, &
              bytf=0.0036608580560606877, tf=273.16)
    !lhe and lhs are the latent heats of evaporation and sublimation

    logical, save     :: firsttime = .true.
    integer           :: num_iter

    real(fp)    cplevels(9) ! cosmic ray pressure levels (for interp)
    data        cplevels  / 959., 894., 786., 634., 468., &
                            321., 201., 103., 27. /

    integer     lev
    real(fp)    weight

    double precision soil_ions(9) ! ion pairs cm-3 s-1 from radioactive elements in soil
    data             soil_ions / 5.,3.3,1.8,0.7,0.,0.,0.,0.,0./

    ! For fields from Input_Opt
    LOGICAL :: prtDebug

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! AEROPHYS begins here
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Check that species units are in [kg]. While species units
    ! are now generally [kg/kg] in GEOS-Chem, they are converted to
    ! kg for TOMAS elsewhere in tomas_mod prior to calling this subroutine
    ! (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine AEROPHYS in tomas_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
    ENDIF

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    ! Get logical values from Input_Opt
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Initialize debugging and error-signal switches
    printneg    = .FALSE.
    errorswitch = .FALSE.
    PRINTDEBUG  = .FALSE.
    ERRSPOT     = .FALSE.

    ! Initialize switches for each microphysical process
    COND = .TRUE.
    COAG = .TRUE.
    NUCL = .TRUE.

    ! Initialize nucleation rate arrays
    nucrate(:,:)  = 0.e+0_fp
    nucrate1(:,:) = 0.e+0_fp

    ! First-time setup
    if (firsttime) then

       !====================================================================
       ! Make sure there is access to the H2SO4 production rate array
       ! H2SO4RATE, which saves the H2SO4 production rate for EACH chemistry
       ! timestep.  The prod/loss family has to be set with at least one
       ! with the family name PSO4 and SO4 as its one member. (win, 9/30/08)
       !====================================================================

       DO N = 1, Input_Opt%NFAM
          ! If family name 'PSO4' is found, then skip error-stop
          IF ( Input_Opt%FAM_NAME(N) == 'PSO4') GOTO 1
       ENDDO
       ! Family name 'PSO4' not found... exit with error message
       write(*,*)'-----------------------------------------------'
       write(*,*)' Need to setup ND65 family PSO4 with SO4 as '
       write(*,*)' a member to have H2SO4RATE array '
       write(*,*)'  ... need H2SO4RATE for nucl & cond in TOMAS'
       write(*,*)'-----------------------------------------------'
       CALL ERROR_STOP('AEROPHYS','Enter microphys')
1      CONTINUE

       write(*,*) 'AEROPHYS: This run uses coupled condensation-', &
                  'nucleation scheme with pseudo-steady state H2SO4'
       if(tern_nuc == 1) then
          write(*,*)'  Nucleation: Ternary (Napari ', &
                    'et al. 2002) and Binary (Vehkamaki et al. 2002)'
       else
          write(*,*)'  Nucleation: Binary (Vehkamaki et al. 2002)'
       endif

       firsttime = .false.
    endif

    ! Get chemistry timestep for use below [s]
    ! NOTE: This doesn't have to be !$OMP+PRIVATE (bmy, 2/7/20)
    ADT = GET_TS_CHEM()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: THIS PARALLEL LOOP MAY BE ABLE TO BE REVERSED TO L-J-I
    !%%% WHICH IS MUCH MORE EFFICIENT (bmy, 1/28/14)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED )   &
    !$OMP PRIVATE( I, J, L )  &
    !$OMP PRIVATE( PRES, TEMPTMS, BOXMASS, RHTOMAS, BOXVOL )       &
    !$OMP PRIVATE( printneg, ionrate, lev, weight, GC, N, NK, JC ) &
    !$OMP PRIVATE( MK, H2SO4rate_o, tot_n_1, k, tot_s_1, MPNUM )   &
    !$OMP PRIVATE( Nkd, Mkd, TOT_NK, TOT_MK, TRANSFER )            &
    !$OMP PRIVATE( Nkout,Mkout,Gcout,fn,fn1 )                      &
    !$OMP PRIVATE( num_iter,Nknuc,Mknuc,Nkcond )                   &
    !$OMP PRIVATE( Mkcond, ERRORSWITCH, tot_s_1a, tot_n_1a )       &
    !$OMP PRIVATE( ERR_VAR, ERR_MSG, ERR_IND )                     &
    !$OMP PRIVATE( TRACNUM, NH3_TO_NH4, SURF_AREA )                &
    !$OMP SCHEDULE( DYNAMIC )
    DO I = 1, State_Grid%NX
    DO J = 1, State_Grid%NY
    DO L = 1, State_Grid%NZ

#ifdef BPCH_DIAG
       ! Reset the AD61_INST array used for tracking instantaneous
       ! certain rates.  As of now, tracking nucleation (win, 10/7/08)
       IF ( Input_Opt%ND61 > 0 ) THEN
          AD61_INST(I,J,L,:) = 0e0
       ENDIF
#endif

       ! Skip non-chemgrid boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       !vbn write(890,89)I,J,L,Spc(i,j,l,id_H2SO4)

       !if(printdebug) print *,'+++++',I,J,L,'inside Aerophys'

       ! Get info on this grid box
       PRES    = State_Met%PMID(i,j,l)*100.0 ! in Pa
       TEMPTMS = State_Met%T(I,J,L)
       BOXMASS = State_Met%AD(I,J,L)
       RHTOMAS = State_Met%RH(I,J,L)/ 1.e2
       IF ( RHTOMAS > 0.99 ) RHTOMAS = 0.99
       BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3

       printneg = .FALSE.

       ! determine ion rate
       !ionrate = 10. ! set as constant now !!jrp

       IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN

          if( (pres / 100.) .lt. cplevels(9) ) then
             ionrate = cosmic_ions(i,j,9) * boxmass / boxvol
             if(State_Met%FRCLND(I,J) .gt. 0.2) then
                ionrate = ionrate + soil_ions(9)
             endif
          elseif((pres/100.) .gt. cplevels(1)) then
             ionrate = cosmic_ions(i,j,1) * boxmass / boxvol
             if( State_Met%FRCLND(I,J) .gt. 0.2 ) then
                ionrate = ionrate + soil_ions(1)
             endif
          else
             lev=2
             do while (pres / 100. .lt. cplevels(lev))
                lev=lev+1
             enddo
             weight=( cplevels( lev - 1 ) - pres / 100. ) / &
                    ( cplevels( lev - 1 ) - cplevels(lev) )
             ionrate=( cosmic_ions(i,j,lev  ) * weight + &
                       cosmic_ions(i,j,lev-1) * (1.e+0_fp - weight) ) &
                       * boxmass / boxvol
             if( State_Met%FRCLND(I,J) .gt. 0.2) then
                ionrate=ionrate + ( soil_ions( lev   ) * weight + &
                        soil_ions( lev-1 ) * (1.e+0_fp-weight) )
             endif
          endif

       ELSE
          ionrate = 0.e+0_fp
       ENDIF

       if(ionrate .le. 1.501) ionrate = 1.501

       !print*,'i',i,'j',j,'l',l,'ionrate',ionrate

       ! Initialize all condensible gas values to zero
       ! Gc(srtso4) will remain zero until within cond_nuc where the
       ! pseudo steady state H2SO4 concentration will be put in this place.
       DO JC=1, ICOMP-1
          Gc(JC) = 0.e+0_fp
       ENDDO

       ! Swap Spc into Nk, Mk, Gc arrays
       DO N = 1, IBINS
          NK(N) = Spc(I,J,L,id_NK1-1+N)
          DO JC = 1, ICOMP-IDIAG
             MK(N,JC) = Spc(I,J,L,id_NK1-1+N+JC*IBINS)

             IF( IT_IS_NAN( MK(N,JC) ) ) THEN
                PRINT *,'+++++++ Found NaN in AEROPHYS ++++++++'
                PRINT *,'Location (I,J,L):',I,J,L,'Bin',N,'comp',JC
             ENDIF

          ENDDO
          MK(N,SRTH2O) = Spc(I,J,L,id_AW1-1+N)
       ENDDO

       ! Get NH4 mass from the bulk mass and scale to bin with sulfate
       IF ( SRTNH4 > 0 ) THEN
          CALL NH4BULKTOBIN( MK(:,SRTSO4), Spc(I,J,L,id_NH4), TRANSFER )
          MK(1:IBINS,SRTNH4) = TRANSFER(1:IBINS)
          Gc(SRTNH4) = Spc(I,J,L,id_NH3)
       ENDIF

       ! Give it the pseudo-steady state value instead later (win,9/30/08)
       !GC(SRTSO4) = Spc(I,J,L,id_H2SO4)

       H2SO4rate_o = H2SO4_RATE(I,J,L)  ! [kg s-1]
       ! Pengfei Liu add 2018/04/18, debug
       IF ( H2SO4rate_o .lt. 0.e0 ) THEN
          Print*, 'Debug TOMAS: H2SO4RATE = ', H2SO4rate_o, 'I = ', I, &
               'J = ', J, 'L = ', L
          H2SO4rate_o = 0.e+0_fp
       ENDIF
       !
       IF ( I == 10 .and. J == 10 .and. L == 10 ) THEN
          Print*, 'Debug TOMAS: H2SO4RATE =', H2SO4rate_o
       ENDIF

       ! nitrogen and sulfur mass checks
       ! get the total mass of N
       tot_n_1 = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
       do k=1,ibins
          tot_n_1 = tot_n_1 + Mk(k,srtnh4)*14.e+0_fp/18.e+0_fp
       enddo

       ! get the total mass of S
       tot_s_1 = H2SO4rate_o*adt*32.e+0_fp/98.e+0_fp
       do k=1,ibins
          tot_s_1 = tot_s_1 + Mk(k,srtso4)*32.e+0_fp/96.e+0_fp
       enddo


       !if (printdebug.and.i==iob .and. j==job .and. l==lob ) then
       !   CALL DEBUGPRINT( Nk, Mk, I, J, L, 'Begin aerophys' )
       !   print *,'H2SO4RATE ',H2SO4rate_o
       !endif

       !*********************
       ! Aerosol dynamics
       !*********************

       !Do water eqm at appropriate times
       CALL EZWATEREQM( MK, RHTOMAS )

       !Fix any inconsistencies in M/N distribution (because of advection)
       CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

       !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.

       !print *, 'mnfix in tomas_mod:533'
       CALL MNFIX( NK, MK, ERRORSWITCH )
       IF ( ERRORSWITCH ) THEN
          PRINT *,'Aerophys: MNFIX found error at',I,J,L
          CALL ERROR_STOP('AEROPHYS-MNFIX (1)','Enter microphys')
       ENDIF

       MPNUM = 5
#ifdef BPCH_DIAG
       IF ( ND60 > 0 ) THEN
          CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
       ENDIF
#endif

       !IF ( printdebug.and.i==iob .and. j==job .and. l==lob ) THEN
       !   CALL DEBUGPRINT( Nk, Mk, I, J, L, 'After mnfix before cond/nucl' )
       !ENDIF

       ! Before doing any cond/nucl/coag, check if there's any aerosol in
       ! the current box
       TOT_NK = 0.e+0_fp
       TOT_MK = 0.e+0_fp
       do k = 1, ibins
          TOT_NK = TOT_NK + Nk(K)
          do jc=1, icomp-idiag
             TOT_MK = TOT_MK + Mk(k,jc)
          enddo
       enddo

       if(TOT_NK .lt. 1.e-10_fp) then
          if( .NOT. SPINUP(5.0)) then
             print *,'No aerosol in box ',I,J,L,'-->SKIP'
          endif
          CYCLE
       endif

       !-------------------------------------
       ! Condensation and nucleation (coupled)
       !-------------------------------------
       IF ( COND .AND. NUCL ) THEN

          !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.

          CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

          !sfdebug if(printdebug) then
          !sfdebug    !print*,'Before COND_NUC Gc(srtso4)=',Gc(srtso4)
          !sfdebug    do N = 1,IBINS
          !sfdebug       IF( IT_IS_NAN( NK(N) ) ) THEN
          !sfdebug          print*, "found NAN in nk", n, nk
          !sfdebug       endif
          !sfdebug       DO JC=1, ICOMP-1
          !sfdebug          IF( IT_IS_NAN( Gc(JC) ) ) THEN
          !sfdebug             print*, "found NAN in gc", jc, gc
          !sfdebug          endif
          !sfdebug       ENDDO
          !sfdebug    enddo
          !sfdebug endif

          CALL COND_NUC(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1, &
                        H2SO4rate_o,adt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, &
                        ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, &
                        RHTOMAS, ERRORSWITCH, l)

          !sfdebug if(printdebug) then
          !sfdebug    !print*,'Before COND_NUC Gc(srtso4)=',Gc(srtso4)
          !sfdebug    do N = 1,IBINS
          !sfdebug       IF( IT_IS_NAN( NKout(N) ) ) THEN
          !sfdebug          print*, "found NAN in nkout", n, nkout
          !sfdebug       endif
          !sfdebug       DO JC=1, ICOMP-1
          !sfdebug          IF( IT_IS_NAN( Gcout(JC) ) ) THEN
          !sfdebug             print*, "found NAN in gcout", jc, gcout
          !sfdebug          endif
          !sfdebug       ENDDO
          !sfdebug    enddo
          !sfdebug endif

          IF ( ERRORSWITCH ) THEN
             PRINT *,'Aerophys: found error at',I,J,L
             CALL ERROR_STOP('AEROPHYS','After cond_nuc')
          ENDIF

          ERR_VAR = 'Gcout'
          ERR_MSG = 'After COND_NUC'
          ! check for NaN and Inf (win, 10/4/08)
          do jc = 1, icomp-1
             ERR_IND(1) = I
             ERR_IND(2) = J
             ERR_IND(3) = L
             ERR_IND(4) = 0
             call check_value( Gcout(jc), ERR_IND, ERR_VAR, ERR_MSG )

             !if( IT_IS_FINITE(Gcout(jc))) then
             !   print *,'xxxxxxxxx Found Inf in Gcout xxxxxxxxxxxxxx'
             !   print *,'Location ',I,J,L, 'comp',jc
             !   call debugprint( Nkout, Mkout, i,j,l,'After COND_NUC')
             !   stop
             !endif
          enddo

          !get nucleation diagnostic
          DO N = 1, IBINS
             NK(N) = NKnuc(N)
             DO JC = 1, ICOMP
                MK(N,JC) = MKnuc(N,JC)
             ENDDO
          ENDDO

          MPNUM = 3
#ifdef BPCH_DIAG
          IF ( ND60 > 0 ) THEN
             CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                            State_Grid )
          ENDIF
#endif

          MPNUM = 7
#ifdef BPCH_DIAG
          IF ( ND61 > 0 )  THEN
             CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                            State_Grid )
          ENDIF
#endif

          IF ( printdebug.and.i==iob .and. j==job .and. l==lob )  THEN
             CALL DEBUGPRINT( Nk, Mk, I, J, L,'After nucleation' )
          ENDIF

          !get condensation diagnostic
          DO N = 1, IBINS
             NK(N) = NKcond(N)
             DO JC = 1, ICOMP
                MK(N,JC) = MKcond(N,JC)
             ENDDO
          ENDDO

          Gc(srtnh4)=Gcout(srtnh4)
          Gc(srtso4)=Gcout(srtso4)

          MPNUM = 1
#ifdef BPCH_DIAG
          IF ( ND60 > 0 ) THEN
             CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                            State_Grid )
          ENDIF
#endif

          IF ( printdebug.and.i==iob .and. j==job .and. l==lob ) THEN
             CALL DEBUGPRINT( Nk, Mk, I, J, L,'After condensation' )
          ENDIF

          nucrate(j,l)=nucrate(j,l)+fn
          nucrate1(j,l)=nucrate1(j,l)+fn1

          ! Write nucleation rate to diagnostric ND61 (win, 10/6/08)
#ifdef BPCH_DIAG
          IF ( ND61 > 0 ) THEN
             IF ( L <= LD61 ) AD61(I,J,L,2) = AD61(I,J,L,2) + fn

             ! Tracks nucleation rates instantaneously for planeflight
             AD61_INST(I,J,L,2) = fn
          ENDIF
#endif

          DO N = 1, IBINS
             NK(N) = NKout(N)
             DO JC = 1, ICOMP
                MK(N,JC) = MKout(N,JC)
             ENDDO
          ENDDO

       ENDIF

       ! nitrogen and sulfur mass checks
       ! get the total mass of N
       tot_n_1a = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
       do k=1,ibins
          tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.e+0_fp/18.e+0_fp
       enddo

       ! get the total mass of S
       tot_s_1a = 0.e+0_fp
       do k=1,ibins
          tot_s_1a = tot_s_1a + Mk(k,srtso4)*32.e+0_fp/96.e+0_fp
       enddo

       CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
       !print *, 'mnfix in tomas_mod:677'

       CALL MNFIX( Nk, Mk, ERRORSWITCH )
       IF ( ERRORSWITCH ) THEN
          PRINT *,'Aerophys: MNFIX found error at',I,J,L
          IF( .not. SPINUP(14.0) ) THEN
             CALL ERROR_STOP('AEROPHYS-MNFIX (2)','After cond/nucl')
          ELSE
             PRINT *,'Let error go during spin up'
          ENDIF
       ENDIF

       MPNUM = 5
#ifdef BPCH_DIAG
       IF ( ND60 > 0 ) THEN
          CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                         State_Grid )
       ENDIF
#endif

       !-----------------------------
       ! Coagulation
       !-----------------------------

       if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.
       !if (i==iob .and. j==job .and. l==lob ) &
       !    CALL DEBUGPRINT( Nk, Mk, I, J, L,'Before coagulation' )

       IF( COAG )  THEN
          CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
          CALL MULTICOAG( ADT, Nk, Mk, BOXVOL, PRES, TEMPTMS, errorswitch )

          if ( errorswitch ) &
               CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )
          !if (i==iob .and. j==job .and. l==lob ) &
          !    CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )

          MPNUM = 2
#ifdef BPCH_DIAG
          IF ( ND60 > 0 ) THEN
             CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                            State_Grid )
          ENDIF
#endif

          !Fix any inconsistency after coagulation (win, 4/18/06)
          CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
          if(printdebug .and. i==iob.and.j==job.and.l==lob) &
               ERRORSWITCH=.true. !4/18/06 win

          !print *, 'mnfix in tomas_mod:719'
          CALL MNFIX( NK, MK, ERRORSWITCH )

          IF ( ERRORSWITCH ) THEN
             PRINT *,'MNFIX found error at',I,J,L
             IF( .not. SPINUP(14.0) ) THEN
                CALL ERROR_STOP('AEROPHYS-MNFIX (3)', 'After COAGULATION'  )
             ELSE
                PRINT *,'Let error go during spin up'
             ENDIF
          ENDIF

          MPNUM = 5
#ifdef BPCH_DIAG
          IF ( ND60 > 0 ) THEN
             CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                            State_Grid )
          ENDIF
#endif
       ENDIF  ! Coagulation

       ! Do water eqm at appropriate times
       CALL EZNH3EQM( Gc, Mk )
       CALL EZWATEREQM ( MK, RHTOMAS )

       !****************************
       ! End of aerosol dynamics
       !****************************

       !Fix any inconsistencies in M/N distribution (because of advection)
       CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

       ! Make sure anything that leaves AEROPHYS is free of any error
       ! This MNFIX call could be temporary (?) or just leave it here and
       ! monitor if the error fixed is significantly large meaning some
       ! serious problem needs to be investigated
       if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.true.

       !print *, 'mnfix in tomas_mod:758'
       CALL MNFIX(NK,MK,ERRORSWITCH)
       IF ( ERRORSWITCH ) THEN
          PRINT *,'End of Aerophys: MNFIX found error at',I,J,L
          IF( .not. SPINUP(14.0) ) THEN
             CALL ERROR_STOP('AEROPHYS-MNFIX (4)', 'End of microphysics')
          ELSE
             PRINT *,'Let error go during spin up'
          ENDIF
       ENDIF

       ! Accumulate changes by mnfix to diagnostic (win, 9/8/05)
       MPNUM = 5
#ifdef BPCH_DIAG
       IF ( ND60 > 0 ) THEN
          CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, &
                         State_Grid )
       ENDIF
#endif

       ! Swap Nk, Mk, and Gc arrays back to Spc
       DO N = 1, IBINS
          TRACNUM = id_NK1 - 1 + N
          Spc(I,J,L,TRACNUM) = NK(N)
          DO JC = 1, ICOMP-IDIAG
             TRACNUM = id_NK1 - 1 + N + IBINS * JC
             Spc(I,J,L,TRACNUM) = MK(N,JC)
          ENDDO
          Spc(I,J,L,id_AW1-1+N) = MK(N,SRTH2O)
       ENDDO
       Spc(I,J,L,id_H2SO4) = GC(SRTSO4)

       ! print to file to check mass conserv
       !write(*,77) I,J,L, Spc(I,J,L,id_NH3), Spc(I,J,L,id_NH3)-GC(SRTNH4)

       ! Calculate NH3 gas lost to aerosol phase as NH4
       NH3_to_NH4 = Spc(I,J,L,id_NH3)-GC(SRTNH4)

       ! Update the bulk NH4 aerosol species
       if ( NH3_to_NH4 > 0e+0_fp ) &
            Spc(I,J,L,id_NH4) = Spc(I,J,L,id_NH4) + &
                                NH3_to_NH4/17.e+0_fp*18.e+0_fp

       ! Update NH3 gas species (win, 10/6/08)
       ! plus tiny amount CEPS in case zero causes some problem
       Spc(I,J,L,id_NH3)   = GC(SRTNH4) + CEPS !MUST CHECK THIS!! (win,9/26/08)


       !vbn write(889,89)I,J,L,Spc(i,j,l,id_H2SO4)
89     format(3I3,'Spc(id_H2SO4) kg', E13.5)

    ENDDO                     !L loop
    ENDDO                     !J loop
    ENDDO                     !I loop
    !$OMP END PARALLEL DO

    !WRITE(777,*) '---------------------------'
77  FORMAT(3I4, '  Spc(id_NH3),'E13.5,'  Used', E13.5 )

    IF ( COND .and. prtDebug ) PRINT *,'### AEROPHYS: SO4 CONDENSATION'
    IF ( COAG .and. prtDebug ) PRINT *,'### AEROPHYS: COAGULATION'
    IF ( NUCL .and. prtDebug ) PRINT *,'### AEROPHYS: NUCLEATION'

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE AEROPHYS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cond_nuc
!
! !DESCRIPTION: This subroutine calculates the change in the aerosol size
!  distribution due to so4 condensation and binary/ternary nucleation during
!  the overal microphysics timestep.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!  Put in GEOS-Chem by Win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COND_NUC(Nki,Mki,Gci,Nkf,Mkf,Gcf,fnavg,fn1avg, &
                      H2SO4rate,dti,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, &
                      ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, &
                      RHTOMAS, errswitch, lev)
!
! !INPUT PARAMETERS:
!
    ! Nki(ibins)        - number of particles per size bin in grid cell
    ! Nnuci             - number of nucleation size particles per size bin in
    !                     grid cell
    ! Mnuci             - mass of given species in nucleation pseudo-bin
    !                     (kg/grid cell)
    ! Mki(ibins, icomp) - mass of a given species per size bin/grid cell
    ! Gci(icomp-1)      - amount (kg/grid cell) of all species present in the
    !                     gas phase except water
    ! H2SO4rate         - rate of H2SO4 chemical production [kg s^-1]
    ! dt                - total model time step to be taken (s)
    REAL(fp) Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
    double precision H2SO4rate
    real             dti
!
! !OUTPUT PARAMETERS:
!
    ! Nkf, Mkf, Gcf  - same as above, but final values
    ! Nknuc,  Mknuc  - same as above, final values from just nucleation
    ! Nkcond, Mkcond - same as above, but final values from just condensation
    ! fn, fn1
    REAL(fp) Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
    REAL(fp) Nknuc(ibins), Mknuc(ibins, icomp)
    REAL(fp) Nkcond(ibins),Mkcond(ibins,icomp)
    double precision fnavg        ! nucleation rate of clusters cm-3 s-1
    double precision fn1avg       ! formation rate of particles to first size bin cm-3 s-1
    REAL*4           BOXVOL, BOXMASS, TEMPTMS, RHTOMAS, PRES
    logical          errswitch    ! signal for error
    integer          lev          ! layer of the model
    REAL(fp)   surf_area
    REAL(fp)   ionrate
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision dt
    integer          i,j,k,c      ! counters
    double precision fn           ! nucleation rate of clusters cm-3 s-1
    double precision fn1          ! formation rate of particles to first size bin cm-3 s-1
    double precision pi, R        ! pi and gas constant (J/mol K)
    double precision CSi,CSa      ! intial and average condensation sinks
    double precision CS1,CS2      ! guesses for condensation sink [s^-1]
    double precision CStest       ! guess for condensation sink
    REAL(fp)           Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
    REAL(fp)           Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
    REAL(fp)           Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
    logical          nflg         ! returned from nucleation, says whether nucleation occurred or not
    double precision mcond,mcond1 ! mass to condense [kg]
    double precision tol          ! tolerance
    double precision eps          ! small number
    double precision sinkfrac(ibins) ! fraction of condensation sink coming from bin k
    double precision totmass      ! the total mass of H2SO4 generated during the timestep
    double precision tmass
    double precision CSch         ! fractional change in condensation sink
    double precision CSch_tol     ! tolerance in change in condensation sink
    double precision addt         ! adaptive timestep time
    double precision time_rem     ! time remaining
    integer          num_iter     ! number of iteration
    double precision sumH2SO4     ! used for finding average H2SO4 conc over timestep
    integer          iter         ! number of iteration
    double precision rnuc         ! critical radius [nm]
    double precision gasConc      ! gas concentration [kg]
    double precision mass_change  ! change in mass during nucleation
    double precision total_nh4_1,total_nh4_2
    double precision min_tstep    ! minimum timestep [s]
    integer          nuc_bin      ! the nucleation bin
    double precision sumfn, sumfn1 ! used for getting average nucleation rates
    logical          tempvar,  pdbg
    real(fp)           tnumb
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(eps=1E-40)
    parameter(CSch_tol=0.01)
    parameter(min_tstep=1.0e+0_fp)

    !=================================================================
    ! COND_NUC begins here
    !=================================================================

    pdbg      = errswitch ! transfer the signal to print debug from outside
    errswitch = .false.   ! flag error to outide to terminate program.

    dt = dble(dti)

    ! Initialize values of Nkf, Mkf, Gcf, and time
    do j=1,icomp-1
       Gc1(j)=Gci(j)
       Gcf(j)=Gci(j)
    enddo
    do k=1,ibins
       Nk1(k)=Nki(k)
       Nknuc(k)=Nki(k)
       Nkcond(k)=Nki(k)
       do j=1,icomp
          Mk1(k,j)=Mki(k,j)
          Mknuc(k,j)=Mki(k,j)
          Mkcond(k,j)=Mki(k,j)
       enddo
    enddo

    ! Get initial condensation sink
    CS1 = 0.e+0_fp
    call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac,surf_area,BOXVOL,TEMPTMS,PRES)
    if( pdbg) print*,'CS1', CS1
    !CS1 = max(CS1,eps)

    !Get initial H2SO4 concentration guess (assuming no nucleation)
    !Make sure that H2SO4 concentration doesn't exceed the amount generated
    !during that timestep (this will happen when the condensation sink is very low)

    ! get the steady state H2SO4 concentration
    call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4), &
                      gasConc, ionrate, surf_area, &
                      BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
    if( pdbg) print*,'gasConc',gasConc
    Gc1(srtso4) = gasConc
    addt = min_tstep
    !addt = 3600.e+0_fp
    totmass = H2SO4rate*addt*96.e+0_fp/98.e+0_fp

    tempvar = pdbg

    !Get change size distribution due to nucleation with initial guess
    call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin, &
                    addt, ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, &
                    PRES, RHTOMAS, PDBG, lev)

    if(pdbg) then
       print*,'COND_NUC: Found an error at nucleation --> TERMINATE'
       errswitch = .true.
       return
    endif
    pdbg = tempvar !put the print debug switch back to pdbg
    if(pdbg) call debugprint(Nk2, Mk2, 0,0,0,'After nucleation[1]')

    !print*,'after nucleation'
    !print*,'Nnuc1',Nnuc1
    !print*,'Nnuc2',Nnuc2
    !print*,'Mnuc1',Mnuc1
    !print*,'Mnuc2',Mnuc2

    mass_change = 0.e+0_fp

    do k=1,ibins
       mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
    enddo
    if( pdbg)  print*,'mass_change',mass_change

    mcond = totmass-mass_change ! mass of h2so4 to condense

    if( pdbg) print*,'after nucleation'
    if( pdbg)  print*,'totmass',totmass,'mass_change1',mass_change,'mcond',mcond
    if( pdbg)  print*,'cs1',CS1, Gc1(srtso4)

    if (mcond.lt.0.e+0_fp)then
       tmass = 0.e+0_fp
       do k=1,ibins
          do j=1,icomp-idiag
             tmass = tmass + Mk2(k,j)
          enddo
       enddo
       !if (abs(mcond).gt.tmass*1.0D-8) then
       if (abs(mcond).gt.totmass*1.0e-8_fp) then
          if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
             !if (CS1.gt.1.0D-5)then
             !   print*,'budget fudge 1 in cond_nuc'
             !endif
             tmass = 0.e+0_fp
             do j=1,icomp-idiag
                tmass = tmass + Mk2(nuc_bin,j)
             enddo
             Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
             Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
             mcond = 0.e+0_fp
          else
             print*,'budget fudge 2 in cond_nuc'
             do k=2,ibins
                Nk2(k) = Nk1(k)
                Mk2(k,srtso4) = Mk1(k,srtso4)
             enddo
             Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
             Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
             mcond = 0.e+0_fp
             !print*,'mcond < 0 in cond_nuc', mcond, totmass
             !stop
          endif
       else
          mcond = 0.e+0_fp
       endif
    endif

    !if (mcond.lt.0.e+0_fp)then
    !   print*,'mcond < 0 in cond_nuc', mcond
    !   stop
    !endif
    tmass = 0.e+0_fp
    do k=1,ibins
       do j=1,icomp-idiag
          tmass = tmass + Mk2(k,j)
       enddo
    enddo
    if( pdbg)  print*, 'mcond',mcond,'tmass',tmass,'nuc',Nk2(1)-Nk1(1)
    tempvar = pdbg

    ! Get guess for condensation
    call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area, &
                BOXVOL, TEMPTMS, PRES, pdbg )

    if(pdbg) then
       print*,'COND_NUC: Found an error at EZCOND --> TERMINATE'
       errswitch = .true.
       return
    endif
    pdbg = tempvar
    if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[1]')
    !print*,'after ezcond',Nk2,Nk3
    !jrp mcond1 = 0.e+0_fp
    !jrp do k=1,ibins
    !jrp    do j=1,icomp
    !jrp       mcond1 = mcond1 + (Mk3(k,j)-Mk2(k,j))
    !jrp    enddo
    !jrp enddo
    !print*,'mcond',mcond,'mcond1',mcond1

    Gc3(srtnh4) = Gc1(srtnh4)

    call eznh3eqm(Gc3,Mk3)
    call ezwatereqm(Mk3, RHTOMAS)

    ! check to see how much condensation sink changed
    call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area, &
                     BOXVOL,TEMPTMS, PRES)
    CSch = abs(CS2 - CS1)/CS1

    !if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
    ! get starting adaptive timestep to not allow condensationk sink
    ! to change that much
    ! Avoid div-by-zero (bmy, 1/28/14)
    IF ( ABS( CSch ) > 0e+0_fp ) THEN
       addt = addt*CSch_tol/CSch/2e+0_fp
    ELSE
       addt = 0e+0_fp
    ENDIF
    addt = min(addt,dt)
    addt = max(addt,min_tstep)

    time_rem = dt ! time remaining
    if( pdbg)    print*,'addt',addt,time_rem
    num_iter = 0
    sumH2SO4=0.e+0_fp
    sumfn = 0.e+0_fp
    sumfn1 = 0.e+0_fp
    ! do adaptive timesteps
    do while (time_rem .gt. 0.e+0_fp)
       num_iter = num_iter + 1
       if( pdbg) print*, 'iter', num_iter, ' addt', addt, 'time_rem', time_rem
       ! get the steady state H2SO4 concentration
       if (num_iter.gt.1)then ! no need to recalculate for first step
          call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4), &
                            gasConc, ionrate, surf_area, &
                            BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
          Gc1(srtso4) = gasConc
       endif
       if( pdbg)    print*,'gasConc',gasConc

       sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
       totmass = H2SO4rate*addt*96.e+0_fp/98.e+0_fp
       !call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2, &
       !                Mk2,Gc2,Nnuc2,Mnuc2,nflg,lev)

       !Debug to see what goes in nucleation (win, 10/3/08)
       if(pdbg) then
          print*,'Temperature',TEMPTMS,'RH',RHTOMAS
          print*,'H2SO4',Gc1(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
          print*,'NH3ppt',Gc1(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp
       endif

       tempvar = pdbg
       call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass, &
                       nuc_bin,addt, ionrate, surf_area, BOXVOL, BOXMASS, &
                       TEMPTMS, PRES, RHTOMAS, PDBG, lev)

       if(pdbg) then
          print*,'COND_NUC: Error at nucleation[2] --> TERMINATE'
          errswitch=.true.
          return
       endif
       pdbg = tempvar
       if(pdbg) call debugprint(Nk2, Mk2, 0,0,0, 'After nucleation[2]')
       !print*,'after nucleation iter'
       sumfn = sumfn + fn*addt
       sumfn1 = sumfn1 + fn1*addt

       !total_nh4_1 = Mnuc1(srtnh4)
       !total_nh4_2 = Mnuc2(srtnh4)
       !do i=1,ibins
       !   total_nh4_1 = total_nh4_1 + Mk1(i,srtnh4)
       !   total_nh4_2 = total_nh4_2 + Mk2(i,srtnh4)
       !enddo
       !print*,'total_nh4',total_nh4_1,total_nh4_2

       mass_change = 0.e+0_fp

       do k=1,ibins
          mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
       enddo
       if( pdbg)    print*,'mass_change2',mass_change

       mcond = totmass-mass_change ! mass of h2so4 to condense

       !print*,'after nucleation'
       !print*,'totmass',totmass,'mass_change',mass_change,'mcond',mcond

       !print*,'2 mass_change',mass_change,mcond,totmass
       !print*,'2 cs1',CS1, Gc1(srtso4)

       if (mcond.lt.0.e+0_fp)then
          tmass = 0.e+0_fp
          do k=1,ibins
             do j=1,icomp-idiag
                tmass = tmass + Mk2(k,j)
             enddo
          enddo
          !if (abs(mcond).gt.tmass*1.0D-8) then
          if (abs(mcond).gt.totmass*1.0e-8_fp) then
             if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
                !if (CS1.gt.1.0D-5)then
                !   print*,'budget fudge 1 in cond_nuc'
                !endif
                tmass = 0.e+0_fp
                do j=1,icomp-idiag
                   tmass = tmass + Mk2(nuc_bin,j)
                enddo
                Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                mcond = 0.e+0_fp
             else
                print*,'budget fudge 2 in cond_nuc'
                do k=2,ibins
                   Nk2(k) = Nk1(k)
                   Mk2(k,srtso4) = Mk1(k,srtso4)
                enddo
                Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
                Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
                print*,'mcond < 0 in cond_nuc', mcond, totmass
                mcond = 0.e+0_fp
                ! should I stop or not?? (win, 10/4/08)
                !stop
                ! change from stop here to stop outside with more info (win, 10/4/08)
                print*,'COND_NUC: --> TERMINATE'
                !10/4/08 errswitch = .true.
                !10/4/08 return
             endif
          else
             mcond = 0.e+0_fp
          endif
       endif

       do k=1,ibins
          Nknuc(k) = Nknuc(k)+Nk2(k)-Nk1(k)
          do j=1,icomp-idiag
             Mknuc(k,j)=Mknuc(k,j)+Mk2(k,j)-Mk1(k,j)
          enddo
       enddo

       !Gc2(srtnh4) = Gc1(srtnh4)
       !call eznh3eqm(Gc2,Mk2,Mnuc2)
       !call ezwatereqm(Mk2,Mnuc2)

       !call getCondSink(Nk2,Mk2,Nnuc2,Mnuc2,srtso4,CStest,sinkfrac)

       ! Before entering ezcond, check if there's enough aerosol to
       ! condense onto. After several iteration in the case with high
       ! H2SO4 amount but little existing aerosol and also lack the conditions
       ! for nucleation, the whole size distribution is grown out of our
       ! tracked size bins, so let's exit the loop if there is no aerosol
       ! to condense onto anymore. (win, 10/4/08)
       tmass = 0.e+0_fp
       tnumb = 0.e+0_fp
       do k=1,ibins
          tnumb = tnumb + Nk2(k)
          do j=1,icomp-idiag
             tmass = tmass + Mk2(k,j)
          enddo
       enddo

       if( (tmass+mcond)/tnumb  > Xk(ibins) ) then
          if( .not. SPINUP(10.0) ) then
             print*,'Not enough aerosol for condensation!'
             print*,'  Exiting COND_NUC iteration with '
             print*,time_rem,'sec remaining time'
          endif

          Gc3(srtnh4)=Gc2(srtnh4)
          do k=1,ibins
             Nk3(k)=Nk2(k)
             do j=1,icomp
                Mk3(k,j)=Mk2(k,j)
             enddo
          enddo
          goto 100
       endif

       tempvar = pdbg

       call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area, &
                   BOXVOL, TEMPTMS, PRES, pdbg)
       do k=1,ibins
          Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
          do j=1,icomp-idiag
             Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
          enddo
       enddo
       Gc3(srtnh4) = Gc1(srtnh4)

       if(pdbg) then
          print*,'COND_NUC: Error at EZCOND[2] --> TERMINATE'
          errswitch=.true.
          return
       endif
       pdbg = tempvar

       if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[2]')

       if( pdbg)    print*,'after ezcond iter'
       call eznh3eqm(Gc3,Mk3)
       call ezwatereqm(Mk3, RHTOMAS)

       ! check to see how much condensation sink changed
       call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area, &
                        BOXVOL,TEMPTMS, PRES)

       time_rem = time_rem - addt
       if (time_rem .gt. 0.e+0_fp) then
          CSch = abs(CS2 - CS1)/CS1
          !jrp if (CSch.lt.0.e+0_fp) then
          !jrp    print*,''
          !jrp    print*,'CSch LESS THAN ZERO!!!!!', CS1,CStest,CS2
          !jrp    print*,'Nnuc',Nnuc1,Nnuc2
          !jrp    print*,''
          !jrp
          !jrp    addt = min(addt,time_rem)
          !jrp else

          ! Allow adaptive timestep to change
          ! Avoid div-by-zero error
          IF ( ABS( CSch ) > 0e+0_fp ) THEN
             addt = min(addt*CSch_tol/CSch,addt*1.5e+0_fp)
          ELSE
             addt = 0e+0_fp
          ENDIF

          ! allow adaptive timestep to change again
          addt = min(addt,time_rem)
          addt = max(addt,min_tstep)
          !jrp endif
          if( pdbg)     print*,'CS1',CS1,'CS2',CS2
          CS1 = CS2
          Gc1(srtnh4)=Gc3(srtnh4)
          do k=1,ibins
             Nk1(k)=Nk3(k)
             do j=1,icomp
                Mk1(k,j)=Mk3(k,j)
             enddo
          enddo
       endif
    enddo ! while loop

100 continue

    Gcf(srtso4)=sumH2SO4/dt
    fnavg = sumfn/dt
    fn1avg = sumfn1/dt
    if( pdbg)    print*,'AVERAGE GAS CONC',Gcf(srtso4)

    !jrp else
    !jrp    num_iter = 1
    !jrp    Gcf(srtso4)=Gc1(srtso4)
    !jrp endif

    if( pdbg) print*, 'cond_nuc num_iter =', num_iter
    !T0M(1,1,1,3) = double(num_iter) ! store iterations here

    if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'End of COND_NUC')

    do k=1,ibins
       Nkf(k)=Nk3(k)
       do j=1,icomp
          Mkf(k,j)=Mk3(k,j)
       enddo
    enddo
    Gcf(srtnh4)=Gc3(srtnh4)

    return

  END SUBROUTINE COND_NUC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getcondsink
!
! !DESCRIPTION: This subroutine calculates the condensation sink (first order
!  loss rate of condensing gases) from the aerosol size distribution.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II
!  Put in GEOS-Chem by Win T. (9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getCondSink(Nko, Mko, spec, CS, sinkfrac, surf_area, &
            BOXVOL, TEMPTMS, PRES)
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nk(ibins) - number of particles per size bin in grid cell
    !Nnuc - number of particles per size bin in grid cell
    !Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
    !Mk(ibins, icomp) - mass of a given species per size bin/grid cell
    !spec - number of the species we are finding the condensation sink for
    double precision Nko(ibins), Mko(ibins, icomp)
    REAL*4, INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
    integer spec
!
! !OUTPUT PARAMETERS:
!
    !CS - condensation sink [s^-1]
    !sinkfrac(ibins) - fraction of condensation sink from a bin
    double precision CS, sinkfrac(ibins)
    REAL(fp), INTENT(OUT)    :: surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision pi, R    ! pi and gas constant (J/mol K)
    double precision mu                  !viscosity of air (kg/m s)
    double precision mfp                 !mean free path of air molecule (m)
    double precision l_ab                !mean free path of h2so4 molecule (m)
    real Di       !diffusivity of gas in air (m2/s)
    double precision Neps     !tolerance for number
    real density  !density [kg m^-3]
    double precision mp       !mass per particle [kg]
    double precision Dpk(ibins) !diameter of particle [m]
    double precision Kn       !Knudson number
    double precision beta(ibins) !non-continuum correction factor
    double precision Mktot    !total mass in bin [kg]
    double precision c_a      !average speed of a, h2so4 molecule
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(Neps=1.0e+10_fp)
    double precision alpha(icomp) ! accomodation coef
    !data alpha/0.65,0.,0.,0.,0.,0.,0.,0.,0./
    real Sv(icomp)         !parameter used for estimating diffusivity
    !data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88, &
    !         42.88,42.88/

    !=================================================================
    ! getCondSink begins here
    !=================================================================

    ! have to find a better way to simply assign contants to these array
    ! The problem is I declare the array with ICOMP - its value will be
    ! determined at time of run, so I can't use DATA statement
    DO J=1,ICOMP
       !IF ( J == SRTSO4 ) THEN
       alpha(J) = 0.65
       !ELSE
       !   alpha(J) = 0.
       !ENDIF
       Sv(J) = 42.88
    ENDDO


    ! get some parameters

    !mu=2.5277e-7 * TEMPTMS**0.75302
    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

    Di=gasdiff(TEMPTMS,pres,98.0,Sv(spec))

    c_a  = sqrt(8.0 * TEMPTMS * R / 0.098)
    l_ab = 2.0 * Di / c_a

    ! get size dependent values
    do k=1,ibins
       if (Nko(k) .gt. Neps) then
          Mktot=0.e+0_fp
          do j=1,icomp
             Mktot=Mktot+Mko(k,j)
          enddo
          !kpc  Density should be changed due to more species involed.
          density=aerodens(Mko(k,srtso4),0.e+0_fp, &
                  Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil), &
                  Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob), &
                  Mko(k,srtdust),Mko(k,srth2o)) !assume bisulfate
          mp=Mktot/Nko(k)
       else
          !nothing in this bin - set to "typical value"
          density=1500.
#if  defined( TOMAS12 ) || defined( TOMAS15 )
          mp=sqrt(xk(k)*xk(k+1))
#else
          mp=1.4*xk(k)
#endif
       endif
       Dpk(k)  = ( (mp/density)*(6./pi) )**(0.333)
       !Kn     = 2.0 * mfp  / Dpk(k)     !S&P eqn 11.35 (text)
       Kn      = 2.0 * l_ab / Dpk(k)     !S&Pv2 chapter 12 - Kn for Dahneke correction factor
       beta(k) = ( 1.+Kn )  / ( 1.+2.*Kn*(1.+Kn)/alpha(spec) )   !S&P eqn 11.35
    enddo

    ! get condensation sink
    CS = 0.e+0_fp
    surf_area = 0.e+0_fp
    do k=1,ibins
       CS = CS + Dpk(k)*Nko(k)*beta(k)
       surf_area = surf_area+Nko(k)*pi*(Dpk(k)*1.0e+6_fp)**2
    enddo
    do k=1,ibins
       sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
    enddo
    CS = 2.e+0_fp*pi*dble(Di)*CS/(dble(boxvol)*1e-6_fp)
    surf_area = surf_area/(dble(boxvol))
    
    return

  end subroutine getcondsink
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getH2SO2conc
!
! !DESCRIPTION: This subroutine uses newtons method to solve for the steady
!  state H2SO4 concentration when nucleation is occuring.
!  It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)
!  where P is the production rate of H2SO4, CS is the condensation sink
!  and M(H2SO4) is the loss of mass towards making new particles.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II
!  Put in GEOS-CHEM by Win T. (9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getH2SO4conc(Nk, Mk, H2SO4rate, CS, NH3conc, gasConc, &
                          ionrate, surf_area, BOXVOL, BOXMASS, &
                          TEMPTMS, PRES, RHTOMAS, lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    ! H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
    ! CS - condensation sink [s-1]
    ! NH3conc - ammonium in box [kg box-1]
    REAL(fp)            :: Nk(IBINS)
    REAL(fp)            :: Mk(IBINS, ICOMP)
    double precision       H2SO4rate
    double precision       CS
    double precision       NH3conc
    REAL*4, INTENT(IN)  :: BOXVOL,  BOXMASS, TEMPTMS
    REAL*4, INTENT(IN)  :: PRES,    RHTOMAS
    integer                lev
!
! !OUTPUT PARAMETERS:
!
    ! gasConc - gas H2SO4 [kg/box]
    double precision       gasConc
    REAL(fp)            :: ionrate
    REAL(fp)            :: surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
    double precision mnuc, mnuc1 ! mass of nucleated particle [kg]
    double precision fn1, rnuc1 ! nucleation rate [# cm-3 s-1] and critical radius [nm]
    double precision res      ! d[H2SO4]/dt, need to find the solution where res = 0
    double precision massnuc     ! mass being removed by nucleation [kg s-1 box-1]
    double precision gasConc1 ! perturbed gasConc
    double precision gasConc_hi, gasConc_lo
    double precision res1     ! perturbed res
    double precision res_new  ! new guess for res
    double precision dresdgasConc ! derivative for newtons method
    double precision Gci(icomp-1)      !array to carry gas concentrations
    logical nflg              !says if nucleation occured
    double precision H2SO4min !minimum H2SO4 concentration in parameterizations (molec/cm3)
    double precision pi
    integer iter,iter1
    double precision CSeps    ! low limit for CS
    double precision max_H2SO4conc !maximum H2SO4 concentration in parameterizations (kg/box)
    double precision nh3ppt   !ammonia concentration in ppt
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654)
    !parameter(H2SO4min=1.D4) !molecules cm-3
    parameter(CSeps=1.0e-20_fp)

    !=================================================================
    ! getH2SO4conc begins here
    !=================================================================

    do i=1,icomp-1
       Gci(i)=0.e+0_fp
    enddo
    Gci(srtnh4)=NH3conc

    ! make sure CS doesn't equal zero
    !CS = max(CS,CSeps)

    ! some specific stuff for napari vs. vehk
    if (ion_nuc.eq.1) then
       H2SO4min=1.0e+5_fp
    elseif (ion_nuc.eq.2) then
       H2SO4min=5.0e+5_fp
    else
       H2SO4min=1.0e+4_fp
    endif

    if ((bin_nuc.eq.1).or.(tern_nuc.eq.1).or.(ion_nuc.le.2))then
       nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
  &             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)
       if (ion_nuc.eq.1)then
          max_H2SO4conc=1.0e+8_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif (ion_nuc.eq.2)then
          max_H2SO4conc=5.0e+8_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif ((nh3ppt.gt.1.0e+0_fp).and.(tern_nuc.eq.1))then
          max_H2SO4conc=1.0e+9_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif (bin_nuc.eq.1)then
          max_H2SO4conc=1.0e+11_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       else
          max_H2SO4conc = 1.0e+100_fp
       endif
    else
       max_H2SO4conc = 1.0e+100_fp
    endif

    ! Checks for when condensation sink is very small
    if (CS.gt.CSeps) then
       gasConc = H2SO4rate/CS
    else
       if((bin_nuc.gt.0).or.(tern_nuc.gt.0).or. (ion_nuc.gt.0))then
          gasConc = max_H2SO4conc
       else
          print*,'condesation sink too small in getH2SO4conc'
          STOP
       endif
    endif

    gasConc = min(gasConc,max_H2SO4conc)
    Gci(srtso4) = gasConc
    call getNucRate(Nk, Mk, Gci,fn,mnuc,nflg,ionrate, surf_area, &
                    BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)

    if (fn.gt.0.e+0_fp) then      ! nucleation occured
       !convert to kg/box
       gasConc_lo = H2SO4min*boxvol/(1000.e+0_fp/98.e+0_fp*6.022e+23_fp)

       ! Test to see if gasConc_lo gives a res < 0
       ! (this means ANY nucleation is too high)
       Gci(srtso4) = gasConc_lo*1.000001e+0_fp
       call getNucRate(Nk,Mk,Gci,fn1,mnuc1,nflg,ionrate,surf_area, &
                       BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
       if (fn1.gt.0.e+0_fp) then
          massnuc = mnuc1*fn1*boxvol*98.e+0_fp/96.e+0_fp
          !massnuc = 4.e+0_fp/3.e+0_fp*pi*(rnuc1*1.e-9_fp)**3*1350.*fn1*boxvol*
          !massnuc = 4.e+0_fp/3.e+0_fp*pi*(rnuc1*1.e-9_fp)**3*1800.*fn1*boxvol*%
          !          98.e+0_fp/96.e+0_fp
          !jrp print*,'res',res
          !jrp print*,'H2SO4rate',H2SO4rate
          !jrp print*,'CS*gasConc_lo',CS*gasConc_lo
          !jrp print*,'mnuc',mnuc
          res = H2SO4rate - CS*gasConc_lo - massnuc
          if (res.lt.0.e+0_fp) then ! any nucleation too high
             ! if (.not. spinup(14.0)) print*,'nucleation cuttoff'
             ! have nucleation occur and fix mass balance after
             gasConc = gasConc_lo*1.000001
             return
          endif
       endif

       ! we know this must be the upper limit (since no nucleation)
       gasConc_hi = gasConc
       !take density of nucleated particle to be 1350 kg/m3
       massnuc = mnuc*fn*boxvol*98.e+0_fp/96.e+0_fp
       !print*,'H2SO4rate',H2SO4rate,'CS*gasConc',CS*gasConc,'mnuc',mnuc
       res = H2SO4rate - CS*gasConc - massnuc

       ! check to make sure that we can get solution
       if (res.gt.H2SO4rate*1.e-10_fp) then
          print*,'gas production rate too high in getH2SO4conc'
          print*,H2SO4rate,CS,gasConc,massnuc,res
          return
          !STOP
       endif

       iter = 0
       !jrp print*, 'iter',iter
       !jrp print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
       !jrp print*,'res',res
       do while ((abs(res/H2SO4rate).gt.1.e-4_fp).and.(iter.lt.40))
          iter = iter+1
          if (res .lt. 0.e+0_fp) then ! H2SO4 concentration too high, must reduce
             gasConc_hi = gasConc ! old guess is new upper bound
          elseif (res .gt. 0.e+0_fp) then ! H2SO4 concentration too low, must increase
             gasConc_lo = gasConc ! old guess is new lower bound
          endif
          !print*, 'iter',iter
          !print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
          gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
          Gci(srtso4) = gasConc
          call getNucRate(Nk, Mk,Gci,fn,mnuc,nflg,ionrate,surf_area, &
                          BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
          massnuc = mnuc*fn*boxvol*98.e+0_fp/96.e+0_fp
          res = H2SO4rate - CS*gasConc - massnuc
          !print*,'res',res
          !print*,'H2SO4rate',H2SO4rate,'CS',CS,'gasConc',gasConc
          if (iter.eq.40.and.CS.gt.1.0e-4_fp)then
             print*,'getH2SO4conc iter break'
             print*,'H2SO4rate',H2SO4rate,'CS',CS
             print*,'gasConc',gasConc,'massnuc',massnuc
             print*,'max_H2SO4conc',max_H2SO4conc
             print*,'fn',fn
             print*,'res/H2SO4rate',res/H2SO4rate
          endif
       enddo

       !print*,'IN getH2SO4conc'
       !print*,'fn',fn
       !print*,'H2SO4rate',H2SO4rate
       !print*,'massnuc',massnuc,'CS*gasConc',CS*gasConc

    else  
       ! nucleation didn't occur
    endif

    return

  end SUBROUTINE GETH2SO4CONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getnucrate
!
! !DESCRIPTION: This subroutine calls the Vehkamaki 2002 and Napari 2002
!  nucleation parameterizations and gets the binary and ternary nucleation
!  rates. 
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II
!  Put in GEOS-Chem by win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getNucRate(Nk, Mk, Gci,fn,mnuc,nflg, ionrate,surf_area, &
                        BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    ! Gci(icomp-1) - amount (kg/grid cell) of all species present in the
    !                gas phase except water
    REAL*4,   INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
    REAL*4,   INTENT(IN)       :: PRES,    RHTOMAS
    REAL(fp), INTENT(IN)       :: Gci(icomp-1)
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(INOUT)    :: Nk(IBINS)
    REAL(fp), INTENT(INOUT)    :: Mk(IBINS, ICOMP)
!
! !OUTPUT PARAMETERS:
!
    ! fn - nucleation rate [# cm-3 s-1]
    ! rnuc - radius of nuclei [nm]
    ! nflg - says if nucleation happend
    REAL(fp)                   :: surf_area
    REAL(fp)                   :: ionrate

    integer j,i,k
    double precision fn       ! nucleation rate to first bin cm-3 s-1
    double precision mnuc     !mass of nucleating particle [kg]
    logical nflg
    integer lev
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision nh3ppt   ! gas phase ammonia in pptv
    double precision h2so4    ! gas phase h2so4 in molec cc-1
    double precision gtime    ! time to grow to first size bin [s]
    double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
    double precision Mktot    ! total mass in bin
    double precision neps
    double precision meps
    double precision density  ! density of particle [kg/m3]
    double precision pi
    double precision frac     ! fraction of particles growing into first size bin
    double precision d1,d2    ! diameters of particles [m]
    double precision mp       ! mass of particle [kg]
    double precision mold     ! saved mass in first bin
    double precision rnuc     ! critical nucleation radius [nm]
    double precision sinkfrac(ibins) ! fraction of loss to different size bins
    double precision nadd     ! number to add
    double precision CS       ! kerminan condensation sink [m-2]
    double precision Dpmean   ! the number wet mean diameter of the existing aerosol
    double precision Dp1      ! the wet diameter of bin 1
    double precision dens1    ! density in bin 1 [kg m-3]
    double precision GR       ! growth rate [nm hr-1]
    double precision gamma,eta ! used in kerminen 2004 parameterzation
    double precision drymass,wetmass,WR
    double precision fn_c     ! barrierless nucleation rate
    double precision h1,h2,h3,h4,h5,h6
    double precision dum1,dum2,dum3,dum4   ! dummy variables
    double precision rhin,tempin ! rel hum in

    real(fp)    mydummy
!
! !DEFINED PARAMETERS:
!
    parameter (neps=1E8, meps=1E-8)
    parameter (pi=3.14159)

    !=================================================================
    ! getNucRate begins here
    !=================================================================

    h2so4 = Gci(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
    nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

    fn = 0.e+0_fp
    rnuc = 0.e+0_fp

    !print*,'h2so4',h2so4,'nh3ppt',nh3ppt

    ! if requirements for nucleation are met, call nucleation subroutines
    ! and get the nucleation rate and critical cluster size
    if (h2so4.gt.1.e+4_fp) then
       if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then
          ! print*, 'napari'
          call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc
          if (ion_nuc.eq.1.and.ionrate.ge.1.e+0_fp) then
             call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                           h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          nflg=.true.
       elseif (bin_nuc.eq.1) then
          ! print*, 'vehk'
          call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
          if ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
             call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                           h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          if (fn.gt.1.0e-6_fp)then
             nflg=.true.
          else
             fn = 0.e+0_fp
             nflg=.false.
          endif
       elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
          call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                        h1,h2,h3,h4,h5,h6)
          fn=h1
          rnuc=h5
          nflg=.true.
       elseif(ion_nuc.eq.2) then
          ! Yu Ion nucleation
          !! first we need to calculate the available surface area
          !surf_area = 0.e+0_fp
          !do k=1, ibins
          !   if (Nki(k) .gt. Neps) then
          !      Mktot=0.e+0_fp
          !      do j=1,icomp
          !         Mktot=Mktot+Mki(k,j)
          !      enddo
          !      mp=Mktot/Nki(k)
          !      density=aerodens(Mki(k,srtso4),0.e+0_fp, &
          !                       Mki(k,srtnh4),0.e+0_fp,Mki(k,srth2o))  ! assume bisulfate
          !      ! diameter = ((mass/density)*(6/pi))**(1/3)
          !      d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
          !      ! surface area per particle = pi*diameter**2
          !      surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)* &
          !                  pi*(d2**2.D0) ! (um2 cm-2)
          !   endif
          !enddo
          rhin=dble(RHTOMAS*100.e+0_fp)
          tempin=dble(TEMPTMS)
          call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area, &
                      fn, dum1, rnuc, dum2)
          nflg=.true.
       else
          nflg=.false.
       endif
       if((act_nuc.eq.1).and.(lev.le.7))then
          call bl_nucl(h2so4,fn,rnuc)
          nflg=.true.
       endif
       call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max for ternary
       fn = min(fn,fn_c)
    else
       nflg=.false.
    endif

    if (fn.gt.0.e+0_fp) then
       call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1, &
                             BOXVOL, TEMPTMS, PRES)
       d1 = rnuc*2.e+0_fp*1e-9_fp
       drymass = 0.e+0_fp
       do j=1,icomp-idiag
          drymass = drymass + Mk(1,j)
       enddo
       wetmass = 0.e+0_fp
       do j=1,icomp
          wetmass = wetmass + Mk(1,j)
       enddo
       !prior 10/15/08
       !WR = wetmass/drymass

       ! prevent division by zero (win, 10/15/08)
       if( drymass == 0.e+0_fp ) then
          WR = 1.e+0_fp
       else
          WR = wetmass/drymass
       endif

       !print*,'[getnucrate] Gci',Gci
       !print*,'WR',WR, 'drymass',drymass, 'wetmass',wetmass
       call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS, &
                          boxvol,dens1,gtime)
       GR = (Dp1-d1)*1e+9_fp/gtime*3600.e+0_fp ! growth rate, nm hr-1

       gamma = 0.23e+0_fp*(d1*1.0e+9_fp)**(0.2e+0_fp)* &
               (Dp1*1.0e+9_fp/3.e+0_fp)**0.075e+0_fp* &
               (Dpmean*1.0e+9_fp/150.e+0_fp)** &
               0.048e+0_fp*(dens1*1.0e-3_fp)** &
               (-0.33e+0_fp)*(TEMPTMS/293.e+0_fp) ! equation 5 in kerminen
       eta = gamma*CS/GR
       !print*,'fn1',fn
       if (Dp1.gt.d1)then
          fn = fn*exp(eta/(Dp1*1.0e+9_fp)-eta/(d1*1.0e+9_fp))
       endif
       !print*,'fn2',fn
       if( IT_IS_NAN( fn ) ) then
          print*, '---------------->>> Found NAN in GETNUCRATE'
          print*,'fn',fn
          print*,'eta',eta, 'Dp1',Dp1,'d1',d1
          print*,'gamma',gamma,'CS',CS,'GR',GR,'gtime',gtime
          call ERROR_STOP('Found NaN in fn','getnucrate')
       endif

       mnuc = sqrt(xk(1)*xk(2))
    endif

    return

  end SUBROUTINE GETNUCRATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vehk_nucl
!
! !DESCRIPTION: Subroutine vehk_nucl calculates the binary nucleation rate and
!  radius of the critical nucleation cluster using the parameterization of...
!  .
!    Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck,
!    M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric
!    Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions."
!    Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).
!  .
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak Sep 29,2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VEHK_NUCL (tempi,rhi,cnai,fn,rnuc)
!
! !INPUT PARAMETERS:
!
    real*4,   intent(in)   :: tempi ! temperature of air [K]
    real*4,   intent(in)   :: rhi ! relative humidity of air as a fraction
    real(fp), intent(in)   :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out)  :: fn ! nucleation rate [cm-3 s-1]
    real(fp), intent(out)  :: rnuc ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  :: fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
    REAL(fp)  :: gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
    REAL(fp)  :: temp    ! temperature of air [K]
    REAL(fp)  :: rh      ! relative humidity of air as a fraction
    REAL(fp)  :: cna     ! concentration of gas phase sulfuric acid [molec cm-3]
    REAL(fp)  :: xstar   ! mole fraction sulfuric acid in cluster
    REAL(fp)  :: ntot    ! total number of molecules in cluster
    integer   :: i       ! counter

    ! Nucleation Rate Coefficients
    data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598, &
              2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
    data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796, &
              0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
    data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, &
               0.00864245, -0.000407382, 0.0000377124, 0.000157098, &
              -0.0000198417, -9.24619e-6_fp /
    data fb3 /0.0000722811, 0.0000404196, &
             -4.7758e-6_fp, 3.10801e-7_fp, &
             -0.0000228947, -4.01957e-7_fp, &
              3.21794e-8_fp, 4.00914e-7_fp, &
              2.46048e-8_fp, 5.00427e-9_fp /
    data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985, &
              0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  /

    ! Coefficients of total number of molecules in cluster
    data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323, &
              -0.0125211, -0.038546, -0.0183749, -0.0619974, &
               0.0121827, 0.000320184 /
    data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104, &
               0.00580655, -0.000672316, 0.000172072, 0.000906958, &
              -0.00010665, -0.0000174762 /
    data gb2 /0.00102485, 0.000192654, &
             -0.0000154757, 2.65066e-6_fp, &
             -0.000101674, 2.60288e-6_fp, &
             -3.71766e-7_fp, -9.11728e-7_fp, &
             2.5346e-7_fp, 6.06504e-8_fp /
    data gb3 /-2.18646e-6_fp, -6.7043e-7_fp, &
               5.66661e-8_fp, -3.67471e-9_fp, &
               2.88195e-7_fp, 1.19416e-8_fp, &
              -5.14875e-10_fp, -5.36796e-9_fp, &
              -3.63519e-10_fp, -1.42177e-11_fp /
    data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251, &
               0.0942243, -0.00851515, 0.00026866, -0.00774234, &
               0.000610065, 0.000135751 /

    !=================================================================
    ! VEHK_NUCL begins here!
    !=================================================================
    temp=dble(tempi)
    rh=dble(rhi)
    cna=cnai

    ! Respect the limits of the parameterization
    if (cna .lt. 1.e4_fp) then ! limit sulf acid conc
       fn = 0.
       rnuc = 1.
       !print*,'cna < 1D4', cna
       goto 10
    endif
    if (cna .gt. 1.0e+11_fp) cna=1.0e11 ! limit sulfuric acid conc
    if (temp .lt. 230.15) temp=230.15 ! limit temp
    if (temp .gt. 305.15) temp=305.15 ! limit temp
    if (rh .lt. 1e-4_fp) rh=1e-4_fp ! limit rh
    if (rh .gt. 1.) rh=1. ! limit rh

    ! Mole fraction of sulfuric acid
    xstar=0.740997-0.00266379*temp-0.00349998*log(cna) &
         +0.0000504022*temp*log(cna)+0.00201048*log(rh) &
         -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2. &
         -0.0000179059*temp*(log(rh))**2. &
         +0.000184403*(log(rh))**3. &
         -1.50345e-6_fp*temp*(log(rh))**3.

    ! Nucleation rate coefficients
    do i=1, 10
       fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2. &
              +fb3(i)*temp**3.+fb4(i)/xstar
    enddo

    ! Nucleation rate (1/cm3-s)
    fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2. &
         +fb(4)*(log(rh))**3.+fb(5)*log(cna) &
         +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna) &
         +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2. &
         +fb(10)*(log(cna))**3.)

    !print*,'in vehk_nuc, fn',fn
    !print*,'cna',cna,'rh',rh,'temp',temp
    !print*,'xstar',xstar

    ! Cap at 10^6 particles/s, limit for parameterization
    if (fn.gt.1.0e+6_fp) then
       fn=1.0e+6_fp
    endif

    ! Coefficients of total number of molecules in cluster
    do i=1, 10
       gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2. &
              +gb3(i)*temp**3.+gb4(i)/xstar
    enddo
    ! Total number of molecules in cluster
    ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2. &
         +gb(4)*(log(rh))**3.+gb(5)*log(cna) &
         +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna) &
         +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2. &
         +gb(10)*(log(cna))**3.)

    ! cluster radius
    rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

10  return

  end SUBROUTINE VEHK_NUCL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bl_nucl
!
! !DESCRIPTION: This subroutine calculates a simple binary nucleation rate of
!  1 nm particles.
!  WRITTEN BY Jeff Pierce, April 2007
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE bl_nucl(cnai,fn,rnuc)
!
! !INPUT PARAMETERS:
!
    ! concentration of gas phase sulfuric acid [molec cm-3]
    double precision cnai
!
! !OUTPUT PARAMETERS:
!
    double precision fn                   ! nucleation rate [cm-3 s-1]
    double precision rnuc                 ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision cna ! concentration of gas phase sulfuric acid [molec cm-3]
    double precision A   ! prefactor... empirical
!
! !DEFINED PARAMETERS:
!
    parameter(A=2.0e-6_fp)

    !=================================================================
    ! bl_nucl begins here
    !=================================================================

    cna=cnai

    fn=A*cna
    rnuc=0.5e+0_fp ! particle diameter of 1 nm

10  return

  end SUBROUTINE bl_nucl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: napa_nucl
!
! !DESCRIPTION:  Subroutine NAPA_NUCL calculates the ternary nucleation rate
!  and radius of the critical nucleation cluster using the parameterization of
!  .
!     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of
!     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical
!     Research-Atmospheres 107, no. D19 (2002).
!  .
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak Sep 29, 2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE napa_nucl(tempi,rhi,cnai,nh3ppti,fn,rnuc)
!
! !INPUT PARAMETERS:
!
    real*4,   intent(in) :: tempi ! temperature of air [K]
    real*4,   intent(in) :: rhi ! relative humidity of air as a fraction
    real(fp), intent(in) :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
    real(fp), intent(in) :: nh3ppti ! concentration of gas phase ammonia
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out):: fn  ! nucleation rate [cm-3 s-1]
    real(fp), intent(out):: rnuc ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)    ::  aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
    real(fp)    ::  fnl     ! natural log of nucleation rate
    real(fp)    ::  temp    ! temperature of air [K]
    real(fp)    ::  rh      ! relative humidity of air as a fraction
    real(fp)    ::  cna     ! concentration of gas phase sulfuric acid [molec cm-3]
    real(fp)    ::  nh3ppt  ! concentration of gas phase ammonia
    integer     ::  i       ! counter

    ! Adjustable parameters
    data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916, &
               0.31176, -0.0200738, 0.165536, &
               6.52645, 3.68024, -0.066514, 0.65874, &
               0.0599321, -0.732731, 0.728429, 41.3016, &
               -0.160336, 8.57868, 0.0530167, -2.32736        /

    data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378, &
               1.64009, -0.752115, 3.26623, -0.258002, -0.204098, &
              -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736, &
              -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/

    data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052, &
             0.00939836, -0.00343852, 0.00525813, -0.0489703, &
             0.00143456, 0.00106259, 0.0122938, -0.00165718, &
            -0.0362432, 0.000147186, -0.027422, 0.000904383, &
            -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/

    data a3 /-0.000824007, -9.74576e-06_fp, &
             -4.14699e-07_fp, 7.61229e-06_fp, &
             -1.74927e-05_fp, -1.09753e-05_fp, &
             -8.98038e-06_fp, 0.000146967, &
             -2.02036e-06_fp, -1.2656e-06_fp, &
              6.18554e-05_fp, 3.41744e-06_fp, &
              4.93337e-05_fp, -2.37711e-07_fp, &
              4.93478e-05_fp, -5.73788e-07_fp, &
              8.39522e-08_fp, -6.48365e-07_fp, &
             -2.93564e-05_fp, 8.0459e-08_fp   /

    !=================================================================
    ! NAPA_NUCL begins here!
    !=================================================================
    temp=dble(tempi)
    rh=dble(rhi)
    cna=cnai
    nh3ppt=nh3ppti

    ! Napari's parameterization is only valid within limited area
    if ((cna .lt. 1.e+4_fp).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
       fn = 0.
       rnuc = 1
       goto 10
    endif
    if (cna .gt. 1.0e+9_fp) cna=1.0e+9_fp ! limit sulfuric acid conc
    if (nh3ppt .gt. 100.) nh3ppt=100. ! limit temp
    if (temp .lt. 240.) temp=240. ! limit temp
    if (temp .gt. 300.) temp=300. ! limit temp
    if (rh .lt. 0.05) rh=0.05 ! limit rh
    if (rh .gt. 0.95) rh=0.95 ! limit rh

    do i=1,20
       fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
    enddo

    fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2. &
       +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh &
       +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt) &
       *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna) &
       +fa(12)*rh &
       *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh) &
       *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna) &
       *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt) &
       +fa(18)*rh &
       *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20) &
       *(log(cna))**2.*(log(nh3ppt))**2.

    fn=exp(fnl)

    ! Try scaling down the rate by 1e-5 to see how the param is
    ! doing on the false positive nucleation (win, 12/18/08)
    !sensitivity simulation, change scaling factor down to 1e-4
    fn = fn * 1.e-5

    ! Cap at 10^6 particles/cm3-s, limit for parameterization
    if (fn.gt.1.0e+6_fp) then
       fn=1.0e+6_fp
       fnl=log(fn)
    endif
      
    rnuc=0.141027-0.00122625*fnl-7.82211e-6_fp*fnl**2. &
        -0.00156727*temp-0.00003076*temp*fnl &
        +0.0000108375*temp**2.

10  return

  end subroutine napa_nucl
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getcondsink
!
! !DESCRIPTION: Subroutine GETCONDSINK\_KERM calculates the condensation sink
!  (first order loss rate of condensing gases) from the aerosol size
!  distribution.
!  .
!  This is the cond sink in kerminen et al 2004 Parameterization for
!  new particle formation AS&T Eqn 6.
!  .
!  Written by Jeff Pierce, May 2007 for GISS GCM-II'
!  Introduced to GEOS-Chem by Win Trivitayanurak, Sep 29, 2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getCondSink_kerm(Nko,Mko,CS,Dpmean,Dp1,dens1, &
                              BOXVOL, TEMPTMS, PRES)
!
! !INPUT PARAMETERS:
!
    ! Nk(ibins) - number of particles per size bin in grid cell
    ! Mk(ibins, icomp) - mass of a given species per size bin/grid cell
    REAL(fp), INTENT(IN)        :: Nko(ibins), Mko(ibins, icomp)
    REAL*4,   INTENT(IN)        :: BOXVOL, TEMPTMS, PRES
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT)       :: CS       ! CS - condensation sink [s^-1]
    REAL(fp), INTENT(OUT)       :: Dpmean   ! the number mean diameter [m]
    REAL(fp), INTENT(OUT)       :: Dp1      ! the size of the first size bin [m]
    REAL(fp), INTENT(OUT)       :: dens1    ! the density of the first size bin [kg/m3]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Nnuc - number of particles per size bin in grid cell
    ! Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
    ! spec - number of the species we are finding the condensation sink for
    ! sinkfrac(ibins) - fraction of condensation sink from a bin
    integer        :: i,j,k,c           ! counters
    REAL(fp)       :: pi, R    ! pi and gas constant (J/mol K)
    REAL(fp)       :: mu                  !viscosity of air (kg/m s)
    REAL(fp)       :: mfp                 !mean free path of air molecule (m)
    REAL*4         :: Di       !diffusivity of gas in air (m2/s)
    REAL(fp)       :: Neps     !tolerance for number
    REAL*4         :: density  !density [kg m^-3]
    REAL(fp)       :: mp       !mass per particle [kg]
    REAL(fp)       :: Dpk(ibins) !diameter of particle [m]
    REAL(fp)       :: Kn       !Knudson number
    REAL(fp)       :: beta(ibins) !non-continuum correction factor
    REAL(fp)       :: Mktot    !total mass in bin [kg]
    REAL(fp)       :: Dtot,Ntot ! used on getting the number mean diameter
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(Neps=1.0e+10_fp)

    !=================================================================
    ! GETCONDSINK_KERM  begins here!
    !=================================================================

    ! get some parameters
    mu=2.5277e-7*TEMPTMS**0.75302
    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
    mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
    !Di=gasdiff(temp,pres,98.0,Sv(srtso4))
    !print*,'Di',Di

    ! get size dependent values
    CS = 0.e+0_fp
    Ntot = 0.e+0_fp
    Dtot = 0.e+0_fp
    do k=1,ibins
       if (Nko(k) .gt. Neps) then
          Mktot=0.e+0_fp
          do j=1,icomp
             Mktot=Mktot+Mko(k,j)
          enddo
          !kpc Density should be changed due to more species involed.
          density=aerodens(Mko(k,srtso4),0.e+0_fp, &
                  Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil), &
                  Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob), &
                  Mko(k,srtdust),Mko(k,srth2o))
          mp=Mktot/Nko(k)
       else
          !nothing in this bin - set to "typical value"
          density=1500.
          mp=1.4*xk(k)
       endif
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)      !S&P eqn 11.35 (text)
       CS=CS+0.5e+0_fp*(Dpk(k)*Nko(k)/(dble(boxvol)*1.0e-6_fp)*(1+Kn)) &
            /(1.e+0_fp+0.377e+0_fp*Kn+1.33e+0_fp*Kn*(1+Kn))
       Ntot = Ntot + Nko(k)
       Dtot = Dtot + Nko(k)*Dpk(k)
       if (k.eq.1)then
          Dp1=Dpk(k)
          dens1 = density
       endif
    enddo

    if (Ntot.gt.1e+15_fp)then
       Dpmean = Dtot/Ntot
    else
       Dpmean = 150.e+0_fp
    endif

    return

  END SUBROUTINE GETCONDSINK_KERM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getgrowthtime
!
! !DESCRIPTION: This subroutine calculates the time it takes for a particle to
!  grow from one size to the next by condensation of sulfuric acid (and
!  associated NH3 and water) onto particles.
!  .
!  This subroutine assumes that the growth happens entirely in the kinetic
!  regine such that the dDp/dt is not size dependent.  The time for growth
!  to the first size bin may then be approximated by the time for growth via
!  sulfuric acid (not including nh4 and water) to the size of the first size bin
!  (not including nh4 and water).
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/29/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getGrowthTime (d1,d2,h2so4,temp,boxvol,density,gtime)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    ! d1: intial diameter [m]
    ! d2: final diameter [m]
    ! h2so4: h2so4 ammount [kg]
    ! temp: temperature [K]
    ! boxvol: box volume [cm3]
    REAL(fp), INTENT(IN)  ::  d1,d2    ! initial and final diameters [m]
    REAL(fp), INTENT(IN)  ::  h2so4    ! h2so4 amount [kg]
    real*4,   INTENT(IN)  ::  temp     ! temperature [K]
    real*4,   INTENT(IN)  ::  boxvol  ! box volume [cm3]
    REAL(fp), INTENT(IN)  ::  density  ! density of particles in first bin [kg/m3]
!
! !OUTPUT PARAMETERS:
!
    ! gtime: the time it takes the particle to grow to first size bin [s]
    REAL(fp), INTENT(OUT) ::  gtime    ! the time it will take the particle to
                                       ! grow to first size bin [s]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)     ::  pi, R, MW
    REAL(fp)     ::  csulf    ! concentration of sulf acid [kmol/m3]
    REAL(fp)     ::  mspeed   ! mean speed of molecules [m/s]
    REAL(fp)     ::  alpha    ! accomidation coef
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654e+0_fp, R=8.314e+0_fp) !pi and gas constant (J/mol K)
    parameter(MW=98.e+0_fp) ! density [kg/m3], mol wgt sulf [kg/kmol]
    parameter(alpha=0.65)

    !=================================================================
    ! GETGROWTHTIME begins here!
    !=================================================================
    !print *,'h2so4',h2so4,'MW',MW,'boxvol',boxvol,dble(boxvol)

    csulf = h2so4/MW/(dble(boxvol)*1e-6_fp) ! SA conc. [kmol/m3]
    mspeed = sqrt(8.e+0_fp*R*dble(temp)*1000.e+0_fp/(pi*MW))

    ! Kinetic regime expression (S&P 11.25) solved for T
    gtime = (d2-d1)/(4.e+0_fp*MW/density*mspeed*alpha*csulf)

    if ( IT_IS_NAN(gtime) ) then
       !jrp
       print*,'IN GET GROWTH TIME'
       print*,'d1',d1,'d2',d2
       print*,'h2so4',h2so4
       print*,'boxvol',boxvol
       print*,'csulf',csulf,'mspeed',mspeed
       print*,'density',density,'gtime',gtime
       call ERROR_STOP('Found NaN in fn','getnucrate')
    endif

    RETURN

  END SUBROUTINE GETGROWTHTIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nucleation
!
! !DESCRIPTION: This subroutine calls the Vehkamaki 2002 and Napari 2002
!  nucleation parameterizations and gets the binary and ternary nucleation
!  rates. The number of particles added to the first size bin is calculated
!  by comparing the growth rate of the new particles to the coagulation sink.
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NUCLEATION(Nki,Mki,Gci,Nkf,Mkf,Gcf,fn,fn1,totsulf, &
                        nuc_bin,dt,ionrate, surf_area, BOXVOL, BOXMASS, &
                        TEMPTMS, PRES, RHTOMAS, pdbg,lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nki(ibins) - number of particles per size bin in grid cell
    !Mki(ibins, icomp) - mass of a given species per size bin/grid cell
    !Gci(icomp-1) - amount (kg/grid cell) of all species present in the
    !               gas phase except water
    !dt - total model time step to be taken (s)
    double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
    REAL*4, INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
    REAL*4, INTENT(IN)       :: PRES,    RHTOMAS
!
! !OUTPUT PARAMETERS:
!
    !Nkf, Mkf, Gcf - same as above, but final values
    !fn, fn1
    double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
    integer j,i,k
    double precision totsulf
    integer nuc_bin
    double precision dt
    double precision fn       ! nucleation rate of clusters cm-3 s-1
    double precision fn1      ! formation rate of particles to first size bin cm-3 s-1

    LOGICAL  PDBG             ! Signal print for debug
    integer lev ! layer of model

    REAL(fp)                     ionrate
    REAL(fp)                     surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision nh3ppt   ! gas phase ammonia in pptv
    double precision h2so4    ! gas phase h2so4 in molec cc-1
    double precision rnuc     ! critical nucleation radius [nm]
    double precision gtime    ! time to grow to first size bin [s]
    double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
    double precision Mktot    ! total mass in bin
    double precision neps
    double precision meps
    double precision density  ! density of particle [kg/m3]
    double precision pi
    double precision frac     ! fraction of particles growing into first size bin
    double precision d1,d2    ! diameters of particles [m]
    double precision mp       ! mass of particle [kg]
    double precision mold     ! saved mass in first bin
    double precision mnuc     !mass of nucleation
    double precision sinkfrac(ibins) ! fraction of loss to different size bins
    double precision nadd     ! number to add
    double precision CS       ! kerminan condensation sink [m-2]
    double precision Dpmean   ! the number wet mean diameter of the existing aerosol
    double precision Dp1      ! the wet diameter of bin 1
    double precision dens1    ! density in bin 1 [kg m-3]
    double precision GR       ! growth rate [nm hr-1]
    double precision gamma,eta ! used in kerminen 2004 parameterzation
    double precision drymass,wetmass,WR
    double precision fn_c     ! barrierless nucleation rate
    double precision h1,h2,h3,h4,h5,h6
    double precision dum1,dum2,dum3,dum4   ! dummy variables
    double precision rhin,tempin ! rel hum in

    LOGICAL ERRORSWITCH
!
! !DEFINED PARAMETERS:
!
    parameter (neps=1E8, meps=1E-8)
    parameter (pi=3.14159)

    !=================================================================
    ! NUCLEATION begins here
    !=================================================================

    errorswitch = .false.

    h2so4 = Gci(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
    nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

    fn = 0.e+0_fp
    fn1 = 0.e+0_fp
    rnuc = 0.e+0_fp
    gtime = 0.e+0_fp
    nuc_bin = 1 ! added by Pengfei Liu,initialize  nuc_bin value
    ! if requirements for nucleation are met, call nucleation subroutines
    ! and get the nucleation rate and critical cluster size
    if (h2so4.gt.1.e+4_fp) then
       if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
          call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc
          if (ion_nuc.eq.1.and.ionrate.ge.1.e+0_fp) then
             call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                           h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
       elseif (bin_nuc.eq.1) then
          call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
          if ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
             call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                           h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          if (fn.lt.1.0e-6_fp)then
             fn = 0.e+0_fp
          endif
       elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
          call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
                        h1,h2,h3,h4,h5,h6)
          fn=h1
          rnuc=h5
       elseif(ion_nuc.eq.2) then
          ! Yu Ion nucleation
          !! first we need to calculate the available surface area
          !surf_area = 0.e+0_fp
          !do k=1, ibins
          !   if (Nki(k) .gt. Neps) then
          !      Mktot=0.e+0_fp
          !      do j=1,icomp
          !         Mktot=Mktot+Mki(k,j)
          !      enddo
          !      mp=Mktot/Nki(k)
          !      density=aerodens(Mki(k,srtso4),0.e+0_fp, &
          !                       Mki(k,srtnh4),0.e+0_fp,Mki(k,srth2o))  ! assume bisulfate
          !      ! diameter = ((mass/density)*(6/pi))**(1/3)
          !      d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
          !      ! surface area per particle = pi*diameter**2
          !      surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)* &
          !                              pi*(d2**2.D0) ! (um2 cm-2)
          !   endif
          !enddo
          rhin=dble(RHTOMAS*100.e+0_fp)
          tempin=dble(TEMPTMS)

          call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area, &
                      fn, dum1, rnuc, dum2)
       endif
       if((act_nuc.eq.1).and.(lev.le.7))then
          call bl_nucl(h2so4,fn,rnuc)
       endif
       call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max
       fn = min(fn,fn_c)
       !if (fn.gt.1.0)then
       !   print*, 'fn',fn
       !   print*, 'Yu Yes!'
       !   print*, 'ionrate',ionrate
       !   print*, 'surf_area',surf_area
       !endif
    endif

    if (pdbg) then
       if( bin_nuc == 1 ) then
          print *, 'BINARY cluster form rate : fn',fn
       else
          print *, 'TERNARY cluster form rate: fn',fn
       endif
    endif

    ! if nucleation occured, see how many particles grow to join the first size
    ! section
    if (fn.gt.0.e+0_fp) then

       if(pdbg) print*,'Nki',Nki
       if(pdbg) print*,'Mki',Mki

       call getCondSink_kerm(Nki,Mki,CS,Dpmean,Dp1,dens1,BOXVOL,TEMPTMS,PRES)

       if(pdbg) print*,'CS',CS,'Dpmean',Dpmean,'Dp1',Dp1,'dens1',dens1

       d1 = rnuc*2.e+0_fp*1e-9_fp
       drymass = 0.e+0_fp
       do j=1,icomp-idiag
          drymass = drymass + Mki(1,j)
       enddo
       wetmass = 0.e+0_fp
       do j=1,icomp
          wetmass = wetmass + Mki(1,j)
       enddo

       ! to prevent division by zero (win, 10/1/08)
       if(drymass == 0.e+0_fp) then
          WR = 1.e+0_fp
       else
          WR = wetmass/drymass
       endif

       if(pdbg) print*,'rnuc',rnuc,'WR',WR
       if(pdbg) print*,'d1',d1,'Gci(srtso4)',Gci(srtso4),&
                       'TEMP',temptms,'boxvol',boxvol

       if( IT_IS_NAN( Gci(srtso4) )) then
          print*,'rnuc',rnuc,'WR',WR
          print*,'d1',d1,'Gci(srtso4)',Gci(srtso4)
          call ERROR_STOP('Found NaN in Gci','nucleation')
       endif
       ! print*,'[nucleation] Gci',Gci
       call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS, &
                          boxvol,dens1,gtime)
       if (pdbg) print*,'gtime',gtime

       GR = (Dp1-d1)*1e+9_fp/gtime*3600.e+0_fp ! growth rate, nm hr-1

       gamma = 0.23e+0_fp*(d1*1.0e+9_fp)**(0.2e+0_fp)* &
               (Dp1*1.0d9/3.e+0_fp)**0.075e+0_fp* &
               (Dpmean*1.0e+9_fp/150.e+0_fp)** &
               0.048e+0_fp*(dens1*1.0e-3_fp)** &
               (-0.33e+0_fp)*(TEMPTMS/293.e+0_fp) ! equation 5 in kerminen
       eta = gamma*CS/GR

       if (Dp1.gt.d1)then
          fn1 = fn*exp(eta/(Dp1*1.0e+9_fp)-eta/(d1*1.0e+9_fp))
       else
          fn1 = fn
       endif

       if (pdbg) print*,'eta',eta,'Dp1',Dp1,'d1',d1,'fn1',fn1

       mnuc = sqrt(xk(1)*xk(2))

       nadd = fn1

       nuc_bin = 1

       mold = Mki(nuc_bin,srtso4)
       Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+nadd*mnuc*boxvol*dt
       Nkf(nuc_bin) = Nki(nuc_bin)+nadd*boxvol*dt

       Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
       Gcf(srtnh4) = Gci(srtnh4)

       if (pdbg) then
          print*, 'nadd',nadd
          print *,'Mass add to bin',nuc_bin,'=',nadd*mnuc*boxvol*dt
          print *,'Number added',nadd*boxvol*dt
          print *,'Gcf(srtso4)',Gcf(srtso4)
          print *,'Gcf(srtnh4)',Gcf(srtnh4)
       endif

       do k=1,ibins
          if (k .ne. nuc_bin)then
             Nkf(k) = Nki(k)
             do i=1,icomp
                Mkf(k,i) = Mki(k,i)
             enddo
          else
             do i=1,icomp
                if (i.ne.srtso4) then
                   Mkf(k,i) = Mki(k,i)
                endif
             enddo
          endif
       enddo

       do k=1,ibins
          if (Nkf(k).lt.1.e+0_fp) then
             Nkf(k) = 0.e+0_fp
             do j=1,icomp
                Mkf(k,j) = 0.e+0_fp
             enddo
          endif
       enddo
       !print *, 'mnfix in tomas_mod:2679'

       call mnfix(Nkf,Mkf, ERRORSWITCH)
       pdbg = errorswitch ! carry the error signal from mnfix to outside
       if (errorswitch) print*,'NUCLEATION: Error after mnfix'

       ! there is a chance that Gcf will go less than zero because we are
       ! artificially growing particles into the first size bin. 
       ! don't let it go less than zero.

    else

       do k=1,ibins
          Nkf(k) = Nki(k)
          do i=1,icomp
             Mkf(k,i) = Mki(k,i)
          enddo
       enddo

    endif

    pdbg = errorswitch        ! carry the error signal from mnfix to outside

    RETURN

  END SUBROUTINE NUCLEATION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ezcond
!
! !DESCRIPTION: This subroutine takes a given amount of mass and condenses it
!     across the bins accordingly.
!     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!     Put in GEOS-Chem by Win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZCOND (Nki,Mki,mcondi,spec,Nkf,Mkf,surf_area, &
                     BOXVOL, TEMPTMS, PRES, errswitch)
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nki(ibins) - number of particles per size bin in grid cell
    !Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
    !mcond - mass of species to condense [kg/grid cell]
    !spec - the number of the species to condense
    double precision Nki(ibins), Mki(ibins, icomp)
    double precision mcondi
    REAL*4, INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
    LOGICAL ERRSWITCH   ! signal error to outside

!
! !OUTPUT PARAMETERS:
!
    !Nkf, Mkf - same as above, but final values
    double precision Nkf(ibins), Mkf(ibins, icomp)
    REAL(fp)           surf_area
    integer spec
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision mcond
    double precision pi, R    ! pi and gas constant (J/mol K)
    double precision CS       ! condensation sink [s^-1]
    double precision sinkfrac(ibins+1) ! fraction of CS in size bin
    double precision Nk1(ibins), Mk1(ibins, icomp)
    double precision Nk2(ibins), Mk2(ibins, icomp)
    double precision madd     ! mass to add to each bin [kg]
    double precision maddp(ibins)    ! mass to add per particle [kg]
    double precision mconds ! mass to add per step [kg]
    integer          nsteps            ! number of condensation steps necessary
    integer          my_floor, my_ceil       ! the floor and ceiling (temporary)
    double precision eps     ! small number
    double precision tdt      !the value 2/3
    double precision mpo,mpw  !dry and "wet" mass of particle
    double precision WR       !wet ratio
    double precision tau(ibins) !driving force for condensation
    double precision totsinkfrac ! total sink fraction not including nuc bin
    double precision CSeps    ! lower limit for condensation sink
    double precision tot_m,tot_s    !total mass, total sulfate mass
    double precision ratio    ! used in mass correction
    double precision fracch(ibins,icomp)
    double precision totch

    double precision tot_i,tot_f,tot_fa ! used for conservation of mass check
    LOGICAL          PDBG,  ERRORSWITCH
    real(fp)         zeros(ibins)
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(eps=1.e-40_fp)
    parameter(CSeps=1.e-20_fp)

    !=================================================================
    ! EZCOND begins here
    !=================================================================

    pdbg = errswitch ! take the signal for print debug from outside
    errswitch = .false. !signal to terminate with error. Initialize with .false.

    tdt=2.e+0_fp/3.e+0_fp

    mcond=mcondi

    ! initialize variables
    do k=1,ibins
       Nk1(k)=Nki(k)
       do j=1,icomp
          Mk1(k,j)=Mki(k,j)
       enddo
    enddo

    !print *, 'mnfix in tomas_mod:2804'
    call mnfix(Nk1,Mk1, errorswitch)
    if(errorswitch) then
       print *, 'EZCOND: MNFIX[1] found error --> TERMINATE'
       errswitch=.true.
       return
    endif

    ! get the sink fractions
    ! set Nnuc to zero for this calc
    call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area, &
                     BOXVOL,TEMPTMS, PRES)

    ! make sure that condensation sink isn't too small
    if (CS.lt.CSeps) then     ! just make particles in first bin
       Mkf(1,spec) = Mk1(1,spec) + mcond
       Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
       do j=1,icomp
          if (icomp.ne.spec) then
             Mkf(1,j) = Mk1(1,j)
          endif
       enddo
       do k=2,ibins
          Nkf(k) = Nk1(k)
          do j=1,icomp
             Mkf(k,j) = Mk1(k,j)
          enddo
       enddo
       return
    endif

    if (pdbg) then
       print*,'CS',CS
       print*,'sinkfrac',sinkfrac
       print*,'mcond',mcond
    endif

    ! determine how much mass to add to each size bin
    ! also determine how many condensation steps we need
    totsinkfrac = 0.e+0_fp
    do k=1,ibins
       totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
    enddo
    nsteps = 1
    do k=1,ibins
       if (sinkfrac(k).lt.1.0e-20_fp)then
          madd = 0.e+0_fp
       else
          madd = mcond*sinkfrac(k)/totsinkfrac
       endif
       mpo=0.0
       do j=1,icomp-idiag
          mpo=mpo + Mk1(k,j)
       enddo
       if(mpo == 0.0 ) then  ! prevent division by zero (win, 10/16/08)
          my_floor = 0
       else
          my_floor = int(madd*0.00001/mpo)
       endif
       my_ceil = my_floor + 1
       nsteps = max(nsteps,my_ceil) ! don't let the mass increase by more than 10%
    enddo

    if(pdbg) print*,'nsteps',nsteps

    ! mass to condense each step
    mconds = mcond/nsteps

    ! do steps of condensation
    do i=1,nsteps
       if (i.ne.1) then
          ! set Nnuc to zero for this calculation
          call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area, &
                           BOXVOL,TEMPTMS, PRES)
          totsinkfrac = 0.e+0_fp
          do k=1,ibins
             totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
          enddo
       endif

       tot_m=0.e+0_fp
       tot_s=0.e+0_fp
       do k=1,ibins
          do j=1,icomp-idiag
             tot_m = tot_m + Mk1(k,j)
             if (j.eq.srtso4) then
                tot_s = tot_s + Mk1(k,j)
             endif
          enddo
       enddo

       if (pdbg) print *,'tot_s ',tot_s,' tot_m ',tot_m

       ! change criteria to bigger amount (win, 9/30/08)
       if (mcond.gt.tot_m*5.0e-2_fp) then
          if (pdbg) print *,'Entering TMCOND '

          do k=1,ibins
             mpo=0.0
             mpw=0.0
             !WIN'S CODE MODIFICATION 6/19/06
             !THIS MUST CHANGED WITH THE NEW dmdt_int
             do j=1,icomp-idiag
                mpo = mpo+Mk1(k,j) !accumulate dry mass
             enddo
             do j=1,icomp
                mpw = mpw+Mk1(k,j) ! have wet mass include amso4
             enddo
             if( mpo > 0.0 ) then    ! prevent division by zero (win, 10/16/08)
                WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
             else
                WR = 1.0
             endif
             if (Nk1(k) .gt. 0.e+0_fp) then
                !Change maddp(k) from mass/no. to be just mass (win,10/3/08)
                ! this is because in tmcond here, the moxd argument takes
                ! mass to add for each bin array, not mass/no. array.
                maddp(k) = mconds*sinkfrac(k)/totsinkfrac
                !Prior to 10/3/08 (win)
                !maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
                mpw=mpw/Nk1(k)

                if(pdbg) print*,'mpw',mpw,'maddp',maddp(k),'WR',WR
                !Change the maddp(k) to accordingly -- adding the /Nk1(k) (win, 10/3/08)
                tau(k)=1.5e+0_fp*((mpw+maddp(k)/Nk1(k)*WR)**tdt-mpw**tdt)
                ! Prior to 10/3/08 (win)
                !tau(k)=1.5e+0_fp*((mpw+maddp(k)*WR)**tdt-mpw**tdt) !added WR to moxid term (win, 5/15/06)
                !     tau(k)=0.e+0_fp
                !     maddp(k)=0.e+0_fp
             else
                !nothing in this bin - set tau to zero
                tau(k)=0.e+0_fp
                maddp(k) = 0.e+0_fp
             endif
          enddo
          !print*,'tau',tau
          !print *, 'mnfix in tomas_mod:2942'
          call mnfix(Nk1,Mk1, errorswitch)
          if (errorswitch) then
             print *, 'EZCOND: MNFIX[2] found error --> TERMINATE'
             errswitch=.true.
             return
          endif
          ! do condensation
          errorswitch = pdbg
          !prior to 9/30/08 from Jeff's version
          call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,maddp)

          ! For SO4 condensation, the last argument should be zeroes (win, 9/30/08)
          !zeros(:) = 0.e+0_fp
          !call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,zeros)

          if( errorswitch) then
             errswitch=.true.
             print *,'EZCOND: error after TMCOND --> TERMINATE'
             return
          endif
          errorswitch = pdbg

          !call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
          !jrp totch=0.0
          !jrp do k=1,ibins
          !jrp    do j=1,icomp
          !jrp       fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
          !jrp       totch = totch + (Mk2(k,j)-Mk1(k,j))
          !jrp    enddo
          !jrp enddo
          !print*,'fracch',fracch,'totch',totch

       elseif (mcond.gt.tot_s*1.0e-12_fp) then
          if (pdbg) print *,'Small mcond: distrib w/ sinkfrac '
          if (pdbg) print *, 'maddp(bin) to add to SO4'
          do k=1,ibins
             if (Nk1(k) .gt. 0.e+0_fp) then
                maddp(k) = mconds*sinkfrac(k)/totsinkfrac
             else
                maddp(k) = 0.e+0_fp
             endif
             if(pdbg) print *, maddp(k)
             Mk2(k,srtso4)=Mk1(k,srtso4)+maddp(k)
             do j=1,icomp
                if (j.ne.srtso4) then
                   Mk2(k,j)=Mk1(k,j)
                endif
             enddo
             Nk2(k)=Nk1(k)
          enddo
          if(pdbg) errorswitch = .true.

          !print *, 'mnfix in tomas_mod:2999'
          call mnfix(Nk2,Mk2, errorswitch)
          if(errorswitch) then
             print *, 'EZCOND: MNFIX[3] found error --> TERMINATE'
             errswitch=.true.
             return
          endif
       else ! do nothing
          if (pdbg) print *,'Very small mcond: do nothing!'
          mcond = 0.e+0_fp
          do k=1,ibins
             Nk2(k)=Nk1(k)
             do j=1,icomp
                Mk2(k,j)=Mk1(k,j)
             enddo
          enddo
       endif
       if (i.ne.nsteps)then
          do k=1,ibins
             Nk1(k)=Nk2(k)
             do j=1,icomp
                Mk1(k,j)=Mk2(k,j)
             enddo
          enddo
       endif

    enddo

    do k=1,ibins
       Nkf(k)=Nk2(k)
       do j=1,icomp
          Mkf(k,j)=Mk2(k,j)
       enddo
    enddo

    ! check for conservation of mass
    tot_i = 0.e+0_fp
    tot_fa = mcond
    tot_f = 0.e+0_fp
    do k=1,ibins
       tot_i=tot_i+Mki(k,srtso4)
       tot_f=tot_f+Mkf(k,srtso4)
       tot_fa=tot_fa+Mki(k,srtso4)
    enddo

    if(pdbg) then
       print *,'Check conserv of mass after mcond is distrib'
       print *,' Initial total so4 ',tot_i
       print *,' Final total so4   ',tot_f
       print *,'Percent error=',abs((mcond-(tot_f-tot_i))/mcond)*1e2
    endif

    if ( mcond > 0.0_fp ) then
       if ( abs((mcond-(tot_f-tot_i))/mcond).gt.0.e+0_fp) then
          IF(mcond > 1.e-8_fp .and. tot_i > 5.e-2_fp)  THEN
             !Add a check to check error if mcond is significant (win, 10/2/08)

             IF (abs((mcond-(tot_f-tot_i))/mcond).lt.1.e+0_fp .OR. &
                  spinup(31.0) ) THEN
                !Prior to 10/2/08 (win)   .. original was Jeff's fix
                !! do correction of mass
                !ratio = (tot_f-tot_i)/mcond
                !if(pdbg) print *,'Mk at mass correction '
                !if(pdbg) print *,'  ratio',ratio
                !do k=1,ibins
                !   Mkf(k,srtso4)=Mki(k,srtso4)+
                !   &              (Mkf(k,srtso4)-Mki(k,srtso4))/ratio
                !   if(pdbg) print *,Mkf(k,srtso4)
                !enddo

                ! Do mass correction (win, 10/2/08)
                ratio = (tot_i+mcond)/tot_f
                if(pdbg) print *,'Mk at mass correction apply ratio= ',ratio
                do k=1,ibins
                   Mkf(k,srtso4)=Mkf(k,srtso4) * ratio
                   if(pdbg) print *,Mkf(k,srtso4)
                enddo

                if(pdbg) errorswitch=.true.
                !print *, 'mnfix in tomas_mod:3079'
                call mnfix(Nkf,Mkf, errorswitch)
                if(errorswitch) then
                   print *, 'EZCOND: MNFIX[4] found error --> TERMINATE'
                   errswitch=.true.
                   return
                endif
             else
                print*,'ERROR in ezcond'
                print*,'Condensation error',(mcond-(tot_f-tot_i))/mcond
                print*,'mcond',mcond,'change',tot_f-tot_i
                print*,'tot_i',tot_i,'tot_fa',tot_fa,'tot_f',tot_f
                print*,'Nki',Nki
                print*,'Nkf',Nkf
                print*,'Mki',Mki
                print*,'Mkf',Mkf
                !Prior to 10/2/08 (win)
                !STOP
                ! Send error signal to outside and terminate with more info
                ! (win, 10/2/08)
                !!as of 10/27/08, try comment out this signal to stop the run
                ! (win, 10/27/08)
                !! the problem is that maybe or mostly the mass conservation is
                ! ruined becuase of the fudging inside mnfix.
                !ERRSWITCH=.TRUE.
                !RETURN
             ENDIF
          ENDIF
       endif
    endif

    !jrp if (abs(tot_f-tot_fa)/tot_i.gt.1.0D-8)then
    !jrp    print*,'No S conservation in ezcond'
    !jrp    print*,'initial',tot_fa
    !jrp    print*,'final',tot_f
    !jrp    print*,'mcond',mcond,'change',tot_f-tot_i
    !jrp    print*,'ERROR',(mcond-(tot_f-tot_i))/mcond
    !jrp endif

    ! check for conservation of mass
    tot_i = 0.e+0_fp
    tot_f = 0.e+0_fp
    do k=1,ibins
       tot_i=tot_i+Mki(k,srtnh4)
       tot_f=tot_f+Mkf(k,srtnh4)
    enddo
    if (.not. spinup(14.0)) then
       if (abs(tot_f-tot_i)/tot_i.gt.1.0e-8_fp)then
          if ( tot_i > 1.0e-20_fp ) then
             print*,'No N conservation in ezcond'
             print*,'initial',tot_i
             print*,'final  ',tot_f
          endif
       endif
    endif

    return

  end SUBROUTINE EZCOND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aqoxid
!
! !DESCRIPTION: Subroutine AQOXID takes an amount of SO4 produced via in-cloud
!  oxidation and condenses it onto an existing aerosol size distribution.  It
!  assumes that only particles larger than the critical activation diameter
!  activate and that all of these have grown to roughly the same size.
!  Therefore, the mass of SO4 produced by oxidation is partitioned to the
!  various size bins according to the number of particles in that size bin.
!  Values of tau are calculated for each size bin accordingly and the TMCOND
!  subroutine is called to update Nk and Mk. (win, 7/23/07)
!  Originally written by Peter Adams for TOMAS in GISS GCM-II', June 2000
!  Modified by Win Trivitayanurak (win@cmu.edu), Oct 3, 2005
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AQOXID( MOXID, KMIN, I, J, L, Input_Opt, &
                     State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD             ! ND60
#endif
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    REAL(fp)                      :: MOXID
    INTEGER                       :: KMIN, I, J, L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER,       PARAMETER :: K_MIN = 4
    REAL(fp)                 :: Nact, Mact, MPO, AQTAU(IBINS)
    REAL(fp)                 :: Nko(IBINS), Mko(IBINS, ICOMP)
    REAL(fp)                 :: Nkf(IBINS), Mkf(IBINS, ICOMP)
    REAL(fp),      PARAMETER :: TDT = 2.e+0_fp / 3.e+0_fp
    REAL(fp)                 :: M_OXID(IBINS)
    INTEGER                  :: K, MPNUM, JC, TRACNUM
    INTEGER                  :: NKID
    LOGICAL                  :: PDBG

    REAL(fp)                 :: Nk(IBINS), Nkd(IBINS)
    REAL(fp)                 :: Mk(IBINS, ICOMP)
    REAL(fp)                 :: Mkd(IBINS,ICOMP)
    REAL(fp)                 :: Gc(ICOMP - 1)
    REAL(fp)                 :: Gcd(ICOMP - 1)
    REAL*4                   :: BOXVOL
    REAL*4                   :: thresh
    CHARACTER(LEN=255)       :: MSG, LOC ! (ewl)
    LOGICAL                  :: UNITCHANGE_KGM2

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !=================================================================
    ! AQOXID begins here
    !=================================================================

    ! Assume success
    RC                =  GC_SUCCESS

    ! Check that species units are in [kg] (ewl, 8/13/15)
    ! Convert species concentration units to [kg] if not necessary.
    ! Units are [kg/m2] if AQOXID is called from wet deposition
    ! and are [kg] if called from sulfate_mod since chemistry is
    ! still in [kg]. Since AQOXID is called within an (I,J,L) loop,
    ! only convert units for a single grid box. Otherwise, run will
    ! take too long (ewl, 9/30/15)
    UNITCHANGE_KGM2 = .FALSE.
    IF ( TRIM( State_Chm%Spc_Units ) .eq. 'kg/m2' ) THEN
       UNITCHANGE_KGM2 = .TRUE.
       CALL ConvertBox_Kgm2_to_Kg( I, J, L,  State_Chm, State_Grid, RC )
    ELSE IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine AQOXID in tomas_mod.F90'
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

    ! Point to the chemical species array [kg]
    Spc => State_Chm%Species

    PDBG = .FALSE.            !For print debugging
    !debug IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.

    BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3
    ! Update aerosol water from the current RH
    DO K = 1, IBINS
       CALL EZWATEREQM2( I, J, L, K, State_Met, State_Chm, RC )
    ENDDO


#if   defined ( TOMAS12 ) || defined ( TOMAS15 )
    thresh = 4.0
#else
    thresh = 1.0
#endif

    ! Swap GEOSCHEM variables into aerosol algorithm variables
    DO K = 1, IBINS
       NKID = id_NK1 - 1 + K
       NK(K) = Spc(I,J,L,NKID)
       DO JC = 1, ICOMP-IDIAG
          MK(K,JC) = Spc(I,J,L,NKID+JC*IBINS)
       ENDDO
       MK(K,SRTH2O) = Spc(I,J,L,id_AW1-1+K)
    ENDDO
    !sfarina - initialize Gc to ensure storenm doesn't go NaN on us.
    DO JC=1, ICOMP-1
       Gc(JC) = 0.e+0_fp
    ENDDO

    ! Take the bulk NH4 and allocate to size-resolved NH4
    IF ( SRTNH4 > 0 ) THEN
       CALL NH4BULKTOBIN( MK(:,SRTSO4), Spc(I,J,L,id_NH4), MK(:,SRTNH4) )
    ENDIF

    ! Fix any inconsistencies in M/N distribution
    CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
    !debug IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.

    ! print *, 'mnfix in tomas_mod:3225'
    CALL MNFIX( NK, MK, PDBG )
    IF( PDBG ) THEN
       PRINT *,'AQOXID: MNFIX found error at',I,J,L
       CALL ERROR_STOP('Found bad error in MNFIX', &
                       'Beginning AQOXID after MNFIX' )
    ENDIF
    MPNUM = 5
#ifdef BPCH_DIAG
    IF ( ND60 > 0 ) &
         CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
#endif
    CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

    !debug IF ( I == 46 .AND. J == 59 .AND. L == 9) &
    !      call debugprint(Nk,Mk,I,J,L,'AQOXID after MNFIX_1')

    ! Calculate which particles activate
10  CONTINUE ! Continue here if KMIN has to be lowered
    Nact = 0.e+0_fp
    Mact = 0.e+0_fp
    DO K = KMIN, IBINS
       Nact = Nact + Nk(k)
       DO JC = 1, ICOMP-IDIAG  !accumulate dry mass exclude NH4
          Mact = Mact + Mk(K,JC)
       ENDDO
    ENDDO

    ! No particles to condense on, then just exit AQOXID
    IF ( Nact == 0e+0_fp ) GOTO 20

    ! If condensing mass is too large for the alloted portion of NK
    ! then lower KMIN
    IF ( ( Mact + MOXID )/ Nact > XK(IBINS+1) / thresh ) THEN
       IF ( KMIN > K_MIN ) THEN
          KMIN = KMIN - 1
          GOTO 10
       ELSE
          ! If there is really not enough number to condense onto when lower
          ! KMIN to the threshold K_MIN (set to 4), then
          !  IF current time is within first 2 weeks from initialization
          !    (spin-up), then skip and exit
          !  IF current time is after the first 2 weeks, then terminate
          !    with an error message.
          IF ( .not. SPINUP(14.0) ) THEN
             !WRITE(*,*) 'Location: ',I,J,L
             !WRITE(*,*) 'Kmin/Nact: ',KMIN,NACT
             !WRITE(*,*) 'MOXID/Mact: ',MOXID,Mact
             !DO K = 1, IBINS
             !   WRITE(*,*) 'K, N, MSO4, MH2O: ',K,Nk(k), &
             !        MK(K,SRTSO4),MK(K,SRTH2O)
             !ENDDO
             IF ( MOXID > 5e+0_fp .and. &
                  ( .not. State_Met%InChemGrid(I,J,L) ) ) THEN
                CALL ERROR_STOP( 'Too few number for condensing mass', &
                                 'AQOXID:1'                           )
             ELSE
                WRITE(*,*) 'AQOXID WARNING: SO4 mass is being discarded'
                GOTO 20
             ENDIF
          ELSE
             IF ( PDBG ) print *,'AQOXID: Discard mass (spin-up)'
             GOTO 20
          ENDIF
       ENDIF
    ENDIF

    ! Calculate Tau (driving force) for each size bin
    MOXID = MOXID/ Nact       !Moxid becomes kg SO4 per activated particle
                              !NOTE: NOT using kg of H2SO4
    DO K = 1, IBINS
       IF ( K < KMIN ) THEN
          !too small to activate - no sulfate for this particle
          AQTAU(K) = 0.e+0_fp
          M_OXID(K) = 0.e+0_fp
       ELSE
          !activated particle - calculate appropriate tau
          MPO=0.e+0_fp
          DO JC = 1, ICOMP-IDIAG
             MPO = MPO + Mk(K,JC) !accumulate dry mass
          ENDDO
          M_OXID(K) = MOXID * NK(K)

          IF (Nk(K) > 0.e+0_fp) THEN
             ! Calculate Tau
             MPO = MPO / Nk(K)
             AQTAU(K) = 1.5e+0_fp * ( ( ( MPO + MOXID) ** TDT ) - &
                                      (   MPO          ** TDT )    )

             ! Error checking for negative Tau
             IF ( AQTAU(K) < 0.e+0_fp ) THEN
                IF ( ABS(AQTAU(K)) < 1.e+0_fp ) THEN
                   AQTAU(K)=1.d-50  !0.e+0_fp  !try change to tiny number instead of 0e+0_fp (win, 5/28/06)
                ELSE
                   PRINT *,' ######### aqoxid:  NEGATIVE TAU'
                   PRINT *,'Error at',i,j,l,'bin',k
                   PRINT *,'aqtau(k)',aqtau(k)
                   CALL ERROR_STOP( 'Negative Tau','AQOXID:2' )
                ENDIF
             ENDIF

          ELSE
             ! Nothing in this bin - set tau to zero
             AQTAU(K) = 0.e+0_fp
          ENDIF                  ! Nk>0
       ENDIF                     ! K<kmin
    ENDDO                     ! Loop ibins

    ! Call condensation algorithm

    ! Swap into Nko, Mko
    Mko(:,:) = 0e+0_fp
    DO K = 1, IBINS
       Nko(K) = Nk(K)
       DO JC = 1, ICOMP-IDIAG ! Now do aqoxid "dry" (win, 7/23/07)
          Mko(K,JC) = Mk(K,JC)
       ENDDO

    ENDDO
    !debug IF ( I == 46 .AND. J == 59 .AND. L == 9) PDBG = .TRUE.

    CALL TMCOND( AQTAU, XK, Mko, Nko, Mkf, Nkf, SRTSO4, PDBG, M_OXID )
    IF(.not.SPINUP(60.) .and.  PDBG ) THEN
       write(116,*) 'Error at',i,j,l
    ELSE
       PDBG = .false.
    ENDIF

    ! Swap out of Nkf, Mkf
    DO K = 1, IBINS
       Nk(k)=Nkf(k)
       DO JC = 1, ICOMP-IDIAG
          Mk(K,JC) = Mkf(K,JC)
       ENDDO
    ENDDO

20  CONTINUE ! Continue here if the process is skipped

    ! Save changes to diagnostic
    MPNUM = 4
#ifdef BPCH_DIAG
    IF ( ND60 > 0 ) &
         CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
#endif

    ! Fix any inconsistencies in M/N distribution
    CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
    !print *, 'mnfix in tomas_mod:3371'
    CALL MNFIX( NK, MK, PDBG )
    IF( PDBG ) THEN
       PRINT *,'AQOXID: MNFIX found error at',I,J,L
       CALL ERROR_STOP('Found bad error in MNFIX', &
                       'End of AQOXID after MNFIX' )
    ENDIF
    MPNUM = 5
#ifdef BPCH_DIAG
    IF ( ND60 > 0 ) &
         CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
#endif

    ! Swap Nk and Mk arrays back to Spc
    DO K = 1, IBINS
       TRACNUM = id_NK1 - 1 + K
       Spc(I,J,L,TRACNUM) = Nk(K)
       DO JC = 1, ICOMP-IDIAG
          TRACNUM = id_NK1 - 1 + K + IBINS*JC
          Spc(I,J,L,TRACNUM) = Mk(K,JC)
       ENDDO
       Spc(I,J,L,id_AW1-1+K) = Mk(K,SRTH2O)
    ENDDO

    ! Free pointer memory
    Spc => NULL()

    ! Convert State_Chm%Species units back to original units
    ! if conversion occurred at start of AQOXID (ewl, 9/30/15)
    IF ( UNITCHANGE_KGM2 ) THEN
       CALL ConvertBox_Kg_to_Kgm2( I, J, L, State_Chm, State_Grid, RC )
    ENDIF

    ! Check that species units are as expected (ewl, 9/29/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' .AND. &
         TRIM( State_Chm%Spc_Units ) /= 'kg/m2' ) THEN
       CALL ERROR_STOP('Incorrect final species units:' // State_Chm%Spc_Units,&
                       'Routine AQOXID in tomas_mod.F90')
    ENDIF

  END SUBROUTINE AQOXID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soacond
!
! !DESCRIPTION: Subroutine SOACOND takes the SOA calculated via 10% yeild
!  assumption and condense onto existing aerosol size distribution in a similar
!  manner as in aqoxid. The difference is that SOA condensational driving force
!  is a function of the amount of soluble mass existing in each bin, unlike
!  aqoxid where driving force depends on activated number (proportional to
!  surface) of each bin. (win, 9/25/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SOACOND( MSOA, I, J, L, BOXVOL, TEMPTMS, PRES, &
                      State_Chm, State_Grid, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
#endif
    USE ErrCode_Mod
    USE ERROR_MOD
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    REAL(fp)                      :: MSOA
    INTEGER,        INTENT(IN)    :: I, J, L
    REAL*4,         INTENT(IN)    :: BOXVOL, TEMPTMS, PRES
    TYPE(GrdState), INTENT(IN)    :: State_Grid
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)                 :: surf_area
    REAL(fp)                 :: MPO, OCTAU(IBINS), ntot, mtot
    REAL(fp)                 :: Nko(IBINS), Mko(IBINS, ICOMP)
    REAL(fp)                 :: Nkf(IBINS), Mkf(IBINS, ICOMP)
    REAL(fp),      PARAMETER :: TDT = 2.e+0_fp / 3.e+0_fp
    REAL(fp)                 :: MEDTOT, MED(IBINS), MABS(IBINS)
    REAL(fp)                 :: MKTOT(IBINS), DENSITY, PI
    REAL(fp)                 :: M_NH4
    REAL(fp)                 :: CS, sinkfrac(IBINS)
    REAL(fp)                 :: partfrac(IBINS),avgfrac(IBINS)
    INTEGER                  :: K, MPNUM, JC, TRACNUM
    LOGICAL                  :: PDBG, negvalue ! negvalue added for xSOA
                                                 ! (JKodros, 6/2/15)
    PARAMETER   (PI = 3.141592654e+0_fp)

    !sfarina
    REAL(fp)                 :: Nk(IBINS), Nkd(IBINS)
    REAL(fp)                 :: Mk(IBINS, ICOMP)
    REAL(fp)                 :: Mkd(IBINS,ICOMP)
    REAL(fp)                 :: Gc(ICOMP - 1)
    REAL(fp)                 :: Gcd(ICOMP - 1)
    REAL*4                   :: thresh
    CHARACTER(LEN=255)       :: MSG, LOC ! For species unit check (ewl)
    LOGICAL                  :: ERRORSWITCH = .FALSE.

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    ! For SOACOND warnings
    INTEGER, SAVE       :: SOACOND_WARNING_CT  = -1
    INTEGER, PARAMETER  :: SOACOND_WARNING_MAX = 20

    !=================================================================
    ! SOACOND begins here
    !=================================================================

    ! Assume success
    RC                 = GC_SUCCESS
    SOACOND_WARNING_CT = 0

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine SOACOND in tomas_mod.F90'
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

    ! Point to the chemical species array [kg]
    Spc => State_Chm%Species

    pdbg = .false.

#if   defined ( TOMAS12 ) || defined ( TOMAS15 )
    thresh = 1.0
#else
    thresh = 1.0
#endif

    ! Swap GEOSCHEM variables into TOMAS variables
    DO K = 1, IBINS
       TRACNUM = id_NK1 - 1 + K
       NK(K) = Spc(I,J,L,TRACNUM)
       DO JC = 1, ICOMP-IDIAG  ! do I need aerosol water here?
          TRACNUM = id_NK1 - 1 + K + IBINS*JC
          MK(K,JC) = Spc(I,J,L,TRACNUM)
          IF( IT_IS_NAN( MK(K,JC) ) ) THEN
             PRINT *,'+++++++ Found NaN in SOACOND ++++++++'
             PRINT *,'Location (I,J,L):',I,J,L,'Bin',K,'comp',JC
          ENDIF
       ENDDO
       MK(K,SRTH2O) = Spc(I,J,L,id_AW1-1+K)
    ENDDO

    ! Take the bulk NH4 and allocate to size-resolved NH4
    IF ( SRTNH4 > 0 ) THEN
       CALL NH4BULKTOBIN( MK(:,SRTSO4), Spc(I,J,L,id_NH4), &
                          MK(:,SRTNH4) )
    ENDIF

    CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

    ! Establish an 30-bin array and accculate the total
    ! of the absorbing media.  The choices can be:
    ! organic mass, surface area, organic+inorganic. (win, 3/5/08)

    MEDTOT = 0.e+0_fp
    MED = 0.e+0_fp
    mtot = 0.e+0_fp
    MKTOT(:) = 0.e+0_fp
    ! Accumulate the total absorbing media
    DO K = 1, IBINS

       DO JC = 1, ICOMP
          Mktot(k) = Mktot(k) + Mk(k,jc)
       ENDDO
       mtot = mtot + Mktot(k)

       if (absall.eq.1)then ! partition to all mass
          MED(K) = Mktot(k)
          MEDTOT = MEDTOT + Mktot(k)
       else
          MED(K) = Mk(k,srtocil) ! partition to just hydrophilic organic
          MEDTOT = MEDTOT + Mk(k,srtocil)
       endif

    ENDDO

    ! Fraction to each bin for mass partitioning
    do k = 1,IBINS
       partfrac(k) = MED(K) / MEDTOT ! MSOA (kg SOA) become (kg SOA per
                                     ! total absorbing media)
    enddo

    ! Fraction to each bin for surface condensation
    call getCondSink(Nk,Mk,srtocil,CS,sinkfrac,surf_area, &
                     BOXVOL,TEMPTMS, PRES)

    do k = 1,IBINS
       avgfrac(k)=soaareafrac*sinkfrac(k)+(1.e+0_fp-soaareafrac)*partfrac(k)
    enddo

    !temporary
    ntot = 0.e+0_fp
    do k = 1, ibins
       ntot = ntot + Nk(k)
    enddo

    IF ( ( Mtot + MSOA ) / Ntot > XK(IBINS+1) / thresh ) THEN
       IF ( .not. SPINUP(14.0) ) THEN
          WRITE(*,*) 'Location: ',I,J,L
          WRITE(*,*) 'Mtot_&_Ntot: ',Mtot, Ntot
          IF ( MSOA > 5e+0_fp ) THEN
             CALL ERROR_STOP('Too few no. for SOAcond','SOACOND:1')
          ENDIF
       ELSE
          ! Put a limit on the amount of screen warnings that we get
          ! to keep logfile sizes low (bmy, 9/30/19)
          SOACOND_WARNING_CT = SOACOND_WARNING_CT + 1
          IF ( SOACOND_WARNING_CT < SOACOND_WARNING_MAX ) THEN
             WRITE(*,*) 'SOACOND WARNING: SOA mass is being discarded'
          ENDIF
          GOTO 30
       ENDIF
    ENDIF

    DO K = 1, IBINS
       MPO = 0.e+0_fp
       DO JC = 1, ICOMP-IDIAG
          MPO = MPO + MK(K,JC)  ! Accumulate dry mass
       ENDDO
       MABS(K) = MSOA * avgfrac(K)

       IF ( Nk(K) > 0.e+0_fp ) THEN
          MPO = MPO / Nk(K)
          OCTAU(K) = 1.5e+0_fp * ( ( ( MPO + MABS(K)/Nk(K) ) ** TDT ) - &
                                   (   MPO                   ** TDT )   )

          ! Error checking for negative Tau
          IF ( OCTAU(K) < 0.e+0_fp ) THEN
             IF ( ABS(OCTAU(K)) < 1.e+0_fp ) THEN
                OCTAU(K)=1.e-50_fp  !0.e+0_fp  !try change to tiny number instead of 0e+0_fp (win, 5/28/06)
             ELSE
                PRINT *,' ######### Subroutine SOACOND:  NEGATIVE TAU'
                PRINT *,'Error at',i,j,l,'bin',k
                PRINT *,'octau(k)',octau(k)
                CALL ERROR_STOP( 'Negative Tau','SOACOND:2' )
             ENDIF
          ENDIF

       ELSE
          OCTAU(K) = 0.e+0_fp
       ENDIF
    ENDDO

    ! Call condensation algorithm
    ! Swap into Nko, Mko
    Mko(:,:) = 0.e+0_fp
    DO K = 1, IBINS
       Nko(K) = Nk(K)
       DO JC = 1, ICOMP-IDIAG    ! Now do SOA condensation "dry"
          Mko(K,JC) = Mk(K,JC)  ! dry mass excl. nh4
       ENDDO
    ENDDO
    !debug      if(i==24.and.j==13)       pdbg = .true.
    CALL TMCOND( OCTAU, XK, Mko, Nko, Mkf, Nkf, SRTOCIL, PDBG, MABS )

    ! ----------- JRP ADD MNFIX...This is for xSOA (JKodros 6/2/15) -----------
    if (pdbg) negvalue=.true. !signal received to printdebug (win, 4/8/06)
    call mnfix(Nkf,Mkf,negvalue) !<step5.1> bug fix call argument
    !(win, 4/15/06) !<step4.2> Add call argument to carry tell where mnfix
    !found
    if(negvalue) STOP 'MNFIX terminate' !(win, 9/12/05)
    ! the negative value (win, 9/12/05)
    !---------------------------------------------------------------------------

    IF( PDBG ) THEN
       !print 12, I,J,L
12     FORMAT( 'Error in SOAcond at ', 3I4 )
       if( .not. SPINUP(60.) )write(116,*) 'Error in SOACOND at',i,j,l
    ELSE
       PDBG = .false.
    ENDIF

    ! Swap out of Nkf, Mkf
    DO K = 1, IBINS
       Nk(k)=Nkf(k)
       DO JC = 1, ICOMP-IDIAG
          Mk(K,JC) = Mkf(K,JC)
       ENDDO
    ENDDO

30  CONTINUE

    ! Save changes to diagnostic
    MPNUM = 6
#ifdef BPCH_DIAG
    IF ( ND60 > 0 ) &
       CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
#endif

    ! Fix any inconsistencies in M/N distribution
    !this never happened?

    ! Swap Nk and Mk arrays back to Spc array
    DO K = 1, IBINS
       TRACNUM = id_NK1 - 1 + K
       Spc(I,J,L,TRACNUM) = Nk(K)
       DO JC = 1, ICOMP-IDIAG
          TRACNUM = id_NK1 - 1 + K + IBINS*JC
          Spc(I,J,L,TRACNUM) = Mk(K,JC)
       ENDDO
       Spc(I,J,L,id_AW1-1+K) = Mk(K,SRTH2O)
    ENDDO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE SOACOND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: multicoag
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MULTICOAG( DT, Nk, Mk, BOXVOL, PRES, TEMPTMS, PDBG )
!
! !INPUT PARAMETERS:
!
    REAL*4,    INTENT(IN)     :: DT                ! Time step (s)
    REAL*4,    INTENT(IN)     :: PRES
    REAL*4,    INTENT(IN)     :: TEMPTMS
    REAL*4,    INTENT(IN)     :: BOXVOL
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT)  :: Nk(IBINS)
    REAL(fp),  INTENT(INOUT)  :: Mk(IBINS, ICOMP)
    LOGICAL,   INTENT(INOUT)  :: PDBG              ! For signalling print debug
!
! !REMARKS:
!  Some key variables
!  kij represents the coagulation coefficient (cm3/s) normalized by the
!      volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
!  dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
!     the mass boundaries of the size bins.  xbar is the average mass
!     of a given size bin (it varies with time in this algorithm).  phi
!     and eff are defined in the reference, equations 13a and b.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER     :: K, J, I, JJ
    REAL(fp)    :: dNdt(ibins), dMdt(ibins,icomp-idiag)
    REAL(fp)    :: xbar(ibins), phi(ibins), eff(ibins)
    REAL*4      :: kij(ibins,ibins)
    REAL*4      :: Dpk(ibins)             !diameter (m) of particles in bin k
    REAL*4      :: Dk(ibins)              !Diffusivity (m2/s) of bin k particles
    REAL*4      :: ck(ibins)              !Mean velocity (m/2) of bin k particles
    REAL*4      :: olddiff                !used to iterate to find diffusivity
    REAL*4      :: density                !density (kg/m3) of particles
    REAL*4      :: mu                     !viscosity of air (kg/m s)
    REAL*4      :: mfp                    !mean free path of air molecule (m)
    REAL*4      :: Kn                     !Knudsen number of particle
    REAL(fp)    :: mp                     !particle mass (kg)
    REAL*4      :: beta                   !correction for coagulation coeff.
    !      real(fp), external ::   aerodens  !<tmp> try change to double precision (win, 1/4/06)

    !temporary summation variables
    REAL(fp)    :: k1m(icomp-idiag),k1mx(icomp-idiag)
    REAL(fp)    :: k1mx2(icomp-idiag)
    REAL(fp)    :: k1mtot,k1mxtot
    REAL(fp)    :: sk2mtot, sk2mxtot
    REAL(fp)    :: sk2m(icomp-idiag), sk2mx(icomp-idiag)
    REAL(fp)    :: sk2mx2(icomp-idiag)
    REAL(fp)    :: High_in
    REAL(fp)    :: mtotal, mktot

    REAL*4      :: zeta                      !see reference, eqn 6
    REAL*4      :: tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
    REAL*4      :: dts  !internal time step (<dt for stability)
    REAL*4      :: tsum !time so far
    REAL(fp)    :: Neps !minimum value for Nk
!dbg
    character*12 limit        !description of what limits time step

    REAL(fp)    :: mi, mf   !initial and final masses

#if  defined( TOMAS12 ) || defined( TOMAS15 )
    parameter(zeta=1.28125 , dtlimit=0.25, itlimit=10.)
#else
    parameter(zeta=1.0625, dtlimit=0.25, itlimit=10.)
#endif
    REAL*4      ::pi, kB  !kB is Boltzmann constant (J/K)
    REAL*4      ::R       !gas constant (J/ mol K)
    parameter (pi=3.141592654, kB=1.38e-23, R=8.314, Neps=1.0e-3)

    REAL(fp)      :: M_NH4

    LOGICAL     :: ERRSPOT

    !sfarina
1   format(16E15.3)

    !=================================================================
    ! MULTICOAG begins here!
    !=================================================================
    tsum = 0.0

    ! If any Nk are zero, then set them to a small value to avoid division by zero
    do k=1,ibins
       if (Nk(k) .lt. Neps) then
          Nk(k)=Neps
#if  defined( TOMAS12 ) || defined( TOMAS15 )
          Mk(k,srtso4)=Neps*sqrt( xk(k)*xk(k+1) ) !make the added particles SO4
#else
          Mk(k,srtso4)=Neps*1.4e+0_fp*xk(k) !make the added particles SO4
#endif
       endif
    enddo

    ! Calculate air viscosity and mean free path
    mu=2.5277e-7*temptms**0.75302
    mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temptms)))  !S&P eqn 8.6

    !<temp>
    !write(6,*)'+++ Nk(1:30)    =',Nk(1:30)
    !write(6,*)'+++ Mk(1:30,SO4)=',Mk(1:30,srtso4)
    !write(6,*)'+++ Mk(1:30,H2O)=',Mk(1:30,srth2o)
    if (pdbg) call debugprint(Nk,Mk,0,0,0,'Inside MULTICOAG')
    ! Calculate particle sizes and diffusivities
    do k=1,ibins

       IF ( SRTNH4 > 0 ) THEN
          M_NH4 = Mk(k,SRTNH4)
       ELSE
          M_NH4 = 0.1875e+0_fp*Mk(k,srtso4)  !assume bisulfate
       ENDIF
       !tmp write(6,*)'+++ multicoag:  Mk(',k,'srtso4)=',Mk(k,srtso4)
       density=aerodens(Mk(k,srtso4),0.e+0_fp, M_NH4,        &
               Mk(k,srtnacl), Mk(k,srtecil), Mk(k,srtecob),  &
               Mk(k,srtocil), Mk(k,srtocob), Mk(k,srtdust),  &
               Mk(k,srth2o))     !use Mk for sea salt mass(win, 4/18/06)
       !Update mp calculation to include all species (win, 4/18/06)

       !prior to 9/26/08 (win)
       !Mktot=0.1875e+0_fp*Mk(k,srtso4) !start with NH4 mass

       Mktot = M_NH4         ! start with ammonium (win, 9/26/08)
       Mktot = Mktot + Mk(k,srth2o) ! then include water

       do j=1, icomp-idiag
          Mktot=Mktot+Mk(k,j)
       enddo
       mp=Mktot/Nk(k)
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
       Dk(k)=kB*temptms/(3.0*pi*mu*Dpk(k)) &        !S&P Table 12.1
         *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
       ck(k)=sqrt(8.0*kB*temptms/(pi*mp))           !S&P Table 12.1
    enddo

    ! Calculate coagulation coefficients
    do i=1,ibins
       do j=1,ibins
          Kn=4.0*(Dk(i)+Dk(j)) &
             /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
          beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
          !This is S&P eqn 12.46 with non-continuum correction, beta
          kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
          kij(i,j)=kij(i,j)*1.0e+6_fp/boxvol  !normalize by grid cell volume
       enddo
    enddo

10  continue     !repeat process here if multiple time steps are needed

    if(pdbg) print*,'In the time steps loop +++++++++++++'

    ! Calculate xbar, phi and eff
#if  defined( TOMAS12 ) || defined( TOMAS15 )
    do k=1,ibins

       xbar(k)=0.0
       do j=1,icomp-idiag
          xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
       enddo
       if(k.lt.ibins-1)then !from 1 to 10 bins

          eff(k)=2./9.*Nk(k)/xk(k) *(4.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
          phi(k)=2./9.*Nk(k)/xk(k) *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

          !Constraints in equation 15
          if (xbar(k) .lt. xk(k)) then
             eff(k)=2./3.*Nk(k)/xk(k)
             phi(k)=0.0

          else if (xbar(k) .gt. xk(k+1)) then
             phi(k)=2./3.*Nk(k)/xk(k)
             eff(k)=0.0
          endif
       else                      ! from 11 bins to 12 bins
          eff(k)=2./31./31.*Nk(k)/xk(k) &
                 *(32.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
          phi(k)=2./31./31.*Nk(k)/xk(k) &
                 *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

          !Constraints in equation 15
          if (xbar(k) .lt. xk(k)) then
             eff(k)=2./31.*Nk(k)/xk(k)
             phi(k)=0.0

          else if (xbar(k) .gt. xk(k+1)) then
             phi(k)=2./31.*Nk(k)/xk(k)
             eff(k)=0.0
          endif
       endif

    enddo

#else
    do k=1,ibins

       xbar(k)=0.0
       do j=1,icomp-idiag
          xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
       enddo

       eff(k)=2.*Nk(k)/xk(k)*(2.-xbar(k)/xk(k))    !eqn 13a
       phi(k)=2.*Nk(k)/xk(k)*(xbar(k)/xk(k)-1.)    !eqn 13b

       !Constraints in equation 15
       if (xbar(k) .lt. xk(k)) then
          eff(k)=2.*Nk(k)/xk(k)
          phi(k)=0.0
       else if (xbar(k) .gt. xk(k+1)) then
          phi(k)=2.*Nk(k)/xk(k)
          eff(k)=0.0
       endif
    enddo
#endif

    ! Necessary initializations
    sk2mtot=0.0
    sk2mxtot=0.0
    do j=1,icomp-idiag
       sk2m(j)=0.0
       sk2mx(j)=0.0
       sk2mx2(j)=0.0
    enddo

    ! Calculate rates of change for Nk and Mk
    do k=1,ibins

       !Initialize to zero
       do j=1,icomp-idiag
          k1m(j)=0.0
          k1mx(j)=0.0
          k1mx2(j)=0.0
       enddo
       High_in=0.0
       k1mtot=0.0
       k1mxtot=0.0

       !Calculate sums
#if  defined( TOMAS12 ) || defined( TOMAS15 )
       do j=1,icomp-idiag
          if (k .gt. 1.and.k.lt.ibins) then
             do i=1,k-1
                k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)*zeta
                k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2.*zeta**3.
             enddo
          elseif(k.eq.ibins)then
             k1m(j)= sk2m(j)+kij(k,k-1)*Mk(k-1,j)
             k1mx(j)=sk2mx(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)*4.754
             k1mx2(j)=sk2mx2(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)**2.*107.4365
          endif
          k1mtot=k1mtot+k1m(j)
          k1mxtot=k1mxtot+k1mx(j)
       enddo
#else
       do j=1,icomp-idiag
          if (k .gt. 1) then
             do i=1,k-1
                k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)
                k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2
             enddo
          endif
          k1mtot=k1mtot+k1m(j)
          k1mxtot=k1mxtot+k1mx(j)
       enddo
#endif

       if (k .lt. ibins) then
          do i=k+1,ibins
             High_in=High_in+Nk(i)*kij(k,i)
          enddo
       endif

       !Calculate rates of change
#if  defined( TOMAS12 ) || defined( TOMAS15 )
       if(k.lt.ibins-1)then

          dNdt(k)= -Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.125 &
                   -(phi(k)*k1mtot+(eff(k)-phi(k))/6./xk(k)*k1mxtot) &
                   -kij(k,k)*(phi(k)/3.*xbar(k)*Nk(k)+(eff(k)-phi(k))/18. &
                   /xk(k)*zeta*xbar(k)*xbar(k)*Nk(k))

          if (k .gt. 1) then
             !yhl Nk*low_in changes to -0.5*Kij*Nk**2.
             dNdt(k)=dNdt(k)+0.625*kij(k-1,k-1)*Nk(k-1)**2 &
                     +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1) &
                     *sk2mxtot) &
                     +kij(k-1,k-1)*(phi(k-1)/3.*xbar(k-1)*Nk(k-1)+(eff(k-1) &
                     -phi(k-1))/18./xk(k-1)*zeta*xbar(k-1)*xbar(k-1) &
                     *Nk(k-1))

          endif

          do j=1,icomp-idiag

             dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+ &
                        (eff(k)+2.*phi(k))/6.*k1mx(j) &
                        +(eff(k)-phi(k))/6./xk(k)*k1mx2(j)) & ! term3
                        - kij(k,k)*Nk(k)*Mk(k,j)/3. & ! I assume 1/2Nk and 2/3Mk for half bin
                        - kij(k,k)*(phi(k)*xk(k+1)*Mk(k,j)/3. &
                        +(eff(k)+2.*phi(k))/6.*zeta*xbar(k)*Mk(k,j)/3. &
                        +(eff(k)-phi(k))/6./xk(k)*zeta**3.*xbar(k)**2. &
                        *Mk(k,j)/3.)

             !yhl  Term9(-kij(k,k)*Nk(k)*Mk(k,j)) is cancled out by term6 (k)
             if (k .gt. 1) then
                dMdt(k,j)=dMdt(k,j) &
                          +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1) &
                          +2.*phi(k-1))/6.*sk2mx(j) &
                          +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) & !term1
                          +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3. &
                          +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3. &
                          +(eff(k-1)+2.*phi(k-1))/6.*zeta &
                          *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6. &
                          /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)
             endif
          enddo
       else if (k.eq.ibins-1)then

          dNdt(k)=0.625*kij(k-1,k-1)*Nk(k-1)**2 &
                  +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1) &
                  *sk2mxtot) &
                  +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)/3.*(phi(k-1) &
                  +(eff(k-1) &
                  -phi(k-1))/6./xk(k-1)*zeta*xbar(k-1))

          !yhl updated the following
          dNdt(k)=dNdt(k)-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.02 &
                  -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot) &
                  -kij(k,k)*xbar(k)*Nk(k)*0.484*(phi(k)+(eff(k) &
                  -phi(k))/62./xk(k)*4.754*xbar(k))

          !yhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
          !yhl It results in much closer result as 30 bins. Apr.27.08

          do j=1,icomp-idiag
             dMdt(k,j)= &
                        +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1) &
                        +2.*phi(k-1))/6.*sk2mx(j) &
                        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) & !term1
                        +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3. &
                        +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3. &
                        +(eff(k-1)+2.*phi(k-1))/6.*zeta &
                        *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6. &
                        /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)

             !yhl updated the following
             dMdt(k,j)= dMdt(k,j)+Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k)) &
                        *k1mx(j)+(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) & ! term3
                        -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 & ! I assume 1/2Nk and 2/3Mk for half bin
                        -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62. &
                        +0.484*phi(k))*4.754*xbar(k) &
                        +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
          enddo

       else if (k.eq.ibins)then
          dNdt(k)=-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.103226 &
                  -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot) &
                  -kij(k,k)*0.484*xbar(k)*Nk(k)*(phi(k)+(eff(k)-phi(k)) &
                  /62./xk(k)*4.754*xbar(k)) &
                  +0.52*kij(k-1,k-1)*Nk(k-1)**2 &
                  +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/62./xk(k-1) &
                  *sk2mxtot) &
                  +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)*0.484*(phi(k-1) &
                  +(eff(k-1) &
                  -phi(k-1))/62./xk(k-1)*4.754*xbar(k-1))

          do j=1,icomp-idiag
             dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484 &
                        *phi(k))*k1mx(j) &
                        +(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) & ! term3
                        -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 & ! I assume 1/2Nk and 2/3Mk for half bin
                        -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62. &
                        +0.484*phi(k))*4.754*xbar(k) &
                        +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.) &
                        +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)/62.+0.484 &
                        *phi(k-1)) &
                        *sk2mx(j)+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mx2(j)) & !term1
                        +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)*0.103226 &
                        +kij(k-1,k-1)*Mk(k-1,j)*0.484*(phi(k-1)*xk(k) &
                        +(eff(k-1)/62.+0.484*phi(k-1))*4.754*xbar(k-1) &
                        +(eff(k-1)-phi(k-1))/62./xk(k-1)*107.4365 &
                        *xbar(k-1)**2.)
          enddo
       endif

#else
       dNdt(k)= &
                -kij(k,k)*Nk(k)**2 &
                -phi(k)*k1mtot &
                -zeta*(eff(k)-phi(k))/(2*xk(k))*k1mxtot &
                -Nk(k)*High_in
       if (k .gt. 1) then
          dNdt(k)=dNdt(k)+ &
                  0.5*kij(k-1,k-1)*Nk(k-1)**2 &
                  +phi(k-1)*sk2mtot &
                  +zeta*(eff(k-1)-phi(k-1))/(2*xk(k-1))*sk2mxtot
       endif

       do j=1,icomp-idiag
          dMdt(k,j)= &
                     +Nk(k)*k1m(j) &
                     -kij(k,k)*Nk(k)*Mk(k,j) &
                     -Mk(k,j)*High_in &
                     -phi(k)*xk(k+1)*k1m(j) &
                     -0.5*zeta*eff(k)*k1mx(j) &
                     +zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
          if (k .gt. 1) then
             dMdt(k,j)=dMdt(k,j)+ &
                       kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j) &
                       +phi(k-1)*xk(k)*sk2m(j) &
                       +0.5*zeta*eff(k-1)*sk2mx(j) &
                       -zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
          endif
          !dbg if (j. eq. srtso4) then
          !dbg    if (k. gt. 1) then
          !dbg       write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j), &
          !dbg          Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j), &
          !dbg          0.5*zeta*eff(k)*k1mx(j), &
          !dbg          zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j), &
          !dbg          kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j), &
          !dbg          phi(k-1)*xk(k)*sk2m(j), &
          !dbg          0.5*zeta*eff(k-1)*sk2mx(j), &
          !dbg          zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
          !dbg    else
          !dbg       write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j), &
          !dbg          Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j), &
          !dbg          0.5*zeta*eff(k)*k1mx(j), &
          !dbg          zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
          !dbg    endif
          !dbg endif
       enddo
#endif

       !dbg
       if(pdbg) write(*,*) 'k,dNdt,dMdt: ', k, dNdt(k), dMdt(k,srtso4)

       !Save the summations that are needed for the next size bin
       sk2mtot=k1mtot
       sk2mxtot=k1mxtot
       do j=1,icomp-idiag
          sk2m(j)=k1m(j)
          sk2mx(j)=k1mx(j)
          sk2mx2(j)=k1mx2(j)
       enddo

    enddo  !end of main k loop

    ! Update Nk and Mk according to rates of change and time step

    !If any Mkj are zero, add a small amount to achieve finite
    !time steps
    do k=1,ibins
       do j=1,icomp-idiag
          if (Mk(k,j) .eq. 0.e+0_fp) then
             !add a small amount of mass
             mtotal=0.e+0_fp
             do jj=1,icomp-idiag
                mtotal=mtotal+Mk(k,jj)
             enddo
             Mk(k,j)=1.e-10_fp*mtotal
          endif
       enddo
    enddo

    call mnfix(NK, MK, PDBG)

    !Choose time step
    dts=dt-tsum      !try to take entire remaining time step
    limit='comp'
    do k=1,ibins
       if(pdbg) print*,'At bin ',k
       if (Nk(k) .gt. Neps) then
          !limit rates of change for this bin
          if (dNdt(k) .lt. 0.0) tlimit=dtlimit
          if (dNdt(k) .gt. 0.0) tlimit=itlimit
          if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then
             dts=Nk(k)*tlimit/abs(dNdt(k))
             if(pdbg) print*,'tlimit',tlimit,'Nk(',k,')',Nk(k), &
                             'dNdt',dNdt(k), ' == dts ',dts
             limit='number'
             if(pdbg) write(limit(8:9),'(I2)') k
             if(pdbg) write(*,*) Nk(k), dNdt(k)
          endif
          do j=1,icomp-idiag
             !limit rates of change x(win, 4/22/06)
             if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
             if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
             if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then
                mtotal=0.e+0_fp
                do jj=1,icomp-idiag
                   mtotal=mtotal+Mk(k,jj)
                enddo
                !only use this criteria if this species is significant
                if ((Mk(k,j)/mtotal) .gt. 1.e-5_fp) then
                   dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
                   if(pdbg) print*,'tlimit',tlimit,'Mk(',k,j,')',Mk(k,j), &
                                   'dMdt',dMdt(k,j), ' == dts ',dts
                else
                   if (dMdt(k,j) .lt. 0.0) then
                      !set dmdt to 0 to avoid very small mk going negative
                      dMdt(k,j)=0.0
                      if(pdbg) print*,' dMdt(k,j) < 0 '
                   endif
                endif
                limit='mass'
                if(pdbg) write(limit(6:7),'(I2)') k
                if(pdbg) write(limit(9:9),'(I1)') j
                if(pdbg) write(*,*) Mk(k,j), dMdt(k,j)
             endif
          enddo
       else
          !nothing in this bin - don't let it affect time step
          Nk(k)=Neps
#if  defined( TOMAS12 ) || defined( TOMAS15 )
          Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
#else
          Mk(k,srtso4)=Neps*1.4e+0_fp*xk(k) !make the added particles SO4
#endif
          !make sure mass/number don't go negative
          if (dNdt(k) .lt. 0.0) dNdt(k)=0.0
          if (pdbg) print*,' dNdt(k) < 0 '
          do j=1,icomp-idiag
             if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
          enddo
       endif
    enddo  !loop bin

    if (pdbg .and. dts .lt. 1. ) then
       write(*,*), dts, 'dts < 1. in multicoag'
    endif

    if (dts .eq. 0.) then
       write(*,*) 'time step is 0 in multicoag - inf/nan/tiny error'
       !pause
       do k = 1,ibins
          print *, 'dNdt(k)', dNdt(k)
          print *, 'dMdt(k,j)'
          do j = 1,icomp - idiag
             print *, dMdt(k,j)
          end do
       end do

       call debugprint(nk, mk, 0,0,0,'MULTICOAG before terminate: dts=0')
       PDBG = .true.
       return
       !stop
       !go to 20
    endif

    !Change Nk and Mk
    !dbg
    if(pdbg) write(*,*) 'tsum=',tsum+dts,' ',limit
    do k=1,ibins
       Nk(k)=Nk(k)+dNdt(k)*dts
       do j=1,icomp-idiag
          Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
       enddo
    enddo

    !Update time and repeat process if necessary
    tsum=tsum+dts
    if (tsum .lt. dt) then
       !print*,'tsum',tsum, 'less than 3600. loop again'
       goto 10
    endif

    RETURN

  END SUBROUTINE MULTICOAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: so4cond
!
! !DESCRIPTION: Subroutine SO4COND determines the condensational driving force
!  for mass transfer of sulfate between gas and aerosol phases.  It then calls
!  a mass- and number-conserving algorithm for condensation (/evaporation) of
!  aerosol.
!  .
!  An adaptive time step is used to prevent numerical difficulties.
!  To account for the changing gas phase concentration of sulfuric
!  acid, its decrease during a condensational time step is well-
!  approximated by an exponential decay with a constant, sK (Hz).
!  sK is calculated from the mass and number distribution of the
!  aerosol.  Not only does this approach accurately take into account
!  the changing sulfuric acid concentration, it is also used to
!  predict (and limit) the final sulfuric acid concentration.  This
!  approach is more accurate and faster (allows longer condensational
!  time steps) than assuming a constant supersaturation of sulfate.
!  .
!  Written by Peter Adams, June 2000, based on thermocond.f
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) July 2007

!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SO4COND(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt, &
                     BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, errspot)
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
!
! !INPUT PARAMETERS:
!
    ! Initial values of
    ! Nki(ibins) - number of particles per size bin in grid cell
    ! Mki(ibins, icomp) - mass of a given species per size bin/grid cell
    ! Gci(icomp-1) - amount (kg/grid cell) of all species present in the
    !                gas phase except water
    ! dt - total model time step to be taken (s)
    REAL(fp)           :: Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
    REAL*4             :: dt
    LOGICAL            :: errspot

    REAL*4, INTENT(IN) :: BOXVOL,  BOXMASS, TEMPTMS
    REAL*4, INTENT(IN) :: PRES, RHTOMAS
!
! !OUTPUT PARAMETERS:
!
    ! Nkf, Mkf, Gcf - same as above, but final values
    REAL(fp)      :: Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)    :: dp(ibins, icomp-1)  !Driving force for condensation (Pa)
    REAL(fp)    :: tau(ibins)          !condensation parameter (see cond)
    REAL(fp)    :: atau(ibins, icomp)  !same as tau, but all species
    REAL(fp)    :: atauc(ibins, icomp) !same as atau, but for const dp
    REAL*4      :: time                !amount of time (s) that has been simulated
    REAL*4      :: cdt                 !internal, adaptive time step
    REAL*4      :: mu                  !viscosity of air (kg/m s)
    REAL*4      :: mfp                 !mean free path of air molecule (m)
    REAL*4      :: Kn                  !Knudsen number of particle
    REAL*4      :: Dpk(ibins)          !diameter (m) of particles in bin k
    REAL*4      :: density             !density (kg/m3) of particles
    INTEGER     :: j,k,jj,kk        !counters
    REAL(fp)    :: tj(icomp-1), tk(ibins)  !factors used for calculating tau
    REAL(fp)    :: sK                  !exponential decay const for H2SO4(g)
    REAL(fp)    :: pi, R            !constants
    REAL(fp)    :: zeta13             !from Eqn B30 of Tzivion et al. (1989)
    REAL*4      ::Di                  !diffusivity of gas in air (m2/s)
    REAL*4      ::gmw(icomp-1)        !molecular weight of condensing gas
    REAL*4      ::Sv(icomp-1)         !parameter used for estimating diffusivity
    REAL*4      ::alpha(icomp-1)      !accomodation coefficients
    REAL*4      ::beta                !correction for non-continuum
    REAL(fp)    :: mp         !particle mass (kg)
    REAL(fp)    :: Nko(ibins), Mko(ibins, icomp), Gco(icomp-1) !output of cond routine
    REAL(fp)    :: mi, mf  !initial and final aerosol masses (updates Gc)
    REAL(fp)    :: tr      ! used to calculate time step limits
    REAL(fp)    :: mc, ttr
    REAL(fp)    :: Neps     !value below which Nk is insignificant
    REAL(fp)    :: cthresh  !determines minimum gas conc. for cond.
    !dbg      character*12 limit        !description of what limits time step
    REAL(fp)    :: tdt      !the value 2/3
    REAL(fp)    :: Ntotf, Ntoto, dNerr  !used to track number cons.
    !dbg      integer numcalls          !number of times cond routine is called
    REAL(fp)    :: Mktot        ! total mass (win, 4/14/06)
    REAL(fp)    :: zeros(IBINS)

    LOGICAL     :: negvalue  ! negative check variable
    LOGICAL     :: printdebug ! signal received from aerophys to print values for debug (win, 4/8/06)
    LOGICAL     :: tempvar    ! just a temporary variable (win, 4/12/06)

    !REAL(fp), EXTERNAL :: AERODENS
    !REAL, EXTERNAL ::  GASDIFF

    PARAMETER(PI=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    PARAMETER(Neps=1.0e10_fp, zeta13=0.98483, cthresh=1.e-16_fp)

    !=================================================================
    ! SO4COND begins here!
    !=================================================================

    negvalue = .false.
    printdebug  = .false.
    tempvar  = .false.

    ! Set some constants
    ! Note: Could have declare this using DATA statement but don't want to
    !       keep modifying when changing the multi-component mass species
    DO J = 1, ICOMP-1
       IF( J == SRTSO4 ) THEN
          gmw(J)  = 98.
          Sv(J)   = 42.88
          alpha(J)= 0.65
          !alpha from U. Poschl et al., J. Phys. Chem. A, 102, 10082-10089, 1998
       ELSE IF( J == SRTNACL ) THEN
          gmw(J)  = 0.
          Sv(J)   = 42.88  !use 42.88 for all components following Jeff Pierce's code (win,9/26/08)
          alpha(J)= 0.
       ELSE IF( J == SRTECOB .or. J == SRTECIL .or. &
                J == SRTOCOB .or. J == SRTOCIL ) THEN
          gmw(J)  = 12.          ! check these values with Jeff again (win, 8/22/07)
          Sv(J)   = 42.88
          alpha(J)= 0.
       ELSE IF( J == SRTDUST ) THEN
          gmw(J)  = 0.
          Sv(J)   = 42.88
          alpha(J)= 0.
       ELSE IF( J == SRTNH4 ) THEN
          gmw(J)  = 0.
          Sv(J)   = 42.88
          alpha(J)= 0.
       ELSE
          PRINT *, 'Modify SO4cond for the new species'
          CALL ERROR_STOP('SO4COND','Need values for Gmw, Sv, alpha')
       ENDIF
    ENDDO

    !dbg numcalls=0
    printdebug = errspot !taking the signal to printdebug from aerophys (win, 4/8/06)
    errspot = .true. !<step4.4> Flag for showing error location outside this subroutine (win,9/21/05)

    dNerr=0.0
    tdt=2.e+0_fp/3.e+0_fp

    ! Initialize values of Nkf, Mkf, Gcf, and time
    !--------------------------------------------------
    TIME = 0.0                !subroutine exits when time=dt
    DO J = 1, ICOMP-1
       GCF(J) = GCI(J)
    ENDDO
    DO K = 1, IBINS
       NKF(K) = NKI(K)
       DO J = 1, ICOMP
          MKF(K,J) = MKI(K,J)
       ENDDO
    ENDDO

    !Leave everything the same if nothing to condense
    IF ( GCI(SRTSO4) < CTHRESH * BOXMASS ) GOTO 100

    IF ( PRINTDEBUG ) PRINT*,'COND NOW: H2SO4=',Gci(srtso4)

    ! Repeat from this point if multiple internal time steps are needed
    !------------------------------------------------------------------
10  CONTINUE

    ! Call thermodynamics to get dp forcings for volatile species
    do k=1,ibins
       do j = 1, icomp-1
          dp(k,j)=0.0
       enddo
    enddo

    ! Set dp for nonvolatile species
    do k=1,ibins
       !<step4.5> correct the MW of Gcf(srtso4) to be 98. (win, 10/13/05)
       dp(k,srtso4)=(Gcf(srtso4)/98.)/(boxmass/28.9)*pres
    enddo

    ! Calculate tj and tk factors needed to calculate tau values
    mu=2.5277e-7*temptms**0.75302
    mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temptms)))  !S&P eqn 8.6
    do j=1,icomp-1
       Di=gasdiff(temptms,pres,gmw(j),Sv(j))
       tj(j)=2.*pi*Di*molwt(j)*1.0e-3_fp/(R*temptms)
    enddo
    sK=0.0e+0_fp
    do k=1,ibins
       if (Nkf(k) .gt. Neps) then
          density=aerodens(Mkf(k,srtso4),0.e+0_fp, &
                  0.1875e+0_fp*Mkf(k,srtso4),Mkf(k,srtnacl), &
                  Mkf(k,srtecil), Mkf(k,srtecob), &
                  Mkf(k,srtocil), Mkf(k,srtocob), &
                  Mkf(k,srtdust), &
                  Mkf(k,srth2o))
          !factor of 1.2 assumes ammonium bisulfate
          !(NH4)H has MW of 19 which is = 0.2*96
          !So the Mass of ammonium bisulfate = 1.2*mass sulfate
          !win, 4/14/06             mp=(1.2*Mkf(k,srtso4)+Mkf(k,srth2o))/Nkf(k)
          !Need to include new mass species in mp (win, 4/14/06)
          !Add 0.1875x first for ammonium, and then add 1.0x in the loop (win, 4/14/06)
          Mktot=0.1875e+0_fp*Mkf(k,srtso4)
          do j=1, icomp
             Mktot=Mktot+Mkf(k,j)
          enddo
          mp=Mktot/Nkf(k)

       else
          !nothing in this bin - set to "typical value"
          density=1500.
#if  defined( TOMAS12 ) || defined( TOMAS15 )
          mp=sqrt(xk(k)*xk(k+1))
#else
          mp=1.4*xk(k)
#endif
       endif
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)
       beta=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtso4)) !S&P eqn 11.35
       tk(k)=(6./(pi*density))**(1./3.)*beta
       if (Nkf(k) .gt. 0.0) then
          do kk=1,icomp
             sK=sK+tk(k)*Nkf(k)*(Mkf(k,kk)/Nkf(k))**(1.e+0_fp/3.e+0_fp)  !<step5.1> (win, 4/14/06)
          enddo
       endif
    enddo  !bin loop
    sK=sK*zeta13*tj(srtso4)*R*temptms/(molwt(srtso4)*1.e-3_fp)/(boxvol*1.e-6_fp)

    ! Choose appropriate time step

    !Try to take full time step
    cdt=dt-time
    !dbg limit='complete'

    !Not more than 15 minutes
    if (cdt .gt. 900.) then
       cdt=900.
       !dbg limit='max'
    endif

20  continue   !if time step is shortened, repeat from here

    !Calculate tau values for all species/bins
    do k=1,ibins
       do j=1,icomp
          atauc(k,j)=0.e+0_fp
          atau(k,j)=0.e+0_fp
       enddo
       !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if(printdebug)then
          write(*,*)'+++ k loop at',k !<temp>
          write(*,*)'+++ tj(srtso4)', tj(srtso4) !<temp>
          write(*,*)'+++ dp(k,srtso4)', dp(k,srtso4) !<temp>
          write(*,*)'+++ tk(k)',tk(k) !<temp>
          write(*,*)'+++ cdt',cdt !<temp>
       endif
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
       atauc(k,srtso4)=tj(srtso4)*tk(k)*dp(k,srtso4)*cdt

       if (sK .gt. 0.e+0_fp) then
          atau(k,srtso4)=tj(srtso4)*R*temptms/(molwt(srtso4)*1.e-3_fp) &
                         /(boxvol*1.e-6_fp)*tk(k)*Gcf(srtso4)/sK &
                         *(1.e+0_fp-exp(-1.e+0_fp*sK*cdt))
       else
          !nothing to condense onto
          atau(k,srtso4)=0.e+0_fp
       endif
    enddo

    !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (printdebug) then
       do j=1,icomp
          print *,'atauc(1:ibins,comp) at comp',j
          print *,atauc(1:ibins,j)
          print *,'atau(1:ibins,comp) at comp',j
          print *,atau(1:ibins,j)
       enddo
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !The following sections limit the condensation time step
    !when necessary.  tr is a factor that describes by
    !how much to reduce the time step.
    tr=1.e+0_fp  !make sure tr is double precision (win, 3/20/05)

    !Make sure masses of individual species don't change too much
    do j=1,icomp-1
       do k=1,ibins
          if (Nkf(k) .gt. Neps) then
             mc=0.e+0_fp
             do jj=1,icomp
                mc=mc+Mkf(k,jj)/Nkf(k)
             enddo
             if (mc/xk(k) .gt. 1.0e-3_fp) then
                !species has significant mass in particle - limit change
                if (abs(atau(k,j))/mc**(2.e+0_fp/3.e+0_fp) > 0.1) then
                   ttr=abs(atau(k,j))/mc**(2.e+0_fp/3.e+0_fp)/5.e-2_fp
                   if (ttr .gt. tr) then
                      tr=ttr
                      !dbg limit='amass'
                      !dbg write(limit(7:11),'(I2,X,I2)') k,j
                   endif
                endif
             else
                !species is new to particle - set max time step
                if ((cdt/tr.gt.1.e-1_fp).and.(atau(k,j)>0.e+0_fp)) then
                   tr=cdt/1.e-1_fp !Make sure tr is double precision (win,3/20/05)
                   !dbg limit='nspec'
                   !dbg write(limit(7:11),'(I2,X,I2)') k,j
                endif
             endif
          endif
       enddo
       !Make sure gas phase concentrations don't change too much
       if (exp(-1.e+0_fp*sK*cdt) .lt. 2.5e-1_fp) then
          ttr=-2.e+0_fp*cdt*sK/log(2.5e-1_fp)
          if (ttr .gt. tr) then
             tr=ttr
             !dbg limit='gphas'
             !dbg write(limit(7:8),'(I2)') j
          endif
       endif
    enddo

    !Never shorten timestep by less than half
    if (tr .gt. 1.e+0_fp) tr=max(tr,2.e+0_fp) !make sure tr is double precision (win,3/20/05)

    !Repeat for shorter time step if necessary
    if (tr .gt. 1.e+0_fp) then  !make sure tr is double precision (win,3/20/05)
       cdt=cdt/tr
       goto 20
    endif

    ! Call condensation subroutine to do mass transfer

    do j=1,icomp-1  !Loop over all aerosol components

       !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if(printdebug) print *,'Call condensation at comp',j
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       !Swap tau values for this species into array for cond
       do k=1,ibins
          tau(k)=atau(k,j)
       enddo

       !dbg write(*,*) 'so4cond - time = ', time, ' ',limit
       !dbg if (j .eq. srtso4) then
       !dbg    do k=1,ibins
       !dbg       write(*,'(I3,4E12.4)') &
       !dbg         k,sK,cdt,atauc(k,srtso4),atau(k,srtso4)
       !dbg    enddo
       !dbg endif

       if (printdebug) negvalue=.true. !signal received to printdebug (win, 4/8/06)
       call mnfix(Nkf,Mkf,negvalue) !<step5.1> bug fix call argument (win, 4/15/06) !<step4.2> Add call argument to carry tell where mnfix found
       ! the negative value (win, 9/12/05)
       if ( negvalue ) STOP 'MNFIX terminate' !(win, 9/12/05)

       !Call condensation routine
       Ntotf=0.e+0_fp  !Force double precision (win, 4/20/06)
       do k=1,ibins
          Ntotf=Ntotf+Nkf(k)
       enddo

       !<step5.1> Skip tmcond call if there is absolutely no particle (win, 4/20/06)
       if(Ntotf.gt.0e+0_fp) then

          !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(printdebug) print *,'=== Entering TMCOND ==='
          tempvar = printdebug
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          zeros(:) = 0.e+0_fp
          call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,printdebug,zeros)
          !dbg numcalls=numcalls+1
          errspot = printdebug !receive the error signal from inside tmcond (win,4/12/06)
          printdebug = tempvar !printdebug gets the originally assigned value (win, 4/12/06)
          !8/2/07 if(errspot) goto 100 !Exit so4cond right away when found error from tmcond. (win, 4/13/06)

          !Check for number conservation
          Ntoto=0.0
          do k=1,ibins
             Ntoto=Ntoto+Nko(k)
          enddo
          !dbg write(*,*) 'Time=', time
          !dbg write(*,*) 'Ntoto=', Ntoto
          !dbg write(*,*) 'Ntotf=', Ntotf
          dNerr=dNerr+Ntotf-Ntoto
          if (abs(dNerr/Ntoto) .gt. 1.e-4) then
             write(*,*) 'ERROR in so4cond: Number not conserved'
             write(*,*) 'time=',time
             write(*,*) Ntoto, Ntotf
             write(*,*) (Nkf(k),k=1,ibins)
             errspot = .true. !<step4.4> This flag will trigger printing of location with error (win, 9/21/05)
          endif

       else !(win, 4/20/06)
          if(printdebug) print *,'so4cond: Nk=0 -> skip tmcond'
          do k=1,ibins
             nko(k) = 0e+0_fp
             do jj=1,icomp-1
                Mko(k,jj) = 0e+0_fp
             enddo
          enddo
       endif !(win, 4/20/06)

       if(printdebug) print *,'Initial gas conc:',Gcf(j)  !<temp> (win, 4/11/06)

       !Update gas phase concentration
       mi=0.0
       mf=0.0
       do k=1,ibins
          mi=mi+Mkf(k,j)
          mf=mf+Mko(k,j)
       enddo
       Gcf(j)=Gcf(j)+(mi-mf)*gmw(j)/molwt(j)

       if(printdebug) print *,'Updated gas conc:',Gcf(j) !<temp> (win, 4/11/06)

       !Swap into Nkf, Mkf
       do k=1,ibins
          Nkf(k)=Nko(k)
          do jj=1,icomp-1
             Mkf(k,jj)=Mko(k,jj)
          enddo
       enddo
       
       !Update water concentrations
       call ezwatereqm(Mkf, RHTOMAS)

    enddo

    ! Update time
    time=time+cdt
    !dbg write(*,*) 'so4cond - time = ', time, ' ',limit
    !dbg write(*,*) 'H2SO4(g)= ', Gcf(srtso4)
    if (Gcf(srtso4) .lt. 0.0) then
       if (abs(Gcf(srtso4)) .gt. 1.e-5_fp) then
          !Gcf is substantially less than zero - this is a problem
          write(*,*) 'ERROR in so4cond: H2SO4(g) < 0'
          write(*,*) 'time=',time
          write(*,*) 'Gcf()=',Gcf(srtso4)
          !4/11/06 STOP
          !Let the run STOP outside so4cond so I can know where the run died (win, 4/11/06)
          errspot=.true. !win, 4/11/06
       else
          !Gcf is negligibly less than zero - probably roundoff error
          Gcf(srtso4)=0.0
       endif
    endif

    ! Repeat process if necessary
    if (time .lt. dt) goto 10

    !dbg write(*,*) 'Cond routine called ',numcalls,' times'
    !dbg write(*,*) 'Number cons. error was ', dNerr

100 continue   !skip to here if there is no gas phase to condense

    RETURN
  END SUBROUTINE SO4COND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tmcond
!
! !DESCRIPTION: Subroutine TMCOND do condensation calculation.
!  Original code from Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayaurak (win@cmu.edu)
!  CONDENSATION
!   Based on Tzivion, Feingold, Levin, JAS 1989 and
!   Stevens, Feingold, Cotton, JAS 1996
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,pdbug,moxd)
!
! !INPUT PARAMETERS:
!
    ! TAU(k) ...... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
    ! X(K) ........ Array of bin limits in mass space
    ! AMKD(K,J) ... Input array of mass moments
    ! ANKD(K) ..... Input array of number moments
    ! CSPECIES .... Index of chemical species that is condensing
    REAL(fp)       :: TAU(ibins)
    REAL(fp)       :: X(ibins+1),AMKD(ibins,icomp),ANKD(ibins)
    INTEGER        :: CSPECIES
    LOGICAL        :: pdbug !(win, 4/10/06)
    REAL(fp)       :: moxd(IBINS) ! condensing mass distributed to size bins
                       ! according to the selected absorbing media (win, 3/5/08)
!
! !OUTPUT PARAMETERS:
!
    ! AMK(K,J) .... Output array of mass moments
    ! ANK(K) ...... Output array of number moments
    REAL(fp)       :: AMK(ibins,icomp),ANK(ibins)
!
! !REMARKS:
! The supersaturation is calculated outside of the routine and assumed
! to be constant at its average value over the timestep.
! .
! The method has three basic components:
! (1) first a top hat representation of the distribution is construced
!     in each bin and these are translated according to the analytic
!     solutions
! (2) The translated tophats are then remapped to bins.  Here if a
!     top hat entirely or in part lies below the lowest bin it is
!     not counted.
!     .
! Additional notes (Peter Adams)
!     .
!     I have changed the routine to handle multicomponent aerosols.  The
!     arrays of mass moments are now two dimensional (size and species).
!     Only a single component (CSPECIES) is allowed to condense during
!     a given call to this routine.  Multicomponent condensation/evaporation
!     is accomplished via multiple calls.  Variables YLC and YUC are
!     similar to YL and YU except that they refer to the mass of the
!     condensing species, rather than total aerosol mass.
!     .
!     I have removed ventilation variables (VSW/VNTF) from the subroutine
!     call.  They still exist internally within this subroutine, but
!     are initialized such that they do nothing.
!     .
!     I have created a new variable, AMKDRY, which is the total mass in
!     a size bin (sum of all chemical components excluding water).  I
!     have also created WR, which is the ratio of total wet mass to
!     total dry mass in a size bin.
!     .
!     AMKC(k,j) is the total amount of mass after condensation of species
!     j in particles that BEGAN in bin k.  It is used as a diagnostic
!     for tracking down numerical errors.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: L,I,J,K,IMN
    REAL(fp)       :: DN,DM,DYI,XL,XU,YL,YLC,YU,YUC
    REAL(fp)       :: TEPS,NEPS,NEPS2,EX2,ZERO
    REAL(fp)       :: XI,XX,XP,YM,WTH,W1,W2,WW,AVG
    REAL(fp)       :: VSW,VNTF(ibins)
    REAL(fp)       :: TAU_L, maxtau

    REAL(fp)       :: AMKDRY(ibins), WR(ibins), AMKWET(ibins)
    REAL(fp)       :: AMKDRYSOL(ibins)

    LOGICAL        :: errspot !(win, 4/12/06)

    REAL(fp)       :: c1, c2 !correction factor (win, 5/25/06)
    REAL(fp)       :: madd(ibins) !condensing mass to be added by aqoxid
                                  !or SOAcond. For error fixing (win, 9/27/07)
    REAL(fp)       :: xadd(ibins) !mass per particle to be added by aqoxid
                                  ! or SOAcond. For error fixing (win, 9/27/07)
    REAL(fp)       :: macc !accumulating the condensing mass (win, 7/24/06)
    REAL(fp)       :: delt1,delt2 !the delta = mass not conserved (win, 7/24/06)
    REAL(fp)       :: dummy, xtra,maddtot ! for mass conserv fixing (win, 9/27/07)
    integer        :: kk !counter (wint, 7/24/06)
    REAL(fp)       :: AMKD_tot

    PARAMETER (TEPS=1.0e-40_fp,NEPS=1.0e-20_fp)
    PARAMETER (EX2=2.e+0_fp/3.e+0_fp,ZERO=0.0e+0_fp)
    PARAMETER (NEPS2=1.0e-10_fp)

    !=================================================================
    ! TMCOND begins here!
    !=================================================================

3   format(I4,200E20.11)

    !<step4.5> This first check cause the error of 'number not conserved'
    ! though only with the small amounts because when ANKD(k) = 0.e+0_fp from start,
    ! the original check just give it a value NEPS = 1.d-20, and then undergo
    ! tmcond calculation.   I'm changing the check to if ANKD(k)= 0.e+0_fp,
    ! then keep it that way and make the following calculations skip when
    ! ANKD(k) is zero (win, 10/18/05)

    ! If any ANKD are zero, set them to a small value to avoid division by zero
    !do k=1,ibins
    !   if (ANKD(k) .lt. NEPS) then
    !      ANKD(k)=NEPS
    !      AMKD(k,srtso4)=NEPS*1.4*xk(k) !make the added particles SO4
    !   endif
    !enddo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !<step5.1> Add print for debugging (win, 4/10/06)
    if (pdbug) then
       call debugprint(ANKD, AMKD, 0,0,0,'Entering TMCOND')
       ! print *, 'TMCOND:entering*************************'
       ! print *,'Nk(1:30)'
       ! print *, ANKD(1:30)
       ! print *,'Mk(1:30,comp)'
       ! do j=1,icomp
       ! print *,'comp',j
       ! print *, AMKD(1:30,j)
       ! enddo
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errspot = .false. !initialize error signal as false (Win, 4/12/06)

    !pja Sometimes, after repeated condensation calls, the average bin mass
    !pja can be just above the bin boundary - in that case, transfer a some
    !pja to the next highest bin
    !sfarina this is also true when small particles are growing really fast?
    !sfarina SOACOND throws thousands of errors for XI < 1
    !sfarina what this really means is AVG particle massfor bin k > XK(k+1)
    !sfarina through the debugger I found that mostly the difference is small
    do k=1,ibins-1
       if ( ANKD(k) .lt. NEPS2) goto 300 !<step4.5> (win, 10/18/05)
       ! Modify the check to include all dry mass (win, 10/3/08)
       AMKD_tot = 0.e+0_fp
       do kk=1,icomp-idiag
          AMKD_tot = AMKD_tot + AMKD(k,kk)
       enddo
       if ((AMKD_tot)/ANKD(k).gt.xk(k+1)) then
          !Prior to 10/3/08 (win)
          !if ((AMKD(k,srtso4))/ANKD(k).gt.xk(k+1)) then
          !sfarina: this does noting to help our avg mass per particle
          !         falling outside of bin boundaries:
          !         amkd_tot / ankd(k) = (amkd_tot * 0.9) / (ankd(k) * 0.9)
          !         we need to shift more mass than number.
          !         assuming we have some kind of distributionof particle sizes in bin K
          !         the largest ones will have more mass than average, so we can safly move
          !         more mass than number.
          !         that or we redistribute mass before SOAcond
          !
          do j=1,icomp-idiag
             AMKD(k+1,j)=AMKD(k+1,j)+0.1e+0_fp*AMKD(k,j)
             AMKD(k,j)=AMKD(k,j)*0.9e+0_fp
          enddo
          ANKD(k+1)=ANKD(k+1)+0.1e+0_fp*ANKD(k)
          ANKD(k)=ANKD(k)*0.9e+0_fp
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.1> Add print for debugging (win, 4/10/06)
          if (pdbug) then
             print *, 'Modified at checkpoint1: BIN',k
             print *,'ANKD(k)',ANKD(k),'ANKD(k+1)',ANKD(k+1)
             print *,'Mk(k,comp)       Mk(k+1,comp)'
             do j=1,icomp
                print *,'comp',j
                print *, AMKD(k,j), AMKD(k+1,j)
             enddo
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       endif
300    continue   !<step4.5> If aerosol number is zero (win, 10/18/05)

    enddo

    !pja Initialize ventilation variables so they don't do anything
    VSW=0.0e+0_fp
    DO L=1,ibins
       VNTF(L)=0.0e+0_fp
    ENDDO

    !pja Initialize AMKDRY and WR
    DO L=1,ibins
       AMKDRY(L)=0.e+0_fp
       AMKWET(L)=0.e+0_fp
       AMKDRYSOL(L) = 0.e+0_fp
       DO J=1,icomp-idiag     ! dry mass excl. nh4 (win, 9/26/08)
          AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
          ! Accumulate the absorbing media (win, 3/5/08)
          IF ( J == SRTOCIL  ) &
               AMKDRYSOL(L) = AMKDRYSOL(L) + AMKD(L,J)
       ENDDO
       DO J=1,ICOMP
          AMKWET(L) = AMKWET(L) + AMKD(L,J)
       ENDDO
       if (AMKDRY(L) .gt. 0.e+0_fp) &   !<step4.5> In case there is no mass, then just skip (win, 10/18/05)
            WR(L)= AMKWET(L) / AMKDRY(L)
    ENDDO

    !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(pdbug)then
       print*,'AMKDRY(1:ibins)'
       print *,AMKDRY(1:ibins)
       print *,'WR(1:ibins)'
       print *,WR(1:ibins)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !pja Initialize X() array of particle masses based on xk()
    DO L=1,ibins
       X(L)=xk(L)
    ENDDO

    !
    ! Only solve when significant forcing is available
    !
    maxtau=0.0e+0_fp
    do l=1,ibins
       maxtau=max(maxtau,abs(TAU(l)))
    enddo

    !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(pdbug) then
       print*,'tau(1:ibins)'
       print *,tau(1:ibins)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF(ABS(maxtau).LT.TEPS)THEN
       DO L=1,ibins
          DO J=1,icomp
             AMK(L,J)=AMKD(L,J)
          ENDDO
          ANK(L)=ANKD(L)
       ENDDO
    ELSE
       !<step5.3> Try to fix the error of mass conservation
       ! during aqueous oxidation. Too little mass is used up
       ! (win, 7/24/06)
       IF ( MAXVAL(MOXD(:)) >  0e+0_fp ) THEN
          IF( PDBUG ) PRINT *,'Mass_to_add_by_aqoxid_or_SOAcond'
          maddtot = 0e+0_fp
          DO L = 1, IBINS
             IF(TAU(L) >  0e+0_fp ) THEN
                MADD(L) = MOXD(L)
                XADD(L) = MOXD(L) / ANKD(L)
                !IF( CSPECIES == SRTSO4 ) THEN
                !   MADD(L) = MOXD * ANKD(L)  ! absolute condensing mass
                !   XADD(L) = MOXD            ! mass per particle
                !ELSE IF ( CSPECIES == SRTOCIL ) THEN
                !   MADD(L) = MOXD * AMKDRYSOL(L)
                !   XADD(L) = MADD(L) / ANKD(L)
                !ELSE
                !   PRINT *,'TMCOND ERROR : mass fixing not supported'
                !ENDIF
             ELSE
                MADD(L) = 0e+0_fp
                XADD(L) = 0e+0_fp
             ENDIF
             IF ( PDBUG ) PRINT *,L,madd(L), xadd(L)
             maddtot = maddtot + madd(L)
          ENDDO
       ENDIF

       DO L=1,ibins
          DO J=1,icomp
             AMK(L,J)=0.e+0_fp
          ENDDO
          ANK(L)=0.e+0_fp
       ENDDO
       WW=0.5e+0_fp
       ! IF(TAU.LT.0.)WW=.5e+0_fp
       !
       ! identify tophats and do lagrangian growth
       !
       DO L=1,ibins
          IF(ANKD(L) .LT. NEPS2)GOTO 200 !skip if Number is effectively zero

          !if tau is zero, leave everything in same bin
          IF (TAU(L) .EQ. 0.) THEN
             ANK(L)=ANK(L)+ANKD(L)
             DO J=1,icomp
                AMK(L,J)=AMK(L,J)+AMKD(L,J)
             ENDDO
          ENDIF
          IF (TAU(L) .EQ. 0.) GOTO 200

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.1> Add print for debugging (win, 4/10/06)
          if (pdbug) then
             print *, 'Identify_tophat_and_grow-BIN',L
             print *,'Starting_Nk(1:ibins)'
             print *, ANK(1:ibins)
             print *,'Starting_Mk(1:ibins,comp)'
             do j=1,icomp-1
                print *,'comp',j
                print *, AMK(1:ibins,j)
             enddo
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          !pja Limiting AVG, the average particle size to lie within the size
          !pja bounds causes particles to grow or shrink arbitrarily and is
          !pja wreacking havoc with choosing condensational timesteps and
          !pja conserving mass.  I have turned them off.
          !AVG=MAX(X(L),MIN(X(L+1),AMKDRY(L)/(NEPS+ANKD(L))))
          !try bring the above line back, win 4/10/06
          !win 4/10/06

          AVG=AMKDRY(L)/ANKD(L)
          XX=X(L)/AVG

#if  defined( TOMAS12 ) || defined( TOMAS15 )
          if(l.lt.ibins-1)then ! bin quadrupuling
             XI=.5e+0_fp + XX*(2.5e+0_fp - 2.0e+0_fp*XX)
             !XI<1 means the AVG falls out of bin bounds
             if (XI .LT. 1.e+0_fp) then
                !W1 will have sqrt of negative number
                write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                write(*,*)'AVG is ',AVG
                write(*,*)'Nk is ', ANKD(L)
                write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
                write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                errspot = .true.
                RETURN
             endif
             W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG/4.0e+0_fp ! cyhl 4.0=xk(k+1)/xk(k)
             W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0e+0_fp
          else ! final 2 bins mass*32
             XI=.5e+0_fp + XX*(16.5e+0_fp - 16.0e+0_fp*XX)
             if (XI .LT. 1.e+0_fp) then
                !W1 will have sqrt of negative number
                write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                write(*,*)'lower limit is',X(L)
                write(*,*)'AVG is ',AVG
                write(*,*)'Nk is ', ANKD(L)
                write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
                write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                errspot = .true.
                RETURN
             endif
             W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG/32.0e+0_fp ! cyhl 32.0=xk(k+1)/xk(k)
             W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0e+0_fp
          endif
#else
          XI=.5e+0_fp + XX*(1.5e+0_fp - XX)
          !XI<1 means the AVG falls out of bin bounds

          if (XI .LT. 1.e+0_fp) then
             !W1 will have sqrt of negative number
             write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
             write(*,*)'AVG is ',AVG
             write(*,*)'Nk is ', ANKD(L)
             write(*,*)'Mk are ', (AMKD(L,j),j=1,icomp)
             write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
             errspot = .true.
             RETURN
          endif
          W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG
          W2 =MIN(X(L+1)-AVG,AVG-X(L))
#endif

          WTH=W1*WW+W2*(1.e+0_fp-WW)
          IF(WTH.GT.1.) then
             write(*,*)'WTH>1 in cond, bin #',L
             errspot = .true.
             RETURN
          ENDIF

          XU=AVG+WTH*.5e+0_fp
          XL=AVG-WTH*.5e+0_fp
          ! Ventilation added bin-by-bin
          TAU_L=TAU(l)*MAX(1.e+0_fp,VNTF(L)*VSW)
          IF(TAU_L/TAU(l).GT. 6.) THEN
             PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
          ENDIF
          IF(TAU_L.GT.TAU(l)) THEN
             PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
          ENDIF
          ! prior to 5/25/06 (win)
          !YU=DMDT_INT(XU,TAU_L,WR(L))
          !YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
          !IF (YU .GT. X(ibins+1) ) THEN
          !   YUC=YUC*X(ibins+1)/YU
          !   YU=X(ibins+1)
          !ENDIF
          !YL=DMDT_INT(XL,TAU_L,WR(L))
          !YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
          !add new correction factor to YU and YL (win, 5/25/06)
          YU=DMDT_INT(XU,TAU_L,WR(L))
          YL=DMDT_INT(XL,TAU_L,WR(L))

          ! change to check MOXD of current bin (win, 10/3/08)
          IF( MOXD(L) == 0e+0_fp) THEN
             !Prior to 10/3/08 (win)
             !IF( MAXVAL(MOXD(:)) == 0e+0_fp ) THEN
             C1=1.e+0_fp          !for so4cond call, without correction factor.
          ELSE
             C1 = XADD(L)*2.e+0_fp/(YU+YL-XU-XL)
          ENDIF
          C2 = C1 - ( C1 - 1.e+0_fp ) * ( XU + XL )/( YU + YL )
          !prior to 10/2/08 (win)
          YU = YU * C2
          YL = YL * C2
          ! Run into a problem that YU < XU creating YUC<0
          ! So let's limit the application of C2 to only if
          ! it does not result in YU < XU and YL < XL (win, 10/2/08)
          !IF(TAU_L > 0.e+0_fp) YU = max( YU*C2, XU )
          !IF(TAU_L > 0.e+0_fp) YL = max( YL*C2, XL )

          !end part for fudging to get higher AVG

          YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
          IF (YU .GT. X(ibins+1) ) THEN
             !IF(.not.SPINUP(60.)) write(116,*) &
             !     'YU > Xk(30+1) ++++++++++++' !debug (win, 7/17/06)
             YUC=YUC*X(ibins+1)/YU
             YU=X(ibins+1)
             !errspot=.true.  !just try temp (win, 7/30/07)
          ENDIF
          YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
          DYI=1.e+0_fp/(YU-YL)

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.2> Debug why there is extra mass added when called
          ! by aqoxid. (win, 5/10/06)
          if (pdbug) then
             print *, 'XU',XU,'YU',YU,'YUC',YUC,'c2',c2
             print *, 'XL',XL,'YL',YL,'YLC',YLC
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          !deal with tiny negative (win, 5/28/06)
          if(YUC.lt.0e+0_fp .or. YLC.lt.0e+0_fp)then
             if(YLC.lt.0e+0_fp) YLC=0e+0_fp
             if(YUC.lt.0e+0_fp) then
                YUC = 0e+0_fp
                YLC = 0e+0_fp
             endif
             if(pdbug) print *,'Fudge negative YUC, YLC to zero'
          endif
          !
          ! deal with portion of distribution that lies below lowest gridpoint
          !
          IF(YL.LT.X(1))THEN

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !<step5.2> Debug step-by-step (win, 5/10/06)
             if (pdbug) print *,'YL<X(1)_Just_condensing_to_current_bin'
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !pja Instead of the following, I will just add all new condensed
             !pja mass to the same size bin
             !if ((YL/XL-1.e+0_fp) .LT. 1.e-3_fp) then
             !   !insignificant growth - leave alone
             !   ANK(L)=ANK(L)+ANKD(L)
             !   DO J=1,icomp-1
             !      AMK(L,J)=AMK(L,J)+AMKD(L,J)
             !   ENDDO
             !   GOTO 200
             !else
             !   !subtract out lower portion
             !   write(*,*)'ERROR in cond - low portion subtracted'
             !   write(*,*) 'Nk,Mk: ',ANKD(L),AMKD(L,1),AMKD(L,2)
             !   write(*,*) 'TAU: ', TAU_L
             !   write(*,*) 'XL, YL, YLC: ',XL,YL,YLC
             !   write(*,*) 'XU, YU, YUC: ',XU,YU,YUC
             !   ANKD(L)=ANKD(L)*MAX(ZERO,(YU-X(1)))*DYI
             !   YL=X(1)
             !   YLC=X(1)*AMKD(1,CSPECIES)/AMKDRY(1)
             !   DYI=1.e+0_fp/(YU-YL)
             !endif
             ANK(L)=ANK(L)+ANKD(L)
             do j=1,icomp
                if (J.EQ.CSPECIES) then
                   AMK(L,J)=AMK(L,J)+(YUC+YLC)*.5e+0_fp*ANKD(L)
                else
                   AMK(L,J)=AMK(L,J)+AMKD(L,J)
                endif
             enddo
             GOTO 200
          ENDIF
          IF(YU.LT.X(1))GOTO 200
          !
          ! Begin remapping (start search at present location if condensation)
          !
          IMN=1
          IF(TAU(l).GT.0.)IMN=L
          DO I=IMN,ibins
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !<step5.2> Debug step-by-step (win, 5/10/06)
             if(pdbug) print *,'Now_remapping_in_bin',I
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             IF(YL.LT.X(I+1))THEN
                ![1] lower bound of new tophat in the current I bin
                IF(YU.LE.X(I+1))THEN
                   ![2] upper bound of new tophat also in the current I bin
                   DN=ANKD(L)      ! DN = number from the bin L being remapped
                   do j=1,icomp
                      DM=AMKD(L,J)
                      IF (J.EQ.CSPECIES) THEN
                         !Add mass from new tophat to the existing mass of bin I
                         AMK(I,J)=(YUC+YLC)*.5e+0_fp*DN+AMK(I,J)
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> Debug step-by-step (win, 5/10/06)
                         if (pdbug) then
                            print *,'CASE_1:_New_Tophat_in_a_single_bin'
                            print *,'SO4_from_tophat=',(YUC+YLC)*.5e+0_fp*DN
                         endif
                         !<step5.3> Check mass conservation (win, 7/24/06)
                         if(MAXVAL(moxd(:)).gt.0e+0_fp)then
                            delt1 = (YUC+YLC)*.5e+0_fp*DN-AMKD(L,J)-madd(L)
                            if( abs(delt1)/madd(L).gt.1e-6_fp .and. &
                                 madd(L).gt.1e-4_fp)then
                               ! Just print out this for debugging
                               IF(.not.SPINUP(60.) .and. pdbug ) then
                                  !write(116,*)'CASE1_mass_conserv_fix'
                                  write(116,13) L, madd(L), delt1
13                                FORMAT('CASE_1 Bin ',I2,' moxid ', &
                                         E13.5,' delta ',E13.5 )
                                  !errspot=.true. !just try temp (win, 7/30/07)
                               ENDIF
                               AMK(I,J) = AMK(I,J)-delt1 !fix the error
                            endif
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         !For non-condensing, migrate the mass to bin I
                         AMK(I,J)=AMK(I,J)+DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> Debug step-by-step (win, 5/10/06)
                         if (pdbug) then
                            !print *,' Migrating_mass(',j,')',DM   !use this debugging line if there are more than seasalt+so4
                            print *,'Migrating_mass',DM
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   !Add number of old bin to ANK (which is blank for the first loop of bin I)
                   ANK(I)=ANK(I)+DN
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   !<step5.2> (win, 5/10/06)
                   if(pdbug) print*,'Migrating_number',DN
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ELSE
                   ![3] upper bound of new tophat grow beyond the upper bound of bin I
                   DN=ANKD(L)*(X(I+1)-YL)*DYI !DN= proportion of the number from tophat that still stays in the bin I
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   !<step5.2> (win, 5/10/06)
                   if ( pdbug) then
                      print*,'Case_2:_Tophat_cross_bin_boundary'
                      print *,'Number_that_remain_in_low_bin',DN
                   endif
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
                   !<step5.3> For fixing mass conserv problem (win, 7/24/06)
                   macc=0e+0_fp

                   do j=1,icomp
                      !DM= proporation of the mass that is still in bin I
                      DM=AMKD(L,J)*(X(I+1)-YL)*DYI
                      IF (J.EQ.CSPECIES) THEN
                         !XP= what would have grown to be X(I+1)
                         XP=DMDT_INT(X(I+1),-1.0e+0_fp*TAU_L,WR(L))
                         YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                         !add the condensing mass to the existing sulfate of bin I
                         AMK(I,J)=DN*(YM+YLC)*0.5e+0_fp+AMK(I,J)
                         !<step5.3>Accumulating the condensing mass for error check (win, 7/24/06)
                         macc = macc + DN*(YM+YLC)*0.5e+0_fp
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> (win, 5/10/06)
                         if(pdbug)then
                            print *,'XP',XP,'YM',YM
                            print *,'Cond_TophatLowEnd',DN*(YM+YLC)*0.5e+0_fp
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         !Add DM to AMK (which is blank for the first loop of bin I)
                         AMK(I,J)=AMK(I,J)+DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(pdbug) print*,'Other___in_low_end',DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   ANK(I)=ANK(I)+DN ! Add DN number to ANK (which is blank for the first loop of bin I)
                   ! Remapping loop from bin I+1 to bin30
                   DO K=I+1,ibins
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      if(pdbug) print *,'Spreading_to_bin',K
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      IF(YU.LE.X(K+1))GOTO 100
                      ![4] Found the bin where the high end of the tophat is in --> do the final loop

                      ![5.1] This part for distributing to the bins in between
                      !      the original and the furthest bin that growing occurs

                      !Use width of bin K to proportionate number from old bin wrt. to the top hat (YU-YL)
                      DN=ANKD(L)*(X(K+1)-X(K))*DYI

                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      if(pdbug) then
                         print *,'Number_migrated',DN
                      endif
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      do j=1,icomp
                         !Proportion of old-bin mass that falls in this current bin K
                         DM=AMKD(L,J)*(X(K+1)-X(K))*DYI
                         IF (J.EQ.CSPECIES) THEN
                            XP=DMDT_INT(X(K),-1.0e+0_fp*TAU_L,WR(L)) !what would have grown to be X(k)
                            YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !what would have grown to be X(k) but just for sulfate
                            AMK(K,J)=DN*1.5e+0_fp*YM+AMK(K,J)    ! A factor of 1.5 is from averaging (YM+2*YM)
                            !<step5.3> Accumulating condensing mass for error check (win, 7/24/06)
                            macc = macc+DN*1.5e+0_fp*YM
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            !<step5.2> (win, 5/10/06)
                            if(pdbug)then
                               print *,'XP',XP,'YM',YM
                               print *,'Cond_mass_spread',DN*1.5e+0_fp*YM
                            endif
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         ELSE
                            AMK(K,J)=AMK(K,J)+DM    !Add migrating mass of non-condensing species
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if(pdbug) print*,'No-cond_mass_migrate',DM
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         ENDIF
                      enddo
                      ANK(K)=ANK(K)+DN  !Add migrating number to the exising number of bin K
                   ENDDO
                   !This STOP is for when there's excessive growth over bin30
                   STOP 'Trying to put stuff in bin ibins+1'

100                CONTINUE
                   ![5.2] Final section that the tophat grows to.
                   DN=ANKD(L)*(YU-X(K))*DYI
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   if(pdbug) then
                      print *,'Found_right_edge_for_tophat'
                      print *,'Number_migrated',DN
                   endif
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   do j=1,icomp
                      DM=AMKD(L,J)*(YU-X(K))*DYI  ! proportion of old mass that gets to this furthest bin.
                      IF (J.EQ.CSPECIES) THEN
                         XP=DMDT_INT(X(K),-1.0e+0_fp*TAU_L,WR(L))   !what would have grown to be X(k)
                         YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !=XP for just sulfate
                         AMK(K,J)=DN*(YUC+YM)*0.5e+0_fp+AMK(K,J) !add condensing mass to existing sulfate of bin K
                         !<step5.3>Accumulating condensing mass for error check (win, 7/24/06)
                         macc = macc+DN*(YUC+YM)*0.5e+0_fp
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> (win, 5/10/06)
                         if(pdbug)then
                            print *,'XP',XP,'YM',YM
                            print *,'Cond_mass_spread_final', &
                                     DN*(YM+YUC)*0.5e+0_fp
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         AMK(K,J)=AMK(K,J)+DM  !This adds the migrating mass to the exising mass of non-condensing species
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(pdbug) print*,'No-cond_mass_migrated',DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   ANK(K)=ANK(K)+DN   !This adds the migrating number to the existing number of bin K

                   !<step5.3> Check mass conservation (win, 7/24/06)
                   if(MAXVAL(moxd(:)).gt.0e+0_fp)then
                      delt2 = 0e+0_fp
                      delt2 = macc-AMKD(L,CSPECIES)-madd(L)
                      if(abs(delt2)/ madd(L) > 1e-6)then
                         if( madd(L) > 10.e+0_fp .and. &
                             abs(delt2)/ madd(L) > 15e-2_fp ) then
                            !print *,'TMCOND ERROR: mass condensation', &
                            !  'discrep >15% during aqoxid or SOAcond'
                            IF(.not.SPINUP(60.))  THEN
14                             FORMAT('CASE_2 Bin',I2,' moxid',F7.1, &
                                      ' delta',F7.1 )
                               write(116,14) L, madd(L),delt2
                               !write(116,*)'CASE_2_mass_not_conserve'
                               !write(116,*)'For_bin',L,'moxid',madd(L) &
                               !     ,'delta',delt2
                            ENDIF
                            errspot=.true. !just try temp (win, 7/30/07)
                         endif !significant mass add (10 kg) - then print error.
                         !<step5.3> Fix the problem of mass not conserved
                         !in case of aqueous oxidation by find the missing mass
                         !and spread them equally into the bins that the final
                         !tophat has grown to. (win, 7/24/06)
                         xtra  = 0e+0_fp
                         dummy = 0e+0_fp
                         do kk = I,K
                            !AMK(kk,CSPECIES) = AMK(kk,CSPECIES)-delt2/(K-I+1)
                            dummy = AMK(kk,CSPECIES) - &
                                    ( delt2/(K-I+1) + xtra )
                            if(dummy < 0.e+0_fp )then
                               xtra = xtra + delt2/(K-I+1)
                            else
                               AMK(kk,CSPECIES) = dummy
                               xtra = 0.e+0_fp
                            endif
                         enddo
                      endif   !error>treshold
                   endif      !moxd>0

                ENDIF  !YU.LE.X(I+1)
                GOTO 200
             ELSE    !YL > X(I+1)
                IF(I == IBINS .and.(madd(L)/maddtot)> 1.5e-1_fp) THEN
11                 FORMAT( 'Tophat>Xk(31) at bin ',I3,' loosing ', &
                           E13.5,' kg = ',F5.1,'%')
                   if(MAXVAL(moxd(:)) > 0e+0_fp) then
                      print 11, L, madd(L),(madd(L)/maddtot)*1.e+2_fp
                      !write(116,11) L, madd(L),(madd(L)/maddtot)*1.e+2_fp
                      !write(117,*) madd(L)  !for accumulating mass loss
                      !PRINT *,'Tophat > Xk(31): growth over bin30,Loss%'
                      !if(moxd >0e+0_fp)print *,madd(L),(madd(L)/maddtot)*1.d2
                      !errspot = .true.
                   endif
                ENDIF
             ENDIF   !YL.LT.X(I+1)
          ENDDO !I loop
200       CONTINUE
       ENDDO    !L loop
    ENDIF

    !Signal error out to so4cond so the run can stop in aerophys and show i,j,l (win, 4/12/06)
    pdbug = errspot

    RETURN

  END SUBROUTINE TMCOND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodiag
!
! !DESCRIPTION: Subroutine AERODIAG saves changes to the appropriate diagnostic !  arrays (win, 7/23/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AERODIAG( PTYPE, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG_MOD,       ONLY : AD60_COND, AD60_COAG, AD60_NUCL
    USE DIAG_MOD,       ONLY : AD60_AQOX, AD60_ERROR, AD60_SOA
    USE DIAG_MOD,       ONLY : AD61,      AD61_INST
#endif
    USE ERROR_MOD,      ONLY : IT_IS_NAN
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: PTYPE    ! Number assigned to each dianostic
    INTEGER ,       INTENT(IN) :: I, J, L  ! Grid box indices
    REAL(fp),       INTENT(IN) :: Nk(IBINS)
    REAL(fp),       INTENT(IN) :: Nkd(IBINS)
    REAL(fp),       INTENT(IN) :: Mk(IBINS, ICOMP)
    REAL(fp),       INTENT(IN) :: Mkd(IBINS,ICOMP)
    REAL*4,         INTENT(IN) :: BOXVOL
    TYPE(GrdState), INTENT(IN) :: State_Grid ! Grid State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: K, JS
    REAL*4                :: ADXX(IBINS*(ICOMP-IDIAG+1))
    REAL*4,    SAVE       :: ACCUN, ACCUM(2)
    LOGICAL,   SAVE       :: FIRST = .TRUE.
    real*4                :: tempsum
    REAL*4                :: DTCHEM

    !=================================================================
    ! AERODIAG begins here!
    !=================================================================
#ifdef BPCH_DIAG

    ! PTYPE = 7 is for ND61  --- NOW use for Nucleation at species NK1
    !  Note: This is created to look at 3-D rate for a selected process
    !        Right now (5/21/08) I created this to watch NUCLEATION rate
    !        We can't afford to save all 30-bin and all mass component
    !        in all (I,J,L), thus this is created. (win, 5/21/08)
    IF ( PTYPE == 7 ) THEN
       IF ( L <= LD61 ) THEN
          DTCHEM = GET_TS_CHEM() ! chemistry time step in sec
          AD61(I,J,L,1) = AD61(I,J,L,1)  + ( NK(1) - NKD(1) )/ DTCHEM / BOXVOL  ! no./cm3/sec
          AD61_INST(I,J,L,1) =  ( NK(1) - NKD(1) ) /DTCHEM / BOXVOL ! no./cm3/sec

          !IF(i==39 .and. j==29 ) then
          !if ( AD61_INST(I,J,L) .gt. 1e18)  write(6,*) '*********', &
          !               'AD61_INST(',I,J,L,')', AD61_INST(I,J,L)
          !endif
       ENDIF
    ELSE ! PTYPE = 1-6 is for ND60

       ADXX(:) = 0e+0_fp
       IF ( FIRST ) THEN
          ACCUN = 0e0
          ACCUM(:) = 0e0
          FIRST = .FALSE.
       ENDIF

       ! Debug: check error fixed accumulated at each step
       !IF ( I == 1 .and. J == 1 .and. L == 1 .and. PTYPE == 2) then
       !   print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
       !   print *, '   Number :',ACCUN
       !   print *, '   Sulf   :',ACCUM(1)
       !   print *, '   NaCl   :',ACCUM(2)
       !ENDIF

       IF ( L <= LD60 ) THEN

          SELECT CASE ( PTYPE )
          CASE ( 1 )                ! Condensation diagnostic
             ADXX(:) = AD60_COND(1,J,L,:)

          CASE ( 2 )                ! Coagulation diagnostic
             ADXX(:) = AD60_COAG(1,J,L,:)

          CASE ( 3 )                ! Nucleation diagnostic
             ADXX(:) = AD60_NUCL(1,J,L,:)

          CASE ( 4 )                ! Aqueous oxidation diagnostic
             ADXX(:) = AD60_AQOX(1,J,L,:)

          CASE ( 5 )                ! Error fudging diagnostic
             ADXX(:) = AD60_ERROR(1,J,L,:)

          CASE ( 6 )                ! SOA condensation diagnostic
             ADXX(:) = AD60_SOA(1,J,L,:)
          END SELECT

          ! Change of aerosol number
          DO K = 1, IBINS
             ADXX(K) =  ADXX(K) + NK(K) - NKD(K)
             !IF ( PTYPE == 2 ) ACCUN = ACCUN + NK(K) - NKD(K)
          ENDDO
          IF ( IT_IS_NAN(ACCUN)) print *,'AERODIAG: Nan',I,J,L

          ! Change of aerosol mass
          DO JS = 1, ICOMP-IDIAG
             tempsum = 0e0
             DO K = 1, IBINS
                ADXX(JS*IBINS+K) = ADXX(JS*IBINS+K) + MK(K,JS) - MKD(K,JS)
                !tempsum = tempsum + MK(K,JS) - MKD(K,JS)
                !IF (PTYPE == 2 ) ACCUM(JS) = ACCUM(JS) + MK(K,JS) - MKD(K,JS)
             ENDDO
          ENDDO
          !IF ( IT_IS_NAN(ACCUM(1))) print *,'ADIAG: Nan',I,J,L

          ! Put the updated values back into the diagnostic arrays
          SELECT CASE ( PTYPE )
          CASE ( 1 )                ! Condensation diagnostic
             AD60_COND(1,J,L,:) = ADXX(:)

          CASE ( 2 )                ! Coagulation diagnostic
             AD60_COAG(1,J,L,:) = ADXX(:)

          CASE ( 3 )                ! Nucleation diagnostic
             AD60_NUCL(1,J,L,:) = ADXX(:)

          CASE ( 4 )                ! Aqueous oxidation diagnostic
             AD60_AQOX(1,J,L,:) = ADXX(:)

          CASE ( 5 )                ! Error fudging diagnostic
             AD60_ERROR(1,J,L,:) = ADXX(:)

          CASE ( 6 )                ! SOA condensation diagnostic
             AD60_SOA(1,J,L,:) = ADXX(:)
          END SELECT

       ENDIF
       ! Debug: check error fixed accumulated at each step
       !IF ( I == 3 .and. J == 41 .and. L == 30 .and. PTYPE == 5) then
       !   print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
       !   print *, '   Number :',ACCUN
       !   print *, '   Sulf   :',ACCUM(1)
       !   print *, '   NaCl   :',ACCUM(2)
       !   print *, ' Nk',Nk(:)
       !   print *, ' Nkd',Nkd(:)
       !   print *, ' Mk',Mk(:,1)
       !   print *, ' Mkd',Mkd(:,1)
       !   print *, ' Mk',Mk(:,2)
       !   print *, ' Mkd',Mkd(:,2)
       !ENDIF

       ! Debug: check error fixed accumulated at each step
       IF ( I == State_Grid%NX .and. J == State_Grid%NY .and. &
            L == State_Grid%NZ .and. PTYPE == 2 ) then
          !print *, 'Accumulated diagnostic for ND60 #',PTYPE,' at',i,j,l
          print *, ' Accumulated Coagulation'
          print *, '   Number :',ACCUN
          print *, '   Sulf   :',ACCUM(1)
          print *, '   NaCl   :',ACCUM(2)
       ENDIF

    ENDIF ! If (PTYPE == 7)

#endif

  END SUBROUTINE AERODIAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tomas
!
! !DESCRIPTION: Subroutine INIT_TOMAS intializes variables for TOMAS
!  microphysics based on switches from input.geos, e.g. what aerosol species to
!  simulate.(win, 7/9/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TOMAS( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR, ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE inquireMod,         ONLY : findFreeLUN
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: AS, K, I,J,L, dum1, dum2, dum3, LUN
    REAL(fp)             :: Mo
    CHARACTER(LEN=255)   :: filename
    CHARACTER(LEN=255)   :: fname(4)
    CHARACTER(LEN=255)   :: DATA_DIR

    !=================================================================
    ! INIT_TOMAS begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Define species indices here, these are now saved as
    ! module variables (bmy, 6/20/16)
    id_NK1   = Ind_('NK1'  )
    id_H2SO4 = Ind_('H2SO4')
    id_AW1   = Ind_('AW1'  )
    id_SF1   = Ind_('SF1'  )
    id_SO4   = Ind_('SO4'  )
    id_NH3   = Ind_('NH3'  )
    id_NH4   = Ind_('NH4'  )
    id_SF1   = Ind_('SF1'  )
    id_SS1   = Ind_('SS1'  )
    id_ECIL1 = Ind_('ECIL1')
    id_ECOB1 = Ind_('ECOB1')
    id_OCIL1 = Ind_('OCIL1')
    id_OCOB1 = Ind_('OCOB1')
    id_DUST1 = Ind_('DUST1')

    ! Now read large TOMAS input files from a common disk directory
    ! (bmy, 1/30/14)
    DATA_DIR = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'TOMAS_201402/'

    ! Define subgrid coagulation timescale (win, 10/28/08)
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN
       SGCTSCALE = 10.*3600.  ! 10 hours
    ELSE IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       SGCTSCALE = 5.*3600.
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       SGCTSCALE = 1.*3600.
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       SGCTSCALE = 0.5*3600.
    ENDIF

#if    defined( TOMAS40 )
    Mo = 1.0e-21_fp*2.e+0_fp**(-10)
#elif  defined( TOMAS15 )
    Mo = 1.0e-21_fp*4.e+0_fp**(-3)
#else
    Mo = 1.0e-21_fp
#endif

    ICOMP = 0
    IDIAG = 0
    K = 0
    !IF( LSULF30 )  THEN
    IF (id_SF1 > 0) THEN
       ICOMP = ICOMP + 1
       !SRTSO4 = ICOMP
    ENDIF
    !IF( LSALT30 )  THEN
    IF ( id_SS1 > 0 ) THEN
       ICOMP = ICOMP + 1
       !SRTNACL = ICOMP
    ENDIF
    !IF( LCARB30 )  THEN
    IF ( id_ECIL1 > 0 .AND. id_ECOB1 > 0 .AND. &
         id_OCIL1 > 0 .AND. id_OCOB1 > 0 ) THEN
       ICOMP = ICOMP + 1
       !SRTECIL = ICOMP
       ICOMP = ICOMP + 1
       !SRTECOB = ICOMP
       ICOMP = ICOMP + 1
       !SRTOCIL = ICOMP
       ICOMP = ICOMP + 1
       !SRTOCOB = ICOMP
    ENDIF
    !IF( LDUST30 )  THEN
    IF ( id_DUST1 > 0 ) THEN
       ICOMP = ICOMP + 1
       !SRTDUST = ICOMP
    ENDIF

    ! Have to add one more for aerosol water
    IF( ICOMP > 1 ) THEN
       ICOMP = ICOMP + 1
       IDIAG = IDIAG + 1
       !SRTNH4 = ICOMP

       ICOMP = ICOMP + 1
       IDIAG = IDIAG + 1
       !SRTH2O = ICOMP
    ENDIF
    print *, 'In init_TOMAS, ICOMP = ', ICOMP
    print *, 'In init_TOMAS, IBINS = ', IBINS

    !=================================================================
    ! Allocate arrays
    !=================================================================

    ALLOCATE( MOLWT( ICOMP ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLWT [TOMAS]' )
    MOLWT(:) = 0e+0_fp
    
    ALLOCATE( H2SO4_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2SO4_RATE' )
    H2SO4_RATE = 0.0e+0_fp

    !=================================================================
    ! Calculate aerosol size bin boundaries (dry mass / particle)
    !=================================================================

#if  defined( TOMAS12 ) || defined( TOMAS15 )
    DO K = 1, IBINS + 1
       if(k.lt.ibins)then
          xk(k)=Mo * 4.e+0_fp**(k-1) !mass quadrupling
       else
          xk(k)=xk(k-1) * 32.e+0_fp
       endif
    ENDDO
#else
    DO K = 1, IBINS + 1
       Xk( k ) = Mo * 2.e+0_fp ** ( K-1 )
    ENDDO
#endif

    DO K = 1, IBINS
       AVGMASS( k ) = sqrt(Xk(k)*Xk(k+1))
    ENDDO

    DO J = 1, ICOMP
       IF ( J == SRTSO4 ) THEN
          MOLWT(J) = 98.0
       ELSE IF ( J == SRTNACL ) THEN
          MOLWT(J) = 58.5
       ELSE IF ( J == SRTH2O ) THEN
          MOLWT(J) = 18.0
       ELSE IF ( J == SRTECIL ) THEN
          MOLWT(J) = 12.0
       ELSE IF ( J == SRTECOB ) THEN
          MOLWT(J) = 12.0
       ELSE IF ( J == SRTOCIL ) THEN
          MOLWT(J) = 12.0
       ELSE IF ( J == SRTOCOB ) THEN
          MOLWT(J) = 12.0
       ELSE IF ( J == SRTDUST ) THEN
          MOLWT(J) = 100.0
       ELSE IF ( J == SRTNH4 ) THEN
          MOLWT(J) = 18.0
       ELSE
          PRINT *,'INIT_TOMAS ERROR: Modify code for more species!!'
          CALL ERROR_STOP('INIT_TOMAS','Modify code for new species')
       ENDIF
    ENDDO

    !=================================================================
    ! Create a look-up table for activating bin and scavenging
    ! fraction as a function of chemical composition.
    !=================================================================
    IF ( IBINS == 12 .OR. IBINS == 15 ) THEN
       fname(1) = TRIM( DATA_DIR ) // 'binact02_12.dat'
       fname(2) = TRIM( DATA_DIR ) // 'binact10_12.dat'
       fname(3) = TRIM( DATA_DIR ) // 'fraction02_12.dat'
       fname(4) = TRIM( DATA_DIR ) // 'fraction10_12.dat'
    ELSE IF ( IBINS == 30 .OR. IBINS == 40 ) THEN
       fname(1) = TRIM( DATA_DIR ) // 'binact02.dat'
       fname(2) = TRIM( DATA_DIR ) // 'binact10.dat'
       fname(3) = TRIM( DATA_DIR ) // 'fraction02.dat'
       fname(4) = TRIM( DATA_DIR ) // 'fraction10.dat'
    END IF

    CALL READBINACT  ( fname(1), BINACT1   )
    CALL READBINACT  ( fname(2), BINACT2   )
    CALL READFRACTION( fname(3), FRACTION1 )
    CALL READFRACTION( fname(4), FRACTION2 )

    !initialize yu lookup table
    call READJIMN5D( Input_Opt, RC )  ! yu nucleation inputs

    ! Find a free file LUN
    LUN = findFreeLUN()

    ! Read cosmic ray ion input file
    FILENAME = TRIM( DATA_DIR ) // 'IonPairs1GV.4x5'
    OPEN( unit=LUN, FILE=TRIM( FILENAME ), FORM='FORMATTED', STATUS='OLD' )

    DO L=1,9
    DO J=1,46
    DO I=1,72
       READ( LUN ,'(I5,I5,I5,E10.3)') dum1,dum2,dum3,cosmic_ions(I,J,L)
       if (I.eq.50.and.J.eq.20.and.L.eq.5)then
          print*,'ion test',cosmic_ions(I,J,L)
       endif
    ENDDO
    ENDDO
    ENDDO

    !carbon emission factors:

    ! Close file
    CLOSE( LUN )

  END SUBROUTINE INIT_TOMAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readbinact
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READBINACT( INFILE, BINACT )
!
! !USES:
!
    USE inquireMod, ONLY : findFreeLun
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=255) INFILE
!
! !OUTPUT PARAMETERS:
!
    INTEGER BINACT(101,101,101)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER INNUM, II, JJ, KK, LUN

    ! Find a free file LUN
    INNUM = findFreeLun()

1   FORMAT(I2)
    OPEN(UNIT=INNUM,FILE=INFILE,FORM='FORMATTED',STATUS='OLD')
    DO II=1,101
    DO JJ=1,101
    DO KK=1,101
       READ(INNUM,1) BINACT(KK,JJ,II)
       IF (BINACT(KK,JJ,II).eq.0) BINACT(KK,JJ,II)=IBINS + 1
    ENDDO
    ENDDO
    ENDDO
    CLOSE(INNUM)

    RETURN

  END SUBROUTINE READBINACT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: readfraction
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READFRACTION(INFILE,FRACTION)
!
! !USES:
!
    USE inquireMod, ONLY : findFreeLun
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=255) INFILE
!
! !OUTPUT PARAMETERS:
!
    REAL(fp) FRACTION(101,101,101)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER INNUM, II, JJ, KK

    ! Find a free file LUN
    INNUM = findFreeLun()

1   FORMAT(F6.5)
    OPEN(UNIT=INNUM,FILE=INFILE,FORM='FORMATTED',STATUS='OLD')
    DO II=1,101
    DO JJ=1,101
    DO KK=1,101
       READ(INNUM,1) FRACTION(KK,JJ,II)
       IF (FRACTION(KK,JJ,II).GT.1.) FRACTION(KK,JJ,II)=0.
    ENDDO
    ENDDO
    ENDDO
    CLOSE(INNUM)

    RETURN

  END SUBROUTINE READFRACTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getfraction
!
! !DESCRIPTION: Subroutine GETFRACTION calculate the mass fraction of each
!  soluble component i.e. SO4, sea-salt, hydrophilic OC to use as inputs for a
!  lookup table of activating bin and scavenging fraction. (win, 9/10/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GETFRACTION( I, J, L, N, LS, State_Chm, State_Grid, &
                          State_Met, FRACTION, SOLFRAC )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box index
    INTEGER,        INTENT(IN)    :: N            ! Species ID
    LOGICAL,        INTENT(IN)    :: LS           ! True=LS (stratiform) precip,
                                                  ! False= convective precip
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    !FRACTION         : Scavenging fraction of the given grid box
    !SOLFRAC          : Soluble mass fraction of the aerosol popultion of the
    !                   given grid box
    REAL(fp),       INTENT(OUT)   :: FRACTION, SOLFRAC
!
! !REMARKS:
!  This routine is called from the convection routines (via wetscav_mod.F90
!  routines COMPUTE_F and DO_RAINOUT_ONLY. (bmy, 7/18/16)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL*4                 ::  MECIL, MOCIL, MOCOB, MSO4, MNACL, MTOT
    REAL*4                 ::  MECOB, MDUST
    REAL*4                 ::  XOCIL, XSO4, XNACL
    INTEGER                ::  ISO4, INACL, IOCIL
    INTEGER                ::  GETBINACT
    INTEGER                ::  BIN
    INTEGER                ::  OFFSET
    CHARACTER(LEN=255)     ::  MSG, LOC
    REAL(fp)               ::  UNITFACTOR

    ! Pointers
    REAL(fp), POINTER      :: Spc(:,:,:,:)

    !=================================================================
    ! GETFRACTION begins here
    !=================================================================

    ! Determine factor used to convert Spc to units of [kg] locally
    ! (ewl, 9/29/15)
    IF ( TRIM( State_Chm%Spc_Units ) == 'kg/m2' ) THEN
       UNITFACTOR = State_Grid%AREA_M2(I,J)

    ELSE IF (  TRIM( State_Chm%Spc_Units ) == 'kg/kg dry' ) THEN
       UNITFACTOR = State_Met%AD(I,J,L)

    ELSE
       MSG = 'Unexpected species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine GETFRACTION in tomas_mod.F90'
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

    ! Point to chemical species array
    ! Units are now [kg/m2] in wet deposition and [kg/kg total air] in
    ! convection (ewl, 9/29/15)
    Spc => State_Chm%Species

    BIN = N - id_NK1 + 1
    IF ( BIN > IBINS ) THEN
       BIN = MOD( BIN, IBINS )
       IF ( BIN == 0 ) BIN = IBINS
    ENDIF

    MECIL = 0.E0
    MOCIL = 0.E0
    MOCOB = 0.E0
    MSO4  = 0.E0
    MNACL = 0.E0
    MDUST = 0.E0

    IF ( id_ECIL1 > 0 .AND.id_OCIL1 > 0 .AND. id_OCOB1 > 0 ) THEN
       MECIL = Spc(I,J,L,id_ECIL1-1+BIN) * UNITFACTOR
       MOCIL = Spc(I,J,L,id_OCIL1-1+BIN) * UNITFACTOR
       MOCOB = Spc(I,J,L,id_OCOB1-1+BIN) * UNITFACTOR
    ENDIF
    IF ( id_DUST1 > 0 ) MDUST = Spc(I,J,L,id_DUST1-1+BIN) * UNITFACTOR
    !account for ammonium sulfate
    IF ( id_SF1 > 0 ) MSO4  = Spc(I,J,L,id_SF1-1+BIN) * 1.2 * UNITFACTOR
    IF ( id_SS1 > 0 ) MNACL = Spc(I,J,L,id_SS1-1+BIN) * UNITFACTOR
    MTOT  = MECIL + MOCIL + MOCOB + MSO4 + MNACL + MDUST + 1.e-20
    XOCIL = MOCIL / MTOT
    XSO4  = MSO4  / MTOT
    XNACL = MNACL / MTOT
    ISO4  = MIN(101, INT(XSO4*100)+1)
    INACL = MIN(101, INT(XNACL*100)+1)
    IOCIL = MIN(101, INT(XOCIL*100)+1)

    !==========================================================
    ! subroutine was written considering bin 1 is 10nm
    ! in TOMAS-40, bin 1 is 1nm and bin 11 is 10nm
    !==========================================================
#if    defined( TOMAS40 )
    OFFSET = 10
#elif  defined( TOMAS15 )
    OFFSET = 3
#else
    OFFSET = 0
#endif

    IF ( LS ) THEN
       GETBINACT = BINACT1(ISO4, INACL, IOCIL) + OFFSET
    ELSE
       GETBINACT = BINACT2(ISO4, INACL, IOCIL) + OFFSET
    ENDIF

    if((GETBINACT.lt.0).or.(GETBINACT.gt.50))then
       print*,'BINACT ERROR GETBINACT=',GETBINACT
       stop
    endif

    !print*,'N, BINACT = ',N,GETBINACT

    IF ( GETBINACT > BIN ) THEN
       FRACTION = 0. !NOT ACTIVATED
    ELSE IF ( GETBINACT == BIN ) THEN
       IF ( LS ) THEN
          FRACTION = FRACTION1(ISO4, INACL, IOCIL ) !PARTLY ACTIVATED
       ELSE
          FRACTION = FRACTION2(ISO4, INACL, IOCIL ) !PARTLY ACTIVATED
       ENDIF
    ELSE
       FRACTION = 1. !ALL ACTIVATED
    ENDIF

    ! Calculate the soluble fraction of mass
    MECOB = 0.E0
    IF ( id_ECOB1 > 0 ) MECOB = Spc(I,J,L,id_ECOB1-1+BIN) * UNITFACTOR
    SOLFRAC = MTOT / ( MTOT + MECOB )

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE GETFRACTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getactbin
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GETACTBIN ( I, J, L, N, LS, BINACT, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box index
    INTEGER,        INTENT(IN)    :: N            ! Species ID
    LOGICAL,        INTENT(IN)    :: LS           ! True=LS (stratiform) precip,
                                                  ! False= convective precip
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: BINACT
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*4                 ::  MECIL, MOCIL, MOCOB, MSO4, MNACL, MTOT
    REAL*4                 ::  MECOB, MDUST
    REAL*4                 ::  XOCIL, XSO4, XNACL
    INTEGER                ::  ISO4, INACL, IOCIL
    INTEGER                ::  BIN
    INTEGER                ::  OFFSET
    CHARACTER(LEN=255)     :: MSG, LOC ! For species unit check (ewl)

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    !=================================================================
    ! GETACTBIN begins here
    !=================================================================

    ! Assume success
    RC                =  GC_SUCCESS

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine GETACTBIN in tomas_mod.F90'
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    BIN = N - id_NK1 + 1
    IF ( BIN > IBINS ) THEN
       BIN = MOD( BIN, IBINS )
       IF ( BIN == 0 ) BIN = IBINS
    ENDIF

    MECIL = 0.E0
    MOCIL = 0.E0
    MOCOB = 0.E0
    MSO4  = 0.E0
    MNACL = 0.E0

    IF ( id_ECIL1 > 0 .AND.id_OCIL1 > 0 .AND. id_OCOB1 > 0 ) THEN
       !IF (LCARB30) THEN
       MECIL = Spc(I,J,L,id_ECIL1-1+BIN)
       MOCIL = Spc(I,J,L,id_OCIL1-1+BIN)
       MOCOB = Spc(I,J,L,id_OCOB1-1+BIN)
    ENDIF
    IF ( id_DUST1 > 0 ) MDUST = Spc(I,J,L,id_DUST1-1+BIN)
    !IF (LDUST30) MDUST = Spc(I,J,L,id_DUST1-1+BIN)
    MSO4  = Spc(I,J,L,id_SF1-1+BIN) * 1.2 !account for ammonium sulfate
    MNACL = Spc(I,J,L,id_SS1-1+BIN)

    MTOT  = MECIL + MOCIL + MOCOB + MSO4 + MNACL + MDUST + 1.e-20
    XOCIL = MOCIL / MTOT
    XSO4  = MSO4 / MTOT
    XNACL = MNACL / MTOT
    ISO4  = MIN(101, INT(XSO4*100)+1)
    INACL = MIN(101, INT(XNACL*100)+1)
    IOCIL = MIN(101, INT(XOCIL*100)+1)

    !==========================================================
    ! subroutine was written considering bin 1 is 10nm
    ! in TOMAS-40, bin 1 is 1nm and bin 11 is 10nm
    !==========================================================
#if   defined( TOMAS40 )
    OFFSET = 10
#elif  defined( TOMAS15 )
    OFFSET = 3
#else
    OFFSET = 0
#endif

    IF ( LS ) THEN
       BINACT = BINACT1(ISO4, INACL, IOCIL) + OFFSET
    ELSE
       BINACT = BINACT2(ISO4, INACL, IOCIL) + OFFSET
    ENDIF

    ! Free pointer
    Spc => NULL()

    RETURN

  END SUBROUTINE GETACTBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ezwatereqm
!
! !DESCRIPTION: WRITTEN BY Peter Adams, March 2000
!     .
!     This routine uses the current RH to calculate how much water is
!     in equilibrium with the aerosol.  Aerosol water concentrations
!     are assumed to be in equilibrium at all times and the array of
!     concentrations is updated accordingly.
!     .
!     Introduced to GEOS-CHEM by Win Trivitayanurak. May 8, 2006.
!     This file is replacing the old ezwatereqm that was not compatible
!     with multicomponent aerosols.  The new ezwatereqm use external
!     functions to do ISORROPIA-result curve fitting for each aerosol
!     component.
!     WARNING :
!      *** Watch out for the new aerosol species added in the future!
!     .
!     This version of the routine works for sulfate and sea salt
!     particles.  They are assumed to be externally mixed and their
!     associated water is added up to get total aerosol water.
!     wr is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fits to estimate wr based on the current humidity.
!     The curve fit is based on ISORROPIA results for ammonium bisulfate
!     at 273 K and sea salt at 273 K.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZWATEREQM( Mke, RHTOMAS )
!
! !INPUT PARAMETERS:
!
    REAL*4,  INTENT(IN)    :: RHTOMAS
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),INTENT(INOUT) :: Mke(IBINS,ICOMP)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: k, j
    REAL(fp)            :: so4mass, naclmass, ocilmass
    REAL(fp)            :: wrso4, wrnacl, wrocil
    REAL(fp)            :: rhe

    !========================================================================
    ! EZWATEREQM begins here!
    !========================================================================

    rhe=100.e+0_fp*rhtomas
    if (lowRH == 1) THEN !JKodros RH switch
       if (rhe .gt. 90.e+0_fp) rhe=90.e+0_fp
    ELSE
       if (rhe .gt. 99.e+0_fp) rhe=99.e+0_fp
    END IF !JKodros RH switch
    if (rhe .lt. 1.e+0_fp) rhe=1.e+0_fp

    do k=1,ibins

       so4mass=Mke(k,srtso4)*1.2  !1.2 converts kg so4 to kg nh4hso4
       wrso4=waterso4(rhe)

       ! Add condition for srtnacl in case of running so4 only. (win, 5/8/06)
       if (srtnacl.gt.0) then
          naclmass=Mke(k,srtnacl) !already as kg nacl - no conv necessary
          wrnacl=waternacl(rhe)
       else
          naclmass = 0.e+0_fp
          wrnacl = 1.e+0_fp
       endif

       if (srtocil.gt.0) then
          ocilmass=Mke(k,srtocil) !already as kg ocil - no conv necessary
          wrocil=waterocil(rhe)
       else
          ocilmass = 0.e+0_fp
          wrocil = 1.e+0_fp
       endif

       Mke(k,srth2o)=so4mass*(wrso4-1.e+0_fp)+naclmass &
                     *(wrnacl-1.e+0_fp) &
                     +ocilmass*(wrocil-1.e+0_fp)

    enddo

    RETURN

  END SUBROUTINE EZWATEREQM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ezwatereqm2
!
! !DESCRIPTION: Subroutine EZWATEREQM2 is just like EZWATEREQM but access
!  directly to STT array unlike EZWATEREQM that needs the array Mke to be
!  passed in and out. This subroutine is for calling from outside microphysics
!  module. (win, 8/6/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZWATEREQM2( I, J, L, BIN, State_Met, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L, BIN
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)               :: RHE
    REAL(fp)               :: SO4MASS, NACLMASS, OCILMASS
    REAL(fp)               :: WRSO4, WRNACL, WROCIL
    CHARACTER(LEN=255)     :: MSG, LOC ! For species unit check (ewl)

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    !========================================================================
    ! EZWATEREQM2 begins here!
    !========================================================================

    ! Assume success
    RC  = GC_SUCCESS

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    rhe = State_Met%RH(i,j,l)             !RH [=] percent

    if (rhe .gt. 99.) rhe=99.
    if (rhe .lt. 1.) rhe=1.

    so4mass=Spc(i,j,l,id_SF1-1+bin)*1.2 !1.2 converts kg so4 to kg nh4hso4
    wrso4=waterso4(rhe)       !use external function

    ! Add condition for srtnacl in case of running so4 only. (win, 5/8/06)
    if (id_SS1.gt.0) then
       naclmass=Spc(i,j,l,id_SS1-1+bin) !already as kg nacl - no conv necessary
       wrnacl=waternacl(rhe)  !use external function
    else
       naclmass = 0.e+0_fp
       wrnacl = 1.e+0_fp
    endif

    if (id_OCIL1 > 0) then
       ocilmass=Spc(i,j,l,id_OCIL1-1+bin)  !already as kg ocil - no conv necessary
       wrocil=waterocil(rhe)
    else
       ocilmass = 0.e+0_fp
       wrocil = 1.e+0_fp
    endif

    Spc(i,j,l,id_AW1-1+bin)= so4mass*(wrso4-1.e+0_fp) + &
                             naclmass*(wrnacl-1.e+0_fp) &
                             + ocilmass*(wrocil-1.e+0_fp)

    ! Free pointer
    Spc => NULL()

    RETURN

  END SUBROUTINE EZWATEREQM2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: eznh3eqm
!
! !DESCRIPTION: Subroutine EZNH3REQM2 puts ammonia to the particle phase until
!  there is 2 moles of ammonium per mole of sulfate and the remainder
!  of ammonia is left in the gas phase. (win, 9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZNH3EQM( Gce, Mke )
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT)  :: Gce(icomp - 1) !sfarina - fixed incorrect definition of Gc array
    REAL(fp),  INTENT(INOUT)  :: Mke(ibins,icomp)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer       ::  k
    REAL(fp)        :: tot_nh3  !total kmoles of ammonia
    REAL(fp)        :: tot_so4  !total kmoles of so4
    REAL(fp)        :: sfrac    !fraction of sulfate that is in that bin

    !========================================================================
    ! EZNH3EQM begins here!
    !========================================================================

    ! get the total number of kmol nh3
    tot_nh3 = Gce(srtnh4)/17.e+0_fp
    do k=1,ibins
       tot_nh3 = tot_nh3 + Mke(k,srtnh4)/18.e+0_fp
    enddo

    ! get the total number of kmol so4
    tot_so4 = 0.e+0_fp
    do k=1,ibins
       tot_so4 = tot_so4 + Mke(k,srtso4)/96.e+0_fp
    enddo

    ! see if there is free ammonia
    if (tot_nh3/2.e+0_fp.lt.tot_so4)then  ! no free ammonia
       Gce(srtnh4) = 0.e+0_fp ! no gas phase ammonia
       do k=1,ibins
          sfrac = Mke(k,srtso4)/96.e+0_fp/tot_so4
          Mke(k,srtnh4) = sfrac*tot_nh3*18.e+0_fp ! put the ammonia where the sulfate is
          ! Debug
          !if ( Mke(k,srtnh4) < 0.0 ) then
          !   print *,'negative gas phase ammonia in eznh3eqm!!'
          !   print *,'bin  ', k
          !endif
       enddo
    else ! free ammonia
       do k=1,ibins
          Mke(k,srtnh4) = Mke(k,srtso4)/96.e+0_fp*2.e+0_fp*18.e+0_fp ! fill the particle phase
          ! Debug
          !if ( Mke(k,srtnh4) < 0.0 ) then
          !   print *,'negative gas phase ammonia in eznh3eqm!!'
          !   print *,'bin  ', k
          !endif
       enddo
       Gce(srtnh4) = (tot_nh3 - tot_so4*2.e+0_fp)*17.e+0_fp ! put whats left over in the gas phase
       ! Debug
       !if ( Gce(srtnh4) < 0.0 ) then
       !   print *,'negative gas phase ammonia in eznh3eqm!!'
       !endif

    endif

    RETURN

  END SUBROUTINE EZNH3EQM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aero_diaden
!
! !DESCRIPTION: AERO/_DIADEN calculate the diameter and density by calling
!  external functions GETDP and AERODENS respectively. (win, 7/19/07)
!  Note: This subroutine is created for supplying diameter and density for
!        dry dep velocity calculation in DEPVEL.  Did not want to add much
!        to DEPVEL.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AERO_DIADEN( LEV, Input_Opt, State_Chm, State_Grid, State_Met, &
                          DIA, DENSITY, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE State_Grid_Mod,     ONLY : GrdState
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)    :: LEV
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState),    INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),    INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),          INTENT(OUT)   :: DIA(State_Grid%NX,State_Grid%NY,IBINS)
    REAL(fp),          INTENT(OUT)   :: DENSITY(State_Grid%NX,State_Grid%NY,IBINS)
    INTEGER,           INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I,J, BIN, JC, TRACID, WID
    INTEGER             :: N, NA, nAdvect
    REAL(fp)            :: MSO4, MNACL, MH2O
    REAL(fp)            :: MECIL, MECOB, MOCIL, MOCOB, MDUST
    CHARACTER(LEN=255)  :: MSG, LOC
    CHARACTER(LEN=63)   :: OrigUnit

    ! Arrays
    REAL(fp)            :: AMASS(ICOMP)

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    !========================================================================
    ! AERO_DIADEN begins here!
    !========================================================================

    ! Assume success
    RC =  GC_SUCCESS

    ! Convert species units to [kg] for this routine.
    ! NOTE: For complete area-independence, species units will need to be
    !       mixing ratio or mass per unit area in TOMAS
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Routine AERO_DIADEN in tomas_mod.F90')
       RETURN
    ENDIF

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine AERO_DIADEN in tomas_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
    ENDIF

    ! Point to the chemical species array [kg]
    Spc => State_Chm%Species

    ! Initialize mass mixing ratios
    MSO4  = 0e+0_fp
    MNACL = 0e+0_fp
    MH2O  = 0e+0_fp
    MECIL = 0e+0_fp
    MECOB = 0e+0_fp
    MOCIL = 0e+0_fp
    MOCOB = 0e+0_fp
    MDUST = 0e+0_fp

    CALL CHECKMN( 0, 0, 0, Input_Opt, State_Chm, State_Grid, &
                  State_Met, 'AERO_DIADEN called from DEPVEL', RC )

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( MECIL, MECOB, MOCIL, MOCOB, MDUST ) &
    !$OMP PRIVATE( BIN, I, J, TRACID, WID, MH2O, MSO4, MNACL ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
    DO BIN = 1, IBINS

       TRACID = id_NK1 + BIN - 1
       !print *,"TRACID=",TRACID,"id_NK1=",id_NK1, "BIN=", BIN
       WID    = id_NK1 + (ICOMP - 1)*IBINS - 1 + BIN  !(fixed WID to 281-310. dmw 10/3/09)
       !print *, "wid=", WID, "ICOMP=", ICOMP, "IBINS=", IBINS

       ! Get the diameter from an external function
       DIA(I,J,BIN) = GETDP( I, J, LEV, TRACID, State_Met, State_Chm, RC )

       ! Prepare the mass mixing ratio to call external function
       ! for density
       MH2O = Spc(I,J,1,WID)
       !IF ( LSULF30 ) MSO4  = Spc(I,J,LEV,id_SF1-1+BIN)
       IF ( id_SF1 > 0 ) MSO4 = Spc(I,J,LEV,id_SF1-1+BIN)
       !IF ( LSALT30 ) MNACL = Spc(I,J,LEV,id_SS1-1+BIN)
       IF ( id_SS1 > 0 ) MNACL = Spc(I,J,LEV,id_SS1-1+BIN)
       IF ( id_ECIL1 > 0 .AND.id_ECOB1 > 0 .AND. &
            id_OCIL1 > 0 .AND. id_OCOB1 > 0 ) THEN
          !IF ( LCARB30 ) THEN
          MECIL = Spc(I,J,LEV,id_ECIL1-1+BIN)
          MECOB = Spc(I,J,LEV,id_ECOB1-1+BIN)
          MOCIL = Spc(I,J,LEV,id_OCIL1-1+BIN)
          MOCOB = Spc(I,J,LEV,id_OCOB1-1+BIN)
       ENDIF
       !IF ( LDUST30 ) MDUST = Spc(I,J,LEV,id_DUST1-1+BIN)
       IF ( id_DUST1 > 0 ) MDUST = Spc(I,J,LEV,id_DUST1-1+BIN)

       ! Get density from external function
       DENSITY(I,J,BIN) = AERODENS(MSO4,0.e+0_fp,1.875e-1_fp*MSO4, &
                                   MNACL, MECIL, MECOB, &
                                   MOCIL, MOCOB, MDUST, MH2O  )

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units at end of AERO_DIADEN: ' &
             // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine AERO_DIADEN in tomas_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
    ENDIF

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Routine AERO_DIADEN in tomas_mod.F90')
       RETURN
    ENDIF

    ! Free pointer
    Spc => NULL()

    RETURN

  END SUBROUTINE AERO_DIADEN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_mn
!
! !DESCRIPTION: Subroutine CHECKMN use the subroutine MNFIX to check for error
!  in the aerosol mass and number inconsistencies. (win, 7/24/07)
!\\
!\\
! !INTERFACE:
!

  SUBROUTINE CHECKMN( II, JJ, LL, Input_Opt, State_Chm, &
                      State_Grid, State_Met, LOCATION, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD             ! ND60
#endif
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: II, JJ, LL
    CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, L
    INTEGER             :: I1, I2, J1, J2, L1, L2
    INTEGER             :: K, JC, NKID, TRACNUM, MPNUM
    CHARACTER(LEN=255)  :: MSG, LOC ! For species unit check (ewl)
    LOGICAL             :: ERRORSWITCH
                          ! Make ERRORSWITCH = .TRUE. to get full print
                          ! for debugging

    !sfarina
    REAL(fp)            :: Nk(IBINS), Nkd(IBINS)
    REAL(fp)            :: Mk(IBINS, ICOMP)
    REAL(fp)            :: Mkd(IBINS,ICOMP)
    REAL(fp)            :: Gc(ICOMP - 1)
    REAL(fp)            :: Gcd(ICOMP - 1)
    REAL*4              :: BOXVOL
    REAL(fp)            :: XFER(IBINS)

    ! For values from Input_Opt
    LOGICAL             :: prtDebug

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    !=================================================================
    ! CHECKMN begins here!
    !=================================================================

    ! Assume success
    RC        =  GC_SUCCESS

    ! Copy values from Input_Opt
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Point to chemical species array [kg]
    Spc       => State_Chm%Species

    ERRORSWITCH = .FALSE.

    ! Check that species units are in [kg] (ewl, 8/13/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect species units: ' // TRIM(State_Chm%Spc_Units)
       LOC = 'Routine CHECKMN in tomas_mod.F90'
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

    ! Check throughout all grid boxes
    IF ( II == 0 .and. JJ == 0 .and. LL == 0 ) THEN
       I1 = 1
       I2 = State_Grid%NX
       J1 = 1
       J2 = State_Grid%NY
       L1 = 1
       L2 = State_Grid%NZ
    ELSE ! Check at a single grid
       I1 = II
       I2 = II
       J1 = JJ
       J2 = JJ
       L1 = LL
       L2 = LL
    ENDIF

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I, J, L ) &
    !$OMP PRIVATE( Nk, Nkd, Mk, Mkd, K, TRACNUM, JC, MPNUM, BOXVOL ) &
    !$OMP PRIVATE( GC, GCd, ERRORSWITCH ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO I = I1, I2
    DO J = J1, J2
    DO L = L1, L2

       BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3

       ! Swap GEOSCHEM variables into aerosol algorithm variables
       DO K = 1, IBINS
          TRACNUM = id_NK1 - 1 + K
          ! Check for nan
          IF ( IT_IS_NAN( Spc(I,J,L,TRACNUM) ) ) &
               print *, 'Found NaN at',I, J, L,'Species',TRACNUM
          NK(K) = Spc(I,J,L,TRACNUM)
          DO JC = 1, ICOMP-IDIAG
             TRACNUM = id_NK1 - 1 + K + IBINS*JC
             IF ( IT_IS_NAN( Spc(I,J,L,TRACNUM) ) ) &
                  print *, 'Found NaN at',I, J, L,'Species',TRACNUM
             MK(K,JC) = Spc(I,J,L,TRACNUM)
          ENDDO
          MK(K,SRTH2O) = Spc(I,J,L,id_AW1-1+K)
       ENDDO

       DO JC = 1, ICOMP - 1
          Gc(JC) = 0.0   !sfarina - init Gc so storenm doesn't go NaN on us.
       END DO

       ! Get NH4 mass from the bulk mass and scale to bin with sulfate
       IF ( SRTNH4 > 0 ) THEN
          CALL NH4BULKTOBIN( MK(:,SRTSO4), Spc(I,J,L,id_NH4), XFER  )
          MK(1:IBINS,SRTNH4) = XFER(1:IBINS)
          Gc(SRTNH4) = Spc(I,J,L,id_NH3)
       ENDIF

       !if ( i==26 .and. j==57 .and. l==13 ) &
       !     call debugprint(Nk,Mk,i,j,l,'In CHECKMN')

       CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
       !if(i==47.and.j==10.and.l==7) ERRORSWITCH = .TRUE.
       !if(i==22.and.j==33.and.l==2)
       !   call debugprint( NK, MK, i,j,l,LOCATION)

       CALL MNFIX( NK, MK, ERRORSWITCH )
       IF ( ERRORSWITCH ) THEN
          PRINT *, 'CHECKMN is going to terminate at grid',I,J,L
          !IF( .not. SPINUP(14.0) ) THEN
          CALL ERROR_STOP( 'MNFIX found error', LOCATION )
          !ELSE
          !   PRINT *,'Let error go during spin up'
          !ENDIF
       ENDIF

       ! Save the error fixing to diagnostic AERO-FIX
       MPNUM = 5
#ifdef BPCH_DIAG
       IF ( ND60 > 0 ) &
          CALL AERODIAG( MPNUM, I, J, L, Nk, Nkd, Mk, Mkd, BOXVOL, State_Grid )
#endif

       ! Swap Nk and Mk arrays back to Spc
       DO K = 1, IBINS
          TRACNUM = id_NK1 - 1 + K
          Spc(I,J,L,TRACNUM) = Nk(K)
          DO JC = 1, ICOMP-IDIAG
             TRACNUM = id_NK1 - 1 + K + IBINS*JC
             Spc(I,J,L,TRACNUM) = Mk(K,JC)
          ENDDO
          Spc(I,J,L,id_AW1-1+K) = MK(K,SRTH2O)
       ENDDO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( prtDebug ) WRITE(6,*)' #### CHECKMN: finish at ',LOCATION

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHECKMN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mnfix
!
! !DESCRIPTION: Subroutine MNFIX examines the mass and number distrubution and
!  determine if any bins have an average mass outside their normal range.  This
!  can happen because some process, e.g. advection, seems to treat the mass and
!  number species inconsistently.  If any bins are out of range, I shift some
!  mass and number to a new bin in a way that conserves both. (win, 7/23/07)
!  Originally written by Peter Adams, September 2000
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MNFIX ( NK, MK, ERRORSWITCH )
!
! !USES:
!
    USE ERROR_MOD,    ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT) :: NK(IBINS),  MK(IBINS, ICOMP)
    LOGICAL, INTENT(INOUT) :: ERRORSWITCH
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer             :: K,J,KK !counters
    integer             :: NEWBIN !bin number into which mass is shifted
    REAL(fp)            :: XOLD, XNEW !average masses of old and new bins
    REAL(fp)            :: DRYMASS !dry mass of in a bin
    REAL(fp)            :: AVG !average dry mass of particles in bin
    REAL(fp)            :: NUM_INITIAL !number of particles initially in problem bin
    REAL(fp)            :: NSHIFT  !number to shift to new bin
    REAL(fp)            :: MSHIFT !mass to shift to new bin
    REAL(fp)            :: FJ !fraction of mass that is component j
    REAL(fp)            :: save1,save2,save3,save4,save5
    REAL(fp), PARAMETER :: EPS  = 1.e-20_fp !small number for Nk
    REAL(fp), PARAMETER :: EPS2 = 1.e-28_fp !small number for Mk
    REAL(fp), PARAMETER :: TINY = 1.e-36_fp !small number
    REAL(fp), PARAMETER :: VTINY= 1.e-50_fp !very small number

    LOGICAL             :: FIXERROR
    LOGICAL             :: PRT
    REAL(fp)            :: TOTMAS, TOTNUM !for print debug

    !=================================================================
    ! MNFIX begins here!
    !=================================================================

    FIXERROR = .TRUE.
    PRT = .FALSE.
    PRT = ERRORSWITCH !just carrying a signal to print out value at the observed box - since mnfix does not have any information about I,J,L location. (Win, 9/27/05)
    ERRORSWITCH = .FALSE.
    PRT = .FALSE.             !TO AVOID THE HUGE AMOUNT OF PRINTING (JKodros 6/2/15)
    !xk(1)=xk(2)/2.e+0_fp  ! jrp for some reason xk(1) is changing?!
    save1=xk(1)

    ! Check for any incoming negative values or NaN
    !--------------------------------------------------------------------------
    DO K = 1, IBINS
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'11 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'11 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'11 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
       IF ( NK(K) < 0e+0_fp ) THEN
          IF ( PRT ) THEN
             PRINT *,'MNFIX[0]: FOUND NEGATIVE N'
             PRINT *, 'Bin, N', K, NK(K)
          ENDIF
          IF ( ABS(NK(K)) < 1e+0_fp .and. FIXERROR ) THEN
             NK(K) = 0e+0_fp
             IF ( PRT ) PRINT *,'Negative N > -1.0 Reset to zero'
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(0): Found negative Nk(',k,') <-1e+0_fp'
             GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
          ENDIF
       ENDIF
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( MK(K,J) < 0e+0_fp ) THEN
             IF ( PRT ) THEN
                PRINT *,'MNFIX[0]: FOUND NEGATIVE M'
                PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
             ENDIF
             IF( ABS(MK(K,J)) < 1e-5_fp .and. FIXERROR ) THEN
                MK(K,J) = 0e+0_fp
                IF ( PRT ) PRINT *,'Negative M > -1.d-5 Reset to zero'
             ELSE
                ERRORSWITCH =.TRUE.
                print *,'MNFIX(0): Found negative Mk(',k,',comp',j,')'
                GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
             ENDIF
          ENDIF
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO                   !icomp
    ENDDO                     !ibins
    save2=xk(1)

    ! JRP check for neg numbers
    !DO K = 1,IBINS
    !  IF (NK(K) < 0.e+0_fp) THEN
    !     print*,'1 NK < 0 in MNFIX',K,NK(K)
    !  ENDIF
    !  DO J=1,ICOMP
    !     IF (MK(K,J) < 0.e+0_fp) THEN
    !        print*,'1 MK < 0 in MNFIX',K,J,MK(K,J)
    !     ENDIF
    !  ENDDO
    !  IF ( IT_IS_NAN(NK(K)) ) THEN
    !     PRINT *,'11 Found Nan in Nk at bin',K
    !     ERRORSWITCH = .TRUE.
    !     print *,'11 MNFIX(0): Found NaN in Nk(,',k,')'
    !     GOTO 300
    !  ENDIF
    !  DO J = 1, ICOMP
    !     IF ( IT_IS_NAN(MK(K,J)) ) THEN
    !        PRINT *,'11 Found Nan in Mk at bin',K,'component',J
    !        ERRORSWITCH = .TRUE.
    !        GOTO 300
    !     ENDIF
    !  ENDDO
    !ENDDO

    ! Check if both number and mass are zero, if yes then exit mnfix.
    !----------------------------------------------------------------
    TOTNUM = 0e+0_fp
    TOTMAS = 0e+0_fp
    DO K = 1,IBINS
       TOTNUM = TOTNUM + NK(K)
       DO J=1,ICOMP-IDIAG
          TOTMAS = TOTMAS + MK(K,J)
       ENDDO
    ENDDO
    IF ( TOTNUM == 0e+0_fp .AND. TOTMAS == 0e+0_fp ) THEN
       IF ( PRT ) PRINT *,'MNFIX: Nk=Mk=0. Exit now'
       GOTO 300
    ENDIF

    ! If number is tiny ( < EPS) then set it to zero
    !DO K = 1,IBINS
    !   IF ( NK(K) <= EPS ) THEN
    !      NK(K) = 0e+0_fp
    !      DO J= 1, ICOMP-1
    !         MK(K,J) = 0e+0_fp
    !      ENDDO               !STOP  !original (win, 9/1/05)
    !   ENDIF
    !ENDDO

    ! If N is tiny and M is tiny, set both to zeroes
    !--------------------------------------------------------
    DO K = 1, IBINS
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'22 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'22 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'22 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
       IF ( NK(K) <= EPS .AND. NK(K)>= 0e+0_fp ) THEN
          !print*,'1111'
          !print*,k,EPS,xk(K),xk(K+1)
          !print*,'word up'
          NK(K) = EPS
          !NK(K) = 0.e+0_fp
          !DO J = 1, ICOMP-IDIAG
          DO J = 1, ICOMP
             if (J .eq. 1) then
                !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
                MK(K,J) = EPS*AVGMASS(k)
                !MK(K,J) = 0.e+0_fp
             else
                MK(K,J) = VTINY
             endif
          enddo
          !print*,'allbins',MK(:,1)
       ENDIF ! If tiny number
       TOTMAS = SUM(MK(K,1:ICOMP-IDIAG))
       if (TOTMAS.lt.eps2) then
          !print*,'2222'
          NK(K) = EPS
          !NK(K) = 0.e+0_fp
          DO J = 1, ICOMP
             !DO J = 1, ICOMP-IDIAG
             if (J .eq. 1) then
                !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
                MK(K,J) = EPS*AVGMASS(k)
                !MK(K,J) = 0.e+0_fp
             else
                MK(K,J) = VTINY
             endif
          enddo
       endif
    ENDDO
    save3=xk(1)

    ! JRP check for neg numbers
    DO K = 1,IBINS
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'2 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMP
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'2 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'2 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'2 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'2 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    ! Check to see if any bins are completely out of bounds for min or max bin
    !-------------------------------------------------------------------------
    DO K = 1, IBINS
       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMP-IDIAG
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO

       IF ( NK(k) == 0e+0_fp ) THEN
          !AVG = SQRT( xk(K)* xk(K+1) )
          AVG = SQRT( AVGMASS(k) )
       ELSE
          AVG = DRYMASS/ NK(K)
       ENDIF

       IF ( AVG >  xk(IBINS+1) ) THEN
          IF ( PRT ) PRINT *, 'MNFIX [1]: AVG > Xk(ibins+1) at bin',K
          IF ( FIXERROR ) THEN
             !out of bin range - remove some mass
             MSHIFT = NK(k)* xk(IBINS+1)/ 1.2
             DO J= 1, ICOMP
                MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS+EPS2)
             ENDDO
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(1): AVG>Xk(ibins+1) at bin',K
             GOTO 300
          ENDIF
       ENDIF
       IF ( AVG < xk(1)) THEN
          IF( PRT ) PRINT *,'MNFIX [2]: AVG < Xk(1)'
          IF( FIXERROR ) THEN
             !out of bin range - remove some number
             NK(K) = DRYMASS/ ( xk(1)* 1.2 )
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(1): AVG < Xk(1) at bin',K
             GOTO 300
          ENDIF
       ENDIF
    ENDDO

    ! JRP check for neg numbers
    DO K = 1,IBINS
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'3 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMP
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'3 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'3 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'3 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'3 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO
    save4=xk(1)

    !if (PRT) then !<step5.1-temp>
    !   print *,'After_Check2 ---------------------'
    !   do k=1,ibins
    !      totmas = sum(MK(k,1:icomp-1))
    !      print *, totmas,NK(k), totmas/NK(k)
    !   enddo
    !endif

    !print*,1,NK(1),NK(2)
    !print*,1,MK(1,:)
    !print*,1,MK(2,:)

    ! Check to see if any bins are out of bounds
    !-------------------------------------------------------------------
    DO K = 1, IBINS
       !if (PRT) print *,'Now at bin',k !<step4.4>tmp (win, 9/28/05)

       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMP-IDIAG
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO

       IF ( NK(K) == 0e+0_fp ) THEN
          !AVG = SQRT(xk(K)*xk(K+1)) !set to mid-range value
          AVG = AVGMASS(k) !set to mid-range value
       ELSE
          AVG = DRYMASS/NK(K)
       ENDIF

       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check3---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif

       ! If over boundary of the current bin
       IF ( AVG >  xk(K+1) ) THEN
          IF ( PRT ) PRINT *, 'MNFIX [3]: AVG>Xk(',K+1,')'
          IF ( PRT ) CALL DEBUGPRINT(NK,MK,0,0,0,'inside MNFIX')
          IF ( FIXERROR ) THEN
             !Average mass is too high - shift to higher bin
             !KK = K + 1 ! jrp, this was causing errors
             !ERRORSWITCH=.TRUE.
             KK = K
             XNEW = xk(KK+1)/ 1.1
             if ( PRT ) PRINT *, 'k',k,'AVG',AVG,' XNEW ',XNEW
100          IF ( XNEW <= AVG ) THEN
                IF ( KK < IBINS ) THEN
                   KK = KK + 1
                   XNEW = xk(KK+1)/ 1.1
                   if (PRT) PRINT *, '..move up to bin ',KK,' XNEW ',XNEW
                   GOTO 100
                ELSE
                   ! Already reach highest bin - must remove some mass (win, 8/1/07)
                   ! Updated by jrp 3/1/2012
                   MSHIFT = NK(k)* xk(k+1)/ 1.1
                   if( PRT ) PRINT*,' Mass being discarded: '
                   DO J= 1, ICOMP
                      !if (PRT)
                      !print*,'Removing mass in MNFIX',MSHIFT, DRYMASS
                      MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS)
                   ENDDO
                   ! and recalculate dry mass (win, 8/1/07)
                   DRYMASS = 0.e+0_fp ! jrp fix 2/29/12
                   DO J = 1, ICOMP-IDIAG
                      DRYMASS = DRYMASS + MK(K,J)
                   ENDDO
                   GOTO 111
                ENDIF
             ENDIF

             if(PRT)print*,'Old NK',NK(k),'Old DRYMASS',DRYMASS,'bin',k

             !XOLD = SQRT( xk(K)* xk(K+1) )
             XOLD = AVGMASS(k)
             NUM_INITIAL = NK(K)
             NSHIFT = ( DRYMASS - XOLD * NUM_INITIAL )/ ( XNEW - XOLD )
             MSHIFT = XNEW * NSHIFT
             NK(K) = NK(K) - NSHIFT
             NK(KK) =NK(KK) + NSHIFT

             if(prt) then
                print*,'NSHIFT',NSHIFT, 'MSHIFT',MSHIFT
                print*,'New NK',k,NK(k),' Nk(kk)',kk,NK(kk)
                print*,'Total mass bin',k,sum(MK(k,1:icomp-idiag))
                print*,'SO4 mass bin  ',k,(MK(k,srtso4))
                print*,'Total mass bin',kk,sum(MK(kk,1:icomp-idiag))
                print*,'SO4 mass bin  ' ,kk,(MK(kk,srtso4))
             endif

             DO J = 1, ICOMP-IDIAG
                FJ = MK(K,J)/ DRYMASS
                MK(K,J) = XOLD * NK(K) * FJ
                MK(KK,J) = MK(KK,J) + MSHIFT * FJ
             ENDDO

             if(prt) then
                print*,'After shift mass'
                print*,'Total mass bin',k,sum(MK(k,1:icomp-idiag))
                print*,'SO4 mass bin  ',k,(MK(k,srtso4))
                print*,'Total mass bin',kk,sum(MK(kk,1:icomp-idiag))
                print*,'SO4 mass bin  ',kk,(MK(kk,srtso4))
             endif

          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(3) : AVG>Xk(',K+1,')'
             GOTO 300
          ENDIF    ! Fixerror
       ENDIF       ! AVG > Xk(k+1)

       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check4---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif

       ! If under boundary of the current bin
111    IF ( AVG <  xk(K) ) THEN
          IF ( PRT ) PRINT *,'MNFIX [4]: AVG<Xk(',K,')'
          IF ( FIXERROR ) THEN
             !average mass is too low - shift number to lower bin
             !KK = K - 1 ! jrp potential for errors here
             KK = K
             XNEW = xk(KK)* 1.1
200          IF ( XNEW >= AVG ) THEN
                IF ( KK > 1 ) THEN
                   KK = KK - 1
                   XNEW = xk(KK)* 1.1
                   GOTO 200
                ELSE
                   ! Already reach lowest bin - must remove some number (win, 8/1/07)
                   NK(K) = DRYMASS/ ( xk(1)* 1.2 )
                   GOTO 222
                ENDIF
             ENDIF
             !XOLD = SQRT(xk(K)* xk(K+1))
             XOLD = AVGMASS(k)
             NUM_INITIAL = NK(K)
             NSHIFT = NUM_INITIAL - DRYMASS/XOLD !(win, 10/20/08)
             !Prior to 10/20/08 (win)
             !NSHIFT = (DRYMASS - XOLD * NUMBER)/ ( XNEW - XOLD )
             MSHIFT = XNEW * NSHIFT
             NK(K) = NK(K) - NSHIFT
             NK(KK) = NK(KK) + NSHIFT
             DO J=1,ICOMP
                FJ = MK(K,J)/ DRYMASS
                MK(K,J) = XOLD * NK(K) * FJ
                MK(KK,J) = MK(KK,J) + MSHIFT * FJ
             ENDDO
             
          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(4): AVG < Xk(',k,')'
             GOTO 300
          ENDIF
222    ENDIF
       
       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check5---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif
       !if (PRT) print *,MK(k,1),NK(k), MK(k,1)/NK(k),'Check5'!<step4.4>tmp (win, 9/28/05)

    ENDDO ! loop bin
    save5=xk(1)

    !print*,2,NK(1),NK(2)
    !print*,2,MK(1,:)
    !print*,2,MK(2,:)

    ! JRP check for neg numbers
    DO K = 1,IBINS
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'4 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMP
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'4 MK < 0 in MNFIX',K,J,MK(K,J)
             print*,'saved xk1s',save1,save2,save3,save4,save5
             print*,'xk',xk
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'4 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'4 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'4 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    !if (PRT) then !<step5.1-temp>
    ! Catch any small negative values resulting from fixing
    !--------------------------------------------------------------------------
    DO K = 1, IBINS
       IF ( NK(K) < 0e+0_fp ) THEN
          IF ( PRT ) THEN
             PRINT *,'MNFIX[5]: FOUND NEGATIVE N'
             PRINT *, 'Bin, N', K, NK(K)
          ENDIF
          IF ( ABS(NK(K)) < 1e+0_fp .and. FIXERROR ) THEN
             NK(K) = 0e+0_fp
             IF ( PRT ) PRINT *,'Negative N > -1.0 Reset to zero'
          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(5): Negative N after fixing at bin',k
             GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
          ENDIF
       ENDIF
       DO J = 1, ICOMP
          IF ( MK(K,J) < 0e+0_fp ) THEN
             IF ( PRT ) THEN
                PRINT *,'MNFIX[6]: FOUND NEGATIVE M'
                PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
             ENDIF
             IF( ABS(MK(K,J)) < 1D-5 .and. FIXERROR ) THEN
                MK(K,J) = 0e+0_fp
                IF ( PRT ) PRINT *,'Negative M > -1.d-5 Reset to zero'
             ELSE
                ERRORSWITCH =.TRUE.
                PRINT *, 'MNFIX(6): Negative M after fixing at bin',k
                GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
             ENDIF
          ENDIF
       ENDDO                   !icomp
    ENDDO                     !ibins

    ! JRP check for neg numbers
    DO K = 1,IBINS
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'5 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMP
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'5 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'5 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'5 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMP
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'5 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    ! Check any last inconsistent M=0 or N=0
    !--------------------------------------------------------
    DO K = 1, IBINS
       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMP-IDIAG
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO
       IF ( NK(K) /= 0e+0_fp .AND. DRYMASS == 0e+0_fp .or. &
            NK(K) == 0e+0_fp .AND. DRYMASS /= 0e+0_fp     ) THEN
          PRINT *, '5.5 set nk, mk to ZERO for all bins'
          DO J = 1, ICOMP
             MK(K,J)=0.e+0_fp
             NK(K) = 0.e+0_fp
          ENDDO
          MK(K,ICOMP) = 0.e+0_fp !Set aerosol water to zero too
       ENDIF                  ! If tiny number
    ENDDO

    ! JRP check for neg numbers
    DO K = 1,IBINS
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'6 NK < 0 in MNFIX',K,NK(K)
          STOP
       ENDIF
       DO J=1,ICOMP
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'6 MK < 0 in MNFIX',K,J,MK(K,J)
             STOP
          ENDIF
       ENDDO
    ENDDO

300 CONTINUE

    IF (ERRORSWITCH) THEN
555    FORMAT (3E15.5E2)
       WRITE(6,*)'END OF MNFIX ( WHERE? )'
       WRITE(6,*)'DRYMAS-excl-NH4  NK      DRYMASS/NK'
       DO K = 1,IBINS
          TOTMAS = SUM(MK(K,1:ICOMP-1))
          !PRINT *, TOTMAS,NK(K), TOTMAS/NK(K)
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% NOTE: NK will be IBINS+1 upon exiting the loop, which will cause an
       !%%% out-of-bounds error.  Comment this out for now, unless it should be
       !%%% inserted into the DOloop
       !WRITE(6,555)
       !        TOTMAS, NK(K),
       !        TOTMAS/ NK(K)
       !print*,'-----------'
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       call debugprint( NK, MK, 0,0,0,'End of MNFIX')

       !write(*,*)'Nk'
       !write(*,*) NK(1:30)
       !write(*,*)'Mk(srtso4)'
       !write(*,*) MK(1:30,srtso4)
       !write(*,*)'Mk(srth2o)'
       !write(*,*) MK(1:30,srth2o)
       !STOP 'Negative Nk or Mk at after mnfix'  !comment out this to make it stop outside mnfix so that I can print out the i,j,l (location) of the error (win, 9/1/05)
    ENDIF

    RETURN

  END SUBROUTINE MNFIX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgridcoag
!
! !DESCRIPTION: Subroutine SUBGRIDCOAG determine how much of each size of
!  freshly emitted aerosol will be scavenged by coagulation prior to being
!  completely mixed in the gridbox and will give the new emissions size
!  distribution along with where the mass of coagulated particles should be
!  added.
!  Written by Jeff Pierce, December, 2006
!  Implement in GEOS-Chem by Win Trivitayanurak, 10/4/07
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SUBGRIDCOAG( NDISTINIT, NDIST, MDIST, BOXVOLUME,TEMP, &
                          PRES,TSCALE, NDISTFINAL, MADDFINAL, pdbug)
!
! !ARGUMENTS:
!
    !ndistinit(nbins)   : the number of particles being added to the gridbox
    !                     before subgrid coag
    !ndist(nbins)       : the number of particles in the box
    !mdist(nbins,icomp) : the mass of each component in the box. (kg)
    !boxvolume          : volume of box in cm3
    !tscale             : the scale time for mixing (s)
    !ndistfinal(nbins)  : the number of particles being added to the gridbox
    !                     after subgrid coag
    !maddfinal(nbins)   : the mass that should be added to each bin due to
    !                     coagulation (kg)
    REAL(fp) ndistinit(ibins)
    REAL(fp) ndist(ibins)
    REAL(fp) mdist(ibins,icomp)
    REAL*4 boxvolume, temp , PRES
    REAL*4 tscale
    REAL(fp) ndistfinal(ibins)
    REAL(fp) maddfinal(ibins)
    logical pdbug  ! for pringing out during debugging
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) mp                     ! mass of the particle (kg)
    REAL*4 density                !density (kg/m3) of particles
    REAL(fp) fracdiaml(ibins,ibins) ! fraction of coagulation that occurs with each bin larger
    REAL(fp) kcoag(ibins) ! the coagulation rate for the particles in each bin (s^-1)
    !REAL*4 aerodens
    !external aerodens
    
    REAL*4 mu                     !viscosity of air (kg/m s)
    REAL*4 mfp                    !mean free path of air molecule (m)
    REAL*4 Kn                     !Knudsen number of particle
    REAL*4 beta                   !correction for coagulation coeff.
    REAL(fp) Mktot      !total mass of aerosol
    REAL*4 kij(ibins,ibins)
    REAL*4 Dpk(ibins)             !diameter (m) of particles in bin k
    REAL*4 Dk(ibins)              !Diffusivity (m2/s) of bin k particles
    REAL*4 ck(ibins)              !Mean velocity (m/2) of bin k particles
    REAL*4 neps
    REAL*4 meps
    INTEGER I, J, K, KK
    LOGICAL ERRORSWITCH

    ! Adjustable parameters
    real*4 pi, kB               !kB is Boltzmann constant (J/K)
    real*4 R       !gas constant (J/ mol K)
    parameter (pi=3.141592654, kB=1.38e-23, R=8.314)
    parameter (neps=1E8, meps=1E-8)

    !=================================================================
    ! SUBGRIDCOAG begins here!
    !=================================================================

    if (pdbug) call debugprint(Ndist,Mdist,0,0,0,'SUBDGRIDCOAG: Entering')

    !Before going in to calculation, check and fix Nk and Mk
    ERRORSWITCH = .FALSE.
    !print *, 'mnfix in tomas_mod:7495'
    CALL MNFIX( NDIST, MDIST, ERRORSWITCH )
    IF ( ERRORSWITCH ) THEN
       PRINT *,'SUBGRIDCOAG: MNFIX found error: Entering SUBGRIDCOAG'
       PDBUG = .TRUE.
       GOTO 11
    ENDIF

    if (pdbug) call debugprint(Ndist,Mdist,0,0,0,'SUBDGRIDCOAG: after MNFIX_1')

    mu=2.5277e-7*temp**0.75302
    mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*temp) ) )  !S&P eqn 9.6
    ! Calculate particle sizes and diffusivities
    do k=1,ibins
       Mktot=0.2*mdist(k,srtso4)
       do j=1, icomp
          Mktot=Mktot+mdist(k,j)
       enddo
       if (Mktot.gt.meps)then
          density=aerodens(mdist(k,srtso4),0e+0_fp, &
                  0.1875e+0_fp*mdist(k,srtso4),mdist(k,srtnacl), &
                  mdist(k,srtecil),mdist(k,srtecob), &
                  mdist(k,srtocil),mdist(k,srtocob),mdist(k,srtdust), &
                  mdist(k,srth2o)) !assume bisulfate
       else
          density = 1400.
       endif
       if(ndist(k).gt.neps .and. Mktot.gt.meps)then
          mp=Mktot/ndist(k)
       else
          mp=sqrt(xk(k)*xk(k+1))
       endif
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
       Dk(k)=kB*temp/(3.0*pi*mu*Dpk(k)) &           !S&P Table 12.1
             *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
       ck(k)=sqrt(8.0*kB*temp/(pi*mp))              !S&P Table 12.1
    enddo

    ! Calculate coagulation coefficients

    do i=1,ibins
       do j=1,ibins
          Kn=4.0*(Dk(i)+Dk(j)) &
             /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
          beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
          !This is S&P eqn 12.46 with non-continuum correction, beta
          kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
          kij(i,j)=kij(i,j)*1.0e6/boxvolume  !normalize by grid cell volume
       enddo
    enddo

    ! get the first order loss rate
    kcoag(ibins)=0.0
    !debug
    if(pdbug) print *,'Bin  KCOAG'
    do k=1,ibins-1
       kcoag(k)=0.0
       do kk=k+1,ibins
          kcoag(k)=kcoag(k)+kij(k,kk)*ndist(kk)
       enddo
       !debug
       if(pdbug) print *, k, kcoag(k)
    enddo

    ! get the fraction of the coagulation that occurs from each bin larger
    do k=1,ibins
       do kk=1,ibins
          fracdiaml(k,kk)=0.
       enddo
    enddo
    do k=1,ibins-1
       !debug
       if(pdbug) print *, 'Bin k', k
       !debug
       if(pdbug) print *, 'Bin kk   fracdiaml(k,kk)'
       do kk=k+1,ibins
          if (kcoag(k).gt.0.e+0_fp)then
             fracdiaml(k,kk)=kij(k,kk)*ndist(kk)/kcoag(k)
          else
             fracdiaml(k,kk)=0
          endif
          !debug
          if(pdbug) print *, kk, fracdiaml(k,kk)
       enddo
    enddo

    ! determine the number of new particles left after coagulation
    do k=1,ibins
       ndistfinal(k)=ndistinit(k)*exp(-kcoag(k)*tscale)
    enddo

    ! determine the mass added to each bin coagulation
    do k=1,ibins
       maddfinal(k)=0.
    enddo
    do k=1,ibins-1
       do kk=k+1,ibins
          maddfinal(kk)=maddfinal(kk) + (ndistinit(k)-ndistfinal(k))* &
                        fracdiaml(k,kk)*AVGMASS(k)
          !sfarina - not the slowest part of subgridcoag, but every little bit helps
          !maddfinal(kk)=maddfinal(kk) + (ndistinit(k)-ndistfinal(k))* &
          !     fracdiaml(k,kk)*sqrt(xk(k)*xk(k+1))
       enddo
    enddo

    pdbug = .FALSE.

11  continue
    return

  END SUBROUTINE SUBGRIDCOAG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: storenm
!
! !DESCRIPTION: Subroutine STORENM stores values of Nk and Mk into Nkd and Mkd
!  for diagnostic purposes.  Also do gas phase concentrations. (win, 7/23/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd )
!
! !INPUT PARAMETERS:
!
    REAL(fp),INTENT(IN)    :: Nk(IBINS)
    REAL(fp),INTENT(IN)    :: Mk(IBINS, ICOMP)
    REAL(fp),INTENT(IN)    :: Gc(ICOMP - 1)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),INTENT(OUT)   :: Nkd(IBINS)
    REAL(fp),INTENT(OUT)   :: Mkd(IBINS,ICOMP)
    REAL(fp),INTENT(OUT)   :: Gcd(ICOMP - 1)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: J, K

    !sfarina
    DO J= 1, ICOMP-1
       Gcd(J)=Gc(J)
    ENDDO
    DO K = 1, IBINS
       Nkd(K)=Nk(K)
       DO J= 1, ICOMP
          Mkd(K,J)=Mk(K,J)
       ENDDO
    ENDDO

    RETURN

  END SUBROUTINE STORENM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: debugprint
!
! !DESCRIPTION: Subroutine DEBUGPRINT print out the Nk and Mk values for error
!  checking (win, 9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DEBUGPRINT( NK, MK, I,J,L, LOCATION)
!
! !INPUT PARAMETERS:
!
    REAL(fp),         INTENT(IN) :: Nk(IBINS), MK(IBINS,ICOMP)
    INTEGER,          INTENT(IN) :: I,J,L
    CHARACTER(LEN=*), INTENT(IN) :: LOCATION
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: JC, k, B

    WRITE(*,*) LOCATION, I, J, L
    !write(6,*) 'Nk(1:30)'
    !write(6,*) Nk(1:30)
    !do jc=1,icomp
    !   write(6,*) 'Mk(1:30) comp:',jc
    !   write(6,*) Mk(1:30,jc)
    !enddo
    write(*,111) 'Bin        Num         SO4        NaCl        ', &
                 'ECIL        ECOB        OCIL       OCOB       ', &
                 'Dust        NH4         Water  '
    DO K = 1, IBINS
       write(*,110) k,Nk(k), ( Mk(K,B), B=1,ICOMP )
    ENDDO
    write(*,*) ' '
110 FORMAT ( I2, 10E12.5 )
111 FORMAT (a,a,a)

  END SUBROUTINE DEBUGPRINT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nh4bulktobin
!
! !DESCRIPTION: Subroutine NH4BULKTOBIN takes the bulk ammonium aerosol from
!  GEOS-Chem and fraction it to each bin according to sulfate mole fraction in
!  each bin
!  Written by Win Trivityanurak, Sep 26, 2008
!  .
!  Make sure that we work with mass or mass conc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NH4BULKTOBIN( MSULF, NH4B, MAMMO )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)   :: MSULF(IBINS)  ! size-resolved sulfate [kg]
    REAL(fp),  INTENT(IN)   :: NH4B          ! Bulk NH4 mass [kg]
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT)  :: MAMMO(IBINS)  ! size-resolved NH4 [kg]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: K
    REAL(fp)                :: TOTMASS, NH4TEMP

    !=================================================================
    ! NH4BULKTOBIN begins here
    !=================================================================

    MAMMO(:) = 0.e+0_fp

    ! Sum the total sulfate
    TOTMASS = 0.e+0_fp
    DO K = 1, IBINS
       TOTMASS = TOTMASS + MSULF(K)
    ENDDO

    IF ( TOTMASS .eq. 0.e+0_fp ) RETURN

    ! Limit the amount of NH4 entering TOMAS calculation
    ! if it is very NH4-rich, just limit the amount to balance
    ! existing 30-bin-summed SO4 assuming (NH4)2SO4 in such case
    !  (NH4)2 mass = (SO4)mass / 96. * 2. * 18. = 0.375*(SO4)mass
    ! (win, 9/28/08)
    NH4TEMP = NH4B
    IF ( NH4B/TOTMASS > 0.375e+0_fp ) &  !make sure we use mass ratio
         NH4TEMP = 0.375e+0_fp * TOTMASS

    ! Calculate ammonium aerosol scale to each bin
    DO K = 1, IBINS
       MAMMO(K) = MSULF(K) / TOTMASS * NH4TEMP
    ENDDO

    !write(777,*) NH4B/TOTMASS

    RETURN

  END SUBROUTINE NH4BULKTOBIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodens
!
! !DESCRIPTION: Function AERODENS calculates the density (kg/m3) of a sulfate-
!  nitrate-ammonium-nacl-OC-EC-dust-water mixture.  Inorganic mass (sulfate-
!  nitrate-ammonium-nacl-water) is assumed to be internally mixed.  Then the
!  density of inorg and EC, OC, and dust is combined weighted by mass.
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II' and extened to include
!  carbonaceous aerosol in Jan, 2002.
!\\
!\\
! !INTERFACE:
!
  FUNCTION AERODENS( MSO4, MNO3, MNH4, MNACL, MECIL, MECOB, MOCIL, &
                     MOCOB, MDUST, MH2O )  RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)  ::  MSO4, MNO3, MNH4, MNACL, MH2O
    REAL(fp),  INTENT(IN)  ::  MECIL, MECOB, MOCIL, MOCOB, MDUST
!
! !RETURN VALUE:
!
    REAL(fp)                  :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)                  :: IDENSITY, DEC, DOC, DDUST, MTOT
    parameter(dec=2200., doc=1400., ddust=2650.)

    !=================================================================
    ! AERODENS begins here!
    !=================================================================

    IDENSITY = INODENS( MSO4, MNO3, MNH4, MNACL, MH2O )
    MTOT = MSO4+MNO3+MNH4+MNACL+MH2O+MECIL+MECOB+MOCIL+MDUST+MOCOB
    IF ( MTOT > 0.e+0_fp ) THEN
       VALUE = ( IDENSITY*(MSO4+MNO3+MNH4+MNACL+MH2O) + &
                 DEC*(MECIL+MECOB) + DOC*(MOCIL+MOCOB)+ &
                 DDUST*MDUST                            )/MTOT
    ELSE
       VALUE = 1400.
    ENDIF

  END FUNCTION AERODENS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: inodens
!
! !DESCRIPTION: Function INODENS calculates the density (kg/m3) of a sulfate-
!  nitrate-ammonium-nacl-water mixture that is assumed to be internally mixed.
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II'
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) 8/6/07 first
!  as AERODENS, then change to INODENS on 9/3/07
!\\
!\\
! !INTERFACE:
!
  FUNCTION INODENS( MSO4_, MNO3_, MNH4_, MNACL_, MH2O_ ) &
       RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    ! mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
    ! component.  Since the density is an intensive property,
    ! these may be input in a variety of units (ug/m3, mass/cell, etc.).
    REAL(fp),  INTENT(IN)  ::  MSO4_, MNO3_, MNH4_, MNACL_, MH2O_
!
! !RETURN VALUE:
!
    REAL(fp)               :: VALUE
!
! !REMARKS:
! ----Literature cited----
!     I. N. Tang and H. R. Munkelwitz, Water activities, densities, and
!       refractive indices of aqueous sulfates and sodium nitrate droplets
!       of atmospheric importance, JGR, 99, 18,801-18,808, 1994
!     Ignatius N. Tang, Chemical and size effects of hygroscopic aerosols
!       on light scattering coefficients, JGR, 101, 19,245-19,250, 1996
!     Ignatius N. Tang, Thermodynamic and optical properties of mixed-salt
!       aerosols of atmospheric importance, JGR, 102, 1883-1893, 1997
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp) MSO4, MNO3, MNH4, MNACL, MH2O
    !real(fp) so4temp, no3temp, nh4temp, nacltemp, h2otemp
    real(fp) mwso4, mwno3, mwnh4, mwnacl, mwh2o            !molecular weights
    real(fp) ntot, mtot                          !total number of moles, mass
    real(fp) nso4, nno3, nnh4, nnacl, nh2o       !moles of each species
    real(fp) xso4, xno3, xnh4, xnacl, xh2o       !mole fractions
    real(fp) rso4, rno3, rnh4, rnacl, rh2o       !partial molar refractions
    real(fp) ran, rs0, rs1, rs15, rs2       !same, but for solute species
    real(fp) asr                            !ammonium/sulfate molar ratio
    real(fp) nan, ns0, ns1, ns15, ns2, nss  !moles of dry solutes (nss = sea salt)
    real(fp) xan, xs0, xs1, xs15, xs2, xss  !mass % of dry solutes - Tang (1997) eq. 10
    real(fp) dan, ds0, ds1, ds15, ds2, dss  !binary solution densities - Tang (1997) eq. 10
    real(fp) mwan, mws0, mws1, mws15, mws2  !molecular weights
    real(fp) yan, ys0, ys1, ys15, ys2, yss  !mole fractions of dry solutes
    real(fp) yh2o
    real(fp) d                              !mixture density
    real(fp) xtot

    ! In the lines above, "an" refers to ammonium nitrate, "s0" to
    ! sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
    ! "nacl" or "ss" is sea salt.
    parameter(mwso4=96.e+0_fp, &
              mwno3=62.e+0_fp, &
              mwnh4=18.e+0_fp, &
              mwh2o=18.e+0_fp, &
              mwnacl=58.45e+0_fp)
    parameter(mwan=mwnh4+mwno3,          &
              mws0=mwso4+2.e+0_fp,       &
              mws1=mwso4+1.e+0_fp+mwnh4, &
              mws2=2.e+0_fp*mwnh4+mwso4)

    !=================================================================
    ! INODENS begins here!
    !=================================================================

    ! Pass initial component masses to local variables
    mso4=mso4_
    mno3=mno3_
    mnh4=mnh4_
    mnacl=mnacl_
    mh2o=mh2o_

    !so4temp=mso4
    !no3temp=mno3
    !nh4temp=mnh4
    !h2otemp=mh2o
    !nacltemp=mnacl

    ! [Pengfei Liu, avoid equality test with floating-point real numbers
    !<step4.7> if the aerosol mass is zero - then just return the
    !typical density = 1500 kg/m3 (win, 1/4/06)
    !if (mso4 .eq. 0.e+0_fp .and. mno3 .eq.0.e+0_fp &
    !    .and. mnh4.eq.0.e+0_fp .and. mnacl .eq. 0.e+0_fp ) then
    !   VALUE = 1500.e+0_fp !kg/m3
    !   goto 10
    !endif
    if ((mso4+mno3+mnh4+mnacl) .gt. 0.e+0_fp) then
       CONTINUE
    else
       VALUE = 1500.e+0_fp !kg/m3
       RETURN
    endif
    ! Pengfei Liu, 2018/02/07]

    ! Calculate mole fractions
    mtot  = mso4+mno3+mnh4+mnacl+mh2o
    nso4  = mso4/mwso4
    nno3  = mno3/mwno3
    nnh4  = mnh4/mwnh4
    nnacl = mnacl/mwnacl
    nh2o  = mh2o/mwh2o
    ntot  = nso4+nno3+nnh4+nnacl+nh2o
    xso4  = nso4/ntot
    xno3  = nno3/ntot
    xnh4  = nnh4/ntot
    xnacl = nnacl/ntot
    xh2o  = nh2o/ntot

    ! If there are more moles of nitrate than ammonium, treat unneutralized
    ! HNO3 as H2SO4
    if (nno3 .gt. nnh4) then
       !make the switch
       nso4=nso4+(nno3-nnh4)
       nno3=nnh4
       mso4=nso4*mwso4
       mno3=nno3*mwno3

       !recalculate quantities
       mtot = mso4+mno3+mnh4+mnacl+mh2o
       nso4 = mso4/mwso4
       nno3 = mno3/mwno3
       nnh4 = mnh4/mwnh4
       nnacl = mnacl/mwnacl
       nh2o = mh2o/mwh2o
       ntot = nso4+nno3+nnh4+nnacl+nh2o
       xso4 = nso4/ntot
       xno3 = nno3/ntot
       xnh4 = nnh4/ntot
       xnacl = nnacl/ntot
       xh2o = nh2o/ntot

    endif

    ! Calculate the mixture density
    ! Assume that nitrate exists as ammonium nitrate and that other ammonium
    ! contributes to neutralizing sulfate
    nan=nno3
    if (nnh4 .gt. nno3) then
       !extra ammonium
       asr=(nnh4-nno3)/nso4
    else
       !less ammonium than nitrate - all sulfate is sulfuric acid
       asr=0.e+0_fp
    endif
    if (asr .ge. 2.e+0_fp) asr=2.e+0_fp
    if (asr .ge. 1.e+0_fp) then
       !assume NH4HSO4 and (NH4)2(SO4) mixture
       !NH4HSO4
       ns1=nso4*(2.e+0_fp-asr)
       !(NH4)2SO4
       ns2=nso4*(asr-1.e+0_fp)
       ns0=0.e+0_fp
    else
       !assume H2SO4 and NH4HSO4 mixture
       !NH4HSO4
       ns1=nso4*asr
       !H2SO4
       ns0=nso4*(1.e+0_fp-asr)
       ns2=0.e+0_fp
    endif

    !Calculate weight percent of solutes
    xan=nan*mwan/mtot*100.e+0_fp
    xs0=ns0*mws0/mtot*100.e+0_fp
    xs1=ns1*mws1/mtot*100.e+0_fp
    xs2=ns2*mws2/mtot*100.e+0_fp
    xnacl=nnacl*mwnacl/mtot*100.e+0_fp
    xtot=xan+xs0+xs1+xs2+xnacl

    ! [Pengfei Liu, fix the polynomial issue
    !Calculate binary mixture densities (Tang, eqn 9)
    !dan=0.9971e+0_fp +4.05e-3_fp*xtot +9.0e-6_fp*xtot**2.e+0_fp
    !ds0=0.9971e+0_fp +7.367e-3_fp*xtot -4.934d-5*xtot**2.e+0_fp &
    !     +1.754e-6_fp*xtot**3.e+0_fp - 1.104d-8*xtot**4.e+0_fp
    !ds1=0.9971e+0_fp +5.87e-3_fp*xtot -1.89e-6_fp*xtot**2.e+0_fp &
    !     +1.763e-7_fp*xtot**3.e+0_fp
    !ds2=0.9971e+0_fp +5.92e-3_fp*xtot -5.036e-6_fp*xtot**2.e+0_fp &
    !     +1.024d-8*xtot**3.e+0_fp
    !dss=0.9971e+0_fp +7.41e-3_fp*xtot -3.741d-5*xtot**2.e+0_fp &
    !     +2.252e-6_fp*xtot**3.e+0_fp   -2.06d-8*xtot**4.e+0_fp
    dan=0.9971e+0_fp + xtot * (4.05e-3_fp + 9.0e-6_fp * xtot)
    ds0=0.9971e+0_fp &
        +xtot*(7.367e-3_fp &
        +xtot*(-4.934d-5 &
        +xtot*(1.754e-6_fp &
        +xtot*(-1.104d-8  ))))
    ds1=0.9971e+0_fp &
        +xtot*(5.87e-3_fp &
        +xtot*(-1.89e-6_fp &
        +xtot*(1.763e-7_fp )))
    ds2=0.9971e+0_fp &
        +xtot*(5.92e-3_fp &
        +xtot*(-5.036e-6_fp &
        +xtot*(1.024d-8    )))
    dss=0.9971e+0_fp &
        +xtot*(7.41e-3_fp &
        +xtot*(-3.741d-5 &
        +xtot*(2.252e-6_fp &
        +xtot*(-2.06d-8     ))))
    ! Pengfei Liu, 2018/02/07]

    !Convert x's (weight percent of solutes) to fraction of dry solute (scale to 1)
    xtot=xan+xs0+xs1+xs2+xnacl
    xan=xan/xtot
    xs0=xs0/xtot
    xs1=xs1/xtot
    xs2=xs2/xtot
    xnacl=xnacl/xtot

    !Calculate mixture density
    d=1.e+0_fp/(xan/dan+xs0/ds0+xs1/ds1+xs2/ds2+xnacl/dss)  !Tang, eq. 10

    if ((d .gt. 2.e+0_fp) .or. (d .lt. 0.997e+0_fp)) then
       write(*,*) 'ERROR in aerodens'
       write(*,*) mso4,mno3,mnh4,mnacl,mh2o
       print *, 'xtot',xtot
       print *, 'xs1',xs1, 'ns1',ns1,'mtot',mtot,'asr',asr
       write(*,*) 'density(g/cm3)',d
       STOP
    endif

    ! Restore masses passed
    !mso4=so4temp
    !mno3=no3temp
    !mnh4=nh4temp
    !mnacl=nacltemp
    !mh2o=h2otemp

    ! Return the density
    VALUE = 1000.e+0_fp*d    !Convert g/cm3 to kg/m3

    !<step4.7> negative value check (win, 1/4/06)
    if ( VALUE < 0e+0_fp ) then
       print *, 'ERROR :: aerodens - negative', VALUE
       STOP
    endif

10  CONTINUE

  END FUNCTION INODENS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dmdt_int
!
! !DESCRIPTION: Function DMDT_INT apply the analytic solution to the droplet
!  growth equation in mass space for a given scale length which mimics the
!  inclusion of gas kinetic effects. (win, 7/23/07)
!  Originally written by Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!\\
!\\
! !INTERFACE:
!
  FUNCTION DMDT_INT ( M0, TAU, WR ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    ! M0  initial mass
    ! L0  length scale
    ! Tau forcing from vapor field
    REAL(fp),   INTENT(IN)  ::  M0,  TAU,  WR
!
! !RETURN VALUE:
!
    REAL(fp)                :: VALUE
!
! !REMARKS:
!  Original note from Peter Adams:
!  I have changed the length scale.  Non-continuum effects are
!  assumed to be taken into account in choice of tau (in so4cond subroutine).
!  .
!  I have also added another argument to the function call, WR.  This
!  is the ratio of wet mass to dry mass of the particle.  I use this
!  information to calculate the amount of growth of the wet particle,
!  but then return the resulting dry mass.  This is the appropriate
!  way to implement the condensation algorithm in a moving sectional
!  framework.
!  .
!  Reference: Stevens et al. 1996, Elements of the Microphysical Structure
!           of Numerically Simulated Nonprecipitating Stratocumulus,
!           J. Atmos. Sci., 53(7),980-1006.
! This calculates a solution for m(t+dt) using eqn.(A3) from the reference
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)                ::  X,  L0,  C,  ZERO,  MH2O
    PARAMETER (C=2.e+0_fp/3.e+0_fp,L0=0.0e+0_fp,ZERO=0.0e+0_fp)

    !=================================================================
    ! DMDT_INT begins here!
    !=================================================================

    MH2O = ( WR - 1.e+0_fp ) * M0
    X = ( ( M0 + MH2O ) ** C + L0 )
    X = MAX( ZERO, SQRT(MAX(ZERO,C*TAU+X))-L0 )

    !<step5.3> Do aqueous oxidation dry - so no need to select process (win, 7/14/06)
    !<step5.3> For so4cond condensation, use constant water amount.
    ! For aqueous oxidation, use constant wet ratio. (win, 7/13/06)
    !prior to 10/2/08
    !VALUE = X * X * X - MH2O
    !!DMDT_INT = X*X*X/WR    !<step5.2> change calculation to keep WR constant after condensation/evap (win, 5/14/06)

    !<step6.3> bring back the previously reverted back (win, 10/2/08)
    VALUE = X*X*X/WR
    !pja Perform some numerical checks on dmdt_int
    IF ((TAU > 0.0) .and. (VALUE < M0)) VALUE = M0
    IF ((TAU < 0.0) .and. (VALUE > M0)) VALUE = M0

  END FUNCTION DMDT_INT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gasdiff
!
! !DESCRIPTION: Function GASDIFF returns the diffusion constant of a species in
!  air (m2/s). It uses the method of Fuller, Schettler, and Giddings as
!  described in Perry's Handbook for Chemical Engineers.
!  WRITTEN BY Peter Adams, May 2000
!\\
!\\
! !INTERFACE:
!
  FUNCTION GASDIFF( TEMP, PRES, MW, SV ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    real temp, pres  !temperature (K) and pressure (Pa) of air
    real mw          !molecular weight (g/mol) of diffusing species
    real Sv          !sum of atomic diffusion volumes of diffusing species
!
! !RETURN VALUE:
!
    real VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real mwair, Svair   !same as above, but for air
    real mwf, Svf
    parameter(mwair=28.9, Svair=20.1)

    !========================================================================
    ! GASDIFF begins here!
    !========================================================================

    mwf=sqrt((mw+mwair)/(mw*mwair))
    Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
    VALUE =1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf

  END FUNCTION GASDIFF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getdp
!
! !DESCRIPTION: Function GETDP calculate multi-component aerosol diameter
!  Originally written by Peter Adams in GISS GCM-II'
!  Use in GEOS-CHEM v5-07-08 and later by Win Trivitayanurak (win@cmu.edu)
!\\
!\\
! !INTERFACE:
!
  FUNCTION GETDP( I, J, L, N, State_Met, State_Chm, RC ) &
       RESULT( VALUE )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indices
    INTEGER,        INTENT(IN)    :: N           ! Species index
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !RETURN VALUE:
!
    REAL(fp)                      :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: NUMBIN,  ID,   JC
    REAL(fp)                :: MSO4, MNO3, MNH4, MH2O, MNACL
    REAL(fp)                :: MECIL, MECOB, MOCIL, MOCOB, MDUST
    REAL(fp)                :: DENSITY !density (kg/m3) of current size bin
    REAL(fp)                :: TOTALMASS !(kg)
    REAL(fp)                :: MCONC, NCONC
    CHARACTER(LEN=255)      :: MSG, LOC ! For species unit check (ewl)

    !real(fp), external :: aerodens

    REAL(fp)                ::  pi
    parameter (pi=3.141592654e+0_fp)

    ! Pointers
    ! We need to define local arrays to hold corresponding values
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER   :: Spc(:,:,:,:)

    !=================================================================
    ! GETDP begins here!
    !=================================================================

    ! Assume success
    RC                =  GC_SUCCESS

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    !-------------------------------------------------------------
    ! Calculate bin that we're working with
    !-------------------------------------------------------------
    NUMBIN = MOD(N-id_NK1+1,IBINS)
    IF (NUMBIN==0) NUMBIN = IBINS
    ID = id_NK1-1+NUMBIN   !ID = Species ID of number at current bin

    !-------------------------------------------------------------
    ! Calculate aerosol water in case it has not been initialized elsewhere
    !-------------------------------------------------------------
    CALL EZWATEREQM2( I, J, L, NUMBIN, State_Met, State_Chm, RC )

    !-------------------------------------------------------------
    ! Check negative Spc
    !-------------------------------------------------------------
    ! Significance limit in concentration unit
    ! Treshold for mass concentration 1d-4 ug/m3 = 1d-13 kg/m3
    ! Treshold for numb concentration 1d-1 #/cm3 = 1d5 kg/m3 (fake kg = #)
    IF( Spc(I,J,L,ID) == 0e+0_fp ) GOTO 10

    IF( Spc(I,J,L,ID) < 0e+0_fp ) THEN
       NCONC = ABS( Spc(I,J,L,ID) )/ State_Met%AIRVOL(I,J,L)/1e+6_fp
       IF ( nconc <= 1e-1_fp ) THEN
          Spc(I,J,L,ID) = 0e+0_fp
       ELSE
          PRINT *,'#### GETDP: negative NK at',I,J,L,'bin',NUMBIN
          PRINT *,'Species',N,'Spc=',Spc(I,J,L,ID)
          CALL ERROR_STOP('Negative NK', 'GETDP:1')
       ENDIF
    ENDIF
    IF(IT_IS_NAN(Spc(I,J,L,ID))) PRINT *,'+++++++++ Found Nan in ' , &
          'GETDP at (I,J,L)',I,J,L,'Bin',NUMBIN,': Nk'
    DO JC = 1, ICOMP-IDIAG
       IF( Spc(I,J,L,ID+JC*IBINS) < 0e+0_fp ) THEN
          MCONC = ( ABS(Spc(I,J,L,ID+JC*IBINS)) * 1.e+9_fp / &
                    State_Met%AIRVOL(I,J,L) )
          IF ( MCONC <= 1.e-4_fp ) THEN
             Spc(I,J,L,ID+JC*IBINS) = 0e+0_fp
          ELSE
             PRINT *,'#### GETDP: negative mass at',I,J,L,'bin',NUMBIN
             PRINT *,'Species',N,'Spc=',Spc(I,J,L,ID+JC*IBINS)
             CALL ERROR_STOP('Negative mass','GETDP:2')
          ENDIF
       ENDIF
       IF(IT_IS_NAN(Spc(I,J,L,ID+JC*IBINS))) PRINT *,'+++++++++ ', &
          'Found Nan in GETDP at (I,J,L)',I,J,L,'Bin',NUMBIN,'comp',JC
    ENDDO

    !-------------------------------------------------------------
    ! Begin calculation of diameter
    !-------------------------------------------------------------

    ! Totalmass is the total mass per particle (including water and ammonia)
    ! The factor of 0.1875 is the proportion of nh4 to make the particle
    ! ammonium bisulfate
    MSO4  = 0.e+0_fp
    MNACL = 0.e+0_fp
    MH2O  = 0.e+0_fp
    MECIL = 0.e+0_fp
    MECOB = 0.e+0_fp
    MOCIL = 0.e+0_fp
    MOCOB = 0.e+0_fp
    MDUST = 0.e+0_fp

    ! Get aerosol masses from GEOS-CHEM's Spc array
    DO JC = 1, ICOMP-IDIAG
       IF( JC == SRTSO4  ) MSO4  = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTNACL ) MNACL = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTECIL ) MECIL = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTECOB ) MECOB = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTOCIL ) MOCIL = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTOCOB ) MOCOB = Spc(I,J,L,ID+JC*IBINS)
       IF( JC == SRTDUST ) MDUST = Spc(I,J,L,ID+JC*IBINS)
    ENDDO
    MH2O  = Spc(I,J,L,id_AW1-1+NUMBIN)

    !dbg print *,'mh2o',mh2o,'at',i,j,l

    MNO3 = 0.e+0_fp
    MNH4 = MSO4 * 1.875e-1_fp
    TOTALMASS = ( MSO4 + MNO3 + MNH4 + MNACL + MH2O + &
                  MECIL + MECOB + MOCIL + MOCOB + MDUST)/ &
                  Spc(I,J,L,ID)
    DENSITY = AERODENS( MSO4, MNO3, MNH4, MNACL, MECIL, MECOB, &
                        MOCIL, MOCOB, MDUST, MH2O)

    VALUE = ( TOTALMASS* 6.e+0_fp/ DENSITY/ PI ) **(1.e+0_fp/3.e+0_fp) !getdp [=] meter

    GOTO 20

    !if number and mass is zero - calculate dp based on the density=1500 kg/m3
10  CONTINUE
    TOTALMASS = 1.414e+0_fp * Xk(NUMBIN)  ! Mid-bin mass per particle
    VALUE = ( TOTALMASS* 6.e+0_fp/ 1500.e+0_fp/ PI ) &
            **(1.e+0_fp/3.e+0_fp) !getdp [=] meter

20  CONTINUE

    IF( IT_IS_NAN( VALUE )) &
         CALL ERROR_STOP('Result is NaN', 'GETDP:3')

    ! Free pointer
    Spc => NULL()

  END FUNCTION GETDP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: spinup
!
! !DESCRIPTION: Function SPINUP retuns .TRUE. or .FALSE. whether or not the
!  current time in the run have passed the spin-up period.  This would be used
!  to determine if certain errors should be fixed and let slipped or to stop a
!  run with an error message.  (win, 8/2/07)
!  ====> Be cautious that TIMEBEGIN should be changed according to
!         whatever your spin-up beginning time is
!  Example of TIMEBEGIN (in julian time)
!         2001/07/01 = 144600.0
!         2000/11/01 = 138792.0
!\\
!\\
! !INTERFACE:
!
  FUNCTION SPINUP( DAYS ) RESULT( VALUE )
!
! !USES:
!
    USE TIME_MOD,     ONLY : GET_TAU , GET_TAUb
!
! !INPUT PARAMETERS:
!
    REAL*4,    INTENT(IN) :: DAYS   ! Spin-up duration (day)
!
! !RETURN VALUE:
!
    LOGICAL               :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*4                 :: TIMENOW, TIMEBEGIN, TIMEINIT, HOURS

    !========================================================================
    ! SPINUP begins here!
    !========================================================================

    TIMENOW   = GET_TAU()   ! Current time in the run (Julian time) (hrs)
    TIMEBEGIN = GET_TAUb()  ! Begin time of this run (hrs)
    TIMEINIT  = 141000. !2/1/2001    ! Start time for spin-up (hrs)
    HOURS = DAYS * 24.0     ! Period allow error to pass (hrs)

    ! Criteria to let error go or to terminate the run
    !IF ( TIMENOW > MIN( TIMEBEGIN, TIMEINIT ) + HOURS  ) THEN
    IF ( TIMENOW > TIMEBEGIN + HOURS  ) THEN
       VALUE = .FALSE.
    ELSE
       VALUE = .TRUE.
    ENDIF

  END FUNCTION SPINUP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodens
!
! !DESCRIPTION: Function STRATSCAV is basically a lookup table to get the below-
!  cloud scavenging rate (per mm of rainfall) as a function of particle
!  diameter.  The data are taken from Dana, M. T., and
!  J. M. Hales, Statistical Aspects of the Washout of Polydisperse
!  Aerosols, Atmos. Environ., 10, 45-50, 1976.  I am using the
!  monodisperse aerosol curve from Figure 2 which assumes a
!  lognormal distribution of rain drops with Rg=0.02 cm and a
!  sigma of 1.86, values typical of a frontal rain spectrum
!  (stratiform clouds).
!  WRITTEN BY Peter Adams, January 2001
!  Intoduced to GEOS-Chem by Win Trivitayanurak, 8/6/07
!\\
!\\
! !INTERFACE:
!
  FUNCTION STRATSCAV( DP ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp),   INTENT(IN)  :: DP  !particle diameter (m)
!
! !RETURN VALUE:
!
    REAL*4                  :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer numpts  !number of points in lookup table
    real dpdat      !particle diameter in lookup table (m)
    real scdat      !scavenging rate in lookup table (mm-1)
    integer n1, n2  !indices of nearest data points
    parameter(numpts=37)
    dimension dpdat(numpts), scdat(numpts)

    data dpdat/ 2.0E-09, 4.0E-09, 6.0E-09, 8.0E-09, 1.0E-08, &
                1.2E-08, 1.4E-08, 1.6E-08, 1.8E-08, 2.0E-08, &
                4.0E-08, 6.0E-08, 8.0E-08, 1.0E-07, 1.2E-07, &
                1.4E-07, 1.6E-07, 1.8E-07, 2.0E-07, 4.0E-07, &
                6.0E-07, 8.0E-07, 1.0E-06, 1.2E-06, 1.4E-06, &
                1.6E-06, 1.8E-06, 2.0E-06, 4.0E-06, 6.0E-06, &
                8.0E-06, 1.0E-05, 1.2E-05, 1.4E-05, 1.6E-05, &
                1.8E-05, 2.0E-05/

    data scdat/ 6.99E-02, 2.61E-02, 1.46E-02, 9.67E-03, 7.07E-03, &
                5.52E-03, 4.53E-03, 3.87E-03, 3.42E-03, 3.10E-03, &
                1.46E-03, 1.08E-03, 9.75E-04, 9.77E-04, 1.03E-03, &
                1.11E-03, 1.21E-03, 1.33E-03, 1.45E-03, 3.09E-03, &
                4.86E-03, 7.24E-03, 1.02E-02, 1.36E-02, 1.76E-02, &
                2.21E-02, 2.70E-02, 3.24E-02, 4.86E-01, 8.36E-01, &
                1.14E+00, 1.39E+00, 1.59E+00, 1.75E+00, 1.85E+00, &
                1.91E+00, 1.91E+00/

    !=================================================================
    ! STRATSCAV begins here!
    !=================================================================

    ! If particle diameter is in bounds, interpolate to find value
    if ((dp .gt. dpdat(1)) .and. (dp .lt. dpdat(numpts))) then
       !loop over lookup table points to find nearest values
       n1=1
       do while (dp .gt. dpdat(n1+1))
          n1=n1+1
       enddo
       n2=n1+1
       VALUE=scdat(n1)+(scdat(n2)-scdat(n1)) &
             *(dp-dpdat(n1))/(dpdat(n2)-dpdat(n1))
    endif

    ! If particle diameter is out of bounds, return reasonable value
    if (dp .gt. dpdat(numpts)) VALUE=2.0
    if (dp .lt. dpdat(1))      VALUE=7.0e-2

  END FUNCTION STRATSCAV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: waternacl
!
! !DESCRIPTION: Function WATERNACL uses the current RH to calculate how much
!  water is in equilibrium with the seasalt.  Aerosol water concentrations are
!  assumed to be in equilibrium at all times and the array of concentrations is
!  updated accordingly.
!  WRITTEN BY Peter Adams, November 2001
!  Introduced to GEOS-CHEM by Win Trivitayanurak. 8/6/07
!\\
!\\
! !INTERFACE:
!
  FUNCTION WATERNACL( RHE ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp) :: RHE ! Relative humidity (0-100 scale)
!
! !RETURN VALUE:
!
    REAL(fp) :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! WATERNACL begins here!
    !=================================================================

    if (rhe .gt. 99.) rhe=99.
    if (rhe .lt. 1.) rhe=1.

    if (rhe .gt. 90.) then
       VALUE=5.1667642e-2*rhe**3-14.153121*rhe**2+1292.8377*rhe-3.9373536e4
    else
       if (rhe .gt. 80.) then
          VALUE=1.0629e-3*rhe**3-0.25281*rhe**2+20.171*rhe-5.3558e2
       else
          if (rhe .gt. 50.) then
             VALUE=4.2967e-5*rhe**3-7.3654e-3*rhe**2+.46312*rhe-7.5731
          else
             if (rhe .gt. 20.) then
                VALUE=2.9443e-5*rhe**3-2.4739e-3*rhe**2+7.3430e-2*rhe+1.3727
             else
                VALUE=1.17
             endif
          endif
       endif
    endif

    !check for error
    if (VALUE .gt. 45.) then
       write(*,*) 'ERROR in waternacl'
       write(*,*) rhe,VALUE
       STOP
    endif

  END FUNCTION WATERNACL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: waterocil
!
! !DESCRIPTION: Function WATEROCIL uses the current RH to calculate how much
!  water is in equilibrium with the hydrophillic OA.  Aerosol water
!  concentrations are assumed to be in equilibrium at all times and the array of
!  concentrations is updated accordingly.
!  MODIFIED BY YUNHA LEE, AUG, 2006
!  Bring to GEOS-CHEM by Win Trivitayanurak 9/3/07
!\\
!\\
! !INTERFACE:
!
  FUNCTION WATEROCIL( RHE ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp) :: RHE ! Relative humidity (0-100 scale)
!
! !RETURN VALUE:
!
    REAL(fp) :: VALUE
!
! !REMARKS:
!  waterocil is the ratio of wet mass to dry mass of a particle.  Instead
!  of calling a thermodynamic equilibrium code, this routine uses a
!  simple curve fit to estimate waterocil based on the current humidity.
!  The curve fit is based on observations of Dick et al. JGR D1 1471-1479
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: a, b, c, d, e, f, prefactor, activcoef
    parameter(a=1.0034, b=0.1614, c=1.1693,d=-3.1,e=6.0)

    !=================================================================
    ! WATEROCIL begins here!
    !=================================================================

    if (rhe .gt. 99.) rhe=99.
    if (rhe .lt. 1.) rhe=1.

    if (rhe .gt. 85.) then
       VALUE =d+e*(rhe/100)
       !yhl Growth factor above RH 85% is not available, so it assumes linear
       !yhl growth at above 85%.
    else
       VALUE =a+b*(rhe/100)+c*(rhe/100)**2.
       !yhl This eq is based on the extrapolation curve obtained from
       !yhl Dick et al 2000 figure 5.(High organic,density=1400g/cm3)
    endif

    !check for error
    if (VALUE .gt. 10.) then
       write(*,*) 'ERROR in waterocil'
       write(*,*) rhe, value
       STOP
    endif

    RETURN

  END FUNCTION WATEROCIL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: waterso4
!
! !DESCRIPTION: Function WATERSO4 uses the current RH to calculate how much
!  water is in equilibrium with the sulfate.  Aerosol water concentrations are
!  assumed to be in equilibrium at all times and the array of concentrations is
!  updated accordingly.
!   Introduced to GEOS-CHEM by Win Trivitayanurak. 8/6/07
!   Adaptation of ezwatereqm used in size-resolved sulfate only sim
!   November, 2001
!   ezwatereqm WRITTEN BY Peter Adams, March 2000
!\\
!\\
! !INTERFACE:
!
  FUNCTION WATERSO4( RHE ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp) :: RHE ! Relative humidity (0-100 scale)
!
! !RETURN VALUE:
!
    REAL(fp) :: VALUE

! !REMARKS:
!  waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
!  of calling a thermodynamic equilibrium code, this routine uses a
!  simple curve fit to estimate wr based on the current humidity.
!  The curve fit is based on ISORROPIA results for ammonium bisulfate
!  at 273 K.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! WATERSO4 begins here!
    !=================================================================

    if (rhe .gt. 99.) rhe=99.
    if (rhe .lt. 1.) rhe=1.

    if (rhe .gt. 96.) then
       value=0.7540688*rhe**3-218.5647*rhe**2+21118.19*rhe-6.801999e5
    else
       if (rhe .gt. 91.) then
          value=8.517e-2*rhe**2 -15.388*rhe +698.25
       else
          if (rhe .gt. 81.) then
             value=8.2696e-3*rhe**2 -1.3076*rhe +53.697
          else
             if (rhe .gt. 61.) then
                value=9.3562e-4*rhe**2 -0.10427*rhe +4.3155
             else
                if (rhe .gt. 41.) then
                   value=1.9149e-4*rhe**2 -8.8619e-3*rhe +1.2535
                else
                   value=5.1337e-5*rhe**2 +2.6266e-3*rhe +1.0149
                endif
             endif
          endif
       endif
    endif
    
    !check for error
    if (value .gt. 30.) then
       write(*,*) 'ERROR in waterso4'
       write(*,*) rhe,value
       STOP
    endif

  END FUNCTION WATERSO4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cf_nucl
!
! !DESCRIPTION: This subroutine calculates the barrierless nucleation rate and
!  radius of the critical nucleation cluster using the parameterization of...
!     Clement and Ford (1999) Atmos. Environ. 33:489-499
!     WRITTEN BY Jeff Pierce, April 2007
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cf_nucl(tempi,rhi,cna,nh3ppt,fn)
!
! !INPUT PARAMETERS:
!
    real tempi                ! temperature of air [K]
    real rhi                  ! relative humidity of air as a fraction
    double precision cna      ! concentration of gas phase sulfuric acid [molec cm-3]
    double precision nh3ppt   ! mixing ratio of ammonia in ppt
!
! !OUTPUT PARAMETERS:
!
    double precision fn                   ! nucleation rate [cm-3 s-1]
    double precision rnuc                 ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision temp                 ! temperature of air [K]
    double precision rh                   ! relative humidity of air as a fraction
    double precision alpha1

    temp=dble(tempi)
    rh=dble(rhi)

    if (nh3ppt .lt. 0.1) then
       alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
    else
       alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
    endif
    fn = alpha1*cna**2*3600.
    ! sensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
    if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion

10  return

  end subroutine cf_nucl
!EOC
!--------------------------------------------------------------------
!
!        Ion nucleation Rate calculation
!        from Modgil et al.(2005), JGR, vol. 110, D19205
!
!--------------------------------------------------------------------
! Mathematical Expressions for particle nucleation rate (h1,cm-3 s-1),
! nucleating H2SO4 flux (h2, cm-3 s-1), number of H2SO4 in average
! nucleating cluster (h3), number of H2O in average nucleating cluster (h4),
! radius of average nucleating cluster (h5, nm), and first order loss of
! H2SO4 to particles (h6) are given below,

  subroutine ion_nucl(h2so4i,sai,ti,qi,rhi,h1,h2,h3,h4,h5,h6)

    ! h2so4 = h2so4 concentration [molec cm^-3]
    ! sa = aerosol surface area [um^2 cm^-3]
    ! t = temperature [K]
    ! q = ion formation rate [ion pairs cm^-3 s-1]
    ! rh = relative humidity as a fraction

    real ti, rhi
    double precision h2so4i,qi,sai
    double precision t,rh,h2so4,q,sa
    double precision h1,h2,h3,h4,h6,h5

    t=dble(ti)
    rh=dble(rhi)
    h2so4=h2so4i
    q=qi
    sa=sai

    if (h2so4.lt.1.e+5_fp.or.t.gt.260.e+0_fp) then ! changed to 2E5 because function was diverging
       ! diverges above 260K
       h1=0.e+0_fp
       h2=0.e+0_fp
       h3=0.e+0_fp
       h4=0.e+0_fp
       h5=0.e+0_fp
       h6=0.e+0_fp
       return
    endif
    if (h2so4 .gt. 1e+8_fp) h2so4 = 1e+8_fp
    if (sa .lt. 2.e+0_fp) sa = 2.e+0_fp
    if (sa .gt. 100.e+0_fp) sa = 100.e+0_fp
    if (t .lt. 190.e+0_fp) t = 190.e+0_fp
    if (t .gt. 300.e+0_fp) t = 300.e+0_fp
    if (rh .lt. 0.05e+0_fp) rh = 0.05e+0_fp
    if (rh .gt. 0.95e+0_fp) rh = 0.95e+0_fp
    if (q .lt. 1.e+0_fp) q = 1.e+0_fp
    if (q .gt. 50.e+0_fp) q = 50.e+0_fp

    h6=0.000026859579119003205*SA + 1.7477354270484002e-8_fp*q*SA + &
       1.5718068902491457e-8_fp*SA**2+8.060796806911441e-8_fp*SA*T + &
       3.904048293417882e-7_fp*SA*Log(H2SO4) + &
       2.727259306977938e-7_fp*SA*Log(RH)


    h3=-198.8039518313554 + 3357.132963009284*h6 - 130510.31325149858* &
       h6**2 - 0.7093715033997716*q - 10.713505046150196*h6*q + &
       1103.4737682776713*h6**2*q + 0.0052565148186649*q**2 - &
       0.20195850414426988*h6*q**2 + 10.961027676935213*h6**2*q**2 - &
       26.553841269634976*RH + 2913.499196899548*h6*RH - &
       7558.996305824136*h6**2*RH + 0.050092880471591994*q*RH + &
       0.39840936335061017*h6*q*RH + 16.140386509938388*h6**2*q*RH - &
       0.0008159572217147427*q**2*RH - 0.02492462618304389*h6*q**2*RH + &
       3.2372842210428825*RH**2 + 1709.7485838150235*h6*RH**2 - &
       4016.182638678486*h6**2*RH**2 - 0.022142010235491123*q*RH**2 - &
       1.620063009925805*h6*q*RH**2-0.00028477984814528825*q**2*RH**2 + &
       22.136724153656015*RH**3 + 170.8982375938333*h6*RH**3 - &
       0.01881686723215867*q*RH**3 + 2.6974144456100957*T - &
       96.60591604496483*h6*T + 1772.137264721083*h6**2*T - &
       0.0009432251807207652*q*T + 0.06072064184950673*h6*q*T - &
       2.5196932894429502*h6**2*q*T - 0.000013848768113392552*q**2*T - &
       0.0001948394841164792*h6*q**2*T + 0.1828636512279507*RH*T - &
       55.135341874839185*h6*RH*T + 164.02631709083576*h6**2*RH*T + &
       0.001745921048607296*q*RH*T + 0.035017713828742754*h6*q*RH*T + &
       4.057082638293583e-6_fp*q**2*RH*T - 0.3900441693913758*RH**2*T - &
       8.955078982582657*h6*RH**2*T+0.00021434974336412074*q*RH**2*T - &
       0.14947568974964962*RH**3*T - 0.022748394377623382*T**2 + &
       0.7227721843282815*h6*T**2 - 5.386480671871808*h6**2*T**2 - &
       0.000035250836279611095*q*T**2-0.0003363405774846326*h6*q*T**2+ &
       2.9254973516794257d-8*q**2*T**2 + 0.003994529829421164*RH*T**2+ &
       0.2074067980035454*h6*RH*T**2-5.136172472264946e-6_fp*q*RH*T**2+ &
       0.0020603328018819816*RH**2*T**2+0.000042019279193164354*T**3 - &
       0.002100661388749787*h6*T**3+7.309966632740304e-8_fp*q*T**3 - &
       0.000016323969052607556*RH*T**3+ 12.330627568298462*Log(H2SO4) + &
       768.3961789008589*h6*Log(H2SO4)-11568.47324943553*h6**2* &
       Log(H2SO4)+0.14349043416922366*q*Log(H2SO4)+0.8946851157223353* &
       h6*q*Log(H2SO4) - 68.46004143191098*h6**2*q*Log(H2SO4) - &
       0.0006241121793370407*q**2*Log(H2SO4)+0.011897674833907721*h6* &
       q**2*Log(H2SO4) + 1.7860574934328677*RH*Log(H2SO4) + &
       316.5406316191978*h6*RH*Log(H2SO4) - 2036.825340216443*h6**2*RH* &
       Log(H2SO4) - 0.026323914507434605*q*RH*Log(H2SO4) - &
       0.37505804775954393*h6*q*RH*Log(H2SO4) + 0.00003454867680790666* &
       q**2*RH*Log(H2SO4) + 2.844877302874606*RH**2*Log(H2SO4) - &
       1.5845895178086176*h6*RH**2*Log(H2SO4) + 0.001732608008714275*q* &
       RH**2*Log(H2SO4) + 0.5611003862827533*RH**3*Log(H2SO4) + &
       0.18033151281768975*T*Log(H2SO4) - 7.807090214680351*h6*T* &
       Log(H2SO4) + 52.76241348342321*h6**2*T*Log(H2SO4) + &
       0.0011535134888316242*q*T*Log(H2SO4) + 0.0068466874844708295*h6* &
       q*T*Log(H2SO4) - 4.231224168766194e-7_fp*q**2*T*Log(H2SO4) - &
       0.10505349775719895*RH*T*Log(H2SO4) - 1.9241106452950727*h6*RH* &
       T*Log(H2SO4) + 2.0451815715440337e-6_fp*q*RH*T*Log(H2SO4) - &
       0.015299483302534183*RH**2*T*Log(H2SO4) + 0.000775633115370002* &
       T**2*Log(H2SO4) + 0.04228608723566267*h6*T**2*Log(H2SO4) - &
       9.572221299803945e-7_fp*q*T**2*Log(H2SO4) + &
       0.0002669785990812474* &
       RH*T**2*Log(H2SO4) - 1.7595742533055222e-6_fp*T**3*Log(H2SO4) - &
       2.7165200489046812*Log(H2SO4)**2 + 1.963036665672805*h6* &
       Log(H2SO4)**2 + 33.88545004559797*h6**2*Log(H2SO4)**2 - &
       0.01647722099703982*q*Log(H2SO4)**2 - 0.08498050322218324*h6*q* &
       Log(H2SO4)**2 + 0.00003320475358154802*q**2*Log(H2SO4)**2 + &
       0.50626436977659*RH*Log(H2SO4)**2 + 3.914586690682404*h6*RH* &
       Log(H2SO4)**2 + 0.0006654515705980484*q*RH*Log(H2SO4)**2 - &
       0.025122425208873058*RH**2*Log(H2SO4)**2 - 0.021380868797539664* &
       T*Log(H2SO4)**2 - 0.35405514523597137*h6*T*Log(H2SO4)**2 - &
       0.000020639290758942666*q*T*Log(H2SO4)**2+0.0003587951811915662* &
       RH*T*Log(H2SO4)**2+7.620708111644729e-6_fp*T**2*Log(H2SO4)**2 + &
       0.2196641696573127*Log(H2SO4)**3 + 1.7291708055805226*h6* &
       Log(H2SO4)**3 + 0.00038146321602426414*q*Log(H2SO4)**3 - &
       0.011997306640487447*RH*Log(H2SO4)**3 + 0.0003857955558500776*T* &
       Log(H2SO4)**3 - 0.004779827937902779*Log(H2SO4)**4

    h3=EXP(h3)

    h1=456229.3726785317 - 696754.0061755505/h3 - &
       8.954389043957226e+7_fp*h6 + (1.4677717736521986e+8_fp*h6)/h3 + &
       1867.5296995211318*q - (2798.172491398116*q)/h3 + &
       1500.05530404756*h6*q - (171625.68387665015*h6*q)/h3 - &
       5924.898937400813*T + (7657.054762118453*T)/h3 + &
       1.1074613167489376e+6_fp*h6*T-(1.556115054313954e+6_fp*h6*T)/h3- &
       28.49469750360138*q*T + (50.46283217861981*q*T)/h3 + &
       1753.7251886642075*h6*q*T - (975.4050494746682*h6*q*T)/h3 + &
       25.45577410331681*T**2 - (26.190051353137182*T**2)/h3 - &
       4516.4884641323815*h6*T**2 + (5030.318058537955*h6*T**2)/h3 + &
       0.14403307295844858*q*T**2 - (0.2897062226196713*q*T**2)/h3 - &
       15.561717210345154*h6*q*T**2+(18.061492698402766*h6*q*T**2)/h3- &
       0.03621852394282133*T**3 + (0.02659567895141871*T**3)/h3 + &
       6.0696427040893886*h6*T**3 - (4.585112604562606*h6*T**3)/h3 - &
       0.00024142548506534756*q*T**3 + &
       (0.0005400351895947918*q*T**3)/h3+0.0341853829283136*h6*q*T**3- &
       (0.04459905758423761*h6*q*T**3)/h3-81941.3895777097*Log(H2SO4)+ &
       (90462.54571516594*Log(H2SO4))/h3 + &
       1.6483362741354546e+7_fp*h6*Log(H2SO4) - &
       (2.3542460879418086e+7_fp*h6*Log(H2SO4))/h3 - &
       319.8790300996385*q*Log(H2SO4) + &
       (226.84336756604796*q*Log(H2SO4))/h3 - &
       26455.51148372377*h6*q*Log(H2SO4) + &
       (60326.79915791394*h6*q*Log(H2SO4))/h3 + &
       1063.8621634534027*T*Log(H2SO4) - &
       (891.6122304868614*T*Log(H2SO4))/h3 - &
       203718.41653995324*h6*T*Log(H2SO4) + &
       (243001.1509313233*h6*T*Log(H2SO4))/h3 + &
       4.963508394566835*q*T*Log(H2SO4) - &
       (5.401955915120605*q*T*Log(H2SO4))/h3 + &
       25.811543727271804*h6*q*T*Log(H2SO4) - &
       (355.292262810142*h6*q*T*Log(H2SO4))/h3 - &
       4.5691538346459195*T**2*Log(H2SO4) + &
       (2.468423948121569*T**2*Log(H2SO4))/h3 + &
       830.195816015575*h6*T**2*Log(H2SO4) - &
       (745.2009340266254*h6*T**2*Log(H2SO4))/h3 - &
       0.025450884793884514*q*T**2*Log(H2SO4) + &
       (0.03593327468289806*q*T**2*Log(H2SO4))/h3 + &
       1.3302486952619124*h6*q*T**2*Log(H2SO4) - &
       (0.320851018492291*h6*q*T**2*Log(H2SO4))/h3 + &
       0.006498242231697865*T**3*Log(H2SO4) - &
       (0.0013540108127227987*T**3*Log(H2SO4))/h3 - &
       1.1150374268629073*h6*T**3*Log(H2SO4) + &
       (0.5953715882656012*h6*T**3*Log(H2SO4))/h3 + &
       0.000043146778270078744*q*T**3*Log(H2SO4) - &
       (0.0000735656618351653*q*T**3*Log(H2SO4))/h3 - &
       0.004059155580949069*h6*q*T**3*Log(H2SO4) + &
       (0.0028886014749212705*h6*q*T**3*Log(H2SO4))/h3 + &
       4912.2066512397205*Log(H2SO4)**2 - &
       (3022.0044642322964*Log(H2SO4)**2)/h3 - &
       1.0095170363345986d6*h6*Log(H2SO4)**2 + &
       (1.1998313113698033d6*h6*Log(H2SO4)**2)/h3 + &
       18.018115617541998*q*Log(H2SO4)**2 + &
       (7.561352267717425*q*Log(H2SO4)**2)/h3 + &
       3133.188395383639*h6*q*Log(H2SO4)**2 - &
       (3338.848456178725*h6*q*Log(H2SO4)**2)/h3 - &
       63.77612570995598*T*Log(H2SO4)**2 + &
       (20.100257708705428*T*Log(H2SO4)**2)/h3 + &
       12469.209965937996*h6*T*Log(H2SO4)**2 - &
       (11883.666207080118*h6*T*Log(H2SO4)**2)/h3 - &
       0.2853187730266578*q*T*Log(H2SO4)**2 + &
       (0.03880844916086884*q*T*Log(H2SO4)**2)/h3 - &
       21.708027380830988*h6*q*T*Log(H2SO4)**2 + &
       (28.111166854183463*h6*q*T*Log(H2SO4)**2)/h3 + &
       0.27388018029546984*T**2*Log(H2SO4)**2 + &
       (0.005568448436161079*T**2*Log(H2SO4)**2)/h3 - &
       50.78398945757657*h6*T**2*Log(H2SO4)**2 + &
       (33.25484929904544*h6*T**2*Log(H2SO4)**2)/h3 + &
       0.0014875097131454307*q*T**2*Log(H2SO4)**2 - &
       (0.0008810798671059988*q*T**2*Log(H2SO4)**2)/h3 + &
       0.007045532214236005*h6*q*T**2*Log(H2SO4)**2 - &
       (0.05521643331038819*h6*q*T**2*Log(H2SO4)**2)/h3 - &
       0.00038943745406964085*T**3*Log(H2SO4)**2 - &
       (0.00015308439633013156*T**3*Log(H2SO4)**2)/h3 + &
       0.06817761559886615*h6*T**3*Log(H2SO4)**2 - &
       (0.019520711865759842*h6*T**3*Log(H2SO4)**2)/h3 - &
       2.554240267270117e-6_fp*q*T**3*Log(H2SO4)**2 + &
       (2.526957991768314e-6_fp*q*T**3*Log(H2SO4)**2)/h3 + &
       0.00011978108924107446*h6*q*T**3*Log(H2SO4)**2 - &
       98.09290363960072*Log(H2SO4)**3 + &
       (3.8697214915504805*Log(H2SO4)**3)/h3 + &
       20559.743351544734*h6*Log(H2SO4)**3 - &
       (18687.867717989076*h6*Log(H2SO4)**3)/h3 - &
       0.33118998078526996*q*Log(H2SO4)**3 - &
       (0.6955616584928909*q*Log(H2SO4)**3)/h3 - &
       92.45613590327787*h6*q*Log(H2SO4)**3 + &
       (8.948488896745264*h6*q*Log(H2SO4)**3)/h3 + &
       1.273863980203035*T*Log(H2SO4)**3 + &
       (0.38051947439108685*T*Log(H2SO4)**3)/h3 - &
       253.8228380167286*h6*T*Log(H2SO4)**3 + &
       (171.19178070061787*h6*T*Log(H2SO4)**3)/h3 + &
       0.005376622635136611*q*T*Log(H2SO4)**3 + &
       (0.006658042079999949*q*T*Log(H2SO4)**3)/h3 + &
       0.8230892248689082*h6*q*T*Log(H2SO4)**3 - &
       (0.058021973134188776*h6*q*T*Log(H2SO4)**3)/h3 - &
       0.005471105747328322*T**2*Log(H2SO4)**3 - &
       (0.003699167378319841*T**2*Log(H2SO4)**3)/h3 + &
       1.0332570198673408*h6*T**2*Log(H2SO4)**3 - &
       (0.38600322062471826*h6*T**2*Log(H2SO4)**3)/h3 - &
       0.000028584962254144285*q*T**2*Log(H2SO4)**3 - &
       (0.000016062832568931284*q*T**2*Log(H2SO4)**3)/h3 - &
       0.001819022861961261*h6*q*T**2*Log(H2SO4)**3 + &
       7.779696142330516e-6_fp*T**3*Log(H2SO4)**3 + &
       (8.513832613126146e-6_fp*T**3*Log(H2SO4)**3)/h3 - &
       0.0013866861532058*h6*T**3*Log(H2SO4)**3 + &
       4.981099027173453e-8_fp*q*T**3*Log(H2SO4)**3 + &
       178821.264938151*Log(RH) - (634125.5419575685*Log(RH))/h3 - &
       2.310391119522488e+7_fp*h6*Log(RH) + &
       (3.564836811852359e+7_fp*h6*Log(RH))/h3 - &
       744.9440395010726*q*Log(RH)+(1019.4951846825716*q*Log(RH))/h3- &
       28008.93484108728*h6*q*Log(RH) - &
       (7136.405275749431*h6*q*Log(RH))/h3 - &
       2275.393426900043*T*Log(RH) + (9231.8244493922*T*Log(RH))/h3 + &
       237571.37589421016*h6*T*Log(RH) - &
       (295411.0553218543*h6*T*Log(RH))/h3 + &
       4.451622613453289*q*T*Log(RH) - &
       (4.150562571935793*q*T*Log(RH))/h3 + &
       595.8398944235868*h6*q*T*Log(RH) - &
       (210.09033695512161*h6*q*T*Log(RH))/h3 + &
       9.678343416541734*T**2*Log(RH) - &
       (45.69882758602488*T**2*Log(RH))/h3 - &
       747.9878102876322*h6*T**2*Log(RH) + &
       (649.4424346377774*h6*T**2*Log(RH))/h3 + &
       0.005153649541788936*q*T**2*Log(RH) - &
       (0.0023655719118541533*q*T**2*Log(RH))/h3 - &
       2.848932474900831*h6*q*T**2*Log(RH) + &
       (0.2889680076326506*h6*q*T**2*Log(RH))/h3 - &
       0.013824587420637588*T**3*Log(RH) + &
       (0.07673955581386811*T**3*Log(RH))/h3 + &
       0.6813253055163031*h6*T**3*Log(RH) - &
       (0.2663343569482817*h6*T**3*Log(RH))/h3 - &
       0.00004364745106737275*q*T**3*Log(RH) - &
       (9.434880576275503e-6_fp*q*T**3*Log(RH))/h3 + &
       0.0030430214338531205*h6*q*T**3*Log(RH) + &
       (0.00371157673745342*h6*q*T**3*Log(RH))/h3 - &
       31842.277946766168*Log(H2SO4)*Log(RH) + &
       (68134.57926920953*Log(H2SO4)*Log(RH))/h3 + &
       3.493985612508899e+6_fp*h6*Log(H2SO4)*Log(RH) - &
       (4.587173661378969e+6_fp*h6*Log(H2SO4)*Log(RH))/h3 + &
       190.19310993267845*q*Log(H2SO4)*Log(RH) - &
       (233.68205749312318*q*Log(H2SO4)*Log(RH))/h3 - &
       1547.3248018686509*h6*q*Log(H2SO4)*Log(RH) + &
       (4045.133810971437*h6*q*Log(H2SO4)*Log(RH))/h3 + &
       405.31725032582267*T*Log(H2SO4)*Log(RH) - &
       (1057.655396923827*T*Log(H2SO4)*Log(RH))/h3 - &
       33903.90722480216*h6*T*Log(H2SO4)*Log(RH) + &
       (32965.76708786528*h6*T*Log(H2SO4)*Log(RH))/h3 - &
       1.5034587825111836*q*T*Log(H2SO4)*Log(RH) + &
       (1.106270806420077*q*T*Log(H2SO4)*Log(RH))/h3 - &
       25.549725610666396*h6*q*T*Log(H2SO4)*Log(RH) + &
       (20.09072468984517*h6*q*T*Log(H2SO4)*Log(RH))/h3 - &
       1.7253372988643991*T**2*Log(H2SO4)*Log(RH) + &
       (5.612615323997249*T**2*Log(H2SO4)*Log(RH))/h3 + &
       95.62445916373963*h6*T**2*Log(H2SO4)*Log(RH) - &
       (48.96291394056088*h6*T**2*Log(H2SO4)*Log(RH))/h3 + &
       0.0019300142376745193*q*T**2*Log(H2SO4)*Log(RH) + &
       (0.00015584582354154648*q*T**2*Log(H2SO4)*Log(RH))/h3 + &
       0.19177967246162606*h6*q*T**2*Log(H2SO4)*Log(RH) - &
       (0.16871358892301552*h6*q*T**2*Log(H2SO4)*Log(RH))/h3 + &
       0.00246781658651071*T**3*Log(H2SO4)*Log(RH) - &
       (0.01009444649935961*T**3*Log(H2SO4)*Log(RH))/h3 - &
       0.06589616928988701*h6*T**3*Log(H2SO4)*Log(RH) - &
       (0.01596274185896012*h6*T**3*Log(H2SO4)*Log(RH))/h3 + &
       4.109041091272658e-6_fp*q*T**3*Log(H2SO4)*Log(RH) + &
       (2.041743340493583e-7_fp*q*T**3*Log(H2SO4)*Log(RH))/h3 - &
       0.0001572125066332154*h6*q*T**3*Log(H2SO4)*Log(RH) + &
       1932.4011188760162*Log(H2SO4)**2*Log(RH) - &
       (781.4228432455719*Log(H2SO4)**2*Log(RH))/h3 - &
       181291.39888342103*h6*Log(H2SO4)**2*Log(RH) + &
       (210184.24050416853*h6*Log(H2SO4)**2*Log(RH))/h3 - &
       13.666906354559005*q*Log(H2SO4)**2*Log(RH) + &
       (15.930151624229865*q*Log(H2SO4)**2*Log(RH))/h3 + &
       305.7703887043172*h6*q*Log(H2SO4)**2*Log(RH) - &
       (391.71164595591085*h6*q*Log(H2SO4)**2*Log(RH))/h3 - &
       24.621398404913865*T*Log(H2SO4)**2*Log(RH) + &
       (20.079316822912197*T*Log(H2SO4)**2*Log(RH))/h3 + &
       1658.5845815289838*h6*T*Log(H2SO4)**2*Log(RH) - &
       (1404.5530139788664*h6*T*Log(H2SO4)**2*Log(RH))/h3 + &
       0.11823318365930788*q*T*Log(H2SO4)**2*Log(RH) - &
       (0.08056635668752316*q*T*Log(H2SO4)**2*Log(RH))/h3 - &
       1.0034435154436134*h6*q*T*Log(H2SO4)**2*Log(RH) + &
       (1.6021411957310931*h6*q*T*Log(H2SO4)**2*Log(RH))/h3 + &
       0.10491613716181539*T**2*Log(H2SO4)**2*Log(RH) - &
       (0.1492097929648805*T**2*Log(H2SO4)**2*Log(RH))/h3 - &
       4.0934834940398215*h6*T**2*Log(H2SO4)**2*Log(RH) + &
       (1.896735315742051*h6*T**2*Log(H2SO4)**2*Log(RH))/h3 - &
       0.00022870693934890214*q*T**2*Log(H2SO4)**2*Log(RH) + &
       (0.000013175311622036916*q*T**2*Log(H2SO4)**2*Log(RH))/h3- &
       0.00231560013471937*h6*q*T**2*Log(H2SO4)**2*Log(RH) - &
       0.00015023186988448507*T**3*Log(H2SO4)**2*Log(RH) + &
       (0.0003388679191513823*T**3*Log(H2SO4)**2*Log(RH))/h3 + &
       0.0015784630518057667*h6*T**3*Log(H2SO4)**2*Log(RH) - &
       1.0056424638350849e-7_fp*q*T**3*Log(H2SO4)**2*Log(RH) - &
       39.74503165968021*Log(H2SO4)**3*Log(RH) - &
       (68.29220045487304*Log(H2SO4)**3*Log(RH))/h3 + &
       3236.6794930454566*h6*Log(H2SO4)**3*Log(RH) - &
       (2511.4313482742737*h6*Log(H2SO4)**3*Log(RH))/h3 + &
       0.3044230466506939*q*Log(H2SO4)**3*Log(RH) - &
       (0.321694018076018*q*Log(H2SO4)**3*Log(RH))/h3 - &
       7.83632554091337*h6*q*Log(H2SO4)**3*Log(RH) + &
       (1.2615699211576334*h6*q*Log(H2SO4)**3*Log(RH))/h3 + &
       0.5071743553774433*T*Log(H2SO4)**3*Log(RH) + &
       (0.6993184748948513*T*Log(H2SO4)**3*Log(RH))/h3 - &
       28.064969097031586*h6*T*Log(H2SO4)**3*Log(RH) + &
       (12.010607450658807*h6*T*Log(H2SO4)**3*Log(RH))/h3 - &
       0.0027511166431233884*q*T*Log(H2SO4)**3*Log(RH) + &
       (0.0015662272224653852*q*T*Log(H2SO4)**3*Log(RH))/h3 + &
       0.03815288545975011*h6*q*T*Log(H2SO4)**3*Log(RH) - &
       0.0021639947693271114*T**2*Log(H2SO4)**3*Log(RH) - &
       (0.0017767382579465334*T**2*Log(H2SO4)**3*Log(RH))/h3 + &
       0.059762491440998364*h6*T**2*Log(H2SO4)**3*Log(RH) + &
       6.130388542051968e-6_fp*q*T**2*Log(H2SO4)**3*Log(RH) + &
       3.1017462628108047e-6_fp*T**3*Log(H2SO4)**3*Log(RH) - &
       7689.12786121461*Log(RH)**2-(25856.380361656236*Log(RH)**2)/h3- &
       1.3915437178042033e+6_fp*h6*Log(RH)**2 + &
       (1.9928171556602388e+6_fp*h6*Log(RH)**2)/h3 + &
       53.79216985973448*q*Log(RH)**2 - &
       (134.8787674139501*q*Log(RH)**2)/h3 - &
       1440.4837667461468*h6*q*Log(RH)**2 - &
       (3141.0130344871345*h6*q*Log(RH)**2)/h3 + &
       94.7545268444358*T*Log(RH)**2 + &
       (342.5857682014369*T*Log(RH)**2)/h3 + &
       11341.294070320842*h6*T*Log(RH)**2 - &
       (18070.160441620144*h6*T*Log(RH)**2)/h3 - &
       0.44696864139255954*q*T*Log(RH)**2 + &
       (1.1199451205710647*q*T*Log(RH)**2)/h3 + &
       19.01601691757823*h6*q*T*Log(RH)**2 + &
       (31.604766584775273*h6*q*T*Log(RH)**2)/h3 - &
       0.3470127737474881*T**2*Log(RH)**2 - &
       (1.4877879003131458*T**2*Log(RH)**2)/h3 - &
       18.727169614346867*h6*T**2*Log(RH)**2 + &
       (42.5448759900235*h6*T**2*Log(RH)**2)/h3 + &
       0.0007727701515171791*q*T**2*Log(RH)**2 - &
       (0.0032041999804810098*q*T**2*Log(RH)**2)/h3 - &
       0.055631814243943*h6*q*T**2*Log(RH)**2 - &
       (0.07696588905487615*h6*q*T**2*Log(RH)**2)/h3 + &
       0.00033600378567117826*T**3*Log(RH)**2 + &
       (0.002069337837641118*T**3*Log(RH)**2)/h3 - &
       0.015788853514098967*h6*T**3*Log(RH)**2 - &
       (0.009301939734036398*h6*T**3*Log(RH)**2)/h3 + &
       5.572064421478949e-7_fp*q*T**3*Log(RH)**2 + &
       (5.677455521083783e-6_fp*q*T**3*Log(RH)**2)/h3 - &
       0.00003507762649727548*h6*q*T**3*Log(RH)**2 + &
       1392.9578972108989*Log(H2SO4)*Log(RH)**2 + &
       (1911.734041113941*Log(H2SO4)*Log(RH)**2)/h3 + &
       179329.36146873524*h6*Log(H2SO4)*Log(RH)**2 - &
       (119328.53674704138*h6*Log(H2SO4)*Log(RH)**2)/h3 - &
       6.462104618811676*q*Log(H2SO4)*Log(RH)**2 + &
       (11.294074621589449*q*Log(H2SO4)*Log(RH)**2)/h3 + &
       13.39033080899747*h6*q*Log(H2SO4)*Log(RH)**2 - &
       (8.52831107206456*h6*q*Log(H2SO4)*Log(RH)**2)/h3 - &
       16.8958062365253*T*Log(H2SO4)*Log(RH)**2 - &
       (24.81735534395592*T*Log(H2SO4)*Log(RH)**2)/h3 - &
       1464.394638508739*h6*T*Log(H2SO4)*Log(RH)**2 + &
       (981.4137808275404*h6*T*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.051952833287610724*q*T*Log(H2SO4)*Log(RH)**2 - &
       (0.061982178004231*q*T*Log(H2SO4)*Log(RH)**2)/h3 - &
       0.7507183329339466*h6*q*T*Log(H2SO4)*Log(RH)**2 - &
       (0.05082083499340442*h6*q*T*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.060317519684020285*T**2*Log(H2SO4)*Log(RH)**2 + &
       (0.10238074986938017*T**2*Log(H2SO4)*Log(RH)**2)/h3 + &
       2.654746963853928*h6*T**2*Log(H2SO4)*Log(RH)**2 - &
       (1.9505367684260284*h6*T**2*Log(H2SO4)*Log(RH)**2)/h3 - &
       0.00008828022840243558*q*T**2*Log(H2SO4)*Log(RH)**2 - &
       (7.843133775464917e-6_fp*q*T**2*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.00439346543661198*h6*q*T**2*Log(H2SO4)*Log(RH)**2 - &
       0.00005488176049657228*T**3*Log(H2SO4)*Log(RH)**2 - &
       (0.00012357549613231235*T**3*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.000877040667994904*h6*T**3*Log(H2SO4)*Log(RH)**2 - &
       4.394231627146249e-8_fp*q*T**3*Log(H2SO4)*Log(RH)**2 - &
       72.11406227014747*Log(H2SO4)**2*Log(RH)**2 - &
       (16.86445272085156*Log(H2SO4)**2*Log(RH)**2)/h3 - &
       6535.765336704547*h6*Log(H2SO4)**2*Log(RH)**2 + &
       (941.2813201741736*h6*Log(H2SO4)**2*Log(RH)**2)/h3 + &
       0.2428031989628447*q*Log(H2SO4)**2*Log(RH)**2 - &
       (0.3050220578412721*q*Log(H2SO4)**2*Log(RH)**2)/h3 + &
       3.3649457326647947*h6*q*Log(H2SO4)**2*Log(RH)**2 + &
       (0.9000970248172522*h6*q*Log(H2SO4)**2*Log(RH)**2)/h3 + &
       0.8503624487797692*T*Log(H2SO4)**2*Log(RH)**2 + &
       (0.2527834506442179*T*Log(H2SO4)**2*Log(RH)**2)/h3 + &
       50.73348890297342*h6*T*Log(H2SO4)**2*Log(RH)**2 - &
       (4.356161470448808*h6*T*Log(H2SO4)**2*Log(RH)**2)/h3 - &
       0.0018342607279861292*q*T*Log(H2SO4)**2*Log(RH)**2 + &
       (0.0020365272689109475*q*T*Log(H2SO4)**2*Log(RH)**2)/h3 - &
       0.03276629834613374*h6*q*T*Log(H2SO4)**2*Log(RH)**2 - &
       0.0028420227040081517*T**2*Log(H2SO4)**2*Log(RH)**2 - &
       (0.0009263932662681723*T**2*Log(H2SO4)**2*Log(RH)**2)/h3 - &
       0.08964774960777555*h6*T**2*Log(H2SO4)**2*Log(RH)**2 + &
       3.036637047548649e-6_fp*q*T**2*Log(H2SO4)**2*Log(RH)**2 + &
       2.0594310710626316e-6_fp*T**3*Log(H2SO4)**2*Log(RH)**2 + &
       0.9950338851031777*Log(H2SO4)**3*Log(RH)**2 - &
       (0.622965174549953*Log(H2SO4)**3*Log(RH)**2)/h3 + &
       48.41845482837981*h6*Log(H2SO4)**3*Log(RH)**2 - &
       (1.7738934087048919*h6*Log(H2SO4)**3*Log(RH)**2)/h3 - &
       0.0022234755298572505*q*Log(H2SO4)**3*Log(RH)**2 - &
       (0.002437957119890622*q*Log(H2SO4)**3*Log(RH)**2)/h3 + &
       0.08069009533680993*h6*q*Log(H2SO4)**3*Log(RH)**2 - &
       0.010935011617216193*T*Log(H2SO4)**3*Log(RH)**2 + &
       (0.003456846985810427*T*Log(H2SO4)**3*Log(RH)**2)/h3 - &
       0.24482068250824904*h6*T*Log(H2SO4)**3*Log(RH)**2 + &
       0.000011334503487127534*q*T*Log(H2SO4)**3*Log(RH)**2 + &
       0.000029425270779265584*T**2*Log(H2SO4)**3*Log(RH)**2

    h1=EXP(h1)

    if (h1.gt.q)then
       h1=q
    elseif (h1.lt.0.e+0_fp)then
       h1=0.e+0_fp
    endif

    h2=-32043.03148295406 + 59725.428570008815/h3 + &
       7.128537634261564e+6_fp*h6 - (1.3833467233343722e+7_fp*h6)/h3+ &
       33.63110252227136*q-(48.61215633992165*q)/h3 - 16602.414377611287 &
       *h6*q + (40754.788181739124* h6*q)/h3 -2.3397851800516185*q**2 + &
       (4.426964073992281*q**2)/h3 + 18.971418767591036*h6*q**2 + &
       (60.446718551038344*h6*q**2)/h3+396.33752593131607*T- &
       (650.8601684277011*T)/h3-83753.25337253512*h6*T+ &
       (140932.8771905448*h6*T)/h3 + 1.7779514612590905*q*T - &
       (4.1201474547289845*q*T)/h3 - 14.23100399324848*h6*q*T - &
       (199.64146093004214*h6*q*T)/h3 + 0.017586270137944494*q**2*T - &
       (0.040669111407675304*q**2*T)/h3 - 0.5190592639152767*h6*q**2*T + &
       (0.4152291779457336*h6*q**2*T)/h3 - 1.6239261240452716*T**2 + &
       (2.290370025362425*T**2)/h3 + 323.71145981054207*h6*T**2 - &
       (443.70978129586814*h6*T**2)/h3 - 0.016944490878752067*q*T**2 + &
       (0.03190748579500371*q*T**2)/h3 + 0.9921294086888162*h6*q*T**2 - &
       (0.017807643647025154*h6*q*T**2)/h3-0.000026182356192850867*q**2* &
       T**2+(0.00018413948144606674*q**2*T**2)/h3+0.001418259144142961* &
       h6*q**2*T**2 - (0.0005314099952660705*h6*q**2*T**2)/h3 + &
       0.002203350630670886*T**3-(0.002628197592847762*T**3)/h3 - &
       0.41364399538414304*h6*T**3+(0.4236315079925672*h6*T**3)/h3 + &
       0.000036034470518775974*q*T**3-(0.00005072877953687334*q*T**3)/h3 &
       -0.0025779163615301994*h6*q*T**3-(0.0001801635475946943*h6*q* &
       T**3)/h3 - 1.0397496813871261e-8_fp*q**2*T**3 - &
       (4.794510915054545e-7_fp*q**2*T**3)/h3+ &
       1.0037940606886203e-6_fp*h6* &
       q**2*T**3 +4057.4863055406972*Log(H2SO4) - (4857.4895775507075* &
       Log(H2SO4))/h3 - 873923.0975121224*h6*Log(H2SO4) + &
       (1.3884698205066123e+6_fp*h6*Log(H2SO4))/h3- &
       16.567334471680468*q*Log(H2SO4)+(39.06695157867744*q* &
       Log(H2SO4))/h3 + &
       4033.3161579197213*h6*q*Log(H2SO4) - (4837.3892389479715*h6*q* &
       Log(H2SO4))/h3 + 0.32643050327864925*q**2*Log(H2SO4) - &
       (0.29168894526478634*q**2*Log(H2SO4))/h3 +   3.0197469938934343* &
       h6*q**2*Log(H2SO4) - (13.620934035683936*h6*q**2*Log(H2SO4))/h3 - &
       50.64853826213881*T*Log(H2SO4) + (46.190674056264115*T* &
       Log(H2SO4))/h3 +10309.186386370407*h6*T*Log(H2SO4) - &
       (13393.112078195805*h6*T* &
       Log(H2SO4))/h3 - 0.05988575301198189*q*T*Log(H2SO4) + &
       (0.12991016720908877*q*T*Log(H2SO4))/h3-25.211481437032155*h6*q* &
       T*Log(H2SO4) + (25.76780209283505*h6*q*T*Log(H2SO4))/h3 - &
       0.0025433151017881413*q**2*T*Log(H2SO4)-(0.00003896417606619728* &
       q**2*T*Log(H2SO4))/h3+0.018638367565640454*h6*q**2*T*Log(H2SO4)- &
       (0.009466280650899394*h6*q**2*T*Log(H2SO4))/h3+ &
       0.20947015213812362*T**2*Log(H2SO4)-(0.1282261301435015*T**2* &
       Log(H2SO4))/h3-40.05118528139879*h6*T**2*Log(H2SO4)+ &
       (37.82992227601891* h6*T**2*Log(H2SO4))/h3+0.001426085690256897* &
       q*T**2*Log(H2SO4)-(0.0026110978166799066*q*T**2*Log(H2SO4))/h3 - &
        0.0017889111275040672*h6*q*T**2*Log(H2SO4) - &
       (0.0010479732179162718* h6*q*T**2*Log(H2SO4))/h3 + &
       4.179579131654065e-6_fp * q**2*T**2*Log(H2SO4) + &
       (7.970218003902365e-6_fp * q**2 * &
       T**2*Log(H2SO4))/h3-0.000116249941623789*h6*q**2*T**2*Log(H2SO4)- &
       0.0002870253324673844*T**3*Log(H2SO4)+(0.00009186981667688776* &
       T**3*Log(H2SO4))/h3 + 0.05147918448072947*h6*T**3*Log(H2SO4) - &
       (0.028081229965757286*h6*T**3*Log(H2SO4))/h3 - &
       3.560981861287522e-6_fp* q* &
       T**3*Log(H2SO4)+(4.961243888005889e-6_fp*q*T**3*Log(H2SO4))/h3+ &
       0.0001420993695424454*h6*q*T**3*Log(H2SO4)+ &
       1.2719812342415112e-9_fp* &
       q**2*T**3*Log(H2SO4) - 124.19545742393447*Log(H2SO4)**2 + &
       (41.54331926307027*Log(H2SO4)**2)/h3 + 26490.80092234763*h6* &
       Log(H2SO4)**2- (31462.37707525288*h6*Log(H2SO4)**2)/h3 + &
       0.8500524958176608*q* Log(H2SO4)**2 - (3.2674516205233317*q* &
       Log(H2SO4)**2)/h3 - 170.8115272895333*h6*q*Log(H2SO4)**2 + &
       (120.16514929908757*h6*q* Log(H2SO4)**2)/h3 - &
       0.011268864148405246*q**2*Log(H2SO4)**2 +(0.019897271615564007* &
       q**2*Log(H2SO4)**2)/h3 -  0.2574159414899909* &
       h6*q**2*Log(H2SO4)**2 + (0.5058681955079796*h6*q**2* &
       Log(H2SO4)**2)/h3 +1.5616031047481935*T*Log(H2SO4)**2 + &
       (0.04762580087464638*T* Log(H2SO4)**2)/h3 - 313.57384684195245* &
       h6*T*Log(H2SO4)**2+(270.50302801880633*h6*T*Log(H2SO4)**2)/h3 - &
        0.0026826956178692776* q* &
       T*Log(H2SO4)**2 + (0.02263655181544639*q*T*Log(H2SO4)**2)/h3 + &
       1.4230688700541085*h6*q*T*Log(H2SO4)**2 - (0.6649649390384627* &
        h6* q* T*Log(H2SO4)**2)/h3 + 0.00009215121101630449*q**2*T* &
       Log(H2SO4)**2 - &
       (0.00010215400234653555*q**2*T*Log(H2SO4)**2)/h3 + &
       0.0010136511879910042*h6*q**2*T*Log(H2SO4)**2 - &
       0.006504837589925847*T**2*Log(H2SO4)**2 - (0.002702979699750996* &
       T**2*Log(H2SO4)**2)/h3 + 1.223462360247441*h6*T**2*Log(H2SO4)**2- &
       (0.5623561810520002*h6*T**2*Log(H2SO4)**2)/h3 - &
       0.000024272415265506407*q*T**2*Log(H2SO4)**2 - &
       (0.0000317372505770328*q*T**2*Log(H2SO4)**2)/h3 - &
       0.0028928533476792*h6*q*T**2* &
       Log(H2SO4)**2-1.7589107450560238e-7_fp* &
       q**2*T**2*Log(H2SO4)** &
       2+8.98009092200891e-6_fp*T**3*Log(H2SO4)**2+ &
       (7.173660190761036e-6_fp*T**3*Log(H2SO4)**2)/h3 - &
       0.0015801362599946241*h6*T**3* &
       Log(H2SO4)**2+8.220297383016583e-8_fp* &
       q*T**3*Log(H2SO4)**2 - 12589.220049398413*Log(RH)+ &
       (71533.62210328173*Log(RH))/h3+ &
       1.0678104434003264e+6_fp*h6*Log(RH)- &
       (1.6102655953002474e+6_fp*h6*Log(RH))/h3+ &
       51.82639715490156*q*Log(RH)- &
       (38.71320196661908*q*Log(RH))/h3+3767.8963041701336*h6*q*Log(RH)- &
       (2025.2741146328435*h6*q*Log(RH))/h3+0.43253945689885376*q**2* &
       Log(RH) - (0.5285905724764535*q**2*Log(RH))/h3+2.800088035800146* &
       h6*q**2*Log(RH)+(5.830987444741845*h6*q**2*Log(RH))/h3+ &
       156.76303859222733*T*Log(RH) -  (1050.9788795963977*T* &
       Log(RH))/h3 + 2968.558570914762*h6*T*Log(RH)- &
       (13144.694726104603*h6*T*Log(RH))/h3-0.3939836553978611*q*T* &
       Log(RH) - (0.15149542273585964*q*T*Log(RH))/h3-52.00514647221519* &
       h6*q*T*Log(RH)+(33.50499777824741*h6*q*T*Log(RH))/h3- &
       0.0044939367750510195*q**2*T*Log(RH)+ &
       (0.0037387212699065446*q**2*T*Log(RH))/h3 - &
       0.019659867844446746*h6*q**2*T*Log(RH) - (0.013728505197733986* &
       h6*q**2*T*Log(RH))/h3-0.6541703069347068*T**2*Log(RH) + &
       (5.190500347012039*T**2*Log(RH))/h3- &
       87.05017506628032*h6*T**2*Log(RH) + (179.83498689899525*h6*T**2* &
       Log(RH))/h3 + 0.00042182566738198824*q*T**2*Log(RH) + &
       (0.004271845253967774*q*T**2*Log(RH))/h3+0.12589872273911631*h6* &
       q* T**2*Log(RH) - (0.08117178992783239*h6*q*T**2*Log(RH))/h3 + &
       0.000010942253625927065*q**2*T**2*Log(RH) - &
       (3.4136292496834907e-6_fp*q**2* &
       T**2*Log(RH))/h3+3.7003812518449353e-6_fp*h6*q**2*T**2*Log(RH) + &
       0.0009227426647911644*T**3*Log(RH) - (0.008600685587907114*T**3* &
       Log(RH))/h3 + 0.23093344962838908*h6*T**3*Log(RH) - &
       (0.38790097045592176*h6*T**3*Log(RH))/h3 + &
       1.2941739001871342e-6_fp* &
       q*T**3*Log(RH)-(0.000013459445557723903*q*T**3*Log(RH))/h3 + &
       0.00016881487161048804* h6*q*T**3*Log(RH) + &
       2.609280671396356e-9_fp* &
       q**2*T**3 *Log(RH)+2130.4592878663534*Log(H2SO4)*Log(RH) - &
       (4313.3562495734495*Log(H2SO4)*Log(RH))/h3 - &
       291848.82747478905*h6*Log(H2SO4)*Log(RH) + &
       (488820.18312982953*h6*Log(H2SO4)*Log(RH))/h3-9.405300971012823* &
       q*Log(H2SO4)*Log(RH) + (12.85263348658218*q*Log(H2SO4)* &
       Log(RH))/h3+102.49140755096978*h6*q*Log(H2SO4)*Log(RH) - &
       (289.1347848442807*h6*q*Log(H2SO4)*Log(RH))/h3 - &
       0.01500268773941549*q**2*Log(H2SO4)*Log(RH) + &
       (0.018600134825856034*q**2*Log(H2SO4)*Log(RH))/h3 - &
       0.07717355560425733* &
       h6*q**2*Log(H2SO4)*Log(RH)-(0.161938896012148*h6*q**2*Log(H2SO4)* &
       Log(RH))/h3 - 27.134659704074295*T*Log(H2SO4)*Log(RH) + &
       (69.2483692103918*T*Log(H2SO4)*Log(RH))/h3 + 1882.179099085481* &
       h6*T*Log(H2SO4)*Log(RH)-(3332.724995554323*h6*T*Log(H2SO4)* &
       Log(RH))/h3 + 0.0786126279415837*q* &
       T*Log(H2SO4)*Log(RH) - (0.12353477014379059*q*T*Log(H2SO4)* &
       Log(RH))/h3+2.221848005253306*h6*q*T*Log(H2SO4)*Log(RH) + &
       (0.532770625490143*h6 *q*T*Log(H2SO4)*Log(RH))/h3 + &
       0.00020323051093599427*q**2*T* Log(H2SO4)* &
       Log(RH) - (0.00014374306638210478*q**2*T*Log(H2SO4)*Log(RH))/h3 + &
       0.0006995787862129086*h6*q**2*T*Log(H2SO4)*Log(RH) + &
       0.11535250161710112*T**2*Log(H2SO4)*Log(RH)-(0.37254120769567123* &
       T**2*Log(H2SO4)*Log(RH))/h3 + 0.3210533828417492*h6*T**2* &
       Log(H2SO4)*Log(RH)+(4.899886058456062*h6*T**2*Log(H2SO4)* &
       Log(RH))/h3-0.00013223802083620908*q*T**2*Log(H2SO4)*Log(RH) + &
       (0.00031364246871883094*q*T**2*Log(H2SO4)*Log(RH))/h3 - &
       0.013017818218354784*h6*q*T**2*Log(H2SO4)*Log(RH) - &
       6.285230956812909e-7_fp* &
       q**2*T**2*Log(H2SO4)*Log(RH) -0.0001644893214543928*T**3* &
       Log(H2SO4)*Log(RH)+(0.0006683187277334718*T**3*Log(H2SO4)* &
       Log(RH))/h3-0.012564516901538326*h6*T**3*Log(H2SO4)*Log(RH)- &
       1.2575655284907035e-7_fp*q*T**3*Log(H2SO4)*Log(RH) - &
       79.49703918541012*Log(H2SO4)**2*Log(RH) - &
       (101.22507090731304*Log(H2SO4)**2*Log(RH))/h3 + &
       12719.306036991808*h6*Log(H2SO4)**2*Log(RH)-(7977.5339222745315* &
       h6*Log(H2SO4)**2*Log(RH))/h3 &
       +0.3744372178099974*q*Log(H2SO4)**2*Log(RH)-(0.025682479965146*q* &
       Log(H2SO4)**2*Log(RH))/h3-22.19218619809272*h6*q*Log(H2SO4)**2* &
       Log(RH) + (5.910834971937343*h6*q*Log(H2SO4)**2*Log(RH))/h3 - &
       0.00037981877987097496*q**2*Log(H2SO4)**2*Log(RH) + &
       (0.00032985668220561373*q**2*Log(H2SO4)**2*Log(RH))/h3 - &
       0.0003940881289017262*h6*q**2*Log(H2SO4)**2*Log(RH) + &
       1.0255615837685275*T*Log(H2SO4)**2*Log(RH) + &
       (1.0124457331884233*T* &
       Log(H2SO4)**2*Log(RH))/h3-111.19643387882614*h6*T*Log(H2SO4)**2* &
       Log(RH) + (37.058545749337576*h6*T*Log(H2SO4)**2*Log(RH))/h3 - &
       0.0034094873321621174*q*T*Log(H2SO4)**2*Log(RH) + &
       (0.000015597253696421123*q*T*Log(H2SO4)**2*Log(RH))/h3 + &
       0.10750304600542902*h6*q*T*Log(H2SO4)**2*Log(RH) + &
      1.8808923253076581e-6_fp* &
       q**2*T*Log(H2SO4)**2*Log(RH) - 0.0044066707400722644*T**2*  &
       Log(H2SO4)**2*Log(RH) -  (0.0025163539282379434*T**2* &
       Log(H2SO4)**2*Log(RH))/h3 + 0.24020207161395796*h6*T**2* &
       Log(H2SO4)**2*Log(RH)+7.634986243805219e-6_fp* &
       q*T**2*Log(H2SO4)**2* &
       Log(RH) + 6.328401246898811e-6_fp*T**3*Log(H2SO4)**2*Log(RH) + &
       3630.6862033225625*Log(RH)**2-(1075.4966438125716*Log(RH)**2)/h3+ &
       54546.69751557024*h6*Log(RH)**2-(49710.530231480734*h6* &
       Log(RH)**2)/h3 - 1.5893360668636096*q*Log(RH)**2 + &
       (7.970183065343727*q*Log(RH)**2)/h3 + 400.88335315886525*h6*q* &
       Log(RH)**2 -(219.0329650802132*h6*q*Log(RH)**2)/h3 + &
       0.015009992625400783*q**2*Log(RH)**2 + (0.0036804085191588995* &
       q**2*Log(RH)**2)/h3+0.4203732847045355*h6*q**2*Log(RH)**2 + &
       (0.20391771437862574*h6*q**2*Log(RH)**2)/h3 - 44.579109086853684* &
       T*Log(RH)**2 + (8.53054910480424*T* &
       Log(RH)**2)/h3 + 76.04921125760566*h6*T*Log(RH)**2 - &
       (329.28669019015933*h6*T*Log(RH)**2)/h3 - 0.019940534315381526*q &
       *T*Log(RH)**2 - (0.050640623913744444*q*T*Log(RH)**2)/h3 - &
       1.5765203107981385*h6*q*T*Log(RH)**2 +  (0.3694218342124642*h6* &
       q*T* Log(RH)**2)/h3-0.0001222138921423729*q**2*T*Log(RH)**2 - &
       (0.000049403831741944986*q**2* &
       T*Log(RH)**2)/h3 - 0.0025870767829077484*h6*q**2*T*Log(RH)**2 + &
       0.17463801499094986*T**2*Log(RH)**2 - (0.02046987546015939*T**2* &
       Log(RH)**2)/h3 - 1.3780954635968534*h6*T**2*Log(RH)**2 + &
       (2.0745682824062985*h6*T**2*Log(RH)**2)/h3 + &
       0.00015350888582331981* q*T**2*Log(RH)**2 + &
       (0.00011285744022151752*q*T**2*Log(RH)**2)/h3 + &
       0.00004952110557031714*h6*q*T**2* &
       Log(RH)**2+3.968468979853371e-7_fp* &
       q**2*T**2*Log(RH)**2 - 0.00021374165659919817*T**3*Log(RH)**2 + &
       (0.000023714077438108623*T**3*Log(RH)**2)/h3 - &
       0.0010671420677054659* h6* &
       T**3*Log(RH)**2 - 1.8562790214558605e-7_fp*q*T**3*Log(RH)**2 - &
       279.60556688841575*Log(H2SO4)*Log(RH)**2 +  (44.483779255378025* &
       Log(H2SO4)*Log(RH)**2)/h3 -11085.947063809159*h6*Log(H2SO4)* &
       Log(RH)**2 +(10910.046235786824*h6*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.47924678115681674*q*Log(H2SO4)*Log(RH)**2-(0.3759508982199312* &
       q*Log(H2SO4)*Log(RH)**2)/h3-29.071647637053548*h6*q*Log(H2SO4)* &
       Log(RH)**2 +  (8.887445070958586*h6*q*Log(H2SO4)*Log(RH)**2)/h3- &
       0.00029636660559253115*q**2 * Log(H2SO4)*Log(RH)**2 + &
       (0.00046301434118339234*q**2*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.005479820108166143*h6*q**2*Log(H2SO4)*Log(RH)**2 + &
       3.2937887074139884*T*Log(H2SO4)*Log(RH)**2 - &
       (0.053358185813312024*T*Log(H2SO4)*Log(RH)**2)/h3 + &
       29.117876039650476*h6*T*Log(H2SO4)*Log(RH)**2 - &
       (34.13451610985067*h6*T*Log(H2SO4)*Log(RH)**2)/h3 - &
       0.001150661716592889*q*T*Log(H2SO4)*Log(RH)**2 + &
       (0.0005497364286653244*q*T* Log(H2SO4)* &
       Log(RH)**2)/h3+0.09171149584822953*h6*q*T*Log(H2SO4)*Log(RH)**2- &
       2.9335804996871176e-6_fp*q**2*T*Log(H2SO4)*Log(RH)**2 - &
       0.012024558565450824*T**2*Log(H2SO4)*Log(RH)**2 - &
       (0.0008282043267863183*T**2*Log(H2SO4)*Log(RH)**2)/h3 + &
       0.11350434127907331*h6* T**2* &
       Log(H2SO4)*Log(RH)**2-3.206307709389464e-6_fp*q*T**2*Log(H2SO4)* &
       Log(RH)**2 + 0.000012855779740625904*T**3*Log(H2SO4)*Log(RH)**2 + &
       3.7577104325072703*Log(H2SO4)**2*Log(RH)**2 -(2.6641547962978636* &
       Log(H2SO4)**2*Log(RH)**2)/h3 +  476.38347866428205*h6* &
       Log(H2SO4)**2 *Log(RH)**2 - (127.7453688897455*h6*Log(H2SO4)**2* &
       Log(RH)**2)/h3 -0.01970020856731193*q*Log(H2SO4)**2*Log(RH)**2 + &
       (0.006983199758249699*q*Log(H2SO4)**2*Log(RH)**2)/h3 + &
       0.3243000572140283*h6*q*Log(H2SO4)**2* &
       Log(RH)**2+0.00002903250111944733*q**2*Log(H2SO4)**2*Log(RH)**2- &
       0.037500815576593044*T*Log(H2SO4)**2*Log(RH)**2 + &
       (0.013073560271857765*T*Log(H2SO4)**2*Log(RH)**2)/h3 - &
       2.2753830913797373*h6*T* Log(H2SO4)**2* &
       Log(RH)**2 + 0.00007918460467376976*q*T*Log(H2SO4)**2*Log(RH)**2+ &
       0.00009291148493939081*T**2*Log(H2SO4)**2*Log(RH)**2

    h2=exp(h2)

    h4=-233.3693139924163 + 3711.127600293859*h6 - &
       127375.45943800849*h6**2 - 0.6541599370168311*q - &
       8.950348936875036*h6*q + 1420.4060399615116*h6**2*q + &
       0.006010885721884837*q**2 - 0.2514391282801529*h6*q**2 + &
       11.74107168004114*h6**2*q**2 - 27.242866772851034*RH + &
       3230.6550683739456*h6*RH - 7739.349030052802*h6**2*RH + &
       0.02657310586465451*q*RH + 0.8072083676135904*h6*q*RH + &
       21.19451916114249*h6**2*q*RH - 0.0013789987709190107*q**2*RH - &
       0.02985872690339605*h6*q**2*RH + 10.1858919054768*RH**2 + &
       1831.5638235525591*h6*RH**2 - 3345.6757256829833*h6**2*RH**2 - &
       0.006667965604100408*q*RH**2 - 1.6996949068091805*h6*q*RH**2 - &
       0.0002938711589315329*q**2*RH**2 + 20.8476664473087*RH**3 + &
       153.46841353011587*h6*RH**3 - 0.023573221099724897*q*RH**3 + &
       3.1234153503175706*T - 106.11845318552218*h6*T + &
       1616.0925045006472*h6**2*T - 0.000870039433136904*q*T + &
       0.07485289069770486*h6*q*T - 2.8492194924162093*h6**2*q*T - &
       0.000025041587132070856*q**2*T - &
       0.00009872033200025474*h6*q**2*T + 0.1139999815738393*RH*T - &
       59.3235586131405*h6*RH*T + 188.66428911211182*h6**2*RH*T + &
       0.0021232663933975797*q*RH*T + 0.03525906170060478*h6*q*RH*T + &
       5.50009952133025e-6_fp*q**2*RH*T - 0.40949236112495135*RH**2*T - &
       9.19061496536003*h6*RH**2*T + &
       0.00018605332374157369*q*RH**2*T - &
       0.15228517482883766*RH**3*T - 0.026779347826859368*T**2 + &
       0.7791278560987452*h6*T**2 - 5.635512845808971*h6**2*T**2 - &
       0.00004290197829367679*q*T**2 - &
       0.00034888935617559023*h6*q*T**2 + &
       3.399576950217114e-8_fp*q**2*T**2 + &
       0.0047495543899101315*RH*T**2 + &
       0.21966980586301704*h6*RH*T**2 - &
       5.402900826149558e-6_fp*q*RH*T**2 + &
       0.00209481671640675*RH**2*T**2 + 0.00004712913787763615*T**3- &
       0.0022395082766246705*h6*T**3 + 7.857921629051029e-8_fp*q*T**3 - &
       0.000017637397545586075*RH*T**3 + &
       14.487613189612244*Log(H2SO4) + &
       834.2954132542228*h6*Log(H2SO4) - &
       9818.723117205202*h6**2*Log(H2SO4) + &
       0.1293876064755356*q*Log(H2SO4) + &
       0.5266286653232314*h6*q*Log(H2SO4) - &
       85.4813173795405*h6**2*q*Log(H2SO4) - &
       0.0005441269120155648*q**2*Log(H2SO4) + &
       0.013203398533271429*h6*q**2*Log(H2SO4) + &
       2.9312753776692975*RH*Log(H2SO4) + &
       329.1801633051309*h6*RH*Log(H2SO4) - &
       2384.7570613779067*h6**2*RH*Log(H2SO4) - &
       0.02754902616213419*q*RH*Log(H2SO4) - &
       0.39067534860398334*h6*q*RH*Log(H2SO4) + &
       0.000059670451934966084*q**2*RH*Log(H2SO4) + &
       2.1796062691153963*RH**2*Log(H2SO4) - &
       4.33856569910847*h6*RH**2*Log(H2SO4) + &
       0.0016947748838741473*q*RH**2*Log(H2SO4) + &
       0.7647868504545748*RH**3*Log(H2SO4) + &
       0.22920773390043764*T*Log(H2SO4) - &
       8.113046295543565*h6*T*Log(H2SO4) + &
       68.80414049795816*h6**2*T*Log(H2SO4) + &
       0.0014123733231002104*q*T*Log(H2SO4) + &
       0.0062552686253774335*h6*q*T*Log(H2SO4) + &
       3.138422092020368d-8*q**2*T*Log(H2SO4) - &
       0.11854147612732342*RH*T*Log(H2SO4) - &
       2.0186043272890317*h6*RH*T*Log(H2SO4) - &
       0.000017510448352402354*q*RH*T*Log(H2SO4) - &
       0.014662481263341446*RH**2*T*Log(H2SO4) + &
       0.001037948515905241*T**2*Log(H2SO4) + &
       0.04470222473271311*h6*T**2*Log(H2SO4) - &
       7.507637022804968e-7_fp*q*T**2*Log(H2SO4) + &
       0.0002763185325508037*RH*T**2*Log(H2SO4) - &
       2.046693588795596e-6_fp*T**3*Log(H2SO4) - &
       3.3271552306181786*Log(H2SO4)**2 - &
       0.6732293679771482*h6*Log(H2SO4)**2 - &
       123.7357115914075*h6**2*Log(H2SO4)**2 - &
       0.017682151597098173*q*Log(H2SO4)**2 - &
       0.06932927383915821*h6*q*Log(H2SO4)**2 + &
       0.000028013174553488913*q**2*Log(H2SO4)**2 + &
       0.5385357123752617*RH*Log(H2SO4)**2 + &
       4.405735442137794*h6*RH*Log(H2SO4)**2 + &
       0.0007699683891127761*q*RH*Log(H2SO4)**2 - &
       0.01981622019488726*RH**2*Log(H2SO4)**2 - &
       0.02863075409562691*T*Log(H2SO4)**2 - &
       0.3815368435566718*h6*T*Log(H2SO4)**2 - &
       0.000032173458991045656*q*T*Log(H2SO4)**2 + &
       0.0006610846557652917*RH*T*Log(H2SO4)**2 + &
       5.708155280776268e-6_fp*T**2*Log(H2SO4)**2 + &
       0.2841922065650612*Log(H2SO4)**3 + &
       1.937916853785423*h6*Log(H2SO4)**3 + &
       0.00046645968682534585*q*Log(H2SO4)**3 - &
       0.01392240494812053*RH*Log(H2SO4)**3 + &
       0.0005604898286672238*T*Log(H2SO4)**3 - &
       0.00648946009121241*Log(H2SO4)**4

    H4=EXP(H4)

    h5=68.64045827314231-3277.3575769882523*h6 + 1.0798559249565618*q- &
       25.296110707348316*h6*q + 13.398992645698215*RH + &
       922.4932305036297*h6*RH - 0.27140107873619296*q*RH + &
       20.08312325165439*h6*q*RH + 66.82077511984484*RH**2 + &
       1611.1977384351555*h6*RH**2 - 0.02661518788217287*q*RH**2 + &
       3.0843227537972138*h6*q*RH**2 - 1.4080258142983926*T - &
       1.8568570408634648*h6*T - 0.0037450866352058397*q*T - &
       0.2576980690602505*h6*q*T - 1.3810906781490837*RH*T - &
       17.730890154257356*h6*RH*T + 0.001076745543266636*q*RH*T - &
       0.08909158555166723*h6*q*RH*T - 0.87735596568044*RH**2*T - &
       5.603639080904061*h6*RH**2*T + 0.0002224812194282388*q*RH**2*T + &
       0.015404311250124585*T**2 + 0.009972776386592903*h6*T**2 - &
       7.895537111616416e-7_fp*q*T**2-0.0003918985727565079*h6*q*T**2 + &
       0.014090613376170114*RH*T**2 + 0.029815313352513736*h6*RH*T**2 - &
       3.607688962369381e-6_fp*q*RH*T** &
       2+0.0028248811086430334*RH**2*T**2- &
       0.00010105566309790024*T**3 - 0.00008368737713099654*h6*T**3 + &
       1.0956121068002908e-8_fp*q*T**3-0.000042201374093774165*RH*T**3+ &
       2.429718549962505e-7_fp*T**4 - 0.9150072698367441*Log(H2SO4) + &
       653.1422380210507*h6*Log(H2SO4) - 0.16169777329907886*q* &
       Log(H2SO4)+6.860485142132651*h6*q*Log(H2SO4)+13.912279451959678* &
       RH*Log(H2SO4) - 12.99843208521885*h6*RH*Log(H2SO4) + &
       0.021566098556935497*q*RH*Log(H2SO4) - 0.17455209117051437*h6*q* &
       RH*Log(H2SO4) + 3.4003096116229066*RH**2*Log(H2SO4) - &
       36.48723431979357*h6*RH**2*Log(H2SO4) - 0.001787998461490877*q* &
       RH**2*Log(H2SO4) - 0.11123959974265946*T*Log(H2SO4) + &
       1.4641341095178422*h6*T*Log(H2SO4) + 0.0005023439404053182*q*T* &
       Log(H2SO4) + 0.031171816066622334*h6*q*T*Log(H2SO4) - &
       0.1666443655044969*RH*T*Log(H2SO4) + 0.9274525071251363*h6*RH*T* &
       Log(H2SO4) + 0.00003323858217484082*q*RH*T*Log(H2SO4) - &
       0.021039758605777333*RH**2*T*Log(H2SO4) + 0.0018827939246916288* &
       T**2*Log(H2SO4) + 0.0014616406839659635*h6*T**2*Log(H2SO4) - &
       2.8887052861955783e-7_fp*q*T**2*Log(H2SO4) + &
       0.0007036899839828474* &
       RH*T**2*Log(H2SO4) - 6.077113667352327e-6_fp*T**3*Log(H2SO4) + &
       0.8758449734328851*Log(H2SO4)**2 - 61.130559299849466*h6* &
       Log(H2SO4)**2 + 0.007311434341744119*q*Log(H2SO4)**2 - &
       0.4613762652439876*h6*q*Log(H2SO4)**2 + 0.11048996127121676*RH* &
       Log(H2SO4)**2 - 5.241218416739315*h6*RH*Log(H2SO4)**2 - &
       0.0009849834685920686*q*RH*Log(H2SO4)**2 + 0.034417912409269905* &
       RH**2*Log(H2SO4)**2 - 0.01737690919778605*T*Log(H2SO4)**2 - &
       0.10868884345327394*h6*T*Log(H2SO4)**2-0.000013815430650770167* &
       q*T*Log(H2SO4)**2 - 0.0038903960432723123*RH*T*Log(H2SO4)**2 + &
       0.00005879779401379263*T**2*Log(H2SO4)**2 + 0.03393085910585771* &
       Log(H2SO4)**3 + 2.2221405215918755*h6*Log(H2SO4)**3 - &
       0.00009520679794957912*q*Log(H2SO4)**3 + 0.0166181059051036*RH* &
       Log(H2SO4)**3 - 0.00014027668903237805*T*Log(H2SO4)**3

    if(h5.lt.0.2)then
       h5=0.5
    elseif(h5.gt.5.)then
       h5=5.
    endif

    !if (h1 .gt. 1.d5) then ! take care of weird nuc rate blow up
    !   h1=0.e+0_fp
    !   h2=0.e+0_fp
    !   h3=0.e+0_fp
    !   h4=0.e+0_fp
    !   h5=0.e+0_fp
    !   h6=0.e+0_fp
    !   return
    !endif

    return

  End subroutine ion_nucl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tomas
!
! !DESCRIPTION: Subroutine CLEANUP/_TOMAS deallocates all module arrays
!  (win, 7/9/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_TOMAS
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! CLEANUP_TOMAS begins here!
    !=================================================================
    IF ( ALLOCATED( MOLWT       ) ) DEALLOCATE( MOLWT      )
    IF ( ALLOCATED( H2SO4_RATE  ) ) DEALLOCATE( H2SO4_RATE )

  END SUBROUTINE CLEANUP_TOMAS
!EOC
END MODULE TOMAS_MOD
#endif
