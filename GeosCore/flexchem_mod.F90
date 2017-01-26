!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: flexchem_mod.F90
!
! !DESCRIPTION: Module FlexChem\_Mod contines arrays and routines for the
!  FlexChem chemical solver.
!\\
!\\
! !INTERFACE: 
!
MODULE FlexChem_Mod
!
! !USES:
!
  USE Precision_Mod            ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Do_FlexChem
  PUBLIC :: Init_FlexChem
  PUBLIC :: Cleanup_FlexChem
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  15 Jun 2016 - M. Sulprizio- Remove STTTOCSPEC mapping array. Species and
!                              tracers have a 1:1 mapping currently so mapping
!                              is not required
!  18 Jul 2016 - M. Sulprizio- Remove FAMILIES_KLUDGE routine. Family tracers
!                              have been eliminated.
!  24 Aug 2016 - M. Sulprizio- Rename from flexchem_setup_mod.F90 to
!                              flexchem_mod.F90
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Save the H value for the Rosenbrock solver
  REAL(fp),  ALLOCATABLE :: HSAVE_KPP(:,:,:) 

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_flexchem
!
! !DESCRIPTION: Subroutine Do\_FlexChem is the driver subroutine for
!  full chemistry with KPP.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_FlexChem( am_I_Root, Input_Opt, State_Met, State_Chm, RC  )
!
! !USES:
!
    USE AEROSOL_MOD,          ONLY : SOILDUST, AEROSOL_CONC, RDAER
    USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_CHEMGRID
    USE CMN_FJX_MOD
    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
    USE DIAG_MOD,             ONLY : AD65
    USE DIAG_OH_MOD,          ONLY : DO_DIAG_OH
    USE DIAG20_MOD,           ONLY : DIAG20, POx, LOx
    USE DUST_MOD,             ONLY : RDUST_ONLINE, RDUST_OFFLINE
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FAST_JX_MOD,          ONLY : PHOTRATE_ADJ, FAST_JX
    USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
    USE GCKPP_HetRates,       ONLY : SET_HET
    USE GCKPP_Monitor,        ONLY : SPC_NAMES, FAM_NAMES
    USE GCKPP_Parameters
    USE GCKPP_Integrator,     ONLY : INTEGRATE, NHnew
    USE GCKPP_Function 
    USE GCKPP_Model
    USE GCKPP_Global
    USE GCKPP_Rates,          ONLY : UPDATE_RCONST, RCONST
    USE GCKPP_Initialize,     ONLY : Init_KPP => Initialize
    USE GC_GRID_MOD,          ONLY : GET_YMID
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE PRESSURE_MOD        
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    USE TIME_MOD,             ONLY : GET_MONTH
    USE TIME_MOD,             ONLY : GET_YEAR
    USE UnitConv_Mod
#if defined( UCX )
    USE UCX_MOD,              ONLY : CALC_STRAT_AER
    USE UCX_MOD,              ONLY : SO4_PHOTFRAC
    USE UCX_MOD,              ONLY : UCX_NOX
    USE UCX_MOD,              ONLY : UCX_H2SO4PHOT
#endif
#if   defined( TOMAS )
    USE TOMAS_MOD,            ONLY : H2SO4_RATE
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  18 Dec 2015 - M. Sulprizio- Add calls to OHSAVE and DO_DIAG_OH
!  22 Dec 2015 - M. Sulprizio- Set NUMDEN to State_Met%AIRNUMDEN for conversions
!                              between v/v and molec/cm3
!  28 Mar 2016 - R. Yantosca - Added several fixes for OpenMP parallelization
!  30 Mar 2016 - R. Yantosca - Bug fix, now make sure to copy C back into
!                              State_Chm%Species.  Also block out temp diags.
!  16 May 2016 - M. Sulprizio- Remove call to RURALBOX. The implemementation of
!                              FlexChem has rendered the routine obsolete.
!  31 May 2016 - E. Lundgren - Use species database MW instead of XNUMOL
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  14 Jun 2016 - M. Sulprizio- Replace loops over N_TRACERS with loops over
!                              State_Chm%nSpecies and add checks for Is_Advected
!                              and Is_Kpp to avoid introducing numerical noise
!                              by applying unit conversions to non-KPP species
!  16 Jun 2016 - M. Sulprizio- Now define IDTCH4 locally
!  20 Jun 2016 - R. Yantosca - Renamed IDTCH4 to id_CH4 for consistency
!  18 Jul 2016 - M. Sulprizio- Remove calls to FAMILIES_KLUDGE. Family tracers
!                              have been eliminated. Also simplify code to copy
!                              to/from STT so that unit conversions to/from
!                              molec/cm3 are done in the same step as the copy.
!  02 Aug 2016 - M. Sulprizio- Connect production and loss rates from KPP to
!                              ND65 diagnostic 
!  16 Aug 2016 - E. Lundgren - Remove all references to tracers, including
!                              STT and Input_Opt%N_TRACERS, and use routines in
!                              unitconv_mod.F for species kg <-> molec/cm3 
!  24 Aug 2016 - M. Sulprizio- Replace CSPECTOKPP with State_Chm%Map_KppSpc
!  24 Aug 2016 - M. Sulprizio- Move this subroutine to flexchem_mod.F90 and
!                              rename from FLEX_CHEMDR to Do_FlexChem
!  22 Sep 2016 - R. Yantosca - Add extra debug printout after FAST_JX
!  14 Nov 2016 - E. Lundgren - Move UCX calls to after spc conversion to kg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N, F, M, SpcID
    INTEGER                :: MONTH, YEAR
    INTEGER                :: WAVELENGTH
    INTEGER                :: HCRC
    LOGICAL                :: prtDebug
    LOGICAL                :: LDUST
    CHARACTER(LEN=155)     :: DiagnName

    ! SAVEd variables
    LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
    INTEGER, SAVE          :: CH4_YEAR  = -1
    INTEGER, SAVE          :: id_CH4    = -1

    ! For grid-box latitude
    REAL(fp)               :: YLAT

    ! For interannually-varying Methane
    REAL(fp), SAVE         :: C3090S, C0030S, C0030N, C3090N

    ! Global rate constants
    REAL(dp)               :: GLOB_RCONST(IIPAR,JJPAR,LLPAR,NREACT)

    ! KPP variables
    INTEGER                :: IERR
    REAL(fp)               :: T, TIN, TOUT
    REAL(dp)               :: RCNTRL(20)
    REAL(dp)               :: RSTATE(20)
    INTEGER                :: ICNTRL(20)
    INTEGER                :: ISTATUS(20)

    ! For computing statistics from the integrator
    REAL(fp)               :: Start, Finish, rtim, itim
    INTEGER                :: TotSteps, TotFuncs, TotJacob
    INTEGER                :: TotAccep, TotRejec, TotNumLU
    CHARACTER(LEN=255)     :: ERR_MSG

    ! For passing SO4_PHOTFRAC output to PHOTRATE_ADJ (bmy, 3/28/16)
    REAL(fp)               :: SO4_FRAC

    ! To avoid array temporaries in call to SET_HET (bmy, 3/28/16)
    REAL(fp)               :: SCOEFF(IIPAR,JJPAR,LLPAR,3)
    REAL(fp)               :: SCF(3)

    ! For prod/loss diagnostic
    REAL(fp)               :: FAM(NFAM)
    INTEGER                :: IND
    INTEGER                :: COEF
    CHARACTER(LEN=14)      :: NAME

    ! For testing only, may be removed later (mps, 4/26/16)
    LOGICAL                :: DO_HETCHEM
    LOGICAL                :: DO_PHOTCHEM

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! Do_FlexChem begins here!
    !=================================================================

    ! Initialize pointer
    SpcInfo => NULL()

    ! Turn heterogeneous chemistry and photolysis on/off here
    ! This is for testing only and may be removed later (mps, 4/26/16)
    DO_HETCHEM  = .TRUE.
    DO_PHOTCHEM = .TRUE.

    IF ( FIRSTCHEM ) THEN
       WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       WRITE( 6, '(a,l,a)' ) '# FLEX_CHEMDR: DO_HETCHEM  =', &
                                             DO_HETCHEM,  ' #'
       WRITE( 6, '(a,l,a)' ) '# FLEX_CHEMDR: DO_PHOTCHEM =', &
                                             DO_PHOTCHEM, ' #'
       WRITE( 6, '(a)' ) REPEAT( '#', 32 )
    ENDIF

    ! Initialize variables
    prtDebug  = ( Input_Opt%LPRT .and. am_I_Root )
    LDUST     = Input_Opt%LDUST

    itim      = 0.0_fp
    rtim      = 0.0_fp
    totsteps  = 0
    totfuncs  = 0
    totjacob  = 0
    totaccep  = 0
    totrejec  = 0
    totnumLU  = 0

    ! Get month and year
    MONTH     = GET_MONTH()
    YEAR      = GET_YEAR()

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
    !
    ! Do nothing, since we will call these setup routines from the 
    ! init method of the ESMF interface (bmy, 10/24/12)
    !-----------------------------------------------------------------
#else
    !-----------------------------------------------------------------
    !         %%%%%%% GEOS-Chem CLASSIC (with OpenMP) %%%%%%%
    !
    ! Call the following chemistry setup routines only on the very
    ! first chemistry timestep.  This is current practice in the
    ! std GEOS-Chem. (bmy, 10/24/12)
    !-----------------------------------------------------------------
    IF ( FIRSTCHEM ) THEN

       !---------------------------------
       ! Set global concentration of CH4
       !---------------------------------
       ! Define advected species ID flag
       id_CH4 = Ind_('CH4','A')

       ! Check that CH4 is not an advected species
       ! Check that CH4 is a KPP species (ind_CH4 is from gckpp_Monitor.F90)
       IF ( id_CH4 <= 0 .and. ind_CH4 > 0 .and. ( CH4_YEAR /= YEAR ) ) THEN

          ! If CH4 is a species, then call GET_GLOBAL_CH4
          ! to return the globally-varying CH4 conc. as a function of
          ! year and latitude bin. (bnd, bmy, 7/1/03)
          !
          ! If we are using the future emissions, then get the CH4
          ! concentrations for FUTURE_YEAR.  Otherwise get the CH4
          ! concentrations for the current met field year. 
          ! (swu, havala, bmy, 1/24/08)
          IF ( Input_Opt%LFUTURE ) THEN
             CH4_YEAR = GET_FUTURE_YEAR()
          ELSE
             CH4_YEAR = YEAR
          ENDIF

          ! Get CH4 [ppbv] in 4 latitude bins for each year
          CALL GET_GLOBAL_CH4( CH4_YEAR,  .TRUE.,  C3090S, &
                               C0030S,    C0030N,  C3090N, & 
                               am_I_Root, Input_Opt        )
       ENDIF

    ENDIF
#endif

    !================================================================
    ! Get concentrations of aerosols in [kg/m3] 
    ! for FAST-JX and optical depth diagnostics
    !=================================================================
    IF ( Input_Opt%LSULF .or. Input_Opt%LCARB .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

#if defined( UCX )
       ! Calculate stratospheric aerosol properties (SDE 04/18/13)
       CALL CALC_STRAT_AER( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### FLEX_CHEMDR: after CALC_PSC' )
       ENDIF
#endif

       ! Compute aerosol concentrations
       CALL AEROSOL_CONC( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

    ENDIF

    !================================================================
    ! Zero out certain species:
    !    - isoprene oxidation counter species (dkh, bmy, 6/1/06)
    !    - isoprene-NO3 oxidation counter species (hotp, 6/1/10)
    !    - if SOA or SOA_SVPOA, aromatic oxidation counter species 
    !      (dkh, 10/06/06)
    !    - if SOA_SVPOA, LNRO2H and LNRO2N for NAP (hotp 6/25/09
    !=================================================================
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! isoprene oxidation counter species
       IF ( TRIM( SpcInfo%Name ) == 'LISOPOH' .or. &
            TRIM( SpcInfo%Name ) == 'LISOPNO3' ) THEN
          State_Chm%Species(:,:,:,N) = 0e+0_fp
       ENDIF

       ! aromatic oxidation counter species
       IF ( Input_Opt%LSOA .or. Input_Opt%LSVPOA ) THEN
          SELECT CASE ( TRIM( SpcInfo%Name ) )
             CASE ( 'LBRO2H', 'LBRO2N', 'LTRO2H', 'LTRO2N', &
                    'LXRO2H', 'LXRO2N', 'LNRO2H', 'LNRO2N' )
                State_Chm%Species(:,:,:,N) = 0e+0_fp
          END SELECT
       ENDIF

       ! Temporary fix for CO2
       ! CO2 is a dead species and needs to be set to zero to 
       ! match the old SMVGEAR code (mps, 6/14/16)
       IF ( TRIM( SpcInfo%Name ) == 'CO2' ) THEN
          State_Chm%Species(:,:,:,N) = 0e+0_fp
       ENDIF

       ! Free pointer
       SpcInfo => NULL()

    ENDDO
      
    !=================================================================
    ! Call RDAER -- computes aerosol optical depths
    !=================================================================

    ! Call RDAER to compute AOD for FAST-JX (skim, 02/03/11)
    WAVELENGTH = 0
    CALL RDAER( am_I_Root, Input_Opt, State_Met,  State_Chm, &
                RC,        MONTH,     YEAR,       WAVELENGTH )

    !### Debug
    IF ( prtDebug ) THEN 
       CALL DEBUG_MSG( '### FLEX_CHEMDR: after RDAER' )
    ENDIF

    !=================================================================
    ! If LDUST is turned on, then we have online dust aerosol in
    ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
    ! compute aerosol optical depth for FAST-JX, etc.
    !
    ! If LDUST is turned off, then we do not have online dust aerosol
    ! in GEOS-CHEM...so read monthly-mean dust files from disk.
    ! (rjp, tdf, bmy, 4/1/04)
    !=================================================================
    IF ( LDUST ) THEN
       CALL RDUST_ONLINE( am_I_Root, Input_Opt,  State_Met, State_Chm, &
                          SOILDUST,  WAVELENGTH, RC )
    ELSE
#if !defined( TOMAS )
       ! Don't read dust emissions from disk when using TOMAS,
       ! because TOMAS uses a different set of dust species than the 
       ! std code (win, bmy, 1/25/10)
       CALL RDUST_OFFLINE( am_I_Root, Input_Opt, State_Met, State_Chm, &
                           MONTH,     YEAR,      WAVELENGTH, RC )
#endif
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### FLEX_CHEMDR: after RDUST' )
    ENDIF

    !================================================================
    ! Convert species from [kg] to [molec/cm3] (ewl, 8/16/16)
    !================================================================
    CALL ConvertSpc_Kg_to_MND( am_I_Root, State_Met, State_Chm, RC )
      
    !=================================================================
    ! Call photolysis routine to compute J-Values
    !=================================================================
    CALL FAST_JX( WAVELENGTH, am_I_Root,  Input_Opt, &
                  State_Met,  State_Chm,  RC         )

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### FLEX_CHEMDR: after FAST_JX' )
    ENDIF

    !=================================================================
    ! Set up integration convergence conditions and timesteps
    ! (cf. M. J. Evans)
    !
    ! NOTE: ATOL and RTOL are defined in gckpp_Global.F90 so they
    ! are probably only used as INTENT(IN).  Therefore, it is
    ! probably safe to define them here outside the OpenMP loop.
    ! (bmy, 3/28/16)
    !
    ! The ICNTRL vector specifies options for the solver.  We can
    ! define ICNTRL outside of the "SOLVE CHEMISTRY" parallel loop 
    ! below because it is passed to KPP as INTENT(IN). 
    ! (bmy, 3/28/16)
    !=================================================================

    !%%%%% TIMESTEPS %%%%%

    ! mje Set up conditions for the integration
    ! mje chemical timestep and convert it to seconds.
    DT        = GET_TS_CHEM() * 60d0
    T         = 0d0
    TIN       = T
    TOUT      = T + DT

    !%%%%% CONVERGENCE CRITERIA %%%%%

    ! Absolute tolerance
    ATOL      = 1e-2_dp    

#if defined( UCX )
    ! Relative tolerance, for UCX-based mechanisms
    RTOL      = 2e-2_dp    
#else
    ! Relative tolerance, non-UCX mechanisms
    RTOL      = 2e-1_dp    
#endif

    !%%%%% SOLVER OPTIONS %%%%%

    ! Zero all slots of ICNTRL
    ICNTRL    = 0
      
    ! 0 - non-autonomous, 1 - autonomous
    ICNTRL(1) = 1

    ! 0 - vector tolerances, 1 - scalars	
    ICNTRL(2) = 0	
      
    ! Select Integrator
    ! ICNTRL(3)  -> selection of a particular method.
    ! For Rosenbrock, options are:
    ! = 0 :  default method is Rodas3
    ! = 1 :  method is  Ros2
    ! = 2 :  method is  Ros3 
    ! = 3 :  method is  Ros4 
    ! = 4 :  method is  Rodas3
    ! = 5:   method is  Rodas4
    ICNTRL(3) = 4

    ! 0 - adjoint, 1 - no adjoint
    ICNTRL(7) = 1

    !=================================================================
    ! %%%%% SOLVE CHEMISTRY -- This is the main KPP solver loop %%%%%
    !=================================================================
100 format('No. of function calls:', i6, /,                        &
           'No. of jacobian calls:', i6, /,                        &
           'No. of steps:         ', i6, /,                        &
           'No. of accepted steps:', i6, /,                        &
           'No. of rejected steps ', i6, /,                        &
           '       (except at very beginning)',          /,        &
           'No. of LU decompositions:             ', i6, /,        &
           'No. of forward/backward substitutions:', i6, /,        &
           'No. of singular matrix decompositions:', i6, /,        &
            /,                                                     &
           'Texit, the time corresponding to the      ',        /, &
           '       computed Y upon return:            ', f11.4, /, &
           'Hexit, last accepted step before exit:    ', f11.4, /, &
           'Hnew, last predicted step (not yet taken):', f11.4 )

    !-----------------------------------------------------------------
    ! NOTE: The following variables are held THREADPRIVATE and 
    ! therefore do not need to be included in the !$OMP+PRIVATE 
    ! statements below: C, VAR, FIX, RCONST, TIME, TEMP, NUMDEN, 
    ! H2O, PRESS, PHOTOL, HET, and CFACTOR. (bmy, 3/28/16)
    !-----------------------------------------------------------------
    !$OMP PARALLEL DO                                           &
    !$OMP DEFAULT  ( SHARED                                   ) &
    !$OMP PRIVATE  ( I,     J,        L,          N,     YLAT ) &
    !$OMP PRIVATE  ( SCF,   SO4_FRAC, IERR,       RCNTRL      ) &
    !$OMP PRIVATE  ( START, FINISH,   ISTATUS,    RSTATE      ) &
    !$OMP PRIVATE  ( FAM,   SpcId                             ) &
    !$OMP PRIVATE  ( NAME,  COEF,     IND,        F,     M    ) &
    !$OMP REDUCTION( +:ITIM                                   ) &
    !$OMP REDUCTION( +:RTIM                                   ) &
    !$OMP REDUCTION( +:TOTSTEPS                               ) &
    !$OMP REDUCTION( +:TOTFUNCS                               ) &
    !$OMP REDUCTION( +:TOTJACOB                               ) &
    !$OMP REDUCTION( +:TOTACCEP                               ) &
    !$OMP REDUCTION( +:TOTREJEC                               ) &
    !$OMP REDUCTION( +:TOTNUMLU                               ) &
    !$OMP SCHEDULE ( DYNAMIC, 1                               )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       !==============================================================
       ! For safety's sake, zero variables for each grid box (I,J,L),
       ! whether or not chemistry will be done there.
       !==============================================================
       HET      = 0.0_dp    ! Het chem array
       H2O      = 0.0_dp    ! H2O concentration
       IERR     = 0         ! Success or failure flag
       ISTATUS  = 0.0_dp    ! Rosenbrock output 
       PHOTOL   = 0.0_dp    ! Photolysis array
       PRESS    = 0.0_dp    ! Pressure
       NUMDEN   = 0.0_dp    ! Air number density
       RCNTRL   = 0.0_fp    ! Rosenbrock input
       RSTATE   = 0.0_dp    ! Rosenbrock output
       SCF      = 0.0_dp    ! For hetchem diagnostics
       SO4_FRAC = 0.0_fp    ! Fraction of SO4 available for photolysis
       TEMP     = 0.0_fp    ! Temperature
       YLAT     = 0.0_fp    ! Latitude

       !==============================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !==============================================================

       ! If we are not in the troposphere don't do the chemistry!
       IF ( .not. ITS_IN_THE_CHEMGRID(I,J,L,State_Met) ) CYCLE

       ! Skipping buffer zone (lzh, 08/10/2014)
       IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
          IF ( J <=         Input_Opt%NESTED_J0W ) CYCLE
          IF ( J >  JJPAR - Input_Opt%NESTED_J0E ) CYCLE
          IF ( I <=         Input_Opt%NESTED_I0W ) CYCLE
          IF ( I >  IIPAR - Input_Opt%NESTED_I0E ) CYCLE
       ENDIF

       !==============================================================
       ! Initialize quantities for boxes where chemistry happens
       !==============================================================

       ! Grid-box latitude [degrees]
       YLAT   = GET_YMID( I, J, L )

       ! Temperature [K]
       TEMP   = State_Met%T(I,J,L)

       ! Pressure [hPa]
       PRESS  = GET_PCENTER( I, J, L )

       ! mje Calculate NUMDEN based on ideal gas law (# cm-3)
       NUMDEN = State_Met%AIRNUMDEN(I,J,L)

       ! mje H2O arrives in g/kg needs to be in mol cm-3 
       H2O    = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)

       ! Intialize CFACTOR, VAR, and FIX arrays
       CALL Init_KPP( )

       !===========================================================
       ! Get rates for heterogeneous chemistry
       !===========================================================

!#if defined( DEVEL )
!       ! Get starting time for rate computation
!       CALL CPU_TIME( start )
!#endif

       IF ( DO_HETCHEM ) THEN

          ! Copy SCOEFF to SCF, this will avoid an array temporary
          SCF = SCOEFF(I,J,L,:)

          ! Set hetchem rates
          CALL SET_HET( I, J, L, State_Chm, State_Met, Input_Opt, SCF )

          ! Copy SCF back into SCOEFF
          SCOEFF(I,J,L,:) = SCF

       ENDIF

       !===========================================================
       ! Get photolysis rates (daytime only)
       !
       ! NOTE: The ordering of the photolysis reactions here is
       ! the order in the Fast-J definition file FJX_j2j.dat.
       ! I've assumed that these are the same as in the text files
       ! but this may have been changed.  This needs to be checked
       ! through more thoroughly -- M. Long (3/28/16)
       !===========================================================
       IF ( State_Met%SUNCOSmid(I,J) > 0.e+0_fp ) THEN

          ! Get the fraction of H2SO4 that is available for photolysis
          ! (this is only valid for UCX-enabled mechanisms)
#if defined( UCX )
          SO4_FRAC = SO4_PHOTFRAC( I, J, L )
#else
          SO4_FRAC = 0.0_fp
#endif

          ! Adjust certain photolysis rates:
          ! (1) H2SO4 + hv -> SO2 + OH + OH   (UCX-based mechanisms)
          ! (2) O3    + hv -> O2  + O         (UCX-based mechanisms)
          ! (2) O3    + hv -> OH  + OH        (trop-only mechanisms)
          CALL PHOTRATE_ADJ( am_I_root, I, J, L, NUMDEN, TEMP, &
                             H2O, SO4_FRAC, IERR  )

          IF ( DO_PHOTCHEM ) THEN

             ! Copy ZPJ from fast_jx_mod.F to PHOTOL array
             PHOTOL(1:JVN_) = ZPJ(L,1:JVN_,I,J)

          ENDIF

       ENDIF

       !===========================================================
       ! Initialize species concentrations
       !===========================================================
       ! Loop over KPP Species
       DO N=1,NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Initialize KPP species concentration array
          IF ( SpcID .eq. 0) THEN
             C(N) = 0d0
          ELSE
             C(N) = State_Chm%Species(I,J,L,SpcID)
          ENDIF

       ENDDO

       C(ind_ACTA)  = 0.0_dp
       C(ind_EOH)   = 0.0_dp
       C(ind_HCOOH) = 0.0_dp
#if defined( UCX )
       C(ind_O2) = 0.2095e+0_dp * NUMDEN
       C(ind_N2) = 0.7808e+0_dp * NUMDEN
       C(ind_H2) = 0.5000e-6_dp * NUMDEN
#else
       ! Need to copy H2O to the C array for KPP (mps, 4/25/16)
       ! NOTE: H2O is a tracer in UCX and is obtained from State_Chem%Species
       C(ind_H2O) = H2O
#endif

       ! id_CH4 is the GEOS-Chem advected species index for CH4
       ! ind_CH4 is the KPP species index for CH4
       ! If both are nonzero, then CH4 is an advected species.
       ! If id_CH4 <= 0 but ind_CH4 > 0, then CH4 is a non-advected species.
       ! (bmy, 6/20/16_)
       IF ( id_CH4 <= 0 .and. ind_CH4 > 0 ) THEN
          ! Set CH4 according to latitude
          ! Convert from [ppbv CH4] to [molec CH4/cm3]
          IF ( YLAT < -30.0_fp ) THEN
             C(ind_CH4) = C3090S * 1e-9_dp * NUMDEN
          ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
             C(ind_CH4) = C0030S * 1e-9_dp * NUMDEN
          ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
             C(ind_CH4) = C0030N * 1e-9_dp * NUMDEN
          ELSE
             C(ind_CH4) = C3090N * 1e-9_dp * NUMDEN
          ENDIF
       ENDIF

       !===========================================================
       ! Update KPP's rates
       !===========================================================

       ! VAR and FIX are chunks of array C (mps, 2/24/16)
       !
       ! NOTE: Because VAR and FIX are !$OMP THREADPRIVATE, they
       ! cannot appear in an EQUIVALENCE statement.  Therfore, we
       ! will just copy the relevant elements of C to VAR and FIX
       ! here. (bmy, 3/28/16)
       VAR(1:NVAR) = C(1:NVAR)
       FIX         = C(NVAR+1:NSPEC)

       ! Update the array of rate constants
       CALL Update_RCONST( )

!#if defined( DEVEL )
!       ! Get time when rate computation finished
!       CALL CPU_TIME( finish )
!
!       ! Compute how long it took for KPP to compute rates
!       rtim = rtim + finish - start
!#endif

       !===========================================================
       ! Set options for the KPP Integrator (M. J. Evans)
       !
       ! NOTE: Because RCNTRL(3) is set to an array value that
       ! depends on (I,J,L), we must declare RCNTRL as PRIVATE
       ! within the OpenMP parallel loop and define it just 
       ! before the call to to Integrate. (bmy, 3/24/16)
       !===========================================================

       ! Zero all slots of RCNTRL
       RCNTRL    = 0.0_fp

       ! Starting value for integration time step
       RCNTRL(3) = HSAVE_KPP(I,J,L)

       !===========================================================
       ! Integrate the box forwards
       !===========================================================

!#if defined( DEVEL )
!         ! Get time before integrator starts
!         CALL CPU_TIME( start )
!#endif

       ! Call the KPP integrator
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
                       RCNTRL, ISTATUS, RSTATE, IERR )

       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
       ENDIF

!#if defined( DEVEL )
!       ! Get time when integrator ends
!       CALL CPU_TIME( finish )
!
!       ! Compute how long the integrator took 
!       itim     = itim + finish - start
!
!       ! Compute other statistics from the integrator
!       totfuncs = totfuncs + ISTATUS(1)
!       totjacob = totjacob + ISTATUS(2)
!       totsteps = totsteps + ISTATUS(3)
!       totaccep = totaccep + ISTATUS(4)
!       totrejec = totrejec + ISTATUS(5)
!       totnumLU = totnumLU + ISTATUS(6)
!#endif

       ! Zero certain species
       C(ind_ACTA)  = 0.e0_dp
       C(ind_EOH)   = 0.e0_dp
       C(ind_HCOOH) = 0.e0_dp

       ! Try another time if it failed
       IF ( IERR < 0 ) THEN

          ! Reset first time step and start concentrations
          ! Retry the integration with non-optimized
          ! settings
          RCNTRL(3)  = 0e+0_fp
          CALL Init_KPP( )
          VAR = C(1:NVAR)
          FIX = C(NVAR+1:NSPEC)
          CALL Update_RCONST( )
          CALL Integrate( TIN,    TOUT,    ICNTRL,      &
                          RCNTRL, ISTATUS, RSTATE, IERR )
          IF ( IERR < 0 ) THEN 
             WRITE(6,*) '## INTEGRATE FAILED TWICE !!! '
             WRITE(ERR_MSG,'(a,i3)') 'Integrator error code :',IERR
             CALL ERROR_STOP(ERR_MSG, 'INTEGRATE_KPP')
          ENDIF
            
       ENDIF

       ! Copy VAR and FIX back into C (mps, 2/24/16)
       C(1:NVAR)       = VAR(:)
       C(NVAR+1:NSPEC) = FIX(:)

       ! Save for next integration time step
       HSAVE_KPP(I,J,L) = RSTATE(Nhnew)

       ! Save rate constants in global array
       GLOB_RCONST(I,J,L,:) = RCONST(:)

       !==============================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species
       !==============================================================
       ! Loop over KPP species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Skip if this is not a GEOS-Chem species
          IF ( SpcID .eq. 0 ) CYCLE

          ! Set negative concentrations to zero
          C(N) = MAX( C(N), 0.0E0_dp )

          ! Copy concentrations back into State_Chm%Species
          State_Chm%Species(I,J,L,SpcID) = REAL( C(N), kind=fp )

       ENDDO

       !==============================================================
       ! Prod/loss diagnostic
       !==============================================================
       IF ( Input_Opt%DO_SAVE_PL ) THEN

          ! Obtain prod/loss rates from KPP [molec/cm3]
          CALL ComputeFamilies( VAR, FAM )

          ! Loop over # prod/loss families
          DO F = 1, NFAM

             !--------------------------------------------------------
             ! Add to AD65 array [molec/cm3/s]
             !--------------------------------------------------------
             AD65(I,J,L,F) = AD65(I,J,L,F) + FAM(F) / DT

             !--------------------------------------------------------
             ! Save out P(Ox) and L(Ox) from the fullchem simulation
             ! for a future tagged O3 run
             !--------------------------------------------------------
             IF ( Input_Opt%DO_SAVE_O3 ) THEN
                IF ( TRIM(FAM_NAMES(F)) == 'POx' ) THEN
                   POx(I,J,L) = FAM(F) / DT
                ENDIF
                IF ( TRIM(FAM_NAMES(F)) == 'LOx' ) THEN
                   LOx(I,J,L) = FAM(F) / DT
                ENDIF
             ENDIF

#if defined( TOMAS )
             !-------------------------------------------------------
             ! FOR TOMAS MICROPHYSICS:
             !
             ! Obtain P/L with a unit [kg S] for tracing 
             ! gas-phase sulfur species production (SO2, SO4, MSA)
             ! (win, 8/4/09)
             !-------------------------------------------------------

             ! Calculate H2SO4 production rate [kg s-1] in each
             ! time step (win, 8/4/09)
             IF ( TRIM(FAM_NAMES(F)) == 'PSO4' ) THEN 
                H2SO4_RATE(I,J,L) = FAM(F) / AVO * 98.e-3_fp * & ! Hard-coded MW
                                    State_Met%AIRVOL(I,J,L)  * &
                                    1e+6_fp / DT
             ENDIF
#endif
          ENDDO

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

#if defined( DEVEL )
    write(*,'(a,F10.3)') 'Flex Rate Time     : ', rtim
    write(*,'(a,F10.3)') 'Flex Intg Time     : ', itim
    write(*,'(a,I9)'   ) 'Flex Function Calls: ', totfuncs
    write(*,'(a,I9)'   ) 'Flex Jacobian Calls: ', totjacob
    write(*,'(a,I9)'   ) 'Flex Total Steps   : ', totsteps
    write(*,'(a,I9)'   ) 'Flex Rejected Steps: ', totrejec
    write(*,'(a,I9)'   ) 'Flex LU Decompos.  : ', totnumLU
#endif

    !=================================================================
    ! Call OHSAVE which saves info on OH AND HO2 concentrations
    !=================================================================
    CALL OHSAVE( State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### FLEX_CHEMDR: after OHSAVE' )
    ENDIF

    !=================================================================
    ! Save quantities for computing mean OH lifetime
    !=================================================================
    CALL DO_DIAG_OH( State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### FLEX_CHEMDR: after DO_DIAG_OH' )
    ENDIF

    !=================================================================
    ! Save out P(O3) and L(O3) for a tagged O3 run
    !=================================================================
    IF ( Input_Opt%DO_SAVE_O3 ) THEN
       CALL DIAG20( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### FLEX_CHEMDR: after DIAG20' )
       ENDIF
    ENDIF

    !================================================================
    ! Convert species from [molec/cm3] to [kg] (ewl, 8/16/16)
    !================================================================
    CALL ConvertSpc_MND_to_Kg( am_I_Root, State_Met, State_Chm, RC )

#if defined( UCX )
    ! If using stratospheric chemistry, applying high-altitude
    ! active nitrogen partitioning and H2SO4 photolysis
    ! approximations  outside the chemgrid
    CALL UCX_NOX( Input_Opt, State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### CHEMDR: after UCX_NOX' )
    ENDIF

    CALL UCX_H2SO4PHOT( Input_Opt, State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### CHEMDR: after UCX_H2SO4PHOT' )
    ENDIF
#endif

    ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
    FIRSTCHEM = .FALSE.

  END SUBROUTINE Do_FlexChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_flexchem
!
! !DESCRIPTION: Subroutine Init\_FlexChem is used to allocate arrays for the
!  KPP solver. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE Init_FlexChem( am_I_Root, Input_Opt, RC)  
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE gckpp_Global,       ONLY : NREACT
    USE gckpp_Monitor,      ONLY : EQN_NAMES
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!    
    LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  22 Dec 2015 - M. Sulprizio- Use State_Met%AIRNUMDEN to convert initial
!                              species concentrations from v/v to molec/cm3
!  29 Jan 2016 - M. Sulprizio- Add calls to Register_Tracer and Register_Species
!                              to populate Tracer_Name, Tracer_Id, Species_Name,
!                              and Species_ID fields in State_Chm
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  06 Jun 2016 - M. Sulprizio- Remove calls to Register_Tracer and
!                              Register_Species; these routines were made
!                              obsolete by the species database
!  14 Jun 2016 - M. Sulprizio- Replace Spc_GetIndx with Ind_ (M. Long)
!  25 Jul 2016 - E. Lundgren - Add check that species was not in restart file
!                              prior to v/v -> molec/cm3 conversion
!  02 Aug 2016 - E. Lundgren - Move unit conversion of species background
!                              values to restart_mod
!  24 Aug 2016 - M. Sulprizio- Remove CSPECTOKPP array. State_Chm%Map_KppSpc is
!                              now used instead.
!  20 Sep 2016 - R. Yantosca - Use fixed integer with in WRITE statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: D
    LOGICAL :: prtDebug

    !=================================================================
    ! Init_FlexChem begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Initialize variables
    prtDebug  = Input_Opt%LPRT

    IF ( prtDebug .and. am_I_Root ) WRITE( 6, 100 )
100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

    IF ( prtDebug .and. am_I_Root ) THEN
       write(*,'(a)') ' KPP Reaction Reference '
       DO D = 1,NREACT
          WRITE( 6, '(i8,a3,a85)' ) D,' | ',EQN_NAMES(D)
       END DO
    ENDIF

  END SUBROUTINE Init_FlexChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_flexchem
!
! !DESCRIPTION: Subroutine Cleanup\_FlexChem deallocate module variables.
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE Cleanup_FlexChem()  
!    
! !REVISION HISTORY:
!  24 Aug 2016 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Cleanup_FlexChem begins here!
    !=================================================================
    IF ( ALLOCATED( HSAVE_KPP ) ) DEALLOCATE( HSAVE_KPP )

  END SUBROUTINE Cleanup_FlexChem
!EOC
END MODULE FlexChem_Mod

