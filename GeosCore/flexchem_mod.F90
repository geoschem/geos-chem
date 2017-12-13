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
  PUBLIC  :: Do_FlexChem
  PUBLIC  :: Init_FlexChem
  PUBLIC  :: Cleanup_FlexChem
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Diag_OH_HO2_O1D_O3P
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
!  17 Nov 2017 - R. Yantosca - Now call Diag_OH_HO2_O1D_O3P, which will let
!                              us remove arrays in CMN_O3_SIZE_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags (and logicals to denote if species are present)
  INTEGER                :: id_OH, id_HO2, id_O3P, id_O1D, id_CH4
  LOGICAL                :: ok_OH, ok_HO2, ok_O1D, ok_O3P

  ! Diagnostic flags
  LOGICAL                :: Do_Diag_OH_HO2_O1D_O3P
  LOGICAL                :: Do_ND43
  LOGICAL                :: Archive_OHconcAfterchem
  LOGICAL                :: Archive_HO2concAfterchem
  LOGICAL                :: Archive_O1DconcAfterchem
  LOGICAL                :: Archive_O3PconcAfterchem
  LOGICAL                :: Archive_Prod
  LOGICAL                :: Archive_Loss
  
  ! Diagnostic indices for ND65 bpch diagnostic
  INTEGER,   ALLOCATABLE :: ND65_KPP_Id(:)

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
  SUBROUTINE Do_FlexChem( am_I_Root, Input_Opt,  State_Met,  &
                          State_Chm, State_Diag, RC         )
!
! !USES:
!
    USE AEROSOL_MOD,          ONLY : SOILDUST, AEROSOL_CONC, RDAER
    USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_CHEMGRID
    USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_STRAT
    USE CMN_DIAG_MOD,         ONLY : ND52
    USE CMN_FJX_MOD
    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
    USE DIAG_MOD,             ONLY : AD65, AD52
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
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE Strat_Chem_Mod,       ONLY : SChem_Tend
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    USE TIME_MOD,             ONLY : GET_MONTH
    USE TIME_MOD,             ONLY : GET_YEAR
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
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
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
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
!  10 Mar 2017 - C. Keller   - Make sure ind_CH4 is correctly specified in
!                              ESMF environment.
!  30 May 2017 - M. Sulprizio- Add code for stratospheric chemical tendency
!                              for computing STE in strat_chem_mod.F90
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!  03 Oct 2017 - E. Lundgren - Pass State_Diag as argument
!  03 Nov 2017 - R. Yantosca - Pass State_Diag to SET_HET
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N, NA, F, SpcID, KppID
    INTEGER                :: MONTH, YEAR
    INTEGER                :: WAVELENGTH
    INTEGER                :: HCRC
    LOGICAL                :: prtDebug
    LOGICAL                :: LDUST
    CHARACTER(LEN=155)     :: DiagnName
    CHARACTER(LEN=63)      :: OrigUnit

    ! SAVEd variables
    LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
    INTEGER, SAVE          :: CH4_YEAR  = -1

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

    ! For passing SO4_PHOTFRAC output to PHOTRATE_ADJ (bmy, 3/28/16)
    REAL(fp)               :: SO4_FRAC

    ! To avoid array temporaries in call to SET_HET (bmy, 3/28/16)
    REAL(fp)               :: SCOEFF(IIPAR,JJPAR,LLPAR,3)
    REAL(fp)               :: SCF(3)

    ! For testing only, may be removed later (mps, 4/26/16)
    LOGICAL                :: DO_HETCHEM
    LOGICAL                :: DO_PHOTCHEM

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! Initial mass [kg] for computing stratospheric chemical tendency
    REAL(fp)               :: Before(IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! Do_FlexChem begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_FlexChem (in module GeosCore/flexchem_mod.F)'
    SpcInfo => NULL()

    ! Turn heterogeneous chemistry and photolysis on/off here
    ! This is for testing only and may be removed later (mps, 4/26/16)
    DO_HETCHEM  = .TRUE.
    DO_PHOTCHEM = .TRUE.

    IF ( FIRSTCHEM .AND. am_I_Root ) THEN
       WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       WRITE( 6, '(a,l,a)' ) '# Do_FlexChem: DO_HETCHEM  =', &
                                             DO_HETCHEM,  ' #'
       WRITE( 6, '(a,l,a)' ) '# Do_FlexChem: DO_PHOTCHEM =', &
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

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( Archive_Loss ) State_Diag%Loss = 0.0_f4
    IF ( Archive_Prod ) State_Diag%Prod = 0.0_f4

    !---------------------------------------
    ! Initialize global concentration of CH4
    !---------------------------------------
    IF ( FIRSTCHEM ) THEN

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

    !================================================================
    ! Get concentrations of aerosols in [kg/m3] 
    ! for FAST-JX and optical depth diagnostics
    !=================================================================
    IF ( Input_Opt%LSULF .or. Input_Opt%LCARB .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

#if defined( UCX )
       ! Calculate stratospheric aerosol properties (SDE 04/18/13)
       CALL CALC_STRAT_AER( am_I_Root, Input_Opt, State_Met, State_Chm, RC  )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FlexChem: after CALC_PSC' )
       ENDIF
#endif

       ! Compute aerosol concentrations
       CALL AEROSOL_CONC( am_I_Root, Input_Opt,  State_Met,                  &
                          State_Chm, State_Diag, RC                         )

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
    CALL RDAER( am_I_Root, Input_Opt,  State_Met,  &
                State_Chm, State_Diag, RC,         &
                MONTH,     YEAR,       WAVELENGTH )

    !### Debug
    IF ( prtDebug ) THEN 
       CALL DEBUG_MSG( '### Do_FlexChem: after RDAER' )
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
       CALL RDUST_ONLINE( am_I_Root,  Input_Opt, State_Met,  State_Chm,      &
                          State_Diag, SOILDUST,  WAVELENGTH, RC             )
    ELSE
#if !defined( TOMAS )
       ! Don't read dust emissions from disk when using TOMAS,
       ! because TOMAS uses a different set of dust species than the 
       ! std code (win, bmy, 1/25/10)
       CALL RDUST_OFFLINE( am_I_Root,  Input_Opt, State_Met, State_Chm,      &
                           State_Diag, MONTH,     YEAR,      WAVELENGTH, RC )
#endif
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after RDUST' )
    ENDIF

#if defined( UCX )
    !================================================================
    ! Archive initial species mass for stratospheric tendency
    !================================================================
    ! Loop over advected species
    DO NA = 1, State_Chm%nAdvect

       ! Get the species ID from the advected species ID
       N = State_Chm%Map_Advect(NA)

       ! Archive initial mass [kg]
       Before(:,:,:,N) = State_Chm%Species(:,:,:,N)

    ENDDO
#endif

    !================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, & 
                            State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90')
       RETURN
    ENDIF 
      
    !=================================================================
    ! Call photolysis routine to compute J-Values
    !=================================================================
    CALL FAST_JX( WAVELENGTH, am_I_Root,  Input_Opt, &
                  State_Met,  State_Chm,  State_Diag, RC )

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after FAST_JX' )
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
    !RTOL      = 2e-2_dp    
    !RTOL      = 1e-2_dp    
    RTOL      = 0.5e-2_dp    
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
    !$OMP PRIVATE  ( SpcID, KppID,    F                       ) &
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
          CALL SET_HET( I,         J,         L,         State_Diag,  &
                        State_Chm, State_Met, Input_Opt, SCF         )

#if defined( BPCH_DIAG )
          IF ( ND52 > 0 ) THEN
             ! Archive gamma values
             AD52(I,J,L,1) = AD52(I,J,L,1) + HET(ind_HO2,   1)
             AD52(I,J,L,2) = AD52(I,J,L,2) + HET(ind_IEPOXA,1) &
                                           + HET(ind_IEPOXB,1) &
                                           + HET(ind_IEPOXD,1)
             AD52(I,J,L,3) = AD52(I,J,L,3) + HET(ind_IMAE,  1)
             AD52(I,J,L,4) = AD52(I,J,L,4) + HET(ind_ISOPND,1) &
                                           + HET(ind_ISOPNB,1)
             AD52(I,J,L,5) = AD52(I,J,L,5) + HET(ind_DHDN,  1)
             AD52(I,J,L,6) = AD52(I,J,L,6) + HET(ind_GLYX,  1)
          ENDIF
#endif

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
                             H2O, SO4_FRAC, State_Diag, IERR  )

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
             C(N) = 0.0_dp
          ELSE
             C(N) = State_Chm%Species(I,J,L,SpcID)
          ENDIF

       ENDDO

       IF ( Input_Opt%DO_SAVE_PL ) THEN

          ! Loop over # prod/loss families
          DO F = 1, NFAM

             ! Determine dummy species index in KPP
             KppID = Ind_(TRIM(FAM_NAMES(F)),'K')

             ! Initialize prod/loss rates
             IF ( KppID > 0 ) C(KppID) = 0.0_dp

          ENDDO

       ENDIF

       C(ind_ACTA)  = 0.0_dp
       C(ind_HCOOH) = 0.0_dp
#if defined( UCX )
       C(ind_O2) = 0.2095e+0_dp * NUMDEN
       C(ind_N2) = 0.7808e+0_dp * NUMDEN
       C(ind_H2) = 0.5000e-6_dp * NUMDEN
#else
       ! Need to copy H2O to the C array for KPP (mps, 4/25/16)
       ! NOTE: H2O is a tracer in UCX and is obtained from State_Chm%Species
       C(ind_H2O) = H2O
#endif

       ! id_CH4 is the GEOS-Chem advected species index for CH4
       ! ind_CH4 is the KPP species index for CH4
       ! If both are nonzero, then CH4 is an advected species.
       ! If id_CH4 <= 0 but ind_CH4 > 0, then CH4 is a non-advected species.
       ! (bmy, 6/20/16_)
#if !defined(ESMF_)
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
#endif

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
             WRITE(ERRMSG,'(a,i3)') 'Integrator error code :',IERR
             CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
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

#if defined( BPCH_DIAG )
       !==============================================================
       ! ND65 (bpch) diagnostic
       !
       ! Obtain prod/loss rates from KPP [molec/cm3]
       !==============================================================
       IF ( Input_Opt%DO_SAVE_PL ) THEN

          ! Loop over # prod/loss species
          DO F = 1, nFam

             ! NOTE: KppId is the KPP ID # for each of the prod and loss
             ! diagnostic species.  This is the value used to index the
             ! KPP "VAR" array (in module gckpp_Global.F90).
             KppId         = ND65_Kpp_Id(F)

             ! Archive prod or loss for species or families [molec/cm3/s]
             AD65(I,J,L,F) = AD65(I,J,L,F) + VAR(KppID) / DT

             ! Save out P(Ox) and L(Ox) from the fullchem simulation
             ! for a future tagged O3 run
             ! NOTE: Probably not needed for netCDF diagnostics
             IF ( Input_Opt%DO_SAVE_O3 ) THEN
                IF ( TRIM(FAM_NAMES(F)) == 'POx' ) THEN
                   POx(I,J,L) = VAR(KppID) / DT
                ENDIF
                IF ( TRIM(FAM_NAMES(F)) == 'LOx' ) THEN
                   LOx(I,J,L) = VAR(KppID) / DT
                ENDIF
             ENDIF

          ENDDO
       ENDIF

#if defined( TOMAS )
       !always calculate rate for TOMAS
       DO F = 1, NFAM

          ! Determine dummy species index in KPP
          KppID =  ND65_Kpp_Id(F)

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
             ! Hard-coded MW
             H2SO4_RATE(I,J,L) = VAR(KppID) / AVO * 98.e-3_fp * &
                                 State_Met%AIRVOL(I,J,L)  * &
                                 1e+6_fp / DT

            IF ( H2SO4_RATE(I,J,L) < 0.0d0) THEN
              write(*,*) "H2SO4_RATE negative in flexchem_mod.F90!!", &
                 I, J, L, "was:", H2SO4_RATE(I,J,L), "  setting to 0.0d0"
              H2SO4_RATE(I,J,L) = 0.0d0
            ENDIF
          ENDIF
       ENDDO
#endif
#endif

#if defined( NC_DIAG )
       !==============================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "VAR" array (in module gckpp_Global.F90).
       !==============================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( Archive_Loss ) THEN
          DO F = 1, State_Chm%nLoss
             KppID                    = State_Chm%Map_Loss(F)
             State_Diag%Loss(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( Archive_Prod ) THEN
          DO F = 1, State_Chm%nProd
             KppID                    = State_Chm%Map_Prod(F)
             State_Diag%Prod(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF
#endif


    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

!!!#if defined( DEVEL )
!    write(*,'(a,F10.3)') 'Flex Rate Time     : ', rtim
!    write(*,'(a,F10.3)') 'Flex Intg Time     : ', itim
!    write(*,'(a,I9)'   ) 'Flex Function Calls: ', totfuncs
!    write(*,'(a,I9)'   ) 'Flex Jacobian Calls: ', totjacob
!    write(*,'(a,I9)'   ) 'Flex Total Steps   : ', totsteps
!    write(*,'(a,I9)'   ) 'Flex Rejected Steps: ', totrejec
!    write(*,'(a,I9)'   ) 'Flex LU Decompos.  : ', totnumLU
!!!#endif

    !=================================================================
    ! Archive OH, HO2, O1D, O3P concentrations after FlexChem solver
    ! (NOTE: This was formerly done by ohsave.F and diagoh.F)
    !=================================================================
    IF ( Do_Diag_OH_HO2_O1D_O3P ) THEN

       ! Save OH, HO2, O1D, O3P for the ND43 diagnostic
       ! NOTE: These might not be needed for netCDF, as they will already
       ! have been archived in State_Chm%Species output.
       CALL Diag_OH_HO2_O1D_O3P( am_I_Root, Input_Opt,  State_Met,           &
                                 State_Chm, State_Diag, RC                  )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Diag_OH_HO2_O1D_O3P"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FlexChem: after OHSAVE' )
       ENDIF
    ENDIF

    !=================================================================
    ! Save quantities for computing mean OH lifetime
    !=================================================================
    CALL DO_DIAG_OH( State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after DO_DIAG_OH' )
    ENDIF

    !=================================================================
    ! Save out P(O3) and L(O3) for a tagged O3 run
    !=================================================================
    IF ( Input_Opt%DO_SAVE_O3 ) THEN
       CALL DIAG20( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FlexChem: after DIAG20' )
       ENDIF
    ENDIF

    !================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90' )
       RETURN
    ENDIF  

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

    !================================================================
    ! Compute stratospheric chemical tendency for UCX simulations
    !================================================================
    !$OMP PARALLEL DO &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, NA )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) THEN

          ! Loop over advected species
          DO NA = 1, State_Chm%nAdvect

             ! Get the species ID from the advected species ID
             N = State_Chm%Map_Advect(NA)

             ! Aggregate stratospheric chemical tendency [kg box-1]
             ! for tropchem simulations
             SChem_Tend(I,J,L,N) = SChem_Tend(I,J,L,N) + &
                 ( State_Chm%Species(I,J,L,N) - Before(I,J,L,N) )

          ENDDO

       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
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
! !IROUTINE: Diag_OH_HO2_O1D_O3P
!
! !DESCRIPTION: Archives the chemical production of OH, HO2, O1D, O3P.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diag_OH_HO2_O1D_O3P( am_I_Root, Input_Opt,  State_Met,          &
                                  State_Chm, State_Diag, RC                 )
!
! !USES:
!
    USE ChemGrid_Mod,   ONLY : Its_In_The_NoChemGrid
    USE CMN_SIZE_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    !-----------------------------------------------------------
    ! NOTE: These modules are only needed for BPCH diagnostics
    ! which are slated to be removed in the near future
    USE Diag_Mod,       ONLY : AD43, LTOH, LTHO2, LTO1D, LTO3P
    !-----------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  This routine replaces both OHSAVE and DIAGOH.  Those routines were needed
!  for SMVGEAR, when we had separate arrays for the non-advected species.
!  But now, all species are stored in State_Chm%SPECIES, so the various
!  arrays (SAVEOH, SAVEHO2, etc.) are no longer necessary.  We can now just
!  just get values directly from State_Chm%SPECIES.
!
!  Also note: remove multiplication by LTOH etc arrays.  These are almost
!  always set between 0 and 24.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scalars
    LOGICAL            :: Do_Diag
    INTEGER            :: I,       J,       L

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    REAL(fp), POINTER  :: AirNumDen(:,:,:  )
    REAL(fp), POINTER  :: Spc      (:,:,:,:)
    
    !=======================================================================
    ! Diag_OH_HO2_O1D_O3P begins here!
    !=======================================================================

    ! Point to the array of species concentrations
    AirNumDen => State_Met%AirNumDen
    Spc       => State_Chm%Species

!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J, L ) &
!$OMP SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) THEN

#if defined( NC_DIAG )
            ! Zero the netCDF diagnostic arrays (if activated) above the 
            ! tropopause or mesopause to avoid having leftover values
            ! from previous timesteps
            IF ( Archive_OHconcAfterChem ) THEN
               State_Diag%OHconcAfterChem(I,J,L)  = 0.0_f4
            ENDIF

            IF ( Archive_HO2concAfterChem ) THEN
               State_Diag%HO2concAfterChem(I,J,L) = 0.0_f4
            ENDIF
            
#if defined( UCX ) 
            IF ( Archive_O1DconcAfterChem ) THEN
               State_Diag%O1DconcAfterChem(I,J,L) = 0.0_f4
            ENDIF

            IF ( Archive_O3PconcAfterChem ) THEN
               State_Diag%O3PconcAfterChem(I,J,L) = 0.0_f4
            ENDIF
#endif
#endif           
            ! Skip to next grid box
            CYCLE
         ENDIF

         !------------------------------------------------------------------
         ! OH concentration [molec/cm3]
         !------------------------------------------------------------------
         IF ( ok_OH ) THEN

#if defined( BPCH_DIAG )
            ! ND43 (bpch) diagnostic
            IF ( Do_ND43 ) THEN
               AD43(I,J,L,1) = AD43(I,J,L,1)                                 &
                             + ( Spc(I,J,L,id_OH) * LTOH(I,J)              )
            ENDIF
#endif

#if defined( NC_DIAG )
            ! HISTORY (aka netCDF diagnostics)
            IF ( Archive_OHconcAfterChem ) THEN
               State_Diag%OHconcAfterChem(I,J,L) = ( Spc(I,J,L,id_OH)        &
                                                 *   LTOH(I,J)             )
            ENDIF
#endif


         ENDIF

         !------------------------------------------------------------------
         ! HO2 concentration [v/v] 
         !------------------------------------------------------------------
         IF ( ok_HO2 ) THEN

#if defined( BPCH_DIAG )
            ! ND43 (bpch) diagnostic
            IF ( Do_ND43 ) THEN
               AD43(I,J,L,2) = AD43(I,J,L,2)                                 &
                             + ( Spc(I,J,L,id_HO2) / AirNumDen(I,J,L) )      &
                             * ( LTHO2(I,J)                           )
            ENDIF
#endif

#if defined( NC_DIAG ) 
            ! HISTORY (aka netCDF diagnostics)
            IF ( Archive_HO2concAfterChem ) THEN
               State_Diag%HO2concAfterChem(I,J,L) = ( Spc(I,J,L,id_HO2)      &
                                                  /   AirNumDen(I,J,L)     ) &
                                                  * ( LTHO2(I,J)           )
            ENDIF
#endif
         ENDIF


#if defined( UCX )

         !------------------------------------------------------------------
         ! O1D concentration [molec/cm3]
         !------------------------------------------------------------------
         IF ( ok_O1D ) THEN

#if defined( BPCH_DIAG ) 
            ! ND43 (bpch) diagnostic
            IF ( Do_ND43 ) THEN
               AD43(I,J,L,3) = AD43(I,J,L,3)                                 &
                             + ( Spc(I,J,L,id_O1D) * LTO1D(I,J)            )
            ENDIF
#endif

#if defined( NC_DIAG )
            ! HISTORY (aka netCDFdiagnostics)
            IF ( Archive_O1DconcAfterChem ) THEN
               State_Diag%O1DconcAfterChem(I,J,L) = ( Spc(I,J,L,id_O1D)      &
                                                  *   LTO1D(I,J)           )
            ENDIF
#endif
         ENDIF


         !------------------------------------------------------------------
         ! O3P concentration [molec/cm3]
         !------------------------------------------------------------------
         IF ( ok_O3P ) THEN

#if defined( BPCH_DIAG )
            ! ND43 (bpch) diagnostic
            IF ( Do_ND43 ) THEN
               AD43(I,J,L,4) = AD43(I,J,L,4)                                 &
                             + ( Spc(I,J,L,id_O3P) * LTO3P(I,J)            )
            ENDIF
#endif

#if defined( NC_DIAG )
            ! HISTORY (aka netCDF diagnostics)
            IF ( Archive_O3PconcAfterChem ) THEN
               State_Diag%O3PconcAfterChem(I,J,L) = ( Spc(I,J,L,id_O3P)      &
                                                  *   LTO3P(I,J)           )
            ENDIF
#endif

         ENDIF
#endif

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      AirNumDen => NULL()
      Spc       => NULL()

  END SUBROUTINE Diag_OH_HO2_O1D_O3P
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
  SUBROUTINE Init_FlexChem( am_I_Root, Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Gckpp_Global,     ONLY : nReact
    USE Gckpp_Monitor,    ONLY : Eqn_Names, Fam_Names
    USE Gckpp_Parameters, ONLY : nFam
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Chm_Mod,    ONLY : Ind_
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS:
!    
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Diagnostics State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: KppId,    N
    LOGICAL            :: prtDebug

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,   ThisLoc

    !=================================================================
    ! Init_FlexChem begins here!
    !=================================================================

    !--------------------------------
    ! Initialize variables
    !-------------------------------
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_FlexChem (in module GeosCore/flexchem_mod.F90)'
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )

    IF ( prtDebug ) WRITE( 6, 100 )
100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    CALL GC_CheckVar( 'flexchem_mod.F90:HSAVE_KPP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

    IF ( prtDebug ) THEN
       write( 6 ,'(a)' ) ' KPP Reaction Reference '
       DO N = 1, NREACT
          WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
       END DO
    ENDIF

    !-------------------------------
    ! Initialize species flags
    !-------------------------------

    ! Initialize species flags
    id_CH4                   = Ind_( 'CH4', 'A'     ) ! CH4 advected species
    id_HO2                   = Ind_( 'HO2'          )
    id_O3P                   = Ind_( 'O'            )
    id_O1D                   = Ind_( 'O1D'          )
    id_OH                    = Ind_( 'OH'           ) 

    ! Set flags to denote if each species is defined
    ok_HO2                   = ( id_HO2 > 0         )
    ok_O1D                   = ( id_O1D > 0         )
    ok_O3P                   = ( id_O3P > 0         )
    ok_OH                    = ( id_OH  > 0         )
    
    !-------------------------------
    ! Initialize diagnostics flags
    !-------------------------------

    ! Is the ND43 bpch diagnostic turned on?
    Do_ND43                  = ( Input_Opt%ND43 > 0 ) 

    ! Are the relevant netCDF diagnostics turned on?
    Archive_OHconcAfterChem  = ASSOCIATED( State_Diag%OHconcAfterChem  )
    Archive_HO2concAfterChem = ASSOCIATED( State_Diag%HO2concAfterChem )
    Archive_O1DconcAfterChem = ASSOCIATED( State_Diag%O1DconcAfterChem )
    Archive_O3PconcAfterChem = ASSOCIATED( State_Diag%O3PconcAfterChem )
    Archive_Loss             = ASSOCIATED( State_Diag%Loss             )
    Archive_Prod             = ASSOCIATED( State_Diag%Prod             )

    ! Should we archive OH, HO2, O1D, O3P diagnostics?
    Do_Diag_OH_HO2_O1D_O3P   = ( Do_ND43                  .or.               &
                                 Archive_OHconcAfterChem  .or.               &
                                 Archive_HO2concAfterChem .or.               &
                                 Archive_O1DconcAfterChem .or.               &
                                 Archive_O3PconcAfterChem                   )
    
    ! Pre-store the KPP indices for each KPP prod/loss species or family
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. nFam > 0 ) THEN
             
       ! Allocate mapping array for KPP Id's for ND65 bpch diagnostic
       ALLOCATE( ND65_Kpp_Id( nFam ), STAT=RC )
       CALL GC_CheckVar( 'flexchem_mod.F90:ND65_Kpp_Id', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Loop over all KPP prod/loss species
       DO N = 1, nFam

          ! NOTE: KppId is the KPP ID # for each of the prod and loss
          ! diagnostic species.  This is the value used to index the
          ! KPP "VAR" array (in module gckpp_Global.F90).
          KppID = Ind_( TRIM ( Fam_Names(N) ), 'K' )

          ! Exit if an invalid ID is encountered
          IF ( KppId <= 0 ) THEN
             ErrMsg = 'Invalid KPP ID for prod/loss species: '            // &
                       TRIM( Fam_Names(N) )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! If the species ID is OK, save in ND65_Kpp_Id
          ND65_Kpp_Id(N) = KppId
       ENDDO
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
