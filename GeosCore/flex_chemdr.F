!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: flex_chemdr
!
! !DESCRIPTION: Subroutine FLEX\_CHEMDR is the driver subroutine for
!  full chemistry with KPP.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FLEX_CHEMDR( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC  )
!
! !USES:
!
      USE AEROSOL_MOD,          ONLY : SOILDUST, AEROSOL_CONC, RDAER
      USE COMODE_MOD,           ONLY : JLOP, TAREA, ERADIUS ! <--!!
      USE COMODE_LOOP_MOD,      ONLY : BK, NAMEGAS, CONSVAP!, SMAL2
      USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_CHEMGRID
      USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_FJX_MOD
      USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
      USE DAO_MOD
      USE DIAG_OH_MOD,          ONLY : DO_DIAG_OH
      USE DUST_MOD
      USE ERROR_MOD,            ONLY : DEBUG_MSG
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE FAST_JX_MOD
      USE FLEXCHEM_SETUP_MOD
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_YEAR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GIGC_State_Chm_Mod,   ONLY : ChmState, Get_Indx
      USE GIGC_State_Met_Mod,   ONLY : MetState
      USE GRID_MOD,             ONLY : GET_YMID
      USE HCO_DIAGN_MOD,        ONLY : Diagn_Update
      USE HCO_ERROR_MOD,        ONLY : HCO_SUCCESS
      USE PRECISION_MOD
      USE PRESSURE_MOD        
      USE TIME_MOD,             ONLY : GET_TS_CHEM
      USE TIME_MOD,             ONLY : GET_MONTH
      USE TIME_MOD,             ONLY : GET_YEAR
      USE TRACERID_MOD,         ONLY : IDTCH4
#if defined( UCX )
      USE UCX_MOD,              ONLY : CALC_STRAT_AER
      USE UCX_MOD,              ONLY : SO4_PHOTFRAC
      USE UCX_MOD,              ONLY : UCX_NOX
      USE UCX_MOD,              ONLY : UCX_H2SO4PHOT
#endif

      ! mje
      USE GCKPP_HETRATES,       ONLY : SET_HET
      USE GCKPP_MONITOR,        ONLY : SPC_NAMES
      USE GCKPP_PARAMETERS
      USE GCKPP_INTEGRATOR,     ONLY : INTEGRATE, NHnew
      USE GCKPP_FUNCTION 
      USE GCKPP_MODEL
      USE GCKPP_GLOBAL
      USE GCKPP_RATES,          ONLY : UPDATE_RCONST, RCONST
      USE GCKPP_INITIALIZE,     ONLY : Init_KPP => Initialize

      IMPLICIT NONE
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, L, N, JLOOP, MONTH, YEAR
      INTEGER                :: WAVELENGTH
      INTEGER                :: HCRC
      INTEGER                :: N_TRACERS
      LOGICAL                :: prtDebug
      LOGICAL                :: LDUST
      CHARACTER(LEN=155)     :: DiagnName

      ! SAVEd variables
      LOGICAL, SAVE          :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE          :: CH4_YEAR  = -1

      ! For grid-box latitude
      REAL(fp)               :: YLAT

      ! For interannually-varying Methane
      REAL(fp),SAVE          :: C3090S, C0030S, C0030N, C3090N

      ! mje Variables for holding rate diagnoastics
      ! mje DERIV: Concentrations derivatives for one box d[X]/dt
      ! mje GLOB_DERIV: Store global values of concentrations devivatives
      ! mje GLOB_RATES: Stores global values of fluxes through each reaction
      REAL(kind=dp)          :: DERIV(NSPEC)
      REAL(kind=dp)          :: GLOB_DERIV(IIPAR,JJPAR,LLPAR,NSPEC)
      REAL(kind=dp)          :: GLOB_RATES(IIPAR,JJPAR,LLPAR,NREACT)
      REAL(kind=dp)          :: GLOB_RCONST(IIPAR,JJPAR,LLPAR,NREACT)

      ! mje KPP variables
      INTEGER                :: IERR
      REAL(kind=fp)          :: T, TIN, TOUT
      REAL(kind=dp)          :: RCNTRL(20)
      REAL(kind=dp)          :: RSTATE(20)
      INTEGER                :: ICNTRL(20)
      INTEGER                :: ISTATUS(20)

      ! mje Array to hold the molecular mass of the species
      REAL(kind=fp)          :: TCVV(Input_Opt%N_TRACERS)

      ! mje STT holds the concentrations of Tracers
      REAL(kind=fp), POINTER :: STT(:,:,:,:)

      ! mje Array for bug checking
      REAL(kind=fp)          :: AverageVAR(NSPEC)
      REAL(kind=fp)          :: Mass, Start, Finish, rtim, itim
      INTEGER                :: TotSteps, TotFuncs, TotJacob
      INTEGER                :: TotAccep, TotRejec, TotNumLU

      CHARACTER(LEN=255)     :: ERR_MSG

      ! To avoid array temporaries in call to SET_HET (bmy, 3/28/16)
      REAL(kind=fp)          :: SCOEFF(IIPAR,JJPAR,LLPAR,3)
      REAL(kind=fp)          :: SCF(3)

      !=================================================================
      ! FLEX_CHEMDR begins here!
      !=================================================================

      ! Initialize variables
      prtDebug  = Input_Opt%LPRT
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

      ! mje Copy Number of tracers (advected species) into N_TRACERS
      N_TRACERS = Input_Opt%N_TRACERS

      ! mje Copy the molecular mass of the Tracers into TCVV
      TCVV      = Input_Opt%TCVV

      ! mje link the STT Pointer to the State_Chm%Tracers stucture
      STT       => State_Chm%Tracers

      ! Locate the tropopsheric boxes
      CALL RURALBOX( State_Met )

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
         ! Check that CH4 is not a defined tracer
         IF ( IDTCH4 == 0 .and. IND_CH4 > 0 .and. 
     &      ( CH4_YEAR /= YEAR ) ) THEN

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
            CALL GET_GLOBAL_CH4( CH4_YEAR,  .TRUE.,  C3090S,
     &                           C0030S,    C0030N,  C3090N, 
     &                           am_I_Root, Input_Opt        )
         ENDIF

      ENDIF
#endif

      !================================================================
      ! Get concentrations of aerosols in [kg/m3] 
      ! for FAST-JX and optical depth diagnostics
      !=================================================================
      IF ( Input_Opt%LSULF .or. Input_Opt%LCARB .or. 
     &     Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

#if defined( UCX )
         ! Calculate stratospheric aerosol properties (SDE 04/18/13)
         CALL CALC_STRAT_AER( am_I_Root, Input_Opt, 
     &                        State_Met, State_Chm, RC )
         IF ( prtDebug ) 
     &      CALL DEBUG_MSG( '### FLEX_CHEMDR: after CALC_PSC' )
#endif

         ! Skip this section if all these are turned off
         CALL AEROSOL_CONC
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

      ENDIF
      !=================================================================
      ! Call RDAER -- computes aerosol optical depths
      !=================================================================

      ! Call RDAER to compute AOD for FAST-JX (skim, 02/03/11)
      WAVELENGTH = 0
      CALL RDAER( am_I_Root, Input_Opt, State_Met,  RC,
     &            MONTH,     YEAR,      WAVELENGTH      )

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
         CALL RDUST_ONLINE( am_I_Root, Input_Opt,  State_Met,
     &                      SOILDUST,  WAVELENGTH, RC         )
      ELSE
#if !defined( TOMAS )
         ! Don't read dust emissions from disk when using TOMAS,
         ! because TOMAS uses a different set of dust tracers than the 
         ! std code (win, bmy, 1/25/10)
         CALL RDUST_OFFLINE( am_I_Root, Input_Opt, State_Met,
     &                       MONTH,     YEAR,      WAVELENGTH, RC )
#endif
      ENDIF

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### FLEX_CHEMDR: after RDUST' )
      ENDIF

      !=================================================================
      ! Calculate hetchem rates
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Get 1-D index
         JLOOP = JLOP(I,J,L)

         ! If this is a valid grid box
         IF ( JLOOP > 0 ) THEN

            ! Copy aerosol areas into State_Chm
            !  -- Should be moved to AerosolMod once JLOOP is gone
            State_Chm%AeroArea(I,J,L,:) = TAREA(JLOOP,:)
            State_Chm%AeroRadi(I,J,L,:) = ERADIUS(JLOOP,:)

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Loop over tracers and grid boxes
      ! The following 4-loop block mimics partition.F
      ! -- Families are dealt with in the KLUDGE Pt.1 below
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Convert STT from kg to molec/cm3
         STT(I,J,L,N) = STT(I,J,L,N) * Input_Opt%XNUMOL(N) /
     &                 ( State_Met%AIRVOL(I,J,L) * 1e+6_fp *
     &                   Input_Opt%TRACER_COEFF(N,1) )

         ! Copy STT to State_Chm%Species [molec/cm3]
         IF (STTTOCSPEC(N) .NE. 0) THEN 
            State_Chm%Species(I,J,L,STTTOCSPEC(N)) = STT(I,J,L,N)
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      CALL FAMILIES_KLUDGE(am_I_Root,STT,Input_Opt,State_Chm,1,RC)
 
      !=================================================================
      ! Convert STT from molec/cm3 back to kg (required for 
      ! gckpp_HetRates.F90) 
      ! 
      ! NOTE: When we remove the tracer families we can get rid of this 
      ! step by converting units when we copy STT to State_Chm%Species 
      ! and leaving STT in units kg (mps, 1/28/16)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT(I,J,L,N) = STT(I,J,L,N) / Input_Opt%XNUMOL(N) *
     &                  ( State_Met%AIRVOL(I,J,L) * 1e+6_fp *
     &                  Input_Opt%TRACER_COEFF(N,1) )
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================
      CALL FAST_JX( WAVELENGTH, am_I_Root,  Input_Opt, 
     &              State_Met,  State_Chm,  RC         )

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
 100  format('No. of function calls:', i6, /,
     &       'No. of jacobian calls:', i6, /,
     &       'No. of steps:         ', i6, /,
     &       'No. of accepted steps:', i6, /,
     &       'No. of rejected steps ', i6, /,
     $       '       (except at very beginning)',          /,
     &       'No. of LU decompositions:             ', i6, /,
     &       'No. of forward/backward substitutions:', i6, /,
     &       'No. of singular matrix decompositions:', i6, /,
     &       /,
     &       'Texit, the time corresponding to the      ',        /,
     &       '       computed Y upon return:            ', f11.4, /,
     &       'Hexit, last accepted step before exit:    ', f11.4, /,
     &       'Hnew, last predicted step (not yet taken):', f11.4 )

      !-----------------------------------------------------------------
      ! NOTE: The following variables are held THREADPRIVATE and 
      ! therefore do not need to be included in the !$OMP+PRIVATE 
      ! statements below: C, VAR, FIX, RCONST, TIME, TEMP, NUMDEN, 
      ! H2O, PRESS, PHOTOL, and HET. (bmy, 3/28/16)
      !-----------------------------------------------------------------
!$OMP PARALLEL DO
!$OMP+DEFAULT  ( SHARED )
!$OMP+PRIVATE  ( I,    J,     L,      N                   )
!$OMP+PRIVATE  ( YLAT, SCF,   IERR,   RCNTRL,  AVERAGEVAR )  
!$OMP+PRIVATE  ( MASS, START, FINISH, ISTATUS, RSTATE     )
!$OMP+REDUCTION( +:ITIM                                   )
!$OMP+REDUCTION( +:RTIM                                   ) 
!$OMP+REDUCTION( +:TOTSTEPS                               ) 
!$OMP+REDUCTION( +:TOTFUNCS                               )
!$OMP+REDUCTION( +:TOTJACOB                               ) 
!$OMP+REDUCTION( +:TOTACCEP                               )
!$OMP+REDUCTION( +:TOTREJEC                               )
!$OMP+REDUCTION( +:TOTNUMLU                               )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! mje If we are not in the troposphere don't do the chemistry!
         IF (.not. ITS_IN_THE_CHEMGRID(I,J,L,State_Met)) CYCLE

         ! Skipping buffer zone (lzh, 08/10/2014)
         IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
           IF ( J <=         Input_Opt%NESTED_J0W ) CYCLE
           IF ( J >  JJPAR - Input_Opt%NESTED_J0E ) CYCLE
           IF ( I <=         Input_Opt%NESTED_I0W ) CYCLE
           IF ( I >  IIPAR - Input_Opt%NESTED_I0E ) CYCLE
         ENDIF

         ! Grid-box latitude index
         YLAT   = GET_YMID( I,J,L )

         ! mje Set up the physical conditions for the current box
         ! mje These variables are defined by the KPP definition files gckpp.kpp
         ! mje #INLINE F90_GLOBAL    
         TEMP   = State_Met%T(I,J,L)
         PRESS  = GET_PCENTER(I,J,L)

         ! mje Calculate NUMDEN based on ideal gas law (# cm-3)
         NUMDEN = State_Met%AIRNUMDEN(I,J,L)

         ! mje H2O arrives in g/kg needs to be in mol cm-3 
         H2O    = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)

         ! Initialize the KPP solver
         CALL Init_KPP( )

         !===========================================================
         ! Get rates for heterogeneous chemistry
         !===========================================================

         ! Get starting time for rate computation
         CALL CPU_TIME( start )

         ! Zero the HET array (holds hetchem rates)
         HET = 0.0_dp

         ! Copy SCOEFF to SCF, this will avoid an array temporary
         SCF = SCOEFF(I,J,L,:)

         ! Set hetchem rates
         CALL SET_HET( I,         J,         L, 
     &                 State_Chm, State_Met, Input_Opt, SCF )

         ! Copy SCF back into SCOEFF
         SCOEFF(I,J,L,:) = SCF
         
         !===========================================================
         ! Get photolysis rates
         !===========================================================
         IF ( State_Met%SUNCOSmid(I,J) > 0.e+0_fp ) THEN

            ! Daytime: Set up the photolysis rates
            ! The order here is the order in the Fast-J definition files
            ! I've assumed that these are the same as in the text files but this
            ! may have been changed. 
            ! This needs to be checked through more thoroughly
            PHOTOL(1:JVN_) = 0.0_fp

#if defined( UCX ) 
            !%%%%% NEED TO HANDLE THIS BETTER (bmy, 3/28/16)
!            IF ( .not. Input_Opt%LUCX ) THEN ! If not UCX, then parameterize O(1D)
            CALL PHOTRATE_ADJ( am_I_root, State_Met, 
     &                         State_Chm, Input_Opt, 
     &                         I, J, L, H2O, SO4_PHOTFRAC(I,J,L), IERR )
!            ENDIF
#endif

            ! Copy ZPJ from fast_jx_mod.F to PHOTOL array
            PHOTOL(1:JVN_) = ZPJ(L,1:JVN_,I,J)

!            IF (Input_Opt%LUCX) THEN !Adjust SO4 photolysis rate
!               ! There should be a better, more elegant means of indexing this.
!               PHOTOL(100) = PHOTOL(100) * SO4_PHOTFRAC(I,J,L)
!            ENDIF               !LUCX

         ELSE

            ! Nighttime: Photolysis rates are zero
            PHOTOL(1:JVN_) = 0.0_fp

         ENDIF

         !===========================================================
         ! Initialize species concentrations
         !===========================================================
         ! mje Set up the initial conditions of the variable species
         ! mje The Integrator works in # cm-3
         DO N=1,NSPEC
            IF ( CSPECTOKPP(N) .eq. 0) THEN
               C(N)=0d0
            ELSE
               C(N)=State_Chm%Species(I,J,L,CSPECTOKPP(N))
            ENDIF
         ENDDO

!------------------------------------------------------------------------------
! Prior to 2/10/16:
! Set CH4 according to latitude (mps, 2/10/16)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Prior to 2/12/16:
! Default concentrations of MOH and H2 are now handled in flexchem_setup.F90
!  - MOH is set according to region (continental BL, marine BL, free trop)
!------------------------------------------------------------------------------

         C(ind_ACTA)  = 0.0_dp
         C(ind_EOH)   = 0.0_dp
         C(ind_HCOOH) = 0.0_dp
#if defined( UCX )
         IF (Input_Opt%LUCX) THEN ! Set fixed species concentrations 
            C(ind_O2) = 0.2095e+0_dp * NUMDEN
            C(ind_N2) = 0.7808e+0_dp * NUMDEN
            C(ind_H2) = 0.5000e-6_dp * NUMDEN
         ENDIF !LUCX
#endif

         IF ( IDTCH4 == 0 .and. IND_CH4 > 0 ) THEN
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
         ENDIF   ! IDTCH4, IND_CH4

         !### Debug printout, we can remove it later
         IF (I .eq. 10 .and. J .eq. 10 .and. L .eq. 10) THEN
            DO N = 1,NSPEC
               write(*,'(2a,e12.4)') '<> b ', trim(SPC_NAMES(N))
     &              ,C(N)
            ENDDO
         ENDIF

         !===========================================================
         ! Update KPP's rates
         !===========================================================

         ! VAR and FIX are chunks of array C (mps, 2/24/16)
         !
         ! NOTE: Because VAR and FIX are !$OMP+THREADPRIVATE, they
         ! cannot appear in an EQUIVALENCE statement.  Therfore, we
         ! will just copy the relevant elements of C to VAR and FIX
         ! here. (bmy, 3/28/16)
         VAR = C(1:NVAR)
         FIX = C(NVAR+1:NSPEC)

         ! Update the array of rate constants
         CALL Update_RCONST( )

         ! Get time when rate computation finished
         CALL CPU_TIME( finish )

         ! Compute how long it took for KPP to compute rates
         rtim = rtim + finish - start

         !===========================================================
         ! Set options for the KPP Integrator (M. J. Evans)
         !
         ! NOTE: Because RCNTRL(3) is set to an array value that
         ! depends on (I,J,L), we must declare RCNTRL as PRIVATE
         ! within the OpenMP parallel loop and define it just 
         ! before the call to to Integrate. (bmy, 3/24/16)
         !===========================================================

         ! Zero all slots of RCNTRL
         RCNTRL    = 0.0

         ! Starting value for integration time step
         RCNTRL(3) = HSAVE_KPP(I,J,L)

         !===========================================================
         ! Integrate the box forwards
         !===========================================================

         ! Get time before integrator starts
         call CPU_TIME( start )

         ! Call the KPP integrator
         call Integrate( TIN,    TOUT,    ICNTRL, 
     &                   RCNTRL, ISTATUS, RSTATE, IERR )

         ! Get time when integrator ends
         call CPU_TIME( finish )

         ! Compute how long the integrator took 
         itim     = itim + finish - start

         ! Compute other statistics from the integrator
         totfuncs = totfuncs + ISTATUS(1)
         totjacob = totjacob + ISTATUS(2)
         totsteps = totsteps + ISTATUS(3)
         totaccep = totaccep + ISTATUS(4)
         totrejec = totrejec + ISTATUS(5)
         totnumLU = totnumLU + ISTATUS(6)

         ! Zero certain species
         C(ind_ACTA)  = 0.e0_dp
         C(ind_EOH)   = 0.e0_dp
         C(ind_HCOOH) = 0.e0_dp

         ! Try another time if it failed
!         IF ( IERR < 0 ) THEN
!            write(6,*) ''
!            write(6,100) istatus(1:8),rstate(1:3)
!            write(6,*) ''
!            write(ERR_MSG,'(a,i3)') 'Integrator error code :',IERR
!            write(6,*)  'I, J, L ', I, J, L
!
!            ! Reset first time step and start concentrations
!            ! Retry the integration with non-optimized
!            ! settings
!            RCNTRL(3)  = 0e+0_fp
!!            CALL Init_KPP( )
!            CALL Update_RCONST( )
!            call Integrate( TIN,    TOUT,    ICNTRL, 
!     &                      RCNTRL, ISTATUS, RSTATE, IERR )
!            IF ( IERR < 0 ) THEN 
!               print*, 'failed twice !!! '
!!               CALL ERROR_STOP(ERR_MSG, 'INTEGRATE_KPP')
!            ENDIF
!            
!         ENDIF

         ! mje calculate some diagnostics
         ! mje calculate d[X]/dt for all the species
         !CALL FUN(VAR, FIX, RCONST, DERIV)

         ! mje Copy d[X]/dt for this box into the global array
         ! mje x Copy Fluxes through each reation into the global array
         !GLOB_DERIV(I,J,L,:)=DERIV
         !GLOB_RATES(I,J,L,:)=A

         ! Copy VAR and FIX back into C (mps, 2/24/16)
         C(1:NVAR)       = VAR(:)
         C(NVAR+1:NSPEC) = FIX(:)

         ! Save for next integration time step
         HSAVE_KPP(I,J,L) = RSTATE(Nhnew)

         ! Save rate constants in global array
         GLOB_RCONST(I,J,L,:) = RCONST(:)

         !=============================================================
         ! mje some simple diagnostics to check that things make sense
         ! mje these can be removed once we are happy with everything
         !=============================================================

         ! NOTE: We need to zero these here within the OpenMP loop
         ! and declare these as !$OMP+PRIVATE. (bmy, 3/28/16)
         AVERAGEVAR = 0.0_fp
         MASS       = 0.0_fp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% BLOCK THIS OUT FOR NOW (bmy, 3/28/16)
!         ! mje some simple diagnostics to check that things make sense
!         ! mje these can be removed once we are happy with everything
!         DO N=1,NSPEC
!            AVERAGEVAR(N)=AVERAGEVAR(N)+
!     &                    C(N)*State_Met%AD(I,J,L)
!         ENDDO
!         MASS=MASS+State_Met%AD(I,J,L)
!
!         ! mje Check we have no negative values and copy the concentrations
!         ! mje calculated into State_Chm%Species
!         DO N=1,NSPEC
!            C(N) = MAX( C(N), 0E0_dp )
!            IF ( CSPECTOKPP(N) .eq. 0) CYCLE
!            State_Chm%Species(I,J,L,CSPECTOKPP(N))= real(C(N),kind=fp)
!            IF (I.eq.36.and.J.eq.5.and.L.eq.1.and.N.eq.IND_OH) THEN
!               write(*,'(a,2e12.4)') 'SPC_OH : ',c(N),
!     &              State_Chm%Species(I,J,L,CSPECTOKPP(N))
!            ENDIF
!         ENDDO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Debug
         IF (I .eq. 10 .and. J .eq. 10 .and. L .eq. 10) THEN
            DO N = 1,NSPEC
               write(*,'(2a,e12.4)') '<> a ', trim(SPC_NAMES(N)),
     &              C(N)
            ENDDO
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      write(*,'(a,F10.3)') 'Flex Rate Time: ', rtim
      write(*,'(a,F10.3)') 'Flex Intg Time: ', itim
      write(*,'(a,I9)') 'Flex Function Calls: ', totfuncs
      write(*,'(a,I9)') 'Flex Jacobian Calls: ', totjacob
      write(*,'(a,I9)') 'Flex Total Steps   : ', totsteps
      write(*,'(a,I9)') 'Flex Rejected Steps: ', totrejec
      write(*,'(a,I9)') 'Flex LU Decompos.  : ', totnumLU

!     mje write out tropospheric mass weighted averages
!      DO N=1,NSPEC
!         WRITE (6,*) SPC_NAMES(N), AVERAGEVAR(N)/MASS
!      ENDDO	

      !=================================================================
      ! Convert STT from kg to molec/cm3
      !
      ! NOTE: When we remove the tracer families we can get rid of this 
      ! step because STT will be left in units kg and we can convert 
      ! units from molec/cm3/s to kg when we copy State_Chm%Species to 
      ! STT (mps, 1/28/16)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )!
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         STT(I,J,L,N) = STT(I,J,L,N) * Input_Opt%XNUMOL(N) /
     &                  ( State_Met%AIRVOL(I,J,L) * 1e+6_fp *
     &                  Input_Opt%TRACER_COEFF(N,1) )
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! Kludge for families: species to tracers
      CALL FAMILIES_KLUDGE(am_I_Root,STT,Input_Opt,State_Chm,2,RC)

      !=================================================================
      ! Copy updated tracer concentrations back into STT
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Copy new tracer concentrations from State_Chm%Species to STT
         ! [molec/cm3]
         IF (STTTOCSPEC(N) .NE. 0) THEN
            STT(I,J,L,N) = State_Chm%Species(I,J,L,STTTOCSPEC(N))
         ENDIF
         
         ! Convert STT from molec/cm3 back to kg
         STT(I,J,L,N) = STT(I,J,L,N) / Input_Opt%XNUMOL(N) *
     &                  ( State_Met%AIRVOL(I,J,L) * 1e+6_fp *
     &                  Input_Opt%TRACER_COEFF(N,1) )
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

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

      !=================================================================
      ! Call OHSAVE which saves info on OH AND HO2 concentrations
      !=================================================================
      CALL OHSAVE( State_Met, State_Chm )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### FLEX_CHEMDR: after OHSAVE' )
      ENDIF

      !=================================================================
      ! Save quantities for computing mean OH lifetime
      !=================================================================
      CALL DO_DIAG_OH( State_Met, State_Chm )

      !### Debug
      IF ( prtDebug ) THEN
         CALL DEBUG_MSG( '### FLEX_CHEMDR: after DO_DIAG_OH' )
      ENDIF


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR DIAGNOSTIC PURPOSES ONLY:
#if defined( DEVEL )
      ! mje  STT currently in kg / grid box convert to v/v
      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT)

      DO N=1,N_TRACERS

            ! Construct diagnostic name
            DiagnName = TRIM( Input_Opt%TRACER_NAME(N) )
            
            ! Update the deposition velocity diagostics 
            CALL Diagn_Update( am_I_Root,
     &                         cName   = TRIM( DiagnName ),
     &                         Array3D = STT(:,:,:,N),
     &                         COL     = Input_Opt%DIAG_COLLECTION,
     &                         RC      = HCRC )

            ! Stop with error if the diagnostics were unsuccessful
            IF ( HCRC /= HCO_SUCCESS ) THEN
               CALL ERROR_STOP( 'Cannot update tracer 
     &                           concentration diagnostics',
     &                          'Flex_Chemdr')
            ENDIF
      ENDDO

      ! mje convert the STT array from v/v to kg / grid box
      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT)

      DO N=1,NSPEC

            ! Construct diagnostic name
            DiagnName = 'SPC_'//TRIM( SPC_NAMES(N) )
            
            ! Update the deposition velocity diagostics 
            CALL Diagn_Update( am_I_Root,
     &                 cName   = TRIM( DiagnName ),
     &                 Array3D = State_Chm%Species(:,:,:,CSPECTOKPP(N)),
     &                 COL     = Input_Opt%DIAG_COLLECTION,
     &                 RC      = HCRC )

            ! Stop wieth error if the diagnostics were unsuccessful
            IF ( HCRC /= HCO_SUCCESS ) THEN
               CALL ERROR_STOP( 'Cannot update tracer 
     &                           concentration diagnostics',
     &                          'Flex_Chemdr')
            ENDIF
      ENDDO

      ! Write out Rate values
      DO N=1,NREACT
         ! Construct diagnostic name
         write(DiagnName,'(a)')
     &        'RR_' // TRIM( ADJUSTL(EQN_NAMES(N)) )
            
         ! Update the diagostics 
         CALL Diagn_Update( am_I_Root,
     &        cName   = TRIM( DiagnName ),
     &        Array2D = GLOB_RCONST(:,:,1,N),
     &        COL     = Input_Opt%DIAG_COLLECTION,
     &        RC      = HCRC )

         ! Stop wieth error if the diagnostics were unsuccessful
         IF ( HCRC /= HCO_SUCCESS ) THEN
            CALL ERROR_STOP( 'Cannot update tracer 
     &           concentration diagnostics',
     &          'Flex_Chemdr')
         ENDIF
      ENDDO

      ! Update the diagostics 
      CALL Diagn_Update( am_I_Root,
     &     cName   = TRIM( 'T' ),
     &     Array3D = State_Met%T(:,:,:),
     &     COL     = Input_Opt%DIAG_COLLECTION,
     &     RC      = HCRC )

      ! Stop with error if the diagnostics were unsuccessful
      IF ( HCRC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'Cannot update temperature diagnostic',
     &        'Flex_Chemdr')
      ENDIF

      CALL Diagn_Update( am_I_Root,
     &     cName   = TRIM( 'STK_N2O5' ),
     &     Array3D = SCOEFF(:,:,:,1),
     &     COL     = Input_Opt%DIAG_COLLECTION,
     &     RC      = HCRC )

      CALL Diagn_Update( am_I_Root,
     &     cName   = TRIM( 'STK_HBr' ),
     &     Array3D = SCOEFF(:,:,:,2),
     &     COL     = Input_Opt%DIAG_COLLECTION,
     &     RC      = HCRC )

      CALL Diagn_Update( am_I_Root,
     &     cName   = TRIM( 'STK_HBrICE' ),
     &     Array3D = SCOEFF(:,:,:,3),
     &     COL     = Input_Opt%DIAG_COLLECTION,
     &     RC      = HCRC )
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! mje set rid of the pointers created
      NULLIFY (STT)

      ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
      FIRSTCHEM = .FALSE.

      END SUBROUTINE FLEX_CHEMDR
!EOC


