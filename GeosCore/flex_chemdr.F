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
      USE FLEXCHEM_SETUP_MOD,   ONLY : STTTOCSPEC, CSPECTOKPP, HSAVE_KPP
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
      USE UCX_MOD,              ONLY : CALC_STRAT_AER

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
      LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC ! Success or failure
! 
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  18 Dec 2015 - M. Sulprizio- Add calls to OHSAVE and DO_DIAG_OH
!  22 Dec 2015 - M. Sulprizio- Set NUMDEN to State_Met%AIRNUMDEN for conversions
!                              between v/v and molec/cm3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, N, S1, S2, JLOOP, MONTH, YEAR
      INTEGER             :: WAVELENGTH
      INTEGER             :: HCRC
      INTEGER             :: N_TRACERS
      LOGICAL             :: prtDebug
      LOGICAL             :: LDUST
      CHARACTER(LEN=155)  :: DiagnName

      ! SAVEd variables
      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.
      INTEGER, SAVE       :: CH4_YEAR  = -1

      ! For grid-box latitude
      REAL(fp)            :: YLAT

      ! For interannually-varying Methane
      REAL(fp)            :: C3090S, C0030S, C0030N, C3090N

      ! mje Variables for holding rate diagnoastics
      ! mje DERIV: Concentrations derivatives for one box d[X]/dt
      ! mje GLOB_DERIV: Store global values of concentrations devivatives
      ! mje GLOB_RATES: Stores global values of fluxes through each reaction
      REAL(kind=dp)       :: DERIV(NSPEC)
      REAL(kind=dp)       :: GLOB_DERIV(IIPAR,JJPAR,LLPAR,NSPEC)
      REAL(kind=dp)       :: GLOB_RATES(IIPAR,JJPAR,LLPAR,NREACT)
      REAL(kind=dp)       :: GLOB_RCONST(IIPAR,JJPAR,LLPAR,NREACT)

      ! mje KPP variables
      INTEGER             :: IERR
      REAL(kind=fp)       :: T, TIN, TOUT
      REAL(kind=dp)       :: RCNTRL(20)
      REAL(kind=dp)       :: RSTATE(20)
      INTEGER             :: ICNTRL(20)
      INTEGER             :: ISTATUS(20)

      ! mje Array to hold the molecular mass of the species
      REAL(kind=fp)       :: TCVV(Input_Opt%N_TRACERS)

      ! mje STT holds the concentrations of Tracers
      REAL(kind=fp), POINTER :: STT(:,:,:,:)

      ! mje Array for bug checking
      REAL(kind=fp)       :: AverageVAR(NSPEC)
      REAL(kind=fp)       :: Mass, Start, Finish, rtim, itim

      CHARACTER(LEN=255)  :: ERR_MSG

      ! MSL - Temporary variables for the family kludge (see below)
      REAL(kind=fp)       ::   SUM(IIPAR,JJPAR,LLPAR)
      REAL(kind=fp)       :: QTEMP(IIPAR,JJPAR,LLPAR)

      REAL(kind=fp)       :: SCOEFF(IIPAR,JJPAR,LLPAR,3)

      prtDebug = Input_Opt%LPRT
      LDUST    = Input_Opt%LDUST

      ! Get month and year
      MONTH     = GET_MONTH()
      YEAR      = GET_YEAR()

      ! mje Copy Number of tracers (advected species) into N_TRACERS
      N_TRACERS = Input_Opt%N_TRACERS

      ! mje Copy the molecular mass of the Tracers into TCVV
      TCVV      = Input_Opt%TCVV

      ! mje link the STT Pointer to the State_Chm%Tracers stucture
      STT =>  State_Chm%Tracers

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

      ! MSL  For calculating hetchem rates
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

      END DO
      END DO
      END DO

      ! Loop over tracers and grid boxes
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Convert STT from kg to molec/cm3
         STT(I,J,L,N) = STT(I,J,L,N) * Input_Opt%XNUMOL(N) /
     &                  ( State_Met%AIRVOL(I,J,L) * 1e+6_fp *
     &                  Input_Opt%TRACER_COEFF(N,1) )

         ! Copy STT to State_Chm%Species [molec/cm3]
         IF (STTTOCSPEC(N) .NE. 0) THEN 
            State_Chm%Species(I,J,L,STTTOCSPEC(N)) = STT(I,J,L,N)
         ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      ! K -- L -- U -- D -- G -- E -- C -- O -- D -- E -- ! -- ! -- ! -- !
      ! THIS IS A TEMPORARY FIX AND NEEDS TO BE RESOLVED THROUGHOUT
      ! THE MODEL AS SOON AS IT APPEARS POSSIBLE TO DO SO!
      ! -- The following code ensures that the two remaining tracer
      !    families are dealt with appropriately.
      !    The Tracer MMN is a sum of the species MVKN and MACRN
      !    The Tracer ISOPN is a sum of ISOPND and ISOPNB
      !       The tracer restart file for the pre-flex GEOS-Chem, includes
      !    these families, but the removal of the routines
      !    "lump" and "partition" killed the code resposible for them.
      !       Here, we want to install a hard-code fix with the complete
      !    expectation that we will no longer use species families. Thus,
      !    when possible, the two remaining families need to be removed, 
      !    and this kludge disabled.
      !   M.S.L. - Jan., 5 2016
      !
      ! Part 1: From Tracers to Species
      ! -- Process MMN family
      ! -- Get the tracer and species indices. Currently hardwired.
      N  = 60   ! MMN
      S1 = 122  ! MVKN
      S2 = 121  ! MACRN
      STT(:,:,:,N) = STT(:,:,:,N)*Input_Opt%TRACER_COEFF(N,1) ! A correction...
      ! ... to account for the division by TRACER_COEFF above.
      SUM   = (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &       +(State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
      ! -- -- First, do MVKN
      QTEMP = (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &        / SUM
      State_Chm%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )
      ! -- -- Then, do MACRN
      QTEMP = (State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
     &        / SUM
      State_Chm%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )
      ! -- Process ISOPN family
      ! -- Get the tracer and species indices
      N  = 55  ! ISOPN Tracer
      S1 = 109 ! ISOPND Species Index
      S2 = 110 ! ISOPNB Species Index
      STT(:,:,:,N) = STT(:,:,:,N)*Input_Opt%TRACER_COEFF(N,1) ! A correction...
      ! ... to account for the division by TRACER_COEFF above.
      SUM   = (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &       +(State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
      ! -- -- First, do ISOPND
      QTEMP = (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &        / SUM
      State_Chm%Species(:,:,:,S1) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )
      ! -- -- Then, do ISOPNB
      QTEMP = (State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
     &        / SUM
      State_Chm%Species(:,:,:,S2) = MAX( QTEMP*STT(:,:,:,N), 1.E-99_fp )
      ! E -- N -- D -- O -- F -- K -- L -- U -- D -- G -- E -- -- P -- T -- 1
 
      ! Convert STT from molec/cm3 back to kg (required for gckpp_HetRates.F90)
      ! NOTE: When we remove the tracer families we can get rid of this step
      ! by converting units when we copy STT to State_Chm%Species and leaving
      ! STT in units kg (mps, 1/28/16)
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

      !     imje  Calculate photolysis rates
      !=================================================================
      ! Call photolysis routine to compute J-Values
      !=================================================================
      CALL FAST_JX( WAVELENGTH, am_I_Root, Input_Opt, State_Met, 
     &              State_Chm,  RC )

      ! mje Set up conditions for the integration
      ! mje chemical timestep and convert it to seconds.
      DT   = GET_TS_CHEM() * 60d0
      T    = 0d0
      TIN  = T
      TOUT = T + DT

      ! mje Set up integration convergence conditions
      !RTOL = 1.d-1
      !ATOL = 1.d0
      RTOL = 2e-1_dp
      ATOL = 1e-2_dp
      
      ! mje Set up internal integrator parameters 
      ICNTRL = 0
      RCNTRL = 0.
      
      ICNTRL(1) = 1             ! Autonomous
      ICNTRL(2) = 0             ! Nonautonomous
      
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
      ICNTRL(7) = 1             ! No adjoint

      ! mje setup some internal diagnostics
      AVERAGEVAR(:) = 0.d0
      MASS          = 0.d0

      !=================================================================
      ! Solve Chemistry
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

      write(*,'(a,e12.4,i)') '<> CO2 : ', C(ind_CO2), ind_CO2

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, start, finish )
!!$OMP+PRIVATE( PHOTOL )
!$OMP+PRIVATE( DERIV )
!$OMP+PRIVATE( ISTATUS, RSTATE, IERR )
!!$OMP+PRIVATE( ITEMPK, RO1D, RO1DplH2O, RO1DplH2, H2O )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! mje If we are not in the troposphere don't do the chemistry!
         IF (.not. ITS_IN_THE_CHEMGRID(I,J,L,State_Met)) CYCLE

         call cpu_time(start)
         HET(:,:) = 0.
         CALL SET_HET(I,J,L,State_Chm,State_Met,Input_Opt,
     &        scoeff(I,J,L,1:3))

         ! Grid-box latitude index
         YLAT  = GET_YMID( I,J,L )

         ! mje Set up the physical conditions for the current box
         ! mje These variables are defined by the KPP definition files gckpp.kpp
         ! mje #INLINE F90_GLOBAL    
         TEMP   = State_Met%T(I,J,L)
         PRESS  = GET_PCENTER(I,J,L)

         ! mje Calculate NUMDEN based on ideal gas law (# cm-3)
         NUMDEN = State_Met%AIRNUMDEN(I,J,L)

         ! mje H2O arrives in g/kg needs to be in mol cm-3 
         H2O = State_Met%AVGW(I,J,L)*State_Met%AIRNUMDEN(I,J,L)

!         H2O    = State_Met%SPHU(I,J,L)*
!     &            (0.21*32d0+0.79*28d0)/(18d0)/1000.
!         H2O    = H2O * NUMDEN

         ! Set up te photolysis rates the order here is the order in the Fast-J
         ! definition files
         ! I've assumed that these are the same as in the text files but this
         ! may have been changed. 
         ! This needs to be checked through more thoroughly
         CALL PHOTRATE_ADJ( am_I_root, State_Met, 
     &                      State_Chm, Input_Opt, 
     &                      I, J, L, H2O, NUMDEN, IERR )
         PHOTOL(1:JVN_) =    ZPJ(L,1:JVN_,I,J)

         ! mje Set up the initial conditions of the variable species
         ! mje The Integrator works in # cm-3
         DO N=1,NSPEC
            IF ( CSPECTOKPP(N) .eq. 0) CYCLE
            C(N)=State_Chm%Species(I,J,L,CSPECTOKPP(N))
         ENDDO

         ! Debug
         IF (I .eq. 10 .and. J .eq. 10 .and. L .eq. 10) THEN
            DO N = 1,NSPEC
               write(*,'(2a,e12.4,i)') '<> b ', SPC_NAMES(N), C(N), N
            ENDDO
         ENDIF

!------------------------------------------------------------------------------
! Prior to 2/10/16:
! Set CH4 according to latitude (mps, 2/10/16)
!         C(ind_CH4) = 1770.e-9_dp * NUMDEN
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Prior to 2/12/16:
! Default concentrations of MOH and H2 are now handled in flexchem_setup.F90
!  - MOH is set according to region (continental BL, marine BL, free trop)
!  - H2  is set to default values from globchem.dat (5.0e-7)
! (mps, 2/12/16)
!         C(ind_MOH) = 0.5e-9_dp   * NUMDEN
!         C(ind_H2)  = 0.5e-6_dp   * NUMDEN
!------------------------------------------------------------------------------

         ! Set CH4 according to latitude
         ! Convert from [ppbv CH4] to [molec CH4/cm3]
         IF ( YLAT < -30e+0_fp ) THEN
            C(ind_CH4) = C3090S * 1e-9_dp * NUMDEN
         ELSE IF ( YLAT >= -30e+0_fp .and. YLAT < 0e+0_fp  ) THEN
            C(ind_CH4) = C0030S * 1e-9_dp * NUMDEN
         ELSE IF ( YLAT >=   0e+0_fp .and. YLAT < 30e+0_fp ) THEN
            C(ind_CH4) = C0030N * 1e-9_dp * NUMDEN
         ELSE
            C(ind_CH4) = C3090N * 1e-9_dp * NUMDEN
         ENDIF

         ! starting value for integration time step
         RCNTRL(3) = HSAVE_KPP(I,J,L)

         ! mje  Update KPP's rates
         CALL Update_RCONST( )
         call cpu_time(finish)
         rtim = rtim + finish - start
         GLOB_RCONST(I,J,L,:) = RCONST(:)

         ! mje Integrate the box forwards
         call cpu_time(start)
         call Integrate( TIN,    TOUT,    ICNTRL, 
     &                   RCNTRL, ISTATUS, RSTATE, IERR )
         call cpu_time(finish)
         itim = itim + finish - start


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

         ! mje some simple diagnostics to check that things make sense
         ! mje these can be removed once we are happy with everything
         DO N=1,NSPEC
            AVERAGEVAR(N)=AVERAGEVAR(N)+
     &                    C(N)*State_Met%AD(I,J,L)
         ENDDO
         MASS=MASS+State_Met%AD(I,J,L)

         ! Debug
         IF (I .eq. 10 .and. J .eq. 10 .and. L .eq. 10) THEN
            DO N = 1,NSPEC
               write(*,'(2a,e12.4)') '<> a ', SPC_NAMES(N), C(N)
            ENDDO
         ENDIF

         ! mje Check we have no negative values and copy the concentrations
         ! mje calculated into State_Chm%Species
         DO N=1,NSPEC
            C(N) = MAX( C(N), 0E0_dp )
            IF ( CSPECTOKPP(N) .eq. 0) CYCLE
            State_Chm%Species(I,J,L,CSPECTOKPP(N))= real(C(N),kind=fp)
         ENDDO

         ! save next integration time step
         HSAVE_KPP(I,J,L) = RSTATE(Nhnew)

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      write(*,'(a,F7.3)') 'Flex Rate Time: ', rtim
      write(*,*) 'Flex Intg Time: ', itim

!     mje write out tropospheric mass weighted averages
!      DO N=1,NSPEC
!         WRITE (6,*) SPC_NAMES(N), AVERAGEVAR(N)/MASS
!      ENDDO	

      ! Convert STT from kg to molec/cm3
      ! NOTE: When we remove the tracer families we can get rid of this step
      ! because STT will be left in units kg and we can convert units from
      ! molec/cm3/s to kg when we copy State_Chm%Species to STT (mps, 1/28/16)
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

      ! K -- L -- U -- D -- G -- E -- -- P -- T -- 2
      ! Part 2: From Species to Tracers
      ! -- Process MMN family
      ! -- Get the tracer and species indices
      N  = 60  ! MMN
      S1 = 122 ! MVKN
      S2 = 121 ! MACRN
      STT(:,:,:,N) = 
     &        (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &       +(State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
      ! -- Process ISOPN family
      ! -- Get the tracer and species indices
      N  = 55  ! ISOPN
      S1 = 109 ! ISOPND
      S2 = 110 ! ISOPNB
      STT(:,:,:,N) = 
     &        (State_Chm%Species(:,:,:,S1)*Input_Opt%TRACER_COEFF(N,1))
     &       +(State_Chm%Species(:,:,:,S2)*Input_Opt%TRACER_COEFF(N,2))
      ! E -- N -- D -- O -- F -- K -- L -- U -- D -- G -- E -- -- P -- T -- 2

      ! Loop over tracers
      DO N=1,N_TRACERS
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! Write out Rate values
      DO N=1,NREACT
         ! Construct diagnostic name
         write(DiagnName,'(a)')
     &        'RR_' // TRIM( ADJUSTL(EQN_NAMES(N)) )
            
         ! Update the diagostics 
         CALL Diagn_Update( am_I_Root,
     &        cName   = TRIM( DiagnName ),
     &        Array3D = GLOB_RCONST(:,:,:,N),
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

      ! mje set rid of the pointers created
      NULLIFY (STT)

      ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
      FIRSTCHEM = .FALSE.

      END SUBROUTINE FLEX_CHEMDR
!EOC


