!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
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
!  29 Dec 2017 - C. Keller   - Make HSAVE_KPP public (needed for GEOS-5 restart)
!  24 Jan 2018 - E. Lundgren - Pass error handling up if RC is GC_FAILURE
!  11 Oct 2018 - M. Sulprizio- Move HSAVE_KPP to State_Chm, rename as KPPHvalue
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags (and logicals to denote if species are present)
  INTEGER               :: id_OH, id_HO2, id_O3P, id_O1D, id_CH4
#if defined( MODEL_GEOS )
  INTEGER               :: id_O3
  INTEGER               :: id_A3O2, id_ATO2, id_B3O2, id_BRO2, id_DHPCARP
  INTEGER               :: id_DIBOO,id_ETO2, id_HC5OO, id_IEPOXOO
  INTEGER               :: id_INPN, id_ISNOOA, id_ISNOOB, id_ISNOHOO, id_LIMO2
  INTEGER               :: id_MAOPO2, id_MO2, id_MRO2, id_PIO2, id_PO2
  INTEGER               :: id_PRNI, id_R4NI, id_R4O2, id_RIO2, id_TRO2
  INTEGER               :: id_VRO2, id_XRO2
#endif
  LOGICAL               :: ok_OH, ok_HO2, ok_O1D, ok_O3P

  ! Diagnostic flags
  LOGICAL               :: Do_Diag_OH_HO2_O1D_O3P
  LOGICAL               :: Do_ND43
#if defined( MODEL_GEOS )
  LOGICAL               :: Archive_O3concAfterchem
  LOGICAL               :: Archive_RO2concAfterchem
#endif

  ! SAVEd scalars
  INTEGER,  SAVE        :: PrevDay   = -1
  INTEGER,  SAVE        :: PrevMonth = -1
  
  ! Arrays
  INTEGER,  ALLOCATABLE         :: ND65_KPP_Id(:      )  ! Indices for ND65 bpch diag
  REAL(f4), ALLOCATABLE         :: JvCountDay (:,:,:  )  ! For daily   avg of J-values
  REAL(f4), ALLOCATABLE         :: JvCountMon (:,:,:  )  ! For daily   avg of J-values
  REAL(f4), ALLOCATABLE         :: JvSumDay   (:,:,:,:)  ! For monthly avg of J-values
  REAL(f4), ALLOCATABLE         :: JvSumMon   (:,:,:,:)  ! For monthly avg of J-values

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
    USE CMN_FJX_MOD
    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
#if defined( BPCH_DIAG )
    USE CMN_DIAG_MOD,         ONLY : ND52
    USE DIAG_MOD,             ONLY : AD65,  AD52, ad22
    USE DIAG20_MOD,           ONLY : DIAG20, POx, LOx
#endif
    USE DIAG_OH_MOD,          ONLY : DO_DIAG_OH
    USE DUST_MOD,             ONLY : RDUST_ONLINE, RDUST_OFFLINE
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FAST_JX_MOD,          ONLY : PHOTRATE_ADJ, FAST_JX
    USE GCKPP_HetRates,       ONLY : SET_HET
    USE GCKPP_Monitor,        ONLY : SPC_NAMES, FAM_NAMES
    USE GCKPP_Parameters
    USE GCKPP_Integrator,     ONLY : INTEGRATE, NHnew
    USE GCKPP_Function 
    USE GCKPP_Model
    USE GCKPP_Global
    USE GCKPP_Rates,          ONLY : UPDATE_RCONST, RCONST
    USE GCKPP_Initialize,     ONLY : Init_KPP => Initialize
#if defined( MODEL_GEOS )
    USE GcKPP_Util,           ONLY : Get_OHreactivity
#endif
    USE GC_GRID_MOD,          ONLY : GET_YMID
    USE GEOS_Timers_Mod
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
    USE TIME_MOD,             ONLY : Get_Day
    USE TIME_MOD,             ONLY : Get_Month
    USE TIME_MOD,             ONLY : Get_Year
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
    USE UCX_MOD,              ONLY : CALC_STRAT_AER
    USE UCX_MOD,              ONLY : SO4_PHOTFRAC
    USE UCX_MOD,              ONLY : UCX_NOX
    USE UCX_MOD,              ONLY : UCX_H2SO4PHOT
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
!  21 Dec 2017 - R. Yantosca - Add netCDF diagnostics for J-values, prod/loss
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  18 Jan 2018 - R. Yantosca - Now do photolysis for all levels, so that 
!                              J-values can be saved up to the atm top
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: prtDebug,  IsLocNoon
    INTEGER                :: I,         J,        L,         N
    INTEGER                :: NA,        F,        SpcID,     KppID
    INTEGER                :: P,         MONTH,    YEAR,      WAVELENGTH
    INTEGER                :: TotSteps,  TotFuncs, TotJacob,  TotAccep
    INTEGER                :: TotRejec,  TotNumLU, HCRC,      IERR
    INTEGER                :: Day
    REAL(fp)               :: Start,     Finish,   rtim,      itim
    REAL(fp)               :: SO4_FRAC,  YLAT,     T,         TIN
    REAL(fp)               :: JNoon_Fac, TOUT

    ! Strings
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc

    ! SAVEd scalars
    LOGICAL,  SAVE         :: FIRSTCHEM = .TRUE.
    INTEGER,  SAVE         :: CH4_YEAR  = -1
    REAL(fp), SAVE         :: C3090S,   C0030S,   C0030N,    C3090N

    ! Arrays
    INTEGER                :: ICNTRL     (                  20               )
    INTEGER                :: ISTATUS    (                  20               )
    REAL(dp)               :: RCNTRL     (                  20               )
    REAL(dp)               :: RSTATE     (                  20               )
#if defined( MODEL_GEOS )
    REAL(f4)               :: GLOB_RCONST(IIPAR,JJPAR,LLPAR,NREACT           )
    REAL(f4)               :: GLOB_JVAL  (IIPAR,JJPAR,LLPAR,JVN_             )
#else
    REAL(dp)               :: GLOB_RCONST(IIPAR,JJPAR,LLPAR,NREACT           )
#endif
    REAL(fp)               :: Before     (IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)

    ! For tagged CO saving
    REAL(fp)               :: LCH4, PCO_TOT, PCO_CH4, PCO_NMVOC

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! For testing only, may be removed later (mps, 4/26/16)
    LOGICAL                :: DO_HETCHEM

#if defined( MODEL_GEOS )
    ! OH reactivity
    LOGICAL                :: DoOHreact
    REAL(fp)               :: OHreact
    REAL(dp)               :: Vloc(NVAR), Aout(NREACT)
#endif

    !=======================================================================
    ! Do_FlexChem begins here!
    !=======================================================================

    ! Initialize
    RC        =  GC_SUCCESS
    ErrMsg    =  ''
    ThisLoc   =  ' -> at Do_FlexChem (in module GeosCore/flexchem_mod.F)'
    SpcInfo   => NULL()
    prtDebug  =  ( Input_Opt%LPRT .and. am_I_Root )
    itim      =  0.0_fp
    rtim      =  0.0_fp
    totsteps  =  0
    totfuncs  =  0
    totjacob  =  0
    totaccep  =  0
    totrejec  =  0
    totnumLU  =  0
    Day       =  Get_Day()    ! Current day
    Month     =  Get_Month()  ! Current month
    Year      =  Get_Year()   ! Current year

    ! Turn heterogeneous chemistry and photolysis on/off here
    ! This is for testing only and may be removed later (mps, 4/26/16)
    DO_HETCHEM  = .TRUE.

    ! Remove debug output
    !IF ( FIRSTCHEM .AND. am_I_Root ) THEN
    !   WRITE( 6, '(a)' ) REPEAT( '#', 32 )
    !   WRITE( 6, '(a,l,a)' ) '# Do_FlexChem: DO_HETCHEM  =', &
    !                                         DO_HETCHEM,  ' #'
    !   WRITE( 6, '(a)' ) REPEAT( '#', 32 )
    !ENDIF

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_Loss  ) State_Diag%Loss  = 0.0_f4
    IF ( State_Diag%Archive_Prod  ) State_Diag%Prod  = 0.0_f4
    IF ( State_Diag%Archive_JVal  ) State_Diag%JVal  = 0.0_f4
    IF ( State_Diag%Archive_JNoon ) State_Diag%JNoon = 0.0_f4

#if defined( MODEL_GEOS )
    GLOB_RCONST = 0.0_f4
    GLOB_JVAL   = 0.0_f4
   
    ! testing only
    IF ( Input_Opt%NN_RxnRates > 0 ) State_Diag%RxnRates(:,:,:,:) = 0.0 

    ! OH reactivity
    DoOHreact = .FALSE.
    IF ( State_Diag%Archive_OHreactivity ) THEN
       DoOHreact = .TRUE.
       State_Diag%OHreactivity(:,:,:) = 0.0
    ENDIF
#endif 

    !=======================================================================
    ! Get concentrations of aerosols in [kg/m3] 
    ! for FAST-JX and optical depth diagnostics
    !=======================================================================
    IF ( Input_Opt%LSULF .or. Input_Opt%LCARB .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

       ! Special handling for UCX
       IF ( Input_Opt%LUCX ) THEN

          ! Calculate stratospheric aerosol properties (SDE 04/18/13)
          CALL CALC_STRAT_AER( am_I_Root, Input_Opt,                         &
                               State_Met, State_Chm, RC                     )
          
          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Debutg utput
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### Do_FlexChem: after CALC_PSC' )
          ENDIF

       ENDIF

       ! Compute aerosol concentrations
       CALL AEROSOL_CONC( am_I_Root, Input_Opt,  State_Met,                  &
                          State_Chm, State_Diag, RC                         )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "AEROSOL_CONC"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !=======================================================================
    ! Zero out certain species:
    !    - isoprene oxidation counter species (dkh, bmy, 6/1/06)
    !    - isoprene-NO3 oxidation counter species (hotp, 6/1/10)
    !    - if SOA or SOA_SVPOA, aromatic oxidation counter species 
    !      (dkh, 10/06/06)
    !    - if SOA_SVPOA, LNRO2H and LNRO2N for NAP (hotp 6/25/09
    !=======================================================================
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
      
    !=======================================================================
    ! Call RDAER -- computes aerosol optical depths
    !=======================================================================
#if defined( USE_TIMERS )
    CALL GEOS_Timer_End  ( "=> Gas-phase chem",   RC )
    CALL GEOS_Timer_Start( "=> All aerosol chem", RC )
#endif


    ! Call RDAER to compute AOD for FAST-JX (skim, 02/03/11)
    WAVELENGTH = 0
    CALL RDAER( am_I_Root, Input_Opt,  State_Met,  &
                State_Chm, State_Diag, RC,         &
                MONTH,     YEAR,       WAVELENGTH )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "RDAER"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN 
       CALL DEBUG_MSG( '### Do_FlexChem: after RDAER' )
    ENDIF

    !=======================================================================
    ! If LDUST is turned on, then we have online dust aerosol in
    ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
    ! compute aerosol optical depth for FAST-JX, etc.
    !
    ! If LDUST is turned off, then we do not have online dust aerosol
    ! in GEOS-CHEM...so read monthly-mean dust files from disk.
    ! (rjp, tdf, bmy, 4/1/04)
    !=======================================================================
    IF ( Input_Opt%LDUST ) THEN
       CALL RDUST_ONLINE( am_I_Root,  Input_Opt, State_Met,  State_Chm,      &
                          State_Diag, SOILDUST,  WAVELENGTH, RC             )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "RDUST_ONLINE"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE
#if !defined( TOMAS )
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%%% NOTE: RDUST_OFFLINE STILL HAS BPCH CODE AND THEREFORE   %%%% 
       !%%%% IS PROBABLY NOW OBSOLETE.  THIS WILL BE REMOVED WHEN WE %%%%
       !%%%% GET HIGH_RESOLUTION DUST EMISSIONS (bmy, 1/18/18)       %%%%
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Don't read dust emissions from disk when using TOMAS,
       ! because TOMAS uses a different set of dust species than the 
       ! std code (win, bmy, 1/25/10)
       CALL RDUST_OFFLINE( am_I_Root,  Input_Opt, State_Met, State_Chm,      &
                           State_Diag, MONTH,     YEAR,      WAVELENGTH, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "RDUST_OFFLINE"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
#endif
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after RDUST' )
    ENDIF

#if defined( USE_TIMERS )
    CALL GEOS_Timer_End  ( "=> All aerosol chem", RC )
    CALL GEOS_Timer_Start( "=> Gas-phase chem",   RC )
#endif

    !=======================================================================
    ! Archive initial species mass for stratospheric tendency
    !=======================================================================
    IF ( Input_Opt%LUCX  ) THEN

       ! Loop over advected species
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N      )
       DO NA = 1, State_Chm%nAdvect

          ! Get the species ID from the advected species ID
          N = State_Chm%Map_Advect(NA)

          ! Archive initial mass [kg]
          Before(:,:,:,N) = State_Chm%Species(:,:,:,N)

       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, & 
                            State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90')
       RETURN
    ENDIF 
      
    !=======================================================================
    ! Call photolysis routine to compute J-Values
    !=======================================================================
#if defined( USE_TIMERS )
    CALL GEOS_Timer_End  ( "=> Gas-phase chem",     RC )
    CALL GEOS_Timer_Start( "=> FAST-JX photolysis", RC )
#endif

    ! Do Photolysis
    CALL FAST_JX( WAVELENGTH, am_I_Root,  Input_Opt, &
                  State_Met,  State_Chm,  State_Diag, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "FAST_JX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#if defined( USE_TIMERS )
    CALL GEOS_Timer_End  ( "=> FAST-JX photolysis", RC )
    CALL GEOS_Timer_Start( "=> Gas-phase chem",     RC )
#endif

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after FAST_JX' )
    ENDIF

#if defined( MODEL_GEOS )
    ! Init diagnostics
    IF ( ASSOCIATED(State_Diag%KppError) ) THEN
       State_Diag%KppError(:,:,:) = 0.0
    ENDIF
#endif

    !=======================================================================
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
    !=======================================================================

    !%%%%% TIMESTEPS %%%%%

    ! mje Set up conditions for the integration
    ! mje chemical timestep and convert it to seconds.
    DT        = GET_TS_CHEM() ! [s]
    T         = 0d0
    TIN       = T
    TOUT      = T + DT

    ! Factor that we need to multiply the JNoon netCDF diagnostic by 
    ! to  account for the number of times the diagnostic array will be 
    ! divided each day.
    JNoon_Fac = 86400.0_fp / DT

    !%%%%% CONVERGENCE CRITERIA %%%%%

    ! Absolute tolerance
    ATOL      = 1e-2_dp    

    ! Relative tolerance
    IF ( Input_Opt%LUCX  ) THEN
       ! UCX-based mechanisms
       !RTOL      = 2e-2_dp
       !RTOL      = 1e-2_dp
       RTOL      = 0.5e-2_dp
    ELSE
       ! Non-UCX mechanisms
       RTOL      = 1e-2_dp
    ENDIF

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

    !=======================================================================
    ! %%%%% SOLVE CHEMISTRY -- This is the main KPP solver loop %%%%%
    !=======================================================================
100 format('No. of function calls:', i6, /,                                 &
           'No. of jacobian calls:', i6, /,                                 &
           'No. of steps:         ', i6, /,                                 &
           'No. of accepted steps:', i6, /,                                 &
           'No. of rejected steps ', i6, /,                                 &
           '       (except at very beginning)',          /,                 &
           'No. of LU decompositions:             ', i6, /,                 &
           'No. of forward/backward substitutions:', i6, /,                 &
           'No. of singular matrix decompositions:', i6, /,                 &
            /,                                                              &
           'Texit, the time corresponding to the      ',        /,          &
           '       computed Y upon return:            ', f11.4, /,          &
           'Hexit, last accepted step before exit:    ', f11.4, /,          &
           'Hnew, last predicted step (not yet taken):', f11.4 )

    !-----------------------------------------------------------------------
    ! NOTE: The following variables are held THREADPRIVATE and 
    ! therefore do not need to be included in the !$OMP+PRIVATE 
    ! statements below: C, VAR, FIX, RCONST, TIME, TEMP, NUMDEN, 
    ! H2O, PRESS, PHOTOL, HET, and CFACTOR. (bmy, 3/28/16)
    !-----------------------------------------------------------------------
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT  ( SHARED                                                 )&
    !$OMP PRIVATE  ( I,        J,        L,       N,     YLAT               )&
    !$OMP PRIVATE  ( SO4_FRAC, IERR,     RCNTRL,  START, FINISH, ISTATUS    )&
    !$OMP PRIVATE  ( RSTATE,   SpcID,    KppID,   F,     P                  )&
#if defined( MODEL_GEOS )
    !$OMP PRIVATE  ( Vloc,     Aout, OHreact                                )&
#endif
    !$OMP PRIVATE  ( LCH4,     PCO_TOT,  PCO_CH4, PCO_NMVOC                 ) &
    !$OMP REDUCTION( +:ITIM                                                 )&
    !$OMP REDUCTION( +:RTIM                                                 )&
    !$OMP REDUCTION( +:TOTSTEPS                                             )&
    !$OMP REDUCTION( +:TOTFUNCS                                             )&
    !$OMP REDUCTION( +:TOTJACOB                                             )&
    !$OMP REDUCTION( +:TOTACCEP                                             )&
    !$OMP REDUCTION( +:TOTREJEC                                             )&
    !$OMP REDUCTION( +:TOTNUMLU                                             )&
    !$OMP SCHEDULE ( DYNAMIC,  1                                            )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       !====================================================================
       ! For safety sake, initialize certain variables for each grid
       ! box (I,J,L), whether or not chemistry will be done there.
       !====================================================================
       HET       = 0.0_dp            ! Het chem array
       IERR      = 0                 ! Success or failure flag
       ISTATUS   = 0.0_dp            ! Rosenbrock output 
       PHOTOL    = 0.0_dp            ! Photolysis array
       RCNTRL    = 0.0_fp            ! Rosenbrock input
       RSTATE    = 0.0_dp            ! Rosenbrock output
       SO4_FRAC  = 0.0_fp            ! Fraction of SO4 available for photolysis
       P         = 0                 ! GEOS-Chem photolyis species ID

       ! For tagged CO
       LCH4     = 0.0_fp    ! Methane loss rate
       PCO_TOT  = 0.0_fp    ! Total CO production
       PCO_CH4  = 0.0_fp    ! CO production from CH4
       PCO_NMVOC  = 0.0_fp  ! Total CO from NMVOC

       ! Grid-box latitude [degrees]
       YLAT      = GET_YMID( I, J, L )

       ! Temperature [K]
       TEMP      = State_Met%T(I,J,L)

       ! Pressure [hPa]
       PRESS     = GET_PCENTER( I, J, L )

       ! mje Calculate NUMDEN based on ideal gas law (# cm-3)
       NUMDEN    = State_Met%AIRNUMDEN(I,J,L)

       ! mje H2O arrives in g/kg needs to be in mol cm-3 
       H2O       = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)

       !====================================================================
       ! Get photolysis rates (daytime only)
       !
       ! NOTE: The ordering of the photolysis reactions here is
       ! the order in the Fast-J definition file FJX_j2j.dat.
       ! I've assumed that these are the same as in the text files
       ! but this may have been changed.  This needs to be checked
       ! through more thoroughly -- M. Long (3/28/16)
       !
       ! ALSO NOTE: We moved this section above the test to see if grid 
       ! box (I,J,L) is in the chemistry grid.  This will ensure that
       ! J-value diagnostics are defined for all levels in the column.
       ! This modification was validated by a geosfp_4x5_standard
       ! difference test. (bmy, 1/18/18)
       !
       ! Update SUNCOSmid threshold from 0 to cos(98 degrees) since
       ! fast-jx allows for SZA down to 98 degrees. This is important in
       ! the stratosphere-mesosphere where sunlight still illuminates at 
       ! high altitudes if the sun is below the horizon at the surface
       ! (update submitted by E. Fleming (NASA), 10/11/2018)
       !====================================================================
       IF ( State_Met%SUNCOSmid(I,J) > -0.1391731e+0_fp ) THEN

          ! Get the fraction of H2SO4 that is available for photolysis
          ! (this is only valid for UCX-enabled mechanisms)
          IF ( Input_Opt%LUCX ) THEN 
             SO4_FRAC = SO4_PHOTFRAC( I, J, L )
          ELSE
             SO4_FRAC = 0.0_fp
          ENDIF

          ! Adjust certain photolysis rates:
          ! (1) H2SO4 + hv -> SO2 + OH + OH   (UCX-based mechanisms)
          ! (2) O3    + hv -> O2  + O         (UCX-based mechanisms)
          ! (2) O3    + hv -> OH  + OH        (trop-only mechanisms)
          CALL PHOTRATE_ADJ( am_I_root, Input_Opt, State_Diag,               &
                             I,         J,         L,                        &
                             NUMDEN,    TEMP,      H2O,                      &
                             SO4_FRAC,  IERR                                )

          ! Loop over the FAST-JX photolysis species
          DO N = 1, JVN_

             ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
             PHOTOL(N) = ZPJ(L,N,I,J)
             
#if defined( MODEL_GEOS )
             ! Archive in local array
             GLOB_JVAL(I,J,L,N) = PHOTOL(N)
#endif

             !--------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Instantaneous photolysis rates [s-1] (aka J-values)
             ! and noontime photolysis rates [s-1]
             !
             !    NOTE: Attach diagnostics here instead of in module 
             !    fast_jx_mod.F so that we can get the adjusted photolysis
             !    rates (output from routne PHOTRATE_ADJ above).
             !
             ! The mapping between the GEOS-Chem photolysis species and 
             ! the FAST-JX photolysis species is contained in the lookup 
             ! table in input file FJX_j2j.dat.
             ! 
             !    NOTE: Depending on the simulation, some GEOS-Chem 
             !    species might not map to a of the FAST-JX species 
             !    (e.g. SOA species will not be present in a tropchem run).
             !
             ! Some GEOS-Chem photolysis species may have multiple 
             ! branches for photolysis reactions.  These will be 
             ! represented by multiple entries in the FJX_j2j.dat
             ! lookup table.
             !
             !    NOTE: For convenience, we have stored the GEOS-Chem 
             !    photolysis species index (range: 1..State_Chm%nPhotol) 
             !    for each of the FAST-JX photolysis species (range; 
             !    1..JVN_) in the GC_PHOTO_ID array (located in module 
             !    CMN_FJX_MOD.F).
             !
             ! To match the legacy bpch diagnostic, we archive the sum of 
             ! photolysis rates for a given GEOS-Chem species over all of 
             ! the reaction branches.
             !
             !    NOTE: The legacy ND22 bpch diagnostic divides by the
             !    number of times the when grid box (I,J,L) was between
             !    11am and 1pm local solar time.  Because the HISTORY
             !    component can only divide by the number of times the
             !    diagnostic was updated, we counteract this by multiplying
             !    by the factor 86400 [sec/day] / chemistry timestep [sec].
             !--------------------------------------------------------------
             
             ! GC photolysis species index
             P = GC_Photo_Id(N)
             
             ! If this FAST_JX photolysis species maps to a valid 
             ! GEOS-Chem photolysis species (for this simulation)...
             IF ( P > 0 ) THEN

                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Diag%Archive_JVal ) THEN
                   State_Diag%JVal(I,J,L,P) = State_Diag%JVal(I,J,L,P)    &
                                            + PHOTOL(N)
                ENDIF

                ! Archive the noontime photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Diag%Archive_JNoon .and. &
                     State_Met%IsLocalNoon(I,J) ) THEN
                   State_Diag%JNoon(I,J,L,P) = State_Diag%JNoon(I,J,L,P)  &
                                             + ( PHOTOL(N) * JNoon_Fac )
                ENDIF

             ENDIF
          ENDDO
       ENDIF

       !====================================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !====================================================================

       ! If we are not in the troposphere don't do the chemistry!
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Skipping buffer zone (lzh, 08/10/2014)
       IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
          IF ( J <=         Input_Opt%NESTED_J0W ) CYCLE
          IF ( J >  JJPAR - Input_Opt%NESTED_J0E ) CYCLE
          IF ( I <=         Input_Opt%NESTED_I0W ) CYCLE
          IF ( I >  IIPAR - Input_Opt%NESTED_I0E ) CYCLE
       ENDIF

       !====================================================================
       ! Intialize KPP solver arrays: CFACTOR, VAR, FIX, etc.
       !====================================================================
       CALL Init_KPP( )

       !====================================================================
       ! Get rates for heterogeneous chemistry
       !====================================================================

!#if defined( DEVEL )
!       ! Get starting time for rate computation
!       CALL CPU_TIME( start )
!#endif

       IF ( DO_HETCHEM ) THEN

          ! Set hetchem rates
          CALL SET_HET( I, J, L, Input_Opt, State_Chm, State_Met )

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

       ENDIF

       !====================================================================
       ! Initialize species concentrations
       !====================================================================

       ! Loop over KPP Species
       DO N = 1, NSPEC

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
             KppID =  ND65_Kpp_Id(F)

             ! Initialize prod/loss rates
             IF ( KppID > 0 ) C(KppID) = 0.0_dp

          ENDDO

       ENDIF

       IF ( .not. Input_Opt%LUCX ) THEN
          ! Need to copy H2O to the C array for KPP (mps, 4/25/16)
          ! NOTE: H2O is a tracer in UCX and is obtained from State_Chm%Species
          C(ind_H2O) = H2O
       ENDIF

       !==================================================================
       ! Update KPP rates
       !==================================================================

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

#if defined( MODEL_GEOS )
       ! Archive 
       CALL Fun ( VAR, FIX, RCONST, Vloc, Aout=Aout )
       IF ( Input_Opt%NN_RxnRates > 0 ) THEN
          DO N = 1, Input_Opt%NN_RxnRates
             State_Diag%RxnRates(I,J,L,N) = Aout(Input_Opt%RxnRates_IDs(N))
          ENDDO
       ENDIF
#endif

!#if defined( DEVEL )
!       ! Get time when rate computation finished
!       CALL CPU_TIME( finish )
!
!       ! Compute how long it took for KPP to compute rates
!       rtim = rtim + finish - start
!#endif

       !=================================================================
       ! Set options for the KPP Integrator (M. J. Evans)
       !
       ! NOTE: Because RCNTRL(3) is set to an array value that
       ! depends on (I,J,L), we must declare RCNTRL as PRIVATE
       ! within the OpenMP parallel loop and define it just 
       ! before the call to to Integrate. (bmy, 3/24/16)
       !=================================================================

       ! Zero all slots of RCNTRL
       RCNTRL    = 0.0_fp

       ! Starting value for integration time step
       RCNTRL(3) = State_Chm%KPPHvalue(I,J,L)

       !=================================================================
       ! Integrate the box forwards
       !=================================================================

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

#if defined( MODEL_GEOS )
       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
          IF ( ASSOCIATED(State_Diag%KppError) ) THEN
             State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
          ENDIF
       ENDIF
#endif

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
#if defined( MODEL_GEOS )
             IF ( Input_Opt%KppStop ) THEN
                CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
             ! Revert to start values
             ELSE
                VAR = C(1:NVAR)
                FIX = C(NVAR+1:NSPEC)
             ENDIF
             IF ( ASSOCIATED(State_Diag%KppError) ) THEN
                State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
             ENDIF
#else
             CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
#endif
          ENDIF
            
       ENDIF

       ! Copy VAR and FIX back into C (mps, 2/24/16)
       C(1:NVAR)       = VAR(:)
       C(NVAR+1:NSPEC) = FIX(:)

       ! Save for next integration time step
       State_Chm%KPPHvalue(I,J,L) = RSTATE(Nhnew)

       ! Save rate constants in global array (not used)
       GLOB_RCONST(I,J,L,:) = RCONST(:)

       !====================================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species
       !====================================================================
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
       !====================================================================
       ! ND65 (bpch) diagnostic
       !
       ! Obtain prod/loss rates from KPP [molec/cm3]
       !====================================================================
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

             !--------------------------------------------------------
             ! Save out P(CO) and L(CH4) from the fullchem simulation
             ! for use in tagged CO
             !--------------------------------------------------------
             IF ( Input_Opt%DO_SAVE_PCO ) THEN
                IF ( TRIM(FAM_NAMES(F)) == 'PCO'  ) THEN
                   PCO_TOT = VAR(KppID) / DT
                ENDIF
                IF ( TRIM(FAM_NAMES(F)) == 'LCH4' ) THEN
                   LCH4    = VAR(KppID) / DT
                ENDIF
             ENDIF

          ENDDO

          ! For tagged CO, use LCH4 to get P(CO) contributions from
          ! CH4 and NMVOC
          IF ( Input_Opt%DO_SAVE_PCO ) THEN
             ! P(CO)_CH4 is LCH4. Cap so that it is never greater
             ! than total P(CO) to prevent negative P(CO)_NMVOC
             PCO_CH4 = MIN( LCH4, PCO_TOT )
   
             ! P(CO) from NMVOC is the remaining P(CO)
             PCO_NMVOC = PCO_TOT - PCO_CH4
   
             ! Add to AD65 array [molec/cm3/s]
             AD65(I,J,L,NFAM+1) = AD65(I,J,L,NFAM+1) + PCO_CH4
             AD65(I,J,L,NFAM+2) = AD65(I,J,L,NFAM+2) + PCO_NMVOC

          ENDIF

       ENDIF

#if defined( TOMAS )
       !always calculate rate for TOMAS
       DO F = 1, NFAM

          ! Determine dummy species index in KPP
          KppID =  ND65_Kpp_Id(F)

          !-----------------------------------------------------------------
          ! FOR TOMAS MICROPHYSICS:
          !
          ! Obtain P/L with a unit [kg S] for tracing
          ! gas-phase sulfur species production (SO2, SO4, MSA)
          ! (win, 8/4/09)
          !-----------------------------------------------------------------

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

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "VAR" array (in module gckpp_Global.F90).
       !====================================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Loss ) THEN
          DO F = 1, State_Chm%nLoss
             KppID                    = State_Chm%Map_Loss(F)
             State_Diag%Loss(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Prod ) THEN
          DO F = 1, State_Chm%nProd
             KppID                    = State_Chm%Map_Prod(F)
             State_Diag%Prod(I,J,L,F) = VAR(KppID) / DT
          ENDDO
       ENDIF

#if defined( MODEL_GEOS )
       !==============================================================
       ! Write out OH reactivity
       ! The OH reactivity is defined here as the inverse of its life-
       ! time. In a crude ad-hoc approach, manually add all OH reactants
       ! (ckeller, 9/20/2017)
       !==============================================================
       IF ( DoOHreact ) THEN
          CALL Get_OHreactivity ( C, RCONST, OHreact )
          State_Diag%OHreactivity(I,J,L) = OHreact
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

    !=======================================================================
    ! Archive OH, HO2, O1D, O3P concentrations after FlexChem solver
    ! (NOTE: This was formerly done by ohsave.F and diagoh.F)
    !=======================================================================
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

    !=======================================================================
    ! Save quantities for computing mean OH lifetime
    !=======================================================================
    CALL DO_DIAG_OH( State_Met, State_Chm )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after DO_DIAG_OH' )
    ENDIF

#if defined( BPCH_DIAG )
    !=======================================================================
    ! Save out P(O3) and L(O3) for a tagged O3 run
    !
    ! %%%% NOTE: Currently only works when BPCH_DIAG=y %%%%
    !=======================================================================
    IF ( Input_Opt%DO_SAVE_O3 ) THEN
       CALL DIAG20( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FlexChem: after DIAG20' )
       ENDIF
    ENDIF
#endif

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90' )
       RETURN
    ENDIF  

    !=======================================================================
    ! Special handling for UCX chemistry
    !=======================================================================
    IF ( Input_Opt%LUCX ) THEN

       !--------------------------------------------------------------------
       ! If using stratospheric chemistry, applying high-altitude
       ! active nitrogen partitioning and H2SO4 photolysis
       ! approximations  outside the chemgrid
       !--------------------------------------------------------------------
       CALL UCX_NOX( Input_Opt, State_Met, State_Chm )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMDR: after UCX_NOX' )
       ENDIF

       CALL UCX_H2SO4PHOT( Input_Opt, State_Met, State_Chm )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMDR: after UCX_H2SO4PHOT' )
       ENDIF
     
       !--------------------------------------------------------------------
       ! Compute stratospheric chemical tendency for UCX simulations
       !--------------------------------------------------------------------

       ! Loop over advected species
       !$OMP PARALLEL DO               & 
       !$OMP DEFAULT( SHARED         ) &
       !$OMP PRIVATE( I, J, L, N, NA )
       DO NA = 1, State_Chm%nAdvect

          ! Get the species ID from the advected species ID
          N = State_Chm%Map_Advect(NA)

          ! Loop over grid boxes
          DO L = 1, LLPAR
          DO J = 1, JJPAR
          DO I = 1, IIPAR

             ! Aggregate stratospheric chemical tendency [kg box-1]
             ! for tropchem simulations
             IF ( State_Met%InStratosphere(I,J,L) ) THEN
                SChem_Tend(I,J,L,N) = SChem_Tend(I,J,L,N) + &
                     ( State_Chm%Species(I,J,L,N) - Before(I,J,L,N) )
             ENDIF

          ENDDO
          ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

#if defined( MODEL_GEOS )
    ! Archive all needed reaction rates in state_diag
    IF ( Input_Opt%NN_RxnRconst > 0 ) THEN
       DO N = 1, Input_Opt%NN_RxnRconst
          State_Diag%RxnRconst(:,:,:,N) = GLOB_RCONST(:,:,:,Input_Opt%RxnRconst_IDs(N))
       ENDDO
    ENDIF
    IF ( Input_Opt%NN_Jvals > 0 ) THEN
       DO N = 1, Input_Opt%NN_Jvals
          State_Diag%JValIndiv(:,:,:,N) = GLOB_JVAL(:,:,:,Input_Opt%Jval_IDs(N))
       ENDDO
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
    USE CMN_SIZE_Mod
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
#if defined( BPCH_DIAG )
    USE Diag_Mod,       ONLY : AD43
    USE Diag_Mod,       ONLY : LTOH, LTHO2, LTO1D, LTO3P
#endif
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
!  Also note: for the netCDF diagnostics, we have removed multiplication by 
!  LTOH etc arrays.  These are almost always set between 0 and 24.
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

    ! Assume success
    RC = GC_SUCCESS

    ! Point to the array of species concentrations
    AirNumDen => State_Met%AirNumDen
    Spc       => State_Chm%Species

    ! Zero the netCDF diagnostic arrays (if activated) above the 
    ! tropopause or mesopause to avoid having leftover values
    ! from previous timesteps
#if defined( MODEL_GEOS )
    IF ( State_Diag%Archive_O3concAfterChem  ) THEN
       State_Diag%O3concAfterChem  = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_RO2concAfterChem ) THEN
       State_Diag%RO2concAfterChem = 0.0_f4
    ENDIF
#endif
    IF ( State_Diag%Archive_OHconcAfterChem  ) THEN
       State_Diag%OHconcAfterChem  = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_HO2concAfterChem ) THEN
       State_Diag%HO2concAfterChem = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
       State_Diag%O1DconcAfterChem = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
       State_Diag%O3PconcAfterChem = 0.0_f4
    ENDIF

!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J, L ) &
!$OMP SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( .not. State_Met%InChemGrid(I,J,L) ) THEN
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
                             + ( Spc(I,J,L,id_OH) * LTOH(I,J)               )
            ENDIF
#endif

            ! HISTORY (aka netCDF diagnostics)
            IF ( State_Diag%Archive_OHconcAfterChem ) THEN
               State_Diag%OHconcAfterChem(I,J,L) = Spc(I,J,L,id_OH)
            ENDIF
#if defined( MODEL_GEOS )
            IF ( State_Diag%Archive_O3concAfterChem ) THEN
               State_Diag%O3concAfterChem(I,J,L) = Spc(I,J,L,id_O3)
            ENDIF
            IF ( Archive_RO2concAfterChem ) THEN
               IF ( id_A3O2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) +  Spc(I,J,L,id_A3O2)
               IF ( id_ATO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ATO2)
               IF ( id_B3O2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_B3O2)
               IF ( id_BRO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_BRO2)
               IF ( id_DHPCARP > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_DHPCARP)
               IF ( id_DIBOO > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_DIBOO)
               IF ( id_ETO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ETO2)
               IF ( id_HC5OO > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_HC5OO)
               IF ( id_HO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_HO2)
               IF ( id_IEPOXOO > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IEPOXOO)
               IF ( id_INPN > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_INPN)
               IF ( id_ISNOOA > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ISNOOA)
               IF ( id_ISNOOB > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ISNOOB)
               IF ( id_ISNOHOO > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ISNOHOO)
               IF ( id_LIMO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_LIMO2)
               IF ( id_MAOPO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_MAOPO2)
               IF ( id_MO2 > 0  ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_MO2)
               IF ( id_MRO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_MRO2)
               IF ( id_PIO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PIO2)
               IF ( id_PO2 > 0  ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PO2)
               IF ( id_PRNI > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PRNI)
               IF ( id_R4NI > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_R4NI)
               IF ( id_R4O2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_R4O2)
               IF ( id_RIO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_RIO2)
               IF ( id_TRO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_TRO2)
               IF ( id_VRO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_VRO2)
               IF ( id_XRO2 > 0 ) THEN
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_XRO2)
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

            ! HISTORY (aka netCDF diagnostics)
            IF ( State_Diag%Archive_HO2concAfterChem ) THEN
               State_Diag%HO2concAfterChem(I,J,L) = ( Spc(I,J,L,id_HO2)      &
                                                  /   AirNumDen(I,J,L)      )
            ENDIF

         ENDIF

         IF ( Input_Opt%LUCX ) THEN

            !---------------------------------------------------------------
            ! O1D concentration [molec/cm3]
            !---------------------------------------------------------------
            IF ( ok_O1D ) THEN

#if defined( BPCH_DIAG ) 
               ! ND43 (bpch) diagnostic
               IF ( Do_ND43 ) THEN
                  AD43(I,J,L,3) = AD43(I,J,L,3)                              &
                                + ( Spc(I,J,L,id_O1D) * LTO1D(I,J)          )
               ENDIF
#endif

               ! HISTORY (aka netCDF diagnostics)
               IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
                  State_Diag%O1DconcAfterChem(I,J,L) = Spc(I,J,L,id_O1D)
               ENDIF

            ENDIF


            !---------------------------------------------------------------
            ! O3P concentration [molec/cm3]
            !---------------------------------------------------------------
            IF ( ok_O3P ) THEN

#if defined( BPCH_DIAG )
               ! ND43 (bpch) diagnostic
               IF ( Do_ND43 ) THEN
                  AD43(I,J,L,4) = AD43(I,J,L,4)                              &
                                + ( Spc(I,J,L,id_O3P) * LTO3P(I,J)          )
               ENDIF
#endif

               ! HISTORY (aka netCDF diagnostics)
               IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
                  State_Diag%O3PconcAfterChem(I,J,L) = Spc(I,J,L,id_O3P)
               ENDIF

            ENDIF

         ENDIF

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

    !=======================================================================
    ! Init_FlexChem begins here!
    !=======================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Do the following only if it is a full-chemistry simulation
    ! NOTE: If future specialty simulations use the KPP solver, 
    ! modify the IF statement accordingly to allow initialization
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    !=======================================================================
    ! Initialize variables
    !=======================================================================
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_FlexChem (in module GeosCore/flexchem_mod.F90)'
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )

    ! Debug output
    IF ( prtDebug ) THEN
       WRITE( 6, 100 )
100    FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

       WRITE( 6 ,'(a)' ) ' KPP Reaction Reference '
       DO N = 1, NREACT
          WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
       END DO
    ENDIF

    ! Initialize species flags
    id_CH4                   = Ind_( 'CH4', 'A'     ) ! CH4 advected species
    id_HO2                   = Ind_( 'HO2'          )
    id_O3P                   = Ind_( 'O'            )
    id_O1D                   = Ind_( 'O1D'          )
    id_OH                    = Ind_( 'OH'           ) 

#if defined( MODEL_GEOS )
    ! ckeller
    id_O3                    = Ind_( 'O3'           ) 
    id_A3O2                  = Ind_( 'A3O2'         ) 
    id_ATO2                  = Ind_( 'ATO2'         ) 
    id_BRO2                  = Ind_( 'BRO2'         ) 
    id_DHPCARP               = Ind_( 'DHPCARP'      ) 
    id_DIBOO                 = Ind_( 'DIBOO'        ) 
    id_ETO2                  = Ind_( 'ETO2'         ) 
    id_HC5OO                 = Ind_( 'HC5OO'        ) 
    id_IEPOXOO               = Ind_( 'IEPOXOO'      ) 
    id_INPN                  = Ind_( 'INPN'         ) 
    id_ISNOOA                = Ind_( 'ISNOOA'       ) 
    id_ISNOOB                = Ind_( 'ISNOOB'       ) 
    id_ISNOHOO               = Ind_( 'ISNOHOO'      ) 
    id_LIMO2                 = Ind_( 'LIMO2'        ) 
    id_MAOPO2                = Ind_( 'MAOPO2'       ) 
    id_MO2                   = Ind_( 'MO2'          ) 
    id_MRO2                  = Ind_( 'MRO2'         ) 
    id_PIO2                  = Ind_( 'PIO2'         ) 
    id_PO2                   = Ind_( 'PO2'          ) 
    id_PRNI                  = Ind_( 'PRNI'         ) 
    id_R4NI                  = Ind_( 'R4NI'         ) 
    id_R4O2                  = Ind_( 'R4O2'         ) 
    id_RIO2                  = Ind_( 'RIO2'         ) 
    id_TRO2                  = Ind_( 'TRO2'         ) 
    id_VRO2                  = Ind_( 'VRO2'         ) 
    id_XRO2                  = Ind_( 'XRO2'         ) 
#endif

    ! Set flags to denote if each species is defined
    ok_HO2                   = ( id_HO2 > 0         )
    ok_O1D                   = ( id_O1D > 0         )
    ok_O3P                   = ( id_O3P > 0         )
    ok_OH                    = ( id_OH  > 0         )
    
    ! Is the ND43 bpch diagnostic turned on?
    Do_ND43                  = ( Input_Opt%ND43 > 0 ) 

    ! Throw an error if certain diagnostics for UCX are turned on,
    ! but the UCX mechanism is not used in this fullchem simulation
    ! NOTE: Maybe eventually move this error check to state_diag_mod.F90
    IF ( .not. Input_Opt%LUCX ) THEN  
       
       ! O1D diagnostic is only used w/ UCX
       IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
          ErrMsg = 'The "O1DconcAfterChem" diagnostic is turned on ' //      &
                   'but the UCX mechanism is not being used!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! O3P diagnostic is only used w/ UCX
       IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
          ErrMsg = 'The "O3PconcAfterChem" diagnostic is turned on ' //      &
                   'but the UCX mechanism is not being used!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Should we archive OH, HO2, O1D, O3P diagnostics?
#if defined( MODEL_GEOS )
    Do_Diag_OH_HO2_O1D_O3P      = ( Do_ND43                             .or. &  
                                    State_Diag%Archive_O3concAfterChem  .or. &
                                    State_Diag%Archive_RO2concAfterChem .or. &
                                    State_Diag%Archive_OHconcAfterChem  .or. &
                                    State_Diag%Archive_HO2concAfterChem .or. &
                                    State_Diag%Archive_O1DconcAfterChem .or. &
                                    State_Diag%Archive_O3PconcAfterChem     )
#else
    Do_Diag_OH_HO2_O1D_O3P      = ( Do_ND43                             .or. &  
                                    State_Diag%Archive_OHconcAfterChem  .or. &
                                    State_Diag%Archive_HO2concAfterChem .or. &
                                    State_Diag%Archive_O1DconcAfterChem .or. &
                                    State_Diag%Archive_O3PconcAfterChem     )
#endif

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    !--------------------------------------------------------------------
    ! Pre-store the KPP indices for each KPP prod/loss species or family
    !--------------------------------------------------------------------
    IF ( nFam > 0 ) THEN
             
       ! Allocate mapping array for KPP Ids for ND65 bpch diagnostic
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
  SUBROUTINE Cleanup_FlexChem( am_I_Root, RC )  
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!    
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!    
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!    
! !REVISION HISTORY:
!  24 Aug 2016 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Cleanup_FlexChem begins here!
    !=================================================================

    ! INitialize
    RC = GC_SUCCESS

    IF ( ALLOCATED( JvCountDay ) ) THEN
       DEALLOCATE( JvCountDay, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumDay ) ) THEN
       DEALLOCATE( JvSumDay, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvCountMon ) ) THEN
       DEALLOCATE( JvCountMon, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumMon ) ) THEN
       DEALLOCATE( JvSumMon, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Cleanup_FlexChem
!EOC
END MODULE FlexChem_Mod
