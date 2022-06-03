!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: isorropiaii_mod.F90
!
! !DESCRIPTION: Module ISORROPIAII\_MOD contains the routines that provide
!  the interface between ISORROPIA II and GEOS-Chem.
!\\
!\\
!  The actual ISORROPIA II code which performs Na-SO4-NH3-NO3-Cl-(Ca-K-Mg)
!  aerosol thermodynamic equilibrium is in \texttt{isorropiaIIcode.f}.
!\\
!\\
! !INTERFACE:
!
MODULE ISORROPIAII_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CLEANUP_ISORROPIAII
  PUBLIC  :: DO_ISORROPIAII
  PUBLIC  :: GET_GNO3
#if defined( MODEL_CESM )
  PUBLIC  :: INIT_ISORROPIAII
#else
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: INIT_ISORROPIAII
#endif
  PRIVATE :: SAFELOG10
  PRIVATE :: SET_HNO3
!
! !REMARKS:
!  Original Author:
!  *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
!  *** GEORGIA INSTITUTE OF TECHNOLOGY
!  *** WRITTEN BY ATHANASIOS NENES
!  *** UPDATED BY CHRISTOS FOUNTOUKIS
!                                                                             .
!  Original v1.3 isorropia implementation into GEOS-Chem by
!  Becky Alexander and Bob Yantosca (bec, bmy, 4/12/05, 11/2/05)
!                                                                             .
!  For Ca,K,Mg = 0, ISORROPIA II performs exactly like ISORROPIAv1.7
!  Ca, K, Mg, Na from dust is not currently considered
!                                                                             .
!  To implement ISORROPIA II into GEOS-Chem:
!    * cleanup_isorropiaII needs to be called from cleanup.f
!    * DO_ISORROPIA needs to be replaced with DO_ISORROPIAII in chemistry_mod.f
!    * Change ISORROPIA to ISORROPIAII in sulfate_mod.f
!    * add isorropiaII_mod.f, isorropiaIIcode.f, and irspia.inc to Makefile
!                                                                             .
!  ISORROPIA II implementation notes by Havala O.T. Pye:
!  (1) The original isorropia code from T.Nenes is left as unmodified as
!       possible. Original isorropia code can be found in isorropiaIIcode.f
!       and common blocks can be found in isrpia.inc. For future upgrades
!       to isorropia, replace isrpia.inc and isorropiaIIcode.f with the new
!       version of isorropia and modify the call to ISORROPIA in this module.
!       Please let the original author know of any changes made to ISORROPIA.
!  (2) As of Nov 2007, routines using non-zero Ca, K, and Mg do not always
!       conserve mass. Ca, K, and Mg are set to zero.
!                                                                             .
!  NOTE: ISORROPIA is Greek for "equilibrium", in case you were wondering.
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  ! Array for offline HNO3 (for relaxation of M.M.)
  REAL(fp), ALLOCATABLE :: HNO3_sav(:,:,:)

  ! Array for offline use in sulfate_mod (SEASALT_CHEM)
  REAL(fp), ALLOCATABLE :: GAS_HNO3(:,:,:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Add a C-preprocessor switch to skip calling ISORROPIA if the pressure
!%%% and/or temperature lie outside of the range that will produce a stable
!%%% solution.  This will eliminate the random noise observed in the
!%%% ISORROPIA output.
!%%%
!%%% Leaving this feature deactivated will replicate the prior behavior in
!%%% v11-01 and earlier GEOS-Chem versions.  This will become the default
!%%% setting in a future version, but give the user the choice to activate
!%%% or deactivate this for now.
!%%%
!%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
!%%%
!#define SKIP_IF_P_AND_T_ARE_OUT_OF_RANGE 1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_isorropiaii
!
! !DESCRIPTION: Subroutine DO\_ISORROPIAII is the interface between the
!  GEOS-Chem model and the aerosol thermodynamical equilibrium routine
!  ISORROPIA II.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_ISORROPIAII( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,         ONLY : NDUST
    USE ErrCode_Mod
    USE ERROR_MOD,            ONLY : DEBUG_MSG
    USE ERROR_MOD,            ONLY : ERROR_STOP
    USE ERROR_MOD,            ONLY : SAFE_DIV
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AIRMW, PI
    USE PhysConstants,        ONLY : PI
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_MONTH
    USE TIME_MOD,             ONLY : ITS_A_NEW_MONTH
    USE TIME_MOD,             ONLY : GET_ELAPSED_SEC
    USE IsorropiaII_Main_Mod, ONLY : Isorropia
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
! !REMARKS:
!  Original isorropia v1.3 implmentation: (rjp, bec, bmy, 12/17/01, 8/22/05)
!
! !REVISION HISTORY:
!  24 Aug 2007 - H. O. T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Array dimensions
    INTEGER, PARAMETER       :: NOTHERA  =  9
    INTEGER, PARAMETER       :: NCTRLA   =  2
    INTEGER, PARAMETER       :: NCOMPA   =  8
    INTEGER, PARAMETER       :: NIONSA   = 10
    INTEGER, PARAMETER       :: NGASAQA  =  3
    INTEGER, PARAMETER       :: NSLDSA   = 19

    ! Concentration lower limit [mole/m3]
    REAL(fp),  PARAMETER       :: CONMIN = 1.0e-30_fp
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE            :: FIRST = .TRUE.
    LOGICAL, SAVE            :: IS_HMS = .FALSE.
    INTEGER, SAVE            :: id_HNO3, id_NH3,  id_NH4
    INTEGER, SAVE            :: id_NIT,  id_SALA, id_SO4
    INTEGER, SAVE            :: id_HMS ! jmm 12/5/18
    INTEGER, SAVE            :: id_SALACL, id_HCL, id_SALCCL
    INTEGER, SAVE            :: id_SO4s, id_NITs, id_SALC
    INTEGER, SAVE            :: id_SALAAL, id_SALCAL

    ! Scalars
    INTEGER                  :: I,    J,    L,    N,  NM
    REAL(fp)                 :: ANO3, GNO3, ACL, GCL
    REAL(f8)                 :: RHI,  TEMPI, P_Pa
    REAL(fp)                 :: TCA,  TMG,  TK,   HNO3_DEN
    REAL(fp)                 :: TNA,  TCL,  TNH3, TNH4
    REAL(fp)                 :: TNIT, TNO3, TSO4, VOL
    REAL(fp)                 :: HNO3_UGM3
    REAL(f8)                 :: AERLIQ(NIONSA+NGASAQA+2)
    REAL(f8)                 :: AERSLD(NSLDSA)
    REAL(f8)                 :: GAS(NGASAQA)
    REAL(f8)                 :: OTHER(NOTHERA)
    REAL(f8)                 :: WI(NCOMPA)
    REAL(f8)                 :: WT(NCOMPA)
    REAL(f8)                 :: CNTRL(NCTRLA)
    REAL(f8)                 :: AlkR !Alkalinity % depleted
    REAL(f8)                 :: Qk, PHCl, F_HCl, F_HNO3
    REAL(f8)                 :: Hplus !H+ in SALC,mol/m3
    REAL(f8)                 :: Dcs !SALC diameter, m
    REAL(f8)                 :: n_air !air density, molec/cm3
    REAL(f8)                 :: n_ssc !SALC number concentration, molec/m3

    ! Strings
    CHARACTER(LEN=15)        :: SCASI
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=255)       :: X
    !Temporary variables to check if division is safe
    REAL(fp)                 :: NUM_SAV, DEN_SAV

    ! AEROPH: Temporary variable for pH (hotp 8/11/09)
    REAL(fp)                 :: HPLUSTEMP
    ! Temporary variable for SO4--
    REAL(fp)                 :: SULFTEMP
    ! Temporary variable for HSO4-
    REAL(fp)                 :: BISULTEMP
    ! Temporary variable for NO3-
    REAL(fp)                 :: NITRTEMP
    ! Temporary variable for Cl-
    REAL(fp)                 :: CLTEMP

    ! debug variables
    INTEGER                  :: Itemp, Jtemp, Ltemp
    LOGICAL, SAVE            :: FIRSTCHECK = .TRUE.

    LOGICAL                  :: IT_IS_AN_AEROSOL_SIM
    LOGICAL                  :: IT_IS_A_FULLCHEM_SIM
    LOGICAL, SAVE            :: USE_HNO3_FROM_HEMCO = .FALSE.
    LOGICAL, SAVE            :: USE_HCl_FROM_HEMCO  = .FALSE.
    LOGICAL                  :: prtDebug

    ! Pointers
    TYPE(SpcConc), POINTER   :: Spc(:)

    ! Are we out of the range of valid inputs?
    Logical                  :: OutOfBounds

    ! Local array for HNO3, HCl from HEMCO
    REAL(fp) :: OFFLINE_HNO3(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp) :: OFFLINE_HCl (State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    !=================================================================
    ! DO_ISORROPIAII begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_ISORROPIAII (in module GeosCore/isorropiaII_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    prtDebug             = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Zero State_Chm arrays to avoid leftover values from hanging
    ! around between calls -- especially up near the tropopause
    State_Chm%IsorropAeropH    = 0.0_fp
    State_Chm%IsorropHplus     = 0.0_fp
    State_Chm%IsorropAeroH2O   = 0.0_fp
    State_Chm%IsorropSulfate   = 0.0_fp
    State_Chm%IsorropNitrate   = 0.0_fp
    State_Chm%IsorropBisulfate = 0.0_fp
    State_Chm%IsorropChloride  = 0.0_fp

    ! First-time initialization
    IF ( FIRST ) THEN

       ! Make sure certain tracers are defined
       id_HNO3   = Ind_('HNO3'  )
       id_NH3    = Ind_('NH3'   )
       id_NH4    = Ind_('NH4'   )
       id_NIT    = Ind_('NIT'   )
       id_SALA   = Ind_('SALA'  )
       id_SO4    = Ind_('SO4'   )
       id_HMS    = Ind_('HMS'   )
       id_SALACL = Ind_('SALACL')
       id_HCL    = Ind_('HCl'   )
       id_SALC   = Ind_('SALC'  )
       id_SALCCL = Ind_('SALCCL')
       !id_NH4s  = Ind_('NH4s'  )
       id_NITs   = Ind_('NITs'  )
       id_SO4s   = Ind_('SO4s'  )
       id_SALAAL = Ind_('SALAAL')
       id_SALCAL = Ind_('SALCAL')

       ! Set a flag if HMS is defined
       IS_HMS    = ( id_HMS > 0 )

       ! Make sure certain tracers are defined
       IF ( id_SO4 <= 0 ) THEN
          ErrMsg = 'SO4 is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_HMS <= 0 .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ErrMsg = 'HMS is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_NH3 <= 0 ) THEN
          ErrMsg = 'NH3 is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_NH4 <= 0 ) THEN
          ErrMsg = 'NH4 is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_NIT <= 0 ) THEN
          ErrMsg = 'NIT is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALA <= 0 ) THEN
          ErrMsg = 'SALA is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALACL <= 0 ) THEN
          ErrMsg = 'SALACL is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALC <= 0 ) THEN
          ErrMsg = 'SALC is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALCCL <= 0 ) THEN
          ErrMsg = 'SALCCL is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       !IF ( id_NH4s <= 0 ) THEN
       !   ErrMsg = 'NH4s is an undefined species!'
       !   CALL GC_Error( ErrMsg, RC, ThisLoc )
       !   RETURN
       !ENDIF
       IF ( id_NITs <= 0 ) THEN
          ErrMsg = 'NITs is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SO4s <= 0 ) THEN
          ErrMsg = 'SO4s is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALAAL <= 0 ) THEN
          ErrMsg = 'SALACL is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( id_SALCAL <= 0 ) THEN
          ErrMsg = 'SALCAL is an undefined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

#if !defined( MODEL_CESM )
       ! Initialize arrays
       CALL INIT_ISORROPIAII( State_Grid )
#endif

       ! Check to see if we need to get HNO3 from HEMCO
       IF ( id_HNO3 <= 0 ) THEN

          IF ( IT_IS_A_FULLCHEM_SIM ) THEN

             ! Coupled simulation: stop w/ error since we need HNO3
             ErrMsg = 'HNO3 is an undefined species!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

             ! Offline simulation: get HNO3 from HEMCO (mps, 9/23/14)
             USE_HNO3_FROM_HEMCO = .TRUE.

          ELSE

             ! ISORROPIA is only valid for full-chem or aerosol-only sims
             ErrMsg = 'Invalid simulation type!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ENDIF
       ENDIF

       ! Check to see if we need to get HCl from HEMCO
       IF ( id_HCL <= 0 ) THEN

          IF ( IT_IS_A_FULLCHEM_SIM ) THEN

             ! Coupled simulation: stop w/ error since we need HCl
             ErrMsg = 'HCl is an undefined species!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

             ! Offline simulation: get HCl from HEMCO (mps, 6/11/2020)
             USE_HCl_FROM_HEMCO = .TRUE.

          ELSE

             ! ISORROPIA is only valid for full-chem or aerosol-only sims
             ErrMsg = 'Invalid simulation type!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ENDIF
       ENDIF

       ! Print out
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) REPEAT( '=', 79 )
          WRITE( 6, 110 )
          WRITE( 6, 100 ) REPEAT( '=', 79 )
       ENDIF

100    FORMAT( a                                            )
110    FORMAT( 'Successfully initialized ISORROPIA code II' )

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    ! Initialize for each timestep (bec, bmy, 4/15/05)
    IF ( IT_IS_AN_AEROSOL_SIM ) THEN
       GAS_HNO3 = 0.0_fp
    ENDIF

    ! Evaluate offline global HNO3 from HEMCO is using. Doing this every
    ! timestep allows usage of HEMCO's scaling and masking functionality
    IF ( USE_HNO3_FROM_HEMCO ) THEN
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_HNO3', OFFLINE_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'GLOBAL_HNO3 not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Evaluate offline global HCl from HEMCO is using. Doing this every
    ! timestep allows usage of HEMCO's scaling and masking functionality
    IF ( USE_HCl_FROM_HEMCO ) THEN
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_HCl', OFFLINE_HCl, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'GLOBAL_HCl not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    !========================================================================
    ! Loop over grid boxes and call ISORROPIA (see comments in the
    ! ISORROPIA routine ISORROPIAIICODE.f which describes
    ! the input/output args)
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                  ) &
    !$OMP PRIVATE( I,        J,         L,           N,         WI         ) &
    !$OMP PRIVATE( WT,       GAS,       TEMPI,       RHI,       VOL        ) &
    !$OMP PRIVATE( TSO4,     TNH3,      TNA,         TCL,       ANO3       ) &
    !$OMP PRIVATE( GNO3,     TCA,       TMG,         TK,        CNTRL      ) &
    !$OMP PRIVATE( SCASI,    P_Pa,      TNO3,        AERLIQ,    AERSLD     ) &
    !$OMP PRIVATE( OTHER,    TNH4,      TNIT,        HPLUSTEMP, NUM_SAV    ) &
    !$OMP PRIVATE( GCL,      ACL,       AlkR,        NM,        PHCl       ) &
    !$OMP PRIVATE( Qk,       n_air,     n_ssc,       Hplus,     Dcs        ) &
    !$OMP PRIVATE( DEN_SAV,  HNO3_DEN,  OutOfBounds, F_HCL,     F_HNO3     ) &
    !$OMP PRIVATE( SULFTEMP, BISULTEMP, NITRTEMP,    HNO3_UGM3, CLTEMP     ) &
    !$OMP SCHEDULE( DYNAMIC, 1                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Only applying ISORROPIA II in troposphere
       IF ( State_Met%InStratMeso(I,J,L) ) CYCLE

       ! Zero PRIVATE variables
       Dcs   = 0.0_fp
       HPlus = 0.0_fp
       n_air = 0.0_fp
       n_ssc = 0.0_fp
       PHCl  = 0.0_fp

       ! Initialize PRIVATE meteorological variables
       TEMPI = State_Met%T(I,J,L)                    ! Temperature [K]
       P_Pa  = State_Met%PMid(I,J,L)*100.0e+0_f8     ! Pressure [Pa]
       RHI   = State_Met%RH(I,J,L) * 1.e-2_fp        ! Rel Humidity [1]
       RHI   = MAX( 0.01e+0_fp, RHI )                !  force into range
       RHI   = MIN( 0.98e+0_fp, RHI )                !  0.01 < RH <= 0.98
       VOL   = State_Met%AIRVOL(I,J,L)               ! Grid box volume [m3]

       !-----------------------------------------------------------
       ! Now include full themodynamics for SSA, xnw 11/20/17
       ! Assume coarse SSA is externally mixed with fine aerosols,
       ! fine mode aerosols first reach equilibrium then coarse SSA
       !-----------------------------------------------------------

       ! This is for saving PHCl
       IF ( id_HCl > 0 ) THEN
          PHCl = Spc(id_HCl)%Conc(I,J,L) !initial HCl, kg
       ENDIF

       ! Only do coarse mode cacluation when SALC exists
       ! SET NM = 1 to skip coarse SSA thermodynamic
       IF (Spc(id_SALC)%Conc(I,J,L) .GT. 1.e-20) THEN
          NM = 2
          ! Calculate parameters for mass transfer
          Dcs = State_Chm%AeroRadi(I,J,L,12) * 2.0e-2_fp !(cm->m)
          n_air = State_Met%AIRNUMDEN(I,J,L) !(molec/cm3)
          n_ssc = Spc(id_SALC)%Conc(I,J,L) / 2.2e+3_fp / &
                  ((1.0_fp/6.0_fp) * Pi * Dcs**3 ) / VOL !(kg->#/m3)
       ELSE
          NM = 1
       ENDIF

       !---------------------------------------------------------------------
       ! 1 for fine, 2 for coarse
       ! The coarse SSA thermodynamic does not work for offline simulation
       !---------------------------------------------------------------------
       DO N = 1, NM

          ! Zero PRIVATE variables that get assigned in this loop over NM.
          ! This will prevent values from prior iterations hanging around.
          ACl       = 0.0_fp
          AERLIQ    = 0.0_f8
          AERSLD    = 0.0_f8
          AlkR      = 0.0_fp
          ANO3      = 0.0_fp
          F_HNO3    = 0.0_fp
          F_HCl     = 0.0_fp
          GAS       = 0.0_fp
          GCl       = 0.0_fp
          GNO3      = 0.0_fp
          HNO3_UGM3 = 0.0_fp
          OTHER     = 0.0_f8
          Qk        = 0.0_f8
          SCASI     = ''
          TCa       = 0.0_fp
          TCl       = 0.0_fp
          TK        = 0.0_fp
          TMg       = 0.0_fp
          TNa       = 0.0_fp
          TNH3      = 0.0_fp
          TNH4      = 0.0_fp
          TNIT      = 0.0_fp
          TNO3      = 0.0_fp
          TSO4      = 0.0_fp
          WI        = 0.0_fp
          WT        = 0.0_fp

          !-----------------------------------------------
          ! Compute Alkalinity % consumed in the grid box
          !-----------------------------------------------

          IF (N == 1) THEN
             IF (Spc(id_SALAAL)%Conc(I,J,L) .GT. CONMIN) THEN
                AlkR = Spc(id_SALAAL)%Conc(I,J,L) / Spc(id_SALA)%Conc(I,J,L)
                AlkR = MAX( (1.0_fp-AlkR), CONMIN)
             ELSE
                AlkR = 1.0_fp
             ENDIF
          ELSE
             IF (Spc(id_SALCAL)%Conc(I,J,L) .GT. CONMIN) THEN
                AlkR = Spc(id_SALCAL)%Conc(I,J,L) / Spc(id_SALC)%Conc(I,J,L)
                AlkR = MAX( (1.0_fp-AlkR), CONMIN)
             ELSE
                AlkR = 1.0_fp
             ENDIF
          ENDIF

          !--------------------------------
          ! Compute quantities for ISORROPIA
          !---------------------------------

          IF ( N == 1 ) THEN

             ! Total SO4 [mole/m3], also consider SO4s in SALA
             IF ( IS_HMS ) THEN
                TSO4 = (Spc(id_SO4)%Conc(I,J,L)+Spc(id_SALA)%Conc(I,J,L)*0.08_fp*AlkR) * &
                     1.e+3_fp / ( 96.0_fp * VOL)                             &
                     + Spc(id_HMS)%Conc(I,J,L) * 0.5e+3_fp / ( 111.0_fp * VOL )
             ELSE
                TSO4 = (Spc(id_SO4)%Conc(I,J,L)+Spc(id_SALA)%Conc(I,J,L)*0.08_fp*AlkR) * &
                     1.e+3_fp / ( 96.0_fp * VOL)
             ENDIF


             ! Total NH3 [mole/m3]
             TNH3 = Spc(id_NH4)%Conc(I,J,L) * 1.0e+3_fp / (18.0_fp * VOL) + &
                    Spc(id_NH3)%Conc(I,J,L) * 1.0e+3_fp / (17.0_fp * VOL)

          ELSE

             ! Total SO4 [mole/m3], also consider SO4s in SALC
             TSO4 = Spc(id_SO4s)%Conc(I,J,L) *                               &
                    1.e+3_fp * AlkR / (31.4_fp * VOL) +                      &
                    Spc(id_SALC)%Conc(I,J,L) * 0.08_fp *                     &
                    1.e+3_fp * AlkR / (96.0_fp * VOL)

             ! Total NH3 [mole/m3]
             !TNH3 = Spc(id_NH4s)%Conc(I,J,L)*1.e+3_fp*AlkR/(31.4e+0_fp*VOL)+ &
             !       Spc(id_NH3)%Conc(I,J,L) * 1.e+3_fp / (17.e+0_fp * VOL)
             TNH3  = 0.0_fp

          ENDIF

          IF (N == 1) THEN
             ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
             !TNA = Spc(id_SALA)%Conc(I,J,L) * 0.3061e+0_fp * 1.e+3_fp / &
             !      ( 22.99e+0_fp  * VOL  )

             ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
             ! increase to account for all cations, xnw 11/26/17
             TNA = Spc(id_SALA)%Conc(I,J,L) * 0.397_fp * 1.0e+3_fp           &
                   * AlkR / ( 23.0_fp  * VOL  )

             ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
             !TCL = Spc(id_SALA)%Conc(I,J,L) * 0.5504e+0_fp * 1.e+3_fp / &
             !      ( 35.45e+0_fp  * VOL  )

             ! track chloride in sea salt correctly, xnw 10/12/17
             ! Aerosol phase Cl-, [mole/m3]
             ACL = Spc(id_SALACL)%Conc(I,J,L) * 1.0e+3_fp * AlkR /           &
                   ( 35.45_fp  * VOL  )
          ELSE

             TNA = Spc(id_SALC)%Conc(I,J,L) * 0.378_fp * 1.0e+3_fp           &
                   * AlkR / ( 23.0_fp  * VOL  )
             ACL = Spc(id_SALCCL)%Conc(I,J,L) * 1.0e+3_fp * AlkR /           &
                   ( 35.45_fp  * VOL  )

          ENDIF

          ! Gas phase Cl-, [mole/m3]
          IF ( id_HCl > 0 ) THEN
             GCL = Spc(id_HCL)%Conc(I,J,L) * 1.0e+3_fp /(36.45_fp * VOL)
          ELSE
             ! HCl is in v/v (from HEMCO)
             GCL = OFFLINE_HCl(I,J,L) / VOL
          ENDIF

          ! Total Cl- [mole/m3]
          TCL = ACL + GCL

          !====================================================================
          ! NOTE: As of 11/2007, ISORROPIAII does not conserve mass when Ca,K,Mg
          ! are non-zero. If you would like to consider Ca, K, Mg from seasalt
          ! and dust, isorropiaIIcode.f ISRP4F routines must be debugged.
          ! (hotp, bmy, 2/1/10)
          !
          ! ! Total Ca2+ (1.16% by weight of seasalt) [mole/m3]
          ! TCA      = Spc(id_SALA)%Conc(I,J,L) * 0.0116e+0_fp * 1.d3 /
          !                            ( 40.08e+0_fp  * VOL  )
          !
          ! ! Total K+   (1.1% by weight of seasalt)  [mole/m3]
          ! TK       = Spc(id_SALA)%Conc(I,J,L) * 0.0110e+0_fp * 1.d3 /
          !                            ( 39.102e+0_fp * VOL  )
          !
          ! ! Total Mg+  (3.69% by weight of seasalt) [mole/m3]
          ! TMG      = Spc(id_SALA)%Conc(I,J,L) * 0.0369e+0_fp * 1.d3 /
          !                            ( 24.312e+0_fp * VOL  )
          !====================================================================
          !! Set Ca, K, Mg to zero for time being (hotp, bmy, 2/1/10)
          !TCA      = 0e+0_fp
          !TK       = 0e+0_fp
          !TMG      = 0e+0_fp

          ! Compute gas-phase NO3
          IF ( id_HNO3 > 0 ) THEN

             !---------------------
             ! COUPLED SIMULATION
             !---------------------

             ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
             GNO3 = Spc(id_HNO3)%Conc(I,J,L)
             GNO3 = MAX( GNO3 * 1.e+3_fp / ( 63.0_fp * VOL ), CONMIN )

             ! Aerosol-phase NO3 [mole/m3]
             IF (N == 1) THEN
                ANO3 = Spc(id_NIT)%Conc(I,J,L) * 1.e+3_fp / (62.0_fp * VOL)
             ELSE
                ANO3 = Spc(id_NITs)%Conc(I,J,L) * 1.e+3_fp * AlkR / (31.4_fp * VOL)
             ENDIF

             ! Total NO3 [mole/m3]
             TNO3    = GNO3 + ANO3

          ELSE

             !---------------------
             ! OFFLINE SIMULATION
             !---------------------

             ! Relax to monthly mean HNO3 concentrations every 3 hours
             ! Otherwise just use the concentration in HNO3_sav
             IF ( MOD( GET_ELAPSED_SEC(), 10800 ) == 0 ) THEN
                ! HNO3 is in v/v (from HEMCO), convert to ug/m3
                HNO3_UGM3 = OFFLINE_HNO3(I,J,L) * State_Met%AIRDEN(I,J,L) &
                            * 1.e+9_fp / ( AIRMW / 63.0_fp )
             ELSE
                HNO3_UGM3 = HNO3_sav(I,J,L)
             ENDIF

             ! Convert total inorganic NO3 from [ug/m3] to [mole/m3].
             TNO3  = HNO3_UGM3 * 1.e-6_fp / 63.0_fp

             ANO3 = 0.0_fp
             GNO3 = TNO3

          ENDIF

          !---------------------------------
          ! Call ISORROPIA
          !---------------------------------

          ! set type of ISORROPIA call
          ! Forward problem, do not change this value
          ! 0e+0_fp represents forward problem
          CNTRL(1) = 0.0_fp

          ! Metastable for now
          ! 1e+0_fp represents metastable problem
          CNTRL(2) = 1.0_fp

          ! Insert concentrations [mole/m3] into WI & prevent underflow
          WI(1)    = MAX( TNA,  CONMIN )
          WI(2)    = MAX( TSO4, CONMIN )
          WI(3)    = MAX( TNH3, CONMIN )
          WI(4)    = MAX( TNO3, CONMIN )
          WI(5)    = MAX( TCL,  CONMIN )
          WI(6)    = MAX( TCA,  CONMIN )
          WI(7)    = MAX( TK,   CONMIN )
          WI(8)    = MAX( TMG,  CONMIN )

#if defined( SKIP_IF_P_AND_T_ARE_OUT_OF_RANGE )

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%% If the C-preprocessor switch is activated then check if
          !%%% pressure and temperature are in the range that will result
          !%%% in a stable solution.  If not, then we will skip calling
          !%%% ISORROPIA to avoid random noise in the output.
          !%%%
          !%%% NOTE: Turning this feature on will result in differences
          !%%% with respect to prior GEOS-Chem versions.  So we'll give
          !%%% the user the option to activate it or not.  At some point
          !%%% in the future this will become the default setting.
          !%%%
          !%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! SDE 2017-01-18: Recommendation from Shannon Capps
          ! Skip equilibrium if T < 250 K or P < 200 hPa
          OutOfBounds = ((P_Pa.lt.200.0_f8).or.(TEMPI.lt.250.0e+0_f8))

#else

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%% Always call ISORROPIA, regardless of the values of pressure
          !%%% and temperature.  This will match the prior behavior of
          !%%% when comparing to v11-01 and earlier versions.
          !%%%
          !%%%  -- Seb Eastham and Bob Yantosca (1/25/17)
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ! Never skip calling ISORROPIA
          OutOfBounds = .FALSE.

#endif

          IF ( OutOfBounds ) THEN

             ! %%% Skip equilibrium %%%
             ! Just keep things the way they are
             ! and spoof the other outputs
             WT     = WI
             AERLIQ = 0.0_f8
             AERSLD = 0.0_f8
             GAS    = 0.0_f8
             OTHER  = 0.0_f8

             IF (N == 1) THEN
                ANO3 = Spc(id_NIT)%Conc(I,J,L)  * 1.e+3_fp  / (62.0_fp * VOL)
             ELSE
                ANO3 = Spc(id_NITs)%Conc(I,J,L) * 1.e+3_fp * AlkR / (31.4_fp * VOL)
             ENDIF

          ELSE

             ! %%% Perform aerosol thermodynamic equilibrium %%%
             ! ISORROPIA can be found in ISORROPIAIICODE.F
             ! inputs are WI, RHI, TEMPI, CNTRL
             CALL ISORROPIA( WI,    RHI,  TEMPI,  CNTRL, &
                             WT,    GAS,  AERLIQ, AERSLD, &
                             SCASI, OTHER                 )

             ! Consider mass transfer and acid limitation for coarse
             ! mode calculation
             IF ( N == 2 ) THEN
                !Hplus = AERLIQ(1) !H+ in aerosol, mol/m3
                Hplus = 1.e-5_fp*18.e-3_fp*AERLIQ(8)
                CALL GET_QK(GNO3, GCL, GAS(2), GAS(3), Hplus, Dcs, TEMPI, &
                            n_air, n_ssc, Qk, F_HNO3, F_HCl)
                GAS(2) = GAS(2) * Qk
                GAS(3) = GAS(3) * Qk
                GAS(2) = GNO3 - (GNO3 - GAS(2)) * (1.e0_fp - exp(-F_HNO3))
                GAS(3) = GCL  - (GCL  - GAS(3)) * (1.e0_fp - exp(-F_HCl ))
             ENDIF

             ! Retrieve concentrations in mol/m3
             TSO4 = WT(2)
             TNH3 = GAS(1)
             TNH4 = WT(3) - GAS(1)
             GNO3 = GAS(2)
             TNO3 = WT(4)
             ANO3 = TNO3 - GNO3
             GCL  = GAS(3)
             TCL  = WT(5)
             ACL  = TCL - GCL

          ENDIF

          !---------------------------------
          ! Save back into tracer array
          !---------------------------------
          ! Convert ISORROPIA output from [mole/m3] to [kg]
          IF ( N == 1 ) THEN
             TSO4 = MAX( 96.e-3_fp   * VOL * TSO4, CONMIN )
             TNH4 = MAX( 18.e-3_fp   * VOL * TNH4, CONMIN )
             TNIT = MAX( 62.e-3_fp   * VOL * ANO3, CONMIN )
             ACL  = MAX( 35.45e-3_fp * VOL * ACL,  CONMIN )
             TNH3 = MAX( 17.e-3_fp   * VOL * TNH3, CONMIN )
          ELSE
             TSO4 = MAX( 31.4e-3_fp  * VOL * TSO4, CONMIN )
            !TNH4 = MAX( 31.4e-3_fp  * VOL * TNH4, CONMIN )
             TNIT = MAX( 31.4e-3_fp  * VOL * ANO3, CONMIN )
             ACL  = MAX( 35.45e-3_fp * VOL * ACL,  CONMIN )
          ENDIF
          !TNH3 = MAX( 17.e-3_fp * VOL * TNH3, CONMIN )
          GCL  = MAX( 36.45e-3_fp * VOL * GCL, CONMIN )

          ! Save tracers back into Spc array [kg]
          ! no longer save TSO4 back into Spc. SO4 is all aerosol phase
          ! (hotp 11/7/07)
          ! Spc(id_SO4)%Conc(I,J,L) = TSO4
          !Spc(id_NH3)%Conc(I,J,L) = TNH3
          IF ( id_HCl > 0 ) THEN
             Spc(id_HCL)%Conc(I,J,L) = GCL
          ENDIF

          IF (N == 1) THEN
             Spc(id_NH3   )%Conc(I,J,L) = TNH3
             Spc(id_NH4   )%Conc(I,J,L) = TNH4
             Spc(id_NIT   )%Conc(I,J,L) = TNIT
             Spc(id_SALACL)%Conc(I,J,L) = &
                            Spc(id_SALACL)%Conc(I,J,L)*(1.0_fp-AlkR) + ACL
          ELSE
            !Spc(id_NH4s  )%Conc(I,J,L) = &
            !               Spc(id_NH4s  )%Conc(I,J,L) * (1.0_fp-AlkR) + TNH4
             Spc(id_NITs  )%Conc(I,J,L) = &
                            Spc(id_NITs  )%Conc(I,J,L) * (1.0_fp-AlkR) + TNIT
             Spc(id_SALCCL)%Conc(I,J,L) = &
                            Spc(id_SALCCL)%Conc(I,J,L) * (1.0_fp-AlkR) + ACL
          ENDIF

          ! Special handling for HNO3 [kg]
          IF ( id_HNO3 > 0 ) THEN

             !---------------------
             ! COUPLED SIMULATION
             !---------------------

             ! HNO3 [mole/m3] is in GAS(2); convert & store in Spc [kg]
             Spc(id_HNO3)%Conc(I,J,L) = MAX( 63.0e-3_fp * VOL * GNO3, CONMIN )

             ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
             HNO3_DEN = Spc(id_HNO3)%Conc(I,J,L)

          ELSE

             !---------------------
             ! OFFLINE SIMULATION:
             !---------------------

             ! Convert total inorganic nitrate from [mole/m3] to [ug/m3]
             ! and save for next time
             ! WT(4) is in [mole/m3] -- unit conv is necessary!
             CALL SET_HNO3( I, J, L, 63.0e+6_f8 * TNO3 )

             ! Save for use in sulfate_mod (SEASALT_CHEM) for offline
             ! aerosol simulations (bec, 4/15/05)
             GAS_HNO3(I,J,L) = GNO3

             ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
             HNO3_DEN        = GNO3 * VOL * 63.0e-3_fp

          ENDIF

          !-------------------------
          ! DIAGNOSTICS
          !-------------------------

          ! AEROPH: get pH related info to SAV arrays (hotp 8/11/09)
          ! HPLUSTEMP is H+ in mol/L water, AERLIQ1 is H, AERLIQ8 is H2O
          ! in mol/m3 air --> convert to mol/L water
          IF ( AERLIQ(8) < 1e-18_fp ) THEN
             ! Aerosol is dry so HPLUSTEMP and PH_SAV are undefined
             ! We force HPLUSTEMP to 1d20 (hotp, ccc, 12/18/09)
             ! Force IsorropAeropH to 20e0 (X. Wang, 6/27/19)
             !HPLUSTEMP       = 1e+20_fp
             HPLUSTEMP       = 1.0e-30_fp
             SULFTEMP        = 1.0e-30_fp
             BISULTEMP       = 1.0e-30_fp
             NITRTEMP        = 1.0e-30_fp
             CLTEMP          = 1.0e-30_fp
             !State_Chm%IsorropAeropH(I,J,L,N) = -999e+0_fp
             State_Chm%IsorropAeropH(I,J,L,N) = 20.0_fp
          ELSE
             HPLUSTEMP    = AERLIQ(1) / AERLIQ(8) * 1.0e+3_fp / 18.0_fp
             SULFTEMP     = AERLIQ(5) / AERLIQ(8) * 1.0e+3_fp / 18.0_fp
             BISULTEMP    = AERLIQ(6) / AERLIQ(8) * 1.0e+3_fp / 18.0_fp
             NITRTEMP     = AERLIQ(7) / AERLIQ(8) * 1.0e+3_fp / 18.0_fp
             CLTEMP       = AERLIQ(4) / AERLIQ(8) * 1.0e+3_fp / 18.0_fp

             ! Use SAFELOG10 to prevent NAN
             State_Chm%IsorropAeropH(I,J,L,N)=-1.0_fp*SAFELOG10(HPLUSTEMP)
          ENDIF

          ! Additional Info
          State_Chm%IsorropHplus(I,J,L,N)   = MAX(HPLUSTEMP, 1e-30_fp)
          State_Chm%IsorropAeroH2O(I,J,L,N) = MAX((AERLIQ(8)*18e+6_fp),1e-30_fp) ! mol/m3 -> ug/m3
          State_Chm%IsorropNitrate(I,J,L,N) = MAX(NITRTEMP, 1e-30_fp)
          State_Chm%IsorropChloride(I,J,L,N)= MAX(CLTEMP, 1e-30_fp)
          IF (N==1) THEN
             State_Chm%IsorropSulfate(I,J,L)  = MAX(SULFTEMP, 1e-30_fp)
             State_Chm%IsorropBisulfate(I,J,L)= MAX(BISULTEMP, 1e-30_fp)
             State_Chm%AeroH2O(I,J,L,1+NDUST) = AERLIQ(8) * 18e+0_fp ! mol/m3 -> g/m3

             NUM_SAV = ( Spc(id_NH3 )%Conc(I,J,L)  / 17.0_fp                 &
                     +   Spc(id_NH4 )%Conc(I,J,L)  / 18.0_fp                 &
                     +   Spc(id_SALA)%Conc(I,J,L) * 0.3061_fp / 23.0_fp     )

             ! HMS is only defined for fullchem is simulations,
             ! so skip it if it is not a defined species
             IF ( IS_HMS ) THEN
                DEN_SAV = ( Spc(id_SO4)%Conc(I,J,L)  / 96.0_fp   * 2.0_fp    &
                        +   Spc(id_HMS)%Conc(I,J,L)  / 111.0_fp              &
                        +   Spc(id_NIT)%Conc(I,J,L)  / 62.0_fp               &
                        +   HNO3_DEN           / 63.0_fp                     &
                        +   Spc(id_SALA)%Conc(I,J,L) * 0.55_fp   / 35.45_fp )
             ELSE
                DEN_SAV = ( Spc(id_SO4)%Conc(I,J,L)  / 96.0_fp   * 2.0_fp    &
                        +   Spc(id_NIT)%Conc(I,J,L)  / 62.0_fp               &
                        +   HNO3_DEN           / 63.0_fp                     &
                        +   Spc(id_SALA)%Conc(I,J,L) * 0.55_fp   / 35.45_fp )
             ENDIF

          ENDIF
       ENDDO

       IF ( id_HCl > 0 ) THEN
          PHCl = ( Spc(id_HCL)%Conc(I,J,L) - PHCl ) * 35.45_fp / 36.45_fp
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

    !### Debug
    IF ( prtDebug ) CALL DEBUG_MSG( '### ISORROPIAII: a DO_ISORROPIAII' )

  END SUBROUTINE DO_ISORROPIAII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safelog10
!
! !DESCRIPTION: Calculates the LOG (base 10) of a number X.  Returns a minimum
!  value if X is too small, in order to avoid NaN or Infinity problems.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SAFELOG10( X ) RESULT ( SAFLOG )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X        ! Argument for LOG10 function
!
! !RETURN VALUE:
!
    REAL(fp)             :: SAFLOG   ! LOG10 output --
!
! !REVISION HISTORY:
!  11 Aug 2009 - H. O. T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( X <= 1e-20_fp ) THEN
       SAFLOG = -1e+0_fp*20e+0_fp   ! if X<0, make pH 20
    ELSE
       SAFLOG = LOG10(X)
    ENDIF

  END FUNCTION SAFELOG10
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_hno3
!
! !DESCRIPTION: Subroutine SET\_HNO3 stores the modified HNO3 value back
!  into the HNO3\_sav array for the next timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_HNO3( I, J, L, HNO3_UGM3 )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I           ! GEOS-Chem longitude index
    INTEGER,  INTENT(IN) :: J           ! GEOS-Chem longitude index
    INTEGER,  INTENT(IN) :: L           ! GEOS-Chem longitude index
    REAL(f8), INTENT(IN) :: HNO3_UGM3   ! HNO3 concentration [ug/m3]
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
#if !defined( USE_REAL_8 )
    HNO3_sav(I,J,L) = SNGL(HNO3_UGM3) ! if we are not using real*8
#else
    HNO3_sav(I,J,L) = HNO3_UGM3       ! if we are using real*8
#endif

  END SUBROUTINE SET_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gno3
!
! !DESCRIPTION: Function GET\_GNO3 returns the gas-phase HNO3 [v/v] for
!  calculation of sea-salt chemistry in sulfate\_mod (SEASALT\_CHEM).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_GNO3( I, J, L, HNO3_kg, State_Met )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I          ! GEOS-Chem longitude index
    INTEGER,        INTENT(IN)    :: J          ! GEOS-Chem latitude index
    INTEGER,        INTENT(IN)    :: L          ! GEOS-Chem level index
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),          INTENT(OUT)  :: HNO3_kg    ! Gas-phase HNO3 [kg]
!
! !REVISION HISTORY:
!  15 Apr 2005 - B. Alexander - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Zero variables
    HNO3_kg  = 0.e+0_fp

    ! convert from [mole/m3] to [kg]
    HNO3_kg = GAS_HNO3(I,J,L) * 63.e-3_fp * State_Met%AIRVOL(I,J,L)

  END SUBROUTINE GET_GNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_Qk
!
! !DESCRIPTION: Function GET\_Qk returns the mass transfer correction
!  factor Qk for the thermodynamic equilibrium of coarse SSA. A dynamic
!  method of Pillinis et al. (2000) is used here.
!\\
!\\
! !INTERFACE:
!
      Subroutine GET_QK(GNO3, GCL, GNO3eq, GCLeq, Hplus, d, temp, &
           denair, n_x, Qk, F_HNO3, F_HCl )

      USE PhysConstants,      ONLY : pi
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
      ! Gas concentration (mole m-3) before equilibrium
      REAL(fp),  INTENT(IN) :: GNO3, GCL
      ! Gas concentration (mole m-3) after ISORROPIA equilibrium
      REAL(fp),  INTENT(IN) :: GNO3eq, GCLeq
      ! Aerosol H+ concentration (mole m-3)
      REAL(fp),  INTENT(IN) :: Hplus
      ! Aerosol (coarse seas salt) diameter (m)
      REAL(fp),  INTENT(IN) :: d
      ! Temprature (K)
      REAL(fp),  INTENT(IN) :: temp
      ! Air density (molec cm-3)
      REAL(fp),  INTENT(IN) :: denair
      ! Aerosol (coarse sea salt) number density (m-3)
      REAL(fp),  INTENT(IN) :: n_x
!
! !Return value:
!
      REAL(fp),  INTENT(OUT)  :: Qk, F_HNO3, F_HCl
!
! !REVISION HISTORY:
!  14 Feb 2018 - X. Wang - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
!  !LOCAL VARIABLES:
      REAL(fp) :: aa1, aa2, bb1, bb2, bb3
      REAL(fp) :: D_g, f, Kn, H_flux_NO3, H_flux_Cl
      REAL(fp) :: alpha, Hlim, ab


      ! Initialize
      aa1        = 0.0_fp
      aa2        = 0.0_fp
      bb1        = 0.0_fp
      bb2        = 0.0_fp
      bb3        = 0.0_fp
      D_g        = 0.0_fp
      f          = 0.0_fp
      Kn         = 0.0_fp
      H_flux_NO3 = 0.0_fp
      H_flux_Cl  = 0.0_fp
      Hlim       = 0.0_fp
      alpha      = 0.0_fp
      ab         = 0.0_fp

      ! Based on (Pilinis et al.AST,2000), for  HNO3+HCl case,
      ! Qx = -bb / aa, first calculate bb:
      ! bb = -D_g*f*C(HNO3) - D_g*f*C(HCl) + alpha/(2*pi*d*N)
      ! aa = D_g*f*Ceq(HNO3) + D_g*f*Ceq(HCl)

      ! In the 1st term for HNO3, D_g is the gas diffusion coefficient (m2 s-1)
      D_g = 9.45E+17_fp/denair * SQRT(temp) * 1.e-4_fp * &
            SQRT(3.472E-2_fp + 1.E+0_fp/63.e+0_fp)
      ! f is the the correction for noncontinuum effects and imperfect
      ! accommodation, based on Kn number and accommodation coefficient
      Kn = 0.068e-6_fp / (d / 2.)
      ! Use S&P eq12.42 (Dahneke 1983) to calculate f
      ab = 7.5e-5 * exp(2.1e3 / temp)
      f = (1+Kn) / (1 + 2.*Kn*(1+Kn)/ab)
      ! C is the gas concentration before equilibrium (mol m-3)
      bb1 = -D_g * f * GNO3 !(mole m-1 s-1)
      ! Also calulate availabe H+ flux for later calculation (mol m-3 s-1)
      H_flux_NO3 = 2*pi*D_g*d*n_x*f*(GNO3-GNO3eq)
      F_HNO3 = 2*pi*D_g*d*n_x*f * GET_TS_CHEM()! * 60e+0_fp !save for flux calc
      aa1 = D_g*f*GNO3eq !(mole m-1 s-1)

      ! Do similar for the 2nd term, for HCl
      D_g = 9.45E+17_fp/denair * SQRT(temp) * 1.e-4_fp * &
            SQRT(3.472E-2_fp + 1.E+0_fp/36.45e+0_fp)
      ab = 4.4e-6 * exp(2.898e3 / temp)
      f = (1+Kn) / (1 + 2.*Kn*(1+Kn)/ab)
      bb2 = -D_g * f * GCL !(mole m-1 s-1)
      H_flux_Cl = 2*pi*D_g*d*n_x*f*(GCL-GCLeq)
      F_HCl = 2*pi*D_g*d*n_x*f * GET_TS_CHEM()! * 60e+0_fp !save for flux calc
      aa2 = D_g*f*GCLeq !(mole m-1 s-1)

      ! In the third term, alpha is the limited changes in the acidity
      ! of the particle. The limitation is set to 0.1s-1 following Pilinis, 2000
      alpha = H_flux_Cl+H_flux_NO3 !(mole m-3 s-1)
      Hlim = 1.e-1_fp * Hplus !(mole m-3 s-1)
      IF (abs(alpha) .GT. Hlim) THEN
         alpha = sign(Hlim, alpha) !(mole m-3 s-1)
         bb3 = alpha/(2*pi*d*n_x) !(mole m-1 s-1)
         IF ((aa1+aa2) .LE. 0) THEN
            Qk = 1e+0_fp
         ELSEIF ((bb1 + bb2 + bb3) .GE. 0 ) THEN
            Qk = 1e+0_fp
         ELSE
            Qk = - (bb1 + bb2 + bb3) / (aa1 + aa2)
         ENDIF
      ELSE
         Qk = 1e+0_fp
      ENDIF

      END SUBROUTINE GET_Qk
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_isorropiaII
!
! !DESCRIPTION: Subroutine INIT\_ISORROPIAII initializes all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_ISORROPIAII( State_Grid )
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    !=================================================================
    ! INIT_ISORROPIAII begins here!
    !=================================================================

    ALLOCATE( HNO3_sav( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3_sav' )
    HNO3_sav = 0e+0_fp

    ALLOCATE( GAS_HNO3( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAS_HNO3' )
    GAS_HNO3 = 0e+0_fp

  END SUBROUTINE INIT_ISORROPIAII
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_isorropiaII
!
! !DESCRIPTION: Subroutine CLEANUP\_ISORROPIAII deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_ISORROPIAII
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( HNO3_sav    ) ) DEALLOCATE( HNO3_sav )
    IF ( ALLOCATED( GAS_HNO3    ) ) DEALLOCATE( GAS_HNO3 )

  END SUBROUTINE CLEANUP_ISORROPIAII
!EOC
END MODULE ISORROPIAII_MOD
