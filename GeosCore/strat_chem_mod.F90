!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: strat_chem_mod.F90
!
! !DESCRIPTION: Module STRAT\_CHEM\_MOD contains variables and routines for 
!  performing a simple linearized chemistry scheme.
!
!  If using a tropchem-based simulation, the simple linearized chemistry here 
!  will be applied to the stratosphere using archived 3D monthly production
!  rates and loss frequencies saved from a GEOS-Chem UCX simulation.
!
!  If using a UCX-based simulation, the simple linearized chemistry here will
!  be applied to the mesosphere instead using using archived 3D monthly
!  climatological production rates and loss frequencies from the GMI combo
!  model.
!
!  In the original schem code (schem.F), only the following species
!  were destroyed by photolysis in the stratosphere:
!    PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O, N2O5, HNO4, MP
!  and by reaction with OH:
!    ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
!    PRPE, C3H8, CH2O, C2H6, HNO4, MP
!  
!  The updated code includes at least all of these, and many more. The code
!  is flexible enough to automatically apply the rate to any new species
!  for future simulations that share the name in species\_mod with the
!  UCX or GMI name.  (See Documentation on wiki).
!
!  The prod rates and loss frequencies are now read via HEMCO. They are 
!  stored in a data structure of flexible length (PLVEC). The file containing
!  the prod rates and loss frequencies need to be specified in the HEMCO 
!  configuration file for each species of interest. They are then automatically 
!  read and remapped onto the simulation grid.
! 
!  The field names assigned to the production and loss fields are expected to 
!  be 'UCX\_PROD\_XXX' and 'UCX\_LOSS\_XXX' (or 'GMI\_PROD\_XXX' and
!  'GMI\_LOSS\_XXX'), respectively, where XXX is the species name. Production
!  rates must be given in units of v/v/s, and loss frequencies in s-1.
!
!  The module variable PLMUSTFIND (set below) determines the behavior if no
!  production rates and/or loss frequencies can be found for any of the UCX/GMI
!  species defined in this module. IF PLMUSTFIND is set to TRUE, the code stops
!  with an error if no entry is found. Otherwise, steady-state values are used
!  for all species with no explicitly given values. 
!
!  The (monthly) OH concentrations are also obtained through HEMCO. The field
!  name must be 'STRAT\_OH', and values must be in v/v.
!\\
!\\
! !INTERFACE:
!
MODULE Strat_Chem_Mod
!
! !USES:
!
! for precisions
  USE HCO_Error_Mod
  USE Precision_Mod    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Strat_Chem
  PUBLIC  :: Do_Strat_Chem
  PUBLIC  :: Cleanup_Strat_Chem
!----------------------------------------------------------------------------
! Prior to 8/9/18:
! Disable CALC_STE routine, which can give misleading STE info (bmy, 8/9/18)
!  PUBLIC  :: Calc_STE
!----------------------------------------------------------------------------
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Set_BryPointers
  PRIVATE :: Set_PLVEC
  PRIVATE :: Do_Synoz
!
! !PUBLIC DATA MEMBERS:
!
  PUBLIC  :: SChem_Tend
  PUBLIC  :: Minit_Is_Set
!
! !REMARKS:
!
!  References:
!==============================================================================
!  (1 )
!
! !REVISION HISTORY:
!  01 Feb 2011 - L. Murray   - Initial version
!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!  20 Jul 2012 - R. Yantosca - Correct compilation error in GET_RATES_INTERP
!  07 Aug 2012 - R. Yantosca - Fix parallelization problem in Bry do loop
!  05 Oct 2012 - R. Yantosca - Add bug fix for IFORT 12 compiler in CALC_STE
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  20 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  30 Dec 2014 - C. Keller   - Now read Bry data through HEMCO
!  16 Jan 2015 - C. Keller   - Now read all prod/loss fields and OH conc.
!                              through HEMCO.
!   4 Mar 2015 - R. Yantosca - Declare pointer args for HCO_GetPtr as REAL(f4)
!  06 Apr 2016 - C. Keller   - Add Minit_Is_Set and SET_MINIT.
!  03 Oct 2016 - R. Yantosca - Now dynamically allocate BrPtrDay, BrPtrNight
!  03 Oct 2016 - R. Yantosca - Dynamically allocate BrPtrDay, BrPtrNight
!  30 May 2017 - M. Sulprizio- SChem_Tend is now public so that it can be
!                              updated in flexchem_mod.F90 for UCX simulations
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  06 Mar 2018 - M. Sulprizio- Update module to use prod/loss rates from UCX
!                              instead of GMI
!  23 Mar 2018 - M. Sulprizio- Restore using prod/loss rates from GMI for UCX-
!                              based simulations. Tropchem-based simulations
!                              will use the rates from UCX.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Species index of Bry species in input files
  ! 1:6 = (Br2, Br, BrO, HOBr, HBr, BrNO3) for both
  INTEGER,           PARAMETER :: br_nos(6)   = (/ 44, 45, 46, 47, 48, 50 /)

  ! BrPointers is a derived type to hold pointers to the Bry day 
  ! and night data 
  TYPE :: BrPointers
     REAL(f4), POINTER  :: MR(:,:,:) => NULL()
  END TYPE BrPointers

  ! Vectors holding the mixing ratios (MR) of day/night Bry data.
  ! The mixing ratios will be read and interpolated by HEMCO. The
  ! corresponding HEMCO fields must be specified in the HEMCO 
  ! configuration file. The field names are assumed to be 
  ! 'GEOSCCM_XY_DAY' and 'GEOSCCM_XY_NIGHT', respectively, where 
  ! XY is the Bry species name. It is also assumed that the input
  ! data is in ppt.
  !
  ! NOTE: Need to dynamically allocate these to avoid having to 
  ! declare these with the SAVE attribute (bmy, 10/3/16)
  TYPE(BrPointers), POINTER :: BrPtrDay(:) 
  TYPE(BrPointers), POINTER :: BrPtrNight(:) 

  ! PL_Pointers is a derived type to hold pointers to the production
  ! and loss fields 
  TYPE :: PL_Pointers 
     REAL(f4), POINTER  :: PROD(:,:,:) => NULL() ! Production rate [v/v/s]
     REAL(f4), POINTER  :: LOSS(:,:,:) => NULL() ! Loss frequency [s-1]
  END TYPE PL_Pointers 

  ! Monthly mean OH [v/v]
  REAL(f4), POINTER     :: STRAT_OH(:,:,:) => NULL()

  ! Vector holding the prod/loss arrays of all active strat chem 
  ! species.
  TYPE(PL_Pointers), POINTER :: PLVEC(:) => NULL()

  ! Toggle to specify whether or not production/loss rates must be provided
  ! for every strat chem species. If set to TRUE, the code will return with
  ! an error if the production/loss rate cannot be found (through HEMCO) for
  ! any of the species. If set to .FALSE., only a warning is prompted and a
  ! value of 0.0 is used for every field that cannotbe found. 
  LOGICAL, PARAMETER   :: PLMUSTFIND = .FALSE.

  ! Number of species from UCX or GMI
  INTEGER              :: NTR

  ! Minit_Is_Set indicates whether Minit has been set or not
  LOGICAL              :: Minit_Is_Set = .FALSE. 
!
! !PRIVATE TYPES:
!
  ! Scalars
  REAL(fp)              :: dTchem          ! chemistry time step [s]
  INTEGER               :: NSCHEM          ! Number of species upon which to 
                                           ! apply P's & k's in GEOS-Chem
  ! Arrays
  CHARACTER(LEN=16), ALLOCATABLE :: TrName(:) !Tracer names in UCX/GMI
  INTEGER, ALLOCATABLE  :: Strat_TrID_GC (:)  !Maps 1:NSCHEM to SC%Species
  INTEGER, ALLOCATABLE  :: Strat_TrID_TND(:)  !Maps 1:NSCHEM to SCHEM_TEND
  INTEGER, ALLOCATABLE  :: Strat_TrID(:)      !Maps 1:NSCHEM to UCX/GMI index

  ! Species index of Bry species in GEOS-Chem species DB(may differ from br_nos)
  INTEGER               :: GC_Bry_TrID(6) 

  ! Variables used to calculate the strat-trop exchange flux
  REAL(fp)              :: TauInit             ! Initial time
  INTEGER               :: NymdInit, NhmsInit  ! Initial date
  REAL(fp)              :: TpauseL_Cnt         ! Tropopause counter
  REAL(fp), ALLOCATABLE :: TpauseL(:,:)        ! Tropopause level aggregator
  REAL(f4), ALLOCATABLE :: MInit(:,:,:,:)      ! Init. atm. state for STE period
  REAL(f4), ALLOCATABLE :: SChem_Tend(:,:,:,:) ! Stratospheric chemical tendency
                                               !   (total P - L) [kg period-1]
  ! Species ID flags
  INTEGER               :: id_Br2,   id_Br,     id_BrNO3
  INTEGER               :: id_BrO,   id_CHBr3,  id_CH2Br2
  INTEGER               :: id_CH3Br, id_HOBr,   id_HBr
  INTEGER               :: id_O3,    id_O3Strat 

  !=================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
  !=================================================================
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Strat_Chem
!
! !DESCRIPTION: Function DO\_STRAT\_CHEM is the driver routine for computing
! the simple linearized stratospheric chemistry scheme.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_STRAT_CHEM( am_I_Root, Input_Opt,          &
                            State_Met, State_Chm, errCode )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE LINOZ_MOD,          ONLY : DO_LINOZ
    USE PhysConstants,      ONLY : XNUMOLAIR, AIRMW, AVO
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : TIMESTAMP_STRING
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object 
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: errCode     ! Success or failure
! 
! !REMARKS:
! 
! !REVISION HISTORY: 
!  01 Feb 2011 - L. Murray   - Initial version  
!  18 Jul 2012 - R. Yantosca - For compatibility w/ the GEOS-5/GCM, we cannot
!                              assume a minimum tropopause level anymore
!  18 Jul 2012 - R. Yantosca - Make sure I is the innermost DO loop
!                              wherever expedient 
!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  07 Aug 2012 - R. Yantosca - Make BEFORE a local variable for parallel loop
!  26 Oct 2012 - R. Yantosca - Now pass the Chemistry State object for GIGC
!  09 Nov 2012 - R. Yantosca - Now pass the Input Options object for GIGC
!  15 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  18 Mar 2013 - R. Yantosca - Now pass Input_Opt via the arg list
!  19 Mar 2013 - R. Yantosca - Now only copy Input_Opt%TCVV(1:N_TRACERS)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  30 Dec 2014 - C. Keller   - Now get Bry data through HEMCO.
!  24 Mar 2015 - E. Lundgren - Replace dependency on tracer_mod with
!                              CMN_GTCM_MOD for XNUMOLAIR
!  30 Sep 2015 - E. Lundgren - Now use UNITCONV_MOD for unit conversion
!  05 Mar 2016 - C. Keller   - Allow O3 P/L be done by GMI if both LINOZ and
!                              SYNOZ are disabled. This is primarily for 
!                              testing/data assimilation applications.
!  16 Jun 2016 - M. Yannetti - Replaced TRACERID_MOD.\
!  20 Jun 2016 - R. Yantosca - Now make species ID flags module variables
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  12 Jul 2016 - R. Yantosca - Bug fix: ISBR2 should be held !$OMP PRIVATE
!  18 Jul 2016 - M. Yannetti - Replaced TCVV with spec db and phys constant
!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  19 Oct 2016 - R. Yantosca - Add routine Set_Init_Conc_Strat_Chem for GCHP
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  17 Jan 2018 - R. Yantosca - Replace GET_TPAUSE_LEVEL w/ State_Met%TropLev
!  09 Apr 2018 - R. Yantosca - Now use safe division in expression for k
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd quantities
    LOGICAL, SAVE      :: FIRST      = .TRUE.
    INTEGER, SAVE      :: LASTMONTH  = -99

    ! Flags for simulation types
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_A_TAGO3_SIM
    LOGICAL            :: IT_IS_A_H2HD_SIM

    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: I,     J,       L,   N
    INTEGER            :: NN,    nAdvect, NA
    REAL(fp)           :: dt,    P,       k,   M0,  RC
    REAL(fp)           :: TK,    RDLOSS,  T1L, mOH, BryTmp
    REAL(fp)           :: BOXVL, SpcConc, Num, Den, M
    REAL(fp)           :: MW_g
    LOGICAL            :: LLINOZ
    LOGICAL            :: LSYNOZ
    LOGICAL            :: LPRT
    LOGICAL            :: LUCX
    LOGICAL            :: LCYCLE
    LOGICAL            :: ISBR2
#if defined( MODEL_GEOS ) 
    LOGICAL            :: SKIP  
#endif

    ! Strings
    CHARACTER(LEN=16)  :: STAMP
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    ! Arrays
    REAL(fp)           :: Spc0  (IIPAR,JJPAR,LLPAR,State_Chm%nAdvect)
    REAL(fp)           :: BEFORE(IIPAR,JJPAR,LLPAR                  )

    ! Pointers
    REAL(fp), POINTER  :: Spc(:,:,:,:)
    REAL(fp), POINTER  :: AD (:,:,:  )
    REAL(fp), POINTER  :: T  (:,:,:  )

    !=======================================================================
    ! DO_STRAT_CHEM begins here!
    !=======================================================================

    ! Initialize
    errCode = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Strat_Chem (in module GeosCore/strat_chem_mod.F90)'

    ! Initialize further
    LLINOZ               = Input_Opt%LLINOZ
    LSYNOZ               = Input_Opt%LSYNOZ
    LPRT                 = Input_Opt%LPRT
    LUCX                 = Input_Opt%LUCX
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGO3_SIM    = Input_Opt%ITS_A_TAGO3_SIM  
    IT_IS_A_H2HD_SIM     = Input_Opt%ITS_A_H2HD_SIM
    Spc                  => NULL()
    AD                   => NULL()
    T                    => NULL()

#if defined( MODEL_GEOS )
    ! Skip strat chem if chemistry is over entire vertical domain
    SKIP = .FALSE.
    IF ( LLCHEM == LLPAR ) THEN
       SKIP = .TRUE.
       IF ( FIRST .AND. am_I_Root ) THEN
          WRITE( 6, * ) 'Fullchem up to top of atm - skip linearized strat chemistry'
       ENDIF 
    ENDIF
    IF ( .NOT. SKIP ) THEN
#endif

    ! Set a flag for debug printing
    prtDebug             = ( LPRT .and. am_I_Root )

    STAMP = TIMESTAMP_STRING()
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 ) STAMP
    ENDIF
10  FORMAT( '     - DO_STRAT_CHEM: Linearized strat chemistry at ', a )
    
    !======================-================================================
    ! On first call, establish pointers to data fields read by HEMCO. These 
    ! are the stratospheric Bry fields as well as the production/loss rates.
    ! (ckeller, 12/30/2014)
    !
    ! If we are doing a tagO3 simulation, then we can skip this section,
    ! since tagO3 only uses Linoz or Synoz, but doesn't read in any P/L
    ! fields from disk. (bmy, 7/11/16)
    !======================-================================================
    IF ( FIRST ) THEN
       IF ( .not. IT_IS_A_TAGO3_SIM ) THEN

          ! BrPtrDay and BrPtrNight have to be allocated dynamically
          ! because they are pointers (bmy, 10/3/16)
          ALLOCATE( BrPtrDay  ( 6 ), STAT=errCode )
          ALLOCATE( BrPtrNight( 6 ), STAT=errCode )

          ! Get pointers to Bry fields via HEMCO
          CALL Set_BryPointers ( am_I_Root, Input_Opt,          &
                                 State_Chm, State_Met, errCode )

          ! Trap potential errors
          IF ( errCode /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Set_BryPointers"!'
             CALL GC_Error( ErrMsg, errCode, ThisLoc )
             RETURN
          ENDIF

          ! Get pointers to prod/loss fields via HEMCO
          CALL Set_PLVEC ( am_I_Root, Input_Opt,                &
                           State_Chm, State_Met, errCode )

          ! Trap potential errors
          IF ( errCode /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Set_PlVec"!'
             CALL GC_Error( ErrMsg, errCode, ThisLoc )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### STRAT_CHEM: at DO_STRAT_CHEM' )
    ENDIF

#if defined( ESMF_ )
    ! Eventually set Minit if in ESMF environment (ckeller, 4/6/16)
    IF ( .NOT. Minit_Is_Set ) THEN
       CALL SET_MINIT( am_I_Root, Input_Opt, State_Met, &
                       State_Chm, errCode )

       ! Trap potential errors
       IF ( errCode /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Set_Minit"!'
          CALL GC_Error( ErrMsg, errCode, ThisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    !======================-================================================
    ! FULL CHEMISTRY SIMULATIONS
    !
    ! %%% NOTE: For now, the algorithm assumes that the advected species
    ! %%% are listed first.  We may have to also store the advected ID's
    ! %%% in an array for SCHEM_TEND for tropchem simulations. SCHEM_TEND
    ! %%% for UCX simulations is handled in flexchem_mod.F90.
    !=======================================================================
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       ! Initialize pointers
       Spc         => State_Chm%Species
       AD          => State_Met%AD
       T           => State_Met%T

       ! Advance counter for number of times we've sampled the tropopause level
       TpauseL_CNT = TpauseL_CNT + 1e+0_fp

       ! Make note of inital state for determining tendency later
       BEFORE = Spc(:,:,:,id_O3)

       !--------------------------------------------------------------------
       ! Do chemical production and loss for non-ozone species for
       ! which we have explicit prod/loss rates from UCX/GMI
       !--------------------------------------------------------------------

       ! Timestep [s], can be pulled outside of parallel loop
       dt = DTCHEM

       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED )                                               &
       !$OMP PRIVATE( I, J, L, N, NN, NA, k, P, M0, MW_g, SpcConc, Num, Den )
       DO J=1,JJPAR
          DO I=1,IIPAR

             ! Add to tropopause level aggregator for later determining STE flux
             TpauseL(I,J) = TpauseL(I,J) + State_Met%TropLev(I,J)

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR

                ! For safety's sake, zero variables at each (I,J,L) box
                Den     = 0.0_fp
                k       = 0.0_fp
                M0      = 0.0_fp
                Mw_g    = 0.0_fp
                Num     = 0.0_fp
                P       = 0.0_fp
                SpcConc = 0.0_fp

                ! Only consider boxes above the chemistry grid
                IF ( State_Met%InChemGrid(I,J,L) ) CYCLE

                !-----------------------------------------------------------
                ! Loop over the # of active strat chem species
                !-----------------------------------------------------------
                DO N = 1, NSCHEM

                   ! Species ID (use this for State_Chm%Species)
                   NN = Strat_TrID_GC(N)

                   ! Advected species ID (use this for SCHEM_TEND)
                   NA = Strat_TrID_TND(N)

                   ! Skip O3; we'll always use either Linoz or Synoz
                   ! Use UCX/GMI O3 P/L if both LINOZ and SYNOZ are deactivated.
                   ! This is for test use in assimilation mode (ckeller, 3/7/16)
                   IF ( IT_IS_A_FULLCHEM_SIM .and. NN .eq. id_O3 .AND. &
                        ( LLINOZ .OR. LSYNOZ ) ) CYCLE

                   ! Molecular weight for the species [g]
                   MW_g = State_Chm%SpcData(NN)%Info%emMW_g

                   !--------------------------------------------------------
                   ! Loss freq [s-1] 
                   !--------------------------------------------------------
                   IF ( .NOT. ASSOCIATED(PLVEC(N)%LOSS) ) THEN

                      ! If the PLVEC(N)%LOSS pointer is null, then the
                      ! species doesn't have a loss rate.  Set k to 0.
                      k = 0.0_fp

                   ELSE

                      ! Determine if we are using a UCX-based simulation
                      ! or a tropchem simulation (without UCX)
                      IF ( LUCX ) THEN

                         ! %%%%% SIMULATION USES THE UCX MECHANISM %%%%%
                         ! Loss rates (archived from GMI) are in s-1
                         k = PLVEC(N)%LOSS(I,J,L)

                      ELSE

                         ! %%%%% SIMULATION USES THE TROPCHEM MECHANISM %%%%%
                         ! Loss rates (archived from UCX) are in molec/cm3/s
                         ! Convert to s-1 here
                         SpcConc = Spc(I,J,L,NN)                             &
                                   * ( AVO / ( MW_g * 1.e-3_fp ) )           &
                                   / ( State_Met%AIRVOL(I,J,L) * 1e+6_fp )

                         ! k is the loss rate divided by the concentration
                         ! Return 0 if the division cannot be done
                         Num = PLVEC(N)%LOSS(I,J,L)
                         Den = SpcConc
                         k   = Safe_Div( N         = Num,                    &
                                         D         = Den,                    &
                                         Alt_NaN   = 0.0_fp,                 &
                                         Alt_Over  = 0.0_fp,                 &
                                         Alt_Under = 0.0_fp                 )

                      ENDIF

                   ENDIF

                   !--------------------------------------------------------
                   ! Prod term [kg/s]
                   !--------------------------------------------------------
                   IF ( .NOT. ASSOCIATED(PLVEC(N)%PROD) ) THEN

                      ! If the PLVEC(N)%PROD pointer is null, then the
                      ! species doesn't have a production rate.  Set P to 0.
                      P = 0.0_fp 

                   ELSE

                      ! Determine if we are using a UCX-based simulation
                      ! or a tropchem simulation (without UCX) 
                      IF ( LUCX ) THEN

                         ! %%%%% SIMULATION USES THE UCX MECHANISM %%%%%
                         ! Prod rates (archived from GMI) are in v/v/s
                         ! Convert to kg/s here
                         P = PLVEC(N)%PROD(I,J,L)                            &
                           * AD(I,J,L) / ( AIRMW / MW_g )

                      ELSE

                         ! %%%%% SIMULATION USES THE TROPCHEM MECHANISM %%%%%
                         ! Prod rates (archived from UCX) are in molec/cm3/s
                         ! Convert to kg/s here
                         P = PLVEC(N)%PROD(I,J,L)                            &
                           / ( AVO / ( MW_g  * 1.e-3_fp ) )                  &
                           * ( State_Met%AIRVOL(I,J,L) * 1e+6_fp ) 

                      ENDIF

                   ENDIF

                   !--------------------------------------------------------
                   ! Apply prod and loss
                   !--------------------------------------------------------

                   ! Initial mass [kg]
                   M0 = Spc(I,J,L,NN)

                   ! No prod or loss at all
                   ! NOTE: Bad form to test for equality on zero!
                   ! Replace this later (bmy, 4/9/18)
                   IF ( k .eq. 0e+0_fp .and. P .eq. 0e+0_fp ) CYCLE

                   ! Simple analytic solution to dM/dt = P - kM over [0,t]
                   IF ( k .gt. 0e+0_fp ) then
                      Spc(I,J,L,NN) = M0 * EXP(-k*dt) + &
                                      (P/k)*(1e+0_fp-EXP(-k*dt))
                   ELSE
                      Spc(I,J,L,NN) = M0 + P*dt
                   ENDIF

                   IF ( .not. LUCX ) THEN
                      ! Aggregate stratospheric chemical tendency [kg box-1]
                      ! for tropchem simulations
                      SCHEM_TEND(I,J,L,NA) = SCHEM_TEND(I,J,L,NA) + &
                                             ( Spc(I,J,L,NN) - M0 )
                   ENDIF

                ENDDO ! N
             ENDDO    ! L
          ENDDO       ! I
       ENDDO          ! J
       !$OMP END PARALLEL DO

       !--------------------------------------------------------------------
       ! Ozone
       !--------------------------------------------------------------------

       IF ( LLINOZ .OR. LSYNOZ ) THEN

          ! Put ozone in [v/v] for Linoz or Synoz
          Spc(:,:,:,id_O3) = Spc(:,:,:,id_O3) * ( AIRMW  &
                             / State_Chm%SpcData(id_O3)%Info%emMW_g ) / AD

          ! Do Linoz or Synoz
          IF ( LLINOZ ) THEN
             CALL Do_Linoz( am_I_Root, Input_Opt,             &
                            State_Met, State_Chm, RC=errCode )
          ELSE
             CALL Do_Synoz( am_I_Root, Input_Opt,             &
                            State_Met, State_Chm, RC=errCode )
          ENDIF

          ! Put ozone back to [kg]
          Spc(:,:,:,id_O3) = Spc(:,:,:,id_O3) * AD / ( AIRMW  &
                             / State_Chm%SpcData(id_O3)%Info%emMW_g )

       ENDIF
 
       IF ( .not. LUCX ) THEN
          ! Aggregate stratospheric chemical tendency [kg box-1]
          ! for tropchem simulations
          SCHEM_TEND(:,:,:,id_O3) = SCHEM_TEND(:,:,:,id_O3) + &
                                    ( Spc(:,:,:,id_O3) - BEFORE )
       ENDIF

       !--------------------------------------------------------------------
       ! Reactions with OH
       ! Currently:
       !   (1) CHBr3  
       !   (2) CH2Br2  
       !   (3) CH3Br
       !--------------------------------------------------------------------

       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, M, TK, RC, RDLOSS, T1L, mOH, BOXVL )
       DO J=1,JJPAR
          DO I=1,IIPAR  

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR

                IF ( State_Met%InChemGrid(I,J,L) ) CYCLE

                ! Grid box volume [cm3]
                BOXVL = State_Met%AIRVOL(I,J,L) * 1e+6_fp

                ! Density of air at grid box (I,J,L) in [molec cm-3]
                M = AD(I,J,L) / BOXVL * XNUMOLAIR

                ! OH number density [molec cm-3]
                mOH = M * STRAT_OH(I,J,L)

                ! Temperature at grid box (I,J,L) in K
                TK = T(I,J,L)

                !============!
                ! CH3Br + OH !
                !============!
                IF ( id_CH3Br .gt. 0 ) THEN
                   RC = 2.35e-12_fp * EXP ( - 1300.e+0_fp / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = Spc(I,J,L,id_CH3Br) * RDLOSS
                   Spc(I,J,L,id_CH3Br) = Spc(I,J,L,id_CH3Br) - T1L
                   IF ( .not. LUCX ) THEN
                      ! Aggregate stratospheric chemical tendency [kg box-1]
                      ! for tropchem simulations
                      SCHEM_TEND(I,J,L,id_CH3Br) = &
                        SCHEM_TEND(I,J,L,id_CH3Br) - T1L
                   ENDIF
                ENDIF

                !============!
                ! CHBr3 + OH !
                !============!
                IF ( id_CHBr3 .gt. 0 ) THEN
                   RC = 1.35e-12_fp * EXP ( - 600.e+0_fp / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = Spc(I,J,L,id_CHBr3) * RDLOSS
                   Spc(I,J,L,id_CHBr3) = Spc(I,J,L,id_CHBr3) - T1L
                   IF ( .not. LUCX ) THEN
                      ! Aggregate stratospheric chemical tendency [kg box-1]
                      ! for tropchem simulations
                      SCHEM_TEND(I,J,L,id_CHBr3) = &
                        SCHEM_TEND(I,J,L,id_CHBr3) - T1L
                   ENDIF
                ENDIF

                !=============!
                ! CH2Br2 + OH !
                !=============!
                IF ( id_CH2Br2 .gt. 0 ) THEN
                   RC = 2.00e-12_fp * EXP ( -  840.e+0_fp / TK )
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = Spc(I,J,L,id_CH2Br2) * RDLOSS
                   Spc(I,J,L,id_CH2Br2) = Spc(I,J,L,id_CH2Br2) - T1L
                   IF ( .not. LUCX ) THEN
                      ! Aggregate stratospheric chemical tendency [kg box-1]
                      ! for tropchem simulations
                      SCHEM_TEND(I,J,L,id_CH2Br2) = &
                        SCHEM_TEND(I,J,L,id_CH2Br2) - T1L
                   ENDIF
                ENDIF

             ENDDO ! J
          ENDDO ! I
       ENDDO ! L

       !$OMP END PARALLEL DO

       !--------------------------------------------------------------------
       ! Prescribe Br_y concentrations
       !--------------------------------------------------------------------

       !$OMP PARALLEL DO                                           &
       !$OMP DEFAULT( SHARED                                     ) &
       !$OMP PRIVATE( NN, BEFORE, ISBR2, L, J, I, LCYCLE, BryTmp )
       DO NN = 1,6

          IF ( GC_Bry_TrID(NN) > 0 ) THEN

             ! Make note of inital state for determining tendency later
             ! NOTE: BEFORE has to be made PRIVATE to the DO loop since
             ! it only has IJL scope, but the loop is over IJLN!
             ! (bmy, 8/7/12)
             BEFORE = Spc(:,:,:,GC_Bry_TrID(NN))

             ! Is this Br2?
             ISBR2  = ( Gc_Bry_TrId(NN) == id_Br2 )

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR
             DO J = 1, JJPAR
             DO I = 1, IIPAR  
                 
                LCYCLE = State_Met%InChemGrid(I,J,L)
                IF ( LCYCLE ) CYCLE

                ! Now get Br data through HEMCO pointers (ckeller, 12/30/14).
                IF ( State_Met%SUNCOS(I,J) > 0.e+0_fp ) THEN
                   ! daytime [ppt] -> [kg]
                   BryTmp = BrPtrDay(NN)%MR(I,J,L)   &
                          * 1.e-12_fp                & ! convert from [ppt]
                          * AD(I,J,L)                &
                          / ( AIRMW                  &
                          / State_Chm%SpcData(GC_Bry_TrID(NN))%Info%emMW_g )

                ELSE
                   ! nighttime [ppt] -> [kg]
                   BryTmp = BrPtrNight(NN)%MR(I,J,L) &
                          * 1.e-12_fp                & ! convert from [ppt]
                          * AD(I,J,L)                &
                          /  ( AIRMW                 & 
                          / State_Chm%SpcData(GC_Bry_TrID(NN))%Info%emMW_g )
                ENDIF

                ! Special adjustment for G-C Br2, 
                ! which is BrCl in the strat (ckeller, 1/2/15)
                IF ( ISBR2 ) BryTmp = BryTmp / 2.0_fp

                ! Pass to Spc array
                Spc(I,J,L, GC_Bry_TrID(NN) ) = BryTmp

             ENDDO
             ENDDO
             ENDDO

             IF ( .not. LUCX ) THEN
                ! Aggregate stratospheric chemical tendency [kg box-1]
                ! for tropchem simulations
                SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) = &
                   SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) + &
                   ( Spc(:,:,:,GC_Bry_TrID(NN)) - BEFORE )
             ENDIF
          
          ENDIF

       ENDDO ! NN
       !$OMP END PARALLEL DO

       ! Free pointers
       Spc => NULL()
       AD  => NULL()
       T   => NULL()

    !======================================================================
    ! TAGGED O3 SIMULATION
    !
    ! Tagged O3 only makes use of Synoz or Linoz. We apply either to
    ! the total O3 species, and the stratospheric O3 species.
    !======================================================================
    ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN

       ! Number of advected species
       nAdvect = State_Chm%nAdvect

       ! Loop over only the advected species
       DO NA = 1, nAdvect

          ! Get the species ID from the advected species ID
          N                          = State_Chm%Map_Advect(NA)

          ! Save initial conditions in Spc0 [kg]
          Spc0(:,:,:,NA)             = State_Chm%Species(:,:,:,N)

       ENDDO

       IF ( LLINOZ .OR. LSYNOZ ) THEN

          ! Convert units to [v/v dry air] for Linoz and Synoz (ewl, 10/05/15)
          CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                                  State_Chm, 'v/v dry', errCode,   &
                                  OrigUnit=OrigUnit )

          ! Trap potential errors
          IF ( errCode /= GC_SUCCESS ) THEN
             ErrMsg = 'Unit conversion error (forward) in "Convert_Spc_Units"!'
             CALL GC_Error( ErrMsg, errCode, ThisLoc )
             RETURN
          ENDIF

          ! Do LINOZ or SYNOZ
          IF ( LLINOZ ) THEN
             CALL Do_Linoz( am_I_Root, Input_Opt,             &
                            State_Met, State_Chm, errCode )
          ELSE 
             CALL Do_Synoz( am_I_Root, Input_Opt,             &
                            State_Met, State_Chm, errCode )
          ENDIF

         ! Convert species units back to original unit 
         CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                                 State_Chm, OrigUnit,  errCode )

          ! Trap potential errors
          IF ( errCode /= GC_SUCCESS ) THEN
             ErrMsg = 'Unit conversion error (backward) in "Convert_Spc_Units"!'
             CALL GC_Error( ErrMsg, errCode, ThisLoc )
             RETURN
          ENDIF

       ENDIF

       ! Add to tropopause level aggregator for later determining STE flux
       TpauseL_CNT = TpauseL_CNT + 1e+0_fp

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          TpauseL(I,J) = TpauseL(I,J) + State_Met%TropLev(I,J)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Loop over the # of strat chem species
       DO N = 1, NSCHEM

          ! Species ID (use for State_Chm%Species)
          NN = Strat_TrID_GC(N)

          ! Advected species ID (use for SCHEM_TEND)
          NA = Strat_TrId_TND(N)

          ! Aggregate chemical tendency [kg box-1]
          SCHEM_TEND(:,:,:,NA) = SCHEM_TEND(:,:,:,NA) + &
               ( State_Chm%Species(:,:,:,NN) - Spc0(:,:,:,NA) )
       ENDDO

    !======================================================================
    ! OTHER SIMULATIONS
    !
    ! The code will need to be modified for other tagged simulations 
    ! (e.g., CO). Simulations like CH4, CO2 with standard species names 
    ! should probably just work as is with the full chemistry code above, 
    ! but would need to be tested.
    !======================================================================
    ELSE

       ! Free pointers
       Spc => NULL()
       AD  => NULL()
       T   => NULL() 

       ! Exit with error if strat chemistry is not activated
       ErrMsg = 'Linearized strat chemistry needs to be activated for '   // &
                'activated for your simulation type.  Please see      '   // &
                'module see GeosCore/strat_chem_mod.F90 or disable    '   // &
                'stratospheric chemistry in your "input.geos" file.'
       CALL GC_Error( ErrMsg, errCode, ThisLoc )
       RETURN

    ENDIF
   
#if defined( MODEL_GEOS )
    ! End of SKIP loop
    ENDIF ! SKIP 
#endif

    ! Set first-time flag to false
    FIRST = .FALSE.    

    ! Free pointer
    Spc => NULL()
    AD  => NULL()
    T   => NULL()  

  END SUBROUTINE Do_Strat_Chem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_BryPointers 
!
! !DESCRIPTION: Subroutine SET\_BryPointers gets the Bry stratospheric data 
! read by HEMCO. The pointers only need to be established once. Target data
! is automatically updated through HEMCO. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_BryPointers( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr 
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorological State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=16)   :: ThisName
    CHARACTER(LEN=255)  :: PREFIX, FIELDNAME, ThisLoc
    CHARACTER(LEN=1024) :: ErrMsg
    INTEGER             :: N

    !=================================================================
    ! Set_BryPointers begins here 
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Set_BryPointers (in module GeosCore/strat_chem_mod.F90)'

    ! Do for every Bry species
    DO N = 1, 6

       ! Skip if species is not defined    
       IF ( GC_Bry_TrID(N) <= 0 ) CYCLE

       ! Get Bry name
       ThisName = State_Chm%SpcData(GC_Bry_TrID(N))%Info%Name

       ! Construct field name using Bry name
       PREFIX = 'GEOSCCM_'//TRIM(ThisName)

       ! Get pointer to this field. These are the mixing ratios (pptv).

       ! Day
       FIELDNAME = TRIM(PREFIX) // '_DAY'
       CALL HCO_GetPtr( am_I_Root, HcoState, FIELDNAME, BrPtrDay(N)%MR, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO! Stratospheric Bry data '//&
                   'is expected to be listed in the HEMCO configuration '  //&
                  'file. This error occured when trying to get field '     //&
                  TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Night
       FIELDNAME = TRIM(PREFIX) // '_NIGHT'
       CALL HCO_GetPtr( am_I_Root, HcoState, FIELDNAME, BrPtrNight(N)%MR, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO! Stratospheric Bry data '//&
                   'is expected to be listed in the HEMCO configuration '  //&
                   'file. This error occured when trying to get field '    //&
                   TRIM( FieldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDDO !N

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Set_BryPointers 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Plvec
!
! !DESCRIPTION: Subroutine SET\_PLVEC gets the production and loss terms of
! all strat chem species from HEMCO. The pointers only need to be established 
! once. Target data is automatically updated through HEMCO. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_PLVEC( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr 
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorological State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  16 Jan 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=16)   :: ThisName
    CHARACTER(LEN=255)  :: FIELDNAME
    INTEGER             :: N, TRCID
    LOGICAL             :: FND
    INTEGER             :: I, J, L

    ! Strings
    CHARACTER(LEN=63)   :: OrigUnit
    CHARACTER(LEN=255)  :: ThisLoc
    CHARACTER(LEN=1024) :: ErrMsg

    !=================================================================
    ! Set_PLVEC begins here 
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Set_PLVEC (in module GeosCore/strat_chem_mod.F90)'

    ! Do for every species 
    DO N = 1,NSCHEM

       ! Get GEOS-Chem species index
       TRCID = Strat_TrID_GC(N)

       ! Skip if species is not defined    
       IF ( TRCID <= 0 ) CYCLE

       ! Get species name
       ThisName = State_Chm%SpcData(TRCID)%Info%Name

       ! ---------------------------------------------------------------
       ! Get pointers to fields
       ! ---------------------------------------------------------------

       ! Production rates [v/v/s]
       IF ( Input_Opt%LUCX ) THEN
          FIELDNAME = 'GMI_PROD_'//TRIM(ThisName)
       ELSE
          FIELDNAME = 'UCX_PROD_'//TRIM(ThisName)
       ENDIF

       ! Get pointer from HEMCO
       CALL HCO_GetPtr( am_I_Root,     HcoState, FIELDNAME, &
                        PLVEC(N)%PROD, RC,       FOUND=FND )

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS .OR. ( PLMUSTFIND .AND. .NOT. FND) ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO! Prod/loss ' // &
                   'data is expected to be listed in the HEMCO '          // &
                   'configuration file. This error occured when trying '  // &
                   'to get field ' // TRIM( FIELDNAME )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Warning message
       IF ( .NOT. FND .AND. AM_I_ROOT ) THEN
          ErrMsg = 'Cannot find archived production rates for '       // &
                    TRIM(ThisName) // ' - will use value of 0.0. '        // &
                   'To use archived rates, add the following field '      // &
                   'to the HEMCO configuration file: '// TRIM( FIELDNAME )
          CALL GC_Warning( ErrMsg, RC, ThisLoc )
       ENDIF

       ! Loss frequency [s-1]
       IF ( Input_Opt%LUCX ) THEN
          FIELDNAME = 'GMI_LOSS_'//TRIM(ThisName)
       ELSE
          FIELDNAME = 'UCX_LOSS_'//TRIM(ThisName)
       ENDIF

       ! Get pointer from HEMCO
       CALL HCO_GetPtr( am_I_Root,     HcoState, FIELDNAME, &
                        PLVEC(N)%LOSS, RC,       FOUND=FND )

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS .OR. ( PLMUSTFIND .AND. .NOT. FND) ) THEN
          ErrMsg = 'Cannot get pointer from HEMCO! Prod/loss ' // &
                   'data is expected to be listed in the HEMCO '          // &
                   'configuration file. This error occured when trying '  // &
                   'to get field ' // TRIM( FIELDNAME )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Warning message
       IF ( .NOT. FND .AND. AM_I_ROOT ) THEN
          ErrMsg= 'Cannot find archived loss frequencies for '        // &
                  TRIM(ThisName) // ' - will use value of 0.0. '          // &
                  'To use archived rates, add the following field '       // &
                  'to the HEMCO configuration file: '//TRIM(FIELDNAME)
          CALL GC_Warning( ErrMsg, RC, ThisLoc )
       ENDIF

    ENDDO !N

    ! Get pointer to STRAT_OH
    CALL HCO_GetPtr( am_I_Root, HcoState, 'STRAT_OH', &
                     STRAT_OH,  RC,        FOUND=FND )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS .OR. .NOT. FND ) THEN
       ErrMsg = 'Cannot find monthly archived strat. OH field '           // &
                '`STRAT_OH`. Please add a corresponding entry to '        // &
                'the HEMCO configuration file.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Set_PLVEC
!EOC
!------------------------------------------------------------------------------
! Prior to 8/9/18:
! Disable CALC_STE routine, which can give misleading STE info (bmy, 8/9/18)
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Calc_Ste
!!
!! !DESCRIPTION: Subroutine CALC\_STE estimates what the stratosphere-to-
!!               troposphere exchange flux must have been since the last time
!!               it was reset
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE Calc_STE( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!!
!! !USES:
!!
!    USE CMN_SIZE_MOD
!    USE ErrCode_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
!    USE Species_Mod,        ONLY : Species
!    USE State_Chm_Mod,      ONLY : ChmState
!    USE State_Met_Mod,      ONLY : MetState
!    USE TIME_MOD,   ONLY : GET_TAU, GET_NYMD, GET_NHMS, EXPAND_DATE
!    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
!    IMPLICIT NONE
!!
!! !INPUT PARAMETERS:
!!
!      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!!
!! !OUTPUT PARAMETERS:
!!
!      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!!
!! !REVISION HISTORY: 
!!  28 Apr 2012 - L. Murray   - Initial version
!!  18 Jul 2012 - R. Yantosca - Make sure I is the innermost DO loop
!!                              (wherever expedient)
!!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!!                              running with the traditional driver main.F
!!  05 Oct 2012 - R. Yantosca - Bug fix for IFORT 12: extend the #if statement
!!                              to avoid including code for nested-grid sims
!!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm, RC arguments
!!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!!  10 Aug 2015 - E. Lundgren - Input tracer concentraiton units are now [kg/kg]
!!  25 May 2016 - E. Lundgren - Replace input_opt%TRACER_MW_KG with species
!!                              database field emMW_g (emitted species g/mol)
!!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!!                              species ID from State_Chm%Map_Advect.
!!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!!  10 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!
!    ! Strings
!    CHARACTER(LEN=255)     :: dateStart, dateEnd
!
!    ! Scalars
!    INTEGER                :: N,         I,       J,    L
!    INTEGER                :: nAdvect,   NA,      NN
!    REAL(fp)               :: dStrat,    STE,     Tend, tauEnd, dt
!    CHARACTER(LEN=63)      :: OrigUnit
!
!    ! Arrays
!    INTEGER                :: LTP(IIPAR,JJPAR      )
!    REAL(fp)               :: M1 (IIPAR,JJPAR,LLPAR)
!    REAL(fp)               :: M2 (IIPAR,JJPAR,LLPAR)
!
!
!    ! Pointers
!    REAL(fp),      POINTER :: Spc(:,:,:,:)
!    TYPE(Species), POINTER :: SpcInfo
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    ! By simple mass balance, dStrat/dt = P - L - STE,
!    ! where STE is the net stratosphere-to-troposphere mass exchange. 
!    !
!    ! Therefore, we estimate STE as
!    !   STE = (P-L) - dStrat/dt
!    !
!    ! As the tropopause is dynamic, we use the mean tropopause level during
!    ! the period for determining initial and end stratospheric masses. 
!    ! (ltm, 04/28/2012)
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    ! Assume success
!    RC      = GC_SUCCESS
!
!    ! Number of advected species
!    nAdvect = State_Chm%nAdvect
!
!#if defined( NESTED_NA ) || defined( NESTED_CH ) || defined( NESTED_EU ) || defined( NESTED_AS )
!    ! This method only works for a global domain.
!    ! It could be modified for nested domains if the total mass flux across the
!    ! boundaries during the period is taken into account.
!    RETURN
!#else
!
!    ! Convert State_Chm%SPECIES to [kg] (ewl, 8/10/15)
!    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
!                            State_Chm, 'kg', RC, OrigUnit=OrigUnit )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error('Unit conversion error', RC, &
!                     'Calc_STE in strat_chem_mod.F')
!       RETURN
!    ENDIF  
!
!    ! Point to chemical species array [kg]
!    Spc => State_Chm%Species
!
!    ! Determine mean tropopause level for the period
!    !$OMP PARALLEL DO                               &
!    !$OMP DEFAULT( SHARED )                         &
!    !$OMP PRIVATE( I,  J  )
!    DO J = 1,JJPAR
!    DO I = 1,IIPAR
!       LTP(I,J) = NINT( TPauseL(I,J) / TPauseL_Cnt )
!    ENDDO
!    ENDDO
!    !$OMP END PARALLEL DO
!
!    ! Period over which STE is being determined [a]
!    tauEnd = GET_TAU() ! [h]
!    dt = ( tauEnd - tauInit ) / 24e+0_fp / 365.25e+0_fp
!
!    dateStart = 'YYYY-MM-DD hh:mm'
!    CALL EXPAND_DATE(dateStart,NymdInit,NhmsInit)
!    dateEnd = 'YYYY-MM-DD hh:mm'
!    CALL EXPAND_DATE(dateEnd,GET_NYMD(),GET_NHMS())
!
!    ! Print to output
!    IF ( am_I_Root ) THEN
!       WRITE( 6, * ) ''
!       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!       WRITE( 6, '(a)' ) '  Strat-Trop Exchange'
!       WRITE( 6, '(a)' ) REPEAT( '-', 79 )
!       WRITE( 6, '(a)' ) &
!            '  Global stratosphere-to-troposphere fluxes estimated over'
!       WRITE( 6, 100 ) TRIM(dateStart), TRIM(dateEnd)
!       WRITE( 6, * ) ''
!       WRITE( 6, 110 ) 'Species','[moles a-1]','* [g/mol]','= [Tg a-1]'
!    ENDIF
!100 FORMAT( 2x,a16,' to ',a16 )
!110 FORMAT( 2x,a8,':',4x,a11  ,4x,a9  ,4x,  a11 )
!
!    ! Loop over only the advected species
!    DO NA = 1, nAdvect
!
!       ! Get the species ID from the advected species ID
!       N = State_Chm%Map_Advect(NA)
!    
!       ! Populate before (M1) and after (M2) state for the species [kg]
!       M1 = MInit(:,:,:,NA)             
!       M2 =   Spc(:,:,:,N )
!
!       ! Zero out troposphere and determine total change in the stratospheric
!       ! burden of species N (dStrat) [kg]
!       !$OMP PARALLEL DO &
!       !$OMP DEFAULT( SHARED ) &
!       !$OMP PRIVATE( I,  J  )
!       DO J = 1, JJPAR  
!       DO I = 1, IIPAR
!          M2(I,J,1:LTP(I,J)) = 0e+0_fp
!          M1(I,J,1:LTP(I,J)) = 0e+0_fp
!       ENDDO
!       ENDDO
!       !$OMP END PARALLEL DO
!       dStrat   = SUM(M2)-SUM(M1)
!
!       ! The total chemical tendency (P-L) over the period for species N [kg]
!       ! Schem_tend is computed in this module for tropchem simulations and
!       ! in flexchem_mod.F90 for UCX simulations
!       Tend   = SUM(Schem_tend(:,:,:,N))
!
!       ! Calculate flux as STE = (P-L) - dStrat/dt
!       STE = (Tend-dStrat)/dt ! [kg a-1]
!
!       ! Get info about this species from the species database
!       SpcInfo => State_Chm%SpcData(N)%Info
!
!       ! Print to standard output
!       IF ( am_I_Root ) THEN
!          WRITE(6,120) TRIM( SpcInfo%Name ),                  &
!               STE / (SpcInfo%emMW_g * 1.e-3_fp),             & ! mol/a-1
!               SpcInfo%emMW_g,                                & ! g/mol
!               STE * 1e-9_fp                                  ! Tg a-1
!       ENDIF 
!    ENDDO
!
!120 FORMAT( 2x,a8,':',4x,e11.3,4x,f9.1,4x,f11.4 )
!
!    IF ( am_I_Root ) THEN
!       WRITE( 6, * ) ''
!       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!       WRITE( 6, * ) ''
!    ENDIF
!
!    ! Reset variables for next STE period
!    NymdInit             = GET_NYMD()
!    NhmsInit             = GET_NHMS()
!    TauInit              = GET_TAU()
!    TPauseL_Cnt          = 0e+0_fp
!    TPauseL(:,:)         = 0e+0_fp
!    SChem_tend(:,:,:,:)  = 0e+0_fp
!
!    ! Loop over only the advected species
!    DO NA = 1, nAdvect
!
!       ! Get the species ID from the advected species ID
!       N = State_Chm%Map_Advect(NA)
!       
!       ! Reset MINIT for next STE period
!       MInit(:,:,:,NA)  = Spc(:,:,:,N)
!
!    ENDDO
!
!    ! Convert species units back to original unit 
!    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
!                            State_Chm, OrigUnit,  RC )
!    IF ( RC /= GC_SUCCESS ) THEN
!       CALL GC_Error('Unit conversion error', RC, &
!                     'Calc_STE in strat_chem_mod.F')
!       RETURN
!    ENDIF  
!
!    ! Free pointers
!    Spc     => NULL()
!    SpcInfo => NULL()
!#endif
!  END SUBROUTINE Calc_STE
!!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Strat_Chem
!
! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.  
!  It also opens the necessary rate files.
!\\
!\\
! !INTERFACE:
!      
  SUBROUTINE INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TAU
    USE TIME_MOD,           ONLY : GET_NYMD
    USE TIME_MOD,           ONLY : GET_NHMS
    USE TIME_MOD,           ONLY : GET_TS_CHEM

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  01 Feb 2011 - L. Murray   - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  26 Oct 2012 - R. Yantosca - Now pass Chemistry State object for GIGC
!  09 Nov 2012 - R. Yantosca - Now pass Input Options object for GIGC
!  05 Nov 2013 - R. Yantosca - Now update tracer flags for tagOx simulation
!  03 Apr 2014 - R. Yantosca - PROD, LOSS, STRAT_OH, MINIT, SCHEM_TEND are 
!                              now REAL*4, so use 0e0 to initialize
!  11 Aug 2015 - E. Lundgren - Tracer units are now kg/kg and are converted 
!                              kg for assignment of MInit
!  16 Jun 2016 - M. Yannetti - Replaced TRACERID_MOD.
!  20 Jun 2016 - R. Yantosca - Now save species ID flags as module variables
!                              and only define them in the INIT phase.
!  12 Jul 2016 - R. Yantosca - Now also store advected species ID's
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  20 Sep 2016 - R. Yantosca - Rewrote GMI_TrName statement for Gfortran
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=16)      :: sname
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

    ! Scalars
    INTEGER                :: AS, N, NN, NA, nAdvect
    LOGICAL                :: IT_IS_A_FULLCHEM_SIM
    LOGICAL                :: IT_IS_A_TAGO3_SIM
    LOGICAL                :: LLINOZ
    LOGICAL                :: LUCX

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !================================================================
    ! Initialize 
    !=================================================================
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> Init_Strat_Chem (in Headers/strat_chem_mod.F90)'

    !=================================================================
    ! Get species ID flags
    !=================================================================
    id_Br                    = Ind_('Br'     )
    id_Br2                   = Ind_('Br2'    )
    id_BrNO3                 = Ind_('BrNO3'  ) 
    id_BrO                   = Ind_('BrO'    )
    id_CHBr3                 = Ind_('CHBr3'  )
    id_CH2Br2                = Ind_('CH2Br2' ) 
    id_CH3Br                 = Ind_('CH3Br'  )
    id_HOBr                  = Ind_('HOBr'   )
    id_HBr                   = Ind_('HBr'    )
    id_O3                    = Ind_('O3'     )
    id_O3Strat               = Ind_('O3Strat')   


    ! Save fields from the Input_Opt object to local variables
    LLINOZ                   = Input_Opt%LLINOZ
    LUCX                     = Input_Opt%LUCX
    IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGO3_SIM        = Input_Opt%ITS_A_TAGO3_SIM

    ! Number of advected species
    nAdvect                  = State_Chm%nAdvect

    ! Initialize counters, initial times, mapping arrays
    TpauseL_Cnt              = 0.e+0_fp
    NSCHEM                   = 0
    TauInit                  = GET_TAU()
    NymdInit                 = GET_NYMD()
    NhmsInit                 = GET_NHMS()

    ! Initialize timestep for chemistry [s]
    dTchem = GET_TS_CHEM()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Allocate and initialize arrays for mapping UCX/GMI species
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! At most NTR species could overlap between GC & UCX/GMI
    IF ( LUCX ) THEN
       NTR = 125 ! Number of GMI species
    ELSE
       NTR = 145 ! Number of UCX species
    ENDIF

    ALLOCATE( TrName(NTR), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TrName' )
    TrName(:) = ''

    ALLOCATE( Strat_TrID_GC (NTR), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Strat_TrID_GC' )
    Strat_TrID_GC(:) = 0

    ALLOCATE( Strat_TrId_TND(NTR), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Strat_TrId_TND' )
    Strat_TrID_TND(:) = 0

    ALLOCATE( Strat_TrID(NTR), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Strat_TrID' )
    Strat_TrID(:) = 0

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Determine the mapping for UCX/GMI to GC variables based on
    ! species name, which only needs to be done once per model run.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! List of available species with archived monthly climatological
    ! production rates, loss frequencies, and mixing ratios 
    !
    ! Rewrote the TrName array constructor so that all of the species
    ! names have the same length.  This will prevent Gfortran from choking 
    ! with an error.  This is OK since we trim TrName before using
    ! it in any string comparisons. (bmy, 9/20/16)
    IF ( LUCX ) THEN

       ! Species in GMI Combo model
       TrName = (/ 'A3O2    ', 'ACET    ', 'ACTA    ', 'ALD2    ',   &
                   'ALK4    ', 'ATO2    ', 'B3O2    ', 'Br      ',   &
                   'BrCl    ', 'BrO     ', 'BrONO2  ', 'C2H6    ',   &
                   'C3H8    ', 'CCl4    ', 'CF2Br2  ', 'CF2Cl2  ',   &
                   'CF2ClBr ', 'CF3Br   ', 'CFC113  ', 'CFC114  ',   &
                   'CFC115  ', 'CFCl3   ', 'CH2O    ', 'CH3Br   ',   &
                   'CH3CCl3 ', 'CH3Cl   ', 'CH4     ', 'CO      ',   &
                   'Cl      ', 'Cl2     ', 'Cl2O2   ', 'ClO     ',   &
                   'ClONO2  ', 'EOH     ', 'ETO2    ', 'ETP     ',   &
                   'GCO3    ', 'GLYC    ', 'GLYX    ', 'GP      ',   &
                   'GPAN    ', 'H       ', 'H2      ', 'H2402   ',   &
                   'H2O     ', 'H2O2    ', 'HAC     ', 'HBr     ',   &
                   'HCFC141b', 'HCFC142b', 'HCFC22  ', 'HCOOH   ',   &
                   'HCl     ', 'HNO2    ', 'HNO3    ', 'HNO4    ',   &
                   'HO2     ', 'HOBr    ', 'HOCl    ', 'IALD    ',   &
                   'IAO2    ', 'IAP     ', 'INO2    ', 'INPN    ',   &
                   'ISN1    ', 'ISNP    ', 'ISOP    ', 'KO2     ',   &
                   'MACR    ', 'MAN2    ', 'MAO3    ', 'MAOP    ',   &
                   'MAP     ', 'MCO3    ', 'MEK     ', 'MGLY    ',   &
                   'MO2     ', 'MOH     ', 'MP      ', 'MRO2    ',   &
                   'MRP     ', 'MVK     ', 'MVN2    ', 'N       ',   &
                   'N2O     ', 'N2O5    ', 'NO      ', 'NO2     ',   &
                   'NO3     ', 'NOx     ', 'O       ', 'O1D     ',   &
                   'O3      ', 'OClO    ', 'OH      ', 'Ox      ',   &
                   'PAN     ', 'PMN     ', 'PO2     ', 'PP      ',   &
                   'PPN     ', 'PRN1    ', 'PRPE    ', 'PRPN    ',   &
                   'R4N1    ', 'R4N2    ', 'R4O2    ', 'R4P     ',   &
                   'RA3P    ', 'RB3P    ', 'RCHO    ', 'RCO3    ',   &
                   'RCOOH   ', 'RIO1    ', 'RIO2    ', 'RIP     ',   &
                   'ROH     ', 'RP      ', 'VRO2    ', 'VRP     ',   &
                   'RIPA    ', 'RIPB    ', 'RIPD    ', 'NPMN    ',   &
                   'IPMN    ' /)

    ELSE

       ! Species in GEOS-Chem UCX simulation
       TrName = (/ 'Ox      ', 'NO      ', 'O3      ', 'PAN     ',   & 
                   'CO      ', 'ALK4    ', 'ISOP    ', 'HNO3    ',   &
                   'H2O2    ', 'ACET    ', 'MEK     ', 'ALD2    ',   &
                   'RCHO    ', 'MVK     ', 'MACR    ', 'NPMN    ',   &
                   'PPN     ', 'R4N2    ', 'PRPE    ', 'C3H8    ',   &
                   'CH2O    ', 'C2H6    ', 'N2O5    ', 'HNO4    ',   &
                   'MP      ', 'DMS     ', 'SO2     ', 'SO4     ',   &
                   'MSA     ', 'Br2     ', 'Br      ', 'BrO     ',   &
                   'HOBr    ', 'HBr     ', 'BrNO2   ', 'BrNO3   ',   &
                   'CHBr3   ', 'CH2Br2  ', 'CH3Br   ', 'MPN     ',   &
                   'ISOPND  ', 'ISOPNB  ', 'MOBA    ', 'PROPNN  ',   &
                   'HAC     ', 'GLYC    ', 'MVKN    ', 'MACRN   ',   &
                   'MAP     ', 'NO2     ', 'NO3     ', 'HNO2    ',   &
                   'HNO2    ', 'BENZ    ', 'TOLU    ', 'XYLE    ',   &
                   'MTPA    ', 'LIMO    ', 'EOH     ', 'MGLY    ',   &
                   'GLYX    ', 'ACTA    ', 'HPALD   ', 'DHDN    ',   &
                   'ETHLN   ', 'HCOOH   ', 'IEPOXA  ', 'IEPOXB  ',   &
                   'IEPOXD  ', 'ISN1    ', 'RIPA    ', 'RIPB    ',   &
                   'RIPD    ', 'IMAE    ', 'SOAIE   ', 'SOAME   ',   &
                   'SOAGX   ', 'SOAMG   ', 'LVOC    ', 'LVOCOA  ',   &
                   'ISN1OG  ', 'ISN1OA  ', 'MONITS  ', 'MONITU  ',   &
                   'HONIT   ', 'IONITA  ', 'MONITA  ', 'INDIOL  ',   &
                   'IPMN    ', 'HC187   ', 'N2O     ', 'OCS     ',   &
                   'CH4     ', 'BrCl    ', 'HCl     ', 'CCl4    ',   &
                   'CH3Cl   ', 'CH3CCl3 ', 'CFC113  ', 'CFC114  ',   &
                   'CFC115  ', 'HCFC123 ', 'HCFC141 ', 'HCFC142 ',   &
                   'CFC11   ', 'CFC12   ', 'HCFC22  ', 'H1211   ',   &
                   'H1301   ', 'H2402   ', 'Cl      ', 'ClO     ',   &
                   'HOCl    ', 'ClNO3   ', 'ClNO2   ', 'ClOO    ',   &
                   'OClO    ', 'Cl2     ', 'Cl2O2   ', 'H2O     ',   &
                   'BrSALA  ', 'BrSALC  ', 'CHCl3   ', 'CH2Cl2  ',   &
                   'CH3I    ', 'CH2I2   ', 'CH2ICl  ', 'CH2IBr  ',   &
                   'HOI     ', 'I2      ', 'IBr     ', 'ICl     ',   &
                   'I       ', 'IO      ', 'HI      ', 'OIO     ',   &
                   'INO     ', 'IONO    ', 'IONO2   ', 'I2O2    ',   &
                   'I2O3    ', 'I2O4    ', 'ISALA   ', 'ISALC   ',   &
                   'AERI    ' /)

    ENDIF

    !===========================!
    ! Full chemistry simulation !
    !===========================!
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       DO NN = 1, NTR

          sname = TRIM(TrName(NN))
 
          ! Loop over only the advected species
          DO NA = 1, nAdvect
          
             ! Get the species ID from the advected species ID
             N       = State_Chm%Map_Advect(NA)

             ! Get the corresponding entry in the species database
             SpcInfo => State_Chm%SpcData(N)%Info

             ! For now, guarantee that GMI/UCX prod/loss rates are not used for
             ! any bromine species
             IF ( TRIM( SpcInfo%Name ) .eq.      'Br' .or. &
                  TRIM( SpcInfo%Name ) .eq.    'BrCl' .or. &
                  TRIM( SpcInfo%Name ) .eq.     'BrO' .or. &
                  TRIM( SpcInfo%Name ) .eq.  'BrONO2' .or. &
                  TRIM( SpcInfo%Name ) .eq.  'CF2Br2' .or. &
                  TRIM( SpcInfo%Name ) .eq. 'CF2ClBr' .or. &
                  TRIM( SpcInfo%Name ) .eq.   'CF3Br' .or. &
                  TRIM( SpcInfo%Name ) .eq.   'CH3Br' .or. &
                  TRIM( SpcInfo%Name ) .eq.     'HBr' .or. &
                  TRIM( SpcInfo%Name ) .eq.    'HOBr'        ) CYCLE

             ! SDE 08/28/13: Full strat. has its own mesospheric NOy handling
             IF ( LUCX ) THEN
                IF ( TRIM( SpcInfo%Name ) .eq.    'NO' .or. &
                     TRIM( SpcInfo%Name ) .eq.   'NO2' .or. &
                     TRIM( SpcInfo%Name ) .eq.   'NO3' .or. &
                     TRIM( SpcInfo%Name ) .eq.   'NOx' .or. &
                     TRIM( SpcInfo%Name ) .eq.     'N' .or. &
                     TRIM( SpcInfo%Name ) .eq.   'N2O' ) CYCLE
             ENDIF

             IF ( TRIM( SpcInfo%Name ) .eq. TRIM(sname) ) THEN
                
                IF ( LLINOZ .and. TRIM( SpcInfo%Name ) .eq. 'O3' ) THEN
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM( SpcInfo%Name ) // ' (via Linoz)'
                   ENDIF
                ELSE IF ( Input_Opt%LSYNOZ .AND. TRIM( SpcInfo%Name ) .eq. 'O3' ) THEN
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM( SpcInfo%Name ) // ' (via Synoz)'
                   ENDIF
                ELSE
                   IF ( am_I_Root ) THEN
                    IF ( LUCX ) THEN
                      WRITE( 6, '(a)' ) TRIM( SpcInfo%Name )//' (via GMI rates)'
                    ELSE
                      WRITE( 6, '(a)' ) TRIM( SpcInfo%Name )//' (via UCX rates)'
                    ENDIF
                   ENDIF
                ENDIF

                NSCHEM                 = NSCHEM + 1
                Strat_TrID_GC (NSCHEM) = N  ! Maps 1:NSCHEM to species array
                Strat_TrId_TND(NSCHEM) = NA ! Maps 1:NSCHEM to SCHEM_TEND
                Strat_TrID    (NSCHEM) = NN ! Maps 1:NSCHEM to TrName index
             ENDIF

             ! Free pointer
             SpcInfo => NULL()
          ENDDO
       ENDDO

       ! These are the reactions with which we will use OH fields
       ! to determine stratospheric loss.
       IF( am_I_Root ) THEN
          IF ( id_CHBr3  .gt. 0 ) WRITE(6,*) 'CHBr3 (from GMI OH)'
          IF ( id_CH2Br2 .gt. 0 ) WRITE(6,*) 'CH2Br2 (from GMI OH)'
          IF ( id_CH3Br  .gt. 0 ) WRITE(6,*) 'CH3Br (from GMI OH)'
       ENDIF

    !===========!
    ! Tagged O3 !
    !===========!
    ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN
       IF ( LLINOZ ) THEN
          IF ( am_I_Root ) THEN
             WRITE(6,*) 'Linoz ozone performed on: '
          ENDIF
       ELSE          
          IF ( am_I_Root ) THEN
             WRITE(6,*) 'Synoz ozone performed on: '
          ENDIF
       ENDIF

       ! Loop over only the advected species
       DO NA = 1, nAdvect
          
          ! Get the species ID from the advected species ID
          N       = State_Chm%Map_Advect(NA)

          ! Get the corresponding entry in the species database
          SpcInfo => State_Chm%SpcData(N)%Info

          IF ( TRIM( SpcInfo%Name ) .eq. 'O3'        .or. &
               TRIM( SpcInfo%Name ) .eq. 'O3Strt'    .or. &
               TRIM( SpcInfo%Name ) .eq. 'O3Strat' ) THEN
             NSCHEM                 = NSCHEM + 1    ! Increment count
             Strat_TrID_GC(NSCHEM)  = N             ! Use for Sc%Species
             Strat_TrID_TND(NSCHEM) = NA            ! Use for SCHEM_TEND
             IF ( am_I_Root ) THEN
                WRITE(6,*) TRIM( SpcInfo%Name )
             ENDIF
          ENDIF

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Allocate and initialize prod/loss vector
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ALLOCATE( PLVEC( NSCHEM ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'PLVEC' )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Allocate and initialize arrays for STE calculation !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    ! Array to hold initial state of atmosphere at the beginning
    ! of the period over which to estimate STE. Populate with
    ! initial atm. conditions from restart file converted to [kg/kg].
    ALLOCATE( MInit( IIPAR, JJPAR, LLPAR, nAdvect ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MInit' )
    MInit = 0.0_fp

#if !defined( ESMF_ ) && !defined( MODEL_WRF )
    ! Set MINIT. Ignore in ESMF environment because State_Chm%Species
    ! is not yet filled during initialization. (ckeller, 4/6/16)
    CALL SET_MINIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Set_MInit"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    ! Array to determine the mean tropopause level over the period
    ! for which STE is being estimated.
    ALLOCATE( TPAUSEL( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPAUSEL' )
    TPAUSEL = 0e+0_fp

    ! Array to aggregate the stratospheric chemical tendency [kg period-1]
    ALLOCATE( SCHEM_TEND(IIPAR,JJPAR,LLPAR,nAdvect), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCHEM_TEND' )
    SCHEM_TEND = 0e0

    ! Allocate and initialize bromine arrays
    GC_Bry_TrID(1:6) = (/id_Br2,id_Br,id_BrO,id_HOBr,id_HBr,id_BrNO3/)

  END SUBROUTINE INIT_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Minit 
!
! !DESCRIPTION: Sets the MINIT array to current values in State\_Chm%Species.
!\\
!\\
! !INTERFACE:
!      
  SUBROUTINE SET_MINIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
    USE UnitConv_Mod,  ONLY : Convert_Spc_Units

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  06 Apr 2016 - C. Keller   - Initial version: moved outside of 
!                              INIT_STRAT_CHEM so that it can be called from
!                              within DO_STRAT_CHEM (for ESMF applications).
!  15 Mar 2017 - E. Lundgren - Add unit catch and error handling
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, NA
    LOGICAL            :: UNITCHANGE

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    CHARACTER(LEN=63)  :: OrigUnit

    !=================================================================
    ! SET_MINIT begins here!
    !=================================================================
    
    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = '-> Set_MInit (in GeosCore/strat_chem_mod.F90)'

    ! Assume no unit change is necessary and species are already in kg
    UNITCHANGE = .FALSE.

    ! Safety check that the advected species are non-zero.
    IF ( SUM( State_Chm%Species(:,:,:,1:State_Chm%nAdvect) ) <= 0.0_fp ) THEN
       ErrMsg = 'One or more advected species have negative concentrations!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Convert species to [kg] so that initial state of atmosphere is 
    ! in same units as chemistry (ewl, 8/10/15)
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'kg', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF 

    ! Loop over only the advected species
    DO NA = 1, State_Chm%nAdvect

       ! Get the species ID from the advected species ID
       N = State_Chm%Map_Advect(NA)

       ! Save initial conditions in MINIT
       MInit(:,:,:,NA) = State_Chm%Species(:,:,:,N)

    ENDDO

    ! Convert species units back to original unit 
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF  

    ! Minit is now set
    MInit_Is_Set = .TRUE.

  END SUBROUTINE SET_MINIT 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Strat_Chem
!
! !DESCRIPTION: Subroutine CLEANUP\_STRAT\_CHEM deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_STRAT_CHEM
!
! !USES:
!
    IMPLICIT NONE
!
! !REVISION HISTORY: 
!  01 Feb 2011 - L. Murray   - Initial version
!  03 Oct 2016 - R. Yantosca - Deallocate BrPtrDay and BrPtrNight
!  05 Oct 2016 - R. Yantosca - Now make deallocations more robust
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I

    ! Deallocate arrays
    IF ( ALLOCATED( MInit          ) ) DEALLOCATE( MInit          )
    IF ( ALLOCATED( TPAUSEL        ) ) DEALLOCATE( TPAUSEL        )
    IF ( ALLOCATED( SCHEM_TEND     ) ) DEALLOCATE( SCHEM_TEND     )
    IF ( ALLOCATED( TrName         ) ) DEALLOCATE( TrName         )
    IF ( ALLOCATED( Strat_TrID_GC  ) ) DEALLOCATE( Strat_TrID_GC  )
    IF ( ALLOCATED( Strat_TrID_TND ) ) DEALLOCATE( Strat_TrID_TND )
    IF ( ALLOCATED( Strat_TrID     ) ) DEALLOCATE( Strat_TrID     )

    ! Cleanup BrPtrDay array of pointers
    IF ( ASSOCIATED( BrPtrDay ) ) THEN
       DO I = 1, 6
          BrPtrDay(I)%MR => NULL()
       ENDDO
       DEALLOCATE( BrPtrDay )
    ENDIF

    ! Cleanup BrPtrNight array of pointers
    IF ( ASSOCIATED( BrPtrNight ) ) THEN
       DO I = 1, 6
          BrPtrNight(I)%MR => NULL()
       ENDDO
       DEALLOCATE( BrPtrNight )
    ENDIF

    ! Cleanup PLVEC array of pointers
    IF ( ASSOCIATED( PLVEC ) ) THEN
       DO I = 1, NSCHEM
          PLVEC(I)%PROD => NULL()
          PLVEC(I)%LOSS => NULL()
       ENDDO
       DEALLOCATE( PLVEC )
    ENDIF

    ! Cleanup pointer to strat. OH
    STRAT_OH => NULL()

  END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Synoz
!
! !DESCRIPTION: Subroutine Do\_Synoz establishes the flux boundary condition 
!  for Ozone coming down from the stratosphere, using the Synoz algorithm of
!  McLinden et al, 2000. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Synoz( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_YEAR

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
! !REMARKS:
!  Reference:
!  ============================================================================
!  C. A. McLinden, S. Olsen, B. Hannegan, O. Wild, M. J. Prather, and
!  J. Sundet, "Stratospheric Ozone in 3-D models: A simple chemistry
!  and the cross-tropopause flux".
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%% NOTE: This SYNOZ scheme is now obsolete, replaced by LINOZ   %%%%
!  %%%% We keep this for backwards compatibility w/ older met fields %%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! !REVISION HISTORY:
!  13 Dec 1999 - Q. Li, R. Martin - Initial version
!  (1 ) The parameter Rdg0 from "CMN_GCTM" = R / g0 = 28.97. 
!  (2 ) Pass PW = PS - PTOP to UPBDFLX via "CMN".
!  (3 ) Now pass IORD, JORD, KORD as arguments (bmy, 12/6/00)
!  (4 ) Now compute the proper value of PO3_vmr that will yield 475 Tg O3/yr
!        for various settings of IORD, JORD, KORD (rvm, bey, bmy, 12/5/00)
!                                                                             .
!        **************************************************************
!        ***** You must use this version of UPBDFLX_O3 if you are *****
!        ***** using the Parallel Processor TPCORE v. 7.1         *****
!        **************************************************************
!                                                                             .
!  (5 ) Added to "upbdflx_mod.f".  Also updated comments and made some
!        cosmetic changes. (bmy, 6/28/01)
!  (6 ) Now reference CMN_SETUP for LSPLIT.  Also store strat O3 into
!        tracer #11 for multi-tracer Ox run. (amf, bmy, 7/3/01)
!  (7 ) Removed IREF, JREF -- these are obsolete.  Also T(IREF,JREF,L) is
!        now T(I,J,L). (bmy, 9/27/01)
!  (8 ) Also replace PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (9 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (10) Removed obsolete commented out code from 7/01 (bmy, 11/26/01)
!  (11) Now write file names to stdout (bmy, 4/3/02)
!  (12) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (13) Now use GET_PEDGE and GET_PCENTER from "pressure_mod.f" to compute
!        the pressure at the bottom edge and center of grid box (I,J,L).
!        Also removed obsolete, commented-out code.  Removed G_SIG and
!        G_SIGE from the arg list. (dsa, bdf, bmy, 8/21/02)
!  (14) Now reference BXHEIGHT and T from "dao_mod.f".  Also reference routine
!        ERROR_STOP from "error_mod.f".  Now references IDTOX from F90 module
!        "tracerid_mod.f" instead of from "comtrid.h". (bmy, 11/6/02)
!  (15) Now define J30S and J30N for 1x1 nested grid (bmy, 3/11/03)
!  (16) Make sure to pass AD via "dao_mod.f" for GEOS-1 (bnd, bmy, 4/14/03)
!  (17) On the first timestep, print how much O3 flux is coming down from the 
!        stratosphere in Tg/yr. (mje, bmy, 8/15/03)
!  (18) Change O3 flux to 500 Tg/yr for GEOS-3 (mje, bmy, 9/15/03)
!  (19) Now calls routine ADD_STRAT_POX from "tagged_ox_mod.f" in order to
!        pass stratospheric flux of Ox to the proper tagged tracer w/o
!        resorting to hardwiring w/in this routine. (bmy, 8/18/03)
!  (20) Add GEOS_4 to the #if defined block. (bmy, 1/29/04)
!  (21) Activated parallel DO-loops.  Now made STFLUX a local array
!        in order to facilitate parallelization. (bmy, 4/15/04)
!  (22) Removed IORD, JORD, KORD from the arg list.  Now reference STT and
!        ITS_A_TAGOX_SIM from "tracer_mod.f".  (bmy, 7/20/04)
!  (23) Use an #ifdef block to comment out an EXIT statement from w/in a
!        parallel loop for COMPAQ compiler.  COMPAQ seems to have some 
!        problems with this.  Now supports 1x125 grid. (auvray, bmy, 12/1/04)
!  (24) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (25) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (26) Now set J30S and J30N for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (27) Remove support for COMPAQ compiler (bmy, 7/8/09)
!  (28) Now do not call ADD_STRAT_POx for tagged Ox (dbj, bmy, 10/16/09)
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5 (bmy, 8/13/10)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  10 Feb 2012 - R. Yantosca - Modified for 0.25 x 0.3125 grids
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  28 Apr 2012 - L. Murray - Moved from upbdflx_mod.F to here, modified to 
!                 F90, renamed from UPBDFLX_O3 to DO_SYNOZ. Use chem timestep
!                 now. Also, removed INIT_UPBDFLX, which was last used for 
!                 GEOS-3.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  04 Feb 2013 - M. Payer    - Replace all JJPAR with values for nested grids
!                              since JJPAR is no longer a parameter
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  25 Mar 2013 - R. Yantosca - Now use explicit numbers for J30S, J30N
!  31 May 2013 - R. Yantosca - Now pass Input_Opt, RC as arguments
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocessor switch
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  05 Nov 2013 - R. Yantosca - Rename IDTOxStrt to id_O3Strat
!  23 Jan 2014 - M. Sulprizio- Linoz does not call UPBDFLX_O3. Synoz does. 
!                              Now uncomment ADD_STRAT_POx (jtl,hyl,dbj,11/3/11)
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE and GET_PCENTER with
!                              State_Met%PEDGE and State_Met%PMID. Remove
!                              dependency on pressure_mod.
!  03 Mar 2015 - E. Lundgren - Use virtual temperature in hypsometric eqn
!  12 Aug 2015 - R. Yantosca - Add placeholder values for 0.5 x 0.625 grids
!  16 Jun 2016 - M. Yannetti - Replaced TRACERID_MOD.
!  20 Jun 2016 - R. Yantosca - Now make species ID flags module variables
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  12 Jul 2016 - R. Yantosca - Remove references to ADD_STRAT_POX
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE        :: FIRST = .TRUE.
    INTEGER              :: I, J, L, L70mb, N, NA, nAdvect
    REAL(fp)             :: P1, P2, P3, T1, T2, DZ, ZUP
    REAL(fp)             :: DTCHEM, H70mb, PO3, PO3_vmr
    REAL(fp)             :: STFLUX(IIPAR,JJPAR,LLPAR)

    ! Pointers
    REAL(fp), POINTER    :: Spc(:,:,:,:)
!
! !DEFINED PARAMETERS:
!
    ! Select the grid boxes at the edges of the O3 release region, 
    ! for the proper model resolution (qli, bmy, 12/1/04)
#if defined( GRID4x5 )
    INTEGER, PARAMETER   :: J30S = 16, J30N = 31 

#elif defined( GRID2x25 ) 
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#elif defined( GRID05x0625 )

    !%%% ADD PLACEHOLDER VALUES, THESE AREN'T REALLY USED ANYMORE %%%
#if   defined( NESTED_AS )
    INTEGER, PARAMETER :: J30S = 1,  J30N = 83
#elif defined( NESTED_NA )
    INTEGER, PARAMETER :: J30S = 1,  J30N = 41
#elif defined( NESTED_EU )
    INTEGER, PARAMETER :: J30S = 1,  J30N = 1  ! add later-checked . it is ok Anna Prot
#endif

#elif defined( GRID025x03125 )

#if   defined( NESTED_CH )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 161
#elif defined( NESTED_NA )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 161 !I think it should be 202/Anna Prot
!Anna Prot added 8 May 2015
#elif defined( NESTED_EU )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 115
#endif

#elif defined( EXTERNAL_GRID )
    ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

    ! Lower pressure bound for O3 release (unit: mb)
    ! REAL(fp),  PARAMETER   :: P70mb = 70e+0_fp !PHS
    REAL(fp)             :: P70mb, PTP

    !=================================================================
    ! Do_Synoz begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Chemical timestep [s]
    ! Originally, Synoz was in transport code, and used dynamic dT.
    DTCHEM = GET_TS_CHEM()

    ! For O3 flux printout
    STFLUX = 0e+0_fp

    ! lower pressure !PHS
    P70mb = 70e+0_fp

    ! Point to chemical species array [v/v dry]
    Spc => State_Chm%Species

    !=================================================================
    ! Compute the proper release rate of O3 coming down from the 
    ! stratosphere for the different GEOS model grids. 
    ! (bey, djj, rvm, bmy, 12/5/00).
    !
    ! PO3_vmr is the O3 release rate constant [v/v/s] that will yield 
    ! a strat-to-trop flux of 475 [Tg O3/yr].  Different TPCORE flags 
    ! create different amounts of ozone in the stratosphere.  Flags 
    ! 337F are currently preferred (bey, djj, rvm).  
    !  
    ! For now, provide values for PO3_vmr for two TPCORE flag settings:
    !  (1) IORD = 3, JORD = 3, KORD = 7  (preferred, assumed to 
    !                                     be the default)
    !  (2) IORD = 5, JORD = 5, KORD = 7 
    !=================================================================
#if defined( GISS ) && defined( MODELE )

    ! For Model E, assuming 3,3,7 and 475 Tg N a-1
    PO3_vmr = 4.84610e-14 !/ 2e+0_fp

    ! Scale as necessary for PreIndustrial and Paleo Climates
    ! These values determined by Linoz simulations run for each climate
    if ( GET_YEAR() .ge. 2102 .and. GET_YEAR() .le. 2105 ) then
       ! PIH was 3% higher
       PO3_vmr = PO3_vmr * 1.0273853e+0_fp
    endif
    if ( GET_YEAR() .ge. 2202 .and. GET_YEAR() .le. 2205 ) then
       ! LGM Webb STE was 7% higher
       PO3_vmr = PO3_vmr * 1.0723525e+0_fp
    endif
    if ( GET_YEAR() .ge. 2302 .and. GET_YEAR() .le. 2305 ) then
       ! LGM CLIMAP STE was 3% higher
       PO3_vmr = PO3_vmr * 1.0285232e+0_fp
    endif

#else

    PO3_vmr = 5.14e-14_fp   

#endif

    ! Only initialize on first time step
    IF ( FIRST ) STFLUX = 0e+0_fp

    ! Loop over latitude (30S -> 30N) and longitude
    !$OMP PARALLEL DO                               &
    !$OMP DEFAULT( SHARED )                         &
    !$OMP PRIVATE( I,  J,  L,  P2,  L70mb, P1, P3 ) &
    !$OMP PRIVATE( T2, T1, DZ, ZUP, H70mb, PO3    )
    DO J = J30S, J30N 
       DO I = 1,    IIPAR

          DO L = 1, LLPAR

             ! P2 = pressure [hPa] at the sigma center of level L70mb
             P2 = State_Met%PMID(I,J,L) 

             ! L70mb is the 1st layer where pressure is equal to
             ! or smaller than 70 mb 
             IF ( P2 < P70mb ) THEN
                L70mb = L
                EXIT
             ENDIF
          ENDDO

          ! P1 = pressure [hPa] at the sigma center of level L70mb - 1   
          P1 = State_Met%PMID(I,J,L70mb-1) 

          ! P3 = pressure [hPa] at the lower sigma edge of level L70mb
          P3 = State_Met%PEDGE(I,J,L70mb) 

          !==============================================================
          ! T2 = temperature (K)  at the sigma center of level L70mb
          ! T1 = temperature (K)  at the sigma center of level L70mb-1
          !
          ! DZ is the height from the sigma center of level L70mb-1 
          ! to 70mb.  Therefore, DZ may be found in either the 
          ! (L70mb)th sigma layer or the (L70mb-1)th sigma layer.  
          !
          ! ZUP is the height from the sigma center of the 
          ! (L70mb-1)th layer
          !============================================================== 

          ! Use virtual temperature for hypsometric equation (ewl, 3/3/15)
          T2   = State_Met%TV(I,J,L70mb  )
          T1   = State_Met%TV(I,J,L70mb-1)

          DZ   = Rdg0 * ( (T1 + T2) / 2e+0_fp ) * LOG( P1 / P70mb ) 
          ZUP  = Rdg0 * T1 * LOG( P1 /P3 )

          !==============================================================       
          ! H70mb is height between 70mb and the upper edge of the 
          ! level where DZ is.
          !  
          ! If DZ >= ZUP then DZ is already in level L70mb.
          ! If DZ <  ZUP then DZ is in level L70mb-1.
          !==============================================================       
          IF ( DZ >= ZUP ) THEN
             H70mb = State_Met%BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
          ELSE
             L70mb = L70mb - 1
             H70mb = ZUP - DZ
          ENDIF

          !=========================================================== 
          ! Distribute O3 into the region (30S-30N, 70mb-10mb)
          !=========================================================== 
          DO L = L70mb, LLPAR 

             ! Convert O3 in grid box (I,J,L) from v/v/s to v/v/box
             PO3 = PO3_vmr * DTCHEM 

             ! For both 2 x 2.5 and 4 x 5 GEOS grids, 30S and 30 N are
             ! grid box centers.  However, the O3 release region is 
             ! edged by 30S and 30N.  Therefore, if we are at the 30S
             ! or 30N grid boxes, divide the O3 flux by 2.
             IF ( J == J30S .or. J == J30N ) THEN
                PO3 = PO3 / 2e+0_fp
             ENDIF

             ! If we are in the lower level, compute the fraction
             ! of this level that lies above 70 mb, and scale 
             ! the O3 flux accordingly.
             IF ( L == L70mb ) THEN
                PO3 = PO3 * H70mb / State_Met%BXHEIGHT(I,J,L) 
             ENDIF

             ! Store O3 flux in the proper species number
             Spc(I,J,L,id_O3) = Spc(I,J,L,id_O3) + PO3 

             ! Store O3 flux for strat O3 species (Tagged O3 only)
             IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN
                Spc(I,J,L,id_O3Strat) = Spc(I,J,L,id_O3Strat) + PO3
             ENDIF

             ! Archive stratospheric O3 for printout in [Tg/yr]
             IF ( FIRST ) THEN
                STFLUX(I,J,L) = STFLUX(I,J,L) + &
                 PO3 * State_Met%AD(I,J,L) * 1000.e+0_fp / 28.8e+0_fp / &
                 DTCHEM * 48.e+0_fp * 365.25e+0_fp * 86400e+0_fp / 1e12
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

    !=================================================================
    ! Print amount of stratospheric O3 coming down
    !=================================================================
    IF ( FIRST ) THEN
       IF ( am_I_Root ) THEN
          WRITE( 6, 20 ) SUM( STFLUX )
       ENDIF
20     FORMAT( '     - Do_Synoz: Strat O3 production is', f9.3, ' [Tg/yr]' )
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE Do_Synoz
!EOC
END MODULE STRAT_CHEM_MOD
