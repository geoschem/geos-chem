!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: strat_chem_mod
!
! !DESCRIPTION: Module STRAT\_CHEM\_MOD contains variables and routines for 
!  performing a simple linearized chemistry scheme in the stratosphere, 
!  using archived 3D monthly climatological production rates and loss 
!  frequencies are applied from the GMI combo model.
!
!  In the original schem code (schem.F), only the following species
!  were destroyed by photolysis in the stratosphere:
!    PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O, N2O5, HNO4, MP
!  and by reaction with OH:
!    ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
!    PRPE, C3H8, CH2O, C2H6, HNO4, MP
!  
!  The updated code includes at least all of these, and many more. The code
!  is flexible enough to automatically apply the rate to any new tracers
!  for future simulations that share the name in tracer\_mod with the
!  GMI name.  (See Documentation on wiki).
!
!  The prod rates and loss frequencies are now read via HEMCO. They are 
!  stored in a data structure of flexible length (PLVEC). The file containing
!  the prod rates and loss frequencies need to be specified in the HEMCO 
!  configuration file for each species of interest. They are then automatically 
!  read and remapped onto the simulation grid. 
!  The field names assigned to the production and loss fields are expected to 
!  be 'GMI\_PROD\_XXX' and 'GMI\_LOSS\_XXX', respectively, where XXX is the 
!  species name. Production rates must be given in units of v/v/s, and loss 
!  frequencies in s-1.
!  The module variable PLMUSTFIND (set below) determines the behavior if no
!  production rates and/or loss frequencies can be found for any of the GMI
!  species defined in this module. IF PLMUSTFIND is set to TRUE, the code stops
!  with an error if no entry is found. Otherwise, stead-state values are used
!  for all species with no explicitly given values. 
!
!  The (monthly) OH concentrations are also obtained through HEMCO. The field
!  name must be 'STRAT\_OH', and values must be in v/v.
!\\
!\\
! !INTERFACE:
!
MODULE STRAT_CHEM_MOD
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
  PUBLIC  :: Calc_STE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Set_BryPointers
  PRIVATE :: Set_PLVEC
  PRIVATE :: Do_Synoz
!
! !PUBLIC DATA MEMBERS:
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Tracer index of Bry species in input files
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
  TYPE(BrPointers)      :: BrPtrDay(6) 
  TYPE(BrPointers)      :: BrPtrNight(6) 

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
  ! for every strat chem tracer. If set to TRUE, the code will return with
  ! an error if the production/loss rate cannot be found (through HEMCO) for
  ! any of the species. If set to .FALSE., only a warning is prompted and a
  ! value of 0.0 is used for every field that cannotbe found. 
  LOGICAL, PARAMETER   :: PLMUSTFIND = .FALSE.

  ! Number of species from GMI model
  INTEGER, PARAMETER   :: NTR_GMI   = 120  
!
! !PRIVATE TYPES:
!
  ! Scalars
  REAL(fp)              :: dTchem          ! chemistry time step [s]
  INTEGER               :: NSCHEM          ! Number of species upon which to 
                                           ! apply P's & k's in GEOS-Chem
  ! Arrays
  CHARACTER(LEN=16)     :: GMI_TrName(NTR_GMI)     ! Tracer names in GMI
  INTEGER               :: Strat_TrID_GC(NTR_GMI)  ! Maps 1:NSCHEM to STT index
  INTEGER               :: Strat_TrID_GMI(NTR_GMI) ! Maps 1:NSCHEM to GMI index
                     ! (At most NTR_GMI species could overlap between G-C & GMI)

  ! Tracer index of Bry species in GEOS-Chem STT (may differ from br_nos)
  INTEGER               :: GC_Bry_TrID(6) 

  ! Variables used to calculate the strat-trop exchange flux
  REAL(fp)              :: TauInit             ! Initial time
  INTEGER               :: NymdInit, NhmsInit  ! Initial date
  REAL(fp)              :: TpauseL_Cnt         ! Tropopause counter
  REAL(fp), ALLOCATABLE :: TpauseL(:,:)        ! Tropopause level aggregator
  REAL(f4), ALLOCATABLE :: MInit(:,:,:,:)      ! Init. atm. state for STE period
  REAL(f4), ALLOCATABLE :: SChem_Tend(:,:,:,:) ! Stratospheric chemical tendency
                                               !   (total P - L) [kg period-1]

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
! !IROUTINE: do_strat_chem
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
    USE CHEMGRID_MOD,       ONLY : GET_TPAUSE_LEVEL
    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
    USE CMN_SIZE_MOD
    USE DAO_MOD,            ONLY : CONVERT_UNITS
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE LINOZ_MOD,          ONLY : DO_LINOZ
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : TIMESTAMP_STRING
    USE TRACER_MOD,         ONLY : XNUMOLAIR
    USE TRACERID_MOD,       ONLY : IDTO3
    USE TRACERID_MOD,       ONLY : IDTCHBr3
    USE TRACERID_MOD,       ONLY : IDTCH2Br2
    USE TRACERID_MOD,       ONLY : IDTCH3Br

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd quantities
    LOGICAL, SAVE     :: FIRST      = .TRUE.
    INTEGER, SAVE     :: LASTMONTH  = -99

    ! Flags for simulation types
    LOGICAL           :: IT_IS_A_FULLCHEM_SIM
    LOGICAL           :: IT_IS_A_TAGOX_SIM
    LOGICAL           :: IT_IS_A_H2HD_SIM

    ! Scalars
    LOGICAL           :: prtDebug
    CHARACTER(LEN=16) :: STAMP
    INTEGER           :: I,    J,      L,   N,   NN
    REAL(fp)          :: dt,   P,      k,   M0,  RC,     M
    REAL(fp)          :: TK,   RDLOSS, T1L, mOH, BryTmp
    REAL(fp)          :: BOXVL
    LOGICAL           :: LLINOZ
    LOGICAL           :: LPRT
    LOGICAL           :: LBRGCCM
    LOGICAL           :: LRESET, LCYCLE
    LOGICAL           :: ISBR2 
    INTEGER           :: N_TRACERS

    ! Arrays
    REAL(fp)          :: STT0  (IIPAR,JJPAR,LLPAR,Input_Opt%N_TRACERS)
    REAL(fp)          :: BEFORE(IIPAR,JJPAR,LLPAR)
    REAL(fp)          :: TCVV(Input_Opt%N_TRACERS)

    ! Pointers
    REAL(fp), POINTER :: STT(:,:,:,:)
    REAL(fp), POINTER :: AD(:,:,:)
    REAL(fp), POINTER :: T(:,:,:)

    !===============================
    ! DO_STRAT_CHEM begins here!
    !===============================

    ! Assume Success
    errCode              = GIGC_SUCCESS

    ! Save values from the Input Options object to local variables
    N_TRACERS            = Input_Opt%N_TRACERS
    LLINOZ               = Input_Opt%LLINOZ
    LPRT                 = Input_Opt%LPRT
    LBRGCCM              = Input_Opt%LBRGCCM
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGOX_SIM    = Input_Opt%ITS_A_TAGOX_SIM  
    IT_IS_A_H2HD_SIM     = Input_Opt%ITS_A_H2HD_SIM
    TCVV                 = Input_Opt%TCVV(1:N_TRACERS)

    ! Initialize pointers
    STT                => State_Chm%Tracers
    AD                 => State_Met%AD
    T                  => State_Met%T

    ! Set a flag for debug printing
    prtDebug             = ( LPRT .and. am_I_Root )

    STAMP = TIMESTAMP_STRING()
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 ) STAMP
    ENDIF
10  FORMAT( '     - DO_STRAT_CHEM: Linearized strat chemistry at ', a )
    
    ! On first call, establish pointers to data fields read by HEMCO. These are
    ! the stratospheric Bry fields as well as the production/loss rates.
    ! (ckeller, 12/30/2014)
    IF ( FIRST ) THEN
       CALL Set_BryPointers ( am_I_Root, Input_Opt, State_Chm, State_Met, errCode )
       IF ( errCode /= GIGC_SUCCESS ) RETURN

       CALL Set_PLVEC ( am_I_Root, Input_Opt, State_Chm, State_Met, errCode )
       IF ( errCode /= GIGC_SUCCESS ) RETURN
    ENDIF

    ! SDE 2014-01-14: Allow the user to overwrite stratospheric
    ! concentrations at model initialization if necessary
    LRESET = (FIRST.AND.LBRGCCM)

    ! Set first-time flag to false
    FIRST = .FALSE.    

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### STRAT_CHEM: at DO_STRAT_CHEM' )
    ENDIF

    !======================>==========================================
    ! Full chemistry simulations
    !================================================================
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       ! Advance counter for number of times we've sampled the tropopause level
       TpauseL_CNT = TpauseL_CNT + 1e+0_fp

       !=============================================================
       ! Do chemical production and loss for non-ozone species for
       ! which we have explicit prod/loss rates from GMI
       !=============================================================

       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, N, NN, k, P, dt, M0 )
       DO J=1,JJPAR
          DO I=1,IIPAR

             ! Add to tropopause level aggregator for later determining STE flux
             TpauseL(I,J) = TpauseL(I,J) + GET_TPAUSE_LEVEL( I, J, State_Met )

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR

                IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) CYCLE

                DO N=1,NSCHEM ! Tracer index of active strat chem species
                   NN = Strat_TrID_GC(N) ! Tracer index in STT

                   ! Skip O3; we'll always use either Linoz or Synoz
                   IF ( IT_IS_A_FULLCHEM_SIM .and. NN .eq. IDTO3 ) CYCLE

                   dt = DTCHEM                              ! timestep [s]

                   ! original code:
!                   k = LOSS(I,J,L,N)                        ! loss freq [s-1]
!                   P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(NN) ! prod term [kg s-1]

                   ! loss freq [s-1] 
                   IF ( .NOT. ASSOCIATED(PLVEC(N)%LOSS) ) THEN
                      k = 0.0_fp
                   ELSE
                      k = PLVEC(N)%LOSS(I,J,L)
                   ENDIF

                   ! prod term [v/v/s --> kg/s]
                   IF ( .NOT. ASSOCIATED(PLVEC(N)%PROD) ) THEN
                      P = 0.0_fp 
                   ELSE
                      P = PLVEC(N)%PROD(I,J,L) * AD(I,J,L) / TCVV(NN) 
                   ENDIF

                   M0 = STT(I,J,L,NN)                       ! initial mass [kg]

                   ! No prod or loss at all
                   IF ( k .eq. 0e+0_fp .and. P .eq. 0e+0_fp ) CYCLE

                   ! Simple analytic solution to dM/dt = P - kM over [0,t]
                   IF ( k .gt. 0e+0_fp ) then
                      STT(I,J,L,NN) = M0 * EXP(-k*dt) + (P/k)*(1e+0_fp-EXP(-k*dt))
                   ELSE
                      STT(I,J,L,NN) = M0 + P*dt
                   ENDIF

                   ! Aggregate chemical tendency [kg box-1]
                   SCHEM_TEND(I,J,L,NN) = SCHEM_TEND(I,J,L,NN) + &
                                                       ( STT(I,J,L,NN) - M0 )

                ENDDO ! N
             ENDDO ! L
          ENDDO ! I
       ENDDO ! J
       !$OMP END PARALLEL DO

       !===================================
       ! Ozone
       !===================================

       ! Make note of inital state for determining tendency later
       BEFORE = STT(:,:,:,IDTO3 )

       ! Put ozone in v/v
       STT(:,:,:,IDTO3) = STT(:,:,:,IDTO3) * TCVV( IDTO3 ) / AD

       ! Do Linoz or Synoz
       IF ( LLINOZ ) THEN
          CALL Do_Linoz( am_I_Root, Input_Opt,             &
                         State_Met, State_Chm, RC=errCode )
       ELSE
          CALL Do_Synoz( am_I_Root, Input_Opt,             &
                         State_Met, State_Chm, RC=errCode )
       ENDIF

       ! Put ozone back to kg
       STT(:,:,:,IDTO3) = STT(:,:,:,IDTO3) * AD / TCVV( IDTO3 )

       ! Put tendency into diagnostic array [kg box-1]
       SCHEM_TEND(:,:,:,IDTO3) = SCHEM_TEND(:,:,:,IDTO3) + &
                                                  ( STT(:,:,:,IDTO3) - BEFORE )


       !========================================
       ! Reactions with OH
       ! Currently:
       !   (1) CHBr3  
       !   (2) CH2Br2  
       !   (3) CH3Br
       !========================================

       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, M, TK, RC, RDLOSS, T1L, mOH, BOXVL )
       DO J=1,JJPAR
          DO I=1,IIPAR  

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR

                IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) CYCLE

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
                IF ( IDTCH3Br .gt. 0 ) THEN
                   RC = 2.35e-12_fp * EXP ( - 1300.e+0_fp / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = STT(I,J,L,IDTCH3Br) * RDLOSS
                   STT(I,J,L,IDTCH3Br) = STT(I,J,L,IDTCH3Br) - T1L
                   SCHEM_TEND(I,J,L,IDTCH3Br) = &
                     SCHEM_TEND(I,J,L,IDTCH3Br) - T1L
                ENDIF

                !============!
                ! CHBr3 + OH !
                !============!
                IF ( IDTCHBr3 .gt. 0 ) THEN
                   RC = 1.35e-12_fp * EXP ( - 600.e+0_fp / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = STT(I,J,L,IDTCHBr3) * RDLOSS
                   STT(I,J,L,IDTCHBr3) = STT(I,J,L,IDTCHBr3) - T1L
                   SCHEM_TEND(I,J,L,IDTCHBr3) = &
                     SCHEM_TEND(I,J,L,IDTCHBr3) - T1L
                ENDIF

                !=============!
                ! CH2Br2 + OH !
                !=============!
                IF ( IDTCH2Br2 .gt. 0 ) THEN
                   RC = 2.00e-12_fp * EXP ( -  840.e+0_fp / TK )
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1e+0_fp )
                   T1L    = STT(I,J,L,IDTCH2Br2) * RDLOSS
                   STT(I,J,L,IDTCH2Br2) = STT(I,J,L,IDTCH2Br2) - T1L
                   SCHEM_TEND(I,J,L,IDTCH2Br2) = &
                     SCHEM_TEND(I,J,L,IDTCH2Br2) - T1L
                ENDIF

             ENDDO ! J
          ENDDO ! I
       ENDDO ! L

       !$OMP END PARALLEL DO

       !===============================
       ! Prescribe Br_y concentrations
       !===============================

       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( NN, BEFORE, I, J, L, BryTmp ) &
       !$OMP PRIVATE( LCYCLE )
       DO NN=1,6

          IF ( GC_Bry_TrID(NN) > 0 ) THEN

             ! Make note of inital state for determining tendency later
             ! NOTE: BEFORE has to be made PRIVATE to the DO loop since
             ! it only has IJL scope, but the loop is over IJLN!
             ! (bmy, 8/7/12)
             BEFORE = STT(:,:,:,GC_Bry_TrID(NN))

             ! Is this Br2?
             ISBR2 = ( TRIM(Input_Opt%TRACER_NAME(Strat_TrID_GC(NN))) == 'Br2' )

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR
             DO J = 1, JJPAR
             DO I = 1, IIPAR  
                 
                IF ( LRESET ) THEN
                   LCYCLE = ITS_IN_THE_TROP( I, J, L, State_Met )
                ELSE 
                   LCYCLE = ITS_IN_THE_CHEMGRID( I, J, L, State_Met )
                ENDIF
                IF ( LCYCLE ) CYCLE

                ! Now get Br data through HEMCO pointers (ckeller, 12/30/14).
                IF ( State_Met%SUNCOS(I,J) > 0.e+0_fp ) THEN
                   ! daytime [ppt] -> [kg]
                   BryTmp = BrPtrDay(NN)%MR(I,J,L)   &
                          * 1.e-12_fp                & ! convert from [ppt]
                          * AD(I,J,L)                &
                          / TCVV(GC_Bry_TrID(NN))

                ELSE
                   ! nighttime [ppt] -> [kg]
                   BryTmp = BrPtrNight(NN)%MR(I,J,L) &
                          * 1.e-12_fp                & ! convert from [ppt]
                          * AD(I,J,L)                &
                          / TCVV(GC_Bry_TrID(NN))
                ENDIF

                ! Special adjustment for G-C Br2 tracer, 
                ! which is BrCl in the strat (ckeller, 1/2/15)
                IF ( ISBR2 ) BryTmp = BryTmp / 2.0_fp

                ! Pass to STT array
                STT(I,J,L, GC_Bry_TrID(NN) ) = BryTmp

             ENDDO
             ENDDO
             ENDDO

             ! Put tendency into diagnostic array [kg box-1]
             SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) = &
                SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) + &
                ( STT(:,:,:,GC_Bry_TrID(NN)) - BEFORE )
          
          ENDIF

       ENDDO ! NN
       !$OMP END PARALLEL DO

    !======================================================================
    ! Tagged Ox simulation
    !======================================================================
    ELSE IF ( IT_IS_A_TAGOX_SIM ) THEN

       ! Tagged Ox only makes use of Synoz or Linoz. We apply either to
       ! the total Ox tracer, and the stratospheric Ox tracer.

       ! Intial conditions
       STT0(:,:,:,:) = STT(:,:,:,:)

       CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v

       IF ( LLINOZ ) THEN
          CALL Do_Linoz( am_I_Root, Input_Opt,             &
                         State_Met, State_Chm, RC=errCode )
       ELSE 
          CALL Do_Synoz( am_I_Root, Input_Opt,             &
                         State_Met, State_Chm, RC=errCode )
       ENDIF

       CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg

       ! Add to tropopause level aggregator for later determining STE flux
       TpauseL_CNT = TpauseL_CNT + 1e+0_fp

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          TpauseL(I,J) = TpauseL(I,J) + GET_TPAUSE_LEVEL( I, J, State_Met )
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Aggregate chemical tendency [kg box-1]
       DO N=1,NSCHEM
          NN = Strat_TrID_GC(N)
          SCHEM_TEND(:,:,:,N) = SCHEM_TEND(:,:,:,N) + &
               ( STT(:,:,:,NN) - STT0(:,:,:,NN) )
       ENDDO

    ELSE

       ! The code will need to be modified for other tagged simulations 
       ! (e.g., CO). Simulations like CH4, CO2 with standard tracer names 
       ! should probably just work as is with the full chemistry code above, 
       ! but would need to be tested.
       IF ( am_I_Root ) THEN
          WRITE( 6, '(a)' ) 'Strat chemistry needs to be activated for ' // &
                            'your simulation type.'
          WRITE( 6, '(a)' ) 'Please see GeosCore/strat_chem_mod.F90' // &
                            'or disable in input.geos'
       ENDIF
       CALL GEOS_CHEM_STOP
       
    ENDIF

    ! Free pointer
    NULLIFY( STT )
    NULLIFY( AD  )
    NULLIFY( T   )

  END SUBROUTINE DO_STRAT_CHEM
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
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr 

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
    CHARACTER(LEN=255)  :: PREFIX, FIELDNAME
    CHARACTER(LEN=1023) :: MSG
    INTEGER             :: N

    !=================================================================
    ! Set_BryPointers begins here 
    !=================================================================

    ! Construct error message
    MSG = 'Cannot get pointer from HEMCO! Stratospheric Bry data ' // &
          'is expected to be listed in the HEMCO configuration '   // &
          'file. This error occured when trying to get field'

    ! Do for every Bry species
    DO N = 1,6

       ! Skip if tracer is not defined    
       IF ( GC_Bry_TrID(N) <= 0 ) CYCLE

       ! Get Bry name
       ThisName = Input_Opt%TRACER_NAME( GC_Bry_TrID(N) )

       ! Construct field name using Bry name
       PREFIX = 'GEOSCCM_'//TRIM(ThisName)

       ! Get pointer to this field. These are the mixing ratios (pptv).

       ! Day
       FIELDNAME = TRIM(PREFIX) // '_DAY'
       CALL HCO_GetPtr( am_I_Root, FIELDNAME, BrPtrDay(N)%MR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP ( TRIM(MSG)//' '//TRIM(FIELDNAME), &
                            'Set_BryPointers (start_chem_mod.F90)' )
       ENDIF

       ! Night
       FIELDNAME = TRIM(PREFIX) // '_NIGHT'
       CALL HCO_GetPtr( am_I_Root, FIELDNAME, BrPtrNight(N)%MR, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP ( TRIM(MSG)//' '//TRIM(FIELDNAME), &
                            'Set_BryPointers (start_chem_mod.F90)' )
       ENDIF

    ENDDO !N

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Set_BryPointers 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_PLVEC
!
! !DESCRIPTION: Subroutine SET\_PLVEC gets the production and loss terms of
! all strat chem tracers from HEMCO. The pointers only need to be established 
! once. Target data is automatically updated through HEMCO. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_PLVEC ( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr 

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
    CHARACTER(LEN=1023) :: MSG, ERR
    INTEGER             :: N, TRCID
    LOGICAL             :: FND

    !=================================================================
    ! Set_PLVEC begins here 
    !=================================================================

    ! Assume error until success
    RC = GIGC_FAILURE

    ! Construct error message
    ERR = 'Cannot get pointer from HEMCO! GMI prod/loss data is '  // &
          'expected to be listed in the HEMCO configuration file. '// &
          'This error occured when trying to get field'

    ! Do for every species 
    DO N = 1,NSCHEM

       ! Get GEOS-Chem tracer index
       TRCID = Strat_TrID_GC(N)

       ! Skip if tracer is not defined    
       IF ( TRCID <= 0 ) CYCLE

       ! Get species name
       ThisName = Input_Opt%TRACER_NAME( TRCID )

       ! ---------------------------------------------------------------
       ! Get pointers to fields
       ! ---------------------------------------------------------------

       ! Production rates [v/v/s]
       FIELDNAME = 'GMI_PROD_'//TRIM(ThisName)
       CALL HCO_GetPtr( am_I_Root, FIELDNAME, PLVEC(N)%PROD, RC, FOUND=FND )
       IF ( RC /= HCO_SUCCESS .OR. ( PLMUSTFIND .AND. .NOT. FND) ) THEN
          CALL ERROR_STOP ( TRIM(ERR)//' '//TRIM(FIELDNAME), &
                            'Set_PLVEC (start_chem_mod.F90)' )
          RETURN
       ENDIF
       IF ( .NOT. FND .AND. AM_I_ROOT ) THEN
          MSG = 'Cannot find archived GMI production rates for ' // &
                TRIM(ThisName) // ' - will use value of 0.0. '   // &
                'To use archived rates, add the following field '// &
                'to the HEMCO configuration file: '//TRIM(FIELDNAME)
          WRITE(6,*) TRIM(MSG)
       ENDIF

       ! Loss frequency [s-1]
       FIELDNAME = 'GMI_LOSS_'//TRIM(ThisName)
       CALL HCO_GetPtr( am_I_Root, FIELDNAME, PLVEC(N)%LOSS, RC, FOUND=FND )
       IF ( RC /= HCO_SUCCESS .OR. ( PLMUSTFIND .AND. .NOT. FND) ) THEN
          CALL ERROR_STOP ( TRIM(ERR)//' '//TRIM(FIELDNAME), &
                            'Set_PLVEC (start_chem_mod.F90)' )
          RETURN
       ENDIF
       IF ( .NOT. FND .AND. AM_I_ROOT ) THEN
          MSG = 'Cannot find archived GMI loss frequencies for ' // &
                TRIM(ThisName) // ' - will use value of 0.0. '   // &
                'To use archived rates, add the following field '// &
                'to the HEMCO configuration file: '//TRIM(FIELDNAME)
          WRITE(6,*) TRIM(MSG)
       ENDIF

    ENDDO !N

    ! Get pointer to STRAT_OH
    CALL HCO_GetPtr( am_I_Root, 'STRAT_OH', STRAT_OH, RC, FOUND=FND )
    IF ( RC /= HCO_SUCCESS .OR. .NOT. FND ) THEN
       MSG = 'Cannot find monthly archived strat. OH field '    // &
             '`STRAT_OH`. Please add a corresponding entry to ' // &
             'the HEMCO configuration file.'
       CALL ERROR_STOP ( TRIM(MSG), 'Set_PLVEC (strat_chem_mod.F90)' )
       RETURN
    ENDIF

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Set_PLVEC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_ste
!
! !DESCRIPTION: Subroutine CALC\_STE estimates what the stratosphere-to-
!               troposphere exchange flux must have been since the last time
!               it was reset
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_STE( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE TIME_MOD,   ONLY : GET_TAU, GET_NYMD, GET_NHMS, EXPAND_DATE

    USE CMN_SIZE_MOD

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  28 Apr 2012 - L. Murray   - Initial version
!  18 Jul 2012 - R. Yantosca - Make sure I is the innermost DO loop
!                              (wherever expedient)
!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  05 Oct 2012 - R. Yantosca - Bug fix for IFORT 12: extend the #if statement
!                              to avoid including code for nested-grid sims
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm, RC arguments
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Scalars
    CHARACTER(LEN=255) :: dateStart, dateEnd
    INTEGER            :: N,         I,      J,    L,      NN
    REAL(fp)           :: dStrat,    STE,    Tend, tauEnd, dt

    ! Arrays
    INTEGER            :: LTP(IIPAR,JJPAR      )
    REAL(fp)           :: M1 (IIPAR,JJPAR,LLPAR)
    REAL(fp)           :: M2 (IIPAR,JJPAR,LLPAR)

    ! For fields from Input_Opt
    INTEGER            :: N_TRACERS

    ! Pointers
    ! We need to define local arrays to hold corresponding values 
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER  :: STT(:,:,:,:)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! By simple mass balance, dStrat/dt = P - L - STE,
    ! where STE is the net stratosphere-to-troposphere mass exchange. 
    !
    ! Therefore, we estimate STE as
    !   STE = (P-L) - dStrat/dt
    !
    ! As the tropopause is dynamic, we use the mean tropopause level during
    ! the period for determining initial and end stratospheric masses. 
    ! (ltm, 04/28/2012)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Assume success
    RC = GIGC_SUCCESS

#if defined( NESTED_NA ) || defined( NESTED_CH ) || defined( NESTED_EU )
    ! This method only works for a global domain.
    ! It could be modified for nested domains if the total mass flux across the
    ! boundaries during the period is taken into account.
    RETURN
#else

    ! Copy values from Input_Opt
    N_TRACERS = Input_Opt%N_TRACERS

    ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
    ! (mpayer, 12/6/12)
    STT       => State_Chm%Tracers

    ! Determine mean tropopause level for the period
    !$OMP PARALLEL DO                               &
    !$OMP DEFAULT( SHARED )                         &
    !$OMP PRIVATE( I,  J  )
    DO J = 1,JJPAR
    DO I = 1,IIPAR
       LTP(I,J) = NINT( TPauseL(I,J) / TPauseL_Cnt )
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Period over which STE is being determined [a]
    tauEnd = GET_TAU() ! [h]
    dt = ( tauEnd - tauInit ) / 24e+0_fp / 365.25e+0_fp

    dateStart = 'YYYY-MM-DD hh:mm'
    CALL EXPAND_DATE(dateStart,NymdInit,NhmsInit)
    dateEnd = 'YYYY-MM-DD hh:mm'
    CALL EXPAND_DATE(dateEnd,GET_NYMD(),GET_NHMS())

    ! Print to output
    IF ( am_I_Root ) THEN
       WRITE( 6, * ) ''
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a)' ) '  Strat-Trop Exchange'
       WRITE( 6, '(a)' ) REPEAT( '-', 79 )
       WRITE( 6, '(a)' ) &
            '  Global stratosphere-to-troposphere fluxes estimated over'
       WRITE( 6, 100 ) TRIM(dateStart), TRIM(dateEnd)
       WRITE( 6, * ) ''
       WRITE( 6, 110 ) 'Species','[moles a-1]','* [g/mol]','= [Tg a-1]'
    ENDIF
100 FORMAT( 2x,a16,' to ',a16 )
110 FORMAT( 2x,a8,':',4x,a11  ,4x,a9  ,4x,  a11 )

    ! Loop through each species
    DO N=1,N_TRACERS

       ! Populate before (M1) and after (M2) state for the species [kg]
       M1 = MInit(:,:,:,N)             
       M2 =   STT(:,:,:,N)

       ! Zero out tropopshere and determine total change in the stratospheric
       ! burden of species N (dStrat) [kg]
       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I,  J  )
       DO J = 1, JJPAR  
       DO I = 1, IIPAR
          M2(I,J,1:LTP(I,J)) = 0e+0_fp
          M1(I,J,1:LTP(I,J)) = 0e+0_fp
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
       dStrat   = SUM(M2)-SUM(M1)

       ! The total chemical tendency (P-L) over the period for species N [kg]
       Tend   = SUM(Schem_tend(:,:,:,N))

       ! Calculate flux as STE = (P-L) - dStrat/dt
       STE = (Tend-dStrat)/dt ! [kg a-1]

       ! Print to standard output
       IF ( am_I_Root ) THEN
          WRITE(6,120) TRIM(Input_Opt%TRACER_NAME(N)),      &
               STE/Input_Opt%TRACER_MW_KG(N),               & ! mol/a-1
               Input_Opt%TRACER_MW_KG(N)*1e+3_fp,           & ! g/mol
               STE*1e-9_fp                                    ! Tg a-1
       ENDIF

    ENDDO

120 FORMAT( 2x,a8,':',4x,e11.3,4x,f9.1,4x,f11.4 )

    IF ( am_I_Root ) THEN
       WRITE( 6, * ) ''
       WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
       WRITE( 6, * ) ''
    ENDIF

    ! Reset variables for next STE period
    NymdInit             = GET_NYMD()
    NhmsInit             = GET_NHMS()
    TauInit              = GET_TAU()
    TPauseL_Cnt          = 0e+0_fp
    TPauseL(:,:)         = 0e+0_fp
    SChem_tend(:,:,:,:)  = 0e+0_fp
    MInit(:,:,:,:)       = STT(:,:,:,:)

    ! Free pointer
    NULLIFY( STT )
#endif
  END SUBROUTINE Calc_STE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_strat_chem
!
! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.  
!  It also opens the necessary rate files.
!\\
!\\
! !INTERFACE:
!      
  SUBROUTINE INIT_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ERROR_MOD,          ONLY : ALLOC_ERR
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE TRACERID_MOD,       ONLY : IDTCHBr3, IDTCH2Br2, IDTCH3Br
    USE TRACERID_MOD,       ONLY : IDTBr2,   IDTBr,     IDTBrO
    USE TRACERID_MOD,       ONLY : IDTHOBr,  IDTHBr,    IDTBrNO3
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=16) :: sname
    INTEGER           :: AS, N, NN
    LOGICAL           :: IT_IS_A_FULLCHEM_SIM
    LOGICAL           :: IT_IS_A_TAGOX_SIM
    LOGICAL           :: LLINOZ
    LOGICAL           :: LUCX
    INTEGER           :: N_TRACERS

    ! Arrays
    CHARACTER(LEN=14) :: TRACER_NAME(Input_Opt%N_TRACERS)

    ! Pointers
    ! We need to define local arrays to hold corresponding values 
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER :: STT(:,:,:,:)

    !=================================================================
    ! INIT_STRAT_CHEM begins here!
    !=================================================================

    ! Assume success
    RC                       = GIGC_SUCCESS

    ! Save fields from the Input_Opt object to local variables
    LLINOZ                   = Input_Opt%LLINOZ
    LUCX                     = Input_Opt%LUCX
    N_TRACERS                = Input_Opt%N_TRACERS
    IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGOX_SIM        = Input_Opt%ITS_A_TAGOX_SIM
    TRACER_NAME(1:N_TRACERS) = Input_Opt%TRACER_NAME(1:N_TRACERS)

    ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
    ! (mpayer, 12/6/12)
    STT => State_Chm%Tracers

    ! Initialize counters, initial times, mapping arrays
    TpauseL_Cnt              = 0.e+0_fp
    NSCHEM                   = 0
    TauInit                  = GET_TAU()
    NymdInit                 = GET_NYMD()
    NhmsInit                 = GET_NHMS()
    strat_trID_GC(:)         = 0
    strat_trID_GMI(:)        = 0

    ! Initialize timestep for chemistry [s]
    dTchem = GET_TS_CHEM() * 60e+0_fp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Determine the mapping for the GMI to the GC variables based on
    ! tracer name, which only needs to be done once per model run.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! List of available tracers with archived monthly climatological
    ! production rates, loss frequencies, and mixing ratios from the 
    ! GMI Combo model (tracer names here are as used in GMI).
    GMI_TrName = (/ &
             'A3O2',     'ACET',   'ACTA',   'ALD2',    'ALK4',  'ATO2', &
             'B3O2',       'Br',   'BrCl',    'BrO',  'BrONO2',  'C2H6', &
             'C3H8',     'CCl4', 'CF2Br2', 'CF2Cl2', 'CF2ClBr', 'CF3Br', &
           'CFC113',   'CFC114', 'CFC115',  'CFCl3',    'CH2O', 'CH3Br', &
          'CH3CCl3',    'CH3Cl',    'CH4',     'CO',      'Cl',   'Cl2', &
            'Cl2O2',      'ClO', 'ClONO2',    'EOH',    'ETO2',   'ETP', &
             'GCO3',     'GLYC',   'GLYX',     'GP',    'GPAN',     'H', &
               'H2',    'H2402',    'H2O',   'H2O2',     'HAC',   'HBr', &
         'HCFC141b', 'HCFC142b', 'HCFC22',  'HCOOH',     'HCl',  'HNO2', &
             'HNO3',     'HNO4',    'HO2',   'HOBr',    'HOCl',  'IALD', &
             'IAO2',      'IAP',   'INO2',   'INPN',    'ISN1',  'ISNP', &
             'ISOP',      'KO2',   'MACR',   'MAN2',    'MAO3',  'MAOP', &
              'MAP',     'MCO3',    'MEK',   'MGLY',     'MO2',   'MOH', &
               'MP',     'MRO2',    'MRP',    'MVK',    'MVN2',     'N', &
              'N2O',     'N2O5',     'NO',    'NO2',     'NO3',   'NOx', &
                'O',      'O1D',     'O3',   'OClO',      'OH',    'Ox', &
              'PAN',      'PMN',    'PO2',     'PP',     'PPN',  'PRN1', &
             'PRPE',     'PRPN',   'R4N1',   'R4N2',    'R4O2',   'R4P', &
             'RA3P',     'RB3P',   'RCHO',   'RCO3',   'RCOOH',  'RIO1', &
             'RIO2',      'RIP',    'ROH',     'RP',    'VRO2',   'VRP'    /)

    !===========================!
    ! Full chemistry simulation !
    !===========================!
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       DO NN = 1, NTR_GMI       

          sname = TRIM(GMI_TrName(NN))

          ! Some species names don't exactly match GEOS-Chem names
          !IF ( TRIM(GMI_TrName(NN)) .eq. 'BrONO2' ) sname = 'BrNO3'
          ! Need to line up CFCs correctly
          IF     ( TRIM(GMI_TrName(NN)) .eq. 'CF2Br2'  ) THEN
              sname = 'H1202'
          ELSEIF ( TRIM(GMI_TrName(NN)) .eq. 'CF2Cl2'  ) THEN
              sname = 'CFC12'
          ELSEIF ( TRIM(GMI_TrName(NN)) .eq. 'CF2ClBr' ) THEN
              sname = 'H1211'
          ELSEIF ( TRIM(GMI_TrName(NN)) .eq. 'CF3Br'   ) THEN
              sname = 'H1311'
          ELSEIF ( TRIM(GMI_TrName(NN)) .eq. 'CFCl3'   ) THEN
              sname = 'CFC11'
          ENDIF
 

          DO N = 1, N_TRACERS

             ! For now, guarantee that GMI prod/loss rates are not used for any
             ! bromine species
             IF ( TRIM(TRACER_NAME(N)) .eq.      'Br' .or. &
                  TRIM(TRACER_NAME(N)) .eq.    'BrCl' .or. &
                  TRIM(TRACER_NAME(N)) .eq.     'BrO' .or. &
                  TRIM(TRACER_NAME(N)) .eq.  'BrONO2' .or. &
                  TRIM(TRACER_NAME(N)) .eq.  'CF2Br2' .or. &
                  TRIM(TRACER_NAME(N)) .eq. 'CF2ClBr' .or. &
                  TRIM(TRACER_NAME(N)) .eq.   'CF3Br' .or. &
                  TRIM(TRACER_NAME(N)) .eq.   'CH3Br' .or. &
                  TRIM(TRACER_NAME(N)) .eq.     'HBr' .or. &
                  TRIM(TRACER_NAME(N)) .eq.    'HOBr'        ) CYCLE

             ! SDE 08/28/13: Full strat. has its own mesospheric NOy handling
             IF ( LUCX ) THEN
                IF ( TRIM(TRACER_NAME(N)) .eq.    'NO' .or. &
                     TRIM(TRACER_NAME(N)) .eq.   'NO2' .or. &
                     TRIM(TRACER_NAME(N)) .eq.   'NO3' .or. &
                     TRIM(TRACER_NAME(N)) .eq.   'NOx' .or. &
                     TRIM(TRACER_NAME(N)) .eq.     'N' .or. &
                     TRIM(TRACER_NAME(N)) .eq.   'N2O' ) CYCLE
             ENDIF

             IF ( TRIM(TRACER_NAME(N)) .eq. TRIM(sname) ) THEN
                
                IF ( LLINOZ .and. TRIM(TRACER_NAME(N)) .eq. 'O3' ) THEN
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM(TRACER_NAME(N)) // ' (via Linoz)'
                   ENDIF
                ELSE IF ( TRIM(TRACER_NAME(N)) .eq. 'O3' ) THEN
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM(TRACER_NAME(N)) // ' (via Synoz)'
                   ENDIF
                ELSE
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM(TRACER_NAME(N))//' (via GMI rates)'
                   ENDIF
                ENDIF

                NSCHEM                 = NSCHEM + 1
                Strat_TrID_GC(NSCHEM)  = N  ! Maps 1:NSCHEM to STT index
                Strat_TrID_GMI(NSCHEM) = NN ! Maps 1:NSCHEM to GMI_TrName index

             ENDIF

          ENDDO
       ENDDO

       ! These are the reactions with which we will use OH fields
       ! to determine stratospheric loss.
       IF( am_I_Root ) THEN
          IF ( IDTCHBr3  .gt. 0 ) WRITE(6,*) 'CHBr3 (from GMI OH)'
          IF ( IDTCH2Br2 .gt. 0 ) WRITE(6,*) 'CH2Br2 (from GMI OH)'
          IF ( IDTCH3Br  .gt. 0 ) WRITE(6,*) 'CH3Br (from GMI OH)'
       ENDIF

!       ! Allocate array to hold monthly mean OH mixing ratio
!       ALLOCATE( STRAT_OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
!       IF ( AS /=0 ) CALL ALLOC_ERR( 'STRAT_OH' )
!       STRAT_OH = 0e0

       !===========!
       ! Tagged Ox !
       !===========!
    ELSE IF ( IT_IS_A_TAGOX_SIM ) THEN
       IF ( LLINOZ ) THEN
          IF ( am_I_Root ) THEN
             WRITE(6,*) 'Linoz ozone performed on: '
          ENDIF
       ELSE          
          IF ( am_I_Root ) THEN
             WRITE(6,*) 'Synoz ozone performed on: '
          ENDIF
       ENDIF
       DO N = 1, N_TRACERS
          IF ( TRIM(TRACER_NAME(N)) .eq. 'O3'        .or. &
               TRIM(TRACER_NAME(N)) .eq. 'O3Strt'    .or. &
               TRIM(TRACER_NAME(N)) .eq. 'O3Strat' ) THEN
             NSCHEM = NSCHEM + 1
             Strat_TrID_GC(NSCHEM) = N
             IF ( am_I_Root ) THEN
                WRITE(6,*) TRIM(TRACER_NAME(N))
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    ! Allocate and initialize prod & loss arrays         !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!    ! Allocate PROD -- array for clim. production rates [v/v/s]
!    ALLOCATE( PROD( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
!    IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
!    PROD = 0e0
!
!    ! Allocate LOSS -- array for clim. loss freq [s-1]
!    ALLOCATE( LOSS( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
!    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS' )
!    LOSS = 0e0

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
    ! initial atm. conditions from restart file [kg].
    ALLOCATE( MInit( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MInit' )
    MInit = STT

    ! Array to determine the mean tropopause level over the period
    ! for which STE is being estimated.
    ALLOCATE( TPAUSEL( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPAUSEL' )
    TPAUSEL = 0e+0_fp

    ! Array to aggregate the stratospheric chemical tendency [kg period-1]
    ALLOCATE( SCHEM_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCHEM_TEND' )
    SCHEM_TEND = 0e0

    ! Allocate and initialize bromine arrays
    GC_Bry_TrID(1:6) = (/IDTBr2,IDTBr,IDTBrO,IDTHOBr,IDTHBr,IDTBrNO3/)

    ! Free pointer
    NULLIFY( STT )

  END SUBROUTINE INIT_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_strat_chem
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
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I

    IF ( ALLOCATED( MInit      ) ) DEALLOCATE( MInit      )
    IF ( ALLOCATED( TPAUSEL    ) ) DEALLOCATE( TPAUSEL    )
    IF ( ALLOCATED( SCHEM_TEND ) ) DEALLOCATE( SCHEM_TEND )

    ! Cleanup pointer vectors
    DO I = 1,6
       BrPtrDay(I)%MR   => NULL()
       BrPtrNight(I)%MR => NULL()
    ENDDO

    DO I = 1,NSCHEM
       PLVEC(I)%PROD => NULL()
       PLVEC(I)%LOSS => NULL()
    ENDDO
    IF ( ASSOCIATED(PLVEC) ) DEALLOCATE(PLVEC)

    ! Cleanup pointer to strat. OH
    STRAT_OH => NULL()

  END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_synoz
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
    USE CHEMGRID_MOD,       ONLY : GET_TPAUSE_LEVEL
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE TAGGED_Ox_MOD,      ONLY : ADD_STRAT_POX
    USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_YEAR
    USE TRACERID_MOD,       ONLY : IDTO3,       IDTO3Strt

    USE CMN_SIZE_MOD             ! Size parameters
    USE CMN_GCTM_MOD             ! Rdg0

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
!  05 Nov 2013 - R. Yantosca - Rename IDTOxStrt to IDTO3Strt
!  23 Jan 2014 - M. Sulprizio- Linoz does not call UPBDFLX_O3. Synoz does. 
!                              Now uncomment ADD_STRAT_POx (jtl,hyl,dbj,11/3/11)
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE and GET_PCENTER with
!                              State_Met%PEDGE and State_Met%PMID. Remove
!                              dependency on pressure_mod.
!  03 Mar 2015 - E. Lundgren - Use virtual temperature in hypsometric eqn
!  12 Aug 2015 - R. Yantosca - Add placeholder values for 0.5 x 0.625 grids
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE        :: FIRST = .TRUE.
    INTEGER              :: I, J, L, L70mb
    REAL(fp)               :: P1, P2, P3, T1, T2, DZ, ZUP
    REAL(fp)               :: DTCHEM, H70mb, PO3, PO3_vmr
    REAL(fp)               :: STFLUX(IIPAR,JJPAR,LLPAR)

    ! Select the grid boxes at the edges of the O3 release region, 
    ! for the proper model resolution (qli, bmy, 12/1/04)
#if defined( GRID4x5 ) && defined( GCAP )
    ! GCAP has 45 latitudes, shift by 1/2 grid box (swu, bmy, 5/25/05)
    INTEGER, PARAMETER   :: J30S = 16, J30N = 30 

#elif defined( GRID4x5 )
    INTEGER, PARAMETER   :: J30S = 16, J30N = 31 

#elif defined( GRID2x25 ) 
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#elif defined( GRID1x125 ) 
    INTEGER, PARAMETER   :: J30S = 61, J30N = 121

#elif defined( GRID05x0625 )

    !%%% ADD PLACEHOLDER VALUES, THESE AREN'T REALLY USED ANYMORE %%%
#if   defined( NESTED_CH )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 83
#elif   defined( NESTED_NA )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 41
#elif   defined( NESTED_EU )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 1  ! add later-checked . it is ok Anna Prot
#endif


#elif defined( GRID05x0666 )

! jtl, 10/26/11 
#if   defined( NESTED_CH )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 83
#elif   defined( NESTED_NA )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 41
#elif   defined( NESTED_EU )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = 1  ! add later-checked . it is ok Anna Prot
#endif

#elif defined( GRID025x03125 )

#if defined( NESTED_CH )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 161
#elif defined( NESTED_NA )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 161 !I think it should be 202/Anna Prot
!Anna Prot added 8 May 2015
#elif defined( NESTED_EU )
    INTEGER, PARAMETER   :: J30S = 1, J30N = 115
#endif

#elif defined( GRID1x1 ) 

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined(NESTED_EU)
    INTEGER, PARAMETER   :: J30S = 1,  J30N = JJPAR  ! 1x1 nested grids
#else  
    INTEGER, PARAMETER   :: J30S = 61, J30N = 121    ! 1x1 global grid
#endif

#elif defined( EXTERNAL_GRID )
    ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

    ! Lower pressure bound for O3 release (unit: mb)
    ! REAL(fp),  PARAMETER   :: P70mb = 70e+0_fp !PHS
    REAL(fp)  :: P70mb, PTP

    ! Pointers
    ! We need to define local arrays to hold corresponding values 
    ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
    REAL(fp), POINTER :: STT(:,:,:,:)

    !=================================================================
    ! Do_Synoz begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Chemical timestep [s]
    ! Originally, Synoz was in transport code, and used dynamic dT.
    DTCHEM = GET_TS_CHEM() * 60e+0_fp

    ! For O3 flux printout
    STFLUX = 0e+0_fp

    ! lower pressure !PHS
    P70mb = 70e+0_fp

    ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
    ! (mpayer, 12/6/12)
    STT => State_Chm%Tracers

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
#if   defined( GEOS_4 )
    PO3_vmr = 5.14e-14_fp                                 ! 3,3,7

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP ) || defined( MERRA2 )

    ! For now assume GEOS-5 has same PO3_vmr value 
    ! as GEOS-4; we can redefine later (bmy, 5/25/05)
    PO3_vmr = 5.14e-14_fp   

#elif defined( GCAP )

    ! For GCAP, assuming 3,3,7 (swu, bmy, 5/25/05)
    PO3_vmr = 5.0e-14_fp

#elif defined( GISS ) && defined( MODELE )

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

          !==============================================================
          ! L70mb is the 1st layer where pressure is equal to
          ! or smaller than 70 mb 
          !==============================================================

          !--------------------------------------------------------------
          ! Comment out for now (bmy, 10/2/07)
          ! replace L70mb with Tropopause pressure if the later is 
          ! lower -PHS #### still Beta testing ####
          !IF ( Input_Opt%LVARTROP ) THEN
          !   PTP = State_Met%TROPP(I,J)
          !   IF ( PTP < P70mb ) THEN
          !      P70mb = PTP
          !      !#### TESTING ####
          !      write(6,*)'#### RAISED bottom of O3 release region'
          !      write(6,*)'at ', i, j
          !      first=.true.
          !   ENDIF
          !ENDIF
          !--------------------------------------------------------------

          DO L = 1, LLPAR

             ! P2 = pressure [hPa] at the sigma center of level L70mb
             P2 = State_Met%PMID(I,J,L) 

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

#if   !defined( GCAP ) 
             ! For both 2 x 2.5 and 4 x 5 GEOS grids, 30S and 30 N are
             ! grid box centers.  However, the O3 release region is 
             ! edged by 30S and 30N.  Therefore, if we are at the 30S
             ! or 30N grid boxes, divide the O3 flux by 2.
             IF ( J == J30S .or. J == J30N ) THEN
                PO3 = PO3 / 2e+0_fp
             ENDIF
#endif
             ! If we are in the lower level, compute the fraction
             ! of this level that lies above 70 mb, and scale 
             ! the O3 flux accordingly.
             IF ( L == L70mb ) THEN
                PO3 = PO3 * H70mb / State_Met%BXHEIGHT(I,J,L) 
             ENDIF

             ! Store O3 flux in the proper tracer number
             STT(I,J,L,IDTO3) = STT(I,J,L,IDTO3) + PO3 

             ! Store O3 flux for strat Ox tracer (Tagged Ox only)
             ! UPBDFLX_O3 and thus ADD_STRAT_POX is called only 
             ! when Synoz is used (LLINOZ is FALSE) (jtl, hyl, dbj, 11/3/11)
             IF ( Input_Opt%ITS_A_TAGOX_SIM ) THEN
                CALL ADD_STRAT_POX( I, J, L, PO3, State_Chm )
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
    NULLIFY( STT )

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
