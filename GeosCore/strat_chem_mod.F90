!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!\\
!\\
! !INTERFACE:
!
MODULE STRAT_CHEM_MOD
!
! !USES:
!
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
  PRIVATE :: Get_Rates
  PRIVATE :: Get_Rates_Interp
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Tracer index of Bry species in input files
  ! 1:6 = (Br2, Br, BrO, HOBr, HBr, BrNO3) for both
  INTEGER, PARAMETER   :: br_nos(6) = (/44, 45, 46, 47, 48, 50/)

  ! Number of species from GMI model
  INTEGER, PARAMETER   :: NTR_GMI   = 120  
!
! !PRIVATE TYPES:
!
  ! Scalars
  REAL*8               :: dTchem          ! chemistry time step [s]
  INTEGER              :: NSCHEM          ! Number of species upon which to 
                                          ! apply P's & k's in GEOS-Chem
  ! Arrays
  REAL*8,  ALLOCATABLE, TARGET :: PROD(:,:,:,:)   ! Production rate [v/v/s]
  REAL*8,  ALLOCATABLE, TARGET :: LOSS(:,:,:,:)   ! Loss frequency [s-1]
  REAL*8,  ALLOCATABLE, TARGET :: STRAT_OH(:,:,:) ! Monthly mean OH [v/v]

  CHARACTER(LEN=16)    :: GMI_TrName(NTR_GMI)     ! Tracer names in GMI
  INTEGER              :: Strat_TrID_GC(NTR_GMI)  ! Maps 1:NSCHEM to STT index
  INTEGER              :: Strat_TrID_GMI(NTR_GMI) ! Maps 1:NSCHEM to GMI index
                     ! (At most NTR_GMI species could overlap between G-C & GMI)

  ! Variables for Br strat chemistry, moved here from SCHEM.f
  REAL*4, ALLOCATABLE  :: Bry_temp(:,:,:) 
  REAL*8, ALLOCATABLE  :: Bry_day(:,:,:,:) 
  REAL*8, ALLOCATABLE  :: Bry_night(:,:,:,:)

  ! Tracer index of Bry species in GEOS-Chem STT (may differ from br_nos)
  INTEGER              :: GC_Bry_TrID(6) 

  ! Variables used to calculate the strat-trop exchange flux
  REAL*8               :: TauInit             ! Initial time
  INTEGER              :: NymdInit, NhmsInit  ! Initial date
  REAL*8               :: TpauseL_Cnt         ! Tropopause counter
  REAL*8, ALLOCATABLE  :: TpauseL(:,:)        ! Tropopause level aggregator
  REAL*8, ALLOCATABLE  :: MInit(:,:,:,:)      ! Init. atm. state for STE period
  REAL*8, ALLOCATABLE  :: SChem_Tend(:,:,:,:) ! Stratospheric chemical tendency
                                              !   (total P - L) [kg period-1]

  !=================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
  !=================================================================
CONTAINS
#if defined( ESMF_ )
  !----------------------------------------------------------------------
  ! TEMPORARY VERSION OF DO_STRAT_CHEM -- only invoked with ESMF!
  !
  ! This temporary version of DO_STRAT_CHEM is only called when we are
  ! connecting to the GEOS-5 GCM via the ESMF interface.  This temporary
  ! version skips the normal GMI strat chem routine and calls LINOZ.
  !
  ! This is just a temporary situation, and is intended only for testing.
  ! It is simpler to bring in the LINOZ chemistry to the GEOS-5 GCM
  ! first.  The GMI strat chem package is more complicated due to the
  ! large # of netCDF files that it has to read.  Those files all have
  ! to be read in via ESMF/MAPL and passed as inputs to the GEOS-Chem
  ! Chemistry Component via the Import State. (bmy, 3/18/13)
  !----------------------------------------------------------------------
  SUBROUTINE Do_Strat_Chem( am_I_Root, Input_Opt,     &
                            State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState 
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Linoz_Mod,          ONLY : Do_Linoz
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?

    ! Assume succes
    RC = GIGC_SUCCESS

    ! Do Linoz or Synoz
    IF ( Input_Opt%LLINOZ ) THEN
       write(6,*) '### Shunting GMI strat chem, doing LINOZ instead'
       CALL Do_Linoz( am_I_Root, Input_Opt, State_Met, RC )
    ENDIF

  END SUBROUTINE DO_STRAT_CHEM
#else
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
    USE TRACER_MOD,         ONLY : STT
    USE TRACER_MOD,         ONLY : XNUMOLAIR
    USE TRACERID_MOD,       ONLY : IDTOX
    USE TRACERID_MOD,       ONLY : IDTCHBr3
    USE TRACERID_MOD,       ONLY : IDTCH2Br2
    USE TRACERID_MOD,       ONLY : IDTCH3Br
    USE TROPOPAUSE_MOD,     ONLY : GET_MIN_TPAUSE_LEVEL
    USE TROPOPAUSE_MOD,     ONLY : GET_TPAUSE_LEVEL
    USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

    IMPLICIT NONE

#include "define.h"
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
    INTEGER,        INTENT(OUT)   :: errCode     ! Success or failure
!f
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
!  18 Mar 2013 - R. Yantosca - Now pass Input_Opt via the arg list
!  19 Mar 2013 - R. Yantosca - Now only copy Input_Opt%TCVV(1:N_TRACERS)
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
    REAL*8            :: dt,   P,      k,   M0,  RC,     M
    REAL*8            :: TK,   RDLOSS, T1L, mOH, BryDay, BryNight
    LOGICAL           :: LLINOZ
    LOGICAL           :: LPRT
    INTEGER           :: N_TRACERS

    ! Arrays
    REAL*8            :: STT0  (IIPAR,JJPAR,LLPAR,Input_Opt%N_TRACERS)
    REAL*8            :: BEFORE(IIPAR,JJPAR,LLPAR)
    REAL*8            :: TCVV(Input_Opt%N_TRACERS)

    ! External functions
    REAL*8, EXTERNAL  :: BOXVL

    !===============================
    ! DO_STRAT_CHEM begins here!
    !===============================

    ! Assume Success
    errCode              = GIGC_SUCCESS

    ! Save values from the Input Options object to local variables
    N_TRACERS            = Input_Opt%N_TRACERS
    LLINOZ               = Input_Opt%LLINOZ
    LPRT                 = Input_Opt%LPRT
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGOX_SIM    = Input_Opt%ITS_A_TAGOX_SIM  
    IT_IS_A_H2HD_SIM     = Input_Opt%ITS_A_H2HD_SIM
    TCVV                 = Input_Opt%TCVV(1:N_TRACERS)

    ! Set a flag for debug printing
    prtDebug             = ( LPRT .and. am_I_Root )

    STAMP = TIMESTAMP_STRING()
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 ) STAMP
    ENDIF
10  FORMAT( '     - DO_STRAT_CHEM: Linearized strat chemistry at ', a )
    
    IF ( GET_MONTH() /= LASTMONTH ) THEN

       IF ( prtDebug ) THEN 
          CALL DEBUG_MSG( '### STRAT_CHEM: at GET_RATES' )
       ENDIF

       ! Read rates for this month
       IF ( IT_IS_A_FULLCHEM_SIM ) THEN
#if defined( GRID4x5 ) || defined( GRID2x25 )
          CALL GET_RATES( GET_MONTH(), Input_Opt, State_Chm,  &
                                       am_I_Root, errCode    )
#else
          ! For resolutions finer than 2x2.5, nested, 
          ! or otherwise exotic domains and resolutions
          CALL GET_RATES_INTERP( GET_MONTH(), am_I_Root )
#endif
       ENDIF

       ! Save month for next iteration
       LASTMONTH = GET_MONTH()
    ENDIF

    ! Set first-time flag to false
    FIRST = .FALSE.    

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### STRAT_CHEM: at DO_STRAT_CHEM' )
    ENDIF

    !================================================================
    ! Full chemistry simulations
    !================================================================
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       ! Advance counter for number of times we've sampled the tropopause level
       TpauseL_CNT = TpauseL_CNT + 1d0

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

                IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) CYCLE

                DO N=1,NSCHEM ! Tracer index of active strat chem species
                   NN = Strat_TrID_GC(N) ! Tracer index in STT

                   ! Skip Ox; we'll always use either Linoz or Synoz
                   IF ( IT_IS_A_FULLCHEM_SIM .and. NN .eq. IDTOx ) CYCLE

                   dt = DTCHEM                              ! timestep [s]
                   k = LOSS(I,J,L,N)                        ! loss freq [s-1]
                   P = PROD(I,J,L,N) * State_Met%AD(I,J,L) &! prod term [kg s-1]
                     / TCVV(NN)
                   M0 = STT(I,J,L,NN)                       ! initial mass [kg]

                   ! No prod or loss at all
                   IF ( k .eq. 0d0 .and. P .eq. 0d0 ) CYCLE

                   ! Simple analytic solution to dM/dt = P - kM over [0,t]
                   IF ( k .gt. 0d0 ) then
                      STT(I,J,L,NN) = M0 * EXP(-k*dt) + (P/k)*(1d0-EXP(-k*dt))
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
       BEFORE = STT(:,:,:,IDTOX )

       ! Put ozone in v/v
       STT(:,:,:,IDTOX ) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / &
                           State_Met%AD

       ! Do Linoz or Synoz
       IF ( LLINOZ ) THEN
          CALL Do_Linoz( am_I_Root, Input_Opt, State_Met, RC=errCode )
       ELSE
          CALL Do_Synoz( am_I_Root, State_Met )
       ENDIF

       ! Put ozone back to kg
       STT(:,:,:,IDTOX) = STT(:,:,:,IDTOX) * State_Met%AD / TCVV( IDTOX )

       ! Put tendency into diagnostic array [kg box-1]
       SCHEM_TEND(:,:,:,IDTOX) = SCHEM_TEND(:,:,:,IDTOX) + &
                                                  ( STT(:,:,:,IDTOX) - BEFORE )

       !========================================
       ! Reactions with OH
       ! Currently:
       !   (1) CHBr3  
       !   (2) CH2Br2  
       !   (3) CH3Br
       !========================================

       !$OMP PARALLEL DO &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, M, TK, RC, RDLOSS, T1L, mOH )
       DO J=1,JJPAR
          DO I=1,IIPAR  

             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR

                IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) CYCLE

                ! Density of air at grid box (I,J,L) in [molec cm-3]
                M = State_Met%AD(I,J,L) / BOXVL(I,J,L,State_Met) * XNUMOLAIR

                ! OH number density [molec cm-3]
                mOH = M * STRAT_OH(I,J,L)

                ! Temperature at grid box (I,J,L) in K
                TK = State_Met%T(I,J,L)

                !============!
                ! CH3Br + OH !
                !============!
                IF ( IDTCH3Br .gt. 0 ) THEN
                   RC = 2.35d-12 * EXP ( - 1300.d0 / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1d0 )
                   T1L    = STT(I,J,L,IDTCH3Br) * RDLOSS
                   STT(I,J,L,IDTCH3Br) = STT(I,J,L,IDTCH3Br) - T1L
                   SCHEM_TEND(I,J,L,IDTCH3Br) = &
                     SCHEM_TEND(I,J,L,IDTCH3Br) - T1L
                ENDIF

                !============!
                ! CHBr3 + OH !
                !============!
                IF ( IDTCHBr3 .gt. 0 ) THEN
                   RC = 1.35d-12 * EXP ( - 600.d0 / TK ) 
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1d0 )
                   T1L    = STT(I,J,L,IDTCHBr3) * RDLOSS
                   STT(I,J,L,IDTCHBr3) = STT(I,J,L,IDTCHBr3) - T1L
                   SCHEM_TEND(I,J,L,IDTCHBr3) = &
                     SCHEM_TEND(I,J,L,IDTCHBr3) - T1L
                ENDIF

                !=============!
                ! CH2Br2 + OH !
                !=============!
                IF ( IDTCH2Br2 .gt. 0 ) THEN
                   RC = 2.00d-12 * EXP ( -  840.d0 / TK )
                   RDLOSS = MIN( RC * mOH * DTCHEM, 1d0 )
                   T1L    = STT(I,J,L,IDTCH2Br2) * RDLOSS
                   STT(I,J,L,IDTCH2Br2) = STT(I,J,L,IDTCH2Br2) - T1L
                   SCHEM_TEND(I,J,L,IDTCHBr3) = &
                     SCHEM_TEND(I,J,L,IDTCHBr3) - T1L
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
       !$OMP PRIVATE( NN, BEFORE, I, J, L, BryDay, BryNight )
       DO NN=1,6

          IF ( GC_Bry_TrID(NN) > 0 ) THEN

             ! Make note of inital state for determining tendency later
             ! NOTE: BEFORE has to be made PRIVATE to the DO loop since
             ! it only has IJL scope, but the loop is over IJLN!
             ! (bmy, 8/7/12)
             BEFORE = STT(:,:,:,GC_Bry_TrID(NN))
          
             ! NOTE: For compatibility w/ the GEOS-5 GCM, we can no longer
             ! assume a minimum tropopause level.  Loop from 1,LLPAR instead.
             ! (bmy, 7/18/12)
             DO L = 1, LLPAR
             DO J = 1, JJPAR
             DO I = 1, IIPAR  
                  
                IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) CYCLE
                   
                IF ( State_Met%SUNCOS(I,J) > 0.d0 ) THEN
                   ! daytime [ppt] -> [kg]
                   BryDay = bry_day(I,J,L,NN)     &
                          * 1.d-12                & ! convert from [ppt]
                          * State_Met%AD(I,J,L)   &
                          / TCVV(GC_Bry_TrID(NN))
                   STT(I,J,L, GC_Bry_TrID(NN) ) = BryDay
                ELSE
                   ! nighttime [ppt] -> [kg]
                   BryNight = bry_night(I,J,L,NN)   &
                            * 1.d-12                & ! convert from [ppt]
                            * State_Met%AD(I,J,L)   &
                            / TCVV(GC_Bry_TrID(NN))
                   STT(I,J,L, GC_Bry_TrID(NN) ) = BryNight
                ENDIF
                   
             ENDDO
             ENDDO
             ENDDO

             ! Put tendency into diagnostic array [kg box-1]
             SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) = &
                SCHEM_TEND(:,:,:,GC_Bry_TrID(NN)) + &
                ( STT(:,:,:,GC_Bry_TrID(NN)) - BEFORE )
          
          ENDIF

       ENDDO
       !$OMP END PARALLEL DO
       
    !======================================================================
    ! Tagged Ox simulation
    !======================================================================
    ELSE IF ( IT_IS_A_TAGOX_SIM ) THEN

       ! Tagged Ox only makes use of Synoz or Linoz. We apply either to
       ! the total Ox tracer, and the stratospheric Ox tracer.

       ! Intial conditions
       STT0(:,:,:,:) = STT(:,:,:,:)

       CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT ) ! kg -> v/v
       IF ( LLINOZ ) THEN
          CALL Do_Linoz( am_I_Root, Input_Opt, State_Met, RC=errCode )
       ELSE 
          CALL Do_Synoz( am_I_Root, State_Met )
       ENDIF
       CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT ) ! v/v -> kg

       ! Add to tropopause level aggregator for later determining STE flux
       TpauseL_CNT = TpauseL_CNT + 1d0

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

    !======================================================================
    ! H2-HD Simulation
    !======================================================================
    ELSE IF ( IT_IS_A_H2HD_SIM ) THEN

       ! H2/HD uses upbdflx_H2, which is a modified Synoz.

       ! Intial conditions
       STT0(:,:,:,:) = STT(:,:,:,:)

       CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT ) ! kg -> v/v
       CALL UPBDFLX_HD( State_Met )
       CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT ) ! v/v -> kg

       ! Add to tropopause level aggregator for later determining STE flux
       TpauseL_CNT = TpauseL_CNT + 1d0
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

  END SUBROUTINE DO_STRAT_CHEM
#endif
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rates
!
! !DESCRIPTION: Function GET\_RATES reads from disk the chemical production
!  and loss rates for the species of interest
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_RATES( THISMONTH, Input_Opt, State_Chm, am_I_Root, RC )
!
! !USES:
!
    ! GEOS-Chem routines
    USE BPCH2_MOD,          ONLY : GET_NAME_EXT
    USE BPCH2_MOD,          ONLY : GET_RES_EXT
    USE BPCH2_MOD,          ONLY : GET_TAU0
    USE BPCH2_MOD,          ONLY : READ_BPCH2
    USE CMN_SIZE_MOD
    USE DIRECTORY_MOD,      ONLY : DATA_DIR
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TRANSFER_MOD,       ONLY : TRANSFER_3D

    ! netCDF routines
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_close

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: THISMONTH   ! Current month
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
!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  26 Oct 2012 - R. Yantosca - Now pass Chemistry State object for GIGC
!   9 Nov 2012 - R. Yantosca - Now pass Input Options object for GIGC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: FILENAME, DAYFILE, NIGHTFILE
    INTEGER            :: N,        M,       S
    INTEGER            :: F,        NN,      fileID
    REAL*8             :: XTAU
    INTEGER            :: N_TRACERS

    ! Arrays
    REAL*4             :: ARRAY ( IIPAR, JJPAR, LGLOB )  ! Full vertical res
    REAL*8             :: ARRAY2( IIPAR, JJPAR, LLPAR )  ! Actual vertical res
    CHARACTER(LEN=14)  :: TRACER_NAME(Input_Opt%N_TRACERS)

    !=================================================================
    ! GET_RATES begins here
    !=================================================================

    ! Intialize arrays
    LOSS = 0d0
    PROD = 0d0

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )

    !-----------------------------------------------------------------------
    !         %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
    !
    ! Do nothing, since we would not be reading in data here when using
    ! the ESMF interface.  The strat-chem data would have to come in
    ! via the ESMF Import State. (bmy, 11/7/12)
    !-----------------------------------------------------------------------

    ! Assume success
    RC                       = GIGC_SUCCESS

#else
    !----------------------------------------------------------------------
    !                %%%%% TRADITIONAL GEOS-Chem %%%%%
    !
    ! Current practice in the standard GEOS-Chem is to read strat chem
    ! prod/loss data from netCDF files. (bmy, 10/26/12)
    !----------------------------------------------------------------------

    ! Copy fields from Input_Opt to Llocal variables
    N_TRACERS                = Input_Opt%N_TRACERS
    TRACER_NAME(1:N_TRACERS) = Input_Opt%TRACER_NAME(1:N_TRACERS)

    IF ( am_I_Root ) THEN
       WRITE( 6, 11 ) &
          '       - Getting new strat prod/loss rates for month: ', THISMONTH
    ENDIF
11  FORMAT( a, I2.2 )

    M = THISMONTH

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Get stratospheric OH mixing ratio [v/v] 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILENAME = 'strat_chem_201206/gmi.clim.OH.' // GET_NAME_EXT() //  &
               '.'                              // GET_RES_EXT()  // '.nc'
    FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME )

    IF ( am_I_Root ) THEN
       WRITE( 6, 100 ) TRIM( filename )
100    FORMAT( '         => Reading from file: ', a )
    ENDIF

    call NcOp_Rd( fileID, TRIM( FILENAME ) )
    call NcRd( array, fileID, 'species',                     &
                              (/     1,     1,     1,  m /), & ! Start
                              (/ iipar, jjpar, lglob,  1 /)  ) ! Count
    call NcCl( fileID )

    ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
    call transfer_3D( array, array2 )

    STRAT_OH(:,:,:) = ARRAY2

    DO N=1,NSCHEM
       NN = Strat_TrID_GMI(N)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Open individual species file
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       FILENAME = 'strat_chem_201206/gmi.clim.' // &
                  TRIM( GMI_TrName(NN) ) // '.' // & 
                  GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
       FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME )

       IF ( am_I_Root ) THEN
          WRITE( 6, 100 ) TRIM( filename )
       ENDIF

       call NcOp_Rd( fileID, TRIM( FILENAME ) )

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Read production rate [v/v/s]
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       call NcRd( array, fileID, 'prod',                          &
                                 (/     1,     1,     1,  m  /),  & ! Start 
                                 (/ iipar, jjpar, lglob,  1  /)  )  ! Count

       ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
       call transfer_3D( array, array2 )

       PROD(:,:,:,N) = ARRAY2

       ! Special adjustment for G-C Br2 tracer, which is BrCl in the strat
       IF ( TRIM(TRACER_NAME(Strat_TrID_GC(N))) .eq. 'Br2' ) &
            PROD(:,:,:,N) = PROD(:,:,:,N) / 2d0

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Read loss frequencies [s^-1]
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       call NcRd( array, fileID, 'loss',                          &
                                 (/     1,     1,     1,  m  /),  & ! Start
                                 (/ iipar, jjpar, lglob,  1  /)  )  ! Count

       ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
       call transfer_3D( array, array2 )

       LOSS(:,:,:,N) = ARRAY2

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Close species file
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       call NcCl( fileID )

!#if defined( DEVEL ) 
! NOTE: Comment out the strat chem fields of State_Chem for now (bmy, 11/26/12)
!       !-----------------------------------------------------------------
!       !   %%%%% TESTING GIGC INTERFACE FROM EXISTING GEOS-CHEM %%%%%
!       !
!       ! Call the routine GIGC_Allocate_All (located in module file
!       ! GeosCore/gigc_environment_mod.F90) to allocate all lat/lon
!       ! allocatable arrays used by GEOS-Chem.  
!       !-----------------------------------------------------------------
!       State_Chm%Schm_P(:,:,:,N) = PROD(:,:,:,N)
!       State_Chm%Schm_k(:,:,:,N) = LOSS(:,:,:,N)
!#endif

    ENDDO

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Br_y species are handled differently
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !==============================================================
    ! Read in stored Bry species concentrations for stratosphere.
    ! Stored by Q. Liang using the GEOS CCM. (jpp, 6/27/2011)
    !==============================================================

    ! TAU value at the beginning of this month
    XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

    ! the daytime concentrations
    dayfile = TRIM( DATA_DIR ) // 'bromine_201205/' // &
         'CCM_stratosphere_Bry/Bry_Stratosphere_day.bpch.'// &
         GET_NAME_EXT()   // '.' // GET_RES_EXT()

    ! the nighttime concentrations
    nightfile = TRIM(DATA_DIR) // 'bromine_201205/' // &
         'CCM_stratosphere_Bry/Bry_Stratosphere_night.bpch.'// &
         GET_NAME_EXT()   // '.' // GET_RES_EXT()
    
    DO NN = 1, 6
       
       ! 1. Read daytime data
       CALL READ_BPCH2( DAYFILE, 'IJ-AVG-$', br_nos(NN), &   
            XTAU,    IGLOB,      JGLOB, &  
            LGLOB,   Bry_temp,   QUIET=.TRUE. )
       
       ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
       CALL TRANSFER_3D( Bry_temp(:,:,:), Bry_day(:,:,:,NN) )
       
       ! 2. Read nighttime data
       CALL READ_BPCH2( NIGHTFILE, 'IJ-AVG-$', br_nos(NN), &
            XTAU,      IGLOB,      JGLOB, &   
            LGLOB,     Bry_temp,   QUIET=.TRUE. )
       
       ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR) 
       CALL TRANSFER_3D( Bry_temp(:,:,:), Bry_night(:,:,:,NN) )
       
    ENDDO
#endif

  END SUBROUTINE GET_RATES
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rates_interp
!
!
! !DESCRIPTION: Function GET\_RATES\_INTERP reads from disk the chemical
! production and loss rates for the species of interest to resolutions finer
! than 2 x 2.5 (e.g., nested simluations) via simple nearest-neighbor mapping.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_RATES_INTERP( THISMONTH, am_I_Root )
!
! !USES:
!
    ! GEOS-Chem routines
    USE BPCH2_MOD,       ONLY : GET_NAME_EXT
    USE BPCH2_MOD,       ONLY : GET_RES_EXT
    USE BPCH2_MOD,       ONLY : GET_TAU0
    USE BPCH2_MOD,       ONLY : READ_BPCH2
    USE CMN_SIZE_MOD
    USE DIRECTORY_MOD,   ONLY : DATA_DIR_1x1
    USE GRID_MOD,        ONLY : GET_XMID
    USE GRID_MOD,        ONLY : GET_YMID
    USE LOGICAL_MOD,     ONLY : LLINOZ
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME
    USE TRANSFER_MOD,    ONLY : TRANSFER_3D
    USE TRANSFER_MOD,    ONLY : TRANSFER_3D_Bry

    ! netCDF routines
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_close

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: THISMONTH   ! Current month
    LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REVISION HISTORY: 
!  01 Feb 2011 - L. Murray   - Initial version
!  18 Jul 2012 - R. Yantosca - Make sure that I is the innermost DO loop
!                              (wherever expedient)
!  20 Jul 2012 - R. Yantosca - Now call routine TRANSFER_3D_Bry, which takes
!                              arrays of size (144,91,:) as input & output
!  20 Jul 2012 - R. Yantosca - Reorganized declarations for clarity
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: FILENAME, DAYFILE, NIGHTFILE
    INTEGER            :: N,        M,       S 
    INTEGER            :: F,        I,       J
    INTEGER            :: NN,       fileID
    REAL*8             :: XTAU

    ! Index arrays
    INTEGER            :: II(1)
    INTEGER            :: JJ(1)

    ! Arrays defined on the 2 x 2.5 grid
    REAL*4             :: XMID_COARSE ( 144                    )
    REAL*4             :: YMID_COARSE (        91              )
    REAL*4             :: BryTemp     ( 144,   91,    LGLOB    )
    REAL*8             :: BryDay2x25  ( 144,   91,    LLPAR, 6 )
    REAL*8             :: BryNight2x25( 144,   91,    LLPAR, 6 )

    ! Arrays defined on the nested grid
    ! "f2c" = fine to coarse mapping
    INTEGER            :: I_f2c       ( IIPAR                  )
    INTEGER            :: J_f2c       (        JJPAR           )
    REAL*4             :: COLUMN      (               LGLOB    )
    REAL*4             :: ARRAY       ( IIPAR, JJPAR, LGLOB    )
    REAL*8             :: ARRAY2      ( IIPAR, JJPAR, LLPAR    )

    ! Pointers
    REAL*8, POINTER    :: ptr_3D(:,:,:)

    !=================================================================
    ! GET_RATES_INTERP begins here
    !=================================================================

    ! Intialize arrays
    LOSS = 0d0
    PROD = 0d0

    ! Path to input data, use 2 x 2.5 file
    FILENAME = 'strat_chem_201206/gmi.clim.OH.' // GET_NAME_EXT() // '.2x25.nc'
    FILENAME = TRIM( DATA_DIR_1x1 ) // TRIM( FILENAME )

    ! Echo info
    IF ( am_I_Root ) THEN
       WRITE( 6, 11 ) THISMONTH
    ENDIF
11  FORMAT( '       - Getting new strat prod/loss rates for month: ', I2.2 )

    ! Open the netCDF file containing the rates
    IF ( am_I_Root ) THEN
       WRITE( 6, 12 ) TRIM( filename )
    ENDIF
12  FORMAT( '         => Interpolate to resolution from file: ', a )
    call Ncop_Rd( fileID, TRIM( filename ) )

    ! Get the lat and lon centers of the 2x2.5 GMI climatology
    ! WARNING MAKE 2x25 after testing
    call NcRd( XMID_COARSE, fileID, 'longitude', (/1/),  (/144/) )
    call NcRd( YMID_COARSE, fileID, 'latitude',  (/1/),  (/91/) )

    ! For each fine grid index, determine the closest coarse (2x2.5) index
    ! Note: This doesn't do anything special for the date line, and may 
    ! therefore not pick the exact closest if it is on the other side.
    ! Note: CMN_SIZE_MOD claims in its comments that IIPAR < IGLOB, but 
    ! in actuality, IIPAR = IGLOB and JJPAR = JGLOB, the dimensions of the 
    ! nested region.
    DO I=1,IGLOB
       II = MINLOC( ABS( GET_XMID(I,1,1) - XMID_COARSE ) )
       I_f2c(I) = II(1)
       !IF ( am_I_Root ) print*,'I:',I,'->',II(1)
    ENDDO
    DO J=1,JGLOB
       JJ = MINLOC( ABS( GET_YMID(1,J,1) - YMID_COARSE ) )
       J_f2c(J) = JJ(1)
       !IF ( am_I_Root ) print*,'J:',J,'->',JJ(1)
    ENDDO

    M = THISMONTH

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Get Stratospheric OH concentrations [v/v] 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DO J = 1, JGLOB
    DO I = 1, IGLOB

       call NcRd( column, fileID, 'species',           &
                  (/ I_f2c(I), J_f2c(J),     1, m /),  & ! Start
                  (/        1,        1, lglob, 1 /)  ) ! Count
       array( I, J, : ) = column

    ENDDO
    ENDDO
    call NcCl( fileID )
    ptr_3D => STRAT_OH
    call transfer_3D( array, ptr_3D )
    NULLIFY( ptr_3D )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Get Bry concentrations [ppt]
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !==============================================================
    ! Read in stored Bry species concentrations for stratosphere.
    ! Stored by Q. Liang using the GEOS CCM. (jpp, 6/27/2011)
    !==============================================================

    ! TAU value at the beginning of this month
    XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

    ! the daytime concentrations
    dayfile = TRIM( DATA_DIR_1x1 ) // 'bromine_201205/' // &
         'CCM_stratosphere_Bry/Bry_Stratosphere_day.bpch.'// &
         GET_NAME_EXT()   // '.2x25'
    
    ! the nighttime concentrations
    nightfile = TRIM(DATA_DIR_1x1 ) // 'bromine_201205/' // &
         'CCM_stratosphere_Bry/Bry_Stratosphere_night.bpch.'// &
         GET_NAME_EXT()   // '.2x25'
    
    DO NN = 1, 6
       
       ! 1. Read daytime data on the 2 x 2.5 grid
       CALL READ_BPCH2( DAYFILE, 'IJ-AVG-$', br_nos(NN), &   
                        XTAU,    144,        91, &  
                        LGLOB,   BryTemp,   QUIET=.TRUE. )
       
       ! Cast from REAL*4 to REAL*8 and resize to LLPAR
       CALL TRANSFER_3D_Bry( BryTemp, BryDay2x25(:,:,:,NN) )
       
       ! 2. Read nighttime data on the 2 x 2.5 grid
       CALL READ_BPCH2( NIGHTFILE, 'IJ-AVG-$', br_nos(NN), &
                        XTAU,    144,        91, &   
                        LGLOB,   BryTemp,   QUIET=.TRUE. )
       
       ! Cast from REAL*4 to REAL*8 and resize to LLPAR 
       CALL TRANSFER_3D_Bry( BryTemp, BryNight2x25(:,:,:,NN) )
       
    ENDDO

    ! Cast from global 2x2.5 to fine nested resolution
    DO J = 1, JGLOB
    DO I = 1, IGLOB
       Bry_day  ( I, J, :, : ) = BryDay2x25  ( I_f2c(I), J_f2c(J), :, : )
       Bry_night( I, J, :, : ) = BryNight2x25( I_f2c(I), J_f2c(J), :, : )
    ENDDO
    ENDDO

    DO N=1,NSCHEM
       NN = Strat_TrID_GMI(N)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Open individual species file
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! Path to input data, use 2 x 2.5 file
       FILENAME = 'strat_chem_201206/gmi.clim.' // &
            TRIM( GMI_TrName(NN) ) // '.' // GET_NAME_EXT() // '.2x25.nc'
       FILENAME = TRIM( DATA_DIR_1x1 ) // TRIM( FILENAME )

       IF ( am_I_Root ) THEN
          WRITE( 6, 12 ) TRIM( filename )
       ENDIF

       call NcOp_Rd( fileID, TRIM( FILENAME ) )

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Read production rate [v/v/s]
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       array = 0.0
       
       DO J = 1, JGLOB
       DO I = 1, IGLOB

          call NcRd( column, fileID, 'prod',                       &
                             (/ I_f2c(I), J_f2c(J),     1,  m /),  & ! Start
                             (/        1,        1, lglob,  1 /)  )  ! Count
          array( I, J, : ) = column

       ENDDO
       ENDDO

       ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR if necessary
       ptr_3D => PROD(:,:,:,N)
       call transfer_3D( array, ptr_3D )
       NULLIFY( ptr_3D )

       ! Special adjustment for Br2 tracer, which is BrCl in the strat
       IF ( TRIM(TRACER_NAME(Strat_TrID_GC(N))) .eq. 'Br2' ) &
                                           PROD(:,:,:,N) = PROD(:,:,:,N) / 2d0

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Read loss frequencies [s^-1]
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       array = 0.0

       DO J = 1, JGLOB
       DO I = 1, IGLOB

          call NcRd( column, fileID, 'loss',                       &
                             (/ I_f2c(I), J_f2c(J),     1,  m /),  & ! Start
                             (/        1,        1, lglob,  1 /)  )  ! Count
          array( I, J, : ) = column

       ENDDO
       ENDDO

       ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR if necessary
       ptr_3d => LOSS(:,:,:,N)
       call transfer_3D( array, array2 )
       NULLIFY( ptr_3D )

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Close species file
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       call NcCl( fileID )

    ENDDO

  END SUBROUTINE GET_RATES_INTERP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Calc_STE( am_I_Root )
!
! !USES:
!
    USE TRACER_MOD, ONLY : STT, TRACER_MW_KG, N_TRACERS, TRACER_NAME
    USE TIME_MOD,   ONLY : GET_TAU, GET_NYMD, GET_NHMS, EXPAND_DATE

    USE CMN_SIZE_MOD

    IMPLICIT NONE

#include "define.h"
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Scalars
    CHARACTER(LEN=255) :: dateStart, dateEnd
    INTEGER            :: N,         I,      J,    L,      NN
    REAL*8             :: dStrat,    STE,    Tend, tauEnd, dt

    ! Arrays
    INTEGER            :: LTP(IIPAR,JJPAR      )
    REAL*8             :: M1 (IIPAR,JJPAR,LLPAR)
    REAL*8             :: M2 (IIPAR,JJPAR,LLPAR)

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

#if defined( NESTED_NA ) || defined( NESTED_CH ) || defined( NESTED_EU )
    ! This method only works for a global domain.
    ! It could be modified for nested domains if the total mass flux across the
    ! boundaries during the period is taken into account.
    RETURN
#else

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
    dt = ( tauEnd - tauInit ) / 24d0 / 365.25d0

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
          M2(I,J,1:LTP(I,J)) = 0d0
          M1(I,J,1:LTP(I,J)) = 0d0
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
          WRITE(6,120) TRIM(TRACER_NAME(N)),  &
               STE/TRACER_MW_KG(N),           & ! mol/a-1
               TRACER_MW_KG(N)*1d3,           & ! g/mol
               STE*1d-9                         ! Tg a-1
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
    TPauseL_Cnt          = 0d0
    TPauseL(:,:)         = 0d0
    SChem_tend(:,:,:,:)  = 0d0
    MInit(:,:,:,:)       = STT(:,:,:,:)

#endif
  END SUBROUTINE Calc_STE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
    USE TRACER_MOD,         ONLY : STT
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
    INTEGER           :: N_TRACERS

    ! Arrays
    CHARACTER(LEN=14) :: TRACER_NAME(Input_Opt%N_TRACERS)

    !=================================================================
    ! INIT_STRAT_CHEM begins here!
    !=================================================================

    ! Assume success
    RC                       = GIGC_SUCCESS

    ! Save fields from the Input_Opt object to local variables
    LLINOZ                   = Input_Opt%LLINOZ
    N_TRACERS                = Input_Opt%N_TRACERS
    IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_TAGOX_SIM        = Input_Opt%ITS_A_TAGOX_SIM
    TRACER_NAME(1:N_TRACERS) = Input_Opt%TRACER_NAME(1:N_TRACERS)

    ! Initialize counters, initial times, mapping arrays
    TpauseL_Cnt              = 0.
    NSCHEM                   = 0
    TauInit                  = GET_TAU()
    NymdInit                 = GET_NYMD()
    NhmsInit                 = GET_NHMS()
    strat_trID_GC(:)         = 0
    strat_trID_GMI(:)        = 0

    ! Initialize timestep for chemistry [s]
    dTchem = GET_TS_CHEM() * 60d0

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

             IF ( TRIM(TRACER_NAME(N)) .eq. TRIM(sname) ) THEN
                
                IF ( LLINOZ .and. TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
                   IF ( am_I_Root ) THEN
                      WRITE( 6, '(a)' ) TRIM(TRACER_NAME(N)) // ' (via Linoz)'
                   ENDIF
                ELSE IF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
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

! NOTE: Comment out strat-chem fields of State_Chm for now (bmy, 11/26/12)
!#if defined( DEVEL ) 
!                !---------------------------------------------------------
!                ! %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
!                !
!                !---------------------------------------------------------
!                State_Chm%Schm_Id(NSCHEM)   = Strat_TrID_GC(NSCHEM)
!                State_Chm%Schm_Name(NSCHEM) = TRIM( TRACER_NAME(N) )
!#endif

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

       ! Allocate array to hold monthly mean OH mixing ratio
       ALLOCATE( STRAT_OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
       IF ( AS /=0 ) CALL ALLOC_ERR( 'STRAT_OH' )
       STRAT_OH = 0d0

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
          IF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' .or. &
               TRIM(TRACER_NAME(N)) .eq. 'OxStrt' ) THEN
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

    ! Allocate PROD -- array for clim. production rates [v/v/s]
    ALLOCATE( PROD( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
    PROD = 0d0

    ! Allocate LOSS -- array for clim. loss freq [s-1]
    ALLOCATE( LOSS( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS' )
    LOSS = 0d0

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
    TPAUSEL = 0d0

    ! Array to aggregate the stratospheric chemical tendency [kg period-1]
    ALLOCATE( SCHEM_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCHEM_TEND' )
    SCHEM_TEND = 0d0

    ! Allocate and initialize bromine arrays
    GC_Bry_TrID(1:6) = (/IDTBr2,IDTBr,IDTBrO,IDTHOBr,IDTHBr,IDTBrNO3/)
    ALLOCATE( Bry_temp( IIPAR, JJPAR, LGLOB ) )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Bry_temp' )
    Bry_temp = 0.
    ALLOCATE( Bry_day( IIPAR, JJPAR, LLPAR, 6 ) )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Bry_day' )
    Bry_day = 0.
    ALLOCATE( Bry_night( IIPAR, JJPAR, LLPAR, 6 ) )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'Bry_night' )
    Bry_night = 0.

  END SUBROUTINE INIT_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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

    IF ( ALLOCATED( PROD       ) ) DEALLOCATE( PROD       )
    IF ( ALLOCATED( LOSS       ) ) DEALLOCATE( LOSS       )
    IF ( ALLOCATED( STRAT_OH   ) ) DEALLOCATE( STRAT_OH   )
    IF ( ALLOCATED( MInit      ) ) DEALLOCATE( MInit      )
    IF ( ALLOCATED( TPAUSEL    ) ) DEALLOCATE( TPAUSEL    )
    IF ( ALLOCATED( SCHEM_TEND ) ) DEALLOCATE( SCHEM_TEND )

  END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Do_Synoz( am_I_Root, State_Met )   
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE LOGICAL_MOD,        ONLY : LVARTROP 
    USE PRESSURE_MOD,       ONLY : GET_PEDGE, GET_PCENTER
    USE TAGGED_OX_MOD,      ONLY : ADD_STRAT_POX
    USE TIME_MOD,           ONLY : GET_TS_CHEM, GET_YEAR
    USE TRACER_MOD,         ONLY : STT, ITS_A_TAGOX_SIM
    USE TRACERID_MOD,       ONLY : IDTOX, IDTOxStrt
    USE TROPOPAUSE_MOD,     ONLY : GET_TPAUSE_LEVEL

    USE CMN_SIZE_MOD             ! Size parameters
    USE CMN_GCTM_MOD             ! Rdg0

    IMPLICIT NONE
#include "define.h"
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE        :: FIRST = .TRUE.
    INTEGER              :: I, J, L, L70mb
    INTEGER              :: NTRACER, NTRACE2
    REAL*8               :: P1, P2, P3, T1, T2, DZ, ZUP
    REAL*8               :: DTCHEM, H70mb, PO3, PO3_vmr
    REAL*8               :: STFLUX(IIPAR,JJPAR,LLPAR)

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

#elif defined( GRID05x0666 )
    INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID025x03125 )
    INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID1x1 ) 

#if   defined( NESTED_CH ) || defined( NESTED_NA )
    INTEGER, PARAMETER   :: J30S = 1,  J30N = JJPAR  ! 1x1 nested grids
#else  
    INTEGER, PARAMETER   :: J30S = 61, J30N = 121    ! 1x1 global grid
#endif

#elif defined( EXTERNAL_GRID )
    ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

    ! Lower pressure bound for O3 release (unit: mb)
    ! REAL*8,  PARAMETER   :: P70mb = 70d0 !PHS
    REAL*8  :: P70mb, PTP

    !=================================================================
    ! Do_Synoz begins here!
    !=================================================================

    ! Chemical timestep [s]
    ! Originally, Synoz was in transport code, and used dynamic dT.
    DTCHEM = GET_TS_CHEM() * 60d0

    ! For O3 flux printout
    STFLUX = 0d0

    ! lower pressure !PHS
    P70mb = 70d0

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
    PO3_vmr = 5.14d-14                                 ! 3,3,7

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

    ! For now assume GEOS-5 has same PO3_vmr value 
    ! as GEOS-4; we can redefine later (bmy, 5/25/05)
    PO3_vmr = 5.14d-14   

#elif defined( GCAP )

    ! For GCAP, assuming 3,3,7 (swu, bmy, 5/25/05)
    PO3_vmr = 5.0d-14 

#elif defined( GISS ) && defined( MODELE )

    ! For Model E, assuming 3,3,7 and 475 Tg N a-1
    PO3_vmr = 4.84610e-14 !/ 2d0

    ! Scale as necessary for PreIndustrial and Paleo Climates
    ! These values determined by Linoz simulations run for each climate
    if ( GET_YEAR() .ge. 2102 .and. GET_YEAR() .le. 2105 ) then
       ! PIH was 3% higher
       PO3_vmr = PO3_vmr * 1.0273853d0
    endif
    if ( GET_YEAR() .ge. 2202 .and. GET_YEAR() .le. 2205 ) then
       ! LGM Webb STE was 7% higher
       PO3_vmr = PO3_vmr * 1.0723525d0
    endif
    if ( GET_YEAR() .ge. 2302 .and. GET_YEAR() .le. 2305 ) then
       ! LGM CLIMAP STE was 3% higher
       PO3_vmr = PO3_vmr * 1.0285232d0
    endif

#endif

    ! Store in the proper Ox tracer #
    NTRACER = IDTOX

    ! Only initialize on first time step
    IF ( FIRST ) STFLUX = 0d0

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
          !IF ( LVARTROP ) THEN
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
             P2 = GET_PCENTER(I,J,L) 

             IF ( P2 < P70mb ) THEN
                L70mb = L
                EXIT
             ENDIF
          ENDDO

          ! P1 = pressure [hPa] at the sigma center of level L70mb - 1   
          P1 = GET_PCENTER(I,J,L70mb-1) 

          ! P3 = pressure [hPa] at the lower sigma edge of level L70mb
          P3 = GET_PEDGE(I,J,L70mb) 

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
          T2   = State_Met%T(I,J,L70mb  )
          T1   = State_Met%T(I,J,L70mb-1)

          DZ   = Rdg0 * ( (T1 + T2) / 2d0 ) * LOG( P1 / P70mb ) 
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
                PO3 = PO3 / 2d0
             ENDIF
#endif
             ! If we are in the lower level, compute the fraction
             ! of this level that lies above 70 mb, and scale 
             ! the O3 flux accordingly.
             IF ( L == L70mb ) THEN
                PO3 = PO3 * H70mb / State_Met%BXHEIGHT(I,J,L) 
             ENDIF

             ! Store O3 flux in the proper tracer number
             STT(I,J,L,IDTOX) = STT(I,J,L,IDTOX) + PO3 

             ! Store O3 flux for strat Ox tracer (Tagged Ox only)
             IF ( ITS_A_TAGOX_SIM() ) THEN
                STT(I,J,L,IDTOxStrt) = STT(I,J,L,IDTOxStrt) + PO3
             ENDIF

             ! Archive stratospheric O3 for printout in [Tg/yr]
             IF ( FIRST ) THEN
                STFLUX(I,J,L) = STFLUX(I,J,L) + &
                     PO3 * State_Met%AD(I,J,L) * 1000.d0 / 28.8d0 / &
                     DTCHEM * 48.d0 * 365.25d0 * 86400d0 / 1e12
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

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
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: upbdflx_hd
!
! !DESCRIPTION: Subroutine UPBDFLX\_HD establishes the flux boundary 
! condition for HD coming down from the stratosphere. This is adapted from 
! the UPBDFLX\_O3 routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE UPBDFLX_HD( State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE PRESSURE_MOD,       ONLY : GET_PEDGE, GET_PCENTER
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TRACER_MOD,         ONLY : STT
    USE TRACERID_MOD,       ONLY : IDTHD, IDTH2
    USE GIGC_State_Met_Mod, ONLY : MetState
    
    USE CMN_SIZE_MOD             ! Size parameters
    USE CMN_GCTM_MOD             ! Rdg0
!
! !INPUT PARAMETERS:
!
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  Instead of calculating the fractionation of H2 in the stratosphere 
!  (where we would have to take into account fractionation of CH4),
!  we simply set the HD tracer concentrations in the stratosphere to
!  reproduce observed profiles in the UT/LS.
!                                                                         
!  References:
!  ============================================================================
!  (1) "Global Budget of Molecular Hydrogen and its Deuterium Content: 
!       Constraints from Ground Station, Cruise, and Aircraft Observations" 
!       Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!       submitted to J. Geophys. Res., 2007.
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - L. Jaegle, H. U. Price, P. Le Sager - Initial version
!  (1 ) First adapted from UPBDFLX_O3 (G-C v5-05-03) then merged w/ v7-04-12.
!        Added parallel DO loops. (phs, 9/18/07)
!  (26) Now set J30S and J30N for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (27) Remove support for COMPAQ compiler (bmy, 7/8/09)
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  10 Feb 2012 - R. Yantosca - Modified for 0.25 x 0.3125 grids
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Jun 2012 - L. Murray   - Moved from upbdflx_mod.F to here.
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: I, J, L, L70mb
    INTEGER              :: NTRACER
    REAL*8               :: P1, P2, P3, T1, T2, DZ, ZUP
    REAL*8               :: DTCHEM, H70mb, PO3_vmr!,PO3
    REAL*8               :: PHD, PHD_vmr, SCALE_HD!, HD_AVG

    ! Select the grid boxes at the edges of the HD release region, 
    ! for the proper model resolution 
#if   defined( GRID4x5 ) && defined( GCAP )
    ! GCAP has 45 latitudes, shift by 1/2 grid box (swu, bmy, 5/25/05)
    INTEGER, PARAMETER   :: J30S = 16, J30N = 30 

#elif defined( GRID4x5 )
    INTEGER, PARAMETER   :: J30S = 16, J30N = 31 

#elif defined( GRID2x25 ) 
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#elif defined( GRID1x125 ) 
    INTEGER, PARAMETER   :: J30S = 61, J30N = 121

#elif defined( GRID05x0666 )
    INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID025x03125 )
    INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID1x1 ) 

#if   defined( NESTED_CH ) || defined( NESTED_NA )
    INTEGER, PARAMETER   :: J30S = 1,  J30N = JJPAR  ! 1x1 nested grids
#else  
    INTEGER, PARAMETER   :: J30S = 61, J30N = 121    ! 1x1 global grid
#endif

#elif defined( EXTERNAL_GRID )
    ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
    INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

    ! Lower pressure bound for HD release (unit: mb)
    REAL*8,  PARAMETER   :: P70mb = 70d0

    !=================================================================
    ! UPBDFLX_HD begins here!
    !=================================================================

    ! Chemistry timestep [s]
    DTCHEM = GET_TS_CHEM() * 60d0

    !=================================================================
    ! For now the only HD release rates are for GEOS-3. This will
    ! likely need to be scaled with other met fields (GEOS-4,
    ! GEOS-5...) jaegle, 2/20/2007
    ! PO3_vmr is the release rate for Ozone. This is then scaled by
    ! SCALE_HD in order to obtain the HD profile, to obtain
    ! PHD_vmr [v/v/s]
    ! For now uses the GEOS-3 scale factor for all cases (phs)
    !=================================================================
#if   defined( GEOS_4 )

    PO3_vmr = 5.14d-14                                 ! 3,3,7

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

    ! For now assume GEOS-5 has same PO3_vmr value 
    ! as GEOS-4; we can redefine later (bmy, 5/25/05)
    PO3_vmr = 5.14d-14   

#elif defined( GCAP )

    ! For GCAP, assuming 3,3,7 (swu, bmy, 5/25/05)
    PO3_vmr = 5.0d-14 

#endif

    ! Define scaling factor for HD and scale PO3_vmr
    ! Standard:
    SCALE_HD = 4.0d-5
    PHD_vmr= PO3_vmr * SCALE_HD

    !=================================================================
    ! Select the proper tracer number to store HD into
    !=================================================================
    NTRACER = IDTHD

    !=================================================================
    ! Loop over latitude/longtitude locations (I,J)
    !=================================================================
    !$OMP PARALLEL DO &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,  J,  L,  P2,  L70mb, P1, P3 ) &
    !$OMP PRIVATE( T2, T1, DZ, ZUP, H70mb, PHD    ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = J30S, J30N 
       DO I = 1, IIPAR

          !===========================================================
          ! L70mb is the 1st layer where pressure is equal to
          ! or smaller than 70 mb 
          !
          ! P1 = pressure [ mb ] at the sigma center     of level L70mb - 1
          ! P3 = pressure [ mb ] at the lower sigma edge of level L70mb
          ! P2 = pressure [ mb ] at the sigma center     of level L70mb
          !===========================================================
          DO L = 1, LLPAR
             P2 = GET_PCENTER(I,J,L) 

             IF ( P2 < P70mb ) THEN
                L70mb = L
                EXIT
             ENDIF
          ENDDO

          P1 = GET_PCENTER(I,J,L70mb-1) 
          P3 = GET_PEDGE(I,J,L70mb) 

          !===========================================================
          ! T2 = temperature (K)  at the sigma center of level L70mb
          ! T1 = temperature (K)  at the sigma center of level L70mb-1
          !
          ! DZ is the height from the sigma center of level L70mb-1 
          ! to 70mb.  Therefore, DZ may be found in either the 
          ! (L70mb)th sigma layer or the (L70mb-1)th sigma layer.  
          !
          ! ZUP is the height from the sigma center of the 
          ! (L70mb-1)th layer
          !=========================================================== 
          T2   = State_Met%T(I,J,L70mb  )
          T1   = State_Met%T(I,J,L70mb-1)

          DZ   = Rdg0 * ( (T1 + T2) / 2d0 ) * LOG( P1 / P70mb ) 
          ZUP  = Rdg0 * T1 * LOG( P1 /P3 )

          !===========================================================       
          ! H70mb is height between 70mb and the upper edge of the 
          ! level where DZ is.
          !  
          ! If DZ >= ZUP then DZ is already in level L70mb.
          ! If DZ <  ZUP then DZ is in level L70mb-1.
          !===========================================================       
          IF ( DZ >= ZUP ) THEN
             H70mb = State_Met%BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
          ELSE
             L70mb = L70mb - 1
             H70mb = ZUP - DZ
          ENDIF

          !=========================================================== 
          ! Distribute HD into the region (30S-30N, 70mb-10mb)
          !=========================================================== 
          DO L = L70mb, LLPAR 

             ! Convert HD in grid box (I,J,L) from v/v/s to kg/box
             PHD = PHD_vmr * DTCHEM 

#if   !defined( GCAP ) 
             ! For both 2 x 2.5 and 4 x 5 GEOS grids, 30S and 30 N are
             ! grid box centers.  However, the O3 release region is 
             ! edged by 30S and 30N.  Therefore, if we are at the 30S
             ! or 30N grid boxes, divide the O3 flux by 2.
             IF ( J == J30S .or. J == J30N ) THEN
                PHD = PHD / 2d0
             ENDIF
#endif
             ! If we are in the lower level, compute the fraction
             ! of this level that lies above 70 mb, and scale 
             ! the HD flux accordingly.
             IF ( L == L70mb ) THEN
                PHD = PHD * H70mb / State_Met%BXHEIGHT(I,J,L) 
             ENDIF

             ! Store HD flux in the proper tracer number
             STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER) + PHD

          ENDDO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE UPBDFLX_HD
!EOC
END MODULE STRAT_CHEM_MOD
