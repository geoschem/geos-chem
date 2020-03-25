#ifdef EXCHANGE
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: exchange_mod.F90
!
! !DESCRIPTION: Module EXCHANGE\_MOD contains variables and routines which
!  are used to exchange data between two or more runs (Global Domain and
!  Nested Domains), thus to combine them into TWO-WAY nested run.
!  (yanyy, 6/18/14)
!\\
!\\
! !INTERFACE:
!
MODULE EXCHANGE_MOD
!
! !USES:
!
  USE PRECISION_MOD
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS
!
  PUBLIC   :: EXCHANGE_GLOBAL_INIT
  PUBLIC   :: EXCHANGE_NESTED_INIT
  PUBLIC   :: EXCHANGE_GLOBAL_POST
  PUBLIC   :: EXCHANGE_NESTED_PRE
  PUBLIC   :: EXCHANGE_NESTED_POST
  PUBLIC   :: ITS_TIME_FOR_EXCHANGE
  PUBLIC   :: INIT_EXCHANGE
  PUBLIC   :: CLEANUP_EXCHANGE
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REMARKS:
!  Diagram:
!
!                 Global            Nested
!  Init--------------------------------------------------
!     |
!     |        wait for unlock      wait for unlock
!     |
!     \--------------------------------------------------
!
!  Loop----------------------------------------------
!     |    /---------
!     | phony_ex_global_pre
!     |    \---------
!     |
!     |        GEOSCHEM STUFF
!     |
!     |    /---------
!     |    |   dump for nested
!     |    |   done lock                  ---------------\
!     |    |                        wait for unlock      |
!     |    |                        read LBCs         ex_nested_pre
!     |    |                        (regridded down)        |
!     | ex_global_post              apply LBCs           |
!     |    |                              ---------------/
!     |    |
!     |    |                        GEOSCHEM STUFF
!     |    |
!     |    |                              ---------------\
!     |    |                        dump for global(regridded up)
!     |    |                                          ex_nested_post
!     |    |                        done lock            |
!     |    |   wait for unlock            ---------------/
!     |    |   read & apply feedback
!     |    \---------
!     \------------------------------------------------
!
!
! !REVISION HISTORY:
!  30 Mar 2014 - Y.Y. Yan - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  !=================================================================
  ! MODULE VARIABLES
  !=================================================================
  ! Tmp arrays for exchange. Maybe we only need one array, but the
  ! code would be more difficult to understand.
  ! TMP_WINDOW is for unit convert, etc.
  REAL(fp), ALLOCATABLE        :: TMP_WINDOW(:,:,:)
  !REAL(fp), ALLOCATABLE        :: TMP_WINDOW_CH(:,:,:)
  REAL(fp)        :: TMP_WINDOW_CH(37,35,47)
  !REAL(fp), ALLOCATABLE        :: TMP_WINDOW_NA(:,:,:)
  REAL(fp)        :: TMP_WINDOW_NA(41,31,47)
  REAL(fp), ALLOCATABLE        :: TMP_WINDOW_EU(:,:,:)

  ! Exchange grid size and offset in global domain.
  ! I{LAT/LONG}_EX is the size of the window region.
  ! I{LAT/LONG}_OFF_EX is the window region's LR corner grid in
  ! global domain.
  INTEGER                    :: J0_BC, I0_BC
  INTEGER                    :: J1_BC, I1_BC
  INTEGER                    :: J2_BC, I2_BC
  INTEGER                    :: IM_BC, JM_BC
  INTEGER                    :: I0_BC_CH,I0_BC_NA,I0_BC_EU
  INTEGER                    :: IM_BC_CH,IM_BC_NA,IM_BC_EU
  INTEGER                    :: J0_BC_CH,J0_BC_NA,J0_BC_EU
  INTEGER                    :: JM_BC_CH,JM_BC_NA,JM_BC_EU

  ! R{LAT/LONG}_OFF_EX is nested's start offset in the window region,
  !     not offset from the global region.
  ! R{LAT/LONG}_OFF_EX is nested's end point offset in the window
  !     region. Counted in reversed direction.
  ! RATIO_{LAT/LONG}_EX is the grid size ratio of global/nested
  REAL(fp)                     :: RATIO_LONG_EX, RATIO_LAT_EX
  ! the parameters of boundary when regrid up global domain
  REAL(fp)                     :: OFFSET_LONG, OFFSET_LAT
  LOGICAL                    :: IS_GLOBAL_DOMAIN
  LOGICAL                    :: IS_EXCHANGE
  CHARACTER*80               :: MY_NAME
  CHARACTER*80               :: NESTED_CH_NAME = ''
  CHARACTER*80               :: NESTED_NA_NAME = ''
  CHARACTER*80               :: NESTED_EU_NAME = ''

  ! EXCHANGE_DIR is added into the run dir for data exchange
  CHARACTER*256, PARAMETER   :: EXCHANGE_DIR = './exchange'
  CHARACTER*256              :: DUMP_CSPEC_FNAME = ''
  CHARACTER*256              :: DUMP_TRACER_FNAME = ''
  CHARACTER*256              :: LOAD_NA_CSPEC_FNAME = ''
  CHARACTER*256              :: LOAD_NA_TRACER_FNAME = ''
  CHARACTER*256              :: LOAD_CH_CSPEC_FNAME = ''
  CHARACTER*256              :: LOAD_CH_TRACER_FNAME = ''
  CHARACTER*256              :: LOAD_EU_CSPEC_FNAME = ''
  CHARACTER*256              :: LOAD_EU_TRACER_FNAME = ''

  ! For Regrid
  ! NEW_GRID = REGRID_MAT_L * OLD_GRID * REGRID_MAT_R
  REAL(fp), ALLOCATABLE        :: REGRID_MAT_L(:,:)
  REAL(fp), ALLOCATABLE        :: REGRID_MAT_R(:,:)

  ! For lock utility in the run dir
  CHARACTER*256, PARAMETER   :: LOCK_DIR = './lock'
  CHARACTER*256              :: UNLOCK_FNAME
  CHARACTER*256              :: DONE_FNAME
  CHARACTER*256              :: KILL_FNAME
  CHARACTER*256              :: END_FNAME

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
! !IROUTINE: exchange_global_init
!
! !DESCRIPTION: Global initialization for exchange
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXCHANGE_GLOBAL_INIT()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_GLOBAL_INIT"
    CALL DONE_LOCK_ME_UP()
    CALL WAIT_FOR_UNLOCK()

  END SUBROUTINE EXCHANGE_GLOBAL_INIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exchange_nested_init
!
! !DESCRIPTION: Nested initialization for exchange
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXCHANGE_NESTED_INIT()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_NESTED_INIT"
    CALL DONE_LOCK_ME_UP()
    CALL WAIT_FOR_UNLOCK()

  END SUBROUTINE EXCHANGE_NESTED_INIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exchange_global_post
!
! !DESCRIPTION: Carry out the process to communicate the nested simulated
!  results back to global model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXCHANGE_GLOBAL_POST( Input_Opt, State_Chm, State_Grid, &
                                   State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GC_Error
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_GLOBAL_POST"
    CALL DONE_LOCK_ME_UP()
    CALL WAIT_FOR_UNLOCK()

    IF ( IS_EXCHANGE ) THEN

       ! Convert species units from [kg/kg dry] to [kg] (ewl, 8/13/15)
       CALL ConvertSpc_KgKgDry_to_Kg( State_Met, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error('Unit conversion error', RC, &
                        'EXCHANGE_GLOBAL_POST in exchange_mod.F')
          RETURN
       ENDIF

       ! mol/cm3 -> mol
       CALL EX_CONVERT_UNITS_CSPEC( Input_Opt, State_Chm, State_Grid, &
                                    State_Met, RC, 1)

# if  defined ( EXCHANGE_4x5_CH ) || defined ( EXCHANGE_2x25_CH )
       CALL EX_READ_AND_APPLY_FEEDBACK( Input_Opt, State_Chm, State_Grid, &
                                        State_Met, RC, 1)
# endif

# if  defined ( EXCHANGE_4x5_NA ) || defined ( EXCHANGE_2x25_NA )
       CALL EX_READ_AND_APPLY_FEEDBACK( Input_Opt, State_Chm, State_Grid, &
                                        State_Met, RC, 2)
# endif

# if  defined ( EXCHANGE_4x5_EU ) || defined ( EXCHANGE_2x25_EU )
       CALL EX_READ_AND_APPLY_FEEDBACK( Input_Opt, State_Chm, State_Grid, &
                                        State_Met, RC, 3)
# endif

       ! Convert species units back to kg/kg (ewl, 8/13/15)
       CALL ConvertSpc_Kg_to_KgKgDry( State_Met, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error('Unit conversion error', RC, &
                        'EXCHANGE_GLOBAL_POST in exchange_mod.F')
          RETURN
       ENDIF

       ! mol -> mol/cm3
       CALL EX_CONVERT_UNITS_CSPEC( Input_Opt, State_Chm, State_Grid, &
                                    State_Met, RC, 2 )

    ENDIF

    IS_EXCHANGE = .TRUE.
  END SUBROUTINE EXCHANGE_GLOBAL_POST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exchange_nested_pre
!
! !DESCRIPTION: Before nested grid exchange
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXCHANGE_NESTED_PRE()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_NESTED_PRE"
    CALL WAIT_FOR_UNLOCK()

  END SUBROUTINE EXCHANGE_NESTED_PRE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exchange_nested_post
!
! !DESCRIPTION: After nested-grid exchange
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXCHANGE_NESTED_POST( Input_Opt, State_Chm, State_Grid, &
                                   State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GC_Error
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !RETURN VALUE:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_NESTED_POST"

    ! Convert species units from [kg/kg dry] to [kg] (ewl, 8/13/15)
    CALL ConvertSpc_KgKgDry_to_Kg( State_Met, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'EXCHANGE_NESTED_POST in exchange_mod.F')
       RETURN
    ENDIF

    ! mol/cm3->mol
    CALL EX_CONVERT_UNITS_CSPEC( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC, 1 )

    CALL EX_DUMP_FOR_GLOBAL( Input_Opt, State_Chm, State_Grid, &
                             State_Met, RC )

    ! Convert species units back to kg/kg (ewl, 8/13/15)
    CALL ConvertSpc_Kg_to_KgKgDry( State_Met, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'EXCHANGE_NESTED_POST in exchange_mod.F')
       RETURN
    ENDIF

    ! mol -> mol/cm3
    CALL EX_CONVERT_UNITS_CSPEC( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC, 2 )

    CALL DONE_LOCK_ME_UP()
    CALL WAIT_FOR_UNLOCK()

  END SUBROUTINE EXCHANGE_NESTED_POST
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_exchange
!
! !DESCRIPTION: Returns TRUE if it's time to do the 2-way exchange of data.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ITS_TIME_FOR_EXCHANGE()
!
! !USES:
!
    USE TIME_MOD, ONLY : ITS_TIME_FOR_EXCH
!
! !RETURN VALUE:
!
    LOGICAL :: ITS_TIME_FOR_EXCHANGE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ITS_TIME_FOR_EXCHANGE = ITS_TIME_FOR_EXCH()

  END FUNCTION ITS_TIME_FOR_EXCHANGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_exchange
!
! !DESCRIPTION: Initialization routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_EXCHANGE( Input_Opt, State_Grid )
!
! !USES:
!
    USE ERROR_MOD
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(GrdState), INTENT(IN) :: State_Grid
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: ALLOCATE_RES
    REAL(fp), ALLOCATABLE       :: REGRID_MAT_L_TRANS(:,:)

    IF ( Input_Opt%amIRoot ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@ INIT_EXCHANGE"
    ENDIF

    IF ( .not. State_Grid%Nested_Grid ) THEN

       IS_GLOBAL_DOMAIN = .TRUE.
       IS_EXCHANGE = .FALSE.
       MY_NAME = 'GLOBAL'

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! This needs to be updated for FlexGrid (mps, 4/13/19)
       NESTED_CH_NAME = 'NESTED_CH'
       NESTED_NA_NAME = 'NESTED_NA'
       NESTED_EU_NAME = 'NESTED_EU'
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( LOAD_CH_CSPEC_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.cspec.', TRIM(NESTED_CH_NAME), '.dump'
          WRITE( LOAD_CH_TRACER_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.tracer.', TRIM(NESTED_CH_NAME), '.dump'
          WRITE( LOAD_NA_CSPEC_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.cspec.', TRIM(NESTED_NA_NAME), '.dump'
          WRITE( LOAD_NA_TRACER_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.tracer.', TRIM(NESTED_NA_NAME), '.dump'
          WRITE( LOAD_EU_CSPEC_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.cspec.', TRIM(NESTED_EU_NAME), '.dump'
          WRITE( LOAD_EU_TRACER_FNAME, '(A, A, A, A)' ) &
                 TRIM(EXCHANGE_DIR), &
                 '/dump.tracer.', TRIM(NESTED_EU_NAME), '.dump'
       ENDIF
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ELSE

       IS_GLOBAL_DOMAIN = .FALSE.
       MY_NAME = 'NESTED'

    ENDIF

    IF ( Input_Opt%amIRoot ) THEN
       WRITE ( DUMP_CSPEC_FNAME, '(A, A, A, A)' ) TRIM(EXCHANGE_DIR), &
               '/dump.cspec.', TRIM(MY_NAME), '.dump'
       WRITE ( DUMP_TRACER_FNAME, '(A, A, A, A)' ) TRIM(EXCHANGE_DIR), &
               '/dump.tracer.', TRIM(MY_NAME), '.dump'

       WRITE ( UNLOCK_FNAME, '(A, A, A, A)' ) TRIM(LOCK_DIR), '/key.', &
                                              TRIM(MY_NAME), '.key'
       WRITE ( DONE_FNAME,   '(A, A, A, A)' ) TRIM(LOCK_DIR), '/done.', &
                                              TRIM(MY_NAME), '.done'
       WRITE ( KILL_FNAME,   '(A, A, A, A)' ) TRIM(LOCK_DIR), '/kill.', &
                                              TRIM(MY_NAME), '.kill'
       WRITE ( END_FNAME,    '(A, A, A, A)' ) TRIM(LOCK_DIR), '/end.', &
                                              TRIM(MY_NAME), '.end'
    ENDIF

    ! size information
    ! Eventually need to remove hardcoding for FlexGrid (mps, 4/13/19)
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN
       I0_BC_CH = 50
       IM_BC_CH = 17
       J0_BC_CH = 20
       JM_BC_CH = 17
       I0_BC_NA = 8
       IM_BC_NA = 21
       J0_BC_NA = 25
       JM_BC_NA = 16
       I0_BC_EU = 30
       IM_BC_EU = 17
       J0_BC_EU = 30
       JM_BC_EU = 11
    ELSEIF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       I0_BC_CH = 100
       IM_BC_CH = 33
       J0_BC_CH = 39
       JM_BC_CH = 35
       I0_BC_NA = 16
       IM_BC_NA = 41
       J0_BC_NA = 50
       JM_BC_NA = 31
       I0_BC_EU = 60
       IM_BC_EU = 33
       J0_BC_EU = 60
       JM_BC_EU = 21
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This needs to be updated for FlexGrid (mps, 4/13/19)
    IF ( State_Grid%NestedGrid ) THEN

#if   defined ( EXCHANGE_4x5_CH )
       I0_BC          = 50
       J0_BC          = 20
       I2_BC          = 67
       J2_BC          = 37
       OFFSET_LONG    = 3.25
       OFFSET_LAT     = 1.5
       RATIO_LONG_EX  = 7.5d0
       RATIO_LAT_EX   = 8.0d0

#elif defined ( EXCHANGE_4x5_NA )
       I0_BC          = 8
       J0_BC          = 25
       I2_BC          = 29
       J2_BC          = 41
       OFFSET_LONG    = 3.25
       OFFSET_LAT     = 3.5
       RATIO_LONG_EX  = 7.5d0
       RATIO_LAT_EX   = 8.0d0

#elif defined ( EXCHANGE_4x5_EU )
       I0_BC          = 30
       J0_BC          = 30
       I2_BC          = 47
       J2_BC          = 41
       OFFSET_LONG    = 3.25
       OFFSET_LAT     = 3.5
       RATIO_LONG_EX  = 7.5d0
       RATIO_LAT_EX   = 8.0d0

#elif defined ( EXCHANGE_2x25_CH )
       I0_BC          = 100
       J0_BC          = 39
       I2_BC          = 133
       J2_BC          = 74
       OFFSET_LONG    = 1.375
       OFFSET_LAT     = 3.5
       RATIO_LONG_EX  = 3.75d0
       RATIO_LAT_EX   = 4.0d0

#elif defined ( EXCHANGE_2x25_NA )
       I0_BC          = 16
       J0_BC          = 50
       I2_BC          = 57
       J2_BC          = 81
       OFFSET_LONG    = 1.375
       OFFSET_LAT     = 1.5
       RATIO_LONG_EX  = 3.75d0
       RATIO_LAT_EX   = 4.0d0

#elif defined ( EXCHANGE_2x25_EU )
       I0_BC          = 60
       J0_BC          = 60
       I2_BC          = 93
       J2_BC          = 81
       OFFSET_LONG    = 1.375
       OFFSET_LAT     =  1.5
       RATIO_LONG_EX  = 3.75d0
       RATIO_LAT_EX   = 4.0d0

#endif
       I1_BC = I0_BC + 1
       J1_BC = J0_BC + 1
       IM_BC = I2_BC - I0_BC
       JM_BC = J2_BC - J0_BC

       ! generate regrid map
       ! CHECK: orientation problem...
       ALLOCATE( REGRID_MAT_L(IM_BC,State_Grid%NX), STAT=ALLOCATE_RES )
       IF ( ALLOCATE_RES /= 0 ) THEN
          CALL GEOS_CHEM_STOP
       ENDIF

       ALLOCATE( REGRID_MAT_L_TRANS(State_Grid%NX, IM_BC), STAT=ALLOCATE_RES )
       IF ( ALLOCATE_RES /= 0 ) THEN
          CALL GEOS_CHEM_STOP
       ENDIF

       ALLOCATE( REGRID_MAT_R(State_Grid%NY, JM_BC), STAT=ALLOCATE_RES )
       IF ( ALLOCATE_RES /= 0 ) THEN
          CALL GEOS_CHEM_STOP
       ENDIF

       CALL GEN_REGRID_MAT( State_Grid%NX, IM_BC, RATIO_LONG_EX, &
                            -OFFSET_LONG, REGRID_MAT_L_TRANS)
       CALL GEN_REGRID_MAT( State_Grid%NY, JM_BC, RATIO_LAT_EX, &
                            -OFFSET_LAT,  REGRID_MAT_R )

       REGRID_MAT_L = TRANSPOSE( REGRID_MAT_L_TRANS )

       IF(ALLOCATED(REGRID_MAT_L_TRANS)) DEALLOCATE( REGRID_MAT_L_TRANS )
       ALLOCATE( TMP_WINDOW(IM_BC, JM_BC, State_Grid%NZ), STAT=ALLOCATE_RES )
       IF ( ALLOCATE_RES /= 0 ) THEN
          CALL GEOS_CHEM_STOP
       ENDIF

    ENDIF ! NestedGrid
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF ( .not. State_Grid%NestedGrid) THEN

#if defined ( EXCHANGE_4x5_CH ) || defined ( EXCHANGE_2x25_CH )
       !ALLOCATE( TMP_WINDOW_CH(IM_BC_CH, JM_BC_CH, State_Grid%NZ), &
       !          STAT=ALLOCATE_RES )
       !IF ( ALLOCATE_RES /= 0 ) THEN
       !   CALL GEOS_CHEM_STOP
       !ENDIF
#endif

#if defined ( EXCHANGE_4x5_NA ) || defined ( EXCHANGE_2x25_NA )
       !ALLOCATE( TMP_WINDOW_NA(IM_BC_NA, JM_BC_NA, State_Grid%NZ), &
       !          STAT=ALLOCATE_RES )
       !IF ( ALLOCATE_RES /= 0 ) THEN
       !   CALL GEOS_CHEM_STOP
       !ENDIF
#endif

#if defined ( EXCHANGE_4x5_EU ) || defined ( EXCHANGE_2x25_EU )
       ALLOCATE( TMP_WINDOW_EU(IM_BC_EU, JM_BC_EU, State_Grid%NZ), &
                 STAT=ALLOCATE_RES )
       IF ( ALLOCATE_RES /= 0 ) THEN
          CALL GEOS_CHEM_STOP
       ENDIF
#endif

    ENDIF

  END SUBROUTINE INIT_EXCHANGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_exchange
!
! !DESCRIPTION: Cleanup routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_EXCHANGE( Input_Opt )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    IF ( Input_Opt%amIRoot ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@ CLEANUP_EXCHANGE"
    ENDIF

    IF( ALLOCATED( TMP_WINDOW     ) ) DEALLOCATE( TMP_WINDOW    )
    IF( ALLOCATED( REGRID_MAT_L   ) ) DEALLOCATE( REGRID_MAT_L  )
    IF( ALLOCATED( REGRID_MAT_R   ) ) DEALLOCATE( REGRID_MAT_R  )
    !IF( ALLOCATED( TMP_WINDOW_CH  ) ) DEALLOCATE( TMP_WINDOW_CH )
    !IF( ALLOCATED( TMP_WINDOW_NA  ) ) DEALLOCATE( TMP_WINDOW_NA )
    IF( ALLOCATED( TMP_WINDOW_EU  ) ) DEALLOCATE( TMP_WINDOW_EU )

    CALL NO_MORE_LOCK()

  END SUBROUTINE CLEANUP_EXCHANGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ex_dump_for_global
!
! !DESCRIPTION: ?
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EX_DUMP_FOR_GLOBAL( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )
!
! !USES:
!
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), POINTER :: SPC(:,:,:,:)
    INTEGER           :: L, I,IU_RST
    LOGICAL           :: prtDebug

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! NOTE: Species units are in kg/kg dry which may be incompatible with
    ! functionality prior to tracer removal. Validation needed by 3rd party
    ! developers (ewl, 8/15/16)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SPC      => State_Chm%Species
    prtDebug =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    IU_RST=2
    IF ( Input_Opt%amIRoot ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_DUMP_FOR_GLOBAL"
    ENDIF

    WRITE( 6, 100 ) TRIM(DUMP_CSPEC_FNAME)
100 FORMAT( '     - EXCHANGE: Writing ', A )
    CALL OPEN_BIN_FILE_FOR_WRITE( TRIM(DUMP_CSPEC_FNAME) )

    IF ( Input_Opt%amIRoot ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_REGRID_UP_SPECIES --Editing"
    ENDIF

    DO L = 1, State_Chm%nSpecies
       CALL EX_REGRID_UP( State_Grid%NX, State_Grid%NY, IM_BC, JM_BC, &
                          State_Grid%NZ, SPC(:,:,:,L), TMP_WINDOW )
       CALL WRITE_3D_REAL8_ARRAY( IM_BC, JM_BC, State_Grid%NZ, TMP_WINDOW )
    ENDDO
    CLOSE( 2 )

    CALL OPEN_BIN_FILE_FOR_WRITE( TRIM(DUMP_TRACER_FNAME) )
    IF ( Input_Opt%amIRoot ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_REGRID_UP_TRACER --Editing"
    ENDIF

    DO L = 1, State_Chm%nAdvect
       CALL EX_REGRID_UP( State_Grid%NX, State_Grid%NY, IM_BC, JM_BC, &
                          State_Grid%NZ, SPC(:,:,:,L), TMP_WINDOW )
       CALL WRITE_3D_REAL8_ARRAY( IM_BC, JM_BC, State_Grid%NZ, TMP_WINDOW )
    ENDDO
    CLOSE( IU_RST )

    SPC => NULL()
    IF ( prtDebug ) CALL DEBUG_MSG( '### WRITE_EXCHANGE: wrote file' )

  END SUBROUTINE EX_DUMP_FOR_GLOBAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ex_read_and_apply_feedback
!
! !DESCRIPTION: ?
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EX_READ_AND_APPLY_FEEDBACK( Input_Opt, State_Chm, State_Grid, &
                                         State_Met, RC, DOMAIN )
!
! !USES:
!
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE FILE_MOD
    USE TRANSFER_MOD,       ONLY : TRANSFER_3D_yan
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
    INTEGER,        INTENT(IN)    :: DOMAIN
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), POINTER :: SPC(:,:,:,:)
    REAL(fp), POINTER :: CSPEC_FULL(:,:,:,:)
    REAL(fp)          :: TMP_WINDOW_CHT(IM_BC_CH,JM_BC_CH,State_Grid%NZ)
    REAL(fp)          :: TMP_WINDOW_NAT(IM_BC_NA,JM_BC_NA,State_Grid%NZ)
    REAL(fp)          :: TMP_WINDOW_EUT(IM_BC_EU,JM_BC_EU,State_Grid%NZ)
    REAL*4            :: ARRAYTEMPR(576,361,73)
    INTEGER           :: L,IU_RST,IOS
    INTEGER           :: NI,NJ,NK,I,J,K
    LOGICAL           :: prtDebug

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! NOTE: Species units are in kg/kg dry which may be incompatible with
    ! functionality prior to tracer removal. Validation needed by
    ! 3rd party developers (ewl, 8/15/16)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SPC      => State_Chm%Species
    prtDebug =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    IU_RST=2
    IOS=1
    NI=IM_BC_CH
    NJ=JM_BC_CH
    NK=State_Grid%NZ

    IF ( prtDebug ) THEN
       WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_READ_AND_APPLY_FEEDBACK"
       WRITE(*,*) "!!!!!!!!!!FILLING THE SPECIES BOUNDARY"
    ENDIF

    ! For Asia nested domain
    IF ( DOMAIN == 1 ) THEN
       WRITE( 6, 100 ) TRIM(LOAD_CH_TRACER_FNAME)
100    FORMAT( '     - EXCHANGE: reading ', A )
       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_CH_TRACER_FNAME) )
       DO L = 1,State_Chm%nAdvect
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_CHT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_CH=TMP_WINDOW_CHT

          SPC( I0_BC_CH + 2 : I0_BC_CH + IM_BC_CH-1, &
               J0_BC_CH + 2 : J0_BC_CH + JM_BC_CH-1, &
               State_Grid%NZ, L ) = &
               TMP_WINDOW_CH(2:IM_BC_CH-1,2:JM_BC_CH-1,State_Grid%NZ)
       ENDDO
       CLOSE( IU_RST )

       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_CH_CSPEC_FNAME) )
       DO L = 1,State_Chm%nSpecies
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_CHT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_CH=TMP_WINDOW_CHT
          SPC( I0_BC_CH + 2:I0_BC_CH + IM_BC_CH-1, &
               J0_BC_CH + 2:J0_BC_CH + JM_BC_CH-1, &
               State_Grid%NZ, L ) = &
              TMP_WINDOW_CH(2:IM_BC_CH-1,2:JM_BC_CH-1,State_Grid%NZ)
       ENDDO
       CLOSE( IU_RST )

    ! For North America domain
    ELSE IF ( DOMAIN == 2 ) THEN
       NI=IM_BC_NA
       NJ=JM_BC_NA
       WRITE( 6, 101 ) TRIM(LOAD_NA_CSPEC_FNAME)
101    FORMAT( '     - EXCHANGE: reading ', A )
       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_NA_TRACER_FNAME) )
       DO L = 1,State_Chm%nAdvect
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_NAT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_NA=TMP_WINDOW_NAT
          SPC( I0_BC_NA + 3 : I0_BC_NA + IM_BC_NA-2, &
               J0_BC_NA + 3 : J0_BC_NA + JM_BC_NA-2, &
               State_Grid%NZ, L ) = &
               TMP_WINDOW_NA(3:IM_BC_NA-2,3:JM_BC_NA-2,State_Grid%NZ)
       ENDDO
       CLOSE( IU_RST )

       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_NA_CSPEC_FNAME) )
       DO L = 1,State_Chm%nSpecies
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_NAT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_NA=TMP_WINDOW_NAT
          SPC( I0_BC_NA + 3:I0_BC_NA + IM_BC_NA-2, &
               J0_BC_NA + 3:J0_BC_NA + JM_BC_NA-2, &
               State_Grid%NZ, L ) = &
               TMP_WINDOW_NA(3:IM_BC_NA-2,3:JM_BC_NA-2,State_Grid%NZ)
       ENDDO
       CLOSE( IU_RST )

    ! For Europe domain
    ELSE IF ( DOMAIN == 3 ) THEN
       NI=IM_BC_EU
       NJ=JM_BC_EU
       WRITE( 6, 102 ) TRIM(LOAD_EU_CSPEC_FNAME)
102    FORMAT( '     - EXCHANGE: reading ', A )
       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_EU_TRACER_FNAME) )
       DO L = 1,State_Chm%nAdvect
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_EUT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_EU=TMP_WINDOW_EUT
          SPC( I0_BC_EU + 2 : I0_BC_EU + IM_BC_EU-1, &
               J0_BC_EU + 2 : J0_BC_EU + JM_BC_EU-1, &
               State_Grid%NZ, L ) = &
               REAL(TMP_WINDOW_EU(2:IM_BC_EU-1,2:JM_BC_EU-1,State_Grid%NZ),8)
       ENDDO
       !CLOSE( 2 )

       CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_EU_CSPEC_FNAME) )
       DO L = 1,State_Chm%nSpecies
          READ( IU_RST, IOSTAT=IOS ) &
              ( ( ( ARRAYTEMPR(I,J,K), I=1,NI ), J=1,NJ ), K=1,NK )
          CALL TRANSFER_3D_yan &
               (NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK),TMP_WINDOW_EUT)
          IF ( IOS /= 0 ) THEN
             WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
             CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
          ENDIF
          TMP_WINDOW_EU=TMP_WINDOW_EUT
          SPC( I0_BC_EU + 2:I0_BC_EU + IM_BC_EU-1, &
               J0_BC_EU + 2:J0_BC_EU + JM_BC_EU-1, &
               State_Grid%NZ, L ) = &
               REAL(TMP_WINDOW_EU(2:IM_BC_EU-1,2:JM_BC_EU-1,State_Grid%NZ),8)
       ENDDO
       !CLOSE( 2 )
    ENDIF

    IF ( prtDebug ) CALL DEBUG_MSG( '### read and apply fb: done' )

  END SUBROUTINE EX_READ_AND_APPLY_FEEDBACK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ex_regrid_up
!
! !DESCRIPTION:  This subroutine will be called within EX\_DUMP\_FOR\_GLOBAL,
!  Using module's global vars: REGRID\_MAT\_L, REGRID\_MAT\_R
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EX_REGRID_UP( NX_ORIG, NY_ORIG, NX_NEW, NY_NEW, &
                           NZ, ORIG_ARRAY, NEW_ARRAY )
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER,  INTENT(IN)     :: NX_ORIG, NY_ORIG, NX_NEW, NY_NEW, NZ
    REAL(fp), INTENT(IN)     :: ORIG_ARRAY(NX_ORIG, NY_ORIG, NZ)
    REAL(fp), INTENT(OUT)    :: NEW_ARRAY(NX_NEW, NY_NEW, NZ)
    INTEGER                  :: K

    ! Use direct matrix multiply method for now.
    !yanyy ,2014,06,10
    DO K = 1, NZ
       NEW_ARRAY(:,:,K) = MATMUL( REGRID_MAT_L, &
                          MATMUL( ORIG_ARRAY(:,:,K), REGRID_MAT_R ) )
    ENDDO

  END SUBROUTINE EX_REGRID_UP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ex_convert_units_cspec
!
! !DESCRIPTION: Converts molec/cm3/box to fakemass (e.g. convert cspecfull
!  to a masslike unit)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EX_CONVERT_UNITS_CSPEC( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, RC, FLAG )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    INTEGER,        INTENT(IN)    :: FLAG
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !RETURN VALUE:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), POINTER :: SPC(:,:,:,:)
    REAL(fp)          :: SRC(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    integer           :: L

    SPC => State_Chm%Species

    IF ( FLAG == 1 ) THEN
       DO L = 1, State_Chm%nSpecies
          SRC( :, :, : ) = SPC( :, :, :, L) !* State_Met%AIRVOL
          SPC( :, :, :, L) = SRC(:, :, : )
       ENDDO
    ELSE IF (FLAG == 2 ) THEN
       DO L = 1, State_Chm%nSpecies
          SRC( :, :, : ) = SPC( :, :, :, L) !/ State_Met%AIRVOL
          SPC( :, :, :, L) = SRC(:, :, : )
       ENDDO
    END IF
    NULLIFY(SPC)

  END SUBROUTINE EX_CONVERT_UNITS_CSPEC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wait_for_unlock
!
! !DESCRIPTION: Waits for MPI unlock.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WAIT_FOR_UNLOCK()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER   :: INTERVAL       = 1
    LOGICAL              :: RES            = .FALSE.
    CHARACTER*512        :: EXT_COMMAND    = ""

    WRITE (*,*) "@@@@@@@@@@@@@@@@@@ WAIT_FOR_UNLOCK"
    WRITE (*,*) "         Waiting for file : ", TRIM(UNLOCK_FNAME)
    WRITE ( EXT_COMMAND, '(A, A)' ) "rm -f ", TRIM(UNLOCK_FNAME)

    DO
       CALL ALERT_FOR_KILL_COMMAND()
       INQUIRE( FILE = TRIM(UNLOCK_FNAME), EXIST = RES )
       IF ( RES .EQV. .TRUE. ) THEN
          CALL SYSTEM( TRIM(EXT_COMMAND) )
          EXIT
       ELSE
          CALL SLEEP( INTERVAL )
       ENDIF
    ENDDO

  END SUBROUTINE WAIT_FOR_UNLOCK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: done_lock_me_up
!
! !DESCRIPTION: ?
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DONE_LOCK_ME_UP()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    WRITE (*,*) "@@@@@@@@@@@@@@@@@@@ DONE_LOCK_ME_UP"
    CALL TOUCH_FILE( TRIM(DONE_FNAME) )

  END SUBROUTINE DONE_LOCK_ME_UP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: no_more_lock
!
! !DESCRIPTION: ?
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NO_MORE_LOCK()
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    WRITE (*,*) "@@@@@@@@@@@@@@@@@@@ NO_MORE_LOCK"
    CALL TOUCH_FILE( TRIM(END_FNAME) )

  END SUBROUTINE NO_MORE_LOCK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alert_for_kill_command
!
! !DESCRIPTION: Stops if there is a kill command issued.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ALERT_FOR_KILL_COMMAND()
!
! !USES:
!
    USE ERROR_MOD
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL       :: RES         = .FALSE.
    CHARACTER*512 :: EXT_COMMAND = ''

    WRITE ( EXT_COMMAND, '(A, A)' ) 'rm -f ', TRIM( KILL_FNAME )
    INQUIRE( FILE = TRIM(KILL_FNAME), EXIST = RES )

    IF ( RES .EQV. .TRUE. ) THEN
       WRITE (*,*) "!!!!!!!!!!!!!!!! Got KILL Command : ", &
                   TRIM(KILL_FNAME), " , STOPPING everything."
       CALL SYSTEM( TRIM(EXT_COMMAND) )
       CALL GEOS_CHEM_STOP()
    ENDIF

  END SUBROUTINE ALERT_FOR_KILL_COMMAND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: touch_file
!
! !DESCRIPTION: Touches a file (i.e. updates the timestamp on disk).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TOUCH_FILE( FNAME )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER*512 :: EXT_COMMAND    = ''

    WRITE ( EXT_COMMAND, '(A, A)' ) 'touch ' , TRIM(FNAME)

    CALL SYSTEM( TRIM(EXT_COMMAND) )

  END SUBROUTINE TOUCH_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bin_file_for_read
!
! !DESCRIPTION: Opens a binary file for reading.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPEN_BIN_FILE_FOR_READ( FNAME )
!
! !USES:
!
    USE FILE_MOD
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                          :: IOS,IU_RST
    IU_RST = 2

    OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='OLD', &
          IOSTAT=IOS, FORM='UNFORMATTED')
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error opening filename=',TRIM(FNAME)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IU_RST, 'open_bin_file_for_read:1' )
    ENDIF

  END SUBROUTINE OPEN_BIN_FILE_FOR_READ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bin_file_for_write
!
! !DESCRIPTION: Opens a binary file for writing.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPEN_BIN_FILE_FOR_WRITE( FNAME )
!
! !USES:
!
    USE FILE_MOD
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: IOS,IU_RST

    IU_RST = 2

    OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='REPLACE', &
          IOSTAT=IOS, FORM='UNFORMATTED')
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error opening filename=',TRIM(FNAME)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IU_RST, 'open_bin_file_for_write:1' )
    ENDIF

  END SUBROUTINE OPEN_BIN_FILE_FOR_WRITE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: open_bin_file_for_append
!
! !DESCRIPTION: Opens a binary file for appending.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPEN_BIN_FILE_FOR_APPEND( FNAME )
!
! !USES:
!
    USE FILE_MOD
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IOS,IU_RST

    IU_RST=2


    OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='UNKNOWN', &
          IOSTAT=IOS, FORM='UNFORMATTED')
    IF ( IOS /= 0 ) THEN
       WRITE(6,*)'Error opening filename=',TRIM(FNAME)
       CALL FLUSH(6)
       CALL IOERROR( IOS, IU_RST, 'open_bin_file_for_append:1' )
    ENDIF

  END SUBROUTINE OPEN_BIN_FILE_FOR_APPEND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: close_bin_file
!
! !DESCRIPTION: Closes a binary file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLOSE_BIN_FILE()
!
    USE FILE_MOD
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER         :: IU_RST

    IU_RST   = 2

    CLOSE( IU_RST )

  END SUBROUTINE CLOSE_BIN_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_3d_real8_array
!
! !DESCRIPTION: Writes data to a 3-D REAL*8 array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WRITE_3D_REAL8_ARRAY( NI, NJ, NK, ARRAY )
!
! !USES:
!
    USE FILE_MOD
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: NI, NJ, NK
    REAL(fp), INTENT(IN) :: ARRAY(NI,NJ,NK)
    REAL(f4)             :: ARRAYTEMPW(NI,NJ,NK)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: IU_RST
    INTEGER :: I, J, L

    IU_RST = 2

    ARRAYTEMPW=REAL(ARRAY,4)
    WRITE ( IU_RST ) ( ( ( ARRAYTEMPW(I,J,L), I=1,NI ), J=1,NJ ), L=1,NK )

  END SUBROUTINE WRITE_3D_REAL8_ARRAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_3d_real8_array
!
! !DESCRIPTION: Reads data from 3-D REAL*8 array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_3D_REAL8_ARRAY( NI, NJ, NK, ARRAY )
!
! !USES:
!
    USE FILE_MOD
    USE TRANSFER_MOD, ONLY : TRANSFER_3D_yan
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: NI, NJ, NK
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: ARRAY(NI,NJ,NK)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*4                :: ARRAYTEMPR(576,361,73)
    REAL(fp), ALLOCATABLE :: ARRAYT(:,:,:)
    INTEGER               :: IU_RST, IOS
    INTEGER               :: I, J, L

    IU_RST = 2
    IOS    = 1
    ARRAYTEMPR(:,:,:) = 0e0
    READ( IU_RST, IOSTAT=IOS ) &
        ( ( ( ARRAYTEMPR(I,J,L), I=1,NI ), J=1,NJ ), L=1,NK )
    CALL TRANSFER_3D_yan(NI,NJ,NK,ARRAYTEMPR(1:NI,1:NJ,1:NK), ARRAYT )
    ARRAY(1:NI,1:NJ,1:NK)=ARRAYT(1:NI,1:NJ,1:NK)
    IF ( IOS /= 0 ) THEN
       WRITE(6,*) 'EXCHANGE: READ ERROR' !, TRIM(FNAME)
       CALL IOERROR( IOS, IU_RST, 'EXCHANGE_MOD:READ')
    ENDIF

  END SUBROUTINE READ_3D_REAL8_ARRAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gen_regrid_mat
!
! !DESCRIPTION: Regridding utility.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEN_REGRID_MAT( N1, N2, RATIO, OFFSET, RES )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: N1, N2
    REAL(fp), INTENT(IN)  :: RATIO, OFFSET
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: RES(N1, N2)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Default for Right Array.
    REAL(fp) :: BEG_POS, END_POS
    INTEGER  :: M, N

    RES = 0d0
    DO N = 1, N2
       BEG_POS = (N-1) * RATIO + 1 + OFFSET
       END_POS = N * RATIO + 1 + OFFSET
       DO M = MAX(FLOOR(BEG_POS), 1), MIN(FLOOR(END_POS), N1)
          RES(M, N) = MIN(DBLE(M+1), END_POS) - MAX(DBLE(M), BEG_POS)
       END DO
    END DO

  END SUBROUTINE GEN_REGRID_MAT
!EOC
END MODULE EXCHANGE_MOD
#endif
