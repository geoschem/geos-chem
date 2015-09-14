#if defined( EXCHANGE )
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
!  (1) Establish the frame of the data exchange between models 
!       (kuangye,yanyy, 6/10/12)
!  (2) Now update for version 9-02 (yanyy, 6/21/14)
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
      REAL*8, ALLOCATABLE        :: TMP_WINDOW(:,:,:)
      REAL*8, ALLOCATABLE        :: TMP_WINDOW_CH(:,:,:)
      REAL*8, ALLOCATABLE        :: TMP_WINDOW_NA(:,:,:)
      REAL*8, ALLOCATABLE        :: TMP_WINDOW_EU(:,:,:)

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
      REAL*8                     :: RATIO_LONG_EX, RATIO_LAT_EX
      ! the parameters of boundary when regrid up global domain
      REAL*8                     :: OFFSET_LONG, OFFSET_LAT
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
      REAL*8, ALLOCATABLE        :: REGRID_MAT_L(:,:)
      REAL*8, ALLOCATABLE        :: REGRID_MAT_R(:,:)

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
      SUBROUTINE EXCHANGE_GLOBAL_POST( am_I_Root, Input_Opt,
     &                                 State_Met, State_Chm, RC )
!
! !USES:
!
      USE TPCORE_BC_MOD
      USE DAO_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      INTEGER         :: N_TRACERS
      REAL*8          :: TCVV(Input_Opt%N_TRACERS)
      REAL*8, POINTER :: STT(:,:,:,:)

      N_TRACERS =  Input_Opt%N_TRACERS
      TCVV      =  Input_Opt%TCVV(1:N_TRACERS)
      STT       => State_Chm%Tracers

      WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_GLOBAL_POST"
      CALL DONE_LOCK_ME_UP()
      CALL WAIT_FOR_UNLOCK()

      IF ( IS_EXCHANGE ) THEN

         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT )  !(v/v ->kg)
         CALL EX_CONVERT_UNITS_CSPEC
     & ( am_I_Root, Input_Opt, State_Met, State_Chm, RC,
     &  1, ILONG, ILAT, IPVERT) !mol/cm3 -> mol

# if  defined ( EXCHANGE_4x5_CH ) || defined ( EXCHANGE_2x25_CH )
         CALL EX_READ_AND_APPLY_FEEDBACK
     &  ( am_I_Root, Input_Opt,State_Met, State_Chm, RC, 1 )
# endif

# if  defined ( EXCHANGE_4x5_NA ) || defined ( EXCHANGE_2x25_NA )
          CALL EX_READ_AND_APPLY_FEEDBACK
     &  ( am_I_Root, Input_Opt,State_Met, State_Chm, RC, 2 )
# endif

# if  defined ( EXCHANGE_4x5_EU ) || defined ( EXCHANGE_2x25_EU )
         CALL EX_READ_AND_APPLY_FEEDBACK
     &  ( am_I_Root, Input_Opt,State_Met, State_Chm, RC, 3 )
# endif

         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT )  !(kg -> v/v)
         CALL EX_CONVERT_UNITS_CSPEC
     & ( am_I_Root, Input_Opt, State_Met, State_Chm, RC,
     &   2, ILONG, ILAT, IPVERT) !mol -> mol/cm3

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

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
      SUBROUTINE EXCHANGE_NESTED_POST( am_I_Root, Input_Opt,
     &                                 State_Met, State_Chm, RC )
!
! !USES:
!
      USE DAO_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: N_TRACERS
      REAL*8           :: TCVV(Input_Opt%N_TRACERS)
      REAL*8,  POINTER :: STT(:,:,:,:)

      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)
      STT => State_Chm%Tracers

      WRITE (*,*) "@@@@@@@@@@@@@@ EXCHANGE_NESTED_POST"

      CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, State_Met%AD, STT )  !(V/V -> KG)

      CALL EX_CONVERT_UNITS_CSPEC
     &  ( am_I_Root, Input_Opt, State_Met, State_Chm, RC,
     &    1,ILONG,ILAT,IPVERT) ! mol/cm3->mol

      CALL EX_DUMP_FOR_GLOBAL
     &        (am_I_Root, Input_Opt,State_Met, State_Chm,RC)

      CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, State_Met%AD, STT )  !(KG ->V/V)

      CALL EX_CONVERT_UNITS_CSPEC
     &   ( am_I_Root, Input_Opt, State_Met, State_Chm, RC,
     &     2,ILONG,ILAT,IPVERT) !mol -> mol/cm3

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
      SUBROUTINE INIT_EXCHANGE( am_I_Root )
!
! !USES:
!
      USE ERROR_MOD
      USE COMODE_MOD
      USE TRACER_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)       :: am_I_Root
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Now only print output on the root CPU
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                   :: ALLOCATE_RES
      REAL*8, ALLOCATABLE       :: REGRID_MAT_L_TRANS(:,:)

      IF ( am_I_Root ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@ INIT_EXCHANGE"
      ENDIF

#if defined ( GRID4x5 ) || defined ( GRID2x25 )

      IS_GLOBAL_DOMAIN = .TRUE.
      IS_EXCHANGE = .FALSE.
      MY_NAME = 'GLOBAL'
      NESTED_CH_NAME = 'NESTED_CH'
      NESTED_NA_NAME = 'NESTED_NA'
      NESTED_EU_NAME = 'NESTED_EU'

      IF ( am_I_Root ) THEN
      WRITE ( LOAD_CH_CSPEC_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.cspec.', TRIM(NESTED_CH_NAME), '.dump'
      WRITE ( LOAD_CH_TRACER_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.tracer.', TRIM(NESTED_CH_NAME), '.dump'
      WRITE ( LOAD_NA_CSPEC_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.cspec.', TRIM(NESTED_NA_NAME), '.dump'
      WRITE ( LOAD_NA_TRACER_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.tracer.', TRIM(NESTED_NA_NAME), '.dump'
      WRITE ( LOAD_EU_CSPEC_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.cspec.', TRIM(NESTED_EU_NAME), '.dump'
      WRITE ( LOAD_EU_TRACER_FNAME, '(A, A, A, A)' )TRIM(EXCHANGE_DIR),
     &            '/dump.tracer.', TRIM(NESTED_EU_NAME), '.dump'
      ENDIF

#elif defined ( NESTED_CH )
      IS_GLOBAL_DOMAIN = .FALSE.
      MY_NAME = 'NESTED_CH'

#elif defined ( NESTED_NA )
      IS_GLOBAL_DOMAIN = .FALSE.
      MY_NAME = 'NESTED_NA'

#elif defined ( NESTED_EU )
      IS_GLOBAL_DOMAIN = .FALSE.
      MY_NAME = 'NESTED_EU'
#else

      WRITE(*,*) "You have to define a domain"
      CALL GEOS_CHEM_STOP

#endif

      IF ( am_I_Root ) THEN
      WRITE ( DUMP_CSPEC_FNAME, '(A, A, A, A)' ) TRIM(EXCHANGE_DIR),
     &            '/dump.cspec.', TRIM(MY_NAME), '.dump'
      WRITE ( DUMP_TRACER_FNAME, '(A, A, A, A)' ) TRIM(EXCHANGE_DIR),
     &            '/dump.tracer.', TRIM(MY_NAME), '.dump'


      WRITE ( UNLOCK_FNAME, '(A, A, A, A)' ) TRIM(LOCK_DIR), '/key.',
     &                                       TRIM(MY_NAME), '.key'
      WRITE ( DONE_FNAME,   '(A, A, A, A)' ) TRIM(LOCK_DIR), '/done.',
     &                                       TRIM(MY_NAME), '.done'
      WRITE ( KILL_FNAME,   '(A, A, A, A)' ) TRIM(LOCK_DIR), '/kill.',
     &                                       TRIM(MY_NAME), '.kill'
      WRITE ( END_FNAME,    '(A, A, A, A)' ) TRIM(LOCK_DIR), '/end.',
     &                                       TRIM(MY_NAME), '.end'
      ENDIF

      ! size information
#if   defined ( GRID4x5 )
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

#elif defined ( GRID2x25 )
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

#endif

#if   defined(NESTED_CH) || defined(NESTED_NA) || defined(NESTED_EU)

#if   defined ( EXCHANGE_4x5_CH )
      I0_BC       = 50
      J0_BC       = 20
      I2_BC       = 67
      J2_BC       = 37
      OFFSET_LONG    = 3.25
      OFFSET_LAT     = 1.5
      RATIO_LONG_EX  = 7.5d0
      RATIO_LAT_EX   = 8.0d0

#elif   defined ( EXCHANGE_4x5_NA )
      I0_BC       = 8
      J0_BC       = 25
      I2_BC       = 29
      J2_BC       = 41
      OFFSET_LONG    = 3.25
      OFFSET_LAT     = 3.5
      RATIO_LONG_EX  = 7.5d0
      RATIO_LAT_EX   = 8.0d0

#elif   defined ( EXCHANGE_4x5_EU )
      I0_BC       = 30
      J0_BC       = 30
      I2_BC       = 47
      J2_BC       = 41
      OFFSET_LONG    = 3.25
      OFFSET_LAT     = 3.5
      RATIO_LONG_EX  = 7.5d0
      RATIO_LAT_EX   = 8.0d0

#elif defined ( EXCHANGE_2x25_CH )
      I0_BC   = 100
      J0_BC    = 39
      I2_BC   = 133
      J2_BC    = 74
      OFFSET_LONG = 1.375
      OFFSET_LAT  = 3.5
      RATIO_LONG_EX  = 3.75d0
      RATIO_LAT_EX   = 4.0d0

#elif defined ( EXCHANGE_2x25_NA )
      I0_BC   = 16
      J0_BC    = 50
      I2_BC   = 57
      J2_BC    = 81
      OFFSET_LONG = 1.375
      OFFSET_LAT = 1.5
      RATIO_LONG_EX  = 3.75d0
      RATIO_LAT_EX   = 4.0d0

#elif defined ( EXCHANGE_2x25_EU )
      I0_BC   = 60
      J0_BC    = 60
      I2_BC   = 93
      J2_BC    = 81
      OFFSET_LONG = 1.375
      OFFSET_LAT =  1.5
      RATIO_LONG_EX  = 3.75d0
      RATIO_LAT_EX   = 4.0d0

#endif
      I1_BC = I0_BC + 1
      J1_BC = J0_BC + 1
      IM_BC = I2_BC - I0_BC
      JM_BC = J2_BC - J0_BC

      ! generate regrid map
      ! CHECK: orientation problem...
      ALLOCATE( REGRID_MAT_L(IM_BC,ILONG), STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF

      ALLOCATE( REGRID_MAT_L_TRANS(ILONG, IM_BC), STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF

      ALLOCATE( REGRID_MAT_R(ILAT, JM_BC), STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF

      CALL GEN_REGRID_MAT( ILONG, IM_BC, RATIO_LONG_EX,
     &            -OFFSET_LONG, REGRID_MAT_L_TRANS)
      CALL GEN_REGRID_MAT( ILAT, JM_BC, RATIO_LAT_EX,
     &            -OFFSET_LAT, REGRID_MAT_R )

      REGRID_MAT_L = TRANSPOSE( REGRID_MAT_L_TRANS )

      IF(ALLOCATED(REGRID_MAT_L_TRANS)) DEALLOCATE( REGRID_MAT_L_TRANS )
      ALLOCATE( TMP_WINDOW(IM_BC, JM_BC, IPVERT),
     &          STAT=ALLOCATE_RES )

      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF
#endif

#if   defined ( GRID4x5 ) || defined ( GRID2x25 )

#if   defined ( EXCHANGE_4x5_CH ) || defined ( EXCHANGE_2x25_CH )
      ALLOCATE( TMP_WINDOW_CH(IM_BC_CH, JM_BC_CH, IPVERT),
     &          STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF
#endif

#if   defined ( EXCHANGE_4x5_NA ) || defined ( EXCHANGE_2x25_NA )
      ALLOCATE( TMP_WINDOW_NA(IM_BC_NA, JM_BC_NA, IPVERT),
     &          STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF
#endif

#if   defined ( EXCHANGE_4x5_EU ) || defined ( EXCHANGE_2x25_EU )
      ALLOCATE( TMP_WINDOW_EU(IM_BC_EU, JM_BC_EU, IPVERT),
     &          STAT=ALLOCATE_RES )
      IF ( ALLOCATE_RES /= 0 ) THEN
         CALL GEOS_CHEM_STOP
      ENDIF
#endif
#endif

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
      SUBROUTINE CLEANUP_EXCHANGE( am_I_Root )
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Now only write output on the root CPU
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      IF ( am_I_Root ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@ CLEANUP_EXCHANGE"
      ENDIF

      IF( ALLOCATED( TMP_WINDOW  ) )  DEALLOCATE( TMP_WINDOW )
      IF( ALLOCATED( REGRID_MAT_L ) ) DEALLOCATE( REGRID_MAT_L )
      IF( ALLOCATED( REGRID_MAT_R ) ) DEALLOCATE( REGRID_MAT_R )
      IF( ALLOCATED( TMP_WINDOW_CH  ) )  DEALLOCATE( TMP_WINDOW_CH )
      IF( ALLOCATED( TMP_WINDOW_NA  ) )  DEALLOCATE( TMP_WINDOW_NA )
      IF( ALLOCATED( TMP_WINDOW_EU  ) )  DEALLOCATE( TMP_WINDOW_EU )

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
      SUBROUTINE EX_DUMP_FOR_GLOBAL( am_I_Root, Input_Opt,
     &                               State_Met, State_Chm, RC )
!
! !USES:
!
      USE ERROR_MOD
      USE COMODE_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8,         POINTER       :: STT(:,:,:,:)
      REAL*8,         POINTER       :: CSPEC_FULL(:,:,:,:)
      INTEGER                       :: L, I, N_TRACERS
      LOGICAL                       :: LPRT

      STT        => State_Chm%Tracers
      CSPEC_FULL => State_Chm%Species
      N_TRACERS  =  Input_Opt%N_TRACERS
      LPRT       =  ( am_I_Root .and. Input_Opt%LPRT )
      IF ( am_I_Root ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_DUMP_FOR_GLOBAL"
      ENDIF

      WRITE( 6, 100 ) TRIM(DUMP_CSPEC_FNAME)
 100  FORMAT( '     - EXCHANGE: Writing ', A )
      CALL OPEN_BIN_FILE_FOR_WRITE( TRIM(DUMP_CSPEC_FNAME) )

      IF ( am_I_Root ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_REGRID_UP_SPECIES --Editing"
      ENDIF

      DO L = 1,NTSPEC(1)
         CALL EX_REGRID_UP( ILONG, ILAT,
     &                      IM_BC, JM_BC, IPVERT,
     &                      CSPEC_FULL(:,:,:,L), TMP_WINDOW )
         CALL WRITE_3D_REAL8_ARRAY( IM_BC, JM_BC, IPVERT,
     &                              TMP_WINDOW )
      ENDDO



      CALL OPEN_BIN_FILE_FOR_WRITE( TRIM(DUMP_TRACER_FNAME) )
      IF ( am_I_Root ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_REGRID_UP_TRACER --Editing"
      ENDIF

      DO L = 1,N_TRACERS
         CALL EX_REGRID_UP( ILONG, ILAT,
     &                      IM_BC, JM_BC, IPVERT,
     &                      STT(:,:,:,L), TMP_WINDOW )
         CALL WRITE_3D_REAL8_ARRAY( IM_BC, JM_BC, IPVERT,
     &                              TMP_WINDOW )
      ENDDO


      IF ( LPRT ) CALL DEBUG_MSG( '### WRITE_EXCHANGE: wrote file' )
      ! Return to calling program
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
      SUBROUTINE EX_READ_AND_APPLY_FEEDBACK( am_I_Root, Input_Opt,
     &                         State_Met, State_Chm, RC, DOMAIN )
!
! !USES:
!
      USE ERROR_MOD
      USE TRACER_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!

      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt
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
!  06 Jan 2015 - R. Yantosca - Now get LPRT from Input_Opt
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:


      REAL*8,         POINTER       :: STT(:,:,:,:)
      REAL*8,         POINTER       :: CSPEC_FULL(:,:,:,:)
      INTEGER                       :: L, N_TRACERS
      LOGICAL                       :: LPRT

      STT        => State_Chm%Tracers
      CSPEC_FULL => State_Chm%Species
      N_TRACERS  = Input_Opt%N_TRACERS
      LPRT       =  ( am_I_Root .and. Input_Opt%LPRT ) 

      IF ( LPRT ) THEN
         WRITE (*,*) "@@@@@@@@@@@@@@@@@@ EX_READ_AND_APPLY_FEEDBACK"
         WRITE(*,*) "!!!!!!!!!!FILLING THE SPECIES BOUNDARY"
      ENDIF

      !For Asia nested domain
      IF ( DOMAIN == 1 ) THEN
          WRITE( 6, 100 ) TRIM(LOAD_CH_CSPEC_FNAME)
 100      FORMAT( '     - EXCHANGE: reading ', A )
         CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_CH_TRACER_FNAME) )
         DO L = 1,N_TRACERS
              CALL READ_3D_REAL8_ARRAY( IM_BC_CH, JM_BC_CH, IPVERT,
     &                             TMP_WINDOW_CH )
              STT( I0_BC_CH + 2 : I0_BC_CH + IM_BC_CH-1,
     &             J0_BC_CH + 2 : J0_BC_CH + JM_BC_CH-1,
     &             IPVERT,L) =
     &        TMP_WINDOW_CH(2:IM_BC_CH-1,2:JM_BC_CH-1,IPVERT)
         ENDDO

          CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_CH_CSPEC_FNAME) )
          DO L = 1,NTSPEC(1)
             CALL READ_3D_REAL8_ARRAY( IM_BC_CH, JM_BC_CH, IPVERT,
     &                             TMP_WINDOW_CH )
             CSPEC_FULL(I0_BC_CH + 2:I0_BC_CH + IM_BC_CH-1,
     &                  J0_BC_CH + 2:J0_BC_CH + JM_BC_CH-1,
     &                  IPVERT,L) =
     &       TMP_WINDOW_CH(2:IM_BC_CH-1,2:JM_BC_CH-1,IPVERT)
          ENDDO

      !For North America domain
      ELSE IF ( DOMAIN == 2 ) THEN
          WRITE( 6, 101 ) TRIM(LOAD_NA_CSPEC_FNAME)
 101      FORMAT( '     - EXCHANGE: reading ', A )
         CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_NA_TRACER_FNAME) )
         DO L = 1,N_TRACERS
              CALL READ_3D_REAL8_ARRAY( IM_BC_NA, JM_BC_NA, IPVERT,
     &                             TMP_WINDOW_NA )
              STT( I0_BC_NA + 3 : I0_BC_NA + IM_BC_NA-2,
     &             J0_BC_NA + 3 : J0_BC_NA + JM_BC_NA-2,
     &             IPVERT,L) =
     &        TMP_WINDOW_NA(3:IM_BC_NA-2,3:JM_BC_NA-2,IPVERT)
         ENDDO

          CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_NA_CSPEC_FNAME) )
          DO L = 1,NTSPEC(1)
             CALL READ_3D_REAL8_ARRAY( IM_BC_NA, JM_BC_NA, IPVERT,
     &                             TMP_WINDOW_NA )
             CSPEC_FULL(I0_BC_NA + 3:I0_BC_NA + IM_BC_NA-2,
     &                  J0_BC_NA + 3:J0_BC_NA + JM_BC_NA-2,
     &                  IPVERT,L) =
     &       TMP_WINDOW_NA(3:IM_BC_NA-2,3:JM_BC_NA-2,IPVERT)
          ENDDO

      !For Europe domain
      ELSE IF ( DOMAIN == 3 ) THEN
          WRITE( 6, 102 ) TRIM(LOAD_EU_CSPEC_FNAME)
  102     FORMAT( '     - EXCHANGE: reading ', A )
         CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_EU_TRACER_FNAME) )
         DO L = 1,N_TRACERS
              CALL READ_3D_REAL8_ARRAY( IM_BC_EU, JM_BC_EU, IPVERT,
     &                             TMP_WINDOW_EU )
              STT( I0_BC_EU + 2 : I0_BC_EU + IM_BC_EU-1,
     &             J0_BC_EU + 2 : J0_BC_EU + JM_BC_EU-1,
     &             IPVERT,L) =
     &        TMP_WINDOW_EU(2:IM_BC_EU-1,2:JM_BC_EU-1,IPVERT)
         ENDDO

         CALL OPEN_BIN_FILE_FOR_READ( TRIM(LOAD_EU_CSPEC_FNAME) )
         DO L = 1,NTSPEC(1)
             CALL READ_3D_REAL8_ARRAY( IM_BC_EU, JM_BC_EU, IPVERT,
     &                             TMP_WINDOW_EU )
             CSPEC_FULL(I0_BC_EU + 2:I0_BC_EU + IM_BC_EU-1,
     &                  J0_BC_EU + 2:J0_BC_EU + JM_BC_EU-1,
     &                  IPVERT,L) =
     &       TMP_WINDOW_EU(2:IM_BC_EU-1,2:JM_BC_EU-1,IPVERT)
         ENDDO
      ENDIF

      IF ( LPRT ) CALL DEBUG_MSG( '### read and apply fb: done' )
      !Return to calling program
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
      SUBROUTINE EX_REGRID_UP( NX_ORIG, NY_ORIG, NX_NEW, NY_NEW,
     &                         NZ, ORIG_ARRAY, NEW_ARRAY )
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, INTENT(IN)     :: NX_ORIG, NY_ORIG, NX_NEW, NY_NEW, NZ
      REAL*8,  INTENT(IN)     :: ORIG_ARRAY(NX_ORIG, NY_ORIG, NZ)
      REAL*8,  INTENT(OUT)    :: NEW_ARRAY(NX_NEW, NY_NEW, NZ)
      ! Local vars
      INTEGER                 :: K


      ! Use direct matrix multiply method for now.
      !yanyy ,2014,06,10
      DO K = 1, NZ
         NEW_ARRAY(:,:,K) = MATMUL( REGRID_MAT_L,
     &            MATMUL( ORIG_ARRAY(:,:,K), REGRID_MAT_R ) )
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
      SUBROUTINE EX_CONVERT_UNITS_CSPEC( am_I_Root, Input_Opt,
     &                 State_Met, State_Chm, RC, FLAG, NX, NY, NZ )
!
! !USES:
!
      USE DAO_MOD
      USE COMODE_MOD
      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! Allocatable arrays for SMVGEAR solver
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, POINTER :: CSPEC_FULL(:,:,:,:)
      INTEGER, INTENT(IN)           :: NX, NY, NZ, FLAG
      REAL*8                        :: SRC(NX, NY, NZ)
      integer                       :: L

      CSPEC_FULL => State_Chm%Species

      IF ( FLAG == 1 ) THEN
         DO L = 1, NTSPEC(1)
            SRC( :, :, : ) = CSPEC_FULL( :, :, :, L) * State_Met%AIRVOL
            CSPEC_FULL( :, :, :, L) = SRC(:, :, : )
         ENDDO
      ELSE IF (FLAG == 2 ) THEN
         DO L = 1, NTSPEC(1)
            SRC( :, :, : ) = CSPEC_FULL( :, :, :, L) / State_Met%AIRVOL
            CSPEC_FULL( :, :, :, L) = SRC(:, :, : )
         ENDDO
      END IF
      NULLIFY(CSPEC_FULL)

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

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
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL              :: RES         = .FALSE.
      CHARACTER*512        :: EXT_COMMAND = ''

      WRITE ( EXT_COMMAND, '(A, A)' ) 'rm -f ', TRIM( KILL_FNAME )
      INQUIRE( FILE = TRIM(KILL_FNAME), EXIST = RES )

      IF ( RES .EQV. .TRUE. ) THEN
         WRITE (*,*) "!!!!!!!!!!!!!!!! Got KILL Command : ",
     &               TRIM(KILL_FNAME), " , STOPPING everything."
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
      CHARACTER(LEN=*), INTENT(IN)     :: FNAME
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER*512                    :: EXT_COMMAND    = ''

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
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=*), INTENT(IN)     :: FNAME
      INTEGER                          :: IOS,IU_RST
      IU_RST = 2

      OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='OLD',
     &      IOSTAT=IOS, FORM='UNFORMATTED')
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
      CHARACTER(LEN=*), INTENT(IN)     :: FNAME
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                          :: IOS,IU_RST
      IU_RST = 2

      OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='REPLACE',
     &      IOSTAT=IOS, FORM='UNFORMATTED')
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
      CHARACTER(LEN=*), INTENT(IN)     :: FNAME
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                          :: IOS,IU_RST
      IU_RST=2


      OPEN( IU_RST,     FILE=TRIM(FNAME), STATUS='UNKNOWN',
     &      IOSTAT=IOS, FORM='UNFORMATTED')
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
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
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
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)     :: NI, NJ, NK
      REAL*8,           INTENT(IN)     :: ARRAY(NI,NJ,NK)
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER         :: IU_RST
      IU_RST   = 2

      WRITE ( IU_RST ) ARRAY

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
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)     :: NI, NJ, NK
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,           INTENT(OUT)    :: ARRAY(NI,NJ,NK)
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      INTEGER         :: IU_RST
      IU_RST   = 2

      READ ( IU_RST ) ARRAY

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
      INTEGER,          INTENT(IN)     :: N1, N2
      REAL*8,           INTENT(IN)     :: RATIO, OFFSET
!
! !OUTPUT PARAMETERS: 
!
      REAL*8,           INTENT(OUT)    :: RES(N1, N2)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Default for Right Array.
      ! LOCAL vars
      REAL*8                           :: BEG_POS, END_POS
      INTEGER                          :: M, N

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
