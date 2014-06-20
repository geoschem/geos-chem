!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_ch4_mod
!
! !DESCRIPTION: Module GLOBAL\_CH4\_MOD contains variables and routines for
!  simulating CH4 chemistry in the troposphere.
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_CH4_MOD
!
! !USES:
!
      USE inquireMod, ONLY : findFreeLUN

      IMPLICIT NONE
      PRIVATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Normally we will not save out the CH4_BUDGET diagnostics (esp. when
!%%% using the nested grid simulations) in order to save memory.  If you
!%%% want to use CH4_BUDGET, then uncomment the following line of code:
!%%% (kjw, bmy, 2/12/14)
!#define USE_CH4_BUDGET_DIAG 1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CH4_AVGTP
      PUBLIC :: CHEMCH4
      PUBLIC :: INIT_GLOBAL_CH4
      PUBLIC :: CLEANUP_GLOBAL_CH4
!
! !PUBLIC DATA MEMBERS:
!
      REAL*8, PARAMETER,   PUBLIC :: XNUMOL_CH4 = 6.0221d23 / 16d-3

#if defined( USE_CH4_BUDGET_DIAG )
      REAL*8, ALLOCATABLE, PUBLIC :: TCH4(:,:,:,:)
#endif
! 
! !REVISION HISTORY: 
!  17 Jan 2001- J. Wang, B. Duncan, R. Yantosca -- Initial version
!  (1 ) Merged routines from jsw's CH4 code  into "global_ch4_mod.f" 
!        (bmy, 1/16/01)
!  (2 ) XNUMOL_CH4 and TCH4 have to be public - all other variables can
!        be made private, so as not to conflict with other common-block
!        definitions (bmy, 1/17/01)
!  (3 ) Minor fixes from jsw added (jsw, bmy, 2/17/01)
!  (4 ) Removed some F90 module references from EMISSCH4 (bmy, 3/20/01)
!  (5 ) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) Fixes for binary punch file in READ_COPROD (bmy, 9/26/01)
!  (8 ) Removed obsolete code from READ_COPROD (bmy, 10/24/01)
!  (9 ) Minor bug fixes for compilation on ALPHA (bmy, 11/15/01)
!  (10) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (11) Now eliminate PS from the arg list to CH4_AVGTP (4/11/02)
!  (12) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (13) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (14) Now reference "file_mod.f".  Also removed obsolete code. (bmy, 6/27/02)
!  (15) Now references "pressure_mod.f" (bmy, 8/21/02)
!  (16) Now reference AD and T from "dao_mod.f".  Now reference "error_mod.f".
!        Remove obsolete code from various routines.  Remove reference to
!        header file "comtrid.h" -- it's not used. (bmy, 11/6/02)
!  (17) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (18) Now references "grid_mod.f" and "time_mod.f" (bmy, 3/27/03)
!  (19) Updates to GET_GLOBAL_CH4 (bmy, 7/1/03)
!  (20) Now references "directory_mod.f", "tracer_mod.f", and "diag_oh_mod.f"
!        (bmy, 7/20/04)
!  (21) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (22) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (23) Updated CH4 simulation (kjw, cph, ccarouge, 10/1/09)
!  (24) Added modifications for MERRA (bmy, 8/13/10)
!  08 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  12 Feb 2014 - K. Wecht    - Add modifications for 0.25 x 0.3125 NA grid
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!  20 Jun 2014 - R. Yantosca - Remove obsolete emissions code; we now use HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!     
! !DEFINED PARAMETERS:
!
      !========================================================================
      ! Module Variables:
      ! N_CH4      : Number of budget items in TCH4
      ! BAIRDENS   : Array for air density                      [molec/cm3]
      ! BOH        : Array for OH values                        [molec/cm3]
      ! COPROD     : Array for zonal mean P(CO)                 [v/v/s]
      ! CH4LOSS    : Array for monthly average CH4 loss freq    [1/s]
      ! PAVG       : Array for 24-h avg surface pressure        [mb]
      ! TAVG       : Array for 24-h avg temperature             [K]
      ! TCH4       : Array for CH4 budget (N_CH4 items)         [molec/box]
      ! NCMSALTS   : # of altitudes for CMS climatological OH
      ! NCMSLATS   : # of latitudes for CMS climatological OH
      ! CMSALTS    : Altitude values for CMS climatological OH
      ! CMSLATS    : Latitude values for CMS climatological OH
      ! AVGOH      : Array for CMS climatological OH            [molec/cm3]
      ! FMOL_CH4   : Molecular weight of CH4                    [kg/mole]
      ! XNUMOL_CH4 : Molecules CH4 / kg CH4                     [molec/kg]
      ! CH4_EMIS   : Array for CH4 Emissions                    [molec/cm2/s]
      !========================================================================

      ! Number of CH4 budget types 
      INTEGER, PARAMETER   :: N_CH4    = 12
      ! For Clarisa's Climatological OH
      INTEGER, PARAMETER   :: NSEAS    = 4
      INTEGER, PARAMETER   :: NCMSALTS = 7
      INTEGER, PARAMETER   :: NCMSLATS = 24
      REAL*8,  PARAMETER   :: FMOL_CH4 = 16d-3
!
! !LOCAL VARIABLES:
!
      ! Various arrays      
      REAL*8,  ALLOCATABLE :: BAIRDENS(:,:,:)
      REAL*8,  ALLOCATABLE :: BOH(:,:,:,:)
      REAL*8,  ALLOCATABLE :: COPROD(:,:,:)
      REAL*8,  ALLOCATABLE :: CH4LOSS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: PAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: TAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: CH4_EMIS(:,:,:)
      REAL*8,  ALLOCATABLE :: AVGOH(:,:,:)
      ! Online CH4 only
      REAL*8,  ALLOCATABLE :: SUMANTHRO(:,:)
      REAL*8,  ALLOCATABLE :: SUMBIOG(:,:)

      ! For Clarisa's Climatological OH      
      REAL*8               :: CMSALTS(NCMSALTS) =
     &    (/ 1000d0, 900d0, 800d0, 700d0, 500d0, 300d0, 200d0 /)

      REAL*8               :: CMSLATS(NCMSLATS) =
     &    (/ 90d0,  84d0,  76d0,  68d0,  60d0,  52d0,  44d0,  36d0, 
     &       28d0,  20d0,  12d0,   4d0,  -4d0, -12d0, -20d0, -28d0,
     &      -36d0, -44d0, -52d0, -60d0, -68d0, -76d0, -84d0, -90d0 /)

      REAL*8               :: TROPOCH4
      REAL*8               :: STRATOCH4

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_avgtp
!
! !DESCRIPTION: Subroutine CH4\_AVGTP gets the 24-h average surface pressure
!  and temperature needed for the CH4 simulation. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_AVGTP( State_Met )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_TS_DYN, GET_TS_CHEM
      USE TIME_MOD,           ONLY : GET_ELAPSED_MIN

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry and
!        placed into module "global_ch4_mod.f" by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_AVGTP is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed duplicate definition for NTDT, NMIN (bmy, 11/15/01)
!  (4 ) Removed PS from argument list.  Now use P(I,J)+PTOP instead of
!        PS, this ensures that we have consistency between P and AD.
!        (bmy, 4/11/02)
!  (5 ) Removed obsolete code (bmy, 6/27/02)
!  (6 ) Now uses GET_PCENTER from "pressure_mod.f" to return the pressure
!        at the midpoint of the box (I,J,L).  Also added parallel DO-loops.
!        Updated comments. (dsa, bdf, bmy, 8/21/02)
!  (7 ) Now reference T from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f" (bmy, 10/15/02)
!  (8 ) Removed NTDT, NMIN from the arg list.  Now uses functions GET_TS_DYN,
!        GET_TS_CHEM, and GET_ELAPSED_MIN from "time_mod.f" (bmy, 3/27/03)
!  (9 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NTDT, NMIN
      INTEGER               :: I, J, L, NTIMES, MNDT, K, M, N
      INTEGER               :: NNEW, NNCOUNT 
      REAL*8                :: Ptemp(IIPAR,JJPAR,LLPAR)
      
      !=================================================================
      ! CH4_AVGTP begins here!
      !=================================================================

      ! Get quantities from "time_mod.f"
      NTDT   = GET_TS_DYN() * 60
      NMIN   = GET_ELAPSED_MIN()
      MNDT   = NTDT  / 60 
      NTIMES = GET_TS_CHEM() / MNDT

      ! NTIMES is the number of dynamic timesteps in a chem timestep
      IF ( NMIN <= GET_TS_CHEM() ) NTIMES = NTIMES + 1

      ! At the start of the run...
      IF ( NMIN == 0 ) THEN

         ! Initialize NNEW
	 NNEW = 0

         ! No longer true (kjw, 11/17/11).
         !! Error check --  NCHEM has to be 1440 min
         !IF ( GET_TS_CHEM() /= 1440 ) THEN   
         !   WRITE(*,*) ' '
         !   WRITE(*,*) 'CH4-OH parameterization option (i.e., NSRCX=5)!' 
         !   WRITE(*,*) 'Use a chemistry time step = 24 hours'
         !   WRITE(*,*) '(i.e., NCHEM=1440 min.)'
         !   WRITE(*,*) ' '
         !   CALL GEOS_CHEM_STOP
         !ENDIF

         ! Error check -- need chem timestep (1440) to be divisible by 
         ! dyn timestep
         IF ( mod( GET_TS_CHEM(), MNDT ) /= 0 ) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'CH4-OH parameterization option (i.e., NSRCX=5)!'
            WRITE(*,*) 'The chemistry time step (i.e., 24 hours) is'
            WRITE(*,*) 'not evenly divisible by the meteorological'
            WRITE(*,*) 'data read-in time step (i.e., 6 hours).  This'
            WRITE(*,*) 'will mess up SR avgtp which calculates a 24-'
            WRITE(*,*) 'hour average temperature and pressure to be'
            WRITE(*,*) 'used by SR getinfo.'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF

         ! If NCHEM < NTDT then stop program.
         IF ( GET_TS_CHEM() < MNDT ) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'When using the CH4-OH parameterization'
            WRITE(*,*) 'option (i.e., NSRCX=5), take a 24-hour'
            WRITE(*,*) 'time step (i.e., NCHEM=1440 min.) because'
            WRITE(*,*) 'the OH parameterization produces a 24-hour'
            WRITE(*,*) 'average [OH]'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDIF

      !=================================================================
      ! If a new 24-hr period, set Pavg = 0, and reset NNEW, NCOUNT
      !=================================================================
      IF ( NNEW == 0 ) THEN 
         Pavg(:,:,:) = 0d0
         Tavg(:,:,:) = 0d0
	 NNEW        = 1
	 NNCOUNT     = 0
      ENDIF

      !=================================================================
      ! Archive quantities
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PTEMP )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                  
         ! Archive pressure
         Pavg(I,J,L) = Pavg(I,J,L) + GET_PCENTER(I,J,L)

         ! Archive temperature
         Tavg(I,J,L) = Tavg(I,J,L) + State_Met%T(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !================================================================
      ! Keep track to see if at end of NCHEM time step.
      ! If so, divide PAVG & TAVG by the number of times archived.
      !=================================================================
	
      NNCOUNT = NNCOUNT + 1

      IF ( NNCOUNT == NTIMES ) THEN
         Pavg(:,:,1:LLPAR) = Pavg(:,:,1:LLPAR) / DBLE( NTIMES )
         Tavg(:,:,1:LLPAR) = Tavg(:,:,1:LLPAR) / DBLE( NTIMES )
         NNEW              = 0
      ENDIF
      
      END SUBROUTINE CH4_AVGTP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemch4
!
! !DESCRIPTION: Subroutine CHEMCH4 computes the chemical loss of CH4
!  (sources - sinks). (jsw, bnd, bmy, 6/8/00, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMCH4( am_I_Root, Input_Opt,
     &                    State_Met, State_Chm, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_TAU0, READ_BPCH2, GET_MODELNAME
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT, GET_RES_EXT
      USE CMN_MOD
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD43
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, IT_IS_FINITE
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GLOBAL_OH_MOD,      ONLY : GET_GLOBAL_OH, OH
      USE TIME_MOD,           ONLY : GET_DAY,   GET_MONTH
      USE TIME_MOD,           ONLY : GET_NYMDb, GET_NYMDe
      USE TIME_MOD,           ONLY : GET_TAU,   GET_YEAR
      USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
      USE TRANSFER_MOD,       ONLY : TRANSFER_2D
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass 
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!                                                                             .
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere 
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and 
!        destruction in strat (CH4_strat.f).
!
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CHEMCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!  (4 ) LD43 is already declared in CMN_DIAG; don't redefine it (bmy, 11/15/01)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Now reference AD from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f"  Now make FIRSTCHEM a local SAVEd variable.  Now 
!        reference ALBD from "dao_mod.f".  Now use MONTH and JDATE from "CMN"
!        instead of LMN and LDY. (bmy, 11/15/02)
!  (7 ) Remove NYMDb, NYMDe from the arg list.  Now use functions GET_MONTH,
!        GET_NYMDb, GET_NYMDe, GET_MONTH, GET_DAY from the new "time_mod.f"
!        (bmy, 3/27/03) 
!  (8 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Remove reference to BPCH2_MOD, it's not needed (bmy, 10/3/05)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER      :: WTAIR = 28.966d0
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, K, M, N
      INTEGER                :: IJ, JJ, NPART, III, JJJ
      INTEGER                :: NOHDO
      INTEGER, SAVE          :: NTALDT

      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4)       :: CYEAR
      REAL*4                 :: ARRAY(IIPAR,JJPAR,LLPAR)
      INTEGER                :: TROPP
      REAL*8                 :: XTAU            	      
      INTEGER                :: LMN
      REAL*8                 :: PREVCH4(IIPAR, JJPAR, LLPAR)

      ! Number of days per month
      INTEGER                :: NODAYS(12) = (/ 31, 28, 31, 30, 
     &                                          31, 30, 31, 31, 
     &                                          30, 31, 30, 31 /)

      ! For fields from Input_Opt
      INTEGER                :: N_TRACERS
      LOGICAL                :: LSPLIT
      LOGICAL                :: LCH4BUD

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER        :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:
!
      REAL*8 , EXTERNAL      :: BOXVL  

      !=================================================================
      ! CHEMCH4 begins here!
      !=================================================================
      WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4! ---'

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      LSPLIT    = Input_Opt%LSPLIT
      LCH4BUD   = Input_Opt%LCH4BUD

      !=================================================================
      ! (0) Calculate each box's air density [molec/cm3]
      !        do this for saving mean OH concentrations (kjw, 6/12/09)
      !=================================================================

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BAIRDENS(I,J,L) = State_Met%AD(I,J,L)    * 1000d0   /
     &                     BOXVL(I,J,L,State_Met) * 6.023D23 / WTAIR
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! (1) If the first time step ...
      !        Read monthly average ch4 loss frequencies for strat chem
      !     Read CH4LOSS at beginning of each month
      !        kjw, 2/3/2014
      !=================================================================
      IF ( FIRSTCHEM ) THEN

         ! Counter for total number of timesteps per month for CO budget.
         NTALDT = 1

         ! Now read CH4 loss frequencies instead of CO production
         ! (kjw, 11/19/11)
         ! Zero CO Production array
         !COPROD(:,:,:) = 0d0
         !print*,'READ_COPROD'
         ! Read zonally-averaged CO production [v/v/s]
         !CALL READ_COPROD
         !print*,'READ_COPROD DONE'
         CH4LOSS(:,:,:,:) = 0d0
         CALL READ_CH4LOSS

#if defined( USE_CH4_BUDGET_DIAG )
         ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
         ! (mpayer, 12/6/12)
         STT => State_Chm%Tracers
         
         ! Initialize the CH4 burden TCH4 [molec/box], (ccc, 7/23/09)
         TCH4(:,:,:,1) = STT(:,:,:,1) * XNUMOL_CH4 
         
         ! Free pointer
         NULLIFY( STT )
#endif
      ENDIF

      ! Initialize current month
      LMN = GET_MONTH()

      ! Increment counter of timesteps
      NTALDT = NTALDT + 1


      !================================================================
      ! (2) get parameterized OH fields or monthly mean fields.
      !
      ! Variables of note:
      ! ---------------------------------------------------------------
      ! (1) BOH = storage array for OH fields.
      !
      ! (2) NOHDO = switch
      !       ONLY CASE 1 available as of 7/20/10 (kjw)
      !       = 1 : Get GEOS-Chem OH (v5-07-08) (kjw, 5/28/08)
      !
      !================================================================


      ! Change value of NOHDO as listed above
      NOHDO = 1

      SELECT CASE ( NOHDO )

         ! NOHDO = 1: GEOS-Chem OH v5-07-08
         CASE ( 1 )

            ! If first of month, read monthly mean OH
            IF ( FIRSTCHEM ) THEN

               ! 3D OH Field
               BOH(:,:,:,:) = 0d0               

               ! Loop over each month, reading OH
               DO M=1,12

                  ! Get current month
                  !LMN = GET_MONTH()

                  CALL GET_GLOBAL_OH( am_I_Root, Input_Opt, 
     &                                State_Met, M,       RC )

                  ! Assign to module variable BOH
                  BOH(:,:,:,M) = OH(:,:,:)

               ENDDO

!               ! Loop over each month, reading OH
!               DO M=1,12
!
!                  ! Global OH 
!!------------------------------------------------------------------------------
!! Prior to 10/23/13:
!!                  CALL GET_GLOBAL_OH( M, State_Met )
!!------------------------------------------------------------------------------
!
!
!
!               ENDDO

            ENDIF

         ! Default.  GEOS_CHEM_STOP
         CASE DEFAULT

            WRITE( 6, '(a)' ) 'Invalid selection for NOHDO!'
            WRITE( 6, '(a)' ) 'Halting execution in CHEMCH4!'
            CALL GEOS_CHEM_STOP
            
      END SELECT

      !=================================================================
      ! (3.1) ND43 diagnostics...save [OH] in molecules/cm3
      !=================================================================

      IF ( ND43 > 0 ) THEN
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( L < LPAUSE(I,J) ) THEN
               AD43(I,J,L,1) = AD43(I,J,L,1) + BOH(I,J,L,LMN)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! (4) Save OH concentrations for printing of global mean [OH] and
      !     CH3CCLl3 at end of simulation.
      !=================================================================
      CALL CH4_OHSAVE( State_Met, State_Chm )

      !=================================================================
      ! (5) If multi-CH4 tracers, we store the CH4 total conc. to
      !     distribute the sink after the chemistry. (ccc, 2/10/09)
      !================================================================= 
      IF ( LSPLIT ) THEN

         ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
         ! (mpayer, 12/6/12)
         STT => State_Chm%Tracers

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            PREVCH4(I,J,L) = STT(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Free pointer
         NULLIFY( STT )
      ENDIF

      !=================================================================
      ! (6) calculate rate of decay of CH4 by OH oxidation.
      !=================================================================
      CALL CH4_DECAY( State_Met, State_Chm )

      !=================================================================
      ! (7) calculate CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT( Input_Opt, State_Met, State_Chm )

      !=================================================================
      ! (8) distribute the chemistry sink from total CH4 to other CH4 
      !     tracers. (ccc, 2/10/09)
      !=================================================================
      IF ( LSPLIT ) THEN
         CALL CH4_DISTRIB( PREVCH4, Input_Opt, State_Chm )
      ENDIF

#if defined( USE_CH4_BUDGET_DIAG )
      !=================================================================
      ! (9) write budget (i.e., monthly average fields).
      !
      ! Check to make sure the start and end times are on the
      ! first of a month.  If not the SR CO_budget will not
      ! work properly!
      !=================================================================
      NPART = GET_NYMDb() / 100 

      IF ( LCH4BUD .and. ( GET_NYMDb() - NPART*100 ) /= 1 ) THEN
         print*,'Start date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF


      NPART = GET_NYMDe() /100 

      IF ( LCH4BUD .and. ( GET_NYMDe() - NPART*100 ) /= 1 ) THEN      
         print*,'End date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Call CH4_BUDGET on the last day of the month
      IF ( LCH4BUD .and. GET_DAY() == NODAYS( GET_MONTH() ) ) THEN
         print*,'CALL CH4_BUDGET'
         CALL FLUSH ( 6 )
     
         ! Disable CH4 budget for SEAC4RS code to save memory
         !   kjw, 2/3/2014
         CALL CH4_BUDGET( State_Chm )

	 NTALDT  = 0
         print*,'END CH4_BUDGET'
         call flush(6)
     
      ENDIF
#endif

      print*,'END CHEMCH4'
      call flush(6)

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      END SUBROUTINE CHEMCH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_coprod
!
! !DESCRIPTION: Subroutine READ\_COPROD reads production and destruction rates
!  for CO in the stratosphere. CO destruction rate is assumed equal to CH4
!  production rate for the GEOS-Chem CH4 simulation. (bnd, bmy, 1/17/01,
!  10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_COPROD
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE TRANSFER_MOD,  ONLY : TRANSFER_ZONAL
        
      USE CMN_SIZE_MOD        ! Size parameters

      IMPLICIT NONE
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_COPROD is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JJPAR,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Treat MERRA in the same way as for GEOS-5 (bmy, 8/13/10)
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, L, M
      REAL*4                 :: ARRAY(1,JJPAR,LGLOB)
      REAL*4                 :: DUMMY_IN(JJPAR,LGLOB)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME
      REAL*8                 :: DUMMY_OUT(JJPAR,LGLOB)

      !=================================================================
      ! READ_COPROD begins here!
      ! 
      ! Read P(CO) for all 12 months
      !=================================================================

      DO M = 1, 12 

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Construct filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/'  //
     &              'COprod.GEOS5.05x0666_NA.trim'
#else
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/' //
     &              'COprod.'        // GET_NAME_EXT()    //
     &              '.'              // GET_RES_EXT()
#endif

         WRITE( 6, 93 ) TRIM ( FILENAME )
 93      FORMAT( '     - READ_COPROD: Reading COprod: ', a )
	 CALL FLUSH( 6 )


         ! Read P(CO) in units of [v/v/s].  drevet.
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,     
     &                    XTAU,      1,         JJPAR,     
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% SCHEM PATCH: Do not call TRANSFER_ZONAL when using the GEOS-5 grid
!%%% since the SCHEM data are saved at the reduced resolution.  This will
!%%% be fixed in a subsequent version. (ltm, bmy, 6/2/10)
!%%%-----------------------------------------------------------------------
!%%% Original code here:
!%%%         ! use 2D arrays for TRANSFER ZONAL
!%%%         DUMMY_IN(:,:) = ARRAY(1,:,:)
!%%%
!%%%         ! Copy REAL*4 to REAL*8 data, and resize from (JJPAR,LGLOB) 
!%%%         ! to (JJPAR,LLPAR) -- vertically regrid if necessary
!%%%         CALL TRANSFER_ZONAL( DUMMY_IN, DUMMY_OUT )
!%%%
!%%%         COPROD(:,:,M) = DUMMY_OUT(:,:)
!%%%-----------------------------------------------------------------------
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
         ! Cast from REAL*4 to REAL*8
         COPROD(:,:,M) = ARRAY(1,:,:)
#else
         ! use 2D arrays for TRANSFER ZONAL
         DUMMY_IN(:,:) = ARRAY(1,:,:)

         ! Copy REAL*4 to REAL*8 data, and resize from (JJPAR,LGLOB) 
         ! to (JJPAR,LLPAR) -- vertically regrid if necessary
         CALL TRANSFER_ZONAL( DUMMY_IN, DUMMY_OUT )

         COPROD(:,:,M) = DUMMY_OUT(:,:)
#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ENDDO

      END SUBROUTINE READ_COPROD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_ch4loss
!
! !DESCRIPTION: Subroutine READ\_CH4LOSS reads CH4 loss frequencies in the
!  stratosphere. These values constitute a linearized stratospheric CH4
!  chemistry scheme. Loss frequencies from 4x5 degree output from the GMI 
!  model. Thanks to Lee Murray for the ch4 loss frequencies. (kjw, 11/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CH4LOSS
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE TIME_MOD,      ONLY : GET_MONTH

      USE CMN_SIZE_MOD        ! Size parameters

      IMPLICIT NONE
! 
! !REVISION HISTORY:!
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_CH4LOSS is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JJPAR,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Treat MERRA in the same way as for GEOS-5 (bmy, 8/13/10)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  12 Feb 2014 - K. Wecht    - Add modifications for 0.25 x 0.3125 NA grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, L, M
      REAL*4                 :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME


      !=================================================================
      ! READ_CH4LOSS begins here!
      ! 
      ! Read P(CO) for all 12 months
      !=================================================================

      ! Construct filename
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR )         // 'CH4_201203/' //
     &           'gmi.ch4loss.geos5_47L.' // GET_RES_EXT() //
     &           '_NA.bpch'
#elif defined( NESTED_NA ) && defined( GRID025x03125 ) 
      FILENAME = TRIM( DATA_DIR )         // 'CH4_201402/' //
     &           'gmi.ch4loss.geos5_47L.'                  //
     &           GET_RES_EXT()            //  '_NA.bpch'
#else
      FILENAME = TRIM( DATA_DIR )         // 'CH4_201203/' //
     &           'gmi.ch4loss.geos5_47L.' // GET_RES_EXT() // 
     &           '.bpch'
#endif


      WRITE( 6, 93 ) TRIM ( FILENAME )
 93   FORMAT( '     - READ_CH4LOSS: Reading Ch4loss: ', a )
      CALL FLUSH( 6 )

      ! Read data for each month
      DO M = 1, 12 

      !M = GET_MONTH()

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Read Loss frequencies in units of [1/s].  drevet.
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,     
     &                    XTAU,      IIPAR,     JJPAR,     
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. )

         ! Place array into CH4LOSS module variable
         CH4LOSS(:,:,:,M) = ARRAY(:,:,:)

      ENDDO

      END SUBROUTINE READ_CH4LOSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_decay
!
! !DESCRIPTION: Subroutine CH4\_DECAY calculates the decay rate of CH4 by OH.
!  OH is the only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_DECAY( State_Met, State_Chm )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : AD19   	      
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM, ITS_A_NEW_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH
      USE GIGC_State_Chm_Mod, ONLY : ChmState

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND19
      USE CMN_MOD                  ! LPAUSE
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!                                                                             .
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!                                                                             .
!  We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
!  (bmy, 4/18/00)  
!                                                                             .
!  Monthly loss of CH4 is summed in TCH4(3)
!     TCH4(3)  = CH4 sink by OH
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_DECAY is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL             :: FIRST_DECAY=.TRUE.
      INTEGER             :: I, J, L, M, N, MM, LMN
      REAL*8              :: DT, GCH4, STT2GCH4, KRATE

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER     :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CH4_DECAY begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM() * 60d0

      ! Current month
      LMN = GET_MONTH()

      !=================================================================
      ! Compute decay of CH4 by OH in the troposphere
      !
      ! The decay for CH4 is calculated by:
      !    OH + CH4 -> CH3 + H2O 
      !    k = 2.45E-12 exp(-1775/T)
      !
      ! This is from JPL '97.
      ! JPL '00 & '06 do not revise '97 value. (jsw, kjw)
      !=================================================================
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      MM = MAXVAL( LPAUSE )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( L, J, I, KRATE, STT2GCH4, GCH4 )
!$OMP+REDUCTION( +:TROPOCH4 )
      DO L = 1, MM !MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( L < LPAUSE(I,J) ) THEN 

            ! Use 24-hr avg temperature to calc. rate coeff.
            ! citation needed
            KRATE = 2.45d-12 * EXP( -1775d0 / State_Met%T(I,J,L) )  

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            STT2GCH4 = 1d0 / State_Met%AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4 

            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4
	
#if defined( USE_CH4_BUDGET_DIAG )
            ! Sum loss in TCH4(3) (molecules/box)
            TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
     &                      ( GCH4  * BOXVL(I,J,L,State_Met) *
     &                        KRATE * BOH(I,J,L,LMN)         * DT )
#endif

            TROPOCH4=TROPOCH4 + GCH4 * KRATE * BOH(I,J,L,LMN) * DT / 
     &               STT2GCH4

            ! Modify AD19 Diagnostic
            ! How much CH4 (kg) is lost by reaction with OH
	    IF ( ND19 > 0 ) THEN  ! --> [kg/box]
               IF ( L <= LD19 ) THEN
!	    	AD19(I,J,12) = AD19(I,J,12) + 
                  AD19(I,J,L) = AD19(I,J,L) + 
     &                 ( GCH4 * KRATE * BOH(I,J,L,LMN) * DT ) / STT2GCH4
               ENDIF
	    ENDIF

            ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt) 
            GCH4 = GCH4 * ( 1d0 - KRATE * BOH(I,J,L,LMN) * DT )
		
            ! Convert back from [molec/cm3] --> [kg/box]
            STT(I,J,L,1) = GCH4 / STT2GCH4

         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      print*,'% --- CHEMCH4: CH4_DECAY: TROP DECAY (Tg): ',TROPOCH4/1e9
      print*,'Trop decay should be over 1Tg per day globally'
      print*,'    ~ 500Tg/365d ~ 1.37/d'

      END SUBROUTINE CH4_DECAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_ohsave
!
! !DESCRIPTION: Subroutine CH4\_OHSAVE archives the CH3CCl3 lifetime from the
!  OH used in the CH4 simulation. (bnd, jsw, bmy, 1/16/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_OHSAVE( State_Met, State_Chm )
!
! !USES:
!
      ! References to F90 modules
      USE DIAG_OH_MOD,        ONLY : DO_DIAG_OH_CH4
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_MONTH

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_MOD                  ! LPAUSE
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!                                                                             .
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_OHSAVE is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now call DO_DIAG_OH_CH4 to pass OH diagnostic info to the
!        "diag_oh_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, LMN
      REAL*8               :: MASST,   AREA_CM2
      REAL*8               :: KCLO,    LOSS,       OHMASS  
      REAL*8               :: KCH4,    CH4LOSE,    CH4MASS
      REAL*8               :: CH4EMIS, CH4TROPMASS

      ! Pointers
      REAL*8, POINTER      :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8, EXTERNAL     :: BOXVL

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! (1) Pass OH mass, total air mass, and  to "diag_oh_mod.f"
      ! (2) ND59: Diagnostic for CH3CCl3 calculation
      !=================================================================

      ! Current month
      LMN = GET_MONTH()

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT => State_Chm%Tracers

      ! 1. Calculate OH mass and total air mass
      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR 
      DO I = 1, IIPAR 
         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( L < LPAUSE(I,J) ) THEN

            ! Calculate OH mass [molec / box]
            OHMASS = BOH(I,J,L,LMN)  * BAIRDENS(I,J,L) *
     &               BOXVL(I,J,L,State_Met)

            ! Calculate total air mass [molec / box]
            MASST  = BAIRDENS(I,J,L) * BOXVL(I,J,L,State_Met)

            ! Calculate CH3CCl3 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCLO   = 1.64d-12 * EXP( -1520.d0 / State_Met%T(I,J,L) )

            ! Calculate Loss term [molec / box / s]
            LOSS   = KCLO            * BOH(I,J,L,LMN)  *
     &               BAIRDENS(I,J,L) * BOXVL(I,J,L,State_Met)


            ! Calculate CH4 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCH4   = 2.45d-12 * EXP( -1775d0 / State_Met%T(I,J,L) )

            ! Calculate CH4 mass [molec / box] from [kg / box]
            CH4TROPMASS = STT(I,J,L,1) * XNUMOL_CH4 
            CH4MASS     = STT(I,J,L,1) * XNUMOL_CH4 

            ! Calculate loss term  [molec /box / s]
            CH4LOSE = KCH4            * BOH(I,J,L,LMN)  *
     &                BAIRDENS(I,J,L) * BOXVL(I,J,L,State_Met)

            ! Calculate CH4 emissions [molec / box / s]
            !   Only for surface level
            ! Grid box surface area [cm2]
            IF ( L .GT. 1 ) THEN 
               CH4EMIS = 0d0
            ELSE
               AREA_CM2 = GET_AREA_CM2( I, J, L )
               CH4EMIS  = SUM(CH4_EMIS(I,J,2:10)) + CH4_EMIS(I,J,12)
               CH4EMIS  = CH4EMIS * AREA_CM2 ! [molec/cm2/s] --> [molec/box/s]
            ENDIF

         ELSE

            OHMASS      = 0d0
            MASST       = 0d0
            LOSS        = 0d0
            CH4LOSE     = 0d0
            CH4TROPMASS = 0d0
            CH4EMIS     = 0d0
            CH4MASS     = STT(I,J,L,1) * XNUMOL_CH4 

         ENDIF

         ! Pass OH mass, total mass, and loss to "diag_oh_mod.f",
         ! which calculates mass-weighted mean [OH] and CH3CCl3
         ! lifetime.
         CALL DO_DIAG_OH_CH4( I, J, L, OHMASS, MASST, LOSS, 
     &        CH4LOSE, CH4TROPMASS, CH4EMIS, CH4MASS )

      ENDDO
      ENDDO
      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CH4_OHSAVE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_strat
!
! !DESCRIPTION: Subroutine CH4\_STRAT calculates uses production rates for CH4
!  to  calculate loss of CH4 in above the tropopause. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_STRAT( Input_Opt, State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : GET_MIN_TPAUSE_LEVEL
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,         ONLY : CHECK_STT
      USE TRACERID_MOD,       ONLY : IDTCH4

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_MOD                  ! LPAUSE
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!  Only production by CH4 + OH is considered.
!                                                                             .
!  The annual mean tropopause is stored in the LPAUSE array 
!  (from header file "CMN").  LPAUSE is defined such that: 
!                                                                             .
!  Levels           1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/18/00)
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use 
!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER   :: WTAIR = 28.966d0
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, LMN, BASELEV
      REAL*8              :: DT, GCH4, STT2GCH4, LRATE
      CHARACTER*20        :: STT_TEST
      CHARACTER*20        :: STT2GCH4_CHAR
      LOGICAL,SAVE        :: FIRSTCHEM=.TRUE.

      ! Local variables for quantities from Input_Opt
      LOGICAL             :: LUCX
      LOGICAL             :: ITS_A_FULLCHEM_SIM

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER     :: STT(:,:,:,:)
!
! !EXTERNAL FUNCTIONS:   
!
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CH4_STRAT begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LUCX                = Input_Opt%LUCX
      ITS_A_FULLCHEM_SIM  = Input_Opt%ITS_A_FULLCHEM_SIM

      ! If unified chemistry is active, ignore all of this
      IF ( .not. LUCX ) THEN

         ! Chemistry timestep [s]
         DT  = GET_TS_CHEM() * 60d0

         ! Current month
         LMN = GET_MONTH()

         IF ( ITS_A_FULLCHEM_SIM ) THEN
            
            !=============================================================
            ! (0) Calculate each box's air density [molec/cm3]
            !=============================================================

            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               BAIRDENS(I,J,L) = State_Met%AD(I,J,L)    * 1d3      / 
     &                           BOXVL(I,J,L,State_Met) * 6.023d23 /
     &                           WTAIR
            ENDDO
            ENDDO
            ENDDO

            !=============================================================
            ! (1) If the first time step...
            !=============================================================
            IF ( FIRSTCHEM ) THEN
               ! Zero CO production array
               COPROD(:,:,:) = 0d0
               print*,'READ_COPROD'
               ! Read zonally-averaged CO production [v/v/s]
               CALL READ_COPROD
               print*,'READ_COPROD DONE'
            ENDIF
         ENDIF

         ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
         ! (mpayer, 12/6/12)
         STT => State_Chm%Tracers

         !=================================================================
         ! Loop over stratospheric boxes only
         !=================================================================
         BASELEV = GET_MIN_TPAUSE_LEVEL()
         DO L = BASELEV, LLPAR 
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( ITS_IN_THE_NOCHEMGRID( I, J,  L, State_Met) ) THEN

               ! Conversion factor [kg/box] --> [molec/cm3]
               ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
               STT2GCH4 = 1d0 / State_Met%AIRVOL(I,J,L) / 1d6 *
     &                    XNUMOL_CH4

               ! CH4 in [molec/cm3]
               GCH4 = STT(I,J,L,1) * STT2GCH4

               ! Loss rate [molec/cm3/s]
               LRATE = GCH4 * CH4LOSS( I,J,L,LMN )

               ! Update Methane concentration in this grid box [molec/cm3]
               GCH4 = GCH4 - ( LRATE * DT )

! Update stratospheric chem to use linearized CH4 loss frequencies
!  (kjw, 11/19/11)
!               ! Sum loss in TCH4(3) [molec CH4/box] in the stratosphere
!               ! [molec/cm3] * [v/v/s] * [s] * [cm3/box] = [molec CH4/box]
!               TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
!        &                      ( BAIRDENS(I,J,L) * COPROD(J,L,LMN) *
!        &                        DT              * BOXVL(I,J,L)    )
!
!               ! Calculate new CH4 value [molec CH4/cm3] in the stratosphere
!               ! [v/v/s] * [s] * [molec/cm3] = [molec CH4/cm3] 
!               GCH4 = GCH4 - ( COPROD(J,L,LMN) * DT * BAIRDENS(I,J,L) )

               ! Convert back from [molec CH4/cm3] --> [kg/box] 
               STT(I,J,L,1) = GCH4 / STT2GCH4

! With new linearized chemistry, STT should never be negative
!  (kjw, 11/19/11)
!               ! Prevent negative values in STT
!              IF ( STT(I,J,L,1) < 0 ) THEN
!          	STT(I,J,L,1)=0
!              ENDIF

            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDIF ! LUCX (SDE 03/25/13)

      FIRSTCHEM = .FALSE.

      ! Check STT for negatives and NaN
      CALL CHECK_STT( State_Chm, 'IN global_ch4_mod.F: CH4_STRAT' )

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CH4_STRAT
!EOC
#if defined( USE_CH4_BUDGET_DIAG )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_budget
!
! !DESCRIPTION: Subroutine CH4\_BUDGET calculates the budget of CH4. This SR
!  only works for monthly averages, so be sure to start on the first of the
!  month and run to another first of the month! (jsw, bnd, bmy, 1/16/01,
!  10/3/05)
!
!  Disable CH4 Budget for SEAC4RS code to save memory
!     kjw, 2/3/2014
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_BUDGET( State_Chm )
!
! !USES:
!
      USE BPCH2_MOD,  ONLY : BPCH2,       BPCH2_HDR,   GET_MODELNAME
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : GET_MONTH,   GET_YEAR
      USE TIME_MOD,   ONLY : GET_DIAGb,   GET_CT_DYN
      USE GIGC_State_Chm_Mod, ONLY : ChmState

      USE CMN_SIZE_MOD     ! Size parameters
!
! !INPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = Tropospheric CH4 sink by OH
!  SOURCES
!           ( 4) = Total Sources
!           ( 5) = Industrial  ( Gas + Oil + Mine )
!           ( 6) = Agriculture ( Enteric fermentation + Manure + Rice + Waste
!                              + Waste water )
!           ( 7) = Biomass burning
!           ( 8) = Termites 
!           ( 9) = Wetland
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/13/01)
!  (4 ) Renamed XLABEL to LABEL so as not to conflict w/ "CMN"
!  (5 ) Now use functions GET_MONTH, GET_YEAR, GET_DIAGb, and GET_CT_DYN from 
!        "time_mod.f".  Removed LMN from the arg list and made it a local 
!        variable.  Use functions GET_XOFFSET and GET_YOFFSET from 
!        "grid_mod.f".  (bmy, 3/27/03)
!  (6 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Modified for the run with new emissions (j drevet, 03/06)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, K, L, M, NERROR, UD, LMN

      REAL*8              :: STTCONV, TGS, SCALEDYN
      REAL*8              :: NTP, NTQ, NTP2, NTQ2 
      REAL*8              :: SOURCES, SINKS

      CHARACTER(LEN=17)   :: MERGE
      CHARACTER(LEN=13)   :: MERGE2

      ! For binary punch file, v. 2.0
      REAL*4              :: ARRAY(IIPAR, JJPAR, LLPAR)
      REAL*4              :: LONRES, LATRES

      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER, PARAMETER  :: HALFPOLAR = 1
      INTEGER, PARAMETER  :: CENTER180 = 1

      CHARACTER (LEN=20)  :: MODELNAME 
      CHARACTER (LEN=40)  :: UNIT
      CHARACTER (LEN=40)  :: RESERVED = ''
      CHARACTER (LEN=40)  :: CATEGORY 
      CHARACTER (LEN=80)  :: LABEL

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CH4_BUDGET begins here!
      !
      ! Initialize quantities 
      !=================================================================
      IFIRST    = GET_XOFFSET() + 1
      JFIRST    = GET_YOFFSET() + 1
      LFIRST    = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Current month
      LMN       = GET_MONTH()

      ! Make up a category name for GAMAP (use 8 characters)
      CATEGORY  = 'CH4BUDT'

      ! Get the proper model name for the binary punch file
      MODELNAME = GET_MODELNAME()

      ! Descriptor string
      LABEL    = 'GEOS-CHEM -- CH4 Budget output (jsw, bmy, 1/16/01)'

      ! Unit of quantity being saved
      UNIT      = 'Tg'  !(NOTE: check w/ bnd to get the right units!!!)

      ! Scale factor for dynamic time steps
      SCALEDYN  = FLOAT( GET_CT_DYN() ) + 1D-20

      !=================================================================
      ! Store the final burden of CH4 in TCH4(2) 
      ! Convert kg CH4/box to molecules/box.
      !=================================================================
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      TCH4(:,:,:,2) = 0d0
      TCH4(:,:,:,2) = STT(:,:,:,1) * XNUMOL_CH4

      ! Free pointer
      NULLIFY( STT )

      !=================================================================
      ! Write GLOBAL AVERAGES for all layers to ASCII file
      !=================================================================
      WRITE( MERGE, 2 ) GET_MONTH(), GET_YEAR()
 2    FORMAT( 'CH4budget.', I2.2, '.',I4 )

      OPEN( 189, FILE=MERGE, STATUS='UNKNOWN' )
      REWIND( 189 )
      
      TGS     = 1.D-9
      STTCONV = XNUMOL_CH4/TGS
      SOURCES = 0.D0
      SINKS   = 0.D0
      NERROR  = 0
      
      WRITE(189,18)
      WRITE(189,1801)
 1801 FORMAT('*************************')
      WRITE(189,1800)
 1800 FORMAT('LAYERS 1 - 20')
      WRITE(189,1801)
      WRITE(189,18)

      WRITE(189,18)
      WRITE(189,38)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)
 1990 FORMAT('Tropospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      WRITE(189,20)NTP,NTP/STTCONV

      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,21)NTP2,NTP2/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
 1991 FORMAT('Stratospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks   jsw has checked correctness of code for sinks.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTP+NTQ

      WRITE(189,22) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      WRITE(189,220) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,29) 
      WRITE(189,34) SINKS,SINKS/STTCONV  !Just OH sink 
      WRITE(189,18)
      WRITE(189,30)

      ! Sources
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,4,4,1)
      SOURCES=NTP

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      !jsw Following lines added by jsw.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,10,10,1)
      WRITE(189,35) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      SINKS=SINKS-NTP  !Minus sign because soil absorption is negative.

      WRITE(189,29) 
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS,
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
      
      !=================================================================
      ! Write SOUTHERN HEMISPHERE averages to ASCII file
      ! jsw:  I have not modified the remaining code for CH4.
      !================================================================= 

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,36)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF( NTP > 0d0) SINKS = SINKS + NTP

      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF( NTP2 > 0d0 ) SINKS = SINKS + NTP2      

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
 
      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF

      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,34) SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)
      
      ! Sources
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,5,9,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF( NTP2 > 0d0 ) SINKS = SINKS + NTP2      

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV

      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

      !=================================================================
      ! Write NORTHERN HEMISPHERE averages to ASCII file 
      ! jsw:  I have not modified the remaining code for CH4.
      !================================================================= 

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,37)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,1991)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF

      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,34)SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)

      ! Sources
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,5,9,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,1)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

 18   FORMAT()
 19   FORMAT('                    #Molecules               TG')
 20   FORMAT('  Start of Month  :',E10.3,10x,F10.3)
 21   FORMAT('  End of Month    :',E10.3,10x,F10.3)
 22   FORMAT('  CH4 decay-trop   :',E10.3,2x,F6.1,2x,F10.3)
 220  FORMAT('  CH4 decay-strat  :',E10.3,2x,F6.1,2x,F10.3)
 24   FORMAT('  Industrial      :',E10.3,2x,F6.1,2x,F10.3)
 25   FORMAT('  Biomass Burning :',E10.3,2x,F6.1,2x,F10.3)
 26   FORMAT('  Termites        :',E10.3,2x,F6.1,2x,F10.3)
 27   FORMAT('  Wetland         :',E10.3,2x,F6.1,2x,F10.3)
 270  FORMAT('  N-S Ex.-trop    :',E10.3,2x,F6.1,2x,F10.3)
 2700 FORMAT('  N-S Ex.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 28   FORMAT('Total Sources     :',E10.3,10x,F10.3)
 288  FORMAT('Initial-Final+Sources-Sinks=',E10.3,2x,F10.3)
 289  FORMAT('Net Gain          : ',E10.3,10x,F10.3)
 29   FORMAT('                     ---------')
 30   FORMAT('SOURCES                          %Source')
 31   FORMAT('SINKS                            %Sink')
 34   FORMAT('Total Sinks       :',E10.3,10x,F10.3)
 35   FORMAT('  Soil absorption :',E10.3,2x,F6.1,2x,F10.3)
 39   FORMAT('  Agriculture     :',E10.3,2x,F6.1,2x,F10.3)
      
 36   FORMAT('*****  Southern Hemisphere  *****')
 37   FORMAT('*****  Northern Hemisphere  *****')
 38   FORMAT('*****  Global  *****')
      
      CLOSE(189)

!     !=================================================================
!     ! Also save to binary punch file. Don't save the bpunch file 
!     ! anymore, because it's not used. Keep the code for reference.
!     ! The code creates the bpunch file fort.190. Should use a diag.
!     ! instead. (ccc, 8/14/09)
!     !=================================================================
!     CALL BPCH2_HDR( 190, LABEL )
!
!     DO K = 1, N_CH4
!       
!        ! Cast REAL*8 into REAL*4, convert from molec to Tg
!        ARRAY(:,:,:) = TCH4(:,:,:,K) / STTCONV
!       
!        ! Save the data block 
!        CALL BPCH2( 190,       MODELNAME,   LONRES,      LATRES,
!    &               HALFPOLAR, CENTER180,   CATEGORY,    K,     
!    &               UNIT,      GET_DIAGB(), GET_DIAGb(), RESERVED, 
!    &               IIPAR,     JJPAR,       LLPAR,       IFIRST,  
!    &               JFIRST,    LFIRST,      ARRAY )
!     ENDDO
!
!     CLOSE(190)

      !=================================================================
      ! Final burden at last of month equals initial burden
      ! of next month.  Also set TCH4 = 0 for next month.
      !=================================================================
      TCH4(:,:,:,1      ) = TCH4(:,:,:,2)
      TCH4(:,:,:,2:N_CH4) = 0d0
	
      END SUBROUTINE CH4_BUDGET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sum_ch4
!
! !DESCRIPTION: Function SUM\_CH4 sums a section of the TCH4 array bounded by
!  the input variables I1, I2, J1, J2, L1, L2, K1, K2.  SUM\_CH4 is called by
!  module subroutine CH4\_BUDGET. (jsw, bnd, bmy, 1/16/01)
!
!  Disable CH4 Budget for SEAC4RS code to save memory
!     kjw, 2/3/2014
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION SUM_CH4( I1, I2, J1, J2, L1, L2, K1, K2, UPDOWN )
!
! !USES:
!
      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_MOD          ! LPAUSE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I1, I2 ! Min/max longitude indices of TCH4 to sum
      INTEGER, INTENT(IN) :: J1, J2 ! Min/max latitude  indices of TCH4 to sum 
      INTEGER, INTENT(IN) :: L1, L2 ! Min/max altitude  indices of TCH4 to sum
      INTEGER, INTENT(IN) :: K1, K2 ! Min/max tracer    indices of TCH4 to sum
      INTEGER, INTENT(IN) :: UPDOWN ! Sum in trop (=1) or in strat (=0)
!
! !REMARKS:
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = Tropospheric CH4 sink by OH
!  SOURCES
!           ( 4) = Total Source
!           ( 5) = Industral
!           ( 6) = Agriculture
!           ( 7) = Biomass Burning
!           ( 8) = Termites 
!           ( 9) = Wetland
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!           (12) = ...
!                                                                             .
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/17/00)
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, K, L
      INTEGER             :: LPAUSE_MIN, LPAUSE_MAX

      !=================================================================
      ! SUM_CH4 begins here!
      !=================================================================

      ! Compute the minimum value of LPAUSE once for use in
      ! the DO-loops below (bmy, 4/18/00)
      LPAUSE_MIN = MINVAL( LPAUSE )
      LPAUSE_MAX = MAXVAL( LPAUSE )

      !### Debug
      !print*,'LPAUSE MIN/MAX=',LPAUSE_MIN,LPAUSE_MAX  
      !print*,'L1,L2=',L1,L2
      
      ! Initialize SUM_CH4
      SUM_CH4 = 0d0

      ! Test on UPDOWN
      IF ( UPDOWN == 1 ) THEN

         !=============================================================
         ! UPDOWN = 1: Sum up from the surface to the tropopause
         !=============================================================
         DO K = K1, K2
         DO L = L1, LPAUSE_MAX
         DO J = J1, J2
         DO I = I1, I2
            IF ( L < LPAUSE(I,J) ) THEN 
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ELSE

         !=============================================================
         ! UPDOWN = 0: Sum up from the tropopause to the atm top
         !=============================================================
         DO K = K1,         K2
         DO L = LPAUSE_MIN, L2
         DO J = J1,         J2
         DO I = I1,         I2
            IF ( L >= LPAUSE(I,J) ) THEN 
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF            
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      
      END FUNCTION SUM_CH4
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_distrib
!
! !DESCRIPTION: Subroutine CH4\_DISTRIB allocates the chemistry sink to
!  different emission tracers. (ccc, 10/2/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CH4_DISTRIB( PREVCH4, Input_Opt, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE ERROR_MOD,          ONLY : SAFE_DIV
   
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      REAL*8                        :: PREVCH4(IIPAR,JJPAR,LLPAR)! CH4 bef chem
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm args
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER         :: N_TRACERS
      INTEGER         :: N, I, J, L

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !========================================================================
      ! CH4_DISTRIB begins here
      !========================================================================

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      DO N=2,N_TRACERS

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            STT(I,J,L,N) = SAFE_DIV(STT(I,J,L,N),PREVCH4(I,J,L),0.d0)  
     &                     * STT(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CH4_DISTRIB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ch4_anthro
!
! !DESCRIPTION: Function GET\_CH4\_ANTHRO returns the monthly average CH4
!  emissions at GEOS-Chem grid box (I,J).  Data will be returned in units
!  of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CH4_ANTHRO( I, J, N, Input_Opt ) RESULT( CH4_ANTHRO )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD,       ONLY : IDTCH4
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN) :: J           ! GEOS-Chem latitude index
      INTEGER,        INTENT(IN) :: N           ! GEOS-Chem tracer index
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options!
!
! !RETURN VALUE:
!
      REAL*8 :: CH4_ANTHRO   ! Anthro CH4 emissions [molec/cm2/s]
! 
! !REVISION HISTORY: 
!  24 Jan 2012 - M. Payer    - Initial version adapted from GET_RETRO_ANTHRO
!  25 Mar 2013 - S.D. Eastham- Adapted for unified simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      ! Local variables for quantities from Input_Opt
      LOGICAL :: LCH4EMIS

      !=================================================================
      !  GET_CH4_ANTHRO begins here
      !=================================================================

      ! Copy fields from INPUT_OPT
      LCH4EMIS = Input_Opt%LCH4EMIS

      IF (( LCH4EMIS ).and.( N == IDTCH4 )) THEN
         CH4_ANTHRO = SUMANTHRO(I,J)
      ELSE
         CH4_ANTHRO = -1d0
      ENDIF

      END FUNCTION GET_CH4_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ch4_biog
!
! !DESCRIPTION: Function GET\_CH4\_BIOG returns annual average biogenic
!  CH4 emissions at GEOS-Chem grid-box (I,J). Data will be returned in
!  units of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CH4_BIOG( I, J, N, Input_Opt ) RESULT( CH4_BIOG )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD,       ONLY : IDTCH4
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: I           ! GEOS-Chem longitude index
      INTEGER,        INTENT(IN) :: J           ! GEOS-Chem latitude index
      INTEGER,        INTENT(IN) :: N           ! GEOS-Chem tracer index
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options!
!
! !RETURN VALUE:
!
      REAL*8  :: CH4_BIOG   ! Biogenic CH4 emissions [molec/cm2/s]
! 
! !REVISION HISTORY: 
!  24 Jan 2012 - M. Payer    - Initial version adapted from GET_RETRO_ANTHRO
!  25 Mar 2013 - S.D. Eastham- Adapted for unified simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      ! Local variables for quantities from Input_Opt
      LOGICAL :: LCH4EMIS

      !=================================================================
      !  GET_CH4_BIOG begins here
      !=================================================================

      ! Copy fields from INPUT_OPT
      LCH4EMIS = Input_Opt%LCH4EMIS

      IF (( LCH4EMIS ).and.( N == IDTCH4 )) THEN
         CH4_BIOG = SUMBIOG(I,J)
      ELSE
         CH4_BIOG = -1d0
      ENDIF

      END FUNCTION GET_CH4_BIOG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_ch4
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_CH4 allocates and zeroes module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_CH4( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CMN_DIAG_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is called from GIGC_INIT_EXTRA (in GeosCore/input_mod.f)
! 
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used
!  11 Apr 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: AS

      !=================================================================
      ! INIT_GLOBAL_CH4 begins here!
      !=================================================================

      ! Assume Success
      RC = GIGC_SUCCESS

      ALLOCATE( AVGOH( NSEAS, NCMSLATS, NCMSALTS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVGOH' )
      AVGOH = 0d0

      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BAIRDENS' )
      BAIRDENS = 0d0

      ALLOCATE( BOH( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BOH' )
      BOH = 0d0

      ALLOCATE( COPROD( JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
      COPROD = 0d0

      ALLOCATE( CH4LOSS( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4LOSS' )
      CH4LOSS = 0d0

      ALLOCATE( PAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PAVG' )
      PAVG = 0d0

      ALLOCATE( TAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAVG' )
      TAVG = 0d0

#if defined( USE_CH4_BUDGET_DIAG )
      ALLOCATE( TCH4( IIPAR, JJPAR, LLPAR, N_CH4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCH4' )
      TCH4 = 0d0      
#endif

      ALLOCATE( CH4_EMIS( IIPAR, JJPAR, PD58), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_EMIS' )
      CH4_EMIS = 0d0

      ALLOCATE( SUMANTHRO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMANTHRO' )
      SUMANTHRO = 0d0

      ALLOCATE( SUMBIOG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMBIOG' )
      SUMBIOG = 0d0

      ! Initialize tropoch4 (counts total decay of CH4 due to OH)
      TROPOCH4 = 0d0

      END SUBROUTINE INIT_GLOBAL_CH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_ch4
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_CH4
! 
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  12 Feb 2014 - K. Wecht    - Disable CH4 budget diagnostic (bracket the 
!                              code out with #ifdef blocks so it can be used)
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_GLOBAL_CH4 begins here!
      !=================================================================

      IF ( ALLOCATED( BAIRDENS  ) ) DEALLOCATE( BAIRDENS  )
      IF ( ALLOCATED( BOH       ) ) DEALLOCATE( BOH       )
      IF ( ALLOCATED( COPROD    ) ) DEALLOCATE( COPROD    )
      IF ( ALLOCATED( CH4LOSS   ) ) DEALLOCATE( CH4LOSS   )
#if defined( USE_CH4_BUDGET_DIAG )
      IF ( ALLOCATED( TCH4      ) ) DEALLOCATE( TCH4      )
#endif
      IF ( ALLOCATED( CH4_EMIS  ) ) DEALLOCATE( CH4_EMIS  )
      IF ( ALLOCATED( SUMANTHRO ) ) DEALLOCATE( SUMANTHRO )
      IF ( ALLOCATED( SUMBIOG   ) ) DEALLOCATE( SUMBIOG   )

      END SUBROUTINE CLEANUP_GLOBAL_CH4
!EOC
      END MODULE GLOBAL_CH4_MOD
