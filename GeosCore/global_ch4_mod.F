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
      USE HCO_ERROR_MOD    ! For HEMCO error reporting
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

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
      PUBLIC :: EMISSCH4
      PUBLIC :: CHEMCH4
      PUBLIC :: INIT_GLOBAL_CH4
      PUBLIC :: CLEANUP_GLOBAL_CH4
!
! !PUBLIC DATA MEMBERS:
!
      REAL(fp), PARAMETER,   PUBLIC :: XNUMOL_CH4 = 6.0221d23 / 16d-3

#if defined( USE_CH4_BUDGET_DIAG )
      REAL(fp), ALLOCATABLE, PUBLIC :: TCH4(:,:,:,:)
#endif

      ! Make CH4_EMIS now public so that it can be used by vdiff_mod.F90
      REAL(fp),  ALLOCATABLE, PUBLIC :: CH4_EMIS(:,:,:)
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
!  09 Sep 2014 - C. Keller   - Implemented HEMCO for emissions. 
!  16 Sep 2014 - C. Keller   - Now also get chemical fields (OH and CH4 prod)
!                              through HEMCO. Removed READ_COPROD and 
!                              READ_CH4LOSS (now obsolete)
!  16 Sep 2014 - C. Keller   - Removed all references to Clarisa's OH clima-
!                              tology (not used anywhere)
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  04 Mar 2015 - R. Yantosca - Now use REAL(f4) for HCO_GetPtr pointer args
!------------------------------------------------------------------------------
!BOC
!     
! !DEFINED PARAMETERS:
!
      !========================================================================
      ! Module Variables:
      ! N_CH4      : Number of budget items in TCH4
      ! BAIRDENS   : Array for air density                      [molec/cm3]
      ! BOH        : Array for OH values                        [kg/m3]
      ! COPROD     : Array for zonal mean P(CO)                 [v/v/s]
      ! CH4LOSS    : Array for monthly average CH4 loss freq    [1/s]
      ! PAVG       : Array for 24-h avg surface pressure        [mb]
      ! TAVG       : Array for 24-h avg temperature             [K]
      ! TCH4       : Array for CH4 budget (N_CH4 items)         [molec/box]
      ! FMOL_CH4   : Molecular weight of CH4                    [kg/mole]
      ! XNUMOL_CH4 : Molecules CH4 / kg CH4                     [molec/kg]
      ! CH4_EMIS   : Array for CH4 Emissions                    [kg/m2/s]
      !========================================================================

      ! Number of CH4 budget types 
      INTEGER,  PARAMETER    :: N_CH4     = 12

      REAL(fp), PARAMETER    :: XNUMOL_OH = 6.0221e+23_fp / 17e-3_fp  ! molec OH / kg OH
      REAL(fp), PARAMETER    :: CM3PERM3  = 1.e+6_fp
!
! !LOCAL VARIABLES:
!
      ! Various arrays      
      REAL(fp), ALLOCATABLE :: BAIRDENS(:,:,:)
      REAL(fp), ALLOCATABLE :: PAVG(:,:,:)
      REAL(fp), ALLOCATABLE :: TAVG(:,:,:)

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared as REAL(f4), aka REAL*4.
      REAL(f4), POINTER     :: BOH    (:,:,:) => NULL()
      REAL(f4), POINTER     :: CH4LOSS(:,:,:) => NULL()

      ! Online CH4 only
      REAL(fp), ALLOCATABLE :: SUMANTHRO(:,:)
      REAL(fp), ALLOCATABLE :: SUMBIOG(:,:)

      REAL(fp)              :: TROPOCH4

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
!  17 Sep 2014 - C. Keller   - Removed legacy code (unused variables).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NTDT, NMIN
      INTEGER               :: I, J, L, NTIMES, MNDT
      INTEGER               :: NNEW, NNCOUNT 
      
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
         Pavg(:,:,:) = 0e+0_fp
         Tavg(:,:,:) = 0e+0_fp
         NNEW        = 1
         NNCOUNT     = 0
      ENDIF

      !=================================================================
      ! Archive quantities
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
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
! !IROUTINE: emissch4
!
! !DESCRIPTION: Subroutine EMISSCH4 places emissions of CH4 [kg] into the STT
!  array. (jsw, bnd, bey, bmy, 1/16/01, 10/3/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSCH4( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_XOFFSET, GET_YOFFSET 
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE DIAG_MOD,           ONLY : AD58

      ! HEMCO update
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState, GetHcoDiagn
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
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
!  WARNING: Soil absorption has to be the 11th field in CH4_EMIS
! 
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f" 
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) EMISSCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) GLOBASEAEMIS, GLOBSEAEMIS are diagnostics by jsw.
!  (4 ) Do not multiply CO emissions by 1.28 anymore (jsw, bmy, 2/12/01)
!  (5 ) Renamed input files to CH4_monthly.geos.{RES} and 
!        CH4_aseasonal.geos.{RES}. (bmy, 2/12/01)
!  (6 ) Add reference to "CMN_SETUP" for the DATA_DIR variable (bmy, 2/13/01)
!  (7 ) Removed references to "biofuel_mod.f" and "biomass_mod.f"; these
!        weren't necessary (bmy, 3/20/01)
!  (8 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUNIT as the file unit #. (bmy, 6/27/02)
!  (9 ) Now reference BXHEIGHT and SUNCOS from "dao_mod.f".  Remove reference 
!        to header file "comtrid.h" -- it's not used.  Make FIRSTEMISS a local
!        SAVEd variable.  Also use MONTH from "CMN" instead of the variable
!        LMN. (bmy, 11/15/02)
!  (10) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f".
!        Now use function GET_MONTH and GET_TS_EMIS from "time_mod.f". 
!        Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        I0 and J0 are now local variables. (bmy, 3/27/03)
!  (11) Now reference STT from "tracer_mod.f".  Now reference DATA_DIR from
!        "directory_mod.f". (bmy, 7/20/04)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Add non-local PBL capability (ccc, 8/31/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  11 Apr 2014 - R. Yantosca - Remove call to INIT_GLOBAL_CH4
!  09 Sep 2014 - C. Keller   - Implemented HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                  :: I, I0, IREF, J, J0, JREF, N
      INTEGER                  :: FLAG, ERR
      REAL(fp)                 :: DTSRCE, AREA_M2

      ! Arrays
      REAL(fp), DIMENSION(7)   :: DBGANTHRO
      REAL(fp), DIMENSION(3)   :: DBGBIO

      ! Strings
      CHARACTER(LEN= 63)       :: DgnName
      CHARACTER(LEN=255)       :: MSG
      CHARACTER(LEN=255)       :: LOC = 'EMISSCH4 (global_ch4_mod.F)'

      ! For fields from Input_Opt
      LOGICAL                  :: ITS_A_CH4_SIM
      LOGICAL                  :: ITS_A_FULLCHEM_SIM
      INTEGER                  :: N_TRACERS
      INTEGER                  :: ID_TRACER(Input_Opt%N_TRACERS)

      ! Pointers
      TYPE(HCO_State), POINTER :: HcoState   => NULL()
      REAL(f4),        POINTER :: Ptr2D(:,:) => NULL()

      !=================================================================
      ! EMISSCH4 begins here!
      !=================================================================

      ! Assume success
      RC                  = GIGC_SUCCESS

      ! Get nested-grid offsets
      I0                  = GET_XOFFSET()
      J0                  = GET_YOFFSET()

      ! Copy values from Input_Opt
      N_TRACERS           = Input_Opt%N_TRACERS
      ID_TRACER           = Input_Opt%ID_TRACER(1:N_TRACERS)
      ITS_A_CH4_SIM       = Input_Opt%ITS_A_CH4_SIM
      ITS_A_FULLCHEM_SIM  = Input_Opt%ITS_A_FULLCHEM_SIM

      ! Import emissions from HEMCO (through HEMCO state)
      CALL GetHcoState( HcoState )
      IF ( .NOT. ASSOCIATED(HcoState) ) THEN
         CALL ERROR_STOP ( 'HcoState not defined!', LOC )
      ENDIF

      ! Emission timestep
      DTSRCE = HcoState%TS_EMIS

      IF ( ITS_A_CH4_SIM ) THEN
         print*,'BEGIN SUBROUTINE: EMISSCH4' 
      ENDIF

      ! =================================================================
      ! --> All emission calculations are now done through HEMCO
      ! HEMCO stores emissions of all tracers internally in the HEMCO
      ! state object. Here, we pass these emissions into module array 
      ! CH4_EMIS. These values are then either added to the tracers 
      ! array (full mixing scheme) or used lateron in vdiff_mod.F90 if 
      ! the non-local PBL mixing scheme is used.
      ! The CH4_EMIS array is mostly used for backwards compatibility
      ! (especially the diagnostics). It is also used to ensure that 
      ! in a multi-tracer simulation, tracer 1 (total CH4) is properly
      ! defined. 
      !                                              (ckeller, 9/12/2013)
      ! =================================================================
      CH4_EMIS(:,:,:) = 0e+0_fp

      ! =================================================================
      ! Pass emissions from HEMCO state to CH4_EMIS via HEMCO
      ! diagnostics. Keep units in kg/m2/s.
      ! 
      ! NOTE: the ND60 diagnostics (wetland fraction) is currently not
      ! written. To support this diagnostics, create a manual diagnostics
      ! in the wetland extension (HEMCO/Extensions/hcox_ch4wetland_mod.F90),
      ! activate it in hcoi\_gc\_diagn\_mod.F90 and import it here.
      !                                              (ckeller, 9/12/2013)
      ! =================================================================

      ! Gas and oil
      DgnName = 'CH4_GAS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,2) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Coal 
      DgnName = 'CH4_COAL'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,3) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Livestock
      DgnName = 'CH4_LIVESTOCK'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,4) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Waste 
      DgnName = 'CH4_WASTE'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,5) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Biofuel 
      DgnName = 'CH4_BIOFUEL'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,6) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Rice 
      DgnName = 'CH4_RICE'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,7) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Other anthropogenic 
      DgnName = 'CH4_ANTHROTHER'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,8) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Biomass burning 
      DgnName = 'CH4_BIOMASS'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,9) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Wetland 
      DgnName = 'CH4_WETLAND'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,10) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Soil absorption (those are negative!) 
      DgnName = 'CH4_SOILABSORB'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,11) =  Ptr2D(:,:) * -1.0_hp
      ENDIF
      Ptr2D => NULL()

      ! Other natural 
      DgnName = 'CH4_OTHERNATUR'
      CALL GetHcoDiagn( am_I_Root, DgnName, .FALSE., ERR, Ptr2D=Ptr2D )
      IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP ( TRIM(DgnName), LOC )
      IF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
         CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
      ELSE
         CH4_EMIS(:,:,12) =  Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()


      ! =================================================================
      ! Total emission: sum of all emissions - (2*soil absorption)
      ! We have to substract soil absorption twice because it is added 
      ! to other emissions in the SUM function. (ccc, 7/23/09)
      ! =================================================================
      CH4_EMIS(:,:,1) = SUM(CH4_EMIS, 3) - (2 * CH4_EMIS(:,:,11))

      ! testing only
      write(*,*) 'CH4_EMIS (kg/m2/s):'
      write(*,*) '1 : ', SUM(CH4_EMIS(:,:,1))
      write(*,*) '2 : ', SUM(CH4_EMIS(:,:,2))
      write(*,*) '3 : ', SUM(CH4_EMIS(:,:,3))
      write(*,*) '4 : ', SUM(CH4_EMIS(:,:,4))
      write(*,*) '5 : ', SUM(CH4_EMIS(:,:,5))
      write(*,*) '6 : ', SUM(CH4_EMIS(:,:,6))
      write(*,*) '7 : ', SUM(CH4_EMIS(:,:,7))
      write(*,*) '8 : ', SUM(CH4_EMIS(:,:,8))
      write(*,*) '9 : ', SUM(CH4_EMIS(:,:,9))
      write(*,*) '10: ', SUM(CH4_EMIS(:,:,10))
      write(*,*) '11: ', SUM(CH4_EMIS(:,:,11))
      write(*,*) '12: ', SUM(CH4_EMIS(:,:,12))
!      write(*,*) 'rice 55,29: ',ch4_emis(55,29,7)
!      write(*,*) 'wetl 55,29: ',ch4_emis(55,29,10)

      ! =================================================================
      ! Modify the STT with emissions rates.               (kjw, 5/29/09)
      ! There are 12 tracers in the multi-tracer run.
      ! One tracer for total CH4 and one for each emission category.
      !
      !  1. Total CH4
      !  2. Gas and Oil 
      !  3. Coal
      !  4. Livestock
      !  5. Waste
      !  6. Biofuel
      !  7. Rice
      !  8. Other Anthropogenic
      !  9. Biomass Burning
      !  10. Wetlands
      !  11. Soil Absorption
      !  12. Other Natural
      ! =================================================================
 
      ! Only added to SST for offline CH4 simulations.
      IF ( ITS_A_CH4_SIM ) THEN
       WRITE( 6, '(a)' ) '% EMISSCH4 --- Adding Emissions to STT array.'
      ENDIF

      ! J0 and I0 are global variables, both set = 0.
      DO J = 1, JJPAR
         JREF = J + J0     

      DO I = 1, IIPAR
         IREF = I + I0

         ! Grid box surface area [m2]
         AREA_M2 = HcoState%Grid%AREA_M2%Val(I,J)

         ! The following is now done in mixing_mod.F90 (ckeller, 3/4/15)
!         ! Offline CH4 only
!         ! Convert kg/m2/s to kg.
!         IF ( ITS_A_CH4_SIM ) THEN
!
!            ! Pass directly to tracers array (surface level) if
!            ! non-local scheme is not used.
!            ! For the non-local PBL scheme, CH4_EMIS will
!            ! be used directly in vdiff_mod.F90
!            IF ( .NOT. Input_Opt%LNLPBL ) THEN
!               DO N = 1, N_TRACERS
!                  State_Chm%Tracers(IREF,JREF,1,N)
!     &            = State_Chm%Tracers(IREF,JREF,1,N)
!     &            + CH4_EMIS(IREF,JREF,ID_TRACER(N))
!     &            * DTSRCE * AREA_M2
!               ENDDO
!            ENDIF
!         ENDIF

         ! Update diagnostics
         IF ( ND58 > 0 ) THEN

            ! All emission sources excluding soil absorption
            AD58(IREF,JREF,1) = AD58(IREF,JREF,1) + 
     &         ( CH4_EMIS(IREF,JREF,1) + CH4_EMIS(IREF,JREF,11) )
     &           * DTSRCE * AREA_M2 

            DO N = 2, PD58
               AD58(IREF,JREF,N) = AD58(IREF,JREF,N) + 
     &            CH4_EMIS(IREF,JREF,N) * DTSRCE * AREA_M2
            ENDDO
         ENDIF

      ENDDO
      ENDDO

      IF ( ITS_A_CH4_SIM ) THEN

#if defined( USE_CH4_BUDGET_DIAG )
         !===============================================================
         ! Sum up CH4 budgets
         !
         ! TCH4      - # molecules emitted from different sources
         !===============================================================

         DO J = 1,JJPAR
            JREF = J + J0

         DO I = 1,IIPAR
            IREF = I + I0

            ! Grid box surface area [m2]
            AREA_M2 = HcoState%Grid%AREA_M2%Val(I,J)

            ! Gas, oil, mine
            TCH4(I,J,1,5) = TCH4(I,J,1,5) +
     &           ( ( CH4_EMIS(I,J,2) + CH4_EMIS(I,J,3) ) * 
     &              * XNUMOL_CH4 * AREA_M2 * DTSRCE )
               
            ! agriculture (rice, animals, waste)
            TCH4(I,J,1,6) = TCH4(I,J,1,6) +
     &           ( ( CH4_EMIS(I,J,4) + CH4_EMIS(I,J,7) + 
     &               CH4_EMIS(I,J,5) ) * XNUMOL_CH4 * AREA_M2 * DTSRCE )
            
            ! Biomass burning (and biofuel?)
            TCH4(I,J,1,7) = TCH4(I,J,1,7)+
     &           ( ( CH4_EMIS(I,J,9) + CH4_EMIS(I,J,6) ) *
     &               * XNUMOL_CH4 * AREA_M2 * DTSRCE )

            ! Termites
            TCH4(I,J,1,8) = TCH4(I,J,1,8)+
     &           ( CH4_EMIS(I,J,12) * XNUMOL_CH4 * AREA_M2 * DTSRCE )

            ! Wetlands
            TCH4(I,J,1,9) = TCH4(I,J,1,9)+
     &           ( CH4_EMIS(I,J,10) * XNUMOL_CH4 * AREA_M2 * DTSRCE )

            ! Soil Absorption
            TCH4(I,J,1,10) = TCH4(I,J,1,10)+
     &           ( CH4_EMIS(I,J,11) * XNUMOL_CH4 * AREA_M2 * DTSRCE )


            TCH4(I,J,1,4) = TCH4(I,J,1,4) + 
     &          ( CH4_EMIS(I,J,1)  + ( 2 * CH4_EMIS(IREF,JREF,11) ))
     &          * XNUMOL_CH4 * AREA_M2 * DTSRCE

         ENDDO
         ENDDO
#endif

      ELSEIF ( ITS_A_FULLCHEM_SIM ) THEN

         !===============================================================
         ! Sum up CH4 emissions
         !
         ! SUMANTHRO - Total anthropogenic emissions [molec/cm2/s]
         ! SUMBIOG   - Total biogenic      emissions [molec/cm2/s]
         !===============================================================

         ! Debug counters
         DBGANTHRO = 0e+0_fp
         DBGBIO    = 0e+0_fp

         DO J = 1,JJPAR
            JREF = J + J0

         DO I = 1,IIPAR
            IREF = I + I0

            ! Grid box surface area [m2]
            AREA_M2 = HcoState%Grid%AREA_M2%Val(I,J) 

            ! Total anthropogenic emissions
            ! Gas & oil, mining, livestock, waste, biofuel, rice, other
            SUMANTHRO(I,J) = CH4_EMIS(I,J,2) + CH4_EMIS(I,J,3) +
     &                     + CH4_EMIS(I,J,4) + CH4_EMIS(I,J,5) +
     &                     + CH4_EMIS(I,J,6) + CH4_EMIS(I,J,7) +
     &                     + CH4_EMIS(I,J,8)

            ! kg/s
            DO N=1,7
                DBGANTHRO(N) = DBGANTHRO(N) + 
     &                         (CH4_EMIS(I,J,N+1)*AREA_M2)
            ENDDO

            ! Total biogenic emissions
            ! Wetlands, other natural, (-)soil absorption
            SUMBIOG(I,J) = CH4_EMIS(I,J,10) + CH4_EMIS(I,J,12)
     &                   - CH4_EMIS(I,J,11)

            ! kg/s
            DO N=1,3
                DBGBIO(N) = DBGBIO(N) + 
     &                      (CH4_EMIS(I,J,N+9)*AREA_M2)
            ENDDO

         ENDDO
         ENDDO

         ! Convert from kg/s to kg/yr
         DBGBIO(2)    = DBGBIO(2)*-1e+0_fp ! Soil absorption
         DBGANTHRO(:) = DBGANTHRO(:)*(3600.e-9_fp*24e+0_fp*365.25e+0_fp)
         DBGBIO(:)    = DBGBIO(:)*(3600.e-9_fp*24e+0_fp*365.25e+0_fp)

         ! SDE - Output projected annual emissions based on last output
         WRITE(6,'(a,1x,F10.4,1x,a)') 'CH4 anthrop. emissions rate:',
     &            SUM(DBGANTHRO), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Gas & oil             :',
     &            DBGANTHRO(1), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Mining                :',
     &            DBGANTHRO(2), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Livestock             :',
     &            DBGANTHRO(3), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Waste                 :',
     &            DBGANTHRO(4), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Biofuel               :',
     &            DBGANTHRO(5), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Rice                  :',
     &            DBGANTHRO(6), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Other                 :',
     &            DBGANTHRO(7), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') 'CH4 biogen. emissions rate: ',
     &            SUM(DBGBIO), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Wetlands              :',
     &            DBGBIO(1), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Other natural         :',
     &            DBGBIO(3), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') ' --> Soil absorption       :',
     &            DBGBIO(2), 'Tg/yr'

         WRITE(6,'(a,1x,F10.4,1x,a)') 'CH4 total emissions rate:   ',
     &            (SUM(DBGBIO)+SUM(DBGANTHRO)), 'Tg/yr'
         
      ENDIF

      IF ( ITS_A_CH4_SIM ) THEN
         print*,'END SUBROUTINE: EMISSCH4'
      ENDIF

      END SUBROUTINE EMISSCH4
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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD43
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, ERROR_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, IT_IS_FINITE
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_DAY,   GET_MONTH

      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
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
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  23 Jul 2014 - R. Yantosca - Reference ITS_IN_THE_CHEMGRID (chemgrid_mod.F)
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  17 Sep 2014 - C. Keller   - Now use HEMCO to get CH4 loss and OH conc. field.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Weight of air (taken from "comode.h") 
      REAL(fp), PARAMETER      :: WTAIR = 28.966e+0_fp
!
! !LOCAL VARIABLES:
!
      LOGICAL                :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, K, M, N
      INTEGER                :: NPART

      REAL(fp)                 :: PREVCH4(IIPAR, JJPAR, LLPAR)
      REAL(fp)                 :: BOXVL  

      ! Number of days per month
      INTEGER                :: NODAYS(12) = (/ 31, 28, 31, 30, 
     &                                          31, 30, 31, 31, 
     &                                          30, 31, 30, 31 /)

      ! For fields from Input_Opt
      LOGICAL                :: LSPLIT
      LOGICAL                :: LCH4BUD

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER        :: STT(:,:,:,:)

      ! For error messages
      CHARACTER(LEN=255)    :: LOC = 'CHEMCH4 (global_ch4_mod.F)'

      !=================================================================
      ! CHEMCH4 begins here!
      !=================================================================
      WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4! ---'

      ! testing only
      write(*,*)'Total tracer before chemistry: ',SUM(State_Chm%Tracers)

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy values from Input_Opt
!      N_TRACERS = Input_Opt%N_TRACERS
      LSPLIT    = Input_Opt%LSPLIT
      LCH4BUD   = Input_Opt%LCH4BUD

      !=================================================================
      ! (0) Calculate each box's air density [molec/cm3]
      !        do this for saving mean OH concentrations (kjw, 6/12/09)
      !=================================================================

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Grid box volume [cm3]
         BOXVL           = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         ! Air density [molec/cm3]
         BAIRDENS(I,J,L) = State_Met%AD(I,J,L) * 1000e+0_fp   /
     &                     BOXVL               * 6.023e+23_fp / WTAIR
      ENDDO
      ENDDO
      ENDDO

      !================================================================
      ! (1) Get CH4 loss rates from HEMCO. the target is automatically 
      ! updated by HEMCO (ckeller, 9/16/2014)
      !================================================================
      IF ( FIRSTCHEM ) THEN

         ! Import CH4 loss frequencies from HEMCO. The target will be 
         ! updated automatically by HEMCO (ckeller, 9/16/2014)
         CALL HCO_GetPtr( am_I_Root, 'CH4_LOSS', CH4LOSS, RC )
         IF ( RC /= HCO_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to CH4_LOSS', LOC )

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

      !================================================================
      ! (2) Get OH field from HEMCO. The target will be updated 
      ! automatically by HEMCO (ckeller, 9/16/2014)
      !================================================================
      IF ( FIRSTCHEM ) THEN
         CALL HCO_GetPtr( am_I_Root, 'GLOBAL_OH', BOH, RC )
         IF ( RC /= HCO_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OH', LOC )
      ENDIF

      !=================================================================
      ! (3.1) ND43 diagnostics...save [OH] in molecules/cm3
      ! BOH is in kg/m3 (from HEMCO), convert to molecules/cm3 
      ! (ckeller, 9/16/2014)
      !=================================================================

      IF ( ND43 > 0 ) THEN
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) THEN
!               AD43(I,J,L,1) = AD43(I,J,L,1) + BOH(I,J,L,LMN)
               AD43(I,J,L,1) = AD43(I,J,L,1) + 
     &                       ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 )
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
      CALL CH4_STRAT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

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

         print*,'END CH4_BUDGET'
         call flush(6)
     
      ENDIF
#endif

      print*,'END CHEMCH4'
      call flush(6)

      ! testing only
      write(*,*)'Total tracer after chemistry: ',SUM(State_Chm%Tracers)

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      END SUBROUTINE CHEMCH4
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
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE DIAG_MOD,           ONLY : AD19   	      
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM, ITS_A_NEW_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH
      USE GIGC_State_Chm_Mod, ONLY : ChmState
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
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
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
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L
      REAL(fp)              :: DT, GCH4, STT2GCH4, KRATE
      REAL(fp)              :: BOXVL
      REAL(fp)              :: C_OH

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER     :: STT(:,:,:,:)

      !=================================================================
      ! CH4_DECAY begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM() * 60e+0_fp

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

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( L, J, I, KRATE, STT2GCH4, GCH4, BOXVL, C_OH )
!$OMP+REDUCTION( +:TROPOCH4 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) THEN

            ! Use 24-hr avg temperature to calc. rate coeff.
            ! citation needed
            KRATE = 2.45e-12_fp * EXP( -1775e+0_fp / State_Met%T(I,J,L))

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            STT2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp 
     &                 * XNUMOL_CH4 

            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4

            ! OH in [molec/cm3]. BOH is imported from HEMCO in units of
            ! kg/m3, convert here to molec/cm3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3

#if defined( USE_CH4_BUDGET_DIAG )
            ! Grid box volume [cm3]
            BOXVL         = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! Sum loss in TCH4(3) (molecules/box)
            TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
     &                      ( GCH4  * BOXVL          *
     &                        KRATE * C_OH * DT )
#endif

            TROPOCH4=TROPOCH4 + GCH4 * KRATE * C_OH * DT / STT2GCH4

            ! Modify AD19 Diagnostic
            ! How much CH4 (kg) is lost by reaction with OH
            IF ( ND19 > 0 ) THEN  ! --> [kg/box]
               IF ( L <= LD19 ) THEN
!	    	AD19(I,J,12) = AD19(I,J,12) + 
                  AD19(I,J,L) = AD19(I,J,L) + 
     &                 ( GCH4 * KRATE * C_OH * DT ) / STT2GCH4
               ENDIF
            ENDIF

            ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt) 
            GCH4 = GCH4 * ( 1e+0_fp - KRATE * C_OH * DT )

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
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_CHEMGRID
      USE CMN_SIZE_MOD
      USE DIAG_OH_MOD,        ONLY : DO_DIAG_OH_CH4
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TIME_MOD,           ONLY : GET_MONTH
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
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
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L
      REAL(fp)               :: MASST,   AREA_CM2
      REAL(fp)               :: KCLO,    LOSS,       OHMASS  
      REAL(fp)               :: KCH4,    CH4LOSE,    CH4MASS
      REAL(fp)               :: CH4EMIS, CH4TROPMASS
      REAL(fp)               :: BOXVL
      REAL(fp)               :: C_OH 

      ! Pointers
      REAL(fp), POINTER      :: STT(:,:,:,:)

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! (1) Pass OH mass, total air mass, and  to "diag_oh_mod.f"
      ! (2) ND59: Diagnostic for CH3CCl3 calculation
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT => State_Chm%Tracers

      ! 1. Calculate OH mass and total air mass
      DO L = 1, LLPAR
      DO J = 1, JJPAR 
      DO I = 1, IIPAR 

         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) THEN

            ! Grid box volume [cm3]
            BOXVL  = State_Met%AIRVOL(I,J,L) * 1e+6_fp

            ! OH concentration in molec/cm3. BOH is imported
            ! from HEMCO in kg/m3 (ckeller, 9/16/2014)
            C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3

            ! Calculate OH mass [molec / box]
            OHMASS = C_OH * BAIRDENS(I,J,L) * BOXVL

            ! Calculate total air mass [molec / box]
            MASST  = BAIRDENS(I,J,L) * BOXVL

            ! Calculate CH3CCl3 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCLO = 1.64e-12_fp * EXP( -1520.e+0_fp / State_Met%T(I,J,L))

            ! Calculate Loss term [molec / box / s]
            LOSS   = KCLO            * C_OH  *
     &               BAIRDENS(I,J,L) * BOXVL


            ! Calculate CH4 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCH4 = 2.45e-12_fp * EXP( -1775e+0_fp / State_Met%T(I,J,L) )

            ! Calculate CH4 mass [molec / box] from [kg / box]
            CH4TROPMASS = STT(I,J,L,1) * XNUMOL_CH4 
            CH4MASS     = STT(I,J,L,1) * XNUMOL_CH4 

            ! Calculate loss term  [molec /box / s]
            CH4LOSE = KCH4            * C_OH  *
     &                BAIRDENS(I,J,L) * BOXVL

            ! Calculate CH4 emissions [molec / box / s]
            !   Only for surface level
            ! Grid box surface area [cm2]
            ! HEMCO update: CH4_EMIS now in kg/m2/s (ckeller, 9/12/2014)
            IF ( L .GT. 1 ) THEN 
               CH4EMIS = 0e+0_fp
            ELSE
               AREA_CM2 = GET_AREA_CM2( I, J, L ) / 10000.0e+0_fp         ! now in m2!!
               CH4EMIS  = SUM(CH4_EMIS(I,J,2:10)) + CH4_EMIS(I,J,12)
               CH4EMIS  = CH4EMIS * AREA_CM2                          ! [kg/m2/s]  --> [kg/box/s]
               CH4EMIS  = CH4EMIS * XNUMOL_CH4                        ! [kg/box/s] --> [molec/box/s]
            ENDIF

         ELSE

            OHMASS      = 0e+0_fp
            MASST       = 0e+0_fp
            LOSS        = 0e+0_fp
            CH4LOSE     = 0e+0_fp
            CH4TROPMASS = 0e+0_fp
            CH4EMIS     = 0e+0_fp
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
      SUBROUTINE CH4_STRAT( am_I_Root, Input_Opt, 
     &                      State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CHEMGRID_MOD,       ONLY : GET_MIN_TPAUSE_LEVEL
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD,         ONLY : CHECK_STT
      USE TRACERID_MOD,       ONLY : IDTCH4
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
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
!  Production (mixing ratio/sec) rate provided by Dylan Jones.  
!  Only production by CH4 + OH is considered.
!                                                                             .
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
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
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root and RC
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Weight of air (taken from "comode.h") 
      REAL(fp), PARAMETER   :: WTAIR = 28.966e+0_fp
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, BASELEV
      REAL(fp)              :: DT, GCH4, STT2GCH4, LRATE
      REAL(fp)              :: BOXVL

      ! Local variables for quantities from Input_Opt
      LOGICAL             :: LUCX
      LOGICAL             :: ITS_A_FULLCHEM_SIM

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER     :: STT(:,:,:,:)

!      LOGICAL, SAVE       :: FIRSTCHEM = .TRUE.

      !=================================================================
      ! CH4_STRAT begins here!
      !=================================================================

      ! Assume success
      RC                  =  GIGC_SUCCESS

      ! Copy fields from INPUT_OPT
      LUCX                = Input_Opt%LUCX
      ITS_A_FULLCHEM_SIM  = Input_Opt%ITS_A_FULLCHEM_SIM

      ! If unified chemistry is active, ignore all of this
      IF ( .not. LUCX ) THEN

         ! Chemistry timestep [s]
         DT  = GET_TS_CHEM() * 60e+0_fp

         IF ( ITS_A_FULLCHEM_SIM ) THEN
            
            !=============================================================
            ! (0) Calculate each box's air density [molec/cm3]
            !=============================================================

            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box volume [cm3]
               BOXVL           = State_Met%AIRVOL(I,J,L) * 1e+6_fp

               ! Air density [molec/cm3]
               BAIRDENS(I,J,L) = State_Met%AD(I,J,L) * 1e+3_fp      / 
     &                           BOXVL               * 6.023e+23_fp /
     &                           WTAIR
            ENDDO
            ENDDO
            ENDDO

            !=============================================================
            ! (1) If the first time step...
            ! COPROD doesn't seem to be used anymore at all. Comment out
            ! for now, may be completely removed later. If necessary, we
            ! can also read the COPROD array through HEMCO (as CH4LOSS). 
            ! (ckeller, 9/15/2014)
            !=============================================================
!            IF ( FIRSTCHEM ) THEN
!               ! Zero CO production array
!               COPROD(:,:,:) = 0e+0_fp
!               !print*,'READ_COPROD'
!               ! Read zonally-averaged CO production [v/v/s]
!               CALL READ_COPROD( am_I_Root, Input_Opt, RC )
!               !print*,'READ_COPROD DONE'
!            ENDIF
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
            IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met) ) THEN

               ! Conversion factor [kg/box] --> [molec/cm3]
               ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
               STT2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp *
     &                    XNUMOL_CH4

               ! CH4 in [molec/cm3]
               GCH4 = STT(I,J,L,1) * STT2GCH4

               ! Loss rate [molec/cm3/s]
!               LRATE = GCH4 * CH4LOSS( I,J,L,LMN )
               LRATE = GCH4 * CH4LOSS( I,J,L )

               ! Update Methane concentration in this grid box [molec/cm3]
               GCH4 = GCH4 - ( LRATE * DT )

! Update stratospheric chem to use linearized CH4 loss frequencies
!  (kjw, 11/19/11)
!               BOXVL = State_Met%AIRVOL(I,J,L) * 1d6
!
!               ! Sum loss in TCH4(3) [molec CH4/box] in the stratosphere
!               ! [molec/cm3] * [v/v/s] * [s] * [cm3/box] = [molec CH4/box]
!               TCH4(I,J,L,3) = TCH4(I,J,L,3) + 
!        &                      ( BAIRDENS(I,J,L) * COPROD(J,L,LMN) *
!        &                        DT              * BOXVL  )    )
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

!      FIRSTCHEM = .FALSE.

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
      INTEGER             :: I, J, K, L, M, NERROR, UD

      REAL(fp)              :: STTCONV, TGS, SCALEDYN
      REAL(fp)              :: NTP, NTQ, NTP2, NTQ2 
      REAL(fp)              :: SOURCES, SINKS

      CHARACTER(LEN=17)   :: MERGE

      ! For binary punch file, v. 2.0
!      REAL*4              :: ARRAY(IIPAR, JJPAR, LLPAR)
!      REAL*4              :: LONRES, LATRES
!
!      INTEGER             :: IFIRST, JFIRST, LFIRST
!      INTEGER, PARAMETER  :: HALFPOLAR = 1
!      INTEGER, PARAMETER  :: CENTER180 = 1
!
!      CHARACTER (LEN=20)  :: MODELNAME 
!      CHARACTER (LEN=40)  :: UNIT
!      CHARACTER (LEN=40)  :: RESERVED = ''
!      CHARACTER (LEN=40)  :: CATEGORY 
!      CHARACTER (LEN=80)  :: LABEL

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CH4_BUDGET begins here!
      !
      ! Initialize quantities 
      !=================================================================
!      IFIRST    = GET_XOFFSET() + 1
!      JFIRST    = GET_YOFFSET() + 1
!      LFIRST    = 1
!      LONRES    = DISIZE
!      LATRES    = DJSIZE
!
!      ! Make up a category name for GAMAP (use 8 characters)
!      CATEGORY  = 'CH4BUDT'

!      ! Get the proper model name for the binary punch file
!      MODELNAME = GET_MODELNAME()

!      ! Descriptor string
!      LABEL    = 'GEOS-CHEM -- CH4 Budget output (jsw, bmy, 1/16/01)'

!      ! Unit of quantity being saved
!      UNIT      = 'Tg'  !(NOTE: check w/ bnd to get the right units!!!)

      ! Scale factor for dynamic time steps
      SCALEDYN  = FLOAT( GET_CT_DYN() ) + 1e-20_fp

      !=================================================================
      ! Store the final burden of CH4 in TCH4(2) 
      ! Convert kg CH4/box to molecules/box.
      !=================================================================
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      TCH4(:,:,:,2) = 0e+0_fp
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
      
      TGS     = 1.e-9_fp
      STTCONV = XNUMOL_CH4/TGS
      SOURCES = 0.e+0_fp
      SINKS   = 0.e+0_fp
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

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1,State_Met)
      WRITE(189,20)NTP,NTP/STTCONV

      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1,State_Met)
      WRITE(189,21)NTP2,NTP2/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
 1991 FORMAT('Stratospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0,State_Met)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0,State_Met)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks   jsw has checked correctness of code for sinks.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,1,State_Met)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,0,State_Met)
      SINKS=NTP+NTQ

      WRITE(189,22) NTP,NTP/SINKS*100.e+0_fp,NTP/STTCONV
      WRITE(189,220) NTQ,NTQ/SINKS*100.e+0_fp,NTQ/STTCONV
      WRITE(189,29) 
      WRITE(189,34) SINKS,SINKS/STTCONV  !Just OH sink 
      WRITE(189,18)
      WRITE(189,30)

      ! Sources
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,4,4,1,State_Met)
      SOURCES=NTP

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,5,5,1,State_Met)
      WRITE(189,24) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,9,9,1,State_Met)
      WRITE(189,27) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,6,6,1,State_Met)
      WRITE(189,39) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,7,7,1,State_Met)
      WRITE(189,25) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,8,8,1,State_Met)
      WRITE(189,26) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      !jsw Following lines added by jsw.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,10,10,1,State_Met)
      WRITE(189,35) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      SINKS=SINKS-NTP  !Minus sign because soil absorption is negative.

      WRITE(189,29) 
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1,State_Met)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1,State_Met)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0,State_Met)
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

      SOURCES = 0.e+0_fp
      SINKS   = 0.e+0_fp

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,36)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1,State_Met)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1,State_Met)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0,State_Met)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0,State_Met)
      WRITE(189,21) NTP,NTP/STTCONV
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF( NTP > 0e+0_fp) SINKS = SINKS + NTP

      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/
     &     2+1,1,LLPAR,11,11,0,State_Met)
      IF( NTP2 > 0e+0_fp ) SINKS = SINKS + NTP2      

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0,State_Met)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.e+0_fp,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.e+0_fp,NTQ2/STTCONV
 
      IF(NTP.GE.0.e+0_fp) THEN
         WRITE(189,270) NTP,NTP/SINKS*100.e+0_fp,NTP/STTCONV
      ENDIF

      IF(NTP.GE.0.e+0_fp) THEN
         WRITE(189,2700) NTP,NTP/SINKS*100.e+0_fp,NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,34) SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)
      
      ! Sources
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,5,9,1,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0,State_Met)
      SOURCES=NTQ+NTQ2
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         SOURCES=SOURCES-NTP
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,1,State_Met)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0,State_Met)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,5,5,1,State_Met)
      WRITE(189,24) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,6,6,1,State_Met)
      WRITE(189,39) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,7,7,1,State_Met)
      WRITE(189,25) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,8,8,1,State_Met)
      WRITE(189,26) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,9,9,1,State_Met)
      WRITE(189,27) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         WRITE(189,270) -NTP,-NTP/SOURCES*100.e+0_fp,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,
     &     JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      IF( NTP2 > 0e+0_fp ) SINKS = SINKS + NTP2      

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0,State_Met)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.e+0_fp,NTQ/STTCONV

      IF(NTP.LT.0.e+0_fp) THEN
         WRITE(189,2700) -NTP,-NTP/SOURCES*100.e+0_fp,-NTP/STTCONV
      ENDIF
      
      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)
      
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1,State_Met)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1,State_Met)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0,State_Met)
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

      SOURCES = 0.e+0_fp
      SINKS   = 0.e+0_fp

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,37)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1,State_Met)
      WRITE(189,20) NTP,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1,State_Met)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,1991)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0,State_Met)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0,State_Met)
      WRITE(189,21) NTP,NTP/STTCONV
      
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,1,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,0,State_Met)
      SINKS=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         SINKS=SINKS-NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         SINKS=SINKS-NTP
      ENDIF

      WRITE(189,22) NTQ,NTQ/SINKS*100.e+0_fp,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.e+0_fp,NTQ2/STTCONV
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         WRITE(189,270) -NTP,-NTP/SINKS*100.e+0_fp,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      IF(NTP.LT.0.e+0_fp) THEN
         WRITE(189,2700) -NTP,-NTP/SINKS*100.e+0_fp,-NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,34)SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)

      ! Sources
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,5,9,1,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0,State_Met)
      SOURCES=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.GE.0.e+0_fp) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,,State_Met0)
      IF(NTP.GE.0.e+0_fp) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,1,State_Met)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0,State_Met)
      
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,5,5,1,State_Met)
      WRITE(189,24) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,6,6,1,State_Met)
      WRITE(189,39) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,7,7,1,State_Met)
      WRITE(189,25) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,8,8,1,State_Met)
      WRITE(189,26) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,9,9,1,State_Met)
      WRITE(189,27) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1,State_Met)
      IF(NTP.GE.0.e+0_fp) THEN
         WRITE(189,270) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0,State_Met)
      IF(NTP.GE.0.e+0_fp) THEN
         WRITE(189,2700) NTP,NTP/SOURCES*100.e+0_fp,NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1,State_Met)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1,State_Met)
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0,State_Met)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0,State_Met)
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
!        ! Cast REAL(fp) into REAL*4, convert from molec to Tg
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
      TCH4(:,:,:,2:N_CH4) = 0e+0_fp
	
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
      REAL(fp) FUNCTION SUM_CH4( I1, I2, J1, J2, L1, L2, K1, K2, UPDOWN,
     &                         State_Met )
!
! !USES:
!
      USE CHEMGRID_MOD,      ONLY : ITS_IN_THE_CHEMGRID
      USE CHEMGRID_MOD,      ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_SIZE_MOD
      USE GIGC_State_Met_Mod
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I1, I2 ! Min/max longitude indices of TCH4 to sum
      INTEGER, INTENT(IN) :: J1, J2 ! Min/max latitude  indices of TCH4 to sum 
      INTEGER, INTENT(IN) :: L1, L2 ! Min/max altitude  indices of TCH4 to sum
      INTEGER, INTENT(IN) :: K1, K2 ! Min/max tracer    indices of TCH4 to sum
      INTEGER, INTENT(IN) :: UPDOWN ! Sum in trop (=1) or in strat (=0)

      ! Also need to pass the State_Met for the CHEMGRID functons
      TYPE(MetState), INTENT(IN) :: State_Met
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
!  23 Jul 2014 - R. Yantosca - Now pass the State_Met for chemgrid functions

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, K, L

      !=================================================================
      ! SUM_CH4 begins here!
      !=================================================================

      !### Debug
      !print*,'LPAUSE MIN/MAX=',LPAUSE_MIN,LPAUSE_MAX  
      !print*,'L1,L2=',L1,L2
      
      ! Initialize SUM_CH4
      SUM_CH4 = 0e+0_fp

      ! Test on UPDOWN
      IF ( UPDOWN == 1 ) THEN

         !=============================================================
         ! UPDOWN = 1: Sum up from the surface to the tropopause
         !=============================================================
         DO K = K1, K2
         DO L = 1,  LLPAR
         DO J = J1, J2
         DO I = I1, I2
            IF ( ITS_IN_THE_CHEMGRID( I, J, L, State_Met ) ) THEN
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
         DO K = K1, K2
         DO L = 1,  LLPAR
         DO J = J1, J2
         DO I = I1, I2
            IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) THEN
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
      REAL(fp)                        :: PREVCH4(IIPAR,JJPAR,LLPAR)! CH4 bef chem
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
      REAL(fp), POINTER :: STT(:,:,:,:)

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
           STT(I,J,L,N) = SAFE_DIV(STT(I,J,L,N),PREVCH4(I,J,L),0.e+0_fp)
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
      REAL(fp) :: CH4_ANTHRO   ! Anthro CH4 emissions [molec/cm2/s]
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
         CH4_ANTHRO = -1e+0_fp
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
      REAL(fp)  :: CH4_BIOG   ! Biogenic CH4 emissions [molec/cm2/s]
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
         CH4_BIOG = -1e+0_fp
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

!      ALLOCATE( AVGOH( NSEAS, NCMSLATS, NCMSALTS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVGOH' )
!      AVGOH = 0e+0_fp

      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BAIRDENS' )
      BAIRDENS = 0e+0_fp

!      ALLOCATE( BOH( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BOH' )
!      BOH = 0e+0_fp

!      ALLOCATE( COPROD( JJPAR, LLPAR, 12 ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
!      COPROD = 0e+0_fp

!      ALLOCATE( CH4LOSS( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
!      ALLOCATE( CH4LOSS( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4LOSS' )
!      CH4LOSS = 0e+0_fp

      ALLOCATE( PAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PAVG' )
      PAVG = 0e+0_fp

      ALLOCATE( TAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAVG' )
      TAVG = 0e+0_fp

#if defined( USE_CH4_BUDGET_DIAG )
      ALLOCATE( TCH4( IIPAR, JJPAR, LLPAR, N_CH4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCH4' )
      TCH4 = 0e+0_fp      
#endif

      ALLOCATE( CH4_EMIS( IIPAR, JJPAR, PD58), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_EMIS' )
      CH4_EMIS = 0e+0_fp

      ALLOCATE( SUMANTHRO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMANTHRO' )
      SUMANTHRO = 0e+0_fp

      ALLOCATE( SUMBIOG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUMBIOG' )
      SUMBIOG = 0e+0_fp

      ! Initialize tropoch4 (counts total decay of CH4 due to OH)
      TROPOCH4 = 0e+0_fp

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
!      IF ( ALLOCATED( BOH       ) ) DEALLOCATE( BOH       )
!      IF ( ALLOCATED( COPROD    ) ) DEALLOCATE( COPROD    )
!      IF ( ALLOCATED( CH4LOSS   ) ) DEALLOCATE( CH4LOSS   )
      IF ( ASSOCIATED( BOH      ) ) BOH     => NULL()
      IF ( ASSOCIATED( CH4LOSS  ) ) CH4LOSS => NULL()
#if defined( USE_CH4_BUDGET_DIAG )
      IF ( ALLOCATED( TCH4      ) ) DEALLOCATE( TCH4      )
#endif
      IF ( ALLOCATED( CH4_EMIS  ) ) DEALLOCATE( CH4_EMIS  )
      IF ( ALLOCATED( SUMANTHRO ) ) DEALLOCATE( SUMANTHRO )
      IF ( ALLOCATED( SUMBIOG   ) ) DEALLOCATE( SUMBIOG   )

      END SUBROUTINE CLEANUP_GLOBAL_CH4
!EOC
      END MODULE GLOBAL_CH4_MOD
