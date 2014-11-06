!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tagged_ox_mod
!
! !DESCRIPTION: Module TAGGED\_OX\_MOD contains variables and routines to 
!  perform a tagged Ox simulation.  P(Ox) and L(Ox) rates need to be archived 
!  from a full chemistry simulation before you can run w/ Tagged Ox. 
!\\
!\\
! !INTERFACE: 
!
      MODULE TAGGED_OX_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD   ! HEMCO precision definition and error handling

      IMPLICIT NONE
      PRIVATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% If you want to the EXTENDED SIMULATION with all 13 tagged O3 tracers, 
!%%% then uncomment this #ifdef statement. (bmy, 4/11/14)
!%%% tracers. (bmy, 4/11/14)
!#define USE_ALL_TAGO3_TRACERS 1
!%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% If you want to read archived O3 prod/loss data from a different year
!%%% than the met field dates, then use this #ifdef.  You can also leave this
!%%% flag commented out and set its value by compiling with the TAGO3YR=2004
!%%% (etc. for other years) Makefile option.  (bmy, 4/24/14)
!#define USE_THIS_O3_YEAR 2004
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: ADD_STRAT_POX 
      PUBLIC  :: CHEM_TAGGED_OX 
      PUBLIC  :: CLEANUP_TAGGED_OX 
      PUBLIC  :: INIT_TAGGED_OX
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: GET_REGIONAL_POX
!
! !REMARKS:
!  THE SIMPLE TAGGED O3 SIMULATION (default setting) HAS THESE TRACERS:
!  ----------------------------------------------------------------------------
!  (1 ) O3      : Total O3
!  (2 ) O3Strt  : Stratospheric O3
!                                                                             .
!  THE EXTENDED TAGGED O3 SIMULATION HAS THESE TRACERS:
!  ----------------------------------------------------------------------------
!  (1 ) O3      : Total O3
!  (2 ) O3Strt  : O3 from the Stratosphere      (tropopause - atm top   )
!  (3 ) O3Ut    : O3 produced in Upper Trop     (350 hPa    - tropopause) 
!  (4 ) O3Mt    : O3 produced in Middle Trop    (PBL top    - 350 hPa   )
!  (5 ) O3Row   : O3 produced in Rest of World  (surface    - PBL top   )
!  (6 ) O3PcBl  : O3 produced in Pacific BL     (surface    - PBL top   )
!  (7 ) O3NaBl  : O3 produced in N. American BL (surface    - PBL top   )
!  (8 ) O3AtBl  : O3 produced in Atlantic BL    (surface    - PBL top   )
!  (9 ) O3EuBl  : O3 produced in European BL    (surface    - PBL top   )
!  (10) O3AfBl  : O3 produced in N. African BL  (surface    - PBL top   )
!  (11) O3AsBl  : O3 produced in Asian          (surface    - PBL top   )
!  (12) O3Init  : O3 initial conditions         (all levels             )
!  (13) O3USA   : O3 produced over the USA      (all levels             )
!                                                                             .
!  NOTES:
!  ----------------------------------------------------------------------------
!  (1) The stratospheric O3 tracer must be tracer #2.  This is due to how
!       the Linoz stratospheric O3 chemistry scheme is written.  We have
!       accordingly reorganized the tracer numbers below.
!  (2) The name "tagged_ox_mod.F" is historical.  The Ox tracer in GEOS-Chem 
!       has now been replaced with O3.  O3 usually makes up about 95% of Ox.
!       The nomenclature "tagged Ox" is interchangeable with "tagged O3".
!
! !REVISION HISTORY:
!  20 Aug 2003 - A. Fiore    - Initial version  
!  (1 ) Now accounts for GEOS-4 PBL being in meters (bmy, 1/15/04)
!  (2 ) Bug fix: don't put function call in WRITE statement (bmy, 2/20/04)
!  (3 ) Now bracket AD44 with an !$OMP CRITICAL block (bmy, 3/24/04)
!  (4 ) Now define regions w/ levels in GET_REGIONAL_POX (amf,rch,bmy,5/27/04)
!  (5 ) Bug fix-avoid seg fault if PBLFRAC isn't allocated (bdf, bmy, 10/12/04)
!  (6 ) Now reference "pbl_mix_mod.f" (bmy, 2/17/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (9 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (10) Modified for variable tropopause (phs, bmy, 1/19/07)
!  (11) Now use LLTROP instead of LLTROP_FIX everywhere (bmy, 12/4/07)
!  (12) Now use LD65 instead of LLTROP everywhere (phs, 11/17/08)
!  (13) Updates for LINOZ (dbj, jliu, bmy, 10/26/09)
!  19 Nov 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  23 Jan 2014 - M. Sulprizio- Now use LLTROP instead of LD65 or LLTROP_FIX
!                              everywhere (hyl, bmy, 11/3/11)
!  11 Apr 2014 - R. Yantosca - Add #ifdef USE_ALL_TAGO3_TRACERS to select
!                              between the simple and extended simulations
!  11 Apr 2014 - R. Yantosca - Updated comments
!  11 Apr 2014 - R. Yantosca - Now make the O3Strt tracer #2 for both the
!                              simple and extended tagged O3 simulations
!  11 Apr 2014 - R. Yantosca - Now make INIT_TAGGED_OX a public routine
!  17 Sep 2014 - C. Keller   - Now use HEMCO for prod./loss rates. This makes
!                              subroutine READ_POX_LOX obsolete.
!  06 Nov 2014 - R. Yantosca - Removed code orphaned by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      REAL(hp), POINTER        :: P24H(:,:,:) => NULL() ! O3 production rate
      REAL(hp), POINTER        :: L24H(:,:,:) => NULL() ! O3 loss rate

      ! Emission timestep (will be imported from HEMCO)
      REAL*8                   :: TS_EMIS

      ! To convert m3 to cm3
      REAL*8, PARAMETER        :: CM3PERM3 = 1.0d6
!
! !DEFINED PARAMETERS:
! 
#if defined( USE_ALL_TAGO3_TRACERS )

      !-----------------------------------------------------------------
      ! EXTENDED SIMULATION : Total, strat, and regional O3 tracers
      !-----------------------------------------------------------------
      INTEGER, PARAMETER,   PRIVATE :: N_TAGGED = 13   ! # of tracers
      INTEGER, PARAMETER,   PRIVATE :: N_STRAT  = 2    ! Stratospheric O3

#else

      !-----------------------------------------------------------------
      ! SIMPLE SIMULATION: Total and strat O3 tracers only
      !
      ! %%% THIS IS THE DEFAULT %%%
      !-----------------------------------------------------------------
      INTEGER, PARAMETER,   PRIVATE :: N_TAGGED = 2    ! # of tracers
      INTEGER, PARAMETER,   PRIVATE :: N_STRAT  = 2    ! Stratospheric O3

#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_strat_pox
!
! !DESCRIPTION: Subroutine ADD\_STRAT\_POX adds the stratospheric influx of 
!  Ox to the stratospheric Ox tracer.  This is called from routine 
!  Do\_Synoz, which is applied when the tracer array has units of [v/v].  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ADD_STRAT_POX( I, J, L, POx, State_Chm )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTO3Strt

!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I     ! GEOS-Chem grid box lon index
      INTEGER, INTENT(IN) :: J     ! GEOS-Chem grid box lat index
      INTEGER, INTENT(IN) :: L     ! GEOS-Chem grid box level index
      REAL*8,  INTENT(IN) :: POx   ! P(Ox) in the stratosphere [v/v]
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REMARKS:
!  This routine is only called if you are using the Synoz stratospheric
!  ozone chemistry.  Synoz is now obsolete and has been replaced by Linoz.
!
! !REVISION HISTORY:
!  19 Aug 2003 - R. Yantosca - Initial version
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  10 Jun 2012 - L. Murray - Enable tagged Ox for Synoz.
!  05 Nov 2013 - R. Yantosca - Rename IDTOxStrt to IDTO3Strt
!  05 Nov 2013 - R. Yantosca - Remove STT pointer for simplicity
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! GET_STRAT_POX begins here!
      !=================================================================

      State_Chm%Tracers(I,J,L,IDTO3Strt) = 
     &State_Chm%Tracers(I,J,L,IDTO3Strt) + POx

      END SUBROUTINE ADD_STRAT_POX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_regional_pox
!
! !DESCRIPTION: Subroutine GET\_REGIONAL\_POX returns the P(Ox) for each of 
!  the tagged Ox tracers. Tagged Ox tracers are defined by both geographic 
!  location and altitude. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_REGIONAL_POX( I, J, L, PP, State_Met )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_XMID,  GET_YMID
      USE TIME_MOD,           ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! ND44, ND65, LD65
      USE CMN_GCTM_MOD             ! SCALE_HEIGHT
!
! !INPUT PARAMETERS: 
!
      ! GEOS-Chem grid box indices for lon, lat, alt
      INTEGER,        INTENT(IN)  :: I, J, L

      ! Meteorology State object
      TYPE(MetState), INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS:
!
      ! Array containing P(Ox) for each tagged tracer
      REAL*8,  INTENT(OUT) :: PP(IIPAR,JJPAR,LLTROP,N_TAGGED)
! 
! !REVISION HISTORY:
!  19 Aug 2003 - A. Fiore - Initial version
!  (1 ) Updated from the old routine "chemo3_split.f" (rch, bmy, 8/20/03)
!  (2 ) For GEOS-4, convert PBL from [m] to [hPa] w/ the hydrostatic law.
!        Now references SCALE_HEIGHT from "CMN_GCTM". (bmy, 1/15/04)
!  (3 ) Now uses model levels instead of pressure in order to delineate
!        between PBL, MT, and UT regions (amf, rch, bmy, 5/27/04)
!  (4 ) Now references ITS_IN_THE_TROP from "tropopause_mod.f".  Now remove
!        reference to "CMN", it's obsolete. (bmy, 8/22/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (6 ) Resize the PP array from LLTROP to LLTROP_FIX (phs, 1/19/07)
!  (7 ) Now use LLTROP instead of LLTROP_FIX (bmy, 12/4/07)
!  (8 ) Now use LD65 instead of LLTROP (phs, 11/17/08)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_XMID(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  23 Jan 2014 - M. Sulprizio- Now use LLTROP instead of LD65 (hyl,bmy,11/3/11)
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_IN_TROP, ITS_IN_PBL, ITS_IN_MT
      LOGICAL :: ITS_IN_UT,   ITS_IN_NH,  ITS_IN_ATL
      LOGICAL :: ITS_IN_PAC,  ITS_IN_AS,  ITS_IN_EUR
      LOGICAL :: ITS_IN_NAM,  ITS_IN_NAF, ITS_IN_USA
      INTEGER :: PBLTOP,      MTTOP
      REAL*8  :: PPROD,       X,          Y
      REAL*8  :: BOXVL

      !=================================================================
      ! Initialization
      !=================================================================

      ! Initialize
      PP(I,J,L,:) = 0d0
      
      ! IS TROP is TRUE if we are in the troposphere
      ITS_IN_TROP = ITS_IN_THE_TROP( I, J, L, State_Met )
      
      ! Skip stratospheric boxes
      IF ( .not. ITS_IN_TROP ) RETURN

      ! Longitude and latitude [degrees]
      X          = GET_XMID( I, J, L )   
      Y          = GET_YMID( I, J, L )

      ! PBLTOP is the model level at ~ 750 hPa
      ! MTTOP  is the model level at ~ 350 hPa
#if   defined( GEOS_4 )
      PBLTOP = 5
      MTTOP  = 10
#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      PBLTOP = 5
      MTTOP  = 10
#endif

#if defined( USE_ALL_TAGO3_TRACERS )

      !=================================================================
      ! For the simulation with all tagged O3 tracers: find the 
      ! regions corresponding to this particular longitude & latitude
      !=================================================================

      ! Define flags for various geographic & altitude regions
      ITS_IN_PBL = ( L <= PBLTOP                                       )
      ITS_IN_MT  = ( L >  PBLTOP .and. L <= MTTOP                      )
      ITS_IN_UT  = ( L >  MTTOP  .and. ITS_IN_TROP                     )

      ITS_IN_NH  = ( Y >=   0.0                                        )
      ITS_IN_EUR = ( Y >=  36.0 .and. ( X >  -15.0 .and. X >=   55.0 ) )
      ITS_IN_NAM = ( Y >=  15.0 .and. ( X > -127.5 .and. X <=  -65.0 ) )
      ITS_IN_AS  = ( Y >= -10.0 .and. ( X >   55.0 .and. X <=  145.0 ) )
      ITS_IN_ATL = ( ITS_IN_NH  .and. ( X >  -65.0 .and. X <=  -15.0 ) )
      ITS_IN_PAC = ( ITS_IN_NH  .and. ( X >  145.0  .or. X <= -127.5 ) )

      ITS_IN_NAF = ( ( X >= -15.0 .and. X <=  55.0 ) .and. 
     &               ( Y >=   0.0 .and. Y <   36.0 ) )  

      ITS_IN_USA = ( ( X > -127.5 .and. X <= -65.0 ) .and. 
     &               ( Y >   22.0 .and. Y <=  50.0 ) )

#endif

      !=================================================================
      ! Assign P(Ox) to tagged tracers by geographic/altitude regions
      !=================================================================

      ! Grid box volume [cm3]
      BOXVL = State_Met%AIRVOL(I,J,L) !* 1d6

      ! P(Ox) [kg]
      ! P24H is in kg/m3 per emission time step (ckeller, 9/17/2014).
      PPROD = P24H(I,J,L) * BOXVL * ( GET_TS_CHEM() * 60d0 / TS_EMIS )

      !-----------------------
      ! #1: Total P(Ox)
      !-----------------------
      PP(I,J,L,1) = PPROD

#if defined( USE_ALL_TAGO3_TRACERS )

      !-----------------------
      ! #2: P(Ox) in UT
      !-----------------------
      IF ( ITS_IN_UT ) THEN
         PP(I,J,L,3) = PPROD
         
      !-----------------------
      ! #3: P(Ox) in MT 
      !-----------------------
      ELSE IF ( ITS_IN_MT ) THEN
         PP(I,J,L,4) = PPROD
                                
      !-----------------------
      ! #5: P(Ox) in Pac BL
      !-----------------------
      ELSE IF ( ITS_IN_PAC .and. ITS_IN_PBL ) THEN
         PP(I,J,L,6) = PPROD

      !-----------------------
      ! #6: P(Ox) in NAm BL
      !-----------------------
      ELSE IF ( ITS_IN_NAM .and. ITS_IN_PBL ) THEN     
         PP(I,J,L,7) = PPROD
                  
      !-----------------------
      ! #7: P(Ox) in Atl BL
      !-----------------------
      ELSE IF ( ITS_IN_ATL .and. ITS_IN_PBL ) THEN
         PP(I,J,L,8) = PPROD  
         
      !-----------------------
      ! #8: P(Ox) in Eur BL
      !-----------------------
      ELSE IF ( ITS_IN_EUR .and. ITS_IN_PBL ) THEN
         PP(I,J,L,9) = PPROD
                  
      !-----------------------
      ! #9: P(Ox) in NAfr BL
      !-----------------------
      ELSE IF ( ITS_IN_NAF .and. ITS_IN_PBL ) THEN
         PP(I,J,L,10) = PPROD
 
      !-----------------------
      ! #10: P(Ox) in Asia BL
      !-----------------------          
      ELSE IF ( ITS_IN_AS .and. ITS_IN_PBL ) THEN
         PP(I,J,L,11) = PPROD                   

      !-----------------------
      ! #4: P(Ox) in R.O.W
      !-----------------------
      ELSE 
         PP(I,J,L,5) = PPROD

      ENDIF

      !-------------------------
      ! #13: P(Ox) in USA
      !-------------------------
      IF ( ITS_IN_USA ) THEN
         PP(I,J,L,13) = PPROD               
      ENDIF

#endif

      END SUBROUTINE GET_REGIONAL_POX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_tagged_ox 
!
! !DESCRIPTION: Subroutine CHEM\_TAGGED\_OX performs chemistry for several 
!  Ox tracers which are tagged by geographic and altitude regions. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_TAGGED_OX( am_I_Root, Input_Opt,
     &                           State_Met, State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY : AD44
      USE DIAG_PL_MOD,        ONLY : AD65
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE DRYDEP_MOD,         ONLY : DEPSAV      
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,        ONLY : GET_PBL_MAX_L
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTO3
      USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
      USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
      USE HCO_STATE_MOD,      ONLY : HCO_STATE
      USE ERROR_MOD,          ONLY : ERROR_STOP

      IMPLICIT NONE
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
! !REVISION HISTORY:
!  20 Aug 2003 - R. Hudman   - Initial version
!  (1 ) Updated from the old routine "chemo3_split.f" (rch, bmy, 8/20/03)
!  (2 ) Bug fix: don't put function call in WRITE statement (bmy, 2/20/04)
!  (3 ) Now use ND44_TMP array to store vertical levels of drydep flux, then
!        sum into AD44 array.  This prevents numerical differences when using
!        multiple processors. (bmy, 3/24/04)
!  (4 ) Now references LDRYD from "logical_mod.f".  Now references STT 
!        and N_TRACERS from "tracer_mod.f".  Now references AD65 from 
!        "diag_pl_mod.f".  Now uses ITS_A_NEW_DAY from "time_mod.f".
!        (bmy, 7/20/04)
!  (5 ) Bug fix: Now avoid a SEG FAULT error if PBLFRAC isn't allocated.
!        (bdf, bmy, 10/12/04)
!  (6 ) Replace PBLFRAC from "drydep_mod.f" with GET_FRAC_UNDER_PBLTOP
!        from "pbl_mix_mod.f".  Now only sum ND44 diagnostic up to the
!        maximum tropopsheric level. (bmy, 2/17/05)
!  (7 ) Resize PP, N D44_TMP arrays from LLTROP to LLTROP_FIX.  Now only loop 
!        up to LLTROP_FIX (phs, 1/19/07) 
!  (8 ) Now use LLTROP instead of LLTROP_FIX (bmy, 12/4/07)
!  (9 ) Now use LD65 instead of LLTROP (phs, 11/17/08)
!  (10) Now only compute loss rate in troposphere (dbj, bmy, 10/26/09)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of NOx-Ox
!                              partitioning
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  23 Jan 2014 - M. Sulprizio- Now use LLTROP instead of LD65 (hyl,bmy,11/3/11)
!  11 Apr 2014 - R. Yantosca - Remove call to INIT_TAGGED_OX
!  24 Jul 2014 - R. Yantosca - Now compute BOXVL internally
!  17 Sep 2014 - C. Keller   - Now use HEMCO for prod. and loss arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE     :: FIRST   = .TRUE.
      INTEGER, SAVE     :: LASTDAY = -1
      INTEGER           :: I, J, L, N
      REAL*8            :: PP(IIPAR,JJPAR,LLTROP,N_TAGGED)
      REAL*8            :: ND44_TMP(IIPAR,JJPAR,LLTROP)
      REAL*8            :: DTCHEM,  DT,      FREQ,    FLUX
      REAL*8            :: LL,      PL,      Ox_0
      REAL*8            :: Ox_LOST, PBL_MAX, F_UNDER_TOP
      REAL*8            :: BOXVL
      CHARACTER(LEN=16) :: STAMP

      ! For fields from Input_Opt
      LOGICAL           :: LDRYD
      LOGICAL           :: LNLPBL
      LOGICAL           :: LSPLIT
      INTEGER           :: N_TRACERS
      REAL*8            :: XNUMOL(Input_Opt%N_TRACERS)

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      ! To get emission timestep
      TYPE(HCO_STATE), POINTER :: HcoState => NULL()

      CHARACTER(LEN=255) :: LOC = 'CHEM_TAGGED_OX (tagged_ox_mod.F)'

      !=================================================================
      ! CHEM_TAGGED_OX begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Get fields from Input_Opt
      LDRYD     = Input_Opt%LDRYD
      LNLPBL    = Input_Opt%LNLPBL
      LSPLIT    = Input_Opt%LSPLIT
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT       => State_Chm%Tracers

      ! Chemistry timestep [s]
      DTCHEM    = GET_TS_CHEM() * 60d0


      ! ----------------------------------------------------------
      ! Get production and loss frequencies from HEMCO. The target
      ! will be updated automatically by HEMCO. 
      ! Important: the file data are converted by HEMCO to HEMCO
      ! concentration units, e.g. the imported data is in kg/m3
      ! (production) and 1/m3 (loss), e.g. the original data in
      ! units of kg/m3/s multiplied by the emission timestep.
      !                                     (ckeller, 9/17/2014)
      ! ----------------------------------------------------------
      IF ( FIRST ) THEN

         ! Get pointer to O3 production
         CALL HCO_GetPtr ( am_I_Root, 'O3_PROD', P24H, RC )
         IF ( RC /= HCO_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to O3_PROD', LOC )

         ! Get pointer to O3 loss
         CALL HCO_GetPtr ( am_I_Root, 'O3_LOSS', L24H, RC )
         IF ( RC /= HCO_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to O3_LOSS', LOC )

         ! Get emissions timestep (seconds)
         CALL GetHcoState ( HcoState )
         IF ( .NOT. ASSOCIATED(HcoState) ) 
     &      CALL ERROR_STOP ( 'Cannot get HcoState', LOC )
         TS_EMIS = HcoState%TS_EMIS

         HcoState => NULL()
         FIRST    = .FALSE.
      ENDIF

      ! Maximum extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L()

      ! DT is the ratio of the chemistry and emission time step.
      ! Use this value to convert from kg/m3 or 1/m3 per emission
      ! time step to kg/m3 or 1/m3 per chemistry time step.
      ! (ckeller, 9/17/2014).
      DT = DTCHEM / TS_EMIS

      !=================================================================
      ! Tagged Ox chemistry contains the following terms:
      !
      !   New Ox = Old Ox - Drydep(Ox) + ( P(Ox,region) - L(Ox) )
      !
      ! P(Ox) and L(Ox) are archived from a previous fullchem run using
      ! the ND20 diagnostic.  P(Ox,region) is the P(Ox) for a specific
      ! tagged Ox tracer, as computed by routine GET_REGIONAL_POX.
      !
      ! Tagged Ox tracers are defined by both geographic location and
      ! altitude, as listed below:
      !
      ! EXTENDED SIMUATION: Modification for quicker spinup
      !=================================================================
      DO N = 1, N_TRACERS

         ! Zero ND44_TMP array
         IF ( ND44 > 0 ) ND44_TMP = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,       L,    LL,          PL,   FREQ )
!$OMP+PRIVATE( Ox_0, Ox_LOST, FLUX, F_UNDER_TOP, BOXVL      )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Grid box volume [m3]
            BOXVL = State_Met%AIRVOL(I,J,L) !* 1d6

            !===========================================================
            ! Get P(Ox) and L(Ox) for each tagged tracer in [kg]
            !===========================================================

            ! P(Ox) is a function of geographic & altitude location
            ! NOTE: We call this only when N==1 for optimal looping
            IF ( N == 1 ) CALL GET_REGIONAL_POX( I, J, L, PP,
     &                                           State_Met   )
            
            ! L(Ox) is originally in [1/cm3/s]; convert to [kg]
            ! L(Ox) is now in [1/m3] (ckeller, 9/17/2014)
            IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN
               LL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL * DT
            ELSE
               LL = 0.0d0
            ENDIF 

            !===========================================================
            ! ND65 diagnostic: Chemical prod/loss [kg/s]
            !===========================================================
            IF ( ND65 > 0 ) THEN

               ! Only archive chemical production if this
               ! region has production to begin with [kg/s]
               IF ( PP(I,J,L,N) > 0d0 ) THEN
                  PL            = P24H(I,J,L) * BOXVL / TS_EMIS
                  AD65(I,J,L,N) = AD65(I,J,L,N) + PL
               ENDIF

               ! Archive loss for all tracers [kg/s]
               PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL / TS_EMIS
               AD65(I,J,L,N_TRACERS+N) = AD65(I,J,L,N_TRACERS+N) + PL

            ENDIF

            !===========================================================
            ! Apply drydep of Ox to each tagged tracer.  We need 
            ! to do this using before P(Ox) - L(Ox) is applied.
            !===========================================================
            IF ( LDRYD ) THEN

               ! Fraction of box underneath the PBL top [unitless]
               F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

               IF ( F_UNDER_TOP > 0d0 ) THEN
               
                  ! Ox Drydep frequency [1/s]
                  FREQ = DEPSAV(I,J,1) * F_UNDER_TOP

                  ! Add option for non-local PBL (Lin, 03/31/09)  
                  IF (LNLPBL) FREQ = 0.D0

                  ! Only proceed if drydep frequency is nonzero
                  IF ( FREQ > 0d0 ) THEN

                     ! Initial Ox [kg]
                     Ox_0    = STT(I,J,L,N)

                     ! Amount of Ox LOST to drydep [kg]
                     Ox_LOST = Ox_0 * ( 1d0 - EXP( -FREQ * DTCHEM ) )
                  
                     ! Prevent underflow condition
                     IF ( Ox_LOST < 1d-20 ) Ox_LOST = 0d0
                       
                     ! Subtract Ox lost [kg] 
                     STT(I,J,L,N) = Ox_0 - Ox_LOST 
                  
                     !==================================================
                     ! ND44 diagnostic: Ox lost to drydep [molec/cm2/s]
                     !==================================================
                     IF ( ND44 > 0 .and. Ox_LOST > 0d0 ) THEN

                        ! Convert from [kg] to [molec/cm2/s]
                        FLUX = Ox_LOST             * XNUMOL(IDTO3) / 
     &                         GET_AREA_CM2(I,J,L) / DTCHEM 
                     
                        ! Store dryd flx in ND44_TMP as a placeholder
                        ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
               
            !===========================================================
            ! After removing Ox lost to dry deposition, apply 
            ! chemical P(Ox) - L(Ox) to each tagged tracer
            !===========================================================
            STT(I,J,L,N) = STT(I,J,L,N) + PP(I,J,L,N) - LL
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! ND44: Sum drydep fluxes by level into the AD44 array in
         ! order to ensure that  we get the same results w/ sp or mp 
         !==============================================================
         IF ( ND44 > 0 ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            DO L = 1, PBL_MAX
               AD44(I,J,N,1) = AD44(I,J,N,1) + ND44_TMP(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_TAGGED_OX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tagged_ox
!
! !DESCRIPTION: Subroutine INIT\_TAGGED\_OX allocates and zeroes all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TAGGED_OX( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is now called from GIGC_INIT_EXTRA in GeosCore/input_mod.F.
! 
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  (1 ) Now reference N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now use LD65 instead of LLTROP to dimension P24H, L24H (phs, 11/18/08)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC args
!  23 Jan 2014 - M. Sulprizio- Now use LLTROP instead of LD65 to dimension P24H,
!                              L24H (hyl, bmy, 11/3/11)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_TAGGED_OX begins here
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Safety valve
      IF ( Input_Opt%N_TRACERS > N_TAGGED ) THEN
         CALL ERROR_STOP( 'NTRACE is too large for Tagged Ox!', 
     &                    'INIT_TAGGED_OX (tagged_ox_mod.f)' )
      ENDIF

      END SUBROUTINE INIT_TAGGED_OX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tagged_ox
!
! !DESCRIPTION:Subroutine CLEANUP\_TAGGED\_OX deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TAGGED_OX
! 
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Deallocate module arrays
      IF ( ASSOCIATED( P24H ) ) P24H => NULL() 
      IF ( ASSOCIATED( L24H ) ) L24H => NULL()

      END SUBROUTINE CLEANUP_TAGGED_OX
!EOC
      END MODULE TAGGED_OX_MOD
