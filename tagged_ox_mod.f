! $Id: tagged_ox_mod.f,v 1.23 2007/12/04 16:23:59 bmy Exp $
      MODULE TAGGED_OX_MOD
!
!******************************************************************************
!  Module TAGGED_OX_MOD contains variables and routines to perform a tagged Ox
!  simulation.  P(Ox) and L(Ox) rates need to be archived from a full chemistry
!  simulation before you can run w/ Tagged Ox. (amf,rch,bmy, 8/20/03, 12/4/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) N_TAGGED (INTEGER) : Total number of tagged tracers
!  (2 ) N_STRAT  (INTEGER) : Denotes tracer # of stratospheric Ox
!  (3 ) N_INIT   (INTEGER) : Denotes tracer # of initial condition Ox
!  (4 ) N_USA    (INTEGER) : Denotes tracer # of USA produced Ox
!  (5 ) P24H     (REAL*8 ) : 24-hr avg P(Ox) saved from fullchem run [kg/cm3/s]
!  (6 ) L24H     (REAL*8 ) : 24-hr avg L(Ox) saved from fullchem run [ 1/cm3/s]
! 
!  Module Routines:
!  ============================================================================
!  (1 ) ADD_STRAT_POX      : Adds strat P(Ox) from UPBDFLX_O3 to tracer array
!  (2 ) READ_POX_LOX       : Reads previously archived P(Ox), L(Ox) from disk
!  (3 ) GET_REGIONAL_POX   : Flags tracers by geographic & vertical location
!  (4 ) CHEM_TAGGED_OX     : Performs Ox chem on geographically tagged tracers
!  (5 ) INIT_TAGGED_OX     : Allocates and zeroes all module arrays
!  (6 ) CLEANUP_TAGGED_OX  : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by tagged_ox_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) diag_pl_mod.f      : Module w/ routines for ND65 & ND20 diagnostics
!  (5 ) error_mod.f        : Module w/ I/O error and NaN check routines
!  (6 ) grid_mod.f         : Module w/ horizontal grid information
!  (7 ) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing 
!  (8 ) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (9 ) time_mod.f         : Module w/ routines for computing time & date
!  (10) transfer_mod.f     : Module w/ routines to cast & resize arrays
!  (11) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tagged_ox_mod.f"
      !=================================================================

      ! PRIVATE module routines  
      PRIVATE :: GET_REGIONAL_POX
      PRIVATE :: INIT_TAGGED_OX
      PRIVATE :: READ_POX_LOX

      ! PRIVATE module variables
      PRIVATE :: N_TAGGED, N_INIT, N_STRAT
      PRIVATE :: N_USA,    P24H,   L24H

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      !------------------------------------------------------
      !%%% Modification for quick Ox spinup (bmy, 5/31/07) 
      !%%%INTEGER, PARAMETER   :: N_TAGGED = 13
      !%%%INTEGER, PARAMETER   :: N_STRAT  = 11
      !%%%INTEGER, PARAMETER   :: N_INIT   = 12
      !%%%INTEGER, PARAMETER   :: N_USA    = 13
      !------------------------------------------------------
      INTEGER, PARAMETER   :: N_TAGGED = 3
      INTEGER, PARAMETER   :: N_STRAT  = 2
      INTEGER, PARAMETER   :: N_INIT   = 3
      INTEGER, PARAMETER   :: N_USA    = -1
      REAL*8,  ALLOCATABLE :: P24H(:,:,:)
      REAL*8,  ALLOCATABLE :: L24H(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE ADD_STRAT_POX( I, J, L, POx )
!
!******************************************************************************
!  Subroutine ADD_STRAT_POX adds the stratospheric influx of Ox to the
!  stratospheric Ox tracer.  This is called from routine UPBDFLX_O3, 
!  which is applied when the tracer array has units of [v/v].  
!  (bmy, 8/19/03, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I,J,L (INTEGER) : GEOS-CHEM grid box indices for lon, lat, alt
!  (4  ) POx   (REAL*8 ) : P(Ox) in the stratosphere [v/v]
!
!  NOTES:
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: POx

      !=================================================================
      ! GET_STRAT_POX begins here!
      !=================================================================
      STT(I,J,L,N_STRAT) = STT(I,J,L,N_STRAT) + POx

      ! Return to calling program
      END SUBROUTINE ADD_STRAT_POX

!------------------------------------------------------------------------------

      SUBROUTINE READ_POX_LOX
!
!******************************************************************************
!  Subroutine READ_POX_LOX reads previously-archived Ox production & loss 
!  rates from binary punch file format. (bmy, 8/20/03, 12/4/07)
! 
!  NOTES:
!  (1 ) Updated from the old routine "chemo3_split.f" (rch, bmy, 8/20/03)
!  (2 ) Now references O3PL_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Use LLTROP_FIX to limit array size to case of non-variable tropopause.
!        Also zero ARRAY to avoid numerical problems (phs, 1/19/07)
!  (5 ) Now use LLTROP instead of LLTROP_FIX (phs, bmy, 12/4/07)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : O3PL_DIR 
      USE TIME_MOD,      ONLY : EXPAND_DATE, GET_NYMD, GET_TAU
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP
           
#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! LD65

      ! Local variables
      REAL*4             :: ARRAY(IGLOB,JGLOB,LLTROP)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_POX_LOX begins here!
      !=================================================================

      ! Filename string
      FILENAME = 'rate.YYYYMMDD'
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), 000000 )

      ! Prefix FILENAME w/ the proper directory
      FILENAME = TRIM( O3PL_DIR ) // FILENAME

      ! Echo information
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_POX_LOX: Reading ', a )

      ! Get the TAU0 value for today
      XTAU = GET_TAU()

      !=================================================================
      ! Read P(O3) [kg/cm3/s]
      !=================================================================

      ! Initialize
      ARRAY = 0e0

      ! Limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
      ! tropopause. This is backward compatibility with offline data set.
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 1,      
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_3D_TROP( ARRAY, P24H )

      !=================================================================
      ! Read L(O3) [1/cm3/s]
      !=================================================================

      ! Initialize
      ARRAY = 0e0

      ! read data
      CALL READ_BPCH2( FILENAME, 'PORL-L=$', 2,      
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LLTROP,    ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAY, L24H )

      ! Return to calling program
      END SUBROUTINE READ_POX_LOX

!------------------------------------------------------------------------------

      SUBROUTINE GET_REGIONAL_POX( I, J, L, PP )
!
!******************************************************************************
!  Subroutine GET_REGIONAL_POX returns the P(Ox) for each of the tagged Ox 
!  tracers. Tagged Ox tracers are defined by both geographic location and 
!  altitude. (amf, rch, bmy, 8/19/03, 12/4/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I,J,L (INTEGER) : GEOS-CHEM grid box indices for lon, lat, alt
!
!  Return Value
!  ============================================================================
!  (4  ) PP    (REAL*8)  : Array containing P(Ox) for each tagged tracer
! 
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : PBL
      USE GRID_MOD,       ONLY : GET_XMID,  GET_YMID
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_GCTM"       ! SCALE_HEIGHT

      ! Arguments
      INTEGER, INTENT(IN)     :: I, J, L
      REAL*8,  INTENT(OUT)    :: PP(IIPAR,JJPAR,LLTROP,N_TAGGED)

      ! Local variables
      LOGICAL                 :: ITS_IN_TROP, ITS_IN_PBL, ITS_IN_MT
      LOGICAL                 :: ITS_IN_UT,   ITS_IN_NH,  ITS_IN_ATL
      LOGICAL                 :: ITS_IN_PAC,  ITS_IN_AS,  ITS_IN_EUR
      LOGICAL                 :: ITS_IN_NAM,  ITS_IN_NAF, ITS_IN_USA
      INTEGER                 :: PBLTOP,      MTTOP
      REAL*8                  :: PPROD,       X,          Y
      
      ! External functions
      REAL*8, EXTERNAL        :: BOXVL

      !=================================================================
      ! GET_REGIONAL_POX begins here!
      !=================================================================

      ! Initialize
      PP(I,J,L,:) = 0d0
      
      ! IS TROP is TRUE if we are in the troposphere
      ITS_IN_TROP = ITS_IN_THE_TROP( I, J, L )
      
      ! Skip stratospheric boxes
      IF ( .not. ITS_IN_TROP ) RETURN

      !=================================================================
      ! Get lon, lat, alt coordinates; test for specific regions
      ! NOTE: You can update the lat/lon & alt boundaries as needed!
      !=================================================================      

      ! Longitude [degrees]
      X          = GET_XMID( I )   
      Y          = GET_YMID( J )

      ! PBLTOP is the model level at ~ 750 hPa
      ! MTTOP  is the model level at ~ 350 hPa
#if   defined( GEOS_3 ) 
      PBLTOP = 10
      MTTOP  = 16
#elif defined( GEOS_4 )
      PBLTOP = 5
      MTTOP  = 10
#elif defined( GEOS_5 )
      PBLTOP = 5
      MTTOP  = 10
#endif

!      ! Define flags for various geographic & altitude regions
!      ITS_IN_PBL = ( L <= PBLTOP                                       )
!      ITS_IN_MT  = ( L >  PBLTOP .and. L <= MTTOP                      )
!      ITS_IN_UT  = ( L >  MTTOP  .and. ITS_IN_TROP                     )
!
!      ITS_IN_NH  = ( Y >=   0.0                                        )
!      ITS_IN_EUR = ( Y >=  36.0 .and. ( X >  -15.0 .and. X >=   55.0 ) )
!      ITS_IN_NAM = ( Y >=  15.0 .and. ( X > -127.5 .and. X <=  -65.0 ) )
!      ITS_IN_AS  = ( Y >= -10.0 .and. ( X >   55.0 .and. X <=  145.0 ) )
!      ITS_IN_ATL = ( ITS_IN_NH  .and. ( X >  -65.0 .and. X <=  -15.0 ) )
!      ITS_IN_PAC = ( ITS_IN_NH  .and. ( X >  145.0  .or. X <= -127.5 ) )
!
!      ITS_IN_NAF = ( ( X >= -15.0 .and. X <=  55.0 ) .and. 
!     &               ( Y >=   0.0 .and. Y <   36.0 ) )  
!
!      ITS_IN_USA = ( ( X > -127.5 .and. X <= -65.0 ) .and. 
!     &               ( Y >   22.0 .and. Y <=  50.0 ) )

      !=================================================================
      ! Assign P(Ox) to tagged tracers by geographic/altitude regions
      !=================================================================

      ! P(Ox) [kg]
      PPROD = P24H(I,J,L) * BOXVL(I,J,L) * ( GET_TS_CHEM() * 60d0 )

      !-----------------------
      ! #1: Total P(Ox)
      !-----------------------
      PP(I,J,L,1) = PPROD

!      !-----------------------
!      ! #2: P(Ox) in UT
!      !-----------------------
!      IF ( ITS_IN_UT ) THEN
!         PP(I,J,L,2) = PPROD
!         
!      !-----------------------
!      ! #3: P(Ox) in MT 
!      !-----------------------
!      ELSE IF ( ITS_IN_MT ) THEN
!         PP(I,J,L,3) = PPROD
!                                
!      !-----------------------
!      ! #5: P(Ox) in Pac BL
!      !-----------------------
!      ELSE IF ( ITS_IN_PAC .and. ITS_IN_PBL ) THEN
!         PP(I,J,L,5) = PPROD
!
!      !-----------------------
!      ! #6: P(Ox) in NAm BL
!      !-----------------------
!      ELSE IF ( ITS_IN_NAM .and. ITS_IN_PBL ) THEN     
!         PP(I,J,L,6) = PPROD
!                  
!      !-----------------------
!      ! #7: P(Ox) in Atl BL
!      !-----------------------
!      ELSE IF ( ITS_IN_ATL .and. ITS_IN_PBL ) THEN
!         PP(I,J,L,7) = PPROD  
!         
!      !-----------------------
!      ! #8: P(Ox) in Eur BL
!      !-----------------------
!      ELSE IF ( ITS_IN_EUR .and. ITS_IN_PBL ) THEN
!         PP(I,J,L,8) = PPROD
!                  
!      !-----------------------
!      ! #9: P(Ox) in NAfr BL
!      !-----------------------
!      ELSE IF ( ITS_IN_NAF .and. ITS_IN_PBL ) THEN
!         PP(I,J,L,9) = PPROD
! 
!      !-----------------------
!      ! #10: P(Ox) in Asia BL
!      !-----------------------          
!      ELSE IF ( ITS_IN_AS .and. ITS_IN_PBL ) THEN
!         PP(I,J,L,10) = PPROD                   
!
!      !-----------------------
!      ! #4: P(Ox) in R.O.W
!      !-----------------------
!      ELSE 
!         PP(I,J,L,4) = PPROD
!
!      ENDIF
!
!      !-------------------------
!      ! #13: P(Ox) in USA
!      !-------------------------
!      IF ( ITS_IN_USA ) THEN
!         PP(I,J,L,N_USA) = PPROD               
!      ENDIF

      ! Return to calling program
      END SUBROUTINE GET_REGIONAL_POX

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_TAGGED_OX
!
!******************************************************************************
!  Subroutine CHEM_TAGGED_OX performs chemistry for several Ox tracers which
!  are tagged by geographic and altitude regions. (rch, bmy, 8/20/03, 12/4/07)
! 
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE DIAG_PL_MOD,  ONLY : AD65
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP
      USE DRYDEP_MOD,   ONLY : DEPSAV      
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LDRYD
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM,           ITS_A_NEW_DAY 
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRACER_MOD,   ONLY : STT,                   N_TRACERS, XNUMOL
      USE TRACERID_MOD, ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND44, ND65, LD65

      ! Local variables
      LOGICAL, SAVE     :: FIRST   = .TRUE.
      INTEGER, SAVE     :: LASTDAY = -1
      INTEGER           :: I, J, L, N
      REAL*8            :: PP(IIPAR,JJPAR,LLTROP,N_TAGGED)
      REAL*8            :: ND44_TMP(IIPAR,JJPAR,LLTROP)
      REAL*8            :: DTCHEM,  FREQ,    FLUX
      REAL*8            :: LL,      PL,      Ox_0
      REAL*8            :: Ox_LOST, PBL_MAX, F_UNDER_TOP
      CHARACTER(LEN=16) :: STAMP

      ! External routines
      REAL*8, EXTERNAL :: BOXVL

      !=================================================================
      ! CHEM_TAGGED_OX begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! First-time initialization only
      IF ( FIRST ) THEN 
         CALL INIT_TAGGED_OX
         FIRST = .FALSE.
      ENDIF

      ! Read P(Ox) and L(Ox) if it's a new day
      IF ( ITS_A_NEW_DAY() ) THEN
         CALL READ_POX_LOX
      ENDIF

      ! Maximum extent of the PBL [model layers]
      PBL_MAX = GET_PBL_MAX_L()

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
      !%%% Modification for quicker spinup
      !%%% (1 ) Total Ox
      !%%% (2 ) Ox produced in Upper Trop     (350 hPa    - tropopause) 
      !%%% (3 ) Ox produced in Middle Trop    (PBL top    - 350 hPa   )
      !%%% (4 ) Ox produced in Rest of World  (surface    - PBL top   )
      !%%% (5 ) Ox produced in Pacific BL     (surface    - PBL top   )
      !%%% (6 ) Ox produced in N. American BL (surface    - PBL top   )
      !%%% (7 ) Ox produced in Atlantic BL    (surface    - PBL top   )
      !%%% (8 ) Ox produced in European BL    (surface    - PBL top   )
      !%%% (9 ) Ox produced in N. African BL  (surface    - PBL top   )
      !%%% (10) Ox produced in Asian          (surface    - PBL top   )
      !%%% (11) Ox from the Stratosphere      (tropopause - atm top   )
      !%%% (12) Ox initial conditions         (all levels             )        
      !%%% (13) Ox produced over the USA      (all levels             )
      !
      ! NOTE: MODIFIED FOR QUICKER SPINUP:
      ! (1 ) Total Ox
      ! (2 ) Ox from the Stratosphere      (tropopause - atm top   )
      ! (3 ) Ox initial conditions         (all levels             )        
      !=================================================================
      DO N = 1, N_TRACERS

         ! Zero ND44_TMP array
         IF ( ND44 > 0 ) ND44_TMP = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, LL, PL, FREQ, Ox_0, Ox_LOST, FLUX, F_UNDER_TOP )  
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLTROP !_FIX
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !===========================================================
            ! Get P(Ox) and L(Ox) for each tagged tracer in [kg]
            !===========================================================

            ! P(Ox) is a function of geographic & altitude location
            ! NOTE: We call this only when N==1 for optimal looping
            IF ( N == 1 ) CALL GET_REGIONAL_POX( I, J, L, PP )
            
            ! L(Ox) is originally in [1/cm3/s]; convert to [kg] 
            LL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L) * DTCHEM 
                           
            !===========================================================
            ! ND65 diagnostic: Chemical prod/loss [kg/s]
            !===========================================================
            IF ( ND65 > 0 ) THEN

               ! Only archive chemical production if this
               ! region has production to begin with [kg/s]
               IF ( PP(I,J,L,N) > 0d0 ) THEN
                  PL            = P24H(I,J,L) * BOXVL(I,J,L)
                  AD65(I,J,L,N) = AD65(I,J,L,N) + PL
               ENDIF

               ! Archive loss for all tracers [kg/s]
               PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L)
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
                        FLUX = Ox_LOST         * XNUMOL(IDTOX) / 
     &                         GET_AREA_CM2(J) / DTCHEM 
                     
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

      ! Return to calling program
      END SUBROUTINE CHEM_TAGGED_OX

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TAGGED_OX
!
!******************************************************************************
!  Subroutine INIT_TAGGED_OX allocates and zeroes all module arrays.
!  (bmy, 8/20/03) 
!
!  NOTES:
!  (1 ) Now reference N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : ALLOC_ERR, ERROR_STOP
      USE TRACER_MOD, ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_TAGGED_OX begins here
      !=================================================================

      ! Safety valve
      IF ( N_TRACERS > N_TAGGED ) THEN
         CALL ERROR_STOP( 'NTRACE is too large for Tagged Ox!', 
     &                    'INIT_TAGGED_OX (tagged_ox_mod.f)' )
      ENDIF

      ! Allocate P24H
      ALLOCATE( P24H( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'P24H' )
      P24H = 0d0

      ! Allocate L24H
      ALLOCATE( L24H( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'L24H' ) 
      L24H = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_TAGGED_OX

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TAGGED_OX
!
!******************************************************************************
!  Subroutine CLEANUP_TAGGED_OX deallocates all module arrays (bmy, 8/20/03)
!
!  NOTES:
!******************************************************************************
!
      ! Deallocate module arrays
      IF ( ALLOCATED( P24H ) ) DEALLOCATE( P24H )
      IF ( ALLOCATED( L24H ) ) DEALLOCATE( L24H )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TAGGED_OX

!------------------------------------------------------------------------------

      END MODULE TAGGED_OX_MOD
