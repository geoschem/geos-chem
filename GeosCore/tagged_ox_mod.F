!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: ADD_STRAT_POX 
      PUBLIC  :: CHEM_TAGGED_OX 
      PUBLIC  :: CLEANUP_TAGGED_OX 
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: GET_REGIONAL_POX
      PRIVATE :: INIT_TAGGED_OX
      PRIVATE :: READ_POX_LOX
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      REAL*8,  ALLOCATABLE, PRIVATE :: P24H(:,:,:)
      REAL*8,  ALLOCATABLE, PRIVATE :: L24H(:,:,:)
!
! !DEFINED PARAMETERS:
! 
      !------------------------------------------------------
      !%%% Modification for quick Ox spinup (bmy, 5/31/07) 
      !%%%INTEGER, PARAMETER, PRIVATE :: N_TAGGED = 13
      !%%%INTEGER, PARAMETER, PRIVATE :: N_STRAT  = 11
      !%%%INTEGER, PARAMETER, PRIVATE :: N_INIT   = 12
      !%%%INTEGER, PARAMETER, PRIVATE :: N_USA    = 13
      !------------------------------------------------------
      ! dbj: stratospheric tracers must be tracer 2
      !------------------------------------------------------
      INTEGER, PARAMETER,   PRIVATE :: N_TAGGED = 2
      INTEGER, PARAMETER,   PRIVATE :: N_STRAT  = 2
      INTEGER, PARAMETER,   PRIVATE :: N_INIT   = 3
      INTEGER, PARAMETER,   PRIVATE :: N_USA    = -1

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE ADD_STRAT_POX( I, J, L, POx )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTOxStrt

      USE CMN_SIZE_MOD           ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I     ! GEOS-Chem grid box lon index
      INTEGER, INTENT(IN) :: J     ! GEOS-Chem grid box lat index
      INTEGER, INTENT(IN) :: L     ! GEOS-Chem grid box level index
      REAL*8,  INTENT(IN) :: POx   ! P(Ox) in the stratosphere [v/v]
! 
! !REVISION HISTORY:
!  19 Aug 2003 - R. Yantosca - Initial version
!  (1 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  10 Jun 2012 - L. Murray - Enable tagged Ox for Synoz.
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_STRAT_POX begins here!
      !=================================================================
      STT(I,J,L,IDTOxStrt) = STT(I,J,L,IDTOxStrt) + POx

      END SUBROUTINE ADD_STRAT_POX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_pox_lox
!
! !DESCRIPTION: Subroutine READ\_POX\_LOX reads previously-archived Ox 
!  production and loss rates from binary punch file format. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_POX_LOX
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : O3PL_DIR 
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TIME_MOD,      ONLY : GET_NYMD
      USE TIME_MOD,      ONLY : GET_TAU
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP
      ! JLIU,2008/10/01
      USE CHARPAK_MOD,   ONLY : STRREPL
      USE TIME_MOD,      ONLY : YMD_EXTRACT
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR
      USE TIME_MOD,      ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,      ONLY : GET_YEAR
      USE TIME_MOD,      ONLY : GET_HOUR
      USE DIAG_PL_MOD                                      !dbj
      USE JULDAY_MOD,    ONLY : JULDAY                     !dbj
           
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_DIAG_MOD ! LD65
! 
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  (1 ) Updated from the old routine "chemo3_split.f" (rch, bmy, 8/20/03)
!  (2 ) Now references O3PL_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Use LLTROP_FIX to limit array size to case of non-variable tropopause.
!        Also zero ARRAY to avoid numerical problems (phs, 1/19/07)
!  (5 ) Now use LLTROP instead of LLTROP_FIX (phs, bmy, 12/4/07)
!  (6 ) Now use LD65, since this is the number of levels use to 
!        save diag20 (phs, 11/17/08)
!  (7 ) Updates for LINOZ (dbj, jliu, bmy, 10/16/09)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4             :: ARRAY(IIPAR,JJPAR,LD65)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
      !JLIU, 2008/10/01
      INTEGER            :: YYYY, MM, DD, FIRST_DATE
      INTEGER            :: DAY_OF_YEAR
      CHARACTER(LEN=2)   :: MM_STR, DD_STR

      !=================================================================
      ! READ_POX_LOX begins here!
      !=================================================================

      ! Filename string
      FILENAME = 'rate.YYYYMMDD'
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), 000000 )

      ! Get the TAU0 value for today
      XTAU = GET_TAU()

!------------------------------------------------------------------------------
! Comment out for now (bmy, 2/4/11)
!      ! -----------------------------------------
!      ! dbj 
!      ! JLIU, 2008/10/01
!      ! -----------------------------------------
!      IF (TAGO3_PL_YEAR .GT. 0) THEN
!
!          ! Extract today's date into year, month, and day sections
!          CALL YMD_EXTRACT( GET_NYMD(), YYYY, MM, DD )
!
!          IF ( ITS_A_LEAPYEAR() ) THEN
!             IF ( (.NOT. ITS_A_LEAPYEAR(TAGO3_PL_YEAR)) .AND. 
!     &            (DD .EQ. 29) ) THEN
!                 DD = DD - 1
!             ENDIF
!          ENDIF
!
!#if       defined( LINUX_PGI )
!          ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
!          ENCODE( 2, '(a,i4.4,i2.2,i2.2)', FILENAME ) 
!     &        'rate.',TAGO3_PL_YEAR,MM,DD
!#else
!          ! For other platforms, use an F90 internal write (bmy, 9/29/03)
!          WRITE(FILENAME, '(a,i4.4,i2.2,i2.2)') 
!     &         'rate.',TAGO3_PL_YEAR,MM,DD
!#endif
!
!      ELSE
!          FILENAME = 'rate.YYYYMMDD'
!          CALL EXPAND_DATE( FILENAME, GET_NYMD(), 000000 )
!
!      ENDIF
!------------------------------------------------------------------------------

      ! Prefix FILENAME w/ the proper directory
      FILENAME = TRIM( O3PL_DIR ) // FILENAME

      ! Echo information
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_POX_LOX: Reading ', a )

!------------------------------------------------------------------------------
! Comment out for now (bmy, 2/4/11)
!      !---------------------------------------
!      ! dbj
!      !JLIU, 2008/10/01
!      !---------------------------------------
!
!      IF (TAGO3_PL_YEAR .GT. 0) THEN
!
!         XTAU =
!     &          ( (JULDAY(TAGO3_PL_YEAR,MM,DFLOAT(DD)) - 
!     &             JULDAY(1985,1,DFLOAT(1)))*24.0d0 ) +
!     &             GET_HOUR()
!
!
!      ELSE
!         ! Get the TAU0 value for today
!         XTAU = GET_TAU()
!      ENDIF
!------------------------------------------------------------------------------

      !=================================================================
      ! Read P(O3) [kg/cm3/s]
      !=================================================================

      ! Initialize
      ARRAY = 0e0

      ! Limit array 3d dimension to LLTROP_FIX, i.e, case of annual mean
      ! tropopause. This is backward compatibility with offline data set.
      CALL READ_BPCH2( FILENAME,  'PORL-L=$', 1,      
     &                 XTAU,       IIPAR,     JJPAR,      
     &                 LLTROP_FIX, ARRAY,     QUIET=.false.)

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_3D_TROP( ARRAY, P24H )

      !=================================================================
      ! Read L(O3) [1/cm3/s]
      !=================================================================

      ! Initialize
      ARRAY = 0e0

      ! read data
      CALL READ_BPCH2( FILENAME,  'PORL-L=$', 2,      
     &                 XTAU,       IIPAR,     JJPAR,      
     &                 LLTROP_FIX, ARRAY,     QUIET=.false.)

      ! Cast from REAL*4 to REAL*8 
      CALL TRANSFER_3D_TROP( ARRAY, L24H )

      END SUBROUTINE READ_POX_LOX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
#if defined( DEVEL )
      SUBROUTINE GET_REGIONAL_POX( I, J, L, PP, LOCAL_MET )
#else
      SUBROUTINE GET_REGIONAL_POX( I, J, L, PP )
#endif
!
! !USES:
!
      USE GRID_MOD,       ONLY : GET_XMID,  GET_YMID
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
#if defined( DEVEL )
      USE GC_TYPE_MOD,    ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_DIAG_MOD       ! ND44, ND65, LD65
      USE CMN_GCTM_MOD       ! SCALE_HEIGHT
!
! !INPUT PARAMETERS: 
!
      ! GEOS-Chem grid box indices for lon, lat, alt
      INTEGER, INTENT(IN)  :: I, J, L
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
!
! !OUTPUT PARAMETERS:
!
      ! Array containing P(Ox) for each tagged tracer
      REAL*8,  INTENT(OUT) :: PP(IIPAR,JJPAR,LD65,N_TAGGED)
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL          :: ITS_IN_TROP, ITS_IN_PBL, ITS_IN_MT
      LOGICAL          :: ITS_IN_UT,   ITS_IN_NH,  ITS_IN_ATL
      LOGICAL          :: ITS_IN_PAC,  ITS_IN_AS,  ITS_IN_EUR
      LOGICAL          :: ITS_IN_NAM,  ITS_IN_NAF, ITS_IN_USA
      INTEGER          :: PBLTOP,      MTTOP
      REAL*8           :: PPROD,       X,          Y
      
      ! External functions
      REAL*8, EXTERNAL :: BOXVL

      !=================================================================
      ! GET_REGIONAL_POX begins here!
      !=================================================================

      ! Initialize
      PP(I,J,L,:) = 0d0
      
      ! IS TROP is TRUE if we are in the troposphere
#if defined( DEVEL )
      ITS_IN_TROP = ITS_IN_THE_TROP( I, J, L, LOCAL_MET )
#else
      ITS_IN_TROP = ITS_IN_THE_TROP( I, J, L )
#endif
      
      ! Skip stratospheric boxes
      IF ( .not. ITS_IN_TROP ) RETURN

      !=================================================================
      ! Get lon, lat, alt coordinates; test for specific regions
      ! NOTE: You can update the lat/lon & alt boundaries as needed!
      !=================================================================      

      ! Longitude [degrees]
      X          = GET_XMID( I, J, L )   
      Y          = GET_YMID( I, J, L )

      ! PBLTOP is the model level at ~ 750 hPa
      ! MTTOP  is the model level at ~ 350 hPa
#if   defined( GEOS_4 )
      PBLTOP = 5
      MTTOP  = 10
#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
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
#if defined( DEVEL )
      PPROD = P24H(I,J,L)     * BOXVL(I,J,L,LOCAL_MET) *
     &        ( GET_TS_CHEM() * 60d0 )
#else
      PPROD = P24H(I,J,L) * BOXVL(I,J,L) * ( GET_TS_CHEM() * 60d0 )
#endif

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

      END SUBROUTINE GET_REGIONAL_POX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
#if defined( DEVEL )
      SUBROUTINE CHEM_TAGGED_OX( LOCAL_MET )
#else
      SUBROUTINE CHEM_TAGGED_OX 
#endif
!
! !USES:
!
      USE DIAG_MOD,       ONLY : AD44
      USE DIAG_PL_MOD,    ONLY : AD65
      USE ERROR_MOD,      ONLY : GEOS_CHEM_STOP
      USE DRYDEP_MOD,     ONLY : DEPSAV      
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,    ONLY : LDRYD
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
      USE PBL_MIX_MOD,    ONLY : GET_PBL_MAX_L
      USE TIME_MOD,       ONLY : GET_TS_CHEM
      USE TIME_MOD,       ONLY : ITS_A_NEW_DAY 
      USE TIME_MOD,       ONLY : TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : STT
      USE TRACER_MOD,     ONLY : N_TRACERS
      USE TRACER_MOD,     ONLY : XNUMOL
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE LOGICAL_MOD,    ONLY : LNLPBL
#if defined( DEVEL )
      USE GC_TYPE_MOD,  ONLY   : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_DIAG_MOD  ! ND44, ND65, LD65

      IMPLICIT NONE
#     include "define.h"

!
! !INPUT PARAMETERS: 
!
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE     :: FIRST   = .TRUE.
      INTEGER, SAVE     :: LASTDAY = -1
      INTEGER           :: I, J, L, N
      REAL*8            :: PP(IIPAR,JJPAR,LD65,N_TAGGED)
      REAL*8            :: ND44_TMP(IIPAR,JJPAR,LD65)
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
         DO L = 1, LD65 ! LLTROP !_FIX
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !===========================================================
            ! Get P(Ox) and L(Ox) for each tagged tracer in [kg]
            !===========================================================

#if defined( DEVEL )
            ! P(Ox) is a function of geographic & altitude location
            ! NOTE: We call this only when N==1 for optimal looping
            IF ( N == 1 ) CALL GET_REGIONAL_POX( I, J, L, PP,
     &                                           LOCAL_MET   )
            
            ! L(Ox) is originally in [1/cm3/s]; convert to [kg] 
            IF ( ITS_IN_THE_TROP( I, J, L, LOCAL_MET ) ) THEN
               LL = STT(I,J,L,N)           * L24H(I,J,L) *
     &              BOXVL(I,J,L,LOCAL_MET) * DTCHEM
            ELSE
               LL = 0.0d0
            ENDIF 
#else
            ! P(Ox) is a function of geographic & altitude location
            ! NOTE: We call this only when N==1 for optimal looping
            IF ( N == 1 ) CALL GET_REGIONAL_POX( I, J, L, PP )
            
            ! L(Ox) is originally in [1/cm3/s]; convert to [kg] 
            IF ( ITS_IN_THE_TROP( I, J, L ) ) THEN
               LL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L) * DTCHEM
            ELSE
               LL = 0.0d0
            ENDIF 
#endif
                           
            !===========================================================
            ! ND65 diagnostic: Chemical prod/loss [kg/s]
            !===========================================================
            IF ( ND65 > 0 ) THEN

#if defined( DEVEL )
               ! Only archive chemical production if this
               ! region has production to begin with [kg/s]
               IF ( PP(I,J,L,N) > 0d0 ) THEN
                  PL            = P24H(I,J,L) * BOXVL(I,J,L,LOCAL_MET)
                  AD65(I,J,L,N) = AD65(I,J,L,N) + PL
               ENDIF

               ! Archive loss for all tracers [kg/s]
               PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L,LOCAL_MET)
               AD65(I,J,L,N_TRACERS+N) = AD65(I,J,L,N_TRACERS+N) + PL
#else
               ! Only archive chemical production if this
               ! region has production to begin with [kg/s]
               IF ( PP(I,J,L,N) > 0d0 ) THEN
                  PL            = P24H(I,J,L) * BOXVL(I,J,L)
                  AD65(I,J,L,N) = AD65(I,J,L,N) + PL
               ENDIF

               ! Archive loss for all tracers [kg/s]
               PL = STT(I,J,L,N) * L24H(I,J,L) * BOXVL(I,J,L)
               AD65(I,J,L,N_TRACERS+N) = AD65(I,J,L,N_TRACERS+N) + PL
#endif
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
                        FLUX = Ox_LOST             * XNUMOL(IDTOX) / 
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

      END SUBROUTINE CHEM_TAGGED_OX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE INIT_TAGGED_OX
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ALLOC_ERR
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE TRACER_MOD, ONLY : N_TRACERS

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! ND44, ND65, LD65
! 
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  (1 ) Now reference N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now use LD65 instead of LLTROP to dimension P24H, L24H (phs, 11/18/08)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
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

      ! Safety valve
      IF ( N_TRACERS > N_TAGGED ) THEN
         CALL ERROR_STOP( 'NTRACE is too large for Tagged Ox!', 
     &                    'INIT_TAGGED_OX (tagged_ox_mod.f)' )
      ENDIF

      ! Allocate P24H
      ALLOCATE( P24H( IIPAR, JJPAR, LD65 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'P24H' )
      P24H = 0d0

      ! Allocate L24H
      ALLOCATE( L24H( IIPAR, JJPAR, LD65 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'L24H' ) 
      L24H = 0d0

      END SUBROUTINE INIT_TAGGED_OX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      IF ( ALLOCATED( P24H ) ) DEALLOCATE( P24H )
      IF ( ALLOCATED( L24H ) ) DEALLOCATE( L24H )

      END SUBROUTINE CLEANUP_TAGGED_OX
!EOC
      END MODULE TAGGED_OX_MOD
