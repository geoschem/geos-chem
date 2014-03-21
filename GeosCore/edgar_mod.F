! $Id: edgar_mod.f,v 1.1 2009/09/16 14:06:33 bmy Exp $
      MODULE EDGAR_MOD
!
!******************************************************************************
!  Module EDGAR_MOD contains variables and routines to read anthropogenic 
!  emissions from the EDGAR inventory for NOx, CO and SO2. 
!  (avd, bmy, phs, 7/14/06, 12/20/10)
!
!  Module Routines:
!  ============================================================================
!  (1 ) EMISS_EDGAR           : Driver program for EDGAR emissions
!  (2 ) COMPUTE_EDGAR_NOx     : Computes EDGAR NOx emissions 
!  (3 ) READ_EDGAR_NOx        : Reads EDGAR NOx data from disk
!  (4 ) SEASCL_EDGAR_NOx      : Applies seasonal scale factors to anthro NOx
!  (5 ) COMPUTE_EDGAR_CO      : Computes EDGAR CO emissions 
!  (6 ) READ_EDGAR_CO         : Reads EDGAR CO data from disk
!  (7 ) COMPUTE_EDGAR_SO2     : Computes EDGAR SO2 emissions 
!  (8 ) READ_EDGAR_SO2        : Reads EDGAR SO2 data from disk
!  (9 ) SEASCL_EDGAR_ANTH_SO2 : Applies seasonal scale factors to anthro SO2
!  (10) SEASCL_EDGAR_SHIP_SO2 : Applies monthy scale factors to ship SO2
!  (11) READ_EDGAR_DATA       : Reads an EDGAR data file for a given sector
!  (12) ADD_EDGAR_DATA        : Sums EDGAR data for several sectors
!  (13) GET_EDGAR_NOx         : Returns NOx      emissions at grid box (I,J)
!  (14) GET_EDGAR_CO          : Returns CO       emissions at grid box (I,J)
!  (15) GET_EDGAR_ANTH_SO2    : Returns SOx anth emissions at grid box (I,J)
!  (16) GET_EDGAR_SHIP_SO2    : Returns SOx ship emissions at grid box (I,J)
!  (17) GET_EDGAR_TODN        : Returns NOx time-of-day scale factors at (I,J)
!  (18) READ_AROMATICS        : Reads EDGAR aromatics emissions from disk
!  (19) READ_C2H4             : Reads EDGAR C2H4 emissions from disk
!  (20) READ_C2H2             : Reads EDGAR C2H2 emissions from disk
!  (21) READ_AROMATICS_05x0666: Reads EDGAR aromatics emissions from disk
!                                for 0.5x0.666 resolution
!  (22) READ_C2H4_05x0666     : Reads EDGAR C2H4 emissions from disk
!                                for 0.5x0.666 resolution
!  (23) READ_C2H2_05x0666     : Reads EDGAR C2H2 emissions from disk
!                                for 0.5x0.666 resolution
!  (24) INIT_EDGAR            : Allocates and zeroes module arrays
!  (25) CLEANUP_EDGAR         : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "edgar_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f           : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f       : Module w/ GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f           : Module w/ I/O error and NaN check routines
!  (4 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (5 ) future_emissions_mod.f: Module w/ routines for IPCC future emissions
!  (6 ) grid_mod.f            : Module w/ horizontal grid information
!  (7 ) logical_mod.f         : Module w/ GEOS-CHEM logical switches!  
!  (8 ) regrid_a2a_mod.f      : Module w/ routines to regrid data
!  (9 ) time_mod.f            : Module w/ routines for computing time and date
!  (10) scale_anthro_mod.f    : Module w/ routines to get annual scale factor
!
!  NOTES:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  (2 ) Replaced READ_EDGAR_SCALE with GET_ANNUAL_SCALAR and
!        GET_ANNAUL_SCALAR_1x1 from SCALE_ANTHRO_MOD (amv, phs, 3/10/08)
!  (3 ) Added EDGAR_NOx_SHIP and EDGAR_CO_SHIP variables. Added SHIP flag
!        to GET_EDGAR_NOx and GET_EDGAR_CO, so they can be accessed from
!        outside the module (5/13/08, phs)
!  (4 ) Added subroutines READ_AROMATICS, READ_C2H4 and READ_C2H2 for GLYX
!        chemistry. They are PUBLIC routines as they are used by anthroems.f
!        (ccc, 3/2/09)
!  (5 ) Added subroutines READ_AROMATICS_05x0666, READ_C2H4_05x0666 and 
!        READ_C2H2_05x0666 for 0.5x0.666 resolution emissions. 
!        They are PUBLIC routines as they are used by anthroems.f.(ccc, 3/2/09)
!  (6 ) Fixed typos in READ_C2H2 and READ_C2H4 (tmf, bmy, 12/20/10)
!  01 Mar 2012 - R. Yantosca - Now reference the new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "edgar_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_EDGAR
      PUBLIC :: EMISS_EDGAR
      PUBLIC :: GET_EDGAR_CO
      PUBLIC :: GET_EDGAR_NOx
      PUBLIC :: GET_EDGAR_ANTH_SO2
      PUBLIC :: GET_EDGAR_SHIP_SO2
      PUBLIC :: GET_EDGAR_TODN
      PUBLIC :: READ_AROMATICS
      PUBLIC :: READ_C2H4
      PUBLIC :: READ_C2H2
      PUBLIC :: READ_AROMATICS_05x0666
      PUBLIC :: READ_C2H4_05x0666
      PUBLIC :: READ_C2H2_05x0666

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Variables
      REAL*8               :: SEC_IN_SEASON      
      REAL*8               :: SEC_IN_MONTH
      REAL*8               :: SEASON_TAU0
      REAL*8               :: MONTH_TAU0
      CHARACTER(LEN=3)     :: SEASON_NAME
      CHARACTER(LEN=3)     :: MONTH_NAME

      ! Parameters
      INTEGER, PARAMETER   :: N_HOURS     = 24
      REAL*8,  PARAMETER   :: SEC_IN_DAY  = 86400d0
      REAL*8,  PARAMETER   :: SEC_IN_2000 = SEC_IN_DAY * 366d0  ! leapyear
      REAL*8,  PARAMETER   :: XNUMOL_NO2  = 6.0225d23  / 46d-3
      REAL*8,  PARAMETER   :: XNUMOL_CO   = 6.0225d23  / 28d-3
      REAL*8,  PARAMETER   :: XNUMOL_SO2  = 6.0225d23  / 64d-3

      ! Arrays
      REAL*8,  ALLOCATABLE :: EDGAR_NOx(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_CO(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_SO2(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_SO2_SHIP(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_NOx_SHIP(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_CO_SHIP(:,:)
      REAL*8,  ALLOCATABLE :: EDGAR_TODN(:,:,:)

      ! To shadow logical values from Input_Opt
      LOGICAL              :: LEDGARCO
      LOGICAL              :: LEDGARNOx
      LOGICAL              :: LEDGARSHIP
      LOGICAL              :: LEDGARSOx
      LOGICAL              :: LFUTURE
      LOGICAL              :: LHIST

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_EDGAR( am_I_Root, YEAR,      MONTH,
     &                        Input_Opt, State_Chm, RC     )
!
!******************************************************************************
!  Subroutine EMISS_EDGAR fills emission arrays with emissions based upon 
!  the EDGAR inventory (avd, bmy, 7/14/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR  (INTEGER) : Compute EDGAR emissions for this year ...
!  (2 ) MONTH (INTEGER) :  ... and for this month
!
!  NOTES:
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!******************************************************************************
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState

#if defined( DEVEL )
      USE TRACERID_MOD,   ONLY : IDTCO, IDTSO2, IDTNO
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: YEAR        ! Year to read
      INTEGER,        INTENT(IN)    :: MONTH       ! Month to read
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
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE           :: FIRST      = .TRUE.
      INTEGER, SAVE           :: YEAR_SAVE  = -1
      INTEGER, SAVE           :: MONTH_SAVE = -1
      REAL*8                  :: E_NOX(IIPAR,JJPAR)
      REAL*8                  :: E_NOX_1x1(I1x1,J1x1)
      REAL*8                  :: E_NOX_HRLY_1x1(I1x1,J1x1)
      REAL*8                  :: TEMP(I1x1,J1x1,1)

      ! Days of year 2000
      INTEGER                 :: M(12) = (/ 31, 29, 31, 30, 31, 30,
     &                                      31, 31, 30, 31, 30, 31 /)        

      ! Month names
      CHARACTER(LEN=3)        :: MON(12) = (/ 'JAN', 'FEB', 'MAR', 
     &                                        'APR', 'MAY', 'JUN', 
     &                                        'JUL', 'AUG', 'SEP', 
     &                                        'OCT', 'NOV', 'DEC' /)

#if defined( DEVEL )
      REAL*8  :: SO2_FUTURE_SCALE(IIPAR,JJPAR)
      REAL*8  :: NOX_FUTURE_SCALE(IIPAR,JJPAR)
      REAL*8  :: CO_FUTURE_SCALE(IIPAR,JJPAR)
      INTEGER :: I, J
#endif

      ! For fields from Input_Opt
      INTEGER :: N_TRACERS

      !=================================================================
      ! EMISS_EDGAR begins here!
      !=================================================================

      ! Assume success
      RC         =  GIGC_SUCCESS

      ! Copy values from Input_Opt to module vars for use in other routines
      LEDGARCO   = Input_Opt%LEDGARCO
      LEDGARNOx  = Input_Opt%LEDGARNOx
      LEDGARSHIP = Input_Opt%LEDGARSHIP
      LEDGARSOx  = Input_Opt%LEDGARSOx
      LHIST      = Input_Opt%LHIST
      N_TRACERS  = Input_Opt%N_TRACERS

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_EDGAR( Input_Opt )
         FIRST = .FALSE.
      ENDIF

      ! Get variables for this month
      SEC_IN_MONTH = M(MONTH) * SEC_IN_DAY
      MONTH_NAME   = MON(MONTH)
      MONTH_TAU0   = GET_TAU0( MONTH, 1, 1985 )

      ! Get variables for this season
      SELECT CASE( MONTH )
         CASE( 12,  1,  2 )
            SEC_IN_SEASON = ( M(12) + M(1 ) + M(2 ) ) * SEC_IN_DAY 
            SEASON_NAME   = 'DJF'
            SEASON_TAU0   = -744d0
         CASE(  3,  4,  5 ) 
            SEC_IN_SEASON = ( M(3 ) + M(4 ) + M(5 ) ) * SEC_IN_DAY
            SEASON_NAME   = 'MAM'
            SEASON_TAU0   = 1416d0
         CASE(  6,  7,  8 )
            SEC_IN_SEASON = ( M(6 ) + M(7 ) + M(8 ) ) * SEC_IN_DAY
            SEASON_NAME   = 'JJA'
            SEASON_TAU0   = 3624d0
         CASE(  9, 10, 11 )
            SEC_IN_SEASON = ( M(9 ) + M(10) + M(11) ) * SEC_IN_DAY
            SEASON_NAME   = 'SON'
            SEASON_TAU0   = 5832d0
      END SELECT

      !=================================================================
      ! CO emissions are read once per year
      !=================================================================
      IF ( YEAR /= YEAR_SAVE ) THEN

         ! CO
         IF ( LEDGARCO ) THEN
            IF ( LEDGARSHIP) THEN
               CALL COMPUTE_EDGAR_CO( YEAR, EDGAR_CO, EDGAR_CO_SHIP )
            ELSE
               CALL COMPUTE_EDGAR_CO( YEAR, EDGAR_CO )
            ENDIF
         ENDIF

         ! Reset YEAR_SAVE
         YEAR_SAVE = YEAR
      ENDIF

      !=================================================================
      ! NOx and SOX are annual emissions.  However, we will apply a
      ! seasonal scale factor to NOx and anthropogenic SO2, and a
      ! monthly scale factor to ship SO2.  So read these every month.
      !=================================================================
      IF ( MONTH /= MONTH_SAVE ) THEN

         ! Now always read NOX to get TODN (amv, phs, 3/10/08)
         IF ( LEDGARSHIP ) THEN
            CALL COMPUTE_EDGAR_NOX( YEAR,      MONTH, 
     &                              EDGAR_NOx, EDGAR_NOx_SHIP    )
         ELSE
            CALL COMPUTE_EDGAR_NOX( YEAR,      MONTH, 
     &                              EDGAR_NOx                    )
         ENDIF
         
         ! SO2
         IF ( LEDGARSOx ) THEN
            IF ( LEDGARSHIP ) THEN
               CALL COMPUTE_EDGAR_SO2( YEAR,      MONTH, 
     &                                 EDGAR_SO2, EDGAR_SO2_SHIP )
            ELSE
               CALL COMPUTE_EDGAR_SO2( YEAR,      MONTH, 
     &                                 EDGAR_SO2                 )
            ENDIF
         ENDIF
         
         ! Reset MONTH_SAVE
         MONTH_SAVE = MONTH
      ENDIF

      !=================================================================
      ! Print EDGAR emission totals
      !=================================================================
      CALL EDGAR_TOTAL_Tg( YEAR, MONTH )

      !### For Debug
      !CALL OUTPUT_TOTAL_2D( 'EDGAR NOx', EDGAR_NOx, 'kg/season' )
      !CALL OUTPUT_TOTAL_2D( 'EDGAR CO',  EDGAR_CO , 'kg/yr'     )
      !CALL OUTPUT_TOTAL_2D( 'EDGAR SO2', EDGAR_SO2, 'kg/season' )

#if defined( DEVEL )

      IF ( LFUTURE ) THEN
         DO J=1,JJPAR
         DO I=1,IIPAR
               SO2_FUTURE_SCALE(I,J) = GET_FUTURE_SCALE_SO2ff(I,J)
               CO_FUTURE_SCALE(I,J)  = GET_FUTURE_SCALE_COff(I,J)
               NOX_FUTURE_SCALE(I,J) = GET_FUTURE_SCALE_NOXff(I,J)
            ENDDO
         ENDDO
      ENDIF

      ! CO is yearly
      ! NOx and anthro SOx are seasonal
      ! Ship SOx are monthly
      ! Result is emission tendency in kg/gridbox/s
      ! Replaced NOx emissions with NO emissions (mpayer, 3/14/13)
      IF (LEDGARSHIP) THEN
         IF ( LFUTURE ) THEN
            State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTSO2) +
     &    (((EDGAR_SO2/SEC_IN_SEASON)+(EDGAR_SO2_SHIP/SEC_IN_MONTH)) *
     &      SO2_FUTURE_SCALE )

            State_Chm%TRAC_TEND(:,:,1,IDTCO) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTCO) +
     &    ((EDGAR_CO+EDGAR_CO_SHIP)   * (1.d0/SEC_IN_2000) *
     &      CO_FUTURE_SCALE )

            State_Chm%TRAC_TEND(:,:,1,IDTNO)  = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTNO)  +
     &    ((EDGAR_NOX+EDGAR_NOX_SHIP) * (1.d0/SEC_IN_SEASON) *
     &      NOX_FUTURE_SCALE )
         ELSE
            State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTSO2) +
     &    (((EDGAR_SO2/SEC_IN_SEASON)+(EDGAR_SO2_SHIP/SEC_IN_MONTH)) )

            State_Chm%TRAC_TEND(:,:,1,IDTCO) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTCO) +
     &    ((EDGAR_CO+EDGAR_CO_SHIP)   * (1.d0/SEC_IN_2000) )

            State_Chm%TRAC_TEND(:,:,1,IDTNO)  = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTNO)  +
     &    ((EDGAR_NOX+EDGAR_NOX_SHIP) * (1.d0/SEC_IN_SEASON) )
         ENDIF
      ELSE
         IF ( LFUTURE ) THEN
            State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTSO2) +
     &    ((EDGAR_SO2/SEC_IN_SEASON) * SO2_FUTURE_SCALE )

            State_Chm%TRAC_TEND(:,:,1,IDTCO) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTCO) +
     &    (EDGAR_CO*(1.d0/SEC_IN_2000) * CO_FUTURE_SCALE )

            State_Chm%TRAC_TEND(:,:,1,IDTNO)  = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTNO)  +
     &    (EDGAR_NOX*(1.d0/SEC_IN_SEASON) * NOX_FUTURE_SCALE )
         ELSE
            State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTSO2) +
     &    (EDGAR_SO2/SEC_IN_SEASON)

            State_Chm%TRAC_TEND(:,:,1,IDTCO) = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTCO) +
     &    (EDGAR_CO * (1.d0/SEC_IN_2000) )

            State_Chm%TRAC_TEND(:,:,1,IDTNO)  = 
     &      State_Chm%TRAC_TEND(:,:,1,IDTNO)  +
     &    (EDGAR_NOX * (1.d0/SEC_IN_SEASON) )
         ENDIF
      ENDIF
#endif

      ! Return to calling program
      END SUBROUTINE EMISS_EDGAR

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_EDGAR_NOx( YEAR,  MONTH, 
     &                              E_NOx, E_NOx_SHIP )
!
!******************************************************************************
!  Subroutine COMPUTE_EDGAR_NOx computes the total EDGAR NOx emissions 
!  (summing over several individual sectors) and also the time-of-day 
!  scale factors for NOx.  (avd, bmy, 7/14/06)
!
!  EDGAR NOx is read as [kg NO2/yr], then converted to [molec/cm2/s] when 
!  you call the access function GET_EDGAR_NOx( I, J, MOLEC_CM2_S=.TRUE. )
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR  (INTEGER) : Current year
!  (2 ) MONTH (INTEGER) : Current month
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) E_NOX      (REAL*4 ) : EDGAR NOx emissions      [kg NO2/season]
!  (4 ) E_NOX_SHIP (REAL*4 ) : EDGAR NOx ship emissions [kg NO2/season]
!
!  NOTES:
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              NOx (E. Leibensperger)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a impelementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR
      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER,        INTENT(IN)            :: YEAR, MONTH
      REAL*8,         INTENT(OUT)           :: E_NOX(IIPAR,JJPAR)
      REAL*8,         INTENT(OUT), OPTIONAL :: E_NOX_SHIP(IIPAR,JJPAR)

      ! Local variables
      INTEGER                    :: I, J, H, YYYYMMDD
      REAL*8, TARGET             :: E_NOX_1x1(I1x1,J1x1-1)
      REAL*8, TARGET             :: E_NOX_SHIP_1x1(I1x1,J1x1-1)
      REAL*8, TARGET             :: E_NOX_HRLY_1x1(I1x1,J1x1-1,N_HOURS)
      REAL*4                     :: SC_NOx(IIPAR,JJPAR)
      REAL*8                     :: TEMP_HRLY(IIPAR,JJPAR)
      REAL*8                     :: TEMP_TOT(IIPAR,JJPAR)
      CHARACTER(LEN=255)         :: NAME
      CHARACTER(LEN=255)         :: LLFILENAME
      REAL*8, POINTER            :: INGRID(:,:) => NULL()
      REAL*8                     :: OUTGRID(IIPAR,JJPAR)

      !=================================================================
      ! COMPUTE_EDGAR_NOx begins here!
      !=================================================================

      ! Initialize
      E_NOx_1x1      = 0d0
      E_NOX_HRLY_1x1 = 0d0
      SC_NOx         = 1d0

      !----------------------------------
      ! Read NOx data from disk
      !----------------------------------

      ! Read EDGAR anthro and ship NOx, and hourly NOx [kg NO2/yr]
      CALL READ_EDGAR_NOx( E_NOx_1x1, E_NOX_SHIP_1x1, E_NOx_HRLY_1x1 )

      !----------------------------------
      ! Compute NOx hourly scale factors
      ! (these average out to 1.0)
      !----------------------------------

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'
      
      ! Regrid anthro emissions to current resolution
      INGRID => E_NOx_1x1(:,:)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,     J1x1-1, 
     &                    INGRID,     TEMP_TOT, IS_MASS=1,
     &                    netCDF=.TRUE.                    )
      NULLIFY( INGRID )

      ! Loop over hours
      DO H = 1, N_HOURS

         ! Regrid hourly emissions to current resolution
         INGRID => E_NOx_HRLY_1x1(:,:,H)
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1-1, 
     &                       INGRID,     TEMP_HRLY, IS_MASS=1,
     &                       netCDF=.TRUE.                     )
         NULLIFY( INGRID )

         ! Scale factors are just the normalized hourly emissions
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( TEMP_TOT(I,J) > 0d0 ) THEN
               EDGAR_TODN(I,J,H) = TEMP_HRLY(I,J) / TEMP_TOT(I,J)
            ELSE
               EDGAR_TODN(I,J,H) = 1d0
            ENDIF
         ENDDO
         ENDDO
      ENDDO

      !----------------------------------
      ! Scale NOx from 2000 -> this year
      !----------------------------------
      CALL GET_ANNUAL_SCALAR( 71, 2000, YEAR, SC_NOx )

      !----------------------------------
      ! Do the regridding
      !----------------------------------

      ! Regrid NOx emissions to current model resolution [kg/yr]
      INGRID => E_NOx_1x1(:,:)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=1,
     &                    netCDF=.TRUE.                   )
      NULLIFY( INGRID )

      ! Apply scaling factor
![eml Skip scaling if using historical emissions
      IF ( .not. LHIST ) THEN
         E_NOx = OUTGRID * SC_NOx
      ELSE
         E_NOx = OUTGRID
      ENDIF
!eml]

      !----------------------------------
      ! Scale ship NOX and regrid
      !----------------------------------

      IF ( LEDGARSHIP ) THEN

         ! Regrid NOX emissions to current model resolution [kg/yr]
         INGRID => E_NOX_SHIP_1x1(:,:)
         CALL DO_REGRID_A2A( LLFILENAME,  I1x1,    J1x1-1, 
     &                       INGRID,      OUTGRID, IS_MASS=1,
     &                       netCDF=.TRUE. )

         E_NOX_SHIP = OUTGRID * SC_NOX(:,:) 
         NULLIFY( INGRID )

      ENDIF
    
      ! Return to calling program
      END SUBROUTINE COMPUTE_EDGAR_NOx

!------------------------------------------------------------------------------

      SUBROUTINE READ_EDGAR_NOx( T_1x1, T_SHIP_1x1, T_HRLY_1x1 )
!
!******************************************************************************
!  Subroutine READ_EDGAR_NOx reads EDGAR NOx emissions for the various sectors
!  and returns total and hourly emissions.  The EDGAR emissions are on the 
!  GENERIC 1x1 GRID and are regridded to the GEOS 1x1 grid. 
!  (avd, bmy, 7/14/06, 8/9/06)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_1x1      (REAL*8 ) : Anthro NOx on GEOS 1x1 grid [kg NO2/season]
!  (2 ) E_SHIP_1x1 (REAL*8 ) : Ship   NOx on GEOS 1x1 grid [kg NO2/season]
!  (2 ) E_HRLY_1x1 (REAL*8 ) : Hourly NOx on GEOS 1x1 grid [kg NO2/season]
!
!  NOTES:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              NOx (E. Leibensperger)
!  13 Mar 2012 - M. Cooper   - No longer regridding to the GEOS 1x1 grid.
!                              Regridding is done later
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
![eml
      USE TIME_MOD,       ONLY : GET_HISTYR
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
!eml]

      USE CMN_SIZE_MOD         ! Size parameters

      ! Arguments
      REAL*8,  INTENT(OUT)    :: T_1x1(I1x1,J1x1-1)
      REAL*8,  INTENT(OUT)    :: T_SHIP_1x1(I1x1,J1x1-1)
      REAL*8,  INTENT(OUT)    :: T_HRLY_1x1(I1x1,J1x1-1,N_HOURS)
                              
      ! Local variables       
      INTEGER                 :: H
      REAL*4                  :: E_IN(I1x1,J1x1-1)
      REAL*8                  :: SC(N_HOURS)    
![eml
      REAL*8                  :: TAU0
      CHARACTER(LEN=  4)      :: CHISTYR
      CHARACTER(LEN=255)      :: FILENAME
!eml]

      !=================================================================
      ! READ_EDGAR_NOx begins here!
      !=================================================================

      ! Initialize
      SC         = 0e0
      E_IN       = 0e0
      T_1x1      = 0d0
      T_SHIP_1x1 = 0d0
      T_HRLY_1x1 = 0d0

![eml The older NOx emissions are not sectorized. Read in a lumped diurnal
!     scale factor and those lumped emissions if LHIST.
      IF( .not. LHIST ) THEN
         !-----------------------------------------------------------------
         ! F10 - Industry (fuel combustion)
         !-----------------------------------------------------------------

          ! Hourly scale factors
         SC(:) = (/ 0.75d0, 0.75d0, 0.78d0, 0.82d0, 0.88d0, 0.95d0,
     &              1.02d0, 1.09d0, 1.16d0, 1.22d0, 1.28d0, 1.30d0,
     &              1.22d0, 1.24d0, 1.25d0, 1.16d0, 1.08d0, 1.01d0,
     &              0.95d0, 0.90d0, 0.85d0, 0.81d0, 0.78d0, 0.75d0 /)

         ! Read data
         CALL READ_EDGAR_DATA( 'f1000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! F20 - Power Generation (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = (/ 0.79d0, 0.72d0, 0.72d0, 0.71d0, 0.74d0, 0.80d0,
     &              0.92d0, 1.08d0, 1.19d0, 1.22d0, 1.21d0, 1.21d0,
     &              1.17d0, 1.15d0, 1.14d0, 1.13d0, 1.10d0, 1.07d0,
     &              1.04d0, 1.02d0, 1.02d0, 1.01d0, 0.96d0, 0.88d0 /)

         ! Read data
         CALL READ_EDGAR_DATA( 'f2000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! F30 - Conversion (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1e0 

         ! Read data
         CALL READ_EDGAR_DATA( 'f3000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )
       
         !-----------------------------------------------------------------
         ! F40 - Residential (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = (/ 0.38d0, 0.36d0, 0.36d0, 0.36d0, 0.37d0, 0.50d0,
     &              1.19d0, 1.53d0, 1.57d0, 1.56d0, 1.35d0, 1.16d0,
     &              1.07d0, 1.06d0, 1.00d0, 0.98d0, 0.99d0, 1.12d0,
     &              1.41d0, 1.52d0, 1.39d0, 1.35d0, 1.00d0, 0.42d0 /)

         ! Read data
         CALL READ_EDGAR_DATA( 'f4000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! F51 - Road Transport (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = (/ 0.19d0, 0.09d0, 0.06d0, 0.05d0, 0.09d0, 0.22d0,
     &              0.86d0, 1.84d0, 1.86d0, 1.41d0, 1.24d0, 1.20d0,
     &              1.32d0, 1.44d0, 1.45d0, 1.59d0, 2.03d0, 2.08d0,
     &              1.51d0, 1.06d0, 0.74d0, 0.62d0, 0.61d0, 0.44d0 /)

         ! Read data
         CALL READ_EDGAR_DATA( 'f5100nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! F54 - Non-Road Land Transport (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = (/ 0.19d0, 0.09d0, 0.06d0, 0.05d0, 0.09d0, 0.22d0,
     &              0.86d0, 1.84d0, 1.86d0, 1.41d0, 1.24d0, 1.20d0,
     &              1.32d0, 1.44d0, 1.45d0, 1.59d0, 2.03d0, 2.08d0,
     &              1.51d0, 1.06d0, 0.74d0, 0.62d0, 0.61d0, 0.44d0 /)

         ! Read data
         CALL READ_EDGAR_DATA( 'f5400nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA(  E_IN, T_1x1, SC, T_HRLY_1x1 )

!------------------------------------------------------------------------------
! NOTE: We don't use EDGAR aircraft data, so comment out.  (avd, bmy, 7/14/06)
!      !------------------------------------------------
!      ! F57 - Aircraft (fuel combustion)
!      !------------------------------------------------
!      SC(:) = (/ 0.19d0, 0.09d0, 0.06d0, 0.05d0, 0.09d0, 0.22d0, 
!     &           0.86d0, 1.84d0, 1.86d0, 1.41d0, 1.24d0, 1.20d0, 
!     &           1.32d0, 1.44d0, 1.45d0, 1.59d0, 2.03d0, 2.08d0,
!     &           1.51d0, 1.06d0, 0.74d0, 0.62d0, 0.61d0, 0.44d0 /)
!
!      ! Read data
!      CALL READ_EDGAR_FILE( 'f5700nox', 1, E_IN )!
!
!      ! Add into cumulative arrays
!      CALL ADD_EDGAR_DATA(  E_IN, T_1x1, SC, T_HRLY_1x1 )
!------------------------------------------------------------------------------

         !-----------------------------------------------------------------
         ! F80 - Oil Production (fuel combustion)
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'f8000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! I10 - Iron and Steel Production
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'i1000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! I30 - Chemical Production
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'i3000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! I40 - Cement Production
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'i4000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! I50 - Pulp and Paper Production
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'i5000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         !-----------------------------------------------------------------
         ! W40 - Waste Incineration
         !-----------------------------------------------------------------

         ! Hourly scale factors
         SC(:) = 1d0

         ! Read data
         CALL READ_EDGAR_DATA( 'w4000nox', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

      ELSE

         !-----------------------------------------------------------------
         ! Historical emissions of NOX
         !-----------------------------------------------------------------

         ! Make sure to do Ship nox

         ! Read in lumped diurnal scale factors
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/NOx/'         //
     &              'lumped_fnox_diurnal_scale_factors.dat'

         OPEN(412, FILE=FILENAME, FORM='formatted')
         READ(412,*) SC
         CLOSE(412)

         ! Get year
         WRITE( CHISTYR, '(I4)' ) GET_HISTYR()

         ! Get TAU0 for this year
         TAU0 = GET_TAU0( 1, 1, GET_HISTYR() )

         ! Read fossil fuel combustion NOx [kg/yr]
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/NOx/'         //
     &              'EDGAR.FCNOX.' // CHISTYR // 'generic.1x1'

         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_EDGAR_NOx: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 1,
     &                    TAU0,           I1x1,       J1x1-1,
     &                    1,              E_IN,       QUIET=.TRUE. )

         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

         ! Read industrial NOx [kg/yr]
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/NOx/'         //
     &              'EDGAR.INNOX.' // CHISTYR // 'generic.1x1'

         WRITE( 6, 100 ) TRIM( FILENAME )
 110     FORMAT( '     - READ_EDGAR_NOx: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 1,
     &                    TAU0,            I1x1,      J1x1-1,
     &                    1,               E_IN,      QUIET=.TRUE. )

         SC(:) = 1

         CALL ADD_EDGAR_DATA( E_IN, T_1x1, SC, T_HRLY_1x1 )

      ENDIF
!eml]

      !-----------------------------------------------------------------
      ! F58 - Shipping (fuel combustion) - done regardless of historical 
      !       emissions
      !-----------------------------------------------------------------
      IF ( LEDGARSHIP ) THEN
         
         ! Hourly scale factors
         SC(:) = 1d0 
         
         ! Read data
         CALL READ_EDGAR_DATA( 'f5800nox(IEA)', 1, E_IN )

         ! Add into cumulative arrays
         CALL ADD_EDGAR_DATA( E_IN, T_SHIP_1x1 )
      ENDIF

      !-----------------------------------------------------------------
      ! Force a seasonal variation onto the anthropogenic NOx emissions
      ! by applying seasonal scale factors.  The scale factors are the
      ! ratio of (seasonal GEIA NOx / annual GEIA NOx).
      !
      ! The emissions on which these scale factors are based are 
      ! defined on the GENERIC 1 x 1 GRID, so apply scale factors 
      ! BEFORE regridding!
      !-----------------------------------------------------------------

      ! Convert [kg NO2/yr] to [kg NO2/season]
      CALL SEASCL_EDGAR_NOx( T_1x1, T_SHIP_1x1, T_HRLY_1x1 )

      ! Return to calling program
      END SUBROUTINE READ_EDGAR_NOx

!------------------------------------------------------------------------------

      SUBROUTINE SEASCL_EDGAR_NOx( E_NOx_1x1,     E_NOx_SHIP_1x1, 
     &                             E_NOx_HRLY_1x1 )
!
!******************************************************************************
!  Subroutine SEASCL_EDGAR_NOx applies seasonal scale factors (computed
!  as the ratio of seasonal/total GEIA NOx emissions) to the annual EDGAR
!  anthropogenic NOx emissions.  This is required to impose a seasonality 
!  onto the EDGAR NOx emissions, which are reported as per year. 
!  (avd, bmy, 7/14/06)
!
!  NOTE: NOx scale factors are on the GENERIC 1x1 GRID.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) E_NOx_1x1      (REAL*8 ) : Anthro NOx 1x1        array [kg NO2/yr]
!  (2 ) E_NOx_1x1      (REAL*8 ) : Anthro NOx 1x1 hourly array [kg NO2/yr]
!  (3 ) E_NOx_SHIP_1x1 (REAL*8 ) : Ship   NOx 1x1        array [kg NO2/yr]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_NOx_1x1      (REAL*8 ) : Anthro NOx 1x1        array [kg NO2/season]
!  (2 ) E_NOx_1x1      (REAL*8 ) : Anthro NOx 1x1 hourly array [kg NO2/season]
!  (3 ) E_NOx_SHIP_1x1 (REAL*8 ) : Ship   NOx 1x1        array [kg NO2/yr]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1

      USE CMN_SIZE_MOD        ! Size parameters
 
      ! Arguments
      REAL*8,  INTENT(INOUT) :: E_NOx_1x1(I1x1,J1x1-1)
      REAL*8,  INTENT(INOUT) :: E_NOx_SHIP_1x1(I1x1,J1x1-1)
      REAL*8,  INTENT(INOUT) :: E_NOx_HRLY_1x1(I1x1,J1x1-1,N_HOURS)
      
      ! Local variables
      INTEGER                :: H
      REAL*4                 :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                 :: THIS_TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! SEASCL_EDGAR_NOx begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'EDGAR_200607/NOx/anth_NOx_scale.' // SEASON_NAME //
     &           '.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - SEASCL_EDGAR_NOx: Reading ', a )
 
      ! Read scale factor data [unitless]
      CALL READ_BPCH2( FILENAME,   'EDGAR-2D', 71, 
     &                 SEASON_TAU0, I1x1,      J1x1-1,     
     &                 1,           ARRAY,     QUIET=.TRUE. ) 

      
      ! Apply seasonal scale factors to anthro and ship NOx
      E_NOx_1x1(:,:)      = E_NOx_1x1(:,:)      * ARRAY(:,:,1)
      E_NOx_SHIP_1x1(:,:) = E_NOx_SHIP_1x1(:,:) * ARRAY(:,:,1)

      ! Apply seasonal scale factors to hourly anthro NOx
      DO H = 1, N_HOURS
         E_NOx_HRLY_1x1(:,:,H) = E_NOx_HRLY_1x1(:,:,H) * ARRAY(:,:,1)
      ENDDO

      ! Return to calling program
      END SUBROUTINE SEASCL_EDGAR_NOx

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_EDGAR_CO( YEAR, E_CO, E_CO_SHIP )
!
!******************************************************************************
!  Subroutine COMPUTE_EDGAR_CO computes the total EDGAR CO emissions, summing 
!  over several individual sectors. (avd, bmy, 7/14/06)
!
!  EDGAR CO is read as [kg CO/yr], then converted to [molec/cm2/s] when you 
!  call the access function GET_EDGAR_CO( I, J, MOLEC_CM2_S=.TRUE. )
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR (INTEGER) : Current year
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) E_CO      (REAL*4 ) : EDGAR CO emissions          [kg CO/year]
!  (3 ) E_CO_SHIP (REAL*4 ) : EDGAR CO shipping emissions [kg CO/year]
!
!  NOTES:
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,    ONLY : DATA_DIR_1x1
      USE TIME_MOD,         ONLY : EXPAND_DATE
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      USE REGRID_A2A_MOD,   ONLY : DO_REGRID_A2A
      USE CMN_SIZE_MOD           ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)            :: YEAR
      REAL*8,  INTENT(OUT)           :: E_CO(IIPAR,JJPAR)
      REAL*8,  INTENT(OUT), OPTIONAL :: E_CO_SHIP(IIPAR,JJPAR)

      ! Local variables
      INTEGER                        :: I, J, H, YYYYMMDD
      REAL*8, TARGET                 :: E_CO_1x1(I1x1,J1x1-1)
      REAL*8, TARGET                 :: E_CO_SHIP_1x1(I1x1,J1x1-1)
      REAL*4                         :: SC_CO(IIPAR,JJPAR)
      CHARACTER(LEN=255 )            :: NAME
      CHARACTER(LEN=255)             :: LLFILENAME
      REAL*8                         :: OUTGRID(IIPAR,JJPAR)
      REAL*8,  POINTER               :: INGRID(:,:) => NULL()

      !=================================================================
      ! COMPUTE_EDGAR_CO begins here!
      !=================================================================

      ! Initialize
      SC_CO(:,:) = 1e0

      !----------------------------------
      ! Read CO data from disk
      !----------------------------------

      ! Read EDGAR total NOx and hourly NOx [kg/yr]
      CALL READ_EDGAR_CO( E_CO_1x1, E_CO_SHIP_1x1 )

      !----------------------------------
      ! Do the regridding
      !----------------------------------

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid CO emissions to current model resolution [kg/yr]
      INGRID => E_CO_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=1,
     &                    netCDF=.TRUE.                   )
      E_CO = OUTGRID

      ! Free pointer
      NULLIFY( INGRID )

      !----------------------------------
      ! Scale CO from 2000 -> this year
      !----------------------------------
      CALL GET_ANNUAL_SCALAR( 72, 2000, YEAR, SC_CO )

      ! Apply scale factors
      E_CO = E_CO * SC_CO(:,:)

      !----------------------------------
      ! Scale ship CO and regrid
      !----------------------------------

      IF ( LEDGARSHIP ) THEN

         ! Regrid CO emissions to current model resolution [kg/yr]
         INGRID => E_CO_SHIP_1x1
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                       INGRID,     OUTGRID, IS_MASS=1,
     &                       netCDF=.TRUE.                   )

         ! Apply scale factors
         E_CO_SHIP = OUTGRID * SC_CO
        
         ! Free pointer
         NULLIFY( INGRID )

      ENDIF

      ! Return to calling program
      END SUBROUTINE COMPUTE_EDGAR_CO

!------------------------------------------------------------------------------

      SUBROUTINE READ_EDGAR_CO( T_CO_1x1, T_CO_SHIP_1x1 )
!
!******************************************************************************
!  Subroutine READ_EDGAR_CO reads EDGAR CO emissions for the various sectors
!  and returns total and hourly emissions.  The EDGAR emissions are on the 
!  GENERIC 1x1 GRID and are regridded to the GEOS 1x1 grid. 
!  (avd, bmy, 7/14/06, 8/9/06)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_CO_1x1 (REAL*4) : Total EDGAR CO emissions on GEOS 1x1 grid [kg/yr]
!
!  NOTES:      
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!******************************************************************************
!
      ! Reference to F90 modules
      USE CMN_SIZE_MOD      ! Size parameters

      ! Arguments
      REAL*8, INTENT(OUT)  :: T_CO_1x1(I1x1,J1x1-1)
      REAL*8, INTENT(OUT)  :: T_CO_SHIP_1x1(I1x1,J1x1-1)

      ! Local variables
      REAL*4               :: E_IN(I1x1,J1x1-1)

      !=================================================================
      ! READ_EDGAR_CO begins here!
      !=================================================================

      ! Initialize
      E_IN          = 0e0
      T_CO_1x1      = 0d0
      T_CO_SHIP_1x1 = 0d0

      !------------------------------------------------
      ! Compute total CO for all sectors
      !------------------------------------------------

      ! F10 - Industrial (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f1000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )
    
      ! F20 - Power Generation (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f2000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F30 - Conversion (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f3000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F40 - Residential + Commercial + Other  (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f4000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F51 - Road Transport (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f5100co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F54 - Land (Non-Road) Transport (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f5400co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F57 - Air Transport (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f5700co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! F58 - Shipping (fossil fuel combustion)
      IF ( LEDGARSHIP ) THEN
         CALL READ_EDGAR_DATA( 'f5800co', 4, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_CO_SHIP_1x1 )
      ENDIF

      ! F80 - Oil Production (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f8000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! I10 - Iron and Steel Production
      CALL READ_EDGAR_DATA( 'i1000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! I20 - Non-Ferrous Production
      CALL READ_EDGAR_DATA( 'i2000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! I50 - Pulp and Paper Production
      CALL READ_EDGAR_DATA( 'i5000co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! W40 - Waste Incineration
      CALL READ_EDGAR_DATA( 'w4095co', 4, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_CO_1x1 )

      ! Return to calling program
      END SUBROUTINE READ_EDGAR_CO

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_EDGAR_SO2( YEAR, MONTH, E_SO2, E_SO2_SHIP )
!
!******************************************************************************
!  Subroutine COMPUTE_EDGAR_SO2 computes the total EDGAR SO2 emissions 
!  (summing over several individual sectors) and also the time-of-day scale
!  factors for SO2.  (avd, bmy, 7/14/06)
!
!  EDGAR anthropogenic SO2 is read as [kg SO2/yr], then converted to [kg/s] 
!  when you call the access functions GET_EDGAR_ANTH_SO2( I, J, KG_S=.TRUE. )
!  and GET_EDGAR_ANTH_SO2( I, J, KG_S=.TRUE. )
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR       (INTEGER) : Current year
!  (2 ) MONTH      (INTEGER) : Current month
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) E_SO2      (REAL*4 ) : EDGAR anthropogenic SO2 emissions
!  (4 ) E_SO2_SHIP (REAL*4 ) : EDGAR ship exhaust  SO2 emissions
!
!  NOTES:
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              SO2 (E. Leibensperger)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  25 Mar 2013 - R. Yantosca - Now use logical fields from Input_Opt
!  12 Jun 2013 - R. Yantosca - Remove local declaration of LEDGARSHIP & LHIST,
!                              which obscures the module variable definitions
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A      
      USE BPCH2_MOD,          ONLY : READ_BPCH2
      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)              :: YEAR, MONTH
      REAL*8,  INTENT(OUT)             :: E_SO2(IIPAR,JJPAR)
      REAL*8,  INTENT(INOUT), OPTIONAL :: E_SO2_SHIP(IIPAR,JJPAR)

      ! Local variables
      INTEGER                          :: I, J, H, YYYYMMDD
      REAL*8,  TARGET                  :: E_SO2_1x1(I1x1,J1x1-1)
      REAL*8,  TARGET                  :: E_SO2_SHIP_1x1(I1x1,J1x1-1)
      REAL*4                           :: SC_SO2(IIPAR,JJPAR)
      CHARACTER(LEN=255)               :: NAME
      CHARACTER(LEN=255)               :: LLFILENAME
      REAL*8                           :: OUTGRID(IIPAR,JJPAR)
      REAL*8,  POINTER                 :: INGRID(:,:) => NULL()

      !=================================================================
      ! COMPUTE_EDGAR_SO2 begins here!
      !=================================================================

      ! Initialize
      SC_SO2(:,:) = 1d0

      !----------------------------------
      ! Read SO2 data from disk and regrid
      !----------------------------------

      ! Read EDGAR anthro SO2 and ship SO2
      CALL READ_EDGAR_SO2( E_SO2_1x1, E_SO2_SHIP_1x1 )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid SO2 emissions to current model resolution [kg/yr]
      INGRID => E_SO2_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=1,
     &                    netCDF=.TRUE.                   )

      ! Free pointer
      NULLIFY( INGRID )

      !----------------------------------
      ! Scale SO2 from 2000 -> this year
      !----------------------------------
      CALL GET_ANNUAL_SCALAR( 73, 2000, YEAR, SC_SO2 ) 

      !----------------------------------
      ! Scale anthro SO2
      !----------------------------------

      ! Apply scale factors
![eml Skip scaling if using historical emissions
      IF ( .not. LHIST ) THEN
         E_SO2(:,:) = OUTGRID * SC_SO2
      ELSE
         E_SO2(:,:) = OUTGRID
      ENDIF
!eml]

      !----------------------------------
      ! Scale ship SO2 and regrid
      !----------------------------------

      IF ( LEDGARSHIP ) THEN

         ! Regrid SO2 emissions to current model resolution [kg/yr]
         INGRID => E_SO2_SHIP_1x1
         CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                       INGRID,     OUTGRID, IS_MASS=1,
     &                       netCDF=.TRUE.                   )
                
         ! Apply scale factors
         E_SO2_SHIP = OUTGRID * SC_SO2

         ! Free pointer
         NULLIFY( INGRID )

         ! Apply seasonal scaling
         CALL SEASCL_EDGAR_SHIP_SO2( E_SO2_SHIP )

      ENDIF

      ! Return to calling program
      END SUBROUTINE COMPUTE_EDGAR_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_EDGAR_SO2( T_SO2_1x1, T_SO2_SHIP_1x1 )
!
!******************************************************************************
!  Subroutine READ_EDGAR_SO2 reads EDGAR SO2 emissions for the various sectors
!  and returns both anthropogenic SO2 emissions and ship exhaust SO2 emissions.
!  The EDGAR emissions are on the GENERIC 1x1 GRID and are regridded to the 
!  GEOS 1x1 GRID. (avd, bmy, 7/14/06, 8/9/06)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_SO2_1x1      (REAL*8) : EDGAR anth SO2 on GEOS 1x1 GRID [kg/season]
!  (2 ) E_SO2_SHIP_1x1 (REAL*8) : EDGAR ship SO2 on GEOS 1x1 GRID [kg/season]
!
!  NOTES:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              SO2 (E. Leibensperger)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
![eml
      USE TIME_MOD,       ONLY : GET_HISTYR
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
!eml]

      USE CMN_SIZE_MOD         ! Size parameters

      ! Arguments
      REAL*8,  INTENT(OUT)    :: T_SO2_1x1(I1x1,J1x1-1)
      REAL*8,  INTENT(OUT)    :: T_SO2_SHIP_1x1(I1x1,J1x1-1)

      ! Local variables
      REAL*4                  :: E_IN(I1x1,J1x1-1)
![eml
      REAL*8                  :: TAU0
      CHARACTER(LEN=  4)      :: CHISTYR
      CHARACTER(LEN=255)      :: FILENAME
!eml]

      !=================================================================
      ! READ_EDGAR_SO2 begins here!
      !=================================================================

      ! Initialize
      E_IN           = 0e0
      T_SO2_1x1      = 0d0
      T_SO2_SHIP_1x1 = 0d0

![eml
      IF ( .not. LHIST ) THEN
      
         !-----------------------------------------------------------------
         ! Read anthropogenic SO2 and ship SO2 emissions
         ! (on GENERIC 1x1 GRID)
         !-----------------------------------------------------------------

         ! F10 - industrial (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f1000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F20 - Power Generation (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f2000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F30 - Conversion (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f3000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F40 - Residential + Commercial + Other  (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f4000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F51 - Road Transport (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f5100so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F54 - Land (Non-Road) Transport (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f5400so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! F80 - Oil Production (fossil fuel combustion)
         CALL READ_EDGAR_DATA( 'f8000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! I10 - Iron and Steel Production
         CALL READ_EDGAR_DATA( 'i1000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! I20 - Non-Ferrous Production
         CALL READ_EDGAR_DATA( 'i2000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! I30 - Chemical Production
         CALL READ_EDGAR_DATA( 'i3000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! I40 - Cement Production
         CALL READ_EDGAR_DATA( 'i4000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )
 
         ! I50 - Pulp and Paper Production
         CALL READ_EDGAR_DATA( 'i5000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )
 
         ! W40 - Waste Incineration
         CALL READ_EDGAR_DATA( 'w4000so2', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

      ELSE

         !-----------------------------------------------------------------
         ! Historical emissions of SO2
         !-----------------------------------------------------------------

         WRITE( CHISTYR, '(I4)' ) GET_HISTYR()

         ! Get TAU0 for this year
         TAU0 = GET_TAU0( 1, 1, GET_HISTYR() )

         ! Read fossil fuel combustion SO2
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/SO2/'         //
     &              'EDGAR.FCSO2.' // CHISTYR // 'generic.1x1'

         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_EDGAR_SO2: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 26,
     &                    TAU0,           I1x1,       J1x1-1,
     &                    1,              E_IN,       QUIET=.TRUE. )

         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

         ! Read industrial SO2
         FILENAME = TRIM( DATA_DIR_1x1 )                       // 
     &              'historical_emissions_201203/SO2/'         //
     &              'EDGAR.INSO2.' // CHISTYR // 'generic.1x1'

         WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( '     - READ_EDGAR_SO2: Reading ', a )
         CALL READ_BPCH2( TRIM(FILENAME), 'EDGAR-2D', 26,
     &                    TAU0,           I1x1,       J1x1-1,
     &                    1,              E_IN,       QUIET=.TRUE. )
         
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

      ENDIF

      ! Do aircraft and ship emissions even if using historical
      ! F57 - Air Transport (fossil fuel combustion)
      CALL READ_EDGAR_DATA( 'f5700so2', 26, E_IN )
      CALL ADD_EDGAR_DATA( E_IN, T_SO2_1x1 )

      ! F58 - Shipping (fossil fuel combustion)
      ! NOTE: Add into separate array for compatibility w/ sulfate_mod.f
      IF ( LEDGARSHIP ) THEN
         CALL READ_EDGAR_DATA( 'f5800so2(IEA)', 26, E_IN )
         CALL ADD_EDGAR_DATA( E_IN, T_SO2_SHIP_1x1 )
      ENDIF
!eml]

      !-----------------------------------------------------------------
      ! Force a seasonal variation onto the anthropogenic SO2 emissions
      ! by applying seasonal scale factors.  The scale factors are the
      ! ratio of (seasonal GEIA SO2 / annual GEIA SO2).
      !
      ! The emissions on which these scale factors are based are 
      ! defined on the GENERIC 1 x 1 GRID, so apply scale factors 
      ! BEFORE regridding!
      !-----------------------------------------------------------------

      ! Convert [kg SO2/yr] to [kg SO2/season]
      CALL SEASCL_EDGAR_ANTH_SO2( T_SO2_1x1 )

      !-----------------------------------------------------------------
      ! Force a monthly variation onto the anthropogenic SO2 emissions
      ! by applying monthly scale factors.  The scale factors are the
      ! ratio of (monthly ship SO2 / total ship SO2) as take from the
      ! inventory of Corbett et al. 
      !
      ! The emissions on which these scale factors are based are 
      ! defined on the GEOS 1 x 1 GRID, so apply scale factors 
      ! AFTER regridding!
      !-----------------------------------------------------------------

      ! Convert [kg SO2/yr] to [kg SO2/month]
      ! ->Temporarily moved to compute_edgar_so2, where scaling is done at
      !   model resolution
      !CALL SEASCL_EDGAR_SHIP_SO2( T_SO2_SHIP_1x1 )

      ! Return to calling program
      END SUBROUTINE READ_EDGAR_SO2

!------------------------------------------------------------------------------

      SUBROUTINE SEASCL_EDGAR_ANTH_SO2( E_SO2_1x1 )
!
!******************************************************************************
!  Subroutine SEASCL_EDGAR_ANTH_SO2 applies seasonal scale factors (computed
!  as the ratio of seasonal/total GEIA SO2 emissions) to the annual EDGAR
!  anthropogenic SO2 emissions.  This is required to impose a seasonality onto
!  the EDGAR ship SO2 emissions, which are reported as per year. 
!  (avd, bmy, 7/14/06)
!
!  NOTE: Ship SO2 scale factors are on the GENERIC 1x1 GRID.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) E_SO2_1x1 (REAL*8 ) : Anthro SO2 1x1 array [kg SO2/yr]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_SO2_1x1 (REAL*8 ) : Anthro SO2 1x1 array [kg SO2/season]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1

      USE CMN_SIZE_MOD        ! Size parameters
 
      ! Arguments
      REAL*8,  INTENT(INOUT) :: E_SO2_1x1(I1x1,J1x1-1)
      
      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(I1x1,J1x1-1,1)
      REAL*8                 :: THIS_TAU
      CHARACTER(LEN=3)       :: THIS_SEA
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! SEASCL_EDGAR_ANTH_SO2 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'EDGAR_200607/SOx/anth_SOx_scale.' // SEASON_NAME //
     &           '.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - SEASCL_EDGAR_ANTH_SO2: Reading ', a )
 
      ! Read scale factor data [unitless]
      CALL READ_BPCH2( FILENAME,   'EDGAR-2D', 73, 
     &                 SEASON_TAU0, I1x1,      J1x1-1,     
     &                 1,           ARRAY,     QUIET=.TRUE. ) 

      
      ! Apply seasonal scale factors to anthro SO2
      E_SO2_1x1(:,:) = E_SO2_1x1(:,:) * ARRAY(:,:,1)

      ! Return to calling program
      END SUBROUTINE SEASCL_EDGAR_ANTH_SO2

!------------------------------------------------------------------------------

      SUBROUTINE SEASCL_EDGAR_SHIP_SO2( E_SO2_SHIP )
!
!******************************************************************************
!  Subroutine SEASCL_EDGAR_SHIP_SO2 applies monthly scale factors (which are 
!  computed as the ratio of monthly/total ship SO2 emissions from Corbett et 
!  al) to the annual EDGAR ship SO2 emissions.  This is required to impose a
!  seasonality onto the EDGAR ship SO2 emissions, which are reported as per
!  year. (avd, bmy, 7/14/06)
!
!  NOTE: Ship SO2 scale factors are on the GEOS 1x1 GRID.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) E_SO2_SHIP (REAL*8 ) : Ship SO2 array [kg SO2/yr]
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) E_SO2_SHIP (REAL*8 ) : Ship SO2 array [kg SO2/month]
!
!  NOTES:
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Bug fix for regridding. Changed to PERAREA=0
!                              since data is unitless.
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters
 
      ! Arguments
      REAL*8,  INTENT(INOUT)  :: E_SO2_SHIP(IIPAR,JJPAR) 

      ! Local variables
      INTEGER                 :: I, J
      REAL*4                  :: ARRAY(I1x1,J1x1,1)
      REAL*8                  :: THIS_TAU
      CHARACTER(LEN=3)        :: THIS_MON
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: INGRID (I1x1,J1x1)
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)

      !=================================================================
      ! SEASCL_EDGAR_SHIP_SO2 begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 )               // 
     &           'EDGAR_200607/SOx/ship_SOx_scale.' // MONTH_NAME //
     &           '.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - SEASCL_EDGAR_SHIP_SO2: Reading ', a )
 
      ! Read scale factor data [unitless]
      CALL READ_BPCH2( FILENAME,  'EDGAR-2D', 73, 
     &                 MONTH_TAU0, I1x1,      J1x1,     
     &                 1,          ARRAY,     QUIET=.TRUE. ) 
      
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid to current model resolution
      INGRID = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1, 
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )

      ! Apply monthly scale factors to ship SO2
      E_SO2_SHIP = E_SO2_SHIP * OUTGRID

      ! Return to calling program
      END SUBROUTINE SEASCL_EDGAR_SHIP_SO2

!------------------------------------------------------------------------------

      SUBROUTINE READ_EDGAR_DATA( NAME, TRACER, E_1x1 )
!
!******************************************************************************
!  Subroutine READ_EDGAR_DATA reads EDGAR emissions data for a single sector
!  from disk, in binary punch file format. (avd, bmy, 7/14/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : string with EDGAR inventory filename
!  (2 ) TRACER   (INTEGER  ) : Tracer number
!  (3 ) E_1x1    (REAL*4   ) : Array to hold emissions
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD,        ONLY : DATA_DIR_1x1
      USE BPCH2_MOD,            ONLY : GET_TAU0, READ_BPCH2

      USE CMN_SIZE_MOD               ! Size parameters

      ! Arguments
      INTEGER,          INTENT(IN)  :: TRACER
      REAL*4,           INTENT(OUT) :: E_1x1(I1x1,J1x1-1)
      CHARACTER(LEN=*), INTENT(IN)  :: NAME

      ! Local variables
      REAL*8                        :: TAU0
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_EDGAR_DATA begins here!
      !=================================================================

      ! Filename
      IF ( TRACER == 1 ) THEN

         ! NOx
         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/NOx/EDGAR.' // 
     &              TRIM( NAME         ) // '.generic.1x1'
         
      ELSE IF ( TRACER == 4 ) THEN

         ! CO
         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/CO/EDGAR.'  // 
     &              TRIM( NAME         ) // '.generic.1x1'

      ELSE IF ( TRACER == 26 ) THEN

         ! SO2
         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/SOx/EDGAR.' // 
     &              TRIM( NAME         ) // '.generic.1x1'

      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_EDGAR_DATA: Reading ', a )

      ! Use TAU0 for year 2000
      TAU0 = GET_TAU0( 1, 1, 2000 )

      ! Read data [kg/yr]
      CALL READ_BPCH2( FILENAME, 'EDGAR-2D', TRACER, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         E_1x1,     QUIET=.TRUE. ) 

      ! Return to calling program
      END SUBROUTINE READ_EDGAR_DATA

!------------------------------------------------------------------------------

      SUBROUTINE ADD_EDGAR_DATA( E, E_1x1, SC, E_HRLY_1x1 )
!
!******************************************************************************
!  Subroutine ADD_EDGAR_DATA adds emissions for a given sector to cumulative
!  data arrays.  Computes total emissions arrays and (if requested) totals per
!  hour. (avd, bmy, 7/14/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) E          (REAL*4) : EDGAR emissions for a given sector       [kg/yr]
!  (3 ) SC         (REAL*8) : OPTIONAL - Hourly EDGAR scale factors [unitless]
! 
!  Arguments as Output:
!  ============================================================================
!  (2 ) E_1x1      (REAL*8) : Cumulative total of EDGAR emissions      [kg/yr]
!  (4 ) E_HRLY_1x1 (REAL*8) : OPTIONAL - Hourly cum total of emissions [kg/yr]
!
!  NOTES:
!******************************************************************************
!
      USE CMN_SIZE_MOD               ! Size parameters

      ! Arguments
      REAL*4, INTENT(IN)            :: E(I1x1,J1x1-1)
      REAL*8, INTENT(OUT)           :: E_1x1(I1x1,J1x1-1)
      REAL*8, INTENT(IN),  OPTIONAL :: SC(N_HOURS)
      REAL*8, INTENT(OUT), OPTIONAL :: E_HRLY_1x1(I1x1,J1x1-1,N_HOURS)

      ! Local variables
      LOGICAL                       :: IS_HRLY
      INTEGER                       :: I, J, H

      !=================================================================
      ! ADD_EDGAR_DATA begins here!
      !=================================================================
      
      ! Are we computing cumulative totals for each hour?
      IS_HRLY    = ( PRESENT( E_HRLY_1x1 ) .and. PRESENT( SC ) )

      ! Create total sum
      E_1x1(:,:) = E_1x1(:,:) + E(:,:)
      
      ! Create hourly sum
      IF ( IS_HRLY ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, H )
         DO H = 1, N_HOURS
         DO J = 1, J1x1-1
         DO I = 1, I1x1
            E_HRLY_1x1(I,J,H) = E_HRLY_1x1(I,J,H) + ( SC(H) * E(I,J) )
         ENDDO 
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_EDGAR_DATA

!------------------------------------------------------------------------------

!!              ######### OBSOLETE (phs, 3/10/08) #########  
!!
!      SUBROUTINE READ_EDGAR_SCALE( NAME, TRACER, YEAR, S_1x1 )
!!
!!******************************************************************************
!!  Subroutine READ_EDGAR_SCALE reads interannual scale factor data from disk.
!!  (avd, bmy, 7/14/06)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) FILENAME (CHARACTER) : String with EDGAR inventory filename
!!  (2 ) TRACER   (INTEGER  ) : Tracer number
!!  (3 ) YEAR     (INTEGER  ) : Current year
!!  (4 ) S_1x1    (REAL*4   ) : Array to hold scale factors
!!
!!  NOTES: ######### OBSOLETE (phs, 3/10/08) #########  
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DIRECTORY_MOD,        ONLY : DATA_DIR_1x1
!      USE BPCH2_MOD,            ONLY : GET_TAU0, READ_BPCH2
!
!      USE CMN_SIZE_MOD             ! Size parameters
!
!      ! Arguments
!      CHARACTER(LEN=*), INTENT(IN)  :: NAME
!      INTEGER,          INTENT(IN)  :: TRACER, YEAR
!      REAL*8,           INTENT(OUT) :: S_1x1(I1x1,J1x1)
!
!      ! Local variables
!      REAL*4                        :: T_1x1(I1x1,J1x1)
!      REAL*8                        :: TAU0
!      CHARACTER(LEN=255)            :: FILENAME
!
!      !=================================================================
!      ! READ_EDGAR_SCALE begins here!
!      !=================================================================
!
!      ! Filename
!      IF ( TRACER == 71 ) THEN
!
!         ! NOx
!         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/NOx/EDGAR.' // 
!     &              TRIM( NAME         ) // '.geos.1x1'
!         
!      ELSE IF ( TRACER == 72 ) THEN
!
!         ! CO
!         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/CO/EDGAR.'  // 
!     &              TRIM( NAME         ) // '.geos.1x1'
!
!      ELSE IF ( TRACER == 73 ) THEN
!
!         ! SO2
!         FILENAME = TRIM( DATA_DIR_1x1 ) // 'EDGAR_200607/SOx/EDGAR.' // 
!     &              TRIM( NAME         ) // '.geos.1x1'
!
!      ENDIF
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_EDGAR_SCALE: Reading ', a )
!
!      ! Use TAU0 for the current year
!      TAU0 = GET_TAU0( 1, 1, YEAR )
!
!      ! Read data
!      CALL READ_BPCH2( FILENAME, 'EDGAR-2D', TRACER, 
!     &                 TAU0,      I1x1,      J1x1,     
!     &                 1,         T_1x1,     QUIET=.TRUE. ) 
!
!      ! Convert to REAL*8 and return
!      S_1x1(:,:) = T_1x1(:,:)
!
!      ! Return to calling program
!      END SUBROUTINE READ_EDGAR_SCALE
!
!!------------------------------------------------------------------------------
!
      FUNCTION GET_EDGAR_NOx( I, J, KG_S, MOLEC_CM2_S, SHIP ) 
     &         RESULT( NOx )
!
!******************************************************************************
!  Function GET_EDGAR_NOx returns the EDGAR NOx emissions at grid box (I,J)
!  in units of [kg/season], [kg/s], or [molec/cm2/s]. (avd, bmy, 7/14/06)
!
!  Arguments: 
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude  index
!  (3 ) KG_S        (LOGICAL) : OPTIONAL - Return data in [kg/s]
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL - Return data in [molec/cm2/s]
!
!  NOTES:
!  (1 ) Added the SHIP switch (phs, 5/13/08)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff

      ! Arguments
      INTEGER, INTENT(IN)           :: I,      J
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S,   MOLEC_CM2_S, SHIP 

      ! Local variables
      LOGICAL                       :: DO_KGS, DO_MCS
      REAL*8                        :: NOx
      
      !=================================================================
      ! GET_EDGAR_NOx begins here!
      !=================================================================
    
      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Get NOx [kg NOx/season]
      IF ( PRESENT( SHIP ) ) THEN
         IF ( SHIP ) THEN
            NOx = EDGAR_NOx_SHIP(I,J)
         ELSE
            NOx = EDGAR_NOx(I,J)
         ENDIF
      ELSE
         NOx = EDGAR_NOx(I,J)
      ENDIF

!%%%%%%%%%%%%% KLUDGE TO EMITT SHIP NOX AS NOX %%%%%%%%%%%%%%
!      ! Get NOx [kg NOx/season]
!      IF ( PRESENT( SHIP ) ) THEN
!         NOx = 0d0
!      ELSE
!         NOx = EDGAR_NOx(I,J) + EDGAR_NOx_SHIP(I,J)
!      ENDIF
!%%%%%%%%%%%%% END KLUDGE %%%%%%%%%%%%%%
      
      
      ! Apply scale factor for future emissions (if necessary)
      IF ( LFUTURE ) THEN
         NOx = NOx * GET_FUTURE_SCALE_NOxff( I, J )
      ENDIF

      ! Convert units (if necessary)
      IF ( DO_KGS ) THEN

         ! Convert to [kg NOx/s]
         NOx = NOx / SEC_IN_SEASON

      ELSE IF ( DO_MCS ) THEN

         ! Convert to [molec/cm2/s]
         ! Convert NOx using XNUMOL for NO2 since original NOx
         ! emissions are in mass units of NO2 (mpayer, 4/18/13)
         NOx = NOx * XNUMOL_NO2 
     &       / ( GET_AREA_CM2( I, J, 1 ) * SEC_IN_SEASON )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EDGAR_NOx

!------------------------------------------------------------------------------

      FUNCTION GET_EDGAR_CO( I, J, KG_S, MOLEC_CM2_S, SHIP ) 
     &         RESULT( CO )
!
!******************************************************************************
!  Function GET_EDGAR_CO returns the EDGAR CO emissions at grid box (I,J)
!  in units of [kg/yr], [kg/s], or [molec/cm2/s]. (avd, bmy, 7/14/06)
!
!  Arguments: 
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude  index
!  (3 ) KG_S        (LOGICAL) : OPTIONAL - Return data in [kg/s]
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL - Return data in [molec/cm2/s]
!
!  NOTES:
!  (1 ) Added the SHIP switch (phs, 5/13/08)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff

      ! Arguments
      INTEGER, INTENT(IN)           :: I,      J
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S,   MOLEC_CM2_S, SHIP

      ! Local variables
      LOGICAL                       :: DO_KGS, DO_MCS
      REAL*8                        :: CO
      
      !=================================================================
      ! GET_EDGAR_CO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.

      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Get CO [kg/yr]
      IF ( PRESENT( SHIP ) ) THEN
         IF ( SHIP ) THEN
            CO = EDGAR_CO_SHIP(I,J)
         ELSE
            CO = EDGAR_CO(I,J)
         ENDIF
      ELSE
         CO = EDGAR_CO(I,J)
      ENDIF

      ! Apply scale factor for future emissions (if necessary)
      IF ( LFUTURE ) THEN
         CO = CO * GET_FUTURE_SCALE_COff( I, J )
      ENDIF

      ! Convert units (if necessary)
      IF ( DO_KGS ) THEN

         ! Convert to [kg CO/s]
         CO = CO / SEC_IN_2000

      ELSE IF ( DO_MCS ) THEN 

         ! Convert to [molec/cm2/s]
         CO = CO * XNUMOL_CO 
     &      / ( GET_AREA_CM2( I, J, 1 ) * SEC_IN_2000 )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EDGAR_CO

!------------------------------------------------------------------------------

      FUNCTION GET_EDGAR_ANTH_SO2( I, J, KG_S, MOLEC_CM2_S ) RESULT(SO2)
!
!******************************************************************************
!  Function GET_EDGAR_ANTH_SO2 returns the EDGAR anthropogenic SO2 emissions 
!  at grid box (I,J) in either [kg/yr] or [molec/cm2/s]. (avd, bmy, 7/14/06)
!
!  Arguments: 
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude  index
!  (3 ) KG_S        (LOGICAL) : OPTIONAL - Return data in [kg/s]
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL - Return data in [molec/cm2/s]
!
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

      ! Arguments
      INTEGER, INTENT(IN)           :: I,      J
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S,   MOLEC_CM2_S 

      ! Local variables
      LOGICAL                       :: DO_KGS, DO_MCS
      REAL*8                        :: SO2
      
      !=================================================================
      ! GET_EDGAR_ANTH_SO2 begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.

      ! Return data in [kg SO2/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Get anthropogenic SO2 [kg SO2/season]
      SO2 = EDGAR_SO2(I,J)

      ! Apply scale factor for future emissions (if necessary)
      IF ( LFUTURE ) THEN
         SO2 = SO2 * GET_FUTURE_SCALE_SO2ff( I, J )
      ENDIF

      ! Convert units
      IF ( DO_KGS ) THEN
         
         ! Convert to [kg SO2/s]
         SO2 = SO2 / SEC_IN_SEASON
     
      ELSE IF ( DO_MCS ) THEN 

         ! Convert to [molec/cm2/s]
         SO2 = SO2 * XNUMOL_SO2 
     &       / ( GET_AREA_CM2( I, J, 1 ) * SEC_IN_SEASON )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EDGAR_ANTH_SO2

!------------------------------------------------------------------------------

      FUNCTION GET_EDGAR_SHIP_SO2( I, J, KG_S, MOLEC_CM2_S ) RESULT(SO2)
!
!******************************************************************************
!  Function GET_EDGAR_SHIP_SO2 returns the EDGAR ship exhaust SO2 emissions 
!  at grid box (I,J) in units of [kg/month], [kg/s] or [molec/cm2/s]. 
!  (avd, bmy, 7/14/06)
!
!  Arguments: 
!  ============================================================================
!  (1 ) I           (INTEGER) : GEOS-Chem longitude index
!  (2 ) J           (INTEGER) : GEOS-Chem latitude  index
!  (3 ) KG_S        (LOGICAL) : OPTIONAL - Return data in [kg/s]
!  (4 ) MOLEC_CM2_S (LOGICAL) : OPTIONAL - Return data in [molec/cm2/s]
!
!  NOTES:
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

      ! Arguments
      INTEGER, INTENT(IN)           :: I,      J
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S,   MOLEC_CM2_S 

      ! Function value
      LOGICAL                       :: DO_KGS, DO_MCS
      REAL*8                        :: SO2
      
      !=================================================================
      ! GET_EDGAR_SHIP_SO2 begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.

      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Get ship SO2 [kg SO2/month]
      SO2 = EDGAR_SO2_SHIP(I,J)

      ! Apply scale factor for future emissions (if necessary)
      IF ( LFUTURE ) THEN
         SO2 = SO2 * GET_FUTURE_SCALE_SO2ff( I, J )
      ENDIF

      ! Convert units (if necessary)
      IF ( DO_KGS ) THEN

         ! Convert to [kg SO2/s]
         SO2 = SO2 / SEC_IN_MONTH

      ELSE IF ( DO_MCS ) THEN 

         ! Convert to [molec/cm2/s]
         SO2 = SO2 * XNUMOL_SO2 
     &       / ( GET_AREA_CM2( I, J, 1 ) * SEC_IN_MONTH )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EDGAR_SHIP_SO2

!------------------------------------------------------------------------------

      FUNCTION GET_EDGAR_TODN( I, J, HOUR ) RESULT( TODN )
!
!******************************************************************************
!  Function GET_EDGAR_TODN returns the time-of-day diurnal scale factor for
!  the EDGAR NOx emissions. (avd, bmy, 7/14/06)
!
!  Arguments: 
!  ============================================================================
!  (1 ) I    (INTEGER) : GEOS-Chem longitude index
!  (2 ) J    (INTEGER) : GEOS-Chem latitude  index
!  (3 ) HOUR (INTEGER) : GMT hour of the day (0-23)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, HOUR

      ! Local variables
      INTEGER             :: H
      REAL*8              :: TODN
      
      !=================================================================
      ! GET_EDGAR_TODN begins here!
      !=================================================================

      ! The 1st element of the array is hour 0, so add 1
      H    = HOUR + 1
      
      ! Get time of day factor for NOx
      TODN = EDGAR_TODN(I,J,H)

      ! Return to calling program
      END FUNCTION GET_EDGAR_TODN

!------------------------------------------------------------------------------
      
      SUBROUTINE EDGAR_TOTAL_Tg( YEAR, MONTH )
!
!******************************************************************************
!  Subroutine EDGAR_TOTAL_Tg prints totals of EDGAR emissions species
!  in units of Tg. (avd, bmy, 7/14/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR  (INTEGER) : Current year
!  (2 ) MONTH (INTEGER) : Current month
!
!  NOTES:
!  (1) Bug fix: now add IF (LEDGARSHIP) block before referencing any of the
!       ship emissions arrays.  This will help prevent a seg fault error.
!       (jensberg, bmy, 8/4/10)
!  08 Mar 2012 - M. Payer    - Added modifications for historical emissions of
!                              NOx and SO2 (E. Leibensperger)
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!******************************************************************************
!
![eml
      USE TIME_MOD,    ONLY : GET_HISTYR
!eml]

      ! Arguments
      INTEGER, INTENT(IN) :: YEAR, MONTH

      ! Local variables
      REAL*8              :: T_NOx, T_CO, T_SO2an, T_SO2sh
      REAL*8              :: T_NOxSh, T_COsh

      !=================================================================
      ! EDGAR_TOTAL_Tg begins here!
      !=================================================================

      ! Compute totals [Tg/yr]
      T_NOx      = SUM( EDGAR_NOx      ) * ( 14d0/46d0 ) / 1d9 ! Tg N
      T_CO       = SUM( EDGAR_CO       )                 / 1d9 ! Tg CO
      T_SO2an    = SUM( EDGAR_SO2      ) * ( 32d0/64d0 ) / 1d9 ! Tg S

      IF ( LEDGARSHIP ) THEN
         T_NOxSh = SUM( EDGAR_NOx_SHIP ) * ( 14d0/46d0 ) / 1d9 ! Tg N
         T_COSh  = SUM( EDGAR_CO_SHIP  )                 / 1d9 ! Tg CO
         T_SO2sh = SUM( EDGAR_SO2_SHIP ) * ( 32d0/64d0 ) / 1d9 ! Tg S
      ELSE
         T_NOxSh = 0d0
         T_COSh  = 0d0
         T_SO2sh = 0d0
      ENDIF

      ! Print totals
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'E D G A R   E M I S S I O N S'
![eml
      IF( .not. LHIST ) THEN
         WRITE( 6, 100     ) YEAR, SEASON_NAME, T_NOx
      ELSE
         WRITE( 6, 100     ) GET_HISTYR(), SEASON_NAME, T_NOx
      ENDIF
!eml]
      WRITE( 6, 105     ) YEAR, SEASON_NAME, T_NOxSh
      WRITE( 6, 110     ) YEAR,              T_CO
      WRITE( 6, 115     ) YEAR,              T_COSh
![eml
      IF( .not. LHIST ) THEN
         WRITE( 6, 120     ) YEAR, SEASON_NAME, T_SO2an
      ELSE
         WRITE( 6, 100     ) GET_HISTYR(), SEASON_NAME, T_SO2an
      ENDIF
!eml]

      WRITE( 6, 130     ) YEAR, MONTH_NAME,  T_SO2sh
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! FORMAT statements
 100  FORMAT( 'NO         for year ', i4, ' and season ', a3, 
     &         ' : ', f10.4, ' [Tg N    ]' )
 105  FORMAT( 'NO  Ship   for year ', i4, ' and season ', a3, 
     &         ' : ', f10.4, ' [Tg N    ]' )
 110  FORMAT( 'CO         for year ', i4, ' (annual total)',
     &         ' : ', f10.4, ' [Tg CO/yr]' )
 115  FORMAT( 'CO Ship    for year ', i4, ' (annual total)',
     &         ' : ', f10.4, ' [Tg CO/yr]' )
 120  FORMAT( 'Anthro SO2 for year ', i4, ' and season ', a3, 
     &         ' : ', f10.4, ' [Tg S    ]' )
 130  FORMAT( 'Ship   SO2 for year ', i4, ' and month  ', a3,
     &         ' : ', f10.4, ' [Tg S    ]' )

      ! Return to calling program
      END SUBROUTINE EDGAR_TOTAL_Tg

!------------------------------------------------------------------------------
      SUBROUTINE READ_AROMATICS( E_BENZ, E_TOLU, E_XYLE )
!
!******************************************************************************
!  Subroutine READ_AROMATICS now reads EDGARv2 inventory for Year 1985
!   from a bpch file (tmf, 1/8/08).
!
!  Arguments as Output
!  ============================================================================
!  (1 ) E_BENZ (REAL*4) : EDGARv2 anthro BENZ for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s]  (tmf, 7/8/05)
!  (2 ) E_TOLU (REAL*4) : EDGARv2 anthro TOLU for year 1985 
!                         (no seasonality, 1 level )
!                         [molec/cm2/s]  (tmf, 7/8/05)
!  (3 ) E_XYLE (REAL*4) : EDGARv2 anthro XYLE for year 1985 
!                         (no seasonality, 1 level )
!                         [molec/cm2/s]  (tmf, 7/8/05)
!
!  NOTES:
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT)     :: E_BENZ(IIPAR,JJPAR )
      REAL*4, INTENT(OUT)     :: E_TOLU(IIPAR,JJPAR )
      REAL*4, INTENT(OUT)     :: E_XYLE(IIPAR,JJPAR )

      ! Local variables
      INTEGER                 :: L
      REAL*4                  :: ARRAY(I1x1,J1x1 )
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)
      REAL*8                  :: INGRID(I1x1,J1x1)

      !=================================================================
      ! READ_AROMATICS begins here!
      !=================================================================
      ! Zero emission arrays
      E_BENZ = 0e0      
      E_TOLU = 0e0      
      E_XYLE = 0e0      

      ! Define the binary punch file name

      !=================================================================
      ! Read BENZ (tracer #57): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR_1x1 )//
     &         'EDGAR_200607/BENZ/BENZ_1985_FF_IND_EDGAR2.1x1geos.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',  57, 
     &                 0.d0,     I1x1,        J1x1,     
     &                 1,        ARRAY,       QUIET=.TRUE. )
         
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid BENZ from GEOS 1x1 to current model resolution
      ! [atoms C/cm2/s] (ccc, 3/9/09)
      INGRID = ARRAY
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                    INGRID,     OUTGRID, IS_MASS=0, 
     &                    netCDF=.TRUE.                   )
      E_BENZ = OUTGRID

      !=================================================================
      ! Read TOLU (tracer #58): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR_1x1 )//
     &         'EDGAR_200607/TOLU/TOLU_1985_FF_IND_EDGAR2.1x1geos.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',  58, 
     &                 0.d0,     I1x1,        J1x1,     
     &                 1,        ARRAY,       QUIET=.TRUE. )
         
      ! Regrid TOLU from GEOS 1x1 to current model resolution
      ! [atoms C/cm2/s] (ccc, 3/9/09)
      INGRID = ARRAY
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1, 
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )
      E_TOLU = OUTGRID

      !=================================================================
      ! Read XYLE (tracer #59): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR_1x1 )//
     &         'EDGAR_200607/XYLE/XYLE_1985_FF_IND_EDGAR2.1x1geos.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',  59, 
     &                 0.d0,     I1x1,        J1x1,     
     &                 1,        ARRAY,       QUIET=.TRUE. )
         
      ! Regrid XYLE from GEOS 1x1 to current model resolution
      ! [atoms C/cm2/s] (ccc, 3/9/09)
      INGRID = ARRAY
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )
      E_XYLE = OUTGRID


 100  FORMAT( 'READ_AROMATICS: Reading ', a )

      ! Return to calling program
      END SUBROUTINE READ_AROMATICS

!------------------------------------------------------------------------------

      SUBROUTINE READ_C2H4( E_C2H4 )
!
!******************************************************************************
!  Subroutine READ_C2H4 now reads EDGARv2 inventory for Year 1985 from a bpch file 
!  (tmf, 11/10/2006).
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_C2H4 (REAL*4) : EDGARv2 anthro C2H4 for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s] (tmf, 7/8/05)
!
!  NOTES:
!  (1) Bug fix: EDGAR200607 should be EDGAR_200607 (tmf, bmy, 12/20/10)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT)     :: E_C2H4(IIPAR,JJPAR )

      ! Local variables
      INTEGER                 :: L
      REAL*4                  :: ARRAY(I1x1,J1x1 )
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)
      REAL*8                  :: INGRID(I1x1,J1x1)

      !=================================================================
      ! READ_C2H4 begins here!
      !=================================================================
      ! Zero emission arrays
      E_C2H4 = 0e0      

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &        'EDGAR_200607/C2H4/C2H4_1985_FF_IND_EDGAR2.1x1geos.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_C2H4: Reading ', a )

      !=================================================================
      ! Read C2H4 (tracer #63): aseasonal
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',  63, 
     &                 0.d0,     I1x1,        J1x1,     
     &                 1,        ARRAY,       QUIET=.TRUE. )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid C2H4 from GEOS 1x1 to current model resolution
      ! [atoms C/cm2/s] (ccc, 3/9/09)
      INGRID = ARRAY
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )

      E_C2H4 = OUTGRID

      ! Return to calling program
      END SUBROUTINE READ_C2H4

!------------------------------------------------------------------------------

      SUBROUTINE READ_C2H2( E_C2H2 )
!
!******************************************************************************
!  Subroutine READ_C2H2 now reads EDGARv2 C2H2 emissions for year 1985
!     (tmf, 11/10/2006).
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_C2H2 (REAL*4) : EDGARv2 anthro C2H2 for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s] (tmf, 7/8/05)
!
!  NOTES:
!  (1) Bug fix: EDGAR200607 should be EDGAR_200607 (tmf, bmy, 12/20/10)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a algorithm
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD  
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT)     :: E_C2H2(IIPAR,JJPAR )

      ! Local variables
      INTEGER                 :: L
      REAL*4                  :: ARRAY(I1x1,J1x1 )
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)
      REAL*8                  :: INGRID(I1x1,J1x1)

      !=================================================================
      ! READ_C2H2 begins here!
      !=================================================================
      ! Zero emission arrays
      E_C2H2 = 0e0      

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &        'EDGAR_200607/C2H2/C2H2_1985_FF_IND_EDGAR2.1x1geos.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_C2H2: Reading ', a )

      !=================================================================
      ! Read C2H2 (tracer #64): aseasonal
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',  64, 
     &                 0.d0,     I1x1,        J1x1,     
     &                 1,        ARRAY,       QUIET=.TRUE. )
         
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid C2H2 from GEOS 1x1 to current model resolution
      ! [atoms C/cm2/s] (ccc, 3/9/09)
      INGRID = ARRAY
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1, 
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )

      E_C2H2 = OUTGRID

      ! Return to calling program
      END SUBROUTINE READ_C2H2

!------------------------------------------------------------------------------
      SUBROUTINE READ_AROMATICS_05x0666( E_BENZ, E_TOLU, E_XYLE )
!
!******************************************************************************
!  Subroutine READ_AROMATICS now reads EDGARv2 inventory for Year 1985
!   from a bpch file (tmf, 1/8/08).
!
!  Arguments as Output
!  ============================================================================
!  (1 ) E_BENZ (REAL*4) : EDGARv2 anthro BENZ for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s]  (tmf, 7/8/05)
!  (2 ) E_TOLU (REAL*4) : EDGARv2 anthro TOLU for year 1985 
!                         (no seasonality, 1 level )
!                         [molec/cm2/s]  (tmf, 7/8/05)
!  (3 ) E_XYLE (REAL*4) : EDGARv2 anthro XYLE for year 1985 
!                         (no seasonality, 1 level )
!                         [molec/cm2/s]  (tmf, 7/8/05)
!
!  NOTES:
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT) :: E_BENZ(IIPAR,JJPAR )
      REAL*4, INTENT(OUT) :: E_TOLU(IIPAR,JJPAR )
      REAL*4, INTENT(OUT) :: E_XYLE(IIPAR,JJPAR )

      ! Local variables
      INTEGER                       :: L
      REAL*4                        :: ARRAY(IIPAR,JJPAR )
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_AROMATICS begins here!
      !=================================================================
      ! Zero emission arrays
      E_BENZ = 0e0      
      E_TOLU = 0e0      
      E_XYLE = 0e0      

      ! Define the binary punch file name

      !=================================================================
      ! Read BENZ (tracer #57): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) //
     &        'EDGAR_200905/BENZ/BENZ_1985_FF_IND_EDGAR2.05x0666CH.bpch'
      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    57, 
     &                 0.d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY,        QUIET=.TRUE. )
         
      E_BENZ = ARRAY

      !=================================================================
      ! Read TOLU (tracer #58): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) //
     &        'EDGAR_200905/TOLU/TOLU_1985_FF_IND_EDGAR2.05x0666CH.bpch'
      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    58, 
     &                 0.d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY,        QUIET=.TRUE. )
         
      E_TOLU = ARRAY

      !=================================================================
      ! Read XYLE (tracer #59): aseasonal
      !=================================================================
      FILENAME = TRIM( DATA_DIR ) //
     &        'EDGAR_200905/XYLE/XYLE_1985_FF_IND_EDGAR2.05x0666CH.bpch'
      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )


      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    59, 
     &                 0.d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY,        QUIET=.TRUE. )
         
      E_XYLE = ARRAY

 100  FORMAT( 'READ_AROMATICS: Reading ', a )


      ! Return to calling program
      END SUBROUTINE READ_AROMATICS_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE READ_C2H4_05x0666( E_C2H4 )
!
!******************************************************************************
!  Subroutine READ_C2H4 now reads EDGARv2 inventory for Year 1985 from a bpch file 
!  (tmf, 11/10/2006).
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_C2H4 (REAL*4) : EDGARv2 anthro C2H4 for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s] (tmf, 7/8/05)
!
!  NOTES:
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT) :: E_C2H4(IIPAR,JJPAR )

      ! Local variables
      INTEGER                       :: L
      REAL*4                        :: ARRAY(IIPAR,JJPAR )
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_C2H4 begins here!
      !=================================================================
      ! Zero emission arrays
      E_C2H4 = 0e0      

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR ) //
     &        'EDGAR_200905/C2H4/C2H4_1985_FF_IND_EDGAR2.05x0666CH.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_C2H4: Reading ', a )

      !=================================================================
      ! Read C2H4 (tracer #63): aseasonal
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    63, 
     &                 0.d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY,        QUIET=.TRUE. )
         
      E_C2H4 = ARRAY


      ! Return to calling program
      END SUBROUTINE READ_C2H4_05x0666

!------------------------------------------------------------------------------

      SUBROUTINE READ_C2H2_05x0666( E_C2H2 )
!
!******************************************************************************
!  Subroutine READ_C2H2 now reads EDGARv2 C2H2 emissions for year 1985
!     (tmf, 11/10/2006).
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) E_C2H2 (REAL*4) : EDGARv2 anthro C2H2 for year 1985 
!                         (no seasonality, 1 level ) 
!                         [molec/cm2/s] (tmf, 7/8/05)
!
!  NOTES:
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,      ONLY : GET_MONTH
      USE BPCH2_MOD  
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE CMN_SIZE_MOD            ! Size parameters

      ! Arguments
      REAL*4, INTENT(OUT) :: E_C2H2(IIPAR,JJPAR )

      ! Local variables
      INTEGER                       :: L
      REAL*4                        :: ARRAY(IIPAR,JJPAR )
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_C2H2 begins here!
      !=================================================================
      ! Zero emission arrays
      E_C2H2 = 0e0      

      ! Define the binary punch file name
      FILENAME = TRIM( DATA_DIR ) //
     &        'EDGAR_200905/C2H2/C2H2_1985_FF_IND_EDGAR2.05x0666CH.bpch'

      ! Write file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_C2H2: Reading ', a )

      !=================================================================
      ! Read C2H2 (tracer #64): aseasonal
      !=================================================================
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',    64, 
     &                 0.d0,       IIPAR,        JJPAR,     
     &                    1,         ARRAY,        QUIET=.TRUE. )
         
      E_C2H2 = ARRAY

      ! Return to calling program
      END SUBROUTINE READ_C2H2_05x0666

!------------------------------------------------------------------------------
! NOTE: This should be for debugging... 
!      SUBROUTINE OUTPUT_TOTAL_2D( DESCRIPTION, EMISSIONS, UNITS )
!!
!!******************************************************************************
!!  Subroutine OUTPUT_TOTAL outputs the total emissions for a given emissions
!!  array. (amv 02/14/06)
!!
!!  NOTES:
!!******************************************************************************
!
!      USE TIME_MOD,        ONLY : GET_TS_EMIS
!
!      USE CMN_SIZE_MOD  ! size parameters
!
!      ! Local variables
!      INTEGER                        :: I, J, K, IMAX, IMIN, JMAX, JMIN
!      REAL*8                         :: TOTAL, SCALAR
!      CHARACTER(LEN=255)             :: LOCATION
!      CHARACTER(LEN=255)             :: OUTUNITS
!
!      ! Arguments
!      CHARACTER(LEN=*), INTENT(IN)   :: DESCRIPTION
!      CHARACTER(LEN=*), INTENT(IN)   :: UNITS
!      REAL*8, INTENT(IN)             :: EMISSIONS(IIPAR, JJPAR)
!
!      ! Associate output units with input description
!      IF ( TRIM( DESCRIPTION ) == 'EDGAR NOx' ) THEN
!         !OUTUNITS = '[Tg N/yr]'
!         OUTUNITS = '[Tg N/season]'
!      ELSEIF ( TRIM( DESCRIPTION ) == 'GEIA NOx' ) THEN
!         OUTUNITS = '[Tg N/yr]'
!      ELSEIF ( TRIM( DESCRIPTION ) == 'EDGAR CO' ) THEN
!         OUTUNITS = '[Tg CO/yr]'
!      ELSEIF ( TRIM( DESCRIPTION ) == 'GEIA CO' ) THEN
!         OUTUNITS = '[Tg CO/yr]'
!      ELSEIF ( TRIM( DESCRIPTION ) == 'GEIA SO2' ) THEN
!         OUTUNITS = '[Tg S/yr]'
!      ELSEIF ( TRIM( DESCRIPTION ) == 'EDGAR SO2' ) THEN
!         !OUTUNITS = '[Tg SO2/yr]'
!         OUTUNITS = '[Tg S/season]'
!      ELSE
!         OUTUNITS = 'Unknown'
!      ENDIF
!
!      ! Scalar to convert to kg/yr from received units
!      IF ( TRIM( UNITS ) == 'kg/yr' ) THEN
!         SCALAR = 1.0 / 1.0d9
!      ELSEIF ( TRIM( UNITS ) == 'kg/s' ) THEN
!         SCALAR = 1.d0/ 1.0d9 * 365.25 * 24.0 * 60.0 * 60.0
!      ELSEIF ( TRIM( UNITS ) == 'kg/ts' ) THEN
!         SCALAR = 1.0 / 1.0d9 / GET_TS_EMIS() * 365.25 * 24.0 * 60.0
!      ELSEIF ( TRIM( UNITS ) == 'kg/season' ) THEN
!         SCALAR   = 1.0 / 1d9
!      ELSE
!         SCALAR = 1.0
!      ENDIF
!
!      ! Extra conversion from NO2 -> N
!      ! (NO2 is stored for NOx)
!      IF ( TRIM( DESCRIPTION ) == 'EDGAR NOx' .or. 
!     &     TRIM( DESCRIPTION ) == 'GEIA NOx' ) THEN
!         SCALAR = SCALAR * 14. / 46.
!      ENDIF
!
!      ! Extra conversion from SO2 -> S
!      ! (NO2 is stored for NOx)
!      IF ( TRIM( DESCRIPTION ) == 'EDGAR NOx' .or. 
!     &     TRIM( DESCRIPTION ) == 'GEIA NOx' ) THEN
!         SCALAR = SCALAR * 32. / 64.
!      ENDIF
!
!      ! loop over each region
!      DO K =  1,6
!
!         ! reset total
!         TOTAL = 0.0
!
!         IF ( K == 1 ) THEN
!            LOCATION = '        World'
!            IMIN = 1
!            IMAX = IIPAR
!            JMIN = 1
!            ! avoid anomylous points at pole
!            IF ( TRIM(DESCRIPTION) == 'GEIA SO2') THEN
!               JMAX = NINT( 170. / 180. * JJPAR )
!            ELSE
!               JMAX = JJPAR
!            ENDIF
!         ELSEIF ( K == 2 ) THEN
!            LOCATION = 'North America'
!            IMIN = NINT(  15. / 360. * IIPAR )
!            IMAX = NINT( 140. / 360. * IIPAR )
!            JMIN = NINT( 110. / 180. * JJPAR )
!            JMAX = NINT( 170. / 180. * JJPAR )
!         ELSEIF ( K == 3 ) THEN
!            LOCATION = 'South America'
!            IMIN = NINT(  90. / 360. * IIPAR )
!            IMAX = NINT( 150. / 360. * IIPAR )
!            JMIN = NINT(  30. / 180. * JJPAR )
!            JMAX = NINT( 105. / 180. * JJPAR )
!         ELSEIF ( K == 4 ) THEN
!            LOCATION = '       Europe'
!            IMIN = NINT( 165. / 360. * IIPAR )
!            IMAX = NINT( 240. / 360. * IIPAR )
!            JMIN = NINT( 125. / 180. * JJPAR )
!            JMAX = NINT( 170. / 180. * JJPAR )
!         ELSEIF ( K == 5 ) THEN
!            LOCATION = '         Asia'
!            IMIN = NINT( 240. / 360. * IIPAR )
!            IMAX = NINT( 350. / 360. * IIPAR )
!            JMIN = NINT(  90. / 180. * JJPAR )
!            JMAX = NINT( 170. / 180. * JJPAR )
!         ELSEIF ( K == 6 ) THEN
!            LOCATION = '       Africa'
!            IMIN = NINT( 160. / 360. * IIPAR )
!            IMAX = NINT( 235. / 360. * IIPAR )
!            JMIN = NINT(  50. / 180. * JJPAR )
!            JMAX = NINT( 125. / 180. * JJPAR )
!         ENDIF
!
!         DO I = IMIN,IMAX
!            DO J = JMIN,JMAX
!               TOTAL = TOTAL + EMISSIONS(I,J) * SCALAR
!            ENDDO
!         ENDDO
!
!         WRITE(6,'(a,a,a,a,a,F7.2,1x,a15)') 'Total ', TRIM(DESCRIPTION), 
!     &     ' Emissions, ', TRIM(LOCATION), ': ', TOTAL, 
!     &     TRIM(OUTUNITS)
!
!      ENDDO
!
!      ! Return to calling program
!      END SUBROUTINE OUTPUT_TOTAL_2D
!

!------------------------------------------------------------------------------

      SUBROUTINE INIT_EDGAR( Input_Opt )
!
!******************************************************************************
!  Subroutine INIT_EDGAR allocates and initializes all module arrays.
!  (avd, bmy, 7/14/06)
!
!  NOTES:
!  (1 ) Now allocates SHIP arrays only if LEDGARSHIP is TRUE (5/13/08)
!  01 Mar 2012 - R. Yantosca - Removed ref to GET_AREA_CM2 and A_CM2 array
!  12 Jun 2013 - R. Yantosca - Now use Input_Opt%LEDGARSHIP
!******************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE CMN_SIZE_MOD  ! Size parameters
      
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      TYPE(OptInput), INTENT(IN) :: Input_Opt

      ! Local Variables
      INTEGER            :: AS, J

      !============================================================
      ! INIT_EDGAR begins here!
      !============================================================

      ALLOCATE( EDGAR_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_NOx' )
      EDGAR_NOx = 0d0

      ALLOCATE( EDGAR_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_CO' )
      EDGAR_CO = 0d0

      ALLOCATE( EDGAR_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_SO2' )
      EDGAR_SO2 = 0d0

      IF ( Input_Opt%LEDGARSHIP ) THEN

         ALLOCATE( EDGAR_SO2_SHIP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_SO2_SHIP' )
         EDGAR_SO2_SHIP = 0d0

         ALLOCATE( EDGAR_NOx_SHIP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_NOx_SHIP' )
         EDGAR_NOx_SHIP = 0d0

         ALLOCATE( EDGAR_CO_SHIP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_CO_SHIP' )
         EDGAR_CO_SHIP = 0d0
 
      ENDIF

      ALLOCATE( EDGAR_TODN( IIPAR, JJPAR, N_HOURS ), STAT=AS )
      IF (AS /= 0 ) CALL ALLOC_ERR( 'EDGAR_TODN' )
      EDGAR_TODN = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_EDGAR

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_EDGAR
!
!******************************************************************************
!  Subroutine CLEANUP_EDGAR deallocates all module arrays, (avd, bmy, 7/14/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_EDGAR begins here!
      !=================================================================
      IF ( ALLOCATED( EDGAR_NOx      ) ) DEALLOCATE( EDGAR_NOx      )
      IF ( ALLOCATED( EDGAR_CO       ) ) DEALLOCATE( EDGAR_CO       )
      IF ( ALLOCATED( EDGAR_SO2      ) ) DEALLOCATE( EDGAR_SO2      )
      IF ( ALLOCATED( EDGAR_SO2_SHIP ) ) DEALLOCATE( EDGAR_SO2_SHIP )
      IF ( ALLOCATED( EDGAR_CO_SHIP  ) ) DEALLOCATE( EDGAR_CO_SHIP  )
      IF ( ALLOCATED( EDGAR_NOX_SHIP ) ) DEALLOCATE( EDGAR_NOx_SHIP )
      IF ( ALLOCATED( EDGAR_TODN     ) ) DEALLOCATE( EDGAR_TODN     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_EDGAR

!------------------------------------------------------------------------------

      ! End of module
      END MODULE EDGAR_MOD


