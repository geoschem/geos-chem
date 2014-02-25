!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vistas_anthro_mod
!
! !DESCRIPTION: Module VISTAS\_ANTHRO\_MOD contains variables and routines 
!  to read the VISTAS anthropogenic emissions.
!\\
!\\
! !INTERFACE: 
!
      MODULE VISTAS_ANTHRO_MOD
! 
! !USES:
!
      USE EPA_NEI_MOD, ONLY : GET_USA_MASK

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_VISTAS_ANTHRO
      PUBLIC :: EMISS_VISTAS_ANTHRO
      PUBLIC :: GET_VISTAS_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: INIT_VISTAS_ANTHRO
      PRIVATE :: VISTAS_SCALE_FUTURE
      PRIVATE :: TOTAL_ANTHRO_Tg
!
! !REVISION HISTORY:
!  24 Nov 2008 - A. v. Donkelaar - Initial version  
!  28 Jan 2009 - P. Le Sager     - Initial Version in GEOS-Chem
!  01 Mar 2012 - R. Yantosca     - Now reference new grid_mod.F90
!  14 Mar 2013 - M. Payer        - Replace NOx emissions with NO emissions as
!                                  part of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Arrays for weekday & weekend emissions
      REAL*8,  ALLOCATABLE :: VISTAS_WD_NOx(:,:)
      REAL*8,  ALLOCATABLE :: VISTAS_WE_NOx(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_vistas_anthro
!
! !DESCRIPTION: Function GET\_VISTAS\_ANTHRO returns the VISTAS emission for 
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_VISTAS_ANTHRO( I,        J,           N,  
     &                            WEEKDAY,  MOLEC_CM2_S, KG_S ) 
     &  RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNO
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N

      ! Return weekday or weekend emissions
      LOGICAL, INTENT(IN)           :: WEEKDAY

      ! OPTIONAL -- return emissions in [molec/cm2/s]
      LOGICAL, INTENT(IN), OPTIONAL :: MOLEC_CM2_S

      ! OPTIONAL -- return emissions in [kg/s]
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
!
! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE    
!
! !REVISION HISTORY: 
!  28 Jan 2009 - A. v. Donkelaar, P. Le Sager - Initial Version
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS

      !=================================================================
      ! GET_VISTA_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S

      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( N == IDTNO ) THEN

         ! NOx [molec/cm2/s]
         IF ( WEEKDAY ) THEN
            VALUE = VISTAS_WD_NOx(I,J)
         ELSE
            VALUE = VISTAS_WE_NOx(I,J)
         ENDIF

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no VISTAS emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [molec/c,2/s] to [kg/s]
         VALUE = VALUE * GET_AREA_CM2( I, J, 1 ) / XNUMOL(N)

      ENDIF

      END FUNCTION GET_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_vistas_anthro
!
! !DESCRIPTION: Subroutine EMISS\_VISTAS\_ANTHRO reads the VISTAS emission 
!  fields at 1x1 resolution and regrids them to the current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_VISTAS_ANTHRO( am_I_Root, Input_Opt, 
     &                                State_Chm, RC         )
!
! !USES:
! 
      USE BPCH2_MOD,          ONLY : GET_TAU0,      READ_BPCH2
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR_1x1
#if defined( DEVEL )
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK_LT
      USE TRACERID_MOD,       ONLY : IDTNO
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE ERROR_MOD,          ONLY : ALLOC_ERR
#endif
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
!  28 Jan 2009 - A. v. Donkelaar, P. Le Sager - Initial Version
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fix minor bugs in map_a2a algorithm
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I, J, THISYEAR
      INTEGER               :: MN, SNo, ScNo
      REAL*4                :: ARRAY(I1x1,J1x1,1)
      REAL*8                :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                :: SC_1x1(I1x1,J1x1)
      REAL*8                :: TAU2002, TAU
      CHARACTER(LEN=255)    :: FILENAME, VISTAS_DIR
      CHARACTER(LEN=4)      :: SYEAR, SNAME
      CHARACTER(LEN=2)      :: SMN
      CHARACTER(LEN=1)      :: SSMN
      CHARACTER(LEN=255)    :: LLFILENAME
      REAL*8                :: INGRID(I1x1,J1x1)
#if defined( DEVEL )
      INTEGER               :: DOW_LT
      LOGICAL               :: WEEKDAY
      REAL*8                :: A(IIPAR,JJPAR)
#endif

      ! For fields from Input_Opt
      LOGICAL               :: LFUTURE
      INTEGER               :: N_TRACERS
      REAL*8                :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! EMISS_VISTAS_ANTHRO begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_VISTAS_ANTHRO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      VISTAS_DIR = TRIM( DATA_DIR_1x1 ) // 'VISTAS_200811/'

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! cap maximum scaling year
      IF ( THISYEAR .gt. 2010 ) THEN
         THISYEAR = 2010
      ENDIF

      SNAME = 'NOx'
      SNo = 1
      ScNo = 71
            
      TAU2002 = GET_TAU0( 1, 1, 2002)
      MN = GET_MONTH()

      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'Vistas-' // TRIM(SNAME) // '-'
     &            // SMN // '.1x1'
      ENDIF

      WRITE( SYEAR, '(i4)' ) THISYEAR


      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_VISTAS_ANTHRO: Reading ', a )

      ! Read data [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Load ozone season regulation factors
      IF (MN .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) MN
         FILENAME  = TRIM( VISTAS_DIR )
     &            // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &            // SMN // '.1x1'
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'RATIO-2D', ScNo,
     &                    TAU2002,   I1x1,      J1x1,
     &                    1,         ARRAY,     QUIET=.TRUE. )

      ! Apply Ozone Season Scalars
      GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * ARRAY(:,:,1)

      ! Load/Apply weekend factors
      TAU = GET_TAU0( MN, 1, 1999)
      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkend_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 2,
     &                 TAU,      I1x1,      J1x1,
     &                 1,        ARRAY,     QUIET=.TRUE. )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid to current model resolution [molec/cm2/s]
      INGRID = GEOS_1x1(:,:,1) * ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,          J1x1,
     &                    INGRID,     VISTAS_WE_NOx, IS_MASS=0,
     &                    netCDF=.TRUE.                         )


      ! Load/Apply weekday factors
      FILENAME  = TRIM( VISTAS_DIR )
     &            // 'wkday_an_scalar.nei99.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'WD-WE-$', 1,
     &                 TAU,      I1x1,      J1x1,
     &                 1,        ARRAY,     QUIET=.TRUE. )

 
      ! Regrid to current model resolution [molec/cm2/s]
      INGRID = GEOS_1x1(:,:,1) * ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,          J1x1,
     &                    INGRID,     VISTAS_WD_NOx, IS_MASS=0,
     &                    netCDF=.TRUE.                         )

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL VISTAS_SCALE_FUTURE
      ENDIF

#if defined( DEVEL )
      !--------------------------
      ! Give to global bundle
      !--------------------------
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Determine if we should use weekday or weekend NEI
         ! emissions at grid box (I,J,L).  Since NEI is over
         ! the US, then weekend is Sat/Sun.
         DOW_LT  = GET_DAY_OF_WEEK_LT( I, J, 1 )
         WEEKDAY = ( DOW_LT > 0 .and. DOW_LT < 6 )

         ! Surface area
         A(I,J)  = GET_AREA_CM2( I, J, 1 ) 

         ! Sum into the tracer tendencies
         IF ( WEEKDAY ) THEN
            IF ( IDTNO > 0 ) THEN
               State_Chm%TRAC_TEND(I,J,1,IDTNO) =
     &           VISTAS_WD_NOx(I,J) * A(I,J) / XNUMOL(IDTNO)
            ENDIF
         ELSE
            IF ( IDTNO > 0 ) THEN
               State_Chm%TRAC_TEND(I,J,1,IDTNO) =
     &           VISTAS_WE_NOx(I,J) * A(I,J) / XNUMOL(IDTNO)
            ENDIF
         ENDIF
      ENDDO
      ENDDO
#endif

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR, MN, Input_Opt )

      END SUBROUTINE EMISS_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vistas_scale_future
!
! !DESCRIPTION: Subroutine VISTAS\_SCALE\_FUTURE applies the IPCC future scale 
!  factors to  the VISTAS anthropogenic emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE VISTAS_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff

      USE CMN_SIZE_MOD             ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - A. v. Donkelaar, P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! VISTAS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         VISTAS_WE_NOx(I,J)  = VISTAS_WE_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )
         VISTAS_WD_NOx(I,J)  = VISTAS_WD_NOx(I,J) 
     &                         * GET_FUTURE_SCALE_NOxff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE VISTAS_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR, THISMONTH, Input_Opt )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE TRACERID_MOD,       ONLY : IDTNO
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: YEAR        ! Year to compute totals
      INTEGER,        INTENT(IN) :: THISMONTH   ! Month to compute totals
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  28 Mar 2013 - R. Yantosca - Now use fields from Input_Opt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      REAL*8              :: WD_NOX, WE_NOX, F_NOX, A
      CHARACTER(LEN=3)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================
      WD_NOX  = 0d0
      WE_NOX  = 0d0
      F_NOX   = Input_Opt%TRACER_MW_KG(IDTNO)

      ! Loop over latitudes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Surface area [cm2] * seconds in this month / AVOGADRO's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A       = GET_AREA_CM2( I, J, 1 ) 
     &           * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

         ! Weekday avg emissions
         WD_NOX  = WD_NOX  + VISTAS_WD_NOX (I,J) * A * F_NOX

         ! Weekend avg emissions
         WE_NOX  = WE_NOX  + VISTAS_WE_NOX (I,J) * A * F_NOX

      ENDDO
      ENDDO

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'VISTAS   U S A   E M I S S I O N S', / )

      ! Weekday avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 110   ) 'NO ', THISMONTH, WD_NOX,  '  '
 110  FORMAT( 'Total weekday avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      ! Weekend avg anthro
      WRITE( 6, '(a)' )
      WRITE( 6, 120   ) 'NO ', THISMONTH, WE_NOX,  '  '
 120  FORMAT( 'Total weekend avg anthro ', a4, ' for 1999/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      END SUBROUTINE TOTAL_ANTHRO_Tg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_vistas_anthro 
!
! !DESCRIPTION: Subroutine INIT\_VISTAS\_ANTHRO allocates and zeroes all 
!  module arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_VISTAS_ANTHRO( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!  01 Mar 2012 - R. Yantosca - Remove A_CM2 array, use GET_AREA_CM2 instead
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! INIT_VISTAS_ANTHRO begins here!
      !=================================================================

      ! Assume success
      RC  =  GIGC_SUCCESS

      ! Return if LVISTAS is false
      IF ( .not. Input_Opt%LVISTAS ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( VISTAS_WD_NOx( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'VISTAS_WD_NOx' )
      VISTAS_WD_NOx = 0d0

      ALLOCATE( VISTAS_WE_NOx( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'VISTAS_WE_NOx' )
      VISTAS_WE_NOx = 0d0

      END SUBROUTINE INIT_VISTAS_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_vistas_anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_VISTAS\_ANTHRO deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_VISTAS_ANTHRO
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!  02 Mar 2012 - R. Yantosca - Remove A_CM2 array
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_VISTAS_ANTHRO begins here!
      !=================================================================
      IF ( ALLOCATED( VISTAS_WD_NOx  ) ) DEALLOCATE( VISTAS_WD_NOx  )
      IF ( ALLOCATED( VISTAS_WE_NOx  ) ) DEALLOCATE( VISTAS_WE_NOx  )

      END SUBROUTINE CLEANUP_VISTAS_ANTHRO
!EOC
      END MODULE VISTAS_ANTHRO_MOD
