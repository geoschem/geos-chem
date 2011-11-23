! $Id: icoads_ship_mod.f,v 1.2 2010/02/23 20:55:44 bmy Exp $
!------------------------------------------------------------------------------
!     Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: icoads_ship_mod
!     
! !DESCRIPTION: Module ICOADS\_SHIP\_MOD contains variables and routines to 
!  read the International Comprehensive Ocean-Atmosphere Data Set (ICOADS)
!  ship emissions. Base year is 2002.
!\\
!\\     
!  Source: ICOADS Emissions data for NOx, SOx, and CO were downloaded from 
!  http://coast.cms.udel.edu/GlobalShipEmissions/Inventories/
!\\
!\\
!  Reference: Wang, C., J. J. Corbett, and J. Firestone, \emph{Improving
!  Spatial representation of Global Ship Emissions Inventories},
!  \underline{Environ. Sci. Technol.}, \textbf{42}, (1), 193-199, 2008.
!\\   
!\\   
! !INTERFACE: 
!     
      MODULE ICOADS_SHIP_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_ICOADS_SHIP
      PUBLIC :: EMISS_ICOADS_SHIP
      PUBLIC :: GET_ICOADS_SHIP
      PUBLIC :: INTERPOLATE_LUT
      PUBLIC :: INTERPOLATE_LUT2
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: ICOADS_SCALE_FUTURE
      PRIVATE :: INIT_ICOADS_SHIP
      PRIVATE :: TOTAL_ICOADS_SHIP_TG
!
! !REVISION HISTORY:
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!
      ! Array for surface area
      REAL*8,  ALLOCATABLE :: A_CM2(:)

      ! Arrays for emissions
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)

      CONTAINS
!
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_ICOADS_SHIP
!
! !DESCRIPTION: Function GET\_ICOADS\_SHIP returns the ICOADS ship emissions for
!  GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in units 
!  of [kg/s] or [molec/cm2/s].  (cklee, 7/09/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ICOADS_SHIP( I,    J,     N, 
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNOx, IDTCO, IDTSO2, IDTNH3
      USE TIME_MOD,     ONLY : GET_YEAR, GET_MONTH
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, N

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
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS, DO_MCS
      INTEGER                       :: YEAR, MONTH
      REAL*8                        :: SEC_IN_MONTH

      !=================================================================
      ! GET_ICOADS_SHIP begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      IF ( N == IDTNOx ) THEN

         ! NOx [kg/month]
         VALUE = NOx(I,J)

      ELSE IF ( N == IDTCO ) THEN

         ! CO [kg/month]
         VALUE = CO(I,J)

      ELSE IF ( N == IDTSO2 ) THEN

         ! SO2 [kg/month]
         VALUE = SO2(I,J)

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no CAC emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      ! Get emissions year
      YEAR = GET_YEAR()

      ! Get emissions month      
      MONTH = GET_MONTH()

      IF ( (MONTH == 4) .OR. (MONTH == 6) .OR.
     &   (MONTH == 9) .OR. (MONTH == 11) ) THEN

         SEC_IN_MONTH = 86400D0*30.0D0

      ELSE IF (MONTH == 2) THEN

         ! ICOADS ship emissions for 2002
         IF (MOD(YEAR,4) == 0) THEN
            SEC_IN_MONTH = 86400D0*29.0D0
         ELSE
            SEC_IN_MONTH = 86400D0*28.0D0
         ENDIF

      ELSE

         SEC_IN_MONTH = 86400D0*31.0D0

      ENDIF

      IF ( DO_KGS ) THEN
            
         ! Convert from [kg/box/month] to [kg/box/s]
         VALUE = VALUE / SEC_IN_MONTH

      ELSE IF ( DO_MCS ) THEN

         ! Convert NOx from [kg/month] to [molec/cm2/s]
         VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_MONTH )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISS_ICOADS_SHIP
!
! !DESCRIPTION: Subroutine EMISS\_ICOADS\_SHIP reads the ICOADS emission fields
!  at 1x1 resolution and regrids them to the current model resolution. 
!  (cklee, 7/09/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_ICOADS_SHIP
!
! !USES:
! 
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR,      GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_O3"            ! FSCALYR
!
! !REVISION HISTORY: 
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR, SPECIES, SNo, ScNo
      INTEGER                    :: THISMONTH
      REAL*4                     :: ARRAY(I1x1,J1x1,1)
      REAL*8                     :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                     :: SC_1x1(I1x1,J1x1)
      REAL*8                     :: TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER (LEN=2)          :: SMONTH


      !=================================================================
      ! EMISS_ICOADS_SHIP begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_ICOADS_SHIP
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! Get emissions month      
      THISMONTH = GET_MONTH()

      WRITE( SMONTH, '(i2.2)' ) THISMONTH


      DO SPECIES = 1,3

         IF ( SPECIES .eq. 1 ) THEN
            SNAME = 'NOx'
            SNo = 1
            ScNo = 71
         ELSEIF ( SPECIES .eq. 2 ) THEN
            SNAME = 'CO'
            SNo = 4
            ScNo = 72
         ELSEIF ( SPECIES .eq. 3 ) THEN
            SNAME = 'SOx'
            SNo = 26
            ScNo = 73
         ENDIF
            

         ! TAU values for 2002
         TAU = GET_TAU0( 1, 1, 2002 )

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) //'ICOADS_200907/' //
     &               TRIM( SNAME ) // '_' // SMONTH // '.geos.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - EMISS_ICOADS_SHIP: Reading ', a )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'ICOADS-$', SNo, 
     &                    TAU,      I1x1,       J1x1,     
     &                    1,        ARRAY,      QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEOS_1x1(:,:,1) = ARRAY(:,:,1)

         ! Convert [kg S/month] to [kg SO2/month]
         IF ( SPECIES .eq. 3 ) THEN
            GEOS_1X1 = GEOS_1x1*64.0D0/32.0D0
         ENDIF

         ! Apply annual scalar factor
         CALL GET_ANNUAL_SCALAR_1x1( ScNo, 2002, THISYEAR, SC_1x1 )

         GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)



         ! Regrid from GEOS 1x1 --> current model resolution
         IF ( SPECIES .eq. 1 ) THEN

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, NOx )

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, CO )

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * 0.986

            CALL DO_REGRID_1x1( 'kg/month', GEOS_1x1, SO2 )

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL ICOADS_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ICOADS_SHIP_TG( THISYEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ICOADS_SCALE_FUTURE
!
! !DESCRIPTION: applies the IPCC future scale factors
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE ICOADS_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

#     include "CMN_SIZE"             ! Size parameters
!
! !REVISION HISTORY: 
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! ICOADS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/month]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /month]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/month] 
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ICOADS_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_ICOADS_SHIP_TG
!
! !DESCRIPTION: Subroutine TOTAL\_ICOADS\_SHIP\_TG prints the totals for  
!   ship emissions of NOx, CO, and SO2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ICOADS_SHIP_TG( MONTH )
!
! !USES:
! 

#     include "CMN_SIZE"            ! Size parameters

!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: MONTH   ! Month of data to compute totals
!
! !REVISION HISTORY: 
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      REAL*8              :: T_NOX, T_CO,  T_SO2
      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ICOADS_SHIP_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'I. C. O. A. D. S.   S H I P   E M I S S I O N S', / )


      ! Total NOx [Tg N]
      T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )

      ! Total CO  [Tg CO]
      T_CO  = SUM( CO  ) * 1d-9

      ! Total SO2 [Tg S]
      T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

      ! Print totals in [kg]
      WRITE( 6, 110 ) 'NOx ', MONTH, T_NOx,  '[Tg N  ]'
      WRITE( 6, 110 ) 'CO  ', MONTH, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'SO2 ', MONTH, T_SO2,  '[Tg S  ]'

      ! Format statement
 110  FORMAT( 'ICOADS ship ', a5, 
     &        'for month ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_ICOADS_SHIP_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_ICOADS_SHIP
!
! !DESCRIPTION: Subroutine INIT\_ICOADS\_SHIP allocates and zeroes all 
!  module arrays. (cklee, 7/09/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ICOADS_SHIP
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LICOADSSHIP

#     include "CMN_SIZE"    ! Size parameters
!
! !REVISION HISTORY: 
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_ICOADS_SHIP begins here!
      !=================================================================

      ! Return if LICOADSSHIP is false
      IF ( .not. LICOADSSHIP ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ALLOCATE( NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      !---------------------------------------------------
      ! Pre-store array for grid box surface area in cm2
      !---------------------------------------------------

      ! Allocate array
      ALLOCATE( A_CM2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

      ! Fill array
      DO J = 1, JJPAR
         A_CM2(J) = GET_AREA_CM2( J )
      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_ICOADS_SHIP
!
! !DESCRIPTION:  Subroutine CLEANUP\_ICOADS\_SHIP deallocates all module arrays. 
!  (cklee, 7/09/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ICOADS_SHIP
!
! !REVISION HISTORY: 
!   21 Jul 2009 - Chulkyu Lee & P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( A_CM2          ) ) DEALLOCATE( A_CM2          )
      IF ( ALLOCATED( NOx            ) ) DEALLOCATE( NOx            )
      IF ( ALLOCATED( CO             ) ) DEALLOCATE( CO             )
      IF ( ALLOCATED( SO2            ) ) DEALLOCATE( SO2            )

      ! Return to calling program
      END SUBROUTINE CLEANUP_ICOADS_SHIP
!EOC
!------------------------------------------------------------------------------

      subroutine INTERPOLATE_LUT(I,J,fraction_nox,int_ope) 
    !=======================================================================
    !
    ! INTERPOLATE_LUT:    Return FracNOx or IntOPE from the lookup table
    !
    ! temp   : model temperature
    ! jno2   : J(NO2) value
    ! cao3   : concentration O3 in ambient air
    ! alfa0 : solar zenith angle 5 hours ago
    ! alfa5 : solar zenith angle at this time
    ! jo1d   : ratio J(O1D)/J(NO2)
    ! caco   : concentration CO in ambient air 
    !
    !
    !                   G.C.M. Vinken, KNMI, june 2010
    !=======================================================================
      USE DAO_MOD, ONLY : TS, AD, SUNCOS, SUNCOS5
      USE TRACERID_MOD, ONLY : IDO3,      IDTOX,      IDTCO
      USE TRACER_MOD, ONLY : STT, TCVV
      USE TIME_MOD, ONLY   : GET_LOCALTIME

! already called in calling program apparently.. #     include "CMN_SIZE"               ! Size parameters
#     include "cmn_fj.h"               ! Photolysis parameters
#     include "CMN_O3"                 ! fracnox, intope, jvalues

      INTEGER,    INTENT(IN)        :: I, J
      REAL,       INTENT(OUT)       :: fraction_nox,int_ope

      ! Local Variables
      INTEGER                       :: IJLOOP
      integer,parameter             :: ntemp  = 4
      integer,parameter             :: njno2  = 4
      integer,parameter             :: ncao3  = 4
      integer,parameter             :: nalfa0 = 12
      integer,parameter             :: nalfa5 = 12
      integer,parameter             :: njo1d  = 4
      integer,parameter             :: ncaco  = 4 
 
      real,dimension(ntemp)         :: templev
      real,dimension(njno2)         :: jno2lev
      real,dimension(ncao3)         :: cao3lev
      real,dimension(nalfa0)        :: alfa0lev
      real,dimension(nalfa5)        :: alfa5lev      
      real,dimension(njo1d)         :: jo1dlev
      real,dimension(ncaco)         :: cacolev

      real              :: temp_tmp,jno2_tmp,cao3_tmp                    ! Temporary variable storage
      REAL              :: alfa0_tmp,alfa5_tmp,jo1d_tmp,caco_tmp
      
      real,dimension(2) :: xtemp,xjno2,xcao3,xalfa0                     ! Interpolation parameters
      real,dimension(2) :: xalfa5,xjo1d,xcaco                           ! Interpolation parameters
      
      ! For loops
      integer           :: itemp, ijno2, icao3, ialfa0
      integer           :: ialfa5, ijo1d, icaco
      integer           :: i0,i1,i2,i3,i4,i5,i6,i7
 
      REAL,DIMENSION(7) :: var_array                                    ! array contain temp, jno2, cao3, alfa_0, alfa_5, jo1d, caco

      ! Set the levels that were chosen in the look up table
      templev = (/ 275. , 280. , 285. , 300. /)
      jno2lev = (/ 5.e-4, 0.0025, 0.0050, 0.012 /)
      cao3lev = (/ 5., 20., 35., 75. /)
      alfa0lev = (/ -90., -60., -45., -30., -15., 0., 15., 30.,
     $              45., 60., 75., 90. /)
      alfa5lev = (/ -90., -60., -45., -30., -15., 0., 15., 30.,
     $              45., 60., 75., 90. /)
      jo1dlev = (/ 5.e-4, 0.0015, 0.0025, 0.0055 /)
      cacolev = (/ 50., 100., 150., 1200. /) 

c      print*,"Temperature levels are: ",templev
c      print*,"This is grid cell: ",I,J
      
      ! Temperature
c      print*,"Temperature here is: ",TS(I,J)
c      print*,"USA: ",TS(32,64)

      ! Tracer concentrations in v/v
c      print*,"[O3] is: ",STT(I,J,1,IDTOX)/ AD(I,J,1) * TCVV(IDTOX)
c      print*,"[CO] is: ",STT(I,J,1,IDTCO)/ AD(I,J,1) * TCVV(IDTCO)
c      print*,"IDTOX is: ", IDTOX
c      print*,"IDO3 is: ", IDO3
c      print*,"In USA: ",STT(32,64,1,IDTOX)/ AD(32,64,1) * TCVV(IDTOX)
      
      ! SOLAR ZENITH ANGLES IN DEGREES
c      IJLOOP = ( (J-1) * IIPAR ) + I
c      print*,"Local Time: ",GET_LOCALTIME(I)
c      print*,"Solar Zenith Angle at this location: ", 
c     $            ASIND(SUNCOS(IJLOOP))
c      IJLOOP = ( (64-1) * IIPAR ) + 32
c      print*,"Local USA time: ", GET_LOCALTIME(32)
c      print*,"Solar Zenith Angle at USA: ", ASIND(SUNCOS(IJLOOP))
c      print*,"Solar Zenith Angle at USA - 5: ",ASIND(SUNCOS5(IJLOOP))
      
      ! Set the variables
      IJLOOP = ( (J-1) * IIPAR ) + I
      var_array(1) = TS(I,J)                                            ! Temperature
      var_array(2) = jvalues(I,J,1)                                     ! J(NO2)
      var_array(3) = STT(I,J,1,IDTOX)/ AD(I,J,1) * TCVV(IDTOX) * 1.E9   ! [O3] in ppbv
      var_array(4) = ASIND(SUNCOS5(IJLOOP))                             ! alfa0
      var_array(5) = ASIND(SUNCOS(IJLOOP))                              ! alfa5
      var_array(6) = jvalues(I,J,2) / jvalues(I,J,1)                    ! J(O1D)/J(NO2)
      var_array(7) = STT(I,J,1,IDTCO)/ AD(I,J,1) * TCVV(IDTCO) * 1.E9   ! [CO] in ppbv
      
      if (jvalues(I,J,1) .eq. 0.) var_array(6) = 0.                     ! prevent NaN when jvalues are 0.
      
      ! First some error checking
c     ########### MAYBE CHECK HERE FOR NEGATIVE VALUES?##########
 
      !
      ! Determine reference index (itemp,ijno2,icao3,ialfa0,ialfa5,ialfajo1d,icaco)
      !
      !===========================================================================
      ! Find smallest temperature reference level (i) for which actual temperature
      ! is smaller, then do
      !
      ! x(1) = (  temperature_level(i+1) - actual temperature  )
      !        ------------------------------------------------- 
      !        ( temperature_level(i+1) - temperature_level(i) )
      !
      ! then x(2) = 1.0 - x(1)
      !
      !===========================================================================

      ! Temperature:
      temp_tmp = var_array(1)
      if ( var_array(1) > templev(ntemp) ) temp_tmp = templev(ntemp)        ! If temperature larger than largest in LUT, assign largest temp
      if ( var_array(1) < templev(1) ) temp_tmp = templev(1)        ! If temp smaller, assign smallest temp level
      do i0=1,ntemp-1
         itemp = i0
         if( templev(itemp+1) > temp_tmp ) exit 
      end do
      xtemp(1)=(templev(itemp+1)-temp_tmp)/
     $                        (templev(itemp+1)-templev(itemp))
      xtemp(2)=1.0-xtemp(1)       

      ! J(NO2):
      jno2_tmp = var_array(2)
      if ( var_array(2) > jno2lev(njno2) ) jno2_tmp = jno2lev(njno2)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(2) < jno2lev(1) ) jno2_tmp = jno2lev(1)        ! If smaller, assign smallest level value
      do i0=1,njno2-1
         ijno2 = i0
         if( jno2lev(ijno2+1) > jno2_tmp ) exit 
      end do
      xjno2(1)=(jno2lev(ijno2+1)-jno2_tmp)/
     $                        (jno2lev(ijno2+1)-jno2lev(ijno2))
      xjno2(2)=1.0-xjno2(1)       

      ! [O3]:
      cao3_tmp = var_array(3)
      if ( var_array(3) > cao3lev(ncao3) ) cao3_tmp = cao3lev(ncao3)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(3) < cao3lev(1) ) cao3_tmp = cao3lev(1)                ! If smaller, assign smallest level value
      do i0=1,ncao3-1
         icao3 = i0
         if( cao3lev(icao3+1) > cao3_tmp ) exit 
      end do
      xcao3(1)=(cao3lev(icao3+1)-cao3_tmp)/
     $                        (cao3lev(icao3+1)-cao3lev(icao3))
      xcao3(2)=1.0-xcao3(1)       

      ! alfa0:
      alfa0_tmp = var_array(4)
      if ( var_array(4) > alfa0lev(nalfa0) ) alfa0_tmp = 
     $                                          alfa0lev(nalfa0)   ! If larger than largest in LUT, assign largest level values
      if ( var_array(4) < alfa0lev(1) ) alfa0_tmp = alfa0lev(1)             ! If smaller, assign smallest level value
      do i0=1,nalfa0-1
         ialfa0 = i0
         if( alfa0lev(ialfa0+1) > alfa0_tmp ) exit 
      end do
      xalfa0(1)=(alfa0lev(ialfa0+1)-alfa0_tmp)/
     $                        (alfa0lev(ialfa0+1)-alfa0lev(ialfa0))
      xalfa0(2)=1.0-xalfa0(1)       

      ! alfa5:
      alfa5_tmp = var_array(5)
      if ( var_array(5) > alfa5lev(nalfa5) ) alfa5_tmp = 
     $                                          alfa5lev(nalfa5)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(5) < alfa5lev(1) ) alfa5_tmp = alfa5lev(1)                ! If smaller, assign smallest level value
      do i0=1,nalfa5-1
         ialfa5 = i0
         if( alfa5lev(ialfa5+1) > alfa5_tmp ) exit 
      end do
      xalfa5(1)=(alfa5lev(ialfa5+1)-alfa5_tmp)/
     $                        (alfa5lev(ialfa5+1)-alfa5lev(ialfa5))
      xalfa5(2)=1.0-xalfa5(1)       

      ! jo1d:
      jo1d_tmp = var_array(6)
      if ( var_array(6) > jo1dlev(njo1d) ) jo1d_tmp = jo1dlev(njo1d)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(6) < jo1dlev(1) ) jo1d_tmp = jo1dlev(1)                ! If smaller, assign smallest level value
      do i0=1,njo1d-1
         ijo1d = i0
         if( jo1dlev(ijo1d+1) > jo1d_tmp ) exit 
      end do
      xjo1d(1)=(jo1dlev(ijo1d+1)-jo1d_tmp)/
     $                        (jo1dlev(ijo1d+1)-jo1dlev(ijo1d))
      xjo1d(2)=1.0-xjo1d(1)       

      ! [CO]:
      caco_tmp = var_array(7)
      if ( var_array(7) > cacolev(ncaco) ) caco_tmp = cacolev(ncaco)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(7) < cacolev(1) ) caco_tmp = cacolev(1)                ! If smaller, assign smallest level value
      do i0=1,ncaco-1
         icaco = i0
         if( cacolev(icaco+1) > caco_tmp ) exit 
      end do
      xcaco(1)=(cacolev(icaco+1)-caco_tmp)/
     $                        (cacolev(icaco+1)-cacolev(icaco))
      xcaco(2)=1.0-xcaco(1)       

c      print*,"The i-values are:", itemp, ijno2, icao3, ialfa0,
c     $                            ialfa5, ijo1d, icaco
c      print*,"Variables are: ", var_array
c      print*,"For testing, xtemp: ",xtemp
      
      !
      ! Linear interpolation
      !
      fraction_nox=0.0
      int_ope = 0.0
      do i1=1,2
       do i2=1,2
        do i3=1,2
         do i4=1,2
          do i5=1,2
           do i6=1,2
            do i7=1,2 
              IF (fracnox(itemp+i1-1,ijno2+i2-1,icao3+i3-1,ialfa0+i4-1,     !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!       
     $        ialfa5+i5-1,ijo1d+i6-1,icaco+i7-1) .LT. 0.) print*,"##E#"
              IF ((ialfa0+i4-1 .eq. 6) .and. (ialfa5+i5-1 .eq. 6) ) 
     $                        CYCLE   !Cycle if both angles are 0   
              fraction_nox=fraction_nox+xtemp(i1)*xjno2(i2)*
     $           xcao3(i3)*xalfa0(i4)*xalfa5(i5)*xjo1d(i6)*xcaco(i7)* 
     $                  fracnox(itemp+i1-1,ijno2+i2-1,icao3+i3-1,          ! fracnox is the array with the actual lut data
     $                  ialfa0+i4-1,ialfa5+i5-1,ijo1d+i6-1,icaco+i7-1)
              int_ope=int_ope+xtemp(i1)*xjno2(i2)*xcao3(i3)* 
     $                  xalfa0(i4)*xalfa5(i5)*xjo1d(i6)*xcaco(i7)* 
     $                  intope(itemp+i1-1,ijno2+i2-1,icao3+i3-1,           ! intope is the array with the actual lut data
     $                  ialfa0+i4-1,ialfa5+i5-1,ijo1d+i6-1,icaco+i7-1)
            end do
           end do
          end do
         end do
        end do
       end do
      end do
      
      if ((I .eq. 108) .and. (J .eq. 49)) then
           print*,"----INTERPOLATE_LUT-----"
           print*,"fraction_nox and int_OPE: ",fraction_nox,
     &                                            int_ope
           print*,"Jvalues are: ",jvalues(I, J,:)
           print*,"Vars are: ",var_array
           print*,"[O3] in interpolate_lut: ",var_array(3)
           print*,"[CO] in interpolate_lut: ",var_array(7)
           print*,"The i-values are:", itemp, ijno2, icao3, ialfa0,
     $                            ialfa5, ijo1d, icaco
       endif
               
c      print*,"fraction_nox is: ",fraction_nox
c      print*,"integrated OPE: ",int_ope
      
      end subroutine INTERPOLATE_LUT

!EOC
!------------------------------------------------------------------------------

      subroutine INTERPOLATE_LUT2(I,J,o3,no,no2,dens,fraction_nox,
     &                            int_ope) 
    !=======================================================================
    !
    ! INTERPOLATE_LUT:    Return FracNOx or IntOPE from the lookup table
    !
    ! temp   : model temperature
    ! jno2   : J(NO2) value
    ! cao3   : concentration O3 in ambient air
    ! alfa0 : solar zenith angle 5 hours ago
    ! alfa5 : solar zenith angle at this time
    ! jo1d   : ratio J(O1D)/J(NO2)
    ! canox  : concentration NOx in ambient air 
    ! 
    ! o3     : incoming o3 concentration
    ! no     : incoming no
    ! no2    : incoming no2
    ! dens   : incoming air density
    !
    !
    !                   G.C.M. Vinken, KNMI, june 2010
    ! 02-21-2011: Updated for NOx in LUT
    !=======================================================================
      USE DAO_MOD, ONLY : TS, SUNCOS, SUNCOS5
      USE TIME_MOD, ONLY   : GET_LOCALTIME

! already called in calling program apparently.. #     include "CMN_SIZE"               ! Size parameters
#     include "cmn_fj.h"               ! Photolysis parameters
#     include "CMN_O3"                 ! fracnox, intope, jvalues

      INTEGER,    INTENT(IN)        :: I, J
      REAL*8,       INTENT(IN)        :: o3,no,no2,dens
      REAL,       INTENT(OUT)       ::  fraction_nox,int_ope

      ! Local Variables
      INTEGER                       :: IJLOOP
      integer,parameter             :: ntemp  = 4
      integer,parameter             :: njno2  = 4
      integer,parameter             :: ncao3  = 4
      integer,parameter             :: nalfa0 = 12
      integer,parameter             :: nalfa5 = 12
      integer,parameter             :: njo1d  = 4
      integer,parameter             :: ncanox = 5 
 
      real,dimension(ntemp)         :: templev
      real,dimension(njno2)         :: jno2lev
      real,dimension(ncao3)         :: cao3lev
      real,dimension(nalfa0)        :: alfa0lev
      real,dimension(nalfa5)        :: alfa5lev      
      real,dimension(njo1d)         :: jo1dlev
      real,dimension(ncanox)         :: canoxlev

      real              :: temp_tmp,jno2_tmp,cao3_tmp                    ! Temporary variable storage
      REAL              :: alfa0_tmp,alfa5_tmp,jo1d_tmp,canox_tmp
      
      real,dimension(2) :: xtemp,xjno2,xcao3,xalfa0                     ! Interpolation parameters
      real,dimension(2) :: xalfa5,xjo1d,xcanox                           ! Interpolation parameters
      
      ! For loops
      integer           :: itemp, ijno2, icao3, ialfa0
      integer           :: ialfa5, ijo1d, icanox
      integer           :: i0,i1,i2,i3,i4,i5,i6,i7
 
      REAL,DIMENSION(7) :: var_array                                    ! array contain temp, jno2, cao3, alfa_0, alfa_5, jo1d, canox

      ! Set the levels that were chosen in the look up table
      templev = (/ 275. , 280. , 285. , 310. /)
      jno2lev = (/ 5.e-4, 0.0025, 0.0050, 0.012 /)
      cao3lev = (/ 5., 20., 35., 75. /)
      alfa0lev = (/ -90., -60., -45., -30., -15., 0., 15., 30.,
     $              45., 60., 75., 90. /)
      alfa5lev = (/ -90., -60., -45., -30., -15., 0., 15., 30.,
     $              45., 60., 75., 90. /)
      jo1dlev = (/ 5.e-4, 0.0015, 0.0025, 0.0055 /)
      canoxlev = (/ 10., 200., 1000., 2000.,6000. /) 

c      print*,"Temperature levels are: ",templev
c      print*,"This is grid cell: ",I,J
      
      ! Temperature
c      print*,"Temperature here is: ",TS(I,J)
c      print*,"USA: ",TS(32,64)

      ! Tracer concentrations in v/v
c      print*,"[O3] is: ",STT(I,J,1,IDTOX)/ AD(I,J,1) * TCVV(IDTOX)
c      print*,"[CO] is: ",STT(I,J,1,IDTCO)/ AD(I,J,1) * TCVV(IDTCO)
c      print*,"IDTOX is: ", IDTOX
c      print*,"IDO3 is: ", IDO3
c      print*,"In USA: ",STT(32,64,1,IDTOX)/ AD(32,64,1) * TCVV(IDTOX)
      
      ! SOLAR ZENITH ANGLES IN DEGREES
c      IJLOOP = ( (J-1) * IIPAR ) + I
c      print*,"Local Time: ",GET_LOCALTIME(I)
c      print*,"Solar Zenith Angle at this location: ", 
c     $            ASIND(SUNCOS(IJLOOP))
c      IJLOOP = ( (64-1) * IIPAR ) + 32
c      print*,"Local USA time: ", GET_LOCALTIME(32)
c      print*,"Solar Zenith Angle at USA: ", ASIND(SUNCOS(IJLOOP))
c      print*,"Solar Zenith Angle at USA - 5: ",ASIND(SUNCOS5(IJLOOP))
      
      ! Set the variables
      IJLOOP = ( (J-1) * IIPAR ) + I
      var_array(1) = TS(I,J)                                            ! Temperature
      var_array(2) = jvalues(I,J,1)                                     ! J(NO2)
      var_array(3) = o3 / dens * 1.E9                                   ! [O3] in ppbv
      var_array(4) = ASIND(SUNCOS5(IJLOOP))                             ! alfa0
      var_array(5) = ASIND(SUNCOS(IJLOOP))                              ! alfa5
      var_array(6) = jvalues(I,J,2) / jvalues(I,J,1)                    ! J(O1D)/J(NO2)
      var_array(7) = (no + no2) / dens * 1.E12                          ! [NOx] in pptv
      
      if (jvalues(I,J,1) .eq. 0.) var_array(6) = 0.                     ! prevent NaN when jvalues are 0.

      ! First some error checking
c     ########### MAYBE CHECK HERE FOR NEGATIVE VALUES?##########
 
      !
      ! Determine reference index (itemp,ijno2,icao3,ialfa0,ialfa5,ialfajo1d,icaco)
      !
      !===========================================================================
      ! Find smallest temperature reference level (i) for which actual temperature
      ! is smaller, then do
      !
      ! x(1) = (  temperature_level(i+1) - actual temperature  )
      !        ------------------------------------------------- 
      !        ( temperature_level(i+1) - temperature_level(i) )
      !
      ! then x(2) = 1.0 - x(1)
      !
      !===========================================================================

      ! Temperature:
      temp_tmp = var_array(1)
      if ( var_array(1) > templev(ntemp) ) temp_tmp = templev(ntemp)        ! If temperature larger than largest in LUT, assign largest temp
      if ( var_array(1) < templev(1) ) temp_tmp = templev(1)        ! If temp smaller, assign smallest temp level
      do i0=1,ntemp-1
         itemp = i0
         if( templev(itemp+1) > temp_tmp ) exit 
      end do
      xtemp(1)=(templev(itemp+1)-temp_tmp)/
     $                        (templev(itemp+1)-templev(itemp))
      xtemp(2)=1.0-xtemp(1)       

      ! J(NO2):
      jno2_tmp = var_array(2)
      if ( var_array(2) > jno2lev(njno2) ) jno2_tmp = jno2lev(njno2)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(2) < jno2lev(1) ) jno2_tmp = jno2lev(1)        ! If smaller, assign smallest level value
      do i0=1,njno2-1
         ijno2 = i0
         if( jno2lev(ijno2+1) > jno2_tmp ) exit 
      end do
      xjno2(1)=(jno2lev(ijno2+1)-jno2_tmp)/
     $                        (jno2lev(ijno2+1)-jno2lev(ijno2))
      xjno2(2)=1.0-xjno2(1)       

      ! [O3]:
      cao3_tmp = var_array(3)
      if ( var_array(3) > cao3lev(ncao3) ) cao3_tmp = cao3lev(ncao3)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(3) < cao3lev(1) ) cao3_tmp = cao3lev(1)                ! If smaller, assign smallest level value
      do i0=1,ncao3-1
         icao3 = i0
         if( cao3lev(icao3+1) > cao3_tmp ) exit 
      end do
      xcao3(1)=(cao3lev(icao3+1)-cao3_tmp)/
     $                        (cao3lev(icao3+1)-cao3lev(icao3))
      xcao3(2)=1.0-xcao3(1)       

      ! alfa0:
      alfa0_tmp = var_array(4)
      if ( var_array(4) > alfa0lev(nalfa0) ) alfa0_tmp = 
     $                                          alfa0lev(nalfa0)   ! If larger than largest in LUT, assign largest level values
      if ( var_array(4) < alfa0lev(1) ) alfa0_tmp = alfa0lev(1)             ! If smaller, assign smallest level value
      do i0=1,nalfa0-1
         ialfa0 = i0
         if( alfa0lev(ialfa0+1) > alfa0_tmp ) exit 
      end do
      xalfa0(1)=(alfa0lev(ialfa0+1)-alfa0_tmp)/
     $                        (alfa0lev(ialfa0+1)-alfa0lev(ialfa0))
      xalfa0(2)=1.0-xalfa0(1)       

      ! alfa5:
      alfa5_tmp = var_array(5)
      if ( var_array(5) > alfa5lev(nalfa5) ) alfa5_tmp = 
     $                                          alfa5lev(nalfa5)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(5) < alfa5lev(1) ) alfa5_tmp = alfa5lev(1)                ! If smaller, assign smallest level value
      do i0=1,nalfa5-1
         ialfa5 = i0
         if( alfa5lev(ialfa5+1) > alfa5_tmp ) exit 
      end do
      xalfa5(1)=(alfa5lev(ialfa5+1)-alfa5_tmp)/
     $                        (alfa5lev(ialfa5+1)-alfa5lev(ialfa5))
      xalfa5(2)=1.0-xalfa5(1)       

      ! jo1d:
      jo1d_tmp = var_array(6)
      if ( var_array(6) > jo1dlev(njo1d) ) jo1d_tmp = jo1dlev(njo1d)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(6) < jo1dlev(1) ) jo1d_tmp = jo1dlev(1)                ! If smaller, assign smallest level value
      do i0=1,njo1d-1
         ijo1d = i0
         if( jo1dlev(ijo1d+1) > jo1d_tmp ) exit 
      end do
      xjo1d(1)=(jo1dlev(ijo1d+1)-jo1d_tmp)/
     $                        (jo1dlev(ijo1d+1)-jo1dlev(ijo1d))
      xjo1d(2)=1.0-xjo1d(1)       

      ! [NOx]:
      canox_tmp = var_array(7)
      if ( var_array(7) > canoxlev(ncanox) ) canox_tmp = 
     $                                         canoxlev(ncanox)        ! If larger than largest in LUT, assign largest level values
      if ( var_array(7) < canoxlev(1) ) canox_tmp = canoxlev(1)                ! If smaller, assign smallest level value
      do i0=1,ncanox-1
         icanox = i0
         if( canoxlev(icanox+1) > canox_tmp ) exit 
      end do
      xcanox(1)=(canoxlev(icanox+1)-canox_tmp)/
     $                        (canoxlev(icanox+1)-canoxlev(icanox))
      xcanox(2)=1.0-xcanox(1)       

c      print*,"The i-values are:", itemp, ijno2, icao3, ialfa0,
c     $                            ialfa5, ijo1d, icanox
c      print*,"Variables are: ", var_array
c      print*,"For testing, xtemp: ",xtemp
      
      !
      ! Linear interpolation
      !
      fraction_nox=0.0
      int_ope = 0.0
      do i1=1,2
       do i2=1,2
        do i3=1,2
         do i4=1,2
          do i5=1,2
           do i6=1,2
            do i7=1,2 
              IF (fracnox(itemp+i1-1,ijno2+i2-1,icao3+i3-1,ialfa0+i4-1,     !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!       
     $      ialfa5+i5-1,ijo1d+i6-1,icanox+i7-1) .LT. 0.) print*,"##E#"
              fraction_nox=fraction_nox+xtemp(i1)*xjno2(i2)
     $             *xcao3(i3)*xalfa0(i4)*xalfa5(i5)*xjo1d(i6)*
     $             xcanox(i7)*fracnox(itemp+i1-1,ijno2+i2-1,icao3+i3-1,          ! fracnox is the array with the actual lut data
     $             ialfa0+i4-1,ialfa5+i5-1,ijo1d+i6-1,icanox+i7-1)
              int_ope=int_ope+xtemp(i1)*xjno2(i2)*xcao3(i3)* 
     $                  xalfa0(i4)*xalfa5(i5)*xjo1d(i6)*xcanox(i7)* 
     $                  intope(itemp+i1-1,ijno2+i2-1,icao3+i3-1,           ! intope is the array with the actual lut data
     $                  ialfa0+i4-1,ialfa5+i5-1,ijo1d+i6-1,icanox+i7-1)
            end do
           end do
          end do
         end do
        end do
       end do
      end do
      !
      
      if ((I .eq. 108) .and. (J .eq. 49)) then
           print*,"-----INTERPOLATE_LUT2, for 108,49-----"
           print*,"Fraction_nox and int_OPE: ",fraction_nox,
     &                                            int_ope
           print*,"Jvalues are: ",jvalues(I, J,:)
           print*,"Vars are: ",var_array
           print*,"[O3] in interpolate_lut: ",var_array(3)
           print*,"[NOx] in interpolate_lut: ",var_array(7)
           print*,"J(O1D)/J(NO2) : ",var_array(6)
           print*,"The i-values are:", itemp, ijno2, icao3, ialfa0,
     $                            ialfa5, ijo1d, icanox
           print*,"Interpolation parameters: ",xtemp,xjno2,xcao3,
     $                  xalfa0,xalfa5,xjo1d,xcanox 
           print*,"---------------------------------"
      endif
       
      if ((I .eq. 73) .and. (J .eq. 76)) then
           print*,"-----INTERPOLATE_LUT2, for 73,76-----"
           print*,"Fraction_nox and int_OPE: ",fraction_nox,
     &                                            int_ope
           print*,"Jvalues are: ",jvalues(I, J,:)
           print*,"Vars are: ",var_array
           print*,"[O3] in interpolate_lut: ",var_array(3)
           print*,"[NOx] in interpolate_lut: ",var_array(7)
           print*,"J(O1D)/J(NO2) : ",var_array(6)
           print*,"The i-values are:", itemp, ijno2, icao3, ialfa0,
     $                            ialfa5, ijo1d, icanox
           print*,"Interpolation parameters: ",xtemp,xjno2,xcao3,
     $                  xalfa0,xalfa5,xjo1d,xcanox 
           print*,"------------------------------------"
      endif
      
      end subroutine INTERPOLATE_LUT2

!EOC
!------------------------------------------------------------------------------

      ! End of module
      END MODULE ICOADS_SHIP_MOD
