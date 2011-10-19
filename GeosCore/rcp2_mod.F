
      MODULE RCP2_MOD


      IMPLICIT NONE

      !Make everything PRIVATE...
      PRIVATE

      !...except these routines
      PUBLIC :: CLEANUP_RCP
      PUBLIC :: EMISS_RCP
      PUBLIC :: GET_RCP_EMISSION
      PUBLIC :: RCPNAME, RCPYEAR
      PUBLIC :: RCP_AIREMISS

      ! Anthropogenic land surface emissions.
      ! includes fossil fuel and biofuel emissions
      ! from energy production and distribution, 
      ! residential and commercial combustion, 
      ! industry, transportation, waste treatment and disposal
      ! solvent production and use, 
      ! agriculture, and agricultural waste burning
      REAL*4, ALLOCATABLE :: RCP_NO(:,:)
      REAL*4, ALLOCATABLE :: RCP_CO(:,:)
      REAL*4, ALLOCATABLE :: RCP_OC(:,:)
      REAL*4, ALLOCATABLE :: RCP_BC(:,:)
      REAL*4, ALLOCATABLE :: RCP_SO2(:,:)
      REAL*4, ALLOCATABLE :: RCP_NH3(:,:)
      REAL*4, ALLOCATABLE :: RCP_ALK4(:,:)
      REAL*4, ALLOCATABLE :: RCP_ACET(:,:)
      REAL*4, ALLOCATABLE :: RCP_MEK(:,:)
      REAL*4, ALLOCATABLE :: RCP_ALD2(:,:)
      REAL*4, ALLOCATABLE :: RCP_PRPE(:,:)
      REAL*4, ALLOCATABLE :: RCP_C3H8(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H6(:,:)
      REAL*4, ALLOCATABLE :: RCP_CH2O(:,:)
      REAL*4, ALLOCATABLE :: RCP_BENZ(:,:)
      REAL*4, ALLOCATABLE :: RCP_TOLU(:,:)
      REAL*4, ALLOCATABLE :: RCP_XYLE(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H4(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H2(:,:)
      REAL*4, ALLOCATABLE :: RCP_HCOOH(:,:)

      ! Aircraft emissions
      REAL*4, ALLOCATABLE :: RCP_NO_AIR(:,:,:)
      REAL*4, ALLOCATABLE :: RCP_BC_AIR(:,:,:)

      ! Ship emissions
      REAL*4, ALLOCATABLE :: RCP_NO_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_CO_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_OC_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_BC_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_SO2_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_NH3_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_ALK4_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_ACET_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_MEK_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_ALD2_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_PRPE_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_C3H8_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H6_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_CH2O_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_BENZ_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_TOLU_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_XYLE_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H4_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_C2H2_SHIP(:,:)
      REAL*4, ALLOCATABLE :: RCP_HCOOH_SHIP(:,:)
      
      REAL*4, ALLOCATABLE :: RCP_LAND(:,:,:)
      REAL*4, ALLOCATABLE :: RCP_SHIP(:,:,:)

      CHARACTER(LEN=20)   :: RCPNAME, RCPYEAR
      INTEGER             :: IDTRCP(20), IDTRCP_SHIP(20)

      !========================================
      !  MODULE ROUTINES -- follow "CONTAINS"
      !========================================
      CONTAINS

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      SUBROUTINE EMISS_RCP

!***********************************************************************
!  Subroutine EMISS_RCP reads all RCP emissions at the beginning of
!  each month. (cdh, 10/14/11)
!***********************************************************************

      USE TIME_MOD,               ONLY : GET_MONTH
      USE ERROR_MOD,              ONLY : GEOS_CHEM_STOP
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE LOGICAL_MOD,            ONLY : LRCPAIR
      USE BPCH2_MOD,              ONLY : GET_TAU0, GET_RES_EXT
      USE TRACERID_MOD
      USE TRACER_MOD, ONLY      : TRACER_NAME
      USE COMODE_LOOP_MOD, ONLY : IDEMS

      !Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: THISMONTH
      CHARACTER(LEN=20)             :: RCPSPECIES
      INTEGER :: I
      CHARACTER(LEN=255)            :: FILENAME
      REAL*8                        :: XTAU


      !=================================================================
      !  EMISS_RCP begins here
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays
         CALL INIT_RCP

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      ! Get month
      THISMONTH = GET_MONTH()      

      ! Echo info
      WRITE(6, '(a)' ) REPEAT( '=', 79)
      WRITE(6, 100   ) RCPNAME, RCPYEAR
 100  FORMAT( 'R C P    E M I S S I O N S',
     &        '  --  Scenario: ', 2A10 , / )

      ! Anthro file name
      FILENAME = TRIM( DATA_DIR )        // 'RCP_201110/'     //
     &                  trim( RCPNAME )  // '_anthropogenic_' // 
     &                  trim( RCPYEAR )  // '.'               //
     &                  GET_RES_EXT()    // '.bpch' 

      XTAU = GET_TAU0( 1, 1, 2000 ) 
      ! Read data
      CALL READ_RCP( FILENAME, SHIP=.FALSE. )

      ! Anthro file name
      FILENAME = TRIM( DATA_DIR )        // 'RCP_201110/'     //
     &                  trim( RCPNAME )  // '_ships_'         // 
     &                  trim( RCPYEAR )  // '.'               //
     &                  GET_RES_EXT()    // '.bpch' 

      XTAU = GET_TAU0( THISMONTH, 1, 2000 )  

      ! Read data
      CALL READ_RCP( FILENAME, SHIP=.TRUE., TAU0=XTAU )

        
      !=================================================================
      ! Aircraft emissions
      !=================================================================

      IF (LRCPAIR) THEN
         FILENAME = TRIM( DATA_DIR ) // 'RCP_201110/'     //
     &              TRIM( RCPNAME )  // '_NO_aircraft_'   // 
     &              TRIM( RCPYEAR )  // '.nc'

         CALL READ_RCP_AIRCRAFT( FILENAME, THISMONTH )

      ENDIF

      !=================================================================
      ! Print totals to log
      !=================================================================

      CALL TOTAL_ANTHRO_RCP( THISMONTH )

      ! Fancy output
      WRITE(6, '(a)' ) REPEAT( '=', 79)

!      PRINT*,'ALL OK'
!      CALL GEOS_CHEM_STOP

      ! Return to calling program
      END SUBROUTINE EMISS_RCP

!----------------------------------------------------------------------

      SUBROUTINE READ_RCP( FILENAME, SHIP, TAU0 )

!***********************************************************************
!  Subroutine READ_RCP reads a BPCH file created from RCP data.
!***********************************************************************

      ! Reference to F90 modules
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      
      USE TIME_MOD,         ONLY : GET_YEAR
      USE ERROR_MOD,        ONLY : ERROR_STOP

      USE COMODE_LOOP_MOD,      ONLY : IDEMS
      USE TRACER_MOD, ONLY      : TRACER_NAME

      USE CMN_SIZE_MOD           ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      LOGICAL,          INTENT(IN) :: SHIP
      REAL*8,OPTIONAL,  INTENT(IN) :: TAU0


      ! Local variables
      INTEGER                         :: I, J, L, N, IOS, K
      INTEGER                         :: NI, NJ, NL
      INTEGER                         :: IFIRST, JFIRST, LFIRST
      INTEGER                         :: NTRACER, NSKIP
      INTEGER                         :: HALFPOLAR, CENTER180
      INTEGER                         :: SCALEYEAR, BASEYEAR
      REAL*4                          :: LONRES, LATRES
      REAL*4                          :: ARRAY(IGLOB,JGLOB,1)
      REAL*4                          :: SC(IIPAR,JJPAR)
      REAL*8                          :: ZTAU0, ZTAU1
      CHARACTER(LEN=20)               :: MODELNAME
      CHARACTER(LEN=40)               :: CATEGORY
      CHARACTER(LEN=40)               :: UNIT
      CHARACTER(LEN=40)               :: RESERVED

      !=================================================================
      !  READ_RCP begins here
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RCP: Reading ', a )

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME)

      ! Initialize
      K = 0

      ! Read the entire file in one pass
      DO

         ! Read 1st data block header
         READ( IU_FILE, IOSTAT=IOS )
     &   MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! Check for EOF or errors
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:2' )

         ! Read 2nd data block header line
         READ (IU_FILE, IOSTAT=IOS )
     &   CATEGORY, NTRACER, UNIT, ZTAU0, ZTAU1, RESERVED,
     &   NI, NJ, NL, IFIRST, JFIRST, LFIRST, NSKIP

         IF ( CATEGORY /= 'ANTHSRCE' ) 
     &        CALL ERROR_STOP( 'ANTHSRCE not found', 'READ_RCP' ) 
         
         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:3' )

         ! Read data
         READ( IU_FILE, IOSTAT=IOS )
     &        ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'read_data:4' )

         !==============================================================
         !  Save into tracer arrays
         !==============================================================

         ! Select date, if this argument is present
         IF ( PRESENT( TAU0 ) ) THEN
            IF (ZTAU0 /= TAU0) CYCLE
         ENDIF
            
         ! Tracer number must be positive
         IF (NTRACER > 0 .AND. NTRACER <= SIZE(TRACER_NAME) ) THEN

            ! Increment tracer counter
            K = K + 1 

            ! Save emissions and tracer number
            IF (SHIP) THEN
               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_SHIP(:,:,K) )
               IDTRCP_SHIP(K) = NTRACER 
            ELSE
               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_LAND(:,:,K) )
               IDTRCP(K) = NTRACER 
            ENDIF
               
         ENDIF

      END DO

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE READ_RCP

!-----------------------------------------------------------------------

      SUBROUTINE READ_RCP_AIRCRAFT( FILENAME, MONTH )

!***********************************************************************
!  Subroutine READ_RCP_AIRCRAFT reads a netCDF file with RCP data.
!***********************************************************************

      ! Reference to F90 modules
      USE GRID_MOD,               ONLY : GET_XEDGE, GET_YEDGE
      USE RCP_MOD,                ONLY : READ_AND_REGRID_VOLUME
      USE PRESSURE_MOD,           ONLY : GET_AP, GET_BP
      USE DAO_MOD,                ONLY : AIRVOL
      USE CMN_SIZE_MOD               ! Size parameters

      ! Arguments
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      INTEGER,         INTENT(IN) :: MONTH

      !Local variables
      LOGICAL, SAVE                         :: FIRST = .TRUE.
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: tmpEmission

      REAL*8, SAVE :: XEDGE(IIPAR+1)
      REAL*8, SAVE :: YEDGE(JJPAR+1)
      REAL*8, SAVE :: PEDGE(LLPAR+1)
      INTEGER :: I
  
      !=================================================================
      !  READ_RCP_AIRCRAFT begins here
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Create array of longitude edges
         DO I=1, IIPAR
            XEDGE(I) = get_XEDGE( I )
         ENDDO
         XEDGE(IIPAR+1) = XEDGE(1)

         ! Create array of latitude edges
         DO I=1, JJPAR
            yedge(I) = get_yedge( I )
         ENDDO
         yedge(JJPAR+1) = 90d0
         
         ! Create mean pressure, 
         ! Assume 1013 hPa surface, which is consistent with netCDF files
         DO I=1, LLPAR+1
            pedge(I) = GET_AP(I) + GET_BP(I) * 1013.25
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RCP_AIRCRAFT: Reading ', a )


      ! Read aircraft NO, molec/box/s
      CALL READ_AND_REGRID_VOLUME( 
     &        FILENAME=FILENAME, 
     &        varname1='emiss_air',
     &        outLonEdge=XEDGE,
     &        outLatEdge=YEDGE,
     &        outLevEdge=PEDGE, 
     &        nMoments=1,
     &        outEmission=tmpEmission )

      ! Select only the desired month, molec/box/s
      RCP_NO_AIR = tmpEmission(:,:,:,1,MONTH) 

      ! Convert molec/box/s -> molec/cm3/s
      RCP_NO_AIR = RCP_NO_AIR / ( AIRVOL * 1D6 )

      ! Deallocate temporary array
      IF ( ALLOCATED(tmpEmission) ) DEALLOCATE( tmpEmission )

      ! Return to calling program
      END SUBROUTINE READ_RCP_AIRCRAFT

!-----------------------------------------------------------------------

      SUBROUTINE RCP_AIREMISS

!***********************************************************************
!  Populate EMIS_AC_NOx with aircraft NOx emissions. Also do diagnostics
!***********************************************************************

      ! References to F90 modules
      USE DAO_MOD,                ONLY : BXHEIGHT
      USE AIRCRAFT_NOX_MOD,       ONLY : EMIS_AC_NOx
      USE DIAG_MOD,               ONLY : AD32_AC
      USE CMN_DIAG_MOD  ! Diagnostic switches

      !=================================================================
      !  RCP_AIREMISS begins here
      !=================================================================

      ! Convert  -> molec/cm3/s
      EMIS_AC_NOx = RCP_NO_AIR

      ! ND32 -- save NOx in [molec/cm2], will convert to
      ! [molec/cm2/s] in subroutine "diag3.f" (bmy, 3/16/00)
      IF ( ND32 > 0 ) THEN
         AD32_ac = AD32_ac + ( EMIS_AC_NOx * BXHEIGHT * 1d2 )
      ENDIF

      ! Return to calling program
      END SUBROUTINE RCP_AIREMISS

!-----------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_RCP( THISMONTH )

!***********************************************************************
!  Subroutine to print total RCP anthropogenic emissions each month
!  in TgC
!***********************************************************************

      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACER_MOD, ONLY      : TRACER_NAME
      USE ERROR_MOD, ONLY :GEOS_CHEM_STOP

      USE CMN_SIZE_MOD       ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J, K
      REAL*8              :: A, TOTAL, TOTAL_SHIP

      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      !  TOTAL_ANTHRO_RCP begins here
      !=================================================================

      !==============================================================
      ! RCP Land emissions
      !==============================================================

      WRITE( 6, '(a)' )
      DO K=1, SIZE(IDTRCP)
         
         IF (IDTRCP(K) < 1) CYCLE

         !==============================================================
         ! Global total emission
         !==============================================================

         TOTAL = 0d0
         
         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2] * seconds in the month / Avogadro's number
            ! Also multiply by the factor 1d-9 to convert kg to Tg
            A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) 
     &           / 6.0225d23
            
            ! Anthro emissions
            TOTAL = TOTAL + SUM(RCP_LAND(:,J,K)) * A * 
     &           TRACER_MW_KG(IDTRCP(K)) 
            
         ENDDO

         !==============================================================
         !  Units
         !==============================================================

         SELECT CASE ( TRACER_NAME(IDTRCP(K)) )
         CASE ( 'NOx' ) 
            UNIT='NO2'
         CASE ( 'CO'  ) 
            UNIT=' CO'
         CASE ( 'SO2' ) 
            UNIT='SO2'
         CASE ( 'NH3' ) 
            UNIT='NH3'
         CASE DEFAULT
            UNIT='  C'
         END SELECT

         !==============================================================
         !  Print info
         !==============================================================

         IF ( TOTAL > 0d0 ) THEN
            WRITE( 6, 100 ) TRACER_NAME(IDTRCP(K)), THISMONTH, 
     &           TOTAL, UNIT
         ENDIF
 100     FORMAT( 'Anthro Land ', a4, ' for month ',
     &           i2.2, ': ',  f13.6, ' Tg ', a3 )

      ENDDO

      !==============================================================
      ! RCP Ship emissions
      !==============================================================

      WRITE( 6, '(a)' )
      DO K=1, SIZE(IDTRCP_SHIP)
         
         IF (IDTRCP_SHIP(K) < 1) CYCLE

         !==============================================================
         ! Global total emission
         !==============================================================

         TOTAL = 0d0
         
         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2] * seconds in the month / Avogadro's number
            ! Also multiply by the factor 1d-9 to convert kg to Tg
            A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) 
     &           / 6.0225d23
            
            ! Anthro emissions
            TOTAL = TOTAL + SUM(RCP_SHIP(:,J,K)) * A * 
     &           TRACER_MW_KG(IDTRCP_SHIP(K)) 
            
         ENDDO

         !==============================================================
         !  Units
         !==============================================================

         SELECT CASE ( TRACER_NAME(IDTRCP_SHIP(K)) )
         CASE ( 'NOx' ) 
            UNIT='NO2'
         CASE ( 'CO'  ) 
            UNIT=' CO'
         CASE ( 'SO2' ) 
            UNIT='SO2'
         CASE ( 'NH3' ) 
            UNIT='NH3'
         CASE DEFAULT
            UNIT='  C'
         END SELECT

         !==============================================================
         !  Print info
         !==============================================================

         IF ( TOTAL > 0d0 ) THEN
            WRITE( 6, 101 ) TRACER_NAME(IDTRCP_SHIP(K)), THISMONTH, 
     &           TOTAL, UNIT
         ENDIF
 101     FORMAT( 'Anthro Ship ', a4, ' for month ',
     &           i2.2, ': ',  f13.6, ' Tg ', a3 )

      ENDDO

!      WRITE(6,*)'ALL OK'
!      CALL GEOS_CHEM_STOP

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_RCP

!-----------------------------------------------------------------------

      SUBROUTINE READ_RCP_SURFACE( RCPNAME, YEAR, MONTH, SPECIES, 
     &     LAND, SHIP, EMISSION )

      USE GRID_MOD,               ONLY : GET_XEDGE, GET_YEDGE
      USE RCP_MOD,                ONLY : READ_AND_REGRID_SURFACE
      USE CMN_SIZE_MOD               ! Size parameters

      ! Arguments
      CHARACTER(LEN=*),INTENT(IN)   :: RCPNAME, YEAR, SPECIES
      INTEGER,INTENT(IN)            :: MONTH
      LOGICAL,INTENT(IN)            :: LAND, SHIP
      REAL*4,DIMENSION(:,:),INTENT(OUT) :: EMISSION
      

      !Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmpEmission

      REAL*8, SAVE :: XEDGE(IIPAR+1)
      REAL*8, SAVE :: YEDGE(JJPAR+1)
      INTEGER :: I
      CHARACTER(LEN=255) :: FNAME

      !=================================================================
      !  READ_RCP_SURFACE begins here
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Create array of longitude edges
         DO I=1, IIPAR
            XEDGE(I) = get_XEDGE( I )
         ENDDO
         XEDGE(IIPAR+1) = XEDGE(1)

         ! Create array of latitude edges
         DO I=1, JJPAR
            yedge(I) = get_yedge( I )
         ENDDO
         yedge(JJPAR+1) = 90d0
         

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      IF ( LAND ) THEN

         FNAME = '/gdata/prather/winds/RCP/' // trim( RCPNAME ) // 
     &        '/' //
     &        trim( YEAR ) // '/' // trim( RCPNAME ) // '_' //
     &        trim( Species ) // '_surface_' // trim( YEAR ) // '.nc'

      ELSE
         
         FNAME = '/gdata/prather/winds/RCP/' // trim( RCPNAME ) // '/' 
     &        // trim( YEAR ) // '/' // trim( RCPNAME ) // '_' //
     &        trim( Species ) // '_ships_' // trim( YEAR ) // '.nc'

      ENDIF


      IF ( LAND ) THEN
         
         CALL READ_AND_REGRID_SURFACE( 
     &        FILENAME=FNAME, 
     &        varname1='emiss_nrd',
     &        varname2='emiss_tra',
     &        outLonEdge=XEDGE,
     &        outLatEdge=YEDGE,
     &        nMoments=1,
     &        outEmission=tmpEmission )

      ELSE 

         CALL READ_AND_REGRID_SURFACE( 
     &        FILENAME=FNAME, 
     &        varname1='emiss_shp',
     &        outLonEdge=XEDGE,
     &        outLatEdge=YEDGE,
     &        nMoments=1,
     &        outEmission=tmpEmission )
         
      ENDIF

      EMISSION = tmpEmission(:,:,1,MONTH) ! molec/cm2/s



      END SUBROUTINE READ_RCP_SURFACE


!----------------------------------------------------------------------

      FUNCTION GET_RCP_EMISSION( I, J, N, LAND, SHIP ) 
     &         RESULT( EMISS )

      USE TRACERID_MOD
      USE ERROR_MOD,            ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J, N !GEOS-Chem advected tracer index
      LOGICAL, INTENT(IN), OPTIONAL :: SHIP, LAND 

      ! Local variables
      REAL*8                        :: EMISS
      CHARACTER(LEN=20)             :: STR
      LOGICAL                       :: DOLAND, DOSHIP
      INTEGER                       :: K

      !=================================================================
      ! GET_RCP_EMISSION begins here!
      !=================================================================

      IF ( PRESENT( LAND ) ) THEN
         DOLAND = LAND
      ELSE
         DOLAND = .FALSE.
      ENDIF

      IF ( PRESENT( SHIP ) ) THEN
         DOSHIP = SHIP
      ELSE
         DOSHIP = .FALSE.
      ENDIF

      IF ( .NOT. (DOLAND .OR. DOSHIP) ) THEN
         WRITE( STR, '(I4)' ) N
         CALL ERROR_STOP( 'No land/ship emissions, tracer '//trim(STR),
     &        'GET_RCP_EMISSION' )
      ENDIF
      
      EMISS = 0d0

      IF ( DOLAND ) THEN
         DO K=1, SIZE(IDTRCP)
            IF (N == IDTRCP(K)) THEN
               EMISS = EMISS + RCP_LAND(I,J,K)
            ENDIF
         ENDDO
      ENDIF

      IF ( DOSHIP ) THEN
         DO K=1, SIZE(IDTRCP_SHIP)
            IF (N == IDTRCP_SHIP(K)) THEN
               EMISS = EMISS + RCP_SHIP(I,J,K)
            ENDIF
         ENDDO
      ENDIF

      
      ! Return to calling program
      END FUNCTION GET_RCP_EMISSION

!------------------------------------------------------------------------------

      SUBROUTINE INIT_RCP

!***********************************************************************
!  Subroutine allocates and zeroes all module arrays
!***********************************************************************

      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LRCP

      USE CMN_SIZE_MOD      ! Size parameters

      ! Local variables
      INTEGER :: AS

      !=================================================================
      !  INIT_RCP begins here
      !=================================================================

      ! Return if we LRCP = .FALSE.
      IF (.not. LRCP ) RETURN

      IDTRCP = 0d0
      IDTRCP_SHIP = 0d0

      ALLOCATE( RCP_NO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NO' )
      RCP_NO = 0e0

      ALLOCATE( RCP_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_CO' )
      RCP_CO = 0e0

      ALLOCATE( RCP_OC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_OC' )
      RCP_OC = 0e0

      ALLOCATE( RCP_BC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BC' )
      RCP_BC = 0e0

      ALLOCATE( RCP_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_SO2' )
      RCP_SO2 = 0e0

      ALLOCATE( RCP_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NH3' )
      RCP_NH3 = 0e0

      ALLOCATE( RCP_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ALK4' )
      RCP_ALK4 = 0e0

      ALLOCATE( RCP_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ACET' )
      RCP_ACET = 0e0

      ALLOCATE( RCP_MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_MEK' )
      RCP_MEK = 0e0

      ALLOCATE( RCP_ALD2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ALD2' )
      RCP_ALD2 = 0e0

      ALLOCATE( RCP_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_PRPE' )
      RCP_PRPE = 0e0

      ALLOCATE( RCP_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C3H8' )
      RCP_C3H8 = 0e0

      ALLOCATE( RCP_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_CH2O' )
      RCP_CH2O = 0e0

      ALLOCATE( RCP_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H6' )
      RCP_C2H6 = 0e0

      ALLOCATE( RCP_BENZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BENZ' )
      RCP_BENZ = 0e0

      ALLOCATE( RCP_TOLU( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_TOLU' )
      RCP_TOLU = 0e0

      ALLOCATE( RCP_XYLE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_XYLE' )
      RCP_XYLE = 0e0

      ALLOCATE( RCP_C2H4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H4' )
      RCP_C2H4 = 0e0

      ALLOCATE( RCP_C2H2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H2' )
      RCP_C2H2 = 0e0

      ALLOCATE( RCP_HCOOH( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_HCOOH' )
      RCP_HCOOH = 0e0

      ! Aircraft emissions
      ALLOCATE( RCP_NO_AIR( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NO_AIR' )
      RCP_NO_AIR = 0e0

      ALLOCATE( RCP_BC_AIR( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BC_AIR' )
      RCP_BC_AIR = 0e0

      ! Ship emissions
      ALLOCATE( RCP_NO_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NO_SHIP' )
      RCP_NO_SHIP = 0e0

      ALLOCATE( RCP_CO_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_CO_SHIP' )
      RCP_CO_SHIP = 0e0

      ALLOCATE( RCP_OC_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_OC_SHIP' )
      RCP_OC_SHIP = 0e0

      ALLOCATE( RCP_BC_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BC_SHIP' )
      RCP_BC_SHIP = 0e0

      ALLOCATE( RCP_SO2_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_SO2_SHIP' )
      RCP_SO2_SHIP = 0e0

      ALLOCATE( RCP_NH3_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NH3_SHIP' )
      RCP_NH3_SHIP = 0e0

      ALLOCATE( RCP_ALK4_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ALK4_SHIP' )
      RCP_ALK4_SHIP = 0e0

      ALLOCATE( RCP_ACET_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ACET_SHIP' )
      RCP_ACET_SHIP = 0e0

      ALLOCATE( RCP_MEK_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_MEK_SHIP' )
      RCP_MEK_SHIP = 0e0

      ALLOCATE( RCP_ALD2_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_ALD2_SHIP' )
      RCP_ALD2_SHIP = 0e0

      ALLOCATE( RCP_PRPE_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_PRPE_SHIP' )
      RCP_PRPE_SHIP = 0e0

      ALLOCATE( RCP_C3H8_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C3H8_SHIP' )
      RCP_C3H8_SHIP = 0e0

      ALLOCATE( RCP_CH2O_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_CH2O_SHIP' )
      RCP_CH2O_SHIP = 0e0

      ALLOCATE( RCP_C2H6_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H6_SHIP' )
      RCP_C2H6_SHIP = 0e0

      ALLOCATE( RCP_BENZ_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BENZ_SHIP' )
      RCP_BENZ_SHIP = 0e0

      ALLOCATE( RCP_TOLU_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_TOLU_SHIP' )
      RCP_TOLU_SHIP = 0e0

      ALLOCATE( RCP_XYLE_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_XYLE_SHIP' )
      RCP_XYLE_SHIP = 0e0

      ALLOCATE( RCP_C2H4_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H4_SHIP' )
      RCP_C2H4_SHIP = 0e0

      ALLOCATE( RCP_C2H2_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_C2H2_SHIP' )
      RCP_C2H2_SHIP = 0e0

      ALLOCATE( RCP_HCOOH_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_HCOOH_SHIP' )
      RCP_HCOOH_SHIP = 0e0


      ! Alternative method
      ALLOCATE( RCP_LAND( IIPAR, JJPAR, SIZE(IDTRCP) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_LAND' )
      RCP_LAND = 0e0

      ALLOCATE( RCP_SHIP( IIPAR, JJPAR, SIZE(IDTRCP_SHIP) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_SHIP' )
      RCP_SHIP = 0e0

      ! Return to calling program
      END SUBROUTINE INIT_RCP

!----------------------------------------------------------------------

      SUBROUTINE CLEANUP_RCP
!**********************************************************************
!  Subroutine deallocates all module arrays
!**********************************************************************

      !=================================================================
      !  CLEANUP_RCP begins here
      !=================================================================

      ! Fossil fuel
      IF ( ALLOCATED( RCP_NO ) ) DEALLOCATE( RCP_NO )
      IF ( ALLOCATED( RCP_CO ) ) DEALLOCATE( RCP_CO )
      IF ( ALLOCATED( RCP_OC ) ) DEALLOCATE( RCP_OC )
      IF ( ALLOCATED( RCP_BC ) ) DEALLOCATE( RCP_BC )
      IF ( ALLOCATED( RCP_SO2 ) ) DEALLOCATE( RCP_SO2 )
      IF ( ALLOCATED( RCP_NH3 ) ) DEALLOCATE( RCP_NH3 )
      IF ( ALLOCATED( RCP_ALK4 ) ) DEALLOCATE( RCP_ALK4 )
      IF ( ALLOCATED( RCP_ACET ) ) DEALLOCATE( RCP_ACET )
      IF ( ALLOCATED( RCP_MEK  ) ) DEALLOCATE( RCP_MEK  )
      IF ( ALLOCATED( RCP_ALD2 ) ) DEALLOCATE( RCP_ALD2 )
      IF ( ALLOCATED( RCP_PRPE ) ) DEALLOCATE( RCP_PRPE )
      IF ( ALLOCATED( RCP_C3H8 ) ) DEALLOCATE( RCP_C3H8 )
      IF ( ALLOCATED( RCP_CH2O ) ) DEALLOCATE( RCP_CH2O )
      IF ( ALLOCATED( RCP_C2H6 ) ) DEALLOCATE( RCP_C2H6 )
      IF ( ALLOCATED( RCP_BENZ ) ) DEALLOCATE( RCP_BENZ )
      IF ( ALLOCATED( RCP_TOLU ) ) DEALLOCATE( RCP_TOLU )
      IF ( ALLOCATED( RCP_XYLE ) ) DEALLOCATE( RCP_XYLE )
      IF ( ALLOCATED( RCP_C2H4 ) ) DEALLOCATE( RCP_C2H4 )
      IF ( ALLOCATED( RCP_C2H2 ) ) DEALLOCATE( RCP_C2H2 )
      IF ( ALLOCATED( RCP_HCOOH ) ) DEALLOCATE( RCP_HCOOH )

      IF ( ALLOCATED( RCP_NO_AIR ) ) DEALLOCATE( RCP_NO_AIR )
      IF ( ALLOCATED( RCP_BC_AIR ) ) DEALLOCATE( RCP_BC_AIR )

      IF ( ALLOCATED( RCP_NO_SHIP   ) ) DEALLOCATE( RCP_NO_SHIP )
      IF ( ALLOCATED( RCP_CO_SHIP   ) ) DEALLOCATE( RCP_CO_SHIP )
      IF ( ALLOCATED( RCP_OC_SHIP   ) ) DEALLOCATE( RCP_OC_SHIP )
      IF ( ALLOCATED( RCP_BC_SHIP   ) ) DEALLOCATE( RCP_BC_SHIP )
      IF ( ALLOCATED( RCP_SO2_SHIP  ) ) DEALLOCATE( RCP_SO2_SHIP )
      IF ( ALLOCATED( RCP_NH3_SHIP  ) ) DEALLOCATE( RCP_NH3_SHIP )
      IF ( ALLOCATED( RCP_ALK4_SHIP ) ) DEALLOCATE( RCP_ALK4_SHIP )
      IF ( ALLOCATED( RCP_ACET_SHIP ) ) DEALLOCATE( RCP_ACET_SHIP )
      IF ( ALLOCATED( RCP_MEK_SHIP  ) ) DEALLOCATE( RCP_MEK_SHIP )
      IF ( ALLOCATED( RCP_ALD2_SHIP ) ) DEALLOCATE( RCP_ALD2_SHIP )
      IF ( ALLOCATED( RCP_PRPE_SHIP ) ) DEALLOCATE( RCP_PRPE_SHIP )
      IF ( ALLOCATED( RCP_C3H8_SHIP ) ) DEALLOCATE( RCP_C3H8_SHIP )
      IF ( ALLOCATED( RCP_C2H6_SHIP ) ) DEALLOCATE( RCP_C2H6_SHIP )
      IF ( ALLOCATED( RCP_CH2O_SHIP ) ) DEALLOCATE( RCP_CH2O_SHIP )
      IF ( ALLOCATED( RCP_BENZ_SHIP ) ) DEALLOCATE( RCP_BENZ_SHIP )
      IF ( ALLOCATED( RCP_TOLU_SHIP ) ) DEALLOCATE( RCP_TOLU_SHIP )
      IF ( ALLOCATED( RCP_XYLE_SHIP ) ) DEALLOCATE( RCP_XYLE_SHIP )
      IF ( ALLOCATED( RCP_C2H4_SHIP ) ) DEALLOCATE( RCP_C2H4_SHIP )
      IF ( ALLOCATED( RCP_C2H2_SHIP ) ) DEALLOCATE( RCP_C2H2_SHIP )
      IF ( ALLOCATED( RCP_HCOOH_SHIP ) ) DEALLOCATE( RCP_HCOOH_SHIP )

      IF ( ALLOCATED( RCP_LAND ) ) DEALLOCATE( RCP_LAND )
      IF ( ALLOCATED( RCP_SHIP ) ) DEALLOCATE( RCP_SHIP )

      ! Return to calling program
      END SUBROUTINE CLEANUP_RCP

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE EMISS_RCP_OLD

!***********************************************************************
!  Subroutine EMISS_RCP reads all RCP emissions at the beginning of
!  each month. (cdh, 10/14/11)
!***********************************************************************

      USE TIME_MOD,               ONLY : GET_MONTH
      USE ERROR_MOD,              ONLY : GEOS_CHEM_STOP
      USE DIRECTORY_MOD,          ONLY : DATA_DIR

      USE BPCH2_MOD,              ONLY : GET_TAU0, GET_RES_EXT
      USE TRACERID_MOD
      USE TRACER_MOD, ONLY      : TRACER_NAME
      USE COMODE_LOOP_MOD, ONLY : IDEMS

      !Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: THISMONTH
      CHARACTER(LEN=20)             :: RCPSPECIES
      INTEGER :: I
      CHARACTER(LEN=255)            :: FILENAME
      REAL*8                        :: XTAU


      !=================================================================
      !  EMISS_RCP begins here
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays
         CALL INIT_RCP

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      ! Get month
      THISMONTH = GET_MONTH()      

      ! Echo info
      WRITE(6, '(a)' ) REPEAT( '=', 79)
      WRITE(6, 100   ) RCPNAME, RCPYEAR
 100  FORMAT( 'R C P    E M I S S I O N S',
     &        '  --  Scenario: ', 2A10 , / )

      ! Anthro file name
      FILENAME = TRIM( DATA_DIR )        // 'RCP_201110/'     //
     &                  trim( RCPNAME )  // '_anthropogenic_' // 
     &                  trim( RCPYEAR )  // '.'               //
     &                  GET_RES_EXT()    // '.bpch' 

      XTAU = GET_TAU0( 1, 1, 2000 ) 
      ! Read data
      CALL READ_RCP( FILENAME, SHIP=.FALSE. )!, TAU0=XTAU )

      ! Anthro file name
      FILENAME = TRIM( DATA_DIR )        // 'RCP_201110/'     //
     &                  trim( RCPNAME )  // '_ships_' // 
     &                  trim( RCPYEAR )  // '.'               //
     &                  GET_RES_EXT()    // '.bpch' 

      XTAU = GET_TAU0( THISMONTH, 1, 2000 )  

      ! Read data
      CALL READ_RCP( FILENAME, SHIP=.TRUE., TAU0=XTAU )

        
      !=================================================================
      ! Aircraft emissions
      !=================================================================

      

      FILENAME = TRIM( DATA_DIR )        // 'RCP_201110/'     //
     &                  trim( RCPNAME )  // '_NO_aircraft_'   // 
     &                  trim( RCPYEAR )  // '.nc'

      CALL READ_RCP_AIRCRAFT( FILENAME, THISMONTH )

      ! Print totals to log
      CALL TOTAL_ANTHRO_RCP( THISMONTH )

      ! Fancy output
      WRITE(6, '(a)' ) REPEAT( '=', 79)

!      PRINT*,'ALL OK'
!      CALL GEOS_CHEM_STOP

      ! Return to calling program
      END SUBROUTINE EMISS_RCP_OLD

!----------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_RCP_OLD( THISMONTH )

!***********************************************************************
!  Subroutine to print total RCP anthropogenic emissions each month
!  in TgC
!***********************************************************************

      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACER_MOD, ONLY      : TRACER_NAME
      USE COMODE_LOOP_MOD, ONLY : IDEMS

      USE CMN_SIZE_MOD       ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J, K
      REAL*8              :: A, TOTAL, TOTAL_SHIP
      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      !  TOTAL_ANTHRO_RCP begins here
      !=================================================================

      !=================================================================
      !  Sum anthropogenic emissions
      !=================================================================

      WRITE( 6, '(a)' )
      DO K=1, SIZE(IDEMS)
         
         IF (IDEMS(K) == 0) CYCLE

         TOTAL = 0d0
         TOTAL_SHIP = 0d0
         
         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2] * seconds in the month / Avogadro's number
            ! Also multiply by the factor 1d-9 to convert kg to Tg
            A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) 
     &           / 6.0225d23
            
            ! Anthro emissions
            TOTAL = TOTAL + SUM(RCP_LAND(:,J,K)) * A * 
     &           TRACER_MW_KG(IDEMS(K)) 
            TOTAL_SHIP = TOTAL_SHIP + SUM(RCP_SHIP(:,J,K)) * A * 
     &           TRACER_MW_KG(IDEMS(K)) 
            
         ENDDO

         SELECT CASE ( TRACER_NAME(IDEMS(K)) )
         CASE ( 'NOx' ) 
            UNIT='NO2'
         CASE ( 'CO'  ) 
            UNIT=' CO'
         CASE ( 'SO2' ) 
            UNIT='SO2'
         CASE ( 'NH3' ) 
            UNIT='NH3'
         CASE DEFAULT
            UNIT='  C'
         END SELECT

         IF ( (TOTAL > 0d0) .OR. (TOTAL_SHIP > 0d0) ) THEN
            WRITE( 6, 100 ) TRACER_NAME(IDEMS(K)), THISMONTH, 
     &           TOTAL, TOTAL_SHIP, UNIT
         ENDIF
 100     FORMAT( 'Anthro Land/Ship ', a4, ' for 2000/',
     &           i2.2, ': ', 2f13.6, ' Tg', a3 )

      ENDDO


      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_RCP_OLD

!----------------------------------------------------------------------

      FUNCTION GET_RCP_EMISSION_OLD( I, J, N, LAND, SHIP ) 
     &         RESULT( EMISS )

      USE TRACERID_MOD
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE COMODE_LOOP_MOD,      ONLY : IDEMS

      ! Arguments
      INTEGER, INTENT(IN)           :: I, J, N !GEOS-Chem advected tracer index
      LOGICAL, INTENT(IN), OPTIONAL :: SHIP, LAND 

      ! Local variables
      REAL*8                        :: EMISS
      CHARACTER(LEN=20)             :: STR
      LOGICAL                       :: DOLAND, DOSHIP
      INTEGER                       :: K

      !=================================================================
      ! GET_RCP_EMISSION begins here!
      !=================================================================

      IF ( PRESENT( LAND ) ) THEN
         DOLAND = LAND
      ELSE
         DOLAND = .FALSE.
      ENDIF

      IF ( PRESENT( SHIP ) ) THEN
         DOSHIP = SHIP
      ELSE
         DOSHIP = .FALSE.
      ENDIF

      IF ( .NOT. (DOLAND .OR. DOSHIP) ) THEN
         WRITE( STR, '(I4)' ) N
         CALL ERROR_STOP( 'No land/ship emissions, tracer '//trim(STR),
     &        'GET_RCP_EMISSION' )
      ENDIF
      
      EMISS = 0d0

      DO K=1, SIZE(IDEMS)

         IF (N == IDEMS(K)) THEN

            IF ( DOLAND ) THEN
               EMISS = RCP_LAND(I,J,K)
            ENDIF

            IF ( DOSHIP ) THEN
               EMISS = EMISS + RCP_SHIP(I,J,K)
            ENDIF

         ENDIF

      ENDDO

      
      ! Return to calling program
      END FUNCTION GET_RCP_EMISSION_OLD

!------------------------------------------------------------------------------


      END MODULE RCP2_MOD
