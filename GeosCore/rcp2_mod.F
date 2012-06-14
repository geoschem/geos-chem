
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
      REAL*4, ALLOCATABLE :: RCP_LAND(:,:,:)

      ! Aircraft emissions
      REAL*4, ALLOCATABLE :: RCP_NO_AIR(:,:,:)
      REAL*4, ALLOCATABLE :: RCP_BC_AIR(:,:,:)

      ! Ship emissions
      REAL*4, ALLOCATABLE :: RCP_SHIP(:,:,:)

      CHARACTER(LEN=20)   :: RCPNAME, RCPYEAR
      INTEGER             :: IDTRCP_LAND(20), IDTRCP_SHIP(20)

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
      USE ERROR_MOD,        ONLY : ERROR_STOP
      USE TRACERID_MOD           ! tracer ID numbers
      USE CMN_SIZE_MOD           ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      LOGICAL,          INTENT(IN) :: SHIP
      REAL*8,OPTIONAL,  INTENT(IN) :: TAU0


      ! Local variables
      INTEGER                         :: I, J, L, N, IOS, K, IDT 
      INTEGER                         :: NI, NJ, NL
      INTEGER                         :: IFIRST, JFIRST, LFIRST
      INTEGER                         :: NTRACER, NSKIP
      INTEGER                         :: HALFPOLAR, CENTER180
      INTEGER                         :: SCALEYEAR, BASEYEAR
      REAL*4                          :: LONRES, LATRES
      REAL*4                          :: ARRAY(IGLOB,JGLOB,1)
      REAL*4                          :: TMP(IIPAR,JJPAR)
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
            
         IDT = 0

         ! Find GEOS-Chem tracer ID for each species in file
         ! In most cases these ID numbers will be the same, but we do this
         ! in case the GEOS-Chem tracer numbers change in the future
         SELECT CASE ( NTRACER )
         CASE ( 1  )
            IDT = IDTNOx
         CASE ( 4  )
            IDT = IDTCO
         CASE ( 5  )
            IDT = IDTALK4
         CASE ( 9  ) 
            ! We expect ACET to be lumped with MEK, as explained below
            ! and in RETRO implementation
            CALL ERROR_STOP( 'RCP file unexpectely contains ACET: ' // 
     &           FILENAME, 'READ_RCP ' )
!            IDT = IDTACET
         CASE ( 10 )
            IDT = IDTMEK
         CASE ( 11 ) 
            IDT = IDTALD2
         CASE ( 18 ) 
            IDT = IDTPRPE
         CASE ( 19 ) 
            IDT = IDTC3H8
         CASE ( 20 ) 
            IDT = IDTCH2O
         CASE ( 21 ) 
            IDT = IDTC2H6
         CASE ( 26 ) 
            IDT = IDTSO2
         CASE ( 30 ) 
            IDT = IDTNH3
         CASE ( 36 ) 
            IDT = IDTBCPO
         CASE ( 37 )
            IDT = IDTOCPO
         CASE ( 59 )
            IDT = IDTBENZ
         CASE ( 60 ) 
            IDT = IDTTOLU
         CASE ( 61 ) 
            IDT = IDTXYLE
         CASE ( 65 ) 
            IDT = IDTC2H4
         CASE ( 66 )
            IDT = IDTC2H2
         CASE ( 76 )
            IDT = IDTHCOOH
         CASE DEFAULT
            ! DO NOTHING
         END SELECT

         ! Tracer number must be positive, 
         ! otherwise it's not used for this simulation type
         IF ( IDT > 0 ) THEN

            ! Increment tracer counter
            K = K + 1 

            ! Save emissions and tracer number
            IF (SHIP) THEN
               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_SHIP(:,:,K) )
               IDTRCP_SHIP(K) = IDT 
            ELSE
               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_LAND(:,:,K) )
               IDTRCP_LAND(K) = IDT
            ENDIF

            !==============================================================
            ! Special case for MEK
            ! Partition ketones into 75% acetone and 25% MEK
            ! In the file, MEK contains all ketones.
            ! As done for RETRO (cdh, 10/18/11; dbm, 8/18/2011)
            !==============================================================
            
            IF (IDT == IDTMEK) THEN
               
               ! Reduce MEK emissions
               IF (SHIP) THEN  
                  RCP_SHIP(:,:,K) = RCP_SHIP(:,:,K) * 0.25D0
               ELSE
                  RCP_LAND(:,:,K) = RCP_LAND(:,:,K) * 0.25D0
               ENDIF

               IF (IDTACET > 0d0) THEN

                  ! Increment tracer counter
                  K = K + 1 
                  
                  ! Save ACET emissions (75% of original MEK = 3*25%
                  IF (SHIP) THEN  
                     RCP_SHIP(:,:,K) = RCP_SHIP(:,:,K-1) * 3d0
                     IDTRCP_SHIP(K)  = IDTACET 
                  ELSE
                     RCP_LAND(:,:,K) = RCP_LAND(:,:,K-1) * 3d0
                     IDTRCP_LAND(K)  = IDTACET
                  ENDIF

               ENDIF


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
         XEDGE(IIPAR+1) = XEDGE(1) + 360d0
         XEDGE = XEDGE + 360d0

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
      USE AIRCRAFT_NOX_MOD,       ONLY : EMIS_AC_NOx, READAIR
      USE DIAG_MOD,               ONLY : AD32_AC
      USE CMN_DIAG_MOD       ! Diagnostic switches
      USE CMN_SIZE_MOD       ! Size parameters      

      ! Local variables
      INTEGER          :: I, J, L
      LOGICAL, SAVE    :: FIRST=.TRUE.

      !=================================================================
      !  RCP_AIREMISS begins here
      !=================================================================

      ! Allocate and initialize arrays
      IF ( FIRST ) THEN 
         CALL READAIR ! use this only because init_aircraft_nox is private
         FIRST = .FALSE.
      ENDIF

      ! Convert  -> molec/cm3/s
      EMIS_AC_NOx = RCP_NO_AIR

      ! ND32 -- save NOx in [molec/cm2], will convert to
      ! [molec/cm2/s] in subroutine "diag3.f" (bmy, 3/16/00)
      IF ( ND32 > 0 ) THEN
         DO L=1, LLPAR
         DO J=1, JJPAR
         DO I=1, IIPAR
            AD32_ac(I,J,L) = AD32_ac(I,J,L) + 
     &           ( EMIS_AC_NOx(I,J,L) * BXHEIGHT(I,J,L) * 1d2 )
         ENDDO
         ENDDO
         ENDDO
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
      DO K=1, SIZE(IDTRCP_LAND)
         
         IF (IDTRCP_LAND(K) < 1) CYCLE

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
     &           TRACER_MW_KG(IDTRCP_LAND(K)) 
            
         ENDDO

         !==============================================================
         !  Units
         !==============================================================

         SELECT CASE ( TRACER_NAME(IDTRCP_LAND(K)) )
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

         WRITE( 6, 100 ) TRACER_NAME(IDTRCP_LAND(K)), THISMONTH, 
     &        TOTAL, UNIT
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

         WRITE( 6, 101 ) TRACER_NAME(IDTRCP_SHIP(K)), THISMONTH, 
     &           TOTAL, UNIT
 101     FORMAT( 'Anthro Ship ', a4, ' for month ',
     &           i2.2, ': ',  f13.6, ' Tg ', a3 )

      ENDDO

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
         XEDGE(IIPAR+1) = XEDGE(1) + 360d0
         XEDGE = XEDGE + 360d0

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
      LOGICAL                       :: DOLAND, DOSHIP, TRACERFOUND
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
      TRACERFOUND = .FALSE.

      IF ( DOLAND ) THEN
         DO K=1, SIZE(IDTRCP_LAND)
            IF (N == IDTRCP_LAND(K)) THEN
               EMISS = EMISS + RCP_LAND(I,J,K)
               TRACERFOUND=.TRUE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      IF ( DOSHIP ) THEN
         DO K=1, SIZE(IDTRCP_SHIP)
            IF (N == IDTRCP_SHIP(K)) THEN
               EMISS = EMISS + RCP_SHIP(I,J,K)
               TRACERFOUND=.TRUE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! Return -1 if there are no emissions for tracer N
      IF (.NOT. TRACERFOUND) EMISS = -1d0
      
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

      IDTRCP_LAND = 0d0
      IDTRCP_SHIP = 0d0

      ! Aircraft emissions
      ALLOCATE( RCP_NO_AIR( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_NO_AIR' )
      RCP_NO_AIR = 0e0

      ALLOCATE( RCP_BC_AIR( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_BC_AIR' )
      RCP_BC_AIR = 0e0

      ! Anthropogenic land surface emissions
      ALLOCATE( RCP_LAND( IIPAR, JJPAR, SIZE(IDTRCP_LAND) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_LAND' )
      RCP_LAND = 0e0

      ! Shipping 
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
      IF ( ALLOCATED( RCP_NO_AIR ) ) DEALLOCATE( RCP_NO_AIR )
      IF ( ALLOCATED( RCP_BC_AIR ) ) DEALLOCATE( RCP_BC_AIR )
      IF ( ALLOCATED( RCP_LAND ) ) DEALLOCATE( RCP_LAND )
      IF ( ALLOCATED( RCP_SHIP ) ) DEALLOCATE( RCP_SHIP )

      ! Return to calling program
      END SUBROUTINE CLEANUP_RCP

!----------------------------------------------------------------------


      END MODULE RCP2_MOD
