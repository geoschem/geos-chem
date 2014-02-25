!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: rcp_mod
!
! !DESCRIPTION: Module RCP\_MOD provides access to the RCP emission inventories
!  that were prepared for IPCC AR5. The inventory includes anthropogenic
!  emissions from land, ships, and aircraft. Species include trace gases 
!  (NOx, CO, NH3, SO2, various VOCs) and aerosols (BC, OC). Land emissions
!  include fossil fuel and biofuel use, energy production and distribution, 
!  residential and commercial combustion, industry, transportation, waste
!  treatment and disposal, solvent production and use, agriculture, and
!  agricultural waste burning. Data sources are documented in the data
!  directories. 
!\\
!\\
! !INTERFACE: 
!
      MODULE RCP_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      !NONE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: CLEANUP_RCP
      PUBLIC :: LOAD_RCP_EMISSIONS
      PUBLIC :: GET_RCP_EMISSION
      PUBLIC :: RCPNAME, RCPYEAR
      PUBLIC :: RCP_AIREMISS
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to pass aircraft NOx emissions to SMVGEAR
      REAL*8,  PUBLIC, ALLOCATABLE :: RCP_AC_NOx(:,:,:)
!
! !PRIVATE DATA MEMBERS:
!
      REAL*4, ALLOCATABLE :: RCP_LAND(:,:,:)
      REAL*4, ALLOCATABLE :: RCP_AIR(:,:,:,:)
      REAL*4, ALLOCATABLE :: RCP_SHIP(:,:,:)
      CHARACTER(LEN=20)   :: RCPNAME
      INTEGER             :: RCPYEAR
      INTEGER             :: IDTRCP_LAND(20), IDTRCP_SHIP(20),
     &     IDTRCP_AIR(3)
!
! !REVISION HISTORY:
!  14 Jun 2012 - C. Holmes   - Initial version
!  31 Jul 2013 - M. Sulprizio- Added new array RCP_AC_NOx to pass aircraft
!                              emissions to SMVGEAR. The array EMIS_AC_NOx
!                              was removed with the obsolete aircraft_nox_mod.F
!                              when AEIC aircraft emissions were added.
!EOP
!------------------------------------------------------------------------------
!BOC

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_rcp_emissions
!
! !DESCRIPTION: Subroutine LOAD\_RCP\_EMISSIONS reads all RCP emissions at the
!  beginning of each month. (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LOAD_RCP_EMISSIONS( Input_Opt )
!
! !USES:
!
      USE BPCH2_MOD,              ONLY : GET_TAU0, GET_RES_EXT
      USE DIRECTORY_MOD,          ONLY : DATA_DIR
      USE ERROR_MOD,              ONLY : GEOS_CHEM_STOP
      USE GIGC_Input_Opt_Mod,     ONLY : OptInput
      USE TIME_MOD,               ONLY : GET_MONTH
      USE TRACERID_MOD
      USE TRACER_MOD,             ONLY : TRACER_NAME
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP, LRCPSHIP, LRCPAIR from Input_Opt
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: THISMONTH, I
      CHARACTER(LEN=20)             :: RCPSPECIES, YEARSTR
      CHARACTER(LEN=255)            :: FILENAME
      REAL*8                        :: XTAU

      ! For fields from Input_Opt
      LOGICAL                       :: LRCP
      LOGICAL                       :: LRCPSHIP
      LOGICAL                       :: LRCPAIR

      !=================================================================
      !  LOAD_RCP_EMISSIONS begins here
      !=================================================================

      ! Copy values from Input_Opt
      LRCP     = Input_Opt%LRCP
      LRCPSHIP = Input_Opt%LRCPSHIP
      LRCPAIR  = Input_Opt%LRCPAIR

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays
         CALL INIT_RCP( Input_Opt )

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      ! Get month
      THISMONTH = GET_MONTH()      
      
      ! Convert to string
      WRITE( YEARSTR, '(I4)' ) RCPYEAR

      !=================================================================
      ! Land and ship emissions
      !=================================================================

      IF( LRCP .OR. LRCPSHIP) THEN

         ! Land file name
         FILENAME = TRIM( DATA_DIR )        // 'RCP_201206/'     //
     &                     trim( RCPNAME )  // '/'               //
     &                     trim( RCPNAME )  // '_anthropogenic_' // 
     &                     trim( YEARSTR )  // '.'               //
     &                     GET_RES_EXT()    // '.bpch' 
         ! Date for emissions
         ! Land emissions dated Jan 1 because all months are the same
         XTAU = GET_TAU0( 1, 1, RCPYEAR ) 

         ! Read data (LAND -> TYPE=1)
         CALL READ_RCP_BPCH( FILENAME, TYPE=1, TAU0=XTAU )

         ! Ship file name
         FILENAME = TRIM( DATA_DIR )        // 'RCP_201206/'     //
     &                     trim( RCPNAME )  // '/'               //
     &                     trim( RCPNAME )  // '_ships_'         // 
     &                     trim( YEARSTR )  // '.'               //
     &                     GET_RES_EXT()    // '.bpch' 

         ! Date for emissions
         XTAU = GET_TAU0( THISMONTH, 1, RCPYEAR )  

         ! Read data (SHIP -> TYPE=2)
         CALL READ_RCP_BPCH( FILENAME, TYPE=2, TAU0=XTAU )
         
      ENDIF
        
      !=================================================================
      ! Aircraft emissions
      !=================================================================

      IF (LRCPAIR) THEN
         FILENAME = TRIM( DATA_DIR )     // 'RCP_201206/'     //
     &                  trim( RCPNAME )  // '/'               //
     &                  trim( RCPNAME )  // '_aircraft_'      // 
     &                  trim( YEARSTR )  // '.'               //
     &                  GET_RES_EXT()    // '.bpch' 

         ! Date for emissions
         XTAU = GET_TAU0( THISMONTH, 1, RCPYEAR ) 

         ! Read data (AIRCRAFT -> TYPE=3)
         CALL READ_RCP_BPCH( FILENAME, TYPE=3, TAU0=XTAU )

      ENDIF

      !=================================================================
      ! Print totals to log
      !=================================================================

      CALL TOTAL_ANTHRO_RCP( THISMONTH )

      ! Fancy output
      WRITE(6, '(a)' ) REPEAT( '=', 79)

      END SUBROUTINE LOAD_RCP_EMISSIONS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_rcp_bpch
!
! !DESCRIPTION: Subroutine READ\_RCP\_BPCH reads a BPCH file containing RCP
!  data. (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_RCP_BPCH( FILENAME, TYPE, TAU0 )
!
! !USES:
!
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE ERROR_MOD,        ONLY : ERROR_STOP
      USE TRACERID_MOD           ! tracer ID numbers
      USE CMN_SIZE_MOD           ! Size parameters
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      INTEGER,          INTENT(IN) :: TYPE ! 1=LAND, 2=SHIP, 3=AIRCRAFT
      REAL*8,OPTIONAL,  INTENT(IN) :: TAU0
!
! !REVISION HISTORY:
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                         :: I, J, L, N, IOS, K, IDT 
      INTEGER                         :: NI, NJ, NL
      INTEGER                         :: IFIRST, JFIRST, LFIRST
      INTEGER                         :: NTRACER, NSKIP
      INTEGER                         :: HALFPOLAR, CENTER180
      INTEGER                         :: SCALEYEAR, BASEYEAR
      REAL*4                          :: LONRES, LATRES
      REAL*4                          :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4                          :: TMP(IIPAR,JJPAR)
      REAL*8                          :: ZTAU0, ZTAU1
      CHARACTER(LEN=20)               :: MODELNAME
      CHARACTER(LEN=40)               :: CATEGORY
      CHARACTER(LEN=40)               :: UNIT
      CHARACTER(LEN=40)               :: RESERVED
      CHARACTER(LEN=20)               :: STR

      !=================================================================
      !  READ_RCP_BPCH begins here
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RCP_BPCH: Reading ', a )

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
     &        CALL ERROR_STOP( 'ANTHSRCE not found', 'READ_RCP_BPCH' ) 
         
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
         ! These ID numbers will be the same as the ID numbers
         ! stored in the files, but we do this in case the GEOS-Chem tracer
         ! numbers change in the future
         SELECT CASE ( NTRACER )
         CASE ( 1  )
            IDT = IDTNO
         CASE ( 4  )
            IDT = IDTCO
         CASE ( 5  )
            IDT = IDTALK4
         CASE ( 9  ) 
            ! We expect ACET to be lumped with MEK, as explained below
            ! and in RETRO implementation
            CALL ERROR_STOP( 'RCP file unexpectely contains ACET: ' // 
     &           FILENAME, 'READ_RCP_BPCH ' )
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
            SELECT CASE ( TYPE )
            CASE ( 1 )

               ! Error check
               IF (K > SIZE( IDTRCP_LAND )) THEN
                  WRITE( STR, '(I4)' ) K
                  CALL ERROR_STOP( 'TOO MANY SPECIES FOR RCP_LAND '//
     &                 TRIM(STR), 'READ_RCP_BPCH' )
               ENDIF

               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_LAND(:,:,K) )
               IDTRCP_LAND(K) = IDT

            CASE ( 2 )

               ! Error check
               IF (K > SIZE( IDTRCP_SHIP )) THEN
                  WRITE( STR, '(I4)' ) K
                  CALL ERROR_STOP( 'TOO MANY SPECIES FOR RCP_SHIP '//
     &                 TRIM(STR), 'READ_RCP_BPCH' )
               ENDIF

               CALL TRANSFER_2D( ARRAY(:,:,1), RCP_SHIP(:,:,K) )
               IDTRCP_SHIP(K) = IDT 

            CASE ( 3 )

               ! Error check
               IF (K > SIZE( IDTRCP_AIR )) THEN
                  WRITE( STR, '(I4)' ) K
                  CALL ERROR_STOP( 'TOO MANY SPECIES FOR RCP_AIR '//
     &                 TRIM(STR), 'READ_RCP_BPCH' )
               ENDIF

               ! Transfer, 
               DO L=1, LLPAR
                  CALL TRANSFER_2D( ARRAY(:,:,L), RCP_AIR(:,:,L,K) )
               ENDDO
               IDTRCP_AIR(K) = IDT 

            CASE DEFAULT
            END SELECT

            !==============================================================
            ! Special case for MEK
            ! Partition ketones into 75% acetone and 25% MEK
            ! In the file, MEK contains all ketones.
            ! As done for RETRO (cdh, 10/18/11; dbm, 8/18/2011)
            !==============================================================
            
            IF (IDT == IDTMEK) THEN
               
               ! Reduce MEK emissions
               SELECT CASE ( TYPE )
               CASE ( 1 )
                  RCP_LAND(:,:,K) = RCP_LAND(:,:,K) * 0.25D0
               CASE ( 2 )
                  RCP_SHIP(:,:,K) = RCP_SHIP(:,:,K) * 0.25D0
               CASE DEFAULT
                  ! No MEK emissions expected for aircraft
               END SELECT

               IF (IDTACET > 0d0) THEN

                  ! Increment tracer counter
                  K = K + 1 
                  
                  ! Save ACET emissions (75% of original MEK = 3*25%)
                  SELECT CASE ( TYPE )
                  CASE ( 1 )
                     RCP_LAND(:,:,K) = RCP_LAND(:,:,K-1) * 3d0
                     IDTRCP_LAND(K)  = IDTACET
                  CASE ( 2 )
                     RCP_SHIP(:,:,K) = RCP_SHIP(:,:,K-1) * 3d0
                     IDTRCP_SHIP(K)  = IDTACET 
                  CASE DEFAULT
                     ! No MEK emissions expected for aircraft
                  END SELECT
 
               ENDIF


            ENDIF
               
         ENDIF

      END DO

      ! Close file
      CLOSE( IU_FILE )

      END SUBROUTINE READ_RCP_BPCH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcp_airemiss
!
! !DESCRIPTION: Subroutine RCP\_AIREMISS populates EMIS\_AC\_NOx with aircraft
!  NOx emissions. Also does diagnostics. (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RCP_AIREMISS( State_Met )
!
! !USES:
!
      USE DIAG_MOD,               ONLY : AD32_AC
      USE ERROR_MOD,              ONLY : ERROR_STOP
      USE GIGC_State_Met_Mod,     ONLY : MetState
      USE TRACERID_MOD,           ONLY : IDTNO
      USE CMN_DIAG_MOD                 ! Diagnostic switches
      USE CMN_SIZE_MOD                 ! Size parameters      
!
! !INPUT PARAMETERS:
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J, L, K
      LOGICAL, SAVE    :: FIRST=.TRUE.
      LOGICAL          :: TRACERFOUND

      !=================================================================
      !  RCP_AIREMISS begins here
      !=================================================================

      ! Initialized
      TRACERFOUND = .FALSE.
      
      ! Locate the NOx tracer in the emission array
      DO K=1, SIZE(IDTRCP_AIR)
         IF (IDTNO == IDTRCP_AIR(K)) THEN
            TRACERFOUND=.TRUE.
            EXIT
         ENDIF
      ENDDO

      ! Error if there are no NOx emissions
      IF (.NOT. TRACERFOUND) 
     &     CALL ERROR_STOP('RCP AVIATION NOX HAS NOT BEEN READ', 
     &                     'RCP_AIREMISS' )
      
      ! Convert molec/cm2/s  -> molec/cm3/s
      RCP_AC_NOx = RCP_AIR(:,:,:,K) / ( State_Met%BXHEIGHT * 1D2 )

      ! ND32 -- save NOx in [molec/cm2], will convert to
      ! [molec/cm2/s] in subroutine "diag3.f" (bmy, 3/16/00)
      IF ( ND32 > 0 ) THEN
         DO L=1, LLTROP
         DO J=1, JJPAR
         DO I=1, IIPAR
            AD32_ac(I,J,L) = AD32_ac(I,J,L) + ( RCP_AC_NOx(I,J,L) * 
     &                       State_Met%BXHEIGHT(I,J,L) * 1d2 )
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE RCP_AIREMISS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_rcp
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_RCP prints total RCP anthropogenic
!  emissions each month. (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_RCP( THISMONTH )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACER_MOD,   ONLY : TRACER_NAME
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: THISMONTH
!
! !REVISION HISTORY: 
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, K
      REAL*8              :: A, TOTAL, TOTAL_SHIP

      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      !  TOTAL_ANTHRO_RCP begins here
      !=================================================================

      ! Echo info
      WRITE(6, '(a)' ) REPEAT( '=', 79)
      WRITE(6, 100   ) RCPNAME, RCPYEAR
 100  FORMAT( 'R C P    E M I S S I O N S',
     &        '  --  Scenario: ', A10, I6, / )

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
            A = GET_AREA_CM2( 1, J, 1 ) * ( D(THISMONTH) * 86400d-9 ) 
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
            ! Convert to Tg(N)
            TOTAL = TOTAL * 0.014 / TRACER_MW_KG(IDTRCP_LAND(K))
            UNIT='N'
         CASE ( 'SO2' ) 
            ! Convert to Tg(S)
            TOTAL = TOTAL * 0.032 / TRACER_MW_KG(IDTRCP_LAND(K))
            UNIT='S'
         CASE ( 'NH3' ) 
            ! Convert to Tg(N)
            TOTAL = TOTAL * 0.014 / TRACER_MW_KG(IDTRCP_LAND(K))
            UNIT='N'
         CASE ( 'CO'  ) 
            UNIT='CO'
         CASE DEFAULT
            UNIT='C'
         END SELECT

         !==============================================================
         !  Print info
         !==============================================================

         WRITE( 6, 101 ) 'Land', TRACER_NAME(IDTRCP_LAND(K)), THISMONTH, 
     &        TOTAL, UNIT
 101     FORMAT( 'Anthro ',a5, ' ', a4, ' for month ',
     &           i2.2, ': ',  f13.6, ' Tg  ', a3 )

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
            A = GET_AREA_CM2( 1, J, 1 ) * ( D(THISMONTH) * 86400d-9 ) 
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
            ! Convert to Tg(N)
            TOTAL = TOTAL *  0.014 / TRACER_MW_KG(IDTRCP_SHIP(K))
            UNIT='N'
         CASE ( 'SO2' ) 
            ! Convert to Tg(S)
            TOTAL = TOTAL * 0.032 / TRACER_MW_KG(IDTRCP_SHIP(K))
            UNIT='S'
         CASE ( 'NH3' ) 
            ! Convert to Tg(N)
            TOTAL = TOTAL * 0.014 / TRACER_MW_KG(IDTRCP_SHIP(K))
            UNIT='N'
         CASE ( 'CO'  ) 
            UNIT='CO'
         CASE DEFAULT
            UNIT='C'
         END SELECT

         !==============================================================
         !  Print info
         !==============================================================

         WRITE( 6, 101 ) 'Ship', TRACER_NAME(IDTRCP_SHIP(K)), THISMONTH, 
     &           TOTAL, UNIT

      ENDDO

      !==============================================================
      ! RCP Aircraft emissions
      !==============================================================

      WRITE( 6, '(a)' )
      DO K=1, SIZE(IDTRCP_AIR)
         
         IF (IDTRCP_AIR(K) < 1) CYCLE

         !==============================================================
         ! Global total emission
         !==============================================================

         TOTAL = 0d0
         
         ! Loop over latitudes
         DO J = 1, JJPAR

            ! Surface area [cm2] * seconds in the month / Avogadro's number
            ! Also multiply by the factor 1d-9 to convert kg to Tg
            A = GET_AREA_CM2( 1, J, 1 ) * ( D(THISMONTH) * 86400d-9 ) 
     &           / 6.0225d23
            
            ! Anthro emissions
            TOTAL = TOTAL + SUM(RCP_AIR(:,J,:,K)) * A * 
     &           TRACER_MW_KG(IDTRCP_AIR(K)) 
            
         ENDDO

         !==============================================================
         !  Units
         !==============================================================

         SELECT CASE ( TRACER_NAME(IDTRCP_AIR(K)) )
         CASE ( 'NOx' )
            ! Convert to Tg(N)
            TOTAL = TOTAL * 0.014 / TRACER_MW_KG(IDTRCP_AIR(K))
            UNIT='N'
         CASE ( 'SO2' ) 
            ! Convert to Tg(S)
            TOTAL = TOTAL * 0.032 / TRACER_MW_KG(IDTRCP_AIR(K))
            UNIT='S'
         CASE ( 'NH3' ) 
            ! Convert to Tg(N)
            TOTAL = TOTAL * 0.014 / TRACER_MW_KG(IDTRCP_AIR(K))
            UNIT='N'
         CASE ( 'CO'  ) 
            UNIT='CO'
         CASE DEFAULT
            UNIT='C'
         END SELECT

         !==============================================================
         !  Print info
         !==============================================================

         WRITE( 6, 101 ) 'Air', TRACER_NAME(IDTRCP_AIR(K)), THISMONTH, 
     &           TOTAL, UNIT

      ENDDO

      END SUBROUTINE TOTAL_ANTHRO_RCP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rcp_emission
!
! !DESCRIPTION: Function GET\_RCP\_EMISSION retrieves the emissions of tracer N
! at grid location (I,J). Use LAND=.TRUE. or SHIP=.TRUE. or both to retrieve
!  either land anthropogenic emissions, ship emissions, or their sum. 
!  "N" is the advected tracer index, i.e. the tracer index for STT.
!  The function will return -1 if no emissions are found for that species.
!  (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RCP_EMISSION( I, J, N, LAND, SHIP ) 
     &         RESULT( EMISS )
!
! !USES:
!
      USE TRACERID_MOD
      USE ERROR_MOD,            ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)           :: I, J
      INTEGER, INTENT(IN)           :: N     !GEOS-Chem advected tracer index
      LOGICAL, INTENT(IN), OPTIONAL :: SHIP
      LOGICAL, INTENT(IN), OPTIONAL :: LAND 
!
! !REVISION HISTORY: 
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                        :: EMISS
      CHARACTER(LEN=20)             :: STR
      LOGICAL                       :: DOLAND, DOSHIP, TRACERFOUND
      INTEGER                       :: K

      !=================================================================
      ! GET_RCP_EMISSION begins here!
      !=================================================================
      
      ! Are we getting land emissions?
      IF ( PRESENT( LAND ) ) THEN
         DOLAND = LAND
      ELSE
         DOLAND = .FALSE.
      ENDIF

      ! Are we getting ship emissions?
      IF ( PRESENT( SHIP ) ) THEN
         DOSHIP = SHIP
      ELSE
         DOSHIP = .FALSE.
      ENDIF

      ! Throw error if neither emission type is requested
      IF ( .NOT. (DOLAND .OR. DOSHIP) ) THEN
         WRITE( STR, '(I4)' ) N
         CALL ERROR_STOP( 'No land/ship emissions, tracer '//trim(STR),
     &        'GET_RCP_EMISSION' )
      ENDIF
      
      ! Initialize
      EMISS = 0d0
      TRACERFOUND = .FALSE.

      ! Find tracer number for land emissions
      IF ( DOLAND ) THEN
         ! Loop over all the species we have land emissions for
         DO K=1, SIZE(IDTRCP_LAND)
            IF (N == IDTRCP_LAND(K)) THEN
               ! We found the desired tracer, so add it up and exit loop
               EMISS = EMISS + RCP_LAND(I,J,K)
               TRACERFOUND=.TRUE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! Find tracer number for ship emissions
      IF ( DOSHIP ) THEN
         ! Loop over all the species we have ship emissions for
         DO K=1, SIZE(IDTRCP_SHIP)
            IF (N == IDTRCP_SHIP(K)) THEN
               ! We found the desired tracer, so add it up and exit loop
               EMISS = EMISS + RCP_SHIP(I,J,K)
               TRACERFOUND=.TRUE.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! Return -1 if there are no emissions for tracer N
      IF (.NOT. TRACERFOUND) EMISS = -1d0
      
      END FUNCTION GET_RCP_EMISSION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rcp
!
! !DESCRIPTION: Subroutine INIT\_RCP allocates and zeroes all module arrays
!  (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_RCP( Input_Opt )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      USE CMN_SIZE_MOD             ! Size parameters

!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY:
!  22 Jul 2013 - M. Sulprizio- Now copy LRCP, LRCPSHIP, LRCPAIR from Input_Opt
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER    :: AS

      ! For fields from Input_Opt
      LOGICAL    :: LRCP
      LOGICAL    :: LRCPSHIP
      LOGICAL    :: LRCPAIR

      !=================================================================
      !  INIT_RCP begins here
      !=================================================================

      ! Copy values from Input_Opt
      LRCP        = Input_Opt%LRCP
      LRCPSHIP    = Input_Opt%LRCPSHIP
      LRCPAIR     = Input_Opt%LRCPAIR

      ! Return if we LRCP = .FALSE.
      IF ( .not. (LRCP .OR. LRCPSHIP .OR. LRCPAIR) ) RETURN

      IDTRCP_LAND = 0d0
      IDTRCP_SHIP = 0d0
      IDTRCP_AIR  = 0d0

      ! Anthropogenic land surface emissions
      ALLOCATE( RCP_LAND( IIPAR, JJPAR, SIZE(IDTRCP_LAND) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_LAND' )
      RCP_LAND = 0e0

      ! Shipping 
      ALLOCATE( RCP_SHIP( IIPAR, JJPAR, SIZE(IDTRCP_SHIP) ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_SHIP' )
      RCP_SHIP = 0e0

      ! Aircraft
      ALLOCATE( RCP_AIR( IIPAR, JJPAR, LLPAR, SIZE(IDTRCP_AIR) ),
     &     STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RCP_AIR' )
      RCP_AIR = 0e0

      ! Array to pass aircraft NOx to SMVGEAR
      ALLOCATE( RCP_AC_NOx( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIREMS' )
      RCP_AC_NOx = 0d0

      END SUBROUTINE INIT_RCP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_rcp
!
! !DESCRIPTION: Subroutine CLEANUP\_RCP deallocates all module arrays
!  (cdh, 10/14/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_RCP
!
! !REVISION HISTORY: 
!  22 Jul 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      !  CLEANUP_RCP begins here
      !=================================================================

      IF ( ALLOCATED( RCP_LAND   ) ) DEALLOCATE( RCP_LAND   )
      IF ( ALLOCATED( RCP_SHIP   ) ) DEALLOCATE( RCP_SHIP   )
      IF ( ALLOCATED( RCP_AIR    ) ) DEALLOCATE( RCP_AIR    )
      IF ( ALLOCATED( RCP_AC_NOx ) ) DEALLOCATE( RCP_AC_NOx )

      END SUBROUTINE CLEANUP_RCP
!EOC
      END MODULE RCP_MOD
