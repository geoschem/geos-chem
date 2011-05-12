!-----------------------------------------------------------------------
!          University of Minnesota Atmospheric Chemistry Group
!-----------------------------------------------------------------------

      MODULE RETRO_MOD

!***********************************************************************
!  Module RETRO_MOD reads emissions from the RETRO emissions inventory
!  (wfr 3/8/11)
!***********************************************************************

      IMPLICIT NONE

      !Make everything PRIVATE...
      PRIVATE

      !...except these routines
      PUBLIC :: CLEANUP_RETRO
      PUBLIC :: EMISS_RETRO
      PUBLIC :: GET_RETRO_ANTHRO
      
      REAL*4, ALLOCATABLE :: RETRO_ALK4(:,:)
      REAL*4, ALLOCATABLE :: RETRO_ACET(:,:)
      REAL*4, ALLOCATABLE :: RETRO_MEK(:,:)
      REAL*4, ALLOCATABLE :: RETRO_ALD2(:,:)
      REAL*4, ALLOCATABLE :: RETRO_PRPE(:,:)
      REAL*4, ALLOCATABLE :: RETRO_C3H8(:,:)
      REAL*4, ALLOCATABLE :: RETRO_C2H6(:,:)
      REAL*4, ALLOCATABLE :: RETRO_CH2O(:,:)
      REAL*4, ALLOCATABLE :: RETRO_BENZ(:,:)
      REAL*4, ALLOCATABLE :: RETRO_TOLU(:,:)
      REAL*4, ALLOCATABLE :: RETRO_XYLE(:,:)
      REAL*4, ALLOCATABLE :: RETRO_C2H4(:,:)
      REAL*4, ALLOCATABLE :: RETRO_C2H2(:,:)

      !========================================
      !  MODULE ROUTINES -- follow "CONTAINS"
      !========================================
      CONTAINS

!-----------------------------------------------------------------------
      SUBROUTINE EMISS_RETRO
!***********************************************************************
!  Subroutine EMISS_RETRO reads all RETRO emissions at the beginning of
!  each month. (wfr, 3/8/11)
!***********************************************************************

      !References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE FILE_MOD,             ONLY : IU_FILE, IOERROR
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C3H8ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff
      USE LOGICAL_MOD,          ONLY : LFUTURE
      USE TIME_MOD,             ONLY : EXPAND_DATE, GET_MONTH

#     include "CMN_SIZE"             ! Size parameters
      
      !Local variables
      LOGICAL, SAVE                 :: FIRST = .TRUE.
      INTEGER                       :: I, J, THISMONTH, YYYYMMDD
      REAL*8                        :: ALK4ff, PRPEff, C3H8ff
      REAL*8                        :: C2H6ff, VOCff
      
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      !  EMISS_RETRO begins here
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         print*,'### RETRO: initializing'

         ! Allocate arrays
         CALL INIT_RETRO

         ! Reset first-time flag
         FIRST = .FALSE.

      ENDIF

      ! Get month
      THISMONTH = GET_MONTH()

      ! Get date for 2000 emissions
      YYYYMMDD = 20000000 + ( THISMONTH * 100 ) + 01

      ! Echo info
      WRITE(6, '(a)' ) REPEAT( '=', 79)
      WRITE(6, 100   )
 100  FORMAT( 'R E T R O    E M I S S I O N S',
     &        '  --  Baseline Year: 2000', / )

      !=================================================================
      !  Read RETRO average annual anthropogenic emissions
      !=================================================================

      ! Anthro file name
      FILENAME = TRIM( DATA_DIR ) // 'RETRO_201103/' //
     &                  'YYYYMM.' // GET_RES_EXT()

      ! Replace date in filename
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Read data
      CALL READ_RETRO( FILENAME,   RETRO_ALK4, RETRO_ACET, RETRO_MEK,
     &                 RETRO_ALD2, RETRO_PRPE, RETRO_C3H8, RETRO_C2H6,
     &                 RETRO_CH2O, RETRO_BENZ, RETRO_TOLU, RETRO_XYLE,
     &                 RETRO_C2H4, RETRO_C2H2                         )

      DO J = 1, JJPAR
      DO I = 1, IIPAR


         !-----------------------------------
         !  Calculate IPCC future emissions
         !-----------------------------------
         IF ( LFUTURE ) THEN

            ! Future anthro scale factors
            ALK4ff  = GET_FUTURE_SCALE_ALK4ff( I, J )
            VOCff   = GET_FUTURE_SCALE_VOCff(  I, J )
            PRPEff  = GET_FUTURE_SCALE_PRPEff( I, J )
            C3H8ff  = GET_FUTURE_SCALE_C3H8ff( I, J )
            C2H6ff  = GET_FUTURE_SCALE_C2H6ff( I, J )

            ! Apply scale factors
            RETRO_ALK4 (I,J) = RETRO_ALK4 (I,J) * ALK4ff
            RETRO_ACET (I,J) = RETRO_ACET (I,J) * VOCff
            RETRO_MEK  (I,J) = RETRO_MEK  (I,J) * VOCff
            RETRO_ALD2 (I,J) = RETRO_ALD2 (I,J) * VOCff
            RETRO_PRPE (I,J) = RETRO_PRPE (I,J) * PRPEff
            RETRO_C3H8 (I,J) = RETRO_C3H8 (I,J) * C3H8ff
            RETRO_C2H6 (I,J) = RETRO_C2H6 (I,J) * C2H6ff
            RETRO_CH2O (I,J) = RETRO_CH2O (I,J) * VOCff
            RETRO_BENZ (I,J) = RETRO_BENZ (I,J) * VOCff
            RETRO_TOLU (I,J) = RETRO_TOLU (I,J) * VOCff
            RETRO_XYLE (I,J) = RETRO_XYLE (I,J) * VOCff
            RETRO_C2H4 (I,J) = RETRO_C2H4 (I,J) * VOCff
            RETRO_C2H2 (I,J) = RETRO_C2H2 (I,J) * VOCff
         ENDIF
      ENDDO
      ENDDO

      ! Print totals to log
      CALL TOTAL_ANTHRO_TG( THISMONTH )

      ! Fancy output
      WRITE(6, '(a)' ) REPEAT( '=', 79)

      ! Return to calling program
      END SUBROUTINE EMISS_RETRO

!-----------------------------------------------------------------------

      SUBROUTINE READ_RETRO( FILENAME, ALK4, ACET, MEK, ALD2, PRPE,
     &                       C3H8, C2H6, CH2O, BENZ, TOLU, XYLE, C2H4,
     &                       C2H2 )                       

!***********************************************************************
!  Subroutine READ_RETRO reads a BPCH file created from RETRO data.
!***********************************************************************

      ! Reference to F90 modules
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      USE TIME_MOD,         ONLY : GET_YEAR

#     include "CMN_SIZE"             ! Size parameters
#     include "CMN_O3"               ! FSCLYR

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      REAL*4,           INTENT(INOUT) :: ALK4(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: ACET(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: MEK(IIPAR,  JJPAR)
      REAL*4,           INTENT(INOUT) :: ALD2(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: PRPE(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: C3H8(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: CH2O(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H6(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: BENZ(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: TOLU(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: XYLE(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H4(IIPAR, JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H2(IIPAR, JJPAR)

      ! Local variables
      INTEGER                         :: I, J, L, N, IOS
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
      !  READ_RETRO begins here
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_RETRO: Reading ', a )

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME)

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
         SELECT CASE ( NTRACER )
            CASE( 5  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ALK4 )
            CASE( 9  )
               CALL TRANSFER_2D( ARRAY(:,:,1), ACET )
            CASE( 10 )
               CALL TRANSFER_2D( ARRAY(:,:,1), MEK  )
            CASE( 11 )
               CALL TRANSFER_2D( ARRAY(:,:,1), ALD2 )
            CASE( 18 )
               CALL TRANSFER_2D( ARRAY(:,:,1), PRPE )
            CASE( 19 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C3H8 )
            CASE( 20 )
               CALL TRANSFER_2D( ARRAY(:,:,1), CH2O )
            CASE( 21 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C2H6 )
            CASE( 59 )
               CALL TRANSFER_2D( ARRAY(:,:,1), BENZ )
            CASE( 60 )
               CALL TRANSFER_2D( ARRAY(:,:,1), TOLU )
            CASE( 61 )
               CALL TRANSFER_2D( ARRAY(:,:,1), XYLE )
            CASE( 65 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C2H4 )
            CASE( 66 )
               CALL TRANSFER_2D( ARRAY(:,:,1), C2H2 )
            CASE DEFAULT
               ! Nothing
         END SELECT
      END DO

      ! Close file
      CLOSE( IU_FILE )

      ! Apply annual scalar factor
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      BASEYEAR=2000

      CALL GET_ANNUAL_SCALAR( 71, BASEYEAR, SCALEYEAR, SC)
      ALK4(:,:) = ALK4(:,:) * SC(:,:)

      ! Return to calling program
      END SUBROUTINE READ_RETRO

!-----------------------------------------------------------------------

      SUBROUTINE TOTAL_ANTHRO_TG( THISMONTH )

!***********************************************************************
!  Subroutine to print total RETRO anthropogenic emissions each month
!  in TgC
!***********************************************************************

      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTALK4, IDTMEK,  IDTPRPE, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTC2H6, IDTBENZ, IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE, IDTC2H4, IDTC2H2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALD2

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH

      ! Local variables
      INTEGER             :: I, J
      REAL*8              :: ALK4, MEK,  ALD2, PRPE, C3H8, CH2O
      REAL*8              :: BENZ, TOLU, XYLE, C2H4, C2H2, C2H6, ACET
      REAL*8              :: F_ALK4, F_MEK,  F_PRPE, F_C3H8, F_CH2O
      REAL*8              :: F_BENZ, F_TOLU, F_XYLE, F_C2H4, F_C2H2
      REAL*8              :: F_C2H6, F_ALD2, F_ACET
      REAL*8              :: A
      CHARACTER(LEN=6)    :: UNIT

      ! Days per month
      INTEGER             :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                                  31, 31, 30, 31, 30, 31 /)

      !=================================================================
      !  TOTAL_ANTHRO_TG begins here
      !=================================================================

      ! Summing variables for anthro
      ALK4 = 0d0
      ACET = 0d0
      MEK  = 0d0
      ALD2 = 0d0
      PRPE = 0d0
      C3H8 = 0d0
      CH2O = 0d0
      C2H6 = 0d0
      BENZ = 0d0
      TOLU = 0d0
      XYLE = 0d0
      C2H4 = 0d0
      C2H2 = 0d0

      ! Molecular weights
      F_ALK4 = 0d0
      F_ACET = 0d0
      F_MEK  = 0d0
      F_ALD2 = 0d0
      F_PRPE = 0d0
      F_C3H8 = 0d0
      F_CH2O = 0d0
      F_C2H6 = 0d0
      F_BENZ = 0d0
      F_TOLU = 0d0
      F_XYLE = 0d0
      F_C2H4 = 0d0
      F_C2H2 = 0d0

      ! Prevent array out of bounds error for undefined tracers
      IF ( IDTALK4 > 0 ) F_ALK4 = TRACER_MW_KG(IDTALK4)
      IF ( IDTACET > 0 ) F_ACET = TRACER_MW_KG(IDTACET)
      IF ( IDTMEK  > 0 ) F_MEK  = TRACER_MW_KG(IDTMEK )
      IF ( IDTALD2 > 0 ) F_ALD2 = TRACER_MW_KG(IDTALD2)
      IF ( IDTPRPE > 0 ) F_PRPE = TRACER_MW_KG(IDTPRPE)
      IF ( IDTC2H6 > 0 ) F_C2H6 = TRACER_MW_KG(IDTC2H6)
      IF ( IDTC3H8 > 0 ) F_C3H8 = TRACER_MW_KG(IDTC3H8)
      IF ( IDTCH2O > 0 ) F_CH2O = TRACER_MW_KG(IDTCH2O)
      IF ( IDTBENZ > 0 ) F_BENZ = TRACER_MW_KG(IDTBENZ)
      IF ( IDTTOLU > 0 ) F_TOLU = TRACER_MW_KG(IDTTOLU)
      IF ( IDTXYLE > 0 ) F_XYLE = TRACER_MW_KG(IDTXYLE)
      IF ( IDTC2H4 > 0 ) F_C2H4 = TRACER_MW_KG(IDTC2H4)
      IF ( IDTC2H2 > 0 ) F_C2H2 = TRACER_MW_KG(IDTC2H2)

      !=================================================================
      !  Sum anthropogenic emissions
      !=================================================================
      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Surface area [cm2] * seconds in the month / Avogadro's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A = GET_AREA_CM2( J ) * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23
      
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Anthro emissions
            ALK4  = ALK4  + RETRO_ALK4(I,J)  * A * F_ALK4
            ACET  = ACET  + RETRO_ACET(I,J)  * A * F_ACET
            MEK   = MEK   + RETRO_MEK(I,J)   * A * F_MEK
            ALD2  = ALD2  + RETRO_ALD2(I,J)  * A * F_ALD2
            PRPE  = PRPE  + RETRO_PRPE(I,J)  * A * F_PRPE
            C3H8  = C3H8  + RETRO_C3H8(I,J)  * A * F_C3H8
            CH2O  = CH2O  + RETRO_CH2O(I,J)  * A * F_CH2O
            C2H6  = C2H6  + RETRO_C2H6(I,J)  * A * F_C2H6
            BENZ  = BENZ  + RETRO_BENZ(I,J)  * A * F_BENZ
            TOLU  = TOLU  + RETRO_TOLU(I,J)  * A * F_TOLU
            XYLE  = XYLE  + RETRO_XYLE(I,J)  * A * F_XYLE
            C2H4  = C2H4  + RETRO_C2H4(I,J)  * A * F_C2H4
            C2H2  = C2H2  + RETRO_C2H2(I,J)  * A * F_C2H2

            ENDDO
         ENDDO

         !==============================================================
         !  Print info
         !==============================================================

         ! Anthro
         WRITE( 6, '(a)' )
         WRITE( 6, 100   ) 'ALK4', THISMONTH, ALK4, ' C'
         WRITE( 6, 100   ) 'ACET', THISMONTH, ACET, ' C'
         WRITE( 6, 100   ) 'MEK',  THISMONTH, MEK,  ' C'
         WRITE( 6, 100   ) 'ALD2', THISMONTH, ALD2, ' C'
         WRITE( 6, 100   ) 'PRPE', THISMONTH, PRPE, ' C'
         WRITE( 6, 100   ) 'C3H8', THISMONTH, C3H8, ' C'
         WRITE( 6, 100   ) 'CH2O', THISMONTH, CH2O, ' C'
         WRITE( 6, 100   ) 'C2H6', THISMONTH, C2H6, ' C'
         WRITE( 6, 100   ) 'BENZ', THISMONTH, BENZ, ' C'
         WRITE( 6, 100   ) 'TOLU', THISMONTH, TOLU, ' C'
         WRITE( 6, 100   ) 'XYLE', THISMONTH, XYLE, ' C'
         WRITE( 6, 100   ) 'C2H4', THISMONTH, C2H4, ' C'
         WRITE( 6, 100   ) 'C2H2', THISMONTH, C2H2, ' C'
 100     FORMAT( 'Total anthro ', a4, ' for 2000/',
     &           i2.2, ': ', f13.6, ' Tg', a2 )

         ! Return to calling program
         END SUBROUTINE TOTAL_ANTHRO_TG

!-----------------------------------------------------------------------

      FUNCTION GET_RETRO_ANTHRO( I, J, N ) RESULT( RETRO )

!***********************************************************************
!  Returns the monthly avg anthropogenic emissions at (I,J)
!***********************************************************************

      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDTALK4, IDTMEK,  IDTPRPE, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTC2H6, IDTBENZ, IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE, IDTC2H4, IDTC2H2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALD2

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, N

      ! Local variables
      REAL*8              :: RETRO

      !=================================================================
      !  GET_RETRO_ANTHRO begins here
      !=================================================================

      IF ( N == IDTALK4 )      THEN
         RETRO = RETRO_ALK4(I,J)
      ELSE IF ( N == IDTACET ) THEN
         RETRO = RETRO_ACET(I,J)
      ELSE IF ( N == IDTMEK )  THEN
         RETRO = RETRO_MEK(I,J)
      ELSE IF ( N == IDTALD2 ) THEN
         RETRO = RETRO_ALD2(I,J)
      ELSE IF ( N == IDTPRPE ) THEN
         RETRO = RETRO_PRPE(I,J)
      ELSE IF ( N == IDTC3H8 ) THEN
         RETRO = RETRO_C3H8(I,J)
      ELSE IF ( N == IDTCH2O ) THEN
         RETRO = RETRO_CH2O(I,J)
      ELSE IF ( N == IDTC2H6 ) THEN
         RETRO = RETRO_C2H6(I,J)
      ELSE IF ( N == IDTBENZ ) THEN
         RETRO = RETRO_BENZ(I,J)
      ELSE IF ( N == IDTTOLU ) THEN
         RETRO = RETRO_TOLU(I,J)
      ELSE IF ( N == IDTXYLE ) THEN
         RETRO = RETRO_XYLE(I,J)
      ELSE IF ( N == IDTC2H4 ) THEN
         RETRO = RETRO_C2H4(I,J)
      ELSE IF ( N == IDTC2H2 ) THEN
         RETRO = RETRO_C2H2(I,J)
      ELSE
         RETRO = -1d0
      ENDIF

      ! Return to calling program
      END FUNCTION GET_RETRO_ANTHRO

!----------------------------------------------------------------------

      SUBROUTINE INIT_RETRO

!***********************************************************************
!  Subroutine allocates and zeroes all module arrays
!***********************************************************************

      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LRETRO

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: AS

      !=================================================================
      !  INIT_RETRO begins here
      !=================================================================

      ! Return if we LRETRO = .FALSE.
      IF (.not. LRETRO ) RETURN

      ALLOCATE( RETRO_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_ALK4' )
      RETRO_ALK4 = 0e0

      ALLOCATE( RETRO_ACET( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_ACET' )
      RETRO_ACET = 0e0

      ALLOCATE( RETRO_MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_MEK' )
      RETRO_MEK = 0e0

      ALLOCATE( RETRO_ALD2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_ALD2' )
      RETRO_ALD2 = 0e0

      ALLOCATE( RETRO_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_PRPE' )
      RETRO_PRPE = 0e0

      ALLOCATE( RETRO_C3H8( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_C3H8' )
      RETRO_C3H8 = 0e0

      ALLOCATE( RETRO_CH2O( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_CH2O' )
      RETRO_CH2O = 0e0

      ALLOCATE( RETRO_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_C2H6' )
      RETRO_C2H6 = 0e0

      ALLOCATE( RETRO_BENZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_BENZ' )
      RETRO_BENZ = 0e0

      ALLOCATE( RETRO_TOLU( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_TOLU' )
      RETRO_TOLU = 0e0

      ALLOCATE( RETRO_XYLE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_XYLE' )
      RETRO_XYLE = 0e0

      ALLOCATE( RETRO_C2H4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_C2H4' )
      RETRO_C2H4 = 0e0

      ALLOCATE( RETRO_C2H2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RETRO_C2H2' )
      RETRO_C2H2 = 0e0

      ! Return to calling program
      END SUBROUTINE INIT_RETRO

!----------------------------------------------------------------------

      SUBROUTINE CLEANUP_RETRO
!**********************************************************************
!  Subroutine deallocates all module arrays
!**********************************************************************

      !=================================================================
      !  CLEANUP_RETRO begins here
      !=================================================================

      ! Fossil fuel
      IF ( ALLOCATED( RETRO_ALK4 ) ) DEALLOCATE( RETRO_ALK4 )
      IF ( ALLOCATED( RETRO_ACET ) ) DEALLOCATE( RETRO_ACET )
      IF ( ALLOCATED( RETRO_MEK  ) ) DEALLOCATE( RETRO_MEK  )
      IF ( ALLOCATED( RETRO_ALD2 ) ) DEALLOCATE( RETRO_ALD2 )
      IF ( ALLOCATED( RETRO_PRPE ) ) DEALLOCATE( RETRO_PRPE )
      IF ( ALLOCATED( RETRO_C3H8 ) ) DEALLOCATE( RETRO_C3H8 )
      IF ( ALLOCATED( RETRO_CH2O ) ) DEALLOCATE( RETRO_CH2O )
      IF ( ALLOCATED( RETRO_C2H6 ) ) DEALLOCATE( RETRO_C2H6 )
      IF ( ALLOCATED( RETRO_BENZ ) ) DEALLOCATE( RETRO_BENZ )
      IF ( ALLOCATED( RETRO_TOLU ) ) DEALLOCATE( RETRO_TOLU )
      IF ( ALLOCATED( RETRO_XYLE ) ) DEALLOCATE( RETRO_XYLE )
      IF ( ALLOCATED( RETRO_C2H4 ) ) DEALLOCATE( RETRO_C2H4 )
      IF ( ALLOCATED( RETRO_C2H2 ) ) DEALLOCATE( RETRO_C2H2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_RETRO

!----------------------------------------------------------------------

      END MODULE RETRO_MOD
