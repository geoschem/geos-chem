!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: retro_mod 
!
! !DESCRIPTION: Module RETRO\_MOD reads emissions from the RETRO emissions 
!  inventory
!\\
!\\
! !INTERFACE:
!
      MODULE RETRO_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!      
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
!
! !PUBLIC MEMBER FUNCTIONS:
!  
      PUBLIC  :: CLEANUP_RETRO
      PUBLIC  :: EMISS_RETRO
      PUBLIC  :: GET_RETRO_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_RETRO
      PRIVATE :: READ_RETRO
      PRIVATE :: TOTAL_ANTHRO_Tg
!
! !REVISION HISTORY:
!  08 Mar 2011 - W. Reinhart - Initial version
!  18 Aug 2011 - D. Millet   - Partition ketones into 25% MEK and 75% ACET
!  18 Aug 2011 - D. Millet   - Remove call to GET_ANNUAL_SCALAR
!  22 Aug 2011 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  22 Mar 2012 - M. Payer    - RETRO C2H6 emissions are too low. Use
!                              Yaping Xiao's C2H6 emissions instead.
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS

!----------------------------------------------------------------------- 

      SUBROUTINE EMISS_RETRO( am_I_Root, Input_Opt, State_Chm, RC )
!
!***********************************************************************
!  Subroutine EMISS_RETRO reads all RETRO emissions at the beginning of
!  each month. (wfr, 3/8/11)
!***********************************************************************

!
! !USES:
!
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,            ONLY : GET_RES_EXT
      USE FILE_MOD,             ONLY : IOERROR
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C3H8ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GIGC_State_Chm_Mod,   ONLY : ChmState
      USE TIME_MOD,             ONLY : EXPAND_DATE
      USE TIME_MOD,             ONLY : GET_MONTH
      USE CMN_SIZE_MOD

#if defined( DEVEL )
      USE TRACERID_MOD,         ONLY : IDTALK4, IDTACET, IDTMEK,
     & IDTALD2, IDTPRPE, IDTC3H8, IDTC2H6, IDTCH2O, IDTBENZ,
     & IDTTOLU, IDTXYLE, IDTC2H4, IDTC2H2
      USE GRID_MOD,             ONLY : GET_AREA_CM2
      USE ERROR_MOD,            ONLY : ALLOC_ERR
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
!  08 Mar 2011 - W. Reinhart - Initial version
!  22 Aug 2011 - R. Yantosca - Added ProTeX header
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      LOGICAL, SAVE      :: FIRST = .TRUE.
      INTEGER            :: I,      J,      THISMONTH, YYYYMMDD
      REAL*8             :: ALK4ff, PRPEff, C3H8ff
      REAL*8             :: C2H6ff, VOCff
      CHARACTER(LEN=255) :: FILENAME

#if defined( DEVEL )
      REAL*8, ALLOCATABLE   :: A(:,:)
      INTEGER AS
#endif

      ! For fields from Input_Opt
      LOGICAL         :: LFUTURE
      INTEGER         :: N_TRACERS
      REAL*8          :: TRACER_MW_KG(Input_Opt%N_TRACERS)

      !=================================================================
      !  EMISS_RETRO begins here
      !=================================================================

      ! Assume success
      RC           =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFUTURE      = Input_Opt%LFUTURE
      N_TRACERS    = Input_Opt%N_TRACERS
      TRACER_MW_KG = Input_Opt%TRACER_MW_KG(1:N_TRACERS)

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays
#if defined( DEVEL )
         ALLOCATE( A( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMISS_EPA_NEI:A' )
         A = 0d0
#endif         

         CALL INIT_RETRO( am_I_Root, Input_Opt, RC )

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

#if defined( DEVEL )
      write(6,*)' ##### DEVEL BLOCK'
      call flush(6)
      DO I=1,IIPAR
      DO J=1,JJPAR
         A(I,J)  = GET_AREA_CM2( I, J, 1 ) 
      ENDDO
      ENDDO


      IF ( IDTALK4 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTALK4) = RETRO_ALK4(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTALK4)

      IF ( IDTACET > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTACET) = RETRO_ACET(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTACET)

      IF ( IDTMEK  > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTMEK)  = RETRO_MEK(:,:)  * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTMEK)

      IF ( IDTALD2 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTALD2) = RETRO_ALD2(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTALD2)

      IF ( IDTPRPE > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTPRPE) = RETRO_PRPE(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTPRPE)

      IF ( IDTC3H8 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTC3H8) = RETRO_C3H8(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTC3H8)

      IF ( IDTC2H6 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTC2H6) = RETRO_C2H6(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTC2H6)

      IF ( IDTCH2O > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTCH2O) = RETRO_CH2O(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTCH2O)

      IF ( IDTBENZ > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTBENZ) = RETRO_BENZ(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTBENZ)

      IF ( IDTTOLU > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTTOLU) = RETRO_TOLU(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTTOLU)

      IF ( IDTXYLE > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTXYLE) = RETRO_XYLE(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTXYLE)

      IF ( IDTC2H4 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTC2H4) = RETRO_C2H4(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTC2H4)

      IF ( IDTC2H2 > 0 ) 
     & State_Chm%TRAC_TEND(:,:,1,IDTC2H2) = RETRO_C2H2(:,:) * 
     & A * 6.0225d-23 * TRACER_MW_KG(IDTC2H2)

#endif

      ! Return to calling program
      END SUBROUTINE EMISS_RETRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_retro
!
! !DESCRIPTION: Subroutine READ\_RETRO reads a BPCH file created from 
!  RETRO data.  The data has units of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_RETRO( FILENAME, ALK4, ACET, MEK,  ALD2, PRPE,
     &                       C3H8,     C2H6, CH2O, BENZ, TOLU, XYLE, 
     &                       C2H4,     C2H2                          )
!
! !USES:
!
      USE BPCH2_MOD,        ONLY : OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IOERROR
      USE inquireMod,       ONLY : findFreeLUN
      USE TRANSFER_MOD,     ONLY : TRANSFER_2D
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR
      USE TIME_MOD,         ONLY : GET_YEAR
      USE CMN_SIZE_MOD           ! Size parameters
      USE CMN_O3_MOD             ! FSCLYR
!
! !INPUT PARAMETERS: 
!
      ! Name of file to read
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME 
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! RETRO emissions for various VOC species [molec/cm2/s]
      REAL*4,           INTENT(INOUT) :: ALK4(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: ACET(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: MEK (IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: ALD2(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: PRPE(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: C3H8(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: CH2O(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H6(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: BENZ(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: TOLU(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: XYLE(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H4(IIPAR,JJPAR)
      REAL*4,           INTENT(INOUT) :: C2H2(IIPAR,JJPAR)
! 
! !REVISION HISTORY: 
!  08 Mar 2011 - W. Reinhart - Initial Version 
!  18 Aug 2011 - D. Millet   - Remove call to GET_ANNUAL_SCALAR
!  22 Aug 2011 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, IU_FILE, J, L, N, IOS
      INTEGER           :: NI, NJ, NL
      INTEGER           :: IFIRST, JFIRST, LFIRST
      INTEGER           :: NTRACER, NSKIP
      INTEGER           :: HALFPOLAR, CENTER180
      INTEGER           :: SCALEYEAR  !, BASEYEAR (dbm, 8/18/11)
      REAL*4            :: LONRES, LATRES
      REAL*4            :: ARRAY(IIPAR,JJPAR,1)
      REAL*4            :: SC(IIPAR,JJPAR)
      REAL*8            :: ZTAU0, ZTAU1
      CHARACTER(LEN=20) :: MODELNAME
      CHARACTER(LEN=40) :: CATEGORY
      CHARACTER(LEN=40) :: UNIT
      CHARACTER(LEN=40) :: RESERVED

      !=================================================================
      !  READ_RETRO begins here
      !=================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME ), IU_FILE
 100  FORMAT( 'READ_RETRO: Reading ', a, ' on unit ', i4 )

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

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

      END SUBROUTINE READ_RETRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_ANTHRO_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_Tg to print total RETRO 
!  anthropogenic VOC emissions each month in [Tg C].
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_Tg( THISMONTH )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : TRACER_MW_KG
      USE TRACERID_MOD, ONLY : IDTALK4, IDTMEK,  IDTPRPE, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTC2H6, IDTBENZ, IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE, IDTC2H4, IDTC2H2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALD2
      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Current month
! 
! !REVISION HISTORY: 
!  08 Mar 2011 - W. Reinhart - Initial Version 
!  22 Aug 2011 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  22 Mar 2012 - M. Payer    - Remove print for C2H6 emissions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J
      REAL*8           :: ALK4, MEK,  ALD2, PRPE, C3H8, CH2O
      REAL*8           :: BENZ, TOLU, XYLE, C2H4, C2H2, C2H6, ACET
      REAL*8           :: F_ALK4, F_MEK,  F_PRPE, F_C3H8, F_CH2O
      REAL*8           :: F_BENZ, F_TOLU, F_XYLE, F_C2H4, F_C2H2
      REAL*8           :: F_C2H6, F_ALD2, F_ACET
      REAL*8           :: A
      CHARACTER(LEN=6) :: UNIT

      ! Days per month
      INTEGER          :: D(12) = (/ 31, 28, 31, 30, 31, 30,
     &                               31, 31, 30, 31, 30, 31 /)

      !=================================================================
      ! TOTAL_ANTHRO_Tg begins here
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

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Surface area [cm2] * seconds in the month / Avogadro's number
         ! Also multiply by the factor 1d-9 to convert kg to Tg
         A     = GET_AREA_CM2( I, J, 1 ) 
     &         * ( D(THISMONTH) * 86400d-9 ) / 6.0225d23

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
      WRITE( 6, '(a)' )
      WRITE( 6, 100   ) 'ALK4', THISMONTH, ALK4, ' C'
      WRITE( 6, 100   ) 'ACET', THISMONTH, ACET, ' C'
      WRITE( 6, 100   ) 'MEK',  THISMONTH, MEK,  ' C'
      WRITE( 6, 100   ) 'ALD2', THISMONTH, ALD2, ' C'
      WRITE( 6, 100   ) 'PRPE', THISMONTH, PRPE, ' C'
      WRITE( 6, 100   ) 'C3H8', THISMONTH, C3H8, ' C'
      WRITE( 6, 100   ) 'CH2O', THISMONTH, CH2O, ' C'
      WRITE( 6, 100   ) 'BENZ', THISMONTH, BENZ, ' C'
      WRITE( 6, 100   ) 'TOLU', THISMONTH, TOLU, ' C'
      WRITE( 6, 100   ) 'XYLE', THISMONTH, XYLE, ' C'
      WRITE( 6, 100   ) 'C2H4', THISMONTH, C2H4, ' C'
      WRITE( 6, 100   ) 'C2H2', THISMONTH, C2H2, ' C'
 100  FORMAT( 'Total anthro ', a4, ' for 2000/',
     &         i2.2, ': ', f13.6, ' Tg', a2 )

      WRITE( 6, '(/,a,/)' ) 'RETRO_MOD: RETRO C2H6 anthro ' //
     &   'emissions are too low. Using offline C2H6 '       //
     &   'emissions from Yaping Xiao.'

      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_retro_anthro
!
! !DESCRIPTION: Function GET\_RETRO\_ANTHRO returns the monthly average
!  anthropogenic VOC emissions at GEOS-Chem grid box (I,J).  Data will
!  be returned in units of [atoms C/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RETRO_ANTHRO( I, J, N ) RESULT( RETRO )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTALK4, IDTMEK,  IDTPRPE, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTC2H6, IDTBENZ, IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE, IDTC2H4, IDTC2H2
      USE TRACERID_MOD, ONLY : IDTACET, IDTALD2
      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I   ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J   ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: N   ! GEOS-Chem tracer index
!
! !RETURN VALUE:
!
      REAL*8              :: RETRO   ! RETRO emissions [mole
! 
! !REVISION HISTORY: 
!  08 Mar 2011 - W. Reinhart - Initial Version 
!  18 Aug 2011 - D. Millet   - Partition RETRO ketones into 75% acetone 
!                              and 25% MEK
!  22 Mar 2012 - M. Payer    - RETRO C2H6 emissions are too low. Use
!                              Yaping Xiao's C2H6 emissions instead.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!

      !=================================================================
      !  GET_RETRO_ANTHRO begins here
      !=================================================================

      IF ( N == IDTALK4 )      THEN
         RETRO = RETRO_ALK4(I,J)
!------------------------------------------------------------------------------
!      ! Partition RETRO ketones into 75% acetone and 25% MEK (dbm, 8/18/2011)
!      ELSE IF ( N == IDTACET ) THEN
!         RETRO = RETRO_ACET(I,J)
!      ELSE IF ( N == IDTMEK )  THEN
!         RETRO = RETRO_MEK(I,J)
!------------------------------------------------------------------------------
      ELSE IF ( N == IDTACET ) THEN
         RETRO = 0.75d0*RETRO_MEK(I,J)   ! RETRO ketones --> 75% ACET
      ELSE IF ( N == IDTMEK )  THEN
         RETRO = 0.25d0*RETRO_MEK(I,J)   ! RETRO ketones --> 25% MEK 
      ELSE IF ( N == IDTALD2 ) THEN
         RETRO = RETRO_ALD2(I,J)
      ELSE IF ( N == IDTPRPE ) THEN
         RETRO = RETRO_PRPE(I,J)
      ELSE IF ( N == IDTC3H8 ) THEN
         RETRO = RETRO_C3H8(I,J)
      ELSE IF ( N == IDTCH2O ) THEN
         RETRO = RETRO_CH2O(I,J)
      ELSE IF ( N == IDTC2H6 ) THEN
         RETRO = -1d0
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

      END FUNCTION GET_RETRO_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_retro
!
! !DESCRIPTION: Subroutine INIT\_RETRO allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_RETRO( am_I_Root, Input_Opt, RC )
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
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  08 Mar 2011 - W. Reinhart - Initial Version 
!  22 Aug 2011 - R. Yantosca - Added ProTeX headers

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: AS

      !=================================================================
      !  INIT_RETRO begins here
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS
      
      ! Return if we LRETRO = .FALSE.
      IF (.not. Input_Opt%LRETRO ) RETURN

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

      END SUBROUTINE INIT_RETRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_retro
!
! !DESCRIPTION: Subroutine CLEANUP\_RETRO deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_RETRO
! 
! !REVISION HISTORY: 
!  08 Mar 2011 - W. Reinhart - Initial Version 
!  22 Aug 2011 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      !  CLEANUP_RETRO begins here
      !=================================================================
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

      END SUBROUTINE CLEANUP_RETRO
!EOC
      END MODULE RETRO_MOD
