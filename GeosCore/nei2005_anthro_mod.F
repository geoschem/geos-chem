!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: nei2005_anthro_mod
!
! !DESCRIPTION: Module NEI2005\_ANTHRO\_MOD contains variables and routines to 
!  read the NEI2005 anthropogenic emissions.   
!\\
!\\
! !INTERFACE: 
!
      MODULE NEI2005_ANTHRO_MOD
! 
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      REAL*8, PUBLIC, ALLOCATABLE :: USA_MASK(:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_NEI2005_ANTHRO
      PUBLIC  :: EMISS_NEI2005_ANTHRO
      PUBLIC  :: EMISS_NEI2005_ANTHRO_05x0666
      PUBLIC  :: GET_NEI2005_ANTHRO
      !--------------------------------------
      ! Leave for future use (bmy, 12/3/09)
      !PUBLIC  :: GET_NEI2005_MASK
      !--------------------------------------
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: NEI2005_SCALE_FUTURE
      PRIVATE :: INIT_NEI2005_ANTHRO
      PRIVATE :: TOTAL_ANTHRO_TG
      PRIVATE :: READ_NEI2005_MASK
      PRIVATE :: GET_NEI99_SEASON
      PRIVATE :: GET_NEI99_SEASON_05x0666
      PRIVATE :: GET_VISTAS_SEASON
      PRIVATE :: GET_VISTAS_SEASON_05x0666
      PRIVATE :: GET_NEI99_WKSCALE
      PRIVATE :: GET_NEI99_WKSCALE_05x0666
!
! !REMARKS:
!   (1) NIT is available in the data file but not read here (it is not
!        emitted in GEOS-Chem). 
!   (2) The algorithms in routines EMISS_NEI2005_ANTHRO and 
!        EMISS_NEI2005_ANTHRO_05x0666 may cause the code to die
!        when running offline simulations.  We will add a fix later.
!     
! !REVISION HISTORY:
!  07 Oct 2009 - A. van Donkelaar - initial version
!  20 Oct 2009 - P. Le Sager - added handling of VOC & masks
!  02 Nov 2009 - A. van Donkelaar - added seasonality, weekday factors
!  02 Dec 2009 - R. Yantosca - Added GET_NEI2005_MASK function
!  02 Dec 2009 - R. Yantosca - Updated comments etc.
!  10 Dec 2009 - D. Millet - Fix scaling, which is by ozone season
!  11 Dec 2009 - L. Zhang, A. Van Donkelaar - Add seasonality for NH3 
!  21 Dec 2009 - R. Yantosca - Added support for 0.5 x 0.666 nested grids
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  27 Jul 2011 - R. Yantosca - Removed typo in EMISS_NEI2005_ANTHRO_05x0666
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
      ! Array for surface area
      REAL*8,  ALLOCATABLE :: A_CM2(:)

      ! Arrays for emissions
      REAL*8,  ALLOCATABLE :: NOx(:,:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:,:)
      REAL*8,  ALLOCATABLE :: SO4(:,:,:)
      REAL*8,  ALLOCATABLE :: NH3(:,:,:)
      REAL*8,  ALLOCATABLE :: OC(:,:,:)
      REAL*8,  ALLOCATABLE :: BC(:,:,:)

      REAL*8,  ALLOCATABLE :: ALK4(:,:,:) ! 105
      REAL*8,  ALLOCATABLE :: ACET(:,:,:) ! 109
      REAL*8,  ALLOCATABLE :: MEK (:,:,:) ! 110
      REAL*8,  ALLOCATABLE :: ALD2(:,:,:) ! 111
      REAL*8,  ALLOCATABLE :: PRPE(:,:,:) ! 118
      REAL*8,  ALLOCATABLE :: C2H6(:,:,:) ! 121
      REAL*8,  ALLOCATABLE :: C3H8(:,:,:) ! 119
      REAL*8,  ALLOCATABLE :: CH2O(:,:,:) ! 120
      
      REAL*8,  ALLOCATABLE :: NOx_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: CO_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: SO2_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: SO4_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: NH3_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: OC_WKEND(:,:,:)
      REAL*8,  ALLOCATABLE :: BC_WKEND(:,:,:)

      REAL*8,  ALLOCATABLE :: ALK4_WKEND(:,:,:) ! 105
      REAL*8,  ALLOCATABLE :: ACET_WKEND(:,:,:) ! 109
      REAL*8,  ALLOCATABLE :: MEK_WKEND(:,:,:)  ! 110
      REAL*8,  ALLOCATABLE :: ALD2_WKEND(:,:,:) ! 111
      REAL*8,  ALLOCATABLE :: PRPE_WKEND(:,:,:) ! 118
      REAL*8,  ALLOCATABLE :: C2H6_WKEND(:,:,:) ! 121
      REAL*8,  ALLOCATABLE :: C3H8_WKEND(:,:,:) ! 119
      REAL*8,  ALLOCATABLE :: CH2O_WKEND(:,:,:) ! 120
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei2005_anthro
!
! !DESCRIPTION: Function GET\_NEI2005\_ANTHRO returns the NEI2005
!  emission for GEOS-Chem grid box (I,J,L) and tracer N.  Emissions can be 
!  returned in units of [kg/s] or [molec/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_NEI2005_ANTHRO( I,    J,     L, N, WEEKDAY,
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
! !USES:
!
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
      USE TRACERID_MOD, ONLY : IDTNOx,  IDTCO,   IDTSO2,  IDTNH3
      USE TRACERID_MOD, ONLY : IDTSO4
!
! !INPUT PARAMETERS: 
!
      ! Longitude, latitude, and tracer indices
      INTEGER, INTENT(IN)           :: I, J, L, N

      ! OPTIONAL -- return emissions in [molec/cm2/s]
      LOGICAL, INTENT(IN), OPTIONAL :: WEEKDAY, MOLEC_CM2_S  

      ! OPTIONAL -- return emissions in [kg/s] or [kg C/s]
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S
!
! !RETURN VALUE:
!     
      ! Emissions output
      REAL*8                        :: VALUE     
!
! !REVISION HISTORY: 
!    7 Oct 2009 - A. van Donkelaar - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_NEI2005_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      IF ( WEEKDAY ) THEN
         IF ( N == IDTNOx ) THEN

            ! NOx [kg/yr]
            VALUE = NOx(I,J,L)

         ELSE IF ( N == IDTCO ) THEN

            ! CO [kg/yr]
            VALUE = CO(I,J,L)

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2(I,J,L)

         ELSE IF ( N == IDTSO4 ) THEN

            ! SO4 [kg/yr]
            VALUE = SO4(I,J,L)

         ELSE IF ( N == IDTNH3 ) THEN

            ! NH3 [kg/yr]
            VALUE = NH3(I,J,L)

         ELSE IF ( N == IDTALK4 ) THEN

            ! [kg C/yr]
            VALUE = ALK4(I,J,L)

         ELSE IF ( N == IDTACET ) THEN

            ! [kg C/yr]
            VALUE = ACET(I,J,L)

         ELSE IF ( N == IDTMEK ) THEN

            ! [kg C/yr]
            VALUE = MEK(I,J,L)

         ELSE IF ( N == IDTPRPE ) THEN

            ! [kg C/yr]
            VALUE = PRPE(I,J,L)

         ELSE IF ( N == IDTC3H8 ) THEN

            ! [kg C/yr]
            VALUE = C3H8(I,J,L)

         ELSE IF ( N == IDTCH2O ) THEN

            ! [kg C/yr]
            VALUE = CH2O(I,J,L)

         ELSE IF ( N == IDTC2H6 ) THEN

            ! [kg C/yr]
            VALUE = C2H6(I,J,L)

         ELSE IF ( N == IDTALD2 ) THEN

            ! [kg C/yr]
            VALUE = ALD2(I,J,L)

         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no NEI2005 emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ELSE
         IF ( N == IDTNOx ) THEN

            ! NOx [kg/yr]
            VALUE = NOx_WKEND(I,J,L)

         ELSE IF ( N == IDTCO ) THEN
   
            ! CO [kg/yr]
            VALUE = CO_WKEND(I,J,L)

         ELSE IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2_WKEND(I,J,L)

         ELSE IF ( N == IDTSO4 ) THEN

            ! SO4 [kg/yr]
            VALUE = SO4_WKEND(I,J,L)

         ELSE IF ( N == IDTNH3 ) THEN
   
            ! NH3 [kg/yr]
            VALUE = NH3_WKEND(I,J,L)

         ELSE IF ( N == IDTALK4 ) THEN

            ! [kg C/yr]
            VALUE = ALK4_WKEND(I,J,L)

         ELSE IF ( N == IDTACET ) THEN

            ! [kg C/yr]
            VALUE = ACET_WKEND(I,J,L)

         ELSE IF ( N == IDTMEK ) THEN

            ! [kg C/yr]
            VALUE = MEK_WKEND(I,J,L)

         ELSE IF ( N == IDTPRPE ) THEN

            ! [kg C/yr]
            VALUE = PRPE_WKEND(I,J,L)

         ELSE IF ( N == IDTC3H8 ) THEN

            ! [kg C/yr]
            VALUE = C3H8_WKEND(I,J,L)

         ELSE IF ( N == IDTCH2O ) THEN

            ! [kg C/yr]
            VALUE = CH2O_WKEND(I,J,L)

         ELSE IF ( N == IDTC2H6 ) THEN

            ! [kg C/yr]
            VALUE = C2H6_WKEND(I,J,L)
 
         ELSE IF ( N == IDTALD2 ) THEN

            ! [kg C/yr]
            VALUE = ALD2_WKEND(I,J,L)

         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no NEI2005 emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
            
         ! Convert from [kg/yr] to [kg/s] or from [kgC/yr] to [kgC/s]
         VALUE = VALUE / SEC_IN_YEAR

      ELSE IF ( DO_MCS ) THEN

         ! Convert NOx from [kg/yr] to [molec/cm2/s] or from
         ! [kg C/yr] to [atom C/cm2/s] 
         VALUE = VALUE * XNUMOL(N) / ( A_CM2(J) * SEC_IN_YEAR )

      ENDIF

      ! Return to calling program
      END FUNCTION GET_NEI2005_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_nei2005_anthro
!
! !DESCRIPTION: Subroutine EMISS\_NEI2005\_ANTHRO reads the NEI2005
!  emission fields at 1x1 resolution and regrids them to the 
!  current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_NEI2005_ANTHRO
!
! !USES:
! 
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_1x1
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
      USE TRACERID_MOD, ONLY : IDTNOx,  IDTCO,   IDTSO2,  IDTNH3
      USE TRACERID_MOD, ONLY : IDTSO4,  IDTOCPI, IDTBCPI

      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
!
! !REVISION HISTORY: 
!  07 Oct 2009 - A. van Donkelaar - initial version
!  20 Oct 2009 - P. Le Sager - added VOC, account for mask to get better total
!  12 Jul 2010 - R. Yantosca - Now point to NEI2005_201007 directory, to read
!                              in updated files (by Aaron van Donkelaar) to
!                              fix a problem in the VOC emissions.
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J,   THISYEAR,       SNo, ScNo
      INTEGER                    :: L, KLM, SPECIES_ID(15), ID,  MN
      INTEGER                    :: OFFLINE_ID(15)
      REAL*4                     :: ARRAY(I1x1,J1x1,5)
      REAL*8                     :: GEOS_1x1(I1x1,J1x1,5)
      REAL*8                     :: SC_1x1(I1x1,J1x1)
      REAL*8                     :: TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR
      CHARACTER(LEN=5)           :: SNAME
      CHARACTER(LEN=1)           :: SSMN
      CHARACTER(LEN=2)           :: SMN

      !=================================================================
      ! EMISS_NEI2005_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_NEI2005_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

#if   defined( GEOS_5 ) || defined( MERRA )
      SNAME = 'GEOS5'
#elif defined( GEOS_4 )
      SNAME = 'GEOS4'
#elif defined( GEOS_3 )
      SNAME = 'GEOS3'
#endif

      ! list of ID of available species
      SPECIES_ID = (/ IDTNOX,  IDTCO,   IDTSO2,  IDTSO4, IDTNH3,
     &                IDTACET, IDTALK4, IDTC2H6, IDTC3H8,
     &                IDTOCPI, IDTBCPI,
     &                IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
     &     /)

      ! ID #'s for that are not tied to IDTxxxx flags
      OFFLINE_ID = (/ 1,       4,       26,      27,      30,
     &                9,       5,       21,      19,      35,           
     &                34,      11,      20,      18,      10       /)

      ! Loop over species
      DO KLM = 1, SIZE( SPECIES_ID )

         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            SNo = SPECIES_ID( KLM ) 
         ELSE
            SNo = OFFLINE_ID( KLM )
         ENDIF

         ! Skip undefined tracers
         IF ( SNo == 0 ) CYCLE

         ! corresponding annual scale factor # if any
         ScNo = 0
         IF ( SNo == 1                 ) ScNo = 71
         IF ( SNo == 4                 ) ScNo = 72
         IF ( SNo == 26 .or. SNo == 27 ) ScNo = 73
         
         ! TAU values for 2005
         TAU2005 = GET_TAU0( 1, 1, 2005 )

         ! File name
         FILENAME  = TRIM( DATA_DIR_1x1 ) // 'NEI2005_201007/' //
     &               'NEI2005.' // TRIM( SNAME ) // '.1x1.AVG.bpch'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - EMISS_NEI2005_ANTHRO: Reading ', a )

         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo, 
     &                    TAU2005,  I1x1,       J1x1,     
     &                    5,        ARRAY,      QUIET=.TRUE. ) 

         ! Cast to REAL*8 before regridding
         GEOS_1x1(:,:,:) = ARRAY(:,:,:)

         ! Apply annual scalar factor. Available for 1985-2005,
         ! and NOx, CO and SO2 only.
         IF ( ScNo .ne. 0 ) THEN

            CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2005, 
     &                                  THISYEAR, SC_1x1 )
            
            DO L = 1, 5
               GEOS_1x1(:,:,L) = GEOS_1x1(:,:,L) * SC_1x1(:,:)
            ENDDO

         ENDIF

         ! Apply Seasonality
         IF ( SNo .eq. 1 ) THEN
            CALL GET_VISTAS_SEASON( ARRAY )
         ELSE
            CALL GET_NEI99_SEASON( SNo, ARRAY )
         ENDIF
         GEOS_1x1(:,:,:) = GEOS_1x1(:,:,:) * ARRAY(:,:,:)

         
         ! Get Weekday/Weekend scaling
         CALL GET_NEI99_WKSCALE( SNo, ARRAY )

         
         ! Regrid from GEOS 1x1 --> current model resolution
         IF ( SNo .eq. 1 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, NOx )
            CALL DO_REGRID_1x1( 5, 'kg/yr', 
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), NOx_WKEND )
            DO L = 1, 5
               NOx(:,:,L) = NOx(:,:,L) * USA_MASK(:,:)
               NOx_WKEND(:,:,L) = 
     &                NOx_WKEND(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 4 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, CO )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), CO_WKEND )
            DO L = 1, 5
               CO(:,:,L) = CO(:,:,L) * USA_MASK(:,:)
               CO_WKEND(:,:,L) = 
     &                CO_WKEND(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 26 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, SO2 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), SO2_WKEND )
            DO L = 1, 5
               SO2_WKEND(:,:,L) =
     &                SO2_WKEND(:,:,L) * USA_MASK(:,:)
               SO2(:,:,L) = SO2(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 27 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, SO4 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), SO4_WKEND )
            DO L = 1, 5
               SO4_WKEND(:,:,L) =
     &                SO4_WKEND(:,:,L) * USA_MASK(:,:)
               SO4(:,:,L) = SO4(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 30 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, NH3 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), NH3_WKEND )
            DO L = 1, 5
               NH3_WKEND(:,:,L) =
     &                NH3_WKEND(:,:,L) * USA_MASK(:,:)
               NH3(:,:,L) = NH3(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 35 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, OC )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), OC_WKEND )
            DO L = 1, 5
               OC_WKEND(:,:,L) =
     &                OC_WKEND(:,:,L) * USA_MASK(:,:)
               OC(:,:,L) = OC(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 34 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, BC )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), BC_WKEND )
            DO L = 1, 5
               BC_WKEND(:,:,L) =
     &               BC_WKEND(:,:,L) * USA_MASK(:,:)
               BC(:,:,L) = BC(:,:,L) * USA_MASK(:,:)
            ENDDO

         !--VOC   
         ELSEIF ( SNo == 5 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, ALK4 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), ALK4_WKEND )
            DO L = 1, 5
               ALK4_WKEND(:,:,L) =
     &                ALK4_WKEND(:,:,L) * USA_MASK(:,:)
               ALK4(:,:,L) = ALK4(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 9 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, ACET )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), ACET_WKEND )
            DO L = 1, 5
               ACET_WKEND(:,:,L) =
     &                ACET_WKEND(:,:,L) * USA_MASK(:,:)
               ACET(:,:,L) = ACET(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 10 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, MEK )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), MEK_WKEND )
            DO L = 1, 5
               MEK_WKEND(:,:,L) =
     &                MEK_WKEND(:,:,L) * USA_MASK(:,:)
               MEK(:,:,L) = MEK(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 18 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, PRPE )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), PRPE_WKEND )
            DO L = 1, 5
               PRPE_WKEND(:,:,L) =
     &                PRPE_WKEND(:,:,L) * USA_MASK(:,:)
               PRPE(:,:,L) = PRPE(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 19 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, C3H8 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), C3H8_WKEND )
            DO L = 1, 5
               C3H8_WKEND(:,:,L) =
     &                C3H8_WKEND(:,:,L) * USA_MASK(:,:)
               C3H8(:,:,L) = C3H8(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 20 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, CH2O )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), CH2O_WKEND )
            DO L = 1, 5
               CH2O_WKEND(:,:,L) =
     &                CH2O_WKEND(:,:,L) * USA_MASK(:,:)
               CH2O(:,:,L) = CH2O(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 21 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, C2H6 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), C2H6_WKEND )
            DO L = 1, 5
               C2H6_WKEND(:,:,L) =
     &                C2H6_WKEND(:,:,L) * USA_MASK(:,:)
               C2H6(:,:,L) = C2H6(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 11 ) THEN

            CALL DO_REGRID_1x1( 5, 'kg/yr', GEOS_1x1, ALD2 )
            CALL DO_REGRID_1x1( 5, 'kg/yr',
     &           GEOS_1x1(:,:,:) * ARRAY(:,:,:), ALD2_WKEND )
            DO L = 1, 5
               ALD2_WKEND(:,:,L) =
     &                ALD2_WKEND(:,:,L) * USA_MASK(:,:)
               ALD2(:,:,L) = ALD2(:,:,L) * USA_MASK(:,:)
            ENDDO

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL NEI2005_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_NEI2005_ANTHRO
!EOC
!------------------------------------------------------------------------------
!       Dalhousie University Atmospheric Compositional Analysis Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_nei2005_anthro_05x0666
!
! !DESCRIPTION: Subroutine EMISS\_NEI2005\_ANTHRO reads the NEI2005
!  emission fields at 1/2 x 2.3 resolution
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_NEI2005_ANTHRO_05x0666
!
! !USES:
!
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR
      USE LOGICAL_MOD,       ONLY : LFUTURE
      USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH
      USE SCALE_ANTHRO_MOD,  ONLY : GET_ANNUAL_SCALAR_05x0666_NESTED
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
      USE TRACERID_MOD, ONLY : IDTNOx,  IDTCO,   IDTSO2,  IDTNH3
      USE TRACERID_MOD, ONLY : IDTSO4,  IDTOCPI, IDTBCPI

      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
!
! !REVISION HISTORY:
!  03 Nov 2009 - A. van Donkelaar - initial version
!  12 Jul 2010 - R. Yantosca - Now point to NEI2005_201007 directory, to read
!                              in updated files (by Aaron van Donkelaar) to
!                              fix a problem in the VOC emissions.
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5 (leave for future use)
!  27 Jul 2011 - R. Yantosca - Fixed typo: now *really* point to the NEI2005
!                              data directory NEI2005_101007/
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J,   THISYEAR,       SNo, ScNo
      INTEGER                    :: L, KLM, SPECIES_ID(15), ID,  MN
      INTEGER                    :: OFFLINE_ID(15)
      REAL*4                     :: ARRAY(IIPAR,JJPAR,5)
      REAL*8                     :: GEOS_05x0666(IIPAR,JJPAR,5)
      REAL*4                     :: SC_05x0666(IIPAR,JJPAR)
      REAL*8                     :: TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR
      CHARACTER(LEN=5)           :: SNAME
      CHARACTER(LEN=1)           :: SSMN
      CHARACTER(LEN=2)           :: SMN

      !=================================================================
      ! EMISS_NEI2005_ANTHRO begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_NEI2005_ANTHRO
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

#if   defined( GEOS_5 ) || defined( MERRA )
      SNAME = 'GEOS5'
#elif defined( GEOS_4 )
      SNAME = 'GEOS4'
#elif defined( GEOS_3 )
      SNAME = 'GEOS3'
#endif

      ! list of ID of available species
      SPECIES_ID = (/ IDTNOX,  IDTCO,   IDTSO2,  IDTSO4, IDTNH3,
     &                IDTACET, IDTALK4, IDTC2H6, IDTC3H8,
     &                IDTOCPI, IDTBCPI,
     &                IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
     &     /)

      ! ID #'s for that are not tied to IDTxxxx flags
      OFFLINE_ID = (/ 1,       4,       26,      27,      30,
     &                9,       5,       21,      19,      35,           
     &                34,      11,      20,      18,      10       /)

      ! Loop over species
      DO KLM = 1, SIZE( SPECIES_ID )

         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            SNo = SPECIES_ID( KLM ) 
         ELSE
            SNo = OFFLINE_ID( KLM )
         ENDIF

         ! Skip undefined tracers
         IF ( SNo == 0 ) CYCLE

         ! corresponding annual scale factor # if any
         ScNo = 0
         IF ( SNo == 1                 ) ScNo = 71
         IF ( SNo == 4                 ) ScNo = 72
         IF ( SNo == 26 .or. SNo == 27 ) ScNo = 73

         ! TAU values for 2005
         TAU2005 = GET_TAU0( 1, 1, 2005 )

         ! File name
         FILENAME  = TRIM( DATA_DIR ) // 'NEI2005_201007/' //
     &            'NEI2005.' // TRIM( SNAME ) 
     &             // '.1t2x2t3.AVG.na.bpch'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - EMISS_NEI2005_ANTHRO_05x0666: 
     &                   Reading ', a )

         CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                    TAU2005,  IIPAR,       JJPAR,
     &                    5,        ARRAY,      QUIET=.TRUE. )

         GEOS_05x0666(:,:,:) = ARRAY(:,:,:)

         ! Apply annual scalar factor. Available for 1985-2005,
         ! and NOx, CO and SO2 only.
         IF ( ScNo .ne. 0 ) THEN

            CALL GET_ANNUAL_SCALAR_05x0666_NESTED( ScNo,2005,
     &                             THISYEAR, SC_05x0666 )

            DO L = 1, 5
               GEOS_05x0666(:,:,L) = GEOS_05x0666(:,:,L) 
     &                               * SC_05x0666(:,:)
            ENDDO

         ENDIF

         ! Apply Seasonality
         IF ( SNo .eq. 1 ) THEN
            CALL GET_VISTAS_SEASON_05x0666( ARRAY )
         ELSE
            CALL GET_NEI99_SEASON_05x0666( SNo, ARRAY )
         ENDIF
         GEOS_05x0666(:,:,:) = GEOS_05x0666(:,:,:) 
     &                         * ARRAY(:,:,:)

         CALL GET_NEI99_WKSCALE_05x0666( SNo, ARRAY )

         IF ( SNo .eq. 1 ) THEN

            NOx(:,:,:) = GEOS_05x0666
            NOx_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               NOx(:,:,L) = NOx(:,:,L) * USA_MASK(:,:)
               NOx_WKEND(:,:,L) =
     &                NOx_WKEND(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 4 ) THEN

            CO(:,:,:) = GEOS_05x0666
            CO_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               CO(:,:,L) = CO(:,:,L) * USA_MASK(:,:)
               CO_WKEND(:,:,L) =
     &                CO_WKEND(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 26 ) THEN

            SO2(:,:,:) = GEOS_05x0666
            SO2_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               SO2_WKEND(:,:,L) =
     &                SO2_WKEND(:,:,L) * USA_MASK(:,:)
               SO2(:,:,L) = SO2(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 27 ) THEN

            SO4(:,:,:) = GEOS_05x0666
            SO4_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               SO4_WKEND(:,:,L) =
     &                SO4_WKEND(:,:,L) * USA_MASK(:,:)
               SO4(:,:,L) = SO4(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 30 ) THEN

            NH3(:,:,:) = GEOS_05x0666
            NH3_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               NH3_WKEND(:,:,L) =
     &                NH3_WKEND(:,:,L) * USA_MASK(:,:)
               NH3(:,:,L) = NH3(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 35 ) THEN

            OC(:,:,:) = GEOS_05x0666
            OC_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               OC_WKEND(:,:,L) =
     &                OC_WKEND(:,:,L) * USA_MASK(:,:)
               OC(:,:,L) = OC(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSEIF ( SNo .eq. 34 ) THEN

            BC(:,:,:) = GEOS_05x0666
            BC_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               BC_WKEND(:,:,L) =
     &               BC_WKEND(:,:,L) * USA_MASK(:,:)
               BC(:,:,L) = BC(:,:,L) * USA_MASK(:,:)
            ENDDO

         !--VOC
         ELSEIF ( SNo == 5 ) THEN

            ALK4(:,:,:) = GEOS_05x0666
            ALK4_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               ALK4_WKEND(:,:,L) =
     &                ALK4_WKEND(:,:,L) * USA_MASK(:,:)
               ALK4(:,:,L) = ALK4(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 9 ) THEN

            ACET(:,:,:) = GEOS_05x0666
            ACET_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               ACET_WKEND(:,:,L) =
     &                ACET_WKEND(:,:,L) * USA_MASK(:,:)
               ACET(:,:,L) = ACET(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 10 ) THEN

            MEK(:,:,:) = GEOS_05x0666
            MEK_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               MEK_WKEND(:,:,L) =
     &                MEK_WKEND(:,:,L) * USA_MASK(:,:)
               MEK(:,:,L) = MEK(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 18 ) THEN

            PRPE(:,:,:) = GEOS_05x0666
            PRPE_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               PRPE_WKEND(:,:,L) =
     &                PRPE_WKEND(:,:,L) * USA_MASK(:,:)
               PRPE(:,:,L) = PRPE(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 19 ) THEN

            C3H8(:,:,:) = GEOS_05x0666
            C3H8_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               C3H8_WKEND(:,:,L) =
     &                C3H8_WKEND(:,:,L) * USA_MASK(:,:)
               C3H8(:,:,L) = C3H8(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 20 ) THEN

            CH2O(:,:,:) = GEOS_05x0666
            CH2O_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               CH2O_WKEND(:,:,L) =
     &                CH2O_WKEND(:,:,L) * USA_MASK(:,:)
               CH2O(:,:,L) = CH2O(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 21 ) THEN

            C2H6(:,:,:) = GEOS_05x0666
            C2H6_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               C2H6_WKEND(:,:,L) =
     &                C2H6_WKEND(:,:,L) * USA_MASK(:,:)
               C2H6(:,:,L) = C2H6(:,:,L) * USA_MASK(:,:)
            ENDDO

         ELSE IF ( SNo == 11 ) THEN

            ALD2(:,:,:) = GEOS_05x0666
            ALD2_WKEND(:,:,:) = GEOS_05x0666 * ARRAY(:,:,:)
            DO L = 1, 5
               ALD2_WKEND(:,:,L) =
     &                ALD2_WKEND(:,:,L) * USA_MASK(:,:)
               ALD2(:,:,L) = ALD2(:,:,L) * USA_MASK(:,:)
            ENDDO

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN
         CALL NEI2005_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

      ! Return to calling program
      END SUBROUTINE EMISS_NEI2005_ANTHRO_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei99_season
!
! !DESCRIPTION: Subroutine GET\_NEI99\_SEASON returns monthly scale
!  factors from EPA 1999
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NEI99_SEASON( TRACER, AS )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_TAU0,    READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1
      USE TIME_MOD,      ONLY : GET_MONTH
      USE TRACERID_MOD,  ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD,  ONLY : IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
      USE TRACERID_MOD,  ONLY : IDTNOx,  IDTCO,   IDTSO2,  IDTNH3
      USE TRACERID_MOD,  ONLY : IDTSO4

      USE CMN_SIZE_MOD                         ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)    :: TRACER   ! Tracer number
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(OUT)   :: AS(I1x1,J1x1,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  30 Oct 2009 - A. van Donkelaar - Initial Version
!   3 Nov 2009 - P. Le Sager      - update handling of boxes w/ zero emissions
!  10 Dec 2009 - D. Millet        - Now scale to August, not an annual average
!  11 Dec 2009 - L. Zhang, A. van Donkelaar - Add seasonality for NH3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                 :: ARRAY(I1x1,J1x1,1)
!      REAL*4                 :: ANNUAL(I1x1,J1x1,1) dbm, 12/9/2009
      REAL*4                 :: AUGUST(I1x1,J1x1,1) ! dbm, 12/9/2009
      REAL*4                 :: MONTHLY(I1x1,J1x1,1)
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=6)       :: MYEAR
      REAL*8                 :: TAU
      INTEGER                :: MN, ThisMN, L

      ! seasonal scalar for NH3 emission (lzh, amv, 12/11/2009)
      REAL*8, PARAMETER :: NH3_SCALE(12) = (/ 
     &     0.426d0, 0.445d0, 0.526d0, 0.718d0, 1.179d0, 1.447d0, 
     &     1.897d0, 1.884d0, 1.577d0, 0.886d0, 0.571d0, 0.445d0 /)
 
      !=================================================================
      ! GET_NEI99_SEASON begins here!
      !=================================================================

      ARRAY(:,:,1) = 0.d0
!      ANNUAL(:,:,1) = 0.d0 dbm, 12/9/2009
      AUGUST(:,:,1) = 0.d0 ! dbm, 12/9/2009
      MONTHLY(:,:,1) = 0.d0

      ThisMN = GET_MONTH()

      ! lzh, amv, 12/11/2009 add NH3 emission seasonality
      IF ( TRACER == 11 .or. TRACER == 20 ) THEN
         AS = 1.d0
         RETURN
      ELSEIF ( TRACER == 30 ) THEN
         AS = NH3_SCALE(ThisMN) / NH3_SCALE(8)   ! Normalize to August
         RETURN
      ENDIF 

      ! Echo info
      WRITE( 6, 100 ) TRACER
 100  FORMAT( '     - GET_NEI99_SEASON: Reading TRACER: ', i )

      !---------------------------------
      ! Read in data for August
      !---------------------------------

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &     'EPA_NEI_200708/wkday_avg_an.199908.geos.1x1'

      ! TAU0 for 1999/08/01
      TAU = GET_TAU0( 8, 1, 1999 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', TRACER,
     &                 TAU,       I1x1,      J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      AUGUST(:,:,1) = ARRAY(:,:,1)

      !---------------------------------
      ! Read in data for current month
      !---------------------------------

      WRITE(MYEAR, '(i6)') 199900 + ThisMN
      
      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &     'EPA_NEI_200708/wkday_avg_an.' // MYEAR // '.geos.1x1'

      ! TAU for this month of 1999
      TAU = GET_TAU0( ThisMN, 1, 1999 )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', TRACER,
     &                 TAU,       I1x1,      J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )
      
      MONTHLY(:,:,1) = ARRAY(:,:,1)

      !---------------------------------
      ! Normalize
      !------------- -------------------

      WHERE ( AUGUST == 0d0 )
         ARRAY = 1d0
      ELSEWHERE
         ARRAY = MONTHLY / AUGUST
      ENDWHERE

      DO L = 1, SIZE(AS,3)
         AS(:,:,L) = ARRAY(:,:,1)
      ENDDO
         
      END SUBROUTINE GET_NEI99_SEASON
!EOC
!------------------------------------------------------------------------------
!      Dalhousie University Atmospheric Composition Analysis Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei99_season_05x0666
!
! !DESCRIPTION: Subroutine GET\_NEI\_SEASON returns monthly scale
!  factors from EPA 1999, for the 0.5 x 0.666 nested grids.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NEI99_SEASON_05x0666( TRACER, AS )
!
! !USES:
!
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1

      USE CMN_SIZE_MOD                         ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)       :: TRACER   ! Tracer number
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(INOUT) :: AS(IIPAR,JJPAR,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  30 Oct 2009 - A. van Donkelaar - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                     :: ARRAY(I1x1,J1x1,5)
      REAL*8                     :: ARRAY_R8(IIPAR,JJPAR,5)

      !=================================================================
      ! GET_NEI99_SEASON_05x0666 begins here!
      !=================================================================

      ARRAY(:,:,:) = 0.d0

      CALL GET_NEI99_SEASON( TRACER, ARRAY )

      CALL DO_REGRID_1x1( 5, 'unitless', ARRAY, ARRAY_R8 )
      AS(:,:,:) = ARRAY_R8(:,:,:)

      END SUBROUTINE GET_NEI99_SEASON_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_vistas_season
!
! !DESCRIPTION: Subroutine GET\_VISTAS\_SEASON returns monthly scale
!  factors to account for monthly variations in NOx emissions
!  on 1x1 resolution grid (amv, 11/02/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_VISTAS_SEASON( AS )
!
! !USES:
!
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE TIME_MOD,          ONLY : GET_MONTH,     GET_YEAR

      USE CMN_SIZE_MOD                         ! Size parameters
      USE CMN_O3_MOD            ! FSCALYR
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(INOUT) :: AS(I1x1,J1x1,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  30 Oct 2009 - A. van Donkelaar - Initial Version
!   3 Nov 2009 - P. Le Sager      - update handling of boxes w/ zero emissions
!  10 Dec 2009 - D. Millet        - Now scale to August, not an annual average
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                     :: ARRAY(I1x1,J1x1,1)
      REAL*4                     :: AUGUST(I1x1,J1x1,1) ! dbm, 12/9/2009
      REAL*4                     :: MONTHLY(I1x1,J1x1,1)
      REAL*4                     :: O3SEASON(I1x1,J1x1,1)
      REAL*4                     :: O3SEASON_AUGUST(I1x1,J1x1,1) ! dbm, 12/9/2009
      CHARACTER(LEN=255)         :: FILENAME, VISTAS_DIR
      CHARACTER(LEN=4)           :: SYEAR
      CHARACTER(LEN=1)           :: SSMN
      CHARACTER(LEN=2)           :: SMN
      REAL*8                     :: TAU2002
      INTEGER                    :: MN, THISMONTH, LEV
      INTEGER                    :: THISYEAR

      !=================================================================
      ! GET_NEI99_SEASON begins here!
      !=================================================================

      ARRAY(:,:,1) = 0.d0
      AUGUST(:,:,1) = 0.d0 ! dbm, 12/9/2009
      MONTHLY(:,:,1) = 0.d0
      O3SEASON(:,:,1) = 0.d0
      O3SEASON_AUGUST(:,:,1) = 0.d0 ! dbm, 12/9/2009

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF

      ! cap maximum scaling year
      IF ( THISYEAR .gt. 2007 ) THEN
         THISYEAR = 2007
      ENDIF

      VISTAS_DIR = TRIM( DATA_DIR_1x1 ) // 'VISTAS_200811/'

      TAU2002 = GET_TAU0( 1, 1, 2002)
      THISMONTH = GET_MONTH()

      ! -------------------
      ! Read in data for August
      ! -------------------

      FILENAME  = TRIM( VISTAS_DIR )
     &     // 'Vistas-NOx-8.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - GET_VISTAS_SEASON: Reading ', a )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,
     &                 TAU2002,   I1x1,      J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )

      AUGUST(:,:,1) = ARRAY(:,:,1)

      ! -------------------
      ! Read in data for current month
      ! -------------------

      IF (THISMONTH .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) THISMONTH
         FILENAME  = TRIM( VISTAS_DIR )
     &        // 'Vistas-NOx-' // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) THISMONTH
         FILENAME  = TRIM( VISTAS_DIR )
     &        // 'Vistas-NOx-' // SMN // '.1x1'
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1,
     &                 TAU2002,   I1x1,      J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )
  
      MONTHLY(:,:,1) = ARRAY(:,:,1)
      
      WRITE( SYEAR, '(i4)') THISYEAR

      ! Load ozone season regulation factors
      IF (THISMONTH .lt. 10) THEN
         WRITE( SSMN, '(i1)' ) THISMONTH
         FILENAME  = TRIM( VISTAS_DIR )
     &         // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &         // SSMN // '.1x1'
      ELSE
         WRITE( SMN, '(i2)' ) THISMONTH
         FILENAME  = TRIM( VISTAS_DIR )
     &         // 'ARP-SeasonalVariation-' // SYEAR // '-'
     &         // SMN // '.1x1'
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'RATIO-2D', 71,
     &                 GET_TAU0(1,1,2002),  I1x1,    J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )
      O3SEASON(:,:,1) = ARRAY(:,:,1)

      ! August ozone season regulation factors
      FILENAME  = TRIM( VISTAS_DIR )
     &     // 'ARP-SeasonalVariation-' // SYEAR // '-8.1x1'
     
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data
      CALL READ_BPCH2( FILENAME, 'RATIO-2D', 71,
     &                 GET_TAU0(1,1,2002),  I1x1,    J1x1,
     &                 1,         ARRAY,     QUIET=.TRUE. )
      O3SEASON_AUGUST(:,:,1) = ARRAY(:,:,1)

      ! First do seasonal scaling according to VISTAS
      WHERE ( AUGUST == 0d0 )
         ARRAY = 1d0
      ELSEWHERE
         ARRAY = MONTHLY / AUGUST
      ENDWHERE
      
      ! Now scale for summertime NOx reductions
      ARRAY = ARRAY * O3SEASON / O3SEASON_AUGUST

      DO LEV = 1, SIZE(AS,3)
         AS(:,:,LEV) = ARRAY(:,:,1)
      ENDDO
      
      
      END SUBROUTINE GET_VISTAS_SEASON
!EOC
!------------------------------------------------------------------------------
!      Dalhousie University Atmospheric Composition Analysis Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_vistas_season_05x0666
!
! !DESCRIPTION: Subroutine GET\_VISTAS\_SEASON\_05x0666 returns monthly scale
!  factors to account for monthly variations in NOx emissions
!  for the 0.5 x 0.666 nested grids. (amv, 11/02/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_VISTAS_SEASON_05x0666( AS )
!
! !USES:
!
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1

      USE CMN_SIZE_MOD                         ! Size parameters
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(INOUT) :: AS(IIPAR,JJPAR,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  03 Nov 2009 - A. van Donkelaar - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                     :: ARRAY(I1x1,J1x1,5)
      REAL*8                     :: ARRAY_R8(IIPAR,JJPAR,5)

      !=================================================================
      ! GET_VISTAS_SEASON_05x0666 begins here!
      !=================================================================

      ARRAY(:,:,:) = 0.d0

      CALL GET_VISTAS_SEASON( ARRAY )

      CALL DO_REGRID_1x1( 5, 'unitless', ARRAY, ARRAY_R8 )
      AS(:,:,:) = ARRAY_R8(:,:,:)

      END SUBROUTINE GET_VISTAS_SEASON_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei99_wkscale
!
! !DESCRIPTION: Subroutine GET\_NEI99\_WKSCALE returns the scale
!  factors to convert weekday to weekend emissions based 
!  on the NEI99.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NEI99_WKSCALE( TRACER, AS )
!
! !USES:
!
      USE BPCH2_MOD,         ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
      USE TIME_MOD,          ONLY : GET_MONTH
      USE TRACERID_MOD, ONLY : IDTACET, IDTALK4, IDTC2H6, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTALD2, IDTCH2O, IDTPRPE, IDTMEK
      USE TRACERID_MOD, ONLY : IDTNOx,  IDTCO,   IDTSO2,  IDTNH3
      USE TRACERID_MOD, ONLY : IDTSO4

      USE CMN_SIZE_MOD                         ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)       :: TRACER   ! Tracer number
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(INOUT) :: AS(I1x1,J1x1,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  30 Oct 2009 - A. van Donkelaar - Initial Version
!   3 Nov 2009 - P. Le Sager - update handling of boxes w/ zero emissions
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                     :: WEEKDAY(I1x1,J1x1,1)
      REAL*4                     :: WEEKEND(I1x1,J1x1,1)
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=6)           :: MYEAR
      REAL*8                     :: TAU
      INTEGER                    :: MN, L

      !=================================================================
      ! GET_NEI99_WKSCALE begins here!
      !=================================================================

      WEEKDAY(:,:,1) = 0.d0
      WEEKEND(:,:,1) = 0.d0

      MN = GET_MONTH()

      ! NH3/ALD2/ISOP not available
      IF (( TRACER .eq. 30 ) .or. (TRACER .eq. 11 ) .or.
     &    ( TRACER .eq. 20 )) THEN
         AS(:,:,:) = 1.d0
         RETURN
      ENDIF

      ! Echo info
      WRITE( 6, 100 ) TRACER
 100  FORMAT( '     - GET_NEI99_WKSCALE: Reading TRACER: ', i )

      WRITE(MYEAR, '(i6)') 199900 + MN

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) //
     &   'EPA_NEI_200708/wkday_avg_an.' // MYEAR // '.geos.1x1'

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', TRACER,
     &                 GET_TAU0(MN,1,1999), I1x1,   J1x1,
     &                 1,  WEEKDAY,   QUIET=.TRUE. )

      WRITE(MYEAR, '(i6)') 199900 + MN

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) //
     &   'EPA_NEI_200708/wkend_avg_an.' // MYEAR // '.geos.1x1'

      ! Read data
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', TRACER,
     &                 GET_TAU0(MN,1,1999), I1x1,   J1x1,
     &                 1,  WEEKEND,   QUIET=.TRUE. )

!---see below
!     ! avoid 0 / 0
!      WEEKDAY(:,:,1) = WEEKDAY(:,:,1) + 1.d0
!      WEEKEND(:,:,1) = WEEKEND(:,:,1) + 1.d0
!
!      DO L = 1,5
!         AS(:,:,L) = WEEKEND(:,:,1) / WEEKDAY(:,:,1)
!      ENDDO


      ! --Get scalings
      WHERE ( WEEKDAY == 0d0 )
         WEEKEND = 1d0
      ELSEWHERE
         WEEKEND = WEEKEND / WEEKDAY
      ENDWHERE

      DO L = 1, SIZE(AS,3)
         AS(:,:,L) = WEEKEND(:,:,1)
      ENDDO

      
      END SUBROUTINE GET_NEI99_WKSCALE
!EOC
!------------------------------------------------------------------------------
!      Dalhousie University Atmospheric Composition Analysis Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_nei99_wkscale_05x0666
!
! !DESCRIPTION: Subroutine GET\_NEI99\_WKSCALE\_05x0666 returns the scale
!  factors (for 0.5 x 0.666 nested grids) to convert weekday to weekend 
!  emissions based on the NEI99.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_NEI99_WKSCALE_05x0666( TRACER, AS )
!
! !USES:
!
      USE REGRID_1x1_MOD,    ONLY : DO_REGRID_1x1

      USE CMN_SIZE_MOD                         ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)       :: TRACER   ! Tracer number
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4,  INTENT(INOUT) :: AS(IIPAR,JJPAR,5)  ! Scale factor array
!
! !REVISION HISTORY:
!  30 Oct 2009 - A. van Donkelaar - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                     :: ARRAY(I1x1,J1x1,5)
      REAL*8                     :: ARRAY_R8(IIPAR,JJPAR,5)

      !=================================================================
      ! GET_NEI99_SEASON_05x0666 begins here!
      !=================================================================

      ARRAY(:,:,:) = 0.d0

      CALL GET_NEI99_WKSCALE( TRACER, ARRAY )

      CALL DO_REGRID_1x1( 5, 'unitless', ARRAY, ARRAY_R8 )
      AS(:,:,:) = ARRAY_R8(:,:,:)

      END SUBROUTINE GET_NEI99_WKSCALE_05x0666
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nei2005_mask
!
! !DESCRIPTION: Subroutine READ\_NEI2005\_MASK reads the mask for NEI data  
!\\
!\\
! !INTERFACE:
      
      SUBROUTINE READ_NEI2005_MASK
!
! !USES:
!     
      ! Reference to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE LOGICAL_MOD,    ONLY : LCAC,            LBRAVO
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD        ! Size parameters
!
! !REMARKS:
!     temporary mask: same as EPA 99
!     
! !REVISION HISTORY: 
!   20 Oct 2009 - P. Le Sager - init
!   26 Oct 2009 - P. Le Sager - new masks
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4             :: ARRAY2(I1x1,J1x1,1)
      REAL*8             :: XTAU, GEOS_1x1(I1x1,J1x1,1)
      CHARACTER(LEN=255) :: FILENAME, SNAME

      !=================================================================
      ! Mask specific to NEI2005 data
      !=================================================================
      
      SNAME = 'usa.'

      ! NEI2005 covers CANADA if we do not use CAC     
      IF ( .NOT. LCAC ) SNAME = TRIM( SNAME ) // 'can.'

      ! NEI2005 covers Mexico if we do not use BRAVO      
      IF ( .NOT. LBRAVO ) SNAME = TRIM( SNAME ) // 'mex.'

      
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 'NEI2005_200910/' //      
     &     TRIM( SNAME ) // 'mask.nei2005.geos.1x1'

      ! Echo info
      WRITE( 6, 200 ) TRIM( FILENAME )
200   FORMAT( '     - READ_NEI2005_MASK: Reading ', a )
     

      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,      I1x1,      J1x1,     
     &                 1,        ARRAY2,    QUIET=.TRUE. ) 

      
      ! Cast to REAL*8 before regridding
      GEOS_1x1(:,:,:) = ARRAY2(:,:,:)

      ! Regrid from GEOS 1x1 --> current model resolution
      CALL DO_REGRID_1x1( 'unitless', GEOS_1x1, USA_MASK )

      WHERE ( USA_MASK /= 0D0 ) USA_MASK = 1D0

      
      ! Return to calling program
      END SUBROUTINE READ_NEI2005_MASK
!------------------------------------------------------------------------------
! Prior to 12/3/09:
! Leave for future use (bmy, 12/3/09)
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: get_nei2005_mask
!!
!! !DESCRIPTION: Subroutine GET\_NEI2005\_MASK returns the value of the 
!!  NEI 2005 mask to the calling program.  Values of 1 denote grid boxes 
!!  within the EPA/NEI2005 emission region.!  
!!\\
!!\\
!! !INTERFACE:
!      
!      FUNCTION GET_NEI2005_MASK( I, J ) RESULT ( USA )
!!
!! !INPUT PARAMETERS:
!!     
!      INTEGER, INTENT(IN) :: I, J   ! GEOS-Chem lon & lat indices
!!
!! !RETURN VALUE:
!!
!      REAL*8              :: USA    ! Value of the mask  
!!
!! !REMARKS:
!!  This is entended to encapsulate the USA_MASK variable.
!!     
!! !REVISION HISTORY: 
!!  02 Dec 2009 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      USA = USA_MASK(I,J)
!
!      END FUNCTION GET_NEI2005_MASK
!------------------------------------------------------------------------------
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nei2005_scale_future
!
! !DESCRIPTION: Subroutine NEI2005\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the NEI2005 anthropogenic emissions.
!\\
!\\
! !INTERFACE:

      SUBROUTINE NEI2005_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_OCff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_BCff

      USE CMN_SIZE_MOD             ! Size parameters
!
! !REMARKS:
!    VOC are not scaled, however scale factors are available (see
!     epa_nei_mod.f for procedure)
!     
! !REVISION HISTORY: 
!    7 Oct 2009 - A. van Donkelaar - initial version
!   20 Oct 2009 - P. Le Sager - set L OpenMP private, put L loop first  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J, L

      !=================================================================
      ! NEI2005_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, 5
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         NOx(I,J,L)  = NOx(I,J,L) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /yr]
         CO(I,J,L)   = CO(I,J,L)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/yr] 
         SO2(I,J,L)  = SO2(I,J,L) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future SO4 [kg SO4/yr]
         SO4(I,J,L)  = SO4(I,J,L) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future NH3 [kg NH3/yr] 
         NH3(I,J,L)  = NH3(I,J,L) * GET_FUTURE_SCALE_NH3an( I, J )

         ! Future OC [kg NH3/yr]
         OC(I,J,L)  = OC(I,J,L) * GET_FUTURE_SCALE_OCff( I, J )

         ! Future BC [kg NH3/yr]
         BC(I,J,L)  = BC(I,J,L) * GET_FUTURE_SCALE_BCff( I, J )

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE NEI2005_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx, CO, SO2 and NH3.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR )
!
! !USES:
! 
      USE CMN_SIZE_MOD            ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: YEAR   ! Year of data to compute totals
!
! !REVISION HISTORY: 
!    7 Oct 2009 - A. van Donkelaar - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L
      REAL*8              :: T_NOX,  T_CO,  T_SO2,  T_NH3
      REAL*8              :: T_SO4,  T_OC,  T_BC,   T_ALK4
      REAL*8              :: T_ACET, T_MEK, T_PRPE, T_C3H8
      REAL*8              :: T_CH2O, T_C2H6,T_ALD2
      CHARACTER(LEN=3)    :: UNIT
                                                   
      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'N. E. I. 2005 U. S. A.   E M I S S I O N S', / )


      ! Total NOx [Tg N]
      T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )

      ! Total CO  [Tg CO]
      T_CO  = SUM( CO  ) * 1d-9

      ! Total SO2 [Tg S]
      T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

      ! Total SO4 [Tg S]
      T_SO4 = SUM( SO4 ) * 1d-9 * ( 32d0 / 96d0 )

      ! Total NH3 [Tg NH3]
      T_NH3 = SUM( NH3 ) * 1d-9

      ! Total OC [Tg]
      T_OC = SUM( OC ) * 1d-9

      ! Total OC [Tg]
      T_BC = SUM( BC ) * 1d-9 

      ! Total ALK4 [Tg C]
      T_ALK4 = SUM( ALK4 ) * 1d-9

      ! Total ACET [Tg C]
      T_ACET = SUM( ACET ) * 1d-9

      ! Total MEK [Tg C]
      T_MEK = SUM( MEK ) * 1d-9

      ! Total PRPE [Tg C]
      T_PRPE = SUM( PRPE ) * 1d-9

      ! Total C3H8 [Tg C]
      T_C3H8 = SUM( C3H8 ) * 1d-9

      ! Total CH2O [Tg C]
      T_CH2O = SUM( CH2O ) * 1d-9

      ! Total C2H6 [Tg C]
      T_C2H6 = SUM( C2H6 ) * 1d-9

      ! Total ALD2 [Tg C]
      T_ALD2 = SUM( ALD2 ) * 1d-9


      
      ! Print totals in [Tg]
      WRITE( 6, 110 ) 'NOx ',  YEAR, T_NOx,  '[Tg N  ]'
      WRITE( 6, 110 ) 'CO  ',  YEAR, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'SO2 ',  YEAR, T_SO2,  '[Tg S  ]'
      WRITE( 6, 110 ) 'SO4 ',  YEAR, T_SO4,  '[Tg S  ]'
      WRITE( 6, 110 ) 'NH3 ',  YEAR, T_NH3,  '[Tg NH3]'
      WRITE( 6, 110 ) 'OC ' ,  YEAR, T_OC,   '[Tg C]'
      WRITE( 6, 110 ) 'BC ' ,  YEAR, T_BC,   '[Tg C]'
      WRITE( 6, 110 ) 'ALK4 ', YEAR, T_ALK4, '[Tg C]'
      WRITE( 6, 110 ) 'ACET ', YEAR, T_ACET, '[Tg C]'
      WRITE( 6, 110 ) 'MEK ' , YEAR, T_MEK,  '[Tg C]'
      WRITE( 6, 110 ) 'PRPE ', YEAR, T_PRPE, '[Tg C]'
      WRITE( 6, 110 ) 'C3H8 ', YEAR, T_C3H8, '[Tg C]'
      WRITE( 6, 110 ) 'CH2O ', YEAR, T_CH2O, '[Tg C]'
      WRITE( 6, 110 ) 'C2H6 ', YEAR, T_C2H6, '[Tg C]'
      WRITE( 6, 110 ) 'ALD2 ', YEAR, T_ALD2, '[Tg C]'

      ! Format statement
 110  FORMAT( 'NEI2005 anthro ', a5, 
     &        'for year ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_Tg
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_nei2005_anthro
!
! !DESCRIPTION: Subroutine INIT\_NEI2005\_ANTHRO allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_NEI2005_ANTHRO
!
! !USES:
! 
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE LOGICAL_MOD, ONLY : LNEI05

      USE CMN_SIZE_MOD    ! Size parameters
!
! !REVISION HISTORY: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_NEI2005_ANTHRO begins here!
      !=================================================================

      ! Return if LNEI05 is false
      IF ( .not. LNEI05 ) RETURN
      
      !--------------------------------------------------
      ! Allocate and zero arrays for emissions
      !--------------------------------------------------

      ! allocate and read USA Mask
      ALLOCATE( USA_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'USA_MASK' )
      USA_MASK = 0d0

      CALL READ_NEI2005_MASK

      ALLOCATE( NOx( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx' )
      NOx = 0d0

      ALLOCATE( CO( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
      CO = 0d0

      ALLOCATE( SO2( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
      SO2 = 0d0

      ALLOCATE( SO4( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4' )
      SO4 = 0d0

      ALLOCATE( NH3( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0

      ALLOCATE( OC( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OC' )
      OC = 0d0 

      ALLOCATE( BC( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
      BC = 0d0 

      ALLOCATE( ALK4( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK4' )
      ALK4 = 0d0

      ALLOCATE( ACET( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACET' )
      ACET = 0d0

      ALLOCATE( MEK( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MEK' )
      MEK = 0d0

      ALLOCATE( ALD2( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALD2' )
      ALD2 = 0d0

      ALLOCATE( PRPE( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRPE' )
      PRPE = 0d0

      ALLOCATE( C2H6( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C2H6' )
      C2H6 = 0d0 

      ALLOCATE( C3H8( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C3H8' )
      C3H8 = 0d0 

      ALLOCATE( CH2O( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH2O' )
      CH2O = 0d0 

      ALLOCATE( NOx_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOx_WKEND' )
      NOx_WKEND = 0d0

      ALLOCATE( CO_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO_WKEND' )
      CO_WKEND = 0d0

      ALLOCATE( SO2_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2_WKEND' )
      SO2_WKEND = 0d0

      ALLOCATE( SO4_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4_WKEND' )
      SO4_WKEND = 0d0

      ALLOCATE( NH3_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3_WKEND' )
      NH3_WKEND = 0d0

      ALLOCATE( OC_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OC_WKEND' )
      OC_WKEND = 0d0 

      ALLOCATE( BC_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC_WKEND' )
      BC_WKEND = 0d0 

      ALLOCATE( ALK4_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK4_WKEND' )
      ALK4_WKEND = 0d0

      ALLOCATE( ACET_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACET_WKEND' )
      ACET_WKEND = 0d0

      ALLOCATE( MEK_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MEK_WKEND' )
      MEK_WKEND = 0d0

      ALLOCATE( ALD2_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALD2_WKEND' )
      ALD2_WKEND = 0d0

      ALLOCATE( PRPE_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRPE_WKEND' )
      PRPE_WKEND = 0d0

      ALLOCATE( C2H6_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C2H6_WKEND' )
      C2H6_WKEND = 0d0 

      ALLOCATE( C3H8_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'C3H8_WKEND' )
      C3H8_WKEND = 0d0 

      ALLOCATE( CH2O_WKEND( IIPAR, JJPAR, 5 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH2O_WKEND' )
      CH2O_WKEND = 0d0 

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
      END SUBROUTINE INIT_NEI2005_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_nei2005_anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_NEI2005\_ANTHRO deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_NEI2005_ANTHRO
!
! !REVISION HISTORY: 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_NEIO2005_ANTHRO begins here!
      !=================================================================
      ! USA mask
      IF ( ALLOCATED( USA_MASK) ) DEALLOCATE( USA_MASK )
      IF ( ALLOCATED( A_CM2   ) ) DEALLOCATE( A_CM2 )
      IF ( ALLOCATED( NOx     ) ) DEALLOCATE( NOx   )
      IF ( ALLOCATED( CO      ) ) DEALLOCATE( CO    )
      IF ( ALLOCATED( SO2     ) ) DEALLOCATE( SO2   )
      IF ( ALLOCATED( SO4     ) ) DEALLOCATE( SO4   )
      IF ( ALLOCATED( NH3     ) ) DEALLOCATE( NH3   )
      IF ( ALLOCATED( OC      ) ) DEALLOCATE( OC    )
      IF ( ALLOCATED( BC      ) ) DEALLOCATE( BC    )
      IF ( ALLOCATED( ALK4    ) ) DEALLOCATE( ALK4  )
      IF ( ALLOCATED( ACET    ) ) DEALLOCATE( ACET  )
      IF ( ALLOCATED( MEK     ) ) DEALLOCATE( MEK   )
      IF ( ALLOCATED( ALD2    ) ) DEALLOCATE( ALD2  )
      IF ( ALLOCATED( PRPE    ) ) DEALLOCATE( PRPE  )
      IF ( ALLOCATED( C2H6    ) ) DEALLOCATE( C2H6  )
      IF ( ALLOCATED( C3H8    ) ) DEALLOCATE( C3H8  )
      IF ( ALLOCATED( CH2O    ) ) DEALLOCATE( CH2O  )
      IF (ALLOCATED(NOx_WKEND) ) DEALLOCATE(NOx_WKEND )
      IF (ALLOCATED(CO_WKEND  )) DEALLOCATE(CO_WKEND  )
      IF (ALLOCATED(SO2_WKEND )) DEALLOCATE(SO2_WKEND )
      IF (ALLOCATED(SO4_WKEND )) DEALLOCATE(SO4_WKEND )
      IF (ALLOCATED(NH3_WKEND )) DEALLOCATE(NH3_WKEND )
      IF (ALLOCATED(OC_WKEND  )) DEALLOCATE(OC_WKEND  )
      IF (ALLOCATED(BC_WKEND  )) DEALLOCATE(BC_WKEND  )
      IF (ALLOCATED(ALK4_WKEND)) DEALLOCATE(ALK4_WKEND)
      IF (ALLOCATED(ACET_WKEND)) DEALLOCATE(ACET_WKEND)
      IF (ALLOCATED(MEK_WKEND )) DEALLOCATE(MEK_WKEND )
      IF (ALLOCATED(ALD2_WKEND)) DEALLOCATE(ALD2_WKEND)
      IF (ALLOCATED(PRPE_WKEND)) DEALLOCATE(PRPE_WKEND)
      IF (ALLOCATED(C2H6_WKEND)) DEALLOCATE(C2H6_WKEND)
      IF (ALLOCATED(C3H8_WKEND)) DEALLOCATE(C3H8_WKEND)
      IF (ALLOCATED(CH2O_WKEND)) DEALLOCATE(CH2O_WKEND)
      
      END SUBROUTINE CLEANUP_NEI2005_ANTHRO
!EOC
      END MODULE NEI2005_ANTHRO_MOD
