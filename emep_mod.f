! $Id: emep_mod.f,v 1.9 2009/01/29 15:37:23 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: EMEP_MOD
!
! !DESCRIPTION: \subsection*{Overview}
!  Module EMEP\_MOD contains variables and routines to read the 
!  EMEP European anthropogenic emission inventory for CO, NOz, and some 
!  NMVOCs.  The EMEP files come from Marion Auvray and Isabelle Bey at EPFL. 
!  (bdf, bmy, amv, phs, 11/1/05, 1/28/09)
!
!\subsection*{References}
! \begin{enumerate}
!     \item Vestreng, V., and H. Klein (2002), \emph{Emission data reported 
!           to UNECE/EMEP: Quality insurance and trend analysis and 
!           presentation of Web-Dab}, \underline{MSC-W Status Rep}. 2002:, 
!           101 pp., Norw. Meteorol. Inst., Oslo, Norway.  This paper is 
!           on the EMEP web site:
!           http://www.emep.int/mscw/mscw\_publications.html
!           http://www.emep.int/publ/reports/2002/mscw\_note\_1\_2002.pdf
!     \item Auvray, M., and I. Bey, \emph{Long-Range Transport to Europe: 
!           Seasonal Variations and Implications for the European Ozone 
!           Budget}, \underline{J. Geophys. Res.}, \textbf{110}, D11303, 
!           doi: 10.1029/2004JD005503, 2005.
! \end{enumerate}
!
! !INTERFACE: 
!
      MODULE EMEP_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMISS_EMEP
      PUBLIC  :: CLEANUP_EMEP
      PUBLIC  :: GET_EUROPE_MASK
      PUBLIC  :: GET_EMEP_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: EMEP_SCALE_FUTURE
      PRIVATE :: READ_EUROPE_MASK 
      PRIVATE :: INIT_EMEP        
!
! !REVISION HISTORY:
!  (1 ) Now only print totals for defined tracers (bmy, 2/6/06)
!  (2 ) Now modified for IPCC future emissions (swu, bmy, 5/30/06)
!  (3 ) Now yearly scale factors can be applied (phs, amv, 3/17/08)
!  (4 ) Now include emep SOx and emep emissions to 2005 (amv, 06/08)
!  (5 ) Modify to access SHIP emissions from outside (phs, 06/08)
!  (6 ) Account for monthly variations (amv, 12/9/08)
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!
      ! Array for geographic mask
      REAL*8,  ALLOCATABLE :: EUROPE_MASK(:,:)

      ! Arrays for ground-based emissions
      REAL*8,  ALLOCATABLE :: EMEP_NOx(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_CO(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_SO2(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_NH3(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_ALK4(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_MEK(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_ALD2(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_PRPE(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_C2H6(:,:)

      ! Arrays for ship emissions
      REAL*8,  ALLOCATABLE :: EMEP_CO_SHIP(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_SO2_SHIP(:,:)
      REAL*8,  ALLOCATABLE :: EMEP_NOx_SHIP(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_EUROPE_MASK
!
! !DESCRIPTION: Function GET\_EUROPE\_MASK returns the value of the EUROPE 
!  mask for EMEP emissions at grid box (I,J).  MASK=1 if (I,J) is in the 
!  European region, or MASK=0 otherwise. (bdf, bmy, 11/1/05)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EUROPE_MASK( I, J ) RESULT( EUROPE )
!
! !INPUT ARGUMENTS:
!
      INTEGER, INTENT(IN) :: I        ! Longitude index
      INTEGER, INTENT(IN) :: J        ! Latitude  index
!
! RETURN VALUE:
! 
      REAL*8              :: EUROPE   ! Returns the mask value @ (I,J)
!
! !REVISION HISTORY: 
!  01 Nov 2005 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! GET_EUROPE_MASK begins here!
      !=================================================================
      EUROPE = EUROPE_MASK(I,J)

      ! Return to calling program
      END FUNCTION GET_EUROPE_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_EMEP_ANTHRO
!
! !DESCRIPTION: Function GET\_EMEP\_ANTHRO returns the EMEP emission for 
!  GEOS-CHEM grid box (I,J) and tracer N. (bdf, bmy, 11/1/05)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EMEP_ANTHRO( I, J, N, KG_S, SHIP ) RESULT( EMEP )
!
! !USES:
! 
      USE TRACERID_MOD, ONLY : IDTNOX,  IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD, ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD, ONLY : IDTNH3
      USE TRACER_MOD,   ONLY : XNUMOL
      USE GRID_MOD,     ONLY : GET_AREA_CM2
!
! !INPUT ARGUMENTS:
!
      INTEGER, INTENT(IN)           :: I       ! Longitude index
      INTEGER, INTENT(IN)           :: J       ! Latitude index
      INTEGER, INTENT(IN)           :: N       ! Tracer number
      LOGICAL, INTENT(IN), OPTIONAL :: KG_S    ! Return emissions in [kg/s]
      LOGICAL, INTENT(IN), OPTIONAL :: SHIP    ! Return ship emissions
!
! RETURN VALUE:
! 
      REAL*8                        :: EMEP    ! Returns emissions at (I,J)
!
! !REVISION HISTORY: 
!  (1 ) added SOx, SOx ship and NH3 emissions, plus optional kg/s output
!       (amv, 06/2008)
!  (2 ) Now returns ship emissions if requested (phs, 6/08)
!  (3 ) Added checks to avoid calling unavailable ship emissions (phs, 6/08)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS, IS_SHIP
      INTEGER                       :: NN,     HAS_SHIP(3)
      
      !=================================================================
      ! GET_EMEP_ANTHRO begins here!
      !=================================================================

      ! Initialize
      NN      = N
      IS_SHIP = .FALSE.
      DO_KGS  = .FALSE.
         
      IF ( PRESENT( KG_S ) ) DO_KGS = KG_S
      IF ( PRESENT( SHIP ) ) IS_SHIP = SHIP

      ! check SHIP availability
      HAS_SHIP = (/ IDTNOX, IDTCO, IDTSO2 /)

      IF ( IS_SHIP .AND. .NOT. ANY( HAS_SHIP == N) ) THEN
         WRITE(6,*)'WARNING: EMEP SHIP emissions not available for'//
     $             'tracer #',N
         EMEP = 0D0
         RETURN
      ENDIF

      ! NOx
      IF ( N  == IDTNOX ) THEN
         IF ( IS_SHIP ) THEN
            EMEP = EMEP_NOx_SHIP(I,J)
         ELSE 
            EMEP = EMEP_NOx(I,J)
         ENDIF

      ! CO
      ELSE IF ( N == IDTCO ) THEN
         IF ( IS_SHIP ) THEN
            EMEP = EMEP_CO_SHIP(I,J)
         ELSE 
            EMEP = EMEP_CO(I,J)
         ENDIF

      ! ALK4 (>= C4 alkanes)
      ELSE IF ( N == IDTALK4 ) THEN
         EMEP = EMEP_ALK4(I,J)

      ! MEK
      ELSE IF ( N == IDTMEK ) THEN
         EMEP = EMEP_MEK(I,J)

      ! ALD2 (acetaldehyde)
      ELSE IF ( N == IDTALD2 ) THEN
         EMEP = EMEP_ALD2(I,J)

      ! PRPE (>= C3 alkenes)
      ELSE IF ( N == IDTPRPE ) THEN
         EMEP = EMEP_PRPE(I,J)

      ! C2H6 
      ELSE IF ( N == IDTC2H6 ) THEN
         EMEP = EMEP_C2H6(I,J)

      ! SO2
      ELSE IF ( N == IDTSO2 ) THEN
         IF ( IS_SHIP ) THEN
            EMEP = EMEP_SO2_SHIP(I,J)
         ELSE 
            EMEP = EMEP_SO2(I,J)
         ENDIF

      ! NH3
      ELSE IF ( N == IDTNH3 ) THEN
         EMEP = EMEP_NH3(I,J)

      ! Otherwise return a negative value to indicate
      ! that there are no EMEP emissions for tracer N
      ELSE
         EMEP = -1d0

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN

         EMEP = EMEP * GET_AREA_CM2(J) / XNUMOL(NN)

      ENDIF

      ! Return to calling program
      END FUNCTION GET_EMEP_ANTHRO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SUBROUTINE EMISS_EMEP
!
! !DESCRIPTION: Subroutine EMISS\_EMEP reads the EMEP emission fields at 
!  1x1 resolution and regrids them to the current model resolution.
!  (bdf, bmy, 11/1/05, 5/30/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EMEP
!
! !USES:
! 
      USE BPCH2_MOD,        ONLY : GET_TAU0,     OPEN_BPCH2_FOR_READ
      USE FILE_MOD,         ONLY : IU_FILE,      IOERROR
      USE DIRECTORY_MOD,    ONLY : DATA_DIR_1x1 
      USE LOGICAL_MOD,      ONLY : LFUTURE
      USE REGRID_1x1_MOD,   ONLY : DO_REGRID_1x1
      USE TIME_MOD,         ONLY : EXPAND_DATE,  GET_YEAR
      USE TIME_MOD,         ONLY : GET_MONTH
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_O3"         ! SCALEYEAR
!
! !REVISION HISTORY: 
!  (1 ) Modified for IPCC future emissions.  Now references LFUTURE from
!        "logical_mod.f". (bmy, 5/30/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: EMEP_NYMD, EMEP_YEAR
      REAL*8                  :: EMEP_TAU,  TAU0
      CHARACTER(LEN=255)      :: FILENAME

      ! For bpch file format
      INTEGER                 :: I,  J,  L,  N,  IOS
      INTEGER                 :: NTRACER,   NSKIP
      INTEGER                 :: HALFPOLAR, CENTER180
      INTEGER                 :: NI,        NJ,        NL
      INTEGER                 :: IFIRST,    JFIRST,    LFIRST
      INTEGER                 :: SCALEYEAR
      REAL*4                  :: ARRAY(I1x1,J1x1,1)
      REAL*4                  :: LONRES,    LATRES
      REAL*4                  :: Sc(IIPAR,JJPAR)
      REAL*8                  :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)       :: MODELNAME
      CHARACTER(LEN=40)       :: CATEGORY
      CHARACTER(LEN=40)       :: UNIT     
      CHARACTER(LEN=40)       :: RESERVED

      !=================================================================
      ! EMISS_EMEP begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_EMEP
         FIRST = .FALSE.
      ENDIF

      ! 1x1 file name for EMEP 2000
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'EMEP_200510/EMEP.geos.1x1.YYYY'

      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      ! EMEP 2000 data is only defined from 1985-2000
      EMEP_YEAR = MAX( MIN( SCALEYEAR, 2000 ), 1985 )

      ! YYYYMMDD value for 1st day of EMEP_YEAR
      EMEP_NYMD = ( EMEP_YEAR * 10000 ) + 0101 

      ! TAU0 value corresponding to EMEP_NYMD
      EMEP_TAU  = GET_TAU0( 1, 1, EMEP_YEAR )
         
      ! Expand filename
      CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )
         
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_EMEP: Reading ', a )

      !=================================================================
      ! Read data at 1x1 resolution and regrid to current grid size
      !=================================================================

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      ! Read the entire file in one pass (for I/O optimization)
      DO 

         ! Read 1st data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! Check for EOF or errors
         IF ( IOS < 0 ) EXIT
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'emiss_emep:2' )

         ! Read 2nd data block header line
         READ( IU_FILE, IOSTAT=IOS ) 
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST, NSKIP

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emiss_emep:3' )

         ! Read data [molec/cm2/s] or [atoms C/cm2/s]
         READ( IU_FILE, IOSTAT=IOS ) 
     &        ( ( ( ARRAY(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'emiss_emep:4' )

         ! Regrid data from 1x1
         SELECT CASE ( NTRACER )

            ! NOx [molec/cm2/s]
            CASE( 1  )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_NOx  )

            ! CO [molec/cm2/s]
            CASE( 4  )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_CO   )

            ! ALK4 [atoms C/cm2/s]
            CASE( 5  )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_ALK4 )

            ! MEK [atoms C/cm2/s]
            CASE( 10 )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_MEK  )

            ! ALD2 [atoms C/cm2/s]
            CASE( 11 )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_ALD2 )

            ! PRPE [atoms C/cm2/s]
            CASE( 18 )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_PRPE )

            ! C2H6 [atoms C/cm2/s]
            CASE( 21 )
               CALL DO_REGRID_1x1( UNIT, ARRAY, EMEP_C2H6 )

            CASE DEFAULT
               ! Nothing

         END SELECT

      ENDDO

      ! Close file
      CLOSE( IU_FILE )

      !=================================================================
      ! Get and apply annual emissions factors (amv, phs, 3/17/08)
      !=================================================================

      !=================================================================
      ! If we are at or above 1990, can apply updated EMEP emissions for
      ! NOx, CO, NH3 and include SOx (amv, 06/04/08)
      !=================================================================

      IF ( SCALEYEAR > 1989 ) THEN

         ! new EMEP data is only defined from 1990-2005
         EMEP_YEAR = MIN( SCALEYEAR, 2005 )

         CALL READ_EMEP_UPDATED(  1, EMEP_YEAR, EMEP_NOx, 0 )
         CALL READ_EMEP_UPDATED(  4, EMEP_YEAR, EMEP_CO, 0 )
         CALL READ_EMEP_UPDATED( 26, EMEP_YEAR, EMEP_SO2, 0 )
         CALL READ_EMEP_UPDATED( 30, EMEP_YEAR, EMEP_NH3, 1 )

         CALL READ_EMEP_UPDATED(  1, EMEP_YEAR, EMEP_NOx_SHIP, 2 )
         CALL READ_EMEP_UPDATED(  4, EMEP_YEAR, EMEP_CO_SHIP, 2 )
         CALL READ_EMEP_UPDATED( 26, EMEP_YEAR, EMEP_SO2_SHIP, 2 )
 
      ! Need to use for SOx/NH3 anyways, but SOx scale back further
      ELSE

         CALL READ_EMEP_UPDATED( 26, 1990, EMEP_SO2, 0 )
         CALL READ_EMEP_UPDATED( 26, 1990, EMEP_SO2_SHIP, 2 )
         CALL READ_EMEP_UPDATED( 30, 1990, EMEP_NH3, 1 )

         CALL GET_ANNUAL_SCALAR( 73, 1990, SCALEYEAR, Sc )
         EMEP_SO2(:,:) = EMEP_SO2(:,:) * Sc(:,:)
!         EMEP_SO2_SHIP = EMEP_SO2_SHIP * Sc  ! do not scale SHIP

      ENDIF

      
      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================
      IF ( LFUTURE ) THEN 
         CALL EMEP_SCALE_FUTURE
      ENDIF

      !=================================================================
      ! Print emission totals
      !=================================================================

      ! Print totals for EMEP_YEAR
      CALL TOTAL_ANTHRO_TG( EMEP_YEAR, SCALEYEAR, GET_MONTH() )

      ! Return to calling program
      END SUBROUTINE EMISS_EMEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMEP_SCALE_FUTURE
!
! !DESCRIPTION: Subroutine EMEP\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the EMEP anthropogenic emissions. (swu, bmy, 5/30/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMEP_SCALE_FUTURE
!
! !USES:
! 
      ! References to F90 modules
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_TONEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff

#     include "CMN_SIZE"             ! Size parameters
!
! !REVISION HISTORY:
!  30 May 2006 - S. Wu & R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! EMEP_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [molec/cm2/s]
         EMEP_NOx(I,J)  = EMEP_NOx(I,J)                   * 
     &                    GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO [molec/cm2/s]
         EMEP_CO(I,J)   = EMEP_CO(I,J)                    *
     &                    GET_FUTURE_SCALE_COff( I, J )

         ! Future ALK4 [atoms C/cm2/s]
         EMEP_ALK4(I,J) = EMEP_ALK4(I,J)                  *
     &                    GET_FUTURE_SCALE_ALK4ff( I, J )
         
         ! Future MEK [atoms C/cm2/s]
         EMEP_MEK(I,J)  = EMEP_MEK(I,J)                   *
     &                    GET_FUTURE_SCALE_TONEff( I, J )     

         ! Future ALD2 [atoms C/cm2/s]
         EMEP_ALD2(I,J) = EMEP_ALD2(I,J)                  *
     &                    GET_FUTURE_SCALE_VOCff( I, J )
     
         ! Future PRPE [atoms C/cm2/s]
         EMEP_PRPE(I,J) = EMEP_PRPE(I,J)                  *
     &                    GET_FUTURE_SCALE_PRPEff( I, J )

         ! Future C2H6 [atoms C/cm2/s]
         EMEP_C2H6(I,J) = EMEP_C2H6(I,J)                  *
     &                    GET_FUTURE_SCALE_C2H6ff( I, J )
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMEP_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_ANTHRO_TG
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the amount of EMEP 
!  anthropogenic emissions that are emitted each month in Tg or Tg C. 
!  (rch, bmy, 11/10/04, 2/6/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( EMEP_YEAR, EMISS_YEAR, EMEP_MONTH )
!
! !USES:
! 
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LEMEPSHIP
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNOX,  IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD, ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD, ONLY : IDTNH3

#     include "CMN_SIZE"     ! Size parameters
!
! !INPUT ARGUMENTS:
!
      INTEGER, INTENT(IN)   :: EMEP_YEAR    ! EMEP base year
      INTEGER, INTENT(IN)   :: EMISS_YEAR   ! Current simulated year
      INTEGER, INTENT(IN)   :: EMEP_MONTH   ! Current simulated month
!
! !REVISION HISTORY: 
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now replace FMOL with TRACER_MW_KG (bmy, 10/25/05) 
!  (3 ) Now only print totals of defined tracers; other totals will be
!        printed as zeroes. (bmy, 2/6/06)
!  (4 ) Now emissions and base year are arguments. Output in Tg/month
!        since this is called monthly (phs, 12/9/08)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J
      REAL*8                :: A,   B(9), NOX,  CO,  ALK4
      REAL*8                :: MEK, ALD2, PRPE, C2H6, SO2
      REAL*8                :: NH3
      CHARACTER(LEN=3)      :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100 )
 100  FORMAT( 'M O N T H L Y   E M E P   E U R O P E A N
     $     E M I S S I O N S', / )
      
      ! indicate if we include ship emissions (automatic before 1990)
      IF ( LEMEPSHIP .OR. ( EMISS_YEAR < 1990 )) WRITE( 6, 101 )
 101  FORMAT( '( INCL. SHIP )', / )
      
      WRITE( 6, 102 ) EMEP_YEAR
 102  FORMAT( 'Base Year :', i4 )

      !----------------
      ! Sum emissions
      !----------------
      
      ! Define conversion factors for kg/molec
      ! (Undefined tracers will be zero)
      B(:) = 0d0
      IF ( IDTNOx  > 0 ) B(1) = 1d0 / ( 6.0225d23 / 14d-3 )  ! Tg N
      IF ( IDTCO   > 0 ) B(2) = 1d0 / XNUMOL(IDTCO  )
      IF ( IDTALK4 > 0 ) B(3) = 1d0 / XNUMOL(IDTALK4)
      IF ( IDTMEK  > 0 ) B(4) = 1d0 / XNUMOL(IDTMEK )
      IF ( IDTALD2 > 0 ) B(5) = 1d0 / XNUMOL(IDTALD2)
      IF ( IDTPRPE > 0 ) B(6) = 1d0 / XNUMOL(IDTPRPE)
      IF ( IDTC2H6 > 0 ) B(7) = 1d0 / XNUMOL(IDTC2H6)
      IF ( IDTSO2  > 0 ) B(8) = 1d0 / ( 6.0225d23 / 32d-3 )  ! Tg S
      IF ( IDTNH3  > 0 ) B(9) = 1d0 / XNUMOL(IDTNH3)


      ! Summing variables
      NOX      = 0d0   
      CO       = 0d0 
      ALK4     = 0d0 
      MEK      = 0d0 
      ALD2     = 0d0 
      PRPE     = 0d0 
      C2H6     = 0d0 
      SO2      = 0d0 
      NH3      = 0d0 

      ! Loop over latitudes
      DO J = 1, JJPAR
            
         ! Surface area [cm2] * seconds in this year 
         ! Multiply by 1d-9 to convert from [kg] to [Tg]
         A = GET_AREA_CM2( J ) * 365.25d0 * 86400d0 * 1d-9 
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Sum emissions (list NOx as Tg N)
            NOX  = NOX  + ( EMEP_NOX (I,J) + EMEP_NOX_SHIP(I,J) )
     $           * A * B(1)
            CO   = CO   + ( EMEP_CO  (I,J) + EMEP_CO_SHIP(I,J) )
     $           * A * B(2) 
            SO2  = SO2  + ( EMEP_SO2 (I,J) + EMEP_SO2_SHIP(I,J) )
     $           * A * B(8) 

            ALK4 = ALK4 + EMEP_ALK4(I,J) * A * B(3) 
            MEK  = MEK  + EMEP_MEK (I,J) * A * B(4) 
            ALD2 = ALD2 + EMEP_ALD2(I,J) * A * B(5) 
            PRPE = PRPE + EMEP_PRPE(I,J) * A * B(6) 
            C2H6 = C2H6 + EMEP_C2H6(I,J) * A * B(7) 
            NH3  = NH3  + EMEP_NH3 (I,J) * A * B(9) 
         ENDDO
      ENDDO
 
      !----------------
      ! Print sums
      !----------------

      ! Print totals in [kg/month]
      WRITE( 6, 110   ) 'NOx ', EMISS_YEAR, EMEP_MONTH, NOx/12d0,  ' N'
      WRITE( 6, 110   ) 'CO  ', EMISS_YEAR, EMEP_MONTH, CO/12d0,   '  '
      WRITE( 6, 110   ) 'SO2 ', EMISS_YEAR, EMEP_MONTH, SO2/12d0,  ' S'
      WRITE( 6, 110   ) 'NH3 ', EMISS_YEAR, EMEP_MONTH, NH3/12d0,  '  '
      WRITE( 6, 110   ) 'ALK4', EMISS_YEAR, EMEP_MONTH, ALK4/12d0, ' C'
      WRITE( 6, 110   ) 'MEK ', EMISS_YEAR, EMEP_MONTH, MEK/12d0,  ' C'
      WRITE( 6, 110   ) 'ALD2', EMISS_YEAR, EMEP_MONTH, ALD2/12d0, ' C'
      WRITE( 6, 110   ) 'PRPE', EMISS_YEAR, EMEP_MONTH, PRPE/12d0, ' C'
      WRITE( 6, 110   ) 'C2H6', EMISS_YEAR, EMEP_MONTH, C2H6/12d0, ' C'
 110  FORMAT( 'EMEP anthropogenic ', a4, ' for ', i4, '/', i2,
     &        ': ', f13.6, ' Tg', a2 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READ_EUROPE_MASK
!
! !DESCRIPTION: Subroutine READ\_EUROPE\_MASK reads and regrids the 
!  Europe mask for the EMEP anthropogenic emissions.  (bmy, 10/18/06)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EUROPE_MASK
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1

#     include "CMN_SIZE"       ! Size parameters
!
! !REVISION HISTORY: 
!  (1 ) Now read the Europe mask from a disk file instead of defining it as 
!        a rectangular box (bmy, 10/18/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_EUROPE_MASK begins here!
      !=================================================================

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'EMEP_200510/EMEP_mask.geos.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_EUROPE_MASK: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,       I1x1,     J1x1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Regrid from GEOS 1x1 GRID to current model resolution
      CALL DO_REGRID_1x1( 'unitless', ARRAY, EUROPE_MASK )

      ! Return to calling program
      END SUBROUTINE READ_EUROPE_MASK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READ_EMEP_UPDATED
!
! !DESCRIPTION: Subroutine READ\_EMEP\_UPDATED reads updated EMEP emissions 
!  from the year 1990 including SOx emissions.  These are regridded to the 
!  simulation resolution. Ship emissions can also be included. 
!  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EMEP_UPDATED( TRACER, EMEP_YEAR, ARRAY, wSHIP )
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : READ_BPCH2, GET_TAU0
      USE TIME_MOD,       ONLY : EXPAND_DATE, GET_MONTH
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1 
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE LOGICAL_MOD,    ONLY : LEMEPSHIP

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN_O3"         ! SCALEYEAR
!
! !INPUT ARGUMENTS:
!
      INTEGER, INTENT(IN)  :: TRACER              ! Tracer number
      INTEGER, INTENT(IN)  :: EMEP_YEAR           ! Year of emissions to read
      INTEGER, INTENT(IN)  :: wSHIP               ! Use ground, ship, or both?
!
! !OUTPUT ARGUMENTS:
!
      REAL*8,  INTENT(OUT) :: ARRAY(IIPAR,JJPAR)  ! Output array
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Now account for LEMEPSHIP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                        :: ARRAY_1x1(I1x1,J1x1,1)
      REAL*4                        :: ARRAY_1x1_SHIP(I1x1,J1x1,1)
      REAL*4                        :: ARRAY_1x1_LAND(I1x1,J1x1,1)
      CHARACTER(LEN=255)            :: FILENAME, DIR
      REAL*8                        :: EMEP_TAU, TAU2005
      INTEGER                       :: EMEP_NYMD, MN
      CHARACTER(LEN=2)              :: SMN
      CHARACTER(LEN=1)              :: SSMN

      ARRAY_1x1_SHIP(:,:,:) = 0.d0
      ARRAY_1x1_LAND(:,:,:) = 0.d0

      ! YYYYMMDD value for 1st day of EMEP_YEAR
      EMEP_NYMD = ( EMEP_YEAR * 10000 ) + 0101

      ! TAU0 value corresponding to EMEP_NYMD
      EMEP_TAU  = GET_TAU0( 1, 1, EMEP_YEAR )

      ! Expand filename
      DIR = TRIM( DATA_DIR_1x1 ) // 'EMEP_200806/'

      ! wSHIP = 0 means no ship emissions included
      ! wSHIP = 1 means include ships emissions
      ! wSHIP = 2 means only ship emissions

      IF ( wSHIP .lt. 2 ) THEN

         IF ( TRACER .eq. 1 ) THEN
            ! NOx
            FILENAME = TRIM( DIR ) // 'NOx/'
     &                // 'EMEP-NOx-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 4 ) THEN
            ! CO
            FILENAME = TRIM( DIR ) // 'CO/'
     &                // 'EMEP-CO-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 26 ) THEN
            ! SOx
            FILENAME = TRIM( DIR ) // 'SOx/'
     &                // 'EMEP-SOx-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 30 ) THEN
            ! NH3
            FILENAME = TRIM( DIR ) // 'NH3/'
     &                // 'EMEP-NH3-YYYY.geos.1x1'
         ENDIF

         CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )

         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_EMEP_UPDATED: Reading ', a )
 
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',      TRACER,
     &                    EMEP_TAU,  I1x1,           J1x1,
     &                    1,         ARRAY_1x1_LAND, QUIET=.TRUE. )

      ENDIF

      IF ( ( wSHIP .gt. 0 ) .AND. LEMEPSHIP ) THEN

         IF ( TRACER .eq. 1 ) THEN
            ! NOx
            FILENAME = TRIM( DIR ) // 'NOx/'
     &                // 'EMEP-NOx-SHIP-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 4 ) THEN
            ! CO
            FILENAME = TRIM( DIR ) // 'CO/'
     &                // 'EMEP-CO-SHIP-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 26 ) THEN
            ! SOx
            FILENAME = TRIM( DIR ) // 'SOx/'
     &                // 'EMEP-SOx-SHIP-YYYY.geos.1x1'
         ELSE IF ( TRACER .eq. 30 ) THEN
            ! NH3
            FILENAME = TRIM( DIR ) // 'NH3/'
     &                // 'EMEP-NH3-SHIP-YYYY.geos.1x1'

         ENDIF

         CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )

         WRITE( 6, 101 ) TRIM( FILENAME ) 
 101     FORMAT( '     - READ_EMEP_UPDATED: Reading ', a )

         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',      TRACER,
     &                    EMEP_TAU,  I1x1,           J1x1,
     &                    1,         ARRAY_1x1_SHIP, QUIET=.TRUE. )

      ENDIF

      ! Apply monthly variation (courtesy of the GENEMIS project 
      ! coordinated by the Institute of Energy Economics and the 
      ! Rational Use of Energy (IER) at the University of 
      ! Stuttgart) (amv, 11/24/2008)
      IF (( TRACER .eq. 1 ) .and. ( wSHIP .lt. 2 )) THEN

         ! Apply Monthly Factors over land
         TAU2005 = GET_TAU0( 1, 1, 2005)
         MN = GET_MONTH()

         IF (MN .lt. 10) THEN
            WRITE( SSMN, '(i1)' ) MN
            FILENAME = TRIM( DIR ) // 'SeasonalVariation/'
     &               // 'EMEP-SeasonalVariation-'
     &               // SSMN // '.1x1'
         ELSE
            WRITE( SMN, '(i2)' ) MN
            FILENAME = TRIM( DIR ) // 'SeasonalVariation/'
     &               // 'EMEP-SeasonalVariation-'
     &               // SMN // '.1x1'
         ENDIF

         ! Echo info
         WRITE( 6, 101 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', 71,
     &                    TAU2005,   I1x1,      J1x1,
     &                    1,         ARRAY_1x1,     QUIET=.TRUE. )

         ARRAY_1x1_LAND(:,:,1) = ARRAY_1x1_LAND(:,:,1) 
     &                           * ARRAY_1x1(:,:,1)

      ENDIF

      IF ( wSHIP .eq. 0 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_LAND(:,:,1)
      IF ( wSHIP .eq. 1 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_LAND(:,:,1) + 
     &                                       ARRAY_1x1_SHIP(:,:,1)
      IF ( wSHIP .eq. 2 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_SHIP(:,:,1)

      CALL DO_REGRID_1x1('molec/cm2/s', ARRAY_1x1, ARRAY)

      ! Convert SOx to SO2 assuming a SOx is 95% SO2 over Europe, as used
      ! throughout GEOS-Chem, and as per Chin et al, 2000
      IF ( TRACER .eq. 26 ) ARRAY(:,:) = ARRAY(:,:) * 0.95d0

      END SUBROUTINE READ_EMEP_UPDATED
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_EMEP
!
! !DESCRIPTION: Subroutine INIT\_EMEP allocates and zeroes EMEP module 
!  arrays, and also creates the mask which defines the European region.
!  (bdf, bmy, 11/1/05, 10/18/06) 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_EMEP
!
! !USES:
! 
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE GRID_MOD,    ONLY : GET_XMID, GET_YMID
      USE LOGICAL_MOD, ONLY : LEMEP

#     include "CMN_SIZE"    ! Size parameters
!
! !REVISION HISTORY: 
!  (1 ) Now call READ_EUROPE_MASK to read & regrid EUROPE_MASK from disk 
!        instead of just defining it as a rectangular box. (bmy, 10/18/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, I, J, X, Y

      !=================================================================
      ! INIT_EMEP begins here!
      !=================================================================

      ! Return if LEMEP is false
      IF ( .not. LEMEP ) RETURN
      
      !--------------------------------
      ! Allocate and zero arrays
      !--------------------------------
      ALLOCATE( EMEP_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_NOx' )
      EMEP_NOx = 0d0

      ALLOCATE( EMEP_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_CO' )
      EMEP_CO = 0d0

      ALLOCATE( EMEP_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2' )
      EMEP_SO2 = 0d0

      ALLOCATE( EMEP_SO2_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2_SHIP' )
      EMEP_SO2_SHIP = 0d0

      ALLOCATE( EMEP_CO_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2_SHIP' )
      EMEP_SO2_SHIP = 0d0

      ALLOCATE( EMEP_NOx_SHIP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2_SHIP' )
      EMEP_SO2_SHIP = 0d0

      ALLOCATE( EMEP_NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_NH3' )
      EMEP_NH3 = 0d0

      ALLOCATE( EMEP_ALK4( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_ALK4' )
      EMEP_ALK4 = 0d0

      ALLOCATE( EMEP_MEK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_MEK' )
      EMEP_MEK = 0d0

      ALLOCATE( EMEP_ALD2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_ALD2' )
      EMEP_ALD2 = 0d0

      ALLOCATE( EMEP_PRPE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_PRPE' )
      EMEP_PRPE = 0d0

      ALLOCATE( EMEP_C2H6( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP_C2H6' )
      EMEP_C2H6 = 0d0

      ALLOCATE( EUROPE_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EUROPE_MASK' )
      EUROPE_MASK = 0d0

      ! Read and regrid the European mask
      CALL READ_EUROPE_MASK

      ! Return to calling program
      END SUBROUTINE INIT_EMEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_EMEP
!
! !DESCRIPTION: Subroutine CLEANUP\_EMEP deallocates all module arrays 
!   (bmy, 11/1/05)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_EMEP
!
! !REVISION HISTORY: 
!  1 Nov 2005 - R. Yantosca - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_EMEP begins here!
      !=================================================================
      IF ( ALLOCATED( EMEP_NOx      ) ) DEALLOCATE( EMEP_NOx      )
      IF ( ALLOCATED( EMEP_CO       ) ) DEALLOCATE( EMEP_CO       )
      IF ( ALLOCATED( EMEP_SO2      ) ) DEALLOCATE( EMEP_SO2      )
      IF ( ALLOCATED( EMEP_SO2_SHIP ) ) DEALLOCATE( EMEP_SO2_SHIP )
      IF ( ALLOCATED( EMEP_CO_SHIP  ) ) DEALLOCATE( EMEP_CO_SHIP  )
      IF ( ALLOCATED( EMEP_NOx_SHIP ) ) DEALLOCATE( EMEP_NOx_SHIP )
      IF ( ALLOCATED( EMEP_NH3      ) ) DEALLOCATE( EMEP_NH3      )
      IF ( ALLOCATED( EMEP_ALK4     ) ) DEALLOCATE( EMEP_ALK4     )
      IF ( ALLOCATED( EMEP_MEK      ) ) DEALLOCATE( EMEP_MEK      )
      IF ( ALLOCATED( EMEP_ALD2     ) ) DEALLOCATE( EMEP_ALD2     )
      IF ( ALLOCATED( EMEP_PRPE     ) ) DEALLOCATE( EMEP_PRPE     )
      IF ( ALLOCATED( EMEP_C2H6     ) ) DEALLOCATE( EMEP_C2H6     )
      IF ( ALLOCATED( EUROPE_MASK   ) ) DEALLOCATE( EUROPE_MASK   )  

      ! Return to calling program
      END SUBROUTINE CLEANUP_EMEP

!------------------------------------------------------------------------------

      ! End of module
      END MODULE EMEP_MOD
!EOC
