!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emep_mod
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
!\begin{verbatim}
!  http://www.emep.int/mscw/mscw\_publications.html
!  http://www.emep.int/publ/reports/2002/mscw\_note\_1\_2002.pdf
!\end{verbatim}
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
      USE ERROR_MOD, ONLY: ALLOC_ERR

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EMISS_EMEP
      PUBLIC  :: EMISS_EMEP_05x0666
      PUBLIC  :: CLEANUP_EMEP
      PUBLIC  :: GET_EUROPE_MASK
      PUBLIC  :: GET_EMEP_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: EMEP_SCALE_FUTURE
      PRIVATE :: READ_EMEP_UPDATED
      PRIVATE :: READ_EMEP_UPDATED_05x0666
      PRIVATE :: READ_EUROPE_MASK 
      PRIVATE :: READ_EUROPE_MASK_05x0666
      PRIVATE :: INIT_EMEP        
!
! !REVISION HISTORY:
!  01 Nov 2005 - B. Field, R. Yantosca - Initial version
!  (1 ) Now only print totals for defined tracers (bmy, 2/6/06)
!  (2 ) Now modified for IPCC future emissions (swu, bmy, 5/30/06)
!  (3 ) Now yearly scale factors can be applied (phs, amv, 3/17/08)
!  (4 ) Now include emep SOx and emep emissions to 2005 (amv, 06/08)
!  (5 ) Modify to access SHIP emissions from outside (phs, 06/08)
!  (6 ) Account for monthly variations (amv, 12/9/08)
!  18 Dec 2009 - Aaron van D - Created routine EMISS_EMEP_05x0666
!  18 Dec 2009 - Aaron van D - Created routine READ_EMEP_UPDATED_05x0666
!  18 Dec 2009 - Aaron van D - Created routine READ_EUROPE_MASK_05x0666
!  11 Jan 2010 - Aaron van D - Max scale year is now 2007, for consistency
!  11 Jan 2010 - Aaron van D - Extend 1x1 emission files to 2007.  Routine
!                              READ_EMEP_UPDATED now mimics routine
!                              READ_EMEP_UPDATED_05x0666.
!  26 Jan 2010 - R. Yantosca - Minor bug fix in INIT_EMEP
!  31 Aug 2010 - R. Yantosca - Updated comments
!  24 Nov 2010 - G. Vinken   - Updated EMEP mask file
!  13 Mar 2012 - M. Cooper   - Changed gridding algorithm to map_a2a
!  22 Mar 2012 - M. Payer    - C2H6 emissions are too low. Use Yaping
!                              Xiao's C2H6 emissions instead.
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!
      ! Array for geographic mask
      REAL*8,  ALLOCATABLE, TARGET  :: EUROPE_MASK(:,:)

      ! Arrays for ground-based emissions
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_NOx(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_CO(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_SO2(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_NH3(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_ALK4(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_MEK(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_ALD2(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_PRPE(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_C2H6(:,:)

      ! Arrays for ship emissions
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_CO_SHIP(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_SO2_SHIP(:,:)
      REAL*8,  ALLOCATABLE, TARGET  :: EMEP_NOx_SHIP(:,:)

      ! Shadow variable for XNUMOL
      REAL*8,  ALLOCATABLE          :: s_XNUMOL(:)

      ! Shadow variables for logicals from Input_Opt
      LOGICAL                       :: LEMEPSHIP
      LOGICAL                       :: LFUTURE

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_europe_mask
!
! !DESCRIPTION: Function GET\_EUROPE\_MASK returns the value of the EUROPE 
!  mask for EMEP emissions at grid box (I,J).  MASK=1 if (I,J) is in the 
!  European region, or MASK=0 otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EUROPE_MASK( I, J ) RESULT( EUROPE )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I        ! Longitude index
      INTEGER, INTENT(IN) :: J        ! Latitude  index
!
! !RETURN VALUE:
! 
      REAL*8              :: EUROPE   ! Returns the mask value @ (I,J)
!
! !REVISION HISTORY: 
!  01 Nov 2005 - B. Field, R. Yantosca - Initial version
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

      END FUNCTION GET_EUROPE_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_emep_anthro
!
! !DESCRIPTION: Function GET\_EMEP\_ANTHRO returns the EMEP emission for 
!  GEOS-CHEM grid box (I,J) and tracer N.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_EMEP_ANTHRO( I, J, N, KG_S, SHIP ) RESULT( EMEP )
!
! !USES:
! 
      USE TRACERID_MOD, ONLY : IDTNO,   IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD, ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD, ONLY : IDTNH3,  IDTNO2
      USE GRID_MOD,     ONLY : GET_AREA_CM2
!
! !INPUT PARAMETERS:
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
!  01 Nov 2005 - B. Field, R. Yantosca - Initial version
!  (1 ) added SOx, SOx ship and NH3 emissions, plus optional kg/s output
!       (amv, 06/2008)
!  (2 ) Now returns ship emissions if requested (phs, 6/08)
!  (3 ) Added checks to avoid calling unavailable ship emissions (phs, 6/08)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  22 Mar 2012 - M. Payer    - C2H6 emissions are too low. Use Yaping
!                              Xiao's C2H6 emissions instead.
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
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
      HAS_SHIP = (/ IDTNO, IDTCO, IDTSO2 /)

      IF ( IS_SHIP .AND. .NOT. ANY( HAS_SHIP == N) ) THEN
         WRITE(6,*)'WARNING: EMEP SHIP emissions not available for'//
     $             'tracer #',N
         EMEP = 0D0
         RETURN
      ENDIF

      ! NOx
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( N  == IDTNO ) THEN
         IF ( IS_SHIP ) THEN
            EMEP = EMEP_NOx_SHIP(I,J)
         ELSE 
            EMEP = EMEP_NOx(I,J)
         ENDIF
!%%%%%%%%%%%%%KLUDGE TO EMITT SHIP NOX AS NOX %%%%%%%%%%%%%%
!         IF ( IS_SHIP ) THEN
!            EMEP = 0d0
!         ELSE 
!            EMEP = EMEP_NOx(I,J) + EMEP_NOx_SHIP(I,J)
!         ENDIF
!%%%%%%%%%%%%% END KLUDGE %%%%%%%%%%%%%%

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
         EMEP = -1d0

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

         IF ( N == IDTNO ) THEN
            ! Convert NOx using XNUMOL for NO2 since original NOx emissions
            ! are in mass units of NO2 (mpayer 4/18/13)
            EMEP = EMEP * GET_AREA_CM2( I, J, 1 ) / s_XNUMOL(IDTNO2)
         ELSE
            EMEP = EMEP * GET_AREA_CM2( I, J, 1 ) / s_XNUMOL(NN)
         ENDIF

      ENDIF

      END FUNCTION GET_EMEP_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_emep
!
! !DESCRIPTION: Subroutine EMISS\_EMEP reads the EMEP emission fields at 
!  1x1 resolution and regrids them to the current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EMEP( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
! 
      USE BPCH2_MOD,          ONLY : GET_TAU0, OPEN_BPCH2_FOR_READ
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE FILE_MOD,           ONLY : IOERROR
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE inquireMod,         ONLY : findFreeLUN
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE TIME_MOD,           ONLY : EXPAND_DATE,  GET_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR
#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNO2,  IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD,       ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD,       ONLY : IDTNH3,  IDTNO
      USE GRID_MOD,           ONLY : GET_AREA_CM2
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
!
! !REVISION HISTORY: 
!  01 Nov 2005 - B. Field, R. Yantosca - Initial version
!  (1 ) Modified for IPCC future emissions.  Now references LFUTURE from
!        "logical_mod.f". (bmy, 5/30/06)
!  13 Mar 2012 - M. Cooper   - Changed gridding algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
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
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: EMEP_NYMD, EMEP_YEAR
      REAL*8              :: EMEP_TAU,  TAU0
      CHARACTER(LEN=255)  :: FILENAME

      ! For bpch file format
      INTEGER             :: I,  IU_FILE, J,  L,  N,  IOS
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      INTEGER             :: NI,        NJ,        NL
      INTEGER             :: IFIRST,    JFIRST,    LFIRST
      INTEGER             :: SCALEYEAR
      REAL*4              :: ARRAY(I1x1,J1x1,1)
      REAL*4              :: LONRES,    LATRES
      REAL*4              :: Sc(IIPAR,JJPAR)
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED
      CHARACTER(LEN=255)  :: LLFILENAME
      REAL*8              :: INGRID (I1x1, J1x1 )

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      REAL*8              :: A(IIPAR,JJPAR)
#endif

      ! For fields from Input_Opt
      INTEGER             :: N_TRACERS
      REAL*8              :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! EMISS_EMEP begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LEMEPSHIP = Input_Opt%LEMEPSHIP
      LFUTURE   = Input_Opt%LFUTURE
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! First-time initialization
      IF ( FIRST ) THEN
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
         A = 0d0
#endif
         
         CALL INIT_EMEP( am_I_Root, Input_Opt, RC )
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
         
      ! Find a free file LUN
      IU_FILE = findFreeLun()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME ), IU_FILE
 100  FORMAT( '     - EMISS_EMEP: Reading ', a, ' on unit ', i4 )

      !=================================================================
      ! Read data at 1x1 resolution and regrid to current grid size
      !=================================================================

      ! Open file
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

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
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_NOx,  IS_MASS=0,
     &                             netCDF=.TRUE.                     )
            ! CO [molec/cm2/s]
            CASE( 4  )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_CO,   IS_MASS=0,
     &                             netCDF=.TRUE.                     )
               
            ! ALK4 [atoms C/cm2/s]
            CASE( 5  )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_ALK4, IS_MASS=0,
     &                             netCDF=.TRUE.                     )

            ! MEK [atoms C/cm2/s]
            CASE( 10 )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_MEK,  IS_MASS=0,
     &                             netCDF=.TRUE.                     )

            ! ALD2 [atoms C/cm2/s]
            CASE( 11 )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_ALD2, IS_MASS=0,
     &                             netCDF=.TRUE.                     )

            ! PRPE [atoms C/cm2/s]
            CASE( 18 )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_PRPE, IS_MASS=0,
     &                             netCDF=.TRUE.                     )

            ! C2H6 [atoms C/cm2/s]
            CASE( 21 )
               INGRID = ARRAY(:,:,1)
               CALL DO_REGRID_A2A( LLFILENAME, I1x1,      J1x1,
     &                             INGRID,     EMEP_C2H6, IS_MASS=0,
     &                             netCDF=.TRUE.                     )

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

      IF ( SCALEYEAR > 1989) THEN

         ! new EMEP data is only defined from 1990-2007
         EMEP_YEAR = MIN( SCALEYEAR, 2007 )

         CALL READ_EMEP_UPDATED(  1, EMEP_YEAR, EMEP_NOx, 0 )
         IF ( SCALEYEAR >2007 .AND. SCALEYEAR <=2010 ) THEN
         CALL GET_ANNUAL_SCALAR( 71, 2007, SCALEYEAR, Sc )
         EMEP_NOx(:,:) = EMEP_NOx(:,:) * Sc(:,:)
         END IF

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

#if defined( DEVEL )
      DO I=1,IIPAR
      DO J=1,JJPAR
         A(I,J)  = GET_AREA_CM2( I, J, 1 ) 
      ENDDO
      ENDDO

      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTNO) = EMEP_NOx_SHIP(:,:)
     &        * A / XNUMOL(IDTNO2)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTNO) = EMEP_NOx(:,:)
     &        * A / XNUMOL(IDTNO2)
      ENDIF

      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTCO) = EMEP_CO_SHIP(:,:)
     &        * A / XNUMOL(IDTCO)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTCO) = EMEP_CO(:,:)
     &        * A / XNUMOL(IDTCO)
      ENDIF

      ! ALK4 (>= C4 alkanes)
      State_Chm%TRAC_TEND(:,:,1,IDTALK4) = EMEP_ALK4(:,:)
     &        * A / XNUMOL(IDTALK4)

      ! MEK
      State_Chm%TRAC_TEND(:,:,1,IDTMEK) = EMEP_MEK(:,:)
     &        * A / XNUMOL(IDTMEK)

      ! ALD2 (acetaldehyde)
      State_Chm%TRAC_TEND(:,:,1,IDTALD2) = EMEP_ALD2(:,:)
     &        * A / XNUMOL(IDTALD2)

      ! PRPE (>= C3 alkenes)
      State_Chm%TRAC_TEND(:,:,1,IDTPRPE) = EMEP_PRPE(:,:)
     &        * A / XNUMOL(IDTPRPE)

      ! C2H6 
      State_Chm%TRAC_TEND(:,:,1,IDTC2H6) = EMEP_C2H6(:,:)
     &        * A / XNUMOL(IDTC2H6)
         
      ! SO2
      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTSO2) = EMEP_SO2_SHIP(:,:)
     &        * A / XNUMOL(IDTSO2)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTSO2) = EMEP_SO2(:,:)
     &        * A / XNUMOL(IDTSO2)
      ENDIF

      ! NH3
      State_Chm%TRAC_TEND(:,:,1,IDTNH3) = EMEP_NH3(:,:)
     &        * A / XNUMOL(IDTNH3)
#endif

      !=================================================================
      ! Print emission totals
      !=================================================================

      ! Print totals for EMEP_YEAR
      CALL TOTAL_ANTHRO_TG( EMEP_YEAR,   SCALEYEAR, 
     &                      GET_MONTH(), Input_Opt  )

      END SUBROUTINE EMISS_EMEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_emep_05x0666
!
! !DESCRIPTION: Subroutine EMISS\_EMEP reads the EMEP emission fields at
!  05x0666 resolution and regrids them to the current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_EMEP_05x0666( am_I_Root, Input_Opt,
     &                               State_Chm, RC         )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_TAU0,     READ_BPCH2
      USE CMN_O3_MOD
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : EXPAND_DATE,  GET_YEAR
      USE TIME_MOD,           ONLY : GET_MONTH
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR_05x0666_NESTED
#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNO2,  IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD,       ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD,       ONLY : IDTNH3,  IDTNO
      USE GRID_MOD,           ONLY : GET_AREA_CM2
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
!  23 Oct 2006 - A. v. Donkelaar - Initial version, modified from EMISS_EMEP
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: EMEP_NYMD, EMEP_YEAR
      REAL*8              :: EMEP_TAU,  TAU0
      CHARACTER(LEN=255)  :: FILENAME

      ! For bpch file format
      INTEGER             :: I,  J,  L,  N,  IOS
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      INTEGER             :: NI,        NJ,        NL
      INTEGER             :: IFIRST,    JFIRST,    LFIRST
      INTEGER             :: SCALEYEAR
      REAL*4              :: ARRAY(IIPAR,JJPAR,1)
      REAL*4              :: LONRES,    LATRES
      REAL*4              :: Sc(IIPAR,JJPAR)
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED

#if defined( DEVEL )
      REAL*8, ALLOCATABLE :: A(:,:)
      INTEGER AS
#endif

      ! For fields from Input_Opt
      INTEGER             :: N_TRACERS
      REAL*8              :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! EMISS_EMEP begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! First-time initialization
      IF ( FIRST ) THEN
#if defined( DEVEL )
         ALLOCATE( A( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMEP:A' )
         A = 0d0
#endif         
         
         CALL INIT_EMEP( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      ! 1x1 file name for EMEP 2000
      FILENAME  = TRIM( DATA_DIR ) //
     &            'EMEP_200510/EMEP.geos.05x0666.YYYY'

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
 100  FORMAT( '     - EMISS_EMEP_05x0666: Reading ', a )

      !=================================================================
      ! Read data at 05x0666 resolution
      !=================================================================

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, EMEP_TAU, 
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_CO(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_NOx(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',18, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_PRPE(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 5, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_ALK4(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',21, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_C2H6(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',11, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_ALD2(:,:) = ARRAY(:,:,1)

      CALL READ_BPCH2( FILENAME, 'ANTHSRCE',10, EMEP_TAU,
     &          IIPAR, JJPAR, 1, ARRAY, QUIET=.TRUE.)
      EMEP_MEK(:,:) = ARRAY(:,:,1)

      !=================================================================
      ! Get and apply annual emissions factors (amv, phs, 3/17/08)
      !=================================================================

      !=================================================================
      ! If we are at or above 1990, can apply updated EMEP emissions for
      ! NOx, CO, NH3 and include SOx (amv, 06/04/08)
      !=================================================================

      IF ( SCALEYEAR > 1989 ) THEN

         ! new EMEP data is only defined from 1990-2007
         EMEP_YEAR = MIN( SCALEYEAR, 2007 )

         CALL READ_EMEP_UPDATED_05x0666(  1, EMEP_YEAR, EMEP_NOx, 0 )
         CALL READ_EMEP_UPDATED_05x0666(  4, EMEP_YEAR, EMEP_CO, 0 )
         CALL READ_EMEP_UPDATED_05x0666( 26, EMEP_YEAR, EMEP_SO2, 0 )
         CALL READ_EMEP_UPDATED_05x0666( 30, EMEP_YEAR, EMEP_NH3, 1 )

         CALL READ_EMEP_UPDATED_05x0666( 1,EMEP_YEAR, EMEP_NOx_SHIP, 2)
         CALL READ_EMEP_UPDATED_05x0666( 4,EMEP_YEAR, EMEP_CO_SHIP, 2)
         CALL READ_EMEP_UPDATED_05x0666( 26,EMEP_YEAR, EMEP_SO2_SHIP, 2)

      ! Need to use for SOx/NH3 anyways, but SOx scale back further
      ELSE

         CALL READ_EMEP_UPDATED_05x0666( 26, 1990, EMEP_SO2, 0 )
         CALL READ_EMEP_UPDATED_05x0666( 26, 1990, EMEP_SO2_SHIP, 2 )
         CALL READ_EMEP_UPDATED_05x0666( 30, 1990, EMEP_NH3, 1 )

         CALL GET_ANNUAL_SCALAR_05x0666_NESTED(73,1990,SCALEYEAR,Sc)
         EMEP_SO2(:,:) = EMEP_SO2(:,:) * Sc(:,:)
!         EMEP_SO2_SHIP = EMEP_SO2_SHIP * Sc  ! do not scale SHIP
      ENDIF


      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================
      IF ( LFUTURE ) THEN
         CALL EMEP_SCALE_FUTURE
      ENDIF

#if defined( DEVEL )
      DO I=1,IIPAR
      DO J=1,JJPAR
         A(I,J)  = GET_AREA_CM2( I, J, 1 ) 
      ENDDO
      ENDDO

      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTNO)  = EMEP_NOx_SHIP(:,:)
     &        * A / XNUMOL(IDTNO2)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTNO)  = EMEP_NOx(:,:)
     &        * A / XNUMOL(IDTNO2)
      ENDIF

      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTCO) = EMEP_CO_SHIP(:,:)
     &        * A / XNUMOL(IDTCO)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTCO) = EMEP_CO(:,:)
     &        * A / XNUMOL(IDTCO)
      ENDIF

      ! ALK4 (>= C4 alkanes)
      State_Chm%TRAC_TEND(:,:,1,IDTALK4) = EMEP_ALK4(:,:)
     &        * A / XNUMOL(IDTALK4)

      ! MEK
      State_Chm%TRAC_TEND(:,:,1,IDTMEK) = EMEP_MEK(:,:)
     &        * A / XNUMOL(IDTMEK)

      ! ALD2 (acetaldehyde)
      State_Chm%TRAC_TEND(:,:,1,IDTALD2) = EMEP_ALD2(:,:)
     &        * A / XNUMOL(IDTALD2)

      ! PRPE (>= C3 alkenes)
      State_Chm%TRAC_TEND(:,:,1,IDTPRPE) = EMEP_PRPE(:,:)
     &        * A / XNUMOL(IDTPRPE)

      ! C2H6 
      State_Chm%TRAC_TEND(:,:,1,IDTC2H6) = EMEP_C2H6(:,:)
     &        * A / XNUMOL(IDTC2H6)
         
      ! SO2
      IF ( LEMEPSHIP ) THEN
         State_Chm%TRAC_TEND(:,:,1,IDTSO2) = EMEP_SO2_SHIP(:,:)
     &        * A / XNUMOL(IDTSO2)
      ELSE 
         State_Chm%TRAC_TEND(:,:,1,IDTSO2) = EMEP_SO2(:,:)
     &        * A / XNUMOL(IDTSO2)
      ENDIF

      ! NH3
      State_Chm%TRAC_TEND(:,:,1,IDTNH3) = EMEP_NH3(:,:)
     &        * A / XNUMOL(IDTNH3)
#endif

      !=================================================================
      ! Print emission totals
      !=================================================================

      ! Print totals for EMEP_YEAR
      CALL TOTAL_ANTHRO_TG( EMEP_YEAR,   SCALEYEAR, 
     &                      GET_MONTH(), Input_Opt )

      END SUBROUTINE EMISS_EMEP_05x0666
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emep_scale_future
!
! !DESCRIPTION: Subroutine EMEP\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the EMEP anthropogenic emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMEP_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_ALK4ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_C2H6ff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_PRPEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_TONEff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_VOCff

      USE CMN_SIZE_MOD             ! Size parameters
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

      END SUBROUTINE EMEP_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the amount of EMEP 
!  anthropogenic emissions that are emitted each month in Tg or Tg C. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( EMEP_YEAR,  EMISS_YEAR, 
     &                            EMEP_MONTH, Input_Opt )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TIME_MOD,           ONLY : ITS_A_LEAPYEAR
      USE TRACERID_MOD,       ONLY : IDTNO,   IDTCO,   IDTALK4, IDTMEK
      USE TRACERID_MOD,       ONLY : IDTALD2, IDTPRPE, IDTC2H6, IDTSO2
      USE TRACERID_MOD,       ONLY : IDTNH3
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: EMEP_YEAR    ! EMEP base year
      INTEGER,        INTENT(IN) :: EMISS_YEAR   ! Current simulated year
      INTEGER,        INTENT(IN) :: EMEP_MONTH   ! Current simulated month
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  10 Nov 2004 - R. Hudman, R. Yantosca - Initial version  
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now replace FMOL with TRACER_MW_KG (bmy, 10/25/05) 
!  (3 ) Now only print totals of defined tracers; other totals will be
!        printed as zeroes. (bmy, 2/6/06)
!  (4 ) Now emissions and base year are arguments. Output in Tg/month
!        since this is called monthly (phs, 12/9/08)
!  (5 ) Bug fix, now print out correct monthly EMEP totals (bmy, 1/30/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  22 Mar 2012 - M. Payer    - Remove print for C2H6 emissions.
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J
      REAL*8           :: A,   B(9), NOX,  CO,  ALK4
      REAL*8           :: MEK, ALD2, PRPE, C2H6, SO2
      REAL*8           :: NH3
      CHARACTER(LEN=3) :: UNIT

      ! Days per month
      REAL*8           :: DAYS_IN_MONTH
      REAL*8           :: DMON(12) = (/ 31d0, 28d0, 31d0, 30d0,
     &                                 31d0, 30d0, 31d0, 31d0,
     &                                 30d0, 31d0, 30d0, 31d0 /) 

      ! For fields from Input_Opt
      LOGICAL          :: LEMEPSHIP
      INTEGER          :: N_TRACERS
      REAL*8           :: XNUMOL(Input_Opt%N_TRACERS)

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)


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
      
      ! Get the proper # of days in the month for totaling
      IF ( EMEP_MONTH == 2 .and. ITS_A_LEAPYEAR( EMISS_YEAR ) ) THEN
         DAYS_IN_MONTH = DMON(EMEP_MONTH) + 1
      ELSE
         DAYS_IN_MONTH = DMON(EMEP_MONTH)
      ENDIF
      
      ! Define conversion factors for kg/molec
      ! (Undefined tracers will be zero)
      B(:) = 0d0
      IF ( IDTNO   > 0 ) B(1) = 14d-3 / 6.0225d23        ! Tg N
      IF ( IDTCO   > 0 ) B(2) = 1d0   / XNUMOL(IDTCO  )
      IF ( IDTALK4 > 0 ) B(3) = 1d0   / XNUMOL(IDTALK4)
      IF ( IDTMEK  > 0 ) B(4) = 1d0   / XNUMOL(IDTMEK )
      IF ( IDTALD2 > 0 ) B(5) = 1d0   / XNUMOL(IDTALD2)
      IF ( IDTPRPE > 0 ) B(6) = 1d0   / XNUMOL(IDTPRPE)
      IF ( IDTC2H6 > 0 ) B(7) = 1d0   / XNUMOL(IDTC2H6)
      IF ( IDTSO2  > 0 ) B(8) = 32d-3 / 6.0225d23        ! Tg S
      IF ( IDTNH3  > 0 ) B(9) = 1d0   / XNUMOL(IDTNH3)

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

      ! Loop over grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Surface area [cm2] * seconds in this year 
         ! Multiply by 1d-9 to convert from [kg] to [Tg]
         A = GET_AREA_CM2( I, J, 1 ) * DAYS_IN_MONTH * 86400d0 * 1d-9 
         
         ! Sum emissions (list NOx as Tg N)
         NOX  = NOX  + ( EMEP_NOX (I,J) + EMEP_NOX_SHIP(I,J) )
     &        * A * B(1)
         CO   = CO   + ( EMEP_CO  (I,J) + EMEP_CO_SHIP(I,J)  )
     &        * A * B(2) 
         SO2  = SO2  + ( EMEP_SO2 (I,J) + EMEP_SO2_SHIP(I,J) )
     &        * A * B(8) 

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
      WRITE( 6, 110   ) 'NO  ', EMISS_YEAR, EMEP_MONTH, NOx,  ' N'
      WRITE( 6, 110   ) 'CO  ', EMISS_YEAR, EMEP_MONTH, CO,   '  '
      WRITE( 6, 110   ) 'SO2 ', EMISS_YEAR, EMEP_MONTH, SO2,  ' S'
      WRITE( 6, 110   ) 'NH3 ', EMISS_YEAR, EMEP_MONTH, NH3,  '  '
      WRITE( 6, 110   ) 'ALK4', EMISS_YEAR, EMEP_MONTH, ALK4, ' C'
      WRITE( 6, 110   ) 'MEK ', EMISS_YEAR, EMEP_MONTH, MEK,  ' C'
      WRITE( 6, 110   ) 'ALD2', EMISS_YEAR, EMEP_MONTH, ALD2, ' C'
      WRITE( 6, 110   ) 'PRPE', EMISS_YEAR, EMEP_MONTH, PRPE, ' C'
 110  FORMAT( 'EMEP anthropogenic ', a4, ' for ', i4, '/', i2.2,
     &        ': ', f13.6, ' Tg', a2 )

      WRITE( 6, '(/,a,/)' ) 'EMEP_MOD: EMEP C2H6 anthro '       //
     &   'emissions are too low. Using offline C2H6 emissions ' //
     &   'from Yaping Xiao.'

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_europe_mask
!
! !DESCRIPTION: Subroutine READ\_EUROPE\_MASK reads and regrids the 
!  Europe mask for the EMEP anthropogenic emissions.
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
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE CMN_SIZE_MOD         ! Size parameters
!
! !REVISION HISTORY: 
!  18 Oct 2006 - R. Yantosca - Initial version
!  (1 ) Now read the Europe mask from a disk file instead of defining it as 
!        a rectangular box (bmy, 10/18/06)
!  (2 ) Updated the mask file to correspond with the 200911 EMEP emissions
!        (gvinken, 11/24/10)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  15 Aug 2012 - M. Payer    - Fixed minor bugs in regridding of mask; Also set
!                              mask to 1 if greater than 0 (L. Murray)
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY(I1x1,J1x1,1)
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8                  :: INGRID (I1x1,J1x1)
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)

      !=================================================================
      ! READ_EUROPE_MASK begins here!
      !=================================================================

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'EMEP_200911/EMEP_mask.geos.1x1'

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_EUROPE_MASK: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 0d0,       I1x1,     J1x1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 


      ! Regrid from GEOS 1x1 GRID to current model resolution
      INGRID = ARRAY(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )  

      EUROPE_MASK = OUTGRID 

      WHERE ( EUROPE_MASK > 0D0 ) EUROPE_MASK = 1D0

      END SUBROUTINE READ_EUROPE_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_europe_mask_05x0666
!
! !DESCRIPTION: Subroutine READ\_EUROPE\_MASK reads and regrids the
!  Europe mask for the EMEP anthropogenic emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EUROPE_MASK_05x0666
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR

      USE CMN_SIZE_MOD       ! Size parameters
!
! !REVISION HISTORY:
!  18 Oct 2006 - R. Yantosca - Initial version
!  (1 ) Now read the Europe mask from a disk file instead of defining it as
!        a rectangular box (bmy, 10/18/06)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_EUROPE_MASK begins here!
      !=================================================================

      ! File name
      FILENAME  = TRIM( DATA_DIR ) //
     &            'EMEP_200510/EMEP_mask.geos.05x0666'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_EUROPE_MASK: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 0d0,       IIPAR,     JJPAR,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      EUROPE_MASK(:,:) = ARRAY(:,:,1)

      END SUBROUTINE READ_EUROPE_MASK_05x0666
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_emep_updated
!
! !DESCRIPTION: Subroutine READ\_EMEP\_UPDATED reads updated EMEP emissions 
!  from the year 1990 including SOx emissions.  These are regridded to the 
!  simulation resolution. Ship emissions can also be included. 
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
      USE GRID_MOD,       ONLY : GET_AREA_CM2
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters
      USE CMN_O3_MOD           ! SCALEYEAR
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: TRACER              ! Tracer number
      INTEGER, INTENT(IN)  :: EMEP_YEAR           ! Year of emissions to read
      INTEGER, INTENT(IN)  :: wSHIP               ! Use ground, ship, or both?
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: ARRAY(IIPAR,JJPAR)  ! Output array
!
! !REVISION HISTORY: 
!  28 Jan 2009 - A. v. Donkelaar, P. Le Sager - Initial version
!  28 Jan 2009 - P. Le Sager - Now account for LEMEPSHIP
!  29 Oct 2009 - Added multi-species seasonality (amv)
!  04 Jan 2010 - Extended to 2007, changed input format (amv)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
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
      REAL*8                        :: EMEP_TAU, TAU, A, B
      INTEGER                       :: EMEP_NYMD, MN, RATIOID, I, J
      CHARACTER(LEN=255)            :: LLFILENAME
      REAL*8                        :: INGRID (I1x1, J1x1 )
      REAL*8                        :: OUTGRID(IIPAR,JJPAR)

      ARRAY_1x1_SHIP(:,:,:) = 0.d0
      ARRAY_1x1_LAND(:,:,:) = 0.d0

      ! YYYYMMDD value for 1st day of EMEP_YEAR
      EMEP_NYMD = ( EMEP_YEAR * 10000 ) + 0101

      ! TAU0 value corresponding to EMEP_NYMD
      EMEP_TAU  = GET_TAU0( 1, 1, EMEP_YEAR )

      ! Expand filename
      DIR = TRIM( DATA_DIR_1x1 ) // 'EMEP_200911/'

      ! wSHIP = 0 means no ship emissions included
      ! wSHIP = 1 means include ships emissions
      ! wSHIP = 2 means only ship emissions

      IF ( wSHIP .lt. 2 ) THEN

         FILENAME = TRIM(DIR) //  
     &              'EMEP-YYYY.geos.1x1' 
 
         CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )
 
         WRITE( 6, 100 ) TRIM( FILENAME ) 
 100     FORMAT( '     - READ_EMEP_UPDATED: Reading ', a ) 
 
         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',      TRACER,
     &                    EMEP_TAU,  I1x1,           J1x1,
     &                    1,         ARRAY_1x1_LAND, QUIET=.TRUE. )
 
      ENDIF

      IF ( ( wSHIP .gt. 0 ) .AND. LEMEPSHIP ) THEN

         FILENAME = TRIM(DIR) //
     &              'EMEP-SHIP-YYYY.geos.1x1'

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
      IF ( wSHIP .lt. 2 ) THEN

         ! Apply Monthly Factors over land
         TAU = GET_TAU0( GET_MONTH(), 1, 2005)

         ! Select proper filename & tracer #
         SELECT CASE( TRACER )

            ! NOx
            CASE( 1 )
               RATIOID = 71
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'NOx-EMEP-SeasonalScalar.geos.1x1'

            ! CO
            CASE( 4 )
               RATIOID = 72
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'CO-EMEP-SeasonalScalar.geos.1x1'

            ! SO2
            CASE( 26 )
               RATIOID = 73
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'SOx-EMEP-SeasonalScalar.geos.1x1'

            ! SO2
            CASE( 30 )
               RATIOID = 74
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'NH3-EMEP-SeasonalScalar.geos.1x1'

         END SELECT

         ! Echo info
         WRITE( 6, 101 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', RATIOID,
     &                    TAU,      I1x1,       J1x1,
     &                    1,        ARRAY_1x1,  QUIET=.TRUE. )

         ARRAY_1x1_LAND(:,:,1) = ARRAY_1x1_LAND(:,:,1) 
     &                           * ARRAY_1x1(:,:,1)

      ENDIF

      IF ( wSHIP .eq. 0 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_LAND(:,:,1)
      IF ( wSHIP .eq. 1 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_LAND(:,:,1) + 
     &                                       ARRAY_1x1_SHIP(:,:,1)
      IF ( wSHIP .eq. 2 ) ARRAY_1x1(:,:,1) = ARRAY_1x1_SHIP(:,:,1)

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Regrid to current model resolution [kg/yr]
      INGRID = ARRAY_1x1(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                    INGRID,     OUTGRID, IS_MASS=1,
     &                    netCDF=.TRUE.                   )
      ARRAY  = OUTGRID

      ! Convert SOx to SO2 assuming a SOx is 95% SO2 over Europe, as used
      ! throughout GEOS-Chem, and as per Chin et al, 2000
      IF ( TRACER .eq. 26 ) ARRAY(:,:) = ARRAY(:,:) * 0.95d0

      ! convert to molec/cm2 for consistency with previous
      ! emissions
      B = 0d0
      IF ( TRACER .eq.  1 ) B = 1.d3 / 46d0 * 6.0225d23
      IF ( TRACER .eq.  4 ) B = 1.d3 / 28d0 * 6.0225d23
      IF ( TRACER .eq. 26 ) B = 1.d3 / 64d0 * 6.0225d23
      IF ( TRACER .eq. 30 ) B = 1.d3 / 17d0 * 6.0225d23

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Surface area [cm2] * sec per year
         A = GET_AREA_CM2( I, J, 1 ) * 365d0 * 86400d0

         ! Unit conversion
         ARRAY(I,J) = ARRAY(I,J) / A * B

      ENDDO
      ENDDO

      END SUBROUTINE READ_EMEP_UPDATED
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_emep_updated_05x0666
!
! !DESCRIPTION: Subroutine READ\_EMEP\_UPDATED reads updated EMEP emissions
!  from the year 1990 including SOx emissions.  These are regridded to the
!  simulation resolution. Ship emissions can also be included.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_EMEP_UPDATED_05x0666( TRACER, EMEP_YEAR, ARRAY, 
     &                   wSHIP )
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : READ_BPCH2, GET_TAU0
      USE TIME_MOD,       ONLY : EXPAND_DATE, GET_MONTH
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE GRID_MOD,       ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD
      USE CMN_O3_MOD
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: TRACER              ! Tracer number
      INTEGER, INTENT(IN)  :: EMEP_YEAR           ! Year of emissions to read
      INTEGER, INTENT(IN)  :: wSHIP               ! Use ground, ship, or both?
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: ARRAY(IIPAR,JJPAR)  ! Output array
!
! !REVISION HISTORY:
!  28 Jan 2009 - A. v. Donkelaar, P. Le Sager - Initial version
!  28 Jan 2009 - P. Le Sager - Now account for LEMEPSHIP
!  29 Oct 2009 - Added multi-species seasonality (amv)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8               :: ARRAY_05x0666(IIPAR,JJPAR,1)
      REAL*4               :: ARRAY_05x0666_R4(IIPAR,JJPAR,1)
      REAL*4               :: ARRAY_05x0666_SHIP(IIPAR,JJPAR,1)
      REAL*4               :: ARRAY_05x0666_LAND(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)   :: FILENAME, DIR
      REAL*8               :: EMEP_TAU, TAU, A, B
      INTEGER              :: EMEP_NYMD, MN, RATIOID, I, J
      CHARACTER(LEN=2)     :: SMN
      CHARACTER(LEN=1)     :: SSMN

      ARRAY_05x0666_SHIP(:,:,:) = 0.d0
      ARRAY_05x0666_LAND(:,:,:) = 0.d0

      ! YYYYMMDD value for 1st day of EMEP_YEAR
      EMEP_NYMD = ( EMEP_YEAR * 10000 ) + 0101

      ! TAU0 value corresponding to EMEP_NYMD
      EMEP_TAU  = GET_TAU0( 1, 1, EMEP_YEAR )

      ! Expand filename
      DIR = TRIM( DATA_DIR ) // 'EMEP_200911/'

      ! wSHIP = 0 means no ship emissions included
      ! wSHIP = 1 means include ships emissions
      ! wSHIP = 2 means only ship emissions

      IF ( wSHIP .lt. 2 ) THEN

         FILENAME = TRIM(DIR) // 
     &              'EMEP-YYYY.1p2x2p3.eu.bpch'

         CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )

         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - READ_EMEP_UPDATED_05x0666
     &               : Reading ', a )

         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',      TRACER,
     &                    EMEP_TAU,  IIPAR,       JJPAR,
     &                    1,      ARRAY_05x0666_LAND, QUIET=.TRUE. )

      ENDIF

      IF ( ( wSHIP .gt. 0 ) .AND. LEMEPSHIP ) THEN

         FILENAME = TRIM(DIR) //
     &              'EMEP-SHIP-YYYY.1p2x2p3.eu.bpch'

         CALL EXPAND_DATE( FILENAME, EMEP_NYMD, 000000 )

         WRITE( 6, 101 ) TRIM( FILENAME )
 101     FORMAT( '     - READ_EMEP_UPDATED_05x0666
     &         : Reading ', a )

         CALL READ_BPCH2( FILENAME, 'ANTHSRCE',      TRACER,
     &                    EMEP_TAU,  IIPAR,       JJPAR,
     &                    1,      ARRAY_05x0666_SHIP, QUIET=.TRUE. )

      ENDIF

      ! Apply monthly variation (courtesy of the GENEMIS project
      ! coordinated by the Institute of Energy Economics and the
      ! Rational Use of Energy (IER) at the University of
      ! Stuttgart) (amv, 11/24/2008)

      ! Expand filename
      DIR = TRIM( DATA_DIR ) // 'EMEP_200806/'

      IF ( wSHIP .lt. 2 ) THEN

         ! Apply Monthly Factors over land
         TAU = GET_TAU0( GET_MONTH(), 1, 2005)

         ! Select proper filename & tracer #
         SELECT CASE( TRACER )

            ! NOx
            CASE( 1 )
               RATIOID = 71
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'NOx-EMEP-SeasonalScalar.geos.05x0666'

            ! CO
            CASE( 4 )
               RATIOID = 72
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'CO-EMEP-SeasonalScalar.geos.05x0666'

            ! SO2
            CASE( 26 )
               RATIOID = 73
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'SOx-EMEP-SeasonalScalar.geos.05x0666'

            ! SO2
            CASE( 30 )
               RATIOID = 74
               FILENAME = TRIM( DIR )          //
     &                    'SeasonalVariation/' //
     &                    'NH3-EMEP-SeasonalScalar.geos.05x0666'

         END SELECT

         ! Echo info
         WRITE( 6, 101 ) TRIM( FILENAME )

         ! Read data
         CALL READ_BPCH2( FILENAME, 'RATIO-2D', RATIOID,
     &                    TAU,   IIPAR,      JJPAR,
     &                    1,    ARRAY_05x0666_R4,     QUIET=.TRUE. )

         ARRAY_05x0666_LAND(:,:,1) = ARRAY_05x0666_LAND(:,:,1)
     &                           * ARRAY_05x0666_R4(:,:,1)

      ENDIF

      IF ( wSHIP .eq. 0 ) ARRAY_05x0666(:,:,1) = 
     &            ARRAY_05x0666_LAND(:,:,1)
      IF ( wSHIP .eq. 1 ) ARRAY_05x0666(:,:,1) = 
     &            ARRAY_05x0666_LAND(:,:,1) + ARRAY_05x0666_SHIP(:,:,1)
      IF ( wSHIP .eq. 2 ) ARRAY_05x0666(:,:,1) = 
     &            ARRAY_05x0666_SHIP(:,:,1)

      ARRAY(:,:) = ARRAY_05x0666(:,:,1)

      ! Convert SOx to SO2 assuming a SOx is 95% SO2 over Europe, as used
      ! throughout GEOS-Chem, and as per Chin et al, 2000
      IF ( TRACER .eq. 26 ) ARRAY(:,:) = ARRAY(:,:) * 0.95d0

      ! convert to molec/cm2 for consistency with previous
      ! emissions
      B = 0d0
      IF ( TRACER .eq.  1 ) B = 1.d3 / 46d0 * 6.0225d23
      IF ( TRACER .eq.  4 ) B = 1.d3 / 28d0 * 6.0225d23
      IF ( TRACER .eq. 26 ) B = 1.d3 / 64d0 * 6.0225d23
      IF ( TRACER .eq. 30 ) B = 1.d3 / 17d0 * 6.0225d23

      ! Loop over surface boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Surface area [cm2] * sec per year
         A = GET_AREA_CM2( I, J, 1 ) * 365d0 * 86400d0

         ! Unit conversion
         ARRAY(I,J) = ARRAY(I,J) / A * B

      ENDDO
      ENDDO

      END SUBROUTINE READ_EMEP_UPDATED_05x0666
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_emep
!
! !DESCRIPTION: Subroutine INIT\_EMEP allocates and zeroes EMEP module 
!  arrays, and also creates the mask which defines the European region.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_EMEP( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
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
!  01 Nov 2005 - B. Field, R. Yantosca - Initial version
!  (1 ) Now call READ_EUROPE_MASK to read & regrid EUROPE_MASK from disk 
!        instead of just defining it as a rectangular box. (bmy, 10/18/06)
!  26 Jan 2010 - R. Yantosca - Fixed cut-n-paste error.  Now make sure to zero 
!                              EMEP_CO_SHIP and EMEP_NOx_SHIP. 
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC 
!  04 Apr 2013 - R. Yantosca - Now allocate a shadow variable of XNUMOL
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! INIT_EMEP begins here!
      !=================================================================

      ! Return if LEMEP is false
      IF ( .not. Input_Opt%LEMEP ) RETURN
      
      !--------------------------------
      ! Allocate and zero arrays
      !--------------------------------
      ALLOCATE( EMEP_NOx( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_NOx' )
      EMEP_NOx = 0d0

      ALLOCATE( EMEP_CO( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_CO' )
      EMEP_CO = 0d0

      ALLOCATE( EMEP_SO2( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2' )
      EMEP_SO2 = 0d0

      ALLOCATE( EMEP_SO2_SHIP( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_SO2_SHIP' )
      EMEP_SO2_SHIP = 0d0

      ALLOCATE( EMEP_CO_SHIP( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_CO_SHIP' )
      EMEP_CO_SHIP = 0d0

      ALLOCATE( EMEP_NOx_SHIP( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_NOx_SHIP' )
      EMEP_NOx_SHIP = 0d0

      ALLOCATE( EMEP_NH3( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_NH3' )
      EMEP_NH3 = 0d0

      ALLOCATE( EMEP_ALK4( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_ALK4' )
      EMEP_ALK4 = 0d0

      ALLOCATE( EMEP_MEK( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_MEK' )
      EMEP_MEK = 0d0

      ALLOCATE( EMEP_ALD2( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_ALD2' )
      EMEP_ALD2 = 0d0

      ALLOCATE( EMEP_PRPE( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_PRPE' )
      EMEP_PRPE = 0d0

      ALLOCATE( EMEP_C2H6( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EMEP_C2H6' )
      EMEP_C2H6 = 0d0

      ALLOCATE( EUROPE_MASK( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'EUROPE_MASK' )
      EUROPE_MASK = 0d0

      ! Shadow variable for XNUMOL
      ALLOCATE( s_XNUMOL( Input_Opt%N_TRACERS ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'p_XNUMOL' )
      s_XNUMOL = Input_Opt%XNUMOL(1:Input_Opt%N_TRACERS)

      ! Read and regrid the European mask
#if   defined(GRID05x0666)
         CALL READ_EUROPE_MASK_05x0666
#else
         CALL READ_EUROPE_MASK
#endif

      END SUBROUTINE INIT_EMEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_emep
!
! !DESCRIPTION: Subroutine CLEANUP\_EMEP deallocates all module arrays.
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

      END SUBROUTINE CLEANUP_EMEP
!EOC
      END MODULE EMEP_MOD
