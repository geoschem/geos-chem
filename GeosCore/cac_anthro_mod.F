!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cac_anthro_mod
!
! !DESCRIPTION: Module CAC\_ANTHRO\_MOD contains variables and routines to 
!  read the  Criteria Air Contaminant Canadian anthropogenic emissions 
!  (amv, phs, 1/28/2009)
!\\
!\\
! !INTERFACE: 
!
      MODULE CAC_ANTHRO_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_CAC_ANTHRO
      PUBLIC :: EMISS_CAC_ANTHRO
      PUBLIC :: EMISS_CAC_ANTHRO_05x0666
      PUBLIC :: GET_CANADA_MASK
      PUBLIC :: GET_CAC_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: CAC_SCALE_FUTURE
      PRIVATE :: READ_CANADA_MASK
      PRIVATE :: READ_CANADA_MASK_05x0666
      PRIVATE :: INIT_CAC_ANTHRO
      PRIVATE :: TOTAL_ANTHRO_TG
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!  18 Dec 2009 - Aaron van D - Added EMISS_CAC_ANTHRO_05x0666 routine
!  18 Dec 2009 - Aaron van D - Added READ_CANADA_MASK_05x0666 routine
!  01 Mar 2012 - R. Yantosca - Remove A_CM2 array, use the new function
!                              GET_AREA_CM2( I, J, L ) from grid_mod.F90 
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!

      ! Arrays for data masks
      INTEGER, ALLOCATABLE :: MASK_CANADA_1x1(:,:)
      REAL*8,  ALLOCATABLE :: MASK_CANADA(:,:)
      
      ! Arrays for emissions
      REAL*8,  ALLOCATABLE :: NOx(:,:)
      REAL*8,  ALLOCATABLE :: CO(:,:)
      REAL*8,  ALLOCATABLE :: SO2(:,:)
      REAL*8,  ALLOCATABLE :: NH3(:,:)
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_canada_mask
!
! !DESCRIPTION: Function GET\_CANADA\_MASK returns the value of the Canadian 
!  geographic mask at grid box (I,J).  MASK=1 if (I,J) is within Canada, 
!  MASK=0 otherwise. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CANADA_MASK( I, J ) RESULT( THISMASK )
!
! !INPUT PARAMETERS: 
!
      ! Longitude and latitude indices
      INTEGER, INTENT(IN) :: I, J
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL*8              :: THISMASK

      !=================================================================
      ! GET_CANADA_MASK begins here!
      !=================================================================
      THISMASK = MASK_CANADA(I,J)

      END FUNCTION GET_CANADA_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_cac_anthro
!
! !DESCRIPTION: Function GET\_CAC\_ANTHRO returns the Critical Air Contaminants
!  emission for GEOS-Chem grid box (I,J) and tracer N.  Emissions can be 
!  returned in units of [kg/s] or [molec/cm2/s].  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_CAC_ANTHRO( I,    J,     N, 
     &                         MOLEC_CM2_S, KG_S ) RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTNO, IDTCO, IDTSO2, IDTNH3, IDTNO2
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
!  28 Jan 2009 - P. Le Sager - Initial Version
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                       :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_CAC_ANTHRO begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.
      
      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( N == IDTNO ) THEN

         ! NOx [kg/yr]
         VALUE = NOx(I,J)

      ELSE IF ( N == IDTCO ) THEN

         ! CO [kg/yr]
         VALUE = CO(I,J)

      ELSE IF ( N == IDTSO2 ) THEN

         ! SO2 [kg/yr]
         VALUE = SO2(I,J)

      ELSE IF ( N == IDTNH3 ) THEN

         ! NH3 [kg/month]
         VALUE = NH3(I,J)

      ELSE

         ! Otherwise return a negative value to indicate
         ! that there are no CAC emissions for tracer N
         VALUE = -1d0
         RETURN

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN
         IF ( N == IDTNH3 ) THEN
            ! Use 30 days per month (actual number of
            ! days may be required for the future)
            ! 2592000 = 30days*24hrs*60min*60sec
            VALUE = VALUE / 2592000
         ELSE
            
            ! Convert from [kg/yr] to [kg/s]
            VALUE = VALUE / SEC_IN_YEAR
         ENDIF

      ELSE IF ( DO_MCS ) THEN

         ! Convert from [kg/yr] to [molec/cm2/s]
         ! Updated on May 3, 2012 by Wai-Ho Lo: Not only NOx, but 
         ! also NH3.
         IF ( N == IDTNH3 ) THEN
            ! Use 30 days per month (actual number of
            ! days may be required for the future)
            ! 2592000 = 30days*24hrs*60min*60sec
            VALUE = VALUE * XNUMOL(N) / ( GET_AREA_CM2( I, J, 1 )
     &            * 2592000 )
         ELSEIF ( N == IDTNO ) THEN
            ! Convert NOx using XNUMOL for NO2 since original NOx
            ! emissions are in mass units of NO2 (mpayer, 4/18/13)
            VALUE = ( VALUE       * XNUMOL(IDTNO2)          )
     &            / ( SEC_IN_YEAR * GET_AREA_CM2( I, J, 1 ) )
         ELSE
            VALUE = ( VALUE       * XNUMOL(N)               )
     &            / ( SEC_IN_YEAR * GET_AREA_CM2( I, J, 1 ) )
         ENDIF
      ENDIF

      END FUNCTION GET_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_cac_anthro
!
! !DESCRIPTION: Subroutine EMISS\_CAC\_ANTHRO reads the Critical Air 
!  Contaminants emission fields at 1x1 resolution and regrids them to the 
!  current model resolution. (amv, phs, 1/28/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_CAC_ANTHRO( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
! 
      USE BPCH2_MOD,          ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_YEAR,  GET_MONTH
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR_1x1
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A
      USE CMN_SIZE_MOD
      USE CMN_O3_MOD

#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNO, IDTCO, IDTSO2, IDTNH3
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
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?!
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!
! !REMARKS:
!  (1 ) Emissions are read for a year b/w 2002-2005, and scaled
!        (except NH3) between 1985-2003 if needed (phs, 3/10/08)
!  (2 ) Now accounts for FSCALYR (phs, 3/17/08)
!  18 Dec 2009 - Aaron van D - Use 2005 scale factors for years beyond 2005
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bug in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  16 Jul 2013 - M. Payer    - Bug fix: Prevent negative emissions when scaling
!                              between 2002 and 2005 (C. Keller)
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
      REAL*8, TARGET             :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                     :: GEOS_1x1_2002(I1x1,J1x1,1)
      REAL*8                     :: GEOS_1x1_2005(I1x1,J1x1,1)
      REAL*8                     :: SC_1x1(I1x1,J1x1)
      REAL*8                     :: TAU2002, TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER(LEN=2)           :: THISMONTHCHAR
      REAL*8                     :: NH3_SCALE(12)
      CHARACTER(LEN=255)         :: LLFILENAME
      REAL*8, POINTER            :: INGRID(:,:) => NULL()

      ! For fields from Input_Opt
      LOGICAL                    :: LFUTURE
      INTEGER                    :: N_TRACERS

      !  seasonal scalar for NH3 emission (lzh, amv, 12/11/2009)
      ! Updated on May 13, 2012 by Wai-Ho Lo, since Agriculture Canada's
      ! NH3 emission inventory is used, monthly scalars are not used

      NH3_SCALE = (/
     &       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
      !&     0.426d0, 0.445d0, 0.526d0, 0.718d0, 1.179d0, 1.447d0,
      !&     1.897d0, 1.884d0, 1.577d0, 0.886d0, 0.571d0, 0.445d0 /)

      !=================================================================
      ! EMISS_CAC_ANTHRO begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      N_TRACERS = Input_Opt%N_TRACERS

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_CAC_ANTHRO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF
      THISMONTH = GET_MONTH()
      
      WRITE( THISMONTHCHAR, '(i2.2)' ) THISMONTH
      THISMONTHCHAR = ADJUSTL( THISMONTHCHAR )

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      DO SPECIES = 1,4

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
         ELSEIF ( SPECIES .eq. 4 ) THEN
            SNAME = 'NH3'
            SNo = 30
            ScNo = 0
         ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%% Prior to 6/28/12:
!%%% The CAC update goes only up to 2008.  Replace 2010 with 2008 in this
!%%% IF statement for the time being. (bmy, 6/28/12)
!%%%         IF ( ( THISYEAR .le. 2002 ) .OR. 
!%%%     &        ( THISYEAR .ge. 2010 ) ) THEN
!%%% 
!%%%            ! TAU values for 2002/2005
!%%%            TAU = GET_TAU0( 1, 1, MIN( MAX( THISYEAR, 2002 ), 2010 ) )
!%%%            WRITE( SYEAR, '(i4)' ) MIN( MAX( THISYEAR, 2002 ), 2010 )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

         IF ( ( ( THISYEAR .le. 2002 ) .OR.
     &        ( THISYEAR .ge. 2005 ) ) .and.
     &          ( THISYEAR .ne. 2006 ) ) THEN

            ! TAU values for 2002/2005
            TAU = GET_TAU0( 1, 1, MIN( MAX( THISYEAR, 2002 ), 2008 ) )
            WRITE( SYEAR, '(i4)' ) MIN( MAX( THISYEAR, 2002 ), 2008 )

            ! File name
            IF (SPECIES .eq. 4 ) THEN
               FILENAME = TRIM( DATA_DIR_1x1 )  //
     &                    'CAC_201307/CAC2008-' // 
     &                    TRIM( SNAME )         // '-' // 
     &                    TRIM( THISMONTHCHAR ) // '.geos.1x1'
            ELSE
               FILENAME  = TRIM( DATA_DIR_1x1 ) // 'CAC_200801/CAC' // 
     &                     SYEAR // '-' // TRIM( SNAME ) // '.geos.1x1'
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_CAC_ANTHRO: Reading ', a )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2005 data is read, a monthly
               ! TAU value has to be read for 2008 for NH3 emissions
               TAU = GET_TAU0( THISMONTH, 1, 2008 )
            ENDIF

            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo, 
     &                       TAU,      I1x1,       J1x1,
     &                       1,        ARRAY,      QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEOS_1x1(:,:,1) = ARRAY(:,:,1)

            ! Apply annual scalar factor. Available for 1985-2006,
            ! and NOx, CO and SO2 only.
            IF ( ( THISYEAR .lt. 2002 ) .and. SPECIES .ne. 4 ) THEN

               CALL GET_ANNUAL_SCALAR_1x1( ScNo,     2002,
     &                                     THISYEAR, SC_1x1 )

               GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

            ELSE IF ((THISYEAR .gt. 2005) .and. SPECIES .ne. 4) THEN

               CALL GET_ANNUAL_SCALAR_1x1( ScNo, MIN(THISYEAR, 2008),
     &                                     THISYEAR, SC_1x1 )

               GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * SC_1x1(:,:)

            ENDIF
            
         ELSE

            IF ( THISYEAR .ne. 2006 ) THEN
               TAU2002 = GET_TAU0( 1, 1, 2002)
               TAU2005 = GET_TAU0( 1, 1, 2005)
            ELSE
               TAU2002 = GET_TAU0( 1, 1, 2005)
               TAU2005 = GET_TAU0( 1, 1, 2007)
            ENDIF

            ! File name for 2002 data
            IF (SPECIES .eq. 4) THEN
               FILENAME = TRIM(DATA_DIR_1x1 )   //
     &                    'CAC_201307/CAC2008-' //
     &                    TRIM( SNAME )         // '-' //
     &                    TRIM( THISMONTHCHAR ) // '.geos.1x1'
            ELSE
               IF ( THISYEAR .ne. 2006 ) THEN
                  FILENAME  = TRIM( DATA_DIR_1x1 )       //
     &                        'CAC_200801/CAC2002-'      //
     &                        TRIM(SNAME) // '.geos.1x1'
               ELSE
                  FILENAME  = TRIM( DATA_DIR_1x1 )       //
     &                        'CAC_200801/CAC2005-'      //
     &                        TRIM(SNAME) // '.geos.1x1'
               ENDIF
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2002 or 2005 data is read, a
               ! monthly TAU value has to be read for 2008 for NH3
               ! emissions
               TAU = GET_TAU0( THISMONTH, 1, 2008 )
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU,       I1x1,      J1x1,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ELSE
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU2002,   I1x1,      J1x1,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ENDIF

            ! Cast to REAL*8 before regridding
            GEOS_1x1_2002(:,:,1) = ARRAY(:,:,1)

            ! File name for 2005 data
            IF (SPECIES .eq. 4) THEN
               FILENAME = TRIM(DATA_DIR_1x1 )   // 
     &                    'CAC_201307/CAC2008-' //
     &                    TRIM( SNAME )         // '-' //
     &                    TRIM( THISMONTHCHAR ) // '.geos.1x1'
            ELSE
               IF ( THISYEAR .ne. 2006 ) THEN
                  FILENAME  = TRIM( DATA_DIR_1x1 )       //
     &                        'CAC_200801/CAC2005-'      //
     &                        TRIM(SNAME) // '.geos.1x1'
               ELSE
                  FILENAME  = TRIM( DATA_DIR_1x1 )       //
     &                        'CAC_200801/CAC2007-'      //
     &                        TRIM(SNAME) // '.geos.1x1'
               ENDIF
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2002 or 2005 data is read, a
               ! monthly TAU value has to be read for 2008 for NH3
               ! emissions
               TAU = GET_TAU0( THISMONTH, 1, 2008 )
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU,       I1x1,      J1x1,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ELSE
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU2005,   I1x1,      J1x1,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ENDIF

            ! Cast to REAL*8 before regridding
            GEOS_1x1_2005(:,:,1) = ARRAY(:,:,1)

            ! Scale b/w years
            IF ( THISYEAR .ne. 2006 ) THEN
               GEOS_1x1(:,:,1) = GEOS_1x1_2002(:,:,1)
     &                   + ( THISYEAR - 2002.) / 3. *
     &                   ( GEOS_1x1_2005(:,:,1) - GEOS_1x1_2002(:,:,1) )
            ELSE
               GEOS_1x1(:,:,1) = GEOS_1x1_2002(:,:,1)
     &                   + ( THISYEAR - 2005.) / 2. *
     &                   ( GEOS_1x1_2005(:,:,1) - GEOS_1x1_2002(:,:,1) )
            ENDIF

            ! Bug fix: Prevent negative emissions (ckeller, mpayer, 7/16/13)
            WHERE ( GEOS_1x1 < 0d0 )
               GEOS_1x1 = 0d0
            END WHERE

         ENDIF

         ! Regrid from GEOS 1x1 --> current model resolution [kg/yr]
         INGRID => GEOS_1x1(:,:,1)

         IF ( SPECIES .eq. 1 ) THEN

            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1, 
     &                          INGRID,     NOx,     IS_MASS=1,
     &                          netCDF=.TRUE.                   )

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                          INGRID,     CO,      IS_MASS=1,
     &                          netCDF=.TRUE.                   )

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            !GEOS_1x1(:,:,1) = GEOS_1x1(:,:,1) * 0.986
            INGRID = INGRID * 0.986

            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                          INGRID,     SO2,     IS_MASS=1,
     &                          netCDF=.TRUE.                   )

         ELSEIF ( SPECIES .eq. 4 ) THEN

            ! Apply seasonality
            ! Using Agriculture Canada's NH3 emission inventory, 
            ! no seasonality scalars are required
            !GEOS_1x1(:,:,1) = NH3_SCALE(THISMONTH) * GEOS_1x1(:,:,1)

            CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1,
     &                          INGRID,     NH3,     IS_MASS=1,
     &                          netCDF=.TRUE.                   )

         ENDIF

         ! Free pointer
         NULLIFY( INGRID )

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN 
         CALL CAC_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

#if defined( DEVEL )
      State_Chm%TRAC_TEND(:,:,1,IDTNH3) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTNH3) + NH3

      State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTSO2) + SO2

      State_Chm%TRAC_TEND(:,:,1,IDTNO ) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTNO ) + NOx

      State_Chm%TRAC_TEND(:,:,1,IDTCO)  = 
     &State_Chm%TRAC_TEND(:,:,1,IDTCO)  + CO
#endif

      END SUBROUTINE EMISS_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_cac_anthro_05x0666
!
! !DESCRIPTION: Subroutine EMISS\_CAC\_ANTHRO\_05x0666 reads the Critical Air
!  Contaminants emission fields at nested NA resolution (1/2 x 2/3)
!  (amv, phs, 11/03/2009)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_CAC_ANTHRO_05x0666( am_I_Root, Input_Opt, 
     &                                     State_Chm, RC         )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TIME_MOD,           ONLY : GET_YEAR,  GET_MONTH
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR_05x0666_NESTED
      USE CMN_SIZE_MOD
      USE CMN_O3_MOD

#if defined( DEVEL )
      USE TRACERID_MOD,       ONLY : IDTNO, IDTCO, IDTSO2, IDTNH3
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
!  03 Nov 2009 - A. van Donkelaar - Initial Version
!
! !REMARKS:
!  (1 ) Emissions are read for a year b/w 2002-2005, and scaled
!        (except NH3) between 1985-2003 if needed (phs, 3/10/08)
!  (2 ) Now accounts for FSCALYR (phs, 3/17/08)
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  16 Jul 2013 - M. Payer    - Bug fix: Prevent negative emissions when scaling
!                              between 2002 and 2005 (C. Keller)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE              :: FIRST = .TRUE.
      INTEGER                    :: I, J, THISYEAR, SPECIES, SNo, ScNo
      INTEGER                    :: THISMONTH
      REAL*4                     :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                     :: GEOS_05x0666(IIPAR,JJPAR,1)
      REAL*8                     :: GEOS_05x0666_2002(IIPAR,JJPAR,1)
      REAL*8                     :: GEOS_05x0666_2005(IIPAR,JJPAR,1)
      REAL*4                     :: SC_05x0666(IIPAR,JJPAR)
      REAL*8                     :: TAU2002, TAU2005, TAU
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=4)           :: SYEAR, SNAME
      CHARACTER(LEN=2)           :: THISMONTHCHAR
      REAL*8                     :: NH3_SCALE(12)

      ! For fields from Input_Opt
      LOGICAL         :: LFUTURE
      INTEGER         :: N_TRACERS

      !  seasonal scalar for NH3 emission (lzh, amv, 12/11/2009)
      ! Updated on May 13, 2012 by Wai-Ho Lo, since Agriculture Canada's
      ! NH3 emission inventory is used, monthly scalars are not used
      NH3_SCALE = (/
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
      !&     0.426d0, 0.445d0, 0.526d0, 0.718d0, 1.179d0, 1.447d0,
      !&     1.897d0, 1.884d0, 1.577d0, 0.886d0, 0.571d0, 0.445d0 /)

      !=================================================================
      ! EMISS_CAC_ANTHRO begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      N_TRACERS = Input_Opt%N_TRACERS

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_CAC_ANTHRO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         THISYEAR = GET_YEAR()
      ELSE
         THISYEAR = FSCALYR
      ENDIF
      THISMONTH = GET_MONTH()

      WRITE( THISMONTHCHAR, '(i2.2)' ) THISMONTH
      THISMONTHCHAR = ADJUSTL( THISMONTHCHAR )


      DO SPECIES = 1,4

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
         ELSEIF ( SPECIES .eq. 4 ) THEN
            SNAME = 'NH3'
            SNo = 30
            ScNo = 0
         ENDIF

         IF ( ( ( THISYEAR .le. 2002 ) .OR.
     &        ( THISYEAR .ge. 2005 ) ) .and.
     &        ( THISYEAR .ne. 2006 ) ) THEN

            ! TAU values for 2002/2005
            TAU = GET_TAU0( 1, 1, MIN( MAX( THISYEAR, 2002 ), 2008 ) )
            WRITE( SYEAR, '(i4)' ) MIN( MAX( THISYEAR, 2002 ), 2008 )

            ! File name
            IF (SPECIES .eq. 4 ) THEN
               FILENAME = TRIM( DATA_DIR )      // 
     &                    'CAC_201307/CAC2008-' // 
     &                    TRIM( SNAME )         // '-' //
     &                    TRIM( THISMONTHCHAR ) // '.geos.1t2x2t3'
            ELSE
               FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC' //
     &                     SYEAR // '-' // TRIM( SNAME ) // 
     &                     '.geos.na.1t2x2t3'
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - EMISS_CAC_ANTHRO_05x0666: Reading ', a )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2002 or 2005 data is read, a 
               ! monthly TAU value has to be read for 2008 for NH3 
               ! emissions
               TAU = GET_TAU0(THISMONTH, 1, 2008 )
            ENDIF

            CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                       TAU,      IIPAR,       JJPAR,
     &                       1,        ARRAY,      QUIET=.TRUE. )

            ! Cast to REAL*8 before regridding
            GEOS_05x0666(:,:,1) = ARRAY(:,:,1)

            ! Apply annual scalar factor. Available for 1985-2006,
            ! and NOx, CO and SO2 only.
            IF ( ( THISYEAR .lt. 2002 ) .and. SPECIES .ne. 4 ) THEN

               CALL GET_ANNUAL_SCALAR_05x0666_NESTED( ScNo, 2002,
     &                                     THISYEAR, SC_05x0666 )

               GEOS_05x0666(:,:,1) = GEOS_05x0666(:,:,1)
     &                               * SC_05x0666(:,:)

            ELSE IF ((THISYEAR .gt. 2005) .and. SPECIES .ne. 4) THEN

               CALL GET_ANNUAL_SCALAR_05x0666_NESTED( ScNo,
     &                    MIN(THISYEAR,2008), THISYEAR, SC_05x0666 )

              GEOS_05x0666(:,:,1) = GEOS_05x0666(:,:,1)
     &                               * SC_05x0666(:,:)

            ENDIF
	

         ELSE

            IF ( THISYEAR .ne. 2006 ) THEN
               TAU2002 = GET_TAU0( 1, 1, 2002)
               TAU2005 = GET_TAU0( 1, 1, 2005)
            ELSE
               TAU2002 = GET_TAU0( 1, 1, 2005)
               TAU2005 = GET_TAU0( 1, 1, 2007)
            ENDIF

            ! File name for 2002 data
            IF (SPECIES .eq. 4 ) THEN
               FILENAME = TRIM( DATA_DIR )      //
     &                    'CAC_201307/CAC2008-' // 
     &                    TRIM( SNAME )         // '-' //
     &                    TRIM( THISMONTHCHAR ) // '.geos.1t2x2t3'
            ELSE
               IF ( THISYEAR .ne. 2006 ) THEN
                  FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2002-'
     &                     // TRIM(SNAME)      // '.geos.na.1t2x2t3'
               ELSE
                  FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2005-'
     &                     // TRIM(SNAME)      // '.geos.na.1t2x2t3'
               ENDIF
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2002 or 2005 data is read, a
               ! monthly TAU value has to be read for 2008 for NH3
               ! emissions
               TAU = GET_TAU0( THISMONTH, 1, 2008 )
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU,       IIPAR,     JJPAR,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ELSE
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU2002,   IIPAR,      JJPAR,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ENDIF

            ! Cast to REAL*8 before regridding
            GEOS_05x0666_2002(:,:,1) = ARRAY(:,:,1)

            ! File name for 2005 data
            IF (SPECIES .eq. 4 ) THEN
               FILENAME = TRIM( DATA_DIR )      //
     &                    'CAC_201307/CAC2008-' // 
     &                    TRIM( SNAME )         // '-' //
     &                    TRIM( THISMONTHCHAR ) // '.geos.1t2x2t3'
            ELSE
               IF ( THISYEAR .ne. 2006 ) THEN
                  FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2005-'
     &                      // TRIM(SNAME)     // '.geos.na.1t2x2t3'
               ELSE
                  FILENAME  = TRIM( DATA_DIR ) // 'CAC_200911/CAC2007-'
     &                      // TRIM(SNAME)     // '.geos.na.1t2x2t3'
               ENDIF
            ENDIF

            ! Echo info
            WRITE( 6, 100 ) TRIM( FILENAME )

            ! Read data
            IF (SPECIES .eq. 4 ) THEN
               ! Since currently the 2002 or 2005 data is read, a
               ! monthly TAU value has to be read for 2008 for NH3
               ! emissions
               TAU = GET_TAU0( THISMONTH, 1, 2008 )
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU,      IIPAR,     JJPAR,
     &                          1,        ARRAY,     QUIET=.TRUE. )
            ELSE
               CALL READ_BPCH2( FILENAME, 'ANTHSRCE', SNo,
     &                          TAU2005,   IIPAR,    JJPAR,
     &                          1,         ARRAY,     QUIET=.TRUE. )
            ENDIF

            ! Cast to REAL*8 before regridding
            GEOS_05x0666_2005(:,:,1) = ARRAY(:,:,1)

            ! Scale b/w 2002-2005
            IF ( THISYEAR .ne. 2006 ) THEN
               GEOS_05x0666(:,:,1) = GEOS_05x0666_2002(:,:,1) +
     &              ( THISYEAR - 2002.) / 3. *
     &              ( GEOS_05x0666_2005(:,:,1)
     &                - GEOS_05x0666_2002(:,:,1) )
            ELSE
               GEOS_05x0666(:,:,1) = GEOS_05x0666_2002(:,:,1) +
     &              ( THISYEAR - 2005.) / 2. *
     &              ( GEOS_05x0666_2005(:,:,1)
     &                - GEOS_05x0666_2002(:,:,1) )
            ENDIF

            ! Prevent negative emissions (ckeller, mpayer, 7/16/13)
            WHERE ( GEOS_05x0666 < 0d0 )
               GEOS_05x0666 = 0d0
            END WHERE

         ENDIF

         IF ( SPECIES .eq. 1 ) THEN

            NOx(:,:) =  GEOS_05x0666(:,:,1)

         ELSEIF ( SPECIES .eq. 2 ) THEN

            CO(:,:) = GEOS_05x0666(:,:,1)

         ELSEIF ( SPECIES .eq. 3 ) THEN

            ! Convert SOx to SO2, where SOx is assumed to be 1.4% SO4 and
            ! 98.6% SO2 over NA, based upon Chin et al, 2000, and as
            ! utilized in sulfate_mod.f
            SO2(:,:) = GEOS_05x0666(:,:,1) * 0.986

         ELSEIF ( SPECIES .eq. 4 ) THEN

            ! Apply seasonality
            !GEOS_05X0666(:,:,1) = NH3_SCALE(THISMONTH) 
            !    * GEOS_05X0666(:,:,1)
            NH3(:,:) = GEOS_05x0666(:,:,1)

         ENDIF

      ENDDO

      !--------------------------
      ! Compute future emissions
      !--------------------------
      IF ( LFUTURE ) THEN
         CALL CAC_SCALE_FUTURE
      ENDIF

      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_ANTHRO_Tg( THISYEAR )

#if defined( DEVEL )
      State_Chm%TRAC_TEND(:,:,1,IDTNH3) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTNH3) + NH3
      
      State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTSO2) + SO2

      State_Chm%TRAC_TEND(:,:,1,IDTNO)  = 
     &State_Chm%TRAC_TEND(:,:,1,IDTNO)  + NOx
      
      State_Chm%TRAC_TEND(:,:,1,IDTCO)  = 
     &State_Chm%TRAC_TEND(:,:,1,IDTCO)  + CO
#endif

      END SUBROUTINE EMISS_CAC_ANTHRO_05x0666
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cac_scale_future
!
! !DESCRIPTION: Subroutine CAC\_SCALE\_FUTURE applies the IPCC future scale 
!  factors to the Criteria Air Contaminant anthropogenic emissions.
!  (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:

      SUBROUTINE CAC_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NH3an 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

      USE CMN_SIZE_MOD             ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: I, J

      !=================================================================
      ! STREETS_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [kg NO2/yr]
         NOx(I,J)  = NOx(I,J) * GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO  [kg CO /yr]
         CO(I,J)   = CO(I,J)  * GET_FUTURE_SCALE_COff(  I, J )

         ! Future SO2 [kg SO2/yr] 
         SO2(I,J)  = SO2(I,J) * GET_FUTURE_SCALE_SO2ff( I, J )

         ! Future NH3 [kg NH3/yr] 
         NH3(I,J)  = NH3(I,J) * GET_FUTURE_SCALE_NH3an( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE CAC_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the 
!  anthropogenic emissions of NOx, CO, SO2 and NH3. (amv, phs, 1/28/09)
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
!  28 Jan 2009 - P. Le Sager - Initial Version
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,     J
      REAL*8              :: T_NOX, T_CO,  T_SO2,  T_NH3
      CHARACTER(LEN=3)    :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'C. A. C.   C A N A D I A N   E M I S S I O N S', / )


      ! Total NOx [Tg N]
      T_NOX = SUM( NOx ) * 1d-9 * ( 14d0 / 46d0 )

      ! Total CO  [Tg CO]
      T_CO  = SUM( CO  ) * 1d-9

      ! Total SO2 [Tg S]
      T_SO2 = SUM( SO2 ) * 1d-9 * ( 32d0 / 64d0 )

      ! Total NH3 [Tg NH3]
      T_NH3 = SUM( NH3 ) * 1d-9

      ! Print totals in [kg]
      WRITE( 6, 110 ) 'NO  ', YEAR, T_NOx,  '[Tg N  ]'
      WRITE( 6, 110 ) 'CO  ', YEAR, T_CO,   '[Tg CO ]'
      WRITE( 6, 110 ) 'SO2 ', YEAR, T_SO2,  '[Tg S  ]'
      WRITE( 6, 110 ) 'NH3 ', YEAR, T_NH3,  '[Tg NH3]'

      ! Format statement
 110  FORMAT( 'C.A.C. Canadian anthro ', a5, 
     &        'for year ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      END SUBROUTINE TOTAL_ANTHRO_Tg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_canada_mask
!
! !DESCRIPTION: Subroutine READ\_CANADA\_MASK reads and regrids the Canadian 
!  geographic mask from disk. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CANADA_MASK
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD         ! Size parameters
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  07 Jun 2012 - M. Payer    - Fixed minor bugs in map_a2a implementation
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
      REAL*8, TARGET          :: GEOS_1x1(I1x1,J1x1,1)
      REAL*8                  :: TAU2000
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8, POINTER         :: INGRID(:,:) => NULL()

      !=================================================================
      ! READ_CANADA_MASK begins here!
      !=================================================================

      TAU2000 = GET_TAU0(1,1,2000)

      ! File name
      FILENAME  = TRIM( DATA_DIR_1x1 ) // 
     &            'CAC_200801/CanadaMask.geos.1x1'

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos1x1.nc'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CANADA_MASK: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 TAU2000,   I1x1,     J1x1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEOS_1x1(:,:,1) = ARRAY(:,:,1)

      ! Save the 1x1 China mask for future use
      MASK_CANADA_1x1(:,:) = GEOS_1x1(:,:,1)

      ! Regrid from GENERIC 1x1 to current model resolution [unitless]
      INGRID => GEOS_1x1(:,:,1)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,        J1x1, 
     &                    INGRID,     MASK_CANADA, IS_MASS=0,
     &                    netCDF=.TRUE.                       )

      ! Free pointer
      NULLIFY( INGRID )

      WHERE ( MASK_CANADA > 0D0 ) MASK_CANADA = 1D0

      END SUBROUTINE READ_CANADA_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_canada_mask_05x0666
!
! !DESCRIPTION: Subroutine READ\_CANADA\_MASK\_05x0666 reads the Canadian
!  geographic mask from disk. (amv, phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_CANADA_MASK_05x0666
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR


      USE CMN_SIZE_MOD       ! Size parameters
!
! !REVISION HISTORY:
!  11 Nov 2009 - A. van Donkelaar - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                  :: GEOS_05x0666(IIPAR,JJPAR,1)
      REAL*8                  :: TAU2000
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! READ_CANADA_MASK begins here!
      !=================================================================

      TAU2000 = GET_TAU0(1,1,2000)

      ! File name
      FILENAME  = TRIM( DATA_DIR ) //
     &            'CAC_200911/CanadaMask.geos.na.1t2x2t3'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CANADA_MASK_05x0666: Reading ', a )

      ! Read data [unitless]
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2,
     &                 TAU2000,   IIPAR,     JJPAR,
     &                 1,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 before regridding
      MASK_CANADA(:,:) = ARRAY(:,:,1)

      END SUBROUTINE READ_CANADA_MASK_05x0666
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cac_anthro
!
! !DESCRIPTION: Subroutine INIT\_CAC\_ANTHRO allocates and zeroes all 
!  module arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CAC_ANTHRO( am_I_Root, Input_Opt, RC )
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
!  01 Mar 2012 - R. Yantosca - Delete the A_CM2 array, we will now just
!                              use the function directly
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS, J

      !=================================================================
      ! INIT_CAC_ANTHRO begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Return if LCAC is false
      IF ( .not. Input_Opt%LCAC ) RETURN
      
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

      ALLOCATE( NH3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
      NH3 = 0d0

      !---------------------------------------------------
      ! Read & Regrid masks for CAC emissions
      !---------------------------------------------------

      ALLOCATE( MASK_CANADA_1x1( I1x1, J1x1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CANADA_1x1' )
      MASK_CANADA_1x1 = 0

      ALLOCATE( MASK_CANADA( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASK_CANADA' )
      MASK_CANADA = 0d0

      ! Read China & SE Asia masks from disk
      CALL READ_CANADA_MASK

      END SUBROUTINE INIT_CAC_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cac_anthro
!
! !DESCRIPTION: Subroutine CLEANUP\_CAC\_ANTHRO deallocates all module 
!  arrays. (phs, 1/28/09)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CAC_ANTHRO
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CAC_ANTHRO begins here!
      !=================================================================
      IF ( ALLOCATED( MASK_CANADA_1x1) ) DEALLOCATE( MASK_CANADA_1x1) 
      IF ( ALLOCATED( MASK_CANADA    ) ) DEALLOCATE( MASK_CANADA    ) 
      IF ( ALLOCATED( NOx            ) ) DEALLOCATE( NOx            )
      IF ( ALLOCATED( CO             ) ) DEALLOCATE( CO             )
      IF ( ALLOCATED( SO2            ) ) DEALLOCATE( SO2            )
      IF ( ALLOCATED( NH3            ) ) DEALLOCATE( NH3            )

      END SUBROUTINE CLEANUP_CAC_ANTHRO
!EOC
      END MODULE CAC_ANTHRO_MOD

