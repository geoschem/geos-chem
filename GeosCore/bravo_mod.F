!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bravo_mod
!
! !DESCRIPTION: \subsection*{Overview}
!  Module BRAVO\_MOD contains variables and routines to read the BRAVO 
!  Mexican anthropogenic emission inventory for NOx, CO, and SO2. 
!  (rjp, kfb, bmy, 6/22/06, 1/30/09)
!
! \subsection*{References}
! \begin{enumerate}
! \item Kuhns, H., M. Green, and Etyemezian, V, \emph{Big Bend Regional 
!       Aerosol and Visibility Observational (BRAVO) Study Emissions 
!       Inventory}, Desert Research Institute, 2003.
! \end{enumerate}
!
! !INTERFACE: 
!
      MODULE BRAVO_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_BRAVO
      PUBLIC  :: EMISS_BRAVO
      PUBLIC  :: GET_BRAVO_MASK
      PUBLIC  :: GET_BRAVO_ANTHRO
!
! !PRIVATE MEMBER FUNCTIONS:
!     
      PRIVATE :: BRAVO_SCALE_FUTURE
      PRIVATE :: INIT_BRAVO 
      PRIVATE :: READ_BRAVO_MASK
!
! !REVISION HISTORY:
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  (2 ) Now scale emissions using int-annual scale factors (amv, 08/24/07)
!  (3 ) Now accounts for FSCLYR (phs, 3/17/08)
!  (4 ) Added ProTeX headers (bmy, 1/30/09)
!  31 Aug 2010 - R. Yantosca - Updated comments
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
! 
      ! Arrays
      REAL*8,  ALLOCATABLE :: BRAVO_MASK(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_NOx(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_CO(:,:)
      REAL*8,  ALLOCATABLE :: BRAVO_SO2(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_bravo_mask
!
! !DESCRIPTION: Function GET\_BRAVO\_MASK returns the value of the Mexico 
!  mask for BRAVO emissions at grid box (I,J).  MASK=1 if (I,J) is in the 
!  BRAVO Mexican region, or MASK=0 otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BRAVO_MASK( I, J ) RESULT( MASK )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I        ! Longitude index
      INTEGER, INTENT(IN) :: J        ! Latitude  index
!
! !RETURN VALUE:
! 
      REAL*8              :: MASK     ! Returns the mask value @ (I,J)
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_BRAVO_MASK begins here!
      !=================================================================
      MASK = BRAVO_MASK(I,J)

      END FUNCTION GET_BRAVO_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_bravo_anthro
!
! !DESCRIPTION: Function GET\_BRAVO\_ANTHRO returns the BRAVO emission 
!  for GEOS-Chem grid box (I,J) and tracer N.  Units are [molec/cm2/s]. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BRAVO_ANTHRO( I, J, N ) RESULT( BRAVO )
!
! !USES:
! 
      USE TRACERID_MOD, ONLY : IDTNO, IDTCO, IDTSO2
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I       ! Longitude index
      INTEGER, INTENT(IN) :: J       ! Latitude index
      INTEGER, INTENT(IN) :: N       ! Tracer number
!
! RETURN VALUE:
! 
      REAL*8              :: BRAVO   ! Returns emissions at (I,J)
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!  (1 ) added SOx, SOx ship and NH3 emissions, plus optional kg/s output
!       (amv, 06/2008)
!  (2 ) Now returns ship emissions if requested (phs, 6/08)
!  (3 ) Added checks to avoid calling unavailable ship emissions (phs, 6/08)
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_BRAVO_ANTHRO begins here!
      !=================================================================

      ! NOx
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( N  == IDTNO ) THEN
         BRAVO = BRAVO_NOx(I,J)

      ! CO
      ELSE IF ( N == IDTCO ) THEN
         BRAVO = BRAVO_CO(I,J)

      ! SO2 
      ELSE IF ( N == IDTSO2 ) THEN
         BRAVO = BRAVO_SO2(I,J)

      ! Otherwise return a negative value to indicate
      ! that there are no BRAVO emissions for tracer N
      ELSE
         BRAVO = -1d0

      ENDIF

      END FUNCTION GET_BRAVO_ANTHRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_bravo
!
! !DESCRIPTION: Subroutine EMISS\_BRAVO reads the BRAVO emission fields at 1x1 
!  resolution and regrids them to the current model resolution. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_BRAVO( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
! 
      USE BPCH2_MOD,          ONLY : GET_TAU0,      READ_BPCH2
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE SCALE_ANTHRO_MOD,   ONLY : GET_ANNUAL_SCALAR
      USE TIME_MOD,           ONLY : GET_YEAR
      USE REGRID_A2A_MOD,     ONLY : DO_REGRID_A2A

      USE CMN_SIZE_MOD           ! Size parameters
      USE CMN_O3_MOD             ! 

#if defined( DEVEL )
      USE TRACERID_MOD, ONLY : IDTNO, IDTCO, IDTSO2
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
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!  (1 ) Now pass the unit string to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: SCALEYEAR
      REAL*4                  :: ARRAY(I1x1,J1x1-1,1)
      REAL*8, TARGET          :: GEN_1x1(I1x1,J1x1-1)
      REAL*4                  :: SC(IIPAR,JJPAR)
      REAL*8                  :: TAU0
      CHARACTER(LEN=255)      :: FILENAME
      CHARACTER(LEN=255)      :: LLFILENAME
      REAL*8, POINTER         :: INGRID(:,:) => NULL()
      REAL*8                  :: OUTGRID(IIPAR,JJPAR)

      ! For fields from Input_Opt
      LOGICAL                 :: LFUTURE
      INTEGER                 :: N_TRACERS

      !=================================================================
      ! EMISS_BRAVO begins here!
      !=================================================================

      ! Copy values from Input_Opt
      LFUTURE   = Input_Opt%LFUTURE
      N_TRACERS = Input_Opt%N_TRACERS

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_BRAVO( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read data from disk
      !=================================================================

      ! Use 1999 for BRAVO emission files (BASE YEAR)
      TAU0  = GET_TAU0( 1, 1, 1999 )
        
      ! Get emissions year
      IF ( FSCALYR < 0 ) THEN
         SCALEYEAR = GET_YEAR()
      ELSE
         SCALEYEAR = FSCALYR
      ENDIF

      !---------------------
      ! Read and regrid NOx
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.NOx.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - EMISS_BRAVO: Reading ', a )
      
      ! Read NOx [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 1, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid to NOx [molec/cm2/s] current model resolution
      INGRID => GEN_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )

      ! Get & apply scaling factor
      CALL GET_ANNUAL_SCALAR( 71, 1999, SCALEYEAR, SC )
      BRAVO_NOx = OUTGRID * SC

      ! Free pointer
      NULLIFY( INGRID )

      !---------------------
      ! Read and regrid CO
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.CO.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read CO [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 4, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid CO [molec/cm2/s] to current model resolution
      INGRID => GEN_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=0, 
     &                    netCDF=.TRUE.                   )


      ! Get & apply scaling factor
      CALL GET_ANNUAL_SCALAR( 72, 1999, SCALEYEAR, SC )
      BRAVO_CO = OUTGRID * SC

      ! Free pointer
      NULLIFY( INGRID )

      !---------------------
      ! Read and regrid SO2
      !---------------------

      ! 1x1 file name
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'BRAVO_200607/BRAVO.SO2.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
      
      ! Read SO2 [molec/cm2/s] on GENERIC 1x1 GRID
      CALL READ_BPCH2( FILENAME, 'ANTHSRCE', 26, 
     &                 TAU0,      I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1)

      ! Regrid SO2 [molec/cm2/s] to current model resolution
      INGRID => GEN_1x1
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=0,
     &                    netCDF=.TRUE.                   )

      ! Get & apply scaling factor
      CALL GET_ANNUAL_SCALAR( 73, 1999, SCALEYEAR, SC )
      BRAVO_SO2 = OUTGRID * SC

      ! Nullify pointer
      NULLIFY( INGRID )

      !=================================================================
      ! Compute IPCC future emissions (if necessary)
      !=================================================================
      IF ( LFUTURE ) THEN 
         CALL BRAVO_SCALE_FUTURE
      ENDIF

      !=================================================================
      ! Print emission totals
      !=================================================================
      CALL TOTAL_ANTHRO_TG( SCALEYEAR )

#if defined( DEVEL )
      State_Chm%TRAC_TEND(:,:,1,IDTNO)  =
     &State_Chm%TRAC_TEND(:,:,1,IDTNO)  + BRAVO_NOx

      State_Chm%TRAC_TEND(:,:,1,IDTCO)  = 
     &State_Chm%TRAC_TEND(:,:,1,IDTCO)  + BRAVO_CO

      State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &State_Chm%TRAC_TEND(:,:,1,IDTSO2) + BRAVO_SO2
#endif      

      END SUBROUTINE EMISS_BRAVO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: bravo_scale_future
!
! !DESCRIPTION: Subroutine BRAVO\_SCALE\_FUTURE applies the IPCC future 
!  scale factors to the BRAVO anthropogenic emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BRAVO_SCALE_FUTURE
!
! !USES:
! 
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_COff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_NOxff
      USE FUTURE_EMISSIONS_MOD, ONLY : GET_FUTURE_SCALE_SO2ff

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
      ! BRAVO_SCALE_FUTURE begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Future NOx [molec/cm2/s]
         BRAVO_NOx(I,J) = BRAVO_NOx(I,J)                * 
     &                    GET_FUTURE_SCALE_NOxff( I, J )

         ! Future CO [molec/cm2/s]
         BRAVO_CO(I,J)  = BRAVO_CO(I,J)                 *
     &                    GET_FUTURE_SCALE_COff( I, J )

         ! Future ALK4 [atoms C/cm2/s]
         BRAVO_SO2(I,J) = BRAVO_SO2(I,J)                *
     &                    GET_FUTURE_SCALE_SO2ff( I, J )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE BRAVO_SCALE_FUTURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: total_anthro_Tg
!
! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the amount of BRAVO 
!  anthropogenic emissions that are emitted each year.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_ANTHRO_TG( YEAR )
!
! !USES:
! 
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACERID_MOD, ONLY : IDTNO, IDTCO, IDTSO2

      USE CMN_SIZE_MOD     ! Size parameters
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)   :: YEAR
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!  (1 ) Now YEAR is input to reflect scaling factors applied (phs, 3/17/08)
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90 
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J
      REAL*8                :: A, B(3), NOx, CO, SO2
      CHARACTER(LEN=3)      :: UNIT

      !=================================================================
      ! TOTAL_ANTHRO_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )
 100  FORMAT( 'B R A V O   M E X I C A N   E M I S S I O N S', /,
     &        'Base year : 1999' )
      
      !----------------
      ! Sum emissions
      !----------------
      
      ! Define conversion factors for kg/molec
      ! (Undefined tracers will be zero)
      B(:) = 0d0
      IF ( IDTNO  > 0 ) B(1) = 1d0 / ( 6.0225d23 / 14d-3 )  ! Tg N
      IF ( IDTCO  > 0 ) B(2) = 1d0 / ( 6.0225d23 / 28d-3 )  ! Tg CO
      IF ( IDTSO2 > 0 ) B(3) = 1d0 / ( 6.0225d23 / 32d-3 )  ! Tg S

      ! Summing variables
      NOX = 0d0   
      CO  = 0d0 
      SO2 = 0d0 

      ! Loop over latitudes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Convert [molec/cm2/s] to [Tg]
         ! (Multiply by 1d-9 to convert from [kg] to [Tg])
         A = GET_AREA_CM2( I, J, 1 ) * 365.25d0 * 86400d0 * 1d-9 

         ! Sum emissions (list NOx as Tg N)
         NOX = NOX + ( BRAVO_NOX(I,J) * A * B(1) )
         CO  = CO  + ( BRAVO_CO (I,J) * A * B(2) )
         SO2 = SO2 + ( BRAVO_SO2(I,J) * A * B(3) )
      ENDDO
      ENDDO
 
      !----------------
      ! Print sums
      !----------------

      ! Print totals in [kg]
      WRITE( 6, 110   ) 'NO  ', YEAR, NOx, ' N'
      WRITE( 6, 110   ) 'CO  ', YEAR, CO,  '  '
      WRITE( 6, 110   ) 'SO2 ', YEAR, SO2, ' S'

 110  FORMAT( 'BRAVO anthropogenic ', a4, 
     &        'for year ', i4, ': ', f9.4, ' Tg', a2 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return to calling program
      END SUBROUTINE TOTAL_ANTHRO_TG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_bravo_mask
!
! !DESCRIPTION: Subroutine READ\_BRAVO\_MASK reads the Mexico mask from 
!  disk.  The Mexico mask is the fraction of the grid box (I,J) which lies 
!  w/in the BRAVO Mexican emissions region.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_BRAVO_MASK
!
! !USES:
! 
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD         ! Size parameters
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!  (1 ) Now pass UNIT to DO_REGRID_G2G_1x1 (bmy, 8/9/06)
!  06 Apr 2012 - M. Payer    - Changed regrid algorithm to map_a2a (M. Cooper)
!  24 May 2012 - R. Yantosca - Fixed minor bugs in map_a2a implementation
!  15 Aug 2012 - M. Payer    - Set mask to 1 if greater than 0 (L. Murray)
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4                  :: ARRAY  (I1x1,J1x1-1,1)
      REAL*8                  :: GEN_1x1(I1x1,J1x1-1  )
      REAL*8                  :: XTAU
      CHARACTER(LEN=255)      :: FILENAME 
      CHARACTER(LEN=255)      :: LLFILENAME

      !=================================================================
      ! READ_BRAVO_MASK begins here!
      !=================================================================

      ! File name
      FILENAME = TRIM( DATA_DIR_1x1 ) //
     &           'BRAVO_200607/BRAVO.MexicoMask.generic.1x1'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_MEXICO_MASK: Reading ', a )

      ! Get TAU0 for Jan 1985
      XTAU  = GET_TAU0( 1, 1, 1999 )

      ! Mask is stored in the bpch file as #2
      CALL READ_BPCH2( FILENAME, 'LANDMAP', 2, 
     &                 XTAU,      I1x1,     J1x1-1,     
     &                 1,         ARRAY,    QUIET=.TRUE. ) 

      ! Cast from REAL*4 to REAL*8
      GEN_1x1(:,:) = ARRAY(:,:,1) 

      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid to current model resolution [unitless]
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,       J1x1-1, 
     &                    GEN_1x1,    BRAVO_MASK, IS_MASS=0,
     &                    netCDF=.TRUE.                      )

      WHERE ( BRAVO_MASK > 0D0 ) BRAVO_MASK = 1D0

      END SUBROUTINE READ_BRAVO_MASK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_bravo
!
! !DESCRIPTION: Subroutine INIT\_BRAVO allocates and zeroes BRAVO module 
!  arrays, and also creates the mask which defines the Mexico region 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_BRAVO( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS

      !=================================================================
      ! INIT_BRAVO begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Return if LBRAVO is false
      IF ( .not. Input_Opt%LBRAVO ) RETURN
      
      !--------------------------
      ! Allocate and zero arrays
      !--------------------------

      ALLOCATE( BRAVO_NOx( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_NOx' )
      BRAVO_NOx = 0d0

      ALLOCATE( BRAVO_CO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_CO' )
      BRAVO_CO = 0d0

      ALLOCATE( BRAVO_SO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_SO2' )
      BRAVO_SO2 = 0d0

      !--------------------------
      ! Read Mexico mask
      !--------------------------
     
      ALLOCATE( BRAVO_MASK( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BRAVO_MASK' )
      BRAVO_MASK = 0d0
      
      ! Read the mask
      CALL READ_BRAVO_MASK

      END SUBROUTINE INIT_BRAVO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_BRAVO
!
! !DESCRIPTION: Subroutine CLEANUP\_BRAVO deallocates all BRAVO module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_BRAVO
!
! !REVISION HISTORY:
!  22 Jun 2006 - R. Park, F. Boersma, R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_BRAVO begins here!
      !=================================================================
      IF ( ALLOCATED( BRAVO_NOx  ) ) DEALLOCATE( BRAVO_NOx  )
      IF ( ALLOCATED( BRAVO_CO   ) ) DEALLOCATE( BRAVO_CO   )
      IF ( ALLOCATED( BRAVO_SO2  ) ) DEALLOCATE( BRAVO_SO2  )
      IF ( ALLOCATED( BRAVO_MASK ) ) DEALLOCATE( BRAVO_MASK )

      END SUBROUTINE CLEANUP_BRAVO
!EOC
      END MODULE BRAVO_MOD
