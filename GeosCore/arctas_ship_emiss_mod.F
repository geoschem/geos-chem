!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: arctas_ship_emiss_mod
!
! !DESCRIPTION: Module ARCTAS\_SHIP\_EMISS\_MOD contains variables and 
!  routines to read the Arctas Ship emissions.
!\\
!\\
! !INTERFACE: 
!
      MODULE ARCTAS_SHIP_EMISS_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_ARCTAS_SHIP
      PUBLIC  :: EMISS_ARCTAS_SHIP
      PUBLIC  :: GET_ARCTAS_SHIP
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_ARCTAS_SHIP
      PRIVATE :: READ_ARCTAS_SHIP
      PRIVATE :: TOTAL_EMISS_TG
!
! !REVISION HISTORY:
!  28 Jan 2009 - P. Le Sager - Initial Version
!  31 Aug 2010 - R. Yantosca - Updated comments
!  01 Mar 2012 - R. Yantosca - Remove A_CM2 array, use the new function
!                              GET_AREA_CM2( I, J, L ) from grid_mod.F90 
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!
! !REMARKS:
!  (1) Part of the ARCTAS pre-campaign composite inventory distributed by
!       David Streets.
!  (2) Only SO2 differs from existing EDGAR inventory. All other
!       species are disregarded for now, except CO2 that we did not have
!       before.
!  (3) The ship emission is based on the work by Eyring et al., JGR 2005,
!      which estimates the total international ship emissions for 1985, 1990,
!      2001, and 2020 (projection). The ship emission for each individual year
!      is interpreted based on the above years, and the spatial pattern
!      (gridded) is mapped based on the EDGAR gridded ship emission for 2000
!      (total amount from EDGAR is scaled to Eyring-based number).
!
!      If you want to reference the work on publication or website, you may
!      either mention "Diehl et al., manuscript in preparation, 2009" or
!      refer to the AeroCom readme document for hindcast emissions for ship
!      (prepared by Diehl):
!      http://www-lscedods.cea.fr/aerocom/AEROCOM_HC/readme_ship/.
!
!      If you have further questions, please contact Thomas directly
!      (thomas.diehl@nasa.gov). 
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE DATA MEMBERS:
!
      ! Anthro emiss arrays
      REAL*8,  TARGET, ALLOCATABLE :: SO2_SHIP(:,:)
      REAL*8,  TARGET, ALLOCATABLE :: CO2_SHIP(:,:)
!
! !DEFINED PARAMETERS:
!
      REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

      CONTAINS

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_arctas_ship
!
! !DESCRIPTION: Function GET\_ARCTAS\_SHIP returns the ARCTAS\_SHIP emission 
!  for GEOS-Chem grid box (I,J) and tracer N.  Emissions can be returned in
!  units of [kg/s] or [molec/cm2/s].
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ARCTAS_SHIP( I, J, N, MOLEC_CM2_S, KG_S ) 
     &     RESULT( VALUE )
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,   ONLY : XNUMOL
      USE TRACERID_MOD, ONLY : IDTSO2
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL                       :: DO_KGS, DO_MCS

      !=================================================================
      ! GET_ARCTAS_SHIP begins here!
      !=================================================================

      ! Initialize
      DO_KGS = .FALSE.
      DO_MCS = .FALSE.

      ! Return data in [kg/s] or [molec/cm2/s]?
      IF ( PRESENT( KG_S        ) ) DO_KGS = KG_S
      IF ( PRESENT( MOLEC_CM2_S ) ) DO_MCS = MOLEC_CM2_S

      ! Test for simulation type
      IF ( ITS_A_CO2_SIM() ) THEN
         
         !-------------------
         ! CO2 simulation
         !-------------------
         VALUE  = CO2_SHIP(I,J)

      ELSE

         !-------------------
         ! Other simulations
         !-------------------
         IF ( N == IDTSO2 ) THEN

            ! SO2 [kg/yr]
            VALUE = SO2_SHIP(I,J)

         ELSE

            ! Otherwise return a negative value to indicate
            ! that there are no ARCTAS_SHIP emissions for tracer N
            VALUE = -1d0
            RETURN

         ENDIF

      ENDIF

      !------------------------------
      ! Convert units (if necessary)
      !------------------------------
      IF ( DO_KGS ) THEN

         ! Convert from [kg/yr] to [kg/s]         
         VALUE = VALUE / SEC_IN_YEAR              

      ELSE IF ( DO_MCS ) THEN

         ! Convert from [kg/yr] to [molec/cm2/s]            
         VALUE = ( VALUE       * XNUMOL(N)               )
     &         / ( SEC_IN_YEAR * GET_AREA_CM2( I, J, 1 ) )  

      ENDIF

      END FUNCTION GET_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_arctas_ship
!
! !DESCRIPTION: Subroutine EMISS\_ARCTAS\_SHIP reads the ARCTAS\_SHIP emissions
!  from disk.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISS_ARCTAS_SHIP( am_I_Root, YEAR, 
     &                              Input_Opt, State_Chm, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE DIRECTORY_MOD,      ONLY : DATA_DIR_1x1 
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTSO2
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: YEAR        ! Year of data to read
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
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: FIRST = .TRUE.
      CHARACTER(LEN=255) :: FILENAME, DIR

      ! For fields from Input_Opt
      INTEGER            :: N_TRACERS

      !=================================================================
      ! EMISS_ARCTAS_SHIP begins here!
      !=================================================================

      ! Copy values from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_ARCTAS_SHIP( am_I_Root, Input_Opt, RC )
         FIRST = .FALSE.
      ENDIF

 100  FORMAT( '     - EMISS_ARCTAS_SHIP: Reading ', a )

      ! Data directory
      DIR = TRIM( DATA_DIR_1x1 ) // 'ARCTAS_SHIP_2008/'
      

      IF ( Input_Opt%ITS_A_CO2_SIM ) THEN

         !--------------------------
         ! Read CO2 and regrid
         !--------------------------
         FILENAME  = TRIM( DIR ) // 'Arctas_CO2_ship_2008.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_ARCTAS_SHIP( FILENAME, 'CO2-SRCE', 1, CO2_SHIP,
     $        YEAR )
            
#if defined( DEVEL )
         State_Chm%TRAC_TEND(:,:,1,N_TRACERS+1) = 
     &   State_Chm%TRAC_TEND(:,:,1,N_TRACERS+1) + CO2_SHIP(:,:)
#endif
      ELSE

         !--------------------------
         ! Read SO2
         !--------------------------
         FILENAME  = TRIM( DIR ) // 'Arctas_SO2_ship_2008.generic.1x1'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! Read data
         CALL READ_ARCTAS_SHIP( FILENAME, 'ANTHSRCE', 26, SO2_SHIP,
     $        YEAR )
         
#if defined( DEVEL )
         State_Chm%TRAC_TEND(:,:,1,IDTSO2) = 
     &   State_Chm%TRAC_TEND(:,:,1,IDTSO2) + SO2_SHIP(:,:)
#endif

      ENDIF
     
      !--------------------------
      ! Print emission totals
      !--------------------------
      CALL TOTAL_EMISS_Tg( Input_Opt )

      END SUBROUTINE EMISS_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_arctas_ship
!
! !DESCRIPTION: Subroutine READ\_ARCTAS\_SHIP reads data from one ARCTAS\_SHIP 
!  data file from disk, at GENERIC 1x1 resolution and regrids them to the 
!  current model resolution.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_ARCTAS_SHIP( FILENAME, CATEGORY, 
     &                             TRACERN,  ARR,      YEAR )

!
! !USES:
! 
      USE BPCH2_MOD,        ONLY : GET_TAU0,      READ_BPCH2
      USE REGRID_A2A_MOD,   ONLY : DO_REGRID_A2A
      USE DIRECTORY_MOD,    ONLY : DATA_DIR_1x1      
      USE SCALE_ANTHRO_MOD, ONLY : GET_ANNUAL_SCALAR

      USE CMN_SIZE_MOD           ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! Year of data to read
      INTEGER,          INTENT(IN)    :: YEAR

      ! Tracer number
      INTEGER,          INTENT(IN)    :: TRACERN

      ! Filename & category under which data is stored in bpch file
      CHARACTER(LEN=*), INTENT(IN)    :: FILENAME, CATEGORY
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Array containing output data
      REAL*8,           INTENT(INOUT) :: ARR(IIPAR,JJPAR)
!
!
! !REMARKS:
!  (1) Even though the inventory was prepared for Arctas 2008 campaign, CO2 
!       base year is 2000, and SO2 base year is 2006. Input YEAR is used to 
!       scale SO2 into 1985-2005 
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!  13 Mar 2012 - M. Cooper   - Changed regrid algorithm to map_a2a
!  14 May 2012 - R. Yantosca - Bug fix: SC should be defined w/ IIPAR,JJPAR
!  24 May 2012 - R. Yantosca - Fix minor bugs in map_a2a implementation
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4              :: ARRAY(I1x1,J1x1-1,1)
      REAL*8, TARGET      :: GEN_1x1(I1x1,J1x1-1)
      REAL*4              :: SC(IIPAR,JJPAR)
      REAL*8              :: TAU2008
      CHARACTER(LEN=255)  :: LLFILENAME
      REAL*8, POINTER     :: INGRID (:,:) => NULL()
      REAL*8              :: OUTGRID(IIPAR,JJPAR)

      ! TAU0 values for 2008
      TAU2008 = GET_TAU0( 1, 1, 2008 )

      ! Initialize
      SC = 1d0

      ! Read data [kg/yr]
      CALL READ_BPCH2( FILENAME,  CATEGORY,  TRACERN, 
     &                 TAU2008,   I1x1,      J1x1-1,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8 before regridding
      GEN_1x1(:,:) = ARRAY(:,:,1)
  
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/MAP_A2A_latlon_generic1x1.nc'

      ! Regrid from GENERIC 1x1 --> Current Model Resolution
      INGRID => GEN_1x1(:,:)
      CALL DO_REGRID_A2A( LLFILENAME, I1x1,    J1x1-1, 
     &                    INGRID,     OUTGRID, IS_MASS=1, 
     &                    netCDF=.TRUE. )

      ! Get & apply scaling factor
      IF ( TRACERN == 26 )
     $     CALL GET_ANNUAL_SCALAR( 73, 2000, YEAR, SC )
      ARR = OUTGRID * SC

      ! Free pointer
      NULLIFY( INGRID )

      END SUBROUTINE READ_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TOTAL_EMISS_TG
!
! !DESCRIPTION: Subroutine TOTAL\_EMISS\_TG prints the totals for the 
!  anthropogenic or biomass emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TOTAL_EMISS_TG( Input_Opt )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8     :: T_SO2, T_CO2

      !=================================================================
      ! TOTAL_EMISS_TG begins here!
      !=================================================================

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100  )  
 100  FORMAT( 'A R C T A S   S H I P   E M I S S I O N S', / )


      ! Test for simulation type
      IF ( Input_Opt%ITS_A_CO2_SIM ) THEN

         !-----------------------
         ! CO2 simulation
         !-----------------------

         ! Total CO2 [Tg CO2]
         T_CO2 = SUM( CO2_SHIP ) * 1d-9

         ! Print totals
         WRITE( 6, 110 ) 'CO2 ', 2008, T_CO2,  ' CO2'

      ELSE

         !-----------------------
         ! Other simulations
         !-----------------------

         ! Total SO2 [Tg S]
         T_SO2 = SUM( SO2_SHIP ) * 1d-9 * ( 32d0 / 64d0 )

         ! Print totals in [Tg]
         WRITE( 6, 110 ) 'SO2 ', 2008, T_SO2, '[Tg S  ]'

      ENDIF

      ! Format statement
 110  FORMAT( 'ARCTAS SHIP ', a5, 
     &        'for base year ', i4, ': ', f11.4, 1x, a8 )

      ! Fancy output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      
      END SUBROUTINE TOTAL_EMISS_Tg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine INIT\_ARCTAS\_SHIP allocates and zeroes all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ARCTAS_SHIP( am_I_Root, Input_Opt, RC )
!
! !USES:
! 
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
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
!  28 Jan 2009 - P. Le Sager - Initial Version
!  01 Mar 2012 - R. Yantosca - Remove A_CM2 array
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Assume success
      RC =  GIGC_SUCCESS

      !=================================================================
      ! INIT_ARCTAS_SHIP begins here!
      !=================================================================

      ! Allocate ANTHRO arrays if LARCTAS_SHIP is TRUE
      IF ( Input_Opt%LARCSHIP ) THEN
      
         !--------------------------------------------------
         ! Allocate and zero arrays for SHIP emissions
         !--------------------------------------------------
         ! Test for simulation type
         IF ( Input_Opt%ITS_A_CO2_SIM ) THEN

            !-----------------------
            ! CO2 simulation
            !-----------------------
            ALLOCATE( CO2_SHIP( IIPAR, JJPAR ), STAT=RC )
            IF ( RC /= 0 ) CALL ALLOC_ERR( 'CO2_SHIP' )
            CO2_SHIP = 0d0

         ELSE

            !-----------------------
            ! Other simulations
            !-----------------------
            ALLOCATE( SO2_SHIP( IIPAR, JJPAR ), STAT=RC )
            IF ( RC /= 0 ) CALL ALLOC_ERR( 'SO2_SHIP' )
            SO2_SHIP = 0d0

         ENDIF

      ENDIF

      END SUBROUTINE INIT_ARCTAS_SHIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_ARCTAS_SHIP
!
! !DESCRIPTION: Subroutine CLEANUP\_ARCTAS\_SHIP deallocates all module 
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ARCTAS_SHIP
!
! !REVISION HISTORY: 
!  28 Jan 2009 - P. Le Sager - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_ARCTAS_SHIP begins here!
      !=================================================================
      IF ( ALLOCATED( SO2_SHIP     ) ) DEALLOCATE( SO2_SHIP     )
      IF ( ALLOCATED( CO2_SHIP     ) ) DEALLOCATE( CO2_SHIP     )

      END SUBROUTINE CLEANUP_ARCTAS_SHIP
!EOC
      END MODULE ARCTAS_SHIP_EMISS_MOD
