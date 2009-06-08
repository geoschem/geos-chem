! $Id: co2_mod.f,v 1.5 2009/06/08 14:09:32 ccarouge Exp $
      MODULE CO2_MOD
!
!******************************************************************************
!  Module CO2_MOD contains variables and routines used for the CO2 simulation.
!  (pns, bmy, 8/16/05, 9/27/06) 
!
!  Module Variables:
!  ============================================================================
!  (1 ) EMFOSSCO2    : Array for fossil fuel CO2 emissions (annual mean)
!  (2 ) EMBIOCO2     : Balanced biosphere CO2 (CASA) emissions (daily)   
!  (3 ) EMOCCO2      : Ocean CO2 emissions (annual mean)  
!  (4 ) EMBIOBRNCO2  : Biomass burning emissions  
!  (5 ) EMBIOFUELCO2 : Biofuel emissions 
!  (6 ) EMBIONETCO2  : Net terrestrial CO2 emissions
!  (7 ) FMOL_CO2     : Molecular weight of CO2 
!  (8 ) LFOSSCO2     : Flag for switching on/off FOSSIL FUEL CO2 emissions
!  (9 ) LBIOCO2      : Flag for switching on/off BALANCED BIOSPHERE CO2 emiss. 
!  (10) LOCCO2       : Flag for switching on/off OCEAN CO2 emissions     
!  (11) LBIOBRNCO2   : Flag for switching on/off BIOMASS BURNING CO2 emissions 
!  (12) LBIOFUELCO2  : Flag for switching on/off BIOFUEL CO2 emissions 
!  (13) LBIONETCO2   : Flag for switching on/off NET BIOSPHERE EXCHANGE of CO2
!  (14) LUSECASANEP  : Flag for reading daily CASA NEP w/ diurnal cycle
!  (15) XNUMOL_CO2   : molec CO2 / kg CO2 
!
!  Module Procedures:
!  ============================================================================
!  (1 ) EMISSCO2               : Emits CO2 into individual tracers
!  (2 ) READ_ANNUAL_FOSSILCO2  : Reads annual mean emission fields for CO2
!  (3 ) READ_ANNUAL_OCEANCO2   : Reads annual mean CO2 ocean emissions
!  (4 ) READ_BBIO_DAILYAVERAGE : Reads daily mean CASA Bal Bio CO2 (no diurnal)
!  (5 ) READ_BBIO_DIURNALCYCLE : Reads CASA NEP fluxes w/ imposed diurnal cycle
!  (6 ) READ_MONTH_BIOBRN_CO2  : Read monthly biomass burning emissions
!  (7 ) READ_ANNUAL_BIOFUELCO2 : Read annual mean biofuel emissions
!  (8 ) READ_ANNUAL_BIONET_CO2 : Read annual net terrestrial exchange
!  (9 ) INIT_CO2               : Allocates and initializes module arrays
!  (10) CLEANUP_CO2            : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by "co2_mod.f"
!  ============================================================================
!  (1 ) biomass_mod.f          : Module w/ routines for biomass burning
!  (2 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (3 ) diag04_mod.f           : Module w/ routines for CO2 diagnostics
!  (4 ) directory_mod.f        : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f            : Module w/ I/O error and NaN check routines
!  (6 ) file_mod.f             : Module w/ file unit numbers and error checks
!  (7 ) grid_mod.f             : Module w/ horizontal grid information
!  (8 ) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (9 ) time_mod.f             : Module w/ routines for computing time & date
!  (10) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.
!  (11) transfer_mod.f         : Module w/ routines to cast & resize arrays 
!
!  CO2 tracers:
!  ============================================================================
!  (1 ) Total CO2
!  (2 ) CO2 from oceans
!  (3 ) CO2 from fossil fuel 
!  (4 ) CO2 from balanced biosphere
!  (5 ) CO2 from biomass burning
!  (6 ) CO2 from biofuels
!
!  References:
!  ============================================================================
!  (1 ) Andres, R.J, G. Marland, I. Fung, and E. Matthews, "A 1x1 distribution
!        of carbon dioxide emissions from fossil fuel consumption and
!        cement manufacture", Glob. Biogeochem. Cycles, Vol 10, 419-429, 1996.
!  (2 ) Randerson, J.T, M.V. Thompson, T.J.Conway, I.Y. Fung, and C.B. Field,
!        "The contribution of terrestrial sources and sinks to trends in the
!        seasonal cycle of atmospheric carbon dioxide", Glob. Biogeochem. 
!        Cycles, Vol 11, 535-560, 1997.
!  (3 ) Takahashi, T, R. Feely, R. Weiss, R. Wanninkof, D. Chipman, 
!        S. Sutherland, and T. Takahashi, "Global air-sea flux of CO2: An
!        estimate based on measurements of sea-air pCO2 difference", 
!        Proceedings of the National Academy of Sciences, 94, 8292-8299,
!        1997
!  (4 ) Yevich, R. and J. A. Logan, "An assesment of biofuel use and burning 
!        of agricultural waste in the developing world", Glob. Biogeochem. 
!        Cycles, Vol 17, 1095, doi:10.1029/2002GB001952, 2003
!
!  NOTES: 
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Now references biomass_mod.f (bmy, 9/27/06)
!******************************************************************************
!
      IMPLICIT NONE 

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "co2_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 

      ! ... except these routines
      PUBLIC :: CLEANUP_CO2
      PUBLIC :: EMISSCO2

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Logical switches
      LOGICAL, PARAMETER   :: LFOSSCO2     = .TRUE.
      LOGICAL, PARAMETER   :: LBIOCO2      = .TRUE.
      LOGICAL, PARAMETER   :: LOCCO2       = .TRUE.
      LOGICAL, PARAMETER   :: LBIOBRNCO2   = .TRUE.
      LOGICAL, PARAMETER   :: LBIOFUELCO2  = .TRUE.
      LOGICAL, PARAMETER   :: LBIONETCO2   = .FALSE. 
      LOGICAL, PARAMETER   :: LUSECASANEP  = .FALSE.

      ! Arrays
      REAL*8,  ALLOCATABLE :: EMFOSSCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMOCCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOBRNCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIOFUELCO2(:,:)
      REAL*8,  ALLOCATABLE :: EMBIONETCO2(:,:)

      ! FMOL_CO2     - kg CO2 / mole CO2 
      REAL*8,  PARAMETER   :: FMOL_CO2   = 44d-3

      ! XNUMOL_CO2   - molecules CO2 / kg CO2 
      REAL*8,  PARAMETER   :: XNUMOL_CO2 = 6.022d+23 / FMOL_CO2

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!-----------------------------------------------------------------------------
   
      SUBROUTINE EMISSCO2
!
!******************************************************************************
!  Subroutine EMISSCO2 is the driver routine for CO2 emissions. 
!  (pns, bmy, 8/16/05, 9/27/06)
!
!  The initial condition for CO2 has to be at least 50 ppm or higher or else
!  the balanced biosphere fluxes will make STT negative. (pns, bmy, 8/16/05)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) We now get CO2 biomass emissions from biomass_mod.f.  This allows us 
!        to use either GFED2 or default Duncan et al biomass emissions. 
!        (bmy, 9/27/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,   ONLY : BIOMASS,       IDBCO2
      USE DIAG04_MOD,    ONLY : AD04,          ND04
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_DAY,       GET_DAY_OF_YEAR
      USE TIME_MOD,      ONLY : GET_HOUR,      GET_MONTH
      USE TIME_MOD,      ONLY : GET_YEAR,      GET_TS_CHEM 
      USE TIME_MOD,      ONLY : ITS_A_NEW_DAY, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS,     STT
      
#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      LOGICAL, SAVE          :: FIRST  = .TRUE.

      ! Local variables
      INTEGER                :: I,     IJLOOP, J,    L,     N
      INTEGER                :: DAY,   DOY,    HOUR, MONTH, YEAR   
      REAL*8                 :: A_CM2, DTSRCE, E_CO2

      !=================================================================
      ! EMISSCO2 begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN 

         ! Allocate arrays and read annual-mean data
         CALL INIT_CO2

         ! Set first-time flag to false
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read in monthly and daily emissions fields
      !=================================================================      

      ! Emission timestep
      DTSRCE = 60d0 * GET_TS_CHEM()
      
      ! Time variables
      DAY    = GET_DAY()
      DOY    = GET_DAY_OF_YEAR()
      HOUR   = GET_HOUR()
      MONTH  = GET_MONTH()
      YEAR   = GET_YEAR()

      ! Check if Balanced Biosphere emissions are required  
      IF ( LBIOCO2 ) THEN  

         ! If LUSECASANEP is TRUE ...
         IF ( LUSECASANEP ) THEN

            ! ... then use 3-hourly NEP emissions for Bal Bio ...
            IF ( MOD( HOUR, 3 ) == 0 ) THEN
               CALL READ_BBIO_DIURNALCYCLE( MONTH, DAY, HOUR )
            ENDIF

         ELSE

            ! ... otherwise use constant daily emissions of NEP for Bal Bio
            IF ( ITS_A_NEW_DAY() ) THEN
               CALL READ_BBIO_DAILYAVERAGE( MONTH, DAY, DOY, YEAR ) 
            ENDIF

         ENDIF

      ENDIF

      !=================================================================
      ! Process emissions and save diagnostics
      !=================================================================

      ! Loop over latitudes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, A_CM2, E_CO2 )
      DO J = 1, JJPAR

         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

      ! Loop over longitudes
      DO I = 1, IIPAR

         !-------------------------------------------
         ! #1: Total CO2
         ! #2: CO2 from fossil fuel emissions 
         !-------------------------------------------
         IF ( LFOSSCO2 ) THEN

            ! Fossil fuel emissions of CO2 [molec/cm2/s]
            E_CO2          = EMFOSSCO2(I,J)
            
            ! ND04 diag: Fossil Fuel CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,1) = AD04(I,J,1) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #2: Fossil CO2 [kg]
            STT(I,J,1,2)   = STT(I,J,1,2) + E_CO2
         ENDIF

         !-------------------------------------------
         ! #3: CO2 from ocean emissions
         !-------------------------------------------
         IF ( LOCCO2 ) THEN

            ! Ocean CO2 emissions in [molec/cm2/s]
            E_CO2          = EMOCCO2(I,J)

            ! ND04 diag: Ocean CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,2) = AD04(I,J,2) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #3: Ocean CO2 [kg]
            STT(I,J,1,3)   = STT(I,J,1,3) + E_CO2   
         ENDIF

         !-------------------------------------------
         ! #4: CO2 from balanced biosphere emissions
         !-------------------------------------------
         IF ( LBIOCO2 ) THEN

            ! Balanced biosphere CO2 [molec/cm2/s]
            E_CO2         = EMBIOCO2(I,J)

            ! ND04 diag: Bal Bio CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,3) = AD04(I,J,3) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 

            ! Add BalBio CO2 to Tracer #1 -- total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add BB CO2 to Tracer #4 -- Bal Bio CO2 [kg CO2]
            STT(I,J,1,4)   = STT(I,J,1,4) + E_CO2
         ENDIF

         !-------------------------------------------
         ! #5: CO2 from biomass burning emissions
         !-------------------------------------------
         IF ( LBIOBRNCO2 ) THEN

            ! Biomass burning emissions [molec/cm2/s]
            E_CO2          = BIOMASS(I,J,IDBCO2)

            ! ND04 diag: Biomass burning CO2 [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,4) = AD04(I,J,4) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2 

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            
            ! Add to Tracer #5: Biomass burning CO2 [kg CO2]
            STT(I,J,1,5)   = STT(I,J,1,5) + E_CO2
         ENDIF

         !-------------------------------------------
         ! #6: CO2 from biofuel emissions
         !-------------------------------------------
         IF ( LBIOFUELCO2 ) THEN

            ! Biofuel CO2 emissions [molec/cm2/s]
            E_CO2          = EMBIOFUELCO2(I,J)

            ! ND04 diag: Biofuel CO2 [molec/cm2/s] 
            IF ( ND04 > 0 ) THEN
               AD04(I,J,5) = AD04(I,J,5) + E_CO2
            ENDIF

            ! Convert E_CO2 from [molec CO2/cm2/s] to [kg CO2]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! Add to Tracer #1: Total CO2 [kg CO2]
            STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2

            ! Add to Tracer #6: Biofuel CO2 [kg CO2]
            STT(I,J,1,6)   = STT(I,J,1,6) + E_CO2
         ENDIF

         !-------------------------------------------
         ! #7: CO2 from net terrestrial exchange
         !-------------------------------------------
         IF ( LBIONETCO2 ) THEN

            ! CO2 from net terrestrial exchange [molec/cm2/s]
            E_CO2          = EMBIONETCO2(I,J)

            ! ND04 diag: net terrestrial exchange [molec/cm2/s]
            IF ( ND04 > 0 ) THEN
               AD04(I,J,6) = AD04(I,J,6) + E_CO2
            ENDIF

            ! Convert from [molec/cm2/s] to [kg]
            E_CO2          = E_CO2 * A_CM2 * DTSRCE / XNUMOL_CO2

            ! DO NOT ADD NET TERR EXCHANGE to Tracer #1: total CO2 [kg]
            ! REMOVE Comment out below if need to add to Tracer #1
            !STT(I,J,1,1)   = STT(I,J,1,1) + E_CO2
            
            ! Add to Tracer #7: Net Terr exchange CO2 [kg]
            STT(I,J,1,7)   = STT(I,J,1,7) + E_CO2
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMISSCO2 

!-----------------------------------------------------------------------------

        SUBROUTINE READ_ANNUAL_FOSSILCO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_FOSSILCO2 reads in annual mean fossil CO2 emissions 
!  from a binary punch file. (pns, bmy, 8/16/05, 10/3/05)
!
!  References:
!  ============================================================================
!  (1 ) CDIAC gridded (1x1) dataset for 1995 (Andres et al.)
!
!  NOTES:
!  (1 ) Emissions read in from directory : DATA_DIR/CO2XXX
!  (2 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_FOSSILCO2 begins here!
      !=================================================================

      ! File contaning fossil fuel CO2 data
      FILENAME = TRIM( DATA_DIR )           // 
     &           'CO2_200508/fossil95_CO2.' // GET_NAME_EXT_2D() //
     &           '.'                        // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_FOSSCO2: Reading ', a ) 
      
      ! TAU value for start of "generic" year 1985 
      TAU = GET_TAU0( 1, 1, 1985 )

      ! Read fossil fuel CO2 [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 1, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), EMFOSSCO2 )

      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_FOSSILCO2

!-----------------------------------------------------------------------------

       SUBROUTINE READ_ANNUAL_OCEANCO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_OCEANCO2 reads in annual mean oceanic CO2 exchange  
!  from a binary punch file. (pns, bmy, 8/16/05, 10/3/05)
!
!  References:
!  ============================================================================
!  (1 ) Takahashi et al. (1997)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_OCEANCO2 begins here!
      !=================================================================

      ! File contaning ocean CO2 data
      FILENAME = TRIM( DATA_DIR )        // 
     &           'CO2_200508/ocean_CO2.' // GET_NAME_EXT_2D() //
     &           '.'                     // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_OCEANCO2: Reading ', a ) 

      ! TAU value for start of "generic" year 1985 
      TAU = GET_TAU0( 1, 1, 1985 )

      ! Read ocean CO2 data [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 2, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMOCCO2 )

      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_OCEANCO2

!------------------------------------------------------------------------------

        SUBROUTINE READ_ANNUAL_BIOFUELCO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_BIOFUELCO2 reads in annual mean biofuel CO2 
!  emissions from a binary punch file (pns, bmy, 8/16/05, 10/3/05)
!
!  References:
!  ============================================================================
!  (1 ) Yevich and Logan 2001 gridded (1x1) dataset in combination with 
!        emission factors for CO2 per kg drymatter burned
!  (2 ) See routines in /users/trop/pns/GEOSCHEM/EMISSIONS/BIOFUEL
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_BIOFUELCO2 begins here!
      !=================================================================

      ! File contaning biofuel CO2 data
      FILENAME = TRIM( DATA_DIR )          // 
     &           'CO2_200508/biofuel_CO2.' // GET_NAME_EXT_2D() //
     &           '.'                       // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_BIOFUELCO2: Reading ', a ) 

      ! Read biofuel CO2 emissions [molec/cm2/s]
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 5, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Cast from REAL*4 to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOFUELCO2 )

      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_BIOFUELCO2

!------------------------------------------------------------------------------
      
      SUBROUTINE READ_ANNUAL_BIONET_CO2
!
!******************************************************************************
!  Subroutine READ_ANNUAL_BIONET_CO2 reads in annual mean values of for Net 
!  Terrestrial exchange from a binary punch file. (pns, bmy, 8/16/05)
!
!  References:
!  ============================================================================
!  (1 ) Currently only exchange for Year 2000
!        Source : Transcom 3 Interannual inversion results
!        From David Baker (pers. comm.)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR 
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER                :: I, J, IOS
      REAL*4                 :: ARRAY(IIPAR,JJPAR,1)
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_ANNUAL_BIONET_CO2 begins here!
      !=================================================================
      
      ! Filename
      FILENAME = TRIM( DATA_DIR )                //
     &           'CO2_200508/net_terr_exch_CO2.' // GET_NAME_EXT_2D() //
     &           '.'                             // GET_RES_EXT()     //
     &           '.txt'

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_ANNUAL_BIONETCO2: Reading ', a ) 

      ! Initialize ARRAY
      ARRAY = 0e0

      ! Open file
      OPEN( IU_FILE,          FILE=TRIM( FILENAME ), 
     &      FORM='FORMATTED', IOSTAT=IOS )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_ann_bionet:1' )      

      ! Read data
      READ( IU_FILE, '(7e13.6)', IOSTAT=IOS )
     &  ( ( ARRAY(I,J,1), I=1,IIPAR), J=1,JJPAR )
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'read_ann_bionet:2' )

      ! Close file
      CLOSE( IU_FILE )

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIONETCO2 )

      ! Return to calling program
      END SUBROUTINE READ_ANNUAL_BIONET_CO2 

!------------------------------------------------------------------------------

      SUBROUTINE READ_BBIO_DAILYAVERAGE( MONTH, DAY, DOY, YEAR ) 
!
!******************************************************************************
!  Subroutine READ_DAILY_BBIO_CO2 reads in daily values for balanced 
!  biospheric exchange from a binary punch file.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month (1-31)
!  (3 ) DOY   (INTEGER) : Current day of year (0-366)
!
!  Data Sources:
!  ============================================================================
!  (1 ) CASA gridded (1x1) dataset for from M. Thompson
!        Monthly values interpolated to daily values : 365 daily files 
!        NB : These files DO NOT have the diurnal cycle in daily emissions
!        See routine ' ' to read in files with diurnal cycle imposed
!
!  References
!  ============================================================================
!  (1 ) Randerson et al. [1997]
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Modify DOY for non leap-year to correspond to DOY for leap-year (2000).
!       (ccc, 6/4/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE TIME_MOD,      ONLY : ITS_A_LEAPYEAR

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: MONTH, DAY, DOY, YEAR

      ! Local variables
      INTEGER                :: I, J
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3  )     :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_BBIO_DAILYAVERAGE begins here!
      !=================================================================

      ! Warning: DOY may be different depending on leap-year or not.
      ! (ccc, 6/4/09)
      IF ( .NOT.ITS_A_LEAP_YEAR( YEAR ) .AND. DOY > 59 ) THEN
         DOY = DOY + 1
      ENDIF

      ! Make a string from DOY
      WRITE( SDOY, '(i3.3)' ) DOY 

      ! Name of file with Balanced Bio CO2 data
      FILENAME = TRIM( DATA_DIR )                      //
     &           'CO2_200508/BBIO_DAILYAVG/CO2.daily.' //
     &           GET_NAME_EXT_2D() // '.'              //
     &           GET_RES_EXT()     // '.'              // SDOY

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DAILYAVERAGE: Reading ', a ) 

      ! Get TAU value corresponding to DOY in year 2000
      TAU = GET_TAU0( MONTH, DAY, 2000 )

      ! Read balanced biosphere CO2 [molec/cm2/s] from disk
      CALL READ_BPCH2( FILENAME, 'CO2-SRCE', 3, 
     &                 TAU,       IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )


      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )

      ! Return to calling program
      END SUBROUTINE READ_BBIO_DAILYAVERAGE

!------------------------------------------------------------------------------

      SUBROUTINE READ_BBIO_DIURNALCYCLE( MONTH, DAY, HOUR )
!
!******************************************************************************
!  Subroutine READ_BBIO_DIURNALCYCLE reads CASA daily Net Ecosystem Production
!  (NEP) fluxes but with a diurnal cycle imposed.  (pns, bmy, 8/16/05, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) MONTH (INTEGER) : Current month of year (1-12)
!  (2 ) DAY   (INTEGER) : Current day of month  (1-31)
!  (3 ) HOUR  (INTEGER) : Current hour of day   (0-23)
!
!  References:
!  ============================================================================
!  (1 ) Olsen et al. [2004]. Fluxes are Net Ecosystem Production (NEP) 
!        from CASA model
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D 

      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: MONTH, DAY, HOUR 

      ! Local variables
      INTEGER                :: DOY
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: TAU
      CHARACTER(LEN=3 )      :: SDOY
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_BBIO_DIURNALCYCLE begins here!
      !=================================================================
      
      ! Create string for day of year
      WRITE( SDOY, '(i3.3)' ) DOY

      ! File name
      FILENAME = TRIM( DATA_DIR )                    //
     &           'CO2_200508/BBIO_DIURNALCYCLE/nep.' // 
     &           GET_NAME_EXT_2D() // '.'            //
     &           GET_RES_EXT()     // '.'            // SDOY

      ! Echo file name to stdout
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_BBIO_DIURNALCYCLE: Reading ', a )

      ! Get TAU of this month & day in "generic" year 1985
      TAU = GET_TAU0( MONTH, DAY, 1985, HOUR ) 
         
      ! Read Net Ecosytem Productivity [molec CO/cm2/s] from disk
      ! The CASA fluxes use atmospheric convention: 
      ! positive = into atm; negative = into biosphere
      CALL READ_BPCH2( FILENAME, 'GLOB-NPP', 2, 
     &                 TAU,       IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )
       
      ! Cast from REAL*4 to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), EMBIOCO2 )
         
      ! Return to calling program
      END SUBROUTINE READ_BBIO_DIURNALCYCLE

!------------------------------------------------------------------------------

      SUBROUTINE TOTAL_BIOMASS_TG( BBARRAY, MOLWT, NAME )
!
!******************************************************************************
!  Subroutine TOTAL_BIOMASS_TG prints the amount of biomass burning emissions 
!  that are emitted each month in Tg or Tg C. (bmy, 3/20/01, 3/14/03)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Biomass burning CO emissions [molec/cm2/month]
!
!  NOTES:
!  (1 ) BBARRAY is now dimensioned (IIPAR,JJPAR).  Also, DXYP is dimensioned
!        as JGLOB, so use J+J0 to reference it. (bmy, 9/28/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (3 ) Now use function GET_AREA_CM2 from "grid_mod.f" to compute grid
!        box surface area in cm2.  Removed reference to CMN header file.
!        Cosmetic changes. (bmy, 3/14/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_CM2

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      REAL*8,           INTENT(IN) :: BBARRAY(IIPAR,JJPAR) 
      REAL*8,           INTENT(IN) :: MOLWT
      CHARACTER(LEN=*), INTENT(IN) :: NAME

      ! Local variables
      INTEGER                      :: I, J
      REAL*8                       :: TOTAL, A_CM2
      CHARACTER(LEN=6)             :: UNIT

      !=================================================================
      ! TOTAL_BIOMASS_TG begins here!
      !=================================================================

      ! Initialize summing variable
      TOTAL = 0d0

      ! Convert from [molec  /cm2/month] to [kg  /month]
      ! or      from [molec C/cm2/month] to [kg C/month]
      DO J = 1, JJPAR
         
         ! Grid box surface area [cm2]
         A_CM2 = GET_AREA_CM2( J )

         DO I = 1, IIPAR
            TOTAL = TOTAL + BBARRAY(I,J) * A_CM2 * ( MOLWT / 6.023d23 )
         ENDDO
      ENDDO
     
      ! Convert from kg --> Tg
      TOTAL = TOTAL * 1d-9

      ! Define unit string
      IF ( NAME == 'NOx' .or. NAME == 'CO' .or. NAME == 'CH2O' ) THEN
         UNIT = '[Tg  ]'
      ELSE
         UNIT = '[Tg C]'
      ENDIF

      ! Write totals
      WRITE( 6, 100 ) NAME, TOTAL, UNIT
 100  FORMAT( 'Sum Biomass ', a4, 1x, ': ', f9.3, 1x, a9  )
      ! Return to calling program
      END SUBROUTINE TOTAL_BIOMASS_TG

!------------------------------------------------------------------------------

      SUBROUTINE INIT_CO2 
!
!******************************************************************************
!  Subroutine INIT_CO2 allocates memory to module arrays and reads in annual
!  mean emissions. (pns, bmy, 8/16/05)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      LOGICAL, SAVE :: IS_INIT = .FALSE.
      INTEGER       :: AS 

      !=================================================================
      ! INIT_CO2 begins here!
      !=================================================================
           
      ! Exit if we have already intialized 
      IF ( IS_INIT ) RETURN

      ! Array for Fossil fuel CO2
      ALLOCATE( EMFOSSCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMFOSSCO2' )         
      EMFOSSCO2 = 0d0 

      ! Array for CO2 from ocean exchange
      ALLOCATE( EMOCCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMOCCO2' )         
      EMOCCO2 = 0d0 

      ! Array for Balanced Bio CO2
      ALLOCATE( EMBIOCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOCO2' )         
      EMBIOCO2 = 0d0 

      ! Array for Biomass burning CO2
      ALLOCATE( EMBIOBRNCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOBRNCO2' )
      EMBIOBRNCO2 = 0d0

      ! Array for Biofuel CO2
      ALLOCATE( EMBIOFUELCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIOFUELCO2' )         
      EMBIOFUELCO2 = 0d0 

      ! Array for NET BIO CO2
      ALLOCATE( EMBIONETCO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMBIONETCO2' )
      EMBIONETCO2  = 0d0

      !=================================================================
      ! Read in annual mean emissions
      !=================================================================

      ! Fossil fuel emissions
      IF ( LFOSSCO2    ) CALL READ_ANNUAL_FOSSILCO2

      ! Oceanic exchange
      IF ( LOCCO2      ) CALL READ_ANNUAL_OCEANCO2

      ! Biofuel emissions
      IF ( LBIOFUELCO2 ) CALL READ_ANNUAL_BIOFUELCO2

      ! Net terrestrial exchange
      IF ( LBIONETCO2  ) CALL READ_ANNUAL_BIONET_CO2

      ! Reset IS_INIT flag
      IS_INIT = .TRUE.
      
      ! Return to calling program
      END SUBROUTINE INIT_CO2 

!------------------------------------------------------------------------------
  
      SUBROUTINE CLEANUP_CO2 
!
!******************************************************************************
!  Subroutine CLEANUP_CO2 deallocates all module arrays (pns, bmy, 8/16/05)
! 
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_CO2 begins here!
      !=================================================================
      IF ( ALLOCATED( EMFOSSCO2    ) ) DEALLOCATE( EMFOSSCO2    )
      IF ( ALLOCATED( EMOCCO2      ) ) DEALLOCATE( EMOCCO2      )
      IF ( ALLOCATED( EMBIOCO2     ) ) DEALLOCATE( EMBIOCO2     )
      IF ( ALLOCATED( EMBIOBRNCO2  ) ) DEALLOCATE( EMBIOBRNCO2  )
      IF ( ALLOCATED( EMBIOFUELCO2 ) ) DEALLOCATE( EMBIOFUELCO2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_CO2

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CO2_MOD
