! $Id: future_emissions_mod.f,v 1.1 2009/09/16 14:06:28 bmy Exp $
      MODULE FUTURE_EMISSIONS_MOD
!
!******************************************************************************
!  Module FUTURE_EMISSIONS_MOD contains variables and routines for returning
!  scale factors for IPCC A1 & B1 emissions scenarios for future years
!  such as 2030, 2050, etc. (swu, bmy, 5/30/06)
!
!  The baseline year for the IPCC scale factors is 1995.  In other words, we
!  compute 1995 emissions in GEOS-Chem and then multiply e.g. 2030/1995 scale 
!  factors to compute the emissions for the future year.
!
!  Module Variables:
!  ============================================================================
!  (1 ) FUTURE_YEAR              : Year of future emissions (e.g. 2030, 2050)
!  (2 ) SCENARIO                 : IPCC emissions scenario (e.g. A1, B1)
!  (3 ) ALK4ff                   : Future scale factor array for fos-fuel ALK4
!  (4 ) BCbb                     : Future scale factor array for biomass  BC
!  (5 ) BCbf                     : Future scale factor array for biofuel  BC
!  (6 ) BCff                     : Future scale factor array for fos-fuel BC
!  (7 ) C2H6ff                   : Future scale factor array for fos-fuel C2H6
!  (8 ) C3H8ff                   : Future scale factor array for fos-fuel C3H8
!  (9 ) CObb                     : Future scale factor array for biomass  CO
!  (10) CObf                     : Future scale factor array for biofuel  CO
!  (11) COff                     : Future scale factor array for fos-fuel CO
!  (12) NH3an                    : Future scale factor array for anthro   NH3
!  (13) NH3bb                    : Future scale factor array for biomass  NH3
!  (14) NH3bf                    : Future scale factor array for biofuel  NH3
!  (15) NOxbb                    : Future scale factor array for biomass  NOx
!  (16) NOxbf                    : Future scale factor array for biofuel  NOx
!  (17) NOxff                    : Future scale factor array for fos-fuel NOx
!  (18) NOxft                    : Future scale factor array for fertiliz NOx
!  (19) OCbb                     : Future scale factor array for biomass  OC
!  (20) OCbf                     : Future scale factor array for biofuel  OC
!  (21) OCff                     : Future scale factor array for fos-fuel OC
!  (22) PRPEff                   : Future scale factor array for fossil   PRPE
!  (23) SO2bb                    : Future scale factor array for biomass  SO2
!  (24) SO2bf                    : Future scale factor array for biofuel  SO2
!  (25) SO2ff                    : Future scale factor array for fos-fuel SO2
!  (26) TONEff                   : Future scale factor array for fos-fuel TONE
!  (27) VOCbb                    : Future scale factor array for biomass  VOC
!  (28) VOCbf                    : Future scale factor array for biofuel  VOC
!  (29) VOCff                    : Future scale factor array for fos-fuel VOC
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_FUTURE_EMISSIONS      : Driver routine
!  (2 ) READ_GROWTH_FACTORS      : Reads future scale factors from disk
!  (3 ) GET_FUTURE_YEAR          : Returns the future emissions year
!  (4 ) GET_FUTURE_SCENARIO      : Returns the future emissions scenario
!  (5 ) GET_FUTURE_SCALE_ALK4ff  : Returns future fos-fuel ALK4 scale factors
!  (6 ) GET_FUTURE_SCALE_BCbb    : Returns future biomass  BC   scale factors
!  (7 ) GET_FUTURE_SCALE_BCbf    : Returns future biofuel  BC   scale factors
!  (8 ) GET_FUTURE_SCALE_BCff    : Returns future fos-fuel BC   scale factors
!  (9 ) GET_FUTURE_SCALE_C2H6ff  : Returns future fos-fuel C2H6 scale factors
!  (10) GET_FUTURE_SCALE_C3H8ff  : Returns future fos-fuel C3H8 scale factors
!  (11) GET_FUTURE_SCALE_CObb    : Returns future biomass  CO   scale factors
!  (12) GET_FUTURE_SCALE_CObf    : Returns future biofuel  CO   scale factors
!  (13) GET_FUTURE_SCALE_COff    : Returns future fos-fuel CO   scale factors
!  (14) GET_FUTURE_SCALE_NH3an   : Returns future anthro   NH3  scale factors
!  (15) GET_FUTURE_SCALE_NH3bb   : Returns future biomass  NH3  scale factors
!  (16) GET_FUTURE_SCALE_NH3bf   : Returns future biofuel  NH3  scale factors
!  (17) GET_FUTURE_SCALE_NOxbb   : Returns future biomass  NOx  scale factors
!  (18) GET_FUTURE_SCALE_NOxbf   : Returns future biofuel  NOx  scale factors
!  (19) GET_FUTURE_SCALE_NOxff   : Returns future fos-fuel NOx  scale factors
!  (20) GET_FUTURE_SCALE_NOxft   : Returns future fertiliz NOx  scale factors
!  (21) GET_FUTURE_SCALE_OCbb    : Returns future biomass  OC   scale factors
!  (22) GET_FUTURE_SCALE_OCbf    : Returns future biofuel  OC   scale factors
!  (23) GET_FUTURE_SCALE_OCff    : Returns future fos-fuel OC   scale factors
!  (24) GET_FUTURE_SCALE_PRPEff  : Returns future fos-fuel PRPE scale factors
!  (25) GET_FUTURE_SCALE_SO2bb   : Returns future biomass  SO2  scale factors
!  (26) GET_FUTURE_SCALE_SO2bf   : Returns future biofuel  SO2  scale factors
!  (27) GET_FUTURE_SCALE_SO2ff   : Returns future fos-fuel SO2  scale factors
!  (28) GET_FUTURE_SCALE_TONEff  : Returns future fos-fuel TONE scale factors
!  (29) GET_FUTURE_SCALE_VOCbb   : Returns future biomass  VOC  scale factors
!  (30) GET_FUTURE_SCALE_VOCbf   : Returns future biofuel  VOC  scale factors
!  (31) GET_FUTURE_SCALE_VOCff   : Returns future fos-fuel VOC  scale factors
!  (32) INIT_FUTURE_EMISSIONS    : Initializes and allocates module arrays
!  (33) CLEANUP_FUTURE_EMISSIONS : Deallocates all module arrays
!
!  GEOS-Chem modules referenced by "future_emissions_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f              : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f          : Module w/ GEOS-CHEM met field and data dirs
!  (3 ) error_mod.f              : Module w/ I/O error and NaN check routines
!  (4 ) file_mod.f               : Module w/ file unit numbers & error checks
!  (5 ) transfer_mod.f           : Module w/ routines to cast & resize arrays
!
!  References:
!  ============================================================================
!
!  NOTES:
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "future_emissions_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE
      
      ! ... except these routines
      PUBLIC :: CLEANUP_FUTURE_EMISSIONS
      PUBLIC :: DO_FUTURE_EMISSIONS
      PUBLIC :: GET_FUTURE_YEAR
      PUBLIC :: GET_FUTURE_SCENARIO
      PUBLIC :: GET_FUTURE_SCALE_ALK4ff  
      PUBLIC :: GET_FUTURE_SCALE_BCbb  
      PUBLIC :: GET_FUTURE_SCALE_BCbf  
      PUBLIC :: GET_FUTURE_SCALE_BCff  
      PUBLIC :: GET_FUTURE_SCALE_C2H6ff
      PUBLIC :: GET_FUTURE_SCALE_C3H8ff
      PUBLIC :: GET_FUTURE_SCALE_CObb  
      PUBLIC :: GET_FUTURE_SCALE_CObf  
      PUBLIC :: GET_FUTURE_SCALE_COff  
      PUBLIC :: GET_FUTURE_SCALE_NH3an 
      PUBLIC :: GET_FUTURE_SCALE_NH3bb 
      PUBLIC :: GET_FUTURE_SCALE_NH3bf 
      PUBLIC :: GET_FUTURE_SCALE_NOxbb 
      PUBLIC :: GET_FUTURE_SCALE_NOxbf 
      PUBLIC :: GET_FUTURE_SCALE_NOxff 
      PUBLIC :: GET_FUTURE_SCALE_NOxft 
      PUBLIC :: GET_FUTURE_SCALE_OCbb  
      PUBLIC :: GET_FUTURE_SCALE_OCbf  
      PUBLIC :: GET_FUTURE_SCALE_OCff  
      PUBLIC :: GET_FUTURE_SCALE_PRPEff
      PUBLIC :: GET_FUTURE_SCALE_SO2bb 
      PUBLIC :: GET_FUTURE_SCALE_SO2bf 
      PUBLIC :: GET_FUTURE_SCALE_SO2ff 
      PUBLIC :: GET_FUTURE_SCALE_TONEff 
      PUBLIC :: GET_FUTURE_SCALE_VOCbb 
      PUBLIC :: GET_FUTURE_SCALE_VOCbf 
      PUBLIC :: GET_FUTURE_SCALE_VOCff 
      PUBLIC :: INIT_FUTURE_EMISSIONS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER             :: FUTURE_YEAR
      CHARACTER(LEN=2)    :: SCENARIO

      ! Arrays
      REAL*8, ALLOCATABLE :: ALK4ff(:,:)  
      REAL*8, ALLOCATABLE :: BCbb(:,:) 
      REAL*8, ALLOCATABLE :: BCbf(:,:) 
      REAL*8, ALLOCATABLE :: BCff(:,:) 
      REAL*8, ALLOCATABLE :: C2H6ff(:,:) 
      REAL*8, ALLOCATABLE :: C3H8ff(:,:) 
      REAL*8, ALLOCATABLE :: CObb(:,:) 
      REAL*8, ALLOCATABLE :: CObf(:,:) 
      REAL*8, ALLOCATABLE :: COff(:,:) 
      REAL*8, ALLOCATABLE :: NH3an(:,:)
      REAL*8, ALLOCATABLE :: NH3bb(:,:)
      REAL*8, ALLOCATABLE :: NH3bf(:,:)
      REAL*8, ALLOCATABLE :: NOxbb(:,:) 
      REAL*8, ALLOCATABLE :: NOxbf(:,:) 
      REAL*8, ALLOCATABLE :: NOxff(:,:) 
      REAL*8, ALLOCATABLE :: NOxft(:,:) 
      REAL*8, ALLOCATABLE :: OCbb(:,:) 
      REAL*8, ALLOCATABLE :: OCbf(:,:) 
      REAL*8, ALLOCATABLE :: OCff(:,:) 
      REAL*8, ALLOCATABLE :: PRPEff(:,:)
      REAL*8, ALLOCATABLE :: TONEff(:,:) 
      REAL*8, ALLOCATABLE :: SO2bb(:,:)
      REAL*8, ALLOCATABLE :: SO2bf(:,:)
      REAL*8, ALLOCATABLE :: SO2ff(:,:)
      REAL*8, ALLOCATABLE :: VOCbb(:,:) 
      REAL*8, ALLOCATABLE :: VOCbf(:,:) 
      REAL*8, ALLOCATABLE :: VOCff(:,:) 

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_FUTURE_EMISSIONS( THIS_YEAR, THIS_SCEN )
!
!******************************************************************************
!  Subroutine DO_FUTURE_EMISSIONS reads future emission growth factors
!  into module arrays.  This can be done once at the beginning of the
!  GEOS-Chem simulation. (swu, bmy, 5/30/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THIS_FY   (INTEGER)   : Year for future emission growth factors
!  (2 ) THIS_SCEN (CHARACTER) : Emissions scenario (e.g. "A1", "B1", etc.)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER,          INTENT(IN) :: THIS_YEAR
      CHARACTER(LEN=*), INTENT(IN) :: THIS_SCEN

      ! Local variables
      LOGICAL                      :: FIRST = .TRUE.

      !=================================================================
      ! DO_FUTURE_EMISSIONS begins here!
      !=================================================================
      
      ! Save module variables
      FUTURE_YEAR = THIS_YEAR
      SCENARIO    = TRIM( THIS_SCEN )

      ! First-time initialization
      IF ( FIRST ) THEN 
         CALL INIT_FUTURE_EMISSIONS
         FIRST = .FALSE.
      ENDIF

      !----------------------
      ! Read growth factors
      !----------------------

      ! ALK4
      CALL READ_GROWTH_FACTORS( 'ALK4FF', 5,  ALK4ff )  

      ! BC
      CALL READ_GROWTH_FACTORS( 'BCBB',   34, BCbb   ) 
      CALL READ_GROWTH_FACTORS( 'BCBF',   34, BCbf   ) 
      CALL READ_GROWTH_FACTORS( 'BCFF',   34, BCff   ) 

      ! C2H6
      CALL READ_GROWTH_FACTORS( 'C2H6FF', 21, C2H6ff ) 

      ! C3H8
      CALL READ_GROWTH_FACTORS( 'C3H8FF', 19, C3H8ff ) 

      ! CO
      CALL READ_GROWTH_FACTORS( 'COBB',   4,  CObb   ) 
      CALL READ_GROWTH_FACTORS( 'COBF',   4,  CObf   ) 
      CALL READ_GROWTH_FACTORS( 'COFF',   4,  COff   ) 

      ! NH3
      CALL READ_GROWTH_FACTORS( 'NH3AN',  30, NH3an  )
      CALL READ_GROWTH_FACTORS( 'NH3BB',  30, NH3bb  )
      CALL READ_GROWTH_FACTORS( 'NH3BF',  30, NH3bf  )

      ! NOx
      CALL READ_GROWTH_FACTORS( 'NOxBB',  1,  NOxbb  )
      CALL READ_GROWTH_FACTORS( 'NOxBF',  1,  NOxbf  )
      CALL READ_GROWTH_FACTORS( 'NOxFF',  1,  NOxff  )
      CALL READ_GROWTH_FACTORS( 'NOxFT',  1,  NOxft  )

      ! OC
      CALL READ_GROWTH_FACTORS( 'OCBB',   35, OCbb   ) 
      CALL READ_GROWTH_FACTORS( 'OCBF',   35, OCbf   ) 
      CALL READ_GROWTH_FACTORS( 'OCFF',   35, OCff   ) 

      ! PRPE
      CALL READ_GROWTH_FACTORS( 'PRPEFF', 18, PRPEff )

      ! TONE (Ketones > C3; use for ACET, MEK)
      CALL READ_GROWTH_FACTORS( 'TONEFF', 9,  TONEff ) 

      ! SO2
      CALL READ_GROWTH_FACTORS( 'SO2BB',  26, SO2bb  )
      CALL READ_GROWTH_FACTORS( 'SO2BF',  26, SO2bf  )
      CALL READ_GROWTH_FACTORS( 'SO2FF',  26, SO2ff  )

      ! VOC
      CALL READ_GROWTH_FACTORS( 'VOCBB',  90, VOCbb  ) 
      CALL READ_GROWTH_FACTORS( 'VOCBF',  90, VOCbf  ) 
      CALL READ_GROWTH_FACTORS( 'VOCFF',  90, VOCff  ) 

      !----------------------
      ! Print ranges 
      !----------------------

      ! Write header
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F U T U R E   E M I S S I O N   ' //
     &                  'G R O W T H   F A C T O R S'
      
      ! Write year and scenario
      WRITE( 6, 100 ) FUTURE_YEAR, SCENARIO

      ! Write totals
      WRITE( 6, 110 ) 'ALK4ff', MINVAL( ALK4ff ), MAXVAL( ALK4ff )
      WRITE( 6, 110 ) 'BCbb',   MINVAL( BCbb   ), MAXVAL( BCbb   )
      WRITE( 6, 110 ) 'BCbf',   MINVAL( BCbf   ), MAXVAL( BCbf   )
      WRITE( 6, 110 ) 'BCff',   MINVAL( BCff   ), MAXVAL( BCff   )
      WRITE( 6, 110 ) 'C2H6ff', MINVAL( C2H6ff ), MAXVAL( C2H6ff )
      WRITE( 6, 110 ) 'C3H8ff', MINVAL( C3H8ff ), MAXVAL( C3H8ff )
      WRITE( 6, 110 ) 'CObb',   MINVAL( CObb   ), MAXVAL( CObb   )
      WRITE( 6, 110 ) 'CObf',   MINVAL( CObf   ), MAXVAL( CObf   )
      WRITE( 6, 110 ) 'COff',   MINVAL( COff   ), MAXVAL( COff   )
      WRITE( 6, 110 ) 'NH3an',  MINVAL( NH3an  ), MAXVAL( NH3an  )
      WRITE( 6, 110 ) 'NH3bb',  MINVAL( NH3bb  ), MAXVAL( NH3bb  )
      WRITE( 6, 110 ) 'NH3bf',  MINVAL( NH3bf  ), MAXVAL( NH3bf  )
      WRITE( 6, 110 ) 'NOxbb',  MINVAL( NOxbb  ), MAXVAL( NOxbb  )
      WRITE( 6, 110 ) 'NOxbf',  MINVAL( NOxbf  ), MAXVAL( NOxbf  )
      WRITE( 6, 110 ) 'NOxff',  MINVAL( NOxff  ), MAXVAL( NOxff  )
      WRITE( 6, 110 ) 'NOxft',  MINVAL( NOxft  ), MAXVAL( NOxft  )
      WRITE( 6, 110 ) 'OCbb',   MINVAL( OCbb   ), MAXVAL( OCbb   )
      WRITE( 6, 110 ) 'OCbf',   MINVAL( OCbf   ), MAXVAL( OCbf   )
      WRITE( 6, 110 ) 'OCff',   MINVAL( OCff   ), MAXVAL( OCff   )
      WRITE( 6, 110 ) 'PRPEff', MINVAL( PRPEff ), MAXVAL( PRPEff )
      WRITE( 6, 110 ) 'TONEff', MINVAL( TONEff ), MAXVAL( TONEff )
      WRITE( 6, 110 ) 'SO2bb',  MINVAL( SO2bb  ), MAXVAL( SO2bb  )
      WRITE( 6, 110 ) 'SO2bf',  MINVAL( SO2bf  ), MAXVAL( SO2bf  )
      WRITE( 6, 110 ) 'SO2ff',  MINVAL( SO2ff  ), MAXVAL( SO2ff  )
      WRITE( 6, 110 ) 'VOCbb',  MINVAL( VOCbb  ), MAXVAL( VOCbb  )
      WRITE( 6, 110 ) 'VOCbf',  MINVAL( VOCbf  ), MAXVAL( VOCbf  )
      WRITE( 6, 110 ) 'VOCff',  MINVAL( VOCff  ), MAXVAL( VOCff  )

      ! Write footer
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )

      ! FORMAT statements
 100  FORMAT( 'for year ', i4, ' and emissions scenario ', a2, / )
 110  FORMAT( a6, ' growth factors range from ', f8.3, ' to ', f8.3 )

      ! Return to calling program
      END SUBROUTINE DO_FUTURE_EMISSIONS

!------------------------------------------------------------------------------

      SUBROUTINE READ_GROWTH_FACTORS( SPECIES, TRACER, GRFACTORS )
!
!******************************************************************************
!  Subroutine READ_GROWTH_FACTORS reads the future growth factors for one
!  species, future year, and scenario from disk. (swu, bmy, 5/30/06)
!
!  If growth factors for a particular species do not exist for a given future
!  year and scenario, then the READ_GROWTH_FACTORS will return and the
!  GRFACTORS array will be set to 1 everywhere.
!
!  The baseline year for the IPCC scale factors is 1995.  In other words, we
!  compute 1995 emissions in GEOS-Chem and then multiply e.g. 2030/1995 scale 
!  factors to compute the emissions for the future year.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SPECIES   (CHARACTER) : Species name to read (e.g. "NOxbb, CObb", etc)
!  (2 ) TRACER    (INTEGER  ) : Tracer number 
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) GRFACTORS (REAL*8   ) : Array of growth factors
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,            ONLY : GET_RES_EXT, GET_NAME_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0,    READ_BPCH2
      USE DIRECTORY_MOD,        ONLY : DATA_DIR
      USE FILE_MOD,             ONLY : FILE_EXISTS
      USE TRANSFER_MOD,         ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD             ! Size parameters

      INTEGER,          INTENT(IN)  :: TRACER
      CHARACTER(LEN=*), INTENT(IN)  :: SPECIES
      REAL*8,           INTENT(OUT) :: GRFACTORS(IIPAR,JJPAR)

      ! Local variables
      INTEGER                       :: IOS, IUNIT
      REAL*4                        :: ARRAY(IIPAR,JJPAR,1)
      REAL*8                        :: TAU0
      CHARACTER(LEN=4)              :: YSTR
      CHARACTER(LEN=255)            :: FILENAME

      !=================================================================
      ! READ_GROWTH_FACTORS begins here!
      !=================================================================

      ! Initialize
      GRFACTORS(:,:) = 1d0

      ! Create a string for the 4-digit year
      WRITE( YSTR, '(i4)' ) FUTURE_YEAR

      ! File name
      FILENAME = TRIM( DATA_DIR )           // 
     &           'future_emissions_200605/' // YSTR           //
     &           '/'                        // SCENARIO       //
     &           '/'                        // SPECIES        // 
     &           '_'                        // SCENARIO       // 
     &           '.'                        // GET_NAME_EXT() // 
     &           '.'                        // GET_RES_EXT()  // 
     &           '.'                        // YSTR

      ! Return if file is not found (growth factors array = 1)
      IF ( .not. FILE_EXISTS( FILENAME ) ) RETURN

      ! TAU0 value for the future year
      TAU0 = GET_TAU0( 1, 1, FUTURE_YEAR )

      ! ACET is stored in the biomass file as tracer #9
      CALL READ_BPCH2( FILENAME, 'FUTURE-E', TRACER, 
     &                 TAU0,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. ) 

      ! Cast to REAL*8
      CALL TRANSFER_2D( ARRAY(:,:,1), GRFACTORS )

      ! Return to calling program
      END SUBROUTINE READ_GROWTH_FACTORS

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_YEAR() RESULT( THIS_YEAR )
!
!******************************************************************************
!  Function GET_FUTURE_YEAR returns the year for future emissions
!  to the calling program. (swu, bmy, 5/30/06)
!
!  NOTES:
!******************************************************************************
! 
      ! Function value
      REAL*8 :: THIS_YEAR

      !=================================================================
      ! GET_FUTURE_YEAR begins here!
      !=================================================================
      THIS_YEAR = FUTURE_YEAR 
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_YEAR

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCENARIO() RESULT( THIS_SCEN )
!
!******************************************************************************
!  Function GET_FUTURE_SCENARIO returns the IPCC future emissions scenario
!  (e.g. A1, B1) for future emissions to the calling program. 
!  (swu, bmy, 5/30/06)
!
!  NOTES:
!******************************************************************************
! 
      ! Function value
      CHARACTER(LEN=255) :: THIS_SCEN

      !=================================================================
      ! GET_FUTURE_SCENARIO begins here!
      !=================================================================
      THIS_SCEN = TRIM( SCENARIO )
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCENARIO

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_ALK4ff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_ALK4ff returns the future scale factor for
!  Fossil Fuel ALK4 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_ALK4ff begins here!
      !=================================================================
      SCALEFAC = ALK4ff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_ALK4ff  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_BCbb( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_BCbb returns the future scale factor for
!  biomass burning BC for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_BCbb begins here!
      !=================================================================
      SCALEFAC = BCbb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_BCbb   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_BCbf( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_BCbf returns the future scale factor for
!  biofuel BC for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_BCbf begins here!
      !=================================================================
      SCALEFAC = BCbf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_BCbf   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_BCff( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_BCff returns the future scale factor for
!  Fossil Fuel BC for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_BCff begins here!
      !=================================================================
      SCALEFAC = BCff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_BCff   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_C2H6ff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_C2H6ff returns the future scale factor for
!  Fossil Fuel C2H6 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_C2H6ff begins here!
      !=================================================================
      SCALEFAC = C2H6ff( I, J )
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_C2H6ff 

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_C3H8ff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_C3H8ff returns the future scale factor for
!  Fossil Fuel C3H8 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_C3H8ff begins here!
      !=================================================================
      SCALEFAC = C3H8ff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_C3H8ff 

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_CObb( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_CObb returns the future scale factor for
!  biomass burning CO for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_CObb begins here!
      !=================================================================
      SCALEFAC = CObb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_CObb   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_CObf( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_CObf returns the future scale factor for
!  biofuel CO for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_CObf begins here!
      !=================================================================
      SCALEFAC = CObf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_CObf   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_COff( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_COff returns the future scale factor for
!  Fossil Fuel CO for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_COff begins here!
      !=================================================================
      SCALEFAC = COff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_COff   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NH3an( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NH3an returns the future scale factor for
!  anthropogenic NH3 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NH3an begins here!
      !=================================================================
      SCALEFAC = NH3an(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NH3an  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NH3bb( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NH3bb returns the future scale factor for
!  biomass burning NH3 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NH3bb begins here!
      !=================================================================
      SCALEFAC = NH3bb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NH3bb  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NH3bf( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NH3bf returns the future scale factor for
!  biofuel NH3 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NH3bf begins here!
      !=================================================================
      SCALEFAC = NH3bf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NH3bf  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NOxbb( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NOxbb returns the future scale factor for
!  biomass burning NOx for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NOxbb begins here!
      !=================================================================
      SCALEFAC = NOXbb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NOxbb  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NOxbf( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NOXbf returns the future scale factor for
!  biofuel NOx for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NOxbf begins here!
      !=================================================================
      SCALEFAC = NOxbf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NOxbf  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NOxff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NOxff returns the future scale factor for
!  Fossil Fuel NOx for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NOxff begins here!
      !=================================================================
      SCALEFAC = NOxff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NOxff  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_NOxft( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_NOxft returns the future scale factor 
!  for NOx from the free tropoposphere the GEOS-Chem grid box (I,J) 
!  (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_NOxft begins here!
      !=================================================================
      SCALEFAC = NOxft(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_NOxft

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_OCbb( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_OCbb returns the future scale factor for
!  biomass burning OC for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_OCbb begins here!
      !=================================================================
      SCALEFAC = OCbb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_OCbb   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_OCbf( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_OCbf returns the future scale factor for
!  biofuel OC for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_OCbf begins here!
      !=================================================================
      SCALEFAC = OCbf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_OCbf   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_OCff( I, J ) RESULT( SCALEFAC )   
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_OCff returns the future scale factor for
!  Fossil Fuel ACET for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_OCff begins here!
      !=================================================================
      SCALEFAC = OCff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_OCff   

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_PRPEff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_PRPEff returns the future scale factor for
!  Fossil Fuel PRPE for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_PRPEff begins here!
      !=================================================================
      SCALEFAC = PRPEff( I, J )
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_PRPEff 

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_TONEff( I, J ) RESULT( SCALEFAC )
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_TONEff returns the future scale factor for
!  Fossil Fuel TONE (Ketones > C3, such as ACET, MEK) for the GEOS-Chem grid 
!  box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_TONEff begins here!
      !=================================================================
      SCALEFAC = TONEff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_TONEff

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_SO2bb( I, J ) RESULT( SCALEFAC )
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_SO2bb returns the future scale factor for
!  Fossil Fuel ACET for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_SO2bb begins here!
      !=================================================================
      SCALEFAC = SO2bb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_SO2bb 

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_SO2bf( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_SO2bf returns the future scale factor for
!  biofuel SO2 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_SO2bf begins here!
      !=================================================================
      SCALEFAC = SO2bf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_SO2bf  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_SO2ff( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_SO2ff returns the future scale factor for
!  Fossil Fuel SO2 for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_SO2ff begins here!
      !=================================================================
      SCALEFAC = SO2ff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_SO2ff  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_VOCbb( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_VOCbb returns the future scale factor for
!  biomass burning VOC's for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_VOCbb begins here!
      !=================================================================
      SCALEFAC = VOCbb(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_VOCbb  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_VOCbf( I, J ) RESULT( SCALEFAC )  
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_VOCbf returns the future scale factor for
!  biofuel VOC's for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_VOCbf begins here!
      !=================================================================
      SCALEFAC = VOCbf(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_VOCbf  

!------------------------------------------------------------------------------

      FUNCTION GET_FUTURE_SCALE_VOCff( I, J ) RESULT( SCALEFAC ) 
!
!******************************************************************************
!  Function GET_FUTURE_SCALE_VOCff returns the future scale factor for
!  Fossil Fuel VOC's for the GEOS-Chem grid box (I,J) (swu, bmy, 5/30/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) :: GEOS-Chem longitude index
!  (2 ) J (INTEGER) :: GEOS-Chem latitude  index
!
!  NOTES:
!******************************************************************************
! 
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Function variable
      REAL*8              :: SCALEFAC

      !=================================================================
      ! GET_FUTURE_SCALE_VOCff begins here!
      !=================================================================
      SCALEFAC = VOCff(I,J)
      
      ! Return to calling program
      END FUNCTION GET_FUTURE_SCALE_VOCff  

!------------------------------------------------------------------------------

      SUBROUTINE INIT_FUTURE_EMISSIONS
!
!******************************************************************************
!  Subroutine CLEANUP_FUTURE_EMISSIONS allocates and initializes all module 
!  arrays. (swu, bmy, 5/30/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD  ! Size parameters

      ! Local variables
      INTEGER            :: AS

      !=================================================================
      ! INIT_FUTURE_EMISSIONS begins here!
      !=================================================================

      ALLOCATE( ALK4ff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'ALK4ff' )
      ALK4ff = 1d0

      ALLOCATE( BCbb( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'BCbb' )
      BCbb = 1d0

      ALLOCATE( BCbf( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( '' )
      BCbf = 1d0

      ALLOCATE( BCff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'BCff' )
      BCff = 1d0

      ALLOCATE( C2H6ff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'C2H6ff' )
      C2H6ff = 1d0

      ALLOCATE( C3H8ff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'CH38ff' )
      C3H8ff = 1d0

      ALLOCATE( CObb( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'CObb' )
      CObb = 1d0

      ALLOCATE( CObf( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'CObf' )
      CObf = 1d0

      ALLOCATE( COff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'COff' )
      COff = 1d0

      ALLOCATE( NH3an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NH3an' )
      NH3an = 1d0

      ALLOCATE( NH3bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NH3bb' )
      NH3bb = 1d0

      ALLOCATE( NH3bf( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NH3bf' )
      NH3bf = 1d0

      ALLOCATE( NOxbb( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NOxbb' )
      NOXbb = 1d0

      ALLOCATE( NOxbf( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NOxbf' )
      NOxbf = 1d0

      ALLOCATE( NOxff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NOxff' )
      NOxff = 1d0

      ALLOCATE( NOxft ( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'NOxft' )
      NOxft = 1d0

      ALLOCATE( OCbb( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'OCbb' )
      OCbb = 1d0

      ALLOCATE( OCbf( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'OCbf' )
      OCbf = 1d0

      ALLOCATE( OCff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'OCff' )
      OCff = 1d0

      ALLOCATE( PRPEff( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'PRPEff' )
      PRPEff = 1d0

      ALLOCATE( TONEff( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'TONEff' )
      TONEff = 1d0

      ALLOCATE( SO2bb( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'SO2bb' )
      SO2bb = 1d0

      ALLOCATE( SO2bf( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'SO2bf' )
      SO2bf = 1d0

      ALLOCATE( SO2ff( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /=0 ) CALL ALLOC_ERR( 'SO2ff' )
      SO2ff = 1d0

      ALLOCATE( VOCbb( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'VOCbb' )
      VOCbb = 1d0

      ALLOCATE( VOCbf ( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'VOCbf' )
      VOCbf = 1d0

      ALLOCATE( VOCff ( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /=0 ) CALL ALLOC_ERR( 'VOCff' )
      VOCff = 1d0

      ! Return to calling program
      END SUBROUTINE INIT_FUTURE_EMISSIONS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_FUTURE_EMISSIONS
!
!******************************************************************************
!  Subroutine CLEANUP_FUTURE_EMISSIONS deallocates all module arrays.
!  (swu, bmy, 5/30/06)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_FUTURE_EMISSIONS begins here!
      !=================================================================
      IF ( ALLOCATED( ALK4ff ) ) DEALLOCATE( ALK4ff )  
      IF ( ALLOCATED( BCbb   ) ) DEALLOCATE( BCbb   ) 
      IF ( ALLOCATED( BCbf   ) ) DEALLOCATE( BCbf   ) 
      IF ( ALLOCATED( BCff   ) ) DEALLOCATE( BCff   ) 
      IF ( ALLOCATED( C2H6ff ) ) DEALLOCATE( C2H6ff ) 
      IF ( ALLOCATED( C3H8ff ) ) DEALLOCATE( C3H8ff ) 
      IF ( ALLOCATED( CObb   ) ) DEALLOCATE( CObb   ) 
      IF ( ALLOCATED( CObf   ) ) DEALLOCATE( CObf   ) 
      IF ( ALLOCATED( COff   ) ) DEALLOCATE( COff   ) 
      IF ( ALLOCATED( NH3an  ) ) DEALLOCATE( NH3an  )
      IF ( ALLOCATED( NH3bb  ) ) DEALLOCATE( NH3bb  )
      IF ( ALLOCATED( NH3bf  ) ) DEALLOCATE( NH3bf  )
      IF ( ALLOCATED( NOxbb  ) ) DEALLOCATE( NOxbb  ) 
      IF ( ALLOCATED( NOxbf  ) ) DEALLOCATE( NOxbf  ) 
      IF ( ALLOCATED( NOxff  ) ) DEALLOCATE( NOxff  ) 
      IF ( ALLOCATED( NOxft  ) ) DEALLOCATE( NOxft  ) 
      IF ( ALLOCATED( OCbb   ) ) DEALLOCATE( OCbb   ) 
      IF ( ALLOCATED( OCbf   ) ) DEALLOCATE( OCbf   ) 
      IF ( ALLOCATED( OCff   ) ) DEALLOCATE( OCff   ) 
      IF ( ALLOCATED( PRPEff ) ) DEALLOCATE( PRPEff )
      IF ( ALLOCATED( TONEff ) ) DEALLOCATE( TONEff ) 
      IF ( ALLOCATED( SO2bb  ) ) DEALLOCATE( SO2bb  )
      IF ( ALLOCATED( SO2bf  ) ) DEALLOCATE( SO2bf  )
      IF ( ALLOCATED( SO2ff  ) ) DEALLOCATE( SO2ff  )
      IF ( ALLOCATED( VOCbb  ) ) DEALLOCATE( VOCbb  ) 
      IF ( ALLOCATED( VOCbf  ) ) DEALLOCATE( VOCbf  ) 
      IF ( ALLOCATED( VOCff  ) ) DEALLOCATE( VOCff  ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_FUTURE_EMISSIONS

!------------------------------------------------------------------------------

      ! End of module
      END MODULE FUTURE_EMISSIONS_MOD
