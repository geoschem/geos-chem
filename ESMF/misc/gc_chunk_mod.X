!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains 
!  the GEOS-Chem chunk code init, run and finalize methods.
!\\
!\\
! !INTERFACE: 
!      
      MODULE GC_CHUNK_MOD
!
! !USES:
!      
      USE GC_TYPE_MOD          

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!      PUBLIC  :: GC_CHUNK_INIT
      PUBLIC  :: GC_CHUNK_RUN
      PUBLIC  :: GC_CHUNK_FINAL
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: ITS_TIME
      PRIVATE :: CONVERT_UNITS
      PRIVATE :: PRINT_DIAG_COL
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  28 Oct 2009 - R. Yantosca - Various updates in module routines
!  14 Dec 2009 - R. Yantosca - GC_CHUNK_RUN now accepts TRACER_1d in [mol/mol]
!                              mixing ratio.  Apropos unit conversions are
!                              now done w/in GC_CHUNK_RUN internally.
!  14 Apr 2010 - R. Yantosca - Removed hardwireing from GC_CHUNK_INIT
!  08 Jul 2010 - R. Yantosca - Add GC_DIAG private type for diag printout
!  08 Jul 2010 - R. Yantosca - Add DIAG_COL object for diag printout
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Derived type for chunk diagnostic output (for code validation)
      TYPE GC_DIAG
         LOGICAL                    :: DO_PRINT     ! Should we print out?
         INTEGER                    :: N_DIAG       ! # of diag quantities
         INTEGER                    :: COUNT        ! Counter for averaging
         CHARACTER(LEN=10), POINTER :: NAME(:)      ! Tracer names
         REAL*8,            POINTER :: TRACER(:,:)  ! Tracer concentrations
         CHARACTER(LEN=40)          :: FILENAME     ! File name for output
         INTEGER                    :: LUN          ! File unit # for output
      END TYPE GC_DIAG

      ! Derived type object for saving concentration diagnostics
      TYPE(GC_DIAG)      :: DIAG_COL

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_chunk_run
!
! !DESCRIPTION: Routine GC\_CHUNK\_RUN is the driver for the following 
! operations:
! \begin{itemize}
! \item Planetary boundary layer mixing
! \item Cloud convection
! \item Dry deposition
! \item Emissions
! \item Chemistry
! \item Wet Depositon
! \end{itemize}
!
! !INTERFACE:
!
      SUBROUTINE GC_CHUNK_RUN( )

!      CALL GEOS_CHEM_
      END SUBROUTINE GC_CHUNK_RUN
!EOC
!v0.1!------------------------------------------------------------------------------
!v0.1!          Harvard University Atmospheric Chemistry Modeling Group            !
!v0.1!------------------------------------------------------------------------------
!v0.1!BOP
!v0.1!
!v0.1! !IROUTINE: gc_chunk_init
!v0.1!
!v0.1! !DESCRIPTION: Subroutine GC\_CHUNK\_INIT calls the various 
!v0.1!  initialization routines that read the setup files for the GEOS-Chem
!v0.1!  chunk code.  Also, ID flags for advected tracers, chemical species,
!v0.1!  dry deposition species and wet deposition species are defined.
!v0.1!\\
!v0.1!\\
!v0.1! !INTERFACE:
!v0.1!
!v0.1      SUBROUTINE GC_CHUNK_INIT( IDENT,       OPTIONS,    L_COLUMN,    
!v0.1     &                           N_TRACERS,   N_MEMBERS,  N_DRYDEP,    
!v0.1     &                           N_WETDEP,    TS_CHEM,    TRACER_NAME, 
!v0.1     &                           ID_SPECIES,  ID_TRACERS, COEF,        
!v0.1     &                           ID_DRYDEP,   ID_WETDEP,  N_SPECIES,
!v0.1     &                           N_REACTIONS, N_JV,       RC )   
!v0.1!
!v0.1! !USES:
!v0.1!
!v0.1      USE DRYDEP_MOD,  ONLY : INIT_DRYDEP
!v0.1      USE WETDEP_MOD,  ONLY : INIT_WETDEP
!v0.1      
!v0.1      IMPLICIT NONE
!v0.1
!v0.1#     include "smv_dimension.h"  ! Dimensions for common blocks
!v0.1#     include "comode_loop.h"    ! SMVGEAR common blocks
!v0.1#     include "smv_errcode.h"    ! Error codes
!v0.1#     include "smv_physconst.h"  ! Physical constants
!v0.1!
!v0.1! !INPUT PARAMETERS:
!v0.1!
!v0.1      ! Object with logical flags
!v0.1      TYPE(GC_OPTIONS),  INTENT(IN)    :: OPTIONS
!v0.1
!v0.1      ! Number of boxes in the atmospheric chunk
!v0.1      INTEGER,           INTENT(IN)    :: L_COLUMN
!v0.1
!v0.1      ! Number of advected tracers
!v0.1      INTEGER,           INTENT(IN)    :: N_TRACERS
!v0.1
!v0.1      ! Max # of species per chemical family
!v0.1      INTEGER,           INTENT(IN)    :: N_MEMBERS
!v0.1
!v0.1      ! Chemistry timestep [minutes]
!v0.1      REAL*8,            INTENT(IN)    :: TS_CHEM
!v0.1
!v0.1      ! Names of advected tracers
!v0.1      CHARACTER(LEN=*),  INTENT(IN)    :: TRACER_NAME(N_TRACERS)
!v0.1!
!v0.1! !INPUT/OUTPUT PARAMETERS:
!v0.1!
!v0.1      ! Object for Identification info from the Gridded Component
!v0.1      TYPE(GC_IDENT),    INTENT(INOUT) :: IDENT
!v0.1
!v0.1      ! Object for ID flags for SMVGEAR chemical species
!v0.1      TYPE(ID_SPEC),     INTENT(INOUT) :: ID_SPECIES
!v0.1
!v0.1      ! Object for ID flags for advected tracers
!v0.1      TYPE(ID_TRAC),     INTENT(INOUT) :: ID_TRACERS
!v0.1
!v0.1      ! Object for # of species per advected tracer, etc.
!v0.1      TYPE(SPEC_2_TRAC), INTENT(INOUT) :: COEF
!v0.1
!v0.1      ! Object for ID flags for dry-deposited and wet-deposited tracers
!v0.1      TYPE(ID_DRYD),     INTENT(INOUT) :: ID_DRYDEP
!v0.1      TYPE(ID_WETD),     INTENT(INOUT) :: ID_WETDEP
!v0.1!
!v0.1! !OUTPUT PARAMETERS:
!v0.1!
!v0.1      ! Number of dry deposited and wet-deposited species
!v0.1      INTEGER,           INTENT(OUT)   :: N_DRYDEP
!v0.1      INTEGER,           INTENT(OUT)   :: N_WETDEP
!v0.1
!v0.1      ! Number of chemical species, reactions, and photolysis reactions
!v0.1      INTEGER,           INTENT(OUT)   :: N_SPECIES
!v0.1      INTEGER,           INTENT(OUT)   :: N_REACTIONS
!v0.1      INTEGER,           INTENT(OUT)   :: N_JV
!v0.1
!v0.1      ! Return code
!v0.1      INTEGER,           INTENT(OUT)   :: RC 
!v0.1!
!v0.1! !REVISION HISTORY: 
!v0.1!  22 Jun 2009 - R. Yantosca - Initial version
!v0.1!  15 Jul 2009 - R. Yantosca - Now call init routines for drydep & wetdep
!v0.1!  14 Apr 2010 - R. Yantosca - Now de-hardwire definition of ID_TRACERS
!v0.1!  16 Apr 2010 - R. Yantosca - Now pass back N_SPECIES, N_REACTIONS as outputs
!v0.1!  23 Apr 2010 - R. Yantosca - Added IDENT to the argument list
!v0.1!  23 Apr 2010 - R. Yantosca - Redirect stdout to a log file
!v0.1!  23 Apr 2010 - R. Yantosca - Write advected tracer info to stdout
!v0.1!  03 May 2010 - R. Yantosca - Remove references to SEASALT_MOD, DIAG_OH_MOD
!v0.1!  06 May 2010 - R. Yantosca - Now make sure that ID_SPECIES contains the
!v0.1!                              indices in CSPEC for active & inactive species
!v0.1!  08 Jul 2010 - R. Yantosca - Archive tracers, OH in DIAG_COL for printout
!v0.1!EOP
!v0.1!------------------------------------------------------------------------------
!v0.1!BOC
!v0.1!
!v0.1! !LOCAL VARIABLES:
!v0.1!
!v0.1      ! Scalars
!v0.1      INTEGER            :: C, L, N
!v0.1      INTEGER            :: IU_LOG
!v0.1      CHARACTER(LEN=4)   :: PETSTR
!v0.1      CHARACTER(LEN=255) :: LOGFILE
!v0.1
!v0.1      ! Arrays
!v0.1      CHARACTER(LEN=14)  :: DRYDEP_NAME(MAX_TRACERS)
!v0.1
!v0.1      !=================================================================
!v0.1      ! Initialization
!v0.1      !=================================================================
!v0.1
!v0.1      ! Allocate pointer fields
!v0.1      ALLOCATE( COEF%SPEC_COEF( N_TRACERS, N_MEMBERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%SPEC_COEF = 0d0
!v0.1
!v0.1      ALLOCATE( COEF%SPEC_ID( N_TRACERS, N_MEMBERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%SPEC_ID = 0
!v0.1
!v0.1      ALLOCATE( COEF%SPEC_EMITTED( N_TRACERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%SPEC_EMITTED = 0
!v0.1
!v0.1      ALLOCATE( COEF%SPEC_PER_TRAC( N_TRACERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%SPEC_PER_TRAC = 0
!v0.1
!v0.1      ALLOCATE( COEF%TRAC_COEF( N_TRACERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%TRAC_COEF = 0d0
!v0.1
!v0.1      ALLOCATE( COEF%MOLWT_KG( N_TRACERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%MOLWT_KG = 0d0
!v0.1
!v0.1      ALLOCATE( COEF%XNUMOL( N_TRACERS ), STAT=RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      COEF%XNUMOL = 0d0
!v0.1
!v0.1      ! Initialize drydep name
!v0.1      DRYDEP_NAME(:)          = ''
!v0.1
!v0.1      ! Put this routine name on error trace stack
!v0.1      IDENT%LEV               = IDENT%LEV + 1
!v0.1      IDENT%I_AM( IDENT%LEV ) = 'GC_CHUNK_INIT'
!v0.1
!v0.1      ! Unit for logfile redirect
!v0.1      IU_LOG                  = IDENT%STDOUT_LUN
!v0.1
!v0.1      !=================================================================
!v0.1      ! Read chemistry mechanism
!v0.1      !=================================================================
!v0.1      
!v0.1      ! Read from data file mglob.dat
!v0.1      CALL READER( IDENT, .TRUE., TS_CHEM, RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1
!v0.1      ! Read "globchem.dat" chemistry mechanism
!v0.1      CALL READCHEM( IDENT, RC )
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1
!v0.1      ! Pass # of species and reactions back as output
!v0.1      N_SPECIES   = NGAS
!v0.1      N_REACTIONS = NTRATES(NCS) 
!v0.1
!v0.1      !=================================================================
!v0.1      ! Return w/ error if the LISOPOH species is undefined
!v0.1      !=================================================================
!v0.1      IF ( OPTIONS%USE_SEC_ORG_AEROSOLS ) THEN
!v0.1         IF ( ILISOPOH == 0 ) THEN
!v0.1            RC           = SMV_FAILURE
!v0.1            IDENT%ERRMSG = 'LISOPOH must be defined for SOA!'
!v0.1            RETURN
!v0.1         ENDIF
!v0.1      ENDIF
!v0.1
!v0.1      !=================================================================
!v0.1      ! Set ID flags for advected tracers
!v0.1      !=================================================================
!v0.1      DO N = 1, N_TRACERS
!v0.1
!v0.1         SELECT CASE( TRIM( TRACER_NAME(N) ) )
!v0.1            CASE( 'NOx', 'NOX' )
!v0.1               ID_TRACERS%NOx  = N
!v0.1            CASE( 'Ox',  'OX'  )
!v0.1               ID_TRACERS%Ox   = N
!v0.1            CASE( 'PAN'  )
!v0.1               ID_TRACERS%PAN  = N
!v0.1            CASE( 'CO'   )
!v0.1               ID_TRACERS%CO   = N
!v0.1            CASE( 'ALK4' )
!v0.1               ID_TRACERS%ALK4 = N
!v0.1            CASE( 'ISOP' )
!v0.1               ID_TRACERS%ISOP = N
!v0.1            CASE( 'HNO3' )
!v0.1               ID_TRACERS%HNO3 = N
!v0.1            CASE( 'H2O2' )
!v0.1               ID_TRACERS%H2O2 = N
!v0.1            CASE( 'ACET' )
!v0.1               ID_TRACERS%ACET = N
!v0.1            CASE( 'MEK'  )
!v0.1               ID_TRACERS%MEK  = N
!v0.1            CASE( 'ALD2' )
!v0.1               ID_TRACERS%ALD2 = N
!v0.1            CASE( 'RCHO' )
!v0.1               ID_TRACERS%RCHO = N
!v0.1            CASE( 'MVK'  )
!v0.1               ID_TRACERS%MVK  = N
!v0.1            CASE( 'MACR' )
!v0.1               ID_TRACERS%MACR = N
!v0.1            CASE( 'PMN'  )
!v0.1               ID_TRACERS%PMN  = N
!v0.1            CASE( 'PPN'  )
!v0.1               ID_TRACERS%PPN  = N
!v0.1            CASE( 'R4N2' )
!v0.1               ID_TRACERS%R4N2 = N
!v0.1            CASE( 'PRPE' )
!v0.1               ID_TRACERS%PRPE = N
!v0.1            CASE( 'C3H8' )
!v0.1               ID_TRACERS%C3H8 = N
!v0.1            CASE( 'CH2O' )
!v0.1               ID_TRACERS%CH2O = N
!v0.1            CASE( 'C2H6' )
!v0.1               ID_TRACERS%C2H6 = N
!v0.1            CASE( 'N2O5' )
!v0.1               ID_TRACERS%N2O5 = N
!v0.1            CASE( 'HNO4' )
!v0.1               ID_TRACERS%HNO4 = N
!v0.1            CASE( 'MP'   )
!v0.1               ID_TRACERS%MP   = N
!v0.1            CASE( 'DMS'  )
!v0.1               ID_TRACERS%DMS  = N
!v0.1            CASE( 'SO2'  )
!v0.1               ID_TRACERS%SO2  = N
!v0.1            CASE( 'SO4'  )
!v0.1               ID_TRACERS%SO4  = N
!v0.1            CASE( 'SO4s' )
!v0.1               ID_TRACERS%SO4s = N
!v0.1            CASE( 'MSA'  )
!v0.1               ID_TRACERS%MSA  = N
!v0.1            CASE( 'NH3'  )
!v0.1               ID_TRACERS%NH3  = N
!v0.1            CASE( 'NH4'  )
!v0.1               ID_TRACERS%NH4  = N
!v0.1            CASE( 'NIT'  )
!v0.1               ID_TRACERS%NIT  = N
!v0.1            CASE( 'NITs' )
!v0.1               ID_TRACERS%NITs = N
!v0.1            CASE( 'BCPI' )
!v0.1               ID_TRACERS%BCPI = N
!v0.1            CASE( 'OCPI' )
!v0.1               ID_TRACERS%OCPI = N
!v0.1            CASE( 'BCPO' )
!v0.1               ID_TRACERS%BCPO = N
!v0.1            CASE( 'OCPO' )
!v0.1               ID_TRACERS%OCPO = N
!v0.1            CASE( 'ALPH' )
!v0.1               ID_TRACERS%ALPH = N
!v0.1            CASE( 'LIMO' )
!v0.1               ID_TRACERS%LIMO = N
!v0.1            CASE( 'ALCO' )
!v0.1               ID_TRACERS%ALCO = N
!v0.1            CASE( 'SOG1' )
!v0.1               ID_TRACERS%SOG1 = N
!v0.1            CASE( 'SOG2' )
!v0.1               ID_TRACERS%SOG2 = N
!v0.1            CASE( 'SOG3' )
!v0.1               ID_TRACERS%SOG3 = N
!v0.1            CASE( 'SOG4' )
!v0.1               ID_TRACERS%SOG4 = N
!v0.1            CASE( 'SOA1' )
!v0.1               ID_TRACERS%SOA1 = N
!v0.1            CASE( 'SOA2' )
!v0.1               ID_TRACERS%SOA2 = N
!v0.1            CASE( 'SOA3' )
!v0.1               ID_TRACERS%SOA3 = N
!v0.1            CASE( 'SOA4' )
!v0.1               ID_TRACERS%SOA4 = N
!v0.1            CASE( 'DST1' )
!v0.1               ID_TRACERS%DST1 = N
!v0.1            CASE( 'DST2' )
!v0.1               ID_TRACERS%DST2 = N
!v0.1            CASE( 'DST3' )
!v0.1               ID_TRACERS%DST3 = N
!v0.1            CASE( 'DST4' )
!v0.1               ID_TRACERS%DST4 = N
!v0.1            CASE( 'SALA' )
!v0.1               ID_TRACERS%SALA = N
!v0.1            CASE( 'SALC' )
!v0.1               ID_TRACERS%SALC = N
!v0.1            CASE DEFAULT
!v0.1               ! Nothing
!v0.1         END SELECT
!v0.1      ENDDO
!v0.1
!v0.1      !=================================================================
!v0.1      ! Set ID flags for chemistry species (active & inactive)
!v0.1      !=================================================================
!v0.1      DO N = 1, NTSPEC(NCS)
!v0.1
!v0.1         SELECT CASE( TRIM( NAMEGAS(N) ) )
!v0.1            CASE( 'A3O2'      )  
!v0.1               ID_SPECIES%A3O2     = N
!v0.1            CASE( 'ACET'      )    
!v0.1               ID_SPECIES%ACET     = N
!v0.1            CASE( 'ACTA'      )    
!v0.1               ID_SPECIES%ACTA     = N
!v0.1            CASE( 'ALD2'      )    
!v0.1               ID_SPECIES%ALD2     = N
!v0.1            CASE( 'ALK4'      )    
!v0.1               ID_SPECIES%ALK4     = N
!v0.1            CASE( 'ATO2'      )    
!v0.1               ID_SPECIES%ATO2     = N
!v0.1            CASE( 'B3O2'      )    
!v0.1               ID_SPECIES%B3O2     = N
!v0.1            CASE( 'C2H6'      )    
!v0.1               ID_SPECIES%C2H6     = N
!v0.1            CASE( 'C3H8'      )    
!v0.1               ID_SPECIES%C3H8     = N
!v0.1            CASE( 'CH2O'      )    
!v0.1               ID_SPECIES%CH2O     = N
!v0.1            CASE( 'CH4'       )    
!v0.1               ID_SPECIES%CH4      = N
!v0.1            CASE( 'CO'        )    
!v0.1               ID_SPECIES%CO       = N
!v0.1            CASE( 'CO2'       )    
!v0.1               ID_SPECIES%CO2      = N
!v0.1            CASE( 'DMS'       )    
!v0.1               ID_SPECIES%DMS      = N
!v0.1            CASE( 'DRYCH2O'   )    
!v0.1               ID_SPECIES%DRYCH2O  = N
!v0.1            CASE( 'DRYDEP'    )    
!v0.1               ID_SPECIES%DRYDEP   = N
!v0.1            CASE( 'DRYH2O2'   )    
!v0.1               ID_SPECIES%DRYH2O2  = N
!v0.1            CASE( 'DRYHNO3'   )    
!v0.1               ID_SPECIES%DRYHNO3  = N
!v0.1            CASE( 'DRYN2O5'   )    
!v0.1               ID_SPECIES%DRYN2O5  = N
!v0.1            CASE( 'DRYNO2'    )    
!v0.1               ID_SPECIES%DRYNO2   = N
!v0.1            CASE( 'DRYO3'     )    
!v0.1               ID_SPECIES%DRYO3    = N
!v0.1            CASE( 'DRYPAN'    )    
!v0.1               ID_SPECIES%DRYPAN   = N
!v0.1            CASE( 'DRYPMN'    )    
!v0.1               ID_SPECIES%DRYPMN   = N
!v0.1            CASE( 'DRYPPN'    )    
!v0.1               ID_SPECIES%DRYPPN   = N
!v0.1            CASE( 'DRYR4N2'   )
!v0.1               ID_SPECIES%DRYR4N2  = N
!v0.1            CASE( 'EMISSION'  )
!v0.1               ID_SPECIES%EMISSION = N
!v0.1            CASE( 'EOH'       )
!v0.1               ID_SPECIES%EOH      = N
!v0.1            CASE( 'ETO2'      )
!v0.1               ID_SPECIES%ETO2     = N
!v0.1            CASE( 'ETP'       )
!v0.1               ID_SPECIES%ETP      = N
!v0.1            CASE( 'GCO3'      )
!v0.1               ID_SPECIES%GCO3     = N
!v0.1            CASE( 'GLCO3'     )
!v0.1               ID_SPECIES%GLCO3    = N
!v0.1            CASE( 'GLP'       )
!v0.1               ID_SPECIES%GLP      = N
!v0.1            CASE( 'GLPAN'     )
!v0.1               ID_SPECIES%GLPAN    = N
!v0.1            CASE( 'GLYC'      )
!v0.1               ID_SPECIES%GLYC     = N
!v0.1            CASE( 'GLYX'      )
!v0.1               ID_SPECIES%GLYX     = N
!v0.1            CASE( 'GP'        )
!v0.1               ID_SPECIES%GP       = N
!v0.1            CASE( 'GPAN'      )
!v0.1               ID_SPECIES%GPAN     = N
!v0.1            CASE( 'H'         )
!v0.1               ID_SPECIES%H        = N
!v0.1            CASE( 'H2'        )
!v0.1               ID_SPECIES%H2       = N
!v0.1            CASE( 'H2O'       )
!v0.1               ID_SPECIES%H2O      = N
!v0.1            CASE( 'H2O2'      )
!v0.1               ID_SPECIES%H2O2     = N
!v0.1            CASE( 'HAC'       )
!v0.1               ID_SPECIES%HAC      = N
!v0.1            CASE( 'HCOOH'     )
!v0.1               ID_SPECIES%HCOOH    = N
!v0.1            CASE( 'HNO2'      )
!v0.1               ID_SPECIES%HNO2     = N
!v0.1            CASE( 'HNO3'      )
!v0.1               ID_SPECIES%HNO3     = N
!v0.1            CASE( 'HNO4'      )
!v0.1               ID_SPECIES%HNO4     = N
!v0.1            CASE( 'HO2'       )
!v0.1               ID_SPECIES%HO2      = N
!v0.1            CASE( 'IALD'      )
!v0.1               ID_SPECIES%IALD     = N
!v0.1            CASE( 'IAO2'      )
!v0.1               ID_SPECIES%IAO2     = N
!v0.1            CASE( 'IAP'       )
!v0.1               ID_SPECIES%IAP      = N
!v0.1            CASE( 'INO2'      )
!v0.1               ID_SPECIES%INO2     = N
!v0.1            CASE( 'INPN'      )
!v0.1               ID_SPECIES%INPN     = N
!v0.1            CASE( 'ISN1'      )
!v0.1               ID_SPECIES%ISN1     = N
!v0.1            CASE( 'ISNO3'     )
!v0.1               ID_SPECIES%ISNO3    = N
!v0.1            CASE( 'ISNP'      )
!v0.1               ID_SPECIES%ISNP     = N
!v0.1            CASE( 'ISOP'      )
!v0.1               ID_SPECIES%ISOP     = N
!v0.1            CASE( 'KO2'       )
!v0.1               ID_SPECIES%KO2      = N
!v0.1            CASE( 'LISOPOH'   )
!v0.1               ID_SPECIES%LISOPOH  = N
!v0.1            CASE( 'M'         )
!v0.1               ID_SPECIES%M        = N
!v0.1            CASE( 'MACR'      )
!v0.1               ID_SPECIES%MACR     = N
!v0.1            CASE( 'MAN2'      )
!v0.1               ID_SPECIES%MAN2     = N
!v0.1            CASE( 'MAO3'      )
!v0.1               ID_SPECIES%MAO3     = N
!v0.1            CASE( 'MAOP'      )
!v0.1               ID_SPECIES%MAOP     = N
!v0.1            CASE( 'MAP'       )
!v0.1               ID_SPECIES%MAP      = N
!v0.1            CASE( 'MCO3'      )
!v0.1               ID_SPECIES%MCO3     = N
!v0.1            CASE( 'MEK'       )
!v0.1               ID_SPECIES%MEK      = N
!v0.1            CASE( 'MGLY'      )
!v0.1               ID_SPECIES%MGLY     = N
!v0.1            CASE( 'MNO3'      )
!v0.1               ID_SPECIES%MNO3     = N
!v0.1            CASE( 'MO2'       )
!v0.1               ID_SPECIES%MO2      = N
!v0.1            CASE( 'MOH'       )
!v0.1               ID_SPECIES%MOH      = N
!v0.1            CASE( 'MP'        )
!v0.1               ID_SPECIES%MP       = N
!v0.1            CASE( 'MRO2'      )
!v0.1               ID_SPECIES%MRO2     = N
!v0.1            CASE( 'MRP'       )
!v0.1               ID_SPECIES%MRP      = N
!v0.1            CASE( 'MSA'       )
!v0.1               ID_SPECIES%MSA      = N
!v0.1            CASE( 'MVK'       )
!v0.1               ID_SPECIES%MVK      = N
!v0.1            CASE( 'MVN2'      )
!v0.1               ID_SPECIES%MVN2     = N
!v0.1            CASE( 'N2'        )
!v0.1               ID_SPECIES%N2       = N
!v0.1            CASE( 'N2O'       )
!v0.1               ID_SPECIES%N2O      = N
!v0.1            CASE( 'N2O5'      )
!v0.1               ID_SPECIES%N2O5     = N
!v0.1            CASE( 'NH2'       )
!v0.1               ID_SPECIES%NH2      = N
!v0.1            CASE( 'NH3'       )
!v0.1               ID_SPECIES%NH3      = N
!v0.1            CASE( 'NO'        )
!v0.1               ID_SPECIES%NO       = N
!v0.1            CASE( 'NO2'       )
!v0.1               ID_SPECIES%NO2      = N
!v0.1            CASE( 'NO3'       )
!v0.1               ID_SPECIES%NO3      = N
!v0.1            CASE( 'O'         )
!v0.1               ID_SPECIES%O        = N
!v0.1            CASE( 'O1D'       )
!v0.1               ID_SPECIES%O1D      = N
!v0.1            CASE( 'O2'        )
!v0.1               ID_SPECIES%O2       = N
!v0.1            CASE( 'O2CH2OH'   )
!v0.1               ID_SPECIES%O2CH2OH  = N
!v0.1            CASE( 'O3'        )
!v0.1               ID_SPECIES%O3       = N
!v0.1            CASE( 'OH'        )
!v0.1               ID_SPECIES%OH       = N
!v0.1            CASE( 'PAN'       )
!v0.1               ID_SPECIES%PAN      = N
!v0.1            CASE( 'PMN'       )
!v0.1               ID_SPECIES%PMN      = N
!v0.1            CASE( 'PO2'       )
!v0.1               ID_SPECIES%PO2      = N
!v0.1            CASE( 'PP'        )
!v0.1               ID_SPECIES%PP       = N
!v0.1            CASE( 'PPN'       )
!v0.1               ID_SPECIES%PPN      = N
!v0.1            CASE( 'PRN1'      )
!v0.1               ID_SPECIES%PRN1     = N
!v0.1            CASE( 'PRPE'      )
!v0.1               ID_SPECIES%PRPE     = N
!v0.1            CASE( 'PRPN'      )
!v0.1               ID_SPECIES%PRPN     = N
!v0.1            CASE( 'R4N1'      )
!v0.1               ID_SPECIES%R4N1     = N
!v0.1            CASE( 'R4N2'      )
!v0.1               ID_SPECIES%R4N2     = N
!v0.1            CASE( 'R4O2'      )
!v0.1               ID_SPECIES%R4O2     = N
!v0.1            CASE( 'R4P'       )
!v0.1               ID_SPECIES%R4P      = N
!v0.1            CASE( 'RA3P'      )
!v0.1               ID_SPECIES%RA3P     = N
!v0.1            CASE( 'RB3P'      )
!v0.1               ID_SPECIES%RB3P     = N
!v0.1            CASE( 'RCHO'      )
!v0.1               ID_SPECIES%RCHO     = N
!v0.1            CASE( 'RCO3'      )
!v0.1               ID_SPECIES%RCO3     = N
!v0.1            CASE( 'RCOOH'     )
!v0.1               ID_SPECIES%RCOOH    = N
!v0.1            CASE( 'RIO1'      )
!v0.1               ID_SPECIES%RIO1     = N
!v0.1            CASE( 'RIO2'      )
!v0.1               ID_SPECIES%RIO2     = N
!v0.1            CASE( 'RIP'       )
!v0.1               ID_SPECIES%RIP      = N
!v0.1            CASE( 'ROH'       )
!v0.1               ID_SPECIES%ROH      = N
!v0.1            CASE( 'RP'        )
!v0.1               ID_SPECIES%RP       = N
!v0.1            CASE( 'SO2'       )
!v0.1               ID_SPECIES%SO2      = N
!v0.1            CASE( 'SO4'       )
!v0.1               ID_SPECIES%SO4      = N
!v0.1            CASE( 'VRO2'      )
!v0.1               ID_SPECIES%VRO2     = N
!v0.1            CASE( 'VRP'       )
!v0.1               ID_SPECIES%VRP      = N
!v0.1            CASE DEFAULT
!v0.1               ! Nothing
!v0.1         END SELECT
!v0.1      ENDDO
!v0.1
!v0.1      !=================================================================
!v0.1      ! Save tracer number in the ID_TRACERS object
!v0.1      ! Create the coefficient arrays that link tracers & species
!v0.1      !=================================================================
!v0.1
!v0.1      ! Write header text
!v0.1      WRITE( IU_LOG, '(/,a)' ) REPEAT( '=', 79 )
!v0.1      WRITE( IU_LOG, 200     )
!v0.1      WRITE( IU_LOG, '(a,/)' ) REPEAT( '=', 79 )
!v0.1      WRITE( IU_LOG, 210     ) 
!v0.1      WRITE( IU_LOG, '(  a)' ) REPEAT( '-', 30 )
!v0.1
!v0.1      ! Formats
!v0.1 200  FORMAT(
!v0.1     &  'ADVECTED TRACERS (==> denotes emitted constituent species)' )
!v0.1 210  FORMAT( '  # Tracer          g/mole' )
!v0.1
!v0.1      ! Loop over all advected tracers
!v0.1      DO N = 1, N_TRACERS
!v0.1         
!v0.1         IF ( N == ID_TRACERS%NOx ) THEN
!v0.1            COEF%SPEC_COEF    (N,1:4) = (/ 1d0, 1d0, 1d0, 1d0 /)
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%NO2
!v0.1            COEF%SPEC_ID      (N,2  ) = ID_SPECIES%NO
!v0.1            COEF%SPEC_ID      (N,3  ) = ID_SPECIES%NO3
!v0.1            COEF%SPEC_ID      (N,4  ) = ID_SPECIES%HNO2
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 4
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 2     ! NO is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 46d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%Ox ) THEN
!v0.1            COEF%SPEC_COEF    (N,1:3) = (/ 1d0, 1d0, 2d0 /)
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%O3
!v0.1            COEF%SPEC_ID      (N,2  ) = ID_SPECIES%NO2
!v0.1            COEF%SPEC_ID      (N,3  ) = ID_SPECIES%NO3  
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 3
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! O3 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 48d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%PAN ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PAN
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 121d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%CO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%CO
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1    ! CO is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 28d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ALK4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 4d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ALK4
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 4d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! ALK4 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ISOP ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 5d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ISOP
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 5d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! ISOP is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%HNO3 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%HNO3
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! HNO3 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 63d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%H2O2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%H2O2
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 34d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ACET ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 3d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ACET
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 3d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! ACET is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%MEK ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 4d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MEK
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 4d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! MEK is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ALD2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 2d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ALD2
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 2d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! ALD2 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%RCHO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%RCHO
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 58d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%MVK ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MVK
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 70d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%MACR ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MACR
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 70d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%PMN ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PMN
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 147d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%PPN ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PPN
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 135d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%R4N2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%R4N2
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 119d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%PRPE ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 3d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PRPE
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 3d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! PRPE is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%C3H8 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 3d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%C3H8
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 3d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! C3H8 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%CH2O ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%CH2O
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! CH2O is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 30d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%C2H6 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 2d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%C2H6
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 2d0
!v0.1            COEF%SPEC_EMITTED (N    ) = 1     ! C2H6 is emitted
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%N2O5 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%N2O5
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 105d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%HNO4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%HNO4
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 79d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%MP ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MP
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 48d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%DMS ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%DMS
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 62d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SO2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%SO2
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 64d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SO4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%SO4
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 96d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SO4s ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 96d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%MSA ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MSA
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 96d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%NH3 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 17d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%NH4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 18d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%NIT ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 62d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%NITs ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 62d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%BCPI ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%OCPI ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%BCPO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%OCPO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 12d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ALPH ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 136.23d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%LIMO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 136.23d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%ALCO ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 142d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOG1 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 150d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOG2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 160d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOG3 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 220d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOG4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 130d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOA1 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 150d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOA2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 160d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOA3 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 220d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SOA4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 130d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%DST1 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 29d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%DST2 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 29d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%DST3 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 29d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%DST4 ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 29d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SALA ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 36d-3
!v0.1
!v0.1         ELSE IF ( N == ID_TRACERS%SALC ) THEN
!v0.1            COEF%SPEC_COEF    (N,1  ) = 1d0
!v0.1            COEF%SPEC_ID      (N,1  ) = 0
!v0.1            COEF%SPEC_PER_TRAC(N    ) = 1
!v0.1            COEF%TRAC_COEF    (N    ) = 1d0
!v0.1            COEF%MOLWT_KG     (N    ) = 36d-3
!v0.1
!v0.1         ELSE
!v0.1
!v0.1            ! Invalid tracer -- return with error status
!v0.1            WRITE( IDENT%ERRMSG, 300 ) N
!v0.1            RC = SMV_FAILURE
!v0.1            RETURN
!v0.1            
!v0.1         ENDIF         
!v0.1
!v0.1         ! XNUMOL = AVO / MOLWT = ratio of molec/kg for each tracer
!v0.1         COEF%XNUMOL(N) = AVO / COEF%MOLWT_KG(N)
!v0.1
!v0.1         ! Write tracer number, name, & mol wt
!v0.1         WRITE( IU_LOG, 310 ) N, TRACER_NAME(N), 
!v0.1     &                            COEF%MOLWT_KG(N) * 1000d0
!v0.1
!v0.1         ! If the tracer ishas constituent species ...
!v0.1         IF ( COEF%SPEC_PER_TRAC(N) > 1  .or. 
!v0.1     &        COEF%SPEC_EMITTED(N)  > 0 ) THEN
!v0.1         
!v0.1            ! ... print out constituent species names
!v0.1            DO C = 1, COEF%SPEC_PER_TRAC(N) 
!v0.1
!v0.1               ! Highlight emitted species with ==> arrow
!v0.1               IF ( COEF%SPEC_EMITTED(N) == C ) THEN
!v0.1                  WRITE( IU_LOG, 320 ) COEF%SPEC_COEF(N,C), 
!v0.1     &                                  NAMEGAS( COEF%SPEC_ID(N,C) )
!v0.1               ELSE
!v0.1                  WRITE( IU_LOG, 330 ) COEF%SPEC_COEF(N,C),
!v0.1     &                                  NAMEGAS( COEF%SPEC_ID(N,C) )
!v0.1               ENDIF
!v0.1
!v0.1            ENDDO 
!v0.1         ENDIF
!v0.1
!v0.1      ENDDO
!v0.1
!v0.1      ! Formats
!v0.1 300  FORMAT( 'Tracer ', i4, ' is not defined!' )
!v0.1 310  FORMAT( I3, 1x, A10, 6x, F6.1 )
!v0.1 320  FORMAT( 5x, '===> ', f4.1, 1x, A4  )
!v0.1 330  FORMAT( 5x, '---> ', f4.1, 1x, A4  )
!v0.1
!v0.1      !=================================================================
!v0.1      ! Initialize dry deposition
!v0.1      !=================================================================
!v0.1      IF ( OPTIONS%USE_DRYDEP ) THEN
!v0.1
!v0.1         ! Call init method from "drydep_mod.f"
!v0.1         CALL INIT_DRYDEP( IDENT, ID_TRACERS, ID_DRYDEP,   
!v0.1     &                     COEF,  N_DRYDEP,   DRYDEP_NAME, RC )
!v0.1
!v0.1         ! Return w/ error if necessary
!v0.1         IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      ENDIF
!v0.1
!v0.1      !=================================================================
!v0.1      ! Initialize wet deposition
!v0.1      !
!v0.1      ! NOTE: Turning off wet deposition will also turn off the 
!v0.1      !       wet scavenging of tracer in cloud updrafts. 
!v0.1      !=================================================================
!v0.1      IF ( OPTIONS%USE_WETDEP .or. OPTIONS%USE_CONVECTION ) THEN
!v0.1
!v0.1         ! Call init method from "wetdep_mod.f"
!v0.1         CALL INIT_WETDEP( IDENT,     L_COLUMN,  ID_TRACERS, 
!v0.1     &                     ID_WETDEP, N_WETDEP,  RC )
!v0.1
!v0.1         ! Return w/ error if necessary
!v0.1         IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1      ENDIF
!v0.1
!v0.1      !=================================================================
!v0.1      ! Flag emission & dry deposition reactions w/in SMVGEAR
!v0.1      !=================================================================
!v0.1      CALL SETEMDEP( IDENT,       N_TRACERS,   N_DRYDEP, 
!v0.1     &               TRACER_NAME, DRYDEP_NAME, COEF,     RC )
!v0.1
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1
!v0.1      !=================================================================
!v0.1      ! Initialize the FAST-J photolysis mechanism
!v0.1      !=================================================================
!v0.1      CALL INPHOT( IDENT, L_COLUMN, NPHOT, N_JV, RC )
!v0.1
!v0.1      IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1
!v0.1      !=================================================================
!v0.1      ! Initialize the DIAG_COL object for diagnostic output
!v0.1      !=================================================================
!v0.1
!v0.1      ! Set a flag if we should print the diagnostic output to a file
!v0.1      DIAG_COL%DO_PRINT = OPTIONS%USE_DEBUG_PRINT
!v0.1
!v0.1      ! Initialize the other fields of DIAG_COL
!v0.1      IF ( DIAG_COL%DO_PRINT .and. IDENT%PET == 0  ) THEN
!v0.1         
!v0.1         ! Allocate the TRACER fields
!v0.1         IF ( .not. ASSOCIATED( DIAG_COL%TRACER ) ) THEN
!v0.1            ALLOCATE( DIAG_COL%TRACER(MAX_CHUNK,N_TRACERS+1), STAT=RC )
!v0.1            IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1         ENDIF
!v0.1         
!v0.1         ! Allocate the tracer name fields
!v0.1         IF ( .not. ASSOCIATED( DIAG_COL%NAME ) ) THEN
!v0.1            ALLOCATE( DIAG_COL%NAME( N_TRACERS+1 ), STAT=RC )
!v0.1            IF ( RC /= SMV_SUCCESS ) RETURN
!v0.1         ENDIF
!v0.1
!v0.1         ! Init fields
!v0.1         DIAG_COL%N_DIAG   = N_TRACERS + 1
!v0.1         DIAG_COL%TRACER   = 0d0
!v0.1         DIAG_COL%COUNT    = 0
!v0.1         DIAG_COL%FILENAME = 'diag_col.txt'
!v0.1         DIAG_COL%LUN      = 750
!v0.1
!v0.1         ! Initialize tracer names (OH is NTRACERS+1)
!v0.1         DO N = 1, N_TRACERS
!v0.1            DIAG_COL%NAME(N) = TRIM( TRACER_NAME(N) )
!v0.1         ENDDO
!v0.1         DIAG_COL%NAME(N_TRACERS+1) = 'OH'
!v0.1      ENDIF
!v0.1
!v0.1      !================================================================
!v0.1      ! Succesful return!
!v0.1      !================================================================
!v0.1
!v0.1      ! Set error code to success
!v0.1      RC                      = SMV_SUCCESS
!v0.1
!v0.1      ! Remove this routine name from error trace stack
!v0.1      IDENT%I_AM( IDENT%LEV ) = ''
!v0.1      IDENT%LEV               = IDENT%LEV - 1
!v0.1
!v0.1      END SUBROUTINE GC_CHUNK_INIT
!v0.1
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_chunk_final
!
! !DESCRIPTION: Subroutine GC\_CHUNK\_FINAL deallocates pointers and
!  arrays used in the chemistry. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GC_CHUNK_FINAL( IDENT, COEF, RC )
!
! !USES:
!
#     include "smv_errcode.h"           ! Error codes
!
! !INPUT/OUTPUT PARAMETERS:
!     
      ! Object w/ info from ESMF and traceback stack
      TYPE(GC_IDENT),    INTENT(INOUT) :: IDENT

      ! Object w/ info for mapping species <--> tracers
      TYPE(SPEC_2_TRAC), INTENT(INOUT) :: COEF      
!
! !OUTPUT PARAMETERS:
!
      ! Return code
      INTEGER,           INTENT(OUT)   :: RC
!
! !REVISION HISTORY: 
!  30 Apr 2009 - R. Yantosca - Initial version
!  05 May 2009 - P. Le Sager - now use module variables; remove call to
!                              cleanup_dust
!  05 Jun 2009 - R. Yantosca - Now deallocate COEF%MOLEC_KG & COEF%XNUMOL
!  30 Jun 2009 - R. Yantosca - Moved here from "chemistry_mod.f" 
!  30 Apr 2010 - R. Yantosca - Now pass IDENT via the arg list
!  30 Apr 2010 - R. Yantosca - Now call CLEANUP_SCHEM outside this routine
!  03 Jun 2010 - R. Yantosca - Removed calls to CLEANUP_* routines.  These
!                              referenced 3-D arrays that only need to be
!                              used in GEOS-Chem.
!  08 Jul 2010 - R. Yantosca - Archive tracers, OH in DIAG_COL for printout
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: IU_LOG

      !=================================================================
      ! Initialization
      !=================================================================

      ! Put this routine name on error trace stack
      IDENT%LEV               = IDENT%LEV + 1
      IDENT%I_AM( IDENT%LEV ) = 'GC_CHUNK_FINAL'

      ! Unit for logfile redirect
      IU_LOG                  = IDENT%STDOUT_LUN

      !=================================================================
      ! Deallocate pointer fields of the COEF object
      !=================================================================
      
      IF ( ASSOCIATED( COEF%SPEC_COEF ) ) THEN
         DEALLOCATE( COEF%SPEC_COEF )
      ENDIF

      IF ( ASSOCIATED( COEF%SPEC_ID ) ) THEN
         DEALLOCATE( COEF%SPEC_ID )
      ENDIF

      IF ( ASSOCIATED( COEF%SPEC_EMITTED ) ) THEN
         DEALLOCATE( COEF%SPEC_EMITTED )
      ENDIF

      IF ( ASSOCIATED( COEF%SPEC_PER_TRAC ) ) THEN
         DEALLOCATE( COEF%SPEC_PER_TRAC )
      ENDIF

      IF ( ASSOCIATED( COEF%TRAC_COEF ) ) THEN
         DEALLOCATE( COEF%TRAC_COEF )
      ENDIF

      IF ( ASSOCIATED( COEF%MOLWT_KG ) ) THEN
         DEALLOCATE( COEF%MOLWT_KG )
      ENDIF

      IF ( ASSOCIATED( COEF%XNUMOL ) ) THEN
         DEALLOCATE( COEF%XNUMOL )
      ENDIF

      !=================================================================
      ! Print out diagnostic output and finalize the DIAG_COL object
      !=================================================================
      IF ( DIAG_COL%DO_PRINT .and. IDENT%PET == 0 ) THEN

         ! Write to file
         CALL PRINT_DIAG_COL( IDENT, RC )
         IF ( RC /= SMV_SUCCESS ) RETURN

         ! Deallocate tracer field array
         IF ( ASSOCIATED( DIAG_COL%TRACER ) ) THEN
            DEALLOCATE( DIAG_COL%TRACER )
         ENDIF

         ! Deallocate tracer names array
         IF ( ASSOCIATED( DIAG_COL%NAME ) ) THEN
            DEALLOCATE( DIAG_COL%NAME )
         ENDIF

         ! Blank other fields
         DIAG_COL%DO_PRINT = .FALSE.
         DIAG_COL%COUNT    = 0
         DIAG_COL%FILENAME = ''
         DIAG_COL%LUN      = 0
         DIAG_COL%N_DIAG   = 0
      ENDIF

      !================================================================
      ! Succesful return!
      !================================================================

      ! Set error code to success
      RC                      = SMV_SUCCESS

      ! Remove this routine name from error trace stack
      IDENT%I_AM( IDENT%LEV ) = ''
      IDENT%LEV               = IDENT%LEV - 1

      END SUBROUTINE GC_CHUNK_FINAL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time
!
! !DESCRIPTION: Function ITS\_TIME returns TRUE if the elapsed time is a 
!  multiple of the timestep, or FALSE otherwise.  This is used to check
!  if it is time to do emissions or chemistry, for example.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME( TIME_ELAPSED, TIMESTEP ) RESULT( IT_IS_TIME )
!
! !INPUT PARAMETERS:
!
      REAL*8,  INTENT(IN) :: TIME_ELAPSED
      REAL*8,  INTENT(IN) :: TIMESTEP
!
! !RETURN_CODE
!
      ! Return code
      LOGICAL             :: IT_IS_TIME
!
! !REVISION HISTORY: 
!     30 Jun 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      IT_IS_TIME = ( MOD( TIME_ELAPSED, TIMESTEP ) == 0d0 )

      END FUNCTION ITS_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_units
!
! !DESCRIPTION: Subroutine CONVERT\_UNITS converts the tracer concentration 
!  array from [kg] to [mol/mol] mixing ratio, or vice versa.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CONVERT_UNITS( IDENT, DIMINFO, IFLAG, 
     &                          COEF,  AD,      TRACER, RC ) 
!
! !USES:
!
      USE ERROR_MOD, ONLY : IS_SAFE_DIV   ! From "GeosUtil" library

#     include "smv_errcode.h"             ! Error codes
#     include "smv_physconst.h"           ! Physical constants
!
! !INPUT PARAMETERS:
!
      ! Object w/ dimension information
      TYPE(GC_DIMS),     INTENT(IN)    :: DIMINFO

      ! If IFLAG==1, then convert from [kg     ] --> [mol/mol]
      ! If IFLAG==2, then convert from [mol/mol] --> [kg     ]
      INTEGER,           INTENT(IN)    :: IFLAG

      ! Object containing molecular weight info
      TYPE(SPEC_2_TRAC), INTENT(IN)    :: COEF

      ! Air mass in grid box [kg]
      REAL*8,            INTENT(IN)    :: AD(:)
!
! !INPUT/OUTPUT PARAMETERS:
! 
      ! Object with info from ESMF and traceback info
      TYPE(GC_IDENT),    INTENT(INOUT) :: IDENT

      ! Tracer concentration in [kg] or [mol/mol]
      REAL*8,            INTENT(INOUT) :: TRACER(:,:)
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,           INTENT(OUT)   :: RC
!
! !REMARKS:
!  The conversion from [kg] to [mol/mol] mixing ratio is as follows:
!                                                                             .
!     kg tracer(N)       1        Air mol wt       moles tracer
!     -----------  * -------- *  -------------  =  ------------
!          1          kg air     tracer mol wt      moles air  
!                                                                             .
!  Therefore, with:
!                                                                             .
!     TCVV(N) = 28.97d-3 / molecular weight of tracer [kg]
!             = mol. wt. of air (AMU) / mol. wt. of tracer (AMU)
!                                                                             .
!     AD(L)   = mass of air (kg) in grid box L  (L = vertical index)
!                                                                             .
!  The conversion is:
!                                                                             .
!     TRACER(L,N) [kg] * TCVV(N) / AD(L) = TRACER(L,N) [mol/mol]
!                                                                             .
!  And the inverse conversion [mol/mol] to [kg] is:
!                                                                             .
!     TRACER(L,N) [mol/mol] / TCVV(N) * AD(L) = TRACER(L,N) [kg]
!
! !REVISION HISTORY: 
!  07 Jul 2009 - R. Yantosca - Initial version
!  14 Dec 2009 - R. Yantosca - Now use [mol/mol] in comments instead of [v/v]
!  14 Dec 2009 - R. Yantosca - Now get molecular weight of air from MW_AIR
!                              in header file "smv_physconst.h"
!  14 Dec 2009 - R. Yantosca - Now make sure division can be performed,
!                              otherwise return an error condition
!  04 May 2010 - R. Yantosca - Add IDENT to the argument list
!  01 Jun 2010 - R. Yantosca - Add DIMINFO to the argument list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local Variables
      INTEGER :: L, N
      REAL*8  :: TCVV

      !=================================================================
      ! Initialization
      !=================================================================

      ! Put this routine name on error trace stack
      IDENT%LEV               = IDENT%LEV + 1
      IDENT%I_AM( IDENT%LEV ) = 'CONVERT_UNITS'

      ! Case statement
      SELECT CASE ( IFLAG )

         !----------------------------------------------
         ! IFLAG = 1: Convert from [kg] -> [mol/mol] 
         !----------------------------------------------
         CASE ( 1 )

            DO N = 1, DIMINFO%N_TRACERS
            DO L = 1, DIMINFO%L_COLUMN

               ! Ratio of air/tracer mol wts
               TCVV        = MW_AIR / COEF%MOLWT_KG(N)

               ! Convert units if the division can be performed
               ! Otherwise return w/ error
               IF ( IS_SAFE_DIV( TCVV, AD(L) ) ) THEN
                  TRACER(L,N) = TRACER(L,N) * TCVV / AD(L)
               ELSE
                  WRITE( IDENT%ERRMSG, 200 )
 200              FORMAT( 'Error: Cannot do division, IFLAG=1' )
                  RC = SMV_FAILURE
                  RETURN
               ENDIF

            ENDDO
            ENDDO

         !----------------------------------------------
         ! IFLAG = 2: Convert from [mol/mol] -> [kg] 
         !---------------------------------------------- 
         CASE ( 2 )

            DO N = 1, DIMINFO%N_TRACERS
            DO L = 1, DIMINFO%L_COLUMN

               ! Ratio of air/tracer mol wts
               TCVV        = MW_AIR / COEF%MOLWT_KG(N)
               
               ! Convert units if the division can be performed
               ! Otherwise return w/ error
               IF ( IS_SAFE_DIV( AD(L), TCVV ) ) THEN
                  TRACER(L,N) = TRACER(L,N) * AD(L) / TCVV
               ELSE
                  WRITE( IDENT%ERRMSG, 210 )
 210              FORMAT( 'Error: Cannot do division, IFLAG=1' )
                  RC = SMV_FAILURE
                  RETURN
               ENDIF

            ENDDO     
            ENDDO

      END SELECT

      !================================================================
      ! Succesful return!
      !================================================================

      ! Set error code to success
      RC                      = SMV_SUCCESS

      ! Remove this routine name from error trace stack
      IDENT%I_AM( IDENT%LEV ) = ''
      IDENT%LEV               = IDENT%LEV - 1

      END SUBROUTINE CONVERT_UNITS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_diag_avg
!
! !DESCRIPTION: Function PRINT\_DIAG\_AVG computes the time-averaged tracer
!  and OH concentrations and writes them to a file for debug output.  This
!  facilitates comparing similar quantities when evaluating the performance
!  of the GEOS-Chem chunk chemistry within the ESMF/MAPL environment.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PRINT_DIAG_COL( IDENT, RC )
!
! !USES:
!
#     include "smv_errcode.h"                  ! Error codes
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(GC_IDENT), INTENT(INOUT) :: IDENT   ! Obj w/ info from ESMF etc
!
! !OUTPUT PARAMETERS:
!
      ! Return code
      INTEGER,        INTENT(OUT)   :: RC      ! Return code
!
! !REVISION HISTORY: 
!  08 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: N, LUN
      REAL*8  :: CT

      !=================================================================
      ! Initialization
      !=================================================================

      ! Put this routine name on error trace stack
      IDENT%LEV               = IDENT%LEV + 1
      IDENT%I_AM( IDENT%LEV ) = 'PRINT_DIAG_COL'

      !===================================================================
      ! Print output for validation at test box
      !===================================================================

      ! Number of timesteps to average
      CT  = DBLE( DIAG_COL%COUNT )

      ! Open file
      OPEN( DIAG_COL%LUN,     FILE=TRIM( DIAG_COL%FILENAME ), 
     &      STATUS='UNKNOWN', IOSTAT=RC ) 
      
      ! Return if error opening file
      IF ( RC /= SMV_SUCCESS ) THEN
         IDENT%ERRMSG = 'Could not open file ' // 
     &                  TRIM( DIAG_COL%FILENAME )
         RC           = SMV_FAILURE
         RETURN
      ENDIF

      ! Print number of timesteps
      WRITE( DIAG_COL%LUN, 100 ) DIAG_COL%COUNT
 100  FORMAT( 'Average concentrations after ', i5, ' timesteps', / )

      ! Make sure we have at least one timestep before printout
      IF ( DIAG_COL%COUNT <= 0 ) THEN
         IDENT%ERRMSG = 'Need at least one timestep for output!'
         RC           = SMV_FAILURE     
         RETURN
      ENDIF

      ! Only print out 1st 36 levels, because we start lumping together
      ! the stratospheric levels starting with level 37.
      DO N = 1, DIAG_COL%N_DIAG
         WRITE( DIAG_COL%LUN, '(a)' ) TRIM( DIAG_COL%NAME(N) )
         WRITE( DIAG_COL%LUN, 110   ) DIAG_COL%TRACER(1:36,N) / CT
 110     FORMAT( 4( es19.12, 1x ) )
      ENDDO

      ! Close file
      CLOSE( DIAG_COL%LUN )

      !================================================================
      ! Succesful return!
      !================================================================

      ! Set error code to success
      RC                      = SMV_SUCCESS

      ! Remove this routine name from error trace stack
      IDENT%I_AM( IDENT%LEV ) = ''
      IDENT%LEV               = IDENT%LEV - 1

      END SUBROUTINE PRINT_DIAG_COL
!EOC
      END MODULE GC_CHUNK_MOD



