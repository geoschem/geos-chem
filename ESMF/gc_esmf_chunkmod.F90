#if defined (ESMF_)
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
      MODULE GC_ESMF_CHUNKMOD
!
! !USES:
!      
      USE GC_TYPE_MOD          

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GC_CHUNK_INIT
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
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_chunk_init
!
! !DESCRIPTION: Subroutine GC\_CHUNK\_INIT calls the various 
!  initialization routines that read the setup files for the GEOS-Chem
!  chunk code.  Also, ID flags for advected tracers, chemical species,
!  dry deposition species and wet deposition species are defined.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GC_CHUNK_INIT( IDENT,       OPTIONS,    L_COLUMN,     &
     &                           N_TRACERS,   N_MEMBERS,  N_DRYDEP,    &
     &                           N_WETDEP,    TS_CHEM,    TRACER_NAME, &
     &                           ID_SPECIES,  ID_TRACERS, COEF,        &
     &                           ID_DRYDEP,   ID_WETDEP,  N_SPECIES,   &
     &                           N_REACTIONS, N_JV,       RC )   
!
! !USES:
!
!v0.1      USE DRYDEP_MOD,  ONLY : INIT_DRYDEP
!v0.1      USE WETDEP_MOD,  ONLY : INIT_WETDEP
      
      IMPLICIT NONE

#     include "smv_dimension.h"  ! Dimensions for common blocks
#     include "comode_loop.h"    ! SMVGEAR common blocks
#     include "smv_errcode.h"    ! Error codes
#     include "smv_physconst.h"  ! Physical constants
!
! !INPUT PARAMETERS:
!
      ! Object with logical flags
      TYPE(GC_OPTIONS),  INTENT(IN)    :: OPTIONS

      ! Number of boxes in the atmospheric chunk
      INTEGER,           INTENT(IN)    :: L_COLUMN

      ! Number of advected tracers
      INTEGER,           INTENT(IN)    :: N_TRACERS

      ! Max # of species per chemical family
      INTEGER,           INTENT(IN)    :: N_MEMBERS

      ! Chemistry timestep [minutes]
      REAL*8,            INTENT(IN)    :: TS_CHEM

      ! Names of advected tracers
      CHARACTER(LEN=*),  INTENT(IN)    :: TRACER_NAME(N_TRACERS)
!
! !INPUT/OUTPUT PARAMETERS:
!
      ! Object for Identification info from the Gridded Component
      TYPE(GC_IDENT),    INTENT(INOUT) :: IDENT

      ! Object for ID flags for SMVGEAR chemical species
      TYPE(ID_SPEC),     INTENT(INOUT) :: ID_SPECIES

      ! Object for ID flags for advected tracers
      TYPE(ID_TRAC),     INTENT(INOUT) :: ID_TRACERS

      ! Object for # of species per advected tracer, etc.
      TYPE(SPEC_2_TRAC), INTENT(INOUT) :: COEF

      ! Object for ID flags for dry-deposited and wet-deposited tracers
      TYPE(ID_DRYD),     INTENT(INOUT) :: ID_DRYDEP
      TYPE(ID_WETD),     INTENT(INOUT) :: ID_WETDEP
!
! !OUTPUT PARAMETERS:
!
      ! Number of dry deposited and wet-deposited species
      INTEGER,           INTENT(OUT)   :: N_DRYDEP
      INTEGER,           INTENT(OUT)   :: N_WETDEP

      ! Number of chemical species, reactions, and photolysis reactions
      INTEGER,           INTENT(OUT)   :: N_SPECIES
      INTEGER,           INTENT(OUT)   :: N_REACTIONS
      INTEGER,           INTENT(OUT)   :: N_JV

      ! Return code
      INTEGER,           INTENT(OUT)   :: RC 
!
! !REVISION HISTORY: 
!  22 Jun 2009 - R. Yantosca - Initial version
!  15 Jul 2009 - R. Yantosca - Now call init routines for drydep & wetdep
!  14 Apr 2010 - R. Yantosca - Now de-hardwire definition of ID_TRACERS
!  16 Apr 2010 - R. Yantosca - Now pass back N_SPECIES, N_REACTIONS as outputs
!  23 Apr 2010 - R. Yantosca - Added IDENT to the argument list
!  23 Apr 2010 - R. Yantosca - Redirect stdout to a log file
!  23 Apr 2010 - R. Yantosca - Write advected tracer info to stdout
!  03 May 2010 - R. Yantosca - Remove references to SEASALT_MOD, DIAG_OH_MOD
!  06 May 2010 - R. Yantosca - Now make sure that ID_SPECIES contains the
!                              indices in CSPEC for active & inactive species
!  08 Jul 2010 - R. Yantosca - Archive tracers, OH in DIAG_COL for printout
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: C, L, N
      INTEGER            :: IU_LOG
      CHARACTER(LEN=4)   :: PETSTR
      CHARACTER(LEN=255) :: LOGFILE

      ! Arrays
      CHARACTER(LEN=14)  :: DRYDEP_NAME(MAX_TRACERS)

      !=================================================================
      ! Initialization
      !=================================================================

      ! Allocate pointer fields
      ALLOCATE( COEF%SPEC_COEF( N_TRACERS, N_MEMBERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%SPEC_COEF = 0d0

      ALLOCATE( COEF%SPEC_ID( N_TRACERS, N_MEMBERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%SPEC_ID = 0

      ALLOCATE( COEF%SPEC_EMITTED( N_TRACERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%SPEC_EMITTED = 0

      ALLOCATE( COEF%SPEC_PER_TRAC( N_TRACERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%SPEC_PER_TRAC = 0

      ALLOCATE( COEF%TRAC_COEF( N_TRACERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%TRAC_COEF = 0d0

      ALLOCATE( COEF%MOLWT_KG( N_TRACERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%MOLWT_KG = 0d0

      ALLOCATE( COEF%XNUMOL( N_TRACERS ), STAT=RC )
      IF ( RC /= SMV_SUCCESS ) RETURN
      COEF%XNUMOL = 0d0

      ! Initialize drydep name
      DRYDEP_NAME(:)          = ''

      ! Put this routine name on error trace stack
      IDENT%LEV               = IDENT%LEV + 1
      IDENT%I_AM( IDENT%LEV ) = 'GC_CHUNK_INIT'

      ! Unit for logfile redirect
      IU_LOG                  = IDENT%STDOUT_LUN

      !=================================================================
      ! Read chemistry mechanism
      !=================================================================
      
      ! Read from data file mglob.dat
      CALL READER( IDENT, .TRUE., TS_CHEM, RC )
      IF ( RC /= SMV_SUCCESS ) RETURN

      ! Read "globchem.dat" chemistry mechanism
      CALL READCHEM( IDENT, RC )
      IF ( RC /= SMV_SUCCESS ) RETURN

      ! Pass # of species and reactions back as output
      N_SPECIES   = NGAS
      N_REACTIONS = NTRATES(NCS) 

      !=================================================================
      ! Return w/ error if the LISOPOH species is undefined
      !=================================================================
      IF ( OPTIONS%USE_SEC_ORG_AEROSOLS ) THEN
         IF ( ILISOPOH == 0 ) THEN
            RC           = SMV_FAILURE
            IDENT%ERRMSG = 'LISOPOH must be defined for SOA!'
            RETURN
         ENDIF
      ENDIF

      !=================================================================
      ! Set ID flags for advected tracers
      !=================================================================
      DO N = 1, N_TRACERS

         SELECT CASE( TRIM( TRACER_NAME(N) ) )
            CASE( 'NOx', 'NOX' )
               ID_TRACERS%NOx  = N
            CASE( 'Ox',  'OX'  )
               ID_TRACERS%Ox   = N
            CASE( 'PAN'  )
               ID_TRACERS%PAN  = N
            CASE( 'CO'   )
               ID_TRACERS%CO   = N
            CASE( 'ALK4' )
               ID_TRACERS%ALK4 = N
            CASE( 'ISOP' )
               ID_TRACERS%ISOP = N
            CASE( 'HNO3' )
               ID_TRACERS%HNO3 = N
            CASE( 'H2O2' )
               ID_TRACERS%H2O2 = N
            CASE( 'ACET' )
               ID_TRACERS%ACET = N
            CASE( 'MEK'  )
               ID_TRACERS%MEK  = N
            CASE( 'ALD2' )
               ID_TRACERS%ALD2 = N
            CASE( 'RCHO' )
               ID_TRACERS%RCHO = N
            CASE( 'MVK'  )
               ID_TRACERS%MVK  = N
            CASE( 'MACR' )
               ID_TRACERS%MACR = N
            CASE( 'PMN'  )
               ID_TRACERS%PMN  = N
            CASE( 'PPN'  )
               ID_TRACERS%PPN  = N
            CASE( 'R4N2' )
               ID_TRACERS%R4N2 = N
            CASE( 'PRPE' )
               ID_TRACERS%PRPE = N
            CASE( 'C3H8' )
               ID_TRACERS%C3H8 = N
            CASE( 'CH2O' )
               ID_TRACERS%CH2O = N
            CASE( 'C2H6' )
               ID_TRACERS%C2H6 = N
            CASE( 'N2O5' )
               ID_TRACERS%N2O5 = N
            CASE( 'HNO4' )
               ID_TRACERS%HNO4 = N
            CASE( 'MP'   )
               ID_TRACERS%MP   = N
            CASE( 'DMS'  )
               ID_TRACERS%DMS  = N
            CASE( 'SO2'  )
               ID_TRACERS%SO2  = N
            CASE( 'SO4'  )
               ID_TRACERS%SO4  = N
            CASE( 'SO4s' )
               ID_TRACERS%SO4s = N
            CASE( 'MSA'  )
               ID_TRACERS%MSA  = N
            CASE( 'NH3'  )
               ID_TRACERS%NH3  = N
            CASE( 'NH4'  )
               ID_TRACERS%NH4  = N
            CASE( 'NIT'  )
               ID_TRACERS%NIT  = N
            CASE( 'NITs' )
               ID_TRACERS%NITs = N
            CASE( 'BCPI' )
               ID_TRACERS%BCPI = N
            CASE( 'OCPI' )
               ID_TRACERS%OCPI = N
            CASE( 'BCPO' )
               ID_TRACERS%BCPO = N
            CASE( 'OCPO' )
               ID_TRACERS%OCPO = N
            CASE( 'ALPH' )
               ID_TRACERS%ALPH = N
            CASE( 'LIMO' )
               ID_TRACERS%LIMO = N
            CASE( 'ALCO' )
               ID_TRACERS%ALCO = N
            CASE( 'SOG1' )
               ID_TRACERS%SOG1 = N
            CASE( 'SOG2' )
               ID_TRACERS%SOG2 = N
            CASE( 'SOG3' )
               ID_TRACERS%SOG3 = N
            CASE( 'SOG4' )
               ID_TRACERS%SOG4 = N
            CASE( 'SOA1' )
               ID_TRACERS%SOA1 = N
            CASE( 'SOA2' )
               ID_TRACERS%SOA2 = N
            CASE( 'SOA3' )
               ID_TRACERS%SOA3 = N
            CASE( 'SOA4' )
               ID_TRACERS%SOA4 = N
            CASE( 'DST1' )
               ID_TRACERS%DST1 = N
            CASE( 'DST2' )
               ID_TRACERS%DST2 = N
            CASE( 'DST3' )
               ID_TRACERS%DST3 = N
            CASE( 'DST4' )
               ID_TRACERS%DST4 = N
            CASE( 'SALA' )
               ID_TRACERS%SALA = N
            CASE( 'SALC' )
               ID_TRACERS%SALC = N
            CASE DEFAULT
               ! Nothing
         END SELECT
      ENDDO

      !=================================================================
      ! Set ID flags for chemistry species (active & inactive)
      !=================================================================
      DO N = 1, NTSPEC(NCS)

         SELECT CASE( TRIM( NAMEGAS(N) ) )
            CASE( 'A3O2'      )  
               ID_SPECIES%A3O2     = N
            CASE( 'ACET'      )    
               ID_SPECIES%ACET     = N
            CASE( 'ACTA'      )    
               ID_SPECIES%ACTA     = N
            CASE( 'ALD2'      )    
               ID_SPECIES%ALD2     = N
            CASE( 'ALK4'      )    
               ID_SPECIES%ALK4     = N
            CASE( 'ATO2'      )    
               ID_SPECIES%ATO2     = N
            CASE( 'B3O2'      )    
               ID_SPECIES%B3O2     = N
            CASE( 'C2H6'      )    
               ID_SPECIES%C2H6     = N
            CASE( 'C3H8'      )    
               ID_SPECIES%C3H8     = N
            CASE( 'CH2O'      )    
               ID_SPECIES%CH2O     = N
            CASE( 'CH4'       )    
               ID_SPECIES%CH4      = N
            CASE( 'CO'        )    
               ID_SPECIES%CO       = N
            CASE( 'CO2'       )    
               ID_SPECIES%CO2      = N
            CASE( 'DMS'       )    
               ID_SPECIES%DMS      = N
            CASE( 'DRYCH2O'   )    
               ID_SPECIES%DRYCH2O  = N
            CASE( 'DRYDEP'    )    
               ID_SPECIES%DRYDEP   = N
            CASE( 'DRYH2O2'   )    
               ID_SPECIES%DRYH2O2  = N
            CASE( 'DRYHNO3'   )    
               ID_SPECIES%DRYHNO3  = N
            CASE( 'DRYN2O5'   )    
               ID_SPECIES%DRYN2O5  = N
            CASE( 'DRYNO2'    )    
               ID_SPECIES%DRYNO2   = N
            CASE( 'DRYO3'     )    
               ID_SPECIES%DRYO3    = N
            CASE( 'DRYPAN'    )    
               ID_SPECIES%DRYPAN   = N
            CASE( 'DRYPMN'    )    
               ID_SPECIES%DRYPMN   = N
            CASE( 'DRYPPN'    )    
               ID_SPECIES%DRYPPN   = N
            CASE( 'DRYR4N2'   )
               ID_SPECIES%DRYR4N2  = N
            CASE( 'EMISSION'  )
               ID_SPECIES%EMISSION = N
            CASE( 'EOH'       )
               ID_SPECIES%EOH      = N
            CASE( 'ETO2'      )
               ID_SPECIES%ETO2     = N
            CASE( 'ETP'       )
               ID_SPECIES%ETP      = N
            CASE( 'GCO3'      )
               ID_SPECIES%GCO3     = N
            CASE( 'GLCO3'     )
               ID_SPECIES%GLCO3    = N
            CASE( 'GLP'       )
               ID_SPECIES%GLP      = N
            CASE( 'GLPAN'     )
               ID_SPECIES%GLPAN    = N
            CASE( 'GLYC'      )
               ID_SPECIES%GLYC     = N
            CASE( 'GLYX'      )
               ID_SPECIES%GLYX     = N
            CASE( 'GP'        )
               ID_SPECIES%GP       = N
            CASE( 'GPAN'      )
               ID_SPECIES%GPAN     = N
            CASE( 'H'         )
               ID_SPECIES%H        = N
            CASE( 'H2'        )
               ID_SPECIES%H2       = N
            CASE( 'H2O'       )
               ID_SPECIES%H2O      = N
            CASE( 'H2O2'      )
               ID_SPECIES%H2O2     = N
            CASE( 'HAC'       )
               ID_SPECIES%HAC      = N
            CASE( 'HCOOH'     )
               ID_SPECIES%HCOOH    = N
            CASE( 'HNO2'      )
               ID_SPECIES%HNO2     = N
            CASE( 'HNO3'      )
               ID_SPECIES%HNO3     = N
            CASE( 'HNO4'      )
               ID_SPECIES%HNO4     = N
            CASE( 'HO2'       )
               ID_SPECIES%HO2      = N
            CASE( 'IALD'      )
               ID_SPECIES%IALD     = N
            CASE( 'IAO2'      )
               ID_SPECIES%IAO2     = N
            CASE( 'IAP'       )
               ID_SPECIES%IAP      = N
            CASE( 'INO2'      )
               ID_SPECIES%INO2     = N
            CASE( 'INPN'      )
               ID_SPECIES%INPN     = N
            CASE( 'ISN1'      )
               ID_SPECIES%ISN1     = N
            CASE( 'ISNO3'     )
               ID_SPECIES%ISNO3    = N
            CASE( 'ISNP'      )
               ID_SPECIES%ISNP     = N
            CASE( 'ISOP'      )
               ID_SPECIES%ISOP     = N
            CASE( 'KO2'       )
               ID_SPECIES%KO2      = N
            CASE( 'LISOPOH'   )
               ID_SPECIES%LISOPOH  = N
            CASE( 'M'         )
               ID_SPECIES%M        = N
            CASE( 'MACR'      )
               ID_SPECIES%MACR     = N
            CASE( 'MAN2'      )
               ID_SPECIES%MAN2     = N
            CASE( 'MAO3'      )
               ID_SPECIES%MAO3     = N
            CASE( 'MAOP'      )
               ID_SPECIES%MAOP     = N
            CASE( 'MAP'       )
               ID_SPECIES%MAP      = N
            CASE( 'MCO3'      )
               ID_SPECIES%MCO3     = N
            CASE( 'MEK'       )
               ID_SPECIES%MEK      = N
            CASE( 'MGLY'      )
               ID_SPECIES%MGLY     = N
            CASE( 'MNO3'      )
               ID_SPECIES%MNO3     = N
            CASE( 'MO2'       )
               ID_SPECIES%MO2      = N
            CASE( 'MOH'       )
               ID_SPECIES%MOH      = N
            CASE( 'MP'        )
               ID_SPECIES%MP       = N
            CASE( 'MRO2'      )
               ID_SPECIES%MRO2     = N
            CASE( 'MRP'       )
               ID_SPECIES%MRP      = N
            CASE( 'MSA'       )
               ID_SPECIES%MSA      = N
            CASE( 'MVK'       )
               ID_SPECIES%MVK      = N
            CASE( 'MVN2'      )
               ID_SPECIES%MVN2     = N
            CASE( 'N2'        )
               ID_SPECIES%N2       = N
            CASE( 'N2O'       )
               ID_SPECIES%N2O      = N
            CASE( 'N2O5'      )
               ID_SPECIES%N2O5     = N
            CASE( 'NH2'       )
               ID_SPECIES%NH2      = N
            CASE( 'NH3'       )
               ID_SPECIES%NH3      = N
            CASE( 'NO'        )
               ID_SPECIES%NO       = N
            CASE( 'NO2'       )
               ID_SPECIES%NO2      = N
            CASE( 'NO3'       )
               ID_SPECIES%NO3      = N
            CASE( 'O'         )
               ID_SPECIES%O        = N
            CASE( 'O1D'       )
               ID_SPECIES%O1D      = N
            CASE( 'O2'        )
               ID_SPECIES%O2       = N
            CASE( 'O2CH2OH'   )
               ID_SPECIES%O2CH2OH  = N
            CASE( 'O3'        )
               ID_SPECIES%O3       = N
            CASE( 'OH'        )
               ID_SPECIES%OH       = N
            CASE( 'PAN'       )
               ID_SPECIES%PAN      = N
            CASE( 'PMN'       )
               ID_SPECIES%PMN      = N
            CASE( 'PO2'       )
               ID_SPECIES%PO2      = N
            CASE( 'PP'        )
               ID_SPECIES%PP       = N
            CASE( 'PPN'       )
               ID_SPECIES%PPN      = N
            CASE( 'PRN1'      )
               ID_SPECIES%PRN1     = N
            CASE( 'PRPE'      )
               ID_SPECIES%PRPE     = N
            CASE( 'PRPN'      )
               ID_SPECIES%PRPN     = N
            CASE( 'R4N1'      )
               ID_SPECIES%R4N1     = N
            CASE( 'R4N2'      )
               ID_SPECIES%R4N2     = N
            CASE( 'R4O2'      )
               ID_SPECIES%R4O2     = N
            CASE( 'R4P'       )
               ID_SPECIES%R4P      = N
            CASE( 'RA3P'      )
               ID_SPECIES%RA3P     = N
            CASE( 'RB3P'      )
               ID_SPECIES%RB3P     = N
            CASE( 'RCHO'      )
               ID_SPECIES%RCHO     = N
            CASE( 'RCO3'      )
               ID_SPECIES%RCO3     = N
            CASE( 'RCOOH'     )
               ID_SPECIES%RCOOH    = N
            CASE( 'RIO1'      )
               ID_SPECIES%RIO1     = N
            CASE( 'RIO2'      )
               ID_SPECIES%RIO2     = N
            CASE( 'RIP'       )
               ID_SPECIES%RIP      = N
            CASE( 'ROH'       )
               ID_SPECIES%ROH      = N
            CASE( 'RP'        )
               ID_SPECIES%RP       = N
            CASE( 'SO2'       )
               ID_SPECIES%SO2      = N
            CASE( 'SO4'       )
               ID_SPECIES%SO4      = N
            CASE( 'VRO2'      )
               ID_SPECIES%VRO2     = N
            CASE( 'VRP'       )
               ID_SPECIES%VRP      = N
            CASE DEFAULT
               ! Nothing
         END SELECT
      ENDDO

      !=================================================================
      ! Save tracer number in the ID_TRACERS object
      ! Create the coefficient arrays that link tracers & species
      !=================================================================

      ! Write header text
      WRITE( IU_LOG, '(/,a)' ) REPEAT( '=', 79 )
      WRITE( IU_LOG, 200     )
      WRITE( IU_LOG, '(a,/)' ) REPEAT( '=', 79 )
      WRITE( IU_LOG, 210     ) 
      WRITE( IU_LOG, '(  a)' ) REPEAT( '-', 30 )

      ! Formats
 200  FORMAT(  &
     &  'ADVECTED TRACERS (==> denotes emitted constituent species)' )
 210  FORMAT( '  # Tracer          g/mole' )

      ! Loop over all advected tracers
      DO N = 1, N_TRACERS
         
         IF ( N == ID_TRACERS%NOx ) THEN
            COEF%SPEC_COEF    (N,1:4) = (/ 1d0, 1d0, 1d0, 1d0 /)
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%NO2
            COEF%SPEC_ID      (N,2  ) = ID_SPECIES%NO
            COEF%SPEC_ID      (N,3  ) = ID_SPECIES%NO3
            COEF%SPEC_ID      (N,4  ) = ID_SPECIES%HNO2
            COEF%SPEC_PER_TRAC(N    ) = 4
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%SPEC_EMITTED (N    ) = 2     ! NO is emitted
            COEF%MOLWT_KG     (N    ) = 46d-3

         ELSE IF ( N == ID_TRACERS%Ox ) THEN
            COEF%SPEC_COEF    (N,1:3) = (/ 1d0, 1d0, 2d0 /)
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%O3
            COEF%SPEC_ID      (N,2  ) = ID_SPECIES%NO2
            COEF%SPEC_ID      (N,3  ) = ID_SPECIES%NO3  
            COEF%SPEC_PER_TRAC(N    ) = 3
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%SPEC_EMITTED (N    ) = 1     ! O3 is emitted
            COEF%MOLWT_KG     (N    ) = 48d-3

         ELSE IF ( N == ID_TRACERS%PAN ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PAN
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 121d-3

         ELSE IF ( N == ID_TRACERS%CO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%CO
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%SPEC_EMITTED (N    ) = 1    ! CO is emitted
            COEF%MOLWT_KG     (N    ) = 28d-3

         ELSE IF ( N == ID_TRACERS%ALK4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 4d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ALK4
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 4d0
            COEF%SPEC_EMITTED (N    ) = 1     ! ALK4 is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%ISOP ) THEN
            COEF%SPEC_COEF    (N,1  ) = 5d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ISOP
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 5d0
            COEF%SPEC_EMITTED (N    ) = 1     ! ISOP is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%HNO3 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%HNO3
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%SPEC_EMITTED (N    ) = 1     ! HNO3 is emitted
            COEF%MOLWT_KG     (N    ) = 63d-3

         ELSE IF ( N == ID_TRACERS%H2O2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%H2O2
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 34d-3

         ELSE IF ( N == ID_TRACERS%ACET ) THEN
            COEF%SPEC_COEF    (N,1  ) = 3d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ACET
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 3d0
            COEF%SPEC_EMITTED (N    ) = 1     ! ACET is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%MEK ) THEN
            COEF%SPEC_COEF    (N,1  ) = 4d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MEK
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 4d0
            COEF%SPEC_EMITTED (N    ) = 1     ! MEK is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%ALD2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 2d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%ALD2
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 2d0
            COEF%SPEC_EMITTED (N    ) = 1     ! ALD2 is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%RCHO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%RCHO
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 58d-3

         ELSE IF ( N == ID_TRACERS%MVK ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MVK
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 70d-3

         ELSE IF ( N == ID_TRACERS%MACR ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MACR
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 70d-3

         ELSE IF ( N == ID_TRACERS%PMN ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PMN
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 147d-3

         ELSE IF ( N == ID_TRACERS%PPN ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PPN
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 135d-3

         ELSE IF ( N == ID_TRACERS%R4N2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%R4N2
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 119d-3

         ELSE IF ( N == ID_TRACERS%PRPE ) THEN
            COEF%SPEC_COEF    (N,1  ) = 3d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%PRPE
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 3d0
            COEF%SPEC_EMITTED (N    ) = 1     ! PRPE is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%C3H8 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 3d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%C3H8
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 3d0
            COEF%SPEC_EMITTED (N    ) = 1     ! C3H8 is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%CH2O ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%CH2O
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%SPEC_EMITTED (N    ) = 1     ! CH2O is emitted
            COEF%MOLWT_KG     (N    ) = 30d-3

         ELSE IF ( N == ID_TRACERS%C2H6 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 2d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%C2H6
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 2d0
            COEF%SPEC_EMITTED (N    ) = 1     ! C2H6 is emitted
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%N2O5 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%N2O5
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 105d-3

         ELSE IF ( N == ID_TRACERS%HNO4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%HNO4
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 79d-3

         ELSE IF ( N == ID_TRACERS%MP ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MP
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 48d-3

         ELSE IF ( N == ID_TRACERS%DMS ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%DMS
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 62d-3

         ELSE IF ( N == ID_TRACERS%SO2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%SO2
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 64d-3

         ELSE IF ( N == ID_TRACERS%SO4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%SO4
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 96d-3

         ELSE IF ( N == ID_TRACERS%SO4s ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 96d-3

         ELSE IF ( N == ID_TRACERS%MSA ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = ID_SPECIES%MSA
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 96d-3

         ELSE IF ( N == ID_TRACERS%NH3 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 17d-3

         ELSE IF ( N == ID_TRACERS%NH4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 18d-3

         ELSE IF ( N == ID_TRACERS%NIT ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 62d-3

         ELSE IF ( N == ID_TRACERS%NITs ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 62d-3

         ELSE IF ( N == ID_TRACERS%BCPI ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%OCPI ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%BCPO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%OCPO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 12d-3

         ELSE IF ( N == ID_TRACERS%ALPH ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 136.23d-3

         ELSE IF ( N == ID_TRACERS%LIMO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 136.23d-3

         ELSE IF ( N == ID_TRACERS%ALCO ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 142d-3

         ELSE IF ( N == ID_TRACERS%SOG1 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 150d-3

         ELSE IF ( N == ID_TRACERS%SOG2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 160d-3

         ELSE IF ( N == ID_TRACERS%SOG3 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 220d-3

         ELSE IF ( N == ID_TRACERS%SOG4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 130d-3

         ELSE IF ( N == ID_TRACERS%SOA1 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 150d-3

         ELSE IF ( N == ID_TRACERS%SOA2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 160d-3

         ELSE IF ( N == ID_TRACERS%SOA3 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 220d-3

         ELSE IF ( N == ID_TRACERS%SOA4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 130d-3

         ELSE IF ( N == ID_TRACERS%DST1 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 29d-3

         ELSE IF ( N == ID_TRACERS%DST2 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 29d-3

         ELSE IF ( N == ID_TRACERS%DST3 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 29d-3

         ELSE IF ( N == ID_TRACERS%DST4 ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 29d-3

         ELSE IF ( N == ID_TRACERS%SALA ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 36d-3

         ELSE IF ( N == ID_TRACERS%SALC ) THEN
            COEF%SPEC_COEF    (N,1  ) = 1d0
            COEF%SPEC_ID      (N,1  ) = 0
            COEF%SPEC_PER_TRAC(N    ) = 1
            COEF%TRAC_COEF    (N    ) = 1d0
            COEF%MOLWT_KG     (N    ) = 36d-3

         ELSE

            ! Invalid tracer -- return with error status
            WRITE( IDENT%ERRMSG, 300 ) N
            RC = SMV_FAILURE
            RETURN
            
         ENDIF         

         ! XNUMOL = AVO / MOLWT = ratio of molec/kg for each tracer
         COEF%XNUMOL(N) = AVO / COEF%MOLWT_KG(N)

         ! Write tracer number, name, & mol wt
         WRITE( IU_LOG, 310 ) N, TRACER_NAME(N),                         &
     &                            COEF%MOLWT_KG(N) * 1000d0

         ! If the tracer ishas constituent species ...
         IF ( COEF%SPEC_PER_TRAC(N) > 1  .or.                            &
     &        COEF%SPEC_EMITTED(N)  > 0 ) THEN
         
            ! ... print out constituent species names
            DO C = 1, COEF%SPEC_PER_TRAC(N) 

               ! Highlight emitted species with ==> arrow
               IF ( COEF%SPEC_EMITTED(N) == C ) THEN
                  WRITE( IU_LOG, 320 ) COEF%SPEC_COEF(N,C),              &
     &                                  NAMEGAS( COEF%SPEC_ID(N,C) )
               ELSE
                  WRITE( IU_LOG, 330 ) COEF%SPEC_COEF(N,C),              &
     &                                  NAMEGAS( COEF%SPEC_ID(N,C) )
               ENDIF

            ENDDO 
         ENDIF

      ENDDO

      ! Formats
 300  FORMAT( 'Tracer ', i4, ' is not defined!' )
 310  FORMAT( I3, 1x, A10, 6x, F6.1 )
 320  FORMAT( 5x, '===> ', f4.1, 1x, A4  )
 330  FORMAT( 5x, '---> ', f4.1, 1x, A4  )

      !=================================================================
      ! Initialize dry deposition
      !=================================================================
      IF ( OPTIONS%USE_DRYDEP ) THEN

         ! Call init method from "drydep_mod.f"
!v0.1         CALL INIT_DRYDEP( IDENT, ID_TRACERS, ID_DRYDEP,                 &
!v0.1     &                     COEF,  N_DRYDEP,   DRYDEP_NAME, RC )
         N_DRYDEP = 7 !v0.1 Temp fix

         ! Return w/ error if necessary
         IF ( RC /= SMV_SUCCESS ) RETURN
      ENDIF

      !=================================================================
      ! Initialize wet deposition
      !
      ! NOTE: Turning off wet deposition will also turn off the 
      !       wet scavenging of tracer in cloud updrafts. 
      !=================================================================
      IF ( OPTIONS%USE_WETDEP .or. OPTIONS%USE_CONVECTION ) THEN

         ! Call init method from "wetdep_mod.f"
!v0.1         CALL INIT_WETDEP( IDENT,     L_COLUMN,  ID_TRACERS,             &
!v0.1     &                     ID_WETDEP, N_WETDEP,  RC )

         N_WETDEP = 7 !v0.1 Temp fix

         ! Return w/ error if necessary
         IF ( RC /= SMV_SUCCESS ) RETURN
      ENDIF

      !=================================================================
      ! Flag emission & dry deposition reactions w/in SMVGEAR
      !=================================================================
      CALL SETEMDEP( IDENT,       N_TRACERS,   N_DRYDEP,                 &
     &               TRACER_NAME, DRYDEP_NAME, COEF,     RC )

      IF ( RC /= SMV_SUCCESS ) RETURN

      !=================================================================
      ! Initialize the FAST-J photolysis mechanism
      !=================================================================
      CALL INPHOT( IDENT, L_COLUMN, NPHOT, N_JV, RC )

      IF ( RC /= SMV_SUCCESS ) RETURN

      !=================================================================
      ! Initialize the DIAG_COL object for diagnostic output
      !=================================================================

      ! Set a flag if we should print the diagnostic output to a file
      DIAG_COL%DO_PRINT = OPTIONS%USE_DEBUG_PRINT

      ! Initialize the other fields of DIAG_COL
      IF ( DIAG_COL%DO_PRINT .and. IDENT%PET == 0  ) THEN
         
         ! Allocate the TRACER fields
         IF ( .not. ASSOCIATED( DIAG_COL%TRACER ) ) THEN
            ALLOCATE( DIAG_COL%TRACER(MAX_COLUMN,N_TRACERS+1), STAT=RC )
            IF ( RC /= SMV_SUCCESS ) RETURN
         ENDIF
         
         ! Allocate the tracer name fields
         IF ( .not. ASSOCIATED( DIAG_COL%NAME ) ) THEN
            ALLOCATE( DIAG_COL%NAME( N_TRACERS+1 ), STAT=RC )
            IF ( RC /= SMV_SUCCESS ) RETURN
         ENDIF

         ! Init fields
         DIAG_COL%N_DIAG   = N_TRACERS + 1
         DIAG_COL%TRACER   = 0d0
         DIAG_COL%COUNT    = 0
         DIAG_COL%FILENAME = 'diag_col.txt'
         DIAG_COL%LUN      = 750

         ! Initialize tracer names (OH is NTRACERS+1)
         DO N = 1, N_TRACERS
            DIAG_COL%NAME(N) = TRIM( TRACER_NAME(N) )
         ENDDO
         DIAG_COL%NAME(N_TRACERS+1) = 'OH'
      ENDIF

      !================================================================
      ! Succesful return!
      !================================================================

      ! Set error code to success
      RC                      = SMV_SUCCESS

      ! Remove this routine name from error trace stack
      IDENT%I_AM( IDENT%LEV ) = ''
      IDENT%LEV               = IDENT%LEV - 1

      END SUBROUTINE GC_CHUNK_INIT

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
      SUBROUTINE CONVERT_UNITS( IDENT, DIMINFO, IFLAG,                   &
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
      OPEN( DIAG_COL%LUN,     FILE=TRIM( DIAG_COL%FILENAME ),           &
     &      STATUS='UNKNOWN', IOSTAT=RC ) 
      
      ! Return if error opening file
      IF ( RC /= SMV_SUCCESS ) THEN
         IDENT%ERRMSG = 'Could not open file ' //                       &
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
      END MODULE GC_ESMF_CHUNKMOD
#endif
