! $Id: tracer_mod.f,v 1.2 2004/12/16 16:52:46 bmy Exp $
      MODULE TRACER_MOD
!
!******************************************************************************
!  Module TRACER_MOD contains the GEOS-CHEM tracer array STT plus various
!  other related quantities.  TRACER_MOD also contains inquiry functions that
!  can be used to determine the type of GEOS-CHEM simulation.
!  (bmy, 7/20/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) SIM_TYPE               : Number denoting simulation type
!  (2 ) N_TRACERS              : Number of GEOS-CHEM tracers
!  (3 ) N_MEMBERS              : Max # of constituents a tracer can have
!  (4 ) ID_TRACER              : Array of tracer numbers
!  (5 ) ID_EMITTED             : Index of which constituent has the emissions
!  (6 ) STT                    : GEOS-CHEM Tracer array [kg] 
!  (7 ) TCVV                   : Molecular weight air / molecular weight tracer
!  (8 ) TRACER_COEFF           : Coefficient of each tracer constituent
!  (9 ) TRACER_MW_G            : Tracer molecular weight [g/mole]
!  (10) TRACER_MW_KG           : Tracer molecular weight [kg/mole]
!  (11) TRACER_N_CONST         : Array of number of constituents per tracer
!  (12) TRACER_NAME            : Array of tracer names
!  (13) TRACER_CONST           : Array of names for tracer constituents
!  (14) SALA_REDGE_um          : Accum mode seasalt radii bin edges [um]
!  (15) SALC_REDGE_um          : Coarse mode seasalt radii bin edges [um]
!
!  Module Routines:
!  ============================================================================
!  (1 ) ITS_A_RnPbBe_SIM       : Returns TRUE if it's a  Rn-Pb-Be  simulation
!  (2 ) ITS_A_CH3I_SIM         : Returns TRUE if it's a  CH3I      simulation
!  (3 ) ITS_A_FULLCHEM_SIM     : Returns TRUE if it's a  fullchem  simulation
!  (4 ) ITS_A_HCN_SIM          : Returns TRUE if it's a  HCN       simulation
!  (5 ) ITS_A_COPARAM_SIM      : Returns TRUE if it's a  CO param  simulation
!  (6 ) ITS_A_TAGOX_SIM        : Returns TRUE if it's a  Tagged Ox simulation
!  (7 ) ITS_A_TAGCO_SIM        : Returns TRUE if it's a  Tagged CO simulation
!  (8 ) ITS_A_C2H6_SIM         : Returns TRUE if it's a  C2H6      simulation
!  (9 ) ITS_A_CH4_SIM          : Returns TRUE if it's a  CH4       simulation
!  (10) ITS_AN_AEROSOL_SIM     : Returns TRUE if it's an aerosol   simulation
!  (11) ITS_A_MERCURY_SIM      : Returns TRUE if it's a  mercury   simulation
!  (12) ITS_NOT_COPARAM_OR_CH4 : Returns TRUE if it's not CO param or CH4 
!  (13) CHECK_STT              : Checks STT array for NaN, Inf, or negatives
!  (14) INIT_TRACER            : Allocates and zeroes all module arrays
!  (15) CLEANUP_TRACER         : Deallocates all module arrays
!
!  Module Routines:
!  ============================================================================
!  (1 ) error_mod.f : Module containing I/O error and NaN check routines 
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER                        :: SIM_TYPE
      INTEGER                        :: N_TRACERS
      INTEGER,           PARAMETER   :: N_MEMBERS = 10
      INTEGER,           ALLOCATABLE :: ID_TRACER(:)
      INTEGER,           ALLOCATABLE :: ID_EMITTED(:)
      INTEGER,           ALLOCATABLE :: TRACER_N_CONST(:)
      REAL*8,            ALLOCATABLE :: STT(:,:,:,:)
      REAL*8,            ALLOCATABLE :: TCVV(:)
      REAL*8,            ALLOCATABLE :: TRACER_COEFF(:,:)
      REAL*8,            ALLOCATABLE :: TRACER_MW_G(:)
      REAL*8,            ALLOCATABLE :: TRACER_MW_KG(:)
      CHARACTER(LEN=14), ALLOCATABLE :: TRACER_NAME(:)
      CHARACTER(LEN=14), ALLOCATABLE :: TRACER_CONST(:,:)

      ! Define seasalt radii bin edges [um] here since these
      ! need to be used both in "seasalt_mod.f" and "drydep_mod.f"
      REAL*8                         :: SALA_REDGE_um(2)
      REAL*8                         :: SALC_REDGE_um(2)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      FUNCTION ITS_A_RnPbBe_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_RnPbBe_SIM returns TRUE if we are doing a GEOS-CHEM
!  Rn-Pb-Be simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 1 )

      ! Return to calling program
      END FUNCTION ITS_A_RnPbBe_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_CH3I_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_CH3I_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH3I (Methyl Iodide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH3I_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 2 )

      ! Return to calling program
      END FUNCTION ITS_A_CH3I_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_FULLCHEM_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_FULLCHEM_SIM returns TRUE if we are doing a GEOS-CHEM
!  full chemistry/aerosol simulation (i.e. via SMVGEAR). (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_FULLCHEM_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 3 )

      ! Return to calling program
      END FUNCTION ITS_A_FULLCHEM_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_HCN_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_HCN_SIM returns TRUE if we are doing a GEOS-CHEM
!  HCN (Hydrogen Cyanide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_HCN_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 4 )

      ! Return to calling program
      END FUNCTION ITS_A_HCN_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_COPARAM_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_COPARAM_SIM returns TRUE if we are doing a GEOS-CHEM CO 
!  simulation with parameterized OH simulation (cf. B. Duncan). (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CO_OH_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 5 )

      ! Return to calling program
      END FUNCTION ITS_A_COPARAM_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_TAGOX_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_TAGOX_SIM returns TRUE if we are doing a GEOS-CHEM
!  Tagged Ox simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGOX_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 6 )

      ! Return to calling program
      END FUNCTION ITS_A_TAGOX_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_TAGCO_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_TAGCO_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH3I (Methyl Iodide) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_TAGCO_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 7 )

      ! Return to calling program
      END FUNCTION ITS_A_TAGCO_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_C2H6_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_C2H6_SIM returns TRUE if we are doing a GEOS-CHEM
!  C2H6 (Ethane) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_RnPbBe_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 8 )

      ! Return to calling program
      END FUNCTION ITS_A_C2H6_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_CH4_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_CH4_SIM returns TRUE if we are doing a GEOS-CHEM
!  CH4 (Methane) simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_CH4_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 9 )

      ! Return to calling program
      END FUNCTION ITS_A_CH4_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_AN_AEROSOL_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_SULFATE_SIM returns TRUE if we are doing a GEOS-CHEM
!  offline Sulfate/Carbon/dust/seasalt aerosol simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_AN_AEROSOL_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 10 )

      ! Return to calling program
      END FUNCTION ITS_AN_AEROSOL_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_A_MERCURY_SIM() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_A_MERCURY_SIM returns TRUE if we are doing a GEOS-CHEM
!  Hg0/Hg2/HgP offline mercury simulation. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_A_MERCURY_SIM begins here!
      !=================================================================
      VALUE = ( SIM_TYPE == 11 )

      ! Return to calling program
      END FUNCTION ITS_A_MERCURY_SIM

!------------------------------------------------------------------------------

      FUNCTION ITS_NOT_COPARAM_OR_CH4() RESULT( VALUE )
!
!******************************************************************************
!  Function ITS_NOT_COPARAM_OR_CH4 returns TRUE if we are doing a GEOS-CHEM
!  simulation other than CO w/ Parameterized OH or CH4. (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! Local variables
      LOGICAL :: VALUE
      
      !=================================================================
      ! ITS_NOT_COPARAM_OR_CH4 begins here!
      !=================================================================
      VALUE = ( SIM_TYPE /= 5 .and. SIM_TYPE /= 9 )

      ! Return to calling program
      END FUNCTION ITS_NOT_COPARAM_OR_CH4

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_STT( LOCATION )
!
!******************************************************************************
!  Subroutine CHECK_STT checks the STT tracer array for negative values,
!  NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message. (bmy, 3/8/01, 7/15/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1) LOCATION (CHARACTER) : String describing location of error in code
!
!  NOTES:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!  (4 ) Moved from "error_mod.f" to "tracer_mod.f" (bmy, 7/15/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD

#     include "CMN_SIZE"           ! Size parameters

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION

      ! Local variables
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N
      
      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
         IF ( STT(I,J,L,N) < 0d0 ) THEN 
!$OMP CRITICAL
            LNEG = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNEG .or. LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a )

      ! Return to calling program
      END SUBROUTINE CHECK_STT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRACER
!
!******************************************************************************
!  Subroutine INIT_TRACER initializes all module arrays (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! Allocate arrays 
      !=================================================================
      ALLOCATE( ID_TRACER( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_TRACER' )
      ID_TRACER = 0

      ALLOCATE( ID_EMITTED( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ID_EMITTED' )
      ID_EMITTED = 0

      ALLOCATE( TCVV( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCVV' )
      TCVV = 0d0

      ALLOCATE( TRACER_NAME( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_NAME' )
      TRACER_NAME = ''

      ALLOCATE( TRACER_MW_G( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_G' )
      TRACER_MW_G = 0d0

      ALLOCATE( TRACER_MW_KG( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_MW_KG' )
      TRACER_MW_KG = 0d0

      ALLOCATE( TRACER_COEFF( N_TRACERS, N_MEMBERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_COEFF' )
      TRACER_COEFF = 0d0

      ALLOCATE( TRACER_CONST( N_TRACERS, N_MEMBERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_CONST' )
      TRACER_CONST = ''

      ALLOCATE( TRACER_N_CONST( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_N_CONST' )
      TRACER_N_CONST = 0

      ALLOCATE( STT( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT' )
      STT = 0d0      

      ! Return to calling program
      END SUBROUTINE INIT_TRACER

!-----------------------------------------------------------------------------
      
      SUBROUTINE CLEANUP_TRACER
!
!******************************************************************************
!  Subroutine CLEANUP_TRACER deallocates all module arrays (bmy, 7/15/04)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRACER begins here!
      !=================================================================
      IF ( ALLOCATED( ID_TRACER      ) ) DEALLOCATE( ID_TRACER      )
      IF ( ALLOCATED( ID_EMITTED     ) ) DEALLOCATE( ID_EMITTED     )
      IF ( ALLOCATED( TCVV           ) ) DEALLOCATE( TCVV           )
      IF ( ALLOCATED( TRACER_NAME    ) ) DEALLOCATE( TRACER_NAME    )
      IF ( ALLOCATED( TRACER_COEFF   ) ) DEALLOCATE( TRACER_COEFF   )
      IF ( ALLOCATED( TRACER_CONST   ) ) DEALLOCATE( TRACER_CONST   )
      IF ( ALLOCATED( TRACER_N_CONST ) ) DEALLOCATE( TRACER_N_CONST )
      IF ( ALLOCATED( TRACER_MW_G    ) ) DEALLOCATE( TRACER_MW_G    )
      IF ( ALLOCATED( TRACER_MW_KG   ) ) DEALLOCATE( TRACER_MW_KG   )
      IF ( ALLOCATED( STT            ) ) DEALLOCATE( STT            )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRACER

!-----------------------------------------------------------------------------

      ! End of module
      END MODULE TRACER_MOD
