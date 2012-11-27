#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_state_chm_mod
!
! !DESCRIPTION: Module GIGC\_STATE\_CHM\_MOD contains the derived type
!  used to define the Chemistry State object for the Grid-Independent 
!  GEOS-Chem implementation (abbreviated "GIGC").
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory 
!  to the Chemistry State object.  The chemistry state object is not defined
!  in this module.  It must be be declared as variable in the top-level 
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_State_Chm_Mod
!
! USES:
!
  IMPLICIT NONE
# include "define.h"
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_Indx
  PUBLIC :: Register_Species
  PUBLIC :: Init_GIGC_State_Chm
  PUBLIC :: Cleanup_GIGC_State_Chm
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Chemistry State
  !=========================================================================
  TYPE, PUBLIC :: ChmState

     ! Values from "tracer_mod.F"
     INTEGER,           POINTER :: Trac_Id   (:      )  ! Tracer ID flags
     CHARACTER(LEN=14), POINTER :: Trac_Name (:      )  ! Tracer names
     INTEGER,           POINTER :: Spec_Id   (:      )  ! Chemical species flags
     CHARACTER(LEN=14), POINTER :: Spec_Name (:      )  ! Chemical species name

     ! Concentrations & tendencies
     REAL*8,            POINTER :: Trac_Tend (:,:,:,:)  ! Tracer tendency
     REAL*8,            POINTER :: Trac_Btend(:,:,:,:)  ! Biomass tendency
     REAL*8,            POINTER :: Tracers   (:,:,:,:)  ! Tracer conc [kg]
     REAL*8,            POINTER :: Species   (:,:,:,:)  ! Species [molec/cm3]

  END TYPE ChmState

  !=========================================================================
  ! Other variables
  !=========================================================================

  ! Position value used for registering CSPEC parameters in the chemical state
  INTEGER,  SAVE :: POSITION = 1 
!
! !REMARKS:
!  -----------------------------------------------------------------------
!   FULLCHEM Simulation Emissions
!   Done | Units? | Routine
!   Y/P  --> Yes or Partially(fix needed)
!           good  --> Verifed that they are in Kg/s
!  -----------------------------------------------------------------------
!   Y    |        | CALL COMPUTE_BIOMASS_EMISSIONS
!   Y    |        | CALL EMISS_STREETS_ANTHRO_05x0666
!   Y    |        | CALL EMISS_STREETS_ANTHRO
!   Y    |        | CALL EMISS_EDGAR( YEAR, MONTH )
!   Y    |  good  | CALL EMISS_RETRO
!   Y    |  good* | CALL EMISS_EPA_NEI
!   Y    |        | CALL EMISS_VISTAS_ANTHRO
!   Y    |        | CALL EMISS_BRAVO
!   Y    |        | CALL EMISS_EMEP_05x0666
!   Y    |        | CALL EMISS_EMEP
!   Y    |        | CALL EMISS_CAC_ANTHRO_05x0666
!   Y    |        | CALL EMISS_CAC_ANTHRO
!   Y    |        | CALL EMISS_EPA_NEI
!   Y    |        | CALL EMISS_NEI2005_ANTHRO_05x0666
!   Y    |        | CALL EMISS_NEI2005_ANTHRO
!   Y    |        | CALL EMISS_ARCTAS_SHIP( YEAR )
!   Y    |        | CALL EMISS_ICOADS_SHIP
!        |        | CALL EMISSDR
!   Y    |        | CALL EMISSSEASALT
!   Y    |        | CALL EMISSSULFATE --> Be sure there's no PBL mixing
!   Y    |        | CALL EMISSCARBON --> Be sure there's no PBL mixing
!   Y    |        | CALL EMISSDUST --> Be sure there's no PBL mixing
!   Y    |        | AIRCRAFT_NOX
!   Y    |        | LIGHTNING_NOX
!   Y    |        | SOIL_NOX
!        |        | BIOFUEL_BURN (NOx and CO)
!  -----------------------------------------------------------------------
!   Notes:
!   LNLPBL Switch --> NEEDS TO BE ON (>=1)
!                 --> But, does VDIFF need to be turned off? 
!   NOT ALL EMISSIONS ARE JUST AT SFC, e.g. SO2
!   LFUTURE --> How to deal with this? e.g. EDGAR emissions
!   REGIONAL EMISSIONS OVERWRITE GLOBAL!!!!! DEAL WITH THIS!
!                                                                             .
!   STT<-->CSPEC mapped in PARTITION
!                                                                             .
!   KEEP EMISSIONS FROM UPDATING STT DIRECTLY
!                                                                             .
!   NEI EMISSIONS: BIOFUEL EMISSIONS ARE NOT 'REALLY' BIOFUEL.
!                  AS IN THERE'S NO IDBF'SPEC' INDEX.
!  
!  -----------------------------------------------------------------------
!   FULLCHEM Simulation Chemistry Routines
!   Done | Units? | Routine
!   Y/P  --> Yes or Partially(fix needed)
!           good  --> Verifed that they are in Kg
!  -----------------------------------------------------------------------
!        |        | CALL CHEMDR
!        |        | CALL CHEMSEASALT
!        |        | CALL CHEMSULFATE
!        |        | CALL DO_ISOROPIAII
!        |        | CALL DO_RPMARES
!        |        | CALL CHEMCARBON
!        |        | CALL CHEMDUST
!        |        | CALL DRYFLX     
!        |        | CALL DIAGOH
!        |        | CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
!  -----------------------------------------------------------------------
!   Notes:
!   (1) If STT IS TIGHTLY LINKED TO CHEM_STATE, THEN THE ONLY CHANGE
!       NEEDED IN "DO_CHEMISTRY" IS TO PASS "CHEM_STATE" IN AND OUT
!   (2) 1st PASS, WORKS ONLY WITH FULLCHEM, NOT APM OR ANY ADD-ON SIM
!       OPTIONS.
!  -----------------------------------------------------------------------
!                                                                             
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on "gc_type2_mod.F90"
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_indx 
!
! !DESCRIPTION: Function GET\_INDX returns the index of an advected tracer or 
!  chemical species contained in the chemistry state object by name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Indx( name, allIds, allNames ) RESULT( Indx )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: name         ! Species or tracer name
    CHARACTER(LEN=*), INTENT(IN) :: allNames(:)  ! List of species/tracer names
    INTEGER,          INTENT(IN) :: allIDs(:)    ! List of species/tracer IDs 
!
! !RETURN VALUE:
!
    INTEGER                      :: Indx         ! Index of this species 
!
! !REVISION HISTORY: 
!  09 Oct 2012 - M. Long     - Initial version, based on gc_esmf_utils_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    ! Initialize
    Indx= -1

    ! Loop over all species names
    DO M = 1, SIZE( allNames )

       ! Return the index of the sought-for species
       IF( TRIM( name ) == TRIM( allNames(M) ) ) THEN
          Indx = allIDs(M)
          EXIT
       ENDIF

    ENDDO

  END FUNCTION Get_Indx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Species
!
! !DESCRIPTION: Routine REGISTER\_SPECIES registers the 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_Species( Name, Id, State_Chm, Status )
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN)    :: Name       ! Name of desired species
    INTEGER,          INTENT(IN)    :: Id         ! ID flag of desired species
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: Status     ! Success or failure
!
! !REMARKS:
!   This routine is called from SETTRACE in tracerid_mod.F.
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version, based on gc_esmf_type_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    !write(*,*) 'POSITION:', POSITION
    
    ! We have not found the desired species yet
    Status                          = -1
    
    ! Locate the species name and ID 
    State_Chm%Spec_Name( POSITION ) = Name
    State_Chm%Spec_Id  ( POSITION ) = ID
    
    ! Return status
    Status                          = POSITION

    ! Increment for next species
    POSITION                        = POSITION + 1
   
  END SUBROUTINE Register_Species
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gigc_state_chm
!
! !DESCRIPTION: Routine INIT\_GIGC\_STATE\_CHM allocates and initializes the 
!  pointer fields of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_GIGC_State_Chm( am_I_Root, IM,        JM,       &
                                  LM,        nTracers,  nBioMax,  &
                                  nSpecies,  State_Chm, RC       )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: nTracers    ! # advected tracers
    INTEGER,        INTENT(IN)    :: nBioMax     ! # biomass burning tracers
    INTEGER,        INTENT(IN)    :: nSpecies    ! # chemical species    
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  In the near future we will put some error trapping on the allocations
!  so that we can stop the simulation if the allocations cannot be made.
! 
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Renamed from gc_type2_mod.F90
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! Assume success until otherwise
    RC = GIGC_SUCCESS

    !--------------------------------
    ! Allocate ID fields
    !--------------------------------
    ALLOCATE( State_Chm%Trac_Id  ( nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Name( nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Spec_Id  ( nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Spec_Name( nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !--------------------------------
    ! Allocate concentration fields
    !--------------------------------
    ALLOCATE( State_Chm%Tracers   ( IM, JM, LM, nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Species   ( IM, JM, LM, nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !--------------------------------
    ! Allocate tendency fields
    !--------------------------------
    ALLOCATE( State_Chm%Trac_Tend ( IM, JM, LM, nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Btend( IM, JM, LM, nBiomax    ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !--------------------------------
    ! Initialize fields
    !--------------------------------
    State_Chm%Trac_Id    = 0
    State_Chm%Trac_name  = ''
    State_Chm%Spec_Id    = 0
    State_Chm%Spec_Name  = ''
    State_Chm%Trac_Tend  = 0d0
    State_Chm%Trac_Btend = 0d0
    State_Chm%Tracers    = 0d0
    State_Chm%Species    = 0d0

  END SUBROUTINE Init_GIGC_State_Chm
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gigc_state_chm
!
! !DESCRIPTION: Routine CLEANUP\_GIGC\_STATE\_CHM deallocates the fields 
!  of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_GIGC_State_Chm( am_I_Root, State_Chm, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod 
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Return code
!
! !REMARKS:
!  For now the am_I_Root and RC arguments are not used.  We include these
!  for consistency and also to facilitate future expansion. (bmy, 10/16/12)
!
! !REVISION HISTORY: 
!  15 Oct 2012 - R. Yantosca -  Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    IF ( ASSOCIATED( State_Chm%Trac_Id ) ) THEN
       DEALLOCATE( State_Chm%Trac_Id )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Trac_Name  ) ) THEN
       DEALLOCATE( State_Chm%Trac_Name  )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Spec_Id ) ) THEN
       DEALLOCATE( State_Chm%Spec_Id )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Spec_Name) ) THEN
       DEALLOCATE( State_Chm%Spec_Name )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Trac_Tend ) ) THEN
       DEALLOCATE( State_Chm%Trac_Tend  )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Trac_Btend ) ) THEN
       DEALLOCATE( State_Chm%Trac_Btend )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Tracers ) ) THEN
       DEALLOCATE( State_Chm%Tracers )
    ENDIF

    IF ( ASSOCIATED( State_Chm%Species ) ) THEN
       DEALLOCATE( State_Chm%Species )
    ENDIF

  END SUBROUTINE Cleanup_GIGC_State_Chm
!EOC
END MODULE GIGC_State_Chm_Mod
#endif
