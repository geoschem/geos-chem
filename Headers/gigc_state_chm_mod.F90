!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_Indx
  PUBLIC :: Register_Species
  PUBLIC :: Register_Tracer
  PUBLIC :: Init_GIGC_State_Chm
  PUBLIC :: Cleanup_GIGC_State_Chm
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Chemistry State
  !=========================================================================
  TYPE, PUBLIC :: ChmState

     ! Advected tracers
     INTEGER,           POINTER :: Trac_Id    (:      )  ! Tracer ID #'s
     CHARACTER(LEN=14), POINTER :: Trac_Name  (:      )  ! Tracer names
     REAL*8,            POINTER :: Tracers    (:,:,:,:)  ! Tracer conc [kg]
     REAL*8,            POINTER :: Trac_Tend  (:,:,:,:)  ! Tracer tendency
     REAL*8,            POINTER :: Trac_Btend (:,:,:,:)  ! Biomass tendency

     ! Chemical species
     INTEGER,           POINTER :: Spec_Id    (:      )  ! Species ID # 
     CHARACTER(LEN=14), POINTER :: Spec_Name  (:      )  ! Species names
     REAL*8,            POINTER :: Species    (:,:,:,:)  ! Species [molec/cm3]

     ! Chemical rates & rate parameters
     REAL*8,            POINTER :: DepSav     (:,:,:  )  ! Drydep freq [1/s]

     ! Stratospheric chemistry 
     INTEGER,           POINTER :: Schm_Id    (:      )  ! Strat Chem ID #'s
     CHARACTER(LEN=14), POINTER :: Schm_Name  (:      )  ! Strat Chem Names
     REAL*8,            POINTER :: Schm_P     (:,:,:,:)  ! Strat prod [v/v/s]
     REAL*8,            POINTER :: Schm_k     (:,:,:,:)  ! Strat loss [1/s]
     INTEGER,           POINTER :: Schm_BryId (:      )  ! Bry tracer #'s
     CHARACTER(LEN=14), POINTER :: Schm_BryNam(:      )  ! Bry Names
     REAL*8,            POINTER :: Schm_BryDay(:,:,:,:)  ! Bry, Day
     REAL*8,            POINTER :: Schm_BryNit(:,:,:,:)  ! Bry, Night
 
  END TYPE ChmState
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
!  26 Oct 2012 - R. Yantosca - Add fields for stratospheric chemistry
!  26 Feb 2013 - M. Long     - Add DEPSAV to derived type ChmState
!  07 Mar 2013 - R. Yantosca - Add Register_Tracer subroutine
!  07 Mar 2013 - R. Yantosca - Now make POSITION a locally SAVEd variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: register_species
!
! !DESCRIPTION: Routine REGISTER\_SPECIES stores the names of GEOS-Chem 
!  chemical species in fields of the Chemistry State (aka State\_Chm) object.
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
!  07 Mar 2013 - R. Yantosca - Now make POSITION a locally saved variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Position index
    INTEGER,  SAVE :: POSITION = 1
    
    !======================================================================
    ! REGISTER_SPECIES begins here!
    !======================================================================

    ! We have not found the desired species yet
    Status                          = -1
    
    ! Locate the species name and ID
    State_Chm%Spec_Name( POSITION ) = TRIM( Name )
    State_Chm%Spec_Id  ( POSITION ) = Id
    
    ! Return status
    Status                          = POSITION

    ! Increment for next species
    POSITION                        = POSITION + 1
   
  END SUBROUTINE Register_Species
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Tracer
!
! !DESCRIPTION: Routine REGISTER\_TRACER stores the names of GEOS-Chem
!  advected tracers in fields of the Chemistry State (aka State\_Chm) object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_Tracer( Name, Id, State_Chm, Status )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: Name       ! Name of desired tracer
    INTEGER,          INTENT(IN)    :: Id         ! ID flag of desired tracer
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: Status     ! Success or failure
!
! !REVISION HISTORY:
!   7 Mar 2013 - R. Yantosca - Initial version, based on Register_SPecies
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Position index
    INTEGER, SAVE :: POSITION = 1

    !======================================================================
    ! REGISTER_TRACER begins here!
    !======================================================================

    ! We have not found the desired species yet
    Status                          = -1

    ! Locate the tracer name and ID
    State_Chm%Trac_Name( POSITION ) = TRIM( Name )
    State_Chm%Trac_Id  ( POSITION ) = ID

    ! Return status
    Status                          = POSITION

    ! Increment for next species
    POSITION                        = POSITION + 1

  END SUBROUTINE Register_Tracer
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE Init_GIGC_State_Chm( am_I_Root, IM,        JM,        LM,     &  
                                  nTracers,  nBioMax,   nSpecies,  nSchm,  &    
                                  nSchmBry,  Input_Opt, State_Chm, RC      )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
    USE GIGC_Input_Opt_Mod, ONLY   : OptInput    ! Derived type
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
    INTEGER,        INTENT(IN)    :: nSchm       ! # of strat chem species
    INTEGER,        INTENT(IN)    :: nSchmBry    ! # of Bry species, strat chm
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  26 Oct 2012 - R. Yantosca - Now allocate Strat_P, Strat_k fields
!  26 Oct 2012 - R. Yantosca - Add nSchem, nSchemBry as arguments
!  01 Nov 2012 - R. Yantosca - Don't allocate strat chem fields if nSchm=0
!                              and nSchmBry=0 (i.e. strat chem is turned off)
!  26 Feb 2013 - M. Long     - Now pass Input_Opt via the argument list
!  26 Feb 2013 - M. Long     - Now allocate the State_Chm%DEPSAV field
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                       :: MAX_DEP

    ! Assume success until otherwise
    RC = GIGC_SUCCESS

    ! Maximum # of drydep species
    MAX_DEP = Input_Opt%MAX_DEP

    !=====================================================================
    ! Allocate advected tracer fields
    !=====================================================================

    ALLOCATE( State_Chm%Trac_Id       (             nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Name     (             nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Tracers       ( IM, JM, LM, nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Tend     ( IM, JM, LM, nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Btend    ( IM, JM, LM, nBiomax    ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !=====================================================================
    ! Allocate chemical species fields
    !=====================================================================

    ALLOCATE( State_Chm%Spec_Id       (             nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Spec_Name     (             nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Species       ( IM, JM, LM, nSpecies   ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !=====================================================================
    ! Allocate chemical rate fields
    !=====================================================================

    ! DEPSAV is allocated in drydep_mod
    ALLOCATE( State_Chm%DEPSAV        ( IM, JM,     Max_Dep    ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

! NOTE: Comment out for now, leave for future expansion (bmy, 11/20/12)
!    !=====================================================================
!    ! Allocate stratospheric chemistry fields
!    !=====================================================================
!
!    ! Only allocate if strat chem is turned on
!    IF ( nSchm > 0 ) THEN
!
!       ALLOCATE( State_Chm%Schm_Id    (             nSchm      ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!       
!       ALLOCATE( State_Chm%Schm_Name  (             nSchm      ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!       
!       ALLOCATE( State_Chm%Schm_P     ( IM, JM, LM, nSchm      ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!       
!       ALLOCATE( State_Chm%Schm_k     ( IM, JM, LM, nSchm      ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!    
!    ENDIF
!
!    ! Only allocate if strat chem is turned on
!    IF ( nSchmBry > 0 ) THEN
!   
!       ALLOCATE( State_Chm%Schm_BryId (             nSchmBry   ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!       
!       ALLOCATE( State_Chm%Schm_BryNam(             nSchmBry   ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!
!       ALLOCATE( State_Chm%Schm_BryDay( IM, JM, LM, nSchmBry   ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!
!       ALLOCATE( State_Chm%Schm_BryNit( IM, JM, LM, nSchmBry   ), STAT=RC )
!       IF ( RC /= GIGC_SUCCESS ) RETURN
!
!    ENDIF

    !=====================================================================
    ! Initialize fields
    !=====================================================================

    ! Advected tracers
    State_Chm%Trac_Id     = 0
    State_Chm%Trac_name   = ''
    State_Chm%Tracers     = 0d0
    State_Chm%Trac_Tend   = 0d0
    State_Chm%Trac_Btend  = 0d0

    ! Dry deposition
    State_Chm%DepSav      = 0d0

    ! Chemical species
    State_Chm%Spec_Id     = 0
    State_Chm%Spec_Name   = ''
    State_Chm%Species     = 0d0

! NOTE: Comment out for now, leave for future expansion (bmy, 11/20/12)
!    ! Stratospheric chemistry    
!    IF ( nSchm > 0 ) THEN
!       State_Chm%Schm_Id     = 0
!       State_Chm%Schm_Name   = ''
!       State_Chm%Schm_P      = 0d0
!       State_Chm%Schm_k      = 0d0
!    ENDIF
!    IF ( nSchmBry > 0 ) THEN
!       State_Chm%Schm_BryId  = 0
!       State_Chm%Schm_BryNam = ''
!       State_Chm%Schm_BryDay = 0d0
!       State_Chm%Schm_BryNit = 0d0
!    ENDIF

  END SUBROUTINE Init_GIGC_State_Chm
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  15 Oct 2012 - R. Yantosca - Initial version
!  26 Oct 2012 - R. Yantosca - Now deallocate Strat_P, Strat_k fields
!  26 Feb 2013 - M. Long     - Now deallocate State_Chm%DEPSAV
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    ! Deallocate fields
    IF ( ASSOCIATED(State_Chm%Trac_Id    ) ) DEALLOCATE(State_Chm%Trac_Id    )
    IF ( ASSOCIATED(State_Chm%Trac_Name  ) ) DEALLOCATE(State_Chm%Trac_Name  )
    IF ( ASSOCIATED(State_Chm%Spec_Id    ) ) DEALLOCATE(State_Chm%Spec_Id    )
    IF ( ASSOCIATED(State_Chm%Spec_Name  ) ) DEALLOCATE(State_Chm%Spec_Name  )
    IF ( ASSOCIATED(State_Chm%Trac_Tend  ) ) DEALLOCATE(State_Chm%Trac_Tend  )
    IF ( ASSOCIATED(State_Chm%Trac_Btend ) ) DEALLOCATE(State_Chm%Trac_Btend )
    IF ( ASSOCIATED(State_Chm%Tracers    ) ) DEALLOCATE(State_Chm%Tracers    )
    IF ( ASSOCIATED(State_Chm%Species    ) ) DEALLOCATE(State_Chm%Species    )
    IF ( ASSOCIATED(State_Chm%DepSav     ) ) DEALLOCATE(State_Chm%DepSav     )

    ! NOTE: Comment out for now, leave for future expansion (bmy, 11/26/12)
    !IF ( ASSOCIATED(State_Chm%Schm_Id    ) ) DEALLOCATE(State_Chm%Schm_Id    )
    !IF ( ASSOCIATED(State_Chm%Schm_Name  ) ) DEALLOCATE(State_Chm%Schm_Name  )
    !IF ( ASSOCIATED(State_Chm%Schm_P     ) ) DEALLOCATE(State_Chm%Schm_P     )
    !IF ( ASSOCIATED(State_Chm%Schm_k     ) ) DEALLOCATE(State_Chm%Schm_k     )
    !IF ( ASSOCIATED(State_Chm%Schm_BryId ) ) DEALLOCATE(State_Chm%Schm_BryId )
    !IF ( ASSOCIATED(State_Chm%Schm_BryNam) ) DEALLOCATE(State_Chm%Schm_BryNam)
    !IF ( ASSOCIATED(State_Chm%Schm_BryDay) ) DEALLOCATE(State_Chm%Schm_BryDay)
    !IF ( ASSOCIATED(State_Chm%Schm_BryNit) ) DEALLOCATE(State_Chm%Schm_BryNit)

  END SUBROUTINE Cleanup_GIGC_State_Chm
!EOC
END MODULE GIGC_State_Chm_Mod
