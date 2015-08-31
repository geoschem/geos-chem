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
  USE Precision_Mod    ! GEOS-Chem precision types 
  USE Species_Mod      ! For species database object

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

     ! Number of species
     INTEGER                    :: nSpecies             ! # of species
     INTEGER                    :: nAdvect              ! # of advected species
     INTEGER                    :: nDrydep              ! # of drydep species
     INTEGER                    :: nWetDep              ! # of wetdep species

     ! Physical properties about tracers & species
     TYPE(SpcPtr),      POINTER :: SpcData(:)           ! Species database

     ! Advected tracers
     INTEGER,           POINTER :: Trac_Id    (:      ) ! Tracer ID #'s
     CHARACTER(LEN=14), POINTER :: Trac_Name  (:      ) ! Tracer names
     REAL(fp),          POINTER :: Tracers    (:,:,:,:) ! Tracer conc [kg]

     ! Chemical species
     INTEGER,           POINTER :: Spec_Id    (:      ) ! Species ID # 
     CHARACTER(LEN=14), POINTER :: Spec_Name  (:      ) ! Species names
     REAL(fp),          POINTER :: Species    (:,:,:,:) ! Species [molec/cm3]

#if defined( ESMF_ )
     ! Chemical rates & rate parameters
     INTEGER,           POINTER :: JLOP       (:,:,:  ) ! 1-D SMVGEAR index
     INTEGER,           POINTER :: JLOP_PREV  (:,:,:  ) ! JLOP, prev timestep
#endif

  END TYPE ChmState
!
! !REMARKS:
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
!  19 May 2014 - C. Keller   - Removed Trac_Btend. DepSav array covers now
!                              all species.
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Dec 2014 - R. Yantosca - Keep JLOP and JLOP_PREV for ESMF runs only
!  28 Aug 2015 - R. Yantosca - Remove strat chemistry fields, these are now
!                              handled by the HEMCO component
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
  SUBROUTINE Init_GIGC_State_Chm( am_I_Root, IM,        JM,        &   
                                  LM,        nTracers,  nSpecies,  &
                                  Input_Opt, State_Chm, RC        )
!
! !USES:
!
    USE Comode_Loop_Mod,      ONLY : ILONG, ILAT, IPVERT
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod,   ONLY : OptInput
    USe Species_Database_Mod, ONLY : Init_Species_Database
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: nTracers    ! # advected tracers
    INTEGER,        INTENT(IN)    :: nSpecies    ! # chemical species
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
!  11 Dec 2014 - R. Yantosca - Remove TRAC_TEND and DEPSAV fields
!  28 Aug 2015 - R. Yantosca - Remove stratospheric chemistry fields; 
!                              these are all now read in via HEMCO
!  28 Aug 2015 - R. Yantosca - Also initialize the species database object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Assume success until otherwise
    RC = GIGC_SUCCESS

    !=====================================================================
    ! Allocate advected tracer fields
    !=====================================================================

    ALLOCATE( State_Chm%Trac_Id       (             nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Trac_Name     (             nTracers+1 ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ALLOCATE( State_Chm%Tracers       ( IM, JM, LM, nTracers+1 ), STAT=RC )
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

#if defined( ESMF_ )
    !=====================================================================
    ! Allocate chemical rate fields
    !=====================================================================

    ! Keep this here for now -- FLEXCHEM will remove this (bmy, 12/11/14)
    ALLOCATE( State_Chm%JLOP( ILONG, ILAT, IPVERT ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Chm%JLOP = 0
    
    ! Keep this here for now -- FLEXCHEM will remove this (bmy, 12/11/14)
    ALLOCATE( State_Chm%JLOP_PREV( ILONG, ILAT, IPVERT ), STAT=RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN
    State_Chm%JLOP_PREV = 0
#endif

    !=====================================================================
    ! Initialize fields
    !=====================================================================

    ! Number of species
    State_Chm%nSpecies    = 0
    State_Chm%nAdvect     = 0
    State_Chm%nDryDep     = 0
    State_Chm%nWetDep     = 0

    ! Advected tracers
    State_Chm%Trac_Id     = 0
    State_Chm%Trac_name   = ''
    State_Chm%Tracers     = 0e+0_fp

    ! Chemical species
    State_Chm%Spec_Id     = 0
    State_Chm%Spec_Name   = ''
    State_Chm%Species     = 0e+0_fp

    !=====================================================================
    ! Populate the species database object field
    ! (assumes Input_Opt has already been initialized)
    !=====================================================================
    CALL Init_Species_Database( am_I_Root = am_I_Root,          &
                                Input_Opt = Input_Opt,          &
                                SpcData   = State_Chm%SpcData,  &
                                RC        = RC                 )

    !=====================================================================
    ! Determine the number of advected, drydep, wetdep, and total species
    !=====================================================================
    State_Chm%nSpecies = SIZE( State_Chm%SpcData )


!    State_Chm%nAdvect  = MAXVAL( State_Chm%SpcData%Info%AdvectID )
!    State_Chm%nDryDep  = MAXVAL( State_Chm%SpcData%Info%DryDepID )
!    State_Chm%nWetDep  = MAXVAL( State_Chm%SpcData%Info%WetDepID )

    print*, 'nSpecies : ', State_Chm%nSpecies
    print*, 'nAdvect  : ', State_Chm%nAdvect
    print*, 'nDryDep  : ', State_Chm%nDryDep
    print*, 'nWetDep  : ', State_Chm%nWetDep
    
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
    USE Species_Database_Mod, ONLY : Cleanup_Species_Database
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
!  11 Dec 2014 - R. Yantosca - Remove TRAC_TEND and DEPSAV fields
!  28 Aug 2015 - R. Yantosca - Remove stratospheric chemistry fields; 
!                              these are all now read in via HEMCO
!  28 Aug 2015 - R. Yantosca - Also initialize the species database object
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
    IF ( ASSOCIATED(State_Chm%Tracers    ) ) DEALLOCATE(State_Chm%Tracers    )
    IF ( ASSOCIATED(State_Chm%Species    ) ) DEALLOCATE(State_Chm%Species    )

#if defined( ESMF_ )
    ! Keep these here for now, FLEXCHEM will remove these (bmy, 12/11/14)
    IF ( ASSOCIATED(State_Chm%JLOP       ) ) DEALLOCATE(State_Chm%JLOP       )
    IF ( ASSOCIATED(State_Chm%JLOP_PREV  ) ) DEALLOCATE(State_Chm%JLOP_PREV  )
#endif    

    ! Deallocate the species database object field
    CALL Cleanup_Species_Database( am_I_Root, State_Chm%SpcData, RC )

  END SUBROUTINE Cleanup_GIGC_State_Chm
!EOC
END MODULE GIGC_State_Chm_Mod
