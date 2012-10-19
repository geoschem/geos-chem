#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_type2_mod
!
! !DESCRIPTION: Module GC\_TYPE2\_MOD defines the derived types that will 
!  represent the "socket" with external modules. These types will be 
!  implemented thru GEOS-Chem in a top-down format to tightly control memory
!  allocation and access both within and outside of GEOS-Chem.
!\\
!\\
!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. 
!  It will remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
MODULE GC_Type2_Mod
!
! USES:
!
  IMPLICIT NONE
# include "define.h"
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_Tracer_Id
  PUBLIC :: Init_Chemistry_State
  PUBLIC :: Cleanup_Chemistry_State
  PUBLIC :: Register_Species
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for chemical state
  !=========================================================================
  PUBLIC :: ChemState
  TYPE   :: ChemState

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

  END TYPE ChemState

  !=========================================================================
  ! Derived type for physics state
  !=========================================================================
  PUBLIC :: PhyState
  TYPE      PhyState

     ! Scalars
     INTEGER             :: BEGIN_I
     INTEGER             :: END_I
     INTEGER             :: BEGIN_J
     INTEGER             :: END_J
     INTEGER             :: NINDX
     INTEGER             :: NJNDX

     ! 1-D arrays
     REAL*8, ALLOCATABLE :: LAT     (:  )  ! LATITUDE (DEG)
     REAL*8, ALLOCATABLE :: LON     (:  )  ! LONGITUDE (DEG)
     REAL*8, ALLOCATABLE :: PS      (:  )  ! SURFACE PRESSURE
     REAL*8, ALLOCATABLE :: PSDRY   (:  )  ! DRY SURFACE PRESSURE
     REAL*8, ALLOCATABLE :: PHIS    (:  )  ! SURFACE GEOPOTENTIAL
     REAL*8, ALLOCATABLE :: ULAT    (:  )  ! UNIQUE LATITUDES  (DEG)
     REAL*8, ALLOCATABLE :: ULON    (:  )  ! UNIQUE LONGITUDES (DEG)

     ! 2-D arrays
     REAL*8, ALLOCATABLE :: T       (:,:)   ! TEMPERATURE (K)
     REAL*8, ALLOCATABLE :: U       (:,:)   ! ZONAL WIND (M/S)
     REAL*8, ALLOCATABLE :: V       (:,:)   ! MERIDIONAL WIND (M/S)
     REAL*8, ALLOCATABLE :: OMEGA   (:,:)   ! VERTICAL PRESSURE VELOCITY (PA/S) 
     REAL*8, ALLOCATABLE :: PMID    (:,:)   ! MIDPOINT PRESSURE (PA) 
     REAL*8, ALLOCATABLE :: PMIDDRY (:,:)   ! MIDPOINT PRESSURE DRY (PA) 
     REAL*8, ALLOCATABLE :: PDEL    (:,:)   ! LAYER THICKNESS (PA)
     REAL*8, ALLOCATABLE :: PDELDRY (:,:)   ! LAYER THICKNESS DRY (PA)
     REAL*8, ALLOCATABLE :: RPDEL   (:,:)   ! 1/LAYER THICKNESS (PA)
     REAL*8, ALLOCATABLE :: RPDELDRY(:,:)   ! 1/LAYER THICKNESS DRY (PA)
     REAL*8, ALLOCATABLE :: UZM     (:,:)   ! ZONAL WIND FOR QBO (M/S)
     REAL*8, ALLOCATABLE :: ZM      (:,:)   ! GEOPOTENTIAL HEIGHT ABOVE 
                                            !  SURFACE AT MIDPOINTS (M)

! Leave commented out for now         
!         REAL*8, DIMENSION(PCOLS,PVER+1)           :: &
!              PINT,    &  ! INTERFACE PRESSURE (PA)
!              PINTDRY, &  ! INTERFACE PRESSURE DRY (PA) 
!              LNPINT,  &  ! LN(PINT)
!              LNPINTDRY,& ! LOG INTERFACE PRESSURE DRY (PA) 
!              ZI          ! GEOPOTENTIAL HEIGHT ABOVE SURFACE AT INTERFACES (M)
!         
!         INTEGER, DIMENSION(PCOLS) :: &
!              LATMAPBACK, &! MAP FROM COLUMN TO UNIQUE LAT FOR THAT COLUMN
!              LONMAPBACK, &! MAP FROM COLUMN TO UNIQUE LON FOR THAT COLUMN
!              CID         ! UNIQUE COLUMN ID
!         INTEGER :: ULATCNT, & ! NUMBER OF UNIQUE LATS IN CHUNK
!              ULONCNT     ! NUMBER OF UNIQUE LONS IN CHUNK

  END TYPE PhyState

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
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers
!  15 Oct 2012 - R. Yantosca - Added routine CLEANUP_CHEMSTATE
!  16 Oct 2012 - R. Yantosca - Add CSPEC to Chemical state fields
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
! !IROUTINE: get_tracer_id
!
! !DESCRIPTION: Routine GET\_TRACER\_ID returns a tracer's ID number given
!  its name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Tracer_Id( SPC_NAME ) RESULT( TRACER_ID )
!
! !USES:
!
    USE TRACER_MOD, ONLY : N_TRACERS 
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*), INTENT(IN) :: SPC_NAME
!
! !RETURN VALUE:
!
    INTEGER                      :: TRACER_ID
! 
! !REMARKS:
!  Not sure why this is commented out, maybe this has been moved to
!  gc_utils_mod.F90.
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers
!  15 Oct 2012 - R. Yantosca - Now pass return value via RESULT, which is
!                              required by some F90 compilers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    TRACER_ID = -1
    DO M = 1,N_TRACERS
!         IF( TRIM( SPC_NAME ) == TRIM( CHEM_STATE%TRAC_NAME(M) ) ) THEN
!            GET_TRACER_ID = M
!            EXIT
!         END IF
    ENDDO

  END FUNCTION Get_Tracer_Id
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Species
!
! !DESCRIPTION: Routine REGISTER\_SPECIES finds the
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
    TYPE(CHEMSTATE),  INTENT(INOUT) :: State_Chm  ! Object for chemical state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: Status     ! Success or failure
!
! !REMARKS:
!   This routine is called from SETTRACE in tracerid_mod.F.
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed chemical state object to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed this subroutine to Register_Species
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
! !IROUTINE: init_chemistry_state
!
! !DESCRIPTION: Routine INIT\_CHEMISTRY\_STATE allocates and initializes the 
!  pointer fields of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Chemistry_State( State_Chm, am_I_Root, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,    ONLY : IIPAR, JJPAR, LLPAR, NBIOMAX
    USE COMODE_LOOP_MOD, ONLY : IGAS
    USE TRACER_MOD,      ONLY : N_TRACERS
    USE SMV_ErrCode_Mod 
!
! !INPUT PARAMETERS:
! 
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(CHEMSTATE), INTENT(INOUT) :: State_Chm    ! Obj for chemistry state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC           ! Return code
!
! !REMARKS:
!  In the near future we will put some error trapping on the allocations
!  so that we can stop the simulation if the allocations cannot be made.
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - R. Yantosca - Initial version
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  15 Oct 2012 - R. Yantosca - Add STAT=AS to the allocation calls
!  15 Oct 2012 - R. Yantosca - Declare CHEM_STATE with INTENT(INOUT)
!  15 Oct 2012 - R. Yantosca - Add am_I_Root, RC parameters
!  16 Oct 2012 - R. Yantosca - Rename local chem state obj to State_Chm
!  16 Oct 2012 - R. Yantosca - Add CSPEC to Chemical state fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! Assume success until otherwise
    RC = SMV_SUCCESS

    ! Allocate 1-D fields
    ALLOCATE( State_Chm%Trac_Id  ( N_TRACERS+1 ), STAT=AS )
    ALLOCATE( State_Chm%Trac_Name( N_TRACERS+1 ), STAT=AS )
    ALLOCATE( State_Chm%Spec_Id  ( IGAS        ), STAT=AS )
    ALLOCATE( State_Chm%Spec_Name( IGAS        ), STAT=AS )

    ! Allocate 3-D fields
    ALLOCATE( State_Chm%Trac_Tend ( IIPAR,JJPAR,LLPAR,N_TRACERS+1 ), STAT=AS )
    ALLOCATE( State_Chm%Trac_Btend( IIPAR,JJPAR,LLPAR,NBIOMAX     ), STAT=AS )
    ALLOCATE( State_Chm%Tracers   ( IIPAR,JJPAR,LLPAR,N_TRACERS+1 ), STAT=AS )
    ALLOCATE( State_Chm%Species   ( IIPAR,JJPAR,LLPAR,IGAS        ), STAT=AS )

    ! Initialize fields of the chemical state object
    State_Chm%Trac_Id    = 0
    State_Chm%Trac_name  = ''
    State_Chm%Spec_Id    = 0
    State_Chm%Spec_Name  = ''
    State_Chm%Trac_Tend  = 0d0
    State_Chm%Trac_Btend = 0d0
    State_Chm%Tracers    = 0d0
    State_Chm%Species    = 0d0

  END SUBROUTINE Init_Chemistry_State
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_chemstate
!
! !DESCRIPTION: Routine CLEANUP\_CHEMISTRY_STATE deallocates the fields of the
!  chemical state object
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Chemistry_State( State_Chm, am_I_Root, RC )
!
! !USES:
!
    USE Smv_ErrCode_Mod 
!
! !INPUT PARAMETERS:
! 
    LOGICAL,         INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(CHEMSTATE), INTENT(INOUT) :: State_Chm    ! Obj for chemistry state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC           ! Return code
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
    RC = SMV_SUCCESS

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

  END SUBROUTINE Cleanup_Chemistry_State
!EOC
END MODULE GC_Type2_Mod
#endif
