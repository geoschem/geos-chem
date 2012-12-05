#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  run, and finalize methods for the ESMF interface to the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Chunk_Mod
!
! !USES:
!      
  USE Mapping_Mod, ONLY : MapWeight

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Chunk_Init
  PUBLIC :: GIGC_Chunk_Run
  PUBLIC :: GIGC_Chunk_Final
!
! !REMARKS:
!  The routines in this module execute only when GEOS-Chem is connected
!  to the GEOS-5 GCM via the ESMF/MAPL interface.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  09 Oct 2012 - R. Yantosca - Now pass am_I_Root to all routines
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chunk_mod.F90
!  01 Nov 2012 - R. Yantosca - Now pass Input Options object to routines
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
  TYPE(GC_DIAG)                 :: DIAG_COL

  ! Derived type objects
  TYPE(MapWeight),      POINTER :: mapping(:,:) => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_init
!
! !DESCRIPTION: Subroutine GIGC\_CHUNK\_INIT is the ESMF init method for
!  the Grid-Independent GEOS-Chem (aka "GIGC").  This routine calls lower-
!  level routines to allocate arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Init( am_I_Root, I_LO,      J_LO,      I_HI,       &
                              J_HI,      IM,        JM,        LM,         &
                              IM_WORLD,  JM_WORLD,  LM_WORLD,  nymd,       &
                              nhms,      tsChem,    lonCtr,    latCtr,     &
                              latEdg,    Input_Opt, State_Chm, State_Met,  &
                              RC                                          )
!
! !USES:
!
    USE ESMF_Mod,                ONLY : ESMF_KIND_R4
    USE GIGC_Initialization_Mod, ONLY : GIGC_Init_Simulation
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod,      ONLY : OptInput
    USE GIGC_State_Chm_Mod,      ONLY : ChmState
    USE GIGC_State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,            INTENT(IN)    :: I_LO        ! Min lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_LO        ! Min lat index, this CPU
    INTEGER,            INTENT(IN)    :: I_HI        ! Max lon index, this CPU
    INTEGER,            INTENT(IN)    :: J_HI        ! Max lat index, this CPU
    INTEGER,            INTENT(IN)    :: IM          ! # lons, this CPU
    INTEGER,            INTENT(IN)    :: JM          ! # lats, this CPU
    INTEGER,            INTENT(IN)    :: LM          ! # levs, this CPU
    INTEGER,            INTENT(IN)    :: IM_WORLD    ! # lons, global grid
    INTEGER,            INTENT(IN)    :: JM_WORLD    ! # lats, global grid
    INTEGER,            INTENT(IN)    :: LM_WORLD    ! # levs, global grid
    INTEGER,            INTENT(IN)    :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER,            INTENT(IN)    :: nhms        ! GMT time (hh:mm:ss)
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep
    REAL(ESMF_KIND_R4), INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: latEdg(:,:) ! Lat centers [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),     INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState),     INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState),     INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  28 Mar 2012 - M. Long     - Rewrite per structure of BCC init interface
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Init
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  01 Nov 2012 - R. Yantosca - Reordered arguments for clarit
!  28 Nov 2012 - M. Long     - Now pass lonCtr, latCtr, latEdg as arguments
!                              to routine GIGC_Init_Simulation
!  03 Dec 2012 - R. Yantosca - Now call Init_CMN_SIZE (in CMN_SIZE_mod.F)
!                              instead of GIGC_Init_Dimensions to initialize
!                              the size parameters.
!  03 Dec 2012 - R. Yantosca - Rename NI, NJ, NL to IM, JM, LM for clarity
!  03 Dec 2012 - R. Yantosca - Now pass I_LO, J_LO, I_HI, J_HI, IM_WORLD, 
!                              JM_WORLD, LM_WORLD via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    !======================================================================
    ! Initialize the G-C simulation and chemistry mechanism
    USE DAO_Mod,            ONLY : Convert_Units
    USE GIGC_ChemDr,        ONLY : GIGC_Do_Chem
    USE GIGC_DryDep_Mod,    ONLY : GIGC_Do_DryDep
   !USE DryDep_Mod,         ONLY : Do_DryDep
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE Pressure_Mod,       ONLY : EXTERNAL_PEDGE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # of lons on this CPU
    INTEGER,        INTENT(IN)    :: JM          ! # of lats on this CPU
    INTEGER,        INTENT(IN)    :: LM          ! # of levs on this CPU
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  Met field inputs from the GC_MET object have SI units.  Some GEOS-Chem 
!  lower-level routines require nonstandard units.  Units are converted and 
!  stored in local variables within this module.
!                                                                             .
!  Need to add better error-handling.
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added extra comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  17 Oct 2012 - R. Yantosca - Need to call AIRQNT before chemistry
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Run
!  25 Oct 2012 - R. Yantosca - Now pass RC to GIGC_DO_CHEM
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt to GIGC_Do_Chem
!  13 Nov 2012 - M. Long     - Added Dry Deposition method
!  29 Nov 2012 - R. Yantosca - Now block off calls to GIGC_DO_DRYDEP and
!                              GIGC_DO_CHEM w/ the appropriate logical flags
!  04 Dec 2012 - R. Yantosca - Now convert units of State_Chm%TRACERS here
!                              instead of in lower-level routines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NC

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC             = GIGC_SUCCESS

    ! Exit if we are doing neither drydep nor chemistry
    IF ( ( .not. Input_Opt%LDRYD ) .and. ( .not. Input_Opt%LCHEM ) ) THEN
       RETURN
    ENDIF

    ! Pressure @ level edges [hPa]
    EXTERNAL_PEDGE = State_Met%PEDGE     

    ! # of tracers
    NC             = Input_Opt%N_TRACERS

    !=======================================================================
    ! The tracer concentrations in State_Chm%TRACERS state have units 
    ! of [v/v] (since these come from the ESMF internal state).  We need
    ! to convert these to [kg] before calling the GEOS-Chem drydep and
    ! chemistry routines.
    !=======================================================================
    CALL Convert_Units    ( IFLAG      = 2,                    &
                            N_TRACERS  = Input_Opt%N_TRACERS,  &
                            TCVV       = Input_Opt%TCVV,       &
                            AD         = State_Met%AD,         &
                            STT        = State_Chm%Tracers    )

    !=======================================================================
    ! Call the run method of the GEOS-Chem dry deposition package
    !
    ! %%%% NOTE: Eventually we can call DO_DRYDEP directly from here
    !=======================================================================

    IF ( Input_Opt%LDRYD ) THEN
       CALL GIGC_Do_DryDep( am_I_Root  = am_I_Root,  &   ! Are we on root CPU?
                            NI         = IM,         &   ! # lons on this CPU
                            NJ         = JM,         &   ! # lats on this CPU
                            NL         = LM,         &   ! # levs on this CPU
                            Input_Opt  = Input_Opt,  &   ! Input Options
                            State_Chm  = State_Chm,  &   ! Chemistry State
                            State_Met  = State_Met,  &   ! Meteorology State
                            RC         = RC         )    ! Success or failure
    ENDIF

    !=======================================================================
    ! Call the run method of the GEOS-Chem dry chemistry package
    ! 
    ! %%%% NOTE: Eventually we can call DO_CHEMISTRY directly from here
    !=======================================================================
    IF ( Input_Opt%LCHEM ) THEN
       CALL GIGC_Do_Chem(   am_I_Root  = am_I_Root,  &   ! Are we on root CPU?
                            NI         = IM,         &   ! # lons on this CPU
                            NJ         = JM,         &   ! # lats on this CPU
                            NL         = LM,         &   ! # levels on this CPU
                            NCNST      = NC,         &   ! # of advected tracers
                            Input_Opt  = Input_Opt,  &   ! Input Options
                            State_Chm  = State_Chm,  &   ! Chemistry State
                            State_Met  = State_Met,  &   ! Meteorology State
                            RC         = RC         )    ! Success or failure
    ENDIF


    !=======================================================================
    ! Convert the tracer concentrations in State_Chm%TRACERS back to
    ! [v/v] so that they can be stored in the ESMF internal state
    ! for the next chemistry timestep.
    !=======================================================================
    CALL Convert_Units    ( IFLAG     = 1,                    &
                            N_TRACERS = Input_Opt%N_TRACERS,  &
                            TCVV      = Input_Opt%TCVV,       &
                            AD        = State_Met%AD,         &
                            STT       = State_Chm%Tracers    )

  END SUBROUTINE GIGC_Chunk_Run
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_chunk_final
!
! !DESCRIPTION: Subroutine GC\_CHUNK\_FINAL deallocates pointers and
!  arrays used in the chemistry. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Chunk_Final( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Finalization_Mod
    USE GIGC_Input_Opt_Mod,    ONLY : OptInput
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE GIGC_State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Chunk_Final
!  01 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume succes
    RC = GIGC_SUCCESS

    ! Finalize GEOS-Chem
    CALL GIGC_Finalize( am_I_Root = am_I_Root,  &  ! Are we on the root CPU?
                        Input_Opt = Input_Opt,  &  ! Input Options
                        State_Chm = State_Chm,  &  ! Chemistry State
                        State_Met = State_Met,  &  ! Meteorology State
                        RC        = RC         )   ! Success or failure?

  END SUBROUTINE GIGC_Chunk_Final
!EOC
END MODULE GIGC_Chunk_Mod
#endif
