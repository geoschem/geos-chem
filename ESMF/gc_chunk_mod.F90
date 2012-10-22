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
MODULE GC_Chunk_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GC_Chunk_Init
  PUBLIC :: GC_Chunk_Run
  PUBLIC :: GC_Chunk_Final
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
! !IROUTINE: gc_chunk_init
!
! !DESCRIPTION: Subroutine GC\_CHUNK\_INIT calls the various 
!  initialization routines that read the setup files for the GEOS-Chem
!  chunk code.  Also, ID flags for advected tracers, chemical species,
!  dry deposition species and wet deposition species are defined.
!\\
!\\
! !INETRFACE:
!
  SUBROUTINE GC_Chunk_Init( NI,     NJ,   NL,   State_Met, State_Chm,  &
                            tsChem, nymd, nhms, am_I_Root, RC        )
!
! !USES:
!
    USE GC_Initialization_Mod
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE GIGC_State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
    INTEGER,        INTENT(IN)    :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER,        INTENT(IN)    :: nhms        ! GMT time (hh:mm:ss)
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    REAL,           INTENT(IN)    :: tsChem      ! Chemistry timestep
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State Object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State Object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  28 Mar 2012 - M. Long     - Rewrite per structure of BCC init interface
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize GEOS-Chem dimensions w/ the dimensions on this PET
    CALL GC_Init_Dimensions( NI, NJ, NL )

    ! Initialize the G-C simulation and chemistry mechanism
    CALL GC_InitRun( State_Met, State_Chm, tsChem, nymd, nhms, am_I_Root, RC )

  END SUBROUTINE GC_Chunk_Init
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
  SUBROUTINE GC_Chunk_Run( am_I_Root, State_Met, State_Chm, NL, NI, NJ, RC )
!
! !USES:
!
    USE Gc_ChemDr
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE GIGC_State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on root CPU?
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
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
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initialf Version
!  09 Oct 2012 - R. Yantosca - Added extra comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  17 Oct 2012 - R. Yantosca - Need to call AIRQNT before chemistry
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NC

    ! Assume success
    RC = GIGC_SUCCESS

    ! Number of advected tracers
    NC = SIZE( State_Chm%Tracers, 4 )

    ! Call the chemistry run method
    CALL Do_GC_Chem( State_Chm, State_Met, am_I_Root, NI, NJ, NL, NC )

    ! Return code
    RC = SMV_SUCCESS

  END SUBROUTINE GC_Chunk_Run

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
  SUBROUTINE GC_Chunk_Final( State_Met, State_Chm, am_I_Root, RC )
!
! !USES:
!
    USE GC_Finalization_Mod
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod,    ONLY : ChmState
    USE GIGC_State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  09 Oct 2012 - R. Yantosca - Added comments & cosmetic changes
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!   write(*,*) 'FINALIZING'

    ! Assume succes
    RC = GIGC_SUCCESS

    ! Finalize GEOS-Chem
    CALL GC_Finalize( State_Met, State_Chm, am_I_Root, RC )

  END SUBROUTINE GC_Chunk_Final
!EOC
END MODULE GC_Chunk_Mod
#endif
