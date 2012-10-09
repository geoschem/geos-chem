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
  PUBLIC  :: GC_CHUNK_INIT
  PUBLIC  :: GC_CHUNK_RUN
  PUBLIC  :: GC_CHUNK_FINAL
!
! !REMARKS:
!  The routines in this module execute only when GEOS-Chem is connected
!  to the GEOS-5 GCM via the ESMF/MAPL interface.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  09 Oct 2012 - R. Yantosca - Now pass am_I_Root to all routines
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
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
  SUBROUTINE GC_Chunk_Run( am_I_Root, GC_MET, GC_STATE, NL, NI, NJ, RC )
!
! !USES:
!
    USE GC_CHEMDR
    USE SMV_ERRCODE_MOD
    USE GC_TYPE2_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on root CPU?
    INTEGER,            INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,            INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,            INTENT(IN)    :: NL          ! # of levels
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(CHEMSTATE),    INTENT(INOUT) :: GC_STATE    ! Chemical state
    TYPE(GC_MET_LOCAL), INTENT(INOUT) :: GC_MET      ! Meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  Met field inputs from the GC_MET object have SI units.  Some GEOS-Chem 
!  lower-level routines require nonstandard units.  Units are converted and 
!  stored in local variables within this module.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Columnized & cleaned up.
!  09 Oct 2012 - R. Yantosca - Added extra comments & cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NC

    ! Number of advected tracers
    NC = SIZE(GC_STATE%TRACERS,4)

    ! Call the chemistry run method
    CALL DO_GC_CHEM( GC_STATE, GC_MET, am_I_Root, NI, NJ, NL, NC )

    ! Return code
    RC = SMV_SUCCESS

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
! !INETRFACE:
!
  SUBROUTINE GC_Chunk_Init( NI,     NJ,   NL,   GC_MET,    GC_STATE,  &
                            tsChem, nymd, nhms, am_I_Root, RC        )
!
! !USES:
!
    USE GC_INITIALIZATION_MOD
    USE GC_TYPE2_MOD
    USE SMV_ERRCODE_MOD
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: NI          ! # of longitudes
    INTEGER, INTENT(IN)  :: NJ          ! # of latitudes
    INTEGER, INTENT(IN)  :: NL          ! # of levels
    INTEGER, INTENT(IN)  :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER, INTENT(IN)  :: nhms        ! GMT time (hh:mm:ss)
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    REAL,    INTENT(IN)  :: tsChem      ! Chemistry timestep
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(CHEMSTATE)      :: GC_STATE    ! Object for chemistry state
    TYPE(GC_MET_LOCAL)   :: GC_MET      ! Object for meteorology state
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Error return code
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long     - Initial Version
!  28 Mar 2012 - M. Long     - Rewrite per structure of BCC init interface
!  09 Oct 2012 - R. Yantosca - Added comments, cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize dimensions
    CALL GC_INIT_DIMENSIONS( NI, NJ, NL )

    ! Initialize the G-C simulation and chemistry mechanism
    CALL GC_INITRUN( GC_MET, GC_STATE, tsChem, nymd, nhms, am_I_Root)

    ! Return code
    RC = SMV_SUCCESS

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
  SUBROUTINE GC_CHUNK_FINAL( am_I_Root, RC )
!
! !USES:
!
    USE SMV_ERRCODE_MOD
!
! !INPUT PARAMETERS
!
    LOGICAL, INTENT(IN) :: am_I_Root
!
! !REMARKS:
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%% NOTE: As of 10/9/12, the finalize is not done; this is a stub. %%%%
! %%%% Need to add the finalize & deallocation code ASAP.             %%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  09 Oct 2012 - R. Yantosca - Added comments & cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)   :: RC         !  Error return code
!
! !LOCAL VARIABLES:
!
!   write(*,*) 'FINALIZING'

    ! Set error code to success
    RC                      = SMV_SUCCESS

  END SUBROUTINE GC_CHUNK_FINAL
!EOC
END MODULE GC_CHUNK_MOD
#endif
