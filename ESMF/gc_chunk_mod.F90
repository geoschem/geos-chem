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
! !PRIVATE MEMBER FUNCTIONS:
!

!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!
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
      SUBROUTINE GC_CHUNK_RUN( GC_MET, GC_STATE, NL, NI, NJ, RC )
!      SUBROUTINE GC_CHUNK_RUN( NLx, NIx, NJx, RC )
!
! !USES:
!
        USE GC_CHEMDR
        USE SMV_ERRCODE_MOD
!        USE GC_INITIALIZATION_MOD
        USE GC_TYPE2_MOD
!
! !REMARKS:
!  Met field inputs from the MET_1d object have SI units.  Some GEOS-Chem 
!  lower-level routines require nonstandard units.  Units are converted and 
!  stored in local variables within this module.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Columnized & cleaned up.

!EOP
!------------------------------------------------------------------------------
!BOC
        INTEGER, INTENT(OUT) :: RC
        INTEGER, INTENT(IN)  :: NL, NI, NJ
        TYPE(CHEMSTATE),    INTENT(INOUT) :: GC_STATE
        TYPE(GC_MET_LOCAL), INTENT(INOUT) :: GC_MET

        INTEGER, PARAMETER :: NCNST=1
        REAL*8 :: DLON(NI,NJ), DLAT(NI,NJ)
        REAL*8 :: LON(NI,NJ), LAT(NI,NJ)

!
! !LOCAL VARIABLES:
!

      DLAT = 1.
      DLON = DLAT

        CALL DO_GC_CHEM(GC_STATE, GC_MET, NI, NJ, NL, NCNST=SIZE(GC_STATE%TRACERS,4))

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
      SUBROUTINE GC_CHUNK_INIT( NL, NI, NJ, GC_MET, GC_STATE, RC )
!
! !USES:
!
        USE GC_INITIALIZATION_MOD
        USE GC_TYPE2_MOD
        USE SMV_ERRCODE_MOD
!
! !REVISION HISTORY: 
!  18 Jul 2011 - M. Long - Initial Version
!  28 Mar 2012 - M. Long - Rewrite per structure of BCC initialization interface
!EOP
!------------------------------------------------------------------------------
!BOC
        INTEGER,            INTENT(OUT) :: RC          ! Error return code
        INTEGER,            INTENT(IN)  :: NI,NJ,NL
        TYPE(CHEMSTATE)                 :: GC_STATE
        TYPE(GC_MET_LOCAL)              :: GC_MET
!
! !LOCAL VARIABLES:
!

        CALL GC_INIT_DIMENSIONS(NI,NJ,NL)

        CALL GC_INITRUN(GC_MET, GC_STATE) !cheminit is in here

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
      SUBROUTINE GC_CHUNK_FINAL( RC )
!
! !USES:
!
        USE SMV_ERRCODE_MOD
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
! !OUTPUT PARAMETERS:
!
        INTEGER,             INTENT(OUT)   :: RC         !  Error return code
!
! !LOCAL VARIABLES:
!
!        write(*,*) 'FINALIZING'

      ! Set error code to success
      RC                      = SMV_SUCCESS

      END SUBROUTINE GC_CHUNK_FINAL
!EOC
      END MODULE GC_CHUNK_MOD
#endif
