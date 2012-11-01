#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_environment_mod
!
! !DESCRIPTION: Module GIGC\_ENVIRONMENT\_MOD establishes the runtime 
!  environment for the Grid-Independent GEOS-Chem (aka "GIGC") model.  It is 
!  designed to receive model parameter and geophysical environment information 
!  and allocate memory based upon it.
!\\
!\\
!  It provides routines to do the following:
!
! \begin{itemize}
! \item Allocate geo-spatial arrays
! \item Initialize met. field derived type.
! \item Initialize CHEM, PHYS, and EMISSIONS states
! \end{itemize}

!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. 
!  It will remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_Environment_Mod
!
! !USES
!        
  IMPLICIT NONE
# include "define.h"

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Allocate_All
  PUBLIC :: GIGC_Init_All
!
! !REMARKS:
!  For consistency, we should probably move the met state initialization
!  to the same module where the met state derived type is contained.
!
! !REVISION HISTORY:
!  26 Jan 2012 - M. Long     - Created module file
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  19 Oct 2012 - R. Yantosca - Removed routine INIT_LOCAL_MET, this is now
!                              handled in Headers/gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_environment_mod.F90
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
! !IROUTINE: gigc_allocate_all
!
! !DESCRIPTION: Subroutine GIGC\_ALLOCATE\_ALL allocates all LAT/LON 
!  ALLOCATABLE arrays for global use by the GEOS-Chem either as a standalone 
!  program or module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Allocate_All( am_I_Root, RC )
!
! !USES:
!
    USE CMN_DEP_MOD,       ONLY : SET_CMN_DEP_MOD
    USE CMN_NOX_MOD,       ONLY : SET_CMN_NOX_MOD
    USE CMN_O3_MOD,        ONLY : SET_CMN_O3_MOD
    USE CMN_MOD,           ONLY : SET_CMN_MOD
    USE CMN_FJ_MOD,        ONLY : SET_CMN_FJ_MOD
    USE CMN_SIZE_MOD,      ONLY : SET_CMN_SIZE_MOD
    USE CMN_DIAG_MOD,      ONLY : SET_CMN_DIAG_MOD
    USE COMODE_LOOP_MOD,   ONLY : SET_COMODE_LOOP_MOD
    USE COMMSOIL_MOD,      ONLY : SET_COMMSOIL_MOD
    USE GIGC_ERRCODE_MOD
    USE JV_CMN_MOD,        ONLY : SET_JV_CMN_MOD
    USE VDIFF_PRE_MOD,     ONLY : SET_VDIFF_PRE_MOD

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking and exit upon failure.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  17 Oct 2012 - R. Yantosca - Add am_I_Root, RC as arguments
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Allocate_All
!  30 Oct 2012 - R. Yantosca - Now pass am_I_Root, RC to SET_COMMSOIL_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
    CALL SET_CMN_SIZE_MOD
    CALL SET_CMN_DEP_MOD
    CALL SET_CMN_DIAG_MOD
    CALL SET_CMN_NOX_MOD
    CALL SET_CMN_O3_MOD
    CALL SET_CMN_MOD
    CALL SET_CMN_FJ_MOD
    CALL SET_COMMSOIL_MOD   ( am_I_Root, RC )
    CALL SET_COMODE_LOOP_MOD( am_I_Root, RC )
    CALL SET_JV_CMN_MOD
    
    CALL SET_VDIFF_PRE_MOD
          
  END SUBROUTINE GIGC_Allocate_All
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_init_all
!
! !DESCRIPTION: Subroutine GIGC\_INIT\_ALL initializes the top-level data 
!  structures that are either passed to/from GC or between GC components 
!  (emis->transport->chem->etc)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Init_All( State_Met, State_Chm, am_I_Root, RC ) 
!
! !USES:
!
    USE CMN_Size_Mod,       ONLY : IIPAR, JJPAR, LLPAR, NBIOMAX
    USE Comode_Loop_Mod,    ONLY : IGAS
    USE GIGC_ErrCode_Mod
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
    USE Tracer_Mod,         ONLY : N_TRACERS
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Need to add better error checking, currently we just return upon error.
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  16 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE  argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Call Init_Chemistry_State (in gc_type2_mod.F90,
!                              which was renamed from INIT_CHEMSTATE)
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_errcode_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference IGAS in Headers/comode_loop_mod.F
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Init_All
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! Initialize object for met fields
    !=======================================================================
    CALL Init_GIGC_State_Met( am_I_Root  = am_I_Root,   &
                              IM         = IIPAR,       &
                              JM         = JJPAR,       &
                              LM         = LLPAR,       &
                              State_Met  = State_Met,   &
                              RC         = RC          )

    ! Return upon error
    IF ( RC /= GIGC_SUCCESS ) RETURN

    !=======================================================================
    ! Initialize object for chemical state
    !=======================================================================
    CALL Init_GIGC_State_Chm( am_I_Root  = am_I_Root,   &
                              IM         = IIPAR,       &
                              JM         = JJPAR,       &
                              LM         = LLPAR,       &
                              nTracers   = N_TRACERS,   &
                              nBioMax    = NBIOMAX,     &
                              nSpecies   = IGAS,        &
                              State_Chm  = State_Chm,   &
                              RC         = RC          )
    
    ! Return upon error
    IF ( RC /= GIGC_SUCCESS ) RETURN

  END SUBROUTINE GIGC_Init_All
!EOC           
END MODULE GIGC_Environment_Mod
#endif
