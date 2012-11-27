#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_finalization_mod
!
! !DESCRIPTION: Module GIGC\_FINALIZATION\_MOD is the module that contains 
!  the GEOS-Chem finalize methods for the ESMF framework.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Finalization_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  PUBLIC :: GIGC_Finalize
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_finalization_mod.F90
!  02 Nov 2012 - R. Yantosca - Now cleanup the Input Options object
!  27 Nov 2012 - R. Yantosca - Remove GIGC_CleanUp_GeosChem

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
! !IROUTINE: gigc_finalize
!
! !DESCRIPTION: GIGC\_Finalize calls the cleanup routines which deallocate
!  memory from the State objects of the Grid-Independent GEOS-Chem (aka
!  "GIGC") code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Finalize( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - R. Yantosca - Initial version
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Finalize
!  02 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!  27 Nov 2012 - R. Yantosca - Now call subroutine CLEANUP to deallocate
!                              the various GEOS-Chem module arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Assume success
    RC = GIGC_SUCCESS

    ! Deallocate fields of the Input Options object
    CALL Cleanup_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_GIGC_State_Chm( am_I_Root, State_Chm, RC )

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_GIGC_State_Met( am_I_Root, State_Met, RC )

    ! Deallocate all other GEOS-Chem allocatable arrays
    CALL CleanUp( am_I_Root, RC )

  END SUBROUTINE GIGC_Finalize
!EOC
END MODULE GIGC_Finalization_Mod
#endif
