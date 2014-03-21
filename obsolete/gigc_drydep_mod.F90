#if defined( ESMF_ )
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_chemdr
!
! !DESCRIPTION: Module GC\_CHEMDR is the "bridge" between the ESMF interface
!  to the GEOS-5 GCM and the GEOS-Chem chemistry routines.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_DryDep_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Do_DryDep
!
! !REMARKS:
!  This mo
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers, F90 indentation
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET to State_Met
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chemdr.F90
!  08 Nov 2012 - R. Yantosca - Deleted obsolete, commented-out stuff
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_do_drydep
!
! !DESCRIPTION: Routine GIGC\_DO\_DRYDEP is the chemistry driver for the
!  Grid-Independent GEOS-Chem (aka "GIGC") model.  This routine is called by
!  the Run method from the ESMF interface to the GEOS-5 GCM.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Do_DryDep( am_I_Root, NI,        NJ,        NL, &
                             Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState

! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    INTEGER,        INTENT(IN)    :: NI          ! # of lons on this CPU
    INTEGER,        INTENT(IN)    :: NJ          ! # of lats on this CPU
    INTEGER,        INTENT(IN)    :: NL          ! # of levs on this CPU
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Nov 2012 - M. Long     - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS

    !### Debug, print values in v/v before chem
    IF ( am_I_Root ) THEN
       WRITE(6,*) '##### GIGC_Do_DryDep, LDRYD'
       WRITE(6,*) Input_Opt%LDRYD
    ENDIF
    
    !=======================================================================
    ! Set 3-D variables
    !=======================================================================

    !=======================================================================
    ! Call the GEOS-Chem dry deposition routines
    !=======================================================================
    IF ( Input_Opt%LDRYD ) THEN

       ! Call the GEOS-Chem Dry Deposition routines
       CALL Do_DryDep( am_I_Root = am_I_Root,            &
                       NI        = NI,                   &
                       NJ        = NJ,                   &
                       NL        = NL,                   &
                       Input_Opt = Input_Opt,            &
                       State_Met = State_Met,            &
                       State_Chm = State_Chm,            &
                       RC        = RC                    &
                     )
    ENDIF
    !=======================================================================
   ! Reset 3-D variables
   !=======================================================================
    !### Debug
   !IF ( am_I_Root ) THEN
   !   WRITE(6,*) '##### GIGC_Do_DryDep, TRC_OX after chem'
   !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
   !ENDIF

    END SUBROUTINE GIGC_Do_DryDep
!EOC
END MODULE GIGC_DryDep_Mod
#endif
