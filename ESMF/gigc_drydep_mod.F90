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
! !IROUTINE: gigc_do_chem
!
! !DESCRIPTION: Routine GIGC\_DO\_CHEM is the chemistry driver for the
!  Grid-Independent GEOS-Chem (aka "GIGC") model.  This routine is called by
!  the Run method from the ESMF interface to the GEOS-5 GCM.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GIGC_Do_DryDep( am_I_Root, NI,        NJ,         &
                               NL,        Input_Opt,             &
                               State_Chm, State_Met, RC            )
!
! !USES:
!

      USE LOGICAL_MOD
      USE DRYDEP_MOD,         ONLY : DO_DRYDEP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState

! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
!>    INTEGER,        INTENT(IN)    :: NCNST       ! # of constituents
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
!
! !LOCAL VARIABLES:
!

! TESTING SECTION
!<><><><><><><><><><><><><><><><><><><><><><>
     !CALL GIGC_Dump_Config( am_I_Root )


!<><><><><><><><><><><><><><><><><><><><><><>
! TEMPORARY INITIALIZATION SECTION
!>      IIPAR      = NI
!>      JJPAR      = NJ
!>      LLPAR      = NL
!>      LLTROP     = NL
!>      LLTROP_FIX = NL

      !======================================================================
      ! Set 2-D variables
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: This is a stopgap measure for testing.  Eventually we   %%%
      ! %%% will carry the Meteorology State and Chemistry State objects  %%%
      ! %%% down to all G-C routines.   In order to continue testing the  %%%
      ! %%% GIGC without disrupting existing workflow, we must populate   %%%
      ! %%% G-C module arays from the Meteorology and Chemistry States.   %%%
      ! %%% (bmy, 10/25/12)                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      ! Met fields

      !======================================================================
      ! Set 3-D variables
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: This is a stopgap measure for testing.  Eventually we   %%%
      ! %%% will carry the Meteorology State and Chemistry State objects  %%%
      ! %%% down to all G-C routines.   In order to continue testing the  %%%
      ! %%% GIGC without disrupting existing workflow, we must populate   %%%
      ! %%% G-C module arays from the Meteorology and Chemistry States.   %%%
      ! %%% (bmy, 10/25/12)                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      ! Met fields

      ! Constituents

      !======================================================================
      ! Call the GEOS-Chem Dry Deposition routines
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%%                                                               %%%
      ! %%%                                                               %%%
      ! %%%                                                               %%%
      ! %%%                                                               %%%
      ! %%%                                                               %%%
      ! %%%                                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      !### Debug, print values in v/v before chem
      IF ( am_I_Root ) THEN
         WRITE(6,*) '##### GIGC_Do_DryDep, LDRYD'
         WRITE(6,*) LDRYD
      ENDIF

      ! If we are doing chemistry
      IF ( LDRYD ) THEN

         ! Convert Units

         ! Call the GEOS-Chem Dry Deposition routines
         CALL Do_DryDep( am_I_Root, NI,        NJ,        NL,  &
                         Input_Opt, State_Met, State_Chm, RC )
         RC = 1

         ! Convert Units
           
      ENDIF

      !======================================================================
      ! Reset 3-D variables
      !
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! %%% NOTE: This is a stopgap measure for testing.  Eventually we   %%%
      ! %%% will carry the Meteorology State and Chemistry State objects  %%%
      ! %%% down to all G-C routines.   In order to continue testing the  %%%
      ! %%% GIGC without disrupting existing workflow, we must populate   %%%
      ! %%% G-C module arays from the Meteorology and Chemistry States.   %%%
      ! %%% (bmy, 10/25/12)                                               %%%
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !======================================================================

      ! Save chemistry output for next timestep
      ! State_Chm%Species = CSPEC_FULL

      !### Debug
      !IF ( am_I_Root ) THEN
      !   WRITE(6,*) '##### GIGC_Do_DryDep, TRC_OX after chem'
      !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
      !ENDIF

    END SUBROUTINE GIGC_Do_DryDep
!EOC
END MODULE GIGC_DryDep_Mod
#endif
