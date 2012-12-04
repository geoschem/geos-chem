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
MODULE GIGC_ChemDr
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GIGC_Do_Chem
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
  SUBROUTINE GIGC_Do_Chem( am_I_Root, NI,        NJ,         &
                           NL,        NCNST,     Input_Opt,  &
                           State_Chm, State_Met, RC         )
!
! !USES:
!
    USE CHEMISTRY_MOD,      ONLY : DO_CHEMISTRY
    USE CMN_SIZE_MOD
    USE COMODE_MOD,         ONLY : CSPEC_FULL
    USE COMODE_LOOP_MOD
    USE DAO_MOD
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_Test_Utils
    USE GIGC_Chem_Utils
    USE PRESSURE_MOD,       ONLY : EXTERNAL_PEDGE
    USE TRACER_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    INTEGER,        INTENT(IN)    :: NI          ! # of longitudes
    INTEGER,        INTENT(IN)    :: NJ          ! # of latitudes
    INTEGER,        INTENT(IN)    :: NL          ! # of levels
    INTEGER,        INTENT(IN)    :: NCNST       ! # of constituents
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
!  NOTE: Eventually we will replace all met field arrays with the
!  values from the meteorology state.  For now we need to copy these
!  over. (bmy, 10/15/12)
!
! !REVISION HISTORY:
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX headers, F90 indentation
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE to State_Chm
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET to State_Met
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Do_Chem
!  25 Oct 2012 - R. Yantosca - Now convert units of State_Chm%TRACERS from 
!                              [v/v] -> [kg] before calling G-C chemistry 
!                              (and from [kg] -> [v/v] after chemistry)
!  08 Nov 2012 - R. Yantosca - Now pass the Input_Opt object to DO_CHEMISTRY
!  27 Nov 2012 - R. Yantosca - Now use met fields directly from State_Met
!  04 Dec 2012 - R. Yantosca - Now do unit conversion one level higher
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !=======================================================================
    ! Initialize
    !=======================================================================
!-----------------------------------------------------------------------------
! Prior to 12/4/12:
! This is now done when we initialize the grid (bmy, 12/4/12)
!!<><><><><><><><><><><><><><><><><><><><><><>
!! TEMPORARY INITIALIZATION SECTION
!      IIPAR      = NI
!      JJPAR      = NJ
!      LLPAR      = NL
!      LLTROP     = NL
!      LLTROP_FIX = NL
!-----------------------------------------------------------------------------
    RRATE      = 0.E0
    TRATE      = 0.E0

    !=======================================================================
    ! Set 3-D variables
    !=======================================================================

    ! Met fields
    EXTERNAL_PEDGE     = State_Met%PEDGE     ! Pressure @ level edges [hPa]

    ! Constituents
    CSPEC_FULL         = State_Chm%Species   ! Chemical species

    !=======================================================================
    ! Call the GEOS-Chem Chemistry routines
    !=======================================================================

    !### Debug, print values in v/v before chem
    !IF ( am_I_Root ) THEN
    !   WRITE(6,*) '##### GIGC_DO_CHEM, TRC_OX before chem [v/v]'
    !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
    !ENDIF

    ! If we are doing chemistry
    IF ( Input_Opt%LCHEM ) THEN

!-----------------------------------------------------------------------------
! Prior to 12/4/12:
! Now convert units in the higher level routine GIGC_Chunk_Run (bmy, 12/4/12)
!         ! The tracer concentrations in State_Chm%TRACERS state have units 
!         ! of [v/v] (since these come from the ESMF internal state).  We need
!         ! to convert these to [kg] before calling the GEOS-Chem chemistry.
!         CALL Convert_Units( IFLAG     = 2,                    &
!                             N_TRACERS = Input_Opt%N_TRACERS,  &
!                             TCVV      = Input_Opt%TCVV,       &
!                             AD        = State_Met%AD,         &
!                             STT       = State_Chm%Tracers    )
!-----------------------------------------------------------------------------

         ! Call the GEOS-Chem chemistry routines
         CALL Do_Chemistry ( am_I_Root = am_I_Root,            &
                             NI        = NI,                   &
                             NJ        = NJ,                   &
                             NL        = NL,                   &
                             Input_Opt = Input_Opt,            &
                             State_Chm = State_Chm,            &
                             State_Met = State_Met,            &    
                             RC        = RC                   )

!-----------------------------------------------------------------------------
! Prior to 12/4/12:
! Now convert units in the higher level routine GIGC_Chunk_Run (bmy, 12/4/12)
!         ! Convert the tracer concentrations in State_Chm%TRACERS back to
!         ! [v/v] so that they can be stored in the ESMF internal state
!         ! for the next chemistry timestep.
!         CALL Convert_Units( IFLAG     = 1,                    &
!                             N_TRACERS = Input_Opt%N_TRACERS,  &
!                             TCVV      = Input_Opt%TCVV,       &
!                             AD        = State_Met%AD,         &
!                             STT       = State_Chm%Tracers    )
!-----------------------------------------------------------------------------

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
      State_Chm%Species = CSPEC_FULL

      !### Debug
      !IF ( am_I_Root ) THEN
      !   WRITE(6,*) '##### GIGC_Do_Chem, TRC_OX after chem'
      !   WRITE(6,*) State_Chm%Tracers(1,1,:,2)
      !ENDIF

    END SUBROUTINE GIGC_Do_Chem
!EOC
END MODULE GIGC_ChemDr
#endif
