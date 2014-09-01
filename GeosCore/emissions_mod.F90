!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emissions_mod.F90
!
! !DESCRIPTION: Module emissions\_mod.F90 is a wrapper module to interface
! GEOS-Chem and HEMCO. It basically just calls the GEOS-Chem - HEMCO interface
! routines. For some specialty sims, a few additional steps are required that
! are also executed here.
!\\
!\\
! !INTERFACE:
!
MODULE EMISSIONS_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: EMISSIONS_INIT
  PUBLIC :: EMISSIONS_RUN
  PUBLIC :: EMISSIONS_FINAL
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_INIT
!
! !DESCRIPTION: Subroutine EMISSIONS\_INIT calls the HEMCO - GEOS-Chem
! interface initialization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_INIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_INIT
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
    IF ( RC/=GIGC_SUCCESS ) RETURN 

  END SUBROUTINE EMISSIONS_INIT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_RUN
!
! !DESCRIPTION: Subroutine EMISSIONS\_RUN calls the HEMCO - GEOS-Chem
! interface run routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_RUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_RUN
    USE CO2_MOD,            ONLY : EMISSCO2
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Run HEMCO
    CALL HCOI_GC_RUN( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
    IF ( RC/=GIGC_SUCCESS ) RETURN 

    ! For CO2 simulation, emissions are not added to Trac_Tend and hence
    ! not passed to the Tracers array during PBL mixing. Thus, need to add 
    ! emissions explicitly to the Tracers array here.
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
       CALL EMISSCO2 ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( RC/=GIGC_SUCCESS ) RETURN 
    ENDIF

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE EMISSIONS_RUN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EMISSIONS_FINAL
!
! !DESCRIPTION: Subroutine EMISSIONS\_FINAL calls the HEMCO - GEOS-Chem
! interface finalization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSIONS_FINAL( am_I_Root )
!
! !USES:
!
    USE HCOI_GC_MAIN_MOD, ONLY : HCOI_GC_FINAL
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    CALL HCOI_GC_Final( am_I_Root )

  END SUBROUTINE EMISSIONS_FINAL
!EOC
END MODULE EMISSIONS_MOD
