!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_diagn_mod.F90
!
! !DESCRIPTION: Module HCOI\_GC\_Diagn\_Mod.F90 is the GEOS-Chem interface 
! module for the HEMCO diagnostics. For every GEOS-Chem emissions diagnostics,
! a corresponding HEMCO diagnostics is created. The HEMCO diagnostics become
! (automatically) filled and updated when calling HEMCO. They are passed
! back to GEOS-Chem when writing the diagnostics (e.g. in diag3.F).
!\\
!\\
! Notes:
! \begin{itemize}
! \item The category specific diagnostics (anthropogenic, aircraft, etc.)
!  explicitly assume certain category numbers in the HEMCO configuration 
!  file (e.g. Cat=1 for anthropogenic, Cat=20 for aircraft, etc.).
!  Diagnostics will not represent what they should if these category numbers
!  get changed!
! \item In HEMCO, ocean sinks are treated as drydep and the calculated 
!  deposition velocities are passed to drydep\_mod.F. Hence, no Acetone or ALD2
!  ocean sink is calculated by HEMCO and the DMS diagnostics only includes
!  the ocean flux (this is NOT the net flux!!). 
!  If needed, we can build a simple wrapper in hcoi\_gc\_main\_mod.F90 that
!  explicitly calculates oceanic fluxes.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCOI_GC_Diagn_Mod
!
! !USES:
!
  ! GEOS-Chem diagnostic switches and arrays
  USE CMN_SIZE_Mod
#if defined( BPCH_DIAG )
  USE CMN_DIAG_Mod
  USE DIAG_Mod
  USE DIAG53_Mod
  USE DIAG56_Mod
#endif
  USE HCO_Diagn_Mod
  USE HCO_Error_Mod
  USE HCO_Interface_Mod

  IMPLICIT NONE
  PRIVATE

  ! Get parameters that define the different categories
#include "hcoi_gc_diagn_include.H"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOI_GC_Diagn_Init
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REMARKS:
!  This is currently a "bridge" module to provide backwards compatibility
!  with existing GEOS-Chem diagnostics.  We will eventually write all
!  diagnostics to netCDF format, but we are not quite there yet.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version. 
!  See the Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Diagn_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Diagn\_Init initializes the HEMCO 
! diagnostics in GEOS-Chem. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Diagn_Init( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  The category numbers must correspond to those in the HEMCO_Config.rc file.
!  We will have to come up with a better way of making sure that these
!  are consistent in the future.

!  CO emissions (ND29) 
!  --> Anthropogenic, biogenic, and biomass emissions are 
!      all covered in the respective sections. 
!  --> CO produced from methanol doesn't seem to be written anymore?!
!      Not filled for now.
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  See the Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: YesOrNo
    INTEGER            :: I, J,  HcoID, N,    AS
    INTEGER            :: ExtNr, Cat, Hier
    CHARACTER(LEN=31)  :: SpcName, DiagnName, Unit
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
 
    !=======================================================================
    ! HCOI_GC_DIAGN_INIT begins here!
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
       ' -> at HCOI_GC_Diagn_Init (in module GeosCore/hcoi_gc_diagn_mod.F90)'

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Emissions for CH4 specialty simulations are passed to  %%%%
    !%%%%  global_ch4_mod.F90 via HEMCO diagnostics, and not directly   %%%%
    !%%%%  from the HEMCO state pointer.  Therefore, we need to make    %%%%
    !%%%%  sure that routine DIAGN_CH4 is outside the BPCH_DIAG #if     %%%%
    !%%%%  block.  -- Bob Yantosca (25 Jan 2018)                        %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_CH4( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_CH4"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Emissions for Hg specialty simulations are passed to   %%%%
    !%%%%  mercury_mod.F90 via HEMCO diagnostics, and not directly      %%%%
    !%%%%  from the HEMCO state pointer.  Therefore, we need to make    %%%%
    !%%%%  sure that routine DIAGN_Hg is outside the BPCH_DIAG #if      %%%%
    !%%%%  block.  -- Bob Yantosca (25 Jan 2018)                        %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_Hg( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_Hg"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%  NOTE: Some diagnostics for the POPs simulation do not have   %%%%
    !%%%%  any species associated with them, and thus need to be        %%%%
    !%%%%  declared as manual diagnostics.  We need to move the call    %%%%
    !%%%%  to Diagn_POPs outside of the #ifdef BPCH_DIAG block.         %%%%
    !%%%%    -- Bob Yantosca (09 Oct 2018)                              %%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL Diagn_POPs( am_I_Root, Input_Opt, HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diagn_POPs"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#ifdef BPCH_DIAG

    !=======================================================================
    ! Define manual diagnostics
    !=======================================================================
    CALL Diagn_Dust    ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_LFlash  ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

#ifdef TOMAS
    CALL Diagn_TOMAS   ( am_I_Root, Input_Opt, HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
#endif
#endif
    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE HCOI_GC_Diagn_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_CH4
!
! !DESCRIPTION: Subroutine Diagn\_CH4 initializes diagnostics for the
!  CH4 simulation (ND58).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_CH4( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
!  Biomass diagnostics are defined in routine Diagn_Biomass.
!\\
!\\
!  CH4 diagnostics need to be defined even if ND58 is turned off because
!  the diagnostics are also being used to write CH4 emissions from the 
!  individual sources (gas, coal, etc.) into STT (in global\_ch4\_mod.F).
!  The categories defined here need to match the ones specified in the 
!  HEMCO configuration file.
!
! !REVISION HISTORY: 
!  13 Sep 2014 - C. Keller   - Initial version
!  16 Jun 2016 - C. Miller   - Now define species ID's with Ind_ funciton
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, id_CH4, Cat, HcoID, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_CH4 (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND58 diagnostics (CH4 emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the CH4 simulation is not selected
    IF ( .NOT. ( Input_Opt%ITS_A_CH4_SIM .OR. Ind_('CH4','A') > 0 ) ) RETURN

    ! Get default HEMCO species ID for CH4 
    id_CH4 = HCO_GetHcoID( 'CH4', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    !-----------------------------------------------------------------
    ! %%%%% CH4 from oil (Category 1 or species CH4_OIL)  %%%%%
    !-----------------------------------------------------------------

    ! Check if tagged CH4 simulation
    ! Otherwise, use CH4 category 1 emissions
    Cat   = 1
    HcoID = HCO_GetHcoID( 'CH4_OIL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_OIL'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! %%%%% CH4 from natural gas (Category 2 or species CH4_GAS)  %%%%%
    !-----------------------------------------------------------------

    ! Check if tagged CH4 simulation
    ! Otherwise, use CH4 category 1 emissions
    Cat   = 2
    HcoID = HCO_GetHcoID( 'CH4_GAS', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_GAS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !----------------------------------------------------------
    ! %%%%% CH4 from coal (Category 3 or species CH4_COL)  %%%%%
    !----------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 3
    HcoID = HCO_GetHcoID( 'CH4_COL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_COAL'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !---------------------------------------------------------------
    ! %%%%% CH4 from livestock (Category 4 or species CH4_LIV)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 4
    HcoID = HCO_GetHcoID( 'CH4_LIV', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LIVESTOCK'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from landfills (Category 5 or species CH4_LDF)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 5
    HcoID = HCO_GetHcoID( 'CH4_LDF', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LANDFILLS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF  

    !---------------------------------------------------------------
    ! %%%%% CH4 from wastewater (Category 6 or species CH4_WST)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 6
    HcoID = HCO_GetHcoID( 'CH4_WST', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_WASTEWATER'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF   
 
    !---------------------------------------------------------------
    ! %%%%% CH4 from rice (Category 7 or species CH4_RIC)  %%%%%
    !---------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 7
    HcoID = HCO_GetHcoID( 'CH4_RIC', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_RICE'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from other anth. sources (Category 8 or species CH4_OTA)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 8
    HcoID = HCO_GetHcoID( 'CH4_OTA', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_ANTHROTHER'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from biomass burning (Category 9 or species CH4_BBN)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 9
    HcoID = HCO_GetHcoID( 'CH4_BBN', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_BIOMASS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from wetlands (Category 10 or species CH4_WTL)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 10
    HcoID = HCO_GetHcoID( 'CH4_WTL', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_WETLAND'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from seeps (Category 11 or species CH4_SEE)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 11
    HcoID = HCO_GetHcoID( 'CH4_SEE', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_SEEPS'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from lakes (Category 12 or species CH4_LAK)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 12
    HcoID = HCO_GetHcoID( 'CH4_LAK', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_LAKES'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from termites (Category 13 or species CH4_TER)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 13
    HcoID = HCO_GetHcoID( 'CH4_TER', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_TERMITES'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-------------------------------------------------------------------------
    ! %%%%% CH4 from soil absorption (Category 14 or species CH4_SAB)  %%%%%
    !-------------------------------------------------------------------------

    ! Check if tagged CH4 simulation
    Cat   = 14
    HcoID = HCO_GetHcoID( 'CH4_SAB', HcoState )
    IF ( HcoID <= 0 ) THEN
       HcoID = id_CH4
    ENDIF

    IF ( HcoID > 0 ) THEN 

       ! Create diagnostic container
       DiagnName = 'CH4_SOILABSORB'
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

  END SUBROUTINE Diagn_CH4
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Dust
!
! !DESCRIPTION: Subroutine Diagn\_Dust initializes diagnostics for the
!  mineral dust aerosols (ND06).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Dust( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use dust
!  30 Sep 2014 - R. Yantosca - Update for TOMAS dust species
!  25 Oct 2016 - R. Yantosca - Make sure to cast INTEGER to LOGICAL values
!                              before comparing them in an IF statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: Is_DustDead
    LOGICAL            :: Is_DustGinoux
    INTEGER            :: ExtNr, Cat, HcoID, I, N
    CHARACTER(LEN=1)   :: ISTR1
    CHARACTER(LEN=2)   :: ISTR2
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_DUST (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_DUST begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#ifdef BPCH_DIAG

    ! Exit if we are doing a specialty simulation w/o dust
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       RETURN
    ENDIF

    ! Now use local LOGICAL variables to save ExtState%DustDead and
    ! ExtState%DustGinoux.  This will make sure these variables are 
    ! cast to LOGICAL, so that we can compare them in the same IF
    ! statement.  Otherwise GNU Fortran will choke. (bmy, 10/25/16)
    Is_DustDead   = ( ExtState%DustDead   > 0 )
    Is_DustGinoux = ( ExtState%DustGinoux > 0 )

    ! Define diagnostics if dust is used
    IF ( Input_Opt%ND06 > 0 ) THEN

       ! Get Ext. Nr
       IF ( Is_DustDead ) THEN
          ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustDead' )
          Cat   = -1
       ELSEIF ( Is_DustGinoux ) THEN
          ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustGinoux' )
          Cat   = -1
       ELSE
          ! Use offline dust emissions
          ExtNr = 0
          Cat   = CATEGORY_NATURAL
       ENDIF

       ! Do for each dust bin
       DO I = 1, NDSTBIN

#ifdef TOMAS

          ! Get species name (i.e. DUST1 .. DUST40) for TOMAS simulatiosn
          IF ( I < 10 )  THEN
             WRITE( ISTR1,'(i1)' ) I
             SpcName = 'DUST'   // ISTR1
          ELSE
             WRITE( ISTR2,'(i2)' ) I
             SpcName = 'DUST'   // ISTR2
          ENDIF
#else

          ! Get species name (i.e. DST1 .. DST4) for non TOMAS simualtions
          WRITE( ISTR1,'(i1)' ) I
          SpcName   = 'DST'   // ISTR1

#endif

          DiagnName = 'AD06_' // TRIM( SpcName )         

          ! HEMCO species ID 
          HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg',              &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                             AutoFill  = 1,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

       ENDDO 

       ! Add diagnostics for dust alkalinity
       IF ( Input_Opt%LDSTUP ) THEN

          ! Get Ext. Nr of used extension
          ExtNr = GetExtNr( HcoState%Config%ExtList, 'DustAlk' )
          IF ( ExtNr <= 0 ) THEN
             CALL HCO_Error( 'Cannot find dust alk extension', RC, &
                              THISLOC=LOC )
             RETURN      
          ENDIF

          ! Do for each dust bin
          DO I = 1, NDSTBIN

             ! Get species name (i.e. DSTAL1 .. DSTAL4)
             WRITE( ISTR1,'(i1)' ) I
             SpcName   = 'DSTAL' // ISTR1
             DiagnName = 'AD06_' // TRIM( SpcName )

             ! HEMCO species ID 
             HcoID = GetHemcoId( TRIM( SpcName ), HcoState, LOC, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Create diagnostic container
             CALL Diagn_Create( am_I_Root,                     & 
                                HcoState  = HcoState,          &
                                cName     = TRIM( DiagnName ), &
                                ExtNr     = ExtNr,             &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = HcoID,             &
                                SpaceDim  = 2,                 &
                                LevIDx    = -1,                &
                                OutUnit   = 'kg',              &
                                COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                                AutoFill  = 1,                 &
                                RC        = RC                  ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDDO
       ENDIF

    ENDIF
#endif

  END SUBROUTINE Diagn_Dust

!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_LFlash
!
! !DESCRIPTION: Subroutine Diagn\_LFlash initializes diagnostics for the
!  lightning flash diagnostics (ND56).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_LFlash( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Lightning flashes (ND56)
!  ==> Manually set in hcox_lightning_mod.F90
!  ==> The original code diagnoses the 'cumulative mean', i.e. it
!      sums all quantities and simply divides by the number of 
!      met. time steps when writing the diagnostics. Here, we only 
!      use the average since the last writeout. To simulate the
!      behavior of the original ND56 diagnostics, set 'OutOper' to
!      'Cumsum' (cumulative sum) and manually divide by the # of
!      timesteps when writing data to the output file!
!
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!  21 Aug 2014 - R. Yantosca - Exit for simulations that don't use lightning=
!  19 Feb 2015 - C. Keller   - Added entry for MTYPE
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, I, N
    CHARACTER(LEN=1)   :: ISTR
    CHARACTER(LEN=15)  :: SpcName
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_LFLASH (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! DIAGN_LFLASH begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if we are doing a specialty simulation w/o lightning
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    ! Define diagnostics
    IF ( Input_Opt%ND56 > 0 ) THEN

       ! Extension number
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'LightNOx')

       ! Exit if LightNOx is not turned on - the lightning NOx extension
       ! must be enabled in the HEMCO configuration file.
       IF ( ExtNr <= 0 ) THEN
          MSG = 'Lightning NOx is not enabled - cannot write diagnostics ND56!'
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! HEMCO species ID
       HcoID = GetHemcoId( 'NO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Loop over lighthing flash quantities
       DO I = 1, 3

          ! Pick the proper diagnostic name
          SELECT CASE( I )
             CASE( 1 )
                DiagnName = 'LIGHTNING_TOTAL_FLASHRATE'
             CASE( 2 )
                DiagnName = 'LIGHTNING_INTRACLOUD_FLASHRATE'
             CASE( 3 )
                DiagnName = 'LIGHTNING_CLOUDGROUND_FLASHRATE'
          END SELECT

          ! Define diagnostics ID
          N = 56000 + I

          ! Create diagnostic container
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             cID       = N,                 &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = HcoID,             &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'flashes/min/km2', &
                             OutOper   = 'Mean',            &
                             COL       = HcoState%Diagn%HcoDiagnIDDefault, &
                             AutoFill  = 0,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
 
       ! ---------------------------------------------------------- 
       ! Diagnostics for convective cloud top height.
       ! ---------------------------------------------------------- 

       ! Define diagnostics name and ID
       DiagnName = 'LIGHTNING_CLOUD_TOP'
       N         = 56004 

       ! Create diagnostic container
       CALL Diagn_Create( am_I_Root,                     & 
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          cID       = N,                 & 
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = '1',               &
                          OutOper   = 'Mean',            &
                          COL       = HcoState%Diagn%HcoDiagnIDDefault, &
                          AutoFill  = 0,                 &
                          RC        = RC                  ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF ! ND56 
#endif

  END SUBROUTINE Diagn_LFlash
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_POPs
!
! !DESCRIPTION: Subroutine Diagn\_POPs initializes several HEMCO manual 
!  diagnostics for the POPs simulation.  These diagnostics are updated
!  in the HEMCO extensions module hcox\_gc\_POPs\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_POPs( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Split off code from HCOI_GC_Diagn_Init into smaller routines in order to
!  make the code more manageable.
!
! !REVISION HISTORY: 
!  26 Aug 2014 - M. Sulprizio- Initial version
!  15 Oct 2018 - R. Yantosca - Changed "AD53_*" to "GCPOPS_*"
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, N, I
    CHARACTER(LEN=31)  :: DiagnName, OutOper, OutUnit
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_POPs (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define ND53 diagnostics (POPs emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the POPs simulation is not selected
    IF ( .not. Input_Opt%ITS_A_POPS_SIM ) RETURN

    ! For the HISTORY netCDF diagnostics, we want to get the instantaneous
    ! values archived by HEMCO and then let HISTORY do the averaging.
    OutOper = 'Instantaneous'

#if defined( BPCH_DIAG )
    ! Exit if ND53 diagnostics aren't turned on
    IF ( ND53 <= 0 ) RETURN

    ! For the bpch diagnostics, change units to kg/s to help in validating
    ! the netCDF diagnostics.  But select "Mean" diagnostics since these
    ! are only ouptut at the end.
    OutOper = 'Mean'
#endif

    ! Define diagnostics
    IF ( ExtState%GC_POPs > 0 ) THEN

       ! HEMCO extension # for POPs
       ExtNr = GetExtNr( HcoState%Config%ExtList, 'GC_POPs' )
       IF ( ExtNr <= 0 ) THEN
          MSG = 'Cannot find the POPs extension for HEMCO!'
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN      
       ENDIF

       !-------------------------------------------
       ! %%%%% Gas-phase POP emissions %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPG', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPG_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% OC-phase POP emissions%%%%%
       !-------------------------------------------
 
       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPOCPO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPPOCPO_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM( DiagnName ),                     &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &  
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% BC-phase POP emissions %%%%%
       !-------------------------------------------

       ! HEMCO species ID
       HcoID = GetHemcoId( 'POPPBCPO', HcoState, LOC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Create diagnostic container
       DiagnName = 'GCPOPS_POPPBCPO_SOURCE'
       CALL Diagn_Create( am_I_Root,                                         &
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = -1,                                    &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          OutOper   = OutOper,                               &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Could not define POPs diagnostic '// TRIM( DiagnName )
          CALL HCO_Error( MSG, RC, THISLOC=LOC )
          RETURN 
       ENDIF

       !-------------------------------------------
       ! %%%%% Manual diagnostics %%%%%
       !-------------------------------------------

       DO I = 1,12

          ! Define diagnostic names. These have to match the names
          ! in module HEMCO/Extensions/hcox_gc_POPs_mod.F90.
          IF ( I == 1 ) THEN
             DiagnName = 'GCPOPS_POPG_SOIL'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 2  ) THEN
             DiagnName = 'GCPOPS_POPG_LAKE'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 3  ) THEN
             DiagnName = 'GCPOPS_POPG_LEAF'
             OutUnit   = 'kg/m2/s'
          ELSEIF ( I == 4  ) THEN
             DiagnName = 'GCPOPS_SOIL2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 5  ) THEN
             DiagnName = 'GCPOPS_AIR2SOIL'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 6  ) THEN
             DiagnName = 'GCPOPS_LAKE2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 7  ) THEN
             DiagnName = 'GCPOPS_AIR2LAKE'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 8  ) THEN
             DiagnName = 'GCPOPS_LEAF2AIR'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 9  ) THEN
             DiagnName = 'GCPOPS_AIR2LEAF'
             OutUnit   = 'ng/m2/day'
          ELSEIF ( I == 10 ) THEN
             DiagnName = 'GCPOPS_SOILAIR_FUG'
             OutUnit   = '1'
          ELSEIF ( I == 11 ) THEN
             DiagnName = 'GCPOPS_LAKEAIR_FUG'
             OutUnit   = '1'
          ELSEIF ( I == 12 ) THEN
             DiagnName = 'GCPOPS_LEAFAIR_FUG'
             OutUnit   = '1'
          ENDIF

          ! Create manual diagnostics
          CALL Diagn_Create( am_I_Root,                                      & 
                             HcoState  = HcoState,                           &
                             cName     = TRIM( DiagnName ),                  &
                             ExtNr     = ExtNr,                              &
                             Cat       = -1,                                 &
                             Hier      = -1,                                 &
                             HcoID     = -1,                                 &
                             SpaceDim  = 2,                                  &
                             LevIDx    = -1,                                 &
                             OutUnit   = OutUnit,                            &
                             OutOper   = OutOper,                            &
                             COL       = HcoState%Diagn%HcoDiagnIDManual,    &
                             AutoFill  = 1,                                  &
                             RC        = RC                                 ) 

          ! Trap potential errors
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Could not define POPs diagnostic: '// TRIM( DiagnName )
             CALL HCO_Error( Msg, RC, THISLOC=LOC )
             RETURN 
          ENDIF
       ENDDO

    ENDIF

  END SUBROUTINE Diagn_POPs
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Hg
!
! !DESCRIPTION: Subroutine Diagn\_Hg initializes diagnostics for the
!  Hg specialty simulation. For now, this is just a placeholder. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Hg( am_I_Root, Input_Opt, HcoState, ExtState, RC ) 
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  23 Sep 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, HcoID, Cat, N
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_Hg (hcoi_gc_diagn_mod.F90)'

    !=======================================================================
    ! Define diagnostics (POPs emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Exit if the mercury simulation is not selected
    IF ( .not. Input_Opt%ITS_A_MERCURY_SIM ) RETURN

    ! HEMCO species ID
    HcoID = GetHemcoId( 'Hg0', HcoState, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Extension number is always zero (HEMCO core)
    ExtNr = 0

    !-------------------------------------------
    ! %%%%% ARTISANAL HG (Hg0) %%%%%
    !-------------------------------------------
 
    ! Create diagnostic container
    DiagnName = 'HG0_ARTISANAL'
    Cat       = 8
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-------------------------------------------
    ! %%%%% NATURAL HG (Hg0) %%%%%
    !-------------------------------------------
 
    ! Create diagnostic container
    DiagnName = 'HG0_NATURAL'
    Cat       = CATEGORY_NATURAL
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       )
    IF ( RC /= HCO_SUCCESS ) RETURN


    !---------------------------------------------
    ! %%%%% ANTHROPOGENIC HG (Hg0, Hg2, HgP) %%%%%
    !---------------------------------------------
 
    Cat = CATEGORY_ANTHRO

    ! Create diagnostic container
    DiagnName = 'HG0_ANTHRO'
    CALL Diagn_Create( am_I_Root,                                            & 
                       HcoState  = HcoState,                                 &
                       cName     = TRIM(DiagnName),                          &
                       ExtNr     = ExtNr,                                    &
                       Cat       = Cat,                                      &
                       Hier      = -1,                                       &
                       HcoID     = HcoID,                                    &
                       SpaceDim  = 2,                                        &
                       LevIDx    = -1,                                       &
                       OutUnit   = 'kg/m2/s',                                &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,          &
                       AutoFill  = 1,                                        &
                       RC        = RC                                       ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Hg2
    HcoID = GetHemcoId( 'Hg2', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HG2_ANTHRO'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! HgP
    HcoID = GetHemcoId( 'HgP', HcoState, LOC, RC, ERR=.FALSE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create diagnostic container
    IF ( HcoID > 0 ) THEN
       DiagnName = 'HGP_ANTHRO'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                     ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !-------------------------------------------
    ! %%%%% BIOMASS BURNING HG %%%%%
    ! ==> defined in Diagn_Biomass
    !-------------------------------------------

    ! Do only if Hg0 is defined ... 
    HcoID = HCO_GetHcoID( 'Hg0', HcoState )
    IF ( HcoID > 0 ) THEN

        Cat   = -1
        ExtNr = GetExtNr( HcoState%Config%ExtList, 'GFED' )
        IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'FINN' )
        IF ( ExtNr <= 0 ) THEN
           ExtNr = 0
           Cat   = CATEGORY_BIOMASS
        ENDIF


       ! Create diagnostic container
       DiagnName = 'BIOMASS_HG0'
       CALL Diagn_Create( am_I_Root,                                         & 
                          HcoState  = HcoState,                              &
                          cName     = TRIM(DiagnName),                       &
                          ExtNr     = ExtNr,                                 &
                          Cat       = Cat,                                   &
                          Hier      = -1,                                    &
                          HcoID     = HcoID,                                 &
                          SpaceDim  = 2,                                     &
                          LevIDx    = -1,                                    &
                          OutUnit   = 'kg/m2/s',                             &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,       &
                          AutoFill  = 1,                                     &
                          RC        = RC                                    ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Diagn_Hg
!EOC
#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_TOMAS
!
! !DESCRIPTION: This creates diagnostics for bulk emissions that will be called
! to scale into TOMAS bins. May not even be necessary. (JKodros 6/2/15)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_TOMAS( am_I_Root, Input_Opt, HcoState, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_State_Mod,      ONLY : HCO_State
    USE HCO_State_Mod,      ONLY : HCO_GetHcoID
    USE HCOX_State_Mod,     ONLY : Ext_State
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(EXT_State),  POINTER        :: ExtState   ! Extensions state object
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success

!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ExtNr, Cat, HcoID, N
    INTEGER            :: id_BCPI, id_BCPO, id_OCPI, id_OCPO
    INTEGER            :: IDSO4
    INTEGER            :: IDCO
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGN_TOMAS (hcoi_gc_diagn_mod.F90)'
    !=======================================================================
    ! Define ND?? diagnostics (TOMAS-related emissions)
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

#if defined( BPCH_DIAG )

    ! Exit if the CH4 simulation is not selected
    !IF ( .NOT. ( Input_Opt%ITS_A_CH4_SIM .OR. id_CH4 > 0 ) ) RETURN
    ! SOME SORT OF IF DEFINED TOMAS HERE

    ! Get default HEMCO species ID for BC/OC
    id_BCPI = HCO_GetHcoID( 'BCPI', HcoState )
    id_BCPO = HCO_GetHcoID( 'BCPO', HcoState )
    id_OCPI = HCO_GetHcoID( 'OCPI', HcoState )
    id_OCPO = HCO_GetHcoID( 'OCPO', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0
    !-----------------------------------------------------------------
    ! %%%%% BCPI from anthro (Category 1 or species BCPI_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_BCPI
       ! Create diagnostic container
       DiagnName = 'BCPI_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BCPO from anthro (Category 1 or species BCPO_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_BCPO
       ! Create diagnostic container
       DiagnName = 'BCPO_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPI from anthro (Category 1 or species OCPI_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_OCPI
       ! Create diagnostic container
       DiagnName = 'OCPI_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPO from anthro (Category 1 or species OCPO_ANTH)  %%%%%
    !-----------------------------------------------------------------

       Cat = CATEGORY_ANTHRO
       HcoID = id_OCPO
       ! Create diagnostic container
       DiagnName = 'OCPO_ANTH'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 3,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ! ------------ NOW DEAL WITH BIOMASS BURNING --------------------
    ! First test if GFED is used.  If not, then test if FINN is used.
    ! If not, then use extension # 0 and the default biomass category.
      Cat   = -1
      ExtNr = GetExtNr( HcoState%Config%ExtList, 'GFED' )
      IF ( ExtNr <= 0 ) ExtNr = GetExtNr( HcoState%Config%ExtList, 'FINN' )
      IF ( ExtNr <= 0 ) THEN
         ExtNr = 0
         Cat   = CATEGORY_BIOMASS
      ENDIF

    !-----------------------------------------------------------------
    ! %%%%% BPCI from BIOB (Category ? or species BCPI_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_BCPI
       ! Create diagnostic container
       DiagnName = 'BCPI_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% BPCO from BIOB (Category ? or species BCPO_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_BCPO
       ! Create diagnostic container
       DiagnName = 'BCPO_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPI from BIOB (Category ? or species OCPI_bb)  %%%%%
    !-----------------------------------------------------------------

       HcoID = id_OCPI
       ! Create diagnostic container
       DiagnName = 'OCPI_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% OCPO from BIOB (Category ? or species OCPI_bb)  %%%%%
    !-----------------------------------------------------------------


       HcoID = id_OCPO
       ! Create diagnostic container
       DiagnName = 'OCPO_BB'
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          ExtNr     = ExtNr,             &
                          Cat       = Cat,               &
                          Hier      = -1,                &
                          HcoID     = HcoID,             &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/m2/s',         &
                          COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                          AutoFill  = 1,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% SO4 from ANTRHO (Category ? or species SO4_ANTH)  %%%%%
    !-----------------------------------------------------------------
    ! ID for SO4
    IDSO4 = HCO_GetHcoID( 'SO4', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    Cat = CATEGORY_ANTHRO
    HcoID = IDSO4
    ! Create diagnostic container
    DiagnName = 'SO4_ANTH'
    CALL Diagn_Create( am_I_Root,                     &
                       HcoState  = HcoState,          &
                       cName     = TRIM( DiagnName ), &
                       ExtNr     = ExtNr,             &
                       Cat       = Cat,               &
                       Hier      = -1,                &
                       HcoID     = HcoID,             &
                       SpaceDim  = 3,                 &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg/m2/s',         &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                       AutoFill  = 1,                 &
                       RC        = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! %%%%% CO from ANTRHO (Category ? or species CO_ANTH)  %%%%%
    !-----------------------------------------------------------------
    ! ID for CO
    IDCO = HCO_GetHcoID( 'CO', HcoState )

    ! Extension number is zero (HEMCO core) until defined otherwise
    ExtNr = 0

    Cat = CATEGORY_ANTHRO
    HcoID = IDCO
    ! Create diagnostic container
    DiagnName = 'CO_ANTH'
    CALL Diagn_Create( am_I_Root,                     &
                       HcoState  = HcoState,          &
                       cName     = TRIM( DiagnName ), &
                       ExtNr     = ExtNr,             &
                       Cat       = Cat,               &
                       Hier      = -1,                &
                       HcoID     = HcoID,             &
                       SpaceDim  = 3,                 &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg/m2/s',         &
                       COL       = HcoState%Diagn%HcoDiagnIDManual,  &
                       AutoFill  = 1,                 &
                       RC        = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN

#endif

  END SUBROUTINE Diagn_TOMAS
!EOC
#endif
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetHemcoId
!
! !DESCRIPTION: Function GetHemcoId returns the HEMCO species ID number 
!  corresponding to a given HEMCO species name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetHemcoId( HcoName, HcoState, Loc, RC, ERR ) RESULT( HcoID )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_State_Mod, ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS: 
!
    CHARACTER(LEN=*),  INTENT(IN)  :: HcoName    ! HEMCO species name
    TYPE(HCO_State),   POINTER     :: HcoState   ! HEMCO State object
    CHARACTER(LEN=*),  INTENT(IN)  :: Loc        ! Calling routine
    LOGICAL, OPTIONAL, INTENT(IN)  :: Err        ! Return error if not found
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC         ! Success or failure?
!
! !RETURN VALUE:
!
    INTEGER                        :: HcoID      ! HEMCO species ID #
!
! !REMARKS:
!  This is a wrapper function to simplify the code above.   Calls to 
!  HCO_GetHcoId and HCO_Error are made from here. 
! 
! !REVISION HISTORY: 
!  20 Aug 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG
    LOGICAL            :: ERROR = .TRUE.

    !=======================================================================
    ! GetHemcoId begins here!
    !=======================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Prompt error?
    IF ( PRESENT(ERR) ) ERROR = ERR

    ! Get the HEMCO species ID from the name
    HcoID = HCO_GetHcoID( HcoName, HcoState )

    ! Exit with error if the species is not valid
    ! (HCO_Error will set RC = HCO_FAIL)
    IF ( HcoID <= 0 .AND. ERROR ) THEN
       MSG = 'This is not a HEMCO species: ' // HcoName
       CALL HCO_Error( MSG, RC, THISLOC=Loc )
    ENDIF

  END FUNCTION GetHemcoId
!EOC
END MODULE HCOI_GC_Diagn_Mod
