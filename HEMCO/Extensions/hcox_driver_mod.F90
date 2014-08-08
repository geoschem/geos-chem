!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_driver_mod.F90 
!
! !DESCRIPTION: Module hcox\_driver\_mod.F90 contains the driver routines 
! (INIT, RUN, FINAL) for the HEMCO extensions. It determines the extensions 
! to be used (based on the settings specified in the configuration file) 
! and invokes the respective extension module calls.
!\\
!\\
! Call this module at the HEMCO - model interface level to execute the
! HEMCO extensions.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_Driver_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State 

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Init
  PUBLIC :: HCOX_Run
  PUBLIC :: HCOX_Final
!
! !REMARKS: 
! (1) The extension option objects (e.g. meteorological variables) are 
!     defined in the HEMCO - model interface module and passed to this 
!     module. 
! (2) To add/remove HEMCO extensions from a model application, just 
!     add/remove the corresponding initialize, run, and finalize calls 
!     in the respective driver routines!
!
! !REVISION HISTORY:
!  15 Dec 2013 - C. Keller   - Initial version 
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
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
! !IROUTINE: HCOX_Init
!
! !DESCRIPTION: Subroutine HCOX\_Init is the driver routine to initialize 
!  all enabled HEMCO extensions. 
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Init( amIRoot, HcoState, ExtState, RC )
!
! !USES:
!
    ! HCO routines
    USE HCOX_State_MOD,        ONLY : ExtStateInit
    USE HCOX_Custom_Mod,       ONLY : HCOX_Custom_Init
    USE HCOX_SeaFlux_Mod,      ONLY : HCOX_SeaFlux_Init
    USE HCOX_PARANOX_Mod,      ONLY : HCOX_PARANOX_Init
    USE HCOX_LightNox_Mod,     ONLY : HCOX_LightNox_Init
    USE HCOX_SoilNox_Mod,      ONLY : HCOX_SoilNox_Init
    USE HCOX_DustDead_Mod,     ONLY : HCOX_DustDead_Init
    USE HCOX_DustGinoux_Mod,   ONLY : HCOX_DustGinoux_Init
    USE HCOX_SeaSalt_Mod,      ONLY : HCOX_SeaSalt_Init
    USE HCOX_GFED3_Mod,        ONLY : HCOX_GFED3_Init
    USE HCOX_Megan_Mod,        ONLY : HCOX_Megan_Init
    USE HcoX_FINN_Mod,         ONLY : HcoX_FINN_Init
    USE HCOX_GC_RnPbBe_Mod,    ONLY : HCOX_GC_RnPbBe_Init
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: amIRoot    ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(Ext_State),  POINTER        :: ExtState   ! HEMCO extension object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  07 Jul 2014 - R. Yantosca - Now init GEOS-Chem Rn-Pb-Be emissions module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! HCOX_INIT begins here!
    !=================================================================

    ! Error handling 
    CALL HCO_ENTER ('HCOX_INIT (hcox_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=================================================================
    ! Initialize extensions 
    !=================================================================

    ! Initialize extension object
    CALL ExtStateInit ( ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Custom 
    !-----------------------------------------------------------------
    CALL HCOX_Custom_Init( amIRoot, HcoState, 'Custom', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! SeaFlux
    !-----------------------------------------------------------------
    CALL HCOX_SeaFlux_Init( amIRoot, HcoState, 'SeaFlux', ExtState, RC )
    IF ( RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! ParaNox
    !-----------------------------------------------------------------
    CALL HCOX_PARANOX_INIT( amIRoot, HcoState, 'ParaNOx', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! LightNox 
    !-----------------------------------------------------------------
    CALL HCOX_LightNox_Init ( amIRoot,  HcoState, 'LightNOx', &
                              ExtState, RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! SoilNox 
    !-----------------------------------------------------------------
    CALL HCOX_SoilNox_Init( amIRoot, HcoState, 'SoilNOx', ExtState, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Dust emissions (DEAD model) 
    !-----------------------------------------------------------------
    CALL HCOX_DustDead_Init( amIRoot, HcoState, 'DustDead', &
                             ExtState,  RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Dust Ginoux emissions 
    !-----------------------------------------------------------------
    CALL HCOX_DustGinoux_Init( amIRoot, HcoState, 'DustGinoux', &
                               ExtState,  RC                       )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! SeaSalt aerosol extension
    !-----------------------------------------------------------------
    CALL HCOX_SeaSalt_Init( amIRoot, HcoState, 'SeaSalt', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! MEGAN extension
    !-----------------------------------------------------------------
    CALL HCOX_Megan_Init( amIRoot, HcoState, 'MEGAN', ExtState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! GFED3 extension
    !-----------------------------------------------------------------
    CALL HCOX_GFED3_Init( amIRoot, HcoState, 'GFED3', ExtState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! FINN biomass burning emissions
    !-----------------------------------------------------------------
    CALL HcoX_FINN_Init( amIRoot, HcoState, 'FINN', ExtState, RC ) 
    IF( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Extension for GEOS-Chem Rn-Pb-Be specialty simulation
    !-----------------------------------------------------------------
    CALL HCOX_GC_RnPbBe_Init( amIRoot, HcoState, 'GC_Rn-Pb-Be', &
                                       ExtState,  RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Add extensions here ...
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Sanity checks 
    !-----------------------------------------------------------------
    IF ( ExtState%DustDead .AND. ExtState%DustGinoux ) THEN
       MSG = 'Ginoux and DEAD dust emissions switched on!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Leave w/ success
    !-----------------------------------------------------------------
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Run
!
! !DESCRIPTION: Subroutine HCOX\_Run is the driver routine to run the HEMCO
! extensions. All enabled emission extensions are executed, and the
! emissions calculated therein become added to the respective flux arrays 
! in HcoState.\\
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Run( amIRoot, HcoState, ExtState, RC )
!
! !USES:
!
    ! HCO extensions
    USE HCOX_Custom_Mod,     ONLY : HCOX_Custom_Run
    USE HCOX_SeaFlux_Mod,    ONLY : HCOX_SeaFlux_Run 
    USE HCOX_ParaNox_Mod,    ONLY : HCOX_ParaNox_Run 
    USE HCOX_LightNox_Mod,   ONLY : HCOX_LightNox_Run 
    USE HCOX_SoilNox_Mod,    ONLY : HCOX_SoilNox_Run 
    USE HCOX_DustDead_Mod,   ONLY : HCOX_DustDead_Run 
    USE HCOX_DustGinoux_Mod, ONLY : HCOX_DustGinoux_Run 
    USE HCOX_SeaSalt_Mod,    ONLY : HCOX_SeaSalt_Run 
    USE HCOX_Megan_Mod,      ONLY : HCOX_Megan_Run 
    USE HCOX_GFED3_Mod,      ONLY : HCOX_GFED3_Run 
    USE HcoX_FINN_Mod,       ONLY : HcoX_FINN_Run 
    USE HCOX_GC_RnPbBe_Mod,  ONLY : HCOX_GC_RnPbBe_Run
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: amIRoot    ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    TYPE(Ext_State),  POINTER        :: ExtState   ! Extension options object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !NOTE:
!
! The ExtState object contains all extension option objects. In particular, 
! it contains the pointers to all met fields used by the extensions. These
! pointers have to be set in the HEMCO-model interface module beforehand!
!
! !REVISION HISTORY: 
!  15 Dec 2013 - C. Keller   - Initial version 
!  07 Jul 2014 - R. Yantosca - Now Run GEOS-Chem Rn-Pb-Be emissions module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG

    !=======================================================================
    ! HCOX_RUN begins here!
    !=======================================================================

    ! For error handling
    CALL HCO_ENTER ('HCOX_RUN (hcox_driver_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------------
    ! Customized emissions 
    !-----------------------------------------------------------------------
    IF ( ExtState%Custom ) THEN
       CALL HCOX_Custom_Run( amIRoot, ExtState, HcoState, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! SeaFlx (Air-sea exchange)
    !-----------------------------------------------------------------------
    IF ( ExtState%SeaFlux) THEN
       CALL HCOX_SeaFlux_Run( amIRoot, ExtState, HcoState, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! ParaNox (Ship NO emissions) 
    !-----------------------------------------------------------------------
    IF (ExtState%ParaNOx ) THEN
       CALL HCOX_ParaNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Lightning NOx  
    !-----------------------------------------------------------------------
    IF ( ExtState%LightNOx ) THEN
       CALL HCOX_LightNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! SoilNOx  
    !-----------------------------------------------------------------------
    IF ( ExtState%SoilNOx ) THEN
       CALL HCOX_SoilNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Dust emissions (DEAD model) 
    !-----------------------------------------------------------------------
    IF ( ExtState%DustDead ) THEN
       CALL HCOX_DustDead_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Dust emissions (Ginoux)
    !-----------------------------------------------------------------------
    IF ( ExtState%DustGinoux ) THEN
       CALL HCOX_DustGinoux_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Sea salt aerosols
    !-----------------------------------------------------------------------
    IF ( ExtState%SeaSalt ) THEN
       CALL HCOX_SeaSalt_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! ----------------------------------------------------------------
    ! Megan biogenic emissions 
    ! ----------------------------------------------------------------
    IF ( ExtState%Megan ) THEN
       CALL HCOX_Megan_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! GFED3 biomass burning emissions 
    !-----------------------------------------------------------------
    IF ( ExtState%GFED3 ) THEN
       CALL HCOX_GFED3_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! ----------------------------------------------------------------
    ! FINN biomass burning emissions 
    ! ----------------------------------------------------------------
    IF ( ExtState%FINN ) THEN
       CALL HcoX_FINN_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! Emissions for GEOS-Chem Rn-Pb-Be specialty simulation
    !-----------------------------------------------------------------
    IF ( ExtState%GC_RnPbBe ) THEN
       CALL HCOX_GC_RnPbBe_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! Add extensions here ...
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Return w/ success 
    !-----------------------------------------------------------------
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_RUN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Final
!
! !DESCRIPTION: Subroutine HCOX\_Final finalizes all HEMCO extensions. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Final( ExtState )
!
! !USES:
!
    USE HCOX_State_Mod,      ONLY : ExtStateFinal
    USE HCOX_Custom_Mod,     ONLY : HCOX_Custom_Final
    USE HCOX_SeaFlux_Mod,    ONLY : HCOX_SeaFlux_Final
    USE HCOX_PARANOX_Mod,    ONLY : HCOX_PARANOX_Final
    USE HCOX_LightNox_Mod,   ONLY : HCOX_LightNox_Final
    USE HCOX_SoilNox_Mod,    ONLY : HCOX_SoilNox_Final
    USE HCOX_DustDead_Mod,   ONLY : HCOX_DustDead_Final
    USE HCOX_DustGinoux_Mod, ONLY : HCOX_DustGinoux_Final
    USE HCOX_SeaSalt_Mod,    ONLY : HCOX_SeaSalt_Final
    USE HCOX_Megan_Mod,      ONLY : HCOX_Megan_Final
    USE HCOX_GFED3_Mod,      ONLY : HCOX_GFED3_Final
    USE HcoX_FINN_Mod,       ONLY : HcoX_FINN_Final
    USE HCOX_GC_RnPbBe_Mod,  ONLY : HCOX_GC_RnPbBe_Final
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER :: ExtState     ! Extension options object 
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  07 Jul 2014 - R. Yantosca - Now finalize GEOS-Chem Rn-Pb-Be emissions pkg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    !=================================================================
    ! HCOX_FINAL begins here!
    !=================================================================

    IF ( ASSOCIATED( ExtState ) ) THEN 

       ! Nullify all ExtState object pointers
       CALL ExtStateFinal( ExtState ) 

       ! Call individual cleanup routines
       IF ( ExtState%Custom     ) CALL HCOX_Custom_Final()
       IF ( ExtState%SeaFlux    ) CALL HCOX_SeaFlux_Final()
       IF ( ExtState%ParaNOx    ) CALL HCOX_PARANOX_Final()
       IF ( ExtState%LightNOx   ) CALL HCOX_LIGHTNOX_Final()
       IF ( ExtState%DustDead   ) CALL HCOX_DustDead_Final()
       IF ( ExtState%DustGinoux ) CALL HCOX_DustGinoux_Final()
       IF ( ExtState%SeaSalt    ) CALL HCOX_SeaSalt_Final()
       IF ( ExtState%Megan      ) CALL HCOX_Megan_Final()
       IF ( ExtState%GFED3      ) CALL HCOX_GFED3_Final()
       IF ( ExtState%SoilNOx    ) CALL HCOX_SoilNox_Final()  
       IF ( ExtState%FINN       ) CALL HcoX_FINN_Final
       IF ( ExtState%GC_RnPbBe  ) CALL HCOX_GC_RnPbBe_Final()
       
       ! Deallocate ExtState object
       DEALLOCATE( ExtState )
       ExtState => NULL()
    ENDIF
 
  END SUBROUTINE HCOX_FINAL
!EOC
END MODULE HCOX_DRIVER_MOD
