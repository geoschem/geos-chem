!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_driver_mod.F90 
!
! !DESCRIPTION: Module HCOX\_DRIVER\_MOD.F90 contains the driver routines 
! (INIT, RUN, FINAL) for the HEMCO extensions. It determines the extensions 
! to be used (based on the settings specified in the configuration file) 
! and invokes the respective extension module calls.\\
! Call this module at the HEMCO - model interface level to execute the
! HEMCO extensions.
! \\
! !INTERFACE:
!
MODULE HCOX_DRIVER_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_STATE_MOD,   ONLY : HCO_State
  USE HCOX_State_MOD,  ONLY : Ext_State 

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_INIT
  PUBLIC :: HCOX_RUN
  PUBLIC :: HCOX_FINAL
!
! !NOTES: 
!\\ 
! - The extension option objects (e.g. meteorological variables) are 
!   defined in the HEMCO - model interface module and passed to this 
!   module. 
!\\
! - To add/remove HEMCO extensions from a model application, just 
!   add/remove the corresponding initialize, run, and finalize calls in 
!   the respective driver routines!
!\\
!\\
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
! !IROUTINE: hcox_init
!
! !DESCRIPTION: Subroutine HCOX\_INIT is the driver routine to initialize 
!  all enabled HEMCO extensions. 
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_Init( amIRoot, HcoState, ExtState, RC )
!
! !USES:
!
    ! HCO routines
    USE HCOX_State_MOD,        ONLY : ExtStateInit
    USE HcoX_Custom_Mod,       ONLY : HcoX_Custom_Init
    USE HcoX_SeaFlux_Mod,      ONLY : HcoX_SeaFlux_Init
    USE HcoX_PARANOX_Mod,      ONLY : HcoX_PARANOX_Init
    USE HcoX_LightNox_Mod,     ONLY : HcoX_LightNox_Init
    USE HcoX_SoilNox_Mod,      ONLY : HcoX_SoilNox_Init
    USE HcoX_DustDead_Mod,     ONLY : HcoX_DustDead_Init
    USE HcoX_DustGinoux_Mod,   ONLY : HcoX_DustGinoux_Init
    USE HcoX_SeaSalt_Mod,      ONLY : HcoX_SeaSalt_Init
    USE HcoX_GFED3_Mod,        ONLY : HcoX_GFED3_Init
    USE HcoX_Megan_Mod,        ONLY : HcoX_Megan_Init
    USE HcoX_GC_RnPbBe_Mod,    ONLY : HcoX_GC_RnPbBe_Init
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
    ! Arguments needed to read soil NOx restart file
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
    CALL HcoX_Custom_Init( amIRoot, HcoState, 'Custom', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! SeaFlux
    !-----------------------------------------------------------------
    CALL HcoX_SeaFlux_Init( amIRoot, HcoState, 'SeaFlux', ExtState, RC )
    IF ( RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! ParaNox
    !-----------------------------------------------------------------
    CALL HCOX_PARANOX_INIT( amIRoot, HcoState, 'ParaNOx', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! LightNox 
    !-----------------------------------------------------------------
    CALL HcoX_LightNox_Init ( amIRoot,  HcoState, 'LightNOx', &
                              ExtState, RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! SoilNox 
    !-----------------------------------------------------------------
    CALL HcoX_SoilNox_Init( amIRoot, HcoState, 'SoilNOx', ExtState, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Dust emissions (DEAD model) 
    !-----------------------------------------------------------------
    CALL HcoX_DustDead_Init( amIRoot, HcoState, 'DustDead', &
                             ExtState,  RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Dust Ginoux emissions 
    !-----------------------------------------------------------------
    CALL HcoX_DustGinoux_Init( amIRoot, HcoState, 'DustGinoux', &
                               ExtState,  RC                       )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! SeaSalt aerosol extension
    !-----------------------------------------------------------------
    CALL HcoX_SeaSalt_Init( amIRoot, HcoState, 'SeaSalt', ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! MEGAN extension
    !-----------------------------------------------------------------
    CALL HcoX_Megan_Init( amIRoot, HcoState, 'MEGAN', ExtState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! GFED3 extension
    !-----------------------------------------------------------------
    CALL HcoX_GFED3_Init( amIRoot, HcoState, 'GFED3', ExtState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Extension for GEOS-Chem Rn-Pb-Be specialty simulation
    !-----------------------------------------------------------------
    CALL HcoX_GC_RnPbBe_Init( amIRoot, HcoState, 'GC_Rn-Pb-Be', &
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

  END SUBROUTINE HcoX_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_run
!
! !DESCRIPTION: Subroutine HCOX\_RUN is the driver routine to run the HEMCO
! extensions. All enabled emission extensions are executed, and the
! emissions calculated therein become added to the respective flux arrays 
! in HcoState.\\
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_RUN( amIRoot, HcoState, ExtState, RC )
!
! !USES:
!
    ! HCO extensions
    USE HcoX_Custom_Mod,     ONLY : HcoX_Custom_Run
    USE HcoX_SeaFlux_Mod,    ONLY : HcoX_SeaFlux_Run 
    USE HcoX_ParaNox_Mod,    ONLY : HcoX_ParaNox_Run 
    USE HcoX_LightNox_Mod,   ONLY : HcoX_LightNox_Run 
    USE HcoX_SoilNox_Mod,    ONLY : HcoX_SoilNox_Run 
    USE HcoX_DustDead_Mod,   ONLY : HcoX_DustDead_Run 
    USE HcoX_DustGinoux_Mod, ONLY : HcoX_DustGinoux_Run 
    USE HcoX_SeaSalt_Mod,    ONLY : HcoX_SeaSalt_Run 
    USE HcoX_Megan_Mod,      ONLY : HcoX_Megan_Run 
    USE HcoX_GFED3_Mod,      ONLY : HcoX_GFED3_Run 
    USE HcoX_GC_RnPbBe_Mod,  ONLY : HcoX_GC_RnPbBe_Run
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
       CALL HcoX_Custom_Run( amIRoot, ExtState, HcoState, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! SeaFlx (Air-sea exchange)
    !-----------------------------------------------------------------------
    IF ( ExtState%SeaFlux) THEN
       CALL HcoX_SeaFlux_Run( amIRoot, ExtState, HcoState, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! ParaNox (Ship NO emissions) 
    !-----------------------------------------------------------------------
    IF (ExtState%ParaNOx ) THEN
       CALL HcoX_ParaNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Lightning NOx  
    !-----------------------------------------------------------------------
    IF ( ExtState%LightNOx ) THEN
       CALL HcoX_LightNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! SoilNOx  
    !-----------------------------------------------------------------------
    IF ( ExtState%SoilNOx ) THEN
       CALL HcoX_SoilNox_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Dust emissions (DEAD model) 
    !-----------------------------------------------------------------------
    IF ( ExtState%DustDead ) THEN
       CALL HcoX_DustDead_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Dust emissions (Ginoux)
    !-----------------------------------------------------------------------
    IF ( ExtState%DustGinoux ) THEN
       CALL HcoX_DustGinoux_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------------
    ! Sea salt aerosols
    !-----------------------------------------------------------------------
    IF ( ExtState%SeaSalt ) THEN
       CALL HcoX_SeaSalt_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    ! ----------------------------------------------------------------
    ! Megan biogenic emissions 
    ! ----------------------------------------------------------------
    IF ( ExtState%Megan ) THEN
       CALL HcoX_Megan_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! GFED3 biomass burning emissions 
    !-----------------------------------------------------------------
    IF ( ExtState%GFED3 ) THEN
       CALL HcoX_GFED3_Run( amIRoot, ExtState, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
    ENDIF

    !-----------------------------------------------------------------
    ! Emissions for GEOS-Chem Rn-Pb-Be specialty simulation
    !-----------------------------------------------------------------
    IF ( ExtState%GC_RnPbBe ) THEN
       CALL HcoX_GC_RnPbBe_Run( amIRoot, ExtState, HcoState, RC )
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
! !IROUTINE: hcox_final
!
! !DESCRIPTION: Subroutine HCOX\_FINAL finalizes all HEMCO extensions. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_FINAL( ExtState )
!
! !USES:
!
    USE HcoX_State_Mod,      ONLY : ExtStateFinal
    USE HcoX_Custom_Mod,     ONLY : HcoX_Custom_Final
    USE HcoX_SeaFlux_Mod,    ONLY : HcoX_SeaFlux_Final
    USE HcoX_PARANOX_Mod,    ONLY : HcoX_PARANOX_Final
    USE HcoX_LightNox_Mod,   ONLY : HcoX_LightNox_Final
    USE HcoX_SoilNox_Mod,    ONLY : HcoX_SoilNox_Final
    USE HcoX_DustDead_Mod,   ONLY : HcoX_DustDead_Final
    USE HcoX_DustGinoux_Mod, ONLY : HcoX_DustGinoux_Final
    USE HcoX_SeaSalt_Mod,    ONLY : HcoX_SeaSalt_Final
    USE HcoX_Megan_Mod,      ONLY : HcoX_Megan_Final
    USE HcoX_GFED3_Mod,      ONLY : HcoX_GFED3_Final
    USE HcoX_GC_RnPbBe_Mod,  ONLY : HcoX_GC_RnPbBe_Final
    USE Hcox_ExtList_Mod,    ONLY : ExtFinal
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

    ! Call individual cleanup routines
    IF ( ExtState%Custom     ) CALL HcoX_Custom_Final()
    IF ( ExtState%SeaFlux    ) CALL HcoX_SeaFlux_Final()
    IF ( ExtState%ParaNOx    ) CALL HcoX_PARANOX_Final()
    IF ( ExtState%LightNOx   ) CALL HcoX_LIGHTNOX_Final()
    IF ( ExtState%DustDead   ) CALL HcoX_DustDead_Final()
    IF ( ExtState%DustGinoux ) CALL HcoX_DustGinoux_Final()
    IF ( ExtState%SeaSalt    ) CALL HcoX_SeaSalt_Final()
    IF ( ExtState%Megan      ) CALL HcoX_Megan_Final()
    IF ( ExtState%GFED3      ) CALL HcoX_GFED3_Final()
    IF ( ExtState%SoilNOx    ) CALL HcoX_SoilNox_Final()  
    IF ( ExtState%GC_RnPbBe  ) CALL HcoX_GC_RnPbBe_Final()
     
    ! Remove ExtState object 
    CALL ExtStateFinal( ExtState ) 

    ! Remove extensions list 
    CALL ExtFinal()

  END SUBROUTINE HCOX_FINAL
!EOC
END MODULE HCOX_DRIVER_MOD
