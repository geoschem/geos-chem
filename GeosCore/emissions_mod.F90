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
MODULE Emissions_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Emissions_Init
  PUBLIC :: Emissions_Run
  PUBLIC :: Emissions_Final
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version. 
!  20 Jun 2016 - R. Yantosca - Declare species ID flags as module variables
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES: 
!
  ! Species ID flags
  INTEGER :: id_BrO, id_CH4, id_CH3Br

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_init
!
! !DESCRIPTION: Subroutine EMISSIONS\_INIT calls the HEMCO - GEOS-Chem
! interface initialization routines. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Init( am_I_Root, Input_Opt, State_Met,                &
                             State_Chm, RC,        HcoConfig                ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOI_GC_Main_Mod,   ONLY : HCoi_GC_Init
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )          :: am_I_Root  ! root CPU?
    TYPE(MetState),  INTENT(IN   )          :: State_Met  ! Met state
    TYPE(ChmState),  INTENT(IN   )          :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT)          :: Input_Opt  ! Input opts
    TYPE(ConfigObj), POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
    INTEGER,         INTENT(INOUT)          :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller   - Initial version 
!  16 Jun 2016 - J. Sheng    - Added tracer index retriever
!  20 Jun 2016 - R. Yantosca - Now define species IDs only in the INIT phase
!  22 Jan 2018 - R. Yantosca - Return error code to calling routine
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Init (in module GeosCore/emissions_mod.F90)'

    ! Define species ID flags for use in routines below
    id_BrO   = Ind_('BrO'  )
    id_CH4   = Ind_('CH4'  )
    id_CH3Br = Ind_('CH3Br')    

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( am_I_Root, Input_Opt, State_Met,                      &
                       State_Chm, RC,        HcoConfig=HcoConfig            ) 

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Emissions_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_run
!
! !DESCRIPTION: Subroutine EMISSIONS\_RUN calls the HEMCO - GEOS-Chem
! interface run routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Run( am_I_Root, Input_Opt, State_Met,                 &
                            State_Chm, EmisTime,  Phase,     RC             ) 
!
! !USES:
!
    USE BROMOCARB_MOD,      ONLY : SET_BRO
    USE BROMOCARB_MOD,      ONLY : SET_CH3BR
    USE CARBON_MOD,         ONLY : EMISSCARBON
    USE CO2_MOD,            ONLY : EMISSCO2
    USE ErrCode_Mod
    USE GLOBAL_CH4_MOD,     ONLY : EMISSCH4
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_RUN
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
#if defined ( TOMAS )
    USE CARBON_MOD,         ONLY : EMISSCARBONTOMAS !jkodros
    USE SULFATE_MOD,        ONLY : EMISSSULFATETOMAS !jkodros
#endif

    ! Setting other surface VMRs
    Use sfcVMR_Mod,         Only : fixSfcVMR

    ! Use old mercury code for now (ckeller, 09/23/2014)
    USE MERCURY_MOD,        ONLY : EMISSMERCURY

    ! For UCX, use Seb's routines for now
    USE UCX_MOD,            ONLY : EMISS_BASIC
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   )  :: am_I_Root  ! root CPU?
    LOGICAL,        INTENT(IN   )  :: EmisTime   ! Emissions in this time step
    INTEGER,        INTENT(IN   )  :: Phase      ! Run phase
 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,        INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller   - Initial version 
!  13 Nov 2014 - C. Keller   - Added EMISSCARBON (for SESQ and POA)
!  21 Nov 2014 - C. Keller   - Added EMISSVOC to prevent VOC build-up
!                              above tropopause
!  22 Sep 2016 - R. Yantosca - Don't call EMISSCARBON unless we are doing
!                              a fullchem or aerosol simulation
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  22 Jan 2018 - R. Yantosca - Return error code to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Run (in module GeosCore/emissions_mod.F90)'

    ! Run HEMCO. Phase 1 will only update the HEMCO clock and the 
    ! HEMCO data list, phase 2 will perform the emission calculations.
    CALL HCOI_GC_Run( am_I_Root, Input_Opt, State_Met, State_Chm, & 
                      EmisTime,  Phase,     RC                     ) 

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Run"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! The following only needs to be done in phase 2
    IF ( Phase /= 1 ) THEN 

       ! Call carbon emissions module to make sure that sesquiterpene
       ! emissions calculated in HEMCO (SESQ) are passed to the internal
       ! species array in carbon, as well as to ensure that POA emissions
       ! are correctly treated.
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN 
          CALL EmissCarbon( am_I_Root, Input_Opt, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "EmissCarbon"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ! Call TOMAS emission routines (JKodros 6/2/15)
#if defined ( TOMAS )
       CALL EmissCarbonTomas( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCarbonTomas"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL EmissSulfateTomas( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissSulfateTomas"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

#endif
   
       ! For CO2 simulation, emissions are not added to STT in mixing_mod.F90 
       ! because the HEMCO CO2 species are not GEOS-Chem tracers. The emissions
       ! thus need to be added explicitly, which is done in EMISSCO2.
       IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
          CALL EmissCO2( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "EmissCO2"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF
   
       ! For CH4 simulation or if CH4 is defined, call EMISSCH4. 
       ! This will get the individual CH4 emission terms (gas, coal, wetlands, 
       ! ...) and write them into the individual emissions arrays defined in
       ! global_ch4_mod (CH4_EMIS). Emissions are all done in mixing_mod, the
       ! call to EMISSCH4 is for backwards consistency, in particular for the
       ! ND58 diagnostics.
       IF ( Input_Opt%ITS_A_CH4_SIM .OR.            &
          ( id_CH4 > 0 .and. Input_Opt%LCH4EMIS ) ) THEN
          CALL EmissCh4( am_I_Root, Input_Opt, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "EmissCH4"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF
   
       ! For UCX, use Seb's routines for stratospheric species for now.
       IF ( Input_Opt%LUCX .and. Input_Opt%LBASICEMIS ) THEN
          CALL Emiss_Basic( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Emiss_Basic"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! For mercury, use old emissions code for now
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
          CALL EmissMercury( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "EmissMercury"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Prescribe some concentrations if needed
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
  
          !========================================================
          !jpp, 2/12/08: putting a call to SET_CH3Br
          !              which is in bromocarb_mod.f
          !       ***** Fix CH3Br Concentration in PBL *****
          ! Kludge: eventually I want to keep the concentration
          !         entirely fixed! Ask around on how to...
          !========================================================
          IF ( Input_Opt%LEMIS .AND. ( id_CH3Br > 0 ) ) THEN
             CALL Set_CH3Br( am_I_Root, Input_Opt, State_Met, &
                             State_Chm, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Set_CH3BR"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
   
          ! ----------------------------------------------------
          ! If selected in input.geos, then set the MBL
          ! concentration of BrO equal to 1 pptv during daytime.
          ! ----------------------------------------------------
          IF ( Input_Opt%LEMIS .AND. ( id_BrO > 0 ) ) THEN
             CALL Set_BrO( am_I_Root, Input_Opt, State_Met, & 
                           State_Chm, RC          )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Set_BrO"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
   
          ! Set other (non-UCX) fixed VMRs
          If ( Input_Opt%LEMIS ) Then
             CALL FixSfcVMR( am_I_Root, Input_Opt, State_Met, & 
                             State_Chm, RC          )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "FixSfcVmr"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          endif
 
       ENDIF
    ENDIF ! Phase/=1  
    
    ! Return w/ success
    RC = GC_SUCCESS
   
   END SUBROUTINE EMISSIONS_RUN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_final
!
! !DESCRIPTION: Subroutine EMISSIONS\_FINAL calls the HEMCO - GEOS-Chem
! interface finalization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Final( am_I_Root, Error, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: am_I_Root  ! root CPU?
    LOGICAL, INTENT(IN)  :: Error      ! Cleanup arrays after crash? 
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCOI_GC_Final (in module GeosCore/hcoi_gc_final_mod.F90)'

    ! Finalize HEMCO
    CALL HCOI_GC_Final( am_I_Root, Error, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Emissions_Final
!EOC
END MODULE Emissions_Mod
