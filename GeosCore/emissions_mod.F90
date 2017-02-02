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
  SUBROUTINE EMISSIONS_INIT( am_I_Root, Input_Opt, State_Met, State_Chm, &
                             RC,        HcoConfig ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_INIT
    USE HCO_TYPES_MOD,      ONLY : ConfigObj
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )          :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )          :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )          :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)          :: Input_Opt  ! Input opts
    TYPE(ConfigObj),  POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
    INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller   - Initial version 
!  16 Jun 2016 - J. Sheng    - Added tracer index retriever
!  20 Jun 2016 - R. Yantosca - Now define species IDs only in the INIT phase
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Define species ID flags for use in routines below
    id_BrO   = Ind_('BrO'  )
    id_CH4   = Ind_('CH4'  )
    id_CH3Br = Ind_('CH3Br')    

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( am_I_Root, Input_Opt, State_Met, State_Chm, &
                       RC,        HcoConfig=HcoConfig ) 
    IF ( RC/=GC_SUCCESS ) RETURN 

  END SUBROUTINE EMISSIONS_INIT
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
  SUBROUTINE EMISSIONS_RUN( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, EmisTime,  Phase,     RC ) 
!
! !USES:
!
    USE BROMOCARB_MOD,      ONLY : SET_BRO
    USE BROMOCARB_MOD,      ONLY : SET_CH3BR
    USE CARBON_MOD,         ONLY : EMISSCARBON
    USE CO2_MOD,            ONLY : EMISSCO2
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GC_Error
    USE GLOBAL_CH4_MOD,     ONLY : EMISSCH4
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_RUN
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
#if defined ( TOMAS )
    USE CARBON_MOD,         ONLY : EMISSCARBONTOMAS !jkodros
    USE SULFATE_MOD,        ONLY : EMISSSULFATETOMAS !jkodros
#endif

    ! Use old mercury code for now (ckeller, 09/23/2014)
    USE MERCURY_MOD,        ONLY : EMISSMERCURY

    ! For UCX, use Seb's routines for now
#if defined( UCX )
    USE UCX_MOD,            ONLY : EMISS_BASIC
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    LOGICAL,          INTENT(IN   )  :: EmisTime   ! Emissions in this time step
    INTEGER,          INTENT(IN   )  :: Phase      ! Run phase
 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!  13 Nov 2014 - C. Keller    - Added EMISSCARBON (for SESQ and POA)
!  21 Nov 2014 - C. Keller    - Added EMISSVOC to prevent VOC build-up
!                               above tropopause
!  22 Sep 2016 - R. Yantosca  - Don't call EMISSCARBON unless we are doing
!                               a fullchem or aerosol simulation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Run HEMCO. Phase 1 will only update the HEMCO clock and the 
    ! HEMCO data list, phase 2 will perform the emission calculations.
    CALL HCOI_GC_RUN( am_I_Root, Input_Opt, State_Met, State_Chm, & 
                      EmisTime,  Phase,     RC                     ) 
    IF ( RC /= GC_SUCCESS ) RETURN 

    ! The following only needs to be done in phase 2
    IF ( Phase /= 1 ) THEN 

       ! Call carbon emissions module to make sure that sesquiterpene
       ! emissions calculated in HEMCO (SESQ) are passed to the internal
       ! species array in carbon, as well as to ensure that POA emissions
       ! are correctly treated.
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
            Input_Opt%ITS_AN_AEROSOL_SIM ) THEN 
          CALL EMISSCARBON( am_I_Root, Input_Opt, State_Met, RC )
          IF ( RC /= GC_SUCCESS ) RETURN 
       ENDIF

    ! Call TOMAS emission routines (JKodros 6/2/15)
#if defined ( TOMAS )
       CALL EMISSCARBONTOMAS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    
       CALL EMISSSULFATETOMAS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
#endif
   
       ! For CO2 simulation, emissions are not added to STT in mixing_mod.F90 
       ! because the HEMCO CO2 species are not GEOS-Chem tracers. The emissions
       ! thus need to be added explicitly, which is done in EMISSCO2.
       IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
          CALL EMISSCO2( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
          IF ( RC /= GC_SUCCESS ) RETURN 
       ENDIF
   
       ! For CH4 simulation or if CH4 is defined, call EMISSCH4. 
       ! This will get the individual CH4 emission terms (gas, coal, wetlands, 
       ! ...) and write them into the individual emissions arrays defined in
       ! global_ch4_mod (CH4_EMIS). Emissions are all done in mixing_mod, the
       ! call to EMISSCH4 is for backwards consistency, in particular for the
       ! ND58 diagnostics.
       IF ( Input_Opt%ITS_A_CH4_SIM .OR.            &
          ( id_CH4 > 0 .and. Input_Opt%LCH4EMIS ) ) THEN
          CALL EMISSCH4( am_I_Root, Input_Opt, State_Met, RC )
          IF ( RC /= GC_SUCCESS ) RETURN 
       ENDIF
   
       ! For UCX, use Seb's routines for stratospheric species for now.
#if defined( UCX )
       IF ( Input_Opt%LBASICEMIS ) THEN
          CALL EMISS_BASIC( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       ENDIF
#endif

       ! For mercury, use old emissions code for now
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
          CALL EMISSMERCURY ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
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
             CALL SET_CH3BR( am_I_Root, Input_Opt, State_Met, &
                             State_Chm, RC )
          ENDIF
   
          ! ----------------------------------------------------
          ! If selected in input.geos, then set the MBL
          ! concentration of BrO equal to 1 pptv during daytime.
          ! ----------------------------------------------------
          IF ( Input_Opt%LEMIS .AND. ( id_BrO > 0 ) ) THEN
             CALL SET_BRO( am_I_Root, Input_Opt, State_Met, & 
                           State_Chm, RC          )
          ENDIF
   
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
  SUBROUTINE EMISSIONS_FINAL( am_I_Root, ERROR )
!
! !USES:
!
    USE HCOI_GC_MAIN_MOD, ONLY : HCOI_GC_FINAL
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    LOGICAL,          INTENT(IN   )  :: ERROR      ! Cleanup after crash? 
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    CALL HCOI_GC_Final( am_I_Root, ERROR )

  END SUBROUTINE EMISSIONS_FINAL
!EOC
END MODULE EMISSIONS_MOD
