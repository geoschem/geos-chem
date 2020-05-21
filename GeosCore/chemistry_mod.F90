!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod.F90
!
! !DESCRIPTION: Module CHEMISTRY\_MOD is used to call the proper chemistry
!  subroutine for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE:
!
MODULE Chemistry_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)
  USE Timers_Mod       ! For GEOS-Chem timers (optional)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: INIT_CHEMISTRY
  PUBLIC  :: DO_CHEMISTRY
  PUBLIC  :: RECOMPUTE_OD
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHEM_PASSIVE_SPECIES
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

  INTEGER :: id_DST1, id_NK1   ! Species ID flags

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_chemistry
!
! !DESCRIPTION: Subroutine DO\_CHEMISTRY is the driver routine which calls
!  the appropriate chemistry subroutine for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_CHEMISTRY( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE AEROSOL_MOD,     ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,     ONLY : RDAER
    USE AEROSOL_MOD,     ONLY : SOILDUST
    USE CARBON_MOD,      ONLY : CHEMCARBON
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE DUST_MOD,        ONLY : CHEMDUST
    USE DUST_MOD,        ONLY : RDUST_ONLINE
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FlexChem_Mod,    ONLY : Do_FlexChem
    USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
    USE Input_Opt_Mod,   ONLY : OptInput
    USE ISORROPIAII_MOD, ONLY : DO_ISORROPIAII
    USE MERCURY_MOD,     ONLY : CHEMMERCURY
    USE POPS_MOD,        ONLY : CHEMPOPS
    USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
    USE RPMARES_MOD,     ONLY : DO_RPMARES
    USE SEASALT_MOD,     ONLY : CHEMSEASALT
    USE SULFATE_MOD,     ONLY : CHEMSULFATE
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Chm_Mod,   ONLY : Ind_
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE STRAT_CHEM_MOD,  ONLY : DO_STRAT_CHEM
    USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
    USE TAGGED_O3_MOD,   ONLY : CHEM_TAGGED_O3
    USE TIME_MOD,        ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,        ONLY : GET_TS_CHEM
    USE UCX_MOD,         ONLY : CALC_STRAT_AER
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units
#ifdef APM
    USE APM_INIT_MOD,    ONLY : APMIDS
    USE APM_DRIV_MOD,    ONLY : PSO4GAS
    USE APM_DRIV_MOD,    ONLY : AERONUM
    USE APM_DRIV_MOD,    ONLY : APM_DRIV
#endif
#ifdef TOMAS
    USE TOMAS_MOD,       ONLY : DO_TOMAS  !(win, 7/14/09)
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N_TROP, N
    INTEGER            :: MONTH
    INTEGER            :: YEAR
    INTEGER            :: WAVELENGTH
    LOGICAL            :: IT_IS_A_CH4_SIM
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_A_MERCURY_SIM
    LOGICAL            :: IT_IS_A_RnPbBe_SIM
    LOGICAL            :: IT_IS_A_TAGCO_SIM
    LOGICAL            :: IT_IS_A_TAGO3_SIM
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: IT_IS_NOT_COPARAM_OR_CH4
    LOGICAL            :: IT_IS_A_POPS_SIM
    LOGICAL            :: LCARB
    LOGICAL            :: LCHEM
    LOGICAL            :: LDUST
    LOGICAL            :: LSCHEM
    LOGICAL            :: LSSALT
    LOGICAL            :: LSULF
    LOGICAL            :: LSOA
    LOGICAL            :: LNLPBL
    LOGICAL            :: LUCX
    LOGICAL            :: prtDebug
    REAL(fp)           :: DT_Chem
#ifdef APM
    INTEGER            :: I,J,L
    REAL*8             :: CONCTMPSO4(State_Grid%NX,                         &
                                     State_Grid%NY,                         &
                                     State_Grid%NZ)
#endif

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! DO_CHEMISTRY begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                    = Input_Opt%LCARB
    LCHEM                    = Input_Opt%LCHEM
    LDUST                    = Input_Opt%LDUST
    LSCHEM                   = Input_Opt%LSCHEM
    LSSALT                   = Input_Opt%LSSALT
    LSULF                    = Input_Opt%LSULF
    LSOA                     = Input_Opt%LSOA
    LNLPBL                   = Input_Opt%LNLPBL
    LUCX                     = Input_Opt%LUCX
    IT_IS_A_CH4_SIM          = Input_Opt%ITS_A_CH4_SIM
    IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_A_MERCURY_SIM      = Input_Opt%ITS_A_MERCURY_SIM
    IT_IS_A_RnPbBe_SIM       = Input_Opt%ITS_A_RnPbBe_SIM
    IT_IS_A_TAGCO_SIM        = Input_Opt%ITS_A_TAGCO_SIM
    IT_IS_A_TAGO3_SIM        = Input_Opt%ITS_A_TAGO3_SIM
    IT_IS_A_POPS_SIM         = Input_Opt%ITS_A_POPS_SIM
    IT_IS_AN_AEROSOL_SIM     = Input_Opt%ITS_AN_AEROSOL_SIM
    prtDebug                 = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Save species ID"s on first call
    IF ( FIRST ) THEN
       id_DST1 = Ind_('DST1')
       id_NK1  = Ind_('NK1' )
    ENDIF

    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN
       ! Get initial column masses
       CALL Compute_Column_Mass( Input_Opt,                              &
                                 State_Chm, State_Grid, State_Met,       &
                                 State_Chm%Map_Advect,                   &
                                 State_Diag%Archive_BudgetChemistryFull, &
                                 State_Diag%Archive_BudgetChemistryTrop, &
                                 State_Diag%Archive_BudgetChemistryPBL,  &
                                 State_Diag%BudgetMass1,                 &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Convert species units to [kg] for chemistry (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> kg)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! If LCHEM=T then call the chemistry subroutines
    !=======================================================================
    IF ( LCHEM ) THEN

       !====================================================================
       ! Full-chemistry simulations:
       !
       ! (1) Benchmark; (2) Standard; (3) SimpleSOA; (4) complexSOA,
       ! (5) complexSOA-SVPOA; (6) aciduptake; (7) marinePOA
       !====================================================================
       IF ( IT_IS_A_FULLCHEM_SIM ) THEN

          !----------------------------------------
          ! Dry-run sulfate chem to get cloud pH
          !----------------------------------------
          IF ( LSULF ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "=> Aerosol chem", RC )
             ENDIF

             ! Dry run only
             CALL ChemSulfate( Input_Opt, State_Chm, State_Diag, State_Grid, &
                               State_Met, .FALSE.,   RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "=> Aerosol chem", RC )
             ENDIF

          ENDIF

#ifdef APM
          ! Save SO4 concentration before chemistry
          N          = APMIDS%id_SO4
          CONCTMPSO4 = State_Chm%Species(:,:,:,N)

          CALL AERONUM( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemSulfate"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

          !---------------------------
          ! Call gas-phase chemistry
          !---------------------------
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> FlexChem", RC )
          ENDIF

          CALL Do_FlexChem( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after FLEX_CHEMDR!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Do_FlexChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> FlexChem", RC )
          ENDIF

          !----------------------------------------
          ! Call linearized stratospheric scheme
          !----------------------------------------
          IF ( LSCHEM ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "=> Linearized strat chem", RC )
             ENDIF

             ! Do linearized chemistry for the stratosphere (tropchem)
             ! or the mesosphere (UCX)
             CALL Do_Strat_Chem( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_STRAT_CHEM!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "=> Linearized strat chem", RC )
             ENDIF

          ENDIF

#ifdef APM
          ! Obtain SO4 production after chemistry
          N = APMIDS%id_SO4
          !$OMP PARALLEL DO         &
          !$OMP DEFAULT( SHARED   ) &
          !$OMP PRIVATE( I, J, L  ) &
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             IF ( State_Chm%Species(I,J,L,N) > CONCTMPSO4(I,J,L) ) THEN
                PSO4GAS(I,J,L) = State_Chm%Species(I,J,L,N)                  &
                               - CONCTMPSO4(I,J,L)
             ELSE
                PSO4GAS(I,J,L) = 0.D0
             ENDIF
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Aerosol chem", RC )
          ENDIF

          !--------------------------------
          ! Do seasalt aerosol chemistry
          !--------------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------------------
          ! Recalculate PSC properties
          !-------------------------------
          IF ( LUCX ) THEN

             ! Recalculate PSC
             CALL Calc_Strat_Aer( Input_Opt, State_Chm, State_Grid, &
                                  State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ENDIF

          !--------------------------------
          ! Also do sulfate chemistry
          !--------------------------------
          IF ( LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( Input_Opt, State_Chm, State_Diag, State_Grid, &
                               State_Met, .TRUE.,     RC )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg =  'Incorrect species units after CHEMSULFATE!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !-----------------------------------------
             ! Do aerosol thermodynamic equilibrium
             !-----------------------------------------
             IF ( LSSALT ) THEN

#ifndef APM
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( Input_Opt,  State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_ISORROPIAII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
#endif

             ELSE

#ifdef APM
                ! Exit with error if RPMARES + APM is selected
                ErrMsg = 'Warning: APM does not want to use DO_RPMARES'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
#endif

                ! RPMARES does not take Na+, Cl- into account
                CALL Do_RPMARES( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )
             ENDIF

          ENDIF

          !-----------------------------------
          ! Do carbonaceous aerosol chemistry
          !-----------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( Input_Opt,  State_Chm, State_Diag, &
                              State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemCarbon"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------------------
          ! Do dust aerosol chemistry/removal
          !------------------------------------
          IF ( LDUST .AND. id_DST1 > 0 ) THEN
             CALL ChemDust( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

#ifdef APM
          !--------------------------------------------
          ! Do APM aerosol microphysics
          !--------------------------------------------
          CALL APM_DRIV( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "APM_DRIV"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
#endif

#ifdef TOMAS
          !--------------------------------------------
          ! Do TOMAS aerosol microphysics and dry dep
          !--------------------------------------------
          IF ( id_NK1 > 0 ) THEN
             CALL Do_TOMAS( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

             ! Check units (ewl, 10/5/15)
             IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
                ErrMsg = 'Incorrect species units after DO_TOMAS!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
             ENDIF

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_TOMAS"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF
#endif

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Aerosol chem", RC )
          ENDIF

       !====================================================================
       ! Aerosol-only simulation
       !====================================================================
       ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "=> Aerosol chem", RC )
          ENDIF

          !-------------------------------------------------------
          ! Compute aerosol & dust concentrations [kg/m3]
          ! (NOTE: SOILDUST in "aerosol_mod.F90" is computed here)
          !-------------------------------------------------------
          CALL Aerosol_Conc( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, State_Met, RC )

          ! Check units (ewl, 10/5/15)
          IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
             ErrMsg = 'Incorrect species units after AEROSOL_CONC'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
          ENDIF

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Aerosol_Conc"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-------------------------------------------
          ! Compute AOD's and surface areas at 999 nm
          !-------------------------------------------
          MONTH      = 0
          YEAR       = 0
          WAVELENGTH = 0
          CALL RdAer( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC,         &
                      MONTH,      YEAR,      WAVELENGTH  )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "RdAer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !--------------------------------------------
          ! Aerosol Thermodynamic Equilibrium
          !--------------------------------------------
          IF ( LSULF ) THEN
             IF ( LSSALT ) THEN

#ifndef APM
                ! ISORROPIA takes Na+, Cl- into account
                CALL Do_IsorropiaII( Input_Opt,  State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )
#endif

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_IsorropiaII"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

             ELSE

#ifdef APM
                ! Exit with error if RPMARES + APM is selected
                ErrMsg = 'Warning: APM does not want to use DO_RPMARES'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
#endif

                ! RPMARES does not take Na+, Cl- into account
                ! (skip for crystalline & aqueous offline run)
                CALL Do_RPMARES( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Do_RPMARES"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF
          ENDIF

          !-----------------------------
          ! Seasalt Aerosols
          !-----------------------------
          IF ( LSSALT ) THEN
             CALL ChemSeaSalt( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSeaSalt"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-------------------
          ! Sulfate aerosols
          !-------------------
          IF ( LSULF ) THEN

             ! Do sulfate chemistry
             CALL ChemSulfate( Input_Opt, State_Chm, State_Diag, State_Grid, &
                               State_Met, .TRUE.,    RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemSulfate"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !-----------------------------------------
          ! Carbon and Secondary Organic Aerosols
          !-----------------------------------------
          IF ( LCARB ) THEN
             CALL ChemCarbon( Input_Opt,  State_Chm, State_Diag, &
                              State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in ""!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          !------------------------
          ! Mineral Dust Aerosols
          !------------------------
          IF ( LDUST ) THEN

             ! Do dust aerosol chemistry
             CALL ChemDust( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "ChemDust"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             ! Compute dust OD's & surface areas
             WAVELENGTH = 0
             CALL Rdust_Online( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, SOILDUST,   &
                                WAVELENGTH, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Rdust_Online"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF
          ENDIF

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "=> Aerosol chem", RC )
          ENDIF

       !====================================================================
       ! Rn-Pb-Be
       !====================================================================
       ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN

          ! Do Rn-Pb-Be chemistry
          CALL ChemRnPbBe( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemRnPbBe"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! Tagged O3
       !====================================================================
       ELSE IF ( IT_IS_A_TAGO3_SIM ) THEN

          !-----------------------------------------------
          ! Do Tagged O3 chemistry
          !-----------------------------------------------
          CALL Chem_Tagged_O3( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_O3"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-----------------------------------------------
          ! Call linearized stratospheric scheme (LINOZ)
          !-----------------------------------------------
          IF ( LSCHEM ) THEN

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_Start( "=> Linearized strat chem", RC )
             ENDIF

             ! Do LINOZ for Ozone
             CALL Do_Strat_Chem( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "Do_Strat_Chem"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             IF ( Input_Opt%useTimers ) THEN
                CALL Timer_End( "=> Linearized strat chem", RC )
             ENDIF

          ENDIF

       !====================================================================
       ! Tagged CO
       !====================================================================
       ELSE IF ( IT_IS_A_TAGCO_SIM ) THEN

          ! Do tagged CO chemistry
          CALL Chem_Tagged_CO( Input_Opt,  State_Chm, State_Diag, &
                               State_Grid, State_Met, RC        )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Tagged_CO"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! CH4
       !====================================================================
       ELSE IF ( IT_IS_A_CH4_SIM ) THEN

          CALL ChemCh4( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemCh4"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! Mercury (only used when compiled with BPCH_DIAG=y)
       !====================================================================
       ELSE IF ( IT_IS_A_MERCURY_SIM ) THEN

          ! Do Hg chemistry
          CALL ChemMercury( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemMercury"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       !====================================================================
       ! POPs (only used when compiled with BPCH_DIAG=y)
       !====================================================================
       ELSE IF ( IT_IS_A_POPS_SIM ) THEN

          ! Do POPS chemistry
          CALL ChemPOPs( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ChemPOPs"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDIF

       !====================================================================
       ! PASSIVE SPECIES
       !
       ! This performs a simple loss chemistry on passive species.  Call
       ! this routine for all simulation types since passive species can
       ! be defined for various simulations (as additional species to the
       ! default! ones). ckeller, 09/04/15
       !
       ! NOTE: To speed up execution, only call Chem_Passive_Species
       ! if there is at least one passive species with a finite
       ! atmospheric lifetime.  There is no reason to apply a loss rate
       ! of unity to those passive species whose lifetime is infinity.
       ! This will speed up GEOS-Chem simulations. (bmy, 12/13/17)
       !====================================================================
       IF ( Input_Opt%NPassive_Decay > 0 ) THEN

          ! Apply loss rate to passive species with finite lifetimes
          CALL Chem_Passive_Species( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_Passive_Species"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !### Debug
          IF ( prtDebug ) THEN
             CALL Debug_Msg( '### MAIN: a CHEMISTRY' )
          ENDIF

       ENDIF

    ENDIF

    !=======================================================================
    ! Convert species units back to original unit (ewl, 8/12/15)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Chemistry timestep [s]
    DT_Chem = Get_Ts_Chem()

    !----------------------------------------------------------
    ! Chemistry budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetChemistry ) THEN
       ! Get final column masses and compute diagnostics
       CALL Compute_Column_Mass( Input_Opt,                              &
                                 State_Chm, State_Grid, State_Met,       &
                                 State_Chm%Map_Advect,                   &
                                 State_Diag%Archive_BudgetChemistryFull, &
                                 State_Diag%Archive_BudgetChemistryTrop, &
                                 State_Diag%Archive_BudgetChemistryPBL,  &
                                 State_Diag%BudgetMass2,                 &
                                 RC )
       CALL Compute_Budget_Diagnostics( State_Grid,                          &
                                     State_Chm%Map_Advect,                   &
                                     DT_Chem,                                &
                                     State_Diag%Archive_BudgetChemistryFull, &
                                     State_Diag%Archive_BudgetChemistryTrop, &
                                     State_Diag%Archive_BudgetChemistryPBL,  &
                                     State_Diag%BudgetChemistryFull,         &
                                     State_Diag%BudgetChemistryTrop,         &
                                     State_Diag%BudgetChemistryPBL,          &
                                     State_Diag%BudgetMass1,                 &
                                     State_Diag%BudgetMass2,                 &
                                     RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Chemistry budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE DO_CHEMISTRY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: recompute_od
!
! !DESCRIPTION: Subroutine RECOMPUTE\_OD will update the optical depth values
!  before accumulating or writing the diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RECOMPUTE_OD( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE AEROSOL_MOD,    ONLY : AEROSOL_CONC
    USE AEROSOL_MOD,    ONLY : RDAER
    USE AEROSOL_MOD,    ONLY : SOILDUST
    USE DUST_MOD,       ONLY : RDUST_ONLINE
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : Debug_Msg
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_MONTH
    USE TIME_MOD,       ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  03 Fev 2011 - Adapted from chemdr.f by skim
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: IT_IS_A_FULLCHEM_SIM
    LOGICAL            :: IT_IS_AN_AEROSOL_SIM
    LOGICAL            :: LCARB, LCHEM,  LDUST
    LOGICAL            :: LSSALT, LSULF, LSOA
    LOGICAL            :: prtDebug
    INTEGER            :: MONTH, YEAR,   WAVELENGTH

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! RECOMPUTE_OD begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Recompute_OD  (in module GeosCore/chemistry_mod.F90)'

    ! Get month and year
    MONTH                = GET_MONTH()
    YEAR                 = GET_YEAR()

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                = Input_Opt%LCARB
    LCHEM                = Input_Opt%LCHEM
    LDUST                = Input_Opt%LDUST
    LSSALT               = Input_Opt%LSSALT
    LSULF                = Input_Opt%LSULF
    LSOA                 = Input_Opt%LSOA
    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
    IT_IS_AN_AEROSOL_SIM = Input_Opt%ITS_AN_AEROSOL_SIM
    prtDebug             = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! First make sure chemistry is turned on
    IF ( LCHEM ) THEN

       ! Then make sure that the simulations use aerosol species
       IF ( IT_IS_A_FULLCHEM_SIM .or. IT_IS_AN_AEROSOL_SIM ) THEN

          ! And then make sure that the aersol species are defined
          IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

             ! Skip this section if all of these are turned off
             CALL AEROSOL_CONC( Input_Opt,  State_Chm, State_Diag, &
                                State_Grid, State_Met, RC )

             !==============================================================
             ! Call RDAER -- computes aerosol optical depths
             !==============================================================

             ! Calculate the AOD at the wavelength specified in jv_spec_aod
             WAVELENGTH = 1
             CALL RDAER( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, State_Met, RC,         &
                         MONTH,     YEAR,       WAVELENGTH  )

             ! Trap potential errors
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error encountered in "RdAer"!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

             !### Debug
             IF ( prtDebug ) THEN
                CALL Debug_Msg( '### RECOMPUTE_OD: after RDAER' )
             ENDIF

             !==============================================================
             ! If LDUST is turned on, then we have online dust aerosol in
             ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order
             ! to compute aerosol optical depth for FAST-JX, etc.
             !
             ! If LDUST is turned off, then we don't have online dust
             ! aerosol in GEOS-CHEM...so read monthly-mean dust files
             ! from disk. (rjp, tdf, bmy, 4/1/04)
             !==============================================================
             IF ( LDUST ) THEN
                CALL RDUST_ONLINE( Input_Opt,  State_Chm, State_Diag, &
                                   State_Grid, State_Met, SOILDUST,   &
                                   WAVELENGTH, RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Rdust_Online"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF

             !### Debug
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDUST' )
             ENDIF
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE RECOMPUTE_OD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_passive_species
!
! !DESCRIPTION: Subroutine CHEM\_PASSIVE\_SPECIES performs loss chemistry
!  on passive species with finite atmospheric lifetimes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Passive_Species( Input_Opt, State_Chm, State_Grid, &
                                   State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt   ! Input options object
    TYPE(GrdState), INTENT(IN   ) :: State_Grid  ! Grid state object
    TYPE(MetState), INTENT(IN   ) :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtDebug
    LOGICAL             :: Is_HalfLife
    INTEGER             :: I,       J,      L
    INTEGER             :: N,       GCId,   Id
    REAL(fp)            :: DT,      Decay,  Rate

    ! SAVEd scalars
    LOGICAL,  SAVE      :: First = .TRUE.

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg,  ThisLoc
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=======================================================================
    ! Chem_Passive_Species begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at Chem_Passive_Species (in module GeosCore/chemistry_mod.F90)'

    DT       = GET_TS_CHEM() ! timestep in seconds

    ! For now, always compute decay using e-folding time
    Is_HalfLife = .FALSE.

    !=======================================================================
    ! Apply decay loss rate only to those passive species that have a
    ! finite atmospheric lifetime (this speeds up execution)
    !=======================================================================

    ! Loop over all decaying passive species
    DO N = 1, Input_Opt%NPassive_Decay

       !----------------------------------
       ! Find the GEOS-Chem species Id
       !----------------------------------

       ! Get the Id of the species in the passive decay menu
       Id   = Input_Opt%Passive_DecayID(N)

       ! Convert this to a GEOS-Chem species Id number
       GcId = Ind_( TRIM( Input_Opt%PASSIVE_NAME(Id) ) )

       ! Make sure the model ID is valid
       IF ( GcId < 0 ) THEN
          ErrMsg = 'Could not find the GEOS-Chem species ID # '        // &
                   'for passive species : '                            // &
                   TRIM( Input_Opt%PASSIVE_NAME(Id) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !----------------------------------
       ! Compute the decay rate
       !----------------------------------

       ! Compute the decay rate for each passive species
       IF ( Is_HalfLife ) THEN
          Decay = ln2 / Input_Opt%PASSIVE_TAU(Id)
       ELSE
          Decay = 1.0 / Input_Opt%PASSIVE_TAU(Id)
       ENDIF
       Rate  = EXP( - DT * Decay )

       !### Debug output
       IF ( First ) THEN
          IF ( prtDebug ) THEN
             WRITE( 6,100 ) ADJUSTL( Input_Opt%PASSIVE_NAME(Id) ),           &
                            GcId, Rate
 100         FORMAT( '     -  Pass. species name, Id, loss rate: ',&
                      a15, i5, 1x, es13.6 )
          ENDIF
       ENDIF

       !----------------------------------
       ! Apply loss
       !----------------------------------

       !$OMP PARALLEL DO                  &
       !$OMP DEFAULT( SHARED            ) &
       !$OMP PRIVATE( I, J, L           )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%Species(I,J,L,GcId) = State_Chm%Species(I,J,L,GcId)      &
                                        * Rate
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO

    ! Reset after the first time
    IF ( First) First = .FALSE.

  END SUBROUTINE Chem_Passive_Species
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_chemistry
!
! !DESCRIPTION: Subroutine INIT\_CHEMISTRY initializes chemistry
! variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Chemistry( Input_Opt,  State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE FAST_JX_MOD,    ONLY : Init_FJX
    USE FlexChem_Mod,   ONLY : Init_FlexChem
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT)  :: State_Diag  ! Diagnostics State object
    INTEGER,        INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 May 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! INIT_CHEMISTRY begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_Chemistry  (in module GeosCore/chemistry_mod.F90)'

    ! Skip if we have already done this
    IF ( FIRST ) THEN

       ! Adjust first flag
       FIRST  = .FALSE.

       ! Define species ID's
       id_DST1 = Ind_( 'DST1' )
       id_NK1  = Ind_( 'NK1'  )

       !--------------------------------------------------------------------
       ! Initialize FlexChem (skip if it is a dry-run)
       !--------------------------------------------------------------------
       IF ( .not. Input_Opt%DryRun ) THEN
          CALL Init_FlexChem( Input_Opt, State_Chm, State_Diag, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Init_FlexChem"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Initialize Fast-JX photolysis
       ! NOTE: we need to call this for a dry-run so that we can get
       ! a list of all of the lookup tables etc. that FAST-JX reads
       !--------------------------------------------------------------------
       CALL Init_FJX( Input_Opt, State_Chm, State_Diag, State_Grid, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Init_FJX"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Init_Chemistry
!EOC
END MODULE Chemistry_Mod
