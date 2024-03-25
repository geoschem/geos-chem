!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aerosol_mod.F90
!
! !DESCRIPTION: Module AEROSOL\_MOD contains variables and routines for
!  computing optical properties for aerosols which are needed for both the
!  FAST-J photolysis and ND21 optical depth diagnostics. (bmy, 7/20/04,
!  2/10/09)
!\\
!\\
! !INTERFACE:
!
MODULE AEROSOL_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_AEROSOL
  PUBLIC :: AEROSOL_CONC
  PUBLIC :: RDAER
  PUBLIC :: RD_AOD   ! Public for use in legacy FAST-JX initialization
  PUBLIC :: CALC_AOD ! needed for Hg simulation
!
! !PUBLIC DATA MEMBERS:
!
  ! Growth factors
  REAL(fp), PUBLIC  :: SIA_GROWTH
  REAL(fp), PUBLIC  :: ORG_GROWTH
  REAL(fp), PUBLIC  :: SSA_GROWTH

  ! Logical flags
  LOGICAL,  PUBLIC  :: IS_OCPI
  LOGICAL,  PUBLIC  :: IS_OCPO
  LOGICAL,  PUBLIC  :: IS_BC
  LOGICAL,  PUBLIC  :: IS_SO4
  LOGICAL,  PUBLIC  :: IS_HMS
  LOGICAL,  PUBLIC  :: IS_NH4
  LOGICAL,  PUBLIC  :: IS_NIT
  LOGICAL,  PUBLIC  :: IS_DST
  LOGICAL,  PUBLIC  :: IS_SAL
  LOGICAL,  PUBLIC  :: IS_POA
  LOGICAL,  PUBLIC  :: IS_OPOA
  LOGICAL,  PUBLIC  :: IS_TSOA
  LOGICAL,  PUBLIC  :: IS_ASOA
  LOGICAL,  PUBLIC  :: IS_SOAGX
  LOGICAL,  PUBLIC  :: IS_SimpleSOA
  LOGICAL,  PUBLIC  :: IS_ComplexSOA
!
! !DEFINED PARAMETERS:
!
  ! For SOAGX, assume the total aerosol mass/glyoxal mass = 1.d0
  ! for now (tmf, 1/7/09)
  REAL(fp), PARAMETER,   PUBLIC :: OCFG = 1.e+0_fp
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Pye, H.O.T., and J.H. Seinfeld, "A global perspective on aerosol from
!        low-volatility organic compounds", Atmos. Chem. & Phys., Vol 10, pp
!        4377-4401, 2010.
!  (2 ) Philip, S., R.V. Martin, J.R. Pierce, J.L. Jimenez, Q. Zhang, M.R.
!       Canagaratna, D.V. Spracklen, C.R. Nowlan, L.N. Lamsal, M.J. Cooper, and
!       N.A. Krotkov, "Spatially and seasonally resolved estimate of the ratio
!       of global organic mass to organic carbon", Atmospheric Environment, 87,
!       34-40, doi:10.1016/j.atmosenv.2013.11.065, 2014
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Add tracer ID flags as module variables (bmy, 6/16/16)
  INTEGER :: id_BCPI,  id_BCPO,  id_DST1,  id_DST2
  INTEGER :: id_DST3,  id_DST4,  id_NH4,   id_NIT
  INTEGER :: id_OCPO,  id_OCPI,  id_SALA,  id_SALC
  INTEGER :: id_SO4,   id_SO4s,  id_NITs,  id_NH4s
  INTEGER :: id_POA1,  id_POA2,  id_OPOA1, id_OPOA2
  INTEGER :: id_TSOA1, id_TSOA2, id_TSOA3, id_TSOA0
  INTEGER :: id_ASOAN, id_ASOA1, id_ASOA2, id_ASOA3
  INTEGER :: id_DUST01, id_SOAS,  id_SALACL, id_HMS   ! (jmm, 06/29/18)
  INTEGER :: id_SOAGX, id_SOAIE
  INTEGER :: id_INDIOL,id_LVOCOA

  ! Index to map between NRHAER and species database hygroscopic species
  ! NOTE: Increasing value of NRHAER in CMN_SIZE_Mod.F90 (e.g. if there is
  ! a new hygroscopic species) requires manual update of this mapping
  ! (ewl, 1/23/17)
  INTEGER  :: Map_NRHAER(5)

  
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerosol_conc
!
! !DESCRIPTION: Subroutine AEROSOL\_CONC computes aerosol concentrations in
!  kg/m3 from the tracer mass in kg in the Species array.  These are needed to
!  compute optical properties for photolysis, for the optical depth diagnostics,
!  and for the SOA concentration diagnostics. (bmy, 7/20/04, 2/10/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AEROSOL_CONC( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
#if !defined( MODEL_CESM )
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
#endif
    USE Input_Opt_Mod,     ONLY : OptInput
    USE Species_Mod,       ONLY : SpcConc
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE UnitConv_Mod
    USE TIME_MOD,          ONLY : GET_MONTH
    USE Timers_Mod,        ONLY : Timer_End, Timer_Start
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd variables
    LOGICAL,  SAVE      :: FIRST = .TRUE.

    ! Non-SAVEd variables
    INTEGER             :: I, J, L, N, NA, ND, K, IBINS
    INTEGER             :: k_SO4
    INTEGER             :: k_ORG
    INTEGER             :: k_SSA
    REAL(fp)            :: Rad_wet, Rad_dry
    REAL(fp)            :: Rho_wet, Rho_dry
    REAL(fp)            :: REFF

    ! Logical flags
    LOGICAL             :: LCARB
    LOGICAL             :: LDUST
    LOGICAL             :: LSSALT
    LOGICAL             :: LSULF

    ! Pointers
    TYPE(SpcConc), POINTER   :: Spc(:)
    REAL*8,        POINTER   :: REAA(:,:,:)
    REAL(fp),      POINTER   :: AIRVOL(:,:,:)
    REAL(fp),      POINTER   :: PMID(:,:,:)
    REAL(fp),      POINTER   :: T(:,:,:)
    REAL(fp),      POINTER   :: SOILDUST(:,:,:,:)
    REAL(fp),      POINTER   :: KG_STRAT_AER(:,:,:,:)

    ! Other variables
    INTEGER             :: previous_units


    ! For spatially and seasonally varying OM/OC
    CHARACTER(LEN=255)  :: FIELDNAME
    INTEGER             :: MONTH
    LOGICAL             :: FND

    ! For errors
    CHARACTER(LEN=255)  :: ThisLoc
    CHARACTER(LEN=1023) :: ErrMsg

    !=================================================================
    ! AEROSOL_CONC begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at AEROSOL_CONC (in module GeosCore/aerosol_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB   = Input_Opt%LCARB
    LDUST   = Input_Opt%LDUST
    LSSALT  = Input_Opt%LSSALT
    LSULF   = Input_Opt%LSULF

#ifdef TOMAS
    ! Number of size bins for TOMAS microphysics
    IBINS   = State_Chm%nTomasBins
#endif
    ! Set pointers
    REAA => State_Chm%Phot%REAA

    ! Stop aerosol chem timer (so that unit conv can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "=> Aerosol chem", RC )
    ENDIF

    ! Convert species to [kg] for this routine
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = KG_SPECIES,                                        &
         mapping        = State_Chm%Map_Advect,                              &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error at start of AEROSOL_CONC!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Start aerosol chem timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "=> Aerosol chem", RC )
    ENDIF

    ! Initialize pointers
    Spc      => State_Chm%Species
    AIRVOL   => State_Met%AIRVOL
    PMID     => State_Met%PMID
    T        => State_Met%T
    SOILDUST => State_Chm%SoilDust
    KG_STRAT_AER => State_Chm%KG_AER

    !=================================================================
    ! OM/OC ratio
    !
    ! Get spatial and seasonally varying OM/OC from Philip et al. (2014)
    ! or use default global mean values recommended by Aerosols WG
    !=================================================================

    ! Get data for OM/OC for current month from HEMCO
    MONTH = GET_MONTH()
    IF      ( MONTH == 12 .or. MONTH == 1  .or. MONTH == 2  ) THEN
       FieldName = 'OMOC_DJF'
    ELSE IF ( MONTH == 3  .or. MONTH == 4  .or. MONTH == 5  ) THEN
       FieldName = 'OMOC_MAM'
    ELSE IF ( MONTH == 6  .or. MONTH == 7  .or. MONTH == 8  ) THEN
       FieldName = 'OMOC_JJA'
    ELSE IF ( MONTH == 9  .or. MONTH == 10 .or. MONTH == 11 ) THEN
       FieldName = 'OMOC_SON'
    ENDIF

    IF ( RC /= GC_SUCCESS ) RETURN
#if !defined( MODEL_CESM )
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, Trim(FieldName), State_Chm%OMOC, RC, FOUND=FND )
#else
    FND = .True.
    RC  = GC_SUCCESS
#endif

    IF ( RC == GC_SUCCESS .AND. FND ) THEN

       ! Set OM/OC using spatially and seasonally varying data from
       ! Philip et al. (2014)
       State_Chm%AerMass%OCFPOA(:,:)  = State_Chm%OMOC(:,:) ! OM/OC for POA
       State_chm%AerMass%OCFOPOA(:,:) = State_Chm%OMOC(:,:) ! OM/OC for OPOA, OCPI, and OCPO

    ELSE

       ! Use default global mean OM/OC recommended by the Aerosols WG
       State_Chm%AerMass%OCFPOA(:,:)  = 1.4e+0_fp ! OM/OC for POA
       State_chm%AerMass%OCFOPOA(:,:) = 2.1e+0_fp ! OM/OC for OPOA, OCPI, and OCPO

    ENDIF

    ! Save OM/OC
    State_Chm%OMOC_POA(:,:) = State_Chm%AerMass%OCFPOA(:,:)
    State_Chm%OMOC_OPOA(:,:) = State_chm%AerMass%OCFOPOA(:,:)

    !=================================================================
    ! Compute growth factors at 35% RH
    !
    ! GF = 1 + [ ( r_wet / r_dry )^3 -1 ] * [ rho_wet / rho_dry ]
    !
    ! and use rho_wet = 1000 kg/m3
    !=================================================================
    IF ( FIRST ) THEN

       ! Species index of REAA from RD_AOD (in fast_jx_mod.F90)
       k_SO4      = 1
       k_ORG      = 3
       k_SSA      = 4

       ! Density of H2O [kg/m3]
       Rho_wet    = 1000e+0_fp

       ! Growth factor for SO4 + NIT + NH4
       Rad_dry    = REAA(1,k_SO4,State_Chm%Phot%DRg) ! DRg = 6. choice of dry size doesn't affect volume growth ratio (hzhu)
       Rad_wet    = REAA(1,k_SO4,State_Chm%Phot%DRg) + 35e+0_fp * &
                  ( REAA(2,k_SO4,State_Chm%Phot%DRg) - REAA(1,k_SO4,State_Chm%Phot%DRg) ) / 50e+0_fp
       Rho_dry    = State_Chm%SpcData(id_SO4)%Info%Density
       SIA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Force SIA growth to 1.1 to treat as partially crystalline
       SIA_GROWTH = 1.1_fp

       ! Growth factor for OCPI + SOA
       Rad_dry    = REAA(1,k_ORG,State_Chm%Phot%DRg)
       Rad_wet    = REAA(1,k_ORG,State_Chm%Phot%DRg) + 35e+0_fp * &
                  ( REAA(2,k_ORG,State_Chm%Phot%DRg) - REAA(1,k_ORG,State_Chm%Phot%DRg) ) / 50e+0_fp
       IF ( IS_POA ) THEN
          Rho_dry    = State_Chm%SpcData(id_POA1)%Info%Density
       ELSE IF ( IS_OCPI ) THEN
          Rho_dry    = State_Chm%SpcData(id_OCPI)%Info%Density
       ENDIF
       ORG_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Growth factor for SALA
       Rad_dry    = REAA(1,k_SSA,State_Chm%Phot%DRg)
       Rad_wet    = REAA(1,k_SSA,State_Chm%Phot%DRg) + 35e+0_fp * &
                  ( REAA(2,k_SSA,State_Chm%Phot%DRg) - REAA(1,k_SSA,State_Chm%Phot%DRg) ) / 50e+0_fp
       Rho_dry    = State_Chm%SpcData(id_SALA)%Info%Density
       SSA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Print values to log file
       IF ( Input_Opt%Verbose ) THEN
          WRITE( 6,'(a)') 'Growth factors at 35% RH:'
          WRITE( 6, 100 ) SIA_GROWTH, ' for SO4, NIT, and NH4'
          WRITE( 6, 100 ) ORG_GROWTH, ' for OCPI and SOA'
          WRITE( 6, 100 ) SSA_GROWTH, ' for SALA'
100       FORMAT(F5.2,A)
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.

    ENDIF

    !$OMP PARALLEL DO               &
    !$OMP DEFAULT( SHARED )         &
    !$OMP PRIVATE( I, J, L, N, K )  &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !==============================================================
       ! S U L F A T E   A E R O S O L S
       !
       ! Dump hydrophilic aerosols into one array that will be passed
       ! to RDAER and then used for heterogeneous chemistry as well
       ! as photolysis rate calculations interatively.
       !
       ! For the full-chemistry run, If LSULF=F, then we read these
       ! aerosol data from Mian's simulation.  If LSULF=T then we use
       ! the online tracers.
       !
       ! Now assume that all sulfate, ammonium, and nitrate are
       ! hydrophilic but sooner or later we can pass only hydrophilic
       ! aerosols from the thermodynamic calculations for this
       ! purpose.  This dumping should be done before calling INITGAS,
       ! which converts the unit of Spc from kg/box to molec/cm3.
       !
       ! Units of SO4_NH4_NIT are [kg/m3].  (rjp, bmy, 3/23/03)
       !==============================================================
       IF ( LSULF ) THEN

          ! If we are using the full stratospheric chemistry mechanism,
          ! stratospheric NH4 is ignored, stratospheric NIT is taken
          ! as available for NAT formation and stratospheric SO4 is
          ! taken as sulfuric acid
          IF ( State_Met%InTroposphere(I,J,L) ) THEN

             ! Tropospheric - keep as normal
             ! now including sulfate and nitrate associated with sea-salt
             ! NOTE: these should be treated as having a sea-salt size
             ! distribution but are currently simply treated in the same
             ! way (size and optics) as all other sulfate aerosol (DAR
             ! 2013)

             IF ( IS_HMS ) THEN

                !%%%%% Fullchem simulations: add contribution from HMS
                State_Chm%AerMass%SO4_NH4_NIT(I,J,L) = ( Spc(id_SO4)%Conc(I,J,L)     &
                                   +   Spc(id_HMS)%Conc(I,J,L)     &
                                   +   Spc(id_NH4)%Conc(I,J,L)     &
                                   +   Spc(id_NIT)%Conc(I,J,L) )   &
                                   / AIRVOL(I,J,L)

                State_Chm%AerMass%HMS(I,J,L) = Spc(id_HMS)%Conc(I,J,L) / AIRVOL(I,J,L)

             ELSE

                !%%%%% Aerosol-only simulations: Skip contribution from HMS
                State_Chm%AerMass%SO4_NH4_NIT(I,J,L) = ( Spc(id_SO4)%Conc(I,J,L)   &
                                   +   Spc(id_NH4)%Conc(I,J,L)   &
                                   +   Spc(id_NIT)%Conc(I,J,L) ) &
                                   / AIRVOL(I,J,L)

                State_Chm%AerMass%HMS(I,J,L) = 0.0_fp
             ENDIF

             State_Chm%AerMass%SO4(I,J,L) = Spc(id_SO4)%Conc(I,J,L) / AIRVOL(I,J,L)
             State_Chm%AerMass%NH4(I,J,L) = Spc(id_NH4)%Conc(I,J,L) / AIRVOL(I,J,L)
             State_Chm%AerMass%NIT(I,J,L) = Spc(id_NIT)%Conc(I,J,L) / AIRVOL(I,J,L)
             State_Chm%AerMass%SLA(I,J,L) = 0.0_fp
             State_Chm%AerMass%SPA(I,J,L) = 0.0_fp


          ELSE

             ! Tropospheric sulfate is zero in stratosphere
             State_Chm%AerMass%SO4_NH4_NIT(I,J,L) = 0.0_fp
             State_Chm%AerMass%SO4(I,J,L) = 0.0_fp
             State_Chm%AerMass%HMS(I,J,L) = 0.0_fp ! (jmm, 06/30/18)
             State_Chm%AerMass%NH4(I,J,L) = 0.0_fp
             State_Chm%AerMass%NIT(I,J,L) = 0.0_fp
             State_Chm%AerMass%SLA(I,J,L) = KG_STRAT_AER(I,J,L,1) / AIRVOL(I,J,L)
             State_Chm%AerMass%SPA(I,J,L) = KG_STRAT_AER(I,J,L,2) / AIRVOL(I,J,L)
          ENDIF

          ! Add error check for safe division (bmy, 4/7/15)
          IF ( State_Chm%AerMass%SO4_NH4_NIT(I,J,L) > 0e+0_fp ) THEN

             ! Save these fractions for partitioning of optics
             ! until later when these may be treated independently
             ! Only use HMS if it is defined (for fullchem sims)
             IF ( IS_HMS ) THEN
                State_Chm%AerMass%FRAC_SNA(I,J,L,1) = ( ( Spc(id_SO4)%Conc(I,J,L) +         &
                                        Spc(id_HMS)%Conc(I,J,L) )         &
                                  /   AIRVOL(I,J,L)              )        &
                                  / State_Chm%AerMass%SO4_NH4_NIT(I,J,L)
             ELSE
                State_Chm%AerMass%FRAC_SNA(I,J,L,1) = ( Spc(id_SO4)%Conc(I,J,L) / AIRVOL(I,J,L) )&
                                    / State_Chm%AerMass%SO4_NH4_NIT(I,J,L)
             ENDIF


             State_Chm%AerMass%FRAC_SNA(I,J,L,2) = ( Spc(id_NIT)%Conc(I,J,L) / AIRVOL(I,J,L) ) &
                               / State_Chm%AerMass%SO4_NH4_NIT(I,J,L)

             State_Chm%AerMass%FRAC_SNA(I,J,L,3) = ( Spc(id_NH4)%Conc(I,J,L) / AIRVOL(I,J,L) ) &
                               / State_Chm%AerMass%SO4_NH4_NIT(I,J,L)

          ELSE

             ! If SO4_NH4_NIT(I,J,L) is zero, then avoid a div-by-zero
             ! error.  Set all of these to zero because the division
             ! cannot be done.
             State_Chm%AerMass%FRAC_SNA(I,J,L,1) = 0e+0_fp
             State_Chm%AerMass%FRAC_SNA(I,J,L,2) = 0e+0_fp
             State_Chm%AerMass%FRAC_SNA(I,J,L,3) = 0e+0_fp

          ENDIF

       ENDIF

       !==============================================================
       ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
       !
       ! Compute hydrophilic and hydrophobic BC and OC in [kg/m3]
       ! Also add online 2ndary organics if necessary
       !==============================================================
       IF ( LCARB ) THEN

          ! Hydrophilic BC [kg/m3]
          State_Chm%AerMass%BCPI(I,J,L) = Spc(id_BCPI)%Conc(I,J,L) / AIRVOL(I,J,L)

          ! Hydrophobic BC [kg/m3]
          State_Chm%AerMass%BCPO(I,J,L) = Spc(id_BCPO)%Conc(I,J,L) / AIRVOL(I,J,L)

          ! Hydrophobic OC [kg/m3]
          ! SOAupdate: Treat either OCPO (x2.1) or POA (x1.4)
          IF ( IS_POA ) THEN
             State_Chm%AerMass%OCPO(I,J,L) = ( Spc(id_POA1)%Conc(I,J,L)     &
                             + Spc(id_POA2)%Conc(I,J,L) ) &
                           * State_Chm%AerMass%OCFPOA(I,J) / AIRVOL(I,J,L)
          ELSE IF ( IS_OCPO ) THEN
             State_Chm%AerMass%OCPO(I,J,L) = Spc(id_OCPO)%Conc(I,J,L) &
                           * State_chm%AerMass%OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

          ! Hydrophilic OC [kg/m3]
          IF ( IS_OCPI ) THEN
             State_Chm%AerMass%OCPI(I,J,L) = Spc(id_OCPI)%Conc(I,J,L) &
                           * State_chm%AerMass%OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

          ! Now avoid division by zero (bmy, 4/20/04)
          State_Chm%AerMass%BCPI(I,J,L)    = MAX( State_Chm%AerMass%BCPI(I,J,L), 1e-35_fp )
          State_Chm%AerMass%OCPI(I,J,L)    = MAX( State_Chm%AerMass%OCPI(I,J,L), 1e-35_fp )
          State_Chm%AerMass%BCPO(I,J,L)    = MAX( State_Chm%AerMass%BCPO(I,J,L), 1e-35_fp )
          State_Chm%AerMass%OCPO(I,J,L)    = MAX( State_Chm%AerMass%OCPO(I,J,L), 1e-35_fp )

       ENDIF ! LCARB

       !===========================================================
       ! M I N E R A L   D U S T   A E R O S O L S
       !
       ! NOTE: We can do better than this! Currently we carry 4
       ! dust tracers...but het. chem and fast-j use 7 dust bins
       ! hardwired from Ginoux.
       !
       ! Now, I apportion the first dust tracer into four smallest
       ! dust bins equally in mass for het. chem and fast-j.
       !
       ! Maybe we need to think about chaning our fast-j and het.
       ! chem to use just four dust bins or more flexible
       ! calculations depending on the number of dust bins.
       ! (rjp, 03/27/04)
       !
       ! Now splitting mass into bins in fractions derived from
       ! Highwood et al. (2003).  Data is from log-normal fit of
       ! PCASP measurements of Saharan dust (Solid line in Fig.4b)
       ! (dar, 04/25/10) [see Ridley et al., 2012, JGR]
       !
       ! Updated for TOMAS (Jeffrey Pierce, 6/17/14)
       !
       ! Now get dust radius from species database (bmy, 3/16/17)
       !===========================================================
#ifdef TOMAS

       !-----------------------------------------------------------
       ! TOMAS simulations only
       !-----------------------------------------------------------
       IF ( LDUST ) THEN

          ! Zero SOILDUST
          SOILDUST(I,J,L,:) = 0.e0_fp

          ! Loop over the # of TOMAS dust bins
          DO K = 1, State_Chm%nTomasBins

             ! Get the overall species index for species K
             N    = id_DUST01 + K - 1

             ! Effective aerosol radius [m]
             REFF = State_Chm%SpcData(N)%Info%Radius

             ! Bin #1
             IF ( REFF < 0.2e-6_fp ) THEN
                SOILDUST(I,J,L,1) = SOILDUST(I,J,L,1) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #2
             ELSE IF ( REFF < 0.325e-6_fp ) THEN
                SOILDUST(I,J,L,2) = SOILDUST(I,J,L,2) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #3
             ELSE IF ( REFF < 0.6e-6_fp ) THEN
                SOILDUST(I,J,L,3) = SOILDUST(I,J,L,3) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #4
             ELSE IF ( REFF < 1.15e-6_fp ) THEN
                SOILDUST(I,J,L,4) = SOILDUST(I,J,L,4) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #5
             ELSE IF ( REFF < 2.0e-6_fp ) THEN
                SOILDUST(I,J,L,5) = SOILDUST(I,J,L,5) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #6
             ELSE IF ( REFF < 3.25e-6_fp ) THEN
                SOILDUST(I,J,L,6) = SOILDUST(I,J,L,6) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ! Bin #7
             ELSE
                SOILDUST(I,J,L,7) = SOILDUST(I,J,L,7) &
                                    + Spc(N)%Conc(I,J,L) / AIRVOL(I,J,L)

             ENDIF
          ENDDO

       ENDIF
#else

       !-----------------------------------------------------------
       ! Preserve original code for non-TOMAS simulations
       !-----------------------------------------------------------
       IF ( LDUST ) THEN

          ! Lump 1st dust tracer for het chem
          ! Now use dust size distribution scheme to improve PM2.5
          ! surface dust conc over western U.S. (L. Zhang, 6/25/15)
          SOILDUST(I,J,L,1) = 0.007e+0_fp  * Spc(id_DST1)%Conc(I,J,L) &
                              / AIRVOL(I,J,L)
          SOILDUST(I,J,L,2) = 0.0332e+0_fp * Spc(id_DST1)%Conc(I,J,L) &
                              / AIRVOL(I,J,L)
          SOILDUST(I,J,L,3) = 0.2487e+0_fp * Spc(id_DST1)%Conc(I,J,L) &
                              / AIRVOL(I,J,L)
          SOILDUST(I,J,L,4) = 0.7111e+0_fp * Spc(id_DST1)%Conc(I,J,L) &
                              / AIRVOL(I,J,L)

          ! Other hetchem bins
          SOILDUST(I,J,L,5) = Spc(id_DST2)%Conc(I,J,L) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,6) = Spc(id_DST3)%Conc(I,J,L) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,7) = Spc(id_DST4)%Conc(I,J,L) / AIRVOL(I,J,L)

       ENDIF

#endif

       !===========================================================
       ! S E A S A L T   A E R O S O L S
       !
       ! Compute accumulation & coarse mode concentration [kg/m3]
       !===========================================================
       IF ( LSSALT ) THEN

          ! Accumulation mode seasalt aerosol [kg/m3]
          State_Chm%AerMass%SALA(I,J,L) = Spc(id_SALA)%Conc(I,J,L) / AIRVOL(I,J,L)

          ! Coarse mode seasalt aerosol [kg/m3]
          State_Chm%AerMass%SALC(I,J,L) = Spc(id_SALC)%Conc(I,J,L) / AIRVOL(I,J,L)

          ! Fine mode Cl-/sulfate interal mixed [kg/m3]
          State_Chm%AerMass%ACL(I,J,L) = ( Spc(id_SALACL)%Conc(I,J,L) + &
                         Spc(id_SALA)%Conc(I,J,L)*0.45e0_fp)/AIRVOL(I,J,L)

          ! Avoid division by zero
          State_Chm%AerMass%SALA(I,J,L) = MAX( State_Chm%AerMass%SALA(I,J,L), 1e-35_fp )
          State_Chm%AerMass%SALC(I,J,L) = MAX( State_Chm%AerMass%SALC(I,J,L), 1e-35_fp )
          State_Chm%AerMass%ACL(I,J,L) = MAX( State_Chm%AerMass%ACL(I,J,L), 1e-35_fp )

       ENDIF

       !===========================================================
       ! S E C O N D A R Y   O R G A N I C   A E R O S O L S
       !
       ! Compute SOA concentration [kg/m3]
       !===========================================================

       !--------------------------------------------------------
       ! Simple SOA scheme
       !--------------------------------------------------------
       IF ( Is_SimpleSOA ) THEN

          ! Simple SOA [kg/m3]
          State_Chm%AerMass%SOAS(I,J,L) = Spc(id_SOAS)%Conc(I,J,L) / AIRVOL(I,J,L)

       ENDIF

       !--------------------------------------------------------
       ! Complex SOA scheme
       !--------------------------------------------------------
       IF ( Is_ComplexSOA ) THEN

          ! TSOA (terpene SOA) [kg/m3]
          IF ( IS_TSOA ) THEN
             State_Chm%AerMass%TSOA(I,J,L) = ( Spc(id_TSOA1)%Conc(I,J,L)    &
                           + Spc(id_TSOA2)%Conc(I,J,L)    &
                           + Spc(id_TSOA3)%Conc(I,J,L)    &
                           + Spc(id_TSOA0)%Conc(I,J,L) )  &
                           / AIRVOL(I,J,L)
          ENDIF

          ! ASOA (benz, tolu, xyle, + NAP/IVOC SOA) [kg/m3]
          IF ( IS_ASOA ) THEN
             State_Chm%AerMass%ASOA(I,J,L) = ( Spc(id_ASOAN)%Conc(I,J,L)   &
                           + Spc(id_ASOA1)%Conc(I,J,L)   &
                           + Spc(id_ASOA2)%Conc(I,J,L)   &
                           + Spc(id_ASOA3)%Conc(I,J,L) ) &
                           / AIRVOL(I,J,L)
          ENDIF

          ! OPOA [kg/m3]
          IF ( IS_OPOA ) THEN
             State_Chm%AerMass%OPOA(I,J,L) = ( Spc(id_OPOA1)%Conc(I,J,L)    &
                           + Spc(id_OPOA2)%Conc(I,J,L) )  &
                           * State_chm%AerMass%OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF
       ENDIF

       !-------------------------------------------------------
       ! Mass loading of isoprene SOA (ISOAAQ) [kg/m3]
       !-------------------------------------------------------

       ! Glyoxal
       IF ( id_SOAGX > 0 ) THEN
          State_Chm%AerMass%ISOAAQ(I,J,L) = Spc(id_SOAGX)%Conc(I,J,L) / AIRVOL(I,J,L)
       ENDIF

       ! IEPOX
       IF ( id_SOAIE > 0 ) THEN
          State_Chm%AerMass%ISOAAQ(I,J,L) = State_Chm%AerMass%ISOAAQ(I,J,L) &
                          + Spc(id_SOAIE)%Conc(I,J,L) / AIRVOL(I,J,L)
       ENDIF

       !-----------------------------------------------------------------------
       ! Exclude INDIOL from AOD and aerosol mass calculations. This results in
       ! lost mass. As noted in Fisher et al. (2016, ACP), this is a source of
       ! uncertainty and would benefit from an update when more information
       ! about this process becomes available. (eam, jaf, mps, 3/5/18)
       !! SOA from alkyl nitrates (some contribution
       !! from non-isoprene sources)
       !IF ( id_INDIOL > 0 ) THEN
       !   State_Chm%AerMass%ISOAAQ(I,J,L) = State_Chm%AerMass%ISOAAQ(I,J,L) + Spc(id_INDIOL)%Conc(I,J,L) / AIRVOL(I,J,L)
       !ENDIF
       !-----------------------------------------------------------------------

       ! SOA from ISOPOOH oxidation product
       IF ( id_LVOCOA > 0 ) THEN
          State_Chm%AerMass%ISOAAQ(I,J,L) = State_Chm%AerMass%ISOAAQ(I,J,L) &
                          + Spc(id_LVOCOA)%Conc(I,J,L) / AIRVOL(I,J,L)
       ENDIF

       !-------------------------------------------------------
       ! Hydrophilic primary OC plus SOA [kg/m3].
       !
       ! We need to multiply by OCF to account for the mass of
       ! other components which are attached to the OC aerosol.
       ! (rjp, bmy, 7/15/04)
       !
       ! SOAupdate: use 2.1 (OCFOPOA) (hotp 7/21/10)
       !
       ! sfarina - add SOA-Simplified to primary OC.
       !         - IDTSOAS is already mass basis, so only apply
       !           OCFOPOA to IDTOCPI
       !
       ! SOAupdate: Update traditional SOA (hotp 7/21/10)
       ! for new mtp + isop + lumparomivoc (hotp 5/20/10)
       !
       ! %%% IMPORTANT %%%
       ! Note that if complex SOA is used then PM2.5 includes all
       ! the SOA formed in both the Marais et al. and Pye et al.
       ! schemes and may include some double-counting of isoprene SOA.
       ! (Aerosol WG)
       !-------------------------------------------------------

       ! Use simple SOA by default over complex SOA in calculations
       IF ( Is_SimpleSOA ) THEN
          State_Chm%AerMass%OCPISOA(I,J,L) = State_Chm%AerMass%OCPI(I,J,L) + &
                                             State_Chm%AerMass%SOAS(I,J,L)

       ELSEIF ( Is_ComplexSOA ) THEN

          State_Chm%AerMass%OCPISOA(I,J,L) = State_Chm%AerMass%TSOA(I,J,L) + &
                                             State_Chm%AerMass%ASOA(I,J,L)

          IF ( IS_OCPI ) THEN  ! hotp 7/28/10
             State_Chm%AerMass%OCPISOA(I,J,L) = State_Chm%AerMass%OCPISOA(I,J,L) + &
                                                State_Chm%AerMass%OCPI(I,J,L)
          ENDIF

          IF ( IS_OPOA ) THEN ! hotp 7/28/10
             State_Chm%AerMass%OCPISOA(I,J,L) = State_Chm%AerMass%OCPISOA(I,J,L) + &
                                                State_Chm%AerMass%OPOA(I,J,L)
          ENDIF

          ! Add mechanistic isoprene OA (eam, 08/2015)
          ! Skip adding this for Simple SOA (jaf, clh, bmy, 5/17/18)
          ! benchmark OCPISOA follows simpleSOA and
          ! should exculde ISOAAQ to avoid double-counting
          ! (yuanjianz, 8 Jun 2024)
          State_Chm%AerMass%OCPISOA(I,J,L) = State_Chm%AerMass%OCPISOA(I,J,L) + State_Chm%AerMass%ISOAAQ(I,J,L)

       ENDIF

       ! Now avoid division by zero (bmy, 4/20/04)
       State_Chm%AerMass%OCPISOA(I,J,L) = MAX( State_Chm%AerMass%OCPISOA(I,J,L), 1e-35_fp )

       !===========================================================
       ! SOAGX [kg/m3]
       !===========================================================
       IF ( IS_SOAGX ) THEN
          State_Chm%AerMass%SOAGX(I,J,L) = Spc(id_SOAGX)%Conc(I,J,L) * OCFG / AIRVOL(I,J,L)
       ENDIF

       !==============================================================
       ! P A R T I C U L A T E   M A T T E R
       !
       ! See this GEOS-Chem wiki page for the most up-to-date
       ! definitions of PM2.5 and PM10 used in GEOS-Chem:
       !
       ! http://wiki.geos.chem.org/Particulate_Matter_in_GEOS-Chem
       !==============================================================

       ! Particulate matter < 2.5um [kg/m3]
       State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%NH4(I,J,L)        * SIA_GROWTH + &
                     State_Chm%AerMass%NIT(I,J,L)        * SIA_GROWTH + &
                     State_Chm%AerMass%SO4(I,J,L)        * SIA_GROWTH + &
                     State_Chm%AerMass%HMS(I,J,L)        * SIA_GROWTH + &   ! (jmm, 06/30/18)
                     State_Chm%AerMass%BCPI(I,J,L)                    + &
                     State_Chm%AerMass%BCPO(I,J,L)                    + &
                     State_Chm%AerMass%OCPO(I,J,L)                    + &
                     State_Chm%AerMass%SALA(I,J,L)       * SSA_GROWTH + &
                     SOILDUST(I,J,L,1)              + &
                     SOILDUST(I,J,L,2)              + &
                     SOILDUST(I,J,L,3)              + &
                     SOILDUST(I,J,L,4)              + &
                     SOILDUST(I,J,L,5) * 0.3_fp           ! + 30%  of DST2
       ! OCPI is not present in SVPOA simulation
       ! OCPO represents all POA intead (factor*POA)
       IF ( Is_OCPI ) THEN
          State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%PM25(I,J,L) + &
                                          State_Chm%AerMass%OCPI(I,J,L) * ORG_GROWTH
       ENDIF

       ! Include either simple SOA (default) or Complex SOA in
       ! PM2.5 calculation.  In simulations where both Simple SOA and
       ! Complex SOA species are carried (i.e. "benchmark"), then
       ! only the Simple SOA will be added to PM2.5 and PM10, in order
       ! to avoid double-counting. (bmy, 03 Nov 2021)
       IF ( Is_SimpleSOA ) THEN
          State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%PM25(I,J,L) + ( State_Chm%AerMass%SOAS(I,J,L) * ORG_GROWTH )

       ELSE IF ( Is_ComplexSOA ) THEN
          State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%PM25(I,J,L)                 + &
                        State_Chm%AerMass%TSOA(I,J,L)   * ORG_GROWTH  + &
                        State_Chm%AerMass%ASOA(I,J,L)   * ORG_GROWTH  + &
                        State_Chm%AerMass%ISOAAQ(I,J,L) * ORG_GROWTH        ! Includes SOAGX

          ! Need to add OPOA to PM2.5 for complexSOA_SVPOA simulations
          ! -- Maggie Marvin (15 Jul 2020)
          IF ( Is_OPOA ) THEN
             State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%PM25(I,J,L) + ( State_Chm%AerMass%OPOA(I,J,L) * ORG_GROWTH )
          ENDIF
       ENDIF

       ! Particulate matter < 10um [kg/m3]
       State_Chm%AerMass%PM10(I,J,L) = State_Chm%AerMass%PM25(I,J,L) +                    &   ! PM2.5
                     SOILDUST(I,J,L,5) * 0.7_fp     + &   ! + 70%  of DST2
                     SOILDUST(I,J,L,6)              + &   ! + 100% of DST3
                     SOILDUST(I,J,L,7) * 0.9_fp     + &   ! + 90%  of DST4
                     State_Chm%AerMass%SALC(I,J,L)       * SSA_GROWTH

       ! Apply STP correction factor based on ideal gas law
       State_Chm%AerMass%PM25(I,J,L) = State_Chm%AerMass%PM25(I,J,L) * ( 1013.25_fp / PMID(I,J,L) ) * &
                     ( T(I,J,L)   / 298.0_fp    )

       State_Chm%AerMass%PM10(I,J,L) = State_Chm%AerMass%PM10(I,J,L) * ( 1013.25_fp / PMID(I,J,L) ) * &
                     ( T(I,J,L)   / 298.0_fp    )


      !===========================================================
      ! PDER [um] ! (hzhu, 04/05/2024)
      ! Parameterized dry effective radius for SNA and OM
      !===========================================================
      IF ( State_Chm%AerMass%SO4_NH4_NIT(I,J,L) > 0e+0_fp ) THEN
         ! dry SNA and OM mass, in unit of ug/m3
         State_Chm%AerMass%SNAOM(I,J,L) = ( State_Chm%AerMass%SO4_NH4_NIT(I,J,L) + &
                                            State_Chm%AerMass%OCPO(I,J,L) + &
                                            State_Chm%AerMass%OCPISOA(I,J,L) ) * 1.0e+9_fp

         ! ratio between OM and SNA, unitless
         State_Chm%AerMass%R_OMSNA(I,J,L) = ( State_Chm%AerMass%OCPO(I,J,L) + &
                                              State_Chm%AerMass%OCPISOA(I,J,L) ) / &
                                              State_Chm%AerMass%SO4_NH4_NIT(I,J,L)

         ! Parameterized dry effective radius, in unit of um
         State_Chm%AerMass%PDER(I,J,L) = (exp( 4.36_fp + 0.20_fp*log(State_Chm%AerMass%SNAOM(I,J,L)) + 0.065_fp*log(State_Chm%AerMass%R_OMSNA(I,J,L)) ) *0.001_fp )/0.9_fp ;  
         
         IF (State_Chm%AerMass%PDER(I,J,L) == 0.0_fp) THEN
            State_Chm%AerMass%PDER(I,J,L) = 0.005_fp ! give it a small value to avoid divided by 0
         ENDIF

      ELSE
         State_Chm%AerMass%SNAOM(I,J,L) = 0.0_fp;
         State_Chm%AerMass%R_OMSNA(I,J,L) = 0.0_fp;
         State_Chm%AerMass%PDER(I,J,L) = 0.005_fp;
         
      ENDIF
      !===========================================================


#ifdef MODEL_GEOS
       ! PM2.5 sulfates
       IF ( State_Diag%Archive_PM25su ) THEN
          State_Diag%PM25su(I,J,L) = ( State_Chm%AerMass%SO4(I,J,L) * SIA_GROWTH  ) &
                                   * ( 1013.25_fp / PMID(I,J,L) ) &
                                   * ( T(I,J,L)   / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 nitrates
       IF ( State_Diag%Archive_PM25ni ) THEN
          State_Diag%PM25ni(I,J,L) = ( State_Chm%AerMass%NH4(I,J,L) * SIA_GROWTH    &
                                   +   State_Chm%AerMass%NIT(I,J,L) * SIA_GROWTH  ) &
                                   * ( 1013.25_fp / PMID(I,J,L) ) &
                                   * ( T(I,J,L)   / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 BC
       IF ( State_Diag%Archive_PM25bc  ) THEN
          State_Diag%PM25bc(I,J,L) = ( State_Chm%AerMass%BCPI(I,J,L) + State_Chm%AerMass%BCPO(I,J,L) ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 OC
       IF ( State_Diag%Archive_PM25oc  ) THEN
          State_Diag%PM25oc(I,J,L) = ( State_Chm%AerMass%OCPO(I,J,L)                 &
                                   +   State_Chm%AerMass%OCPI(I,J,L) * ORG_GROWTH  ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 dust
       IF ( State_Diag%Archive_PM25du  ) THEN
          State_Diag%PM25du(I,J,L) = ( SOILDUST(I,J,L,1)           &
                                   +   SOILDUST(I,J,L,2)           &
                                   +   SOILDUST(I,J,L,3)           &
                                   +   SOILDUST(I,J,L,4)           &
                                   +   SOILDUST(I,J,L,5) * 0.38  ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 sea salt
       IF ( State_Diag%Archive_PM25ss  ) THEN
          State_Diag%PM25ss(I,J,L) = ( State_Chm%AerMass%SALA(I,J,L) * SSA_GROWTH  ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 SOA
       IF ( State_Diag%Archive_PM25soa ) THEN
          State_Diag%PM25soa(I,J,L) = ( State_Chm%AerMass%TSOA(I,J,L)   * ORG_GROWTH    &
                                    +   State_Chm%AerMass%ASOA(I,J,L)   * ORG_GROWTH    &
                                    +   State_Chm%AerMass%SOAS(I,J,L)   * ORG_GROWTH    &
                                    +   State_Chm%AerMass%ISOAAQ(I,J,L) * ORG_GROWTH  ) &
                                    * ( 1013.25_fp    / PMID(I,J,L) ) &
                                    * ( T(I,J,L)      / 298.0_fp    ) &
                                    * 1.0e+9_fp
       ENDIF

       ! PM2.5 nitrate 
       IF ( State_Diag%Archive_PM25nit ) THEN
          State_Diag%PM25nit(I,J,L) = ( State_Chm%AerMass%NIT(I,J,L) * SIA_GROWTH  ) &
                                    * ( 1013.25_fp / PMID(I,J,L) ) &
                                    * ( T(I,J,L)   / 298.0_fp    ) &
                                    * 1.0e+9_fp
       ENDIF

       ! PM2.5 ammonium 
       IF ( State_Diag%Archive_PM25nh4 ) THEN
          State_Diag%PM25nh4(I,J,L) = ( State_Chm%AerMass%NH4(I,J,L) * SIA_GROWTH  ) &
                                    * ( 1013.25_fp / PMID(I,J,L) ) &
                                    * ( T(I,J,L)   / 298.0_fp    ) &
                                    * 1.0e+9_fp
       ENDIF
#endif

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Stop aerosol chem timer (so that unit conv can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "=> Aerosol chem", RC )
    ENDIF

    ! Convert species back to original unit
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         mapping    = State_Chm%Map_Advect,                                  &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'End of AEROSOL_CONC in aerosol_mod.F90')
       RETURN
    ENDIF

    ! Start aerosol chem timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "=> Aerosol chem", RC )
    ENDIF

    ! Free pointers
    Spc      => NULL()
    REAA     => NULL()
    AIRVOL   => NULL()
    PMID     => NULL()
    T        => NULL()
    SOILDUST => NULL()

  END SUBROUTINE AEROSOL_CONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rdaer
!
! !DESCRIPTION: Subroutine RDAER reads global aerosol concentrations as
!  determined by Mian Chin.  Calculates optical depth at each level for
!  "set\_prof". Also calculates surface area for heterogeneous chemistry. It
!  uses aerosol parameters in FAST-J input file "jv\_spec.dat" for these
!  calculations. (rvm, rjp, tdf, bmy, 11/04/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RDAER( Input_Opt, State_Chm, State_Diag, State_Grid, State_Met, &
                    RC,        MONTH,     YEAR,       ODSWITCH )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NAER, NRH, NDUST, NRHAER, NSTRATAER
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : ERROR_STOP, Safe_Div
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : CONSVAP
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
    USE TIME_MOD,       ONLY : SYSTEM_TIMESTAMP
    USE UCX_MOD,        ONLY : GET_STRAT_OPT
    USE Species_Mod,    ONLY : Species

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        OPTIONAL      :: MONTH       ! # of current month
    INTEGER,        OPTIONAL      :: YEAR        ! 4-digit year
    INTEGER,        OPTIONAL      :: ODSWITCH    ! Logical indicator
                                                 !  = 0: AOD computed
                                                 !       at 999 nm
                                                 !  = 1: AOD computed
                                                 !       at wavelength set
                                                 !       in Radiation Menu
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
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: FIRST = .TRUE.
    LOGICAL             :: LINTERP
    CHARACTER(LEN=16)   :: STAMP
    INTEGER             :: I, J, L, N, R, IRH, W, IRHN, NA, SpcID, g
    INTEGER             :: AA, IWV, IIWV, NWVS, IR, NRT, S
    REAL*4              :: TEMP( State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp)            :: TEMP2(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp)            :: MSDENS(NAER), DRYAREA, VDRY, VH2O
    REAL(f8)            :: XTAU
    REAL*8              :: BCSCAT_AE  !(xnw, 8/24/15)

    ! Variables for speed diagnostics
    INTEGER             :: ITIMEVALS(8)
    REAL*8              :: OLDSECS, NEWSECS
    REAL*8              :: OLDSECST, NEWSECST

    ! Effective radius at RH bins read in from "FJX_spec.dat"
    REAL(fp)            :: RW(NRH)

    ! Effective radius at RH after interpolation
    REAL(fp)            :: REFF

    ! Q at different RH bins read in from "FJX_spec.dat"
    REAL(fp)            :: AW(NRH)
    REAL(fp)            :: QW(NRH)
    REAL(fp)            :: SSW(NRH)
    REAL(fp)            :: ASYW(NRH)
    REAL(fp)            :: AW0
    REAL(fp)            :: QW0
    REAL(fp)            :: SSW0
    REAL(fp)            :: ASYW0
    REAL(fp)            :: DENOM,NUMER

    ! Used to interpolate between sizes
    REAL(fp)            :: FRAC

    ! Change in Q (extinction efficiency)
    REAL(fp)            :: SCALEQ

    ! Change in Radius with RH
    REAL(fp)            :: SCALER

    ! Change in SSA with RH
    REAL(fp)            :: SCALESSA

    ! Change in asym parameter with RH
    REAL(fp)            :: SCALEASY

    ! Change in Optical properties vs RH
    REAL(fp)            :: SCALEA
    REAL(fp)            :: SCALEOD

    ! Change in Vol vs RH
    REAL(fp)            :: SCALEVOL

    ! Convert AbsHum to RelHum
    REAL(fp)            :: TK,CONSEXP,VPRESH2O,RELHUM

    ! Relative Humidities
    REAL(fp),  SAVE     :: RH(NRH)   = (/0e+0_fp,0.5e+0_fp, &
                                         0.7e+0_fp,0.8e+0_fp,0.9e+0_fp/)

    ! Temporary variables
    REAL(fp)            :: RAER, SADSTRAT, RHOSTRAT, XSASTRAT
    INTEGER             :: ISTRAT

    ! Aqueous aerosol volume (cm3/cm3):
    REAL(fp)            :: TAERVOL

    ! Local variables for quantities from Input_Opt
    LOGICAL             :: LCARB
    LOGICAL             :: LSSALT
    LOGICAL             :: LSULF
    LOGICAL             :: LSTRATOD
    LOGICAL             :: LRAD
    LOGICAL             :: LBCAE  ! (xnw, 8/24/15)
    REAL(fp)            :: GF_RH
    REAL(fp)            :: BCAE_1, BCAE_2

    ! Pointers to State_Chm%Phot
    INTEGER,  POINTER   :: IWVREQUIRED(:)
    INTEGER,  POINTER   :: IWVSELECT  (:,:)
    INTEGER,  POINTER   :: IRHARR     (:,:,:)
    REAL*8,   POINTER   :: ACOEF_WV (:)
    REAL*8,   POINTER   :: BCOEF_WV (:)
    REAL*8,   POINTER   :: REAA     (:,:,:)
    REAL*8,   POINTER   :: QQAA     (:,:,:,:)
    REAL*8,   POINTER   :: ALPHAA   (:,:,:,:)
    REAL*8,   POINTER   :: SSAA     (:,:,:,:)
    REAL*8,   POINTER   :: ASYMAA   (:,:,:,:)
    REAL*8,   POINTER   :: ISOPOD   (:,:,:,:)
    REAL*8,   POINTER   :: ODAER    (:,:,:,:,:)
#ifdef RRTMG
    REAL*8,   POINTER   :: RTODAER  (:,:,:,:,:)
    REAL*8,   POINTER   :: RTSSAER  (:,:,:,:,:)
    REAL*8,   POINTER   :: RTASYMAER(:,:,:,:,:)
#endif

    ! Other pointers
    REAL(fp), POINTER   :: BXHEIGHT (:,:,:)
    REAL(fp), POINTER   :: ERADIUS  (:,:,:,:)
    REAL(fp), POINTER   :: TAREA    (:,:,:,:)
    REAL(fp), POINTER   :: WERADIUS (:,:,:,:)
    REAL(fp), POINTER   :: WTAREA   (:,:,:,:)
    REAL(fp), POINTER   :: ACLRADIUS(:,:,:)
    REAL(fp), POINTER   :: ACLAREA  (:,:,:)

    ! For diagnostics
    LOGICAL                :: IsWL1
    LOGICAL                :: IsWL2
    LOGICAL                :: IsWL3
    LOGICAL                :: IsSLA
    LOGICAL                :: IsPSC
    CHARACTER(LEN=255)     :: ErrMsg
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! RDAER begins here!
    !=================================================================

    ! speed diagnostic
    !CALL DATE_AND_TIME( VALUES=ITIMEVALS )
    !OLDSECS=real(ITIMEVALS(5))*3600.0+real(ITIMEVALS(6))*60.0+ &
    !        real(ITIMEVALS(7))+real(ITIMEVALS(8))/1000.0

    ! Assume success
    RC                   = GC_SUCCESS

    ! Copy fields from INPUT_OPT to local variables for use below
    LCARB                = Input_Opt%LCARB
    LSSALT               = Input_Opt%LSSALT
    LSULF                = Input_Opt%LSULF
    LSTRATOD             = Input_Opt%LSTRATOD
    LRAD                 = Input_Opt%LRAD
    LBCAE                = Input_Opt%LBCAE !(xnw, 8/24/15)
    BCAE_1               = Input_Opt%BCAE_1
    BCAE_2               = Input_Opt%BCAE_2

    ! Initialize pointers
    IWVREQUIRED => State_Chm%Phot%IWVREQUIRED ! WL indexes for interpolation
    IWVSELECT   => State_Chm%Phot%IWVSELECT   ! Indexes of requested WLs
    IRHARR      => State_Chm%Phot%IRHARR      ! Relative humidity indexes
    ACOEF_WV    => State_Chm%Phot%ACOEF_WV    ! Coeffs for WL interpolation
    BCOEF_WV    => State_Chm%Phot%BCOEF_WV    ! Coeffs for WL interpolation
    REAA        => State_Chm%Phot%REAA
    QQAA        => State_Chm%Phot%QQAA
    ALPHAA      => State_Chm%Phot%ALPHAA
    SSAA        => State_Chm%Phot%SSAA
    ASYMAA      => State_Chm%Phot%ASYMAA
    ISOPOD      => State_Chm%Phot%ISOPOD      ! Isoprene optical depth
    ODAER       => State_Chm%Phot%ODAER       ! Aerosol optical depth
#ifdef RRTMG
    RTODAER     => State_Chm%Phot%RTODAER     ! Optical dust
    RTSSAER     => State_Chm%Phot%RTSSAER
    RTASYMAER   => State_Chm%Phot%RTASYMAER
#endif
    BXHEIGHT    => State_Met%BXHEIGHT    ! Grid box height [m]
    ERADIUS     => State_Chm%AeroRadi    ! Aerosol Radius [cm]
    TAREA       => State_Chm%AeroArea    ! Aerosol Area [cm2/cm3]
    WERADIUS    => State_Chm%WetAeroRadi ! Wet Aerosol Radius [cm]
    WTAREA      => State_Chm%WetAeroArea ! Wet Aerosol Area [cm2/cm3]
    ACLRADIUS   => State_Chm%AClRadi     ! Fine Cl- Radius [cm]
    ACLAREA     => State_Chm%AClArea     ! Fine Cl- Area [cm2/cm3]

    !=================================================================
    ! S U L F A T E   A E R O S O L S
    !
    ! If LSULF = TRUE, then take the lumped SO4, NH4, NIT
    ! concentrations [kg/m3] computed by AEROSOL_CONC, and save
    ! into WAERSL(:,:,:,1) for use w/ FAST-J and hetchem.  This is
    ! updated every timestep.  (For fullchem and offline runs)
    !
    ! If LSULF = FALSE, then read monthly mean offline sulfate aerosol
    ! concentrations [kg/m3] from disk at the start of each month.
    ! (For fullchem simulations only)
    !=================================================================
    IF ( LSULF ) THEN

       !-----------------------------------
       ! Use online aerosol concentrations
       !-----------------------------------
       IF ( FIRST ) THEN
          IF ( Input_Opt%amIRoot .and. Input_Opt%Verbose ) WRITE( 6, 100 )
100       FORMAT( '     - RDAER: Using online SO4 NH4 NIT!' )
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%AerMass%WAERSL(I,J,L,1) = State_Chm%AerMass%SO4_NH4_NIT(I,J,L)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    !=================================================================
    ! C A R B O N  &  2 n d A R Y   O R G A N I C   A E R O S O L S
    !
    ! If LCARB = TRUE, then take Hydrophilic OC, Hydrophobic OC,
    ! Hydropilic BC, and Hydrophobic BC, and 2ndary organic aerosol
    ! concentrations [kg/m3] that have been computed by AEROSOL_CONC.
    ! Save these into DAERSL and WAERSL for use w/ FAST-J and hetchem.
    ! These fields are updated every chemistry timestep.
    ! (For both fullchem and offline simulations)
    !
    ! If LCARB = FALSE, then read monthly mean carbon aerosol
    ! concentrations [kg/m3] from disk at the start of each month.
    ! (For full chemistry simulations only)
    !=================================================================
    IF ( LCARB ) THEN

       !-----------------------------------
       ! Use online aerosol concentrations
       !-----------------------------------
       IF ( FIRST ) THEN
          IF ( Input_Opt%amIRoot .and. Input_Opt%Verbose ) WRITE( 6, 110 )
110       FORMAT( '     - RDAER: Using online BCPI OCPI BCPO OCPO!' )
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Hydrophilic BC (a.k.a EC) [kg/m3]
          State_Chm%AerMass%WAERSL(I,J,L,2) = State_Chm%AerMass%BCPI(I,J,L)

          ! Hydrophilic OC [kg/m3]
          State_Chm%AerMass%WAERSL(I,J,L,3) = State_Chm%AerMass%OCPISOA(I,J,L)

          ! Hydrophobic BC (a.k.a EC) [kg/m3]
          State_Chm%AerMass%DAERSL(I,J,L,1) = State_Chm%AerMass%BCPO(I,J,L)

          ! Hydrophobic OC [kg/m3]
          State_Chm%AerMass%DAERSL(I,J,L,2) = State_Chm%AerMass%OCPO(I,J,L)

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    !=================================================================
    ! S E A S A L T   A E R O S O L S
    !
    ! If LSSALT = TRUE, then take accumulation and coarse mode
    ! seasalt aerosol concentrations [kg/m3] that are passed from
    ! KPP.  Save these into WAERSL for use w/ FAST-J and
    ! hetchem.  These fields are updated every chemistry timestep.
    ! (For both fullchem and offline simulations)
    !
    ! If LSSALT = FALSE, then read monthly-mean coarse sea-salt
    ! aerosol concentrations [kg/m3] from the binary punch file.
    ! Also merge the coarse sea salt aerosols into a combined bin
    ! rather than carrying them separately.
    ! (For fullchem simulations only)
    !=================================================================
    IF ( LSSALT ) THEN

       !-----------------------------------
       ! Use online aerosol concentrations
       !-----------------------------------
       IF ( FIRST ) THEN
          IF ( Input_Opt%amIRoot .and. Input_Opt%Verbose ) THEN
             WRITE( 6, 120 )
120          FORMAT( '     - RDAER: Using online SALA SALC' )
          ENDIF
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Accumulation mode seasalt aerosol [kg/m3]
          State_Chm%AerMass%WAERSL(I,J,L,4) = State_Chm%AerMass%SALA(I,J,L)

          ! Coarse mode seasalt aerosol [kg/m3]
          State_Chm%AerMass%WAERSL(I,J,L,5) = State_Chm%AerMass%SALC(I,J,L)

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! Transfer stratospheric aerosol data
    ! SDE 04/17/13
    State_Chm%AerMass%WAERSL(:,:,:,NRHAER+1) = State_Chm%AerMass%SLA
    State_Chm%AerMass%WAERSL(:,:,:,NRHAER+2) = State_Chm%AerMass%SPA

    !=================================================================
    ! Calculate optical depth and surface area at each timestep
    ! to account for the change in relative humidity
    !
    ! For the optical depth calculation, this involves carrying the
    ! optical depth at each RH as separate aerosols since OPMIE
    ! treats the phase functions and single scattering albedos
    ! separately. (An alternative would be to rewrite OPMIE)
    !
    ! Scaling is sufficient for the surface area calculation
    !=================================================================
    ! Representative aerosol densities (kg/m3):
    MSDENS(1) = State_Chm%SpcData(id_SO4)%Info%Density
    MSDENS(2) = State_Chm%SpcData(id_BCPI)%Info%Density
    IF ( IS_POA ) THEN
       MSDENS(3) = State_Chm%SpcData(id_POA1)%Info%Density
    ELSE IF ( IS_OCPI ) THEN
       MSDENS(3) = State_Chm%SpcData(id_OCPI)%Info%Density
    ENDIF
    MSDENS(4) = State_Chm%SpcData(id_SALA)%Info%Density
    MSDENS(5) = State_Chm%SpcData(id_SALC)%Info%Density

    ! These default values unused (actively retrieved from ucx_mod)
    MSDENS(NRHAER+1) = 1700.0d0 ! SSA/STS
    MSDENS(NRHAER+2) = 1000.0d0 ! NAT/ice PSC

    !set default values for RH index array
    IRHARR(:,:,:)=1

    ! Empty ODAER before refilling. This is required to make sure that
    ! all gridboxes outside the chemistry grid are zero and do not carry
    ! over values from previous time steps. (ckeller, 10/15/15)
    ODAER(:,:,:,:,:) = 0.0d0

    ! DAR 09/2013
    ! There are two ways RDAER can be called:
    ! (1) When Fast-J requires aerosol optics at 1000nm (ODSWITCH=0)
    ! (2) Before diags are accumulated, from RECOMPUTE_AOD (ODSWITCH=1)
    ! for (1) we just need the optics stored at a single wavelength
    !     not for all user specified wavelengths, hence NWVS=1, IWV=IWV1000
    ! for (2) we need to determine if RRTMG is switched on
    !     if LRAD=true, calculation is over total minus standard wavelengths
    !     in optics dat files (NWVAA-NWVAA0)
    !     if LRAD=false, calculation is for the wavelengths required for
    !     user-requested wavelength output. These are determined in CALC_AOD
    !     within RD_AOD and stored in IWVREQUIRED. The coefficients to
    !     interpolate from the LUT wavelengths to the user-requested
    !     waveelenths (in CALC_AOD) are used here.

    ! Select number of wavelengths required to loop over
    IF (ODSWITCH .EQ. 0) THEN !this is the call for Fast_JX at 1000nm
       NWVS   = 1
    ELSE
       IF ( LRAD ) THEN
          !Loop over all RT wavelengths (30)
          ! plus any required for calculating the AOD
          NWVS = State_Chm%Phot%NWVAA - State_Chm%Phot%NWVAA0 + &
                 State_Chm%Phot%NWVREQUIRED
       ELSE
          !Loop over wavelengths needed for
          !interpolation to those requested in geoschem_config.yml
          !(determined in RD_AOD)
          NWVS = State_Chm%Phot%NWVREQUIRED
       ENDIF
    ENDIF

    DO IIWV = 1, NWVS
       !now select the correct LUT wavelength
       IF (ODSWITCH .EQ. 0) THEN
          ! only doing for 1000nm (IWV1000 is set in RD_AOD)
          ! N.B. NWVS is fixed to 1 above - only one wavelength
          IWV=State_Chm%Phot%IWV1000
       ELSE
          IF ( LRAD ) THEN
             ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
             ! but add on any others required
             IF (IIWV.LE.30) THEN
                !index of RRTMG wavelengths starts after the standard NWVAA0
                !(currently NWVAA0=11, hard-coded in phot_container_mod based
                ! on the .dat LUT)
                IWV = IIWV + State_Chm%Phot%NWVAA0
             ELSE
                !now we calculate at wvs for the requested AOD
                IWV = IWVREQUIRED(IIWV-30)
             ENDIF
          ELSE
             ! IWVREQUIRED lists the index of requires standard wavelengths
             IWV = IWVREQUIRED(IIWV)
          ENDIF
       ENDIF

       ! Loop over types of aerosol with hygroscopic growth
       DO NA = 1, NAER

          ! Get ID following ordering of aerosol densities in RD_AOD
          IF ( NA <= NRHAER) THEN
             N = Map_NRHAER(NA)
          ELSE
             N = NA
          ENDIF

          !index for strat aerosol (only >0 for strat aero)
          ISTRAT=N-NRHAER

          ! NRT is subscript for RT arrays that contain SNA separately
          ! so their optics can be treated separately in future
          IF (N.GT.1) THEN
             NRT=N+2
          ELSE
             NRT=N
          ENDIF

          !==============================================================
          ! Determine aerosol growth rates from the relative
          ! humidity in each box
          !
          ! The optical depth scales with the radius and Q alone
          ! since SCALEDENS cancels as follows
          !
          !    SCALER 	= RW / RDRY
          !    SCALEDENS = DENSWET / DENSDRY
          !    SCALEM 	= SCALEDENS * SCALER**3
          !    SCALEOD 	= (SCALEQ * SCALEM) / (SCALEDENS * SCALER)
          !          	= SCALEQ * SCALER**2
          !
          ! Cap aerosol values at 90% relative humidity since
          ! aerosol growth at that point becomes highly nonlinear and
          ! relative humidities above this value essentially mean
          ! there is a cloud in that grid box
          !
          ! Q is the extinction efficiency
          !
          ! Each grid box (I,J,L) will fall into one of the RH bins,
          ! since each grid box will have a different RH value.  So,
          ! for SCALEOD(I,J,L,:), only one of the IRH bins will contain
          ! nonzero data, while the other IRH bins will all be zero.
          !==============================================================
          ! We loop over all selected wavelengths, referenced by IWV
          ! if FAST_J is calling then IWV will be for 1000nm only
          ! if RRTMG is on then IWV will be 30 wavelengths + AOD wavs,
          ! otherwise IWV will be at user input specified wavelengths

          ! Loop over grid boxes
          !$OMP PARALLEL DO                                                 &
          !$OMP PRIVATE( I,        J,       L,        R,        IRH       ) &
          !$OMP PRIVATE( RW,       QW,      AW,       SSW,      ASYW      ) &
          !$OMP PRIVATE( AW0,      QW0,     SSW0,     ASYW0,    REFF      ) &
          !$OMP PRIVATE( SCALEA,   SCALEQ,  SCALESSA, SCALEASY, FRAC      ) &
          !$OMP PRIVATE( SCALER,   SCALEOD, SCALEVOL, DRYAREA,  TAERVOL   ) &
          !$OMP PRIVATE( TK,       CONSEXP, VPRESH2O, RELHUM,   BCSCAT_AE ) &
#ifdef RRTMG
          !$OMP PRIVATE( IR                                               ) &
#endif
          !$OMP PRIVATE( RHOSTRAT, RAER,    SADSTRAT, XSASTRAT            ) &
          !$OMP PRIVATE( VDRY,     VH2O,    S,        g                   ) &
          !$OMP SCHEDULE( DYNAMIC, 8                                      ) &
          !$OMP COLLAPSE( 3                                               )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Zero private loop variables
             g         = 0
             IRH       = 0
             S         = 0
             RW        = 0.0_fp
             QW        = 0.0_fp
             AW        = 0.0_fp
             SSW       = 0.0_fp
             ASYW      = 0.0_fp 
             FRAC      = 0.0_fp
             AW0       = 0.0_fp
             QW0       = 0.0_fp
             SSW0      = 0.0_fp
             ASYW0     = 0.0_fp
             REFF      = 0.0_fp
             SCALEA    = 0.0_fp
             SCALEQ    = 0.0_fp
             SCALESSA  = 0.0_fp
             SCALEASY  = 0.0_fp
             FRAC      = 0.0_fp
             SCALER    = 0.0_fp
             SCALEOD   = 0.0_fp
             SCALEVOL  = 0.0_fp
             DRYAREA   = 0.0_fp
             TAERVOL   = 0.0_fp
             TK        = 0.0_fp
             CONSEXP   = 0.0_fp
             VPRESH2O  = 0.0_fp
             RELHUM    = 0.0_fp
             RHOSTRAT  = 0.0_fp
             RAER      = 0.0_fp
             SADSTRAT  = 0.0_fp
             XSASTRAT  = 0.0_fp
             VDRY      = 0.0_fp
             VH2O      = 0.0_fp
             BCSCAT_AE = 0.0_fp

             ! Loop over relative humidity bins
             IF (N == 1 .or. N == 3) THEN ! (hzhu, 08/2023)
                ! For SNA or Organics
                g = 1
                DO WHILE ( State_Chm%AerMass%PDER(I,J,L) > REAA(1,N,g) .and. &
                           g < State_Chm%Phot%NDRg )
                   ! REAA(1,N,g) is the upper limit of REFF
                   g = g + 1
                END DO

                IF (g == 1) THEN
                   DO R = 1, NRH
                      ! Wet radius in aerosol LUT files
                      RW(R) = REAA(R,N,g)

                      ! Extinction efficiency for Q for each RH bin
                      QW(R)   = QQAA(IWV,R,N,g)
                      AW(R)   = ALPHAA(IWV,R,N,g)
                      SSW(R)  = SSAA(IWV,R,N,g)
                      ASYW(R) = ASYMAA(IWV,R,N,g)
                   ENDDO

                ELSE
                   FRAC = (State_Chm%AerMass%PDER(I,J,L) - REAA(1,N,g-1))/  &
                          (REAA(1,N,g) - REAA(1,N,g-1))
                   IF ( FRAC > 1.0d0 ) FRAC = 1.0d0
                   DO R = 1, NRH
                      RW(R)  = FRAC*REAA(R,N,g) + (1.d0-FRAC)*REAA(R,N,g-1)

                      QW(R)  = FRAC*QQAA(IWV,R,N,g) + &
                              (1.d0-FRAC)*QQAA(IWV,R,N,g-1)

                      AW(R)  = FRAC*ALPHAA(IWV,R,N,g)+ &
                               (1.d0-FRAC)*ALPHAA(IWV,R,N,g-1)

                      SSW(R) = FRAC*SSAA(IWV,R,N,g) + &
                               (1.d0-FRAC)*SSAA(IWV,R,N,g-1)
                        
                      ASYW(R)= FRAC*ASYMAA(IWV,R,N,g)+ &
                               (1.d0-FRAC)*ASYMAA(IWV,R,N,g-1)
                   END DO
                END IF

             ELSE
               ! For other species
               DO R = 1, NRH
                  ! Wet radius in aerosol LUT files
                  RW(R) = REAA(R,N,State_Chm%Phot%DRg)

                  ! Extinction efficiency for Q for each RH bin
                  QW(R)   = QQAA(IWV,R,N,State_Chm%Phot%DRg)
                  AW(R)   = ALPHAA(IWV,R,N,State_Chm%Phot%DRg)
                  SSW(R)  = SSAA(IWV,R,N,State_Chm%Phot%DRg)
                  ASYW(R) = ASYMAA(IWV,R,N,State_Chm%Phot%DRg)
               ENDDO

            ENDIF

             ! Skip non-chemistry boxes
             IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

             ! Calculate RH. Not clear why the result of this calc is
             ! slightly different than State_Met%RH
             RELHUM   = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L)
             TK       = State_Met%T(I,J,L)
             CONSEXP  = 17.2693882e+0_fp * (TK - 273.16e+0_fp) / &
                        (TK - 35.86e+0_fp)
             VPRESH2O = CONSVAP * EXP(CONSEXP) / TK
             RELHUM   = RELHUM / VPRESH2O

             ! Sort into relative humidity bins
             ! Currently uses the 0, 50, 70, 80, 90% RH bins
             ! (95 and 99% also stored in data files but not used)
             IF (      RELHUM <= RH(2) ) THEN
                IRH = 1
             ELSE IF ( RELHUM <= RH(3) ) THEN
                IRH = 2
             ELSE IF ( RELHUM <= RH(4) ) THEN
                IRH = 3
             ELSE IF ( RELHUM <= RH(5) ) THEN
                IRH = 4
             ELSE
                IRH = 5
             ENDIF

             ! save the index of the relative humidity into array that is
             ! used by the photolysis to pull the right optics
             ! information from the FJX_spec.dat
             ! Previously the ODAER was organized by aerosol species and
             ! each RH bin, but this leaves a lot of the array redundant
             ! because only one RH bin is used at once
             ! This becomes a waste of memory with multiple wavelengths
             IRHARR(I,J,L) = IRH

             ! For the NRHth bin, we don't have to interpolate
             ! For the other bins, we have to interpolate
             ! Sometimes values can be zero, don't want to divide by 0...
             IF (AW(1).EQ.0.0) THEN
                AW0=1.0
             ELSE
                AW0=AW(1)
             ENDIF
             IF (QW(1).EQ.0.0) THEN
                QW0=1.0
             ELSE
                QW0=QW(1)
             ENDIF
             IF (SSW(1).EQ.0.0) THEN
                SSW0=1.0
             ELSE
                SSW0=SSW(1)
             ENDIF
             IF (ASYW(1).EQ.0.0) THEN
                ASYW0=1.0
             ELSE
                ASYW0=ASYW(1)
             ENDIF

             IF ( IRH == NRH ) THEN
                REFF     = RW(NRH)
                SCALEA   = AW(NRH)   / AW0  !QW(1) is dry extinction eff.
                SCALEQ   = QW(NRH)   / QW0
                SCALESSA = SSW(NRH)  / SSW0
                SCALEASY = ASYW(NRH) / ASYW0
             ELSE
                ! Interpolate between different RH
                FRAC = (RELHUM-RH(IRH)) / (RH(IRH+1)-RH(IRH))
                IF ( FRAC > 1.0d0 ) FRAC = 1.0d0
                REFF    = FRAC*RW(IRH+1)  + (1.d0-FRAC)*RW(IRH)
                SCALEA  =(FRAC*AW(IRH+1)  + (1.d0-FRAC)*AW(IRH))/AW0
                SCALEQ  =(FRAC*QW(IRH+1)  + (1.d0-FRAC)*QW(IRH))/QW0
                SCALESSA=(FRAC*SSW(IRH+1) + (1.d0-FRAC)*SSW(IRH))/SSW0
                SCALEASY=(FRAC*ASYW(IRH+1)+ (1.d0-FRAC)*ASYW(IRH))/ASYW0
             ENDIF

             SCALER  = REFF / RW(1)
             SCALEOD = SCALEQ * SCALER * SCALER


             IF ( N.LE.NRHAER ) THEN

                !--------------------------------------------------------
                ! %%%%%%% Aerosols undergoing hygroscopic growth %%%%%%%
                !--------------------------------------------------------

                !calculate optics for hyrdophillic aerosol here
                !However MDENS in LUT was in g/cm3 not kg/m3 so x1e3
                ODAER(I,J,L,IWV,N) = SCALEOD * BXHEIGHT(I,J,L) * 0.75d0 * &
                                     State_Chm%AerMass%WAERSL(I,J,L,N) * QW(1)   / &
                                     ( MSDENS(N) * RW(1) * 1.0D-6 )

                !Include BC absorption enhancement (xnw, 8/24/15)
                IF (N.eq.2) THEN

                   IF (LBCAE) THEN
                      BCSCAT_AE = ODAER(I,J,L,IWV,N)*SCALESSA*SSAA(IWV,1,N,State_Chm%Phot%DRg)
                      ODAER(I,J,L,IWV,N) = ODAER(I,J,L,IWV,N) * &
                                ( BCAE_1 + SCALESSA*SSAA(IWV,1,N,State_Chm%Phot%DRg) - &
                                  SCALESSA*SSAA(IWV,1,N,State_Chm%Phot%DRg)*BCAE_1 )

                      !now combine with hydrophilic OD as before
                      BCSCAT_AE = BCSCAT_AE + SSAA(IWV,1,N,State_Chm%Phot%DRg) * &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  State_Chm%AerMass%DAERSL(I,J,L,N-1) * QW(1)   / &
                                  ( MSDENS(N) * REAA(1,N,State_Chm%Phot%DRg) * 1.0D-6 )
                      ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                           (BCAE_2+SSAA(IWV,1,N,State_Chm%Phot%DRg) - SSAA(IWV,1,N,State_Chm%Phot%DRg)*BCAE_2) * &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  State_Chm%AerMass%DAERSL(I,J,L,N-1) * QW(1)   / &
                                  ( MSDENS(N) * REAA(1,N,State_Chm%Phot%DRg) * 1.0D-6 )

                   ELSE
                      !now combine with hydrophilic OD as before
                      ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  State_Chm%AerMass%DAERSL(I,J,L,N-1) * QW(1)   / &
                                  ( MSDENS(N) * REAA(1,N,State_Chm%Phot%DRg) * 1.0D-6 )
                   ENDIF

                ENDIF

                IF (N.eq.3) THEN
                   !now combine with hydrophilic OD as before
                   ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                                   0.75d0 * BXHEIGHT(I,J,L) * &
                                   State_Chm%AerMass%DAERSL(I,J,L,N-1) * QW(1)  / &
                                   ( MSDENS(N) * State_Chm%AerMass%PDER(I,J,L) * 1.0D-6 )
                ENDIF

                ! Get the AOD contribution from isoprene SOA only (eam, 2014)
                IF ( N == 3 .and. Is_ComplexSOA ) THEN
                   ISOPOD(I,J,L,IWV) = SCALEOD*BXHEIGHT(I,J,L)*0.75d0 &
                                   * State_Chm%AerMass%ISOAAQ(I,J,L) * QW(1)  / &
                                   ( MSDENS(N) * State_Chm%AerMass%PDER(I,J,L) * 1.0D-6 )
                ENDIF

             ELSE

                !--------------------------------------------------------
                !        %%%%%%% Stratospheric aerosols %%%%%%%
                !--------------------------------------------------------

                ! Get aerosol effective radius
                CALL GET_STRAT_OPT(State_Chm, I, J, L, ISTRAT, RAER, REFF, &
                                   SADSTRAT, XSASTRAT)

                ! SDE 2014-02-04
                ! The calculation used for the aerosols above
                ! is essentially a roundabout way of deriving
                ! the cross-sectional area. For log-normally
                ! distributed aerosols, this is much easier,
                ! and a direct query prevents the possibility
                ! of dividing a small mass by a small calculated
                ! radius and blowing up

                ! Aerosol optical depth
                ODAER(I,J,L,IWV,N) = BXHEIGHT(I,J,L) * XSASTRAT * QW(1) 

             ENDIF

#ifdef RRTMG
             !SNA currently treated as one with optics but considered
             !separately for RT, so we split them by mass here
             IF (N.EQ.1) THEN
                DO IR=1,3
                   RTODAER(I,J,L,IWV,N+IR-1)= ODAER(I,J,L,IWV,N)* &
                                              State_Chm%AerMass%FRAC_SNA(I,J,L,IR)
                   RTSSAER(I,J,L,IWV,N+IR-1)   = SCALESSA*SSAA(IWV,1,N,State_Chm%Phot%DRg)
                   RTASYMAER(I,J,L,IWV,N+IR-1) = SCALEASY*ASYMAA(IWV,1,N,State_Chm%Phot%DRg)
                ENDDO
             ELSE
                !RT arrays now offset from NAER by 2 (NRT=N+2 for N>1)
                !This will automatically be added after the standard aerosol
                !(NRHAER+1,2) but before dust
                RTODAER(I,J,L,IWV,NRT)     = ODAER(I,J,L,IWV,N)
                RTSSAER(I,J,L,IWV,NRT)     = SCALESSA*SSAA(IWV,1,N,State_Chm%Phot%DRg)
                !for BC SSA with absorption enhancement (xnw 8/24/15)
                IF ((N .EQ. 2) .AND. (LBCAE)) THEN
                   RTSSAER(I,J,L,IWV,NRT)  = BCSCAT_AE / &
                                             ODAER(I,J,L,IWV,N)
                ENDIF
                RTASYMAER(I,J,L,IWV,NRT)   = SCALEASY*ASYMAA(IWV,1,N,State_Chm%Phot%DRg)
             ENDIF
#endif

             ! Only need to do hyg once, not for each wavelength
             ! Skip stratospheric aerosols
             IF ((IIWV.EQ.1).and.(ISTRAT.le.0)) THEN

                !----------------------------------------------------
                ! Netcdf diagnostics computed here:
                !  Hygroscopic growth of SO4                [unitless]
                !  Hygroscopic growth of Black Carbon       [unitless]
                !  Hygroscopic growth of Organic Carbon     [unitless]
                !  Hygroscopic growth of Sea Salt (accum)   [unitless]
                !  Hygroscopic growth of Sea Salt (coarse)  [unitless]
                IF ( State_Diag%Archive_AerHygGrowth .AND. &
                     L <= State_Grid%MaxChemLev      .AND. &
                     ODSWITCH.EQ.1 ) THEN
                   S = State_Diag%Map_AerHygGrowth%id2slot(NA)
                   IF ( S > 0 ) THEN
                      State_Diag%AerHygGrowth(I,J,L,S) = SCALEOD
                   ENDIF
                ENDIF

                !=======================================================
                !now calulate the surface areas
                !==============================================================
                !  Calculate Aerosol Surface Area
                !
                !  Units ==> AERSL    [ kg aerosol m^-3 air ]
                !            MSDENS   [ kg aerosol m^-3 aerosol ]
                !            ERADIUS  [ cm      ]
                !            TAREA    [ cm^2 dry aerosol/cm^3 air ]
                !
                !  Note: first find volume of aerosol (cm^3 arsl/cm^3 air), then
                !        multiply by 3/radius to convert to surface area in cm^2
                !
                !  Wet Volume = AERSL * SCALER**3 / MSDENS
                !  Wet Surface Area = 3 * (Wet Volume) / ERADIUS
                !
                !  Effective radius for surface area and optical depths
                !  are identical.
                !==============================================================
                !========================================================
                ! NOTES:
                !    WAERSL   [ kg dry mass of wet aerosol m^-3 air ]
                !    ERADIUS  [ cm wet aerosol radius ]
                !    MSDENS   [ kg dry mass of aerosol m^-3 dry volume of aer]
                !    TAREA    [ cm^2 wet sfc area of aerosol cm^-3 air ]
                !    WTAREA   : same as TAREA, but excludes dry dust, BCPO, OCPO
                !               use same units as TAREA    (tmf, 4/18/07)
                !    WERADIUS : same as ERADIUS, but excludes dry dust,BCPO,OCPO
                !               use same units as ERADIUS  (tmf, 4/18/07)
                ! Wet dust WTAREA and WERADIUS are archived in dust_mod.F90.
                !========================================================

                !get scaling for R and VOL
                SCALER                 = REFF / RW(1)
                SCALEVOL               = SCALER**3
                ERADIUS(I,J,L,N+NDUST) = 1.0D-4 * REFF

                ! Store aerosol surface areas in TAREA, and be sure
                ! to list them following the dust surface areas
                TAREA(I,J,L,N+NDUST)   = 3.D0 * State_Chm%AerMass%WAERSL(I,J,L,N) * SCALEVOL / &
                                         ( ERADIUS(I,J,L,N+NDUST) * MSDENS(N) )

                WTAREA(I,J,L,N+NDUST)   = TAREA(I,J,L,N+NDUST)
                WERADIUS(I,J,L,N+NDUST) = ERADIUS(I,J,L,N+NDUST)
                ! For SO4-NIT-NH4-fine sea salt aerosol, re-calculate the wet
                ! effective
                ! radius using the water content from ISORROPIA/HETP.
                ! This new effective radius will be used for surface area
                ! used in heterogeneous chemistry. We don't use this
                ! effective radius in the optics above (OD, scattering,
                ! absorption) because the index of refraction, phase
                ! function, and Q must all be consistent with the radius and
                ! composition.
                ! (cdholmes, 5/17/2019, with update by XW, 5/28/2020)
                IF (N == 1) THEN

                   ! Volume of water, m3(H2O)/m3(air)
                   ! AeroH2O has units g/m3
                   VH2O = State_Chm%AeroH2O(I,J,L,NDUST+1) / 1e6
                   ! Volume of dry aerosol, m3(aerosol)/m3(air)
                   VDry = State_Chm%AerMass%WAERSL(I,J,L,1)/MSDENS(1) + State_Chm%AerMass%WAERSL(I,J,L,4)/MSDENS(4)

                   ! Notes on REFF derivation
                   ! Volume of wet aerosol: VWet = VDry + VH2O [note:
                   ! this is incorrect but has the correct limits for
                   ! VH2O/VDry << 1 and VH2O/VDry >> 1. It would be
                   ! better to use an empirical function for density.]
                   ! Volume of one dry particle v1dry = 4/3*pi*RDry**3
                   ! [note: RW(1) = RDry]
                   ! Number of aerosol particles: n = VDry / v1dry
                   ! Volume of wet aerosol is also: VWet = 4/3*pi * RWet**3 * n
                   ! So RWet = ( 3*VWet / (4 pi n) )**(1/3)
                   ! RWet = RDry * ( 1 + VH2O/Vdry )**(1/3)

                   ! Wet effective radius, um
                   ! Here assume the dry radius of the mixture = SNA
                   REFF = RW(1) * min( 3d0, &
                          ( 1d0 + safe_div( VH2O, VDry, 0d0 ) )**(1d0/3d0))

                   ACLRADIUS(I,J,L) = 1.0D-4 * REFF
                   ACLAREA(I,J,L) = 3.D0*(VH2O + VDry) / ACLRADIUS(I,J,L)
                ENDIF

                ! Save aerosol water content. Assume that the increase in volume
                ! equals the volume of pure water added, m3(H2O)/m3(air),
                ! then convert to g/m3
                ! Don't update SNA, keep ISORROPIA/HETP values
                IF (N.ne.1) THEN
                   State_Chm%AeroH2O(I,J,L,N+NDUST) = 1e+6_fp * &
                       State_Chm%AerMass%WAERSL(I,J,L,N) / MSDENS(N) * (ScaleVol - 1d0)
                ENDIF

                !include hydrophobic BC and OC
                !stored separate to hydrophillic in RT variables
                IF ((N.eq.2).or.(N.eq.3)) THEN

                   ! Dry surface area
                   ! SDE 2015-10-27: RW is in um, but everything
                   ! else is in terms of cm. Correct with 10^-4 factor
                   DRYAREA = 3.D0 * State_Chm%AerMass%DAERSL(I,J,L,N-1) / ( RW(1) * &
                             1.0D-4 * MSDENS(N) )

                   ! Add surface area to TAREA array
                   TAREA(I,J,L,N+NDUST) = WTAREA(I,J,L,N+NDUST) + DRYAREA

                   ! Define a new effective radius that accounts
                   ! for the hydrophobic aerosol
                   ERADIUS(I,J,L,NDUST+N) = ( WERADIUS(I,J,L,NDUST+N) *   &
                                              WTAREA(I,J,L,N+NDUST)   +   &
                                              RW(1) * 1.0D-4 * DRYAREA )/ &
                                            ( WTAREA(I,J,L,N+NDUST)   +   &
                                              DRYAREA )
                ENDIF !Hydrophobic aerosol surface area

                !----------------------------------------------------
                ! Netcdf diagnostics computed here:
                !  Aqueous aerosol volume
                IF ( State_Diag%Archive_AerAqVol .AND. N.EQ.1 ) THEN
                   State_Diag%AerAqVol(I,J,L) = &
                      ( ERADIUS(I,J,L,NDUST+N) * TAREA(I,J,L, N+NDUST) ) / 3.D0

                ENDIF

             ENDIF !Surface area calcs (1st wavelength only, non-strat)

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDDO !Loop over NAER
    ENDDO !End loop over NWVS

    !==============================================================
    ! Account for stratospheric aerosols (SDE 04/17/13)
    !==============================================================
    ! Loop over stratospheric aerosols
    DO ISTRAT = 1,NSTRATAER

       ! Index for combination of aerosol type and RH
       N = NRHAER + ISTRAT

       ! Assume ISTRAT == 1 is strat liquid aerosol and
       ! ISTRAT == 2 is polar strat cloud
       IsSLA = .FALSE.
       IsPSC = .FALSE.
       SELECT CASE ( ISTRAT )
       CASE ( 1 )
          IsSLA = .TRUE.
       CASE ( 2 )
          IsPSC = .TRUE.
       CASE DEFAULT
          WRITE( 6,'(a)') "WARNING: aerosol diagnostics not defined " // &
               "for NSTRATAER greater than 2!"
       END SELECT

       !$OMP PARALLEL DO                    &
       !$OMP DEFAULT( SHARED )              &
       !$OMP PRIVATE( I, J, L, RAER, REFF ) &
       !$OMP PRIVATE( SADSTRAT, XSASTRAT )  &
       !$OMP SCHEDULE( DYNAMIC )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Get aerosol effective radius
          CALL GET_STRAT_OPT(State_Chm, I,J,L,ISTRAT,RAER,REFF,SADSTRAT,XSASTRAT)

          ! Moved this from a separate loop for clarity
          IF ( State_Met%InChemGrid(I,J,L) ) THEN

             ! Add surface area to TAREA array
             TAREA(I,J,L,NDUST+N)  = SADSTRAT
             WTAREA(I,J,L,NDUST+N) = SADSTRAT

             ! Store radius
             ERADIUS(I,J,L,NDUST+N)  = RAER
             WERADIUS(I,J,L,NDUST+N) = RAER
          ENDIF

          !----------------------------------------------------
          ! Netcdf diagnostics computed here:
          !  Strat aerosol liquid surface area
          !  Polar strat cloud surface area
          IF ( State_Diag%Archive_AerSurfAreaSLA .AND. ISTRAT==1 ) THEN
             State_Diag%AerSurfAreaSLA(I,J,L) = SADSTRAT
          ENDIF
          IF ( State_Diag%Archive_AerSurfAreaPSC .AND. ISTRAT==2 ) THEN
             State_Diag%AerSurfAreaPSC(I,J,L) = SADSTRAT
          ENDIF

          !----------------------------------------------------
          ! Netcdf diagnostics computed here:
          !  Strat. liquid aerosol optical depth     [unitless]
          !  Strat. liquid aerosol number density    [#/cm3   ]
          !  Strat. particulate aerosol opt. depth   [unitless]
          !  Strat. particulate aerosol num. density [#/cm3   ]
          ! NOTE: Diagnostic output for each wavelength
          IF ( State_Diag%Archive_AODStrat .and. ODSWITCH .EQ. 1 ) THEN
             IF ( IsSLA ) THEN
                IF ( State_Diag%Archive_AerNumDenSLA ) THEN
                   State_Diag%AerNumDenSLA(I,J,L) = &
                        State_Chm%NDENS_AER(I,J,L,ISTRAT)*1.d-6
                ENDIF
                IF ( State_Diag%Archive_AODSLAWL1 ) THEN
                   State_Diag%AODSLAWL1(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,1),N)
                ENDIF
                IF ( State_Diag%Archive_AODSLAWL2 ) THEN
                   State_Diag%AODSLAWL2(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,2),N)
                ENDIF
                IF ( State_Diag%Archive_AODSLAWL3 ) THEN
                   State_Diag%AODSLAWL3(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,3),N)
                ENDIF
             ELSEIF ( IsPSC ) THEN
                IF ( State_Diag%Archive_AerNumDenPSC ) THEN
                   State_Diag%AerNumDenPSC(I,J,L) = &
                        State_Chm%NDENS_AER(I,J,L,ISTRAT)*1.d-6
                ENDIF
                IF ( State_Diag%Archive_AODPSCWL1 ) THEN
                   State_Diag%AODPSCWL1(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,1),N)
                ENDIF
                IF ( State_Diag%Archive_AODPSCWL2 ) THEN
                   State_Diag%AODPSCWL2(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,2),N)
                ENDIF
                IF ( State_Diag%Archive_AODPSCWL3 ) THEN
                   State_Diag%AODPSCWL3(I,J,L) = &
                        ODAER(I,J,L,IWVSELECT(1,3),N)
                ENDIF
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO ! end loop over stratospheric aerosols

    !==============================================================
    ! Non-stratospheric cloud diagnostics
    ! NOTE: The cloud optical depths are actually recorded at
    !       1000 nm, but vary little with wavelength.
    !==============================================================
    IF ( State_Diag%Archive_AOD .and. ODSWITCH .EQ. 1 ) THEN

       ! Loop over aerosol types (dust handled in dust_mod.F90)
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                               ) &
       !$OMP PRIVATE( I, J, L, N, W, LINTERP, IsWL1, IsWL2, IsWL3, S       ) &
       !$OMP SCHEDULE( DYNAMIC                                             )
       ! Loop over hydroscopic aerosols
       DO NA = 1, NRHAER

          ! Get ID following ordering of aerosol densities in RD_AOD
          N = Map_NRHAER(NA)

          ! Loop over wavelengths set in geoschem_config.yml
          DO W = 1, Input_Opt%NWVSELECT

             ! Set wavelength logical
             IsWL1 = .FALSE.
             IsWL2 = .FALSE.
             IsWL3 = .FALSE.
             SELECT CASE ( W )
             CASE ( 1 )
                IsWL1 = .TRUE.
             CASE ( 2 )
                IsWL2 = .TRUE.
             CASE ( 3 )
                IsWL3 = .TRUE.
             END SELECT

             ! no interpolation required if these are the same
             IF (IWVSELECT(1,W).EQ.IWVSELECT(2,W)) THEN
                LINTERP=.FALSE.
             ELSE
                LINTERP=.TRUE.
             ENDIF

             ! Loop over grid boxes
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX

                !----------------------------------------------------
                ! Netcdf diagnostics computed here:
                !  Sulfate Optical Depth (lambda1,2,3 nm)       [.]
                !  Black Carbon Optical Depth (lambda1,2,3 nm)  [.]
                !  Organic Carbon Optical Depth (lambda1,2,3 nm)[.]
                !  Sea Salt (accum) Opt Depth (lambda1,2,3 nm)  [.]
                !  Sea Salt (coarse) Opt Depth(lambda1,2,3 nm)  [.]
                IF ( .not. LINTERP ) THEN
                   IF ( State_Diag%Archive_AODHygWL1 .AND. IsWL1 ) THEN
                      S = State_Diag%Map_AODHygWL1%id2slot(NA)
                      IF ( S > 0 ) THEN
                         State_Diag%AODHygWL1(I,J,L,S) = &
                              ODAER(I,J,L,IWVSELECT(1,W),N)
                      ENDIF

                   ELSE IF ( State_Diag%Archive_AODHygWL2 .AND. IsWL2 ) THEN
                      S = State_Diag%Map_AODHygWL2%id2slot(NA)
                      IF ( S > 0 ) THEN
                         State_Diag%AODHygWL2(I,J,L,S) = &
                              ODAER(I,J,L,IWVSELECT(1,W),N)
                      ENDIF

                   ELSE IF ( State_Diag%Archive_AODHygWL3 .AND. IsWL3 ) THEN
                      S = State_Diag%Map_AODHygWL3%id2slot(NA)
                      IF ( S > 0 ) THEN
                         State_Diag%AODHygWL3(I,J,L,S) = &
                           ODAER(I,J,L,IWVSELECT(1,W),N)
                      ENDIF
                   ENDIF
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,W),N).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,W),N).GT.0)) THEN
                      IF ( State_Diag%Archive_AODHygWL1 .AND. IsWL1 ) THEN
                         S = State_Diag%Map_AODHygWL1%id2slot(NA)
                         IF ( S > 0 ) THEN
                            State_Diag%AODHygWL1(I,J,L,S) =                  &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**    &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/&
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
                         ENDIF

                      ELSE IF ( State_Diag%Archive_AODHygWL2 .AND. IsWL2 ) THEN
                         S = State_Diag%Map_AODHygWL2%id2slot(NA)
                         IF ( S > 0 ) THEN
                            State_Diag%AODHygWL2(I,J,L,S) =                  &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**    &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/&
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
                         ENDIF
                      ELSE IF ( State_Diag%Archive_AODHygWL3 .AND. IsWL3 ) THEN
                         S = State_Diag%Map_AODHygWL3%id2slot(NA)
                         IF ( S > 0 ) THEN
                            State_Diag%AODHygWL3(I,J,L,S) =                  &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**    &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/&
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF

                !----------------------------------------------------
                ! Netcdf diagnostics computed here:
                !  AOD for SOA from aq isoprene (lambda1,2,3 nm) [unitless]
                IF ( ( N == 3 ) .and. Is_ComplexSOA ) THEN
                   IF ( State_Diag%Archive_AODSOAfromAqIsopWL1 .AND. &
                        IsWL1 ) THEN
                      State_Diag%AODSOAfromAqIsopWL1(I,J,L) = &
                           ISOPOD(I,J,L,IWVSELECT(1,W))
                   ENDIF
                   IF ( State_Diag%Archive_AODSOAfromAqIsopWL2 .AND. &
                        IsWL2 ) THEN
                      State_Diag%AODSOAfromAqIsopWL2(I,J,L) = &
                           ISOPOD(I,J,L,IWVSELECT(1,W))
                   ENDIF
                   IF ( State_Diag%Archive_AODSOAfromAqIsopWL3 .AND. &
                        IsWL3 ) THEN
                      State_Diag%AODSOAfromAqIsopWL3(I,J,L) = &
                           ISOPOD(I,J,L,IWVSELECT(1,W))
                   ENDIF
                ENDIF

             ENDDO
             ENDDO
             ENDDO

          ENDDO ! end loop over wavelengths
       ENDDO ! end of loop over hygroscopic aerosols
       !$OMP END PARALLEL DO

    ENDIF

    !------------------------------------
    ! Aerosol Surface Areas
    !------------------------------------
    IF ( State_Diag%Archive_AerSurfAreaHyg .AND. ODSWITCH .EQ. 1) THEN

       !$OMP PARALLEL DO              &
       !$OMP DEFAULT( SHARED        ) &
       !$OMP PRIVATE( I, J, L, N, S ) &
       !$OMP SCHEDULE( DYNAMIC      )
       ! Loop over hydroscopic aerosols
       DO NA = 1, NRHAER

          ! Get ID following ordering of aerosol densities in RD_AOD
          N = Map_NRHAER(NA)

          !----------------------------------------------------
          ! Netcdf diagnostics computed here:
          !  Sulfate Surface Area                       [cm2/cm3]
          !  Black Carbon (hydrophilic) Surface Area    [cm2/cm3]
          !  Organic Carbon (hydrophilic) Surface Area  [cm2/cm3]
          !  Sea Salt (accum) Surface Area              [cm2/cm3]
          !  Sea Salt (coarse) Surface Area             [cm2/cm3]
          !----------------------------------------------------
          S = State_Diag%Map_AerSurfAreaHyg%id2slot(NA)
          IF ( S > 0 ) THEN
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX
                State_Diag%AerSurfAreaHyg(I,J,L,S) = &
                     TAREA(I,J,L,N+NDUST)
             ENDDO
             ENDDO
             ENDDO
          ENDIF

       ENDDO ! end of loop over hygroscopic aerosols
       !$OMP END PARALLEL DO

    ENDIF

    ! Turn off radiative effects of stratospheric aerosols?
    IF ( .not. LSTRATOD ) THEN
       ODAER(:,:,:,:,NRH+1) = 0.d0
       ODAER(:,:,:,:,NRH+2) = 0.d0
    ENDIF

    !=================================================================
    ! To turn off the radiative effects of different aerososl
    ! uncomment the following lines
    !=================================================================
    !DO R = 1,NRH
    !  ODAER(:,:,:,R)       = 0.d0  !sulfate
    !  ODAER(:,:,:,R+NRH)   = 0.d0  !BC
    !  ODAER(:,:,:,R+2*NRH) = 0.d0  !OC
    !  ODAER(:,:,:,R+3*NRH) = 0.d0  !SS(accum)
    !  ODAER(:,:,:,R+4*NRH) = 0.d0  !SS(coarse)
    !ENDDO
    !ODAER(:,:,:,NRHAER*NRH+1) = 0.d0   !SLA
    !ODAER(:,:,:,NRHAER*NRH+2) = 0.d0   !SPA

    !=================================================================
    ! To turn off heterogeneous chemistry on different aerosols
    ! uncomment the following lines
    !=================================================================
    !TAREA(:,NDUST+1) = 0.d0	!Sulfate
    !TAREA(:,NDUST+2) = 0.d0	!BC
    !TAREA(:,NDUST+3) = 0.d0	!OC
    !TAREA(:,NDUST+4) = 0.d0	!SS (accum)
    !TAREA(:,NDUST+5) = 0.d0	!SS (coarse)
    !TAREA(:,NDUST+NRHAER+1) = 0.d0 !SLA
    !TAREA(:,NDUST+NRHAER+2) = 0.d0 !SPA

    ! Free pointers
    IWVREQUIRED => NULL()
    IWVSELECT   => NULL()
    IRHARR      => NULL()
    ACOEF_WV    => NULL()
    BCOEF_WV    => NULL()
    REAA        => NULL()
    QQAA        => NULL()
    ALPHAA      => NULL()
    SSAA        => NULL()
    ASYMAA      => NULL()
    ISOPOD      => NULL()
    ODAER       => NULL()
#ifdef RRTMG
    RTODAER     => NULL()
    RTSSAER     => NULL()
    RTASYMAER   => NULL()
#endif
    BXHEIGHT    => NULL()
    ERADIUS     => NULL()
    TAREA       => NULL()
    WERADIUS    => NULL()
    WTAREA      => NULL()
    ACLRADIUS   => NULL()
    ACLAREA     => NULL()

    ! Reset first-time flag
    FIRST = .FALSE.

  END SUBROUTINE RDAER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_aerosol
!
! !DESCRIPTION: Subroutine INIT\_AEROSOL initializes module variables
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Aerosol( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,   ONLY : NAER, NDUST, NRHAER
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N, SpcID
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    TYPE(Species), POINTER :: SpcInfo
    
    !=================================================================
    ! INIT_AEROSOL begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_Aerosol (in module GeosCore/aerosol_mod.F90)'

    ! Exit immediately if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    ! Add tracer ID flags as module variables (bmy, 6/16/16)
    id_BCPI   = Ind_( 'BCPI'   )
    id_BCPO   = Ind_( 'BCPO'   )
    id_DST1   = Ind_( 'DST1'   )
    id_DST2   = Ind_( 'DST2'   )
    id_DST3   = Ind_( 'DST3'   )
    id_DST4   = Ind_( 'DST4'   )
    id_DUST01 = Ind_( 'DUST01' )
    id_NH4    = Ind_( 'NH4'    )
    id_NIT    = Ind_( 'NIT'    )
    id_OCPO   = Ind_( 'OCPO'   )
    id_OCPI   = Ind_( 'OCPI'   )
    id_SOAS   = Ind_( 'SOAS'   )
    id_SALA   = Ind_( 'SALA'   )
    id_SALC   = Ind_( 'SALC'   )
    id_SALACL = Ind_( 'SALACL' )
    id_SO4    = Ind_( 'SO4'    )
    id_SO4s   = Ind_( 'SO4s'   )
    id_HMS    = Ind_( 'HMS'    )
    id_NITs   = Ind_( 'NITs'   )
    id_POA1   = Ind_( 'POA1'   )
    id_POA2   = Ind_( 'POA2'   )
    id_OPOA1  = Ind_( 'OPOA1'  )
    id_OPOA2  = Ind_( 'OPOA2'  )
    id_TSOA1  = Ind_( 'TSOA1'  )
    id_TSOA2  = Ind_( 'TSOA2'  )
    id_TSOA3  = Ind_( 'TSOA3'  )
    id_TSOA0  = Ind_( 'TSOA0'  )
    id_ASOAN  = Ind_( 'ASOAN'  )
    id_ASOA1  = Ind_( 'ASOA1'  )
    id_ASOA2  = Ind_( 'ASOA2'  )
    id_ASOA3  = Ind_( 'ASOA3'  )
    id_SOAGX  = Ind_( 'SOAGX'  )
    id_SOAIE  = Ind_( 'SOAIE'  )
    id_INDIOL = Ind_( 'INDIOL' )
    id_LVOCOA = Ind_( 'LVOCOA' )

    ! Define logical flags
    IS_OCPI    = ( id_OCPI  > 0 )
    IS_OCPO    = ( id_OCPO  > 0 )
    IS_BC      = ( id_BCPI  > 0 .AND. id_BCPO  > 0 )
    IS_SO4     = ( id_SO4   > 0 )
    IS_HMS     = ( id_HMS   > 0 )
    IS_NH4     = ( id_NH4   > 0 )
    IS_NIT     = ( id_NIT   > 0 )
    IS_DST     = ( id_DST1  > 0 .AND. id_DST2  > 0 )
    IS_SAL     = ( id_SALA  > 0 .AND. id_SALC  > 0 )
    IS_POA     = ( id_POA1  > 0 .AND. id_POA2  > 0 )
    IS_OPOA    = ( id_OPOA1 > 0 .AND. id_OPOA2 > 0 )
    IS_TSOA    = ( id_TSOA1 > 0 .AND. id_TSOA2 > 0 .AND. &
                   id_TSOA3 > 0 .AND. id_TSOA0 > 0 )
    IS_ASOA    = ( id_ASOAN > 0 .AND. id_ASOA1 > 0 .AND. &
                   id_ASOA2 > 0 .AND. id_ASOA3 > 0 )
    IS_SOAGX   = ( id_SOAGX > 0 )
    Is_SimpleSOA  = ( id_SOAS > 0 )
    Is_ComplexSOA = Input_Opt%LSOA

    ! Set logicals also used for diagnostics
    IF ( IS_POA  .AND. State_Diag%Archive_AerMassPOA  ) State_Diag%isPOA  = .TRUE.
    IF ( IS_OPOA .AND. State_Diag%Archive_AerMassOPOA ) State_Diag%isOPOA = .TRUE.

    ! Initialize the mapping between hygroscopic species in the
    ! species database and the species order in NRHAER
    DO N = 1, NRHAER

       ! Get the species database index from the species database
       ! mapping array for hygroscopic growth species
       SpcID = State_Chm%Map_HygGrth(N)

       ! Point to the Species Database entry for species N
       SpcInfo => State_Chm%SpcData(SpcID)%Info

       ! Set the mapping to the ordering of aerosol densities in RD_AOD
       SELECT CASE ( TRIM(SpcInfo%Name) )
       CASE ( 'SO4' )
          Map_NRHAER(N) = 1
       CASE ( 'BCPI' )
          Map_NRHAER(N) = 2
       CASE ( 'OCPI', 'POA1' )
          Map_NRHAER(N) = 3
       CASE ( 'SALA' )
          Map_NRHAER(N) = 4
       CASE ( 'SALC' )
          Map_NRHAER(N) = 5
       CASE DEFAULT
          ErrMsg = 'WARNING: aerosol diagnostics not defined' // &
                   ' for NRHAER greater than 5!'
          CALL GC_ERROR( ErrMsg, RC, 'Init_Aerosol in aerosol_mod.F90' )
       END SELECT

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    !------------------------------------------------------------------------
    ! Read in AOD data
    !------------------------------------------------------------------------
    CALL RD_AOD( Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "RD_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF (Input_Opt%amIRoot) WRITE(6,*) 'Wavelength optics read successfully'

    !------------------------------------------------------------------------
    ! Compute the required wavelengths in the LUT to calculate requested AOD
    !------------------------------------------------------------------------
    CALL CALC_AOD( Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "CALC_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE INIT_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_AOD( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE InquireMod,    ONLY : FindFreeLUN
    USE State_Chm_Mod, ONLY : ChmState
#if defined( MODEL_CESM )
    USE UNITS,         ONLY : freeUnit
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The .dat files for each species contain the optical properties
!  at multiple wavelengths to be used in the online calculation of the aerosol
!  optical depth diagnostics.
!  These properties have been calculated using the same size and optical
!  properties as the FJX_spec.dat file used for the FAST-J photolysis
!  calculations (which is now redundant for aerosols, the values in the .dat
!  files here are now used). The file currently contains 11 wavelengths
!  for Fast-J and other commonly used wavelengths for satellite and
!  AERONET retrievals. 30 wavelengths follow that map onto RRTMG
!  wavebands for radiaitive flux calculations (not used if RRTMG is off).
!  A complete set of optical properties from 250-2000 nm for aerosols is
!  available at:
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!     -- David A. Ridley, 05/10/13 (update for new optics files)
!
! !REVISION HISTORY:
!  10 May 2010 - C. Heald      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER            :: I, J, K, N, g
    INTEGER            :: IOS, NJ1
    LOGICAL            :: LBRC, FileExists

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: THISFILE
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! String arrays
    CHARACTER(LEN=30)  :: SPECFIL(8)

    ! Pointers
    REAL*8, POINTER :: WVAA  (:,:)
    REAL*8, POINTER :: RHAA  (:,:)
    REAL*8, POINTER :: RDAA  (:,:,:)
    REAL*8, POINTER :: RWAA  (:,:,:)
    REAL*8, POINTER :: SGAA  (:,:)
    REAL*8, POINTER :: REAA  (:,:,:)
    REAL*8, POINTER :: NCMAA (:,:,:)
    REAL*8, POINTER :: NRLAA (:,:,:)
    REAL*8, POINTER :: QQAA  (:,:,:,:)
    REAL*8, POINTER :: ALPHAA(:,:,:,:)
    REAL*8, POINTER :: SSAA  (:,:,:,:)
    REAL*8, POINTER :: ASYMAA(:,:,:,:)
    REAL*8, POINTER :: PHAA  (:,:,:,:,:)

    !================================================================
    ! RD_AOD begins here!
    !================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at RD_AOD (in module GeosCore/photolysis_mod.F90)'
    LBRC     = Input_Opt%LBRC
    DATA_DIR = TRIM( Input_Opt%AER_OPTICS_DIR )

    ! Set Pointers
    WVAA   => State_Chm%Phot%WVAA
    RHAA   => State_Chm%Phot%RHAA
    RDAA   => State_Chm%Phot%RDAA
    RWAA   => State_Chm%Phot%RWAA
    SGAA   => State_Chm%Phot%SGAA
    REAA   => State_Chm%Phot%REAA
    NRLAA  => State_Chm%Phot%NRLAA
    NCMAA  => State_Chm%Phot%NCMAA
    QQAA   => State_Chm%Phot%QQAA
    ALPHAA => State_Chm%Phot%ALPHAA
    SSAA   => State_Chm%Phot%SSAA
    ASYMAA => State_Chm%Phot%ASYMAA
    PHAA   => State_Chm%Phot%PHAA

    ! Get a free LUN
    NJ1 = findFreeLUN()

    ! IMPORTANT: aerosol_mod.F and dust_mod.F expect aerosols in this order
    !
    ! Treating strat sulfate with GADS data but modified to match
    ! the old Fast-J values size (r=0.09um, sg=0.6) - I think there's
    ! evidence that this is too smale and narrow e.g. Deshler et al. 2003
    ! NAT should really be associated with something like cirrus cloud
    ! but for now we are just treating the NAT like the sulfate... limited
    ! info but ref index is similar e.g. Scarchilli et al. (2005)
    !(DAR 05/2015)
    DATA SPECFIL /"so4.dat","soot.dat","org.dat", &
                  "ssa.dat","ssc.dat",            &
                  "h2so4.dat","h2so4.dat",        &
                  "dust.dat"/

    ! Loop over the array of filenames
    DO k = 1, State_Chm%Phot%NSPAA

       ! Choose different set of input files for standard (trop+strat chenm)
       ! and tropchem (trop-only chem) simulations
       THISFILE = TRIM( DATA_DIR ) // '/' // TRIM( SPECFIL(k) )

       !--------------------------------------------------------------
       ! In dry-run mode, print file path to dryrun log and cycle.
       ! Otherwise, print file path to stdout and continue.
       !--------------------------------------------------------------

       ! Test if the file exists
       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )

       ! Test if the file exists and define an output string
       IF ( FileExists ) THEN
          FileMsg = 'PHOTOLYSIS (RD_AOD): Opening'
       ELSE
          FileMsg = 'PHOTOLYSIS (RD_AOD): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write to stdout for both regular and dry-run simulations
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
300       FORMAT( a, ' ', a )
       ENDIF

       ! For dry-run simulations, cycle to next file.
       ! For regular simulations, throw an error if we can't find the file.
       IF ( Input_Opt%DryRun ) THEN
          CYCLE
       ELSE
          IF ( .not. FileExists ) THEN
             WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! If not a dry-run, read data from each species file
       !--------------------------------------------------------------

       ! Open file
       OPEN( NJ1, FILE=TRIM( THISFILE ), STATUS='OLD', IOSTAT=RC )

       ! Error check
       IF ( RC /= 0 ) THEN
          ErrMsg = 'Error opening file: ' // TRIM( ThisFile )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read header lines
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       ! Second header line added for more info
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       READ(  NJ1, '(A)' ) TITLE0
110    FORMAT( 3x, a20 )

       IF (k == 1 .OR. k == 3) THEN
       ! for SO4 and ORGANICS, dry aerosol size varies, therefore all
       ! opt properties vary.
       DO g = 1, State_Chm%Phot%NDRg
       DO i = 1, State_Chm%Phot%NRAA
       DO j = 1, State_Chm%Phot%NWVAA

          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
                      RDAA(i,k,g),RWAA(i,k,g),SGAA(i,k),QQAA(j,i,k,g),   &
                      ALPHAA(j,i,k,g),REAA(i,k,g),SSAA(j,i,k,g),         &
                      ASYMAA(j,i,k,g),(PHAA(j,i,k,n,g),n=1,8)

          ! make note of where 1000nm is for FAST-J calcs
          IF (WVAA(j,k).EQ.1000.0) State_Chm%Phot%IWV1000=J

       ENDDO
       ENDDO
       ENDDO

       ELSE
       ! For other species, keep g = default Rg (DRg)
       g = State_Chm%Phot%DRg
       DO i = 1, State_Chm%Phot%NRAA
       DO j = 1, State_Chm%Phot%NWVAA

          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
                      RDAA(i,k,g),RWAA(i,k,g),SGAA(i,k),QQAA(j,i,k,g),   &
                      ALPHAA(j,i,k,g),REAA(i,k,g),SSAA(j,i,k,g),         &
                      ASYMAA(j,i,k,g),(PHAA(j,i,k,n,g),n=1,8)

          ! make note of where 1000nm is for FAST-J calcs
          IF (WVAA(j,k).EQ.1000.0) State_Chm%Phot%IWV1000=J

       ENDDO
       ENDDO

       ENDIF

       ! Close file
       CLOSE( NJ1 )

    ENDDO

#if defined( MODEL_CESM )
   CALL freeUnit(NJ1)
#endif

  ! Free pointers
    WVAA   => NULL()
    RHAA   => NULL()
    RDAA   => NULL()
    RWAA   => NULL()
    SGAA   => NULL()
    REAA   => NULL()
    NCMAA  => NULL()
    NRLAA  => NULL()
    QQAA   => NULL()
    ALPHAA => NULL()
    SSAA   => NULL()
    ASYMAA => NULL()
    PHAA   => NULL()

  END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_aod
!
! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
! in the optics LUT wavelengths and the coefficients required to
! calculate the angstrom exponent for interpolating optics to the requested
! wavelength. If the wavelength requested matches a standard wavelength
! in the LUT then we skip the interpolation (DAR 09/2013)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_AOD( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
#ifdef RRTMG
    USE PARRRTM,       ONLY : NBNDLW
#endif
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: RC
!
! !REMARKS:
!  Now the user is able to select any 3 wavelengths for optics
!  output in the geoschem_config.yml file we need to be able to interpolate
!  to those wavelengths based on what is available in the optics
!  look-up table.
!                                                                             .
!  The standard lookup table currently has values for
!  11 common wavelengths followed by 30 that are required by RRTMG.
!  Only those required to interpolate to user requested
!  wavelengths are selected from the standard wavelengths. RRTMG
!  wavelengths are not used in the interpolation for AOD output
!  (DAR 10/2013)
!                                                                             .
!   UPDATE: because the RT optics output doesnt have access to the
!   standard wavelengths we now calculate two sets of values: one
!   for the ND21 and diag3 outputs that use the standard wavelengths
!   and one for RRTMG diagnostics that interpolate the optics from RRTMG
!   wavelengths. Perhaps a switch needs adding to switch off the RT
!   optics output (and interpolation) if this ends up costing too
!   much and is not used, but it is ideal to have an optics output
!   that matches exactly what RRTMG uses to calculate the fluxes
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
    INTEGER             :: NWVAA, NWVAA0, NWVREQUIRED, NRTWVREQUIRED
    REAL(fp)            :: WVDIF

    ! Pointers
    INTEGER, POINTER :: IWVREQUIRED  (:)
    INTEGER, POINTER :: IRTWVREQUIRED(:)
    INTEGER, POINTER :: IWVSELECT    (:,:)
    INTEGER, POINTER :: IRTWVSELECT  (:,:)
    REAL*8,  POINTER :: ACOEF_WV     (:)
    REAL*8,  POINTER :: BCOEF_WV     (:)
    REAL*8,  POINTER :: CCOEF_WV     (:)
    REAL*8,  POINTER :: ACOEF_RTWV   (:)
    REAL*8,  POINTER :: BCOEF_RTWV   (:)
    REAL*8,  POINTER :: CCOEF_RTWV   (:)
    REAL*8,  POINTER :: WVAA         (:,:)

    !================================================================
    ! CALC_AOD begins here!
    !================================================================

    ! Constants State_Chm%Phot
    NWVAA         = State_Chm%Phot%NWVAA
    NWVAA0        = State_Chm%Phot%NWVAA0

    ! Scalars in State_Chm%Phot that will be set in this subroutine
    NWVREQUIRED   = State_Chm%Phot%NWVREQUIRED
    NRTWVREQUIRED = State_Chm%Phot%NRTWVREQUIRED

    ! Set pointers
    IWVREQUIRED   => State_Chm%Phot%IWVREQUIRED
    IRTWVREQUIRED => State_Chm%Phot%IRTWVREQUIRED
    IWVSELECT     => State_Chm%Phot%IWVSELECT
    IRTWVSELECT   => State_Chm%Phot%IRTWVSELECT
    ACOEF_WV      => State_Chm%Phot%ACOEF_WV
    BCOEF_WV      => State_Chm%Phot%BCOEF_WV
    CCOEF_WV      => State_Chm%Phot%CCOEF_WV
    ACOEF_RTWV    => State_Chm%Phot%ACOEF_RTWV
    BCOEF_RTWV    => State_Chm%Phot%BCOEF_RTWV
    CCOEF_RTWV    => State_Chm%Phot%CCOEF_RTWV
    WVAA          => State_Chm%Phot%WVAA

    !cycle over standard wavelengths
    N0=1
    N1=NWVAA0
    NSTEP=1
    NWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP ! 1 to 11
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested wavelength
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(1,1)
          IWVSELECT(1,W)=1
          IWVSELECT(2,W)=1 !300nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0,1)
          IWVSELECT(1,W)=NWVAA0
          IWVSELECT(2,W)=NWVAA0 !1020nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
          ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/ &
                                  WVAA(IWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_WV(W) =(Input_Opt%WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/ &
                      (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
       ENDIF
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
          !write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
          !write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
          !write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#ifdef RRTMG
    !repeat for RRTMG wavelengths to get the closest wavelength
    !indices and the interpolation coefficients
    !Indices are relative to all wavelengths in the LUT i.e. the RRTMG
    !wavelengths start at NWVAA0+1
    N0=NWVAA0+1
    N1=NWVAA
    NSTEP=1
    NRTWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IRTWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IRTWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IRTWVSELECT(2,W).EQ.IRTWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested
          !wavelength
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA-1,1)
          IRTWVSELECT(1,W)=NWVAA-1
          IRTWVSELECT(2,W)=NWVAA-1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0+1,1)
          IRTWVSELECT(1,W)=NWVAA0+1
          IRTWVSELECT(2,W)=NWVAA0+1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IRTWVSELECT(1,W).NE.IRTWVSELECT(2,W)) THEN
          ACOEF_RTWV(W) = WVAA(IRTWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_RTWV(W) =1.0d0/(LOG(WVAA(IRTWVSELECT(2,W),1)/ &
                                    WVAA(IRTWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_RTWV(W) =(Input_Opt%WVSELECT(W)-WVAA(IRTWVSELECT(1,W),1))/ &
                      (WVAA(IRTWVSELECT(2,W),1)-WVAA(IRTWVSELECT(1,W),1))
       ENDIF
       !convert wavelength index to that required by rrtmg_rad_transfer
       !i.e. without the standard and LW wavelengths
       IRTWVSELECT(1,W) = IRTWVSELECT(1,W) - NWVAA0 - NBNDLW
       IRTWVSELECT(2,W) = IRTWVSELECT(2,W) - NWVAA0 - NBNDLW
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N RT WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'RT WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'RT WAVELENGTH REQUIRED:', NRTWVREQUIRED
          write(6,*) IRTWVSELECT(1,W),WVAA(IRTWVSELECT(1,W)+NWVAA0+NBNDLW,1)
          write(6,*) IRTWVSELECT(2,W),WVAA(IRTWVSELECT(2,W)+NWVAA0+NBNDLW,1)
          write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#endif

    ! Copy values back into State_Chm
    State_Chm%Phot%NWVREQUIRED   = NWVREQUIRED
    State_Chm%Phot%NRTWVREQUIRED = NRTWVREQUIRED

    ! Free pointers
    IWVREQUIRED   => NULL()
    IRTWVREQUIRED => NULL()
    IWVSELECT     => NULL()
    IRTWVSELECT   => NULL()
    ACOEF_WV      => NULL()
    BCOEF_WV      => NULL()
    CCOEF_WV      => NULL()
    ACOEF_RTWV    => NULL()
    BCOEF_RTWV    => NULL()
    CCOEF_RTWV    => NULL()
    WVAA          => NULL()

  END SUBROUTINE CALC_AOD
!EOC
END MODULE AEROSOL_MOD
