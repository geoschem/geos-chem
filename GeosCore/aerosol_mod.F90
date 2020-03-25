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
  PUBLIC :: AEROSOL_CONC
  PUBLIC :: CLEANUP_AEROSOL
  PUBLIC :: INIT_AEROSOL
  PUBLIC :: RDAER
  PUBLIC :: Set_AerMass_Diagnostic
!
! !PUBLIC DATA MEMBERS:
!
  !========================================================================
  ! BCPI        : Hydrophilic black carbon aerosol   [kg/m3]
  ! BCPO        : Hydrophobic black carbon aerosol   [kg/m3]
  ! OCPI        : Hydrophilic organic carbon aerosol [kg/m3]
  ! OCPO        : Hydrophobic organic carbon aerosol [kg/m3]
  ! OCPISOA     : Hydrophilic OC + SOA aerosol       [kg/m3]
  ! SALA        : Accumulation mode seasalt aerosol  [kg/m3]
  ! SALC        : Coarse mode seasalt aerosol        [kg/m3]
  ! SO4_NH4_NIT : Lumped SO4-NH4-NIT aerosol         [kg/m3]
  ! SO4         : Sulfate aerosol                    [kg/m3]
  ! NH4         : Ammonium aerosol                   [kg/m3]
  ! NIT         : Inorganic nitrate aerosol          [kg/m3]
  ! SOILDUST    : Mineral dust aerosol from soils    [kg/m3]
  ! SLA         : Stratospheric liquid aerosol       [kg/m3]
  ! SPA         : Stratospheric particulate aerosol  [kg/m3]
  ! TSOA        : Terpene SOA                        [kg/m3]
  ! ASOA        : Aromatic + IVOC SOA                [kg/m3]
  ! OPOA        : Aerosol product of SVOC oxidation  [kg/m3]
  ! SOAGX       : SOA product of GLYX                [kg/m3]
  ! SOAIE       : SOA product of IEPOX & HMML        [kg/m3]
  ! PM25        : Particulate matter < 2.5 um        [kg/m3]
  ! ISOAAQ      : Isoprene SOA (aqueous formation)   [kg/m3]
  ! SOAS        : Simple SOA                         [kg/m3]
  ! OCFPOA      : OM/OC for POA                      [unitless]
  ! OCFOPOA     : OM/OC for OPOA, OCPI, OCPO         [unitless]
  !========================================================================
  REAL(fp), ALLOCATABLE, PUBLIC :: BCPI(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: BCPO(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OCPI(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OCPO(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OCPISOA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SALA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SALC(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SO4_NH4_NIT(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SO4(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: NH4(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: NIT(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: FRAC_SNA(:,:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SOILDUST(:,:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SLA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SPA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: TSOA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: ASOA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OPOA(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SOAGX(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: PM25(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: ISOAAQ(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: SOAS(:,:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OCFPOA(:,:)
  REAL(fp), ALLOCATABLE, PUBLIC :: OCFOPOA(:,:)

  ! Growth factors
  REAL(fp),              PUBLIC :: SIA_GROWTH
  REAL(fp),              PUBLIC :: ORG_GROWTH
  REAL(fp),              PUBLIC :: SSA_GROWTH
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
  ! Mass of hydrophobic aerosol from Mian Chin
  REAL(fp), ALLOCATABLE, SAVE   :: DAERSL(:,:,:,:)

  ! Mass of hydrophilic aerosol from Mian Chin
  REAL(fp), ALLOCATABLE, SAVE   :: WAERSL(:,:,:,:)

  ! Add tracer ID flags as module variables (bmy, 6/16/16)
  INTEGER :: id_BCPI,  id_BCPO,  id_DST1,  id_DST2
  INTEGER :: id_DST3,  id_DST4,  id_NH4,   id_NIT
  INTEGER :: id_OCPO,  id_OCPI,  id_SALA,  id_SALC
  INTEGER :: id_SO4,   id_SO4s,  id_NITs
  INTEGER :: id_POA1,  id_POA2,  id_OPOA1, id_OPOA2
  INTEGER :: id_TSOA1, id_TSOA2, id_TSOA3, id_TSOA0
  INTEGER :: id_ASOAN, id_ASOA1, id_ASOA2, id_ASOA3
  INTEGER :: id_DUST1, id_SOAS
  INTEGER :: id_SOAGX, id_SOAIE
  INTEGER :: id_INDIOL,id_LVOCOA

  ! Index to map between NRHAER and species database hygroscopic species
  ! NOTE: Increasing value of NRHAER in CMN_SIZE_Mod.F90 (e.g. if there is
  ! a new hygroscopic species) requires manual update of this mapping
  ! (ewl, 1/23/17)
  INTEGER :: Map_NRHAER(5)

  ! Diagnostic switches
  LOGICAL :: Is_POA

  ! Conversionf factors to ugC/m3 for Total Organic Carbon diagnostic
  REAL(fp) :: Fac_INDIOL
  REAL(fp) :: Fac_LVOCOA
  REAL(fp) :: Fac_SOAGX
  REAL(fp) :: Fac_SOAIE

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
    USE CMN_FJX_MOD,       ONLY : REAA
    USE ErrCode_Mod
    USE ERROR_MOD
    USE HCO_EMISLIST_MOD,  ONLY : HCO_GetPtr
    USE HCO_Error_Mod
    USE HCO_INTERFACE_MOD, ONLY : HcoState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE UCX_MOD,           ONLY : KG_STRAT_AER
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units
    USE TIME_MOD,          ONLY : GET_MONTH
#ifdef TOMAS
    USE TOMAS_MOD,         ONLY : IBINS
#endif
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
    INTEGER             :: I, J, L, N, NA, ND, K
    INTEGER             :: k_SO4
    INTEGER             :: k_ORG
    INTEGER             :: k_SSA
    REAL(fp)            :: Rad_wet, Rad_dry
    REAL(fp)            :: Rho_wet, Rho_dry
    REAL(fp)            :: REFF

    ! Logical flags
    LOGICAL             :: prtDebug
    LOGICAL             :: LCARB
    LOGICAL             :: LDUST
    LOGICAL             :: LSSALT
    LOGICAL             :: LSULF
    LOGICAL             :: LUCX
    LOGICAL             :: IS_OCPO,  IS_OCPI,  IS_BC
    LOGICAL             :: IS_SO4,   IS_NH4,   IS_NIT
    LOGICAL             :: IS_SAL,   IS_DST
    LOGICAL             :: IS_TSOA,  IS_ASOA
    LOGICAL             :: IS_POA,   IS_OPOA
    LOGICAL             :: IS_SOAGX
    LOGICAL             :: Is_SimpleSOA
    LOGICAL             :: Is_ComplexSOA

    ! Pointers
    REAL(fp), POINTER   :: Spc(:,:,:,:)
    REAL(fp), POINTER   :: AIRVOL(:,:,:)
    REAL(fp), POINTER   :: PMID(:,:,:)
    REAL(fp), POINTER   :: T(:,:,:)

    ! Other variables
    CHARACTER(LEN=63)   :: OrigUnit

    ! For spatially and seasonally varying OM/OC
    CHARACTER(LEN=255)  :: FIELDNAME
    INTEGER             :: MONTH
    LOGICAL             :: FND
    REAL(f4), POINTER   :: OMOC(:,:) => NULL()

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
    LUCX    = Input_Opt%LUCX

    ! Do we have to print debug output?
    prtDebug   = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Define logical flags
    IS_OCPI    = ( id_OCPI  > 0 )
    IS_OCPO    = ( id_OCPO  > 0 )
    IS_BC      = ( id_BCPI  > 0 .AND. id_BCPO  > 0 )
    IS_SO4     = ( id_SO4   > 0 )
    IS_NH4     = ( id_NH4   > 0 )
    IS_NIT     = ( id_NIT   > 0 )
    IS_DST     = ( id_DST1  > 0 .AND. id_DST2  > 0 )
    IS_SAL     = ( id_SALA  > 0 .AND. id_SALC  > 0 )
    IS_TSOA    = ( id_TSOA1 > 0 .AND. id_TSOA2 > 0 .AND. &
                   id_TSOA3 > 0 .AND. id_TSOA0 > 0 )
    IS_ASOA    = ( id_ASOAN > 0 .AND. id_ASOA1 > 0 .AND. &
                   id_ASOA2 > 0 .AND. id_ASOA3 > 0 )
    IS_POA     = ( id_POA1  > 0 .AND. id_POA2  > 0 )
    IS_OPOA    = ( id_OPOA1 > 0 .AND. id_OPOA2 > 0 )
    IS_SOAGX   = ( id_SOAGX > 0 )

    ! Logical flags for SOA scheme
    Is_SimpleSOA  = ( id_SOAS > 0 )
    Is_ComplexSOA = Input_Opt%LSOA

    ! Convert species to [kg] for this routine
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error at start of AEROSOL_CONC!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize pointers
    Spc    => State_Chm%Species
    AIRVOL => State_Met%AIRVOL
    PMID   => State_Met%PMID
    T      => State_Met%T

    !=================================================================
    ! OM/OC ratio
    !
    ! Get spatial and seasonally varying OM/OC from Philip et al. (2014)
    ! or use default global mean values recommended by Aerosols WG
    !=================================================================

    ! Attenot to get pointer to OM/OC for current month from HEMCO
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
    CALL HCO_GetPtr( HcoState, FIELDNAME, OMOC, RC, FOUND=FND )

    IF ( RC == HCO_SUCCESS .AND. FND ) THEN

       ! Set OM/OC using spatially and seasonally varying data from
       ! Philip et al. (2014)
       OCFPOA(:,:)  = OMOC(:,:) ! OM/OC for POA
       OCFOPOA(:,:) = OMOC(:,:) ! OM/OC for OPOA, OCPI, and OCPO

    ELSE

       ! Use default global mean OM/OC recommended by the Aerosols WG
       OCFPOA(:,:)  = 1.4e+0_fp ! OM/OC for POA
       OCFOPOA(:,:) = 2.1e+0_fp ! OM/OC for OPOA, OCPI, and OCPO

    ENDIF

    ! Save OM/OC
    State_Chm%OMOC_POA(:,:) = OCFPOA(:,:)
    State_Chm%OMOC_OPOA(:,:) = OCFOPOA(:,:)

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
       Rad_dry    = REAA(1,k_SO4)
       Rad_wet    = REAA(1,k_SO4) + 35e+0_fp * &
                  ( REAA(2,k_SO4) - REAA(1,k_SO4) ) / 50e+0_fp
       Rho_dry    = State_Chm%SpcData(id_SO4)%Info%Density
       SIA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Growth factor for OCPI + SOA
       Rad_dry    = REAA(1,k_ORG)
       Rad_wet    = REAA(1,k_ORG) + 35e+0_fp * &
                  ( REAA(2,k_ORG) - REAA(1,k_ORG) ) / 50e+0_fp
       IF ( IS_POA ) THEN
          Rho_dry    = State_Chm%SpcData(id_POA1)%Info%Density
       ELSE IF ( IS_OCPI ) THEN
          Rho_dry    = State_Chm%SpcData(id_OCPI)%Info%Density
       ENDIF
       ORG_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Growth factor for SALA
       Rad_dry    = REAA(1,k_SSA)
       Rad_wet    = REAA(1,k_SSA) + 35e+0_fp * &
                  ( REAA(2,k_SSA) - REAA(1,k_SSA) ) / 50e+0_fp
       Rho_dry    = State_Chm%SpcData(id_SALA)%Info%Density
       SSA_GROWTH = 1 + ( ( ( Rad_wet / Rad_dry ) ** 3 - 1 ) * &
                            ( Rho_wet / Rho_dry ) )

       ! Print debug info
       IF ( prtDebug ) THEN
          WRITE( 6,'(a)') 'Growth factors at 35% RH:'
          WRITE( 6, 100 ) SIA_GROWTH, ' for SO4, NIT, and NH4'
          WRITE( 6, 100 ) ORG_GROWTH, ' for OCPI and SOA'
          WRITE( 6, 100 ) SSA_GROWTH, ' for SALA'
100       FORMAT(F5.2,A)
       ENDIF

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

          IF ( LUCX ) THEN

             !--------------------------------------------------------
             !          %%%%%%% UCX-based mechanisms %%%%%%%
             !--------------------------------------------------------

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
                SO4_NH4_NIT(I,J,L) = ( Spc(I,J,L,id_SO4)    + &
                                       Spc(I,J,L,id_NH4)    + &
                                       Spc(I,J,L,id_NIT)    + &
                                       Spc(I,J,L,id_SO4s)   + &
                                       Spc(I,J,L,id_NITs) ) / &
                                       AIRVOL(I,J,L)
                SO4(I,J,L) = Spc(I,J,L,id_SO4) / AIRVOL(I,J,L)
                NH4(I,J,L) = Spc(I,J,L,id_NH4) / AIRVOL(I,J,L)
                NIT(I,J,L) = Spc(I,J,L,id_NIT) / AIRVOL(I,J,L)
                SLA(I,J,L) = 0e+0_fp
                SPA(I,J,L) = 0e+0_fp

             ELSE

                ! Tropospheric sulfate is zero in stratosphere
                SO4_NH4_NIT(I,J,L) = 0e+0_fp
                SO4(I,J,L) = 0e+0_fp
                NH4(I,J,L) = 0e+0_fp
                NIT(I,J,L) = 0e+0_fp
                SLA(I,J,L) = KG_STRAT_AER(I,J,L,1) / AIRVOL(I,J,L)
                SPA(I,J,L) = KG_STRAT_AER(I,J,L,2) / AIRVOL(I,J,L)

             ENDIF

          ELSE

             !--------------------------------------------------------
             !        %%%%%%% Tropchem only mechanism %%%%%%%
             !--------------------------------------------------------

             ! Compute SO4 aerosol concentration [kg/m3]
             ! now including sulfate and nitrate associated with sea-salt
             ! NOTE: these should be treated as having a sea-salt size
             ! distribution but are currently simply treated in the same
             ! way (size and optics) as all other sulfate aerosol (DAR
             ! 2013)
             SO4_NH4_NIT(I,J,L) = ( Spc(I,J,L,id_SO4)    + &
                                    Spc(I,J,L,id_NH4)    + &
                                    Spc(I,J,L,id_NIT)    + &
                                    Spc(I,J,L,id_SO4s)   + &
                                    Spc(I,J,L,id_NITs) ) / &
                                    AIRVOL(I,J,L)
             SO4(I,J,L) = Spc(I,J,L,id_SO4) / AIRVOL(I,J,L)
             NH4(I,J,L) = Spc(I,J,L,id_NH4) / AIRVOL(I,J,L)
             NIT(I,J,L) = Spc(I,J,L,id_NIT) / AIRVOL(I,J,L)

          ENDIF

          ! Add error check for safe division (bmy, 4/7/15)
          IF ( SO4_NH4_NIT(I,J,L) > 0e+0_fp ) THEN

             ! Save these fractions for partitioning of optics
             ! until later when these may be treated independently
             FRAC_SNA(I,J,L,1) = ( ( Spc(I,J,L,id_SO4 ) + Spc(I,J,L,id_SO4s) ) &
                                 / AIRVOL(I,J,L) ) / SO4_NH4_NIT(I,J,L)

             FRAC_SNA(I,J,L,2) = ( ( Spc(I,J,L,id_NIT ) + Spc(I,J,L,id_NITs) ) &
                                 / AIRVOL(I,J,L) ) / SO4_NH4_NIT(I,J,L)

             FRAC_SNA(I,J,L,3) = ( Spc(I,J,L,id_NH4) &
                                 / AIRVOL(I,J,L) ) / SO4_NH4_NIT(I,J,L)

          ELSE

             ! If SO4_NH4_NIT(I,J,L) is zero, then avoid a div-by-zero
             ! error.  Set all of these to zero because the division
             ! cannot be done.
             FRAC_SNA(I,J,L,1) = 0e+0_fp
             FRAC_SNA(I,J,L,2) = 0e+0_fp
             FRAC_SNA(I,J,L,3) = 0e+0_fp

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
          BCPI(I,J,L) = Spc(I,J,L,id_BCPI) / AIRVOL(I,J,L)

          ! Hydrophobic BC [kg/m3]
          BCPO(I,J,L) = Spc(I,J,L,id_BCPO) / AIRVOL(I,J,L)

          ! Hydrophobic OC [kg/m3]
          ! SOAupdate: Treat either OCPO (x2.1) or POA (x1.4)
          IF ( IS_POA ) THEN
             OCPO(I,J,L) = ( Spc(I,J,L,id_POA1) + Spc(I,J,L,id_POA2) ) &
                           * OCFPOA(I,J) / AIRVOL(I,J,L)
          ELSE IF ( IS_OCPO ) THEN
             OCPO(I,J,L) = Spc(I,J,L,id_OCPO) * OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

          ! Hydrophilic OC [kg/m3]
          IF ( IS_OCPI ) THEN
             OCPI(I,J,L) = Spc(I,J,L,id_OCPI) * OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

          ! Now avoid division by zero (bmy, 4/20/04)
          BCPI(I,J,L)    = MAX( BCPI(I,J,L), 1e-35_fp )
          OCPI(I,J,L)    = MAX( OCPI(I,J,L), 1e-35_fp )
          BCPO(I,J,L)    = MAX( BCPO(I,J,L), 1e-35_fp )
          OCPO(I,J,L)    = MAX( OCPO(I,J,L), 1e-35_fp )

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
          DO K = 1, IBINS

             ! Get the overall species index for species K
             N    = id_DUST1 + K - 1

             ! Effective aerosol radius [m]
             REFF = State_Chm%SpcData(N)%Info%Radius

             ! Bin #1
             IF ( REFF < 0.2e-6_fp ) THEN
                SOILDUST(I,J,L,1) = SOILDUST(I,J,L,1) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #2
             ELSE IF ( REFF < 0.325e-6_fp ) THEN
                SOILDUST(I,J,L,2) = SOILDUST(I,J,L,2) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #3
             ELSE IF ( REFF < 0.6e-6_fp ) THEN
                SOILDUST(I,J,L,3) = SOILDUST(I,J,L,3) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #4
             ELSE IF ( REFF < 1.15e-6_fp ) THEN
                SOILDUST(I,J,L,4) = SOILDUST(I,J,L,4) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #5
             ELSE IF ( REFF < 2.0e-6_fp ) THEN
                SOILDUST(I,J,L,5) = SOILDUST(I,J,L,5) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #6
             ELSE IF ( REFF < 3.25e-6_fp ) THEN
                SOILDUST(I,J,L,6) = SOILDUST(I,J,L,6) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

             ! Bin #7
             ELSE
                SOILDUST(I,J,L,7) = SOILDUST(I,J,L,7) &
                                    + Spc(I,J,L,N) / AIRVOL(I,J,L)

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
          SOILDUST(I,J,L,1) = 0.007e+0_fp  * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,2) = 0.0332e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,3) = 0.2487e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,4) = 0.7111e+0_fp * Spc(I,J,L,id_DST1) / AIRVOL(I,J,L)

          ! Other hetchem bins
          SOILDUST(I,J,L,5) = Spc(I,J,L,id_DST2) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,6) = Spc(I,J,L,id_DST3) / AIRVOL(I,J,L)
          SOILDUST(I,J,L,7) = Spc(I,J,L,id_DST4) / AIRVOL(I,J,L)

       ENDIF

#endif

       !===========================================================
       ! S E A S A L T   A E R O S O L S
       !
       ! Compute accumulation & coarse mode concentration [kg/m3]
       !===========================================================
       IF ( LSSALT ) THEN

          ! Accumulation mode seasalt aerosol [kg/m3]
          SALA(I,J,L) = Spc(I,J,L,id_SALA) / AIRVOL(I,J,L)

          ! Coarse mode seasalt aerosol [kg/m3]
          SALC(I,J,L) = Spc(I,J,L,id_SALC) / AIRVOL(I,J,L)

          ! Avoid division by zero
          SALA(I,J,L) = MAX( SALA(I,J,L), 1e-35_fp )
          SALC(I,J,L) = MAX( SALC(I,J,L), 1e-35_fp )

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
          SOAS(I,J,L) = Spc(I,J,L,id_SOAS) / AIRVOL(I,J,L)

       ENDIF

       !--------------------------------------------------------
       ! Complex SOA scheme
       !--------------------------------------------------------
       IF ( Is_ComplexSOA ) THEN

          ! TSOA (terpene SOA) [kg/m3]
          IF ( IS_TSOA ) THEN
             TSOA(I,J,L) = ( Spc(I,J,L,id_TSOA1) + Spc(I,J,L,id_TSOA2) + &
                             Spc(I,J,L,id_TSOA3) + Spc(I,J,L,id_TSOA0) ) &
                             / AIRVOL(I,J,L)
          ENDIF

          ! ASOA (benz, tolu, xyle, + NAP/IVOC SOA) [kg/m3]
          IF ( IS_ASOA ) THEN
             ASOA(I,J,L) = ( Spc(I,J,L,id_ASOAN) + Spc(I,J,L,id_ASOA1) + &
                             Spc(I,J,L,id_ASOA2) + Spc(I,J,L,id_ASOA3) ) &
                             / AIRVOL(I,J,L)
          ENDIF

          ! OPOA [kg/m3]
          IF ( IS_OPOA ) THEN
             OPOA(I,J,L) = ( Spc(I,J,L,id_OPOA1) + Spc(I,J,L,id_OPOA2) ) &
                            * OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF
       ENDIF

       !-------------------------------------------------------
       ! Mass loading of isoprene SOA (ISOAAQ) [kg/m3]
       !-------------------------------------------------------

       ! Glyoxal
       IF ( id_SOAGX > 0 ) THEN
          ISOAAQ(I,J,L) = Spc(I,J,L,id_SOAGX) / AIRVOL(I,J,L)
       ENDIF

       ! IEPOX
       IF ( id_SOAIE > 0 ) THEN
          ISOAAQ(I,J,L) = ISOAAQ(I,J,L) + Spc(I,J,L,id_SOAIE) / AIRVOL(I,J,L)
       ENDIF

       !-----------------------------------------------------------------------
       ! Exclude INDIOL from AOD and aerosol mass calculations. This results in
       ! lost mass. As noted in Fisher et al. (2016, ACP), this is a source of
       ! uncertainty and would benefit from an update when more information
       ! about this process becomes available. (eam, jaf, mps, 3/5/18)
       !! SOA from alkyl nitrates (some contribution
       !! from non-isoprene sources)
       !IF ( id_INDIOL > 0 ) THEN
       !   ISOAAQ(I,J,L) = ISOAAQ(I,J,L) + Spc(I,J,L,id_INDIOL) / AIRVOL(I,J,L)
       !ENDIF
       !-----------------------------------------------------------------------

       ! SOA from ISOPOOH oxidation product
       IF ( id_LVOCOA > 0 ) THEN
          ISOAAQ(I,J,L) = ISOAAQ(I,J,L) + Spc(I,J,L,id_LVOCOA) / AIRVOL(I,J,L)
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
          OCPISOA(I,J,L) = ( Spc(I,J,L,id_OCPI) * OCFOPOA(I,J) + &
                             Spc(I,J,L,id_SOAS) ) / AIRVOL(I,J,L)

       ELSEIF ( Is_ComplexSOA ) THEN

          OCPISOA(I,J,L) = ( Spc(I,J,L,id_TSOA1) + Spc(I,J,L,id_TSOA2) + &
                             Spc(I,J,L,id_TSOA3) + Spc(I,J,L,id_TSOA0) + &
                             Spc(I,J,L,id_ASOAN) + Spc(I,J,L,id_ASOA1) + &
                             Spc(I,J,L,id_ASOA2) + Spc(I,J,L,id_ASOA3) ) &
                             / AIRVOL(I,J,L)

          IF ( IS_OPOA ) THEN ! hotp 7/28/10
             OCPISOA(I,J,L) = OCPISOA(I,J,L) + &
                              ( Spc(I,J,L,id_OPOA1) + Spc(I,J,L,id_OPOA2) ) &
                              * OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

          IF ( IS_OCPI ) THEN  ! hotp 7/28/10
             OCPISOA(I,J,L) = OCPISOA(I,J,L) + Spc(I,J,L,id_OCPI) &
                              * OCFOPOA(I,J) / AIRVOL(I,J,L)
          ENDIF

       ENDIF

       ! Add mechanistic isoprene OA (eam, 08/2015)
       ! Skip adding this for Simple SOA (jaf, clh, bmy, 5/17/18)
       IF ( Is_ComplexSOA ) THEN
          OCPISOA(I,J,L) = OCPISOA(I,J,L) + ISOAAQ(I,J,L)
       ENDIF

       ! Now avoid division by zero (bmy, 4/20/04)
       OCPISOA(I,J,L) = MAX( OCPISOA(I,J,L), 1e-35_fp )

       !===========================================================
       ! SOAGX [kg/m3]
       !===========================================================
       IF ( IS_SOAGX ) THEN
          SOAGX(I,J,L) = Spc(I,J,L,id_SOAGX) * OCFG / AIRVOL(I,J,L)
       ENDIF

       !==============================================================
       ! P A R T I C U L A T E   M A T T E R
       !
       ! Compute PM2.5 concentration [kg/m3]
       !
       ! PM25 = 1.33 (NH4 + NIT  + SO4) + BCPI + BCPO +
       !        2.10 (OCPO + 1.16 OCPI) + 1.16 SOA*   +
       !        DST1 + 0.38 DST2 + 1.86 SALA
       !
       !   * If using simple  SOA, SOA = SOAS;
       !     If using complex SOA, SOA = TSOA + ASOA + ISOAAQ
       !
       ! NOTES:
       ! - We apply growth factors at 35% RH (computed above):
       !    1.33 for SO4, NIT, and NH4
       !    1.16 for OCPI and SOA
       !    1.86 for SALA
       ! - Ratio of OM/OC = 2.1 is applied to OCPI and OCPO above
       ! - Aerosol WG recommends including 38% of DST2 in PM2.5
       ! - Use either simple SOA or complex SOA in PM2.5 calculation.
       !   By default simple SOA will be used.
       !
       ! %%% IMPORTANT %%%
       ! Note that if complex SOA is used then PM2.5 includes all
       ! the SOA formed in both the Marais et al. and Pye et al.
       ! schemes and may include some double-counting of isoprene SOA.
       ! (Aerosol WG)
       !==============================================================

       ! Units: [kg/m3]
       PM25(I,J,L) = NH4(I,J,L)        * SIA_GROWTH + &
                     NIT(I,J,L)        * SIA_GROWTH + &
                     SO4(I,J,L)        * SIA_GROWTH + &
                     BCPI(I,J,L)                    + &
                     BCPO(I,J,L)                    + &
                     OCPO(I,J,L)                    + &
                     OCPI(I,J,L)       * ORG_GROWTH + &
                     SALA(I,J,L)       * SSA_GROWTH + &
                     SOILDUST(I,J,L,1)              + & ! DST1
                     SOILDUST(I,J,L,2)              + & ! DST1
                     SOILDUST(I,J,L,3)              + & ! DST1
                     SOILDUST(I,J,L,4)              + & ! DST1
                     SOILDUST(I,J,L,5) * 0.38           ! 38% of DST2

       ! Include either simple SOA (default) or Complex SOA in
       ! PM2.5 calculation.  In simulations where both Simple SOA and
       ! Complex SOA species are carried (i.e. "benchmark"), then
       ! only the Simple SOA will be added to PM2.5, in order to avoid
       ! double-counting. (bmy, 5/11/18)
       IF ( Is_SimpleSOA ) THEN
          PM25(I,J,L) = PM25(I,J,L) + SOAS(I,J,L) * ORG_GROWTH
       ELSEIF ( Is_ComplexSOA ) THEN
          PM25(I,J,L) = PM25(I,J,L)                + &
                        TSOA(I,J,L)   * ORG_GROWTH + &
                        ASOA(I,J,L)   * ORG_GROWTH + &
                        ISOAAQ(I,J,L) * ORG_GROWTH    ! Includes SOAGX
       ENDIF

       ! Apply STP correction factor based on ideal gas law
       PM25(I,J,L) = PM25(I,J,L) * ( 1013.25_fp / PMID(I,J,L) ) * &
                     ( T(I,J,L)   / 298.0_fp    )

#ifdef MODEL_GEOS
       ! PM2.5 sulfates
       IF ( State_Diag%Archive_PM25su ) THEN
          State_Diag%PM25su(I,J,L) = ( SO4(I,J,L) * SIA_GROWTH  ) &
                                   * ( 1013.25_fp / PMID(I,J,L) ) &
                                   * ( T(I,J,L)   / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 nitrates
       IF ( State_Diag%Archive_PM25ni ) THEN
          State_Diag%PM25ni(I,J,L) = ( NH4(I,J,L) * SIA_GROWTH    &
                                   +   NIT(I,J,L) * SIA_GROWTH  ) &
                                   * ( 1013.25_fp / PMID(I,J,L) ) &
                                   * ( T(I,J,L)   / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 BC
       IF ( State_Diag%Archive_PM25bc  ) THEN
          State_Diag%PM25bc(I,J,L) = ( BCPI(I,J,L) + BCPO(I,J,L) ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 OC
       IF ( State_Diag%Archive_PM25oc  ) THEN
          State_Diag%PM25oc(I,J,L) = ( OCPO(I,J,L)                 &
                                   +   OCPI(I,J,L) * ORG_GROWTH  ) &
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
          State_Diag%PM25ss(I,J,L) = ( SALA(I,J,L) * ORG_GROWTH  ) &
                                   * ( 1013.25_fp  / PMID(I,J,L) ) &
                                   * ( T(I,J,L)    / 298.0_fp    ) &
                                   * 1.0e+9_fp
       ENDIF

       ! PM2.5 SOA
       IF ( State_Diag%Archive_PM25soa ) THEN
          State_Diag%PM25soa(I,J,L) = ( TSOA(I,J,L)   * ORG_GROWTH    &
                                    +   ASOA(I,J,L)   * ORG_GROWTH    &
                                    +   SOAS(I,J,L)   * ORG_GROWTH    &
                                    +   ISOAAQ(I,J,L) * ORG_GROWTH  ) &
                                    * ( 1013.25_fp    / PMID(I,J,L) ) &
                                    * ( T(I,J,L)      / 298.0_fp    ) &
                                    * 1.0e+9_fp
       ENDIF
#endif

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Convert species back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'End of AEROSOL_CONC in aerosol_mod.F90')
       RETURN
    ENDIF

    ! Free pointers
    Spc    => NULL()
    AIRVOL => NULL()
    OMOC   => NULL()

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
    USE CMN_SIZE_Mod,   ONLY : NRH, NRHAER, NSTRATAER
    USE CMN_FJX_MOD
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
    USE UCX_MOD,        ONLY : NDENS_AER
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
    INTEGER             :: I, J, L, N, R, IRH, W, IRHN
    INTEGER             :: AA, IWV, IIWV, NWVS, IR, NRT
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
    LOGICAL             :: Is_ComplexSOA
    LOGICAL             :: LSTRATOD
    LOGICAL             :: LRAD
    LOGICAL             :: LBCAE  ! (xnw, 8/24/15)
    LOGICAL             :: LUCX
    LOGICAL             :: IS_POA, IS_OCPI
    REAL(fp)            :: GF_RH
    REAL(fp)            :: BCAE_1, BCAE_2

    ! Pointers
    REAL(fp), POINTER   :: BXHEIGHT(:,:,:)
    REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
    REAL(fp), POINTER   :: TAREA(:,:,:,:)
    REAL(fp), POINTER   :: WERADIUS(:,:,:,:)
    REAL(fp), POINTER   :: WTAREA(:,:,:,:)

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
    Is_ComplexSOA        = Input_Opt%LSOA
    LSTRATOD             = Input_Opt%LSTRATOD
    LRAD                 = Input_Opt%LRAD
    LBCAE                = Input_Opt%LBCAE !(xnw, 8/24/15)
    BCAE_1               = Input_Opt%BCAE_1
    BCAE_2               = Input_Opt%BCAE_2
    LUCX                 = Input_Opt%LUCX

    ! Define logical flags
    IS_OCPI              = ( id_OCPI > 0 )
    IS_POA               = ( id_POA1 > 0 .AND. id_POA2 > 0 )

    ! Initialize pointers
    BXHEIGHT            => State_Met%BXHEIGHT    ! Grid box height [m]
    ERADIUS             => State_Chm%AeroRadi    ! Aerosol Radius [cm]
    TAREA               => State_Chm%AeroArea    ! Aerosol Area [cm2/cm3]
    WERADIUS            => State_Chm%WetAeroRadi ! Wet Aerosol Radius [cm]
    WTAREA              => State_Chm%WetAeroArea ! Wet Aerosol Area [cm2/cm3]

    ! Initialize the mapping between hygroscopic species in the
    ! species database and the species order in NRHAER
    IF ( FIRST ) THEN
       DO I = 1, NRHAER

          ! Get the species database index from the species database
          ! mapping array for hygroscopic growth species
          N  = State_Chm%Map_HygGrth(I)

          ! Point to the Species Database entry for species N
          SpcInfo => State_Chm%SpcData(N)%Info

          ! Set the mapping to the local ordering of aerosol densities
          ! in RDAER
          SELECT CASE ( TRIM(SpcInfo%Name) )
          CASE ( 'SO4' )
             Map_NRHAER(I) = 1
          CASE ( 'BCPI' )
             Map_NRHAER(I) = 2
          CASE ( 'OCPI', 'POA1' )
             Map_NRHAER(I) = 3
          CASE ( 'SALA' )
             Map_NRHAER(I) = 4
          CASE ( 'SALC' )
             Map_NRHAER(I) = 5
          CASE DEFAULT
             ErrMsg = 'WARNING: aerosol diagnostics not defined' // &
                      ' for NRHAER greater than 5!'
             CALL GC_ERROR( ErrMsg, RC, 'RDAER in aerosol_mod.F90' )
          END SELECT

       ENDDO
    ENDIF

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
          IF ( Input_Opt%amIRoot ) WRITE( 6, 100 )
100       FORMAT( '     - RDAER: Using online SO4 NH4 NIT!' )
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          WAERSL(I,J,L,1) = SO4_NH4_NIT(I,J,L)
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
          IF ( Input_Opt%amIRoot ) WRITE( 6, 110 )
110       FORMAT( '     - RDAER: Using online BCPI OCPI BCPO OCPO!' )
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Hydrophilic BC (a.k.a EC) [kg/m3]
          WAERSL(I,J,L,2) = BCPI(I,J,L)

          ! Hydrophilic OC [kg/m3]
          WAERSL(I,J,L,3) = OCPISOA(I,J,L)

          ! Hydrophobic BC (a.k.a EC) [kg/m3]
          DAERSL(I,J,L,1) = BCPO(I,J,L)

          ! Hydrophobic OC [kg/m3]
          DAERSL(I,J,L,2) = OCPO(I,J,L)

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
          IF ( Input_Opt%amIRoot ) WRITE( 6, 120 )
120       FORMAT( '     - RDAER: Using online SALA SALC' )
       ENDIF

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Accumulation mode seasalt aerosol [kg/m3]
          WAERSL(I,J,L,4) = SALA(I,J,L)

          ! Coarse mode seasalt aerosol [kg/m3]
          WAERSL(I,J,L,5) = SALC(I,J,L)

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    IF ( LUCX ) THEN
       ! For UCX-based mechanisms transfer stratospheric aerosol data
       ! SDE 04/17/13
       WAERSL(:,:,:,NRHAER+1) = SLA
       WAERSL(:,:,:,NRHAER+2) = SPA
    ENDIF

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

    IF ( LUCX ) THEN
       ! These default values unused (actively retrieved from ucx_mod)
       MSDENS(NRHAER+1) = 1700.0d0 ! SSA/STS
       MSDENS(NRHAER+2) = 1000.0d0 ! NAT/ice PSC
    ENDIF

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
          NWVS = NWVAA-NWVAA0+NWVREQUIRED
       ELSE
          !Loop over wavelengths needed for
          !interpolation to those requested in input.geos
          !(determined in RD_AOD)
          NWVS = NWVREQUIRED
       ENDIF
    ENDIF

    DO IIWV = 1, NWVS
       !now select the correct LUT wavelength
       IF (ODSWITCH .EQ. 0) THEN
          ! only doing for 1000nm (IWV1000 is set in RD_AOD)
          ! N.B. NWVS is fixed to 1 above - only one wavelength
          IWV=IWV1000
       ELSE
          IF ( LRAD ) THEN
             ! RRTMG wavelengths begin after NWVAA0 standard wavelengths
             ! but add on any others required
             IF (IIWV.LE.30) THEN
                !index of RRTMG wavelengths starts after the standard NWVAA0
                !(currently NWVAA0=11, set in CMN_FJX_mod based on the
                ! .dat LUT)
                IWV = IIWV+NWVAA0
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
       ! (this will include strat aerosol when using UCX)
       DO N = 1, NAER

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

          ! Loop over relative humidity bins
          DO R = 1, NRH

             ! Wet radius in aerosol LUT files
             RW(R) = REAA(R,N)

             ! Extinction efficiency for Q for each RH bin
             QW(R)   = QQAA(IWV,R,N)
             AW(R)   = ALPHAA(IWV,R,N)
             SSW(R)  = SSAA(IWV,R,N)
             ASYW(R) = ASYMAA(IWV,R,N)

          ENDDO

          ! Loop over grid boxes
          !$OMP PARALLEL DO                                                 &
          !$OMP PRIVATE( I,        J,       L,        IRH                 ) &
          !$OMP PRIVATE( AW0,      QW0,     SSW0,     ASYW0,    REFF      ) &
          !$OMP PRIVATE( SCALEA,   SCALEQ,  SCALESSA, SCALEASY, FRAC      ) &
          !$OMP PRIVATE( SCALER,   SCALEOD, SCALEVOL, DRYAREA,  TAERVOL   ) &
          !$OMP PRIVATE( TK,       CONSEXP, VPRESH2O, RELHUM              ) &
#ifdef RRTMG
          !$OMP PRIVATE( IR                                               ) &
#endif
          !$OMP PRIVATE( RHOSTRAT, RAER,    SADSTRAT, XSASTRAT            ) &
          !$OMP PRIVATE( VDRY,     VH2O                                   ) &
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

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

             IF ( .not. LUCX .or. ( LUCX .and. (N.LE.NRHAER)) ) THEN

                !--------------------------------------------------------
                ! %%%%%%% Aerosols undergoing hygroscopic growth %%%%%%%
                !--------------------------------------------------------

                !calculate optics for hyrdophillic aerosol here
                !However MDENS in LUT was in g/cm3 not kg/m3 so x1e3
                ODAER(I,J,L,IWV,N) = SCALEOD * BXHEIGHT(I,J,L) * 0.75d0 * &
                                     WAERSL(I,J,L,N) * QQAA(IWV,1,N)    / &
                                     ( MSDENS(N) * REAA(1,N) * 1.0D-6 )

                !Include BC absorption enhancement (xnw, 8/24/15)
                IF (N.eq.2) THEN

                   IF (LBCAE) THEN
                      BCSCAT_AE = ODAER(I,J,L,IWV,N)*SCALESSA*SSAA(IWV,1,N)
                      ODAER(I,J,L,IWV,N) = ODAER(I,J,L,IWV,N) * &
                                ( BCAE_1 + SCALESSA*SSAA(IWV,1,N) - &
                                  SCALESSA*SSAA(IWV,1,N)*BCAE_1 )

                      !now combine with hydrophilic OD as before
                      BCSCAT_AE = BCSCAT_AE + SSAA(IWV,1,N) * &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  DAERSL(I,J,L,N-1) * QQAA(IWV,1,N)  / &
                                  ( MSDENS(N) * REAA(1,N) * 1.0D-6 )
                      ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                           (BCAE_2+SSAA(IWV,1,N) - SSAA(IWV,1,N)*BCAE_2) * &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  DAERSL(I,J,L,N-1) * QQAA(IWV,1,N)  / &
                                  ( MSDENS(N) * REAA(1,N) * 1.0D-6 )

                   ELSE
                      !now combine with hydrophilic OD as before
                      ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                                  0.75d0 * BXHEIGHT(I,J,L) * &
                                  DAERSL(I,J,L,N-1) * QQAA(IWV,1,N)  / &
                                  ( MSDENS(N) * REAA(1,N) * 1.0D-6 )
                   ENDIF

                ENDIF

                IF (N.eq.3) THEN
                   !now combine with hydrophilic OD as before
                   ODAER(I,J,L,IWV,N)= ODAER(I,J,L,IWV,N) + &
                                   0.75d0 * BXHEIGHT(I,J,L) * &
                                   DAERSL(I,J,L,N-1) * QQAA(IWV,1,N)  / &
                                   ( MSDENS(N) * REAA(1,N) * 1.0D-6 )
                ENDIF

                ! Get the AOD contribution from isoprene SOA only (eam, 2014)
                IF ( N == 3 .and. Is_ComplexSOA ) THEN
                   ISOPOD(I,J,L,IWV) = SCALEOD*BXHEIGHT(I,J,L)*0.75d0 &
                                   * ISOAAQ(I,J,L) * QQAA(IWV,1,N) / &
                                   ( MSDENS(N) * REAA(1,N) * 1.0D-6 )
                ENDIF

             ELSE

                !--------------------------------------------------------
                !        %%%%%%% Stratospheric aerosols %%%%%%%
                !--------------------------------------------------------

                ! Get aerosol effective radius
                CALL GET_STRAT_OPT(I, J, L, ISTRAT, RAER, REFF, &
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
                ODAER(I,J,L,IWV,N) = BXHEIGHT(I,J,L) * XSASTRAT * QQAA(IWV,1,N)

             ENDIF

#ifdef RRTMG
             !SNA currently treated as one with optics but considered
             !separately for RT, so we split them by mass here
             IF (N.EQ.1) THEN
                DO IR=1,3
                   RTODAER(I,J,L,IWV,N+IR-1)= ODAER(I,J,L,IWV,N)* &
                                              FRAC_SNA(I,J,L,IR)
                   RTSSAER(I,J,L,IWV,N+IR-1)   = SCALESSA*SSAA(IWV,1,N)
                   RTASYMAER(I,J,L,IWV,N+IR-1) = SCALEASY*ASYMAA(IWV,1,N)
                ENDDO
             ELSE
                !RT arrays now offset from NAER by 2 (NRT=N+2 for N>1)
                !If strat aerosol switched on (UCX) then this will
                !automatically be added after the standard aerosol
                !(NRHAER+1,2) but before dust
                RTODAER(I,J,L,IWV,NRT)     = ODAER(I,J,L,IWV,N)
                RTSSAER(I,J,L,IWV,NRT)     = SCALESSA*SSAA(IWV,1,N)
                !for BC SSA with absorption enhancement (xnw 8/24/15)
                IF ((N .EQ. 2) .AND. (LBCAE)) THEN
                   RTSSAER(I,J,L,IWV,NRT)  = BCSCAT_AE / &
                                             ODAER(I,J,L,IWV,N)
                ENDIF
                RTASYMAER(I,J,L,IWV,NRT)   = SCALEASY*ASYMAA(IWV,1,N)
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
                   State_Diag%AerHygGrowth(I,J,L,Map_NRHAER(N)) = SCALEOD
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

                ! For SO4-NIT-NH4 aerosol, re-calculate the wet effective
                ! radius using the water content from ISORROPIA.
                ! This new effective radius will be used for surface area
                ! used in heterogeneous chemistry. We don't use this
                ! effective radius in the optics above (OD, scattering,
                ! absorption) because the index of refraction, phase
                ! function, and Q must all be consistent with the radius and
                ! composition.
                ! Note: ISORROPIA water includes fine sea salt aerosol,
                ! which we are assigning all to SNA here without decreasing
                ! the sea salt volume. This double counts the fine SSA
                ! volume. (cdholmes, 5/17/2019)
                IF (N == 1) THEN

                   ! Volume of water, m3(H2O)/m3(air)
                   ! AeroH2O has units g/m3
                   VH2O = State_Chm%AeroH2O(I,J,L,NDUST+1) / 1e6
                   ! Volume of dry aerosol, m3(aerosol)/m3(air)
                   VDry = WAERSL(I,J,L,N) / MSDENS(N)

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
                   REFF = RW(1) * &
                          ( 1d0 + safe_div( VH2O, VDry, 0d0 ) )**(1d0/3d0)

                ENDIF

                !get scaling for R and VOL
                SCALER                 = REFF / RW(1)
                SCALEVOL               = SCALER**3
                ERADIUS(I,J,L,N+NDUST) = 1.0D-4 * REFF

                ! Store aerosol surface areas in TAREA, and be sure
                ! to list them following the dust surface areas
                TAREA(I,J,L,N+NDUST)   = 3.D0 * WAERSL(I,J,L,N) * SCALEVOL / &
                                         ( ERADIUS(I,J,L,N+NDUST) * MSDENS(N) )

                WTAREA(I,J,L,N+NDUST)   = TAREA(I,J,L,N+NDUST)
                WERADIUS(I,J,L,N+NDUST) = ERADIUS(I,J,L,N+NDUST)

                ! Save aerosol water content. Assume that the increase in volume
                ! equals the volume of pure water added, m3(H2O)/m3(air),
                ! then convert to g/m3
                State_Chm%AeroH2O(I,J,L,N+NDUST) = 1e+6_fp * &
                   WAERSL(I,J,L,N) / MSDENS(N) * (ScaleVol - 1d0)

                !include hydrophobic BC and OC
                !stored separate to hydrophillic in RT variables
                IF ((N.eq.2).or.(N.eq.3)) THEN

                   ! Dry surface area
                   ! SDE 2015-10-27: RW is in um, but everything
                   ! else is in terms of cm. Correct with 10^-4 factor
                   DRYAREA = 3.D0 * DAERSL(I,J,L,N-1) / ( RW(1) * &
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
    ! Account for stratospheric aerosols in UCX mechanisms (SDE 04/17/13)
    !==============================================================
    IF ( LUCX ) THEN

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
             CALL GET_STRAT_OPT(I,J,L,ISTRAT,RAER,REFF,SADSTRAT,XSASTRAT)

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
                           NDENS_AER(I,J,L,ISTRAT)*1.d-6
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
                           NDENS_AER(I,J,L,ISTRAT)*1.d-6
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

    ENDIF    ! end if UCX

    !==============================================================
    ! Non-stratospheric cloud diagnostics
    ! NOTE: The cloud optical depths are actually recorded at
    !       1000 nm, but vary little with wavelength.
    !==============================================================
    IF ( State_Diag%Archive_AOD .and. ODSWITCH .EQ. 1 ) THEN

       ! Loop over aerosol types (dust handled in dust_mod.F90)
       !$OMP PARALLEL DO                                            &
       !$OMP DEFAULT( SHARED )                                      &
       !$OMP PRIVATE( I, J, L, N, W, LINTERP, IsWL1, IsWL2, IsWL3 ) &
       !$OMP SCHEDULE( DYNAMIC )
       ! Loop over hydroscopic aerosols
       DO N = 1, NRHAER

          ! Loop over wavelengths set in input.geos radiation menu
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
                      State_Diag%AODHygWL1(I,J,L,Map_NRHAER(N)) = &
                           ODAER(I,J,L,IWVSELECT(1,W),N)
                   ELSEIF ( State_Diag%Archive_AODHygWL2 .AND. IsWL2 ) THEN
                      State_Diag%AODHygWL2(I,J,L,Map_NRHAER(N)) = &
                           ODAER(I,J,L,IWVSELECT(1,W),N)
                   ELSEIF ( State_Diag%Archive_AODHygWL3 .AND. IsWL3 ) THEN
                      State_Diag%AODHygWL3(I,J,L,Map_NRHAER(N)) = &
                           ODAER(I,J,L,IWVSELECT(1,W),N)
                   ENDIF
                ELSE
                   ! Interpolated using angstrom exponent between
                   ! Closest available wavelengths
                   ! (coefs pre-calculated in CALC_AOD (RD_AOD)
                   !catch any zero values before interpolation
                   IF ((ODAER(I,J,L,IWVSELECT(2,W),N).GT.0).AND. &
                       (ODAER(I,J,L,IWVSELECT(1,W),N).GT.0)) THEN
                      IF ( State_Diag%Archive_AODHygWL1 .AND. IsWL1 ) THEN
                         State_Diag%AODHygWL1(I,J,L,Map_NRHAER(N)) =          &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**     &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/ &
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
                      ELSEIF ( State_Diag%Archive_AODHygWL2 .AND. IsWL2 ) THEN
                         State_Diag%AODHygWL2(I,J,L,Map_NRHAER(N)) =          &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**     &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/ &
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
                      ELSEIF ( State_Diag%Archive_AODHygWL3 .AND. IsWL3 ) THEN
                         State_Diag%AODHygWL3(I,J,L,Map_NRHAER(N)) =          &
                              ODAER(I,J,L,IWVSELECT(2,W),N)*ACOEF_WV(W)**     &
                              (BCOEF_WV(W)*LOG(ODAER(I,J,L,IWVSELECT(1,W),N)/ &
                              ODAER(I,J,L,IWVSELECT(2,W),N)))
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

       !$OMP PARALLEL DO           &
       !$OMP DEFAULT( SHARED )     &
       !$OMP PRIVATE( I, J, L, N ) &
       !$OMP SCHEDULE( DYNAMIC )
       ! Loop over hydroscopic aerosols
       DO N = 1, NRHAER

          !----------------------------------------------------
          ! Netcdf diagnostics computed here:
          !  Sulfate Surface Area                       [cm2/cm3]
          !  Black Carbon (hydrophilic) Surface Area    [cm2/cm3]
          !  Organic Carbon (hydrophilic) Surface Area  [cm2/cm3]
          !  Sea Salt (accum) Surface Area              [cm2/cm3]
          !  Sea Salt (coarse) Surface Area             [cm2/cm3]
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             State_Diag%AerSurfAreaHyg(I,J,L,Map_NRHAER(N)) = &
                  TAREA(I,J,L,N+NDUST)
          ENDDO
          ENDDO
          ENDDO

       ENDDO ! end of loop over hygroscopic aerosols
       !$OMP END PARALLEL DO

    ENDIF

    ! Turn off radiative effects of stratospheric aerosols?
    IF ( LUCX .and. .not.(LSTRATOD)) THEN
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
    NULLIFY( BXHEIGHT, ERADIUS, TAREA, WERADIUS, WTAREA )

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
! !DESCRIPTION: Subroutine INIT\_AEROSOL allocates and zeroes module arrays
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Aerosol( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,   ONLY : NAER, NDUST
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
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
    INTEGER            :: AS
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_AEROSOL begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_Aerosol (in module GeosCore/aerosol_mod.F90)'

    ! Add tracer ID flags as module variables (bmy, 6/16/16)
    id_BCPI   = Ind_( 'BCPI'   )
    id_BCPO   = Ind_( 'BCPO'   )
    id_DST1   = Ind_( 'DST1'   )
    id_DST2   = Ind_( 'DST2'   )
    id_DST3   = Ind_( 'DST3'   )
    id_DST4   = Ind_( 'DST4'   )
    id_DUST1  = Ind_( 'DUST1'  )
    id_NH4    = Ind_( 'NH4'    )
    id_NIT    = Ind_( 'NIT'    )
    id_OCPO   = Ind_( 'OCPO'   )
    id_OCPI   = Ind_( 'OCPI'   )
    id_SOAS   = Ind_( 'SOAS'   )
    id_SALA   = Ind_( 'SALA'   )
    id_SALC   = Ind_( 'SALC'   )
    id_SO4    = Ind_( 'SO4'    )
    id_SO4s   = Ind_( 'SO4s'   )
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

    !=================================================================
    ! Allocate arrays
    !=================================================================
    ALLOCATE( BCPI( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BCPI = 0.0_fp

    ALLOCATE( BCPO( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:BCPO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    BCPO = 0.0_fp

    ALLOCATE( OCPI( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OCPI', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OCPI = 0.0_fp

    ALLOCATE( OCPO( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OCPO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OCPO = 0.0_fp

    ALLOCATE( OCPISOA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OCPISOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OCPISOA = 0.0_fp

    ALLOCATE( SALA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SALA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SALA = 0.0_fp

    ALLOCATE( SALC( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SALC', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SALC = 0.0_fp

    ALLOCATE( SO4_NH4_NIT(State_Grid%NX,State_Grid%NY,State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SO4_NH4_NIT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SO4_NH4_NIT = 0.0_fp

    ALLOCATE( SO4( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SO4', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SO4 = 0.0_fp

    ALLOCATE( NH4( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F:NH4', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NH4 = 0.0_fp

    ALLOCATE( NIT( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:NIT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    NIT = 0.0_fp

    ALLOCATE( FRAC_SNA( State_Grid%NX, State_Grid%NY, State_Grid%NZ, 3 ), &
              STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:FRAC_SNA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    FRAC_SNA = 0.0_fp

    ALLOCATE( SOILDUST( State_Grid%NX, State_Grid%NY, State_Grid%NZ, NDUST ), &
              STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SOILDUST', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SOILDUST = 0.0_fp

    IF ( Input_Opt%LUCX ) THEN
       ALLOCATE( SLA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
       CALL GC_CheckVar( 'aerosol_mod.F90:SLA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SLA = 0.0_fp

       ALLOCATE( SPA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
       CALL GC_CheckVar( 'aerosol_mod.F90:SPA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SPA   = 0.0_fp
    ENDIF

    ALLOCATE( TSOA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:TSOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TSOA = 0.0_fp

    ALLOCATE( ASOA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:ASOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ASOA = 0.0_fp

    ALLOCATE( OPOA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OPOA = 0.0_fp

    ALLOCATE( PM25( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:PM25', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PM25 = 0.0_fp

    ALLOCATE( SOAGX( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SOAGX', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SOAGX = 0.0_fp

    ! Mass of hydrophobic aerosol from Mian Chin
    ALLOCATE( DAERSL(State_Grid%NX,State_Grid%NY,State_Grid%NZ,2), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:DAERSL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    DAERSL = 0.0_fp

    ! Mass of hydrophilic aerosol from Mian Chin
    ALLOCATE( WAERSL(State_Grid%NX,State_Grid%NY,State_Grid%NZ,NAER), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:WAERSL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    WAERSL = 0.0_fp

    ! Mechanistic isoprene SOA (eam, 2014):
    ALLOCATE( ISOAAQ(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:ISOAAQ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ISOAAQ = 0.0_fp

    ! Simple SOA
    ALLOCATE( SOAS(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:SOAS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SOAS = 0.0_fp

    ! OM/OC for POA
    ALLOCATE( OCFPOA(State_Grid%NX,State_Grid%NY), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OCFPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OCFPOA = 0.0_fp

    ! OM/OC for OPOA, OCPI, OCPO
    ALLOCATE( OCFOPOA(State_Grid%NX,State_Grid%NY), STAT=RC )
    CALL GC_CheckVar( 'aerosol_mod.F90:OCFOPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    OCFOPOA = 0.0_fp

  END SUBROUTINE INIT_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_aerosol
!
! !DESCRIPTION: Subroutine CLEANUP\_AEROSOL deallocates all module arrays
!  (bmy, 7/20/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_AEROSOL
!
! !REVISION HISTORY:
!  20 Jul 2004 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    IF ( ALLOCATED( BCPI        ) ) DEALLOCATE( BCPI        )
    IF ( ALLOCATED( BCPO        ) ) DEALLOCATE( BCPO        )
    IF ( ALLOCATED( OCPI        ) ) DEALLOCATE( OCPI        )
    IF ( ALLOCATED( OCPO        ) ) DEALLOCATE( OCPO        )
    IF ( ALLOCATED( OCPISOA     ) ) DEALLOCATE( OCPISOA     )
    IF ( ALLOCATED( SALA        ) ) DEALLOCATE( SALA        )
    IF ( ALLOCATED( SALC        ) ) DEALLOCATE( SALC        )
    IF ( ALLOCATED( SO4_NH4_NIT ) ) DEALLOCATE( SO4_NH4_NIT )
    IF ( ALLOCATED( SO4         ) ) DEALLOCATE( SO4         )
    IF ( ALLOCATED( NH4         ) ) DEALLOCATE( NH4         )
    IF ( ALLOCATED( NIT         ) ) DEALLOCATE( NIT         )
    IF ( ALLOCATED( FRAC_SNA    ) ) DEALLOCATE( FRAC_SNA    )
    IF ( ALLOCATED( SOILDUST    ) ) DEALLOCATE( SOILDUST    )
    IF ( ALLOCATED( SLA         ) ) DEALLOCATE( SLA         )
    IF ( ALLOCATED( SPA         ) ) DEALLOCATE( SPA         )
    IF ( ALLOCATED( TSOA        ) ) DEALLOCATE( TSOA        )
    IF ( ALLOCATED( ASOA        ) ) DEALLOCATE( ASOA        )
    IF ( ALLOCATED( OPOA        ) ) DEALLOCATE( OPOA        )
    IF ( ALLOCATED( SOAGX       ) ) DEALLOCATE( SOAGX       )
    IF ( ALLOCATED( PM25        ) ) DEALLOCATE( PM25        )
    IF ( ALLOCATED( WAERSL      ) ) DEALLOCATE( WAERSL      )
    IF ( ALLOCATED( DAERSL      ) ) DEALLOCATE( DAERSL      )
    IF ( ALLOCATED( ISOAAQ      ) ) DEALLOCATE( ISOAAQ      )
    IF ( ALLOCATED( SOAS        ) ) DEALLOCATE( SOAS        )
    IF ( ALLOCATED( OCFPOA      ) ) DEALLOCATE( OCFPOA      )
    IF ( ALLOCATED( OCFOPOA     ) ) DEALLOCATE( OCFOPOA     )

  END SUBROUTINE CLEANUP_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aermass_diagnostic
!
! !DESCRIPTION: Computes the aerosol mass diagnostic (formerly ND42 bpch
!  diagnostic).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE PhysConstants,  ONLY : MwCarb
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState),   INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  NOTE: This diagnostic mimics the bpch diagnostic routine "DIAG42".
!
! !REVISION HISTORY:
!  05 Feb 2018 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL                  :: First = .TRUE.

    ! Scalars
    INTEGER                  :: I, J, L

    ! Strings
    CHARACTER(LEN=63)        :: OrigUnit
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg

    ! Pointers
    REAL(fp),      POINTER   :: AirDen(:,:,:  )
    REAL(fp),      POINTER   :: Spc   (:,:,:,:)
    TYPE(Species), POINTER   :: SpcInfo
!
! !DEFINED PARAMETERS:
!
    ! Convert [kg/m3] to [ug/m3]
    REAL(fp),      PARAMETER :: kgm3_to_ugm3 = 1.0e+9_fp

    ! Define number of carbon atoms in each irreversible isoprene
    ! SOA tracer species. Named according to the parent HC (same
    ! number of carbons):
    REAL(fp),      PARAMETER :: NCIMAE   = 4e+0_fp
    REAL(fp),      PARAMETER :: NCIEPOX  = 5e+0_fp
    REAL(fp),      PARAMETER :: NCINDIOL = NCIEPOX
    REAL(fp),      PARAMETER :: NCGLYX   = 2e+0_fp
    REAL(fp),      PARAMETER :: NCGLYC   = NCGLYX
    REAL(fp),      PARAMETER :: NCMGLY   = 3e+0_fp
    REAL(fp),      PARAMETER :: NCLVOC   = NCIEPOX
    REAL(fp),      PARAMETER :: NCISN1   = NCIEPOX

    !=======================================================================
    ! Set_AerMass_Diagnostic begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Set_AerMass_Diagnostic (in module GeosCore/aerosol_mod.F90)'

    ! Check that species units are kg/kg dry air
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       CALL GC_Error( 'Incorrect species units: ' // &
                      State_Chm%Spc_Units, RC, ThisLoc )
       RETURN
    ENDIF

    ! Define species ID flags for the aerosol mass diagnostics
    IF ( First ) THEN

       !--------------------------------------------------------------------
       ! Look up species indices in State_Chm%SPECIES
       !--------------------------------------------------------------------
       id_INDIOL  = Ind_( 'INDIOL' )
       id_LVOCOA  = Ind_( 'LVOCOA' )
       id_POA1    = Ind_( 'POA1'   )
       id_POA2    = Ind_( 'POA2'   )
       id_SOAGX   = Ind_( 'SOAGX'  )
       id_SOAIE   = Ind_( 'SOAIE'  )
       Is_POA     = ( id_POA1 > 0 .and. id_POA2 > 0 )

       ! Initialize conversion factors for total OC diagnostic
       Fac_INDIOL = 0.0_fp
       Fac_LVOCOA = 0.0_fp
       Fac_SOAGX  = 0.0_fp
       Fac_SOAIE  = 0.0_fp

       !--------------------------------------------------------------------
       ! Set conversion factors for certain isoprene SOA species,
       ! or, if they aren't present, disable their diagnostics
       !--------------------------------------------------------------------
       IF ( id_INDIOL > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_INDIOL)%Info
          Fac_INDIOL =  ( NCINDIOL  * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( State_Diag%Archive_AerMassINDIOL ) THEN
             State_Diag%Archive_AerMassINDIOL = .FALSE.
             ErrMsg = 'Disabling AerMassINDIOL diagnostic. ' // &
                      'INDIOL is not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_LVOCOA > 0  ) THEN
          SpcInfo    => State_Chm%SpcData(id_LVOCOA)%Info
          Fac_LVOCOA = ( NCLVOC * MwCarb / ( SpcInfo%Mw_G * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( State_Diag%Archive_AerMassLVOCOA ) THEN
             State_Diag%Archive_AerMassLVOCOA = .FALSE.
             ErrMsg = 'Disabling AerMassLVOCOA diagnostic. ' // &
                      'LVOCOA is not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAGX > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAGX)%Info
          Fac_SOAGX  = ( NCGLYX * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( State_Diag%Archive_AerMassSOAGX ) THEN
             State_Diag%Archive_AerMassSOAGX = .FALSE.
             ErrMsg = 'Disabling AerMassSOAGX diagnostic.' // &
                      'SOAGX is not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAIE > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAIE)%Info
          Fac_SOAIE  =  ( NCIEPOX * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( State_Diag%Archive_AerMassSOAIE ) THEN
             State_Diag%Archive_AerMassSOAIE = .FALSE.
             ErrMsg = 'Disabling AerMassSOAIE diagnostic. ' // &
                      'SOAIE is not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       ! Reset first-time flag
       First = .FALSE.
    ENDIF

    !=======================================================================
    ! Compute Aerosol mass and PM2.5 diagnostics using concentrations
    ! from the end of the chemistry timestep, which should be more
    ! consistent with the legacy ND42 bpch diagnostics
    !=======================================================================

    ! Point to fielss of State_Chm and State_Met
    Spc    => State_Chm%Species
    AirDen => State_Met%AIRDEN

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED   ) &
    !$OMP PRIVATE( I, J, L  )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !--------------------------------------
       ! AerMassASOA [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassASOA ) THEN
          State_Diag%AerMassASOA(I,J,L) = ASOA(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassBC [ug C/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassBC ) THEN
          State_Diag%AerMassBC(I,J,L) = ( BCPI(I,J,L) + BCPO(I,J,L) ) * &
                                          kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassINDIOL [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassINDIOL ) THEN
          State_Diag%AerMassINDIOL(I,J,L) = Spc(I,J,L,id_INDIOL) * &
                                            kgm3_to_ugm3 * AirDen(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassLVOCOA [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassLVOCOA ) THEN
          State_Diag%AerMassLVOCOA(I,J,L) = Spc(I,J,L,id_LVOCOA) * &
                                            kgm3_to_ugm3 * AirDen(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassNH4 [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassNH4 ) THEN
          State_Diag%AerMassNH4(I,J,L) = NH4(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassNIT [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassNIT ) THEN
          State_Diag%AerMassNIT(I,J,L) = NIT(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassOPOA [ug/m3], OA:OC=2.1
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassOPOA ) THEN
          State_Diag%AerMassOPOA(I,J,L) = OPOA(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassPOA [ug/m3], OA:OC=2.1
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassPOA ) THEN
          IF ( Is_POA ) THEN
             State_Diag%AerMassPOA(I,J,L) = OCPO(I,J,L) * kgm3_to_ugm3
          ELSE
             State_Diag%AerMassPOA(I,J,L) = ( OCPI(I,J,L) + OCPO(I,J,L) ) * &
                                              kgm3_to_ugm3
          ENDIF
       ENDIF

       !--------------------------------------
       ! AerMassSAL [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassSAL ) THEN
          State_Diag%AerMassSAL(I,J,L) = ( SALA(I,J,L) + SALC(I,J,L) ) * &
                                           kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSO4 [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassSO4 ) THEN
          State_Diag%AerMassSO4(I,J,L) = SO4(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSOAGX [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassSOAGX ) THEN
          State_Diag%AerMassSOAGX(I,J,L) = SOAGX(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSOAIE [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassSOAIE ) THEN
          State_Diag%AerMassSOAIE(I,J,L) = Spc(I,J,L,id_SOAIE) * &
                                           kgm3_to_ugm3 * AirDen(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassTSOA [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_AerMassTSOA ) THEN
          State_Diag%AerMassTSOA(I,J,L) = TSOA(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! PM25 [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_PM25 ) THEN
          State_Diag%PM25(I,J,L) = PM25(I,J,L) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! Sum of all biogenic organic aerosol
       !--------------------------------------
       IF ( State_Diag%Archive_TotalBiogenicOA ) THEN
          State_Diag%TotalBiogenicOA(I,J,L) = ( TSOA(I,J,L) + ISOAAQ(I,J,L) ) &
                                                * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! Sum of all organic aerosol [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_TotalOA ) THEN
          State_Diag%TotalOA(I,J,L) = ( TSOA(I,J,L) + &
                                        ASOA(I,J,L) + &
                                        OCPO(I,J,L) + &
                                        OCPI(I,J,L) + &
                                        OPOA(I,J,L) + &
                                        ISOAAQ(I,J,L) ) * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! Sum of all organic carbon [ug/m3]
       !--------------------------------------
       IF ( State_Diag%Archive_TotalOC ) THEN

          IF ( Is_POA ) THEN
             State_Diag%TotalOC(I,J,L) = &
                  ( ( TSOA(I,J,L) + ASOA(I,J,L) &
                    + OCPI(I,J,L) + OPOA(I,J,L) ) / OCFOPOA(I,J) &
                    + OCPO(I,J,L) / OCFPOA(I,J) ) * kgm3_to_ugm3

          ELSE
             State_Diag%TotalOC(I,J,L) = &
                  ( ( TSOA(I,J,L) + ASOA(I,J,L) &
                    + OCPO(I,J,L) + OCPI(I,J,L) + OPOA(I,J,L) ) &
                    / OCFOPOA(I,J) ) * kgm3_to_ugm3
          ENDIF

          IF ( Input_Opt%LSOA ) THEN
             State_Diag%TotalOC(I,J,L) =  State_Diag%TotalOC(I,J,L) + &
                  ( ( Spc(I,J,L,id_SOAIE)  * Fac_SOAIE  ) + &
                    ( Spc(I,J,L,id_INDIOL) * Fac_INDIOL ) + &
                    ( Spc(I,J,L,id_SOAGX)  * Fac_SOAGX  ) + &
                    ( Spc(I,J,L,id_LVOCOA) * Fac_LVOCOA ) ) &
                    * AirDen(I,J,L) * kgm3_to_ugm3
          ENDIF

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    Spc      => NULL()
    AirDen   => NULL()

  END SUBROUTINE Set_AerMass_Diagnostic
!EOC
END MODULE AEROSOL_MOD
