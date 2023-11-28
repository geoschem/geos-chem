!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cldj_interface_mod.F90
!
! !DESCRIPTION: Module CLDJ\_INTERFACE\_MOD contains routines and variables
!  for interfacing with the Cloud-J scheme (Prather et al) that calculates
!  photolysis rates.
!\\
!\\
! !INTERFACE:
!
MODULE CLDJ_INTERFACE_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_CloudJ
  PUBLIC  :: Run_CloudJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Set_Clim_Profiles
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version, adapted from fast_jx_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: int_cloudj
!
! !DESCRIPTION: Subroutine INIT\_CLOUDJ initializes Cloud-J variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_CLOUDJ( Input_Opt, State_Grid, State_Diag, State_Chm, RC )
!
! !USES:
!

! ewl: Use, inputs/outputs, and local vars could be slimmed down
    ! ewl: if these are in cloud-j, why do I need to pass them???
    USE Cldj_Cmn_Mod,   ONLY : JVN_, NRatJ, W_
    USE Cldj_Init_Mod,  ONLY : Init_CldJ
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)     :: State_Grid  ! Grid State object
    TYPE(DgnState), INTENT(IN)     :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry State object

!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version, adapted from fast_jx_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: NJXX ! ewl: what is this?
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_) ! ewl: do we do anything with this??

    !=================================================================
    ! INIT_CLOUDJ begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_CloudJ (in module GeosCore/cldj_interface_mod.F90)'

    ! Skip these operations when running in dry-run mode
    IF ( .NOT. Input_Opt%DryRun ) THEN

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Cloud-J'

          ! ewl: can this be put into the initialization???
          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
             ErrMsg =  ' INIT_CLOUDJ: invalid no. wavelengths'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          endif
       ENDIF

    ENDIF

    ! Initialize Cloud-J. Includes reading input data files
    ! FJX_spec.dat (RD_XXX), FJX_scat-aer.dat (RD_MIE), and 
    ! FJX_j2j.dat (RD_JS_JX)
    CALL Init_CldJ(Input_Opt%CloudJ_Dir, State_Grid%NZ, TITLEJXX, JVN_, NJXX)

    ! Store # of photolysis reactions in State_Chm object
    State_Chm%Phot%nPhotRxns = NRatJ

  END SUBROUTINE INIT_CLOUDJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_cloudj
!
! !DESCRIPTION: Subroutine RUN\_CLOUDJ loops over horizontal grid boxes to call
!  Cloud-J subroutine CLOUD\_JX for computation of J-Values for each column.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_CloudJ( Input_Opt, State_Chm, State_Diag, &
                         State_Grid, State_Met, RC )
!
! !USES:
!
    USE Aerosol_Mod,    ONLY : SO4_NH4_NIT, BCPI, BCPO, OCPO, OCPISOA
    USE Aerosol_Mod,    ONLY : SALA, SALC, SLA, SPA
    ! ewl: if these are in cloud-j, why do we need them here??
    USE Cldj_Cmn_Mod,   ONLY : L_, L1_, W_, S_, LWEPAR
    USE Cldj_Cmn_Mod,   ONLY : JVN_, AN_, NQD_, W_r
    USE Cldj_Cmn_Mod,   ONLY : JIND, JFACTA, FL, QAA, RAA, SAA
    USE Cld_Sub_Mod,    ONLY : Cloud_JX
    USE Cmn_Size_Mod,   ONLY : NRHAER, NRH, NDUST
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AVO, H2OMW, AIRMW, G0_100, PI, PI_180
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
!ewl    USE TOMS_MOD,       ONLY : GET_OVERHEAD_O3

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm     ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid    ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    INTEGER            :: A, I, J, L, K, N, S, MaxLev, RH_ind
    INTEGER            :: SO4_ind, BC_ind, OC_ind, SALA_ind, SALC_ind
    INTEGER            :: S_rh0, S_rhx, K_rh0, K_rhx, ind_1000
    REAL(8)            :: MW_g, BoxHt, Delta_P, IWC, LWC
    REAL(8)            :: FRAC, RAA_eff, QAA_eff, SAA_eff
    REAL(8)            :: dry_to_wet_factor
    REAL(8)            :: R_interp_factor, Q_interp_factor
    REAL(fp)           :: RH_lut(NRH)
    LOGICAL, SAVE      :: FIRST = .true.

    !------------------------------------------------------------------------
    ! Solar_JX inputs
    !------------------------------------------------------------------------
    INTEGER  :: DAY_OF_YEAR               ! simulation day of year
    REAL(fp) :: U0                        ! cosine of SZA

    !------------------------------------------------------------------------
    ! Set_Clim_Profiles inputs
    !------------------------------------------------------------------------
    INTEGER  :: MONTH, DAY                ! simulation month and day
    REAL(fp) :: T_CTM   (State_Grid%NZ+1) ! temperature profile [K]
    REAL(fp) :: P_CTM   (State_Grid%NZ+2) ! pressure profile (edges) [hPa]
    REAL(fp) :: O3_CTM  (State_Grid%NZ+1) ! ozone profile [molec/cm3]

    !------------------------------------------------------------------------
    ! Cloud_JX inputs
    !------------------------------------------------------------------------

    ! Scalars
    LOGICAL  :: LPRTJ   ! Debug prints
    INTEGER  :: IRAN
    REAL(fp) :: SZA             ! Computed in Solar_JX. Should this be real8?
    REAL(fp) :: SOLF            ! Computed in Solar_JX. Should this be real8?
    REAL(8)  :: CLDCOR

    ! 1D arrays
    INTEGER  :: CLDIW   (L1_  )
    REAL(fp) :: T_CLIM  (L1_  ) ! Computed in Set_Prof_CloudJ, should be real8?
    REAL(fp) :: O3_CLIM (L1_  ) ! Computed in Set_Prof_CloudJ, should be real8?
    REAL(fp) :: AIR_CLIM(L1_  ) ! Computed in Set_Prof_CloudJ, should be real8?
    REAL(fp) :: Z_CLIM  (L1_+1) ! Computed in Set_Prof_CloudJ, should be real8?
    REAL(8)  :: HHH     (L1_  )
    REAL(8)  :: RRR     (L1_  )
    REAL(8)  :: CCC     (L1_  )
    REAL(8)  :: LWP     (L1_  )
    REAL(8)  :: IWP     (L1_  )
    REAL(8)  :: REFFL   (L1_  )
    REAL(8)  :: REFFI   (L1_  )
    REAL(8)  :: CLDF    (L1_  )

    ! 2D arrays
    INTEGER  :: NDXAER (L1_, AN_   )
    REAL(8)  :: AERSP  (L1_, AN_   )
    REAL(8)  :: RFL    (5  , W_+W_r)

    !------------------------------------------------------------------------
    ! Cloud_JX outputs
    !------------------------------------------------------------------------

    ! Scalars
    LOGICAL  :: LDARK
    INTEGER  :: NICA
    INTEGER  :: JCOUNT

    ! 1D arrays
    REAL(8)  :: SWMSQ (6   )
    REAL(8)  :: OD18  (L1_ )
    REAL(8)  :: WTQCA (NQD_)

    ! 2D arrays
    REAL(8)  :: SKPERD(S_+2, L1_)

    ! Which of the below is correct???
    REAL(8)  :: VALJXX(L_,JVN_)

    !------------------------------------------------------------------------
    ! For diagnostics
    !------------------------------------------------------------------------

    ! These are currently never set. Should they be output from Cloud-J?
    REAL(fp) :: FJBOT(W_)
    REAL(fp) :: FSBOT(W_)
    REAL(fp) :: FLXD(L1_,W_)
    REAL(fp) :: FJFLX(L_,W_)

    ! For UVFlux* diagnostics
    REAL(fp) :: FDIRECT (L1_)
    REAL(fp) :: FDIFFUSE(L1_)
    REAL(fp) :: UVX_CONST

    ! Species ids
    INTEGER, SAVE :: id_O3
    INTEGER, SAVE :: id_SO4

    ! Index for Cloud-J prints if GEOS-Chem verbose is on
    INTEGER :: I_PRT, J_PRT

    ! Debugging logicals to turn optical depth sources on/off
    LOGICAL :: use_liqcld
    LOGICAL :: use_icecld
    LOGICAL :: use_dust
    LOGICAL :: use_so4
    LOGICAL :: use_bc
    LOGICAL :: use_oc
    LOGICAL :: use_sala
    LOGICAL :: use_salc
    LOGICAL :: use_stratso4
    LOGICAL :: use_psc

    !=================================================================
    ! Run_CloudJ begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Run_CloudJ (in module GeosCore/cldj_interface_mod.F90)'

    ! Set debugging logicals to turn optical depth sources on/off. If using,
    ! uncomment where they are appear later in this file. They are commented out by
    ! default to avoid unnecessary slow-down.
    use_liqcld   = .true.
    use_icecld   = .true.
    use_dust     = .true.
    use_so4      = .true.
    use_bc       = .true.
    use_oc       = .true.
    use_sala     = .true.
    use_salc     = .true.
    use_stratso4 = .true.
    use_psc      = .true.

    ! Aerosol indexes (must match mapping set in RD_AOD)
    SO4_ind  = 1
    BC_ind   = 2
    OC_ind   = 3
    SALA_ind = 4
    SALC_ind = 5

    ! Relative humidities in FJX_spec-aer.dat
    RH_lut(1) = 0.d0
    RH_lut(2) = 50.d0
    RH_lut(3) = 70.d0
    RH_lut(4) = 80.d0
    RH_lut(5) = 90.d0

    ! Index for wavelength 1000 in optical property LUT
    ind_1000 = 5

    ! Diagnostic initialization
    IF ( State_Diag%Archive_UVFluxDiffuse ) State_Diag%UVFluxDiffuse = 0.0_f4
    IF ( State_Diag%Archive_UVFluxDirect  ) State_Diag%UVFluxDirect  = 0.0_f4
    IF ( State_Diag%Archive_UVFluxNet     ) State_Diag%UVFluxNet     = 0.0_f4
    IF ( State_Diag%Archive_OD600         ) State_Diag%OD600         = 0.0_f4
    IF ( State_Diag%Archive_TCOD600       ) State_Diag%TCOD600       = 0.0_f4
#if defined( MODEL_GEOS )
    ! TODO: implement these
    IF ( State_Diag%Archive_EXTRALNLEVS ) State_Diag%EXTRALNLEVS = 0.0
    IF ( State_Diag%Archive_EXTRALNITER ) State_Diag%EXTRALNITER = 0.0
#endif

    ! Set species ids for use in diagnostics
    IF ( FIRST ) THEN
       id_O3   = Ind_('O3')
       id_SO4  = Ind_('SO4')
       IF ( id_O3 < 0 ) THEN
          ErrMsg = 'O3 is not a defined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
 
    ! Constant values across grid boxes
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    DAY_OF_YEAR = GET_DAY_OF_YEAR()


    ! ewl: set NDXAER to MIEDX duplicated for all levels. If this works
    ! will want to store this elsewhere. Don't want to do this computation
    ! every timestep.
    ! Why isn't this integer???
    ! Should change MIEDX to be NDXAER for cloud-j?
    NDXAER(:,:) = 0.d0
    DO N = 1, AN_
    DO L = 1, L1_
       NDXAER(L,N) = State_Chm%Phot%MIEDX(N)
    ENDDO
    ENDDO

    !=================================================================
    ! For each column compute Cloud-J inputs and call Cloud_JX to compute J-values
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( A, I, J, L, K, N, S, MW_g, BoxHt, RH_ind                    ) &
    !$OMP PRIVATE( S_rh0, S_rhx, K_rh0, K_rhx, FRAC, RAA_eff, QAA_eff, SAA_eff ) &
    !$OMP PRIVATE( dry_to_wet_factor                                           ) &
    !$OMP PRIVATE( R_interp_factor, Q_interp_factor                            ) &
    !$OMP PRIVATE( U0, SZA, SOLF, T_CTM, P_CTM,  O3_CTM                        ) &
    !$OMP PRIVATE( T_CLIM, O3_CLIM, AIR_CLIM, Z_CLIM                           ) &
    !$OMP PRIVATE( CLDIW, CLDF, IWP, LWP, REFFI, REFFL, IWC, LWC, DELTA_P      ) &
    !$OMP PRIVATE( AERSP, RFL, RRR, LPRTJ, IRAN, CLDCOR, HHH, CCC              ) &
    !$OMP PRIVATE( LDARK, NICA, JCOUNT, SWMSQ, OD18, WTQCA, SKPERD, VALJXX     ) &
    !$OMP PRIVATE( FJBOT, FSBOT, FLXD, FJFLX, FDIRECT, FDIFFUSE, UVX_CONST     ) &
    !$OMP SCHEDULE( DYNAMIC )

    ! Loop over all latitudes and all longitudes
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Debug prints in Cloud-J. Limit to one grid cell so not excessive.
       ! Use this for debugging purposes only.
       LPRTJ = .false.
       I_PRT = 20
       J_PRT = 20
       IF ( I .eq. I_PRT .and. J .eq. J_PRT .and. Input_Opt%Verbose ) THEN
          print *, "cldj_interface_mod.F90: Cloud-J prints on for lat, lon: ", &
                   State_Grid%GlobalYMid(I,J), State_Grid%GlobalXMid(I,J)
          LPRTJ = .true.
       ENDIF

       !-----------------------------------------------------------------
       ! Solar zenith angle
       !-----------------------------------------------------------------

       ! Cosine of solar zenith angle [unitless]
       U0 = State_Met%SUNCOSmid(I,J)

       ! Solar zenith angle
       SZA = ACOS( MIN( MAX( U0, -1._fp ), 1._fp ) ) / PI_180

       ! Skip if dark conditions SZA > 98.0 deg => tan ht = 63 km
       if (SZA .gt. 98.e+0_fp) cycle

       ! Offset used for GEOS-Chem with fast-jx is 186; 172 in Cloud-J
       SOLF  = 1.e+0_fp - ( 0.034e+0_fp   &
               * cos( dble( DAY_OF_YEAR - 172 ) * 2.e+0_fp * PI / 365.e+0_fp ) )

       !-----------------------------------------------------------------
       ! Vertical climatology profiles
       !-----------------------------------------------------------------
       ! Cloud-J requires climatology vertical profiles for:
       !    Temperature               [K]
       !    Ozone                     [# O3 molec/cm2]
       !    Edge altitude             [cm]
       !    Path density (air column) [# molec/cm2]
       ! We compute these from CTM values per timestep

       ! Temperature profile [K]
       T_CTM(1:State_Grid%NZ) = State_Met%T(I,J,1:State_Grid%NZ)
       T_CTM(State_Grid%NZ+1) = T_CTM(State_Grid%NZ)

       ! Pressure profile [hPa]
       P_CTM(1:State_Grid%NZ+1) = State_Met%PEDGE(I,J,1:State_Grid%NZ+1)
       P_CTM(State_Grid%NZ+2) = State_Met%PEDGE(I,J,State_Grid%NZ+1) / 10.d0

       ! Ozone profile [molec/cm3]
       MaxLev = State_Met%ChemGridLev(I,J)
       O3_CTM = 0e+0_fp
       O3_CTM(1:MaxLev) = State_Chm%Species(id_O3)%Conc(I,J,1:MaxLev)

       ! Compute climatology. This subroutine is analogous to Cloud-J ACLIM_FJX.
       CALL Set_Clim_Profiles (I,          J,        MONTH,    DAY,      &
                               T_CTM,      P_CTM,    O3_CTM,             &
                               T_CLIM,     O3_CLIM,  Z_CLIM,   AIR_CLIM, &
                               Input_Opt,  State_Grid, State_Chm )

       !-----------------------------------------------------------------
       ! Clouds and humidity
       !-----------------------------------------------------------------

       CLDIW(:) = 0     ! Cloud type flag [0=none, 1=water, 2=ice, 3=both]
       CLDF(:)  = 0.d0  ! Cloud fraction [unitless]
       IWP(:)   = 0.d0  ! Ice cloud mass   [g/m2]
       LWP(:)   = 0.d0  ! Water cloud mass [g/m2]
       REFFI(:) = 0.d0  ! Ice cloud effective radius   [microns]
       REFFL(:) = 0.d0  ! Water cloud effective radius [microns]

       ! Set cloud fraction from input meteorology field
       CLDF(1:State_Grid%NZ) = State_Met%CLDF(I,J,1:State_Grid%NZ)
       CLDF(State_Grid%NZ+1) = CLDF(State_Grid%NZ)

       ! Loop over # layers in cloud-j (layers with clouds)
       DO L = 1, LWEPAR

          ! Get in-cloud liquid and ice water content from met-fields [kg/kg]
          LWC = State_Met%QL(I,J,L)
          IWC = State_Met%QI(I,J,L)

          ! Compute cloud type flag and reset cloud fraction if below threshold
          IF ( CLDF(L) .GT. 0.005d0 ) THEN
             IF ( LWC .GT. 1.d-11 ) CLDIW(L) = 1
             IF ( IWC .GT. 1.d-11 ) CLDIW(L) = CLDIW(L) + 2
          ELSE
             CLDF(L) = 0.d0
          ENDIF

          ! NOTES ON EFFECTIVE RADIUS FROM M. PRATHER:
          ! Compute effective radius [microns] of liquid water cloud and liquid ice cloud
          ! based on met-fields for in-cloud water content and in-cloud optical depth.
          !
          ! Note: The approach used here is consistent with the cloud optical depth
          ! calculation within Cloud-J but makes an assumption of extinction efficiency
          ! Q = 2.06. This works because all cloud Reffs are much bigger
          ! (2*pi*Reff >> 500 nm) so that Q (extinction efficiency = optical cross
          ! section / pi*r*r) is nearly constant at 2.06 (see the cloud scattering
          ! tables in Cloud-J).

          ! Compute liquid water path [g/m2] and effective radius [microns]
          DELTA_P = P_CTM(L) - P_CTM(L+1)
          IF ( State_Met%TAUCLW(I,J,L) .GT. 0.d0 ) THEN
             LWP(L) = 1000.d0 * LWC * DELTA_P * g0_100
             REFFL(L) = LWP(L) * 0.75d0 * 2.06d0 / ( State_Met%TAUCLW(I,J,L) * 1.d0 )
          ENDIF

          ! Compute ice water path [g/m2] and effective radius [microns]
          IF ( State_Met%TAUCLI(I,J,L) .GT. 0.d0 ) THEN
             IWP(L) = 1000.d0 * IWC * DELTA_P * g0_100
             REFFI(L) = IWP(L) * 0.75d0 * 2.06d0 / ( State_Met%TAUCLI(I,J,L) * 0.917d0 )
          ENDIF

          ! Convert relative humidity [unitless] from percent to fraction
          RRR(L) = State_Met%RH(I,J,L) / 100.d0

       ENDDO

       ! Set top of atmosphere relative humidity to 10% of layer below
       RRR(State_Grid%NZ+1) = RRR(State_Grid%NZ) * 1.d-1

       !-----------------------------------------------------------------
       ! Compute concentration per aerosol [g/m2]
       !-----------------------------------------------------------------
       ! AERSP is column concentration in g/m2 for each aerosol. The array currently
       ! includes entries for clouds but these are not used in Cloud-J and can be
       ! left as zero.
       ! Size is (L1_, AN_) where,
       !    AN_ is # of separate aerosols per layer (=37 for GEOS-Chem)
       !    L1_ is # of layer of edges (=73)

       ! Initialize conentration array to zero
       AERSP(:,:) = 0.d0

       ! Set values in loop over levels
       DO L= 1, State_Grid%NZ

          ! Layer height [m]
          BoxHt = State_Met%BXHEIGHT(I,J,L)   


          !---------------------------------------
          !  Non-aerosols in array
          !--------------------------------------
          ! Leave AERSP(L,1:3) as zero since non-aerosols (black carbon
          ! absorber, water cloud, and irregular ice cloud)

          !---------------------------------------
          !  Mineral dust [kg/m3] -> [g/m2]
          !--------------------------------------
          DO K = 4, 10
             AERSP(L,K) = State_Chm%SOILDUST(I,J,L,K-3) * BoxHt * 1.d3
          ENDDO

          !---------------------------------------
          ! Aerosols undergoing hydroscopic growth
          !--------------------------------------
          IF ( State_Met%InChemGrid(I,J,L) ) THEN

             ! For aerosols undergoing hygroscopic growth we need to pass the
             ! concentration that will be used in Cloud-J with humidity-dependent
             ! parameters. This means we need to convert from dry to wet concentration
             ! depending on this grid cell's humidity. We do this by computing an
             ! effective radius by linearly interpolating between LUT radius values
             ! for each aerosol given the current relative humidity and then applying
             ! a wet to dry conversion factor of ( Reff / Rdry )^^3 to the dry
             ! concentration. This is only done for the hydrophilic concentrations.
             ! Hydrophobic black and organic carbon are not converted to wet.
             !
             ! In addition to this we apply two conversion factors to take into account
             ! that Cloud-J does not interpolate extinction and effective radius between
             ! relative humidity entries in FJX_spec-aer.dat when computing optical
             ! depth. A conversion factor based on linear interpolation of the parameters
             ! is computed and applied to the concentration passed into Cloud-J. This
             ! effectively scales the extinction and radius within the optical depth
             ! calculation within Cloud-J. For consistency with Fast-JX previously
             ! implemented in GEOS-Chem we use values at 1000 nm for extinction interpolation.
             !
             ! We also separate the concentration array into 5 different arrays,
             ! one for each relative humidity entry in FJX_spec-aer.dat. Values are
             ! zero in each array except where current relative humidity falls within the
             ! pre-defined relative humidity range for each entry. For example,
             ! AERSP(L,11-15) are 5 values of sulfate for the same grid box. If RH is
             ! between 0 and 50 then only AERSP(L,11) is non-zero. AERSP(L,12:15) are
             ! all zero values.

             ! Humidity bin for aerosols (1:<=50, 2:<=70, 3:<=80; 4:<=90, else 5)
             IF ( State_Met%RH(I,J,L) <= 50 ) THEN
                RH_ind = 1
             ELSE IF ( State_Met%RH(I,J,L) <= 70 ) THEN
                RH_ind = 2
             ELSE IF ( State_Met%RH(I,J,L) <= 80 ) THEN
                RH_ind = 3
             ELSE IF ( State_Met%RH(I,J,L) <= 90 ) THEN
                RH_ind = 4
             ELSE
                RH_ind = 5
             ENDIF

             !----------------------------------------------------
             ! Sulfate [dry molec/cm3] -> [wet g/m2] (troposphere only)
             !----------------------------------------------------

             IF ( Input_Opt%LSULF .AND. State_Met%InTroposphere(I,J,L) ) THEN

                ! Get indexes to optical property LUT
                S_rh0 = 3 + NDUST + NRHAER*(SO4_ind-1) + 1  ! SO4 index for RH=0 in NDXAER
                S_rhx = S_rh0 + RH_ind - 1  ! Sulfate index for this RH
                K_rh0 = NDXAER(L,S_rh0)     ! index for RH=0 in FJX_spec-aer.dat
                K_rhx = NDXAER(L,S_rhx)     ! index for this RH in FJX_spec-aer.dat

                ! Get interpolated effective radius and extinction for RH in this grid box
                IF ( RH_ind == NRH ) THEN
                   RAA_eff = RAA(K_rhx)
                   QAA_eff = QAA(ind_1000,K_rhx)
                ELSE
                   FRAC = ( State_Met%RH(I,J,L) - RH_lut(RH_ind) ) &
                        / ( RH_lut(RH_ind+1) - RH_lut(RH_ind) )
                   RAA_eff = RAA(K_rhx) + FRAC * ( RAA(K_rhx+1) - RAA(K_rhx) )
                   QAA_eff = QAA(ind_1000,K_rhx) &
                        + FRAC * ( QAA(ind_1000,K_rhx+1) - QAA(ind_1000,K_rhx) )
                ENDIF
                dry_to_wet_factor = ( RAA_eff / RAA(K_rh0) )**3
                R_interp_factor = RAA_eff / RAA(K_rhx)
                Q_interp_factor = QAA_eff / QAA(ind_1000,K_rhx)

                ! Set concentration, converting [dry kg/m3] -> [wet g/m2]
                AERSP(L,S_rhx) = SO4_NH4_NIT(I,J,L) * BoxHt * 1.d3 * dry_to_wet_factor &
                     * Q_interp_factor / R_interp_factor

             ENDIF

             !----------------------------------------------------
             ! Carbon
             !----------------------------------------------------
             IF ( Input_Opt%LCARB ) THEN

                !----------------------------------------------------
                ! Black carbon
                !----------------------------------------------------

                ! Get indexes to optical property LUT
                S_rh0 = 3 + NDUST + NRHAER*(BC_ind-1) + 1  ! BC index for RH=0 in NDXAER
                S_rhx = S_rh0 + RH_ind - 1  ! BC index for this RH
                K_rh0 = NDXAER(L,S_rh0)     ! index for RH=0 in FJX_spec-aer.dat
                K_rhx = NDXAER(L,S_rhx)     ! index for this RH in FJX_spec-aer.dat

                ! Get interpolated effective radius and extinction for RH in this grid box
                IF ( RH_ind == NRH ) THEN
                   RAA_eff = RAA(K_rhx)          ! effective radius
                   QAA_eff = QAA(ind_1000,K_rhx) ! scattering phase function
                   SAA_eff = SAA(ind_1000,K_rhx) ! single scattering albedo
                ELSE
                   FRAC = ( State_Met%RH(I,J,L) - RH_lut(RH_ind) ) &
                        / ( RH_lut(RH_ind+1) - RH_lut(RH_ind) )
                   RAA_eff = RAA(K_rhx) + FRAC * ( RAA(K_rhx+1) - RAA(K_rhx) )
                   QAA_eff = QAA(ind_1000,K_rhx) + FRAC * ( QAA(ind_1000,K_rhx+1) - QAA(ind_1000,K_rhx) )
                   SAA_eff = SAA(ind_1000,K_rhx) + FRAC * ( SAA(ind_1000,K_rhx+1) - SAA(ind_1000,K_rhx) )
                ENDIF
                dry_to_wet_factor = ( RAA_eff / RAA(K_rh0) )**3
                R_interp_factor = RAA_eff / RAA(K_rhx)
                Q_interp_factor = QAA_eff / QAA(ind_1000,K_rhx)

                ! Set concentration
                IF ( Input_Opt%LBCAE ) THEN

                   ! Apply BC absorption enhancement (if using) first for hydrophilic BC
                   AERSP(L,S_rhx) = BCPI(I,J,L)                                      &
                        * ( Input_Opt%BCAE_1 + SAA_eff * (1.d0 - Input_Opt%BCAE_1) ) &
                        * dry_to_wet_factor * Q_interp_factor / R_interp_factor

                   ! Now apply hydrophobic using single scattering albedo for zero humidity
                   AERSP(L,S_rhx) = AERSP(L,S_rhx) + BCPO(I,J,L) &
                        * ( Input_Opt%BCAE_2 + SAA(ind_1000,K_rh0) * (1.d0 - Input_Opt%BCAE_2) )
                ELSE

                   ! No BC absorption enhancement
                   AERSP(L,S_rhx) = BCPO(I,J,L) &
                        + BCPI(I,J,L) * dry_to_wet_factor * Q_interp_factor / R_interp_factor

                ENDIF

                ! Convert to [dry kg/m3] -> [wet g/m2]
                AERSP(L,S_rhx) = AERSP(L,S_rhx) * 1.d3 * BoxHt

                !----------------------------------------------------
                ! Organic carbon
                !----------------------------------------------------

                ! Get indexes to optical property LUT
                S_rh0 = 3 + NDUST + NRHAER*(OC_ind-1) + 1  ! OC index for RH=0 in NDXAER
                S_rhx = S_rh0 + RH_ind - 1  ! OC index for this RH
                K_rh0 = NDXAER(L,S_rh0)     ! index for RH=0 in FJX_spec-aer.dat
                K_rhx = NDXAER(L,S_rhx)     ! index for this RH in FJX_spec-aer.dat

                ! Get interpolated effective radius and extinction for RH in this grid box
                IF ( RH_ind == NRH ) THEN
                   RAA_eff = RAA(K_rhx)
                   QAA_eff = QAA(ind_1000,K_rhx)
                ELSE
                   FRAC = ( State_Met%RH(I,J,L) - RH_lut(RH_ind) ) &
                        / ( RH_lut(RH_ind+1) - RH_lut(RH_ind) )
                   RAA_eff = RAA(K_rhx) + FRAC * ( RAA(K_rhx+1) - RAA(K_rhx) )
                   QAA_eff = QAA(ind_1000,K_rhx) + FRAC * ( QAA(ind_1000,K_rhx+1) - QAA(ind_1000,K_rhx) )
                ENDIF
                dry_to_wet_factor = ( RAA_eff / RAA(K_rh0) )**3
                R_interp_factor = RAA_eff / RAA(K_rhx)
                Q_interp_factor = QAA_eff / QAA(ind_1000,K_rhx)

                ! Set concentration, converting [dry kg/m3] -> [wet g/m2]
                AERSP(L,S_rhx) = ( OCPO(I,J,L)                                                      &
                     + ( OCPISOA(I,J,L) * dry_to_wet_factor * Q_interp_factor / R_interp_factor ) ) &
                     * 1.d3 * BoxHt

             ENDIF

             !----------------------------------------------------
             ! Seasalt [dry molec/cm3] -> [wet g/m2]
             !----------------------------------------------------

             IF ( Input_Opt%LSSALT ) THEN

                !----------------------------------------------------
                ! Accumulation mode seasalt
                !----------------------------------------------------

                ! Get indexes to optical property LUT
                S_rh0 = 3 + NDUST + NRHAER*(SALA_ind-1) + 1  ! SALA index for RH=0 in NDXAER
                S_rhx = S_rh0 + RH_ind - 1  ! SALA index for this RH
                K_rh0 = NDXAER(L,S_rh0)     ! index for RH=0 in FJX_spec-aer.dat
                K_rhx = NDXAER(L,S_rhx)     ! index for this RH in FJX_spec-aer.dat

                ! Get interpolated effective radius and extinction for RH in this grid box
                IF ( RH_ind == NRH ) THEN
                   RAA_eff = RAA(K_rhx)
                   QAA_eff = QAA(ind_1000,K_rhx)
                ELSE
                   FRAC = ( State_Met%RH(I,J,L) - RH_lut(RH_ind) ) &
                        / ( RH_lut(RH_ind+1) - RH_lut(RH_ind) )
                   RAA_eff = RAA(K_rhx) + FRAC * ( RAA(K_rhx+1) - RAA(K_rhx) )
                   QAA_eff = QAA(ind_1000,K_rhx) + FRAC * ( QAA(ind_1000,K_rhx+1) - QAA(ind_1000,K_rhx) )
                ENDIF
                dry_to_wet_factor = ( RAA_eff / RAA(K_rh0) )**3
                R_interp_factor = RAA_eff / RAA(K_rhx)
                Q_interp_factor = QAA_eff / QAA(ind_1000,K_rhx)

                ! Set concentration, converting [dry kg/m3] -> [wet g/m2]
                AERSP(L,S_rhx) = SALA(I,J,L) * BoxHt * 1.d3 * dry_to_wet_factor &
                     * Q_interp_factor / R_interp_factor

                !----------------------------------------------------
                ! Coarse seasalt
                !----------------------------------------------------

                ! Get indexes to optical property LUT
                S_rh0 = 3 + NDUST + NRHAER*(SALC_ind-1) + 1  ! SALC index for RH=0 in NDXAER
                S_rhx = S_rh0 + RH_ind - 1  ! SALC index for this RH
                K_rh0 = NDXAER(L,S_rh0)     ! index for RH=0 in FJX_spec-aer.dat
                K_rhx = NDXAER(L,S_rhx)     ! index for this RH in FJX_spec-aer.dat

                ! Get interpolated effective radius and extinction for RH in this grid box
                IF ( RH_ind == NRH ) THEN
                   RAA_eff = RAA(K_rhx)
                   QAA_eff = QAA(ind_1000,K_rhx)
                ELSE
                   FRAC = ( State_Met%RH(I,J,L) - RH_lut(RH_ind) ) &
                        / ( RH_lut(RH_ind+1) - RH_lut(RH_ind) )
                   RAA_eff = RAA(K_rhx) + FRAC * ( RAA(K_rhx+1) - RAA(K_rhx) )
                   QAA_eff = QAA(ind_1000,K_rhx) + FRAC * ( QAA(ind_1000,K_rhx+1) - QAA(ind_1000,K_rhx) )
                ENDIF
                dry_to_wet_factor = ( RAA_eff / RAA(K_rh0) )**3
                R_interp_factor = RAA_eff / RAA(K_rhx)
                Q_interp_factor = QAA_eff / QAA(ind_1000,K_rhx)

                ! Set concentration, converting [dry molec/cm3] -> [wet g/m2]
                AERSP(L,S_rhx) = SALC(I,J,L) * BoxHt  * 1.d3 * dry_to_wet_factor &
                     * Q_interp_factor / R_interp_factor

             ENDIF

          ENDIF

          !------------------------
          ! Stratospheric aerosols
          !------------------------

          MW_g = State_Chm%SpcData(id_SO4)%Info%MW_g

          ! Use sulfate concentration for stratospheric aerosols. Only set if the optical
          ! depth computed in GEOS-Chem is non-zero.

          !  SSA/LBS/STS
          IF ( State_Chm%Phot%ODAER(I,J,L,State_Chm%Phot%IWV1000,6) > 0._fp ) THEN
             AERSP(L,36) = State_Chm%Species(id_SO4)%Conc(I,J,L) &
                  * MW_g / AVO * BoxHt * 1e+6_fp
          ENDIF

          !  NAT/ice PSCs
          IF ( State_Chm%Phot%ODAER(I,J,L,State_Chm%Phot%IWV1000,7) > 0._fp ) THEN
             AERSP(L,37) = State_Chm%Species(id_SO4)%Conc(I,J,L) &
                  * MW_g / AVO * BoxHt * 1e+6_fp
          ENDIF

       ENDDO ! levels

       ! Set TOA equal to concentration in top level
       AERSP(State_Grid%NZ+1,:) = AERSP(State_Grid%NZ,:)

       ! Debugging option to set contributions from different sources to zero.
       ! Uncomment if using.
       !IF ( .NOT. use_liqcld   ) LWP(:)         = 0.d0
       !IF ( .NOT. use_icecld   ) IWP(:)         = 0.d0
       !IF ( .NOT. use_dust     ) AERSP(:,4:10)  = 0.d0
       !IF ( .NOT. use_so4      ) AERSP(:,11:15) = 0.d0
       !IF ( .NOT. use_bc       ) AERSP(:,16:20) = 0.d0
       !IF ( .NOT. use_oc       ) AERSP(:,21:25) = 0.d0
       !IF ( .NOT. use_sala     ) AERSP(:,26:30) = 0.d0
       !IF ( .NOT. use_salc     ) AERSP(:,31:35) = 0.d0
       !IF ( .NOT. use_stratso4 ) AERSP(:,36)    = 0.d0
       !IF ( .NOT. use_psc      ) AERSP(:,37)    = 0.d0

       !-----------------------------------------------------------------
       ! Set remaining inputs needed for Cloud_JX
       !-----------------------------------------------------------------

       ! UV surface albedo [unitless]
       ! Use same value for all levels and wavelengths
       RFL(1:5,:) = State_Met%UVALBEDO(I,J)

       ! Cloud correlation coefficient
       CLDCOR = 0.33

       ! Only used for CLDFLAG = 5
       IRAN = 1

       ! Required variables that are not used
       HHH = 0.d0
       CCC = 0.d0

       !-----------------------------------------------------------------
       ! Call Cloud_JX
       !-----------------------------------------------------------------
       
       ! Cloud_JX output list for easy reference:
       ! SKPERD, SWMSQ, OD18, NICA, JCOUNT, LDARK, WTQCA

       ! ewl debug
       IF ( LPRTJ ) THEN
          print *, "Calling Cloud_JX with the following inputs: "
          print *, " -> U0       : ", U0
          print *, " -> SZA      : ", SZA
          print *, " -> RFL      : ", RFL
          print *, " -> SOLF     : ", SOLF
          print *, " -> P_CTM    : ", P_CTM
          print *, " -> Z_CLIM   : ", Z_CLIM
          print *, " -> T_CLIM   : ", T_CLIM
          print *, " -> HHH      : ", HHH
          print *, " -> AIR_CLIM : ", AIR_CLIM
          print *, " -> RRR      : ", RRR
          print *, " -> O3_CLIM  : ", O3_CLIM
          print *, " -> CCC      : ", CCC
          print *, " -> LWP      : ", LWP
          print *, " -> IWP      : ", IWP
          print *, " -> REFFL    : ", REFFL
          print *, " -> REFFI    : ", REFFI
          print *, " -> CLDF     : ", CLDF
          print *, " -> CLDCOR   : ", CLDCOR
          print *, " -> CLDIW    : ", CLDIW
          print *, " -> AERSP    : ", AERSP
          print *, " -> IRAN     : ", IRAN
       ENDIF

       CALL Cloud_JX( U0,       SZA,      RFL,      SOLF,     LPRTJ,     &
                      P_CTM,    Z_CLIM,   T_CLIM,   HHH,      AIR_CLIM,  &
                      RRR,      O3_CLIM,  CCC,      LWP,      IWP,       &
                      REFFL,    REFFI,    CLDF,     CLDCOR,   CLDIW,     &
                      AERSP,    NDXAER,   L1_,      AN_,      JVN_,      &
                      VALJXX,   SKPERD,   SWMSQ,    OD18,     IRAN,      &
                      NICA,     JCOUNT,   LDARK,    WTQCA               )

       !-----------------------------------------------------------------
       ! Fill GEOS-Chem array ZPJ with J-values
       !-----------------------------------------------------------------
       DO L=1,State_Grid%MaxChemLev
          DO K=1,State_Chm%Phot%nPhotRxns
             IF (JIND(K).gt.0) THEN
                State_Chm%Phot%ZPJ(L,K,I,J) = VALJXX(L,JIND(K))*JFACTA(K)
             ELSE
                State_Chm%Phot%ZPJ(L,K,I,J) = 0.e+0_fp
             ENDIF
          ENDDO
       ENDDO

       ! Set J-rates outside the chemgrid to zero
       IF (State_Grid%MaxChemLev.lt.L_) THEN
          DO L=State_Grid%MaxChemLev+1,L_
             DO K=1,State_Chm%Phot%nPhotRxns
                State_Chm%Phot%ZPJ(L,K,I,J) = 0.e+0_fp
             ENDDO
          ENDDO
       ENDIF

       !-----------------------------------------------------------------
       ! Diagnostics for 600 nm optical depth computed in Cloud-J
       !-----------------------------------------------------------------
       IF ( State_Diag%Archive_OD600 ) THEN
          State_Diag%OD600(I,J,1:State_Grid%NZ) = OD18(1:State_Grid%NZ)
       ENDIF
       IF ( State_Diag%Archive_TCOD600 ) THEN
          State_Diag%TCOD600(I,J) = SUM(OD18(:))
       ENDIF

       !-----------------------------------------------------------------
       ! UV radiative flux diagnostics (direct, diffuse, net) [W/m2]
       ! Convention: negative is downwards
       !-----------------------------------------------------------------
       IF ( State_Diag%Archive_UVFluxDiffuse .or. &
            State_Diag%Archive_UVFluxDirect .or. &
            State_Diag%Archive_UVFluxNet ) THEN
       
          ! Loop over wavelength bins
          DO K = 1, W_
       
             ! Initialize
             FDIRECT  = 0.0_fp
             FDIFFUSE = 0.0_fp
       
             ! ewl: this is messed up. FSBOT and FJBOT aren't set.

             ! Direct & diffuse fluxes at each level
             FDIRECT(1)  = FSBOT(K)                    ! surface
             FDIFFUSE(1) = FJBOT(K)                    ! surface
             DO L = 2, State_Grid%NZ
                FDIRECT(L) = FDIRECT(L-1) + FLXD(L-1,K)
                FDIFFUSE(L) = FJFLX(L-1,K)
             ENDDO
       
             ! Constant to multiply UV fluxes at each wavelength bin
             UVX_CONST = SOLF * FL(K) * State_Chm%Phot%UVXFACTOR(K)
       
             ! Archive into diagnostic arrays
             DO L = 1, State_Grid%NZ
       
                IF ( State_Diag%Archive_UVFluxNet ) THEN
                   S = State_Diag%Map_UvFluxNet%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxNet(I,J,L,S) =  &
                      State_Diag%UVFluxNet(I,J,L,S) +  &
                           ( ( FDIRECT(L) + FDIFFUSE(L) ) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDirect ) THEN
                   S = State_Diag%Map_UvFluxDirect%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDirect(I,J,L,S) =  &
                      State_Diag%UVFluxDirect(I,J,L,S) +  &
                           ( FDIRECT(L) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
                   S = State_Diag%Map_UvFluxDiffuse%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDiffuse(I,J,L,S) =  &
                      State_Diag%UVFluxDiffuse(I,J,L,S) +  &
                           ( FDIFFUSE(L) * UVX_CONST )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Reset first-time flag
    FIRST=.FALSE.

  END SUBROUTINE Run_CloudJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_clim_profiles
!
! !DESCRIPTION: Subroutine SET\_CLIM_\PROFILES sets vertical climatology profiles
!  for a given latitude and longitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Clim_Profiles( ILON,       ILAT,       MONTH,   DAY,       &
                                T_CTM,      P_CTM,      O3_CTM,             &
                                T_CLIM,     O3_CLIM,    Z_CLIM,  AIR_CLIM,  &
                                Input_Opt,  State_Grid, State_Chm )
!
! !USES:
!
    USE Cldj_Cmn_Mod,    ONLY : L_, L1_, ZZHT
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PhysConstants,   ONLY : AIRMW, AVO, g0, BOLTZ
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Chm_Mod,   ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: ILON              ! Longitude index
    INTEGER,        INTENT(IN) :: ILAT              ! Latitude index
    INTEGER,        INTENT(IN) :: MONTH             ! Month
    INTEGER,        INTENT(IN) :: DAY               ! Day *of month*
    REAL(fp),       INTENT(IN) :: T_CTM(L1_)        ! CTM temperatures (K)
    REAL(fp),       INTENT(IN) :: P_CTM(L1_+1)      ! CTM edge pressures (hPa)
    REAL(fp),       INTENT(IN) :: O3_CTM(L1_)       ! CTM ozone (molec/cm3)
    TYPE(OptInput), INTENT(IN) :: Input_Opt         ! Input Options object
    TYPE(GrdState), INTENT(IN) :: State_Grid        ! Grid State object
    TYPE(ChmState), INTENT(IN) :: State_Chm         ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT)      :: T_CLIM(L1_)       ! Clim. temperatures (K)
    REAL(fp), INTENT(OUT)      :: Z_CLIM(L1_+1)     ! Edge altitudes (cm)
    REAL(fp), INTENT(OUT)      :: O3_CLIM(L1_)      ! O3 column depth (#/cm2)
    REAL(fp), INTENT(OUT)      :: AIR_CLIM(L1_)     ! O3 column depth (#/cm2)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - Adapted from Set_Prof_FastJX for use with Cloud-J
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: J, K, L, M
    REAL(fp) :: DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH
    REAL(fp) :: PSTD(52),OREF2(51),TREF2(51)
    REAL(fp) :: PROFCOL

    !=================================================================
    ! Set_Clim_Profiles begins here!
    !=================================================================

    !=================================================================
    ! Set up pressure levels for O3/T climatology - assume that value
    ! given for each 2 km z* level applies from 1 km below to 1 km
    ! above, so select pressures at these boundaries. Surface level
    ! values at 1000 mb are assumed to extend down to the actual
    ! surface pressure for this lat/lon.
    !=================================================================
    PSTD(1)  = MAX(P_CTM(1),1000.e+0_fp)
    PSTD(2)  = 1000.e+0_fp * 10.e+0_fp ** (-1.e+0_fp/16.e+0_fp)
    DLOGP    = 10.e+0_fp**(-2.e+0_fp/16.e+0_fp)
    DO K=3,51
       PSTD(K) = PSTD(K-1) * DLOGP
    ENDDO
    PSTD(52) = 0.e+0_fp

    ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
    MASFAC = 100.e+0_fp * AVO / ( AIRMW * g0 * 10.e+0_fp )

    ! Select appropriate monthly and latitudinal profiles
    ! Now use State_Grid%YMid instead of Oliver's YDGRD(NSLAT)
    M = MAX( 1, MIN( 12, MONTH ) )
    J = MAX( 1, MIN( 18, ( INT(State_Grid%YMid(ILON,ILAT)) + 99 ) / 10 ) )

    ! Temporary arrays for climatology data
    DO K = 1, 51
       OREF2(K) = State_Chm%Phot%OREF(K,J,M)
       TREF2(K) = State_Chm%Phot%TREF(K,J,M)
    ENDDO

    ! Apportion O3 and T on supplied climatology z* levels onto CTM levels
    ! with mass (pressure) weighting, assuming constant mixing ratio and
    ! temperature half a layer on either side of the point supplied.
    DO L = 1, L1_
       F0 = 0.e+0_fp
       T0 = 0.e+0_fp
       DO K = 1, 51
          PC = MIN( P_CTM(L),   PSTD(K)   )
          PB = MAX( P_CTM(L+1), PSTD(K+1) )
          IF ( PC .GT. PB ) THEN
             XC = ( PC - PB ) / ( P_CTM(L) - P_CTM(L+1) )
             F0 = F0 + OREF2(K)*XC
             T0 = T0 + TREF2(K)*XC
          ENDIF
       ENDDO
       T_CLIM(L)  = T0
       O3_CLIM(L) = F0 * 1.e-6_fp
    ENDDO

    !=================================================================
    ! Calculate effective altitudes using scale height at each level
    !=================================================================
    Z_CLIM(1) = 0.e+0_fp
    DO L = 1, L_
       SCALEH = BOLTZ * 1.e+4_fp * MASFAC * T_CLIM(L)
       Z_CLIM(L+1) = Z_CLIM(L) - ( LOG( P_CTM(L+1) / P_CTM(L) ) * SCALEH )
    ENDDO
    Z_CLIM(L1_+1)=Z_CLIM(L1_) + ZZHT

    !=================================================================
    ! Calculate column quantities for Cloud-J
    !=================================================================
    PROFCOL = 0e+0_fp

    DO L = 1, L1_

       ! Monthly mean air Column [molec/cm2]
       AIR_CLIM(L)  = ( P_CTM(L) - P_CTM(L+1) ) * MASFAC

       ! Monthly mean O3 column [molec/cm2]
       O3_CLIM(L) = O3_CLIM(L) * AIR_CLIM(L)

       ! Monthly mean O3 column [DU]
       PROFCOL = PROFCOL + ( O3_CLIM(L) / 2.69e+16_fp )
    ENDDO

!ewl: is this needed?
    !! Top values are special (do not exist in CTM data)
    !AIR_CLIM(L1_)     = P_CTM(L1_) * MASFAC
    !O3_CLIM(L1_) = O3_CLIM(L1_) * AIR_CLIM(L1_)

    ! Scale monthly O3 profile to the daily O3 profile (if available)
    DO L = 1, L1_

       ! Use online O3 values in the chemistry grid if selected; otherwise use
       ! O3 values from the met fields or TOMS/SBUV
       IF ( ( Input_opt%Use_Online_O3 )             &
            .AND. ( L <= State_Grid%MaxChemLev )    &
            .AND. ( O3_CTM(L) > 0e+0_fp ) ) THEN

          ! Convert from molec/cm3 to molec/cm2
          O3_CLIM(L) = O3_CTM(L) * (Z_CLIM(L+1)-Z_CLIM(L))

       ELSEIF (State_Chm%TO3_Daily(ILON,ILAT) > 0e+0_fp) THEN

          ! NOTE: replaced O3_TOMS with State_Chm%TO3_Daily since is the equivalent
          O3_CLIM(L) = O3_CLIM(L) * ( State_Chm%TO3_Daily(ILON,ILAT) / PROFCOL )

       ENDIF

    ENDDO

  END SUBROUTINE Set_Clim_Profiles
!EOC
END MODULE CLDJ_INTERFACE_MOD
