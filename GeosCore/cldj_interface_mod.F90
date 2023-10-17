#ifdef CLOUDJ
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
#if defined( MODEL_CESM ) && defined( SPMD )
  USE MPISHORTHAND
  USE SPMD_UTILS
#endif

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
  SUBROUTINE INIT_CLOUDJ( Input_Opt, State_Diag, State_Chm, RC )
!
! !USES:
!

! ewl: Use, inputs/outputs, and local vars could be slimmed down
    ! ewl: if these are in cloud-j, why do I need to pass them???
    USE Cldj_Cmn_Mod,   ONLY : JVN_, NRatJ, W_
!ewl    USE Cldj_Cmn_Mod,   ONLY : JVN_, NJX, NRATJ, W_, WL
!ewl    USE Cldj_Cmn_Mod,   ONLY : TITLEJX, JLABEL, JFACTA, RNAMES
!ewl    USE Cldj_Cmn_Mod,   ONLY : JIND, BRANCH     ! ewl debugging
    USE Cldj_Init_Mod,  ONLY : Init_CldJ
    USE ErrCode_Mod
!ewl    USE ERROR_MOD,      ONLY : ERROR_STOP
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)     :: Input_Opt   ! Input Options object
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
    CALL Init_CldJ(Input_Opt%CloudJ_Dir, TITLEJXX, JVN_, NJXX)

    ! Store # of photolysis reactions in State_Chm object
    State_Chm%Phot%nPhotRxns = NRatJ

! ewl: do we need this?? state_grid not passed here.
!    ! Random error handling??? Is this needed?
!    if (State_Grid%NZ+1 .gt. JXL1_) then
!       CALL ERROR_STOP('Init_CloudJ: not enough levels in JX', &
!                       'cldj_interface_mod.F90' )
!    endif

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
    ! ewl: if these are in cloud-j, why do we need them here??
    USE Cldj_Cmn_Mod,   ONLY : L_, L1_, W_, S_, LWEPAR
    USE Cldj_Cmn_Mod,   ONLY : JVN_, JXL_, JXL1_, AN_, NQD_, W_r
    USE Cldj_Cmn_Mod,   ONLY : JIND, JFACTA, FL 
    USE Cld_Sub_Mod,    ONLY : Cloud_JX
    USE Cmn_Size_Mod,   ONLY : NRHAER
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
    INTEGER            :: A, I, J, L, K, N, S, MaxLev, RH_bin
    REAL(8)            :: MW_kg, BoxHt, Delta_P, IWC, LWC
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
    REAL(fp) :: FLXD(JXL1_,W_)
    REAL(fp) :: FJFLX(JXL_,W_)

    ! For UVFlux* diagnostics
    REAL(fp) :: FDIRECT (JXL1_)
    REAL(fp) :: FDIFFUSE(JXL1_)
    REAL(fp) :: UVX_CONST

    ! Species ids
    INTEGER, SAVE :: id_O3
    INTEGER, SAVE :: id_SO4
    INTEGER, SAVE :: id_NH4
    INTEGER, SAVE :: id_NIT
    INTEGER, SAVE :: id_BCPI
    INTEGER, SAVE :: id_BCPO
    INTEGER, SAVE :: id_OCPI
    INTEGER, SAVE :: id_OCPO
    INTEGER, SAVE :: id_SALA
    INTEGER, SAVE :: id_SALC

    ! Index for Cloud-J prints if GEOS-Chem verbose is on
    INTEGER :: I_PRT, J_PRT

    !=================================================================
    ! Run_CloudJ begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Run_CloudJ (in module GeosCore/cldj_interface_mod.F90)'

    ! Diagnostic initialization
    IF ( State_Diag%Archive_UVFluxDiffuse ) State_Diag%UVFluxDiffuse = 0.0_f4
    IF ( State_Diag%Archive_UVFluxDirect  ) State_Diag%UVFluxDirect  = 0.0_f4
    IF ( State_Diag%Archive_UVFluxNet     ) State_Diag%UVFluxNet     = 0.0_f4
    IF ( State_Diag%Archive_OD600         ) State_Diag%OD600         = 0.0_f4
    IF ( State_Diag%Archive_TCOD600       ) State_Diag%TCOD600       = 0.0_f4
#if defined( MODEL_GEOS )
    ! ewl: should these diags be set later? They are not.
    IF ( State_Diag%Archive_EXTRALNLEVS ) State_Diag%EXTRALNLEVS = 0.0
    IF ( State_Diag%Archive_EXTRALNITER ) State_Diag%EXTRALNITER = 0.0
#endif

    ! Set species ids for use in diagnostics. Require ozone.
    IF ( FIRST ) THEN
       id_O3   = Ind_('O3')
       id_SO4  = Ind_('SO4')
       id_NH4  = Ind_('NH4')
       id_NIT  = Ind_('NIT')
       id_BCPI = Ind_('BCPI')
       id_BCPO = Ind_('BCPO')
       id_OCPI = Ind_('OCPI')
       id_OCPO = Ind_('OCPO')
       id_SALA = Ind_('SALA')
       id_SALC = Ind_('SALC')
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
    !$OMP PRIVATE( A, I, J, L, K, N, S, MW_kg, BoxHt, RH_bin               ) &
    !$OMP PRIVATE( U0, SZA, SOLF, T_CTM, P_CTM,  O3_CTM                    ) &
    !$OMP PRIVATE( T_CLIM, O3_CLIM, AIR_CLIM, Z_CLIM                       ) &
    !$OMP PRIVATE( CLDIW, CLDF, IWP, LWP, REFFI, REFFL, IWC, LWC, DELTA_P  ) &
    !$OMP PRIVATE( AERSP, RFL, RRR, LPRTJ, IRAN, CLDCOR, HHH, CCC          ) &
    !$OMP PRIVATE( LDARK, NICA, JCOUNT, SWMSQ, OD18, WTQCA, SKPERD, VALJXX ) &
    !$OMP PRIVATE( FJBOT, FSBOT, FLXD, FJFLX, FDIRECT, FDIFFUSE, UVX_CONST ) &
    !$OMP SCHEDULE( DYNAMIC )

    ! Loop over all latitudes and all longitudes
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Debug prints in Cloud-J. Limit to one grid cell so not excessive.
       ! Use this for debugging purposes only.
       LPRTJ = .false.
       I_PRT = 1
       J_PRT = 10
       IF ( I .eq. I_PRT .and. J .eq. J_PRT .and. Input_Opt%Verbose ) THEN
          print *, "cldj_interface_mod.F90: Dloud-J prints on for lat, lon: ", &
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
       DO L=1,LWEPAR

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
       ! AERSP is concentration in g/m2. Size is L1_, AN_. Will need to
       ! compute this for each lat/lon in this loop. Could be a separate
       ! subroutine, or in set_prof. Best put in set_prof.
       ! AN_ is # of separate aerosols per layer, = 37.
       ! L1_ is # of layer of edges, so 73. Make the top one the same as 72?
       ! Just do that for now, but need to check if okay.

       ! Initialize to zero
       AERSP(:,:) = 0.d0

       ! Set values in loop over levels
       DO L=1,State_Grid%NZ

          ! Layer height [m]
          BoxHt = State_Met%BXHEIGHT(I,J,L)   

          ! Humidity bin for aerosols (1:<=50, 2:<=70, 3:<=80; 4:<=90, else 5
          IF ( State_Met%RH(I,J,L) <= 50 ) THEN
             RH_bin = 1
          ELSE IF ( State_Met%RH(I,J,L) <= 70 ) THEN
             RH_bin = 2
          ELSE IF ( State_Met%RH(I,J,L) <= 80 ) THEN
             RH_bin = 3
          ELSE IF ( State_Met%RH(I,J,L) <= 90 ) THEN
             RH_bin = 4
          ELSE
             RH_bin = 5
          ENDIF
          ! Leave AERSP(L,1:3) as zero since non-aerosols (black carbon
          ! absorber, water cloud, and irregular ice cloud)

          ! Mineral dust [kg/m3] -> [g/m2]
          DO K = 4, 10
             AERSP(L,K)  = State_Chm%SOILDUST(I,J,L,K-3)*BoxHt*1.d3 
          ENDDO

          ! Sulfate [molec/cm3] -> [g/m2]
          MW_kg = State_Chm%SpcData(id_SO4)%Info%MW_g * 1.e-3_fp
          AERSP(L,11:15) =                                            &
              ( State_Chm%Species(id_SO4)%Conc(I,J,L)                 &
              * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(I,J,L) ) &
              * BoxHt

          ! Black carbon [molec/cm3] -> [g/m2]
          MW_kg = State_Chm%SpcData(id_BCPI)%Info%MW_g * 1.e-3_fp
          AERSP(L,16:20) =                                            &
              ( ( State_Chm%Species(id_BCPI)%Conc(I,J,L) +            &
                State_Chm%Species(id_BCPO)%Conc(I,J,L) )              &
              * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(I,J,L) ) &
              * BoxHt

          ! Organic carbon [molec/cm3] -> [g/m2]
          MW_kg = State_Chm%SpcData(id_OCPI)%Info%MW_g * 1.e-3_fp
          AERSP(L,21:25) =                                            &
              ( ( State_Chm%Species(id_OCPI)%Conc(I,J,L) +            &
                State_Chm%Species(id_OCPO)%Conc(I,J,L) )              &
              * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(I,J,L) ) &
              * BoxHt

          ! Seasalt (accum) [molec/cm3] -> [g/m2]
          MW_kg = State_Chm%SpcData(id_SALA)%Info%MW_g * 1.e-3_fp
          AERSP(L,26:30) =                                            &
              ( State_Chm%Species(id_SALA)%Conc(I,J,L)                &
              * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(I,J,L) ) &
              * BoxHt

          ! Seasalt (coarse) [molec/cm3] -> [g/m2]
          MW_kg = State_Chm%SpcData(id_SALC)%Info%MW_g * 1.e-3_fp
          AERSP(L,31:35) =                                            &
              ( State_Chm%Species(id_SALC)%Conc(I,J,L)                &
              * 1e+6_fp / ( AVO / MW_kg ) / State_Met%AIRDEN(I,J,L) ) &
              * BoxHt
          
          ! Use sulfate concentration for sulfate stratospheric aerosols
          AERSP(L,36) = AERSP(L,11)
          AERSP(L,37) = AERSP(L,11)

       ENDDO ! levels

       ! Set TOA equal to concentration in top level. Not sure if this is
       ! right. We set TOA optical depth to zero when passing to photo_jx
       ! in old fast-jx. 
       AERSP(State_Grid%NZ+1,:) = AERSP(State_Grid%NZ,:)

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
       ! SKPERD
       ! SWMSQ
       ! OD18
       ! NICA
       ! JCOUNT
       ! LDARK
       ! WTQCA

       ! ewl debug
       IF ( LPRTJ ) THEN
          print *, "ewl: Calling Cloud_JX with the following inputs: "
          print *, " -> U0       : ", U0         ! cldj 1.0
          print *, " -> SZA      : ", SZA        ! cldj 0.0
          print *, " -> RFL      : ", RFL        ! looks ok 5e-2 vs 6e-2
          print *, " -> SOLF     : ", SOLF       ! looks ok 1 vs 0.97
          print *, " -> P_CTM    : ", P_CTM      ! looks ok
          print *, " -> Z_CLIM   : ", Z_CLIM     ! looks ok. gc sfc is 0, vs 171
          print *, " -> T_CLIM   : ", T_CLIM     ! looks ok
          print *, " -> HHH      : ", HHH        ! not used
          print *, " -> AIR_CLIM : ", AIR_CLIM   ! looks ok
          print *, " -> RRR      : ", RRR        ! FIXED
          print *, " -> O3_CLIM  : ", O3_CLIM    ! looks ok
          print *, " -> CCC      : ", CCC        ! not used
          print *, " -> LWP      : ", LWP        ! LOOKS WRONG! cldj 0-15, gc 1e10 where not zero.??
          print *, " -> IWP      : ", IWP        ! SAME ISSUE AS LWP!!!
          print *, " -> REFFL    : ", REFFL      ! TOTALLY OFF!! GC is huge!
          print *, " -> REFFI    : ", REFFI      ! SAME ISSUE AS REFFL!!!
          print *, " -> CLDF     : ", CLDF       ! presumably ok. highly variable
          print *, " -> CLDCOR   : ", CLDCOR     ! ok
          print *, " -> CLDIW    : ", CLDIW      ! looks like bug in cldj???
          print *, " -> AERSP    : ", AERSP      ! all zeros in cldj??
          print *, " -> IRAN     : ", IRAN       ! ok
       ENDIF

       ! ewl: deleted arguments CLDFLAG, NRANDO, and LNRG since globally
       ! set in cloud-j init via reading CJ77_inputs file
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

       ! Use online O3 values in the chemistry grid if selected
       IF ( ( Input_opt%Use_Online_O3 ) .and. &
            (L <= State_Grid%MaxChemLev) .and. &
            (O3_CTM(L) > 0e+0_fp) ) THEN

          ! Convert from molec/cm3 to molec/cm2
          O3_CLIM(L) = O3_CTM(L) * (Z_CLIM(L+1)-Z_CLIM(L))

       ! Otherwise, use O3 values from the met fields or TOMS/SBUV
       ELSEIF (State_Chm%TO3_Daily(ILON,ILAT) > 0e+0_fp) THEN

!ewl: note I replaced O3_TOMS with State_Chm%TO3_Daily since should equiv.
          O3_CLIM(L) = O3_CLIM(L) * ( State_Chm%TO3_Daily(ILON,ILAT) / PROFCOL )

       ENDIF

    ENDDO

  END SUBROUTINE Set_Clim_Profiles
!EOC
END MODULE CLDJ_INTERFACE_MOD
#endif
