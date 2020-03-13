!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: convection_mod.F90
!
! !DESCRIPTION: Module CONVECTION\_MOD contains routines which select the
!  proper convection code for different met field data sets.
!\\
!\\
! !INTERFACE:
!
MODULE CONVECTION_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DO_CONVECTION
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DO_CLOUD_CONVECTION
!
! !REVISION HISTORY:
!  27 Jan 2004 - R. Yantosca - Initial version
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
! !IROUTINE: do_convection
!
! !DESCRIPTION: Subroutine DO\_CONVECTION calls the appropriate convection q
!  driver program for different met field data sets.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_CONVECTION( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
!
! !USES:
!
    USE Diagnostics_Mod, ONLY : Compute_Column_Mass
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PhysConstants
    USE Species_Mod,     ONLY : Species
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : GET_TS_CONV
    USE UnitConv_Mod
    USE WETSCAV_MOD,     ONLY : COMPUTE_F
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  25 May 2005 - S. Wu - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, NA, nAdvect, NW, EC, ISOL
    INTEGER            :: I, J, L, NN, TS_DYN
    REAL(fp)           :: AREA_M2, DT
    LOGICAL            :: DO_ND14, DoConvFlux
    LOGICAL            :: DO_ND38, DoWetLoss
    REAL(fp)           :: DT_Conv

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Arrays
    REAL(fp)           :: PEDGE (State_Grid%NZ+1)
    REAL(fp)           :: DIAG14(State_Grid%NZ,State_Chm%nAdvect)
    REAL(fp)           :: DIAG38(State_Grid%NZ,State_Chm%nWetDep)
    REAL(fp)           :: F     (State_Grid%NZ,State_Chm%nAdvect)
    REAL(fp), TARGET   :: FSOL  (State_Grid%NX,State_Grid%NY,&
                                 State_Grid%NZ,State_Chm%nAdvect)

    ! Pointers
    REAL(fp), POINTER  :: p_FSOL (:,:,:)

    !=================================================================
    ! DO_CONVECT begins here!
    !=================================================================

    !-----------------------------------------------------------------
    ! Initialize
    !----------------------------------------------------------------
    RC      = GC_SUCCESS
    p_FSOL  => NULL()
    ErrMsg  = ''
    ThisLoc = ' -> at Do_Convection (in module GeosCore/convection_mod.F)'

    !----------------------------------------------------------
    ! Convection budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetConvection ) THEN
       ! Get initial column masses
       CALL Compute_Column_Mass( Input_Opt,                               &
                                 State_Chm,                               &
                                 State_Grid,                              &
                                 State_Met,                               &
                                 State_Chm%Map_Advect,                    &
                                 State_Diag%Archive_BudgetConvectionFull, &
                                 State_Diag%Archive_BudgetConvectionTrop, &
                                 State_Diag%Archive_BudgetConvectionPBL,  &
                                 State_Diag%BudgetMass1,                  &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Convection budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! More initializations
    !-----------------------------------------------------------------
    TS_DYN     = GET_TS_DYN()                           ! Dyn timestep [sec]
    DT         = DBLE( TS_DYN )                         ! Dyn timestep [sec]
    FSOL       = 0e+0_fp                                ! Zero the FSOL array
    DoConvFlux = State_Diag%Archive_CloudConvFlux       ! Save mass flux?
    DoWetLoss  = State_Diag%Archive_WetLossConv         ! Save wet loss?

    ! Number of advected species
    nAdvect = State_Chm%nAdvect

    !=================================================================
    ! Compute fraction of soluble species lost in conv. updrafts
    !=================================================================

    ! Loop over advected species
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( NA, N, p_FSOL, EC, ISOL )
    DO NA = 1, nAdvect

       ! Species ID
       N = State_Chm%Map_Advect(NA)

       ! Now point to a 3D slice of the FSOL array
       p_FSOL => FSOL(:,:,:,NA)

       ! Fraction of soluble species
       CALL COMPUTE_F( N, p_FSOL, ISOL, Input_Opt, State_Chm, &
                       State_Grid, State_Met, RC=EC )

       ! Trap potential errors (we can't exit an OpenMP loop)
       IF ( EC /= GC_SUCCESS ) THEN
          RC = EC
       ENDIF

       !--------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Fraction of soluble species lost in convective updrafts
       !--------------------------------------------------------------
       IF ( State_Diag%Archive_WetLossConvFrac .and. ISOL > 0 ) THEN
          State_Diag%WetLossConvFrac(:,:,:,ISOL) = p_FSOL
       ENDIF

       ! Free pointer memory
       p_FSOL => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Return if COMPUTE_F returned an error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Compute_F"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Do convection column by column
    !
    ! NOTE: Later on, consider moving the I,J loops within the call
    ! to DO_CLOUD_CONVECTION, to gain computational efficiency
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( J,      I,      AREA_M2, L, F,  EC ) &
    !$OMP PRIVATE( DIAG14, DIAG38, RC,      N, NA, NW ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! PRIVATE error trapping variable
       EC = GC_SUCCESS

       ! Skip buffer zone for nested grids (lzh, 4/1/15)
       IF ( State_Grid%NestedGrid ) THEN
          IF ( J <=                 State_Grid%SouthBuffer ) CYCLE
          IF ( J >  State_Grid%NY - State_Grid%NorthBuffer ) CYCLE
          IF ( I <=                 State_Grid%WestBuffer  ) CYCLE
          IF ( I >  State_Grid%NX - State_Grid%EastBuffer  ) CYCLE
       ENDIF

       ! Grid box surface area [m2]
       AREA_M2 =  State_Grid%Area_M2(I,J)

       ! NOTE: For some reason the code chokes when we try
       ! to use a pointer, so we'll use a 2D array here.
       F       =  FSOL(I,J,:,:)

       !--------------------------
       ! Do the cloud convection
       !--------------------------
       CALL DO_CLOUD_CONVECTION( Input_Opt  = Input_Opt,  &
                                 State_Chm  = State_Chm,  &
                                 State_Diag = State_Diag, &
                                 State_Grid = State_Grid, &
                                 State_Met  = State_Met,  &
                                 I          = I,          &
                                 J          = J,          &
                                 AREA_M2    = AREA_M2,    &
                                 F          = F,          &
                                 TS_DYN     = DT,         &
                                 USE_DIAG14 = DoConvFlux, &
                                 DIAG14     = DIAG14,     &
                                 USE_DIAG38 = DoWetLoss,  &
                                 DIAG38     = DIAG38,     &
                                 RC         = EC          )

       ! Trap potential errors (we can't exit an OpenMP loop)
       IF ( EC /= GC_SUCCESS ) THEN
          RC = EC
       ENDIF

       !--------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Convective mass flux [kg/s]
       ! NOTE: May be replaced soon with better flux diagnostics
       !--------------------------------------------------------------
       IF ( State_Diag%Archive_CloudConvFlux ) THEN
          DO N = 1, State_Chm%nAdvect
          DO L = 1, State_Grid%NZ
             State_Diag%CloudConvFlux(I,J,L,N) = Diag14(L,N)
          ENDDO
          ENDDO
       ENDIF

       !--------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Loss of soluble species in convective mass flux [kg/s]
       ! NOTE: May be replaced soon with better flux diagnostics
       !--------------------------------------------------------------
       IF ( State_Diag%Archive_WetLossConv ) THEN
          DO NW = 1, State_Chm%nWetDep
          DO L  = 1, State_Grid%NZ
             State_Diag%WetLossConv(I,J,L,NW) = Diag38(L,NW)
          ENDDO
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Return if COMPUTE_F returned an error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Do_Cloud_Convection"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Convection timestep [s]
    DT_Conv = GET_TS_CONV()

    !----------------------------------------------------------
    ! Convection budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetConvection ) THEN
       ! Get final masses and compute diagnostics
       CALL Compute_Column_Mass( Input_Opt,                               &
                                 State_Chm,                               &
                                 State_Grid,                              &
                                 State_Met,                               &
                                 State_Chm%Map_Advect,                    &
                                 State_Diag%Archive_BudgetConvectionFull, &
                                 State_Diag%Archive_BudgetConvectionTrop, &
                                 State_Diag%Archive_BudgetConvectionPBL,  &
                                 State_Diag%BudgetMass2,                  &
                                 RC )
       CALL Compute_Budget_Diagnostics( State_Grid,                       &
                                 State_Chm%Map_Advect,                    &
                                 DT_Conv,                                 &
                                 State_Diag%Archive_BudgetConvectionFull, &
                                 State_Diag%Archive_BudgetConvectionTrop, &
                                 State_Diag%Archive_BudgetConvectionPBL,  &
                                 State_Diag%BudgetConvectionFull,         &
                                 State_Diag%BudgetConvectionTrop,         &
                                 State_Diag%BudgetConvectionPBL,          &
                                 State_Diag%BudgetMass1,                  &
                                 State_Diag%BudgetMass2,                  &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Convection budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE DO_CONVECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_cloud_convection
!
! !DESCRIPTION: Subroutine DO\_CLOUD\_CONVECTION (formerly called NFCLDMX)
!  is S-J Lin's cumulus transport module for 3D GSFC-CTM, modified for the 
!  GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_CLOUD_CONVECTION( Input_Opt,  &
                                  State_Chm,  &
                                  State_Diag, &
                                  State_Grid, &
                                  State_Met,  &
                                  I,          &
                                  J,          &
                                  AREA_M2,    &
                                  F,          &
                                  TS_DYN,     &
                                  USE_DIAG14, &
                                  DIAG14,     &
                                  USE_DIAG38, &
                                  DIAG38,     &
                                  RC          )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : IT_IS_NAN
    USE ERROR_MOD,          ONLY : IT_IS_FINITE
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE Species_Mod,        ONLY : Species
    USE WETSCAV_MOD,        ONLY : WASHOUT
    USE WETSCAV_MOD,        ONLY : LS_K_RAIN
    USE WETSCAV_MOD,        ONLY : LS_F_PRIME
    USE WETSCAV_MOD,        ONLY : CONV_F_PRIME
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_SNOWPACK
    USE DEPO_MERCURY_MOD,   ONLY : ADD_Hg2_WD
    USE DEPO_MERCURY_MOD,   ONLY : ADD_HgP_WD
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: I, J        ! Lon & lat indices
    REAL(fp),       INTENT(IN)    :: AREA_M2     ! Surface area [m2]
    REAL(fp),       INTENT(IN)    :: F(:,:)      ! Fraction of soluble species
                                                 !  for updraft scavenging
                                                 !  [unitless].  Computed by
                                                 !  routine  COMPUTE_F.
    REAL(fp),       INTENT(IN)    :: TS_DYN      ! Dynamic timestep [sec]
    LOGICAL,        INTENT(IN)    :: USE_DIAG14  ! Archive DIAG14?
    LOGICAL,        INTENT(IN)    :: USE_DIAG38  ! Archive DIAG38?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: DIAG14(:,:) ! Array for ND14 diagnostic
    REAL(fp),       INTENT(OUT)   :: DIAG38(:,:) ! Array for ND38 diagnostic
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  Reference:
!  ============================================================================
!  Lin, SJ.  "Description of the parameterization of cumulus transport
!     in the 3D Goddard Chemistry Transport Model, NASA/GSFC, 1996.
!                                                                             .
!  Unit conversion for BMASS:
!
!      Ps - Pt (mb)| P2 - P1 | 100 Pa |  s^2  | 1  |  1 kg        kg
!     -------------+---------+--------+-------+----+--------  =  -----
!                  | Ps - Pt |   mb   | 9.8 m | Pa | m^2 s^2      m^2
!
!                                                                             .
!  NOTE: We are passing I & J down to this routine so that it can call the
!  proper code from "mercury_mod.f".  Normally, we wouldn't pass I & J as
!  arguments to columnized code.  This prevents rewriting the mercury_mod.f
!  routines ADD_Hg2_
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER    :: TINYNUM = 1e-14_fp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: AER,         IS_Hg
    INTEGER                :: IC,          ISTEP,     K
    INTEGER                :: KTOP,        NC,        NDT
    INTEGER                :: NLAY,        NS,        CLDBASE
    INTEGER                :: Hg_Cat,      NA,        nAdvect, NW
    REAL(fp)               :: CMFMC_BELOW, ALPHA,     ALPHA2
    REAL(fp)               :: CMOUT,       DELQ,      DQ
    REAL(fp)               :: DNS,         ENTRN,     QC
    REAL(fp)               :: QC_PRES,     QC_SCAV,   SDT
    REAL(fp)               :: T0,          T0_SUM,    T1
    REAL(fp)               :: T2,          T3,        T4
    REAL(fp)               :: TSUM,        LOST,      GAINED
    REAL(fp)               :: WETLOSS,     MASS_WASH, MASS_NOWASH
    REAL(fp)               :: QDOWN,       DT,        F_WASHOUT
    REAL(fp)               :: K_RAIN,      WASHFRAC,  WET_Hg2
    REAL(fp)               :: WET_HgP,     MB,        QB
    REAL(fp)               :: QB_NUM,      DELP_DRY_NUM

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

    ! Arrays
    REAL(fp)               :: BMASS    (State_Grid%NZ)
    REAL(fp)               :: PDOWN    (State_Grid%NZ)

    ! Pointers
    REAL(fp),      POINTER :: BXHEIGHT (:        )
    REAL(fp),      POINTER :: CMFMC    (:        )
    REAL(fp),      POINTER :: DQRCU    (:        )
    REAL(fp),      POINTER :: DTRAIN   (:        )
    REAL(fp),      POINTER :: PFICU    (:        )
    REAL(fp),      POINTER :: PFLCU    (:        )
    REAL(fp),      POINTER :: REEVAPCN (:        )
    REAL(fp),      POINTER :: DELP_DRY (:        )
    REAL(fp),      POINTER :: T        (:        )
    REAL(fp),      POINTER :: H2O2s    (:        )
    REAL(fp),      POINTER :: SO2s     (:        )
    REAL(fp),      POINTER :: Q        (:,:      )
    TYPE(Species), POINTER :: SpcInfo

    !========================================================================
    ! (0)  I n i t i a l i z a t i o n
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Do_Cloud_Convection (in convection_mod.F90)'

    ! Point to columns of derived-type object fields
    BXHEIGHT => State_Met%BXHEIGHT(I,J,:        ) ! Box height [m]
    CMFMC    => State_Met%CMFMC   (I,J,2:State_Grid%NZ+1) ! Cloud mass flux
                                                          ! [kg/m2/s]
    DQRCU    => State_Met%DQRCU   (I,J,:        ) ! Precip production rate:
    DTRAIN   => State_Met%DTRAIN  (I,J,:        ) ! Detrainment flux [kg/m2/s]
    REEVAPCN => State_Met%REEVAPCN(I,J,:        ) ! Evap of precip'ing conv.
    DELP_DRY => State_Met%DELP_DRY(I,J,:        ) ! Edge dry P diff [hPa]
    T        => State_Met%T       (I,J,:        ) ! Air temperature [K]
    H2O2s    => State_Chm%H2O2AfterChem(I,J,:   ) ! H2O2s from sulfate_mod
    SO2s     => State_Chm%SO2AfterChem (I,J,:   ) ! SO2s from sulfate_mod
    Q        => State_Chm%Species (I,J,:,:      ) ! Chemical species
                                                  ! [mol/mol dry air]
    SpcInfo  => NULL()                            ! Species database entry

    ! PFICU and PFLCU are on level edges
    PFICU    => State_Met%PFICU   (I,J,2:State_Grid%NZ+1) ! Dwnwd flx of conv
                                                          !  ice precip
                                                          !  [kg/m2/s]
    PFLCU    => State_Met%PFLCU   (I,J,2:State_Grid%NZ+1) ! Dwnwd flux of conv
                                                          !  liquid precip
                                                          !  [kg/m2/s]

    ! # of levels and # of species
    NLAY     = State_Grid%NZ
    NC       = State_Chm%nAdvect

    ! Top level for convection
    KTOP     = NLAY - 1

    ! Convection timestep [s]
    NDT      = TS_DYN

    ! Internal time step for convective mixing is 300 sec.
    ! Doug Rotman (LLNL) says that 450 sec works just as well.
    NS       = NDT / 300                ! Num internal timesteps (int)
    NS       = MAX( NS, 1 )             ! Set lower bound to 1
    DNS      = DBLE( NS )               ! Num internal timesteps (real)
    SDT      = DBLE( NDT ) / DBLE( NS ) ! seconds in internal timestep

    !-----------------------------------------------------------------
    ! Determine location of the cloud base, which is the level where
    ! we start to have non-zero convective precipitation formation
    !-----------------------------------------------------------------

    ! Minimum value of cloud base is the surface level
    CLDBASE = 1

    ! Find the cloud base
    DO K = 1, NLAY
       IF ( DQRCU(K) > 0e+0_fp ) THEN
          CLDBASE = K
          EXIT
       ENDIF
    ENDDO

    !-----------------------------------------------------------------
    ! Compute PDOWN and BMASS
    !-----------------------------------------------------------------
    DO K = 1, NLAY

       ! PDOWN is the convective precipitation leaving each
       ! box [cm3 H2O/cm2 air/s]. This accounts for the
       ! contribution from both liquid & ice precip.
       ! PFLCU and PFICU are converted from kg/m2/s to m3/m2/s
       ! using water and ice densities, respectively.
       ! m3/m2/s becomes cm3/cm2/s using a factor of 100.
       PDOWN(K) = ( ( PFLCU(K) / 1000e+0_fp ) &
                +   ( PFICU(K) /  917e+0_fp ) ) * 100e+0_fp

       ! BMASS is the dry air mass per unit area for the grid box
       ! bounded by level K and K+1 [kg/m2]
       ! BMASS is equivalent to deltaP (dry) * 100 / g
       ! This is done to keep BMASS in the same units as CMFMC * SDT
       BMASS(K) = DELP_DRY(K) * G0_100

    ENDDO

    !-----------------------------------------------------------------
    ! Compute MB, the mass per unit area of dry air below the cloud
    ! base [kg/m2]. Calculate MB by looping over levels below the
    ! cloud base.
    !-----------------------------------------------------------------
    MB = 0e+0_fp
    DO K = 1, CLDBASE-1
       MB = MB + BMASS(K)
    ENDDO

    ! Is this a Hg simulation?
    IS_Hg = Input_Opt%ITS_A_MERCURY_SIM

    !========================================================================
    ! (1)  A d v e c t e d   S p e c i e s   L o o p
    !========================================================================

    ! Loop over only the advected species
    DO NA = 1, NC

       ! Get the species ID (modelID) from the advected species ID
       IC       =  State_Chm%Map_Advect(NA)

       ! Look up the corresponding entry in the species database
       SpcInfo  => State_Chm%SpcData(IC)%Info

       ! Also get the corresponding wetdep ID
       NW       =  SpcInfo%WetDepId

       ! Zero the DIAG14 diagnostic array
       DIAG14(:,NA) = 0.0_fp

       ! Zero the DIAG38 diagnostic array
       IF ( NW > 0 ) THEN
          DIAG38(:,NW) = 0.0_fp
       ENDIF

       !=====================================================================
       ! (2)  I n t e r n a l   T i m e   S t e p   L o o p
       !=====================================================================
       DO ISTEP = 1, NS

          ! Initialize
          QC     = 0e+0_fp    ! [kg species/kg dry air]
          T0_SUM = 0e+0_fp    ! [kg species/m2/timestep]

          !----------------------------------------------------------
          ! B e l o w   C l o u d   B a s e   (K < CLDBASE)
          !
          ! QB is the "weighted avg" mixing ratio below the cloud
          ! base [kg/kg dry air].
          ! QC is the mixing ratio of the air that moved in cumulus
          ! transport up to the next level [kg/kg dry air].
          ! MB is the dry mass of air below the cloud base per
          ! unit area [kg/m2] (see calculation before loop).
          !-----------------------------------------------------------

          ! We need to make this a nested IF statement so that we don't
          ! get an out-of-bounds error when CLDBASE=1 (bmy, 11/18/10)
          IF ( CLDBASE > 1 ) THEN

             IF ( CMFMC(CLDBASE-1) > TINYNUM ) THEN

                !-----------------------------------------------------
                ! %%% Non-negligible Cloud mass flux %%%
                !-----------------------------------------------------

                ! Calculate QB_NUM, the numerator for QB. QB is the
                ! weighted average mixing ratio below the cloud base.
                ! QB_NUM is equal to the grid box species concentrations 
                ! [kg/kg dry air] weighted by the adjacent level pressure 
                ! differences and summed over all levels up to just
                ! below the cloud base (ewl, 6/22/15)
                QB_NUM  = 0e+0_fp
                DELP_DRY_NUM = 0e+0_fp

                DO K  = 1, CLDBASE-1
                   QB_NUM = QB_NUM + Q(K,IC) * DELP_DRY(K)
                   DELP_DRY_NUM = DELP_DRY_NUM + DELP_DRY(K)
                ENDDO

                ! Compute QB, the weighted avg mixing ratio below
                ! the cloud base [kg/kg dry air]
                QB = QB_NUM / DELP_DRY_NUM

                ! Compute QC, the mixing ratio of the air that moved
                ! in cumulus transport up to the next level [kg/kg]
                !
                !        Dry mass of species below cloud base  +
                !        Subsidence into cloud base from above
                ! QC =  --------------------------------------------
                !            Dry air mass below cloud base
                !
                QC = ( MB*QB + CMFMC(CLDBASE-1) * Q(CLDBASE,IC) * SDT  ) / &
                     ( MB    + CMFMC(CLDBASE-1) * SDT  )

                ! Copy QC to all levels of the species array Q
                ! that are below the cloud base level [kg/kg]
                Q(1:CLDBASE-1,IC) = QC

             ELSE

                !-----------------------------------------------------
                ! %%% Negligible cloud mass flux %%%
                !-----------------------------------------------------

                ! When CMFMC is negligible, then set QC to the species
                ! concentration at the cloud base level [kg/kg]
                QC = Q(CLDBASE,IC)

             ENDIF

          ELSE

             !-----------------------------------------------------
             ! If the cloud base happens at level 1, then just
             ! set QC to the species concentration at the surface
             ! level [kg/kg]
             !-----------------------------------------------------
             QC = Q(CLDBASE,IC)

          ENDIF

          !==================================================================
          ! (3)  A b o v e   C l o u d   B a s e
          !==================================================================
          DO K = CLDBASE, KTOP

             ! Initialize
             ALPHA   = 0e+0_fp
             ALPHA2  = 0e+0_fp
             CMOUT   = 0e+0_fp
             ENTRN   = 0e+0_fp
             QC_PRES = 0e+0_fp
             QC_SCAV = 0e+0_fp

             ! CMFMC_BELOW is the air mass [kg/m2/s] coming into the
             ! grid box (K) from the box immediately below (K-1).
             IF ( K == 1 ) THEN
                CMFMC_BELOW = 0e+0_fp
             ELSE
                CMFMC_BELOW = CMFMC(K-1)
             ENDIF

             ! If we have a nonzero air mass flux coming from
             ! grid box (K-1) into (K) ...
             IF ( CMFMC_BELOW > TINYNUM ) THEN

                !------------------------------------------------------------
                ! (3.1)  M a s s   B a l a n c e   i n   C l o u d
                !
                ! F(K,NA) = fraction of species IC in level K that is
                !           available for wet-scavenging by cloud updrafts.  
                !
                ! If ENTRN > 0 then compute the new value of QC:
                !
                !      species mass from below      (i.e. level K-1) +
                !      species mass from this level (i.e. level K)
                !  = -----------------------------------------------------
                !             dry mass coming into cloud
                !
                ! Otherwise, preserve the previous value of QC.  This will
                ! ensure that TERM1 - TERM2 is not a negative quantity (see
                ! below).
                !
                ! Entrainment must be >= 0 (since we cannot have a negative 
                ! flux of air into the cloud).  This condition is strong 
                ! enough to ensure that CMOUT > 0 and will prevent floating-
                ! point exception.
                !------------------------------------------------------------

                ! Air mass flowing out of cloud at grid box (K) [kg/m2/s]
                CMOUT   = CMFMC(K) + DTRAIN(K)

                ! Air mass flowing into cloud at grid box (K) [kg/m2/s]
                ENTRN   = CMOUT - CMFMC_BELOW

                ! Amount of QC preserved against scavenging [kg/kg]
                !QC_PRES = QC * ( 1e+0_fp - F(K,IC) )
                QC_PRES = QC * ( 1e+0_fp - F(K,NA) )

                ! Amount of QC lost to scavenging [kg/kg]
                ! QC_SCAV = 0 for non-soluble species
                !QC_SCAV = QC * F(K,IC)
                QC_SCAV = QC * F(K,NA)

                ! - - - - - - - - FOR SOLUBLE SPECIES ONLY - - - - - - - - - 
                IF ( QC_SCAV > 0e+0_fp ) THEN

                   ! The fraction ALPHA is the fraction of raindrops that
                   ! will re-evaporate soluble species while falling from
                   ! grid box K+1 down to grid box K.  Avoid div-by-zero.

                   ! Initialize
                   ALPHA = 0e+0_fp

                   IF ( PDOWN(K+1)  > TINYNUM ) THEN

                      ! %%%% CASE 1 %%%%
                      ! Partial re-evaporation. Less precip is leaving
                      ! the grid box then entered from above.
                      IF ( PDOWN(K+1) > PDOWN(K) .AND. &
                           PDOWN(K)   > TINYNUM        ) THEN

                         ! Define ALPHA, the fraction of raindrops that
                         ! re-evaporate when falling from grid box
                         ! (I,J,L+1) to (I,J,L):
                         ! NOTE:
                         !   REEVAPCN is in units of [kg/kg/s]
                         !   Now use BMASS [kg/m2] instead of AD/area to
                         !   remove area dependency
                         !   PDOWN is in units of [cm3/cm2/s]
                         !   Factor of 10 in denom for unit conversion
                         !     1000 kg/m3 * 0.01 m/cm = 10 kg/m2/cm
                         ALPHA = REEVAPCN(K) * BMASS(K) &
                                 / ( PDOWN(K+1)  * 10e+0_fp )

                         ! Restrict ALPHA to be less than 1
                         ! (>1 is unphysical)  (hma, 24-Dec-2010)
                         IF ( ALPHA > 1e+0_fp ) THEN
                            ALPHA = 1e+0_fp
                         ENDIF

                         ! We assume that 1/2 of the soluble species w/in
                         ! the raindrops actually gets resuspended into
                         ! the atmosphere
                         ALPHA2   = ALPHA * 0.5e+0_fp

                      ENDIF

                      ! %%%% CASE 2 %%%%
                      ! Total re-evaporation. Precip entered from above,
                      ! but no precip is leaving grid box (ALPHA = 2 so
                      ! that  ALPHA2 = 1)
                      IF ( PDOWN(K) < TINYNUM ) THEN
                         ALPHA2 = 1e+0_fp
                      ENDIF

                   ENDIF

                   ! The resuspension takes 1/2 the amount of the scavenged
                   ! aerosol (QC_SCAV) and adds that back to QC_PRES
                   QC_PRES  = QC_PRES + ( ALPHA2 * QC_SCAV )

                   ! ... then we decrement QC_SCAV accordingly
                   QC_SCAV  = QC_SCAV * ( 1e+0_fp    - ALPHA2     )

                ENDIF
                !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                ! Update QC taking entrainment into account [kg/kg]
                ! Prevent div by zero condition
                IF ( ENTRN >= 0e+0_fp .and. CMOUT > 0e+0_fp ) THEN
                   QC   = ( CMFMC_BELOW * QC_PRES   + &
                          ENTRN       * Q(K,IC) ) / CMOUT
                ENDIF

                !------------------------------------------------------------
                ! (3.2)  M a s s   B a l a n c e   i n   L e v e l  ==> Q
                !
                ! Terminology:
                !
                !  C_k-1   = cloud air mass flux from level k-1 to level k
                !  C_k     = cloud air mass flux from level k   to level k+1
                !  QC_k-1  = mixing ratio of species INSIDE CLOUD at level k-1
                !  QC_k    = mixing ratio of species INSIDE CLOUD at level k
                !  Q_k     = mixing ratio of species in level k
                !  Q_k+1   = mixing ratio of species in level k+1
                !
                ! For convenience we denote:
                !
                !  QC_SCAV = Amount of species wet-scavenged in updrafts
                !          = QC_k-1 * F(k,IC)    [kg/kg]
                !
                !  QC_PRES = Amount of species preserved against
                !            wet-scavenging in updrafts [kg/kg]
                !          = QC_k-1 * ( 1 - F(k,IC) )
                !
                ! Where F(k,IC) is the fraction of species IC in level k
                ! that is available for wet-scavenging by cloud updrafts.
                ! F(k,IC) is computed by routine COMPUTE_UPDRAFT_FSOL
                ! and passed to this routine as an argument.
                !
                ! The cumulus transport above the cloud base is done as
                ! follows:
                !
                !                 ||///////////////////||
                !                 ||//// C L O U D ////||
                !                 ||                   ||
                !   k+1     ^     ||         ^         ||3)   C_k * Q_k+1
                !           |     ||         |         ||         |
                !   --------|-----++---------|---------++---------|--------
                !           |     ||         |         ||         |
                !   k      C_k    ||2)   C_k * QC_k    ||         V
                !                 ||                   ||
                !                 ||                   ||
                !           ^     ||         ^         ||4)   C_k-1 * Q_k
                !           |     ||         |         ||         |
                !   --------|-----++---------|---------++---------|--------
                !           |     ||         |         ||         |
                !   k-1   C_k-1   ||1) C_k-1 * QC_k-1  ||         V
                !                 ||         * (1 - F) ||
                !                 ||                   ||
                !                 ||//// C L O U D ////||
                !                 ||///////////////////||
                !
                ! There are 4 terms that contribute to mass flow in
                ! and out of level k:
                !
                ! 1) C_k-1 * QC_PRES = species convected from k-1 to k
                ! 2) C_k   * QC_k    = species convected from k   to k+1
                ! 3) C_k   * Q_k+1   = species subsiding from k+1 to k
                ! 4) C_k-1 * Q_k     = species subsiding from k   to k-1
                !
                ! Therefore the change in species concentration is given by
                !
                !    DELQ = (Term 1) - (Term 2) + (Term 3) - (Term 4)
                !
                ! and Q(K,IC) = Q(K,IC) + DELQ.
                !
                ! The term T0 is the amount of species that is scavenged
                ! out of the box.
                !
                ! Units of T0, T1, T2, T3, T4, and TSUM are
                ! [kg/m2/s * kg species / kg dry air]
                !------------------------------------------------------------
                T0      =  CMFMC_BELOW * QC_SCAV
                T1      =  CMFMC_BELOW * QC_PRES
                T2      = -CMFMC(K  )  * QC
                T3      =  CMFMC(K  )  * Q(K+1,IC)
                T4      = -CMFMC_BELOW * Q(K,  IC)

                TSUM    = T1 + T2 + T3 + T4

                DELQ    = ( SDT / BMASS(K) ) * TSUM    ! change in [kg/kg]

                ! If DELQ > Q then do not make Q negative!!!
                IF ( Q(K,IC) + DELQ < 0 ) THEN
                   DELQ = -Q(K,IC)
                ENDIF

                ! Increment the species array [kg/kg]
                Q(K,IC) = Q(K,IC) + DELQ

                ! Return if we encounter NInf
                IF ( .not. IT_IS_FINITE( Q(K,IC) ) ) THEN
                   WRITE( 6, 200 )
200                FORMAT( 'Infinity in DO_CLOUD_CONVECTION!' )
                   WRITE( 6, 255 ) K, IC, Q(K,IC)
                   RC = GC_FAILURE
                   RETURN
                ENDIF

                ! Return if we encounter NaN
                IF ( IT_IS_NAN( Q(K,IC) ) ) THEN
                   WRITE( 6, 250 )
                   WRITE( 6, 255 ) K, IC, Q(K,IC)
250                FORMAT( 'NaN encountered in DO_CLOUD_CONVECTION!' )
255                FORMAT( 'K, IC, Q(K,IC): ', 2i4, 1x, es13.6 )
                   RC = GC_FAILURE
                   RETURN
                ENDIF

                ! Pass T0_SUM in units of [kg species/m2/timestep].
                ! Converting kg dry air to kg species requires use
                ! of the molecular weight of air including moisture
                ! (ewl, 6/5/15)
                T0_SUM = T0_SUM + T0 * SDT

                !------------------------------------------------------------
                ! (3.3)  N D 1 4   D i a g n o s t i c
                !
                ! Archive upward mass flux due to wet convection
                ! [kg/sec] in the box (I,J), for the species IC going out 
                ! of the top of the layer K to the layer above (K+1)
                ! (bey, 11/10/99). We must divide by DNS, the # of internal 
                ! timesteps so that the sum represents the average loss rate 
                ! across all internal timesteps.
                !------------------------------------------------------------
                IF ( USE_DIAG14 ) THEN
                   DIAG14(K,NA) = DIAG14(K,NA) + ( -T2-T3 ) * AREA_M2 / DNS
                ENDIF

                !------------------------------------------------------------
                ! (3.4)  N D 3 8   D i a g n o s t i c
                !
                ! Archive the loss of soluble species to wet scavenging in
                ! cloud updrafts [kg/s].  We must divide by DNS, the # of
                ! internal timesteps so that the sum represents the
                ! average loss rate across all internal timesteps.
                !------------------------------------------------------------
                IF ( USE_DIAG38 .and. NW > 0 ) THEN
                   DIAG38(K,NW) = DIAG38(K,NW) + ( T0 * AREA_M2 / DNS )

                   ! check for infinity (added by hma, 20101117)
                   IF ( .not. IT_IS_FINITE( DIAG38(K,NW) ) ) THEN
                      WRITE( ErrMsg, 300 ) K
300                   FORMAT( 'DIAG38 is infinity (K,NW)= ', i6,' #1')
                      CALL GC_Error( ErrMsg, RC, ThisLoc )
                      RETURN
                   ENDIF
                ENDIF

             ELSE

                !------------------------------------------------------------
                ! (3.5)  N o   C l o u d   M a s s   F l u x   B e l o w
                !------------------------------------------------------------

                ! If there is no cloud mass flux coming from below, set
                ! QC to the species concentration at this level [kg/kg]
                QC = Q(K,IC)

                ! Bug fix for the cloud base layer, which is not necessarily
                ! in the boundary layer, and there could be
                ! "secondary convection" plumes - one in the PBL and another 
                ! one not.  NOTE: T2 and T3 are the same terms as described 
                ! in the above section.  (swu, 08/13/2007)
                IF ( CMFMC(K) > TINYNUM ) THEN

                   ! Species convected from K -> K+1
                   ! [kg/m2/s * kg species/kg dry air]
                   T2   = -CMFMC(K) * QC

                   ! Species subsiding from K+1 -> K [kg/m2/s]
                   ! [kg/m2/s * kg species/kg dry air]
                   T3   =  CMFMC(K) * Q(K+1,IC)

                   ! Change in species concentration [kg/kg]
                   DELQ = ( SDT / BMASS(K) ) * (T2 + T3)

                   ! If DELQ > Q then do not make Q negative!!!
                   IF ( Q(K,IC) + DELQ < 0.0e+0_fp ) THEN
                      DELQ = -Q(K,IC)
                   ENDIF

                   ! Add change in species to Q array [kg/kg]
                   Q(K,IC) = Q(K,IC) + DELQ

                ENDIF
             ENDIF
          ENDDO     ! End of loop over levels above cloud base

          !==================================================================
          ! (4)  B e l o w   C l o u d   B a s e
          !==================================================================
          DO K = CLDBASE-1, 1, -1

             ! Initialize
             QDOWN       = 0e+0_fp
             F_WASHOUT   = 0e+0_fp
             WASHFRAC    = 0e+0_fp
             ALPHA       = 0e+0_fp
             ALPHA2      = 0e+0_fp
             GAINED      = 0e+0_fp
             WETLOSS     = 0e+0_fp
             LOST        = 0e+0_fp
             MASS_WASH   = 0e+0_fp
             MASS_NOWASH = 0e+0_fp
             K_RAIN      = 0e+0_fp

             ! Check if...
             ! (1) there is precip coming into box (I,J,K) from (I,J,K+1)
             ! (2) there is re-evaporation happening in grid box (I,J,K)
             ! (3) there is species to re-evaporate
             IF ( PDOWN(K+1)  > 0 .and. &
                  REEVAPCN(K) > 0 .and. &
                  T0_SUM      > 0        ) THEN

                ! Compute F_WASHOUT, the fraction of grid box (I,J,L)
                ! experiencing washout. First, convert units of PDOWN, 
                ! the downward flux of precip leaving grid box (K+1)
                ! from [cm3 H20/cm2 area/s] to [cm3 H20/cm3 air/s]
                ! by dividing by box height in cm
                QDOWN = PDOWN(K+1) / ( BXHEIGHT(K+1) * 100e+0_fp  )

                ! Compute K_RAIN and F_WASHOUT based on the flux of precip 
                ! leaving grid box (K+1).
#ifdef LUO_WETDEP
                ! Luo et al scheme: Use COND_WATER_CONTENT = 2e-6 [cm3/cm3]
                K_RAIN   = LS_K_RAIN( QDOWN, 2.0e-6_fp )
                F_WASHOUT= CONV_F_PRIME( QDOWN, K_RAIN, SDT )
#else
                ! Default scheme: Use COND_WATER_CONTENT = 1e-6 [cm3/cm3]
                ! (which was recommended by Qiaoqiao Wang et al [2014])
                K_RAIN   = LS_K_RAIN(  QDOWN,         1.0e-6_fp )
                F_WASHOUT= LS_F_PRIME( QDOWN, K_RAIN, 1.0e-6_fp )
#endif

                ! Call WASHOUT to compute the fraction of species lost
                ! to washout in grid box (I,J,K)
                CALL WASHOUT( I,         J,                         &
                              K,         IC,         BXHEIGHT(K),   &
                              T(K),      QDOWN,      SDT,           &
                              F_WASHOUT, H2O2s(K),   SO2s(K),       &
                              WASHFRAC,  AER,        Input_Opt,     &
                              State_Chm, State_Grid, State_Met,  RC )

                ! Trap potential errors
                IF ( RC /= GC_SUCCESS ) THEN
                   ErrMsg = 'Error encountered in "Washout"!'
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF

                ! Check if the species is an aerosol or not
                IF ( AER ) THEN

                   !---------------------------------------------------------
                   ! Washout of aerosol species
                   ! This is modeled as a kinetic process
                   !---------------------------------------------------------

                   ! Define ALPHA, the fraction of raindrops that
                   ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                   ! NOTE:
                   !   REEVAPCN is in units of [kg/kg/s]
                   !   Now use BMASS [kg/m2] instead of AD/area to
                   !   remove area dependency
                   !   PDOWN is in units of [cm3/cm2/s]
                   !   Factor of 10 in denom for unit conversion
                   !     1000 kg/m3 * 0.01 m/cm = 10 kg/m2/cm

                   ! %%%% CASE 1 %%%%
                   ! Partial re-evaporation. Less precip is leaving
                   ! the grid box then entered from above (V. Shah, 9/14/15)
                   IF ( PDOWN(K+1) > PDOWN(K) .AND. &
                        PDOWN(K)   > TINYNUM        ) THEN

                      ! Define ALPHA, the fraction of raindrops that
                      ! re-evaporate when falling from grid box
                      ! (I,J,L+1) to (I,J,L)
                      ALPHA = REEVAPCN(K) * BMASS(K) &
                              / ( PDOWN(K+1) * 10e+0_fp  )

                      ! For safety
                      ALPHA = MIN( ALPHA, 1e+0_fp )

                      ! ALPHA2 is the fraction of the rained-out aerosols
                      ! that gets resuspended in grid box (I,J,L)
                      ALPHA2  = 0.5e+0_fp * ALPHA

                   ENDIF

                   ! %%%% CASE 2 %%%%
                   ! Total re-evaporation. Precip entered from above,
                   ! but no precip is leaving grid box (ALPHA = 2 so
                   ! that  ALPHA2 = 1) (V. Shah, 9/14/15)
                   IF ( PDOWN(K) < TINYNUM ) THEN
                      ALPHA2 = 1e+0_fp
                   ENDIF

                   ! GAINED is the rained out aerosol coming down from
                   ! grid box (I,J,L+1) that will evaporate and re-enter
                   ! the atmosphere in the gas phase in grid box (I,J,L)
                   ! [kg species/m2/timestep]
                   GAINED  = T0_SUM * ALPHA2

                   ! Amount of aerosol lost to washout in grid box [kg/m2]
                   ! (V. Shah, 9/14/15)
                   WETLOSS = ( Q(K,IC) * BMASS(K) + GAINED ) * &
                             WASHFRAC - GAINED

                   ! LOST is the rained out aerosol coming down from
                   ! grid box (I,J,L+1) that will remain in the liquid
                   ! phase in grid box (I,J,L) and will NOT re-evaporate
                   ! [kg/m2/timestep]
                   LOST    = T0_SUM - GAINED

                   ! Update species concentration (V. Shah, mps, 5/20/15)
                   ! [kg/kg]
                   Q(K,IC) = Q(K,IC) - WETLOSS / BMASS(K)

                   ! Update T0_SUM, the total amount of scavenged
                   ! species that will be passed to the grid box below
                   ! [kg/m2/timestep]
                   T0_SUM = T0_SUM + WETLOSS

                ELSE

                   !---------------------------------------------------------
                   ! Washout of non-aerosol species
                   ! This is modeled as an equilibrium process
                   !---------------------------------------------------------

                   ! MASS_NOWASH is the amount of non-aerosol species in
                   ! grid box (I,J,L) that is NOT available for washout.
                   ! Calculate in units of [kg/kg]
                   MASS_NOWASH = ( 1e+0_fp - F_WASHOUT ) * Q(K,IC)

                   ! MASS_WASH is the total amount of non-aerosol species
                   ! that is available for washout in grid box (I,J,L).
                   ! It consists of the mass in the precipitating
                   ! part of box (I,J,L), plus the previously rained-out
                   ! species coming down from grid box (I,J,L+1).
                   ! (Eq. 15, Jacob et al, 2000)
                   ! Units are [kg species/m2/timestep]
                   MASS_WASH = ( F_WASHOUT * Q(K,IC) ) * BMASS(K) + T0_SUM

                   ! WETLOSS is the amount of species mass in
                   ! grid box (I,J,L) that is lost to washout.
                   ! (Eq. 16, Jacob et al, 2000)
                   ! [kg species/m2/timestep]
                   WETLOSS     = MASS_WASH * WASHFRAC - T0_SUM

                   ! The species left in grid box (I,J,L) is what was
                   ! originally in the non-precipitating fraction
                   ! of the box, plus MASS_WASH, less WETLOSS.
                   ! [kg/kg]
                   Q(K,IC) = Q(K,IC) - WETLOSS / BMASS(K)

                   ! Update T0_SUM, the total scavenged species
                   ! that will be passed to the grid box below
                   ! [kg species/m2/timestep]
                   T0_SUM      = T0_SUM + WETLOSS

                ENDIF

                !------------------------------------------------------------
                ! N D 1 4   D i a g n o s t i c
                !
                ! Archive upward mass flux due to wet convection.
                ! [kg/sec] in the box (I,J), for the species IC going
                ! out of the top of the layer K to the layer above (K+1)  
                ! (bey, 11/10/99). We must divide by DNS, the # of internal 
                ! timesteps so that the sum represents the average loss 
                ! rate across all internal timesteps.
                !------------------------------------------------------------
                IF ( USE_DIAG14 ) THEN
                   DIAG14(K,NA) = DIAG14(K,NA) + ( -T2-T3 ) * AREA_M2 / DNS
                ENDIF

                !------------------------------------------------------------
                !  N D 3 8   D i a g n o s t i c
                !
                ! Archive the loss of soluble species to wet scavenging in
                ! cloud updrafts [kg/s].  We must divide by NDT, the # of
                ! seconds in the convective timestep, equal to DNS * SDT,
                ! in order to make diag38 represent the average loss rate
                ! across all internal timesteps. Note that the units of
                ! WETLOSS are [kg/m2/timestep].
                !------------------------------------------------------------
                !%%% NOTE: SHOULD TEST FOR NW > 0 BUT IF WE DO THAT WE
                !%%% NO LONGER GET IDENTICAL RESULTS WITH THE REF CODE.
                !%%% LOOK INTO THIS LATER.  (bmy, 7/7/16)
                IF ( USE_DIAG38 .and. F(K,NA) > 0.0_fp ) THEN
                   DIAG38(K,NW) = DIAG38(K,NW) + ( WETLOSS * AREA_M2 / NDT )
                ENDIF

                ! check for infinity (added by hma, 20101117)
                IF ( .not. IT_IS_FINITE( DIAG38(K,NW) ) ) THEN
                   WRITE( ErrMsg, 310 ) K, NW
310                FORMAT( 'DIAG38 is infinity (K,NW)= ', 2i6, ' #3' )
                   CALL GC_Error( ErrMsg, RC, ThisLoc )
                   RETURN
                ENDIF
             ENDIF
          ENDDO     ! End of loop over levels below cloud base

          !==================================================================
          ! (5)  M e r c u r y   O c e a n   M o d e l   A r c h i v a l
          !
          ! Pass the amount of Hg2 and HgP lost in wet  scavenging [kg] to 
          ! "ocean_mercury_mod.f" via ADD_Hg2_WET and ADD_HgP_WET.   We must 
          ! also divide  by DNS, the # of internal timesteps.
          ! (sas, bmy, eck, eds, 1/19/05, 1/6/06, 7/30/08)
          !
          ! NOTE: Reorder
          !==================================================================
          IF ( Is_Hg ) THEN

             !--------------------------------------
             ! Hg2
             !--------------------------------------
             IF ( SpcInfo%IS_Hg2 ) THEN

                ! Wet scavenged Hg(II) in [kg]
                WET_Hg2 = ( T0_SUM * AREA_M2 )

                ! Category # for this Hg2 species
                Hg_Cat  = SpcInfo%Hg_Cat

                ! Pass to "ocean_mercury_mod.f"
                CALL ADD_Hg2_WD      ( I, J, Hg_Cat, WET_Hg2  )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_Hg2, &
                                       State_Met, State_Chm, State_Diag )
             ENDIF

             !--------------------------------------
             ! HgP
             !--------------------------------------
             IF ( SpcInfo%Is_HgP ) THEN

                ! Wet scavenged Hg(P) in [kg]
                WET_HgP = ( T0_SUM * AREA_M2 )

                ! Category # for this Hg2 species
                Hg_Cat  = SpcInfo%Hg_Cat

                ! Pass to "ocean_mercury_mod.f"
                CALL ADD_HgP_WD      ( I, J, Hg_Cat, WET_HgP  )
                CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat, WET_HgP, &
                                       State_Met, State_Chm, State_Diag )
             ENDIF
          ENDIF
       ENDDO               ! End internal timestep loop

       ! Free pointer
       SpcInfo => NULL()
    ENDDO                  ! End loop over advected species

    !================================================================
    ! Succesful return!
    !================================================================

    ! Nullify pointers
    NULLIFY( BXHEIGHT )
    NULLIFY( CMFMC    )
    NULLIFY( DQRCU    )
    NULLIFY( DTRAIN   )
    NULLIFY( PFICU    )
    NULLIFY( PFLCU    )
    NULLIFY( REEVAPCN )
    NULLIFY( DELP_DRY )
    NULLIFY( T        )
    NULLIFY( H2O2s    )
    NULLIFY( SO2s     )
    NULLIFY( Q        )

    ! Set error code to success
    RC                      = GC_SUCCESS

  END SUBROUTINE DO_CLOUD_CONVECTION
!EOC
END MODULE CONVECTION_MOD
