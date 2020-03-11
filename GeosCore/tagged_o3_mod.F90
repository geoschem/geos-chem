!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tagged_o3_mod.F90
!
! !DESCRIPTION: Module TAGGED\_O3\_MOD contains variables and routines to
!  perform a tagged O3 simulation.  P(O3) and L(O3) rates need to be archived
!  from a full chemistry simulation before you can run w/ Tagged O3.
!\\
!\\
! !INTERFACE:
!
MODULE TAGGED_O3_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% If you want to the EXTENDED SIMULATION with all 13 tagged O3 species,
!%%% then uncomment this #ifdef statement. (bmy, 4/11/14)
!#define USE_ALL_TAGO3_SPECIES 1
!%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEM_TAGGED_O3
  PUBLIC  :: CLEANUP_TAGGED_O3
  PUBLIC  :: INIT_TAGGED_O3
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GET_REGIONAL_PO3
!
! !REMARKS:
!  THE SIMPLE TAGGED O3 SIMULATION (default setting) HAS THESE ADVECTED SPECIES:
!  ----------------------------------------------------------------------------
!  (1 ) O3      : Total O3
!  (2 ) O3Strt  : Stratospheric O3
!                                                                             .
!  THE EXTENDED TAGGED O3 SIMULATION HAS THESE ADVECTED SPECIES:
!  ----------------------------------------------------------------------------
!  (1 ) O3      : Total O3
!  (2 ) O3Strt  : O3 from the Stratosphere      (tropopause - atm top   )
!  (3 ) O3Ut    : O3 produced in Upper Trop     (350 hPa    - tropopause)
!  (4 ) O3Mt    : O3 produced in Middle Trop    (PBL top    - 350 hPa   )
!  (5 ) O3Row   : O3 produced in Rest of World  (surface    - PBL top   )
!  (6 ) O3PcBl  : O3 produced in Pacific BL     (surface    - PBL top   )
!  (7 ) O3NaBl  : O3 produced in N. American BL (surface    - PBL top   )
!  (8 ) O3AtBl  : O3 produced in Atlantic BL    (surface    - PBL top   )
!  (9 ) O3EuBl  : O3 produced in European BL    (surface    - PBL top   )
!  (10) O3AfBl  : O3 produced in N. African BL  (surface    - PBL top   )
!  (11) O3AsBl  : O3 produced in Asian          (surface    - PBL top   )
!  (12) O3Init  : O3 initial conditions         (all levels             )
!  (13) O3USA   : O3 produced over the USA      (all levels             )
!                                                                             .
!  NOTES:
!  ----------------------------------------------------------------------------
!  (1) The stratospheric O3 species must be species #2.  This is due to how
!       the Linoz stratospheric O3 chemistry scheme is written.  We have
!       accordingly reorganized the species numbers below.
!                                                                             .
!  (3) If using all of the tagged O3 species, note that that the sum of
!       production of species #2 - #12 sums together to equal the production
!       of species #1.  In other words:
!                                                                             .
!         PP(I,J,L,1) will equal SUM( PP(I,J,L,2:12) )
!                                                                             .
!       ALSO NOTE: The O3USA species is defined to include all of the
!       O3 produced over the US (at all levels).  Therefore, it should
!       be treated separately from species 2 - 12.
!                                                                             .
!  (4) When starting a long tagged O3 simulation, we recommend that you use
!      a restart file where all species concentrations are set to zero.
!      Then spin up for as many years as it takes to get into steady-state.
!      This will ensure that the sum of tagged O3 species (#2  - #12) will
!      equal the total O3 species (#1).
!
! !REVISION HISTORY:
!  20 Aug 2003 - A. Fiore    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! These are pointers to fields in the HEMCO data structure.
  ! Declare these with REAL(f4), aka REAL*4. (bmy, 3/4/15)
  REAL(f4), POINTER            :: P24H(:,:,:) => NULL() ! O3 production rate
  REAL(f4), POINTER            :: L24H(:,:,:) => NULL() ! O3 loss rate

  ! Emission timestep (will be imported from HEMCO)
  REAL(fp)                     :: TS_EMIS

  ! Species ID flags
  INTEGER                      :: id_O3Strat
!
! !DEFINED PARAMETERS:
!
  ! To convert m3 to cm3
  REAL(fp), PARAMETER          :: CM3PERM3 = 1.0e+6_fp

#if defined( USE_ALL_TAGO3_SPECIES )
  !-----------------------------------------------------------------
  ! EXTENDED SIMULATION : Total, strat, and regional O3 species
  !-----------------------------------------------------------------
  INTEGER,  PARAMETER, PRIVATE :: N_TAGGED = 13   ! # of species
  INTEGER,  PARAMETER, PRIVATE :: N_STRAT  = 2    ! Stratospheric O3
#else
  !-----------------------------------------------------------------
  ! SIMPLE SIMULATION: Total and strat O3 species only
  !
  ! %%% THIS IS THE DEFAULT %%%
  !-----------------------------------------------------------------
  INTEGER,  PARAMETER, PRIVATE :: N_TAGGED = 2    ! # of species
  INTEGER,  PARAMETER, PRIVATE :: N_STRAT  = 2    ! Stratospheric O3
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_regional_po3
!
! !DESCRIPTION: Subroutine GET\_REGIONAL\_PO3 returns the P(O3) for each of
!  the tagged O3 species. Tagged O3 species are defined by both geographic
!  location and altitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_REGIONAL_PO3( I, J, L, PP, State_Grid, State_Met )
!
! !USES:
!
    USE PhysConstants            ! SCALE_HEIGHT
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I, J, L     ! Grid box indices
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    ! Array containing P(O3) for each tagged species
    REAL(fp),  INTENT(OUT) :: PP(State_Grid%NX,State_Grid%NY,State_grid%NZ,&
                                 N_TAGGED)
!
! !REVISION HISTORY:
!  19 Aug 2003 - A. Fiore - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL  :: ITS_IN_TROP, ITS_IN_PBL, ITS_IN_MT
    LOGICAL  :: ITS_IN_UT,   ITS_IN_NH,  ITS_IN_ATL
    LOGICAL  :: ITS_IN_PAC,  ITS_IN_AS,  ITS_IN_EUR
    LOGICAL  :: ITS_IN_NAM,  ITS_IN_NAF, ITS_IN_USA
    INTEGER  :: PBLTOP,      MTTOP
    REAL(fp) :: PPROD,       X,          Y
    REAL(fp) :: BOXVL

    !=================================================================
    ! Initialization
    !=================================================================

    ! Initialize
    PP(I,J,L,:) = 0e+0_fp

    ! IS TROP is TRUE if we are in the troposphere
    ITS_IN_TROP = State_Met%InTroposphere(I,J,L)

    ! Skip stratospheric boxes
    IF ( .not. ITS_IN_TROP ) RETURN

    ! Longitude and latitude [degrees]
    X = State_Grid%XMid(I,J)
    Y = State_Grid%YMid(I,J)

    ! PBLTOP is the model level at ~ 750 hPa
    ! MTTOP  is the model level at ~ 350 hPa
    PBLTOP = 16
    MTTOP  = 27

#if defined( USE_ALL_TAGO3_SPECIES )

    !=================================================================
    ! For the simulation with all tagged O3 species: find the
    ! regions corresponding to this particular longitude & latitude
    !=================================================================

    ! Define flags for various geographic & altitude regions
    ITS_IN_PBL = ( L <= PBLTOP                                       )
    ITS_IN_MT  = ( L >  PBLTOP .and. L <= MTTOP                      )
    ITS_IN_UT  = ( L >  MTTOP  .and. ITS_IN_TROP                     )

    ITS_IN_NH  = ( Y >=   0.0                                        )
    ITS_IN_EUR = ( Y >=  36.0 .and. ( X >  -15.0 .and. X <=   55.0 ) )
    ITS_IN_NAM = ( Y >=  15.0 .and. ( X > -127.5 .and. X <=  -65.0 ) )
    ITS_IN_AS  = ( Y >= -10.0 .and. ( X >   55.0 .and. X <=  145.0 ) )
    ITS_IN_ATL = ( ITS_IN_NH  .and. ( X >  -65.0 .and. X <=  -15.0 ) )
    ITS_IN_PAC = ( ITS_IN_NH  .and. ( X >  145.0  .or. X <= -127.5 ) )

    ITS_IN_NAF = ( ( X >= -15.0 .and. X <=  55.0 ) .and. &
                   ( Y >=   0.0 .and. Y <   36.0 ) )

    ITS_IN_USA = ( ( X > -127.5 .and. X <= -65.0 ) .and. &
                   ( Y >   22.0 .and. Y <=  50.0 ) )

#endif

    !=================================================================
    ! Assign P(O3) to tagged species by geographic/altitude regions
    !
    ! NOTE: If using all of the tagged O3 species, note that the
    ! sum of production of species #2 - #12 sums together to equal
    ! the production of species #1.  In other words:
    !
    !    PP(I,J,L,1) will equal SUM( PP(I,J,L,2:12) )
    !
    ! ALSO NOTE: The O3USA species is defined to include all of the
    ! O3 produced over the US (at all levels).  Therefore, it should
    ! be treated separately from species #2 - #12.
    !=================================================================

    ! Grid box volume [cm3]
    BOXVL = State_Met%AIRVOL(I,J,L) !* 1d6

    ! P(O3) [kg]
    ! P24H is in kg/m3 per emission time step (ckeller, 9/17/2014).
    PPROD = P24H(I,J,L) * BOXVL * ( GET_TS_CHEM()/TS_EMIS )

    !-----------------------
    ! #1: Total P(O3)
    !-----------------------
    PP(I,J,L,1) = PPROD

#if defined( USE_ALL_TAGO3_SPECIES )

    !-----------------------
    ! #2: P(O3) in UT
    !-----------------------
    IF ( ITS_IN_UT ) THEN
       PP(I,J,L,3) = PPROD

    !-----------------------
    ! #3: P(O3) in MT
    !-----------------------
    ELSE IF ( ITS_IN_MT ) THEN
       PP(I,J,L,4) = PPROD

    !-----------------------
    ! #5: P(O3) in Pac BL
    !-----------------------
    ELSE IF ( ITS_IN_PAC .and. ITS_IN_PBL ) THEN
       PP(I,J,L,6) = PPROD

    !-----------------------
    ! #6: P(O3) in NAm BL
    !-----------------------
    ELSE IF ( ITS_IN_NAM .and. ITS_IN_PBL ) THEN
       PP(I,J,L,7) = PPROD

    !-----------------------
    ! #7: P(O3) in Atl BL
    !-----------------------
    ELSE IF ( ITS_IN_ATL .and. ITS_IN_PBL ) THEN
       PP(I,J,L,8) = PPROD

    !-----------------------
    ! #8: P(O3) in Eur BL
    !-----------------------
    ELSE IF ( ITS_IN_EUR .and. ITS_IN_PBL ) THEN
       PP(I,J,L,9) = PPROD

    !-----------------------
    ! #9: P(O3) in NAfr BL
    !-----------------------
    ELSE IF ( ITS_IN_NAF .and. ITS_IN_PBL ) THEN
       PP(I,J,L,10) = PPROD

    !-----------------------
    ! #10: P(O3) in Asia BL
    !-----------------------
    ELSE IF ( ITS_IN_AS .and. ITS_IN_PBL ) THEN
       PP(I,J,L,11) = PPROD

    !-----------------------
    ! #4: P(O3) in R.O.W
    !-----------------------
    ELSE
       PP(I,J,L,5) = PPROD

    ENDIF

    !-------------------------
    ! #13: P(O3) in USA
    !-------------------------
    IF ( ITS_IN_USA ) THEN
       PP(I,J,L,13) = PPROD
    ENDIF

#endif

  END SUBROUTINE GET_REGIONAL_PO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_tagged_o3
!
! !DESCRIPTION: Subroutine CHEM\_TAGGED\_O3 performs chemistry for several
!  O3 species which are tagged by geographic and altitude regions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_TAGGED_O3( Input_Opt, State_Chm, State_Diag, State_Grid, &
                             State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM

    IMPLICIT NONE
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
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMAKRS:
!  Dry deposition is now applied in mixing_mod.F90.  We have the application
!  of Ox dry deposition from this routine, as well as the archival of the
!  ND44 drydep diagnostic. (bmy, 6/15/15)
!
! !REVISION HISTORY:
!  20 Aug 2003 - R. Hudman   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE            :: FIRST = .TRUE.

    ! Scalars
    INTEGER                  :: I,     J,      L,  N,  NA, nAdvect
    REAL(fp)                 :: BOXVL, DTCHEM, DT, LL, PL

    ! Arrays
    REAL(fp)                 :: PP(State_Grid%NX,State_Grid%NY,State_grid%NZ, &
                                   N_TAGGED)

    ! Pointers
    REAL(fp),        POINTER :: Spc(:,:,:,:)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! CHEM_TAGGED_O3 begins here!
    !=================================================================

    ! Initialize
    RC        =  GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at CHEM_TAGGED_O3 (in GeosCore/tagged_o3_mod.F90)'

    ! Number of advected species
    nAdvect   =  State_Chm%nAdvect

    ! Pointers
    Spc       => State_Chm%Species   ! Points to chemical species [kg]

    ! Chemistry timestep [s]
    DTCHEM    =  GET_TS_CHEM()

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_Loss ) State_Diag%Loss = 0.0_f4
    IF ( State_Diag%Archive_Prod ) State_Diag%Prod = 0.0_f4

    !=================================================================
    ! Get production and loss frequencies from HEMCO. The target
    ! will be updated automatically by HEMCO.
    !
    ! Important: the file data are converted by HEMCO to HEMCO
    ! concentration units, e.g. the imported data is in kg/m3
    ! (production) and 1/m3 (loss), e.g. the original data in
    ! units of kg/m3/s multiplied by the emission timestep.
    !                                     (ckeller, 9/17/2014)
    !=================================================================
    IF ( FIRST ) THEN

       ! Make sure the HEMCO State object has been first allocated
       IF ( .NOT. ASSOCIATED(HcoState) ) THEN
          ErrMsg = 'The HcoState object is not allocated!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get pointer to O3 production
       CALL HCO_GetPtr( HcoState, 'O3_PROD', P24H, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to O3_PROD!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get pointer to O3 loss
       CALL HCO_GetPtr( HcoState, 'O3_LOSS', L24H, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to O3_LOSS!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get emissions timestep (seconds)
       TS_EMIS = HcoState%TS_EMIS

       ! Reset first-time variable
       ! so that this won't be executed again
       FIRST    = .FALSE.

    ENDIF

    ! DT is the ratio of the chemistry and emission time step.
    ! Use this value to convert from kg/m3 or 1/m3 per emission
    ! time step to kg/m3 or 1/m3 per chemistry time step.
    ! (ckeller, 9/17/2014).
    DT = DTCHEM / TS_EMIS

    !=================================================================
    ! Tagged O3 chemistry contains the following terms:
    !
    !   New O3 = Old O3 + ( P(O3,region) - L(O3) )
    !
    ! P(O3) and L(O3) are archived from a previous fullchem run using
    ! the ND20 diagnostic.  P(O3,region) is the P(O3) for a specific
    ! tagged O3 species, as computed by routine GET_REGIONAL_PO3.
    !
    ! Tagged O3 species are defined by both geographic location and
    ! altitude, as defined in GET_REGIONAL_PO3.  If you are running
    ! the
    !=================================================================

    ! Loop over the # of advected species
    DO NA = 1, nAdvect

       ! Advected species ID
       N = State_Chm%Map_Advect(NA)

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, BOXVL, LL, PL ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Grid box volume [m3]
          BOXVL = State_Met%AIRVOL(I,J,L) !* 1d6

          !===========================================================
          ! Get P(O3) and L(O3) for each tagged species in [kg/m3/s]
          !===========================================================

          ! P(O3) is a function of geographic & altitude location
          ! NOTE: We call this only when N==1 for optimal looping
          ! ALSO NOTE: PP is 4-D so it doesn't have to be PRIVATE.
          IF ( N == 1 ) THEN
             CALL GET_REGIONAL_PO3(I, J, L, PP, State_Grid, State_Met)
          ENDIF

          ! L(O3) is now in [1/m3] (ckeller, 9/17/2014)
          IF ( State_Met%InTroposphere(I,J,L) ) THEN
             LL = Spc(I,J,L,N) * L24H(I,J,L) * BOXVL * DT
          ELSE
             LL = 0.0e+0_fp
          ENDIF

          !===========================================================
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Chemical production and loss of tagged O3 species
          !===========================================================

          ! Production of tagged O3 species [kg/s]
          IF ( State_Diag%Archive_Prod ) THEN
             IF ( PP(I,J,L,N) > 0e+0_fp ) THEN
                State_Diag%Prod(I,J,L,N) = P24H(I,J,L) * BOXVL / TS_EMIS
             ENDIF
          ENDIF

          ! Loss of tagged O3 species [kg/s]
          IF ( State_Diag%Archive_Loss ) THEN
             State_Diag%Loss(I,J,L,N) = Spc(I,J,L,N) * L24H(I,J,L) &
                                        * BOXVL / TS_EMIS
          ENDIF

          !===========================================================
          ! Apply chemical P(O3) - L(O3) to each tagged species
          !===========================================================
          Spc(I,J,L,N) = Spc(I,J,L,N) + PP(I,J,L,N) - LL

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_TAGGED_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tagged_o3
!
! !DESCRIPTION: Subroutine INIT\_TAGGED\_O3 allocates and zeroes all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TAGGED_O3( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_TAGGED_O3 begins here
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Tagged_O3 (in module GeosCore/tagged_o3_mod.F90)'

    ! Define species ID flag
    id_O3Strat = Ind_('O3Strat')

    ! Add error check to make sure O3Strt is defined (bmy, 6/20/16)
    IF ( id_O3Strat <= 0 ) THEN
       ErrMsg = 'O3Strat is an undefined species!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Safety valve
    IF ( State_Chm%nAdvect > N_TAGGED ) THEN
       ErrMsg = 'State_Chm%nAdvect is too large for Tagged O3!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE INIT_TAGGED_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tagged_o3
!
! !DESCRIPTION:Subroutine CLEANUP\_TAGGED\_O3 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_TAGGED_O3()
!
! !REVISION HISTORY:
!  20 Aug 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ASSOCIATED( P24H ) ) P24H => NULL()
    IF ( ASSOCIATED( L24H ) ) L24H => NULL()

  END SUBROUTINE CLEANUP_TAGGED_O3
!EOC
END MODULE TAGGED_O3_MOD
