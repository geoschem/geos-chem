!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tagged_o3_mod.F90
!
! !DESCRIPTION: Contains variables and routines to perform a tagged O3
!  simulation.  P(O3) and L(O3) rates need to be archived from a full
!  chemistry simulation (via the "ProdLoss" History diagnostics collection).
!\\
!\\
! !INTERFACE:
!
MODULE Tagged_O3_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Chem_Tagged_O3
  PUBLIC  :: Init_Tagged_O3
!
! !REMARKS:
!  THE EXTENDED TAGGED O3 SIMULATION (DEFAULT) HAS THESE ADVECTED SPECIES
!  ----------------------------------------------------------------------------
!  O3      : Total O3
!  O3Strat : O3 from the Stratosphere      (tropopause - atm top   )
!  O3ut    : O3 produced in Upper Trop     (350 hPa    - tropopause)
!  O3mt    : O3 produced in Middle Trop    (PBL top    - 350 hPa   )
!  O3row   : O3 produced in Rest of World  (surface    - PBL top   )
!  O3pcbl  : O3 produced in Pacific BL     (surface    - PBL top   )
!  O3nabl  : O3 produced in N. American BL (surface    - PBL top   )
!  O3atbl  : O3 produced in Atlantic BL    (surface    - PBL top   )
!  O3eubl  : O3 produced in European BL    (surface    - PBL top   )
!  O3afbl  : O3 produced in N. African BL  (surface    - PBL top   )
!  O3asbl  : O3 produced in Asian          (surface    - PBL top   )
!  O3init  : O3 initial conditions         (all levels             )
!  O3usa   : O3 produced over the USA      (all levels             )
!                                                                             .
!  THE SIMPLE TAGGED O3 SIMULATION HAS THESE ADVECTED SPECIES:
!  -----------------------------------------------------------------------------
!  O3      : Total O3
!  O3Strat : Stratospheric O3
!                                                                             .
!  NOTES:
!  ----------------------------------------------------------------------------
!  When starting a long tagged O3 simulation, we recommend that you use
!  a restart file where all species concentrations are set to zero.
!  Then spin up for as many years as it takes to get into steady-state.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID's
  INTEGER  :: id_O3,      id_O3Strat,  id_O3ut,    id_O3mt,    id_O3row
  INTEGER  :: id_O3pcbl,  id_O3nabl,   id_O3atbl,  id_O3eubl,  id_O3afbl
  INTEGER  :: id_O3asbl,  id_O3init,   id_O3usa

  ! Global variables
  REAL(fp) :: molcm3_to_kgm3   ! Conversion factor [molec/cm3] -> [kg/m3]

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TagSpcId
!
! !DESCRIPTION: Returns the tagged species ID corresponding to the given
!  (X,Y) horizontal position and level L in the atmosphere.
!\\
!\\
! !INTERFACE:
!
  FUNCTION TagSpcId( X, Y, L, InTrop ) RESULT( tagId )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X        ! Longitude value [deg]
    REAL(fp), INTENT(IN) :: Y        ! Latitude  value [deg]
    INTEGER,  INTENT(IN) :: L        ! Level index [1]
    LOGICAL,  INTENT(IN) :: InTrop   ! =T if we are in the troposphere
                                     ! =F otherwise
!
! !RETURN VALUE:
!
    INTEGER              :: tagId    ! # of the region
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL  :: Its_In_TROP, Its_In_PBL, Its_In_MT
    LOGICAL  :: Its_In_UT,   Its_In_NH,  Its_In_ATL
    LOGICAL  :: Its_In_PAC,  Its_In_AS,  Its_In_EUR
    LOGICAL  :: Its_In_NAM,  Its_In_NAF, Its_In_USA
    INTEGER  :: PblTop,      MtTop

    !========================================================================
    ! Region begins here!
    !========================================================================

    ! Return 0 by default
    tagId = 0

    ! Return if we are not in the troposphere
    IF ( .not. InTrop ) RETURN

    !========================================================================
    ! Set up logicals that define each of the extended tagged regions
    !========================================================================

    ! PBLTOP is the model level at ~ 750 hPa
    ! MTTOP  is the model level at ~ 350 hPa
    PblTop     = 16
    MtTop      = 27

    ! Define flags for various geographic & altitude regions
    Its_In_PBL = ( L <= PblTop                                              )
    Its_In_MT  = ( L >  PblTop  .and. L <= MtTop                            )
    Its_In_UT  = ( L >  MtTop   .and. InTrop                                )
    Its_In_NH  = ( Y >=   0.0                                               )
    Its_In_EUR = ( Y >=  36.0   .and. ( X >  -15.0 .and. X <=   55.0 )      )
    Its_In_NAM = ( Y >=  15.0   .and. ( X > -127.5 .and. X <=  -65.0 )      )
    Its_In_AS  = ( Y >= -10.0   .and. ( X >   55.0 .and. X <=  145.0 )      )
    Its_In_ATL = ( Its_In_NH    .and. ( X >  -65.0 .and. X <=  -15.0 )      )
    Its_In_PAC = ( Its_In_NH    .and. ( X >  145.0  .or. X <= -127.5 )      )
    Its_In_NAF = ( ( X >= -15.0 .and. X <=  55.0 ) .and.                     &
                   ( Y >=   0.0 .and. Y <   36.0 )                          )

    !========================================================================
    ! Return the tagged species ID corresponding to the given location
    !========================================================================
    IF ( Its_In_UT ) THEN
       tagId = id_O3ut                              ! Upper trop

    ELSE IF ( Its_In_MT ) THEN
       tagId = id_O3mt                              ! Middle trop

    ELSE IF ( Its_In_PAC .and. Its_In_PBL ) THEN
       tagId = id_O3pcbl                            ! Pacific PBL

    ELSE IF ( Its_In_NAM .and. Its_In_PBL ) THEN
       tagId = id_O3nabl                            ! N. Am. PBL

    ELSE IF ( Its_In_ATL .and. Its_In_PBL ) THEN
       tagId = id_O3atbl                            ! Atlantic PBL

    ELSE IF ( Its_In_EUR .and. Its_In_PBL ) THEN
       tagId = id_O3eubl                            ! European PBL

    ELSE IF ( Its_In_NAF .and. Its_In_PBL ) THEN
       tagId = id_O3afbl                            ! N. African PBL

    ELSE IF ( Its_In_AS  .and. Its_In_PBL ) THEN
       tagId = id_O3asbl                            ! Asian PBL

    ELSE
       tagId = id_O3row                             ! Rest of world

    ENDIF

  END FUNCTION TagSpcId
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConUS
!
! !DESCRIPTION: Indicates if a lon/lat position is over the Continental USA.
!  If so, returns the species ID of the USA tagged species (id_O3usa)
!\\
!\\
! !INTERFACE:
!
  FUNCTION ConUS( X, Y, inTrop ) RESULT( isConUS )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: X         ! Longitude value [deg]
    REAL(fp), INTENT(IN) :: Y         ! Latitude value [deg]
    LOGICAL,  INTENT(IN) :: inTrop    ! =T if we are in the troposphere
                                      ! =F otherwise
!
! !RETURN VALUE:
!
    INTEGER              :: isConUS   ! = 1 if over Continental US
!                                     ! = 0 otherwise
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Return 0 by default
    isConUS = 0

    ! Exit if we are not in the tropopshere
    IF ( .not. InTrop ) RETURN

    ! If we are over CONUS, return the ID of the O3usa species
    IF ( ( X > -127.5_fp .and. X <= -65.0_fp )   .and. &
         ( Y >   22.0_fp .and. Y <=  50.0_fp ) ) THEN
       isConUS = id_O3usa
    ENDIF

  END FUNCTION ConUS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_tagged_o3
!
! !DESCRIPTION: Performs chemistry (by applying archived prod and loss rates)
!  for several O3 species tagged by geographic and altitude regions.  This
!  is useful for attributing where O3 is being produced in the atmosphere.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Tagged_O3( Input_Opt,  State_Chm, State_Diag,              &
                             State_Grid, State_Met, RC                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE Hco_Utilities_GC_Mod, ONLY : Hco_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Chem
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I,      J,      L,        N
    REAL(fp)           :: dtChem, LO3_kg, LO3_kgps, PO3_kg, PO3_kgps

    ! Arrays
    INTEGER            :: GeoMask(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    INTEGER            :: UsaMask(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp)           :: PO3_hco(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp)           :: LO3_hco(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Strings
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Chem_Tagged_O3 begins here!
    !========================================================================

    ! Initialize
    RC      =  GC_SUCCESS            ! Success or failure
    Spc     => State_Chm%Species     ! Ptr to State_
    geoMask = 0                      ! Geo-tagged O3 species numbers
    usaMask = 0                      ! Mask for CONUS-produced O3
    PO3_hco =  0.0_fp                ! Array for P(O3) from HEMCO
    LO3_hco =  0.0_fp                ! Array for L(O3) from HEMCO
    dtChem  =  Get_Ts_Chem()         ! Chemistry timestep [s]
    errMsg  =  ''
    thisLoc =  &
     ' -> at Chem_Tagged_O3 (in module GeosCore/tagged_o3_mod.F90)'

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_Loss ) State_Diag%Loss = 0.0_f4
    IF ( State_Diag%Archive_Prod ) State_Diag%Prod = 0.0_f4

    !========================================================================
    ! Get production and loss frequencies from HEMCO. These are read
    ! from the "ProdLoss" History collection (archived from a fullchem
    ! simulation) and have units [molec/cm3/s].
    !========================================================================

    ! P(O3) from HEMCO [molec/cm3/s]
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'O3_PROD', PO3_hco, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get O3_PROD [molec/cm3/s]!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! L(O3) from HEMCO [molec/cm3/s]
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'O3_LOSS', LO3_hco, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get O3_LOSS [molec/cm3/s]!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Tagged O3 chemistry contains the following terms:
    !
    !   New O3 = Old O3 + ( P(O3,region) - L(O3) )
    !
    ! P(O3) and L(O3) are archived from a previous fullchem run using
    ! the ProdLoss collection from the History diagnostics.
    !========================================================================

    ! Loop over the # of advected species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N, LO3_kg, LO3_kgps, PO3_kg, PO3_kgps           )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 4                                              )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=====================================================================
       ! Convert P(O3) and L(O3) from [molec/cm3/s] -> [kg]
       ! cf. https://github.com/geoschem/geos-chem/issues/1109 by Xingpei Ye
       !=====================================================================

       ! P(O3) in [kg/s] and [kg]
       PO3_kgps = PO3_hco(I,J,L)              & ! in molec/cm3/s
                * molcm3_to_kgm3              & ! molec/cm3/s -> kg/m3/s
                * State_Met%AIRVOL(I,J,L)       ! kg/m3/s     -> kg/s
       PO3_kg   = PO3_kgps * dtChem             ! kg/s        -> kg

       ! L(O3) in [kg/s] and [kg]
       LO3_kgps = LO3_hco(I,J,L)              & ! in molec/cm3/s
                * molcm3_to_kgm3              & ! molec/cm3/s -> kg/m3/s
                * State_Met%AIRVOL(I,J,L)       ! kg/m3/s     -> kg/s
       LO3_kg  = LO3_kgps * dtchem              ! kg/s        -> kg

       ! Prevent denormal values
       PO3_kg   = MAX( PO3_kg,   1.0e-30_fp )
       PO3_kgps = MAX( PO3_kgps, 1.0e-30_fp )
       LO3_kg   = MAX( LO3_kg,   1.0e-30_fp )
       LO3_kgps = MAX( LO3_kgps, 1.0e-30_fp )

       !------------------------------------------------------------------
       ! Find species IDs corresponding to geographic location & altitude
       ! (Extended tagged O3 simulaton only)
       !------------------------------------------------------------------
       IF ( Input_Opt%LSplit ) THEN

          ! O3usa
          usaMask(I,J,L) = ConUS(                                            &
             X      = State_Grid%XMid(I,J),                                  &
             Y      = State_Grid%YMid(I,J),                                  &
             inTrop = State_Met%InTroposphere(I,J,L)                        )

          ! All other species except O3init
          geoMask(I,J,L) = TagSpcId(                                         &
             X      = State_Grid%XMid(I,J),                                  &
             Y      = State_Grid%YMid(I,J),                                  &
             L      = L,                                                     &
             inTrop = State_Met%InTroposphere(I,J,L)                        )

       ENDIF

       !=====================================================================
       ! Apply chemical production of ozone (only where it is produced)
       !=====================================================================

       ! Add P(O3) [kg] to the total O3 species
       Spc(id_O3)%Conc(I,J,L) = Spc(id_O3)%Conc(I,J,L) + PO3_kg

       ! Add P(O3) [kg] to the stratospheric O3 species
       IF ( State_Met%InStratosphere(I,J,L) ) THEN
          Spc(id_O3Strat)%Conc(I,J,L) = Spc(id_O3Strat)%Conc(I,J,L) + PO3_kg
       ENDIF

       ! Add P(O3) to extended tagged O3 species
       IF ( Input_Opt%LSPLIT ) THEN

          ! Add P(O3) [kg] to the tagged species corresponding to the
          ! geographic/altitude region at this grid box (I,J,L)
          ! These regions only are defined in the troposphere (N > 0).
          ! Also, do not apply
          N = GeoMask(I,J,L)
          IF ( N > 0 ) THEN
             Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) + PO3_kg
          ENDIF

          ! Add P(O3) [kg] to the O3usa species, if we are within
          ! the continental USA and below the tropopause.
          N = UsaMask(I,J,L)
          IF ( N == id_O3usa ) THEN
             Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) + PO3_kg
          ENDIF
       ENDIF

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Archive chemical loss of tagged O3 species [kg/s]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_Prod ) THEN

          ! Total P(O3) and Stratospheric P(O3)
          State_Diag%Prod(I,J,L,id_O3     ) = PO3_kgps
          State_Diag%Prod(I,J,L,id_O3Strat) = PO3_kgps

          ! Archve P(O3) for extended tagged O3 species to History
          IF ( Input_Opt%LSplit ) THEN

             ! P(O3) over the continental US
             IF ( usaMask(I,J,L) == id_O3usa ) THEN
                State_Diag%Prod(I,J,L,id_O3usa) = PO3_kgps
             ENDIF

             ! P(O3) over the continental USA
             N  = GeoMask(I,J,L)
             IF ( N > 0 ) THEN
                State_Diag%Prod(I,J,L,N) = PO3_kgps
             ENDIF
          ENDIF
        ENDIF

       !=====================================================================
       ! Apply chemical loss of ozone (everywhere)
       !=====================================================================
       DO N = 1, State_Chm%nAdvect

          ! Do not apply loss to the O3init species,
          ! which preserves the initial conditions.
          IF ( N /= id_O3init ) THEN

             ! Apply chemical loss [kg]
             Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) - LO3_kg

             ! Prevent denormal values
             IF ( Spc(N)%Conc(I,J,L) < 1.0e-30_fp ) Spc(N)%Conc(I,J,L) = 0.0_fp

             !---------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Archive chemical loss of tagged O3 species [kg/s]
             !---------------------------------------------------------------
             IF ( State_Diag%Archive_Loss ) THEN
                State_Diag%Loss(I,J,L,N) = LO3_kgps
             ENDIF
          ENDIF
       ENDDO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE Chem_Tagged_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tagged_o3
!
! !DESCRIPTION: Gets species indices defines the conversion molcm3_to_kgm3
!  conversion factor.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Tagged_O3( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AVO
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
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Init_Tagged_O3 begins here
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Init_Tagged_O3 (in module GeosCore/tagged_o3_mod.F90)'

    ! Exit immediately if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    !------------------------------------------------------------------------
    ! Get O3 and O3Strat species indices
    !------------------------------------------------------------------------
    id_O3 = Ind_('O3')
    IF ( id_O3 < 0 ) THEN
       errMsg = 'O3 is an undefined species!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Define species ID flag for strat O3 (should be 2)
    id_O3Strat = Ind_('O3Strat')
    IF ( id_O3Strat < 0 ) THEN
       errMsg = 'O3Strat is an undefined species!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Conversion factor from molec/cm3 to kg/m3
    ! cf https://github.com/geoschem/geos-chem/issues/1109
    !------------------------------------------------------------------------
    molcm3_to_kgm3 = ( State_Chm%SpcData(id_O3)%Info%MW_g * 1000.0_fp ) / AVO

    !------------------------------------------------------------------------
    ! Define ID's for extended tagged O3 species
    !------------------------------------------------------------------------
    IF ( Input_Opt%LSplit ) THEN

       id_O3ut = Ind_('O3ut')
       IF ( id_O3ut < 0 ) THEN
          errMsg = 'O3ut is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3mt = Ind_('O3mt')
       IF ( id_O3ut < 0 ) THEN
          errMsg = 'O3mt is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3row = Ind_('O3row')
       IF ( id_O3row < 0 ) THEN
          errMsg = 'O3row is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3pcbl = Ind_('O3pcbl')
       IF ( id_O3pcbl < 0 ) THEN
          errMsg = 'O3pcbl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3nabl = Ind_('O3nabl')
       IF ( id_O3nabl < 0 ) THEN
          errMsg = 'O3nabl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3atbl = Ind_('O3atbl')
       IF ( id_O3atbl < 0 ) THEN
          errMsg = 'O3atbl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3eubl = Ind_('O3eubl')
       IF ( id_O3eubl < 0 ) THEN
          errMsg = 'O3eubl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3afbl = Ind_('O3afbl')
       IF ( id_O3afbl < 0 ) THEN
          errMsg = 'O3afbl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3asbl = Ind_('O3asbl')
       IF ( id_O3asbl < 0 ) THEN
          errMsg = 'O3asbl is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3init = Ind_('O3init')
       IF ( id_O3init < 0 ) THEN
          errMsg = 'O3init is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       id_O3usa = Ind_('O3usa')
       IF ( id_O3usa < 0 ) THEN
          errMsg = 'O3usa is an undefined species!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Init_Tagged_O3
!EOC
END MODULE Tagged_O3_Mod
