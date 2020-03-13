!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pops_mod.F90
!
! !DESCRIPTION: Module POPS\_MOD contains variables and routines for the
!  GEOS-Chem peristent organic pollutants (POPs) simulation.
!\\
!\\
! !INTERFACE:
!
MODULE POPS_MOD
!
! !USES:
!
  USE PhysConstants    ! For physical constants
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GetPopsDiagsFromHemco
  PUBLIC :: ChemPOPs
  PUBLIC :: Init_POPs
  PUBLIC :: Cleanup_POPs
!
! !REMARKS:
!  POPs Tracers
!  ============================================================================
!  (1 ) POPG   : Gaseous POP - total tracer
!  (2 ) POPPOC : OC-sorbed POP  - total tracer
!  (3 ) POPPBC : BC-sorbed POP  - total tracer
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin    - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL(fp), PARAMETER  :: SMALLNUM  = 1e-20_fp

  !=================================================================
  ! MODULE VARIABLES
  !=================================================================

  ! Arrays
  !=================================================================
  ! TCOSZ     : Sum of COS(Solar Zenith Angle ) [unitless]
  ! TTDAY     : Total daylight time at location (I,J) [minutes]
  ! ZERO_DVEL : Array with zero dry deposition velocity [cm/s]
  ! COSZM     : Max daily value of COS(S.Z. angle) [unitless]
  !=================================================================
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)
  REAL(fp), ALLOCATABLE :: TTDAY(:,:)
  REAL(fp), ALLOCATABLE :: ZERO_DVEL(:,:)
  REAL(fp), ALLOCATABLE :: COSZM(:,:)

  ! Pointers to fields in the HEMCO data structure.
  ! These need to be declared REAL(f4), aka REAL*4.
  ! NOTE: These are globally SAVEd pointers, so we can
  ! nullify them in the declaration here (bmy, 4/29/16)
  REAL(f4), POINTER     :: C_OC(:,:,:) => NULL()
  REAL(f4), POINTER     :: C_BC(:,:,:) => NULL()
  REAL(f4), POINTER     :: O3(:,:,:)   => NULL()
  REAL(f4), POINTER     :: OH(:,:,:)   => NULL()
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER,  PRIVATE     :: id_POPG
  INTEGER,  PRIVATE     :: id_POPPBCPI
  INTEGER,  PRIVATE     :: id_POPPBCPO
  INTEGER,  PRIVATE     :: id_POPPOCPI
  INTEGER,  PRIVATE     :: id_POPPOCPO

  ! Species drydep ID flags
  INTEGER,  PRIVATE     :: dd_POPG
  INTEGER,  PRIVATE     :: dd_POPP_BCPI
  INTEGER,  PRIVATE     :: dd_POPP_BCPO
  INTEGER,  PRIVATE     :: dd_POPP_OCPI
  INTEGER,  PRIVATE     :: dd_POPP_OCPO

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetPOPsDiagsfromHemco
!
! !DESCRIPTION: Copies manually-archived diagnostic values for the POPs
!  specialty simulation from the HEMCO state object into the relevant fields
!  of the State\_Diag object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetPopsDiagsFromHemco( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_Mod,  ONLY : GetHcoDiagn
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine should be called from EMISSIONS_RUN.  It does not add any
!  emissions for the POPS simulation per se.  But because several POPs
!  diagnostics are archived via HEMCO, we need to copy information out
!  of the HEMCO state and into the relevant fields of State_Diag.
!
! !REVISION HISTORY:
!  15 Oct 2018 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=63)  :: DgnName
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Pointers
    REAL(f4), POINTER  :: Ptr2D(:,:)

    !=================================================================
    ! GetPopsDiagsFromHemco begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = '-> at GetPopsDiagsFromHemco (in module GeosCore/pops_mod.F90)'

    !=================================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Get pointers to the POPs diagnostics that are tracked by HEMCO
    ! and then attach them to fields of the State_Diag array.
    !
    ! NOTE: The diagnostic names (which start with GCPOPS) have
    ! to match those in HEMCO/Extensions/hcox_gc_POPS_mod.F90.
    !=================================================================

    !-----------------------------------------------------------------
    ! Primary POPPOCPO emissions
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPPOCPO ) THEN

       ! Get pointer from HEMCO diagnostics
       DgnName = 'GCPOPS_POPPOCPO_SOURCE'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPPOCPO = Ptr2D

       ! Free Pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Primary POPPBCPO emissions
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPPBCPO ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_POPPBCPO_SOURCE'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPPBCPO = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Primary POPG emissions
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPG ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_POPG_SOURCE'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPG = Ptr2D

       ! Free Pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Secondary POPG emissions from soil
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPGfromSoil ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_POPG_SOIL'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPGfromSoil = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Secondary POPG emissions from lakes
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPGfromLake ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_POPG_LAKE'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPGfromLake = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Secondary POPG emissions from leaves
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_EmisPOPGfromLeaf ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_POPG_LEAF'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%EmisPOPGfromLeaf= Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Positive POPG soil flux (from soil to air)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromSoilToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_SOIL2AIR'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromSoilToAir = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Negative POPG soil flux (from air to soil)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromAirToSoil ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_AIR2SOIL'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromAirToSoil = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Positive POPG lake flux (from lake to air)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromLakeToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_LAKE2AIR'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromLakeToAir = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Negative POPG lake flux (from air to soil)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromAirToLake ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_AIR2LAKE'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromAirToLake = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Positive POPG leaf flux (from leaves to air)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromLeafToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_LEAF2AIR'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromLeafToAir = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Negative POPG leaf flux (from air to leaves)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxPOPGfromAirToLeaf ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_AIR2LEAF'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FluxPOPGfromAirToLeaf = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Fugacity ratio: soil/air
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FugacitySoilToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_SOILAIR_FUG'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FugacitySoilToAir = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Fugacity ratio: lake/air
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FugacityLakeToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_LAKEAIR_FUG'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FugacityLakeToAir = Ptr2D

       ! Free pointer
       Ptr2D => NULL()
    ENDIF

    !-----------------------------------------------------------------
    ! Fugacity ratio: leaf/air
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FugacityLeafToAir ) THEN

       ! Get pointer from HEMCO
       DgnName = 'GCPOPS_LEAFAIR_FUG'
       CALL GetHcoDiagn( DgnName, .TRUE., RC, Ptr2D=Ptr2D )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find HEMCO field: ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Diag
       State_Diag%FugacityLeafToAir = Ptr2D

       ! Free pointer
       Ptr2D   => NULL()
    ENDIF

  END SUBROUTINE GetPopsDiagsFromHemco
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  chempops
!
! !DESCRIPTION: Subroutine CHEMPOPS is the driver routine for POPs chemistry
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMPOPS( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,          ONLY : DEBUG_MSG
    USE HCO_Interface_Mod,  ONLY : GetHcoDiagn
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
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
!  20 September 2010 - N.E. Selin - Initial Version based on CHEMMERCURY
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg

    ! Pointers
    REAL(fp),      POINTER  :: DEPSAV(:,:,:  )

    !=================================================================
    ! CHEMPOPS begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ErrMsg   = ''
    ThisLoc  = ' -> at ChemPops (in module GeosCore/pops_mod.F90)'

    ! Point to columns of derived-type object fields (hplin, 12/1/18)
    DEPSAV            => State_Chm%DryDepSav

    !-----------------------------------------------------------------
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Because we need to sum the drydep from all levels in the PBL
    ! into State_Diag%DryDepChm, we need to zero it at the start
    ! of this routine.  This will avoid data from the last timestep
    ! from being accumulated in the averaging. (bmy, 10/23/18)
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_DryDepChm .or. &
         State_Diag%Archive_DryDep    ) THEN
       State_Diag%DryDepChm = 0.0_f4
    ENDIF

    !=================================================================
    ! Get pointers to fields that are read in by HEMCO.
    ! These are global concentrations of OH, O3, OC, and BC.
    !
    ! NOTE: The HEMCO pointers will update with time, so we only need
    ! to make this assignment on the first call to CHEMPOPS.
    !=================================================================
    IF ( FIRST ) THEN

       ! OC
       CALL HCO_GetPtr( HcoState, 'GLOBAL_OC', C_OC, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OC!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! BC
       CALL HCO_GetPtr( HcoState, 'GLOBAL_BC', C_BC, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_BC!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! OH
       CALL HCO_GetPtr( HcoState, 'GLOBAL_OH', OH,   RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OC!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! O3
       CALL HCO_GetPtr( HcoState, 'GLOBAL_O3', O3,   RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OC!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Perform chemistry on POPs tracers
    !=================================================================

    ! Compute diurnal scaling for OH
    CALL OHNO3TIME( State_Grid )
    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMPOPS: a OHNO3TIME' )

    !-------------------------
    ! GAS AND PARTICLE PHASE chemistry
    !-------------------------
    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMPOPS: b CHEM_GASPART' )

    ! Add option for non-local PBL (cdh, 08/27/09)
    IF ( Input_Opt%LNLPBL ) THEN

       ! Dry deposition occurs with PBL mixing,
       ! pass zero deposition frequency
       CALL CHEM_POPGP( Input_Opt,  &
                        State_Chm,  &
                        State_Diag, &
                        State_Grid, &
                        State_Met,  &
                        ZERO_DVEL,  &
                        ZERO_DVEL,  &
                        ZERO_DVEL,  &
                        ZERO_DVEL,  &
                        ZERO_DVEL,  &
                        RC )

    ELSE

       ! For addition of hydrophilic OC/BC POP tracers: for now, assume that
       ! if dry deposition of hydrophoblic tracers is active, then dry
       ! deposition of hydrophilic tracers is also active
       ! (ie., they're coupled). clf, 2/12/2012
       IF ( dd_POPG      > 0 .and. &
            dd_POPP_OCPO > 0 .and. &
            dd_POPP_BCPO > 0 ) THEN

          ! Dry deposition active for both POP-Gas and POP-Particle;
          ! pass drydep frequency to CHEM_POPGP (NOTE: DEPSAV has units 1/s)
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           DEPSAV(:,:,dd_POPG),      &
                           DEPSAV(:,:,dd_POPP_OCPO), &
                           DEPSAV(:,:,dd_POPP_OCPI), &
                           DEPSAV(:,:,dd_POPP_BCPO), &
                           DEPSAV(:,:,dd_POPP_BCPI), &
                           RC)

       ELSEIF ( dd_POPG      >  0 .and. &
                dd_POPP_OCPO >  0 .and. &
                dd_POPP_BCPO <= 0 ) THEN

          ! Only POPG and POPP_OC dry deposition are active
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           DEPSAV(:,:,dd_POPG),      &
                           DEPSAV(:,:,dd_POPP_OCPO), &
                           DEPSAV(:,:,dd_POPP_OCPI), &
                           ZERO_DVEL,                &
                           ZERO_DVEL,                &
                           RC )

       ELSEIF ( dd_POPG      >  0 .and. &
                dd_POPP_OCPO <= 0 .and. &
                dd_POPP_BCPO >  0 ) THEN

          ! Only POPG and POPP_BC dry deposition are active
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           DEPSAV(:,:,dd_POPG),      &
                           ZERO_DVEL,                &
                           ZERO_DVEL,                &
                           DEPSAV(:,:,dd_POPP_BCPO), &
                           DEPSAV(:,:,dd_POPP_BCPI), &
                           RC )

       ELSEIF ( dd_POPG      >  0 .and. &
                dd_POPP_OCPO <= 0 .and. &
                dd_POPP_BCPO <= 0 ) THEN

          ! Only POPG dry deposition is active
          CALL CHEM_POPGP( Input_Opt,           &
                           State_Chm,           &
                           State_Diag,          &
                           State_Grid,          &
                           State_Met,           &
                           DEPSAV(:,:,dd_POPG), &
                           ZERO_DVEL,           &
                           ZERO_DVEL,           &
                           ZERO_DVEL,           &
                           ZERO_DVEL,           &
                           RC )

       ELSEIF ( dd_POPG      <= 0 .and. &
                dd_POPP_OCPO >  0 .and. &
                dd_POPP_BCPO >  0 ) THEN

          ! Only POPP dry deposition is active
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           ZERO_DVEL,                &
                           DEPSAV(:,:,dd_POPP_OCPO), &
                           DEPSAV(:,:,dd_POPP_OCPI), &
                           DEPSAV(:,:,dd_POPP_BCPO), &
                           DEPSAV(:,:,dd_POPP_BCPI), &
                           RC )

       ELSEIF ( dd_POPG      <= 0 .and. &
                dd_POPP_OCPO >  0 .and. &
                dd_POPP_BCPO <= 0 ) THEN

          ! Only POPP_OC dry deposition is active
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           ZERO_DVEL,                &
                           DEPSAV(:,:,dd_POPP_OCPO), &
                           DEPSAV(:,:,dd_POPP_OCPI), &
                           ZERO_DVEL,                &
                           ZERO_DVEL,                &
                           RC )

       ELSEIF ( dd_POPG      <= 0 .and. &
                dd_POPP_OCPO <= 0 .and. &
                dd_POPP_BCPO >  0 ) THEN

          ! Only POPP_OC dry deposition is active
          CALL CHEM_POPGP( Input_Opt,                &
                           State_Chm,                &
                           State_Diag,               &
                           State_Grid,               &
                           State_Met,                &
                           ZERO_DVEL,                &
                           ZERO_DVEL,                &
                           ZERO_DVEL,                &
                           DEPSAV(:,:,dd_POPP_BCPO), &
                           DEPSAV(:,:,dd_POPP_BCPI), &
                           RC )

       ELSE

          ! No dry deposition, pass zero deposition frequency
          CALL CHEM_POPGP( Input_Opt,  &
                           State_Chm,  &
                           State_Diag, &
                           State_Grid, &
                           State_Met,  &
                           ZERO_DVEL,  &
                           ZERO_DVEL,  &
                           ZERO_DVEL,  &
                           ZERO_DVEL,  &
                           ZERO_DVEL,  &
                           RC )
       ENDIF

    ENDIF

    IF ( prtDebug ) CALL DEBUG_MSG( 'CHEMPOPS: a CHEM_GASPART' )

    ! Nullify pointers
    NULLIFY( DEPSAV )

  END SUBROUTINE CHEMPOPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  chem_popgp
!
! !DESCRIPTION: Subroutine CHEM\_POPGP is the chemistry subroutine for the
!  oxidation, gas-particle partitioning, and deposition of POPs.
!  (eck, clf, 1/4/2011)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_POPGP ( Input_Opt,                  &  
                          State_Chm,    State_Diag,   & 
                          State_Grid,   State_Met,    &
                          V_DEP_G,                    &
                          V_DEP_P_OCPO, V_DEP_P_OCPI, &
                          V_DEP_P_BCPO, V_DEP_P_BCPI, &
                          RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG53_MOD
#endif
    USE ErrCode_Mod
    USE ERROR_MOD,      ONLY : DEBUG_MSG
    USE ERROR_MOD,      ONLY : SAFE_DIV
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

    ! Dry deposition frequencies [/s]
    REAL(fp),       INTENT(IN)    :: V_DEP_G(State_Grid%NX,State_Grid%NY)
    REAL(fp),       INTENT(IN)    :: V_DEP_P_OCPO(State_Grid%NX,State_Grid%NY)
    REAL(fp),       INTENT(IN)    :: V_DEP_P_BCPO(State_Grid%NX,State_Grid%NY)
    REAL(fp),       INTENT(IN)    :: V_DEP_P_OCPI(State_Grid%NX,State_Grid%NY)
    REAL(fp),       INTENT(IN)    :: V_DEP_P_BCPI(State_Grid%NX,State_Grid%NY)
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
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) For OH rate constant: Brubaker & Hites. 1998. OH reaction kinetics of
!  PAHs and PCDD/Fs. J. Phys. Chem. A. 102:915-921.
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on CHEM_HG0_HG2
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! DENS_OCT = density of octanol, needed for partitioning into OC
    ! 820 [kg/m^3]
    REAL(fp), PARAMETER :: DENS_OCT   = 82e+1_fp

    ! DENS_BC = density of BC, needed for partitioning onto BC
    ! 1 [kg/L] or 1000 [kg/m^3]
    ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
    REAL(fp), PARAMETER :: DENS_BC    = 1e+3_fp

    ! k for reaction POPP + NO3 taken from Liu et al. EST 2012 "Kinetic
    ! studies of heterogeneous reactions of PAH aerosols with NO3 radicals",
    ! for now. CLF, 1/24/2012
    ! For PYR, 6.4 x 10^-12 [cm3 / molec / s]
    REAL(fp), PARAMETER :: K_POPP_NO3 = 0e+0_fp ! 6.4d-12

    ! OC/BC hydrophobic POP lifetime before folding to hydrophilic
    REAL(fp), PARAMETER :: OCBC_LIFE  = 1.15e+0_fp
!
! !LOCAL VARIABLES:
!
    INTEGER           :: I,             J,             L
    INTEGER           :: N
    REAL(fp)          :: DTCHEM,        SUM_F
    REAL(fp)          :: KOA_T,         KBC_T
    REAL(fp)          :: KOC_BC_T,      KBC_OC_T
    REAL(fp)          :: AREA_CM2
    REAL(fp)          :: F_PBL
    REAL(fp)          :: C_OH,          C_OC_CHEM,     C_BC_CHEM
    REAL(fp)          :: C_OC_CHEM1,    C_BC_CHEM1
    REAL(fp)          :: K_OH
    REAL(fp)          :: K_OX,          C_O3
    REAL(fp)          :: E_KOX_T
    REAL(fp)          :: K_DEPG,        K_DEPP_OCPO,   K_DEPP_BCPO
    REAL(fp)          :: K_DEPP_OCPI,   K_DEPP_BCPI
    REAL(fp)          :: OLD_POPG,      OLD_POPP_OCPO, OLD_POPP_BCPO
    REAL(fp)          :: OLD_POPP_OCPI, OLD_POPP_BCPI
    REAL(fp)          :: NEW_POPG,      NEW_POPP_OCPO, NEW_POPP_BCPO
    REAL(fp)          :: NEW_POPP_OCPI, NEW_POPP_BCPI
    REAL(fp)          :: POPG_BL,       POPP_OCPO_BL,  POPP_BCPO_BL
    REAL(fp)          :: POPP_OCPI_BL,  POPP_BCPI_BL
    REAL(fp)          :: POPG_FT,       POPP_OCPO_FT,  POPP_BCPO_FT
    REAL(fp)          :: POPP_OCPI_FT,  POPP_BCPI_FT
    REAL(fp)          :: TMP_POPG,      TMP_OX
    REAL(fp)          :: GROSS_OX,      GROSS_OX_OH,   NET_OX
    REAL(fp)          :: DEP_POPG,      DEP_POPP_OCPO, DEP_POPP_BCPO
    REAL(fp)          :: DEP_POPP_OCPI, DEP_POPP_BCPI
    REAL(fp)          :: DEP_POPG_DRY,  DEP_POPP_OCPO_DRY
    REAL(fp)          :: DEP_POPP_BCPO_DRY, DEP_POPP_OCPI_DRY
    REAL(fp)          :: DEP_POPP_BCPI_DRY
    REAL(fp)          :: DEP_DRY_FLXG,  DEP_DRY_FLXP_OCPO
    REAL(fp)          :: DEP_DRY_FLXP_BCPO, DEP_DRY_FLXP_OCPI
    REAL(fp)          :: DEP_DRY_FLXP_BCPI
    REAL(fp)          :: OLD_POP_T
    REAL(fp)          :: VR_OC_AIR,     VR_BC_AIR
    REAL(fp)          :: VR_OC_BC,      VR_BC_OC
    REAL(fp)          :: F_POP_OC,      F_POP_BC
    REAL(fp)          :: F_POP_G
    REAL(fp)          :: MPOP_OCPO,     MPOP_BCPO,   MPOP_G
    REAL(fp)          :: MPOP_OCPI,     MPOP_BCPI
    REAL(fp)          :: MPOP_OC,       MPOP_BC
    REAL(fp)          :: DIFF_G,        DIFF_OC,     DIFF_BC
    REAL(fp)          :: OC_AIR_RATIO,  OC_BC_RATIO, BC_AIR_RATIO
    REAL(fp)          :: BC_OC_RATIO,   SUM_DIFF
    REAL(fp)          :: FOLD, KOCBC,   NEW_OCPI, NEW_BCPI
    REAL(fp)          :: GROSS_OX_OCPO, GROSS_OX_OCPI
    REAL(fp)          :: GROSS_OX_BCPO, GROSS_OX_BCPI
    REAL(fp)          :: GROSS_OX_O3_OCPO, GROSS_OX_O3_OCPI
    REAL(fp)          :: GROSS_OX_O3_BCPO, GROSS_OX_O3_BCPI
    REAL(fp)          :: GROSS_OX_NO3_OCPO, GROSS_OX_NO3_OCPI
    REAL(fp)          :: GROSS_OX_NO3_BCPO, GROSS_OX_NO3_BCPI
    REAL(fp)          :: TMP_POPP_OCPO, TMP_POPP_OCPI
    REAL(fp)          :: E_KOX_T_BC
    REAL(fp)          :: TMP_OX_P_OCPO, TMP_OX_P_OCPI
    REAL(fp)          :: TMP_OX_P_BCPO, TMP_OX_P_BCPI
    REAL(fp)          :: TMP_POPP_BCPO, TMP_POPP_BCPI
    REAL(fp)          :: NET_OX_OCPO,   NET_OX_OCPI
    REAL(fp)          :: NET_OX_BCPO,   NET_OX_BCPI
    REAL(fp)          :: K_O3_BC, C_NO3, K_OX_P, K_NO3_BC

    ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
    ! from gas phase to OC. For now we use Delta H for phase transfer
    ! from the gas phase to the pure liquid state.
    ! For PHENANTHRENE:
    ! this is taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
    ! For PYRENE:
    ! this is taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].
    ! For BENZO[a]PYRENE:
    ! this is also taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
    REAL(fp)            :: DEL_H

    ! KOA_298 for partitioning of gas phase POP to atmospheric OC
    ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K
    ! For PHENANTHRENE:
    ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
    ! For PYRENE:
    ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
    ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
    REAL(fp)            :: KOA_298

    ! KBC_298 for partitioning of gas phase POP to atmospheric BC
    ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
    ! For PHENANTHRENE:
    ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
    ! For PYRENE:
    ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
    ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
    REAL(fp)            :: KBC_298

    ! K for reaction POPG + OH  [cm3 /molecule /s]
    ! For PHENANTHRENE: 2.70d-11
    ! (Source: Brubaker & Hites, J. Phys Chem A 1998)
    ! For PYRENE: 5.00d-11
    ! Calculated by finding the ratio between kOH of phenanthrene and
    ! kOH of pyrene using structure-activity relationships (Schwarzenback,
    ! Gschwend, Imboden, pg 680) and scaling the experimental kOH for
    ! phenanthrene from Brubaker and Hites
    ! For BENZO[a]PYRENE: 5.68d-11
    ! Calculated by finding the ratio between kOH of phenanthrene and
    ! kOH of pyrene using structure-activity relationships (Schwarzenback,
    ! Gschwend, Imboden, pg 680) and scaling the experimental kOH for
    ! phenanthrene from Brubaker and Hites
    ! Could potentially set this to change with temperature

    ! Using EPA AOPWIN values:
    ! For PHENANTHRENE: 13d-12
    ! For PYRENE: 50d-12
    ! For BaP: 50d-12
    REAL(fp)            :: K_POPG_OH !(Gas phase)

    ! k for reaction POPP + O3 [/s] depends on fitting parameters A and B.
    ! A represents the maximum number of surface sites available to O3, and B
    ! represents the ratio of desorption/adsorption rate coefficients for
    ! both bulk phases (Ref: Kahan et al Atm Env 2006, 40:3448)
    ! k(obs) = A x [O3(g)] / (B + [O3(g)])
    ! For PHENANTHRENE: A = 0.5 x 10^-3 s^-1, B = 2.15 x 10^15 molec/cm3
    ! For PYRENE: A = 0.7 x 10^-3 s^-1, B = 3 x 10^15 molec/cm3
    ! for BaP: A = 5.5 x 10^-3 s^-1, B = 2.8 x 10^15 molec/cm3
    REAL(fp)            :: AK  ! s^-1
    REAL(fp)            :: BK  ! molec/cm3

    ! For fields from Input_Opt
    LOGICAL             :: LNLPBL
    LOGICAL             :: LGTMM
    LOGICAL             :: Archive_Drydep

    ! For SAFE_DIV
    REAL(fp)            :: DENOM

    ! Pointers
    REAL(fp),      POINTER :: Spc(:,:,:,:)
    TYPE(Species), POINTER :: ThisSpc

    !=================================================================
    ! CHEM_POPGP begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy values from Input_Opt
    DEL_H     = Input_Opt%POP_DEL_H
    KOA_298   = Input_Opt%POP_KOA
    KBC_298   = Input_Opt%POP_KBC
    K_POPG_OH = Input_Opt%POP_K_POPG_OH
    AK        = Input_Opt%POP_K_POPP_O3A
    BK        = Input_Opt%POP_K_POPP_O3B
    LNLPBL    = Input_Opt%LNLPBL
    LGTMM     = Input_Opt%LGTMM

    ! Point to the chemical species array [kg]
    Spc       => State_Chm%Species

    ! Pointer for the species database object
    ThisSpc   => NULL()

    ! Chemistry timestep [s]
    DTCHEM = GET_TS_CHEM()

    !================================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Zero out all diagnostic arrays to avoid leftover values from
    ! from the last timestep from being archived again (bmy,10/15/18)
    !================================================================

    !----------------------------------------------------------------
    ! Initialize State_Diag fields for prod/loss diagnostics
    !----------------------------------------------------------------
    IF ( State_Diag%Archive_LossPOPPOCPObyGasPhase ) THEN
       State_Diag%LossPOPPOCPObyGasPhase = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPOCPOfromGasPhase ) THEN
       State_Diag%ProdPOPPOCPOfromGasPhase = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossPOPPBCPObyGasPhase ) THEN
       State_Diag%LossPOPPBCPObyGasPhase = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPBCPOfromGasPhase ) THEN
       State_Diag%ProdPOPPBCPOfromGasPhase = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPGfromOH ) THEN
       State_Diag%ProdPOPGfromOH = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPOCPOfromO3 ) THEN
       State_Diag%ProdPOPPOCPOfromO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPOCPIfromO3 ) THEN
       State_Diag%ProdPOPPOCPIfromO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPBCPOfromO3 ) THEN
       State_Diag%ProdPOPPBCPOfromO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPBCPIfromO3 ) THEN
       State_Diag%ProdPOPPBCPIfromO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPOCPOfromNO3 ) THEN
       State_Diag%ProdPOPPOCPOfromNO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPOCPIfromNO3 ) THEN
       State_Diag%ProdPOPPOCPIfromNO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPBCPOfromNO3 ) THEN
       State_Diag%ProdPOPPBCPOfromNO3 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdPOPPBCPIfromNO3 ) THEN
       State_Diag%ProdPOPPBCPIfromNO3 = 0.0_f4
    ENDIF

    !----------------------------------------------------------------
    ! Determine if we need to save drydep to netCDF
    !----------------------------------------------------------------
    Archive_Drydep = ( State_Diag%Archive_DryDepChm .or. &
                       State_Diag%Archive_DryDep         )

    ! Eventually should parallelize this ...
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NZ

       ! Zero concentrations in loop
       MPOP_G            = 0e+0_fp
       MPOP_OCPO         = 0e+0_fp
       MPOP_BCPO         = 0e+0_fp
       MPOP_OCPI         = 0e+0_fp
       MPOP_BCPI         = 0e+0_fp
       OLD_POPG          = 0e+0_fp
       OLD_POPP_OCPO     = 0e+0_fp
       OLD_POPP_BCPO     = 0e+0_fp
       OLD_POPP_OCPI     = 0e+0_fp
       OLD_POPP_BCPI     = 0e+0_fp
       OLD_POP_T         = 0e+0_fp
       NEW_POPG          = 0e+0_fp
       NEW_POPP_OCPO     = 0e+0_fp
       NEW_POPP_BCPO     = 0e+0_fp
       NEW_POPP_OCPI     = 0e+0_fp
       NEW_POPP_BCPI     = 0e+0_fp
       POPG_BL           = 0e+0_fp
       POPP_OCPO_BL      = 0e+0_fp
       POPP_BCPO_BL      = 0e+0_fp
       POPP_OCPI_BL      = 0e+0_fp
       POPP_BCPI_BL      = 0e+0_fp
       POPG_FT           = 0e+0_fp
       POPP_OCPO_FT      = 0e+0_fp
       POPP_BCPO_FT      = 0e+0_fp
       POPP_OCPI_FT      = 0e+0_fp
       POPP_BCPI_FT      = 0e+0_fp
       DIFF_G            = 0e+0_fp
       DIFF_OC           = 0e+0_fp
       DIFF_BC           = 0e+0_fp
       NET_OX            = 0e+0_fp
       TMP_POPG          = 0e+0_fp
       TMP_OX            = 0e+0_fp
       GROSS_OX          = 0e+0_fp
       GROSS_OX_OH       = 0e+0_fp
       DEP_POPG          = 0e+0_fp
       DEP_POPP_OCPO     = 0e+0_fp
       DEP_POPP_BCPO     = 0e+0_fp
       DEP_POPP_OCPI     = 0e+0_fp
       DEP_POPP_BCPI     = 0e+0_fp
       DEP_POPG_DRY      = 0e+0_fp
       DEP_POPP_OCPO_DRY = 0e+0_fp
       DEP_POPP_BCPO_DRY = 0e+0_fp
       DEP_POPP_OCPI_DRY = 0e+0_fp
       DEP_POPP_BCPI_DRY = 0e+0_fp
       DEP_DRY_FLXG      = 0e+0_fp
       DEP_DRY_FLXP_OCPO = 0e+0_fp
       DEP_DRY_FLXP_BCPO = 0e+0_fp
       DEP_DRY_FLXP_OCPI = 0e+0_fp
       DEP_DRY_FLXP_BCPI = 0e+0_fp
       E_KOX_T           = 0e+0_fp
       K_OX              = 0e+0_fp
       K_O3_BC           = 0e+0_fp
       GROSS_OX_OCPO     = 0e+0_fp
       GROSS_OX_BCPO     = 0e+0_fp
       GROSS_OX_OCPI     = 0e+0_fp
       GROSS_OX_BCPI     = 0e+0_fp
       GROSS_OX_O3_OCPO  = 0e+0_fp
       GROSS_OX_O3_BCPO  = 0e+0_fp
       GROSS_OX_O3_OCPI  = 0e+0_fp
       GROSS_OX_O3_BCPI  = 0e+0_fp
       GROSS_OX_NO3_OCPO = 0e+0_fp
       GROSS_OX_NO3_BCPO = 0e+0_fp
       GROSS_OX_NO3_OCPI = 0e+0_fp
       GROSS_OX_NO3_BCPI = 0e+0_fp
       TMP_POPP_OCPO     = 0e+0_fp
       TMP_POPP_OCPI     = 0e+0_fp
       TMP_POPP_BCPO     = 0e+0_fp
       TMP_POPP_BCPI     = 0e+0_fp
       E_KOX_T_BC        = 0e+0_fp
       TMP_OX_P_OCPO     = 0e+0_fp
       TMP_OX_P_OCPI     = 0e+0_fp
       TMP_OX_P_BCPO     = 0e+0_fp
       TMP_OX_P_BCPI     = 0e+0_fp
       NET_OX_OCPO       = 0e+0_fp
       NET_OX_OCPI       = 0e+0_fp
       NET_OX_BCPO       = 0e+0_fp
       NET_OX_BCPI       = 0e+0_fp

       ! Get monthly mean OH concentrations
       C_OH  = GET_OH( I, J, L, State_Met )

       ! Get monthly mean O3 concentrations
       ! O3 is in v/v (from HEMCO), convert to molec/cm3 (mps, 9/10/14)
       C_O3  = O3(I,J,L) * ( AVO / ( AIRMW * 1.e-3_fp ) ) * &
               State_Met%AIRDEN(I,J,L) * 1e-6_fp

       ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
       F_PBL = State_Met%F_UNDER_PBLTOP(I,J,L)

       ! Define K for the oxidation reaction with POPG [/s]
       K_OH  = K_POPG_OH * C_OH

       ! Define K for the oxidation reaction with POPPOC and POPPBC [/s]
       ! Kahan:
       K_O3_BC  = (AK * C_O3) / (BK + C_O3)

       ! Bug fix: Zero C_NO3 to avoid a floating-point exception error.
       ! K_POPP_NO3 is already set to zero, so this will already cause
       ! K_NO3_BC to be set to zero.  We think that C_NO3 once was
       ! assigned a value but was later taken out. (bmy, 10/6/16)
       C_NO3 = 0.0_fp

       ! Define K for oxidation of POPP by NO3  [/s]
       K_NO3_BC = K_POPP_NO3 * C_NO3

       ! Total K for gas phase oxidation [/s]
       K_OX   = K_OH !+ ...

       ! Total K for particle phase oxidation [/s]
       K_OX_P = K_O3_BC + K_NO3_BC

       ! Define Ks for dry deposition of gas phase POP [/s]
       K_DEPG = V_DEP_G(I,J)

       ! Define Ks for dry deposition of hydrophoblic particle phase POP [/s]
       K_DEPP_OCPO = V_DEP_P_OCPO(I,J)
       K_DEPP_BCPO = V_DEP_P_BCPO(I,J)

       ! Define Ks for dry deposition of hydrophilic particle phase POP [/s]
       K_DEPP_OCPI = V_DEP_P_OCPI(I,J)
       K_DEPP_BCPI = V_DEP_P_BCPI(I,J)

       ! Precompute exponential factors [dimensionless]
       ! For gas phase
       E_KOX_T    = EXP( -K_OX   * DTCHEM )
       ! For particle phase
       E_KOX_T_BC = EXP( -K_OX_P * DTCHEM )

       !==============================================================
       ! GAS-PARTICLE PARTITIONING
       !==============================================================

       OLD_POPG      = MAX( Spc(I,J,L,id_POPG),     SMALLNUM ) ![kg]
       OLD_POPP_OCPO = MAX( Spc(I,J,L,id_POPPOCPO), SMALLNUM ) ![kg]
       OLD_POPP_BCPO = MAX( Spc(I,J,L,id_POPPBCPO), SMALLNUM ) ![kg]
       OLD_POPP_OCPI = MAX( Spc(I,J,L,id_POPPOCPI), SMALLNUM ) ![kg]
       OLD_POPP_BCPI = MAX( Spc(I,J,L,id_POPPBCPI), SMALLNUM ) ![kg]

       ! Total POPs in box I,J,L
       OLD_POP_T = OLD_POPG      + &
                   OLD_POPP_OCPO + OLD_POPP_BCPO + &
                   OLD_POPP_OCPI + OLD_POPP_BCPI

       ! Define temperature-dependent partition coefficients
       ! NOTE: State_Met%T is in units of [K]
       KOA_T = KOA_298 * EXP( ( -DEL_H / RSTARG ) * &
                            ( ( 1e+0_fp / State_Met%T(I,J,L) ) - &
                            ( 1e+0_fp / 298e+0_fp ) ) )

       ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
       ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
       ! NOTE: State_Met%T is in units of [K]
       KBC_T = KBC_298 * EXP( ( -DEL_H / RSTARG ) * &
                            ( ( 1e+0_fp / State_Met%T(I,J,L) ) - &
                            ( 1e+0_fp / 298e+0_fp ) ) )

       ! Define KOC_BC_T, theoretical OC-BC part coeff at temp T [unitless]
       KOC_BC_T = KOA_T / KBC_T

       ! Define KBC_OC_T, theoretical BC_OC part coeff at temp T [unitless]
       KBC_OC_T = 1e+0_fp / KOC_BC_T

       ! Get monthly mean OC and BC concentrations [kg/box]
       C_OC_CHEM = C_OC(I,J,L)
       C_BC_CHEM = C_BC(I,J,L)

       ! Make sure OC is not negative
       C_OC_CHEM = MAX( C_OC_CHEM, 0e+0_fp )

       ! Convert to units of volume per box [m^3 OC or BC/box]
       C_OC_CHEM1 = C_OC_CHEM / DENS_OCT
       C_BC_CHEM1 = C_BC_CHEM / DENS_BC

       ! Define volume ratios:
       ! VR_OC_AIR = volume ratio of OC to air [unitless]
       VR_OC_AIR     = C_OC_CHEM1 / State_Met%AIRVOL(I,J,L) ! could be zero

       ! VR_OC_BC = volume ratio of OC to BC [unitless]
       VR_OC_BC      = SAFE_DIV( C_OC_CHEM1, C_BC_CHEM1, 0e+0_fp )

       ! VR_BC_AIR = volume ratio of BC to air [unitless]
       VR_BC_AIR     = SAFE_DIV( VR_OC_AIR,  VR_OC_BC,   0e+0_fp )

       ! VR_BC_OC = volume ratio of BC to OC [unitless]
       VR_BC_OC      = SAFE_DIV( 1e+0_fp,    VR_OC_BC,   0e+0_fp )

       ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase,
       ! BC-phase, and gas phase with new time step (should only change if
       ! temp changes or OC/BC concentrations change)
       DENOM         = KOA_T * VR_OC_AIR
       OC_AIR_RATIO  = SAFE_DIV( 1e+0_fp,  DENOM,      0e+0_fp )

       DENOM         = KOC_BC_T * VR_OC_BC
       OC_BC_RATIO   = SAFE_DIV( 1e+0_fp,  DENOM,      0e+0_fp )

       DENOM         = KBC_T * VR_BC_AIR
       BC_AIR_RATIO  = SAFE_DIV( 1e+0_fp,  DENOM,      0e+0_fp )

       DENOM         = KBC_OC_T * VR_BC_OC
       BC_OC_RATIO   = SAFE_DIV( 1e+0_fp,  DENOM,      0e+0_fp )

       ! If there are zeros in OC or BC concentrations, make sure they
       ! don't cause problems with phase fractions
       IF ( C_OC_CHEM > SMALLNUM .and. C_BC_CHEM > SMALLNUM ) THEN
          F_POP_OC  = 1e+0_fp / (1e+0_fp + OC_AIR_RATIO + OC_BC_RATIO)
          F_POP_BC  = 1e+0_fp / (1e+0_fp + BC_AIR_RATIO + BC_OC_RATIO)

       ELSE IF (C_OC_CHEM >  SMALLNUM .and. &
                C_BC_CHEM <= SMALLNUM ) THEN
          F_POP_OC  = 1e+0_fp / (1e+0_fp + OC_AIR_RATIO)
          F_POP_BC  = SMALLNUM

       ELSE IF ( C_OC_CHEM <= SMALLNUM .and. &
                 C_BC_CHEM >  SMALLNUM ) THEN
          F_POP_OC  = SMALLNUM
          F_POP_BC  = 1e+0_fp / (1e+0_fp + BC_AIR_RATIO)

       ELSE IF ( C_OC_CHEM <= SMALLNUM .and. &
                 C_BC_CHEM <= SMALLNUM) THEN
          F_POP_OC = SMALLNUM
          F_POP_BC = SMALLNUM
       ENDIF

       ! Gas-phase:
       F_POP_G = 1e+0_fp - F_POP_OC - F_POP_BC

       ! Check that sum equals 1
       SUM_F   = F_POP_OC + F_POP_BC + F_POP_G

       ! Calculate new masses of POP in each phase [kg]
       ! OC-phase:
       MPOP_OC = F_POP_OC * OLD_POP_T

       ! BC-phase
       MPOP_BC = F_POP_BC * OLD_POP_T

       ! Gas-phase
       MPOP_G  = F_POP_G  * OLD_POP_T

       ! Ensure new masses of POP in each phase are positive
       MPOP_OC = MAX(MPOP_OC, SMALLNUM)
       MPOP_BC = MAX(MPOP_BC, SMALLNUM)
       MPOP_G  = MAX(MPOP_G,  SMALLNUM)

       ! Calculate differences in masses in each phase from previous time
       ! step for storage in ND53 diagnostic
       DIFF_G  = MPOP_G  - OLD_POPG
       DIFF_OC = MPOP_OC - OLD_POPP_OCPO - OLD_POPP_OCPI
       DIFF_BC = MPOP_BC - OLD_POPP_BCPO - OLD_POPP_BCPI

       ! Sum of differences should equal zero
       SUM_DIFF = MAX(DIFF_G + DIFF_OC + DIFF_BC, SMALLNUM)

#ifdef BPCH_DIAG
       !==============================================================
       ! %%%%% ND53 (bpch) DIAGNOSTIC %%%%%
       ! ND53 diagnostic: Differences in distribution of gas and
       ! particle phases between time steps [kg/s]
       !==============================================================
       IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

          IF (DIFF_OC .lt. 0) THEN
             AD53_PG_OC_NEG(I,J,L) = AD53_PG_OC_NEG(I,J,L)  + DIFF_OC / DTCHEM
          ELSE IF (DIFF_OC .eq. 0 .or. DIFF_OC .gt. 0) THEN
             AD53_PG_OC_POS(I,J,L) = AD53_PG_OC_POS(I,J,L)  + DIFF_OC / DTCHEM
          ENDIF

          IF (DIFF_BC .lt. 0) THEN
             AD53_PG_BC_NEG(I,J,L) = AD53_PG_BC_NEG(I,J,L)  + DIFF_BC / DTCHEM
          ELSE IF (DIFF_BC .eq. 0 .or. DIFF_BC .gt. 0) THEN
             AD53_PG_BC_POS(I,J,L) = AD53_PG_BC_POS(I,J,L)  + DIFF_BC / DTCHEM
          ENDIF

       ENDIF
#endif

       !==============================================================
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Prod and loss of total POPPOC and total POPPBC species
       ! by gas-phase reactions [kg/s]
       !==============================================================

       IF ( DIFF_OC < 0.0_fp ) THEN

          ! Loss of POPPOCPO by gas-phase reactions
          IF ( State_Diag%Archive_LossPOPPOCPObyGasPhase ) THEN
             State_Diag%LossPOPPOCPObyGasPhase(I,J,L)   = DIFF_OC / DTCHEM
          ENDIF

       ELSE

          ! Prod of POPPOCPO by gas-phase reactions
          IF ( State_Diag%Archive_ProdPOPPOCPOfromGasPhase ) THEN
             State_Diag%ProdPOPPOCPOfromGasPhase(I,J,L) = DIFF_OC / DTCHEM
          ENDIF

       ENDIF

       IF ( DIFF_BC < 0.0_fp ) THEN

          ! Loss of POPPBCPO by gas-phase reactions
          IF ( State_Diag%Archive_LossPOPPBCPObyGasPhase ) THEN
             State_Diag%LossPOPPBCPObyGasPhase(I,J,L)   = DIFF_BC / DTCHEM
          ENDIF

       ELSE

          ! Prod of POPPBCPO by gas-phase reactions
          IF ( State_Diag%Archive_ProdPOPPOCPOfromGasPhase ) THEN
             State_Diag%ProdPOPPBCPOfromGasPhase(I,J,L) = DIFF_BC / DTCHEM
          ENDIF

       ENDIF

       !==============================================================
       ! HYDROPHOBIC PARTICULATE POPS DECAY TO HYDROPHILIC (clf, 2/12/2012)
       !==============================================================

       ! Define the lifetime and e-folding time to be the same as hydrophobic
       ! to  hydrophilic aerosols
       KOCBC = 1.e+0_fp / (86400e+0_fp * OCBC_LIFE)

       ! Send hydrophobic to hydrophilic
       FOLD = KOCBC * DTCHEM

       ! Amount of hydrophobic particulate POP left after folding to
       ! hydrophilic
       MPOP_OCPO = MPOP_OC * EXP( -FOLD )
       MPOP_BCPO = MPOP_BC * EXP( -FOLD )

       ! Hydrophilic particulate POP already existing
       MPOP_OCPI = MAX( Spc(I,J,L,id_POPPOCPI), SMALLNUM )  ![kg]
       MPOP_BCPI = MAX( Spc(I,J,L,id_POPPBCPI), SMALLNUM )  ![kg]

       ! Hydrophilic POP that used to be hydrophobic
       NEW_OCPI = MPOP_OC - MPOP_OCPO
       NEW_BCPI = MPOP_BC - MPOP_BCPO

       ! Add new hydrophilic to old hydrophilic
       ! Don't do this - already added into total for redistribution (clf 3/27/2012)
       ! MPOP_OCPI = MPOP_OCPI + NEW_OCPI
       ! MPOP_BCPI = MPOP_BCPI + NEW_BCPI
       MPOP_OCPI = NEW_OCPI
       MPOP_BCPI = NEW_BCPI

       !==============================================================
       ! CHEMISTRY AND DEPOSITION REACTIONS
       !==============================================================
       IF ( F_PBL < 0.05e+0_fp .OR. K_DEPG < SMALLNUM ) THEN

          !==============================================================
          ! Entire box is in the free troposphere
          ! or deposition is turned off, so use RXN without deposition
          ! for gas phase POPs
          ! For particle POPs, no rxn and no deposition
          !==============================================================

          ! OH:
          CALL RXN_OX_NODEP( MPOP_G, K_OX, &
                             E_KOX_T, NEW_POPG, GROSS_OX )

          ! O3 and NO3:
          CALL RXN_OX_NODEP( MPOP_OCPO, K_OX_P, &
                             E_KOX_T_BC, NEW_POPP_OCPO, GROSS_OX_OCPO)

          CALL RXN_OX_NODEP( MPOP_OCPI, K_OX_P, &
                             E_KOX_T_BC, NEW_POPP_OCPI, GROSS_OX_OCPI)

          CALL RXN_OX_NODEP( MPOP_BCPO, K_OX_P, &
                             E_KOX_T_BC, NEW_POPP_BCPO, GROSS_OX_BCPO)

          CALL RXN_OX_NODEP( MPOP_BCPI, K_OX_P, &
                             E_KOX_T_BC, NEW_POPP_BCPI, GROSS_OX_BCPI)

          ! No deposition occurs [kg]
          DEP_POPG      = 0e+0_fp
          DEP_POPP_OCPO = 0e+0_fp
          DEP_POPP_BCPO = 0e+0_fp
          DEP_POPP_OCPI = 0e+0_fp
          DEP_POPP_BCPI = 0e+0_fp

       ELSE IF ( F_PBL > 0.95e+0_fp ) THEN

          !==============================================================
          ! Entire box is in the boundary layer
          ! so use RXN with deposition for gas phase POPs
          ! Deposition only (no rxn) for particle phase POPs
          !==============================================================

          CALL RXN_OX_WITHDEP( MPOP_G,        K_OX,          &
                               K_DEPG,        DTCHEM,        &
                               E_KOX_T,       NEW_POPG,      &
                               GROSS_OX,      DEP_POPG       )

          CALL RXN_OX_WITHDEP( MPOP_OCPO,     K_OX_P,        &
                               K_DEPP_OCPO,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_OCPO, &
                               GROSS_OX_OCPO, DEP_POPP_OCPO  )

          CALL RXN_OX_WITHDEP( MPOP_OCPI,     K_OX_P,        &
                               K_DEPP_OCPI,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_OCPI, &
                               GROSS_OX_OCPI, DEP_POPP_OCPI  )

          CALL RXN_OX_WITHDEP( MPOP_BCPO,     K_OX_P,        &
                               K_DEPP_BCPO,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_BCPO, &
                               GROSS_OX_BCPO, DEP_POPP_BCPO  )

          CALL RXN_OX_WITHDEP( MPOP_BCPI,     K_OX_P,        &
                               K_DEPP_BCPI,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_BCPI, &
                               GROSS_OX_BCPI, DEP_POPP_BCPI  )

       ELSE

          !==============================================================
          ! Box spans the top of the boundary layer
          ! Part of the mass is in the boundary layer and subject to
          ! deposition while part is in the free troposphere and
          ! experiences no deposition.
          !
          ! We apportion the mass between the BL and FT according to the
          ! volume fraction of the box in the boundary layer.
          ! Arguably we should assume uniform mixing ratio, instead of
          ! uniform density but if the boxes are short, the air density
          ! doesn't change much.
          ! But assuming uniform mixing ratio across the inversion layer
          ! is a poor assumption anyway, so we are just using the
          ! simplest approach.
          !==============================================================

          ! Boundary layer portion of POPG [kg]
          POPG_BL = MPOP_G * F_PBL

          ! Boundary layer portion of POPP_OCPO [kg]
          POPP_OCPO_BL = MPOP_OCPO * F_PBL

          ! Boundary layer portion of POPP_BCPO [kg]
          POPP_BCPO_BL = MPOP_BCPO * F_PBL

          ! Boundary layer portion of POPP_OCPI [kg]
          POPP_OCPI_BL = MPOP_OCPI * F_PBL

          ! Boundary layer portion of POPP_BCPI [kg]
          POPP_BCPI_BL = MPOP_BCPI * F_PBL

          ! Free troposphere portion of POPG [kg]
          POPG_FT = MPOP_G - POPG_BL

          ! Free troposphere portion of POPP_OCPO [kg]
          POPP_OCPO_FT = MPOP_OCPO - POPP_OCPO_BL

          ! Free troposphere portion of POPP_BCPO [kg]
          POPP_BCPO_FT = MPOP_BCPO - POPP_BCPO_BL

          ! Free troposphere portion of POPP_OCPI [kg]
          POPP_OCPI_FT = MPOP_OCPI - POPP_OCPI_BL

          ! Free troposphere portion of POPP_BCPI [kg]
          POPP_BCPI_FT = MPOP_BCPI - POPP_BCPI_BL

          ! Do chemistry with deposition on BL fraction for gas phase
          CALL RXN_OX_WITHDEP( POPG_BL,  K_OX,     &
                               K_DEPG,   DTCHEM,   &
                               E_KOX_T,  NEW_POPG, &
                               GROSS_OX, DEP_POPG  )

          ! Do chemistry without deposition on the FT fraction for gas phase
          CALL RXN_OX_NODEP(   POPG_FT,  K_OX,               &
                               E_KOX_T,  TMP_POPG, TMP_OX    )

          ! Now do the same with the OC and BC phase:

          ! Do chemistry with deposition on BL fraction for OCPO phase
          CALL RXN_OX_WITHDEP( POPP_OCPO_BL,  K_OX_P,        &
                               K_DEPP_OCPO,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_OCPO, &
                               GROSS_OX_OCPO, DEP_POPP_OCPO  )

          ! Do chemistry without deposition on the FT fraction for OCPO phase
          CALL RXN_OX_NODEP(   POPP_OCPO_FT,  K_OX_P,        &
                               E_KOX_T_BC,    TMP_POPP_OCPO, &
                               TMP_OX_P_OCPO  )

          ! Do chemistry with deposition on BL fraction for OCPI phase
          CALL RXN_OX_WITHDEP( POPP_OCPI_BL,  K_OX_P,        &
                               K_DEPP_OCPI,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_OCPI, &
                               GROSS_OX_OCPI, DEP_POPP_OCPI  )

          ! Do chemistry without deposition on the FT fraction for OCPI phase
          CALL RXN_OX_NODEP(   POPP_OCPI_FT,  K_OX_P,        &
                               E_KOX_T_BC,    TMP_POPP_OCPI, &
                               TMP_OX_P_OCPI  )

          ! Do chemistry with deposition on BL fraction for BCPO phase
          CALL RXN_OX_WITHDEP( POPP_BCPO_BL,  K_OX_P,        &
                               K_DEPP_BCPO,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_BCPO, &
                               GROSS_OX_BCPO, DEP_POPP_BCPO  )

          ! Do chemistry without deposition on the FT fraction for BCPO phase
          CALL RXN_OX_NODEP(   POPP_BCPO_FT,  K_OX_P,        &
                               E_KOX_T_BC,    TMP_POPP_BCPO, &
                               TMP_OX_P_BCPO  )

          ! Do chemistry with deposition on BL fraction for BCPI phase
          CALL RXN_OX_WITHDEP( POPP_BCPI_BL,  K_OX_P,        &
                               K_DEPP_BCPI,   DTCHEM,        &
                               E_KOX_T_BC,    NEW_POPP_BCPI, &
                               GROSS_OX_BCPI, DEP_POPP_BCPI  )

          ! Do chemistry without deposition on the FT fraction for BCPI phase
          CALL RXN_OX_NODEP(   POPP_BCPI_FT,  K_OX_P,        &
                               E_KOX_T_BC,    TMP_POPP_BCPI, &
                               TMP_OX_P_BCPI  )

          ! Do deposition (no chemistry) on BL fraction for particulate phase
          ! No deposition (and no chem) on the FT fraction
          ! for the particulate phase
          !CALL NO_RXN_WITHDEP(POPP_OCPO_BL, K_DEPP_OCPO, DTCHEM, &
          !                    NEW_POPP_OCPO, DEP_POPP_OCPO)

          !CALL NO_RXN_WITHDEP(POPP_BCPO_BL, K_DEPP_BCPO, DTCHEM, &
          !                    NEW_POPP_BCPO, DEP_POPP_BCPO)

          !CALL NO_RXN_WITHDEP(POPP_OCPI_BL, K_DEPP_OCPI, DTCHEM, &
          !                    NEW_POPP_OCPI, DEP_POPP_OCPI)

          !CALL NO_RXN_WITHDEP(POPP_BCPI_BL, K_DEPP_BCPI, DTCHEM, &
          !                    NEW_POPP_BCPI, DEP_POPP_BCPI)

          ! Recombine the boundary layer and free troposphere parts [kg]
          NEW_POPG      = NEW_POPG      + TMP_POPG
          NEW_POPP_OCPO = NEW_POPP_OCPO + TMP_POPP_OCPO
          NEW_POPP_BCPO = NEW_POPP_BCPO + TMP_POPP_BCPO
          NEW_POPP_OCPI = NEW_POPP_OCPI + TMP_POPP_OCPI
          NEW_POPP_BCPI = NEW_POPP_BCPI + TMP_POPP_BCPI

          ! Total gross oxidation of gas phase in the BL and FT [kg]
          GROSS_OX      = GROSS_OX      + TMP_OX

          ! Total gross oxidation of particle phases in the BL and FT [kg]
          GROSS_OX_OCPO = GROSS_OX_OCPO + TMP_OX_P_OCPO
          GROSS_OX_OCPI = GROSS_OX_OCPI + TMP_OX_P_OCPI
          GROSS_OX_BCPO = GROSS_OX_BCPO + TMP_OX_P_BCPO
          GROSS_OX_BCPI = GROSS_OX_BCPI + TMP_OX_P_BCPI

       ENDIF

       ! Ensure positive concentration [kg]
       NEW_POPG      = MAX( NEW_POPG, SMALLNUM )
       NEW_POPP_OCPO = MAX( NEW_POPP_OCPO, SMALLNUM )
       NEW_POPP_BCPO = MAX( NEW_POPP_BCPO, SMALLNUM )
       NEW_POPP_OCPI = MAX( NEW_POPP_OCPI, SMALLNUM )
       NEW_POPP_BCPI = MAX( NEW_POPP_BCPI, SMALLNUM )

       ! Archive new POPG and POPP values [kg]
       Spc(I,J,L,id_POPG)     = NEW_POPG
       Spc(I,J,L,id_POPPOCPO) = NEW_POPP_OCPO
       Spc(I,J,L,id_POPPBCPO) = NEW_POPP_BCPO
       Spc(I,J,L,id_POPPOCPI) = NEW_POPP_OCPI
       Spc(I,J,L,id_POPPBCPI) = NEW_POPP_BCPI

       ! Net oxidation [kg] (equal to gross ox for now)
       NET_OX      = MPOP_G    - NEW_POPG      - DEP_POPG
       NET_OX_OCPO = MPOP_OCPO - NEW_POPP_OCPO - DEP_POPP_OCPO
       NET_OX_OCPI = MPOP_OCPI - NEW_POPP_OCPI - DEP_POPP_OCPI
       NET_OX_BCPO = MPOP_BCPO - NEW_POPP_BCPO - DEP_POPP_BCPO
       NET_OX_BCPI = MPOP_BCPI - NEW_POPP_BCPI - DEP_POPP_BCPI

       ! Error check on gross oxidation [kg]
       IF ( GROSS_OX < 0e+0_fp ) &
            CALL DEBUG_MSG('CHEM_POPGP: negative gross gas oxid')

       IF ( GROSS_OX_OCPO < 0e+0_fp ) &
            CALL DEBUG_MSG('CHEM_POPGP: negative gross OCPO oxid')

       IF ( GROSS_OX_OCPI < 0e+0_fp ) &
            CALL DEBUG_MSG('CHEM_POPGP: negative gross OCPI oxid')

       IF ( GROSS_OX_BCPO < 0e+0_fp ) &
            CALL DEBUG_MSG('CHEM_POPGP: negative gross BCPO oxid')

       IF ( GROSS_OX_BCPI < 0e+0_fp ) &
            CALL DEBUG_MSG('CHEM_POPGP: negative gross BCPI oxid')

       ! Apportion gross oxidation between OH (and no other gas-phase
       ! oxidants considered now) [kg]
       IF ( (K_OX < SMALLNUM) .OR. (GROSS_OX < SMALLNUM) ) THEN
          GROSS_OX_OH = 0e+0_fp
       ELSE
          GROSS_OX_OH = GROSS_OX * K_OH / K_OX
       ENDIF

       ! Small number check for particulate O3 oxidation
       ! Now apportion total particulate oxidation between O3 and NO3
       IF ( (K_OX_P < SMALLNUM) .OR. (GROSS_OX_OCPO < SMALLNUM) ) THEN
          GROSS_OX_OCPO     = 0e+0_fp
       ELSE
          GROSS_OX_O3_OCPO  = GROSS_OX_OCPO * K_O3_BC / K_OX_P
          GROSS_OX_NO3_OCPO = GROSS_OX_OCPO * K_NO3_BC / K_OX_P
       ENDIF

       IF ( (K_OX_P < SMALLNUM) .OR. (GROSS_OX_OCPI < SMALLNUM) ) THEN
          GROSS_OX_OCPI     = 0e+0_fp
       ELSE
          GROSS_OX_O3_OCPI  = GROSS_OX_OCPI * K_O3_BC / K_OX_P
          GROSS_OX_NO3_OCPI = GROSS_OX_OCPI * K_NO3_BC / K_OX_P
       ENDIF

       IF ( (K_OX_P < SMALLNUM) .OR. (GROSS_OX_BCPO < SMALLNUM) ) THEN
          GROSS_OX_BCPO     = 0e+0_fp
       ELSE
          GROSS_OX_O3_BCPO  = GROSS_OX_BCPO * K_O3_BC / K_OX_P
          GROSS_OX_NO3_BCPO = GROSS_OX_BCPO * K_NO3_BC / K_OX_P
       ENDIF

       IF ( (K_OX_P < SMALLNUM) .OR. (GROSS_OX_BCPI < SMALLNUM) ) THEN
          GROSS_OX_BCPI     = 0e+0_fp
       ELSE
          GROSS_OX_O3_BCPI  = GROSS_OX_BCPI * K_O3_BC / K_OX_P
          GROSS_OX_NO3_BCPI = GROSS_OX_BCPI * K_NO3_BC / K_OX_P
       ENDIF

       ! Apportion deposition [kg]
       ! Right now only using dry deposition (no sea salt) (clf, 1/27/11)
       ! If ever use dep with sea salt aerosols,
       ! will need to multiply DEP_POPG by the ratio
       ! of K_DRYG (rate of dry dep) to K_DEPG (total dep rate).
       IF ( (K_DEPG < SMALLNUM) .OR. (DEP_POPG < SMALLNUM) ) THEN
          DEP_POPG_DRY  = 0e+0_fp
       ELSE
          DEP_POPG_DRY  = DEP_POPG
       ENDIF

       IF ( (K_DEPP_OCPO < SMALLNUM) .OR. (DEP_POPP_OCPO < SMALLNUM) ) THEN
          DEP_POPP_OCPO_DRY  = 0e+0_fp
       ELSE
          DEP_POPP_OCPO_DRY  = DEP_POPP_OCPO
       ENDIF

       IF ( (K_DEPP_BCPO < SMALLNUM) .OR. (DEP_POPP_BCPO < SMALLNUM) ) THEN
          DEP_POPP_BCPO_DRY  = 0e+0_fp
       ELSE
          DEP_POPP_BCPO_DRY  = DEP_POPP_BCPO
       ENDIF

       IF ( (K_DEPP_OCPI < SMALLNUM) .OR. (DEP_POPP_OCPI < SMALLNUM) ) THEN
          DEP_POPP_OCPI_DRY  = 0e+0_fp
       ELSE
          DEP_POPP_OCPI_DRY  = DEP_POPP_OCPI
       ENDIF

       IF ( (K_DEPP_BCPI < SMALLNUM) .OR. (DEP_POPP_BCPI < SMALLNUM) ) THEN
          DEP_POPP_BCPI_DRY  = 0e+0_fp
       ELSE
          DEP_POPP_BCPI_DRY  = DEP_POPP_BCPI
       ENDIF

#ifdef BPCH_DIAG
       !==============================================================
       ! %%%% ND53 (bpch) DIAGNOSTIC %%%%
       ! ND53 diagnostic: Oxidized POPG (OH-POPG) production [kg/s]
       !==============================================================
       IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

          ! OH:
          AD53_POPG_OH(I,J,L) = AD53_POPG_OH(I,J,L) + &
                                GROSS_OX_OH / DTCHEM

          ! O3:
          AD53_POPP_OCPO_O3(I,J,L) = AD53_POPP_OCPO_O3(I,J,L) + &
                                     GROSS_OX_O3_OCPO / DTCHEM
          AD53_POPP_OCPI_O3(I,J,L) = AD53_POPP_OCPI_O3(I,J,L) + &
                                     GROSS_OX_O3_OCPI / DTCHEM
          AD53_POPP_BCPO_O3(I,J,L) = AD53_POPP_BCPO_O3(I,J,L) + &
                                     GROSS_OX_O3_BCPO / DTCHEM
          AD53_POPP_BCPI_O3(I,J,L) = AD53_POPP_BCPI_O3(I,J,L) + &
                                     GROSS_OX_O3_BCPI / DTCHEM

          ! NO3:
          AD53_POPP_OCPO_NO3(I,J,L) = AD53_POPP_OCPO_NO3(I,J,L) + &
                                      GROSS_OX_NO3_OCPO / DTCHEM
          AD53_POPP_OCPI_NO3(I,J,L) = AD53_POPP_OCPI_NO3(I,J,L) + &
                                      GROSS_OX_NO3_OCPI / DTCHEM
          AD53_POPP_BCPO_NO3(I,J,L) = AD53_POPP_BCPO_NO3(I,J,L) + &
                                      GROSS_OX_NO3_BCPO / DTCHEM
          AD53_POPP_BCPI_NO3(I,J,L) = AD53_POPP_BCPI_NO3(I,J,L) + &
                                      GROSS_OX_NO3_BCPI / DTCHEM

       ENDIF
#endif

       !==============================================================
       ! %%%%% HISTORY (netCDF DIAGNOSTICS) %%%%%
       !
       ! If we are using full PBL mixing, then archive dry dep flux
       ! due to chemistry and the total dry dep flux [molec/cm2/s]
       !
       ! NOTE: Only archive drydep within the PBL (bmy, 10/23/18)
       ! ALSO NOTE: Need to sum contributions from all layers
       ! into the State_Diag%DryDepChm array.
       !==============================================================
       IF ( ( Archive_Drydep         )  .and. &
            ( .not. Input_Opt%LNLPBL )  .and. &
            ( F_PBL > 0.0_fp         ) ) THEN

          ! Grid box surface area [cm2]
          AREA_CM2 = State_Grid%Area_M2(I,J) * 1e+4_fp

          !----------------------------------
          ! Save drydep flux for POPG
          !----------------------------------
          IF ( id_POPG > 0 ) THEN

             ! Point to POPG entry in the species database
             ThisSpc => State_Chm%SpcData(id_POPG)%Info

             ! Amt of POPG lost to drydep [molec/cm2/s]
             DEP_DRY_FLXG = DEP_POPG_DRY * AVO &
                            / ( 1.e-3_fp * ThisSpc%EmMW_g ) &
                            / ( AREA_CM2 * DTCHEM         )

             ! Save into State_Diag%DryDepChm
             State_Diag%DryDepChm(I,J,id_POPG) = &
                  State_Diag%DryDepChm(I,J,id_POPG) + DEP_DRY_FLXG

             ! Free pointer
             ThisSpc => NULL()
          ENDIF

          !----------------------------------
          ! Save drydep flux for POPPOCPO
          !----------------------------------
          IF ( id_POPPOCPO > 0 ) THEN

             ! Point to POPPOCPO entry in the species database
             ThisSpc => State_Chm%SpcData(id_POPPOCPO)%Info

             ! Amt of POPPOCPO lost to drydep [molec/cm2/s]
             DEP_DRY_FLXP_OCPO = DEP_POPP_OCPO_DRY * AVO &
                                 / ( 1.e-3_fp * ThisSpc%EmMW_g ) &
                                 / ( AREA_CM2 * DTCHEM         )

             ! Save into State_Diag%DryDepChm
             State_Diag%DryDepChm(I,J,id_POPPOCPO) = &
                  State_Diag%DryDepChm(I,J,id_POPPOCPO) + DEP_DRY_FLXP_OCPO

             ! Free pointer
             ThisSpc => NULL()
          ENDIF

          !----------------------------------
          ! Save drydep flux for POPPOCPI
          !----------------------------------
          IF ( id_POPPOCPI > 0 ) THEN

             ! Point to POPPOCPO entry in the species database
             ThisSpc => State_Chm%SpcData(id_POPPOCPI)%Info

             ! Amt of POPPOCPO lost to drydep [molec/cm2/s]
             DEP_DRY_FLXP_OCPI = DEP_POPP_OCPI_DRY * AVO &
                                 / ( 1.e-3_fp * ThisSpc%EmMW_g ) &
                                 / ( AREA_CM2 * DTCHEM         )

             ! Save into State_Diag%DryDepChm
             State_Diag%DryDepChm(I,J,id_POPPOCPI) = &
                  State_Diag%DryDepChm(I,J,id_POPPOCPI) + DEP_DRY_FLXP_OCPI

             ! Free pointer
             ThisSpc => NULL()
          ENDIF

          !----------------------------------
          ! Save drydep flux for POPPBCPO
          !----------------------------------
          IF ( id_POPPBCPO > 0 ) THEN

             ! Point to POPPOCPO entry in the species database
             ThisSpc => State_Chm%SpcData(id_POPPBCPO)%Info

             ! Amt of POPPBCPO lost to drydep [molec/cm2/s]
             DEP_DRY_FLXP_BCPO = DEP_POPP_BCPO_DRY * AVO &
                                 / ( 1.e-3_fp * ThisSpc%EmMW_g ) &
                                 / ( AREA_CM2 * DTCHEM         )

             ! Save into State_Diag%DryDepChm
             State_Diag%DryDepChm(I,J,id_POPPBCPO) = &
                  State_Diag%DryDepChm(I,J,id_POPPBCPO) + DEP_DRY_FLXP_BCPO

             ! Free pointer
             ThisSpc => NULL()
          ENDIF

          !----------------------------------
          ! Save drydep flux for POPPBCPI
          !----------------------------------
          IF ( id_POPPBCPI > 0 ) THEN

             ! Point to POPPOCPO entry in the species database
             ThisSpc => State_Chm%SpcData(id_POPPBCPI)%Info

             ! Amt of POPPBCPI lost to drydep [molec/cm2/s]
             DEP_DRY_FLXP_BCPI = DEP_POPP_BCPI_DRY * AVO &
                                 / ( 1.e-3_fp * ThisSpc%EmMW_g ) &
                                 / ( AREA_CM2 * DTCHEM         )

             ! Save into State_Diag%DryDepChm
             State_Diag%DryDepChm(I,J,id_POPPBCPI) = &
                  State_Diag%DryDepChm(I,J,id_POPPBCPI) + DEP_DRY_FLXP_BCPI

             ! Free pointer
             ThisSpc => NULL()
          ENDIF
       ENDIF

       !==============================================================
       ! %%%%% HISTORY (netCDF) DIAGNOSTICS %%%%%
       !
       ! POPs prod/loss diagnostics [kg/s]
       !==============================================================
       IF ( State_Diag%Archive_ProdPOPGfromOH ) THEN
          State_Diag%ProdPOPGfromOH(I,J,L)      = GROSS_OX_OH / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPOCPOfromO3 ) THEN
          State_Diag%ProdPOPPOCPOfromO3(I,J,L)  = GROSS_OX_O3_OCPO / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPOCPIfromO3 ) THEN
          State_Diag%ProdPOPPOCPIfromO3(I,J,L)  = GROSS_OX_O3_OCPI / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPBCPOfromO3 ) THEN
          State_Diag%ProdPOPPBCPOfromO3(I,J,L)  = GROSS_OX_O3_BCPO / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPBCPIfromO3 ) THEN
          State_Diag%ProdPOPPBCPIfromO3(I,J,L)  = GROSS_OX_O3_BCPI / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPOCPOfromNO3 ) THEN
          State_Diag%ProdPOPPOCPOfromNO3(I,J,L) = GROSS_OX_NO3_OCPO / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPOCPIfromNO3 ) THEN
          State_Diag%ProdPOPPOCPIfromNO3(I,J,L) = GROSS_OX_NO3_OCPI / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPBCPOfromNO3 ) THEN
          State_Diag%ProdPOPPBCPOfromNO3(I,J,L) = GROSS_OX_O3_BCPO / DTCHEM
       ENDIF

       IF ( State_Diag%Archive_ProdPOPPBCPIfromNO3 ) THEN
          State_Diag%ProdPOPPBCPIfromNO3(I,J,L) = GROSS_OX_NO3_BCPI / DTCHEM
       ENDIF

    ENDDO
    ENDDO
    ENDDO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_POPGP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  rxn_ox_nodep
!
! !DESCRIPTION: Subroutine RXN\_OX\_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on
! RXN\_REDOX\_NODEP in mercury\_mod.F90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RXN_OX_NODEP( OLD_POPG, K_OX, E_KOX_T, NEW_POPG, GROSS_OX )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)  :: OLD_POPG
    REAL(fp),  INTENT(IN)  :: K_OX
    REAL(fp),  INTENT(IN)  :: E_KOX_T
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT) :: NEW_POPG
    REAL(fp),  INTENT(OUT) :: GROSS_OX
!
! !REVISION HISTORY:
!  27 January 2011 - CL Friedman - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! RXN_OX_NODEP begins here!
    !=================================================================

    ! Oxidation
    IF (K_OX < SMALLNUM ) THEN

       GROSS_OX = 0e+0_fp
       NEW_POPG = OLD_POPG

    ELSE

       ! New concentration of POPG
       NEW_POPG = OLD_POPG * E_KOX_T

       ! Gross oxidation
       GROSS_OX = OLD_POPG - NEW_POPG
       GROSS_OX = MAX( GROSS_OX, 0e+0_fp )

    ENDIF

  END SUBROUTINE RXN_OX_NODEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  rxn_ox_withdep
!
! !DESCRIPTION: Subroutine RXN\_OX\_WITHDEP calculates new mass of POPG for
!  given rates of oxidation and deposition. This is for the boundary layer.
!  (clf, 1/27/11, based on RXN\_REDOX\_NODEP in mercury\_mod.F90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RXN_OX_WITHDEP( OLD_POPG, K_OX, K_DEPG, DT, E_KOX_T, &
                             NEW_POPG, GROSS_OX, DEP_POPG )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)  :: OLD_POPG
    REAL(fp),  INTENT(IN)  :: DT
    REAL(fp),  INTENT(IN)  :: K_OX
    REAL(fp),  INTENT(IN)  :: K_DEPG
    REAL(fp),  INTENT(IN)  :: E_KOX_T
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT) :: NEW_POPG
    REAL(fp),  INTENT(OUT) :: GROSS_OX
    REAL(fp),  INTENT(OUT) :: DEP_POPG
!
! !REVISION HISTORY:
!  27 January 2011 - CL Friedman - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)               :: E_KDEPG_T
    REAL(fp)               :: NEWPOPG_OX
    REAL(fp)               :: NEWPOPG_DEP

    !=================================================================
    ! RXN_OX_WITHDEP begins here!
    !=================================================================

    ! Precompute exponential factor for deposition [dimensionless]
    E_KDEPG_T = EXP( -K_DEPG * DT )

    IF (K_OX < SMALLNUM) THEN

       !=================================================================
       ! No Chemistry, Deposition only
       !=================================================================

       ! New mass of POPG [kg]
       NEW_POPG = OLD_POPG * E_KDEPG_T

       ! Oxidation of POPG [kg]
       GROSS_OX = 0e+0_fp

       ! Deposited POPG [kg]
       DEP_POPG = OLD_POPG - NEW_POPG

    ELSE

       !=================================================================
       ! Oxidation and Deposition
       !=================================================================

       ![POPG](t) = [POPG](0) exp( -(kOx + kDPOPG) t)
       !Ox(t)     = ( [POPG](0) - [POPG](t) ) * kOx / ( kOx + kDPOPG )
       !Dep_POPG(t)   = ( [POPG](0) - [POPG](t) - Ox(t) )

       ! New concentration of POPG [kg]
       NEW_POPG = OLD_POPG * E_KOX_T * E_KDEPG_T

       ! Gross oxidized gas phase mass [kg]
       GROSS_OX = ( OLD_POPG - NEW_POPG ) * K_OX / ( K_OX + K_DEPG )
       GROSS_OX = MAX( GROSS_OX, 0e+0_fp )

       ! POPG deposition [kg]
       DEP_POPG = ( OLD_POPG - NEW_POPG - GROSS_OX )
       DEP_POPG = MAX( DEP_POPG, 0e+0_fp )

    ENDIF

  END SUBROUTINE RXN_OX_WITHDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  no_rxn_withdep
!
! !DESCRIPTION: Subroutine NO\_RXN\_WITHDEP calculates new mass of POPP for
!  given rate of deposition. No oxidation of POPP. This is for the boundary
!  layer. (clf, 2/9/11)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NO_RXN_WITHDEP( OLD_POPP, K_DEPP, DT, NEW_POPP, DEP_POPP )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)  :: OLD_POPP
    REAL(fp),  INTENT(IN)  :: K_DEPP
    REAL(fp),  INTENT(IN)  :: DT
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT) :: NEW_POPP
    REAL(fp),  INTENT(OUT) :: DEP_POPP
!
! !REVISION HISTORY:
!  9 February 2011 - CL Friedman - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)               :: E_KDEPP_T

    !=================================================================
    ! NO_RXN_WITHDEP begins here!
    !=================================================================

    ! Precompute exponential factors [dimensionless]
    E_KDEPP_T = EXP( -K_DEPP * DT )

    !=================================================================
    ! No Chemistry, Deposition only
    !=================================================================

    ! New mass of POPP [kg]
    NEW_POPP = OLD_POPP * E_KDEPP_T

    ! POPP deposition [kg]
    DEP_POPP = OLD_POPP - NEW_POPP
    DEP_POPP = MAX( DEP_POPP, 0e+0_fp )

  END SUBROUTINE NO_RXN_WITHDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_oh
!
! !DESCRIPTION: Function GET\_OH returns monthly mean OH and imposes a diurnal
! variation.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_OH( I, J, L, State_Met ) RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I, J, L
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  03 Feb 2011 - CL Friedman - Initial Version, copied from mercury_mod.f
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)        :: OH_MOLEC_CM3

    !=================================================================
    ! GET_OH begins here!
    !=================================================================

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN

       ! Impose a diurnal variation on OH during the day
       ! OH from HEMCO is in kg/m3
       OH_MOLEC_CM3 = OH(I,J,L) * 1e-6_fp * ( AVO / 0.017) * &
                      ( State_Met%SUNCOS(I,J) / TCOSZ(I,J)    ) * &
                      ( 86400e+0_fp           / GET_TS_CHEM() )

       ! Make sure OH is not negative
       OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

    ELSE

       ! At night, OH goes to zero
       OH_MOLEC_CM3 = 0e+0_fp

    ENDIF

  END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight.
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version for POPS_MOD
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: I, J, L, N, NT, NDYSTEP
    REAL(fp)            :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)            :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)            :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! OHNO3TIME begins here!
    !=================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

       ! Zero arrays
       TTDAY(:,:) = 0e+0_fp
       TCOSZ(:,:) = 0e+0_fp
       COSZM(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R( I, J )

             TIMLOC = real(LHR0) + real(NT)/3600.0 + &
                      State_Grid%XMid( I, J )/15.0

             DO WHILE (TIMLOC .lt. 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC .gt. 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() ) * 60e+0_fp
             ENDIF
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

          ! Increment elapsed time [sec]
          NT = NT + GET_TS_CHEM()
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_pops
!
! !DESCRIPTION: Subroutine INIT\_POPS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_POPS( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on INIT_MERCURY
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Assume success
    RC    = GC_SUCCESS

    !=================================================================
    ! Allocate and initialize arrays
    ! NOTE: These might have to go into state_chm_mod.F90 eventually
    !=================================================================
    ALLOCATE( COSZM( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pops_mod.F90:COSZM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    COSZM = 0.0_fp

    ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pops_mod.F90:TCOSZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOSZ = 0.0_fp

    ALLOCATE( TTDAY( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pops_mod.F90:TTDAY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TTDAY = 0.0_fp

    ! Now always allocate ZERO_DVEL.  This is needed because the calls
    ! to CHEM_POPGP have been modified in the CHEMPOPS routine above.
    ! (bmy, 3/10/15)
    ALLOCATE( ZERO_DVEL( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'pops_mod.F90:ZERO_DVEL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ZERO_DVEL = 0.0_fp

    !=================================================================
    ! Initialize species ID and drydep ID flags
    !=================================================================
    id_POPG      = Ind_('POPG'        )
    dd_POPG      = Ind_('POPG',    'D')

    id_POPPOCPO  = Ind_('POPPOCPO'    )
    dd_POPP_OCPO = Ind_('POPPOCPO','D')

    id_POPPBCPO  = Ind_('POPPBCPO'    )
    dd_POPP_BCPO = Ind_('POPPBCPO','D')

    id_POPPOCPI  = Ind_('POPPOCPI'    )
    dd_POPP_OCPI = Ind_('POPPOCPI','D')

    id_POPPBCPI  = Ind_('POPPBCPI'    )
    dd_POPP_BCPI = Ind_('POPPBCPI','D')

  END SUBROUTINE INIT_POPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cleanup_pops
!
! !DESCRIPTION: Subroutine CLEANUP\_POPS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_POPS( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  20 September 2010 - N.E. Selin - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate variables
    IF ( ALLOCATED( COSZM ) ) THEN
       DEALLOCATE( COSZM, STAT=RC )
       CALL GC_CheckVar( 'pops_mod.F90:COSZM', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( TCOSZ ) ) THEN
       DEALLOCATE( TCOSZ, STAT=RC )
       CALL GC_CheckVar( 'pops_mod.F90:TCOSZ', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( TTDAY  ) ) THEN
       DEALLOCATE( TTDAY, STAT=RC )
       CALL GC_CheckVar( 'pops_mod.F90:TTDAY', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( ZERO_DVEL ) ) THEN
       DEALLOCATE( ZERO_DVEL, STAT=RC  )
       CALL GC_CheckVar( 'pops_mod.F90:ZERO_DVEL', 2, RC )
       RETURN
    ENDIF

    ! Free pointers
    IF ( ASSOCIATED( C_OC ) ) C_OC => NULL()
    IF ( ASSOCIATED( C_BC ) ) C_BC => NULL()
    IF ( ASSOCIATED( O3   ) ) O3   => NULL()
    IF ( ASSOCIATED( OH   ) ) OH   => NULL()

  END SUBROUTINE CLEANUP_POPS
!EOC
END MODULE POPS_MOD
