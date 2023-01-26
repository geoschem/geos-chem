!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mercury_mod.F90
!
! !DESCRIPTION: Contains variables and routines for the GEOS-Chem mercury
!  simulation.  Many choices of reaction mechanism and model processes can
!  be selected with logical switches located in INIT\_MERCURY.
!\\
!\\
! !INTERFACE:
!
MODULE Mercury_Mod
!
! !USES:
!
  USE Depo_Mercury_Mod,  ONLY : ADD_Hg2_SNOWPACK
  USE Depo_Mercury_Mod,  ONLY : LHgSNOW
  USE Ocean_Mercury_Mod, ONLY : LDYNSEASALT
  USE Ocean_Mercury_Mod, ONLY : LPOLARBr
  USE Ocean_Mercury_Mod, ONLY : L_ADD_MBL_Br
  USE Ocean_Mercury_Mod, ONLY : LGEIA05
  USE Ocean_Mercury_Mod, ONLY : LVEGEMIS
  USE Ocean_Mercury_Mod, ONLY : LBrCHEM
  USE Ocean_Mercury_Mod, ONLY : LRED_CLOUDONLY
  USE Ocean_Mercury_Mod, ONLY : LHALOGENCHEM
  USE Ocean_Mercury_Mod, ONLY : LHgAQCHEM
  USE Ocean_Mercury_Mod, ONLY : LHg2HalfAerosol
  USE Ocean_Mercury_Mod, ONLY : STRAT_Br_FACTOR
  USE Ocean_Mercury_Mod, ONLY : LAnthroHgOnly
  USE Ocean_Mercury_Mod, ONLY : LOHO3CHEM
  USE Ocean_Mercury_Mod, ONLY : LGCBrOMINE
  USE Ocean_Mercury_Mod, ONLY : LnoUSAemis
  USE Ocean_Mercury_Mod, ONLY : LBrOCHEM
  USE Ocean_Mercury_Mod, ONLY : LNEI2005
  USE Ocean_Mercury_Mod, ONLY : LInPlume
  USE Ocean_Mercury_Mod, ONLY : LOCEANCOEF
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Cleanup_Mercury
  PUBLIC :: ChemMercury
  PUBLIC :: EmissMercury
  PUBLIC :: Init_Mercury
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Shah, V. et al (2021), "Improved mechanistic model of the atmospheric
!        redox chemistry of mercury", Environ. Sci. Technol., 55, 14445-14456/
!  (2 ) Saiz-Lopez, A. et al (2020), "Photochemistry of oxidized Hg(I) and
!        Hg(II) species suggests missing mercury oxidation in the
!         troposphere", PNAS, 117, 30949-3095, 2020.
!  (3 ) Parrella, J. et al. (2012), Tropospheric bromine chemistry:
!        implications for present and pre-industrial ozone and mercury, ACP.
!  (4 ) Prados-Roman, C. et al. (2011), Airborne DOAS limb measurements of
!        tropospheric trace gas profiles: case studies on the profile retrieval
!        of O4 and BrO, Atmos. Meas. Tech., 4: 1241-1260.
!  (5 ) Pohler, D. et al. (2010), Observation of halogen species in the Amundsen
!        Gulf, Arctic, by active long-path differential optical absorption
!        spectroscopy, Proc. Natl. Acad. Sci, 107(15): 6528-6587.
!  (6 ) Holmes, C.D., et al. (2010) Global atmospheric model for mercury
!        including oxidation by bromine atoms, AC&P, 10, 12,037-12,057.
!  (7 ) Streets, D.G. et al. (2009), Projections of global mercury emissions
!        in 2050, Environ. Sci. Technol., 43, 2983-2988.
!  (8 ) Corbitt, E.S. et al. (2011), Global source-receptor relationships for
!        mercury deposition under present-day and 2050 emissions scenarios,
!        Environ. Sci. Technol., 45, 10477-10484.
!  (8 ) Soerensen, A. et al. (2010), An improved global model for air-sea
!        exchange of mercury: High concentrations over the North Atlantic,
!        Environ. Sci. Technol., 44, 8574-8580.
!  (9 ) Mintz, Y and G.K. Walker (1993). "Global fields of soil moisture
!        and land surface evapotranspiration derived from observed
!        precipitation and surface air temperature." J. Appl. Meteorol. 32 (8),
!        1305-1334.
!  (10) Allison, J.D. and T.L. Allison (2005) "Partition coefficients for
!        metals in surface water, soil and waste." Rep. EPA/600/R-05/074,
!        US EPA, Office of Research and Development, Washington, D.C.
!  (11) Selin, N., et al. (2008). "Global 3-D land-ocean-atmospehre model
!        for mercury: present-day versus preindustrial cycles and
!        anthropogenic enrichment factors for deposition." Global
!        Biogeochemical Cycles 22: GB2011.
!  (12) Selin, N., et al. (2007). "Chemical cycling and deposition of
!        atmospheric mercury: Global constraints from observations."
!        J. Geophys. Res. 112.
!  (13) Sommar, J., et al. (2001). "A kinetic study of the gas-phase
!        reaction between the hydroxyl radical and atomic mercury."
!        Atmospheric Environment 35: 3049-3054.
!  (14) Hall, B. (1995). "The gas phase oxidation of elemental mercury by
!        ozone.", Water, Air, and Soil Pollution 80: 301-315.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !--------------------------------------------------------------------------
  ! Scalars
  !--------------------------------------------------------------------------
  LOGICAL  :: Failed2x
  INTEGER  :: N_Hg_CATS
  INTEGER  :: id_Hg0,         id_Hg2,      id_HgP
  INTEGER  :: id_phot_NO2,    id_phot_BrO, id_phot_ClO
  INTEGER  :: id_phot_Hg2Org, id_O3,       id_OH
  INTEGER  :: id_HO2,         id_ClO,      id_Cl
  INTEGER  :: id_NO2,         id_NO,       id_Br
  INTEGER  :: id_BrO,         id_HgBrNO2,  id_HgBrHO2
  INTEGER  :: id_HgBrOH,      id_HgBrBrO,  id_HgBrClO
  INTEGER  :: id_HgBr2,       id_HgClNO2,  id_HgClHO2
  INTEGER  :: id_HgClOH,      id_HgClBrO,  id_HgClClO
  INTEGER  :: id_HgClBr,      id_HgOHNO2,  id_HgOHHO2
  INTEGER  :: id_HgOHOH,      id_HgOHBrO,  id_HgOHClO
  INTEGER  :: id_HgCl2,       id_Hg2Clp,   id_Hg2ORGp
  INTEGER  :: id_Hg2STRP,     id_HgBr,     id_HgCl
  INTEGER  :: id_HgOH,        id_HgBrO,    id_HgClO
  INTEGER  :: id_HgOHO,       nHg2gasSpc,  n_Aer
  INTEGER  :: n_Dust
  REAL(fp) :: srMw_HgCl2

  !--------------------------------------------------------------------------
  ! Arrays
  !--------------------------------------------------------------------------
  INTEGER               :: Map_Hg2gas(25)
  INTEGER,  ALLOCATABLE :: PL_Kpp_ID(:)
  REAL(fp), ALLOCATABLE :: EHg0_an(:,:)
  REAL(fp), ALLOCATABLE :: EHg0_dist(:,:)
  REAL(fp), ALLOCATABLE :: EHg0_ln(:,:)
  REAL(fp), ALLOCATABLE :: EHg0_oc(:,:)
  REAL(fp), ALLOCATABLE :: EHg0_so(:,:)
  REAL(fp), ALLOCATABLE :: EHg0_snow(:,:)
  REAL(fp), ALLOCATABLE :: EHg2_an(:,:)
  REAL(fp), ALLOCATABLE :: COSZM(:,:)               ! Max daily SZA
  REAL(fp), ALLOCATABLE :: srMw(:)
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)               ! Sum of SZA
  REAL(fp), ALLOCATABLE :: TTDAY(:,:)               ! Total daylight time [min]
  REAL(fp), ALLOCATABLE :: ZERO_DVEL(:,:)           ! Zero drydep vel [cm/s]
  REAL(fp), ALLOCATABLE :: Hg2_SEASALT_LOSSRATE(:,:)
  CHARACTER(LEN=8),                                                          &
            ALLOCATABLE :: AerSpcNames(:)

  ! For now, we need an emission array for the HG simulation
  ! that can be passed to vdiff_mod.F90 since Trac_Tend does
  ! not exist anymore (ckeller, 10/21/2014).
  REAL(fp), ALLOCATABLE, PUBLIC :: HG_EMIS(:,:,:)

  !--------------------------------------------------------------------------
  ! Pointers to fields in the HEMCO data structure, which must be REAL*4.
  ! (NOTE: We can set them to NULL here because hey are globally
  ! SAVEd variables (bmy, 4/29/16)
  !--------------------------------------------------------------------------
  REAL(f4), POINTER :: O3(:,:,:)         => NULL()
  REAL(f4), POINTER :: OH(:,:,:)         => NULL()
  REAL(f4), POINTER :: JNO2(:,:,:)       => NULL()
  REAL(f4), POINTER :: NO2(:,:,:)        => NULL()
  REAL(f4), POINTER :: NO(:,:,:)         => NULL()
  REAL(f4), POINTER :: HOCl(:,:,:)       => NULL()
  REAL(f4), POINTER :: HO2(:,:,:)        => NULL()
  REAL(f4), POINTER :: ClO(:,:,:)        => NULL()
  REAL(f4), POINTER :: Cl(:,:,:)         => NULL()
  REAL(f4), POINTER :: OA(:,:,:)         => NULL()
  REAL(f4), POINTER :: OCEAN_CONC(:,:,:) => NULL()
  REAL(f4), POINTER :: GLOB_PM25(:,:,:)  => NULL()
  REAL(f4), POINTER :: GLOB_fOA (:,:,:)  => NULL()
  REAL(f4), POINTER :: GLOB_RH(:,:,:)    => NULL()

  !--------------------------------------------------------------------------
  ! Derived types and derived-type arrays
  !--------------------------------------------------------------------------

  ! For oxidants and related quantities read from HEMCO
  TYPE :: ConcPtrObj
     REAL(f4), POINTER :: Data(:,:,:) => NULL()  ! [molec/cm3]
  END TYPE ConcPtrObj

  ! For AOD quantities read from HEMCO
  TYPE :: AeroPtrObj
     REAL(f4), POINTER :: AOD(:,:,:)  => NULL()   ! [unitless]
     REAL(f4), POINTER :: Area(:,:,:) => NULL()   ! [cm2/cm3]
     REAL(f4), POINTER :: Radi(:,:,:) => NULL()   ! [cm]
  END TYPE AeroPtrObj

  ! Vector of type ConcPtrObj
  TYPE(ConcPtrObj), POINTER :: FixSpcPtr(:)

  ! Vector of type AeroPtrObj
  TYPE(AeroPtrObj), POINTER :: AeroPtr(:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissmercury
!
! !DESCRIPTION: Subroutine EMISSMERCURY is the driver routine for mercury
!  emissions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmissMercury( Input_Opt,  State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
!
! !USES:
!
    USE DEPO_MERCURY_MOD,     ONLY : RESET_HG_DEP_ARRAYS
    USE ErrCode_Mod
    USE ERROR_MOD
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetPtr
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Land_Mercury_Mod,     ONLY : Land_Mercury_Flux
    USE Land_Mercury_Mod,     ONLY : SOILEMIS
    USE Land_Mercury_Mod,     ONLY : SNOWPACK_MERCURY_FLUX
    USE Ocean_Mercury_Mod,    ONLY : OCEAN_MERCURY_FLUX
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : GET_MONTH, ITS_A_NEW_MONTH
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
    
    ! Added for GTMM (ccc, 11/19/09)
    !USE LAND_MERCURY_MOD,   ONLY : GTMM_DR
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
! !REMARKS:
!
!
! !REVISION HISTORY:
!  03 Jun 2013 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRST = .TRUE.
    INTEGER            :: THISMONTH, I, J
    LOGICAL            :: prtDebug
    
    ! Pointers
    REAL(f4),  POINTER :: Ptr2D(:,:)

    ! Strings
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! EMISSMERCURY begins here!
    !=================================================================

    ! Assume success
    RC       = GC_SUCCESS
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ErrMsg   = ''
    ThisLoc  = ' -> at EMISSMERCURY (in module GeosCore/mercury_mod.F90)'

    ! Convert species units to [kg] for EMISSMERCURY (ewl, 8/12/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units" #1!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Get data pointers from HEMCO on the first call
    !=================================================================
    IF ( FIRST ) THEN

       ! Soil distribution
       CALL HCO_GC_GetPtr( Input_Opt, State_Grid, 'HG0_SOILDIST', Ptr2D, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to HEMCO field HG0_SOILDIST!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       EHg0_dist =  Ptr2D(:,:)
       Ptr2D     => NULL()

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !--------------------------------------------------------------
    ! Here we are NOT using the Global Terrstrial Mercury Model...
    !--------------------------------------------------------------

    ! Read offline ocean Hg0
    EHg0_oc = 0.0_fp
    CALL OFFLINEOCEAN_READMO( State_Chm, State_Diag, State_Grid, &
                              State_Met, EHg0_oc,    RC )

    ! Get land mercury flux of Hg0
    CALL LAND_MERCURY_FLUX( LHgSNOW    = LHgSNOW,                            &
                            State_Grid = State_Grid,                         &
                            State_Met  = State_Met,                          &
                            LFLUX      = EHg0_ln                            )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a LAND_FLUX' )

    ! Get soil mercury flux of Hg0
    CALL SOILEMIS( EHg0_dist  = EHg0_dist,                                   &
                   State_Grid = State_Grid,                                  &
                   State_Met  = State_Met,                                   &
                   EHg0_so    = EHg0_so                                     )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SOILEMIS' )

    ! Get snow mercury flux of Hg0
    CALL SNOWPACK_MERCURY_FLUX( LHgSNOW    = LHgSNOW,                        &
                                State_Chm  = State_Chm,                      &
                                State_Grid = State_Grid,                     &
                                State_Met  = State_Met,                      &
                                FLUX       = EHg0_snow                      )
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SNOW_FLUX' )

    ! If we are using the non-local PBL mixing,
    ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
    ! EMIS_SAVE is now HG_EMIS (ckeller, 10/21/2014)
    IF ( Input_Opt%LNLPBL ) HG_EMIS = 0.0e+0_fp

    ! Zero arrays for Hg deposition
    CALL RESET_HG_DEP_ARRAYS()
    IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: ' // &
         'a RESET_HG_DEP_ARRAYS' )

    ! Add Hg(0) source into State_Chm%Species [kg]
    CALL SRCHg0( Input_Opt,  State_Chm, State_Diag, State_Grid, State_Met, RC )
   IF ( prtDebug ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Routine EMISSMERCURY in mercury_mod.F90')
       RETURN
    ENDIF

  END SUBROUTINE EMISSMERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emithg
!
! !DESCRIPTION: Subroutine EMITHG directs emission either to the chemical
!  species array (State\_Chm%Species) directly or to Hg\_EMIS for use by the
!  non-local PBL mixing. This is a programming convenience.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMITHG( I, J, L, N, E_HG, Input_Opt, State_Chm, State_Grid )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_TS_EMIS
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L, N  ! Grid boxes + species #
    REAL(fp),       INTENT(IN)    :: E_Hg        ! Hg emissions
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  27 Aug 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)               :: AM2, TS

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !=================================================================
    ! EMITHG begins here!
    !=================================================================

    IF ( Input_Opt%LNLPBL ) THEN

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       !
       ! Save emissions for non-local PBL mixing or emit directly.
       ! Make sure that emitted mass is non-negative
       ! This is hear only for consistency with old code which warned
       ! of underflow error (cdh, 08/27/09)
       ! EMIS_SAVE is now HG_EMIS array. Convert kg to kg/m2/s
       ! (ckeller, 09/23/2014)
       !--------------------------------------------------------------

       ! Surface area [m2]
       AM2            = State_Grid%Area_M2(I,J)

       ! Emission timestep
       TS             = GET_TS_EMIS()

       ! Save emissions as [kg/m2/s].  These will be added
       ! to the chemical species array in routine DO_TEND
       ! (in mixing_mod.F90).
       Hg_EMIS(I,J,N) = Hg_EMIS(I,J,N) + ( MAX( E_Hg, 0.0_fp ) / AM2 / TS )

    ELSE

       !--------------------------------------------------------------
       ! We are using FULL PBL MIXING (routine TURBDAY)
       ! so add directly to the State_Chm%Species array
       !--------------------------------------------------------------

       ! Point to the chemical spcies array [kg]
       Spc                => State_Chm%Species

       ! Add emissions
       Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) + MAX( E_Hg, 0.0_fp )

       ! Free pointer
       Spc                => NULL()

    ENDIF

  END SUBROUTINE EMITHG
!EOP
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcHg0
!
! !DESCRIPTION: Subroutine SRCHg0 is the subroutine for Hg(0) emissions.
!  Emissions of Hg(0) will be distributed throughout the boundary layer.
!  (eck, cdh, bmy, 1/21/05, 4/6/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SRCHg0( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : GET_TS_EMIS
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
!  21 Jan 2005 - N. (Eckley) Selin, C. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      J,    L
    REAL(fp) :: DTSRCE, E_Hg, T_Hg

    !=================================================================
    ! SRCHg0 begins here!
    !=================================================================
    
    ! Initialize
    RC        = GC_SUCCESS       
    DTSRCE    = GET_TS_EMIS()        ! Timestep [s]
    
    ! Zero ocean and snow if the "anthro Hg only" option is selected
    IF ( LAnthroHgOnly ) THEN
       EHg0_oc   = 0.0_fp
       EHg0_snow = 0.0_fp
    ENDIF

    ! Loop over grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, T_Hg, E_Hg                                      )&
    !$OMP COLLAPSE( 2                                                       )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Compute total Hg(0) emissions not computed by HEMCO
       ! (oceans + land + natural + snow)
       T_Hg = EHg0_oc(I,J) + EHg0_ln(I,J) + EHg0_so(I,J) + EHg0_snow(I,J)

       !=====================================================================
       ! Partition Hg0 throughout PBL; store into State_Chm%Species [kg]
       ! Now make sure State_Chm%Species does not underflow (cdh, bmy, 4/6/06)
       !=====================================================================
       DO L = 1, State_Met%PBL_MAX_L
          E_Hg = State_Met%F_OF_PBL(I,J,L) * T_Hg * DTSRCE
          CALL EmitHg( I, J, L, id_Hg0, E_Hg, Input_Opt, State_Chm, State_Grid )
       ENDDO

       !=====================================================================
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Save the various Hg0 emissions (land, ocean, snow, soil)
       ! in units of [kg/s].
       !
       ! NOTE: All other emission categories are archived via HEMCO.
       !=====================================================================

       ! Land Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0land ) THEN
          State_Diag%EmisHg0land(I,J) = EHg0_ln(I,J)
       ENDIF

       ! Oceanic Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0ocean ) THEN
          State_Diag%EmisHg0ocean(I,J) = EHg0_oc(I,J)
       ENDIF

       ! Snow Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0snow ) THEN
          State_Diag%EmisHg0snow(I,J) = EHg0_snow(I,J)
       ENDIF

       ! Soil Hg0 emissions [kg/s]
       IF ( State_Diag%Archive_EmisHg0soil ) THEN
          State_Diag%EmisHg0soil(I,J) = EHg0_so(I,J)
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SRCHg0
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemmercury
!
! !DESCRIPTION: Subroutine CHEMMERCURY is the driver routine for mercury
!  chemistry in the GEOS-CHEM module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ChemMercury( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, RC                         )
!
! !USES:
!
    USE Depo_Mercury_Mod,   ONLY : ADD_Hg2_DD
    USE Depo_Mercury_Mod,   ONLY : ADD_HgP_DD
    USE FAST_JX_MOD,        ONLY : FAST_JX
    USE CMN_FJX_MOD
    USE GcKpp_Monitor,      ONLY : SPC_NAMES, FAM_NAMES
    USE GcKpp_Parameters
    USE GcKpp_Integrator,   ONLY : Integrate
    USE GcKpp_Function
    USE GcKpp_Model
    USE Gckpp_Global
    USE GcKpp_Rates,        ONLY : UPDATE_RCONST, RCONST
    USE Timers_Mod
    USE PhysConstants,      ONLY : AVO
    USE State_Chm_Mod,      ONLY : Ind_
    USE PRESSURE_MOD
    USE Species_Mod,        ONLY : Species
    USE Time_Mod,           ONLY : Get_Ts_Chem
    USE Time_Mod,           ONLY : Get_Day
    USE Time_Mod,           ONLY : Get_Month
    USE Time_Mod,           ONLY : Get_Year
    USE Time_Mod,           ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_DAY
    USE Time_Mod,           ONLY : ITS_TIME_FOR_A3
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP, DEBUG_MSG, SAFE_DIV
    USE HCO_STATE_GC_MOD,   ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species, SpcConc
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE          :: FIRST = .TRUE.

    ! Scalars
    LOGICAL                :: prtDebug
    INTEGER                :: I,         J,        L,         K
    INTEGER                :: N,         NN,       CN,        Hg_Cat
    INTEGER                :: NA,        F,        SpcID,     KppID
    INTEGER                :: P,         MONTH,    YEAR,      IRH
    INTEGER                :: TotSteps,  TotFuncs, TotJacob,  TotAccep
    INTEGER                :: TotRejec,  TotNumLU, HCRC,      IERR
    INTEGER                :: Day,       S
    REAL(fp)               :: REL_HUM,   Start,     Finish,   rtim
    REAL(fp)               :: itim,      TOUT,      T,        TIN

    ! Strings
    CHARACTER(LEN=16)      :: thisName
    CHARACTER(LEN=63)      :: origUnit
    CHARACTER(LEN=255)     :: errMsg
    CHARACTER(LEN=255)     :: thisLoc

    ! Arrays
    INTEGER                :: ICNTRL(20)
    INTEGER                :: ISTATUS(20)
    REAL(dp)               :: RCNTRL(20)
    REAL(dp)               :: RSTATE(20)
    REAL(dp)               :: Vloc(NVAR)
    REAL(dp)               :: Aout(NREACT)
#ifdef MODEL_GEOS
    REAL(dp)               :: localC(NSPEC)
#endif

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(fp),      POINTER :: TK(:,:,:   )

    ! Objects
    TYPE(Species), POINTER :: SpcInfo
!
! !DEFINED PARAMETERS:
!
    ! Toggle hetchem or photolysis on/off for testing (default=on)
    LOGICAL,  PARAMETER :: DO_HETCHEM  = .TRUE.
    LOGICAL,  PARAMETER :: DO_PHOTCHEM = .TRUE.

    ! Relative Humidities (to be passed to FAST_JX)
    REAL(fp), PARAMETER :: RH(5) = (/0.0_fp, 0.5_fp, 0.7_fp, 0.8_fp, 0.9_fp/)

    !========================================================================
    ! CHEMMERCURY begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    errMsg   = ''
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    thisLoc  = ' -> at ChemMercury (in GeosCore/mercury_mod.F90)'
    itim     =  0.0_fp            ! For KPP timing
    rtim     =  0.0_fp            ! For KPP timing
    totsteps =  0                 ! Total # of KPP timesteps
    totfuncs =  0                 ! Total # of integrator function calls
    totjacob =  0                 ! Total # of jacobian calls
    totaccep =  0                 ! Total # of KPP calls that finished OK
    totrejec =  0                 ! Total # of KPP calls that didn't converge
    totnumLU =  0                 ! Total # of LU decomposition calls
    Day      =  Get_Day()         ! Current day
    Month    =  Get_Month()       ! Current month
    Year     =  Get_Year()        ! Current year
    Spc      => State_Chm%Species ! Chemical species array [kg]
    TK       => State_Met%T       ! Temperature [K]
    SpcInfo  => NULL()            ! Pointer to GEOS-Chem species database
    Failed2x =  .FALSE.           ! Flag for graceful exit of simulation

    !========================================================================
    ! Set chemistry options and pointers to chemical inputs from HEMCO
    !========================================================================
    IF ( FIRST ) THEN
       IF ( .not. DO_HETCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Heterogeneous chemistry'    // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
       IF ( .not. DO_PHOTCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Photolysis chemistry'       // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
    ENDIF

    !========================================================================
    ! Recompute AOD and related properties when it's a new month
    ! from the data read in via HEMCO
    !========================================================================
    IF ( ITS_A_NEW_MONTH() ) THEN

       ! Get pointers to fields read via HEMCO
       CALL Set_HCOPointers ( Input_Opt, State_Chm, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
           errMsg = 'Error encountered in "Set_HCOPointers"!'
           CALL GC_Error( errMsg, RC, thisLoc )
           RETURN
       ENDIF

       ! Set AOD fields to pass to FastJX
       IRHARR  = 1
       ODAER   = 0.0_fp
       ODMDUST = 0.0_fp

       !$OMP PARALLEL DO                                                    &
       !$OMP DEFAULT( SHARED                                               )&
       !$OMP PRIVATE( I, J, L, N, REL_HUM                                  )&
       !$OMP COLLAPSE( 3                                                   )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Dust OD
          DO N = 1, N_Dust
             ODMDUST(I,J,L,1,N) = AeroPtr(N)%AOD(I,J,L)
          ENDDO

          ! Aerosol OD
          DO N = 1, N_Aer
             ODAER(I,J,L,1,N) = AeroPtr(N_Dust+N)%AOD(I,J,L)
          ENDDO

          ! Save IRHARR
          REL_HUM =  GLOB_RH(I,J,L)
          IF (      REL_HUM <= RH(2) ) THEN
             IRHARR(I,J,L) = 1
          ELSE IF ( REL_HUM <= RH(3) ) THEN
             IRHARR(I,J,L) = 2
          ELSE IF ( REL_HUM <= RH(4) ) THEN
             IRHARR(I,J,L) = 3
          ELSE IF ( REL_HUM <= RH(5) ) THEN
             IRHARR(I,J,L) = 4
          ELSE
             IRHARR(I,J,L) = 5
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF (State_Diag%Archive_Loss           ) State_Diag%Loss           = 0.0_f4
    IF (State_Diag%Archive_Prod           ) State_Diag%Prod           = 0.0_f4
    IF (State_Diag%Archive_SatDiagnLoss   ) State_Diag%SatDiagnLoss   = 0.0_f4
    IF (State_Diag%Archive_SatDiagnProd   ) State_Diag%SatDiagnProd   = 0.0_f4
    IF (State_Diag%Archive_JVal           ) State_Diag%JVal           = 0.0_f4
    IF (State_Diag%Archive_SatDiagnJVal   ) State_Diag%SatDiagnJVal   = 0.0_f4
    IF (State_Diag%Archive_JNoon          ) State_Diag%JNoon          = 0.0_f4
    IF (State_Diag%Archive_OHreactivity   ) State_Diag%OHreactivity   = 0.0_f4
    IF (State_Diag%Archive_RxnRate        ) State_Diag%RxnRate        = 0.0_f4
    IF (State_Diag%Archive_HgBrAfterChem  ) State_Diag%HgBrAfterChem  = 0.0_f4
    IF (State_Diag%Archive_HgClAfterChem  ) State_Diag%HgClAfterChem  = 0.0_f4
    IF (State_Diag%Archive_HgOHAfterChem  ) State_Diag%HgOHAfterChem  = 0.0_f4
    IF (State_Diag%Archive_HgBrOAfterChem ) State_Diag%HgBrOAfterChem = 0.0_f4
    IF (State_Diag%Archive_HgClOAfterChem ) State_Diag%HgClOAfterChem = 0.0_f4
    IF (State_Diag%Archive_HgOHOAfterChem ) State_Diag%HgOHOAfterChem = 0.0_f4
    IF (State_Diag%Archive_KppDiags) THEN
       IF (State_Diag%Archive_KppIntCounts) State_Diag%KppIntCounts   = 0.0_f4
       IF (State_Diag%Archive_KppJacCounts) State_Diag%KppJacCounts   = 0.0_f4
       IF (State_Diag%Archive_KppTotSteps ) State_Diag%KppTotSteps    = 0.0_f4
       IF (State_Diag%Archive_KppAccSteps ) State_Diag%KppAccSteps    = 0.0_f4
       IF (State_Diag%Archive_KppRejSteps ) State_Diag%KppRejSteps    = 0.0_f4
       IF (State_Diag%Archive_KppLuDecomps) State_Diag%KppLuDecomps   = 0.0_f4
       IF (State_Diag%Archive_KppSubsts   ) State_Diag%KppSubsts      = 0.0_f4
       IF (State_Diag%Archive_KppSmDecomps) State_Diag%KppSmDecomps   = 0.0_f4
    ENDIF

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Unit conversion error!'
       CALL GC_Error( errMsg, RC, 'mercury_mod.F90')
       RETURN
    ENDIF

    !========================================================================
    ! Call photolysis routine to compute J-Values
    !========================================================================
    IF ( DO_PHOTCHEM ) THEN

        !Compute J values
        CALL Fast_JX( 0,          Input_Opt,  State_Chm,                     &
                      State_Diag, State_Grid, State_Met, RC                 )

        ! Trap potential errors
        IF ( RC /= GC_SUCCESS ) THEN
           errMsg = 'Error encountered in "FAST_JX"!'
           CALL GC_Error( errMsg, RC, thisLoc )
           RETURN
        ENDIF

        !### Debug
        IF ( prtDebug ) THEN
           CALL DEBUG_MSG( '### ChemMercury: after FAST_JX' )
        ENDIF
    ENDIF

    !========================================================================
    ! Set instantaneous oxidant concentrations (molec cm-3)
    !========================================================================
    CALL Set_HgOxidConc( Input_Opt, State_Chm, State_Grid, State_Met, RC    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Set_HgOxidConc"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### ChemMercury: after Set_HgOxidConc' )
    ENDIF

    !========================================================================
    ! Set up integration convergence conditions and timesteps
    ! (cf. M. J. Evans)
    !========================================================================

    !%%%%% TIMESTEPS %%%%%
    DT        = GET_TS_CHEM() ! [s]
    T         = 0d0
    TIN       = T
    TOUT      = T + DT

    !%%%%% CONVERGENCE CRITERIA %%%%%

    ! Absolute tolerance
    ATOL      = 1e-2_dp

    ! Relative tolerance
    RTOL      = 1e-2_dp

    !%%%%% SOLVER OPTIONS %%%%%

    ! Zero all slots of ICNTRL
    ICNTRL    = 0

    ! 0 - non-autonomous, 1 - autonomous
    ICNTRL(1) = 1

    ! 0 - vector tolerances, 1 - scalars
    ICNTRL(2) = 0

    ! Select Integrator
    ICNTRL(3) = 4 ! Rodas3

    ! 0 - adjoint, 1 - no adjoint
    ICNTRL(7) = 1

    ! Turn off calling Update_SUN, Update_RCONST, Update_PHOTO from within
    ! the integrator.  Rate updates are done before calling KPP.
    !  -- Bob Yantosca (03 May 2022)
    ICNTRL(15) = -1

    !=======================================================================
    ! %%%%% SOLVE CHEMISTRY -- This is the main KPP solver loop %%%%%
    !=======================================================================
100 format('No. of function calls:', i6, /,                                 &
           'No. of jacobian calls:', i6, /,                                 &
           'No. of steps:         ', i6, /,                                 &
           'No. of accepted steps:', i6, /,                                 &
           'No. of rejected steps ', i6, /,                                 &
           '       (except at very beginning)',          /,                 &
           'No. of LU decompositions:             ', i6, /,                 &
           'No. of forward/backward substitutions:', i6, /,                 &
           'No. of singular matrix decompositions:', i6, /,                 &
            /,                                                              &
           'Texit, the time corresponding to the      ',        /,          &
           '       computed Y upon return:            ', f11.4, /,          &
           'Hexit, last accepted step before exit:    ', f11.4, /,          &
           'Hnew, last predicted step (not yet taken):', f11.4 )

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,        J,        L,       N                           )&
    !$OMP PRIVATE( IERR,     RCNTRL,   START,   FINISH, ISTATUS             )&
    !$OMP PRIVATE( RSTATE,   SpcID,    KppID,   F,      P                   )&
    !$OMP PRIVATE( Vloc,     Aout,     NN                                   )&
    !$OMP REDUCTION( +:ITIM                                                 )&
    !$OMP REDUCTION( +:RTIM                                                 )&
    !$OMP REDUCTION( +:TOTSTEPS                                             )&
    !$OMP REDUCTION( +:TOTFUNCS                                             )&
    !$OMP REDUCTION( +:TOTJACOB                                             )&
    !$OMP REDUCTION( +:TOTACCEP                                             )&
    !$OMP REDUCTION( +:TOTREJEC                                             )&
    !$OMP REDUCTION( +:TOTNUMLU                                             )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE ( DYNAMIC,  24                                           )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !====================================================================
       ! For safety sake, initialize certain variables for each grid
       ! box (I,J,L), whether or not chemistry will be done there.
       !====================================================================
       IERR      = 0         ! Success or failure flag
       P         = 0         ! GEOS-Chem photolyis species ID
       ISTATUS   = 0.0_dp    ! Rosenbrock output
       PHOTOL    = 0.0_dp    ! Photolysis array
       RCNTRL    = 0.0_fp    ! Rosenbrock input
       RSTATE    = 0.0_dp    ! Rosenbrock output
       C         = 0.0_dp    ! KPP species conc's
       RCONST    = 0.0_dp    ! KPP rate constants
       CFACTOR   = 1.0_dp    ! KPP conversion factor (not really needed)

       !====================================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !====================================================================

       ! If we are not below the stratopause don't do the chemistry!
       !IF ( L > State_Grid%MaxStratLev ) CYCLE
       IF ( .not. State_Met%InChemGrid( I, J, L ) ) CYCLE

       ! Skipping buffer zone (lzh, 08/10/2014)
       IF ( State_Grid%NestedGrid ) THEN
          IF ( J <=                 State_Grid%SouthBuffer ) CYCLE
          IF ( J >  State_Grid%NY - State_Grid%NorthBuffer ) CYCLE
          IF ( I <=                 State_Grid%EastBuffer  ) CYCLE
          IF ( I >  State_Grid%NX - State_Grid%WestBuffer  ) CYCLE
       ENDIF

       !====================================================================
       ! Get photolysis rates (daytime only)
       !====================================================================
       IF ( State_Met%SUNCOSmid(I,J) > -0.1391731e+0_fp ) THEN

          ! Loop over the FAST-JX photolysis species
          DO N = 1, nRatJ

             ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
             IF ( DO_PHOTCHEM ) THEN
                PHOTOL(N) = ZPJ(L,N,I,J)
             ENDIF

             !---------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Instantaneous photolysis rates [s-1] (aka J-values)
             ! and noontime photolysis rates [s-1]
             !---------------------------------------------------------------
             IF ( State_Diag%Archive_JVal ) THEN

                ! GC photolysis species index
                P = GC_Photo_Id(N)

                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                State_Diag%JVal(I,J,L,P) = State_Diag%JVal(I,J,L,P)          &
                                         + PHOTOL(N)
             ENDIF

             ! J_values for satellite diagnostics
             IF ( State_Diag%Archive_SatDiagnJVal ) THEN

                ! GC photolysis species index
                P = GC_Photo_Id(N)

                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                State_Diag%SatDiagnJVal(I,J,L,P) =                           &
                State_Diag%SatDiagnJVal(I,J,L,P) + PHOTOL(N)
             ENDIF

          ENDDO
       ENDIF

       !====================================================================
       ! Copy values at each gridbox into variables in gckpp_Global.F90
       ! This includes e.g. temperature, air density, and quantities
       ! needed for heterogeneous chemistry
       !====================================================================
       CALL Set_Kpp_GridBox_Values( I          = I,                          &
                                    J          = J,                          &
                                    L          = L,                          &
                                    Input_Opt  = Input_Opt,                  &
                                    State_Chm  = State_Chm,                  &
                                    State_Grid = State_Grid,                 &
                                    State_Met  = State_Met,                  &
                                    RC         = RC                         )

       !=====================================================================
       ! Update KPP rates
       !=====================================================================

       ! Zero out dummy species index in KPP
       DO F = 1, NFAM
          KppID = PL_Kpp_Id(F)
          IF ( KppID > 0 ) C(KppID) = 0.0_dp
       ENDDO

       ! Update the array of rate constants
       CALL Update_RCONST( )

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Archive KPP equation rates (Aout).  For GEOS-Chem in GEOS, also
       ! archive the time derivative of variable species (Vdot).
       !
       ! NOTE: Replace VAR with C(1:NVAR) and FIX with C(NVAR+1:NSPEC),
       ! because VAR and FIX are now local to the integrator
       !  -- Bob Yantosca (03 May 2022)
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_RxnRate ) THEN

#ifdef MODEL_GEOS
          !------------------------------------------
          ! GEOS-Chem in GEOS: Get Aout and Vdotout
          !------------------------------------------
          CALL Fun ( V       = C(1:NVAR),                                    &
                     F       = C(NVAR+1:NSPEC),                              &
                     RCT     = RCONST,                                       &
                     Vdot    = Vloc,                                         &
                     Aout    = Aout )
#else
          !------------------------------------------
          ! All other contexts: Get Aout only
          !------------------------------------------
          CALL Fun ( V       = C(1:NVAR),                                    &
                     F       = C(NVAR+1:NSPEC),                              &
                     RCT     = RCONST,                                       &
                     Vdot    = Vloc,                                         &
                     Aout    = Aout                                         )
#endif

          ! Only save requested equation rates
          DO S = 1, State_Diag%Map_RxnRate%nSlots
             N = State_Diag%Map_RxnRate%slot2Id(S)
             State_Diag%RxnRate(I,J,L,S) = Aout(N)
          ENDDO

       ENDIF

       !=====================================================================
       ! Set options for the KPP Integrator (M. J. Evans)
       !=====================================================================

       ! Zero all slots of RCNTRL
       RCNTRL    = 0.0_fp

       ! Initialize Hstart (the starting value of the integration step
       ! size with the value of Hnew (the last predicted but not yet 
       ! taken timestep)  saved to the the restart file.
       RCNTRL(3) = State_Chm%KPPHvalue(I,J,L)

       !=====================================================================
       ! Integrate the box forwards
       !=====================================================================
       CALL Integrate( TIN,    TOUT,    ICNTRL,                              &
                       RCNTRL, ISTATUS, RSTATE, IERR                        )

       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
       ENDIF

       !--------------------------------------------------------------------
       ! Try another time if it failed
       !--------------------------------------------------------------------
       IF ( IERR < 0 ) THEN

#ifdef MODEL_GEOS
          ! Save a backup copy of C (for GEOS-Chem in GEOS only)
          localC    = C
#endif

          ! Reset first time step and start concentrations
          ! Retry the integration with non-optimized
          ! settings
          RCNTRL(3) = 0e+0_fp
          C         = 0.0_dp

          ! Update rate constants again
          CALL Update_RCONST( )

          ! Integrate again
          CALL Integrate( TIN,    TOUT,    ICNTRL,                           &
                          RCNTRL, ISTATUS, RSTATE, IERR                     )

          !------------------------------------------------------------------
          ! Exit upon the second failure
          !------------------------------------------------------------------
          IF ( IERR < 0 ) THEN
             WRITE(6,     '(a   )' ) '## INTEGRATE FAILED TWICE !!! '
             WRITE(ERRMSG,'(a,i3)' ) 'Integrator error code :', IERR
#ifdef MODEL_GEOS
             IF ( Input_Opt%KppStop ) THEN
                CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
                ! Revert to start values
             ELSE
                C = localC
             ENDIF
             IF ( ASSOCIATED(State_Diag%KppError) ) THEN
                State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
             ENDIF
#else
             ! Set a flag to break out of loop gracefully
             ! NOTE: You can set a GDB breakpoint here to examine the error
             !$OMP CRITICAL
             Failed2x = .TRUE.

             ! Print concentrations at trouble box KPP error
             PRINT*, REPEAT( '###', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Species concentrations at problem box ', I, J, L
             PRINT*, REPEAT( '###', 79 )
             DO N = 1, NSPEC
                PRINT*, '### ', C(N), TRIM( ADJUSTL( SPC_NAMES(N) ) )
             ENDDO

             ! Print rate constants at trouble box KPP error
             PRINT*, REPEAT( '###', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Reaction rates at problem box ', I, J, L
             PRINT*, REPEAT( '###', 79 )
             DO N = 1, NREACT
                PRINT*, RCONST(N), TRIM( ADJUSTL( EQN_NAMES(N) ) )
             ENDDO
             !$OMP END CRITICAL
#endif
          ENDIF
       ENDIF

       !=====================================================================
       ! Continue upon successful return...
       !=====================================================================

       ! Save Hnew (the last predicted but not taken step) from the 3rd slot
       ! of RSTATE into State_Chm so that it can be written to the restart
       ! file.  For simulations that are broken into multiple stages,
       ! Hstart will be initialized to the value of Hnew from the restart
       ! file at startup (see above).
       State_Chm%KPPHvalue(I,J,L) = RSTATE(3)

       !=====================================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species
       !=====================================================================
       ! Loop over KPP species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Skip if this is not a GEOS-Chem species
          IF ( SpcID <= 0 ) CYCLE

          ! Set negative concentrations to zero
          C(N) = MAX( C(N), 0.0_dp )

          ! Copy concentrations back into species concentration array
          State_Chm%Species(SpcID)%Conc(I,J,L) = REAL( C(N), kind=fp )

       ENDDO

       !=====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "VAR" array (in module GcKpp_Global.F90).
       !=====================================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Loss ) THEN
          DO S = 1, State_Diag%Map_Loss%nSlots
             KppId = State_Diag%Map_Loss%slot2Id(S)
             State_Diag%Loss(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Prod ) THEN
          DO S = 1, State_Diag%Map_Prod%nSlots
             KppID = State_Diag%Map_Prod%slot2Id(S)
             State_Diag%Prod(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Satellite diagnostic: chemical loss of species/families [molec/cm3/s]
       IF ( State_Diag%Archive_SatDiagnLoss ) THEN
          DO S = 1, State_Diag%Map_SatDiagnLoss%nSlots
             KppId = State_Diag%Map_SatDiagnLoss%slot2Id(S)
             State_Diag%SatdiagnLoss(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Satellite diagnostic: chemical prod of species/families [molec/cm3/s]
       IF ( State_Diag%Archive_SatDiagnProd ) THEN
          DO S = 1, State_Diag%Map_SatDiagnProd%nSlots
             KppID = State_Diag%Map_SatDiagnProd%slot2Id(S)
             State_Diag%SatDiagnProd(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       !=====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Archive concetration of short-lived radicals [mol/mol]
       !=====================================================================
       IF ( State_Diag%Archive_HgBrAfterChem  ) THEN
          State_Diag%HgBrAfterChem(I,J,L)  = Spc(id_HgBr)%Conc(I,J,L)         &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

       IF ( State_Diag%Archive_HgClAfterChem  ) THEN
          State_Diag%HgClAfterChem(I,J,L)  = Spc(id_HgCl)%Conc(I,J,L)         &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

       IF ( State_Diag%Archive_HgOHAfterChem  ) THEN
          State_Diag%HgOHAfterChem(I,J,L)  = Spc(id_HgOH)%Conc(I,J,L)         &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

       IF ( State_Diag%Archive_HgBrOAfterChem ) THEN
          State_Diag%HgBrOAfterChem(I,J,L) = Spc(id_HgBrO)%Conc(I,J,L)        &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

       IF ( State_Diag%Archive_HgClOAfterChem ) THEN
          State_Diag%HgClOAfterChem(I,J,L) = Spc(id_HgClO)%Conc(I,J,L)        &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

       IF ( State_Diag%Archive_HgOHOAfterChem ) THEN
          State_Diag%HgOHOAfterChem(I,J,L) = Spc(id_HgOHO)%Conc(I,J,L)        &
                                           / State_Met%AirNumDen(I,J,L)
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=======================================================================
    ! Return gracefully if integration failed 2x anywhere
    ! (as we cannot break out of a parallel DO loop!)
    !=======================================================================
    IF ( Failed2x ) THEN
       ErrMsg = 'KPP failed to converge after 2 iterations!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Partition Hg2 between gas and aerosol phase
    !========================================================================
    CALL PartitionHg2( Input_Opt,  State_Chm, State_Diag,                    &
                       State_Grid, State_Met, RC                            )

    !========================================================================
    ! Hg2 uptake by seasalt aerosols in the MBL
    !========================================================================
    CALL SeaSaltUptake( Input_Opt,  State_Chm, State_Diag,                   &
                        State_Grid, State_Met, RC                           )

    !========================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !========================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Unit conversion error!'
       CALL GC_Error( errMsg, RC, 'mercury_mod.F90' )
       RETURN
    ENDIF

    !### Debug output (uncomment if needed)
    !IF ( ITS_A_NEW_DAY() ) THEN
    !    WRITE(*,*) 'Total Hg0 mass [Mg]: ', &
    !                SUM ( State_Chm%Species(id_Hg0)%Conc(:,:,:) )
    !    Hg2Sum = 0.0_fp
    !    !$OMP PARALLEL DO       &
    !    !$OMP DEFAULT( SHARED ) &
    !    !$OMP PRIVATE( N )
    !    DO N = 2, 25
    !       Hg2Sum = Hg2Sum + SUM(State_Chm%Species(N)%Conc(:,:,:))
    !    ENDDO
    !    !$OMP END PARALLEL DO
    !    WRITE(*,*) 'Total Hg2 mass [Mg]: ', Hg2Sum
    !ENDIF

    ! Free pointer memory
    TK      => NULL()
    SpcInfo => NULL()

  END SUBROUTINE ChemMercury
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_HcoPointers
!
! !DESCRIPTION: Subroutine Set_HCOPointers gets the offline chemistry data
! read by HEMCO. The pointers only need to be established once. Target data
! is automatically updated through HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_HcoPointers( Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_State_GC_Mod, ONLY : HcoState
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_ExtList_Mod,  ONLY : HCO_GetOpt
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorological State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (based on set_brypointers)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: is_BrOx_GC
    INTEGER             :: N, spcID

    ! Strings
    CHARACTER(LEN=16)   :: thisOpt
    CHARACTER(LEN=255)  :: fieldName
    CHARACTER(LEN=255)  :: thisLoc
    CHARACTER(LEN=1024) :: errMsg

    !=================================================================
    ! Set_HCOPointers begins here
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Set_HCOPointers (in module GeosCore/mercury_mod.F90)'

    ! Check if the Br_Ox switch is set
    thisOpt    = HCO_GetOpt( HcoState%Config%ExtList, 'BrOx_GC', ExtNr=0 )
    is_BrOx_GC = ( INDEX( thisOpt, 'true' ) > 0 )

    ! Do for each fixed KPP species
    DO N = 1, State_Chm%nKppFix

       ! Get species ID
       spcID = State_Chm%Map_KppFix(N)

       ! Construct field name using species name
       fieldName = 'GLOBAL_' // TRIM( State_Chm%SpcData(SpcID)%Info%Name )

       ! Get pointer to oxidant field [molec/cm3]
       CALL HCO_GetPtr( HcoState, fieldName, FixSpcPtr(N)%Data, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot get pointer from HEMCO! Oxidant data '//&
                   'is expected to be listed in the HEMCO configuration '  //&
                   'file. This error occured when trying to get field '     //&
                  TRIM( FieldName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDDO !N

    !----------------------------------
    ! PM2.5 mass concentration (ug m-3)
    !----------------------------------
    ! Get pointer to this field.
    CALL HCO_GetPtr( HcoState, 'GLOBAL_PM25', GLOB_PM25, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get pointer from HEMCO for Global PM2.5'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !----------------------------------
    ! OA fraction
    !----------------------------------
    ! Get pointer to this field.
    CALL HCO_GetPtr( HcoState, 'GLOBAL_fOA', GLOB_fOA, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get pointer from HEMCO for Global fOA'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !----------------------------------
    ! Aerosol fields
    !----------------------------------

    ! Do for each aerosol species
    DO N = 1, N_Aer + N_Dust

       !------------------------------
       ! AOD
       !------------------------------

       ! Get aerosol species name
       fieldName = 'AOD_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, fieldName, AeroPtr(N)%AOD, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot get pointer from HEMCO for ' // TRIM( fieldName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !------------------------------
       ! Area (cm2 cm-3)
       !------------------------------

       ! Get aerosol species name
       fieldName = 'Area_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, fieldName, AeroPtr(N)%Area, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( fieldName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !------------------------------
       ! Radi (cm)
       !------------------------------

       ! Get aerosol species name
       fieldName = 'Radi_' // TRIM( AerSpcNames(N) )

       ! Get pointer to this field. These are AODs.
       CALL HCO_GetPtr( HcoState, fieldName, AeroPtr(N)%Radi, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot get pointer from HEMCO for ' //&
                  TRIM( FieldName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDDO !N

    !----------------------------------
    ! Relative humidity
    !----------------------------------
    ! Get pointer to this field.
    CALL HCO_GetPtr( HcoState, 'GLOBAL_RH', GLOB_RH, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get pointer from HEMCO for Global RH'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Set_HCOPointers
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_HgOxidConc
!
! !DESCRIPTION: Transfers oxidant concentration fields
!               to State_Chm after applying diurnal variation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_HgOxidConc( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE HCO_State_GC_Mod,   ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH


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
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: I, J, L, N         ! lon, lat, lev, indexes
    INTEGER            :: SpcID
    CHARACTER(LEN=60)  :: Prefix             ! utility string
    CHARACTER(LEN=255) :: thisLoc            ! routine location
    CHARACTER(LEN=255) :: errMsg             ! message

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize pointers
    SpcInfo     => NULL()

    !================================================================
    ! Get_HgOxConc begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    errMsg    = ''
    thisLoc   = ' -> at Get_HgOxConc (in GeosCore/mercury_mod.F90)'

    ! Copy values from Input_Opt
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    !=================================================================
    ! Set instantaneous species concentrations
    !=================================================================

    ! Set species concentration to monthly mean value
    ! Do for each fixed KPP species
    DO N = 1, State_Chm%nKppFix

       ! Get species ID
       SpcID = State_Chm%Map_KppFix(N)

       ! Get value from pointer to monthly mean field
       State_Chm%Species(SpcID)%Conc(:,:,:) = FixSpcPtr(N)%Data(:,:,:)
    ENDDO

    ! Impose diurnal cycle
    ! Compute sum of cosine of the solar zenith angle over a 24 hour day
    CALL OhNO3time( State_Grid )

    ! Calculate instantaneous HOx
    CALL DiurnalHOx( State_Chm, State_Grid, State_Met )

    ! Partition NOx based on NO-NO2-O3 photochemical steady state
    CALL PartNOx( State_Chm, State_Grid, State_Met )

    ! Apply dirunal cycle and partition XOx based on O3, NO and J_XO
    CALL PartXOx( State_Chm, State_Grid, State_Met )

    ! Add BrOx from springtime polar bromine explosion events
    IF ( LPOLARBr ) CALL PolarBrOx( State_Chm, State_Grid, State_Met )

  END SUBROUTINE Set_HgOxidConc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
  SUBROUTINE DiurnalHOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Time_Mod,       ONLY : Get_Ts_Chem
    USE Species_Mod,    ONLY : SpcConc
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from GET_OH)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, J, L
    REAL(fp) :: diurnalFac, dtChem, C_OH, C_HO2

    ! Pointer to Species array
    TYPE(SpcConc),  POINTER :: Spc(:)

    !========================================================================
    ! DirunalHOx begins here!
    !========================================================================

    ! Chemistry timestep [s]
    dtChem = Get_Ts_Chem()

    ! Loop over gridcells
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, diurnalFac, C_OH, C_HO2                         )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       diurnalFac = 0.0_fp
       C_OH       = State_Chm%Species(id_OH )%Conc(I,J,L)
       C_HO2      = State_Chm%Species(id_HO2)%Conc(I,J,L)

       ! If sunlight, compute diurnalFac
       ! If nighttime, skip (diurnalFac = 0)
       IF ( State_Met%SUNCOS(I,J) > 0.0_fp .and.  TCOSZ(I,J) > 0.0_fp ) THEN

          diurnalFac = ( State_Met%SUNCOS(I,J)  / TCOSZ(I,J)               ) &
                     * ( 86400.0_fp             / dtChem                   )

          ! Make sure factor is not negative
          diurnalFac = MAX( diurnalFac, 0.0_fp )

        ENDIF

        ! Scale the HOx species by the diurnal factor
        State_Chm%Species(id_OH )%Conc(I,J,L) = C_OH  * DiurnalFac
        State_Chm%Species(id_HO2)%Conc(I,J,L) = C_HO2 * DiurnalFac

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE DiurnalHOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PolarBrOx
!
! !DESCRIPTION: Calculates BrOx during bromine explosion events.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PolarBrOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Species_Mod,    ONLY : SpcConc
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Time_Mod,       ONLY : Get_Month
    USE Cmn_FJX_Mod,    ONLY : ZPJ
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from cdh's and jaf's
!                            GET_BR)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: IS_MOSTLY_ICE
    INTEGER             :: I,              J,             L,        month
    REAL(fp)            :: BrO_POLAR_PPTV, Br_POLAR_PPTV, O3_POLAR, SWRAD
    REAL(fp)            :: BrO_POLAR_CONC, Br_POLAR_CONC, BrO_CONC, Br_CONC
    REAL(fp)            :: FPBL,           JBrO
!
! !DEFINED PARAMETERS:
!
    ! Assume 5 ppb O3 when BrO present. Pohler et al. (2010) shows 5 ppb O3
    ! can exist for 3 ppt < [BrO] < 40 ppt. (jaf, 11/29/11)
    REAL(fp), PARAMETER :: O3_POLAR_PPBV = 5e+0_fp

    ! Rate coefficient BrO + NO -> Br + NO2, cm3/molec/s
    REAL(fp), PARAMETER :: K_BrO_NO = 2.1e-11_fp

    ! Rate coefficient Br + O3 -> BrO + O2, cm3/molec/s
    REAL(fp), PARAMETER :: K_Br_O3  = 1.2e-12_fp

    ! Concentration of NO, based on 10pptv, molec/cm3
    REAL(fp), PARAMETER :: C_NO     = 2.5e+8_fp

    !=================================================================
    ! PolarBrOx begins here!
    !=================================================================

    ! Initialize
    month = Get_Month()

    ! Loop over gridcells and add BrOx in polar regions
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,              J,             L                         )&
    !$OMP PRIVATE( BrO_POLAR_PPTV, Br_POLAR_PPTV, O3_POLAR, SWRAD           )&
    !$OMP PRIVATE( BrO_POLAR_CONC, Br_POLAR_CONC, BrO_CONC, Br_CONC         )&
    !$OMP PRIVATE( FPBL,           IS_MOSTLY_ICE, JBrO                      )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       IS_MOSTLY_ICE   = .FALSE.
       Br_CONC         = 0.0_fp
       Br_POLAR_CONC   = 0.0_fp
       BrO_POLAR_PPTV  = 0.0_fp
       BrO_CONC        = 0.0_fp
       BrO_POLAR_CONC  = 0.0_fp
       Br_POLAR_PPTV   = 0.0_fp
       FPBL            = 0.0_fp
       JBrO            = 0.0_fp
       O3_POLAR        = 0.0_fp
       SWRAD           = 0.0_fp

       ! Skip if we aren't poleward of 50N/S lat in the wintertime hemisphere
       IF ( ( ( State_Grid%YMid(I,J) >  50.0_fp                  )     .and. &
              ( month >= 2                                       )     .and. &
              ( month <= 6                                       ) )   .or.  &
            ( ( State_Grid%YMid(I,J) < -50.0_fp                  )     .and. &
              ( month >= 8                                       )     .and. &
              ( month <= 12                                      ) ) ) THEN

            !----------------------------------------------------------------
            ! Add Br in the polar PBL
            !----------------------------------------------------------------
            FPBL = State_Met%F_UNDER_PBLTOP(I,J,L)

            !----------------------------------------------------------------
            ! Bromine in the polar boundary layer requires the following
            ! criteria be met:
            ! - In the PBL
            ! - Downward shortwave radiation > 100 W/m2 (Pohler et al. 2010)
            ! - Sea ice exists (>50% native boxes have >10% ice cover)
            ! - Breaks in sea ice exist (<100% native boxes have >90% ice cover)
            ! - Month is between Feb & June (Arctic) or Aug & Dec (Antarctic)
            !   based on http://bro.aeronomie.be/level3_monthly.php?cmd=map
            ! - Temperature is less than 0C
            !
            ! If these criteria are met, BrO is a function of ambient temp.
            ! with [BrO] based on findings from Pohler et al. (2010) and
            ! Prados-Roman et al. (2011). O3 used to convert BrO to Br is 5
            ! ppb, based on data from Pohler et al. (2010).
            !----------------------------------------------------------------
            SWRAD         = State_Met%SWGDN(I,J)
            IS_MOSTLY_ICE = ( State_Met%SEAICE00(I,J) <= 0.50_fp       .and. &
                              State_Met%SEAICE90(I,J) <  1.0_fp    )

            IF ( ( FPBL              >  0.0_fp                     )   .and. &
                 ( IS_MOSTLY_ICE                                   )   .and. &
                 ( SWRAD             >  1.0e+2_fp                  )   .and. &
                 ( State_Met%TS(I,J) <= 273.0_fp                   ) ) THEN

              ! Get BrOx concentration from species array
              BrO_CONC = State_Chm%Species(id_BrO)%Conc(I,J,L)
              Br_CONC  = State_Chm%Species(id_Br )%Conc(I,J,L)

              ! Get JBrO
              JBrO     = ZPJ(L,id_phot_BrO,I,J)

              ! [BrO] is a linear function of temperature derived based on
              ! results from Pohler et al. (2010), Prados-Roman et al. (2011)
              ! and ability to match Hg0 seasonal cycle at Alert. (jaf,
              ! 12/24/11)
              IF ( State_Met%TS(I,J) <= 253.0_fp ) THEN
                 BrO_POLAR_PPTV = 20.0_fp
              ELSE
                 BrO_POLAR_PPTV = -1.0_fp                                    &
                                *  ( State_Met%TS(I,J) - 253.0_fp )          &
                                +  20.0_fp
              ENDIF

              ! Convert O3 to molec/cm3
              O3_POLAR       = O3_POLAR_PPBV                                 &
                             * 1.0e-9_fp                                     &
                             * State_Met%AIRNUMDEN(I,J,L)

              ! Polar BrO concentration [pptv]
              BrO_POLAR_PPTV = BrO_POLAR_PPTV                                &
                             * FPBL

              ! Polar Br concentration [pptv]
              Br_POLAR_PPTV  = BrO_POLAR_PPTV                                &
                             * ( JBrO + K_BrO_NO * C_NO  )                   &
                             / ( K_Br_O3 * O3_POLAR      )

              ! Convert BrO conc from [pptv] to molec/cm3
              BrO_POLAR_CONC = BrO_POLAR_PPTV                                &
                             * 1.0e-12_fp                                    &
                             * State_Met%AIRNUMDEN(I,J,L)

              ! Convert Br conc from [pptv] to [molec/cm3]
              Br_POLAR_CONC  = Br_POLAR_PPTV                                 &
                             * 1.0e-12_fp                                    &
                             * State_Met%AIRNUMDEN(I,J,L)

              ! Replace concentrations in species array
              State_Chm%Species(id_BrO)%Conc(I,J,L) = BrO_CONC + BrO_POLAR_CONC
              State_Chm%Species(id_Br )%Conc(I,J,L) = Br_CONC  + Br_POLAR_CONC
              State_Chm%Species(id_O3 )%Conc(I,J,L) = O3_POLAR

            ENDIF
        ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE PolarBrOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PartXOx
!
! !DESCRIPTION: Subroutine PartXOx partitions halogen radicals based on
!               photochemical steady state with ozone and NO
!\\
!\\
! !INTERFACE:
  !
  SUBROUTINE PartXOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Error_Mod,      ONLY : Safe_Div
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Cmn_FJX_Mod,    ONLY : ZPJ
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from hmh's GET_NO2)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,        J,         L
    REAL(fp) :: C_Br,     C_BrO,     C_O3,     C_BrOx,    C_NO
    REAL(fp) :: C_Cl,     C_ClO,     C_ClOx,   k_Br_O3,   k_BrO_NO
    REAL(fp) :: F_Br_BrO, J_BrO,     k_Cl_O3,  k_ClO_NO,  F_Cl_ClO
    REAL(fp) :: J_ClO,    diurnalFac
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: A_BrO_NO     =    8.8e-12_fp
    REAL(fp), PARAMETER :: EdivR_BrO_NO = -260.0_fp
    REAL(fp), PARAMETER :: A_ClO_NO     =    6.4e-12_fp
    REAL(fp), PARAMETER :: EdivR_ClO_NO = -290.0_fp
    REAL(fp), PARAMETER :: A_Br_O3      =    1.6e-11_fp
    REAL(fp), PARAMETER :: EdivR_Br_O3  =  780.0_fp
    REAL(fp), PARAMETER :: A_Cl_O3      =    2.3e-11_fp
    REAL(fp), PARAMETER :: EdivR_Cl_O3  =  200.0_fp

    !========================================================================
    ! PartXOx begins here!
    !
    ! NOTES:
    ! (R1) XO + hv -> X + O3,  j
    ! (R2) XO + NO -> X + NO2, k_XO_NO
    ! (R2) O3 + X ->  XO + O2, k_X_O3
    !
    ! XOx steady state:
    ! [X]/[XO] = (j+k_XO_NO * C_NO) / (k_X_O3 * C_O3)
    !========================================================================

    ! Loop over gridcells
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,         J,          L,        C_Br                    )&
    !$OMP PRIVATE( C_BrO,     C_BrOx,     C_Cl,     C_ClO                   )&
    !$OMP PRIVATE( C_ClOx,    C_O3,       C_NO,     diurnalFac              )&
    !$OMP PRIVATE( F_Br_BrO,  F_Cl_ClO,   J_BrO,    J_ClO                   )&
    !$OMP PRIVATE( k_Br_O3,   k_BrO_NO ,  k_Cl_O3,  k_ClO_NO                )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !---------------------------------------------------------------------
       ! Initialize loop variables
       !---------------------------------------------------------------------
       C_Br       = 0.0_fp
       C_BrO      = 0.0_fp
       C_BrOx     = 0.0_fp
       C_Cl       = 0.0_fp
       C_ClO      = 0.0_fp
       C_ClOx     = 0.0_fp
       C_O3       = 0.0_fp
       C_NO       = 0.0_fp
       diurnalFac = 0.0_fp
       F_Br_BrO   = 0.0_fp
       F_Cl_ClO   = 0.0_fp
       J_BrO      = 0.0_fp
       J_ClO      = 0.0_fp
       k_Br_O3    = 0.0_fp
       k_BrO_NO   = 0.0_fp
       k_Cl_O3    = 0.0_fp
       k_ClO_NO   = 0.0_fp

       !---------------------------------------------------------------------
       ! Only do the following if there is sunlight at this box
       !---------------------------------------------------------------------
       IF ( ( State_Met%SUNCOS(I,J) > 0.0_fp                       )   .and. &
            ( TTDAY(I,J)            > 0.0_fp                       ) ) THEN

          ! Use a constant function for XOx
          diurnalFac = Safe_Div( 1440.0_fp, TTDAY(I,J), 0.0_fp )

          ! Species concentrations
          C_Br       = State_Chm%Species(id_Br )%Conc(I,J,L)
          C_BrO      = State_Chm%Species(id_BrO)%Conc(I,J,L)
          C_Cl       = State_Chm%Species(id_Cl )%Conc(I,J,L)
          C_ClO      = State_Chm%Species(id_ClO)%Conc(I,J,L)
          C_NO       = State_Chm%Species(id_NO )%Conc(I,J,L)
          C_O3       = State_Chm%Species(id_O3 )%Conc(I,J,L)

          ! BrOx and ClOx concentrations
          C_BrOx     = ( C_Br + C_BrO ) * diurnalFac
          C_ClOx     = ( C_Cl + C_ClO ) * diurnalFac

          ! Calculate temperature dependent reaction rates for partitioning
          k_BrO_NO   = A_BrO_NO * EXP( -EdivR_BrO_NO / State_Met%T(I,J,L) )
          k_ClO_NO   = A_ClO_NO * EXP( -EdivR_ClO_NO / State_Met%T(I,J,L) )
          k_Br_O3    = A_Br_O3  * EXP( -EdivR_Br_O3  / State_Met%T(I,J,L) )
          k_Cl_O3    = A_Cl_O3  * EXP( -EdivR_Cl_O3  / State_Met%T(I,J,L) )

          ! Instantaneous J-values [1/s]
          J_BrO      = ZPJ(L,id_phot_BrO,I,J)
          J_ClO      = ZPJ(L,id_phot_ClO,I,J)

          ! Fraction of [X]/[XO]
          F_Br_BrO   = Safe_Div( J_BrO+(k_BrO_NO*C_NO), k_Br_O3*C_O3, 0.0_fp )
          F_Cl_ClO   = Safe_Div( J_ClO+(k_ClO_NO*C_NO), k_Cl_O3*C_O3, 0.0_fp )

          ! Species concentrations, adjusted
          C_Br       =  C_BrOx * F_Br_BrO / ( 1.0_fp + F_Br_BrO )
          C_BrO      =  C_BrOx - C_Br
          C_Cl       =  C_ClOx * F_Cl_ClO / ( 1.0_fp + F_Cl_ClO )
          C_ClO      =  C_ClOx - C_Cl

       ENDIF

       !---------------------------------------------------------------------
       ! Store back into the species array
       ! If nighttime, these will have already been set to zero
       !---------------------------------------------------------------------
       State_Chm%Species(id_Br )%Conc(I,J,L) = C_Br
       State_Chm%Species(id_Cl )%Conc(I,J,L) = C_Cl
       State_Chm%Species(id_BrO)%Conc(I,J,L) = C_BrO
       State_Chm%Species(id_ClO)%Conc(I,J,L) = C_ClO

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE PartXOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PartNOx
!
! !DESCRIPTION: Subroutine PartNOx partitions NOx based on PSS with ozone
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PartNOx( State_Chm, State_Grid, State_Met )
!
! !USES:
!
    USE Error_MoD,      ONLY : Safe_Div
    Use Species_Mod,    ONLY : SpcConc
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Cmn_Fjx_Mod,    ONLY : ZPJ
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  15 Jul 2020 - V. Shah   - Initial version (modified from hmh's GET_NO2)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I,    J,     L
    REAL(fp)            :: C_O3, C_NOx, F_NO2, J_NO2, k3

    ! Pointers
    TYPE(SpcConc), POINTER   :: Spc(:)
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: A     = 3.0e-12_fp
    REAL(fp), PARAMETER :: EdivR = 1.5e+3_fp

    !========================================================================
    ! PartNOx begins here!
    !
    ! NOTES:
    ! (R1) NO2 + hv -> NO  + O  , j
    ! (R2) O   + O2 -> O3       , k2
    ! (R3) O3  + NO -> NO2 + O2 , k3
    !
    ! [NOx] = [NO] + [NO2]
    !
    ! NOx steady state:
    ! j[NO2] = k3[O3][NO]
    ! j[NO2] = k3[O3]([NOx]-[NO2])
    ! [NO2](j+k3[O3] = k3[O3][NOx]
    ! [NO2]/[NOx] = k3[O3]/(j+k3[O3])
    !
    !  k3 = A exp(-E / RT) Arrhenius Equation
    !  A = 3.0e-12 Seinfeld & Pandis
    !  E/R = 1500
    !========================================================================

    ! Point to the chemical spcies array
    Spc => State_Chm%Species

    ! Loop over gridcells
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, C_NOx, C_O3, F_NO2, J_NO2, k3                   )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Zero/initialize PRIVATE loop variables
       C_NOx = Spc(id_NO)%Conc(I,J,L) + Spc(id_NO2)%Conc(I,J,L) ! molec/cm3
       C_O3  = Spc(id_O3)%Conc(I,J,L)                           ! molec/cm3
       J_NO2 = 0.0_fp                                           ! 1/s
       k3    = 0.0_fp                                           ! 1/s
       F_NO2 = 1.0_fp                                           ! unitless

       ! Test for sunlight...
       IF ( State_Met%SUNCOS(I,J) > 0.0_fp .and. TTDAY(I,J) > 0.0_fp ) THEN

          ! Reaction rate
          k3    = A * EXP( -EdivR / State_Met%T(I,J,L) )

          ! Instantaneous JNO2
          J_NO2 = ZPJ(L,id_phot_NO2,I,J)

          ! Fraction of NO2/NOx
          F_NO2 = SAFE_DIV( k3*C_O3, J_NO2+k3*C_O3, 0.0_fp )

       ENDIF

       ! Partition NOx into NO and NO2
       Spc(id_NO2)%Conc(I,J,L) = F_NO2 * C_NOx
       Spc(id_NO )%Conc(I,J,L) = ( 1.0_fp - F_NO2 ) * C_NOx

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE PartNOx
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the
!  solar zenith angle over a 24 hour day, as well as the total length of
!  daylight.  This is needed to scale the offline OH and NO3 concentrations.
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
!  16 Dec 2002 - R. Park & R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER       :: I, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(State_Grid%NX,State_Grid%NY)

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
          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR                          )&
          !$OMP COLLAPSE( 2                                                 )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R( I, J )

             TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMid(I,J)/15.0

             DO WHILE (TIMLOC < 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC > 24.0)
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
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0.0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() ) / 60e+0_fp
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
! !IROUTINE: calc_hg2_seasalt_lossrate
!
! !DESCRIPTION: Calculates the loss rate of RGM [s-1] by uptake of RGM into
!  sea salt aerosol for each model grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Hg2_SeaSalt_LossRate( State_Grid, State_Met )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  The formula used here is a least-squares fit to the full-physics model of
!  sea-salt aerosol emissions, hydroscopic growth, mass-transport limited
!  uptake of Hg(II), and aerosol deposition presented by Holmes et al. (2009)
!  See Holmes et al. 2010 for evaluation of this parameterization.
!  (cdh, 11/25/09)
!
! !REVISION HISTORY:
!  25 Nov 2009 - C. Holmes   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,         J
    REAL(fp) :: LOSS_FREQ, S, SFCWINDSQR, U10M

    !========================================================================
    ! Hg2_SeaSalt_LossRate begins here!
    !========================================================================

    ! Loop over gridboxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, U10M, S, LOSS_FREQ, SFCWINDSQR                     )&
    !$OMP COLLAPSE( 2                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       LOSS_FREQ  = 0.0_fp
       S          = 0.0_fp
       SFCWINDSQR = 0.0_fp
       U10M       = 0.0_fp

       ! Only calculate deposition via sea salt over water
       IF ( State_Met%IsWater(I,J) ) THEN

          ! Wind speed at 10m altitude [m/s]
          ! Don't allow wind >20 m/s (limit of this parameterization)
          SFCWINDSQR = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2
          U10M       = SQRT( SFCWINDSQR )
          U10M       = MAX( MIN( U10M, 20.0_fp ), 1.0_fp )

          ! Relative humidity as a saturation ratio
          ! Use the relative humidity of the lowest layer, although this is
          ! lower than the higher layers of the MBL
          !
          ! Don't allow supersaturation, as [Cl-] is undefined for RH>=1
          ! Cap RH at 99%, Don't allow RH < 75% as happens in coastal areas
          S = MAX( MIN( State_Met%RH(I,J,1), 99.0_fp ), 75.0_fp ) * 1.0e-2_fp

          ! Seasalt loss requency [1/s]
          LOSS_FREQ = 1.0e-10_fp                                             &
                    * ( 1.0_fp - EXP( -57.758_fp * ( 1.0_fp - S ) ) )        &
                    * EXP( ( -1.9351_fp  * U10M                   )          &
                         + (  9.0047_fp  * SQRT( U10M )           )          &
                         + (  0.14788_fp * U10M**1.5_fp           ) )

          ! Loss frequency must be positive
          LOSS_FREQ = MAX( LOSS_FREQ, 1e-10_fp )
       ENDIF

       ! Save loss frequency to an array
       ! If not over water, LOSS_FREQ is already set to zero
       Hg2_SEASALT_LOSSRATE(I,J) = LOSS_FREQ
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Calc_Hg2_SeaSalt_LossRate
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SeasaltUptake
!
! !DESCRIPTION: Subroutine SeasaltUptake does the uptake of Hg2 by seasalt
! aerosols in the MBL.
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SeaSaltUptake( Input_Opt,  State_Chm, State_Diag,               &
                            State_Grid, State_Met, RC                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE Species_Mod,     ONLY : Species
    USE Time_Mod,        ONLY : Get_Ts_Chem
    USE Time_Mod,        ONLY : Its_Time_For_A3
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
!  10-Dec-2021 - V. Shah     - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I,       J,        L,  N,     S
    REAL(fp)           :: gasConc, dGasConc, dt, f_PBL, k

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    ! Boolean
    LOGICAL            :: prtDebug

    !========================================================================
    ! SeaSaltUptake begins here!
    !========================================================================

    ! Assume success
    RC        =  GC_SUCCESS                                ! Success?
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) ! Debug prints?
    DT        = Get_Ts_Chem()                              ! Chem timestep [s]
    errMsg    = ''
    thisLoc   = ' -> at SeaSaltUptake (in GeosCore/mercury_mod.F90)'

    ! Initialize diagnostic arrays
    IF ( State_Diag%Archive_Hg2GasToSSA ) State_Diag%Hg2GasToSSA = 0.0_f4

    ! Calculate loss rate by seasalt uptake
    IF ( ITS_TIME_FOR_A3() ) THEN
       CALL Calc_Hg2_SeaSalt_LossRate( State_Grid, State_Met )
    ENDIF

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N, S, k, f_PBL, gasConc, dGasConc               )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       k        = 0.0_fp
       f_PBL    = 0.0_fp
       gasConc  = 0.0_fp
       dGasConc = 0.0_fp

       ! Proceed only if gridbox below the PBL top
       IF ( L > State_Met%PBL_TOP_L(I,J) ) CYCLE

       ! RGM uptake on sea-salt aerosol, 1/s
       ! (based on SSA production rate (wind speed)
       ! and salinity (RH)) already zero over land
       k = Hg2_SeaSalt_LossRate(I,J)

       ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
       f_PBL = State_Met%F_Under_PblTop(I,J,L)

       ! Do seasalt uptake only for the fraction of the box in the PBL
       IF ( f_PBL > 0.1_fp ) k = f_PBL * k

       ! Calculate loss of each Hg2(g) species
       DO N = 1, nHg2gasSpc

          ! Index for State_Chm%Species array
          S = Map_Hg2gas(N)

          ! Initial Hg(II) gas concentration [molec/cm3]
          gasConc = State_Chm%Species(S)%Conc(I,J,L)

          ! Remove Hg2 lost to sea salt aerosol from the gas phase [molec/cm3]
          dGasConc = gasConc * ( 1.0_fp - EXP( -k * dt ) )
          gasConc  = gasConc - dGasConc

          ! Final Hg2 gas concentration [molec/cm3]
          State_Chm%Species(S)%Conc(I,J,L) = gasConc

          !------------------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Archive Hg2 lost from gas to sea salt aerosol [molec/cm3/s]
          !------------------------------------------------------------------
          IF ( State_Diag%Archive_Hg2GasToSSA ) THEN
             State_Diag%Hg2GasToSSA(I,J,L) =                                 &
             State_Diag%Hg2GasToSSA(I,J,L) + ( dGasConc / DT )
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SeaSaltUptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: partitionhg2
!
! !DESCRIPTION: Subroutine PartitionHg2 calculates the uptake, speciation,
!               and volatilization of Hg2 gas in aerosols.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PartitionHg2( Input_Opt,  State_Chm, State_Diag,                &
                           State_Grid, State_Met, RC                        )
!
! !USES:
!
    USE ErrCode_Mod
    USE gckpp_Global,     ONLY : NUMDEN, SR_TEMP, TEMP
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE Species_Mod,      ONLY : Species, SpcConc
    USE TIME_MOD,         ONLY : GET_TS_CHEM
    USE rateLawUtilFuncs, ONLY : Ars_L1K
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
!  01-Oct-2020 - V. Shah     - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! !LOCAL VARIABLES:
    !
    ! Scalars
    INTEGER             :: I,         J,             L
    INTEGER             :: N,         S
    REAL(fp)            :: aerConc,   aerConcInorg,  aerConcOrg
    REAL(fp)            :: dGasConc,  dt,            Fgas
    REAL(fp)            :: fracOA,    gasConc
    REAL(fp)            :: k,         Kp,            gasAerTot
    REAL(fp)            :: gasTot,    gasTot_eq,     pm25

    ! Arrays
    REAL(fp)            :: xArea(N_Dust+N_Aer)
    REAL(fp)            :: xRadi(N_Dust+N_Aer)
    REAL(fp)            :: xVol (N_Dust+N_Aer)

    ! Strings
    CHARACTER(LEN=255)  :: errMsg
    CHARACTER(LEN=255)  :: thisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
!
! !DEFINED PARAMETERS:
!
    ! Indices for aerosol type
    INTEGER,  PARAMETER :: DU1       = 1  ! Dust (Reff = 0.151 um)
    INTEGER,  PARAMETER :: DU2       = 2  ! Dust (Reff = 0.253 um)
    INTEGER,  PARAMETER :: DU3       = 3  ! Dust (Reff = 0.402 um)
    INTEGER,  PARAMETER :: DU4       = 4  ! Dust (Reff = 0.818 um)
    INTEGER,  PARAMETER :: SUL       = 8  ! Tropospheric Sulfate
    INTEGER,  PARAMETER :: BKC       = 9  ! Black Carbon
    INTEGER,  PARAMETER :: ORC       = 10 ! Organic Carbon
    INTEGER,  PARAMETER :: SSA       = 11 ! Accum-mode sea salt
    INTEGER,  PARAMETER :: SLA       = 13 ! Strat sulfate liq aer

    ! Hg(II) mass accommodation coefficient
    REAL(fp), PARAMETER :: ALPHA_Hg2 = 0.1_fp

    !========================================================================
    ! PARTITIONHg2 begins here!
    !=======================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    DT      =  Get_Ts_Chem()         ! Chemistry timestep [s]
    Spc     => State_Chm%Species     ! Point to species array [molec/cm3]
    errMsg  =  ''
    thisLoc =  ' -> at PartitionHg2 (in module GeosCore/mercury_mod.F90)'

    ! Initialize diagnostic arrays
    IF (State_Diag%Archive_Hg2GToHg2P     ) State_Diag%Hg2GToHg2P      = 0.0_f4
    IF (State_Diag%Archive_Hg2PToHg2G     ) State_Diag%Hg2PToHg2G      = 0.0_f4
    IF (State_Diag%Archive_Hg2GasToHg2StrP) State_Diag%Hg2GasToHg2StrP = 0.0_f4

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,         J,          L,             N                  )&
    !$OMP PRIVATE( S,         aerConc,    aerConcInorg,  aerConcOrg         )&
    !$OMP PRIVATE( dGasConc,  Fgas,       fracOA,        gasAerTot          )&
    !$OMP PRIVATE( gasConc,   gasTot,     gasTot_eq,     k                  )&
    !$OMP PRIVATE( Kp,        pm25,       xArea,         xRadi              )&
    !$OMP PRIVATE( xVol                                                     )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Zero/initialize loop variables
       aerConc      = 0.0_fp
       aerConcInorg = 0.0_fp
       aerConcOrg   = 0.0_fp
       dGasConc     = 0.0_fp
       Fgas         = 0.0_fp
       fracOA       = 0.0_fp
       gasAerTot    = 0.0_fp
       gasConc      = 0.0_fp
       gasTot       = 0.0_fp
       gasTot_eq    = 0.0_fp
       k            = 0.0_fp
       Kp           = 0.0_fp
       pm25         = 0.0_fp
       xArea        = 0.0_fp
       xRadi        = 0.0_fp
       xVol         = 0.0_fp

       !======================================================================
       ! Proceed only if gridbox is below the stratopause
       !======================================================================
       IF ( L > State_Grid%MaxStratLev ) CYCLE

       !----------------------------------------------------------------------
       ! Copy values into THREADPRIVATE variables in gckpp_Global.F90
       ! This is needed in order for Ars_L1K to have the proper inputs
       !---------------------------------------------------------------------
       NUMDEN  = State_Met%AIRNUMDEN(I,J,L)   ! Air density [molec/cm3]
       TEMP    = State_Met%T(I,J,L)           ! Temperature [K]
       SR_TEMP = SQRT( TEMP )                 ! Square root of temperaure

       !---------------------------------------------------------------------
       ! Get aerosol physical properties
       ! xArea = Aerosol specific surface area [cm2/cm3 air]
       ! xRadi = Aerosol effective radius      [cm         ]
       ! xVol  = Aerosol specific volume       [cm3/cm3 air]
       !---------------------------------------------------------------------
       DO N = 1, N_Dust + N_Aer
          xArea(N) = AeroPtr(N)%Area(I,J,L)
          xRadi(N) = AeroPtr(N)%Radi(I,J,L)
          xVol(N)  = xArea(N) * xRadi(N) / 3.0_fp
       ENDDO

       IF ( State_Met%InTroposphere(I,J,L) ) THEN

          !==================================================================
          ! IN THE TROPOSPHERE:
          !
          ! Perform gas-particle partitioning on fine mode aerosols.
          ! Begin by calculating equilibrium concentrations
          ! following Amos et al. (2012)
          !==================================================================

          ! Get PM2.5 concentrations [ug m-3] (skip if too small)
          pm25 = GLOB_PM25(I,J,L)
          IF ( pm25 < 1.0e-3_fp ) CYCLE

          ! Calculate partitioning coefficient (m-3/ug)
          ! This is from Amos et al. (2012)
          Kp = 10.0_fp**( ( 2.5e+3_fp / TEMP ) - 10.0_fp )

          ! Gas fraction
          Fgas = 1.0_fp / ( 1.0_fp + ( Kp * pm25 ) )

          ! Initial Hg2 gas concentration [molec/cm3]
          ! NOTE: gasTot is zeroed at the top of the loop
          DO N = 1, nHg2gasSpc
             S      = Map_Hg2gas(N)
             gasTot = gasTot + Spc(S)%Conc(I,J,L)
          ENDDO

          ! Concentration of Hg2 on aerosols (inorganic + organic),
          ! include any Hg2+ transported from stratosphere
          aerConcInorg = Spc(id_Hg2ClP)%Conc(I,J,L) + Spc(id_Hg2STRP)%Conc(I,J,L)
          aerConcOrg   = Spc(id_Hg2ORGP)%Conc(I,J,L)

          ! Zero stratospheic Hg2
          Spc(id_Hg2STRP)%Conc(I,J,L) = 0.0_fp

          ! Total HgP concentration [molec/cm3]
          aerConc  =  aerConcInorg + aerConcOrg

          ! Add particle-bound species [molec/cm3]
          gasAerTot = gasTot + aerConc

          ! Total Hg2Gas at equilibrium [molec/cm3]
          gasTot_eq = gasAerTot * Fgas

          !-----------------------------------------------------------------
          ! Mass transfer from gas to particles
          !-----------------------------------------------------------------

          ! Loop over all Hg2 gas spcies
          DO N = 1, nHg2gasSpc

             ! Index for State_Chm%Species
             S = Map_Hg2gas(N)

             ! Initial gas concentration [molec/cm3]
             gasConc = Spc(S)%Conc(I,J,L)

             ! Mass transfer rate [1/s] onto dust, sulfate, BC/OC and fine SS
             k = 0.0_fp
             k = k + Ars_L1k( xArea(DU1), xRadi(DU1), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(DU2), xRadi(DU2), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(DU3), xRadi(DU3), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(DU4), xRadi(DU4), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(SUL), xRadi(SUL), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(BKC), xRadi(BKC), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(ORC), xRadi(ORC), ALPHA_Hg2, srMw(N) )
             k = k + Ars_L1k( xArea(SSA), xRadi(SSA), ALPHA_Hg2, srMw(N) )

             ! Amount of mass [molec/cm3] transferred from gas to aerosol
             dGasConc = gasConc * ( 1.0_fp - EXP( -k * DT ) )

             ! Remove transferred mass from gas and add to aerosol
             Spc(S)%Conc(I,J,L) = gasConc - dGasConc
             aerConc            = aerConc + dGasConc

             !---------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Hg2 mass transferred from gas to aerosol [molec/cm3/s]
             !---------------------------------------------------------------
             IF ( State_Diag%Archive_Hg2GToHg2P ) THEN
                State_Diag%Hg2GToHg2P(I,J,L) =                      &
                State_Diag%Hg2GToHg2P(I,J,L) + ( dGasConc / DT )
             ENDIF
          ENDDO

          !------------------------------------------------------------------
          ! Mass transfer from particle to gas
          !------------------------------------------------------------------

          ! Mass transfer rate [1/s] from dust, sulfate, BC/OC and fine SS
          k = 0.0_fp
          k = k + Ars_L1k( xArea(DU1), xRadi(DU1), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(DU2), xRadi(DU2), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(DU3), xRadi(DU3), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(DU4), xRadi(DU4), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(SUL), xRadi(SUL), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(BKC), xRadi(BKC), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(ORC), xRadi(ORC), ALPHA_Hg2, srMw_HgCl2 )
          k = k + Ars_L1k( xArea(SSA), xRadi(SSA), ALPHA_Hg2, srMw_HgCl2 )

          ! Amount of mass [molec/cm3] transferred from aerosol to gas
          ! Limit transferred mass to the amount of HgP present
          dGasConc = gasTot_eq * ( 1.0_fp - EXP( -k * DT ) )
          dGasConc = MIN( dGasConc, aerConc )

          ! Remove transferred mass from aerosol HgCl2 and add to gaseous HgCl2
          aerConc                   = aerConc - dGasConc
          Spc(id_HgCl2)%Conc(I,J,L) = Spc(id_HgCl2)%Conc(I,J,L) + dGasConc

          !------------------------------------------------------------------
          ! Partition aerosol concentration between org and inorg
          !------------------------------------------------------------------

          ! Fraction of organic aerosol in the grid box [unitless]
          fracOA = MIN( GLOB_fOA(I,J,L), 1.0_fp )

          ! Organic and inorganic HgIIP [molec/cm3]
          Spc(id_Hg2OrgP)%Conc(I,J,L) = aerConc * fracOA                ! Org
          Spc(id_Hg2ClP)%Conc(I,J,L)  = aerConc * ( 1.0_fp - fracOA )   ! Inorg

          !------------------------------------------------------------------
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Hg2 mass transferred from aerosol to gas [molec/cm3/s]
          !------------------------------------------------------------------
          IF ( State_Diag%Archive_Hg2PToHg2G ) THEN
             State_Diag%Hg2PToHg2G(I,J,L) =                                  &
             State_Diag%Hg2PToHg2G(I,J,L) + ( dGasConc / DT )
          ENDIF

       ELSE

          !==================================================================
          ! IN THE STRATOSPHERE:
          !
          ! Calculate heterogeneous uptake on stratospheric aqueous aerosols
          !==================================================================

          ! Concentration of Hg2 on aerosols [molec/cm3]
          aerConc = Spc(id_Hg2STRP)%Conc(I,J,L)                              &
                  + Spc(id_Hg2ORGP)%Conc(I,J,L)                              &
                  + Spc(id_Hg2ClP)%Conc(I,J,L)

          ! Zero organic Hg2 aerosol and Hg2Cl aerosol in stratosphere
          Spc(id_Hg2ORGP)%Conc(I,J,L) = 0.0_fp
          Spc(id_Hg2ClP)%Conc(I,J,L)  = 0.0_fp

          !------------------------------------------------------------------
          ! Perform mass transfer between gas and stratopsheric aerosol
          !------------------------------------------------------------------
          DO N = 1, nHg2gasSpc

             ! Index for State_Chm%Species
             S  = Map_Hg2gas(N)

             ! Mass transfer rate [s-1]
             k = Ars_L1K( xArea(SLA), xRadi(SLA), ALPHA_Hg2, srMw(N) )

             ! Initial Hg(II) gas [molec/cm3]
             gasConc = Spc(S)%Conc(I,J,L)

             ! Amount of mass [molec/cm3] transferred from gas to aerosol
             dGasConc = gasConc * ( 1.0_fp - EXP( -k * DT ) )

             ! Remove transferred mass from gas and add to aerosol
             gasConc  = gasConc - dGasConc
             aerConc  = aerConc + dGasConc

             ! Final Hg2 gas concentration [molec/cm3]
             Spc(S)%Conc(I,J,L)  = gasConc

             !---------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Hg2 mass transferred from gas to strat aerosol [molec/cm3/s]
             !---------------------------------------------------------------
             IF ( State_Diag%Archive_Hg2GasToHg2StrP ) THEN
                State_Diag%Hg2GasToHg2StrP(I,J,L) =                         &
                State_Diag%Hg2GasToHg2StrP(I,J,L) + ( dGasConc / DT )
             ENDIF
          ENDDO

          ! Update Hg2 stratospheric aerosol concentration [molec/cm3]
          Spc(id_Hg2STRP)%Conc(I,J,L) = aerConc

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer memory
    Spc => NULL()

  END SUBROUTINE PartitionHg2
!!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: red_inplume_grid
!
! !DESCRIPTION: Subroutine RED\_INPLUME\_GRID conducts in plume reduction of
!  Hg2 for selected grids
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Red_InPlume_Grid( I, J, E_plant )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J
    REAL(fp), INTENT(IN) :: E_plant
!
! !REVISION HISTORY:
!  11 Jan 2011 - Y. Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: E_deg
!
! !DEFINED PARAMETERS:
!
    ! The source species profile (0.2 Hg0)
    ! 0.78 Hg2, 0.02 HgP, Streets, 2005), we use 75%, very high end
    REAL(fp),  PARAMETER :: K_Red_InPlume = 7.5e-1_fp

    !=================================================================
    ! RED_INPLUME_GRID begins here!
    !=================================================================

    ! Calculate the mass of Hg2 been degraded in plume
    E_deg = K_Red_InPlume * E_plant

    ! Subtract this part of Hg2 from the emission
    EHg2_an(I,J) = EHg2_an(I,J) - E_deg

    ! Degraded to Hg0
    EHg0_an(I,J) = EHg0_an(I,J) + E_deg

  END SUBROUTINE Red_InPlume_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_red_inplume
!
! !DESCRIPTION: Subroutine DO\_RED\_INPLUME conducts in plume reduction of
!  Hg2 for selected grids.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Red_InPlume( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_State_GC_Mod, ONLY : HcoState
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Grid_Mod,   ONLY : GrdState
    USE Time_Mod,         ONLY : Expand_Date
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  11 Jan 2011 - Y. Zhang    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I, J
    REAL(fp)            :: E_plant

    ! Strings
    CHARACTER(LEN=255)  :: thisLoc
    CHARACTER(LEN=512)  :: errMsg

    ! Pointers
    REAL(f4), POINTER   :: E_ELEC_Hg2(:,:)
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: SEC_PER_YR = 365.25_fp * 86400_fp

    !=================================================================
    ! DO_RED_INPLUME begins here!
    !=================================================================

    ! Initialize
    RC         =  GC_SUCCESS
    errMsg     =  ''
    thisLoc    =  ' -> at DO_RED_INPLUME (in GeosCore/mercury_mod.F90)'
    E_ELEC_Hg2 => NULL()

    ! Get a pointer to the monthly mean OH from HEMCO (bmy, 3/11/15)
    CALL HCO_GetPtr( HcoState, 'CFPP_NEI2005_Hg2', E_ELEC_Hg2, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Cannot get pointer to HEMCO field CFPP_NEI2005_Hg2!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Do the reduction
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Hg emission in CFPP sector decreased for 45.14% during 2005-2010
       ! Convert [kg/m2/s] --> [kg/s]
       E_plant = E_ELEC_Hg2(I,J) * State_Grid%Area_M2(I,J)

       ! Reduce the Hg2 from plume
       CALL RED_INPLUME_GRID( I, J, E_plant )

    ENDDO
    ENDDO

    ! Free npointer
    E_ELEC_Hg2 => NULL()

  END SUBROUTINE Do_Red_InPlume
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: offlineocean_readmo
!
! !DESCRIPTION: Subroutine OFFLINEOCEAN\_READMO reads the monthly varying
!     offline ocean evasion emissions if LDYNOCEAN is FALSE. Will not actually
!     need mixed layer depth when i get stuff from Yanxu
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OfflineOcean_ReadMo( State_Chm, State_Diag, State_Grid,         &
                                  State_Met, FLUX, RC                       )
!
! !USES:
!
    USE ErrCode_Mod
    USE TIME_MOD,          ONLY : EXPAND_DATE,GET_YEAR, GET_TS_EMIS
    USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH, GET_MONTH
    USE State_Chm_Mod,     ONLY : ChmState
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE HCO_State_GC_Mod,  ONLY : HcoState
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State Object
!
! !OUTPUT PARAMETERS:
!
    REAL*8,         INTENT(OUT)   :: FLUX(State_Grid%NX,State_Grid%NY)
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  12 Aug 2015 - H. Horowitz - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: I, J, M, DN(12) ! M is mon., DN is days in mon
    INTEGER                 :: NNN
    INTEGER                 :: THISYEAR
    INTEGER                 :: THISMONTH

    LOGICAL, SAVE           :: FIRST = .TRUE.
    LOGICAL                 :: IS_OCEAN_BOX
    INTEGER                 :: NN, N
    REAL(fp)                :: A_M2,     DTSRCE
    REAL(fp)                :: CHg0aq,   CHg0,     vi, Hg0aqtemp
    REAL(fp)                :: TC,       TK,       Kw
    REAL(fp)                :: Sc,       ScCO2,    USQ
    REAL(fp)                :: FRAC_L,   FRAC_O,   H, D
    REAL(fp)                :: FUP(State_Grid%NX,State_Grid%NY) 
    REAL(fp)                :: FDOWN(State_Grid%NX,State_Grid%NY)
    REAL(fp)                :: Hg0aq(State_Grid%NX,State_Grid%NY)
    REAL(fp)                :: MHg0_air

    ! Conversion factor from [cm/h * ng/L] --> [kg/m2/s]
    REAL(fp),  PARAMETER    :: TO_KGM2S = 1.0e-11_fp / 3600.0e+0_fp

    ! Small numbers to avoid dividing by zero
    REAL(fp),  PARAMETER    :: SMALLNUM = 1.0e-32_fp

    REAL(fp)                :: SFCWINDSQR

    ! Characters
    CHARACTER(LEN=255)      :: thisLoc
    CHARACTER(LEN=512)      :: errMsg

    !=================================================================
    ! OFFLINEOCEAN_READMO begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at OFFLINEOCEAN_READMO (in GeosCore/mercury_mod.F90)'

    ! Get month
    THISMONTH = GET_MONTH()
    M         = THISMONTH

    ! Days in each month (will use later) 9/16/15 hmh
    DN =  (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    !-----------------------------------------------------------------
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Zero flux arrays of State_Diag to prevent leftover values
    ! from the last timestep from being included in the averaging
    !-----------------------------------------------------------------
    IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
       State_Diag%FluxHg0fromOceanToAir = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
       State_Diag%FluxHg0fromAirToOcean = 0.0_f4
    ENDIF

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! read monthly ocean evasion  !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF ( ITS_A_NEW_MONTH() ) THEN

       CALL HCO_GetPtr( HcoState, 'GLOBAL_OCEAN', OCEAN_CONC, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Cannot get pointer to HEMCO field GLOBAL_OCEAN!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Only doing Hg0 overall, should add trap for LSPLIT (cpt - 2017)
    Hg0aq = OCEAN_CONC(:,:,1)

    ! Emission timestep [s]
    DTSRCE = GET_TS_EMIS()

    ! Loop over surface boxes
    ! NOTE: Remove the loop over NN -- the tagged Hg simulation is not used
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,            A_M2,    vi,     ScCO2                     )&
    !$OMP PRIVATE( J,            TC,      TK                                )&
    !$OMP PRIVATE( N,            CHg0,    FRAC_L, FRAC_O                    )&
    !$OMP PRIVATE( H,            Kw,      CHg0aq, Hg0aqtemp, MHg0_air       )&
    !$OMP PRIVATE( IS_OCEAN_BOX, Sc,      Usq,    D                         )&
    !$OMP COLLAPSE( 2                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Zero/initialize loop varaibles
       A_M2 = State_Grid%Area_M2( I, J )   ! Area [m2]
       Kw   = 0.0_fp
       TK   = 0.0_fp
       TC   = 0.0_fp

       ! Get fractions of land and ocean in the grid box [unitless]
       ! Use fractional land type information in MERRA. Also make sure
       ! we do not use boxes that are mostly sea ice for consistency
       ! FROCEAN is a constant, so to get correct ocean fraction we
       ! need to subtract the sea ice fraction. Don't let the fraction
       ! be less than zero (jaf, 4/26/11)
       FRAC_L       = State_Met%FRLAND(I,J)
       FRAC_O       = MAX( State_Met%FROCEAN(I,J) - &
                           State_Met%FRSEAICE(I,J), 0.0_fp )
       IS_OCEAN_BOX = ( ( FRAC_O > 0.0_fp ) .and. &
                        ( State_Met%SEAICE00(I,J)  > 0.5_fp ) )

       IF ( IS_OCEAN_BOX ) THEN

          !--------------------------------------------------------------
          ! Sea surface temperature in both [K] and [C]
          !--------------------------------------------------------------
          ! where TSKIN is the temperature (K) at the ground/sea surface
          ! (Use as surrogate for SST, cap at freezing point)
          TK     = MAX( State_Met%TSKIN(I,J), 273.15_fp )
          TC     = TK - 273.15_fp

          !==============================================================
          ! Volatilisation of Hg0
          !==============================================================

          ! Henry's law constant (gas->liquid) [unitless] [L water/L air]
          ! (ref: Andersson et al. 2008)
          H      = EXP( ( -2404.3_fp / TK ) + 6.92_fp )

          ! Viscosity as a function of changing temperatures
          ! (ref: Loux 2001)
          ! The paper says the viscosity is given in cP but us really P
          ! and we therefor multiply with 100 to get cP.
          vi    = ( 10**( ( 1301.0_fp / ( 998.333_fp + 8.1855_fp             &
                  * ( TC - 20.0_fp )  + 0.00585_fp                           &
                  * ( TC - 20.0_fp )**2 ) ) - 3.30233_fp ) ) * 100.0_fp

          ! Schmidt # for Hg [unitless]
          ! Sc = v/D = kinematic viscosity/diffusivity
          ! (ref: Poissant et al 2000; Wilke and Chang 1995)
          ! to correct for seawater D0 is decreased by 6% as suggested
          ! by Wanninkhof (1992)
          D  = 7.4e-8_fp * SQRT( 2.26_fp * 18.0_fp   )                       &
             * TK        / ( ( 14.8_fp**0.6_fp ) *vi )
          Sc = ( 0.017_fp * EXP( -0.025_fp * TC ) ) / D

          ! Schmidt # of CO2 [unitless] for CO2 in seawater at 20 degrees C
          ! The value is set to a constant based on other ocean studies
          ! (Gardfeld et al. 2003, Rolfhus & Fitzgerald2004,Mason et al.2001)
          !
          ! Correction of the Schmidt # with temperature based on Poissant
          ! et al. (2000) (for freshwatersystems).
          ScCO2  = 644.7_fp + TC * ( -6.16_fp + TC * ( 0.11_fp ) )

          ! Square of surface (actually 10m) wind speed [m2/s2]
          Usq    = State_Met%U10M(I,J)**2 + State_Met%V10M(I,J)**2

          !------------------------------------------------------
          ! Parameterizations for calculating water side mass trasfer
          ! coefficient
          !------------------------------------------------------
          ! Mass transfer coefficient [cm/h], from Nightingale et al. 2000
          Kw     = ( 0.25_fp * Usq ) / SQRT( Sc / ScCO2 )

          ! Hg0 tracer number (for Spc)
          N = id_Hg0 !Hg0_Id_List(NN)

          !--------------------------------------------------------
          ! Calculate oceanic and gas-phase concentration of Hg(0)
          !--------------------------------------------------------

          ! Concentration of Hg(0) in the ocean [ng/L]
          ! now converting from Hg0aq in mol/m3 to ng/L
          CHg0aq = Hg0aq(I,J) * 200.59_fp * 1.0e9_fp / 1.0e3_fp

          ! Gas phase Hg(0) concentration: convert [kg] -> [ng/L]
          MHg0_air = State_Chm%Species(N)%Conc(I,J,1)
          CHg0     = MHg0_air *  1.0e9_fp /State_Met%AIRVOL(I,J,1)

          !--------------------------------------------------------
          ! Compute flux of Hg(0) from the ocean to the air
          !--------------------------------------------------------

          ! Compute ocean flux of Hg0 [cm/h*ng/L]
          FLUX(I,J)     = Kw * ( CHg0aq - ( CHg0 / H ) )

          !Extra diagnostic: compute flux up and flux down
          FUP(I,J)   = ( Kw * CHg0aq )
          FDOWN(I,J) = ( Kw * CHg0 / H )

          !--------------------------------------------------
          ! Convert [cm/h*ng/L] --> [kg/m2/s] --> [kg/s]
          ! Also account for ocean fraction of grid box
          FLUX(I,J)  = FLUX(I,J) * TO_KGM2S * A_M2 * FRAC_O

          ! hmh 5/11/16 reverting to old version and uncommenting here
          FUP(I,J)   = FUP(I,J)   * TO_KGM2S * A_M2 * FRAC_O
          FDOWN(I,J) = FDOWN(I,J) * TO_KGM2S * A_M2 * FRAC_O
          !--------------------------------------------------

          !--------------------------------------------------------
          ! Flux limited by ocean and atm Hg(0)
          !--------------------------------------------------------

          ! Cap the flux w/ the available Hg(0) ocean mass
          Hg0aqtemp = CHg0aq * A_M2 * FRAC_O *1.0e-8_fp

          IF ( FLUX(I,J) * DTSRCE > Hg0aqtemp ) THEN
             FLUX(I,J) = Hg0aqtemp / DTSRCE
             FUP(I,J)  = FLUX(I,J)-FDOWN(I,J)
          ENDIF

          ! Cap the neg flux w/ the available Hg(0) atm mass
          IF ( (-FLUX(I,J) * DTSRCE ) > MHg0_air ) THEN
             FLUX(I,J) = -MHg0_air / DTSRCE
          ENDIF

          ! make sure Fup and Fdown do not underflow either
          ! debug 2x2.5 diagnostic?
          FUP(I,J)   = MAX( FUP(I,J),   SMALLNUM )
          FDOWN(I,J) = MAX( FDOWN(I,J), SMALLNUM )

          !--------------------------------------------------------
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Fluxes of Hg0 from air to ocean and ocean to air [kg/s]
          ! NOTE: Implement for total Hg species at ths time
          !--------------------------------------------------------

          ! Flux of Hg0 from ocean to air [kg/s]
          IF ( State_Diag%Archive_FluxHg0fromOceanToAir ) THEN
             State_Diag%FluxHg0fromOceanToAir(I,J) = FUP(I,J)
          ENDIF

          IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
             State_Diag%FluxHg0fromAirToOcean(I,J) = FDOWN(I,J)
          ENDIF

       ELSE

          FLUX(I,J)  = 0.0_fp
          FUP(I,J)   = 0.0_fp
          FDOWN(I,J) = 0.0_fp

       ENDIF
    ENDDO
    ENDDO
   !$OMP END PARALLEL DO

  END SUBROUTINE OfflineOcean_ReadMo
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mercury
!
! !DESCRIPTION: Subroutine INIT\_MERCURY allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Mercury( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE Cmn_FJX_Mod
    USE Cmn_Size_Mod,     ONLY : nAer, nDust
    USE ErrCode_Mod
    USE Fast_JX_Mod,      ONLY : Init_FJX
    USE GcKpp_Monitor,    ONLY : Eqn_Names, Fam_Names
    USE GcKpp_Parameters, ONLY : nFam, nReact
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Species_Mod,      ONLY : Species
    USE State_Chm_Mod,    ONLY : Ind_
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm     ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid    ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag    ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  02 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE          :: IS_INIT = .FALSE.

    ! Scalars
    INTEGER                :: I,     KppId, N,         P
    INTEGER                :: p_BrO, p_ClO, p_Hg2ORGP, p_NO2

    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=512)     :: errMsg

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !========================================================================
    ! INIT_MERCURY begins here!
    !========================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Return if we have already allocated arrays
    IF ( IS_INIT ) RETURN

    ! Continue initialization
    SpcInfo   => NULL()
    errMsg    = ''
    thisLoc   = '-> at DEFINE_TAGGED_Hg (in GeosCore/mercury_mod.F90)'

    ! Write reactions
    WRITE( 6 ,'(a)' ) ' KPP Reaction Reference '
    DO N = 1, NREACT
        WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
    END DO

    !========================================================================
    ! Pre-store the KPP indices for each KPP prod/loss species or family
    !========================================================================
    IF ( nFam > 0 ) THEN

        ! Allocate mapping array for KPP Ids for ND65 bpch diagnostic
        ALLOCATE( PL_Kpp_Id( nFam ), STAT=RC )
        CALL GC_CheckVar( 'mercury_mod.F90:PL_Kpp_Id', 0, RC )
        IF ( RC /= GC_SUCCESS ) RETURN

        ! Loop over all KPP prod/loss species
        DO N = 1, nFam
            ! NOTE: KppId is the KPP ID # for each of the prod and loss
            ! diagnostic species.  This is the value used to index the
            ! KPP "VAR" array (in module gckpp_Global.F90).
            KppID = Ind_( TRIM ( Fam_Names(N) ), 'K' )

            ! If the species ID is OK, save in ND65_Kpp_Id
            PL_Kpp_Id(N) = KppId
        ENDDO

    ENDIF

    ! Set oxidant species IDs
    id_O3       = Ind_( 'O3'      )
    id_OH       = Ind_( 'OH'      )
    id_HO2      = Ind_( 'HO2'     )
    id_ClO      = Ind_( 'ClO'     )
    id_Cl       = Ind_( 'Cl'      )
    id_NO2      = Ind_( 'NO2'     )
    id_NO       = Ind_( 'NO'      )
    id_Br       = Ind_( 'Br'      )
    id_BrO      = Ind_( 'BrO'     )

    ! Locate Hg gas species
    id_Hg0      = Ind_( 'Hg0'     )
    id_HgBrNO2  = Ind_( 'HgBrNO2' )
    id_HgBrHO2  = Ind_( 'HgBrHO2' )
    id_HgBrOH   = Ind_( 'HgBrOH ' )
    id_HgBrBrO  = Ind_( 'HgBrBrO' )
    id_HgBrClO  = Ind_( 'HgBrClO' )
    id_HgBr2    = Ind_( 'HgBr2  ' )
    id_HgClNO2  = Ind_( 'HgClNO2' )
    id_HgClHO2  = Ind_( 'HgClHO2' )
    id_HgClOH   = Ind_( 'HgClOH ' )
    id_HgClBrO  = Ind_( 'HgClBrO' )
    id_HgClClO  = Ind_( 'HgClClO' )
    id_HgClBr   = Ind_( 'HgClBr'  )
    id_HgOHNO2  = Ind_( 'HgOHNO2' )
    id_HgOHHO2  = Ind_( 'HgOHHO2' )
    id_HgOHOH   = Ind_( 'HgOHOH ' )
    id_HgOHBrO  = Ind_( 'HgOHBrO' )
    id_HgOHClO  = Ind_( 'HgOHClO' )
    id_HgCl2    = Ind_( 'HgCl2'   )
    id_Hg2ClP   = Ind_( 'Hg2ClP'  )
    id_Hg2ORGP  = Ind_( 'Hg2ORGP' )
    id_Hg2STRP  = Ind_( 'Hg2STRP' )
    id_HgBr     = Ind_( 'HgBr'    )
    id_HgCl     = Ind_( 'HgCl'    )
    id_HgOH     = Ind_( 'HgOH'    )
    id_HgBrO    = Ind_( 'HgBrO'   )
    id_HgClO    = Ind_( 'HgClO'   )
    id_HgOHO    = Ind_( 'HgOHO'   )

    ! Initialize variables
    nHg2gasSpc = 0
    Map_Hg2gas = 0

    IF ( id_HGBrNO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBrNO2
    ENDIF
    IF ( id_HGBrHO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBrHO2
    ENDIF
    IF ( id_HGBrOH  > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBrOH
    ENDIF
    IF ( id_HGBrBrO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBrBrO
    ENDIF
    IF ( id_HGBrClO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBrClO
    ENDIF
    IF ( id_HGBr2 > 0   ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGBr2
    ENDIF
    IF ( id_HGClNO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGClNO2
    ENDIF
    IF ( id_HGClHO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGClHO2
    ENDIF
    IF ( id_HGClOH  > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HgClOH
    ENDIF
    IF ( id_HGClBrO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGClBrO
    ENDIF
    IF ( id_HGClClO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGClClO
    ENDIF
    IF ( id_HGClBr > 0  ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGClBr
    ENDIF
    IF ( id_HGOHNO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGOHNO2
    ENDIF
    IF ( id_HGOHHO2 > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGOHHO2
    ENDIF
    IF ( id_HGOHOH  > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HgOHOH
    ENDIF
    IF ( id_HGOHBrO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGOHBrO
    ENDIF
    IF ( id_HGOHClO > 0 ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGOHClO
    ENDIF
    IF ( id_HGCl2  > 0  ) THEN
       nHg2gasSpc           = nHg2gasSpc + 1
       Map_Hg2gas(nHg2gasSpc) = id_HGCl2
    ENDIF

    !========================================================================
    ! Allocate module arrays
    !========================================================================
    ALLOCATE( COSZM( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:COSZM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    COSZM = 0e+0_fp

    ALLOCATE( EHg0_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_an', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_an = 0e+0_fp

    ALLOCATE( EHg0_dist( State_Grid%NX, State_Grid%NY), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_dist', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_dist = 0e+0_fp

    ALLOCATE( EHg0_ln( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_ln', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_ln = 0e+0_fp

    ALLOCATE( EHg0_oc( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_oc', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_oc = 0e+0_fp

    ALLOCATE( EHg0_so( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_snow', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_so = 0e+0_fp

    ALLOCATE( EHg0_snow( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg0_snow', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg0_snow = 0e+0_fp

    ALLOCATE( EHg2_an( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:EHg2_an', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    EHg2_an = 0e+0_fp

    ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:TCOSZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOSZ = 0e+0_fp

    ALLOCATE( TTDAY( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:TTDAY', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TTDAY = 0e+0_fp

    ! Allocate ZERO_DVEL if we use non-local PBL mixing or
    ! if drydep is turned off
    IF ( Input_Opt%LNLPBL .OR. ( .not. Input_Opt%LDRYD ) ) THEN
       ALLOCATE( ZERO_DVEL( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'mercury_mod.F90:ZERO_DVEL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ZERO_DVEL = 0e+0_fp
    ENDIF

    ALLOCATE( HG2_SEASALT_LOSSRATE( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:HG2_SEASALT_LOSSRATE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HG2_SEASALT_LOSSRATE = 0e+0_fp

    ALLOCATE( srMw( nHg2GasSpc ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:srMw', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    srMw = 0.0_fp

    ! HG_EMIS is needed for non-local PBL mixing
    ALLOCATE( HG_EMIS( State_Grid%NX, State_Grid%NY, State_Chm%nAdvect ),    &
              STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:HG_EMIS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HG_EMIS = 0e+0_fp

    !========================================================================
    ! Allocate and initialize oxidant concentration pointer
    !========================================================================
    ALLOCATE( FixSpcPtr( State_Chm%nKppFix ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:FixSpcPtr', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Number of aerosol & dust speices
    N_Aer  =  NAER
    N_Dust =  NDUST

    ! Aerosol species name
    ALLOCATE( AerSpcNames ( N_Dust + N_Aer ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:AerSpcNames', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    AerSpcNames = (/'DST1  ','DST2  ','DST3  ','DST4  ','DST5  ','DST6  ',   &
                    'DST7  ','SO4   ','BC    ','OC    ','SSA   ','SSC   ',   &
                    'BGSULF','ICEI  '                                      /)

    ALLOCATE( AeroPtr( N_Dust + N_Aer ), STAT=RC )
    CALL GC_CheckVar( 'mercury_mod.F90:AeroPtr', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !========================================================================
    ! Pre-save square root of Hg2 species molecular weights, which are
    ! needed within PartitionHg2.  This saves unnecesary CPU cycles.
    !========================================================================
    DO N = 1, nHg2gasSpc
       P       = Map_Hg2gas(N)
       srMw(N) = SQRT( State_Chm%SpcData(P)%Info%MW_g )

       ! Also save the sqrt(MW) for HgCl2
       IF ( State_Chm%SpcData(P)%Info%Name(1:5) == 'HgCl2' ) THEN
          srMw_HgCl2 = srMw(N)
       ENDIF
    ENDDO

    !========================================================================
    ! Various Settings (not sure how many of these still work)
    !========================================================================

    ! Switch uses ocean rate coefficients from parameter inversion,
    ! ref. Song et al. 2015 ACP
    LOCEANCOEF=.FALSE.

    ! Switch determines whether uptake of Hg2 by sea-salt aerosol
    ! is calculated dynamically (TRUE) or uses a constant rate (FALSE)
    LDYNSEASALT = .TRUE.

    ! Use GEIA 2005 inventory
    LGEIA05=.FALSE.

    ! Switch use NEI2005 and NPRI2005 emission inventories
    LNEI2005=.TRUE.

    ! Switch modifying the speciation profile of Hg emission
    LInPlume=.FALSE.

    ! no Hg emitted through transpiration (VEGEMIS off)
    LVEGEMIS=.FALSE.

    ! Switch adds bromine in marine boundary layer
    L_ADD_MBL_Br=.FALSE.

    ! Switch adds bromine explosion in Northern springtime
    LPOLARBr=.TRUE.

    ! Switch for only doing reduction in cloud water
    LRED_CLOUDONLY = .TRUE.

    ! Switch for using GEOS-Chem tropospheric bromine fields,
    ! ref. Parrella et al. 2012, instead of older TOMCAT fields
    LGCBROMINE = .TRUE.

    ! Switch specifies that Hg2 is 50% bound to aerosol and 50% in
    ! gas phase (TRUE). If FALSE, then use temperature dependent
    ! partitioning as described in Amos et al. (2011, ACPD)
    LHg2HalfAerosol = .FALSE.

    ! Switch turns on snowpack Hg storage until snowmelt
    LHGSNOW = .TRUE.

    ! Multiplicative factor for increasing stratospheric Br and BrO
    STRAT_Br_FACTOR = 1e+0_fp

    ! Switch turns off all emissions except direct anthropogenic
    LAnthroHgOnly = .FALSE.

    ! Switch turns off all anthropogenic emissions from contiguous USA
    LnoUSAemis = .FALSE.

    !========================================================================
    ! Initialize FAST-JX photolysis
    !========================================================================
    CALL Init_FJX( Input_Opt, State_Chm, State_Diag, State_Grid, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Init_FJX"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Find where certain species are stored in the FAST-JX rate array ZPJ
    !========================================================================

    ! Initialize photolysis indices from species database
    p_BrO          = Ind_( 'BrO',     'P' )
    p_ClO          = Ind_( 'ClO',     'P' )
    p_Hg2ORGP      = Ind_( 'Hg2ORGP', 'P' )
    p_NO2          = Ind_( 'NO2',     'P' )

    ! Initialize variables for slots of ZPJ
    id_phot_BrO    = 0
    id_phot_ClO    = 0
    id_phot_Hg2Org = 0
    id_phot_NO2    = 0

    ! Loop over all photolysis reactions
    DO N = 1, nRatJ

       ! GC photolysis species index (skip if not present)
       P = GC_Photo_Id(N)
       IF ( P <= 0 ) CYCLE

       ! Define the slots in the ZPJ array for several species.
       ! We will use this in the ChemMercury routine above.
       IF ( P == p_BrO     ) id_phot_BrO    = N
       IF ( P == p_ClO     ) id_phot_ClO    = N
       IF ( P == p_Hg2ORGP ) id_phot_Hg2Org = N
       IF ( P == p_NO2     ) id_phot_NO2    = N
    ENDDO

    ! Error checks
    IF ( id_phot_BrO <= 0 .or. id_phot_BrO > nRatJ ) THEN
       errMsg = 'Invalid photolysis index for BrO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( id_phot_ClO <= 0 .or. id_phot_ClO > nRatJ ) THEN
       errMsg = 'Invalid photolysis index for ClO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( id_phot_Hg2Org <= 0 .or. id_phot_Hg2Org > nRatJ ) THEN
       errMsg = 'Invalid photolysis index for HG2ORGP!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( id_phot_NO2 <= 0 .or. id_phot_NO2 > nRatJ ) THEN
       errMsg = 'Invalid photolysis index for NO2!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Done!  Reset IS_INIT, since we have already allocated arrays
    !=================================================================
    IS_INIT = .TRUE.

  END SUBROUTINE Init_Mercury
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_options_from_hemco
!
! !DESCRIPTION: Overrides some of the Hg simulation settings depending on
!  the inputs that are specified in the HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Options_From_Hemco( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ERROR_MOD
    USE HCO_State_GC_Mod,  ONLY : HcoState
    USE HCO_ExtList_Mod,   ONLY : GetExtOpt
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Grid_Mod,    ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
    TYPE(GrdState), INTENT(IN)  :: State_Grid
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC   ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL :: LGC
    LOGICAL :: LTOMCAT
    LOGICAL :: LPREINDHG
    LOGICAL :: FOUND

    ! Strings
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !-----------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------
    RC      = HCO_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at SET_OPTIONS_FROM_HEMCO (in GeosCore/mercury_mod.F90)'

    !-----------------------------------------------------------------
    ! Set the value of chemistry flags depending on whether or not
    ! the HEMCO collection LFLAGNAME is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'LHALOGENCHEM', &
                    OptValBool=LHALOGENCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LHALOGENCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LHALOGENCHEM = .TRUE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LHGAQCHEM', &
                    OptValBool=LHGAQCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LHGAQCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LHGAQCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LBRCHEM', &
                    OptValBool=LBRCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LBRCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LBRCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LBROCHEM', &
                    OptValBool=LBROCHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LBROCHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LBROCHEM = .FALSE.
    ENDIF

    CALL GetExtOpt( HcoState%Config, -999, 'LOHO3CHEM', &
                    OptValBool=LOHO3CHEM, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LOHO3CHEM not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LOHO3CHEM = .FALSE.
    ENDIF

    !-----------------------------------------------------------------
    ! Set the value of LGCBROMINE depending on the values of the
    ! HEMCO collections BrOx_GC and BrOx_TOMCAT
    !-----------------------------------------------------------------

    ! First look for BrOx_GC
    CALL GetExtOpt( HcoState%Config, -999, 'BrOx_GC', &
                    OptValBool=LGC, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'BrOx_GC not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LGC = .FALSE.
    ENDIF

    ! Set LGCBROMINE = .TRUE. if BrOx_GC is true
    LGCBROMINE = LGC

    ! Set BrOx_TOMCAT to be the opposite of LGCBROMINE
    LTOMCAT   = ( .not. LGCBROMINE )

    ! Are we doing a preindustrial simulation?
    LPREINDHG = Input_Opt%LPREINDHG

    !-----------------------------------------------------------------
    ! Set the value of LNEI2005 depending on whether or not
    ! the HEMCO collection NEI2005 is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'NEI2005', &
                    OptValBool=LNEI2005, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'NEI2005 not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LNEI2005 = .FALSE.
    ENDIF

    !-----------------------------------------------------------------
    ! Set the value of LInPlume depending on whether or not
    ! the HEMCO collection NEI2005 is activated
    !-----------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'LRED_INPLUME', &
                    OptValBool=LInPlume, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       errMsg = 'LRED_INPLUME not found in HEMCO_Config.rc file!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       LInPlume = .FALSE.
    ENDIF

    ! In plume degradation of Hg2 by SO2 in U.S. and Canada at CFPPs,
    ! (yzh 11/1/2011).  Move this here so that the HEMCO_Config file
    ! will have been already read by this point. (bmy, 10/11/16)
    IF ( LInPlume .AND. .NOT. LPREINDHG .AND. LNEI2005 ) THEN
       CALL DO_RED_INPLUME( Input_Opt, State_Grid, RC )
    ENDIF

    !-----------------------------------------------------------------
    ! Echo output
    !-----------------------------------------------------------------
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 100   )
    WRITE( 6, 120   ) LGCBROMINE
    WRITE( 6, 130   ) LNEI2005
    WRITE( 6, 140   ) LInPlume
    WRITE( 6, 150   ) LHALOGENCHEM
    WRITE( 6, 160   ) LHGAQCHEM
    WRITE( 6, 170   ) LBRCHEM
    WRITE( 6, 180   ) LBROCHEM
    WRITE( 6, 190   ) LOHO3CHEM
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
100 FORMAT( 'Adjusting Hg simulation settings from HEMCO inputs' )
120 FORMAT( 'LGCBROMINE is set to ', L1                          )
130 FORMAT( 'LNEI2005   is set to ', L1                          )
140 FORMAT( 'LInPlume   is set to ', L1                          )
150 FORMAT( 'LHALOGENCHEM   is set to ', L1                          )
160 FORMAT( 'LHGAQCHEM  is set to ', L1                          )
170 FORMAT( 'LBRCHEM    is set to ', L1                          )
180 FORMAT( 'LBROCHEM   is set to ', L1                          )
190 FORMAT( 'LOHO3CHEM  is set to ', L1                          )

  END SUBROUTINE Set_Options_From_Hemco
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_MERCURY deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Mercury
!
! !REVISION HISTORY:
!  06 Dec 2004 - N. (Eckley) Selin - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Deallocate module arrays
    IF ( ALLOCATED( COSZM                ) ) DEALLOCATE( COSZM                )
    IF ( ALLOCATED( EHg0_an              ) ) DEALLOCATE( EHg0_an              )
    IF ( ALLOCATED( EHg0_ln              ) ) DEALLOCATE( EHg0_ln              )
    IF ( ALLOCATED( EHg0_oc              ) ) DEALLOCATE( EHg0_oc              )
    IF ( ALLOCATED( EHg0_so              ) ) DEALLOCATE( EHg0_so              )
    IF ( ALLOCATED( EHg0_snow            ) ) DEALLOCATE( EHg0_snow            )
    IF ( ALLOCATED( EHg2_an              ) ) DEALLOCATE( EHg2_an              )
    IF ( ALLOCATED( srMw                 ) ) DEALLOCATE( srMw                 )
    IF ( ALLOCATED( TCOSZ                ) ) DEALLOCATE( TCOSZ                )
    IF ( ALLOCATED( TTDAY                ) ) DEALLOCATE( TTDAY                )
    IF ( ALLOCATED( ZERO_DVEL            ) ) DEALLOCATE( ZERO_DVEL            )
    IF ( ALLOCATED( HG2_SEASALT_LOSSRATE ) ) DEALLOCATE( HG2_SEASALT_LOSSRATE )
    IF ( ALLOCATED( HG_EMIS              ) ) DEALLOCATE( HG_EMIS              )

    ! Free pointers to HEMCO fields
    O3          => NULL()
    OH          => NULL()
    JNO2        => NULL()
    NO          => NULL()
    NO2         => NULL()
    CLO         => NULL()
    CL          => NULL()
    OA          => NULL()
    HOCl        => NULL()
    HO2         => NULL()
    OCEAN_CONC  => NULL()

  END SUBROUTINE Cleanup_Mercury
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Kpp_GridBox_Values
!
! !DESCRIPTION: Populates KPP variables in the gckpp_Global.F90 module
!  for a particular (I,J,L) grid box.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Kpp_GridBox_Values( I,         J,         L,                &
                                     Input_Opt, State_Chm, State_Grid,       &
                                     State_Met, RC                          )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE GcKpp_Global,     ONLY : State_Het
    USE Hg_HetStateFuncs, ONLY : Hg_SetStateHet
    USE Input_Opt_Mod,    ONLY : OptInput
    USE PhysConstants,    ONLY : AVO, CONSVAP, PI, RGASLATM, RSTARG
    USE Pressure_Mod,     ONLY : Get_Pcenter
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,       S
    REAL(f8)           :: CONSEXP, VPRESH2O

    ! Characters
    CHARACTER(LEN=255) :: errMsg,  thisLoc

    !========================================================================
    ! Set_Kpp_GridBox_Values begins here!
    !========================================================================

    ! Initialization
    RC      =  GC_SUCCESS
    errMsg  = ''
    thisLoc = &
       ' -> at Set_Kpp_GridBox_Values (in module GeosCore/mercury_mod.F90'

    !========================================================================
    ! Copy species concentrations into gckpp_Global variables
    !========================================================================
    DO N = 1, NSPEC
       S = State_Chm%Map_KppSpc(N)
       IF ( S > 0 ) THEN
          C(N) = State_Chm%Species(S)%Conc(I,J,L)
       ELSE
          C(N) = 0.0_f8
       ENDIF
    ENDDO

    !========================================================================
    ! Populate global variables in gckpp_Global.F90
    !========================================================================

    ! Pressure and density quantities
    NUMDEN          = State_Met%AIRNUMDEN(I,J,L)
    H2O             = State_Met%AVGW(I,J,L) * NUMDEN
    PRESS           = Get_Pcenter( I, J, L )

    ! Temperature quantities
    TEMP            = State_Met%T(I,J,L)
    INV_TEMP        = 1.0_dp   / TEMP
    TEMP_OVER_K300  = TEMP     / 300.0_dp
    K300_OVER_TEMP  = 300.0_dp / TEMP
    SR_TEMP         = SQRT( TEMP )

    ! Relative humidity quantities
    CONSEXP         = 17.2693882_f8 * (TEMP - 273.16_f8) / (TEMP - 35.86_f8)
    VPRESH2O        = CONSVAP * EXP( CONSEXP ) / TEMP

    !========================================================================
    ! Populate variables in the HetChem state object
    !========================================================================
    CALL Hg_SetStateHet(                                                     &
         I          = I,                                                     &
         J          = J,                                                     &
         L          = L,                                                     &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Met  = State_Met,                                             &
         H          = State_Het,                                             &
         fracOrgAer = GLOB_FOA(I,J,L),                                       &
         RC         = RC                                                     )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in routine "fullchem_SetStateHet"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Set_Kpp_GridBox_Values
!EOC
END MODULE Mercury_Mod
