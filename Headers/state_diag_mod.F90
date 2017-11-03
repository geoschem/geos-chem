!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: state_diag_mod.F90
!
! !DESCRIPTION: Module STATE\_DIAG\_MOD contains the derived type
!  used to define the Diagnostics State object for GEOS-Chem.
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory
!  to the Diagnostics State object.  The Diagnostics State object is not
!  defined in this module.  It must be be declared as variable in the top-level
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE:
!
MODULE State_Diag_Mod
!
! USES:
!
  USE Diagnostics_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod
  USE Species_Mod,   ONLY : Species
  USE State_Chm_Mod, ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Register_DiagField
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Diagnostics State
  !=========================================================================
  TYPE, PUBLIC :: DgnState

     !----------------------------------------------------------------------
     ! Standard Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     ! Concentrations
     REAL(f8),  POINTER :: SpeciesConc     (:,:,:,:) ! Spc Conc for diag output

     ! Dry deposition
     REAL(f4),  POINTER :: DryDepChm       (:,:,:,:) ! Drydep flux in chemistry
     REAL(f4),  POINTER :: DryDepMix       (:,:,:,:) ! Drydep flux in mixing
     REAL(f4),  POINTER :: DryDep          (:,:,:,:) ! Total drydep flux
     REAL(f4),  POINTER :: DryDepVel       (:,:,:  ) ! Dry deposition velocity
    !REAL(f4),  POINTER :: DryDepRst_RA    (:,:,:  ) ! Aerodynamic resistance
    !REAL(f4),  POINTER :: DryDepRst_RB    (:,:,:  ) ! Aerodynamic resistance
    !REAL(f4),  POINTER :: DryDepRst_RC    (:,:,:  ) ! Total drydep resistance
    !REAL(f4),  POINTER :: DryDepRst_RI    (:,:    ) ! Stomatal resistance
     ! Waiting for inputs on new resistance diagnostics commented out above

     ! Chemistry
     REAL(f4),  POINTER :: JValues         (:,:,:,:) ! Photolysis rates
     REAL(f4),  POINTER :: RxnRates        (:,:,:,:) ! Reaction rates from KPP
     REAL(f4),  POINTER :: UVFluxDiffuse   (:,:,:  ) ! Diffuse UV flux per bin
     REAL(f4),  POINTER :: UVFluxDirect    (:,:,:  ) ! Direct UV flux per bin
     REAL(f4),  POINTER :: UVFluxNet       (:,:,:  ) ! Net UV flux per bin
     
     ! Aerosols
     ! *** Need to add AOD ***
     ! Waiting for input on rest of list from Aerosol WG?

     ! Advection
     REAL(f4),  POINTER :: AdvFluxZonal    (:,:,:,:) ! EW Advective Flux
     REAL(f4),  POINTER :: AdvFluxMerid    (:,:,:,:) ! NW Advective Flux
     REAL(f4),  POINTER :: AdvFluxVert     (:,:,:,:) ! Vertical Advective Flux

     ! Mixing
     REAL(f4),  POINTER :: PBLMixFrac      (:,:,:  ) ! Frac of BL occupied by lev
     REAL(f4),  POINTER :: PBLFlux         (:,:,:,:) ! BL mixing mass flux

     ! Convection and wet deposition
     REAL(f4),  POINTER :: CloudConvFlux   (:,:,:,:) ! Cloud conv. mass flux
     REAL(f4),  POINTER :: WetLossConv     (:,:,:,:) ! Loss in convect. updraft
     REAL(f4),  POINTER :: PrecipFracConv  (:,:,:  ) ! Frac convective precip
     REAL(f4),  POINTER :: RainFracConv    (:,:,:,:) ! Frac lost to conv rainout
     REAL(f4),  POINTER :: WashFracConv    (:,:,:,:) ! Frac lost to conv washout
     REAL(f4),  POINTER :: WetLossLS       (:,:,:,:) ! Loss in LS rainout/washout
     REAL(f4),  POINTER :: PrecipFracLS    (:,:,:  ) ! Frac large scale precip
     REAL(f4),  POINTER :: RainFracLS      (:,:,:,:) ! Frac lost to LS rainout
     REAL(f4),  POINTER :: WashFracLS      (:,:,:,:) ! Frac lost to LS washout
     
     ! Carbon aerosols
     REAL(f4),  POINTER :: ProdBCPIfromBCPO(:,:,:  ) ! Prod BCPI from BCPO
     REAL(f4),  POINTER :: ProdOCPIfromOCPO(:,:,:  ) ! Prod OCPI from OCPO

     ! For isoprene SOA (moved here from State_Chm)
     REAL(fp),  POINTER :: pHSav           (:,:,:  ) ! ISORROPIA aerosol pH
     REAL(fp),  POINTER :: HplusSav        (:,:,:  ) ! H+ concentration [M]
     REAL(fp),  POINTER :: WaterSav        (:,:,:  ) ! ISORROPIA aerosol H2O
     REAL(fp),  POINTER :: SulRatSav       (:,:,:  ) ! Sulfate conc [M]
     REAL(fp),  POINTER :: NaRatSav        (:,:,:  ) ! Nitrate conc [M]
     REAL(fp),  POINTER :: AcidPurSav      (:,:,:  ) !
     REAL(fp),  POINTER :: BiSulSav        (:,:,:  ) ! Bisulfate conc [M]
 
     !----------------------------------------------------------------------
     ! Specialty Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     ! Radon / Lead / Beryllium specialty simulation
     REAL(f4),  POINTER :: PbFromRnDecay   (:,:,:  ) ! Pb emitted from Rn decay
     REAL(f4),  POINTER :: RadDecay        (:,:,:,:) ! Radioactive decay

     ! TOMAS aerosol microphysics specialty simulation
 
     ! CO2 specialty simulation

     ! CH4 specialty simulation
 
     ! Persistent Organic Pollutants specialty simulation

     ! Hg specialty simulation

     ! Radiation simulation
     REAL(f4),  POINTER :: RadAllSkyLWSurf(:,:,:  ) ! All-sky LW rad @ surface
     REAL(f4),  POINTER :: RadAllSkyLWTOA (:,:,:  ) ! All-sky LW rad @ atm top
     REAL(f4),  POINTER :: RadAllSkySWSurf(:,:,:  ) ! All-sky SW rad @ surface
     REAL(f4),  POINTER :: RadAllSkySWTOA (:,:,:  ) ! All-sky SW rad @ atm top
     REAL(f4),  POINTER :: RadClrSkyLWSurf(:,:,:  ) ! Clear-sky SW rad @ surface
     REAL(f4),  POINTER :: RadClrSkyLWTOA (:,:,:  ) ! Clear-sky LW rad @ atm top
     REAL(f4),  POINTER :: RadClrSkySWSurf(:,:,:  ) ! Clear-sky SW rad @ surface
     REAL(f4),  POINTER :: RadClrSkySWTOA (:,:,:  ) ! Clear-sky SW rad @ atm top
    
     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Diag
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'DIAG'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object

  END TYPE DgnState
!
! !REMARKS:
!  TBD
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca - Initial version
!  22 Sep 2017 - E. Lundgren - Fill in content to allocate State_Diag; add 
!                              subroutines to get metadata and interface to
!                              register fields
!  26 Sep 2017 - E. Lundgren - Remove Lookup_State_Diag and Print_State_Diag
!  05 Oct 2017 - R. Yantosca - Add separate drydep fields for chem & mixing
!  06 Oct 2017 - R. Yantosca - Declare SpeciesConc as an 8-byte real field
!  02 Nov 2017 - R. Yantosca - Update wetdep and convection diagnostic names
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Register_DiagField
     MODULE PROCEDURE Register_DiagField_R4_2D
     MODULE PROCEDURE Register_DiagField_R4_3D
     MODULE PROCEDURE Register_DiagField_R4_4D
     MODULE PROCEDURE Register_DiagField_Rfp_3D
     MODULE PROCEDURE Register_DiagField_R8_4D
  END INTERFACE Register_DiagField
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_State_Diag
!
! !DESCRIPTION: Subroutine INIT\_STATE\_DIAG allocates all fields of
!  the diagnostics state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Diag( am_I_Root, IM, JM, LM, Input_Opt, &
                              State_Chm, Diag_List,  State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # latitudes
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes
    INTEGER,        INTENT(IN)    :: LM          ! # levels
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry state object
    TYPE(DgnList),  INTENT(IN)    :: Diag_List   ! Diagnostics list object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca -  Initial version
!  22 Sep 2017 - E. Lundgren -  Fill in content
!  06 Oct 2017 - R. Yantosca -  State_Diag%SpeciesConc is now an 8-byte real
!  11 Oct 2017 - R. Yantosca -  Bug fix: nAdvect is now defined properly
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc
    CHARACTER(LEN=255)     :: arrayID,  diagID
    INTEGER                :: nSpecies, nAdvect, nDryDep
    INTEGER                :: nKppSpc,  nWetDep, N
    LOGICAL                :: EOF,      Found

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC       =  GC_SUCCESS
    arrayID  = ''
    diagID   = ''
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Found    = .FALSE.

    ! Number of species per category
    nSpecies = State_Chm%nSpecies
    nAdvect  = State_Chm%nAdvect
    nDryDep  = State_Chm%nDryDep
    nKppSpc  = State_Chm%nKppSpc
    nWetDep  = State_Chm%nWetDep

    ! Free pointers
    State_Diag%SpeciesConc         => NULL()
    State_Diag%DryDepChm           => NULL()
    State_Diag%DryDepMix           => NULL()
    State_Diag%DryDep              => NULL()
    State_Diag%DryDepVel           => NULL()
    State_Diag%JValues             => NULL()
    State_Diag%RxnRates            => NULL()
    State_Diag%UVFluxDiffuse       => NULL()
    State_Diag%UVFluxDirect        => NULL()
    State_Diag%UVFluxNet           => NULL()
    State_Diag%AdvFluxZonal        => NULL()
    State_Diag%AdvFluxMerid        => NULL()
    State_Diag%AdvFluxVert         => NULL()
    State_Diag%PBLMixFrac          => NULL()
    State_Diag%PBLFlux             => NULL()
    State_Diag%CloudConvFlux       => NULL()
    State_Diag%WetLossConv         => NULL()
    State_Diag%PrecipFracConv      => NULL()
    State_Diag%RainFracConv        => NULL()
    State_Diag%WashFracConv        => NULL()
    State_Diag%WetLossLS           => NULL()
    State_Diag%PrecipFracLS        => NULL()
    State_Diag%RainFracLS          => NULL()
    State_Diag%WashFracLS          => NULL()
    State_Diag%PbFromRnDecay       => NULL()
    State_Diag%RadDecay            => NULL()
    State_Diag%RadAllSkyLWSurf     => NULL()
    State_Diag%RadAllSkyLWTOA      => NULL()
    State_Diag%RadAllSkySWSurf     => NULL()
    State_Diag%RadAllSkySWTOA      => NULL()
    State_Diag%RadClrSkyLWSurf     => NULL()
    State_Diag%RadClrSkyLWTOA      => NULL()
    State_Diag%RadClrSkySWSurf     => NULL()
    State_Diag%RadClrSkySWTOA      => NULL()
    State_Diag%ProdBCPIfromBCPO    => NULL()
    State_Diag%ProdOCPIfromOCPO    => NULL()
    State_Diag%pHSav              => NULL()
    State_Diag%HplusSav           => NULL()
    State_Diag%WaterSav           => NULL()
    State_Diag%SulRatSav          => NULL()
    State_Diag%NaRatSav           => NULL()
    State_Diag%AcidPurSav         => NULL()
    State_Diag%BisulSav           => NULL()

#if defined( NC_DIAG )

    ! Write header
    WRITE( 6, 10 )
 10 FORMAT( /, 'Allocating the following fields of the State_Diag object:' )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    !--------------------------------------------
    ! Species Concentration
    !--------------------------------------------
    arrayID = 'State_Diag%SpeciesConc'
    diagID  = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc = 0.0_f8
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%SpeciesConc, &
                                State_Chm, State_Diag, RC )
    ENDIF

    !--------------------------------------------
    ! Dry deposition flux from chemistry
    !--------------------------------------------
    arrayID = 'State_Diag%DryDepChm'
    diagID  = 'DryDepChm'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepChm( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( ArrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepChm = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepChm, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Dry deposition flux from mixing
    !--------------------------------------------
    arrayID = 'State_Diag%DryDepMix'
    diagID  = 'DryDepMix'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepMix( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepMix = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepMix, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Total dry deposition flux
    !--------------------------------------------
    arrayID = 'State_Diag%DryDep'
    diagID  = 'DryDep'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDep( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDep = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDep, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Dry deposition velocity
    !-------------------------------------------- 
    arrayID = 'State_Diag%DryDepVel'
    diagID  = 'DryDepVel'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepVel( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepVel, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! J-Values
    !-------------------------------------------- 
    ! TODO: Mapping needs work
    arrayID = 'State_Diag%JValues'
    diagID  = 'JVal'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%JValues( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JValues = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%JValues, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! KPP Reaction Rates
    !-------------------------------------------- 
    arrayID = 'State_Diag%RxnRates'
    diagID  = 'RxnRates'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RxnRates( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RxnRates = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RxnRates, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Diffuse UV flux per wavelength bin
    !-------------------------------------------- 
    arrayID = 'State_Diag%UVFluxDiffuse'
    diagID  = 'UVFluxDiffuse'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%UVFluxDiffuse( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxDiffuse = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%UVFluxDiffuse, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Direct UV flux per wavelength bin
    !-------------------------------------------- 
    arrayID = 'State_Diag%UVFluxDirect'
    diagID  = 'UVFluxDirect'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%UVFluxDirect( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxDirect = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%UVFluxDirect, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Net UV flux per wavelength bin
    !-------------------------------------------- 
    arrayID = 'State_Diag%UVFluxNet'
    diagID  = 'UVFluxNet'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%UVFluxNet( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%UVFluxNet = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%UVFluxNet, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Zonal Advective Flux (east positive)
    !-------------------------------------------- 
    arrayID = 'State_Diag%AdvFluxZonal'
    diagID  = 'AdvFluxZonal'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxZonal( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxZonal = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxZonal, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Meridional Advective Flux (south positive)
    !-------------------------------------------- 
    arrayID = 'State_Diag%AdvFluxMerid'
    diagID  = 'AdvFluxMerid'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxMerid( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxMerid = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxMerid, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !---------------------------------------------
    ! Vertical Advective Flux (downwards positive)
    !--------------------------------------------- 
    arrayID = 'State_Diag%AdvFluxVert'
    diagID  = 'AdvFluxVert'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxVert( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxVert = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxVert, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of BL occupied by level L
    !-------------------------------------------- 
    arrayID = 'State_Diag%PBLMixFrac'
    diagID  = 'PBLMixFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLMixFrac( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLMixFrac = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PBLMixFrac, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Mass change due to boundary layer mixing
    !-------------------------------------------- 
    arrayID = 'State_Diag%PBLFlux'
    diagID  = 'PBLFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PBLFlux, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Mass change due to cloud convection
    !-------------------------------------------- 
    arrayID = 'State_Diag%CloudConvFlux'
    diagID  = 'CloudConvFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%CloudConvFlux( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CloudConvFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%CloudConvFlux, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------
    ! Loss of soluble species in convective updrafts
    !-----------------------------------------------
    arrayID = 'State_Diag%WetLossConv'
    diagID  = 'WetLossConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossConv( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConv = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WetLossConv, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !---------------------------------------------------------
    ! Fraction of grid box undergoing convective precipitation
    !---------------------------------------------------------
    arrayID = 'State_Diag%PrecipFracConv'
    diagID  = 'PrecipFracConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PrecipFracConv( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracConv = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PrecipFracConv, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Fraction of soluble species lost to rainout in convective precipitation
    !------------------------------------------------------------------------ 
    arrayID = 'State_Diag%RainFracConv'
    diagID  = 'RainFracConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RainFracConv( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracConv = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RainFracConv, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of soluble species lost to rainout in convective precipitation
    !-------------------------------------------- 
    arrayID = 'State_Diag%WashFracConv'
    diagID  = 'WashFracConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WashFracConv( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracConv = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WashFracConv, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Loss of solutble species in large-scale rainout/washout
    !-------------------------------------------- 
    arrayID = 'State_Diag%WetLossLS'
    diagID  = 'WetLossLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WetLossLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of grid box undergoing large-scale precipitation
    !-------------------------------------------- 
    arrayID = 'State_Diag%PrecipFracLS'
    diagID  = 'PrecipFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PrecipFracLS( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PrecipFracLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of soluble species lost to rainout in large-scale precipitation
    !-------------------------------------------- 
    arrayID = 'State_Diag%RainFracLS'
    diagID  = 'RainFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RainFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RainFracLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of soluble species lost to washout in large-scale precipitation
    !-------------------------------------------- 
    arrayID = 'State_Diag%WashFracLS'
    diagID  = 'WashFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WashFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WashFracLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Emission of Pb210 from Rn222 decay
    !-------------------------------------------- 
    arrayID = 'State_Diag%PbFromRnDecay'
    diagID  = 'PbFromRnDecay'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PbFromRnDecay( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PbFromRnDecay = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PbFromRnDecay, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Radioactive decay of Rn, Pb, and Be7
    ! (separate into 3 different arrays??)
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadDecay'
    diagID  = 'RadDecay'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadDecay( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadDecay = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadDecay, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! All-sky LW rad @ surface
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadAllSkyLWSurf'
    diagID  = 'RadAllSkyLWSurf'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadAllSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkyLWSurf = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadAllSkyLWSurf, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! All-sky LW rad @ atm top
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadAllSkyLWTOA'
    diagID  = 'RadAllSkyLWTOA'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadAllSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkyLWTOA = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadAllSkyLWTOA, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! All-sky SW rad @ surface
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadAllSkySWSurf'
    diagID  = 'RadAllSkySWSurf'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadAllSkySWSurf( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkySWSurf = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadAllSkySWSurf, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! All-sky SW rad @ atm top 
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadAllSkySWTOA'
    diagID  = 'RadAllSkySWTOA'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadAllSkySWTOA( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadAllSkySWTOA = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadAllSkySWTOA, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Clear-sky SW rad @ surface
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadClrSkyLWSurf'
    diagID  = 'RadClrSkyLWSurf'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadClrSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkyLWSurf = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadClrSkyLWSurf, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Clear-sky LW rad @ atm top
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadClrSkyLWTOA'
    diagID  = 'RadClrSkyLWTOA'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadClrSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkyLWTOA = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadClrSkyLWTOA, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Clear-sky SW rad @ surface 
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadClrSkySWSurf'
    diagID  = 'RadClrSkySWSurf'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadClrSkySWSurf( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkySWSurf = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadClrSkySWSurf, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    !  Clear-sky SW rad @ atm top
    !-------------------------------------------- 
    arrayID = 'State_Diag%RadClrSkySWTOA'
    diagID  = 'RadClrSkySWTOA'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RadClrSkySWTOA( IM, JM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RadClrSkySWTOA = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadClrSkySWTOA, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Production of Hydrophilic BC (aka BCPI)
    ! from Hydrophobic BC (aka BCPO)
    !--------------------------------------------
    arrayID = 'State_Diag%ProdBCPIfromBCPO'
    diagID  = 'ProdBCPIfromBCPO'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ProdBCPIfromBCPO( IM, JM, LM ), STAT=RC ) 
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdBCPIfromBCPO = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID,               &
                                State_Diag%ProdBCPIfromBCPO,     &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Production of Hydrophilic OC (aka OCPI)
    ! from Hydrophobic OC (aka OCPO)
    !--------------------------------------------
    arrayID = 'State_Diag%ProdOCPIfromOCPO'
    diagID  = 'ProdOCPIfromOCPO'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ProdOCPIfromOCPO( IM, JM, LM ), STAT=RC ) 
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ProdOCPIfromOCPO = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID,               &
                                State_Diag%ProdOCPIfromOCPO,     &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! The following quantities are only relevant for either fullchem
    ! or aerosol-only simulations.
    !
    ! If any of these quantites is listed in the HISTORY.rc file, it will
    ! be archived to a netCDF diagnostic collection.  Otherwise, it will
    ! be used internally in the chemistry modules but not saved to netCDF.
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !-----------------------------------------
       ! phSav
       !-----------------------------------------
       arrayID = 'State_Diag%phSav'
       diagID  = 'pHSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%phSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%phSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%phSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !-----------------------------------------
       ! HplusSav
       !-----------------------------------------
       arrayID = 'State_Diag%HplusSav'
       diagID  = 'HplusSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%HplusSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%HplusSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%HplusSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !-----------------------------------------
       ! WaterSav
       !-----------------------------------------
       arrayID = 'State_Diag%WaterSav'
       diagID  = 'WaterSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WaterSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WaterSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WaterSav, &
                                State_Chm, State_Diag, RC )

       !-----------------------------------------
       ! SulRatSav
       !-----------------------------------------
       arrayID = 'State_Diag%SulRatSav'
       diagID  = 'SulRatSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SulRatSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SulRatSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%SulRatSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !-----------------------------------------
       ! NaRatSav
       !-----------------------------------------
       arrayID = 'State_Diag%NaRatSav'
       diagID  = 'NaRatSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%NaRatSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%NaRatSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%NaRatSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !-----------------------------------------
       ! AcidPurSav
       !-----------------------------------------
       arrayID = 'State_Diag%AcidPurSav'
       diagID  = 'AcidPurSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AcidPurSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AcidPurSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AcidPurSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !-----------------------------------------
       ! BisulSav
       !-----------------------------------------
       arrayID = 'State_Diag%BisulSav'
       diagID  = 'BisulSav'
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%BisulSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%BisulSav = 0.0_fp
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%BisulSav, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ELSE
       
       !----------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! full-chemistry or aerosol-only.  
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !----------------------------------------------------------------
       DO N = 1, 7
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1 )
                diagID = 'pHSav'
             CASE( 2 )                
                diagID = 'HplusSav'
             CASE( 3 )
                diagID = 'WaterSav'
             CASE( 4 )
                diagID = 'SulRatSav'
             CASE( 5 ) 
                diagID = 'NaRatSav'
             CASE( 6 )
                diagID = 'AcidPurSav'
             CASE( 7 ) 
                diagID = 'BisulSav'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for full-chemistry '  // &
                      'or aerosol-only simulations.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !!--------------------------------------------
    !! Template for adding more diagnostics arrays
    !! Search and replace 'xxx' with array name
    !!-------------------------------------------- 
    !arrayID = 'State_Diag%xxx'
    !diagID  = 'xxx'
    !CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    !IF ( Found ) THEN
    !   WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
    !   ALLOCATE( State_Diag%xxx( IM, JM, LM, n ), STAT=RC ) ! Edits dims
    !   CALL GC_CheckVar( arrayID, 0, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Diag%xxx = 0.0_f4
    !   CALL Register_DiagField( am_I_Root, diagID, State_Diag%xxx, &
    !                            State_Chm, State_Diag, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !ENDIF

    ! Format statement
 20 FORMAT( 1x, a32, ' is registered as: ', a )

    !=======================================================================
    ! Print information about the registered fields (short format)
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, 30 )
 30    FORMAT( /, &
            'Registered variables contained within the State_Diag object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( am_I_Root   = am_I_Root,            &
                         Registry    = State_Diag%Registry,  &
                         ShortFormat = .TRUE.,               &
                         RC          = RC                      )
#endif

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Init_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Diag
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_DIAG deallocates all fields 
!  of the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Diag( am_I_Root, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca - Initial version
!   5 Oct 2017 - R. Yantosca - Now put error trapping on deallocations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Cleanup_State_Diag (in Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Deallocate module variables
    !=======================================================================
    IF ( ASSOCIATED( State_Diag%SpeciesConc ) ) THEN
       DEALLOCATE( State_Diag%SpeciesConc, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%SpeciesConc"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF
 
    IF ( ASSOCIATED( State_Diag%DryDepChm ) ) THEN
       DEALLOCATE( State_Diag%DryDepChm, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepChm"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepMix ) ) THEN
       DEALLOCATE( State_Diag%DryDepMix, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepMix"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDep ) ) THEN
       DEALLOCATE( State_Diag%DryDep, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDep"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepVel ) ) THEN
       DEALLOCATE( State_Diag%DryDepVel, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepVel"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%JValues ) ) THEN
       DEALLOCATE( State_Diag%JValues, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%JValues"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RxnRates ) ) THEN
       DEALLOCATE( State_Diag%RxnRates, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RxnRates"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxDiffuse ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDiffuse, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%UVFluxDiffuse"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxDirect ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDirect, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%UVFluxDirect"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxNet ) ) THEN
       DEALLOCATE( State_Diag%UVFluxNet, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%UVFluxNet"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxZonal ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxZonal, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%AdvFluxZonal"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxMerid ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxMerid, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%AdvFluxMerid"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxVert ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxVert, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%AdvFluxVert"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLMixFrac ) ) THEN
       DEALLOCATE( State_Diag%PBLMixFrac, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%PBLMixFrac"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLFlux ) ) THEN
       DEALLOCATE( State_Diag%PBLFlux, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%PBLFlux"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%CloudConvFlux ) ) THEN
       DEALLOCATE( State_Diag%CloudConvFlux, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%CloudConvFlux"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossConv ) ) THEN
       DEALLOCATE( State_Diag%WetLossConv, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%WetLossConv"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%PrecipFracConv ) ) THEN
       DEALLOCATE( State_Diag%PrecipFracConv, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%PrecipFracConv"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RainFracConv ) ) THEN
       DEALLOCATE( State_Diag%RainFracConv, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RainFracConv"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%WashFracConv ) ) THEN
       DEALLOCATE( State_Diag%WashFracConv, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%WashFracConv"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossLS ) ) THEN
       DEALLOCATE( State_Diag%WetLossLS, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%WetLossLS"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%PrecipFracLS ) ) THEN
       DEALLOCATE( State_Diag%PrecipFracLS, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%PrecipFracLS"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RainFracLS ) ) THEN
       DEALLOCATE( State_Diag%RainFracLS, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RainFracLS"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%WashFracLS ) ) THEN
       DEALLOCATE( State_Diag%WashFracLS, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%WashFracLS"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%PbFromRnDecay ) ) THEN
       DEALLOCATE( State_Diag%PbFromRnDecay, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%PbFromRnDecay"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadDecay ) ) THEN
       DEALLOCATE( State_Diag%RadDecay, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadDecay"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWSurf, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadAllSkyLWSurf"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWTOA, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadAllSkyLWTOA"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWSurf, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadAllSkySWSurf"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWTOA, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadAllSkySWTOA"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWSurf, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadClrSkyLWSurf"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWTOA, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadClrSkyLWTOA"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWSurf, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadClrSkySWSurf"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWTOA, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%RadClrSkySWTOA"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdBCPIfromBCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdBCPIfromBCPO, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%%ProdBCPIfromBCPO"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdOCPIfromOCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdOCPIfromOCPO, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%ProdOCPIfromOCPO"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%phSav ) ) THEN
       DEALLOCATE( State_Diag%phSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%phSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%HplusSav ) ) THEN
       DEALLOCATE( State_Diag%HplusSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%HplusSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%WaterSav ) ) THEN
       DEALLOCATE( State_Diag%WaterSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%WaterSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%SulRatSav ) ) THEN
       DEALLOCATE( State_Diag%SulRatSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%SulRatSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%NaRatSav ) ) THEN
       DEALLOCATE( State_Diag%NaRatSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%Na_RatSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%AcidPurSav ) ) THEN
       DEALLOCATE( State_Diag%AcidPurSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%AcidPurSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%BisulSav ) ) THEN
       DEALLOCATE( State_Diag%BisulSav, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%BisulSav"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
    !-----------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Diag%xxx ) ) THEN
    !   DEALLOCATE( State_Diag%xxx, STAT=RC  )
    !   IF ( RC /= GC_SUCCESS ) THEN
    !      ErrMsg = 'Could not deallocate "State_Diag%xxx"!'
    !      CALL GC_Error( ErrMsg, RC, ThisLoc )
    !      RETURN
    !   ENDIF
    !ENDIF

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, State_Diag%Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Diag%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Metadata_State_Diag
!
! !DESCRIPTION: Subroutine GET\_METADATA\_STATE\_DIAG retrieves basic 
!  information about each State_Diag field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root,  metadataID, Found,    &
                                      RC,         Desc,       Units,    &
                                      PerSpecies, Rank,       Type,     &
                                      VLoc )
!
! !USES:
!
    USE Charpak_Mod,         ONLY: To_UpperCase
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
    CHARACTER(LEN=*),    INTENT(IN)  :: metadataID ! State_Diag field ID
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT)           :: Found      ! Item found?
    INTEGER,             INTENT(OUT)           :: RC         ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc       ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units      ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: PerSpecies ! Max spc wildcard
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank       ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type       ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc       ! Vert placement
!
! !REMARKS:
!  If a diagnostic cannot use a wildcard, then set PerSpecies=''.
!
! !REVISION HISTORY: 
!  20 Sep 2017 - E. Lundgren - Initial version
!  06 Oct 2017 - R. Yantosca - State_Diga%SpeciesConc is now an 8-byte real
!  01 Nov 2017 - R. Yantosca - Now get To_UpperCase from charpak_mod.F90
!  02 Nov 2017 - R. Yantosca - Update metadata to be consistent w/ arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Name_AllCaps
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc, isSpecies
    
    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Diag (in Headers/state_diag_mod.F90)'
    Found = .TRUE.

    ! Optional arguments present?
    isDesc    = PRESENT( Desc       )
    isUnits   = PRESENT( Units      )
    isRank    = PRESENT( Rank       )
    isType    = PRESENT( Type       )
    isVLoc    = PRESENT( VLoc       )
    isSpecies = PRESENT( PerSpecies )

    ! Set defaults for optional arguments. Assume type and vertical 
    ! location are real (flexible precision) and center unless specified 
    ! otherwise
    IF ( isUnits   ) Units = ''
    IF ( isDesc    ) Desc  = ''              
    IF ( isRank    ) Rank  = -1              ! Initialize # dims as bad value 
    IF ( isType    ) Type  = KINDVAL_F4      ! Assume real*4 for diagnostics
    IF ( isVLoc    ) VLoc  = VLocationCenter ! Assume vertically centered
    IF ( isSpecies ) PerSpecies = ''         ! Assume not per species

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM( Name_AllCaps ) )

       CASE ( 'SPECIESCONC' )
          IF ( isDesc    ) Desc       = 'Dry mixing ratio of species'
          IF ( isUnits   ) Units      = 'mol mol-1 dry'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ALL'
          IF ( isType    ) Type       = KINDVAL_F8

       CASE ( 'DRYDEPCHM' )
          IF ( isDesc    ) Desc       = 'Dry deposition flux of species, from chemistry'
          IF ( isUnits   ) Units      = 'molec cm-2 s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'DRYDEPMIX' )
          IF ( isDesc    ) Desc       = 'Dry deposition flux of species, from mixing'
          IF ( isUnits   ) Units      = 'molec cm-2 s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'DRYDEP' )
          IF ( isDesc    ) Desc       = 'Dry deposition flux of species'
          IF ( isUnits   ) Units      = 'molec cm-2 s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'DRYDEPVEL' )
          IF ( isDesc    ) Desc       = 'Dry deposition velocity of species'
          IF ( isUnits   ) Units      = 'cm s-1'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'JVAL' )
          IF ( isDesc    ) Desc       = 'Photolysis rate' !TODO: append to this?
          IF ( isUnits   ) Units      = 's-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'JVN' ! TODO: fix species mapping

       CASE ( 'RXNRATES' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE ( 'UVFLUXDIFFUSE' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 3

       CASE ( 'UVFLUXDIRECT' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 3

       CASE ( 'UVFLUXNET' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 3

       CASE ( 'ADVFLUXZONAL' )
          IF ( isDesc    ) Desc       = 'Advection of species in zonal direction'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'ADVFLUXMERID' )
          IF ( isDesc    ) Desc       = 'Advection of species in meridional direction'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'
     
       CASE ( 'ADVFLUXVERT' )
          IF ( isDesc    ) Desc       = 'Advection of species in vertical direction'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'PBLMIXFRAC' )
          IF ( isDesc    ) Desc       = 'Fraction of boundary layer occupied by each level'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'PBLFLUX' )
          IF ( isDesc    ) Desc       = 'Species mass change due to boundary-layer mixing'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'CLOUDCONVFLUX' )
          IF ( isDesc    ) Desc       = 'Mass change due to cloud convection'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'WETLOSSCONV' )
          IF ( isDesc    ) Desc       = 'Loss of soluble species in convective updrafts'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'PRECIPFRACCONV' )
          IF ( isDesc    ) Desc       = 'Fraction of grid box undergoing convective precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3

       CASE ( 'RAINFRACCONV' )
          IF ( isDesc    ) Desc       = 'Fraction of soluble species lost to rainout in convective precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'WASHFRACCONV' )
          IF ( isDesc    ) Desc       = 'Fraction of soluble species lost to washout in convective precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'WETLOSSLS' )
          IF ( isDesc    ) Desc       = 'Loss of soluble species in large-scale precipitation'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'PRECIPFRACLS' )
          IF ( isDesc    ) Desc       = 'Fraction of grid box undergoing large-scale precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3

       CASE ( 'RAINFRACLS' )
          IF ( isDesc    ) Desc       = 'Fraction of soluble species lost to rainout in large-scale precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'WASHFRACLS' )
          IF ( isDesc    ) Desc       = 'Fraction of soluble species lost to washout in large-scale precipitation'
          IF ( isUnits   ) Units      = '1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'PBFROMRNDECAY' )
          IF ( isDesc    ) Desc       = 'Pb210 created from radioactive decay of Rn222'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3

       CASE ( 'RADDECAY' )
          IF ( isDesc    ) Desc       = 'Radioactive decay of radionuclide species'
          IF ( isUnits   ) Units      = 'kg s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'RADALLSKYLWSURF' )
          IF ( isDesc    ) Desc       = 'All-sky long-wave radiation at surface'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYLWTOA' )
          IF ( isDesc    ) Desc       = 'All-sky long-wave radiation at top of atmosphere'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYSWSURF' )
          IF ( isDesc    ) Desc       = 'All-sky short-wave radiation at surface'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYSWTOA ' )
          IF ( isDesc    ) Desc       = 'All-sky short-wave radiation at top of atmosphere'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYLWSURF' )
          IF ( isDesc    ) Desc       = 'Clear-sky long-wave radiation at surface'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYLWTOA ' )
          IF ( isDesc    ) Desc       = 'Clear-sky long-wave radiation at top of atmosphere'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYSWSURF' )
          IF ( isDesc    ) Desc       = 'Clear-sky short-wave radiation at surface'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYSWTOA' )
          IF ( isDesc    ) Desc       = 'Clear-sky short-wave radiation at top of atmosphere'
          IF ( isUnits   ) Units      = 'W m-2'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'PRODBCPIFROMBCPO' )
          IF ( isDesc    ) Desc       = 'Production of hydrophilic black carbon from hydrophobic black carbon'
          IF ( isUnits   ) Units      = 'kg'
          IF ( isRank    ) Rank       = 3

       CASE ( 'PRODOCPIFROMOCPO' )
          IF ( isDesc    ) Desc       = 'Production of hydrophilic organic carbon from hydrophobic organic carbon'
          IF ( isUnits   ) Units      = 'kg'
          IF ( isRank    ) Rank       = 3

       CASE ( 'PHSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA aerosol pH'
          IF ( isUnits   ) Units = '1'
          IF ( isRank    ) Rank  = 3

       CASE ( 'HPLUSSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA H+ concentration'
          IF ( isUnits   ) Units = 'mol L-1'
          IF ( isRank    ) Rank  = 3

       CASE ( 'WATERSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA aerosol water concentration'
          IF ( isUnits   ) Units = 'ug m-3'
          IF ( isRank    ) Rank  = 3

       CASE ( 'SULRATSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits   ) Units = 'M'
          IF ( isRank    ) Rank  = 3

       CASE ( 'NARATSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits   ) Units = 'M'
          IF ( isRank    ) Rank  = 3

       CASE ( 'ACIDPURSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA ACIDPUR'
          IF ( isUnits   ) Units = 'M'
          IF ( isRank    ) Rank  = 3

       CASE ( 'BISULSAV' )
          IF ( isDesc    ) Desc  = 'ISORROPIA Bisulfate (general acid)' &
                                 // ' concentration'
          IF ( isUnits   ) Units = 'M'
          IF ( isRank    ) Rank  =  3

       CASE DEFAULT
          Found = .False.
          ErrMsg = 'Metadata not found for State_Diag field ID: '            &
                   // TRIM( metadataID ) // '. If the name in HISTORY.rc '   &
                   // 'has species appended, make sure the species name '    &
                   // 'is preceded by double underscores. Otherwise, '       &
                   // 'check that the name is listed with all caps in '      &
                   // 'subroutine Get_Metadata_State_Diag '                  &
                   // '(Headers/state_diag_mod.F90).'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

    END SELECT

   END SUBROUTINE Get_Metadata_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 4-byte real field of State_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_2D( am_I_Root, metadataID, Ptr2Data, &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(f4),          POINTER       :: Ptr2Data(:,:)   ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    TYPE(DgnState),    INTENT(IN)    :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type, vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( ( ( perSpecies == '' ) .AND. ( rank /= 2 ) )  &
         .OR. ( ( perSpecies /= '' ) .AND. ( rank /= 1 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
          
    IF ( perSpecies /= '' ) THEN

       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL' )
             nSpecies = State_Chm%nSpecies
          CASE ( 'ADV' )
             nSpecies = State_Chm%nAdvect
          CASE ( 'DRY' )
             nSpecies = State_Chm%nDryDep
          CASE ( 'WET' )
             nSpecies = State_Chm%nWetDep
          CASE DEFAULT
             ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                      ' is not implemented for this combo of data type and size'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       DO N = 1, nSpecies          
          ! TODO: add more cases as needed
          SELECT CASE ( perSpecies )
             CASE ( 'ALL', 'ADV' )
                D = N
             CASE ( 'DRY' )
                D =  State_Chm%Map_DryDep(N)
             CASE ( 'WET' )
                D =  State_Chm%Map_WetDep(N)
          END SELECT
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSE
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = MetadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data2d_4     = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where diagnostics is not tied to species'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_DiagField_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 4-byte real field of State_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_3D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(f4),          POINTER       :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    TYPE(DgnState),    INTENT(INOUT) :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( ( ( perSpecies == '' ) .AND. ( rank /= 3 ) )  &
         .OR. ( ( perSpecies /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( perSpecies /= '' ) THEN

       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL' )
             nSpecies = State_Chm%nSpecies
          CASE ( 'ADV' )
             nSpecies = State_Chm%nAdvect
          CASE ( 'DRY' )
             nSpecies = State_Chm%nDryDep
          CASE ( 'WET' )
             nSpecies = State_Chm%nWetDep
          CASE DEFAULT
             ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                      ' is not implemented for this combo of data type and size'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       DO N = 1, nSpecies          
          ! TODO: add more cases as needed
          SELECT CASE ( perSpecies )
             CASE ( 'ALL', 'ADV' )
                D = N
             CASE ( 'DRY' )
                D =  State_Chm%Map_DryDep(N)
             CASE ( 'WET' )
                D =  State_Chm%Map_WetDep(N)
          END SELECT
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),      &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSE
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = metadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where diagnostics is not tied to species'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_DiagField_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_4D
!
! !DESCRIPTION: Registers a 4-dimensional, 4-byte real field of State_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_4D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root         ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID        ! Name
    REAL(f4),          POINTER       :: Ptr2Data(:,:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm         ! Obj for chem state
    TYPE(DgnState),    INTENT(IN)    :: State_Diag        ! Obj for diag state
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC                ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    CHARACTER(LEN=10)      :: JNames(37) ! Temp for J-Values until in specdb
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_4D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Assume always tied to a species
    ! TODO: add more cases as needed
    SELECT CASE ( perSpecies )
       CASE ( 'ALL' )
          nSpecies = State_Chm%nSpecies
       CASE ( 'ADV' )
          nSpecies = State_Chm%nAdvect
       CASE ( 'DRY' )
          nSpecies = State_Chm%nDryDep
       CASE ( 'WET' )
          nSpecies = State_Chm%nWetDep
       CASE ( 'JVN' ) 
          ! TODO: # of J-Values. For now, hard-code them based on AD22.
          !       In the near future, build photol into species database.
          nSpecies = 37
       CASE DEFAULT
          ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                   ' is not implemented for this combo of data type and size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

    DO N = 1, nSpecies          
       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL', 'ADV' )
             D = N
          CASE ( 'DRY' )
             D =  State_Chm%Map_DryDep(N)
          CASE ( 'WET' )
             D =  State_Chm%Map_WetDep(N)
          CASE ( 'JVN' )
             ! TODO: For now, hard-code the names based on AD22. RNAMES and
             !  JLABEL are not yet initialized at this point. Put index mapping
             !  in  species database in near future if possible.
             JNAMES = [ 'NO2       ', 'HNO3      ', 'H2O2      ',            &
                        'CH2O      ', 'x         ', 'x         ',            &
                        'GLYX      ', 'MGLY      ', 'BrO       ',            &
                        'HOBr      ', 'BrNO2     ', 'BrNO3     ',            &
                        'CHBr3     ', 'Br2       ', 'O2inadj   ',            &
                        'N2O       ', 'NO        ', 'NO3       ',            &
                        'CFC11     ', 'CFC12     ', 'CCl4      ',            &
                        'CH3Cl     ', 'ACET      ', 'ALD2      ',            &
                        'MVK       ', 'MACR      ', 'HAC       ',            &
                        'GLYC      ', 'PIP       ', 'IPMN      ',            &
                        'ETHLN     ', 'DHDN      ', 'HPALD     ',            &
                        'ISN1      ', 'MONITS    ', 'MONITU    ',            &
                        'HONIT     ' ]
#if defined( UCX )
             JNAMES(5) = 'O3_O1D'
             JNAMES(6) = 'O3_O3P'
#else
             JNAMES(5) = 'O3'
             JNAMES(6) = 'O3_POH'
#endif
       END SELECT
       IF ( perSpecies /= 'JVN' ) THEN
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
       ELSE
          ! Temporary implementation for J-Values. Add to specdb in near future
          ! to avoid having to do special treatment.
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( JNAMES(N) )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( JNAMES(N) )
       ENDIF
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = thisSpcName,          &
                               Description  = thisSpcDesc,          &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data(:,:,:,N),    &
                               RC           = RC                   )
       SpcInfo => NULL()
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE Register_DiagField_R4_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_Rfp_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 4-byte real field of State_Diag,
!  so that we can include it in the netCDF diagnostic output archive.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_Rfp_3D( am_I_Root, metadataID, Ptr2Data,  &
                                        State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(fp),          POINTER       :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    TYPE(DgnState),    INTENT(INOUT) :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  03 Nov 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( ( ( perSpecies == '' ) .AND. ( rank /= 3 ) )  &
         .OR. ( ( perSpecies /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( perSpecies /= '' ) THEN

       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL' )
             nSpecies = State_Chm%nSpecies
          CASE ( 'ADV' )
             nSpecies = State_Chm%nAdvect
          CASE ( 'DRY' )
             nSpecies = State_Chm%nDryDep
          CASE ( 'WET' )
             nSpecies = State_Chm%nWetDep
          CASE DEFAULT
             ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                      ' is not implemented for this combo of data type and size'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       DO N = 1, nSpecies          
          ! TODO: add more cases as needed
          SELECT CASE ( perSpecies )
             CASE ( 'ALL', 'ADV' )
                D = N
             CASE ( 'DRY' )
                D =  State_Chm%Map_DryDep(N)
             CASE ( 'WET' )
                D =  State_Chm%Map_WetDep(N)
          END SELECT
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data2d       = Ptr2Data(:,:,N),      &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSE
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = metadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d       = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where diagnostics is not tied to species'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Register_DiagField_Rfp_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R8_4D
!
! !DESCRIPTION: Registers a 4-dimensional, 8-byte real field of State_Diag,
!  so that we can include it in the netCDF diagnostic output archive.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R8_4D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root         ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID        ! Name
    REAL(f8),          POINTER       :: Ptr2Data(:,:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm         ! Obj for chem state
    TYPE(DgnState),    INTENT(IN)    :: State_Diag        ! Obj for diag state
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC                ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    CHARACTER(LEN=10)      :: JNames(37) ! Temp for J-Values until in specdb
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R8_4D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Assume always tied to a species
    ! TODO: add more cases as needed
    SELECT CASE ( perSpecies )
       CASE ( 'ALL' )
          nSpecies = State_Chm%nSpecies
       CASE ( 'ADV' )
          nSpecies = State_Chm%nAdvect
       CASE ( 'DRY' )
          nSpecies = State_Chm%nDryDep
       CASE ( 'WET' )
          nSpecies = State_Chm%nWetDep
       CASE ( 'JVN' ) 
          ! TODO: # of J-Values. For now, hard-code them based on AD22.
          !       In the near future, build photol into species database.
          nSpecies = 37
       CASE DEFAULT
          ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                   ' is not implemented for this combo of data type and size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

    DO N = 1, nSpecies          
       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL', 'ADV' )
             D = N
          CASE ( 'DRY' )
             D =  State_Chm%Map_DryDep(N)
          CASE ( 'WET' )
             D =  State_Chm%Map_WetDep(N)
          CASE ( 'JVN' )
             ! TODO: For now, hard-code the names based on AD22. RNAMES and
             !  JLABEL are not yet initialized at this point. Put index mapping
             !  in  species database in near future if possible.
             JNAMES = [ 'NO2       ', 'HNO3      ', 'H2O2      ',            &
                        'CH2O      ', 'x         ', 'x         ',            &
                        'GLYX      ', 'MGLY      ', 'BrO       ',            &
                        'HOBr      ', 'BrNO2     ', 'BrNO3     ',            &
                        'CHBr3     ', 'Br2       ', 'O2inadj   ',            &
                        'N2O       ', 'NO        ', 'NO3       ',            &
                        'CFC11     ', 'CFC12     ', 'CCl4      ',            &
                        'CH3Cl     ', 'ACET      ', 'ALD2      ',            &
                        'MVK       ', 'MACR      ', 'HAC       ',            &
                        'GLYC      ', 'PIP       ', 'IPMN      ',            &
                        'ETHLN     ', 'DHDN      ', 'HPALD     ',            &
                        'ISN1      ', 'MONITS    ', 'MONITU    ',            &
                        'HONIT     ' ]
#if defined( UCX )
             JNAMES(5) = 'O3_O1D'
             JNAMES(6) = 'O3_O3P'
#else
             JNAMES(5) = 'O3'
             JNAMES(6) = 'O3_POH'
#endif
       END SELECT
       IF ( perSpecies /= 'JVN' ) THEN
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
       ELSE
          ! Temporary implementation for J-Values. Add to specdb in near future
          ! to avoid having to do special treatment.
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( JNAMES(N) )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( JNAMES(N) )
       ENDIF
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = thisSpcName,          &
                               Description  = thisSpcDesc,          &
                               Units        = units,                &
                               Data3d_8     = Ptr2Data(:,:,:,N),    &
                               RC           = RC                   )
       SpcInfo => NULL()
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE Register_DiagField_R8_4D
!EOC
END MODULE State_Diag_Mod
