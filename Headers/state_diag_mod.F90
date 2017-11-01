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
     REAL(f8),  POINTER :: SpeciesConc    (:,:,:,:) ! Spc Conc for diag output

     ! Dry deposition
     REAL(f4),  POINTER :: DryDepFlux_Chm (:,:,:,:) ! Drydep flux from chemistry
     REAL(f4),  POINTER :: DryDepFlux_Mix (:,:,:,:) ! Drydep flux from mixing
     REAL(f4),  POINTER :: DryDepFlux     (:,:,:,:) ! Total drydep flux
     REAL(f4),  POINTER :: DryDepVel      (:,:,:  ) ! Dry deposition velocity
     !REAL(f4),  POINTER :: DryDepRst_RA   (:,:,:  ) ! Aerodynamic resistance
     !REAL(f4),  POINTER :: DryDepRst_RB   (:,:,:  ) ! Aerodynamic resistance
     !REAL(f4),  POINTER :: DryDepRst_RC   (:,:,:  ) ! Total drydep resistance
     !REAL(f4),  POINTER :: DryDepRst_RI   (:,:    ) ! Stomatal resistance
     ! Waiting for inputs on new resistance diagnostics commented out above

     ! Chemistry
     REAL(f4),  POINTER :: JValues        (:,:,:,:) ! Photolysis rates
     REAL(f4),  POINTER :: RxnRates       (:,:,:,:) ! Reaction rates from KPP
     REAL(f4),  POINTER :: UVFluxDiffuse  (:,:,:  ) ! Diffuse UV flux per bin
     REAL(f4),  POINTER :: UVFluxDirect   (:,:,:  ) ! Direct UV flux per bin
     REAL(f4),  POINTER :: UVFluxNet      (:,:,:  ) ! Net UV flux per bin
     
     ! Aerosols
     ! *** Need to add AOD ***
     ! Waiting for input on rest of list from Aerosol WG?

     ! Advection
     REAL(f4),  POINTER :: AdvFluxZonal   (:,:,:,:) ! EW Advective Flux
     REAL(f4),  POINTER :: AdvFluxMerid   (:,:,:,:) ! NW Advective Flux
     REAL(f4),  POINTER :: AdvFluxVert    (:,:,:,:) ! Vertical Advective Flux

     ! Mixing
     REAL(f4),  POINTER :: PBLMixFrac     (:,:,:  ) ! Frac of BL occupied by lev
     REAL(f4),  POINTER :: PBLFlux        (:,:,:,:) ! BL mixing mass flux

     ! Convection and wet deposition
     REAL(f4),  POINTER :: CloudConvFlux  (:,:,:,:) ! cloud convection mass flux
     REAL(f4),  POINTER :: ConvLoss       (:,:,:,:) ! Loss in convective updraft
     REAL(f4),  POINTER :: ConvPrecipFrac (:,:,:  ) ! Frac convective precip
     REAL(f4),  POINTER :: ConvRainFrac   (:,:,:,:) ! Frac lost to conv rainout
     REAL(f4),  POINTER :: ConvWashFrac   (:,:,:,:) ! Frac lost to conv washout
     REAL(f4),  POINTER :: WetLossLS      (:,:,:,:) ! Loss in LS rainout/washout
     REAL(f4),  POINTER :: PrecipFracLS   (:,:,:  ) ! Frac large scale precip
     REAL(f4),  POINTER :: RainFracLS     (:,:,:,:) ! Frac lost to LS rainout
     REAL(f4),  POINTER :: WashFracLS     (:,:,:,:) ! Frac lost to LS washout
     
     !----------------------------------------------------------------------
     ! Specialty Simulation Diagnostic Arrays
     !----------------------------------------------------------------------

     ! Radon / Lead / Beryllium specialty simulation
     REAL(f4),  POINTER :: PbFromRnDecay  (:,:,:  ) ! Pb emitted from Rn decay
     REAL(f4),  POINTER :: RadDecay       (:,:,:,:) ! Radioactive decay

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
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc, nWetDep
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
    State_Diag%SpeciesConc     => NULL()
    State_Diag%DryDepFlux_Chm  => NULL()
    State_Diag%DryDepFlux_Mix  => NULL()
    State_Diag%DryDepVel       => NULL()
    State_Diag%JValues         => NULL()
    State_Diag%DryDepFlux_Chm  => NULL()
    State_Diag%DryDepFlux_Mix  => NULL()
    State_Diag%DryDepFlux      => NULL()
    State_Diag%DryDepVel       => NULL()
    State_Diag%JValues         => NULL()
    State_Diag%RxnRates        => NULL()
    State_Diag%UVFluxDiffuse   => NULL()
    State_Diag%UVFluxDirect    => NULL()
    State_Diag%UVFluxNet       => NULL()
    State_Diag%AdvFluxZonal    => NULL()
    State_Diag%AdvFluxMerid    => NULL()
    State_Diag%AdvFluxVert     => NULL()
    State_Diag%PBLMixFrac      => NULL()
    State_Diag%PBLFlux         => NULL()
    State_Diag%CloudConvFlux   => NULL()
    State_Diag%ConvLoss        => NULL()
    State_Diag%ConvPrecipFrac  => NULL()
    State_Diag%ConvRainFrac    => NULL()
    State_Diag%ConvWashFrac    => NULL()
    State_Diag%WetLossLS       => NULL()
    State_Diag%PrecipFracLS    => NULL()
    State_Diag%RainFracLS      => NULL()
    State_Diag%WashFracLS      => NULL()
    State_Diag%PbFromRnDecay   => NULL()
    State_Diag%RadDecay        => NULL()
    State_Diag%RadAllSkyLWSurf => NULL()
    State_Diag%RadAllSkyLWTOA  => NULL()
    State_Diag%RadAllSkySWSurf => NULL()
    State_Diag%RadAllSkySWTOA  => NULL()
    State_Diag%RadClrSkyLWSurf => NULL()
    State_Diag%RadClrSkyLWTOA  => NULL()
    State_Diag%RadClrSkySWSurf => NULL()
    State_Diag%RadClrSkySWTOA  => NULL()

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
    arrayID = 'State_Diag%DryDepFlux_Chm'
    diagID  = 'DryDepFlux_Chm'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepFlux_Chm( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( ArrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepFlux_Chm = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepFlux_Chm, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Dry deposition flux from mixing
    !--------------------------------------------
    arrayID = 'State_Diag%DryDepFlux_Mix'
    diagID  = 'DryDepFlux_Mix'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepFlux_Mix( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepFlux_Mix = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepFlux_Mix, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Total dry deposition flux
    !--------------------------------------------
    arrayID = 'State_Diag%DryDepFlux'
    diagID  = 'DryDepFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepFlux( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepFlux, &
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
    arrayID = 'State_Diag%ConvLoss'
    diagID  = 'ConvLoss'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ConvLoss( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConvLoss = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%ConvLoss, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !---------------------------------------------------------
    ! Fraction of grid box undergoing convective precipitation
    !---------------------------------------------------------
    arrayID = 'State_Diag%ConvPrecipFrac'
    diagID  = 'ConvPrecipFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ConvPrecipFrac( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConvPrecipFrac = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%ConvPrecipFrac, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Fraction of soluble species lost to rainout in convective precipitation
    !------------------------------------------------------------------------ 
    arrayID = 'State_Diag%ConvRainFrac'
    diagID  = 'ConvRainFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ConvRainFrac( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConvRainFrac = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%ConvRainFrac, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Fraction of soluble species lost to rainout in convective precipitation
    !-------------------------------------------- 
    arrayID = 'State_Diag%ConvWashFrac'
    diagID  = 'ConvWashFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%ConvWashFrac( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ConvWashFrac = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%ConvWashFrac, &
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
 
    IF ( ASSOCIATED( State_Diag%DryDepFlux_Chm ) ) THEN
       DEALLOCATE( State_Diag%DryDepFlux_Chm, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepFlux_Chm"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepFlux_Mix ) ) THEN
       DEALLOCATE( State_Diag%DryDepFlux_Mix, STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepFlux_Mix"!'
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

    IF ( ASSOCIATED( State_Diag%DryDepFlux ) ) THEN
       DEALLOCATE( State_Diag%DryDepFlux, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%DryDepFlux"!'
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

    IF ( ASSOCIATED( State_Diag%ConvLoss ) ) THEN
       DEALLOCATE( State_Diag%ConvLoss, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%ConvLoss"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%ConvPrecipFrac ) ) THEN
       DEALLOCATE( State_Diag%ConvPrecipFrac, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%ConvPrecipFrac"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%ConvRainFrac ) ) THEN
       DEALLOCATE( State_Diag%ConvRainFrac, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%ConvRainFrac"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    IF ( ASSOCIATED( State_Diag%ConvWashFrac ) ) THEN
       DEALLOCATE( State_Diag%ConvWashFrac, STAT=RC  )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not deallocate "State_Diag%ConvWashFrac"!'
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
!
! !REVISION HISTORY: 
!  20 Sep 2017 - E. Lundgren - Initial version
!  06 Oct 2017 - R. Yantosca - State_Diga%SpeciesConc is now an 8-byte real
!  01 Nov 2017 - R. Yantosca - Now get To_UpperCase from charpak_mod.F90
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

       CASE ( 'DRYDEPFLUX_CHM' )
          IF ( isDesc    ) Desc       = 'Dry deposition flux of species, from chemistry'
          IF ( isUnits   ) Units      = 'molec cm-2 s-1'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'DRYDEPFLUX_MIX' )
          IF ( isDesc    ) Desc       = 'Dry deposition flux of species, from mixing'
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

       CASE ( 'DRYDEPFLUX' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'RXNRATES' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE ( 'UVFLUXDIFFUSE' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'UVFLUXDIRECT' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'UVFLUXNET' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'ADVFLUXZONAL' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'ADVFLUXMERID' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'
     
       CASE ( 'ADVFLUXVERT' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'PBLMIXFRAC' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'PBLFLUX' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'CLOUDCONVFLUX' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'CONVLOSS' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'CONVPRECIPFRAC' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'CONVRAINFRAC' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'CONVWASHFRAC' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'WETLOSSLS' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'PRECIPFRACLS' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'RAINFRACLS' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'WASHFRACLS' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'WET'

       CASE ( 'PBFROMRNDECAY' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3

       CASE ( 'RADDECAY' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ADV'

       CASE ( 'RADALLSKYLWSURF' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYLWTOA' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYSWSURF' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADALLSKYSWTOA ' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYLWSURF' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYLWTOA ' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYSWSURF' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE ( 'RADCLRSKYSWTOA' )
          IF ( isDesc    ) Desc       = 'placeholder'
          IF ( isUnits   ) Units      = 'placeholder'
          IF ( isRank    ) Rank       = 2
          IF ( isSpecies ) PerSpecies = 'placeholder'

       CASE DEFAULT
          Found = .False.
          ErrMsg = 'Metadata not found for State_Diag field ID: ' &
                   // TRIM( metadataID ) // '. If the name in HISTORY.rc ' &
                   // 'has species appended, make sure the species name ' &
                   // 'is preceded by double underscores. Otherwise, ' &
                   // 'check that the name is listed with all caps in ' &
                   // 'subroutine Get_Metadata_State_Diag ' &
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
