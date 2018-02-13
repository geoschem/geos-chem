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

  USE CMN_Size_Mod,    ONLY : IIPAR, JJPAR, LLPAR, NDUST
  USE CMN_FJX_Mod,     ONLY : STRWVSELECT
  USE Diagnostics_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod
  USE Species_Mod,     ONLY : Species
  USE State_Chm_Mod,   ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
  PUBLIC :: Get_TagInfo
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
     REAL(f4),  POINTER :: JVal            (:,:,:,:) ! J-values, instantaneous
     REAL(f4),  POINTER :: JNoonDailyAvg   (:,:,:,:) ! Noon J-values, day avg
     REAL(f4),  POINTER :: JNoonMonthlyAvg (:,:,:,:) ! Noon J-values, mon avg
     REAL(f4),  POINTER :: RxnRates        (:,:,:,:) ! Reaction rates from KPP
     REAL(f4),  POINTER :: UVFluxDiffuse   (:,:,:  ) ! Diffuse UV flux per bin
     REAL(f4),  POINTER :: UVFluxDirect    (:,:,:  ) ! Direct UV flux per bin
     REAL(f4),  POINTER :: UVFluxNet       (:,:,:  ) ! Net UV flux per bin
     REAL(f4),  POINTER :: OHconcAfterChem (:,:,:  ) ! OH, HO2, O1D, and O3P
     REAL(f4),  POINTER :: HO2concAfterChem(:,:,:  ) !  concentrations 
     REAL(f4),  POINTER :: O1DconcAfterChem(:,:,:  ) !  upon exiting the
     REAL(f4),  POINTER :: O3PconcAfterChem(:,:,:  ) !  FlexChem solver 
     REAL(f4),  POINTER :: Loss            (:,:,:,:) ! Chemical loss of species
     REAL(f4),  POINTER :: Prod            (:,:,:,:) ! Chemical prod of species

     ! Aerosols
     ! *** Need to add AOD ***
     ! Waiting for input on rest of list from Aerosol WG?
     REAL(f4),  POINTER :: ODDust          (:,:,:  ) ! Dust optical depth
     REAL(f4),  POINTER :: ODDustBinsWL1   (:,:,:,:) ! All bins 1st WL dust OD
     REAL(f4),  POINTER :: ODDustBinsWL2   (:,:,:,:) ! All bins 2nd WL dust OD
     REAL(f4),  POINTER :: ODDustBinsWL3   (:,:,:,:) ! All bins 3rd WL dust OD
     REAL(f4),  POINTER :: PM25            (:,:,:  ) ! Part. matter, r < 2.5 um
     REAL(f4),  POINTER :: TerpeneSOA      (:,:,:  )
     REAL(f4),  POINTER :: IsopreneSOA     (:,:,:  )
     REAL(f4),  POINTER :: AromaticSOA     (:,:,:  )

     ! Advection
     REAL(f4),  POINTER :: AdvFluxZonal    (:,:,:,:) ! EW Advective Flux
     REAL(f4),  POINTER :: AdvFluxMerid    (:,:,:,:) ! NW Advective Flux
     REAL(f4),  POINTER :: AdvFluxVert     (:,:,:,:) ! Vertical Advective Flux

     ! Mixing
     REAL(f4),  POINTER :: PBLMixFrac      (:,:,:  ) ! Frac of BL occupied by lev
     REAL(f4),  POINTER :: PBLFlux         (:,:,:,:) ! BL mixing mass flux

     ! Convection
     REAL(f4),  POINTER :: CloudConvFlux   (:,:,:,:) ! Cloud conv. mass flux
     REAL(f4),  POINTER :: WetLossConv     (:,:,:,:) ! Loss in convect. updraft

     ! Wet deposition
     REAL(f4),  POINTER :: WetLossLS       (:,:,:,:) ! Loss in LS rainout/washout
     REAL(f4),  POINTER :: PrecipFracLS    (:,:,:  ) ! Frac of box in LS precip
     REAL(f4),  POINTER :: RainFracLS      (:,:,:,:) ! Frac lost to LS rainout
     REAL(f4),  POINTER :: WashFracLS      (:,:,:,:) ! Frac lost to LS washout
     
     ! Carbon aerosols
     REAL(f4),  POINTER :: ProdBCPIfromBCPO(:,:,:  ) ! Prod BCPI from BCPO
     REAL(f4),  POINTER :: ProdOCPIfromOCPO(:,:,:  ) ! Prod OCPI from OCPO

     ! Sulfur aerosols prod & loss
     REAL(f4),  POINTER :: ProdSO2fromDMSandOH        (:,:,:) 
     REAL(f4),  POINTER :: ProdSO2fromDMSandNO3       (:,:,:) 
     REAL(f4),  POINTER :: ProdSO2fromDMS             (:,:,:)   
     REAL(f4),  POINTER :: ProdMSAfromDMS             (:,:,:) 
     REAL(f4),  POINTER :: ProdNITfromHNO3uptakeOnDust(:,:,:)
     REAL(f4),  POINTER :: ProdSO4fromGasPhase        (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromH2O2inCloud     (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromO3inCloud       (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromO2inCloudMetal  (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromO3inSeaSalt     (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromOxidationOnDust (:,:,:) 
     REAL(f4),  POINTER :: ProdSO4fromUptakeOfH2SO4g  (:,:,:)
     REAL(f4),  POINTER :: ProdSO4fromHOBrInCloud     (:,:,:)
     REAL(f4),  POINTER :: ProdSO4fromSRO3            (:,:,:)
     REAL(f4),  POINTER :: ProdSO4fromSRHOBr          (:,:,:)
     REAL(f4),  POINTER :: ProdSO4fromO3s             (:,:,:)
     REAL(f4),  POINTER :: LossOHbyDMS                (:,:,:) 
     REAL(f4),  POINTER :: LossNO3byDMS               (:,:,:) 
     REAL(f4),  POINTER :: LossHNO3onSeaSalt          (:,:,:) 
     
     ! Stratospheric aerosols 
     REAL(f4),  POINTER :: StratSADS      (:,:,:) ! Surface area density, solid aerosols 
     REAL(f4),  POINTER :: StratSADL      (:,:,:) ! Surface area density, liquid aerosols

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

     ! Radiation simulation (RRTMG)
     REAL(f4),  POINTER :: RadAllSkyLWSurf (:,:,:  ) ! All-sky LW rad @ surface
     REAL(f4),  POINTER :: RadAllSkyLWTOA  (:,:,:  ) ! All-sky LW rad @ atm top
     REAL(f4),  POINTER :: RadAllSkySWSurf (:,:,:  ) ! All-sky SW rad @ surface
     REAL(f4),  POINTER :: RadAllSkySWTOA  (:,:,:  ) ! All-sky SW rad @ atm top
     REAL(f4),  POINTER :: RadClrSkyLWSurf (:,:,:  ) ! Clr-sky SW rad @ surface
     REAL(f4),  POINTER :: RadClrSkyLWTOA  (:,:,:  ) ! Clr-sky LW rad @ atm top
     REAL(f4),  POINTER :: RadClrSkySWSurf (:,:,:  ) ! Clr-sky SW rad @ surface
     REAL(f4),  POINTER :: RadClrSkySWTOA  (:,:,:  ) ! Clr-sky SW rad @ atm top

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
     ! real(fp) and real(8) 2D not implemented
     MODULE PROCEDURE Register_DiagField_R4_3D
     MODULE PROCEDURE Register_DiagField_Rfp_3D
     ! real(8) 3D not implemented
     MODULE PROCEDURE Register_DiagField_R4_4D
     ! real(fp) 4D not implemented
     MODULE PROCEDURE Register_DiagField_R8_4D
  END INTERFACE Register_DiagField
!
! !PRIVATE TYPES:
!
  ! Shadow variables from Input_Opt
  LOGICAL :: Is_UCX

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
  SUBROUTINE Init_State_Diag( am_I_Root, Input_Opt,  State_Chm, &
                              Diag_List, State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
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
!  05 Jul 2017 - R. Yantosca - Initial version
!  22 Sep 2017 - E. Lundgren - Fill in content
!  06 Oct 2017 - R. Yantosca - State_Diag%SpeciesConc is now an 8-byte real
!  11 Oct 2017 - R. Yantosca - Bug fix: nAdvect is now defined properly  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc
    CHARACTER(LEN=255)     :: arrayID,  diagID
    INTEGER                :: N, IM, JM, LM
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc
    INTEGER                :: nWetDep,  nPhotol, nProd,   nLoss
    LOGICAL                :: EOF,      Found

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC        =  GC_SUCCESS
    arrayID   = ''
    diagID    = ''
    ErrMsg    = ''
    ThisLoc   = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Found     = .FALSE.

    ! Save shadow variables from Input_Opt
    Is_UCX    = Input_Opt%LUCX
    
    ! Shorten grid parameters for readability
    IM        = IIPAR ! # latitudes
    JM        = JJPAR ! # longitudes
    LM        = LLPAR ! # levels

    ! Number of species per category
    nSpecies  = State_Chm%nSpecies
    nAdvect   = State_Chm%nAdvect
    nDryDep   = State_Chm%nDryDep
    nKppSpc   = State_Chm%nKppSpc
    nLoss     = State_Chm%nLoss
    nPhotol   = State_Chm%nPhotol
    nProd     = State_Chm%nProd
    nWetDep   = State_Chm%nWetDep

    ! Free pointers
    State_Diag%AdvFluxZonal               => NULL()
    State_Diag%AdvFluxMerid               => NULL()
    State_Diag%AdvFluxVert                => NULL()

    State_Diag%DryDep                     => NULL()
    State_Diag%DryDepChm                  => NULL()
    State_Diag%DryDepMix                  => NULL()
    State_Diag%DryDepVel                  => NULL()

    State_Diag%CloudConvFlux              => NULL()
    State_Diag%WetLossConv                => NULL()

    State_Diag%PBLMixFrac                 => NULL()
    State_Diag%PBLFlux                    => NULL()

    State_Diag%WetLossLS                  => NULL()
    State_Diag%PrecipFracLS               => NULL()
    State_Diag%RainFracLS                 => NULL()
    State_Diag%WashFracLS                 => NULL()

    State_Diag%SpeciesConc                => NULL()
    State_Diag%JVal                       => NULL()
    State_Diag%JNoonDailyAvg              => NULL()
    State_Diag%JNoonMonthlyAvg            => NULL()
    State_Diag%RxnRates                   => NULL()
    State_Diag%UVFluxDiffuse              => NULL()
    State_Diag%UVFluxDirect               => NULL()
    State_Diag%UVFluxNet                  => NULL()
    State_Diag%ProdBCPIfromBCPO           => NULL()
    State_Diag%ProdOCPIfromOCPO           => NULL()
    State_Diag%OHconcAfterChem            => NULL()
    State_Diag%HO2concAfterChem           => NULL()
    State_Diag%O1DconcAfterChem           => NULL()
    State_Diag%O3PconcAfterChem           => NULL()
    State_Diag%Loss                       => NULL()
    State_Diag%Prod                       => NULL()
    State_Diag%PM25                       => NULL()
    State_Diag%TerpeneSOA                 => NULL()
    State_Diag%IsopreneSOA                => NULL()
    State_Diag%AromaticSOA                => NULL()

    State_Diag%PbFromRnDecay              => NULL()
    State_Diag%RadDecay                   => NULL()

    State_Diag%ProdMSAfromDMS             => NULL() 
    State_Diag%ProdSO2fromDMSandOH        => NULL() 
    State_Diag%ProdSO2fromDMSandNO3       => NULL() 
    State_Diag%ProdSO2fromDMS             => NULL()   
    State_Diag%ProdSO4fromGasPhase        => NULL() 
    State_Diag%ProdSO4fromH2O2inCloud     => NULL() 
    State_Diag%ProdSO4fromO3inCloud       => NULL() 
    State_Diag%ProdSO4fromHOBrInCloud     => NULL()
    State_Diag%ProdSO4fromO2inCloudMetal  => NULL() 
    State_Diag%ProdSO4fromO3inSeaSalt     => NULL() 
    State_Diag%ProdSO4fromOxidationOnDust => NULL() 
    State_Diag%ProdSO4fromUptakeOfH2SO4g  => NULL()
    State_Diag%ProdSO4fromSRO3            => NULL()
    State_Diag%ProdSO4fromSRHOBr          => NULL()
    State_Diag%ProdSO4fromO3s             => NULL()
    State_Diag%ProdNITfromHNO3uptakeOnDust=> NULL()
    State_Diag%LossOHbyDMS                => NULL() 
    State_Diag%LossNO3byDMS               => NULL() 
    State_Diag%LossHNO3onSeaSalt          => NULL() 

    State_Diag%StratSADS                  => NULL() 
    State_Diag%StratSADL                  => NULL() 

    State_Diag%RadAllSkyLWSurf            => NULL()
    State_Diag%RadAllSkyLWTOA             => NULL()
    State_Diag%RadAllSkySWSurf            => NULL()
    State_Diag%RadAllSkySWTOA             => NULL()
    State_Diag%RadClrSkyLWSurf            => NULL()
    State_Diag%RadClrSkyLWTOA             => NULL()
    State_Diag%RadClrSkySWSurf            => NULL()
    State_Diag%RadClrSkySWTOA             => NULL()

#if defined( NC_DIAG )

    ! Write header
    IF ( am_I_Root ) THEN
    WRITE( 6, 10 )
 10 FORMAT( /, 'Allocating the following fields of the State_Diag object:' )
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    !------------------------------------------------------------------------
    ! Species Concentration
    !------------------------------------------------------------------------
    arrayID = 'State_Diag%SpeciesConc'
    diagID  = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc = 0.0_f8
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%SpeciesConc, &
                                State_Chm, State_Diag, RC )
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition flux from chemistry
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepChm'
    diagID  = 'DryDepChm'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepChm( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( ArrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepChm = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepChm, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition flux from mixing
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepMix'
    diagID  = 'DryDepMix'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepMix( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepMix = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepMix, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Total dry deposition flux
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDep'
    diagID  = 'DryDep'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDep( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDep = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDep, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Dry deposition velocity
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%DryDepVel'
    diagID  = 'DryDepVel'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%DryDepVel( IM, JM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepVel, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Zonal Advective Flux (east positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxZonal'
    diagID  = 'AdvFluxZonal'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxZonal( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxZonal = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxZonal, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Meridional Advective Flux (south positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxMerid'
    diagID  = 'AdvFluxMerid'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxMerid( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxMerid = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxMerid, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Vertical Advective Flux (downwards positive)
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%AdvFluxVert'
    diagID  = 'AdvFluxVert'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%AdvFluxVert( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%AdvFluxVert = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%AdvFluxVert, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of BL occupied by level L
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PBLMixFrac'
    diagID  = 'PBLMixFrac'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLMixFrac( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLMixFrac = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PBLMixFrac, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to boundary layer mixing
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PBLFlux'
    diagID  = 'PBLFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PBLFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PBLFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PBLFlux, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Mass change due to cloud convection
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%CloudConvFlux'
    diagID  = 'CloudConvFlux'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%CloudConvFlux( IM, JM, LM, nAdvect ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%CloudConvFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%CloudConvFlux, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of soluble species in convective updrafts
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WetLossConv'
    diagID  = 'WetLossConv'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossConv( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossConv = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WetLossConv, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss of solutble species in large-scale rainout/washout
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WetLossLS'
    diagID  = 'WetLossLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WetLossLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WetLossLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WetLossLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of grid box undergoing large-scale precipitation
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%PrecipFracLS'
    diagID  = 'PrecipFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%PrecipFracLS( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%PrecipFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%PrecipFracLS, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost to rainout in large-scale precip
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%RainFracLS'
    diagID  = 'RainFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%RainFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%RainFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%RainFracLS,   &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Fraction of soluble species lost to washout in large-scale precip
    !-----------------------------------------------------------------------
    arrayID = 'State_Diag%WashFracLS'
    diagID  = 'WashFracLS'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
       ALLOCATE( State_Diag%WashFracLS( IM, JM, LM, nWetDep ), STAT=RC )
       CALL GC_CheckVar( arrayID, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%WashFracLS = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%WashFracLS,   &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for the 
    ! Rn-Pb-Be radionuclide simulation
    !=======================================================================
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN

       !--------------------------------------------------------------------
       ! Emission of Pb210 from Rn222 decay
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PbFromRnDecay'
       diagID  = 'PbFromRnDecay'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PbFromRnDecay( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PbFromRnDecay = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%PbFromRnDecay,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Radioactive decay of Rn, Pb, and Be7
       ! (separate into 3 different arrays??)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadDecay'
       diagID  = 'RadDecay'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadDecay( IM, JM, LM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadDecay = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%RadDecay,  &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! the Rn-Pb-Be-PASV simulation.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 2
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1 ) 
                diagID = 'PbFromRnDecay'
             CASE( 2 ) 
                diagID = 'RadDecay'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for Rn-Pb-Be-PASV '   // &
                      'simulations.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !=======================================================================
    ! The following diagnostic quantities are only relevant for simulations
    ! with the RRTMG radiative transfer model turned on.
    !=======================================================================
    IF ( Input_Opt%LRAD ) THEN

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkyLWSurf'
       diagID  = 'RadAllSkyLWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWSurf = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadAllSkyLWSurf,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky LW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkyLWTOA'
       diagID  = 'RadAllSkyLWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkyLWTOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadAllSkyLWTOA,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkySWSurf'
       diagID  = 'RadAllSkySWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWSurf = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadAllSkySWSurf,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: All-sky SW rad @ atm top 
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadAllSkySWTOA'
       diagID  = 'RadAllSkySWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadAllSkySWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadAllSkySWTOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadAllSkySWTOA,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkyLWSurf'
       diagID  = 'RadClrSkyLWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWSurf = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadClrSkyLWSurf,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky LW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkyLWTOA'
       diagID  = 'RadClrSkyLWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkyLWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkyLWTOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadClrSkyLWTOA,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ surface 
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkySWSurf'
       diagID  = 'RadClrSkySWSurf'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWSurf( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWSurf = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadClrSkySWSurf,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! RRTMG: Clear-sky SW rad @ atm top
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RadClrSkySWTOA'
       diagID  = 'RadClrSkySWTOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RadClrSkySWTOA( IM, JM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RadClrSkySWTOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%RadClrSkySWTOA, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! the RRTMG radiatve transfer model.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 8
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1 ) 
                diagID = 'RadAllSkyLWSurf'
             CASE( 2 ) 
                diagID = 'RadAllSkyLWTOA'
             CASE( 3 ) 
                diagID = 'RadAllSkySWSurf'
             CASE( 4 ) 
                diagID = 'RadAllSkySWTOA'
             CASE( 5 ) 
                diagID = 'RadClrSkyLWSurf'
             CASE( 6 ) 
                diagID = 'RadClrSkyLWTOA'
             CASE( 7 ) 
                diagID = 'RadClrSkySWSurf'
             CASE( 8 ) 
                diagID = 'RadClrSkySWTOA'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for simulations '     // &
                      'with the RRTMG radiative transfer model.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !=======================================================================
    ! The following quantities are only relevant for fullchem simulations
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !--------------------------------------------------------------------
       ! KPP Reaction Rates
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%RxnRates'
       diagID  = 'RxnRates'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%RxnRates( IM, JM, LM, nSpecies ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%RxnRates = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%RxnRates,   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! J-Values (instantaneous values)
       !
       ! NOTE: Dimension array nPhotol+2 to archive special photolysis
       ! reactions for O3_O1D, O3_O3P (with UCX) or O3, POH (w/o UCX)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%JVal'
       diagID  = 'JVal'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JVal( IM, JM, LM, nPhotol+2 ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JVal = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%JVal,       &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Noontime J-values, daily average
       !
       ! NOTE: Dimension array nPhotol+2 to archive special photolysis
       ! reactions for O3_O1D, O3_O3P (with UCX) or O3, POH (w/o UCX)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%JNoonDailyAvg'
       diagID  = 'JNoonDailyAvg'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JNoonDailyAvg(IM,JM,LM,nPhotol+2), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JNoonDailyAvg = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%JNoonDailyAvg,                 &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Noontime J-values, monthly average
       !
       ! NOTE: Dimension array nPhotol+2 to archive special photolysis
       ! reactions for O3_O1D, O3_O3P (with UCX) or O3, POH (w/o UCX)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%JNoonMonthlyAvg'
       diagID  = 'JNoonMonthlyAvg'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%JNoonMonthlyAvg(IM,JM,LM,nPhotol+2), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%JNoonMonthlyAvg = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%JNoonMonthlyAvg,               &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Diffuse UV flux per wavelength bin
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%UVFluxDiffuse'
       diagID  = 'UVFluxDiffuse'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxDiffuse( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDiffuse = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%UVFluxDiffuse,                 &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Direct UV flux per wavelength bin
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%UVFluxDirect'
       diagID  = 'UVFluxDirect'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxDirect( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxDirect = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%UVFluxDirect,                  &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Net UV flux per wavelength bin
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%UVFluxNet'
       diagID  = 'UVFluxNet'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%UVFluxNet( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%UVFluxNet = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%UVFluxNet, &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! HO2 concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%HO2concAfterChem'
       diagID  = 'HO2concAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%HO2concAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%HO2concAfterChem = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%HO2concAfterChem,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! O1D concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%O1DconcAfterChem'
       diagID  = 'O1DconcAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O1DconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O1DconcAfterChem = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        & 
                                   State_Diag%O1DconcAfterChem,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! O3P concentration upon exiting the FlexChem solver
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%O3PconcAfterChem'
       diagID  = 'O3PconcAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%O3PconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%O3PconcAfterChem = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%O3PconcAfterChem,              &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! full-chemistry simulations.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 10
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  ) 
                diagID = 'RxnRates'
             CASE( 2  ) 
                diagID = 'JVal'
             CASE( 3  ) 
                diagID = 'JNoonDailyAvg'
             CASE( 4  ) 
                diagID = 'JNoonMonthlyAvg'
             CASE( 5  ) 
                diagID = 'UvFluxDiffuse'
             CASE( 6  ) 
                diagID = 'UvFluxDirect'
             CASE( 7  ) 
                diagID = 'UvFluxNet'
             CASE( 8  )
                diagID = 'HO2concAfterChem'
             CASE( 9  )
                diagID = 'O1DconcAfterChem'
             CASE( 10 )
                diagID = 'O3PconcAfterChem'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for full-chemistry '  // &
                      'simulations.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ENDIF

    !=======================================================================
    ! The following quantities are only relevant for either fullchem
    ! or CH4 simulations.
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_A_CH4_SIM ) THEN

       !--------------------------------------------------------------------
       ! OH concentration upon exiting the FlexChem solver (fullchem
       ! simulations) or the CH4 specialty simulation chemistry routine
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%OHconcAfterChem'
       diagID  = 'OHconcAfterChem'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%OHconcAfterChem( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%OHconcAfterChem = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%OHconcAfterChem,               &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE
       
       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! full-chemistry or CH4 simulations.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       diagID  = 'OHconcAfterChem'

       ! Exit if any of the above are in the diagnostic list
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '      // &
                   'but this is only appropriate for full-chemistry '    // &
                   'or CH4 simulations.' 
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !=======================================================================
    ! The following quantities are only relevant for fullchem simulations
    ! or aerosol-only simulations
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !--------------------------------------------------------------------
       ! Production of Hydrophilic BC (aka BCPI)
       ! from Hydrophobic BC (aka BCPO)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdBCPIfromBCPO'
       diagID  = 'ProdBCPIfromBCPO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdBCPIfromBCPO( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdBCPIfromBCPO = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdBCPIfromBCPO,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of Hydrophilic OC (aka OCPI)
       ! from Hydrophobic OC (aka OCPO)
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdOCPIfromOCPO'
       diagID  = 'ProdOCPIfromOCPO'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdOCPIfromOCPO( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdOCPIfromOCPO = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdOCPIfromOCPO,             &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ODDust'
       diagID  = 'OpticalDepthDust'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ODDust( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ODDust = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%ODDust,    &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 1st wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ODDustBinsWL1'
       diagID  = 'OpticalDepthDust' // TRIM(RadWL(1))
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ODDustBinsWL1( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ODDustBinsWL1 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ODDustBinsWL1,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 2nd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ODDustBinsWL2'
       diagID  = 'OpticalDepthDust' // TRIM(RadWL(2))
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ODDustBinsWL2( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ODDustBinsWL2 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ODDustBinsWL2,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Dust Optical Depth per bin at 3rd wavelength
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ODDustBinsWL3'
       diagID  = 'OpticalDepthDust' // TRIM(RadWL(3))
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ODDustBinsWL3( IM, JM, LM, NDUST ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ODDustBinsWL3 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ODDustBinsWL3,                &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 in gas phase
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromGasPhase'
       diagID  = 'ProdSO4fromGasPhase'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromGasPhase( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromGasPhase = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       & 
                                   State_Diag%ProdSO4fromGasPhase,            &
                                   State_Chm, State_Diag, RC               )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of H2O2 in cloud
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromH2O2inCloud'
       diagID  = 'ProdSO4fromH2O2inCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromH2O2inCloud( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromH2O2inCloud = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromH2O2inCloud,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of O3 in cloud
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromO3inCloud'
       diagID  = 'ProdSO4fromO3inCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO3inCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO3inCloud = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromO3inCloud,          &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by aqueous oxidation of HOBr in cloud
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromHOBrInCloud'
       diagID  = 'ProdSO4fromHOBrInCloud'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromHOBrInCloud( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromHOBrInCloud = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        & 
                                   State_Diag%ProdSO4fromHOBrInCloud,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from aqueous oxidation of O2 metal-catalyzed
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromO2inCloudMetal'
       diagID  = 'ProdSO4fromO2inCloudMetal'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          WRITE( 6, 20 ) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO2inCloudMetal(IM, JM, LM), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO2inCloudMetal = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromO2inCloudMetal,     &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 from O3 in sea salt aerosols
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSo4fromO3inSeaSalt'
       diagID  = 'ProdSo4fromO3inSeaSalt'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSo4fromO3inSeaSalt( IM, JM, LM ), STAT=RC ) 
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSo4fromO3inSeaSalt = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSo4fromO3inSeaSalt,        &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by SRO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromSRO3'
       diagID  = 'ProdSO4fromSRO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromSRO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromSRO3 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdSO4fromSRO3,                &
                                   State_Chm, State_Diag, RC               )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by SRHOBr
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromSRHOBr'
       diagID  = 'ProdSO4fromSRHOBr'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromSRHOBr( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromSRHOBr = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromSRHOBr,             &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO4 by O3s
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO4fromO3s'
       diagID  = 'ProdSO4fromO3s'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO4fromO3s( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO4fromO3s = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdSO4fromO3s,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! PM25, particulate matter with radius < 2.5 um
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%PM25'
       diagID  = 'PM25'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%PM25( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%PM25 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,     State_Diag%PM25,    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Terpene SOA
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%TerpeneSOA'
       diagID  = 'TerpeneSOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%TerpeneSOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%TerpeneSOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%TerpeneSOA,                    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Isoprene SOA
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%IsopreneSOA'
       diagID  = 'IsopreneSOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%IsopreneSOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%IsopreneSOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%IsopreneSOA,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Aromatic SOA
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%AromaticSOA'
       diagID  = 'AromaticSOA'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%AromaticSOA( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%AromaticSOA = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,     State_Diag%AromaticSOA,    &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! full-chemistry simulations or aerosol-only simulations.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 19
          
          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  ) 
                diagID = 'ProdBCPIfromBCPO'
             CASE( 2  ) 
                diagID = 'ProdOCPIfromOCPO'
             CASE( 3  ) 
                diagID = 'OpticalDepthDust'
             CASE( 4  ) 
                diagID = 'OpticalDepthDust' // TRIM( RadWL(1) )
             CASE( 5  ) 
                diagID = 'OpticalDepthDust' // TRIM( RadWL(2) )
             CASE( 6  )                
                diagID = 'OpticalDepthDust' // TRIM( RadWL(3) )
             CASE( 7  )                
                diagID = 'ProdSO4fromGasPhase'
             CASE( 8  )                
                diagID = 'ProdSO4fromH2O2inCloud'
             CASE( 9  )                
                diagID = 'ProdSO4fromO3inCloud'
             CASE( 10 )                
                diagID = 'ProdSO4fromHOBrInCloud'
             CASE( 11  )                
                diagID = 'ProdSO4fromO2inCloudMetal'
             CASE( 12 )                
                diagID = 'ProdSO4fromO3inSeaSalt'
             CASE( 13 )                
                diagID = 'ProdSO4fromSRO3'
             CASE( 14 )                
                diagID = 'ProdSO4fromSRHOBr'
             CASE( 15 )                
                diagID = 'ProdSO4fromO3s'
             CASE( 16 )                
                diagID = 'PM25'
             CASE( 17 )                
                diagID = 'TerpeneSOA'
             CASE( 18 )                
                diagID = 'IsopreneSOA'
             CASE( 19 )
                diagID = 'AromaticSOA'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for full-chemistry '  // &
                      'simulations or aerosol-only simulations.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO   
    ENDIF

    !=======================================================================
    ! The following quantities are only relevant for 
    ! aerosol-only simulations
    !=======================================================================
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !--------------------------------------------------------------------
       ! Loss of NO3 by DMS
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossNO3byDMS'
       diagID  = 'LossNO3byDMS'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossNO3byDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossNO3byDMS = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%LossNO3byDMS,                  &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Loss of OH by DMS
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%LossOHbyDMS'
       diagID  = 'LossOHbyDMS'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%LossOHbyDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%LossOHbyDMS = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%LossOHbyDMS,                   &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of MSA from DMS
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdMSAfromDMS'
       diagID  = 'ProdMSAfromDMS'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdMSAfromDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdMSAfromDMS = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                        &
                                   State_Diag%ProdMSAfromDMS,                &
                                   State_Chm, State_Diag, RC                )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Total production of SO2 from DMS
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO2fromDMS'
       diagID  = 'ProdSO2fromDMS'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMS( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMS = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdSO2fromDMS,               &
                                   State_Chm, State_Diag, RC               )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and NO3
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO2fromDMSandNO3'
       diagID  = 'ProdSO2fromDMSandNO3'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandNO3( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandNO3 = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdSO2fromDMSandNO3,         &
                                   State_Chm, State_Diag, RC               )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Production of SO2 from DMS and OH
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%ProdSO2fromDMSandOH'
       diagID  = 'ProdSO2fromDMSandOH'
       CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%ProdSO2fromDMSandOH( IM, JM, LM ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%ProdSO2fromDMSandOH = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID,                       &
                                   State_Diag%ProdSO2fromDMSandOH,          &
                                   State_Chm, State_Diag, RC               )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! aerosol-only.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 6

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1  ) 
                diagID = 'LossNO3byDMS'
             CASE( 2  ) 
                diagID = 'LossOHbyDMS'
             CASE( 3  ) 
                diagID = 'ProdMSAfromDMS'
             CASE( 4  ) 
                diagID = 'ProdSO2fromDMS'
             CASE( 5  ) 
                diagID = 'ProdSO2fromDMSandNO3'
             CASE( 6  ) 
                diagID = 'ProdSO2fromDMSandOH'
          END SELECT

          ! Exit if any of the above are in the diagnostic list
          CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
          IF ( Found ) THEN
             ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '    // &
                      'but this is only appropriate for aerosol-only '    // &
                      'simulations.' 
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO
 
    ENDIF

    !=======================================================================
    ! The production and loss diagnostics are only relevant for:
    ! (1) All full-chemistry simulations with FlexChem/KP
    ! (2) Tagged CO simulations
    ! (3) Tagged O3 simulations
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or.                                  &
         Input_Opt%ITS_A_TAGCO_SIM    .or. Input_Opt%ITS_A_TAGO3_SIM ) THEN

       !--------------------------------------------------------------------
       ! Chemical loss for selected species or families
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%Loss'
       diagID  = 'Loss'
       
       ! NOTE: Use "Loss_" as the search string so that other diagnostics
       ! such as "LossCH4byOH" won't be confused with this diagnostic.
       CALL Check_DiagList( am_I_Root, Diag_List, 'Loss_', Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Loss( IM, JM, LM, nLoss ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Loss = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%Loss,      &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

       !--------------------------------------------------------------------
       ! Chemical production for selected species or families
       !--------------------------------------------------------------------
       arrayID = 'State_Diag%Prod'
       diagID  = 'Prod'

       ! NOTE: Use "Prod_" as the search string so that other diagnostics
       ! such as "ProdBCPIfromBCPO" won't be confused with this diagnostic.
       CALL Check_DiagList( am_I_Root, Diag_List, 'Prod_', Found, RC )
       IF ( Found ) THEN
          IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
          ALLOCATE( State_Diag%Prod( IM, JM, LM, nProd ), STAT=RC )
          CALL GC_CheckVar( arrayID, 0, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
          State_Diag%Prod = 0.0_f4
          CALL Register_DiagField( am_I_Root, diagID, State_Diag%Prod,      &
                                   State_Chm, State_Diag, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF

    ELSE

       !-------------------------------------------------------------------
       ! Halt with an error message if any of the following quantities
       ! have been requested as diagnostics in simulations other than
       ! full-chemistry, tagged CO, or tagged O3 simulations.
       !
       ! This will prevent potential errors caused by the quantities
       ! being requested as diagnostic output when the corresponding
       ! array has not been allocated.
       !-------------------------------------------------------------------
       DO N = 1, 2

          ! Select the diagnostic ID
          SELECT CASE( N )
             CASE( 1 )
                diagID  = 'Loss_'
             CASE( 2 )
                diagID  = 'Prod_'
           END SELECT

           ! Exit if any of the above are in the diagnostic list
           CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
           IF ( Found ) THEN
              ErrMsg = TRIM( diagId ) // ' is a requested diagnostic, '  // &
                      'but this is only appropriate for full-chemistry, '// &
                      'tagged CO, or tagged O3 simulations.' 
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
           ENDIF
        ENDDO

    ENDIF



    ! Format statement
20  FORMAT( 1x, a32, ' is registered as: ', a )

    !-----------------------------------------------------------------
    ! TODO:
    ! 1. Hydroscopic growth - (:,:,:,N) where N is one of five hygro spc
    ! 2. Optical depth for each of five hygro spc, for each wavelength
    ! 3+ UCX-only strat diags - 5 or 7 total (hard-code)
    ! 4? isoprene optical depth??? check if AD21(:,:,:,58) is actually set
    !-----------------------------------------------------------------

    !!-------------------------------------------------------------------
    !! Template for adding more diagnostics arrays
    !! Search and replace 'xxx' with array name
    !!-------------------------------------------------------------------
    !arrayID = 'State_Diag%xxx'
    !diagID  = 'xxx'
    !CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    !IF ( Found ) THEN
    !   IF(am_I_Root) WRITE(6,20) ADJUSTL( arrayID ), TRIM( diagID )
    !   ALLOCATE( State_Diag%xxx( IM, JM, LM, n ), STAT=RC ) ! Edits dims
    !   CALL GC_CheckVar( arrayID, 0, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Diag%xxx = 0.0_f4
    !   CALL Register_DiagField( am_I_Root, diagID, State_Diag%xxx, &
    !                            State_Chm, State_Diag, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !ENDIF

    !=======================================================================
    ! Print information about the registered fields (short format)
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, 30 )
 30    FORMAT( /, &
            'Registered variables contained within the State_Diag object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( am_I_Root   = am_I_Root,                            &
                         Registry    = State_Diag%Registry,                  &
                         ShortFormat = .TRUE.,                               &
                         RC          = RC                                   )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Registry_Print"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#endif

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
!  05 Oct 2017 - R. Yantosca - Now put error trapping on deallocations
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
       CALL GC_CheckVar( 'State_Diag%SpeciesConc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
 
    IF ( ASSOCIATED( State_Diag%DryDepChm ) ) THEN
       DEALLOCATE( State_Diag%DryDepChm, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepChm', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepMix ) ) THEN
       DEALLOCATE( State_Diag%DryDepMix, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepMix', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDep ) ) THEN
       DEALLOCATE( State_Diag%DryDep, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDep', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%DryDepVel ) ) THEN
       DEALLOCATE( State_Diag%DryDepVel, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%DryDepVel', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%JVal ) ) THEN
       DEALLOCATE( State_Diag%JVal, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%Jval', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%JNoonDailyAvg ) ) THEN
       DEALLOCATE( State_Diag%JNoonDailyAvg, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%JNoonDailyAvg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%JNoonMonthlyAvg ) ) THEN
       DEALLOCATE( State_Diag%JNoonMonthlyAvg, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%JNoonMonthlyAvg', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RxnRates ) ) THEN
       DEALLOCATE( State_Diag%RxnRates, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RxnRates', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxDiffuse ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDiffuse, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%UvFluxDiffuse', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxDirect ) ) THEN
       DEALLOCATE( State_Diag%UVFluxDirect, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%UvFluxDirect', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%UVFluxNet ) ) THEN
       DEALLOCATE( State_Diag%UVFluxNet, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%UvFluxNet', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxZonal ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxZonal, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AdvFluxZonal', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxMerid ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxMerid, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AdvFluxMerid', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AdvFluxVert ) ) THEN
       DEALLOCATE( State_Diag%AdvFluxVert, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AdvFluxVert', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLMixFrac ) ) THEN
       DEALLOCATE( State_Diag%PBLMixFrac, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PBLMixFrac', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PBLFlux ) ) THEN
       DEALLOCATE( State_Diag%PBLFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%PBLFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%CloudConvFlux ) ) THEN
       DEALLOCATE( State_Diag%CloudConvFlux, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%CloudConvFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossConv ) ) THEN
       DEALLOCATE( State_Diag%WetLossConv, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%WetLossConv', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%WetLossLS ) ) THEN
       DEALLOCATE( State_Diag%WetLossLS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%WetLossLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PrecipFracLS ) ) THEN
       DEALLOCATE( State_Diag%PrecipFracLS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PrecipFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RainFracLS ) ) THEN
       DEALLOCATE( State_Diag%RainFracLS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RainFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%WashFracLS ) ) THEN
       DEALLOCATE( State_Diag%WashFracLS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%WashFracLS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PbFromRnDecay ) ) THEN
       DEALLOCATE( State_Diag%PbFromRnDecay, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PbFromRnDecay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadDecay ) ) THEN
       DEALLOCATE( State_Diag%RadDecay, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadDecay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWSurf, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadAllSkyLWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkyLWTOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadAllSkyLWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWSurf, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadAllSkySWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadAllSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadAllSkySWTOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadAllSkySWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWSurf, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadClrSkyLWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkyLWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkyLWTOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadClrSkyLWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWSurf ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWSurf, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadClrSkySWSurf', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%RadClrSkySWTOA ) ) THEN
       DEALLOCATE( State_Diag%RadClrSkySWTOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%RadClrSkySWTOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdBCPIfromBCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdBCPIfromBCPO, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdBCPIfromBCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdOCPIfromOCPO ) ) THEN
       DEALLOCATE( State_Diag%ProdOCPIfromOCPO, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdOCPIfromOCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%OHconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%OhconcAfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%OHconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 

    IF ( ASSOCIATED( State_Diag%HO2concAfterChem ) ) THEN
       DEALLOCATE( State_Diag%HO2concAfterChem, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%HO2concAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 

    IF ( ASSOCIATED( State_Diag%O1DconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O1DconcAfterChem, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%O1DconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%O3PconcAfterChem ) ) THEN
       DEALLOCATE( State_Diag%O3PconcAfterChem, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%O3PconcAfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ODDust ) ) THEN
       DEALLOCATE( State_Diag%ODDust, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ODDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ODDustBinsWL1 ) ) THEN
       DEALLOCATE( State_Diag%ODDustBinsWL1, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ODDustBinsWL1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ODDustBinsWL2 ) ) THEN
       DEALLOCATE( State_Diag%ODDustBinsWL2, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ODDustBinsWL2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ODDustBinsWL3 ) ) THEN
       DEALLOCATE( State_Diag%ODDustBinsWL3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ODDustBinsWL3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%Loss ) ) THEN
       DEALLOCATE( State_Diag%Loss, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%Loss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%Prod ) ) THEN
       DEALLOCATE( State_Diag%Prod, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%Prod', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMSandOH ) ) THEN
       DEALLOCATE( State_Diag%ProdSO2fromDMSandOH, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO2fromDMSandOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMSandNO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdSO2fromDMSandNO3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO2fromDMSandNO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO2fromDMS ) ) THEN
       DEALLOCATE(  State_Diag%ProdSO2fromDMS, STAT=RC  )
       CALL GC_CheckVar( ' State_Diag%ProdSO2fromDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdMSAfromDMS ) ) THEN
       DEALLOCATE( State_Diag%ProdMSAfromDMS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdMSAfromDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdNITfromHNO3uptakeOnDust ) ) THEN
       DEALLOCATE( State_Diag%ProdNITfromHNO3uptakeOnDust, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdNITfromHNO3uptakeOnDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromGasPhase ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromGasPhase, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromGasPhase', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromH2O2inCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromH2O2inCloud, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromH2O2inCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3inCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3inCloud, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3inCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromHOBrInCloud ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromHOBrInCloud, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromHOBrInCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO2inCloudMetal ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO2inCloudMetal, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO2inCloudMetal', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3inSeaSalt ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3inSeaSalt, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3inSeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromOxidationOnDust  ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromOxidationOnDust, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromOxidationOnDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromUptakeOfH2SO4g ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromUptakeOfH2SO4g, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromUptakeOfH2SO4g', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromSRO3 ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromSRO3, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromSRO3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromSRHOBr ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromSRHOBr, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromSRHOBr', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%ProdSO4fromO3s ) ) THEN
       DEALLOCATE( State_Diag%ProdSO4fromO3s, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%ProdSO4fromO3s', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossOHbyDMS ) ) THEN
       DEALLOCATE( State_Diag%LossOHbyDMS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%LossOHbyDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossNO3byDMS ) ) THEN
       DEALLOCATE( State_Diag%LossNO3byDMS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%LossNO3byDMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%LossHNO3onSeaSalt  ) ) THEN
       DEALLOCATE( State_Diag%LossHNO3onSeaSalt, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%LossHNO3onSeaSalt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%StratSADS  ) ) THEN
       DEALLOCATE( State_Diag%StratSADS, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%StratSADS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%StratSADL  ) ) THEN
       DEALLOCATE( State_Diag%StratSADL, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%StratSADL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%PM25 ) ) THEN
       DEALLOCATE( State_Diag%PM25, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%PM25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%TerpeneSOA ) ) THEN
       DEALLOCATE( State_Diag%TerpeneSOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%TerpeneSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%IsopreneSOA ) ) THEN
       DEALLOCATE( State_Diag%IsopreneSOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%IsopreneSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Diag%AromaticSOA ) ) THEN
       DEALLOCATE( State_Diag%AromaticSOA, STAT=RC  )
       CALL GC_CheckVar( 'State_Diag%AromaticSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF


    !-----------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
    !-----------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Diag%xxx ) ) THEN
    !   DEALLOCATE( State_Diag%xxx, STAT=RC  )
    !   CALL GC_CheckVar( 'State_Diag%xxx', 2, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
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
                                      TagId,      Rank,       Type,     &
                                      VLoc                             )
!
! !USES:
!
    USE Charpak_Mod,         ONLY: StrSplit, To_UpperCase
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
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: TagId      ! Tag wildcard (wc)
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank       ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type       ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc       ! Vert placement
!
! !REMARKS:
!  If a diagnostic cannot use a wildcard, then set Tag=''.
!
! !REVISION HISTORY: 
!  20 Sep 2017 - E. Lundgren - Initial version
!  06 Oct 2017 - R. Yantosca - State_Diag%SpeciesConc is now an 8-byte real
!  01 Nov 2017 - R. Yantosca - Now get To_UpperCase from charpak_mod.F90
!  02 Nov 2017 - R. Yantosca - Update metadata to be consistent w/ arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: isDesc,  isUnits,  isRank
    LOGICAL            :: isVLoc,  isTagged, isType

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Name_AllCaps
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC        =  GC_SUCCESS
    Found     = .TRUE.
    ErrMsg    = ''
    ThisLoc   =  &
         ' -> at Get_Metadata_State_Diag (in Headers/state_diag_mod.F90)'

    ! Optional arguments present?
    isDesc    = PRESENT( Desc    )
    isUnits   = PRESENT( Units   )
    isRank    = PRESENT( Rank    )
    isType    = PRESENT( Type    )
    isVLoc    = PRESENT( VLoc    )
    isTagged  = PRESENT( TagID   ) 

    ! Set defaults for optional arguments. Assume type and vertical 
    ! location are real (flexible precision) and center unless specified 
    ! otherwise
    IF ( isUnits  ) Units  = ''
    IF ( isDesc   ) Desc   = ''              
    IF ( isRank   ) Rank   = -1 
    IF ( isType   ) Type   = KINDVAL_F4      ! Assume real*4
    IF ( isVLoc   ) VLoc   = VLocationCenter ! Assume vertically centered
    IF ( isTagged ) TagID  = '' 

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    IF ( TRIM( Name_AllCaps ) == 'SPECIESCONC' ) THEN
       IF ( isDesc    ) Desc  = 'Dry mixing ratio of species'
       IF ( isUnits   ) Units = 'mol mol-1 dry'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'
       IF ( isType    ) Type  = KINDVAL_F8

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPCHM' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species, from chemistry'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPMIX' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species, from mixing'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEP' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition flux of species'
       IF ( isUnits   ) Units = 'molec cm-2 s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'DRYDEPVEL' ) THEN
       IF ( isDesc    ) Desc  = 'Dry deposition velocity of species'
       IF ( isUnits   ) Units = 'cm s-1'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'DRY'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JVAL' ) THEN
       IF ( isDesc    ) Desc  = 'Photolysis rate for species' 
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOONDAILYAVG' ) THEN
       IF ( isDesc    ) Desc  = 'Noontime photolysis rate for species' 
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'JNOONMONTHLYAVG' ) THEN
       IF ( isDesc    ) Desc  = 'Noontime photolysis rate for species' 
       IF ( isUnits   ) Units = 's-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PHO'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RXNRATES' ) THEN
       IF ( isDesc    ) Desc  = 'placeholder'
       IF ( isUnits   ) Units = 'placeholder'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ALL'

    ELSE IF ( TRIM( Name_AllCaps ) == 'UVFLUXDIFFUSE' ) THEN
       IF ( isDesc    ) Desc  = 'placeholder'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'UVFLUXDIRECT' ) THEN
       IF ( isDesc    ) Desc  = 'placeholder'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3

    ELSEIF ( TRIM( Name_AllCaps ) == 'UVFLUXNET' ) THEN
       IF ( isDesc    ) Desc  = 'placeholder'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'ADVFLUXZONAL' ) THEN
       IF ( isDesc    ) Desc  = 'Advection of species in zonal direction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'ADVFLUXMERID' ) THEN
       IF ( isDesc    ) Desc  = 'Advection of species in meridional direction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'
     
    ELSE IF ( TRIM( Name_AllCaps ) == 'ADVFLUXVERT' ) THEN
       IF ( isDesc    ) Desc  = 'Advection of species in vertical direction'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBLMIXFRAC' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of boundary layer occupied by each level'
       IF ( isUnits   ) Units = 'placeholder'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBLFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'Species mass change due to boundary-layer mixing'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'CLOUDCONVFLUX' ) THEN
       IF ( isDesc    ) Desc  = 'Mass change due to cloud convection'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of soluble species in convective updrafts'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRECIPFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of grid box undergoing convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RAINFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to rainout in convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WASHFRACCONV' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to washout in convective precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WETLOSSLS' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of soluble species in large-scale precipitation'
       IF ( isUnits   ) Units = 'placeholder'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRECIPFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of grid box undergoing large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RAINFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to rainout in large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'WASHFRACLS' ) THEN
       IF ( isDesc    ) Desc  = 'Fraction of soluble species lost to washout in large-scale precipitation'
       IF ( isUnits   ) Units = '1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'WET'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PBFROMRNDECAY' ) THEN
       IF ( isDesc    ) Desc  = 'Pb210 created from radioactive decay of Rn222'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADDECAY' ) THEN
       IF ( isDesc    ) Desc  = 'Radioactive decay of radionuclide species'
       IF ( isUnits   ) Units = 'kg s-1'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'ADV'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYLWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky long-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYLWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky long-wave radiation at top of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSEIF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADALLSKYSWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'All-sky short-wave radiation at top of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYLWTOA ' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky long-wave radiation at top of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWSURF' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at surface'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId =  'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'RADCLRSKYSWTOA' ) THEN
       IF ( isDesc    ) Desc  = 'Clear-sky short-wave radiation at top of atmosphere'
       IF ( isUnits   ) Units = 'W m-2'
       IF ( isRank    ) Rank  = 2
       IF ( isTagged  ) TagId = 'placeholder'

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODBCPIFROMBCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of hydrophilic black carbon from hydrophobic black carbon'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODOCPIFROMOCPO' ) THEN
       IF ( isDesc    ) Desc  = 'Production of hydrophilic organic carbon from hydrophobic organic carbon'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'OHCONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'OH concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'HO2CONCAFTERCHEM' )  THEN
       IF ( isDesc    ) Desc  = 'HO2 concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'mol mol-1'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'O1DCONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'O1D concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'O3PCONCAFTERCHEM' ) THEN
       IF ( isDesc    ) Desc  = 'O3P concentration immediately after chemistry'
       IF ( isUnits   ) Units = 'molec cm-3'
       IF ( isRank    ) Rank  = 3

    ELSE IF ( TRIM( Name_AllCaps ) == 'OPTICALDEPTHDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Optical depth for mineral dust'
       IF ( isUnits   ) Units = 'unitless'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'OPTICALDEPTHDUST' // TRIM(RadWL(1)) ) THEN 
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(1)) // ' nm'
       IF ( isUnits   ) Units   = 'unitless'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'OPTICALDEPTHDUST' // TRIM(RadWL(2)) ) THEN
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(2)) // ' nm'
       IF ( isUnits   ) Units   = 'unitless'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'OPTICALDEPTHDUST' // TRIM(RadWL(3)) ) THEN
       IF ( isDesc    ) Desc    = 'Optical depth for dust at ' // &
                                   TRIM(RadWL(3)) // ' nm'
       IF ( isUnits   ) Units   = 'unitless'
       IF ( isRank    ) Rank    =  3
       IF ( isTagged  ) TagId   = 'DUSTBIN'

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSS' ) THEN
       IF ( IsDesc    ) Desc  = 'Chemical loss of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'LOS'

       ! NOTE: Units are different depending on simulation, due to historical
       ! baggage.  Maybe clean this up at a later point to use the same units
       ! regardless of simulation type. (bmy, 12/4/17)
       IF ( isUnits   ) THEN
          IF ( IsFullChem ) THEN
             Units = 'molec/cm3/s'
          ELSE
             Units = 'kg/s'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PROD' ) THEN
       IF ( isDesc    ) Desc  = 'Chemical production of'
       IF ( isRank    ) Rank  = 3
       IF ( isTagged  ) TagId = 'PRD'

       ! NOTE: Units are different depending on simulation, due to historical
       ! baggage.  Maybe clean this up at a later point to use the same units
       ! regardless of simulation type. (bmy, 12/4/17)
       IF ( isUnits   ) THEN
          IF ( IsFullChem ) THEN
             Units = 'molec/cm3/s'
          ELSE
             Units = 'kg/s'
          ENDIF
       ENDIF

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMSANDOH' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO2 from DMS+OH reaction'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMSANDNO3' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO2 from DMS+NO3 reaction'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO2FROMDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Total production of SO2 from DMS'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODMSAFROMDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Production of MSA from DMS'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMGASPHASE' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from gas phase reactions'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMH2O2INCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous oxidation of H2O2 in clouds'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3INCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous oxidation of O3 in clouds'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMHOBRINCLOUD' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous oxidation of HOBr in clouds'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO2INCLOUDMETAL' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from aqueous oxidation of O2 metal-catalyzed'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3INSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from O3 in sea salt aerosols'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMOXIDATIONONDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from oxidation on dust aerosols'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODNITFROMHNO3UPTAKEONDUST' ) THEN
       IF ( isDesc    ) Desc  = 'Production of NIT from HNO3 uptake on dust aerosols'
       IF ( isUnits   ) Units = 'kg N'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMUPTAKEOFH2SO4G' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from uptake of H2SO4(g)'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMSRO3' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 by SRO3'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMSRHOBR' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from SRHOBr'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PRODSO4FROMO3S' ) THEN
       IF ( isDesc    ) Desc  = 'Production of SO4 from O3s'
       IF ( isUnits   ) Units = 'kg S'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSOHBYDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of OH by DMS'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSNO3BYDMS' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of NO3 by DMS'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'LOSSHNO3ONSEASALT' ) THEN
       IF ( isDesc    ) Desc  = 'Loss of HNO3 on sea salt aerosols'
       IF ( isUnits   ) Units = 'kg'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'STRATSADL' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric liquid aerosol surface aera density'
       IF ( isUnits   ) Units = 'cm2/cm3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'STRATSADS' ) THEN
       IF ( isDesc    ) Desc  = 'Stratospheric solid aerosol surface aera density'
       IF ( isUnits   ) Units = 'cm2/cm3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'PM25' ) THEN
       IF ( isDesc    ) Desc  = 'Particulate matter with radii < 2.5 um'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'TERPENESOA' ) THEN
       IF ( isDesc    ) Desc  = 'Monoterpene and sesqiterpene SOA'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'ISOPRENESOA' ) THEN
       IF ( isDesc    ) Desc  = 'Isoprene (biogenic) SOA from either semivolatile partitioning or irreversible uptake'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE IF ( TRIM( Name_AllCaps ) == 'AROMATICSOA' ) THEN
       IF ( isDesc    ) Desc  = 'Aromatic and intermediate volatility (anthropogenic) SOA'
       IF ( isUnits   ) Units = 'ug m-3'
       IF ( isRank    ) Rank  =  3

    ELSE
       
       !--------------------------------------------------------------------
       ! Could not find metadata, so exit with error message
       !--------------------------------------------------------------------
       Found = .False.
       ErrMsg = 'Metadata not found for State_Diag field ID: '            // &
                 TRIM( metadataID ) // '. If the name in HISTORY.rc '     // &
                'has species appended, make sure the species name '       // &
                'is preceded by a single underscore. Otherwise, '         // &
                'check that the name is listed with all capitals in '     // &
                'subroutine Get_Metadata_State_Diag '                     // &
                '(Headers/state_diag_mod.F90).'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

   END SUBROUTINE Get_Metadata_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_TagInfo
!
! !DESCRIPTION: Subroutine GET\_TAGINFO retrieves basic information about 
! tags given a wildcard string.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_TagInfo( am_I_Root, tagID, State_Chm, Found,                &
                          RC,        N,     tagName,   nTags                )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,            INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
    CHARACTER(LEN=*),   INTENT(IN)  :: tagID       ! ID of tag (e.g. wildcard)
    TYPE(ChmState),     INTENT(IN)  :: State_Chm   ! Chemistry State object
    INTEGER,            OPTIONAL    :: N           ! index (1 to # tags)
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,            INTENT(OUT) :: Found       ! Item found?
    INTEGER,            INTENT(OUT) :: RC          ! Return code
    CHARACTER(LEN=255), OPTIONAL    :: tagName     ! tag name for index N 
    INTEGER,            OPTIONAL    :: nTags       ! # tags
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  16 Nov 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: D,         numTags
    LOGICAL            :: isNumTags, isTagName, isN

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,    ThisLoc,   Nstr

    !=======================================================================
    ! Get_TagName begins here
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    ErrMsg     = ''
    ThisLoc    = ' -> at Get_TagInfo (in Headers/state_diag_mod.F90)'
    Found      = .TRUE.
    numTags    = 0
    
    ! Optional arguments present?
    isN        = PRESENT( N       )
    isTagName  = PRESENT( TagName )
    isNumTags  = PRESENT( nTags   )

    ! Exit with error if getting tag name but index not specified
    IF ( isTagName .AND. .NOT. isN ) THEN
       ErrMsg = 'Index must be specified if retrieving an individual tag name'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Get number of tags
    !=======================================================================
    SELECT CASE( TRIM( tagId ) )
       CASE( 'ALL'     )
          numTags = State_Chm%nSpecies
       CASE( 'ADV'     )
          numTags = State_Chm%nAdvect
       CASE( 'AER'     )
          numTags = State_Chm%nAero
       CASE( 'DRY'     )
          numTags = State_Chm%nDryDep
       CASE( 'DUSTBIN' )
          numTags = NDUST
       CASE( 'FIX'     )
          numTags = State_Chm%nKppFix
       CASE( 'GAS'     )
          numTags = State_Chm%nGasSpc
       CASE( 'HYG'     )
          numTags = State_Chm%nHygGrth
       CASE( 'KPP'     )
          numTags = State_Chm%nKppSpc
       CASE( 'LOS'     )
          numTags = State_Chm%nLoss
       CASE( 'PHO'     )
          numTags = State_Chm%nPhotol+2  ! NOTE: Extra slots for diagnostics
       CASE( 'PRD'     )
          numTags = State_Chm%nProd
       CASE( 'VAR'     )
          numTags = State_Chm%nKppVar
       CASE( 'WET'     )
          numTags = State_Chm%nWetDep
       CASE DEFAULT
          FOUND = .FALSE.
          ErrMsg = 'Handling of tagId ' // TRIM(tagId) // &
                   ' is not implemented for getting number of tags'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

    !=======================================================================
    ! Sanity checks -- exit under certain conditions
    !=======================================================================

    ! If not getting tag name then set nTags and exit
    IF ( .NOT. isTagName ) THEN 
       nTags = numTags
       RETURN
    ENDIF

    ! Exit with error if index exceeds number of tags for this wildcard
    IF ( isTagName .AND. .NOT. isN ) THEN
       ErrMsg = 'Index must be greater than total number of tags for wildcard' &
                // TRIM(tagId)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Get mapping index
    !=======================================================================
    SELECT CASE( TRIM( tagID ) )
       CASE( 'ALL', 'ADV', 'DUSTBIN', 'PRD', 'LOS' )
          D = N
       CASE( 'AER'  )
          D = State_Chm%Map_Aero(N)
       CASE( 'DRY'  )
          D = State_Chm%Map_DryDep(N)
       CASE( 'GAS'  )
          D = State_Chm%Map_GasSpc(N)
       CASE( 'HYG'  )
          D = State_Chm%Map_HygGrth(N)
       CASE( 'VAR'  )
          D = State_Chm%Map_KppVar(N)
       CASE( 'FIX'  )
          D = State_Chm%Map_KppFix(N)
       CASE( 'KPP'  )
          D = State_Chm%Map_KppSpc(N)
       CASE( 'PHO'  )
          IF ( N > State_Chm%nPhotol ) THEN  ! NOTE: To denote the nPhotol+1
             D = 5000 + N                    ! and nPhotol+2 slots, add a
          ELSE                               ! large offset, so that we don't
             D = State_Chm%Map_Photol(N)     ! clobber any existing species #'s
          ENDIF
       CASE( 'WET'  )
          D = State_Chm%Map_WetDep(N)
       CASE DEFAULT
          FOUND = .FALSE.
          ErrMsg = 'Handling of tagId ' // TRIM( tagId ) // &
                   ' is not implemented for getting tag name'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

    !=======================================================================
    ! Return the tag name
    !=======================================================================

    ! Initialize
    tagName = ''  

    ! Special handling for certain tagID's
    SELECT CASE( TRIM( tagID ) )

       ! Dust bins
       CASE( 'DUSTBIN' )
          WRITE ( Nstr, "(I1)" ) D
          tagName = 'BIN' // TRIM(Nstr)

       ! Loss species
       CASE( 'LOS' ) 
          tagName = State_Chm%Name_Loss(N)
          D       = INDEX( tagName, '_' )
          tagName = tagName(D+1:)

       ! Prod species
       CASE( 'PRD' )
          tagName = State_Chm%Name_Prod(N)
          D       = INDEX( tagName, '_' )
          tagName = tagName(D+1:)

       ! Photolysis species
       CASE( 'PHO' )

          ! Save O1_O3D (UCX) or O3 (non-UCX) in the nPhotol+1 slot
          IF ( D == 5001 + State_Chm%nPhotol ) THEN
             IF ( Is_UCX ) THEN
                tagName = 'O3O1D'
             ELSE
                tagName = 'O3O1Da'
             ENDIF

          ! Save O3_O3P (UCX) or POH (non UCX) in the nPhotol+2 slot
          ELSE IF ( D == 5002 + State_Chm%nPhotol ) THEN
             IF ( Is_UCX ) THEN
                tagName = 'O3O3P'
             ELSE
                tagName = 'O3O1Db'
             ENDIF
          
          ! For all other photolysis species, get the name
          ! from the GEOS-Chem species database
          ELSE
             tagName = State_Chm%SpcData(D)%Info%Name
             
          ENDIF

       ! Default tag name is the name in the species database
       CASE DEFAULT
          tagName = State_Chm%SpcData(D)%Info%Name

    END SELECT

  END SUBROUTINE Get_TagInfo
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
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: found

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    ! Get metadata for this diagnostic
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,   tagId=tagId      )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata dimensions consistent with data pointer
    IF ( ( ( tagId == '' ) .AND. ( rank /= 2 ) )  &
         .OR. ( ( tagId /= '' ) .AND. ( rank /= 1 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
          
    ! Special handling if there are tags (wildcard)
    IF ( tagId /= '' ) THEN
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                         nTags=nTags )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags          
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )
          IF ( RC /= GC_SUCCESS ) THEN
             CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
             RETURN
          ENDIF
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc ) // ' ' // TRIM( tagName )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = diagName,             &
                                  Description  = diagDesc,             &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where tagID is ' // TRIM(tagID)
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
    CHARACTER(LEN=255)     :: desc, units, tagID, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type,  vloc
    LOGICAL                :: Found

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    ! Get metadata for this diagnostic
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  tagID=tagID                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata dimensions consistent with data pointer
    IF ( ( ( tagID == '' ) .AND. ( rank /= 3 ) )  &
         .OR. ( ( tagID /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Special handling if there are tags
    IF ( tagID /= '' ) THEN
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                         nTags=nTags )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags          
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )
          IF ( RC /= GC_SUCCESS ) THEN
             CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
             RETURN
          ENDIF
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc ) // ' ' // TRIM( tagName )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = diagName,             &
                                  Description  = diagDesc,             &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),      &
                                  RC           = RC                   )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where tagID is ' // TRIM(tagID)
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
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: found

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_4D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    ! Get metadata for this diagnostic
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,   tagId=tagId  )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata dimensions consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Assume always tagged if 4D
    CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                      nTags=nTags )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Register each tagged name as a separate diagnostic
    DO N = 1, nTags          
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                         N=N, tagName=tagName )

       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF
       diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
       diagDesc = TRIM( Desc ) // ' ' // TRIM( tagName )
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = diagName,             &
                               Description  = diagDesc,             &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data(:,:,:,N),    &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where tagId is ' // TRIM(tagId)
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
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type, vloc
    LOGICAL                :: Found

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    ! Get metadata for this diagnostic
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  tagId=tagId                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata dimensions consistent with data pointer
    IF ( ( ( tagId == '' ) .AND. ( rank /= 3 ) )  &
         .OR. ( ( tagId /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Special handling if there are tags
    IF ( tagId /= '' ) THEN
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                         nTags=nTags )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Register each tagged name as a separate diagnostic
       DO N = 1, nTags          
          CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                            N=N, tagName=tagName )
          IF ( RC /= GC_SUCCESS ) THEN
             CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
             RETURN
          ENDIF
          diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
          diagDesc = TRIM( Desc ) // ' ' // TRIM( tagName )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = diagName,             &
                                  Description  = diagDesc,             &
                                  Units        = units,                &
                                  Data2d       = Ptr2Data(:,:,N),      &
                                  RC           = RC                   )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where tagId is ' // TRIM(tagId)
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
    CHARACTER(LEN=255)     :: desc, units, tagId, tagName
    CHARACTER(LEN=255)     :: diagName, diagDesc
    INTEGER                :: N, nTags, rank, type,  vloc
    LOGICAL                :: found

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R8_4D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    ! Get metadata for this diagnostic
    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  tagId=tagId                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata dimensions consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Assume always tagged. Get number of tags.
    CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                      nTags=nTags )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Register each tagged name as a separate diagnostic
    DO N = 1, nTags          
       CALL Get_TagInfo( am_I_Root, tagId, State_Chm, Found, RC, &
                         N=N, tagName=tagName )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF
       diagName = TRIM( metadataID ) // '_' // TRIM( tagName )
       diagDesc = TRIM( Desc ) // ' ' // TRIM( tagName )
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = diagName,             &
                               Description  = diagDesc,             &
                               Units        = units,                &
                               Data3d_8     = Ptr2Data(:,:,:,N),    &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where tagId is ' // TRIM(tagId)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE Register_DiagField_R8_4D
!EOC
END MODULE State_Diag_Mod
